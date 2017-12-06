#include "dna_builder.h"

using namespace std;

DNABuilder::DNABuilder() :
	curr_sequence_index(0),
	curr_zero_index(0),
	curr_one_index(0),
	buff_number_of_zeroes(0),
	buff_number_of_ones(0) {}

DNABuilder::DNABuilder(const vector<uchar> * file, const SequenceOptions * options) :
	input(file),
	options(options),
	curr_sequence_index(0),
	curr_zero_index(0),
	curr_one_index(0),
	buff_number_of_zeroes(0),
	buff_number_of_ones(0) {}

string DNABuilder::construct_sequence()
{
	string sequence;
	make_nucleotide_map();
	add_pam(sequence);
	add_header(sequence);
	add_content(sequence);
	add_padding(sequence);
	add_marker_bit(sequence);
	sequence += options->get_trailing_invariant();
	return sequence;
}

vector<string> & DNABuilder::construct_sequences()
{
	uint input_length = input->size(); // Input length in bytes
	uint content_length = options->get_content_length(); // Content length in bits
	vector<uchar> bin_buffer; // Binary buffer of content
	uint i = 0; // Byte index
	while (i < input_length)
	{
		get_uchar_binary(input->at(i), bin_buffer); // Get the binary digits of the current byte
		i++;
		if (bin_buffer.size() >= content_length) // Check if there are enough binary digits to construct a sequence
		{
			buffer.insert(buffer.end(), bin_buffer.begin(), bin_buffer.begin() + content_length);
			bin_buffer.erase(bin_buffer.begin(), bin_buffer.begin() + content_length);
			dna_sequences.push_back(construct_sequence());
			reset_current_sequence();
			curr_sequence_index++;
		}
	}
	if (bin_buffer.size() > 0) // If the buffer is not empty (it means our input length in bits was not a multiple of <content_length>
	{
		buffer.insert(buffer.end(), bin_buffer.begin(), bin_buffer.begin() + bin_buffer.size());
		dna_sequences.push_back(construct_sequence());
		reset_current_sequence();
		curr_sequence_index++;
	}
	return dna_sequences;
}

void DNABuilder::reset_current_sequence()
{
	nucleotide_to_digit.clear();
	zero_digit_to_nucleotide.clear();
	one_digit_to_nucleotide.clear();
	buffer.clear();
	curr_zero_index = 0;
	curr_one_index = 0;
	buff_number_of_zeroes = 0;
	buff_number_of_ones = 0;
}

void DNABuilder::reset_builder()
{
	reset_current_sequence();
	input = nullptr;
	options = nullptr;
	curr_sequence_index = 0;
	dna_sequences.clear();
}

void DNABuilder::make_nucleotide_map()
{
	// Analyze current buffer (count the number of binary 0s and 1s in the content)
	for (uchar c : buffer)
	{
		(c == '0') ? buff_number_of_zeroes++ : buff_number_of_ones++;
	}
	/**
	 *	- Make nucleotide map -
	 *	In order to keep GC% content >= 50%, assign GC to 
	 *	represent the binary digits which appear the most.
	*/
	if (buff_number_of_zeroes < buff_number_of_ones)
	{
		nucleotide_to_digit['A'] = 0;
		nucleotide_to_digit['T'] = 0;
		nucleotide_to_digit['G'] = 1;
		nucleotide_to_digit['C'] = 1;
		zero_digit_to_nucleotide.push_back('A');
		zero_digit_to_nucleotide.push_back('T');
		one_digit_to_nucleotide.push_back('G');
		one_digit_to_nucleotide.push_back('C');
		curr_marker_bit = 'G'; // Not definite - can be modified
	}
	else
	{
		nucleotide_to_digit['A'] = 1;
		nucleotide_to_digit['T'] = 1;
		nucleotide_to_digit['G'] = 0;
		nucleotide_to_digit['C'] = 0;
		zero_digit_to_nucleotide.push_back('G');
		zero_digit_to_nucleotide.push_back('C');
		one_digit_to_nucleotide.push_back('A');
		one_digit_to_nucleotide.push_back('T');
		curr_marker_bit = 'A'; // Not definite - can be modified
	}
}

void DNABuilder::add_pam(string & sequence)
{
	string pam = options->get_pam_invariant();
	// Analyze PAM (update 0's and 1's indexes)
	uchar prev_AT = 'X';
	uchar prev_GC = 'X';
	for (uchar c : pam)
	{
		if (c == 'A' || c == 'T') // check for repeats
		{
			if (c != prev_AT)
			{
				prev_AT = c;
			}
			else
			{
				continue;
			}
		}
		else
		{
			if (c != prev_GC)
			{
				prev_GC = c;
			}
			else
			{
				continue;
			}
		}
		if (nucleotide_to_digit[c] == 0)
		{
			curr_zero_index++;
		}
		else
		{
			curr_one_index++;
		}
	}
	sequence += pam;
}

void DNABuilder::add_header(string & sequence)
{
	string header;
	// Get binary digits of sequence index
	vector<uchar> binary_seq_index;
	get_uint_binary(curr_sequence_index, binary_seq_index);
	// Remove redundance (only keep <sequence_index_length> bits)
	binary_seq_index.erase(binary_seq_index.begin(), binary_seq_index.begin() + binary_seq_index.size() - options->get_sequence_index_length());
	uint seq_index_zeroes = 0, seq_index_ones = 0;
	// Count binary 0s and 1s
	for (uchar c : binary_seq_index)
	{
		(c == '0') ? seq_index_zeroes++ : seq_index_ones++;
	}
	// Add parity bit
	uchar parity_bit;
	if (nucleotide_to_digit['G'] == 0)
	{
		if ((buff_number_of_zeroes + seq_index_zeroes) % 2 == 0)
		{
			parity_bit = zero_digit_to_nucleotide[curr_zero_index % 2];
			curr_zero_index++;
		}
		else
		{
			parity_bit = one_digit_to_nucleotide[curr_one_index % 2];
			curr_one_index++;
		}
	}
	else
	{
		if ((buff_number_of_ones + seq_index_ones) % 2 == 0)
		{
			parity_bit = one_digit_to_nucleotide[curr_one_index % 2];
			curr_one_index++;
		}
		else
		{
			parity_bit = zero_digit_to_nucleotide[curr_zero_index % 2];
			curr_zero_index++;
		}
	}
	header += parity_bit;
	// Add sequence index
	string dna_seq_index;
	for (uchar digit : binary_seq_index)
	{
		if (digit == '0')
		{
			dna_seq_index += zero_digit_to_nucleotide[curr_zero_index % 2];
			curr_zero_index++;
		}
		else
		{
			dna_seq_index += one_digit_to_nucleotide[curr_one_index % 2];
			curr_one_index++;
		}
	}
	header += dna_seq_index;
	sequence += header;
}

void DNABuilder::add_content(string & sequence)
{
	// Translate bits to nucleotides
	for (uchar c : buffer)
	{
		if (c == '0')
		{
			sequence += zero_digit_to_nucleotide[curr_zero_index % 2];
			curr_zero_index++;
		}
		else
		{
			sequence += one_digit_to_nucleotide[curr_one_index % 2];
			curr_one_index++;
		}
	}
}

void DNABuilder::add_padding(string & sequence)
{
	// Check if padding is necessary
	if (buffer.size() < options->get_content_length())
	{
		/**
		* Re-add the last nucleotide as a repeat for padding detection.
		* Note: This algorithm is designed not to permit any mononucleotide repeats.
		* Using this method, there's no need to specify input length within the sequence.
		*/
		uchar repeat = sequence[sequence.size() - 1];
		sequence += repeat;
		uint padding_index = 0;
		for (uint i = buffer.size() + 1; i < options->get_content_length(); i++)
		{
			// Add further GC-rich padding
			if (nucleotide_to_digit['G'] == 0)
			{
				sequence += zero_digit_to_nucleotide[(curr_zero_index + padding_index) % 2];
			}
			else
			{
				sequence += one_digit_to_nucleotide[(curr_one_index + padding_index) % 2];
			}
			padding_index++;
		}
	}
}

void DNABuilder::add_marker_bit(string & sequence)
{
	/**
	* In case the last nucleotide of the content is the same as the bit marker,
	* in order to avoid 5'-3' internal PAM, A gets swapped with T and G gets swapped with C.
	* Note: Canonical PAM is guaranteed to contain mononucleotide repeats.
	*/
	if (sequence[sequence.size() - 1] == curr_marker_bit)
	{
		if (nucleotide_to_digit['G'] == 0)
		{
			curr_marker_bit = 'T';
		}
		else
		{
			curr_marker_bit = 'C';
		}
	}
	sequence += curr_marker_bit;
}

void DNABuilder::count_binary_digits(const uchar c, uint & zeroes, uint & ones)
{
	for (uint i = 0; i < 8; i++)
	{
		((c >> i) & 0x01) ? ones++ : zeroes++;
	}
}

void DNABuilder::count_binary_digits(const uint c, uint & zeroes, uint & ones)
{
	for (uint i = 0; i < 8; i++)
	{
		((c >> i) & 0x01) ? ones++ : zeroes++;
	}
}

void DNABuilder::get_uchar_binary(const uchar c, vector<uchar> & binary)
{
	for (int i = 7; i >= 0; i--)
	{
		binary.push_back((c & (1 << i)) ? '1' : '0');
	}
}

void DNABuilder::get_uint_binary(const uint a, vector<uchar> & binary)
{
	for (int i = 7; i >= 0; i--)
	{
		binary.push_back((a & (1 << i)) ? '1' : '0');
	}
}
