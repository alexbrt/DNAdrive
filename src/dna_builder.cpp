#include "dna_builder.h"

using namespace std;

DNABuilder::DNABuilder() :
	input(nullptr),
	options(nullptr) {}

DNABuilder::DNABuilder(const vector<uchar> * file, const SequenceOptions * options) :
	input(file),
	options(options) {}

void DNABuilder::set_input(const vector<uchar> * file)
{
	input = file;
}

void DNABuilder::set_options(const SequenceOptions * opt)
{
	options = opt;
}

vector<string> DNABuilder::construct_sequences()
{
	uint input_length = input->size(); // Input length in bytes
	uint content_length = options->get_content_length(); // Content length in bits
	vector<uchar> to_be_encoded; // Buffer containing binary digits to be encoded
	vector<uchar> bit_buffer; // Buffer which transfers binary digits from input file to to_be_encoded in order for the encoder to build a sequence
	vector<string> sequences;
	uint i = 0; // Byte index
	uint sequence_index = 0;
	while (i < input_length)
	{
		get_uchar_binary(input->at(i), bit_buffer); // Get the binary digits of the current byte
		i++;
		if (bit_buffer.size() >= content_length) // Check if there are enough binary digits to construct a sequence
		{
			to_be_encoded.insert(to_be_encoded.end(), bit_buffer.begin(), bit_buffer.begin() + content_length); // Transfer <content_length> bits
			bit_buffer.erase(bit_buffer.begin(), bit_buffer.begin() + content_length); // Erase bits <content_length> bits from bit_buffer
			string encoded;
			construct_sequence(to_be_encoded, sequence_index, encoded);
			sequences.push_back(encoded);
			to_be_encoded.clear();
			sequence_index++;
		}
	}
	if (bit_buffer.size() > 0) // If the buffer is not empty (it means our input length in bits was not a multiple of <content_length>
	{
		to_be_encoded.insert(to_be_encoded.end(), bit_buffer.begin(), bit_buffer.begin() + bit_buffer.size()); // Transfer remaining bits
		bit_buffer.clear();
		string encoded;
		construct_sequence(to_be_encoded, sequence_index, encoded); // Padding is automatically added if bit_buffer.size() < content_length
		sequences.push_back(encoded);
		to_be_encoded.clear();
		sequence_index++;
	}
	return sequences;
}

void DNABuilder::reset_sequence()
{
	nucleotide_to_digit.clear();
	zero_digit_to_nucleotide.clear();
	one_digit_to_nucleotide.clear();
}

void DNABuilder::construct_sequence(vector<uchar> input, uint sequence_index, string & sequence)
{
	reset_sequence();
	uint zero_index = 0, one_index = 0;
	// 1. Analyze current buffer (count the number of binary 0s and 1s in the content)
	uint buff_number_of_zeroes = 0, buff_number_of_ones = 0;
	for (uchar c : input)
	{
		(c == '0') ? buff_number_of_zeroes++ : buff_number_of_ones++;
	}
	/**
	*	- Make nucleotide map -
	*	In order to keep GC% content >= 50%, assign GC to
	*	represent the binary digits which appear the most.
	*/
	uchar marker_bit;
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
		marker_bit = 'G'; // Not definite - can be modified
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
		marker_bit = 'A'; // Not definite - can be modified
	}

	// 2. Analyze PAM (update 0's and 1's indexes)
	string pam = options->get_pam_invariant();
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
			zero_index++;
		}
		else
		{
			one_index++;
		}
	}
	sequence += pam;

	// 3. Add header
	string header;

	// 3.1. Get binary digits of sequence index
	vector<uchar> binary_seq_index;
	get_uint_binary(sequence_index, binary_seq_index);

	// 3.2. Remove redundancy (only keep <sequence_index_length> bits)
	binary_seq_index.erase(binary_seq_index.begin(), binary_seq_index.begin() + binary_seq_index.size() - options->get_sequence_index_length());
	uint seq_index_zeroes = 0, seq_index_ones = 0;

	// 3.3. Count binary 0s and 1s
	for (uchar c : binary_seq_index)
	{
		(c == '0') ? seq_index_zeroes++ : seq_index_ones++;
	}

	// 3.4. Add parity bit
	uchar parity_bit;
	if (nucleotide_to_digit['G'] == 0)
	{
		if ((buff_number_of_zeroes + seq_index_zeroes) % 2 == 0)
		{
			parity_bit = zero_digit_to_nucleotide[zero_index % 2];
			zero_index++;
		}
		else
		{
			parity_bit = one_digit_to_nucleotide[one_index % 2];
			one_index++;
		}
	}
	else
	{
		if ((buff_number_of_ones + seq_index_ones) % 2 == 0)
		{
			parity_bit = zero_digit_to_nucleotide[zero_index % 2];
			zero_index++;
		}
		else
		{
			parity_bit = one_digit_to_nucleotide[one_index % 2];
			one_index++;
		}
	}
	header += parity_bit;

	// 3.5. Add sequence index
	string dna_seq_index;
	for (uchar digit : binary_seq_index)
	{
		if (digit == '0')
		{
			dna_seq_index += zero_digit_to_nucleotide[zero_index % 2];
			zero_index++;
		}
		else
		{
			dna_seq_index += one_digit_to_nucleotide[one_index % 2];
			one_index++;
		}
	}
	header += dna_seq_index;
	sequence += header;

	// 4. Add content
	for (uchar c : input)	// Translate bits to nucleotides
	{
		if (c == '0')
		{
			sequence += zero_digit_to_nucleotide[zero_index % 2];
			zero_index++;
		}
		else
		{
			sequence += one_digit_to_nucleotide[one_index % 2];
			one_index++;
		}
	}

	// 5. Check if padding is necessary
	if (input.size() < options->get_content_length())
	{
		/**
		* Re-add the last nucleotide as a repeat for padding detection.
		* Note: This algorithm is designed not to permit any mononucleotide repeats.
		* Using this method, there's no need to specify input length within the sequence.
		*/
		uchar repeat = sequence[sequence.size() - 1];
		sequence += repeat;
		uint padding_index = 0;
		for (uint i = input.size() + 1; i < options->get_content_length(); i++)
		{
			// Add further GC-rich padding
			if (nucleotide_to_digit['G'] == 0)
			{
				sequence += zero_digit_to_nucleotide[(zero_index + padding_index) % 2];
			}
			else
			{
				sequence += one_digit_to_nucleotide[(one_index + padding_index) % 2];
			}
			padding_index++;
		}
	}

	// 6. Add marker bit
	/**
	* In case the last nucleotide of the content is the same as the bit marker,
	* in order to avoid 5'-3' internal PAM, A gets swapped with T and G gets swapped with C.
	* Note: Canonical PAM is guaranteed to contain mononucleotide repeats.
	*/
	if (sequence[sequence.size() - 1] == marker_bit)
	{
		if (nucleotide_to_digit['G'] == 0)
		{
			marker_bit = 'T';
		}
		else
		{
			marker_bit = 'C';
		}
	}
	sequence += marker_bit;

	// 7. Add trailing invariant
	sequence += options->get_trailing_invariant();
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
