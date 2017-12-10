#include "dna_decoder.h"

using namespace std;
using namespace boost;

DNADecoder::DNADecoder() :
	sequences(nullptr),
	options(nullptr) {}

DNADecoder::DNADecoder(const std::vector<std::string> * sequences, const SequenceOptions * options) :
	number_of_sequences(sequences->size()),
	number_of_output_bytes(sequences->size() * options->get_content_length() / 8),
	sequences(sequences),
	options(options) {}

void DNADecoder::set_input(const vector<string>* seq)
{
	sequences = seq;
}

void DNADecoder::set_options(const SequenceOptions * opt)
{
	options = opt;
}

vector<DNADecoder::uchar> DNADecoder::decode_sequences()
{
	uint length = options->get_content_length() * sequences->size();
	dynamic_bitset<uchar> bit_buffer(length); // All 0's by default
	uint real_length = 0;
	for (uint i = 0; i < sequences->size(); i++)
	{
		uint sequence_index = 0;
		string decoded_sequence;
		decode_sequence(sequences->at(i), i, sequence_index, decoded_sequence);
		uint j = 0;
		for (j = 0; j < decoded_sequence.size(); j++)
		{
			if (decoded_sequence[j] == '1')
			{
				bit_buffer[length - 1 - (sequence_index * options->get_content_length() + j)] = 1; // Append the bits in reverse order
			}
		}
		if (j < options->get_content_length()) // Padding was detected
		{
			real_length = (length - (options->get_content_length() - j)) / 8; // Compute real content length in bytes (without padding)
		}
	}
	deque<uchar> decoded;
	to_block_range(bit_buffer, front_inserter(decoded));
	while (real_length > 0 && decoded.size() > real_length)
	{
		decoded.pop_back(); // Pop padded bytes
	}
	return vector<uchar>(decoded.begin(), decoded.end()); // Return vector
}

void DNADecoder::decode_sequence(string sequence, uint sequence_number, uint & sequence_index, string & decoded)
{
	// Remove invariants
	sequence.erase(0, options->get_pam_invariant().size()); // Remove PAM
	sequence.erase(sequence.size() - options->get_trailing_invariant().size(), options->get_trailing_invariant().size()); // Remove trailing invariant

	// Make nucleotide map
	uchar marker_bit = sequence[sequence.size() - 1];
	if (marker_bit == 'A' || marker_bit == 'T')
	{
		nucleotide_to_digit['A'] = 1;
		nucleotide_to_digit['T'] = 1;
		nucleotide_to_digit['G'] = 0;
		nucleotide_to_digit['C'] = 0;
	}
	else
	{
		nucleotide_to_digit['A'] = 0;
		nucleotide_to_digit['T'] = 0;
		nucleotide_to_digit['G'] = 1;
		nucleotide_to_digit['C'] = 1;
	}
	sequence.erase(sequence.end() - 1); // Remove marker bit

	// Decode header

	// Step 1 - Check for data corruption
	uchar parity_bit = sequence[0];
	uint GC_count = 0;
	uchar prev_c = 'X';
	uint padding_begin = options->get_content_length(); // Index of padding start | Default = no padding
	for (uint i = 1; i < sequence.size(); i++)
	{
		if (sequence[i] == prev_c) // If nucleotide repetition is detected within content, break and set padding_begin
		{
			padding_begin = i;
			break;
		}
		if (sequence[i] == 'G' || sequence[i] == 'C')
		{
			GC_count++;
		}
		prev_c = sequence[i];
	}
	if ((GC_count % 2) != nucleotide_to_digit[parity_bit])
	{
		cout << "Corrupted data detected in sequence #" << sequence_number << endl;
	}
	// Step 2 - Get sequence index
	string index_bits;
	for (uint i = 1; i < options->get_header_length(); i++)
	{
		index_bits += (nucleotide_to_digit[sequence[i]] + '0');
	}
	get_uint_from_binary(index_bits, sequence_index);
	sequence.erase(0, options->get_header_length()); // Remove header

	// Decode content
	string decoded_content;
	for (uint i = 0; i < padding_begin; i++)
	{
		decoded_content += (nucleotide_to_digit[sequence[i]] + '0'); // Get correponding bit in its char representation
	}
	decoded = decoded_content;
}

void DNADecoder::get_uchar_from_binary(string & bits, uchar & out)
{
	out = 0x00;
	for (uint i = 0, j = bits.size() - 1; i < bits.size(); i++, j--)
	{
		out |= ((bits[i] - '0') << j);
	}
}

void DNADecoder::get_uint_from_binary(std::string & bits, uint & out)
{
	out = 0;
	for (uint i = 0, j = bits.size() - 1; i < bits.size(); i++, j--)
	{
		out |= ((bits[i] - '0') << j);
	}
}
