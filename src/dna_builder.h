#pragma once

#include<string>
#include<utility>
#include<map>
#include<vector>
#include"sequence_options.h"

class DNABuilder
{
	typedef unsigned int uint;
	typedef unsigned char uchar;

	const std::vector<uchar> * input;
	const SequenceOptions * options;

	std::map<uchar, uint> nucleotide_to_digit;
	std::vector<uchar> zero_digit_to_nucleotide;
	std::vector<uchar> one_digit_to_nucleotide;
	std::vector<uchar> buffer; // Current sequence buffer (contains binary digits)
	uchar curr_marker_bit; // Indicates what G/C represent as binary digits (A or T means G/C => 0, G or C means G/C => 1)
	uint curr_sequence_index;
	uint buff_number_of_zeroes;
	uint buff_number_of_ones;
	uint curr_zero_index;
	uint curr_one_index;

	std::vector<std::string> dna_sequences;
public:
	DNABuilder();
	DNABuilder(const std::vector<uchar> * file, const SequenceOptions * options);
	std::string construct_sequence();
	std::vector<std::string> & construct_sequences();
	void reset_current_sequence();
	void reset_builder();
private:
	void make_nucleotide_map();
	void add_pam(std::string & sequence);
	void add_header(std::string & sequence);
	void add_content(std::string & sequence);
	void add_padding(std::string & sequence);
	void add_marker_bit(std::string & sequence);
	void count_binary_digits(const uchar c, uint & zeroes, uint & ones);
	void count_binary_digits(const uint c, uint & zeroes, uint & ones);
	void get_uchar_binary(const uchar c, std::vector<uchar> & binary);
	void get_uint_binary(const uint a, std::vector<uchar> & binary);
};
