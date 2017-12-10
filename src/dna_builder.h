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
public:
	DNABuilder();
	DNABuilder(const std::vector<uchar> * file, const SequenceOptions * options);
	void set_input(const std::vector<uchar> * file);
	void set_options(const SequenceOptions * opt);
	std::vector<std::string> construct_sequences();
	void reset_sequence();
private:
	void construct_sequence(std::vector<uchar> input, uint sequence_index, std::string & sequence);
	void get_uchar_binary(const uchar c, std::vector<uchar> & binary);
	void get_uint_binary(const uint a, std::vector<uchar> & binary);
};
