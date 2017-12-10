#pragma once

#include<vector>
#include<deque>
#include<string>
#include<map>
#include<cmath>
#include<iostream>
#include<boost/dynamic_bitset.hpp>
#include"sequence_options.h"

class DNADecoder
{
	typedef unsigned int uint;
	typedef unsigned char uchar;

	const std::vector<std::string> * sequences;
	const SequenceOptions * options;

	uint number_of_sequences;
	uint number_of_output_bytes;
	std::map<uchar, uint> nucleotide_to_digit;
public:
	DNADecoder();
	DNADecoder(const std::vector<std::string> * sequences, const SequenceOptions * options);
	void set_input(const std::vector<std::string> * seq);
	void set_options(const SequenceOptions * opt);
	std::vector<uchar> decode_sequences();
private:
	void decode_sequence(std::string sequence, uint sequence_number, uint & sequence_index, std::string & decoded);
	void get_uchar_from_binary(std::string & bits, uchar & out);
	void get_uint_from_binary(std::string & bits, uint & out);
};
