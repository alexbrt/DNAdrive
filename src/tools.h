#pragma once

#include <string>
#include <algorithm>

unsigned char get_nucleotide_complement(unsigned char c);

std::string get_complement(const std::string & seq);
std::string get_reverse_complement(const std::string & seq);

std::string protospacer_to_hairpin(const std::string & seq, unsigned int overhang_length, unsigned int loop_length);
