#pragma once

#include<string>
#include<algorithm>

unsigned char get_nucleotide_complement(unsigned char c);

std::string get_complement(std::string seq);
std::string get_reverse_complement(std::string seq);
