#include "tools.h"

using namespace std;

unsigned char get_nucleotide_complement(unsigned char c)
{
	switch (c)
	{
	case 'A': return 'T';
	case 'T': return 'A';
	case 'G': return 'C';
	case 'C': return 'G';
	default: return '\0';
	}
}

string get_complement(const string & seq)
{
	string complement;
	for (auto c : seq)
	{
		complement.push_back(get_nucleotide_complement(c));
	}
	return complement;
}

string get_reverse_complement(const string & seq)
{
	string complement = get_complement(seq);
	reverse(complement.begin(), complement.end());
	return complement;
}

string protospacer_to_hairpin(const string & seq, unsigned int overhang_length, unsigned int loop_length)
{
	string out;
	string copy_1 = seq;
	string copy_2 = seq;
	copy_1.erase(0, overhang_length);
	copy_2.erase(copy_2.size() - loop_length, loop_length);
	reverse(copy_2.begin(), copy_2.end());
	for (unsigned int i = 0; i < copy_2.size(); i++)
	{
		switch (copy_2[i])
		{
		case 'A':
			copy_2[i] = 'T';
			break;
		case 'T':
			copy_2[i] = 'A';
			break;
		case 'G':
			copy_2[i] = 'C';
			break;
		case 'C':
			copy_2[i] = 'G';
			break;
		}
	}
	return copy_1 + copy_2;
}
