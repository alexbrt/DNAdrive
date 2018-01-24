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

string get_complement(string seq)
{
	string complement;
	for (auto c : seq)
	{
		complement.push_back(get_nucleotide_complement(c));
	}
	return complement;
}

string get_reverse_complement(string seq)
{
	string complement = get_complement(seq);
	reverse(complement.begin(), complement.end());
	return complement;
}
