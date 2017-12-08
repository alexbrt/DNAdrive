#include<iostream>
#include<vector>
#include<fstream>
#include<iterator>
#include"dna_builder.h"

using namespace std;

int main()
{
	ifstream infile("test2.txt", ios::binary);
	vector<unsigned char> content((istreambuf_iterator<char>(infile)),
									istreambuf_iterator<char>());

	SequenceOptions options(4, 24, "AAG", "GA");
	DNABuilder builder(&content, &options);
	vector<string> protospacers = builder.construct_sequences();
	builder.reset_builder();

	for (string & p : protospacers)
	{
		cout << p << endl;
	}

	return 0;
}
