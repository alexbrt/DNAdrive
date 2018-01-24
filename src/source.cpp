#include<iostream>
#include<vector>
#include<fstream>
#include<iterator>
#include<random> // Testing
#include<algorithm> // Testing
#include<assert.h> // Testing
#include"dna_builder.h"
#include"dna_decoder.h"

using namespace std;

void randomize(vector<unsigned char> & in, unsigned int min_length, int max_length)
{
	random_device rd_1;
	mt19937 mt_1(rd_1());
	uniform_real_distribution<double> dist_1(min_length, nextafter(max_length, DBL_MAX));
	unsigned int rand_len = dist_1(mt_1);
	for (unsigned int i = 0; i < rand_len; i++)
	{
		random_device rd_2;
		mt19937 mt_2(rd_2());
		uniform_real_distribution<double> dist_2(0, nextafter(255, DBL_MAX));
		unsigned int rand_value = dist_2(mt_2);
		in.push_back(rand_value);
	}
}

int main()
{
	SequenceOptions options(4, 24, "AAG", "GA");

	DNABuilder builder;
	DNADecoder decoder;

	// "hello."
	ifstream hello_file("hello.txt", ios::binary);
	vector<unsigned char> hello_stream((istreambuf_iterator<char>(hello_file)),
		istreambuf_iterator<char>());
	builder.set_input(&hello_stream);
	builder.set_options(&options);
	vector<string> protospacers_hello = builder.construct_sequences();
	cout << "Protospacer Set Hello:" << endl;
	for (string & p : protospacers_hello)
	{
		cout << p << endl;
	}
	decoder.set_input(&protospacers_hello);
	decoder.set_options(&options);
	vector<unsigned char> decoded_hello = decoder.decode_sequences();
	cout << "Decoded Set Hello:" << endl;
	for (char c : decoded_hello)
	{
		cout << c;
	}
	cout << endl << endl;

	/*unsigned int iterations = 1000000;
	unsigned int i = 0;
	while (i < iterations)
	{
		vector<unsigned char> random;
		randomize(random, 5, 16 * 3);

		builder.set_input(&random);
		builder.set_options(&options);

		vector<string> seqs = builder.construct_sequences();
		random_shuffle(seqs.begin(), seqs.end());

		decoder.set_input(&seqs);
		decoder.set_options(&options);

		vector<unsigned char> decoded = decoder.decode_sequences();

		assert(decoded == random);
		i++;
	}*/
	return 0;
}
