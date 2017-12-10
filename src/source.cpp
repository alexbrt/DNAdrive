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

	//// Test 1
	//ifstream infile1("test1.txt", ios::binary);
	//vector<unsigned char> content1((istreambuf_iterator<char>(infile1)),
	//								istreambuf_iterator<char>());
	//builder.set_input(&content1);
	//builder.set_options(&options);
	//vector<string> protospacers1 = builder.construct_sequences();
	//cout << "Protospacer Set 1:" << endl;
	//for (string & p : protospacers1)
	//{
	//	cout << p << endl;
	//}
	//decoder.set_input(&protospacers1);
	//decoder.set_options(&options);
	//vector<unsigned char> decoded_1 = decoder.decode_sequences();
	//cout << "Decoded Set 1:" << endl;
	//for (char c : decoded_1)
	//{
	//	cout << c;
	//}
	//cout << endl << endl;

	//// Test 2
	//ifstream infile2("test2.txt", ios::binary);
	//vector<unsigned char> content2((istreambuf_iterator<char>(infile2)),
	//	istreambuf_iterator<char>());
	//builder.set_input(&content2);
	//builder.set_options(&options);
	//vector<string> protospacers2 = builder.construct_sequences();
	//cout << "Protospacer Set 2:" << endl;
	//for (string & p : protospacers2)
	//{
	//	cout << p << endl;
	//}
	//decoder.set_input(&protospacers2);
	//decoder.set_options(&options);
	//vector<unsigned char> decoded_2 = decoder.decode_sequences();
	//cout << "Decoded Set 2:" << endl;
	//for (char c : decoded_2)
	//{
	//	cout << c;
	//}
	//cout << endl << endl;

	unsigned int iterations = 100000;
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
	}
	return 0;
}
