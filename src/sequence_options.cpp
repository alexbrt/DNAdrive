#include "sequence_options.h"

using namespace std;

namespace // Default strings and values for use in Escherichia Coli B strains (such as BL21-AI or BL21-DE3)
{
	unsigned int default_packet_index_length = 4;
	unsigned int default_content_length = 24;
	string default_pam_invariant = "AAG";
	string default_trailing_invariant = "GA";
}

SequenceOptions::SequenceOptions() :
	packet_index_length(default_packet_index_length),
	content_length(default_content_length),
	pam_invariant(default_pam_invariant),
	trailing_invariant(default_trailing_invariant) {}

SequenceOptions::SequenceOptions(uint packet_index_length, uint content_length, string pam_invariant, string trailing_invariant) :
	packet_index_length(packet_index_length),
	content_length(content_length),
	pam_invariant(pam_invariant),
	trailing_invariant(trailing_invariant) {}

SequenceOptions::uint SequenceOptions::get_sequence_length() const
{
	return pam_invariant.size() + (packet_index_length + 1) + content_length + (trailing_invariant.size() + 1);
}

SequenceOptions::uint SequenceOptions::get_header_length() const
{
	return (packet_index_length + 1);
}

void SequenceOptions::set_sequence_index_length(uint len)
{
	packet_index_length = len;
}

SequenceOptions::uint SequenceOptions::get_sequence_index_length() const
{
	return packet_index_length;
}

void SequenceOptions::set_content_length(uint len)
{
	content_length = len;
}

SequenceOptions::uint SequenceOptions::get_content_length() const
{
	return content_length;
}

void SequenceOptions::set_pam_invariant(string pam)
{
	pam_invariant = pam;
}

string SequenceOptions::get_pam_invariant() const
{
	return pam_invariant;
}

void SequenceOptions::set_trailing_invariant(string trailing)
{
	trailing_invariant = trailing;
}

string SequenceOptions::get_trailing_invariant() const
{
	return trailing_invariant;
}
