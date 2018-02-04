#pragma once

#include <string>

class SequenceOptions
{
	typedef unsigned int uint;
	typedef unsigned int uchar;

	uint packet_index_length; // length of packet index in bits
	uint content_length; // length of content (packet length - header length - trailing length)
	std::string pam_invariant;
	std::string trailing_invariant;
public:
	SequenceOptions();
	SequenceOptions(uint packet_index_length, uint content_length, std::string pam_invariant, std::string trailing_invariant);
	uint get_sequence_length() const;
	uint get_header_length() const;
	void set_sequence_index_length(uint len);
	uint get_sequence_index_length() const;
	void set_content_length(uint len);
	uint get_content_length() const;
	void set_pam_invariant(std::string pam);
	std::string get_pam_invariant() const;
	void set_trailing_invariant(std::string trailing);
	std::string get_trailing_invariant() const;
};
