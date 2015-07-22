#ifndef FFEA_INPUT_READER
#define FFEA_INPUT_READER

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"

using namespace std;

class FFEA_input_reader {
	
	public:
		FFEA_input_reader();
		~FFEA_input_reader();

		/* Set of functions to read FFEA files (pseudo-xml)*/

		// Get all lines from ffea, strip them of whitespace and return as a vector object
		int file_to_lines(string script_fname, vector<string> *output);

		// Extract any block from the current block
		int extract_block(string block_title, int block_index, vector<string> input, vector<string> *output);

		// Get rvalue from block
		int parse_tag(string input, string *output);

		// Specifically return map data
		int parse_map_tag(string input, int *map_indices, string *map_fname);

		// Split string around delim and return as strings
		int split_string(string input, string *output, string delim);

		// Split string around delim and return as ints
		int split_string(string input, int *output, string delim);

		// Split string around delim and return as scalars
		int split_string(string input, scalar *output, string delim);

	private:

		int max_buf_size;
		char *buf;
		string buf_string;
		int copying;
		vector<string>::iterator string_it;
		
};

#endif
