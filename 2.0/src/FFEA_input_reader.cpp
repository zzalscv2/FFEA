#include "FFEA_input_reader.h"

FFEA_input_reader::FFEA_input_reader() {
	max_buf_size = 255;
	buf = new char[max_buf_size];
	buf_string = "";
	copying = 0;
}

FFEA_input_reader::~FFEA_input_reader() {
	max_buf_size = 0;
	delete[] buf;
	buf = NULL;
	buf_string = "";
	copying = 0;
}

int FFEA_input_reader::file_to_lines(string script_fname, vector<string> *script_vector) {

	// Open script
	ifstream fin;
	fin.open(script_fname.c_str());
	if(fin.fail()) {
		FFEA_error_text();
		cout << script_fname << " does not exist.\n" << endl;
		return FFEA_ERROR;
	}

	// Copy entire script into string
	script_vector->clear();

	while(!fin.eof()) {
		fin.getline(buf, max_buf_size);
		buf_string = string(buf);
		boost::trim(buf_string);
		if(buf_string != "\n" && buf_string != "") {
			script_vector->push_back(buf_string);
		}
	}
	fin.close();
	return FFEA_OK;
}

int FFEA_input_reader::extract_block(string block_title, int block_index, vector<string> input, vector<string> *output) {

	// Immediate error checking
	if(block_title != "param" && block_title != "system" && block_title != "blob" && block_title != "conformation" && block_title != "kinetics" && block_title != "maps" && block_title != "interactions" && block_title != "springs") {
		FFEA_error_text();
		cout << "Unrecognised block. Block structure is:" << endl;
		cout << "<param>\n</param>\n<system>\n\t<blob>\n\t\t<conformation>\n\t\t</conformation>\n\t\t\t.\n\t\t\t.\n\t\t\t.\n\t\t<kinetics>\n\t\t\t<maps>\n\t\t\t</maps>\n\t\t</kinetics>\n\t</blob>\n\t\t.\n\t\t.\n\t\t.\n\t<interactions>\n\t\t<springs>\n\t\t</springs>\n\t</interactions>\n</system>" << endl;
		return FFEA_ERROR;
	}

	int count = -1;
	output->clear();
	for(string_it = input.begin(); string_it != input.end(); ++string_it) {

		buf_string = boost::erase_last_copy(boost::erase_first_copy(*string_it, "<"), ">");
		boost::trim(buf_string);
		if(buf_string == block_title) {
			
			if(copying == 1) {
				FFEA_error_text();
				cout << "Shouldn't have found '" << buf_string << "' within " << block_title << " block." << endl;
				return FFEA_ERROR;
			}

			count++;
			if(count == block_index) {
				copying = 1;
				continue;
			}
		}

		if(copying == 1) {
			if(buf_string == "/" + block_title) {
				copying = 0;
				return FFEA_OK;
			} else {
				output->push_back(*string_it);
			}
		}
	}

	FFEA_error_text();
	if(copying == 1) {
		cout << "Never found closing tag '/" << block_title << "'." << endl;
		return FFEA_ERROR;
	} else {
		cout << "Specified block_index " << block_index << " for block '" << block_title << "' not found." << endl;
		return FFEA_OK;
	}
	
}

int FFEA_input_reader::parse_tag(string input, string *output) {

	string tag;
	vector<string> lrvalvec;
	tag = boost::erase_last_copy(boost::erase_first_copy(input, "<"), ">");
	boost::trim(tag);

	// Split around "=", trim and return
	boost::split(lrvalvec, tag, boost::is_any_of("="));
	for(int i = 0; i < lrvalvec.size(); ++i) {
		output[i] = lrvalvec.at(i);
		boost::trim(output[i]);
	}

	return FFEA_OK;
}

int FFEA_input_reader::parse_map_tag(string input, int *map_indices, string *map_fname) {

	// Parse whole tag
	string lrvalue[2], indexlrvalue[2];
	parse_tag(input, lrvalue);

	// Check if map
	split_string(lrvalue[0], indexlrvalue, "(");
	if(indexlrvalue[0] != "map") {
		cout << indexlrvalue[0] << endl;
		FFEA_error_text();
		cout << "Expected '<map (from,to) = fname>' but got " << input << endl;
		return FFEA_ERROR;
	}

	// Assign map!
	*map_fname = lrvalue[1];

	// Get indices
	indexlrvalue[1] = boost::erase_last_copy(indexlrvalue[1], ")");
	split_string(indexlrvalue[1], map_indices, ",");
	return FFEA_OK;
}

int FFEA_input_reader::split_string(string input, string *output, string delim) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		output[i] = *it;
		boost::trim(output[i++]);
	}
	return FFEA_OK;
}

int FFEA_input_reader::split_string(string input, int *output, string delim) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		output[i++] = atoi((*it).c_str());
	}
	return FFEA_OK;
}

int FFEA_input_reader::split_string(string input, scalar *output, string delim) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		output[i++] = atof((*it).c_str());
	}
	return FFEA_OK;
}