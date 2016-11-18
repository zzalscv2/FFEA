#include <iostream>	
#include <sstream>
#include "FFEA_input_reader.h"
#include "SimulationParams.h"
#include <boost/filesystem.hpp>

int main (void)
{
	boost::filesystem::path path( boost::filesystem::current_path() );
	string file_to_open = path.string();
	file_to_open += string("/2r5u_8ang.ffea");
	// const int MAX_BUF_SIZE = 255;
	string buf_string;
	SimulationParams params;
	FFEA_input_reader *ffeareader;
	ffeareader = new FFEA_input_reader();

	// Copy entire script into string
	vector<string> script_vector;
	if(ffeareader->file_to_lines(file_to_open, &script_vector) == FFEA_ERROR) {
		return(1);
	}

	// Get params section
	cout << "Extracting Parameters..." << endl;
	params.FFEA_script_filename = file_to_open;  // includes absolute path.
	if(params.extract_params(script_vector) != 0) {
		FFEA_error_text();
		printf("Error parsing parameters in SimulationParams::extract_params()\n");
		return(1);
	}
	cout << "...done!" << endl;

	// Check for consistency
	cout << "\nVerifying Parameters..." << endl;
	if(params.validate() != 0) {
		FFEA_error_text();
		printf("Parameters found to be inconsistent in SimulationParams::validate()\n");
		return(1);
	}
	
	return(0);
}
