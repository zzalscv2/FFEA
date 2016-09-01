#include "FFEA_input_reader.h"

int main (void)
{
	FFEA_input_reader *ffeareader;
	ffeareader = new FFEA_input_reader();
	cout << "if you're seeing this it built correctly" << endl;
	vector<string> script_vector;

	if(ffeareader->file_to_lines("2r5u_8ang.ffea", &script_vector) == FFEA_ERROR) {
		return FFEA_ERROR;
	}
	cout << "if you're seeing this it read the filename correctly" << endl;
    return(0);
}
