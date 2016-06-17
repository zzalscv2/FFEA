#include "FFEA_user_info.h"

using namespace std;

string userInfo::log_out_fname;
int userInfo::verblevel;
FILE * userInfo::log_out;

void set_verbosity_level(int l) {
	userInfo::verblevel = (l < 4) ? l : 3;
}

void set_log_fname(string s) {
	userInfo::log_out_fname = s;
}

void print_mid(string s) {
	if(userInfo::verblevel > 0) {
		cout << s << endl;
	}
}

void print_high(string s) {
	if(userInfo::verblevel > 1) {
		cout << s << endl;
	}
}

void print_mania(string s) {
	if(userInfo::verblevel > 2) {
		cout << s << endl;
	}
}
