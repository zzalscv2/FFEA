#include "FFEA_user_info.h"

using namespace std;

void set_verbosity_level(int l) {
	verblevel = (l < 4) ? l : 3;
}

void print_normal(string s) {
	cout << s << endl;
}

void print_mid(string s) {
	if(verblevel > 0) {
		cout << s << endl;
	}
}

void print_high(string s) {
	if(verblevel > 1) {
		cout << s << endl;
	}
}

void print_mania(string s) {
	if(verblevel > 2) {
		cout << s << endl;
	}
}