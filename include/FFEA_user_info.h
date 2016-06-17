#ifndef FFEA_USER_INFO_H_INCLUDED
#define FFEA_USER_INFO_H_INCLUDED

#include <iostream>
#include <string>

using namespace std;

// Variables
namespace userInfo {
	extern string log_out_fname;
	extern int verblevel;
	extern FILE *log_out;
}

// Methods
void set_verbosity_level(int l);
void set_log_fname(string s);

void print_normal(string s);
void print_mid(string s);
void print_high(string s);
void print_mania(string s);


#endif
