#ifndef FFEA_USER_INFO_H_INCLUDED
#define FFEA_USER_INFO_H_INCLUDED

#include <iostream>
#include <string>

using namespace std;
static int verblevel;

void set_verbosity_level(int l);
void print_normal(string s);
void print_mid(string s);
void print_high(string s);
void print_mania(string s);

#endif