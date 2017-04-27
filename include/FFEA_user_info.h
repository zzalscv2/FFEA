// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

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

void print_ffea_version();
void print_ffea_compilation_details();


#endif
