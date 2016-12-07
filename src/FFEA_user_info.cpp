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
