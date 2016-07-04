#ifndef MATRICES_H
#define MATRICES_H

#include <cmath>
#include "Vectors.hpp"

using namespace std;

class matrix33
{
	public:

		// Constructors/Destructors
		matrix33() {
			for(int i = 0; i < 3; ++i) {
				for(int j = 0; j < 3; ++j) {
					val[i][j] = 0.0;
				}
			}
		}

		matrix33(vector3 a, vector3 b, vector3 c) {
			val[0][0] = a.x;
			val[1][0] = a.y;
			val[2][0] = a.z;
			val[0][1] = b.x;
			val[1][1] = b.y;
			val[2][1] = b.z;
			val[0][2] = c.x;
			val[1][2] = c.y;
			val[2][2] = c.z;
		}

		~matrix33() {
			for(int i = 0; i < 3; ++i) {
				for(int j = 0; j < 3; ++j) {
					val[i][j] = 0.0;
				}
			}
		}

		// Overloaded operators
		matrix33 operator*(double rvalue) {
			matrix33 result;
			for(int i = 0; i < 3; ++i) {
				for(int j = 0; j < 3; ++j) {
					result.val[i][j] = val[i][j] * rvalue;
				}
			}
			return result;
		}

		matrix33 operator*=(double rvalue) {
			for(int i = 0; i < 3; ++i) {
				for(int j = 0; j < 3; ++j) {
					val[i][j] *= rvalue;
				}
			}
		}

		// Member functions
		void set(vector3 a, vector3 b, vector3 c) {
			val[0][0] = a.x;
			val[1][0] = a.y;
			val[2][0] = a.z;
			val[0][1] = b.x;
			val[1][1] = b.y;
			val[2][1] = b.z;
			val[0][2] = c.x;
			val[1][2] = c.y;
			val[2][2] = c.z;
			return;
		}
	
		double get_determinant() {
			return val[0][0] * (val[1][1] * val[2][2] - val[2][1] * val[1][2]) + val[0][1] * (val[2][0] * val[1][2] - val[1][0] * val[2][2]) + val[0][2] * (val[1][0] * val[2][1] - val[2][0] * val[1][1]);
		}

		matrix33 get_inverse() {
			double det;
			matrix33 inv;

			det = get_determinant();

			inv.val[0][0] = val[1][1] * val[2][2] - val[2][1] * val[1][2];
			inv.val[1][0] = val[1][2] * val[2][0] - val[1][0] * val[2][2];
			inv.val[2][0] = val[1][0] * val[2][1] - val[1][1] * val[2][0];
			inv.val[0][1] = val[0][2] * val[2][1] - val[0][1] * val[2][2];
			inv.val[1][1] = val[0][0] * val[2][2] - val[0][2] * val[2][0];
			inv.val[2][1] = val[0][1] * val[2][0] - val[0][0] * val[2][1];
			inv.val[0][2] = val[0][1] * val[1][2] - val[0][2] * val[1][1];
			inv.val[1][2] = val[0][2] * val[1][0] - val[0][0] * val[1][2];
			inv.val[2][2] = val[0][0] * val[1][1] - val[0][1] * val[1][0];

			inv *= 1.0 / det;
			return inv;
		}

		matrix33 apply(matrix33 a) {
			int i, j, k;
			matrix33 b;

			for(i = 0; i < 3; ++i) {
				for(j = 0; j < 3; ++j) {
					for(k = 0; k < 3; ++k) {
						b.val[i][j] += val[i][k] * a.val[k][j];
					}	
				}
			}
			return b;
		}

		vector3 apply(vector3 a) {
			int i;
			vector3 b;

			b.x = val[0][0] * a.x + val[0][1] * a.y + val[0][2] * a.z;
			b.y = val[1][0] * a.x + val[1][1] * a.y + val[1][2] * a.z;
			b.z = val[2][0] * a.x + val[2][1] * a.y + val[2][2] * a.z;

			return b;
		}

		void print() {

			int i, j;
			cout << "Matrix:\n" << endl;
			for(i = 0; i < 3; ++i) {
				for(j = 0; j < 3; ++j) {
					cout << val[i][j] << " ";				
				}
				cout << endl;
			}

			return;
		}

		// Data members
		double val[3][3];
};

#endif
