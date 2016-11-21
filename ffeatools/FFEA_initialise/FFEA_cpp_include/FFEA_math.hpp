#ifndef FFEA_MATH_INCLUDED
#define FFEA_MATH_INCLUDED

#include <iostream>
using namespace std;

class vector3i {

	public:

		// Constructor / Destructor
		vector3i();
		~vector3i();

		// Overload Operators
		vector3i operator+(vector3i rvalue);
		vector3i operator+=(vector3i rvalue);
		vector3i operator-(vector3i rvalue);
		vector3i operator-=(vector3i rvalue);
		vector3i operator*(float rvalue);
		vector3i operator*=(float rvalue);

		// Data access functions
		void set_components(int a, int b, int c);
		float get_magnitude();

		// Operation functions
		void normalise();

		// Member variables
		int comp[3];
		
	private:

		float magnitude;
};

class vector3f {

	public:

		// Constructor / Destructor
		vector3f();
		~vector3f();

		// Overload Operators
		vector3f operator+(vector3f rvalue);
		vector3f operator+=(vector3f rvalue);
		vector3f operator-(vector3f rvalue);
		vector3f operator-=(vector3f rvalue);
		vector3f operator*(float rvalue);
		vector3f operator*=(float rvalue);

		// Data access functions
		void set_components(float a, float b, float c);
		float get_magnitude();

		// Operation functions
		void normalise();

		// Member variables
		float comp[3];
		
	private:

		float magnitude;
};
#endif
