#ifndef MAT_VEC_TYPES_H_INCLUDED
#define MAT_VEC_TYPES_H_INCLUDED

#include <cstddef>
#include <limits>
#include <cmath> 

/*
 * Defines what is meant by a scalar (essentially sets the precision of
 * the code between float or double).
 */
#ifdef USE_DOUBLE
typedef double scalar;
typedef double geoscalar;
#else
typedef float scalar;
typedef double geoscalar;
#endif 
//typedef long double scalar;

////////  Constants and scalar functions ///////
namespace ffea_const {
   const scalar threeErr = 3.0*std::numeric_limits<scalar>::epsilon();
   const scalar threeGeoErr = 3.0*std::numeric_limits<geoscalar>::epsilon();
   const scalar mOne = -1.0;
   const scalar zero = 0.0;
   const scalar half = 0.5;
   const scalar one = 1.0;
   const scalar two = 2.0;
   const scalar eight = 8.0;
   const geoscalar ten = 10.0;
   const scalar twentyfour = 24.0;
   const scalar oneOverThree = 0.33333333333333333;
   const scalar oneOverSix = 0.166666666666666667;
   const scalar oneOverEight = 0.12500000000000000; 
   const scalar sphereFactor = 4.0 * 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148 / 3.0;
}

/**
 * A simple 3 dimensional vector (x, y, z)
 */
typedef struct {
    scalar x, y, z;
} vector3;

/** arr3 will hopefully substitute vector3 one day */ 
typedef scalar arr3[3]; 
typedef geoscalar grr3[3]; 

/*  // The following class will be a temporary replacement
    //   for the current vector3 struct so that we have some
    //   time to change from Whatever.x into Whatever[0] 
    // Then vector3 will become arr3.
    // Therefore, don't ever populate this class with methods!!
class vector3 {
public:
    scalar data[3]; 
    scalar& x = data[0]; 
    scalar& y = data[0]; 
    scalar& z = data[0]; 
    scalar& operator [](std::size_t i) { return data[i]; }
};*/

/** arr3_view has the spirit of "span" or "array_view" */
template <class t_scalar, class brr3> class arr3_view
{
public:
    arr3_view(brr3 (&arr) ) : data(arr) {}
    arr3_view(t_scalar* data, std::size_t size) : data(data) {
       // if (size != 3) // simple check removed in benefit of performance.
       //  throw std::runtime_error("arr3 - wrong size of data: " + std::to_string(size)); 
    }

    t_scalar* begin() { return data; }
    t_scalar* end() { return data + 3; }
    t_scalar& operator [](std::size_t i) { return data[i]; }
    // and similar stuff as above for const versions

private:
    t_scalar* data;
}; 



/**
 * Defines a 12 vector (just for ease of use, clarity and compiler type checking)
 */
typedef scalar vector12[12];

/**
 * Defines a 12x12 matrix
 */
typedef scalar matrix12[12][12];

/**
 * Defines a 3x3 matrix
 */
typedef scalar matrix3[3][3];

/**
 * Defines a 4x4 matrix
 */
typedef scalar matrix4[4][4];

/**
 * A useful type for holding the upper triangular part of symmetric 4x4 matrices
 */
typedef struct {
    scalar u00, u01, u02, u03,
    u11, u12, u13,
    u22, u23,
    u33;
} upper_triangular_matrix4;
#endif
