#ifndef DIMENSIONS_H_INCLUDED
#define DIMENSIONS_H_INCLUDED

#include "mat_vec_types.h"


namespace mesoDimensions {
   const scalar length = 1.7e-10; ///< length = C atom VdW radius
   const scalar Energy = 4.141945559999999e-21; ///< Energy = KbT, T=300K
   const scalar mass = 1.994307387553024e-26; ///< mass = C atom mass
   const scalar charge = 1.602176565e-19;  ///< charge = electron charge
   const scalar area = 2.89e-20; ///< area = l**2
   const scalar volume = 4.913e-30; ///< volume = l**3
   const scalar force = 2.4364705882352941e-11; ///< force = E/l
   const scalar time = 3.7302915560886126e-13; ///< time = l*sqrt(m/E)
   const scalar pressure = 8.4305832688784826e8; ///< pressure = E/l**3
   const scalar velocity = 4.5572845297447225e2; ///< velocity = sqrt(E/m)
}

/// atomic units
namespace atomicDimensions {
   const scalar length = 5.2917721092e-11; /* length */
   const scalar Energy = 4.35974417e-18; /* Energy */
   const scalar mass = 9.10938291e-31; /* mass */
   const scalar charge = 1.602176565e-19; /* charge */
   const scalar area = 2.8002852055707021e-21; /* area */ 
   const scalar volume = 1.4818471148644432e-31; /* volume */ 
   const scalar force = 8.2387224544692213e-08; /* force */
   const scalar time = 2.418884326505e-17; /* time */
   const scalar pressure = 2.9421912e13; /* pressure */
   const scalar velocity = 2.1876912633e6; /* velocity */
}

#endif
