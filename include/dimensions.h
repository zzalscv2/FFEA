#ifndef DIMENSIONS_H_INCLUDED
#define DIMENSIONS_H_INCLUDED

#include "mat_vec_types.h"


typedef struct {
   const scalar length;
   const scalar Energy;
   const scalar mass;
   const scalar charge;
   const scalar area; 
   const scalar volume; 
   const scalar force;
   const scalar time;
   const scalar pressure;
   const scalar velocity;
} dimset;

class Dimensions{
  public:
    dimset meso {
          1.7e-10, /* length = C atom VdW radius */
          4.142e-21, /* Energy = KbT */
          1.994307387553024e-26, /* mass = C atom mass */
          1.602176565e-19,  /* charge = electron charge */
          2.89e-20, /* area */ 
          4.913e-30, /* volume */ 
          2.4364705882352941e-11, /* force = E/l */
          3.7302670416342907e-13, /* time = sqrt(m*l/f) */
          8.4306940769387329e8, /* pressure */
          4.5573144791671604e2 /* velocity */
    };

    dimset atomic {
          5.2917721092e-11, /* length */ 
          4.35974417e-18, /* Energy */
          9.10938291e-31, /* mass */
          1.602176565e-19, /* charge */
          2.8002852055707021e-21, /* area */
          1.4818471148644432e-31, /* volume */
          8.2387224544692213e-08, /* force */
          2.418884326505e-17, /* time */
          2.9421912e13, /* pressure */
          2.1876912633e6 /* velocity */
    };
    
    dimset is {
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1
    };
    
    
};

#endif
