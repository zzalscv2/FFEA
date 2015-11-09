#ifndef DIMENSIONS_H_INCLUDED
#define DIMENSIONS_H_INCLUDED


typedef struct {
   double length;
   double Energy;
   double time;
   double pressure;
   double Temperature;
   double velocity;
   double mass;
   double permittivity;
   double charge;
} dimset;

class Dimensions{
  public:
    dimset atomic {
          5.2917721092e-11,
          4.35974417e-18,
          2.418884326505e-17,
          2.9421912e13,
          3.1577464e5,
          2.1876912633e6,
          9.10938291e-31,
          8.8541878176204667e-12,
          1.602176565e-19
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
          1
    };
    
    
};

#endif
