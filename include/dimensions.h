#ifndef DIMENSIONS_H_INCLUDED
#define DIMENSIONS_H_INCLUDED


typedef struct {
   const double length;
   const double area; 
   const double volume; 
   const double Energy;
   const double force;
   const double time;
   const double pressure;
   const double Temperature;
   const double velocity;
   const double mass;
   const double permittivity;
   const double charge;
} dimset;

class Dimensions{
  public:
    dimset atomic {
          5.2917721092e-11,
          2.8002852055707021e-21,
          1.4818471148644432e-31,
          4.35974417e-18,
          8.2387224544692213e-08,
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
          1,
          1,
          1,
          1
    };
    
    
};

#endif
