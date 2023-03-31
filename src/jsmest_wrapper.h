#ifndef JSMEST_WRAPPER
#define JSMEST_WRAPPER

#include<math.h>
#include<iostream>
#include<cstdlib>
#include<vector>

using namespace std;

typedef double real;

extern "C"{
  real stemp2_(const real *p);                // Saturation temperature [K]
  real htp2_(const real *t,const real *p);    // Specific enthalpy [J/kg]
  real vtp2_(const real *t,const real *p);    // Specific volume [m3/kg]
  real dhtptp2_(const real *t,const real *p);
  real dhpttp2_(const real *t,const real *p);
  real dvtp2_(const real *t,const real *p);
  //real kvtp2_(const real *t,const real *p);  
};


#endif
