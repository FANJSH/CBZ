#ifndef GEOMVECTOR
#define GEOMVECTOR

#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "Numeric.h"

class GeomVector{
  real x,y;
 public:
  GeomVector(){};
  GeomVector(real a, real b=0.){Initialize(a,b);};
  void Initialize(real a,real b);
  GeomVector operator+(const GeomVector &second);
  GeomVector operator-(const GeomVector &second);
  GeomVector operator*(real a);
  real getx(){return x;};
  real gety(){return y;};
  real CalcDistance(GeomVector second);
  void show_self();
  void put_data(real a=0,real b=0){x=a; y=b;};
  void unity(){NormalizedToUnity();};
  void NormalizedToUnity();
  int SameOrNot(GeomVector sec,real Judge=0.00001);
  bool BeIdentical(GeomVector sec,real Judge=1e-5);
  GeomVector OperatingElementaryReflector(GeomVector u);
  GeomVector OperatingRotator(real rot);
  void Shift(real dx, real dy);
};

#endif
