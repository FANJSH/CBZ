#ifndef GEOMPLANE
#define GEOMPLANE

#include <fstream>
#include "GeomVector.h"

class GeomPlane{
 protected:
  GeomVector p1,p2;
  GeomVector pver; // Vertical vector
  real a,b,c;     // (ax+by+c=0)
  real a1,b1,c1;  // Vertical line (a1x+b1y+c1=0)
  bool Exist;
 public:
  GeomPlane(){Exist=false;};
  GeomPlane(GeomVector i,GeomVector j){PutVector(i,j);};
  GeomPlane(real x1,real y1,real x2,real y2){PutVector(x1,y1,x2,y2);};
  void PutVector(GeomVector i,GeomVector j);
  void PutVector(real x1,real y1,real x2,real y2);
  void put_line();  // Calc. a,b,c,a1,b1,c1
  void show_self();
  real CalcDistance(GeomVector &p);
  real CalcDistance(real x,real y);
  GeomVector& get_p1(){return p1;};
  GeomVector& get_p2(){return p2;};
  GeomVector CrossPlane(GeomPlane &second);
  bool BeIdentical(GeomPlane sec,real judge=1e-5);
  bool BePartial(GeomPlane sec,real judge=1e-5);
  bool BeCrossing(GeomPlane sec,real judge=1e-5);
  bool VectorOnPlane(GeomVector p,real judge=1e-4);
  real get_xsi(GeomVector r);
  real get_long(){return p1.CalcDistance(p2);};
  void put_pver(GeomVector r){pver=r;};
  real get_a1(){return a1;};
  real get_b1(){return b1;};
  real get_c1(){return c1;};
  real get_a(){return a;};
  real get_b(){return b;};
  real get_c(){return c;};
  bool GetExist(){return Exist;};
  GeomVector GetVecXsi(real xsi);
  GeomPlane OperatingElementaryReflector(GeomVector u);
  GeomPlane OperatingRotator(real rot);
  void Shift(real dx, real dy);
};

#endif
