#ifndef FUNCTIONTABLE
#define FUNCTIONTABLE

#include "Numeric.h"

// *Caution*
// fk(4) is already multiplied by 1.5
// for anisotropic pij calculation

class fktab{
 protected:
  real fk_a[5500];
  real fk_b[5500];
  real fk_c[5500];
  real errx;
  real aa[24];
  real bb[30];
  real cc[6];
  real px[6];
 public:
  fktab(){};
  void init(real e);
  real fkin(int n,real x);
  real xkini(real x,real y,int n);
  real simpsk(real err,real y,int n);
  real fkin0(real y,int n,real err);
};

class exptab
{
  float p0[1001];
  float p1[1000];
 public:
  exptab();
  real e(real x);
};

class funcmoc
{
  int num, num2;
  real wid; 
  vector<real> p0,p1;
 public:
  funcmoc();
  real get(real x);
};

#endif
