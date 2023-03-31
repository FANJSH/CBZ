#ifndef SNRZQUADRATURE
#define SNRZQUADRATURE

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "SphericalHarmonics.h"

using namespace std;

class SNRZQuadrature{
 protected:
  int sn, snnum;
  int pl, plnum;
  vector<real> omega, mu, xi;
  vector<real> mc; 
  vector< vector<real> > val;
  vector<int> xref, yref, xyref;
  vector<int> sizex;
  bool exist;
 public:
  SNRZQuadrature(int pl=0);
  ~SNRZQuadrature();
  // Access function for input
  void PutSnnum(int sninp);
  // Access function for output
  bool Exist(){return exist;};
  real GetMC(int i){return mc[i];};
  int GetSn(){return sn;};
  int GetSnnum(){return snnum;};
  int GetPlnum(){return plnum;};
  real GetOmega(int i){return omega[i];};
  real GetMu(int i){return mu[i];};
  real GetXi(int i){return xi[i];};
  int GetSizex(int i){return sizex[i];};
  int GetXref(int i){return xref[i];};
  int GetYref(int i){return yref[i];};
  int GetXYref(int i){return xyref[i];};

  real GetMoment(int m,int is);
  real GetMoment(int m,real i);
  real GetMoment(int l,int m,real mu,real xi);
  real GetMoment(int l,int m,int is);

  void PutLevelSymmetric(int sn);
  void PutS2();
  void PutS4();
  void PutS6();
  void PutS8();
  void PutS12();
  void PutS16();
  void PutS20();
  void PutRectangularDPnTn(int dpn,int tn);
  void PutBiasedRectangularDPnTn(int dpn1,int dpn2,int tn,real mu_b);
  void PutBiasedRectangularDPnTn(int dpn1,int dpn2,int dpn3,int tn,real mu_b);
  // dpn1 for [-1,mu_b], dpn2 for [mu_b,+1]
  void ArrangeDiscretedPoint(real *p,real *w,int *wm);
  void CalMC();
  void CalValue();
  void WeightNormalize();
  void CalXYref();
  void CheckOrthogonalityTo00(real eps=1e-4);
  void show_self();
};


#endif
