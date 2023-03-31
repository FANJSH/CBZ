#ifndef SNRQUADRATURE
#define SNRQUADRATURE

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "SphericalHarmonics.h"

using namespace std;

class SNRQuadrature{
 protected:
  int sn, pl;
  vector<real> omega, mu, mc;
  vector< vector<real> > val;
  vector<int> xref;
  bool exist;
 public:
  SNRQuadrature(int pl);
  ~SNRQuadrature();
  bool Exist(){return exist;};
  void PutSN(int sninp);
  void CalWeightGaussian();
  void CalValue();
  void WeightNormalize();
  void SymmetricMU();
  void PutGaussian(int s,bool odd_permit=false);
  void PutDoubleGaussian(int s);
  void CalMC();
  real GetMC(int i){return mc[i];};
  int GetSN(){return sn;};
  real GetOM(int i){return omega[i]*mu[i];};
  real GetOmega(int i){return omega[i];};
  real GetMu(int i){return mu[i];};
  real GetMoment(int m,int is);
  real GetMoment(int m,real i);
  void CalXYZref();
  int GetXref(int i){return xref[i];};
  void show_self();
  void PutMu(int i,real m){mu[i]=m;};
};


#endif
