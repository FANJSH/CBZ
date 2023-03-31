#ifndef GENERALOPTION
#define GENERALOPTION

#include <time.h>
#include <string>
#include <math.h>
#include <vector>
#include "Numeric.h"

using namespace std;

class GeneralOption{
  bool forward;
  real omegao;
  vector<real> omegai;
  int grp;
  real epsf,epsk,epss;
  bool print;
  int outitermax;
  int itcmfd,itmin_cmfd;
 public:
  GeneralOption();
  void PutGrp(int i);
  real GetOmegao(){return omegao;};
  real GetOmegai(int i){return omegai[i];};
  void PutOmegai(int i,real j){omegai[i]=j;};
  void PutAdjointCal(){forward=false;};
  void PutEpsf(real i) {epsf=i;};
  void PutEpsk(real i) {epsk=i;};
  void PutEpss(real i) {epss=i;};
  void PutOutitermax(int i){outitermax=i;};
  void PutItmin_cmfd(int i){itmin_cmfd=i;};
  void PutItcmfd(int i){itcmfd=i;};
  real GetEpsf() {return epsf;};
  real GetEpsk() {return epsk;};
  real GetEpss() {return epss;};
  int GetOutitermax(){return outitermax;};
  int GetItcmfd(){return itcmfd;};
  int GetItmin_cmfd(){return itmin_cmfd;};
  bool Forward(){return forward;};
  bool Print(){return print;};
  bool Converged(real errf,real errk,real errs);
};

#endif
