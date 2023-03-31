#ifndef TRAJECTORY
#define TRAJECTORY

#include <fstream>
#include <vector>
#include <iostream>
#include <cstdio>
#include"FunctionTable.h"
#include"GeomVector.h"
#include"Numeric.h"

using namespace std;

class Trajectory{
 protected:
  int numreg;
  vector<int>  reg;
  vector<real> dist;
  vector<real> opt;
  vector<real> xsinv;
  real Weight; // Width*(PI/angle_div);
  real TotalOpt;
 public:
  Trajectory(){numreg=-1;};
  Trajectory(int i){Initialize(i);};
  void Initialize(int i);
  void PutData(int *ireg, real *idist);
  void PutData(vector<int> ireg, vector<real> idist);
  void PutDist(int i,real in){dist[i]=in;};
  void PutXS(real *xs);
  void CalVol(int regnum,real *vol){GetRegionwiseLength(regnum,vol);};
  void GetRegionwiseLength(int regnum,real *length);
  // Pij calculation for single lattice
  void CalPij(int matnum,fktab &fin,real *retpij,bool aniso=false);
  void CalPijSphere(int matnum,real *retpij);
  void CalPijSphereWithExpTable(int matnum,exptab &etab,real *retpij);
  void CalPijAniso(int matnum,fktab &fin,real *retpij);
  // Pij calculation for periodic lattice
  void CalPij
    (int matnum,Trajectory &sec,real Optb,real *pij,fktab &fin,bool aniso,bool oppo1,bool oppo2);
  void CalPijIso
    (int matnum,Trajectory &sec,real Optb,real *pij,fktab &fin,bool oppo1,bool oppo2);
  void CalPijAniso
    (int matnum,Trajectory &sec,real Optb,real *pij,fktab &fin,bool oppo1,bool oppo2);
     // oppo=1 - Opposite direction
  // +++ for CBG/MEC
  real SolveCFormTransport(exptab &e,real fin,real *src,real *avef);
  real SolveCFormTransport(exptab &e,real fin,real *src,real *avef,real sin_inv);
  void SolveCFormTransport(exptab &e,real fin,real *src,real *avef,real *outf,real sin_inv);
  void SolveCFormTransportOpposite(exptab &e,real fin,real *src,real *avef,real *outf,real sin_inv);
  void SolveCFormTransport(funcmoc &fmoc,real fin,real *src,real *avef,real *outf,real sin_inv);
  void SolveCFormTransport(funcmoc &fmoc,real fin,real *src,real *avef,real &outf,real sin_inv);
  void SolveCFormTransport(funcmoc &fmoc,real fin,vector<real> &src,vector<real> &avef,real &outf,real sin_inv);
  void SolveCFormTransport(funcmoc &fmoc,real fin,real *src,real *avef,real *outf,real *ep,real sin_inv);
  void SolveCFormTransportOpposite(funcmoc &fmoc,real fin,real *src,real *avef,real &outf,real sin_inv);
  void SolveCFormTransportOpposite(funcmoc &fmoc,real fin,real *src,real *avef,real *outf,real sin_inv);
  void SolveCFormTransportOpposite(funcmoc &fmoc,real fin,real *src,real *avef,real *outf,real *ep,real sin_inv);
  void SolveCFormTransportDD(real flxin,real *src,real *aveflx,real *outflx,real sin_inv);
  void SolveCFormTransportDoubleDirection(funcmoc &fmoc, real fin1, real fin2, real *src, real *avef1, real *avef2, real &outf1, real &outf2, real sin_inv);
  void SolveCFormTransportDoubleDirection(funcmoc &fmoc, real fin1, real fin2, real *src, real *avef1, real *avef2, real *outf1, real *outf2, real sin_inv);
  void SolveCFormTransportDoubleDirection(funcmoc &fmoc, real fin1, real fin2, real *src, real *avef1, real *avef2, vector<real> &outf1, vector<real> &outf2, real sin_inv);
  void show_self();
  void AdjustTrajectoryLength(Trajectory &sec,Trajectory &thi);
  void LengthFactorize(real factor);
  // Access function
  int GetNumreg(){return numreg;};
  int GetReg(int i){return reg[i];};
  int GetLastReg(){return reg[numreg-1];};
  real GetDist(int i){return dist[i];};
  real GetOpt(int i){return opt[i];};
  real GetWeight(){return Weight;};
  real GetTotalOpt(){return TotalOpt;};
  bool GetExist(){if(numreg>0)return true; return false;};
  void PutWeight(real i){Weight=i;};
};

#endif
