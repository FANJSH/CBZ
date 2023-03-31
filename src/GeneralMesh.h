#ifndef GENERALMESH
#define GENERALMESH

#include <string>
#include <math.h>
#include <vector>
#include "Medium.h"

using namespace std;

class GeneralMesh{
 protected:
  int dim, grp, pl, maxpl;
  real Len[3], Sur[6], Volume;
  vector<GroupData1D> Flux; // Point flux angular-moment
  vector<GroupData1D> SSrc; // Scattering source (volume-integrated) 
  vector<GroupData1D> UpSSrc; // Up scattering source 
  vector<real> srcin;      // External source for inner iteration (volume-integrated)
  real FSrc;  // Fission source (angle- and volume-integrated)
  Medium *Med;
 public:
  //Constructor & Destructor
  GeneralMesh();
  GeneralMesh(int grp){PutGrp(grp);};
  ~GeneralMesh();

  void PutDim(int i);
  void PutDim(int n, real *l);
  void PutDimSphere(real l,real ll,real lr);
  void PutDimCylinder(real l,real ll,real lr,real z);

  void SetVectorSizeForFlux();
  void PutGrp(int g);
  void PutPL(int pl=0);

  real GetDif(int g,int flag);

  void SetZeroScatSrc(int g);
  void SetZeroScalarFlux(int g);
  void SetZeroUpScatSrc(int g);
  void AddDownScatSrc(int DepGrp,int pls=0,bool selfscat=false);
  void AddDownScatSrcReactionWise(int DepGrp,enum xstype xst,int pls=0,bool selfscat=false);
  void AddUpScatSrc(int DepGrp,int pls=0);
  void AddDownScatSrcAdjoint(int DepGrp,int pls=0);
  void AddDownScatSrcAdjointArv(int DepGrp,int ArvMinGrp,int pls=0);
  // The above method is to reduce calculation time for adjoint-thermal iteration
  void AddUpScatSrcAdjoint(int DepGrp,int pls=0);
  real GetScatSrc(int g,int ind=0);
  real GetUpScatSrc(int g,int ind=0);
  void PutScatSrc(int g,real inp,int ind=0){SSrc[ind].put_data(g,inp);};
  void AddScatSrc(int g,real inp,int ind=0){SSrc[ind].add_data(g,inp);};

  void PutFissionSrc(real i){FSrc=i;};
  real GetFissionSrc(){return FSrc;};
  void CalFissionSrc();
  void CalFissionSrcAdjoint();
  real GetFissionSrcChi(int g);
  real GetFissionSrcKai(int g){return GetFissionSrcChi(g);};
  real GetFissionSrcSigf(int g);
  GroupData1D GetFissionSrcSigf();

  real GetReactionRate(enum xstype ss);
  real GetReactionRate(int matnum,enum xstype ss);
  real GetVolumeFlux();

  void PutMedium(Medium *inp){Med=inp;};
  Medium *GetMed(){return Med;};
  void PutVolume(real inp){Volume=inp;}
  real GetVolume(){return Volume;};
  void PutLen(int i,real j){Len[i]=j;};
  real GetLen(int i){return Len[i];};
  real GetSurR(int i){return Sur[i*2+1];};
  real GetSurL(int i){return Sur[i*2];};
  int GetPL(){return pl;};
  int GetMaxPL(){return maxpl;};
  GroupData1D &GetFlux(int i=0){return Flux[i];};
  real GetEnergyIntegratedFlux(real etop,real elow,int mom=0);
  void PutSrcin(real inp,int ind=0){srcin[ind]=inp;};
  void AddSrcin(real inp,int ind=0){srcin[ind]+=inp;};
  real GetSrcin(int ind=0){return srcin[ind];};

  real GetFluxData(int pl, int group){return Flux[pl].get_dat(group);}
  void PutFluxData(int pl, int group, real val){Flux[pl].put_data(group, val);}
  void AddFluxData(int pl, int group, real val){Flux[pl].add_data(group, val);}
  real GetMediumData1D(int group, enum xstype rname){return Med->GetData1D(group, rname);}

  void InitializeUpScatSrc(int wgrp);

  void SrcDataClear();
};

#endif
