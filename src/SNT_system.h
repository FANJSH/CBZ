#ifndef SYSTEM
#define SYSTEM

#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>

#include "GroupData.h"
#include "Medium.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"
#include "GeneralMesh.h"
#include "GeneralSystem.h"
#include "Quadrature.h"
#include "PLOSE_system.h"
#include "PLOS_system.h"
#include "FunctionTable.h"
#include "SensitivityData.h"

using namespace std;

class SNTSystem:public GeneralSystem{

  int nquad;
  int BC[6];
  PLOSESystem psys;
  vector<int> quadid;
  vector<Quadrature*> quad;
  vector< vector<GroupData1D> > aflux;
  // Storing angular flux
  vector<bool> wrtflx;
  int wrtflx_endgrp;
  bool transport; // transport approximation
  bool dsa; // DSA acceleration
  //exptab etab;
  // for CMFD
  vector< vector< vector< vector< vector<real> > > > > CurFF;
  vector<real> curxp;
  vector<real> curxn;
  vector<real> curyp;
  vector<real> curyn;
  vector<real> curzp;
  vector<real> curzn;
  int itmax_inner;
 public:
  SNTSystem(){};
  SNTSystem(int Ndiminp,int grp,int MAX_MED=20){Init(Ndiminp,grp,MAX_MED);};
  ~SNTSystem();
  void Init(int Ndiminp,int grp,int MAX_MED);

  void NoTransportApprox(){transport=false;};
  void PutCartMeshInfo(CartMeshInfo inp,string geom="Cartesian");

  void SetQuadratureNum(int num);
  void SetQuadratureID(int *gbnd, int *id);
  void SetQuadrature(Quadrature *qinp, int id=0);
  void SetArray();
  Quadrature* GetQuadrature(int g){return quad[quadid[g]];};

  real CalFluxXYZ(int g,int oiter,real epsif=1e-5);
  real CalFluxXYZAllAbsorption(int g,int oiter,real epsif=1e-5);
  real CalFluxXY(int g,int oiter,real epsif=1e-5,bool angular_dependent_source=false);
  void PutSigmaForInnerIteration(int g,real *sigtvol, real *sigsself);
  void SetInitialFlux();
  void InitialFluxGuessInnerIteration(int g,real *flxmom);
  void PutSourceInnerIteration(real *Src);
  void CalSelfScatteringSource(real *Src,real *ssrc,real *flxmom,real *ssc,int qid,int is);
  void RenewFluxMomentInnerIteration(int sn,int qid,real *flxmom,real *anflx);
  void AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself,real *anflx,bool cmfd_on);
  real PutFluxAfterInnerIteration(int g,real *flxmom);
  void WriteAngularFlux(int sn,int g,real *anflx);

  void PutWriteFlux(int med=-2);
  void PutWriteFlux(int med1,int med2);
  bool GetWrtflx(int m){return wrtflx[m];};
  void PutWriteFluxEndGrp(int i);
  void NoPrint(){print=false; psys.NoPrint();};

  GroupData1D &GetAFlux(int m,int g){return aflux[m][g];};

  void NoDSAAcceleration(){dsa=false;};
  void DSAAcceleration(){dsa=true;};

  PLOSESystem &GetPLOSESystem(){return psys;};
  // For power iteration
  real CalFluxGeneral(int g,real cin,int it);
  void DoAcceleration(int it,real errf,real k);

  // For CMFD
  void DoCMFDAcceleration(real delk);
  void SetZeroCurFF();
  void CalCoarseCur();

  void AddMedium(Medium inp);

  // For perturbation calculation
  real CalReactivity(SNTSystem *sec,real kunp,real kp,bool pr=true);
  void CalSensitivity(SNTSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  SensitivityData CalSensitivityNew(SNTSystem *sec,real keff,int nucnum,int *nucid);
  void CalPerturbLeakScat(SNTSystem *sec,bool *flag,int g,real *sct,real *leak);
  
  void MemoryReductionForPerturbation();
  void AddFlux(SNTSystem &sec);
  void NegFlux(SNTSystem &sec);
  void NegFlux(SNTSystem &sec,real fact);
  void PutItmaxInner(int i){itmax_inner=i;};
};

#endif
