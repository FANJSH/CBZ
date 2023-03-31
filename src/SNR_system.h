#ifndef SNRSYSTEM
#define SNRSYSTEM

#include <time.h>
#include <string.h>
#include <math.h>
#include "GroupData.h"
#include "Medium.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"
#include "GeneralMesh.h"
#include "GeneralSystem.h"

#include "SNR_quadrature.h"
#include "PLOSE_system.h"
#include "SensitivityData.h"
#include<iostream>
#include <vector>

using namespace std;

class SNRSystem:public GeneralSystem{

  int nquad;
  int BC[6];
  PLOSESystem psys;
  vector<int> quadid;
  vector<SNRQuadrature*> quad;
  vector<bool> wrtflx;
  vector< vector<GroupData1D> > aflux;
  bool transport; // transport approximation
  bool etransport;
  bool sphere;
  bool dsa;
  // for total leakage
  GroupData1D leakage;
  // for leakage sensitivity calculation
  bool leak_sens;
  int leak_mesh,top_grp,btm_grp;
  bool step_d;
  bool surface_aflx_store;
  // +++ for white boundary condition
  vector<real> inflx_white;

 public:
  SNRSystem(int Ndiminp,int grp,int MAX_MED=20);
  ~SNRSystem();

  void NoTransportApprox(){transport=false;};
  void OnETransportApprox(){etransport=true;};
  void PutCartMeshInfo(CartMeshInfo inp,string geom="Sphere");

  void StepDifferencing(){step_d=true;};
  void SurfaceFluxStore(){surface_aflx_store=true;};

  void SetQuadratureNum(int num);
  void SetQuadratureID(int *gbnd, int *id);
  void SetQuadrature(SNRQuadrature *qinp, int id=0);
  void SetArray();
  SNRQuadrature* GetQuadrature(int g){return quad[quadid[g]];};

  real CalFlux(int g,real epsif=1e-5);
  real CalFluxSlab(int g,real epsif=1e-5,bool angular_dependent_src=false);
  //real CalFluxSlabDF(int g,real epsif,vector<real> inflx_df,vector< vector<real> > outflx_df_l,vector< vector<real> > outflx_df_r);
  real CalFluxSlabDF(int g,real epsif,vector<real> inflx_df,vector<real> outflx_df_l,vector<real> outflx_df_r);
  real CalFluxSlabDF(int g,real epsif,vector<real> outflx_df_l,vector<real> outflx_df_r);
  void CalFluxInSource(real *ssrc,real *flxmom,real *ssc,int qid,int sn);
  void CalFluxInSourceSP(real *ssrc,real *flxmom,real *ssc,int qid,int sn);
  //void CalFluxInSourceSP(real *ssrc,int g);
  void AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself);

  void PutWriteFlux(int med=-1);
  bool GetWrtflx(int m){return wrtflx[m];};

  GroupData1D &GetAFlux(int m,int g){return aflux[m][g];};

  void NoDSAAcceleration(){dsa=false;};

  void ClearFlux();
  void AddFlux(SNRSystem &sec);
  void NegFlux(SNRSystem &sec);
  void NegFlux(SNRSystem &sec,real fact);

  void FluxAngleReverse();
  // For power iteration
  real CalFluxGeneral(int ginp,real cin,int iter);
  void SetInitialFlux();
  void DoAcceleration(int it,real errf,real k){};

  // For perturbation calculation
  real CalReactivity(SNRSystem *sec,real kunp,real kp,vector<real> &rho,bool pr=true);
  real CalReactivity(SNRSystem *sec,real kunp,real kp,bool pr=true);
  real CalReactivity(SNRSystem *sec,real kunp,real kp,int g);
  void CalSensitivity(SNRSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  void CalSensitivityScatteringMatrix(SNRSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  void CalSensitivityFixedSource(SNRSystem *sec,int nucnum,int *nucid);
  SensitivityData CalSensitivityNewFixedSource(SNRSystem *sec,real val,int nucnum,int *nucid);
  void CalPerturbLeakScat(SNRSystem *sec,bool *flag,int g,real *scat,real *leak);
  real CalPerturbLeakScatSimplifiedNew(SNRSystem *sec, bool *flag, int g,int gg=-1);
  SensitivityData CalSensitivityNew(SNRSystem *sec,real keff,int nucnum,int *nucid,bool fiss_matrix=false,bool ipcal=true);

  SensitivityData CalSensitivityRRR(SNRSystem &lat,int nume_id,enum xstype nume_xs,int denom_id,enum xstype denom_xs,int mesh_pos,int nucnum,int *nucid,real keff);
  /*
  void CalGPT_SNR(real keff,int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh,real rr,int itermax=20);
  void CalGPT_SNR(real keff,real eps_f=1e-4,int itermax=20);
  */
  // *** Note in 2018/3/6
  //   Implementation of GPT capability into SNR was failed 
  //   because angular flux and high angular moment of generalized adjoint flux
  //   should be stored.  In MEC or PJI, P0 transport treatment is usually adopted,
  //   so this issue is NOT problematic.  But small critical assembly calculations
  //   with SNR requires high-order angular treatment, so this issue becomes critical.

  void AddMedium(Medium inp);
  GroupData1D &GetLeakage(){return leakage;};
  //
  // For sensitivity
  void PutLeakSens(){leak_sens=true;};
  void PutLeakMesh(int i){leak_mesh=i;};
  void PutTopGrp(int i){top_grp=i;};
  void PutBtmGrp(int i){btm_grp=i;};
};

#endif
