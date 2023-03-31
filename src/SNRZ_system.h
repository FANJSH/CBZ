#ifndef SNRZSYSTEM
#define SNRZSYSTEM

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

#include "SNRZ_quadrature.h"
#include "PLOSE_system.h"
#include "PLOS_system.h"
#include "SensitivityData.h"

using namespace std;

class SNRZSystem:public GeneralSystem{

  int nquad;
  int BC[6];
  PLOSESystem psys;
  vector<int> quadid;
  vector<SNRZQuadrature*> quad;
  vector<bool> wrtflx;
  vector< vector<GroupData1D> > aflux;
  bool transport; // transport approximation
  bool etransport;
  bool dsa; 
  bool step_d;
  bool array_set;
  // for CMFD
  vector< vector< vector< vector< vector<real> > > > > CurFF;
  vector<real> curxp;
  vector<real> curxn;
  vector<real> curyp;
  vector<real> curyn;
  // for total leakage
  GroupData1D leakage;
  int itmax_inner;
 public:
  SNRZSystem(){};
  SNRZSystem(int Ndiminp,int grp,int MAX_MED=20){Init(Ndiminp,grp,MAX_MED);};
  void Init(int Ndiminp,int grp,int MAX_MED);
  ~SNRZSystem();

  void NoTransportApprox(){transport=false;};
  void OnETransportApprox(){etransport=true;};
  void PutCartMeshInfo(CartMeshInfo inp,string geom="Cylinder");

  void StepDifferencing(){step_d=true;};

  void SetQuadratureNum(int num);
  void SetQuadratureID(int *gbnd, int *id);
  void SetQuadrature(SNRZQuadrature *qinp, int id=0);
  void SetArray();
  SNRZQuadrature* GetQuadrature(int g){return quad[quadid[g]];};

  real CalFlux(int g,int oiter,real epsif=1e-5);
  real CalFluxAllAbsorption(int g,int oiter,real epsif=1e-5);
  real CalFluxFirstCollision(int g,int oiter,real epsif=1e-5);
  void PutSigmaForInnerIteration(int g,real *sigtvol,real *sigsself,real *sigt);
  void SetInitialFlux();
  void InitialFluxGuessInnerIteration(int g,real *flxmom);
  void PutSourceInnerIteration(real *Src);
  void CalFluxInSource(real *src,real *ssrc,real *flxmom,real *ssc,int qid,int is);
  void CalFluxInSourceSP(real *src,real *ssrc,real *flxmom,real *ssc,int qid,real xi,real mu);
  void RenewFluxMomentInnerIteration(int sn,int qid,real *flxmom,real *anflx);
  void AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself,bool cmfd_on);
  real PutFluxAfterInnerIteration(int g,real *flxmom);
  real CalFluxMomentFromAngularFlux(int g,int mom,int meshid);
  GroupData1D GetIntegratedFlux(int medid,int mom=0);
  GroupData1D GetIntegratedFlux(int x1,int x2,int y1,int y2)
  {return GeneralSystem::GetIntegratedFlux(x1,x2,y1,y2);};

  void PutWriteFlux(int med=-1);
  bool GetWrtflx(int m){return wrtflx[m];};
  void NoPrint(){print=false; psys.NoPrint();};

  GroupData1D &GetAFlux(int m,int g){return aflux[m][g];};
  GroupData1D &GetAFlux(int r,int z,int g){return aflux[meshid[0][z][r]][g];};

  void NoDSAAcceleration(){dsa=false;};

  PLOSESystem &GetPLOSESystem(){return psys;};

  // For power iteration
  real CalFluxGeneral(int ginp,real cin,int iter);

  void DoAcceleration(int it,real errf,real k);

  void AddFlux(SNRZSystem &sec);

  // For CMFD
  void DoCMFDAcceleration(real delk);
  void SetZeroCurFF();
  void CalCoarseCur();

  void AddMedium(Medium inp);
  void AddMedium(string mdir,string ss,int plt=1,bool upscat=false);

  GroupData1D &GetLeakage(){return leakage;};
  void ShowKeffContribution();

  // For perturbation calculation
  real CalReactivity(SNRZSystem *sec,real kunp,real kp,bool pr=true);
  real CalReactivity(SNRZSystem *sec,real kunp,real kp,bool* flag,bool pr=true);
  void CalSensitivity(SNRZSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  void CalPerturbLeakScat(SNRZSystem *sec, bool *flag, int g, real *sct,real *leak);
  real CalPerturbLeakScatSimplified(SNRZSystem *sec, bool *flag, int g);
  real CalPerturbLeakScatSimplifiedNew(SNRZSystem *sec, bool *flag, int g,int gg=-1);
  SensitivityData CalSensitivityNew(SNRZSystem *sec,real keff,int nucnum,int *nucid,bool fission_matrix=false);
  SensitivityData CalSensitivityNewFixedSource(SNRZSystem *sec,real val,int nucnum,int *nucid);

  real CalZCurrent(int r,int z,int g,bool pos=true);
  void PutItmaxInner(int i){itmax_inner=i;};

  void CalFixedSourceUpScat(real epsf=1e-4, int oiter=1, bool print=true);
};

#endif
