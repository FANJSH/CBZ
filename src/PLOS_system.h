#ifndef PLOSSYSTEM
#define PLOSSYSTEM

#include <iostream>
#include <time.h>
#include <string>
#include <vector>
#include <math.h>
#include "GroupData.h"
#include "CartCore.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"
#include "GeneralSystem.h"
#include "GeneralMesh.h"
#include "Medium.h"
#include "SensitivityData.h"

using namespace std;

class PLOSSystem:public GeneralSystem{
  int BC[6];
  int difc[3];
  vector< vector<real> > coefs;
  vector< vector<real> > coefx;
  vector< vector<real> > coefy;
  vector< vector<real> > coefz;
  bool Cylinder,Sphere;
  int size_coef;
  // (for heat conductivity calculation)
  real xr_bc; 
 public:
  PLOSSystem(){};
  PLOSSystem(int Ndiminp,int grp,int MAX_MED=20){Init(Ndiminp,grp,MAX_MED);};
  ~PLOSSystem(){};
  void Init(int Ndiminp,int grp,int MAX_MED);
  void CalCoef();
  void CalCoef1DWithDF();
  void CopyCoef(PLOSSystem &sec);
  void PutCartMeshInfo(CartMeshInfo inp, string geom="Cartesian");

  real CalFlux(int grp,int iter,real epsif=1e-7);
  real CalFluxOmega(int grp,real epsif=1e-7);
  real CalFluxModifiedLeakage(int g,real epsif,
    vector< vector< vector< vector<real> > > >  &d1);
  real RenewFluxInnerIteration(int g,real *Fl);
  bool SweepInnerIteration(int itmax,real epsif,int g,real omega,real *Fl,real *Src);
  void GetFluxAndSrcInnerIteration(int g,real *Fl,real *Src);
  void DetermineSizeCoef();
  void CalInnerCoef(int g, real *tmparray);

  void SetDifc(string d1,string d2,string d3);
  int GetDifc(int i){return difc[i];};
  bool cylinder(){return Cylinder;};
  bool sphere(){return Sphere;};

  // For power iteration
  real CalFluxGeneral(int g,real cin,int it);
  void SetInitialFlux();
  void DoAcceleration(int it,real errf,real k);

  // For CMFD
  void DoCMFDAcceleration(real delk);
  void CalCoarseCur(int g,vector< vector< vector< vector<real> > > > &cur, bool cumulative=true);

  real CalSP3(bool rigourous_bc=false,real alpha=2.);
  real CalSP3(PLOSSystem &sys2, bool rigourous_bc=false,real alpha=2.);
  real CalSP3Adjoint(PLOSSystem &sys2, bool rigourous_bc=false,real alpha=2.);
  real CalSP3Adjoint(bool rigourous_bc=false,real alpha=2.);
  real CalOSP3Adjoint(PLOSSystem &sys2, bool rigourous_bc=false);
  real CalOSP3Adjoint(bool rigourous_bc=false);
  real CalSP3AdjointOld(PLOSSystem &sys2, bool rigourous_bc=false,real alpha=2.);

  void CalCorDif
  (vector< vector< vector< vector< vector<real> > > > > &curc,
   vector< vector< vector< vector< vector<real> > > > > &curf); // Calculate correction for diffusion coefficient
  void ModCoef(vector< vector< vector< vector< vector<real> > > > > &moddif); // Modify coefs from corrected diffision coefficient
  void PreIgenCoarse(vector< vector< vector< vector< vector<real> > > > > &curf,
                     vector< vector<real> > &flxf,
                     vector< vector< vector< vector< vector<real> > > > > &moddif);
  real CalCurrent(int x,int y,int z,int g,int dir); 
  real CalEdgeCurZeroFlux(int m,int dir,int g,int flag);
  real CalEdgeCurVacuum(int m,int dir,int g,int flag);
  real CalIgenCoarse(vector< vector< vector< vector< vector<real> > > > > &moddif);

  // For perturbation calc.
  real CalPerturbLeakageTerm(int dir,PLOSSystem *sec,int g,bool *flag);
  real CalReactivity(PLOSSystem *sec,real kunp,real kp,bool pr=true);
  real CalReactivity(PLOSSystem *sec,real kunp,real kp,bool* flag,bool pr=true);
  real CalReactivitySP3(PLOSSystem *adj_p2, PLOSSystem *fwd_p0, PLOSSystem *fwd_p2, real kunp, real kup,bool pr=true);
  real CalReactivityOSP3(PLOSSystem *adj_p2, PLOSSystem *fwd_p0, PLOSSystem *fwd_p2, real kunp, real kup,bool pr=true);  
  void CalSensitivity(PLOSSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  SensitivityData CalSensitivityNew(PLOSSystem *sec,real keff,int nucnum,int *nucid,bool fiss_matrix=false);
  //void PutPerturbationDenominator(real val);

  void MemoryReductionForPerturbation();
  void MemoryReductionForHighMomentCal();
  void PrintFluxDistributionForCylinder();
  real GetCoefs(int i,int j){return coefs[i][j];};
  real GetCoefx(int i,int j){return coefx[i][j];};
  real GetCoefy(int i,int j){return coefy[i][j];};
  real GetCoefz(int i,int j){return coefz[i][j];};
  void PutCoefs(int i,int j,real k){coefs[i][j]=k;};
  void PrintFluxDistribution(int z=0);

  // (for heat conductivity calculation)
  void PutXrBoundaryCondition(real inp){xr_bc=inp;}; 
  
};

#endif
