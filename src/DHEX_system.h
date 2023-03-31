#ifndef DHEXSYSTEM
#define DHEXSYSTEM

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

class DHEXSystem:public GeneralSystem{
  int BC[6];
  int difc[3];
  vector< vector<real> > coefs;
  vector< vector<real> > coefx;
  vector< vector<real> > coefyl;
  vector< vector<real> > coefyr;
  vector< vector<real> > coefz;
  // hexagonal-z
  real pitch,half_pitch,inv_pitch,hex_area,len_plane;
  vector< vector<bool> > ledge_yl;
  vector< vector<bool> > ledge_yr;
  vector< vector<bool> > redge_yl;
  vector< vector<bool> > redge_yr;
  // For modified coarse mesh correction
  vector< vector<real> > mcor;
  bool cal_mcor;
  // For triangular mesh correction
  //
  // (index for triangular mesh)
  //            [y=0] 
  //    5 0
  //   4   1
  //    3 2
  //            [y=Y]

  bool cal_tri;
  bool flx_store;
  vector< vector<real> > tricor1; // leakage from outer-bound
  vector< vector<real> > tricor2; // leakage between adjacent meshes in the same assembly
  vector< vector<real> > outer_fsrc_cal;
  vector< vector<real> > outer_fsrc; // [m,6]
  vector< vector< vector<real> > > outer_ssrc;
  vector< vector< vector<real> > > flx_triang; // [g,m,6]
  real keff_tmp;
 public:
  DHEXSystem(){};
  DHEXSystem(int Ndiminp,int grp,int MAX_MED=20){Init(Ndiminp,grp,MAX_MED);};
  ~DHEXSystem(){};

  void Init(int Ndiminp,int grp,int MAX_MED);
  void PutCartMeshInfo(CartMeshInfo inp, real pitchinp);
  void CalCoef(bool mcor=false);
  void CopyCoef(DHEXSystem &sec);
  real CalEdgeCurVacuum(int m,int dir,int g,int flag,bool mcor=false);
  void EdgeCalculationForHex();
  void SetInitialFlux();
  void SetDifc(string d1,string d2,string d3);
  int GetDifc(int i){return difc[i];};
  real CalCurrentX(int x,int y,int z,int g); 
  real CalCurrentZ(int x,int y,int z,int g); 
  real CalCurrentYL(int x,int y,int z,int g); 
  real CalCurrentYR(int x,int y,int z,int g); 
  real CalCurrentX_Tri(int x,int y,int z,int g); 
  real CalCurrentZ_Tri(int x,int y,int z,int g); 
  real CalCurrentYL_Tri(int x,int y,int z,int g); 
  real CalCurrentYR_Tri(int x,int y,int z,int g); 
  // For power iteration
  real CalFluxGeneral(int g,real cin,int it);
  real CalFlux(int grp,int iter,real epsif=1e-7);
  real CalFluxOmega(int grp,real epsif);
  real CalFluxOmegaTriangle(int grp,real epsif);
  bool SweepInnerIteration(int itmax,real epsif,int g,real omega,real *Fl,real *Src);
  bool SweepInnerIterationTriangle(int itmax,real epsif,int g,real omega,real *Fl,real *Src);
  real RenewFluxInnerIteration(int g,real *Fl);
  void GetFluxAndSrcInnerIteration(int g,real *Fl,real *Src);
  // For perturbation calc.
  real CalPerturbDenominatorTriangular(DHEXSystem *sec);
  real CalPerturbLeakageTermInsideAssembly(DHEXSystem *sec,bool *flag,int g);
  void CalPerturbLeakageTermInsideAssembly(DHEXSystem *sec,vector<real> &sol,int g);
  real CalPerturbLeakageTerm(DHEXSystem *sec,bool *flag,int g,int dir);
  void CalPerturbLeakageTerm(DHEXSystem *sec,vector<real> &sol,int g,int dir);
  real CalReactivity(DHEXSystem *sec,real kunp,real kp,bool pr=true);
  real CalReactivity(DHEXSystem *sec,real kunp,real kp,bool* flag,bool pr=true);
  real CalReactivityTriangular(DHEXSystem *sec,real kunp,real kp,bool* flag,bool pr=true);
  void CalSensitivity(DHEXSystem *sec,real kunp,real kp,int nucnum,int *nucid);
  SensitivityData CalSensitivityNew(DHEXSystem *sec,real keff,int nucnum,int *nucid);
  // For CMFD acceleration
  void DoAcceleration(int iter,real errs,real fiss);
  void DoCMFDAcceleration(real delk);
  real CalIgenCoarse(vector< vector< vector<real> > > &moddif_x,
                     vector< vector< vector<real> > > &moddif_yl,
                     vector< vector< vector<real> > > &moddif_yr,
                     vector< vector< vector<real> > > &moddif_z);
  real CalFluxModifiedLeakage(int g,real epsif,
                 vector< vector< vector<real> > > &moddif_x,
                 vector< vector< vector<real> > > &moddif_yl,
                 vector< vector< vector<real> > > &moddif_yr,
  	         vector< vector< vector<real> > > &moddif_z);
  void CalCoarseCur(int g,vector< vector< vector<real> > > &curx,
		    vector< vector< vector<real> > > &curyl,
		    vector< vector< vector<real> > > &curyr,
		    vector< vector< vector<real> > > &curz, bool cumulative);
  void PreIgenCoarse(vector< vector< vector<real> > > &fcur_x,
                     vector< vector< vector<real> > > &fcur_yl,
                     vector< vector< vector<real> > > &fcur_yr,
                     vector< vector< vector<real> > > &fcur_z,
                     vector< vector< vector<real> > > &moddif_x,
                     vector< vector< vector<real> > > &moddif_yl,
                     vector< vector< vector<real> > > &moddif_yr,
                     vector< vector< vector<real> > > &moddif_z,
                     vector< vector<real> > &flxf);
  void CalCorDif(vector< vector< vector<real> > > &ccur_x,
                 vector< vector< vector<real> > > &ccur_yl,
                 vector< vector< vector<real> > > &ccur_yr,
                 vector< vector< vector<real> > > &ccur_z,
                 vector< vector< vector<real> > > &fcur_x,
                 vector< vector< vector<real> > > &fcur_yl,
                 vector< vector< vector<real> > > &fcur_yr,
                 vector< vector< vector<real> > > &fcur_z);
  void ModCoef(vector< vector< vector<real> > > &moddif_x,
               vector< vector< vector<real> > > &moddif_yl,
               vector< vector< vector<real> > > &moddif_yr,
               vector< vector< vector<real> > > &moddif_z);
  void PrintFluxDistribution(int z=0);
  void MemoryReductionForHighMomentCal();

  void ModifiedCoarseMeshCorrectionCal(real k);
  real GetMCMCor(int m,int g){return mcor[m][g];};
  real CalIgenMCMCor(int pl=0);

  real GetPitch(){return pitch;};

  void PutTriangularMeshCorrection(bool flx_store_inp=true);
  real GetFluxTriangularMesh(int g,int m,int pos){return flx_triang[g][m][pos];};
  real GetFissionSourceTriangularMesh(int m,int pos){return outer_fsrc[m][pos];};
  bool IsTriangularCorrection(){return cal_tri;};
  //
  real GetCoefs(int i,int j){return coefs[i][j];};
  real GetCoefx(int i,int j){return coefx[i][j];};
  real GetCoefyl(int i,int j){return coefyl[i][j];};
  real GetCoefyr(int i,int j){return coefyr[i][j];};
  real GetCoefz(int i,int j){return coefz[i][j];};

  void Get24PointFlux(int g,int m,real *res);
  void Get24PointFlux(int g,int x,int y,int z,real *res);
  /*
  void DoAcceleration(int it,real errf,real k);
  real CalEdgeCurZeroFlux(int m,int dir,int g,int flag);
  */
};

#endif
