#ifndef MECSYSTEM
#define MECSYSTEM

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Medium.h"
#include "CartMeshInfo.h"
#include "Numeric.h"
#include "GeneralMesh.h"
#include "GeneralSystem.h"
#include "PJI_trajectoryset.h"
#include "PJI_slabpij.h"
#include "FunctionTable.h"
#include "Quadrature.h"
#include "PLOS_system.h"

using namespace std;

class MECSystem:public GeneralSystem{
 protected:
  exptab etab;
  funcmoc fmoc;
  TrajectorySet *tset;
  int max_segments;
  Quadrature quad;

  vector< vector< vector<real> > > influx; // Angular flux from out of region
  int polar_angle_division;
  vector<real> polar_angle_sin_inv;
  vector<real> polar_angle_weight;

  int azimuthal_angle_division;
  vector<real> azimuthal_angle_weight;

  bool transport; // transport approximation
  bool set_array;

  bool Vacuum; // 
  bool degree360;
  vector<int> medid_fmesh;
  vector< vector<GroupData1D> > aflux;
  bool wrtflx;
  // For coarse-mesh rebalancing
  bool cmr;
  int cmesh,cmeshx,cmeshy;
  vector<real> cmxl,cmyl;
  vector<int> cmesh_bound;
  vector<int> cmeshid;
  vector<int> fmesh_par_cmesh; // Number of fine meshes in each coarse mesh
  vector<int> regid_45;
  vector<int> regid_90;
  vector<int> regid_180;
  // For CMFD
  vector< vector< vector< vector< vector<real> > > > > CurFF;

  bool already_run;
 public:
  MECSystem(int grpi,int MAX_MED=20);
  //void PutTrajectorySet(TrajectorySet *tsetinp,bool deg360=false,bool deg720=false);
  void PutTrajectorySet(TrajectorySet *tsetinp,bool deg720=false);
  void PutRelationRegionMedium(int *inp);
  void PutRelationRegionMedium(vector<int> &inp);
  void PutRegMed(int *inp){PutRelationRegionMedium(inp);};
  void PutRegMed(vector<int> &inp){PutRelationRegionMedium(inp);};
  void SetInitialFlux();
  void SetArray();
  void PutWriteFlux(){wrtflx=true;};
  bool WriteFlux(){return wrtflx;};
  GroupData1D &GetAFlux(int m,int g){return aflux[m][g];};
  Quadrature &GetQuad(){return quad;};

  void NoTransportApprox(){transport=false;};
  void NoTransportApproximation(){transport=false;};

  real CalFluxGeneral(int ng,real cin,int iter);
  real CalFluxGeneralDetail(int ng,real cin,int iter);
  real CalFluxAngle360Old(int ng,real cin,int iter,bool angular_dependent_source=false);
  real CalFluxAngle360(int ng,real cin,int iter);
  real CalFluxDegree360(int ng,real cin,int iter,bool angular_dependent_source=false){return CalFluxAngle360(ng,cin,iter);};
  real CalFluxAngle360AngularDependentSource(int ng,real cin,int iter,bool benoist_consistent=false);
  real CalFluxDegree720(int ng,real cin,int iter);
  void PutSigmaForInnerIteration(int g,real *sigt,real *sigsself);
  void InitialFluxGuessInnerIteration(int g,real *flx);
  void PutSourceInnerIteration(real *src);

  void PutPolarAngleDivision(int i);
  int GetPolarAngleDivision(){return polar_angle_division;};
  void PutAzimuthalAngleDivision(int i);

  void PutVacuum(){Vacuum=true;};
  // Perturbation calculation
  real CalReactivity(MECSystem *sec,real kunp,real kp,bool pr=true,bool ipcal=true);
  real CalReactivity(MECSystem *sec,real kunp,real kp,vector<real> &rho,bool pr,bool ipcal=true);
  real CalPerturbScatteringTerm(MECSystem *sec,int g,bool *flag);
  real CalPerturbP1LeakageTerm(MECSystem *sec,int g,bool *flag);
  SensitivityData CalSensitivityNew(MECSystem *sec,real keff,int nucnum,int *nucid,bool ipcal=true);
  SensitivityData CalSensitivityRRR(MECSystem &lat, int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh,int nucnum,int *nucid,real keff);
  void CalGPT_MEC(real keff,real eps_f,int itermax);
  // For coarse-mesh rebalancing
  void PutCMesh(int i,int j);
  void PutCMeshBound(int *i);
  void PutCMeshLength(real *xin,real *yin);
  void NoCMRAcceleration(){cmr=false;};
  void PutCellSymmetricData(int i,int *i1,int *i2,int *i3);
  void PutCellSymmetricDataForAllOctant(int sqr=1);
  void PutCellSymmetricDataForAllQuarter();
  void PutCellSymmetricDataForRectangular(int m);
  void PutCellSymmetricDataNoSymmetry();
  int GetCMeshBound(int cid){return cmesh_bound[cid];};
  // For CMFD for outer iteration
  void DoAcceleration(int it,real errf,real k);
  void DoCMFDAcceleration(real delk);
  void SetZeroCurFF();
  GroupData1D GetIntegratedFlux(int medid,int mom=0);
  GroupData1D GetIntegratedReactionRate(enum xstype react, int medid);
  real GetVolumeParMedium(int medid);
};


#endif
