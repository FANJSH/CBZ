#ifndef PJISYSTEM
#define PJISYSTEM

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
#include "SensitivityData.h"
#include "GData.h"

using namespace std;

class PJISystem:public GeneralSystem{
 protected:
  TrajectorySet *cpij;
  PJISlabPij slabpij;
  vector<GroupData2D> pij;
  vector<GroupData2D> pijr; // perpendicular
  vector<GroupData2D> pijz; // parallel
  bool slab;
  bool write_proc;
  bool put_pijr;
  vector<bool> pij_transform;
  vector<int> medium_id_par_region;

  string LeakTreat;
  enum xstype SigmaCol;
  real b2, b2r, b2z;

  void CalFluxPij(int g,real b2);
  void CalFluxPijNew(int g,real b2,int iter);
  void CalFluxPijModifiedSource(int g);
  void CalFluxTibere(int g);
  void SetInitialFlux();
  real BucklingSearchPij();
  real BucklingSearchMSS();
  real BucklingSearchTibere();
  void PutPij(enum xstype ss,real b2=0.);
  void PutPijr(enum xstype ss,real b2=0.);
 public:
  PJISystem(int grpi,int MAX_MED=20);
  void PutRegMed(vector<int> inp);
  void PutRegMed(int *inp);
  void PutRelationRegionMedium(int *inp){PutRegMed(inp);};
  void PutLeakageTreatment(string inp);
  void PutSigmaCol(enum xstype inp);
  void PutTrajectorySet(TrajectorySet *cinp);
  void PutCartMeshInfo(CartMeshInfo cminp);
  void PutApproximatedCurrent(real b2);
  void PutFluxAsCurrent();
  real CalIgenPij(bool chi_mat=false);
  real CalIgenTibere();
  real CalFluxGeneral(int g,real cin,int iter);
  void PutBuckling(real bi){b2=bi;};
  real BucklingSearch(string opt="PseudoAbsorption",enum xstype sigt=sigtr);
  real CalHomoSigt(int g,enum xstype ss=sigt);
  GroupData2D &GetPij(int grp){return pij[grp];};
  GroupData2D &GetPijr(int grp){return pijr[grp];};
  void PutPij(){PutPij(SigmaCol,b2);};
  void PutPijSphere();
  void PutPijr(){PutPijr(SigmaCol,b2);};
  void PijTransformTrue();
  void PijTransformFalse();
  Medium Homogenize(int hreg, int*nreg,bool benoist_d=true, bool micxs=false); // nreg:mesh-id
  Medium HomogenizeAll(bool benoist_d=true, bool micxs=false);
  void TransformForAdjointFlux();
  // Perturbation calculation
  GData CalReactivityGData(PJISystem *sec,real kunp,real kp,bool pr=true,bool ipcal=true);
  real CalReactivity(PJISystem *sec,real kunp,real kp,bool pr=true,bool ipcal=true);  
  real CalReactivity(PJISystem *sec,real kunp,real kp,vector<real> &rho,bool pr,bool ipcal=true);
  real CalReactivity(PJISystem *sec,real kunp,real kp,bool pr, int meshid,bool ipcal=true);
  real CalPerturbScatteringTerm(PJISystem *sec,int g,bool *flag);
  void CalSensitivity(PJISystem *sec,real kunp,real kp,int nucnum,int *nucid);
  SensitivityData CalSensitivityNew(PJISystem *sec,real keff,int nucnum,int *nucid,bool ipcal=true);
  SensitivityData CalSensitivityRRR(PJISystem &lat,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh,int nucnum,int *nucid,real keff);
  GroupData1D GetIntegratedFlux(int medid,int mom=0);
  real GetVolumePerMedium(int medid);
  void CopyPij(PJISystem &sys2);

  void WriteProcOff(){write_proc=false;};

  //void AddFlux(PJISystem &sec);
  //void NegFlux(PJISystem &sec);
  //void NegFlux(PJISystem &sec,real fact);
  void ClearFlux();
};

#endif
