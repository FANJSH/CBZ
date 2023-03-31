#ifndef PLOSESYSTEM
#define PLOSESYSTEM

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

using namespace std;

class PLOSESystem:public GeneralSystem{
  int BC[6];
  int difc[3];
  vector< vector<real> > coefs;
  vector< vector<real> > coefx;
  vector< vector<real> > coefy;
  vector< vector<real> > coefz; // o
  bool Cylinder,Sphere; // o
  int size_coef;
  int TotME;
  vector<int> xedgel_e; // o
  vector<int> xedger_e; // o
  vector<int> yedgel_e; // o
  vector<int> yedger_e; // o
  vector< vector< vector< vector<int> > > > edge_e;
 public:
  PLOSESystem(){};
  PLOSESystem(int Ndiminp,int grp,int MAX_MED=20){Init(Ndiminp,grp,MAX_MED);};
  void Init(int Ndiminp,int grp,int MAX_MED);

  void CalCoef();
  void CalCoefSphere();
  void CalCoefCylinder();
  void PutCartMeshInfo(CartMeshInfo inp, string geom="Cartesian");

  real CalFlux(int grp,int it,real epsif=1e-7);  
  real CalFluxOmega(int grp,real epsif=1e-7);
  real RenewFluxInnerIteration(int g,real *FlE);
  real RenewFluxInnerIterationSphere(int g,real *FlE);
  real RenewFluxInnerIterationCylinder(int g,real *FlE);
  void SetEdgeSource(real *SrcE);
  void SetEdgeSourceSphere(real *SrcE);
  void SetEdgeSourceCylinder(real *SrcE);
  void SetInitialFluxInnerIteration(int g,real *FlE);
  void SetInitialFluxInnerIterationSphere(int g,real *FlE);
  void SetInitialFluxInnerIterationCylinder(int g,real *FlE);
  bool SweepInnerIteration(int itmax,real epsif,int g,real* FlE,real* SrcE);
  void DetermineSizeCoef();
  void CalInnerCoef(int g,real *tmparray);

  void SetDifc(string d1, string d2, string d3);
  int GetDifc(int i){return difc[i];};

  bool cylinder(){return Cylinder;}; // o
  bool sphere(){return Sphere;}; // o

  void CalOmegaInFixedSource();

  // For power iteration
  real CalFluxGeneral(int grp,real epsif,int iter);
  void SetInitialFlux();
  void DoAcceleration(int iter,real errf,real fiss){};
  // For DSA of SNT
  void CalFluxAutoConv(int grp,real epsif=1e-7);

  void AllVectorClear();
};

#endif
