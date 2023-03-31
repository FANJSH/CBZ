#ifndef PNXSYSTEM
#define PNXSYSTEM

#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>

#include "GroupData.h"
#include "Medium.h"
#include "SphericalHarmonics.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"
#include "GeneralMesh.h"
#include "GeneralSystem.h"

using namespace std;

class PNXSystem:public GeneralSystem{
  int BC[2];
  int plexp,plexp1;
  int sz;
  vector< vector<real> > pl_pl; // for Marshak BC
  vector<int> pos;
  vector< vector<real> > amatg;
 public:
  PNXSystem(int Ndiminp,int grp,int MAX_MED=20);
  ~PNXSystem();

  void PutPLexp(int i);
  void PutCartMeshInfo(CartMeshInfo inp);
  void SetArray();
  void SetInitialFlux();
  real CalFluxGeneral(int ginp,real cin,int iter);
  real CalFluxReflectiveReflective(int ginp);
  real CalFluxReflectiveVacuum(int ginp);
  real CalFluxReflectiveVacuumCullen(int ginp);
  real CalFluxGeneral2(int ginp,real cin,int iter);
  real GetAngularFlux(int m,int g,real mu,real pl=-1);
  void CalAmat();

  void ShowAngularFlux(int m,int g);
};

#endif
