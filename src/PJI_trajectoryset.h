#ifndef TRAJECTORYSET
#define TRAJECTORYSET

#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "FunctionTable.h"
#include "PJI_trajectory.h"
#include "GeomVector.h"
#include "GeomPlane.h"
#include "Numeric.h"
#include "IrregularGeometryInformation.h"

using namespace std;

enum BCondition{Black,White,Periodic,Reflective};

class TrajectorySet{
 protected:
  int NumTrajectory, regnum;

  int angular_division;
  vector<int> trajectory_on_angular_boundary;
  vector<real> discreted_angle;
  vector<GeomVector> vdir; 

  real surface;
  fktab fkk;
  exptab etab;
  vector<Trajectory> line;
  vector<int> NextTrajectory;
  vector<int> PreTrajectory;
  vector<real> volume;
  BCondition BoundaryCondition;

  bool cal_vol_sph;
  bool angle360;
 public:
  TrajectorySet();
  void PutXS(real *xsinp);
  void CalculationPij(real *xsinp,real *pijinp,bool aniso=false);
  void CalculationPijSphere(real *xsinp,real *pijinp);
  void CalVolume();
  void CalVol(){CalVolume();};
  void CalVolumeSphere(int r,int *rid,real *rinp);
  real GetVol(int i){return volume[i];};
  real GetVolume(int i){return volume[i];};
  void Cal_White(real *xsinp,real *pijinp,bool black=false,bool aniso=false);
  //black=1 : black boundary
  void Cal_Periodic(real *xsinp,real *pijinp,bool aniso=false);
  //aniso=1 : Anisotropic pij
  void CalTrajectory(IrregularGeometryInformation &ginp,int AngleDiv,real dl,real angle=90.);
  void CalTrajectoryMB(IrregularGeometryInformation &ginp,int AngleDiv,real dl,real angle=90.);
  Trajectory PlaneToTrajectory(GeomPlane &inp, GeomVector& p1, GeomVector& p2, IrregularGeometryInformation &gset);
  // instance of 'GeomPlane' -> instance of 'Trajectory'
  void PutBoundaryCondition(BCondition i);
  BCondition GetBoundaryCondition(){return BoundaryCondition;};
  bool IsFixedBoundaryCondition();
  bool IsAngle360(){return angle360;};
  int GetRegnum(){return regnum;};
  int GetNextTrajectory(int i){return NextTrajectory[i];};
  int GetPreTrajectory(int i){return PreTrajectory[i];};

  Trajectory &GetTrajectory(int i){return line[i];};
  void PutAngularDivision(int i);
  int GetAngularDivision(){return angular_division;};
  int GetTrajectoryOnAngularBoundary(int i){return trajectory_on_angular_boundary[i];};
  real GetDiscretedAngle(int i){return discreted_angle[i];};
  int GetNumTrajectory(){return NumTrajectory;};
  // I/O
  void WriteFile(string mdir,string ss);
  void ReadFile(string mdir,string ss);

  void LengthFactorize(real factor);
  int GetMaximumSegments();
};

#endif
