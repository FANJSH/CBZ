#ifndef PJISLABPIJ
#define PJISLABPIJ

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Numeric.h"

using namespace std;

class PJISlabPij{
  int mesh;
  vector<int> regid;
  vector<real> width;
 public:
  PJISlabPij(){};
  void PutMesh(int i);
  void PutRegionID(int *i);
  void PutRegionID(vector<int> i){regid=i;};
  void PutWidth(real *i);
  void PutWidth(vector<real> i){width=i;};
  void CalculationPij(real *xsinp,real *pij,bool aniso=false,bool black=false);
  void CalculationPijBlack(real *xsinp,real *pij,bool aniso=false){CalculationPij(xsinp,pij,aniso,true);};
  int GetMesh(){return mesh;};
  real GetWidth(int i){return width[i];};
};

#endif
