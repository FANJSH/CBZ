#ifndef NUCLIDE
#define NUCLIDE

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include "Numeric.h"
#include "GroupData.h"
#include "GroupDataSet.h"
#include "LibData.h"

using namespace std;

class Nuclide{
 private:
  int grp;
  int matnum;   // Nuclide ID
  real density; // Nuclide number density
  real temperature;
  GroupDataSet Micxs; // Microscopic cross section
 public:
  Nuclide(bool simple=false);
  ~Nuclide();
  void PutGrp(int i);
  int GetGrp(){return grp;};
  void PutMatnum(int i);
  void PutDensity(real i){density=i;};
  void AddDensity(real i){density+=i;};
  void PutTemperature(real i){temperature=i;};
  int GetMatnum(){return matnum;};
  int GetID(){return GetMatnum();};
  GroupDataSet &GetMicxs(){return Micxs;};
  void PutMicxs(GroupDataSet micinp);
  real GetDensity(){return density;};
  real GetTemperature(){return temperature;};
  Nuclide Cond(int ngrp,int *bgrp,GroupData1D fl,GroupData1D cu);
  void CalMicroFromF(int g,LibData &lib,vector<real> &f,bool matrix=true);
  bool IsFissionable();
  void PrintMicroSectionTable(GroupData1D &enband);
};

#endif

