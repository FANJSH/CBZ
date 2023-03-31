#ifndef GAMMATOOL
#define GAMMATOOL

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <fstream>
#include <map>
#include <string>
#include "Numeric.h"
#include "GammaLibData.h"
#include "Medium.h"
#include "GeneralSystem.h"
#include "Burnup.h"

using namespace std;

class DecayGammaSpectrum{
 private:
  int ggroup;
  GroupData1D ebnd;
  map<int,GroupData1D> data;
 public:
  DecayGammaSpectrum(){};
  void PutGgroup(int i);
  void ReadFile(string cbglibdir, string filename);
  int GetGgroup(){return ggroup;};
  bool ExistData(int mat);
  map<int,GroupData1D> &GetData(){return data;};
  GroupData1D &GetData(int mat);
  GroupData1D &GetEbnd(){return ebnd;};
  void ShowGammaSpectrum(GroupData1D &spec);
};

class GammaTool{
 public:
  GammaTool(){};
  void GammaXSOUT(GammaXSLibrary &gxslib,Medium &med,string mdir,string ss,bool print=true);
  GroupData2D ReadYieldDataFromFile(string mdir,string filename);
  void ReadDataFromFile(string mdir,string filename,Medium &med,GroupData2D &yld);
  void HeatingCalculation(GeneralSystem &sys,GeneralSystem &gam,Burnup &bu,int *reg_med,bool decay=false);
  void GammaSourceCalculation(GeneralSystem &sys,GeneralSystem &gam,vector<GroupData2D> &yld,int *reg_med);
  void AddGammaSourceFromDecay(GeneralSystem &sys,GeneralSystem &gam,map<int,GroupData1D> &spec,Burnup &bu);
  void GammaTransportCalculation(GeneralSystem &gam);
  void SetGammaEnergyAsProductionXS(GeneralSystem &gam);
};

#endif


