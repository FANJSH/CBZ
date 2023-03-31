#ifndef DELAYEDNEUTRONDATA
#define DELAYEDNEUTRONDATA

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Numeric.h"
#include "GroupData.h"
#include "MATIDTranslator.h"

using namespace std;

class DelayedNeutronData{
 private:
  int nucnum,famnum;
  vector<int> nucid; // Mat No. in JFA format
  vector< vector<GroupData1D> > chi;
  vector<GroupData1D> yield;
  vector< vector<real> > fraction; 
  vector< vector<real> > lambda;
  MATIDTranslator midt;
 public:
  DelayedNeutronData();
  ~DelayedNeutronData(){};
  void PutFamilyNumber(int famin);
  int GetNumFromNucid(int nucid);
  void PutDataFromFile(string mdir,string ss);
  int GetFamnum(){return famnum;};
  int GetNucnum(){return nucnum;};
  int GetNucid(int i){return nucid[i];};
  GroupData1D &GetChi(int i,int j){return chi[i][j];};
  GroupData1D &GetKai(int i,int j){return GetChi(i,j);};
  GroupData1D GetAveragedChi(int nucid);
  GroupData1D GetAveragedKai(int nucid){return GetAveragedChi(nucid);};
  real GetFraction(int i,int j){return fraction[i][j];};
  real GetLambda(int i,int j){return lambda[i][j];};
  GroupData1D &GetYield(int i){return yield[i];};
  GroupData1D &GetYieldFromMatNumber(int i){return yield[GetNumFromNucid(i)];};
  // 
  void ShiftChi(GroupData1D &ebnd, real factor);
  void SetConstantYield(int g_set);
  void ShowSelf();
};

#endif


