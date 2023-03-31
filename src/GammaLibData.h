#ifndef GAMMALIBDATA
#define GAMMALIBDATA

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <fstream>
#include <map>
#include <string>
#include "Numeric.h"
#include "GroupData.h"
#include "GroupDataSet.h"
#include "MATIDTranslator.h"

using namespace std;

class GammaLibData{
  int mat;
  int ngroup;
  int ggroup;
  vector<GroupData2D> yield; // (n,2n), (n,f), (n,c), (n,n'), (n,non-e)
  GroupDataSet xsdata;
 public:
  GammaLibData();
  void PutNGroup(int i);
  void PutGGroup(int i);
  int GetNGroup(){return ngroup;};
  int GetGGroup(){return ggroup;};
  GroupData2D &GetYield(int i){return yield[i];};
  void ReadFile(string mdir,string ss,MATIDTranslator &midt);
  void WriteFile(string mdir,string ss);
  GroupDataSet &GetXSData(){return xsdata;};
  void ShowSelf();
  int GetMat(){return mat;};
};

class GammaXSLibrary{
  map<int,GammaLibData> data;
  GroupData1D enband;
  GroupData1D wtflux;
  MATIDTranslator midt;
 public:
  GammaXSLibrary(string dir,string ss){ReadGEnergy(dir,ss);};
  void ReadGEnergy(string dir,string ss);
  GammaLibData& GetGammaLibData(int mat);
  bool ExistGammaLibData(int mat);
  GroupData1D &GetEnband(){return enband;};
  GroupData1D &GetWtflux(){return wtflux;};
  int GetGroup(){return enband.get_imax()-1;};
  void ReadFile(int nucnum,string mdir,string *filename);
  void ReadFile(int nucnum,string mdir,string *filename,int *nucid);  
  void AddGammaLibData(int mat,GammaLibData &libinp){data[mat]=libinp;};
  void AddDelayedFissionGamma(string mdir,string fname,real fraction);
};


#endif


