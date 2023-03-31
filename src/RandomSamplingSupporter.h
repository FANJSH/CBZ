#ifndef RANDOMSAMPLINGSUPPORTER
#define RANDOMSAMPLINGSUPPORTER

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "LibData.h"
#include "UNC_CrossSection.h"
#include "UNC.h"
#include "UNC_Covariance.h"
#include "UNC_Parameters.h"
#include "UNC_Sensitivity.h"
#include "SensitivityData.h"

using namespace std;

class RandomSamplingSupporter{
 private:
  XSLibrary xslib;
  LibraryCovariance xscov;
  int rs_nuc_data;
  vector<int> rs_matid;
  vector<int> rs_mtid;
  int sample;
  int grp;
  vector< vector<GroupData1D> > sampled_data;
  vector<GroupData1D> sampled_avg;
  vector< vector<real> > orgxs;
  int outp_num;
  vector< vector<real> > outp_data;
  vector<real> output_avg;
 public:
  RandomSamplingSupporter(){rs_nuc_data=0;};
  ~RandomSamplingSupporter(){};
  void ReadXSLibrary(int nucnum,string libdir,string *fname);
  void ReadCovarianceFromFile(int nucnum,string covdir,string *fname);
  void PutSampledNuclearDataInfo(int i1,int* i2, int *i3);
  void PutSampleNum(int i);
  XSLibrary &GetXSLib(){return xslib;};
  void PrepareSamplingData(int sample,real factor);
  void SampleAverageCalculation(real factor);
  void PrepareOriginalData();
  void PerturbXSLib(int id,real factor=1.0);
  void PutOutputNum(int i);
  void PutOutputData(int i,int j,real k);
  void PutOutputData(int i,real k);
  void ShowOutputStatics(int i);
  void CalOutputAverage();
  SensitivityData SensitivityCalculation(real factor);
  void OutputPartialUncertaintyCalculation();
};

#endif


