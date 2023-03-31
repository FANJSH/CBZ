#ifndef COVARIANCE
#define COVARIANCE

#include <fstream>
#include <iostream>
#include <cstdio>

#include "GroupData.h"

class Covariance{
 protected:
  int size1,size2;
  vector<GroupData1D> val;
  GroupData2D cov; // (absolute value)
  bool Variance;
 public:
  Covariance(){};
  Covariance(string type){PutType(type);};
  void PutType(string type="variance");
  void PutSize(int inp);
  void PutSize(int inp1, int inp2);
  void PutVal(real *inp);
  void PutVal(int i,real *inp);
  void PutVal(int i,GroupData1D &inp);
  void PutVal(int i,real inp);
  void PutCov(real *inp, string type="Absolute");
  void PutCov(real *dev, real *corr, string type="Absolute");
  void PutCov(real *dev, GroupData2D corr, string type="Absolute");
  void PutCov(GroupData2D cov, string type="Absolute");
  void PutStandardDeviation(int i,real inp,string type="Relative");
  void CheckPosition(int i,int j=-1);
  int GetSize(){return size1;};
  int GetSize1(){return size1;};
  int GetSize2(){return size2;};
  GroupData2D GetCovariance(string type="Relative");
  GroupData2D &GetCov(){return cov;};
  GroupData1D GetStandardDeviation(string type="Relative");
  GroupData2D GetCorrelationMatrix();
  GroupData1D GetVal(int i=0){return val[i];};
  void Multiply(real v){Factorize(v);};
  void Factorize(real v);
  void Factorize(vector<real>& v);
  void Factorize(vector<real>& v1, vector<real>& v2);
  void SetNoCorrelation();
  //
  void PutGrp(int i){PutSize(i);};
  void PutGrp(int i,int j){PutSize(i,j);};
  int GetGrp(){return GetSize();};
  void NormalizeWeightFunction(GroupData1D &wgt,int bgrp,int* bndgrp);
  Covariance Cond(GroupData1D &wgt,int bgrp,int* bndgrp);
  Covariance CondRelative(GroupData1D &wgt,int bgrp,int* bndgrp);
  Covariance GetCovarianceLognormal(real factor=1.);
  bool GetVariance(){return Variance;};
  // 
  void DoRandomSampling(int num, vector<GroupData1D> &sample,string type);
};

#endif
