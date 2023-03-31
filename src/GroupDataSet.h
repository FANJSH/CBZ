#ifndef GROUPDATASET
#define GROUPDATASET

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "Numeric.h"
#include "GroupData.h"

using namespace std;

class GroupDataSet{
 private:
  string name;
  int num1d, num2d, grp;
  vector<int> dim1d;
  vector<int> dim2d;
  map<enum xstype,int> map1d;
  map<enum xstype,int> map2d;
  vector< vector<GroupData1D> > data1d;
  vector< vector<GroupData2D> > data2d;
 public:
  GroupDataSet(){grp=-1; num1d=0; num2d=0;};
  GroupDataSet(string naminp){grp=-1; Init(naminp);};
  ~GroupDataSet();
  void Init(string naminp);
  void PutGrp(int i);
  void PutNum1d(int i);
  void PutNum2d(int i);
  void PutDim1d(int i,int j);
  void PutDim2d(int i,int j);
  void PutDim1d(enum xstype name,int i);
  void PutDim2d(enum xstype name,int i);
  int GetNum1d(){return num1d;};
  int GetNum2d(){return num2d;};
  int GetDim1d(int i){return dim1d[i];};
  int GetDim2d(int i){return dim2d[i];};
  int GetDim1d(enum xstype name){return dim1d[GetPos1d(name)];};
  int GetDim2d(enum xstype name){return dim2d[GetPos2d(name)];};
  int GetGrp(){return grp;};
  string GetName(){return name;};
  GroupData1D &GetData1d(enum xstype name, int l=0);
  GroupData2D &GetData2d(enum xstype name, int l=0);
  GroupData1D &GetData1d(int i,int l);
  GroupData2D &GetData2d(int i,int l);
  GroupDataSet Cond(int ngrp,int *bgrp,GroupData1D fl,GroupData1D cu);
  void ShowSelf();
  void ShowData1D(GroupData1D &enband);
  int GetPos1d(enum xstype name);
  int GetPos2d(enum xstype name);
  bool CheckSameType(GroupDataSet &sec);
  void DataCopy(GroupDataSet &sec);
  void DataCopyPL(GroupDataSet &sec,int plmax_2d);
  void DataCopy(GroupDataSet &sec,int lwgrp);
  void DataCopy(GroupDataSet &sec,int upgrp,int lwgrp);
  void MultiplyExceptForChi(real factor,int g);
  void Multiply(real factor,int g);
  void Multiply(real factor);
  void Add(GroupDataSet &sec);
  void AllVectorClear();
  void AllVectorClear2DData();

  void XSMatrixNormalizationTo1DXS();
  void P1MatrixNormalizationToMubar();
  // fast access to macroscopic cross sections
  GroupData2D &GetSigs(int l=0){return data2d[0][l];};
  GroupData1D &GetSigt(int l=0){return data1d[0][l];};
  GroupData1D &GetSigtr(int l=0){return data1d[3][l];};
  GroupData1D &GetNusigf(){return data1d[1][0];};
  GroupData1D &GetChi(){return data1d[2][0];};
  GroupData1D &GetKai(){return data1d[2][0];};
  GroupData1D &GetSiga(){return data1d[5][0];};
  GroupData1D &GetD(){return data1d[6][0];};
  GroupData1D &GetDr(){return data1d[7][0];};
  GroupData1D &GetDz(){return data1d[8][0];};
  GroupData1D &GetSigr(){return data1d[4][0];};
  GroupData1D &GetSign2n(){return data1d[9][0];};
  // fast access to microscopic cross sections
  GroupData1D &GetMicSigf(){return data1d[0][0];};
  GroupData1D &GetMicSigc(){return data1d[1][0];};
  GroupData1D &GetMicNu(){return data1d[2][0];};
  GroupData1D &GetMicChi(){return data1d[5][0];};
  GroupData1D GetMicNusigf(){return data1d[2][0].mult(data1d[0][0]);};
  GroupData2D &GetMicSigel(int l){return data2d[0][l];};
  GroupData2D &GetMicSiginel(int l){return data2d[1][l];};
  GroupData2D &GetMicSign2n(int l){return data2d[2][l];};

  // For NFI2020 joint research
  // Implemented by Yanagihara-kun
  // Source was sent from Yanagohara-kun in 2021/12/26
  GroupDataSet Caldsdx1(GroupDataSet &a, GroupDataSet &b, real c);
  GroupDataSet Caldsdx2(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d, bool nu=false);
  GroupDataSet Caldsdx3(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e);
  GroupDataSet Caldsdx4(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, real *f);
  GroupDataSet Caldsdx5(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, real *g);
  vector<GroupDataSet> Caldsdx1Quadr(GroupDataSet &a, GroupDataSet &b, real c);
  vector<GroupDataSet> Caldsdx2Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d);
  vector<GroupDataSet> Caldsdx3Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e);
  vector<GroupDataSet> Caldsdx5Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, real *g);
  vector<GroupDataSet> Caldsdx6Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, GroupDataSet &g, real *h);
  vector<GroupDataSet> Caldsdx7Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, GroupDataSet &g, GroupDataSet &h, real *i);
  //vector<GroupDataSet> CaldsdxQuadr(int num, vector<GroupDataSet> &a, real *b);
  vector<GroupDataSet> Caldsdx2Cubic(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d);
  vector<GroupDataSet> Caldsdx3Cubic(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e);
  GroupDataSet CaldeltaXS(GroupDataSet &baseXS, GroupDataSet &coeff, real delta);
  GroupDataSet CaldeltaXSQuadr(GroupDataSet &baseXS, vector<GroupDataSet> &coeff, real delta);
  GroupDataSet CaldeltaXSCubic(GroupDataSet &baseXS, vector<GroupDataSet> &coeff, real delta);
};

#endif


