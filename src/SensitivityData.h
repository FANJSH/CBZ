#ifndef SENSITIVITIYDATA
#define SENSITIVITIYDATA

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

class SensitivityData{
 private:
  string name1; // core/experimental name
  string name2; // parameter name
  string name3; // used library
  real val;  // value of parameter
  int group; // number of energy groups
  GroupData1D enband; // energy boundary in group structure
  vector<real> sens0d;
  vector<GroupData1D> sens1d;
  vector<GroupData2D> sens2d;
  vector<int> mat_list0d;
  vector<int> mt_list0d;
  vector<int> mat_list1d;
  vector<int> mt_list1d;
  vector<int> mat_list2d;
  vector<int> mt_list2d;

 public:
  SensitivityData();

  void PutName(string i1,string i2,string i3){name1=i1; name2=i2; name3=i3;};
  void PutValue(real i){val=i;};
  void PutGroup(int i);
  void PutEnergyBoundaryData(int i,real j);

  GroupData1D &GetEnband(){return enband;};
  int GetGroup(){return group;};
  real GetValue(){return val;};
  string GetName1(){return name1;};
  string GetName2(){return name2;};
  string GetName3(){return name3;};

  int FindData0D(int mat,int mt);
  int FindData1D(int mat,int mt);
  int FindData2D(int mat,int mt);
  real& GetSensitivity0D(int mat,int mt);
  real& GetSensitivity0D(int i){return sens0d[i];};
  GroupData1D &GetSensitivity1D(int mat,int mt);
  GroupData1D &GetSensitivity1D(int i){return sens1d[i];};
  GroupData2D &GetSensitivity2D(int mat,int mt);
  GroupData2D &GetSensitivity2D(int i){return sens2d[i];};
  void PutSensitivity0D(int mat,int mt,real dat);
  void PutSensitivity1D(int mat,int mt,GroupData1D dat);
  void PutSensitivity2D(int mat,int mt,GroupData2D dat);
  int GetMatList0D(int i){return mat_list0d[i];};
  int GetMtList0D(int i){return mt_list0d[i];};
  int GetMatList1D(int i){return mat_list1d[i];};
  int GetMtList1D(int i){return mt_list1d[i];};
  int GetMatList2D(int i){return mat_list2d[i];};
  int GetMtList2D(int i){return mt_list2d[i];};
  int GetSize0D(){int tmp=sens0d.size(); return tmp;};
  int GetSize1D(){int tmp=sens1d.size(); return tmp;};
  int GetSize2D(){int tmp=sens2d.size(); return tmp;};

  void AddSensitivityData(SensitivityData &sens2,bool replace=false);
  void AverageSensitivityData(vector<SensitivityData> &sens);
  void WithdrawSensitivityData(SensitivityData &sens2);
  SensitivityData CalSensitivityDataForReactivity(SensitivityData &sns2);
  SensitivityData CalSensitivityDataForBeta(SensitivityData &sns2);
  void CalSensitivity1DFrom2D(int mat,int mt);

  void ShowSensitivity1D(int mat,int mt);
  void ShowSensitivity1DExcel(int mat,int mt);
  void ShowSelf(real minumum=0.);
  void ShowSelf(MATIDTranslator &midt,real minumum=0.);

  real GetFPYieldSensitivityFissileWise(int fisid);
  real GetFPYieldSensitivityFPWise(int fpid);
  void ShowFPYieldSummary(MATIDTranslator &midt, real minimum=0.);


  void WriteFile(string mdir,string filename);
  void ReadFile(string mdir,string filename);
  void ReadFileOldFormat(string mdir,string filename);

  void ChangeMatMt(int mato,int mto,int matn,int mtn);

  void TransformToConstrainedSensitivity(int mat,int mt,GroupData1D &xs);
  void Factorize(real factor);
  void Diff(SensitivityData &sns2,real dif_threshold=0.01,real sns_threshold=1e-5);

  SensitivityData EnergyGroupTranslation(GroupData1D& ebnd);
  SensitivityData Compressing(int num, string *nuclist);

  void CompareAnotherSensitivityData(SensitivityData &sd2);
  real CalEuclideanNormForAbsorptionXS();

};

#endif


