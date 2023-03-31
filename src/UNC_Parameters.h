#ifndef UNCPARAMETERS
#define UNCPARAMETERS

#include <fstream>
#include <list>
#include <vector>
#include "GroupData.h"
#include "UNC_Covariance.h"

using namespace std;

class ParameterList
{
  vector<string> core_tag;
  vector<string> chara_tag;
  vector<int> step_tag;
 public:
  ParameterList();
  int GetSize();
  int FindData(string core,string chara,int step=0);
  void ShowSelf();
  void AddNewData(string core,string chara,int step=0);
  ParameterList operator+(ParameterList &other_list);
  string GetCoreTag(int i){return core_tag[i];};
  string GetCharaTag(int i){return chara_tag[i];};
  int GetStepTag(int i){return step_tag[i];};
};

class Parameters{
  ParameterList *plist;
  GroupData1D val;
 public:
  Parameters(ParameterList* plist);
  void PutValue(string core,string chara,real value,int step=0);
  void PutValue(int i,real value){val.put_data(i,value);};
  ParameterList *GetParameterList(){return plist;};
  real GetValue(int i){return val.get_dat(i);};
  GroupData1D &GetValue(){return val;};
  void ShowSelf();
};

class ParameterCovariance:public Covariance{
  ParameterList* plist;
 public:
  ParameterCovariance(Parameters& paras);
  ParameterList *GetParameterList(){return plist;};
  void PutValue(string core,string chara,real value,int step=0);
  void PutStandardDeviation(string core,string chara,real value,int step=0,string type="Relative");
  real GetStandardDeviation(string core,string chara,int step=0,string type="Relative");
  GroupData1D GetStandardDeviation(string type="Relative"){return Covariance::GetStandardDeviation(type);};
  void PutCorrelation(string core1,string chara1,string core2,string chara2,real cor);
  void PutCorrelation(string core,string chara,real cor);
  void PutCorrelation(real cor);
  real GetCorrelation(string core1,string chara1,string core2,string chara2,int step1=0,int step2=0);
};

class ParametersContainer{
  vector<Parameters> paras;
  vector<string> name;
 public:
  ParametersContainer(){};
  int GetSize(){int tmp=paras.size(); return tmp;};
  void PutParameters(Parameters inp,string nameinp);
  int FindData(string nameinp);
  Parameters &GetParameters(string nameinp);
  Parameters &GetParameters(int i){return paras[i];};
  string GetName(int i){return name[i];};
};


#endif
