#ifndef SENSITIVITIES
#define SENSITIVITIES

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Numeric.h"
#include "GroupData.h"
#include "UNC_Parameters.h"
#include "UNC_CrossSection.h"
#include "SensitivityData.h"

using namespace std;

class SensitivityContainer{
  ParameterList* plist;
  vector<SensitivityData> sens;
 public:
  SensitivityContainer(ParameterList* plist);
  ParameterList *GetParameterList(){return plist;};
  void ReadSensitivityDataFromFile(string dirname,string filename,string core,string chara,int step=0);
  real &GetSensitivity0D(int mat,int mt,string core,string chara,int step=0);
  GroupData1D &GetSensitivity1D(int mat,int mt,string core,string chara,int step=0);
  GroupData2D &GetSensitivity2D(int mat,int mt,string core,string chara,int step=0);
  SensitivityData &GetSensitivityData(string core,string chara,int step=0);
  SensitivityData &GetSensitivityData(int i){return sens[i];};
  void SetRelativeSensitivity(Library &lib);
  void Produce1DSensitivity();
  void GetDiscreteInelasticSensitivity(string mdir, string filename, int matid);
  void TransformToConstrainedSensitivityForChi(Library &lib);
  int GetSize(){return sens.size();};
  void ChangeMTNumber(int mat,int mt1,int mt2);
};


#endif


