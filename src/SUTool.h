#ifndef SUTOOL
#define SUTOOL

#include "SensitivityData.h"
#include "LibData.h"
#include "GroupData.h"
#include "Numeric.h"

class SUTool{
 protected:
  vector<SensitivityData> sns_data;
  vector<string> sns_name;
 public:
  SUTool(){};
  int GetSensitivityIndex(string name);
  void ReadSensitivityDataFromFile(string mdir,string filename,string data_name);
  void CalLibraryEffect(string data_name,XSLibrary &lib1,XSLibrary &lib2,bool grp_wise_prt=false);
  SensitivityData &GetSensitivityData(string name);
  enum xstype GetXSType(int mt);
};

#endif
