#ifndef OKTAVIANTOOL
#define OKTAVIANTOOL

#include "GroupData.h"
#include "Numeric.h"
#include "LibData.h"
#include "Medium.h"
#include "SensitivityData.h"

class OktavianTool{
 protected:
  int grp;
  vector<GroupData1D> sns1d;
  vector<int> sns1d_mat;
  vector<int> sns1d_mt;
  vector<GroupData2D> sns2d;
  vector<int> sns2d_mat;
  vector<int> sns2d_mt;
  SensitivityData sens;
 public:
  OktavianTool(int ginp);
  void PutGroup(int ginp);
  // (sensitivity data)
  void ReadSensFile(string mdir,string filename);
  void PutSensitivityData(SensitivityData &snsin){sens=snsin;};
  void ShowSensData();
  GroupData1D &GetSens1D(int mati,int mti);
  GroupData2D &GetSens2D(int mati,int mti);
  // (library effect calculation)
  void LibraryEffectCal(XSLibrary &lib1,XSLibrary &lib2,int mati,bool pring_eng=false);
  void LibraryEffectCal(Medium &med1,Medium &med2,int mati,bool print_eng=false);
  void LibraryEffectCalInelastic(string lib1,string lib2,string fname,int mat);
  void ShowSecondaryDistribution(string lib1,string fname,int srcg,GroupData1D& ebnd);
};

#endif
