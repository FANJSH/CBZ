#ifndef FMG_TOOL
#define FMG_TOOL

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "GroupData.h"
#include "MATIDTranslator.h"

using namespace std;

class FMG_Tool{
 private:
  string file_grp;
  string dir_acef;
  bool slowing_down_calc;
  real sd_wgt, fehi;
  int grp;
  vector<float> ebnd;

  int num_bgxs;
  vector<real> bgxs;

  int pnt;
  vector<real> epnt;
  vector<real> wgt;
  bool wgt_inverse_e;

  int nucnum;
  vector<string> aceid;
  vector<real> den;
  vector<string> nucname;
  
 public:
  FMG_Tool();
  void PutWeightFunctionEnergyInverse(){wgt_inverse_e=true;};
  void SetBackgroundCrossSection(int numbg, real *xsinp);
  void SetGroupStructureFile(string inp);
  void SetAceFileDirectory(string inp);
  void SlowingDownCalculationOn(real wgt, real fehi_inp=1e4);
  void ReadEnergyGroupStructure();
  void ReadWeightFunction(string inp);
  void PutMaterialInfo(int i, string *sinp, real *dinp, string *nameinp);
  void input_generator_lanlmix(string filenameadd);
  void input_generator(string aceid, string nucname, string temp);  
  void input_generator_tsl(string aceid1,string aceid2, string nuc,string acefname,string sab_type);  
  void delayed_neutron_data_reader(int famnum, int group, string nucname, string frendy_dir, string cbglib_dir);
  void nenergy_generator(int group, string frendy_dir, string cbglib_dir, string nucname);    
};


#endif


