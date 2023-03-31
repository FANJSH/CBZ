#ifndef GENERALBURNER
#define GENERALBURNER

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

#include "LibData.h"
#include "MATIDTranslator.h"
#include "OnePointCalculator.h"
#include "SelfShieldingCalculator.h"
#include "Medium.h"
#include "FRDesignTool.h"


/*
#include "Numeric.h"
#include "Medium.h"
#include "OnePointCalculator.h"
#include "IrregularGeometryInformation.h"
#include "PJI_trajectoryset.h"
#include "PJI_system.h"
#include "Burnup.h"
#include "ENDFIDTranslator.h"
#include "MATIDTranslator.h"
#include "BurnupChainGenerator.h"
#include "FRDesignTool.h"
#include "DelayedNeutronData.h"
#include "SensitivityData.h"
#include "SelfShieldingCalculator.h"
*/

using namespace std;

class GeneralBurner{
 protected:
  int group;
  XSLibrary xslib;
  bool up_scattering;

  MATIDTranslator midt;
  OnePointCalculator opc;
  SelfShieldingCalculator ssc;
  // +++ 

  // +++ burnup condition
  bool burn_time_calc;
  int burn_step;
  int sub_step_org;
  vector<real> power_density_list;
  vector<real> flux_level_list;
  vector<real> burn_time;
  vector<real> burn_time_gwd;
  vector<int> sub_step_list;
  bool burn_time_GWd_t;
  bool burn_time_accumulate;
  bool input_flux_level;
  string input_power_unit;
  vector<real> keff;
  vector<real> abs_frac; // absorption reaction rate fraction of fuel media to total
  // +++
  TrajectorySet sys, sys_f;
  int mednum;
  int mednum_fuel;
  int med_clad;
  int med_water;
  vector<Medium> med;
  real hm_weight_init; // [g/cm] (NOT [g/cm3])
  bool pl0_calc;
  // +++
  bool sensitivity;
 public:
  GeneralBurner();
  ~GeneralBurner(){};

  void PutGroup(int i){group=i;};

  void SensitivityOff(){sensitivity=false;};

  // +++ Multigroup library treatment
  void SetLibrary(string cbglibdir,string lib,int fp_num,string *fp_nuc);  
  void SetLibrary(string cbglibdir,string lib="jendl-4.0");  
  void SetLibraryTest(string cbglibdir,string lib="jendl-4.0");  
  void SetLibrary_172g(string cbglibdir,string lib="jendl-4.0");
  void SetLibrary_361g(string cbglibdir,string lib="jendl-4.0");  
  void SetLibrary(XSLibrary& libinp){up_scattering=true; PutGroup(107); xslib=libinp;};
  void SetLibraryFP103(string cbglibdir,string lib="jendl-4.0");
  void SetLibraryAll(string cbglibdir,string lib="jendl-4.0");
  void SetLibraryForThorium(string cbglibdir,string lib="jendl-4.0");
  void SetLibraryForThorium_361g(string cbglibdir,string lib="jendl-4.0");  
  void AddLibraryForVariousCoolant(string cbglibdir,string lib="jendl-4.0");
  void AddLibraryForMinorFP(string cbglibdir);
  void SetLibrary70g(string cbglibdir);
  void SetLibraryFP103_70g(string cbglibdir);
  XSLibrary &GetXSLibrary(){return xslib;};

  // +++ Burnup condition
  void PutBurnStep(int i);
  void PutSubstep(int i){sub_step_org=i;};
  void PutSubstepList(int *inp);
  void PutPowerDensityList(real *inp,string type="W_cm"); // [W_cm] or [MW_t]
  void PutPowerDensity(real inp,string type="W_cm");
  void PutFluxLevelList(real *inp);
  void PutFluxLevel(real inp);
  void PutBurnTime(real *inp,bool GWd_t=true,bool accumulate=true,int div=1);
  void PutBurnTime(real inp,bool GWd_t=true,bool accumulate=true);
  void PreCalculation_bt();
  void ShowBurnupConditionSetting();

  // +++ Medium data
  void PutFuelData(FRDTFuelComposition &fcom,real temp);
  void PutFuelData(int nuc,int *mat,real *den,real temp,int medid=0);
  void PutNonfuelData(int medid, int nuc,int *mat,real *den,real temp);
  void PutCladData(int nuc,int *mat,real *den,real temp);
  void PutModeratorData(int nuc,int *mat,real *den,real temp);
  Medium& GetMedium(int i){return med[i];};
  void PutHeavyMetalInitialWeight(real i){hm_weight_init=i;};
  real GetHeavyMetalInitialWeight(){return hm_weight_init;};
  void PutTemperature(int medid, real temp);
  void PutFuelTemperature(real temp);
  void PutModeratorTemperature(real temp){PutTemperature(med_water, temp);};    

  // +++ Calculation condition
  void PutWhiteBoundary();
  
  // +++ Access function
  real GetKeff(int st);

};

#endif


