#ifndef RESONANCACALCULATION_FAN
#define RESONANCACALCULATION_FAN

#include "Medium.h"
#include "Numeric.h"
//#include "MATIDTranslator.h"
#include "iomanip"
#include "vector"
using namespace std;

class ResonanceCalculationInHomoSystemFan{
  // Class: CamelCase;
  // parameter & array: snake_case;
  // function: lowerCamelCase;
  public:
  // Accessable parameters and functions. 
  int group,pl,nucnum;
  real temperature;
  string libdir;
  XSLibrary xslib;
  Medium this_medium;
  vector<string> NuclideNameSaving;
  vector<int> NuclideIDSaving;
  vector<real> NumberDensitySaving;
  vector<GroupData1D> BackgroundXS;
  vector<GroupData1D> FluxWeightedTotalXSSaving;
  vector<vector<GroupData1D>> EffectiveXS1D;
  vector<vector<real>> EffectiveTotalXSSaving;
  vector<vector<GroupData2D>> EffectiveN2NMatrixPL;
 
  /*
  vector<GroupData1D> MicroFissionXS;
  vector<GroupData2D> MicroAbsorptionXS;
  vector<GroupData2D> MicroTotalXS;
  vector<GroupData2D> MicroScatteringXS;
  */
  void WithSelfShieldingCalculation(int arg_grp, int arg_pl, int arg_nucnum, string arg_libdir, string *arg_filename, int *arg_nuclideID, real *arg_numberdensity, real arg_temperature,bool sig0state);
  int testvalue;
  void AssignInfiniteDilutionCrossSection(int arg_grp, int arg_pl, int arg_nucnum, string arg_libdir, string *arg_filename, int *arg_nuclideID, real *arg_numberdensity, real arg_temperature);

};



#endif 