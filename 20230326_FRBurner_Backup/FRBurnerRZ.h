#ifndef FR_BURNER_RZ
#define FR_BURNER_RZ

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <chrono> // to logging the running time 

#include "PLOSE_system.h"
#include "PLOS_system.h"
#include "SNRZ_system.h"
#include "Medium.h"
#include "CartCore.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"
#include "Burnup.h"
#include "OnePointCalculator.h"
#include "SelfShieldingCalculator.h"
#include "FRDesignTool.h"
#include "MATIDTranslator.h"
#include "PJI_system.h"
#include "BurnupChainGenerator.h"


using namespace std;

class FRBurnerRZ{
 private:
  int dim, group;
  int pl, sn;
  int mednum, med_burn;
  vector<Medium> med;
  vector<Medium> med_voided;
  vector<Medium> med_doppler;
  //vector<Medium> med_high_temperature;
  vector<bool> med_data_inp;
  int burn_nuc;
  vector<int> matno;
  CartMeshInfo cmi;
  real power; // Wth
  real cycle_length; // day
  real refuel_day; // days for refueling
  real cooling_day; 
  int trace_cycle, cycle_div;
  int batch_ic, batch_oc, batch_rb, batch_max;
  vector< vector< vector<real> > > density; // [med][batch][nuc] (Batch-wise ND)
  vector< vector<real> > density_ave;       // [med][nuc] (Batch-averaged ND)
  vector< vector< vector<GroupData1D> > > fwd_nuc; // [cyc][step][med][nuc]
  vector< vector< vector< vector<real> > > > sigc_1g; // [med][cyc][step][nuc]
  vector< vector< vector< vector<real> > > > sigf_1g;
  vector< vector< vector< vector<real> > > > sign2n_1g;
  vector< vector<int> > ir_cycle; // number of irradiated cycles
  vector< vector<real> > ac_burn; // accumulated burnup
  vector< vector< vector<real> > > flux_level; // flux history [1/s/cm2]
  vector< vector<real> > keff;    
  vector<real> vol_per_med;
  vector<real> thm_per_med; // initial tHM for each medium [g]
  OnePointCalculator opc;
  MATIDTranslator midt;

  // material zone management
  int num_mzone;
  vector<int> mzone_per_med;
  vector<int> num_med_mzone; 
  vector< vector<real> > iden_mzone;
  vector<int> batch_mzone;
  vector<string> name_mzone;
  vector<int> mzone_rep; // representative medium for each mzone
  vector<int> num_pin_mzone;
  vector<real> area_assembly_mzone;

  bool void_cal;
  real void_ratio;
  bool dop_cal;
  real delta_t;
  bool xscalc_for_reactivity;
  real sodium_density_in_voided_case;
  real delta_t_in_doppler_case;

  bool print_linepower_map;
  bool cmfd_on;
  // For Ogawa calculations
  bool ogawa_cal;
  bool show_fission_info;
  vector<real> ogawa_year;
  int ogawa_mat1,ogawa_mat2;
  string dirname_boc, dirname_eoc;

  void PutSADataToMZoneHomogeneous(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib);
  void PutSADataToMZoneHeterogeneous1D(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib);
  void PutSADataToMZoneHeterogeneous2D(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib);
  void PutSADataToMZoneDefiningMediumData(Medium &med_ref, int mzone, int ii);
 public:
  FRBurnerRZ(int group_inp, int tnn, int* mninp);
  void PutCartMeshInfo(CartMeshInfo &cmii);
  Medium &GetMedium(int i){return med[i];};
  void CalMacroFromMicroBurnupMedium();

  void PutNumberofMedium(int i);
  void PutMedium(Medium &min,int id);
  void PutNonfuelInfo(int medid, int nucnum, int *nuc, real *den, real temp);
  void PutNonfuelInfo(int medid, FRDTComposition &inp, real temp);

  void AddActinideDecayHeatDataToBurnupChain(string cbglibdir,string fname,Burnup &bu);

  // material zone assignment
  void PutNumberMZone(int i);
  void PutMZoneInfo(int *mzone);
  void PutNameMZone(string *nameinp);
  void PutBatchMZone(int *batchinp);
  void PutMediumToMZone(Medium &min, int mzone);
  void PutSADataToMZone(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib, int model=0);
  void PutSADataToNonfuel(FRDTSubAssembly &sa, int medid);
  GroupData1D GetIntegratedFluxPerMZone(int i);


  void PutReactorPower(real i){power=i;};
  void PutTraceCycle(int i){trace_cycle=i;};
  void PutCycleDiv(int i){cycle_div=i;};
  void PutCycleLength(real i){cycle_length=i;};
  void PutRefuelDay(real i){refuel_day=i;};   // for ADS calculation
  void PutCoolingDay(real i){cooling_day=i;}; // for FR calculation
  void PutBatchNumber(int i,int j,int k);
  int GetBatchNumberFromMediumID(int medid);
  void PutBurnNucNum(int i){burn_nuc=i;};
  void PreCalculation(XSLibrary &xslib, Burnup &bu);
  void InitializeNumberDensity();
  void Run(XSLibrary &xslib, Burnup &bu);
  void CalRegionAveragedDensity();
  void StoreRegionAveragedDensity(int cyc,int step);
  //void CalNeutronFluxEnergySpectrum();
  void CMFDOff(){cmfd_on=false;};
  real MaximumPowerDensityCalculation(vector< vector<real> > &power_map, int &maxp_x,int &maxp_y);
  real CalTotalPower(Burnup &bu, vector< vector<real> > &pow_store, bool heat_by_nonfuel, bool heat_by_capture);
  void PrintPowerDensity(vector< vector<real> > &power_map);


  void VoidCalculationOn(real ratio=0.){void_cal=true; void_ratio=ratio;};
  void DopplerCalculationOn(real dt=10.){dop_cal=true; delta_t=dt;};
  void XSCalculationForReactivityOff(){xscalc_for_reactivity=false;};
  void PutSodiumDensityInVoidedCase(real den){sodium_density_in_voided_case=den;};
  void PutDeltaTInDopplerCase(real dt){delta_t_in_doppler_case=dt;};
  real CalVoidReactivity(int meds=-1,int mede=-1);
  real CalVoidReactivity(int vmednum, int *vmedlist);
  real CalDopplerReactivity(int meds=-1,int mede=-1);

  // +++ for ADS calculations
  void RunADS(XSLibrary &xslib, Burnup &bu,bool print=false,bool denout=true);
  // Removed in 2018/6/23
  //void RunADSSensitivityCalculation(XSLibrary &xslib, Burnup &bu,int target_med, int target_id, int cycle_id=-1, bool print=false);
  void ExternalSourceReading(vector< vector<GroupData1D> > &esrc);
  void ExternalSourceSetting(GeneralSystem &sys, vector< vector<GroupData1D> > &esrc);
  void PutSAGeometryData(FRDTSubAssemblyGeometry &sa_geom, int mzone);

  // [for Ogawa calculations]
  void OgawaCalculationOn(real y1, real y2, real y3, int mat1=96, int mat2=96);
  void ShowFissionInfoOn(){show_fission_info=true;};
  void ShowPlotExternalSource(string filename);
  real GetTonHeavyMetalPerMedium(int i){return thm_per_med[i];};
  void PutDirname(string a1,string a2){dirname_boc=a1; dirname_eoc=a2;};

  // +++ post-burnup calculation
  void PutNumberDensity(int cyc,int step);
  void CalNeutronFluxEnergySpectrum(bool adjoint=false,int meds=-1,int mede=-1);
  void Cal1groupXS();
  void ShowNeutronMultiplicationInfo();
  SensitivityData CalKeffSensitivity();
  SensitivityData CalKeffSensitivity(XSLibrary &xslib);  
  SensitivityData CalVoidReactivitySensitivity();  
  real CalDelayedNeutronParameters(DelayedNeutronData &dnd);
  void WriteFileMediumData(string dirname,string filename);
  void WriteFileNDData(string dirname,bool init=false,string addname="");
  vector<real> GetHMNDData(int mzone_max=-1);
  void DischargedFuelAfterCooling(Burnup &bu, real cooling_year, int num, string *nuc_list, vector<real>& atom_ratio, bool weight=false, bool xs=false);
  void DischargedFuelAfterCooling(int mzones, int mzonee, Burnup &bu, real cooling_year, int num, string *nuc_list, vector<real>& atom_ratio, bool weight=false, bool xs=false);
  void DischargedFuelAfterCooling(Burnup &bu, real cooling_year, int num, vector<string> &nuc_list, vector<real>& atom_ratio, bool weight=false, bool xs=false);
  void DischargedFuelAfterCooling(int mzones, int mzonee, Burnup &bu, real cooling_year, int num, vector<string> &nuc_list, vector<real>& atom_ratio, bool weight=false, bool xs=false);
  void GetPrintingNuclideList(int prt_nuc, string *prt_nuc_nam, vector< vector<int> > &each_nuc_turn);

  real GetKeff(int cyc, int substep);

  void PerturbInitialND(int mzone_inp, int nuclide_id, real rel_change);
  
  // +++ printing
  void PrintMacroXS();
  void PrintMicroXS(int nucid);
  void PrintMacroXSMediumWise(int inp);
  void PrintInitialND();
  void PrintInitialNDForPhits();
  void PrintNuclideWeightPerBatch(Burnup &bu,bool init=false);
  void PrintND(Burnup &bu,string opt="nd");
  void PrintND(int nuc,string *nucname,int medi=0,int mede=-1);
  void PrintBurnup(int medi=0,int mede=-1);
  void PrintFluxLevelHistory();
  void Print1groupXS(string nucname,int medid);
  void PrintLinepowerMap(){print_linepower_map=true;};
  real GetDensity(int med,int batch,int nucid);
  real GetVolumePerMedium(int i){return vol_per_med[i];};
};


#endif


