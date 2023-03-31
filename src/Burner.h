#ifndef BURNER
#define BURNER

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

#include "Numeric.h"
#include "LibData.h"
#include "Medium.h"

#include "IrregularGeometryInformation.h"
#include "PJI_trajectoryset.h"
#include "PJI_system.h"
#include "Burnup.h"
#include "BurnupChainGenerator.h"
#include "DelayedNeutronData.h"
#include "SensitivityData.h"

#include "GeneralBurner.h"
#include "GData.h"
#include "RMforFR.h"

using namespace std;

class Burner:public GeneralBurner{
 private:
  int nucn;
  //
  int *region_medium;
  real fuel_r,clad_r,fuel_vol,clad_vol,mod_vol;
  int mesh;
  int mesh_fuel; // The number of meshes in fuel region
  int mesh_clad; // Mesh ID of cladding region
  //
  vector< vector<real> > density_data;
  vector< vector<real> > xsc_1g;
  vector< vector<real> > xsn2n_1g;
  vector< vector<real> > xsf_1g;
  //vector< vector<real> > xsnsf_1g;
  vector<real> acday;
  vector<real> acburn;
  vector<GroupData1D> fuel_flux;
  vector<GroupData1D> clad_flux;
  vector< vector<real> > power_factor;
  vector< vector<real> > delt;
  vector< vector<real> > total_flux; // normalized neutron flux in fuel region
  vector< vector<real> > total_flux_unnormalized; // normalized neutron flux in fuel region
  vector< vector<real> > total_flux_clad; // normalized nxeutron flux in clad region
  //
  vector< vector<real> > four_factor; 
  //
  vector<int> nuclide_info; // 0:no cross section, 1:fissile, 2:other
  //
  real hm_weight_init_u, hm_weight_init_tru;
  real h_hm_ratio;
  //
  GroupData1D dancoff;

  // Group collapsing during burnup
  bool collapsing;
  string collapse_dir, collapse_filename;
  bool collapsing_micro;
  int ngrp;
  vector<int> bgrp;

  //
  bool limit_burnup;
  real k_limit;
  int batch;
  real bumax;
  //
  // (for burnup matrix extraction)
  bool matrix_extract;
  string mat_ext_filename;
  int mat_ext_prc;
  // (for sensitivity calculation)
  bool adj_nuc_data;
  real target_val;
  vector< vector<GroupData1D> > fwd_nuc;
  vector< vector<GroupData1D> > fwd_nuc_int;
  vector< vector<GroupData1D> > adj_nuc;
  vector< vector<GroupData1D> > adj_nuc_e;
  vector< vector<GroupData1D> > adj_nuc_s;
  vector< vector<GroupData2D> > pij_store; // [step][group](TotM*TotM)
  vector< vector<real> > pow_adj;
  vector< vector<GroupData1D> > volflx_mesh; // volume-integrated flux per mesh
  vector< GroupDataSet > macxs;
  vector< vector<GroupData1D> > mic_sigf;
  vector< vector<GroupData1D> > mic_sigc;
  vector< vector<GroupData1D> > mic_sign2n;
  //vector< vector<GroupData2D> > mic_sigs;
  vector< vector<GroupData1D> > gpt_flx;
  vector< vector< vector<GroupData1D> > > dmdf_nuc;  // (dM/dphi) * (NUC_FWD)
  GroupData2D trmat_flxindep;
  // (for cooling calculation)
  bool cooling_cal;
 public:
  Burner();
  ~Burner(){};
  //~Burner(){delete [] region_medium;};

  void CalHeavyMetalInitialWeight(Burnup &bu);

  void CalHHMRatio();

  void AddActinideDecayHeatDataToBurnupChain(string cbglibdir,string fname,Burnup &bu);

  void PreCalculation(bool burn_time_cal=true);

  void ConstructingFuelNuclide(Burnup &bu);
  void PutNucn(int i);
  void SelfShieldingCalculation(int step,GroupData1D &bell_default);
  void Calculation(Burnup &bu,bool adjoint=false);
  void CoolingCalculation(Burnup &bu,int sub_step=1);
  void NuclideWiseReactivityEffectCalculationDuringCooling(Burnup &bu,int nuc,string *nuc_name);
  void CalculationHomogeneous(Burnup &bu);
  void dkdN_Calculation(Burnup &bu,string filename="out.dkdn");
  void BurnupDependentKeffSensitivityCalculation(Burnup &bu,string filename="out.dkdn"){dkdN_Calculation(bu,filename);};
  real CalMacroscopicReactionRate(int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh);
  void GroupCollapsingDuringBurnup(string dirname, string filename, bool micro=false);
  void SetDefaultCollapsingInfo();

  // (burnup-sensitivity calculation)
  void SensitivityCalculation(Burnup &bu,int nucnum,string *nucname,bool fileout=false,string prefilename="");
  void SensitivityCalculationDirect(Burnup &bu, int nucnum,string *nucname,int matid,int mt,string snsname);
  void SensitivityCalculationDecayHeat(Burnup &bu,bool beta=true,bool gamma=true,bool alpha=true,string snsname="");
  void SensitivityCalculationKeffEOC(Burnup &bu,bool vid=false,string prefilename="");
  void SensitivityCalculationKeffBOC(Burnup &bu,string prefilename="");
  void SensitivityCalculationRRRBOC(Burnup &bu,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh);
  void SensitivityCalculationRRREOC(Burnup &bu,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh);
  void GetHMIndex(int *nucid,int &nucnum);
  void SetArrayForSensitivityCalculation();
  void AveragingForwardNumberDensity();
  void IntegratingForwardNumberDensity();
  void AdjointBurnupCalculation(Burnup &bu);
  void SensitivityPrinting(Burnup &bu,real response,string sys_name,string para_name,string lib_name,string filename);
  void SensitivityPrintingForDecayHeat(Burnup &bu,real response,vector<int> &nuc,vector<vector<real> > &energy_sens, string snsname);//kawamoto
  void SensitivityPrintingForDecayEnergy(vector<int> &nuc,vector<vector<real> > &energy_sens, real end_decayheat);//kawamoto
  void WriteFileNumberDensity(string filename);
  void WriteFileAdjointNumberDensity(string filename,int targetID);
  void WriteFileContributionFunction(Burnup &bu,string filename,int targetID,real end_nuc);

  void SensitivityRead(string filename,string nucname);
  void BeffCalculation(DelayedNeutronData &dnd);
  real EigenvalueCalculation(Burnup &bu,bool vid=false);
  GroupData1D GetFuelNeutronFluxInitialState(Burnup &bu, real ebnd=-1.);
  Medium GetHomogenizedCollapsedMedium(Burnup &bu, bool micro=false);
  Medium GetHomogenizedCollapsedMediumWOSSC(Burnup &bu, bool micro=false);  
  Medium GetHomogenizedCollapsedMedium(Burnup &bu, Medium &med_mod, bool micro=false);    
  Medium GetCollapsedFuelMedium(Burnup &bu, bool micro=false);  
  GroupData1D GetFuelFlux(int i){return fuel_flux[i];};
  void Cal1groupBranchingRatio(Burnup &bu, string cbglibdir,string filename,string outfile);

  real GetFuelVol(){return fuel_vol;};
  real GetCladVol(){return clad_vol;};

  void PutGeometryData(real pitch,int ring,real *rri,int *rmedi);
  void PutGeometryDataCylinder(real pitch,int ring,real *rri,int *rmedi);
  void PutThreeRegionCircularGeometry(real *rri);
  void PutThreeRegionHexagonalGeometry(real pitch,real *rri);
  //void PutWhiteBoundary();


  GroupData2D GetTransitionMatrix(Burnup &bu);  

  void WriteFileFuelNumberDensity(string mdir,string ss){med[0].WriteFileNumberDensity(mdir,ss);};
  void WriteFileFuelNumberDensity(string mdir,string ss,int cyc);
  void ReadFileFuelNumberDensity(string mdir,string ss){med[0].ReadFileNumberDensity(mdir,ss);};
  void CheckNumberDensityToRef(string mdir,string ss,int nucnum,string *nucname);
  // +++ (Functions for users' editting)
  void GetNuclideDensity(vector<string> &nuc,vector<real> &den,int step);
  real GetNuclideDensity(string nuc, int step){return GetNuclideDensity(midt.ID(nuc),step);};
  real GetNuclideDensity(int id, int step);
  // +++ (Functions for back-end quantity evaluation)
  void CobaltActivationCalculation(real den_init,real flux_factor=1.);// initial number density of Cobalt-59
  void Mn54ActivationCalculation(real den_init,real flux_factor=1.);// initial number density of Fe-54
  void Reprocessing(real tc, real ru_rh, real fp, real ac);
  void NumberDensityReset(int mat,real den);
  void NumberDensityReset(int mat_st,int mat_ed,real den);
  // +++ (Access function)
  int GetNucn(){return nucn;};
  real GetBurnup(int i){return acburn[i];};
  real GetEigenvalue(int i){return keff[i];};
  real GetFuelFlux(int i,int j){return total_flux[i][j];};
  real Get0Flux(int i){return GetFuelFlux(i,0);};
  int GetBurnStep(){return burn_step;};
  real GetXSC1g(int i, int j){return xsc_1g[i][j];};
  // +++ (For limited burnup calculation)
  void LimitBurnupOn(real kl=1.03, int batchinp=3);
  real GetEOCBurnup();
  real GetAtomWiseWeightAtBOC(Burnup &bu,int atm);
  real GetAtomWiseWeightAtEOC(Burnup &bu,int atm);
  // +++ (For lumped-FP creation)
  void CreatePseudoFP(Burnup &bu,string fissle_name,int nucnum,string *nucname,int outid,string outname);
  void LFP(Burnup &bu,string fissle_name,int nucnum,string *nucname);
  void ShowPseudoFPCrossSection(Burnup &bu,int nucnum,string *nucname);
  // +++ (For burnup matrix extraction)
  void SetMatrixExtraction(string filename="matdat",int prc=16){matrix_extract=true; mat_ext_filename=filename; mat_ext_prc=prc;};
  // +++ (Printing calculation result)
  void ShowEigenvalue();
  void ShowFourFactor();
  void ShowNeutronFluxHistory();
  void ShowNuclideList();
  //void ShowDecayHeat(Burnup &bu,bool detail_print=true,real factor=1.); 
  //void ShowRadioactivity(Burnup &bu,bool detail_print=true,real factor=1.); 
  //void ShowDecayHeatAtomWise(Burnup &bu,string atom_name,real factor=1.);
  void ShowNumberDensityChange(real limit=0.);
  //void ShowNumberDensityChange(Burnup &bu, real limit=0.);
  void ShowNumberDensityChange(int nucnum, string *nuclist);
  void ShowWeightChange(Burnup &bu, real limit=0.);
  void ShowNumberDensityHistoryAtomWise(string atom_name,Burnup &bu,string opt="nd_per_vol",bool sum_only=false);
  void ShowNumberDensityHistoryOld(int prt_nuc, string *prt_nuc_nam,Burnup &bu,string opt="nd_per_vol",bool sum_only=false);
  GData ShowNumberDensityHistory(int prt_nuc, string *prt_nuc_nam,Burnup &bu,string opt="nd_per_vol",bool sum_only=false);
  void ShowPlutonium(Burnup &bu);
  void ShowMA(Burnup &bu);
  void ShowCrossSection(int prt_nuc,string *prt_nuc_nam,bool reaction_rate=false, int digit=4);
  void ShowCrossSection(int step);
  void ShowGroupDependentCrossSection(string nucname,bool excel=false);
  void ShowGroupDependentCrossSectionExcel(string nucname){ShowGroupDependentCrossSection(nucname,true);};
  void ShowCaptureToDecayRatio(int prt_nuc,string *prt_nuc_nam,Burnup &bu);
  void ShowAdjointNumberDensity(int prt_nuc,string *prt_nuc_nam);
  void ShowFuelNeutronFlux(int step,bool excel=false);  
  void ShowFuelNeutronFluxExcel(int step){ShowFuelNeutronFlux(step,true);};
  void ShowHeavyMetalWeightRatio(Burnup &bu);
  void ShowHeavyMetalNumberDensityRatio();
  void GetHeavyMetalWeightRatio(Burnup &bu, vector<real> &output);
  void ShowNuclideWiseContributionForFission(int prt_nuc,string *prt_nuc_nam);
  void ShowRadioactivityRatio(Burnup &bu,string nuc1,string nuc2);
  real GetRadioactivityRatio(Burnup &bu,string nuc1,string nuc2);
  real GetRadioactivity(Burnup &bu,string nuc1);
  void GetPrintNuclide(int prt_nuc,string *prt_nuc_nam,vector<int> &prt_nuc_turn);
  // +++ (Writing calculation results on external file)
  void WriteFileNumberDensity(int cyc_num, int *cyc, string mdir, string filename, int digit=7);
  // +++ (Animation)
  void ShowRadioactivityHistoryForAnimation(Burnup &bu);
  void ShowDataForAnimation(Burnup &bu,string filename,string type="nd");
  void ShowDataForAnimationNew(Burnup &bu,string filename,string type="nd");
  // +++ (For NEL-2019 funding research)
  void UseLibraryDataForNonExistingNuclide(Burnup &bu, string cbglibdir, string nucname);
  void ShowDominantCaptureReactionRate(real threshold=1e-2);
  // +++ (estimation by RMforFR)
#if 1  
  void EstimateTRUcompositionByROM(Burnup &bu, real burnup, real coolyear, bool pwr=true, string filename="");
  real EstimateUenrichment(Burnup &bu);
  vector<real> EstimateTRUcompositionByROMForUO2(Burnup &bu, real burnup, real coolyear, bool pwr=true);
  vector<real> EstimateTRUcompositionByROMForMOX(Burnup &bu, real burnup, real coolyear, real pu_enrichment, bool pwr=true);
#endif  
  // The following was replaced by the above in 2022/5/21.
  /*  
  void EstimateTRUcompositionByROM(Burnup &bu, real burnup, bool pwr=true, string filename="");
  vector<real> EstimateTRUcompositionByROMForUO2(Burnup &bu, real burnup, bool pwr=true);
  vector<real> EstimateTRUcompositionByROMForMOX(Burnup &bu, real burnup, real pu_enrichment, bool pwr=true);  
  */
};


#endif


