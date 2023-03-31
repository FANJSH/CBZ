#ifndef MULTICELLBURNER
#define MULTICELLBURNER

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

#include "Numeric.h"
#include "SensitivityData.h"
#include "Medium.h"
#include "IrregularGeometryInformation.h"
#include "PJI_trajectoryset.h"
#include "PJI_system.h"
#include "MEC_system.h"
#include "Burnup.h"
#include "SelfShieldingCalculator.h"
#include "OnePointCalculator.h"
#include "LibData.h"

#include "GeneralBurner.h"

using namespace std;

class MulticellBurner:public GeneralBurner{
 protected:
  TrajectorySet sys_f2;  // for NEL2020
  int nucn;
  real fuel_r;
  int totm;
  int mednum_nonfuel;
  vector<int> region_medium;
  vector<real> vol_med;
  // Number density information at initial state
  vector<real> acday;
  vector<real> acburn;

  vector<int> init_nucnum;
  vector< vector<int> > init_nucid;
  vector< vector<real> > init_nucden;
  vector<real> init_temperature;
  vector<int> nuclide_info; // 0:no cross section, 1:fissile, 2:other

  vector< vector< vector<GroupData1D> > > fwd_nuc;     // [step][substep][medid](nucn)
  vector< vector< vector<GroupData1D> > > fwd_nuc_int; // [step][substep][medid](nucn)

  // (in the case of PC calculations, result of corrector step is stored.) 
  vector< vector< vector<real> > > xsc_1g; // [step][medid][nucn]
  vector< vector< vector<real> > > xsn2n_1g;
  vector< vector< vector<real> > > xsf_1g;
  vector< vector< vector<real> > > total_flux;    // [step][substep][medid] 

  vector< vector< vector<real> > > xsc_1g_p;
  vector< vector< vector<real> > > xsn2n_1g_p;
  vector< vector< vector<real> > > xsf_1g_p;
  vector< vector< vector<real> > > total_flux_p;  // [step][substep][medid]

  vector< vector<real> > delt;
  vector< vector< vector<GroupData1D> > > mic_sigf; //[step][medid][nucn](group)
  vector< vector< vector<GroupData1D> > > mic_sigc;
  vector< vector< vector<GroupData1D> > > mic_sign2n; 
  vector< vector<GroupData1D> > volflx_mesh; // [step][TotM](group)
  vector< vector<real> > power_factor;

  bool pij_storing; 
  vector< vector< vector<GroupData2D> > > pij_store; // [step][predict/corrector][group](TotM*TotM)
  vector< vector<real> > pow_adj;
  vector< vector< vector<GroupData1D> > > adj_nuc;
  vector< vector<GroupDataSet> > macxs;
  vector< vector<GroupData1D> > bilinear_flx; // [step][medid](group)
  vector< vector<GroupData1D> > bilinear_aflx; // [step][medid](group)
  vector< vector<GroupData1D> > chi_gpt_fwdflx; // [step][medid](group) 
  vector< vector< vector<GroupData1D> > > volaflx_mesh; // [step][TotM][group](totsn)


  // GPT with WPC
  vector< vector< vector<GroupData1D> > > fwd_nuc_p;     // [step][substep][medid](nucn)
  vector< vector< vector<GroupData1D> > > fwd_nuc_p_int; // [step][substep][medid](nucn)
  vector< vector<GroupData1D> > volflx_mesh_p; // [step][TotM](group)
  vector< vector< vector<GroupData1D> > > volaflx_mesh_p; // [step][TotM][group](totsn)
  vector< vector<GroupDataSet> > macxs_p;

  vector<real> keff_p;
  vector< vector<real> > power_factor_p;
  vector< vector<real> > pow_adj_c;
  vector< vector< vector<GroupData1D> > > adj_nuc_c;
  vector< vector<GroupData1D> > bilinear_flx_c; // [step][medid](group)
  vector< vector<GroupData1D> > bilinear_aflx_c; // [step][medid](group)
  vector< vector<GroupData1D> > chi_gpt_fwdflx_c; // [step][medid](group) 
  real wc_gpt_wpc;

  bool corrector_calc;
  bool wpc_direct_calc;

  // (Spherical harmonics expansion for angular flux)
  Quadrature quad;
  int aflx_legendre;
  vector< vector< vector<GroupData1D> > > volaflx_pl;   // [step][TotM][group](l,m)
  vector< vector< vector<GroupData1D> > > volaflx_pl_p; // [step][TotM][group](l,m)

  vector<real> hm_weight_init_per_medium;
  vector< vector<real> > acburn_per_medium;
  
  // Dancoff factor
  bool dancoff_input;
  vector<real> dancoff_factor;

  // Group collapsing during burnup
  bool collapsing;
  string collapse_dir;


  void ForwardCalculation(Burnup &bu, int med_target, bool adjoint=false);
  void ForwardCalculationOld(Burnup &bu, int med_target, bool adjoint=false);
  void ForwardCalculation20210108(Burnup &bu, int med_target, bool adjoint=false);
  
  void ForwardCalculationNEL2020(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp);
  void ForwardCalculationADTS_toOWPC(Burnup &bu, int med_target, int ngrp, int *bgrp); // fuga additon
  void ForwardCalculationNEL2020_3x3(Burnup &bu, int med_target,int sm_step, int ngrp, int *bgrp);
  void ForwardCalculationNEL2020_3x3_pc(Burnup &bu, int med_target,int sm_step, int ngrp, int *bgrp); // fuga addition 
  void ForwardCalculationADTS_SamSys(Burnup &bu,int med_target,int sm_step,int ngrp,int *bgrp); // fuga addition
  void ForwardCalculationNEL2020Sasuga(Burnup &bu, int med_target);
  void ForwardCalculationNEL2020SasugaFinal(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp);
  void ForwardCalculationADTS_DifSys(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp); // fuga addition   
  void SensitivityCalculationPC(Burnup &bu, int med_normalize, int num_med, int *med_target, int num, int *nucid,string filename);

 public:
  MulticellBurner();
  ~MulticellBurner(){};
  void PutIGI(IrregularGeometryInformation &igia,IrregularGeometryInformation &igib,bool reflective=false);
  void PutIGIQuarter(IrregularGeometryInformation &igia,IrregularGeometryInformation &igib);
  void PutIGI_NEL2020_3x3(IrregularGeometryInformation &igia,IrregularGeometryInformation &igib,int num_divAng, real dist_trac, bool reflective=false);
  void PutIGI_adts_SameSystem(IrregularGeometryInformation &igia,IrregularGeometryInformation &igib,int num_divAng, real dist_trac, bool reflective=false);  
  void PutIGI_NEL2020(IrregularGeometryInformation &igi_33, int num_divAng, real dist_trac);
  void PutIGI_adts_DifSystem(IrregularGeometryInformation &igi_33, int num_divAng, real dist_trac);  // NEL2020
  void PutMednumFuel(int i);
  void AddNonfuelData(int i,int *id,real *den,real temp);
  void PutFuelDataMB(int i,int *id,real *den,real temp,int medid);
  void PutRelationRegionMedium(int *inp);
  void PutRelationRegionMedium(vector<int> &inp);
  void PutNuclideDataToMedium(GroupData1D &den, int medid_med);
  void PutMicroXSDataToMedium(GroupData1D &sf,GroupData1D &sc,GroupData1D &s2n,int nuc,int medid_med,int nuclide_info);
  //void Run(Burnup &bu,string filename);
  void PreCalculation(Burnup &bu);
  void Calculation(Burnup &bu, int med_target, bool adjoint=false);
  void CalculationADTS_toOWPC(Burnup &bu, int med_target, int ngrp, int *bgrp);
  void CalculationNEL2020(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp);
  void CalculationNEL2020_3x3(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp);  
  real CalculationPinPower(Burnup &bu, int st, int sst, int medid, real vol_flx);
  void CalculationPinBurnup(Burnup &bu, int st, int sst, int medid, vector<real> &xsf1g, vector<real> &xsc1g, vector<real> &xsn2n1g, real tflx, real period, bool mstep=false);
  void PijStoringOff(){pij_storing=false;};
  void PredictorCorrector(){corrector_calc=true;};
  void WPCDirectCalculation(){wpc_direct_calc=true; wc_gpt_wpc=0.6;};
  void PutDancoffFactor(real *dancoff_inp);
  void GroupCollapsingDuringBurnup(string dirname);
  // +++ For OWPC study
  void PreEigenvalueCalculation(Burnup &bu);
  real EigenvalueCalculation(Burnup &bu);
  void PerturbNumberDensity(int nucid, int medid, real factor=1.01);
  real GetMicroscopicReactionRate(int nucid, int medid);
  // +++ Sensitivity Calculation
  void SensitivityCalculation(Burnup &bu, int med_normalize, int med_target, int num, int *nucid, string filename);
  void SensitivityCalculation(Burnup &bu, int med_normalize, int num, int *nucid, string filename);
  void SensitivityCalculation(Burnup &bu, int med_normalize, int num_med, int *med_target, int num, int *nucid, string filename);
  void SensitivityCalculation(Burnup &bu, int med_normalize, int num_med, int *med_target, int num, string *nucname, string filename);
  void SensitivityCalculation(Burnup &bu, int med_normalize, int med_target, int num, string *nucname, string filename);
  void SensitivityCalculation(Burnup &bu, int med_normalize, int num, string *nucname, string filename);
  void IntegratingForwardNumberDensity(vector< vector< vector<GroupData1D> > > &fwd_nuc, vector< vector< vector<GroupData1D> > > &fwd_nuc_int);
  void CalNadjDMNfwd(int inuc, int matnum, int bc_channel, Burnup &bu, vector< vector< vector<GroupData1D> > > &fwd_nuc, vector< vector< vector<GroupData1D> > > &adj_nuc, vector< vector< vector<real> > > &nadj_dm_nfwd);
  //void SensitivityCalculationDirect(Burnup &bu, int medid, int num, int *nucid,int matid, int mt,string filename);
  void SensitivityCalculationDirect(Burnup &bu, int med_normalize, int matid, int mtid, int target_num, int *medid, int *nucid, string filename,bool sigt_preserve=false);
  void SensitivityCalculationKeffEOC(Burnup &bu, int med_normalize, string filename); // Predictor-corrector calculation is expected.
  SensitivityData CalInitialAdjNDKeffEOC(real &response, vector<GroupData1D> &adj_nuc_e);
  // +++ Output printing
  void ShowEigenvalue();
  void ShowCrossSection(int step, int medid);
  void ShowCrossSectionNuclideWise(int nucid, bool reactionrate=false);
  void ShowAbsorptionRate(int medid=-1);
  void ShowNumberDensity(int medid=-1);
  void ShowNumberDensityReactionRate(int nucid, int nucid2=-1);
  void ShowNumberDensityChange(int medid, real limit=0.);
  void ShowNeutronFluxHistory();
  void ShowNumberDensityHistory(int prt_nuc, string *prt_nuc_nam,Burnup &bu,string opt="nd_per_vol",bool sum_only=false);
  void ShowBurnupHistory();
  void ShowFuelNeutronFlux(int st, bool excel=false);
  // +++ (Writing calculation results on external file)
  void WriteFileNumberDensityEOC(string mdir, string filename, int digit=7);
  void WriteFileNumberDensity(int cyc_num, int *cyc, string mdir, string filename, int medid, int digit=7);  
  void ReadFileNumberDensity(int pos, string mdir, string filename, int medid);
  void ReadFileNumberDensityForInitialState(Burnup &bu, string mdir, string filename);
  void ReadFileNumberDensityForInitialState(Burnup &bu, string mdir, string filename, int medid);    
};

#endif


