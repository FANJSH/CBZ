#ifndef UNCERTAINTYCALCULATION
#define UNCERTAINTYCALCULATION

#include "UNC_CrossSection.h"
#include "UNC_Sensitivity.h"
#include "UNC_Parameters.h"
#include "GroupData.h"
#include "Numeric.h"
#include "MATIDTranslator.h"
#include "BurnupChain.h"
#include "YieldDecayCovariance.h"

const int MAX_DAT=20000;

class UncertaintyCalculation{
 protected:
  LibraryCovariance *cov_con;
  SensitivityContainer *sens_con;
  LibraryContainer *lib_con;
  ParameterCovariance *cval_cov;
  ParameterCovariance *eval_cov;
  ParametersContainer *cval_con;
  Parameters *eval;
  GroupData1D ebnd;
  int grp;
  MATIDTranslator midt;
  // For Uncertainty Quantification
  string filename_uq;
  bool uq_on;
  int count;
 public:
  UncertaintyCalculation
  (int g,real *bnd,LibraryCovariance *cinp,
   SensitivityContainer *sinp,LibraryContainer *linp,
   ParametersContainer *cvalinp, Parameters *evalinp,
   ParameterCovariance *ccovinp, ParameterCovariance *ecovinp);
  UncertaintyCalculation(int g,real *bnd);
  UncertaintyCalculation(XSLibrary &xslib);
  UncertaintyCalculation();

  void PutEnergyGroupInfo(int g,real *bnd);
  GroupData1D GetEbnd(){return ebnd;};

  void GetCEValue(ParametersContainer &cval_con,Parameters &eval,string libname);
  void GetCEValue(ParametersContainer &cval_con,Parameters &eval);

  void ShowCrossSectionStandardDeviation(LibraryCovariance &libcov,int mat,int mt,bool excel=false);
  void ShowCrossSectionStandardDeviationExcel(LibraryCovariance &libcov,int mat,int mt)
  {ShowCrossSectionStandardDeviation(libcov,mat,mt,true);};
  void ShowCorrelationForXYPlot(LibraryCovariance &libcov,int mat,int mt);

  void PrintingOn(string filename){filename_uq=filename; uq_on=true;};
  real CalCrossSectionUncertainty
    (SensitivityContainer &sens_con,LibraryCovariance &lib_cov,
     string core,string chara,int step=0,real mean=0.0001,bool print=true);
  void CalCrossSectionUncertaintyDetail
    (SensitivityContainer &sens_con,LibraryCovariance &lib_cov,
     int mat, int mt, string core,string chara,int step=0);
  void ShowCrossSectionUncertaintyComponentXYPlot
    (SensitivityContainer &sens_con,LibraryCovariance &lib_cov,
     string core,string chara,int step,int mat1,int mt1,int mat2,int mt2,
     string filename="out");
  void ShowCrossSectionUncertaintyComponentXYPlot
    (GroupData1D &s1,GroupData1D &s2, GroupData2D &cov, real factor, string filename="out");
  real CalG1MG2
    (SensitivityContainer &sens_con,LibraryCovariance &lib_cov,
     string core,string chara,int step,string core2,string chara2,int step2);
  real CalG1MG2
    (SensitivityData &sens1,SensitivityData &sens2,LibraryCovariance &lib_cov);
  real CalUncertaintyWithoutAnyInformation(string core,string chara,int step=0);

  real CalCrossSectionUncertaintyBiasMethod
    (string core1,string chara1,string core2,string chara2,
     int step1=0,int step2=0,bool print=true);
  real CalUncertaintyBiasMethod
    (string core1,string chara1,string core2,string chara2,
     int step1=0,int step2=0,bool print=true);
  void SearchBestMockup(string core,string chara,int step=0);

  // +++ Library effect calculation
  real CalLibraryEffect
    (SensitivityContainer &sens_con,LibraryContainer &lib_con,
     string core,string chara,int step,string libname1,string libname2,real mean=0.0001,bool energy=false,bool print=true);
  void CalLibraryEffect
    (SensitivityContainer &sens_con, BurnupChain &chain1, BurnupChain &chain2,
     string core,string chara,int step=0,real mean=0.0001);
  void EstimateCEValueFromSensitivity
    (SensitivityContainer &sens_con,LibraryContainer &lib_con, string libname1,string libname2,
     Parameters &cal_val, Parameters &exp_val, ParameterCovariance &exp_cov);
  void EstimateCEValueFromSensitivity
    (SensitivityContainer &sens_con,LibraryContainer &lib_con, string libname1,string libname2,
     Parameters &cal_val, Parameters &exp_val, ParameterCovariance &exp_cov,string *name);

  ParameterCovariance GetGMG
  (Library &lib,LibraryCovariance &xscov_cov,Parameters &cal_val,SensitivityContainer &sens);
  GroupData2D GetGMG
    (LibraryCovariance &xscov_cov,Parameters &cal_val1,Parameters &cal_val2,SensitivityContainer &sens1,SensitivityContainer &sens2);

  void DoUncertaintyLibraryAdjustment
  (SensitivityContainer &sens_con, ParameterCovariance &cal_cov,
   ParameterCovariance &exp_cov, Parameters &cal_val,
   Parameters &exp_val, Library &lib,LibraryCovariance &lib_cov);
  void DoUncertaintyLibraryAdjustmentWithNewIntegralDataPrediction
  (SensitivityContainer &sens_con, ParameterCovariance &cal_cov,
   ParameterCovariance &exp_cov, Parameters &cal_val,
   Parameters &exp_val, Library &lib,LibraryCovariance &lib_cov,
   SensitivityContainer &sens_con_im, Parameters &cal_val_im);


  // Decay heat
  void CalCEwithUncertainty
  (SensitivityContainer &sens_con, ParameterCovariance &cal_cov, ParameterCovariance &exp_cov,
   Parameters &cal_val, Parameters &exp_val,
   IndependentYieldCovariance &icov, HalfLifeCovariance &hcov,
   DecayEnergyCovariance &dcov, BranchingRatioCovariance &bcov,bool adjustment=false);
  void DoUncertaintyLibraryAdjustmentForBurn
  (SensitivityContainer &sens_con, SensitivityContainer &sens_con_im, 
   ParameterCovariance &cal_cov, ParameterCovariance &exp_cov,
   Parameters &cal_val, Parameters &cal_val_im,
   Parameters &exp_val, Library &lib,
   LibraryCovariance &lib_cov, IndependentYieldCovariance &icov,
   HalfLifeCovariance &hcov, DecayEnergyCovariance &dcov,
   BranchingRatioCovariance &bcov);//kawamoto
  //void CondSensitivity(int cgrp, int *bgrp,string name);


  void PutCvalueFromSensitivity(SensitivityContainer &sens_con,Parameters &cal_val);

  void MIMP(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,string core,string chara,int step,real mean=0.01)
  {CalVarianceReductionFactor(sens_con,lib_cov,core,chara,step,mean);};
  void CalVarianceReductionFactor(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,string core,string chara,int step,real mean);

};

  //(kawamoto)
class UncertaintyCalculationForYieldDecay{ //kawamoto
 protected:
 public:
  UncertaintyCalculationForYieldDecay(){};
  ~UncertaintyCalculationForYieldDecay(){};

  void CalFissionYieldUncertainty(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,string core,string chara,int step);
  void CalFissionYieldUncertaintyDetail(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,string core,string chara,int step,real mean);
  void CalDecayEnergyUncertainty(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step);
  real CalDecayEnergyUncertaintyNuclideWise(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step, string nucname);
  void CalDecayEnergyUncertaintyDetail(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step,real mean);
  void CalHalfLifeUncertainty(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step);
  real CalHalfLifeUncertaintyNuclideWise(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step, string nucname);
  void CalHalfLifeUncertaintyDetail(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step,real mean);
  void CalBranchingRatioUncertainty(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,string core,string chara,int step);
  void CalBranchingRatioUncertaintyDetail(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,string core,string chara,int step,real mean,bool decay=true);
  GroupData2D GetFissionYieldGMG(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,int chara_num,string *core,string *chara,int *step);
  GroupData2D GetFissionYieldGMG(SensitivityContainer &sens_con,IndependentYieldCovariance &cov);
  GroupData2D GetFissionYieldGMG(SensitivityContainer &sens_con1,IndependentYieldCovariance &cov,SensitivityContainer &sens_con2);
  GroupData2D GetDecayEnergyGMG(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,int chara_num,string *core,string *chara,int *step);
  GroupData2D GetDecayEnergyGMG(SensitivityContainer &sens_con,DecayEnergyCovariance &cov);
  GroupData2D GetDecayEnergyGMG(SensitivityContainer &sens_con1,DecayEnergyCovariance &cov,SensitivityContainer &sens_con2);
  GroupData2D GetHalfLifeGMG(SensitivityContainer &sens_con,HalfLifeCovariance &cov,int chara_num,string *core,string *chara,int *step);
  GroupData2D GetHalfLifeGMG(SensitivityContainer &sens_con,HalfLifeCovariance &cov);
  GroupData2D GetHalfLifeGMG(SensitivityContainer &sens_con1,HalfLifeCovariance &cov,SensitivityContainer &sens_con2);
  GroupData2D GetBranchingRatioGMG(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,int chara_num,string *core,string *chara,int *step,bool decay=true);
  GroupData2D GetBranchingRatioGMG(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,bool decay=true);
  GroupData2D GetBranchingRatioGMG(SensitivityContainer &sens_con1,BranchingRatioCovariance &cov,SensitivityContainer &sens_con2,bool decay=true);
  // okumura
  GroupData2D GetFissionYieldCorrGMG(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &fbcov,int chara_num,string *core,string *chara,int *step);
  GroupData2D GetFissionYieldCorrGMG(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &fbcov); // okumura
  GroupData2D GetFissionYieldCorrGMG(SensitivityContainer &sens_con1,FissionYieldCovarianceCorr &fbcov,SensitivityContainer &sens_con2); // okumura
  real CalFissionYieldCorrUncertainty(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &cov,string core,string chara,int step);
  void CalFissionYieldCorrUncertaintyDetail(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &cov,string core,string chara,int step,real mean);
};

#endif
