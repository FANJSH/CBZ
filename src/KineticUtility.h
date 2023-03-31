#ifndef KINETICUTILITY
#define KINETICUTILITY

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "BurnupChainGenerator.h"
#include "GroupData.h"
#include "Burnup.h"
#include "BurnupChain.h"
#include "SensitivityData.h"

using namespace std;

class OnePointSolverExplicitFPModel{
 private:
  bool chain_effect;
  vector<string> name;
  vector<int> matid;
  vector<int> atomn;
  vector<int> massn;
  vector<int> level;
  vector<real> lmd;
  vector<real> dnt; // number of emitted delayed neutrons per one decay
  vector<real> dnt_dtr; // number of total emitted DN including daughter nuclides
  vector<real> yield;
  vector<real> delta_yield;
  vector<int> pos;  // BCGmanager index to EFP index
  vector<int> pos2; // EFP index to BCGManager index
  real rho; 
  real beta; 
  real ge;
  real nu;
  real k_p;
  real nu_p, nu_d;
  real l;
  // (for sensitivity calculation)
  int step;
  real delta_t;

  vector<GroupData1D> n_fwd;
  vector<GroupData1D> n_adj;
  vector<real> second; // elapsed time
  vector<real> k_p_hist;

  vector<real> sens1; // for fission yield sensitivity
  vector<real> sens2; // for decay constant sensitivity
  vector< vector<real> > sens3; // for branching ratio sensitivity
 public:
  OnePointSolverExplicitFPModel(bool chain_effect_inp);
  ~OnePointSolverExplicitFPModel(){};
  void PutReactivity(real i){rho=i;};
  void PutBeta(real i){beta=i;};
  void PutGenerationTime(real i){ge=i;};
  void PutNu(real i){nu=i;};
  void PutNu_p(real i){nu_p=i;};
  void PutK_p(real i){k_p=i;};
  void PutKandBeta(real k, real b){PutK_p(k*(1.-b));};
  void PutRhoandBeta(real r, real b){PutKandBeta(1./(1.-r), b);};
  void PutLifetime(real i){l=i;};
  real GetYield(int i){return yield[i];};
  real GetLambda(int i){return lmd[i];};
  int GetMatID(int i){return matid[i];};
  real GetDntdtr(int i){return dnt_dtr[i];};
  GroupData1D GetNfwd(int i);

  void DataInputFromChainData(BCGManager& man, string fisnuc);
  void DataInputFromChainDataAll(BCGManager& man, string fisnuc);
  void ChainModification(BCGManager& man);
  void YieldDataSetting(BCGManager& man, string fisnuc);
  GroupData2D CalDecayMatrix(BCGManager &man);
  GroupData2D CalDecayMatrixKp(BCGManager &man);
  int GetTreatedNuclideNumber(){return name.size();};

  void SetDefaultInitialCondition();
  void SetStationaryStateAsInitialCondition(BCGManager &man, real k_p);    
  void ShowNumberOfNeutrons();
  void ShowDecayHeat(BCGManager &man);
  void CalTransientNew(BCGManager &man,int step, real delta_t);
  void Restarting();

  
  void CalTransient(BCGManager &man,int step, real delta_t);  
  void CalTransientFromBurstFission(BCGManager &man,int step, real delta_t);
  void CalTransientNeutronEmission(string srcnuc,BCGManager &man,int step, real delta_t);
  void CalTransientWithPlotXY(BCGManager &man,int step, real delta_t);
  void YieldFactorize(real fact);
  void CalTransientWithAdjoint(BCGManager &man,int step,real delta_t);
  void CalTransientWithAdjointForDNEmission(BCGManager &man,int step,real delta_t);
  SensitivityData SensitivityCalculation(BCGManager &man, int step1,real nu_d,string tagname,int fisid);
  SensitivityData SensitivityCalculationForDNEmission(BCGManager &man,string tagname,int fisid);
  real CalDecayConstantSensitivityPart(BCGManager &man,GroupData1D &adj,int ii);
  real CalDecayConstantSensitivityPart2(BCGManager &man,GroupData1D &adj,int ii,int id,int jj);
  real CalDecayConstantSensitivityPart3(BCGManager &man,GroupData1D &adj,int ii);
  real CalDecayConstantSensitivityPart4(BCGManager &man,GroupData1D &adj,int ii,int id,int jj);
  void BranchingRatioCheck(BCGManager &man);
  void DensityDataPrinting(int nucp,string* nuc,string type="fwd");
  void DensityDataPrintingForXYPlot();
  void UncertaintyCalculation(BCGManager &man);

  void MakeStationaryState(BCGManager &man);
};

#endif


