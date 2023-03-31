#ifndef BURNUP
#define BURNUP

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "GroupData.h"
#include "GeneralSystem.h"
#include "BurnupChain.h"
#include "MATIDTranslator.h"

using namespace std;

class AtomicMassData{
 private:
  map<int,real> atomic_weight;
  MATIDTranslator midt;
 public:
  AtomicMassData();
  real GetData(int id){return atomic_weight[id];};
  void ReadFile(string cbglibdir,string filename);
  void ShowSelf();
};

class ReactionEnergyData{
 private:
  map<int,real> fis; // Unit: J/reaction
  map<int,real> cap;
  MATIDTranslator midt;
 public:
  ReactionEnergyData();
  void PutGammaEnergySeparatedData();
  real GetFissionEnergy(int mat){return fis[mat];};
  real GetCaptureEnergy(int mat){return cap[mat];};
  void ReadFile(string mdir,string filename);
  // +++ MATIDTranslator
  int ID(string name){return midt.ID(name);};
  string Name(int id){return midt.Name(id);};
};

class Burnup{
 private:
  BurnupChain bc;
  AtomicMassData amd;
  ReactionEnergyData red;

  int nucnum;
  vector<int> nucid;
  vector<real> sgf;
  vector<real> sgc;
  vector<real> sgn2n;
  vector<real> dens;

  int order_krylov;
  GroupData2D trmat_flxdep;
  GroupData2D trmat_flxindep;

  MATIDTranslator midt;
 public:
  Burnup();
  //
  void SetOKUMURAChainFine(bool fr=false){bc.SetOKUMURAChainFine(fr);};
  void SetDefaultChain(){bc.SetDefault();};
  //
  void OverWritingChainData(string mdir,string name){bc.OverWritingChainData(mdir,name);};
  //void OverWritingBranchingRatioDataJENDLDD(){bc.OverWritingBranchingRatioDataJENDLDD();};
  //void OverWritingBranchingRatioDataENDFDD(){bc.OverWritingBranchingRatioDataENDFDD();};
  //void OverWritingBranchingRatioDataJEFFDD(){bc.OverWritingBranchingRatioDataJEFFDD();};
  //
  void PutNucnum(int i);
  int GetNucnum(){return nucnum;};
  void PutNuclideData(int i,int id,real den,real sf,real sc,real sn2n);
  int GetNuclideID(int i){return nucid[i];};
  int SearchNuclide(int id);

  real GetDensity(int i){return dens[i];};
  vector<real> GetDensity(){return dens;};
  void PutDensity(int i,real j){dens[i]=j;};
  void PutDensity(vector<real> &j);
  void SetZeroDensity();
  void ShowDensity();

  real GetAtomicWeight(int mat){return amd.GetData(mat);};
  //real GetAtomicWeight(int mat){return atomic_weight[mat];};

  void ReadReactionEnergyFromFile(string cbglibdir,string filename){red.ReadFile(cbglibdir+"CBGLIB_BURN/ReactionEnergy/",filename);};
  void PutGammaEnergySeparatedData(){red.PutGammaEnergySeparatedData();};

  void CalTransitionMatrixFluxDependentPart();
  void CalTransitionMatrixFluxInDependentPart();
  void UpdateTransitionMatrixFluxInDependentPartForSpontaneousFission(int id, real half_life);

  void AddNuclideToMediumFromBurnupChain(Medium &med);
  GroupData2D CalTransitionMatrix(real flx,bool decay=true);
  GroupData2D CalTransitionMatrix(Medium &med,real flx);
  GroupData2D CalTransitionMatrixByGroup(Medium &med,int grp,real flx);
  GroupData2D CaldTMatdFlux(Medium &med,int grp,int nn=-1);
  GroupData1D CalMatrixExponentialByKrylov(GroupData2D &inp,GroupData1D &inp2,real delt);
  //GroupData2D CalMatrixExponentialByChebyshev(GroupData2D &inp,real delt,int order);
  void BurnupCalculation(real flx, real delt,string opt="pade",bool print=true);

  void BurnupCalculation(Medium &med,real flx,real delt,bool putmed=true);
  void BurnupCalculationByKrylov(Medium &med,real flx,real delt,bool putmed=true);
  void BurnupCalculationByPade(Medium &med,real flx,real delt,bool putmed=true);
  void BurnupCalculationByChebyshev(Medium &med,real flx,real delt,bool putmed=true);
  void BurnupCalculationByMMPA(Medium &med,real flx,real delt,bool putmed=true);
  void BurnupCalculationByMultiStepCalc(Medium &med,GroupData1D &nuc_int,real flx_level,real delt,bool putmed);

  void BurnupCalculationNew(real flx, real delt,bool print=true);
  void BurnupCalculationByKrylov(real flx,real delt);
  //
  void PutOrderKrylov(int i){order_krylov=i;};
  void PutOrder(int i){PutOrderKrylov(i);};
  void PutMediumData(Medium &med);
  void CalOneGroupCrossSection(Medium &med,GroupData1D &flx);
  void PutMediumDataWithZeroCrossSection(Medium &med);
  // Unit of flux /cm2/sec :  Unit of delt : sec

  real CalPowerNormalizationFactor(GeneralSystem &sys,real power);
  // [Unit : W]
  real GetIntegratedPower(Medium &med,GroupData1D &flx,bool capture=true);
  real GetIntegratedPower(Medium &med,bool capture=true);
  real GetIntegratedPower(real total_flux,bool capture=true);
  real GetDecayPower(Medium &med,int type);
  real GetIntegratedPowerParUnitVolume(Medium &med,bool capture=true){return GetIntegratedPower(med,capture);};
  real GetIntegratedPowerParMesh(GeneralSystem &sys,int m,bool cap=true);
  real GetIntegratedPowerParMesh(GeneralSystem &sys,int x,int y,int z=0,bool cap=true);
  real GetIntegratedPowerParMedium(GeneralSystem &sys,int m);
  real CalTotalPower(GeneralSystem &sys,bool cap=true);
  real CalMaximumFastNeutronFlux(GeneralSystem &sys,real e=1e5);
  real CalConversionRatio(GeneralSystem &sys);
  real CalWeightOfHeavyNuclideParUnitVolume(Medium &med, string keyword="all");
  real CalWeightOfAllNuclideParUnitVolume(Medium &med);
  // [Unit : g]
  real CalWeightOfHeavyNuclide(GeneralSystem &sys,int medid1,int medid2);
  real CalWeightOfAllNuclide(GeneralSystem &sys,int medid1,int medid2);
  ReactionEnergyData &GetReactionEnergyData(){return red;};
  // 
  real GetSigf(int i){return sgf[i];};
  real GetSigc(int i){return sgc[i];};
  real GetSign2n(int i){return sgn2n[i];};
  void PutSigf(int i, real j){sgf[i]=j;};
  void PutSigc(int i, real j){sgc[i]=j;};
  void PutSign2n(int i, real j){sgn2n[i]=j;};
  GroupData2D &GetTrmatFlxDep(){return trmat_flxdep;};
  GroupData2D &GetTrmatFlxInDep(){return trmat_flxindep;};
  real GetDecayConstant(int i){return bc.GetDecayConstant(i);};
  //
  BurnupChain &GetBC(){return GetBurnupChain();};
  BurnupChain &GetBurnupChain(){return bc;};
  void CalKermaFactor(Medium &med,bool cap=true);
  // +++ MATIDTranslator
  int ID(string name){return midt.ID(name);};
  string Name(int id){return midt.Name(id);};
  //
  void ReadAtomicMassDataFromFile(string cbglibdir,string filename){amd.ReadFile(cbglibdir, filename);};
  //
  void SetFissionYieldAsInitialNumberDensity(int matid);
  //void SetFissionYieldAsInitialNumberDensity(string fisnuc);
  void PutNuclideDataFromBurnupChain();
  void CoolingCalculation(real init_dt,int step,int sub_step=1);
  void CoolingCalculationByKrylov(real init_dt,int step);
};

#endif


