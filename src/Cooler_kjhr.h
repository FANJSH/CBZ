#ifndef COOLER
#define COOLER

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

#include "Cooler.h"
#include "Burnup.h"
#include "GroupData.h"
#include "ENDFIDTranslator.h"
#include "MATIDTranslator.h"
#include "BurnupChainGenerator.h"
#include "GammaTool.h"

using namespace std;

class Cooler{
 private:
  Burnup bu;
  ENDFIDTranslator eidt;
  MATIDTranslator midt;
  int nucnum;
  int total_step;
  bool beta;
  bool alpha;
  bool gamma;
  vector<GroupData1D> density;
  vector<GroupData1D> adj_density;
  vector<GroupData1D> int_density;
  vector<GroupData1D> int_adj_density;
  vector<GroupData1D> ave_density;
  vector<GroupData1D> ave_adj_density;
  vector<real> time;
  //
  bool spontaneous_fission;
  int id_sf;
  real hl_sf;
 public:
  Cooler(Burnup &bui);
  ~Cooler(){};
  void PutTotalStep(int i);
  void PutNucnum(int i){nucnum=i;};
  void CoolingCalculation(real init_dt,int step,int sub_step=1);
  void SpontaneousFissionConsideration(int matid, real hl);
  void IrradiationAndCoolingCalculationArbitralTimeStep(int mat,real duration,int step_num,real *time_step);
  void IrradiationAndCoolingCalculation(int mat,real duration, real init_dt, int step, int sub_step=1);
  void CoolingCalculationArbitralTimeStep(int step_num,real *time_step);
  void Adjoint(int step_num,real *time_step,real *power);
  void FiniteIrradiationCalculation(int mat,real dt); // flux level is unity (1.0)
  void FiniteIrradiationCalculation(int mat,int step_num,real *time_step,real *power);
  void ABGSwitch(bool alpha, bool beta, bool gamma);
  void IntegrationForward();
  void IntegrationAdjoint();
  void SensitivityCalculation(int mat,int step_num,real *time_step,real *power);

  void ShowDecayHeat(real tcor=0.);
  void ShowRadioactivity(int prt_nuc, string *prt_nuc_nam);
  void ShowDecayHeatJoule();
  void ShowDecayHeatJoule(int num,string *name);
  void ShowRadioactivityRatio(string nuc1,string nuc2);
  void ShowNumberDensityForXYPlot(int step=-1);
  void ShowDecayEnergyForXYPlot(int step=-1);
  void ShowDecayGammaSpectrum(DecayGammaSpectrum &dgs);

  int GetNucpos(int id);

  void Calculation(int mat,int step_num,real *time_step,real *power);  

};


#endif


