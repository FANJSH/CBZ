#ifndef HEAT_CONDUCTION_CALCULATOR
#define HEAT_CONDUCTION_CALCULATOR

#include "GroupData.h"
#include "Medium.h"
#include "PLOS_system.h"
#include "CartMeshInfo.h"
#include "GeneralOption.h"

#include<string>
#include<iostream>
#include<fstream>

using namespace std;

class HeatConductionCalculator{
 protected:
  //Medium med_pellet, med_gas, med_cladding;
  PLOSSystem test;
  CartMeshInfo cmi;
  int mesh_p, mesh_g, mesh_c, totm;
  real rad_p, rad_g, rad_c;
  real temp_boundary;
  real v_pellet;
  vector<real> r;
  bool xs_corr;
  real theta;
 public:
  HeatConductionCalculator(int mp, int mg, int mc){Initialize(mp,mg,mc);};
  HeatConductionCalculator(){};  
  ~HeatConductionCalculator(){};
  void Initialize(int mp, int mg, int mc);
  void SetTemperatureBoundary(real tin);
  void Run();
  void RunTransient(real temp_b_in, real line_power); // [K] and [W/cm]
  void StationaryStateCalculation(real temp_b_in, real line_power); // [K] and [W/cm]
  void TemperatureDistributionIteration(vector<real> &src_ex);
  void ResultPrinting();
  vector<real> GetTemperatureData();
  real GetMeshWidth(int i);
  void CrossSectionCorrectionForDynamicCalculation(real theta_inp, real dt);
  real GetHeatFluxAtBoundaryPerUnitHeight(); // [W/cm]
  real GetVolumeAveragedTemperature(); // [K]
  real GetCenterTemperature(); // [K]
  real GetOuterTemperature(); // [K]
};

#endif
