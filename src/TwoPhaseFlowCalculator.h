#ifndef TWO_PHASE_FLOW_CALCULATOR
#define TWO_PHASE_FLOW_CALCULATOR

#include "GroupData.h"
#include "jsmest_wrapper.h"
#include "Numeric.h"


#include<string>
#include<iostream>
#include<fstream>

using namespace std;

class TwoPhaseFlowCalculator{
 protected:
  real length;    // [m]
  real diameter;  // [m]
  real area;      // [m2]
  real dz;        // [m] The length per one mesh 
  int num_node;   // The number of node points (boundaries)

  real c_f; // Friction coefficient
  real g;   // [m/s2] gravity acceleration
  int itermax_rho;  // Maximum number for density iteration
  real eps_rho;     // Criteria for density iteration
  int itermax_vel;

  vector<real> p;    // pressure [Pa] at mesh-center + at exit
  vector<real> rho;  // density [kg/m3] at mesh-center 
  vector<real> h;    // specific enthalpy [J/kg] at mesh-center 
  vector<real> x;    // quality at mesh-center 
  vector<real> vr;   // void rate at mesh-center 
  vector<real> temp; // temperature [K] at mesh-center
  vector<real> q_ex; // Power density [J/(m3s)] at mesh-center
  vector<real> v;    // velocity [m/s] at node point
 public:
  TwoPhaseFlowCalculator(real lin, real din);
  ~TwoPhaseFlowCalculator(){};
  void PutNumNode(int i);
  void PutConstantPowerDensity(real q_in);
  void PutPowerDensity(int i, real q);
  void TransientFixedOutletPressure();
  void Run(real flow_rate_in, real pressure_in, real t_in); // [kg/s] [Pa] [K]
  real RunTransient(real flow_rate_in, real pressure_in, real t_in, real dt, bool conv); // [kg/s] [Pa] [K] [sec]
  void Printing(int digit=4);

  vector<real> GetP(){return p;};
  vector<real> GetRho(){return rho;};
  vector<real> GetH(){return h;};
  vector<real> GetX(){return x;};
  vector<real> GetVr(){return vr;};
  vector<real> GetTemp(){return temp;};
  vector<real> GetV(){return v;};

  vector<real> GetHeatTransferCoefficient();
  
  real get_temp(real p, real h);
  real quality(real p, real h);
  real void_rate(real p, real x, real h, real &density);
  real get_specific_volume(real t, real p);  // Specific volume [m3/kg]
  real get_specific_enthalpy(real t, real p);  // Specific enthalpy [J/kg]
  real get_saturation_temperature(real p);
  real get_dynamic_viscocity(real t, real p); // dynamic viscocity, mu [Pa sec]
  //real get_kinetic_viscocity(real t, real p); // kinetic viscocity, [Pa sec]?    
};

#endif
