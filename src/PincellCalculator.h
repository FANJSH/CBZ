#ifndef PINCELLCALCULATOR
#define PINCELLCALCULATOR

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "PJI_trajectoryset.h"
#include "Medium.h"
#include "LibData.h"
#include "IrregularGeometryInformation.h"
#include "OnePointCalculator.h"
#include "GeneralOption.h"
#include "SelfShieldingCalculator.h"
#include "PJI_system.h"
#include "MEC_system.h"
#include "SensitivityData.h"
#include "GData.h"

using namespace std;

class PincellCalculator{
 private:
  int group,mednum;

  real pin_pitch;
  int ring;
  vector<real> radius;
  vector<int> medium_id;
  vector<Medium> med;

  bool hexagonal;
  bool geom_print;

  TrajectorySet sys,sys_f;
  enum BCondition bc_ssc; // (self-shielding calculation)
  enum BCondition bc_flx; // (eigenvalue calculation)

  bool background_xs_printing;

  real thermal_cutoff_energy_coolant;
  real thermal_cutoff_energy_fuel;

  bool cw_correction; // CW-total correction
  int weight_transport_xs_fuel;

 public:
  PincellCalculator(bool hex_in=false);
  void BackgroundXSPrinting(){background_xs_printing=true;};  
  void PutGroup(int i){group=i;};
  void PutMednum(int i);
  void PutPinPitch(real i){pin_pitch=i;};
  void PutRing(int i);
  void PutRadius(real *i);
  void PutMediumID(int *i);
  void PutBoundaryCondition(enum BCondition bci_ssc,enum BCondition bci_flx);
  void PutMedium(vector<Medium> &minp);
  void GeometoryPrinting(){geom_print=true;};
  void SetTrajectorySet();
  void SelfShieldingCalculation(XSLibrary &xslib);
  real EigenvalueCalculation(XSLibrary &xslib, bool transport=true);
  real EigenvalueCalculationByMEC(XSLibrary &xslib);
  void ReactivityCalculation(XSLibrary &xslib, vector<Medium> &minp_pert);
  GData ReactivityCalculation(XSLibrary &xslib, XSLibrary &xslib2,vector<Medium> &minp_pert);  
  void ReactivityCalculationByMEC(XSLibrary &xslib, vector<Medium> &minp_pert);  
  void ReactivityCalculation(XSLibrary &xslib, XSLibrary &xslib2);
  void CurrentWeightTotalCalculation(XSLibrary &xslib,string filename,vector<real> &err);
  Medium GetHomogenizedCrossSection(XSLibrary &xslib,bool ssf_cal=true);
  Medium GetHomogenizedCrossSection(XSLibrary &xslib,int tmesh,int *meshid,bool ssf_cal=true);
  Medium GetHomogenizedCrossSectionBilinear(XSLibrary &xslib,bool ssf_cal=true);
  Medium GetHomogenizedCrossSectionBilinear(XSLibrary &xslib,int tmesh,int *meshid,bool ssf_cal=true);
  Medium &GetMedium(int i){return med[i];};
  SensitivityData CalKinfSensitivity(XSLibrary &xslib,int nucnum,int* nucid,bool ssf_cal=true);
  SensitivityData CalKinfSensitivityByMEC(XSLibrary &xslib,int nucnum,int* nucid,bool ssf_cal=true);
  SensitivityData CalRRRSensitivity
    (int numenum,int* numeid, enum xstype* numetype,
     int denomnum,int* denomid, enum xstype* denomtype,
     int ss_nuc, int* matid, XSLibrary &xslib);
  SensitivityData CalRRRSensitivityByMEC
    (int numenum,int* numeid, enum xstype* numetype,
     int denomnum,int* denomid, enum xstype* denomtype,
     int ss_nuc, int* matid, XSLibrary &xslib);
  SensitivityData CalRRRSensitivityDirect
    (int numenum,int* numeid, enum xstype* numetype,
     int denomnum,int* denomid, enum xstype* denomtype,
     int matid,enum xstype xst,XSLibrary &xslib);
  void FissionSpectrumVectorCalculation(XSLibrary &xslib);
  void ShowNeutronFluxInFuel(XSLibrary &xslib, bool adj=false, bool ssf_cal=true);
  GroupData1D CalNeutronFluxInFuel(XSLibrary &xslib, bool adj=false, bool ssf_cal=true);
  TrajectorySet& GetTrajectorySetFlux(){return sys_f;};
  vector<int> & GetMediumID(){return medium_id;};
  void ReadFile(string mdir,string filename);
  void PutThermalCutoffEnergyFuel(real inp){thermal_cutoff_energy_fuel=inp;};
  void CWCorrectionOff(){cw_correction=false;};
  void FluxWeightedTotalForTransportXS(){weight_transport_xs_fuel=0;};
};


#endif


