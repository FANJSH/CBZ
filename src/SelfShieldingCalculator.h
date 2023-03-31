#ifndef SELFSHIELDINGCALCULATOR
#define SELFSHIELDINGCALCULATOR

#include "Numeric.h"
#include "Medium.h"
#include "LibData.h"
#include "OnePointCalculator.h"
#include "PJI_trajectoryset.h"
#include "PJI_slabpij.h"

#include<string>
#include<iostream>
#include<fstream>

using namespace std;

enum pij_calculator{slab,general};

class SelfShieldingCalculator{
  int num_medium,region;
  vector<Medium> medium;
  PJISlabPij *pijcalculator_slab;
  TrajectorySet *pijcalculator_general;
  enum pij_calculator pij_c;
  vector<int> medium_id;
  vector<GroupData1D> dancoff;
  bool background_xs_printing;
  bool cw_correction; // CW-total correction (limited to 107- or 172-group)  
 public:
  SelfShieldingCalculator();
  void BackgroundXSPrinting(){background_xs_printing=true;};  
  void PutPijCalculator(TrajectorySet *tset);
  void PutPijCalculator(PJISlabPij *pspij);
  void SetNumberOfMedium(int i);
  void PutMedium(int number,Medium &medinp);
  void PutMediumID(int *reg_med);
  void GiveInfiniteDilutionCrossSection(XSLibrary &xslib);
  void WithToneMethod(XSLibrary &xslib, bool print=true);
  int GetRegion();
  real GetVolume(int i);
  void CalculationPij(real *xs,real *pij);
  Medium &GetMedium(int i){return medium[i];};
  void CalMixture();
  //
  void ThreeRegionDancoffMethod
    (XSLibrary &xslib, TrajectorySet &tset, Medium &med0, Medium &med1, Medium &med2, bool clad=true,real eng=-1.);
  void FourRegionDancoffMethod
    (XSLibrary &xslib, TrajectorySet &tset, Medium &med0, Medium &med1, Medium &med2, Medium &med3,
     bool clad=true,real eng=-1.);
  void DancoffMethod(XSLibrary &xslib, TrajectorySet &tset, vector<Medium> &med, real eng=-1.);
  GroupData1D &GetDancoff(int mat){return dancoff[mat];};
  //void CalDancoffCorrection(real average_chord);
  //void CalDancoffCorrection(TrajectorySet &tset,int group,int nucid,real average_chord);
  void WithDancoffMethod(XSLibrary &xslib,TrajectorySet &tset);
  //void WriteDancoffCorrection();
  void CWCorrectionOff(){cw_correction=false;};  
};

#endif
