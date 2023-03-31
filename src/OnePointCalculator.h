#ifndef ONEPOINTCALCULATOR
#define ONEPOINTCALCULATOR

#include "Numeric.h"
#include "Medium.h"
#include "LibData.h"

#include<string>
#include<iostream>
#include<fstream>

using namespace std;

class OnePointCalculator{
 private:
  bool background_xs_printing;
  bool cw_correction_org; // CW-total correction (limited to 107- or 172-group)
 public:
  //OnePointCalculator(){background_xs_printing=false;};
  OnePointCalculator();
  void BackgroundXSPrinting(){background_xs_printing=true;};
  void B1SpectrumCalculationWithFixedSource(Medium &inp,GroupData1D &src,real square_buckling=1e-10);
  void B1SpectrumCalculation(Medium &inp,real square_buckling=1e-10);
  real BucklingSearchB1Method(Medium &inp);
  void GiveInfiniteDillutionCrossSection(Medium &inp,XSLibrary &xslib);
  void GiveInfiniteDilutionCrossSection(Medium &inp,XSLibrary &xslib){GiveInfiniteDillutionCrossSection(inp,xslib);};
  void CalSelfShieldingWithDancoffCorrection
  (Medium &inp,XSLibrary &xslib,real average_chord,GroupData1D bell_factor,GroupData1D dancoff_correction);
  void CalSelfShieldingWithDancoffCorrection
    (Medium &inp,XSLibrary &xslib,real average_chord,GroupData1D bell_factor,vector<real> dancoff_correction);
  // Energy-independent nuclide-wise Dancoff
  void CalSelfShieldingWithDancoffCorrection
    (Medium &inp,XSLibrary &xslib,real average_chord,GroupData1D bell_factor,vector<real> dancoff_correction,int g);
  // Energy-independent nuclide-wise Dancoff
  void CalSelfShieldingWithDancoffCorrection
  (Medium &inp,XSLibrary &xslib,real average_chord,real bell_factor,GroupData1D dancoff_correction);
  void CalSelfShieldingInfiniteSystemWithSig0Correction(Medium &inp,XSLibrary &xslib,real cor);
  void CalSelfShieldingInfiniteSystem(Medium &inp,XSLibrary &xslib);
  void CalSelfShieldingWithHeterogeneousCorrection
  (Medium &inp,XSLibrary &xslib,real average_chord,real bell_factor=1.20);

  void CoarseGroupCorrectionForCurrentWeightCrossSection(Medium &inp);

  void CalIncidentEnergyDependentChi(Medium &inp,XSLibrary &xslib);
  void CalIncidentEnergyDependentChi(Medium &inp,XSLibrary &xslib,GroupData1D &flux);
  void CalFissionSpectrumMatrix(Medium &inp,XSLibrary &xslib);
  void CalEquivalentFissionSpectrumVector(Medium &inp,GroupData1D &flux);
  bool IsFissionMaterial(Medium &inp,XSLibrary &xslib);

  void CalThermalScatteringMatrix(Medium &inp,XSLibrary &xslib,real cutoff);
  void CheckEnergyGroupConsistency(Medium &inp,XSLibrary &xslib);
  void CWCorrectionOff(){cw_correction_org=false;};
  bool CWCorrection(){return cw_correction_org;};
};

#endif
