#ifndef FRDESIGNTOOL
#define FRDESIGNTOOL

#include "Numeric.h"
#include "Medium.h"
#include "MATIDTranslator.h"
#include "Burnup.h"

class FRDTComposition{
 protected:
  int nucnum;
  vector<int> nucid;
  vector<real> density;

 public:
  FRDTComposition();
  int GetNucnum(){return nucnum;};
  int GetNucid(int i){return nucid[i];};
  real GetDensity(int i){return density[i];};
  real GetDensityFromID(int matid);
  void PutNucnum(int i);
  void ShowSelf(real factor=1.);
  void Multiply(real factor);
  void ShowSelf(MATIDTranslator &midt);
  void AddData(int id,real den);
  void AddData(int inpnum, int* id,real* den);
  void PutNumberDensity(Medium &med);
  void PutOxide(real density_oxide, real num_comp, real num_oxigen, int num_iso_comp, int *nucid_comp, vector<real> &wgt_ratio_comp, AtomicMassData &amd);
  void PutOxide(real density_oxide, real num_comp, real num_oxigen, int num_iso_comp, int *nucid_comp, real *wgt_ratio_comp, AtomicMassData &amd);
  void PutX2O3(real density_oxide, int num_iso_comp, int *nucid_comp, real *wgt_ratio_comp, AtomicMassData &amd)
   {PutOxide(density_oxide, 2., 3., num_iso_comp, nucid_comp, wgt_ratio_comp, amd);};
  void PutX2O3(real density_oxide, int num_iso_comp, int *nucid_comp, vector<real> &wgt_ratio_comp, AtomicMassData &amd)
   {PutOxide(density_oxide, 2., 3., num_iso_comp, nucid_comp, wgt_ratio_comp, amd);};
  void PutXO2(real density_oxide, int num_iso_comp, int *nucid_comp, real *wgt_ratio_comp, AtomicMassData &amd)
   {PutOxide(density_oxide, 1., 2., num_iso_comp, nucid_comp, wgt_ratio_comp, amd);};
  void PutXO2(real density_oxide, int num_iso_comp, int *nucid_comp, vector<real> &wgt_ratio_comp, AtomicMassData &amd)
   {PutOxide(density_oxide, 1., 2., num_iso_comp, nucid_comp, wgt_ratio_comp, amd);};

  FRDTComposition operator*(real a);
  FRDTComposition operator+(const FRDTComposition &second);
};

class FRDTFuelComposition:public FRDTComposition{
 protected:
  int nucnum_hm,nucnum_u;
  real pu_fissile_enrichment; // [wt] (P9&P1)
  real pu_enrichment; // [wt]
  real u5_enrichment; // [wt]
  real o_m; // (O/M)
  real pellet_theoretical_density;  // [Rel.]
  //
  vector<int> uo2_comp_id;
  vector<real> uo2_comp_wgt; // weight ratio of U isotope
  vector<int> puo2_comp_id;
  vector<real> puo2_comp_wgt;
  AtomicMassData amd;
 public:
  FRDTFuelComposition(bool special=false);
  void PutPuFissileEnrichment(real a){pu_fissile_enrichment=a*0.01;}; // wt%
  void PutPuEnrichment(real a){pu_enrichment=a*0.01;}; // wt%
  void PutU5Enrichment(real a){u5_enrichment=a*0.01;}; // wt%
  void PutOM(real a){o_m=a;};
  void PutPelletTheoreticalDensity(real a){pellet_theoretical_density=a*0.01;}; // %->Rel.
  void PutTRUComposition(real rat_pu, real rat_ma, int nuc_pu, string *id_pu, real *inp_pu, int nuc_ma, string *id_ma, real *inp_ma, MATIDTranslator &midt);
  void PutTRUComposition(real rat_pu, real rat_ma, int nuc_pu, int *id_pu, real *inp_pu, int nuc_ma, int *id_ma, real *inp_ma);
  void PutTRUComposition(int nuc,string *id,real *inp,MATIDTranslator &midt); // wt%
  void PutTRUComposition(int nuc,int *id,real *inp); // wt%    
  real GetPelletTheoreticalDensity(){return pellet_theoretical_density;};
  void ShowPuEnrichment();
  void ShowTRUWeightInfo();
  //
  void PutUO2Composition(int nuc,int *id,real *inp);
  void PutPuO2Composition(int nuc,int *id,real *inp);
  void CalNumberDensity(real uo2_wgt,real puo2_wgt,real tho2_wgt);
};

class FRDTMetalFuelComposition:public FRDTComposition{
 protected:
 public:
  FRDTMetalFuelComposition(real smear_density, real pu_enrichment, real zr_composition, int num_pu_iso, string *pu_name, real *pu_vector, AtomicMassData &amd);
  ~FRDTMetalFuelComposition(){};
};

class FRDTSUSComposition:public FRDTComposition{
 protected:
  AtomicMassData amd;
 public:
  FRDTSUSComposition();
  void PutDensityAndRatio(real dinp,real *inp);
};

class FRDTB4CComposition:public FRDTComposition{
 protected:
 public:
  FRDTB4CComposition(real b10c_wt,real theo_den); 
  // b10c: B-10 concentration (0-1.)
  // theo_den:theoretical density (0-1.)
};

class FRDTADSFuelComposition:public FRDTComposition{
 protected:
  real zrn_wgt_density;
  real puman_wgt_density;
  real smear_density;
  int pu_nuc;
  vector<int> pu_composition;
  vector<real> pu_isotopic_ratio;
  int ma_nuc;
  vector<int> ma_composition;
  vector<real> ma_isotopic_ratio;
  int num_non_fuel_data;
  AtomicMassData amd;
 public:
  FRDTADSFuelComposition();
  void PutPuData(int n, int *id, real *rat);
  void PutMAData(int n, int *id, real *rat);
  void PutNonfuelData(int n,int *id, real *den);
  void PutFuelWeightRatio(real zrn,real pun,real man,real pellet_fraction);
};

// ++++++++++++++++++++++++++++++++++++++++

class FRDTBaseGeometry{
 protected:
  // Unit of size [mm]
  real assembly_pitch;
  real pin_outer_diameter;
  real pin_thickness;
  int number_of_pin;
  real spacerwire_diameter;
  real spacerwire_pitch;
  real pellet_outer_diameter;
  real pellet_invoid_diameter;
  real pin_pitch;
  int pin_layer;
 public:
  FRDTBaseGeometry();
  void PutAssemblyPitch(real i){assembly_pitch=i;};
  void PutPinOuterDiameter(real i){pin_outer_diameter=i;};
  void PutPinThickness(real i){pin_thickness=i;};
  void PutNumberOfPin(int i){number_of_pin=i;};
  void PutSpacerwireDiameter(real i){spacerwire_diameter=i;};
  void PutSpacerwirePitch(real i){spacerwire_pitch=i;};
  void PutPelletOuterDiameter(real i){pellet_outer_diameter=i;};
  void PutPelletInvoidDiameter(real i){pellet_invoid_diameter=i;};
  void PutPinPitch(real i){pin_pitch=i;};
  void PutPinLayer(int i){pin_layer=i;};

  int GetNumberOfPin(){return number_of_pin;};

  real GetAssemblyVolume();
  real GetCladVolume();
  real GetPinVolume(){return GetCladVolume();};
  real GetSpacerwireVolume();
  real GetPelletVolume();
  real GetPelletGapVolume();
  real GetWholePinVolume();
  real GetPelletVoidVolume();
  int GetPinLayer(){return pin_layer;};
  real GetPinPitch(){return pin_pitch;};
  real GetAssemblyPitch(){return assembly_pitch;};
  real GetPinOuterDiameter(){return pin_outer_diameter;};
  real GetPinThickness(){return pin_thickness;};
};

class FRDTSubAssemblyGeometry:public FRDTBaseGeometry{
 protected:
  // Unit of size [mm]
  real duct_outersize;
  real duct_thickness;
 public:
  FRDTSubAssemblyGeometry();
  void PutDuctOutersize(real i){duct_outersize=i;};
  void PutDuctThickness(real i){duct_thickness=i;};

  real GetDuctVolume();
  real GetCoolantVolume();
  real GetOutDuctVolume();
  void ShowVolumeInformation();
  void CheckValidity();
  real GetDuctOuterSize(){return duct_outersize;};
  real GetDuctThickness(){return duct_thickness;};

  void MonjuFuelSubAssembly();
  void MonjuRadialBlanketSubAssembly();
  void JSFR1500FuelSubAssembly();
  void JSFR1500RadialBlanketSubAssembly();
};

class FRDTControlRodGeometry:public FRDTBaseGeometry{
 protected:
  // Unit of size [mm]
  real guidetube_outer_diameter;
  real guidetube_thickness;
  real shieldtube_outer_diameter;
  real shieldtube_thickness;
 public:
  FRDTControlRodGeometry(){};

  void PutGuidetubeOuterDiameter(real i){guidetube_outer_diameter=i;};
  void PutGuidetubeThickness(real i){guidetube_thickness=i;};
  void PutShieldtubeOuterDiameter(real i){shieldtube_outer_diameter=i;};
  void PutShieldtubeThickness(real i){shieldtube_thickness=i;};

  real GetGuidetubeVolume();
  real GetShieldtubeVolume();
  real GetCoolantVolume();
  real GetShieldtubeInsideVolume();
};

// +++++++++++++++++++++++++++++++++++++++

class FRDTSubAssembly{
 protected:
  FRDTSubAssemblyGeometry *fuelsa_geometry;
  //FRDTFuelComposition *fuel_composition;
  FRDTComposition *fuel_composition;
  FRDTSUSComposition *sus_composition; // Material information for both cladding and duct
  FRDTSUSComposition *duct_composition; // Material information for duct
  FRDTSUSComposition *clad_composition; // Material information for cladding
  bool exist_fuel;
  bool exist_sus, exist_clad, exist_duct;
  bool exist_geom;
  real sodium_number_density;
  real temperature;
  real sodium_temperature;
  AtomicMassData amd;
 public:
  FRDTSubAssembly();
  void PutSubAssemblyGeometry(FRDTSubAssemblyGeometry *i){fuelsa_geometry=i; i->CheckValidity(); exist_geom=true;};
  void PutFuelComposition(FRDTComposition *i){fuel_composition=i; exist_fuel=true;};
  //void PutFuelComposition(FRDTFuelComposition *i){fuel_composition=i; exist_fuel=true;};
  //void PutSUSComposition(FRDTSUSComposition *i){sus_composition=i; exist_sus=true;};
  void PutSUSComposition(FRDTSUSComposition *i){PutDuctComposition(i); PutCladComposition(i);};
  void PutDuctComposition(FRDTSUSComposition *i){duct_composition=i; exist_duct=true;};
  void PutCladComposition(FRDTSUSComposition *i){clad_composition=i; exist_clad=true;};

  void PutSodiumDensity(real i); // g/cm3
  real GetSodiumNumberDensity(){return sodium_number_density;};

  void PutTemperature(real temp1,real temp2=0.); // non-sodium and sodium
  void PutSodiumTemperature(real temp){sodium_temperature=temp;};
  real GetTemperature(){return temperature;};
  real GetSodiumTemperature(){return sodium_temperature;};

  void CheckExistAllData();
  void ShowHomogenizedNumberDensity();
  void ShowHeterogeneousNumberDensity();
  void PutHomogenizedNumberDensity(Medium &med);

  FRDTSubAssemblyGeometry *GetFuelSAGeometry(){return fuelsa_geometry;};
  //FRDTFuelComposition *GetFuelComposition(){return fuel_composition;};
  FRDTComposition *GetFuelComposition(){return fuel_composition;};
  FRDTSUSComposition *GetSUSComposition(){return sus_composition;};
  FRDTSUSComposition *GetCladComposition(){return clad_composition;};
  FRDTSUSComposition *GetDuctComposition(){return duct_composition;};
};

class FRDTControlRod{
 protected:
  FRDTControlRodGeometry *cr_geometry;
  FRDTB4CComposition *b4c_composition;
  FRDTSUSComposition *sus_composition;
  real sodium_number_density;
  AtomicMassData amd;
 public:
  FRDTControlRod(){};
  void PutControlRodGeometry(FRDTControlRodGeometry *i){cr_geometry=i;};
  void PutB4CComposition(FRDTB4CComposition *i){b4c_composition=i;};
  void PutSUSComposition(FRDTSUSComposition *i){sus_composition=i;};

  void PutSodiumDensity(real i); // g/cm3
  real GetSodiumNumberDensity(){return sodium_number_density;};

  void ShowHomogenizedNumberDensity();
  void PutHomogenizedNumberDensity(Medium &med);
  void ShowShieldtubeInsideHomogenizedNumberDensity();
  void ShowGapSmearedCladNumberDensity();
  void ShowSpacerwireSmearedCoolantNumberDensity();
  void ShowSUSSmearedCoolantNumberDensity();
  void ShowSelf();
};

// ++++++++++++++++++++++++++++++++++++++++++++++

class FRDTCoreDataManager{
 private:
 public:
  FRDTCoreDataManager(){};
  void CylinderDataGeneration(int x, int *map, real pitch);
};

class FRDTAssypow{
 private:
  real pitch;
  real coef[10];
  real xpos[24];
  real ypos[24];
  vector< vector<real> > xyf;
  real flx[24];
 public:
  FRDTAssypow(real ipitch);
  void DoFitting();
};

class FRDTAssypow4th{
 private:
  real pitch;
  vector<real> coef;
  vector<real> xpos;
  vector<real> ypos;
  vector< vector<real> > xyf;
  vector<real> flx;
  int num_c, num_p;
 public:
  FRDTAssypow4th(real ipitch,bool xy=false);
  void PutNump(int i);
  void PutNumc(int i);
  void DoFitting();
  void PutFlux(real *inp);
  void PutFlux(int i,real inp){flx[i]=inp;};
  real GetValue(real x,real y);
  real GetFlux(int i){return flx[i];};
  real GetXpos(int i){return xpos[i];};
  real GetYpos(int i){return ypos[i];};
};

void GetAssemblyPosition(int layer,real pin_pitch,real *x,real *y);

#endif
