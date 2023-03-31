#include <cstdlib>
#include "Burner.h"

Burner::Burner():GeneralBurner()
{
  mednum=3;
  mednum_fuel=1;
  med_clad=1;
  med_water=2;
  med.resize(3);

  adj_nuc_data=false;
  hm_weight_init=0.;

  //sensitivity=true;
  sensitivity=false;

  limit_burnup=false;
  bumax=0.;

  matrix_extract=false;
  mat_ext_filename="";

  cooling_cal=false;

  collapsing=false;
  collapse_dir="";
  collapsing_micro=false;

  four_factor.resize(4);
};

void Burner::PutGeometryData(real pitch,int ring,real *rri,int *rmedi)
{
  // rri   : outer radius array from inner meshes
  // rmedi : medium ID array from inner meshes

  if(rmedi[0]!=0){
    cout<<"# Error in Burner::PutGeometryData.\n";
    cout<<"# Fuel region should be assigned to the center of the pincell.\n";
    exit(0);
  };

  int mid_max=0;
  for(int i=1;i<ring;i++){
    if(rmedi[i]>mid_max)mid_max=rmedi[i];
  };

  if(rmedi[ring-1]!=mid_max){
    cout<<"# Error in Burner::PutGeometryData.\n";
    cout<<"# Moderator region should be assigned to the most peripheral region in the ring.\n";
    exit(0);
  };

  mesh=ring+1;
  vector<real> rr(ring);

  region_medium=new int[mesh];

  fuel_r=0.;
  clad_r=0.;
  real non_mod_vol=0.;
  for(int i=0;i<ring;i++){
    region_medium[i]=rmedi[i];
    if(rmedi[i]!=0&&fuel_r==0.){
      fuel_r=rri[i-1];
      mesh_fuel=i;
      fuel_vol=fuel_r*fuel_r*PI;
    };
    if(rmedi[i]==1){
      if(clad_r!=0.){
	cout<<"# Error in Burner::PutGeometryData.\n";
	cout<<"# Cladding region CANNOT be divided into multiple spatial meshes.\n";
	exit(0);
      };
      clad_r=rri[i];
      mesh_clad=i;
      clad_vol=(clad_r*clad_r-rri[i-1]*rri[i-1])*PI;
    };
    if(rmedi[i]>1&&rmedi[i]<mid_max){
      non_mod_vol+=(rri[i]*rri[i]-rri[i-1]*rri[i-1])*PI;
    };
  };

  mod_vol=pitch*pitch-fuel_vol-clad_vol-non_mod_vol;

  

  if(clad_r==0.){
    cout<<"# Error in Burner::PutGeometryData.\n";
    cout<<"# Cladding region is not defined.\n";
    exit(0);
  };
  region_medium[ring]=region_medium[ring-1];

  real pin_pitch=pitch;

  // +++ boundary condition
  enum BCondition bc_ssc=Periodic;  // (self-shielding calculation)
  enum BCondition bc_flx=Periodic;  // (eigenvalue calculation)

  // +++ IGI for self-shielding calculation
  int mesh_ssecal=0; // The number of meshes in self-shielding calculation
  vector<real> r_ssecal;
  vector<int> mid_ssecal;
  for(int i=1;i<ring;i++){
    if(rmedi[i]!=rmedi[i-1]){
      mesh_ssecal++;
      r_ssecal.push_back(rri[i-1]);
      mid_ssecal.push_back(rmedi[i-1]);
    };
  };

  IrregularGeometryInformation igi;
  GeomPolygon pol;
  pol.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol.PutRegionID(mid_max);
  igi.AddGeom(pol);
  for(int i=0;i<mesh_ssecal;i++){
    GeomCircle cir(0.,0.,r_ssecal[mesh_ssecal-i-1]);
    cir.PutRegionID(mid_ssecal[mesh_ssecal-i-1]);
    igi.AddGeom(cir);
  };

  /*
  GeomPolygon pol;
  pol.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol.PutRegionID(2);
  igi.AddGeom(pol);
  GeomCircle cir1(0.,0.,clad_r);
  cir1.PutRegionID(1);
  igi.AddGeom(cir1);
  GeomCircle cir2(0.,0.,fuel_r);
  cir2.PutRegionID(0);
  igi.AddGeom(cir2);
  */

  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,8,0.02,45.);

  // +++ IGI for flux calculation
  IrregularGeometryInformation igi_f;
  GeomPolygon pol2;
  pol2.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol2.PutRegionID(mesh-1);
  igi_f.AddGeom(pol2);
  int *rid=new int[ring];
  for(int i=0;i<ring;i++){
    int tmp=ring-1-i;
    rr[i]=rri[tmp];
    rid[i]=tmp;
  };
  igi_f.AddCircleRing(ring,rr,rid);
  delete [] rid;

  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi_f,8,0.02,45.);

  igi_f.WriteGnuplotFile(0.001);
};

void Burner::PutGeometryDataCylinder(real pitch,int ring,real *rri,int *rmedi)
{
  mesh=ring+1;
  vector<real> rr(ring);

  mesh_fuel=0;
  mesh_clad=0;
  for(int i=0;i<ring;i++){
    if(rmedi[i]==0)mesh_fuel++;
    if(rmedi[i]==1)mesh_clad++;
  };

  region_medium=new int[mesh];
  fuel_r=0.;
  clad_r=0.;
  for(int i=0;i<ring;i++){
    region_medium[i]=rmedi[i];
    rr[ring-1-i]=rri[i];
    if(rmedi[i]==1&&fuel_r==0.){
      fuel_r=rri[i-1];
    };
    if(rmedi[i]==2&&clad_r==0.)clad_r=rri[i-1];
  };
  if(clad_r==0.)clad_r=rri[ring-1];
  region_medium[ring]=2;

  real pin_pitch=pitch;
  fuel_vol=fuel_r*fuel_r*PI;
  clad_vol=clad_r*clad_r*PI-fuel_vol;

  real outr=sqrt(pin_pitch*pin_pitch/PI);

  // +++ boundary condition
  enum BCondition bc_ssc=White;     // (self-shielding calculation)
  enum BCondition bc_flx=White;  // (eigenvalue calculation)

  // +++ IGI for self-shielding calculation
  IrregularGeometryInformation igi;
  GeomCircle pol(0.,0.,outr);
  pol.PutRegionID(2);
  igi.AddGeom(pol);
  GeomCircle cir1(0.,0.,clad_r);
  cir1.PutRegionID(1);
  igi.AddGeom(cir1);
  GeomCircle cir2(0.,0.,fuel_r);
  cir2.PutRegionID(0);
  igi.AddGeom(cir2);

  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,1,0.02,45.);

  // +++ IGI for flux calculation
  IrregularGeometryInformation igi_f;
  GeomCircle pol2(0.,0.,outr);
  pol2.PutRegionID(mesh-1);
  igi_f.AddGeom(pol2);
  int *rid=new int[ring];
  for(int i=0;i<ring;i++){
    rid[i]=ring-1-i;
  };
  igi_f.AddCircleRing(ring,rr,rid);
  delete [] rid;

  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi_f,1,0.02,45.);

  igi_f.WriteGnuplotFile(0.001);
};

void Burner::PutThreeRegionCircularGeometry(real *rri)
{
  mesh=3;

  region_medium=new int[mesh];
  for(int i=0;i<3;i++){
    region_medium[i]=i;
  };

  fuel_r=rri[0];
  clad_r=rri[1];
  mesh_fuel=1;
  mesh_clad=1;

  fuel_vol=fuel_r*fuel_r*PI;
  clad_vol=clad_r*clad_r*PI-fuel_vol;

  // +++ boundary condition
  enum BCondition bc_ssc=Periodic;     // (self-shielding calculation)
  enum BCondition bc_flx=Periodic;  // (eigenvalue calculation)

  // +++ IGI for self-shielding calculation
  IrregularGeometryInformation igi;
  GeomCircle cir1(0.,0.,rri[2]);
  cir1.PutRegionID(2);
  GeomCircle cir2(0.,0.,clad_r);
  cir2.PutRegionID(1);
  GeomCircle cir3(0.,0.,fuel_r);
  cir3.PutRegionID(0);
  igi.AddGeom(cir1);
  igi.AddGeom(cir2);
  igi.AddGeom(cir3);

  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,1,0.01,45.);

  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi,1,0.01,45.);
};

void Burner::PutThreeRegionHexagonalGeometry(real pitch,real *rri)
{
  mesh=3;

  region_medium=new int[mesh];
  for(int i=0;i<3;i++){
    region_medium[i]=i;
  };

  fuel_r=rri[0];
  clad_r=rri[1];
  mesh_fuel=1;
  mesh_clad=1;

  fuel_vol=fuel_r*fuel_r*PI;
  clad_vol=clad_r*clad_r*PI-fuel_vol;

  // +++ boundary condition
  enum BCondition bc_ssc=Periodic;     // (self-shielding calculation)
  enum BCondition bc_flx=Periodic;  // (eigenvalue calculation)

  // +++ IGI for self-shielding calculation

  IrregularGeometryInformation igi;
  real bar=pitch/sqrt(3.);
  GeomPolygon hex;
  hex.PutHexagon(0.,0.,bar);
  hex.PutRegionID(2);
  GeomCircle cir2(0.,0.,clad_r);
  cir2.PutRegionID(1);
  GeomCircle cir3(0.,0.,fuel_r);
  cir3.PutRegionID(0);
  igi.AddGeom(hex);
  igi.AddGeom(cir2);
  igi.AddGeom(cir3);

  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,12,0.01,60.);

  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi,12,0.01,60.);
};

/*
void Burner::PutWhiteBoundary()
{
  sys.PutBoundaryCondition(White);
  sys_f.PutBoundaryCondition(White);
};
*/

GroupData2D Burner::GetTransitionMatrix(Burnup &bu)
{
  int sz=total_flux.size();
  return bu.CalTransitionMatrix(total_flux[sz-1][0],true);
};

//

void Burner::CalHeavyMetalInitialWeight(Burnup &bu)
{
  hm_weight_init=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*fuel_vol; // [g] (NOT [g/cm3])
  hm_weight_init_u=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0],"U")*fuel_vol; // [g] (NOT [g/cm3])
  hm_weight_init_tru=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0],"TRU")*fuel_vol; // [g] (NOT [g/cm3])    

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  /*
    cout<<"# Power density [W/cm]\n";
    for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<power_density_list[i]<<"\n";
    };
  */
};

void Burner::CalHHMRatio()
{
  int nucn=med[0].GetNucnum();
  real den_hm=0.;
  for(int i=0;i<nucn;i++){
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    if(matno>=900000)den_hm+=med[0].GetNuclideInTurn(i).GetDensity();
  };
  //cout<<den_hm<<"\n";
  den_hm*=fuel_vol;

  int nucn2=med[2].GetNucnum();
  real den_h=0.;
  for(int i=0;i<nucn2;i++){
    int matno=med[2].GetNuclideInTurn(i).GetMatnum();
    if(matno==10010)den_h+=med[2].GetNuclideInTurn(i).GetDensity();
  };
  //cout<<mod_vol<<"\n";
  //cout<<den_h<<"\n";
  den_h*=mod_vol;

  //cout<<fuel_vol<<"\n";
  //cout<<den_hm<<"\n";
  //cout<<den_h<<"\n";
  //exit(0);

  h_hm_ratio=den_h/den_hm;  
  //cout<<"#\n# H/HM = "<<h_hm_ratio<<"\n#\n";
};

void Burner::PreCalculation(bool burn_time_cal)
{
  // +++ Pre self-shielding calculation +++++++++++++++++++++++++++++++++++++++++++++++++
  opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[1],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[2],xslib);
  //if(up_scattering)opc.CalThermalScatteringMatrix(med[2],xslib, 3.93); // 3.93 : thermal cut-off energy
  if(up_scattering)opc.CalThermalScatteringMatrix(med[2],xslib, 4.048); // ESAB from MVP library
  //if(up_scattering)opc.CalThermalScatteringMatrix(med[2],xslib, 10.);
  med[2].CalSigtr(0);

  for(int i=0;i<nucn;i++){
    nuclide_info[i]=0;
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      nuclide_info[i]=2;
      bool fissile=false;
      for(int g=0;g<group;g++){
	if(med[0].GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(g)){
	  fissile=true;
	  break;
	};
      };
      if(fissile)nuclide_info[i]=1;
    };
  };

  if(burn_time_cal)PreCalculation_bt();
};

void Burner::SelfShieldingCalculation(int step,GroupData1D &bell_default)
{
  cout<<"#\n# ... self-shielding calculation ...\n";

  if(step==0){
    ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[1],med[2],true);
    dancoff=ssc.GetDancoff(0);
  }else{
    // At the non-initial cycles, Dancoff correction calculated at the initial cycle is used.
    opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell_default,dancoff);
    //med[0].CalSigtr(0);
  };
  //med[0].CalSigtr(1);

  /*
  for(int g=0;g<group;g++){
    real e=med[0].GetEnband().get_dat(g);
    real xs1=med[0].GetNuclide(962440).GetMicxs().GetData1d(sigc).get_dat(g);
    cout<<e<<" "<<xs1<<"\n";
  };
  */

  // Thermal scattering matrices are over-written
  //if(up_scattering)opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
  if(up_scattering)opc.CalThermalScatteringMatrix(med[0],xslib,4.048);   // ESAB from MVP library
  //if(up_scattering)opc.CalThermalScatteringMatrix(med[0],xslib,10.);

  //med[0].CalMacroFromMicro();
  //med[0].CalSigtr(0);
  //med[0].CalSigtr(1);
  
  //opc.CalIncidentEnergyDependentChi(med[0],xslib);
  //opc.CalIncidentEnergyDependentChi(med[0],xslib,xslib.GetWtflux());
  opc.CalFissionSpectrumMatrix(med[0],xslib);

  cout<<"#      ... end\n";
};

void Burner::BeffCalculation(DelayedNeutronData &dnd)
{
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PJISystem lata(group,3);
  lata.PutTrajectorySet(&sys_f);
  lata.AddMedium(med[0]);
  lata.AddMedium(med[1]);
  lata.AddMedium(med[2]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  lata.CalIgenPij();

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij();

  lata.CalBetaEffective(&lat,dnd);
};

real Burner::EigenvalueCalculation(Burnup &bu, bool vid)
{
  ConstructingFuelNuclide(bu);

  Medium med2_copy=med[2];
  if(vid){
    for(int i=0;i<med[2].GetNucnum();i++){
      int mat=med[2].GetNuclideInTurn(i).GetMatnum();
      real den=med[2].GetNuclideInTurn(i).GetDensity();
      if(mat==10010||mat==80160)den*=0.9;
      med[2].GetNuclideInTurn(i).PutDensity(den);
    };
    med[2].CalMacroFromMicro();
  };
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);
  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  med[2]=med2_copy;
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  return lat.CalIgenPij(true);
};

GroupData1D Burner::GetFuelNeutronFluxInitialState(Burnup &bu, real ebnd)
{
  ConstructingFuelNuclide(bu);
  
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij(true);

  GroupData1D flx=lat.GetIntegratedFluxMeshID(0,mesh_fuel-1);

  if(ebnd>0.){
    real flx_fast=0.;
    real flx_thermal=0.;
    for(int g=0;g<group;g++){
      real e_top=med[0].GetEnband().get_dat(g);
      real e_btm=med[0].GetEnband().get_dat(g+1);
      real gflx=flx.get_dat(g);
      if(e_btm>ebnd){
        flx_fast+=gflx;
      }else if(e_top<ebnd){
	flx_thermal+=gflx;
      }else{
	real letwid=log(e_top/e_btm);
	real letwid_f=log(e_top/ebnd);
	real letwid_t=letwid-letwid_f;
	flx_fast+=gflx*letwid_f/letwid;
	flx_thermal+=gflx*letwid_t/letwid;	
      };
    };
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"#\n# Spectrum index in fuel region at initial state (cut-off energy : "<<ebnd<<" eV)\n#\n";
    cout<<"#  "<<flx_fast/flx_thermal<<"\n#\n";

  };
  
  return flx;
};

Medium Burner::GetHomogenizedCollapsedMedium(Burnup &bu, bool micro)
{
  // Hard-coded for NFI-2020 joint research
  
  ConstructingFuelNuclide(bu);
  
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);  
  lat.PutPij();
  lat.CalIgenPij(true); // Fission spectrum is treated as VECTOR

  /*
  for(int g=0;g<group;g++){
    real e0=med[0].GetEnband().get_dat(g);
    real e1=med[0].GetEnband().get_dat(g+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<lat.GetIntegratedFlux(0).get_dat(g)/letwid<<"\n";
  };
  */

  lat.PutFluxAsCurrent(); // Current is approximated by Flux
  Medium hmed=lat.HomogenizeAll(false, micro);
  GroupData1D flx_wgt=hmed.GetFlux();
  
  hmed.CalHomoB1(1e-20);  // Neutron spectrum for group-collapsing is calculated here.

  cout.setf(ios::scientific);
  cout.precision(6);
  cout<<"# K-inf in the homogenized system : "<<hmed.CalKinf()<<"\n";

  int ngrp=9;
  int bgrp[]={21,46,91,115,125,134,151,159,171};
  Medium bmed=hmed.Cond(ngrp,bgrp);
  bmed.CalHomoB1(1e-20);    // Neutron flux spectrum is re-calculated for this medium.
  bmed.CalMacroFromMicro(); // Macroscopic XS is re-calculated with neutron flux calculated above to revise fission spectrum
  GroupData1D bflx=flx_wgt.CondSum(ngrp,bgrp);
  bmed.GetMacxs().GetData1d(d).copy(bflx);

  return bmed;
};

Medium Burner::GetHomogenizedCollapsedMediumWOSSC(Burnup &bu, bool micro)
{
  // Hard-coded for NFI-2020 joint research
  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);  
  lat.PutPij();
  lat.CalIgenPij(true); // Fission spectrum is treated as VECTOR

  /*
  for(int g=0;g<group;g++){
    real e0=med[0].GetEnband().get_dat(g);
    real e1=med[0].GetEnband().get_dat(g+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<lat.GetIntegratedFlux(0).get_dat(g)/letwid<<"\n";
  };
  */

  lat.PutFluxAsCurrent(); // Current is approximated by Flux
  Medium hmed=lat.HomogenizeAll(false, micro);
  GroupData1D flx_wgt=hmed.GetFlux();
  
  hmed.CalHomoB1(1e-20);  // Neutron spectrum for group-collapsing is calculated here.

  cout.setf(ios::scientific);
  cout.precision(6);
  cout<<"# K-inf in the homogenized system : "<<hmed.CalKinf()<<"\n";

  int ngrp=9;
  int bgrp[]={21,46,91,115,125,134,151,159,171};
  Medium bmed=hmed.Cond(ngrp,bgrp);
  bmed.CalHomoB1(1e-20);    // Neutron flux spectrum is re-calculated for this medium.
  bmed.CalMacroFromMicro(); // Macroscopic XS is re-calculated with neutron flux calculated above to revise fission spectrum
  GroupData1D bflx=flx_wgt.CondSum(ngrp,bgrp);
  bmed.GetMacxs().GetData1d(d).copy(bflx);

  return bmed;
};

Medium Burner::GetHomogenizedCollapsedMedium(Burnup &bu, Medium &med_mod, bool micro)
{
  // Hard-coded for NFI-2020 joint research
  
  ConstructingFuelNuclide(bu);
  
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  //lat.AddMedium(med[2]); // 
  lat.AddMedium(med_mod);  // Only this part is different from the other [GetHomogenizedCollapsedMedium]
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij(false); // Fission spectrum is treated as VECTOR

  lat.GetMedium(2)=med[2]; // replaced by original
  
  lat.PutFluxAsCurrent(); // Current is approximated by Flux
  Medium hmed=lat.HomogenizeAll(false, micro);
  GroupData1D flx_wgt=hmed.GetFlux();
  
  hmed.CalHomoB1(1e-20);  // Neutron spectrum for group-collapsing is calculated here.
  int ngrp=9;
  int bgrp[]={21,46,91,115,125,134,151,159,171};
  Medium bmed=hmed.Cond(ngrp,bgrp);

  GroupData1D bflx=flx_wgt.CondSum(ngrp,bgrp);
  bmed.GetMacxs().GetData1d(d).copy(bflx);

  return bmed;
};




Medium Burner::GetCollapsedFuelMedium(Burnup &bu, bool micro)
{
  // Hard-coded for NFI-2020 joint research
  
  ConstructingFuelNuclide(bu);
  
  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij(true);

  lat.PutFluxAsCurrent(); // Current is approximated by Flux
  GroupData1D flx_wgt=lat.GetIntegratedFlux(0);
  
  int ngrp=9;
  int bgrp[]={21,46,91,115,125,134,151,159,171};
  Medium bmed=med[0].Cond(ngrp,bgrp,flx_wgt,flx_wgt);

  return bmed;
};

void Burner::Cal1groupBranchingRatio(string cbglibdir,string filename,string outname)
{
  // Eigenvalue calculation to obtain neutron energy spectrum

  PreCalculation(false);
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  SelfShieldingCalculation(0,b);

  GeneralOption opt;

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij();

  GroupData1D flx=lat.GetIntegratedFluxMeshID(0,mesh_fuel-1);
  GroupData1D dat(group);
  int bgrp[]={0};
  bgrp[0]=group-1;

  ifstream fin;
  string libdir=cbglibdir+"CBGLIB_BURN/CBG_Chain/";
  libdir.append(filename);
  fin.open(libdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<libdir<<"\n";
    exit(0);
  };

  int tmp;
  fin>>tmp;
  if(tmp!=group){
    cout<<"# Error in Burner::OverWritingNGBranchingRatioData\n";
    cout<<"# The number of groups is inconsistent.\n";
    exit(0);
  };

  ofstream fout;
  fout.open(outname.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file in Burner::Cal1groupNGBranchingRatio.\n";
    exit(1);
  };

  int mat=999;
  while(mat!=-1){
    fin>>mat;
    if(mat==-1){
      fout<<"-1";
      fin.close();
      fout.close();
      return;
    };

    int flag;
    fin>>flag;

    int ich;
    fin>>ich;

    bool process=true;
    GroupData1D ngxs;
    if(xslib.ExistLibData(mat)){
      if(flag==0)ngxs=xslib.GetLibData(mat).GetXSData().GetData1d(sigc);
      if(flag==1)ngxs=xslib.GetLibData(mat).GetXSData().GetData1d(sign2n);
      real sum=ngxs.get_sum();
      if(sum==0.)process=false;
    }else{
      process=false;
    };

    if(process){
      fout.setf(ios::showpoint);
      fout.precision(9);
      ngxs=ngxs.mult(flx);
      fout<<"   "<<mat<<"\n";
      fout<<"   "<<flag<<"\n";
      fout<<"   "<<ich<<"\n";
      vector<real> br(ich);
      for(int i=0;i<ich;i++){
	fin>>tmp;
        for(int j=0;j<group;j++){
	  real tmpr;
	  fin>>tmpr;
	  dat.put_data(j,tmpr);
	};
        br[i]=dat.Cond(ngxs,1,bgrp).get_dat(0);
      };
      if(flag==1){
	real sum=0.;
	for(int i=0;i<ich;i++){
	  sum+=br[i];
	};
	for(int i=0;i<ich;i++){
	  br[i]/=sum;
	};
      };
      for(int i=0;i<ich;i++){
	fout<<"   "<<br[i]<<"\n";
      };
    }else{
      for(int i=0;i<ich;i++){
        fin>>tmp;
        for(int j=0;j<group;j++){
          real tmpr;
          fin>>tmpr;
	};
      };
    };
  };

  fin.close();
};

void Burner::ConstructingFuelNuclide(Burnup &bu)
{
  //
  // Instances of the nuclide class are stored in acsending MATID order
  // to perform efficient burnup calculation with CRAM.
  //
  bu.AddNuclideToMediumFromBurnupChain(med[0]);
  PutNucn(med[0].GetNucnum());
};

void Burner::PutNucn(int i)
{
  nucn=i;
  nuclide_info.resize(nucn);
};

void Burner::Calculation(Burnup &bu,bool adjoint)
{
  setlinebuf(stdout);

  ConstructingFuelNuclide(bu);
  if(hm_weight_init==0.)CalHeavyMetalInitialWeight(bu);
  CalHHMRatio();

  PreCalculation();

  //
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  GeneralOption opt;

  // +++ Number density data storing ++++++++++++++++++++++++++++++++++++++++++++++++++++
  density_data.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  fuel_flux.resize(burn_step+1);
  clad_flux.resize(burn_step+1);

  for(int i=0;i<4;i++){
    four_factor[i].resize(burn_step+1);
  };

  power_factor.resize(burn_step);
  delt.resize(burn_step);
  total_flux.resize(burn_step);
  total_flux_unnormalized.resize(burn_step);
  total_flux_clad.resize(burn_step);
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    total_flux[i].resize(sub_step);
    total_flux_unnormalized[i].resize(sub_step);
    total_flux_clad[i].resize(sub_step);
  };

  Burnup bu_dmdf;

  /*
  // +++ Preparing to plot
  ofstream fout;
  fout.open("./out_flx",ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    exit(1);
  ;}
  */


  // +++ for sensitivity calculation +++++++++++++++++++++++++++++++++++++++
  if(adjoint){
    bu_dmdf=bu;
    SetArrayForSensitivityCalculation();
    adj_nuc_data=true;
  };
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  cout<<"#\n";
  
  real bu0=0.;
  real accumulated_day=0.;
  real accumulated_burn=0.;
  for(int st=0;st<burn_step+1;st++){

    //acday[st]=accumulated_day;
    //acburn[st]=accumulated_burn;
    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);

    cout<<"#\n# +++ Burnup step : "<<st<<" ";
    cout<<" (accumulated day : "<<accumulated_day<<" )\n";

    density_data[st].resize(nucn);
    for(int i=0;i<nucn;i++){
      density_data[st][i]=med[0].GetNuclideInTurn(i).GetDensity();
    };

    SelfShieldingCalculation(st,b);

    //med[0].ShowMacroCrossSection1D();

    /*
    real org=med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc).get_dat(56);
    med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc).add_data(56,org*0.001);
    med[0].GetNuclide(922380).GetMicxs().GetData1d(sigt).add_data(56,org*0.001);
    med[0].GetNuclide(922380).GetMicxs().GetData1d(sigt,1).add_data(56,org*0.001);
    med[0].CalMacroFromMicro();
    */

    // +++ Eigenvalue calculation
    PJISystem lat(group,3);
    lat.PutTrajectorySet(&sys_f);
    lat.AddMedium(med[0]);
    lat.AddMedium(med[1]);
    lat.AddMedium(med[2]);
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutSigmaCol(sigtr);
    lat.PutPij();
    keff[st]=lat.CalIgenPij(true); // `true' means matrix chi

    /*
    cout.setf(ios::scientific);
    cout.precision(6);
    real totrr=0.;
    for(int g=0;g<group;g++){
      real en=med[0].GetEnband().get_dat(g);
      real en1=med[0].GetEnband().get_dat(g+1);
      real letwid=log(en/en1);
      real xs=med[0].GetNuclide(551330).GetMicxs().GetData1d(sigc).get_dat(g);
      //real xs=med[1].GetMacxs().GetData1d(sign2n).get_dat(g);      
      real flx1=lat.GetIntegratedFlux(0).get_dat(g);
      totrr+=xs*flx1;
      cout<<en<<" "<<xs<<" "<<xs*flx1/letwid<<" "<<totrr<<"\n";
    };
    exit(0);
    */
    /*
    cout.setf(ios::scientific);
    cout.precision(6);
    real sum_prod=0.;
    real sum_absf=0.;
    real sum_absc=0.;
    real sum_absm=0.;        
    for(int g=0;g<group;g++){
      real flx0=lat.GetIntegratedFlux(0).get_dat(g);
      real flx1=lat.GetIntegratedFlux(1).get_dat(g);
      real flx2=lat.GetIntegratedFlux(2).get_dat(g);
      sum_prod+=flx0*med[0].GetMacxs().GetData1d(nusigf).get_dat(g);
      sum_absf+=flx0*med[0].GetMacxs().GetData1d(siga).get_dat(g);
      sum_absc+=flx1*med[1].GetMacxs().GetData1d(siga).get_dat(g);
      sum_absm+=flx2*med[2].GetMacxs().GetData1d(siga).get_dat(g);            
    };
    cout<<sum_prod<<"\n";
    cout<<sum_absf<<"\n";
    cout<<sum_absc<<"\n";
    cout<<sum_absm<<"\n";
    cout<<sum_prod/(sum_absf+sum_absc+sum_absm)<<"\n";
    exit(0);
    */

    /*
    // +++ Adjoint calculation
    GeneralOption opta;
    opta.PutAdjointCal();

    PJISystem lata(group,3);
    lata.PutTrajectorySet(&sys_f);
    lata.AddMedium(med[0]);
    lata.AddMedium(med[1]);
    lata.AddMedium(med[2]);
    lata.PutRegMed(region_medium);
    lata.PutGeneralOption(opta);
    lata.PutSigmaCol(sigtr);
    lata.PutPij();
    real kdum=lata.CalIgenPij();
    */

    /*
    // +++ Flux plotting
    for(int i=0;i<group;i++){
      real e0=med[0].GetEnband().get_dat(i);
      real e1=med[0].GetEnband().get_dat(i+1);
      real letwid=log(e0/e1);
      fout<<e0<<" ";
      fout<<lat.GetIntegratedFlux(0).get_dat(i)/letwid<<"\n";
    };
    fout<<"\n\n";
    */

    // --- (Generating medium instances include collapsing cross section data) 
    if(collapsing){
      lat.PutFluxAsCurrent(); // Current is approximated by Flux
      //Medium hmed=lat.HomogenizeAll(false,collapsing_micro); // [benoist_d, micxs]
      Medium hmed=lat.HomogenizeAll(false,true); // [benoist_d, micxs]      
      hmed.CalHomoB1(1e-20);  // Neutron spectrum for group-collapsing is calculated here.
      cout.setf(ios::scientific);
      cout.precision(6);
      cout<<"# K-inf in the homogenized system : "<<hmed.CalKinf()<<"\n";
      int ngrp=9;
      int bgrp[]={21,46,91,115,125,134,151,159,171};
      Medium bmed=hmed.Cond(ngrp,bgrp);
      bmed.CalHomoB1(1e-20);    // Neutron flux spectrum is re-calculated for this medium.
      bmed.CalMacroFromMicro(); // Macroscopic XS is re-calculated with neutron flux calculated above to revise fission spectrum
      string filename=collapse_filename+"_st"+IntToString(st);
      bmed.WriteFile(collapse_dir,filename,collapsing_micro);
    };
    // --------------------------------------------------------------------------------



    // +++ four-factor calculation ++++++++++++++++++++++++
    /*
    real nusigf_t=lat.GetIntegratedReactionRate(nusigf,0).get_sum();
    real siga_t=lat.GetIntegratedReactionRate(siga).get_sum();
    GroupData1D mic_nusigf_u8=med[0].GetNuclide(922380).GetMicxs().GetData1d(sigf).mult(med[0].GetNuclide(922380).GetMicxs().GetData1d(nu));
    GroupData1D mic_nusigf_p0=med[0].GetNuclide(942400).GetMicxs().GetData1d(sigf).mult(med[0].GetNuclide(942400).GetMicxs().GetData1d(nu));
    real nusigf_u8=(lat.GetIntegratedFlux(0)*mic_nusigf_u8)*med[0].GetNuclide(922380).GetDensity();
    real nusigf_p0=(lat.GetIntegratedFlux(0)*mic_nusigf_p0)*med[0].GetNuclide(942400).GetDensity();
    real sigf_u8=(lat.GetIntegratedFlux(0)*med[0].GetNuclide(922380).GetMicxs().GetData1d(sigf))*med[0].GetNuclide(922380).GetDensity();
    real sigc_u8=(lat.GetIntegratedFlux(0)*med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc))*med[0].GetNuclide(922380).GetDensity();
    //real sigc_o=(lat.GetIntegratedFlux(0)*med[0].GetNuclide(80160).GetMicxs().GetData1d(sigc))*med[0].GetNuclide(80160).GetDensity();
    real siga_clad=lat.GetIntegratedReactionRate(siga,1).get_sum();
    real siga_mod=lat.GetIntegratedReactionRate(siga,2).get_sum();
    real siga_fp=0.;
    int nucnum=med[0].GetNucnum();
    for(int i=0;i<nucnum;i++){
      int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
      if(nucid>10000&&nucid<900000&&nuclide_info[i]!=0){
        siga_fp+=(lat.GetIntegratedFlux(0)*med[0].GetNuclideInTurn(i).GetMicxs().GetData1d(sigc))*med[0].GetNuclideInTurn(i).GetDensity();
      };
    };
    four_factor[0][st]=nusigf_t/(nusigf_t-nusigf_u8-nusigf_p0);
    four_factor[1][st]=(siga_t-sigf_u8-sigc_u8)/siga_t;
    four_factor[2][st]=(siga_t-sigf_u8-sigc_u8-siga_clad-siga_mod-siga_fp)/(siga_t-sigf_u8-sigc_u8);
    four_factor[3][st]=(nusigf_t-nusigf_u8-nusigf_p0)/(siga_t-sigf_u8-sigc_u8-siga_clad-siga_mod-siga_fp);
    */
    /*
    cout.setf(ios::showpoint);
    cout.precision(5);
    cout<<four_factor[0][st]<<" ";
    cout<<four_factor[1][st]<<" ";
    cout<<four_factor[2][st]<<" ";
    cout<<four_factor[3][st]<<" ";
    cout<<four_factor[0][st]*four_factor[1][st]*four_factor[2][st]*four_factor[3][st]<<" ";
    cout<<"\n";
    */
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++

    real vol_inv=1./fuel_vol;
    fuel_flux[st]=lat.GetIntegratedFluxMeshID(0,mesh_fuel-1)*vol_inv;
    //clad_flux[st]=lat.GetIntegratedFluxMeshID(mesh_fuel,mesh_fuel+mesh_clad-1)/clad_vol;
    clad_flux[st]=lat.GetIntegratedFluxMeshID(mesh_clad,mesh_clad)/clad_vol;

    // Total-flux per unit volume in fuel region
    med[0].GetFlux().copy(fuel_flux[st]);

    opc.CalEquivalentFissionSpectrumVector(med[0],fuel_flux[st]);

    //real org=med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc).get_dat(56);
    //med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc).put_data(56,org*1.001);
    bu.PutMediumData(med[0]); // burnup data
    //med[0].GetNuclide(922380).GetMicxs().GetData1d(sigc).put_data(56,org);

    // +++ one-group cross section storing
    xsc_1g[st].resize(nucn);
    xsn2n_1g[st].resize(nucn);
    xsf_1g[st].resize(nucn);
    for(int j=0;j<nucn;j++){
      xsf_1g[st][j]=bu.GetSigf(j);
      xsc_1g[st][j]=bu.GetSigc(j);
      xsn2n_1g[st][j]=bu.GetSign2n(j);
    };

    if(limit_burnup){
      if(bu0==0.){
	if(keff[st]<k_limit){
          bu0=acburn[st-1]+(acburn[st]-acburn[st-1])/(keff[st-1]-keff[st])*(k_limit-keff[st]);
	};
      }else{
        bumax=bu0*2*batch/(batch+1);
	if(acburn[st]>bumax){
          cout<<"#\n# Burnup calculation is forced to be terminated in a limited burnup.\n#\n";
          burn_step=st;
          return;
	};
      };
    };

    if(adjoint){
      for(int g=0;g<group;g++){
	pij_store[st][g].copy(lat.GetPij(g));
      };
    };

    if(st!=burn_step){
      // +++ (for sensitivity calculation) +++++++++++++++++++++++++
      if(adjoint){
        macxs[st].DataCopyPL(med[0].GetMacxs(),0);
        if(st==0)trmat_flxindep=bu.GetTrmatFlxInDep();
        for(int i=0;i<mesh_fuel;i++){
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*lat.GetMesh(i).GetVolume();
        };
        for(int k=0;k<nucn;k++){
	  if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1)mic_sigf[st][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            mic_sigc[st][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n[st][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
	    /*
	    GroupData2D tmp=med[0].GetNuclideInTurn(k).GetMicxs().GetData2d(sigel)
	                   +med[0].GetNuclideInTurn(k).GetMicxs().GetData2d(siginel)
	      +med[0].GetNuclideInTurn(k).GetMicxs().GetData2d(sign2n);
	    mic_sigs[st][k].copy(tmp);
	    */
	  };
        };
      };
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day

      int sub_step=sub_step_list[st];
      burn_span/=sub_step;
    
      // +++ Burnup calculation
      if(!input_flux_level)accumulated_burn+=burn_time_gwd[st];
      //cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
      real flx1g=fuel_flux[st].get_sum();
      real flx1g_clad=clad_flux[st].get_sum();

      for(int j=0;j<sub_step;j++){
	//cout<<"#... burnup calculation : "<<j<<"/"<<sub_step<<"\n";   

	// +++ (for sensitivity calculation) ++++++++++++++++++++++++++++++
	if(adjoint){
          for(int i=0;i<nucn;i++){
	    fwd_nuc[st][j].put_data(i,med[0].GetNuclideInTurn(i).GetDensity());
          };
	};
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        power_factor[st][j]=power_density/(bu.GetIntegratedPower(med[0])*fuel_vol);
        total_flux[st][j]=flx1g*power_factor[st][j];
        total_flux_unnormalized[st][j]=flx1g;
        total_flux_clad[st][j]=flx1g_clad*power_factor[st][j];
	if(input_flux_level){
          real ratio=flx1g_clad/flx1g;
          total_flux[st][j]=flux_level_list[st];
	  power_factor[st][j]=total_flux[st][j]/flx1g;
	  total_flux_clad[st][j]=total_flux[st][j]*ratio;
          real power_org=bu.GetIntegratedPower(med[0])*fuel_vol;
          power_density=power_org*total_flux[st][j]/flx1g;
          accumulated_burn+=(power_density*1e-9)/(hm_weight_init*1e-6)*burn_span;
	};
	accumulated_day+=burn_span;
	delt[st][j]=burn_span*24*60*60;
	// +++ (for burnup-matrix extraction)+++++++++++++++
	if(j==0&&matrix_extract){
	  ofstream fout;
	  string tmp="./"+mat_ext_filename+"_step"+IntToString(st);
	  fout.open(tmp.data(),ios::out);
	  if(fout.fail()){
	    cout<<"# Failed to open the file.\n";
	    exit(1);
	  };
	  fout.setf(ios::scientific);
	  fout.precision(mat_ext_prc);
          fout<<nucn<<"\n";
	  for(int i=0;i<nucn;i++){
	    int mn=med[0].GetNuclideInTurn(i).GetMatnum();
	    fout<<mn<<"\n";
  	    fout<<"   "<<midt.Name(mn)<<"\n";
	  };
	  GroupData2D mat=bu.CalTransitionMatrix(total_flux[st][j]);
	  for(int i=0;i<nucn;i++){
	    for(int j=0;j<nucn;j++){
	      fout<<mat.get_dat(i,j)<<"\n";
	    };
	  };
	  fout.close();

	};
	// ++++++++++++++++++++++++++++++++++++++++++++++++++

	if(adjoint==false){
	  bu.BurnupCalculationByMMPA(med[0],total_flux[st][j],delt[st][j],false);
	  //bu.BurnupCalculationByChebyshev(med[0],total_flux[st][j],delt[st][j],false);
	}else if(adjoint==true){
	  bu.BurnupCalculationByMultiStepCalc(med[0],fwd_nuc_int[st][j],total_flux[st][j],delt[st][j],false);
	};

      };

      // +++ (for sensitivity calculation)
      if(adjoint){
        for(int i=0;i<nucn;i++){
  	  fwd_nuc[st][sub_step].put_data(i,med[0].GetNuclideInTurn(i).GetDensity());
        };

        for(int jj=0;jj<group;jj++){
          GroupData2D tmp=bu_dmdf.CaldTMatdFlux(med[0],jj);
          for(int j=0;j<sub_step+1;j++){
            dmdf_nuc[st][j][jj]=tmp*fwd_nuc[st][j];
	  };
  	};
      };

      // ++++++++++++++++++++++++++++

      //cout<<"#    ... terminated.\n";
    };
  };
};

void Burner::dkdN_Calculation(Burnup &bu,string filename)
{
  setlinebuf(stdout);

  ConstructingFuelNuclide(bu);
  if(hm_weight_init==0.)CalHeavyMetalInitialWeight(bu);
  PreCalculation();
  //
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  power_factor.resize(burn_step);
  delt.resize(burn_step);
  total_flux.resize(burn_step);
  fuel_flux.resize(burn_step+1);
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    total_flux[i].resize(sub_step);
  };

  // +++ dkdN_storing
  vector< vector<real> > dkdN(burn_step+1);
  for(int i=0;i<burn_step+1;i++){
    dkdN[i].resize(nucn,0.);
  };

  real accumulated_day=0.;
  real accumulated_burn=0.;
  for(int st=0;st<burn_step+1;st++){

    //acday[st]=accumulated_day;
    //acburn[st]=accumulated_burn;
    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    SelfShieldingCalculation(st,b);

    med[0].CalMacroFromMicroSimple();

    // +++ Eigenvalue calculation
    PJISystem lata(group,3);
    lata.PutTrajectorySet(&sys_f);
    lata.AddMedium(med[0]);
    lata.AddMedium(med[1]);
    lata.AddMedium(med[2]);
    lata.PutRegMed(region_medium);
    lata.PutGeneralOption(opta);
    lata.PutSigmaCol(sigtr);
    lata.PutPij();
    real k_adj=lata.CalIgenPij();

    PJISystem lat(group,3);
    lat.PutTrajectorySet(&sys_f);
    lat.AddMedium(med[0]);
    lat.AddMedium(med[1]);
    lat.AddMedium(med[2]);
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutSigmaCol(sigtr);
    lat.PutPij();
    keff[st]=lat.CalIgenPij();
 
    for(int i=0;i<nucn;i++){
      cout<<"  dkdn calculation : "<<i<<"/"<<nucn<<"\n";
      real org=med[0].GetNuclideInTurn(i).GetDensity();
      real sns=0.;
      if(org!=0.){
        lat.GetMedium(0).GetNuclideInTurn(i).PutDensity(org*1.01);
        lat.GetMedium(0).CalMacroFromMicroSimple();
        sns=lata.CalReactivity(&lat,k_adj,keff[st],false)*k_adj*100.;
        lat.GetMedium(0).GetNuclideInTurn(i).PutDensity(org);
      };
      dkdN[st][i]=sns;
    };

    real vol_inv=1./fuel_vol;
    fuel_flux[st]=lat.GetIntegratedFluxMeshID(0,mesh_fuel-1)*vol_inv;
    // Total-flux per unit volume in fuel region
    med[0].GetFlux().copy(fuel_flux[st]);

    bu.PutMediumData(med[0]); // burnup data

    if(st!=burn_step){

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day

      int sub_step=sub_step_list[st];
      burn_span/=sub_step;
    
      // +++ Burnup calculation
      if(!input_flux_level)accumulated_burn+=burn_time_gwd[st];
      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
      real flx1g=fuel_flux[st].get_sum();
      for(int j=0;j<sub_step;j++){
        power_factor[st][j]=power_density/(bu.GetIntegratedPower(med[0])*fuel_vol);
        total_flux[st][j]=flx1g*power_factor[st][j];
	if(input_flux_level){
          total_flux[st][j]=flux_level_list[st];
          real power_org=bu.GetIntegratedPower(med[0])*fuel_vol;
          power_density=power_org*total_flux[st][j]/flx1g;
          accumulated_burn+=(power_density*1e-9)/(hm_weight_init*1e-6)*burn_span;
	};
	accumulated_day+=burn_span;
	delt[st][j]=burn_span*24*60*60;
        bu.BurnupCalculationByChebyshev(med[0],total_flux[st][j],delt[st][j],false);
      };

      cout<<"#    ... terminated.\n";
    };
  };

  // +++ Editing
  ofstream fout;
  string fname="./"+filename;
  fout.open(fname.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file in Burner::dkdN_Calculation.\n";
    exit(1);
  };

  fout.setf(ios::showpoint);
  fout.precision(5);

  int judge_num=7;
  real judge[]={1., 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
  for(int jj=0;jj<judge_num;jj++){
    fout<<"#\n# Maximum (dk/k)/(dn/n) through burnup\n";
    fout<<"#  [criteria : "<<judge[jj]<<"]\n#\n";
    for(int i=0;i<nucn;i++){
      int mmt=med[0].GetNuclideInTurn(i).GetMatnum();
      real maxval=0.;
      for(int j=0;j<burn_step+1;j++){
        if(abs(dkdN[j][i])>maxval)maxval=abs(dkdN[j][i]);
      };
      if(maxval>judge[jj]){
        fout<<midt.Name(mmt)<<" "<<maxval<<"\n";
      };
    };
    fout<<"\n\n";
  };

  fout<<"#\n";
  fout<<"# Burnup-dependent keff sensitivity\n";
  fout<<"#\n";
  for(int i=0;i<burn_step+1;i++){
    fout<<"# Burnup step : "<<i<<"\n";
    fout<<"# Accumulated day : "<<acday[i]<<"\n";
    fout<<"# Burnup : "<<acburn[i]<<"\n";
    fout<<"#\n";
    for(int j=0;j<nucn;j++){
      if(dkdN[i][j]!=0.){
        fout<<"  "<<midt.Name(med[0].GetNuclideInTurn(j).GetMatnum())<<" ";
        fout<<dkdN[i][j]<<"\n";
      };
    };
    fout<<"#\n#\n";
  };

  vector<real> dkdn_max(nucn,0.);
  fout<<"#\n# Burnup dependent (dk/k)/(dn/n) per each nuclide\n#\n";
  for(int i=0;i<nucn;i++){
    int mmt=med[0].GetNuclideInTurn(i).GetMatnum();
    fout<<"# "<<midt.Name(mmt)<<"\n";
    fout<<"# [day] [GWd/t] [(dk/k)/(dn/n)]\n";
    for(int j=0;j<burn_step+1;j++){
      fout<<acday[j]<<" "<<acburn[j]<<" "<<dkdN[j][i]<<"\n";
      if(fabs(dkdN[j][i])>dkdn_max[i])dkdn_max[i]=dkdN[j][i];
    };
    fout<<"\n\n";
  };

  fout<<"# Maximum dk/dn for each nuclide\n";
  for(int i=0;i<nucn;i++){
    fout<<"#  "<<midt.Name(med[0].GetNuclideInTurn(i).GetMatnum())<<" ";
    fout<<dkdn_max[i]<<"\n";
  };

};

void Burner::NuclideWiseReactivityEffectCalculationDuringCooling(Burnup &bu, int nuc, string *nuc_name)
{
  // substep is fixed to 1.
  int sub_step=1;

  vector< vector<real> > react(nuc);
  vector<int> nuc_turn(nuc);
  GetPrintNuclide(nuc,nuc_name,nuc_turn);

  //ConstructingFuelNuclide(bu); 
  for(int i=0;i<burn_step;i++){
    power_density_list[i]=0.;
  };
  PreCalculation_bt();

  // forward/adjoint calculation at t=0
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  med[0].CalMacroFromMicro();
  Medium med_org=med[0];

  PJISystem lata(group,3);
  lata.PutTrajectorySet(&sys_f);
  lata.AddMedium(med_org);
  lata.AddMedium(med[1]);
  lata.AddMedium(med[2]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  real k_adj=lata.CalIgenPij();

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med_org);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real k_fwd=lat.CalIgenPij();

  // +++ Number density data storing ++++++++++++++++++++++++++++++++++++++++++++++++++++
  delt.resize(burn_step);
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
  };

  acday.clear();
  acburn.clear();
  density_data.clear();

  real accumulated_day=0.;
  for(int st=0;st<burn_step+1;st++){

    acday.push_back(accumulated_day);
    acburn.push_back(0.);

    vector<real> density_tmp(nucn);
    for(int i=0;i<nucn;i++){
      density_tmp[i]=med[0].GetNuclideInTurn(i).GetDensity();
    };
    density_data.push_back(density_tmp);

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    bu.PutMediumDataWithZeroCrossSection(med[0]); // burnup data

    if(st!=burn_step){

      real burn_span=burn_time[st]; // day
      burn_span/=sub_step;
      // +++ Burnup calculation
      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
      for(int j=0;j<sub_step;j++){
	accumulated_day+=burn_span;
	delt[st][j]=burn_span*24*60*60;
        bu.BurnupCalculationByChebyshev(med[0],0.,delt[st][j],false);
      };

      // perturbation calculation
      for(int i=0;i<nuc;i++){
        real den_org=med_org.GetNuclideInTurn(nuc_turn[i]).GetDensity();
        real den_new=med[0].GetNuclideInTurn(nuc_turn[i]).GetDensity();
        med_org.GetNuclideInTurn(nuc_turn[i]).PutDensity(den_new);
        med_org.CalMacroFromMicro();
        lat.GetMedium(0).GetMacxs().DataCopy(med_org.GetMacxs());
        real tmp=lata.CalReactivity(&lat,k_adj,k_fwd,false);
        react[i].push_back(tmp);
        med_org.GetNuclideInTurn(nuc_turn[i]).PutDensity(den_org);
      };
      cout<<"#    ... terminated.\n";
    };
  };

  cout<<"# Nuclide-wise reactivity effect during cooling period.\n";
  cout<<"#  Day     ";
  for(int i=0;i<nuc;i++){
    WriteOut(nuc_name[i],11);
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<burn_step;i++){
    cout<<acday[i+1]<<" ";
    for(int j=0;j<nuc;j++){
      if(react[j][i]>=0.)cout<<" ";
      cout<<react[j][i]<<" ";
    };
    cout<<"\n";
  };

};

void Burner::CoolingCalculation(Burnup &bu,int sub_step)
{

  ConstructingFuelNuclide(bu); 
  for(int i=0;i<burn_step;i++){
    power_density_list[i]=0.;
  };
  PreCalculation_bt();
  
  // +++ Number density data storing ++++++++++++++++++++++++++++++++++++++++++++++++++++
  //density_data.resize(burn_step+1);
  delt.resize(burn_step);
  for(int i=0;i<burn_step;i++){
    //int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
  };

  if(!cooling_cal){
    acday.clear();
    acburn.clear();
    density_data.clear();
  };

  real accumulated_day=0.;
  if(cooling_cal){
    int sz=acday.size();
    accumulated_day=acday[sz-1];
  };

  for(int st=0;st<burn_step+1;st++){

    if(!cooling_cal||st!=0){
      acday.push_back(accumulated_day);
      acburn.push_back(0.);

      vector<real> density_tmp(nucn);
      for(int i=0;i<nucn;i++){
        density_tmp[i]=med[0].GetNuclideInTurn(i).GetDensity();
      };
      density_data.push_back(density_tmp);
    };

    //cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    //cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    //density_data[st].resize(nucn);
    //for(int i=0;i<nucn;i++){
    //density_data[st][i]=med[0].GetNuclideInTurn(i).GetDensity();
    //};

    bu.PutMediumDataWithZeroCrossSection(med[0]); // burnup data

    if(st!=burn_step){

      real burn_span=burn_time[st]; // day
      burn_span/=sub_step;
      // +++ Burnup calculation
      //cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
      for(int j=0;j<sub_step;j++){
	accumulated_day+=burn_span;
	delt[st][j]=burn_span*24*60*60;
        //bu.BurnupCalculationByKrylov(med[0],0.,delt[st][j],false);
        //bu.BurnupCalculationByChebyshev(med[0],0.,delt[st][j],false);
        bu.BurnupCalculationByMMPA(med[0],0.,delt[st][j],false);
      };

      //cout<<"#    ... terminated.\n";
    };
  };

  cooling_cal=true;
  int sz=acday.size();
  burn_step=sz-1;
};

real Burner::CalMacroscopicReactionRate(int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh)
{
  GeneralOption opt;
  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  lat.CalIgenPij();
  return lat.CalMacroscopicReactionRate(nuc,nuc_id,xs,on_mesh);
};

void Burner::CalculationHomogeneous(Burnup &bu)
{
  fuel_vol=1.;
  setlinebuf(stdout);

  //if(hm_weight_init==0.)hm_weight_init=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*fuel_vol; // [g] (NOT [g/cm3])
  ConstructingFuelNuclide(bu);
  if(hm_weight_init==0.)CalHeavyMetalInitialWeight(bu);
  PreCalculation_bt();

  // +++ Number density data storing ++++++++++++++++++++++++++++++++++++++++++++++++++++
  density_data.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  fuel_flux.resize(burn_step+1);

  real accumulated_day=0.;
  real accumulated_burn=0.;
  for(int st=0;st<burn_step+1;st++){

    //acday[st]=accumulated_day;
    //acburn[st]=accumulated_burn;
    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    density_data[st].resize(nucn);
    for(int i=0;i<nucn;i++){
      density_data[st][i]=med[0].GetNuclideInTurn(i).GetDensity();
    };

    opc.CalSelfShieldingInfiniteSystem(med[0],xslib);
    med[0].CalMacroFromMicro();

    med[0].CalHomoB1(0.);
    keff[st]=med[0].CalKeff();
    fuel_flux[st]=med[0].GetFlux();

    bu.PutMediumData(med[0]); // burnup data

    // +++ one-group cross section printing
    xsc_1g[st].resize(nucn);
    xsn2n_1g[st].resize(nucn);
    xsf_1g[st].resize(nucn);
    for(int j=0;j<nucn;j++){
      xsf_1g[st][j]=bu.GetSigf(j);
      xsc_1g[st][j]=bu.GetSigc(j);
      xsn2n_1g[st][j]=bu.GetSign2n(j);
    };

    if(st!=burn_step){

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day

      int sub_step=sub_step_list[st];
      burn_span/=sub_step;
    
      // +++ Burnup calculation
      if(!input_flux_level)accumulated_burn+=burn_time_gwd[st];
      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
      real flx1g=med[0].GetFlux().get_sum();
      for(int j=0;j<sub_step;j++){
	//cout<<"#... burnup calculation : "<<j<<"/"<<sub_step<<"\n";   
        real power_factor=power_density/bu.GetIntegratedPower(med[0]);
        real total_flux=flx1g*power_factor;
	if(input_flux_level){
          total_flux=flux_level_list[st];
          real power_org=bu.GetIntegratedPower(med[0]);
          power_density=power_org*total_flux/flx1g;
          accumulated_burn+=(power_density*1e-9)/(hm_weight_init*1e-6)*burn_span;
	};
	accumulated_day+=burn_span;
	real delt=burn_span*24*60*60;
        bu.BurnupCalculationByKrylov(med[0],total_flux,delt,false);
	cout<<"#  Flux level    : "<<total_flux<<" [/s/cm2]\n";
	//cout<<"#  Power density : "<<power_density<<" [W/cm]\n";
      };
      cout<<"#    ... terminated.\n";
    };
  };
};

void Burner::SensitivityCalculation(Burnup &bu,int target_nuclide_num,string *target_nuclide_snscal,bool fileout, string prefilename)
{
  sensitivity=true;
  if(prefilename=="")prefilename="sns.";

  // +++ Forward calculation
  Calculation(bu,true);

  // +++ for forward number density printing +++++++++++++++++++++++
  if(fileout)WriteFileNumberDensity("nuc_fwd");

  // (Sensitivity Calculation)
  //AveragingForwardNumberDensity();
  IntegratingForwardNumberDensity();
  for(int iii=0;iii<target_nuclide_num;iii++){

    cout<<"# Sensitivity calculation : "<<iii<<"/"<<target_nuclide_num<<"\n";
    int sub_step=sub_step_list[burn_step-1];
    // (initialization for final condition)
    for(int i=0;i<nucn;i++){
      adj_nuc_e[burn_step-1][sub_step-1].put_data(i,0.);
    };

    // (Set-up for final condition)
    int target_nuclideID=midt.ID(target_nuclide_snscal[iii]);
    int target_nuclideID_inturn=-1;
    for(int i=0;i<nucn;i++){
      if(med[0].GetNuclideInTurn(i).GetMatnum()==target_nuclideID){
        target_nuclideID_inturn=i;
        adj_nuc_e[burn_step-1][sub_step-1].put_data(i,1.);
      };
    };
    if(target_nuclideID_inturn==-1){
      cout<<"# Error in Burner::Sensitivity Calculation.\n";
      cout<<"# Targert nuclide "<<target_nuclide_snscal[iii]<<" cannot be found.\n";
      exit(0);
    };

    AdjointBurnupCalculation(bu);

    real end_nuc=med[0].GetNuclideInTurn(target_nuclideID_inturn).GetDensity();
    target_val=end_nuc;
    SensitivityPrinting(bu,end_nuc,"burner_pincell","number_density","unknown",prefilename+target_nuclide_snscal[iii]);

    // +++ adjoint number density printing +++++++++++++++++++++++++++
    if(fileout){
      string fname="nuc_adj_"+target_nuclide_snscal[iii];
      WriteFileAdjointNumberDensity(fname,target_nuclideID);
      string fname2="cf_"+target_nuclide_snscal[iii];
      WriteFileContributionFunction(bu,fname2,target_nuclideID,end_nuc);    
    };

  }; // loop-end for nuclide

};

void Burner::SensitivityCalculationDecayHeat(Burnup &bu,bool beta,bool gamma,bool alpha,string snsname)
{
  real ev_to_j=1.60219e-19;
  sensitivity=true;

  // +++ Forward calculation
  Calculation(bu,true);

  // (Sensitivity Calculation)
  //AveragingForwardNumberDensity();
  IntegratingForwardNumberDensity();

  cout<<"# Sensitivity calculation for decay heat\n";
  int sub_step=sub_step_list[burn_step-1];

  // (Set-up for final condition)
  real tmp=0.;
  real tmp1;
  real end_decayheat=0.;//[W/cm3]
  vector<int> nuc;
  vector<vector<real> > energy_sens;
  nuc.resize(nucn);
  energy_sens.resize(nucn);
  vector<bool> decay;
  decay.push_back(beta);
  decay.push_back(gamma);
  decay.push_back(alpha);

  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    tmp1=0.;
    real decay_const=bu.GetBurnupChain().GetDecayConstant(id);
    real dens=med[0].GetNuclideInTurn(i).GetDensity();
    nuc[i]=id;
    energy_sens[i].resize(3);
    for(int j=0;j<3;j++){
      if(decay[j]==true){
	real energy=bu.GetBurnupChain().GetDecayEnergy(id,j);
	tmp1=tmp1+energy;
	energy_sens[i][j]=energy*decay_const*dens;
      }else if(decay[j]==false){
	energy_sens[i][j]=0.;
      };
    };
    tmp=tmp1*decay_const*ev_to_j*1e+24;
    end_decayheat=end_decayheat+tmp*dens;
    adj_nuc_e[burn_step-1][sub_step-1].put_data(i,tmp);
  };

  for(int i=0;i<nucn;i++){
    for(int j=0;j<3;j++){
      real tmp=energy_sens[i][j];
      energy_sens[i][j]=tmp/(end_decayheat/(ev_to_j*1e+24));
    };
  };

  AdjointBurnupCalculation(bu);
  SensitivityPrintingForDecayHeat(bu,end_decayheat,nuc,energy_sens,snsname);

};

void Burner::WriteFileNumberDensity(string filename)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file in Burner::WriteFileNumberDensity.\n";
    exit(1);
  };
  fout<<"   "<<burn_step<<"\n";
  fout<<"   "<<nucn<<"\n";
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<burn_step;i++){
    fout<<burn_time_gwd[i]<<"\n";
  };
  for(int i=0;i<nucn;i++){
    fout<<med[0].GetNuclideInTurn(i).GetMatnum()<<"\n";
    for(int j=0;j<burn_step;j++){
      int ss=sub_step_list[j];
      fout<<ss<<"\n";
      for(int k=0;k<ss;k++){
        fout<<(fwd_nuc[j][k].get_dat(i)+fwd_nuc[j][k+1].get_dat(i))*0.5<<"\n";
      };
    };
  };
  fout.close();
};

void Burner::WriteFileAdjointNumberDensity(string filename,int target_nuclideID)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file in Burner::WriteFileAdjointNumberDensity.\n";
    exit(1);
  };

  fout<<"   "<<target_nuclideID<<"\n";
  fout<<"   "<<burn_step<<"\n";
  fout<<"   "<<nucn<<"\n";
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<nucn;i++){
    fout<<med[0].GetNuclideInTurn(i).GetMatnum()<<"\n";
    for(int j=0;j<burn_step;j++){
      int ss=sub_step_list[j];
      fout<<ss<<"\n";
      for(int k=0;k<ss;k++){
        fout<<adj_nuc[j][k].get_dat(i)<<"\n";
      };
    };
  };

  fout.close();
};

void Burner::WriteFileContributionFunction(Burnup &bu,string filename,int target_nuclideID,real end_nuc)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file in Burner::WriteFileContributionFunction.\n";
    exit(1);
  };

  fout<<"   "<<target_nuclideID<<"\n";
  fout<<"   "<<burn_step<<"\n";
  fout<<"   "<<nucn<<"\n";
  fout.setf(ios::scientific);
  fout.precision(5);

  for(int i=0;i<nucn;i++){
    fout<<"    "<<med[0].GetNuclideInTurn(i).GetMatnum()<<"\n";
    for(int j=burn_step-1;j>=0;j--){
      int ss=sub_step_list[j];
      fout<<ss<<"\n";
      for(int k=ss-1;k>=0;k--){
        real pf=power_factor[j][k];
	real cf=adj_nuc[j][k].get_dat(i)*fwd_nuc[j][k].get_dat(i)/end_nuc;
	if(nuclide_info[i]!=0){ // no cross section data
	  real tmp=0.;
	  for(int g=0;g<group;g++){
	    // [flux term]
	    for(int m=0;m<mesh_fuel;m++){
	      tmp+=volflx_mesh[j][m].get_dat(g)*pf*mic_sigc[j][i].get_dat(g)*gpt_flx[j][m].get_dat(g);
	      if(nuclide_info[i]==1){ // fissile
		real src=0.;
		for(int gg=0;gg<group;gg++){
		  src+=volflx_mesh[j][m].get_dat(gg)*pf*mic_sigf[j][i].get_dat(gg)
	            *med[0].GetNuclideInTurn(i).GetMicxs().GetData1d(nu).get_dat(gg);
		};
		tmp-=macxs[j].GetData1d(chi).get_dat(g)*src/keff[j]*gpt_flx[j][m].get_dat(g);
	      };
	    };
	    // [power normalization term]
	    if(nuclide_info[i]==1){
	      int iid=med[0].GetNuclideInTurn(i).GetMatnum();
	      real flx=fuel_flux[j].get_dat(g)*pf*fuel_vol;
	      tmp-=pow_adj[j][k]*flx*mic_sigf[j][i].get_dat(g)*bu.GetReactionEnergyData().GetFissionEnergy(iid);
	    };
	  };
	  cf-=tmp*fwd_nuc[j][k].get_dat(i)/end_nuc;
	}; // 
	fout<<"  "<<cf<<"\n"; // contribution function
      };
    };
  };
};

void Burner::SensitivityCalculationKeffEOC(Burnup &bu,bool vid,string prefilename)
{
  if(prefilename==""){
    prefilename="sns.k_EOC";
  };
  string filename_dir=prefilename+"_dir";
  string filename_indir=prefilename+"_indir";

  sensitivity=true;
  // +++ Forward calculation
  Calculation(bu,true);

  med[0].CalMacroFromMicro();

  Medium med2_copy=med[2];
  if(vid){
    for(int i=0;i<med[2].GetNucnum();i++){
      int mat=med[2].GetNuclideInTurn(i).GetMatnum();
      real den=med[2].GetNuclideInTurn(i).GetDensity();
      if(mat==10010||mat==80160)den*=0.9;
      med2_copy.GetNuclideInTurn(i).PutDensity(den);
    };
    med2_copy.CalMacroFromMicro();
  };

  // +++ EOC calculation
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PJISystem lata(group,3);
  lata.PutTrajectorySet(&sys_f);
  lata.AddMedium(med[0]);
  lata.AddMedium(med[1]);
  lata.AddMedium(med2_copy);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  real k_adj=lata.CalIgenPij();

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med2_copy);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real k_fwd=lat.CalIgenPij();

  // (Direct term)
  // (Sensitivity only for HM is calculated)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      //real den=med[0].GetNuclideInTurn(i).GetDensity();
      //if(id>=900000&&den>0.){
        //if((id>=900000||(id>=620000&&id<630000))&&den>0.){
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sens=lata.CalSensitivityNew(&lat,k_fwd,nucnum,nucid);
  delete [] nucid;
  sens.PutName("burner_pincell","k_EOC","unknown");
  if(vid){
    sens.WriteFile("./","sns.k_vid_EOC_dir");
  }else{
    //sens.WriteFile("./","sns.k_EOC_dir");
    sens.WriteFile("./",filename_dir);
  };

  int sub_step=sub_step_list[burn_step-1];
  // (initialization for final condition)
  for(int i=0;i<nucn;i++){
    adj_nuc_e[burn_step-1][sub_step-1].put_data(i,0.);
  };

  for(int i=0;i<nucn;i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      real org=med[0].GetNuclideInTurn(i).GetDensity();
      if(org>1e-20){
      lat.GetMedium(0).GetNuclideInTurn(i).PutDensity(org*1.01);
      lat.GetMedium(0).CalMacroFromMicro();
      real dk=lata.CalReactivity(&lat,k_adj,k_fwd,false)*k_adj*k_fwd;
      lat.GetMedium(0)=med[0];
      //int mmt=med[0].GetNuclideInTurn(i).GetMatnum();
      //cout<<mmt<<" "<<midt.Name(mmt)<<" "<<dk*100.<<"\n";
      real src=0.;
      if(org!=0.)src=dk/(org*0.01);
      adj_nuc_e[burn_step-1][sub_step-1].put_data(i,src);
      };
    };
  };
  //exit(0);

  cout<<"# Sensitivity calculation for Keff EOC\n";
  AdjointBurnupCalculation(bu);
  AveragingForwardNumberDensity();

  SensitivityData sens2;
  if(vid){
    SensitivityPrinting(bu,k_fwd,"burner_pincell","k_EOC","unknown","sns.k_vid_EOC_indir");
    sens2.ReadFile("./","sns.k_vid_EOC_indir");
    sens.AddSensitivityData(sens2);
    sens.WriteFile("./","sns.k_vid_EOC");
  }else{
    SensitivityPrinting(bu,k_fwd,"burner_pincell","k_EOC","unknown",filename_indir);
    //sens2.ReadFile("./","sns.k_EOC_indir");
    sens2.ReadFile("./",filename_indir);
    sens.AddSensitivityData(sens2);
    //sens.WriteFile("./","sns.k_EOC");
    sens.WriteFile("./",prefilename);
  };


};

void Burner::SensitivityCalculationKeffBOC(Burnup &bu,string prefilename)
{
  if(prefilename==""){
    prefilename="sns.k_BOC";
  };

  setlinebuf(stdout);

  ConstructingFuelNuclide(bu);
  if(hm_weight_init==0.)CalHeavyMetalInitialWeight(bu);

  PreCalculation();

  //
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};

  //

  SelfShieldingCalculation(0,b);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PJISystem lata(group,3);
  lata.PutTrajectorySet(&sys_f);
  lata.AddMedium(med[0]);
  lata.AddMedium(med[1]);
  lata.AddMedium(med[2]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  real k_adj=lata.CalIgenPij();

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real k_fwd=lat.CalIgenPij();

  // (Direct term)
  // (Sensitivity only for HM is calculated)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sens=lata.CalSensitivityNew(&lat,k_fwd,nucnum,nucid);
  delete [] nucid;
  sens.PutName("burner_pincell","k_BOC","unknown");
  sens.WriteFile("./",prefilename);

};

void Burner::SensitivityCalculationRRREOC(Burnup &bu,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh)
{
  sensitivity=true;
  // +++ Forward calculation
  Calculation(bu,true);

  med[0].CalMacroFromMicro();

  // +++ EOC calculation
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PJISystem lat(group,3);
  lat.PutTrajectorySet(&sys_f);
  lat.AddMedium(med[0]);
  lat.AddMedium(med[1]);
  lat.AddMedium(med[2]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  PJISystem lata(group,3);
  lata.PutTrajectorySet(&sys_f);
  lata.AddMedium(med[0]);
  lata.AddMedium(med[1]);
  lata.AddMedium(med[2]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  lata.PutPL(0); // PJI uses PutPL method in the "CalIgen" method.

  real nume=lat.CalMacroscopicReactionRate(nume_nuc,nume_id,nume_xs,on_mesh);
  real denom=lat.CalMacroscopicReactionRate(denom_nuc,denom_id,denom_xs,on_mesh);
  real rrr=nume/denom;

  // (Direct term)
  // (Sensitivity only for HM is calculated)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    real den=med[0].GetNuclideInTurn(i).GetDensity();
    if(id>=900000&&den>0.){
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sens=lata.CalSensitivityRRR(lat,nume_nuc,nume_id,nume_xs,denom_nuc,denom_id,denom_xs,on_mesh,nucnum,nucid,keff);
  delete [] nucid;
  sens.PutName("burner_pincell","RRR_EOC","unknown");
  sens.WriteFile("./","sns.RRR_EOC_dir");

  int sub_step=sub_step_list[burn_step-1];
  // (initialization for final condition)
  for(int i=0;i<nucn;i++){
    adj_nuc_e[burn_step-1][sub_step-1].put_data(i,0.);
  };

  // +++ Initial adjoint number density calculation
  for(int i=0;i<nucn;i++){
    real org=med[0].GetNuclideInTurn(i).GetDensity();
    lat.GetMedium(0).GetNuclideInTurn(i).PutDensity(org*1.01);
    // (indirect term calculation)
    lat.GetMedium(0).CalMacroFromMicro();
    real dr=lata.CalReactivity(&lat,keff,keff,false,false)*rrr;
    // (direct term calculation)
    real nume=lat.CalMacroscopicReactionRate(nume_nuc,nume_id,nume_xs,on_mesh);
    real denom=lat.CalMacroscopicReactionRate(denom_nuc,denom_id,denom_xs,on_mesh);
    real rrr_pert=nume/denom;
    dr+=(rrr_pert-rrr);
    //
    lat.GetMedium(0)=med[0];
    //int mmt=med[0].GetNuclideInTurn(i).GetMatnum();
    //cout<<mmt<<" "<<midt.Name(mmt)<<" "<<rrr<<" "<<dr<<"\n";
    real src=0.;
    if(org!=0.)src=dr/(org*0.01);
    adj_nuc_e[burn_step-1][sub_step-1].put_data(i,src);
  };
  //exit(0);

  cout<<"# Sensitivity calculation for RRR EOC\n";
  AdjointBurnupCalculation(bu);
  AveragingForwardNumberDensity();
  SensitivityPrinting(bu,rrr,"burner_pincell","RRR_EOC","unknown","sns.RRR_EOC_indir");
};

void Burner::SetArrayForSensitivityCalculation()
{
  fwd_nuc.resize(burn_step);
  fwd_nuc_int.resize(burn_step);
  dmdf_nuc.resize(burn_step);  // (dM/dphi) * (NUC_FWD)
  macxs.resize(burn_step);
  mic_sigf.resize(burn_step);
  mic_sigc.resize(burn_step);
  mic_sign2n.resize(burn_step);
  //mic_sigs.resize(burn_step);
  volflx_mesh.resize(burn_step);
  adj_nuc.resize(burn_step);
  adj_nuc_e.resize(burn_step);
  adj_nuc_s.resize(burn_step);
  pow_adj.resize(burn_step);
  gpt_flx.resize(burn_step);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    fwd_nuc_int[i].resize(sub_step);
    dmdf_nuc[i].resize(sub_step+1);
    volflx_mesh[i].resize(mesh_fuel);
    macxs[i].Init("MacroCrossSection");
    mic_sigf[i].resize(nucn);
    mic_sigc[i].resize(nucn);
    mic_sign2n[i].resize(nucn);
    //mic_sigs[i].resize(nucn);
    adj_nuc[i].resize(sub_step);
    adj_nuc_e[i].resize(sub_step);
    adj_nuc_s[i].resize(sub_step);
    pow_adj[i].resize(sub_step+1);
    gpt_flx[i].resize(mesh_fuel);
    for(int ii=0;ii<mesh_fuel;ii++){
      gpt_flx[i][ii].put_imax(group);
    };
    for(int j=0;j<sub_step;j++){
      adj_nuc[i][j].put_imax(nucn);
      adj_nuc_e[i][j].put_imax(nucn);
      adj_nuc_s[i][j].put_imax(nucn);
      fwd_nuc_int[i][j].put_imax(nucn);
    };
    for(int j=0;j<sub_step+1;j++){
      dmdf_nuc[i][j].resize(group);
      fwd_nuc[i][j].put_imax(nucn);
    };
  };

  pij_store.resize(burn_step+1);
  for(int i=0;i<burn_step+1;i++){
    pij_store[i].resize(group);
  };
};

void Burner::AdjointBurnupCalculation(Burnup &bu)
{
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PJISystem lat_gpt(group,3);
  lat_gpt.NoPrint();
  lat_gpt.PutTrajectorySet(&sys_f);
  lat_gpt.AddMedium(med[0]);
  lat_gpt.AddMedium(med[1]);
  lat_gpt.AddMedium(med[2]);
  lat_gpt.PutRegMed(region_medium);
  lat_gpt.PutGeneralOption(opt);
  lat_gpt.PutSigmaCol(sigtr);
  lat_gpt.PutBuckling(0.);
  lat_gpt.WriteProcOff();

  GroupData1D gpt_src(group);  // positive
  GroupData1D gpt_src2(group); // negative

  real fuel_vol_inv=1./fuel_vol;

  // (sub-sub-time step division)
  int ssv=40;
  //int ssv=1;

  // (Adjoint burnup calculation)
  for(int st=burn_step-1;st>=0;st--){

    real power_density=power_density_list[st];
    int sub_step=sub_step_list[st];

    lat_gpt.GetMed(0).GetMacxs().DataCopy(macxs[st]);

    // +++ Old
    lat_gpt.PutPij();

    // +++ New
    /*
    for(int g=0;g<group;g++){
      lat_gpt.GetPij(g).copy(pij_store[st][g]);
    };
    lat_gpt.PijTransformTrue();
    */

    real keff_step=lat_gpt.CalIgenPij();

    for(int i=0;i<nucn;i++){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      bu.PutNuclideData(i,id,0.,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i]);
    };
    bu.CalTransitionMatrixFluxDependentPart();
    GroupData2D trmat_flxdep=bu.GetTrmatFlxDep();

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int n=0;n<nucn;n++){
	if(fabs(adj_nuc_e[st][j].get_dat(n))>1e-20){
	  zero_adj=false;
	  break;
	};
      };
      if(zero_adj){
	cout<<"# Error in Burner::AdjointBurnupCalculation.\n";
	cout<<"# All adjoint number density is zero.\n";
	cout<<"# Step : "<<st<<"  /  Sub-step : "<<j<<"\n";
	adj_nuc_e[st][j].show_self();
	exit(0);
      };

      // +++ Adjoint nuclide density calculation
      GroupData2D mmat1=trmat_flxdep*(total_flux[st][j]*1e-24);
      GroupData2D mmat2=trmat_flxindep+mmat1;

      mmat2.Transposition();
      GroupData1D ttt2=adj_nuc_e[st][j];

      real dts=0.;
      adj_nuc[st][j]=ttt2*(0.5/ssv);
      vector<GroupData1D> ans(ssv);
      mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
      for(int k=0;k<ssv-1;k++){
        adj_nuc[st][j]=adj_nuc[st][j]+ans[k]/ssv;
      };
      adj_nuc[st][j]=adj_nuc[st][j]+ans[ssv-1]*(0.5/ssv);
      adj_nuc_s[st][j]=ans[ssv-1];
      
      // +++ Adjoint power calculation 
      real factor=1./power_density;
      if(power_density==0.)factor=1e10;
      if(input_flux_level)factor=1./(flux_level_list[st]*fuel_vol);

      pow_adj[st][j]=adj_nuc[st][j]*(mmat1*(fwd_nuc[st][j+1]+fwd_nuc[st][j]))*0.5*delt[st][j]*factor;

      // (Jump condition by adjoint power)
      GroupData1D newnuc=adj_nuc_s[st][j];
      if(!input_flux_level){
        for(int k=0;k<nucn;k++){
          real xsf1g=xsf_1g[st][k];
          if(xsf1g>0.){
            int id=med[0].GetNuclideInTurn(k).GetMatnum();
            real tmp1=xsf1g*total_flux[st][j]*fuel_vol
  	      *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
            newnuc.add_data(k,-tmp1);
	  };
	};
      };
      if(j!=0){
        adj_nuc_e[st][j-1]=newnuc;
        //newnuc.show_self();
      }else if(st!=0){
        int sst=sub_step_list[st-1];
        adj_nuc_e[st-1][sst-1]=newnuc;
      };
    };

    // +++ Generalized adjoint flux calculation
    // (source calculation)
    for(int g=0;g<group;g++){
      real tmp=0.;
      if(!input_flux_level){
        // (power term)
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]==1){
            int id=med[0].GetNuclideInTurn(k).GetMatnum();
            for(int j=sub_step-1;j>=0;j--){
               tmp+=mic_sigf[st][k].get_dat(g)*fwd_nuc[st][j].get_dat(k)
               *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j]
               *power_factor[st][j];
 	    };
	  };
        };
      }else{
        // (total-flux term)
        for(int j=sub_step-1;j>=0;j--){
          tmp+=pow_adj[st][j]*power_factor[st][j];
        };
      };
      // (number density term)
      real tmp2=0.;
      for(int j=sub_step-1;j>=0;j--){
        tmp2+=(adj_nuc[st][j]*(dmdf_nuc[st][j][g]+dmdf_nuc[st][j+1][g]))*0.5*delt[st][j]*power_factor[st][j];
      };
      tmp2*=fuel_vol_inv;
      tmp-=tmp2;

      if(tmp>0.){
        gpt_src.put_data(g,tmp);
        gpt_src2.put_data(g,0.);
      }else{
        gpt_src.put_data(g,0.);
        gpt_src2.put_data(g,-tmp);
      };
    };

    if(total_flux[st][0]>0.){
	//if(power_density>1e-5){

      lat_gpt.PutGeneralOption(opta);
        
      // (positive source)
      lat_gpt.SetZeroScatSrc();
      for(int m=0;m<mesh_fuel;m++){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src);
      };
      lat_gpt.CalGPT(keff_step,1e-4,10);
      for(int rr=0;rr<mesh_fuel;rr++){
        gpt_flx[st][rr]=lat_gpt.GetMesh(rr).GetFlux();
      };

      // (negative source)
      lat_gpt.SetZeroScatSrc();
      for(int m=0;m<mesh_fuel;m++){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2);
      };
      lat_gpt.CalGPT(keff_step,1e-4,10);
      for(int rr=0;rr<mesh_fuel;rr++){
        gpt_flx[st][rr]=gpt_flx[st][rr]-lat_gpt.GetMesh(rr).GetFlux();
      };

      // (Jump condition for generalized adjoint flux)
      if(st!=0){
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            real xsf1g=xsf_1g[st][k];
	    // (scattering term)
            real tmp1=0.;
	    /*
	    for(int m=0;m<mesh_fuel;m++){
  	      for(int g=0;g<group;g++){
		real gpt1=gpt_flx[st][m].get_dat(g);
		real flx1=volflx_mesh[st][m].get_dat(g);
	        for(int g2=g+1;g2<group;g2++){
	  	  tmp1+=(gpt_flx[st][m].get_dat(g2)-gpt1)*mic_sigs[st][k].get_dat(g,g2)*flx1;
		};
	      };
	    };
	    */
            // (absorption term)
            real tmp2=0.;
            for(int g=0;g<group;g++){
              real xsa=mic_sigc[st][k].get_dat(g)-mic_sign2n[st][k].get_dat(g); 
              if(xsf1g>0.)xsa+=mic_sigf[st][k].get_dat(g);
              for(int m=0;m<mesh_fuel;m++){
                tmp2+=volflx_mesh[st][m].get_dat(g)*gpt_flx[st][m].get_dat(g)*xsa;
              };
            };
            // (yield term)
            real tmp3=0.;
            if(xsf1g>0.){
              for(int m=0;m<mesh_fuel;m++){
                real fsrc=0.;
                for(int g=0;g<group;g++){
                  fsrc+=volflx_mesh[st][m].get_dat(g)*mic_sigf[st][k].get_dat(g)
	               *med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(nu).get_dat(g);
  	        };
		real fsrc_adj=0.;
	        for(int g=0;g<group;g++){
	          fsrc_adj+=gpt_flx[st][m].get_dat(g)*macxs[st].GetData1d(chi).get_dat(g);
		};
		tmp3+=fsrc*fsrc_adj;
	      };
	      tmp3/=keff_step;
	    };
            int sst=sub_step_list[st-1];
            adj_nuc_e[st-1][sst-1].add_data(k,-tmp1+tmp2-tmp3);
	    //cout<<"Nuc: "<<k<<" "<<tmp1<<" "<<tmp2<<" "<<tmp3<<"\n";
	  };
	};
      };
    };

  };
};

void Burner::AveragingForwardNumberDensity()
{
  // (forward nuclide density is averaged in each time sub-step)
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      for(int k=0;k<nucn;k++){
	//fwd_nuc[i][j].put_data(k,(fwd_nuc[i][j].get_dat(k)+fwd_nuc[i][j+1].get_dat(k))*0.5);
	fwd_nuc[i][j].put_data(k,sqrt(fwd_nuc[i][j].get_dat(k)*fwd_nuc[i][j+1].get_dat(k)));
      };
    };
  };
};

void Burner::IntegratingForwardNumberDensity()
{
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      real dt=delt[i][j];
      for(int k=0;k<nucn;k++){
	real n1=fwd_nuc[i][j].get_dat(k);
	real n2=fwd_nuc_int[i][j].get_dat(k);
	real n3=fwd_nuc[i][j+1].get_dat(k);
	real tmp;
	if(n1==n3){
	  tmp=n1;
	}else if(2*n2==(n1+n3)){
	  tmp=n2;
	}else if(n1>n3&&n2>=n1){
	  tmp=n1;
	}else if(n1>n3&&n2<=n3){
	  tmp=n3;
	}else if(n1<n3&&n2>=n3){
	  tmp=n3;
	}else if(n1<n3&&n2<=n1){
	  tmp=n1;
	}else {
	  real a=-(n1-n2)*(n1-n2)/(2*n2-n1-n3);
	  real b=log((n2-n3)/(n1-n2))/(dt*0.5);
	  real c=(n2*n2-n1*n3)/(2*n2-n1-n3);
	  tmp=(a*exp(b*dt)/b+c*dt-a/b)/dt;
	};
	fwd_nuc[i][j].put_data(k,tmp);
      };
    };
  };
};

void Burner::SensitivityPrinting(Burnup &bu,real response,string sys_name,string para_name,string lib_name,string filename)
{
  SensitivityData sns;
  sns.PutName(sys_name,para_name,lib_name);
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

  real end_nuc=response;
  real en_inverse=1./end_nuc;

  GroupData1D sns1d(group);
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    int rmax=3;
    if(matnum<900000)rmax=2;
    if(nuclide_info[i]==0)rmax=0; // no cross section data
    for(int r=0;r<rmax;r++){

      enum xstype sigxx=sigc;
      int mt=102;
      int bc_channel=1;
      if(r==1){
	sigxx=sign2n;
	mt=16;
	bc_channel=2;
      }else if(r==2){
        sigxx=sigf;
        mt=18;
        bc_channel=0;
      };

      // (pre-calculation for number density term)
      vector< vector<real> > nadj_dm_nfwd(burn_step);
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        nadj_dm_nfwd[k].resize(sub_step);
        for(int l=0;l<sub_step;l++){
          real val=0.;
          val+=-fwd_nuc[k][l].get_dat(i)*adj_nuc[k][l].get_dat(i);
          int tmp=bu.GetBC().GetNdiv(matnum,bc_channel);
          for(int j=0;j<tmp;j++){
            int id2=bu.GetBC().GetNextID(matnum,bc_channel,j);
            int pos=bu.SearchNuclide(id2);
            if(pos!=-1){
              real rat=bu.GetBC().GetRatio(matnum,bc_channel,j);
              val+=adj_nuc[k][l].get_dat(pos)*rat*fwd_nuc[k][l].get_dat(i);
            };
          };
          nadj_dm_nfwd[k][l]=val*delt[k][l];
        };
      };

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
          real xs=0.;
          if(r==0)xs=mic_sigc[k][i].get_dat(j);
          if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
          if(r==2)xs=mic_sigf[k][i].get_dat(j);
          real den0=fwd_nuc[k][0].get_dat(i);
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l].get_dat(i);
            // --- Number density term
            // (dM)*(flx*1e-24) = (dsig_j*phi_j/flx)*(flx*1e-24) = dsig_j*phi_j*1e-24
            real dsig=xs*(fuel_flux[k].get_dat(j)*power_factor[k][l])*1e-24;
            sum+=dsig*nadj_dm_nfwd[k][l];
            // --- Power normalization term (fission case)
	    if(!input_flux_level){
              if(sigxx==sigf){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real flx=fuel_flux[k].get_dat(j)*power_factor[k][l]*fuel_vol;
                sum-=pow_adj[k][l]*flx*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
	  };
          // --- flux term [(n,2n) reaction is not well treated yet.]

          real nu_value=0.;
          if(sigxx==sigf)nu_value=med[0].GetNuclideInTurn(i).GetMicxs().GetData1d(nu).get_dat(j);
          for(int m=0;m<mesh_fuel;m++){
            real tmp2=volflx_mesh[k][m].get_dat(j)*den0*xs;
            // (absorption)
            sum+=tmp2*gpt_flx[k][m].get_dat(j);
            // (yield)(only for fission case)
            if(sigxx==sigf){
              tmp2*=nu_value/keff[k];
	      for(int gg=0;gg<group;gg++){
  	        sum-=tmp2*gpt_flx[k][m].get_dat(gg)*macxs[k].GetData1d(chi).get_dat(gg);
	      };
            };
	  };

        };
        sns1d.put_data(j,sum*en_inverse);
      };
      sns.PutSensitivity1D(matnum,mt,sns1d);
    };
  };

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
  //int idfisn=5;
  //int idfisorg[]={922350,922380,942380,942390,942410};
  int idfisn=21;//kawamoto
  int idfisorg[]={
    922340,922350,922360,922370,922380,
    932370,932390,
    942380,942390,942400,942410,942420,
    952410,952420,952421,952430,
    962420,962430,962440,962450,962460,
  };

  for(int ii=0;ii<idfisn;ii++){
    int idfis=idfisorg[ii];
    int pos0=bu.SearchNuclide(idfis);
    int nuct=bu.GetBC().GetNdivFission(idfis);
    for(int i=0;i<nuct;i++){
      int id=bu.GetBC().GetNextIDFission(idfis,i);
      real rat=bu.GetBC().GetRatioFission(idfis,i);
      int pos=bu.SearchNuclide(id);
      if(pos!=-1){
        real val=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
          real factor=mic_sigf[k][pos0]*fuel_flux[k]*1e-24*rat;
          for(int l=0;l<sub_step;l++){
            val+=(adj_nuc[k][l].get_dat(pos)*power_factor[k][l]*factor*fwd_nuc[k][l].get_dat(pos0))*delt[k][l];
	  };
	};
        sns.PutSensitivity0D(id,18000000+idfis,val*en_inverse);
      };
    };
  };

  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda
      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int l=0;l<sub_step;l++){
          real den=fwd_nuc[k][l].get_dat(i);
          real val=-decay_c*den*adj_nuc[k][l].get_dat(i);
          int tmp=bu.GetBC().GetNdivDecay(matnum);
          for(int j=0;j<tmp;j++){
            int id2=bu.GetBC().GetNextIDDecay(matnum,j);
            int pos=bu.SearchNuclide(id2);
            if(pos!=-1){
   	      real rat=bu.GetBC().GetRatioDecay(matnum,j);
	      val+=rat*decay_c*den*adj_nuc[k][l].get_dat(pos);
	    };
	  };
	  val*=delt[k][l];
	  sum+=val;
	};
      };
      sum*=100.;// because dT=0.01T
      sns.PutSensitivity0D(matnum,8888,sum*en_inverse);
    };
  };

  // +++ For Decay Branching Ratio +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(matnum);
      if(channel>1){
	vector<real> sns_tmp1;
	vector<real> ratio;
	sns_tmp1.resize(channel);
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(matnum,j);
	  ratio[j]=rat;
	  rat*=0.01;
	  real time=0.;
	  for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
	    for(int l=0;l<sub_step;l++){
	      real den=fwd_nuc[k][l].get_dat(i);
	      real val=0.;
	      int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	      int pos=bu.SearchNuclide(id2);
	      if(pos!=-1){
		val+=adj_nuc[k][l].get_dat(pos)*decay_c*rat*den;
	      };
	      time+=delt[k][l];
	      val*=delt[k][l];
	      sum+=val;
	    };
	  };
	  sum*=100.;// because dr=0.01r
	  sum=sum*en_inverse;
	  sns_tmp1[j]=sum;
	};

	// Make constrained sensitivity
	vector<real> sns_tmp2;
	sns_tmp2.resize(channel);
	for(int j=0;j<channel;j++){
	  real tmps=0.;
	  for(int jj=0;jj<channel;jj++){
	    tmps+=sns_tmp1[jj];
	  };
	  sns_tmp2[j]=sns_tmp1[j]-tmps*ratio[j];
	};

	for(int j=0;j<channel;j++){
	  int mt=88880+j;
	  sns.PutSensitivity0D(matnum,mt,sns_tmp2[j]);
	};
      };
    };
  };

  // +++ For Capture Branching Ratio +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    int channel=bu.GetBC().GetNdivCapture(matnum);
    if(channel>1){
      vector<real> sns_tmp1;
      vector<real> ratio;
      sns_tmp1.resize(channel);
      ratio.resize(channel);
      for(int j=0;j<channel;j++){
	real sum=0.;
	real rat=bu.GetBC().GetRatioCapture(matnum,j);
	ratio[j]=rat;
	rat*=0.01;
	real time=0.;
	for(int k=0;k<burn_step;k++){
	  int sub_step=sub_step_list[k];
          real factor=mic_sigc[k][i]*fuel_flux[k]*1e-24;
	  for(int l=0;l<sub_step;l++){
	    real den=fwd_nuc[k][l].get_dat(i);
	    real val=0.;
	    int id2=bu.GetBC().GetNextIDCapture(matnum,j);
	    int pos=bu.SearchNuclide(id2);
	    if(pos!=-1){
	      val+=adj_nuc[k][l].get_dat(pos)*power_factor[k][l]*factor*rat*den;
	    };
	    time+=delt[k][l];
	    val*=delt[k][l];
	    sum+=val;
	  };
	};
	sum*=100.;// because dr=0.01r
	sum=sum*en_inverse;
	sns_tmp1[j]=sum;
      };

      // Make constrained sensitivity
      vector<real> sns_tmp2;
      sns_tmp2.resize(channel);
      for(int j=0;j<channel;j++){
        real tmps=0.;
        for(int jj=0;jj<channel;jj++){
          tmps+=sns_tmp1[jj];
        };
        sns_tmp2[j]=sns_tmp1[j]-tmps*ratio[j];
      };

      for(int j=0;j<channel;j++){
        int mt=1020+j;
        sns.PutSensitivity0D(matnum,mt,sns_tmp2[j]);
      };
    };
  };
      
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sns.WriteFile("./",filename);
};

void Burner::SensitivityPrintingForDecayHeat(Burnup &bu,real response,vector<int> &nuc,vector<vector<real> > &energy_sens,string snsname)
{
  /*
  int target_nucnum=2;
  int target_nucid[]={561361,551360};
  int target_nucturn[]={-1,-1};
  for(int ii=0;ii<target_nucnum;ii++){
    for(int i=0;i<nucn;i++){
      if(med[0].GetNuclideInTurn(i).GetMatnum()==target_nucid[ii]){
	target_nucturn[ii]=i;
      };
    };
    if(target_nucturn[ii]==-1){
      cout<<"# Error in Burner::SensitivityPrintingForDecayHeat.\n";
      exit(0);
    };
  };

  cout<<"# ";
  for(int i=0;i<target_nucnum;i++){
    WriteOut(target_nucid[i],7);
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  for(int k=0;k<burn_step;k++){
    int sub_step=sub_step_list[k];
    real del_day=(acday[k+1]-acday[k])/sub_step;
    for(int l=0;l<sub_step;l++){
      cout.precision(10);
      cout<<acday[k]+del_day*l<<" ";
      cout.precision(4);
      for(int i=0;i<target_nucnum;i++){
	cout<<adj_nuc[k][l].get_dat(target_nucturn[i])<<" ";
      };
      cout<<"\n";
    };
  };
  exit(0);
  */

  string sys_name="burner_pincell";
  string para_name="decayheat";
  string lib_name="unknown";
  string addname=snsname;

  SensitivityData sns;
  sns.PutName(sys_name,para_name,lib_name);
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

  real end_nuc=response;
  real en_inverse=1./end_nuc;

  GroupData1D sns1d(group);
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    int rmax=3;
    if(matnum<900000)rmax=2;
    if(nuclide_info[i]==0)rmax=0; // no cross section data
    for(int r=0;r<rmax;r++){
      enum xstype sigxx=sigc;
      int mt=102;
      int bc_channel=1;
      if(r==1){
	sigxx=sign2n;
	mt=16;
	bc_channel=2;
      }else if(r==2){
        sigxx=sigf;
        mt=18;
        bc_channel=0;
      };

      // (pre-calculation for number density term)
      vector< vector<real> > nadj_dm_nfwd(burn_step);
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        nadj_dm_nfwd[k].resize(sub_step);
        for(int l=0;l<sub_step;l++){
          real val=0.;
          val+=-fwd_nuc[k][l].get_dat(i)*adj_nuc[k][l].get_dat(i);
          int tmp=bu.GetBC().GetNdiv(matnum,bc_channel);
          for(int j=0;j<tmp;j++){
            int id2=bu.GetBC().GetNextID(matnum,bc_channel,j);
            int pos=bu.SearchNuclide(id2);
            if(pos!=-1){
              real rat=bu.GetBC().GetRatio(matnum,bc_channel,j);
              val+=adj_nuc[k][l].get_dat(pos)*rat*fwd_nuc[k][l].get_dat(i);
            };
          };
          nadj_dm_nfwd[k][l]=val*delt[k][l];
        };
      };

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
          real xs=0.;
          if(r==0)xs=mic_sigc[k][i].get_dat(j);
          if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
          if(r==2)xs=mic_sigf[k][i].get_dat(j);
          real den0=fwd_nuc[k][0].get_dat(i);
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l].get_dat(i);
            // --- Number density term
            // (dM)*(flx*1e-24) = (dsig_j*phi_j/flx)*(flx*1e-24) = dsig_j*phi_j*1e-24
            real dsig=xs*(fuel_flux[k].get_dat(j)*power_factor[k][l])*1e-24;
            sum+=dsig*nadj_dm_nfwd[k][l];
            // --- Power normalization term (fission case)
	    if(!input_flux_level){
              if(sigxx==sigf){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real flx=fuel_flux[k].get_dat(j)*power_factor[k][l]*fuel_vol;
                sum-=pow_adj[k][l]*flx*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
	      };
            };
	  };
          // --- flux term [(n,2n) reaction is not well treated yet.]
          real nu_value=0.;
	  if(sigxx==sigf)nu_value=med[0].GetNuclideInTurn(i).GetMicxs().GetData1d(nu).get_dat(j);
          for(int m=0;m<mesh_fuel;m++){
            real tmp2=volflx_mesh[k][m].get_dat(j)*den0*xs;
            // (absorption)
            sum+=tmp2*gpt_flx[k][m].get_dat(j);
            // (yield)(only for fission case)
            if(sigxx==sigf){
              tmp2*=nu_value/keff[k];
              for(int gg=0;gg<group;gg++){
  	        sum-=tmp2*gpt_flx[k][m].get_dat(gg)*macxs[k].GetData1d(chi).get_dat(gg);
	      };
	    };
	  };
        };
        sns1d.put_data(j,sum*en_inverse);
      };
      sns.PutSensitivity1D(matnum,mt,sns1d);
    };
  };

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
  //int idfisn=5;
  //int idfisorg[]={922350,922380,942380,942390,942410};
  //int idfisn=21;//kawamoto
  int idfisn=1+7+2+5+4+5;
  int idfisorg[]={
    902320,
    922320,922330,922340,922350,922360,922370,922380,
    932370,932390,
    942380,942390,942400,942410,942420,
    952410,952420,952421,952430,
    962420,962430,962440,962450,962460,
  };

  for(int ii=0;ii<idfisn;ii++){
    int idfis=idfisorg[ii];
    int pos0=bu.SearchNuclide(idfis);
    int nuct=bu.GetBC().GetNdivFission(idfis);
    for(int i=0;i<nuct;i++){
      int id=bu.GetBC().GetNextIDFission(idfis,i);
      real rat=bu.GetBC().GetRatioFission(idfis,i);
      int pos=bu.SearchNuclide(id);
      if(pos!=-1){
        real val=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
          real factor=mic_sigf[k][pos0]*fuel_flux[k]*1e-24*rat;
          for(int l=0;l<sub_step;l++){
            val+=(adj_nuc[k][l].get_dat(pos)*power_factor[k][l]*factor*fwd_nuc[k][l].get_dat(pos0))*delt[k][l];
	  };
	};
        sns.PutSensitivity0D(id,18000000+idfis,val*en_inverse);
      };
    };
  };

  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda
      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int l=0;l<sub_step;l++){
          real den=fwd_nuc[k][l].get_dat(i);
          real val=-decay_c*den*adj_nuc[k][l].get_dat(i);
          int tmp=bu.GetBC().GetNdivDecay(matnum);
          for(int j=0;j<tmp;j++){
            int id2=bu.GetBC().GetNextIDDecay(matnum,j);
            int pos=bu.SearchNuclide(id2);
            if(pos!=-1){
   	      real rat=bu.GetBC().GetRatioDecay(matnum,j);
	      val+=rat*decay_c*den*adj_nuc[k][l].get_dat(pos);
	    };
	  };
	  val*=delt[k][l];
	  sum+=val;
	};
      };
      sum*=100.;// because dT=0.01T
      sum=sum*en_inverse;

      //+++Add direct sensitivity term+++

      real sns_direct=0.;
      for(int j=0;j<nuc.size();j++){
	if(nuc[j]==matnum){
	  real tmpsum=0.;
	  for(int jj=0;jj<3;jj++){
	    tmpsum+=energy_sens[j][jj];
	  };
	  sns_direct=tmpsum;
	};
      };
      sum=sum-sns_direct;

      //+++++++++++++++++++++++++++++++++

      sns.PutSensitivity0D(matnum,8888,sum);
    };
  };

  // +++ For Branching Ratio +++++++++++++++++++++++++++++++++++++++++++

  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(matnum);
      if(channel>1){
	vector<real> sns_tmp1;
	vector<real> ratio;
	sns_tmp1.resize(channel);
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(matnum,j);
	  ratio[j]=rat;
	  rat*=0.01;
	  real time=0.;
	  for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
	    for(int l=0;l<sub_step;l++){
	      real den=fwd_nuc[k][l].get_dat(i);
	      real val=0.;
	      int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	      int pos=bu.SearchNuclide(id2);
	      if(pos!=-1){
		val+=adj_nuc[k][l].get_dat(pos)*decay_c*rat*den;
	      };
	      time+=delt[k][l];
	      val*=delt[k][l];
	      sum+=val;
	    };
	  };
	  sum*=100.;// because dr=0.01r
	  sum=sum*en_inverse;
	  sns_tmp1[j]=sum;
	};

	//Make constrained sensitivity
	vector<real> sns_tmp2;
	sns_tmp2.resize(channel);
	for(int j=0;j<channel;j++){
	  real tmps=0.;
	  for(int jj=0;jj<channel;jj++){
	    tmps+=sns_tmp1[jj];
	  };
	  sns_tmp2[j]=sns_tmp1[j]-tmps*ratio[j];
	};
	//

	for(int j=0;j<channel;j++){
	  int mt=88880+j;
	  sns.PutSensitivity0D(matnum,mt,sns_tmp2[j]);
	};
      };
    };
  };
      
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ For Decay Energy +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nuc.size();i++){
    for(int j=0;j<3;j++){
      int mt=99990+j;
      sns.PutSensitivity0D(nuc[i],mt,energy_sens[i][j]);
    };
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  string filename="sns."+addname;
  sns.WriteFile("./",filename);
};

void Burner::SensitivityRead(string filename,string prtname)
{
  int prtid=midt.ID(prtname);

  bool nucname=true;

  int maxn=10;
  int maxfis=5;

  real sens_low=1e-2;

  vector<int> target_id(maxn);

  vector< vector<int> > xssns_mat(maxn);
  vector< vector<int> > xssns_mt(maxn);
  vector< vector<GroupData1D> > xssns(maxn);

  vector< vector<int> > yldsns_fis(maxn);
  vector< vector< vector<int> > > yldsns_mat(maxn);
  vector< vector< vector<real> > > yldsns(maxn);
  for(int i=0;i<maxn;i++){
    yldsns_mat[i].resize(maxfis);
    yldsns[i].resize(maxfis);
  };

  vector< vector<int> > hlsns_mat(maxn);
  vector< vector<real> > hlsns(maxn);
  // 

  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(0);
  };

  int target_nuc=0;

  int itmp=0;
  real tmp;
  while(itmp!=-1){
    if(target_nuc==10){
      cout<<"Error !\n";
      cout<<"Please increase maxn\n";
      exit(0);
    };

    fin>>itmp;
    if(itmp==-1)break;
    target_id[target_nuc]=itmp;

    int idflag=0;
    while(idflag!=-1){
      fin>>idflag;
      if(idflag==3){ // XS sensitivity read
        int grp;
        fin>>grp;
        int mat;
        fin>>mat;
        int mt;
        fin>>mt;
        xssns_mat[target_nuc].push_back(mat);
        xssns_mt[target_nuc].push_back(mt);
        GroupData1D dat(grp);
        for(int g=0;g<grp;g++){
	  fin>>tmp;
	  dat.put_data(g,tmp);
        };
        xssns[target_nuc].push_back(dat);
      }else if(idflag==18){
	int sz=yldsns_fis[target_nuc].size();
	if(sz==maxfis){
	  cout<<"Error !\n";
	  cout<<"Please increase maxfis.\n";
	  exit(0);
	};
	fin>>itmp;
	yldsns_fis[target_nuc].push_back(itmp);
	int inuc;
	fin>>inuc;
	for(int i=0;i<inuc;i++){
	  fin>>itmp;
	  yldsns_mat[target_nuc][sz].push_back(itmp);
	  fin>>tmp;
	  yldsns[target_nuc][sz].push_back(tmp);
	};
      }else if(idflag==88){
	int idf=0;
	while(idf!=-1){
  	  fin>>idf;
	  if(idf!=-1){
  	    fin>>tmp;
	    hlsns_mat[target_nuc].push_back(idf);
	    hlsns[target_nuc].push_back(tmp);
	  };
	};
      };
    };
    target_nuc++;
  };

  // ++++++++++++++

  int tid=-1;
  for(int i=0;i<target_nuc;i++){
    if(prtid==target_id[i])tid=i;
  };
  if(tid==-1){
    cout<<"+ No data "<<prtid<<"\n";
  };

  // +++ XS sensitivity print +++
  cout<<"\n+++ 1-group cross section sensitivity +++\n\n";
  int xs_size=xssns_mat[tid].size();
  for(int i=0;i<xs_size;i++){
    int mat=xssns_mat[tid][i];
    int mt=xssns_mt[tid][i];
    real sns1g=xssns[tid][i].get_sum();
    if(fabs(sns1g)>sens_low){
      cout.setf(ios::scientific);
      cout.precision(5);
      if(nucname){
        string tmp=midt.Name(mat);
        cout<<tmp<<" "<<mt<<" "<<sns1g<<"\n";
      }else{
        cout<<mat<<" "<<mt<<" "<<sns1g<<"\n";
      };
    };
  };

  // +++ Fission yield print +++
  cout<<"\n+++ fission yield sensitivity +++\n";
  cout<<"   (fis-nuc / FP-nuc / sens )\n\n";
  //int fy_size=yldsns_fis[tid].size()+1;
  int fy_size=yldsns_fis[tid].size();
  real sum=0.;
  for(int i=0;i<fy_size;i++){
    int fismat=yldsns_fis[tid][i];
    //int nuc=yldsns_mat[tid][i].size()+1;
    int nuc=yldsns_mat[tid][i].size();
    for(int j=0;j<nuc;j++){
      sum+=yldsns[tid][i][j];
      if(yldsns[tid][i][j]>sens_low){
	if(nucname){
          string tmp=midt.Name(fismat);
          string tmp2=midt.Name(yldsns_mat[tid][i][j]);
	  cout<<tmp<<" "<<tmp2<<" "<<yldsns[tid][i][j]<<"\n";
	}else{
  	  cout<<fismat<<" "<<yldsns_mat[tid][i][j]<<" "<<yldsns[tid][i][j]<<"\n";
	};
      };
    };
  };
  cout<<"#  (Sum) : "<<sum<<"\n\n";

  // +++ Half-life 
  cout<<"\n+++ half-life sensitivity +++\n\n";
  int hl_size=hlsns_mat[tid].size();
  for(int j=0;j<hl_size;j++){
    int mat=hlsns_mat[tid][j];
    if(hlsns[tid][j]>sens_low){
      if(nucname){
        string tmp=midt.Name(mat);
        cout<<tmp<<" "<<hlsns[tid][j]<<"\n";
      }else{
        cout<<mat<<" "<<hlsns[tid][j]<<"\n";
      };
    };
  };

};

void Burner::AddActinideDecayHeatDataToBurnupChain(string cbglibdir,string fname,Burnup &bu)
{
  cbglibdir=cbglibdir+"CBGLIB_BURN/DDfile/"+fname;

  BCGManager man;
  man.ReadDecayDataFromFile(cbglibdir);

  // Heavy isotopes than Tl are considered.
  int sz=man.GetNuclideNumber();
  for(int i=0;i<sz;i++){
    int atm=man.GetAtomicNumber(i);
    if(atm>=81){
      int id=man.GetID(i);
      if(bu.GetBurnupChain().GetNuclideChainData(id).GetID()==id){
	//cout<<atm<<" "<<id<<"\n";
        for(int j=0;j<3;j++){
  	  real ee=man.GetNuclide(i).GetDecayEnergy(j);
	  bu.GetBurnupChain().GetNuclideChainData(id).PutDecayEnergy(j,ee);
	};
      };
    };
  };
  /*
  int iill=8;
  int ialist[]={92,94,94,94,94,95,96,96};
  int izlist[]={237,238,239,240,241,241,242,244};
  int matno[8];
  for(int i=0;i<8;i++){
    matno[i]=ialist[i]*10000+izlist[i]*10;
  };
  for(int i=0;i<iill;i++){
    for(int j=0;j<3;j++){
      real ee=man.GetNuclide(ialist[i],izlist[i],0).GetDecayEnergy(j);
      bu.GetBurnupChain().GetNuclideChainData(matno[i]).PutDecayEnergy(j,ee);
    };
  };
  */
};

void Burner::GetNuclideDensity(vector<string> &nuc,vector<real> &den,int step)
{
  if(step<0||step>burn_step){
    cout<<"# Error in Burner::GetNuclideDensity.\n";
    cout<<"# Chosen burnup step is inappropriate.\n";
    cout<<"# You chose the step "<<step<<"\n";
    cout<<"# Maximum step is "<<burn_step<<"\n";
    exit(0);
  };

  nuc.clear();
  den.clear();
  nuc.resize(nucn);
  den.resize(nucn);
  for(int i=0;i<nucn;i++){
    nuc[i]=midt.Name(med[0].GetNuclideInTurn(i).GetMatnum());
    den[i]=density_data[step][i];
  };
};

/*
  real Burner::GetNuclideDensity(string nuc,int step)
  {
  if(step<0||step>burn_step){
  cout<<"# Error in Burner::GetNuclideDensity.\n";
  cout<<"# Chosen burnup step is inappropriate.\n";
  cout<<"# You chose the step "<<step<<"\n";
  cout<<"# Maximum step is "<<burn_step<<"\n";
  exit(0);
  };

  int nucid=midt.ID(nuc);

  for(int i=0;i<nucn;i++){
  int tmp=med[0].GetNuclideInTurn(i).GetMatnum();
  if(tmp==nucid)return density_data[step][i];
  };

  cout<<"# Error in Burner::GetNuclideDensity.\n";
  cout<<"# No material in fuel region.\n";
  cout<<"# You chose the nuclide "<<nuc<<"\n";
  exit(0);
  };
*/

real Burner::GetNuclideDensity(int nucid,int step)
{
  if(step<0||step>burn_step){
    cout<<"# Error in Burner::GetNuclideDensity.\n";
    cout<<"# Chosen burnup step is inappropriate.\n";
    cout<<"# You chose the step "<<step<<"\n";
    cout<<"# Maximum step is "<<burn_step<<"\n";
    exit(0);
  };

  for(int i=0;i<nucn;i++){
    int tmp=med[0].GetNuclideInTurn(i).GetMatnum();
    if(tmp==nucid)return density_data[step][i];
  };

  cout<<"# Error in Burner::GetNuclideDensity.\n";
  cout<<"# No material in fuel region.\n";
  cout<<"# You chose the nuclide : "<<nucid<<"\n";
  exit(0);
};

// +++++++++++++++++++++++++++++++++++++++++++
//
// +++++++++++++++++++++++++++++++++++++++++++

void Burner::CobaltActivationCalculation(real den_init, real flux_factor)
// (implemented in 2012/09/01, see the notebook)
{
  real half_life_co60=5.271*365*24*60*60; // (s)
  //real br_co60g=0.444; // branching ratio of Co-59 (n,g) to ground (from JEFF-3.1A)
  real br_co60g=1.; // since half-life of Co-60 meta is about 10 minutes
  real decay_e_co60=9.67734e4+2.50384e6; // (eV) from JEFF-3.1.1 decay data file
  real ev_to_j=1.60219e-19;

  real dc=0.6931471818056/half_life_co60;

  vector<real> n_co59(burn_step+1);
  vector<real> n_co60(burn_step+1);

  vector<int> bgrp(1);
  bgrp[0]=group;

  GroupData1D sigc_co59=xslib.GetLibData(270590).GetXSData().GetData1d(sigc);

  cout<<"#\n# +++ Co-60 activation calculation +++\n#\n";
  real n59=den_init;
  real n60=0.;
  for(int i=0;i<burn_step+1;i++){
    n_co59[i]=n59;
    n_co60[i]=n60;
    if(i!=burn_step){
      real sig=sigc_co59.Cond(clad_flux[i],1,bgrp).get_dat(0)*1e-24;
      if(i==0)cout<<"#   Co-59 (n,g) cross section at initial state : "<<sig*1e24<<"\n#\n";
      for(int j=0;j<sub_step_list[i];j++){
        real flx=total_flux_clad[i][j]*flux_factor;
        real dt=delt[i][j];
        real coef=br_co60g*n59*sig*flx/(dc-sig*flx);
        n59=n59*exp(-sig*flx*dt);
        n60=coef*(exp(-sig*flx*dt)-exp(-dc*dt))+n60*exp(-dc*dt);
      };
    };
  };

  cout<<"#\n# +++ Information of Co-60 in clad region +++\n";
  if(flux_factor!=1.)cout<<"#     (flux is multiplied by "<<flux_factor<<")\n";
  cout<<"# (day)   (GWd/t)   N.D.[/cm3]  N.D.[/cm]   Heat[W/cm] Heat[W/tHM] RAD[Bq/tHM]\n#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    cout<<n_co60[i]*1e24<<" ";
    cout<<n_co60[i]*1e24*clad_vol<<" ";
    cout<<n_co60[i]*1e24*clad_vol*dc*decay_e_co60*ev_to_j<<" ";
    cout<<n_co60[i]*1e24*clad_vol*dc*decay_e_co60*ev_to_j/(hm_weight_init*1e-6)<<" ";
    cout<<n_co60[i]*1e24*clad_vol*dc/(hm_weight_init*1e-6)<<" ";
    cout<<"\n";
  };

  Nuclide co60;
  co60.PutMatnum(270600);
  real density_org=n_co60[burn_step];
  density_org*=clad_vol/fuel_vol;
  co60.PutDensity(density_org);
  med[0].AddNuclide(co60);
};

void Burner::Mn54ActivationCalculation(real den_init,real flux_factor)
{
  if(group!=70&&group!=107){
    cout<<"# Error in Burner::Mn54ActivationCalculation.\n";
    cout<<"# 70- and 107-group calculations have been implemented.\n";
    cout<<"# You are calculating with "<<group<<"-group.\n";
    exit(0);
  };

  if(den_init==0.){
    if(!med[1].ExistNuclide(260000)){
      cout<<"# Mn-54 activation calculation is not performed\n";
      cout<<"# because Fe-54 is not included in clad regin.\n";
      return;
    }else{
      den_init=med[1].GetNuclide(260000).GetDensity()*0.058;// (5.8% in natural isotope)
    };
  };

  real half_life_mn54=312.05*24*60*60; // (s)
  real decay_e_mn54=8.36e5+4.03e3; // (eV) from JEFF-3.1.1 decay data file
  real ev_to_j=1.60219e-19;

  real dc=0.6931471818056/half_life_mn54;

  vector<real> n_fe54(burn_step+1);
  vector<real> n_mn54(burn_step+1);

  vector<int> bgrp(1);
  bgrp[0]=group;

  GroupData1D sigp_fe54(group);
  sigp_fe54.set_zero();
  real datinp[]={
    5.003E-01, 5.020E-01, 4.483E-01, 2.967E-01, 1.831E-01,
    7.703E-02, 1.964E-02, 2.661E-03, 3.949E-04, 1.126E-05
  }; // JENDL-4.0 multigroup data
  for(int i=0;i<10;i++){
    sigp_fe54.put_data(i,datinp[i]);
  };

  cout<<"#\n# +++ Mn-54 activation calculation +++\n#\n";
  real n0=den_init;
  real n1=0.;
  for(int i=0;i<burn_step+1;i++){
    n_fe54[i]=n0;
    n_mn54[i]=n1;
    if(i!=burn_step){
      real sig=sigp_fe54.Cond(clad_flux[i],1,bgrp).get_dat(0)*1e-24;
      if(i==0)cout<<"#   Fe-54 (n,p) cross section at initial state : "<<sig*1e24<<"\n#\n";
      for(int j=0;j<sub_step_list[i];j++){
        real flx=total_flux_clad[i][j]*flux_factor;
        real dt=delt[i][j];
        real coef=n0*sig*flx/(dc-sig*flx);
        n0=n0*exp(-sig*flx*dt);
        n1=coef*(exp(-sig*flx*dt)-exp(-dc*dt))+n1*exp(-dc*dt);
      };
    };
  };

  cout<<"#\n# +++ Information of Mn-54 in clad region +++\n";
  if(flux_factor!=1.)cout<<"#     (flux is multiplied by "<<flux_factor<<")\n";
  cout<<"# (day)   (GWd/t)   N.D.[/cm3]  N.D.[/cm]   Heat[W/cm] Heat[W/tHM] RAD[Bq/tHM]\n#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    cout<<n_mn54[i]*1e24<<" ";
    cout<<n_mn54[i]*1e24*clad_vol<<" ";
    cout<<n_mn54[i]*1e24*clad_vol*dc*decay_e_mn54*ev_to_j<<" ";
    cout<<n_mn54[i]*1e24*clad_vol*dc*decay_e_mn54*ev_to_j/(hm_weight_init*1e-6)<<" ";
    cout<<n_mn54[i]*1e24*clad_vol*dc/(hm_weight_init*1e-6)<<" ";
    cout<<"\n";
  };
};

void Burner::Reprocessing(real tc, real ru_rh, real fp, real ac)
{
  int nucnum=med[0].GetNucnum();
  for(int i=0;i<nucnum;i++){
    int mat=med[0].GetNuclideInTurn(i).GetMatnum();
    real den=med[0].GetNuclideInTurn(i).GetDensity();
    real ratio=fp;
    if(mat>=900000){
      ratio=ac;
    }else if(mat>=430000&&mat<=439999){
      ratio=tc;
    }else if(mat>=440000&&mat<=459999){
      ratio=ru_rh;
    };
    med[0].GetNuclideInTurn(i).PutDensity(den*ratio);
  };
};

void Burner::NumberDensityReset(int mat_st,int mat_ed,real den)
{
  if(mat_st>mat_ed){
    cout<<"# Error in Burner::NumberDensityReset.\n";
    cout<<"# Please check start/end points of material ID.\n";
    exit(0);
  };

  int nucnum=med[0].GetNucnum();
  for(int i=0;i<nucnum;i++){
    int mat=med[0].GetNuclideInTurn(i).GetMatnum();
    if(mat>=mat_st&&mat<=mat_ed){
      real newden=den;
      if(den<0.){
	newden=med[0].GetNuclideInTurn(i).GetDensity()*den*-1;
      };
      med[0].GetNuclideInTurn(i).PutDensity(newden);
    };
  };
};

void Burner::NumberDensityReset(int mat,real den)
{
  if(!med[0].ExistNuclide(mat)){
    cout<<"# Error in Burner::NumberDensityReset.\n";
    cout<<"# Nuclide "<<mat<<" does NOT exist.\n";
    exit(0);
  };
  med[0].GetNuclide(mat).PutDensity(den);
};

// +++++++++++++++++++++++++++++++++++++++++++
//  LIMITED BURNUP CALCULATION
// +++++++++++++++++++++++++++++++++++++++++++

void Burner::LimitBurnupOn(real kl, int batchinp)
{
  limit_burnup=true;
  k_limit=kl;
  batch=batchinp;
};

real Burner::GetEOCBurnup(){
  if(!limit_burnup){
    cout<<"# Error in Burner::GetEOCBurnup\n";
    cout<<"# This method is effective only for limited burnup calculation.\n";
    exit(0);
  };
  return bumax;
};

real Burner::GetAtomWiseWeightAtBOC(Burnup &bu,int atm)
{
  real avo=0.60221367;

  real sum=0.;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=atm*100&&id<(atm+1)*100){
      real dd=density_data[0][i];
      real aw=bu.GetAtomicWeight(id);
      real mol=(dd*fuel_vol)/avo;
      real wt=aw*mol*1e-3;
      sum+=wt/(hm_weight_init*1e-6);
    };
  };
  return sum;
};

real Burner::GetAtomWiseWeightAtEOC(Burnup &bu,int atm)
{
  real avo=0.60221367;

  real bu0=acburn[burn_step-1];
  real bu1=acburn[burn_step];

  real sum=0.;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=atm*100&&id<(atm+1)*100){
      real dd0=density_data[burn_step-1][i];
      real dd1=density_data[burn_step][i];
      real dd=dd0+(dd1-dd0)/(bu1-bu0)*(bumax-bu0);
      real aw=bu.GetAtomicWeight(id);
      real mol=(dd*fuel_vol)/avo;
      real wt=aw*mol*1e-3;
      sum+=wt/(hm_weight_init*1e-6);
    };
  };
  return sum;
};

// +++++++++++++++++++++++++++++++++++++++++++
//  Lumped FP creation
// +++++++++++++++++++++++++++++++++++++++++++

void Burner::CreatePseudoFP(Burnup &bu,string fissle_name,int nucnum,string *nucname,int outid,string outname)
{
  int fisid=midt.ID(fissle_name);
  
  vector<real> sc(group,0.);
  int fpnum=bu.GetBurnupChain().GetNdivFission(fisid);

  vector<int> nucid(nucnum);
  for(int i=0;i<nucnum;i++){
    nucid[i]=midt.ID(nucname[i]);
  };

  real sum=0.;
  for(int i=0;i<fpnum;i++){
    int fpid=bu.GetBurnupChain().GetNextIDFission(fisid,i);
    real fprat=bu.GetBurnupChain().GetRatioFission(fisid,i);
    bool same=false;
    for(int j=0;j<nucnum;j++){
      if(fpid==nucid[j])same=true;
    };
    if(!same){
      real dc=bu.GetBurnupChain().GetDecayConstant(fpid);
      real hl=0.;
      if(dc>0.)hl=0.693/dc;
      if(xslib.GetLibData(fpid).GetMat()!=-1){
	//if(xslib.GetLibData(fpid).GetMat()!=-1){
	//if(xslib.GetLibData(fpid).GetMat()!=-1&&(hl==0.||hl>60*60*24*10)){
        for(int g=0;g<group;g++){
  	  real sigval=xslib.GetLibData(fpid).GetXSData().GetData1d(sigc).get_dat(g);
	  sc[g]+=sigval*fprat;
        };
        sum+=fprat;
      };
    };
  };

  sum=1./sum;
  for(int g=0;g<group;g++){
    sc[g]*=sum;
  };

  LibData pfp;
  pfp.PutGroup(group);
  pfp.PutMat(outid);
  for(int g=0;g<group;g++){
    pfp.GetXSData().GetData1d(sigc).put_data(g,sc[g]);
    pfp.GetXSData().GetData1d(sigt).put_data(g,sc[g]);
    //cout<<xslib.GetEnband().get_dat(g)<<" "<<sc[g]<<"\n";
  };

  pfp.GetFtable().PutNomtft(0);
  pfp.GetFtable().PutNsig0(0);
  pfp.GetFtable().PutMaxtemp(0);
  pfp.GetFtable().PutMaxnr(0);
  pfp.WriteFile("./",outname);

};

void Burner::LFP(Burnup &bu,string fissle_name,int nucnum,string *nucname)
{
  int fisid=midt.ID(fissle_name);
  
  vector<real> sc(group,0.);
  int fpnum=bu.GetBurnupChain().GetNdivFission(fisid);

  vector<int> nucid(nucnum);
  for(int i=0;i<nucnum;i++){
    nucid[i]=midt.ID(nucname[i]);
  };

  real sum=0.;
  for(int i=0;i<fpnum;i++){
    int fpid=bu.GetBurnupChain().GetNextIDFission(fisid,i);
    real fprat=bu.GetBurnupChain().GetRatioFission(fisid,i);
    bool same=false;
    for(int j=0;j<nucnum;j++){
      if(fpid==nucid[j])same=true;
    };
    if(!same){
      if(xslib.GetLibData(fpid).GetMat()!=-1){
        for(int g=0;g<group;g++){
  	  real sigval=xslib.GetLibData(fpid).GetXSData().GetData1d(sigc).get_dat(g);
	  sc[g]+=sigval*fprat;
        };
      };
      sum+=fprat;
    };
  };

  sum=1./sum;
  for(int g=0;g<group;g++){
    sc[g]*=sum;
  };

  LibData pfp;
  pfp.PutGroup(group);
  pfp.PutMat(9999);
  for(int g=0;g<group;g++){
    pfp.GetXSData().GetData1d(sigc).put_data(g,sc[g]);
    pfp.GetXSData().GetData1d(sigt).put_data(g,sc[g]);
  };
  pfp.GetFtable().PutNomtft(0);
  pfp.GetFtable().PutNsig0(0);
  pfp.GetFtable().PutMaxtemp(0);
  pfp.GetFtable().PutMaxnr(0);
  pfp.WriteFile("./","pfp");

  /*
  // (judging from relative error in all/some groups)
  real errmax=1e10;
  int eid=-1;
  for(int i=0;i<fpnum;i++){
  int fpid=bu.GetBurnupChain().GetNextIDFission(fisid,i);
  if(xslib.GetLibData(fpid).GetMat()!=-1){
  real errsum=0.;
  real errmax2=0.;
  for(int g=0;g<group;g++){
  //for(int g=59;g<group;g++){ // (below 3.9279eV)
  real err=fabs(xslib.GetLibData(fpid).GetXSData().GetData1d(sigc).get_dat(g)/sc[g]-1.);
  errsum=err*err;
  if(errsum>errmax2)errmax2=errsum;
  };
  if(errsum<errmax){
  errmax=errsum;
  eid=fpid;
  };

  };
  };
  */
  
  // (judging from thermal capture)
  /*
    real errmax=1e10;
    int eid=-1;
    for(int i=0;i<fpnum;i++){
    int fpid=bu.GetBurnupChain().GetNextIDFission(fisid,i);
    if(xslib.GetLibData(fpid).GetMat()!=-1){
    int g=98;
    real err=fabs(xslib.GetLibData(fpid).GetXSData().GetData1d(sigc).get_dat(g)-sc[g]);
    if(err<errmax){
    errmax=err;
    eid=fpid;
    };
    };
    };
  */

  /*
    cout<<"# Candidate is "<<midt.Name(eid)<<"\n";
    for(int g=0;g<group;g++){
    real e=xslib.GetEnband().get_dat(g);
    cout<<e<<" "<<sc[g]<<" ";
    real xs=xslib.GetLibData(eid).GetXSData().GetData1d(sigc).get_dat(g);
    cout<<xs<<"\n";
    };
  */
};

void Burner::ShowPseudoFPCrossSection(Burnup &bu,int nucnum,string *nucname)
{
  vector<int> nucid(nucnum);
  for(int i=0;i<nucnum;i++){
    nucid[i]=midt.ID(nucname[i]);
  };

  vector< vector<real> > cont(burn_step+1);
  vector<real> contsum(burn_step+1,0.);

  cout<<"################################################\n";
  cout<<"# Cross sections of non-important FP nuclides. #\n";
  cout<<"################################################\n";

  for(int i=1;i<=burn_step;i++){

    cout<<"# Step   : "<<i<<"\n";
    cout<<"# Burnup : "<<acburn[i]<<" [GWd/t]\n";

    cont[i].resize(nucn);
    vector<real> sc(group,0.);
    real sum=0.;
    for(int j=0;j<nucn;j++){
      int id=med[0].GetNuclideInTurn(j).GetMatnum();
      bool same=false;
      for(int k=0;k<nucnum;k++){
        if(id==nucid[k])same=true;
      };
      if(!same){
	if(xslib.ExistLibData(id)&&id<900000&&id>=100000){
	  real den=density_data[i][j];
	  real tmp=den*xslib.GetLibData(id).GetXSData().GetData1d(sigc).get_dat(99);
	  cont[i][j]=tmp;
	  contsum[i]+=tmp;
	  for(int g=0;g<group;g++){
	    real sigval=xslib.GetLibData(id).GetXSData().GetData1d(sigc).get_dat(g);
	    sc[g]+=sigval*den;
	  };
	  sum+=den;
	};
      };
    };
 
    if(sum!=0)sum=1./sum;
    for(int g=0;g<group;g++){
      sc[g]*=sum;
    };

    cout<<"\n\n# [Averaged-microscopic cross section]\n";
    // cross section printing
    for(int g=0;g<group;g++){
      cout.setf(ios::scientific);
      cout.precision(5);
      real e=xslib.GetEnband().get_dat(g);
      cout<<e<<" "<<sc[g]<<"\n";
    };
    cout<<"\n\n";

  };

  // contribution printing for group 100
  cout<<" Nuc ";
  for(int i=1;i<burn_step;i++){
    cout<<acburn[i]<<" ";
  };
  cout<<"\n";

  for(int j=0;j<nucn;j++){
    int id=med[0].GetNuclideInTurn(j).GetMatnum();
    bool same=false;
    for(int k=0;k<nucnum;k++){
      if(id==nucid[k])same=true;
    };
    if(!same){
      if(xslib.ExistLibData(id)&&id<900000&&id>99999){
	bool large=false;
	for(int i=1;i<burn_step;i++){
          if(cont[i][j]/contsum[i]>0.01)large=true;
 	};
	if(large){
          cout<<midt.Name(id)<<" ";
	  for(int i=1;i<burn_step;i++){
	    cout<<cont[i][j]/contsum[i]<<" ";
	  };
    	  cout<<"\n";
	};
      };
    };
  };
};

// +++++++++++++++++++++++++++++++++++++++++++
//  OUTPUT EDITING
// +++++++++++++++++++++++++++++++++++++++++++

/*
  void Burner::ShowDecayHeat(Burnup &bu,bool detail_print,real factor)
  {
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"#\n# Time-dependent decay heat (day, W/cm , W/tHM";
  if(detail_print)cout<<", ([HM][FP])";
  cout<<")\n";
  if(factor!=1.){
  cout<<"#     factorized by "<<factor<<"\n";
  };
  cout<<"#\n";
  real ev_to_j=1.60219e-19;
  ev_to_j*=factor;
  for(int i=0;i<burn_step+1;i++){
  cout<<acday[i]<<" ";
  real en=0.;
  real en_hm=0.;
  real en_fp=0.;
  for(int j=0;j<nucn;j++){
  int id=med[0].GetNuclideInTurn(j).GetMatnum();
  real den=density_data[i][j]*1e24;
  real lambda=bu.GetDecayConstant(id);
  real e=0.;
  for(int k=0;k<3;k++){
  e+=bu.GetBurnupChain().GetDecayEnergy(id,k);
  };
  real enuc=e*(den*fuel_vol)*lambda;
  en+=enuc;
  if(id>=900000){en_hm+=enuc;}else{en_fp+=enuc;};
  };
  en*=ev_to_j; // J
  en_hm*=ev_to_j; 
  en_fp*=ev_to_j;
  cout<<en<<" ";
  if(detail_print){
  cout<<"( "<<en_hm<<" ";
  cout<<en_fp<<" ) ";
  };
  cout<<en/(hm_weight_init*1e-6)<<" "; // J/tHM
  if(detail_print){
  cout<<"( "<<en_hm/(hm_weight_init*1e-6)<<" ";
  cout<<en_fp/(hm_weight_init*1e-6)<<" ) ";
  };
  cout<<"\n";
  };

  };
*/

/*
  void Burner::ShowRadioactivity(Burnup &bu,bool detail_print,real factor)
  {
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"#\n# Time-dependent radioactivity (day, year, Bq/cm , Bq/tHM";
  if(detail_print)cout<<", ([HM][FP])";
  cout<<")\n";
  if(factor!=1.){
  cout<<"#     factorized by "<<factor<<"\n";
  };
  cout<<"#\n";
  for(int i=0;i<burn_step+1;i++){
  cout<<acday[i]<<" "<<acday[i]/365.<<" ";
  real en=0.;
  real en_hm=0.;
  real en_fp=0.;
  for(int j=0;j<nucn;j++){
  int id=med[0].GetNuclideInTurn(j).GetMatnum();
  real den=density_data[i][j]*1e24;
  real lambda=bu.GetDecayConstant(id);
  real enuc=(den*fuel_vol)*lambda;
  en+=enuc;
  if(id>=900000){en_hm+=enuc;}else{en_fp+=enuc;};
  };
  cout<<en<<" ";
  if(detail_print){
  cout<<"  "<<en_hm<<" ";
  cout<<en_fp<<"   ";
  };
  cout<<en/(hm_weight_init*1e-6)<<" "; // Bq/tHM
  if(detail_print){
  cout<<" "<<en_hm/(hm_weight_init*1e-6)<<" ";
  cout<<en_fp/(hm_weight_init*1e-6)<<"  ";
  };
  cout<<"\n";
  };

  };
*/
/*
  void Burner::ShowDecayHeatAtomWise(Burnup &bu,string nucname,real factor)
  {
  int z=midt.AtomicNumberFromName(nucname);
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"#\n# Time-dependent decay heat of "<<nucname<<" (day, W/cm , W/tHM,)\n";
  if(factor!=1.){
  cout<<"#     factorized by "<<factor<<"\n";
  };
  
  real ev_to_j=1.60219e-19*factor;
  
  for(int i=0;i<burn_step+1;i++){
  cout<<acday[i]<<" ";
  real en=0.;
  for(int j=0;j<nucn;j++){
  int id=med[0].GetNuclideInTurn(j).GetMatnum();
  int nucz=(int)id/10000;
  if(nucz==z){
  real den=density_data[i][j]*1e24;
  real lambda=bu.GetDecayConstant(id);
  real e=0.;
  for(int k=0;k<3;k++){
  e+=bu.GetBurnupChain().GetDecayEnergy(id,k);
  };
  real enuc=e*(den*fuel_vol)*lambda;
  en+=enuc;
  };
  };
  en*=ev_to_j; // J
  cout<<en<<" ";
  cout<<en/(hm_weight_init*1e-6); // J/tHM
  cout<<"\n";
  };
  };
*/

void Burner::GetPrintNuclide(int prt_nuc,string *prt_nuc_nam,vector<int> &prt_nuc_turn)
{
  for(int i=0;i<prt_nuc;i++){
    int id=midt.ID(prt_nuc_nam[i]);
    for(int j=0;j<nucn;j++){
      if(med[0].GetNuclideInTurn(j).GetMatnum()==id){
	prt_nuc_turn[i]=j;
	break;
      };
    };
    if(prt_nuc_turn[i]==-1){
      cout<<"# Error in Burner::GetPrintNuclide.\n";
      cout<<"# Nuclide "<<prt_nuc_nam[i]<<" is NOT included in the burnup chain.\n";
      exit(0);
    };
  };
};

void Burner::ShowEigenvalue()
{
  // (for conversion ratio calculation)
  int fert_num=5;
  int fert_mat[]={902320,922340,922380,942400,-912330};
  int fiss_num=4;
  int fiss_mat[]={922330,922350,942390,942410};

  cout<<"#\n# +++ Time-dependent eigenvalue +++\n";
  cout<<"#  (day)       ";
  cout<<"(GWd/t)     ";
  cout<<"(k-inf)     ";
  cout<<"(C.R.)      ";
  cout<<"(flux[/cm2/s])";
  cout<<"\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<burn_step+1;i++){
    if(i<10)cout<<" ";
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" "<<keff[i]<<" ";
    // (conversion ratio)
    real fert_cap=0.;
    for(int j=0;j<fert_num;j++){
      int mat=abs(fert_mat[j]);
      for(int k=0;k<nucn;k++){
	if(med[0].GetNuclideInTurn(k).GetMatnum()==mat){
  	  real factor=1.;
	  if(fert_mat[j]<0)factor=-1.;
          fert_cap+=density_data[i][k]*xsc_1g[i][k]*factor;
	};
      };
    };
    real fiss_abs=0.;
    for(int j=0;j<fiss_num;j++){
      int mat=abs(fiss_mat[j]);
      for(int k=0;k<nucn;k++){
	if(med[0].GetNuclideInTurn(k).GetMatnum()==mat){
  	  real factor=1.;
	  if(fiss_mat[j]<0)factor=-1.;
          fiss_abs+=density_data[i][k]*(xsc_1g[i][k]+xsf_1g[i][k]+xsn2n_1g[i][k])*factor;
	};
      };
    };
    cout<<fert_cap/fiss_abs<<" ";
    //
    if(i!=burn_step)cout<<total_flux[i][0]<<" ";
    cout<<"\n";
  };
};

void Burner::ShowFourFactor()
{
  cout<<"#\n# +++ Time-dependent eigenvalue and four-factor +++\n";
  cout<<"#  (day)       ";
  cout<<"(GWd/t)     ";
  cout<<"(k-inf)     ";
  cout<<"(Epsilon)   ";
  cout<<"(p)         ";
  cout<<"(f)         ";
  cout<<"(eta)       ";
  cout<<"\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<burn_step+1;i++){
    if(i<10)cout<<" ";
    cout<<i<<" "<<acday[i]<<" "<<acburn[i];
    cout<<" "<<keff[i]<<" ";
    for(int j=0;j<4;j++){
      cout<<four_factor[j][i]<<" ";
    };
    cout<<"\n";
  };
};

void Burner::ShowNeutronFluxHistory()
{
  cout<<"#\n# +++ Time-dependent neutron flux [/cm2/s] +++\n";
  cout<<"# (day)       (GWd/t)     (fuel)      (clad)\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  //cout.precision(10);
  for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
    cout<<total_flux[i][0]<<" ";
    cout<<total_flux_clad[i][0]<<" ";
    cout<<"\n";
  };
};

void Burner::ShowNuclideList()
{
  cout<<"#\n# List of nuclides included in the fuel region\n#\n";
  for(int i=0;i<nucn;i++){
    WriteOut(i,5);
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    string name=midt.Name(matno);
    WriteOut(matno,8);
    cout<<"  "<<name<<"\n";
  };
};

//void Burner::ShowNumberDensityChange(Burnup &bu, real limit)
void Burner::ShowNumberDensityChange(real limit)
{
  cout<<"#\n#+++ Number density before/after burnup +++\n#\n";
  cout.setf(ios::scientific);
  cout.precision(7);
  //cout.precision(15);
  for(int ii=0;ii<2;ii++){
    for(int i=0;i<nucn;i++){
      real den_b=density_data[0][i];
      real den_a=density_data[burn_step][i];
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      if(((ii==0&&matno>=900000)||(ii==1&&matno<900000))&&(den_a>limit||den_b>limit)){
        //real dc=bu.GetBurnupChain().GetDecayConstant(matno);
	cout<<"# "<<matno<<" "<<midt.Name(matno)<<"   "<<den_b<<"   "<<den_a<<"\n";
	//cout<<" "<<matno<<" "<<den_a<<"\n";
	//cout<<" "<<matno<<" "<<dc<<"\n";
	//cout<<"   "<<den_b<<"   "<<den_a<<"  "<<den_a*(xsc_1g[burn_step][i]+xsf_1g[burn_step][i])<<"\n";
	//cout<<den_a<<"\n";
      };
    };
  };
};

void Burner::ShowWeightChange(Burnup &bu, real limit)
{
  real avo=0.60221367;

  cout<<"#\n#+++ Weight before/after burnup [g/cm] +++\n#\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  real sum_hm_a=0.;
  real sum_hm_b=0.;
  real sum_fp_a=0.;
  real sum_fp_b=0.;
  for(int ii=0;ii<2;ii++){
    for(int i=0;i<nucn;i++){
      real den_b=density_data[0][i];
      real den_a=density_data[burn_step][i];
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      if(((ii==0&&matno>=900000)||(ii==1&&matno<900000))&&(den_a>limit||den_b>limit)){
	real aw=bu.GetAtomicWeight(matno);
        real wgt_b=den_b*fuel_vol/avo*aw;
        real wgt_a=den_a*fuel_vol/avo*aw;
	if(matno>=900000){
	  sum_hm_a+=wgt_a;
	  sum_hm_b+=wgt_b;
	}else if(matno>=300000){
	  sum_fp_a+=wgt_a;
	  sum_fp_b+=wgt_b;
	};
        cout<<"# "<<matno<<" "<<midt.Name(matno);
        cout<<"   "<<wgt_b<<"   "<<wgt_a<<"\n";
      };
    };
  };
  cout<<"#\n";
  cout<<"# Total weight of heavy metal     (IZ>=90)    : "<<sum_hm_b<<" "<<sum_hm_a<<"\n";
  cout<<"# Total weight of fission product (90>IZ>=30) : "<<sum_fp_b<<" "<<sum_fp_a<<"\n";
};

void Burner::ShowNumberDensityChange(int nucnum, string *nuclist)
{
  cout<<"#\n#+++ Number density before/after burnup +++\n#\n";
  cout.setf(ios::scientific);
  cout.precision(7);
  for(int ii=0;ii<2;ii++){
    for(int j=0;j<nucnum;j++){
      real den_init=0.;
      real den_final=0.;
      int sz=nuclist[j].size();
      real hm=true;
      if(sz>2){
        int idt=midt.ID(nuclist[j]);
        for(int i=0;i<nucn;i++){
          int matno=med[0].GetNuclideInTurn(i).GetMatnum();
  	  if(matno==idt){
	    den_init+=density_data[0][i];
	    den_final+=density_data[burn_step][i];
	  };
        };
	if(idt<900000)hm=false;
      }else{
        for(int i=0;i<nucn;i++){
          int matno=med[0].GetNuclideInTurn(i).GetMatnum();
          string name=midt.Name(matno);
  	  if(nuclist[j]==name.substr(0,sz)){
	    den_init+=density_data[0][i];
	    den_final+=density_data[burn_step][i];
	    if(matno<900000)hm=false;
	  };
        };
      };
      if(((ii==0&&hm)||(ii==1&&!hm))&&(den_init!=0.||den_final!=0.)){
        cout<<" ";
        WriteOut(nuclist[j],7);
        cout<<"   "<<den_init<<"   "<<den_final<<"\n";
      };
    };
  };
};

void Burner::ShowCrossSection(int step)
{
  cout<<"#\n# One-group cross section list (cap/fis/n2n/[sum])\n#\n";
  cout<<"#   burnup-step : "<<step<<"\n#\n";
  cout<<"#\n#       ";
  cout<<"[CAPTURE]   ";   
  cout<<"[FISSION]   ";   
  cout<<"[(N,2N)]    ";   
  cout<<"[ALL]\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int ii=0;ii<2;ii++){
    for(int i=0;i<nucn;i++){
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      if((ii==0&&matno>=900000)||(ii==1&&matno<900000)){
	string tmp=midt.Name(matno);
	cout<<tmp;
	for(int k=0;k<8-tmp.size();k++){cout<<" ";};
	cout<<xsc_1g[step][i]<<"  ";
	cout<<xsf_1g[step][i]<<"  ";
	cout<<xsn2n_1g[step][i]<<"  ";
	cout<<xsc_1g[step][i]+xsf_1g[step][i]+xsn2n_1g[step][i]<<"  ";
	cout<<"\n";
      };
    };
  };
};

void Burner::ShowCrossSection(int prt_nuc,string *prt_nuc_nam, bool reaction_rate, int digit)
{
  vector<int> prt_nuc_turn(prt_nuc);
  GetPrintNuclide(prt_nuc,prt_nuc_nam,prt_nuc_turn);

  cout.setf(ios::scientific);
  cout.precision(digit);
  for(int ii=0;ii<3;ii++){
    cout<<"#\n# ";
    if(!reaction_rate)cout<<"One-group ";
    if(ii==0){
      cout<<"CAPTURE";
    }else if(ii==1){
      cout<<"FISSION";
    }else{
      cout<<"(N,2N)";
    };
    if(reaction_rate){
      cout<<" reaction rate\n#\n";
    }else{
      cout<<" cross section\n#\n";
    };
    cout<<"#\n#(day)     (GWd/t)    ";    
    for(int j=0;j<prt_nuc;j++){
      cout<<prt_nuc_nam[j];
      int tmp=prt_nuc_nam[j].size();
      for(int k=0;k<11-tmp;k++){cout<<" ";};
    };
    cout<<"     (sum)\n";
    int maxstep=burn_step+1;
    if(reaction_rate)maxstep=burn_step;
    for(int i=0;i<maxstep;i++){
      cout<<acday[i]<<" "<<acburn[i]<<" ";
      real sum=0.;
      for(int j=0;j<prt_nuc;j++){
        real val=0.;
        if(ii==0){
          val=xsc_1g[i][prt_nuc_turn[j]];
	}else if(ii==1){
          val=xsf_1g[i][prt_nuc_turn[j]];
	}else{
          val=xsn2n_1g[i][prt_nuc_turn[j]];
	};
        real factor=1.;
        if(reaction_rate){
          //factor=total_flux[i][0]*fuel_vol*density_data[i][prt_nuc_turn[j]];
          //factor=total_flux[i][0]/power_factor[i][0]*fuel_vol*density_data[i][prt_nuc_turn[j]];
          factor=total_flux_unnormalized[i][0]*fuel_vol*density_data[i][prt_nuc_turn[j]];
        };
	cout<<val*factor<<" ";
	sum+=val*factor;
      };
      cout<<sum<<"\n";
    };
  };

};

void Burner::ShowCaptureToDecayRatio(int prt_nuc,string *prt_nuc_nam,Burnup &bu)
{
  cout<<"#\n# [Microscopic capture reaction rate]/[Decay constant]\n#\n";

  vector<int> prt_nuc_turn(prt_nuc);
  GetPrintNuclide(prt_nuc,prt_nuc_nam,prt_nuc_turn);

  vector<real> dc(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    dc[i]=bu.GetBurnupChain().GetDecayConstant(midt.ID(prt_nuc_nam[i]));
    if(dc[i]==0.){
      cout<<"# Warning in Burner::ShowCaptureToDecayRatio.\n";
      cout<<"# Decay constant is zero : "<<prt_nuc_nam[i]<<"\n";
      return;
    };
  };

  cout<<"#(day)     (GWd/t)    ";
  for(int j=0;j<prt_nuc;j++){
    cout<<prt_nuc_nam[j];
    int tmp=prt_nuc_nam[j].size();
    for(int k=0;k<11-tmp;k++){cout<<" ";};
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<prt_nuc;j++){
      cout<<xsc_1g[i][prt_nuc_turn[j]]*1e-24*total_flux[i][0]/dc[j]<<" ";
    };
    cout<<"\n";
  };
};

void Burner::ShowFuelNeutronFlux(int st,bool excel)
{
  if(st>burn_step){
    cout<<"# Error in Burner::ShowFuelNeutronFlux\n";
    exit(0);
  };

  cout<<"#\n# Neutron flux per lethargy\n#\n";
  cout<<"#   burnup step : "<<st<<"\n";
  cout<<"#\n";
  if(excel){
    cout<<"# Eenrgy      Fuel          Clad\n";
    cout<<"# [eV]        region        region\n";
  }else{
    cout<<"# Upper       Fuel          Clad\n";
    cout<<"# energy[eV]  region        region\n";
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<group;i++){
    real e0=med[0].GetEnband().get_dat(i);
    real e1=med[0].GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    cout<<e0<<"   "<<fuel_flux[st].get_dat(i)/letwid<<"   ";
    cout<<clad_flux[st].get_dat(i)/letwid<<"\n";
    if(excel){
      cout<<e1<<"   "<<fuel_flux[st].get_dat(i)/letwid<<"   ";
      cout<<clad_flux[st].get_dat(i)/letwid<<"\n";
    };
  };
};

void Burner::ShowAdjointNumberDensity(int prt_nuc,string *prt_nuc_nam)
{
  if(!adj_nuc_data){
    cout<<"# Warning in Burner::ShowAdjointNumberDensity.\n";
    cout<<"# Adjoint number density is not calculated.\n";
    return;
  };

  cout<<"# Time step-averaged adjoint number density\n";

  vector<int> prt_nuc_turn(prt_nuc,-1);
  GetPrintNuclide(prt_nuc,prt_nuc_nam,prt_nuc_turn);
  cout.setf(ios::scientific);
  cout.precision(5);

  /*
    cout<<"# COARSE time mesh\n";
    cout<<"#  ";
    for(int k=0;k<prt_nuc;k++){
    cout<<prt_nuc_nam[k]<<" ";
    };
    cout<<"\n";
    for(int i=0;i<burn_step;i++){
    cout<<i<<" ";
    for(int k=0;k<prt_nuc;k++){
    cout<<adj_nuc[i][0].get_dat(prt_nuc_turn[k])<<" ";
    };
    cout<<"\n";
    };
  */

  
  cout<<"#\n# FINE time mesh\n";
  cout<<"# [Adjoint] ";
  for(int k=0;k<prt_nuc;k++){
    cout<<prt_nuc_nam[k]<<" ";
  };
  cout<<"\n";
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      cout<<acday[i]+(acday[i+1]-acday[i])/sub_step*j<<" ";
      cout<<acburn[i]+(acburn[i+1]-acburn[i])/sub_step*j<<"   ";
      //cout<<i<<" "<<j<<" ";
      for(int k=0;k<prt_nuc;k++){
	cout<<adj_nuc[i][j].get_dat(prt_nuc_turn[k])<<" ";
      };
      cout<<"\n";
    };
  };

  cout<<"\n\n";
  cout<<"# [Contribution] ";
  for(int k=0;k<prt_nuc;k++){
    cout<<prt_nuc_nam[k]<<" ";
  };
  cout<<"\n";
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      cout<<acday[i]+(acday[i+1]-acday[i])/sub_step*j<<" ";
      cout<<acburn[i]+(acburn[i+1]-acburn[i])/sub_step*j<<"   ";
      //cout<<i<<" "<<j<<" ";
      for(int k=0;k<prt_nuc;k++){
	cout<<adj_nuc[i][j].get_dat(prt_nuc_turn[k])*fwd_nuc[i][j].get_dat(prt_nuc_turn[k])/target_val<<" ";
      };
      cout<<"\n";
    };
  };

};

void Burner::ShowNumberDensityHistoryAtomWise(string atom_name,Burnup &bu,string opt,bool sum_only)
{
  int sz=atom_name.size();

  vector<string> prt_nuc_nam_org;
  for(int i=0;i<nucn;i++){
    string name=midt.Name(med[0].GetNuclideInTurn(i).GetMatnum());
    int leng=name.size();
    string endchr=name.substr(leng-1,1);
    if(endchr=="M"||endchr=="m")leng--;
    leng-=3;
    if(name.substr(0,sz)==atom_name&&sz==leng)prt_nuc_nam_org.push_back(name);
  };

  int prt_nuc=prt_nuc_nam_org.size();
  string *prt_nuc_nam=new string[prt_nuc];
  for(int i=0;i<prt_nuc;i++){
    prt_nuc_nam[i]=prt_nuc_nam_org[i];
  };

  if(sum_only)cout<<"#\n# Atom-wise information : "<<atom_name<<"\n";
  ShowNumberDensityHistory(prt_nuc,prt_nuc_nam,bu,opt,sum_only);
  cout<<"\n\n";

  delete [] prt_nuc_nam;
};

void Burner::ShowNumberDensityHistoryOld(int prt_nuc,string *prt_nuc_nam,Burnup &bu,string opt,bool sum_only)
{
  real avo=0.60221367;
  real ev_to_j=1.60219e-19;

  vector<int> prt_nuc_turn(prt_nuc,-1);
  vector<int> prt_nuc_id(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    prt_nuc_id[i]=midt.ID(prt_nuc_nam[i]);
  };
  GetPrintNuclide(prt_nuc,prt_nuc_nam,prt_nuc_turn);

  cout<<"#\n#(day)     (GWd/t)      ";    
  if(opt=="nd_per_vol"){
    cout<<"N.D. [1e24/cm3]\n";
  }else if(opt=="nd"){
    cout<<"N.D. [1e24])\n";
  }else if(opt=="bq"){
    cout<<"Radioactivity [Bq]\n";
  }else if(opt=="bq_per_thm"){
    cout<<"Radioactivity [Bq/tHM]\n";
  }else if(opt=="kg_per_thm"){
    cout<<"Weight [kg/tHM]\n";
  }else if(opt=="w_per_thm"){
    cout<<"Heat [W/tHM]\n";
  }else if(opt=="w"){
    cout<<"Heat [W]\n";
  };
  cout<<"#                       ";
  cout.setf(ios::scientific);
  cout.precision(4);
  if(!sum_only){
    for(int j=0;j<prt_nuc;j++){
      cout<<prt_nuc_nam[j];
      int tmp=prt_nuc_nam[j].size();
      for(int k=0;k<11-tmp;k++){cout<<" ";};
    };

  };
  cout<<" (sum)\n";
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<"   ";
    real sum=0.;
    for(int j=0;j<prt_nuc;j++){
      real dd=density_data[i][prt_nuc_turn[j]];
      real dc=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[j]);
      real val=0.;
      if(opt=="nd_per_vol"){
        val=dd;
      }else if(opt=="nd"){
        val=dd*fuel_vol;
      }else if(opt=="bq"){
        val=dd*fuel_vol*dc*1e24;
      }else if(opt=="bq_per_thm"){
        val=dd*fuel_vol*dc*1e24/(hm_weight_init*1e-6);
      }else if(opt=="kg_per_thm"){
        real aw=bu.GetAtomicWeight(prt_nuc_id[j]);
	real mol=(dd*fuel_vol)/avo;
        real wt=aw*mol*1e-3;
	val=wt/(hm_weight_init*1e-6);
      }else if(opt=="w_per_thm"){
        real e=0.;
        for(int k=0;k<3;k++){
          e+=bu.GetBurnupChain().GetDecayEnergy(prt_nuc_id[j],k);
        };
        val=e*dd*1e24*fuel_vol*dc*ev_to_j/(hm_weight_init*1e-6);
      }else if(opt=="w"){
        real e=0.;
        for(int k=0;k<3;k++){
          e+=bu.GetBurnupChain().GetDecayEnergy(prt_nuc_id[j],k);
        };
        val=e*dd*1e24*fuel_vol*dc*ev_to_j;
      };
      if(!sum_only)cout<<val<<" ";
      sum+=val;
    };
    cout<<"  "<<sum<<"\n";
  };
};

GData Burner::ShowNumberDensityHistory(int prt_nuc,string *prt_nuc_nam,Burnup &bu,string opt,bool sum_only)
{
  GData gdata(prt_nuc);
  gdata.PutTagX("day");
  for(int i=0;i<prt_nuc;i++){
    string tmp=opt+"_"+prt_nuc_nam[i];
    gdata.PutTagY(i,tmp);
  };
  
  
  //int sz=density_data.size();
  //cout<<sz<<"\n"; exit(0);

  real avo=0.60221367;
  real ev_to_j=1.60219e-19;

  vector< vector<int> > each_nuc_turn(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    if(prt_nuc_nam[i]=="HM"){
      for(int j=0;j<nucn;j++){
	int matid=med[0].GetNuclideInTurn(j).GetMatnum();
	if(matid>=800000){
	  each_nuc_turn[i].push_back(j);          
	};
      };
    }else if(prt_nuc_nam[i]=="FP"){
      for(int j=0;j<nucn;j++){
	int matid=med[0].GetNuclideInTurn(j).GetMatnum();
	if(matid<800000&&matid>=310000){
	  each_nuc_turn[i].push_back(j);          
	};
      };
    }else if(prt_nuc_nam[i]=="ALL"||prt_nuc_nam[i]=="All"){
      for(int j=0;j<nucn;j++){
        each_nuc_turn[i].push_back(j);          
      };
    }else{
      int sz=prt_nuc_nam[i].size();
      if(sz<=2){
	// (atom-wise)
        for(int j=0;j<nucn;j++){
          string name=midt.Name(med[0].GetNuclideInTurn(j).GetMatnum());
          int leng=name.size();
          string endchr=name.substr(leng-1,1);
          if(endchr=="M"||endchr=="m")leng--;
          leng-=3;
          if(name.substr(0,sz)==prt_nuc_nam[i]&&sz==leng)each_nuc_turn[i].push_back(j);
        };
      }else{
	// (nuclide-wise)
	for(int j=0;j<nucn;j++){
          string name=midt.Name(med[0].GetNuclideInTurn(j).GetMatnum());
	  if(name==prt_nuc_nam[i])each_nuc_turn[i].push_back(j);
	};
      };
    };
  };

  cout<<"#\n#(day)     (GWd/t)      ";    
  if(opt=="nd_per_vol"){
    cout<<"N.D. [1e24/cm3]\n";
  }else if(opt=="nd"){
    cout<<"N.D. [1e24])\n";
  }else if(opt=="bq"){
    cout<<"Radioactivity [Bq]\n";
  }else if(opt=="bq_per_thm"){
    cout<<"Radioactivity [Bq/tHM]\n";
  }else if(opt=="kg_per_thm"){
    cout<<"Weight [kg/tHM]\n";
  }else if(opt=="w"){
    cout<<"Heat [W]\n";
  }else if(opt=="w_per_thm"){
    cout<<"Heat [W/tHM]\n";
  }else if(opt=="w_gamma_per_thm"){
    cout<<"Gamma heat [W/tHM]\n";
  }else if(opt=="w_beta_per_thm"){
    cout<<"Beta Heat [W/tHM]\n";
  }else if(opt=="sv_per_thm_ing"){
    cout<<"Ingestion toxicity [Sv/tHM]\n";
  }else if(opt=="sv_per_thm_inh"){
    cout<<"Inhalation toxicity [Sv/tHM]\n";
  }else if(opt=="sv_ing"){
    cout<<"Ingestion toxicity [Sv]\n";
  }else if(opt=="sv_inh"){
    cout<<"Inhalation toxicity [Sv]\n";
  };
  cout<<"#                       ";
  cout.setf(ios::scientific);
  cout.precision(4);
  if(!sum_only){
    for(int j=0;j<prt_nuc;j++){
      cout<<"  (";
      WriteOut(j+3,3);
      cout<<")    ";
    };
    cout<<"\n";
    cout<<"#                       ";
    for(int j=0;j<prt_nuc;j++){
      cout<<prt_nuc_nam[j];
      int tmp=prt_nuc_nam[j].size();
      for(int k=0;k<11-tmp;k++){cout<<" ";};
    };
  };
  cout<<"  (sum)\n";
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<"   ";
    gdata.push_back_x(acday[i]);
    real sum=0.;
    for(int j=0;j<prt_nuc;j++){
      int sz=each_nuc_turn[j].size();
      real val=0.;
      for(int k=0;k<sz;k++){
        int tt=each_nuc_turn[j][k];
        int id=med[0].GetNuclideInTurn(tt).GetMatnum();
        real dd=density_data[i][tt];
        real dc=bu.GetBurnupChain().GetDecayConstant(id);
        if(opt=="nd_per_vol"){
          val+=dd;
        }else if(opt=="nd"){
          val+=dd*fuel_vol;
        }else if(opt=="bq"){
          val+=dd*fuel_vol*dc*1e24;
        }else if(opt=="bq_per_thm"){
          val+=dd*fuel_vol*dc*1e24/(hm_weight_init*1e-6);
        }else if(opt=="kg_per_thm"){
          real aw=bu.GetAtomicWeight(id); 
    	  real mol=(dd*fuel_vol)/avo;
          real wt=aw*mol*1e-3;
	  val+=wt/(hm_weight_init*1e-6);
        }else if(opt=="w_per_thm"){
          real e=0.;
          for(int k=0;k<3;k++){
            e+=bu.GetBurnupChain().GetDecayEnergy(id,k);
          };
          val+=e*dd*1e24*fuel_vol*dc*ev_to_j/(hm_weight_init*1e-6);
        }else if(opt=="w"){
          real e=0.;
          for(int k=0;k<3;k++){
            e+=bu.GetBurnupChain().GetDecayEnergy(id,k);
          };
          val+=e*dd*1e24*fuel_vol*dc*ev_to_j;
        }else if(opt=="w_gamma_per_thm"){
          real e=bu.GetBurnupChain().GetDecayEnergy(id,1);
          val+=e*dd*1e24*fuel_vol*dc*ev_to_j/(hm_weight_init*1e-6);
        }else if(opt=="w_beta_per_thm"){
          real e=bu.GetBurnupChain().GetDecayEnergy(id,0);
          val+=e*dd*1e24*fuel_vol*dc*ev_to_j/(hm_weight_init*1e-6);
        }else if(opt=="sv_per_thm_ing"){
          real bq=dd*fuel_vol*dc*1e24/(hm_weight_init*1e-6);
	  val+=bq*bu.GetBurnupChain().GetDoseCoefficientIngestion(id);
        }else if(opt=="sv_per_thm_inh"){
          real bq=dd*fuel_vol*dc*1e24/(hm_weight_init*1e-6);
	  val+=bq*bu.GetBurnupChain().GetDoseCoefficientInhalation(id);
        }else if(opt=="sv_ing"){
          real bq=dd*fuel_vol*dc*1e24;
	  val+=bq*bu.GetBurnupChain().GetDoseCoefficientIngestion(id);
        }else if(opt=="sv_inh"){
          real bq=dd*fuel_vol*dc*1e24;
	  val+=bq*bu.GetBurnupChain().GetDoseCoefficientInhalation(id);
        };
      };
      if(!sum_only)cout<<val<<" ";
      gdata.push_back_y(j,val);
      sum+=val;
    };
    cout<<"  "<<sum<<"\n";
  };

  return gdata;
};

void Burner::ShowRadioactivityRatio(Burnup &bu,string nuc1,string nuc2)
{
  string *prt_nuc_nam=new string[2];
  prt_nuc_nam[0]=nuc1;
  prt_nuc_nam[1]=nuc2;
  vector<int> prt_nuc_turn(2,-1);
  vector<int> prt_nuc_id(2);
  prt_nuc_id[0]=midt.ID(nuc1);
  prt_nuc_id[1]=midt.ID(nuc2);
  GetPrintNuclide(2,prt_nuc_nam,prt_nuc_turn);

  delete [] prt_nuc_nam;

  int idt1=prt_nuc_turn[0];
  int idt2=prt_nuc_turn[1];

  cout<<"#\n#(day)     (GWd/t)      (Bq-ratio : "<<nuc1<<"/"<<nuc2<<")\n#\n";    
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step+1;i++){
    real dd1=density_data[i][idt1];
    real dd2=density_data[i][idt2];
    real dc1=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[0]);
    real dc2=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[1]);
    if(dd2>0.){
      real val=(dd1*dc1)/(dd2*dc2);
      cout<<acday[i]<<" "<<acburn[i]<<"   ";
      cout<<val<<"\n";
    };
  };
};

real Burner::GetRadioactivityRatio(Burnup &bu,string nuc1,string nuc2)
{
  string *prt_nuc_nam=new string[2];
  prt_nuc_nam[0]=nuc1;
  prt_nuc_nam[1]=nuc2;
  vector<int> prt_nuc_turn(2,-1);
  vector<int> prt_nuc_id(2);
  prt_nuc_id[0]=midt.ID(nuc1);
  prt_nuc_id[1]=midt.ID(nuc2);
  GetPrintNuclide(2,prt_nuc_nam,prt_nuc_turn);

  delete [] prt_nuc_nam;

  int idt1=prt_nuc_turn[0];
  int idt2=prt_nuc_turn[1];

  int i=burn_step;
  real dd1=density_data[i][idt1];
  real dd2=density_data[i][idt2];
  real dc1=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[0]);
  real dc2=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[1]);
  return  (dd1*dc1)/(dd2*dc2);
};

real Burner::GetRadioactivity(Burnup &bu,string nuc1)
{
  string *prt_nuc_nam=new string[1];
  prt_nuc_nam[0]=nuc1;
  vector<int> prt_nuc_turn(1,-1);
  vector<int> prt_nuc_id(1);
  prt_nuc_id[0]=midt.ID(nuc1);
  GetPrintNuclide(1,prt_nuc_nam,prt_nuc_turn);
  delete [] prt_nuc_nam;

  int idt1=prt_nuc_turn[0];
  int i=burn_step;
  real dd1=density_data[i][idt1];
  real dc1=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[0]);
  return  dd1*dc1;
};

void Burner::ShowHeavyMetalWeightRatio(Burnup &bu)
{
  real avo=0.60221367;

  int nucnam_u=3;
  int nucnam_tru=1+5+3+5;
  int nucnam=nucnam_u+nucnam_tru;
  string hm_name[]={
    "U235","U236","U238",
    "Np237",
    "Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242m","Am243",
    "Cm242","Cm243","Cm244","Cm245","Cm246"
  };
  int p9_pos=5;

  vector<int> prt_nuc_turn(nucnam);
  vector<int> prt_nuc_id(nucnam);
  for(int i=0;i<nucnam;i++){
    int id=midt.ID(hm_name[i]);
    prt_nuc_id[i]=id;
    for(int j=0;j<nucn;j++){
      if(med[0].GetNuclideInTurn(j).GetMatnum()==id){
	prt_nuc_turn[i]=j;
	break;
      };
    };
  };

  vector<real> wgt(nucnam);
  cout<<"#\n# Heavy metal weight [g/cm]\n";
  cout<<"#    (before / after burnup)\n#\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  real sum_b=0.;
  real sum_a=0.;
  real sumfp_b=0.;
  real sumfp_a=0.;
  real del_p9=0.;
  real sum=0.;
  for(int i=0;i<nucnam;i++){
    real aw=bu.GetAtomicWeight(prt_nuc_id[i]);
    cout<<hm_name[i]<<" ";
    real wgt_b=density_data[0][prt_nuc_turn[i]]*fuel_vol/avo*aw;
    real wgt_a=density_data[burn_step][prt_nuc_turn[i]]*fuel_vol/avo*aw;
    wgt[i]=wgt_a;
    sum_b+=wgt_b;
    sum_a+=wgt_a;
    sumfp_b+=wgt_b;
    sumfp_a+=wgt_a;
    if(i>=nucnam_u)sum+=wgt_a;
    if(i==p9_pos)del_p9=wgt_b-wgt_a;
    cout<<wgt_b<<" "<<wgt_a<<" ";
    cout<<"\n";
  };
  cout<<"(Sum)   "<<sum_b<<" "<<sum_a<<"\n";

  cout.setf(ios::showpoint);
  cout.precision(5);
  cout<<"\n\n# TRU concentration [wt%]\n";
  cout<<"#   (left  : numerical result)\n";
  cout<<"#   (right : Pu-239 is added to preserve Pu-239 inventory)\n#\n";
  for(int i=nucnam_u;i<nucnam;i++){
    cout<<hm_name[i]<<"  "<<wgt[i]/sum*100<<" ";
    if(i!=p9_pos){
      cout<<wgt[i]/(sum+del_p9)*100<<" ";
    }else{
      cout<<(wgt[i]+del_p9)/(sum+del_p9)*100<<" ";
    };
    cout<<"\n";
  };

  //
  vector<real> hm_weight(nucn,0.);
  real wgtsum=0.;
  for(int i=0;i<nucn;i++){
    int matid=med[0].GetNuclideInTurn(i).GetMatnum();
    //if(matid>810000){
    if(matid>810000&&matid!=922350&&matid!=922380){
      real aw=bu.GetAtomicWeight(matid);
      if(aw==0.){
	int iz,ia,il;
	midt.GetParameter(matid,iz,ia,il);
	aw=real(ia);
      };
      hm_weight[i]=density_data[burn_step][i]*fuel_vol/avo*aw;
      wgtsum+=hm_weight[i];
      //cout<<matid<<" "<<aw<<" "<<density_data[burn_step][i]<<"\n";
    };
  };

  cout.setf(ios::showpoint);
  cout.precision(5);
  cout<<"\n\n# HM weight rate [wt%]\n";
  cout<<"#   (Tl, Pb, Bi, ...)\n";
  cout<<"#\n";
  for(int i=0;i<nucn;i++){
    if(hm_weight[i]>0.){
      int matid=med[0].GetNuclideInTurn(i).GetMatnum();
      cout<<matid<<" "<<midt.Name(matid)<<" "<<hm_weight[i]/wgtsum<<"\n";
    };
  };

  /*
  int cnt=0;
  cout<<"\n\n";
  for(int i=0;i<nucn;i++){
    if(hm_weight[i]>0.){
      int matid=med[0].GetNuclideInTurn(i).GetMatnum();
      if(matid!=922350&&matid!=922380){
	//cout<<midt.Name(matid)<<", ";
	cout<<matid<<", ";
	cnt++;
	if(cnt==5){
	  cnt=0;
	  cout<<"\n";
	};
      };
    };
  };
  cout<<"\n\n";

  cnt=0;
  for(int i=0;i<nucn;i++){
    if(hm_weight[i]>0.){
      int matid=med[0].GetNuclideInTurn(i).GetMatnum();
      if(matid!=922350&&matid!=922380){
	cout<<hm_weight[i]/wgtsum<<", ";
	cnt++;
	if(cnt==5){
	  cnt=0;
	  cout<<"\n";
	};
      };
    };
  };
  cout<<"\n\n";
  */
};

void Burner::ShowHeavyMetalNumberDensityRatio()
{
  cout<<"#\n# Heavy metal number density ratio\n#\n";

  vector<int> nucid;
  vector<real> den;
  real densum=0.;
  int nucn=med[0].GetNucnum();
  for(int i=0;i<nucn;i++){
    int matid=med[0].GetNuclideInTurn(i).GetMatnum();
    if(matid>900000&&matid<990000){
      nucid.push_back(matid);
      real tmp=med[0].GetNuclideInTurn(i).GetDensity();
      den.push_back(tmp);
      densum+=tmp;      
    };
  };

  cout.setf(ios::showpoint);
  cout.precision(6);
  int nn=nucid.size();
  for(int i=0;i<nn;i++){
    cout<<nucid[i]<<" "<<den[i]/densum<<"\n";
  };
};

void Burner::GetHeavyMetalWeightRatio(Burnup &bu, vector<real> &output)
{
  real avo=0.60221367;

  vector<real> hm_weight(nucn,0.);
  real wgtsum=0.;
  for(int i=0;i<nucn;i++){
    int matid=med[0].GetNuclideInTurn(i).GetMatnum();
    if(matid>810000&&matid!=922350&&matid!=922380){
      real aw=bu.GetAtomicWeight(matid);
      if(aw==0.){
	int iz,ia,il;
	midt.GetParameter(matid,iz,ia,il);
	aw=real(ia);
      };
      hm_weight[i]=density_data[burn_step][i]*fuel_vol/avo*aw;
      wgtsum+=hm_weight[i];
    };
  };

  for(int i=0;i<nucn;i++){
    if(hm_weight[i]>0.){
      output.push_back(hm_weight[i]/wgtsum);
    };
  };

};

void Burner::ShowNuclideWiseContributionForFission(int prt_nuc,string *prt_nuc_nam)
{
  vector<int> prt_nuc_turn(prt_nuc,-1);
  vector<int> prt_nuc_id(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    prt_nuc_id[i]=midt.ID(prt_nuc_nam[i]);
  };
  GetPrintNuclide(prt_nuc,prt_nuc_nam,prt_nuc_turn);

  vector<real> fis(prt_nuc);

  cout<<"#\n";
  cout<<"# Time dependent nuclide-wise contribution to fission reaction\n";
  cout<<"#\n#(day)    (GWd/t)   ";    
  for(int j=0;j<prt_nuc;j++){
    cout<<prt_nuc_nam[j]<<" ";
    int tmp=prt_nuc_nam[j].size();
    for(int k=0;k<9-tmp;k++){cout<<" ";};
  };
  cout<<" (sum)\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<prt_nuc;j++){
      fis[j]=0.;
    };
    real tot=0.;
    real totp=0.;
    for(int j=0;j<nucn;j++){
      real rrf=xsf_1g[i][j]*density_data[i][j];
      tot+=rrf;
      for(int k=0;k<prt_nuc;k++){
	if(prt_nuc_turn[k]==j){
          fis[k]=rrf;
	  totp+=rrf;
	};
      };
    };
    for(int j=0;j<prt_nuc;j++){
      cout<<fis[j]/tot<<" ";
    };
    cout<<" ("<<totp/tot<<")";
    cout<<"\n";
  };
};

void Burner::ShowPlutonium(Burnup &bu)
{
  real avo=0.60221367;

  cout<<"\n\n#  (day, GWd/t,   Ratio_Annihilation,  Amount_Annihilation,  Amount_Annihilation/BU,  MA-Amount\n";
  cout.setf(ios::scientific);
  cout.precision(4);

  vector<real> pu_w(burn_step+1);
  vector<real> ma_w(burn_step+1);
  for(int i=0;i<burn_step+1;i++){
    real sum=0.;
    real sum2=0.;
    for(int j=0;j<nucn;j++){
      int id=med[0].GetNuclideInTurn(j).GetMatnum();
      real dd=density_data[i][j];
      real val=0.;
      real aw=bu.GetAtomicWeight(id);
      real mol=(dd*fuel_vol)/avo;
      real wt=aw*mol*1e-3;
      val=wt/(hm_weight_init*1e-6);
      if(id>=940000&&id<950000){
        sum+=val;
      }else if(id>=930000){
        sum2+=val;
      };
    };
    pu_w[i]=sum;
    ma_w[i]=sum2;
  };

  for(int i=0;i<burn_step+1;i++){
    cout<<"# ";
    cout<<acday[i]<<" "<<acburn[i]<<"   ";
    real pu2=pu_w[i]-pu_w[0];
    real pu1=pu2/pu_w[0];
    real pu3=0.;
    if(acburn[i]>0.)pu3=pu2/acburn[i];
    cout<<pu1<<" "<<pu2<<" "<<pu3<<" "<<ma_w[i]<<"\n";
  };
};


void Burner::ShowMA(Burnup &bu)
{
  real avo=0.60221367;

  cout<<"\n\n#  (day, GWd/t,   Amount,  Amount/BU,  Amount/Pu_Annihilation\n";
  cout.setf(ios::scientific);
  cout.precision(4);

  vector<real> pu_w(burn_step+1);
  vector<real> ma_w(burn_step+1);
  for(int i=0;i<burn_step+1;i++){
    real sum_ma=0.;
    real sum_pu=0.;
    for(int j=0;j<nucn;j++){
      int id=med[0].GetNuclideInTurn(j).GetMatnum();
      real dd=density_data[i][j];
      real val=0.;
      real aw=bu.GetAtomicWeight(id);
      real mol=(dd*fuel_vol)/avo;
      real wt=aw*mol*1e-3;
      val=wt/(hm_weight_init*1e-6);
      if(id>=940000&&id<950000){
        sum_pu+=val;
      }else if(id>=930000){
        sum_ma+=val;
      };
    };
    pu_w[i]=sum_pu;
    ma_w[i]=sum_ma;
  };

  for(int i=0;i<burn_step+1;i++){
    cout<<"# ";
    cout<<acday[i]<<" "<<acburn[i]<<"   ";
    real ma1=ma_w[i]-ma_w[0];
    real ma2=0.;
    if(acburn[i]>0.)ma2=ma1/acburn[i];
    real ma3=0.;
    if(acburn[i]>0.)ma3=ma1/abs(pu_w[i]-pu_w[0]);
    cout<<ma1<<" "<<ma2<<" "<<ma3<<"\n";
  };
};


void Burner::ShowGroupDependentCrossSection(string nucname,bool excel)
{
  int matid=midt.ID(nucname);

  if(!med[0].ExistNuclide(matid)){
    cout<<"# Burner::ShowGroupDependentCrossSection.\n";
    cout<<"# No nuclide : "<<nucname<<"\n";
    exit(0);
  };

  cout<<"#\n# Group-dependent cross section of "<<nucname<<"\n#\n";
  if(excel){
    cout<<"# (energy)  (fission)   (capture)   (nu-value)\n"; 
  }else{
    cout<<"# (upper)   (fission)   (capture)   (nu-value)\n"; 
    cout<<"# (energy)\n";
  };
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int g=0;g<group;g++){
    real e0=med[0].GetEnband().get_dat(g);
    real e1=med[0].GetEnband().get_dat(g+1);
    real fis=med[0].GetNuclide(matid).GetMicxs().GetData1d(sigf).get_dat(g);
    real cap=med[0].GetNuclide(matid).GetMicxs().GetData1d(sigc).get_dat(g);
    real vnu=xslib.GetLibData(matid).GetXSData().GetData1d(nu).get_dat(g);
    cout<<e0<<" "<<fis<<" "<<cap<<" "<<vnu<<"\n";
    if(excel)cout<<e1<<" "<<fis<<" "<<cap<<" "<<vnu<<"\n";
  };
  
};

// +++ animation

void Burner::ShowRadioactivityHistoryForAnimation(Burnup &bu)
{
  ofstream fout;
  fout.setf(ios::scientific);
  fout.precision(5);

  int graph=0; 
  int z=0;
  int a=0;
  int n=0;
  int l=0;
  int z_x=0;
  int a_x=0;
  int l_x=0;
  int z_e=0;
  int a_e=0;
  int l_e=0;
  fout.open("animation.data",ios::out);
  for(int i=0;i<nucn;i++){
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    midt.GetParameter(matno,z,a,l);
    real dc=bu.GetBurnupChain().GetDecayConstant(matno);
    n=a-z;
    if(l==0)fout<<std::setfill('0')<<std::setw(3)<<std::right<<z<<" "<<std::setfill('0')<<std::setw(3)<<std::right<<n<<" ";
    for(int j=0;j<burn_step+1;j++){
      if(l==0){
	real rad=density_data[j][i]*fuel_vol*dc*1e24;
	for(int k=1;k<3;k++){
	  for(int p=0;p<nucn;p++){
	    int matno_x=med[0].GetNuclideInTurn(p).GetMatnum();
	    midt.GetParameter(matno_x,z_x,a_x,l_x);
	    if(l_x==k&&z_x==z&&a_x==a)rad+=density_data[j][p]*fuel_vol*dc*1e24;
	  };
	};
	fout<<rad+1.0<<" ";
	if(i==0)graph++;
      };
    };
    if(l==0)fout<<"\n";
  };
  for(int ii=1;ii<101;ii++){
    for(int jj=1;jj<161;jj++){
      for(int kk=0;kk<nucn;kk++){
	int matno_e=med[0].GetNuclideInTurn(kk).GetMatnum();
	midt.GetParameter(matno_e,z_e,a_e,l_e);
	if(z_e==ii&&(a_e-z_e)==jj)goto OUT;
      };
      fout<<std::setfill('0')<<std::setw(3)<<std::right<<ii<<" "<<std::setfill('0')<<std::setw(3)<<std::right<<jj<<" ";
      for(int pp=0;pp<burn_step+1;pp++){
	fout<<1.0<<" ";
      };
      fout<<"\n";
    OUT:;
    };
  };
  fout<<"#Graphs: "<<graph<<"\n";	
  fout.close();
  cout<<"Graphs: "<<graph<<"\n";
};

void Burner::ShowDataForAnimation(Burnup &bu,string filename,string type)
{
  int zmin=92;
  int nmin=141;
  int zmax=96;
  int nmax=150;

  int z,a,l;

  ofstream fout;
  fout.open(filename.data(),ios::out);
  real minv=1e30;
  real maxv=0.;
  for(int st=0;st<burn_step+1;st++){

    vector< vector<real> > dat(zmax+1,vector<real>(nmax+1,0.));
    for(int i=0;i<nucn;i++){
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      real dc=bu.GetBurnupChain().GetDecayConstant(matno);
      midt.GetParameter(matno,z,a,l);
      int n=a-z;
      if(z<=zmax&&n<=nmax){
	real tmp=0.;
	if(type=="nd"){
	  tmp=density_data[st][i];
	}else if(type=="bq"){
	  tmp=density_data[st][i]*dc;
	}else if(type=="w"||type=="en"){
	  real e=0.;
	  for(int k=0;k<3;k++){
	    e+=bu.GetBurnupChain().GetDecayEnergy(matno,k);
	  };
	  tmp=density_data[st][i]*dc*e;
	};
	dat[z][n]+=tmp;
      };
    };


    for(int i=nmin;i<=nmax;i++){
      for(int jj=0;jj<2;jj++){
        for(int j=zmin;j<=zmax;j++){
          real tmp=dat[j][i];
          if(tmp>maxv)maxv=tmp;
	  if(tmp<minv&&tmp>0.)minv=tmp;
          fout<<i-0.5+jj<<" "<<j-0.5<<" "<<tmp<<"\n";
          fout<<i-0.5+jj<<" "<<j+0.5<<" "<<tmp<<"\n";
        };
        fout<<"\n";
      };
    };
    fout<<"\n\n";

  };

  cout<<"# step : "<<burn_step+1<<"\n";
  cout<<"# max. : "<<maxv<<"\n";
  cout<<"# min. : "<<minv<<"\n";

  fout.close();
};

void Burner::ShowDataForAnimationNew(Burnup &bu,string filename,string type)
{
  int zmin=92;
  int zmax=96;
  int amin=234;
  int amax=246;

  int z,a,l;

  ofstream fout;
  fout.open(filename.data(),ios::out);
  real minv=1e30;
  real maxv=0.;
  for(int st=0;st<burn_step+1;st++){

    vector< vector<real> > dat(zmax+1,vector<real>(amax+1,0.));
    for(int i=0;i<nucn;i++){
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      real dc=bu.GetBurnupChain().GetDecayConstant(matno);
      midt.GetParameter(matno,z,a,l);
      if(z<=zmax&&a<=amax){
	real tmp=0.;
	if(type=="nd"){
	  tmp=density_data[st][i];
	}else if(type=="bq"){
	  tmp=density_data[st][i]*dc;
	}else if(type=="w"||type=="en"){
	  real e=0.;
	  for(int k=0;k<3;k++){
	    e+=bu.GetBurnupChain().GetDecayEnergy(matno,k);
	  };
	  tmp=density_data[st][i]*dc*e;
	};
	dat[z][a]+=tmp;
      };
    };


    for(int i=amin;i<=amax;i++){
      for(int jj=0;jj<2;jj++){
        for(int j=zmin;j<=zmax;j++){
          real tmp=dat[j][i];
          if(tmp>maxv)maxv=tmp;
	  if(tmp<minv&&tmp>0.)minv=tmp;
          fout<<i-0.5+jj<<" "<<j-0.5<<" "<<tmp<<"\n";
          fout<<i-0.5+jj<<" "<<j+0.5<<" "<<tmp<<"\n";
        };
        fout<<"\n";
      };
    };
    fout<<"\n\n";

  };

  cout<<"# step : "<<burn_step+1<<"\n";
  cout<<"# max. : "<<maxv<<"\n";
  cout<<"# min. : "<<minv<<"\n";

  fout.close();
};

void Burner::CheckNumberDensityToRef(string mdir,string ss,int nucnum,string *nucname)
{
  Medium med_ref;
  med_ref.ReadFileNumberDensity(mdir,ss);
  
  real sum=0.;
  real sum2=0.;
  real max=0.;
  string maxname="";
  for(int i=0;i<nucnum;i++){
    real den=GetNuclideDensity(nucname[i],burn_step);
    real den_ref=med_ref.GetNuclide(midt.ID(nucname[i])).GetDensity();
    real dif=(den-den_ref)/den_ref;
    if(fabs(dif)>max){
      max=fabs(dif);
      maxname=nucname[i];
    };
    sum+=fabs(dif);
    sum2+=dif*dif;
    cout<<nucname[i]<<" "<<dif<<" "<<den<<" "<<den_ref<<"\n";
  };

  cout<<"# Average of absolute error      : "<<sum/nucnum<<"\n";
  cout<<"# Square root of squared average : "<<sqrt(sum2/nucnum)<<"\n";
  cout<<"# Maximum                        : "<<max<<" ("<<maxname<<")\n";
};

void Burner::SensitivityCalculationDirect(Burnup &bu, int target_num, string *target, int matid, int mt,string snsname)
{
  vector<SensitivityData> sns(target_num);

  int nstep=burn_step;
  if(mt!=102&&mt!=18){
    cout<<"# Error in SensitivityCalculationDirect.\n";
    cout<<"# Not yet be coded for MT="<<mt<<"\n";
    exit(0);
  };
  enum xstype xst=sigf;
  if(mt==102)xst=sigc;

  real factor=0.01;

  vector< vector<real> > response(target_num);
  for(int i=0;i<target_num;i++){
    response[i].resize(group+1);
  };

  // (initial ND storing)
  int nuc0=med[0].GetNucnum();
  int *mat0=new int[nuc0];
  real *den0=new real[nuc0];
  real temp0=med[0].GetNuclideInTurn(0).GetTemperature();
  for(int i=0;i<nuc0;i++){
    mat0[i]=med[0].GetNuclideInTurn(i).GetMatnum();
    den0[i]=med[0].GetNuclideInTurn(i).GetDensity();
  };

  for(int g=0;g<group+1;g++){
    real dxs;
    if(g!=group){
      dxs=GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).get_dat(g)*factor;
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).add_data(g,dxs);
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigt).add_data(g,dxs);
    };
    PutFuelData(nuc0,mat0,den0,temp0);
    Calculation(bu);
    for(int i=0;i<target_num;i++){
      response[i][g]=GetNuclideDensity(target[i],nstep);
    };
    if(g!=group){
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).add_data(g,-dxs);
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigt).add_data(g,-dxs);
    };
  };

  for(int i=0;i<target_num;i++){
    sns[i].PutName("dummy","dummy","dummy");
    sns[i].PutValue(response[i][group]);
    sns[i].PutGroup(group);
    sns[i].GetEnband().copy(GetXSLibrary().GetEnband());
    GroupData1D sns1d(group);
    for(int g=0;g<group;g++){
      real val=(response[i][g]-response[i][group])/(response[i][group]*factor);
      sns1d.put_data(g,val);
    };
    sns[i].PutSensitivity1D(matid,mt,sns1d);
  };

  for(int i=0;i<target_num;i++){
    sns[i].WriteFile("./",snsname+"."+target[i]);
  };

  delete[] mat0;
  delete[] den0;

};

void Burner::WriteFileFuelNumberDensity(string mdir, string ss, int cyc)
{
  if(cyc<0||cyc>burn_step){
    cout<<"# Error in Burner::WriteFileFuelNumberDensity.\n";
    cout<<"# The target cycle number "<<cyc<<" is inappropriate.\n";
    exit(0);
  };
  
  for(int i=0;i<nucn;i++){
    med[0].GetNuclideInTurn(i).PutDensity(density_data[cyc][i]);
  };

  med[0].WriteFileNumberDensity(mdir,ss);

  for(int i=0;i<nucn;i++){
    med[0].GetNuclideInTurn(i).PutDensity(density_data[burn_step][i]);
  };
};
void Burner::WriteFileNumberDensity(int cyc_num, int *cyc, string mdir, string filename, int digit)
{
  ofstream fout;
  mdir.append(filename);
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Error in Burner::WriteFileNumberDensity.\n";
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(0);
  };

  for(int i=0;i<cyc_num;i++){
    int cc=cyc[i];
    if(cc<0||cc>burn_step){
      cout<<"# Error in Burner::WriteFileNumberDensity.\n";
      cout<<"# Number density data at cycle "<<cc<<" do NOT exist.\n";
      exit(0);
    };
  };

  for(int i=0;i<cyc_num;i++){
    int cc=cyc[i];
    fout<<"      "<<cc<<"\n";
    fout.setf(ios::scientific);
    fout.precision(digit);
    fout<<"      "<<keff[cc]<<"\n";
    for(int i=0;i<nucn;i++){
      fout<<" "<<density_data[cc][i]<<"\n";
    };
  };

  fout.close();
};

void Burner::UseLibraryDataForNonExistingNuclide(Burnup &bu, string cbglibdir, string nucname)
{
  string libdir=cbglibdir+"CBGLIB/j4.107g.iwt7/";
  string namelist[]={"dummy"};
  int matidlist[]={0};

  ConstructingFuelNuclide(bu);

  for(int i=0;i<nucn;i++){
    int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
    bool exist=xslib.ExistLibData(nucid);
    if(!exist){
      namelist[0]=nucname;
      matidlist[0]=nucid;
      xslib.ReadFile(1,libdir,namelist,matidlist);
    };
  };
};

void Burner::ShowDominantCaptureReactionRate(real threshold)
{
  cout<<"#\n";
  cout<<"# Dominant capture reaction rate\n";
  cout<<"#    (threshold is "<<threshold<<")\n";
  cout<<"#\n";

  vector<bool> nuc_dom(nucn,false);

  cout.setf(ios::scientific);
  cout.precision(6);
  for(int i=0;i<burn_step;i++){
    for(int j=0;j<nucn;j++){
      real xs=xsc_1g[i][j];
      real den=density_data[i][j];
      real flx=total_flux_unnormalized[i][0]*fuel_vol;
      real reaction_rate=xs*den*flx;
      if(reaction_rate>threshold)nuc_dom[j]=true;
    };
  };

  int dom_num=0;
  for(int i=0;i<nucn;i++){
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    string tmp=midt.Name(matno);
    if(nuc_dom[i]){
      cout<<"# "<<tmp<<" "<<matno<<"\n";
      dom_num+=1;
    };
  };
  cout<<"#\n";
  cout<<"# The number of dominant nuclides : "<<dom_num<<"\n";
  cout<<"#\n";

  int count=0;
  for(int i=0;i<nucn;i++){

    if(nuc_dom[i]){
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      string tmp=midt.Name(matno);
      cout<<"#\n";
      cout<<"#   count : "<<count<<"\n";
      cout<<"# Nuclide : "<<tmp<<" ("<<matno<<")\n";
      cout<<"#\n";
      cout<<"# (Day)  (Burnup)\n";

      for(int j=0;j<burn_step;j++){
        cout<<acday[j]<<" "<<acburn[j]<<" ";
        real xs=xsc_1g[j][i];
        real den=density_data[j][i];
        real flx=total_flux_unnormalized[j][0]*fuel_vol;
	real reaction_rate=xs*den*flx;
	cout<<reaction_rate<<"\n";
      };
      cout<<"\n\n";
      count+=1;
    };

  };


  /*


  for(int i=0;i<burn_step;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<nucn;j++){
      int matno=med[0].GetNuclideInTurn(j).GetMatnum();
      string tmp=midt.Name(matno);
      real xs=xsc_1g[i][j];
      real den=density_data[i][j];
      real flx=total_flux_unnormalized[i][0]*fuel_vol;
    };
  };
  */

};

void Burner::GroupCollapsingDuringBurnup(string dirname, string filename, bool micro)
{
  if(group!=172){
    cout<<"# Error in Burner::GroupCollapsingDuringBurnup.\n";
    cout<<"# This method is hard-coded only for 172-group calculations.\n";
    exit(0);
  };
  collapsing=true;
  collapse_dir="./"+dirname+"/";
  collapse_filename=filename;
  collapsing_micro=micro;
};

void Burner::EstimateTRUcompositionByROM(Burnup &bu, real burnup, real coolyear, bool pwr, string filename)
{
  real avo=0.60221367;
  
  int num_tru=16;
  int tru_list[]={
    932370,932390,942380,942390,942400,    942410,942420,952410,952420,952421,
    952430,962420,962430,962440,962450,    962460
  };

  // Initial composition
  vector<real> tru_comp_ref_init(num_tru,0.);
  for(int i=0;i<nucn;i++){
    int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
    for(int j=0;j<num_tru;j++){
      if(nucid==tru_list[j]){
	real mol=(density_data[0][i]*fuel_vol)/avo;
	tru_comp_ref_init[j]=bu.GetAtomicWeight(nucid)*mol; // [g]
      };
    };
  };
  
  // After long-term cooling period
  vector<real> tru_comp_ref(num_tru,0.);
  for(int i=0;i<nucn;i++){
    int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
    for(int j=0;j<num_tru;j++){
      if(nucid==tru_list[j]){
	real mol=(density_data[burn_step][i]*fuel_vol)/avo;
	tru_comp_ref[j]=bu.GetAtomicWeight(nucid)*mol; // [g]
	//tru_comp_ref[j]=density_data[burn_step][i];		
      };
    };
  };
  vector<real> tru_comp_ref_org=tru_comp_ref;
  Normalize(tru_comp_ref);

  // Before long-term cooling period  
  vector<real> tru_comp_refb(num_tru,0.);
  for(int i=0;i<nucn;i++){
    int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
    for(int j=0;j<num_tru;j++){
      if(nucid==tru_list[j]){
	real mol=(density_data[burn_step-1][i]*fuel_vol)/avo;
	tru_comp_refb[j]=bu.GetAtomicWeight(nucid)*mol; // [g]
	//tru_comp_refb[j]=density_data[burn_step-1][i];
      };
    };
  };
  vector<real> tru_comp_refb_org=tru_comp_refb;  
  Normalize(tru_comp_refb);

  string reactor_type="PWR";
  if(!pwr)reactor_type="BWR";
  cout<<"#\n";
  cout<<"# <input information>\n";
  cout<<"#\n";
  cout<<"# Reactor type                  : "<<reactor_type<<"\n";
  cout<<"# Fuel burnup                   : "<<burnup<<" [GWD/t]\n";
  //cout<<"# Cooling year                  : "<<cooling_year<<" [year]\n";
  cout<<"# H/HM number density datio     : "<<h_hm_ratio<<"\n";
  cout<<"# Cooling-year                  : "<<coolyear<<" year\n";
  cout<<"# Initial weight of heavy metal : "<<hm_weight_init<<" [g]\n";
  cout<<"# Initial weight of U           : "<<hm_weight_init_u<<" [g]\n";
  cout<<"# Initial weight of TRU         : "<<hm_weight_init_tru<<" [g]\n";

  vector<real> tru;
  vector<real> trub;
  real u5=EstimateUenrichment(bu);
  cout<<"# Uranium-235 enrichment        : "<<u5<<" [wt%]\n";
  
  if(hm_weight_init_tru>1e-5){
    trub=EstimateTRUcompositionByROMForMOX(bu, burnup, 1., hm_weight_init_tru/hm_weight_init*100., pwr);
    tru=EstimateTRUcompositionByROMForMOX(bu, burnup, coolyear, hm_weight_init_tru/hm_weight_init*100., pwr);
  }else{
    trub=EstimateTRUcompositionByROMForUO2(bu, burnup, 1., pwr);
    tru=EstimateTRUcompositionByROMForUO2(bu, burnup, coolyear, pwr);
  };

  cout<<"#\n";
  cout<<"# <output information>\n";
  cout<<"# Cooling-year                  : 1year\n";
  cout<<"#\n";
  cout<<"#        [ROM]        [Original]   [ROM]-[Original]\n";
  int sz=tru.size();
  cout.setf(ios::scientific);
  cout.precision(6);
  for(int i=0;i<sz;i++){
    cout<<tru_list[i]<<" : ";
    cout<<trub[i]<<" "<<tru_comp_refb[i]<<" "<<trub[i]-tru_comp_refb[i]<<"\n"; // TRU composition ratio
  };

    cout<<"# Cooling-year                  : "<<int(coolyear)<<"year\n";
    cout<<"#\n";
    cout<<"#        [ROM]        [Original]   [ROM]-[Original]\n";
    cout.setf(ios::scientific);
    cout.precision(6);
    for(int i=0;i<sz;i++){
      cout<<tru_list[i]<<" : ";
      cout<<tru[i]<<" "<<tru_comp_ref[i]<<" "<<tru[i]-tru_comp_ref[i]<<"\n"; // TRU composition ratio
    };


  if(filename!=""){
    ofstream fout;
    fout.open(filename.data(),ios::out);
    if(fout.fail()){
      cout<<"# Error in Burner::EstimateTRUcompositionByROM.\n";
      exit(1);
    };
    fout.setf(ios::scientific);
    fout.precision(6);

      fout<<"# <input information>\n";
      fout<<"#\n";
      fout<<"# Reactor type                  : "<<reactor_type<<"\n";
      fout<<"# Fuel burnup                   : "<<burnup<<" [GWD/t]\n";
      //cout<<"# Cooling year                  : "<<cooling_year<<" [year]\n";
      fout<<"# H/HM number density datio     : "<<h_hm_ratio<<"\n";
      fout<<"# Initial weight of heavy metal : "<<hm_weight_init<<" [g]\n";
      fout<<"# Initial weight of U           : "<<hm_weight_init_u<<" [g]\n";
      fout<<"# Initial weight of TRU         : "<<hm_weight_init_tru<<" [g]\n";
      fout<<"# Uranium-235 enrichment        : "<<u5<<" [wt%]\n";
      fout<<"# Cooling-year                  : 1year\n";
      for(int i=0;i<sz;i++){
        fout<<tru_list[i]<<" ";
        fout<<trub[i]<<" "<<tru_comp_refb[i]<<" "<<trub[i]-tru_comp_refb[i]<<"\n"; // TRU composition ratio
      };
      fout<<"# Cooling-year                  : "<<int(coolyear)<<" year\n";
      for(int i=0;i<sz;i++){
        fout<<tru_list[i]<<" ";
        fout<<tru[i]<<" "<<tru_comp_ref[i]<<" "<<tru[i]-tru_comp_ref[i]<<"\n"; // TRU composition ratio
      };
      fout<<"# Absolute_weight_of_TRU_by_the_reference_model_[g]\n";
      fout<<"#    (initial)(after_burnup_with_1year_cooling)(after_burnup_with_long-term_cooling)\n";
      for(int i=0;i<sz;i++){
        fout<<tru_list[i]<<" ";
        fout<<tru_comp_ref_init[i]<<" ";
        fout<<tru_comp_refb_org[i]<<" ";	
        fout<<tru_comp_ref_org[i]<<" ";
	fout<<"\n";
      };
      fout.close();
  };
  
};

real Burner::EstimateUenrichment(Burnup &bu)
{
  // ... U-235 enrichment calculation
  real u_sum=0.;
  real u5;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=920000&&id<930000){
      u_sum+=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==922350)u5=density_data[0][i]*bu.GetAtomicWeight(id);
    };
  };
  u5=u5/u_sum*100; // [at%]
  //cout<<"# Uranium-235 enrichment        : "<<u5<<" [wt%]\n";

    return u5;
};

vector<real> Burner::EstimateTRUcompositionByROMForUO2(Burnup &bu, real burnup, real coolyear, bool pwr)
{
  // ... U-235 enrichment calculation
  real u_sum=0.;
  real u5;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=920000&&id<930000){
      u_sum+=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==922350)u5=density_data[0][i]*bu.GetAtomicWeight(id);
    };
  };
  u5=u5/u_sum*100; // [at%]

  // ... ROM calculations with RMforFR
  
  int input_data  = 4;
  int output_data = 16;
  int num_vib     = 7;

  RMforFR_LowOrderModel WR(input_data, output_data, input_data, num_vib);
    
  vector<string> filename_rm;
  for(int i=0;i<num_vib;i++){
    filename_rm.push_back("./LIB/polynomial0_func" + IntToString(i));
  };
  WR.ReadRegressionModel(filename_rm);

  WR.ReadPostProcessor("./LIB/SVD0_func");

  vector<double> xinp(input_data);

  xinp[0]=1.0;          // PWR
  if(!pwr)xinp[0]=0.0;  // BWR
  xinp[1]=u5;           // U-235 enrichment [wt%]
  xinp[2]=burnup;       // Burnup [GWD/t]  
  xinp[3]=h_hm_ratio;   // H/HM in mol ratio


  vector<double> tru=WR.GetRegressionResult(xinp);

  if(coolyear > 1.0){
    RMforFR_Coolingperiod uo2cooling(output_data, coolyear);
    uo2cooling.Readfilecp("./LIB/Coolingperiod");
    tru=uo2cooling.Coolingcal(tru);
  }

  return tru;
};


vector<real> Burner::EstimateTRUcompositionByROMForMOX(Burnup &bu, real burnup, real coolyear, real pu_enrichment, bool pwr)
{
    cout<<"# Pu enrichment                 : "<<pu_enrichment<<" [wt%]\n";

    // ... Pu-239 & -240 composition calculation
    real tru_sum=0.;
    real p9,p0,p1,p2;
    for(int i=0;i<nucn;i++){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      if(id>=930000&&id<990000){
        tru_sum+=density_data[0][i]*bu.GetAtomicWeight(id);
        if(id==942390)p9=density_data[0][i]*bu.GetAtomicWeight(id);
        if(id==942400)p0=density_data[0][i]*bu.GetAtomicWeight(id);
        if(id==942410)p1=density_data[0][i]*bu.GetAtomicWeight(id);
        if(id==942420)p2=density_data[0][i]*bu.GetAtomicWeight(id);
      };
    };

    p9=p9/tru_sum; // [at]
    p0=p0/tru_sum; // [at]
    p1=p1/tru_sum; // [at]
    p2=p2/tru_sum; // [at]
    cout<<"# Pu-239 ratio in TRU           : "<<p9<<" [-]\n";
    cout<<"# Pu-240 ratio in TRU           : "<<p0<<" [-]\n";
    cout<<"# Pu-241 ratio in TRU           : "<<p1<<" [-]\n";
    cout<<"# Pu-242 ratio in TRU           : "<<p2<<" [-]\n";

  int input_data  = 8;
  int output_data = 16;
  int num_vib     = 8;

  RMforFR_LowOrderModel WR(input_data, output_data, input_data, num_vib);
    
  vector<string> filename_rm;
  for(int i=0;i<num_vib;i++){
    filename_rm.push_back("./LIB/polynomial1_func" + IntToString(i));
  };
  WR.ReadRegressionModel(filename_rm);

  WR.ReadPostProcessor("./LIB/SVD1_func");

  vector<double> xinp(input_data);

  xinp[0]=1.0;           // PWR
  if(!pwr)xinp[0]=0.0;   // BWR

  xinp[1]=pu_enrichment; // Pu enrichment [wt%]
  xinp[2]=burnup;        // Burnup [GWD/t]
  xinp[3]=h_hm_ratio;    // H/HM in mol ratio
  xinp[4]=p9;       // Pu-239 atomic ratio in TRU
  xinp[5]=p0;       // Pu-240 atomic ratio in TRU
  xinp[6]=p1;       // Pu-241 atomic ratio in TRU
  xinp[7]=p2;       // Pu-242 atomic ratio in TRU

    int sz_xinp=xinp.size();
    for(int i=0;i<sz_xinp;i++){
        cout<<xinp[i]<<endl;
    }

  vector<double> tru=WR.GetRegressionResult(xinp);
    if(coolyear > 1.0){
        RMforFR_Coolingperiod moxcooling(output_data, coolyear);
        moxcooling.Readfilecp("./LIB/Coolingperiod");
        tru=moxcooling.Coolingcal(tru);
    }
  return tru;
};

// The following was replaced by one modified by Okuyama-san in 2022/5/21.
#if 0
void Burner::EstimateTRUcompositionByROM(Burnup &bu, real burnup, bool pwr, string filename)
{
  int num_tru=16;
  int tru_list[]={
    932370,932390,942380,942390,942400,    942410,942420,952410,952420,952421,
    952430,962420,962430,962440,962450,    962460
  };

  vector<real> tru_comp_ref(num_tru,0.);
  for(int i=0;i<nucn;i++){
    int nucid=med[0].GetNuclideInTurn(i).GetMatnum();
    for(int j=0;j<num_tru;j++){
      if(nucid==tru_list[j])tru_comp_ref[j]=density_data[burn_step][i];
    };
  };
  Normalize(tru_comp_ref);

  string reactor_type="PWR";
  if(!pwr)reactor_type="BWR";

  cout.setf(ios::scientific);
  cout.precision(6);
  
  cout<<"#\n";
  cout<<"# <input information>\n";
  cout<<"#\n";
  cout<<"# Reactor type                  : "<<reactor_type<<"\n";
  cout<<"# Fuel burnup                   : "<<burnup<<" [GWD/t]\n";
  //cout<<"# Cooling year                  : "<<cooling_year<<" [year]\n";
  cout<<"# H/HM number density datio     : "<<h_hm_ratio<<"\n";
  cout<<"# Initial weight of heavy metal : "<<hm_weight_init<<" [g]\n";
  cout<<"# Initial weight of U           : "<<hm_weight_init_u<<" [g]\n";
  cout<<"# Initial weight of TRU         : "<<hm_weight_init_tru<<" [g]\n";

  vector<real> tru;
  if(hm_weight_init_tru>1e-5){
    tru=EstimateTRUcompositionByROMForMOX(bu, burnup, hm_weight_init_tru/hm_weight_init*100., pwr);    
  }else{
    tru=EstimateTRUcompositionByROMForUO2(bu, burnup, pwr);        
  };

  cout<<"#\n";
  cout<<"# <output information>\n";
  cout<<"#\n";
  cout<<"#        [ROM]        [Original]   [ROM]-[Original]\n";
  int sz=tru.size();
  for(int i=0;i<sz;i++){
    cout<<tru_list[i]<<" : ";
    cout<<tru[i]<<" "<<tru_comp_ref[i]<<" "<<tru[i]-tru_comp_ref[i]<<"\n"; // TRU composition ratio
  };

  if(filename!=""){
    ofstream fout;
    fout.open(filename.data(),ios::out);
    if(fout.fail()){
      cout<<"# Error in Burner::EstimateTRUcompositionByROM.\n";
      exit(1);
    };
    fout.setf(ios::scientific);
    fout.precision(6);
    for(int i=0;i<sz;i++){
      fout<<tru_list[i]<<" ";
      fout<<tru[i]<<" "<<tru_comp_ref[i]<<" "<<tru[i]-tru_comp_ref[i]<<"\n"; // TRU composition ratio
    };
    fout.close();
  };
  
};

vector<real> Burner::EstimateTRUcompositionByROMForUO2(Burnup &bu, real burnup, bool pwr)
{
  // ... U-235 enrichment calculation
  real u_sum=0.;
  real u5;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=920000&&id<930000){
      u_sum+=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==922350)u5=density_data[0][i]*bu.GetAtomicWeight(id);
    };
  };
  u5=u5/u_sum*100; // [wt%]
  cout<<"# Uranium-235 enrichment        : "<<u5<<" [wt%]\n";

  // ... ROM calculations with RMforFR
  
  int input_data  = 4;
  int output_data = 16;
  int num_vib     = 8;

  RMforFR_LowOrderModel WR(input_data, output_data, input_data, num_vib);
    
  vector<string> filename_rm;
  for(int i=0;i<num_vib;i++){
    filename_rm.push_back("./LIB/polynomial0_func" + IntToString(i));
  };
  WR.ReadRegressionModel(filename_rm);

  WR.ReadPostProcessor("./LIB/SVD0_func");

  vector<double> xinp(input_data);

  xinp[0]=1.0;          // PWR
  if(!pwr)xinp[0]=0.0;  // BWR
  xinp[1]=u5;           // U-235 enrichment [wt%]
  xinp[2]=burnup;       // Burnup [GWD/t]  
  xinp[3]=h_hm_ratio;   // H/HM in mol ratio  

  vector<double> tru=WR.GetRegressionResult(xinp);
  return tru;
};


vector<real> Burner::EstimateTRUcompositionByROMForMOX(Burnup &bu, real burnup, real pu_enrichment, bool pwr)
{
  fout<<"# TRU enrichment                : "<<pu_enrichment<<" [wt%]\n";  

  // ... Pu-239 & -240 composition calculation
  real tru_sum=0.;
  real p9,p0,p1,p2;
  for(int i=0;i<nucn;i++){
    int id=med[0].GetNuclideInTurn(i).GetMatnum();
    if(id>=930000&&id<990000){
      tru_sum+=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==942390)p9=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==942400)p0=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==942410)p1=density_data[0][i]*bu.GetAtomicWeight(id);
      if(id==942420)p2=density_data[0][i]*bu.GetAtomicWeight(id);
    };
  };
  p9=p9/tru_sum; // [at]
  p0=p0/tru_sum; // [at]
  p1=p1/tru_sum; // [at]
  p2=p2/tru_sum; // [at]
  fout<<"# Pu-239 ratio in TRU           : "<<p9<<" [-]\n";
  fout<<"# Pu-240 ratio in TRU           : "<<p0<<" [-]\n";
  fout<<"# Pu-241 ratio in TRU           : "<<p1<<" [-]\n";
  fout<<"# Pu-242 ratio in TRU           : "<<p2<<" [-]\n";

  int input_data  = 8;
  int output_data = 16;
  int num_vib     = 8;

  RMforFR_LowOrderModel WR(input_data, output_data, input_data, num_vib);
    
  vector<string> filename_rm;
  for(int i=0;i<num_vib;i++){
    filename_rm.push_back("./LIB/polynomial1_func" + IntToString(i));
  };
  WR.ReadRegressionModel(filename_rm);

  WR.ReadPostProcessor("./LIB/SVD1_func");

  vector<double> xinp(input_data);

  xinp[0]=1.0;           // PWR
  if(!pwr)xinp[0]=0.0;   // BWR

    xinp[1]=pu_enrichment; // Pu enrichment [wt%]
    xinp[2]=burnup;        // Burnup [GWD/t]
    xinp[3]=h_hm_ratio;    // H/HM in mol ratio
    xinp[4]=p9;       // Pu-239 atomic ratio in TRU
    xinp[5]=p0;       // Pu-240 atomic ratio in TRU
    xinp[6]=p1;       // Pu-240 atomic ratio in TRU
    xinp[7]=p2;       // Pu-240 atomic ratio in TRU

  vector<double> tru=WR.GetRegressionResult(xinp);
  return tru;
};
#endif
