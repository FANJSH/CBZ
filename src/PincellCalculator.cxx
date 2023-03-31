#include <cstdlib>
#include "PincellCalculator.h"

PincellCalculator::PincellCalculator(bool hex_in)
{
  hexagonal=hex_in;
  geom_print=false;
  background_xs_printing=false;
  thermal_cutoff_energy_coolant=10.;
  thermal_cutoff_energy_fuel=10.;
  cw_correction=true;
  weight_transport_xs_fuel=1;
};

void PincellCalculator::PutMednum(int i)
{
  mednum=i;
  med.resize(mednum);
};

void PincellCalculator::PutRing(int i)
{
  ring=i;
  radius.resize(ring);
  medium_id.resize(ring+1);
};

void PincellCalculator::PutRadius(real *inp)
{
  for(int i=0;i<ring;i++){
    radius[i]=inp[i];
    if(i!=0){
      if(radius[i]<=radius[i-1]){
	cout<<"# Error in PincellCalculator::PutRadius.\n";
        cout<<"# Radius values should be in ascending order.\n";
	exit(0);
      };
    };
  };
};

void PincellCalculator::PutMediumID(int *inp)
{
  for(int i=0;i<ring;i++){
    medium_id[i]=inp[i];
    if(i!=0){
      if(medium_id[i]<medium_id[i-1]){
	cout<<"# Error in PincellCalculator::PutMediumID.\n";
        cout<<"# Medium ID should be in ascending order.\n";
	exit(0);
      };
    };
  };
  medium_id[ring]=inp[ring-1];
};

void PincellCalculator::PutBoundaryCondition(enum BCondition bci_ssc,enum BCondition bci_flx)
{
  bc_ssc=bci_ssc;
  bc_flx=bci_flx;
};

void PincellCalculator::PutMedium(vector<Medium> &minp)
{
  for(int i=0;i<mednum;i++){
    med[i]=minp[i];
  };
};

void PincellCalculator::SetTrajectorySet()
{
  real *rr=new real[ring];
  for(int i=0;i<ring;i++){
    rr[i]=radius[ring-1-i];
  };

  vector<real> macro_radius(mednum-1,0.);
  int mesh_fuel=0; // number of mesh in fuel region
  for(int i=0;i<ring;i++){
    if(medium_id[i]==0)mesh_fuel++;
    for(int j=0;j<mednum-1;j++){
      if(medium_id[i]==j+1&&macro_radius[j]==0.)macro_radius[j]=radius[i-1];
    };
  };

  /*
  real fuel_r=macro_radius[0];
  real fuel_vol=fuel_r*fuel_r*PI;
  real r=fuel_r;
  */

  // TrajectorySet for self-shielding calculation
  IrregularGeometryInformation igi;
  GeomPolygon pol;
  if(hexagonal){
    real bar=pin_pitch/sqrt(3.);
    pol.PutHexagon(0.,0.,bar);
  }else{
    pol.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  };
  pol.PutRegionID(mednum-1);
  igi.AddGeom(pol);

  for(int i=mednum-2;i>=0;i--){
    GeomCircle cir(0.,0.,macro_radius[i]);
    cir.PutRegionID(i);
    igi.AddGeom(cir);
  };

  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,8,0.02,45.); 

  // TrajectorySet for flux distribution calculation
  IrregularGeometryInformation igi_f;
  GeomPolygon pol2;
  if(hexagonal){
    real bar=pin_pitch/sqrt(3.);
    pol2.PutHexagon(0.,0.,bar);
  }else{
    pol2.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  };
  pol2.PutRegionID(ring);
  igi_f.AddGeom(pol2);
  int *rid=new int[ring];
  for(int i=0;i<ring;i++){rid[i]=ring-1-i;};
  igi_f.AddCircleRing(ring,rr,rid);

  if(geom_print)igi_f.WriteGnuplotFile(0.001);

  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi_f,8,0.02,45.);
  //sys_f.CalTrajectory(igi_f,64,0.02,360.);  

  delete [] rr;
  delete [] rid;
};

void PincellCalculator::SelfShieldingCalculation(XSLibrary &xslib)
{
  OnePointCalculator opc;
  if(!cw_correction)opc.CWCorrectionOff();
  
  for(int i=0;i<mednum;i++){
    opc.GiveInfiniteDillutionCrossSection(med[i],xslib);
  };

  //opc.CalThermalScatteringMatrix(med[mednum-1],xslib, 3.93); // 3.93 : thermal cut-off energy
  opc.CalThermalScatteringMatrix(med[mednum-1],xslib, thermal_cutoff_energy_coolant);// thermal cut-off energy
  med[mednum-1].CalMacroFromMicro();
  med[mednum-1].CalSigtr(0);

  GroupData1D c(group); // Dancoff correction
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};

  SelfShieldingCalculator ssc;
  if(background_xs_printing)ssc.BackgroundXSPrinting();
  if(!cw_correction)ssc.CWCorrectionOff();

  ssc.DancoffMethod(xslib,sys,med);
  //ssc.GetDancoff(0).show_self();

  // Thermal scattering matrices are over-written
  //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
  opc.CalThermalScatteringMatrix(med[0],xslib,thermal_cutoff_energy_fuel); // thermal cut-off energy  
  med[0].CalMacroFromMicro();
  med[0].CalSigtr(weight_transport_xs_fuel);

  opc.CalFissionSpectrumMatrix(med[0],xslib);
};

Medium PincellCalculator::GetHomogenizedCrossSection(XSLibrary &xslib, bool ssf_cal)
{
  int tmesh=0;
  int tmeshid[]={0};
  return GetHomogenizedCrossSection(xslib,tmesh,tmeshid,ssf_cal);
};

Medium PincellCalculator::GetHomogenizedCrossSection(XSLibrary &xslib,int tmesh,int *tmeshid,bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt;

  // +++ Eigenvalue calculation
  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigt);

  lat.PutPij();
  real keff=lat.CalIgenPij();

  lat.PutFluxAsCurrent();
  //lat.PutPijr();

  Medium homoxs;
  if(tmesh==0){
    homoxs=lat.HomogenizeAll(false);
  }else{
    homoxs=lat.Homogenize(tmesh,tmeshid,false);
  };
  
  return homoxs;
};

Medium PincellCalculator::GetHomogenizedCrossSectionBilinear(XSLibrary &xslib, bool ssf_cal)
{
  int tmesh=0;
  int tmeshid[]={0};
  return GetHomogenizedCrossSectionBilinear(xslib,tmesh,tmeshid,ssf_cal);
};

Medium PincellCalculator::GetHomogenizedCrossSectionBilinear(XSLibrary &xslib,int tmesh,int *tmeshid,bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation
  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigt);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  PJISystem lata(group,mednum);
  lata.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lata.AddMedium(med[i]);
  };
  lata.PutRegMed(medium_id);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigt);
  lata.PutPij();
  real keffa=lata.CalIgenPij();

  int ttm=lat.GetTotM();
  for(int i=0;i<ttm;i++){
    GroupData1D flx=lat.GetMesh(i).GetFlux();
    GroupData1D flxa=lata.GetMesh(i).GetFlux();
    flx=flx.mult(flxa);
    lat.GetMesh(i).GetFlux().copy(flx);
  };

  lat.PutFluxAsCurrent();
  //lat.PutPijr();

  Medium homoxs;
  if(tmesh==0){
    homoxs=lat.HomogenizeAll(false);
  }else{
    homoxs=lat.Homogenize(tmesh,tmeshid,false);
  };
  
  return homoxs;
};

SensitivityData PincellCalculator::CalKinfSensitivity(XSLibrary &xslib,int nucnum,int *nucid,bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opta);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  PJISystem lat2(group,mednum);
  lat2.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat2.AddMedium(med[i]);
  };
  lat2.PutRegMed(medium_id);
  lat2.PutGeneralOption(opt);
  lat2.PutSigmaCol(sigtr);
  lat2.PutPij();
  real k2=lat2.CalIgenPij();

  SensitivityData sens=lat.CalSensitivityNew(&lat2,keff,nucnum,nucid);
  return sens;
};

SensitivityData PincellCalculator::CalKinfSensitivityByMEC(XSLibrary &xslib,int nucnum,int *nucid,bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      med[i].CalSigt();
      lat.AddMedium(med[i]);
    };
    lat.PutRelationRegionMedium(medium_id);
    lat.PutGeneralOption(opta);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    lat.PutWriteFlux();
    lat.SetArray();
    real keff=lat.CalIgen();

    MECSystem lat2(group,mednum);
    lat2.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      lat2.AddMedium(med[i]);
    };
    lat2.PutRelationRegionMedium(medium_id);
    lat2.PutGeneralOption(opt);
    lat2.PutPL(0);
    lat2.NoCMRAcceleration();
    //lat2.NoTransportApprox();
    lat2.PutWriteFlux();
    lat2.SetArray();
    real k2=lat2.CalIgen();

  SensitivityData sens=lat.CalSensitivityNew(&lat2,keff,nucnum,nucid);
  return sens;
};

SensitivityData PincellCalculator::CalRRRSensitivity
(int nume_nuc,int* nume_id, enum xstype* nume_xs,
 int denom_nuc,int* denom_id, enum xstype* denom_xs, 
 int ss_nuc, int* mat_id, XSLibrary &xslib)
{
  SetTrajectorySet();

  GeneralOption opt,opta;
  opta.PutAdjointCal();
    
  SelfShieldingCalculation(xslib);

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  // (adjoint system)
  PJISystem lata(group,mednum);
  lata.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lata.AddMedium(med[i]);
  };
  lata.PutRegMed(medium_id);
  lata.PutGeneralOption(opta);
  lata.PutSigmaCol(sigtr);
  lata.PutPij();
  lata.PutPL(0);  // !!! This method is important in the NON-eigenvalue calculation mode !!!

  int ttm=lat.GetTotM();
  bool *on_mesh=new bool[ttm];
  for(int i=0;i<ttm;i++){
    on_mesh[i]=false;
    if(medium_id[i]==0)on_mesh[i]=true;
  };

  SensitivityData sens=lata.CalSensitivityRRR(lat,nume_nuc,nume_id,nume_xs,denom_nuc,denom_id,denom_xs,on_mesh,ss_nuc,mat_id,keff);

  delete [] on_mesh;

  return sens;
};

SensitivityData PincellCalculator::CalRRRSensitivityByMEC
(int nume_nuc,int* nume_id, enum xstype* nume_xs,
 int denom_nuc,int* denom_id, enum xstype* denom_xs, 
 int ss_nuc, int* mat_id, XSLibrary &xslib)
{
  // Transport correction is NOT adopted
  // because no convergence is obtained.
  //
  //

  SetTrajectorySet();

  GeneralOption opt,opta;
  opta.PutAdjointCal();
    
  SelfShieldingCalculation(xslib);

  // +++ Eigenvalue calculation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      med[i].CalSigt();
      lat.AddMedium(med[i]);
    };
    lat.PutRelationRegionMedium(medium_id);
    lat.PutGeneralOption(opt);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.NoTransportApprox();
    lat.PutWriteFlux();
    real keff=lat.CalIgen();

    MECSystem lata(group,mednum);
    lata.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      lata.AddMedium(med[i]);
    };
    lata.PutRelationRegionMedium(medium_id);
    lata.PutGeneralOption(opta);
    lata.PutPL(0);
    lata.NoCMRAcceleration();
    lata.NoTransportApprox();
    lata.PutWriteFlux();
    lata.SetArray();

  int ttm=lat.GetTotM();
  bool *on_mesh=new bool[ttm];
  for(int i=0;i<ttm;i++){
    on_mesh[i]=false;
    if(medium_id[i]==0)on_mesh[i]=true;
  };

  SensitivityData sens=lata.CalSensitivityRRR(lat,nume_nuc,nume_id,nume_xs,denom_nuc,denom_id,denom_xs,on_mesh,ss_nuc,mat_id,keff);

  delete [] on_mesh;

  return sens;
};

SensitivityData PincellCalculator::CalRRRSensitivityDirect
(int nume_nuc,int* nume_id, enum xstype* nume_xs,
 int denom_nuc,int* denom_id, enum xstype* denom_xs,
 int matid,enum xstype xst,XSLibrary &xslib)
{
  real factor=0.01;

  int mtid=0;
  if(xst==sigc){
    mtid=102;
  }else if(xst==sigf){
    mtid=18;
  }else{
    cout<<"# Error in PincellCalculator::CalRRRSensitivityDirect.\n";
    cout<<"# The reaction cannot be coded yet.\n";
    exit(0);
  };

  SetTrajectorySet();

  GeneralOption opt;

  vector<real> val(group+1);

  for(int ii=0;ii<group+1;ii++){

    real dxs;
    if(ii!=group){
      dxs=xslib.GetLibData(matid).GetXSData().GetData1d(xst).get_dat(ii)*factor;
      xslib.GetLibData(matid).GetXSData().GetData1d(xst).add_data(ii,dxs);
      xslib.GetLibData(matid).GetXSData().GetData1d(sigt).add_data(ii,dxs);
    };
    
    SelfShieldingCalculation(xslib);

    // +++ Eigenvalue calculation

    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      med[i].CalSigt();
      lat.AddMedium(med[i]);
    };
    lat.PutRelationRegionMedium(medium_id);
    lat.PutGeneralOption(opt);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.NoTransportApprox();
    lat.PutWriteFlux();
    lat.SetArray();
    real keff=lat.CalIgen();

    /*
    PJISystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      lat.AddMedium(med[i]);
    };
    lat.PutRegMed(medium_id);
    lat.PutGeneralOption(opt);
    lat.PutSigmaCol(sigtr);
    lat.PutPij();
    real keff=lat.CalIgenPij();
    */

    int ttm=lat.GetTotM();
    bool *on_mesh=new bool[ttm];
    for(int i=0;i<ttm;i++){
      on_mesh[i]=false;
      if(medium_id[i]==0)on_mesh[i]=true;
    };
    real nume=lat.CalMacroscopicReactionRate(nume_nuc,nume_id,nume_xs,on_mesh);
    real denom=lat.CalMacroscopicReactionRate(denom_nuc,denom_id,denom_xs,on_mesh);
    delete [] on_mesh;

    val[ii]=nume/denom;

    if(ii!=107){
      xslib.GetLibData(matid).GetXSData().GetData1d(xst).add_data(ii,-dxs);
      xslib.GetLibData(matid).GetXSData().GetData1d(sigt).add_data(ii,-dxs);
    };

  };

  SensitivityData sens;
  sens.PutName("dummy","dummy","dummy");
  sens.PutValue(val[group]);
  sens.PutGroup(group);
  sens.GetEnband().copy(xslib.GetEnband());
  GroupData1D sns1d(group);
  for(int g=0;g<group;g++){
    real dd=(val[g]-val[group])/(val[group]*factor);
    sns1d.put_data(g,dd);
  };
  sens.PutSensitivity1D(matid,mtid,sns1d);

  return sens;
};

void PincellCalculator::FissionSpectrumVectorCalculation(XSLibrary &xslib)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt;

  // +++ Eigenvalue calculation
  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  med[0].GetFlux().copy(lat.GetIntegratedFlux(0));
  med[0].CalMacroFromMicro(); // Chi-vector calculation
};

void PincellCalculator::ShowNeutronFluxInFuel(XSLibrary &xslib, bool adj, bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt;
  if(adj)opt.PutAdjointCal();

  // +++ Eigenvalue calculation
  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  cout<<"#\n# Eup       Flux";
  if(!adj)cout<<"        Flux per lethargy";
  cout<<"\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<group;i++){
    real e0=med[0].GetEnband().get_dat(i);
    real e1=med[0].GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    cout<<e0<<" ";
    real tmp=lat.GetIntegratedFlux(0).get_dat(i);
    cout<<tmp<<" ";
    if(!adj){
      cout<<tmp/letwid;
    };
    cout<<"\n";
    //cout<<e0<<" "<<lat.GetIntegratedFlux(0).get_dat(i)*factor<<"\n";
  };

};

GroupData1D PincellCalculator::CalNeutronFluxInFuel(XSLibrary &xslib, bool adj, bool ssf_cal)
{
  SetTrajectorySet();
  if(ssf_cal)SelfShieldingCalculation(xslib);

  GeneralOption opt;
  if(adj)opt.PutAdjointCal();

  // +++ Eigenvalue calculation
  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();
  
  return lat.GetIntegratedFlux(0);

};

real PincellCalculator::EigenvalueCalculation(XSLibrary &xslib,bool transport)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);
  GeneralOption opt;

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opt);
  if(!transport){
    lat.PutSigmaCol(sigt);
  }else{
    lat.PutSigmaCol(sigtr);
  };
  lat.PutPij();
  real keff=lat.CalIgenPij(true); // fission spectrum matrix is used.

  return keff;
};

real PincellCalculator::EigenvalueCalculationByMEC(XSLibrary &xslib)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt;
  //opt.PutEpsf(1e-6);
  //opt.PutEpsk(1e-7);

  // +++ Eigenvalue calculation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);
    for(int i=0;i<mednum;i++){
      med[i].CalSigt();
      lat.AddMedium(med[i]);
    };
    lat.PutRelationRegionMedium(medium_id);
    lat.PutGeneralOption(opt);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    //lat.PutWriteFlux();
    real keff=lat.CalIgen();

  return keff;
};

void PincellCalculator::ReactivityCalculation(XSLibrary &xslib, vector<Medium> &minp_pert)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opta);
  //lat.PutSigmaCol(sigt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  PutMedium(minp_pert);
  SelfShieldingCalculation(xslib);

  PJISystem lat2(group,mednum);
  lat2.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat2.AddMedium(med[i]);
  };
  lat2.PutRegMed(medium_id);
  lat2.PutGeneralOption(opt);
  //lat.PutSigmaCol(sigt);
  lat2.PutSigmaCol(sigtr);
  lat2.PutPij();
  real keff2=lat2.CalIgenPij();

  lat.CalReactivity(&lat2,keff,keff2);
};

GData PincellCalculator::ReactivityCalculation(XSLibrary &xslib, XSLibrary &xslib2, vector<Medium> &minp_pert)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opta);
  //lat.PutSigmaCol(sigt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  PutMedium(minp_pert);
  SelfShieldingCalculation(xslib2);

  PJISystem lat2(group,mednum);
  lat2.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat2.AddMedium(med[i]);
  };
  lat2.PutRegMed(medium_id);
  lat2.PutGeneralOption(opt);
  //lat.PutSigmaCol(sigt);
  lat2.PutSigmaCol(sigtr);
  lat2.PutPij();
  real keff2=lat2.CalIgenPij();

  GData ret;
  ret=lat.CalReactivityGData(&lat2,keff,keff2);
  return ret;
};

void PincellCalculator::ReactivityCalculationByMEC(XSLibrary &xslib, vector<Medium> &minp_pert)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation
  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    med[i].CalSigt();
    lat.AddMedium(med[i]);
  };
  lat.PutRelationRegionMedium(medium_id);
  lat.PutGeneralOption(opta);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  //lat.NoTransportApprox();
  lat.PutWriteFlux();
  real keff=lat.CalIgen();
  
  PutMedium(minp_pert);
  SelfShieldingCalculation(xslib);

  MECSystem lat2(group,mednum);
  lat2.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    med[i].CalSigt();
    lat2.AddMedium(med[i]);
  };
  lat2.PutRelationRegionMedium(medium_id);
  lat2.PutGeneralOption(opt);
  lat2.PutPL(0);
  lat2.PutWriteFlux();
  //lat2.NoTransportApprox();  
  lat2.NoCMRAcceleration();
  real keff2=lat2.CalIgen();

  lat.CalReactivity(&lat2,keff,keff2);
};

void PincellCalculator::ReactivityCalculation(XSLibrary &xslib, XSLibrary &xslib2)
{
  SetTrajectorySet();
  SelfShieldingCalculation(xslib);

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Eigenvalue calculation

  PJISystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat.AddMedium(med[i]);
  };
  lat.PutRegMed(medium_id);
  lat.PutGeneralOption(opta);
  //lat.PutSigmaCol(sigt);
  lat.PutSigmaCol(sigtr);
  lat.PutPij();
  real keff=lat.CalIgenPij();

  SelfShieldingCalculation(xslib2);

  PJISystem lat2(group,mednum);
  lat2.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum;i++){
    lat2.AddMedium(med[i]);
  };
  lat2.PutRegMed(medium_id);
  lat2.PutGeneralOption(opt);
  //lat.PutSigmaCol(sigt);
  lat2.PutSigmaCol(sigtr);
  lat2.PutPij();
  real keff2=lat2.CalIgenPij();

  lat.CalReactivity(&lat2,keff,keff2);
};

void PincellCalculator::ReadFile(string mdir, string ss)
{
  mdir.append(ss);

  ifstream fin;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(1);
  };

  int tmp;
  real ftmp;

  fin>>tmp; // # of mediums
  PutMednum(tmp);

  fin>>ftmp;
  PutPinPitch(ftmp); // Pitch

  fin>>tmp;
  PutRing(tmp); // # of rings

  real *radin=new real[ring];
  for(int i=0;i<ring;i++){
    fin>>radin[i];
  };
  PutRadius(radin);
  delete [] radin;

  int *medidinp=new int[ring];
  for(int i=0;i<ring;i++){
    fin>>medidinp[i];
  };
  PutMediumID(medidinp);
  delete [] medidinp;

  vector<Medium> medin(mednum);
  for(int i=0;i<mednum;i++){
    medin[i].PutImax(group);
    medin[i].PutPL(1);
    int nucnum;
    fin>>nucnum;
    int *idin=new int[nucnum];
    real *denin=new real[nucnum];
    for(int i=0;i<nucnum;i++){
      fin>>idin[i];
    };
    for(int i=0;i<nucnum;i++){
      fin>>denin[i];
    };
    medin[i].PutNuclide(nucnum,idin,denin);
    delete [] idin;
    delete [] denin;
    fin>>ftmp; // temperature
    medin[i].PutTemperatureForAllNuclide(ftmp);
  };

  PutMedium(medin);

};

void PincellCalculator::CurrentWeightTotalCalculation(XSLibrary &xslib,string filename,vector<real> &err)
{
  /*
  int group=707;
  int bgroup=107;
  int *bgrp=new int[bgroup];
  int tg=-1;
  int delta_tg=1;
  for(int i=0;i<107;i++){
    if(i==37)delta_tg=25;
    if(i==62)delta_tg=1;
    tg+=delta_tg;
    bgrp[i]=tg;
  };
  */

  int group=907;
  int bgroup=91;
  int *bgrp=new int[bgroup];
  int tg=-1;
  int delta_tg=1;
  for(int i=0;i<91;i++){
    if(i==56)delta_tg=25;
    if(i==90)delta_tg=1;
    tg+=delta_tg;
    bgrp[i]=tg;
  };

  if(xslib.GetGroup()!=group){
    cout<<"# Error in PincellCalculator::CurrentWeightTotalCalculation.\n";
    cout<<"# The number of energy groups in the library should be "<<group<<".\n";
    cout<<"# The number of energy groups used is "<<xslib.GetGroup()<<"\n";
    exit(0);
  };

  if(mednum!=3){
    cout<<"# Error in PincellCalculator::CurrentWeightTotalCalculation.\n";
    cout<<"# The number of mediums should be 3.\n";
    cout<<"# The number of mediums is "<<mednum<<"\n";
    exit(0);
  };

  int reg=3;
  real hpitch=pin_pitch*0.5;

  vector<real> macro_radius(mednum-1,0.);
  int mesh_fuel=0; // number of mesh in fuel region
  for(int i=0;i<ring;i++){
    if(medium_id[i]==0)mesh_fuel++;
    for(int j=0;j<mednum-1;j++){
      if(medium_id[i]==j+1&&macro_radius[j]==0.)macro_radius[j]=radius[i-1];
    };
  };
  real fuel_r=macro_radius[0];
  real clad_r=macro_radius[1];

  // IrregularGeometryInformation
  IrregularGeometryInformation igi;
  GeomPolygon pol;
  pol.PutRectangular(0.,0.,hpitch,hpitch);
  pol.PutRegionID(2);
  igi.AddGeom(pol);
  GeomCircle cir1(0.,0.,clad_r);
  cir1.PutRegionID(1);
  igi.AddGeom(cir1);
  GeomCircle cir2(0.,0.,fuel_r);
  cir2.PutRegionID(0);
  igi.AddGeom(cir2);

  // TrajectorySet
  TrajectorySet sys;
  sys.PutBoundaryCondition(Periodic);
  sys.CalTrajectory(igi,15,0.01,45.);

  OnePointCalculator opc;
  opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[1],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[2],xslib);

  SelfShieldingCalculator ssc;
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[1],med[2],false); // [false] is for cladding SSC calculation.

  // +++ finemesh calculation
  IrregularGeometryInformation ginp;

  GeomPolygon back;
  back.PutRectangular(0.,0.,hpitch,hpitch);
  back.PutRegionID(0);
  ginp.AddGeom(back);

  real line=hpitch*0.5;
  real x=line;
  real y=line;

  int quaddd[]={2,1,3,4};

  real xtmp=-hpitch*0.5;
  real ytmp=hpitch*0.5;
  // water
  int meshid=0;
  for(int i=0;i<4;i++){
    GeomPolygon pol;
    pol.PutRectangular(xtmp,ytmp,hpitch*0.5,hpitch*0.5);
    pol.PutRegionID(meshid++);
    ginp.AddGeom(pol);
    if(i==0||i==2)xtmp+=hpitch;
    if(i==1){
      xtmp-=hpitch;
      ytmp-=hpitch;
    };
  };
  // clading
  for(int i=0;i<4;i++){
    GeomDividedCircle cir;
    cir.PutQuarterCircle(0.,0.,clad_r,quaddd[i]);
    cir.PutRegionID(meshid++);
    ginp.AddGeom(cir);
  };
  // fuel
  for(int i=0;i<4;i++){
    GeomDividedCircle cir;
    cir.PutQuarterCircle(0.,0.,fuel_r,quaddd[i]);
    cir.PutRegionID(meshid++);
    ginp.AddGeom(cir);
  };

  //ginp.WriteGnuplotFile(0.01);

  TrajectorySet test;
  test.PutBoundaryCondition(Periodic);
  test.CalTrajectory(ginp,8,0.01,360.); 

  // +++ MEC_system
  GeneralOption opt;

  MECSystem lat(group,3);
  lat.PutTrajectorySet(&test);
  for(int i=0;i<3;i++){
    lat.AddMedium(med[i]);
  };
  int region_medium[]={
    2,2,2,2,    1,1,1,1,    0,0,0,0,
  };
  lat.PutRelationRegionMedium(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  lat.NoTransportApproximation(); // P0-XS
  real k_ref=lat.CalIgen("none"); 


  vector<real> k_cond(2);
  for(int jj=0;jj<2;jj++){

  vector<Medium> bmed(3);
  int regid[]={8,4,0};
  for(int ii=0;ii<3;ii++){
    GroupData1D current=lat.GetMesh(regid[ii]).GetFlux(1);
    GroupData1D flux=lat.GetMesh(regid[ii]).GetFlux(0);
    if(jj==0)bmed[ii]=med[ii].Cond(bgroup,bgrp,flux,current,false);
    if(jj==1)bmed[ii]=med[ii].Cond(bgroup,bgrp,flux,flux,false);
    for(int i=0;i<bgroup;i++){
      real st0=bmed[ii].GetMacxs().GetData1d(sigt).get_dat(i);
      real st1=bmed[ii].GetMacxs().GetData1d(sigt,1).get_dat(i);
      bmed[ii].GetMacxs().GetData2d(sigs,0).add_data(i,i,st1-st0);
    };
    bmed[ii].GetMacxs().GetData1d(sigt).copy(bmed[ii].GetMacxs().GetData1d(sigt,1));
  };

  MECSystem blat(bgroup,3);
  blat.PutTrajectorySet(&test);
  blat.AddMedium(bmed[0]);
  blat.AddMedium(bmed[1]);
  blat.AddMedium(bmed[2]);
  blat.PutRelationRegionMedium(region_medium);
  blat.PutGeneralOption(opt);
  blat.PutPL(0);
  blat.NoCMRAcceleration();
  blat.NoTransportApproximation();
  k_cond[jj]=blat.CalIgen("none");

  };

  ofstream fout;
  string mdir="./"+filename;
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Error in PincellCalculator::CalWeightTotalCalculation.\n";
    cout<<"# Failed to open the file : "<<mdir<<"\n";
    exit(0);
  };


  GroupData1D fl0=lat.GetMesh(8).GetFlux();
  GroupData1D fl1=lat.GetMesh(8).GetFlux(1);
  GroupData1D con(group);
  for(int i=0;i<group;i++){
    con.put_data(i,1.);
  };

  GroupData1D flux=lat.GetMesh(0).GetFlux(0);
  Medium mtmp=med[0].Cond(bgroup,bgrp,flux,flux,false);
  GroupData1D st0mac=med[0].GetMacxs().GetData1d(sigt).Cond(fl0,bgroup,bgrp);
  GroupData1D st1mac=med[0].GetMacxs().GetData1d(sigt).Cond(fl1,bgroup,bgrp);
  fout.setf(ios::showpoint);
  fout.precision(5);
  /*
  cout<<"# P1-total/P0-total in macroscopic cross section.\n";
  for(int g=0;g<70;g++){
    cout<<mtmp.GetEnband().get_dat(g)<<" ";
    cout<<st1mac.get_dat(g)/st0mac.get_dat(g)<<"\n";
  };
  */

  fout<<"# P1-total/P0-total in uranium-238 cross section.\n";
  int matid_cw=922380;
  //int matid_cw=942420;
  GroupData1D st0=med[0].GetNuclide(matid_cw).GetMicxs().GetData1d(sigt).Cond(fl0,bgroup,bgrp);
  GroupData1D st1=med[0].GetNuclide(matid_cw).GetMicxs().GetData1d(sigt).Cond(fl1,bgroup,bgrp);
  GroupData1D stcon=med[0].GetNuclide(matid_cw).GetMicxs().GetData1d(sigt).Cond(con,bgroup,bgrp);
  for(int g=0;g<bgroup;g++){
    fout<<mtmp.GetEnband().get_dat(g)<<" ";
    fout<<st1.get_dat(g)/st0.get_dat(g)<<"\n";
    //fout<<st1.get_dat(g)/stcon.get_dat(g)<<"\n";
  };


  fout<<"# K_eff\n";
  fout.setf(ios::showpoint);
  fout.precision(6);
  fout<<"# Reference : "<<k_ref<<"\n";
  fout<<"# Current-weighted : "<<k_cond[0]<<"\n";
  fout<<"#  (Error)         : "<<(k_cond[0]-k_ref)/k_ref<<"\n";
  fout<<"# Flux-weighted    : "<<k_cond[1]<<"\n";
  fout<<"#  (Error)         : "<<(k_cond[1]-k_ref)/k_ref<<"\n";

  fout.close();

  err[0]=(k_cond[0]-k_ref)/k_ref;
  err[1]=(k_cond[1]-k_ref)/k_ref;

  delete [] bgrp;
};
