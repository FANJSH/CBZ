#include <cstdlib>
#include "HeatConductionCalculator.h"

void HeatConductionCalculator::Initialize(int mp, int mg, int mc)
{
  //temp_boundary=300.; // [K]
  
  // ... Geometric data specification

  rad_p=0.41;
  rad_g=0.42;
  rad_c=0.475;

  /*
  mesh_p=100;
  mesh_g=10;
  mesh_c=20;
  */
  mesh_p=mp;
  mesh_g=mg;
  mesh_c=mc;

  v_pellet=PI*rad_p*rad_p; // [cm3]  
  
  // *** CartMeshInfo ***
  real *xl=new real[mesh_p+2];
  int *xm=new int[mesh_p+2];
  for(int i=0;i<mesh_p;i++){
    xl[i]=rad_p/mesh_p;
    xm[i]=1;
  };
  xl[mesh_p]=rad_g-rad_p;
  xl[mesh_p+1]=rad_c-rad_g;
  xm[mesh_p]=mesh_g;
  xm[mesh_p+1]=mesh_c;  

  real yl[]={1.};
  int ym[]={1};
  
  int *mat=new int[mesh_p+2];
  for(int i=0;i<mesh_p+2;i++){
    mat[i]=i;
  };

  cmi.PutMeshInfo(mesh_p+2,1,xm,ym,xl,yl,mat,"width");  
  cmi.PutBoundaryCondition("Reflective","Zeroflux","Reflective","Reflective");
  delete [] xl;
  delete [] xm;
  delete [] mat;  
  
  
  // ... Material data specification

  Medium med_pellet, med_gas, med_cladding;
  int group=1;

  med_pellet.PutImax(group);
  med_gas.PutImax(group);
  med_cladding.PutImax(group);  

  med_pellet.PutPL(0);
  med_gas.PutPL(0);
  med_cladding.PutPL(0);  

  real xs[1];

  // (pellet)
  real t_f=1000; // [oC]
  xs[0]=1./(11.75+0.0235*t_f);
  med_pellet.PutData(d,xs); // heat conductivity [W/cm/K] as diffusion coefficient
  xs[0]=1e-10;  
  med_pellet.PutData(siga,xs);
  med_pellet.PutData(sigt,xs);  
  xs[0]=0.0;
  med_pellet.PutData(nusigf,xs);
  xs[0]=0.;
  med_pellet.PutDataSigs(0,xs);
  xs[0]=10.97*(0.052*4.184); // density [g/cm3] x  specific heat capacity ([cal/g/K] x [J/cal])
  med_pellet.PutData(dr,xs);

  // (gas)
  xs[0]=0.3e-2;
  med_gas.PutData(d,xs);
  xs[0]=1e-10;
  med_gas.PutData(siga,xs);
  med_gas.PutData(sigt,xs);  
  xs[0]=0.0;
  med_gas.PutData(nusigf,xs);
  xs[0]=0.;
  med_gas.PutDataSigs(0,xs);
  //xs[0]=10.97*(0.052*4.184); // density [g/cm3] x  specific heat capacity ([cal/g/K] x [J/cal])
  xs[0]=0.4e-6*1.12; // density [g/cm3] x  specific heat capacity ([J/g/K])  
  // From Haraguchi-kun in 2021/5/17
  med_gas.PutData(dr,xs);

  // (cladding)
  xs[0]=12.1e-2;
  med_cladding.PutData(d,xs);
  xs[0]=1e-10;
  med_cladding.PutData(siga,xs);
  med_cladding.PutData(sigt,xs);  
  xs[0]=0.0;
  med_cladding.PutData(nusigf,xs);
  xs[0]=0.;
  med_cladding.PutDataSigs(0,xs);
  //xs[0]=10.97*(0.052*4.184); // density [g/cm3] x  specific heat capacity ([cal/g/K] x [J/cal])
  xs[0]=6.52*(0.302); // density [g/cm3] x  specific heat capacity ([J/g/K])  
  // From Haraguchi-kun in 2021/5/17  
  med_cladding.PutData(dr,xs);

  // ... PLOS instance preparation
  int mednum=mesh_p+2;  

  GeneralOption opt;
  //opt.PutEpsf(1e-6);

  test.Init(2,1,mednum);
  for(int i=0;i<mesh_p;i++){
    test.AddMedium(med_pellet);
  };
  test.AddMedium(med_gas);
  test.AddMedium(med_cladding);  
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.NoPrint();  

  totm=test.GetTotM();

  r.resize(totm,0.);
  xs_corr=false;


};

void HeatConductionCalculator::SetTemperatureBoundary(real tin)
{
  temp_boundary=tin;
  test.PutXrBoundaryCondition(temp_boundary);
};

void HeatConductionCalculator::TemperatureDistributionIteration(vector<real> &src_ex)
{
  int iteration_heat_conductivity_max=10;
  real eps_ihc=1.; // absolute residual in temperature
  
  vector<real> temp_ihc(totm,0.);  
  GroupData1D inp(1);
  
  for(int i=0;i<iteration_heat_conductivity_max;i++){
    
    test.CalCoef();
    test.SetZeroScatSrc();
    for(int i=0;i<test.GetTotM();i++){
      inp.put_data(0,src_ex[i]);
      test.PutIsotropicSourceParVolume(i,inp);
    };

    test.CalFixedSource();

    // ... heat conductivity reset
    real errmax=0.;
    for(int i=0;i<mesh_p;i++){
      real temp=test.GetMesh(i).GetFlux().get_dat(0);
      // ... default 
      real heat_conductivity_new=1./(11.75+0.0235*(temp-273.15));
      // ... optional      
      //real heat_conductivity_new=1./(11.75+0.0235*((1453+733)*0.5-273.15));
      //real heat_conductivity_new=1./(11.75+0.0235*(1453-273.15));
      //real heat_conductivity_new=1./(11.75+0.0235*(715-273.15));      
      //real heat_conductivity_new=5.066e-2;
      test.GetMedium(i).GetMacxs().GetData1d(d).put_data(0,heat_conductivity_new);
      real err=fabs(temp-temp_ihc[i]);
      if(err>errmax)errmax=err;
      temp_ihc[i]=temp;
    };

    //cout<<i<<" "<<errmax<<"\n";
    if(errmax<eps_ihc&&i!=0){
      //cout<<i<<" "<<errmax<<"\n"; exit(0);
      break;
    };

  };

};

void HeatConductionCalculator::StationaryStateCalculation(real temp_b_in, real line_power)
{
  SetTemperatureBoundary(temp_b_in);
  
  real pow=line_power/v_pellet; // [W/cm3]  
  
  vector<real> src_ex(totm);
  for(int i=0;i<totm;i++){
    real inp=0.;
    if(i<mesh_p)src_ex[i]=pow;
  };
  TemperatureDistributionIteration(src_ex);
};


void HeatConductionCalculator::RunTransient(real temp_b_in, real line_power)
{
  SetTemperatureBoundary(temp_b_in);  
  real pow2=line_power/v_pellet;

  vector<real> q_ex(totm);
  vector<real> src_ex(totm);    

  // ... Calculation of Q_ex
  for(int i=0;i<totm;i++){
    q_ex[i]=(1./theta-1.)*r[i]+test.GetMesh(i).GetFlux().get_dat(0)*(test.GetMesh(i).GetMed()->GetMacxs().GetData1d(siga).get_dat(0));
  };

  for(int i=0;i<totm;i++){
    real inp=q_ex[i];
    if(i<mesh_p)inp+=pow2;
    src_ex[i]=inp;
  };

  TemperatureDistributionIteration(src_ex);
  
  // calculation R
  for(int i=0;i<totm;i++){
    r[i]=test.GetMesh(i).GetFlux().get_dat(0)*(test.GetMesh(i).GetMed()->GetMacxs().GetData1d(siga).get_dat(0))-q_ex[i];
  };
};

void HeatConductionCalculator::Run()  
{
};

void HeatConductionCalculator::CrossSectionCorrectionForDynamicCalculation(real theta_inp, real dt)
{
  if(xs_corr){
    cout<<"# Error in HeatConductionCalculator::CrossSectionCorrectionForDynamicCalculation.\n";
    cout<<"# Correction has been already done.\n";
    exit(0);
  };

  xs_corr=true;

  theta=theta_inp;
  real factor=1./(theta*dt);
  // ... pseudo absorption adding
  for(int i=0;i<mesh_p+2;i++){
    real xs_corr=test.GetMed(i).GetMacxs().GetData1d(dr).get_dat(0)*factor;
    test.GetMed(i).GetMacxs().GetData1d(siga).put_data(0,xs_corr);
    test.GetMed(i).GetMacxs().GetData1d(sigt).put_data(0,xs_corr);
  };
};
  


void HeatConductionCalculator::ResultPrinting()
{
  cout<<"# R (at mesh-left-boundary) [cm]   temperature [K]\n";
  real x=0.;
  for(int i=0;i<totm+1;i++){
    int ii=i;
    if(ii==totm)ii=i-1;
    cout<<x<<" "<<test.GetMesh(ii).GetFlux().get_dat(0)<<"\n";    
    x+=cmi.GetFMeshL(0,i);
  };
  cout<<"\n\n";
};

vector<real> HeatConductionCalculator::GetTemperatureData()
{
  vector<real> ret(totm);
  for(int i=0;i<totm;i++){
    ret[i]=test.GetMesh(i).GetFlux().get_dat(0);    
  };
  return ret;
};

real HeatConductionCalculator::GetMeshWidth(int i)
{
  if(i<mesh_p){
    return rad_p/mesh_p;
  }else if(i<mesh_p+mesh_g){
    return (rad_g-rad_p)/mesh_g;
  }else if(i<mesh_p+mesh_g+mesh_c){
    return (rad_c-rad_g)/mesh_c;
  }else{
    cout<<"# Error in HeatConductionCalculator::GetMeshWidth.\n";
    exit(0);
  };
  return 0.;
};

real HeatConductionCalculator::GetHeatFluxAtBoundaryPerUnitHeight()
{
  real temp=test.GetMesh(totm-1).GetFlux().get_dat(0);
  real half_width=(rad_c-rad_g)/mesh_c*0.5;
  real lambda=test.GetMesh(totm-1).GetMed()->GetMacxs().GetData1d(d).get_dat(0);
  real surface_area=2.*PI*rad_c;
  return lambda*(temp-temp_boundary)/half_width*surface_area; 
};

real HeatConductionCalculator::GetVolumeAveragedTemperature()
{
  real sum1=0.;
  real sum2=0.;
  for(int i=0;i<mesh_p;i++){
    real vol=test.GetMesh(i).GetVolume();
    real temp=test.GetMesh(i).GetFlux().get_dat(0);
    sum1+=vol;
    sum2+=temp*vol;
  };
  return sum2/sum1;
};

real HeatConductionCalculator::GetCenterTemperature()
{
  return test.GetMesh(0).GetFlux().get_dat(0);
};

real HeatConductionCalculator::GetOuterTemperature()
{
  return test.GetMesh(mesh_p-1).GetFlux().get_dat(0);
};
