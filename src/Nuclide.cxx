#include <cstdlib>
#include "Nuclide.h"

Nuclide::Nuclide(bool simple)
{
  if(simple){
    Micxs.Init("MicroCrossSectionSimple");
  }else{
    Micxs.Init("MicroCrossSection");
  };

  grp=-1;
  temperature=300.; // Room temperature as default
};

Nuclide::~Nuclide()
{
  Micxs.AllVectorClear();
};

void Nuclide::PutGrp(int i)
{
  grp=i;
  Micxs.PutGrp(i);
};

void Nuclide::PutMatnum(int i)
{
  matnum=i;
};

void Nuclide::PutMicxs(GroupDataSet micinp)
{
  if(micinp.GetName()!="MicroCrossSection"){
    cout<<"# Error in Nuclide class.\n";
    cout<<"# You are trying to put GroupDataSet into Nuclide class.\n";
    cout<<"# The name of GroupDataSet should be MicroCrossSection.\n";
    cout<<"# Your case is "<<micinp.GetName()<<"\n";
    exit(0);
  };
  Micxs=micinp;
};

Nuclide Nuclide::Cond(int ngrp,int *bgrp,GroupData1D fl,GroupData1D cu)
{
  Nuclide ret;
  ret.PutGrp(ngrp);
  ret.PutMatnum(matnum);
  ret.PutDensity(density);
  ret.PutMicxs(Micxs.Cond(ngrp,bgrp,fl,cu));
  return ret;
};

void Nuclide::CalMicroFromF(int g,LibData &lib,vector<real> &f,bool matrix)
{
  real f_fis=f[0];
  real f_cap=f[1];
  real f_ela=f[2];
  real f_tot=f[3];
  real f_rem=f[4];
  real f_ine=f[5];

  real s_fis=lib.GetXSData().GetData1d(sigf).get_dat(g);
  real s_cap=lib.GetXSData().GetData1d(sigc).get_dat(g);
  real s_ela=lib.GetXSData().GetData1d(sigel).get_dat(g);
  real s_ine=lib.GetXSData().GetData1d(siginel).get_dat(g);
  real s_n2n=lib.GetXSData().GetData1d(sign2n).get_dat(g);

  real eff_fis=f_fis*s_fis;
  real eff_cap=f_cap*s_cap;
  real eff_ela=f_ela*s_ela;
  real eff_ine=f_ine*s_ine;

  Micxs.GetData1d(sigf).put_data(g,eff_fis);
  Micxs.GetData1d(sigc).put_data(g,eff_cap);
  Micxs.GetData1d(sigel).put_data(g,eff_ela);
  Micxs.GetData1d(siginel).put_data(g,eff_ine);

  real s_tot=s_fis+s_cap+s_ela+s_ine+s_n2n;
  real eff_sigt0=eff_fis+eff_cap+eff_ela+eff_ine+s_n2n;
  real eff_sigt1=f_tot*s_tot;
  if(f_tot<0.)eff_sigt1=-f_tot*eff_sigt0;
  //if(f_tot<0.)eff_sigt1=-f_tot*s_tot;
  if(f_tot==0.)eff_sigt1=eff_sigt0;
  Micxs.GetData1d(sigt,0).put_data(g,eff_sigt0);
  Micxs.GetData1d(sigt,1).put_data(g,eff_sigt1);
  //real corr_p1=eff_sigt1/eff_sigt0;

  // In the following, P1 correction for high-order scattering matrix is ignored.
  if(matrix){

    real f_self=f_ela;

    // ... self-shielding factor for self-scattering cross section
    /*
    real el_self=lib.GetXSData().GetData2d(sigel,0).get_dat(g,g);
    real el_rem=0.;
    for(int j=g+1;j<grp;j++){
      el_rem+=lib.GetXSData().GetData2d(sigel,0).get_dat(g,j);
    };
    f_self=(eff_ela-f_rem*el_rem)/el_self;
    */
    
    for(int l=0;l<Micxs.GetDim2d(0);l++){
      real tmp=lib.GetXSData().GetData2d(sigel,l).get_dat(g,g);
      //if(l!=0&&f_ela!=1.)tmp*=corr_p1;
      //Micxs.GetData2d(sigel,l).put_data(g,g,tmp*f_ela);  
      Micxs.GetData2d(sigel,l).put_data(g,g,tmp*f_self);      
      real ftmp=f_rem;
      //if(l!=0&&ftmp!=1.)ftmp*=corr_p1;
      for(int j=g+1;j<grp;j++){
        real tmp=lib.GetXSData().GetData2d(sigel,l).get_dat(g,j);
        Micxs.GetData2d(sigel,l).put_data(g,j,tmp*ftmp);
      };
    };
    if(s_ine>0.){
    for(int l=0;l<Micxs.GetDim2d(1);l++){
      real ftmp=f_ine;
      for(int j=g;j<grp;j++){
        real tmp=lib.GetXSData().GetData2d(siginel,l).get_dat(g,j);
        Micxs.GetData2d(siginel,l).put_data(g,j,tmp*ftmp);
      };
    };
    };
  };
};

bool Nuclide::IsFissionable()
{
  for(int i=0;i<grp;i++){
    if(Micxs.GetData1d(sigf).get_dat(i)>0.)return true;
  };
  return false;
};

void Nuclide::PrintMicroSectionTable(GroupData1D &enband)
{
  cout.setf(ios::scientific);
  cout.precision(4);

  cout<<"# Upper    Nu         Fission    Capture    (n,2n)     Chi    total    Sigel(0~5)\n";
  cout<<"# energy                                                                         \n";
  for(int g=0;g<grp;g++){
    cout<<enband.get_dat(g)<<" ";
    cout<<Micxs.GetData1d(nu).get_dat(g)<<" ";
    cout<<Micxs.GetData1d(sigf).get_dat(g)<<" ";
    cout<<Micxs.GetData1d(sigc).get_dat(g)<<" ";
    cout<<Micxs.GetData1d(sign2n).get_dat(g)<<" ";
    cout<<Micxs.GetData1d(chi).get_dat(g)<<" ";
    cout<<Micxs.GetData1d(sigt).get_dat(g)<<" ";
    cout<<Micxs.GetSigs(0).get_dat(g)<<" ";
    cout<<"\n";
  };
};

