#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "FRDesignTool.h"

using namespace std;

// +++ FRDTComposition +++

FRDTComposition::FRDTComposition()
{
  nucnum=0;
};

void FRDTComposition::PutNucnum(int inp)
{
  nucnum=inp;
  nucid.resize(nucnum);
  density.resize(nucnum,0.);
};

void FRDTComposition::ShowSelf(real factor)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<nucnum;i++){
    cout<<nucid[i]<<"  "<<density[i]*factor<<"\n";
  };
};

real FRDTComposition::GetDensityFromID(int matid)
{
  for(int i=0;i<nucnum;i++){
    if(nucid[i]==matid)return density[i];
  };

  cout<<"# Error in FRDTComposition::GetDensityFromID.\n";
  cout<<"# Materiai ID "<<matid<<" cannot be found.\n";
  exit(0);
};

void FRDTComposition::Multiply(real factor)
{
  for(int i=0;i<nucnum;i++){
    density[i]*=factor;
  };
};

void FRDTComposition::ShowSelf(MATIDTranslator &midt)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<nucnum;i++){
    cout<<midt.Name(nucid[i])<<"  "<<density[i]<<"\n";
  };
};

void FRDTComposition::AddData(int id,real den)
{
  for(int i=0;i<nucnum;i++){
    if(nucid[i]==id){
      density[i]+=den;
      return;
    };
  };
  nucnum++;
  nucid.push_back(id);
  density.push_back(den);
};

void FRDTComposition::AddData(int inpnum, int* id,real* den)
{
  for(int i=0;i<inpnum;i++){
    AddData(id[i],den[i]);
  };
};

void FRDTComposition::PutNumberDensity(Medium &med)
{
  for(int i=0;i<nucnum;i++){
    if(nucid[i]>900000){
      if(med.ExistNuclide(nucid[i])){
	med.GetNuclide(nucid[i]).PutDensity(density[i]);
      }else{
        Nuclide tmp;
        tmp.PutMatnum(nucid[i]);
        tmp.PutDensity(density[i]);
        med.AddNuclide(tmp);
      };
    };
  };

  for(int i=0;i<nucnum;i++){
    if(nucid[i]<=900000){
      if(med.ExistNuclide(nucid[i])){
	med.GetNuclide(nucid[i]).PutDensity(density[i]);
      }else{
        Nuclide tmp;
        tmp.PutMatnum(nucid[i]);
        tmp.PutDensity(density[i]);
        med.AddNuclide(tmp);
      };
    };
  };
};

FRDTComposition FRDTComposition::operator*(real a)
{
  FRDTComposition ret;
  for(int i=0;i<nucnum;i++){
    ret.AddData(nucid[i],density[i]*a);
  };
  return ret;
};

FRDTComposition FRDTComposition::operator+(const FRDTComposition &second)
{
  FRDTComposition ret;
  for(int i=0;i<nucnum;i++){
    ret.AddData(nucid[i],density[i]);
  };
  for(int i=0;i<second.nucnum;i++){
    ret.AddData(second.nucid[i],second.density[i]);
  };
  return ret;
};

void FRDTComposition::PutOxide(real density_oxide, real num_comp, real num_oxigen, int num_iso_comp, int *nucid_comp, vector<real> &wgt_ratio_comp_vec, AtomicMassData &amd)
{
  real *wgt_ratio_comp=new real[num_iso_comp];
  for(int i=0;i<num_iso_comp;i++){
    wgt_ratio_comp[i]=wgt_ratio_comp_vec[i];
  };

  PutOxide(density_oxide, num_comp, num_oxigen, num_iso_comp, nucid_comp, wgt_ratio_comp, amd);

  delete [] wgt_ratio_comp;
};

void FRDTComposition::PutOxide(real density_oxide, real num_comp, real num_oxigen, int num_iso_comp, int *nucid_comp, real *wgt_ratio_comp, AtomicMassData &amd)
{
  // If [wgt_ratio_comp] is negative, number density ratio is transferred.

  /*
  for(int i=0;i<num_iso_comp;i++){
    cout<<nucid_comp[i]<<" "<<wgt_ratio_comp[i]<<"\n";
  };
  */


  for(int i=1;i<num_iso_comp;i++){
    if(wgt_ratio_comp[0]*wgt_ratio_comp[i]<0.){
      cout<<"# Error in FRDTComposition::PutOxide.\n";
      cout<<"# Inconsistent isotopic ratio input is detected.\n";
      exit(0);
    };
  };

  vector<real> atom_ratio_comp(num_iso_comp);

  real total_rat_comp=0.;
  real total_wgt_comp=0.;
  for(int i=0;i<num_iso_comp;i++){
    real wgt=amd.GetData(nucid_comp[i]);
    if(wgt<=0.){
      cout<<"# Error in FRDTComposition::PutCompoundInformation.\n";
      cout<<"# Atomic mass data of nuclide "<<nucid_comp[i]<<" is not defined.\n";
      exit(0);
    };
    if(wgt_ratio_comp[i]<0.){
      atom_ratio_comp[i]=-wgt_ratio_comp[i];
    }else{
      atom_ratio_comp[i]=wgt_ratio_comp[i]/wgt;
    };
    total_rat_comp+=atom_ratio_comp[i];
    total_wgt_comp+=wgt*atom_ratio_comp[i];
  };
  total_wgt_comp/=total_rat_comp;

  real total_wgt=num_oxigen*amd.GetData(80160)+num_comp*total_wgt_comp;

  real mol=density_oxide/total_wgt; // [mol/cm3]
  real nd=mol*avo;

  for(int i=0;i<num_iso_comp;i++){
    AddData(nucid_comp[i], atom_ratio_comp[i]/total_rat_comp*nd*num_comp);
  };
  AddData(80160, nd*num_oxigen);
};

// +++ FRDTFuelComposition +++

FRDTFuelComposition::FRDTFuelComposition(bool special)
{
  if(!special){

    /*
  nucnum=3+1+5+3+5+1;
  nucnum_u=3;
  nucnum_hm=3+1+5+3+5;

  nucid.resize(nucnum);
  density.resize(nucnum,0.);
  int nucid_initial[]={
    922350,922360,922380,
    932370,
    942380,942390,942400,942410,942420,
    952410,952421,952430,
    962420,962430,962440,962450,962460,
    80160,
  };
    */

  nucnum=3+2+5+4+5+1;
  nucnum_u=3;
  nucnum_hm=3+2+5+4+5;

  nucid.resize(nucnum);
  density.resize(nucnum,0.);
  int nucid_initial[]={
    922350,922360,922380,
    932370,932390,
    942380,942390,942400,942410,942420,
    952410,952420,952421,952430,
    962420,962430,962440,962450,962460,
    80160,
  };
    
  for(int i=0;i<nucnum;i++){
    nucid[i]=nucid_initial[i];
  };

  }else{


  nucnum_u=2;
  nucnum_hm=2+60;
  nucnum=nucnum_hm+1;

  nucid.resize(nucnum);
  density.resize(nucnum,0.);
  int nucid_initial[]={
    922350,922380,
812080, 812090, 822080, 822090, 822120, 
832090, 832120, 832130, 832150, 842120, 
842130, 842160, 852170, 852190, 862170, 
862190, 862200, 862220, 872210, 872230, 
882230, 882240, 882250, 882260, 892250, 
892270, 902270, 902280, 902290, 902300, 
902310, 902320, 912310, 922320, 922330, 
922340, 922360, 922370, 932360, 932370, 
932390, 942360, 942380, 942390, 942400, 
942410, 942420, 952410, 952420, 952421, 
952430, 962420, 962430, 962440, 962450, 
962460, 962470, 962480, 962500, 972490, 
    80160,
  };

  for(int i=0;i<nucnum;i++){
    nucid[i]=nucid_initial[i];
  };

  };



  pu_fissile_enrichment=0.;
  pu_enrichment=0.;

  u5_enrichment=0.2;
  o_m=1.97;
  pellet_theoretical_density=85.;
};

void FRDTFuelComposition::PutTRUComposition(int nuc,string *id,real *inp,MATIDTranslator &midt)
{
  int *idi=new int[nuc];
  for(int i=0;i<nuc;i++){
    idi[i]=midt.ID(id[i]);
  };
  PutTRUComposition(nuc,idi,inp);
  delete [] idi;
};

void FRDTFuelComposition::PutTRUComposition(real rat_pu, real rat_ma, int nuc_pu, string *id_pu, real *inp_pu, int nuc_ma, string *id_ma, real *inp_ma, MATIDTranslator &midt)
{
  int *idi_pu=new int[nuc_pu];
  int *idi_ma=new int[nuc_ma];  
  for(int i=0;i<nuc_pu;i++){
    idi_pu[i]=midt.ID(id_pu[i]);
  };
  for(int i=0;i<nuc_ma;i++){
    idi_ma[i]=midt.ID(id_ma[i]);
  };
  PutTRUComposition(rat_pu,rat_ma,nuc_pu,idi_pu,inp_pu,nuc_ma,idi_ma,inp_ma);
  delete [] idi_pu;
  delete [] idi_ma;  
};
  
void FRDTFuelComposition::PutTRUComposition(real rat_pu, real rat_ma, int nuc_pu, int *id_pu, real *inp_pu, int nuc_ma, int *id_ma, real *inp_ma)
{
  real rat_sum=rat_pu+rat_ma;
  rat_pu/=rat_sum;
  rat_ma/=rat_sum;

  for(int i=0;i<nuc_pu;i++){
    int tmp=id_pu[i];
    for(int j=0;j<nuc_ma;j++){
      if(tmp==id_ma[j]){
	cout<<"# Error in FRDTFuelComposition::PutTRUComposition.\n";
	cout<<"# Nuclide "<<tmp<<" is included in both Pu and MA.\n";
	exit(0);
      };
    };
  };

  real inp_pu_sum=0.;
  for(int i=0;i<nuc_pu;i++){
    inp_pu_sum+=inp_pu[i];
  };
  for(int i=0;i<nuc_pu;i++){
    inp_pu[i]/=inp_pu_sum;
  };
  /*
  if(inp_pu_sum<99.999||inp_pu_sum>100.0001){
    cout<<"# Error in FRDTFuelComposition::PutTRUComposition.\n";
    cout<<"# The sum of weight ratios of Pu is NOT unity.\n";
    exit(0);
  };
  */
  
  real inp_ma_sum=0.;
  for(int i=0;i<nuc_ma;i++){
    inp_ma_sum+=inp_ma[i];
  };
  for(int i=0;i<nuc_ma;i++){
    inp_ma[i]/=inp_ma_sum;
  };
  /*
  if(inp_ma_sum<99.999||inp_ma_sum>100.0001){
    cout<<"# Error in FRDTFuelComposition::PutTRUComposition.\n";
    cout<<"# The sum of weight ratios of MA is NOT unity.\n";
    exit(0);
  };
  */
  
  int nuc=nuc_pu+nuc_ma;
  int *id=new int[nuc];
  real *inp=new real[nuc];

  for(int i=0;i<nuc_pu;i++){
    id[i]=id_pu[i];
    inp[i]=inp_pu[i]*rat_pu*100.;
  };
  for(int i=0;i<nuc_ma;i++){
    id[nuc_pu+i]=id_ma[i];
    inp[nuc_pu+i]=inp_ma[i]*rat_ma*100.;
  };

  PutTRUComposition(nuc,id,inp);
  
  delete [] id;
  delete [] inp;
};

void FRDTFuelComposition::PutTRUComposition(int nuc,int *id,real *inp)
{
  // (unit of input data is %)
  for(int i=0;i<nuc;i++){
    inp[i]*=0.01;
  };

  // (weight-ratio is normalized to 1)
  real sum=0.;
  for(int i=0;i<nuc;i++){
    sum+=inp[i];
  };
  for(int i=0;i<nuc;i++){
    inp[i]*=1./sum;
  };

  // 
  vector<real> wgt_rat(nucnum_hm,0.);

  real pu_per_tru=0.;
  real pu_fissile=0.;
  for(int i=0;i<nuc;i++){
    if(id[i]>=940000&&id[i]<950000)pu_per_tru+=inp[i];
    if(id[i]==942390||id[i]==942410)pu_fissile+=inp[i];
  };

  if(pu_enrichment==0.){
    pu_enrichment=pu_fissile_enrichment*pu_per_tru/pu_fissile;
    if(pu_enrichment>1.){
      cout<<"# Error in FRDTFuelComposition::PutTRUComposition.\n";
      cout<<"# Pu-fissile-enrichment is too high.\n";
      exit(0);
    };
  };
  
  // weight ratio
  real total_hm=1.0;
  real total_u=1.0-pu_enrichment/pu_per_tru;
  real total_tru=1.0-total_u;

  wgt_rat[0]=total_u*u5_enrichment; // (U-235)
  if(nucid[1]==922380){
    wgt_rat[1]=total_u-wgt_rat[0]; // (U-238)
  }else{
    wgt_rat[2]=total_u-wgt_rat[0]; // (U-238)
  };
  real tru_enrichment=(total_hm-total_u)/total_hm; // wt
  real u_enrichment=total_u/total_hm; // wt
  cout<<"# TRU (HM) enrichment : "<<tru_enrichment*100<<" wt%\n";

  for(int i=0;i<nuc;i++){
    bool flag=false;
    for(int j=0;j<nucnum_hm;j++){
      if(id[i]==nucid[j]){
	wgt_rat[j]=total_tru*inp[i];
	flag=true;        
      };
    };
    if(!flag){
      cout<<"# Error in FRDTFuelComposition::PutTRUComposition.\n";
      cout<<"# Nuclide "<<id[i]<<" cannot be treated.\n";
      exit(0);
    };
  };

  vector<real> den_rat(nucnum_hm,0.);
  MATIDTranslator midt;
  for(int i=0;i<nucnum_hm;i++){
    //den_rat[i]=wgt_rat[i]/atomic_weight_hm[i];
    //real aw=bu.GetAtomicWeight(nucid[i]);
    real aw=amd.GetData(nucid[i]);
      if(aw==0.){
      int iz,ia,il;
      midt.GetParameter(nucid[i],iz,ia,il);
      aw=real(ia);
    };
    den_rat[i]=wgt_rat[i]/aw;
    //cout<<i<<" "<<nucid[i]<<" "<<den_rat[i]<<"\n";
  };

  // (normalization to 1.0)
  real wgt_sum=0.;
  real den_sum=0.;
  for(int i=0;i<nucnum_hm;i++){
    wgt_sum+=wgt_rat[i];
    den_sum+=den_rat[i];
  };
  for(int i=0;i<nucnum_hm;i++){
    wgt_rat[i]/=wgt_sum;
    den_rat[i]/=den_sum;
  }; 

  /*
  for(int i=0;i<nucnum_hm;i++){
    cout<<nucid[i]<<" "<<wgt_rat[i]<<" "<<den_rat[i]<<" ";
    cout<<"\n";
  };
  */

  // (TRU averaged weight)
  real avg_wgt_tru=0.;
  if(total_tru>=0.){
    for(int i=nucnum_u;i<nucnum_hm;i++){
      //avg_wgt_tru+=wgt_rat[i]*atomic_weight_hm[i];
      //avg_wgt_tru+=wgt_rat[i]*bu.GetAtomicWeight(nucid[i]);
      avg_wgt_tru+=wgt_rat[i]*amd.GetData(nucid[i]);
    };
    avg_wgt_tru/=tru_enrichment;
  };

  // (U averaged weight)
  real avg_wgt_u=0.;
  for(int i=0;i<nucnum_u;i++){
    //avg_wgt_u+=wgt_rat[i]*atomic_weight_hm[i];
    //avg_wgt_u+=wgt_rat[i]*bu.GetAtomicWeight(nucid[i]);
    avg_wgt_u+=wgt_rat[i]*amd.GetData(nucid[i]);
  };
  if(avg_wgt_u!=0.)avg_wgt_u/=(1.-tru_enrichment);

  // (TRU-O2 enrichment)
  real tru_o2_er=0.;
  if(total_tru>0.){
    real tmp=0.;
    real tru_enrichment_at=0.;
    if(tru_enrichment>0.){
      tmp+=tru_enrichment/avg_wgt_tru;
      tru_enrichment_at=tru_enrichment/avg_wgt_tru;
    };
    if(u_enrichment>0.)tmp+=u_enrichment/avg_wgt_u;
    tru_enrichment_at/=tmp;
    //tru_o2_er=tru_enrichment_at*(avg_wgt_tru+o_m*atomic_weight_o);
    tru_o2_er=tru_enrichment_at*(avg_wgt_tru+o_m*amd.GetData(80160));
    //tru_o2_er/=tru_o2_er+(1.-tru_enrichment_at)*(avg_wgt_u+o_m*atomic_weight_o);
    tru_o2_er/=tru_o2_er+(1.-tru_enrichment_at)*(avg_wgt_u+o_m*amd.GetData(80160));
  };

  real tru_ox_density=0.48*tru_o2_er+10.96-2.5*(2.-o_m); 

  // [given in JNC TN8410 2000-011, p.123]

  real tru_density=0.;
  real u_density=0.;
  //if(tru_o2_er>0.)tru_density=tru_ox_density*tru_o2_er/(avg_wgt_tru+o_m*atomic_weight_o)*avo;
  if(tru_o2_er>0.)tru_density=tru_ox_density*tru_o2_er/(avg_wgt_tru+o_m*amd.GetData(80160))*avo;
  //if(tru_o2_er<1.)u_density=tru_ox_density*(1.-tru_o2_er)/(avg_wgt_u+o_m*atomic_weight_o)*avo;
  if(tru_o2_er<1.)u_density=tru_ox_density*(1.-tru_o2_er)/(avg_wgt_u+o_m*amd.GetData(80160))*avo;

  real den_rat_u=0.;
  for(int i=0;i<nucnum_u;i++){
    den_rat_u+=den_rat[i];
  };
  real den_rat_tru=1.-den_rat_u;

  for(int i=0;i<nucnum_u;i++){
    density[i]=0.;
    if(den_rat_u>0.)density[i]=u_density*den_rat[i]/den_rat_u*pellet_theoretical_density;
  };
  if(total_tru>0.){
    for(int i=nucnum_u;i<nucnum_hm;i++){
      density[i]=tru_density*den_rat[i]/den_rat_tru*pellet_theoretical_density;
    };
  };

  density[nucnum-1]=(tru_density+u_density)*o_m*pellet_theoretical_density; // O-16
};

void FRDTFuelComposition::ShowPuEnrichment()
{
  cout<<"\n\n";
  cout<<"  Pu   enrichment : "<<pu_enrichment<<"\n";
  cout<<"  Pu-f enrichment : "<<pu_fissile_enrichment<<"\n";
};

void FRDTFuelComposition::ShowTRUWeightInfo()
{
  real wgt_sum_hm=0.;
  real wgt_sum_tru=0.;
  real wgt_sum_ma=0.;    
  for(int i=0;i<nucnum;i++){
    int id=nucid[i];
    real tmp=density[i]*amd.GetData(id);    
    if(id>900000&&id<999990){
      wgt_sum_hm+=tmp;
      if(id>930000&&id<999990){
        wgt_sum_tru+=tmp;
        if(id<940000||id>950000){
          wgt_sum_ma+=tmp;
	};
      };
    };
  };

  cout<<"#---------------------------------------------------------------------\n";
  cout<<"# Actinoids weight ratio in fuel pellet ([in HM] / [in TRU])\n";
  cout<<"#\n";
  for(int i=0;i<nucnum;i++){
    int id=nucid[i];
    if(id>900000&&id<999990&&density[i]>0){
      cout<<"#   "<<id<<" : "<<density[i]*amd.GetData(id)/wgt_sum_hm;
      if(id>930000){
        cout<<" / "<<density[i]*amd.GetData(id)/wgt_sum_tru;
      };
      cout<<"\n";
    };
  };
  cout<<"#\n";
  cout<<"# Weight ratio of MA (non-Pu) ([in HM]/[in TRU]) : "<<wgt_sum_ma/wgt_sum_hm<<" / "<<wgt_sum_ma/wgt_sum_tru<<"\n";
  cout<<"#---------------------------------------------------------------------\n";  
};

void FRDTFuelComposition::PutUO2Composition(int nuc,int *id,real *inp)
{
  real sum=0.;
  for(int i=0;i<nuc;i++){
    sum+=inp[i];
  };

  uo2_comp_id.clear();
  uo2_comp_wgt.clear();
  for(int i=0;i<nuc;i++){
    if(id[i]<920000||id[i]>=930000){
      cout<<"# Error in FRDTFuelComposition::PutUO2Composition.\n";
      cout<<"# Nuclide "<<id[i]<<" is NOT accepted in the UO2 fuel.\n";
      exit(0);
    };
    uo2_comp_id.push_back(id[i]);
    uo2_comp_wgt.push_back(inp[i]/sum);
  };
};

void FRDTFuelComposition::PutPuO2Composition(int nuc,int *id,real *inp)
{
  real sum=0.;
  for(int i=0;i<nuc;i++){
    sum+=inp[i];
  };

  puo2_comp_id.clear();
  puo2_comp_wgt.clear();
  for(int i=0;i<nuc;i++){
    if(id[i]<940000||id[i]>=960000){
      cout<<"# Error in FRDTFuelComposition::PutPuO2Composition.\n";
      cout<<"# Nuclide "<<id[i]<<" is NOT accepted in the PuO2 fuel.\n";
      exit(0);
    };
    puo2_comp_id.push_back(id[i]);
    puo2_comp_wgt.push_back(inp[i]/sum);
  };
};

void FRDTFuelComposition::CalNumberDensity(real uo2_wgt,real puo2_wgt,real tho2_wgt)
{
  cout<<"# WARNING !!!\n";
  cout<<"# FRDTFuelComposition::CalNumberDensity.\n";
  cout<<"# Map variable <atomic_weight> should be replaced by instance of the burnup class.\n";
  exit(0);

  map<int,real> atomic_weight;
  // Atomic weight (atomic mass unit)
  real anmu=1.008665; // From SRAC burnup-chain data
  // (from neutron mass unit to atomic mass unit)
  // Th
  atomic_weight[902320]=230.045*anmu;
  // Pa
  atomic_weight[912310]=229.051*anmu;
  atomic_weight[912330]=229.051*anmu;
  // U
  atomic_weight[922320]=230.044*anmu;
  atomic_weight[922330]=231.038*anmu;
  atomic_weight[922340]=232.03*anmu;
  atomic_weight[922350]=233.025*anmu;
  atomic_weight[922360]=234.018*anmu;
  atomic_weight[922370]=235.012*anmu;
  atomic_weight[922380]=236.006*anmu;
  // Np
  atomic_weight[932360]=233.973*anmu;
  atomic_weight[932370]=235.012*anmu;
  atomic_weight[932390]=236.999*anmu;
  // Pu
  atomic_weight[942360]=234.018*anmu;
  atomic_weight[942380]=236.005*anmu;
  atomic_weight[942390]=236.999*anmu;
  atomic_weight[942400]=237.992*anmu;
  atomic_weight[942410]=238.986*anmu;
  atomic_weight[942420]=239.979*anmu;
  // Am
  atomic_weight[952410]=238.986*anmu;
  atomic_weight[952420]=239.98*anmu;
  atomic_weight[952421]=239.98*anmu;
  atomic_weight[952430]=240.973*anmu;
  // Cm
  atomic_weight[962420]=239.980*anmu;
  atomic_weight[962430]=240.972*anmu;
  atomic_weight[962440]=241.966*anmu;
  atomic_weight[962450]=242.961*anmu;
  atomic_weight[962460]=243.953*anmu;

  nucnum=0;
  nucid.clear();
  density.clear();

  // ++ Calculation of average weight for each atom
  real atomic_weight_u=0.;
  int uo2_nuc=uo2_comp_id.size();
  for(int i=0;i<uo2_nuc;i++){
    int nid=uo2_comp_id[i];
    atomic_weight_u+=atomic_weight[nid]*uo2_comp_wgt[i];
  };

  real atomic_weight_pu=0.;
  int puo2_nuc=puo2_comp_id.size();
  for(int i=0;i<puo2_nuc;i++){
    int nid=puo2_comp_id[i];
    atomic_weight_pu+=atomic_weight[nid]*puo2_comp_wgt[i];
  };

  //cout<<atomic_weight_u<<"\n";
  //cout<<atomic_weight_pu<<"\n";
  //exit(0);

  // ++ Weight ratio translation from HM to molecule
  //    (U:Pu:Th -> UO2:PuO2:ThO2)
  uo2_wgt*=(atomic_weight_u+amd.GetData(80160)*o_m)/atomic_weight_u;
  puo2_wgt*=(atomic_weight_pu+amd.GetData(80160)*o_m)/atomic_weight_pu;
  tho2_wgt*=(atomic_weight[902320]+amd.GetData(80160)*o_m)/atomic_weight[902320];

  real sum=uo2_wgt+puo2_wgt+tho2_wgt;
  uo2_wgt/=sum;
  puo2_wgt/=sum;
  tho2_wgt/=sum;

  // ++ Weight ratio translation from HM to molecule
  sum=0.;
  for(int i=0;i<uo2_nuc;i++){
    int nid=uo2_comp_id[i];
    uo2_comp_wgt[i]*=(atomic_weight[nid]+amd.GetData(80160)*o_m)/atomic_weight[nid];
    sum+=uo2_comp_wgt[i];
  };
  for(int i=0;i<uo2_nuc;i++){
    uo2_comp_wgt[i]/=sum;
  };

  sum=0.;
  for(int i=0;i<puo2_nuc;i++){
    int nid=puo2_comp_id[i];
    puo2_comp_wgt[i]*=(atomic_weight[nid]+amd.GetData(80160)*o_m)/atomic_weight[nid];
    sum+=puo2_comp_wgt[i];
  };
  for(int i=0;i<puo2_nuc;i++){
    puo2_comp_wgt[i]/=sum;
  };

  // Density [g/cm3]
  real den_uo2=10.96;
  real den_puo2=11.44;
  real den_tho2=9.87;

  real uo2_v=uo2_wgt/den_uo2;
  real puo2_v=puo2_wgt/den_puo2;
  real tho2_v=tho2_wgt/den_tho2;
  real tot_v=uo2_v+puo2_v+tho2_v;
  uo2_v/=tot_v;
  puo2_v/=tot_v;
  tho2_v/=tot_v;

  // Weight per unit volume [g/cm3]
  real uo2_w=den_uo2*uo2_v;
  real puo2_w=den_puo2*puo2_v;
  real tho2_w=den_tho2*tho2_v;
  
  real den_o=0.;

  // UO2
  for(int i=0;i<uo2_nuc;i++){
    int nid=uo2_comp_id[i];
    nucid.push_back(nid);
    nucnum++;
    real tmp=uo2_w*uo2_comp_wgt[i]/(atomic_weight[nid]+amd.GetData(80160)*o_m);
    // molecular of UO2 [mol/cm3] 
    tmp*=avo; // number density of molecular of UO2 [#/cm3]
    density.push_back(tmp);
    den_o+=tmp*o_m;
  };

  // PuO2
  for(int i=0;i<puo2_nuc;i++){
    int nid=puo2_comp_id[i];
    nucid.push_back(nid);
    nucnum++;
    real tmp=puo2_w*puo2_comp_wgt[i]/(atomic_weight[nid]+amd.GetData(80160)*o_m);
    // molecular of PuO2 [mol/cm3] 
    tmp*=avo; // number density of molecular of UO2 [#/cm3]
    density.push_back(tmp);
    den_o+=tmp*o_m;
  };

  // ThO2
  nucid.push_back(902320);
  nucnum++;
  real tmp=tho2_w/(atomic_weight[902320]+amd.GetData(80160)*o_m);
  // molecular of ThO2 [mol/cm3]
  tmp*=avo; // ND of ThO2
  density.push_back(tmp);
  den_o+=tmp*o_m;

  // O
  nucid.push_back(80160);
  nucnum++;
  density.push_back(den_o);

  for(int i=0;i<nucnum;i++){
    density[i]*=pellet_theoretical_density;
  };

};

// +++ FRDTMetalFuelComposition +++
FRDTMetalFuelComposition::FRDTMetalFuelComposition(real smear_density, real pu_enrichment, real zr_composition, int num_pu_iso, string *pu_name, real *pu_vector, AtomicMassData &amd)
{
  //real pellet_density=15.8; // [g/cm3] from CRIEPI review No.37, p.88
                            // U-15%Pu-10%Zr (text by Tohoku Univ. & PNC)

  real pellet_density=-1.; // [g/cm3] adjusted for IC
  if(zr_composition<10.1&&zr_composition>9.9)pellet_density=15.7; // [g/cm3] adjusted for IC
  if(zr_composition<6.1&&zr_composition>5.9)pellet_density=16.8; // [g/cm3] adjusted for OC  
  if(pellet_density<0.){
    cout<<"# Error in FRDTMetalFuelComposition::FRDTMetalFuelComposition.\n";
    cout<<"# Assined Zr composition rate is not yet coded.\n";
    exit(0);
  };

  real u5_enrichment=0.3;  // [wt%]

  MATIDTranslator midt;  

  real pu_ratio_in_tru=0.;
  for(int i=0;i<num_pu_iso;i++){
    if(midt.ID(pu_name[i])>940000&&midt.ID(pu_name[i])<950000)pu_ratio_in_tru+=pu_vector[i];
  };
  pu_enrichment/=(pu_ratio_in_tru*0.01);


  // ... Zr
  int zr_iso_num=5;
  int zr_iso[]={400900,400910,400920,400940,400960};
  real zr_ratio[]={0.5145, 0.1122, 0.1715, 0.1738, 0.028};
  real aw_zr_nat=0.;
  for(int i=0;i<zr_iso_num;i++){
    aw_zr_nat+=zr_ratio[i]*amd.GetData(zr_iso[i]);
  };
  //cout<<"# Atomic weight of Zr natural abundance : "<<aw_zr_nat<<"\n";

  real wgt_zr=pellet_density*(smear_density*0.01)*(zr_composition*0.01);
  real nd_zr=wgt_zr/aw_zr_nat*avo;
  AddData(400000,nd_zr);

  // ... U
  int u_iso_num=2;
  int u_iso[]={922350,922380};
  real u_ratio[]={u5_enrichment*0.01, 1.-u5_enrichment*0.01}; // weight ratio

  real tmpsum=0.;
  for(int i=0;i<u_iso_num;i++){
    u_ratio[i]/=amd.GetData(u_iso[i]); // ND ratio
    tmpsum+=u_ratio[i];
  };
  for(int i=0;i<u_iso_num;i++){
    u_ratio[i]/=tmpsum;
    //cout<<i<<" "<<u_ratio[i]<<"\n";
  };
  
  real aw_u_nat=0.;
  for(int i=0;i<u_iso_num;i++){
    aw_u_nat+=u_ratio[i]*amd.GetData(u_iso[i]);
  };

  real wgt_u=pellet_density*(smear_density*0.01)*(1.-zr_composition*0.01)*(1.0-pu_enrichment*0.01);
  real nd_u=wgt_u/aw_u_nat*avo;
  AddData(922350, nd_u*u_ratio[0]);
  AddData(922380, nd_u*u_ratio[1]);  

  // ... Pu
  vector<real> pu_vector_nd(num_pu_iso);
  {
  real tmpsum=0.;
  for(int i=0;i<num_pu_iso;i++){
    pu_vector_nd[i]=pu_vector[i]/amd.GetData(midt.ID(pu_name[i])); // ND ratio
    tmpsum+=pu_vector_nd[i];
  };
  for(int i=0;i<num_pu_iso;i++){
    pu_vector_nd[i]/=tmpsum;
    pu_vector_nd[i]*=100.; // -> [%]
    //cout<<i<<" "<<pu_vector[i]<<"\n";    
  };
  };
  
  real aw_pu_nat=0.;
  for(int i=0;i<num_pu_iso;i++){
    aw_pu_nat+=(pu_vector_nd[i]*0.01)*amd.GetData(midt.ID(pu_name[i]));
  };

  real wgt_pu=pellet_density*(smear_density*0.01)*(1.-zr_composition*0.01)*(pu_enrichment*0.01);  
  real nd_pu=wgt_pu/aw_pu_nat*avo;  
  for(int i=0;i<num_pu_iso;i++){
    AddData(midt.ID(pu_name[i]),nd_pu*(pu_vector_nd[i]*0.01));
  };

  /*  
  // ... Bond sodium
  real na_density=0.838008; // g/cm3
  real wgt_na=na_density*(1.-smear_density*0.01);
  real nd_na=wgt_na/(amd.GetData(110230))*avo;
  AddData(110230,nd_na);
  */  

  //ShowSelf();
  //exit(0);
};


// +++ FRDTSUSComposition +++

FRDTSUSComposition::FRDTSUSComposition()
{
  PutNucnum(7);
  nucid[0]=250550;
  nucid[1]=280000;
  nucid[2]=240000;
  nucid[3]=420000;
  nucid[4]=260000;
  nucid[5]=740000;
  nucid[6]=270590;
};

void FRDTSUSComposition::PutDensityAndRatio(real dinp,real *inp)
{
  real sum=0.;
  for(int i=0;i<nucnum;i++){
    sum+=inp[i];
  };
  if(sum<99.999||sum>100.001){
    cout<<"Error in FRDTSUSComposition::PutRatio.\n";
    cout<<"Sum of composition ratio is not 100.\n";
    exit(0);
  };

  for(int i=0;i<nucnum;i++){
    real tmp=amd.GetData(nucid[i]);
    if(tmp<=0.){
      cout<<"# Error in FRDTSUSComposition.\n";
      cout<<"# Atomic mass of nuclide "<<nucid[i]<<" is not defined.\n";
      exit(0);
    };
    density[i]=dinp*inp[i]*0.01/tmp*avo;
  };
};

// +++ FRDTB4CComposition +++

FRDTB4CComposition::FRDTB4CComposition(real b10c_wt,real theo_den)
{
  PutNucnum(3);
  nucid[0]=50100;
  nucid[1]=50110;
  nucid[2]=60000;

  real atw_b10=10.0129380;
  real atw_b11=11.0093050;
  real atw_c=12.;

  real sum=b10c_wt/atw_b10+(1.-b10c_wt)/atw_b11;
  real b10c=b10c_wt/atw_b10/sum; // wt->at

  real atw_b=b10c*atw_b10+(1.-b10c)*atw_b11;
  real b4c_den=-1.869300912e-3*b10c*100.+2.553556231;

  real b_den=b4c_den*4./(atw_b*4.+atw_c)*avo;
  b_den*=theo_den;
  density[0]=b_den*b10c;
  density[1]=b_den*(1.-b10c);
  density[2]=b_den/4.;
};

// +++ FRDTADSFuelComposition

FRDTADSFuelComposition::FRDTADSFuelComposition()
{
  zrn_wgt_density=6.38; // [g/cm3]
  puman_wgt_density=12.89; // [g/cm3]
  smear_density=0.85;

  num_non_fuel_data=0;
};

void FRDTADSFuelComposition::PutFuelWeightRatio(real zrn,real pun,real man,real pellet_fraction)
{
  real hmn_density=14.32; // [g/cm3] 
  real zrn_density=7.09; // [g/cm3] 
  real smear_density=0.9*0.85; // theoretical density * smear density

  if(num_non_fuel_data==0){
    cout<<"# Error in FRDTADSFuelComposition::PutFuelWeightRatio.\n";
    cout<<"# Please do [PutNonFuelData] in prior of this method.\n";
    exit(0);
  };

  if(nucnum!=num_non_fuel_data){
    cout<<"# Error in FRDTADSFuelComposition::PutFuelWeightRatio.\n";
    cout<<"# This instance of FRDTADSFuelComposition already has fuel data.\n";
    cout<<"# Please do [PutNonFuelData] again to initialize this instance.\n";
    exit(0);
  };

  /*
  // Fission product nuclides are registered in the main program
  int lump_fp_id[]={9950000, 9980000, 9990000, 9910000};
  for(int i=0;i<4;i++){
    AddData(lump_fp_id[i],0.);
  };
  */

  real sum=zrn+pun+man;
  if(sum<0.999||sum>1.001){
    cout<<"# FRDTADSFuelComposition::PutFuelWeightRatio.\n";
    cout<<"# Ratio is normalized to unity from "<<sum<<"\n";
  };

  //real n15_mass=bu.GetAtomicWeight(70150);
  real n15_mass=amd.GetData(70150);
  real zrn_avg=91.224+n15_mass; // atomic weight of ZrN
 
  real pun_avg=0.;
  for(int i=0;i<pu_nuc;i++){
    //pun_avg+=pu_isotopic_ratio[i]/bu.GetAtomicWeight(pu_composition[i]);
    pun_avg+=pu_isotopic_ratio[i]/amd.GetData(pu_composition[i]);
    //cout<<pu_isotopic_ratio[i]<<" "<<atomic_weight[pu_composition[i]]<<" "<<pun_avg<<"\n";
    //if(bu.GetAtomicWeight(pu_composition[i])<=0.){
    if(amd.GetData(pu_composition[i])<=0.){
      cout<<"# Zero atomic weight is detected.\n";
      cout<<"# Nuclide is "<<pu_composition[i]<<"\n";
      exit(0);
    };
  };
  pun_avg=1./pun_avg;
  pun_avg+=n15_mass;

  real man_avg=0.;
  for(int i=0;i<ma_nuc;i++){
    //man_avg+=ma_isotopic_ratio[i]/bu.GetAtomicWeight(ma_composition[i]);
    man_avg+=ma_isotopic_ratio[i]/amd.GetData(ma_composition[i]);
    //cout<<ma_isotopic_ratio[i]<<" "<<atomic_weight[ma_composition[i]]<<" "<<man_avg<<"\n";
    //if(bu.GetAtomicWeight(ma_composition[i])<=0.){
    if(amd.GetData(ma_composition[i])<=0.){
      cout<<"# Zero atomic weight is detected.\n";
      cout<<"# Nuclide is "<<ma_composition[i]<<"\n";
      exit(0);
    };
  };
  man_avg=1./man_avg;
  man_avg+=n15_mass;

  //cout<<pun_avg<<" "<<man_avg<<"\n"; exit(0);

  // (Volume ratio)
  real zrn_volr=zrn/zrn_density;
  real pun_volr=pun/hmn_density;
  real man_volr=man/hmn_density;
  real tot=zrn_volr+pun_volr+man_volr;
  zrn_volr/=tot;
  pun_volr/=tot;
  man_volr/=tot;
  
  real pellet_density=hmn_density*(pun_volr+man_volr)+zrn_density*zrn_volr;
  pellet_density*=smear_density;
  pellet_density*=pellet_fraction;
  
  // (Mol per 1cm3);
  real zrn_mol=pellet_density*zrn/zrn_avg;
  real pun_mol=pellet_density*pun/pun_avg;
  real man_mol=pellet_density*man/man_avg;

  //cout<<zrn_molr<<" "<<pun_molr<<" "<<man_molr<<"\n";
  //cout<<(zrn_molr+pun_molr+man_molr)*avo*smear_density*pellet_fraction<<"\n"; exit(0);

  real n15_den=(zrn_mol+pun_mol+man_mol)*avo;  

  AddData(70150,n15_den);
  AddData(400000,zrn_mol*avo);
  for(int i=0;i<pu_nuc;i++){
    //AddData(pu_composition[i],pun_mol*avo*(pun_avg-n15_mass)*pu_isotopic_ratio[i]/bu.GetAtomicWeight(pu_composition[i]));
    AddData(pu_composition[i],pun_mol*avo*(pun_avg-n15_mass)*pu_isotopic_ratio[i]/amd.GetData(pu_composition[i]));
  };
  for(int i=0;i<ma_nuc;i++){
    AddData(ma_composition[i],man_mol*avo*(man_avg-n15_mass)*ma_isotopic_ratio[i]/amd.GetData(ma_composition[i]));
  };

};

void FRDTADSFuelComposition::PutPuData(int n, int *id, real *rat)
{
  pu_nuc=n;

  pu_composition.resize(pu_nuc);
  pu_isotopic_ratio.resize(pu_nuc);

  real sum=0.;
  for(int i=0;i<pu_nuc;i++){
    sum+=rat[i];
  };
  if(sum!=1.){
    cout<<"# FRDTADSFuelComposition::PutPuData.\n";
    cout<<"# Isotopic ratio is normalized to unity from "<<sum<<"\n";
  };
  for(int i=0;i<pu_nuc;i++){
    pu_composition[i]=id[i];
    pu_isotopic_ratio[i]=rat[i]/sum;
  };
};

void FRDTADSFuelComposition::PutMAData(int n, int *id, real *rat)
{
  ma_nuc=n;

  ma_composition.resize(ma_nuc);
  ma_isotopic_ratio.resize(ma_nuc);

  real sum=0.;
  for(int i=0;i<ma_nuc;i++){
    sum+=rat[i];
  };
  if(sum!=1.){
    cout<<"# FRDTADSFuelComposition::PutMAData.\n";
    cout<<"# Isotopic ratio is normalized to unity from "<<sum<<"\n";
  };
  for(int i=0;i<ma_nuc;i++){
    ma_composition[i]=id[i];
    ma_isotopic_ratio[i]=rat[i]/sum;
  };
};

void FRDTADSFuelComposition::PutNonfuelData(int n,int *id,real *den)
{
  num_non_fuel_data=n;

  PutNucnum(n);
  for(int i=0;i<nucnum;i++){
    nucid[i]=id[i];
    density[i]=den[i];
  };
};

// +++ FRDTBaseGeometry +++

FRDTBaseGeometry::FRDTBaseGeometry()
{
  pin_pitch=0.;
  pin_layer=0;

  // Default setting
  spacerwire_diameter=0.;
  pellet_invoid_diameter=0.;

  /*
  assembly_pitch=115.6;
  pin_outer_diameter=6.5;
  pin_thickness=0.47;
  number_of_pin=169;
  spacerwire_diameter=1.32;
  spacerwire_pitch=307;
  pellet_outer_diameter=5.4;

  */
};

real FRDTBaseGeometry::GetAssemblyVolume()
{
  return assembly_pitch*assembly_pitch*sqrt(3.)*0.5;
};

real FRDTBaseGeometry::GetCladVolume()
{
  real tmp=pin_outer_diameter*0.5;
  return (pow(tmp,2)-pow(tmp-pin_thickness,2))*PI*number_of_pin;
};

real FRDTBaseGeometry::GetSpacerwireVolume()
{
  if(number_of_pin==0)return 0.;

  if(spacerwire_diameter<1e-20)return 0.;

  real tmp1=pow(spacerwire_diameter*0.5,2);
  real tmp2=sqrt(pow(spacerwire_pitch,2)+pow(pin_outer_diameter*PI,2))/spacerwire_pitch*PI*number_of_pin;
  return tmp1*tmp2;
};

real FRDTBaseGeometry::GetPelletVolume()
{
  real tmp1=pow(pellet_outer_diameter*0.5,2)-pow(pellet_invoid_diameter*0.5,2);
  tmp1*=(PI*number_of_pin);
  return tmp1;
};

real FRDTBaseGeometry::GetPelletVoidVolume()
{
  real tmp1=pow(pellet_invoid_diameter*0.5,2);
  tmp1*=(PI*number_of_pin);
  return tmp1;
};

real FRDTBaseGeometry::GetPelletGapVolume()
{
  real tmp1=pow(pin_outer_diameter*0.5-pin_thickness,2);
  real tmp2=pow(pellet_outer_diameter*0.5,2);
  return (tmp1-tmp2)*PI*number_of_pin;
};

real FRDTBaseGeometry::GetWholePinVolume()
{
  return pow(pin_outer_diameter*0.5,2)*PI*number_of_pin;
};

// +++ FRDTSubAssemblyGeometry +++

FRDTSubAssemblyGeometry::FRDTSubAssemblyGeometry()
{
  duct_outersize=110.6;
  duct_thickness=3.;
};

real FRDTSubAssemblyGeometry::GetDuctVolume()
{
  return (pow(duct_outersize,2)-pow(duct_outersize-2.*duct_thickness,2))*sqrt(3.)*0.5;
};

real FRDTSubAssemblyGeometry::GetOutDuctVolume()
{
  return GetAssemblyVolume()-(pow(duct_outersize,2)*sqrt(3.)*0.5);
};

real FRDTSubAssemblyGeometry::GetCoolantVolume()
{
  return GetAssemblyVolume()-GetDuctVolume()-GetWholePinVolume()-GetSpacerwireVolume();
};

void FRDTSubAssemblyGeometry::CheckValidity()
{
  if(pin_pitch==0.){
    cout<<"# Error in FRDTSubAssembly::CheckValidity.\n";
    cout<<"# Pin pitch is NOT yet determined.\n";
    exit(0);
  };


  int layer=0;
  if(pin_layer==0){
    for(int i=1;i<20;i++){
      int num=3*i*(i-1)+1;
      if(num==number_of_pin)layer=i;
    };
    //pin_layer=layer;
    if(layer==0){
      cout<<"# Error in FRDTSubAssemblyGeometry::CheckValidity.\n";
      cout<<"# The number of fuel layers cannot be determined.\n";
      cout<<"# Please check the number of pins.\n";
      cout<<"# Current setting : "<<number_of_pin<<"\n";
      exit(0);
    };
  }else{
    layer=pin_layer;
    int num=3*pin_layer*(pin_layer-1)+1;
    if(num<number_of_pin){
      cout<<"# Error in FRDTSubAssemblyGeometry::CheckValidity.\n";
      cout<<"# The number of fuel layers is NOT appropriate\n";
      cout<<"# because the number of pins is larger than that calculated from layer.\n";
      exit(0);
    }else if(num>number_of_pin){
      cout<<"# "<<num-number_of_pin<<" fuel pins are removed from regular array.\n";
    };
  };

  real d=duct_outersize-(duct_thickness*2);
  real tmp=(d-pin_outer_diameter)/sqrt(3.)/(layer-1);
  if(tmp<pin_pitch){
    cout<<"# Error in FRDTSubAssemblyGeometry::CheckValidity.\n";
    cout<<"# Assembly geometory is invalid.\n";
    cout<<"# All the fuel pins cannot be located in the duct.\n";
    exit(0);
  };

  real tmp2=(d/sqrt(3.)-pin_pitch*(layer-1))*sqrt(3.)*0.5-pin_outer_diameter*0.5;
  cout<<"# Minimum distance from pin outer surface to duct wall : "<<tmp2<<"[mm]\n";

  if(tmp2>5.){
    cout<<"# Warning in FRDTSubAssemblyGeometry::CheckValidity.\n";
    cout<<"# Minimum distance from pin outer surface to duct wall is too large.\n";
    cout<<"# Minimum distance from pin outer surface to duct wall : "<<tmp2<<"[mm]\n";
  };

  string vol_err="";
  if(GetAssemblyVolume()<0.)vol_err="assembly";
  if(GetPelletVolume()<0.)vol_err="pellet";
  if(GetPelletVoidVolume()<0.)vol_err="void in pellet";
  if(GetPelletGapVolume()<0.)vol_err="void at pellet-clad gap";
  if(GetCladVolume()<0.)vol_err="cladding";
  if(GetSpacerwireVolume()<0.)vol_err="spacer wise";
  if(GetDuctVolume()<0.)vol_err="wrapper tube";
  if(GetCoolantVolume()<0.)vol_err="coolant";
  if(GetCoolantVolume()-GetOutDuctVolume()<0.)vol_err="coolant inside of wrapper tube";
  if(vol_err!=""){
    cout<<"#\n";
    cout<<"# Error in FRDTSubAssemblyGeometry::CheckValidity.\n";
    cout<<"# Negative volume is detected : "<<vol_err<<"\n";
    cout<<"#\n";
    exit(0);
  };
  /*
  real tmp1=pin_outer_diameter*(layer-1);
  real tmp2=(duct_outersize-(duct_thickness*2))/sqrt(3.);
  real tmp3=(tmp2-tmp1)*0.5*sqrt(3.);
  //cout<<tmp1<<" "<<tmp2<<" "<<tmp3<<"\n";
  //cout<<layer<<"\n";
  if(tmp3<pin_outer_diameter*0.5){
    cout<<"# Error in FRDTSubAssemblyGeometry::CheckValidity.\n";
    cout<<"# Assembly geometory is invalid.\n";
    cout<<"# All the fuel pins cannot be located in the duct.\n";
    exit(0);
  };
  */
};

void FRDTSubAssemblyGeometry::ShowVolumeInformation()
{
  real tot=GetAssemblyVolume();
  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"################################################################################\n";
  cout<<"# Region-wise volume of subassembly [mm^2] and fraction to total.\n";
  cout<<"#\n";
  cout<<"#   Total                    : "<<GetAssemblyVolume()<<"\n";
  cout<<"#   Fuel pellet              : "<<GetPelletVolume()<<" "<<GetPelletVolume()/tot<<"\n";
  cout<<"#   Void in pellet           : "<<GetPelletVoidVolume()<<" "<<GetPelletVoidVolume()/tot<<"\n"; 
  cout<<"#   Void at pellet-clad gap  : "<<GetPelletGapVolume()<<" "<<GetPelletGapVolume()/tot<<"\n";
  cout<<"#   Cladding                 : "<<GetCladVolume()<<" "<<GetCladVolume()/tot<<"\n";
  cout<<"#   Spacer wire              : "<<GetSpacerwireVolume()<<" "<<GetSpacerwireVolume()/tot<<"\n";
  cout<<"#   Wrapper tube             : "<<GetDuctVolume()<<" "<<GetDuctVolume()/tot<<"\n";
  cout<<"#   Coolant                  : "<<GetCoolantVolume()<<" "<<GetCoolantVolume()/tot<<"\n";
  cout<<"#   (inside of wrapper tube) : "<<GetCoolantVolume()-GetOutDuctVolume()<<"\n";
  cout<<"###############################################################################\n";
};

// +++ FRDTControlRodGeometry +++

real FRDTControlRodGeometry::GetGuidetubeVolume()
{
  real tmp1=pow(guidetube_outer_diameter*0.5,2);
  real tmp2=pow(guidetube_outer_diameter*0.5-guidetube_thickness,2);
  return (tmp1-tmp2)*PI;
};

real FRDTControlRodGeometry::GetShieldtubeVolume()
{
  real tmp1=pow(shieldtube_outer_diameter*0.5,2);
  real tmp2=pow(shieldtube_outer_diameter*0.5-shieldtube_thickness,2);
  return (tmp1-tmp2)*PI;
};

real FRDTControlRodGeometry::GetCoolantVolume()
{
  real ret=GetAssemblyVolume();
  ret-=GetGuidetubeVolume();
  ret-=GetShieldtubeVolume();
  ret-=GetWholePinVolume();
  ret-=GetSpacerwireVolume();
  return ret;
};

real FRDTControlRodGeometry::GetShieldtubeInsideVolume()
{
  return PI*pow(shieldtube_outer_diameter*0.5-shieldtube_thickness,2);
};

// +++ FRDTSubAssembly

FRDTSubAssembly::FRDTSubAssembly()
{
  sodium_number_density=0.;
  exist_fuel=false;
  exist_sus=false;
  exist_duct=false;
  exist_clad=false;
  exist_geom=false;
};

void FRDTSubAssembly::PutSodiumDensity(real inp)
{
  sodium_number_density=inp/amd.GetData(110230)*avo;
  //cout<<"# Sodium ND : "<<sodium_number_density<<"\n";
};

void FRDTSubAssembly::PutTemperature(real tmp1,real tmp2)
{
  if(tmp2==0.)tmp2=tmp1;
  temperature=tmp1;
  sodium_temperature=tmp2;
};

void FRDTSubAssembly::CheckExistAllData()
{
  if(!exist_fuel){
    cout<<"Error in FRDTSubAssembly::CalHomogenizedNumberDensity.\n";
    cout<<"Fuel data is not yet put to this instance.\n";
    exit(0);
  };
  if(!exist_duct){
    cout<<"Error in FRDTSubAssembly::CalHomogenizedNumberDensity.\n";
    cout<<"Duct data is not yet put to this instance.\n";
    exit(0);
  };
  if(!exist_clad){
    cout<<"Error in FRDTSubAssembly::CalHomogenizedNumberDensity.\n";
    cout<<"Cladding data is not yet put to this instance.\n";
    exit(0);
  };
  if(!exist_geom){
    cout<<"Error in FRDTSubAssembly::CalHomogenizedNumberDensity.\n";
    cout<<"Geometry data is not yet put to this instance.\n";
    exit(0);
  };
};

void FRDTSubAssembly::ShowHomogenizedNumberDensity()
{
  CheckExistAllData();

  cout<<"#\n# Fuel pellet\n#\n";
  real factor=fuelsa_geometry->GetPelletVolume()/fuelsa_geometry->GetAssemblyVolume();
  fuel_composition->ShowSelf(factor);  

  cout<<"#\n# Coolant\n#\n";  
  cout<<"110230  "<<sodium_number_density*fuelsa_geometry->GetCoolantVolume()/fuelsa_geometry->GetAssemblyVolume()<<"\n";

  cout<<"#\n# Cladding & Spacer-wire\n#\n";    
  real factor2=(fuelsa_geometry->GetPinVolume()+fuelsa_geometry->GetSpacerwireVolume())/fuelsa_geometry->GetAssemblyVolume();
  clad_composition->ShowSelf(factor2);

  cout<<"#\n# Duct\n#\n";    
  real factor3=(fuelsa_geometry->GetDuctVolume())/fuelsa_geometry->GetAssemblyVolume();
  duct_composition->ShowSelf(factor3);
};

void FRDTSubAssembly::ShowHeterogeneousNumberDensity()
{
  cout<<"+++ Fuel +++\n";
  fuel_composition->ShowSelf();
  cout<<"+++ SUS +++\n";
  sus_composition->ShowSelf();
  cout<<"+++ Coolant and spacer wire +++\n";
  real cvol=fuelsa_geometry->GetCoolantVolume();
  real wvol=fuelsa_geometry->GetSpacerwireVolume();
  cout<<"1125  "<<sodium_number_density*cvol/(cvol+wvol)<<"\n";
  sus_composition->ShowSelf(wvol/(cvol+wvol));
};

void FRDTSubAssembly::PutHomogenizedNumberDensity(Medium &med)
{
  CheckExistAllData();

  real factor=fuelsa_geometry->GetPelletVolume()/fuelsa_geometry->GetAssemblyVolume();
  // Volume fraction of fuel pellet to whole assembly

  // Homogenized number density of fuel pellet is calculated.
  for(int i=0;i<fuel_composition->GetNucnum();i++){
    //med.GetNuclide(fuel_composition->GetNucid(i)).PutDensity(fuel_composition->GetDensity(i)*factor);
    if(med.ExistNuclide(fuel_composition->GetNucid(i))){
      med.GetNuclide(fuel_composition->GetNucid(i)).PutDensity(fuel_composition->GetDensity(i)*factor);
    };
  };

  real factor3=fuelsa_geometry->GetCoolantVolume()/fuelsa_geometry->GetAssemblyVolume();
  // Volume fraction of sodium region to whole assembly

  // Homogenized number density of sodium is calculated.
  med.GetNuclide(110230).PutDensity(sodium_number_density*factor3);
  //med.GetNuclide(110230).AddDensity(sodium_number_density*factor3);

  /*
  real factor2=(fuelsa_geometry->GetDuctVolume()+fuelsa_geometry->GetPinVolume()+fuelsa_geometry->GetSpacerwireVolume())/fuelsa_geometry->GetAssemblyVolume();
  for(int i=0;i<sus_composition->GetNucnum();i++){
    //med.GetNuclide(sus_composition->GetNucid(i)).PutDensity(sus_composition->GetDensity(i)*factor2);
    med.GetNuclide(sus_composition->GetNucid(i)).AddDensity(sus_composition->GetDensity(i)*factor2);
  };
  */

  // Homogenized number density of cladding & spacer-wire is calculated
  real factor_clad=(fuelsa_geometry->GetPinVolume()+fuelsa_geometry->GetSpacerwireVolume())/fuelsa_geometry->GetAssemblyVolume();
  for(int i=0;i<clad_composition->GetNucnum();i++){
    med.GetNuclide(clad_composition->GetNucid(i)).AddDensity(clad_composition->GetDensity(i)*factor_clad);
  };

  // Homogenized number density of duct is calculated
  real factor_duct=fuelsa_geometry->GetDuctVolume()/fuelsa_geometry->GetAssemblyVolume();
  for(int i=0;i<duct_composition->GetNucnum();i++){
    med.GetNuclide(duct_composition->GetNucid(i)).AddDensity(duct_composition->GetDensity(i)*factor_duct);
  };

};

void FRDTSubAssemblyGeometry::MonjuFuelSubAssembly()
{
  PutAssemblyPitch(115.6);
  PutDuctOutersize(110.6);
  PutDuctThickness(3.);
  PutPinOuterDiameter(6.5);
  PutPinThickness(0.47);
  PutNumberOfPin(169);
  //PutSpacerwireDiameter(0.);
  PutSpacerwireDiameter(1.32);
  PutSpacerwirePitch(307.);
  PutPelletOuterDiameter(5.4);
  PutPelletInvoidDiameter(0.);
  PutPinPitch(7.9);
};

void FRDTSubAssemblyGeometry::MonjuRadialBlanketSubAssembly()
{
  PutAssemblyPitch(115.6);
  PutDuctOutersize(110.6);
  PutDuctThickness(3.);
  PutPinOuterDiameter(11.6);
  PutPinThickness(0.5);
  PutNumberOfPin(61);
  real wrad1=0.9;
  real wrad2=1.5;
  real wrad=(wrad1*24+wrad2*37)/61;
  //PutSpacerwireDiameter(0.);
  PutSpacerwireDiameter(wrad);
  PutSpacerwirePitch(251.);
  PutPelletOuterDiameter(10.4);
  PutPelletInvoidDiameter(0.);
  PutPinPitch(13.0);
};

void FRDTSubAssemblyGeometry::JSFR1500FuelSubAssembly()
{
  // JSFR-1500 assembly model
  //
  // Reference
  // 
  // [1] Ohki, Physor2008
  // [2] Ogawa, JAEA-Research 2007-084
  // [3] FS Phase-II, JAEA-Research 2006-042

  PutAssemblyPitch(206.0); // Table 2-1 in [2]
  PutDuctOutersize(201.6); // Table 2-4 in [2]
  PutDuctThickness(5.5); // Table 2-4 in [2]
  // Original is 5.0, but inner duct component is artificially added.
  PutDuctThickness(5.5); // Table 2-4 in [2]
  PutPinOuterDiameter(10.4); // Table 2-4 in [2]
  PutPinThickness(0.71); // Table 2-4 in [2]
  //PutNumberOfPin(271);// Table 3 in [1]
  PutNumberOfPin(255);// Table 3 in [1]
  // Actually 255 because of duct region
  PutSpacerwireDiameter(1.03); // Table 2-4 in [2]
  PutSpacerwirePitch(307.);
  PutPelletOuterDiameter(10.4-0.71*2);  // No gap
  PutPelletInvoidDiameter(0.);
  PutPinPitch(11.48); // Table 2-4 in [2]
  PutPinLayer(10);
};

void FRDTSubAssemblyGeometry::JSFR1500RadialBlanketSubAssembly()
{
  // JSFR-1500 assembly model
  //
  // Reference
  // 
  // [1] Ohki, Physor2008
  // [2] Ogawa, JAEA-Research 2007-084
  // [3] FS Phase-II, JAEA-Research 2006-042

  PutAssemblyPitch(206.0); // Table 2-1 in [2]
  PutDuctOutersize(201.6); // Table 2-4 in [2]
  PutDuctThickness(5.0); // Table 2-4 in [2]
  PutPinOuterDiameter(11.7); // Table 2-4 in [2]
  PutPinThickness(0.42); // Table 2-4 in [2]
  PutNumberOfPin(217);// Table 2-4 in [2]
  PutSpacerwireDiameter(1.07); // Table 2-4 in [2]
  PutSpacerwirePitch(307.);
  PutPelletOuterDiameter(11.7-0.42*2);  // No gap
  PutPelletInvoidDiameter(0.);
  PutPinPitch(12.82); // Table 2-4 in [2]
  PutPinLayer(9);
};

// +++ FRDTControlRod +++

void FRDTControlRod::PutSodiumDensity(real inp)
{
  sodium_number_density=inp/amd.GetData(110230)*avo;
};

void FRDTControlRod::ShowHomogenizedNumberDensity()
{
  real tv=cr_geometry->GetAssemblyVolume();
  real factor=cr_geometry->GetPelletVolume()/tv;
  b4c_composition->ShowSelf(factor);  

  cout<<"110230  "<<sodium_number_density*cr_geometry->GetCoolantVolume()/tv<<"\n";

  real ssv=cr_geometry->GetGuidetubeVolume();
  ssv+=cr_geometry->GetShieldtubeVolume();
  ssv+=cr_geometry->GetCladVolume();
  ssv+=cr_geometry->GetSpacerwireVolume();
  real factor2=ssv/tv;
  sus_composition->ShowSelf(factor2);
};

void FRDTControlRod::PutHomogenizedNumberDensity(Medium &med)
{
  vector<int> nucid;
  vector<real> den;

  real tv=cr_geometry->GetAssemblyVolume();
  real factor=cr_geometry->GetPelletVolume()/tv;

  int tmp=b4c_composition->GetNucnum();
  for(int i=0;i<tmp;i++){
    nucid.push_back(b4c_composition->GetNucid(i));
    den.push_back(b4c_composition->GetDensity(i)*factor);
  };


  nucid.push_back(110230);
  den.push_back(sodium_number_density*cr_geometry->GetCoolantVolume()/tv);

  real ssv=cr_geometry->GetGuidetubeVolume();
  ssv+=cr_geometry->GetShieldtubeVolume();
  ssv+=cr_geometry->GetCladVolume();
  ssv+=cr_geometry->GetSpacerwireVolume();
  real factor2=ssv/tv;

  int tmp2=sus_composition->GetNucnum();
  for(int i=0;i<tmp2;i++){
    nucid.push_back(sus_composition->GetNucid(i));
    den.push_back(sus_composition->GetDensity(i)*factor2);
  };

  int nucnum=nucid.size();
  med.PutNuclide(nucnum,nucid,den);
  //sus_composition->ShowSelf(factor2);
};

void FRDTControlRod::ShowShieldtubeInsideHomogenizedNumberDensity()
{
  real tv=cr_geometry->GetShieldtubeInsideVolume();
  real factor=cr_geometry->GetPelletVolume()/tv;
  b4c_composition->ShowSelf(factor);  

  real cv=tv-cr_geometry->GetWholePinVolume();
  cv-=cr_geometry->GetSpacerwireVolume();
  cout<<"110230  "<<sodium_number_density*cv/tv<<"\n";

  real ssv=cr_geometry->GetCladVolume();
  ssv+=cr_geometry->GetSpacerwireVolume();
  real factor2=ssv/tv;
  sus_composition->ShowSelf(factor2);
};

void FRDTControlRod::ShowGapSmearedCladNumberDensity()
{
  real v1=cr_geometry->GetWholePinVolume();
  real v2=cr_geometry->GetCladVolume();
  real v3=cr_geometry->GetPelletVolume();
  real factor=v2/(v1-v3);
  sus_composition->ShowSelf(factor);
};

void FRDTControlRod::ShowSpacerwireSmearedCoolantNumberDensity()
{
  real total_v=cr_geometry->GetShieldtubeInsideVolume();
  real wpin_v=cr_geometry->GetWholePinVolume();
  real swire_v=cr_geometry->GetSpacerwireVolume();
  real na_v=total_v-wpin_v-swire_v;
  cout<<"110230  "<<sodium_number_density*na_v/(na_v+swire_v)<<"\n";
  real factor=swire_v/(swire_v+na_v);
  sus_composition->ShowSelf(factor);
};

void FRDTControlRod::ShowSUSSmearedCoolantNumberDensity()
{
  real total_v=cr_geometry->GetShieldtubeInsideVolume();
  real wpin_v=cr_geometry->GetWholePinVolume();
  real swire_v=cr_geometry->GetSpacerwireVolume();
  real pellet_v=cr_geometry->GetPelletVolume()+cr_geometry->GetPelletVoidVolume();
  real sus_v=swire_v+cr_geometry->GetCladVolume();
  real na_v=total_v-wpin_v-swire_v;
  real hom_v=total_v-pellet_v;
  cout<<"1125  "<<sodium_number_density*na_v/hom_v<<"\n";
  real factor=sus_v/hom_v;
  sus_composition->ShowSelf(factor);
};

void FRDTControlRod::ShowSelf()
{
  cout<<"\n# B4C pellet number density\n\n";
  b4c_composition->ShowSelf();

  cout<<"\n# Clad number density \n\n";
  sus_composition->ShowSelf();

  cout<<"\n#   (Gap-smeared density)\n\n";
  ShowGapSmearedCladNumberDensity();

  cout<<"\n# Sodium\n\n";
  cout<<"110230  "<<sodium_number_density<<"\n";

  cout<<"\n#   (Spacer wire-smeared coolant density)\n";
  cout<<"#    (inside shielding tube)\n\n";
  ShowSpacerwireSmearedCoolantNumberDensity();

  cout<<"\n# +++++++++++++++++++++++++++++++++++\n\n";
  cout<<"\n# Homogenized number density\n\n";
  ShowHomogenizedNumberDensity();

  cout<<"\n# Shieldtube inside homogenized number density\n\n";
  ShowShieldtubeInsideHomogenizedNumberDensity();

  cout<<"\n# Coolant-Clad-wire homogenized number density\n";
  cout<<"#  (inside shielding tube)(for ring model)\n\n";
  ShowSUSSmearedCoolantNumberDensity();

  /*
  real tv=cr_geometry->GetShieldtubeInsideVolume();
  tv/=cr_geometry->GetNumberOfPin();
  cout<<sqrt(tv/PI)<<"\n";
  */
};

// +++ FRDTCoreDataManager +++

void FRDTCoreDataManager::CylinderDataGeneration(int siz, int *map, real pitch)
{
  real area=pitch*pitch*sqrt(3.)*0.5;

  int center=(siz-1)/2;

  int maxid=0;
  for(int i=0;i<siz*siz;i++){
    if(map[i]>maxid)maxid=map[i];
  };

  vector<bool> mape(siz*siz,false);
  mape[center*siz+center]=true;

  vector<bool> mape_org(siz*siz,false);

  real area_cum=area;
  real radius=sqrt(area_cum/PI);
  cout<<"Ring 0   ("<<radius<<")\n";
  real radius_old=radius;
  for(int i=0;i<center-1;i++){
    for(int j=0;j<siz*siz;j++){
      mape_org[j]=mape[j];
    };
    vector<int> idnum(maxid+1,0);
    int num=0;
    for(int y=1;y<siz-1;y++){
      for(int x=1;x<siz-1;x++){
	if(!mape[y*siz+x]){
	  int pos=y*siz+x;
	  bool flag=false;
          if(mape_org[pos+1])flag=true;
          if(mape_org[pos-1])flag=true;
          int tmp=0;
	  if(y%2==1)tmp=1;
	  if(mape_org[pos-siz-tmp])flag=true;
	  if(mape_org[pos-siz-tmp+1])flag=true;
	  if(mape_org[pos+siz-tmp])flag=true;
	  if(mape_org[pos+siz-tmp+1])flag=true;
          if(flag){
	    num++;
            idnum[map[pos]]++;
	    mape[pos]=true;
	  };
	};
      };
    };
    area_cum+=num*area;
    radius=sqrt(area_cum/PI);
    cout<<"Ring "<<i+1<<" "<<num<<"  ("<<radius<<" "<<radius-radius_old<<")\n";
    int mat=0;
    for(int j=0;j<maxid;j++){
      if(idnum[j]!=0){
        cout<<"   "<<j<<" "<<idnum[j]<<"\n";
	mat++;
      };
    };
    if(mat>1){
      int ist=-1;
      for(int j=0;j<maxid;j++){
	if(idnum[j]!=0&&ist==-1)ist=j;
      };
      real area2=area_cum-(num-idnum[ist])*area;
      real radius2=sqrt(area2/PI);
      cout<<"\n";
      cout<<"  case A : dominant region locates inside ... radius : "<<radius2<<" "<<radius2-radius_old<<"\n";
      cout<<"\n";

      real area3=area_cum-num*area;
      area3+=idnum[ist]/2*area;
      real radius3=sqrt(area3/PI);
      cout<<"  case B : \n";
      cout<<"   ... half of dominant region ... radius : "<<radius3<<" "<<radius3-radius_old<<"\n";
      area3+=(num-idnum[ist])*area;
      real radius4=sqrt(area3/PI);
      cout<<"   ... non-dominant region     ... radius : "<<radius4<<" "<<radius4-radius3<<"\n";
      cout<<"   ... half of dominant region ... radius : "<<radius<<" "<<radius-radius4<<"\n";
      cout<<"\n";

      area2=area_cum-idnum[ist]*area;
      radius2=sqrt(area2/PI);
      cout<<"  case C : non-dominant region locates inside ... radius : "<<radius2<<" "<<radius2-radius_old<<"\n";
      cout<<"\n";
    };
    radius_old=radius;
  };

};

// ++ FRDTAssypow

FRDTAssypow::FRDTAssypow(real ipitch)
{
  pitch=ipitch;

  real r=pitch*0.3333333333333333;
  real theta=-PI*0.33333333333;
  for(int i=0;i<6;i++){
    xpos[i]=r*cos(theta);
    ypos[i]=r*sin(theta);
    theta+=PI*0.333333333333333;
  };

  int index=6;
  real r2=pitch;
  real theta2=-PI*0.333333333333;
  for(int i=0;i<6;i++){
    real theta3=theta2+PI*0.3333333333333*2.;
    real xc=r2*cos(theta2);
    real yc=r2*sin(theta2);
    for(int j=0;j<3;j++){
      xpos[index]=xc+r*cos(theta3);
      ypos[index]=yc+r*sin(theta3);
      theta3+=PI*0.333333333333333333;
      index++;
    };
    theta2+=PI*0.33333333333333;
  };

  xyf.resize(24);
  for(int i=0;i<24;i++){
    xyf[i].resize(10);
  };

  for(int i=0;i<24;i++){
    xyf[i][0]=1.;
    xyf[i][1]=xpos[i];
    xyf[i][2]=ypos[i];
    xyf[i][3]=xpos[i]*xpos[i];
    xyf[i][4]=xpos[i]*ypos[i];
    xyf[i][5]=ypos[i]*ypos[i];
    xyf[i][6]=xpos[i]*xpos[i]*xpos[i];
    xyf[i][7]=xpos[i]*xpos[i]*ypos[i];
    xyf[i][8]=xpos[i]*ypos[i]*ypos[i];
    xyf[i][9]=ypos[i]*ypos[i]*ypos[i];
  };
};

void FRDTAssypow::DoFitting()
{
  for(int i=0;i<24;i++){
    flx[i]=1.;
    if(i>5)flx[i]=0.;
  };

  real amat[10*10];
  real bmat[10];

  for(int i=0;i<10*10;i++){
    amat[i]=0.;
  };
  for(int i=0;i<10;i++){
    bmat[i]=0.;
  };

  for(int y=0;y<10;y++){
    for(int i=0;i<24;i++){
      for(int j=0;j<10;j++){
        amat[y*10+j]+=xyf[i][y]*xyf[i][j];
      };
      bmat[y]+=xyf[i][y]*flx[i];
    };
  };

  gauss(amat,bmat,10,1);

  /*
  for(int i=0;i<24;i++){
    real sol=0.;
    for(int j=0;j<10;j++){
      sol+=bmat[j]*xyf[i][j];
    };
    cout<<sol<<" "<<flx[i]<<"\n";
  };
  */
};

// ++ FRDTAssypow4th

FRDTAssypow4th::FRDTAssypow4th(real ipitch,bool xy)
{
  pitch=ipitch;

  if(xy){

    PutNumc(15);
    PutNump(25);

    real vol=0.5*sqrt(3.)*pitch*pitch;
    real xlen=pitch*0.5;
    real ylen=vol/(xlen*2.)*0.5;

    real yy=-ylen*2;
    int index=0;
    for(int y=0;y<5;y++){
      real xx=-xlen*2;
      for(int x=0;x<5;x++){
        xpos[index]=xx;
        ypos[index]=yy;
        index++;
        xx+=xlen;
      };
      yy+=ylen;
    };

  }else{

    PutNumc(15);
    PutNump(24);

    real r=pitch*0.3333333333333333;
    real theta=-PI*0.33333333333;
    for(int i=0;i<6;i++){
      xpos[i]=r*cos(theta);
      ypos[i]=r*sin(theta);
      theta+=PI*0.333333333333333;
    };

    int index=6;
    real r2=pitch;
    real theta2=-PI*0.333333333333;
    for(int i=0;i<6;i++){
      real theta3=theta2+PI*0.3333333333333*2.;
      real xc=r2*cos(theta2);
      real yc=r2*sin(theta2);
      for(int j=0;j<3;j++){
        xpos[index]=xc+r*cos(theta3);
        ypos[index]=yc+r*sin(theta3);
        theta3+=PI*0.333333333333333333;
        index++;
      };
      theta2+=PI*0.33333333333333;
    };

  };

  xyf.resize(num_p);
  for(int i=0;i<num_p;i++){
    xyf[i].resize(num_c);
  };

  for(int i=0;i<num_p;i++){
    xyf[i][0]=1.;
    xyf[i][1]=xpos[i];
    xyf[i][2]=ypos[i];
    xyf[i][3]=xpos[i]*xpos[i];
    xyf[i][4]=xpos[i]*ypos[i];
    xyf[i][5]=ypos[i]*ypos[i];
    xyf[i][6]=xpos[i]*xpos[i]*xpos[i];
    xyf[i][7]=xpos[i]*xpos[i]*ypos[i];
    xyf[i][8]=xpos[i]*ypos[i]*ypos[i];
    xyf[i][9]=ypos[i]*ypos[i]*ypos[i];
    xyf[i][10]=xpos[i]*xpos[i]*xpos[i]*xpos[i];
    xyf[i][11]=xpos[i]*xpos[i]*xpos[i]*ypos[i];
    xyf[i][12]=xpos[i]*xpos[i]*ypos[i]*ypos[i];
    xyf[i][13]=xpos[i]*ypos[i]*ypos[i]*ypos[i];
    xyf[i][14]=ypos[i]*ypos[i]*ypos[i]*ypos[i];
  };
};

void FRDTAssypow4th::PutNumc(int i)
{
  num_c=i;
  coef.resize(num_c);
};

void FRDTAssypow4th::PutNump(int i)
{
  num_p=i;
  xpos.resize(num_p);
  ypos.resize(num_p);
  flx.resize(num_p);
};

void FRDTAssypow4th::DoFitting()
{
  /*
  for(int i=0;i<num_p;i++){
    flx[i]=1.;
    if(i>5)flx[i]=0.8;
  };
  */

  real amat[num_c*num_c];
  real bmat[num_c];

  for(int i=0;i<num_c*num_c;i++){
    amat[i]=0.;
  };
  for(int i=0;i<num_c;i++){
    bmat[i]=0.;
  };

  for(int y=0;y<num_c;y++){
    for(int i=0;i<num_p;i++){
      for(int j=0;j<num_c;j++){
        amat[y*num_c+j]+=xyf[i][y]*xyf[i][j];
      };
      bmat[y]+=xyf[i][y]*flx[i];
    };
  };

  gauss(amat,bmat,num_c,1);

  for(int i=0;i<num_c;i++){
    coef[i]=bmat[i];
  };
  /*
  for(int i=0;i<num_p;i++){
    real sol=0.;
    for(int j=0;j<num_c;j++){
      sol+=bmat[j]*xyf[i][j];
    };
    cout<<sol<<" "<<flx[i]<<"\n";
  };
  */
};

void FRDTAssypow4th::PutFlux(real *inp)
{
  for(int i=0;i<num_p;i++){
    flx[i]=inp[i];
  };
};

real FRDTAssypow4th::GetValue(real x,real y)
{
  real xyfn[num_c];
  xyfn[0]=1;
  xyfn[1]=x;
  xyfn[2]=y;
  xyfn[3]=x*x;
  xyfn[4]=x*y;
  xyfn[5]=y*y;
  xyfn[6]=x*x*x;
  xyfn[7]=x*x*y;
  xyfn[8]=x*y*y;
  xyfn[9]=y*y*y;
  xyfn[10]=x*x*x*x;
  xyfn[11]=x*x*x*y;
  xyfn[12]=x*x*y*y;
  xyfn[13]=x*y*y*y;
  xyfn[14]=y*y*y*y;

  real ret=0.;
  for(int i=0;i<num_c;i++){
    ret+=xyfn[i]*coef[i];
  };

  return ret;
};

// +++++++++++++++++++++++++++++++++

void GetAssemblyPosition(int layer,real pin_pitch,real *x,real *y)
{
  x[0]=0.;
  y[0]=0.;
  real xtm=0.;
  real ytm=0.;
  real pnp=pin_pitch;

  real vec_x[6];
  real vec_y[6];
  for(int i=0;i<6;i++){
    real theta=PI*0.16666666666+PI*0.333333333333*(i+3);
    vec_x[i]=cos(theta);
    vec_y[i]=sin(theta);
  };

  int idx=1;
  for(int i=0;i<layer;i++){
    ytm=pnp*(i+1);
    for(int j=0;j<6;j++){
      for(int k=0;k<i+1;k++){
        x[idx]=xtm;
        y[idx]=ytm;
        idx++;
        xtm+=vec_x[j]*pnp;
	ytm+=vec_y[j]*pnp;
      };
    };
  };

};

