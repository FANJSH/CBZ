#include <cstdlib>
#include "MulticellBurner.h"

MulticellBurner::MulticellBurner():GeneralBurner()
{
  pl0_calc=false;
  //pl0_calc=true;
  // If [pl0_calc] is true, PL order is automatically set to be zero because of memory reduction.
  // In this case, transport cross section is defined by ([Sigt] - [mu-bar]*[sigel]).
  // If PL order is 1, [Sigtr] = [Sigt] - [Sigs_1], which is more rigorous definition.
  pij_storing=true;

  sensitivity=true;

  corrector_calc=false;
  wpc_direct_calc=false;
  dancoff_input=false;
  mednum_fuel=0;
  mednum_nonfuel=2;

  collapsing=false;
  collapse_dir="";

  aflx_legendre=1;
};

void MulticellBurner::PutMednumFuel(int i)
{
  init_nucnum.resize(i);
  init_nucid.resize(i);
  init_nucden.resize(i);
  init_temperature.resize(i);

  hm_weight_init_per_medium.resize(i);
  acburn_per_medium.resize(i);

  mednum_fuel=i;

  mednum=i+mednum_nonfuel;

  med_clad=1;
  med_water=2;
  med.resize(1+mednum_nonfuel);
};

void MulticellBurner::PutFuelDataMB(int i, int *id, real *den, real temp, int medid)
{
  if(medid>=mednum_fuel){
    cout<<"# Error in MulticellBurner::PutFuelDataMB.\n";
    cout<<"# Medium ID "<<medid<<" is larger than [mednum_fuel].\n";
    exit(0);
  };

  init_nucnum[medid]=i;
  init_temperature[medid]=temp;

  init_nucid[medid].resize(i);
  init_nucden[medid].resize(i);
  for(int j=0;j<i;j++){
    init_nucid[medid][j]=id[j];
    init_nucden[medid][j]=den[j];
  };

  if(medid==0)PutFuelData(i,id,den,temp,medid);
};

void MulticellBurner::AddNonfuelData(int i, int *id, real *den, real temp)
{
  mednum+=1;
  mednum_nonfuel+=1;

  Medium newmed;
  med.push_back(newmed);

  PutNonfuelData(mednum_nonfuel,i,id,den,temp);
};

void MulticellBurner::PutIGI(IrregularGeometryInformation &igi0,IrregularGeometryInformation &igi,bool reflective)
{
  totm=igi.GetRegion();
  fuel_r=igi0.GetCircle(igi0.GetIcir()-1).GetR(); // [cm]

  // +++ TrajectorySet for self-shielding calculation ++++++++++++++++++++++++
  sys.PutBoundaryCondition(Periodic);
  sys.CalTrajectory(igi0,8,0.02,45.);
  //sys.CalTrajectory(igi0,1,0.08,45.);

  // +++ TrajectorySet for neutron flux calculation
  if(reflective){
    sys_f.PutBoundaryCondition(Reflective);
    sys_f.CalTrajectory(igi,64,0.02,360.); // Default
    //sys_f.CalTrajectory(igi,16,0.05,360.);
    //sys_f.CalTrajectory(igi,16,0.08,360.);
    //sys_f.CalTrajectory(igi,8,0.1,360.);            
    // (for NEL2020_3x3)
    //sys_f2.PutBoundaryCondition(Reflective);
    //sys_f2.CalTrajectory(igi,16,0.05,360.);
  }else{
    sys_f.PutBoundaryCondition(Periodic);
    //sys_f.CalTrajectory(igi,1,0.08,45.);
    //sys_f.CalTrajectory(igi,5,0.03,45.);
    sys_f.CalTrajectory(igi,8,0.02,45.);
    //sys_f.CalTrajectory(igi,16,0.02,45.);
    // (for NEL2020_3x3)
    //sys_f2.PutBoundaryCondition(Periodic);
    //sys_f2.CalTrajectory(igi,2,0.05,45.);
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

void MulticellBurner::PutIGIQuarter(IrregularGeometryInformation &igi0,IrregularGeometryInformation &igi)
{
  // ... Quarter symmetric, but full system

  totm=igi.GetRegion();
  fuel_r=igi0.GetCircle(igi0.GetIcir()-1).GetR(); // [cm]

  // +++ TrajectorySet for self-shielding calculation ++++++++++++++++++++++++
  sys.PutBoundaryCondition(Periodic);
  sys.CalTrajectory(igi0,8,0.02,45.);
  //sys.CalTrajectory(igi0,1,0.08,45.);

  // +++ TrajectorySet for neutron flux calculation  
  sys_f.PutBoundaryCondition(Periodic);
  sys_f.CalTrajectory(igi,16,0.02,90.);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

void MulticellBurner::PutIGI_NEL2020_3x3(IrregularGeometryInformation &igi0,IrregularGeometryInformation &igi,int num_divAng,real dist_trac,bool reflective)
{
  totm=igi.GetRegion();
  fuel_r=igi0.GetCircle(igi0.GetIcir()-1).GetR(); // [cm]

  // +++ TrajectorySet for self-shielding calculation ++++++++++++++++++++++++
  sys.PutBoundaryCondition(Periodic);
  sys.CalTrajectory(igi0,8,0.02,45.);
  //sys.CalTrajectory(igi0,1,0.08,45.);

  // +++ TrajectorySet for neutron flux calculation
  if(reflective){
    sys_f.PutBoundaryCondition(Reflective);
    sys_f.CalTrajectory(igi,64,0.02,360.); // Default
    //sys_f.CalTrajectory(igi,16,0.05,360.);
    //sys_f.CalTrajectory(igi,16,0.08,360.);
    //sys_f.CalTrajectory(igi,8,0.1,360.);            
    // (for NEL2020_3x3)
    sys_f2.PutBoundaryCondition(Reflective);
    //sys_f2.CalTrajectory(igi,64,0.02,360.);
    //sys_f2.CalTrajectory(igi,16,0.05,360.);
    sys_f2.CalTrajectory(igi,num_divAng,dist_trac,360.);        
  }else{
    sys_f.PutBoundaryCondition(Periodic);
    sys_f.CalTrajectory(igi,8,0.02,45.); // Default    
    //sys_f.CalTrajectory(igi,2,0.05,45.);
    //sys_f.CalTrajectory(igi,5,0.03,45.);

    //sys_f.CalTrajectory(igi,16,0.02,45.);
    // (for NEL2020_3x3)
    sys_f2.PutBoundaryCondition(Periodic);
    //sys_f2.CalTrajectory(igi,8,0.02,45.);
    //sys_f2.CalTrajectory(igi,2,0.05,45.);
    sys_f2.CalTrajectory(igi,num_divAng,dist_trac,45.);        
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

void MulticellBurner::PutIGI_NEL2020(IrregularGeometryInformation &igi_33,int num_divAng,real dist_trac)
{
  sys_f2.PutBoundaryCondition(Reflective);
  //sys_f2.CalTrajectory(igi_33,16,0.05,360.);
  sys_f2.CalTrajectory(igi_33,num_divAng,dist_trac,360.);  
};

void MulticellBurner::PutRelationRegionMedium(int *inp)
{
  region_medium.resize(totm);
  for(int i=0;i<totm;i++){
    region_medium[i]=inp[i];
  };
};

void MulticellBurner::PutRelationRegionMedium(vector<int> &inp)
{
  region_medium.resize(totm);
  for(int i=0;i<totm;i++){
    region_medium[i]=inp[i];
  };
};

void MulticellBurner::PutNuclideDataToMedium(GroupData1D &den, int medid_med)
{
  int nucn=med[medid_med].GetNucnum();
  for(int i=0;i<nucn;i++){
    med[medid_med].GetNuclideInTurn(i).PutDensity(den.get_dat(i));
  };
};

void MulticellBurner::PutMicroXSDataToMedium(GroupData1D &sf,GroupData1D &sc,GroupData1D &s2n,int nuc,int medid_med,int nuclide_info)
{
  if(nuclide_info==1)med[medid_med].GetNuclideInTurn(nuc).GetMicxs().GetData1d(sigf).copy(sf);
  med[medid_med].GetNuclideInTurn(nuc).GetMicxs().GetData1d(sigc).copy(sc);
  med[medid_med].GetNuclideInTurn(nuc).GetMicxs().GetData1d(sign2n).copy(s2n);
};

void MulticellBurner::PreCalculation(Burnup &bu)
{
  bu.AddNuclideToMediumFromBurnupChain(med[0]);
  nucn=med[0].GetNucnum();

  // +++ Initial selfshielding calculation for non-fuel medium
  opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[med_clad],xslib);

  /*
  opc.GiveInfiniteDillutionCrossSection(med[med_water],xslib);
  opc.CalThermalScatteringMatrix(med[med_water],xslib,3.93);
  med[med_water].CalSigtr(0);
  */
  for(int i=2;i<=mednum_nonfuel;i++){
    opc.GiveInfiniteDillutionCrossSection(med[i],xslib);
    //opc.CalThermalScatteringMatrix(med[i],xslib,3.93);
    opc.CalThermalScatteringMatrix(med[i],xslib,4.048); // ESAB from MVP library    
    med[i].CalSigtr(0);
  };
  //med[med_water].GetMacxs().GetData1d(sigtr).show_self(); exit(0);
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Medium-wise volume calculation
  vol_med.resize(mednum,0.);
  for(int i=0;i<sys_f.GetRegnum();i++){
    vol_med[region_medium[i]]+=sys_f.GetVolume(i);
  };
  cout<<"#\n# Medium-wise total volume\n";
  for(int i=0;i<mednum;i++){
    cout<<"# "<<i<<" "<<vol_med[i]<<"\n";
  };

  // +++ nuclide info
  nuclide_info.resize(nucn,0); // 0:no cross section, 1:fissile, 2:other
  for(int i=0;i<nucn;i++){
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

};

void MulticellBurner::Calculation(Burnup &bu, int med_target, bool adjoint)
{
  PreCalculation(bu);
  ForwardCalculation(bu,med_target,adjoint);
};

void MulticellBurner::CalculationADTS_toOWPC(Burnup &bu, int med_target, int ngrp, int *bgrp)
{
  PreCalculation(bu);
  ForwardCalculationADTS_toOWPC(bu,med_target,ngrp,bgrp);
};

void MulticellBurner::CalculationNEL2020_3x3(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  PreCalculation(bu);
  // ForwardCalculationNEL2020_3x3(bu,med_target,sm_step,ngrp,bgrp);   
  // ForwardCalculationADTS_toOWPC(bu,med_target,ngrp,bgrp);
  ForwardCalculationADTS_SamSys(bu,med_target,sm_step,ngrp,bgrp);   
};



void MulticellBurner::CalculationNEL2020(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  PreCalculation(bu);
  //ForwardCalculationNEL2020(bu,med_target);
  //ForwardCalculationNEL2020Sasuga(bu,med_target);
  ForwardCalculationNEL2020SasugaFinal(bu,med_target,sm_step,ngrp,bgrp);
  // ForwardCalculationADTS_DifSys(bu,med_target,sm_step,ngrp,bgrp);
};

void MulticellBurner::ForwardCalculation(Burnup &bu, int med_target, bool adjoint)
{
  // This source code was finally prepared in 2021/5/26 by Chiba based on the Sasuga-kun preparing version.
  
  bool owpc_corr=true; // If [false], the (log-averaging) PC is used.
  
  // ... Hard-coded parameters for weighted predictor-corrector
  //real wgt_nc=1.2; // relative weight for corrector
  real wgt_nc=1.0; // relative weight for corrector

  // ... OWPC to store positions of Gd-155 and 157
  int pos_gd155, pos_gd157;
  for(int j=0;j<nucn;j++){
    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    if(nucid==641550)pos_gd155=j;
    if(nucid==641570)pos_gd157=j;
  };

  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculation.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting for forward calculation +++

  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  flx_med.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  // +++ Array setting for predictor-corrector calculation
  //
  // - Generally multi-group microscopic cross section data at every burn steps
  //   are NOT stored because of their large required memory, but those are 
  //   required for GPT (adjoint) calculations, so those at every burnup steps
  //   are stored in the array [mic_sigx].
 
  {
  int tmp=1;
  if(adjoint)tmp=burn_step+1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };

  vector< vector<GroupData1D> > mic_sigf_c;
  vector< vector<GroupData1D> > mic_sigc_c;
  vector< vector<GroupData1D> > mic_sign2n_c;
  vector<GroupData1D> flx_med_c;
  if(corrector_calc){
    //int tmp=1;
    //if(adjoint)tmp=burn_step+1;
    int tmp=burn_step+1;
    xsc_1g_p.resize(tmp);
    xsn2n_1g_p.resize(tmp);
    xsf_1g_p.resize(tmp);
    for(int j=0;j<tmp;j++){
      xsc_1g_p[j].resize(mednum_fuel);
      xsn2n_1g_p[j].resize(mednum_fuel);
      xsf_1g_p[j].resize(mednum_fuel);
      for(int i=0;i<mednum_fuel;i++){
        xsc_1g_p[j][i].resize(nucn,0.);
        xsn2n_1g_p[j][i].resize(nucn,0.);
        xsf_1g_p[j][i].resize(nucn,0.);
      };
    };

    mic_sigf_c.resize(mednum_fuel); 
    mic_sigc_c.resize(mednum_fuel); 
    mic_sign2n_c.resize(mednum_fuel);
    for(int i=0;i<mednum_fuel;i++){
      mic_sigf_c[i].resize(nucn);
      mic_sigc_c[i].resize(nucn);
      mic_sign2n_c[i].resize(nucn);
    };

    flx_med_c.resize(mednum);
    for(int i=0;i<mednum;i++){
      flx_med_c[i].put_imax(group);
    };

    total_flux_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      total_flux_p[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        total_flux_p[i][j].resize(mednum_fuel);
      };
    };

  };

  // ... OWPC
  vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC (Rp)
  vector<real> np_gd5, np_gd7; // Np

  vector<real>  old_n0155;
  vector<real>  old_n0157;
  vector<real>  old_r0155;
  vector<real>  old_r0157;
  vector<real>  old_np155;
  vector<real>  old_rp155;
  vector<real>  old_rp157;
  vector<real>  old_rc155;
  vector<real>  old_rc157;
  vector<real>  old_n0155_2;
  vector<real>  old_r0155_2;
  vector<real>  old_r0157_2;
  vector<real>  old_x0155;
  vector<real>  old_x0157;
  real n0_xx[2][mednum_fuel];
  real np_xx[2][mednum_fuel];
  real rp_xx[2][mednum_fuel];
  real rc_xx[2][mednum_fuel];

  
  real old_xsc5[mednum_fuel];
  real old_xsc7[mednum_fuel];
  real old_xsc5_2[mednum_fuel];
  real old_xsc7_2[mednum_fuel];
  real old_xscc_5[mednum_fuel];
  real old_xscc_7[mednum_fuel];
  real old_rc_5[mednum_fuel];
  real old_rc_7[mednum_fuel];
  real tmp_5[mednum_fuel];
  real tmp_7[mednum_fuel];

  real alpha_5[mednum_fuel];
  real alpha_7[mednum_fuel];
  real time;
  real time_be;
  real time_be_2;
  
  //sasuga addition
  old_n0155.resize(mednum_fuel);
  old_n0157.resize(mednum_fuel);
  old_r0155.resize(mednum_fuel);
  old_r0157.resize(mednum_fuel);
  old_np155.resize(mednum_fuel);
  old_rp155.resize(mednum_fuel);
  old_rp157.resize(mednum_fuel);
  old_rc155.resize(mednum_fuel);
  old_rc157.resize(mednum_fuel);
  old_n0155_2.resize(mednum_fuel);
  old_r0155_2.resize(mednum_fuel);
  old_r0157_2.resize(mednum_fuel);
  old_x0155.resize(mednum_fuel);
  old_x0157.resize(mednum_fuel);

  if(corrector_calc){
    rr_gd5.resize(mednum_fuel);
    rr_gd7.resize(mednum_fuel);
    np_gd5.resize(mednum_fuel);
    np_gd7.resize(mednum_fuel);
  };

  // +++ Array setting for adjoint calculation +++
  //
  //  [fwd_nuc_int] is a number density at time-mesh-center point.
  //  Time-averaged number density during one burnup step is calculated 
  //  from those at beginning, center and end of this burnup step
  if(adjoint){
    fwd_nuc_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
  	  fwd_nuc_int[i][j][k].put_imax(nucn);
	};
      };
    };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();

  // +++ Forward burn-up calculation +++
  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  // +++ Eigenvalue calculation
  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);  

  if(adjoint){

    // +++ array setting for angular flux storing in adjoint calculations.
    int totsn=lat.GetQuad().GetSN();
    // Quadrature setting
    quad.Initialize(2,aflx_legendre);
    quad.PutSN(totsn);
    real *mui=new real[totsn];
    real *eai=new real[totsn];
    real *xii=new real[totsn];
    real *w=new real[totsn];
    for(int i=0;i<totsn;i++){
      mui[i]=lat.GetQuad().GetMu(i);
      eai[i]=lat.GetQuad().GetEata(i);
      xii[i]=lat.GetQuad().GetXi(i);
      w[i]=lat.GetQuad().GetOmega(i);
    };
    quad.PutData(mui,eai,xii,w);
    quad.CalValue();

    if(aflx_legendre==-1){
      volaflx_mesh.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };

  };

  // +++ GPT-PC +++++++++++++++++++++++++++++++++++++++++++++++++
  if(adjoint&&corrector_calc){
    fwd_nuc_p.resize(burn_step+1);
    fwd_nuc_p_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_p[i].resize(sub_step+1);
      fwd_nuc_p_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_p[i][j].resize(mednum_fuel);
        fwd_nuc_p_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
          fwd_nuc_p[i][j][k].put_imax(nucn);
          fwd_nuc_p_int[i][j][k].put_imax(nucn);
        };
      };
    };
    volflx_mesh_p.resize(burn_step);  
    power_factor_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      volflx_mesh_p[i].resize(totm);
      for(int j=0;j<totm;j++){
        if(region_medium[j]<mednum_fuel){
          volflx_mesh_p[i][j].put_imax(group);
        };
      };
      int sub_step=sub_step_list[i];
      power_factor_p[i].resize(sub_step+1);
    };
    int totsn=lat.GetQuad().GetSN();

    if(aflx_legendre==-1){
      volaflx_mesh_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh_p[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh_p[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl_p[i].resize(totm);
        for(int j=0;j<totm;j++){
    	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl_p[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real abs_frac[burn_step+1];
  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;
    if(adjoint)bstmp=st;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    //lat.MedClear();
    for(int i=0;i<mednum_fuel;i++){
      // +++ Self-shielding calculation
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library      

      /*
      // ... neutron flux energy spectrum calculated at the last step is used
      //     to calculate macroscopic fission spectrum
      if(st!=0){
	GroupData1D aveflx;
	aveflx.set_zero();
        for(int ii=0;ii<totm;ii++){
          int medid=region_medium[ii];
          if(medid==i)aveflx=aveflx+volflx_mesh[st-1][ii];
        };
	med[0].GetFlux().copy(aveflx);
      };
      */

      med[0].CalMacroFromMicro();

      /*
      cout<<"# Medium :: "<<i<<"\n";
      med[0].ShowMacroXS1D();
      */

      /*
      if(i==40){
	for(int k=0;k<nucn;k++){
	  cout<<k<<" "<<med[0].GetNuclideInTurn(k).GetMatnum()<<" "<<fwd_nuc[st][0][i].get_dat(k)<<"\n";
	};
      };
      */

      // +++ Macro&Micro data storing
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };
      if(st!=burn_step&&adjoint){
        macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
      };

      if(st==0){
        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();
      }else{
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
      };
    };
    if(st==0){
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };
    }else{
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };
    };
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);

    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    if(adjoint)lat.PutWriteFlux();
    lat.Print();
    keff[st]=lat.CalIgen();

    //cout<<"# Total inner iteration : "<<lat.GetTotalInnerIterationSum()<<"\n"; exit(0);

    // --- (Generating medium instances include collapsing cross section data) 
    if(collapsing){

    for(int i=0;i<mednum_fuel;i++){
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      //med[0].CalMacroFromMicro();

      int ngrp=9;
      int bgrp[]={21,46,91,115,125,134,151,159,171};
      GroupData1D flx=lat.GetIntegratedFlux(i);
      Medium bmed=med[0].Cond(ngrp,bgrp,flx,flx,true);
      bmed.GetMacxs().GetData1d(d).copy(bmed.GetFlux());
      // Neutron flux is stored as diffusion coefficient (temporal treatment)
      string filename="bmed_med"+IntToString(i)+"_st"+IntToString(st);
      bmed.WriteFile(collapse_dir,filename,true);
    };

    };
    // --------------------------------------------------------------------------------


    // ------------------------------

    real sum_tot=0.;
    real sum_fmed=0.;
    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);

      real tmp=flx*lat.GetMedium(i).GetMacxs().GetData1d(siga);
      sum_tot+=tmp;
      if(i<mednum_fuel)sum_fmed+=tmp;

      flx_med[i]=flx*(1./vol_med[i]);
      // +++ One-group cross section storing (at predictor step)
      if(i<mednum_fuel){
        for(int j=0;j<nucn;j++){
          if(nuclide_info[j]!=0){
            if(nuclide_info[j]==1){
   	      xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	    };
	    xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	    xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
          };
	};
      };
    };
    abs_frac[st]=sum_fmed/sum_tot; // absorption rate fraction in fuel media
    cout<<"# Absorption fraction : "<<abs_frac[st]<<"\n";
  
    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
          if(adjoint){
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
  	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };
	  };
	};
      };

      // +++ Burnup calculation
      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
          //power_org=CalculationPinPower(bu,st,j,med_target,sumflx*vol_med[med_target]);
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        accumulated_day+=burn_span;
	//total_flux[st][j].resize(mednum_fuel);
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //   CORRECTOR CALCULATION
      if(corrector_calc){

	cout<<"# Corrector calculation ...\n";

	// (OWPC)	
	for(int i=0;i<mednum_fuel;i++){
	  real tmp=total_flux[st][0][i];
  	  rr_gd5[i]=tmp*xsc_1g[st][i][pos_gd155];  // Reaction rate at BOC (Rp)
  	  rr_gd7[i]=tmp*xsc_1g[st][i][pos_gd157];
	};

	// ++++++++++++++++++++++++

	// (Predictor calculation results are stored in the array of [XXX_p].)
        swap(total_flux_p[st], total_flux[st]);
        swap(xsc_1g_p[st], xsc_1g[st]);
        swap(xsn2n_1g_p[st], xsn2n_1g[st]);
        swap(xsf_1g_p[st], xsf_1g[st]);

	if(adjoint){
          swap(keff_p[st], keff[st]);
          swap(fwd_nuc_p[st], fwd_nuc[st]);
          swap(fwd_nuc_p_int[st], fwd_nuc_int[st]);
          swap(power_factor_p[st], power_factor[st]);
	  swap(volflx_mesh_p[st], volflx_mesh[st]);
	  if(aflx_legendre==-1){
            swap(volaflx_mesh_p[st], volaflx_mesh[st]);
	  }else{
  	    swap(volaflx_pl_p[st], volaflx_pl[st]);
	  };
	  swap(macxs_p[st], macxs[st]);
	  fwd_nuc[st][0]=fwd_nuc_p[st][0];
	};

      for(int i=0;i<mednum_fuel;i++){
        // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
        PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
  	  };
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        }else{
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        };
        //opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library	 
        med[0].CalMacroFromMicro();
        //opc.CalFissionSpectrumMatrix(med[0],xslib);

        // +++ Macro&Micro data storing
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
	if(adjoint){
          macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
	};
      };

      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };

      lat.PutThermalIteration(3);
      lat.PutPL(0);
      lat.NoCMRAcceleration();
      real keff_corr=lat.CalIgen();
      
      // ------------------------------

      if(adjoint)keff[st]=keff_corr;

      for(int i=0;i<mednum;i++){
        GroupData1D flx=lat.GetIntegratedFlux(i);
        flx_med_c[i]=flx*(1./vol_med[i]);
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1)xsf_1g[st][i][j]=mic_sigf_c[i][j].Cond(flx);
	      xsc_1g[st][i][j]=mic_sigc_c[i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n_c[i][j].Cond(flx);
	    };
          };
	};
      };

      if(adjoint){
        for(int i=0;i<totm;i++){
          if(region_medium[i]<mednum_fuel){
            real vol=lat.GetMesh(i).GetVolume();
            volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol;
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };

	  };
        };


      };

      real acburn_pre=accumulated_burn-acburn[st];
      vector<real> acburn_pre_per_medium(mednum_fuel);
      for(int j=0;j<mednum_fuel;j++){
	acburn_pre_per_medium[j]=accumulated_burn_per_medium[j]-acburn_per_medium[j][st];
      };
      // accumulated burnup calculated by the predictor step

      
      // 
      // .....  One-group cross section correction for Rc in A-OWPC .....
      //

      
      if(owpc_corr){
	
      for(int i=0;i<mednum_fuel;i++){

	real n0=fwd_nuc[st][0][i].get_dat(pos_gd157);   // initial 
        real np=fwd_nuc[st+1][0][i].get_dat(pos_gd157); // predictor results      
        if((n0-np)/np>1e-2){
          // This if branch is added in 2021/12/4 since negative cross section is detected in the Gd cross section.
	  // This should be consistent in the following process in the actual correction to the final ND.
	
	real xsc0_5=xsc_1g_p[st][i][pos_gd155]; //rp
	real xsc0_7=xsc_1g_p[st][i][pos_gd157]; //rp
	real xscc_5=xsc_1g[st][i][pos_gd155];//Rcに相当するXc
	real xscc_7=xsc_1g[st][i][pos_gd157];
	if(st>0){
	  real n0_5=fwd_nuc[st][0][i].get_dat(pos_gd155);
	  real np_5=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	  real t_5=old_xsc5[i]+(xsc0_5-old_xsc5[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);//前のステップでの正しいXcを計算
	  real t_7=old_xsc7[i]+(xsc0_7-old_xsc7[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
	  if(st>1){
	    real a=n0_5, b=old_n0155[i], c=old_n0155_2[i];
	    real aa_5=xsc0_5, bb_5=old_xsc5[i], cc_5=old_xsc5_2[i];
	    real aa_7=xsc0_7, bb_7=old_xsc7[i], cc_7=old_xsc7_2[i];
	    real X_5=old_np155[i];
	    t_5=aa_5*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_5*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_5*(X_5-a)*(X_5-b)/(c-a)/(c-b);//前のステップでの正しいXcを計算
	    t_7=aa_7*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_7*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_7*(X_5-a)*(X_5-b)/(c-a)/(c-b);
	  };

	  // correction without any consideration on time step difference
	  /*
	  xsc_1g[st][i][pos_gd155]=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5);
	  xsc_1g[st][i][pos_gd157]=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7);
	  */
	  /*
          // correction with considering the time step length difference
	  xsc_1g[st][i][pos_gd155]=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)/burn_time[st-1]*burn_time[st];
	  xsc_1g[st][i][pos_gd157]=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)/burn_time[st-1]*burn_time[st];
	  */

	  // correction with considering time step length x power density (= burnup) 
	  real tmp1=power_density_list[st-1]*burn_time[st-1];
	  real tmp2=power_density_list[st]*burn_time[st];
	  real factor=tmp2/tmp1;
	  real xsc155=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)*factor;
	  real xsc157=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)*factor;
	  if(factor<10.){ // if the burnup length is too different, the correction is NOT applied.
	    if(xsc155>0.)xsc_1g[st][i][pos_gd155]=xsc155;
	    if(xsc157>0.)xsc_1g[st][i][pos_gd157]=xsc157;
	  };
	  
	  old_xsc5_2[i]=old_xsc5[i];
	  old_xsc7_2[i]=old_xsc7[i];
	};
	old_xsc5[i]=xsc0_5;
	old_xsc7[i]=xsc0_7;
	old_xscc_5[i]=xscc_5;
	old_xscc_7[i]=xscc_7;
      };

      }; // end of [if((n0-np)/np>1e-2)]
      }; // end of [if(owpc_corr)]
      // ............................................................................
      
      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med_c[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med_c[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          //accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        //accumulated_day+=burn_span;
        //delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
	  //CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.

	  /*
	  if(i==40){
	    cout<<"# Medium : "<<i<<"\n";
	    for(int j=0;j<nucn;j++){
              cout<<j<<" "<<med[0].GetNuclideInTurn(j).GetMatnum()<<" "<<xsc_1g[st][i][j]<<" "<<xsn2n_1g[st][i][j]<<"\n";
	    };
	    cout<<" total_flux : "<<total_flux[st][j][i]<<"\n";
	    cout<<" delt       : "<<delt[st][j]<<"\n";	    
	  };
	  */

	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	};
	
      }; // end of sub-step
      
      // (OWPC)
      real dt=burn_time[st]*60*60*24;
      real time_mesh=10000;
      real dt_r=dt/time_mesh;
      real omega_155, omega_157; // various correlation conditions

      for(int i=0;i<mednum_fuel;i++){

	  // OWPC treatment for Gd-155 and -157 with various correlation conditions
	if(owpc_corr){

	  if(st>0){ // ... quadratic model
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	    real x0_155=xsc_1g_p[st][i][pos_gd155];//rpに対応する断面積
	    real x0_157=xsc_1g_p[st][i][pos_gd157];//rpに対応する断面積
	    real xp_155=xsc_1g[st][i][pos_gd155]; //rcに対応する断面積
	    real xp_157=xsc_1g[st][i][pos_gd157]; //rcに対応する断面積
	    real n0_r_155;
	    real n0_r_157;
	    
	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt);//Not consider 154,156
	    	    
	    //quadratic model---------------------------
	    real a_5,b_5,c_5,a_7,b_7,c_7;
	    real aa_5,bb_5,cc_5,aa_7,bb_7,cc_7;
	    a_5=n0_155;
	    b_5=np_155;
	    c_5=old_n0155[i];
	    a_7=n0_157;
	    b_7=np_157;
	    c_7=old_n0157[i];
	    real X_5=np_p_155;
	    real X_7=np_p_157;

	    //-------------------------reaction rate base-------------------------------------
	    aa_5=r0_155;
	    bb_5=rp_155;
	    //cc_5=old_r0155[i];
	    cc_5=old_r0155[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba    
	    
	    aa_7=r0_157;
	    bb_7=rp_157;
	    //cc_7=old_r0157[i];
	    cc_7=old_r0157[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba

	    real rc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    real rc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    //----------------------------------------------------------------------------------------

	    /*
	    //---------------------------------------cross section base-------------------------------------
	    aa_5=x0_155;
	    bb_5=xp_155;
	    cc_5=old_x0155[i];
	    
	    aa_7=x0_157;
	    bb_7=xp_157;
	    cc_7=old_x0157[i];
	        
	    real xc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //xc in toyproblem
	    real xc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    //----------------------------------------------------------------
	    */
	    
	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	    real xp_r_155=x0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	    real xp_r_157=x0_157;

	    //------------------------------------start toy problem calculation----------------------------
	    for(int k=0;k<time_mesh;k++){


	      //-----reaction rate base-------
	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      X_5=np_r_155;
	      X_7=np_r_157;
	      
	      real rc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	      real rc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r); //Not consider 154,156
	    
	      rp_r_155=rc_r_155;
	      rp_r_157=rc_r_157;
	      //----------------------------------------

	      /*
	      //-------cross section base-------
	      real np_r_155=n0_r_155*exp(-xp_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-xp_r_157*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156

	      X_5=np_r_155;
	      X_7=np_r_157;

	      real xc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	      real xc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

	      real x_r_155=(xp_r_155+xc_r_155)*0.5;
	      real x_r_157=(xp_r_157+xc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-x_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-x_r_157*total_flux_p[st][0][i]*1e-24*dt_r); //Not consider 154,156
	    
	      xp_r_155=xc_r_155;
	      xp_r_157=xc_r_157;
	      //----------------------------------------
	      */
	      
	    };
	    //-----------------------------------------end toy problem calculation------------------------------------------
	    
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

	    //-------reaction rate base----------
	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	    //------------------------------------------

	    /*
	    //-------cross section base----------
	    omega_155=(R_reference_155-x0_155*total_flux_p[st][0][i])/(xc_155*total_flux_p[st][0][i]-x0_155*total_flux_p[st][0][i]);
	    omega_157=(R_reference_157-x0_157*total_flux_p[st][0][i])/(xc_157*total_flux_p[st][0][i]-x0_157*total_flux_p[st][0][i]);
	    //--------------------------------------
	    */
	    
	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;
	    
	    old_n0155_2[i]=old_n0155[i];
	    old_r0155_2[i]=old_r0155[i];
	    old_r0157_2[i]=old_r0157[i];
	    
	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;
	    old_x0155[i]=x0_155;
	    old_x0157[i]=x0_157;
	    
	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;
	      	    
	  }else{ // ... linear model for [st]==0
	  
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	
	    real n0_r_155;
	    real n0_r_157;

	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt); //Not consider 154,156
	  
	    // (Correlation to Gd-155 number density)
	    real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	    real alpha_157=(rp_157-r0_157)/(np_155-n0_155);    
	    
	    // (Correlation to Gd-155 number density)	
	    real rc_155=r0_155+(np_p_155-n0_155)*alpha_155; // Rc in toy-problem
	    real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	   
	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	  
	    for(int k=0;k<time_mesh;k++){

	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      // (Correlation to Gd-155 number density)
	      real rc_r_155=r0_155+(np_r_155-n0_155)*alpha_155;
	      real rc_r_157=r0_157+(np_r_155-n0_155)*alpha_157;
	
	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r);//Not consider 154,156

	      rp_r_155=rc_r_155;
	      rp_r_157=rc_r_157;

	    };
	        
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);

	    //cout<<"#   Omega_155 / Omega_157 : "<<omega_155<<" "<<omega_157<<"\n";
	    
	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;

	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;

	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;

	  }; //

	}; // end of [if(owpc_corr)]
	  
	for(int j=0;j<nucn;j++){
	
	  real np=fwd_nuc[st+1][0][i].get_dat(j);      // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
	  real n_next=exp((log(np)+log(nc))*0.5);
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();

	  if(nucid==641550||nucid==641570){

  	    real n0=fwd_nuc[st][0][i].get_dat(j);
            if((n0-np)/np>1e-2){

	      if(owpc_corr){
		if(nucid==641550){
		  n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
		}else if(nucid==641570){
		  n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
    	        };
	      }else{
		/*
		// -- OWPC w/o (N,R) correlation ----------------------- 
		real r0=rr_gd5[i]; //Rp
		if(nucid==641570)r0=rr_gd7[i];

		real rp=xsc_1g[st][i][j]*total_flux[st][0][i]; //Rc
		real alpha=(rp-r0)/(np-n0);
		  
		real np_p=n0*exp(-r0*1e-24*dt);   // Nc in toy-problem
		real rc=alpha*np_p+(r0-n0*alpha); // Rc in toy-problem
	    
		real n0_r=n0;
		real rp_r=r0;
		for(int k=0;k<time_mesh;k++){
		  real np_r=n0_r*exp(-rp_r*1e-24*dt_r);
		  real rc_r=alpha*np_r+(r0-n0*alpha);
		  real r_r=(rp_r+rc_r)/2;
		  n0_r=n0_r*exp(-r_r*1e-24*dt_r);
		  rp_r=rc_r;
		};
		
		real R_reference=(log(n0)-log(n0_r))/dt*1e24;
		real omega=(R_reference-r0)/(rc-r0);
		  
		// -- PPC -------------------------------
		real R_p=-log(np/n0)/dt;
		real R_c=-log(nc/n0)/dt;
		real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
		R_reference=(R_p+R_c_c)/2;
		omega=(R_reference-R_p)/(R_c-R_p);
		// --------------------------------

		n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応
		*/

 	        // ... conventional PC ...
                //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
                //n_next=(np+nc)*0.5;  // Conventional PC (in linear)
                //n_next=nc;          
              };

	    };

	  };

	  /*
          if(!owpc_corr){
	    n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
	    cout<<nucid<<" "<<n_next<<" "<<np<<" "<<nc<<" "<<"\n";
	  };
	  */

          if(adjoint)        n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
          if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC

	  // (For OWPC)
    	  if(nucid==641550)np_gd5[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
  	  if(nucid==641570)np_gd7[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd157);

	  /*
	  if(i==40){
	  cout<<"# Medium : "<<i<<"\n";
  	    cout<<j<<" "<<med[0].GetNuclideInTurn(j).GetMatnum()<<" "<<n_next<<" "<<np<<" "<<nc<<"\n";	     
	  };
	  */
	  
          fwd_nuc[st+1][0][i].put_data(j,n_next);

	}; // end of nuclides iteration

      }; // end of medium iteration

      // Adjustment of accumulated burn
      for(int i=0;i<mednum_fuel;i++){
        real acburn_cor=accumulated_burn_per_medium[i]-acburn_per_medium[i][st]-acburn_pre_per_medium[i];
        accumulated_burn_per_medium[i]=acburn_per_medium[i][st]+(acburn_pre_per_medium[i]+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };

      if(input_flux_level){
        real acburn_cor=accumulated_burn-acburn[st]-acburn_pre;
        accumulated_burn=acburn[st]+(acburn_pre+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };
      cout<<"#     ... terminated.\n";
      }; // END OF CORRECTOR CALCULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    }; // end part of [If (st!=burn_step) ]
 
  }; // loop-end of burnup step

};

void MulticellBurner::ForwardCalculationNEL2020Sasuga(Burnup &bu, int med_target)
{
  // Note that this is hard-coded for VERA 2O benchmark
  
  // Pre-determined parameters/constants
  ofstream fout1("out1_2.dat");
  ofstream fout2("out2_2.dat");

  int num_33_system=2;
  // The number of Gd pins in the system
  // = The number of the 3x3 multicells treated in ADTS
  
  int gd_id_list[]={16,40};
  // The first medium ID of Gd pins.
  // Medium IDs [16-23] and [40-47] are used for Gd pins.
  
  int uo2_id_list[]={24,48};
  // Medium ID for UO2 considered in the 3x3 multicell
  // ID24 and ID48 are at the right-neighboring position of the Gd pins.

  // ... [Region ID]-[Medium ID] relation in the 3x3 multicell ...
  vector<int> region_medium_33;
  int region_medium_num=9+12+8*4*2; // The number of total regions in the 3x3 multicell
  int region_medium_33_inp[]={ // Medium ID for each region
    10,10,10,10,10, 10,10,10,10, // background meshes(9)
    10,10,10,9,0,1,2,3,4,5,6,7,  // center Gd-pin(12)
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // right-neighboring pin
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // top-neighboring pin    
  };
  for(int i=0;i<region_medium_num;i++){
    region_medium_33.push_back(region_medium_33_inp[i]);
  };
  
  
  
  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  vector<GroupData1D> dflux_1; // absolute difference in medium-wise neutron flux for the first 3x3 multicell
  vector<GroupData1D> dflux_2; // absolute difference in medium-wise neutron flux for the second 3x3 multicell
  flx_med.resize(mednum);
  dflux_1.resize(mednum);
  dflux_2.resize(mednum);
  
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };
    
    
  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  // !!! sasuga addition
  double dxsc_5[num_33_system][8];//Gd-155の断面積の差を保存しておく, [8] is the number of medium in Gd-pin
  double dxsc_7[num_33_system][8];//Gd-155の断面積の差を保存しておく
  int st_10; // Time step at which the neutron transport equation of the original system is solved.
  double factor_3_3;
  double factor_as;
  // !!!

  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";


    // +++ SOLVING NEUTRON TRANSPORT EQUATION +++++++++++++++++++++++++++++++++++++++

    // +++ Instance generation of MECSystem to solve neutron transport equation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);  

    // +++ fuel medium data preparation
    for(int i=0;i<mednum_fuel;i++){

      // (Self-shielding calculation)
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
      };
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib, 4.048); // ESAB from MVP library     
      med[0].CalMacroFromMicro();

      // (Macro&Micro data storing)
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      lat.AddMedium(med[0]); 
      lat.GetMedium(i).NuclideClear();

    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.Print();

    int pp=200; //sasuga addition

    // ------------------------------------------------------------------------------
    if((st<50)||((st-50)%pp==0)){ //sasuga addition
      //if(st==0){ //sasuga addition

      // ... Neutron transport equation for the original system
      //     is solved at this time step point
      
      keff[st]=lat.CalIgen(); 

      for(int i=0;i<mednum;i++){
	GroupData1D flx=lat.GetIntegratedFlux(i);
	flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
	// +++ One-group cross section storing
	if(i<mednum_fuel){
	  for(int j=0;j<nucn;j++){
	    if(nuclide_info[j]!=0){
	      if(nuclide_info[j]==1){
		xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	      };
	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
	    };
	  };
	};
      };
      
    }else{

      // ... Neutron transport equation is NOT solved at this time step point,
      //     and one-group cross section is set as those of the previous time step point.
      //     This will be corrected(modified) later.

      // Copying (using) one-group cross section at the preceding step
      // where the original problem is solved.
      for(int i=0;i<mednum;i++){
	//for(int j=0;j<nucn;j++){
	  /*
		xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
		xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
		xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	  */
	  if(i<mednum_fuel){
	    for(int j=0;j<nucn;j++){
	      if(nuclide_info[j]!=0){
		if(nuclide_info[j]==1){
		  xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
		};
		xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
		xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	      };
	    };
	  };

	  //};
      };
    }; // the end of [if((st<50)||((st-50)%pp==0)]
    // ------------------------------------------------------------------------------    


    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    //     Neutron transport calculations for 3x3 multicell including Gd-bearing rod at a center 
    //     This part is conducted at the every burnup step.
    for(int ii=0;ii<num_33_system;ii++){
    
    MECSystem lat2(group,11); // "11" is the number of mediums and hard-coded
    lat2.PutTrajectorySet(&sys_f2);  

    // ... fuel medium data preparation
    for(int i=0;i<9;i++){
      int medid=gd_id_list[ii]+i;
      if(i==8)medid=uo2_id_list[ii];
      lat2.AddMedium(lat.GetMedium(medid));
    };

    // ... non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat2.AddMedium(med[1+jj]);
    };

    // ... eigenvalue calculation
    lat2.PutRegMed(region_medium_33);
    lat2.PutGeneralOption(opt);
    lat2.PutThermalIteration(3);
    lat2.PutPL(0);
    lat2.NoCMRAcceleration();
    lat2.Print();
    real keff_dummy=lat2.CalIgen();

    vector<GroupData1D> flx_med_33(8); // "8" is the number of mediums in Gd-pin
    vector<real> gd5_xsc1g(8);
    vector<real> gd7_xsc1g(8);
    int pos_gd5, pos_gd7;
    for(int i=0;i<8;i++){
      GroupData1D flx=lat2.GetIntegratedFlux(i);
      flx_med_33[i]=flx*(1./vol_med[gd_id_list[ii]+i]); // medium-wise neutron flux data is stored here      
      // +++ One-group cross section calculation
      for(int j=0;j<nucn;j++){
        if(med[0].GetNuclideInTurn(j).GetMatnum()==641550){
	  gd5_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd5=j;
	}else if(med[0].GetNuclideInTurn(j).GetMatnum()==641570){
	  gd7_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd7=j;
	};
      };
    };

    
    //sasuga addition
    if((st<50)||((st-50)%pp==0)){

      // Since neutron transport equation is solved at this time point,
      // difference in one-group cross section between the original system and the simplified system is calculated and stored.
      
      //if(st==0){
      st_10=st; 
      cout<<"ok"<<endl;
      factor_3_3=0; // Total neutron flux of the 3x3 multi-cell in the simplified system
      factor_as=0;  // Total neutron flux of the 3x3 multi-cell in the ORIGINAL system
      for(int i=0;i<8;i++){
	factor_3_3+=flx_med_33[i].get_sum()*vol_med[gd_id_list[ii]+i];
	factor_as+=flx_med[gd_id_list[ii]+i].get_sum()*vol_med[gd_id_list[ii]+i];
	cout<<"# flux"<<endl;
	cout<<flx_med_33[i].get_sum()<<" "<<flx_med[gd_id_list[ii]+i].get_sum()<<endl;
      };

      for(int i=0;i<8;i++){
	// ABSOLUTE difference in one-group cross section is calculated and stored here.
	dxsc_5[ii][i]=gd5_xsc1g[i]-xsc_1g[st][gd_id_list[ii]+i][pos_gd5];
	dxsc_7[ii][i]=gd7_xsc1g[i]-xsc_1g[st][gd_id_list[ii]+i][pos_gd7];

	//flux_med_unity[ii][i]=flx_med[gd_id_list[ii]+i]*factor_3_3/factor_as;//集合体ベースのfluxを3×3体系で規格化
	//dflux[ii][i]=flx_med_33[i]-flux_med_unity[ii][i];
	//flux_10[ii][i]=flux_med_unity[ii][i];

	// Difference in medium-wise neutron flux is stored here
	if(ii==0) dflux_1[i]=flx_med_33[i]-flx_med[gd_id_list[ii]+i]*factor_3_3/factor_as;
	if(ii==1) dflux_2[i]=flx_med_33[i]-flx_med[gd_id_list[ii]+i]*factor_3_3/factor_as;

	fout1<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" "<<gd5_xsc1g[i]<<" ";
	fout2<<flx_med_33[i].get_sum()<<" "<<flx_med_33[i].get_sum()*factor_as/factor_3_3<<" "<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
      };
      fout1<<endl;
      fout2<<endl;
    }else{

      // Since neutron transport equation is NOT solved at this time point,
      // one-group cross section is corrected.
      
      fout1<<"aaa"<<endl;
      fout2<<"aaa"<<endl;
      for(int i=0;i<8;i++){
        cout<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
        cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
        cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
        cout<<endl;

	// One-group cross section in the 3x3 multicell system
	// is corrected by the difference between the original and simplified system.
	// Note that the ABSOLUTE difference is corrected.
	gd5_xsc1g[i]=gd5_xsc1g[i]-dxsc_5[ii][i];
	gd7_xsc1g[i]=gd7_xsc1g[i]-dxsc_7[ii][i];
	// Neutron flux in the 3x3 multicell system
	// is corrected by the difference between the original and simplified system.
	// Note that the ABSOLUTE difference is corrected.
	if(ii==0) flx_med_33[i]=flx_med_33[i]-dflux_1[i];
	if(ii==1) flx_med_33[i]=flx_med_33[i]-dflux_2[i];
	/*
	//補正した値を燃焼計算に利用------
	flx_med[gd_id_list[ii]+i]=flx_med_33[i]*factor_as/factor_3_3;
	xsc_1g[st][gd_id_list[ii]+i][pos_gd5]=gd5_xsc1g[i];
	xsc_1g[st][gd_id_list[ii]+i][pos_gd7]=gd7_xsc1g[i];
	//--------------------------------------------
	*/
	fout1<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" "<<gd5_xsc1g[i]<<" ";
	fout2<<flx_med_33[i].get_sum()<<" "<<flx_med_33[i].get_sum()*factor_as/factor_3_3<<" "<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";

      };
      fout1<<endl;
      fout2<<endl;
	  
    };   // the end of [if((st<50)||((st-50)%pp==0)]
    
    	
    cout<<"# Total neutron flux (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
      cout<<flx_med_33[i].get_sum()<<" ";
      cout<<"\n";
    };
    cout<<"\n";

    cout<<"# Gd-155 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
      cout<<gd5_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    cout<<"# Gd-157 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
      cout<<gd7_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    };
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  }; // loop-end of burnup step
  
};


void MulticellBurner::ForwardCalculationNEL2020SasugaFinal(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  // Pre-determined parameters/constants
  ofstream fout1("sigma5.dat");
  ofstream fout2("sigma7.dat");
  ofstream fout3("out_flux.dat");
  int num_33_system=2;
  int gd_id_list[]={16,40};
  int uo2_id_list[]={24,48};
  vector<int> region_medium_33;
  int region_medium_num=9+12+8*4*2;
  int region_medium_33_inp[]={
    10,10,10,10,10, 10,10,10,10, // background meshes(9)
    10,10,10,9,0,1,2,3,4,5,6,7,  // center pin(12)
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // right pin
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // top pin    
  };
  for(int i=0;i<region_medium_num;i++){
    region_medium_33.push_back(region_medium_33_inp[i]);
  };
  
  
  
  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  vector<GroupData1D> dflux_1;
  vector<GroupData1D> dflux_2;
  flx_med.resize(mednum);
  dflux_1.resize(mednum);
  dflux_2.resize(mednum);
  
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };
    
    
  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);
  //sasuga addition
  double xsc5_pre[num_33_system][8];//Gd-155の断面積の差を保存しておく
  double xsc7_pre[num_33_system][8];//Gd-155の断面積の差を保存しておく
  int st_10;
  double factor_3_3;
  double factor_as;

  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";


    // +++ SOLVING NEUTRON TRANSPORT EQUATION +++++++++++++++++++++++++++++++++++++++

    // +++ Instance generation of MECSystem to solve neutron transport equation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);  

    // +++ fuel medium data preparation
    for(int i=0;i<mednum_fuel;i++){

      // (Self-shielding calculation)
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
      };
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library      
      med[0].CalMacroFromMicro();

      // (Macro&Micro data storing)
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      lat.AddMedium(med[0]); 
      lat.GetMedium(i).NuclideClear();

    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.Print();
    int pp=100; //sasuga addition
    if(st%pp==0){ //sasuga addition
    //if(st==0){ //sasuga addition
      keff[st]=lat.CalIgen(); 
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(int i=0;i<mednum;i++){
	GroupData1D flx=lat.GetIntegratedFlux(i);
	flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
	// +++ One-group cross section storing
	if(i<mednum_fuel){
	  for(int j=0;j<nucn;j++){
	    if(nuclide_info[j]!=0){
	      if(nuclide_info[j]==1){
		xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	      };
	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
	    };
	  };
	};
      };
    }else{
      
      for(int i=0;i<mednum;i++){
	  if(i<mednum_fuel){
	    for(int j=0;j<nucn;j++){
	      if(nuclide_info[j]!=0){
		if(nuclide_info[j]==1){
		  xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
		};
		xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
		xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	      };
	    };
	  };

	  };      
      };
      
    // +++ Neutron transport calculations for 3x3 multicell including Gd-bearing rod at a center
    for(int ii=0;ii<num_33_system;ii++){
    
    MECSystem lat2(group,11); // "11" is the number of mediums and hard-coded
    lat2.PutTrajectorySet(&sys_f2);  

    // +++ fuel medium data preparation
    for(int i=0;i<9;i++){
      int medid=gd_id_list[ii]+i;
      if(i==8)medid=uo2_id_list[ii];
      lat2.AddMedium(lat.GetMedium(medid));
    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat2.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat2.PutRegMed(region_medium_33);
    lat2.PutGeneralOption(opt);
    lat2.PutThermalIteration(3);
    lat2.PutPL(0);
    lat2.NoCMRAcceleration();
    lat2.Print();
    real keff_dummy=lat2.CalIgen();

    vector<GroupData1D> flx_med_33(8);
    vector<real> gd5_xsc1g(8);
    vector<real> gd7_xsc1g(8);
    int pos_gd5, pos_gd7;
    for(int i=0;i<8;i++){
      GroupData1D flx=lat2.GetIntegratedFlux(i);
      flx_med_33[i]=flx*(1./vol_med[gd_id_list[ii]+i]); // medium-wise neutron flux data is stored here      
      // +++ One-group cross section calculation
      for(int j=0;j<nucn;j++){
        if(med[0].GetNuclideInTurn(j).GetMatnum()==641550){
	  gd5_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd5=j;
	}else if(med[0].GetNuclideInTurn(j).GetMatnum()==641570){
	  gd7_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd7=j;
	};
      };
    };
    
    //sasuga addition
    if(st%pp==0){
      //if(st==0){
      st_10=st;
    }else{
      for(int i=0;i<8;i++){
	
	//補正した値を燃焼計算に利用------
	xsc_1g[st][gd_id_list[ii]+i][pos_gd5]=xsc_1g[st-1][gd_id_list[ii]+i][pos_gd5]*gd5_xsc1g[i]/xsc5_pre[ii][i];
	xsc_1g[st][gd_id_list[ii]+i][pos_gd7]=xsc_1g[st-1][gd_id_list[ii]+i][pos_gd7]*gd7_xsc1g[i]/xsc7_pre[ii][i];
	//--------------------------------------------
	
      };
	  
      //----------------------------------
    };
    
      cout<<"# Total neutron flux (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
      cout<<flx_med_33[i].get_sum()<<" ";
      cout<<"\n";
    };
    cout<<"\n";

    cout<<"# Gd-155 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
      cout<<gd5_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    cout<<"# Gd-157 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
      cout<<gd7_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    for(int i=0;i<8;i++){
      fout1<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
      fout1<<gd5_xsc1g[i]<<" ";
      
      fout2<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
      fout2<<gd7_xsc1g[i]<<" ";
      
      fout3<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
      fout3<<flx_med_33[i].get_sum()<<" ";
    };
    fout1<<"\n";
    fout2<<"\n";
    fout3<<"\n";
    
    for(int i=0;i<8;i++){
      xsc5_pre[ii][i]=gd5_xsc1g[i];
      xsc7_pre[ii][i]=gd7_xsc1g[i];
    };
    
    };
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  }; // loop-end of burnup step
  
};

void MulticellBurner::ForwardCalculationADTS_DifSys(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  // fuga addition ------
    // output data
  ofstream fout_1g("xs1g.txt");
  ofstream fout_tf("total_flux.txt");
  ofstream fout_om("micXS_om.txt");
  ofstream fout_sm("micXS_sm.txt");
  ofstream fout_fai("fuelneutronflux.txt");

#if 1
  int each_nuc_turn[]={641550,641570};
#endif

#if 0
  int prt_nuc=2;
  string prt_nuc_nam[]={"Gd155","Gd157"}; // "HM" or "FP" or "ALL" or each nuclide's name
#endif
  // --------------------


  // Pre-determined parameters/constants
  int num_33_system=2;
  int gd_id_list[]={16,40};
  int uo2_id_list[]={24,48};
  vector<int> region_medium_33;
  int mednum_fuel_sm=9; // fuga addition
  int region_medium_num=9+12+8*4*2;
  int region_medium_33_inp[]={
    10,10,10,10,10, 10,10,10,10, // background meshes(9)
    10,10,10,9,0,1,2,3,4,5,6,7,  // center pin(12)
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // right pin
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // top pin    
  };
  for(int i=0;i<region_medium_num;i++){
    region_medium_33.push_back(region_medium_33_inp[i]);
  };
  
  
  
  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  // fuga addition
  vector< vector< vector<real> > > xsc_1g_sm; // [step][medid][nucn]
  vector< vector< vector<real> > > xsn2n_1g_sm;
  vector< vector< vector<real> > > xsf_1g_sm;
  vector< vector< vector<real> > > pre_xsc_1g_sm; // [step][medid][nucn]
  vector< vector< vector<real> > > pre_xsn2n_1g_sm;
  vector< vector< vector<real> > > pre_xsf_1g_sm;
  vector< vector< vector<GroupData1D> > > mic_sigf_sm; //[step][medid][nucn](group)
  vector< vector< vector<GroupData1D> > > mic_sigc_sm;
  vector< vector< vector<GroupData1D> > > mic_sign2n_sm; 
  // -------------

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  vector<GroupData1D> dflux_1;
  vector<GroupData1D> dflux_2;
  flx_med.resize(mednum);
  dflux_1.resize(mednum);
  dflux_2.resize(mednum);
  
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  // fuga addition
  xsf_1g_sm.resize(tmp);
  xsc_1g_sm.resize(tmp);
  xsn2n_1g_sm.resize(tmp);
  pre_xsf_1g_sm.resize(tmp);
  pre_xsc_1g_sm.resize(tmp);
  pre_xsn2n_1g_sm.resize(tmp);
  mic_sigf_sm.resize(tmp);
  mic_sigc_sm.resize(tmp);
  mic_sign2n_sm.resize(tmp);
  // -------------
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    // fuga addition
    xsf_1g_sm[i].resize(mednum_fuel_sm);
    xsc_1g_sm[i].resize(mednum_fuel_sm);
    xsn2n_1g_sm[i].resize(mednum_fuel_sm);
    pre_xsf_1g_sm[i].resize(mednum_fuel_sm);
    pre_xsc_1g_sm[i].resize(mednum_fuel_sm);
    pre_xsn2n_1g_sm[i].resize(mednum_fuel_sm);
    mic_sigf_sm[i].resize(mednum_fuel_sm);
    mic_sigc_sm[i].resize(mednum_fuel_sm);
    mic_sign2n_sm[i].resize(mednum_fuel_sm);
    //-------------
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      // fuga addition
      xsf_1g_sm[i][j].resize(nucn,0.);
      xsc_1g_sm[i][j].resize(nucn,0.);
      xsn2n_1g_sm[i][j].resize(nucn,0.);
      pre_xsf_1g_sm[i][j].resize(nucn,0.);
      pre_xsc_1g_sm[i][j].resize(nucn,0.);
      pre_xsn2n_1g_sm[i][j].resize(nucn,0.);
      mic_sigf_sm[i][j].resize(nucn);
      mic_sigc_sm[i][j].resize(nucn);
      mic_sign2n_sm[i][j].resize(nucn);
      // -------------
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
            mic_sigf_sm[i][j][k].put_imax(ngrp); // fuga addition
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
          mic_sigc_sm[i][j][k].put_imax(ngrp); // fuga addition
          mic_sign2n_sm[i][j][k].put_imax(ngrp); // fuga addition
        };
      };
    };
  };
  };
    
    
  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();

  // ... Preparing medium class instances for ADTS model
  //
  vector< vector<Medium> > med_adts(num_33_system); // fuga addition
  for(int ii=0;ii<num_33_system;ii++){
    med_adts[ii].resize(mednum_fuel_sm);
  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);
  //sasuga addition
  double xsc5_pre[num_33_system][8];//Gd-155の断面積の差を保存しておく
  double xsc7_pre[num_33_system][8];//Gd-155の断面積の差を保存しておく
  int st_10;
  double factor_3_3;
  double factor_as;

  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    if(st%sm_step==0){ // fuga addition
      st_10=st;        // fuga addition

    // +++ SOLVING NEUTRON TRANSPORT EQUATION +++++++++++++++++++++++++++++++++++++++

    // +++ Instance generation of MECSystem to solve neutron transport equation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);  

    // +++ fuel medium data preparation
    for(int i=0;i<mednum_fuel;i++){

      // (Self-shielding calculation)
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
      };
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library
    // fuga addition
    for(int ii=0;ii<num_33_system;ii++){
      for(int j=0;j<8;j++){
  if(i==gd_id_list[ii]+j) med_adts[ii][j]=med[0];
      };
  if(i==uo2_id_list[ii]) med_adts[ii][8]=med[0]; 
    }; 
    // ----------------   
    // med_adts[i]=med[0]; // ... adts
      med[0].CalMacroFromMicro();

      // (Macro&Micro data storing)
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      lat.AddMedium(med[0]); 
      lat.GetMedium(i).NuclideClear();

    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.Print();
    // int pp=100; //sasuga addition
    // if(st%==0){ //sasuga addition
    //if(st==0){ //sasuga addition
      keff[st]=lat.CalIgen(); 
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(int i=0;i<mednum;i++){
	GroupData1D flx=lat.GetIntegratedFlux(i);
	flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
	// +++ One-group cross section storing
	if(i<mednum_fuel){
	  for(int j=0;j<nucn;j++){
	    if(nuclide_info[j]!=0){
	      if(nuclide_info[j]==1){
		xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	      };
	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
	    };
	  };
	};
      };
    }else{
      
      for(int i=0;i<mednum;i++){
	  if(i<mednum_fuel){
	    for(int j=0;j<nucn;j++){
	      if(nuclide_info[j]!=0){
		if(nuclide_info[j]==1){
		  xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
		};
		xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
		xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	      };
	    };
	  };

	  };      
      };
      
    // +++ Neutron transport calculations for 3x3 multicell including Gd-bearing rod at a center
    for(int ii=0;ii<num_33_system;ii++){
    
    MECSystem lat2(group,11); // "11" is the number of mediums and hard-coded
    lat2.PutTrajectorySet(&sys_f2);  

    // +++ fuel medium data preparation
    // for(int i=0;i<9;i++){
    for(int i=0;i<mednum_fuel_sm;i++){ // fuga addition
  /*
      // (Copy Medium in Original-Model)-----------
      int medid=gd_id_list[ii]+i;
      if(i==8)medid=uo2_id_list[ii];
      lat2.AddMedium(lat.GetMedium(medid));
      // -----------------------------------------
  */
    // fuga addition---
      // Use of the Original burnup step results------
      int nucn=med_adts[ii][i].GetNucnum();
      for(int j=0;j<nucn;j++){
        med_adts[ii][i].GetNuclideInTurn(j).PutDensity(fwd_nuc[st][0][i].get_dat(j));
      };
      med_adts[ii][i].CalMacroFromMicro();
      Medium bmed=med_adts[ii][i].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
      // --------------------------------------------

      // Macro&Micro data storing
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf_sm[bstmp][i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc_sm[bstmp][i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n_sm[bstmp][i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      lat2.AddMedium(bmed);
      lat2.GetMedium(i).NuclideClear();
    // ---------------

    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      // lat2.AddMedium(med[1+jj]);

    // fuga addition
      Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel+jj],true); // [true] is for microscopic XS
      lat2.AddMedium(bmed);
    //--------------
    };

    // +++ eigenvalue calculation
    lat2.PutRegMed(region_medium_33);
    lat2.PutGeneralOption(opt);
    lat2.PutThermalIteration(3);
    lat2.PutPL(0);
    lat2.NoCMRAcceleration();
    lat2.Print();
    real keff_dummy=lat2.CalIgen();
    vector<real> medium_wise_flux_change(mednum_fuel_sm); // fuga addition

    // vector<GroupData1D> flx_med_33(8);
    vector<GroupData1D> flx_med_sm(mednum_fuel_sm); // fuga addition
    // vector<real> gd5_xsc1g(8);
    // vector<real> gd7_xsc1g(8);
    // int pos_gd5, pos_gd7;
    // for(int i=0;i<8;i++){
      // GroupData1D flx=lat2.GetIntegratedFlux(i);
      // flx_med_33[i]=flx*(1./vol_med[gd_id_list[ii]+i]); // medium-wise neutron flux data is stored here      

  // fuga addition
    for(int i=0;i<mednum_fuel_sm;i++){
      GroupData1D flx_sm=lat2.GetIntegratedFlux(i);

      // ... Total neutron flux change during the subsequesnt ADTS steps is calculated.
      medium_wise_flux_change[i]=(flx_sm.get_sum()/vol_med[i])/flx_med_sm[i].get_sum();
      flx_med_sm[i]=flx_sm*(1./vol_med[i]); // medium-wise neutron flux data is stored here 

      // +++ One-group cross section calculation
      for(int j=0;j<nucn;j++){
        if(nuclide_info[j]!=0){
          if(nuclide_info[j]==1){
        xsf_1g_sm[bstmp][i][j]=mic_sigf_sm[bstmp][i][j].Cond(flx_sm);
      };
    xsc_1g_sm[bstmp][i][j]=mic_sigf_sm[bstmp][i][j].Cond(flx_sm);
    xsn2n_1g_sm[bstmp][i][j]=mic_sign2n_sm[bstmp][i][j].Cond(flx_sm);
        };
      };
  // -----------

  //     for(int j=0;j<nucn;j++){
  //       if(med[0].GetNuclideInTurn(j).GetMatnum()==641550){
	//   gd5_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	//   pos_gd5=j;
	// }else if(med[0].GetNuclideInTurn(j).GetMatnum()==641570){
	//   gd7_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	//   pos_gd7=j;
	// };
  //     };
    };
    
    //sasuga addition
    // if(st%pp==0){
    if(st%sm_step==0){ // fuga addition
      //if(st==0){
      st_10=st;
      // fuga addition
      for(int i=0;i<mednum_fuel_sm;i++){
    for(int k=0;k<nucn;k++){
      int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
      if(nuclide_info[k]!=0){
        if(nuclide_info[k]==1){
          pre_xsf_1g_sm[bstmp][i][k]=xsf_1g_sm[bstmp][i][k];
        };
        pre_xsc_1g_sm[bstmp][i][k]=xsc_1g_sm[bstmp][i][k];
        pre_xsn2n_1g_sm[bstmp][i][k]=xsn2n_1g_sm[bstmp][i][k];
      };
    };
      };
      // --------------

    }else{
  //     for(int i=0;i<8;i++){
	
	// //補正した値を燃焼計算に利用------
	// // xsc_1g[st][gd_id_list[ii]+i][pos_gd5]=xsc_1g[st-1][gd_id_list[ii]+i][pos_gd5]*gd5_xsc1g[i]/xsc5_pre[ii][i];
	// // xsc_1g[st][gd_id_list[ii]+i][pos_gd7]=xsc_1g[st-1][gd_id_list[ii]+i][pos_gd7]*gd7_xsc1g[i]/xsc7_pre[ii][i];
	// //----------------------------------------
	  
  //   };
     
    // fuga addition
      // ... One-group cross section correction in ADTS
      for(int i=0;i<mednum_fuel_sm;i++){
    for(int k=0;k<nucn;k++){
      int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
      if(nuclide_info[k]!=0){
        if(nuclide_info[k]==1){
          xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_sm[bstmp][i][k]/pre_xsf_1g_sm[bstmp][i][k]); // for neutron fission cross section
        };
        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_sm[bstmp][i][k]/pre_xsc_1g_sm[bstmp][i][k]);
        if(xsn2n_1g[st-1][i][k]>0.){
          xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_sm[bstmp][i][k]/pre_xsn2n_1g_sm[bstmp][i][k]);
        };
      };
    };
      // ... Neutron flux correction in ADTS
      for(int i=0;i<mednum_fuel_sm;i++){
        flx_med[i]=flx_med[i]*medium_wise_flux_change[i];
      };

      };
    // --------------
    };

    
    //   cout<<"# Total neutron flux (original & 3x3)\n";
    // for(int i=0;i<8;i++){
    //   cout<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
    //   cout<<flx_med_33[i].get_sum()<<" ";
    //   cout<<"\n";
    // };
    // cout<<"\n";

    // cout<<"# Gd-155 one-group (n,g) cross section (original & 3x3)\n";
    // for(int i=0;i<8;i++){
    //   cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
    //   cout<<gd5_xsc1g[i]<<" ";
    //   cout<<"\n";
    // };
    // cout<<"\n";
    
    // cout<<"# Gd-157 one-group (n,g) cross section (original & 3x3)\n";
    // for(int i=0;i<8;i++){
    //   cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
    //   cout<<gd7_xsc1g[i]<<" ";
    //   cout<<"\n";
    // };
    // cout<<"\n";
    
    // for(int i=0;i<8;i++){
    //   fout1<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
    //   fout1<<gd5_xsc1g[i]<<" ";
      
    //   fout2<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
    //   fout2<<gd7_xsc1g[i]<<" ";
      
    //   fout3<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
    //   fout3<<flx_med_33[i].get_sum()<<" ";
    // };
    // fout1<<"\n";
    // fout2<<"\n";
    // fout3<<"\n";
    
    // for(int i=0;i<8;i++){
    //   xsc5_pre[ii][i]=gd5_xsc1g[i];
    //   xsc7_pre[ii][i]=gd7_xsc1g[i];
    // };
    
    }; // end of number of Gd-pin
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){

  //     for(int i=0;i<totm;i++){
  //       if(region_medium[i]<mednum_fuel){
  //         real vol=lat.GetMesh(i).GetVolume();
  //         volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	// };
  //     };

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  }; // loop-end of burnup step
  
};


void MulticellBurner::ForwardCalculationNEL2020_3x3(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  // fuga addition ------------------------
    // output data file
    ofstream fout_1g("xs1g_5c7c.txt");
    ofstream fout_tf("total_flux.txt");
    ofstream fout_sm("fuelFluxXs_sm.txt");
    ofstream fout_om("fuelfluxXs_om.txt");

    fout_1g<<"# One-group capture cross section in Gd-155,-157 [barn]\n";
    fout_1g<<"# burnup [GWD/t] : 1 \n";
    for(int i=0;i<mednum_fuel;i++){
      fout_1g<<"# Med "<<IntToString(i)<<" : ";
      for(int j=0;j<4;j++){
        string column=IntToString(4*i+j+2);
        fout_1g<<column<<" ";
        };
      fout_1g<<endl;
    };
    fout_1g<<"# (Burnup)   ";
    for(int i=0;i<mednum_fuel;i++){
      fout_1g<<" (Med "+IntToString(i)+")";
      for(int j=0;j<13*4-8;j++){fout_1g<<" ";};
    };
    fout_1g<<endl;
    fout_1g<<"# [GWd/t]    ";
    for(int i=0;i<mednum_fuel;i++){fout_1g<<" (1g5c.ori)   (1g5c.sim)   (1g7c.ori)   (1g7c.sim)  ";};
    fout_1g<<endl;

    fout_tf<<"# Total neutron flux [/cm2/s] \n";
    for(int i=0;i<mednum_fuel;i++){
      fout_tf<<"# Med "<<IntToString(i)<<" : ";
      for(int j=0;j<2;j++){
        string column=IntToString(2*i+j+2);
        fout_tf<<column<<" ";
      };
      fout_tf<<endl;
    };
    fout_tf<<"# (Burnup)   ";
    for(int i=0;i<mednum_fuel;i++){
      fout_tf<<" (Med "+IntToString(i)+")";
      for(int j=0;j<13*2-8;j++){fout_tf<<" ";};
    };
    fout_tf<<endl;
    fout_tf<<"# [GWd/t]    ";
    for(int i=0;i<mednum_fuel;i++){fout_tf<<" (tFlux.ori)  (tFlux.sim) ";};
    fout_tf<<endl;
    //--------------------------------------


  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };
  
  vector< vector< vector<real> > > xsc_1g_2; // [step][medid][nucn]
  vector< vector< vector<real> > > xsn2n_1g_2;
  vector< vector< vector<real> > > xsf_1g_2;
  vector< vector< vector<GroupData1D> > > mic_sigf_2; //[step][medid][nucn](group)
  vector< vector< vector<GroupData1D> > > mic_sigc_2;
  vector< vector< vector<GroupData1D> > > mic_sign2n_2; 
  
  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  //sasuga addition
  xsc_1g_2.resize(burn_step+1);
  xsn2n_1g_2.resize(burn_step+1);
  xsf_1g_2.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    //sasuga adition
    xsc_1g_2[i].resize(mednum_fuel);
    xsn2n_1g_2[i].resize(mednum_fuel);
    xsf_1g_2[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
       //sasuga addition
      xsc_1g_2[i][j].resize(nucn,0.);
      xsn2n_1g_2[i][j].resize(nucn,0.);
      xsf_1g_2[i][j].resize(nucn,0.);
   };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  flx_med.resize(mednum);
  vector<GroupData1D> flx_med_2;
  flx_med_2.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
    flx_med_2[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  //sasuga addition
  mic_sigf_2.resize(tmp); 
  mic_sigc_2.resize(tmp); 
  mic_sign2n_2.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    //sasuga addiiton
    mic_sigf_2[i].resize(mednum_fuel);
    mic_sigc_2[i].resize(mednum_fuel);
    mic_sign2n_2[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      //sasuga addition
      mic_sigf_2[i][j].resize(nucn);
      mic_sigc_2[i][j].resize(nucn);
      mic_sign2n_2[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
            // mic_sigf_2[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
          // mic_sigc_2[i][j][k].put_imax(group);
          // mic_sign2n_2[i][j][k].put_imax(group);
        };
      };
    };
  };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();


  // ... Preparing medium class instances for ADTS model
  //
  //
  vector<Medium> med_adts(mednum_fuel);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  double dxsf[mednum_fuel][nucn];
  double dxsc[mednum_fuel][nucn];
  vector<GroupData1D> pre_xsf_1g_smlux2;
  pre_xsf_1g_smlux2.resize(mednum);
  double dxsn2n[mednum_fuel][nucn];
  int st_10;
  for(int st=0;st<burn_step+1;st++){
    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    
    if(st==0||st%sm_step==0){

      st_10=st;      

      // +++ Instance generation of MECSystem to solve neutron transport equation
      MECSystem lat(group,mednum);
      lat.PutTrajectorySet(&sys_f);  

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel;i++){

        // (Self-shielding calculation)
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library
	med_adts[i]=med[0];  // ... adts
        med[0].CalMacroFromMicro();

        // (Macro&Micro data storing)
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();

      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };

      // +++ eigenvalue calculation
      lat.PutRegMed(region_medium);
      lat.PutGeneralOption(opt);
      lat.PutThermalIteration(3);
      lat.PutPL(0);
      lat.NoCMRAcceleration();
      lat.Print();
      
      keff[st]=lat.CalIgen();

      for(int i=0;i<mednum;i++){
        GroupData1D flx=lat.GetIntegratedFlux(i);
        flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	      };
	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
            };
	  };
        };
      };

    }else{

      // One-group cross section at the last burnup step is assumed
      for(int i=0;i<mednum;i++){
        if(i<mednum_fuel){
	  for(int j=0;j<nucn;j++){
	    if(nuclide_info[j]!=0){
	      if(nuclide_info[j]==1){
	        xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
	      };
	      xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
	      xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	    };
	  };
	};
      };
      
    };
      
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){


    #if 0
      // ... No group collapsing ...
      // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
      MECSystem lat2(group,mednum);
      lat2.PutTrajectorySet(&sys_f2);        

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel;i++){

	/*
        // (Self-shielding calculation) -----------------------------------------------------------
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
        med[0].CalMacroFromMicro();
	// ----------------------------------------------------------------
	*/

	// (Use of the original burnup step results) -------------------------------
        int nucn=med_adts[i].GetNucnum();
        for(int ii=0;ii<nucn;ii++){
           med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
        };
	med_adts[i].CalMacroFromMicro();
	med[0]=med_adts[i];
	// -------------------------------------------------------------------------
	
	
        // (Macro&Micro data storing)
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat2.AddMedium(med[0]); 
        lat2.GetMedium(i).NuclideClear();

      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat2.AddMedium(med[1+jj]);
      };

      // +++ eigenvalue calculation
      lat2.PutRegMed(region_medium);
      lat2.PutGeneralOption(opt);
      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      lat2.Print();

      real keff_dummy=lat2.CalIgen();

      vector<real> medium_wise_flux_change(mednum);
      
      for(int i=0;i<mednum;i++){
        GroupData1D flx_2=lat2.GetIntegratedFlux(i);

	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
	medium_wise_flux_change[i]=(flx_2.get_sum()/vol_med[i])/flx_med_2[i].get_sum();

        flx_med_2[i]=flx_2*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g_2[st][i][j]=mic_sigf_2[bstmp][i][j].Cond(flx_2);
  	      };
	      xsc_1g_2[st][i][j]=mic_sigc_2[bstmp][i][j].Cond(flx_2);
	      xsn2n_1g_2[st][i][j]=mic_sign2n_2[bstmp][i][j].Cond(flx_2);
            };
  	  };
        };
      };
    #endif      


    #if 1
      // ... with group collapsing from 107-group to fewer-group ...
      // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
      

      //int ngrp=21; // The number of energy groups after the group collapsing
      //int bgrp[]={4,9,14,19,24, 29,34,39,44,49, 54,59,64,69,74, 79,84,89,94,99, 106};
      // Energy group boundary in the collapsed structure.
      // The final one should be [106] when the number of groups in the original structure is the 107-group.
      
      // int ngrp=11;
      // int bgrp[]={9,19,29,39,49, 59,69,79,89,99, 106};
      /*
      int ngrp=107;
      int *bgrp=new int[ngrp];
      for(int i=0;i<ngrp;i++){
	bgrp[i]=i;
      };
      */

      vector< vector<GroupData1D> > mic_sigf_3(mednum_fuel);  // [medid][nucn](group)
      vector< vector<GroupData1D> > mic_sigc_3(mednum_fuel);
      vector< vector<GroupData1D> > mic_sign2n_3(mednum_fuel); 
	
      MECSystem lat2(ngrp,mednum);
      lat2.PutTrajectorySet(&sys_f2);        

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel;i++){

	/*
        // (Self-shielding calculation) ---------------------------------------
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
        med[0].CalMacroFromMicro();

        Medium bmed=med[0].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
        // ------------------------------------------------------------------------
	*/

	// (Use of the original burnup step results) -------------------------------
        int nucn=med_adts[i].GetNucnum();
        for(int ii=0;ii<nucn;ii++){
           med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
        };
	med_adts[i].CalMacroFromMicro();
	Medium bmed=med_adts[i].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
	// -------------------------------------------------------------------------


        // (Macro&Micro data storing)
	mic_sigf_3[i].resize(nucn);
	mic_sigc_3[i].resize(nucn);
	mic_sign2n_3[i].resize(nucn);	
        for(int k=0;k<nucn;k++){
	  
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat2.AddMedium(bmed); 
        lat2.GetMedium(i).NuclideClear();

      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel*jj],true); // [true] is for microscopic XS
        lat2.AddMedium(bmed);
      };

      // +++ eigenvalue calculation
      lat2.PutRegMed(region_medium);
      lat2.PutGeneralOption(opt);
      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      lat2.Print();

      real keff_dummy=lat2.CalIgen();

      vector<real> medium_wise_flux_change(mednum);
      
      for(int i=0;i<mednum;i++){
        GroupData1D flx_2=lat2.GetIntegratedFlux(i);

	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
	medium_wise_flux_change[i]=(flx_2.get_sum()/vol_med[i])/flx_med_2[i].get_sum();

        flx_med_2[i]=flx_2*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g_2[st][i][j]=mic_sigf_3[i][j].Cond(flx_2);
  	      };
	      xsc_1g_2[st][i][j]=mic_sigc_3[i][j].Cond(flx_2);
	      xsn2n_1g_2[st][i][j]=mic_sign2n_3[i][j].Cond(flx_2);
            };
  	  };
        };
      };
    #endif      



      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      if(st==0||st%sm_step==0){	
        //st_10=st;
      }else{

        // ... One-group cross section correction in ADTS	
	
        for(int i=0;i<mednum_fuel;i++){
	  for(int k=0;k<nucn;k++){
	  
  	    int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
    #if 1	  
            // ... Copy of the ADTS step data
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
  	        xsf_1g[st][i][k]=xsf_1g_2[st][i][k];
	      };
	      xsc_1g[st][i][k]=xsc_1g_2[st][i][k];
	      if(xsn2n_1g[st-1][i][k]>0.){
  	        xsn2n_1g[st][i][k]=xsn2n_1g_2[st][i][k];
	      };
	    };
    #endif	  

    #if 0	  
            // ... Correcting all nuclides in all burnup regions
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
  	        xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
	      };
	      xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
	      if(xsn2n_1g[st-1][i][k]>0.){
  	        xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
	      };
	    };
    #endif	  

    #if 0	  
	    if(i<8){
              // ... Correcting all nuclides in the central Gd-bearing pin
	      if(nuclide_info[k]!=0){
	        if(nuclide_info[k]==1){
  	          xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
	        };
	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
  	        if(xsn2n_1g[st-1][i][k]>0.){	      
	          xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
	        };
	      };
	    };
    #endif	  

    #if 0	  
	    if(i<8){
	      // ... Correcting Gd-155 and -157 capture cross sections in the central Gd-bearing pin
	      if(nucid==641550||nucid==641570){
	        if(nucid==641550) cout<<"Gd-155 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        if(nucid==641570) cout<<"Gd-157 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        //xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
	        //xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
	      };
	    };
    #endif	  

	  };
        };

	// ... Neutron flux correction in ADTS

	// for(int i=0;i<mednum_fuel;i++){ // All the medium included in the system
	//   //for(int i=0;i<8;i++){ // Only for Gd-bearing rod
	//   flx_med[i]=flx_med[i]*medium_wise_flux_change[i];
	// };
      };

    // fuga addition -------------------------------------
      //output data

      // Original-Model
      fout_om<<"#\n# Neutron flux per lethargy\n#\n";
      fout_om<<"#  burnup step : "<<st<<endl;
      fout_om<<"#  burnup[GWD/t] : "<<acburn[st]<<endl;
      fout_om<<"#\n";
      fout_om<<"# Energy    Med ID\n";
      fout_om<<"# [eV]     ";
      for(int i=0;i<mednum_fuel;i++){
        string strid=IntToString(i);
        fout_om<<" ("+strid+")";
        for(int j=0;j<(10*3-3-strid.size());j++){fout_om<<" ";};
      };
      fout_om<<endl;
      fout_om<<"#         ";
      for(int i=0;i<mednum_fuel;i++){
        fout_om<<" (flux)   (sigma5c) (sigma7c) ";
      };
      fout_om<<endl;

      fout_om.setf(ios::scientific);
      fout_om.precision(3);
      for(int i=0;i<group;i++){
        real e0=med[0].GetEnband().get_dat(i);
        real e1=med[0].GetEnband().get_dat(i+1);
        real letwid=log(e0/e1);
        fout_om<<e0<<" ";

        for(int j=0;j<mednum_fuel;j++){
          fout_om<<flx_med[j].get_dat(i)/letwid<<" ";
          for(int k=0;k<nucn;k++){
            int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
      if(nucid==641550) fout_om<<mic_sigc[bstmp][j][k].get_dat(i)<<" ";
      if(nucid==641570) fout_om<<mic_sigc[bstmp][j][k].get_dat(i)<<" ";    
          };
        };
        fout_om<<endl;
      };
      fout_om<<"\n\n";
      
      // Simplified-Model  
      fout_sm<<"#\n# Neutron flux per lethargy\n#\n";
      fout_sm<<"#  burnup step : "<<st<<endl;
      fout_sm<<"#  burnup[GWD/t] : "<<acburn[st]<<endl;
      fout_sm<<"#\n";
      fout_sm<<"# Energy\n";
      fout_sm<<"# [eV]    ";
      for(int i=0;i<mednum_fuel;i++){
        string strid=IntToString(i);
        fout_sm<<" ("+strid+")";
        for(int j=0;j<(10*3-3-strid.size());j++){fout_sm<<" ";};
      };
      fout_sm<<endl;
      fout_sm<<"#         ";
      for(int i=0;i<mednum_fuel;i++){
        fout_sm<<" (flux)   (sigma5c) (sigma7c) ";
      };
      fout_sm<<endl;

      fout_sm.setf(ios::scientific);
      fout_sm.precision(3);
      for(int i=0;i<ngrp;i++){
        Medium bmed=med_adts[0].Cond(ngrp,bgrp,flx_med[0],flx_med[0],true);
        real e0=bmed.GetEnband().get_dat(i);
        real e1=bmed.GetEnband().get_dat(i+1);
        real letwid=log(e0/e1);
        fout_sm<<e0<<" ";

        for(int j=0;j<mednum_fuel;j++){
          fout_sm<<flx_med_2[j].get_dat(i)/letwid<<" ";
          Medium bmed=med_adts[j].Cond(ngrp,bgrp,flx_med[j],flx_med[j],true);
          for(int k=0;k<nucn;k++){
            int nucid=bmed.GetNuclideInTurn(k).GetMatnum();
      if(nucid==641550) fout_sm<<mic_sigc_3[j][k].get_dat(i)<<" ";
      if(nucid==641570) fout_sm<<mic_sigc_3[j][k].get_dat(i)<<" ";    
          };
        };
        fout_sm<<endl;
      };
      fout_sm<<"\n\n";

      fout_1g.setf(ios::scientific);
      fout_1g.precision(6);
      fout_1g<<acburn[st]<<" ";
      for(int i=0;i<mednum_fuel;i++){
        for(int j=0;j<nucn;j++){
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    if(nucid==641550||nucid==641570) fout_1g<<xsc_1g[st][i][j]<<" "<<xsc_1g_2[st][i][j]<<" ";
        };
      };
      fout_1g<<endl;
      
      fout_tf.setf(ios::scientific); 
      fout_tf.precision(6); 
      fout_tf<<acburn[st]<<" ";
      for(int i=0;i<mednum_fuel;i++){
    fout_tf<<flx_med[i].get_sum()<<" "<<flx_med_2[i].get_sum()<<" ";
      };
      fout_tf<<endl;

    //---------------------------------------------------------



      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    /*      
      // Commented out in this method since it is unnecessary.
      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };
    */

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
    
      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          //real tmp=flx_med_2[i].get_sum();	  
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
  	  //sumflx=flx_med_2[med_target].get_sum();	  
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
	
	for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
	  //total_flux[st][j][i]=flx_med_2[i].get_sum()*power_factor[st][j];	  
	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	};
	
      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  }; // loop-end of burnup step
};

// void MulticellBurner::ForwardCalculationNEL2020_3x3_pc(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
// {
//   if(input_flux_level&&med_target==-1){
//     cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
//     cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
//     exit(0);
//   };
  
//   vector< vector< vector<real> > > xsc_1g_2; // [step][medid][nucn]
//   vector< vector< vector<real> > > xsn2n_1g_2;
//   vector< vector< vector<real> > > xsf_1g_2;
//   vector< vector< vector<GroupData1D> > > mic_sigf_2; //[step][medid][nucn](group)
//   vector< vector< vector<GroupData1D> > > mic_sigc_2;
//   vector< vector< vector<GroupData1D> > > mic_sign2n_2; 
  
//   GeneralOption opt;

//   // +++ Pre-calculation of Dancoff factor in pincell model +++
//   SelfShieldingCalculator ssc;  
//   ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
//   GroupData1D dancoff=ssc.GetDancoff(0);

//   // +++ Initialization of Bell factor +++
//   GroupData1D bell(group);
//   for(int i=0;i<group;i++){bell.put_data(i,1.2);};

//   // +++ Array setting +++
//   fwd_nuc.resize(burn_step+1);
//   xsc_1g.resize(burn_step+1);
//   xsn2n_1g.resize(burn_step+1);
//   xsf_1g.resize(burn_step+1);
//   //sasuga addition
//   xsc_1g_2.resize(burn_step+1);
//   xsn2n_1g_2.resize(burn_step+1);
//   xsf_1g_2.resize(burn_step+1);
//   total_flux.resize(burn_step);
//   delt.resize(burn_step);
//   power_factor.resize(burn_step);

//   for(int i=0;i<burn_step+1;i++){
//     int sub_step=sub_step_list[i];
//     fwd_nuc[i].resize(sub_step+1);
//     for(int j=0;j<sub_step+1;j++){
//       fwd_nuc[i][j].resize(mednum_fuel);
//       for(int k=0;k<mednum_fuel;k++){
// 	fwd_nuc[i][j][k].put_imax(nucn);
//       };
//     };
//     xsc_1g[i].resize(mednum_fuel);
//     xsn2n_1g[i].resize(mednum_fuel);
//     xsf_1g[i].resize(mednum_fuel);
//     //sasuga adition
//     xsc_1g_2[i].resize(mednum_fuel);
//     xsn2n_1g_2[i].resize(mednum_fuel);
//     xsf_1g_2[i].resize(mednum_fuel);
//     for(int j=0;j<mednum_fuel;j++){
//       xsc_1g[i][j].resize(nucn,0.);
//       xsn2n_1g[i][j].resize(nucn,0.);
//       xsf_1g[i][j].resize(nucn,0.);
//        //sasuga addition
//       xsc_1g_2[i][j].resize(nucn,0.);
//       xsn2n_1g_2[i][j].resize(nucn,0.);
//       xsf_1g_2[i][j].resize(nucn,0.);
//    };
//   };

//   for(int i=0;i<burn_step;i++){
//     int sub_step=sub_step_list[i];
//     delt[i].resize(sub_step);
//     total_flux[i].resize(sub_step);
//     power_factor[i].resize(sub_step+1);
//     for(int j=0;j<sub_step;j++){
//       total_flux[i][j].resize(mednum_fuel);
//     };
//   };

//   vector<GroupData1D> flx_med;
//   flx_med.resize(mednum);
//   vector<GroupData1D> flx_med_2;
//   flx_med_2.resize(mednum);
//   for(int j=0;j<mednum;j++){
//     flx_med[j].put_imax(group);
//     flx_med_2[j].put_imax(group);
//   };

//   volflx_mesh.resize(burn_step); 
//   for(int i=0;i<burn_step;i++){
//     volflx_mesh[i].resize(totm);
//     for(int j=0;j<totm;j++){
//       if(region_medium[j]<mednum_fuel){
//         volflx_mesh[i][j].put_imax(group);
//       };
//     };
//   };

//   {
//   int tmp=1;
//   mic_sigf.resize(tmp); 
//   mic_sigc.resize(tmp); 
//   mic_sign2n.resize(tmp);
//   //sasuga addition
//   mic_sigf_2.resize(tmp); 
//   mic_sigc_2.resize(tmp); 
//   mic_sign2n_2.resize(tmp);
//   for(int i=0;i<tmp;i++){
//     mic_sigf[i].resize(mednum_fuel);
//     mic_sigc[i].resize(mednum_fuel);
//     mic_sign2n[i].resize(mednum_fuel);
//     //sasuga addiiton
//     mic_sigf_2[i].resize(mednum_fuel);
//     mic_sigc_2[i].resize(mednum_fuel);
//     mic_sign2n_2[i].resize(mednum_fuel);
//     for(int j=0;j<mednum_fuel;j++){
//       mic_sigf[i][j].resize(nucn);
//       mic_sigc[i][j].resize(nucn);
//       mic_sign2n[i][j].resize(nucn);
//       //sasuga addition
//       mic_sigf_2[i][j].resize(nucn);
//       mic_sigc_2[i][j].resize(nucn);
//       mic_sign2n_2[i][j].resize(nucn);
//       for(int k=0;k<nucn;k++){
//         if(nuclide_info[k]!=0){
//           if(nuclide_info[k]==1){
//             mic_sigf[i][j][k].put_imax(group);
//             // mic_sigf_2[i][j][k].put_imax(group);
//           };
//           mic_sigc[i][j][k].put_imax(group);
//           mic_sign2n[i][j][k].put_imax(group);
//           // mic_sigc_2[i][j][k].put_imax(group);
//           // mic_sign2n_2[i][j][k].put_imax(group);
//         };
//       };
//     };
//   };
//   };

//   vector< vector<GroupData1D> > mic_sigf_c;
//   vector< vector<GroupData1D> > mic_sigc_c;
//   vector< vector<GroupData1D> > mic_sign2n_c;
//   vector<GroupData1D> flx_med_c;
//   if(corrector_calc){
//     //int tmp=1;
//     //if(adjoint)tmp=burn_step+1;
//     int tmp=burn_step+1;
//     xsc_1g_p.resize(tmp);
//     xsn2n_1g_p.resize(tmp);
//     xsf_1g_p.resize(tmp);
//     for(int j=0;j<tmp;j++){
//       xsc_1g_p[j].resize(mednum_fuel);
//       xsn2n_1g_p[j].resize(mednum_fuel);
//       xsf_1g_p[j].resize(mednum_fuel);
//       for(int i=0;i<mednum_fuel;i++){
//         xsc_1g_p[j][i].resize(nucn,0.);
//         xsn2n_1g_p[j][i].resize(nucn,0.);
//         xsf_1g_p[j][i].resize(nucn,0.);
//       };
//     };

//     mic_sigf_c.resize(mednum_fuel); 
//     mic_sigc_c.resize(mednum_fuel); 
//     mic_sign2n_c.resize(mednum_fuel);
//     for(int i=0;i<mednum_fuel;i++){
//       mic_sigf_c[i].resize(nucn);
//       mic_sigc_c[i].resize(nucn);
//       mic_sign2n_c[i].resize(nucn);
//     };

//     flx_med_c.resize(mednum);
//     for(int i=0;i<mednum;i++){
//       flx_med_c[i].put_imax(group);
//     };

//     total_flux_p.resize(burn_step);
//     for(int i=0;i<burn_step;i++){
//       int sub_step=sub_step_list[i];
//       total_flux_p[i].resize(sub_step);
//       for(int j=0;j<sub_step;j++){
//         total_flux_p[i][j].resize(mednum_fuel);
//       };
//     };

//   };

  
//   // ... OWPC
//   vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC (Rp)
//   vector<real> np_gd5, np_gd7; // Np

//   vector<real>  old_n0155;
//   vector<real>  old_n0157;
//   vector<real>  old_r0155;
//   vector<real>  old_r0157;
//   vector<real>  old_np155;
//   vector<real>  old_rp155;
//   vector<real>  old_rp157;
//   vector<real>  old_rc155;
//   vector<real>  old_rc157;
//   vector<real>  old_n0155_2;
//   vector<real>  old_r0155_2;
//   vector<real>  old_r0157_2;
//   vector<real>  old_x0155;
//   vector<real>  old_x0157;
//   real n0_xx[2][mednum_fuel];
//   real np_xx[2][mednum_fuel];
//   real rp_xx[2][mednum_fuel];
//   real rc_xx[2][mednum_fuel];

  
//   real old_xsc5[mednum_fuel];
//   real old_xsc7[mednum_fuel];
//   real old_xsc5_2[mednum_fuel];
//   real old_xsc7_2[mednum_fuel];
//   real old_xscc_5[mednum_fuel];
//   real old_xscc_7[mednum_fuel];
//   real old_rc_5[mednum_fuel];
//   real old_rc_7[mednum_fuel];
//   real tmp_5[mednum_fuel];
//   real tmp_7[mednum_fuel];

//   real alpha_5[mednum_fuel];
//   real alpha_7[mednum_fuel];
//   real time;
//   real time_be;
//   real time_be_2;
  
//   //sasuga addition
//   old_n0155.resize(mednum_fuel);
//   old_n0157.resize(mednum_fuel);
//   old_r0155.resize(mednum_fuel);
//   old_r0157.resize(mednum_fuel);
//   old_np155.resize(mednum_fuel);
//   old_rp155.resize(mednum_fuel);
//   old_rp157.resize(mednum_fuel);
//   old_rc155.resize(mednum_fuel);
//   old_rc157.resize(mednum_fuel);
//   old_n0155_2.resize(mednum_fuel);
//   old_r0155_2.resize(mednum_fuel);
//   old_r0157_2.resize(mednum_fuel);
//   old_x0155.resize(mednum_fuel);
//   old_x0157.resize(mednum_fuel);

//   if(corrector_calc){
//     rr_gd5.resize(mednum_fuel);
//     rr_gd7.resize(mednum_fuel);
//     np_gd5.resize(mednum_fuel);
//     np_gd7.resize(mednum_fuel);
//   };


//   // +++ Initial number density setting +++
//   for(int i=0;i<mednum_fuel;i++){
//     fwd_nuc[0][0][i].set_zero();
//     for(int j=0;j<init_nucnum[i];j++){
//       int idtmp=init_nucid[i][j];
//       real dtmp=init_nucden[i][j];
//       int idpos=med[0].SearchNuclide(idtmp);
//       if(idpos==-1){
//         cout<<"# Error !!\n";
//         exit(0);
//       };
//       fwd_nuc[0][0][i].put_data(idpos,dtmp);
//     };
//   };


//   // +++ Initial heavy metal weight calculation +++
//   for(int i=0;i<mednum_fuel;i++){
//     PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
//     hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
//   };

//   if(med_target!=-1){
//     hm_weight_init=hm_weight_init_per_medium[med_target];
//   }else{
//     hm_weight_init=0.;
//     for(int i=0;i<mednum_fuel;i++){
//       hm_weight_init+=hm_weight_init_per_medium[i];
//     };
//   };

//   cout.setf(ios::scientific);
//   cout.precision(5);
//   cout<<"#\n# Initial heavy metal weight [g]\n#\n";
//   cout<<"#     Total : "<<hm_weight_init<<"\n";
//   for(int i=0;i<mednum_fuel;i++){
//     cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
//   };
//   cout<<"#\n";

//   if(input_power_unit=="MW_t"){
//     input_power_unit="W_cm";
//     for(int i=0;i<burn_step;i++){
//       power_density_list[i]*=hm_weight_init;
//     };
//   };

//   // +++ Burnup calculation condition setting +++
//   PreCalculation_bt();


//   // ... Preparing medium class instances for ADTS model
//   //
//   //
//   vector<Medium> med_adts(mednum_fuel);

//   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   //
//   //  RUNNING BURNUP CALCULATION 
//   //
//   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//   real accumulated_day=0.;
//   real accumulated_burn=0.;
//   vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

//   double dxsf[mednum_fuel][nucn];
//   double dxsc[mednum_fuel][nucn];
//   vector<GroupData1D> pre_xsf_1g_smlux2;
//   pre_xsf_1g_smlux2.resize(mednum);
//   double dxsn2n[mednum_fuel][nucn];
//   int st_10;
//   for(int st=0;st<burn_step+1;st++){
//     int bstmp=0;

//     acday.push_back(accumulated_day);
//     acburn.push_back(accumulated_burn);
//     for(int i=0;i<mednum_fuel;i++){
//       acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
//     };

//     cout<<"#\n# +++ Burnup step : "<<st<<"\n";
//     cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    
//     if(st==0||st%sm_step==0){

//       st_10=st;      

//       // +++ Instance generation of MECSystem to solve neutron transport equation
//       MECSystem lat(group,mednum);
//       lat.PutTrajectorySet(&sys_f);  

//       // +++ fuel medium data preparation
//       for(int i=0;i<mednum_fuel;i++){

//         // (Self-shielding calculation)
//         PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
//         opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
//         if(dancoff_input){
//           for(int g=0;g<group;g++){
//             dancoff.put_data(g,1.-dancoff_factor[i]);
// 	  };
//         };
//         opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
//         opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library
// 	med_adts[i]=med[0];  // ... adts
//         med[0].CalMacroFromMicro();

//         // (Macro&Micro data storing)
//         for(int k=0;k<nucn;k++){
//           if(nuclide_info[k]!=0){
//             if(nuclide_info[k]==1){
//               mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
//             };
//             mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
//             mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
//           };
//         };

//         lat.AddMedium(med[0]); 
//         lat.GetMedium(i).NuclideClear();

//       };

//       // +++ non-fuel medium data preparation
//       for(int jj=0;jj<mednum_nonfuel;jj++){
//         lat.AddMedium(med[1+jj]);
//       };

//       // +++ eigenvalue calculation
//       lat.PutRegMed(region_medium);
//       lat.PutGeneralOption(opt);
//       lat.PutThermalIteration(3);
//       lat.PutPL(0);
//       lat.NoCMRAcceleration();
//       lat.Print();
      
//       keff[st]=lat.CalIgen();

//       for(int i=0;i<mednum;i++){
//         GroupData1D flx=lat.GetIntegratedFlux(i);
//         flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
//         // +++ One-group cross section storing
//         if(i<mednum_fuel){
//           for(int j=0;j<nucn;j++){
//             if(nuclide_info[j]!=0){
//               if(nuclide_info[j]==1){
//    	        xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
// 	      };
// 	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
// 	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
//             };
// 	  };
//         };
//       };

//     }else{

//       // One-group cross section at the last burnup step is assumed
//       for(int i=0;i<mednum;i++){
//         if(i<mednum_fuel){
// 	  for(int j=0;j<nucn;j++){
// 	    if(nuclide_info[j]!=0){
// 	      if(nuclide_info[j]==1){
// 	        xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
// 	      };
// 	      xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
// 	      xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
// 	    };
// 	  };
// 	};
//       };
      
//     };
      
    
//     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//     // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     if(st!=burn_step){


//     #if 0
//       // ... No group collapsing ...
//       // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
//       MECSystem lat2(group,mednum);
//       lat2.PutTrajectorySet(&sys_f2);        

//       // +++ fuel medium data preparation
//       for(int i=0;i<mednum_fuel;i++){

// 	/*
//         // (Self-shielding calculation) -----------------------------------------------------------
//         PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
//         opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
//         if(dancoff_input){
//           for(int g=0;g<group;g++){
//             dancoff.put_data(g,1.-dancoff_factor[i]);
// 	  };
//         };
//         opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
//         opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
//         med[0].CalMacroFromMicro();
// 	// ----------------------------------------------------------------
// 	*/

// 	// (Use of the original burnup step results) -------------------------------
//         int nucn=med_adts[i].GetNucnum();
//         for(int ii=0;ii<nucn;ii++){
//            med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
//         };
// 	med_adts[i].CalMacroFromMicro();
// 	med[0]=med_adts[i];
// 	// -------------------------------------------------------------------------
	
	
//         // (Macro&Micro data storing)
//         for(int k=0;k<nucn;k++){
//           if(nuclide_info[k]!=0){
//             if(nuclide_info[k]==1){
//               mic_sigf_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
//             };
//             mic_sigc_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
//             mic_sign2n_2[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
//           };
//         };

//         lat2.AddMedium(med[0]); 
//         lat2.GetMedium(i).NuclideClear();

//       };

//       // +++ non-fuel medium data preparation
//       for(int jj=0;jj<mednum_nonfuel;jj++){
//         lat2.AddMedium(med[1+jj]);
//       };

//       // +++ eigenvalue calculation
//       lat2.PutRegMed(region_medium);
//       lat2.PutGeneralOption(opt);
//       lat2.PutThermalIteration(3);
//       lat2.PutPL(0);
//       lat2.NoCMRAcceleration();
//       lat2.Print();

//       real keff_dummy=lat2.CalIgen();

//       vector<real> medium_wise_flux_change(mednum);
      
//       for(int i=0;i<mednum;i++){
//         GroupData1D flx_2=lat2.GetIntegratedFlux(i);

// 	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
// 	medium_wise_flux_change[i]=(flx_2.get_sum()/vol_med[i])/flx_med_2[i].get_sum();

//         flx_med_2[i]=flx_2*(1./vol_med[i]); // medium-wise neutron flux data is stored here
//         // +++ One-group cross section storing
//         if(i<mednum_fuel){
//           for(int j=0;j<nucn;j++){
//             if(nuclide_info[j]!=0){
//               if(nuclide_info[j]==1){
//    	        xsf_1g_2[st][i][j]=mic_sigf_2[bstmp][i][j].Cond(flx_2);
//   	      };
// 	      xsc_1g_2[st][i][j]=mic_sigc_2[bstmp][i][j].Cond(flx_2);
// 	      xsn2n_1g_2[st][i][j]=mic_sign2n_2[bstmp][i][j].Cond(flx_2);
//             };
//   	  };
//         };
//       };
//     #endif      


//     #if 1
//       // ... with group collapsing from 107-group to fewer-group ...
//       // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
      

//       //int ngrp=21; // The number of energy groups after the group collapsing
//       //int bgrp[]={4,9,14,19,24, 29,34,39,44,49, 54,59,64,69,74, 79,84,89,94,99, 106};
//       // Energy group boundary in the collapsed structure.
//       // The final one should be [106] when the number of groups in the original structure is the 107-group.
      
//       // int ngrp=11;
//       // int bgrp[]={9,19,29,39,49, 59,69,79,89,99, 106};
//       /*
//       int ngrp=107;
//       int *bgrp=new int[ngrp];
//       for(int i=0;i<ngrp;i++){
// 	bgrp[i]=i;
//       };
//       */

//       vector< vector<GroupData1D> > mic_sigf_3(mednum_fuel);  // [medid][nucn](group)
//       vector< vector<GroupData1D> > mic_sigc_3(mednum_fuel);
//       vector< vector<GroupData1D> > mic_sign2n_3(mednum_fuel); 
	
//       MECSystem lat2(ngrp,mednum);
//       lat2.PutTrajectorySet(&sys_f2);     


//   // +++ GPT-PC +++++++++++++++++++++++++++++++++++++++++++++++++
//   if(adjoint&&corrector_calc){
//     fwd_nuc_p.resize(burn_step+1);
//     fwd_nuc_p_int.resize(burn_step+1);
//     for(int i=0;i<burn_step+1;i++){
//       int sub_step=sub_step_list[i];
//       fwd_nuc_p[i].resize(sub_step+1);
//       fwd_nuc_p_int[i].resize(sub_step+1);
//       for(int j=0;j<sub_step+1;j++){
//         fwd_nuc_p[i][j].resize(mednum_fuel);
//         fwd_nuc_p_int[i][j].resize(mednum_fuel);
//         for(int k=0;k<mednum_fuel;k++){
//           fwd_nuc_p[i][j][k].put_imax(nucn);
//           fwd_nuc_p_int[i][j][k].put_imax(nucn);
//         };
//       };
//     };
//     volflx_mesh_p.resize(burn_step);  
//     power_factor_p.resize(burn_step);
//     for(int i=0;i<burn_step;i++){
//       volflx_mesh_p[i].resize(totm);
//       for(int j=0;j<totm;j++){
//         if(region_medium[j]<mednum_fuel){
//           volflx_mesh_p[i][j].put_imax(group);
//         };
//       };
//       int sub_step=sub_step_list[i];
//       power_factor_p[i].resize(sub_step+1);
//     };
//     int totsn=lat.GetQuad().GetSN();

//     if(aflx_legendre==-1){
//       volaflx_mesh_p.resize(burn_step);
//       for(int i=0;i<burn_step;i++){
//         volaflx_mesh_p[i].resize(totm);
//         for(int j=0;j<totm;j++){
// 	  if(region_medium[j]<mednum_fuel){
// 	    volaflx_mesh_p[i][j].resize(group);
// 	    for(int k=0;k<group;k++){
// 	      volaflx_mesh_p[i][j][k].put_imax(totsn);
// 	    };
// 	  };
//         };
//       };
//     }else{
//       volaflx_pl_p.resize(burn_step);
//       for(int i=0;i<burn_step;i++){
//         volaflx_pl_p[i].resize(totm);
//         for(int j=0;j<totm;j++){
//     	  if(region_medium[j]<mednum_fuel){
// 	    volaflx_pl_p[i][j].resize(group);
// 	    for(int k=0;k<group;k++){
// 	      volaflx_pl_p[i][j][k].put_imax(quad.GetPlnum());
// 	    };
// 	  };
//         };
//       };
//     };
//   };
//   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

//       // +++ fuel medium data preparation
//       for(int i=0;i<mednum_fuel;i++){

// 	/*
//         // (Self-shielding calculation) ---------------------------------------
//         PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
//         opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
//         if(dancoff_input){
//           for(int g=0;g<group;g++){
//             dancoff.put_data(g,1.-dancoff_factor[i]);
// 	  };
//         };
//         opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
//         opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
//         med[0].CalMacroFromMicro();

//         Medium bmed=med[0].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
//         // ------------------------------------------------------------------------
// 	*/

// 	// (Use of the original burnup step results) -------------------------------
//         int nucn=med_adts[i].GetNucnum();
//         for(int ii=0;ii<nucn;ii++){
//            med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
//         };
// 	med_adts[i].CalMacroFromMicro();
// 	Medium bmed=med_adts[i].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
// 	// -------------------------------------------------------------------------


//         // (Macro&Micro data storing)
// 	mic_sigf_3[i].resize(nucn);
// 	mic_sigc_3[i].resize(nucn);
// 	mic_sign2n_3[i].resize(nucn);	
//         for(int k=0;k<nucn;k++){
	  
//           if(nuclide_info[k]!=0){
//             if(nuclide_info[k]==1){
//               mic_sigf_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
//             };
//             mic_sigc_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
//             mic_sign2n_3[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
//           };
//         };

//         lat2.AddMedium(bmed); 
//         lat2.GetMedium(i).NuclideClear();

//       };

//       // +++ non-fuel medium data preparation
//       for(int jj=0;jj<mednum_nonfuel;jj++){
//         Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel*jj],true); // [true] is for microscopic XS
//         lat2.AddMedium(bmed);
//       };

//       // +++ eigenvalue calculation
//       lat2.PutRegMed(region_medium);
//       lat2.PutGeneralOption(opt);
//       lat2.PutThermalIteration(3);
//       lat2.PutPL(0);
//       lat2.NoCMRAcceleration();
//       lat2.Print();

//       real keff_dummy=lat2.CalIgen();

//       vector<real> medium_wise_flux_change(mednum);
      
//       for(int i=0;i<mednum;i++){
//         GroupData1D flx_2=lat2.GetIntegratedFlux(i);

// 	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
// 	medium_wise_flux_change[i]=(flx_2.get_sum()/vol_med[i])/flx_med_2[i].get_sum();

//         flx_med_2[i]=flx_2*(1./vol_med[i]); // medium-wise neutron flux data is stored here
//         // +++ One-group cross section storing
//         if(i<mednum_fuel){
//           for(int j=0;j<nucn;j++){
//             if(nuclide_info[j]!=0){
//               if(nuclide_info[j]==1){
//    	        xsf_1g_2[st][i][j]=mic_sigf_3[i][j].Cond(flx_2);
//   	      };
// 	      xsc_1g_2[st][i][j]=mic_sigc_3[i][j].Cond(flx_2);
// 	      xsn2n_1g_2[st][i][j]=mic_sign2n_3[i][j].Cond(flx_2);
//             };
//   	  };
//         };
//       };
//     #endif     

//     abs_frac[st]=sum_fmed/sum_tot; // absorption rate fraction ini fuel media
//     cout<<"# Absorption fraction : "<<abs_frac[st]<<"\n"; 

//     /*      
//       // Commented out in this method since it is unnecessary.
//       for(int i=0;i<totm;i++){
//         if(region_medium[i]<mednum_fuel){
//           real vol=lat.GetMesh(i).GetVolume();
//           volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
// 	};
//       };
//     */

//       real power_density=power_density_list[st];
//       real burn_span=burn_time[st]; // day
//       int sub_step=sub_step_list[st];
//       burn_span/=sub_step;   

//       cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
    
//       for(int j=0;j<sub_step;j++){

// 	// (Line power of target medium is calculated)
// 	real sumflx=0.;
// 	real power_org=0.;

//         vector<real> power_per_medium(mednum_fuel);
// 	for(int i=0;i<mednum_fuel;i++){
//           real tmp=flx_med[i].get_sum();
//           //real tmp=flx_med_2[i].get_sum();	  
//           power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
// 	};

// 	if(med_target!=-1){
//   	  sumflx=flx_med[med_target].get_sum();
//   	  //sumflx=flx_med_2[med_target].get_sum();	  
//           power_org=power_per_medium[med_target];
// 	}else{
//           power_org=0.;
// 	  for(int i=0;i<mednum_fuel;i++){
//             power_org+=power_per_medium[i];
// 	  };
// 	};

// 	if(input_flux_level){
//           power_factor[st][j]=flux_level_list[st]/sumflx;
// 	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
// 	}else{
//   	  power_factor[st][j]=power_density/power_org;
//           accumulated_burn+=burn_time_gwd[st]/sub_step;
// 	};

//         for(int i=0;i<mednum_fuel;i++){
// 	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
// 	};
//         accumulated_day+=burn_span;
//         delt[st][j]=burn_span*24*60*60;
	
// 	for(int i=0;i<mednum_fuel;i++){
// 	  total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
// 	  //total_flux[st][j][i]=flx_med_2[i].get_sum()*power_factor[st][j];	  
// 	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
// 	};
	
//       }; // end of sub-step loop

//       fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

//     // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //   CORRECTOR CALCULATION
//       if(corrector_calc){

// 	cout<<"# Corrector calculation in ADTS model ...\n";

// 	// (OWPC)	
// 	for(int i=0;i<mednum_fuel;i++){
// 	  real tmp=total_flux[st][0][i];
//   	  rr_gd5[i]=tmp*xsc_1g[st][i][pos_gd155];  // Reaction rate at BOC (Rp)
//   	  rr_gd7[i]=tmp*xsc_1g[st][i][pos_gd157];
// 	};

// 	// ++++++++++++++++++++++++

// 	// (Predictor calculation results are stored in the array of [XXX_p].)
//         swap(total_flux_p[st], total_flux[st]);
//         swap(xsc_1g_p[st], xsc_1g[st]);
//         swap(xsn2n_1g_p[st], xsn2n_1g[st]);
//         swap(xsf_1g_p[st], xsf_1g[st]);

// 	if(adjoint){
//           swap(keff_p[st], keff[st]);
//           swap(fwd_nuc_p[st], fwd_nuc[st]);
//           swap(fwd_nuc_p_int[st], fwd_nuc_int[st]);
//           swap(power_factor_p[st], power_factor[st]);
// 	  swap(volflx_mesh_p[st], volflx_mesh[st]);
// 	  if(aflx_legendre==-1){
//             swap(volaflx_mesh_p[st], volaflx_mesh[st]);
// 	  }else{
//   	    swap(volaflx_pl_p[st], volaflx_pl[st]);
// 	  };
// 	  swap(macxs_p[st], macxs[st]);
// 	  fwd_nuc[st][0]=fwd_nuc_p[st][0];
// 	};

//       for(int i=0;i<mednum_fuel;i++){
//         // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
//         PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
//         opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
//         if(dancoff_input){
//           for(int g=0;g<group;g++){
//             dancoff.put_data(g,1.-dancoff_factor[i]);
//   	  };
//           opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         }else{
//           opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         };
//         //opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
//         //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
//         opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library	 
//         med[0].CalMacroFromMicro();
//         //opc.CalFissionSpectrumMatrix(med[0],xslib);

//         // +++ Macro&Micro data storing
//         for(int k=0;k<nucn;k++){
//           if(nuclide_info[k]!=0){
//             if(nuclide_info[k]==1){
//               mic_sigf_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
//             };
//             mic_sigc_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
//             mic_sign2n_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
//           };
//         };
// 	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
// 	if(adjoint){
//           macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
// 	};
//       };

//       for(int jj=0;jj<mednum_nonfuel;jj++){
//         lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
//       };

//       lat.PutThermalIteration(3);
//       lat.PutPL(0);
//       lat.NoCMRAcceleration();
//       real keff_corr=lat.CalIgen();
      
//       // ------------------------------

//       if(adjoint)keff[st]=keff_corr;

//       for(int i=0;i<mednum;i++){
//         GroupData1D flx=lat.GetIntegratedFlux(i);
//         flx_med_c[i]=flx*(1./vol_med[i]);
//         // +++ One-group cross section storing
//         if(i<mednum_fuel){
//           for(int j=0;j<nucn;j++){
//             if(nuclide_info[j]!=0){
//               if(nuclide_info[j]==1)xsf_1g[st][i][j]=mic_sigf_c[i][j].Cond(flx);
// 	      xsc_1g[st][i][j]=mic_sigc_c[i][j].Cond(flx);
// 	      xsn2n_1g[st][i][j]=mic_sign2n_c[i][j].Cond(flx);
// 	    };
//           };
// 	};
//       };

//       //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//       if(st==0||st%sm_step==0){	
//         //st_10=st;
//       }else{

//         // ... One-group cross section correction in ADTS	
	
//         for(int i=0;i<mednum_fuel;i++){
// 	  for(int k=0;k<nucn;k++){
	  
//   	    int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
//     #if 1	  
//             // ... Copy of the ADTS step data
// 	    if(nuclide_info[k]!=0){
// 	      if(nuclide_info[k]==1){
//   	        xsf_1g[st][i][k]=xsf_1g_2[st][i][k];
// 	      };
// 	      xsc_1g[st][i][k]=xsc_1g_2[st][i][k];
// 	      if(xsn2n_1g[st-1][i][k]>0.){
//   	        xsn2n_1g[st][i][k]=xsn2n_1g_2[st][i][k];
// 	      };
// 	    };
//     #endif	  

//     #if 0	  
//             // ... Correcting all nuclides in all burnup regions
// 	    if(nuclide_info[k]!=0){
// 	      if(nuclide_info[k]==1){
//   	        xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
// 	      };
// 	      xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
// 	      if(xsn2n_1g[st-1][i][k]>0.){
//   	        xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
// 	      };
// 	    };
//     #endif	  

//     #if 0	  
// 	    if(i<8){
//               // ... Correcting all nuclides in the central Gd-bearing pin
// 	      if(nuclide_info[k]!=0){
// 	        if(nuclide_info[k]==1){
//   	          xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
// 	        };
// 	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
//   	        if(xsn2n_1g[st-1][i][k]>0.){	      
// 	          xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
// 	        };
// 	      };
// 	    };
//     #endif	  

//     #if 0	  
// 	    if(i<8){
// 	      // ... Correcting Gd-155 and -157 capture cross sections in the central Gd-bearing pin
// 	      if(nucid==641550||nucid==641570){
// 	        if(nucid==641550) cout<<"Gd-155 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
// 	        if(nucid==641570) cout<<"Gd-157 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
// 	        //xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
// 	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_2[st][i][k]/xsc_1g_2[st-1][i][k]);//for neutron capture cross section
// 	        //xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
// 	      };
// 	    };
//     #endif	  

// 	  };
//         };

// 	// ... Neutron flux correction in ADTS

// 	// for(int i=0;i<mednum_fuel;i++){ // All the medium included in the system
// 	//   //for(int i=0;i<8;i++){ // Only for Gd-bearing rod
// 	//   flx_med[i]=flx_med[i]*medium_wise_flux_change[i];
// 	// };
//       };

//       //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//       if(adjoint){
//         for(int i=0;i<totm;i++){
//           if(region_medium[i]<mednum_fuel){
//             real vol=lat.GetMesh(i).GetVolume();
//             volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol;
// 	    if(aflx_legendre==-1){
//               for(int g=0;g<group;g++){
// 	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
// 	      };
// 	    }else{
//               for(int g=0;g<group;g++){
//                 int sntot=quad.GetSN();	      
// 	        for(int lm=0;lm<quad.GetPlnum();lm++){
//                   real tmp=0.;
//   	          for(int w=0;w<sntot;w++){
// 	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
// 		  };
// 		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
// 	        };
// 	      };
// 	    };

// 	  };
//         };


//       };

//       real acburn_pre=accumulated_burn-acburn[st];
//       vector<real> acburn_pre_per_medium(mednum_fuel);
//       for(int j=0;j<mednum_fuel;j++){
// 	acburn_pre_per_medium[j]=accumulated_burn_per_medium[j]-acburn_per_medium[j][st];
//       };
//       // accumulated burnup calculated by the predictor step

      
//       // 
//       // .....  One-group cross section correction for Rc in A-OWPC .....
//       //

      
//       if(owpc_corr){
	
//       for(int i=0;i<mednum_fuel;i++){

// 	real n0=fwd_nuc[st][0][i].get_dat(pos_gd157);   // initial 
//         real np=fwd_nuc[st+1][0][i].get_dat(pos_gd157); // predictor results      
//         if((n0-np)/np>1e-2){
//           // This if branch is added in 2021/12/4 since negative cross section is detected in the Gd cross section.
// 	  // This should be consistent in the following process in the actual correction to the final ND.
	
// 	real xsc0_5=xsc_1g_p[st][i][pos_gd155]; //rp
// 	real xsc0_7=xsc_1g_p[st][i][pos_gd157]; //rp
// 	real xscc_5=xsc_1g[st][i][pos_gd155];//Rcに相当するXc
// 	real xscc_7=xsc_1g[st][i][pos_gd157];
// 	if(st>0){
// 	  real n0_5=fwd_nuc[st][0][i].get_dat(pos_gd155);
// 	  real np_5=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
// 	  real t_5=old_xsc5[i]+(xsc0_5-old_xsc5[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);//前のステップでの正しいXcを計算
// 	  real t_7=old_xsc7[i]+(xsc0_7-old_xsc7[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
// 	  if(st>1){
// 	    real a=n0_5, b=old_n0155[i], c=old_n0155_2[i];
// 	    real aa_5=xsc0_5, bb_5=old_xsc5[i], cc_5=old_xsc5_2[i];
// 	    real aa_7=xsc0_7, bb_7=old_xsc7[i], cc_7=old_xsc7_2[i];
// 	    real X_5=old_np155[i];
// 	    t_5=aa_5*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_5*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_5*(X_5-a)*(X_5-b)/(c-a)/(c-b);//前のステップでの正しいXcを計算
// 	    t_7=aa_7*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_7*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_7*(X_5-a)*(X_5-b)/(c-a)/(c-b);
// 	  };

// 	  // correction without any consideration on time step difference
// 	  /*
// 	  xsc_1g[st][i][pos_gd155]=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5);
// 	  xsc_1g[st][i][pos_gd157]=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7);
// 	  */
// 	  /*
//           // correction with considering the time step length difference
// 	  xsc_1g[st][i][pos_gd155]=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)/burn_time[st-1]*burn_time[st];
// 	  xsc_1g[st][i][pos_gd157]=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)/burn_time[st-1]*burn_time[st];
// 	  */

// 	  // correction with considering time step length x power density (= burnup) 
// 	  real tmp1=power_density_list[st-1]*burn_time[st-1];
// 	  real tmp2=power_density_list[st]*burn_time[st];
// 	  real factor=tmp2/tmp1;
// 	  real xsc155=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)*factor;
// 	  real xsc157=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)*factor;
// 	  if(factor<10.){ // if the burnup length is too different, the correction is NOT applied.
// 	    if(xsc155>0.)xsc_1g[st][i][pos_gd155]=xsc155;
// 	    if(xsc157>0.)xsc_1g[st][i][pos_gd157]=xsc157;
// 	  };
	  
// 	  old_xsc5_2[i]=old_xsc5[i];
// 	  old_xsc7_2[i]=old_xsc7[i];
// 	};
// 	old_xsc5[i]=xsc0_5;
// 	old_xsc7[i]=xsc0_7;
// 	old_xscc_5[i]=xscc_5;
// 	old_xscc_7[i]=xscc_7;
//       };

//       }; // end of [if((n0-np)/np>1e-2)]
//       }; // end of [if(owpc_corr)]
//       // ............................................................................
      
//       for(int j=0;j<sub_step;j++){

// 	// (Line power of target medium is calculated)
// 	real sumflx=0.;
// 	real power_org=0.;

//         vector<real> power_per_medium(mednum_fuel);
// 	for(int i=0;i<mednum_fuel;i++){
//           real tmp=flx_med_c[i].get_sum();
//           power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
// 	};

// 	if(med_target!=-1){
//   	  sumflx=flx_med_c[med_target].get_sum();
//           power_org=power_per_medium[med_target];
// 	}else{
//           power_org=0.;
// 	  for(int i=0;i<mednum_fuel;i++){
//             power_org+=power_per_medium[i];
// 	  };
// 	};

// 	if(input_flux_level){
//           power_factor[st][j]=flux_level_list[st]/sumflx;
// 	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
// 	}else{
//   	  power_factor[st][j]=power_density/power_org;
//           //accumulated_burn+=burn_time_gwd[st]/sub_step;
// 	};

//         for(int i=0;i<mednum_fuel;i++){
// 	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
// 	};

//         //accumulated_day+=burn_span;
//         //delt[st][j]=burn_span*24*60*60;
//         for(int i=0;i<mednum_fuel;i++){
// 	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
// 	  //CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.

// 	  /*
// 	  if(i==40){
// 	    cout<<"# Medium : "<<i<<"\n";
// 	    for(int j=0;j<nucn;j++){
//               cout<<j<<" "<<med[0].GetNuclideInTurn(j).GetMatnum()<<" "<<xsc_1g[st][i][j]<<" "<<xsn2n_1g[st][i][j]<<"\n";
// 	    };
// 	    cout<<" total_flux : "<<total_flux[st][j][i]<<"\n";
// 	    cout<<" delt       : "<<delt[st][j]<<"\n";	    
// 	  };
// 	  */

// 	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
// 	};
	
//       }; // end of sub-step
      
//       // (OWPC)
//       real dt=burn_time[st]*60*60*24;
//       real time_mesh=10000;
//       real dt_r=dt/time_mesh;
//       real omega_155, omega_157; // various correlation conditions

//       for(int i=0;i<mednum_fuel;i++){

// 	  // OWPC treatment for Gd-155 and -157 with various correlation conditions
// 	if(owpc_corr){

// 	  if(st>0){ // ... quadratic model
// 	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
// 	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
// 	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
// 	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
// 	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
// 	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
// 	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
// 	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
// 	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
// 	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
// 	    real x0_155=xsc_1g_p[st][i][pos_gd155];//rpに対応する断面積
// 	    real x0_157=xsc_1g_p[st][i][pos_gd157];//rpに対応する断面積
// 	    real xp_155=xsc_1g[st][i][pos_gd155]; //rcに対応する断面積
// 	    real xp_157=xsc_1g[st][i][pos_gd157]; //rcに対応する断面積
// 	    real n0_r_155;
// 	    real n0_r_157;
	    
// 	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
// 	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt);//Not consider 154,156
	    	    
// 	    //quadratic model---------------------------
// 	    real a_5,b_5,c_5,a_7,b_7,c_7;
// 	    real aa_5,bb_5,cc_5,aa_7,bb_7,cc_7;
// 	    a_5=n0_155;
// 	    b_5=np_155;
// 	    c_5=old_n0155[i];
// 	    a_7=n0_157;
// 	    b_7=np_157;
// 	    c_7=old_n0157[i];
// 	    real X_5=np_p_155;
// 	    real X_7=np_p_157;

// 	    //-------------------------reaction rate base-------------------------------------
// 	    aa_5=r0_155;
// 	    bb_5=rp_155;
// 	    //cc_5=old_r0155[i];
// 	    cc_5=old_r0155[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba    
	    
// 	    aa_7=r0_157;
// 	    bb_7=rp_157;
// 	    //cc_7=old_r0157[i];
// 	    cc_7=old_r0157[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba

// 	    real rc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
// 	    real rc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
// 	    //----------------------------------------------------------------------------------------

// 	    /*
// 	    //---------------------------------------cross section base-------------------------------------
// 	    aa_5=x0_155;
// 	    bb_5=xp_155;
// 	    cc_5=old_x0155[i];
	    
// 	    aa_7=x0_157;
// 	    bb_7=xp_157;
// 	    cc_7=old_x0157[i];
	        
// 	    real xc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //xc in toyproblem
// 	    real xc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
// 	    //----------------------------------------------------------------
// 	    */
	    
// 	    n0_r_155=n0_155;
// 	    real rp_r_155=r0_155;
// 	    real xp_r_155=x0_155;
	  
// 	    n0_r_157=n0_157;
// 	    real rp_r_157=r0_157;
// 	    real xp_r_157=x0_157;

// 	    //------------------------------------start toy problem calculation----------------------------
// 	    for(int k=0;k<time_mesh;k++){


// 	      //-----reaction rate base-------
// 	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
// 	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

// 	      X_5=np_r_155;
// 	      X_7=np_r_157;
	      
// 	      real rc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
// 	      real rc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

// 	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
// 	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

// 	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
// 	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r); //Not consider 154,156
	    
// 	      rp_r_155=rc_r_155;
// 	      rp_r_157=rc_r_157;
// 	      //----------------------------------------

// 	      /*
// 	      //-------cross section base-------
// 	      real np_r_155=n0_r_155*exp(-xp_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
// 	      real np_r_157=n0_r_157*exp(-xp_r_157*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156

// 	      X_5=np_r_155;
// 	      X_7=np_r_157;

// 	      real xc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
// 	      real xc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

// 	      real x_r_155=(xp_r_155+xc_r_155)*0.5;
// 	      real x_r_157=(xp_r_157+xc_r_157)*0.5;

// 	      n0_r_155=n0_r_155*exp(-x_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
// 	      n0_r_157=n0_r_157*exp(-x_r_157*total_flux_p[st][0][i]*1e-24*dt_r); //Not consider 154,156
	    
// 	      xp_r_155=xc_r_155;
// 	      xp_r_157=xc_r_157;
// 	      //----------------------------------------
// 	      */
	      
// 	    };
// 	    //-----------------------------------------end toy problem calculation------------------------------------------
	    
// 	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
// 	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

// 	    //-------reaction rate base----------
// 	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
// 	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
// 	    //------------------------------------------

// 	    /*
// 	    //-------cross section base----------
// 	    omega_155=(R_reference_155-x0_155*total_flux_p[st][0][i])/(xc_155*total_flux_p[st][0][i]-x0_155*total_flux_p[st][0][i]);
// 	    omega_157=(R_reference_157-x0_157*total_flux_p[st][0][i])/(xc_157*total_flux_p[st][0][i]-x0_157*total_flux_p[st][0][i]);
// 	    //--------------------------------------
// 	    */
	    
// 	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
// 	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;
	    
// 	    old_n0155_2[i]=old_n0155[i];
// 	    old_r0155_2[i]=old_r0155[i];
// 	    old_r0157_2[i]=old_r0157[i];
	    
// 	    old_n0155[i]=n0_155;
// 	    old_n0157[i]=n0_157;
// 	    old_r0155[i]=r0_155;
// 	    old_r0157[i]=r0_157;
// 	    old_x0155[i]=x0_155;
// 	    old_x0157[i]=x0_157;
	    
// 	    old_np155[i]=np_155;
// 	    old_rc155[i]=rp_155;
// 	    old_rc157[i]=rp_157;
	      	    
// 	  }else{ // ... linear model for [st]==0
	  
// 	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
// 	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
// 	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
// 	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
// 	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
// 	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
// 	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
// 	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
// 	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
// 	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	
// 	    real n0_r_155;
// 	    real n0_r_157;

// 	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
// 	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt); //Not consider 154,156
	  
// 	    // (Correlation to Gd-155 number density)
// 	    real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
// 	    real alpha_157=(rp_157-r0_157)/(np_155-n0_155);    
	    
// 	    // (Correlation to Gd-155 number density)	
// 	    real rc_155=r0_155+(np_p_155-n0_155)*alpha_155; // Rc in toy-problem
// 	    real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	   
// 	    n0_r_155=n0_155;
// 	    real rp_r_155=r0_155;
	  
// 	    n0_r_157=n0_157;
// 	    real rp_r_157=r0_157;
	  
// 	    for(int k=0;k<time_mesh;k++){

// 	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
// 	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

// 	      // (Correlation to Gd-155 number density)
// 	      real rc_r_155=r0_155+(np_r_155-n0_155)*alpha_155;
// 	      real rc_r_157=r0_157+(np_r_155-n0_155)*alpha_157;
	
// 	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
// 	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

// 	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
// 	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r);//Not consider 154,156

// 	      rp_r_155=rc_r_155;
// 	      rp_r_157=rc_r_157;

// 	    };
	        
// 	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
// 	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

// 	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
// 	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);

// 	    //cout<<"#   Omega_155 / Omega_157 : "<<omega_155<<" "<<omega_157<<"\n";
	    
// 	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
// 	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;

// 	    old_n0155[i]=n0_155;
// 	    old_n0157[i]=n0_157;
// 	    old_r0155[i]=r0_155;
// 	    old_r0157[i]=r0_157;

// 	    old_np155[i]=np_155;
// 	    old_rc155[i]=rp_155;
// 	    old_rc157[i]=rp_157;

// 	  }; //

// 	}; // end of [if(owpc_corr)]
	  
// 	for(int j=0;j<nucn;j++){
	
// 	  real np=fwd_nuc[st+1][0][i].get_dat(j);      // predictor results
//   	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
// 	  real n_next=exp((log(np)+log(nc))*0.5);
//           int nucid=med[0].GetNuclideInTurn(j).GetMatnum();

// 	  if(nucid==641550||nucid==641570){

//   	    real n0=fwd_nuc[st][0][i].get_dat(j);
//             if((n0-np)/np>1e-2){

// 	      if(owpc_corr){
// 		if(nucid==641550){
// 		  n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
// 		}else if(nucid==641570){
// 		  n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
//     	        };
// 	      }else{
// 		/*
// 		// -- OWPC w/o (N,R) correlation ----------------------- 
// 		real r0=rr_gd5[i]; //Rp
// 		if(nucid==641570)r0=rr_gd7[i];

// 		real rp=xsc_1g[st][i][j]*total_flux[st][0][i]; //Rc
// 		real alpha=(rp-r0)/(np-n0);
		  
// 		real np_p=n0*exp(-r0*1e-24*dt);   // Nc in toy-problem
// 		real rc=alpha*np_p+(r0-n0*alpha); // Rc in toy-problem
	    
// 		real n0_r=n0;
// 		real rp_r=r0;
// 		for(int k=0;k<time_mesh;k++){
// 		  real np_r=n0_r*exp(-rp_r*1e-24*dt_r);
// 		  real rc_r=alpha*np_r+(r0-n0*alpha);
// 		  real r_r=(rp_r+rc_r)/2;
// 		  n0_r=n0_r*exp(-r_r*1e-24*dt_r);
// 		  rp_r=rc_r;
// 		};
		
// 		real R_reference=(log(n0)-log(n0_r))/dt*1e24;
// 		real omega=(R_reference-r0)/(rc-r0);
		  
// 		// -- PPC -------------------------------
// 		real R_p=-log(np/n0)/dt;
// 		real R_c=-log(nc/n0)/dt;
// 		real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
// 		R_reference=(R_p+R_c_c)/2;
// 		omega=(R_reference-R_p)/(R_c-R_p);
// 		// --------------------------------

// 		n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応
// 		*/

//  	        // ... conventional PC ...
//                 //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
//                 //n_next=(np+nc)*0.5;  // Conventional PC (in linear)
//                 //n_next=nc;          
//               };

// 	    };

// 	  };

// 	  /*
//           if(!owpc_corr){
// 	    n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
// 	    cout<<nucid<<" "<<n_next<<" "<<np<<" "<<nc<<" "<<"\n";
// 	  };
// 	  */

//           if(adjoint)        n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
//           if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC

// 	  // (For OWPC)
//     	  if(nucid==641550)np_gd5[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
//   	  if(nucid==641570)np_gd7[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd157);

// 	  /*
// 	  if(i==40){
// 	  cout<<"# Medium : "<<i<<"\n";
//   	    cout<<j<<" "<<med[0].GetNuclideInTurn(j).GetMatnum()<<" "<<n_next<<" "<<np<<" "<<nc<<"\n";	     
// 	  };
// 	  */
	  
//           fwd_nuc[st+1][0][i].put_data(j,n_next);

// 	}; // end of nuclides iteration

//       }; // end of medium iteration

//       // Adjustment of accumulated burn
//       for(int i=0;i<mednum_fuel;i++){
//         real acburn_cor=accumulated_burn_per_medium[i]-acburn_per_medium[i][st]-acburn_pre_per_medium[i];
//         accumulated_burn_per_medium[i]=acburn_per_medium[i][st]+(acburn_pre_per_medium[i]+acburn_cor*wgt_nc)/(1.+wgt_nc);
//       };

//       if(input_flux_level){
//         real acburn_cor=accumulated_burn-acburn[st]-acburn_pre;
//         accumulated_burn=acburn[st]+(acburn_pre+acburn_cor*wgt_nc)/(1.+wgt_nc);
//       };
//       cout<<"#     ... terminated.\n";
//       }; // END OF CORRECTOR CALCULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
//     }; // end part of [If (st!=burn_step)]
  
//   }; // loop-end of burnup step
// };

void MulticellBurner::ForwardCalculationADTS_toOWPC(Burnup &bu, int med_target, int ngrp, int *bgrp)
{
  // based on ForwardCalculation()

  // fuga addition
  // ofstream fout_1g("xs1g.txt");
  // ofstream fout_tf("total_flux.txt");
  // ofstream fout_p("micXS_p");
  // ofstream fout_c("micXS_c");
  // -------------

  bool owpc_corr=true; // If [false], the (log-averaging) PC is used.
  
  // ... Hard-coded parameters for weighted predictor-corrector
  //real wgt_nc=1.2; // relative weight for corrector
  real wgt_nc=1.0; // relative weight for corrector

  // ... OWPC to store positions of Gd-155 and 157
  int pos_gd155, pos_gd157;
  for(int j=0;j<nucn;j++){
    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    if(nucid==641550)pos_gd155=j;
    if(nucid==641570)pos_gd157=j;
  };

  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculation.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting for forward calculation +++

  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  // vector<GroupData1D> flx_med_sm;
  flx_med.resize(mednum);
  // flx_med_sm.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
    // flx_med_sm[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  // +++ Array setting for predictor-corrector calculation
  //
  // - Generally multi-group microscopic cross section data at every burn steps
  //   are NOT stored because of their large required memory, but those are 
  //   required for GPT (adjoint) calculations, so those at every burnup steps
  //   are stored in the array [mic_sigx].
 
  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };


  vector<GroupData1D> flx_med_c;
  if(corrector_calc){
    //int tmp=1;
    int tmp=burn_step+1;
    xsc_1g_p.resize(tmp);
    xsn2n_1g_p.resize(tmp);
    xsf_1g_p.resize(tmp);
    for(int j=0;j<tmp;j++){
      xsc_1g_p[j].resize(mednum_fuel);
      xsn2n_1g_p[j].resize(mednum_fuel);
      xsf_1g_p[j].resize(mednum_fuel);
      for(int i=0;i<mednum_fuel;i++){
        xsc_1g_p[j][i].resize(nucn,0.);
        xsn2n_1g_p[j][i].resize(nucn,0.);
        xsf_1g_p[j][i].resize(nucn,0.);
      };
    };

    flx_med_c.resize(mednum);
    for(int i=0;i<mednum;i++){
      flx_med_c[i].put_imax(group);
    };

    total_flux_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      total_flux_p[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        total_flux_p[i][j].resize(mednum_fuel);
      };
    };

  };

  // ... OWPC
  vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC (Rp)
  vector<real> np_gd5, np_gd7; // Np

  vector<real>  old_n0155;
  vector<real>  old_n0157;
  vector<real>  old_r0155;
  vector<real>  old_r0157;
  vector<real>  old_np155;
  vector<real>  old_rp155;
  vector<real>  old_rp157;
  vector<real>  old_rc155;
  vector<real>  old_rc157;
  vector<real>  old_n0155_2;
  vector<real>  old_r0155_2;
  vector<real>  old_r0157_2;
  vector<real>  old_x0155;
  vector<real>  old_x0157;
  real n0_xx[2][mednum_fuel];
  real np_xx[2][mednum_fuel];
  real rp_xx[2][mednum_fuel];
  real rc_xx[2][mednum_fuel];

  
  real old_xsc5[mednum_fuel];
  real old_xsc7[mednum_fuel];
  real old_xsc5_2[mednum_fuel];
  real old_xsc7_2[mednum_fuel];
  real old_xscc_5[mednum_fuel];
  real old_xscc_7[mednum_fuel];
  real old_rc_5[mednum_fuel];
  real old_rc_7[mednum_fuel];
  real tmp_5[mednum_fuel];
  real tmp_7[mednum_fuel];

  real alpha_5[mednum_fuel];
  real alpha_7[mednum_fuel];
  real time;
  real time_be;
  real time_be_2;
  
  //sasuga addition
  old_n0155.resize(mednum_fuel);
  old_n0157.resize(mednum_fuel);
  old_r0155.resize(mednum_fuel);
  old_r0157.resize(mednum_fuel);
  old_np155.resize(mednum_fuel);
  old_rp155.resize(mednum_fuel);
  old_rp157.resize(mednum_fuel);
  old_rc155.resize(mednum_fuel);
  old_rc157.resize(mednum_fuel);
  old_n0155_2.resize(mednum_fuel);
  old_r0155_2.resize(mednum_fuel);
  old_r0157_2.resize(mednum_fuel);
  old_x0155.resize(mednum_fuel);
  old_x0157.resize(mednum_fuel);

  if(corrector_calc){
    rr_gd5.resize(mednum_fuel);
    rr_gd7.resize(mednum_fuel);
    np_gd5.resize(mednum_fuel);
    np_gd7.resize(mednum_fuel);
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();

  // ... Preparing medium class instances for ADTS model
  //
  // 
  // vector<Medium> med_adts(mednum_fuel);

  // +++ Forward burn-up calculation +++
  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  // +++ Eigenvalue calculation
  cout<<"#\n# Predictor calculation in Original-Model ...\n";
  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);  

  real abs_frac[burn_step+1];
  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    //lat.MedClear();
    for(int i=0;i<mednum_fuel;i++){
      // +++ Self-shielding calculation
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library      

  // med_adts[i]=med[0]; // fuga addition

      /*
      // ... neutron flux energy spectrum calculated at the last step is used
      //     to calculate macroscopic fission spectrum
      if(st!=0){
	GroupData1D aveflx;
	aveflx.set_zero();
        for(int ii=0;ii<totm;ii++){
          int medid=region_medium[ii];
          if(medid==i)aveflx=aveflx+volflx_mesh[st-1][ii];
        };
	med[0].GetFlux().copy(aveflx);
      };
      */

      med[0].CalMacroFromMicro();

      /*
      cout<<"# Medium :: "<<i<<"\n";
      med[0].ShowMacroXS1D();
      */

      /*
      if(i==40){
	for(int k=0;k<nucn;k++){
	  cout<<k<<" "<<med[0].GetNuclideInTurn(k).GetMatnum()<<" "<<fwd_nuc[st][0][i].get_dat(k)<<"\n";
	};
      };
      */

      // +++ Macro&Micro data storing
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      if(st==0){
        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();
      }else{
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
      };
    };
    if(st==0){
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };
    }else{
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };
    };

    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);

    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    lat.Print();
    keff[st]=lat.CalIgen();

    //cout<<"# Total inner iteration : "<<lat.GetTotalInnerIterationSum()<<"\n"; exit(0);

    // ------------------------------

    real sum_tot=0.;
    real sum_fmed=0.;
    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);

      real tmp=flx*lat.GetMedium(i).GetMacxs().GetData1d(siga);
      sum_tot+=tmp;
      if(i<mednum_fuel)sum_fmed+=tmp;

      flx_med[i]=flx*(1./vol_med[i]);
      // +++ One-group cross section storing (at predictor step)
      if(i<mednum_fuel){
        for(int j=0;j<nucn;j++){
          if(nuclide_info[j]!=0){
            if(nuclide_info[j]==1){
   	      xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	    };
	    xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	    xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
          };
	};
      };
    };


    abs_frac[st]=sum_fmed/sum_tot; // absorption rate fraction in fuel media
    cout<<"# Absorption fraction : "<<abs_frac[st]<<"\n";
  
    if(st!=burn_step){

  // fuga addition
  // ... with group collapsing from 172-group to fewer-group ...
  // +++ Instancee generation of MECSystem to solve neutron transport equation for simplified system

  cout<<"#\n# Predictor calculation in Simplified-Model ...\n";

  vector< vector<GroupData1D> > mic_sigf_sm(mednum_fuel); // [medid][nucn](group)
  vector< vector<GroupData1D> > mic_sigc_sm(mednum_fuel);
  vector< vector<GroupData1D> > mic_sign2n_sm(mednum_fuel);
  vector< vector<real> > xsf_1g_p_sm(mednum_fuel); // [medid][nucn]
  vector< vector<real> > xsc_1g_p_sm(mednum_fuel);
  vector< vector<real> > xsn2n_1g_p_sm(mednum_fuel);

  MECSystem lat2(ngrp,mednum);
  lat2.PutTrajectorySet(&sys_f2);

	  for(int i=0;i<mednum_fuel;i++){
        // (Self-shielding calculation) ---------------------------------------
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
        med[0].CalMacroFromMicro();

        Medium bmed=med[0].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
        // ------------------------------------------------------------------------
	
  
	// (Use of the original burnup step results) -------------------------------
  //       int nucn=med_adts[i].GetNucnum();
  //       for(int ii=0;ii<nucn;ii++){
  //          med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
  //       };
	// med_adts[i].CalMacroFromMicro();
	// Medium bmed=med_adts[i].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
	// -------------------------------------------------------------------------
  

        // (Macro&Micro data storing)
  mic_sigf_sm[i].resize(nucn);
  mic_sigc_sm[i].resize(nucn);
  mic_sign2n_sm[i].resize(nucn);
  xsf_1g_p_sm[i].resize(nucn,0.);
  xsc_1g_p_sm[i].resize(nucn,0.);
  xsn2n_1g_p_sm[i].resize(nucn,0.);
        for(int k=0;k<nucn;k++){
	  
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat2.AddMedium(bmed); 
        lat2.GetMedium(i).NuclideClear();
      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel*jj],true); // [true] is for microscopic XS
        lat2.AddMedium(bmed);
      };

      // +++ eigenvalue calculation
      lat2.PutRegMed(region_medium);
      lat2.PutGeneralOption(opt);
      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      lat2.Print();

      real keff_dummy=lat2.CalIgen();

      vector<real> medium_wise_flux_change(mednum);
      
      for(int i=0;i<mednum;i++){
        GroupData1D flx_sm=lat2.GetIntegratedFlux(i);

	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
	medium_wise_flux_change[i]=flx_med[i].get_sum()/(flx_sm.get_sum()/vol_med[i]);

        // flx_med_sm[i]=flx_sm*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g_p_sm[i][j]=mic_sigf_sm[i][j].Cond(flx_sm);
  	      };
	      xsc_1g_p_sm[i][j]=mic_sigc_sm[i][j].Cond(flx_sm);
	      xsn2n_1g_p_sm[i][j]=mic_sign2n_sm[i][j].Cond(flx_sm);
            };
  	  };
        };
      };
  // ------------------------------

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //output 
  // fout_1g<<"\n#\n# p : "<<st<<endl;
  // fout_p<<"\n#\n# "<<st<<endl;
  // for(int i=0;i<mednum_fuel;i++){
  //   for(int j=0;j<nucn;j++){
  //     int nulid=med[0].GetNuclideInTurn(j).GetMatnum();

  //     if(nuclide_info[j]!=0){
  //       if(nuclide_info[j]==1){
  //   fout_1g<<xsf_1g[st][i][j]<<" "<<xsf_1g_p_sm[i][j]<<" ";
  //   // fout_p<<mic_sigf[bstmp][i][j]<<" "<<mic_sigf_sm[i][j]<<" ";
  //       };
  //   fout_1g<<xsc_1g[st][i][j]<<" "<<xsc_1g_p_sm[i][j]<<" ";
  //   fout_1g<<xsn2n_1g[st][i][j]<<" "<<xsn2n_1g_p_sm[i][j]<<" \n";

  //   // fout_p<<mic_sigc[bstmp][i][j]<<" "<<mic_sigc_sm[i][j]<<" ";
  //   // fout_p<<mic_sign2n[bstmp][i][j]<<" "<<mic_sign2n_sm[i][j]<<" \n";
  //     };
  //   };
  //   fout_tf<<flx_med[i].get_sum()<<medium_wise_flux_change[i]<<" ";
  //   cout<<flx_med[i].get_sum()<<" "<<medium_wise_flux_change[i]<<endl;
  // };
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };

      // +++ Burnup calculation
      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
          //power_org=CalculationPinPower(bu,st,j,med_target,sumflx*vol_med[med_target]);
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        accumulated_day+=burn_span;
	//total_flux[st][j].resize(mednum_fuel);
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false); // if [adjoint] is true, multistep calculation is done.
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //   CORRECTOR CALCULATION
        //  to solve neutron transport equation for Simplified system
      if(corrector_calc){

	cout<<"#\n# Corrector calculation in Simplified-Model...\n";

  // (adts) 
   // fuga addition
  vector< vector<GroupData1D> > mic_sigf_sm(mednum_fuel); // [medid][nucn](gourp)
  vector< vector<GroupData1D> > mic_sigc_sm(mednum_fuel);
  vector< vector<GroupData1D> > mic_sign2n_sm(mednum_fuel); 
  vector< vector<real> > xsf_1g_c_sm(mednum_fuel); // [medid][nucn]
  vector< vector<real> > xsc_1g_c_sm(mednum_fuel);
  vector< vector<real> > xsn2n_1g_c_sm(mednum_fuel);
  //-------------

	// (OWPC)	
	for(int i=0;i<mednum_fuel;i++){
	  real tmp=total_flux[st][0][i];
  	  rr_gd5[i]=tmp*xsc_1g[st][i][pos_gd155];  // Reaction rate at BOC (Rp)
  	  rr_gd7[i]=tmp*xsc_1g[st][i][pos_gd157];
	};

	// ++++++++++++++++++++++++

	// (Predictor calculation results are stored in the array of [XXX_p].)
        swap(total_flux_p[st], total_flux[st]);
        swap(xsc_1g_p[st], xsc_1g[st]);
        swap(xsn2n_1g_p[st], xsn2n_1g[st]);
        swap(xsf_1g_p[st], xsf_1g[st]);


      for(int i=0;i<mednum_fuel;i++){
        // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
        PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
  	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library	 
      // med_adts[i]=med[0] // ... adts
        med[0].CalMacroFromMicro();
        //opc.CalFissionSpectrumMatrix(med[0],xslib);

        Medium bmed=med[0].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS // fuga addition

        // +++ Macro&Micro data storing
  mic_sigf_sm[i].resize(nucn);
  mic_sigc_sm[i].resize(nucn);
  mic_sign2n_sm[i].resize(nucn);
  xsf_1g_c_sm[i].resize(nucn,0.);
  xsc_1g_c_sm[i].resize(nucn,0.);
  xsn2n_1g_c_sm[i].resize(nucn,0.);
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf)); 
            };
            mic_sigc_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc)); 
            mic_sign2n_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n)); 
          };
        };

	lat2.GetMedium(i).GetMacxs()=bmed.GetMacxs(); // fuga addition
      };

      for(int jj=0;jj<mednum_nonfuel;jj++){
        Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel*jj],true); // [true] is for microscopic XS
        lat2.GetMedium(mednum_fuel+jj).GetMacxs()=bmed.GetMacxs();
      };

      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      real keff_corr=lat2.CalIgen();
      
      // ------------------------------

      for(int i=0;i<mednum;i++){
        GroupData1D flx_sm=lat2.GetIntegratedFlux(i);
        flx_med_c[i]=flx_sm*(1./vol_med[i]);
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1)xsf_1g_c_sm[i][j]=mic_sigf_sm[i][j].Cond(flx_sm);
	      xsc_1g_c_sm[i][j]=mic_sigc_sm[i][j].Cond(flx_sm);
	      xsn2n_1g_c_sm[i][j]=mic_sign2n_sm[i][j].Cond(flx_sm);
	    };
          };
	};
      };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // fout_c<<"\n#\n# "<<st<<endl;
    // fout_1g<<"\n#\n# c : "<<st<<endl;
    // for(int i=0;i<mednum_fuel;i++){
    //   for(int j=0;j<nucn;j++){
    //     if(nuclide_info[j]!=0){
    //       if(nuclide_info[j]==1){
    //   fout_1g<<xsf_1g_c_sm[i][j]<<" ";
    //   // fout_c<<mic_sigf_sm[i][j]<<" ";
    //       };
    //   fout_1g<<xsc_1g_c_sm[i][j]<<" "<<xsn2n_1g_c_sm[i][j]<<" \n";
    //   // fout_c<<mic_sigc_sm[i][j]<<mic_sign2n_sm[i][j]<<" \n";
    //     };
    //   };
    //   fout_tf<<flx_med_c[i].get_sum()<<" \n";
    //   cout<<flx_med_c[i].get_sum()<<endl;
    // };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //++++++++++++++++++++++++++++++++++++++++++++++
    // ... Correction in adts
    cout<<"# Correction ...\n";
      
      // ... one group cross section  
      for(int i=0;i<mednum_fuel;i++){
	      for(int k=0;k<nucn;k++){
  
  	    int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
  #if 0
            // ... Copy of the ADTS step data
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
  	        xsf_1g[st][i][k]=xsf_1g_c_sm[i][k];
	      };
	      xsc_1g[st][i][k]=xsc_1g_c_sm[i][k];
	      if(xsn2n_1g_p_sm[i][k]>0.){
  	        xsn2n_1g[st][i][k]=xsn2n_1g_c_sm[i][k];
	      };
	    };
  #endif	  

  #if 1
            // ... Correcting all nuclides in all burnup regions
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
            // cout<<"# fission\n";
            // cout<<xsf_1g[st][i][k]<<" "<<xsf_1g_p_sm[i][k]<<" "<<xsf_1g_c_sm[i][k]<<endl;
  	        if(xsf_1g_p_sm[i][k]>0.)xsf_1g[st][i][k]=xsf_1g[st][i][k]*xsf_1g_c_sm[i][k]/xsf_1g_p_sm[i][k]; //for neutron fission cross section
            // if(isfinite(xsf_1g[st][i][k]))cout<<"# error in fission.\n";exit(0);
        };
        // cout<<"# capture\n"<<xsc_1g[st][i][k]<<" "<<xsc_1g_p_sm[i][k]<<" "<<xsf_1g_c_sm[i][k]<<endl;
	      if(xsc_1g_p_sm[i][k]>0.)xsc_1g[st][i][k]=xsc_1g[st][i][k]*xsc_1g_c_sm[i][k]/xsc_1g_p_sm[i][k]; // for neutron capture cross section
        // cout<<xsc_1g[st][i][k]<<endl;
        // if(isfinite(xsc_1g[st][i][k]))cout<<"# error in capture.\n";exit(0);
        // if(xsn2n_1g[st][i][k]>0.){
            // cout<<"# n2n\n"<<xsn2n_1g[st][i][k]<<" "<<xsn2n_1g_p_sm[i][k]<<" "<<xsn2n_1g_c_sm[i][k]<<endl;
  	        if(xsn2n_1g_p_sm[i][k]>0.)xsn2n_1g[st][i][k]=xsn2n_1g[st][i][k]*xsn2n_1g_c_sm[i][k]/xsn2n_1g_p_sm[i][k]; // for n_2n reaction cross section
            // if(isfinite(xsn2n_1g[st][i][k]))cout<<"# error in n2n.\n";exit(0);
        // };
	    };
  #endif	  

  #if 0	  
	    if(i<8){
              // ... Correcting all nuclides in the central Gd-bearing pin
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
              xsf_1g[st][i][k]=xsf_1g[st[i][k]*xsf_1g_c_sm[i][k]/xsf_1g_p_sm[i][k]; //for neutron fission cross section
          };
          xsc_1g[st][i][k]=xsc_1g[st][i][k]*xsc_1g_c_sm[i][k]/xsc_1g_p_sm[i][k];//for neutron capture cross section
          if(xsn2n_1g[st-1][i][k]>0.){
              xsn2n_1g[st][i][k]=xsn2n_1g[st][i][k]*xsn2n_1g_c_sm[i][k]/xsn2n_1g_p_sm[i][k]; // for n_2n reaction cross section
          };
        };
	    };
  #endif	  

  #if 0	  
	    if(i<8){
	      // ... Correcting Gd-155 and -157 capture cross sections in the central Gd-bearing pin
	      if(nucid==641550||nucid==641570){
	        // if(nucid==641550) cout<<"Gd-155 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        // if(nucid==641570) cout<<"Gd-157 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        // xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
	        xsc_1g[st][i][k]=xsc_1g[st][i][k]*xsc_1g_c_sm[i][k]/xsc_1g_p_sm[i][k];//for neutron capture cross section
	        // xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
	      };
	    };
  #endif	  

	  };
        };

	// ... Neutron flux correction in ADTS
	for(int i=0;i<mednum_fuel;i++){ // All the medium included in the system
	//for(int i=0;i<8;i++){ // Only for Gd-bearing rod
	  flx_med_c[i]=flx_med_c[i]*medium_wise_flux_change[i];
	};

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real acburn_pre=accumulated_burn-acburn[st];
      vector<real> acburn_pre_per_medium(mednum_fuel);
      for(int j=0;j<mednum_fuel;j++){
	acburn_pre_per_medium[j]=accumulated_burn_per_medium[j]-acburn_per_medium[j][st];
      };
      // accumulated burnup calculated by the predictor step

      
      // 
      // .....  One-group cross section correction for Rc in A-OWPC .....
      //

      
      if(owpc_corr){
	
      for(int i=0;i<mednum_fuel;i++){

	real n0=fwd_nuc[st][0][i].get_dat(pos_gd157);   // initial 
        real np=fwd_nuc[st+1][0][i].get_dat(pos_gd157); // predictor results      
        if((n0-np)/np>1e-2){
          // This if branch is added in 2021/12/4 since negative cross section is detected in the Gd cross section.
	  // This should be consistent in the following process in the actual correction to the final ND.
	
	real xsc0_5=xsc_1g_p[st][i][pos_gd155]; //rp
	real xsc0_7=xsc_1g_p[st][i][pos_gd157]; //rp
	real xscc_5=xsc_1g[st][i][pos_gd155];//Rcに相当するXc
	real xscc_7=xsc_1g[st][i][pos_gd157];
	if(st>0){
	  real n0_5=fwd_nuc[st][0][i].get_dat(pos_gd155);
	  real np_5=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	  real t_5=old_xsc5[i]+(xsc0_5-old_xsc5[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);//前のステップでの正しいXcを計算
	  real t_7=old_xsc7[i]+(xsc0_7-old_xsc7[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
	  if(st>1){
	    real a=n0_5, b=old_n0155[i], c=old_n0155_2[i];
	    real aa_5=xsc0_5, bb_5=old_xsc5[i], cc_5=old_xsc5_2[i];
	    real aa_7=xsc0_7, bb_7=old_xsc7[i], cc_7=old_xsc7_2[i];
	    real X_5=old_np155[i];
	    t_5=aa_5*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_5*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_5*(X_5-a)*(X_5-b)/(c-a)/(c-b);//前のステップでの正しいXcを計算
	    t_7=aa_7*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_7*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_7*(X_5-a)*(X_5-b)/(c-a)/(c-b);
	  };

	  // correction with considering time step length x power density (= burnup) 
	  real tmp1=power_density_list[st-1]*burn_time[st-1];
	  real tmp2=power_density_list[st]*burn_time[st];
	  real factor=tmp2/tmp1;
	  real xsc155=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)*factor;
	  real xsc157=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)*factor;
	  if(factor<10.){ // if the burnup length is too different, the correction is NOT applied.
	    if(xsc155>0.)xsc_1g[st][i][pos_gd155]=xsc155;
	    if(xsc157>0.)xsc_1g[st][i][pos_gd157]=xsc157;
	  };
	  
	  old_xsc5_2[i]=old_xsc5[i];
	  old_xsc7_2[i]=old_xsc7[i];
	};
	old_xsc5[i]=xsc0_5;
	old_xsc7[i]=xsc0_7;
	old_xscc_5[i]=xscc_5;
	old_xscc_7[i]=xscc_7;
      };

      }; // end of [if((n0-np)/np>1e-2)]
      }; // end of [if(owpc_corr)]
      // ............................................................................
      
      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med_c[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med_c[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          //accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        //accumulated_day+=burn_span;
        //delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
	  //CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.

	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false); // if [adjoint] is true, multistep calculation is done.
	};
	
      }; // end of sub-step
      
      // (OWPC)
      real dt=burn_time[st]*60*60*24;
      real time_mesh=10000;
      real dt_r=dt/time_mesh;
      real omega_155, omega_157; // various correlation conditions

      for(int i=0;i<mednum_fuel;i++){

	  // OWPC treatment for Gd-155 and -157 with various correlation conditions
	if(owpc_corr){

	  if(st>0){ // ... quadratic model
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	    real x0_155=xsc_1g_p[st][i][pos_gd155];//rpに対応する断面積
	    real x0_157=xsc_1g_p[st][i][pos_gd157];//rpに対応する断面積
	    real xp_155=xsc_1g[st][i][pos_gd155]; //rcに対応する断面積
	    real xp_157=xsc_1g[st][i][pos_gd157]; //rcに対応する断面積
	    real n0_r_155;
	    real n0_r_157;
	    
	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt);//Not consider 154,156
	    	    
	    //quadratic model---------------------------
	    real a_5,b_5,c_5,a_7,b_7,c_7;
	    real aa_5,bb_5,cc_5,aa_7,bb_7,cc_7;
	    a_5=n0_155;
	    b_5=np_155;
	    c_5=old_n0155[i];
	    a_7=n0_157;
	    b_7=np_157;
	    c_7=old_n0157[i];
	    real X_5=np_p_155;
	    real X_7=np_p_157;

	    //-------------------------reaction rate base-------------------------------------
	    aa_5=r0_155;
	    bb_5=rp_155;
	    //cc_5=old_r0155[i];
	    cc_5=old_r0155[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba    
	    
	    aa_7=r0_157;
	    bb_7=rp_157;
	    //cc_7=old_r0157[i];
	    cc_7=old_r0157[i]*power_density_list[st]/power_density_list[st-1];	 // Correction by chiba

	    real rc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    real rc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    //----------------------------------------------------------------------------------------
	    
	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	    real xp_r_155=x0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	    real xp_r_157=x0_157;

	    //------------------------------------start toy problem calculation----------------------------
	    for(int k=0;k<time_mesh;k++){


	      //-----reaction rate base-------
	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      X_5=np_r_155;
	      X_7=np_r_157;
	      
	      real rc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	      real rc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r); //Not consider 154,156
	    
	      rp_r_155=rc_r_155;
	      rp_r_157=rc_r_157;
	      //----------------------------------------

	
	    };
	    //-----------------------------------------end toy problem calculation------------------------------------------
	    
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

	    //-------reaction rate base----------
	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	    //------------------------------------------

	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;
	    
	    old_n0155_2[i]=old_n0155[i];
	    old_r0155_2[i]=old_r0155[i];
	    old_r0157_2[i]=old_r0157[i];
	    
	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;
	    old_x0155[i]=x0_155;
	    old_x0157[i]=x0_157;
	    
	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;
	      	    
	  }else{ // ... linear model for [st]==0
	  
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	
	    real n0_r_155;
	    real n0_r_157;

	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt); //Not consider 154,156
	  
	    // (Correlation to Gd-155 number density)
	    real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	    real alpha_157=(rp_157-r0_157)/(np_155-n0_155);    
	    
	    // (Correlation to Gd-155 number density)	
	    real rc_155=r0_155+(np_p_155-n0_155)*alpha_155; // Rc in toy-problem
	    real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	   
	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	  
	    for(int k=0;k<time_mesh;k++){

	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      // (Correlation to Gd-155 number density)
	      real rc_r_155=r0_155+(np_r_155-n0_155)*alpha_155;
	      real rc_r_157=r0_157+(np_r_155-n0_155)*alpha_157;
	
	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r);//Not consider 154,156

	      rp_r_155=rc_r_155;
	      rp_r_157=rc_r_157;

	    };
	        
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);

	    //cout<<"#   Omega_155 / Omega_157 : "<<omega_155<<" "<<omega_157<<"\n";
	    
	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;

	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;

	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;

	  }; //

	}; // end of [if(owpc_corr)]
	  
	for(int j=0;j<nucn;j++){
	
	  real np=fwd_nuc[st+1][0][i].get_dat(j);      // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
	  real n_next=exp((log(np)+log(nc))*0.5);
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();

	  if(nucid==641550||nucid==641570){

  	    real n0=fwd_nuc[st][0][i].get_dat(j);
            if((n0-np)/np>1e-2){

	      if(owpc_corr){
		if(nucid==641550){
		  n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
		}else if(nucid==641570){
		  n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
    	        };
	      }else{
		/*
		// -- OWPC w/o (N,R) correlation ----------------------- 
		real r0=rr_gd5[i]; //Rp
		if(nucid==641570)r0=rr_gd7[i];

		real rp=xsc_1g[st][i][j]*total_flux[st][0][i]; //Rc
		real alpha=(rp-r0)/(np-n0);
		  
		real np_p=n0*exp(-r0*1e-24*dt);   // Nc in toy-problem
		real rc=alpha*np_p+(r0-n0*alpha); // Rc in toy-problem
	    
		real n0_r=n0;
		real rp_r=r0;
		for(int k=0;k<time_mesh;k++){
		  real np_r=n0_r*exp(-rp_r*1e-24*dt_r);
		  real rc_r=alpha*np_r+(r0-n0*alpha);
		  real r_r=(rp_r+rc_r)/2;
		  n0_r=n0_r*exp(-r_r*1e-24*dt_r);
		  rp_r=rc_r;
		};
		
		real R_reference=(log(n0)-log(n0_r))/dt*1e24;
		real omega=(R_reference-r0)/(rc-r0);
		  
		// -- PPC -------------------------------
		real R_p=-log(np/n0)/dt;
		real R_c=-log(nc/n0)/dt;
		real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
		R_reference=(R_p+R_c_c)/2;
		omega=(R_reference-R_p)/(R_c-R_p);
		// --------------------------------

		n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応
		*/

 	        // ... conventional PC ...
                //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
                //n_next=(np+nc)*0.5;  // Conventional PC (in linear)
                //n_next=nc;          
              };

	    };

	  };


          // if(adjoint)        n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
          if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC

	  // (For OWPC)
    	  if(nucid==641550)np_gd5[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
  	  if(nucid==641570)np_gd7[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd157);

	  
          fwd_nuc[st+1][0][i].put_data(j,n_next);

	}; // end of nuclides iteration

      }; // end of medium iteration

      // Adjustment of accumulated burn
      for(int i=0;i<mednum_fuel;i++){
        real acburn_cor=accumulated_burn_per_medium[i]-acburn_per_medium[i][st]-acburn_pre_per_medium[i];
        accumulated_burn_per_medium[i]=acburn_per_medium[i][st]+(acburn_pre_per_medium[i]+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };

      if(input_flux_level){
        real acburn_cor=accumulated_burn-acburn[st]-acburn_pre;
        accumulated_burn=acburn[st]+(acburn_pre+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };
      cout<<"#     ... terminated.\n";
      }; // END OF CORRECTOR CALCULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    }; // end part of [If (st!=burn_step) ]
 
  }; // loop-end of burnup step

};

void MulticellBurner::ForwardCalculationADTS_SamSys(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  // fuga addition ------------------------
  
  // enter
    bool nuclide_all=false;
    vector<int> id,newline;
    if(nuclide_all){
      newline={0,1,2,5,11,13,18,22,29,34,40,41,46,47,54,59,69,74,77,80,86,91,102,107,118,125,133,135,138,143,145,150,154,159};
      for(int i=0;i<nucn;i++){
        int matid=med[0].GetNuclideInTurn(i).GetMatnum();
        id.push_back(matid);
      };
    }else{
      newline={0};
      id={641550,641570};
    };
    int sz=id.size();
    vector<string> rct={"ALL"}; // "f", "c", "n2n" or "ALL"
    int react=0;
    bool ncl_f,ncl_c,ncl_n2n;
    ncl_f=false;
    ncl_n2n=false;
    ncl_c=false;
  // -----
   
    // output data file
    ofstream fout_1g("xs1g.txt"); // One-gourp cross section with OM and SM 
    ofstream fout_tf("total_flux.txt"); // Total neutron flux with OM and SM
    ofstream fout_sm("micXS_sm.txt"); // Enerrgy distribution of cross section in Simplified-Model
    ofstream fout_om("micXS_om.txt"); // Enerrgy distribution of cross section in Original-Model 
    // ofstream fout_fai("fuelNflux.txt"); // Neuton energy spectrum with OM and SM 

    for(int i=0;i<rct.size();i++){
      if("f"==rct[i]){
        ncl_f=true;
        react++;
      }else if("c"==rct[i]){
        ncl_c=true;
        react++;
      }else if("n2n"==rct[i]){
        ncl_n2n=true;
        react++;
      }else if("ALL"==rct[i]){
        ncl_f=true;
        ncl_c=true;
        ncl_n2n=true;
        react+=3;
      };
    };
    vector<string> prt_nuc_nam;
    for(int i=0;i<sz;i++){
      for(int j=0;j<nucn;j++){
        int metid=med[0].GetNuclideInTurn(j).GetMatnum();
        if(metid==id[i]){prt_nuc_nam.push_back(midt.Name(metid));};
      };
    };

    
    fout_1g<<"# One-group cross section [barn]\n# ";
    for(int i=0;i<sz;i++){fout_1g<<prt_nuc_nam[i]<<" ";};
    fout_1g<<"\n# ";
    if(ncl_f){fout_1g<<"fission ";};
    if(ncl_c){fout_1g<<"capture ";};
    if(ncl_n2n){fout_1g<<"n2n ";};
    fout_1g<<endl;    
    fout_1g<<"#----------------------------------------------------------------------------\n";
    fout_1g<<"# Number of Medium   : "<<mednum_fuel<<endl;
    fout_1g<<"# Number of nuclides : "<<sz<<endl;
    int num_newline=0;
    fout_1g<<"# ";
    for(int i=0;i<sz;i++){
      fout_1g<<prt_nuc_nam[i]<<":"<<IntToString(i)<<" ";
      if(newline[num_newline]==i){fout_1g<<"\n# ";num_newline++;};
    };
    fout_1g<<"{reaction}           : ";
    if(ncl_f){fout_1g<<"(fission)";};
    if(ncl_c){fout_1g<<"(capture)";};
    if(ncl_n2n){fout_1g<<"(n2n)";};
    fout_1g<<endl;
    fout_1g<<"# Model              : OM, SM\n#\n";
    fout_1g<<"# column = Burnup(1) + Medium("<<IntToString(mednum_fuel)<<")*{flx(1) + nuclide("<<IntToString(sz)<<")*reactio("<<IntToString(react)<<")*Model(2)}\n";
    fout_1g<<"#----------------------------------------------------------------------------\n";
    fout_1g<<"# (Burnup)   ";
    int num=13*sz*react*2;
    for(int i=0;i<mednum_fuel;i++){
      fout_1g<<" (Med "+IntToString(i)+")";   
      for(int j=0;j<num;j++){fout_1g<<" ";};
    };
    fout_1g<<"# [GWd/t]\n";


    fout_tf<<"# Total neutron flux\n#\n";
    fout_tf<<"# ---------------------------------------------------------------------\n";
    fout_tf<<"# Number of Medium : "<<mednum_fuel<<endl;
    fout_tf<<"# Model            : OM, SM\n";
    fout_tf<<"# ---------------------------------------------------------------------\n";
    fout_tf<<"# (Burnup)   ";
    for(int i=0;i<mednum_fuel;i++){
      string numid=IntToString(i);
      fout_tf<<" (Med "+numid+")";
      for(int j=0;j<13*2-7-numid.size();j++){fout_tf<<" ";};
    };
    fout_tf<<endl;
    fout_tf<<"# [GWd/t]    ";
    for(int i=0;i<mednum_fuel;i++){fout_tf<<" (ori)       (sim)        ";};
    fout_tf<<endl;

    
    fout_om<<"# Energy distribution of cross section in Original-Model\n# ";
    for(int i=0;i<sz;i++){fout_om<<prt_nuc_nam[i]<<" ";};
    fout_om<<"\n# ";
    if(ncl_f){fout_om<<"fission ";};
    if(ncl_c){fout_om<<"capture ";};
    if(ncl_n2n){fout_om<<"n2n ";};
    fout_om<<endl;
    fout_om<<"#----------------------------------------------------------------------------\n";
    fout_om<<"# Number of Medium   : "<<mednum_fuel<<endl;
    fout_om<<"# Number of nuclides : "<<sz<<endl;
    num_newline=0;
    fout_om<<"# ";
    for(int i=0;i<sz;i++){
      fout_om<<prt_nuc_nam[i]<<":"<<IntToString(i)<<" ";
      if(newline[num_newline]==i){fout_om<<"\n# ";num_newline++;};
    };
    fout_om<<"{reaction}           : ";
    if(ncl_f){fout_om<<"(fission)";};
    if(ncl_c){fout_om<<"(capture)";};
    if(ncl_n2n){fout_om<<"(n2n)";};
    fout_om<<endl;
    fout_om<<"# Model              : OM, SM\n#\n";
    fout_om<<"# column = Energy(1) + Medium("<<IntToString(mednum_fuel)<<")*{flx(1) + nuclide("<<IntToString(sz)<<")*reactio("<<IntToString(react)<<")*Model(2)}\n";
    fout_om<<"#----------------------------------------------------------------------------\n";

    
    fout_sm<<"# Energy distribution of cross section in Simplified-Model\n# ";
    for(int i=0;i<sz;i++){fout_sm<<prt_nuc_nam[i]<<" ";};
    fout_sm<<"\n# ";
    if(ncl_f){fout_sm<<"fission ";};
    if(ncl_c){fout_sm<<"capture ";};
    if(ncl_n2n){fout_sm<<"n2n ";};
    fout_sm<<endl;
    fout_sm<<"#----------------------------------------------------------------------------\n";
    fout_sm<<"# Number of Medium   : "<<mednum_fuel<<endl;
    fout_sm<<"# Number of nuclides : "<<sz<<endl;
    num_newline=0;
    fout_sm<<"# ";
    for(int i=0;i<sz;i++){
      fout_sm<<prt_nuc_nam[i]<<":"<<IntToString(i)<<" ";
      if(newline[num_newline]==i){fout_sm<<"\n# ";num_newline++;};
    };
    fout_sm<<"{reaction}           : ";
    if(ncl_f){fout_sm<<"(fission)";};
    if(ncl_c){fout_sm<<"(capture)";};
    if(ncl_n2n){fout_sm<<"(n2n)";};
    fout_sm<<endl;
    fout_sm<<"# Model              : OM, SM\n#\n";
    fout_sm<<"# column = Energy(1) + Medium("<<IntToString(mednum_fuel)<<")*{flx(1) + nuclide("<<IntToString(sz)<<")*reactio("<<IntToString(react)<<")*Model(2)}\n";
    fout_sm<<"#----------------------------------------------------------------------------\n";


    // fout_fai<<"# Neutron energy spectrum\n#\n" 
    // fout_fai<<"# ---------------------------------------------------------------------\n";
    // fout_fai<<"# Number of Medium : "<<mednum_fuel<<endl;
    // fout_fai<<"# Model            : OM, SM\n";
    // fout_fai<<"# ---------------------------------------------------------------------\n";
  //--------------------------------------


  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };
  
  // fuga addition
  vector< vector< vector<real> > > xsc_1g_sm; // [step][medid][nucn]
  vector< vector< vector<real> > > xsn2n_1g_sm;
  vector< vector< vector<real> > > xsf_1g_sm;
  // vector< vector< vector<GroupData1D> > > mic_sigf_sm; //[step][medid][nucn](group)
  // vector< vector< vector<GroupData1D> > > mic_sigc_sm;
  // vector< vector< vector<GroupData1D> > > mic_sign2n_sm; 
  int mednum_fuel_sm=mednum_fuel;
  // ------------ 
  
  
  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  // fuga addition
  xsc_1g_sm.resize(burn_step+1);
  xsn2n_1g_sm.resize(burn_step+1);
  xsf_1g_sm.resize(burn_step+1);
  // ------------
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    // fuga adition
    xsc_1g_sm[i].resize(mednum_fuel_sm);
    xsn2n_1g_sm[i].resize(mednum_fuel_sm);
    xsf_1g_sm[i].resize(mednum_fuel_sm);
    // ------------
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
       // fuga addition
      xsc_1g_sm[i][j].resize(nucn,0.);
      xsn2n_1g_sm[i][j].resize(nucn,0.);
      xsf_1g_sm[i][j].resize(nucn,0.);
      // -------------
   };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  vector<GroupData1D> flx_med_sm;
  flx_med.resize(mednum);
  flx_med_sm.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
    flx_med_sm[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  // // fuga addition
  // mic_sigf_sm.resize(tmp); 
  // mic_sigc_sm.resize(tmp); 
  // mic_sign2n_sm.resize(tmp);
  // -----------
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    // // fuga addiiton
    // mic_sigf_sm[i].resize(mednum_fuel_sm);
    // mic_sigc_sm[i].resize(mednum_fuel_sm);
    // mic_sign2n_sm[i].resize(mednum_fuel_sm);
    // ------------
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
            // mic_sigf_sm[i][j][k].put_imax(group); // fuga addition
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
          // mic_sigc_sm[i][j][k].put_imax(group);

        };
      };
    };
  };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();


  // ... Preparing medium class instances for ADTS model
  //
  //
  vector<Medium> med_adts(mednum_fuel);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  int st_10;
  for(int st=0;st<burn_step+1;st++){
    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    
    if(st%sm_step==0){

      st_10=st;      

      // +++ Instance generation of MECSystem to solve neutron transport equation
      MECSystem lat(group,mednum);
      lat.PutTrajectorySet(&sys_f);  

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel;i++){

        // (Self-shielding calculation)
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048); // ESAB from MVP library
	med_adts[i]=med[0];  // ... adts
        med[0].CalMacroFromMicro();

        // (Macro&Micro data storing)
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();

      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };

      // +++ eigenvalue calculation
      lat.PutRegMed(region_medium);
      lat.PutGeneralOption(opt);
      lat.PutThermalIteration(3);
      lat.PutPL(0);
      lat.NoCMRAcceleration();
      lat.Print();
      
      keff[st]=lat.CalIgen();

      for(int i=0;i<mednum;i++){
        GroupData1D flx=lat.GetIntegratedFlux(i);
        flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	      };
	      xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
            };
	  };
        };
      };

    }else{

      // One-group cross section at the last burnup step is assumed
      for(int i=0;i<mednum;i++){
        if(i<mednum_fuel){
	  for(int j=0;j<nucn;j++){
	    if(nuclide_info[j]!=0){
	      if(nuclide_info[j]==1){
	        xsf_1g[st][i][j]=xsf_1g[st_10][i][j];
	      };
	      xsc_1g[st][i][j]=xsc_1g[st_10][i][j];
	      xsn2n_1g[st][i][j]=xsn2n_1g[st_10][i][j];
	    };
	  };
	};
      };
      
    };
      
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){


  #if 0
      // ... No group collapsing ...
      // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
      MECSystem lat2(group,mednum);
      lat2.PutTrajectorySet(&sys_f2);        

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel;i++){

	/*
        // (Self-shielding calculation) -----------------------------------------------------------
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
        med[0].CalMacroFromMicro();
	// ----------------------------------------------------------------
	*/

	// (Use of the original burnup step results) -------------------------------
        int nucn=med_adts[i].GetNucnum();
        for(int ii=0;ii<nucn;ii++){
           med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
        };
	med_adts[i].CalMacroFromMicro();
	med[0]=med_adts[i];
	// -------------------------------------------------------------------------
	
	
        // (Macro&Micro data storing)
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_sm[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_sm[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_sm[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat2.AddMedium(med[0]); 
        lat2.GetMedium(i).NuclideClear();

      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat2.AddMedium(med[1+jj]);
      };

      // +++ eigenvalue calculation
      lat2.PutRegMed(region_medium);
      lat2.PutGeneralOption(opt);
      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      lat2.Print();

      real keff_dummy=lat2.CalIgen();

      vector<real> medium_wise_flux_change(mednum);
      
      for(int i=0;i<mednum;i++){
        GroupData1D flx_sm=lat2.GetIntegratedFlux(i);

	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
	medium_wise_flux_change[i]=(flx_sm.get_sum()/vol_med[i])/flx_med_sm[bstmp][i].get_sum();

        flx_med_sm[i]=flx_sm*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g_sm[st][i][j]=mic_sigf_sm[bstmp][i][j].Cond(flx_sm);
  	      };
	      xsc_1g_sm[st][i][j]=mic_sigc_sm[bstmp][i][j].Cond(flx_sm);
	      xsn2n_1g_sm[st][i][j]=mic_sign2n_sm[bstmp][i][j].Cond(flx_sm);
            };
  	  };
        };
      };
  #endif      


  #if 1
      // ... with group collapsing from 107-group to fewer-group ...
      // +++ Instance generation of MECSystem to solve neutron transport equation for simplified system
      

      //int ngrp=21; // The number of energy groups after the group collapsing
      //int bgrp[]={4,9,14,19,24, 29,34,39,44,49, 54,59,64,69,74, 79,84,89,94,99, 106};
      // Energy group boundary in the collapsed structure.
      // The final one should be [106] when the number of groups in the original structure is the 107-group.
      
      // int ngrp=11;
      // int bgrp[]={9,19,29,39,49, 59,69,79,89,99, 106};
      /*
      int ngrp=107;
      int *bgrp=new int[ngrp];
      for(int i=0;i<ngrp;i++){
	bgrp[i]=i;
      };
      */

      vector< vector<GroupData1D> > mic_sigf_sm(mednum_fuel);  // [medid][nucn](group)
      vector< vector<GroupData1D> > mic_sigc_sm(mednum_fuel);
      vector< vector<GroupData1D> > mic_sign2n_sm(mednum_fuel); 
	
      MECSystem lat2(ngrp,mednum);
      lat2.PutTrajectorySet(&sys_f2);        

      // +++ fuel medium data preparation
      for(int i=0;i<mednum_fuel_sm;i++){

	/*
        // (Self-shielding calculation) ---------------------------------------
        PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
	  };
        };
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);     
        med[0].CalMacroFromMicro();

        Medium bmed=med[0].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
        // ------------------------------------------------------------------------
	*/

	// (Use of the original burnup step results) -------------------------------
        int nucn=med_adts[i].GetNucnum();
        for(int ii=0;ii<nucn;ii++){
           med_adts[i].GetNuclideInTurn(ii).PutDensity(fwd_nuc[st][0][i].get_dat(ii));
        };
	med_adts[i].CalMacroFromMicro();
	Medium bmed=med_adts[i].Cond(ngrp,bgrp,flx_med[i],flx_med[i],true); // [true] is for microscopic XS
	// -------------------------------------------------------------------------


        // (Macro&Micro data storing)
  mic_sigf_sm[i].resize(nucn);
  mic_sigc_sm[i].resize(nucn);
  mic_sign2n_sm[i].resize(nucn);
        for(int k=0;k<nucn;k++){
	  
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_sm[i][k].copy(bmed.GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };

        lat2.AddMedium(bmed); 
        lat2.GetMedium(i).NuclideClear();
      };

      // +++ non-fuel medium data preparation
      for(int jj=0;jj<mednum_nonfuel;jj++){
        Medium bmed=med[1+jj].Cond(ngrp,bgrp,flx_med[mednum_fuel+jj],flx_med[mednum_fuel*jj],true); // [true] is for microscopic XS
        lat2.AddMedium(bmed);
      };

      // +++ eigenvalue calculation
      lat2.PutRegMed(region_medium);
      lat2.PutGeneralOption(opt);
      lat2.PutThermalIteration(3);
      lat2.PutPL(0);
      lat2.NoCMRAcceleration();
      lat2.Print();

      real keff_dummy=lat2.CalIgen();

      vector<real> medium_wise_flux_change(mednum);
      
      for(int i=0;i<mednum;i++){
        GroupData1D flx_sm=lat2.GetIntegratedFlux(i);

	// ... Total neutron flux change during the subsequent ADTS steps is calculated.
	medium_wise_flux_change[i]=(flx_sm.get_sum()/vol_med[i])/flx_med_sm[i].get_sum();

        flx_med_sm[i]=flx_sm*(1./vol_med[i]); // medium-wise neutron flux data is stored here
        // +++ One-group cross section storing
        if(i<mednum_fuel_sm){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1){
   	        xsf_1g_sm[st][i][j]=mic_sigf_sm[i][j].Cond(flx_sm);
  	      };
	      xsc_1g_sm[st][i][j]=mic_sigc_sm[i][j].Cond(flx_sm);
	      xsn2n_1g_sm[st][i][j]=mic_sign2n_sm[i][j].Cond(flx_sm);
            };
  	  };
        };
      };
  #endif      



      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      if(st%sm_step==0){	
        //st_10=st;
      }else{

        // ... One-group cross section correction in ADTS	
             for(int i=0;i<mednum_fuel;i++){
	  for(int k=0;k<nucn;k++){
	  
  	    int nucid=med[0].GetNuclideInTurn(k).GetMatnum();
    #if 0	  
            // ... Copy of the ADTS step data
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
  	        xsf_1g[st][i][k]=xsf_1g_sm[st][i][k];
	      };
	      xsc_1g[st][i][k]=xsc_1g_sm[st][i][k];
	      if(xsn2n_1g[st-1][i][k]>0.){
  	        xsn2n_1g[st][i][k]=xsn2n_1g_sm[st][i][k];
	      };
	    };
    #endif	  

    #if 1  
            // ... Correcting all nuclides in all burnup regions
	    if(nuclide_info[k]!=0){
	      if(nuclide_info[k]==1){
  	        xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_sm[st][i][k]/xsf_1g_sm[st-1][i][k]); //for neutron fission cross section
	      };
	      xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k]);//for neutron capture cross section
	      if(xsn2n_1g[st-1][i][k]>0.){
  	        xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_sm[st][i][k]/xsn2n_1g_sm[st-1][i][k]); // for n_2n reaction cross section
	      };
	    };
    #endif	  

    #if 0	  
	    if(i<8){
              // ... Correcting all nuclides in the central Gd-bearing pin
	      if(nuclide_info[k]!=0){
	        if(nuclide_info[k]==1){
  	          xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_sm[st][i][k]/xsf_1g_sm[st-1][i][k]); //for neutron fission cross section
	        };
	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k]);//for neutron capture cross section
  	        if(xsn2n_1g[st-1][i][k]>0.){	      
	          xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_sm[st][i][k]/xsn2n_1g_sm[st-1][i][k]); // for n_2n reaction cross section
	        };
	      };
	    };
    #endif	  

    #if 0	  
	    if(i<8){
	      // ... Correcting Gd-155 and -157 capture cross sections in the central Gd-bearing pin
	      if(nucid==641550||nucid==641570){
	        if(nucid==641550) cout<<"Gd-155 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        if(nucid==641570) cout<<"Gd-157 "<< xsc_1g[st][i][k]<<" "<<xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k])<<endl;
	        //xsf_1g[st][i][k]=xsf_1g[st-1][i][k]*(xsf_1g_2[st][i][k]/xsf_1g_2[st-1][i][k]); //for neutron fission cross section
	        xsc_1g[st][i][k]=xsc_1g[st-1][i][k]*(xsc_1g_sm[st][i][k]/xsc_1g_sm[st-1][i][k]);//for neutron capture cross section
	        //xsn2n_1g[st][i][k]=xsn2n_1g[st-1][i][k]*(xsn2n_1g_2[st][i][k]/xsn2n_1g_2[st-1][i][k]); // for n_2n reaction cross section
	      };
	    };
    #endif	  

	  };
        };

	// ... Neutron flux correction in ADTS

	for(int i=0;i<mednum_fuel;i++){ // All the medium included in the system
	  //for(int i=0;i<8;i++){ // Only for Gd-bearing rod
	  flx_med[i]=flx_med[i]*medium_wise_flux_change[i];
	};
      };
        

    // fuga addition -------------------------------------
      //output data

      // One-group cross section
      fout_1g.setf(ios::scientific);
      fout_1g.precision(6);
      fout_1g<<acburn[st]<<" ";
      for(int i=0;i<mednum_fuel;i++){
        for(int j=0;j<nucn;j++){
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    for(int ll=0;ll<sz;ll++){
      if(nucid==id[ll]){
        if(ncl_f){fout_1g<<xsf_1g[st][i][j]<<" "<<xsf_1g_sm[st][i][j]<<" ";};
        if(ncl_c){fout_1g<<xsc_1g[st][i][j]<<" "<<xsc_1g_sm[st][i][j]<<" ";};
        if(ncl_n2n){fout_1g<<xsn2n_1g[st][i][j]<<" "<<xsn2n_1g_sm[st][i][j]<<" ";};
      };
    };
        };
      };
      fout_1g<<endl;

      // Total neutron flux
      fout_tf.setf(ios::scientific); 
      fout_tf.precision(6); 
      fout_tf<<acburn[st]<<" ";
      for(int i=0;i<mednum_fuel;i++){
    fout_tf<<flx_med[i].get_sum()<<" "<<flx_med_sm[i].get_sum()<<" ";
      };
      fout_tf<<endl;


      // cross section in Original-Model 
      fout_om<<"#\n# Neutron flux per lethargy\n#\n";
      fout_om<<"#  burnup step : "<<st<<endl;
      fout_om<<"#  burnup[GWD/t] : "<<acburn[st]<<endl;
      fout_om<<"#\n";
      fout_om<<"# Energy    Med ID\n";
      fout_om<<"# [eV]     ";
      for(int i=0;i<mednum_fuel;i++){
        string strid=IntToString(i);
        fout_om<<" ("+strid+")";
        for(int j=0;j<(10*(1+sz*react)-3-strid.size());j++){fout_om<<" ";};
      };
      fout_om<<endl;

      fout_om.setf(ios::scientific);
      fout_om.precision(3);
      for(int i=0;i<group;i++){
        real e0=med[0].GetEnband().get_dat(i);
        real e1=med[0].GetEnband().get_dat(i+1);
        real letwid=log(e0/e1);
        fout_om<<e0<<" ";

        for(int j=0;j<mednum_fuel;j++){
          fout_om<<flx_med[j].get_dat(i)/letwid<<" ";
          for(int k=0;k<nucn;k++){
            int nucid=med[0].GetNuclideInTurn(k).GetMatnum();

            if(nuclide_info[k]!=0){
      for(int ll=0;ll<sz;ll++){

        if(nucid==id[ll]){
          if(ncl_f){
            if(nuclide_info[k]==1){
              fout_om<<mic_sigf[bstmp][j][k].get_dat(i)<<" ";
            }else{
              fout_om<<0.<<" ";
            };
          };
          if(ncl_c){fout_om<<mic_sigc[bstmp][j][k].get_dat(i)<<" ";};
          if(ncl_n2n){fout_om<<mic_sign2n[bstmp][j][k].get_dat(i)<<" ";};
        };

      };  
            };

          };
        };
        fout_om<<endl;

      };
      fout_om<<"\n\n";
      

      // cross section in Simplified-Model  
      fout_sm<<"#\n# Neutron flux per lethargy\n#\n";
      fout_sm<<"#  burnup step : "<<st<<endl;
      fout_sm<<"#  burnup[GWD/t] : "<<acburn[st]<<endl;
      fout_sm<<"#\n";
      fout_sm<<"# Energy    Med ID\n";
      fout_sm<<"# [eV]    ";
      for(int i=0;i<mednum_fuel_sm;i++){
        string strid=IntToString(i);
        fout_sm<<" ("+strid+")";
        for(int j=0;j<(10*(1+sz*react)-3-strid.size());j++){fout_sm<<" ";};
      };
      fout_sm<<endl;

      fout_sm.setf(ios::scientific);
      fout_sm.precision(3);
      for(int i=0;i<ngrp;i++){
        Medium bmed=med_adts[0].Cond(ngrp,bgrp,flx_med[0],flx_med[0],true);
        real e0=bmed.GetEnband().get_dat(i);
        real e1=bmed.GetEnband().get_dat(i+1);
        real letwid=log(e0/e1);
        fout_sm<<e0<<" ";

        for(int j=0;j<mednum_fuel_sm;j++){
          fout_sm<<flx_med_sm[j].get_dat(i)/letwid<<" ";
          // Medium bmed=med_adts[j].Cond(ngrp,bgrp,flx_med[j],flx_med[j],true);
          for(int k=0;k<nucn;k++){
            int nucid=bmed.GetNuclideInTurn(k).GetMatnum();
      
            if(nuclide_info[k]!=0){
      for(int ll=0;ll<sz;ll++){

        if(nucid==id[ll]){
          if(ncl_f){
            if(nuclide_info[k]==1){
              fout_sm<<mic_sigf_sm[j][k].get_dat(i)<<" ";
            }else{
              fout_sm<<0.<<" ";
            };
          };
          if(ncl_c){fout_sm<<mic_sigc_sm[j][k].get_dat(i)<<" ";};
          if(ncl_n2n){fout_sm<<mic_sign2n_sm[j][k].get_dat(i)<<" ";};
        };

      };  
            };

          };
        };
        fout_sm<<endl;
      };
      fout_sm<<"\n\n";

      // fout_fai.setf(ios::scientific);
      // fout_fai.precision(s);
      // int grp_sm=ngrp;
      // grp_sm--;
      // real aveflx_sm;
      // for(int i=group-1;i>=0;i--){
      //   real e0=med[0].GetEnband().get_dat(i);
      //   real e1=med[0].GetEnband().get_dat(i+1);
      //   real letwid=log(e0/e1);
      //   if(bgrp[group_sm]==gourp){
      //       Medium bmed=med_adts[0].Cond(ngrp,bgrp,flx_med[0],flx_med[0],true);
      //       real e0=bmed.GetEnband().get_dat(i);
      //       real e1=bmed.GetEnband().get_dat(i+1);
      //       real letwid=log(e0/e1);
            
      //   }
      //   fout_fai<<e0<<" ";

      //   for(int j=0;j<mednum_fuel;j++){
      //     fout_fai<<flx_med[j].get_dat(i)/letwid<<" ";
          
      //   }


        
        
      


    //---------------------------------------------------------



      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  /*      
      // Commented out in this method since it is unnecessary.
      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };
  */

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";
    
      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          //real tmp=flx_med_2[i].get_sum();	  
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
  	  //sumflx=flx_med_2[med_target].get_sum();	  
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
	
	for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
	  //total_flux[st][j][i]=flx_med_2[i].get_sum()*power_factor[st][j];	  
	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	};
	
      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  }; // loop-end of burnup step
};

void MulticellBurner::ForwardCalculationNEL2020(Burnup &bu, int med_target, int sm_step, int ngrp, int *bgrp)
{
  // Pre-determined parameters/constants
  int num_33_system=2;
  int gd_id_list[]={16,40};
  int uo2_id_list[]={24,48};
  vector<int> region_medium_33;
  int region_medium_num=9+12+8*4*2;
  int region_medium_33_inp[]={
    10,10,10,10,10, 10,10,10,10, // background meshes(9)
    10,10,10,9,0,1,2,3,4,5,6,7,  // center pin(12)
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // right pin
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8,
    10,10,10,9,8,8,8,8,    10,10,10,9,8,8,8,8, // top pin    
  };
  for(int i=0;i<region_medium_num;i++){
    region_medium_33.push_back(region_medium_33_inp[i]);
  };
  
  
  
  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculationNEL2020.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor in pincell model +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // +++ Initialization of Bell factor +++
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting +++
  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  flx_med.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  {
  int tmp=1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //  RUNNING BURNUP CALCULATION 
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";


    // +++ SOLVING NEUTRON TRANSPORT EQUATION +++++++++++++++++++++++++++++++++++++++

    // +++ Instance generation of MECSystem to solve neutron transport equation
    MECSystem lat(group,mednum);
    lat.PutTrajectorySet(&sys_f);  

    // +++ fuel medium data preparation
    for(int i=0;i<mednum_fuel;i++){

      // (Self-shielding calculation)
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
      };
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048);      
      med[0].CalMacroFromMicro();

      // (Macro&Micro data storing)
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };

      lat.AddMedium(med[0]); 
      lat.GetMedium(i).NuclideClear();

    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    lat.Print();
    keff[st]=lat.CalIgen();
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);
      flx_med[i]=flx*(1./vol_med[i]); // medium-wise neutron flux data is stored here
      // +++ One-group cross section storing
      if(i<mednum_fuel){
        for(int j=0;j<nucn;j++){
          if(nuclide_info[j]!=0){
            if(nuclide_info[j]==1){
   	      xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	    };
	    xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	    xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
          };
	};
      };
    };


    // +++ Neutron transport calculations for 3x3 multicell including Gd-bearing rod at a center
    for(int ii=0;ii<num_33_system;ii++){
    
    MECSystem lat2(group,11); // "11" is the number of mediums and hard-coded
    lat2.PutTrajectorySet(&sys_f2);  

    // +++ fuel medium data preparation
    for(int i=0;i<9;i++){
      int medid=gd_id_list[ii]+i;
      if(i==8)medid=uo2_id_list[ii];
      lat2.AddMedium(lat.GetMedium(medid));
    };

    // +++ non-fuel medium data preparation
    for(int jj=0;jj<mednum_nonfuel;jj++){
      lat2.AddMedium(med[1+jj]);
    };

    // +++ eigenvalue calculation
    lat2.PutRegMed(region_medium_33);
    lat2.PutGeneralOption(opt);
    lat2.PutThermalIteration(3);
    lat2.PutPL(0);
    lat2.NoCMRAcceleration();
    lat2.Print();
    real keff_dummy=lat2.CalIgen();

    vector<GroupData1D> flx_med_33(8);
    vector<real> gd5_xsc1g(8);
    vector<real> gd7_xsc1g(8);
    int pos_gd5, pos_gd7;
    for(int i=0;i<8;i++){
      GroupData1D flx=lat2.GetIntegratedFlux(i);
      flx_med_33[i]=flx*(1./vol_med[gd_id_list[ii]+i]); // medium-wise neutron flux data is stored here      
      // +++ One-group cross section calculation
      for(int j=0;j<nucn;j++){
        if(med[0].GetNuclideInTurn(j).GetMatnum()==641550){
	  gd5_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd5=j;
	}else if(med[0].GetNuclideInTurn(j).GetMatnum()==641570){
	  gd7_xsc1g[i]=mic_sigc[bstmp][gd_id_list[ii]+i][j].Cond(flx);
	  pos_gd7=j;
	};
      };
    };

    // normalization factor calculation
    /*
    real tflx_ref=0.;
    real tflx_33=0.;
    for(int i=0;i<8;i++){
      tflx_ref+=flx_med[gd_id_list[ii]+i].get_sum()*vol_med[gd_id_list[ii]+i];
      tflx_33+=flx_med_33[gd_id_list[ii]+i].get_sum()*vol_med[gd_id_list[ii]+i];      
    };
    real norm_factor=tflx_ref/tflx_33;
    */
    
    cout.setf(ios::scientific);
    cout.precision(8);
    cout<<"# Total neutron flux (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<flx_med[gd_id_list[ii]+i].get_sum()<<" ";
      cout<<flx_med_33[i].get_sum()<<" ";
      cout<<"\n";
    };
    cout<<"\n";

    cout<<"# Gd-155 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd5]<<" ";
      cout<<gd5_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    cout<<"# Gd-157 one-group (n,g) cross section (original & 3x3)\n";
    for(int i=0;i<8;i++){
      cout<<xsc_1g[st][gd_id_list[ii]+i][pos_gd7]<<" ";
      cout<<gd7_xsc1g[i]<<" ";
      cout<<"\n";
    };
    cout<<"\n";
    
    };
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    // +++ SOLVING BURNUP EQUTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
	};
      };

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};
        accumulated_day+=burn_span;
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],false);
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  }; // loop-end of burnup step
};

void MulticellBurner::ShowEigenvalue()
{
  cout<<"#\n# +++ Time-dependent eigenvalue +++\n";
  cout<<"#  (day)       ";
  cout<<"(GWd/t)     ";
  cout<<"(keff)      ";
  //cout<<"(C.R.)      ";
  //cout<<"(flux[/cm2/s])";
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(5);
  //cout.precision(10);
  for(int i=0;i<burn_step+1;i++){
    if(i<10)cout<<" ";
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" "<<keff[i]<<" ";
    cout<<"\n";
  };
};

real MulticellBurner::CalculationPinPower(Burnup &bu, int st, int sst, int medid, real vol_flx)
{
  bu.PutNucnum(nucn);
  for(int k=0;k<nucn;k++){
    int id=med[0].GetNuclideInTurn(k).GetMatnum();
    real den=fwd_nuc[st][sst][medid].get_dat(k);
    bu.PutNuclideData(k,id,den,xsf_1g[st][medid][k],xsc_1g[st][medid][k],xsn2n_1g[st][medid][k]);
  };
  return bu.GetIntegratedPower(vol_flx,true);
};

void MulticellBurner::CalculationPinBurnup(Burnup &bu, int st, int sst, int medid, vector<real> &xsf1g, vector<real> &xsc1g, vector<real> &xsn2n1g, real tflx, real period, bool mstep)
{
  bu.PutNucnum(nucn);
  for(int k=0;k<nucn;k++){
    int id=med[0].GetNuclideInTurn(k).GetMatnum();
    real den=fwd_nuc[st][sst][medid].get_dat(k);
    //cout<<k<<" "<<id<<" "<<den<<" "<<xsf1g[k]<<" "<<xsc1g[k]<<" "<<xsn2n1g[k]<<"\n";
    bu.PutNuclideData(k,id,den,xsf1g[k],xsc1g[k],xsn2n1g[k]);
  };
  bu.CalTransitionMatrixFluxDependentPart();
  bu.CalTransitionMatrixFluxInDependentPart();
  if(mstep){
    bu.BurnupCalculationByMultiStepCalc(med[0],fwd_nuc_int[st][sst][medid],tflx,period,false);
  }else{
    bu.BurnupCalculationByMMPA(med[0],tflx,period,false);
  };
  for(int k=0;k<nucn;k++){
    real tmp=bu.GetDensity(k);
    if(tmp<0.){
      cout<<"# Warning ... \n";
      cout<<"#   Negative number density is detected in "<<med[0].GetNuclideInTurn(k).GetMatnum()<<" and reset to zero.\n";
      tmp=1e-20; // Zero-reset in negative ND case (1e-20 for log-PC method)
    };
    //cout<<k<<" "<<tmp<<"\n";
    fwd_nuc[st][sst+1][medid].put_data(k,tmp);
  };
};

void MulticellBurner::PutDancoffFactor(real *dancoff_inp)
{
  dancoff_input=true;
  if(mednum_fuel==0){
    cout<<"# Error in MulticellBurner::PutDancoffFactor.\n";
    cout<<"# No information on medium has not yet been given.\n";
    exit(0);
  };

  dancoff_factor.resize(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    dancoff_factor[i]=dancoff_inp[i];
    if(dancoff_factor[i]<0.||dancoff_factor[i]>1.){
      cout<<"# Error in MulticellBurner::PutDancoffFactor.\n";
      cout<<"# Inpropoer Dancoff factor is assigned.\n";
      cout<<"# Medium ID is "<<i<<" and Dancoff factor is "<<dancoff_factor[i]<<"\n";
      exit(0);
    };
  };
  
};

void MulticellBurner::GroupCollapsingDuringBurnup(string dirname)
{
  collapsing=true;
  collapse_dir="./"+dirname+"/";
};

void MulticellBurner::PreEigenvalueCalculation(Burnup &bu)
{
  PreCalculation(bu);

  GeneralOption opt;

  fwd_nuc.resize(1);
  xsc_1g.resize(1);
  xsn2n_1g.resize(1);
  xsf_1g.resize(1);
  total_flux.resize(1);

  fwd_nuc[0].resize(1);
  fwd_nuc[0][0].resize(mednum_fuel);
  for(int k=0;k<mednum_fuel;k++){
    fwd_nuc[0][0][k].put_imax(nucn);
  };

  total_flux[0].resize(1);
  total_flux[0][0].resize(mednum_fuel);

  xsc_1g[0].resize(mednum_fuel);
  xsn2n_1g[0].resize(mednum_fuel);
  xsf_1g[0].resize(mednum_fuel);
  for(int j=0;j<mednum_fuel;j++){
    xsc_1g[0][j].resize(nucn,0.);
    xsn2n_1g[0][j].resize(nucn,0.);
    xsf_1g[0][j].resize(nucn,0.);
  };

  mic_sigf.resize(1); 
  mic_sigc.resize(1); 
  mic_sign2n.resize(1);
  mic_sigf[0].resize(mednum_fuel);
  mic_sigc[0].resize(mednum_fuel);
  mic_sign2n[0].resize(mednum_fuel);
  for(int j=0;j<mednum_fuel;j++){
    mic_sigf[0][0].resize(nucn);
    mic_sigc[0][0].resize(nucn);
    mic_sign2n[0][0].resize(nucn);
    for(int k=0;k<nucn;k++){
      if(nuclide_info[k]!=0){
        if(nuclide_info[k]==1){
          mic_sigf[0][0][k].put_imax(group);
        };
        mic_sigc[0][0][k].put_imax(group);
        mic_sign2n[0][0][k].put_imax(group);
      };
    };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
  };

};

real MulticellBurner::EigenvalueCalculation(Burnup &bu)
{
  int st=0;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);  

  for(int i=0;i<mednum_fuel;i++){
    // +++ Self-shielding calculation
    PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    if(dancoff_input){
      for(int g=0;g<group;g++){
        dancoff.put_data(g,1.-dancoff_factor[i]);
      };
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    }else{
      opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    };
    //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    opc.CalThermalScatteringMatrix(med[0],xslib,4.048);    
    med[0].CalMacroFromMicro();

    // +++ Micro data storing
    for(int k=0;k<nucn;k++){
      if(nuclide_info[k]!=0){
        if(nuclide_info[k]==1){
          mic_sigf[0][0][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
        };
        mic_sigc[0][0][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
        mic_sign2n[0][0][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
      };
    };

    lat.AddMedium(med[0]); 
    lat.GetMedium(i).NuclideClear();
  };

  GeneralOption opt;

  for(int jj=0;jj<mednum_nonfuel;jj++){
    lat.AddMedium(med[1+jj]);
  };
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);

  lat.PutThermalIteration(3);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  //lat.NoTransportApprox();
  lat.Print();
  real keff=lat.CalIgen();

  for(int i=0;i<mednum_fuel;i++){
    GroupData1D flx=lat.GetIntegratedFlux(i);
    total_flux[0][0][i]=flx.get_sum();
    for(int j=0;j<nucn;j++){
      if(nuclide_info[j]!=0){
        if(nuclide_info[j]==1){
 	  xsf_1g[0][0][j]=mic_sigf[0][0][j].Cond(flx);
	};
        xsc_1g[0][0][j]=mic_sigc[0][0][j].Cond(flx);
        xsn2n_1g[0][0][j]=mic_sign2n[0][0][j].Cond(flx);
      };
    };
  };

  return keff;
};

void MulticellBurner::PerturbNumberDensity(int nucid, int medid, real factor)
{
  int pos=-1;
  for(int i=0;i<nucn;i++){
    if(med[0].GetNuclideInTurn(i).GetMatnum()==nucid)pos=i;
  };

  if(pos==-1){
    cout<<"# Error in MulticellBurner::PerturbNumberDensity.\n";
    exit(0);
  };

  real org=fwd_nuc[0][0][medid].get_dat(pos);
  fwd_nuc[0][0][medid].put_data(pos,org*factor);
};

real MulticellBurner::GetMicroscopicReactionRate(int nucid, int medid)
{
  int pos=-1;
  for(int i=0;i<nucn;i++){
    if(med[0].GetNuclideInTurn(i).GetMatnum()==nucid)pos=i;
  };

  if(pos==-1){
    cout<<"# Error in MulticellBurner::GetMicroscopicReactionRate.\n";
    exit(0);
  };

  return xsc_1g[0][0][pos]*total_flux[0][0][medid];
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize,int med_target,int num,string *nucname_target_list,string filename)
{
  int *nucid_target_list=new int[num];
  for(int i=0;i<num;i++){
    nucid_target_list[i]=midt.ID(nucname_target_list[i]);
  };

  int num_med=1;
  int *med_target_list=new int[num_med];
  for(int i=0;i<num_med;i++){ 
    med_target_list[i]=med_target;
  };
 
  SensitivityCalculation(bu,med_normalize,num_med,med_target_list,num,nucid_target_list,filename);

  delete [] nucid_target_list;
  delete [] med_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize,int num,string *nucname_target_list,string filename)
{
  int *nucid_target_list=new int[num];
  for(int i=0;i<num;i++){
    nucid_target_list[i]=midt.ID(nucname_target_list[i]);
  };

  SensitivityCalculation(bu,med_normalize,num,nucid_target_list,filename);

  delete [] nucid_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize,int med_target,int num,int *nucid_target_list,string filename)
{
  int num_med=1;
  int *med_target_list=new int[num_med];
  for(int i=0;i<num_med;i++){ 
    med_target_list[i]=med_target;
  };
 
  SensitivityCalculation(bu,med_normalize,num_med,med_target_list,num,nucid_target_list,filename);
  delete [] med_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize,int num,int *nucid_target_list,string filename)
{
  int num_med=mednum_fuel;
  int *med_target_list=new int[num_med];
  for(int i=0;i<num_med;i++){ 
    med_target_list[i]=i;
  };
 
  SensitivityCalculation(bu,med_normalize,num_med,med_target_list,num,nucid_target_list,filename);
  delete [] med_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize, int num_med, int *med_target,int num,string *nucname_target_list,string filename)
{
  int *nucid_target_list=new int[num];
  for(int i=0;i<num;i++){
    nucid_target_list[i]=midt.ID(nucname_target_list[i]);
  };

  SensitivityCalculation(bu,med_normalize,num_med,med_target,num,nucid_target_list,filename);

  delete [] nucid_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize, int num_med, int *med_target,int num,int *nucid_target_list,string filename)
{
  if(corrector_calc){
    SensitivityCalculationPC(bu,med_normalize,num_med,med_target,num,nucid_target_list,filename);
    return;
  };

  bool isotropic_approx=false; // Isotropic approximation in perturbation integral

  int ssv=40;  // (sub-sub-time step division)

  // +++

  vector<bool> med_target_bool(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    med_target_bool[i]=false;
  };

  for(int i=0;i<num_med;i++){
    if(med_target[i]>=mednum_fuel){
      cout<<"# Error in MulticellBurner::SensitivityCalculation(PC).\n";
      cout<<"# Medium index "<<med_target[i]<<" is inappopriate.\n";
      exit(0);
    };
    med_target_bool[med_target[i]]=true;
  };

  // +++

  PreCalculation(bu);

  // +++ Array setting for adjoint calculation +++
  pow_adj.resize(burn_step);
  adj_nuc.resize(burn_step);
  macxs.resize(burn_step);

  bilinear_flx.resize(burn_step);
  bilinear_aflx.resize(burn_step);
  chi_gpt_fwdflx.resize(burn_step);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];

    pow_adj[i].resize(sub_step);
    adj_nuc[i].resize(sub_step);
    for(int j=0;j<sub_step;j++){
      adj_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
        adj_nuc[i][j][k].put_imax(nucn); 
      };
    };

    bilinear_flx[i].resize(mednum_fuel);
    bilinear_aflx[i].resize(mednum_fuel);
    chi_gpt_fwdflx[i].resize(mednum_fuel);
    for(int k=0;k<mednum_fuel;k++){
      bilinear_flx[i][k].put_imax(group);
      bilinear_aflx[i][k].put_imax(group);
      chi_gpt_fwdflx[i][k].put_imax(group);
    };

    macxs[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      macxs[i][j].Init("MacroCrossSection");
      macxs[i][j].PutGrp(107);
    };
  };

  // PJI
  /*
  if(pij_storing){
    pij_store.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      pij_store[i].resize(1);
      pij_store[i][0].resize(group);
    };
  };
  */
  // MEC
  // noprocessing

  // +++ Forward burnup calculation +++
  ForwardCalculation(bu,med_normalize,true);

  // +++ Integrating forward number density
  // [fwd_nuc_int] is a number density at time-mesh-center point
  // time-averagted number density is calculated 
  // from those at beginning, center and end of the time mesh

  IntegratingForwardNumberDensity(fwd_nuc,fwd_nuc_int);

  // +++ Adjoint burnup calculation +++

  GeneralOption opta;
  opta.PutAdjointCal();

  for(int ii=0;ii<num;ii++){

    int nucid_target=nucid_target_list[ii];
    int nucid_turn=med[0].SearchNuclide(nucid_target);
    if(nucid_turn==-1){
      cout<<"# Error !\n";
      exit(0);
    };
    //real response=fwd_nuc[burn_step][0][med_target].get_dat(nucid_turn);
    real response=0.;
    real volsum=0.;
    for(int jj=0;jj<mednum_fuel;jj++){
      if(med_target_bool[jj]){
        real nd=fwd_nuc[burn_step][0][jj].get_dat(nucid_turn);
	response+=nd*vol_med[jj];
	volsum+=vol_med[jj];
      };        
    };
    response/=volsum;

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response ("<<midt.Name(nucid_target)<<" ND) : "<<response<<"\n";

  // +++ Adjoint calculation

  MECSystem lat_gpt(group,mednum);

  lat_gpt.NoPrint();
  lat_gpt.PutTrajectorySet(&sys_f);
  int totsn=lat_gpt.GetQuad().GetSN();

  for(int i=0;i<mednum_fuel;i++){
    lat_gpt.AddMedium(med[0]);
    lat_gpt.GetMedium(i).NuclideClear();
  };
  for(int jj=0;jj<mednum_nonfuel;jj++){
    lat_gpt.AddMedium(med[1+jj]);
  };
  lat_gpt.PutRegMed(region_medium);
  lat_gpt.PutGeneralOption(opta);
  lat_gpt.PutPL(0);

  lat_gpt.NoCMRAcceleration();
  //lat_gpt.NoTransportApprox();
  lat_gpt.PutWriteFlux();
  lat_gpt.SetArray();

  vector<GroupData1D> gpt_flx(totm);
  vector< vector<GroupData1D> > gpt_aflx(totm);
  for(int i=0;i<totm;i++){
    gpt_flx[i].put_imax(group);
    gpt_aflx[i].resize(group);
    for(int g=0;g<group;g++){
      gpt_aflx[i][g].put_imax(totsn);
    };
  };

  vector<GroupData1D> gpt_src(mednum_fuel);
  vector<GroupData1D> gpt_src2(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    gpt_src[i].put_imax(group);
    gpt_src2[i].put_imax(group);
  };

  vector<GroupData1D> adj_nuc_e(mednum_fuel);

  // +++ Final condition for adjoint number density

  //real total_decay_energy=0.;
  for(int i=0;i<mednum_fuel;i++){
    adj_nuc_e[i].put_imax(nucn);
    adj_nuc_e[i].set_zero();
    for(int j=0;j<nucn;j++){
      int id=med[0].GetNuclideInTurn(j).GetID();

      // For number density
      /*
      if(id==nucid_target&&i==med_target){
	adj_nuc_e[i].put_data(j,1.0);
      };
      */
      if(id==nucid_target&&med_target_bool[i]){
	adj_nuc_e[i].put_data(j,vol_med[i]/volsum);
      };

      // For decay heat
      /*
      real ev_to_j=1.60219e-19;
      real decay_const=bu.GetBurnupChain().GetDecayConstant(id);
      real dens=fwd_nuc[burn_step][0][i].get_dat(j);
      real decay_energy=0.;
      for(int k=0;k<3;k++){
        decay_energy+=bu.GetBurnupChain().GetDecayEnergy(id,k);
      };
      real tmp=decay_energy*decay_const*ev_to_j*1e+24;
      total_decay_energy+=tmp*dens;
      adj_nuc_e[i].put_data(j,tmp);
      */
    };
  };

  //response=total_decay_energy; // for decay heat

  // +++ Averaging forward number density
  /*
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      for(int m=0;m<mednum_fuel;m++){
        for(int k=0;k<nucn;k++){
  	  //fwd_nuc[i][j].put_data(k,(fwd_nuc[i][j].get_dat(k)+fwd_nuc[i][j+1].get_dat(k))*0.5);
	  fwd_nuc[i][j][m].put_data(k,sqrt(fwd_nuc[i][j][m].get_dat(k)*fwd_nuc[i][j+1][m].get_dat(k))*0.5);
	};
      };
    };
  };
  */


  GroupData2D trmat_flxindep=bu.GetTrmatFlxInDep();
  for(int st=burn_step-1;st>=0;st--){

    real power_density=power_density_list[st];
    int sub_step=sub_step_list[st];

    cout<<"#   Adjoint calculation step : "<<st<<"\n";

    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross sectio data
      // is updated in the preceding lines
    };

    real keff_step=keff[st];

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculation.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g[st][m][i],xsc_1g[st][m][i],xsn2n_1g[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        //pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc[st][j][m]+fwd_nuc[st][j+1][m]))*0.5*delt[st][j];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*fwd_nuc[st][j][m])*delt[st][j];
	// ++++++++
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        //pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
        pow_adj[st][j]/=flux_level_list[st]*volsum;
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
  	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);

        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor[st][j];
	        };
	      };
	    };
	  };
        }else{
	  //if(m==med_target){
	  if(med_target_bool[m]){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc[st][j][m]))*delt[st][j]*power_factor[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    // ++++
    int sntot=lat_gpt.GetQuad().GetSN();
    int pdiv=lat_gpt.GetPolarAngleDivision();
    int sn_per_pang=sntot/pdiv;

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh[st][mm]));
          real tmp=gpt_flx[mm]*macxs[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
              real tmp=omega*volaflx_mesh[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // jump condition for generalized adjoint
    if(st!=0){
      for(int m=0;m<mednum_fuel;m++){
        for(int k=0;k<nucn;k++){
          real xsf1g=xsf_1g[st][m][k];

          if(nuclide_info[k]!=0){
            // (absorption or total term)
            real tmp2=0.;
            for(int g=0;g<group;g++){
              real xsa=mic_sigc[st][m][k].get_dat(g);
              if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
  	      if(isotropic_approx){
                tmp2+=xsa*bilinear_flx[st][m].get_dat(g); // absorption
	      }else{
  	        tmp2+=xsa*bilinear_aflx[st][m].get_dat(g); // total
	      };
	    };
	    // (yield term)
	    real tmp3=0.;
	    if(xsf1g>0.){
	      for(int mm=0;mm<totm;mm++){
	        if(region_medium[mm]==m){
		  real fsrc=0.;
		  for(int g=0;g<group;g++){
		    fsrc+=volflx_mesh[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
		      *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);
  	  	  };
		  for(int g=0;g<group;g++){
		    tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs[st][m].GetData1d(chi).get_dat(g);
		  };
	        };
	      };
  	      tmp3/=keff_step;
  	    };
            adj_nuc_e[m].add_data(k,tmp2-tmp3);
	  };
        };
      };
    };

  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Sensitivity Printing

  SensitivityData sns;
  sns.PutName("dummy","dummy","dummy");
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

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
      vector< vector< vector<real> > > nadj_dm_nfwd(burn_step);
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        nadj_dm_nfwd[k].resize(sub_step);
        for(int l=0;l<sub_step;l++){
	  nadj_dm_nfwd[k][l].resize(mednum_fuel);
	  for(int m=0;m<mednum_fuel;m++){
            real val=0.;
            val+=-fwd_nuc[k][l][m].get_dat(i)*adj_nuc[k][l][m].get_dat(i);
            int tmp=bu.GetBC().GetNdiv(matnum,bc_channel);
            for(int j=0;j<tmp;j++){
              int id2=bu.GetBC().GetNextID(matnum,bc_channel,j);
              int pos=bu.SearchNuclide(id2);
              if(pos!=-1){
                real rat=bu.GetBC().GetRatio(matnum,bc_channel,j);
                val+=adj_nuc[k][l][m].get_dat(pos)*rat*fwd_nuc[k][l][m].get_dat(i);
              };
            };
            nadj_dm_nfwd[k][l][m]=val*delt[k][l];
          };
        };
      };

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real den0=fwd_nuc[k][0][m].get_dat(i);

            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];

            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc[k][l][m].get_dat(i);
              // --- Number density term
              // (dM)*(flx*1e-24) = (dsig_j*phi_j/flx)*(flx*1e-24) = dsig_j*phi_j*1e-24
              real dsig=xs*(flx*power_factor[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor[k][l]*vol_med[m];
                sum-=pow_adj[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff[k];
	      sum-=tmp*chi_gpt_fwdflx[k][m].get_dat(j);
	    };
	  };
        };
        sns1d.put_data(j,sum/response);
      };
      sns.PutSensitivity1D(matnum,mt,sns1d);
    };

  }; // end of nuclide loop for cross section sensitivities

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
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
          for(int m=0;m<mednum_fuel;m++){
            real rr=xsf_1g[k][m][pos0]*rat;
            for(int l=0;l<sub_step;l++){
              val+=(adj_nuc[k][l][m].get_dat(pos)*rr*total_flux[k][l][m]*1e-24*fwd_nuc[k][l][m].get_dat(pos0))*delt[k][l];
	    };
	  };
	};
        sns.PutSensitivity0D(id,18000000+idfis,val/response);
      };
    };
  };


  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  // (For decay heat sensitivity, direct term should be taken into account)
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda
      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int m=0;m<mednum_fuel;m++){
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l][m].get_dat(i);
            real val=-decay_c*den*adj_nuc[k][l][m].get_dat(i);
            int tmp=bu.GetBC().GetNdivDecay(matnum);
            for(int j=0;j<tmp;j++){
              int id2=bu.GetBC().GetNextIDDecay(matnum,j);
              int pos=bu.SearchNuclide(id2);
              if(pos!=-1){
   	        real rat=bu.GetBC().GetRatioDecay(matnum,j);
	        val+=rat*decay_c*den*adj_nuc[k][l][m].get_dat(pos);
	      };
	    };
	    val*=delt[k][l];
	    sum+=val;
	  };
	};
      };
      sum*=100.;// because dT=0.01T
      sns.PutSensitivity0D(matnum,8888,sum/response);
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
            for(int m=0;m<mednum_fuel;m++){
	      for(int l=0;l<sub_step;l++){
	        real den=fwd_nuc[k][l][m].get_dat(i);
	        real val=0.;
	        int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	        int pos=bu.SearchNuclide(id2);
	        if(pos!=-1){
	  	  val+=adj_nuc[k][l][m].get_dat(pos)*decay_c*rat*den;
	        };
	        time+=delt[k][l];
	        val*=delt[k][l];
	        sum+=val;
	      };
	    };
	  };
	  sum*=100.;// because dr=0.01r
	  sum=sum/response;
	  sns_tmp1[j]=sum;
	};

	// (Make constrained sensitivity)
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

  string filename2=filename+"."+midt.Name(nucid_target);
  if(num_med==1){
    filename2=filename2+".med"+IntToString(med_target[0]);
  }else if(num_med==mednum_fuel){
  }else{
    filename2=filename2+".mednum"+IntToString(num_med);
  };
  //sns.WriteFile("./",filename+"."+midt.Name(nucid_target)+".med"+IntToString(med_target));
  sns.WriteFile("./",filename2);

  }; 
};

void MulticellBurner::SensitivityCalculationPC(Burnup &bu,int med_normalize, int num_med, int *med_target,int num,int *nucid_target_list,string filename)
{
  //corrector_calc=true;
  bool isotropic_approx=false; // Isotropic approximation in perturbation integral

  int ssv=40;  // (sub-sub-time step division)

  wc_gpt_wpc=0.6; // !!!

  // +++

  vector<bool> med_target_bool(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    med_target_bool[i]=false;
  };
  for(int i=0;i<num_med;i++){
    if(med_target[i]>=mednum_fuel){
      cout<<"# Error in MulticellBurner::SensitivityCalculation(PC).\n";
      cout<<"# Medium index "<<med_target[i]<<" is inappopriate.\n";
      exit(0);
    };
    med_target_bool[med_target[i]]=true;
  };

  // +++

  PreCalculation(bu);

  // +++ Array setting for adjoint calculation +++
  pow_adj.resize(burn_step);
  adj_nuc.resize(burn_step);
  macxs.resize(burn_step);

  pow_adj_c.resize(burn_step); // !!!
  adj_nuc_c.resize(burn_step); // !!!
  macxs_p.resize(burn_step); // !!!

  bilinear_flx.resize(burn_step);
  bilinear_flx_c.resize(burn_step); // !!!
  bilinear_aflx.resize(burn_step);
  bilinear_aflx_c.resize(burn_step); // !!!
  chi_gpt_fwdflx.resize(burn_step);
  chi_gpt_fwdflx_c.resize(burn_step); // !!!
  keff_p.resize(burn_step+1);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];

    pow_adj[i].resize(sub_step);
    pow_adj_c[i].resize(sub_step); // !!!
    adj_nuc[i].resize(sub_step);
    adj_nuc_c[i].resize(sub_step); // !!!
    for(int j=0;j<sub_step;j++){
      adj_nuc[i][j].resize(mednum_fuel);
      adj_nuc_c[i][j].resize(mednum_fuel); // !!!
      for(int k=0;k<mednum_fuel;k++){
        adj_nuc[i][j][k].put_imax(nucn); 
        adj_nuc_c[i][j][k].put_imax(nucn); // !!! 
      };
    };

    bilinear_flx[i].resize(mednum_fuel);
    bilinear_flx_c[i].resize(mednum_fuel); // !!!
    bilinear_aflx[i].resize(mednum_fuel);
    bilinear_aflx_c[i].resize(mednum_fuel); // !!!
    chi_gpt_fwdflx[i].resize(mednum_fuel);
    chi_gpt_fwdflx_c[i].resize(mednum_fuel); // !!!
    for(int k=0;k<mednum_fuel;k++){
      bilinear_flx[i][k].put_imax(group);
      bilinear_flx_c[i][k].put_imax(group); // !!!
      bilinear_aflx[i][k].put_imax(group);
      bilinear_aflx_c[i][k].put_imax(group); // !!!
      chi_gpt_fwdflx[i][k].put_imax(group);
      chi_gpt_fwdflx_c[i][k].put_imax(group); // !!!
    };

    macxs[i].resize(mednum_fuel);
    macxs_p[i].resize(mednum_fuel); // !!!
    for(int j=0;j<mednum_fuel;j++){
      macxs[i][j].Init("MacroCrossSection");
      macxs_p[i][j].Init("MacroCrossSection"); // !!!
      macxs[i][j].PutGrp(107);
      macxs_p[i][j].PutGrp(107);
    };
  };

  // +++ Forward burnup calculation +++

  ForwardCalculation(bu,med_normalize,true);

  // +++ Integrating forward number density

  IntegratingForwardNumberDensity(fwd_nuc,fwd_nuc_int);
  IntegratingForwardNumberDensity(fwd_nuc_p,fwd_nuc_p_int);

  // +++ Adjoint burnup calculation +++

  GeneralOption opta;
  opta.PutAdjointCal();

  for(int ii=0;ii<num;ii++){

    int nucid_target=nucid_target_list[ii];
    int nucid_turn=med[0].SearchNuclide(nucid_target);
    if(nucid_turn==-1){
      cout<<"# Error !\n";
      exit(0);
    };

    //real response=fwd_nuc[burn_step][0][med_target].get_dat(nucid_turn);
    real response=0.;
    real volsum=0.;
    for(int jj=0;jj<mednum_fuel;jj++){
      if(med_target_bool[jj]){
        real nd=fwd_nuc[burn_step][0][jj].get_dat(nucid_turn);
	response+=nd*vol_med[jj];
	volsum+=vol_med[jj];
      };        
    };
    response/=volsum;

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response ("<<midt.Name(nucid_target)<<" ND) : "<<response<<"\n";

  // +++ Adjoint calculation

  MECSystem lat_gpt(group,mednum);

  lat_gpt.NoPrint();
  lat_gpt.PutTrajectorySet(&sys_f);
  int totsn=lat_gpt.GetQuad().GetSN();

  for(int i=0;i<mednum_fuel;i++){
    lat_gpt.AddMedium(med[0]);
    lat_gpt.GetMedium(i).NuclideClear();
  };
  for(int jj=0;jj<mednum_nonfuel;jj++){
    lat_gpt.AddMedium(med[1+jj]);
  };
  lat_gpt.PutRegMed(region_medium);
  lat_gpt.PutGeneralOption(opta);
  lat_gpt.PutPL(0);

  lat_gpt.NoCMRAcceleration();
  //lat_gpt.NoTransportApprox();
  lat_gpt.PutWriteFlux();
  lat_gpt.SetArray();

  
  vector<GroupData1D> gpt_flx(totm);
  vector< vector<GroupData1D> > gpt_aflx(totm);
  for(int i=0;i<totm;i++){
    gpt_flx[i].put_imax(group);
    gpt_aflx[i].resize(group);
    for(int g=0;g<group;g++){
      gpt_aflx[i][g].put_imax(totsn);
    };
  };

  vector<GroupData1D> gpt_src(mednum_fuel);
  vector<GroupData1D> gpt_src2(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    gpt_src[i].put_imax(group);
    gpt_src2[i].put_imax(group);
  };

  vector<GroupData1D> adj_nuc_e(mednum_fuel);
  vector<GroupData1D> adj_nuc_e_c(mednum_fuel); // for corrector calc

  // +++ Final condition for adjoint number density
  //real ev_to_j=1.60219e-19;
  //real total_decay_energy=0.;
  for(int i=0;i<mednum_fuel;i++){
    adj_nuc_e[i].put_imax(nucn);
    adj_nuc_e[i].set_zero();
    for(int j=0;j<nucn;j++){
      int id=med[0].GetNuclideInTurn(j).GetID();
      /*
      if(id==nucid_target&&i==med_target){
	adj_nuc_e[i].put_data(j,1.0);
      };
      */
      if(id==nucid_target&&med_target_bool[i]){
	adj_nuc_e[i].put_data(j,vol_med[i]/volsum);
      };

    };
  };

  GroupData2D trmat_flxindep=bu.GetTrmatFlxInDep();
  for(int st=burn_step-1;st>=0;st--){

    real power_density=power_density_list[st];
    int sub_step=sub_step_list[st];

    cout<<"#   Adjoint calculation step : "<<st<<"\n";

    adj_nuc_e_c=adj_nuc_e; // to store initial adjoint for following predictor calculation

    // +++ CORRECTOR PART ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real keff_step=keff[st];
    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross section data
      // is updated in the preceding lines in MEC
    };

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g[st][m][i],xsc_1g[st][m][i],xsn2n_1g[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc[st][j][m]+fwd_nuc[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        //pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
        pow_adj[st][j]/=flux_level_list[st]*volsum;
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);
  
        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor[st][j];
		};
	      };
	    };
	  };
        }else{
	  //if(m==med_target){
	  if(med_target_bool[m]){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc[st][j][m]))*delt[st][j]*power_factor[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    // ++++
    int sntot=lat_gpt.GetQuad().GetSN();
    int pdiv=lat_gpt.GetPolarAngleDivision();
    int sn_per_pang=sntot/pdiv;

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh[st][mm]));
          real tmp=gpt_flx[mm]*macxs[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
	      real tmp;
	      if(aflx_legendre==-1){
                tmp=omega*volaflx_mesh[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
	      }else{
                real tmp2=0.;
                int pln=0;
  	        for(int l=0;l<=aflx_legendre;l++){
		  for(int m=0;m<=l;m++){
		    tmp2+=(2.*l+1.)/PI4*quad.GetMoment(pln,sn)*volaflx_pl[st][mm][g].get_dat(pln);
		    pln++;
		  };
	        };
                tmp=omega*tmp2*gpt_aflx[mm][g].get_dat(sn2);
	      };
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // +++ CORRECTOR CALCULATION END ++++++++++++++++++++++++++++++++++++++++

    swap(adj_nuc_e, adj_nuc_e_c);
    swap(pow_adj[st], pow_adj_c[st]);
    swap(adj_nuc[st], adj_nuc_c[st]);
    swap(bilinear_flx[st], bilinear_flx_c[st]);
    swap(bilinear_aflx[st], bilinear_aflx_c[st]);
    swap(chi_gpt_fwdflx[st], chi_gpt_fwdflx_c[st]);
    
    // +++ PREDICTOR CALCULATION +++++++++++++++++++++++++++++++++++++++++++++

    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs_p[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross sectio data
      // is updated in the preceding lines in MEC
    };

    // jump condition for generalized adjoint at the end of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
              tmp2+=xsa*bilinear_flx_c[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx_c[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step; // keff should be one at corrector step
          };
          adj_nuc_e[m].add_data(k,(tmp2-tmp3)*wc_gpt_wpc/(1.-wc_gpt_wpc));
	};
      };
    };

    keff_step=keff_p[st];

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g_p[st][m][i],xsc_1g_p[st][m][i],xsn2n_1g_p[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux_p[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc_p[st][j][m]+fwd_nuc_p[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        //pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
        pow_adj[st][j]/=flux_level_list[st]*volsum;
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g_p[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux_p[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);

        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc_p[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor_p[st][j];
		};
	      };
	    };
	  };
        }else{
	  //if(m==med_target){
	  if(med_target_bool[m]){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor_p[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc_p[st][j][m]))*delt[st][j]*power_factor_p[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh_p[st][mm]));
          real tmp=gpt_flx[mm]*macxs_p[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh_p[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
	      real tmp;
	      if(aflx_legendre==-1){
                tmp=omega*volaflx_mesh_p[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
	      }else{
                real tmp2=0.;
                int pln=0;
	        for(int l=0;l<=aflx_legendre;l++){
		  for(int m=0;m<=l;m++){
		    tmp2+=(2.*l+1.)/PI4*quad.GetMoment(pln,sn)*volaflx_pl_p[st][mm][g].get_dat(pln);
		    pln++;
		  };
	        };
                tmp=omega*tmp2*gpt_aflx[mm][g].get_dat(sn2);
	      };
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // jump condition for generalized adjoint at the beginning of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption or total term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
	      tmp2+=xsa*bilinear_flx[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh_p[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs_p[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step;
          };
          adj_nuc_e[m].add_data(k,tmp2-tmp3);
	};

      };
    };

    // +++ PREDICTOR CALCULATION END +++++++++++++++++++++++++++++++++++++++++

    // Averaging of adjoint number density
    for(int m=0;m<mednum_fuel;m++){
      adj_nuc_e[m]=adj_nuc_e[m]*(1.-wc_gpt_wpc)+adj_nuc_e_c[m]*wc_gpt_wpc;
    };

  }; // the end ot step

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Sensitivity Printing

  vector< vector< vector<real> > > nadj_dm_nfwd;

  SensitivityData sns;
  sns.PutName("dummy","dummy","dummy");
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

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

      // +++ predictor
      GroupData1D sns1d_p(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc_p,adj_nuc,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh_p[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc_p[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor_p[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor_p[k][l]*vol_med[m];
                sum-=pow_adj[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][0][m].get_dat(i);
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff_p[k];
	      sum-=tmp*chi_gpt_fwdflx[k][m].get_dat(j);
	    };
	  };
        };
        sns1d_p.put_data(j,sum/response);
      };

      // +++ corrector
      GroupData1D sns1d_c(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc,adj_nuc_c,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor[k][l]*vol_med[m];
                sum-=pow_adj_c[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][sub_step][m].get_dat(i); // !!! CAUTION
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx_c[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx_c[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff[k];
	      sum-=tmp*chi_gpt_fwdflx_c[k][m].get_dat(j);
	    };
	  };
	};
        sns1d_c.put_data(j,sum/response);
      };

      sns1d=sns1d_p*(1.-wc_gpt_wpc)+sns1d_c*wc_gpt_wpc;
      sns.PutSensitivity1D(matnum,mt,sns1d);

    };
  }; // end of nuclide loop for cross section sensitivities

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
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

        vector<real> snsval(2,0.);
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
        for(int jj=0;jj<2;jj++){ // predictor & corrector

          real val=0.;
          for(int k=0;k<burn_step;k++){
            int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
              real rr=xsf_1g[k][m][pos0]*rat;
              for(int l=0;l<sub_step;l++){
                val+=(adj_nuc[k][l][m].get_dat(pos)*rr*total_flux[k][l][m]*1e-24*fwd_nuc[k][l][m].get_dat(pos0))*delt[k][l];
	      };
	    };
	  };
  	  snsval[jj]=val/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

        }; // end of loop-jj
        swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

        real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
        sns.PutSensitivity0D(id,18000000+idfis,snstot);
        //sns.PutSensitivity0D(id,18000000+idfis,val/response);

      };
    };
  };

  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  // (For decay heat sensitivity, direct term should be taken into account)
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda

      vector<real> snsval(2,0.);
      swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
      for(int jj=0;jj<2;jj++){ // predictor & corrector

      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int m=0;m<mednum_fuel;m++){
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l][m].get_dat(i);
            real val=-decay_c*den*adj_nuc[k][l][m].get_dat(i);
            int tmp=bu.GetBC().GetNdivDecay(matnum);
            for(int j=0;j<tmp;j++){
              int id2=bu.GetBC().GetNextIDDecay(matnum,j);
              int pos=bu.SearchNuclide(id2);
              if(pos!=-1){
   	        real rat=bu.GetBC().GetRatioDecay(matnum,j);
	        val+=rat*decay_c*den*adj_nuc[k][l][m].get_dat(pos);
	      };
	    };
	    val*=delt[k][l];
	    sum+=val;
	  };
	};
      };
      sum*=100.;// because dT=0.01T

      snsval[jj]=sum/response;

      if(jj==0){
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
        swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
      };

      }; // end of loop-jj
      swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

      real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
      sns.PutSensitivity0D(matnum,8888,snstot);
      //sns.PutSensitivity0D(matnum,8888,sum/response);
    };
  };

  // +++ For Branching Ratio +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(matnum);
      if(channel>1){

	vector<real> ratio;
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  ratio[j]=bu.GetBC().GetRatioDecay(matnum,j);
	};
	vector<real> sns_tmp1;
	sns_tmp1.resize(channel);

	for(int j=0;j<channel;j++){

          vector<real> snsval(2,0.);
          swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
          for(int jj=0;jj<2;jj++){ // predictor & corrector

	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(matnum,j);
	  rat*=0.01;
	  real time=0.;
	  for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
	      for(int l=0;l<sub_step;l++){
	        real den=fwd_nuc[k][l][m].get_dat(i);
	        real val=0.;
	        int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	        int pos=bu.SearchNuclide(id2);
	        if(pos!=-1){
	  	  val+=adj_nuc[k][l][m].get_dat(pos)*decay_c*rat*den;
	        };
	        time+=delt[k][l];
	        val*=delt[k][l];
	        sum+=val;
	      };
	    };
	  };
	  sum*=100.;// because dr=0.01r
          snsval[jj]=sum/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

          }; // end of loop-jj
          swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c
          real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
	  sns_tmp1[j]=snstot;

	};

	// (Make constrained sensitivity)
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

 

  //

  string filename2=filename+"."+midt.Name(nucid_target);
  if(num_med==1){
    filename2=filename2+".med"+IntToString(med_target[0]);
  }else if(num_med==mednum_fuel){
  }else{
    filename2=filename2+".mednum"+IntToString(num_med);
  };
  //sns.WriteFile("./",filename+"."+midt.Name(nucid_target)+".med"+IntToString(med_target));
  sns.WriteFile("./",filename2);

  }; // end of reactor physics parameter loop for sensitivity calculation

  
};


void MulticellBurner::SensitivityCalculationKeffEOC(Burnup &bu,int med_normalize, string filename)
{
  if(!corrector_calc){
    cout<<"# Error in MulticellBurner::SensitivityCalculationKeffEOC.\n";
    cout<<"# Predictor-corrector method should be used.\n";
    exit(0);
  };

  bool isotropic_approx=false; // Isotropic approximation in perturbation integral

  int ssv=40;  // (sub-sub-time step division)

  wc_gpt_wpc=0.6; // !!!

  PreCalculation(bu);

  // +++ Array setting for adjoint calculation +++
  pow_adj.resize(burn_step);
  adj_nuc.resize(burn_step);
  macxs.resize(burn_step);

  pow_adj_c.resize(burn_step); // !!!
  adj_nuc_c.resize(burn_step); // !!!
  macxs_p.resize(burn_step); // !!!

  bilinear_flx.resize(burn_step);
  bilinear_flx_c.resize(burn_step); // !!!
  bilinear_aflx.resize(burn_step);
  bilinear_aflx_c.resize(burn_step); // !!!
  chi_gpt_fwdflx.resize(burn_step);
  chi_gpt_fwdflx_c.resize(burn_step); // !!!
  keff_p.resize(burn_step+1);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];

    pow_adj[i].resize(sub_step);
    pow_adj_c[i].resize(sub_step); // !!!
    adj_nuc[i].resize(sub_step);
    adj_nuc_c[i].resize(sub_step); // !!!
    for(int j=0;j<sub_step;j++){
      adj_nuc[i][j].resize(mednum_fuel);
      adj_nuc_c[i][j].resize(mednum_fuel); // !!!
      for(int k=0;k<mednum_fuel;k++){
        adj_nuc[i][j][k].put_imax(nucn); 
        adj_nuc_c[i][j][k].put_imax(nucn); // !!! 
      };
    };

    bilinear_flx[i].resize(mednum_fuel);
    bilinear_flx_c[i].resize(mednum_fuel); // !!!
    bilinear_aflx[i].resize(mednum_fuel);
    bilinear_aflx_c[i].resize(mednum_fuel); // !!!
    chi_gpt_fwdflx[i].resize(mednum_fuel);
    chi_gpt_fwdflx_c[i].resize(mednum_fuel); // !!!
    for(int k=0;k<mednum_fuel;k++){
      bilinear_flx[i][k].put_imax(group);
      bilinear_flx_c[i][k].put_imax(group); // !!!
      bilinear_aflx[i][k].put_imax(group);
      bilinear_aflx_c[i][k].put_imax(group); // !!!
      chi_gpt_fwdflx[i][k].put_imax(group);
      chi_gpt_fwdflx_c[i][k].put_imax(group); // !!!
    };

    macxs[i].resize(mednum_fuel);
    macxs_p[i].resize(mednum_fuel); // !!!
    for(int j=0;j<mednum_fuel;j++){
      macxs[i][j].Init("MacroCrossSection");
      macxs_p[i][j].Init("MacroCrossSection"); // !!!
      macxs[i][j].PutGrp(107);
      macxs_p[i][j].PutGrp(107);
    };
  };

  // +++ Forward burnup calculation +++

  ForwardCalculation(bu,med_normalize,true);

  // +++ Integrating forward number density

  IntegratingForwardNumberDensity(fwd_nuc,fwd_nuc_int);
  IntegratingForwardNumberDensity(fwd_nuc_p,fwd_nuc_p_int);

  // +++ Adjoint burnup calculation +++

  GeneralOption opta;
  opta.PutAdjointCal();

  // !!!

  // +++ EOC calculation
  GeneralOption opt;


  /*
  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  vector<Medium> med_tmp(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[burn_step][0][i],0);
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    med[0].CalMacroFromMicro();
    med_tmp[i]=med[0];
  };


  MECSystem lata(group,mednum);
  //lata.NoPrint();
  lata.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum_fuel;i++){
    lata.AddMedium(med_tmp[i]);
  };
  lata.AddMedium(med[med_clad]);
  lata.AddMedium(med[med_water]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutPL(0);
  lata.NoCMRAcceleration();
  lata.PutWriteFlux();
  real k_adj=lata.CalIgen();

  MECSystem lat(group,mednum);
  //lata.NoPrint();
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum_fuel;i++){
    lat.AddMedium(med_tmp[i]);
  };
  lat.AddMedium(med[med_clad]);
  lat.AddMedium(med[med_water]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  lat.PutWriteFlux();
  real k_fwd=lat.CalIgen();

  // (Direct term)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sns_dir=lata.CalSensitivityNew(&lat,k_fwd,nucnum,nucid);
  delete [] nucid;

  sns_dir.PutName("mburner","k_EOC","unknown");
  sns_dir.WriteFile("./","sns.k_EOC_dir");

  real response=k_fwd;

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response (k-inf) : "<<response<<"\n";
  */

  // +++ Adjoint calculation
  // (MEC)
  MECSystem lat_gpt(group,mednum);

  lat_gpt.NoPrint();
  lat_gpt.PutTrajectorySet(&sys_f);
  int totsn=lat_gpt.GetQuad().GetSN();

  for(int i=0;i<mednum_fuel;i++){
    lat_gpt.AddMedium(med[0]);
    lat_gpt.GetMedium(i).NuclideClear();
  };
  for(int jj=0;jj<mednum_nonfuel;jj++){
    lat_gpt.AddMedium(med[1+jj]);
  };
  //lat_gpt.AddMedium(med[med_clad]);
  //lat_gpt.AddMedium(med[med_water]);
  lat_gpt.PutRegMed(region_medium);
  lat_gpt.PutGeneralOption(opta);
  lat_gpt.PutPL(0);

  lat_gpt.NoCMRAcceleration();
  //lat_gpt.NoTransportApprox();
  lat_gpt.PutWriteFlux();
  lat_gpt.SetArray();

  
  vector<GroupData1D> gpt_flx(totm);
  vector< vector<GroupData1D> > gpt_aflx(totm);
  for(int i=0;i<totm;i++){
    gpt_flx[i].put_imax(group);
    gpt_aflx[i].resize(group);
    for(int g=0;g<group;g++){
      gpt_aflx[i][g].put_imax(totsn);
    };
  };

  vector<GroupData1D> gpt_src(mednum_fuel);
  vector<GroupData1D> gpt_src2(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    gpt_src[i].put_imax(group);
    gpt_src2[i].put_imax(group);
  };

  vector<GroupData1D> adj_nuc_e(mednum_fuel);
  vector<GroupData1D> adj_nuc_e_c(mednum_fuel); // for corrector calc

  //


  real response;
  SensitivityData sns_dir=CalInitialAdjNDKeffEOC(response, adj_nuc_e);

  /*
  // +++ Final condition for adjoint number density
  for(int i=0;i<mednum_fuel;i++){
    adj_nuc_e[i].put_imax(nucn);
    adj_nuc_e[i].set_zero();
    for(int j=0;j<nucn;j++){
      if(med_tmp[i].GetNuclideInTurn(j).GetGrp()!=-1){
        real org=med_tmp[i].GetNuclideInTurn(j).GetDensity();
        if(org>1e-20){
          lat.GetMedium(i).CalMacroFromMicro();
          lata.GetMedium(i).CalMacroFromMicro();
          lat.GetMedium(i).GetNuclideInTurn(j).PutDensity(org*1.01);
          lat.GetMedium(i).CalMacroFromMicro();
          real dk=lata.CalReactivity(&lat,k_adj,k_fwd,false)*k_adj*k_fwd;
          lat.GetMedium(i)=med_tmp[i];
          //int mmt=med[0].GetNuclideInTurn(j).GetMatnum();
          //cout<<mmt<<" "<<midt.Name(mmt)<<" "<<dk*100.<<"\n";
          real src=0.;
          if(org!=0.)src=dk/(org*0.01);
          adj_nuc_e[i].put_data(j,src);
	};
      };
    };
    med_tmp[i].NuclideClear();
  };
  */

  GroupData2D trmat_flxindep=bu.GetTrmatFlxInDep();
  for(int st=burn_step-1;st>=0;st--){

    real power_density=power_density_list[st];
    int sub_step=sub_step_list[st];

    cout<<"#   Adjoint calculation step : "<<st<<"\n";

    adj_nuc_e_c=adj_nuc_e; // to store initial adjoint for following predictor calculation

    // +++ CORRECTOR PART ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real keff_step=keff[st];
    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross section data
      // is updated in the preceding lines in MEC
    };

    // (PJI)
    /*
    if(pij_storing){
      for(int g=0;g<group;g++){
        lat_gpt.GetPij(g).copy(pij_store[st][1][g]);
      };
    }else{
      lat_gpt.PutPij();
    };
    */
    // (MEC)
    // no processing    

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g[st][m][i],xsc_1g[st][m][i],xsn2n_1g[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc[st][j][m]+fwd_nuc[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_normalize];
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);
  
        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor[st][j];
		};
	      };
	    };
	  };
        }else{
	  if(m==med_normalize){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc[st][j][m]))*delt[st][j]*power_factor[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    // ++++
    int sntot=lat_gpt.GetQuad().GetSN();
    int pdiv=lat_gpt.GetPolarAngleDivision();
    int sn_per_pang=sntot/pdiv;

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh[st][mm]));
          real tmp=gpt_flx[mm]*macxs[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
	      real tmp;
	      if(aflx_legendre==-1){
                tmp=omega*volaflx_mesh[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
	      }else{
                real tmp2=0.;
                int pln=0;
  	        for(int l=0;l<=aflx_legendre;l++){
		  for(int m=0;m<=l;m++){
		    tmp2+=(2.*l+1.)/PI4*quad.GetMoment(pln,sn)*volaflx_pl[st][mm][g].get_dat(pln);
		    pln++;
		  };
	        };
                tmp=omega*tmp2*gpt_aflx[mm][g].get_dat(sn2);
	      };
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // +++ CORRECTOR CALCULATION END ++++++++++++++++++++++++++++++++++++++++

    swap(adj_nuc_e, adj_nuc_e_c);
    swap(pow_adj[st], pow_adj_c[st]);
    swap(adj_nuc[st], adj_nuc_c[st]);
    swap(bilinear_flx[st], bilinear_flx_c[st]);
    swap(bilinear_aflx[st], bilinear_aflx_c[st]);
    swap(chi_gpt_fwdflx[st], chi_gpt_fwdflx_c[st]);
    
    // +++ PREDICTOR CALCULATION +++++++++++++++++++++++++++++++++++++++++++++

    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs_p[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross sectio data
      // is updated in the preceding lines in MEC
    };

    // (PJI)
    /*
    if(pij_storing){
      for(int g=0;g<group;g++){
        lat_gpt.GetPij(g).copy(pij_store[st][0][g]);
      };
    }else{
      lat_gpt.PutPij();
    };
    */
    // (MEC)
    // no processing
    
    // jump condition for generalized adjoint at the end of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
              tmp2+=xsa*bilinear_flx_c[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx_c[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step; // keff should be one at corrector step
          };
          adj_nuc_e[m].add_data(k,(tmp2-tmp3)*wc_gpt_wpc/(1.-wc_gpt_wpc));
	};
      };
    };

    keff_step=keff_p[st];

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g_p[st][m][i],xsc_1g_p[st][m][i],xsn2n_1g_p[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux_p[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc_p[st][j][m]+fwd_nuc_p[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_normalize];
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g_p[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux_p[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);

        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc_p[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor_p[st][j];
		};
	      };
	    };
	  };
        }else{
	  if(m==med_normalize){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor_p[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc_p[st][j][m]))*delt[st][j]*power_factor_p[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh_p[st][mm]));
          real tmp=gpt_flx[mm]*macxs_p[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh_p[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
	      if(aflx_legendre==-1){
                tmp=omega*volaflx_mesh_p[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
	      }else{
                real tmp2=0.;
                int pln=0;
	        for(int l=0;l<=aflx_legendre;l++){
		  for(int m=0;m<=l;m++){
		    tmp2+=(2.*l+1.)/PI4*quad.GetMoment(pln,sn)*volaflx_pl_p[st][mm][g].get_dat(pln);
		    pln++;
		  };
	        };
                tmp=omega*tmp2*gpt_aflx[mm][g].get_dat(sn2);
	      };
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // jump condition for generalized adjoint at the beginning of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption or total term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
	      tmp2+=xsa*bilinear_flx[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh_p[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs_p[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step;
          };
          adj_nuc_e[m].add_data(k,tmp2-tmp3);
	};

      };
    };

    // +++ PREDICTOR CALCULATION END +++++++++++++++++++++++++++++++++++++++++

    // Averaging of adjoint number density
    for(int m=0;m<mednum_fuel;m++){
      adj_nuc_e[m]=adj_nuc_e[m]*(1.-wc_gpt_wpc)+adj_nuc_e_c[m]*wc_gpt_wpc;
    };

  }; // the end ot step

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Sensitivity Printing

  vector< vector< vector<real> > > nadj_dm_nfwd;

  SensitivityData sns;
  sns.PutName("dummy","dummy","dummy");
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

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

      // +++ predictor
      GroupData1D sns1d_p(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc_p,adj_nuc,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh_p[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc_p[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor_p[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor_p[k][l]*vol_med[m];
                sum-=pow_adj[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][0][m].get_dat(i);
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff_p[k];
	      sum-=tmp*chi_gpt_fwdflx[k][m].get_dat(j);
	    };
	  };
        };
        sns1d_p.put_data(j,sum/response);
      };

      // +++ corrector
      GroupData1D sns1d_c(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc,adj_nuc_c,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor[k][l]*vol_med[m];
                sum-=pow_adj_c[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][sub_step][m].get_dat(i); // !!! CAUTION
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx_c[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx_c[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff[k];
	      sum-=tmp*chi_gpt_fwdflx_c[k][m].get_dat(j);
	    };
	  };
	};
        sns1d_c.put_data(j,sum/response);
      };

      sns1d=sns1d_p*(1.-wc_gpt_wpc)+sns1d_c*wc_gpt_wpc;
      sns.PutSensitivity1D(matnum,mt,sns1d);

    };
  }; // end of nuclide loop for cross section sensitivities

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
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

        vector<real> snsval(2,0.);
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
        for(int jj=0;jj<2;jj++){ // predictor & corrector

          real val=0.;
          for(int k=0;k<burn_step;k++){
            int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
              real rr=xsf_1g[k][m][pos0]*rat;
              for(int l=0;l<sub_step;l++){
                val+=(adj_nuc[k][l][m].get_dat(pos)*rr*total_flux[k][l][m]*1e-24*fwd_nuc[k][l][m].get_dat(pos0))*delt[k][l];
	      };
	    };
	  };
  	  snsval[jj]=val/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

        }; // end of loop-jj
        swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

        real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
        sns.PutSensitivity0D(id,18000000+idfis,snstot);
        //sns.PutSensitivity0D(id,18000000+idfis,val/response);

      };
    };
  };

  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  // (For decay heat sensitivity, direct term should be taken into account)
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda

      vector<real> snsval(2,0.);
      swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
      for(int jj=0;jj<2;jj++){ // predictor & corrector

      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int m=0;m<mednum_fuel;m++){
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l][m].get_dat(i);
            real val=-decay_c*den*adj_nuc[k][l][m].get_dat(i);
            int tmp=bu.GetBC().GetNdivDecay(matnum);
            for(int j=0;j<tmp;j++){
              int id2=bu.GetBC().GetNextIDDecay(matnum,j);
              int pos=bu.SearchNuclide(id2);
              if(pos!=-1){
   	        real rat=bu.GetBC().GetRatioDecay(matnum,j);
	        val+=rat*decay_c*den*adj_nuc[k][l][m].get_dat(pos);
	      };
	    };
	    val*=delt[k][l];
	    sum+=val;
	  };
	};
      };
      sum*=100.;// because dT=0.01T

      snsval[jj]=sum/response;

      if(jj==0){
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
        swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
      };

      }; // end of loop-jj
      swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

      real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
      sns.PutSensitivity0D(matnum,8888,snstot);
      //sns.PutSensitivity0D(matnum,8888,sum/response);
    };
  };

  // +++ For Branching Ratio +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(matnum);
      if(channel>1){

	vector<real> ratio;
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  ratio[j]=bu.GetBC().GetRatioDecay(matnum,j);
	};
	vector<real> sns_tmp1;
	sns_tmp1.resize(channel);

	for(int j=0;j<channel;j++){

          vector<real> snsval(2,0.);
          swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
          for(int jj=0;jj<2;jj++){ // predictor & corrector

	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(matnum,j);
	  rat*=0.01;
	  real time=0.;
	  for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
	      for(int l=0;l<sub_step;l++){
	        real den=fwd_nuc[k][l][m].get_dat(i);
	        real val=0.;
	        int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	        int pos=bu.SearchNuclide(id2);
	        if(pos!=-1){
	  	  val+=adj_nuc[k][l][m].get_dat(pos)*decay_c*rat*den;
	        };
	        time+=delt[k][l];
	        val*=delt[k][l];
	        sum+=val;
	      };
	    };
	  };
	  sum*=100.;// because dr=0.01r
          snsval[jj]=sum/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

          }; // end of loop-jj
          swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c
          real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
	  sns_tmp1[j]=snstot;

	};

	// (Make constrained sensitivity)
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

  //

  sns.WriteFile("./",filename+"_indir");
  sns.AddSensitivityData(sns_dir);
  sns.WriteFile("./",filename);
  
};

SensitivityData MulticellBurner::CalInitialAdjNDKeffEOC(real &response, vector<GroupData1D> &adj_nuc_e)
{
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  /*
  vector<Medium> med_tmp(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[burn_step][0][i],0);
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    med[0].CalMacroFromMicro();
    med_tmp[i]=med[0];
  };
  */

  MECSystem lata(group,mednum);
  //lata.NoPrint();
  lata.PutTrajectorySet(&sys_f);

  MECSystem lat(group,mednum);
  //lata.NoPrint();
  lat.PutTrajectorySet(&sys_f);

  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[burn_step][0][i],0);
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    opc.CalThermalScatteringMatrix(med[0],xslib,4.048);    
    med[0].CalMacroFromMicro();
    lata.AddMedium(med[0]);
    lata.GetMedium(i).NuclideClear();
    lat.AddMedium(med[0]);
  };
  for(int jj=0;jj<mednum_nonfuel;jj++){
    lata.AddMedium(med[1+jj]);
  };
  //lata.AddMedium(med[med_clad]);
  //lata.AddMedium(med[med_water]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutPL(0);
  lata.NoCMRAcceleration();
  lata.PutWriteFlux();
  real k_adj=lata.CalIgen();

  for(int jj=0;jj<mednum_nonfuel;jj++){
    lat.AddMedium(med[1+jj]);
  };
  //lat.AddMedium(med[med_clad]);
  //lat.AddMedium(med[med_water]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  lat.PutWriteFlux();
  real k_fwd=lat.CalIgen();

  // (Direct term)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sns_dir=lata.CalSensitivityNew(&lat,k_fwd,nucnum,nucid);
  delete [] nucid;

  /*
  for(int i=0;i<mednum_fuel;i++){
    lata.GetMedium(i).CalMacroFromMicro(); // fission spectrum vector re-calculation?
    lata.GetMedium(i).NuclideClear();
  };
  */

  sns_dir.PutName("mburner","k_EOC","unknown");
  sns_dir.WriteFile("./","sns.k_EOC_dir");

  response=k_fwd;

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response (k-inf) : "<<response<<"\n";


  //vector<GroupData1D> adj_nuc_e(mednum_fuel);

  // +++ Final condition for adjoint number density
  for(int i=0;i<mednum_fuel;i++){
    adj_nuc_e[i].put_imax(nucn);
    adj_nuc_e[i].set_zero();
    for(int j=0;j<nucn;j++){
      if(med[0].GetNuclideInTurn(j).GetGrp()!=-1){
        //real org=med_tmp[i].GetNuclideInTurn(j).GetDensity();
        real org=fwd_nuc[burn_step][0][i].get_dat(j);
        if(org>1e-20){
          //lat.GetMedium(i).CalMacroFromMicro();
          //lata.GetMedium(i).CalMacroFromMicro();
          lat.GetMedium(i).GetNuclideInTurn(j).PutDensity(org*1.01);
          lat.GetMedium(i).CalMacroFromMicro();
	  lat.GetMedium(i).TransportApproximation();
          real dk=lata.CalReactivity(&lat,k_adj,k_fwd,false)*k_adj*k_fwd;
          lat.GetMedium(i).GetNuclideInTurn(j).PutDensity(org);
          //lat.GetMedium(i)=lata.GetMedium(i);
          //int mmt=med[0].GetNuclideInTurn(j).GetMatnum();
          //cout<<mmt<<" "<<midt.Name(mmt)<<" "<<dk*100.<<"\n";
          real src=0.;
          if(org!=0.)src=dk/(org*0.01);
          adj_nuc_e[i].put_data(j,src);
	};
      };
    };
    //med_tmp[i].NuclideClear();
  };

  return sns_dir;

};

void MulticellBurner::IntegratingForwardNumberDensity
(vector< vector< vector<GroupData1D> > > &fwd_nuc, vector< vector< vector<GroupData1D> > > &fwd_nuc_int)
{
  // [fwd_nuc_int] is a number density at time-mesh-center point
  // time-averagted number density is calculated 
  // from those at beginning, center and end of the time mesh

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    for(int j=0;j<sub_step;j++){
      real dt=delt[i][j];
      for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
	real n1=fwd_nuc[i][j][m].get_dat(k);
	real n2=fwd_nuc_int[i][j][m].get_dat(k);
	real n3=fwd_nuc[i][j+1][m].get_dat(k);
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
	fwd_nuc[i][j][m].put_data(k,tmp);
      };
      };
    };
  };
};

void MulticellBurner::CalNadjDMNfwd
(int i, int matnum, int bc_channel, Burnup &bu, vector< vector< vector<GroupData1D> > > &fwd_nuc, vector< vector< vector<GroupData1D> > > &adj_nuc, vector< vector< vector<real> > > &nadj_dm_nfwd)
{
  nadj_dm_nfwd.resize(burn_step);
  for(int k=0;k<burn_step;k++){
    int sub_step=sub_step_list[k];
    nadj_dm_nfwd[k].resize(sub_step);
    for(int l=0;l<sub_step;l++){
      nadj_dm_nfwd[k][l].resize(mednum_fuel);
      for(int m=0;m<mednum_fuel;m++){
        real val=0.;
        val+=-fwd_nuc[k][l][m].get_dat(i)*adj_nuc[k][l][m].get_dat(i);
        int tmp=bu.GetBC().GetNdiv(matnum,bc_channel);
        for(int j=0;j<tmp;j++){
          int id2=bu.GetBC().GetNextID(matnum,bc_channel,j);
          int pos=bu.SearchNuclide(id2);
          if(pos!=-1){
            real rat=bu.GetBC().GetRatio(matnum,bc_channel,j);
            val+=adj_nuc[k][l][m].get_dat(pos)*rat*fwd_nuc[k][l][m].get_dat(i);
          };
            };
            nadj_dm_nfwd[k][l][m]=val*delt[k][l];
          };
        };
      };
};

void MulticellBurner::SensitivityCalculationDirect(Burnup &bu, int med_normalize, int matid, int mt, int target_num, int *medid, int *nucid, string snsname,bool sigt_preserve)
{
  PreCalculation(bu);

  vector<SensitivityData> sns(target_num);

  if(mt!=102&&mt!=18){
    cout<<"# Error in SensitivityCalculationDirect.\n";
    cout<<"# Not yet be coded for MT="<<mt<<"\n";
    exit(0);
  };
  enum xstype xst=sigf;
  if(mt==102)xst=sigc;

  real factor=-0.01;
  //real factor=0.01;

  vector< vector<real> > response(target_num);
  for(int i=0;i<target_num;i++){
    response[i].resize(group+1);
  };

  /*
  // +++ Initial ND storing +++
  vector<int> mat0(nucn);
  vector< vector<real> > den0(mednum_fuel);
  vector<real> temp0(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    temp0[i]=med[i].GetNuclideInTurn(0).GetTemperature();
    den0[i].resize(nucn);
    for(int j=0;j<nucn;j++){
      if(j==0)mat0[i]=med[0].GetNuclideInTurn(i).GetMatnum();
      den0[i][j]=med[i].GetNuclideInTurn(j).GetDensity();
    };
  };
  */

  // +++ Calculation with perturbed cross section
  // (Initial iteration calculation is discarded
  //  because there is an effect only in the initial iteration
  //  and it changes sensitivity calculation results)
  for(int gg=0;gg<group+2;gg++){
    int g=gg-1;
    if(g==-1)g=0;

    real dxs=0.;
    if(g!=group){
      dxs=GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).get_dat(g)*factor;
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).add_data(g,dxs);
      if(sigt_preserve){
        //GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(siginel).add_data(g,-dxs);
        //GetXSLibrary().GetLibData(matid).GetXSData().GetData2d(siginel).add_data(g,g,-dxs);
        GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigel).add_data(g,-dxs);
        GetXSLibrary().GetLibData(matid).GetXSData().GetData2d(sigel).add_data(g,g,-dxs);
      }else{
        GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigt).add_data(g,dxs);
      };
    };

    // +++ Initial selfshielding calculation for non-fuel medium
    med[0].GetFlux().copy(xslib.GetWtflux());
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    opc.GiveInfiniteDillutionCrossSection(med[med_clad],xslib);
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //ForwardCalculation(bu,med_target,false);
    ForwardCalculation(bu,med_normalize,false);
    for(int i=0;i<target_num;i++){
      for(int j=0;j<nucn;j++){
	if(med[0].GetNuclideInTurn(j).GetMatnum()==nucid[i]){
          response[i][g]=fwd_nuc[burn_step][0][medid[i]].get_dat(j);
	};
      };
    };

    if(g!=group){
      GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(xst).add_data(g,-dxs);
      if(sigt_preserve){
        GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigel).add_data(g,dxs);
        GetXSLibrary().GetLibData(matid).GetXSData().GetData2d(sigel).add_data(g,g,dxs);
        //GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(siginel).add_data(g,dxs);
        //GetXSLibrary().GetLibData(matid).GetXSData().GetData2d(siginel).add_data(g,g,dxs);
      }else{
        GetXSLibrary().GetLibData(matid).GetXSData().GetData1d(sigt).add_data(g,-dxs);
      };
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
    sns[i].WriteFile("./",snsname+"."+midt.Name(nucid[i])+".med"+IntToString(medid[i]));
  };

};

void MulticellBurner::ShowCrossSection(int step, int medid)
{
  cout<<"#\n# One-group cross section list (cap/fis/n2n/[sum])\n#\n";
  cout<<"#   Burnup-step : "<<step<<"\n";
  cout<<"#   Medium ID   : "<<medid<<"\n#\n";
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
        int tmp2=8-tmp.size();
	for(int k=0;k<tmp2;k++){cout<<" ";};
	cout<<xsc_1g[step][medid][i]<<"  ";
	cout<<xsf_1g[step][medid][i]<<"  ";
	cout<<xsn2n_1g[step][medid][i]<<"  ";
	cout<<xsc_1g[step][medid][i]+xsf_1g[step][medid][i]+xsn2n_1g[step][medid][i]<<"  ";
	cout<<"\n";
      };
    };
  };
};

void MulticellBurner::ShowCrossSectionNuclideWise(int nucid,bool reactionrate)
{
  int target=-1;
  for(int i=0;i<nucn;i++){
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    if(matno==nucid)target=i;
  };
  if(target==-1){
    cout<<"# Error in MulticellBurner::ShowCrossSectionNuclideReactionWise.\n";
    cout<<"# No nuclide of ID "<<nucid<<" exists.\n";
    exit(0);
  };

  for(int ii=0;ii<2;ii++){

    cout<<"#\n";
    if(!reactionrate){
      cout<<"# One-group cross section [barn]";
    }else{
      cout<<"# Microscopic reaction rate [barn/cm2/s]";
    };
    cout<<" of "<<midt.Name(nucid)<<" : ";
    if(ii==0){
      cout<<" capture";
    }else{
      cout<<" fission";
    };
    cout<<"\n#\n";
    cout<<"# (day)      (GWd/t)    ";
    for(int j=0;j<mednum_fuel;j++){
      cout<<"(Med ";
      WriteOut(j,2);
      cout<<")   ";
    };
    cout<<"\n";

    cout.setf(ios::scientific);
    cout.precision(4);
    for(int i=0;i<burn_step;i++){
      cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
      for(int j=0;j<mednum_fuel;j++){
        real tmp=xsc_1g[i][j][target];
        if(ii==1)tmp=xsf_1g[i][j][target];
	if(reactionrate)tmp*=total_flux[i][0][j];
        cout<<tmp<<" ";
      };
      cout<<"\n";
    };

    cout<<"\n\n";

  };

};

void MulticellBurner::ShowAbsorptionRate(int medid)
{
  cout<<"# Absorption rate results\n";
  cout<<"# Nuclide ID list\n#\n";
  for(int j=0;j<nucn;j++){
    int matid=med[0].GetNuclideInTurn(j).GetMatnum();
    cout<<"# "<<j<<" "<<matid<<" "<<midt.Name(matid)<<"\n";
  };

  if(medid!=-1){

    // Medium-wise absorption rate

    // (pre-calculation of medium-wise contribution to total absorption rate)
    cout<<"# Contribution to total absorption rate of medium "<<medid<<"\n#\n#\n";
    cout<<"# Step Day  Burnup  Contribution\n#\n";
    for(int i=0;i<burn_step;i++){
      cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
      real sum=0.;
      real sum_t=0.;
      for(int j=0;j<nucn;j++){
	for(int k=0;k<mednum_fuel;k++){
          real xs=xsc_1g[i][k][j]+xsf_1g[i][k][j];
          real den=fwd_nuc[i][0][k].get_dat(j);
          real volflx=vol_med[k]*total_flux[i][0][k];
	  real tmp=xs*den*volflx;
	  if(k==medid)sum+=tmp;
          sum_t+=tmp;
	};
      };
      cout<<sum/sum_t<<" "<<sum<<" "<<sum_t<<"\n";
    };
    cout<<"\n\n";

  cout<<"# Absorption rate in medium "<<medid<<"\n#\n#\n";
  cout<<"# Step Day  Burnup  Absorption rate\n#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";

    real sum=0.;
    vector<real> rr(nucn);
    for(int j=0;j<nucn;j++){
      real xs=xsc_1g[i][medid][j]+xsf_1g[i][medid][j];
      real den=fwd_nuc[i][0][medid].get_dat(j);
      rr[j]=xs*den;
      sum+=rr[j];
    };

    for(int j=0;j<nucn;j++){
      cout<<rr[j]/sum<<" ";
    };
    cout<<"\n";

  };
  cout<<"\n\n";

  }else{

    // Whole-system absorption rate
    cout<<"# Whole-system absorption rate\n#\n#\n";
    cout<<"# Step Day  Burnup  Absorption rate\n#\n";
    cout.setf(ios::scientific);
    cout.precision(4);
    for(int i=0;i<burn_step;i++){
      cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
            
      real sum=0.;
      vector<real> rr(nucn,0.);
      for(int j=0;j<nucn;j++){
	for(int k=0;k<mednum_fuel;k++){
          real xs=xsc_1g[i][k][j]+xsf_1g[i][k][j];
          real den=fwd_nuc[i][0][k].get_dat(j);
          real volflx=vol_med[k]*total_flux[i][0][k];
	  real tmp=xs*den*volflx;
          rr[j]+=tmp;
          sum+=tmp;
	};
      };

      for(int j=0;j<nucn;j++){
        cout<<rr[j]/sum<<" ";
      };
      cout<<"\n";

    };
    cout<<"\n\n";

  };

};

void MulticellBurner::ShowNumberDensity(int medid)
{
  cout<<"# Number density results\n";
  cout<<"# Nuclide ID list\n#\n";
  for(int j=0;j<nucn;j++){
    int matid=med[0].GetNuclideInTurn(j).GetMatnum();
    cout<<"# "<<j<<" "<<matid<<" "<<midt.Name(matid)<<"\n";
  };

  if(medid!=-1){

    // Medium-wise number density
  cout<<"#  Number density in medium "<<medid<<"\n#\n#\n";

  cout<<"# Step Day  Burnup  ND\n#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";

    for(int j=0;j<nucn;j++){
      real den=fwd_nuc[i][0][medid].get_dat(j);
      cout<<den<<" ";
    };
    cout<<"\n";

  };

  }else{

    // Whole-system number density
    cout<<"# Whole-system ND\n#\n#\n";

    cout<<"# Step Day  Burnup  ND\n#\n";
    cout.setf(ios::scientific);
    cout.precision(4);
    for(int i=0;i<burn_step;i++){
      cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
            

      vector<real> rr(nucn,0.);
      for(int j=0;j<nucn;j++){

        real sum=0.; 
        real sumvol=0.;
	for(int k=0;k<mednum_fuel;k++){
          real den=fwd_nuc[i][0][k].get_dat(j);
	  real tmp=den*vol_med[k];
          sum+=tmp;
	  sumvol+=vol_med[k];
	};
        cout<<sum/sumvol<<" ";

      };
      cout<<"\n";

    };

  };

};

void MulticellBurner::ShowNumberDensityReactionRate(int nucid, int nucid2)
{
  if(nucid2==-1)nucid2=nucid;

  int target=-1;
  int target2=-1;
  for(int i=0;i<nucn;i++){
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    if(matno==nucid)target=i;
    if(matno==nucid2)target2=i;
  };

  if(target==-1){
    cout<<"# Error in MulticellBurner::ShowCrossSectionNuclideReactionWise.\n";
    cout<<"# No nuclide of ID "<<nucid<<" exists.\n";
    exit(0);
  };
  if(target2==-1){
    cout<<"# Error in MulticellBurner::ShowCrossSectionNuclideReactionWise.\n";
    cout<<"# No nuclide of ID "<<nucid2<<" exists.\n";
    exit(0);
  };

  cout<<"#\n";
  cout<<"# A pair of ND [/cm3] and microscopic capture reaction rate [barn/cm2/s]";
  cout<<" of "<<midt.Name(nucid)<<" (ND) and "<<midt.Name(nucid2)<<" (RR) :";
  cout<<"\n#\n";
  cout<<"# (day)     (GWd/t)   ";
  for(int j=0;j<mednum_fuel;j++){
    cout<<"(Med ";
    WriteOut(j,2);
    cout<<")            ";
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<mednum_fuel;j++){
      real tmp;
      if(corrector_calc){
        tmp=xsc_1g_p[i][j][target2]*total_flux_p[i][0][j];
      }else{
        tmp=xsc_1g[i][j][target2]*total_flux[i][0][j];
      };
      cout<<fwd_nuc[i][0][j].get_dat(target)<<" "<<tmp<<" ";
    };
    cout<<"\n";
  };
  cout<<"\n\n";

};

void MulticellBurner::ShowNeutronFluxHistory()
{
  cout<<"#\n# +++ Time-dependent neutron flux [/cm2/s] +++\n";
  cout<<"# (day)     (GWd/t)   ";
  for(int j=0;j<mednum_fuel;j++){
    cout<<"(Med ";
    WriteOut(j,2);
    cout<<")  ";
  };
  cout<<"\n";
  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<burn_step;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<mednum_fuel;j++){
      cout<<total_flux[i][0][j]<<" ";
    };
    cout<<"\n";
  };
};

void MulticellBurner::ShowNumberDensityHistory(int prt_nuc,string *prt_nuc_nam,Burnup &bu,string opt,bool sum_only)
{
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

  // for volume-integrated quantity calculations
  vector< vector<real> > medsum(burn_step+1);
  for(int i=0;i<burn_step+1;i++){
    medsum[i].resize(prt_nuc,0.);
  };

  for(int iii=0;iii<mednum_fuel;iii++){
    cout<<"\n\n";
    cout<<"###################################################################\n";
    cout<<"#\n";
    cout<<"# Fuel medium ID : "<<iii<<"\n";

    cout<<"#\n";
    cout<<"#(Day)        (Burnup in GWd/t)    \n";    
    cout<<"#          [Average]  [Med-wise]   ";    
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
  }else if(opt=="w_gamma_per_thm"){
    cout<<"Gamma heat [W/tHM]\n";
  }else if(opt=="w_beta_per_thm"){
    cout<<"Beta Heat [W/tHM]\n";
  }else if(opt=="sv_per_thm_ing"){
    cout<<"Ingestion toxicity [Sv/tHM]\n";
  }else if(opt=="sv_per_thm_inh"){
    cout<<"Inhalation toxicity [Sv/tHM]\n";
  };
  cout<<"#                                  ";
  cout.setf(ios::scientific);
  cout.precision(4);
  if(!sum_only){
    for(int j=0;j<prt_nuc;j++){
      cout<<"  (";
      WriteOut(j+3,2);
      cout<<")     ";
    };
    cout<<"\n";
    cout<<"#                                  ";
    for(int j=0;j<prt_nuc;j++){
      cout<<prt_nuc_nam[j];
      int tmp=prt_nuc_nam[j].size();
      for(int k=0;k<11-tmp;k++){cout<<" ";};
    };
  };
  cout<<"  (sum)\n";


  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" "<<acburn_per_medium[iii][i]<<"   ";
    real sum=0.;
    for(int j=0;j<prt_nuc;j++){
      int sz=each_nuc_turn[j].size();
      real val=0.;
      real fuel_vol=vol_med[iii];
      for(int k=0;k<sz;k++){
        int tt=each_nuc_turn[j][k];
        int id=med[0].GetNuclideInTurn(tt).GetMatnum();
        //real dd=density_data[i][tt];
        real dd=fwd_nuc[i][0][iii].get_dat(tt);
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
        };
      };
      if(!sum_only)cout<<val<<" ";
      medsum[i][j]+=val*fuel_vol;
      sum+=val;
    };
    cout<<"  "<<sum<<"\n";
  };

  };

  real sumvol=0.;
  for(int i=0;i<mednum_fuel;i++){
    sumvol+=vol_med[i];
  };

  if(opt=="nd_per_vol"){

    cout<<"\n\n";
    cout<<"###################################################################\n";
    cout<<"#\n";
    cout<<"# Whole system-averaging\n";
    cout<<"#\n";
    cout<<"#(Day) (Burnup in GWd/t)\n";    
    cout<<"#          [Average]     ";    
    cout<<"N.D. [1e24/cm3]\n";
    cout<<"#                        ";
    cout.setf(ios::scientific);
    cout.precision(4);
    for(int j=0;j<prt_nuc;j++){
      cout<<"  (";
      WriteOut(j+2,2);
      cout<<")     ";
    };
    cout<<"\n";
    cout<<"#                        ";
    for(int j=0;j<prt_nuc;j++){
      cout<<prt_nuc_nam[j];
      int tmp=prt_nuc_nam[j].size();
      for(int k=0;k<11-tmp;k++){cout<<" ";};
    };
    cout<<"\n";

    for(int i=0;i<burn_step+1;i++){
      cout<<acday[i]<<" "<<acburn[i]<<"   ";
      for(int j=0;j<prt_nuc;j++){
        cout<<medsum[i][j]/sumvol<<" ";
      };
      cout<<"\n";
    };

  };

};

void MulticellBurner::ShowBurnupHistory()
{
  cout<<"#\n# +++ Time-dependent burnup [GWD/t] +++\n";
  cout<<"# (day)     (Target)  ";
  for(int j=0;j<mednum_fuel;j++){
    cout<<"(Med ";
    WriteOut(j,2);
    cout<<")  ";
  };
  cout<<"\n";
  cout<<"#           (or AVG) \n";
  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<burn_step+1;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" ";
    for(int j=0;j<mednum_fuel;j++){
      cout<<acburn_per_medium[j][i]<<" ";
    };
    cout<<"\n";
  };
};

void MulticellBurner::ShowFuelNeutronFlux(int st, bool excel)
{
  if(st>burn_step){
    cout<<"# Error in MulticellBurner::ShowFuelNeutronFlux\n";
    exit(0);
  };

  cout<<"#\n# Neutron flux per lethargy\n#\n";
  cout<<"#   burnup step : "<<st<<"\n";
  cout<<"#\n";
  if(excel)cout<<"# Upper\n";
  cout<<"# Eenrgy\n";
  cout<<"# [eV]     ";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"(";
    WriteOut(i,3);
    cout<<")     ";
  };
  cout<<"\n";

  vector<GroupData1D> aveflx(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    aveflx[i].put_imax(group);
    aveflx[i].set_zero();
  };


  for(int i=0;i<totm;i++){
    int medid=region_medium[i];
    if(medid<mednum_fuel){
      aveflx[medid]=aveflx[medid]+volflx_mesh[st][i];
    };
  };

  for(int i=0;i<mednum_fuel;i++){
    aveflx[i]=aveflx[i]*(1./vol_med[i]);
  };

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<group;i++){
    real e0=med[0].GetEnband().get_dat(i);
    real e1=med[0].GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    cout<<e0<<" ";
    for(int j=0;j<mednum_fuel;j++){
      cout<<aveflx[j].get_dat(i)/letwid<<" ";
    };
    cout<<"\n";
    if(excel){
      cout<<e1<<" ";
      for(int j=0;j<mednum_fuel;j++){
        cout<<aveflx[j].get_dat(i)/letwid<<" ";
      };
      cout<<"\n";
    };
  };
};

void MulticellBurner::ShowNumberDensityChange(int medid, real limit)
{
  cout<<"#\n#+++ Number density before/after burnup +++\n#\n";
  cout<<"#    (Medium : "<<medid<<" )\n";
  cout.setf(ios::scientific);
  //cout.precision(7);
  cout.precision(15);
  for(int ii=0;ii<2;ii++){
    for(int i=0;i<nucn;i++){
      real den_b=fwd_nuc[0][0][medid].get_dat(i);
      real den_a=fwd_nuc[burn_step][0][medid].get_dat(i);
      int matno=med[0].GetNuclideInTurn(i).GetMatnum();
      if(((ii==0&&matno>=900000)||(ii==1&&matno<900000))&&(den_a>limit||den_b>limit)){
	//cout<<"# "<<matno<<" "<<midt.Name(matno)<<"   "<<den_b<<"   "<<den_a<<"\n";
	cout<<" "<<matno<<" "<<den_a<<"\n";
	//cout<<"   "<<den_b<<"   "<<den_a<<"  "<<den_a*(xsc_1g[burn_step][i]+xsf_1g[burn_step][i])<<"\n";
	//cout<<den_a<<"\n";
      };
    };
  };
};

void MulticellBurner::WriteFileNumberDensityEOC(string mdir, string filename, int digit)
{
  int cyc_num=1;
  int cyc[]={burn_step};
  for(int i=0;i<mednum_fuel;i++){
    string filename2=filename+"_"+IntToString(i);
    WriteFileNumberDensity(cyc_num, cyc, mdir, filename2, i);
  };
};

void MulticellBurner::WriteFileNumberDensity(int cyc_num, int *cyc, string mdir, string filename, int medid, int digit)  
{
  ofstream fout;
  mdir.append(filename);
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Error in MulticellBurner::WriteFileNumberDensity.\n";
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(0);
  };

  for(int i=0;i<cyc_num;i++){
    int cc=cyc[i];
    if(cc<0||cc>burn_step){
      cout<<"# Error in MulticellBurner::WriteFileNumberDensity.\n";
      cout<<"# Number density data at cycle "<<cc<<" do NOT exist.\n";
      exit(0);
    };
  };

  if(medid>=mednum_fuel){
    cout<<"# Error in MulticellBurner::WriteFileNumberDensity.\n";
    cout<<"# Medium ID "<<medid<<" does NOT exist.\n";
    exit(0);
  };

  for(int i=0;i<cyc_num;i++){
    int cc=cyc[i];
    fout<<"      "<<cc<<"\n";
    fout.setf(ios::scientific);
    fout.precision(digit);
    fout<<"      "<<keff[cc]<<"\n";
    for(int j=0;j<nucn;j++){
      fout<<" "<<fwd_nuc[cc][0][medid].get_dat(j)<<"\n";
    };
  };

  fout.close();
};

void MulticellBurner::ReadFileNumberDensity(int cyc, string mdir, string filename, int medid)
{
  ifstream fin;
  mdir.append(filename);
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(0);
  };

  if(cyc<0||cyc>burn_step){
    cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
    cout<<"# Number density data at cycle "<<cyc<<" do NOT exist.\n";
    exit(0);
  };

  if(medid>=mednum_fuel){
    cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
    cout<<"# Medium ID "<<medid<<" does NOT exist.\n";
    exit(0);
  };

  bool end=false;
  while(!end){

    int dummy;
    fin>>dummy;
    if(dummy<0||dummy>100){
      cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
      exit(0);
    };
    if(dummy==cyc){
      real dummy2;
      fin>>dummy2;
      for(int i=0;i<nucn;i++){
        fin>>dummy2;
        fwd_nuc[0][0][medid].put_data(i,dummy2);
      };
      fin.close();
      return;
    }else{
      real dummy2;
      fin>>dummy2;
      for(int j=0;j<nucn;j++){
        fin>>dummy2;
      };
    };

  };

};

void MulticellBurner::ReadFileNumberDensityForInitialState(Burnup &bu, string mdir, string filename)
{
  for(int i=0;i<mednum_fuel;i++){
    string filename2=filename+"_"+IntToString(i);
    ReadFileNumberDensityForInitialState(bu, mdir, filename2, i);
  };
};

void MulticellBurner::ReadFileNumberDensityForInitialState(Burnup &bu, string mdir, string filename, int medid)  
{
  if(medid==0){
    bu.AddNuclideToMediumFromBurnupChain(med[0]);
    nucn=med[0].GetNucnum();
  };
  
  ifstream fin;
  mdir.append(filename);
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(0);
  };

  if(medid>=mednum_fuel){
    cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
    cout<<"# Medium ID "<<medid<<" does NOT exist.\n";
    exit(0);
  };

  int sz=med[0].GetNucnum();
  init_nucnum[medid]=sz;
  init_nucid[medid].resize(sz);
  init_nucden[medid].resize(sz);  
  for(int i=0;i<sz;i++){
    init_nucid[medid][i]=med[0].GetNuclideInTurn(i).GetMatnum();
  };
  
  bool end=false;
  while(!end){

    int dummy;
    fin>>dummy;
    if(dummy<0||dummy>100){
      cout<<"# Error in MulticellBurner::ReadFileNumberDensity.\n";
      exit(0);
    };
    if(dummy!=0){
      real dummy2;
      fin>>dummy2;
      for(int i=0;i<sz;i++){
        fin>>dummy2;
	init_nucden[medid][i]=dummy2;
	//cout<<medid<<" "<<init_nucid[medid][i]<<" "<<init_nucden[medid][i]<<"\n";	
      };
      fin.close();
      return;
    }else{
      real dummy2;
      fin>>dummy2;
      for(int j=0;j<nucn;j++){
        fin>>dummy2;
      };
    };

  };

};



void MulticellBurner::ForwardCalculation20210108(Burnup &bu, int med_target, bool adjoint)
{
  // This module is revised by Sasuga-kun to implement the new OWPC method,
  // and the old one is replaced by this one in 2021/01/08.
  //
  // Correlation between ND and reaction rate is replaced by
  // that between ND and one-group cross section.
  //
  // Later several inproper implementation was found by Chiba,
  // and this was replaced by the present [ForwardCalculation] in 2021/5/26.

  // Hard-coded parameters for weighted predictor-corrector
  //real wgt_nc=1.2; // relative weight for corrector
  real wgt_nc=1.0; // relative weight for corrector

  // OWPC to store positions of Gd-155 and 157
  int pos_gd155, pos_gd157;
  for(int j=0;j<nucn;j++){
    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    if(nucid==641550)pos_gd155=j;
    if(nucid==641570)pos_gd157=j;
  };

  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculation.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting for forward calculation +++

  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  flx_med.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  // +++ Array setting for predictor-corrector calculation
  //
  // - Generally multi-group microscopic cross section data at every burn steps
  //   are NOT stored because of their large required memory, but those are 
  //   required for GPT (adjoint) calculations, so those at every burnup steps
  //   are stored in the array [mic_sigx].
 
  {
  int tmp=1;
  if(adjoint)tmp=burn_step+1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };

  vector< vector<GroupData1D> > mic_sigf_c;
  vector< vector<GroupData1D> > mic_sigc_c;
  vector< vector<GroupData1D> > mic_sign2n_c;
  vector<GroupData1D> flx_med_c;
  if(corrector_calc){
    //int tmp=1;
    //if(adjoint)tmp=burn_step+1;
    int tmp=burn_step+1;
    xsc_1g_p.resize(tmp);
    xsn2n_1g_p.resize(tmp);
    xsf_1g_p.resize(tmp);
    for(int j=0;j<tmp;j++){
      xsc_1g_p[j].resize(mednum_fuel);
      xsn2n_1g_p[j].resize(mednum_fuel);
      xsf_1g_p[j].resize(mednum_fuel);
      for(int i=0;i<mednum_fuel;i++){
        xsc_1g_p[j][i].resize(nucn,0.);
        xsn2n_1g_p[j][i].resize(nucn,0.);
        xsf_1g_p[j][i].resize(nucn,0.);
      };
    };

    mic_sigf_c.resize(mednum_fuel); 
    mic_sigc_c.resize(mednum_fuel); 
    mic_sign2n_c.resize(mednum_fuel);
    for(int i=0;i<mednum_fuel;i++){
      mic_sigf_c[i].resize(nucn);
      mic_sigc_c[i].resize(nucn);
      mic_sign2n_c[i].resize(nucn);
    };

    flx_med_c.resize(mednum);
    for(int i=0;i<mednum;i++){
      flx_med_c[i].put_imax(group);
    };

    total_flux_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      total_flux_p[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        total_flux_p[i][j].resize(mednum_fuel);
      };
    };

  };

  // (OWPC)
  vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC (Rp)
  vector<real> np_gd5, np_gd7; // Np
  
    //sasuga addition
  vector<real>  old_n0155;
  vector<real>  old_n0157;
  vector<real>  old_r0155;
  vector<real>  old_r0157;
  vector<real>  old_np155;
  vector<real>  old_rp155;
  vector<real>  old_rp157;
  vector<real>  old_rc155;
  vector<real>  old_rc157;
  vector<real>  old_n0155_2;
  vector<real>  old_r0155_2;
  vector<real>  old_r0157_2;
  vector<real>  old_x0155;
  vector<real>  old_x0157;
  real n0_xx[2][mednum_fuel];
  real np_xx[2][mednum_fuel];
  real rp_xx[2][mednum_fuel];
  real rc_xx[2][mednum_fuel];

  real old_xsc5[mednum_fuel];
  real old_xsc7[mednum_fuel];
  real old_xsc5_2[mednum_fuel];
  real old_xsc7_2[mednum_fuel];
  real old_xscc_5[mednum_fuel];
  real old_xscc_7[mednum_fuel];
  real old_rc_5[mednum_fuel];
  real old_rc_7[mednum_fuel];
  real tmp_5[mednum_fuel];
  real tmp_7[mednum_fuel];

  real alpha_5[mednum_fuel];
  real alpha_7[mednum_fuel];
  real time;
  real time_be;
  real time_be_2;
  real burn_time_be;


  //sasuga addition
  old_n0155.resize(mednum_fuel);
  old_n0157.resize(mednum_fuel);
  old_r0155.resize(mednum_fuel);
  old_r0157.resize(mednum_fuel);
  old_np155.resize(mednum_fuel);
  old_rp155.resize(mednum_fuel);
  old_rp157.resize(mednum_fuel);
  old_rc155.resize(mednum_fuel);
  old_rc157.resize(mednum_fuel);
  old_n0155_2.resize(mednum_fuel);
  old_r0155_2.resize(mednum_fuel);
  old_r0157_2.resize(mednum_fuel);
  old_x0155.resize(mednum_fuel);
  old_x0157.resize(mednum_fuel);
 
  if(corrector_calc){
    rr_gd5.resize(mednum_fuel);
    rr_gd7.resize(mednum_fuel);
    np_gd5.resize(mednum_fuel);
    np_gd7.resize(mednum_fuel);
  };

  // +++ Array setting for adjoint calculation +++
  //
  //  [fwd_nuc_int] is a number density at time-mesh-center point.
  //  Time-averaged number density during one burnup step is calculated 
  //  from those at beginning, center and end of this burnup step
  if(adjoint){
    fwd_nuc_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
  	  fwd_nuc_int[i][j][k].put_imax(nucn);
	};
      };
    };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();

  // +++ Forward burn-up calculation +++
  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  // +++ Eigenvalue calculation
  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);  

  if(adjoint){

    // +++ array setting for angular flux storing in adjoint calculations.
    int totsn=lat.GetQuad().GetSN();
    // Quadrature setting
    quad.Initialize(2,aflx_legendre);
    quad.PutSN(totsn);
    real *mui=new real[totsn];
    real *eai=new real[totsn];
    real *xii=new real[totsn];
    real *w=new real[totsn];
    for(int i=0;i<totsn;i++){
      mui[i]=lat.GetQuad().GetMu(i);
      eai[i]=lat.GetQuad().GetEata(i);
      xii[i]=lat.GetQuad().GetXi(i);
      w[i]=lat.GetQuad().GetOmega(i);
    };
    quad.PutData(mui,eai,xii,w);
    quad.CalValue();

    if(aflx_legendre==-1){
      volaflx_mesh.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };

  };

  // +++ GPT-PC +++++++++++++++++++++++++++++++++++++++++++++++++
  if(adjoint&&corrector_calc){
    fwd_nuc_p.resize(burn_step+1);
    fwd_nuc_p_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_p[i].resize(sub_step+1);
      fwd_nuc_p_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_p[i][j].resize(mednum_fuel);
        fwd_nuc_p_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
          fwd_nuc_p[i][j][k].put_imax(nucn);
          fwd_nuc_p_int[i][j][k].put_imax(nucn);
        };
      };
    };
    volflx_mesh_p.resize(burn_step);  
    power_factor_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      volflx_mesh_p[i].resize(totm);
      for(int j=0;j<totm;j++){
        if(region_medium[j]<mednum_fuel){
          volflx_mesh_p[i][j].put_imax(group);
        };
      };
      int sub_step=sub_step_list[i];
      power_factor_p[i].resize(sub_step+1);
    };
    int totsn=lat.GetQuad().GetSN();

    if(aflx_legendre==-1){
      volaflx_mesh_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh_p[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh_p[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl_p[i].resize(totm);
        for(int j=0;j<totm;j++){
    	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl_p[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real abs_frac[burn_step+1];
  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;
    if(adjoint)bstmp=st;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    //lat.MedClear();
    for(int i=0;i<mednum_fuel;i++){
      // +++ Self-shielding calculation
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048);      
      med[0].CalMacroFromMicro();

      // +++ Macro&Micro data storing
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };
      if(st!=burn_step&&adjoint){
        macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
      };

      if(st==0){
        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();
      }else{
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
      };
    };
    if(st==0){
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };
    }else{
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };
    };
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);

    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    if(adjoint)lat.PutWriteFlux();
    lat.Print();
    keff[st]=lat.CalIgen();

    //cout<<"# Total inner iteration : "<<lat.GetTotalInnerIterationSum()<<"\n"; exit(0);

    // --- (Generating medium instances include collapsing cross section data) 
    if(collapsing){

    for(int i=0;i<mednum_fuel;i++){
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      //med[0].CalMacroFromMicro();

      int ngrp=9;
      int bgrp[]={21,46,91,115,125,134,151,159,171};
      GroupData1D flx=lat.GetIntegratedFlux(i);
      Medium bmed=med[0].Cond(ngrp,bgrp,flx,flx,true);
      bmed.GetMacxs().GetData1d(d).copy(bmed.GetFlux());
      // Neutron flux is stored as diffusion coefficient (temporal treatment)
      string filename="bmed_med"+IntToString(i)+"_st"+IntToString(st);
      bmed.WriteFile(collapse_dir,filename,true);
    };

    };
    // --------------------------------------------------------------------------------


    // ------------------------------

    real sum_tot=0.;
    real sum_fmed=0.;
    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);

      real tmp=flx*lat.GetMedium(i).GetMacxs().GetData1d(siga);
      sum_tot+=tmp;
      if(i<mednum_fuel)sum_fmed+=tmp;

      flx_med[i]=flx*(1./vol_med[i]);
      // +++ One-group cross section storing (at predictor step)
      if(i<mednum_fuel){
        for(int j=0;j<nucn;j++){
          if(nuclide_info[j]!=0){
            if(nuclide_info[j]==1){
   	      xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	    };
	    xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	    xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
          };
	};
      };
    };
    abs_frac[st]=sum_fmed/sum_tot; // absorption rate fraction in fuel media
    cout<<"# Absorption fraction : "<<abs_frac[st]<<"\n";
  
    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
          if(adjoint){
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
  	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };
	  };
	};
      };

      // +++ Burnup calculation
      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
          //power_org=CalculationPinPower(bu,st,j,med_target,sumflx*vol_med[med_target]);
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        accumulated_day+=burn_span;
	//total_flux[st][j].resize(mednum_fuel);
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){ 
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //   CORRECTOR CALCULATION
      if(corrector_calc){

	cout<<"# Corrector calculation ...\n";

	// (For OWPC-2018)	
	for(int i=0;i<mednum_fuel;i++){
	  real tmp=total_flux[st][0][i];
  	  rr_gd5[i]=tmp*xsc_1g[st][i][pos_gd155];
  	  rr_gd7[i]=tmp*xsc_1g[st][i][pos_gd157];
	};

	// ++++++++++++++++++++++++

	// (Predictor calculation results are stored in the array of [XXX_p].)
        swap(total_flux_p[st], total_flux[st]);
        swap(xsc_1g_p[st], xsc_1g[st]);
        swap(xsn2n_1g_p[st], xsn2n_1g[st]);
        swap(xsf_1g_p[st], xsf_1g[st]);

	if(adjoint){
          swap(keff_p[st], keff[st]);
          swap(fwd_nuc_p[st], fwd_nuc[st]);
          swap(fwd_nuc_p_int[st], fwd_nuc_int[st]);
          swap(power_factor_p[st], power_factor[st]);
	  swap(volflx_mesh_p[st], volflx_mesh[st]);
	  if(aflx_legendre==-1){
            swap(volaflx_mesh_p[st], volaflx_mesh[st]);
	  }else{
  	    swap(volaflx_pl_p[st], volaflx_pl[st]);
	  };
	  swap(macxs_p[st], macxs[st]);
	  fwd_nuc[st][0]=fwd_nuc_p[st][0];
	};

      for(int i=0;i<mednum_fuel;i++){
        // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
        PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
  	  };
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        }else{
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        };
        //opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);	
        med[0].CalMacroFromMicro();
        //opc.CalFissionSpectrumMatrix(med[0],xslib);

        // +++ Macro&Micro data storing
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
	if(adjoint){
          macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
	};
      };

      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };

      lat.PutThermalIteration(3);
      lat.PutPL(0);
      lat.NoCMRAcceleration();
      real keff_corr=lat.CalIgen();
      
      // ------------------------------

      if(adjoint)keff[st]=keff_corr;

      for(int i=0;i<mednum;i++){
        GroupData1D flx=lat.GetIntegratedFlux(i);
        flx_med_c[i]=flx*(1./vol_med[i]);
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1)xsf_1g[st][i][j]=mic_sigf_c[i][j].Cond(flx);
	      xsc_1g[st][i][j]=mic_sigc_c[i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n_c[i][j].Cond(flx);
	    };
          };
	};
      };

      if(adjoint){
        for(int i=0;i<totm;i++){
          if(region_medium[i]<mednum_fuel){
            real vol=lat.GetMesh(i).GetVolume();
            volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol;
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };
	  };
        };


      };

      real acburn_pre=accumulated_burn-acburn[st];
      vector<real> acburn_pre_per_medium(mednum_fuel);
      for(int j=0;j<mednum_fuel;j++){
	acburn_pre_per_medium[j]=accumulated_burn_per_medium[j]-acburn_per_medium[j][st];
      };
      // accumulated burnup calculated by the predictor step
      
      //sasuga addition  断面積を相対的に変動させる
      for(int i=0;i<mednum_fuel;i++){	  
	real xsc0_5=xsc_1g_p[st][i][pos_gd155]; //rp
	real xsc0_7=xsc_1g_p[st][i][pos_gd157]; //rp
	real xscc_5=xsc_1g[st][i][pos_gd155];//Rcに相当するXc
	real xscc_7=xsc_1g[st][i][pos_gd157];
	if(st>0){
	    real n0_5=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real np_5=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real t_5=old_xsc5[i]+(xsc0_5-old_xsc5[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);//前のステップでの正しいXcを計算
	    real t_7=old_xsc7[i]+(xsc0_7-old_xsc7[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
	    real xsc155,xsc157;
	    if(st>1){
	      real a=n0_5, b=old_n0155[i], c=old_n0155_2[i];
	      real aa_5=xsc0_5, bb_5=old_xsc5[i], cc_5=old_xsc5_2[i];
	      real aa_7=xsc0_7, bb_7=old_xsc7[i], cc_7=old_xsc7_2[i];
	      real X_5=old_np155[i];
	      t_5=aa_5*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_5*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_5*(X_5-a)*(X_5-b)/(c-a)/(c-b);//前のステップでの正しいXcを計算
	      t_7=aa_7*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_7*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_7*(X_5-a)*(X_5-b)/(c-a)/(c-b);
	    };	          
	    //xsc155=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)/burn_time_be*burn_time[st]; //time
	    //xsc157=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)/burn_time_be*burn_time[st]; //time
	    xsc155=xsc_1g[st][i][pos_gd155]-(old_xscc_5[i]-t_5)/acburn[st-1]*acburn[st]; //burnup
	    xsc157=xsc_1g[st][i][pos_gd157]-(old_xscc_7[i]-t_7)/acburn[st-1]*acburn[st]; //burnup
	  
	  real d_total_flux=(total_flux[st][0][i]-total_flux[st-1][0][i])/total_flux[st-1][0][i];
	  if(fabs(d_total_flux)<0.1){
	    xsc_1g[st][i][pos_gd155]=xsc155;
	    xsc_1g[st][i][pos_gd157]=xsc157;
	  };//if d_total flux
	  
	  old_xsc5_2[i]=old_xsc5[i];
	  old_xsc7_2[i]=old_xsc7[i];
	};
	old_xsc5[i]=xsc0_5;
	old_xsc7[i]=xsc0_7;
	old_xscc_5[i]=xscc_5;
	old_xscc_7[i]=xscc_7;
	
      };
      burn_time_be=burn_time[st];
      //-------------------------- sasuga addition
      
      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med_c[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med_c[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          //accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        //accumulated_day+=burn_span;
        //delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
	  //CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	  /*
	  //sasuga addition  反応率を相対的に変動させる
	  if(j==0){
	    real r0_5=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][j][i]; //rp
	    real r0_7=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][j][i]; //rp
	    real rc_5=xsc_1g[st][i][pos_gd155]*total_flux[st][j][i];
	    real rc_7=xsc_1g[st][i][pos_gd157]*total_flux[st][j][i];
	    cout<<i<<"  before  "<<xsc_1g[st][i][pos_gd155]<<endl;
	    if(st>0){
	      real n0_5=fwd_nuc[st][0][i].get_dat(pos_gd155);
	      real np_5=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	      real t_5=old_r0155[i]+(r0_5-old_r0155[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
	      real t_7=old_r0157[i]+(r0_7-old_r0157[i])/(n0_5-old_n0155[i])*(old_np155[i]-old_n0155[i]);
	      real rcc155,rcc157;
	      if(st>1){
		real a=n0_5, b=old_n0155[i], c=old_n0155_2[i];
		real aa_5=r0_5, bb_5=old_r0155[i], cc_5=old_r0155_2[i];
		real aa_7=r0_7, bb_7=old_r0157[i], cc_7=old_r0157_2[i];
		real X_5=old_np155[i];
		t_5=aa_5*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_5*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_5*(X_5-a)*(X_5-b)/(c-a)/(c-b);
		t_7=aa_7*(X_5-b)*(X_5-c)/(a-b)/(a-c)+bb_7*(X_5-a)*(X_5-c)/(b-a)/(b-c)+cc_7*(X_5-a)*(X_5-b)/(c-a)/(c-b);
	    
		X_5=old_n0155[i];
		real al_5_be=aa_5*(2*X_5-b-c)/(a-b)/(a-c)+bb_5*(2*X_5-a-c)/(b-a)/(b-c)+cc_5*(2*X_5-a-b)/(c-a)/(c-b);
		real al_7_be=aa_7*(2*X_5-b-c)/(a-b)/(a-c)+bb_7*(2*X_5-a-c)/(b-a)/(b-c)+cc_7*(2*X_5-a-b)/(c-a)/(c-b);
	    
		X_5=n0_5;
		real al_5=aa_5*(2*X_5-b-c)/(a-b)/(a-c)+bb_5*(2*X_5-a-c)/(b-a)/(b-c)+cc_5*(2*X_5-a-b)/(c-a)/(c-b);
		real al_7=aa_7*(2*X_5-b-c)/(a-b)/(a-c)+bb_7*(2*X_5-a-c)/(b-a)/(b-c)+cc_7*(2*X_5-a-b)/(c-a)/(c-b);
	    
		//rcc155=xsc_1g[st][i][pos_gd155]*total_flux[st][j][i]-(R_old_5[i]-t_5)/R_old_5[i]*xsc_1g[st][i][pos_gd155]*total_flux[st][j][i];
		//rcc157=xsc_1g[st][i][pos_gd157]*total_flux[st][j][i]-(R_old_7[i]-t_7)/R_old_7[i]*xsc_1g[st][i][pos_gd157]*total_flux[st][j][i]; //ip補正
	    
		rcc155=xsc_1g[st][i][pos_gd155]*total_flux[st][j][i]-(old_rc_5[i]-t_5)/al_5_be*al_5;
		rcc157=xsc_1g[st][i][pos_gd157]*total_flux[st][j][i]-(old_rc_7[i]-t_7)/al_7_be*al_7; //gamma補正
	      }else{
		rcc155=xsc_1g[st][i][pos_gd155]*total_flux[st][j][i]-(old_rc_5[i]-t_5)/old_rc_5[i]*xsc_1g[st][i][pos_gd155]*total_flux[st][j][i];
		rcc157=xsc_1g[st][i][pos_gd157]*total_flux[st][j][i]-(old_rc_7[i]-t_7)/old_rc_7[i]*xsc_1g[st][i][pos_gd157]*total_flux[st][j][i]; //ip補正
	      };	    
	      xsc_1g[st][i][pos_gd155]=rcc155/total_flux[st][j][i];
	      xsc_1g[st][i][pos_gd157]=rcc157/total_flux[st][j][i];

	      old_rc_5[i]=rc_5;
	      old_rc_7[i]=rc_7;
	     };
	  };
	  //-------------------------- sasuga addition
	  */
	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	};
	
      }; // end of sub-step
      
      // (Variables required for OWPC)
      real dt=burn_time[st]*60*60*24;
      real time_mesh=10000;
      real dt_r=dt/time_mesh;

      bool owpc_corr=true; // If [true], OWPC treatment for Gd-155 and -157 with various correlation conditions are used.

      real omega_155, omega_157; // various correlation conditions

      for(int i=0;i<mednum_fuel;i++){

	  // OWPC treatment for Gd-155 and -157 with various correlation conditions
	if(owpc_corr){

	  //consider 154,156--------
	  /*
	  real  n_before_155=fwd_nuc[st][0][i].get_dat(pos_gd154);
	  real  r_before_155=rr_gd4[i];
	  real  n_before_157=fwd_nuc[st][0][i].get_dat(pos_gd156);
	  real  r_before_157=rr_gd6[i];
	  */

	  if(st>0){
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	    real x0_155=xsc_1g_p[st][i][pos_gd155];//rpに対応する断面積
	    real x0_157=xsc_1g_p[st][i][pos_gd157];//rpに対応する断面積
	    real xp_155=xsc_1g[st][i][pos_gd155]; //rcに対応する断面積
	    real xp_157=xsc_1g[st][i][pos_gd157]; //rcに対応する断面積
	    real n0_r_155;
	    real n0_r_157;
	    
	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt);//Not consider 154,156
	    	    
	    //quadratic model---------------------------
	    real a_5,b_5,c_5,a_7,b_7,c_7;
	    real aa_5,bb_5,cc_5,aa_7,bb_7,cc_7;
	    a_5=n0_155;
	    b_5=np_155;
	    c_5=old_n0155[i];
	    a_7=n0_157;
	    b_7=np_157;
	    c_7=old_n0157[i];
	    real X_5=np_p_155;
	    real X_7=np_p_157;

	    /*
	    //-------------------------reaction rate base-------------------------------------
	    aa_5=r0_155;
	    bb_5=rp_155;
	    cc_5=old_r0155[i];
	    
	    aa_7=r0_157;
	    bb_7=rp_157;
	    cc_7=old_r0157[i];


	    real rc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    real rc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    //real rc_157=(aa_7*(X_7-b_7)*(X_7-c_7)/(a_7-b_7)/(a_7-c_7))+(bb_7*(X_7-a_7)*(X_7-c_7)/(b_7-a_7)/(b_7-c_7))+(cc_7*(X_7-a_7)*(X_7-b_7)/(c_7-a_7)/(c_7-b_7));
	    //----------------------------------------------------------------------------------------
	    */

	    //---------------------------------------cross section base-------------------------------------
	    aa_5=x0_155;
	    bb_5=xp_155;
	    cc_5=old_x0155[i];
	    
	    aa_7=x0_157;
	    bb_7=xp_157;
	    cc_7=old_x0157[i];
	        
	    real xc_155=(aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //xc in toyproblem
	    real xc_157=(aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5))+(bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5))+(cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5)); //rc in toyproblem
	    //----------------------------------------------------------------

	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	    real xp_r_155=x0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	    real xp_r_157=x0_157;

	    //------------------------------------start toy problem calculation----------------------------
	    for(int k=0;k<time_mesh;k++){

	      /*
	      //-----reaction rate base-------
	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      X_5=np_r_155;
	      X_7=np_r_157;
	      
	      real rc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	      real rc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	    //real rc_r_157=aa_7*(X_7-b_7)*(X_7-c_7)/(a_7-b_7)/(a_7-c_7)+bb_7*(X_7-a_7)*(X_7-c_7)/(b_7-a_7)/(b_7-c_7)+cc_7*(X_7-a_7)*(X_7-b_7)/(c_7-a_7)/(c_7-b_7);

	    real r_r_155=(rp_r_155+rc_r_155)*0.5;
	    real r_r_157=(rp_r_157+rc_r_157)*0.5;

	    n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	    n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r); //Not consider 154,156
	    
	    rp_r_155=rc_r_155;
	    rp_r_157=rc_r_157;
	    //----------------------------------------
	    */
	    //-------cross section base-------
	      real np_r_155=n0_r_155*exp(-xp_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-xp_r_157*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156

	      X_5=np_r_155;
	      X_7=np_r_157;
	    real xc_r_155=aa_5*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_5*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_5*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);
	    real xc_r_157=aa_7*(X_5-b_5)*(X_5-c_5)/(a_5-b_5)/(a_5-c_5)+bb_7*(X_5-a_5)*(X_5-c_5)/(b_5-a_5)/(b_5-c_5)+cc_7*(X_5-a_5)*(X_5-b_5)/(c_5-a_5)/(c_5-b_5);

	    real x_r_155=(xp_r_155+xc_r_155)*0.5;
	    real x_r_157=(xp_r_157+xc_r_157)*0.5;

	    n0_r_155=n0_r_155*exp(-x_r_155*total_flux_p[st][0][i]*1e-24*dt_r);//Not consider 154,156
	    n0_r_157=n0_r_157*exp(-x_r_157*total_flux_p[st][0][i]*1e-24*dt_r); //Not consider 154,156
	    
	    xp_r_155=xc_r_155;
	    xp_r_157=xc_r_157;
	    //----------------------------------------
	    
	    
	    };
	    //-----------------------------------------end toy problem calculation------------------------------------------
	    
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;
	    /*
	    //-------reaction rate base----------
	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	    //------------------------------------------
	    */
	    //-------cross section base----------
	    omega_155=(R_reference_155-x0_155*total_flux_p[st][0][i])/(xc_155*total_flux_p[st][0][i]-x0_155*total_flux_p[st][0][i]);
	    omega_157=(R_reference_157-x0_157*total_flux_p[st][0][i])/(xc_157*total_flux_p[st][0][i]-x0_157*total_flux_p[st][0][i]);
	    //--------------------------------------
	    
	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;
	    
	    old_n0155_2[i]=old_n0155[i];
	    old_r0155_2[i]=old_r0155[i];
	    old_r0157_2[i]=old_r0157[i];
	    
	    /*
	    old_n0155_2[i]=old_np155[i];
	    old_r0155_2[i]=m_rc_155;
	    old_r0157_2[i]=m_rc_157;
	    */
	  
	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;
	    old_x0155[i]=x0_155;
	    old_x0157[i]=x0_157;
	    
	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;
	      	    
	  }else{
	  
	
	    real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	    real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	    real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	    real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	    real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	    real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	    real r0_155=xsc_1g_p[st][i][pos_gd155]*total_flux_p[st][0][i]; //rp
	    real r0_157=xsc_1g_p[st][i][pos_gd157]*total_flux_p[st][0][i]; //rp
	    real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	    real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	
	    real n0_r_155;
	    real n0_r_157;

	    real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem //Not consider 154,156
	    real np_p_157=n0_157*exp(-r0_157*1e-24*dt);//Not consider 154,156
	  
	    //real np_p_155=(n0_155-r_before_155/r0_155*n_before_155)*exp(-r0_155*1e-24*dt)+r_before_155/r0_155*n_before_155; // consider 154,156
	    //real np_p_157=(n0_157-r_before_157/r0_157*n_before_157)*exp(-r0_157*1e-24*dt)+r_before_157/r0_157*n_before_157; // consider 154,156

	    // (Correlation to Gd-155 number density)
	    real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	    real alpha_157=(rp_157-r0_157)/(np_155-n0_155);    
	    
	    //linear model------------------------------------------------------------
	    // (Correlation to Gd-155 number density)	
	    real rc_155=r0_155+(np_p_155-n0_155)*alpha_155; // Rc in toy-problem
	    real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	   
	    n0_r_155=n0_155;
	    real rp_r_155=r0_155;
	  
	    n0_r_157=n0_157;
	    real rp_r_157=r0_157;
	  
	    for(int k=0;k<time_mesh;k++){

	      real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);//Not consider 154,156
	      real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);//Not consider 154,156

	      //real np_r_155=(n0_r_155-r_before_155/rp_r_155*n_before_155)*exp(-rp_r_155*1e-24*dt_r)+r_before_155/rp_r_155*n_before_155; // consider 154,156
	      //real np_r_157=(n0_r_157-r_before_157/rp_r_157*n_before_157)*exp(-rp_r_157*1e-24*dt_r)+r_before_157/rp_r_157*n_before_157; // consider 154,156
	      
	      //linear model
	      // (Correlation to Gd-155 number density)
	      real rc_r_155=r0_155+(np_r_155-n0_155)*alpha_155;
	      real rc_r_157=r0_157+(np_r_155-n0_155)*alpha_157;
	
	      real r_r_155=(rp_r_155+rc_r_155)*0.5;
	      real r_r_157=(rp_r_157+rc_r_157)*0.5;

	      n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);//Not consider 154,156
	      n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r);//Not consider 154,156

	      rp_r_155=rc_r_155;
	      rp_r_157=rc_r_157;

	      };
	        
	    real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	    real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

	    omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	    omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	    
	    if((0>omega_155)||(1<omega_155)) omega_155=0.5;
	    if((0>omega_157)||(1<omega_157)) omega_157=0.5;
	    /*
	    // omega is defined by number density----------------------------------
	    real nc_c_155=n0_155*exp(-rp_155*1e-24*dt);
	    real nc_c_157=n0_157*exp(-rp_157*1e-24*dt);

	    omega_155=(n0_r_155-np_p_155)/(nc_c_155-np_p_155);
	    omega_157=(n0_r_157-np_p_157)/(nc_c_157-np_p_157);
	    //------------------------------------------------------------------------------------
	    */
	    old_n0155[i]=n0_155;
	    old_n0157[i]=n0_157;
	    old_r0155[i]=r0_155;
	    old_r0157[i]=r0_157;

	    old_np155[i]=np_155;
	    old_rc155[i]=rp_155;
	    old_rc157[i]=rp_157;
	    };
	};
	  
	for(int j=0;j<nucn;j++){
	
	  real np=fwd_nuc[st+1][0][i].get_dat(j);      // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
	  real n_next=exp((log(np)+log(nc))*0.5);

          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  if(nucid==641550||nucid==641570){

  	    real n0=fwd_nuc[st][0][i].get_dat(j);

	    if((n0-np)/np>1e-2){
	    
		if(!owpc_corr){

		  // -- OWPC ----------------------- 
		  real r0=rr_gd5[i]; //Rp
		  //real r0mm=rr_gd4[i]*fwd_nuc[st][0][i].get_dat(pos_gd154); //N Rp
		  if(nucid==641570)r0=rr_gd7[i];
		  //if(nucid==641570)r0mm=rr_gd6[i]*fwd_nuc[st][0][i].get_dat(pos_gd156);

		  real rp=xsc_1g[st][i][j]*total_flux[st][0][i]; //Rc
		  real alpha=(rp-r0)/(np-n0);
		  
		  real np_p=n0*exp(-r0*1e-24*dt);   // Nc in toy-problem
		  real rc=alpha*np_p+(r0-n0*alpha); // Rc in toy-problem
	    
		  real n0_r=n0;
		  real rp_r=r0;
		  for(int k=0;k<time_mesh;k++){
		    real np_r=n0_r*exp(-rp_r*1e-24*dt_r);
		    real rc_r=alpha*np_r+(r0-n0*alpha);
		    real r_r=(rp_r+rc_r)/2;
		    n0_r=n0_r*exp(-r_r*1e-24*dt_r);
		    rp_r=rc_r;
		  };
		  
		  real R_reference=(log(n0)-log(n0_r))/dt*1e24;
		  real omega=(R_reference-r0)/(rc-r0);
		  
		  // --------------------------------
	    	    
		  // -- PPC -------------------------------
		  		  
		    real R_p=-log(np/n0)/dt;
		    real R_c=-log(nc/n0)/dt;
		    real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
		    R_reference=(R_p+R_c_c)/2;
		    omega=(R_reference-R_p)/(R_c-R_p);
		  
		  // --------------------------------
	    
		  n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応

		  //cout<<"# OWPC factor is "<<omega<<" in medium "<<i<<" and nuclide "<<nucid<<"\n";

		};
	    
		// OWPC treatment for Gd-155 and -157 with various correlation conditions
		if(owpc_corr){
		  if(nucid==641550){
		    n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
		    //n_next=np*(1-omega_155)+nc*omega_155;
		  }else if(nucid==641570){
		    n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
		    //n_next=np*(1-omega_157)+nc*omega_157;
		  };
		  
		};
	      }else{
	      
	      //cout<<"# OWPC is not adopted in medium "<<i<<" and nuclide "<<nucid<<"\n";
		
	      };

	  };

	  // ... Conventional PC ...
	  
          //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
          //n_next=(np+nc)*0.5; // Conventional PC (in linear)
          //n_next=nc;

          if(adjoint)        n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
          if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC

  	  // (For OWPC-2018)	
  	  if(nucid==641550)np_gd5[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
  	  if(nucid==641570)np_gd7[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd157);

          fwd_nuc[st+1][0][i].put_data(j,n_next);

	};

      };
      
      // Adjustment of accumulated burn
      for(int i=0;i<mednum_fuel;i++){
        real acburn_cor=accumulated_burn_per_medium[i]-acburn_per_medium[i][st]-acburn_pre_per_medium[i];
        accumulated_burn_per_medium[i]=acburn_per_medium[i][st]+(acburn_pre_per_medium[i]+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };

      if(input_flux_level){
        real acburn_cor=accumulated_burn-acburn[st]-acburn_pre;
        accumulated_burn=acburn[st]+(acburn_pre+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };
      cout<<"#     ... terminated.\n";
      }; // END OF CORRECTOR CALCULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    }; // end part of [If (st!=burn_step) ]
 
  }; // loop-end of burnup step

};

void MulticellBurner::ForwardCalculationOld(Burnup &bu, int med_target, bool adjoint)
{
  bool owpc_corr=true; // If [true], OWPC treatment for Gd-155 and -157 with various correlation conditions are used.  
  
  // ... Hard-coded parameters for weighted predictor-corrector
  //real wgt_nc=1.2; // relative weight for corrector
  real wgt_nc=1.0; // relative weight for corrector

  // ... To store positions of Gd-155 and 157 in A-OWPC
  int pos_gd155, pos_gd157;
  for(int j=0;j<nucn;j++){
    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
    if(nucid==641550)pos_gd155=j;
    if(nucid==641570)pos_gd157=j;
  };

  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculation.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  // +++ Array setting for forward calculation +++

  fwd_nuc.resize(burn_step+1);
  xsc_1g.resize(burn_step+1);
  xsn2n_1g.resize(burn_step+1);
  xsf_1g.resize(burn_step+1);
  total_flux.resize(burn_step);
  delt.resize(burn_step);
  power_factor.resize(burn_step);

  for(int i=0;i<burn_step+1;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].resize(mednum_fuel);
      for(int k=0;k<mednum_fuel;k++){
	fwd_nuc[i][j][k].put_imax(nucn);
      };
    };
    xsc_1g[i].resize(mednum_fuel);
    xsn2n_1g[i].resize(mednum_fuel);
    xsf_1g[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      xsc_1g[i][j].resize(nucn,0.);
      xsn2n_1g[i][j].resize(nucn,0.);
      xsf_1g[i][j].resize(nucn,0.);
    };
  };

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    total_flux[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    for(int j=0;j<sub_step;j++){
      total_flux[i][j].resize(mednum_fuel);
    };
  };

  vector<GroupData1D> flx_med;
  flx_med.resize(mednum);
  for(int j=0;j<mednum;j++){
    flx_med[j].put_imax(group);
  };

  volflx_mesh.resize(burn_step); 
  for(int i=0;i<burn_step;i++){
    volflx_mesh[i].resize(totm);
    for(int j=0;j<totm;j++){
      if(region_medium[j]<mednum_fuel){
        volflx_mesh[i][j].put_imax(group);
      };
    };
  };

  // +++ Array setting for predictor-corrector calculation
  //
  // - Generally multi-group microscopic cross section data at every burn steps
  //   are NOT stored because of their large required memory, but those are 
  //   required for GPT (adjoint) calculations, so those at every burnup steps
  //   are stored in the array [mic_sigx].
 
  {
  int tmp=1;
  if(adjoint)tmp=burn_step+1;
  mic_sigf.resize(tmp); 
  mic_sigc.resize(tmp); 
  mic_sign2n.resize(tmp);
  for(int i=0;i<tmp;i++){
    mic_sigf[i].resize(mednum_fuel);
    mic_sigc[i].resize(mednum_fuel);
    mic_sign2n[i].resize(mednum_fuel);
    for(int j=0;j<mednum_fuel;j++){
      mic_sigf[i][j].resize(nucn);
      mic_sigc[i][j].resize(nucn);
      mic_sign2n[i][j].resize(nucn);
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[i][j][k].put_imax(group);
          };
          mic_sigc[i][j][k].put_imax(group);
          mic_sign2n[i][j][k].put_imax(group);
        };
      };
    };
  };
  };

  vector< vector<GroupData1D> > mic_sigf_c;
  vector< vector<GroupData1D> > mic_sigc_c;
  vector< vector<GroupData1D> > mic_sign2n_c;
  vector<GroupData1D> flx_med_c;
  if(corrector_calc){
    //int tmp=1;
    //if(adjoint)tmp=burn_step+1;
    int tmp=burn_step+1;
    xsc_1g_p.resize(tmp);
    xsn2n_1g_p.resize(tmp);
    xsf_1g_p.resize(tmp);
    for(int j=0;j<tmp;j++){
      xsc_1g_p[j].resize(mednum_fuel);
      xsn2n_1g_p[j].resize(mednum_fuel);
      xsf_1g_p[j].resize(mednum_fuel);
      for(int i=0;i<mednum_fuel;i++){
        xsc_1g_p[j][i].resize(nucn,0.);
        xsn2n_1g_p[j][i].resize(nucn,0.);
        xsf_1g_p[j][i].resize(nucn,0.);
      };
    };

    mic_sigf_c.resize(mednum_fuel); 
    mic_sigc_c.resize(mednum_fuel); 
    mic_sign2n_c.resize(mednum_fuel);
    for(int i=0;i<mednum_fuel;i++){
      mic_sigf_c[i].resize(nucn);
      mic_sigc_c[i].resize(nucn);
      mic_sign2n_c[i].resize(nucn);
    };

    flx_med_c.resize(mednum);
    for(int i=0;i<mednum;i++){
      flx_med_c[i].put_imax(group);
    };

    total_flux_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      total_flux_p[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        total_flux_p[i][j].resize(mednum_fuel);
      };
    };

  };

  // (OWPC)
  vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC (Rp)
  vector<real> np_gd5, np_gd7; // Np
  if(corrector_calc){
    rr_gd5.resize(mednum_fuel);
    rr_gd7.resize(mednum_fuel);
    np_gd5.resize(mednum_fuel);
    np_gd7.resize(mednum_fuel);
  };

  // +++ Array setting for adjoint calculation +++
  //
  //  [fwd_nuc_int] is a number density at time-mesh-center point.
  //  Time-averaged number density during one burnup step is calculated 
  //  from those at beginning, center and end of this burnup step
  if(adjoint){
    fwd_nuc_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
  	  fwd_nuc_int[i][j][k].put_imax(nucn);
	};
      };
    };
  };

  // +++ Initial number density setting +++
  for(int i=0;i<mednum_fuel;i++){
    fwd_nuc[0][0][i].set_zero();
    for(int j=0;j<init_nucnum[i];j++){
      int idtmp=init_nucid[i][j];
      real dtmp=init_nucden[i][j];
      int idpos=med[0].SearchNuclide(idtmp);
      if(idpos==-1){
        cout<<"# Error !!\n";
        exit(0);
      };
      fwd_nuc[0][0][i].put_data(idpos,dtmp);
    };
  };

  // +++ Initial heavy metal weight calculation +++
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[0][0][i],0);
    hm_weight_init_per_medium[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*vol_med[i];
  };

  if(med_target!=-1){
    hm_weight_init=hm_weight_init_per_medium[med_target];
  }else{
    hm_weight_init=0.;
    for(int i=0;i<mednum_fuel;i++){
      hm_weight_init+=hm_weight_init_per_medium[i];
    };
  };

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"#\n# Initial heavy metal weight [g]\n#\n";
  cout<<"#     Total : "<<hm_weight_init<<"\n";
  for(int i=0;i<mednum_fuel;i++){
    cout<<"#       Medium "<<i<<" : "<<hm_weight_init_per_medium[i]<<"\n";
  };
  cout<<"#\n";

  if(input_power_unit=="MW_t"){
    input_power_unit="W_cm";
    for(int i=0;i<burn_step;i++){
      power_density_list[i]*=hm_weight_init;
    };
  };

  // +++ Burnup calculation condition setting +++
  PreCalculation_bt();

  // +++ Forward burn-up calculation +++
  real accumulated_day=0.;
  real accumulated_burn=0.;
  vector<real> accumulated_burn_per_medium(mednum_fuel,0.);

  // +++ Eigenvalue calculation
  MECSystem lat(group,mednum);
  lat.PutTrajectorySet(&sys_f);  

  if(adjoint){

    // +++ array setting for angular flux storing in adjoint calculations.
    int totsn=lat.GetQuad().GetSN();
    // Quadrature setting
    quad.Initialize(2,aflx_legendre);
    quad.PutSN(totsn);
    real *mui=new real[totsn];
    real *eai=new real[totsn];
    real *xii=new real[totsn];
    real *w=new real[totsn];
    for(int i=0;i<totsn;i++){
      mui[i]=lat.GetQuad().GetMu(i);
      eai[i]=lat.GetQuad().GetEata(i);
      xii[i]=lat.GetQuad().GetXi(i);
      w[i]=lat.GetQuad().GetOmega(i);
    };
    quad.PutData(mui,eai,xii,w);
    quad.CalValue();

    if(aflx_legendre==-1){
      volaflx_mesh.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };

  };

  // +++ GPT-PC +++++++++++++++++++++++++++++++++++++++++++++++++
  if(adjoint&&corrector_calc){
    fwd_nuc_p.resize(burn_step+1);
    fwd_nuc_p_int.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      int sub_step=sub_step_list[i];
      fwd_nuc_p[i].resize(sub_step+1);
      fwd_nuc_p_int[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
        fwd_nuc_p[i][j].resize(mednum_fuel);
        fwd_nuc_p_int[i][j].resize(mednum_fuel);
        for(int k=0;k<mednum_fuel;k++){
          fwd_nuc_p[i][j][k].put_imax(nucn);
          fwd_nuc_p_int[i][j][k].put_imax(nucn);
        };
      };
    };
    volflx_mesh_p.resize(burn_step);  
    power_factor_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      volflx_mesh_p[i].resize(totm);
      for(int j=0;j<totm;j++){
        if(region_medium[j]<mednum_fuel){
          volflx_mesh_p[i][j].put_imax(group);
        };
      };
      int sub_step=sub_step_list[i];
      power_factor_p[i].resize(sub_step+1);
    };
    int totsn=lat.GetQuad().GetSN();

    if(aflx_legendre==-1){
      volaflx_mesh_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_mesh_p[i].resize(totm);
        for(int j=0;j<totm;j++){
	  if(region_medium[j]<mednum_fuel){
	    volaflx_mesh_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_mesh_p[i][j][k].put_imax(totsn);
	    };
	  };
        };
      };
    }else{
      volaflx_pl_p.resize(burn_step);
      for(int i=0;i<burn_step;i++){
        volaflx_pl_p[i].resize(totm);
        for(int j=0;j<totm;j++){
    	  if(region_medium[j]<mednum_fuel){
	    volaflx_pl_p[i][j].resize(group);
	    for(int k=0;k<group;k++){
	      volaflx_pl_p[i][j][k].put_imax(quad.GetPlnum());
	    };
	  };
        };
      };
    };
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for(int st=0;st<burn_step+1;st++){

    int bstmp=0;
    if(adjoint)bstmp=st;

    acday.push_back(accumulated_day);
    acburn.push_back(accumulated_burn);
    for(int i=0;i<mednum_fuel;i++){
      acburn_per_medium[i].push_back(accumulated_burn_per_medium[i]);
    };

    cout<<"#\n# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    //lat.MedClear();
    for(int i=0;i<mednum_fuel;i++){
      // +++ Self-shielding calculation
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      opc.CalThermalScatteringMatrix(med[0],xslib,4.048);      
      med[0].CalMacroFromMicro();
      //opc.CalFissionSpectrumMatrix(med[0],xslib);
      //cout<<"# Total cross section of medium : "<<i<<"\n";
      //med[0].GetMacxs().GetData1d(sigt).show_self();
      //med[0].CalSigtr(0);

      // (Output of external file of medium instance)
      //string filename="med_med"+IntToString(i)+"_st"+IntToString(st);
      //med[0].WriteFile("./",filename,true);
      

      // +++ Macro&Micro data storing
      for(int k=0;k<nucn;k++){
        if(nuclide_info[k]!=0){
          if(nuclide_info[k]==1){
            mic_sigf[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
          };
          mic_sigc[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
          mic_sign2n[bstmp][i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
        };
      };
      if(st!=burn_step&&adjoint){
        macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
      };

      if(st==0){
        lat.AddMedium(med[0]); 
        lat.GetMedium(i).NuclideClear();
      }else{
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
      };
    };
    if(st==0){
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.AddMedium(med[1+jj]);
      };
    }else{
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };
    };
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);

    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    if(adjoint)lat.PutWriteFlux();
    lat.Print();
    keff[st]=lat.CalIgen();

    //cout<<"# Total inner iteration : "<<lat.GetTotalInnerIterationSum()<<"\n"; exit(0);

    // --- (Generating medium instances include collapsing cross section data) 
    if(collapsing){

    for(int i=0;i<mednum_fuel;i++){
      PutNuclideDataToMedium(fwd_nuc[st][0][i],0);
      opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
      if(dancoff_input){
        for(int g=0;g<group;g++){
          dancoff.put_data(g,1.-dancoff_factor[i]);
	};
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      }else{
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
      };
      //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      //med[0].CalMacroFromMicro();

      int ngrp=9;
      int bgrp[]={21,46,91,115,125,134,151,159,171};
      GroupData1D flx=lat.GetIntegratedFlux(i);
      Medium bmed=med[0].Cond(ngrp,bgrp,flx,flx,true);
      bmed.GetMacxs().GetData1d(d).copy(bmed.GetFlux());
      // Neutron flux is stored as diffusion coefficient (temporal treatment)
      string filename="bmed_med"+IntToString(i)+"_st"+IntToString(st);
      bmed.WriteFile(collapse_dir,filename,true);
    };

    };
    // --------------------------------------------------------------------------------


    // ------------------------------
    real sum_tot=0.;
    real sum_fmed=0.;
    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);

      real tmp=flx*lat.GetMedium(i).GetMacxs().GetData1d(siga);
      sum_tot+=tmp;
      if(i<mednum_fuel)sum_fmed+=tmp;

      flx_med[i]=flx*(1./vol_med[i]);
      // +++ One-group cross section storing (at predictor step)
      if(i<mednum_fuel){
        for(int j=0;j<nucn;j++){
          if(nuclide_info[j]!=0){
            if(nuclide_info[j]==1){
   	      xsf_1g[st][i][j]=mic_sigf[bstmp][i][j].Cond(flx);
	    };
	    xsc_1g[st][i][j]=mic_sigc[bstmp][i][j].Cond(flx);
	    xsn2n_1g[st][i][j]=mic_sign2n[bstmp][i][j].Cond(flx);
          };
	};
      };
    };
    abs_frac[st]=sum_fmed/sum_tot; // absorption rate fraction in fuel media
    cout<<"# Absorption fraction : "<<abs_frac[st]<<"\n";

    if(st!=burn_step){

      for(int i=0;i<totm;i++){
        if(region_medium[i]<mednum_fuel){
          real vol=lat.GetMesh(i).GetVolume();
          volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol; 
          if(adjoint){
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
  	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };
	  };
	};
      };

      // +++ Burnup calculation
      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;   

      cout<<"#... burnup calculation (total step:"<<sub_step<<")\n";

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med[med_target].get_sum();
          power_org=power_per_medium[med_target];
          //power_org=CalculationPinPower(bu,st,j,med_target,sumflx*vol_med[med_target]);
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        accumulated_day+=burn_span;
	//total_flux[st][j].resize(mednum_fuel);
        delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){
          total_flux[st][j][i]=flx_med[i].get_sum()*power_factor[st][j];
          CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	}; 

      }; // end of sub-step loop

      fwd_nuc[st+1][0]=fwd_nuc[st][sub_step];


      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //   CORRECTOR CALCULATION
      if(corrector_calc){

	cout<<"# Corrector calculation ...\n";

	// ... For A-OWPC
	for(int i=0;i<mednum_fuel;i++){
	  real tmp=total_flux[st][0][i];
  	  rr_gd5[i]=tmp*xsc_1g[st][i][pos_gd155]; // Reaction rate at BOC (Rp)
  	  rr_gd7[i]=tmp*xsc_1g[st][i][pos_gd157];
	};

	// ++++++++++++++++++++++++

	// (Predictor calculation results are stored in the array of [XXX_p].)
        swap(total_flux_p[st], total_flux[st]);
        swap(xsc_1g_p[st], xsc_1g[st]);
        swap(xsn2n_1g_p[st], xsn2n_1g[st]);
        swap(xsf_1g_p[st], xsf_1g[st]);

	if(adjoint){
          swap(keff_p[st], keff[st]);
          swap(fwd_nuc_p[st], fwd_nuc[st]);
          swap(fwd_nuc_p_int[st], fwd_nuc_int[st]);
          swap(power_factor_p[st], power_factor[st]);
	  swap(volflx_mesh_p[st], volflx_mesh[st]);
	  if(aflx_legendre==-1){
            swap(volaflx_mesh_p[st], volaflx_mesh[st]);
	  }else{
  	    swap(volaflx_pl_p[st], volaflx_pl[st]);
	  };
	  swap(macxs_p[st], macxs[st]);
	  fwd_nuc[st][0]=fwd_nuc_p[st][0];
	};

      for(int i=0;i<mednum_fuel;i++){
        // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
        PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        if(dancoff_input){
          for(int g=0;g<group;g++){
            dancoff.put_data(g,1.-dancoff_factor[i]);
  	  };
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        }else{
          opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        };
        //opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        //opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
        opc.CalThermalScatteringMatrix(med[0],xslib,4.048);	
        med[0].CalMacroFromMicro();
        //opc.CalFissionSpectrumMatrix(med[0],xslib);

        // +++ Macro&Micro data storing
        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            if(nuclide_info[k]==1){
              mic_sigf_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
            };
            mic_sigc_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
            mic_sign2n_c[i][k].copy(med[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sign2n));
          };
        };
	lat.GetMedium(i).GetMacxs()=med[0].GetMacxs();
	if(adjoint){
          macxs[st][i].DataCopyPL(med[0].GetMacxs(),0);
	};
      };

      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };

      lat.PutThermalIteration(3);
      lat.PutPL(0);
      lat.NoCMRAcceleration();
      real keff_corr=lat.CalIgen();
      
      // ------------------------------

      if(adjoint)keff[st]=keff_corr;

      for(int i=0;i<mednum;i++){
        GroupData1D flx=lat.GetIntegratedFlux(i);
        flx_med_c[i]=flx*(1./vol_med[i]);
        // +++ One-group cross section storing
        if(i<mednum_fuel){
          for(int j=0;j<nucn;j++){
            if(nuclide_info[j]!=0){
              if(nuclide_info[j]==1)xsf_1g[st][i][j]=mic_sigf_c[i][j].Cond(flx);
	      xsc_1g[st][i][j]=mic_sigc_c[i][j].Cond(flx);
	      xsn2n_1g[st][i][j]=mic_sign2n_c[i][j].Cond(flx);
	    };
          };
	};
      };

      if(adjoint){
        for(int i=0;i<totm;i++){
          if(region_medium[i]<mednum_fuel){
            real vol=lat.GetMesh(i).GetVolume();
            volflx_mesh[st][i]=lat.GetMesh(i).GetFlux()*vol;
	    if(aflx_legendre==-1){
              for(int g=0;g<group;g++){
	        volaflx_mesh[st][i][g]=lat.GetAFlux(i,g)*vol;
	      };
	    }else{
              for(int g=0;g<group;g++){
                int sntot=quad.GetSN();	      
	        for(int lm=0;lm<quad.GetPlnum();lm++){
                  real tmp=0.;
  	          for(int w=0;w<sntot;w++){
	  	    tmp+=quad.GetOmega(w)*quad.GetMoment(lm,w)*lat.GetAFlux(i,g).get_dat(w);
		  };
		  volaflx_pl[st][i][g].put_data(lm,tmp*vol);
	        };
	      };
	    };

	  };
        };


      };

      real acburn_pre=accumulated_burn-acburn[st];
      vector<real> acburn_pre_per_medium(mednum_fuel);
      for(int j=0;j<mednum_fuel;j++){
	acburn_pre_per_medium[j]=accumulated_burn_per_medium[j]-acburn_per_medium[j][st];
      };
      // accumulated burnup calculated by the predictor step

      for(int j=0;j<sub_step;j++){

	// (Line power of target medium is calculated)
	real sumflx=0.;
	real power_org=0.;

        vector<real> power_per_medium(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real tmp=flx_med_c[i].get_sum();
          power_per_medium[i]=CalculationPinPower(bu,st,j,i,tmp*vol_med[i]);
	};

	if(med_target!=-1){
  	  sumflx=flx_med_c[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
          power_org=0.;
	  for(int i=0;i<mednum_fuel;i++){
            power_org+=power_per_medium[i];
	  };
	};

	if(input_flux_level){
          power_factor[st][j]=flux_level_list[st]/sumflx;
	  accumulated_burn+=(power_org*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init*1e-6);
	}else{
  	  power_factor[st][j]=power_density/power_org;
          //accumulated_burn+=burn_time_gwd[st]/sub_step;
	};

        for(int i=0;i<mednum_fuel;i++){
	  accumulated_burn_per_medium[i]+=(power_per_medium[i]*power_factor[st][j]*1e-9)*burn_span/(hm_weight_init_per_medium[i]*1e-6);
	};

        //accumulated_day+=burn_span;
        //delt[st][j]=burn_span*24*60*60;
        for(int i=0;i<mednum_fuel;i++){
	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
	  //CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	  CalculationPinBurnup(bu,st,j,i,xsf_1g[st][i],xsc_1g[st][i],xsn2n_1g[st][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	};

      }; // end of sub-step


      // ... Variables required for OWPC)
      real dt=burn_time[st]*60*60*24;  // [day] -> [sec]
      real time_mesh=10000;            // The number of time meshes in the toy problem
      real dt_r=dt/time_mesh;
      real omega_155, omega_157; // various correlation conditions

      for(int i=0;i<mednum_fuel;i++){

        // ... OWPC treatment for Gd-155 and -157 with various correlation conditions
	if(owpc_corr){

	real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	real r0_155=rr_gd5[i]; //rp
	real r0_157=rr_gd7[i]; //rp
	real rp_155=xsc_1g[st][i][pos_gd155]*total_flux[st][0][i]; //rc
	real rp_157=xsc_1g[st][i][pos_gd157]*total_flux[st][0][i]; //rc
	real n0_r_155;
	real n0_r_157;

	real np_p_155=n0_155*exp(-r0_155*1e-24*dt); // Np in toy-problem
	real np_p_157=n0_157*exp(-r0_157*1e-24*dt);
	  
        // (Correlation to Gd-155 number density)
	real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	real alpha_157=(rp_157-r0_157)/(np_155-n0_155);
	//if(np_155>1e-4)alpha_157=(rp_157-r0_157)/(np_157-n0_157);
        // (Correlation to Gd-157 number density)
        //real alpha_155=(rp_155-r0_155)/(np_157-n0_157);
        //real alpha_157=(rp_157-r0_157)/(np_157-n0_157);

	// (Correlation to Gd-155 number density)
	real rc_155=r0_155+(np_p_155-n0_155)*alpha_155; // Rc in toy-problem
	real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	//if(np_155>1e-4)rc_157=r0_157+(np_p_157-n0_157)*alpha_157;
        // (Correlation to Gd-157 number density)
	//real rc_155=r0_155+(np_p_157-n0_157)*alpha_155;
        //real rc_157=r0_157+(np_p_157-n0_157)*alpha_157;

	n0_r_155=n0_155;
	real rp_r_155=r0_155;
	  
	n0_r_157=n0_157;
	real rp_r_157=r0_157;
	  
	for(int k=0;k<time_mesh;k++){

	  real np_r_155=n0_r_155*exp(-rp_r_155*1e-24*dt_r);
          real np_r_157=n0_r_157*exp(-rp_r_157*1e-24*dt_r);

          // (Correlation to Gd-155 number density)
	  real rc_r_155=r0_155+(np_r_155-n0_155)*alpha_155;
	  real rc_r_157=r0_157+(np_r_155-n0_155)*alpha_157;
	  //if(np_155>1e-4)rc_r_157=r0_157+(np_r_157-n0_157)*alpha_157;
          // (Correlation to Gd-157 number density)
          //real rc_r_155=r0_155+(np_r_157-n0_157)*alpha_155;
          //real rc_r_157=r0_157+(np_r_157-n0_157)*alpha_157;

	  real r_r_155=(rp_r_155+rc_r_155)*0.5;
	  real r_r_157=(rp_r_157+rc_r_157)*0.5;

	  n0_r_155=n0_r_155*exp(-r_r_155*1e-24*dt_r);
	  n0_r_157=n0_r_157*exp(-r_r_157*1e-24*dt_r);

	  rp_r_155=rc_r_155;
	  rp_r_157=rc_r_157;

	};

	real R_reference_155=(log(n0_155)-log(n0_r_155))/dt*1e24; // 1e24 is multiplied to get R in the unit of [burn]
	real R_reference_157=(log(n0_157)-log(n0_r_157))/dt*1e24;

        omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	};

	for(int j=0;j<nucn;j++){

	  real np=fwd_nuc[st+1][0][i].get_dat(j);      // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
	  real n_next=exp((log(np)+log(nc))*0.5);

          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  if(nucid==641550||nucid==641570){

  	    real n0=fwd_nuc[st][0][i].get_dat(j);
	    if((n0-np)/n0>1.e-2){

	      if(!owpc_corr){

	      // -- OWPC ----------------------- 

	      real r0=rr_gd5[i]; //Rp
	      //real r0mm=rr_gd4[i]*fwd_nuc[st][0][i].get_dat(pos_gd154); //N Rp
	      if(nucid==641570)r0=rr_gd7[i];
	      //if(nucid==641570)r0mm=rr_gd6[i]*fwd_nuc[st][0][i].get_dat(pos_gd156);

	      //real rp=xsc_1g[st][i][j]*total_flux_p[st][sub_step-1][i]; //Rc
	      real rp=xsc_1g[st][i][j]*total_flux[st][0][i]; //Rc
	      real alpha=(rp-r0)/(np-n0);

	      /*
	      //if(nucid==641550&&i==0&&st!=0){
		if(nucid==641570&&i==0&&st!=0){
		cout.setf(ios::scientific);
		cout.precision(5);
		cout<<"# OWPC\n";
		cout<<fwd_nuc[st-1][0][i].get_dat(j)<<" "<<xsc_1g_p[st-1][i][j]*total_flux_p[st-1][0][i]<<" ";
		//cout<<np_gd5[i]<<" "<<xsc_1g[st-1][i][j]*total_flux[st-1][0][i]<<" ";
		cout<<np_gd7[i]<<" "<<xsc_1g[st-1][i][j]*total_flux[st-1][0][i]<<" ";
		cout<<n0<<" "<<r0<<" ";
		cout<<np<<" "<<rp<<" ";
		cout<<"\n";
	      };
	      */

	      real np_p=n0*exp(-r0*1e-24*dt);   // Nc in toy-problem
	      real rc=alpha*np_p+(r0-n0*alpha); // Rc in toy-problem
	    
	      real n0_r=n0;
	      real rp_r=r0;
	      for(int k=0;k<time_mesh;k++){
	        real np_r=n0_r*exp(-rp_r*1e-24*dt_r);
	        real rc_r=alpha*np_r+(r0-n0*alpha);
	        real r_r=(rp_r+rc_r)/2;
	        n0_r=n0_r*exp(-r_r*1e-24*dt_r);
	        rp_r=rc_r;
	      };
              real R_reference=(log(n0)-log(n0_r))/dt*1e24;
	      real omega=(R_reference-r0)/(rc-r0);

	      // --------------------------------
	    	    
 	      // -- PPC -------------------------------
	      /*
  	      real R_p=-log(np/n0)/dt;
	      real R_c=-log(nc/n0)/dt;
  	      real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
  	      real R_reference=(R_p+R_c_c)/2;
	      real omega=(R_reference-R_p)/(R_c-R_p);
	      */
	      // --------------------------------
	    
  	      n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応

	      //cout<<"# OWPC factor is "<<omega<<" in medium "<<i<<" and nuclide "<<nucid<<"\n";

	      };
	    
              // OWPC treatment for Gd-155 and -157 with various correlation conditions
	      if(owpc_corr){
	      if(nucid==641550){
	        n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
	      }else if(nucid==641570){
	        n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
	      };
	      };

	    }else{

	      //cout<<"# OWPC is not adopted in medium "<<i<<" and nuclide "<<nucid<<"\n";

	    };

	  };

          //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
          //n_next=(np+nc)*0.5; // Conventional PC (in linear)
          //n_next=nc;
          //n_next=np;

          if(adjoint)        n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
          if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC

  	  // (For OWPC-2018)	
  	  if(nucid==641550)np_gd5[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
  	  if(nucid==641570)np_gd7[i]=fwd_nuc[st+1][0][i].get_dat(pos_gd157);

          fwd_nuc[st+1][0][i].put_data(j,n_next);
	};
      };

      // Adjustment of accumulated burn
      for(int i=0;i<mednum_fuel;i++){
        real acburn_cor=accumulated_burn_per_medium[i]-acburn_per_medium[i][st]-acburn_pre_per_medium[i];
        accumulated_burn_per_medium[i]=acburn_per_medium[i][st]+(acburn_pre_per_medium[i]+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };

      if(input_flux_level){
        real acburn_cor=accumulated_burn-acburn[st]-acburn_pre;
        accumulated_burn=acburn[st]+(acburn_pre+acburn_cor*wgt_nc)/(1.+wgt_nc);
      };
      cout<<"#     ... terminated.\n";
      }; // END OF CORRECTOR CALCULATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    }; // end part of [If (st!=burn_step) ]
 
  }; // loop-end of burnup step
};
