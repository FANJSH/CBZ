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
  //igi0=igia;
  //igi=igib;
  totm=igi.GetRegion();
  fuel_r=igi0.GetCircle(igi0.GetIcir()-1).GetR(); // [cm]

  // +++ TrajectorySet for self-shielding calculation ++++++++++++++++++++++++
  sys.PutBoundaryCondition(Periodic);
  sys.CalTrajectory(igi0,8,0.02,45.);
  //sys.CalTrajectory(igi0,1,0.08,45.);

  // +++ TrajectorySet for neutron flux calculation
  if(reflective){
    sys_f.PutBoundaryCondition(Reflective);
    sys_f.CalTrajectory(igi,64,0.02,360.);
  }else{
    sys_f.PutBoundaryCondition(Periodic);
    //sys_f.CalTrajectory(igi,1,0.08,45.);
    //sys_f.CalTrajectory(igi,5,0.03,45.);
    sys_f.CalTrajectory(igi,8,0.02,45.);
    //sys_f.CalTrajectory(igi,16,0.02,45.);
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    opc.CalThermalScatteringMatrix(med[i],xslib,3.93);
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

void MulticellBurner::ForwardCalculation(Burnup &bu, int med_target, bool adjoint)
{
  // Hard-coded parameters for weighted predictor-corrector
  //real wgt_nc=1.2; // relative weight for corrector
  real wgt_nc=1.0; // relative weight for corrector

  if(input_flux_level&&med_target==-1){
    cout<<"# Error in MulticellBurner::ForwardCalculation.\n";
    cout<<"# [med_target] should NOT be -1 if neutron flux level is posed.\n";
    exit(0);
  };

  GeneralOption opt;
  //opt.PutAdjointCal();
  //opt.PutEpsf(1e-30);
  //opt.PutOutitermax(1000);

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);

  // Dancoff correction printing
  /*
  cout<<"# Dancoff correction\n";
  for(int g=0;g<group;g++){
    cout<<g<<" ";  
    cout<<xslib.GetEnband().get_dat(g)<<" ";
    cout<<dancoff.get_dat(g)<<"\n";
  };
  exit(0);
  */

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
  vector< vector<GroupData1D> > mic_sigf_c;
  vector< vector<GroupData1D> > mic_sigc_c;
  vector< vector<GroupData1D> > mic_sign2n_c;
  vector<GroupData1D> flx_med_c;
  if(corrector_calc){
    int bstmp=1;
    if(adjoint)bstmp=burn_step+1;
    xsc_1g_p.resize(bstmp);
    xsn2n_1g_p.resize(bstmp);
    xsf_1g_p.resize(bstmp);
    for(int j=0;j<bstmp;j++){
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
    flx_med_c.resize(mednum);
    for(int i=0;i<mednum_fuel;i++){
      mic_sigf_c[i].resize(nucn);
      mic_sigc_c[i].resize(nucn);
      mic_sign2n_c[i].resize(nucn);
    };
    for(int i=0;i<mednum;i++){
      flx_med_c[i].put_imax(group);
    };
  };

  int bstmp=1;
  if(adjoint)bstmp=burn_step+1;
  mic_sigf.resize(bstmp); 
  mic_sigc.resize(bstmp); 
  mic_sign2n.resize(bstmp);
  for(int i=0;i<bstmp;i++){
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

  // (for OWPC-2018)
  vector<real> rr_gd5, rr_gd7; // Reaction rate at BOC
  if(corrector_calc){
    rr_gd5.resize(mednum_fuel);
    rr_gd7.resize(mednum_fuel);
  };

  // +++ Array setting for adjoint calculation +++
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
    //for(int k=0;k<nucn;k++){
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
      //};
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
  // (PJI)
  //PJISystem lat(group,mednum);
  // (MEC)
  MECSystem lat(group,mednum);

  lat.PutTrajectorySet(&sys_f);  


  // (PJI)
  // no processing
  // (MEC)
  if(adjoint){

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
    total_flux_p.resize(burn_step);
    power_factor_p.resize(burn_step);
    for(int i=0;i<burn_step;i++){
      volflx_mesh_p[i].resize(totm);
      for(int j=0;j<totm;j++){
        if(region_medium[j]<mednum_fuel){
          volflx_mesh_p[i][j].put_imax(group);
        };
      };
      int sub_step=sub_step_list[i];
      total_flux_p[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        total_flux_p[i][j].resize(mednum_fuel);
      };
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
      opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
      med[0].CalMacroFromMicro();
      //opc.CalFissionSpectrumMatrix(med[0],xslib);
      //cout<<"# Total cross section of medium : "<<i<<"\n";
      //med[0].GetMacxs().GetData1d(sigt).show_self();
      //med[0].CalSigtr(0);

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
      //lat.AddMedium(med[med_clad]);
      //lat.AddMedium(med[med_water]);
    }else{
      for(int jj=0;jj<mednum_nonfuel;jj++){
        lat.GetMedium(mednum_fuel+jj).GetMacxs()=med[1+jj].GetMacxs();
      };
      //lat.GetMedium(mednum_fuel).GetMacxs()=med[med_clad].GetMacxs();
      //lat.GetMedium(mednum_fuel+1).GetMacxs()=med[med_water].GetMacxs();
    };
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);

    // - (PJI)------------------------
    /*
    lat.PutSigmaCol(sigtr);
    lat.PutBuckling(0.);
    lat.PutPij();
    //keff[st]=lat.CalIgenPij(true); // fission spectrum matrix
    keff[st]=lat.CalIgenPij(); 
    if(adjoint&&pij_storing){
      for(int g=0;g<group;g++){
        pij_store[st][0][g].copy(lat.GetPij(g));
      };
    };
    */
    // ------------------------------

    //- (MEC) ------------------------

    lat.PutThermalIteration(3);
    lat.PutPL(0);
    lat.NoCMRAcceleration();
    //lat.NoTransportApprox();
    if(adjoint)lat.PutWriteFlux();
    lat.Print();
    keff[st]=lat.CalIgen();

    //cout<<"# Total inner iteration : "<<lat.GetTotalInnerIterationSum()<<"\n"; exit(0);

    // ------------------------------

    //vector<real> sumflx(mednum);

    for(int i=0;i<mednum;i++){
      GroupData1D flx=lat.GetIntegratedFlux(i);
      flx_med[i]=flx*(1./vol_med[i]);
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
	      /*
	      //cout<<"# Group "<<g<<"\n";
	      for(int w=0;w<sntot;w++){
		real tmp=0.;
                int pln=0;
		for(int l=0;l<=aflx_legendre;l++){
		  for(int m=0;m<=l;m++){
		    tmp+=(2.*l+1.)/PI4*quad.GetMoment(pln,w)*volaflx_pl[st][i][g].get_dat(pln);
		    pln++;
		  };
		};
                //cout<<w<<" "<<quad.GetMu(w)<<" "<<quad.GetEata(w)<<" "<<quad.GetXi(w)<<" "<<tmp<<" "<<lat.GetAFlux(i,g).get_dat(w)<<" "<<tmp/lat.GetAFlux(i,g).get_dat(w)<<"\n";
	      };
	      //cout<<"\n\n";
	      //if(g==group-1)exit(0);
	      */
	      // +++
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

  	// + OWPC  +++++++++++++++++++++;
	real span=burn_time[st]*24*60*60;
	int sub_step=sub_step_list[st];
	vector< vector< real > > omega_opt(mednum_fuel);
	for(int i=0;i<mednum_fuel;i++){
          real flxtmp=total_flux[st][0][i]; // at BOC
	  omega_opt[i].resize(nucn);
  	  for(int j=0;j<nucn;j++){
	    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	    if(nucid==641550)pos_gd155=j;
	    if(nucid==641570)pos_gd157=j;
  	    if(nucid==641550||nucid==641570){
	      real xsc1g=xsc_1g[st][i][j];
	      xsc1g=xsc1g*1e-24;
	      real rrr=xsc1g*flxtmp; //[cm2]*[1/cm2/s]
	      real total_reac=span+((exp(-rrr*span)-1)/rrr);
	      real slide_ratio=2*total_reac/(span*(1-exp(-rrr*span)));
	      omega_opt[i][j]=slide_ratio/(2-slide_ratio);
	      // cout<<"# omega_opt["<<i<<"]["<<j<<"]="<<omega_opt[i][j]<<"\n";
	    }else{
	      omega_opt[i][j]=wgt_nc;
	    };
	  };
	};

	// (For OWPC-2018)	
	for(int i=0;i<mednum_fuel;i++){
	  real tmp=total_flux[st][0][i];
	  for(int j=0;j<nucn;j++){
	    int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
  	    if(nucid==641550)rr_gd5[i]=tmp*xsc_1g[st][i][j];
  	    if(nucid==641570)rr_gd7[i]=tmp*xsc_1g[st][i][j];
	  };
	};

	// ++++++++++++++++++++++++

	if(adjoint){

          swap(keff_p[st], keff[st]);
          swap(fwd_nuc_p[st], fwd_nuc[st]);
          swap(fwd_nuc_p_int[st], fwd_nuc_int[st]);
          swap(total_flux_p[st], total_flux[st]);
          swap(power_factor_p[st], power_factor[st]);
	  swap(volflx_mesh_p[st], volflx_mesh[st]);
	  if(aflx_legendre==-1){
            swap(volaflx_mesh_p[st], volaflx_mesh[st]);
	  }else{
  	    swap(volaflx_pl_p[st], volaflx_pl[st]);
	  };
	  swap(macxs_p[st], macxs[st]);
	  //swap(dmdf_nuc_p[st], dmdf_nuc[st]);
          swap(xsc_1g_p[st], xsc_1g[st]);
          swap(xsn2n_1g_p[st], xsn2n_1g[st]);
          swap(xsf_1g_p[st], xsf_1g[st]);
	  fwd_nuc[st][0]=fwd_nuc_p[st][0];

	};

      for(int i=0;i<mednum_fuel;i++){
        // +++ Self-shielding calculation by the nuclide number densities (NND) at the END of burnup step
        PutNuclideDataToMedium(fwd_nuc[st+1][0][i],0); // fwd_nuc[st+1][0][i] (NND at the beginning of the next step is same as NND at the end of the present step)
        opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
        opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
        opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
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

      //lat.GetMedium(mednum_fuel).GetMacxs()=med[med_clad].GetMacxs();
      //lat.GetMedium(mednum_fuel+1).GetMacxs()=med[med_water].GetMacxs();

      //lat.PutRegMed(region_medium);
      //lat.PutGeneralOption(opt);
      //lat.NoPrint();

      // - (PJI) ---------------------------
      /*
      lat.PutSigmaCol(sigtr);
      lat.PutBuckling(0.);
      lat.PutPij();
      real keff_corr=lat.CalIgenPij();  // Eigenvalue calculation for corrector step
      if(adjoint&&pij_storing){
        for(int g=0;g<group;g++){
          pij_store[st][1][g].copy(lat.GetPij(g));
        };
      };
      */
      // -------------------

      //- (MEC) ------------------------

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
              if(nuclide_info[j]==1)xsf_1g[bstmp][i][j]=mic_sigf_c[i][j].Cond(flx);
	      xsc_1g[bstmp][i][j]=mic_sigc_c[i][j].Cond(flx);
	      xsn2n_1g[bstmp][i][j]=mic_sign2n_c[i][j].Cond(flx);
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
          bu.PutNucnum(nucn);
          for(int k=0;k<nucn;k++){
            int id=med[0].GetNuclideInTurn(k).GetMatnum();
            real den=fwd_nuc[st][j][i].get_dat(k);
            bu.PutNuclideData(k,id,den,xsf_1g[bstmp][i][k],xsc_1g[bstmp][i][k],xsn2n_1g[bstmp][i][k]);
          };
          power_per_medium[i]=bu.GetIntegratedPower(tmp*vol_med[i],true);
	};

	if(med_target!=-1){
  	  sumflx=flx_med_c[med_target].get_sum();
          power_org=power_per_medium[med_target];
	}else{
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
	  real tmp=total_flux[st][j][i];
	  total_flux[st][j][i]=flx_med_c[i].get_sum()*power_factor[st][j];
	  CalculationPinBurnup(bu,st,j,i,xsf_1g[bstmp][i],xsc_1g[bstmp][i],xsn2n_1g[bstmp][i],total_flux[st][j][i],delt[st][j],adjoint); // if [adjoint] is true, multistep calculation is done.
	  //if(!adjoint)total_flux[st][j][i]=(tmp+total_flux[st][j][i]*wgt_nc)/(1.+wgt_nc); // average total flux
	};

      }; // end of sub-step

      //ここから相関式
      /*
      // Calculation of NND at the end of the step by OWPC
      for(int i=0;i<mednum_fuel;i++){
	for(int j=0;j<nucn;j++){
	  //int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  real np=fwd_nuc[st+1][0][i].get_dat(j); // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results

	  // (For OWPC-2018)
	  real n_next;
	  real af_owpc2018;
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  if(nucid==641550||nucid==641570){  

	    real n0=fwd_nuc[st][0][i].get_dat(j);
	    real r0=rr_gd5[i];
	    if(nucid==641570)r0=rr_gd7[i];
	    real rp=xsc_1g[bstmp][i][j]*total_flux[st][0][i];
	    real alpha=(rp-r0)/(np-n0);
	    //cout<<i<<" "<<nucid<<" "<<alpha<<"\n";
	    af_owpc2018=-1.36392e-21*alpha+1.48821;
	    //cout<<"# Fa-new : "<<nucid<<" "<<alpha<<" "<<fa_new<<"\n";
	    real x=log10(-alpha);
	    real dt=burn_time[st]; // day

	    real a=0.000781757*dt+0.0206011;//dt仕様
	    real b=-0.0143447*dt+0.0791623;
	    if(nucid==641570){
	      a=0.000238469*dt+0.0657154;
	      b=-0.00262474*dt-0.854548;
	    };
	    /*
	    real a=0.00106012*dt+0.0146436;
	    real b=-0.0144283*dt+0.0882253;
	    if(nucid==641570){
	      a=0.000262386*dt+0.0648907;
	      b=-0.0031313*dt-0.836972;
	    };
	    // OLD 
	    real omega=a*x+b;


	    // +++ dynamic optimized weight calculation +++++++
	    /*
	    r0=r0*1.e-24; //rp
	    rp=rp*1.e-24; //rc

	    real n0_r=n0;
	    real rp_r=r0;

	    real r_r=0;
	    real rc_r=0;
	    real np_r=0;

            real tmp_np=n0*exp(-(dt*24*60*60)*r0);
	    real tmp_rp=r0+(tmp_np-n0)*alpha*1e-24;
	    real tmp_nc=n0*exp(-(dt*24*60*60)*tmp_rp);
	    
	    real time_mesh=1000;
	    real dt_r=(dt*24*60*60)/time_mesh;
	    for(int k=0;k<time_mesh;k++){
	      // n0_r*=exp(-rp_r*dt_r); // (predictor)

	      rp_r=r0+(np_r-n0)*alpha*1e-24;
	      np_r=n0_r*exp(-rp_r*dt_r);
	      //rc_r=alpha*np_r+(r0-n0*alpha);
	      rc_r=r0+(np_r-n0)*alpha*1e-24;
	      r_r=(rp_r+rc_r)/2;
	      n0_r=n0_r*exp(-r_r*dt_r);
	      //real nc_r=n0_r*exp(-rc_r*dt_r);
	      //n0_r=(np_r+nc_r)*0.5;
	      //rp_r=r_r;

	      //if(k==time_mesh-1)cout<<"           "<<r_r<<" "<<r0<<" "<<rp<<"\n";
	    };
	    cout.setf(ios::scientific);
	    cout.precision(5);
	    //cout<<n0<<" "<<np<<" "<<nc<<" "<<n0_r<<"\n";
	    //cout<<n0<<" "<<tmp_np<<" "<<tmp_nc<<" "<<n0_r<<"\n";
	    //real R_reference=(log(n0)-log(n0_r))/(dt*24*60*60);
	    //real omega_dyn=(R_reference-r0)/(tmp_rp-r0);

	    real omega_dyn=(log(n0_r)-log(tmp_np))/(log(tmp_nc)-log(tmp_np));
	    //omega=(log(n0_r)-log(tmp_np))/(log(tmp_nc)-log(tmp_np));

	    if(nucid==641550||(nucid==641570&&st<=5)){
              omega=omega_dyn;
	    };
	    

	    //real omega_dyn=(n0_r-np)/(nc-np);
	    //cout<<r0<<" "<<rp<<" "<<R_reference<<"\n";

	    if(nucid==641570&&i<8){
              //cout<<i<<" "<<nucid<<" "<<omega<<" "<<omega_dyn<<"\n";
	    };
	    // +++++++++++++++++++

	    //cout<<i<<" "<<nucid<<" "<<omega<<" "<<alpha<<"\n";

	    if(!(omega>0.&&omega<1.))omega=0.5;
	    n_next=exp(log(np)*(1-omega)+log(nc)*omega);

	  }else{
            n_next=exp((log(np)+log(nc))*0.5);
	  };

          //real n_next=(np+nc)*0.5(1.0+1.2)); // Weighted linear-PC
          //real n_next=exp((log(np)*1.0+log(nc)*1.2)/(1.0+1.2)); // Weighted log-PC
          //real n_next=exp((log(np)*1.0+log(nc)*wgt_nc)/(1.0+wgt_nc)); // Weighted log-PC
	  //opt//real n_next=exp((log(np)*1.0+log(nc)*omega_opt[st][i][j])/(1.0+omega_opt[st][i][j])); // Optimized-weighted log-PC
	  //optim//
	  */
      //---ここまで相関式
      
      //ここから究極---
      
      // Calculation of NND at the end of the step by OWPC
      for(int i=0;i<mednum_fuel;i++){
	
	real np_155=fwd_nuc[st+1][0][i].get_dat(pos_gd155);
	real np_157=fwd_nuc[st+1][0][i].get_dat(pos_gd157);
	real nc_155=fwd_nuc[st][sub_step][i].get_dat(pos_gd155);
	real nc_157=fwd_nuc[st][sub_step][i].get_dat(pos_gd157);
	real n0_155=fwd_nuc[st][0][i].get_dat(pos_gd155);
	real n0_157=fwd_nuc[st][0][i].get_dat(pos_gd157);
	real r0_155=rr_gd5[i]; //rpを意味している
	real r0_157=rr_gd7[i]; //rpを意味している
	real rp_155=xsc_1g[bstmp][i][pos_gd155]*total_flux[st][0][i]; //rcを意味している
	real rp_157=xsc_1g[bstmp][i][pos_gd157]*total_flux[st][0][i]; //rcを意味している
	real n0_r_155;
	real n0_r_157;
	real rc_155;
	real rc_157;
	real dt=burn_time[st]*60*60*24;
	real time_mesh=10000;
	real dt_r=dt/time_mesh;
	/*
	if(n0_157-np_157>n0_155-np_155){
	  /*	
	//155のみの相関
	real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	real alpha_157=(rp_157-r0_157)/(np_155-n0_155);
	  	
	  /*	  
	//157のみの相関
	real alpha_155=(rp_155-r0_155)/(np_157-n0_157);
	real alpha_157=(rp_157-r0_157)/(np_157-n0_157);
	  
	
	//それぞれの相関
	real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	real alpha_157=(rp_157-r0_157)/(np_157-n0_157);
	  	  
	r0_155=r0_155*1.e-24;//rp
	r0_157=r0_157*1.e-24;//rp
	  
	real np_p_155=n0_155*exp(-r0_155*dt);
	real np_p_157=n0_157*exp(-r0_157*dt);
	alpha_155=alpha_155*1.e-24;
	alpha_157=alpha_157*1.e-24;
	  	
	  /*
	//155のみ相関
	//real rc_155=alpha_155*np_p_155+(r0_155-n0_155*alpha_155);
	//real rc_157=alpha_157*np_p_155+(r0_157-n0_155*alpha_157);
	real rc_155=r0_155+(np_p_155-n0_155)*alpha_155;
	real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;
	
	/*
	//157のみ相関
	//rc_155=alpha_155*np_p_157+(r0_155-n0_157*alpha_155);
	//rc_157=alpha_157*np_p_157+(r0_157-n0_157*alpha_157);
	real rc_155=r0_155+(np_p_157-n0_157)*alpha_155;
	real rc_157=r0_157+(np_p_157-n0_157)*alpha_157;

	
	//それぞれの相関
	//rc_155=alpha_155*np_p_155+(r0_155-n0_155*alpha_155);
	//rc_157=alpha_157*np_p_157+(r0_157-n0_157*alpha_157);
	real rc_155=r0_155+(np_p_155-n0_155)*alpha_155;
	real rc_157=r0_157+(np_p_157-n0_157)*alpha_157;

	  
	n0_r_155=n0_155;
	real rp_r_155=r0_155;
	real r_r_155=0;
	real rc_r_155=0;
	real np_r_155=0;

	n0_r_157=n0_157;
	real rp_r_157=r0_157;
	real r_r_157=0;
	real rc_r_157=0;
	real np_r_157=0;
	  
	for(int k=0;k<time_mesh;k++){
	  /*
	  //155のみ相関
	  np_r_155=n0_r_155*exp(-rp_r_155*dt_r);
	  rc_r_155=alpha_155*np_r_155+(r0_155-n0_155*alpha_155);
	  rc_r_157=alpha_157*np_r_155+(r0_157-n0_155*alpha_157);
	  r_r_155=(rp_r_155+rc_r_155)/2;
	  r_r_157=(rp_r_157+rc_r_157)/2;
	  n0_r_155=n0_r_155*exp(-r_r_155*dt_r);
	  n0_r_157=n0_r_157*exp(-r_r_157*dt_r);
	  rp_r_155=rc_r_155;
	  rp_r_157=rc_r_157;
	  
	/*    
	  //157のみ相関
	  np_r_157=n0_r_157*exp(-rp_r_157*dt_r);
	  rc_r_155=alpha_155*np_r_157+(r0_155-n0_157*alpha_155);
	  rc_r_157=alpha_157*np_r_157+(r0_157-n0_157*alpha_157);
	  r_r_155=(rp_r_155+rc_r_155)/2;
	  r_r_157=(rp_r_157+rc_r_157)/2;
	  n0_r_155=n0_r_155*exp(-r_r_155*dt_r);
	  n0_r_157=n0_r_157*exp(-r_r_157*dt_r);
	  rp_r_155=rc_r_155;
	  rp_r_157=rc_r_157;
		  
	 
	  //それぞれの相関
	  np_r_157=n0_r_157*exp(-rp_r_157*dt_r);
	  np_r_155=n0_r_155*exp(-rp_r_155*dt_r);
	  rc_r_155=alpha_155*np_r_155+(r0_155-n0_155*alpha_155);
	  rc_r_157=alpha_157*np_r_157+(r0_157-n0_157*alpha_157);
	  r_r_155=(rp_r_155+rc_r_155)/2;
	  r_r_157=(rp_r_157+rc_r_157)/2;
	  n0_r_155=n0_r_155*exp(-r_r_155*dt_r);
	  n0_r_157=n0_r_157*exp(-r_r_157*dt_r);
	  rp_r_155=rc_r_155;
	  rp_r_157=rc_r_157;
	
	};
	
	}else{
	*/  
	  //155のみの相関
	  real alpha_155=(rp_155-r0_155)/(np_155-n0_155);
	  real alpha_157=(rp_157-r0_157)/(np_155-n0_155);
	
	  /* //157のみの相関
	     real alpha_155=(rp_155-r0_155)/(np_157-n0_157);
	     real alpha_157=(rp_157-r0_157)/(np_157-n0_157);
	  */
	r0_155=r0_155*1.e-24;//rp
	r0_157=r0_157*1.e-24;//rp

	  real np_p_155=n0_155*exp(-r0_155*dt);
	  real np_p_157=n0_157*exp(-r0_157*dt);
	  alpha_155=alpha_155*1.e-24;
	  alpha_157=alpha_157*1.e-24;
	  
	  //155のみ相関
	  //rc_155=alpha_155*np_p_155+(r0_155-n0_155*alpha_155);
	  //rc_157=alpha_157*np_p_155+(r0_157-n0_155*alpha_157);
	  real rc_155=r0_155+(np_p_155-n0_155)*alpha_155;
	  real rc_157=r0_157+(np_p_155-n0_155)*alpha_157;

	  /* //157のみ相関
	     //real rc_155=alpha_155*np_p_155+(r0_155-n0_157*alpha_155);
	     //real rc_157=alpha_157*np_p_155+(r0_157-n0_157*alpha_157);
	     real rc_155=r0_155+(np_p_157-n0_157)*alpha_155;
             real rc_157=r0_157+(np_p_157-n0_157)*alpha_157;
	  */
	  
	  n0_r_155=n0_155;
	  real rp_r_155=r0_155;
	  real r_r_155=0;
	  real rc_r_155=0;
	  real np_r_155=0;
	  
	  n0_r_157=n0_157;
	  real rp_r_157=r0_157;
	  real r_r_157=0;
	  real rc_r_157=0;
	  real np_r_157=0;
	  
	  for(int k=0;k<time_mesh;k++){
	    //155のみ相関
	    np_r_155=n0_r_155*exp(-rp_r_155*dt_r);
	    rc_r_155=alpha_155*np_r_155+(r0_155-n0_155*alpha_155);
	    rc_r_157=alpha_157*np_r_155+(r0_157-n0_155*alpha_157);
	    r_r_155=(rp_r_155+rc_r_155)/2;
	    r_r_157=(rp_r_157+rc_r_157)/2;
	    n0_r_155=n0_r_155*exp(-r_r_155*dt_r);
	    n0_r_157=n0_r_157*exp(-r_r_157*dt_r);
	    rp_r_155=rc_r_155;
	    rp_r_157=rc_r_157;
	  
	  
	    /* //157のみ相関
	       np_r_157=n0_r_157*exp(-rp_r_157*dt_r);
	       rc_r_155=alpha_155*np_r_157+(r0_155-n0_157*alpha_155);
	       rc_r_157=alpha_157*np_r_157+(r0_157-n0_157*alpha_157);
	       r_r_155=(rp_r_155+rc_r_155)/2;
	       r_r_157=(rp_r_157+rc_r_157)/2;
	       n0_r_155=n0_r_155*exp(-r_r_155*dt_r);
	       n0_r_157=n0_r_157*exp(-r_r_157*dt_r);
	       rp_r_155=rc_r_155;
	       rp_r_157=rc_r_157;
	    */
	  };
	
	  //};
	
	real R_reference_155=(log(n0_155)-log(n0_r_155))/dt;
	real R_reference_157=(log(n0_157)-log(n0_r_157))/dt;
        omega_155=(R_reference_155-r0_155)/(rc_155-r0_155);
	omega_157=(R_reference_157-r0_157)/(rc_157-r0_157);
	
	for(int j=0;j<nucn;j++){
	  //int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  real np=fwd_nuc[st+1][0][i].get_dat(j); // predictor results
  	  real nc=fwd_nuc[st][sub_step][i].get_dat(j); // corrector results
	  
	  // (For OWPC-2018)
	  real n_next;
	  real af_owpc2018;
          int nucid=med[0].GetNuclideInTurn(j).GetMatnum();
	  real n0=fwd_nuc[st][0][i].get_dat(j);
	  if(((n0-np)/n0>1.e-2)&&(nucid==641550||nucid==641570)){
	    //ここから流石変更
	    real r0=rr_gd5[i]; //rpを意味している
	    if(nucid==641570)r0=rr_gd7[i];
	    real rp=xsc_1g[bstmp][i][j]*total_flux[st][0][i]; //rcを意味している
	    

	    real alpha=(rp-r0)/(np-n0);
	    af_owpc2018=-1.36392e-21*alpha+1.48821;
	    real x=log10(-alpha);
	    
	    real dt=burn_time[st]*60*60*24; // second //流石変更
	    
	    r0=r0*1.e-24;//rp
	    real np_p=n0*exp(-r0*dt);
	    alpha=alpha*1.e-24;
	    real rc=alpha*np_p+(r0-n0*alpha);
	    
	    real n0_r=n0;
	    real rp_r=r0;
	    real time_mesh=10000;
	    real dt_r=dt/time_mesh;
	    real r_r=0;
	    real rc_r=0;
	    real np_r=0;
	    
	    for(int k=0;k<time_mesh;k++){
	      np_r=n0_r*exp(-rp_r*dt_r);
	      rc_r=alpha*np_r+(r0-n0*alpha);
	      r_r=(rp_r+rc_r)/2;
	      n0_r=n0_r*exp(-r_r*dt_r);
	      rp_r=rc_r;
	    };

            real R_reference=(log(n0)-log(n0_r))/dt;
	    real omega=(R_reference-r0)/(rc-r0);
	    
	    	    
	    //ppc---- 
	    real R_p=-log(np/n0)/dt;
	    real R_c=-log(nc/n0)/dt;
	    real R_c_c=(R_p-R_c)/(n0-np)*((np+nc)/2-np)+R_c;
	    R_reference=(R_p+R_c_c)/2;
	    omega=(R_reference-R_p)/(R_c-R_p);
	    //----ppc終わり
	    
	    
	    n_next=exp(log(np)*(1-omega)+log(nc)*omega);//ppc,OWPCの両方に対応

	    
	    //quadraticモデル参照--- Gdの同位体の相関を考慮するときはコメント外す
	    if(nucid==641550){
	      n_next=exp(log(np)*(1-omega_155)+log(nc)*omega_155);
	    }else if(nucid==641570){
	      n_next=exp(log(np)*(1-omega_157)+log(nc)*omega_157);
	    }
	    //---
	    
	    
	  }else{
	    n_next=exp((log(np)+log(nc))*0.5);
	  };
	     
          //real n_next=(np+nc)*0.5(1.0+1.2)); // Weighted linear-PC
          //real n_next=exp((log(np)*1.0+log(nc)*1.2)/(1.0+1.2)); // Weighted log-PC
          //real n_next=exp((log(np)*1.0+log(nc)*wgt_nc)/(1.0+wgt_nc)); // Weighted log-PC
	  //opt//real n_next=exp((log(np)*1.0+log(nc)*omega_opt[st][i][j])/(1.0+omega_opt[st][i][j])); // Optimized-weighted log-PC
	  //optim//
	  //---ここまで究極
	  
          // +++ OWPC 
	  /*
	  real af=2.8;
	  if(af_owpc2018>1.&&af_owpc2018<10.)af=af_owpc2018;
	  omega_opt[i][j]=1+(omega_opt[i][j]-1)*af;
	  */

          //n_next=exp((log(np)+log(nc))*0.5); // Conventional PC (in log)
          //n_next=(np+nc)*0.5; // Conventional PC (in linear)
	  //n_next=exp((log(np)*1.0+log(nc)*omega_opt[i][j])/(1.0+omega_opt[i][j])); // Original OWPC (in log)
          //n_next=nc;


      if(adjoint)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
      if(wpc_direct_calc)n_next=np*(1.-wc_gpt_wpc)+nc*wc_gpt_wpc; // linear-PC
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

  /*
  cout.setf(ios::scientific);
  cout.precision(5);
  //cout.precision(10);
  for(int i=0;i<burn_step+1;i++){
    cout<<" "<<i<<" "<<acday[i]<<" "<<acburn[i]<<" "<<keff[i]<<"\n";
  };
  */

  /*
  for(int i=0;i<nucn;i++){
    if(med[0].GetNuclideInTurn(i).GetMatnum()==962420){
      cout.setf(ios::scientific);
      cout.precision(10);
      cout<<"# Cm-242 number density : "<<fwd_nuc[burn_step][0][0].get_dat(i)<<"\n";
    };
  };
  */
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

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize, int med_target,int num,string *nucname_target_list,string filename)
{
  int *nucid_target_list=new int[num];
  for(int i=0;i<num;i++){
    nucid_target_list[i]=midt.ID(nucname_target_list[i]);
  };
  SensitivityCalculation(bu,med_normalize,med_target,num,nucid_target_list,filename);
  delete [] nucid_target_list;
};

void MulticellBurner::SensitivityCalculation(Burnup &bu,int med_normalize, int med_target,int num,int *nucid_target_list,string filename)
{
  if(corrector_calc){
    SensitivityCalculationPC(bu,med_normalize,med_target,num,nucid_target_list,filename);
    return;
  };

  bool isotropic_approx=false; // Isotropic approximation in perturbation integral

  int ssv=40;  // (sub-sub-time step division)

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
  //ForwardCalculation(bu,med_target,true);
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
    real response=fwd_nuc[burn_step][0][med_target].get_dat(nucid_turn);

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response ("<<midt.Name(nucid_target)<<" ND) : "<<response<<"\n";

  // +++ Adjoint calculation

  // (PJI)
  //PJISystem lat_gpt(group,mednum);
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

  // (PJI)
  /*
  lat_gpt.PutSigmaCol(sigtr);
  lat_gpt.PutBuckling(0.);
  lat_gpt.WriteProcOff();
  lat_gpt.PijTransformTrue();
  */
  // (MEC)

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
      if(id==nucid_target&&i==med_target){
	adj_nuc_e[i].put_data(j,1.0);
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
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
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
	  if(m==med_target){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc[st][j][m][g])*delt[st][j]*power_factor[st][j]/vol_med[m];
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

  sns.WriteFile("./",filename+"."+midt.Name(nucid_target)+".med"+IntToString(med_target));

  }; 
};

void MulticellBurner::SensitivityCalculationPC(Burnup &bu,int med_normalize, int med_target,int num,int *nucid_target_list,string filename)
{
  //corrector_calc=true;
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

  // PJI
  /*
  if(pij_storing){
    pij_store.resize(burn_step+1);
    for(int i=0;i<burn_step+1;i++){
      pij_store[i].resize(2);      
      pij_store[i][0].resize(group);
      pij_store[i][1].resize(group);
    };
  };
  */
  // MEC
  // no processing
  
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
    real response=fwd_nuc[burn_step][0][med_target].get_dat(nucid_turn);

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response ("<<midt.Name(nucid_target)<<" ND) : "<<response<<"\n";

  // +++ Adjoint calculation

  // (PJI)
  //PJISystem lat_gpt(group,mednum);
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

  // (PJI)
  /*
  lat_gpt.PutSigmaCol(sigtr);
  lat_gpt.PutBuckling(0.);
  lat_gpt.WriteProcOff();
  lat_gpt.PijTransformTrue();
  */
  // (MEC)

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
      if(id==nucid_target&&i==med_target){
	adj_nuc_e[i].put_data(j,1.0);
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
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
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
	  if(m==med_target){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc[st][j][m][g])*delt[st][j]*power_factor[st][j]/vol_med[m];
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
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_target];
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
	  if(m==med_target){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor_p[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc_p[st][j][m][g])*delt[st][j]*power_factor_p[st][j]/vol_med[m];
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

  sns.WriteFile("./",filename+"."+midt.Name(nucid_target)+".med"+IntToString(med_target));

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
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc[st][j][m][g])*delt[st][j]*power_factor[st][j]/vol_med[m];
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
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc_p[st][j][m][g])*delt[st][j]*power_factor_p[st][j]/vol_med[m];
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
    opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
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

void MulticellBurner::ShowNumberDensityReactionRate(int nucid)
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

  cout<<"#\n";
  cout<<"# A pair of ND [/cm3] and microscopic capture reaction rate [barn/cm2/s]";
  cout<<" of "<<midt.Name(nucid)<<" : ";
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
      real tmp=xsc_1g[i][j][target]*total_flux[i][0][j];
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
      for(int k=0;k<sz;k++){
        int tt=each_nuc_turn[j][k];
        int id=med[0].GetNuclideInTurn(tt).GetMatnum();
        //real dd=density_data[i][tt];
        real dd=fwd_nuc[i][0][iii].get_dat(tt);
        real dc=bu.GetBurnupChain().GetDecayConstant(id);
        real fuel_vol=vol_med[iii];
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
      sum+=val;
    };
    cout<<"  "<<sum<<"\n";
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

