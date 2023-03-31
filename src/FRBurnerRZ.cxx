#include <cstdlib>
#include "FRBurnerRZ.h"

#define calculation_theory 3 // 0:diffusion, 1:Tranpsort(SN), 2: SP3, 3:OSP3

FRBurnerRZ::FRBurnerRZ(int group_inp, int tnn, int *mninp)
{
  cout<<"# Current calculation theory: "<<calculation_theory<<" (0:diffusion, 1:Tranpsort(SN), 2: SP3, 3:OSP3). \n";
  dim    = 2;
  group  = group_inp;

  pl     = 1;
  sn     = 4;
  //pl     = 3;
  //sn     = 8;

  refuel_day=0.;
  cooling_day=0.;

  void_cal=false;
  dop_cal=false;
  xscalc_for_reactivity=true;

  print_linepower_map=false;
  cmfd_on=true;

  show_fission_info=false;
  ogawa_cal=false;

  dirname_boc="DEN_BOC";
  dirname_eoc="DEN_EOC";

  matno.resize(tnn);
  for(int i=0;i<tnn;i++){
    matno[i]=mninp[i];
  };

  //sodium_density_in_voided_case=1e-6;
  sodium_density_in_voided_case=1e-10;

  delta_t_in_doppler_case=534.+273.15;//MET-1000
  //delta_t_in_doppler_case=1027.+273.15; //MOX-1000
  //delta_t_in_doppler_case=1227.+273.15;//MOX-3600 
  //delta_t_in_doppler_case=987.+273.15;//CAR-3600
};

void FRBurnerRZ::PutNumberofMedium(int i)
{
  mednum=i;
  med.resize(mednum);
  med_voided.resize(mednum);
  med_doppler.resize(mednum);
  med_data_inp.resize(mednum,false);
};

void FRBurnerRZ::PutMedium(Medium &min,int mid)
{
  med[mid]=min;
  med_data_inp[mid]=true;
};

void FRBurnerRZ::PutCartMeshInfo(CartMeshInfo &cmii)
{
  cmi=cmii;
};

int FRBurnerRZ::GetBatchNumberFromMediumID(int i)
{
  return batch_mzone[mzone_per_med[i]];
};

void FRBurnerRZ::PreCalculation(XSLibrary &xslib, Burnup &bu)
{
  med_burn=0;
  for(int i=0;i<mednum;i++){
    if(mzone_per_med[i]!=-1)med_burn++;
  };

  for(int i=med_burn;i<mednum;i++){
    //for(int i=0;i<mednum;i++){
    if(!med_data_inp[i]){
      cout<<"# Error !!\n";
      cout<<"# Actual data for medium (ID="<<i<<") is not yet defined.\n";
      exit(0);
    };
  };

  vol_per_med.resize(mednum);
  thm_per_med.resize(mednum);

  // +++ Initial selfshielding calculation
  for(int i=med_burn;i<mednum;i++){
    opc.CalSelfShieldingInfiniteSystem(med[i],xslib);
  };

  density.resize(med_burn);
  density_ave.resize(med_burn);
  ir_cycle.resize(med_burn);
  ac_burn.resize(med_burn);
  for(int i=0;i<med_burn;i++){
    int tmp=GetBatchNumberFromMediumID(i);
    density[i].resize(tmp);
    ir_cycle[i].resize(tmp);
    ac_burn[i].resize(tmp);
    for(int j=0;j<tmp;j++){
      density[i][j].resize(burn_nuc,0.);
    };
    density_ave[i].resize(burn_nuc,0.);
  };

  sigc_1g.resize(med_burn);
  sigf_1g.resize(med_burn);
  sign2n_1g.resize(med_burn);
  for(int i=0;i<med_burn;i++){
    sigc_1g[i].resize(trace_cycle);
    sigf_1g[i].resize(trace_cycle);
    sign2n_1g[i].resize(trace_cycle);
    for(int j=0;j<trace_cycle;j++){
      sigc_1g[i][j].resize(cycle_div);
      sigf_1g[i][j].resize(cycle_div);
      sign2n_1g[i][j].resize(cycle_div);
      for(int k=0;k<cycle_div;k++){
        sigc_1g[i][j][k].resize(burn_nuc,0.);
        sigf_1g[i][j][k].resize(burn_nuc,0.);
        sign2n_1g[i][j][k].resize(burn_nuc,0.);
      };
    };
  };

  for(int i=0;i<num_mzone;i++){
    iden_mzone[i].resize(burn_nuc);
    for(int j=0;j<burn_nuc;j++){
      iden_mzone[i][j]=med[mzone_rep[i]].GetNuclideInTurn(j).GetDensity();
    };
  };
  
  InitializeNumberDensity();
  CalRegionAveragedDensity();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // External neutron source setting
  /*
  ifstream fin;
  fin.open("./SRC/tmp",ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  vector< vector<GroupData1D> > esrc(54);
  for(int iz=0;iz<54;iz++){
    esrc[iz].resize(38);
    for(int ir=0;ir<38;ir++){
      esrc[iz][ir].put_imax(group);
      esrc[iz][ir].set_zero();
      int tmp=0;
      while(tmp!=-1){
        fin>>tmp;
        if(tmp!=-1){
	  real val;
	  fin>>val;
	  esrc[iz][ir].put_data(tmp-1,val*1e-12);
	};
      };
    };
  };
  */
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  CalMacroFromMicroBurnupMedium();

  // Neutron flux calculation to store flux data in the Medium class
  // because this flux is used to calculate fission spectum vector
  GeneralOption opt;
  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    //cout<<"# Medium : "<<i<<"\n";
    //med[i].ShowMacroCrossSection1D();
    test.AddMedium(med[i]);
    test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();

#if 0
  // (fixed-source calculation)
  test.SetZeroScatSrc();
  for(int iz=0;iz<54;iz++){
    for(int ir=0;ir<38;ir++){
      test.PutIsotropicSourceParVolume(ir,ir,iz,iz,0,0,esrc[iz][ir]);
    };
  };
  real epsf=1e-4;
  int itermax=1000;
  bool high_speed_option=true;

  test.CalFixedSourceWithFission(epsf,itermax,high_speed_option);
#endif

  string optname="cmfd";
  if(!cmfd_on)optname="";
  real k1=test.CalIgen(optname);

  for(int i=0;i<mednum;i++){
    GroupData1D flx=test.GetIntegratedFlux(i);
    med[i].GetFlux().copy(flx);
    if(xscalc_for_reactivity&&i<med_burn){
      med_voided[i].GetFlux().copy(flx);
      med_doppler[i].GetFlux().copy(flx);
    };
  };

  // Total weight calculation

  real wgt=bu.CalWeightOfHeavyNuclide(test,0,mednum-1);
  real wgt_all=bu.CalWeightOfAllNuclide(test,0,mednum-1);
  cout<<"####################################################\n";
  cout<<"# Total weight of all   nuclide : "<<wgt_all*1e-6<<" [t]\n";
  cout<<"# Total weight of heavy nuclide : "<<wgt*1e-6<<" [t]\n";
  cout<<"####################################################\n";

  // Medium-wise weight calculation
  cout<<"# Medium ID   All[t]   Heavy[t]\n";
  for(int i=0;i<mednum;i++){
    real wgt=bu.CalWeightOfHeavyNuclide(test,i,i);
    real wgt_all=bu.CalWeightOfAllNuclide(test,i,i);
    cout<<"#     "<<i<<"      : "<<wgt_all*1e-6<<" "<<wgt*1e-6<<"\n";
  };
  cout<<"####################################################\n";
};

void FRBurnerRZ::InitializeNumberDensity()
{
  for(int i=0;i<mednum;i++){
    int mz=mzone_per_med[i];
    if(mz!=-1){
      for(int b=0;b<batch_mzone[mz];b++){
        for(int j=0;j<burn_nuc;j++){
          density[i][b][j]=iden_mzone[mz][j];
	};
      };
    };
  };
};

void FRBurnerRZ::Run(XSLibrary &xslib, Burnup &bu)
{
  bool heat_by_capture=true;
  bool heat_by_nonfuel=true;

  cout<<"###########################################################################\n";

  if(void_cal||dop_cal){
    cout<<"#  FRBurnerRZ : Calculation condition\n#\n";
    if(void_cal)cout<<"#       number density of sodium in voided case : "<<sodium_density_in_voided_case<<"\n";
    if(dop_cal) cout<<"#       delta temperature in Doppler case       : "<<delta_t_in_doppler_case<<"\n";
    cout<<"#########################################################################\n";
  };

  // (print option)
  bool print_flux_spectra_eoc=false;
  bool print_flux_spectra_boc=false;
  bool print_power_density=false;

  // (for linepower map printing)
  ofstream fout_pmap;
  if(print_linepower_map){
    fout_pmap.open("./linepower_map",ios::out);
    if(fout_pmap.fail()){
      cout<<"# Error in FRBurnerRZ.\n";
      cout<<"# Failed to open the file: linepower_map.\n";
      exit(0);
    };
  };
  int yf=cmi.GetYF();
  int xf=cmi.GetXF();
  vector<real> yf_edge(yf+1);
  vector<real> xf_edge(xf+1);
  xf_edge[0]=0.;
  yf_edge[0]=0.;

  real yy=0.;
  for(int i=0;i<yf;i++){
    yy+=cmi.GetFMeshL(1,i);
    yf_edge[i+1]=yy;
  };
  real xx=0.;
  for(int i=0;i<xf;i++){
    xx+=cmi.GetFMeshL(0,i);
    xf_edge[i+1]=xx;
  };
  int tot_page=trace_cycle*(cycle_div+1);
  fout_pmap<<"# Total number of pages : "<<tot_page<<"\n";

  vector< vector<real> > power_map(yf);
  for(int y=0;y<yf;y++){
    power_map[y].resize(xf);
  };

  // (initialize number of irradiated cycles and accumulated burnup)
  for(int i=0;i<med_burn;i++){
    int sz=ir_cycle[i].size();
    for(int b=0;b<sz;b++){
      ir_cycle[i][b]=0;
      ac_burn[i][b]=0.;
    };
  };

  // +++
  //density_store.resize(trace_cycle);
  fwd_nuc.resize(trace_cycle);  

  // +++
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  vector<int> refuel_id_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    refuel_id_mzone[i]=0;
  };
  /*
  vector<real> fis(burn_nuc);
  vector<real> cap(burn_nuc);
  vector<real> n2n(burn_nuc);
  */
  bu.PutNucnum(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    bu.PutNuclideData(i,med[0].GetNuclideInTurn(i).GetMatnum(),0.,0.,0.,0.);
  };
  bu.CalTransitionMatrixFluxInDependentPart();

  flux_level.resize(trace_cycle);
  keff.resize(trace_cycle);

  cout<<"#C Step  Day  Keff   Max. line   C.R.  Material zone-wise      ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  if(void_cal)cout<<"Void      ";
  if(dop_cal) cout<<"Doppler   ";
  cout<<"\n";
  cout<<"#y                   power[W/cm]       Power Dist. [%]         ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  if(void_cal)cout<<"React.    ";
  if(dop_cal) cout<<"React.    ";
  cout<<"\n";
  cout<<"#                    (pos:r,z)                              ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  if(void_cal)cout<<"[dk/kk']  ";
  if(dop_cal) cout<<"[dk/kk']  ";
  cout<<"\n";
  cout<<"#                                       ";
  for(int i=0;i<num_mzone;i++){
    cout<<"("<<i<<")  ";
  };
  cout<<"\n";

  // (cycle-loop)
  real ac_day=0.;
  for(int cycle=0;cycle<trace_cycle;cycle++){

    flux_level[cycle].resize(cycle_div+1);
    keff[cycle].resize(cycle_div+1);
    //density_store[cycle].resize(cycle_div+1);
    fwd_nuc[cycle].resize(cycle_div+1);

      // +++ Refueling
    if(cycle!=0){
      for(int i=0;i<med_burn;i++){
	int mz=mzone_per_med[i];
        int rid=refuel_id_mzone[mz];
	ir_cycle[i][rid]=0;
        ac_burn[i][rid]=0.;
	for(int j=0;j<burn_nuc;j++){
	  density[i][rid][j]=iden_mzone[mz][j];
	};
      };
      for(int i=0;i<num_mzone;i++){
	refuel_id_mzone[i]++;
	if(refuel_id_mzone[i]==batch_mzone[i])refuel_id_mzone[i]=0;
      };
    };

    // (in cycle-loop)
    for(int tim=0;tim<cycle_div+1;tim++){

      flux_level[cycle][tim].resize(mednum);
      real delt_d=cycle_length/cycle_div; // [d]
      real delt=delt_d*24*60*60; // [s]

      // +++ Calculation of region-averaged density
      CalRegionAveragedDensity();
      StoreRegionAveragedDensity(cycle,tim);

      // +++ Macroscopic cross section calculation
      CalMacroFromMicroBurnupMedium();

      // +++ Eigenvalue calculation ++++++++++++++++++++++++++++++++++++++
#if (calculation_theory !=1)
      // (DIFFUSION THEORY)
      PLOSSystem test(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
	//if(dop_cal)opc.CalSelfShieldingInfiniteSystem(med[i],xslib);
        test.AddMedium(med[i]);
	test.GetMedium(i).MicxsVectorClear2DData();
      };
      test.NoPrint();
      test.PutCartMeshInfo(cmi,"Cylinder");
      test.PutGeneralOption(opt);
      test.CalCoef();
#endif

#if calculation_theory == 1
      // (TRANSPORT THEORY)
      SNRZQuadrature quad(pl);
      quad.PutLevelSymmetric(sn);
      SNRZSystem test(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
	//if(dop_cal)opc.CalSelfShieldingInfiniteSystem(med[i],xslib);
        test.AddMedium(med[i]);
	test.GetMedium(i).MicxsVectorClear2DData();
      };
      test.NoPrint();
      test.PutPL(pl);
      test.SetQuadrature(&quad);
      test.PutCartMeshInfo(cmi,"Cylinder");
      test.PutGeneralOption(opt);
      test.SetArray();
#endif

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      string optname="cmfd";
      if(!cmfd_on)optname="";
      real k1=test.CalIgen(optname);
      //real k1=test.CalSP3(false,2.,true);
      //real k1=test.CalSP3();


      keff[cycle][tim]=k1;


      if(show_fission_info)test.GetNeutronMultiplicationInfo(k1);

      // ++++ Void condition +++++++++++++++++++++++++++++++++++++++++
      real kv=0.;
      if(void_cal){
        PLOSSystem test_v(dim,group,mednum); 
        for(int i=0;i<mednum;i++){
	  if(i<med_burn){
            int mz=mzone_per_med[i];
            int repid=mzone_rep[mz];
            for(int j=0;j<burn_nuc;j++){
  	      med_voided[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
            };
            if(med_voided[repid].ExistNuclide(110230)){
	      med_voided[repid].GetNuclide(110230).PutDensity(sodium_density_in_voided_case);
            };
            med_voided[repid].CalMacroFromMicro();
            if(i!=repid)med_voided[i].GetMacxs()=med_voided[repid].GetMacxs();
            test_v.AddMedium(med_voided[i]);
	  }else{
            test_v.AddMedium(med[i]);
	  };
  	  test_v.GetMedium(i).MicxsVectorClear2DData();
        };
        test_v.NoPrint();
        test_v.PutCartMeshInfo(cmi,"Cylinder");
        test_v.PutGeneralOption(opt);
        test_v.CalCoef();
        kv=test_v.CalIgen(optname);
      };

      // ++++ Doppler condition +++++++++++++++++++++++++++++++++++++++++
      real kd=0.;
      if(dop_cal){

        PLOSSystem test_v(dim,group,mednum); 
        for(int i=0;i<mednum;i++){
          if(i<med_burn){
            int mz=mzone_per_med[i];
            int repid=mzone_rep[mz];
            for(int j=0;j<burn_nuc;j++){
    	      med_doppler[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
            };
            med_doppler[repid].CalMacroFromMicro();
            if(i!=repid)med_doppler[i].GetMacxs()=med_doppler[repid].GetMacxs();
            test_v.AddMedium(med_doppler[i]);
          }else{
            test_v.AddMedium(med[i]);
          };
 	  test_v.GetMedium(i).MicxsVectorClear2DData();
        };
        test_v.NoPrint();
        test_v.PutCartMeshInfo(cmi,"Cylinder");
        test_v.PutGeneralOption(opt);
        test_v.CalCoef();
        kd=test_v.CalIgen("cmfd");

      };

      if(cycle==0&&tim==0){
	for(int i=0;i<mednum;i++){
	  vol_per_med[i]=test.GetVolumePerMedium(i);
          real weight=bu.CalWeightOfHeavyNuclideParUnitVolume(med[i]);
          thm_per_med[i]=weight*vol_per_med[i];
	};
      };

      // Neutron flux printing
      if((print_flux_spectra_boc&&tim==0)||(print_flux_spectra_eoc&&tim==cycle_div)){
	cout<<"\n+++ medium-wise neutron flux spectra (volume-integrated/per-volume) +++\n";
	for(int i=0;i<mednum;i++){
	  cout<<"\n#  Medium : "<<i<<"\n\n";
	  GroupData1D tmpe=med[0].GetEnband();
	  GroupData1D tmp0=test.GetIntegratedFlux(i);
	  GroupData1D tmp1=test.GetIntegratedFluxParVolume(i);
	  for(int g=0;g<group;g++){
	    cout<<tmpe.get_dat(g)<<"  "<<tmp0.get_dat(g)<<"  "<<tmp1.get_dat(g)<<"\n";
	  };
	};
      };

      /*
      real pow_sum=0.;
      for(int i=0;i<test.GetTotM();i++){
	real pow=bu.GetIntegratedPowerParMesh(test,i);
	pow_sum+=pow;
      };
      */

      for(int i=0;i<mednum;i++){
        GroupData1D flx=test.GetIntegratedFlux(i);
        med[i].GetFlux().copy(flx);
      };

      vector< vector<real> > pow_store(med_burn);
      // Volume-integrated power for each medium/batch
      // Note : "Volume" is defined as "total" one for all the batches

      // (power normalization)
      real pow_sum=0.;
      for(int i=mednum-1;i>=0;i--){
	//for(int i=0;i<mednum;i++){
        GroupData1D flx=med[i].GetFlux();
	if(i<med_burn){
          int mz=mzone_per_med[i];
	  int repid=mzone_rep[mz];
	  int sz=ac_burn[i].size(); // batch-wise
	  real pow=0.;
          pow_store[i].resize(sz);
	  for(int j=0;j<sz;j++){ 
	    for(int k=0;k<burn_nuc;k++){
	      med[repid].GetNuclideInTurn(k).PutDensity(density[i][j][k]);
	    };
            real tmp=bu.GetIntegratedPower(med[repid],flx,heat_by_capture);
            pow_store[i][j]=tmp;
	    //cout<<i<<" "<<pow_store[i][j]<<"\n";
	    pow+=tmp;
	  };
	  pow_sum+=pow/sz;
	}else if(heat_by_nonfuel){
	  pow_sum+=bu.GetIntegratedPower(med[i]);
	};
      };
      real factor=power/pow_sum;

      // +++ Power density calculation
      int yf=cmi.GetYF();
      int xf=cmi.GetXF();
      vector< vector<real> > power_map(yf);
      int index=0;
      for(int y=0;y<yf;y++){
	power_map[y].resize(xf);
	for(int x=0;x<xf;x++){
          int medid=cmi.GetFMat(index);
	  if(heat_by_nonfuel||medid<med_burn){
            int mz=mzone_per_med[medid];
  	    real pow=bu.GetIntegratedPowerParMesh(test,x,y,0,heat_by_capture);
 	    real vol=test.GetMesh(x,y,0).GetVolume();
	    //pow=pow/vol; // power density [power/volume]
	    real zl=test.GetMesh(x,y,0).GetLen(1);
	    real area=vol/zl;
            int pinn=num_pin_mzone[mz];
	    if(pinn==0)pinn=num_pin_mzone[0];
            //int pinn=num_pin;
            //if(medid>=rb_s&&medid<=rb_e)pinn=num_pin_rb;
	    //pow=(pow/zl)/(pinn*area/area_assembly); // line-power
	    pow=(pow/zl)/(pinn*area/area_assembly_mzone[mz]); // line-power
	    power_map[y][x]=pow*factor;
	  };
          index++;
	};
      };

      // +++ Maximum line power calculation
      int maxp_x,maxp_y;
      real maxpow=MaximumPowerDensityCalculation(power_map,maxp_x,maxp_y);
      if(print_linepower_map){
        cout<<"#--------------------------------------------------------------------\n";
        //cout<<"# Power density map in fine mesh [W/cm3]\n";
        cout<<"\n# Line power map in fine mesh [W/cm]\n";
        PrintPowerDensity(power_map);
      };

      if(print_linepower_map){
	fout_pmap<<"# New data\n";
	fout_pmap.setf(ios::scientific);
	fout_pmap.precision(5);
        for(int x=0;x<xf;x++){
	  //for(int x=0;x<=xmax;x++){
	  for(int jj=0;jj<2;jj++){
  	    //for(int y=ymin;y<=ymax;y++){
  	    for(int y=0;y<yf;y++){
	      real tmp=power_map[y][x];
	      fout_pmap<<xf_edge[x+jj]<<" "<<yf_edge[y]<<" "<<tmp<<"\n";
	      fout_pmap<<xf_edge[x+jj]<<" "<<yf_edge[y+1]<<" "<<tmp<<"\n";
	    };
	    fout_pmap<<"\n";
	  };
	};
	fout_pmap<<"\n\n";
      };
      //cout<<"   Max. line power : "<<maxpow*factor<<" (W/cm)\n";


      WriteOut(cycle,2);
      WriteOut(tim,4);
      int iac_day=int(ac_day);
      WriteOut(iac_day,6);
      cout<<" ";
      WriteOut(k1,"%7.5f");
      cout<<" ";
      WriteOut(maxpow,"%5.1f");
      cout<<" ("<<maxp_x<<","<<maxp_y<<")";

      // ++ Conversion ratio printing 
      real cr=test.CalConversionRatio();
      cout<<" ";
      WriteOut(cr,"%5.3f");

      // ++ Power distribution printing
      vector<real> pow_reg(num_mzone,0.);
      for(int i=0;i<med_burn;i++){
        int sz=ac_burn[i].size();
	for(int j=0;j<sz;j++){
  	  pow_reg[mzone_per_med[i]]+=pow_store[i][j]/sz;
	};
      };
      for(int i=0;i<num_mzone;i++){
	pow_reg[i]*=(100/pow_sum);
	WriteOut(pow_reg[i],"%5.1f");
      };

      cout.setf(ios::scientific);
      cout.precision(2);
      if(void_cal){
	cout<<" "<<1./k1-1./kv;
      };
      if(dop_cal){
	cout<<" "<<1./k1-1./kd;
      };

      cout<<"\n";

      for(int m=0;m<mednum;m++){
        real sumflx=med[m].GetFlux().get_sum()*factor/vol_per_med[m];
        flux_level[cycle][tim][m]=sumflx;
      };
      
      // Burnup calculation
      if(tim!=cycle_div){

	for(int i=0;i<med_burn;i++){
	  int sz=ac_burn[i].size();
	  for(int j=0;j<sz;j++){
	    ac_burn[i][j]+=((pow_store[i][j]*factor*1e-9)*delt_d)/(thm_per_med[i]*1e-6); // [GWd/t]
	  };
	};

        ac_day+=delt_d;
        for(int m=0;m<med_burn;m++){
          GroupData1D flx=med[m].GetFlux();
          int repid=mzone_rep[mzone_per_med[m]];
	  //bu.CalOneGroupCrossSection(med[m],flx);
	  bu.CalOneGroupCrossSection(med[repid],flx);
	  for(int i=0;i<burn_nuc;i++){
	    sigf_1g[m][cycle][tim][i]=bu.GetSigf(i);
            sigc_1g[m][cycle][tim][i]=bu.GetSigc(i);
            sign2n_1g[m][cycle][tim][i]=bu.GetSign2n(i);

	    /*
	    cout.setf(ios::scientific);
	    cout.precision(5);
	    int idd=med[0].GetNuclideInTurn(i).GetMatnum();
	    if(idd==270590)cout<<m<<" "<<sigc_1g[m][cycle][tim][i]<<"\n";
	    */

	  };
          bu.CalTransitionMatrixFluxDependentPart();

          int tmp=GetBatchNumberFromMediumID(m);
          for(int j=0;j<tmp;j++){
	    bu.PutDensity(density[m][j]);
            bu.BurnupCalculation(flux_level[cycle][tim][m],delt,"chebyshev",false);
	    density[m][j]=bu.GetDensity();
          };
	};

      }else{
        // +++ refueling
        bool decay_cal=false;
        real decay_cal_day=refuel_day;
	if(cycle!=trace_cycle-1&&decay_cal_day>0.){
	  decay_cal=true;
	}else{
	  if(cooling_day>0.){
	    decay_cal=true;
	    decay_cal_day=cooling_day;
	  };
	};
	if(decay_cal){
          ac_day+=refuel_day;
	  for(int m=0;m<med_burn;m++){
	    //cout<<"Decay Calc. : "<<m<<"/"<<med_burn<<"\n";
	    int sz=ac_burn[m].size();
	    for(int i=0;i<sz;i++){
	      for(int j=0;j<burn_nuc;j++){
	        bu.PutDensity(j,density[m][i][j]);
	      };

              int div=1;
	      if(decay_cal_day>100){
		div=20;
	      };
	      for(int j=0;j<div;j++){
	        bu.BurnupCalculationByKrylov(1e-15,decay_cal_day*24*60*60/div);
	      };

	      //bu.BurnupCalculationByKrylov(1e-15,decay_cal_day*24*60*60);
	      //bu.BurnupCalculation(0.,decay_cal_day*24*60*60,"mmpa",false);
	      for(int j=0;j<burn_nuc;j++){
	        density[m][i][j]=bu.GetDensity(j);
	      };
	    };
	  };
	};

      };

    }; // (incycle-loop end)

    // (increase the number of irradiated cycles)
    for(int i=0;i<med_burn;i++){
      int sz=ir_cycle[i].size();
      for(int j=0;j<sz;j++){
	ir_cycle[i][j]+=1;
      };
    };

  }; // (cycle-loop end)

  CalRegionAveragedDensity();
};

real FRBurnerRZ::MaximumPowerDensityCalculation(vector< vector<real> > &power_map, int &maxp_x, int &maxp_y)
{
  maxp_x=-1;
  maxp_y=-1;
  real maxpow=0.;
  int yf=power_map.size();
  int xf=power_map[0].size();
  for(int y=0;y<yf;y++){
    for(int x=0;x<xf;x++){
      if(power_map[y][x]>maxpow){
        maxpow=power_map[y][x];
        maxp_x=x;
        maxp_y=y;
      };
    };
  };
  return maxpow;
};

real FRBurnerRZ::CalTotalPower(Burnup &bu, vector< vector<real> > &pow_store, bool heat_by_nonfuel, bool heat_by_capture)
{
  real pow_sum=0.;
  int med_tmp=med_burn;
  if(heat_by_nonfuel)med_tmp=mednum;
  for(int i=0;i<med_tmp;i++){ 
    GroupData1D flx=med[i].GetFlux();
    if(i<med_burn){
      int sz=ac_burn[i].size(); // batch-wise
      real pow=0.;
      pow_store[i].resize(sz);
      for(int j=0;j<sz;j++){ 
        for(int k=0;k<burn_nuc;k++){
          med[i].GetNuclideInTurn(k).PutDensity(density[i][j][k]);
        };
        real tmp=bu.GetIntegratedPower(med[i],heat_by_capture); 
        pow_store[i][j]=tmp;
        pow+=tmp;
      };
      pow_sum+=pow/sz;
    }else{
      pow_sum+=bu.GetIntegratedPower(med[i]);
    };
  };
  return pow_sum;
};

void FRBurnerRZ::PrintPowerDensity(vector< vector<real> > &power_map)
{
  int yf=power_map.size();
  int xf=power_map[0].size();

  int xright=0;
  int xleft=xf;
  int ytop=-1;
  int ybottom=0;

  for(int y=0;y<yf;y++){
    bool all_zero=true;
    for(int x=0;x<xf;x++){
      if(int(power_map[y][x])>0.){
        all_zero=false;
        if(x>xright)xright=x;
        if(x<xleft)xleft=x;
      };
      if(!all_zero&&ytop==-1)ytop=y;
      if(!all_zero)ybottom=y;
    };
  };
  cout<<"#\n";
  //cout<<"#    cycle: "<<cycle<<" / step: "<<tim<<"\n";
  cout<<"#    Fine mesh range in x : "<<xleft<<" to "<<xright<<"\n";
  cout<<"#                    in y : "<<ytop<<" to "<<ybottom<<"\n";
  cout<<"#\n";
  for(int y=ytop;y<=ybottom;y++){
    for(int x=xleft;x<=xright;x++){
      int tmp=power_map[y][x];
      WriteOut(tmp,3);
      cout<<" ";
    };
    cout<<"\n";
  };
  cout<<"#--------------------------------------------------------------------\n";
};

void FRBurnerRZ::ExternalSourceReading(vector< vector<GroupData1D> > &esrc)
{
  ifstream fin;
  fin.open("./SRC/tmp",ios::in);
  if(fin.fail()){
    cout<<"# Error in FRBurnerRZ::ExternalSourceReading.\n";
    cout<<"# Failed to open the file.\n";
    exit(1);
  };

  esrc.resize(54);
  for(int iz=0;iz<54;iz++){
    esrc[iz].resize(38);
    for(int ir=0;ir<38;ir++){
      esrc[iz][ir].put_imax(group);
      esrc[iz][ir].set_zero();
      int tmp=0;
      while(tmp!=-1){
        fin>>tmp;
        if(tmp!=-1){
	  real val;
	  fin>>val;
	  esrc[iz][ir].put_data(tmp-1,val*1e-12); // [1e-12] is for unit converting
	};
      };
    };
  };
};

void FRBurnerRZ::ExternalSourceSetting(GeneralSystem &test, vector< vector<GroupData1D> > &esrc)
{
  test.SetZeroScatSrc();
  for(int iz=0;iz<54;iz++){
    for(int ir=0;ir<38;ir++){
      test.PutIsotropicSourceParVolume(ir,ir,iz,iz,0,0,esrc[iz][ir]);
    };
  };
};

void FRBurnerRZ::PutSAGeometryData(FRDTSubAssemblyGeometry &sa_geom, int mzone)
{
  if(mzone<0||mzone>=num_mzone){
    cout<<"# Error in FRBurnerRZ::PutSAGeometryData.\n";
    cout<<"# Material zone "<<mzone<<" is NOT defined.\n";
    exit(0);
  };

  num_pin_mzone[mzone]=sa_geom.GetNumberOfPin();
  area_assembly_mzone[mzone]=sa_geom.GetAssemblyVolume()*0.01; // 0.01 is to convert unit from mm2 to cm2 
};

void FRBurnerRZ::CalRegionAveragedDensity()
{
  for(int i=0;i<med_burn;i++){
    int btmp=GetBatchNumberFromMediumID(i);
    real btmp_inv=1./btmp;
    for(int j=0;j<burn_nuc;j++){
      real tmp=0.;
      for(int l=0;l<btmp;l++){
        tmp+=density[i][l][j];
      };
      tmp*=btmp_inv;
      med[i].GetNuclideInTurn(j).PutDensity(tmp);
      density_ave[i][j]=tmp;
    };
  };
};

void FRBurnerRZ::ShowPlotExternalSource(string filename)
{
  string filenametrue="./SRC/"+filename;

  ifstream fin;
  fin.open(filenametrue.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    exit(1);
  };

  vector< vector<GroupData1D> > esrc(54);
  for(int iz=0;iz<54;iz++){
    esrc[iz].resize(38);
    for(int ir=0;ir<38;ir++){
      esrc[iz][ir].put_imax(group);
      esrc[iz][ir].set_zero();
      int tmp=0;
      while(tmp!=-1){
        fin>>tmp;
        if(tmp!=-1){
	  real val;
	  fin>>val;
	  esrc[iz][ir].put_data(tmp-1,val*1e-12);
	};
      };
    };
  };

  GroupData2D esrcsum(54,38);
  for(int iz=0;iz<54;iz++){
    for(int ir=0;ir<38;ir++){
      esrcsum.put_data(iz,ir,esrc[iz][ir].get_sum());
    };
  };
  esrcsum.show_plot();

  fin.close();
};

void FRBurnerRZ::StoreRegionAveragedDensity(int cyc,int step)
{

  fwd_nuc[cyc][step].resize(med_burn);
  for(int i=0;i<med_burn;i++){
    fwd_nuc[cyc][step][i].put_imax(burn_nuc);
    for(int j=0;j<burn_nuc;j++){
      real den=med[i].GetNuclideInTurn(j).GetDensity();
      fwd_nuc[cyc][step][i].put_data(j,den);
    };
  };
  /*
  fwd_nuc[cyc][step].resize(med_burn);
  for(int i=0;i<med_burn;i++){
    fwd_nuc[cyc][step][i].put_imax(burn_nuc);
    int btmp=GetBatchNumberFromMediumID(i);
    real btmp_inv=1./btmp;
    for(int j=0;j<burn_nuc;j++){
      real tmp=0.;
      for(int l=0;l<btmp;l++){
        tmp+=density[i][l][j];
      };
      tmp*=btmp_inv;
      fwd_nuc[cyc][step][i].put_data(j,tmp);
    };
  };
  */

};


/*
void FRBurnerRZ::CalNeutronFluxEnergySpectrum()
{
  CalRegionAveragedDensity();

  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicroSimple();
  };

  GeneralOption opt;

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  real k2=test.CalIgen("cmfd");

  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<mednum;i++){
    cout<<"#\n# Medium-wise neutron flux energy spectrum [/lethargy]\n";
    cout<<"#      (Medium ID : "<<i<<")\n#\n";
    cout<<"# Upper Eng.   Flux\n";
    for(int g=0;g<group;g++){
      real e0=med[0].GetEnband().get_dat(g);
      real e1=med[0].GetEnband().get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<"  "<<e0<<"   ";
      cout<<test.GetIntegratedFlux(i).get_dat(g)/letwid<<"\n";
    };
    cout<<"\n\n";
  };
};
*/

// +++ post-burnup calculation

void FRBurnerRZ::PutNumberDensity(int cyc,int step)
{  
  for(int i=0;i<med_burn;i++){
    for(int j=0;j<burn_nuc;j++){
      real den=fwd_nuc[cyc][step][i].get_dat(j);
      med[i].GetNuclideInTurn(j).PutDensity(den);
      density_ave[i][j]=den;
    };
  };
};

/*
void FRBurnerRZ::EigenvalueCalculation()
{
};
*/

void FRBurnerRZ::CalNeutronFluxEnergySpectrum(bool adjoint,int meds,int mede)
{
  if(mede==-1)mede=meds;
  if(meds==-1){
    meds=0;
    mede=mednum;
  };

  // Macroscopic cross section calculations
  /*
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };
  */
  CalMacroFromMicroBurnupMedium();

  GeneralOption opt;
  if(adjoint)opt.PutAdjointCal();

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
    test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on||adjoint)optname="";
  real k1=test.CalIgen(optname);

  cout.setf(ios::scientific);
  cout.precision(4);

  for(int i=meds;i<=mede;i++){
    cout<<"#\n# Medium-wise neutron flux energy spectrum";
    if(adjoint) cout<<" (adjoint) ";
    if(!adjoint)cout<<" [/lethargy]";
    cout<<"\n";
    cout<<"#      (Medium ID : "<<i<<")\n#\n";
    cout<<"# Upper Eng.   Flux\n";
    for(int g=0;g<group;g++){
      real e0=med[0].GetEnband().get_dat(g);
      real e1=med[0].GetEnband().get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<"  "<<e0<<"   ";
      real tmp=test.GetIntegratedFlux(i).get_dat(g);
      if(!adjoint)tmp/=letwid;
      cout<<tmp<<"\n";
    };
    cout<<"\n\n";
  };

};

void FRBurnerRZ::Cal1groupXS()
{
  // Macroscopic cross section calculations
  /*
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };
  */
  CalMacroFromMicroBurnupMedium();

  GeneralOption opt;

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
    test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  real k1=test.CalIgen("cmfd");

  cout<<"#-------------------------------------------------------------\n";
  cout<<"#\n# MEDIUM-WISE ONE GROUP CROSS SECTION LIST\n#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<mednum;i++){
    cout<<"#-------------------------------------------------------------\n";
    cout<<"#\n# Medium : "<<i<<"\n#\n";
    cout<<"#          FISSION     CAPTURE     (N,2N)\n";
    int sz=med[i].GetNucnum();
    GroupData1D flx=test.GetIntegratedFlux(i);
    real flxsum_inv=1./flx.get_sum();
    for(int j=0;j<sz;j++){
      int id=med[i].GetNuclideInTurn(j).GetMatnum();
      real sigf_1g=med[i].GetNuclideInTurn(j).GetMicxs().GetData1d(sigf)*flx*flxsum_inv;
      real sigc_1g=med[i].GetNuclideInTurn(j).GetMicxs().GetData1d(sigc)*flx*flxsum_inv;
      real sign2n_1g=med[i].GetNuclideInTurn(j).GetMicxs().GetData1d(sign2n)*flx*flxsum_inv;
      cout<<"  ";
      WriteOut(id,7);
      cout<<"  "<<sigf_1g<<"  "<<sigc_1g<<"  "<<sign2n_1g<<"\n";
    };
  };
  cout<<"#-------------------------------------------------------------\n";

};

void FRBurnerRZ::ShowNeutronMultiplicationInfo()
{
  CalMacroFromMicroBurnupMedium();

  GeneralOption opt;

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
    test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real k1=test.CalIgen(optname);
  test.GetNeutronMultiplicationInfo(k1);
};

SensitivityData FRBurnerRZ::CalKeffSensitivity()
{
  //CalMacroFromMicroBurnupMedium();

  // The following is modified version of [CalMacroFromMicroBurnupMedium]
  // This is because the microscopic cross section data are required for all the media.
  for(int i=0;i<num_mzone;i++){
    int id=mzone_rep[i]; // ID for representative medium instance containing microscopic XS
    for(int jj=med_burn-1;jj>=0;jj--){
      if(mzone_per_med[jj]==i){
        for(int j=0;j<burn_nuc;j++){
          med[id].GetNuclideInTurn(j).PutDensity(density_ave[jj][j]);
	};
	med[id].CalMacroFromMicro();
	if(jj!=id){
  	  //med[jj].GetMacxs()=med[id].GetMacxs();  // !!!
  	  med[jj]=med[id]; // !!!
	};
      };
    };
  };

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
    //test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real k1=test.CalIgen(optname);

  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    //testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();
  real ka=testa.CalIgen();

  /*
  int nuc=med[0].GetNucnum();
  int *nucid=new int[nuc];
  for(int i=0;i<nuc;i++){
    nucid[i]=med[0].GetNuclideInTurn(i).GetMatnum();
  };
  */

  int nuc=med[0].GetNucnum();
  vector<int> nucid_tmp;
  for(int i=0;i<nuc;i++){
    int tmp=med[0].GetNuclideInTurn(i).GetMatnum();
    int grp=med[0].GetNuclideInTurn(i).GetGrp();
    if(grp!=-1){
      nucid_tmp.push_back(tmp);
    };
  };

  int sz=nucid_tmp.size();
  int *nucid=new int[sz];
  for(int i=0;i<sz;i++){
    nucid[i]=nucid_tmp[i];
  };
  
  SensitivityData sns=testa.CalSensitivityNew(&test,ka,sz,nucid);
  delete [] nucid;
  return sns;
};

SensitivityData FRBurnerRZ::CalKeffSensitivity(XSLibrary &xslib)
{
  for(int i=0;i<num_mzone;i++){
    int id=mzone_rep[i]; // ID for representative medium instance containing microscopic XS
    for(int jj=med_burn-1;jj>=0;jj--){
      if(mzone_per_med[jj]==i){
        for(int j=0;j<burn_nuc;j++){
          med[id].GetNuclideInTurn(j).PutDensity(density_ave[jj][j]);
	};
	//med[id].CalMacroFromMicro();
	if(jj!=id){
  	  //med[jj].GetMacxs()=med[id].GetMacxs();  // !!!
  	  med[jj]=med[id]; // !!!
	};
      };
    };
  };

  // Fe-Nat. are decomposed into isotope-wise data
  int matin[]={260540,260560,260570,260580};
  real ratio[]={0.05845, 0.91752, 0.02119, 0.00282};
  
  vector<Medium> mednew(mednum);
  for(int i=0;i<mednum;i++){
    mednew[i].PutImax(group);
    int plinp=pl;
    if(plinp==0)plinp=1;
    mednew[i].PutPL(plinp);
    int nucnum=med[i].GetNucnum();
    for(int j=0;j<nucnum;j++){
      int id=med[i].GetNuclideInTurn(j).GetMatnum();
      real den=med[i].GetNuclideInTurn(j).GetDensity();
      real temp=med[i].GetNuclideInTurn(j).GetTemperature();
      if(id!=260000){
        Nuclide nucnew;
        nucnew.PutMatnum(id);
        nucnew.PutDensity(den);
        nucnew.PutTemperature(temp);
        mednew[i].AddNuclide(nucnew);
      }else{
	for(int k=0;k<4;k++){
          Nuclide nucnew;
          nucnew.PutMatnum(matin[k]);
          nucnew.PutDensity(den*ratio[k]);
          nucnew.PutTemperature(temp);
          mednew[i].AddNuclide(nucnew);
	};
      };
    };

  };

  /*
  for(int i=0;i<mednum;i++){
    cout<<"# Medium : "<<i<<"\n";
    int nucnum=med[i].GetNucnum();
    cout<<"#   # of nuclides : "<<nucnum<<"\n";
    for(int j=0;j<nucnum;j++){
      cout<<j<<" "<<med[i].GetNuclideInTurn(j).GetMatnum()<<" "<<med[i].GetNuclideInTurn(j).GetDensity()<<" "<<med[i].GetNuclideInTurn(j).GetTemperature()<<"\n";
    };
  };
  exit(0);
  */

  OnePointCalculator opc;  
  for(int i=0;i<mednum;i++){
    opc.CalSelfShieldingInfiniteSystem(mednew[i],xslib);
    mednew[i].GetMacxs().GetData1d(chi)=med[i].GetMacxs().GetData1d(chi);    
  };


  GeneralOption opt,opta;
  opta.PutAdjointCal();

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(mednew[i]);
    //test.GetMedium(i).MicxsVectorClear2DData();
  };
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real k1=test.CalIgen(optname);

  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(mednew[i]);
    //testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();
  real ka=testa.CalIgen();

  int nuc=mednew[0].GetNucnum();
  vector<int> nucid_tmp;
  for(int i=0;i<nuc;i++){
    int tmp=mednew[0].GetNuclideInTurn(i).GetMatnum();
    int grp=mednew[0].GetNuclideInTurn(i).GetGrp();
    if(grp!=-1){
      nucid_tmp.push_back(tmp);
    };
  };

  int sz=nucid_tmp.size();
  int *nucid=new int[sz];
  for(int i=0;i<sz;i++){
    nucid[i]=nucid_tmp[i];
  };
  SensitivityData sns=testa.CalSensitivityNew(&test,ka,sz,nucid);
  delete [] nucid;
  return sns;
};


real FRBurnerRZ::CalVoidReactivity(int meds, int mede)
{
  if(mede==-1)mede=meds;
  if(meds==-1){
    meds=0;
    mede=med_burn-1;
  };
  
  int vmednum=mede-meds+1;
  int *vmedlist=new int[vmednum];
  for(int i=0;i<vmednum;i++){
    vmedlist[i]=meds+i;
  };
  real rho=CalVoidReactivity(vmednum,vmedlist);
  delete [] vmedlist;
  
#if 0
  if(!xscalc_for_reactivity){
    cout<<"# Error in FRBurnerRZ::CalVoidReactivity.\n";
    cout<<"# Cross section in voided condition is NOT calculated.\n";
    exit(0);
  };

  CalMacroFromMicroBurnupMedium();

  if(mede==-1)mede=meds;
  if(meds==-1){
    meds=0;
    mede=med_burn-1;
  };

  GeneralOption opt,opta;
  opta.PutAdjointCal();
#endif
  
  // --- (DIFFUSION/SP3 THEORY) ------------------------------
#if 0
  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();

  PLOSSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(i>=meds&&i<=mede){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	med_voided[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      if(med_voided[repid].ExistNuclide(110230)){
	med_voided[repid].GetNuclide(110230).PutDensity(sodium_density_in_voided_case);
      };
      med_voided[repid].CalMacroFromMicro();
      if(i!=repid)med_voided[i].GetMacxs()=med_voided[repid].GetMacxs();
      test_v.AddMedium(med_voided[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.CalCoef();

  // (diffusion)

  real k1=testa.CalIgen();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);
  real rho=testa.CalReactivity(&test_v,k1,kv);


  // (SP3)

  PLOSSystem testa_p2,test_v_p2;
  real k1=testa.CalSP3(testa_p2);
  real kv=test_v.CalSP3(test_v_p2);
  testa.CalReactivitySP3(&testa_p2,&test_v,&test_v_p2,k1,kv);

#endif

  // -----------------------------------------------------

  // --- (TRANSPORT THEORY) --------------------------------
#if 0
  SNRZQuadrature quad(pl);
  quad.PutLevelSymmetric(sn);
  SNRZSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutPL(pl);
  testa.SetQuadrature(&quad);
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.PutWriteFlux();
  testa.SetArray();
  real k1=testa.CalIgen();

  SNRZSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(i>=meds&&i<=mede){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	med_voided[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      if(med_voided[repid].ExistNuclide(110230)){
	real org=med_voided[repid].GetNuclide(110230).GetDensity();
	med_voided[repid].GetNuclide(110230).PutDensity(org*void_ratio);
      };
      med_voided[repid].CalMacroFromMicro();
      if(i!=repid)med_voided[i].GetMacxs()=med_voided[repid].GetMacxs();
      test_v.AddMedium(med_voided[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutPL(pl);
  test_v.SetQuadrature(&quad);
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.PutWriteFlux();
  test_v.SetArray();
  real kv=test_v.CalIgen();

  real rho=testa.CalReactivity(&test_v,k1,kv);
#endif
  // -----------------------------------------------------
  return rho;
};


real FRBurnerRZ::CalVoidReactivity(int vmednum, int *vmedlist)
{
  using namespace std::chrono;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  // ... Checking voiding medium
  for(int i=0;i<vmednum;i++){
    if(vmedlist[i]>=med_burn){
      cout<<"# Error in FRBurnerRZ::CalVoidReactivity.\n";
      cout<<"# Medium "<<vmedlist[i]<<" cannot be voided at present implementation.\n";
      exit(0);
    };
  };

  vector<bool> voided_medium(mednum);
  for(int i=0;i<mednum;i++){
    voided_medium[i]=false;
    for(int j=0;j<vmednum;j++){
      if(vmedlist[j]==i)voided_medium[i]=true;
    };
  };

  if(!xscalc_for_reactivity){
    cout<<"# Error in FRBurnerRZ::CalVoidReactivity.\n";
    cout<<"# Cross section in voided condition is NOT calculated.\n";
    exit(0);
  };

  CalMacroFromMicroBurnupMedium();

  GeneralOption opt,opta;
  opta.PutAdjointCal();
cout<<"# 3. Adjoint calculation actived. "<<endl;
// --- (DIFFUSION/SP3 THEORY) ------------------------------
#if (calculation_theory != 1)
  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();

  PLOSSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(voided_medium[i]){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	      med_voided[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      if(med_voided[repid].ExistNuclide(110230)){
	      med_voided[repid].GetNuclide(110230).PutDensity(sodium_density_in_voided_case);
      };
      med_voided[repid].CalMacroFromMicro();
      if(i!=repid)med_voided[i].GetMacxs()=med_voided[repid].GetMacxs();
      test_v.AddMedium(med_voided[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.CalCoef();
#endif 
#if (calculation_theory == 0)
cout<<"# 5. Current calculation method: diffusion."<<endl;
  // (diffusion)
  real k1=testa.CalIgen();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);
  real rho=testa.CalReactivity(&test_v,k1,kv);
#endif

#if(calculation_theory == 2)
cout<<"# 5. Current calculation method: SP3."<<endl;
  // (SP3)
  PLOSSystem testa_p2,test_v_p2;
  real k1=testa.CalSP3(testa_p2);
  real kv=test_v.CalSP3(test_v_p2);
  real rho = testa.CalReactivitySP3(&testa_p2,&test_v,&test_v_p2,k1,kv);
#endif

#if(calculation_theory == 3)
cout<<"# 5. Current calculation method: OSP3."<<endl;
  // (OSP3)
  PLOSSystem testa_p2,test_v_p2;
  cout<<"# Creating two PLOS system instances for phi-2 equation."<<endl;
  real k1=testa.CalOSP3Adjoint(testa_p2);
  cout<<"# Calculate keff of adjoint phi-2."<<endl;
  real kv=test_v.CalSP3(test_v_p2);
  cout<<"# Calculate keff of void system forward equation."<<endl;
  real rho=testa.CalReactivityOSP3(&testa_p2,&test_v,&test_v_p2,k1,kv);
  cout<<"# OSP3 preparation cal finished."<<endl;
#endif

  // -----------------------------------------------------

  // --- (TRANSPORT THEORY) --------------------------------
#if (calculation_theory == 1)
  SNRZQuadrature quad(pl);
  quad.PutLevelSymmetric(sn);
  SNRZSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutPL(pl);
  testa.SetQuadrature(&quad);
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.PutWriteFlux();
  testa.SetArray();
  real k1=testa.CalIgen();

  SNRZSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(voided_medium[i]){    
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	      med_voided[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      if(med_voided[repid].ExistNuclide(110230)){
	      real org=med_voided[repid].GetNuclide(110230).GetDensity();
	      med_voided[repid].GetNuclide(110230).PutDensity(org*void_ratio);
      };
      med_voided[repid].CalMacroFromMicro();
      if(i!=repid)med_voided[i].GetMacxs()=med_voided[repid].GetMacxs();
      test_v.AddMedium(med_voided[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutPL(pl);
  test_v.SetQuadrature(&quad);
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.PutWriteFlux();
  test_v.SetArray();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);
  real rho=testa.CalReactivity(&test_v,k1,kv);
#endif
  // -----------------------------------------------------
  high_resolution_clock::time_point t2 = high_resolution_clock::now();         
  duration<real, std::milli> time_span = t2 - t1;      
  cout<<"# Void reactivity computing time (in second):  "<<time_span.count()/1000<<"\n";
  // -----------------------------------------------------
  return rho;
};


real FRBurnerRZ::CalDopplerReactivity(int meds, int mede)
{
  // ! NOTE !

  if(!xscalc_for_reactivity){
    cout<<"# Error in FRBurnerRZ::CalDopplerReactivity.\n";
    cout<<"# Cross section in Doppler condition is NOT calculated.\n";
    exit(0);
  };

  CalMacroFromMicroBurnupMedium();

  if(mede==-1)mede=meds;
  if(meds==-1){
    meds=0;
    mede=med_burn-1;
  };

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  OnePointCalculator opc;

#if (calculation_theory != 1)
  // --- (DIFFUSION THEORY) ------------------------------
  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();
  real k1=testa.CalIgen();

  PLOSSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(i>=meds&&i<=mede){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	med_doppler[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      med_doppler[repid].CalMacroFromMicro();
      if(i!=repid)med_doppler[i].GetMacxs()=med_doppler[repid].GetMacxs();
      test_v.AddMedium(med_doppler[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);
  // -----------------------------------------------------
#endif

#if calculation_theory == 1
  // --- (TRANSPORT THEORY) --------------------------------
  SNRZQuadrature quad(pl);
  quad.PutLevelSymmetric(sn);

  SNRZSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutPL(pl);
  testa.SetQuadrature(&quad);
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.PutWriteFlux();
  testa.SetArray();
  real k1=testa.CalIgen();

  SNRZSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    if(i>=meds&&i<=mede){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];
      for(int j=0;j<burn_nuc;j++){
	med_doppler[repid].GetNuclideInTurn(j).PutDensity(density_ave[i][j]);
      };
      med_doppler[repid].CalMacroFromMicro();
      if(i!=repid)med_doppler[i].GetMacxs()=med_doppler[repid].GetMacxs();
      test_v.AddMedium(med_doppler[i]);
    }else{
      test_v.AddMedium(med[i]);
    };
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutPL(pl);
  test_v.SetQuadrature(&quad);
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.PutWriteFlux();
  test_v.SetArray();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);
#endif
  
  // -----------------------------------------------------

  real rho=testa.CalReactivity(&test_v,k1,kv);
  return rho;
};


real FRBurnerRZ::CalDelayedNeutronParameters(DelayedNeutronData &dnd)
{
  CalMacroFromMicroBurnupMedium();

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  // --- (DIFFUSION THEORY) ------------------------------
#if (calculation_theory != 1)
  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
  };
  testa.NoPrint();
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();
  real k1=testa.CalIgen();

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test.CalIgen(optname);
#endif  

  // --- (TRANSPORT THEORY) ------------------------------
#if calculation_theory == 1
  SNRZQuadrature quad(pl);
  quad.PutLevelSymmetric(sn);
  SNRZSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
  };
  testa.NoPrint();
  testa.PutPL(pl);
  testa.SetQuadrature(&quad);
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.PutWriteFlux();
  testa.SetArray();
  real k1=testa.CalIgen();

  SNRZSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
  };
  test.NoPrint();
  test.PutPL(pl);
  test.SetQuadrature(&quad);
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.PutWriteFlux();
  test.SetArray();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test.CalIgen();
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

  real beff=testa.CalBetaEffective(&test,dnd);
  testa.CalNeutronLifeTime(&test);
  return beff;
};

void FRBurnerRZ::WriteFileMediumData(string dirname,string filename)
{
  CalMacroFromMicroBurnupMedium();
  
  for(int i=0;i<mednum;i++){
    med[i].WriteFile(dirname,filename+IntToString(i),false); // micro XS cannot be outputed.
  };

  if(xscalc_for_reactivity){
    for(int i=0;i<med_burn;i++){
      int mz=mzone_per_med[i];
      int repid=mzone_rep[mz];

      // ... Doppler
      for(int j=0;j<burn_nuc;j++){
	med_doppler[repid].GetNuclideInTurn(j).PutDensity(med[i].GetNuclideInTurn(j).GetDensity());
      };
      med_doppler[repid].CalMacroFromMicro();
      med_doppler[repid].WriteFile(dirname,filename+IntToString(i)+"_dop",false); // micro XS cannot be outputed.      
      
      // ... Void
      for(int j=0;j<burn_nuc;j++){
	med_voided[repid].GetNuclideInTurn(j).PutDensity(med[i].GetNuclideInTurn(j).GetDensity());
      };
      if(med_voided[repid].ExistNuclide(110230)){
	med_voided[repid].GetNuclide(110230).PutDensity(sodium_density_in_voided_case);
      };
      med_voided[repid].CalMacroFromMicro();
      med_voided[repid].WriteFile(dirname,filename+IntToString(i)+"_vid",false); // micro XS cannot be outputed.            
      
    };
    
  };

};


void FRBurnerRZ::WriteFileNDData(string dirname,bool init,string addname)
// Total nuclide inventories are written on external files.
// If [init] is false, information on discharged fuel is printed out.
{
  //real avo=0.60221367;

  real factor=1.;
  if(cmi.GetBC(2)==1)factor=2.; // when axially-symmetric core is treated

  vector< vector<real> > den_sum((num_mzone+1),vector<real>(burn_nuc,0.)); // den_sum[material zone][nuclide]

  for(int ireg=0;ireg<med_burn;ireg++){
    int ii=mzone_per_med[ireg];
    int batch=GetBatchNumberFromMediumID(ireg);
    real vol=vol_per_med[ireg];
    int sz=density[ireg].size();
    int cyid=-1;
    for(int i=0;i<sz;i++){
      if(ir_cycle[ireg][i]==batch)cyid=i;
    };
    if(!init&&cyid==-1){
      cout<<"# Error in FRBurnerRZ::WriteFileNDData.\n";
      cout<<"# Assembly burnup does NOT reach the fuel exchange criteria.\n";
      cout<<"# In other words, the burnup does NOT reach equiliblium state.\n";
      exit(0);
    };
    for(int i=0;i<burn_nuc;i++){
      //real tmp=density_store[0][0][ireg][i];
      real tmp=fwd_nuc[0][0][ireg].get_dat(i);
      if(!init)tmp=density[ireg][cyid][i];
      den_sum[ii][i]+=tmp*vol/batch*factor;
    };
  };

  for(int i=0;i<burn_nuc;i++){
    for(int j=0;j<num_mzone;j++){
      den_sum[num_mzone][i]+=den_sum[j][i];
    };
  };

  //string filename_add[]={"ic","oc","rb","ab","all"};
  for(int ii=0;ii<=num_mzone;ii++){
    ofstream fout;
    string fname_add="all";
    if(ii!=num_mzone)fname_add=name_mzone[ii];
    string fname=dirname+"nddata."+addname+fname_add;
    fout.open(fname.data(),ios::out);
    if(fout.fail()){
      cout<<"# Failed to open the file.\n";
      exit(1);
    };
    fout<<"  "<<1<<"\n"; // Number of energy groups (dummy)
    fout<<"  "<<burn_nuc<<"\n";
    fout.setf(ios::scientific);
    fout.precision(9);
    for(int i=0;i<burn_nuc;i++){
      fout<<"  "<<med[0].GetNuclideInTurn(i).GetMatnum()<<"  \n"; // Material ID of the i-th nuclide
      fout<<"  "<<den_sum[ii][i]<<"\n";  // Number density of the i-th nuclide
      fout<<"  300\n";  // Temperature of the i-th nuclide (dummy)
      fout<<"  1\n"; // (dummy)
    // Number of energy groups of the i-th nuclide:
    // If this is zero, this nuclide has no cross section data.
    };
    fout.close();
  };

};

vector<real> FRBurnerRZ::GetHMNDData(int mzone_max)
{
  if(mzone_max==-1)mzone_max=num_mzone-1;
  
  bool init=false;
  
  int num_tru=21;
  int tru_list[]={
    922340,922350,922360,922370,922380,    932370,932390,942380,942390,942400,
    942410,942420,952410,952420,952421,    952430,962420,962430,962440,962450,
    962460
  };

  
  real factor=1.;
  if(cmi.GetBC(2)==1)factor=2.; // when axially-symmetric core is treated

  vector< vector<real> > den_sum((num_mzone+1),vector<real>(burn_nuc,0.)); // den_sum[material zone][nuclide]  

  for(int ireg=0;ireg<med_burn;ireg++){
    int ii=mzone_per_med[ireg];
    int batch=GetBatchNumberFromMediumID(ireg);
    real vol=vol_per_med[ireg];
    int sz=density[ireg].size();
    int cyid=-1;
    for(int i=0;i<sz;i++){
      if(ir_cycle[ireg][i]==batch)cyid=i;
    };
    if(!init&&cyid==-1){
      cout<<"# Error in FRBurnerRZ::WriteFileNDData.\n";
      cout<<"# Assembly burnup does NOT reach the fuel exchange criteria.\n";
      cout<<"# In other words, the burnup does NOT reach equiliblium state.\n";
      exit(0);
    };
    for(int i=0;i<burn_nuc;i++){
      //real tmp=density_store[0][0][ireg][i];
      real tmp=fwd_nuc[0][0][ireg].get_dat(i);
      if(!init)tmp=density[ireg][cyid][i];
      den_sum[ii][i]+=tmp*vol/batch*factor;
    };
  };

  
  for(int i=0;i<burn_nuc;i++){
    //for(int j=0;j<num_mzone;j++){
    for(int j=0;j<mzone_max+1;j++){      
      den_sum[num_mzone][i]+=den_sum[j][i];
    };
  };

  /*  
  for(int i=0;i<burn_nuc;i++){
    for(int j=0;j<num_mzone+1;j++){
      cout<<den_sum[j][i]<<" ";
    };
    cout<<"\n";
  };
  exit(0);
  */
  
  vector<real> nd_tru(num_tru,0.);
  for(int i=0;i<num_tru;i++){
    for(int j=0;j<burn_nuc;j++){
      if(med[0].GetNuclideInTurn(j).GetMatnum()==tru_list[i])nd_tru[i]=den_sum[num_mzone][j];
    };
  };

  return nd_tru;
};


void FRBurnerRZ::DischargedFuelAfterCooling(Burnup &bu, real cooling_year, int num, vector<string> &nuc_list_vec, vector<real> &atom_ratio, bool weight, bool xs)
{
  DischargedFuelAfterCooling(0, num_mzone-1, bu, cooling_year, num, nuc_list_vec, atom_ratio, weight, xs);
};


void FRBurnerRZ::DischargedFuelAfterCooling(int mzones, int mzonee, Burnup &bu, real cooling_year, int num, vector<string> &nuc_list_vec, vector<real> &atom_ratio, bool weight, bool xs)
{
  string *nuc_list=new string[num];
  for(int i=0;i<num;i++){
    nuc_list[i]=nuc_list_vec[i];
  };

  DischargedFuelAfterCooling(mzones, mzonee, bu, cooling_year, num, nuc_list, atom_ratio, weight, xs);

  delete [] nuc_list;
};

void FRBurnerRZ::DischargedFuelAfterCooling(Burnup &bu, real cooling_year, int num, string *nuc_list, vector<real> &atom_ratio, bool weight, bool xs)
{
  DischargedFuelAfterCooling(0, num_mzone-1, bu, cooling_year, num, nuc_list, atom_ratio, weight, xs);
};

void FRBurnerRZ::DischargedFuelAfterCooling(int mzones, int mzonee, Burnup &bu, real cooling_year, int num, string *nuc_list, vector<real> &atom_ratio, bool weight, bool xs)
{
  if(mzones<0||mzonee>=num_mzone||mzones>mzonee){
    cout<<"# Error in FRBurnerRZ::DischargedFuelAfterCooling.\n";
    cout<<"# Zone number assignment is inappropriate.\n";
    exit(0);
  };

  vector< vector<int> > each_nuc_turn;
  GetPrintingNuclideList(num, nuc_list, each_nuc_turn);

  real factor=1.;
  if(cmi.GetBC(2)==1)factor=2.; // when axially-symmetric core is treated

  vector< vector<real> > den_sum((num_mzone+1),vector<real>(burn_nuc,0.)); // den_sum[material zone][nuclide]

  for(int ireg=0;ireg<med_burn;ireg++){
    int ii=mzone_per_med[ireg];
    int batch=GetBatchNumberFromMediumID(ireg);
    real vol=vol_per_med[ireg];
    int sz=density[ireg].size();
    int cyid=-1;
    for(int i=0;i<sz;i++){
      if(ir_cycle[ireg][i]==batch)cyid=i;
    };
    if(cyid==-1){
      cout<<"# Error in FRBurnerRZ::DischargedFuelAfterCooling.\n";
      cout<<"# Assembly burnup does NOT reach the fuel exchange criteria.\n";
      cout<<"# In other words, the burnup does NOT reach equiliblium state.\n";
      exit(0);
    };
    for(int i=0;i<burn_nuc;i++){
      real tmp=density[ireg][cyid][i];
      den_sum[ii][i]+=tmp*vol/batch*factor;
    };
  };

  // Whole-core total inventory
  for(int i=0;i<burn_nuc;i++){
    //for(int j=0;j<num_mzone;j++){
    for(int j=mzones;j<=mzonee;j++){
      den_sum[num_mzone][i]+=den_sum[j][i];
    };
  };

  bu.PutDensity(den_sum[num_mzone]);
  bu.BurnupCalculation(0, cooling_year*365*24*60*60,"mmpa",false);
  den_sum[num_mzone]=bu.GetDensity();

  // Averaged 1-group (n,g) cross section calculations
  vector<real> xs1g_ave(num);
  for(int i=0;i<num;i++){
    real sum1=0.;
    real sum2=0.;
    for(int j=0;j<each_nuc_turn[i].size();j++){
      real tmp=den_sum[num_mzone][each_nuc_turn[i][j]];
      sum1+=tmp;
      sum2+=tmp*sigc_1g[0][trace_cycle-1][cycle_div-1][each_nuc_turn[i][j]];
    };
    xs1g_ave[i]=sum2/sum1;
  };

  // Conversion from ND to weight
  if(weight){
    for(int i=0;i<burn_nuc;i++){
      int mid=med[0].GetNuclideInTurn(i).GetMatnum();
      real aw=bu.GetAtomicWeight(mid);
      den_sum[num_mzone][i]*=aw*0.001/avo; // [kg]
    };
  };

  vector<real> val(num);
  for(int i=0;i<num;i++){
    val[i]=0.;
    for(int j=0;j<each_nuc_turn[i].size();j++){
      val[i]+=den_sum[num_mzone][each_nuc_turn[i][j]];
    };
  };

  real total=0.;
  for(int i=0;i<num;i++){
    total+=val[i];
  };

  atom_ratio.resize(num);
  for(int i=0;i<num;i++){
    cout<<"# "<<nuc_list[i]<<" ";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<val[i]<<" ";  
    cout.setf(ios::showpoint);
    cout.precision(5);
    cout<<val[i]/total*100.<<"\n";
    atom_ratio[i]=val[i]/total;
  };
  cout<<"# Total : "<<total<<"\n";

  // overwriting averaged 1-group (n,g) cross section
  // at medium 0 at EOEC
  if(xs){
  for(int i=0;i<num;i++){
    atom_ratio[i]=xs1g_ave[i];
  };
  };

};

void FRBurnerRZ::GetPrintingNuclideList(int prt_nuc, string *prt_nuc_nam, vector< vector<int> > &each_nuc_turn)
{
  // Taken from Burner::ShowNumberDensityHistory
  int nucn=burn_nuc;

  each_nuc_turn.resize(prt_nuc);

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
	if((matid<800000&&matid>=310000)||matid>9900000){
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
};

// +++ printing

void FRBurnerRZ::PrintNuclideWeightPerBatch(Burnup &bu,bool init)
{
  //real avo=0.60221367;

  real factor=1.;
  if(cmi.GetBC(2)==1)factor=2.; // when axially-symmetric core is treated

  vector< vector<real> > den_sum((num_mzone),vector<real>(burn_nuc,0.));
  //int ii=0;
  for(int ireg=0;ireg<med_burn;ireg++){
    int ii=mzone_per_med[ireg];
    int batch=GetBatchNumberFromMediumID(ireg);
    real vol=vol_per_med[ireg];
    int sz=density[ireg].size();
    int cyid=-1;
    for(int i=0;i<sz;i++){
      if(ir_cycle[ireg][i]==batch)cyid=i;
    };
    if(!init&&cyid==-1){
      cout<<"# Error in FRBurnerRZ::PrintNuclideWeightPerBatch.\n";
      cout<<"# Assembly burnup does NOT reach the fuel exchange criteria.\n";
      cout<<"# In other words, the burnup does NOT reach equiliblium state.\n";
      exit(0);
    };
    for(int i=0;i<burn_nuc;i++){
      //real tmp=density_store[0][0][ireg][i];
      real tmp=fwd_nuc[0][0][ireg].get_dat(i);
      if(!init)tmp=density[ireg][cyid][i];
      den_sum[ii][i]+=tmp*vol/batch*factor;
    };
    //if(ireg==ic_e||ireg==oc_e||ireg==rb_e)ii++;
  };

  cout<<"##########################################################\n";
  cout<<"# Nuclide Weight [kg] of\n";
  if(init){
    cout<<"# loaded fuel assembly per batch\n";
  }else{
    cout<<"# discharged fuel assembly per batch\n";
  };
  if(cmi.GetBC(2)==1){
    cout<<"#\n";
    cout<<"# Note that you are treating axially-symmetric core,\n";
    cout<<"# so the following quantities have been multiplied by 2.\n";
    cout<<"#\n";
  };
  for(int i=0;i<num_mzone;i++){
    cout<<"# "<<i<<" : "<<name_mzone[i]<<"\n";
  };
  cout<<"#\n";
  cout<<"##########################################################\n";
  //cout<<"#(Nucl)     (IC)     (OC)     (RB)     (AB)    (Sum)\n";
  cout<<"#(Nucl)    ";
  for(int i=0;i<num_mzone;i++){
    cout<<i<<"         ";
  };
  cout<<"(Sum)\n";

  vector<real> wgtsum(num_mzone,0.);
  vector<real> wgtsum_hm(num_mzone,0.);
  for(int i=0;i<burn_nuc;i++){
    int mid=med[0].GetNuclideInTurn(i).GetMatnum();
    real aw=bu.GetAtomicWeight(mid);
    WriteOut(mid,7);
    cout<<" ";
    real tmp=0.;
    for(int j=0;j<num_mzone;j++){
      real wgt=den_sum[j][i]/avo*aw*1e-3; // [kg]
      WriteOut(wgt,"%9.3f");
      cout<<" ";
      wgtsum[j]+=wgt;
      if(mid>=900000)wgtsum_hm[j]+=wgt;
      tmp+=wgt;
    };
    WriteOut(tmp,"%9.3f");
    cout<<"\n";
  };

  cout<<" (Sum)  ";
  real tmp=0.;
  for(int j=0;j<num_mzone;j++){
    WriteOut(wgtsum[j],"%9.3f");
    tmp+=wgtsum[j];
    cout<<" ";
  };
  WriteOut(tmp,"%9.3f");
  cout<<"\n";

  cout<<" (HM)   ";
  tmp=0.;
  for(int j=0;j<num_mzone;j++){
    WriteOut(wgtsum_hm[j],"%9.3f");
    tmp+=wgtsum_hm[j];
    cout<<" ";
  };
  WriteOut(tmp,"%9.3f");
  cout<<"\n";
  cout<<"##########################################################\n";
};

void FRBurnerRZ::PrintND(Burnup &bu,string opt)
{
  //real avo=0.60221367;
  //real ev_to_j=1.60219e-19;
  int batch_max=0;

  vector< vector< vector<real> > > dsum(num_mzone);
  for(int i=0;i<num_mzone;i++){
    dsum[i].resize(batch_mzone[i]);
    if(batch_mzone[i]>batch_max)batch_max=batch_mzone[i];
    for(int j=0;j<batch_mzone[i];j++){
      dsum[i][j].resize(burn_nuc,0.);
    };
  };

  vector<real> vol_mz(num_mzone,0.);
  for(int i=0;i<med_burn;i++){
    int id=mzone_per_med[i];
    if(id!=-1){
      vol_mz[id]+=vol_per_med[i];
    };
  };

  //vector< vector<real> > density_sum(batch_max);
  vector< vector<real> > density_sum_all(batch_max);  

  for(int i=0;i<batch_max;i++){
    //density_sum[i].resize(burn_nuc,0.);
    density_sum_all[i].resize(burn_nuc,0.);    
  };
  real vol_sum=0.;
  real vol_sum_all=0.;

  cout<<"\n# +++ Batch-wise data +++\n";
  if(opt=="nd"){
    cout<<"#  Number density [#/cm3]\n";
  }else if(opt=="kg_per_thm"){
    cout<<"#  Weight [kg/tHM]\n";
  }else if(opt=="w_per_thm"){
    cout<<"#  Decay heat [W/tHM]\n";
  }else{
    cout<<"# Error stop !\n";
    cout<<"# No keyword : "<<opt<<"\n";
    exit(0);
  };

  cout.setf(ios::scientific);
  cout.precision(4);

  for(int ireg=0;ireg<med_burn;ireg++){
    int mz=mzone_per_med[ireg];
     real vol=vol_per_med[ireg];
     vol_sum+=vol;
     vol_sum_all+=vol;
     cout<<"#\n# Burnup region : "<<ireg;
     cout<<"  (Volume : "<<vol<<" cm3)\n#\n";
     int sz=density[ireg].size();

     cout<<"#(Cycles) ";
     for(int i=0;i<sz;i++){
       cout<<ir_cycle[ireg][i]<<"          ";
     };
     if(sz>1)cout<<"   (Average)";
     cout<<"\n";

     cout<<"#(GWd/t)  ";
     {
     real avg=0.;
     for(int i=0;i<sz;i++){
       cout<<ac_burn[ireg][i]<<" ";
       avg+=ac_burn[ireg][i];
     };
     avg/=sz;
     if(sz>1)cout<<"   "<<avg;
     cout<<"\n";
     };

     cout<<"#---------";
     for(int i=0;i<sz;i++){
       cout<<"-----------";
     };
     cout<<"\n";

     vector<real> sum(sz,0.);

     for(int i=0;i<burn_nuc;i++){
       int mid=med[0].GetNuclideInTurn(i).GetMatnum();
       real aw=bu.GetAtomicWeight(mid);
       real dc=bu.GetBurnupChain().GetDecayConstant(mid);
       real e=0.;
       for(int k=0;k<3;k++){
	 e+=bu.GetBurnupChain().GetDecayEnergy(mid,k);
       };
       WriteOut(midt.Name(mid),8);
       cout<<"  ";
       //cout<<midt.Name(mid)<<"  ";
       real avg=0.;
       for(int j=0;j<sz;j++){
	 real val=0.;
	 if(opt=="nd"){
	   val=density[ireg][j][i];
	 }else if(opt=="kg_per_thm"){
	   val=(density[ireg][j][i]/avo*aw)*1e-3;
	   val/=((thm_per_med[ireg]/vol_per_med[ireg])*1e-6);
	 }else if(opt=="w_per_thm"){
           val=density[ireg][j][i]*1e24*dc*e*ev_to_j;
	   val/=((thm_per_med[ireg]/vol_per_med[ireg])*1e-6);
	 };
	 cout<<val<<" ";
	 sum[j]+=val;
         avg+=val;      
         density_sum_all[j][i]+=val*vol;
	 dsum[mz][j][i]+=val*vol; 
       };
       avg/=sz;
       if(sz>1)cout<<"   "<<avg;
       cout<<"\n";
     };
     cout<<"#(sum)    ";
     for(int i=0;i<sz;i++){
       cout<<sum[i]<<" ";
     };
     cout<<"\n#\n";

  };

  cout<<"#\n# Material zone-averaged\n#\n";
  for(int i=0;i<num_mzone;i++){
    cout<<"#\n#  Material zone "<<i<<" : "<<name_mzone[i];
    cout<<" ( volume : "<<vol_mz[i]<<" cm3 )\n#\n";
    int sz=dsum[i].size();
    for(int j=0;j<burn_nuc;j++){
      int mid=med[0].GetNuclideInTurn(j).GetMatnum();
      WriteOut(midt.Name(mid),8);
      cout<<"  ";
      real sum=0.;
      for(int k=0;k<sz;k++){
	real tmp=dsum[i][k][j]/vol_mz[i];
        cout<<tmp<<" ";
        sum+=tmp;
      };
      sum/=sz;
      if(sz>1)cout<<"   "<<sum;
      cout<<"\n";
    };
  };

  bool sz1=false;
  for(int ireg=0;ireg<med_burn;ireg++){
     int sz=density[ireg].size();
     if(sz!=1)sz1=true;
  };

  if(!sz1){
  cout<<"#\n# +++ Whole core-averaged";
  cout<<" ( volume : "<<vol_sum_all<<" cm3 )\n#\n";
  for(int i=0;i<burn_nuc;i++){
    int mid=med[0].GetNuclideInTurn(i).GetMatnum();
    WriteOut(midt.Name(mid),8);
    cout<<"  ";
    cout<<density_sum_all[0][i]/vol_sum_all<<" ";
    cout<<"\n";
  };
  };
};

void FRBurnerRZ::PrintND(int prt_nuc,string *prt_nuc_nam,int medi,int mede)
{
  cout<<"\n";
  cout<<"# Medium-wise averaged number density after burnup.\n";
  cout<<"#\n";

  if(mede==-1)mede=med_burn-1;

  vector<int> prt_nuc_turn(prt_nuc);
  vector<int> prt_nuc_id(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    int id=midt.ID(prt_nuc_nam[i]);
    prt_nuc_id[i]=id;
    for(int j=0;j<burn_nuc;j++){
      if(med[0].GetNuclideInTurn(j).GetMatnum()==id){
	prt_nuc_turn[i]=j;
	break;
      };
    };
  };

  cout<<"#Med GWd/t ";
  for(int i=0;i<prt_nuc;i++){
    //cout<<prt_nuc_nam[i]<<" ";
    WriteOut(prt_nuc_nam[i],8);
    cout<<"    ";
  };
  cout<<"\n";

  for(int i=medi;i<=mede;i++){
    int sz=density[i].size();
    int ir=0;
    int irmax=0;
    for(int j=0;j<sz;j++){
      if(ir_cycle[i][j]>irmax){
	irmax=ir_cycle[i][j];
	ir=j;
      };
    };
    WriteOut(i,3);
    //cout.setf(ios::scientific);
    //cout.precision(3);
    //cout<<ac_burn[i][ir]<<" ";
    cout<<"  ";
    WriteOut(ac_burn[i][ir],"%5.1f");
    cout<<" ";
    //cout.precision(5);
    cout.setf(ios::scientific);
    cout.precision(5);
    for(int j=0;j<prt_nuc;j++){
      cout<<density[i][ir][prt_nuc_turn[j]]<<" ";
    };
    cout<<"\n";
  }; 
};

void FRBurnerRZ::PrintBurnup(int medi,int mede)
{
  cout<<"#\n# Medium-dependent burnup [GWD/t]\n#\n";
  cout<<"#Med   AVG   MAX\n#\n";
  if(mede==-1)mede=med_burn-1;
  for(int i=medi;i<=mede;i++){
    int sz=density[i].size();
    real max=0;
    real sum=0.;
    for(int j=0;j<sz;j++){
      sum+=ac_burn[i][j];
      if(ac_burn[i][j]>max){
	max=ac_burn[i][j];
      };
    };
    sum/=sz;
    WriteOut(i,3);
    cout<<"  ";
    WriteOut(sum,"%5.1f");
    cout<<"  ";
    WriteOut(max,"%5.1f");
    cout<<"\n";
  };

  // Summary printing
  if(medi==0&&mede==med_burn-1){

    cout<<"#\n# Material zone-dependent burnup [GWD/t]\n#\n";
    cout<<"#MZone   AVG   MAX    (THM)\n#\n";
    vector<real> avgbu(num_mzone,0.);
    vector<real> maxbu(num_mzone,0.);
    vector<real> sumwgt(num_mzone,0.);
    for(int i=0;i<med_burn;i++){
      int sz=density[i].size();
      real thm=thm_per_med[i]/sz;
      real max=0;
      real sum=0.;
      for(int j=0;j<sz;j++){
        sum+=ac_burn[i][j];
        if(ac_burn[i][j]>max)max=ac_burn[i][j];
      };
      sum/=sz;
      int mm=mzone_per_med[i];
      avgbu[mm]+=sum*thm;
      maxbu[mm]+=max*thm;
      sumwgt[mm]+=thm;
    };

    real sum1=0.;
    real sum2=0.;
    real sum3=0.;
    for(int i=0;i<num_mzone;i++){
      cout<<name_mzone[i]<<"  ";
      WriteOut(avgbu[i]/sumwgt[i],"%5.1f");
      cout<<"  ";
      WriteOut(maxbu[i]/sumwgt[i],"%5.1f");
      cout<<"  ";
      WriteOut(sumwgt[i]*1e-6,"%5.1f");
      cout<<"\n";
      sum1+=avgbu[i];
      sum2+=maxbu[i];
      sum3+=sumwgt[i];
    };
    cout<<" All     ";
    WriteOut(sum1/sum3,"%5.1f");
    cout<<"  ";
    WriteOut(sum2/sum3,"%5.1f");
    cout<<"\n";
  };

};

void FRBurnerRZ::PrintMacroXS()
{
  cout<<"###################################\n";
  cout<<"# Macroscopic cross section table #\n";
  cout<<"###################################\n";

  for(int i=0;i<num_mzone;i++){
    cout<<"\n\n# Material zone "<<i<<" : "<<name_mzone[i]<<"\n#\n";
    med[mzone_rep[i]].PrintMacroSectionTable();
  };
};

void FRBurnerRZ::PrintMicroXS(int nucid)
{
  cout<<"#####################################################\n";
  cout<<"# Microscopic cross section table of nuclide "<<nucid<<" #\n";
  cout<<"#####################################################\n";

  int num=1;
  if(xscalc_for_reactivity)num=3;

  for(int ii=0;ii<num;ii++){

    if(ii==0&&num==3){
      cout<<"#\n";
      cout<<"# (Reference state)\n";
      cout<<"#\n";
    }else if(ii==1){
      cout<<"#\n";
      cout<<"# (Voided state)\n";
      cout<<"#\n";
    }else {
      cout<<"#\n";
      cout<<"# (Doppler state)\n";
      cout<<"#\n";
    };

    for(int i=0;i<num_mzone;i++){
      int tmp=med[mzone_rep[i]].SearchNuclide(nucid);
      if(tmp!=0){
        cout<<"\n\n# Material zone "<<i<<" : "<<name_mzone[i]<<"\n#\n";
        GroupData1D ebnd=med[mzone_rep[i]].GetEnband();
	if(ii==0){
          med[mzone_rep[i]].GetNuclideInTurn(tmp).PrintMicroSectionTable(ebnd);
	}else if(ii==1){
          med_voided[mzone_rep[i]].GetNuclideInTurn(tmp).PrintMicroSectionTable(ebnd);
	}else{
          med_doppler[mzone_rep[i]].GetNuclideInTurn(tmp).PrintMicroSectionTable(ebnd);
	};
      };
    };

  };

};

void FRBurnerRZ::PrintMacroXSMediumWise(int inp)
{
  if(inp<0||inp>=mednum){
    cout<<"# Error in FRBurnerRZ::PrintMacroXSMediumWise.\n";
    cout<<"# Medium ID "<<inp<<" does not exist.\n";
    exit(0);
  };

  med[inp].CalMacroFromMicro();
  cout<<"###################################\n";
  cout<<"# Macroscopic cross section table #\n";
  cout<<"# Medium ID : "<<inp<<"\n";
  cout<<"###################################\n";

  med[inp].PrintMacroSectionTable();
};

void FRBurnerRZ::PrintInitialND()
{
  cout<<"# Initial number density [/barn/cm]\n#\n";
  for(int i=0;i<num_mzone;i++){
    cout<<"# "<<i<<" : "<<name_mzone[i]<<"\n";
  };
  cout<<"#\n";

  cout<<"#         ";
  for(int i=0;i<num_mzone;i++){
    cout<<i<<"          ";
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  int tmp=med[0].GetNucnum();
  for(int i=0;i<tmp;i++){
    //for(int i=0;i<burn_nuc;i++){
    int matid=med[0].GetNuclideInTurn(i).GetMatnum();
    string name=midt.Name(matid);
    vector<real> den(num_mzone);
    real sum=0.;
    for(int j=0;j<num_mzone;j++){
      den[j]=med[mzone_rep[j]].GetNuclideInTurn(i).GetDensity();
      sum+=den[j];
    };
    if(sum>0.){
      WriteOut(name,8);
      cout<<"  ";
      for(int j=0;j<num_mzone;j++){
	cout<<den[j]<<" ";
      };
      cout<<"\n";
    };
  };
};

void FRBurnerRZ::PerturbInitialND(int mzone_inp, int nuclide_id, real rel_change)
{
  if(mzone_inp<0||mzone_inp>=num_mzone){
    cout<<"# Error in FRBurnerRZ::PerturbInitialND.\n";
    cout<<"# Material zone setting (ID="<<mzone_inp<<") is NOT appropriate.\n";
    exit(0);
  };

  int rep=mzone_rep[mzone_inp];
  int tmp=med[rep].GetNucnum();
  for(int i=0;i<tmp;i++){
    int matid=med[rep].GetNuclideInTurn(i).GetMatnum();
    if(matid==nuclide_id){
      real org=med[rep].GetNuclideInTurn(i).GetDensity();
      med[rep].GetNuclideInTurn(i).PutDensity(org+org*rel_change);            
    };
  };
};

void FRBurnerRZ::PrintFluxLevelHistory()
{

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<mednum;i++){
    cout<<"#\n";
    cout<<"# Time-dependent flux level [/s/cm2]\n";
    cout<<"# Medium ID : "<<i<<" (Vol.: "<<vol_per_med[i]<<"[cm3])\n";
    cout<<"# [ST] [SST]  [FLUX]\n"; 
    for(int j=0;j<trace_cycle;j++){
      for(int k=0;k<cycle_div+1;k++){
        WriteOut(j,5);
        WriteOut(k,5);
        cout<<"   "<<flux_level[j][k][i]<<"\n";
      };
    };
    cout<<"\n";
  };
};

void FRBurnerRZ::Print1groupXS(string nucname,int medid)
{
  if(medid>=mednum){
    cout<<"# Error in FRBurnerRZ::Print1groupXS.\n";
    cout<<"# Medium "<<medid<<" does not exist.\n";
    exit(0);
  };

  int nucid=midt.ID(nucname);
  int nuct=-1;
  for(int i=0;i<burn_nuc;i++){
    int idtmp=med[medid].GetNuclideInTurn(i).GetMatnum();
    if(idtmp==nucid){
      nuct=i;
      break;
    };
  };
  if(nuct==-1){
    cout<<"# Error in FRBurnerRZ::Print1groupXS.\n";
    cout<<"# Nuclide "<<nucname<<" is not found.\n";
    exit(0);
  };

  if(medid>=med_burn){
    cout<<"# Error in FRBurnerRZ::Print1groupXS.\n";
    cout<<"# Medium "<<medid<<" is NOT fuel medium.\n";
    cout<<"# This method cannot be used for non-fuel medium.\n";
    exit(0);
  };

  cout<<"#\n# 1-group cross section of "<<nucname<<"\n";
  cout<<"#           (medium ID is "<<medid<<")\n";
  cout<<"#\n";
  cout<<"#Cycle  Step  Capture    Fission    (n,2n)\n";
  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<trace_cycle;i++){
    for(int j=0;j<cycle_div;j++){
      WriteOut(i,5);
      cout<<" ";
      WriteOut(j,5);
      cout<<"  ";
      cout<<sigc_1g[medid][i][j][nuct]<<"  ";
      cout<<sigf_1g[medid][i][j][nuct]<<"  ";
      cout<<sign2n_1g[medid][i][j][nuct]<<"  ";
      cout<<"\n";
    };
  };
  cout<<"#\n";
};

void FRBurnerRZ::AddActinideDecayHeatDataToBurnupChain(string cbglibdir,string fname,Burnup &bu)
{
  cbglibdir=cbglibdir+"CBGLIB_BURN/DDfile/"+fname;

  BCGManager man;
  man.ReadDecayDataFromFile(cbglibdir);

  /*
  int iill=14;
  int ialist[]={ 92, 92, 92, 93, 93, 94, 94, 94, 94, 95, 95, 95,  96, 96};
  int izlist[]={232,234,237,237,239,238,239,240,241,241,242,242, 242,244};
  int illist[]={  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,   0,  0};  
  vector<int> matno(iill);
  for(int i=0;i<iill;i++){
    matno[i]=ialist[i]*10000+izlist[i]*10+illist[i];
  };
  for(int i=0;i<iill;i++){
    for(int j=0;j<3;j++){
      real ee=man.GetNuclide(ialist[i],izlist[i],0).GetDecayEnergy(j);
      bu.GetBurnupChain().GetNuclideChainData(matno[i]).PutDecayEnergy(j,ee);
    };
  };
  */
  for(int i=90;i<100;i++){
    for(int j=200;j<300;j++){
      for(int k=0;k<2;k++){
	if(man.GetNuclideIndex(i,j,k)!=-1){
  	  int matno=i*10000+j*10+k;
	  if(bu.GetBurnupChain().GetNuclideChainData(matno).GetID()==matno){
            for(int jj=0;jj<3;jj++){
              real ee=man.GetNuclide(i,j,k).GetDecayEnergy(jj);
              bu.GetBurnupChain().GetNuclideChainData(matno).PutDecayEnergy(jj,ee);
            };
	  };
	};
      };
    };
  };
};

real FRBurnerRZ::GetDensity(int mid, int batch, int nucid)
{
  for(int i=0;i<burn_nuc;i++){
    int id=med[mid].GetNuclideInTurn(i).GetMatnum();
    if(id==nucid)return density[mid][batch][i];
  };
  cout<<"# Error in FRBurnerRZ::GetDensity.\n";
  exit(0);
};

// +++ Ogawa calculation
void FRBurnerRZ::OgawaCalculationOn(real y1, real y2, real y3, int mat1, int mat2)
{
  ogawa_cal=true;
  ogawa_year.resize(3);
  ogawa_year[0]=y1;
  ogawa_year[1]=y2;
  ogawa_year[2]=y3;
  ogawa_mat1=mat1*10000;
  ogawa_mat2=(mat2+1)*10000;
};

// +++ Material zone assignemtn

void FRBurnerRZ::PutNumberMZone(int i)
{
  num_mzone=i;

  num_med_mzone.resize(num_mzone);
  iden_mzone.resize(num_mzone);
  batch_mzone.resize(num_mzone);
  name_mzone.resize(num_mzone);
  num_pin_mzone.resize(num_mzone);
  area_assembly_mzone.resize(num_mzone);
  mzone_rep.resize(num_mzone,-1);
};

void FRBurnerRZ::PutMZoneInfo(int *mzone)
{
  mzone_per_med.resize(mednum);
  for(int i=0;i<mednum;i++){
    mzone_per_med[i]=-1;
  };

  int ir=cmi.GetXC();
  int iz=cmi.GetYC();
  int index=0;
  for(int z=0;z<iz;z++){
    for(int r=0;r<ir;r++){
      int mat=cmi.GetCMat(index); 
      int zone=mzone[index];
      if(zone<num_mzone){
	if(mzone_per_med[mat]==-1){
  	  num_med_mzone[zone]++;
	  mzone_per_med[mat]=zone;
	};
      };
      index++;
    };
  };

  bool negative=false;
  for(int i=0;i<mednum;i++){
    if(mzone_per_med[i]!=-1&&negative){
      cout<<"# Error in FRBurnerRZ::PutMZoneInfo.\n";
      cout<<"# Order of material (Medium) ID is inappropriate.\n";
      exit(0);
    };
    if(mzone_per_med[i]==-1)negative=true;
  };

  /*
  for(int i=0;i<num_mzone;i++){cout<<i<<" "<<num_med_mzone[i]<<"\n";};
  for(int i=0;i<mednum;i++){cout<<i<<" "<<mzone_per_med[i]<<"\n";};
  */
};

void FRBurnerRZ::PutNameMZone(string *nameinp)
{
  for(int i=0;i<num_mzone;i++){
    name_mzone[i]=nameinp[i];
  };
};

void FRBurnerRZ::PutBatchMZone(int *batchinp)
{
  for(int i=0;i<num_mzone;i++){
    batch_mzone[i]=batchinp[i];
  };
};

void FRBurnerRZ::PutMediumToMZone(Medium &min, int mzone)
{
  if(mzone<0||mzone>=num_mzone){
    cout<<"# Error in FRBurnerRZ::PutMediumToZone.\n";
    cout<<"# Material zone "<<mzone<<" is NOT defined.\n";
    exit(0);
  };

  for(int i=0;i<mednum;i++){
    //if(mzone_per_med[i]==mzone)med[i]=min;
    if(mzone_per_med[i]==mzone&&mzone_rep[mzone]==-1){
      //med[i]=min;
      mzone_rep[mzone]=i;
      PutSADataToMZoneDefiningMediumData(min, mzone, 0);
    };
  };

};

void FRBurnerRZ::PutSADataToMZone(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib, int model)
{
  for(int i=0;i<mednum;i++){
    if(mzone_per_med[i]==mzone&&mzone_rep[mzone]==-1){
      mzone_rep[mzone]=i;
    };
  };

  num_pin_mzone[mzone]=sa.GetFuelSAGeometry()->GetNumberOfPin();
  area_assembly_mzone[mzone]=sa.GetFuelSAGeometry()->GetAssemblyVolume()*0.01; // convert from mm^2 to cm^2

  if(model==0){
    PutSADataToMZoneHomogeneous(sa, mzone, xslib);
  }else if(model==1){
    PutSADataToMZoneHeterogeneous1D(sa, mzone, xslib);
  }else if(model==2){
    PutSADataToMZoneHeterogeneous2D(sa, mzone, xslib);
  }else{
    cout<<"# Error in FRBurnerRZ::PutSADataToMZone.\n";
    cout<<"# Assembly model "<<model<<" cannot be assigned.\n";
    exit(0);
  };

};

void FRBurnerRZ::PutSADataToNonfuel(FRDTSubAssembly &sa, int medid)
{
  int fuel_nuc_num=matno.size();
  med[medid].PutImax(group);
  int plinp=pl;
  if(plinp==0)plinp=1;
  med[medid].PutPL(plinp);
  med[medid].PutNuclide(fuel_nuc_num,matno);
  sa.PutHomogenizedNumberDensity(med[medid]);
  med[medid].PutTemperatureForAllNuclide(sa.GetTemperature());
  med[medid].GetNuclide(110230).PutTemperature(sa.GetSodiumTemperature());
  med_data_inp[medid]=true;
};

void FRBurnerRZ::PutSADataToMZoneHomogeneous(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib)
{
  int fuel_nuc_num=matno.size();

  Medium med_ref;
  
  int num_case=3;
  if(!xscalc_for_reactivity)num_case=1;

  for(int ii=0;ii<num_case;ii++){ // ii: 0/ordinary state, 1/voided state, 2/Doppler state

    med_ref.PutImax(group);
    int plinp=pl;
    if(plinp==0)plinp=1;
    med_ref.PutPL(plinp);
    med_ref.PutNuclide(fuel_nuc_num,matno);
    sa.PutHomogenizedNumberDensity(med_ref);
    med_ref.PutTemperatureForAllNuclide(sa.GetTemperature());
    med_ref.GetNuclide(110230).PutTemperature(sa.GetSodiumTemperature());

    if(ii==1){ // sodium-voided system
      med_ref.GetNuclide(110230).PutDensity(sodium_density_in_voided_case);
    };
    if(ii==2){ // doppler system
      int nucnum=med_ref.GetNucnum();
      for(int j=0;j<nucnum;j++){
        int nucid=med_ref.GetNuclideInTurn(j).GetMatnum();
        if(nucid>900000||nucid==420000||nucid==400000){ // default
	//if(nucid>200000){
	//if(nucid>100000){
          real org=med_ref.GetNuclideInTurn(j).GetTemperature();
          med_ref.GetNuclideInTurn(j).PutTemperature(org+delta_t_in_doppler_case);
	};
      };
    };

    opc.CalSelfShieldingInfiniteSystem(med_ref,xslib);

    PutSADataToMZoneDefiningMediumData(med_ref, mzone, ii);

  };

};

void FRBurnerRZ::PutSADataToMZoneHeterogeneous1D(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib)
{
  int fuel_nuc_num=matno.size();

  Medium med_ref;
  
  int num_case=3;
  if(!xscalc_for_reactivity)num_case=1;

  real area_spacer_wire=sa.GetFuelSAGeometry()->GetSpacerwireVolume();
  if(area_spacer_wire>1e-10){
    cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneous1D.\n";
    cout<<"# Not yet coded for cases of non-zero spacer-wire volume.\n";
    exit(0);
  };

  int num_pin=sa.GetFuelSAGeometry()->GetNumberOfPin();
  real pin_pitch=sa.GetFuelSAGeometry()->GetPinPitch();

  // ... single pincell geometric information
  real area_clad=sa.GetFuelSAGeometry()->GetCladVolume()/num_pin;
  real area_pellet=sa.GetFuelSAGeometry()->GetPelletVolume()/num_pin;
  real area_pellet_invoid=sa.GetFuelSAGeometry()->GetPelletVoidVolume()/num_pin;
  real area_pellet_gap=sa.GetFuelSAGeometry()->GetPelletGapVolume()/num_pin;
  real area_total=0.5*sqrt(3.)*pin_pitch*pin_pitch;
  real area_sodium=area_total-(area_clad+area_pellet+area_pellet_invoid+area_pellet_gap);
  //cout<<area_clad<<" "<<area_pellet<<" "<<area_pellet_invoid<<" "<<area_pellet_gap<<" "<<area_sodium<<"\n"; exit(0);

  real area_vpg=area_pellet_invoid+area_pellet+area_pellet_gap;

  real area_sum=0.;
  vector<real> outer_radius_ring;
  vector<int> material_ring; // 0:v+p+g, 1:cladding, 2:sodium, 3:duct

  // Center pin (first layer)
  area_sum+=area_vpg;
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(0);

  area_sum+=area_clad;
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(1);

  area_sum+=area_sodium;
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(2);

  // Non-central pins (second layer, third layer, ...)
  int sum_pin=1;
  for(int i=1;i<20;i++){
    int numpin=i*6;
    sum_pin+=numpin;

    area_sum+=area_sodium*numpin*0.5;
    outer_radius_ring.push_back(sqrt(area_sum/PI));
    material_ring.push_back(2);

    area_sum+=area_clad*numpin*0.5;
    outer_radius_ring.push_back(sqrt(area_sum/PI));
    material_ring.push_back(1);

    area_sum+=area_vpg*numpin;
    outer_radius_ring.push_back(sqrt(area_sum/PI));
    material_ring.push_back(0);

    area_sum+=area_clad*numpin*0.5;
    outer_radius_ring.push_back(sqrt(area_sum/PI));
    material_ring.push_back(1);

    area_sum+=area_sodium*numpin*0.5;
    outer_radius_ring.push_back(sqrt(area_sum/PI));
    material_ring.push_back(2);

    //cout<<i<<" "<<numpin<<" "<<sum_pin<<"\n";
    if(sum_pin==sa.GetFuelSAGeometry()->GetNumberOfPin())break;

    if(i==19){
      cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneousModel.\n";
      cout<<"# Ther number of pin layers cannot be determined from the total number of pins.\n";
      exit(0);
    };
  };

  real area_total_sodium=sa.GetFuelSAGeometry()->GetCoolantVolume();
  real area_outduct_sodium=sa.GetFuelSAGeometry()->GetOutDuctVolume();
  real area_remaining_induct_sodium=area_total_sodium-area_outduct_sodium-area_sodium*num_pin;

  area_sum+=area_remaining_induct_sodium;
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(2);
  
  area_sum+=sa.GetFuelSAGeometry()->GetDuctVolume();
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(3);
  
  area_sum+=area_outduct_sodium;
  outer_radius_ring.push_back(sqrt(area_sum/PI));
  material_ring.push_back(2);

  int sz=outer_radius_ring.size();

  // ... Unit is converted from [mm] to [cm]
  for(int i=0;i<sz;i++){
    outer_radius_ring[i]*=0.1;
  };

#if 0
  cout<<"# Ring information.\n";
  for(int i=0;i<sz;i++){
    cout<<i<<" "<<material_ring[i]<<" "<<outer_radius_ring[i]<<"\n";
  };
  exit(0);
#endif


  // Adjacent-sodium regions are unified
  for(int i=0;i<sz-1;i++){
    if(material_ring[i]==material_ring[i+1]){
      outer_radius_ring[i]=outer_radius_ring[i+1];
      sz--;
      for(int j=i+1;j<sz;j++){
	material_ring[j]=material_ring[j+1];
        outer_radius_ring[j]=outer_radius_ring[j+1];
      };
    };
  };

#if 0
  cout<<"# Ring information.\n";
  for(int i=0;i<sz;i++){
    cout<<i<<" "<<material_ring[i]<<" "<<outer_radius_ring[i]<<"\n";
  };
  exit(0);
#endif


  /*
  vector<real> vol_part(4,0.);
  for(int i=0;i<sz;i++){
    real vol=PI*pow(outer_radius_ring[i],2);
    if(i!=0){
      vol-=PI*pow(outer_radius_ring[i-1],2);
    };
    cout<<i<<" "<<material_ring[i]<<" "<<outer_radius_ring[i]<<" "<<vol<<"\n";
    vol_part[material_ring[i]]+=vol;
  };

  for(int i=0;i<4;i++){
    cout<<i<<" "<<vol_part[i]<<"\n";
  };
  exit(0);
  */

  real factor_vpg=area_pellet/area_vpg; // This factor should be multiplied to pellet number densities

  // Cylinder geometric information is defined here.
  IrregularGeometryInformation ginp;
  real *rinp=new real[sz];
  int *regid=new int[sz];
  for(int i=0;i<sz;i++){
    rinp[i]=outer_radius_ring[sz-i-1];
    regid[i]=sz-i-1;
  };
  ginp.ExCircle(sz,rinp,regid);
  delete [] rinp;
  delete [] regid;


  for(int ii=0;ii<num_case;ii++){ // ii: 0/ordinary state, 1/voided state, 2/Doppler state

    //cout<<" Case : "<<ii<<"\n";
    // (Trajectory-set to calculate collision probabilities)
    TrajectorySet test;
    test.PutBoundaryCondition(White);
    test.CalTrajectory(ginp,1,0.02,45.); // Equidistant
    //ginp.WriteGnuplotFile(0.03); exit(0);
    
    // Material composition (nuclide number density) data transfer
    vector<Medium> med_ring(sz);
    for(int i=0;i<sz;i++){

      med_ring[i].PutImax(group);
      int plinp=pl;
      if(plinp==0)plinp=1;
      med_ring[i].PutPL(plinp);

      if(material_ring[i]==0){// vpg
        int num=sa.GetFuelComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){
  	  idinp[j]=sa.GetFuelComposition()->GetNucid(j);
	  deninp[j]=sa.GetFuelComposition()->GetDensity(j)*factor_vpg;
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==1){ // cladding
        int num=sa.GetCladComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){
	  idinp[j]=sa.GetCladComposition()->GetNucid(j);
	  deninp[j]=sa.GetCladComposition()->GetDensity(j);
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==2){ // sodium
        int num=1;
        int *idinp=new int[num];
        real *deninp=new real[num];
        idinp[0]=110230;
        deninp[0]=sa.GetSodiumNumberDensity();
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetSodiumTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==3){ // duct
        int num=sa.GetDuctComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){ 
	  idinp[j]=sa.GetDuctComposition()->GetNucid(j);
	  deninp[j]=sa.GetDuctComposition()->GetDensity(j);
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else{
        cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneousModel.\n";
        cout<<"# Material ring "<<material_ring[i]<<" cannot not be used.\n";
        exit(0);
      };
    };

    // The first medium instance, med_ring[0], should contain the nuclide data
    // in the same order as [matno] because of burup calculations.
    // So, the number density data for med_ring[0] are re-assigned in the following
    vector<real> den_med0(fuel_nuc_num,0.);
    for(int i=0;i<fuel_nuc_num;i++){
      int num_med0=med_ring[0].GetNucnum();
      for(int j=0;j<num_med0;j++){
        if(matno[i]==med_ring[0].GetNuclideInTurn(j).GetID()){
          den_med0[i]=med_ring[0].GetNuclideInTurn(j).GetDensity();
        };
      };
    };
    med_ring[0].PutNuclide(fuel_nuc_num,matno,den_med0);


    if(ii==1){
      // sodium-voided system
      for(int i=0;i<sz;i++){
        int nucnum=med_ring[i].GetNucnum();
        for(int j=0;j<nucnum;j++){
	  int nucid=med_ring[i].GetNuclideInTurn(j).GetMatnum();
	  if(nucid==110230){
	    med_ring[i].GetNuclideInTurn(j).PutDensity(sodium_density_in_voided_case);
	  };
        };
      };
    };

    if(ii==2){
      // Doppler system
      for(int i=0;i<sz;i++){
        int nucnum=med_ring[i].GetNucnum();
        for(int j=0;j<nucnum;j++){
  	  int nucid=med_ring[i].GetNuclideInTurn(j).GetMatnum();
          if(nucid>900000||nucid==420000||nucid==400000){
	  //if(nucid>900000){
	  //if(nucid>200000){
  	    real org=med_ring[i].GetNuclideInTurn(j).GetTemperature();
	    med_ring[i].GetNuclideInTurn(j).PutTemperature(org+delta_t_in_doppler_case);
	  };
        };
      };
    };

    // ... Self-shielding calculation
    SelfShieldingCalculator ssc;
    ssc.PutPijCalculator(&test);
    ssc.SetNumberOfMedium(sz);
    for(int i=0;i<sz;i++){
      ssc.PutMedium(i,med_ring[i]);
    };
    ssc.WithToneMethod(xslib,false); // [false] is for print option

    // ... Neutron transport calculation 
    int *region_medium=new int[sz];
    for(int i=0;i<sz;i++){region_medium[i]=i;};


#if 0
    if(ii==0){
      for(int i=0;i<sz;i++){
        cout<<"# Medium : "<<i<<"\n";
        ssc.GetMedium(i).ShowMacroCrossSection1D();
      };
      exit(0);
    };
#endif

#if 0
    int medid=4;
    int nnn=ssc.GetMedium(medid).GetNucnum();
    for(int i=0;i<nnn;i++){
      cout<<ssc.GetMedium(medid).GetNuclideInTurn(i).GetMatnum()<<"\n";
      GroupData1D enband=med_ring[0].GetEnband();
      ssc.GetMedium(medid).GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,1).show_self();
    };
    exit(0);
#endif

    // ... Eigenvalue calculation or critical buckling search
    PJISystem sys(group,sz);
    sys.NoPrint();
    sys.PutTrajectorySet(&test);
    for(int i=0;i<sz;i++){
      sys.AddMedium(ssc.GetMedium(i));
    };
    sys.PutRegMed(region_medium);
    GeneralOption opt;
    sys.PutGeneralOption(opt);
    sys.PutSigmaCol(sigtr);
    sys.PutPij();

    //sys.BucklingSearch("ModifiedSource");
    //sys.BucklingSearch("ModifiedSelfScattering");
    sys.BucklingSearch(); // ... dafault
    //real kinf=sys.CalIgen();

    delete [] region_medium;

    // ... Homogenization 
    sys.PutFluxAsCurrent();
    med_ref=sys.HomogenizeAll(false,true); // [false] for Benoist's D calculation option

    /*
    med_ref.CalHomoB1(1e-20);
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<med_ref.CalKinf()<<"\n";
    GroupData1D chitmp=med_ref.GetMacxs().GetData1d(chi);

    med_ref.CalMacroFromMicro();
    //med_ref.GetMacxs().GetData1d(chi)=chitmp;
    med_ref.FissionSpectrumVectorReconstruction();
    med_ref.CalHomoB1(1e-20);
    cout<<med_ref.CalKinf()<<"\n"; exit(0);
    */

    PutSADataToMZoneDefiningMediumData(med_ref, mzone, ii);

  }; // end of ii-loop


};

void FRBurnerRZ::PutSADataToMZoneHeterogeneous2D(FRDTSubAssembly &sa, int mzone, XSLibrary &xslib)
{
  int fuel_nuc_num=matno.size();

  Medium med_ref;
  
  int num_case=3;
  if(!xscalc_for_reactivity)num_case=1;

  real area_spacer_wire=sa.GetFuelSAGeometry()->GetSpacerwireVolume();
  if(area_spacer_wire>1e-10){
    cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneous2D.\n";
    cout<<"# Not yet coded for cases of non-zero spacer-wire volume.\n";
    exit(0);
  };

  int num_pin=sa.GetFuelSAGeometry()->GetNumberOfPin();

  real pin_pitch=sa.GetFuelSAGeometry()->GetPinPitch()*0.1; // [cm]
  real assembly_pitch=sa.GetFuelSAGeometry()->GetAssemblyPitch()*0.1; // [cm]
  real pin_outer_diameter=sa.GetFuelSAGeometry()->GetPinOuterDiameter()*0.1; // [cm]
  real pin_inner_diameter=pin_outer_diameter-sa.GetFuelSAGeometry()->GetPinThickness()*2*0.1; // [cm]
  real duct_outer_distance=sa.GetFuelSAGeometry()->GetDuctOuterSize()*0.1; // [cm]
  real duct_inner_distance=duct_outer_distance-sa.GetFuelSAGeometry()->GetDuctThickness()*2*0.1; // [cm]

  int layer=0;
  for(int i=1;i<20;i++){
    int num=3*i*(i-1)+1;
    if(num==num_pin)layer=i;
  };
  if(layer==0){
    cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneous2D.\n";
    cout<<"# The number of fuel layers cannot be determined.\n";
    exit(0);
  };

  /*
  cout<<num_pin<<"\n";
  cout<<pin_pitch<<"\n";
  cout<<assembly_pitch<<"\n";
  cout<<duct_outer_distance<<"\n";
  cout<<duct_inner_distance<<"\n";
  cout<<pin_outer_diameter<<"\n";
  cout<<pin_inner_diameter<<"\n";
  cout<<layer<<"\n";
  exit(0);
  */

  IrregularGeometryInformation ginp;

  real scon=1./sqrt(3.0);

  real *hexsize=new real[layer+2];
  int *hexreg=new int[layer+2];
  hexsize[0]=assembly_pitch*scon;
  hexsize[1]=duct_outer_distance*scon;
  hexsize[2]=duct_inner_distance*scon;
  for(int i=0;i<layer-1;i++){
    hexsize[3+i]=pin_pitch*(layer-1-i);
  };

  // (simple model: one material per one layer)

  int matnum=layer*2;

  for(int i=0;i<layer+2;i++){
    hexreg[layer+1-i]=matnum+i;
  };

  ginp.AddHexagonLayer(layer+2,hexsize,hexreg);
  delete [] hexsize;
  delete [] hexreg;

  real cring[]={pin_outer_diameter*0.5, pin_inner_diameter*0.5};
  int regj[]={1,0};
  ginp.AddCircleRing(2,cring,regj);
  for(int i=0;i<layer-1;i++){
    regj[0]=i*2+3;
    regj[1]=i*2+2;
    ginp.AddHexagonRing(pin_pitch*(i+1),i+1,2,cring,regj);
  };
  //ginp.WriteGnuplotFile(0.03); exit(0);


  // (rigorous model: several materials per one layer)
  /*
  int matnum=0;
  int mattmp=1;
  for(int i=0;i<layer;i++){
    matnum+=mattmp*2;
    if((i+1)%2==0)mattmp+=1;    
  };

  for(int i=0;i<layer+2;i++){
    hexreg[layer+1-i]=matnum+i;
  };

  ginp.AddHexagonLayer(layer+2,hexsize,hexreg);
  delete [] hexsize;
  delete [] hexreg;

  real cring[]={pin_outer_diameter*0.5, pin_inner_diameter*0.5};
  int regj[]={1,0};
  ginp.AddCircleRing(2,cring,regj);
  int tmp1=1;
  int tmp2=2;
  for(int i=1;i<layer;i++){
    regj[0]=tmp2+1;
    regj[1]=tmp2;
    ginp.AddHexagonRingFine(pin_pitch*i,i,2,cring,regj);
    tmp2+=tmp1*2;
    if((i+1)%2==0)tmp1+=1;
  };
  //ginp.WriteGnuplotFile(0.03); exit(0);
  */


  int sz=layer+2+matnum;

  vector<int> material_ring(sz); // 0:v+p+g, 1:cladding, 2:sodium, 3:duct
  for(int i=0;i<matnum/2;i++){
    material_ring[i*2]=0;
    material_ring[i*2+1]=1;
  };
  for(int i=0;i<layer;i++){
    material_ring[matnum+i]=2;
  };
  material_ring[matnum+layer]=3;
  material_ring[matnum+layer+1]=2;


  real area_pellet=sa.GetFuelSAGeometry()->GetPelletVolume()/num_pin;
  real area_pellet_invoid=sa.GetFuelSAGeometry()->GetPelletVoidVolume()/num_pin;
  real area_pellet_gap=sa.GetFuelSAGeometry()->GetPelletGapVolume()/num_pin;

  real area_vpg=area_pellet_invoid+area_pellet+area_pellet_gap;
  real factor_vpg=area_pellet/area_vpg; // This factor should be multiplied to pellet number densities

  for(int ii=0;ii<num_case;ii++){ // ii: 0/ordinary state, 1/voided state, 2/Doppler state
    
    // (Trajectory-set to calculate collision probabilities)
    TrajectorySet test_ss; // self-shielding
    //test_ss.PutBoundaryCondition(White);
    test_ss.PutBoundaryCondition(Periodic);
    test_ss.CalTrajectory(ginp, 8, 0.02, 60.);

    TrajectorySet test_fc; // flux-calculation
    //test_fc.PutBoundaryCondition(White);
    test_fc.PutBoundaryCondition(Periodic);
    test_fc.CalTrajectory(ginp, 8, 0.02, 60.);

    // Material composition (nuclide number density) data transfer
    vector<Medium> med_ring(sz);
    for(int i=0;i<sz;i++){

      med_ring[i].PutImax(group);
      int plinp=pl;
      if(plinp==0)plinp=1;
      med_ring[i].PutPL(plinp);

      if(material_ring[i]==0){// vpg
        int num=sa.GetFuelComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){
  	  idinp[j]=sa.GetFuelComposition()->GetNucid(j);
	  deninp[j]=sa.GetFuelComposition()->GetDensity(j)*factor_vpg;
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==1){ // cladding
        int num=sa.GetCladComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){
	  idinp[j]=sa.GetCladComposition()->GetNucid(j);
	  deninp[j]=sa.GetCladComposition()->GetDensity(j);
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==2){ // sodium
        int num=1;
        int *idinp=new int[num];
        real *deninp=new real[num];
        idinp[0]=110230;
        deninp[0]=sa.GetSodiumNumberDensity();
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetSodiumTemperature());
        delete [] idinp;
        delete [] deninp;
      }else if(material_ring[i]==3){ // duct
        int num=sa.GetDuctComposition()->GetNucnum();
        int *idinp=new int[num];
        real *deninp=new real[num];
        for(int j=0;j<num;j++){ 
	  idinp[j]=sa.GetDuctComposition()->GetNucid(j);
	  deninp[j]=sa.GetDuctComposition()->GetDensity(j);
        };
        med_ring[i].PutNuclide(num,idinp,deninp);
        med_ring[i].PutTemperatureForAllNuclide(sa.GetTemperature());
        delete [] idinp;
        delete [] deninp;
      }else{
        cout<<"# Error in FRBurnerRZ::PutSADataToMZoneHeterogeneousModel.\n";
        cout<<"# Material ring "<<material_ring[i]<<" cannot not be used.\n";
        exit(0);
      };
    };

    // The first medium instance, med_ring[0], should contain the nuclide data
    // in the same order as [matno] because of burup calculations.
    // So, the number density data for med_ring[0] are re-assigned in the following
    vector<real> den_med0(fuel_nuc_num,0.);
    for(int i=0;i<fuel_nuc_num;i++){
      int num_med0=med_ring[0].GetNucnum();
      for(int j=0;j<num_med0;j++){
        if(matno[i]==med_ring[0].GetNuclideInTurn(j).GetID()){
          den_med0[i]=med_ring[0].GetNuclideInTurn(j).GetDensity();
        };
      };
    };
    med_ring[0].PutNuclide(fuel_nuc_num,matno,den_med0);


    if(ii==1){
      // sodium-voided system
      for(int i=0;i<sz;i++){
        int nucnum=med_ring[i].GetNucnum();
        for(int j=0;j<nucnum;j++){
	  int nucid=med_ring[i].GetNuclideInTurn(j).GetMatnum();
	  if(nucid==110230){
	    med_ring[i].GetNuclideInTurn(j).PutDensity(sodium_density_in_voided_case);
	  };
        };
      };
    };

    if(ii==2){
      // Doppler system
      for(int i=0;i<sz;i++){
        int nucnum=med_ring[i].GetNucnum();
        for(int j=0;j<nucnum;j++){
  	  int nucid=med_ring[i].GetNuclideInTurn(j).GetMatnum();
          if(nucid>900000||nucid==420000||nucid==400000){
	  //if(nucid>900000){
	  //if(nucid>200000){
  	    real org=med_ring[i].GetNuclideInTurn(j).GetTemperature();
	    med_ring[i].GetNuclideInTurn(j).PutTemperature(org+delta_t_in_doppler_case);
	  };
        };
      };
    };

    // ... Self-shielding calculation
    SelfShieldingCalculator ssc;
    ssc.PutPijCalculator(&test_ss);
    ssc.SetNumberOfMedium(sz);
    for(int i=0;i<sz;i++){
      ssc.PutMedium(i,med_ring[i]);
    };
    ssc.WithToneMethod(xslib,false); // [false] is for print option

    // ... Neutron transport calculation 
    int *region_medium=new int[sz];
    for(int i=0;i<sz;i++){region_medium[i]=i;};

    // ... Eigenvalue calculation or critical buckling search
    PJISystem sys(group,sz);
    sys.NoPrint();
    sys.PutTrajectorySet(&test_fc);
    for(int i=0;i<sz;i++){
      sys.AddMedium(ssc.GetMedium(i));
    };
    sys.PutRegMed(region_medium);
    GeneralOption opt;
    sys.PutGeneralOption(opt);
    sys.PutSigmaCol(sigtr);
    sys.PutPij();

    //sys.BucklingSearch("ModifiedSelfScattering");
    sys.BucklingSearch();
    //real kinf=sys.CalIgen();

    delete [] region_medium;

    // ... Homogenization 
    sys.PutFluxAsCurrent();
    med_ref=sys.HomogenizeAll(false,true); // [false] for Benoist's D calculation option

    PutSADataToMZoneDefiningMediumData(med_ref, mzone, ii);

  }; // end of ii-loop


};

void FRBurnerRZ::PutSADataToMZoneDefiningMediumData(Medium &med_ref, int mzone, int ii)
{
    med_ref.CalMuFromSigel_p1();

    for(int i=0;i<mednum;i++){
      if(mzone_per_med[i]==mzone){
        med_data_inp[i]=true;
        if(i==mzone_rep[mzone]){
          if(ii==0)med[i]=med_ref;
          if(ii==1)med_voided[i]=med_ref;
          if(ii==2)med_doppler[i]=med_ref;
        }else{
          if(ii==0)med[i].CopyWithoutMicroscopicScatteringMatrices(med_ref);
          if(ii==1)med_voided[i].CopyWithoutMicroscopicScatteringMatrices(med_ref);
          if(ii==2)med_doppler[i].CopyWithoutMicroscopicScatteringMatrices(med_ref);
	};
      };
    };

};

void FRBurnerRZ::PutNonfuelInfo(int medid, int nucnum, int *nuc, real *den, real temp)
{
  med[medid].PutImax(group);
  int plinp=pl;
  if(plinp==0)plinp=1;
  med[medid].PutPL(plinp);
  med[medid].PutNuclide(nucnum,nuc,den);
  med[medid].PutTemperatureForAllNuclide(temp);

  med_data_inp[medid]=true;
};

void FRBurnerRZ::PutNonfuelInfo(int medid, FRDTComposition &inp, real temp)
{
  int nucinp=inp.GetNucnum();
  int *idinp=new int[nucinp];
  real *deninp=new real[nucinp];
  for(int i=0;i<nucinp;i++){
    idinp[i]=inp.GetNucid(i);
    deninp[i]=inp.GetDensity(i);
  };
  PutNonfuelInfo(medid,nucinp,idinp,deninp,temp);
  delete [] idinp;
  delete [] deninp;
};

void FRBurnerRZ::CalMacroFromMicroBurnupMedium()
{
  /*
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
    //med[i].CalMacroFromMicroSimple();
  };
  */

  for(int i=0;i<num_mzone;i++){
    int id=mzone_rep[i]; // ID for representative medium instance containing microscopic XS
    for(int jj=med_burn-1;jj>=0;jj--){
      if(mzone_per_med[jj]==i){
        for(int j=0;j<burn_nuc;j++){
          med[id].GetNuclideInTurn(j).PutDensity(density_ave[jj][j]);
	};
	med[id].CalMacroFromMicro();
	if(jj!=id){
  	  med[jj].GetMacxs()=med[id].GetMacxs();
	};
      };
    };

  };
};

GroupData1D FRBurnerRZ::GetIntegratedFluxPerMZone(int i)
{
  GroupData1D flx(group);
  flx.set_zero();

  for(int i=0;i<med_burn;i++){
    if(mzone_per_med[i]==i){
      real vol=vol_per_med[i];
      flx=flx+med[i].GetFlux()*vol;
    };
  };

  return flx;
};




void FRBurnerRZ::RunADS(XSLibrary &xslib, Burnup &bu, bool print, bool denout)
{
  // Functions for Ogawa-kun's work were deleted.
  // If necessary, this can be recovered from the old source code in 2020/6/28.

  //bool heat_by_capture=false; 
  //bool heat_by_nonfuel=false;

  bool heat_by_capture=true; 
  bool heat_by_nonfuel=true;

  bool pc_cal=false; // Predictor-corrector calculation

  // ! NOTE !
  //
  // - Heating from Non-fuel regions are ignored.
  // - Heating by capture reactions are ignored.

  // +++ Medium-wise volume and heavy metal weight calculation +++++++

  // This instance is dummy to calculate medium-wise volume
  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
  };
  test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
 
  real vol_fuel_total=0.; // Total volume of fuel mediums
  vector<real> vol_fuel_total_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    vol_fuel_total_mzone[i]=0.;
  };
  for(int i=0;i<med_burn;i++){
    vol_per_med[i]=test.GetVolumePerMedium(i);                   // Volume of medium-i
    real weight=bu.CalWeightOfHeavyNuclideParUnitVolume(med[i]); // Heavy metal weight per unit volume of medium-i
    thm_per_med[i]=weight*vol_per_med[i];                        // Heavy metal weight of medium-i
    vol_fuel_total+=vol_per_med[i];
    vol_fuel_total_mzone[mzone_per_med[i]]+=vol_per_med[i];
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Material zone-wise and Total heavy metal weight calculation +++++++++++++++++++++++++
  vector<real> hm_w_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    hm_w_mzone[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(med[mzone_rep[i]]);
  };

  real hm_w_total=0.;
  vector<real> hm_w_total_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    hm_w_total_mzone[i]=0.;
  };
  for(int i=0;i<med_burn;i++){
    int mz=mzone_per_med[i];
    real tmp=hm_w_mzone[mz]*vol_per_med[i];
    hm_w_total+=tmp;
    hm_w_total_mzone[mz]+=tmp;
  };

  cout<<"#######################################################\n";
  cout<<"# Total weight of heavy metal : "<<hm_w_total<<"[g]\n";
  for(int i=0;i<num_mzone;i++){
    cout<<"#         ("<<name_mzone[i]<<")        : "<<hm_w_total_mzone[i]<<"[g]\n";
  };
  cout<<"#######################################################\n";

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Instance [newfuel] of the medium class is to store number densities of initially-loaded MA nuclides (Np, Am, Cm).
  // Compositions of fuel loaded at the following cycles are based on this number density data. 
  // This composition is assumed to be those taken from LWR spent fuels.
  Medium newfuel; 
  vector<int> matid(burn_nuc);
  vector<real> matden(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    matid[i]=med[0].GetNuclideInTurn(i).GetMatnum();
    matden[i]=0.;
    if((matid[i]>930000&&matid[i]<940000)||
       (matid[i]>950000&&matid[i]<970000)){
      // Number densities of only Np, Am, Cm are stored.
      matden[i]=med[0].GetNuclideInTurn(i).GetDensity();
    };
  };
  newfuel.PutNuclide(burn_nuc,matid,matden);

  real hm_w_newfuel=bu.CalWeightOfHeavyNuclideParUnitVolume(newfuel)*vol_fuel_total;
  vector<real> hm_w_newfuel_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    hm_w_newfuel_mzone[i]=bu.CalWeightOfHeavyNuclideParUnitVolume(newfuel)*vol_fuel_total_mzone[i];
  };
  //  heavy metal weight of MA nuclides in initially-loaded fuel
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Reloading matrix&vector construction +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  int szmat=med_burn*burn_nuc;
  GroupData1D vec_density(szmat);
  GroupData2D mat_reload(szmat,szmat);
  GroupData1D vec_weight_cal(szmat);
  // conversion factor from number density per unit volume to volume-integrated weight
  GroupData1D vec_reload(szmat);

  //FuelExchangeOperatorSetting(mat_reload,vec_weight_cal,vec_reload);

  mat_reload.set_zero();
  vec_weight_cal.set_zero();
  vec_reload.set_zero();

  for(int i=0;i<szmat;i++){
    mat_reload.put_data(i,i,-1.);  // Extraction component
  };

#if 0
  // One-region refueling scheme
  for(int i=0;i<burn_nuc;i++){
    int mid=med[0].GetNuclideInTurn(i).GetMatnum();

    for(int j=0;j<med_burn;j++){
      int pos_y=j*burn_nuc+i; 
      real factor=matden[i]/hm_w_newfuel;
      vec_reload.add_data(pos_y,factor*hm_w_total);
      for(int k=0;k<burn_nuc;k++){
	int mid2=med[0].GetNuclideInTurn(k).GetMatnum();
	if(mid2>=920000&&mid2<=1000000){ // for heavy metal weight calculation
	  for(int l=0;l<med_burn;l++){
	    int pos_x=l*burn_nuc+k;
  	    real aw=bu.GetAtomicWeight(mid2);
            real wgt_factor=1./avo*aw; // [g/#]
	    real val=-factor*vol_per_med[l]*wgt_factor;
            mat_reload.add_data(pos_y,pos_x,val); 
	    if(i==0&&j==0)vec_weight_cal.add_data(pos_x,vol_per_med[l]*wgt_factor);
	  };
	};
      };
    };

    if(mid>=920000&&mid<=1000000){ // Actinoid nuclides are reused in the next cycle
      for(int j=0;j<med_burn;j++){
        int pos_x=j*burn_nuc+i; 
	real val=vol_per_med[j]/vol_fuel_total; // Reloading component
	for(int k=0;k<med_burn;k++){
	  int pos_y=k*burn_nuc+i;
	  mat_reload.add_data(pos_y,pos_x,val);
	};
      };
    };
  };
#endif

#if 1
  // Material zone-wise refueling scheme
  // If radial blanket is considerd, this option cannot be used.

  for(int i=0;i<burn_nuc;i++){
    int mid=med[0].GetNuclideInTurn(i).GetMatnum();

    for(int j=0;j<med_burn;j++){
      int mz=mzone_per_med[j];
      int pos_y=j*burn_nuc+i;
      real factor=matden[i]/hm_w_newfuel_mzone[mz];
      vec_reload.add_data(pos_y,factor*hm_w_total_mzone[mz]);
      for(int k=0;k<burn_nuc;k++){
	int mid2=med[j].GetNuclideInTurn(k).GetMatnum();
	if(mid2>=920000&&mid2<=1000000){ // for heavy metal weight calculation
	  for(int l=0;l<med_burn;l++){
	    int mz2=mzone_per_med[l];
	    if(mz==mz2){
  	      int pos_x=l*burn_nuc+k;
  	      real aw=bu.GetAtomicWeight(mid2);
              real wgt_factor=1./avo*aw; // [g/#]
	      real val=-factor*vol_per_med[l]*wgt_factor;
              mat_reload.add_data(pos_y,pos_x,val); 
	      if(i==0&&j==mzone_rep[mz])vec_weight_cal.add_data(pos_x,vol_per_med[l]*wgt_factor);
	    };
	  };
	};
      };
    };

    if(mid>=920000&&mid<=1000000){ // Actinoid nuclides are reused in the next cycle
      for(int j=0;j<med_burn;j++){
	int mz=mzone_per_med[j];
        int pos_x=j*burn_nuc+i; 
  	real val=vol_per_med[j]/vol_fuel_total_mzone[mz]; // Reloading component
	for(int k=0;k<med_burn;k++){
	  int mz2=mzone_per_med[k];
          if(mz==mz2){
            int pos_y=k*burn_nuc+i;
            mat_reload.add_data(pos_y,pos_x,val);
	  };
	};
      };
    };

  }; // end of loop "i"
#endif


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Fixed-source calculation option
  real epsf=1e-4;
  int itermax=1000;
  bool high_speed_option=true;

  // +++ External neutron source reading
  vector< vector<GroupData1D> > esrc;
  ExternalSourceReading(esrc);

  // (initialize number of irradiated cycles and accumulated burnup)
  for(int i=0;i<med_burn;i++){
    int sz=ir_cycle[i].size();
    for(int b=0;b<sz;b++){
      ir_cycle[i][b]=0;
      ac_burn[i][b]=0.;
    };
  };

  // +++
  //density_store.resize(trace_cycle);
  fwd_nuc.resize(trace_cycle);

  // +++
  GeneralOption opt;

  vector<int> refuel_id_mzone(num_mzone);
  for(int i=0;i<num_mzone;i++){
    refuel_id_mzone[i]=0;
  };

  bu.PutNucnum(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    bu.PutNuclideData(i,med[0].GetNuclideInTurn(i).GetMatnum(),0.,0.,0.,0.);
  };
  bu.CalTransitionMatrixFluxInDependentPart();

  flux_level.resize(trace_cycle);
  keff.resize(trace_cycle);

  string kout_name="Ksub";
  cout<<"#C Step  Day  Keff   Max. power   Material zone-wise      ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  cout<<"\n";
  cout<<"#y                   dens[W/cm3]  Power Dist. [%]         ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  cout<<"\n";
  cout<<"#                    (pos:r,z)                         ";
  if(num_mzone>4){
    for(int i=0;i<num_mzone-4;i++)cout<<"     ";
  };
  cout<<"\n";
  cout<<"#                                  ";
  for(int i=0;i<num_mzone;i++){
    cout<<"("<<i<<")  ";
  };
  cout<<"\n";


  // (cycle-loop)
  real ac_day=0.;
  for(int cycle=0;cycle<trace_cycle;cycle++){

    flux_level[cycle].resize(cycle_div+1);
    keff[cycle].resize(cycle_div+1);    
    //density_store[cycle].resize(cycle_div+1);
    fwd_nuc[cycle].resize(cycle_div+1);

    // +++ Refueling after the second cycle +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(cycle!=0){

      // (Number density data just before refueling are stored in one-dimensional vector [vec_density])
      int index=0;
      for(int i=0;i<med_burn;i++){
	for(int j=0;j<burn_nuc;j++){
	  vec_density.put_data(index,density[i][0][j]);
	  index++;
	};
      };

      real weight_at_present=vec_weight_cal*vec_density;
      // Total heavy metal weight

      real factor=(hm_w_total-weight_at_present)/hm_w_newfuel;
      
      if(print){
        cout<<"#\n# Total weight of newly-loaded MA materials [kg]\n#\n";
        cout<<"#    (Cycle : "<<cycle<<")\n#\n";
        real sum=0.;
        for(int i=0;i<burn_nuc;i++){
	  if(matden[i]>0.){
	    int mid=med[0].GetNuclideInTurn(i).GetMatnum();
            real aw=bu.GetAtomicWeight(mid);
            real wgt=(matden[i]*factor*vol_fuel_total)/avo*aw*1e-3; // [kg]
	    sum+=wgt;
	    cout<<"#  ";
  	    WriteOut(midt.Name(mid),6);
	    cout<<" ";
	    WriteOut(wgt,"%8.2f");
	    cout<<"\n";
	  };
        };
        cout<<"#  Total  ";
        WriteOut(sum,"%8.2f");
        cout<<"\n#\n";

	cout<<"# Total weight of remaining MA materials [kg]\n#\n";
	sum=0.;
	for(int i=0;i<burn_nuc;i++){
	  int mid=med[0].GetNuclideInTurn(i).GetMatnum();
	  if(mid>900000&&mid<1000000){
  	    real mass=0;
            for(int j=0;j<med_burn;j++){
	      mass+=bu.GetAtomicWeight(mid)/avo*vol_per_med[j]*density[j][0][i]*1e-3;
	    };
	    sum+=mass;
	    cout<<"#  "<<mid<<" ";
	    WriteOut(mass,"%8.2f");
	    cout<<" "<<"\n";
	  };
	};
	cout<<"#  Total  ";
        WriteOut(sum,"%8.2f");
        cout<<"\n#\n";
	cout<<"#\n";
      };
      /*
      for(int j=0;j<burn_nuc;j++){
      real sum1=0.;
      real sum2=0.;
      for(int i=0;i<med_burn;i++){
	real vol=vol_per_med[i];
	real den=vec_density.get_dat(i*burn_nuc+j);
	sum1+=vol*den;
	sum2+=vol;
	//cout<<i<<" "<<den<<" "<<vol<<"\n";
      };
      WriteOut(midt.Name(med[0].GetNuclideInTurn(j).GetMatnum()),8);
      cout<<sum1/sum2<<"\n";
      };
      exit(0);
      */

      vec_density=vec_density+vec_reload+mat_reload*vec_density;

      /*
      weight_at_present=vec_weight_cal*vec_density;
      cout<<"# Weight check.\n#\n";
      cout<<"#   Expected : "<<hm_w_total<<"\n";
      cout<<"#   Actual   : "<<weight_at_present<<"\n";
      */
      /*
      real tmp=0.;
      for(int i=0;i<(ic_e+1)*burn_nuc;i++){
	tmp+=vec_weight_cal.get_dat(i)*vec_density.get_dat(i);
      };
      cout<<"        inner core : "<<tmp<<"\n";
      */
      /*
      cout<<"# Newly-loaded fuel composition\n#\n";
      for(int j=0;j<burn_nuc;j++){
	WriteOut(midt.Name(med[0].GetNuclideInTurn(j).GetMatnum()),6);
	cout<<" "<<vec_density.get_dat(j)<<"\n";
      };
      */

      index=0;
      for(int i=0;i<med_burn;i++){
	for(int j=0;j<burn_nuc;j++){
	  density[i][0][j]=vec_density.get_dat(index);
	  index++;
	};
      };

    };
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // (in cycle-loop)
    for(int tim=0;tim<cycle_div+1;tim++){

      flux_level[cycle][tim].resize(mednum);
      real delt_d=cycle_length/cycle_div; // [d]
      real delt=delt_d*24*60*60; // [s]

      // +++ Calculation of region-averaged density
      CalRegionAveragedDensity();
      StoreRegionAveragedDensity(cycle,tim);

#if 0
      if(denout){
        if(cycle!=0&&tim==0){
          med[ic_s].WriteFileNumberDensity("./"+dirname_boc+"/","den_ic_cyc"+IntToString(cycle));
          med[oc_s].WriteFileNumberDensity("./"+dirname_boc+"/","den_oc_cyc"+IntToString(cycle));
          med[rb_s].WriteFileNumberDensity("./"+dirname_boc+"/","den_rb_cyc"+IntToString(cycle));
        };
        if(tim==cycle_div){
	  for(int i=0;i<med_burn;i++){
	    med[i].WriteFileNumberDensity("./"+dirname_eoc+"/","den_cyc"+IntToString(cycle)+"_med"+IntToString(i));
	  };
        };
      };
#endif

      // +++ Macroscopic cross section calculation
      CalMacroFromMicroBurnupMedium();

      // ++ Eigenvalue calculation
      PLOSSystem test_eg(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
        test_eg.AddMedium(med[i]);
        test_eg.GetMedium(i).MicxsVectorClear2DData();
      };
      test_eg.NoPrint();
      test_eg.PutCartMeshInfo(cmi,"Cylinder");
      test_eg.PutGeneralOption(opt);
      test_eg.CalCoef();
      string optname="cmfd";
      if(!cmfd_on)optname="";
      real kout=test_eg.CalIgen(optname);

      // +++ Fixed-source calculation ++++++++++++++++++++++++++++++++++++

      PLOSSystem test(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
        test.AddMedium(med[i]);
	test.GetMedium(i).MicxsVectorClear2DData();
      };
      test.NoPrint();
      test.PutCartMeshInfo(cmi,"Cylinder");
      test.PutGeneralOption(opt);
      test.CopyCoef(test_eg);

      ExternalSourceSetting(test,esrc);

      test.CalFixedSourceWithFission(epsf,itermax,high_speed_option);

      if(show_fission_info)test.GetNeutronMultiplicationInfo(kout);


      for(int i=0;i<mednum;i++){
        GroupData1D flx=test.GetIntegratedFlux(i);
        med[i].GetFlux().copy(flx);
      };

      vector< vector<real> > pow_store(med_burn);
      // Volume-integrated power for each medium/batch
      // Note : "Volume" is defined as "total" one for all the batches

      // (power normalization)
      real pow_sum=0.;
      for(int i=mednum-1;i>=0;i--){
	//for(int i=0;i<mednum;i++){
        GroupData1D flx=med[i].GetFlux();
	if(i<med_burn){
          int mz=mzone_per_med[i];
	  int repid=mzone_rep[mz];
	  int sz=ac_burn[i].size(); // batch-wise
	  real pow=0.;
          pow_store[i].resize(sz);
	  for(int j=0;j<sz;j++){ 
	    for(int k=0;k<burn_nuc;k++){
	      med[repid].GetNuclideInTurn(k).PutDensity(density[i][j][k]);
	    };
            real tmp=bu.GetIntegratedPower(med[repid],flx,heat_by_capture);
            pow_store[i][j]=tmp;
	    //cout<<i<<" "<<pow_store[i][j]<<"\n";
	    pow+=tmp;
	  };
	  pow_sum+=pow/sz;
	}else if(heat_by_nonfuel){
	  pow_sum+=bu.GetIntegratedPower(med[i]);
	};
      };
      real factor=power/pow_sum;

      // +++ Power density calculation
      int yf=cmi.GetYF();
      int xf=cmi.GetXF();
      vector< vector<real> > power_map(yf);
      int index=0;
      for(int y=0;y<yf;y++){
	power_map[y].resize(xf);
	for(int x=0;x<xf;x++){
          int medid=cmi.GetFMat(index);
	  if(heat_by_nonfuel||medid<med_burn){
            int mz=mzone_per_med[medid];
  	    real pow=bu.GetIntegratedPowerParMesh(test,x,y,0,heat_by_capture);
 	    real vol=test.GetMesh(x,y,0).GetVolume();
	    pow=pow/vol; // power density [power/volume]
	    //real zl=test.GetMesh(x,y,0).GetLen(1);
	    //real area=vol/zl;
            //int pinn=num_pin_mzone[mz];
	    //if(pinn==0)pinn=num_pin_mzone[0];
	    //pow=(pow/zl)/(pinn*area/area_assembly_mzone[mz]); // line-power
	    power_map[y][x]=pow*factor;
	  };
          index++;
	};
      };

      // +++ Maximum line power calculation
      int maxp_x,maxp_y;
      real maxpow=MaximumPowerDensityCalculation(power_map,maxp_x,maxp_y);
      if(print_linepower_map){
        cout<<"#--------------------------------------------------------------------\n";
        cout<<"# Power density map in fine mesh [W/cm3]\n";
        //cout<<"\n# Line power map in fine mesh [W/cm]\n";
        PrintPowerDensity(power_map);
      };

      // +++ Printing summary of calculation results
      WriteOut(cycle,2);
      WriteOut(tim,4);
      int iac_day=int(ac_day);
      WriteOut(iac_day,6);
      cout<<" ";
      WriteOut(kout,"%8.6f");
      cout<<" ";
      WriteOut(maxpow,"%5.1f");
      cout<<" ("<<maxp_x<<","<<maxp_y<<")";

      // (power distribution)
      vector<real> pow_reg(num_mzone,0.);
      for(int i=0;i<med_burn;i++){
        int sz=ac_burn[i].size();
	for(int j=0;j<sz;j++){
  	  pow_reg[mzone_per_med[i]]+=pow_store[i][j]/sz;
	};
      };
      for(int i=0;i<num_mzone;i++){
	pow_reg[i]*=(100/pow_sum);
	WriteOut(pow_reg[i],"%5.1f");
      };

      cout<<"\n";

      for(int m=0;m<mednum;m++){
        real sumflx=med[m].GetFlux().get_sum()*factor/vol_per_med[m];
        flux_level[cycle][tim][m]=sumflx;
      };
      real factor_org=factor;

      // +++ Burnup calculation
      if(tim!=cycle_div){

	// +++ Predictor-corrector +++++++++++++++++++++++++++++++++++++
	if(pc_cal){

	vector< vector< vector<real> > > density_tmp=density;
	vector<GroupData1D> flx_tmp(mednum);
	for(int i=0;i<mednum;i++){
	  flx_tmp[i]=med[i].GetFlux(); // initial flux
	};

	// (predictor calculation)
        for(int m=0;m<med_burn;m++){ // fuel medium
          int mz=mzone_per_med[m];
          GroupData1D flx=med[m].GetFlux();
	  bu.CalOneGroupCrossSection(med[mzone_rep[mz]],flx);
          bu.CalTransitionMatrixFluxDependentPart();
          int tmp=GetBatchNumberFromMediumID(m);
          for(int j=0;j<tmp;j++){
	    bu.PutDensity(density[m][j]);
            bu.BurnupCalculation(flux_level[cycle][tim][m],delt,"chebyshev",false);
	    density[m][j]=bu.GetDensity();
          };
	};

	CalRegionAveragedDensity();

        // +++ Macroscopic cross section calculation
        CalMacroFromMicroBurnupMedium();

        PLOSSystem test_cor(dim,group,mednum); 
        for(int i=0;i<mednum;i++){
          test_cor.AddMedium(med[i]);
    	  test_cor.GetMedium(i).MicxsVectorClear2DData();
        };
        test_cor.NoPrint();
        test_cor.PutCartMeshInfo(cmi,"Cylinder");
        test_cor.PutGeneralOption(opt);
        test_cor.CalCoef();

        ExternalSourceSetting(test_cor,esrc);

        test_cor.CalFixedSourceWithFission(epsf,itermax,high_speed_option);

	for(int i=0;i<mednum;i++){
	  GroupData1D flx=test_cor.GetIntegratedFlux(i); // corrector flux
	  flx=flx+flx_tmp[i];
	  flx=flx*0.5;
	  med[i].GetFlux().copy(flx);
	};

        density=density_tmp;

        // (power normalization)
        real pow_sum=CalTotalPower(bu,pow_store,heat_by_nonfuel,heat_by_capture);
        real factor=power/pow_sum;

        for(int m=0;m<mednum;m++){
          real sumflx=med[m].GetFlux().get_sum()*factor/vol_per_med[m];
          flux_level[cycle][tim][m]=sumflx;
        };

	};
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	for(int i=0;i<med_burn;i++){
	  int sz=ac_burn[i].size();
	  for(int j=0;j<sz;j++){
	    ac_burn[i][j]+=((pow_store[i][j]*factor*1e-9)*delt_d)/(thm_per_med[i]*1e-6); // [GWd/t]
	  };
	};

        ac_day+=delt_d;
        for(int m=0;m<med_burn;m++){ // fuel medium

          GroupData1D flx=med[m].GetFlux();
          int repid=mzone_rep[mzone_per_med[m]];
	  bu.CalOneGroupCrossSection(med[repid],flx);
	  for(int i=0;i<burn_nuc;i++){
	    sigf_1g[m][cycle][tim][i]=bu.GetSigf(i);
            sigc_1g[m][cycle][tim][i]=bu.GetSigc(i);
            sign2n_1g[m][cycle][tim][i]=bu.GetSign2n(i);
	  };

          bu.CalTransitionMatrixFluxDependentPart();

          int tmp=GetBatchNumberFromMediumID(m);
          for(int j=0;j<tmp;j++){
	    bu.PutDensity(density[m][j]);
            bu.BurnupCalculation(flux_level[cycle][tim][m],delt,"chebyshev",false);
	    density[m][j]=bu.GetDensity();
          };

	};

      }else{
        // +++ end of cycle: cooling calculation
        bool decay_cal=false;
        real decay_cal_day=refuel_day;
	if(cycle!=trace_cycle-1&&decay_cal_day>0.){
	  decay_cal=true;
	}else{
	  if(cooling_day>0.){
	    decay_cal=true;
	    decay_cal_day=cooling_day;
	  };
	};
	if(decay_cal){
 	  for(int m=0;m<med_burn;m++){
	    int sz=ac_burn[m].size();
	    for(int i=0;i<sz;i++){
	      for(int j=0;j<burn_nuc;j++){
	        bu.PutDensity(j,density[m][i][j]);
	      };
              bu.BurnupCalculation(1e-15,decay_cal_day*24*60*60,"chebyshev",false);
	      for(int j=0;j<burn_nuc;j++){
	        density[m][i][j]=bu.GetDensity(j);
	      };
	    };
	  };
	};

      };

    }; // (incycle-loop end)

    // (increase the number of irradiated cycles)
    for(int i=0;i<med_burn;i++){
      int sz=ir_cycle[i].size();
      for(int j=0;j<sz;j++){
	ir_cycle[i][j]+=1;
      };
    };

  }; // (cycle-loop end)

  CalRegionAveragedDensity();
};

real FRBurnerRZ::GetKeff(int cycin, int substepin)
{
  if(cycin<0||cycin>=trace_cycle){
    cout<<"# Error in FRBurnerRZ::GetKeff\n";
    cout<<"# Inappropriate cycle "<<cycin<<" is requiested.\n";
    exit(0);
  };
  if(substepin<0||substepin>cycle_div){
    cout<<"# Error in FRBurnerRZ::GetKeff\n";
    cout<<"# Inappropriate substep "<<substepin<<" is requiested.\n";
    exit(0);
  };
  return keff[cycin][substepin];
};



