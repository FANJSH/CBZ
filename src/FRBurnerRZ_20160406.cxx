#include <cstdlib>
#include "FRBurnerRZ.h"

FRBurnerRZ::FRBurnerRZ()
{
  dim    = 2;
  group  = 70;
  refuel_day=0.;
  cooling_day=0.;
  
  nic=0;
  noc=0;
  nrb=0;
  nabic=0;
  naboc=0;
 
  ic_s=-1;
  ic_e=-2;
  oc_s=-1;
  oc_e=-2;
  rb_s=-1;
  rb_e=-2;
  abic_s=-1;
  abic_e=-2;
  aboc_s=-1;
  aboc_e=-2;

  void_cal=false;
  dop_cal=false;

  print_linepower_map=false;
  cmfd_on=true;

  show_fission_info=false;
  ogawa_cal=false;
};

void FRBurnerRZ::PutNumberofMedium(int i)
{
  mednum=i;
  med.resize(mednum);
};

void FRBurnerRZ::PutMediumIC(Medium &min)
{
  med[ic_s]=min;
};

void FRBurnerRZ::PutMediumOC(Medium &min)
{
  med[oc_s]=min;
};

void FRBurnerRZ::PutMediumRB(Medium &min)
{
  med[rb_s]=min;
};

void FRBurnerRZ::PutMediumAB(Medium &min)
{
  med[abic_s]=min;
};

void FRBurnerRZ::PutMedium(Medium &min,int mid)
{
  med[mid]=min;
};

void FRBurnerRZ::PutSADataFuel(FRDTSubAssemblyGeometry sa_geom)
{
  num_pin=sa_geom.GetNumberOfPin();
  area_assembly=sa_geom.GetAssemblyVolume();
  area_assembly*=0.01; // convert from mm^2 to cm^2
};

void FRBurnerRZ::PutSADataRadialBlanket(FRDTSubAssemblyGeometry sarb_geom)
{
  num_pin_rb=sarb_geom.GetNumberOfPin();
};

void FRBurnerRZ::PutCartMeshInfo(CartMeshInfo &cmii)
{
  cmi=cmii;
};

void FRBurnerRZ::PutMediumIDForIC(int is,int ie)
{
  ic_s=is;
  ic_e=ie;
  nic=ie-is+1;
};

void FRBurnerRZ::PutMediumIDForOC(int is,int ie)
{
  oc_s=is;
  oc_e=ie;
  noc=ie-is+1;
};

void FRBurnerRZ::PutMediumIDForRB(int is,int ie)
{
  rb_s=is;
  rb_e=ie;
  nrb=ie-is+1;
};

void FRBurnerRZ::PutMediumIDForABIC(int is,int ie)
{
  abic_s=is;
  abic_e=ie;
  nabic=ie-is+1;
};

void FRBurnerRZ::PutMediumIDForABOC(int is,int ie)
{
  aboc_s=is;
  aboc_e=ie;
  naboc=ie-is+1;
};

void FRBurnerRZ::PutBatchNumber(int i,int j,int k)
{
  batch_ic = i;
  batch_oc = j;
  batch_rb = k;

  batch_max=batch_ic;
  if(batch_oc>batch_max)batch_max=batch_oc;
  if(batch_rb>batch_max)batch_max=batch_rb;
};

int FRBurnerRZ::GetBatchNumberFromMediumID(int i)
{
  int tmp=batch_ic;
  if(i>=oc_s&&i<=oc_e)tmp=batch_oc;
  if(i>=rb_s&&i<=rb_e)tmp=batch_rb;
  if(i>=aboc_s&&i<=aboc_e)tmp=batch_oc;
  return tmp;
};

void FRBurnerRZ::PreCalculation(XSLibrary &xslib, Burnup &bu)
{
  med_burn=nic+noc+nrb+nabic+naboc;
 
  vol_per_med.resize(mednum);
  thm_per_med.resize(mednum);

  // +++ Initial selfshielding calculation
  for(int i=med_burn;i<mednum;i++){
    opc.CalSelfShieldingInfiniteSystem(med[i],xslib);
  };
  if(nic!=0){
    opc.CalSelfShieldingInfiniteSystem(med[ic_s],xslib);
    med[ic_s].CalMuFromSigel_p1();
  };
  if(noc!=0){
    opc.CalSelfShieldingInfiniteSystem(med[oc_s],xslib);
    med[oc_s].CalMuFromSigel_p1();
  };
  if(nrb!=0){
    opc.CalSelfShieldingInfiniteSystem(med[rb_s],xslib);
    med[rb_s].CalMuFromSigel_p1();
  };
  if(nabic!=0||naboc!=0){
    opc.CalSelfShieldingInfiniteSystem(med[abic_s],xslib);
    med[abic_s].CalMuFromSigel_p1();
  };

  for(int i=ic_s+1;i<=ic_e;i++){
    med[i]=med[ic_s];
  };
  for(int i=oc_s+1;i<=oc_e;i++){
    med[i]=med[oc_s];
  };
  for(int i=rb_s+1;i<=rb_e;i++){
    med[i]=med[rb_s];
  };
  for(int i=abic_s+1;i<=aboc_e;i++){
    med[i]=med[abic_s];
  };

  density.resize(med_burn);
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
  };

  sigc_1g.resize(med_burn);
  sigf_1g.resize(med_burn);
  sign2n_1g.resize(med_burn);
  for(int i=0;i<med_burn;i++){
    sigc_1g[i].resize(trace_cycle);
    sigf_1g[i].resize(trace_cycle);
    sign2n_1g[i].resize(trace_cycle);
    for(int j=0;j<trace_cycle;j++){
      sigc_1g[i][j].resize(burn_nuc,0.);
      sigf_1g[i][j].resize(burn_nuc,0.);
      sign2n_1g[i][j].resize(burn_nuc,0.);
    };
  };

  iden_ic.resize(burn_nuc);
  iden_oc.resize(burn_nuc);
  iden_rb.resize(burn_nuc);
  iden_ab.resize(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    if(nic!=0)iden_ic[i]=med[ic_s].GetNuclideInTurn(i).GetDensity();
    if(noc!=0)iden_oc[i]=med[oc_s].GetNuclideInTurn(i).GetDensity();
    if(nrb!=0)iden_rb[i]=med[rb_s].GetNuclideInTurn(i).GetDensity();
    if(nabic!=0||naboc!=0)iden_ab[i]=med[abic_s].GetNuclideInTurn(i).GetDensity();
  };
  InitializeNumberDensity();

  // Neutron flux calculation to store flux data in the Medium class
  // because this flux is used to calculate fission spectum vector
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
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real k1=test.CalIgen(optname);
  for(int i=0;i<mednum;i++){
    GroupData1D flx=test.GetIntegratedFlux(i);
    med[i].GetFlux().copy(flx);
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
  for(int i=ic_s;i<=ic_e;i++){
    for(int b=0;b<batch_ic;b++){
      for(int j=0;j<burn_nuc;j++){ 
        density[i][b][j]=iden_ic[j];
      };
    };
  };
  for(int i=oc_s;i<=oc_e;i++){
    for(int b=0;b<batch_oc;b++){
      for(int j=0;j<burn_nuc;j++){
        density[i][b][j]=iden_oc[j];
      };
    };
  };
  for(int i=rb_s;i<=rb_e;i++){
    for(int b=0;b<batch_rb;b++){
      for(int j=0;j<burn_nuc;j++){
        density[i][b][j]=iden_rb[j];
      };
    };
  };
  for(int i=abic_s;i<=abic_e;i++){
    for(int b=0;b<batch_ic;b++){
      for(int j=0;j<burn_nuc;j++){
        density[i][b][j]=iden_ab[j];
      };
    };
  };
  for(int i=aboc_s;i<=aboc_e;i++){
    for(int b=0;b<batch_oc;b++){
      for(int j=0;j<burn_nuc;j++){
        density[i][b][j]=iden_ab[j];
      };
    };
  };
};

void FRBurnerRZ::Run(XSLibrary &xslib, Burnup &bu)
{
  cout<<"################################################################\n";

  if(void_cal||dop_cal){
    cout<<"#  FRBurnerRZ : Calculation condition\n#\n";
    if(void_cal)cout<<"#       void ratio        : "<<void_ratio<<"\n";
    if(dop_cal) cout<<"#       delta temperature : "<<delta_t<<"\n";
    cout<<"################################################################\n";
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
  density_store.resize(trace_cycle);

  // +++
  GeneralOption opt,opta;
  opta.PutAdjointCal();

  int refuel_id_ic=0;
  int refuel_id_oc=0;
  int refuel_id_b=0;
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

  int med_fuel_max=nic+noc+nrb+nabic+naboc-1;

  cout<<"#C Step  Day  Keff   Max. line   C.R.  Power Dist. [%]      ";
  if(void_cal)cout<<"Void      ";
  if(dop_cal) cout<<"Doppler   ";
  cout<<"\n";
  cout<<"#y                   power[W/cm]        IC   OC   RB   AB   ";
  if(void_cal)cout<<"React.    ";
  if(dop_cal) cout<<"React.    ";
  cout<<"\n";
  cout<<"#                    (pos:r,z)                              ";
  if(void_cal)cout<<"[dk/kk']  ";
  if(dop_cal) cout<<"[dk/kk']  ";
  cout<<"\n";
  cout<<"#\n";

  // (cycle-loop)
  real ac_day=0.;
  for(int cycle=0;cycle<trace_cycle;cycle++){

    flux_level[cycle].resize(cycle_div+1);
    density_store[cycle].resize(cycle_div+1);

      // +++ Refueling
    if(cycle!=0){
      for(int i=ic_s;i<=ic_e;i++){
        ir_cycle[i][refuel_id_ic]=0;
        ac_burn[i][refuel_id_ic]=0.;
        for(int j=0;j<burn_nuc;j++){
          density[i][refuel_id_ic][j]=iden_ic[j];
        };
      };
      for(int i=oc_s;i<=oc_e;i++){
        ir_cycle[i][refuel_id_oc]=0;
        ac_burn[i][refuel_id_oc]=0.;
        for(int j=0;j<burn_nuc;j++){
	  density[i][refuel_id_oc][j]=iden_oc[j];
        };
      };
      for(int i=rb_s;i<=rb_e;i++){
        ir_cycle[i][refuel_id_b]=0;
        ac_burn[i][refuel_id_b]=0.;
        for(int j=0;j<burn_nuc;j++){
	  density[i][refuel_id_b][j]=iden_rb[j];
	};
      };
      for(int i=abic_s;i<=abic_e;i++){
        ir_cycle[i][refuel_id_ic]=0;
        ac_burn[i][refuel_id_ic]=0.;
        for(int j=0;j<burn_nuc;j++){
	  density[i][refuel_id_ic][j]=iden_ab[j];
	};
      };
      for(int i=aboc_s;i<=aboc_e;i++){
        ir_cycle[i][refuel_id_oc]=0;
        ac_burn[i][refuel_id_oc]=0.;
        for(int j=0;j<burn_nuc;j++){
	  density[i][refuel_id_oc][j]=iden_ab[j];
	};
      };
      refuel_id_ic++;
      refuel_id_oc++;
      refuel_id_b++;
      if(refuel_id_ic==batch_ic)refuel_id_ic=0;
      if(refuel_id_oc==batch_oc)refuel_id_oc=0;
      if(refuel_id_b==batch_rb)refuel_id_b=0;
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
      for(int i=0;i<med_burn;i++){
        med[i].CalMacroFromMicro();
        //med[i].CalMacroFromMicroSimple();
      };

      // +++ Eigenvalue calculation ++++++++++++++++++++++++++++++++++++++
      PLOSSystem test(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
	if(dop_cal)opc.CalSelfShieldingInfiniteSystem(med[i],xslib);
        test.AddMedium(med[i]);
	test.GetMedium(i).MicxsVectorClear2DData();
      };
      test.NoPrint();
      test.PutCartMeshInfo(cmi,"Cylinder");
      test.PutGeneralOption(opt);
      test.CalCoef();
      string optname="cmfd";
      if(!cmfd_on)optname="";
      real k1=test.CalIgen(optname);

      if(show_fission_info)test.GetNeutronMultiplicationInfo(k1);

      // ++++ Void condition +++++++++++++++++++++++++++++++++++++++++
      real kv=0.;
      if(void_cal){
        PLOSSystem test_v(dim,group,mednum); 
        for(int i=0;i<mednum;i++){
          test_v.AddMedium(med[i]);
	  //if(i<med_burn&&test_v.GetMedium(i).ExistNuclide(110230)){
	  if(test_v.GetMedium(i).ExistNuclide(110230)){
	    real org=test_v.GetMedium(i).GetNuclide(110230).GetDensity();
	    test_v.GetMedium(i).GetNuclide(110230).PutDensity(org*void_ratio);
	    test_v.GetMedium(i).CalMacroFromMicro();
	    //test_v.GetMedium(i).CalMacroFromMicroSimple();
	  };
  	  test_v.GetMedium(i).MicxsVectorClear2DData();
        };
        test_v.NoPrint();
        test_v.PutCartMeshInfo(cmi,"Cylinder");
        test_v.PutGeneralOption(opt);
        test_v.CalCoef();
        kv=test_v.CalIgen(optname);
        //kv=test_v.CalIgen();
        //test_v.CalReactivity(&test,kv,k1);
      };

      // ++++ Doppler condition +++++++++++++++++++++++++++++++++++++++++
      real kd=0.;
      if(dop_cal){
        PLOSSystem test_v(dim,group,mednum); 
        for(int i=0;i<mednum;i++){
          test_v.AddMedium(med[i]);
	  for(int j=0;j<med[i].GetNucnum();j++){
            int id=med[i].GetNuclideInTurn(j).GetMatnum();
	    if(id>900000){
    	      real org=test_v.GetMedium(i).GetNuclideInTurn(j).GetTemperature();
  	      test_v.GetMedium(i).GetNuclideInTurn(j).PutTemperature(org+delta_t);
	    };
	  };
          opc.CalSelfShieldingInfiniteSystem(test_v.GetMedium(i),xslib);
          test_v.GetMedium(i).CalMacroFromMicro();
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

      vector< vector<real> > pow_store(med_burn);
      // Volume-integrated power for each medium/batch
      // Note : "Volume" is defined as "total" one for all the batches

      real pow_sum=0.;
      for(int i=0;i<mednum;i++){
        GroupData1D flx=test.GetIntegratedFlux(i);
        med[i].GetFlux().copy(flx);
	if(i<med_burn){
	  int sz=ac_burn[i].size(); // batch-wise
	  real pow=0.;
          pow_store[i].resize(sz);
	  for(int j=0;j<sz;j++){ 
	    for(int k=0;k<burn_nuc;k++){
	      med[i].GetNuclideInTurn(k).PutDensity(density[i][j][k]);
	    };
            real tmp=bu.GetIntegratedPower(med[i]);
            pow_store[i][j]=tmp;
	    pow+=tmp;
	  };
	  pow_sum+=pow/sz;
	}else{
	  pow_sum+=bu.GetIntegratedPower(med[i]);
	};
      };

      real factor=power/pow_sum;

      /*
      if(print_power_density){
	cout<<"\n+++ Power density ratio to core average value +++\n";
	cout<<"\n";
        real volsum=0.;
        for(int m=0;m<mednum;m++){
          volsum+=test.GetVolumeParMedium(m);
        };
        real pow_den_avg=pow_sum/volsum;
	for(int m=0;m<mednum;m++){
	  real vol=test.GetVolumeParMedium(m);
          real poww=bu.GetIntegratedPowerParMedium(test,m);
	  cout<<m<<" "<<(poww/vol)/pow_den_avg<<"\n";
	};
	cout<<"\n";
      };
      */

      // +++ Maximum line power
      real maxpow=0.;
      int index=0;
      int maxp_x=-1;
      int maxp_y=-1;

      int xmax=0;
      int ymin=1000;
      int ymax=0;
      for(int y=0;y<yf;y++){
        for(int x=0;x<xf;x++){
  	  power_map[y][x]=0.;
          int medid=cmi.GetFMat(index);
          if(medid<=med_fuel_max){
	    real pow=bu.GetIntegratedPowerParMesh(test,x,y,0);
	    real vol=test.GetMesh(x,y,0).GetVolume();
	    real zl=test.GetMesh(x,y,0).GetLen(1);
	    real area=vol/zl;
            int pinn=num_pin;
            if(medid>=rb_s&&medid<=rb_e)pinn=num_pin_rb;
	    pow=(pow/zl)/(pinn*area/area_assembly);
	    if(pow>maxpow){
              maxpow=pow;
              maxp_x=x;
              maxp_y=y;
	    };
	    power_map[y][x]=pow*factor;
	    if(x>xmax)xmax=x;
	    if(y>ymax)ymax=y;
	    if(y<ymin)ymin=y;
	  };
          index++;
	};
      };

      if(print_linepower_map){
        cout<<"\n# Line power map in fine mesh [W/cm]\n";
        cout<<"#    cycle: "<<cycle<<" / step: "<<tim<<"\n#\n";
        for(int y=ymin;y<=ymax;y++){
          for(int x=0;x<=xmax;x++){
	    int tmp=power_map[y][x];
	    WriteOut(tmp,3);
	    cout<<" ";
	  };
	  cout<<"\n";
	};
	cout<<"\n\n";
	fout_pmap.setf(ios::scientific);
	fout_pmap.precision(5);
        for(int x=0;x<=xmax;x++){
	  for(int jj=0;jj<2;jj++){
  	    for(int y=ymin;y<=ymax;y++){
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
      WriteOut(maxpow*factor,"%5.1f");
      cout<<" ("<<maxp_x<<","<<maxp_y<<")";

      // ++ Conversion ratio printing 
      real cr=test.CalConversionRatio();
      cout<<" ";
      WriteOut(cr,"%5.3f");

      // ++ Power distribution printing
      int stt[]={ic_s,oc_s,rb_s,abic_s};
      int edd[]={ic_e,oc_e,rb_e,aboc_e};
      for(int ii=0;ii<4;ii++){
        real pow_reg=0.;
        for(int i=stt[ii];i<=edd[ii];i++){
          int sz=ac_burn[i].size();
	  for(int j=0;j<sz;j++){
            pow_reg+=pow_store[i][j]/sz;
	  };
	};
	pow_reg/=pow_sum;
        pow_reg*=100;
	WriteOut(pow_reg,"%5.1f");
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
          GroupData1D flx=test.GetIntegratedFluxParVolume(m);
	  real flxsum_inv=1./flx.get_sum();
          real sumflx=flx.get_sum()*factor;
	  flux_level[cycle][tim][m]=sumflx;
          for(int i=0;i<burn_nuc;i++){
            sigf_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sigf)*flx*flxsum_inv;
            sigc_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sigc)*flx*flxsum_inv;
            sign2n_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n)*flx*flxsum_inv;
	    //if(m==0)cout<<med[0].GetNuclideInTurn(i).GetMatnum()<<" "<<fis[i]<<" "<<cap[i]<<" "<<n2n[i]<<"\n";
          };
	  //exit(0);

          for(int i=0;i<burn_nuc;i++){
            int id=med[m].GetNuclideInTurn(i).GetMatnum();
            //bu.PutNuclideData(i,id,0.,fis[i],cap[i],n2n[i]);
            bu.PutNuclideData(i,id,0.,sigf_1g[m][cycle][i],sigc_1g[m][cycle][i],sign2n_1g[m][cycle][i]);
          };
          int tmp=GetBatchNumberFromMediumID(m);
          for(int j=0;j<tmp;j++){
            for(int i=0;i<burn_nuc;i++){
              bu.PutDensity(i,density[m][j][i]);
  	    };
            bu.CalTransitionMatrixFluxDependentPart();
            bu.BurnupCalculationByKrylov(sumflx,delt);
            for(int i=0;i<burn_nuc;i++){
              density[m][j][i]=bu.GetDensity(i);
            };
          };
	};
	for(int m=med_burn;m<mednum;m++){
          GroupData1D flx=test.GetIntegratedFluxParVolume(m);
          real sumflx=flx.get_sum()*factor;
	  flux_level[cycle][tim][m]=sumflx;
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
	      //bu.BurnupCalculation(0.,decay_cal_day*24*60*60,"pade",false);
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

void FRBurnerRZ::RunADS(XSLibrary &xslib, Burnup &bu, bool print, bool denout)
{
  real avo=0.60221367;

  Medium newfuel; 
  // !! Ogawa-kun !!
  // This medium is to store number densities of initially-loaded MA nuclides (Np, Am, Cm).
  // It is assumed that these MA nuclides are taken from waste fuel of LWR.
  vector<int> matid(burn_nuc);
  vector<real> matden(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    matid[i]=med[0].GetNuclideInTurn(i).GetMatnum();
    matden[i]=0.;
    if((matid[i]>930000&&matid[i]<940000)||
       (matid[i]>950000&&matid[i]<970000)){
      // !! Ogawa-kun !!
      // Number densities of only Np, Am, Cm are stored.
      matden[i]=med[0].GetNuclideInTurn(i).GetDensity();
    };
  };
  newfuel.PutNuclide(25,matid,matden);
  real hm_w_newfuel=bu.CalWeightOfHeavyNuclideParUnitVolume(newfuel);
  // !! Ogawa-kun !!
  // Heavy metal weight of initially-loaded MAs taken from LWR

  real hm_w_icfuel=bu.CalWeightOfHeavyNuclideParUnitVolume(med[ic_s]);
  //real hm_w_ocfuel=bu.CalWeightOfHeavyNuclideParUnitVolume(med[oc_s]);
  // !! Ogawa-kun !!
  // Heavy metal weight of initially-loaded fuel

  // +++ Fixed-source calculation option
  real epsf=1e-4;
  int itermax=1000;
  bool high_speed_option=true;

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // External neutron source setting

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
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // (initialize number of irradiated cycles and accumulated burnup)
  for(int i=0;i<med_burn;i++){
    int sz=ir_cycle[i].size();
    for(int b=0;b<sz;b++){
      ir_cycle[i][b]=0;
      ac_burn[i][b]=0.;
    };
  };

  // +++
  density_store.resize(trace_cycle);

  // +++
  GeneralOption opt;

  int refuel_id_ic=0;
  int refuel_id_oc=0;
  int refuel_id_b=0;
  bu.PutNucnum(burn_nuc);
  for(int i=0;i<burn_nuc;i++){
    bu.PutNuclideData(i,med[0].GetNuclideInTurn(i).GetMatnum(),0.,0.,0.,0.);
  };
  bu.CalTransitionMatrixFluxInDependentPart();

  flux_level.resize(trace_cycle);

  int med_fuel_max=nic+noc+nrb+nabic+naboc-1;

  cout<<"#C Step EFPD  Keff   Max. line   Power Dist. [%]     \n";
  cout<<"#y                   power[W/cm]  IC   OC   RB   AB  \n";
  cout<<"#                    (pos:r,z)                       \n";
  cout<<"#\n";

  // +++ for ogawa calculation
  int p8_pos=bu.SearchNuclide(942380);
  int p0_pos=bu.SearchNuclide(942400);
  if(p8_pos==-1){
    cout<<"# Error in FRBurnerRZ::RunADS.\n";
    cout<<"# Pu-238 cannot be found.\n";
    exit(0);
  };
  if(p0_pos==-1){
    cout<<"# Error in FRBurnerRZ::RunADS.\n";
    cout<<"# Pu-240 cannot be found.\n";
    exit(0);
  };
  vector<real> den_p8_store(med_burn);
  vector<real> den_p0_store(med_burn);
  // +++++++++++++++++++++++++++++++;



  // (cycle-loop)
  real ac_day=0.;
  for(int cycle=0;cycle<trace_cycle;cycle++){

    flux_level[cycle].resize(cycle_div+1);
    density_store[cycle].resize(cycle_div+1);

    // +++ Refueling
    if(cycle!=0){

      vector<real> density_tmp(burn_nuc);
      real density_tmp_p8=0.;
      real density_tmp_p0=0.;
      // !! Ogawa-kun !!
      // MA Nuclide number densities of whole core after this burnup cycle are stored.
      real volsum=0.;
      for(int i=ic_s;i<=oc_e;i++){ // medium-loop
	real vol=vol_per_med[i];
        for(int j=0;j<burn_nuc;j++){ // nuclide-loop
	  int id=med[0].GetNuclideInTurn(j).GetMatnum();
	  if(id>920000&&id<1000000){
	    if(id==942380&&ogawa_cal){
  	      density_tmp[j]+=(density[i][0][j]-den_p8_store[i])*vol;
  	      density_tmp_p8+=den_p8_store[i]*vol;
	    }else if(id==942400&&ogawa_cal){
  	      density_tmp[j]+=(density[i][0][j]-den_p0_store[i])*vol;
  	      density_tmp_p0+=den_p0_store[i]*vol;
	    }else{
	      // !! Ogawa-kun !!
              // Only actinoids [U~Es] are considered
  	      density_tmp[j]+=density[i][0][j]*vol;
	      // !! Ogawa-kun !!
              // Volume-integrated number density is calculated for each nuclide. [#]
	    };
	  };
	};
	volsum+=vol;
      };

      if(print){
      cout<<"#\n";
      cout<<"# Total heavy nuclides inventory and decay heat at EOC [per cm^3]\n";
      cout<<"#\n";
      cout.setf(ios::scientific);
      cout.precision(3);
      real ensum=0.;
      for(int i=0;i<burn_nuc;i++){
	if(density_tmp[i]>0.){
	int mid=med[0].GetNuclideInTurn(i).GetMatnum();
	cout<<"#   ";
	WriteOut(mid,7);
	cout<<" "<<density_tmp[i]/volsum<<" "; // number density [#/cm3]
	real dc=bu.GetBurnupChain().GetDecayConstant(mid);
	real de=bu.GetBurnupChain().GetDecayEnergy(mid,0);
	de+=bu.GetBurnupChain().GetDecayEnergy(mid,1);
	de+=bu.GetBurnupChain().GetDecayEnergy(mid,2);
	real tmp1=density_tmp[i]/volsum*dc*de; // decay heat [W/cm3]
	ensum+=tmp1;
	cout<<tmp1<<"\n";
	};
      };
      cout<<"#   (Total)           "<<ensum<<"\n"; // Total decay heat
      cout<<"#\n#\n";
      };

      Medium med_dummy; 
      // !! Ogawa-kun !!
      // This instance is to store core-averagd number densities for actinoids
      for(int i=0;i<burn_nuc;i++){
        Nuclide dummy;
	dummy.PutMatnum(med[0].GetNuclideInTurn(i).GetMatnum());
	dummy.PutDensity(density_tmp[i]/volsum);
	med_dummy.AddNuclide(dummy);
      };
      real hm_w_oldfuel=bu.CalWeightOfHeavyNuclideParUnitVolume(med_dummy);
      // !! Ogawa-kun !!
      // Heavy metal weight of the core-averaged number densities for actinoids
      real factor=(hm_w_icfuel-hm_w_oldfuel)/hm_w_newfuel;
      // !! Ogawa-kun !!
      // Heavy metal weight lost during burnup

      vector<real> added_den(burn_nuc,0.);
      for(int i=ic_s;i<=oc_e;i++){
	real vol=vol_per_med[i];
        for(int j=0;j<burn_nuc;j++){
	  density[i][0][j]=density_tmp[j]/volsum+matden[j]*factor;
	  // !! Ogawa-kun !! 
	  // MA from LWR is added to preserve heavy metal weights
	  added_den[j]+=matden[j]*factor*vol; 
        };
      };

      if(print){
      cout<<"#\n# Total weight of newly-loaded MA materials [kg]\n#\n";
      cout<<"#    (Cycle : "<<cycle<<")\n#\n";
      real sum=0.;
      for(int i=0;i<burn_nuc;i++){
	if(added_den[i]>0.){
	int mid=med[0].GetNuclideInTurn(i).GetMatnum();
        real aw=bu.GetAtomicWeight(mid);
        real wgt=added_den[i]/avo*aw*1e-3; // [kg]
	sum+=wgt;
	cout<<"# ";
	WriteOut(mid,7);
	cout<<" ";
	WriteOut(wgt,"%8.2f");
	cout<<"\n";
	};
      };
      cout<<"#  Total  ";
      WriteOut(sum,"%8.2f");
      cout<<"\n#\n";
    };
      if(ogawa_cal){
	cout<<"#\n# Total weight of extracted Pu-238 and -240 [kg]\n#\n";
        cout<<"#    (Cycle : "<<cycle<<")\n#\n";
        real aw_p8=bu.GetAtomicWeight(942380);
        real wgt_p8=(density_tmp_p8/avo)*aw_p8*1e-3; // [kg]
        real aw_p0=bu.GetAtomicWeight(942400);
        real wgt_p0=(density_tmp_p0/avo)*aw_p0*1e-3; // [kg]
	cout<<"#      Pu-238 : ";
        WriteOut(wgt_p8,"%8.2f");
	cout<<"\n";
	cout<<"#      Pu-240 : ";
        WriteOut(wgt_p0,"%8.2f");
	cout<<"\n#\n";
      };
    };

    // (in cycle-loop)
    for(int tim=0;tim<cycle_div+1;tim++){

      flux_level[cycle][tim].resize(med_burn);
      real delt_d=cycle_length/cycle_div; // [d]
      real delt=delt_d*24*60*60; // [s]

      // +++ Calculation of region-averaged density
      CalRegionAveragedDensity();
      StoreRegionAveragedDensity(cycle,tim);

      if(denout){
      if(cycle!=0&&tim==0){
        med[0].WriteFileNumberDensity("./DEN_BOC/","den_cyc"+IntToString(cycle));
      };
      if(tim==cycle_div){
	for(int i=0;i<med_burn;i++){
	  med[i].WriteFileNumberDensity("./DEN_EOC/","den_cyc"+IntToString(cycle)+"_med"+IntToString(i));
	};
      };
      };

      // +++ Macroscopic cross section calculation
      for(int i=0;i<med_burn;i++){
        med[i].CalMacroFromMicro();
      };

      // +++ Eigenvalue calculation ++++++++++++++++++++++++++++++++++++++
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
      if(!cmfd_on)optname="";
      real k1=test.CalIgen(optname);

      // +++ Fixed-source calculation ++++++++++++++++++++++++++++++++++++
      PLOSSystem test2(dim,group,mednum); 
      for(int i=0;i<mednum;i++){
        test2.AddMedium(med[i]);
	test2.GetMedium(i).MicxsVectorClear2DData();
      };
      test2.NoPrint();
      test2.PutCartMeshInfo(cmi,"Cylinder");
      test2.PutGeneralOption(opt);
      test2.CalCoef();
      test2.SetZeroScatSrc();
      for(int iz=0;iz<54;iz++){
	for(int ir=0;ir<38;ir++){
	  test2.PutIsotropicSourceParVolume(ir,ir,iz,iz,0,0,esrc[iz][ir]);
	};
      };
      test2.CalFixedSourceWithFission(epsf,itermax,high_speed_option);

      if(show_fission_info)test2.GetNeutronMultiplicationInfo(k1);

      if(cycle==0&&tim==0){
	for(int i=0;i<med_burn;i++){
	  vol_per_med[i]=test.GetVolumePerMedium(i);
          real weight=bu.CalWeightOfHeavyNuclideParUnitVolume(med[i]);
          thm_per_med[i]=weight*vol_per_med[i];
	};
      };

      vector< vector<real> > pow_store(med_burn);
      // Volume-integrated power for each medium/batch
      // Note : "Volume" is defined as "total" one for all the batches

      real pow_sum=0.;
      for(int i=0;i<mednum;i++){
        GroupData1D flx=test2.GetIntegratedFlux(i);
        med[i].GetFlux().copy(flx);
	if(i<med_burn){
	  int sz=ac_burn[i].size(); // batch-wise
	  real pow=0.;
          pow_store[i].resize(sz);
	  for(int j=0;j<sz;j++){ 
	    for(int k=0;k<burn_nuc;k++){
	      med[i].GetNuclideInTurn(k).PutDensity(density[i][j][k]);
	    };
            real tmp=bu.GetIntegratedPower(med[i]);
            pow_store[i][j]=tmp;
	    pow+=tmp;
	  };
	  pow_sum+=pow/sz;
	}else{
	  pow_sum+=bu.GetIntegratedPower(med[i]);
	};
      };

      real factor=power/pow_sum;

      // +++ Maximum line power
      if(print_linepower_map){
        cout<<"\n# Line power map in fine mesh [W/cm]\n";
        cout<<"#    cycle: "<<cycle<<" / step: "<<tim<<"\n#\n";
      };
      real maxpow=0.;
      int index=0;
      int maxp_x=-1;
      int maxp_y=-1;
      for(int y=0;y<cmi.GetYF();y++){
        for(int x=0;x<cmi.GetXF();x++){
          int medid=cmi.GetFMat(index);
          if(medid<=med_fuel_max){
	    real pow=bu.GetIntegratedPowerParMesh(test2,x,y,0);
	    real vol=test2.GetMesh(x,y,0).GetVolume();
	    real zl=test2.GetMesh(x,y,0).GetLen(1);
	    real area=vol/zl;
            int pinn=num_pin;
            if(medid>=rb_s&&medid<=rb_e)pinn=num_pin_rb;
	    pow=(pow/zl)/(pinn*area/area_assembly);
	    if(pow>maxpow){
              maxpow=pow;
              maxp_x=x;
              maxp_y=y;
	    };
	    if(print_linepower_map){
              int tmp=pow*factor;
              WriteOut(tmp,3);
  	      cout<<" ";
	    };
	  }else{
            if(print_linepower_map)cout<<"    ";
	  };
          index++;
	};
	if(print_linepower_map)cout<<"\n";
      };
      if(print_linepower_map)cout<<"\n\n";
      //cout<<"   Max. line power : "<<maxpow*factor<<" (W/cm)\n";

      WriteOut(cycle,2);
      WriteOut(tim,4);
      int iac_day=int(ac_day);
      WriteOut(iac_day,6);
      cout<<" ";
      WriteOut(k1,"%7.5f");
      cout<<" ";
      WriteOut(maxpow*factor,"%5.1f");
      cout<<" ("<<maxp_x<<","<<maxp_y<<")";

      // ++ Power distribution printing
      int stt[]={ic_s,oc_s,rb_s,abic_s};
      int edd[]={ic_e,oc_e,rb_e,aboc_e};
      for(int ii=0;ii<4;ii++){
        real pow_reg=0.;
        for(int i=stt[ii];i<=edd[ii];i++){
          int sz=ac_burn[i].size();
	  for(int j=0;j<sz;j++){
            pow_reg+=pow_store[i][j]/sz;
	  };
	};
	pow_reg/=pow_sum;
        pow_reg*=100;
	WriteOut(pow_reg,"%5.1f");
      };

      cout<<"\n";

      // Burnup calculation
      if(tim!=cycle_div){

	for(int i=0;i<med_burn;i++){
	  int sz=ac_burn[i].size();
	  for(int j=0;j<sz;j++){
	    ac_burn[i][j]+=((pow_store[i][j]*factor*1e-9)*delt_d)/(thm_per_med[i]*1e-6); // [GWd/t]
	  };
	};

        ac_day+=delt_d;
        for(int m=0;m<med_burn;m++){ // fuel medium
          GroupData1D flx=test2.GetIntegratedFluxParVolume(m);
	  real flxsum_inv=1./flx.get_sum();
          real sumflx=flx.get_sum()*factor;
	  flux_level[cycle][tim][m]=sumflx;
          for(int i=0;i<burn_nuc;i++){
            sigf_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sigf)*flx*flxsum_inv;
            sigc_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sigc)*flx*flxsum_inv;
            sign2n_1g[m][cycle][i]=med[m].GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n)*flx*flxsum_inv;
          };

          for(int i=0;i<burn_nuc;i++){
            int id=med[m].GetNuclideInTurn(i).GetMatnum();
            bu.PutNuclideData(i,id,0.,sigf_1g[m][cycle][i],sigc_1g[m][cycle][i],sign2n_1g[m][cycle][i]);
          };
          int tmp=GetBatchNumberFromMediumID(m);
          for(int j=0;j<tmp;j++){
            for(int i=0;i<burn_nuc;i++){
              bu.PutDensity(i,density[m][j][i]);
  	    };
            bu.CalTransitionMatrixFluxDependentPart();
            bu.BurnupCalculationByKrylov(sumflx,delt);
            for(int i=0;i<burn_nuc;i++){
              density[m][j][i]=bu.GetDensity(i);
            };
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
	  if(!ogawa_cal){
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
	  }else{
	    // (ogawa-calculation)
  	    for(int m=0;m<med_burn;m++){
	      int sz=ac_burn[m].size();
	      if(sz!=1){
		cout<<"# Error in FRBurnerRZ::RunADS.\n";
		cout<<"# Please see the source.\n";
		exit(0);
	      };
	      for(int i=0;i<sz;i++){
	        for(int j=0;j<burn_nuc;j++){
	          bu.PutDensity(j,density[m][i][j]);
	        };
		// The first step
                bu.BurnupCalculation(1e-15,ogawa_year[0]*365*24*60*60,"chebyshev",false);
		// The second step
		// (Number densities of Cm isotopes are set to zero)
		for(int j=0;j<burn_nuc;j++){
                  int idtmp=bu.GetNuclideID(j);
		  if(idtmp>=960000&&idtmp<970000){
                    density[m][i][j]=bu.GetDensity(j);
                    bu.PutDensity(j,0.);
		  };
		};
                bu.BurnupCalculation(1e-15,ogawa_year[1]*365*24*60*60,"chebyshev",false);
	        for(int j=0;j<burn_nuc;j++){
                  int idtmp=bu.GetNuclideID(j);
		  if(idtmp>=960000&&idtmp<970000){
                    real tmp=bu.GetDensity(j);
		    bu.PutDensity(j,density[m][i][j]);
                    density[m][i][j]=tmp;
		  }else{
                    real tmp=bu.GetDensity(j);
                    density[m][i][j]=tmp;
		    bu.PutDensity(j,0.);
		  };
	        };
                bu.BurnupCalculation(1e-15,ogawa_year[1]*365*24*60*60,"chebyshev",false);
		den_p8_store[m]=bu.GetDensity(p8_pos);
		den_p0_store[m]=bu.GetDensity(p0_pos);
	        for(int j=0;j<burn_nuc;j++){
                  real tmp=bu.GetDensity(j);
                  bu.PutDensity(j,tmp+density[m][i][j]);
		};
		//cout<<"   "<<m<<" "<<den_p8_store[m]<<" "<<den_p0_store[m]<<"\n";
		//cout<<"   "<<m<<" "<<den0_pu238<<" "<<den1_pu238<<" "<<den0_pu240<<" "<<den1_pu240<<"\n";
                real tmp_sec=ogawa_year[2]*365.*24*60*60;
                bu.BurnupCalculation(1e-15,tmp_sec,"chebyshev",false);
                den_p8_store[m]*=exp(-bu.GetBurnupChain().GetDecayConstant(942380)*tmp_sec);
                den_p0_store[m]*=exp(-bu.GetBurnupChain().GetDecayConstant(942400)*tmp_sec);
		//cout<<"   "<<m<<" "<<den_p8_store[m]<<" "<<den_p0_store[m]<<"\n";
	        for(int j=0;j<burn_nuc;j++){
	          density[m][i][j]=bu.GetDensity(j);
	        };
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
    };
  };
};

void FRBurnerRZ::StoreRegionAveragedDensity(int cyc,int step)
{
  density_store[cyc][step].resize(med_burn);
  for(int i=0;i<med_burn;i++){
    density_store[cyc][step][i].resize(burn_nuc);
    for(int j=0;j<burn_nuc;j++){
      real den=med[i].GetNuclideInTurn(j).GetDensity();
      density_store[cyc][step][i][j]=den;
    };
  };
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
      real den=density_store[cyc][step][i][j];
      med[i].GetNuclideInTurn(j).PutDensity(den);
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
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };

  GeneralOption opt;
  if(adjoint)opt.PutAdjointCal();

  PLOSSystem test(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test.AddMedium(med[i]);
    test.GetMedium(i).MicxsVectorClear2DData();
  };
  //test.NoPrint();
  test.PutCartMeshInfo(cmi,"Cylinder");
  test.PutGeneralOption(opt);
  test.CalCoef();
  string optname="cmfd";
  if(!cmfd_on||adjoint)optname="";
  real k1=test.CalIgen(optname);

  cout.setf(ios::scientific);
  cout.precision(4);
  //for(int i=0;i<mednum;i++){
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
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };

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
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };

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

void FRBurnerRZ::CalVoidReactivity(real void_ratio, int meds, int mede)
{
  // ! NOTE !
  //
  // Microscopic cross section is invaried after coolant voiding,
  // so yield term should be zero.

  if(mede==-1)mede=meds;
  if(meds==-1){
    meds=0;
    mede=mednum;
  };

  // Macroscopic cross section calculations
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };

  GeneralOption opt,opta;
  opta.PutAdjointCal();

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
    test_v.AddMedium(med[i]);
    if(test_v.GetMedium(i).ExistNuclide(110230)&&i>=meds&&i<=mede){
      real org=test_v.GetMedium(i).GetNuclide(110230).GetDensity();
      test_v.GetMedium(i).GetNuclide(110230).PutDensity(org*void_ratio);
      test_v.GetMedium(i).CalMacroFromMicro();
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

  testa.CalReactivity(&test_v,k1,kv);
};


void FRBurnerRZ::CalDopplerReactivity(XSLibrary &xslib,real dt, int meds, int mede)
{
  // ! NOTE !

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

  GeneralOption opt,opta;
  opta.PutAdjointCal();

  OnePointCalculator opc;

  PLOSSystem testa(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    testa.AddMedium(med[i]);
    opc.CalSelfShieldingInfiniteSystem(testa.GetMedium(i),xslib);
    testa.GetMedium(i).CalMacroFromMicro();
    testa.GetMedium(i).MicxsVectorClear2DData();
  };
  testa.NoPrint();
  testa.PutCartMeshInfo(cmi,"Cylinder");
  testa.PutGeneralOption(opta);
  testa.CalCoef();
  real k1=testa.CalIgen();

  PLOSSystem test_v(dim,group,mednum); 
  for(int i=0;i<mednum;i++){
    test_v.AddMedium(med[i]);
    if(i>=meds&&i<=mede){
      for(int j=0;j<med[i].GetNucnum();j++){
        int id=med[i].GetNuclideInTurn(j).GetMatnum();
        if(id>900000){
          real org=test_v.GetMedium(i).GetNuclideInTurn(j).GetTemperature();
          test_v.GetMedium(i).GetNuclideInTurn(j).PutTemperature(org+dt);
        };
      };
    };
    opc.CalSelfShieldingInfiniteSystem(test_v.GetMedium(i),xslib);
    test_v.GetMedium(i).CalMacroFromMicro();
    test_v.GetMedium(i).MicxsVectorClear2DData();
  };
  test_v.NoPrint();
  test_v.PutCartMeshInfo(cmi,"Cylinder");
  test_v.PutGeneralOption(opt);
  test_v.CalCoef();
  string optname="cmfd";
  if(!cmfd_on)optname="";
  real kv=test_v.CalIgen(optname);

  testa.CalReactivity(&test_v,k1,kv);
};


void FRBurnerRZ::CalDelayedNeutronParameters(DelayedNeutronData &dnd)
{
  for(int i=0;i<med_burn;i++){
    med[i].CalMacroFromMicro();
  };

  GeneralOption opt,opta;
  opta.PutAdjointCal();

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

  testa.CalBetaEffective(&test,dnd);
  testa.CalNeutronLifeTime(&test);
};

void FRBurnerRZ::WriteFileMediumData(string dirname,string filename)
{
  for(int i=0;i<mednum;i++){
    if(i<med_burn)med[i].CalMacroFromMicro();
    med[i].WriteFile(dirname,filename+IntToString(i),true);
  };
};


// +++ printing

void FRBurnerRZ::PrintNuclideWeightPerBatch(Burnup &bu,bool init)
{
  real avo=0.60221367;

  real factor=1.;
  if(cmi.GetBC(2)==1)factor=2.; // when axially-symmetric core is treated

  vector< vector<real> > den_sum((4),vector<real>(burn_nuc,0.));
  int ii=0;
  for(int ireg=0;ireg<med_burn;ireg++){
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
      real tmp=density_store[0][0][ireg][i];
      if(!init)tmp=density[ireg][cyid][i];
      den_sum[ii][i]+=tmp*vol/batch*factor;
    };
    if(ireg==ic_e||ireg==oc_e||ireg==rb_e)ii++;
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
    cout<<"# so the following quantities are multiplied by 2.\n";
    cout<<"#\n";
  };
  cout<<"##########################################################\n";
  cout<<"#(Nucl)     (IC)     (OC)     (RB)     (AB)    (Sum)\n";

  vector<real> wgtsum(4,0.);
  for(int i=0;i<burn_nuc;i++){
    int mid=med[0].GetNuclideInTurn(i).GetMatnum();
    real aw=bu.GetAtomicWeight(mid);
    WriteOut(mid,7);
    cout<<" ";
    real tmp=0.;
    for(int j=0;j<4;j++){
      real wgt=den_sum[j][i]/avo*aw*1e-3; // [kg]
      WriteOut(wgt,"%8.2f");
      cout<<" ";
      wgtsum[j]+=wgt;
      tmp+=wgt;
    };
    WriteOut(tmp,"%8.2f");
    cout<<"\n";
  };

  cout<<" (Sum)  ";
  real tmp=0.;
  for(int j=0;j<4;j++){
    WriteOut(wgtsum[j],"%8.2f");
    tmp+=wgtsum[j];
    cout<<" ";
  };
  WriteOut(tmp,"%8.2f");
  cout<<"\n";
  cout<<"##########################################################\n";
};

void FRBurnerRZ::PrintND(Burnup &bu,string opt)
{
  real avo=0.60221367;
  real ev_to_j=1.60219e-19;

  vector< vector<real> > density_sum(batch_max);

  for(int i=0;i<batch_max;i++){
    density_sum[i].resize(burn_nuc,0.);
  };
  real vol_sum=0.;

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
     real vol=vol_per_med[ireg];
     vol_sum+=vol;
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
         //sum+=tmp;
         density_sum[j][i]+=val*vol;
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

     if(ireg==ic_e||ireg==oc_e||ireg==rb_e||ireg==abic_e||ireg==aboc_e){
       cout<<"#\n# +++ Macro region-averaged";
       if(ireg==ic_e)cout<<"   (inner core)";
       if(ireg==oc_e)cout<<"   (outer core)";
       if(ireg==rb_e)cout<<"   (radial blanket)";
       if(ireg==abic_e)cout<<"   (axial blanket-ic)";
       if(ireg==aboc_e)cout<<"   (axial blanket-oc)";
       cout<<" ( volume : "<<vol_sum<<" cm3 )\n#\n";
       for(int i=0;i<burn_nuc;i++){
         int mid=med[0].GetNuclideInTurn(i).GetMatnum();
         WriteOut(midt.Name(mid),8);
         cout<<"  ";
         //cout<<midt.Name(mid)<<"  ";
         real sum=0.;
         for(int j=0;j<sz;j++){
           cout<<density_sum[j][i]/vol_sum<<" ";
           sum+=density_sum[j][i]/vol_sum;
         };
         sum/=sz;
         if(sz>1)cout<<"   "<<sum;
         cout<<"\n";
       };
    
       for(int i=0;i<batch_max;i++){
         for(int j=0;j<burn_nuc;j++){               
            density_sum[i][j]=0.;
          };
       };
       vol_sum=0.;
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
  cout<<"#\n# Medium-dependent burnup\n#\n";
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
};

void FRBurnerRZ::PrintMacroXS()
{
  string nam[]={"Inner core","Outer core","Radial blanket","Axial blanket"};
  cout<<"###################################\n";
  cout<<"# Macroscopic cross section table #\n";
  cout<<"###################################\n";

  for(int i=0;i<4;i++){
    cout<<"\n\n# Medium : "<<nam[i]<<"\n\n";
    if(i==0)med[ic_s].PrintMacroSectionTable();
    if(i==1)med[oc_s].PrintMacroSectionTable();
    if(i==2)med[rb_s].PrintMacroSectionTable();
    if(i==3)med[abic_s].PrintMacroSectionTable();
  };
};

void FRBurnerRZ::PrintInitialND()
{
  cout<<"# Initial number density\n";
  cout<<"#            (IC)       (OC)       (RB)       (AB)\n";

  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<burn_nuc;i++){
    int matid=med[ic_s].GetNuclideInTurn(i).GetMatnum();
    string name=midt.Name(matid);
    real den0=med[ic_s].GetNuclideInTurn(i).GetDensity();
    real den1=med[oc_s].GetNuclideInTurn(i).GetDensity();
    real den2=med[rb_s].GetNuclideInTurn(i).GetDensity();
    real den3=med[abic_s].GetNuclideInTurn(i).GetDensity();
    if(den0+den1+den2+den3>0.){
      WriteOut(name,8);
      cout<<"  "<<den0<<" "<<den1<<" "<<den2<<" "<<den3<<"\n";
    };
  };
  cout<<"\n";
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
      for(int k=0;k<cycle_div;k++){
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
  cout<<"#Cy  Capture    Fission    (n,2n)\n";
  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<trace_cycle;i++){
    WriteOut(i,3);
    cout<<"  ";
    cout<<sigc_1g[medid][i][nuct]<<"  ";
    cout<<sigf_1g[medid][i][nuct]<<"  ";
    cout<<sign2n_1g[medid][i][nuct]<<"  ";
    cout<<"\n";
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

// +++ Ogawa calculation
void FRBurnerRZ::OgawaCalculationOn(real y1, real y2, real y3)
{
  ogawa_cal=true;
  ogawa_year.resize(3);
  ogawa_year[0]=y1;
  ogawa_year[1]=y2;
  ogawa_year[2]=y3;
};
