#include <cstdlib>
#include "Burner.h"

Burner::Burner()
{
  group = 107;

  med.resize(3);
  for(int i=0;i<3;i++){
    med[i].PutImax(group);
    med[i].PutPL(1);
  };

  rr.resize(6);

  sub_step_org=20; 
};

void Burner::SetLibrarySimple(string cbglibdir)
{
  nucn=21+98+8;

  int matno[]={
    9225,9228,9231,9234,9237,9346,9352,9434,9437,9440,
    9443,9446,9543,9546,9547,9549,9631,9634,9637,9640,
    9643,
    3640,3646,3843,3928,4034, 4040,4043,4131,4234,4240,
    4243,4246,4249,4331,4437, 4440,4443,4446,4449,4452,
    4455,4525,4531,4731,4631, 4634,4637,4640,4643,4837,
    4840,4843,4846,4849,4855, 4725,4931,5067,5137,5140,
    5247,5325,5331,5337,5349, 5446,5449,5452,5455,5458,
    5461,5728,5731,5525,5528, 5531,5537,5837,5840,5849,
    5646,5649,5655,5925,5931, 6025,6028,6031,6034,6037,
    6040,6043,6049,6149,6152, 6155,6161,6234,6237,6240,
    6243,6246,6249,6325,6328, 6331,6334,6337,6340,6343,
    6425,6431,6434,6437,6440, 6443,6449,6153,
      //4126,4534,5141,5647,5934,
    825,4000,125,525,2600, 2400,2800,528
  };

  string libdir=cbglibdir+"/CBGLIB/j4.107g.iwt7/";
  xslib.Initialize(libdir,"N-ENERGY");
  string filename_hm[]={
    // (Heavy nuclides : 21 nuclides)
    "U234","U235.U8","U236","U237","U238.mix",
    "Np237","Np239","Pu238","Pu239.mix","Pu240.mix",
    "Pu241.mix","Pu242.mix","Am241.mix","Am242","Am242m",
    "Am243","Cm242","Cm243","Cm244.T2","Cm245",
    "Cm246",
  };
  string filename_fp[]={
    // (Fission products : 98 nuclides)
    "Kr083","Kr085","Sr090","Y090","Zr093",
    "Zr095","Zr096","Nb095","Mo095","Mo097",
    "Mo098","Mo099","Mo100","Tc099.T2","Ru100",
    "Ru101","Ru102","Ru103","Ru104","Ru105",
    "Ru106","Rh103","Rh105","Ag109","Pd104",
    "Pd105","Pd106","Pd107","Pd108","Cd110",
    "Cd111","Cd112","Cd113","Cd114","Cd116",
    "Ag107","In115","Sn126","Sb125","Sb126",
    "Te127m","I127","I129","I131","I135",
    "Xe131.T2","Xe132","Xe133","Xe134","Xe135",

    "Xe136","La139","La140","Cs133","Cs134",
    "Cs135","Cs137","Ce140","Ce141","Ce144",
    "Ba137","Ba138","Ba140","Pr141","Pr143",
    "Nd142","Nd143","Nd144","Nd145","Nd146",
    "Nd147","Nd148","Nd150","Pm147","Pm148",
    "Pm149","Pm151","Sm147","Sm148","Sm149",
    "Sm150.T2","Sm151.T2","Sm152.T2","Eu151","Eu152",
    "Eu153","Eu154","Eu155","Eu156","Eu157",
    "Gd152","Gd154","Gd155","Gd156","Gd157",
    "Gd158","Gd160","Pm148m",
  };
  string filename_str[]={
    // (light or medium-heavy nuclides : 8 nuclides)
    "O016","Zr000","H001","B010","Fe000",
    "Cr000","Ni000","B011"
  };

  xslib.ReadFile(21,libdir,filename_hm);
  xslib.ReadFile(98,libdir,filename_fp);
  xslib.ReadFile(8,libdir,filename_str);
  // (thermal data)
  string libth=libdir;
  libth.append("Thermal/");
  xslib.GetLibData(125).GetThScat().ReadFile(libth,"H.H2O");
  xslib.GetLibData(825).GetThScat().ReadFile(libth,"O");
  // (Bell factor-optimization)
  string belldir=cbglibdir+"/CBGLIB/cbg-107g/Bell/";
  xslib.ReadBellFactor(belldir,"U235.300K",9228);
  xslib.ReadBellFactor(belldir,"U238.300K",9237);
  xslib.ReadBellFactor(belldir,"Pu239.300K",9437);
  xslib.ReadBellFactor(belldir,"Pu240.300K",9440);
  xslib.ReadBellFactor(belldir,"Pu241.300K",9443);
  xslib.ReadBellFactor(belldir,"Pu242.300K",9446);
  
  // (Create medium data for fuel)
  real *den0in=new real[nucn];
  for(int i=0;i<nucn;i++){
    den0in[i]=0.;
  };
  med[0].PutNuclide(nucn,matno,den0in);

  delete [] den0in;
};

void Burner::SetLibrary(string cbglibdir)
{
  nucn=21+185+13+8; // !! VERY IMPORTANT VALIABLE

  int matno[]={
    9225,9228,9231,9234,9237, 9346,9352,9434,9437,9440,
    9443,9446,9543,9546,9547, 9549,9631,9634,9637,9640,
    9643,
    3234,3237,3243,3325,3431, 3434,3437,3440,3443,3449,
    3531,3637,3640,3643,3646, 3649,3725,3728,3731,3831,
    3834,3837,3840,3843,3925, 3928,3931,4025,4028,4031,
    4034,4037,4040,4043,4125, 4128,4131,4225,4231,4234,
    4237,4240,4243,4246,4249, 4331,4437,4440,4443,4446,
    4449,4452,4455,4525,4531, 4631,4634,4637,4640,4643,
    4649,4725,4731,4735,4837, 4840,4843,4846,4849,4855,
    4925,4931,5037,5040,5043, 5046,5049,5055,5058,5061,
    5067,5125,5131,5134,5137, 5140,5231,5234,5237,5240,
    5243,5247,5249,5253,5255, 5261,5325,5331,5334,5337,
    5349,5431,5437,5440,5443, 5446,5449,5452,5455,5458,
    5461,5525,5528,5531,5534, 5537,5637,5640,5643,5646,
    5649,5655,5728,5731,5837, 5840,5843,5846,5849,5925,
    5931,6025,6028,6031,6034, 6037,6040,6043,6049,6149,
    6152,6153,6155,6161,6234, 6237,6240,6243,6246,6249,
    6252,6255,6325,6328,6331, 6334,6337,6340,6343,6425,
    6431,6434,6437,6440,6443, 6449,6525,6528,6637,6640,
    6643,6646,6649,6725,6825, 6831,6837,6840,6843,6849,
    7231,7234,7237,7240,7243,
    4126,4534,4847,5047,5052, 5053,5141,5235,5241,5647,
    5934,6719,6729,
    825,4000,125,525,2600, 2400,2800,528
  };
 
  string libdir=cbglibdir+"/CBGLIB/j4.107g.iwt7/";
  string libdir2=cbglibdir+"/CBGLIB/tdl11.107g.iwt7/";
  xslib.Initialize(libdir,"N-ENERGY");
  string filename_hm[]={
    // (Heavy nuclides : 21 nuclides)
    "U234","U235.U8","U236","U237","U238.mix",
    "Np237","Np239","Pu238","Pu239.mix","Pu240.mix",
    "Pu241.mix","Pu242.mix","Am241.mix","Am242","Am242m",
    "Am243","Cm242","Cm243","Cm244.T2","Cm245",
    "Cm246",
  };
  string filename_fp1[]={
    // (Fission products : 185 nuclides)
    "Ge073","Ge074","Ge076","As075","Se076",
    "Se077","Se078","Se079","Se080","Se082",
    "Br081","Kr082","Kr083","Kr084","Kr085",
    "Kr086","Rb085","Rb086","Rb087","Sr086",
    "Sr087","Sr088","Sr089","Sr090","Y089",
    "Y090","Y091","Zr090","Zr091","Zr092",
    "Zr093","Zr094","Zr095","Zr096","Nb093",
    "Nb094","Nb095","Mo092","Mo094","Mo095",
    "Mo096","Mo097","Mo098","Mo099","Mo100",
    "Tc099.T2","Ru100","Ru101","Ru102","Ru103",
    "Ru104","Ru105","Ru106","Rh103","Rh105",
    "Pd104","Pd105","Pd106","Pd107","Pd108",
    "Pd110","Ag107","Ag109","Ag110m","Cd110",
    "Cd111","Cd112","Cd113","Cd114","Cd116",
    "In113","In115","Sn116","Sn117","Sn118",
    "Sn119","Sn120","Sn122","Sn123","Sn124",
    "Sn126","Sb121","Sb123","Sb124","Sb125",
    "Sb126","Te122","Te123","Te124","Te125",
    "Te126","Te127m","Te128","Te129m","Te130",
    "Te132","I127","I129","I130","I131",
    "I135","Xe126","Xe128","Xe129","Xe130",
    "Xe131.T2","Xe132","Xe133","Xe134","Xe135",
    "Xe136","Cs133.T2","Cs134","Cs135","Cs136",
    "Cs137","Ba134","Ba135","Ba136","Ba137",
    "Ba138","Ba140","La139","La140","Ce140",
    "Ce141","Ce142","Ce143","Ce144","Pr141",
    "Pr143","Nd142","Nd143","Nd144","Nd145",
    "Nd146","Nd147","Nd148","Nd150","Pm147",
    "Pm148","Pm148m","Pm149","Pm151","Sm147",
    "Sm148","Sm149","Sm150.T2","Sm151.T2","Sm152.T2",
    "Sm153","Sm154","Eu151","Eu152","Eu153",
    "Eu154","Eu155","Eu156","Eu157","Gd152",
    "Gd154","Gd155","Gd156","Gd157","Gd158",
    "Gd160","Tb159","Tb160","Dy160","Dy161",
    "Dy162","Dy163","Dy164","Ho165","Er162",
    "Er164","Er166","Er167","Er168","Er170",
    "Hf176","Hf177","Hf178","Hf179","Hf180",
  };
  string filename_str[]={
    // (light or medium-heavy nuclides : 8 nuclides)
    "O016","Zr000","H001","B010","Fe000",
    "Cr000","Ni000","B011"
  };
  string filename_fp2[]={
    "Nb093m","Rh106","Cd113m","Sn119m","Sn121",
    "Sn121m","Sb126m","Te123m","Te125m","Ba137m",
    "Pr144","Ho163","Ho166m",
  };

  xslib.ReadFile(21,libdir,filename_hm);
  xslib.ReadFile(185,libdir,filename_fp1);
  xslib.ReadFile(8,libdir,filename_str);
  xslib.ReadFile(13,libdir2,filename_fp2);
  // (thermal data)
  string libth=libdir;
  libth.append("Thermal/");
  xslib.GetLibData(125).GetThScat().ReadFile(libth,"H.H2O");
  xslib.GetLibData(825).GetThScat().ReadFile(libth,"O");
  // (Bell factor-optimization)
  string belldir=cbglibdir+"/CBGLIB/cbg-107g/Bell/";
  xslib.ReadBellFactor(belldir,"U235.300K",9228);
  xslib.ReadBellFactor(belldir,"U238.300K",9237);
  xslib.ReadBellFactor(belldir,"Pu239.300K",9437);
  xslib.ReadBellFactor(belldir,"Pu240.300K",9440);
  xslib.ReadBellFactor(belldir,"Pu241.300K",9443);
  xslib.ReadBellFactor(belldir,"Pu242.300K",9446);
  
  // (Create medium data for fuel)
  real *den0in=new real[nucn];
  for(int i=0;i<nucn;i++){
    den0in[i]=0.;
  };
  med[0].PutNuclide(nucn,matno,den0in);

  delete [] den0in;
};

void Burner::PutFuelData(int nuc,int *mat,real *den,real temp)
{
  for(int i=0;i<nuc;i++){
    med[0].GetNuclide(mat[i]).PutDensity(den[i]);
  };

  med[0].PutTemperatureForAllNuclide(temp);
};

void Burner::PutCladData(int nuc,int *mat,real *den,real temp)
{
  med[1].PutNuclide(nuc,mat,den);
  med[1].PutTemperatureForAllNuclide(temp);
};

void Burner::PutModeratorData(int nuc,int *mat,real *den,real temp)
{
  med[2].PutNuclide(nuc,mat,den);
  med[2].PutTemperatureForAllNuclide(temp);
};

void Burner::PutGeometryData(real pitch,real *rri)
{
  pin_pitch=pitch;
  for(int i=0;i<6;i++){rr[i]=rri[i];};
  fuel_r=rr[3]; 
  clad_r=rr[2]; 
  fuel_vol=fuel_r*fuel_r*PI;
};

void Burner::PutBurnStep(int i)
{
  burn_step=i;
  power_density_list.resize(burn_step);
  burn_time.resize(burn_step);
  sub_step_list.resize(burn_step+1);
  keff.resize(burn_step+1);
  acday.resize(burn_step+1);
  acburn.resize(burn_step+1);
};

void Burner::PutPowerDensityList(real *inp)
{
  for(int i=0;i<burn_step;i++){
    power_density_list[i]=inp[i];
  };

  for(int i=0;i<burn_step;i++){
    sub_step_list[i]=sub_step_org;
    if(power_density_list[i]<1e-5){
      power_density_list[i]=0.;
      //sub_step_list[i]=1;
    };
  };

};

void Burner::PutBurnTime(real *inp,bool GWd_t,bool accumulate)
{
  for(int i=0;i<burn_step;i++){
    burn_time[i]=inp[i];
  };
  burn_time_GWd_t=GWd_t;
  burn_time_accumulate=accumulate;
};

void Burner::Calculation(Burnup &bu)
{
  setlinebuf(stdout);

  // +++ boundary condition for Pij calculation
  enum BCondition bc_ssc=White;     // (self-shielding calculation)
  //enum BCondition bc_flx=Periodic;  // (eigenvalue calculation)
  enum BCondition bc_flx=White;  // (eigenvalue calculation)

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ IrregularGeometryInformation & TrajectorySet +++++++++++++++++++++++++++++++++++
  //     (for self-shielding calculation)
  IrregularGeometryInformation igi;
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
  //
  TrajectorySet sys;
  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,8,0.02,45.);
  //     (for flux distribution calculation)
  int mesh_fuel=3; // number of mesh in fuel region
  IrregularGeometryInformation igi_f;
  GeomPolygon pol2;
  pol2.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol2.PutRegionID(6);
  igi_f.AddGeom(pol2);
  int rid[]={5,4,3,2,1,0};
  igi_f.AddCircleRing(6,rr,rid);
  //
  TrajectorySet sys_f;
  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi_f,8,0.02,45.);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Pre self-shielding calculation +++++++++++++++++++++++++++++++++++++++++++++++++
  OnePointCalculator opc;
  opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[1],xslib);
  opc.GiveInfiniteDillutionCrossSection(med[2],xslib);
  opc.CalThermalScatteringMatrix(med[2],xslib, 3.93); // 3.93 : thermal cut-off energy
  med[2].CalMacroFromMicro();
  med[2].CalSigtr(0);
  //
  GroupData1D c(group); // Dancoff correction
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  GeneralOption opt;
  int region_medium[]={0,0,0,1,2,2,2};
  int region_medium2[]={0,1,2};

  // +++ Number density data storing ++++++++++++++++++++++++++++++++++++++++++++++++++++
  density_data.resize(burn_step+1);

  //
  hm_weight_init=bu.CalWeightOfHeavyNuclideParUnitVolume(med[0])*fuel_vol; // [g] (NOT [g/cm3])

  if(burn_time_accumulate){
    real tmp=0.;
    for(int i=burn_step-1;i>0;i--){
      if(power_density_list[i]>0.){
	bool tag=false;
	for(int j=i-1;j>=0;j--){
          if(power_density_list[j]>0.&&!tag){
	    tag=true;
            burn_time[i]-=burn_time[j];
	  };
	};
      };
    };
  };

  vector<real> burn_time_gwd(burn_step);
  if(burn_time_GWd_t){
    for(int i=0;i<burn_step;i++){
      burn_time_gwd[i]=burn_time[i];
      if(power_density_list[i]==0.){
        burn_time_gwd[i]=0.;
      }else{
        burn_time[i]*=hm_weight_init*1e-6/1e-9/power_density_list[i];
      };
    };
  }else{
    for(int i=0;i<burn_step;i++){
      if(power_density_list[i]==0.){
        burn_time_gwd[i]=0.;      
      }else{
        burn_time_gwd[i]=burn_time[i]/(hm_weight_init*1e-6/1e-9/power_density_list[i]);
      };
    };
  };

  real accumulated_day=0.;
  real accumulated_burn=0.;
  for(int st=0;st<burn_step+1;st++){

    acday[st]=accumulated_day;
    acburn[st]=accumulated_burn;

    cout<<"# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";

    // 
    density_data[st].resize(nucn);
    for(int i=0;i<nucn;i++){
      density_data[st][i]=med[0].GetNuclideInTurn(i).GetDensity();
    };

    cout<<"#\n ... self-shielding calculation ...\n";
    // (Self-shielding calculation)
    int matmax=2;
    if(st!=0)matmax=1; // clad region is calculated only at t=0
    for(int mat=0;mat<matmax;mat++){
      real r;
    if(mat==0){r=fuel_r;}else{r=clad_r-fuel_r;};
    for(int i=0;i<group;i++){
      real xs[3];
      for(int j=0;j<3;j++){
        xs[j]=med[j].GetMacxs().GetData1d(sigt).get_dat(i);
      };
      xs[mat]=30000.;
      real pij[3*3];
      sys.CalculationPij(xs,pij,false);
      real pesc=1.-pij[mat*3+mat];
      real dancoff_corr=1.-xs[mat]*2.*r*pesc;
      c.put_data(i,dancoff_corr);
    };
    opc.CalSelfShieldingWithDancoffCorrection(med[mat],xslib,r*2.,b,c);
    med[mat].CalSigtr(0);
    };
    // Thermal scattering matrices are over-written
    opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    cout<<"#      ... end\n";

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
    keff[st]=lat.CalIgenPij();
    real vol_inv=1./fuel_vol;
    GroupData1D flx=lat.GetIntegratedFluxMeshID(0,2)*vol_inv; // Total-flux per unit volume in fuel region
    med[0].GetFlux().copy(flx);

    /*
    cout<<"# Neutron flux energy spectrum \n";
    for(int g=0;g<107;g++){
      real e0=med[0].GetEnband().get_dat(g);
      real e1=med[0].GetEnband().get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<e0<<" "<<flx.get_dat(g)/letwid<<"\n";
    };
    cout<<"\n\n";
    */

    bu.PutMediumData(med[0]); // burnup data

    // +++ one-group cross section printing
    /*
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# One-group XS (fission/capture/n2n)\n";
    for(int j=0;j<nucn;j++){
      cout<<filename[j]<<"   ";
      cout<<bu.GetSigf(j)<<" ";
      cout<<bu.GetSigc(j)<<" ";
      cout<<bu.GetSign2n(j)<<" ";
      cout<<"\n";
    };
    */

    if(st!=burn_step){

      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day

      int sub_step=sub_step_list[st];
      burn_span/=sub_step;
    
      // +++ Burnup calculation
      accumulated_burn+=burn_time_gwd[st];
      for(int j=0;j<sub_step;j++){
	cout<<"#... burnup calculation : "<<j<<"/"<<sub_step<<"\n";   
        real power_factor=power_density/(bu.GetIntegratedPower(med[0])*fuel_vol);
        real total_flux=flx.get_sum()*power_factor;
	accumulated_day+=burn_span;
	real delt=burn_span*24*60*60;
        bu.BurnupCalculationByKrylov(med[0],total_flux,delt,false);
      };
    };
  };
};

void Burner::ShowDecayHeat(Burnup &bu)
{
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"#\n# Time-dependent decay heat (day, W/cm , W/tHM, ([HM][FP]))\n#\n";
  real ev_to_j=1.60219e-19;
  real ent;
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
      if(id>=9000){en_hm+=enuc;}else{en_fp+=enuc;};
    };
    en*=ev_to_j; // J
    en_hm*=ev_to_j; 
    en_fp*=ev_to_j;
    cout<<en<<" ( ";
    cout<<en_hm<<" ";
    cout<<en_fp<<" ) ";
    cout<<en/(hm_weight_init*1e-6)<<" ( "; // J/tHM
    cout<<en_hm/(hm_weight_init*1e-6)<<" ";
    cout<<en_fp/(hm_weight_init*1e-6)<<" ) ";
    cout<<"\n";
    if(i==burn_step)ent=en;
  };

  /*
  int ist=21;
  for(int j=0;j<nucn;j++){
    int id=med[0].GetNuclideInTurn(j).GetMatnum();
    real den=density_data[ist][j]*1e24;
    real lambda=bu.GetDecayConstant(id);
    real e=0.;
    for(int k=0;k<3;k++){
      e+=bu.GetBurnupChain().GetDecayEnergy(id,k);
    };
    real en=e*den*lambda;
    en*=ev_to_j; // J
    if(en/ent>0.01)cout<<id<<" "<<eidt.Name(id)<<" "<<en/ent<<"\n";
  };
  */
};

void Burner::ShowEigenvalue()
{
  cout<<"#\n# +++ Time-dependent eigenvalue +++\n";
  cout.setf(ios::showpoint);
  cout.precision(6);
  cout<<"# (day)   (GWd/t)  (keff)\n";
  for(int i=0;i<burn_step+1;i++){
    cout<<i<<" "<<acday[i]<<" "<<acburn[i]<<" "<<keff[i]<<"\n";
  };
};

void Burner::ShowNumberDensityChange()
{
  cout<<"#+++ Number density before/after burnup +++\n";
  for(int i=0;i<nucn;i++){
    cout.setf(ios::scientific);
    cout.precision(5);
    real den_b=density_data[0][i];
    real den_a=density_data[burn_step][i];
    int matno=med[0].GetNuclideInTurn(i).GetMatnum();
    if(den_a!=0.||den_b!=0.){
      cout<<"# "<<matno<<" "<<eidt.Name(matno);
      cout<<"   "<<den_b<<"   "<<den_a<<"\n";
    };
  };
};

void Burner::ShowNumberDensityHistory(int prt_nuc,string *prt_nuc_nam,Burnup &bu,string opt)
{
  real avo=0.60221367;

  vector<int> prt_nuc_turn(prt_nuc);
  vector<int> prt_nuc_id(prt_nuc);
  for(int i=0;i<prt_nuc;i++){
    int id=eidt.ID(prt_nuc_nam[i]);
    prt_nuc_id[i]=id;
    for(int j=0;j<nucn;j++){
      if(med[0].GetNuclideInTurn(j).GetMatnum()==id){
	prt_nuc_turn[i]=j;
	break;
      };
    };
  };

  if(opt=="nd_per_vol"){
    cout<<"\n\n#  (day, GWd/t, N.D. [1e24/cm3]).\n";
  }else if(opt=="nd"){
    cout<<"\n\n#  (day, GWd/t, N.D. [1e24]).\n";
  }else if(opt=="bq"){
    cout<<"\n\n#  (day, GWd/t, radioactivity [Bq]).\n";
  }else if(opt=="bq_per_thm"){
    cout<<"\n\n#  (day, GWd/t, radioactivity [Bq/tHM]).\n";
  }else if(opt=="kg_per_thm"){
    cout<<"\n\n#  (day, GWd/t, Weight [kg/tHM]).\n";
  };
  cout<<"# ";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int j=0;j<prt_nuc;j++){
    cout<<prt_nuc_nam[j]<<" ";
  };
  cout<<"\n";
  for(int i=0;i<burn_step+1;i++){
    cout<<acday[i]<<" "<<acburn[i]<<" ";
    real sum=0.;
    for(int j=0;j<prt_nuc;j++){
      real dd=density_data[i][prt_nuc_turn[j]];
      real dc=bu.GetBurnupChain().GetDecayConstant(prt_nuc_id[j]);
      real val;
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
      };
      cout<<val<<" ";
      sum+=val;
    };
    cout<<" ("<<sum<<")\n";
  };
};

void Burner::SensitivityCalculation()
{
  //
  // Burnup-sensitivity calculation
  //

  setlinebuf(stdout);

  int group = 107;
 
  // +++ GEOMETRY DATA BLOCK ++++++++++++++++++++++++++
  real pin_pitch=0.6325*2.; // [cm]
  real rr[]={0.6,0.55,0.476,0.412,0.35,0.25};
  //                    r2    r1  
  // r1 : fuel region outer radius [cm]
  // r2 : clading region outer radius [cm]

  // +++ BOUNDARY CONDITION BLOCK +++++++++++++++++++++
  enum BCondition bc_ssc=White;     // (self-shielding calculation)
  //enum BCondition bc_flx=Periodic;  // (eigenvalue calculation)
  enum BCondition bc_flx=White;  // (eigenvalue calculation)

  // +++ MATERIAL DATA BLOCK ++++++++++++++++++++++++++
  // (fuel)
  //   ... MOX
  int nuc0=9;
  int mat0[]={9228,9237,9434,9437,9440, 9443,9446,9543,825};
  real den0[]={
    4.104e-5, 2.022e-2, 4.731e-5, 1.223e-3, 5.586e-4,
    2.069e-4, 1.418e-4, 6.007e-5, 4.500e-2
  };
  /*
  //   ... UO2
  int nuc0=3;
  int mat0[]={9228,9237,825};
  real den0[]={9.349e-4, 2.159e-2, 4.505e-2};
  */
  real temp0=968.8;
  // (clading)
  int nuc1=3;
  int mat1[]={4000,2600,2400};
  real den1[]={3.786e-2, 2.382e-4, 6.770e-5};
  real temp1=604.;
  // (moderator)
  int nuc2=6;
  int mat2[]={125,825,525,2800,2400,2600};
  real den2[]={5.572e-2, 2.786e-2, 4.592e-6, 3.688e-4, 1.609e-4, 1.306e-4};
  real temp2=574.2;

  // +++ BURNUP DATA BLOCK ++++++++++++++++++++++++++++++;
  int burn_step=5;
  //int burn_step=21;
  real power_density_list[]={
    179., 179., 179., 179., 179.,
    179., 179., 179., 179., 179.,
    179., 179., 179., 179., 179.,
    179., 179., 179., 179., 179.,
    0.,
  }; // [W/cm]
  real burn_time[]={
     0.1, 1.0, 2.5, 5.0, 7.5,  10., 12.5, 15., 17.5, 20.,
     22.5, 25., 27.5, 30., 32.5,  35., 37.5, 40., 42.5, 45.,
     45.
  }; // GWd/t
  bool burn_time_GWd_t=true; // default unit is "Day"
  bool burn_time_accumulate=true;
  real cool_year=5.; // It is used for burn_time=0 GWd/t
  int sub_step_org=20;

  // +++ SENSITIVITY CALCULATION OPTION BLOCK +++++++++++++++++++++;
  bool sns_cal=true;   // Sensitivity calculation
  bool snsdcy_cal=true; // Sensitivity calculation for decay constant
  bool snsyld_cal=true; // Sensitivity calculation for fission yield
  int target_nuclide_num=1;
  string target_nuclide_snscal[]={
    //"U234","U235.U8","U236","U237","U238.mix",
    //"Np237","Np239","Pu238","Pu239.mix","Pu240.mix",
    //"Pu241.mix","Pu242.mix","Am241.mix","Am242","Am242m",
    //"Am243","Cm242","Cm243","Cm244.T2","Cm245",
    //"Cm246",
    "Cm242","Cm243","Cm244.T2","Cm245","Cm246"
  };
  int gptitermax=10; 

  // +++ DIRECT PERTURBATION CALCULATION OPTION BLOCK +++++++++++++++
  bool pertxs=false;
  string pert_nuclide="Pu239.mix";
  int pertgrp1=0;
  int pertgrp2=106;
  real amp=0.01;
  enum xstype ttype=sigf;

  //
  //
  // +++ xs library +++++++++++++++++++++++++++++++++++++++++++++++

  string filename[]={
    // (Heavy nuclides : 21 nuclides)
    "U234","U235.U8","U236","U237","U238.mix",
    "Np237","Np239","Pu238","Pu239.mix","Pu240.mix",
    "Pu241.mix","Pu242.mix","Am241.mix","Am242","Am242m",
    "Am243","Cm242","Cm243","Cm244.T2","Cm245",
    "Cm246",
    // (Fission products : 98 nuclides)
    "Kr083","Kr085","Sr090","Y090","Zr093",
    "Zr095","Zr096","Nb095","Mo095","Mo097",

    "Mo098","Mo099","Mo100","Tc099.T2","Ru100",
    "Ru101","Ru102","Ru103","Ru104","Ru105",

    "Ru106","Rh103","Rh105","Ag109","Pd104",
    "Pd105","Pd106","Pd107","Pd108","Cd110",

    "Cd111","Cd112","Cd113","Cd114","Cd116",
    "Ag107","In115","Sn126","Sb125","Sb126",

    "Te127m","I127","I129","I131","I135",
    "Xe131.T2","Xe132","Xe133","Xe134","Xe135",

    "Xe136","La139","La140","Cs133","Cs134",
    "Cs135","Cs137","Ce140","Ce141","Ce144",

    "Ba137","Ba138","Ba140","Pr141","Pr143",
    "Nd142","Nd143","Nd144","Nd145","Nd146",
    "Nd147","Nd148","Nd150","Pm147","Pm148",
    "Pm149","Pm151","Sm147","Sm148","Sm149",
    "Sm150.T2","Sm151.T2","Sm152.T2","Eu151","Eu152",
    "Eu153","Eu154","Eu155","Eu156","Eu157",
    "Gd152","Gd154","Gd155","Gd156","Gd157",
    "Gd158","Gd160","Pm148m",
    // (light or medium-heavy nuclides : 8 nuclides)
    "O016","Zr000","H001","B010","Fe000",
    "Cr000","Ni000","B011"
  };
  int matno[]={
    9225,9228,9231,9234,9237, 9346,9352,9434,9437,9440,
    9443,9446,9543,9546,9547, 9549,9631,9634,9637,9640,
    9643,
    3640,3646,3843,3928,4034, 4040,4043,4131,4234,4240,
    4243,4246,4249,4331,4437, 4440,4443,4446,4449,4452,
    4455,4525,4531,4731,4631, 4634,4637,4640,4643,4837,
    4840,4843,4846,4849,4855, 4725,4931,5067,5137,5140, // Sb126
    5247,5325,5331,5337,5349, 5446,5449,5452,5455,5458, // Xe135
    5461,5728,5731,5525,5528, 5531,5537,5837,5840,5849, // Ce144
    5646,5649,5655,5925,5931, 6025,6028,6031,6034,6037,
    6040,6043,6049,6149,6152, 6155,6161,6234,6237,6240,
    6243,6246,6249,6325,6328, 6331,6334,6337,6340,6343,
    6425,6431,6434,6437,6440, 6443,6449,6153,
    825,4000,125,525,2600, 2400,2800,528
  };
  int nucn=21+98+8;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // (identification of nuclide)
  // (perturbation to library)
  if(pertxs){
  int pertmatid=-1;
    for(int i=0;i<nucn;i++){
      if(filename[i]==pert_nuclide)pertmatid=matno[i];
    };
    if(pertmatid==-1){
      cout<<"ERROR!\n";
      exit(0);
    };
    for(int g=pertgrp1;g<=pertgrp2;g++){
      real org=xslib.GetLibData(pertmatid).GetXSData().GetData1d(ttype).get_dat(g);
      xslib.GetLibData(pertmatid).GetXSData().GetData1d(ttype).add_data(g,org*amp);
    };
  };

  //vector<int> sub_step_list(burn_step+1);
  vector<int> sub_step_list(burn_step);
  for(int i=0;i<burn_step;i++){
    sub_step_list[i]=sub_step_org;
    if(power_density_list[i]<1e-5){
      power_density_list[i]=1e-10;
      if(!sns_cal)sub_step_list[i]=1; 
      // In sensitivity calculation, large number of sub-steps may be necessary
      // because time-integration for N* N is performed.
    };
  };

  real fuel_r=rr[3]; 
  real clad_r=rr[2]; 
  real fuel_vol=fuel_r*fuel_r*PI;

  // +++ lattice data +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Medium pin[3];
  for(int i=0;i<3;i++){
    pin[i].PutImax(group);
    pin[i].PutPL(1);
  };

  real *den0in=new real[nucn];
  for(int i=0;i<nucn;i++){
    den0in[i]=0.;
    for(int j=0;j<nuc0;j++){
      if(matno[i]==mat0[j])den0in[i]=den0[j];
    };
  };
  // (fuel region)
  pin[0].PutNuclide(nucn,matno,den0in);
  pin[0].PutTemperatureForAllNuclide(temp0);
  // (clading region)
  pin[1].PutNuclide(nuc1,mat1,den1);
  pin[1].PutTemperatureForAllNuclide(temp1);
  // (moderator region)
  pin[2].PutNuclide(nuc2,mat2,den2);
  pin[2].PutTemperatureForAllNuclide(temp2);
  delete [] den0in;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ IrregularGeometryInformation & TrajectorySet +++++++++++++++++++++++++++++++++++
  //     (for self-shielding calculation)
  IrregularGeometryInformation igi;
  GeomPolygon pol;
  pol.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol.PutRegionID(2);
  igi.AddGeom(pol);
  GeomCircle cir1(0.,0.,clad_r);
  cir1.PutRegionID(1);
  igi.AddGeom(cir1);
  real r=fuel_r;
  GeomCircle cir2(0.,0.,r);
  cir2.PutRegionID(0);
  igi.AddGeom(cir2);
  //
  TrajectorySet sys;
  sys.PutBoundaryCondition(bc_ssc);
  sys.CalTrajectory(igi,8,0.02,45.);
  //     (for flux distribution calculation)
  int mesh_fuel=3; // number of mesh in fuel region
  IrregularGeometryInformation igi_f;
  GeomPolygon pol2;
  pol2.PutRectangular(0.,0.,pin_pitch*0.5,pin_pitch*0.5);
  pol2.PutRegionID(6);
  igi_f.AddGeom(pol2);
  int rid[]={5,4,3,2,1,0};
  igi_f.AddCircleRing(6,rr,rid);
  //
  TrajectorySet sys_f;
  sys_f.PutBoundaryCondition(bc_flx);
  sys_f.CalTrajectory(igi_f,8,0.02,45.);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Pre self-shielding calculation +++++++++++++++++++++++++++++++++++++++++++++++++
  OnePointCalculator opc;
  opc.GiveInfiniteDillutionCrossSection(pin[0],xslib);
  opc.GiveInfiniteDillutionCrossSection(pin[1],xslib);
  opc.GiveInfiniteDillutionCrossSection(pin[2],xslib);
  opc.CalThermalScatteringMatrix(pin[2],xslib, 3.93); // 3.93 : thermal cut-off energy
  pin[2].CalMacroFromMicro();
  pin[2].CalSigtr(0);
  //
  GroupData1D c(group); // Dancoff correction
  GroupData1D b(group); // Bell factor
  for(int i=0;i<group;i++){b.put_data(i,1.2);};
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  GeneralOption opt,opta;
  opta.PutAdjointCal();
  int region_medium[]={0,0,0,1,2,2,2};
  int region_medium2[]={0,1,2};

  // +++ Burnup Calculation +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Burnup bu; // "true" means usage of SRAC fine burnup chain.
  bu.PutSracTypeReactionEnergyData(); 
  bu.SetOKUMURAChain();
  bu.GetBurnupChain().ReadDecayConstantFromFile("/home/chiba/CBGLIB/decay_constant/","srac_org");
  //
  real hm_weight_init=bu.CalWeightOfHeavyNuclideParUnitVolume(pin[0])*fuel_vol; // g
  if(burn_time_accumulate){
    for(int i=burn_step-1;i>0;i--){
      burn_time[i]-=burn_time[i-1];
    };
  };
  if(burn_time_GWd_t){
    for(int i=0;i<burn_step;i++){
      burn_time[i]*=hm_weight_init*1e-6/1e-9/power_density_list[i];
      if(burn_time[i]==0.)burn_time[i]=cool_year*365.;
    };
  };

  vector<real> keff(burn_step+1);
  vector<real> acday(burn_step+1);
  vector< vector<real> > power_factor(burn_step);
  vector< vector<real> > delt(burn_step);
  vector< vector<real> > total_flux(burn_step);
  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    delt[i].resize(sub_step);
    power_factor[i].resize(sub_step+1);
    total_flux[i].resize(sub_step);
  };

  // +++ for sensitivity calculation +++++++++++++++++++++++++++++++++++++++
  Burnup bu_dmdf;
  bu_dmdf.SetOKUMURAChain(); 
  bu_dmdf.PutSracTypeReactionEnergyData();
  bu_dmdf.GetBurnupChain().ReadDecayConstantFromFile("/home/chiba/CBGLIB/decay_constant/","srac_org");
  //
  vector< vector<GroupData1D> > fwd_nuc(burn_step);
  vector<GroupData2D> trmat_flxdep(burn_step);
  GroupData2D trmat_flxindep;
  vector< vector< vector<GroupData1D> > > dmdf_nuc(burn_step);  // (dM/dphi) * (NUC_FWD)
  vector<GroupData1D> vol_flx(burn_step);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];
    fwd_nuc[i].resize(sub_step+1);
    dmdf_nuc[i].resize(sub_step+1);
    for(int j=0;j<sub_step+1;j++){
      fwd_nuc[i][j].put_imax(nucn);
      dmdf_nuc[i][j].resize(group);
    };
  };
  //
  vector< vector<GroupDataSet> > macxs(burn_step);
  vector< vector<GroupData1D> > mic_sigf(burn_step);
  vector< vector<GroupData1D> > mic_sigc(burn_step);
  vector< vector<GroupData1D> > fwd_vol_flx(burn_step);
  for(int i=0;i<burn_step;i++){
    fwd_vol_flx[i].resize(mesh_fuel);
    int sub_step=sub_step_list[i];
    macxs[i].resize(sub_step);
    mic_sigf[i].resize(nucn);
    mic_sigc[i].resize(nucn);
    for(int j=0;j<sub_step;j++){
      macxs[i][j].Init("MacroCrossSection");
    };
  };
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real accumulated_day=0.;
  for(int st=0;st<burn_step+1;st++){

    acday[st]=accumulated_day;

    cout<<"# +++ Burnup step : "<<st<<"\n";
    cout<<"#     (accumulated day : "<<accumulated_day<<" )\n";
    cout<<"#\n ... self-shielding calculation ...\n";

    // (Self-shielding calculation)
    int matmax=2;
    if(st!=0)matmax=1; // clad region is calculated only at t=0
    for(int mat=0;mat<matmax;mat++){
      if(mat==0){r=fuel_r;}else{r=clad_r-fuel_r;};
      for(int i=0;i<group;i++){
        real xs[3];
        for(int j=0;j<3;j++){
          xs[j]=pin[j].GetMacxs().GetData1d(sigt).get_dat(i);
        };
        xs[mat]=30000.;
        real pij[3*3];
        sys.CalculationPij(xs,pij,false);
        real pesc=1.-pij[mat*3+mat];
        real dancoff_corr=1.-xs[mat]*2.*r*pesc;
        c.put_data(i,dancoff_corr);
      };
      opc.CalSelfShieldingWithDancoffCorrection(pin[mat],xslib,r*2.,b,c);
      pin[mat].CalSigtr(0);
    };
    // Thermal scattering matrices are over-written
    opc.CalThermalScatteringMatrix(pin[0],xslib,3.93);

    // +++ Eigenvalue calculation
    PJISystem lat(group,3);
    lat.PutTrajectorySet(&sys_f);
    lat.AddMedium(pin[0]);
    lat.AddMedium(pin[1]);
    lat.AddMedium(pin[2]);
    lat.PutRegMed(region_medium);
    lat.PutGeneralOption(opt);
    lat.PutSigmaCol(sigtr);
    lat.PutPij();
    keff[st]=lat.CalIgenPij();

    /*
    cout<<"# Neutron flux energy spectrum \n";
    for(int g=0;g<107;g++){
      real e0=pin[0].GetEnband().get_dat(g);
      real e1=pin[0].GetEnband().get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<e0<<" "<<flx.get_dat(g)/letwid<<"\n";
    };
    cout<<"\n\n";
    */

    if(st!=burn_step){

      real vol_inv=1./fuel_vol;
      vol_flx[st]=lat.GetIntegratedFluxMeshID(0,2); // Volume-integrated total flux in fuel region
      pin[0].GetFlux().copy(vol_flx[st]*vol_inv); // flux per unit volume
      bu.PutMediumData(pin[0]); // burnup data
      // +++ (for sensitivity calculation) +++++++++++++++++++++++++
      if(sns_cal){
      trmat_flxdep[st]=bu.GetTrmatFlxDep();
      if(st==0)trmat_flxindep=bu.GetTrmatFlxInDep();
      for(int i=0;i<mesh_fuel;i++){
        fwd_vol_flx[st][i]=lat.GetMesh(i).GetFlux()*lat.GetMesh(i).GetVolume();
      };
      for(int k=0;k<nucn;k++){
        mic_sigf[st][k].copy(pin[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigf));
	mic_sigc[st][k].copy(pin[0].GetNuclideInTurn(k).GetMicxs().GetData1d(sigc));
      };
      };
      vector<GroupData2D> dTMatdFlux(group);
      for(int g=0;g<group;g++){
	dTMatdFlux[g]=bu_dmdf.CaldTMatdFlux(pin[0],g);
      };
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real power_density=power_density_list[st];
      real burn_span=burn_time[st]; // day
      int sub_step=sub_step_list[st];
      burn_span/=sub_step;
      // +++ Burnup calculation
      for(int j=0;j<sub_step+1;j++){
        for(int i=0;i<nucn;i++){
	  fwd_nuc[st][j].put_data(i,pin[0].GetNuclideInTurn(i).GetDensity());
        };
        power_factor[st][j]=power_density/(bu.GetIntegratedPower(pin[0])*fuel_vol);
	if(j!=sub_step){
  	  cout<<"#... burnup calculation : "<<j<<"/"<<sub_step<<"\n";   
	  // +++ (for SNSCAL) ++++++++++++++++++++++++++++
	  if(sns_cal){
          if(j!=0)pin[0].CalMacroFromMicro();
          //macxs[st][j].DataCopy(pin[0].GetMacxs());
          macxs[st][j].DataCopyPL(pin[0].GetMacxs(),0);
          for(int jj=0;jj<group;jj++){
            //dmdf_nuc[st][j][jj]=bu_dmdf.CaldTMatdFlux(pin[0],jj)*fwd_nuc[st][j];
            dmdf_nuc[st][j][jj]=dTMatdFlux[jj]*fwd_nuc[st][j];
          };
	  };
	  // +++++++++++++++++++++++++++++++++++++++++++++
          total_flux[st][j]=vol_flx[st].get_sum()*power_factor[st][j]/fuel_vol;
	  accumulated_day+=burn_span;
	  delt[st][j]=burn_span*24*60*60;
          //bu.BurnupCalculationByPade(pin[0],total_flux[st][j],delt[st][j],false);
          bu.BurnupCalculationByKrylov(pin[0],total_flux[st][j],delt[st][j],false);
	};
      };

    };
  };

  cout.setf(ios::showpoint);
  cout.precision(6);
  cout<<"# (day)   (keff)\n";
  for(int i=0;i<burn_step+1;i++){
    cout<<i<<" "<<acday[i]<<" "<<keff[i]<<"\n";
  };

  // +++ Sensitivity calculation part
  if(sns_cal){

        PJISystem lat_gpt(group,3);
        lat_gpt.PutTrajectorySet(&sys_f);
        lat_gpt.AddMedium(pin[0]);
        lat_gpt.AddMedium(pin[1]);
        lat_gpt.AddMedium(pin[2]);
        lat_gpt.PutRegMed(region_medium);
        lat_gpt.PutGeneralOption(opt);
        lat_gpt.PutSigmaCol(sigtr);
        lat_gpt.PutBuckling(0.);
        lat_gpt.WriteProcOff();
        lat_gpt.NoPrint();

    for(int iii=0;iii<target_nuclide_num;iii++){

    // (identification of nuclide)
    int target_nuclideID_inturn=-1;
    for(int i=0;i<nucn;i++){
      if(filename[i]==target_nuclide_snscal[iii])target_nuclideID_inturn=i;
    };
    if(target_nuclideID_inturn==-1){
      cout<<"ERROR!\n";
      exit(0);
    };

    vector< vector<GroupData1D> > adj_nuc(burn_step);
    vector< vector<real> > pow_adj(burn_step);
    vector< vector< vector<GroupData1D> > > gpt_flx(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      adj_nuc[i].resize(sub_step+1);
      pow_adj[i].resize(sub_step+1);
      for(int j=0;j<sub_step+1;j++){
	adj_nuc[i][j].put_imax(nucn);
      };
      gpt_flx[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
        gpt_flx[i][j].resize(mesh_fuel);
        for(int ii=0;ii<mesh_fuel;ii++){
	  gpt_flx[i][j][ii].put_imax(group);
        };
      };
    };
 
    // (Adjoint burnup calculation)
    for(int st=burn_step-1;st>=0;st--){

      real power_density=power_density_list[st];
      int sub_step=sub_step_list[st];

      if(st==burn_step-1){
        adj_nuc[st][sub_step].put_data(target_nuclideID_inturn,1.);
      }else{
        adj_nuc[st][sub_step].copy(adj_nuc[st+1][0]);
      };
      for(int j=sub_step-1;j>=0;j--){

        // Adjoint nuclide density
	//GroupData2D mmat1=trmat_flxdep[st]*(vol_flx[st].get_sum()/fuel_vol*power_factor[st][j]*1e-24);
	GroupData2D mmat1=trmat_flxdep[st]*total_flux[st][j]*1e-24;
	GroupData2D mmat2=trmat_flxindep+mmat1;
        mmat2=mmat2.GetTransposedMatrix();
        GroupData2D eee=bu.CalMatrixExponentialByPade(mmat2,delt[st][j]);
        //GroupData2D eee=bu.CalMatrixExponentialByKrylov(mmat2,delt[st][j]);
        adj_nuc[st][j]=eee*adj_nuc[st][j+1];
        // Adjoint power calculation
        pow_adj[st][j]=adj_nuc[st][j]*(mmat1*fwd_nuc[st][j])*delt[st][j]/power_density;
        // Generalized adjoint flux calculation
        // (source calculation)
        GroupData1D gpt_src(group); // positive
        GroupData1D gpt_src2(group); // negative
        for(int g=0;g<group;g++){
	  real tmp=0.;
	  // (power term)
	  for(int k=0;k<nucn;k++){
  	    tmp+=mic_sigf[st][k].get_dat(g)*fwd_nuc[st][j].get_dat(k)
                *bu.GetReactionEnergyData().GetFissionEnergy(matno[k]);
	  };
	  tmp*=pow_adj[st][j];
  	  // (number density term)
          tmp-=(adj_nuc[st][j]*dmdf_nuc[st][j][g])*delt[st][j]/fuel_vol;
	  if(tmp>0.){
            gpt_src.put_data(g,tmp);
	  }else{
            gpt_src2.put_data(g,-tmp);
	  };
        };
	// Orthgolality check
	/*
	real sum=0.;
	real sump=0.;
	real sumn=0.;
        for(int g=0;g<group;g++){
	  sum+=vol_flx[st].get_dat(g)*(gpt_src.get_dat(g)-gpt_src2.get_dat(g));
	  sump+=vol_flx[st].get_dat(g)*gpt_src.get_dat(g);
	  sumn+=vol_flx[st].get_dat(g)*gpt_src2.get_dat(g);
	};
	cout<<" Phi * Q : "<<sum<<" "<<sump<<" "<<sumn<<"\n";
	*/

	for(int rr=0;rr<mesh_fuel;rr++){
          gpt_flx[st][j][rr].set_zero();
	};
	if(power_density>1e-5){

	  vector<GroupData1D> gpt_flx_old(mesh_fuel);
	  for(int m=0;m<mesh_fuel;m++){
	    gpt_flx_old[m].put_imax(group);
	  };

        // Initial iteration
	  /*
        PJISystem lat_gpt(group,3);
        lat_gpt.PutTrajectorySet(&sys_f);
        pin[0].GetMacxs().DataCopy(macxs[st][j]);
        lat_gpt.AddMedium(pin[0]);
        lat_gpt.AddMedium(pin[1]);
        lat_gpt.AddMedium(pin[2]);
        lat_gpt.PutRegMed(region_medium);
        lat_gpt.PutGeneralOption(opt);
        lat_gpt.PutSigmaCol(sigtr);
        lat_gpt.PutBuckling(0.);
        lat_gpt.WriteProcOff();
        lat_gpt.NoPrint();
	  */

	lat_gpt.GetMed(0).GetMacxs().DataCopy(macxs[st][j]);
        lat_gpt.PutPij();

	real keff_step=lat_gpt.CalIgenPij();
	lat_gpt.PutGeneralOption(opta);

	// (positive source)
        for(int g=0;g<group;g++){
          lat_gpt.SetZeroScatSrc(g);
        };
	for(int rr=0;rr<mesh_fuel;rr++){
          lat_gpt.PutIsotropicSourcePerUnitVolume(rr,gpt_src);
	};
        lat_gpt.CalFixedSourceUpScat();

        for(int m=0;m<mesh_fuel;m++){
          for(int g=0;g<group;g++){
            real tmp=lat_gpt.GetMesh(m).GetFlux().get_dat(g);
            gpt_flx[st][j][m].add_data(g,tmp);
            gpt_flx_old[m].put_data(g,tmp);
          };
        };

	int gptiter=1;
        for(int k=0;k<gptitermax;k++){
          for(int g=0;g<group;g++){
            lat_gpt.SetZeroScatSrc(g);
          };
          for(int m=0;m<mesh_fuel;m++){
            lat_gpt.GetMesh(m).CalFissionSrcAdjoint();
	  };
          lat_gpt.CalFixedFissionSourceUpScat(keff_step);
          //lat_gpt.CalFixedFissionSourceUpScat(1.);
	  gptiter++;
	  real errmax=0.;
          for(int m=0;m<mesh_fuel;m++){
            for(int g=0;g<group;g++){
              real tmp=lat_gpt.GetMesh(m).GetFlux().get_dat(g);
              real err=fabs(tmp/gpt_flx_old[m].get_dat(g)-1.);
	      if(err>errmax)errmax=err;
              gpt_flx[st][j][m].add_data(g,tmp);
              gpt_flx_old[m].put_data(g,tmp);
   	    };
	  };
	  if(errmax<1e-4)break;
        };
	/*
        for(int m=0;m<mesh_fuel;m++){
          for(int k=0;k<=gptitermax;k++){
            gpt_flx[st][j][m]=gpt_flx[st][j][m]-lat_gpt.GetMesh(m).GetFlux();
          };
        };
	*/
        for(int m=0;m<mesh_fuel;m++){
          gpt_flx[st][j][m]=gpt_flx[st][j][m]-(lat_gpt.GetMesh(m).GetFlux()*gptiter);
          //gpt_flx[st][j][m]=gpt_flx[st][j][m]
	  //-(lat_gpt.GetMesh(m).GetFlux()*pow((1-keff_tmp),gptiter)/(1-keff_tmp));
          //gpt_flx[st][j][m]=gpt_flx[st][j][m]-(lat_gpt.GetMesh(m).GetFlux()*factor);
        };

	// (negative source)
        for(int g=0;g<group;g++){
          lat_gpt.SetZeroScatSrc(g);
        };
	for(int rr=0;rr<mesh_fuel;rr++){
          lat_gpt.PutIsotropicSourcePerUnitVolume(rr,gpt_src2);
	};
        lat_gpt.CalFixedSourceUpScat();

        for(int m=0;m<mesh_fuel;m++){
          for(int g=0;g<group;g++){
            real tmp=lat_gpt.GetMesh(m).GetFlux().get_dat(g);
            gpt_flx[st][j][m].add_data(g,-tmp);
            gpt_flx_old[m].put_data(g,tmp);
          };
        };

	gptiter=1;
        for(int k=0;k<gptitermax;k++){
          for(int g=0;g<group;g++){
            lat_gpt.SetZeroScatSrc(g);
          };
          for(int m=0;m<mesh_fuel;m++){
            lat_gpt.GetMesh(m).CalFissionSrcAdjoint();
	  };
          //lat_gpt.CalFixedFissionSourceUpScat(keff[st]);
          lat_gpt.CalFixedFissionSourceUpScat(keff_step);
	  gptiter++;

	  real errmax=0.;
          for(int m=0;m<mesh_fuel;m++){
            for(int g=0;g<group;g++){
              real tmp=lat_gpt.GetMesh(m).GetFlux().get_dat(g);
              real err=fabs(tmp/gpt_flx_old[m].get_dat(g)-1.);
	      if(err>errmax)errmax=err;
              gpt_flx[st][j][m].add_data(g,-tmp);
              gpt_flx_old[m].put_data(g,tmp);
   	    };
	  };
          //cout<<"Neg.. "<<errmax<<"\n";
	  if(errmax<1e-4)break;
        };
	/*
        for(int m=0;m<mesh_fuel;m++){
          for(int k=0;k<=gptitermax;k++){
            gpt_flx[st][j][m]=gpt_flx[st][j][m]+lat_gpt.GetMesh(m).GetFlux();
          };
        };
	*/
        for(int m=0;m<mesh_fuel;m++){
	  //gpt_flx[st][j][m]=gpt_flx[st][j][m]+(lat_gpt.GetMesh(m).GetFlux()*(gptitermax+1));
          gpt_flx[st][j][m]=gpt_flx[st][j][m]+(lat_gpt.GetMesh(m).GetFlux()*gptiter);
          //gpt_flx[st][j][m]=gpt_flx[st][j][m]
	  //  -(lat_gpt.GetMesh(m).GetFlux()*pow((1-keff_tmp),gptiter)/(1-keff_tmp));
          //gpt_flx[st][j][m]=gpt_flx[st][j][m]-(lat_gpt.GetMesh(m).GetFlux()*factor);
        };

	};

        // Jump condition
        for(int k=0;k<nucn;k++){
          real tmp1=(mic_sigf[st][k]*(vol_flx[st]*power_factor[st][j]))
                   *bu.GetReactionEnergyData().GetFissionEnergy(matno[k]);
	  tmp1*=pow_adj[st][j];
          real tmp2=0.;
	  for(int g=0;g<group;g++){
	    real xsc=mic_sigc[st][k].get_dat(g); 
            real xsf=mic_sigf[st][k].get_dat(g);
            for(int m=0;m<3;m++){
	      tmp2+=fwd_vol_flx[st][m].get_dat(g)*power_factor[st][j]
		*gpt_flx[st][j][m].get_dat(g)*(xsc+xsf);
	    };
	  };
          adj_nuc[st][j].add_data(k,-tmp1+tmp2);
        };
      };
    };

    // +++ adjoint number density at t=0
    /*
    cout<<"\n\n+ adjoint number density at t=0 \n";
    for(int i=0;i<nucn;i++){
      real end_nuc=pin[0].GetNuclideInTurn(iii).GetDensity();
      real ini_nuc=den0[i];
      real adjn=adj_nuc[0][0].get_dat(i);
      //int matnum=pin[0].GetNuclideInTurn(i).GetMatnum();
      cout<<matno[i]<<" "<<adjn<<" "<<adjn/end_nuc<<" "<<adjn/end_nuc*ini_nuc<<"\n;
    };
    cout<<"\n\n";
    */

    // +++ sensitivity calculation
    cout<<"     "<<matno[target_nuclideID_inturn]<<"\n";// target-nuclide ID

    vector< vector< vector<real> > > gptadj_chi(burn_step);
    for(int i=0;i<burn_step;i++){
      int sub_step=sub_step_list[i];
      gptadj_chi[i].resize(sub_step);
      for(int j=0;j<sub_step;j++){
	gptadj_chi[i][j].resize(mesh_fuel);
	for(int m=0;m<mesh_fuel;m++){
	  real sum=0.;
	  for(int g=0;g<group;g++){
	    sum+=gpt_flx[i][j][m].get_dat(g)*macxs[i][j].GetData1d(chi).get_dat(g);
	  };
	  gptadj_chi[i][j][m]=sum;
	};
      };
    };

    real end_nuc=pin[0].GetNuclideInTurn(target_nuclideID_inturn).GetDensity();
    for(int i=0;i<nucn;i++){
      int matnum=pin[0].GetNuclideInTurn(i).GetMatnum();
      int rmax=2;
      if(matnum<9000)rmax=1;
      for(int r=0;r<rmax;r++){
	cout<<" 3\n"; // cross section indicator
	cout<<" "<<group<<"\n";
        cout<<" "<<matno[i]<<"\n";
        enum xstype sigxx=sigc;
	if(r==0){
	  cout<<" 102\n";
	}else{
	  sigxx=sigf;
	  cout<<" 18\n";
	};
	// (pre-calculation for number density term)
        vector< vector<real> > nadj_dm_nfwd(burn_step);
	for(int k=0;k<burn_step;k++){
	  int sub_step=sub_step_list[k];
	  nadj_dm_nfwd[k].resize(sub_step);
	  for(int l=0;l<sub_step;l++){
            real val=0.;
            // absorption
            val+=-fwd_nuc[k][l].get_dat(i)*adj_nuc[k][l].get_dat(i);
            if(sigxx==sigc){
              int tmp=bu.GetBC().GetNdivCapture(matnum);
              for(int j=0;j<tmp;j++){
                int id2=bu.GetBC().GetNextIDCapture(matnum,j);
                int pos=bu.SearchNuclide(id2);
                if(pos!=-1){
                  real rat=bu.GetBC().GetRatioCapture(matnum,j);
                  val+=adj_nuc[k][l].get_dat(pos)*rat*fwd_nuc[k][l].get_dat(i);
                };
              };
	    };
            // fission
            if(sigxx==sigf){
              int tmp2=bu.GetBC().GetNdivFission(matnum);
              for(int j=0;j<tmp2;j++){
                int id2=bu.GetBC().GetNextIDFission(matnum,j);
                int pos=bu.SearchNuclide(id2);
                if(pos!=-1){
                  real rat=bu.GetBC().GetRatioFission(matnum,j);
                  val+=adj_nuc[k][l].get_dat(pos)*rat*fwd_nuc[k][l].get_dat(i);
	        };
  	      };
	    };
  	    nadj_dm_nfwd[k][l]=val*delt[k][l];
	  };
	};

        real tot=0.;
        for(int j=0;j<group;j++){
          real sum=0.;
          for(int k=0;k<burn_step;k++){
            int sub_step=sub_step_list[k];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc[k][l].get_dat(i);
  	      // --- Number density term
              real org=0.;
              if(sigxx==sigc)org=mic_sigc[k][i].get_dat(j);
	      if(sigxx==sigf)org=mic_sigf[k][i].get_dat(j);
	      // (dM)*(flx*1e-24) = (dsig_j*phi_j/flx)*(flx*1e-24) = dsig_j*phi_j*1e-24
              real dsig=org*(vol_flx[k].get_dat(j)*power_factor[k][l]/fuel_vol)*1e-24;
 	      sum+=dsig*nadj_dm_nfwd[k][l];
  	      // --- Power normalization term (fission case)
              if(sigxx==sigf){
 	        real flx=vol_flx[k].get_dat(j)*power_factor[k][l];
                sum-=pow_adj[k][l]*flx*org*den*bu.GetReactionEnergyData().GetFissionEnergy(matno[i]);
	      };
	      // --- flux term
              for(int m=0;m<3;m++){
                real tmp2=fwd_vol_flx[k][m].get_dat(j)*power_factor[k][l];
  	        // (absorption)
                sum+=den*org*tmp2*gpt_flx[k][l][m].get_dat(j);
  	        // (yield)(only for fission case)
                if(sigxx==sigf){
                  real temp=den*org*tmp2*pin[0].GetNuclideInTurn(i).GetMicxs().GetData1d(nu).get_dat(j);
		  /*
    	          for(int g=0;g<group;g++){
	            sum-=temp*gpt_flx[k][l][m].get_dat(g)*macxs[k][l].GetData1d(chi).get_dat(g);
  	          };
		  */
		  sum-=temp*gptadj_chi[k][l][m];
		};
	      };
	    };
          };
          cout.setf(ios::scientific);
          cout.precision(6);
          cout<<"   "<<sum/end_nuc<<"\n";// (group-wise relative sensitivity output)
          tot+=sum;
	};
        //cout.setf(ios::scientific);
        //cout.precision(6);
        //cout<<"    "<<tot/end_nuc<<"\n";
      };
    };

    if(snsyld_cal){
      //cout<<"Sensitivity to fission yield\n";
    // +++ sensitivity calculation
    real end_nuc=pin[0].GetNuclideInTurn(target_nuclideID_inturn).GetDensity();
    int idfisn=3;
    int idfisorg[]={9228,9437,9443};
    for(int ii=0;ii<idfisn;ii++){
      int idfis=idfisorg[ii];
      //cout<<"Fission nuclide : "<<idfis<<"\n";
      cout<<" 18\n"; // yield sensitivity indicator
      cout<<" "<<idfis<<"\n"; // fission nuclide ID
      int pos0=bu.SearchNuclide(idfis);
      int nuct=bu.GetBC().GetNdivFission(idfis);
      cout<<" "<<nuct<<"\n";
      for(int i=0;i<nuct;i++){
        int id=bu.GetBC().GetNextIDFission(idfis,i);
        real rat=bu.GetBC().GetRatioFission(idfis,i);
        int pos=bu.SearchNuclide(id);
	if(pos!=-1){
	  //for(int tt=0;tt<nucn;tt++){
	  //  if(id==matno[tt])cout<<filename[tt]<<"   ";
	  //};
          real val=0.;
          for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
            for(int l=0;l<sub_step;l++){
              real sf=mic_sigf[k][pos0]*vol_flx[k]/vol_flx[k].get_sum();
	      sf*=total_flux[k][l]*1e-24;
              val+=(adj_nuc[k][l].get_dat(pos)*sf*rat*fwd_nuc[k][l].get_dat(pos0))*delt[k][l];
	    };
	  };
	  cout<<" "<<id<<"\n";
          cout.setf(ios::scientific);
          cout.precision(6);
          cout<<" "<<val/end_nuc<<"\n";
	};
      };
    };
    };

    if(snsdcy_cal){
      //cout<<"\n";
      //cout<<"Sensitivity to half-life\n";
      //cout<<"\n";
    // +++ sensitivity calculation
    real end_nuc=pin[0].GetNuclideInTurn(target_nuclideID_inturn).GetDensity();
    cout<<" 88\n"; // half-life sensitivity indicator
    for(int i=0;i<nucn;i++){
      int matnum=pin[0].GetNuclideInTurn(i).GetMatnum();
      real decay_c=bu.GetDecayConstant(matnum);
      if(decay_c!=0.){
        cout<<" "<<matno[i]<<"\n";
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
	sum*=100.;
        cout.setf(ios::scientific);
        cout.precision(6);
        cout<<" "<<sum/end_nuc<<"\n";
      };
    };
    cout<<" -1\n"; // half-life sensitivity end
    };
    cout<<" -1\n"; //
    }; // loop-end for nuclide
    cout<<" -1\n"; // 
  };


  cout<<"#+++ Number density after burnup +++\n";
  for(int i=0;i<nucn;i++){
    cout.setf(ios::scientific);
    cout.precision(5);
    real dens=pin[0].GetNuclideInTurn(i).GetDensity();
    if(dens!=0.)cout<<"# "<<filename[i]<<" "<<dens<<"\n";
  };
};
