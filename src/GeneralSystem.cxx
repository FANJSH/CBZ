#include <cstdlib>
#include "GeneralSystem.h"

using namespace std;

GeneralSystem::~GeneralSystem()
{
  vector<GeneralMesh>().swap(mesh);
  vector<Medium>().swap(med);

};


void GeneralSystem::Init(int n, int g, int i)
{
  PutNdim(n);
  grp = g;
  MAX_MED = i;
  nmed = 0;
  med.resize(MAX_MED);
  cmfd=false;
  print=true;
  IsUpScatterSystem=false;
  ThermalIteration=1;
  cmfd_factor=1.;
  total_inner_iteration.resize(grp,0);
};

void GeneralSystem::AddMedium(Medium &inp)
{
  if(inp.GetImax()!=grp){
    cout<<"# Different energy group!\n";
    cout<<"# between system and Medium.\n";
    cout<<"#  group of medium : "<<inp.GetImax()<<"\n";
    cout<<"#  group of system : "<<grp<<"\n";
    cout<<"#  name of system  : "<<name<<"\n";
    exit(0);
  };

  if(nmed==MAX_MED){
    cout<<"# Error in GeneralSystem::AddMedium.\n";
    cout<<"# You cannot put a instance of medium more.\n";
    cout<<"# Please increase MAX_MED in constructor of GeneralSystem.\n";
    cout<<"# [nmed] is "<<nmed<<" and [MAX_MED] is "<<MAX_MED<<"\n";
    exit(0);
  };

  med[nmed]=inp;
  med[nmed].PutLowestDownScatGrp();
  med[nmed].PutHighestUpScatGrp();
  med[nmed].GetMacxs().GetData1d(sigst).copy(med[nmed].GetMacxs().GetData2d(sigs).get_sumx());

  for(int i=0;i<grp;i++){
    if(inp.GetMacxs().GetData1d(nusigf).get_dat(i)<0.){
      cout<<"# !! Warning !!\n";
      cout<<"# Negative fission yield cross section is detected.\n";
      cout<<"# This cross section is set to zero.\n";
      cout<<"# Medium ID is "<<nmed<<"\n";
      med[nmed].GetMacxs().GetData1d(nusigf).put_data(i,0.);
    };

    if(name!="PJI"&&inp.GetMacxs().GetData1d(d).get_dat(i)<0.){
      cout<<"# !! Warning !!\n";
      cout<<"# Negative diffusion coefficient is detected in group "<<i<<".\n";
      cout<<"# The value is "<<inp.GetMacxs().GetData1d(d).get_dat(i)<<"\n";
      cout<<"# Medium ID is "<<nmed<<"\n";
    };
  };

  nmed++;
}

void GeneralSystem::AddMedium(int id, Medium &inp)
{
  if(id>=nmed){
    cout<<"# Error in GeneralSystem::AddMedium.\n";
    exit(0);
  };

  med[id]=inp;

  med[id].PutLowestDownScatGrp();
  med[id].PutHighestUpScatGrp();
  med[id].GetMacxs().GetData1d(sigst).copy(med[id].GetMacxs().GetData2d(sigs).get_sumx());

  for(int i=0;i<grp;i++){
    if(inp.GetMacxs().GetData1d(nusigf).get_dat(i)<0.){
      cout<<"# !! Warning !!\n";
      cout<<"# Negative fission yield cross section is detected.\n";
      cout<<"# This cross section is set to zero.\n";
      cout<<"# Medium ID is "<<nmed<<"\n";
      med[id].GetMacxs().GetData1d(nusigf).put_data(i,0.);
    };

    if(name!="PJI"&&inp.GetMacxs().GetData1d(d).get_dat(i)<0.){
      cout<<"# !! Warning !!\n";
      cout<<"# Negative diffusion coefficient is detected in group "<<i<<".\n";
      cout<<"# The value is "<<inp.GetMacxs().GetData1d(d).get_dat(i)<<"\n";
      cout<<"# Medium ID is "<<id<<"\n";
    };
  };
}

void GeneralSystem::AddMedium(string mdir,string ss,int plt,bool upscat)
{
  Medium minp;
  minp.ReadPDSFile(mdir,ss,plt,upscat);
  AddMedium(minp);
};

void GeneralSystem::UnifyFissionSpectrum(int i)
{
  if(i<0||i>=nmed){
    cout<<"Error in UnifyFissionSpectrum of GeneralSystem.\n";
    cout<<"Please check i value.\n";
    exit(0);
  };
  GroupData1D kai=med[i].GetMacxs().GetKai();
  for(int j=0;j<nmed;j++){
    med[j].GetMacxs().GetKai().copy(kai);
  };
};

void GeneralSystem::EdgeCalculation()
{
  vector<int> num_fmesh(3,0);
  for(int i=0;i<3;i++){
    int tmp=1;
    if(i<=Ndim){
      tmp=mi.GetFMesh(i);
    };
    num_fmesh[i]=tmp;
  };

  for(int i=0;i<3;i++){

    int axis3=i;
    int axis1=1;
    int axis2=2;
    if(i==1)axis1=0; 
    if(i==2){ 
      axis1=0;
      axis2=1;
    };

    int fic_x=num_fmesh[axis1];
    int fic_y=num_fmesh[axis2];
    int fic_z=num_fmesh[axis3];
    edge[i].resize(fic_x);
    for(int j=0;j<fic_x;j++){
      edge[i][j].resize(fic_y);
      for(int k=0;k<fic_y;k++){
        edge[i][j][k].resize(2,-1);
      };
    };

    for(int j=0;j<fic_x;j++){
      for(int k=0;k<fic_y;k++){

        bool ex1=false;
	bool ex2=false;
	for(int l=0;l<fic_z;l++){
	  int real_x,real_y,real_z;
	  if(i==0){
	    real_x=l;
	    real_y=j;
	    real_z=k;
	  }else if(i==1){
	    real_x=j;
	    real_y=l;
	    real_z=k;
	  }else{
	    real_x=j;
	    real_y=k;
	    real_z=l;
	  };

	  int mid=meshid[real_z][real_y][real_x];
          if(mid!=-1&&!ex1){
            edge[i][j][k][0]=l;
	    ex1=true;
	  };
	  if(mid==-1&&ex1&&!ex2){
	    edge[i][j][k][1]=l-1;
	    ex2=true;
	  };
	  if(mid!=-1&&ex2){
	    if(name!="DHEX"||i!=1){
  	      cout<<"Error in GeneralSystem::EdgeCalculation.\n";
	      cout<<"  position   x : "<<real_x<<"\n";
              cout<<"             y : "<<real_y<<"\n";
              cout<<"             z : "<<real_z<<"\n";
              cout<<"  axis         : "<<i<<"\n";
	      exit(0);
	    };
	  };
	};
	if(ex1&&!ex2){
          edge[i][j][k][1]=fic_z-1;
        };
	if(edge[i][j][k][0]==-1){
	  edge[i][j][k][1]=-2;
	};
      };
    };
  };

};

void GeneralSystem::SetCoarseMeshInfo()
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  cmeshid.resize(xr);
  for(int i=0;i<xr;i++){
    cmeshid[i].resize(yr);
    for(int j=0;j<yr;j++){
      cmeshid[i][j].resize(zr,0);
    };
  };

  int index=0;
  int index2=0;
  for(int z1=0;z1<zr;z1++){
    for(int y1=0;y1<yr;y1++){
      for(int x1=0;x1<xr;x1++){
        int id=mi.GetCMat(index);
        if(id!=-1){
	  cmeshid[x1][y1][z1]=index2;
	  index2++;
	}else{
	  cmeshid[x1][y1][z1]=-1;
	};
	index++;
      };
    };
  };
  CTotM=index2;
};

void GeneralSystem::PutGeneralOption(GeneralOption inp)
{
  opt = inp;
  opt.PutGrp(grp);
};

void GeneralSystem::PutNdim(int i)
{
  Ndim=i;
  edge.resize(3);
};

void GeneralSystem::PutPL(int p)
{
  pl = p;
  if(Ndim==1)plnum=pl+1;
  if(Ndim==2||Ndim==3){
    int ii=0;
    for(int i=0;i<=pl;i++){
      if(Ndim==2)ii+=i+1; 
      if(Ndim==3)ii+=i*2+1;
    };
    plnum=ii;
  };
};

void GeneralSystem::WriteSolverName(string acname)
{
  if(print){
  cout<<"#*************************************************\n";
  cout<<"#* System CBG                                     \n";
  cout<<"#*   Power iteration                              \n";
  cout<<"#*   Solver : "<<name<<"\n";
  cout<<"#*   Acceleration : "<<acname<<"\n";
  if(opt.Forward()){
    cout<<"#*   Forward calculation\n";}
  else{
    cout<<"#*   Adjoint calculation\n";
  };
  cout<<"#*   Convergence condition (k_eff) : "<<opt.GetEpsk()<<"\n";
  cout<<"#*                         (flux)  : "<<opt.GetEpsf()<<"\n";
  cout<<"#*                         (source): "<<opt.GetEpss()<<"\n";
  cout<<"#*************************************************\n";
  cout<<"#*It:   K_eff   : Err in k : Err in   : Err in    \n";
  cout<<"#*er:           :          : flux     : source    \n";
  cout<<"#*************************************************\n";  
  };
};

void GeneralSystem::WritePerturbName()
{
  cout<<"#*************************************************\n";
  cout<<"#* System CBG                                     \n";
  cout<<"#*   Perturbation calculation\n";
  cout<<"#*   Solver : "<<name<<"\n";
  cout<<"#*************************************************\n";  
};


// **********************************************
// * For power iteration                        *
// **********************************************

void GeneralSystem::CheckUpScattering()
{
  int ussrc=grp;
  int ussink=grp;
  for(int i=0;i<nmed;i++){
    int tmp2=med[i].GetUpScatteringGrp();
    if(tmp2!=-1&&tmp2<ussrc)ussrc=tmp2;
    int tmp3=med[i].GetUpScatteringSinkGrp();
    if(tmp3!=-1&&tmp3<ussink)ussink=tmp3;
  };
  UpScatSrcGrp=ussrc;
  UpScatSinkGrp=ussink;

  if(UpScatSrcGrp!=grp)IsUpScatterSystem=true;
  UpScatteringGrp=UpScatSrcGrp;
  if(!opt.Forward())UpScatteringGrp=UpScatSinkGrp;

  if(IsUpScatterSystem){
    int tmp=grp-UpScatteringGrp;
    for(int i=0;i<TotM;i++){
      mesh[i].InitializeUpScatSrc(tmp);
    };
  };
};

real GeneralSystem::CalIgen(string accel,real maxval)
{
  real keff=0;
  if(accel=="cmfd"&&cmfdimp){
    cmfd=true;
    //if(!opt.Forward()){
    //  cout<<"Adjoint calculation is performed without CMFD acceleration.\n";
    //  cmfd=false;
    //};
    keff=CalIgenNormal(pl);
  }else if(accel=="extrapolation"){
    keff=CalIgenExtrapolation(pl);
  }else if(accel=="chebyshev"){
    keff=CalIgenChebychev(pl,maxval);
  }else if(accel=="chebychev"){
    keff=CalIgenChebychev(pl,maxval);
  }else{
    cmfd=false;
    keff=CalIgenNormal(pl);
  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };

  return keff;
};

real GeneralSystem::CalIgenWithFissionSpectrumMatrix()
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  SetInitialFlux();
  NormalizeFissionSrc();
  CheckUpScattering();
  
  real cin=0.1;
  fiss=1.;

  string acname="none"; // o
  if(cmfd)acname="CMFD"; // o
  WriteSolverName(acname);

  // Group-dependent fission source
  vector< vector<real> > group_fiss(TotM);
  for(int i=0;i<TotM;i++){
    group_fiss[i].resize(grp,0.);
  };

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    //real tmp=2e-5;
    //if(opt.GetEpsf()<2e-5)tmp=opt.GetEpsf()*0.5;
    //if(cin<tmp)cin=tmp;

    // Group-dependent fission source calculation
    real invsrc=1./(fiss*PI4);
    for(int i=0;i<TotM;i++){
      for(int j=0;j<grp;j++){
	group_fiss[i][j]=0.;
      };
      real vol=mesh[i].GetVolume();
      if(opt.Forward()){
        for(int j=0;j<grp;j++){
	  real srcj=mesh[i].GetMed()->GetMacxs().GetData1d(nusigf).get_dat(j)
	           *mesh[i].GetFlux().get_dat(j)*vol*invsrc;
	  for(int k=0;k<grp;k++){
	    group_fiss[i][k]+=mesh[i].GetMed()->GetMacxs().GetData2d(chi).get_dat(j,k)*srcj;
	  };
	};
      }else{
	for(int j=0;j<grp;j++){
	  real tmp=0.;
	  for(int k=0;k<grp;k++){
	    tmp+=mesh[i].GetMed()->GetMacxs().GetData2d(chi).get_dat(j,k)
	      *mesh[i].GetFlux().get_dat(k);
	  };
          group_fiss[i][j]+=mesh[i].GetMed()->GetMacxs().GetData1d(nusigf).get_dat(j)
   	                   *tmp*vol*invsrc;
	};
      };
    };

    errf=0.;
    bool InnerConvergence=true;

    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;

      for(int i=0;i<TotM;i++){
	int index=0;
	for(int l=0;l<=pl;l++){
	  int is=0;
	  int ie=0;
	  if(Ndim!=1)ie=l;
	  if(Ndim==3)is=-l;
	  for(int m=is;m<=ie;m++){
	    real tmp=mesh[i].GetScatSrc(ginp,index);
	    if(IsUpScatterSystem)tmp+=mesh[i].GetUpScatSrc(ginp,index);
	    if(index==0)tmp+=group_fiss[i][ginp];
	    mesh[i].PutSrcin(tmp,index);
	    index++;
	  };
	};
      };
      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
      if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
      if(IsUpScatterSystem&&UpScatteringGrp<=ginp)AddUpScatSrc(ginp,pl);
    };

    // +++ Thermal iteration
    real sor_factor=1.2;
    if(IsUpScatterSystem){
      for(int ii=0;ii<ThermalIteration*3;ii++){
        real errf=0.;
        for(int g=0;g<grp;g++){
          int ginp=g;
          if(!opt.Forward())ginp=grp-1-g;
          if(UpScatteringGrp<=ginp){ 

            for(int i=0;i<TotM;i++){
	      int index=0;
   	      for(int l=0;l<=pl;l++){
	        int is=0;
	        int ie=0;
	        if(Ndim!=1)ie=l;
	        if(Ndim==3)is=-l;
	        for(int m=is;m<=ie;m++){
	          real tmp=mesh[i].GetScatSrc(ginp,index);
	          if(IsUpScatterSystem)tmp+=mesh[i].GetUpScatSrc(ginp,index);
	          if(index==0)tmp+=group_fiss[i][ginp];
	          mesh[i].PutSrcin(tmp,index);
	          index++;
	        };
	      };
            };
            vector<real> flxold(TotM);
	    for(int i=0;i<TotM;i++){
	      flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	    };
  	    real err=CalFluxGeneral(ginp,cin,1);
	    //             (important)-------^
	    for(int i=0;i<TotM;i++){
	      real newflx=mesh[i].GetFlux().get_dat(g);
	      mesh[i].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
  	    };
	    if(err>errf)errf=err;

          };
          AddDownScatSrc(ginp,pl);
          SetZeroScatSrc(ginp);
          if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
          if(IsUpScatterSystem&&UpScatteringGrp<=ginp)AddUpScatSrc(ginp,pl);
	};
        if(print)cout<<"#  (Thermal iteration) "<<ii<<" : "<<errf<<"\n";
      };
    };


    if(cmfd){
      DoAcceleration(iter,errs,fiss); // o
      // Correction for up-scattering
      if(IsUpScatterSystem){
        for(int g=0;g<grp;g++){
          int ginp=g;
          if(!opt.Forward())ginp=grp-1-g;
          SetZeroScatSrc(ginp);
          if(UpScatteringGrp<=ginp)AddUpScatSrc(ginp,pl);
	};
      };
    };


    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real tmp2=mesh[i].GetFissionSrc();
	real err=fabs(tmp2/tmp-1.0);
	if(err>errs)errs=err;
      };
    };

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){break;}
      else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };
  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };


  return fiss;
}

real GeneralSystem::CalIgenGroupParallel()
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  SetInitialFlux();
  NormalizeFissionSrc();
  CheckUpScattering();

  
  real cin=0.1;
  fiss=1.;

  string acname="none"; // o
  if(cmfd)acname="CMFD"; // o
  WriteSolverName(acname);

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    real tmp=2e-5;
    if(opt.GetEpsf()<2e-5)tmp=opt.GetEpsf()*0.5;
    if(cin<tmp)cin=tmp;

    errf=0.;
    bool InnerConvergence=true;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,fiss,pl);
      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      //AddDownScatSrc(ginp,pl);
      //SetZeroScatSrc(ginp);
    };

    for(int g=0;g<grp;g++){
      SetZeroScatSrc(g);
    };
    for(int g=0;g<grp;g++){
      AddDownScatSrc(g,pl);
    };

    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real tmp2=mesh[i].GetFissionSrc();
	real err=fabs(tmp2/tmp-1.0);
	if(err>errs)errs=err;
      };
    };

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){break;}
      else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };
  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };

  return fiss;
}

real GeneralSystem::CalIgenNormal(int pl)
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  SetInitialFlux();
  NormalizeFissionSrc();
  CheckUpScattering();
  
  real cin=0.1;
  fiss=1.;

  bool DoThermalIT=false;
  if(IsUpScatterSystem&&ThermalIteration>0&&opt.Forward()){
    DoThermalIT=true;
  };

  vector< vector< vector<real> > > down_scat_src_store;
  if(DoThermalIT){
    int tmp=grp-UpScatSinkGrp;
    down_scat_src_store.resize(tmp);
    for(int i=0;i<tmp;i++){
      down_scat_src_store[i].resize(TotM);
      for(int j=0;j<TotM;j++){
	down_scat_src_store[i][j].resize(plnum,0.);
      };
    };
  };

  if(IsUpScatterSystem&&ThermalIteration>0&&print){
    cout<<"#\n";
    cout<<"# !! Caution !!\n";
    cout<<"# Thermal iteration is customized for 107-group.\n";
    cout<<"#\n";
  };

  string acname="none"; // o
  if(cmfd)acname="CMFD"; // o
  WriteSolverName(acname);

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    //real tmp=1e-5;
    //if(opt.GetEpsf()<1e-5)tmp=opt.GetEpsf()*0.5;
    //if(cin<tmp)cin=tmp;

    errf=0.;
    bool InnerConvergence=true;

    // Thermal iteration for adjoint
    if(iter!=0&&!opt.Forward()&&IsUpScatterSystem&&ThermalIteration>0){
      int tmp=ThermalIteration;
      //if(iter==0)tmp=3;
      for(int itt=0;itt<tmp*3;itt++){
      //for(int itt=0;itt<tmp;itt++){
      real errf=0.;

      for(int g=grp-1;g>=UpScatSinkGrp;g--){

        bool flxcal=false;
	if(itt%3==2){
	  flxcal=true;
	}else if(itt%3==0&&g>grp-1-20){
	  flxcal=true;
	}else if(itt%3==1&&g>grp-1-20){
	  flxcal=true;
	};
	//flxcal=true;

	if(flxcal){

        real sor_factor=1.2;
	CalSrcMultiplySystem(g,fiss,pl);
        vector<real> flxold(TotM);
	for(int i=0;i<TotM;i++){
	  flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	};
	real err=CalFluxGeneral(g,cin,iter);
	for(int i=0;i<TotM;i++){
	  real newflx=GetMesh(i).GetFlux().get_dat(g);
	  //GetMesh(i).GetFlux().put_data(g,flxold[i]+(newflx-flxold[i])*1.1);
	  GetMesh(i).GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
	};
	if(err>errf)errf=err;
	//AddDownScatSrc(g,pl);
	for(int i=0;i<TotM;i++){
	  mesh[i].AddDownScatSrcAdjointArv(g,UpScatSinkGrp,pl);
	};

	};

	SetZeroScatSrc(g);
        SetZeroUpScatSrc(g);
	AddUpScatSrc(g,pl);
      };

      if(print&&itt%3==2)cout<<"#  (Thermal iteration) "<<itt<<" : "<<errf<<"\n";
      //if(print)cout<<"#  (Thermal iteration) "<<itt<<" : "<<errf<<"\n";
      };
      //for(int g=0;g<UpScatSinkGrp;g++){
      //SetZeroScatSrc(g);
      //};
    };

    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;

      if(DoThermalIT&&ginp==UpScatSinkGrp){
        for(int i=UpScatSinkGrp;i<grp;i++){
	  for(int j=0;j<TotM;j++){
	    for(int k=0;k<plnum;k++){
	      down_scat_src_store[i-UpScatSinkGrp][j][k]=mesh[j].GetScatSrc(i,k);
	    };
	  };
	};
      };

      CalSrcMultiplySystem(ginp,fiss,pl);
      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
      if((opt.Forward()&&ginp>=UpScatSinkGrp)||
         (!opt.Forward()&&ginp>=UpScatSrcGrp)) SetZeroUpScatSrc(ginp);
      if((opt.Forward()&&ginp>=UpScatSrcGrp)||
         (!opt.Forward()&&ginp>=UpScatSinkGrp)) AddUpScatSrc(ginp,pl);
    };

    // Thermal iteration
    if(DoThermalIT){
    for(int i=0;i<ThermalIteration*3;i++){
      real errf=0.;

      int grptmp=UpScatSinkGrp;
      if(i%3!=0){
        if(grp==107)grptmp=grp-1-20;
	if(grp==21)grptmp=grp-1-6;
      };
      //for(int ii=UpScatSinkGrp;ii<grp;ii++){
      for(int ii=grptmp;ii<grp;ii++){
        for(int j=0;j<TotM;j++){
	  for(int k=0;k<plnum;k++){
	    mesh[j].PutScatSrc(ii,down_scat_src_store[ii-UpScatSinkGrp][j][k],k);
	  };
	};
      };

      for(int g=UpScatSinkGrp;g<grp;g++){
        int ginp=g;
        if(!opt.Forward())ginp=grp-1-g;
	bool flxcal=false;
	if(i%3==0){
	  flxcal=true;
	}else if(i%3==1&&g>grptmp){
	  flxcal=true;
	}else if(i%3==2&&g>grptmp){
	  flxcal=true;
	};

	if(flxcal){
          real sor_factor=1.2;
  	  CalSrcMultiplySystem(ginp,fiss,pl);
          vector<real> flxold(TotM);
	  for(int i=0;i<TotM;i++){
	    flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	  };
	  real err=CalFluxGeneral(ginp,cin,1);
	  //             (important)-------^
	  for(int i=0;i<TotM;i++){
	    real newflx=mesh[i].GetFlux().get_dat(g);
	    mesh[i].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
	  };
	  if(err>errf)errf=err;
	};
        AddDownScatSrc(ginp,pl);
        SetZeroScatSrc(ginp);
        if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
        AddUpScatSrc(ginp,pl);
      };
      //if(print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
      if(i%3==0&&print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
    };
    };

    if(cmfd){
      DoAcceleration(iter,errs,fiss); // o
      // Correction for up-scattering
      if(IsUpScatterSystem){
        //for(int g=UpScatSinkGrp;g<grp;g++){
        for(int g=0;g<grp;g++){
          int ginp=g;
          if(!opt.Forward())ginp=grp-1-g;
          if(IsUpScatterSystem){
            if((opt.Forward()&&ginp>=UpScatSinkGrp)||
               (!opt.Forward()&&ginp>=UpScatSrcGrp)) SetZeroUpScatSrc(ginp);
            if((opt.Forward()&&ginp>=UpScatSrcGrp)||
               (!opt.Forward()&&ginp>=UpScatSinkGrp)) AddUpScatSrc(ginp,pl);
	  };
	};
      };
    };


    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real tmp2=mesh[i].GetFissionSrc();
	/*
	cout.setf(ios::scientific);
	cout.precision(5);
	cout<<i<<" "<<tmp2-tmp<<"\n";
	*/
	real err=fabs(tmp2/tmp-1.0);
	if(err>errs)errs=err;
      };
    };
    //cout<<"\n\n";

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){
        out_iter_end=iter;
        break;
      }else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };
  };

  return fiss;
}

real GeneralSystem::CalIgenChebychev(int pl,real maxval)
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  SetInitialFlux();
  NormalizeFissionSrc();
  CheckUpScattering();


  real cin=0.1;
  fiss=1.;

  // *** For Chebyshev ***
  real *fsold=new real[TotM];
  for(int i=0;i<TotM;i++){fsold[i]=0.;};
  int itchev=10;
  real cha,chb,chc,ch1,ch2,alp,bet,gam;
  real chp[3], chq[3];
  cha=-maxval;   chb=maxval;
  chc=(2.-chb-cha)/(chb-cha);
  ch1=chc+sqrt(chc*chc-1);
  ch2=chc-sqrt(chc*chc-1);
  chp[0]=1.;  chq[0]=1.;
  chp[1]=ch1;  chq[1]=ch2;
  chp[2]=ch1*ch1;  chq[2]=ch2*ch2;
  // ********************

  WriteSolverName("Chebychev");

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    real tmp=2e-5;
    if(opt.GetEpsf()<2e-5)tmp=opt.GetEpsf()*0.5;
    if(cin<tmp)cin=tmp;

    errf=0.;
    bool InnerConvergence=true;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,fiss,pl);
      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
      if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
      if(UpScatteringGrp<=ginp)AddUpScatSrc(ginp,pl);
    };


    // *** For Chebychev ***
    if(iter==itchev){
      alp=2./(chc*(chb-cha));
      bet=-(cha+chb)/(chc*(chb-cha));
      gam=0.;
    }else if(iter>itchev){
      alp=4./(chb-cha)*(chp[1]+chq[1])/(chp[2]+chq[2]);
      bet=-2*(chb+cha)/(chb-cha)*(chp[1]+chq[1])/(chp[2]+chq[2]);
      gam=(chp[0]+chq[0])/(chp[2]+chq[2]);
    }else{
      alp=1.;   bet=0.;   gam=0.;
    };
    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real fsnew=mesh[i].GetFissionSrc();
	mesh[i].PutFissionSrc(alp*fsnew+bet*tmp-gam*fsold[i]);
	real tmp2=mesh[i].GetFissionSrc();
	real err=fabs(tmp2/tmp-1.);
	if(err>errs)errs=err;
	fsold[i]=tmp;
      };
    };
    if(iter>itchev-1){
      for(int i=0;i<3;i++){
	chp[i]*=ch1;
	chq[i]*=ch2;
      };
    };
    // ****
 
    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){break;}
      else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };
 };

  delete [] fsold; // *** for Chebyshev

  return fiss;
}

real GeneralSystem::CalIgenExtrapolation(int pl)
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  SetInitialFlux();
  NormalizeFissionSrc();
  CheckUpScattering();
  
  real cin=0.1;
  fiss=1.;

  bool DoThermalIT=false;
  if(IsUpScatterSystem&&ThermalIteration>0&&opt.Forward()){
    DoThermalIT=true;
  };

  if(IsUpScatterSystem&&ThermalIteration>0&&print){
    cout<<"#\n";
    cout<<"# !! Caution !!\n";
    cout<<"# Thermal iteration is customized for 107-group.\n";
    cout<<"#\n";
  };

  vector< vector< vector<real> > > down_scat_src_store;
  if(DoThermalIT){
    int tmp=grp-UpScatSinkGrp;
    down_scat_src_store.resize(tmp);
    for(int i=0;i<tmp;i++){
      down_scat_src_store[i].resize(TotM);
      for(int j=0;j<TotM;j++){
	down_scat_src_store[i][j].resize(plnum,0.);
      };
    };
  };

  // *** for use of extrapolation
  vector<real> fsold(TotM);
  vector<real> res(TotM);
  real lamda=0.;
  real lamdaold;
  int expit=0;
  bool exp=false;
  // ***

  WriteSolverName("Source extrapolation");

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    real tmp=2e-5;
    if(opt.GetEpsf()<2e-5)tmp=opt.GetEpsf()*0.5;
    if(cin<tmp)cin=tmp;

    errf=0.;
    bool InnerConvergence=true;

    // +++++ Thermal iteration for adjoint
    if(iter!=0&&!opt.Forward()&&IsUpScatterSystem&&ThermalIteration>0){
      int tmp=ThermalIteration;
      //if(iter==0)tmp=3;
      for(int itt=0;itt<tmp*3;itt++){
	//for(int itt=0;itt<tmp;itt++){
      real errf=0.;

      for(int g=grp-1;g>=UpScatSinkGrp;g--){

        bool flxcal=false;
	if(itt%3==2){
	  flxcal=true;
	}else if(itt%3==0&&g>grp-1-20){
	  flxcal=true;
	}else if(itt%3==1&&g>grp-1-20){
	  flxcal=true;
	};
	//flxcal=true;

	if(flxcal){

        real sor_factor=1.2;
	CalSrcMultiplySystem(g,fiss,pl);
        vector<real> flxold(TotM);
	for(int i=0;i<TotM;i++){
	  flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	};
	real err=CalFluxGeneral(g,cin,iter);
	for(int i=0;i<TotM;i++){
	  real newflx=GetMesh(i).GetFlux().get_dat(g);
	  GetMesh(i).GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
	};
	if(err>errf)errf=err;
	//AddDownScatSrc(g,pl);
	for(int i=0;i<TotM;i++){
	  mesh[i].AddDownScatSrcAdjointArv(g,UpScatSinkGrp,pl);
	};

	};

	SetZeroScatSrc(g);
        SetZeroUpScatSrc(g);
	AddUpScatSrc(g,pl);
      };

      if(print&&itt%3==2)cout<<"#  (Thermal iteration) "<<itt<<" : "<<errf<<"\n";
      //if(print)cout<<"#  (Thermal iteration) "<<itt<<" : "<<errf<<"\n";
      };
      //for(int g=0;g<UpScatSinkGrp;g++){
      //SetZeroScatSrc(g);
      //};
    }; // End of thermal iteration for adjoint +++++

    //int gmax=0;
    for(int g=0;g<grp;g++){

      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;

      if(DoThermalIT&&ginp==UpScatSinkGrp){
        for(int i=UpScatSinkGrp;i<grp;i++){
	  for(int j=0;j<TotM;j++){
	    for(int k=0;k<plnum;k++){
	      down_scat_src_store[i-UpScatSinkGrp][j][k]=mesh[j].GetScatSrc(i,k);
	    };
	  };
	};
      };

      CalSrcMultiplySystem(ginp,fiss,pl);

      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      //if(err>errf)gmax=ginp;
      if(err>errf)errf=err;
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
      if((opt.Forward()&&ginp>=UpScatSinkGrp)||
         (!opt.Forward()&&ginp>=UpScatSrcGrp)) SetZeroUpScatSrc(ginp);
      if((opt.Forward()&&ginp>=UpScatSrcGrp)||
         (!opt.Forward()&&ginp>=UpScatSinkGrp)) AddUpScatSrc(ginp,pl);
    };
    //cout<<gmax<<" !\n";

    // Thermal iteration
    if(DoThermalIT){
      //for(int i=0;i<ThermalIteration;i++){
    for(int i=0;i<ThermalIteration*3;i++){
      real errf=0.;

      int grptmp=UpScatSinkGrp;
      if(i%3!=0){
        if(grp==107)grptmp=grp-1-20;
	if(grp==21)grptmp=grp-1-6;
      };

      //for(int ii=UpScatSinkGrp;ii<grp;ii++){
      for(int ii=grptmp;ii<grp;ii++){
        for(int j=0;j<TotM;j++){
	  for(int k=0;k<plnum;k++){
	    mesh[j].PutScatSrc(ii,down_scat_src_store[ii-UpScatSinkGrp][j][k],k);
	  };
	};
      };

      for(int g=UpScatSinkGrp;g<grp;g++){
        int ginp=g;
        if(!opt.Forward())ginp=grp-1-g;
	bool flxcal=false;
	if(i%3==0){
	  flxcal=true;
	}else if(i%3==1&&g>grptmp){
	  flxcal=true;
	}else if(i%3==2&&g>grptmp){
	  flxcal=true;
	};
        //flxcal=true;

	if(flxcal){
          //real sor_factor=1.;
          real sor_factor=1.2;
  	  CalSrcMultiplySystem(ginp,fiss,pl);
          vector<real> flxold(TotM);
	  for(int i=0;i<TotM;i++){
	    flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	  };
	  real err=CalFluxGeneral(ginp,cin,1);
	  //             (important)-------^
	  for(int i=0;i<TotM;i++){
	    real newflx=mesh[i].GetFlux().get_dat(g);
	    mesh[i].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
	  };
	  if(err>errf)errf=err;
	};
        AddDownScatSrc(ginp,pl);
        SetZeroScatSrc(ginp);
        if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
        AddUpScatSrc(ginp,pl);
      };
      //if(print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
      if(i%3==0&&print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
    };
    };

    // +++ for use of extrapolation
    lamdaold=lamda;
    lamda=SourceAndResidualRevision(fsold,res);
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&errf<0.1&&lamda<1.0&&!exp&&iter>0){
      real omega=1./(1.-lamda);
      exp=true;
      if(print)WriteSourceExtrapolationInfo(lamda,omega);
      SourceExtrapolation(fsold,omega);
    }else if(exp){
      expit=0;
      exp=false;
    };
    errs=GetSourceError(fsold);
    // ****

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){break;}
      else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };
  };

  return fiss;
}

real GeneralSystem::GetFissionSum()
{
  real fiss=0.;
  for(int i=0;i<TotM;i++){
    fiss+=mesh[i].GetFissionSrc();
  };
  return fiss;
}

void GeneralSystem::WriteIterationInfo(int iter,real fiss,real errk,real errf,real errs)
{
  if(print){
  cout<<"# ";
  cout.width(3);
  cout<<iter<<":   ";
  cout.setf(ios::showpoint);
  cout.precision(6);
  cout<<fiss<<"  ";
  cout.unsetf(ios::showpoint);
  cout.setf(ios::scientific);
  cout.precision(3);
  cout<<errk<<"  "<<errf<<"  "<<errs<<"\n";
  cout.unsetf(ios::scientific);
  };
};

void GeneralSystem::CalSrcMultiplySystem(int g,real keff,int pl)
{
  bool ups=false;
  if(IsUpScatterSystem&&g>=UpScatSinkGrp)ups=true;

  real invsrc=1./(keff*PI4);
  for(int i=0;i<TotM;i++){
    real fiss;
    if(opt.Forward()){fiss=mesh[i].GetFissionSrcKai(g);}
    else{fiss=mesh[i].GetFissionSrcSigf(g);};
    fiss*=invsrc;
    int index=0;
    for(int l=0;l<=pl;l++){
      int is=0;
      int ie=0;
      if(Ndim!=1)ie=l;
      if(Ndim==3)is=-l;
      for(int m=is;m<=ie;m++){   
        real tmp=mesh[i].GetScatSrc(g,index);
	if(ups)tmp+=mesh[i].GetUpScatSrc(g,index);
	if(index==0)tmp+=fiss;
	mesh[i].PutSrcin(tmp,index);
	index++;
      };
    };
  };
};

void GeneralSystem::NormalizeFissionSrc()
{
  real fiss=0.;
  for(int i=0;i<TotM;i++){
    if(opt.Forward()){mesh[i].CalFissionSrc();}
    else{mesh[i].CalFissionSrcAdjoint();};
    fiss+=mesh[i].GetFissionSrc();
  };
  if(fiss==0.){
    cout<<"Error in GeneralSystem::NormalizeFissionSrc.\n";
    cout<<"Fission source is zero.\n";
    exit(0);
  };
  real inv_fiss=1./fiss;
  for(int i=0;i<TotM;i++){
    mesh[i].GetFlux()=mesh[i].GetFlux()*inv_fiss;
    if(opt.Forward()){mesh[i].CalFissionSrc();}
    else{mesh[i].CalFissionSrcAdjoint();};
  };
};


void GeneralSystem::SetInitialFlatFlux(real inp)
{
  for(int i=0;i<TotM;i++){
    for(int j=0;j<grp;j++){
      mesh[i].GetFlux().put_data(j,inp);
      /*
      if(i<TotM/2){
        mesh[i].GetFlux().put_data(j,inp);
      }else{
        mesh[i].GetFlux().put_data(j,1e-5);
      };
      */
    };
  };
};

void GeneralSystem::AddDownScatSrcReactionWise(int g,int pl,enum xstype xst,bool selfscat)
{
  for(int i=0;i<TotM;i++){
    mesh[i].AddDownScatSrcReactionWise(g,xst,pl,selfscat);
  };
};

void GeneralSystem::AddDownScatSrc(int g,int pl,bool selfscat)
{
  for(int i=0;i<TotM;i++){
    if(opt.Forward()){mesh[i].AddDownScatSrc(g,pl,selfscat);}
    else{mesh[i].AddDownScatSrcAdjoint(g,pl);};
  };
};

void GeneralSystem::AddUpScatSrc(int g,int pl)
{
  for(int i=0;i<TotM;i++){
    if(opt.Forward()){mesh[i].AddUpScatSrc(g,pl);}
    else{mesh[i].AddUpScatSrcAdjoint(g,pl);};
  };
};

void GeneralSystem::CalFissionSrcAllMesh()
{
  for(int i=0;i<TotM;i++){
    if(opt.Forward()){
      mesh[i].CalFissionSrc();
    }else{
      mesh[i].CalFissionSrcAdjoint();
    };
  };
};

void GeneralSystem::SetZeroFissionSrc()
{
  for(int i=0;i<TotM;i++){
    mesh[i].PutFissionSrc(0.);
  };
};

void GeneralSystem::SetZeroScatSrc()
{
  for(int g=0;g<grp;g++){
    SetZeroScatSrc(g);
  };
};

void GeneralSystem::SetZeroUpScatSrc()
{
  for(int g=0;g<grp;g++){
    SetZeroUpScatSrc(g);
  };
};

void GeneralSystem::SetZeroScatSrc(int g)
{
  for(int i=0;i<TotM;i++){
    mesh[i].SetZeroScatSrc(g);
  };
};

void GeneralSystem::SetZeroUpScatSrc(int g)
{
  for(int i=0;i<TotM;i++){
    mesh[i].SetZeroUpScatSrc(g);
  };
};

void GeneralSystem::SetZeroScalarFlux()
{
  for(int g=0;g<grp;g++){
    SetZeroScalarFlux(g);
  };
};

void GeneralSystem::SetZeroScalarFlux(int g)
{
  for(int i=0;i<TotM;i++){
    mesh[i].SetZeroScalarFlux(g);
  };
};

GroupData1D GeneralSystem::GetIntegratedReactionRate(enum xstype react,int medid)
{
  GroupData1D ret(grp);  

  vector<real> dat(grp,0.);

  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    if(medid==-1||mi.GetFMatParMesh(i)==medid){
    for(int j=0;j<grp;j++){
      dat[j]+=mesh[i].GetMed()->GetData1D(react).get_dat(j)*
	mesh[i].GetFlux().get_dat(j)*vol;
    };
    };
  };
  

  for(int i=0;i<grp;i++){
    ret.put_data(i,dat[i]);
  };

  /*
  ret.set_zero();
  for(int i=0;i<nmed;i++){
    if(medid==-1||medid==i){
      ret=ret+med[i].GetMacxs().GetData1d(react)*GetIntegratedFlux(i);
    };
  };
  */
  
  return ret;


  
};

GroupData1D GeneralSystem::GetIntegratedReactionRatePerMesh(enum xstype react,int meshid)
{
  vector<real> dat(grp,0.);

  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    if(meshid==-1||i==meshid){
    for(int j=0;j<grp;j++){
      dat[j]+=mesh[i].GetMed()->GetMacxs().GetData1d(react).get_dat(j)*
	mesh[i].GetFlux().get_dat(j)*vol;
    };
    };
  };
  GroupData1D ret(grp);
  for(int i=0;i<grp;i++){
    ret.put_data(i,dat[i]);
  };
  return ret;
};

void GeneralSystem::ShowNeutronFluxAlongXAxis(int y,int z,int g)
{
  if(y>=mi.GetFMesh(1)||z>=mi.GetFMesh(2)){
    cout<<"# Error in GeneralSystem::ShowNeutronFluxAlongXAxis.\n";
    cout<<"# y or z position is inappropriate.\n";
    exit(0);
  };
  if(g>=grp){
    cout<<"# Error in GeneralSystem::ShowNeutronFluxAlongXAxis.\n";
    cout<<"# Energy group is inappropriate.\n";
    exit(0);
  };

  int mesh_num=mi.GetFMesh(0);
  real pos=0.;
  cout<<"# Neutron flux distribution along X-axis\n";
  cout<<"#\n";
  cout<<"#    Energy group [0-"<<grp-1<<"] : "<<g<<"\n";
  cout<<"#    Y-position [0-"<<mi.GetFMesh(1)-1<<"] : "<<y<<"\n";
  cout<<"#    Z-position [0-"<<mi.GetFMesh(2)-1<<"] : "<<z<<"\n";
  cout<<"#\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"# Mesh-center   Flux\n";
  for(int i=0;i<mesh_num;i++){
    pos+=mi.GetFMeshL(0,i)*0.5;
    int ind=meshid[z][y][i];
    cout<<" "<<pos<<"   ";
    if(g==-1){
      for(int j=0;j<grp;j++){
	cout<<mesh[ind].GetFlux().get_dat(j)<<" ";
      };
    }else{
      cout<<mesh[ind].GetFlux().get_dat(g);
    };
    cout<<"\n";
    pos+=mi.GetFMeshL(0,i)*0.5;
  };
};

void GeneralSystem::CalReactionRate
(int x1,int x2,int y1,int y2,int z1,int z2,int matnum,enum xstype react,bool normalize)
{
  int tmp=(x2-x1+1)*(y2-y1+1)*(z2-z1+1);
  vector<real> dummy(tmp);
  CalReactionRate(x1,x2,y1,y2,z1,z2,matnum,react,dummy,normalize);
};

void GeneralSystem::CalReactionRate
(int x1,int x2,int y1,int y2,int z1,int z2,int matnum,enum xstype react,vector<real> &rra,bool normalize)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  
  real r1=1;
  real zv=0.;
  int nn=0;
  for(int z=0;z<mi.GetZF();z++){
    zv+=mi.GetFMeshL(2,z)*0.5;
    real yv=0.;
    for(int y=0;y<mi.GetYF();y++){
      yv+=mi.GetFMeshL(1,y)*0.5;
      real xv=0.;
      for(int x=0;x<mi.GetXF();x++){
	xv+=mi.GetFMeshL(0,x)*0.5;
	real rr;
	if(x>=x1&&x<=x2&&y>=y1&&y<=y2&&z>=z1&&z<=z2){
	  int ind=meshid[z][y][x];
	  bool ExistNuc=mesh[ind].GetMed()->ExistNuclide(matnum);
	  if(ExistNuc){
	    rr=mesh[ind].GetReactionRate(matnum,react);
	    if(x==x1&&y==y1&&z==z1)r1=rr;
	    cout<<"  "<<xv<<"  "<<yv<<"  "<<zv<<"    ";
	    cout<<rr/r1<<"\n";
	    //cout<<rr<<" "<<r1<<"\n";	    
	    //cout<<rr<<"\n";
	    if(normalize)rra[nn]=rr/r1;
	    if(!normalize)rra[nn]=rr;
	    nn++;
	  };
	};
	xv+=mi.GetFMeshL(0,x)*0.5;
      };
      yv+=mi.GetFMeshL(1,y)*0.5;
    };
    zv+=mi.GetFMeshL(2,z)*0.5;
  };
};

void GeneralSystem::CalReactionRateRatio(int x1,int x2,int y1,int y2,int z1,int z2,int matnum1,enum xstype react1,int matnum2,enum xstype react2)
{
  real zv=0.;
  for(int z=0;z<mi.GetZF();z++){
    if(z!=0)zv+=mi.GetFMeshL(2,z)*0.5;
    real yv=0.;
    for(int y=0;y<mi.GetYF();y++){
      if(y!=0)yv+=mi.GetFMeshL(1,y)*0.5;
      real xv=0.;
      for(int x=0;x<mi.GetXF();x++){
	if(x!=0)xv+=mi.GetFMeshL(0,x)*0.5;
	if(x>=x1&&x<=x2&&y>=y1&&y<=y2&&z>=z1&&z<=z2){
	  int ind=meshid[z][y][x];
	  bool ExistNuc1=mesh[ind].GetMed()->ExistNuclide(matnum1);
	  bool ExistNuc2=mesh[ind].GetMed()->ExistNuclide(matnum2);
	  if(ExistNuc1&&ExistNuc2){
            real rr1=mesh[ind].GetReactionRate(matnum1,react1);
            real rr2=mesh[ind].GetReactionRate(matnum2,react2);
	    cout.setf(ios::fixed);
	    cout.precision(3);
	    cout<<"  "<<xv<<"  "<<yv<<"  "<<zv<<"    ";
	    cout.precision(5);
	    cout<<rr1/rr2<<"\n";
	  };
	};
	xv+=mi.GetFMeshL(0,x)*0.5;
	if(x==0)xv+=mi.GetFMeshL(0,x)*0.5;
      };
      yv+=mi.GetFMeshL(1,y)*0.5;
      if(y==0)yv+=mi.GetFMeshL(1,y)*0.5;
    };
    zv+=mi.GetFMeshL(2,z)*0.5;
    if(z==0)zv+=mi.GetFMeshL(2,z)*0.5;
  };
};

real GeneralSystem::CalConversionRatio()
{
  int fert_num=4;
  int fert_mat[]={902320,922340,922380,942400};
  enum xstype fert_xs[]={sigc,sigc,sigc,sigc};

  int fert_num_neg=1;
  int fert_mat_neg[]={912330};
  enum xstype fert_xs_neg[]={sigc};

  int fiss_num=4;
  int fiss_mat[]={922330,922350,942390,942410};
  enum xstype fiss_xs_sigc[]={sigc,sigc,sigc,sigc};
  enum xstype fiss_xs_sigf[]={sigf,sigf,sigf,sigf};

  real fert_cap=CalMacroscopicReactionRate(fert_num,fert_mat,fert_xs);
  real fert_cap_neg=CalMacroscopicReactionRate(fert_num_neg,fert_mat_neg,fert_xs_neg);
  real fiss_cap=CalMacroscopicReactionRate(fiss_num,fiss_mat,fiss_xs_sigc); 
  real fiss_fis=CalMacroscopicReactionRate(fiss_num,fiss_mat,fiss_xs_sigf);

  return (fert_cap-fert_cap_neg)/(fiss_cap+fiss_fis);
};

real GeneralSystem::CalMacroscopicReactionRate(int nuc,int *nuc_id,enum xstype *xs)
{
  bool *on_mesh=new bool[TotM];
  for(int i=0;i<TotM;i++){
    on_mesh[i]=true;
  };
  real ret=CalMacroscopicReactionRate(nuc,nuc_id,xs,on_mesh);
  delete [] on_mesh;
  return ret;
};

real GeneralSystem::CalMacroscopicReactionRate(int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh)
{
  real rr=0.;
  for(int i=0;i<TotM;i++){
    if(on_mesh[i]){
      real vol=mesh[i].GetVolume();
      for(int j=0;j<nuc;j++){
        //cout<<i<<" "<<vol<<" "<<j<<"/"<<nuc<<" "<<nuc_id[j]<<" "<<rr<<"\n";
        if(mesh[i].GetMed()->ExistNuclide(nuc_id[j])){
          real volden=vol*mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetDensity();
	  //cout<<volden<<"\n";
          rr+=mesh[i].GetFlux()*mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetMicxs().GetData1d(xs[j])*volden;
	};
	//cout<<" ...\n";
      };
    };
  };
  return rr;
};

real GeneralSystem::CalPointReactionRateRatio
(int x,int matnum1,enum xstype react1,int matnum2,enum xstype react2)
{
  return CalPointReactionRateRatio(x,0,0,matnum1,react1,matnum2,react2);
};

real GeneralSystem::CalPointReactionRateRatio
(int x,int y,int matnum1,enum xstype react1,int matnum2,enum xstype react2)
{
  return CalPointReactionRateRatio(x,y,0,matnum1,react1,matnum2,react2);
};

real GeneralSystem::CalPointReactionRateRatio
(int x,int y,int z,int matnum1,enum xstype react1,int matnum2,enum xstype react2)
{
  int ind=meshid[z][y][x];
  bool ExistNuc1=mesh[ind].GetMed()->ExistNuclide(matnum1);
  bool ExistNuc2=mesh[ind].GetMed()->ExistNuclide(matnum2);
  if(ExistNuc1&&ExistNuc2){
    real rr1=mesh[ind].GetReactionRate(matnum1,react1);
    real rr2=mesh[ind].GetReactionRate(matnum2,react2);
    return rr1/rr2;
  };
  return 0.;
};

void GeneralSystem::CalAxialPowerDistribution(int x,int y)
{
  cout<<"# Axial power distribution at (x/y)="<<x<<"/"<<y<<")\n";
  vector<real> pow(mi.GetFMesh(2));
  real totp=0.;
  for(int z=0;z<mi.GetFMesh(2);z++){
    GroupData1D sigf=GetMesh(x,y,z).GetMed()->GetMacroSigf();
    real v=GetMesh(x,y,z).GetVolume();
    pow[z]=GetMesh(x,y,z).GetFlux()*sigf*v;
    totp+=pow[z];
  };

  cout.setf(ios::showpoint);
  cout.precision(5);
  real zpos=0.;
  for(int z=0;z<mi.GetFMesh(2);z++){
    real zpos2=zpos;
    zpos+=mi.GetFMeshL(2,z);
    cout<<zpos2<<"-"<<zpos<<" : "<<pow[z]/totp<<"\n";    
  };
};

void GeneralSystem::ShowFissionSourceForXYPlot(int zpl)
{
  real xpos=0.;
  for(int x=0;x<mi.GetXF();x++){
    for(int j=0;j<2;j++){
      if(j==1)xpos+=mi.GetFMeshL(0,x);
      real ypos=0.;
      for(int y=0;y<mi.GetYF();y++){
        real src=GetMesh(x,y,zpl).GetFissionSrc();
        cout<<xpos<<" "<<ypos<<" "<<src<<"\n";
        ypos+=mi.GetFMeshL(1,y);
        cout<<xpos<<" "<<ypos<<" "<<src<<"\n";
      };
      cout<<"\n";
    };
  };
};

void GeneralSystem::ShowIntegratedFluxForXYPlot(int gtop, int gend, int zpl)
{
  real xpos=0.;
  for(int x=0;x<mi.GetXF();x++){
    for(int j=0;j<2;j++){
      if(j==1)xpos+=mi.GetFMeshL(0,x);
      real ypos=0.;
      for(int y=0;y<mi.GetYF();y++){
	real flxsum=0.;
	for(int g=gtop;g<=gend;g++){
          flxsum+=GetMesh(x,y,zpl).GetFlux(0).get_dat(g);
	};
	//flxsum*=GetMesh(x,y,zpl).GetVolume();
        cout<<xpos<<" "<<ypos<<" "<<flxsum<<"\n";
        ypos+=mi.GetFMeshL(1,y);
        cout<<xpos<<" "<<ypos<<" "<<flxsum<<"\n";
      };
      cout<<"\n";
    };
  };
};

void GeneralSystem::ShowReactionRateForXYPlot(GroupData1D &wgt, int zpl)
{
  real xpos=0.;
  for(int x=0;x<mi.GetXF();x++){
    for(int j=0;j<2;j++){
      if(j==1)xpos+=mi.GetFMeshL(0,x);
      real ypos=0.;
      for(int y=0;y<mi.GetYF();y++){
        real rr=GetMesh(x,y,zpl).GetFlux(0)*wgt;
        cout<<xpos<<" "<<ypos<<" "<<rr<<"\n";
        ypos+=mi.GetFMeshL(1,y);
        cout<<xpos<<" "<<ypos<<" "<<rr<<"\n";
      };
      cout<<"\n";
    };
  };
};

void GeneralSystem::CalPowerXY(bool fission_source)
{
  int xc=mi.GetXC();
  int yc=mi.GetYC();
  int xyc=xc*yc;
  real *pow=new real[xyc];
  real *vol=new real[xyc];
  for(int i=0;i<xyc;i++){pow[i]=0.; vol[i]=0.;};

  int index=0;
  int index2=0;
  real totpow=0.;
  real totvol=0.;
  for(int z=0;z<mi.GetZC();z++){
    for(int z2=0;z2<mi.GetCMeshF(2,z);z2++){
      for(int y=0;y<mi.GetYC();y++){
        for(int y2=0;y2<mi.GetCMeshF(1,y);y2++){
          for(int x=0;x<mi.GetXC();x++){
            for(int x2=0;x2<mi.GetCMeshF(0,x);x2++){
	      if(mi.GetFMat(index)!=-1){
		GroupData1D sigf=mesh[index2].GetMed()->GetMacxs().GetData1d(nusigf);
		if(!fission_source){
                  sigf=mesh[index2].GetMed()->GetMacroSigf();
		};
                real v=mesh[index2].GetVolume();
	        real p=mesh[index2].GetFlux()*sigf*v;
	        pow[y*mi.GetXC()+x]+=p;
	        totpow+=p;
	        vol[y*mi.GetXC()+x]+=v;
	        if(p>0.)totvol+=v;
		index2++;
	      };
              index++;
	    };
	  };
	};
      };
    };
  };

  /*
  totpow/=totvol;
  for(int i=0;i<xyc;i++){
    pow[i]/=vol[i];
  };
  */

  cout<<"# Radial power distribution (fission reaction rate)\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<yc;i++){
    for(int j=0;j<xc;j++){
      if(pow[i*xc+j]>0.){
        cout<<"(x/y)=("<<j<<"/"<<i<<") "<<pow[i*xc+j]/totpow<<"\n";
      };
    };
  };

  delete [] pow;
  delete [] vol;
}

void GeneralSystem::ShowIntegratedFlux(GroupData1D &flx)
{
  cout<<"# E_up     Flux       Flux/lethargy [normalized]\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  GroupData1D ebnd=med[0].GetEnband();
  real sum=flx.get_sum();
  for(int g=0;g<grp;g++){
    real e0=ebnd.get_dat(g);
    real e1=ebnd.get_dat(g+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<flx.get_dat(g)<<" "<<flx.get_dat(g)/letwid<<"    "<<flx.get_dat(g)/letwid/sum<<"\n";
    if(g==grp-1)cout<<e1<<" "<<flx.get_dat(g)<<" "<<flx.get_dat(g)/letwid<<"    "<<flx.get_dat(g)/letwid/sum<<"\n";
  };
};

void GeneralSystem::ShowIntegratedFluxMeshID(int m1,int m2)
{
  cout<<"#\n# Integrated neutron flux (mesh: "<<m1<<" to "<<m2<<")\n#\n";
  GroupData1D flx=GetIntegratedFluxMeshID(m1,m2);
  ShowIntegratedFlux(flx);
};

void GeneralSystem::ShowIntegratedFlux(int medid,int mom)
{
  cout<<"#\n# Integrated neutron flux (medium ID :  "<<medid<<")\n#\n";
  GroupData1D flx=GetIntegratedFlux(medid,mom);
  ShowIntegratedFlux(flx);
};

GroupData1D GeneralSystem::GetIntegratedFlux(int medid,int mom)
{
  if(name=="PJI"||name=="MEC"){
    cout<<"GetIntegratedFlux is not coded for "<<name<<" with IrregularGeometryInfo.\n";
    exit(0);
  };

  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  for(int i=0;i<TotM;i++){
    //cout<<i<<" "<<medid<<" "<<mi.GetFMatParMesh(i)<<"\n";
    if(mi.GetFMatParMesh(i)==medid){
      real vol=mesh[i].GetVolume();
      for(int g=0;g<grp;g++){
	sum[g]+=fabs(mesh[i].GetFluxData(mom,g))*vol;
	//cout<<vol<<" "<<g<<" "<<mesh[i].GetFluxData(mom,g)<<"\n";
      };
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};

GroupData1D GeneralSystem::GetIntegratedFluxMeshID(int m1,int m2)
{
  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  for(int i=m1;i<=m2;i++){
    real vol=mesh[i].GetVolume();
    for(int g=0;g<grp;g++){
      sum[g]+=fabs(mesh[i].GetFluxData(0,g))*vol;
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};

real GeneralSystem::GetVolumePerMedium(int medid)
{
  //if(name=="PJI"||name=="MEC"){
  if(name=="MEC"){
    cout<<"Error in GeneralSystem::GetVolumePerMedium,\n";
    cout<<"This method cannot be used in "<<name<<" solver.\n";
    exit(0);
  };

  real volsum=0.;
  for(int i=0;i<TotM;i++){
    if(mi.GetFMatParMesh(i)==medid){
      volsum+=mesh[i].GetVolume();
    };
  };
  return volsum;
};

GroupData1D GeneralSystem::GetIntegratedFluxPerVolume(int medid,int mom)
{
  if(name=="PJI"||name=="MEC"){
    cout<<"GetIntegratedFlux is not coded for "<<name<<" with IrregularGeometryInfo.\n";
    exit(0);
  };

  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  real volsum=0.;
  for(int i=0;i<TotM;i++){
    if(mi.GetFMatParMesh(i)==medid){
      real vol=mesh[i].GetVolume();
      volsum+=vol;
      for(int g=0;g<grp;g++){
	sum[g]+=fabs(mesh[i].GetFluxData(mom,g))*vol;
      };
    };
  };
  volsum=1./volsum;
  for(int i=0;i<grp;i++){
    sum[i]*=volsum;
  };

  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};

GroupData1D GeneralSystem::GetIntegratedFlux(int x1,int x2,int y1,int y2,int z1,int z2,int mom)
{
  if(name=="PJI"||name=="MEC"){
    cout<<"GetIntegratedFlux is not coded for "<<name<<" with IrregularGeometryInfo.\n";
    exit(0);
  };

  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};

  int ind=0;
  int z=0;
  for(int i=0;i<mi.GetZC();i++){
    for(int ii=0;ii<mi.GetCMeshF(2,i);ii++){
      int y=0;
      for(int j=0;j<mi.GetYC();j++){
	for(int jj=0;jj<mi.GetCMeshF(1,j);jj++){
	  int x=0;
	  for(int k=0;k<mi.GetXC();k++){
	    for(int kk=0;kk<mi.GetCMeshF(0,k);kk++){
	      if(x>=x1&&x<=x2&&
                 y>=y1&&y<=y2&&
                 z>=z1&&z<=z2){
                real vol=mesh[i].GetVolume();
                for(int g=0;g<grp;g++){
		  sum[g]+=fabs(mesh[ind].GetFluxData(mom,g))*vol;
                };
	      };
	      ind++;
	      x++;
	    };
	  };
	  y++;
	};
      };
      z++;
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};

GroupData1D GeneralSystem::GetIntegratedFluxPerVolume(int x1,int x2,int y1,int y2,int z1,int z2,int mom)
{
  if(name=="PJI"||name=="MEC"){
    cout<<"GetIntegratedFlux is not coded for "<<name<<" with IrregularGeometryInfo.\n";
    exit(0);
  };

  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};

  int ind=0;
  int z=0;
  real volsum=0.;
  for(int i=0;i<mi.GetZC();i++){
    for(int ii=0;ii<mi.GetCMeshF(2,i);ii++){
      int y=0;
      for(int j=0;j<mi.GetYC();j++){
	for(int jj=0;jj<mi.GetCMeshF(1,j);jj++){
	  int x=0;
	  for(int k=0;k<mi.GetXC();k++){
	    for(int kk=0;kk<mi.GetCMeshF(0,k);kk++){
	      if(x>=x1&&x<=x2&&
                 y>=y1&&y<=y2&&
                 z>=z1&&z<=z2){
                real vol=mesh[i].GetVolume();
                volsum+=vol;
                for(int g=0;g<grp;g++){
		  sum[g]+=fabs(mesh[ind].GetFluxData(mom,g))*vol;
                };
	      };
	      ind++;
	      x++;
	    };
	  };
	  y++;
	};
      };
      z++;
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);

  volsum=1./volsum;
  ret=ret*volsum;

  delete []sum;
  return ret;
};

// **********************************************
// * For fixed source calculation               *
// **********************************************

void GeneralSystem::PutIsotropicSourcePerUnitVolume
(int x1,int x2,int y1,int y2,int z1,int z2,GroupData1D &src,int mom)
{
  for(int z=z1;z<=z2;z++){
    for(int y=y1;y<=y2;y++){
      for(int x=x1;x<=x2;x++){
	int tmp=meshid[z][y][x];
	PutIsotropicSourcePerUnitVolume(tmp,src,mom);
	//real vol=mesh[tmp].GetVolume()*INV_PI4;
	//for(int g=0;g<grp;g++){
	//  mesh[tmp].PutScatSrc(g,vol*src.get_dat(g),mom);
	//};
      };
    };
  };
};

void GeneralSystem::PutIsotropicSourcePerUnitVolume
(int mesh_id,GroupData1D &src,int mom)
{
  //src.show_self();
  real vol=mesh[mesh_id].GetVolume()*INV_PI4;
  for(int g=0;g<grp;g++){
    mesh[mesh_id].PutScatSrc(g,vol*src.get_dat(g),mom);
  };
};

void GeneralSystem::CalFixedSource(real epsf,int oiter,bool print_inp)
{
  SetZeroFissionSrc();
  real err;
  for(int g=0;g<grp;g++){
    int ginp=g;
    if(!opt.Forward())ginp=grp-1-g;
    if(print&&print_inp)cout<<"# Fixed source calculation in group "<<ginp<<"\n";
    CalSrcMultiplySystem(ginp,1.,pl);
    err=CalFluxGeneral(ginp,epsf,oiter);
    AddDownScatSrc(ginp,pl);
  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };
};

void GeneralSystem::CalFixedSourceUpScat(real epsf, int oiter, bool print)
{
  // - External source is stored in the scattering source array.
  // - Fission source is set to zero.

  CheckUpScattering();
  SetZeroFissionSrc();
  
  if(!IsUpScatterSystem){
    cout<<"#\n";
    cout<<"# Warning in GeneralSystem::CalFixedSourceUpScat.\n";
    cout<<"# Since this system is NOT up-scattering-included system,\n";
    cout<<"# [CalFixedSource] is runned instead.\n";
    cout<<"#\n";
    CalFixedSource(epsf,oiter,print);
    return;
  };

  //int itermax=50;
  int itermax=5000;
  real err;

  // The new forward part below is verified in 2017/11/2.
  // Anisotropy of stored source data is important.
  if(opt.Forward()){

    // +++ FORWARD CALCULATION +++
    /*
    vector<real> flxold(TotM);    
    real sor_factor=1.5;    
    */
    
    // (initial guess of up-scattering source)
    if(oiter!=0){
      for(int g=0;g<grp;g++){
        if(g>=UpScatSrcGrp)AddUpScatSrc(g,pl);
        SetZeroUpScatSrc(g);
      };
    };

    // +++ Fixed source term calculation for thermal iteration
    for(int g=0;g<UpScatSinkGrp;g++){
      if(print)cout<<"# Fixed source calculation in group "<<g<<"\n";      
      CalSrcMultiplySystem(g,1.,pl);
      err=CalFluxGeneral(g,epsf,0);
      //cout<<"            "<<g<<" "<<err<<"\n";
      AddDownScatSrc(g,pl);
    };

    vector< vector< vector<real> > > src_fixed(TotM);
    for(int i=0;i<TotM;i++){
      src_fixed[i].resize(plnum);
      for(int l=0;l<plnum;l++){
        src_fixed[i][l].resize(grp);
        for(int g=UpScatSinkGrp;g<grp;g++){
          src_fixed[i][l][g]=mesh[i].GetScatSrc(g,l);
	};
      };
    };

    // +++ Thermal iteration
    vector<bool> convergence(grp,false);
    for(int i=0;i<itermax;i++){
      bool conv=true;
      if(i==0)conv=false;
      for(int g=UpScatSinkGrp;g<grp;g++){
	if(i!=0){
          for(int j=0;j<TotM;j++){
	    for(int l=0;l<plnum;l++){
      	      mesh[j].AddScatSrc(g,src_fixed[j][l][g],l);
	    };
	  };
	};	    
	if(!convergence[g]){

          if(print)cout<<"# Fixed source calculation in group "<<g<<"\n";      	  
	  /*
  	  for(int j=0;j<TotM;j++){
	    flxold[j]=GetMesh(j).GetFlux().get_dat(g);
	  };
	  */
	  
  	  CalSrcMultiplySystem(g,1.,pl);
	  err=CalFluxGeneral(g,epsf,i);
	  err=fabs(err);
	  //if(g==100)cout<<i<<"         "<<g<<" : "<<err<<"\n";
	  //cout<<i<<"         "<<g<<" : "<<err<<"\n";
	  //if(err<epsf*0.1){
	  /*
  	  for(int j=0;j<TotM;j++){
	    real newflx=mesh[j].GetFlux().get_dat(g);
	    mesh[j].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[j]);
	  };
	  */
	  
	  if(err<epsf){	    
	    convergence[g]=true;
	  };
	};
        AddDownScatSrc(g,pl);
        AddUpScatSrc(g,pl);
        SetZeroScatSrc(g);
        SetZeroUpScatSrc(g);
        if(err>epsf)conv=false;
      };
      if(conv){
	cout<<"#    Number of iterations : "<<i<<"\n";
        break;
      };
      if(i==itermax-1){
        cout<<"# !! Caution !!\n";
        cout<<"# Fixed source calculation is not converged.\n";
        cout<<"# Maximum residual is "<<err<<"\n";
        cout<<"# Convergence criteria "<<epsf<<"\n";
      };
    };

  };

  // +++ ADJOINT calculation
  if(!opt.Forward()){

    vector<real> flxold(TotM);

    // (external source is stored)
    vector< vector<real> > src_inp(TotM);
    for(int i=0;i<TotM;i++){
      src_inp[i].resize(grp);
      for(int j=0;j<grp;j++){
        src_inp[i][j]=mesh[i].GetScatSrc(j);
	//if(src_inp[i][j]>0)cout<<i<<" "<<j<<" "<<src_inp[i][j]<<"\n";
      };
    };

    // (thermal iteration)
    // The SOR factor of 1.2 is introduced to accelerate adjoint calculation of MulticellBurner,
    // but it fails in Burner calculation, so the SOR factor is reset to 1.0 in 2019/1/31.

    //real sor_factor=1.2;
    real sor_factor=1.0;

    real cin=0.1;  
    real errf=1.0;
    real epsf=1e-4;
    for(int i=0;i<itermax;i++){

      if(errf*0.02<cin)cin=errf*0.02;
      if(cin<epsf)cin=epsf;

      bool conv=true;
      if(i==0)conv=false;
      errf=0.;
      for(int g=grp-1;g>=UpScatSinkGrp;g--){
        if(i!=0){
          for(int j=0;j<TotM;j++){
  	    mesh[j].AddScatSrc(g,src_inp[j][g]);
  	  };
        };
	CalSrcMultiplySystem(g,1.,pl);

	for(int j=0;j<TotM;j++){
	  flxold[j]=GetMesh(j).GetFlux().get_dat(g);
	};

	//real eps_inner=1e-4;
        //err=CalFluxGeneral(g,eps_inner,i);
        err=CalFluxGeneral(g,cin,i);
        err=fabs(err);
	if(err>errf)errf=err;

	for(int j=0;j<TotM;j++){
	  real newflx=mesh[j].GetFlux().get_dat(g);
	  //if(i==0)cout<<i<<" "<<g<<" "<<j<<" "<<flxold[j]<<" "<<newflx<<"\n";
	  mesh[j].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[j]);
	};

        AddDownScatSrc(g,pl);
        AddUpScatSrc(g,pl);
        SetZeroScatSrc(g);
        SetZeroUpScatSrc(g);
	//cout<<g<<" "<<err<<"\n";
        if(err>epsf)conv=false;
      };
      if(conv){
        //cout<<"#  First-Iteration : "<<i<<"\n";
        break;
      }else{
	for(int g=UpScatSinkGrp-1;g>=0;g--){
          SetZeroScatSrc(g);
	};
      };

      if(i==itermax-1){
        cout<<"# !! Caution !!\n";
        cout<<"# Fixed source calculation is not converged.\n";
        cout<<"# Maximum residual is "<<err<<"\n";
        cout<<"# Convergence criteria "<<epsf<<"\n";
      };

    };

    for(int g=UpScatSinkGrp-1;g>=0;g--){
      for(int j=0;j<TotM;j++){
        mesh[j].AddScatSrc(g,src_inp[j][g]);
      };
      CalSrcMultiplySystem(g,1.,pl);
      err=CalFluxGeneral(g,epsf,0);
      AddDownScatSrc(g,pl);
    };

    if(name=="SNR")FluxAngleReverse();

    return;
  };

  // Original
  /*
  vector< vector<real> > src_inp(TotM);
  for(int i=0;i<TotM;i++){
    src_inp[i].resize(grp);
    for(int j=0;j<grp;j++){
      src_inp[i][j]=mesh[i].GetScatSrc(j);
      //cout<<" src_inp : "<<i<<" "<<j<<" "<<src_inp[i][j]<<"\n";
    };
  };

  for(int i=0;i<itermax;i++){

    if(i!=0){
      for(int j=0;j<TotM;j++){
	for(int k=0;k<grp;k++){
	  //mesh[j].PutScatSrc(k,src_inp[j][k]);
	  mesh[j].AddScatSrc(k,src_inp[j][k]);
	};
      };
    };

    bool conv=true;
    if(i==0)conv=false;
    //real errmax=0.;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      if(opt.Forward()||
	 //if((opt.Forward()&&ginp>=UpScatSinkGrp)||
        (!opt.Forward()&&ginp>=UpScatSinkGrp)||
        (!opt.Forward()&&ginp<UpScatSinkGrp&&conv)){
        //if(print)cout<<"# Fixed source calculation in group "<<ginp<<"\n"; 
        if(i!=0){
          for(int i=0;i<TotM;i++){
    	    //mesh[i].AddScatSrc(ginp,src_inp[i][ginp]);
	  };
	};
	CalSrcMultiplySystem(ginp,1.,pl);
        err=CalFluxGeneral(ginp,1e-4,i);
        err=fabs(err);
	//if(err>errmax)errmax=err;
	//if(g==100)cout<<i<<"         "<<ginp<<" : "<<err<<"\n";
        AddDownScatSrc(ginp,pl);
        AddUpScatSrc(ginp,pl);
      };
      SetZeroScatSrc(ginp);
      SetZeroUpScatSrc(ginp);
      if(opt.Forward()){
        if(err>epsf)conv=false;
        //cout<<"    "<<i<<" "<<g<<" "<<err<<"\n";
      }else{
        //cout<<"    "<<g<<" "<<err<<"\n";
        if(ginp>=UpScatSinkGrp&&err>epsf)conv=false;
      };
    };
    //cout<<"     "<<errmax<<"\n";
    if(conv){
      //cout<<"#  First-Iteration : "<<i<<"\n";
      //cout<<"#  Iteration : "<<i<<"\n";
      break;
    };


    if(i==itermax-1){
      cout<<"# !! Caution !!\n";
      cout<<"# Fixed source calculation is not converged.\n";
      cout<<"# Maximum residual is "<<err<<"\n";
      cout<<"# Convergence criteria "<<epsf<<"\n";
    };

  };
  */


  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };

};

void GeneralSystem::CalFixedFissionSourceUpScat(real k)
{
  //CheckUpScattering();

  if(!IsUpScatterSystem){
    cout<<"# Error in GeneralSystem::CalFixedFissionSourceUpScat.\n";
    cout<<"# The method [CheckUpScattering] has not yet been done.\n";
    exit(0);
  };
  for(int g=0;g<grp;g++){
    int ginp=g;
    if(!opt.Forward())ginp=grp-1-g;
    AddUpScatSrc(ginp,pl);
    SetZeroUpScatSrc(ginp);
  };

  int itermax=500;
  real err;
  real epsf=1e-4;  
  real sor_factor=1.2;

  real cin=0.1;
  real errf=1.0;

  vector<real> flxold(TotM);

  for(int i=0;i<itermax;i++){

    if(errf*0.02<cin)cin=errf*0.02;
    if(cin<epsf)cin=epsf;

    bool conv=true;
    if(i==0)conv=false;
    errf=0.;

    /*
    // ...... chibatmp
    real sum_source_pre=0.;
    for(int m=0;m<TotM;m++){
      sum_source_pre+=mesh[m].GetFissionSrc();
    };
    */

    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      if(opt.Forward()||(!opt.Forward()&&ginp>=UpScatSinkGrp)||(!opt.Forward()&&ginp<UpScatSinkGrp&&conv)){
        CalSrcMultiplySystem(ginp,k,pl);

	for(int j=0;j<TotM;j++){
	  flxold[j]=GetMesh(j).GetFlux().get_dat(ginp);
	};

        err=CalFluxGeneral(ginp,cin,i); // Default?
        //err=CalFluxGeneral(ginp,1e-4,i);
        //err=CalFluxGeneral(ginp,1e-4,0);

	if(i!=0&&!opt.Forward()&&ginp>=UpScatSinkGrp){
	for(int j=0;j<TotM;j++){
	  real newflx=mesh[j].GetFlux().get_dat(ginp);
	  mesh[j].GetFlux().put_data(ginp,sor_factor*newflx-(sor_factor-1.)*flxold[j]);
	};
	};

        err=fabs(err);
	if(err>errf)errf=err;
        AddDownScatSrc(ginp,pl);
        AddUpScatSrc(ginp,pl);
      };
      SetZeroScatSrc(ginp);
      SetZeroUpScatSrc(ginp);
      if(opt.Forward()){
        if(err>epsf)conv=false;
      }else{
        //cout<<"    "<<ginp<<" ... "<<err<<"\n";
	if(ginp>=UpScatSinkGrp&&err>epsf)conv=false;
      };
    };
 
    //if(i<20)conv=false;// ! chibatmp
       
    if(conv){
      //cout<<"#  Iteration : "<<i<<"\n";
      break;
    };

    /*
    // ...... chibatmp
    real sum_source_post=0.;
    for(int m=0;m<TotM;m++){
      mesh[m].CalFissionSrcAdjoint();
      sum_source_post+=mesh[m].GetFissionSrc();
    };
    cout.setf(ios::showpoint);
    cout.precision(7);
    cout<<i<<" "<<sum_source_post/sum_source_pre<<" "<<k<<"\n";
    */

  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };

};

void GeneralSystem::CalFixedFissionSource(real k, int oit)
{
  // COMMENT in 2016/06/15
  //
  // I am confused why the outer iteration is implemented in the previous version.
  // Because the fission source is fixed and there is no up-scattering,
  // no outer iteration might be required.
  //
  // Note that [SetZeroScatSrc] is added before the [CalFixedFIssionSource(UpScat)]
  // in the [CalGPT] method of the GeneralSystem class.

  CheckUpScattering();
  if(IsUpScatterSystem){
    cout<<"# Error in GeneralSystem::CalFixedFissionSource.\n";
    cout<<"# This is up-scattered system, so please use [CalFixedFissionSourceUpScat].\n";
    exit(0);
  };

  // (Original part in the update in 2016/06/15)
  /*
  int itermax=50;
  real err;
  real epsf=1e-4;
  for(int i=0;i<itermax;i++){
    bool conv=true;
    if(i==0)conv=false;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,k,pl);
      err=CalFluxGeneral(ginp,1e-4,0);
      //if(i==1)cout<<ginp<<" "<<err<<"\n";
      //cout<<i<<" "<<ginp<<" "<<err<<"\n";
      err=fabs(err);
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
      if(err>epsf)conv=false;
      //for(int m=0;m<TotM;m++){
      //cout<<"      "<<ginp<<" "<<m<<" "<<mesh[m].GetFlux().get_dat(ginp)<<"\n";
      //};
    };
    if(conv){
      //cout<<"# CalFixedFissionSource : "<<i<<"\n";
      break;
    };
  };
  */

  // (New part in the update in 2016/06/15)
  for(int g=0;g<grp;g++){
    int ginp=g;
    if(!opt.Forward())ginp=grp-1-g;
    CalSrcMultiplySystem(ginp,k,pl);
    real err=CalFluxGeneral(ginp,1e-4,oit);
    AddDownScatSrc(ginp,pl);
    SetZeroScatSrc(ginp);
  };

  if(name=="SNR"&&!opt.Forward()){
    FluxAngleReverse();
  };

};


real GeneralSystem::CalFixedSourceWithFission(real epsf, int itermax, bool high_speed_option)
{
  //
  // Subcriticality multiplication rate (k_sub) is returned.
  //

  CheckUpScattering();
  /*
  if(IsUpScatterSystem){
    cout<<"# Error in GeneralSystem::CalFixedSourceWithFission.\n";
    cout<<"# Not coded for up-scattering problems.\n";
    exit(0);
  };
  */

  if(print){
  cout<<"#\n";
  cout<<"# Fixed source calculations with fission\n";
  cout<<"#\n";
  };

  vector< vector<real> > flx_store(TotM);
  vector< vector<real> > flxold(TotM); // to estimate maximum eigenvalue from residual
  for(int m=0;m<TotM;m++){
    flx_store[m].resize(grp,0.);
    flxold[m].resize(grp);
  };

  // +++ Condition to activate high-speed calculation +++
  // 
  //  Original it has been 1e-7, but this was found insufficient 
  //  in kinetic calculations with null-transient conditin
  //  with large time step (extremely-small source under extremely-shallow 
  //  subcritical condition)

  real eps_hso=1e-15; // for high-speed option
  //real eps_hso=1e-7; // for high-speed option

  // initial source
  real ssum=0.;
  for(int m=0;m<TotM;m++){
    for(int g=0;g<grp;g++){
      ssum+=mesh[m].GetScatSrc(g)*PI4;
    };
  };

  if(print)cout<<"# Iter / Source / FlxErr / K\n";
  real k_q=0.;
  real ssum_old=ssum;
  real esrc=ssum;
  real normalization_factor=1.;
  for(int i=0;i<itermax;i++){

    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
	flxold[m][g]=mesh[m].GetFlux().get_dat(g)/normalization_factor;
      };
    };

    //CalFixedSource(1e-4,i,false);
    if(IsUpScatterSystem){
      CalFixedSourceUpScat(epsf,i,false);
      //CalFixedSourceUpScat(1e-4,i,false);
    }else{
      CalFixedSource(epsf,i,false);
      //CalFixedSource(1e-4,i,false);
    };
    
    // estimation of maximum eigenvalue from residual
    real sum=0.;
    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
	flxold[m][g]=mesh[m].GetFlux().get_dat(g)/flxold[m][g]; // lambda
	sum+=flxold[m][g];
      };
    };
    sum/=(TotM*grp); // estimated maximum eigenvalue

    real sum2=0;
    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
	sum2+=pow(flxold[m][g]-sum,2);
      };
    };
    sum2/=(TotM*grp); // variance

    bool end_loop=false;
    if(high_speed_option&&sum2<eps_hso&&i!=0){
      end_loop=true;
      if(print){
      cout<<"# ... higher modes are sufficiently decreased ...\n";
      cout<<"#     (estimated lambda is "<<sum<<" )\n";
      };
    };


    real errmax=0.;
    real factor=1.;
    if(end_loop)factor=1./(1.-sum);
    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
	real org=flx_store[m][g];
        real add=mesh[m].GetFlux().get_dat(g)*factor;
	real err=fabs(add/org);
	if(err>errmax)errmax=err;
	flx_store[m][g]+=add;
      };
    };
    if(print){
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# "<<i<<" "<<ssum<<" "<<errmax<<" ";
    cout.setf(ios::showpoint);
    cout.precision(5);
    cout<<"# "<<ssum/ssum_old<<"\n";
    };
    normalization_factor=ssum/ssum_old;
    //normalization_factor=1.;

    if(errmax<epsf||end_loop)break;

    SetZeroScatSrc();
    CalFissionSrcAllMesh();

    ssum_old=ssum;
    ssum=0.;
    for(int m=0;m<TotM;m++){
      GroupData1D inp(grp);
      real src=mesh[m].GetFissionSrc();
      real vol_inv=1./mesh[m].GetVolume();
      ssum+=src;
      for(int g=0;g<grp;g++){
        real chiv=0.;
	if(opt.Forward()){
          chiv=mesh[m].GetMed()->GetMacxs().GetData1d(chi).get_dat(g);
	}else{
          chiv=mesh[m].GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g);
	};
	inp.put_data(g,chiv*src*vol_inv);
      };
      PutIsotropicSourceParVolume(m,inp);
    };
    if(i==0)k_q=ssum;

    // (flux normalization for accurate estimation of up-scattering source)
    //if(!high_speed_option){
    for(int m=0;m<TotM;m++){
      for(int l=0;l<plnum;l++){
        for(int g=0;g<grp;g++){
  	  real org=mesh[m].GetFlux(l).get_dat(g);
  	  mesh[m].GetFlux(l).put_data(g,org*normalization_factor);
	};
      };
    };
    //};

  };

  // +++ iteration end +++
  for(int m=0;m<TotM;m++){
    for(int g=0;g<grp;g++){
      mesh[m].GetFlux().put_data(g,flx_store[m][g]);
    };
  };
  CalFissionSrcAllMesh();
  real sum_fsrc=0.;
  for(int i=0;i<TotM;i++){  
    sum_fsrc+=mesh[i].GetFissionSrc();
  };

  real ksub=sum_fsrc/(sum_fsrc+esrc);

  if(print&&opt.Forward()){
  cout<<"#\n#+++ External source      : "<<esrc<<"\n";
  cout<<"#+++ Total fission source : "<<sum_fsrc<<"\n";
  k_q/=esrc;
  real k_f=((esrc+sum_fsrc)*ksub-k_q*esrc)/sum_fsrc;
  cout<<"#\n#+++ Subcritical multiplication rate     : "<<ksub<<"\n";
  cout<<"#+++ Source neutron multiplication rate  : "<<k_q<<"\n";
  cout<<"#+++ Fission neutron multiplication rate : "<<k_f<<"\n";
  };

  return ksub;

};

real GeneralSystem::CalFixedSourceWithFissionUpScat()
{
  CheckUpScattering();
  if(!IsUpScatterSystem){
    cout<<"# Error in GeneralSystem::CalFixedSourceWithFissionUpScat.\n";
    cout<<"# This method should NOT used because no up scattering.\n";
    exit(0);
  };

  if(!opt.Forward()){
    cout<<"# Error in GeneralSystem::CalFixedSourceWithFissionUpScat.\n";
    cout<<"# Not yet coded for adjoint problems.\n";
    exit(0);
  };

  // (external source is stored)
  vector< vector<real> > src_inp(TotM);
  for(int i=0;i<TotM;i++){
    src_inp[i].resize(grp);
    for(int j=0;j<grp;j++){
      src_inp[i][j]=mesh[i].GetScatSrc(j);
    };
  };

  /*
  vector< vector<real> > flx_store(TotM);
  vector< vector<real> > flxold(TotM); // to estimate maximum eigenvalue from residual
  for(int m=0;m<TotM;m++){
    flx_store[m].resize(grp,0.);
    flxold[m].resize(grp);
  };

  real eps_hso=1e-7; // for high-speed option

  // initial source
  real ssum=0.;
  for(int m=0;m<TotM;m++){
    for(int g=0;g<grp;g++){
      ssum+=mesh[m].GetScatSrc(g)*PI4;
    };
  };

  if(print)cout<<"# Iter / Source / FlxErr / K\n";
  real k_q=0.;
  real ssum_old=ssum;
  real esrc=ssum;
  */

  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  //SetInitialFlux();

  real cin=0.1; 

  vector< vector< vector<real> > > down_scat_src_store;
  int tmp=grp-UpScatSinkGrp;
  down_scat_src_store.resize(tmp);
  for(int i=0;i<tmp;i++){
    down_scat_src_store[i].resize(TotM);
    for(int j=0;j<TotM;j++){
      down_scat_src_store[i][j].resize(plnum,0.);
    };
  };

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    // (Fixed source is added)
    for(int j=0;j<TotM;j++){
      for(int k=0;k<grp;k++){
        mesh[j].PutScatSrc(k,src_inp[j][k]);
      };
    };

    if(errf*0.02<cin)cin=errf*0.02;
    
    errf=0.;
    bool InnerConvergence=true;

    for(int g=0;g<grp;g++){

      if(g==UpScatSinkGrp){
        for(int i=UpScatSinkGrp;i<grp;i++){
	  for(int j=0;j<TotM;j++){
	    for(int k=0;k<plnum;k++){
	      down_scat_src_store[i-UpScatSinkGrp][j][k]=mesh[j].GetScatSrc(i,k);
	    };
	  };
	};
      };

      CalSrcMultiplySystem(g,1.,pl);
      real err=CalFluxGeneral(g,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      AddDownScatSrc(g,pl);
      SetZeroScatSrc(g);
      if(g>=UpScatSinkGrp)SetZeroUpScatSrc(g);
      if(g>=UpScatSrcGrp)AddUpScatSrc(g,pl);

    };

    // Thermal iteration
    for(int i=0;i<ThermalIteration*3;i++){
      real errf=0.;

      int grptmp=UpScatSinkGrp;
      if(i%3!=0){
        if(grp==107)grptmp=grp-1-20;
      };
      for(int ii=grptmp;ii<grp;ii++){
        for(int j=0;j<TotM;j++){
	  for(int k=0;k<plnum;k++){
	    mesh[j].PutScatSrc(ii,down_scat_src_store[ii-UpScatSinkGrp][j][k],k);
	  };
	};
      };

      for(int g=UpScatSinkGrp;g<grp;g++){
        int ginp=g;
	bool flxcal=false;
	if(i%3==0){
	  flxcal=true;
	}else if(i%3==1&&g>grptmp){
	  flxcal=true;
	}else if(i%3==2&&g>grptmp){
	  flxcal=true;
	};

	if(flxcal){
          real sor_factor=1.2;
  	  CalSrcMultiplySystem(ginp,1.,pl);
          vector<real> flxold(TotM);
	  for(int i=0;i<TotM;i++){
	    flxold[i]=GetMesh(i).GetFlux().get_dat(g);
	  };
	  real err=CalFluxGeneral(ginp,cin,1);
	  //             (important)-------^
	  for(int i=0;i<TotM;i++){
	    real newflx=mesh[i].GetFlux().get_dat(g);
	    mesh[i].GetFlux().put_data(g,sor_factor*newflx-(sor_factor-1.)*flxold[i]);
	  };
	  if(err>errf)errf=err;
	};
        AddDownScatSrc(ginp,pl);
        SetZeroScatSrc(ginp);
        SetZeroUpScatSrc(ginp);
        AddUpScatSrc(ginp,pl);
      };
      //if(print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
      if(i%3==0&&print)cout<<"#  (Thermal iteration) "<<i<<" : "<<errf<<"\n";
    };


    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      mesh[i].CalFissionSrc();
      if(tmp>0.){
	real tmp2=mesh[i].GetFissionSrc();
        real err=fabs(tmp2/tmp-1.0);
        if(err>errs)errs=err;
      };
    };

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss/fissold,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){
        out_iter_end=iter;
        break;
      }else{cout<<"#  ... Inner iteration is not converged ...\n";
      };
    };

  };
};

void GeneralSystem::CalGPT(real keff,int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh,real rr,int itermax)
{
  GroupData1D src(grp);

  SetZeroScatSrc();
  for(int i=0;i<TotM;i++){
    if(on_mesh[i]){
      src.set_zero();
      for(int j=0;j<nuc;j++){
	if(mesh[i].GetMed()->ExistNuclide(nuc_id[j])){
          src=src+mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetMicxs().GetData1d(xs[j])*mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetDensity()/rr;
        };
      };
      PutIsotropicSourcePerUnitVolume(i,src);
    };
  };
  CalGPT(keff,1e-4,itermax);
};

void GeneralSystem::CalGPTForPointMicroReactionRate(real keff,int mesh_id,int nuc_id,enum xstype xs,real rr,int itermax)
{
  if(mesh_id<0||mesh_id>=TotM){
    cout<<"# Error in GeneralSystem::CalGPTForPointMicroReacrionRate.\n";
    cout<<"# Mesh ID "<<mesh_id<<" is incorrect.\n";
    exit(0);
  };

  if(!mesh[mesh_id].GetMed()->ExistNuclide(nuc_id)){
    cout<<"# Error in GeneralSystem::CalGPTForPointMicroReacrionRate.\n";
    cout<<"# Nuclide (ID "<<nuc_id<<") is NOT included in the medium.\n";
    exit(0);
  };

  GroupData1D src(grp);
  src=mesh[mesh_id].GetMed()->GetNuclide(nuc_id).GetMicxs().GetData1d(xs)/rr;
  SetZeroScatSrc();
  PutIsotropicSourcePerUnitVolume(mesh_id,src);

  CalGPT(keff,1e-4,itermax);
};

void GeneralSystem::CalGPTForFissionSource(real keff,real rr,bool subcritical)
{
  GroupData1D src(grp);

  SetZeroScatSrc();
  for(int i=0;i<TotM;i++){
    if(mesh[i].GetMed()->IsFissionable()){
      src.set_zero();
      for(int j=0;j<mesh[i].GetMed()->GetNucnum();j++){
	if(mesh[i].GetMed()->GetNuclideInTurn(j).IsFissionable()){
	  GroupData1D micnusigf=mesh[i].GetMed()->GetNuclideInTurn(j).GetMicxs().GetData1d(nu).mult(mesh[i].GetMed()->GetNuclideInTurn(j).GetMicxs().GetData1d(sigf));
          src=src+micnusigf*mesh[i].GetMed()->GetNuclideInTurn(j).GetDensity()/rr;
        };
      };
      PutIsotropicSourcePerUnitVolume(i,src);
    };
  };
  if(subcritical){
    //CalFixedSourceWithFission(1e-4,20,true);
    CalFixedSourceWithFission(1e-5,1000,false);
  }else{
    CalGPT(keff);
  };
};

void GeneralSystem::CalGPTForFissionReactionRate(real keff,real rr,bool subcritical)
{
  GroupData1D src(grp);

  SetZeroScatSrc();
  for(int i=0;i<TotM;i++){
    if(mesh[i].GetMed()->IsFissionable()){
      src.set_zero();
      for(int j=0;j<mesh[i].GetMed()->GetNucnum();j++){
	if(mesh[i].GetMed()->GetNuclideInTurn(j).IsFissionable()){
	  GroupData1D micsigf=mesh[i].GetMed()->GetNuclideInTurn(j).GetMicxs().GetData1d(sigf);
          src=src+micsigf*mesh[i].GetMed()->GetNuclideInTurn(j).GetDensity()/rr;
        };
      };
      PutIsotropicSourcePerUnitVolume(i,src);
    };
  };
  if(subcritical){
    CalFixedSourceWithFission(1e-5,1000,false);
  }else{
    CalGPT(keff);
  };
};

void GeneralSystem::CalGPTForFissionSourceDelayed(real keff,real rr,DelayedNeutronData& dnd,bool subcritical)
{
  GroupData1D src(grp);

  SetZeroScatSrc();
  for(int i=0;i<TotM;i++){
    if(mesh[i].GetMed()->IsFissionable()){
      src.set_zero();
      for(int j=0;j<mesh[i].GetMed()->GetNucnum();j++){
	if(mesh[i].GetMed()->GetNuclideInTurn(j).IsFissionable()&&mesh[i].GetMed()->GetNuclideInTurn(j).GetDensity()!=0.){
	  int matno=mesh[i].GetMed()->GetNuclideInTurn(j).GetMatnum();
	  if(dnd.GetNumFromNucid(matno)!=-1){
  	    GroupData1D micnusigf=dnd.GetYieldFromMatNumber(matno).mult(mesh[i].GetMed()->GetNuclideInTurn(j).GetMicxs().GetData1d(sigf));
            src=src+micnusigf*mesh[i].GetMed()->GetNuclideInTurn(j).GetDensity()/rr;
	  };
        };
      };
      PutIsotropicSourcePerUnitVolume(i,src);
    };
  };
  if(subcritical){
    //CalFixedSourceWithFission(1e-4,20,true);
    CalFixedSourceWithFission(1e-5,1000,false);
  }else{
    CalGPT(keff);
  };
};

void GeneralSystem::CalGPT(real keff,real eps_f,int itermax)
{
  if(name!="PLOS"&&name!="PJI"&&name!="MEC"&&name!="SNR"&&name!="SNT"){    
    cout<<"# Error in GeneralSystem::CalGPT.\n";
    cout<<"# This method cannot be used in solver "<<name<<"\n";
    exit(0);
  };

  if(name=="SNR"||name=="SNT"){
    cout<<"# Warning in GeneralSystem::CalGPT.\n";
    cout<<"# Angular flux is NOT touched;\n";
    cout<<"# fundamental mode components are NOT extracted.\n";
  };

  CheckUpScattering();

  vector<GroupData1D> gpt_flx_store(TotM);
  vector<GroupData1D> gpt_flx_old(TotM);
  for(int i=0;i<TotM;i++){
    gpt_flx_store[i].put_imax(grp);
    gpt_flx_old[i].put_imax(grp);
  };

  /*
  vector< vector<GroupData1D> > gpt_aflx_store;
  if(name=="MEC"&&WriteFlux()){
gpt_aflx_store
  };
  */

  int iter=0;
  real sumold=0.;
  real sum=0.;
  for(int k=0;k<itermax;k++){

    sumold=sum;
    sum=0.;
    for(int m=0;m<TotM;m++){
      if(k!=0)mesh[m].CalFissionSrcAdjoint();
      sum+=mesh[m].GetFissionSrc();
    };
    //if(k!=0)cout<<"     "<<sumold<<" -> "<<sum<<" : "<<sum/sumold<<"\n";
    
    if(k==0){
      if(IsUpScatterSystem){
        CalFixedSourceUpScat();
      }else{
        CalFixedSource(1e-4,1,false);
      };
    }else{
      SetZeroScatSrc();
      if(IsUpScatterSystem){
        CalFixedFissionSourceUpScat(keff);
      }else{
        CalFixedFissionSource(keff,k);
      };
    };
    iter++;

    real errmax=0.;
    real flxmax=0.;
    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
        real tmp=mesh[m].GetFlux().get_dat(g);
        real flxold=gpt_flx_old[m].get_dat(g);
	//cout<<m<<" "<<g<<" "<<flxold<<" -> "<<tmp<<" "<<tmp/flxold<<"\n";
        if(flxold!=0.&&k!=0){
	  //if(flxold>0.&&k!=0){
          real err=fabs(tmp/flxold-1.);
          if(err>errmax)errmax=err;
        };
	if(fabs(tmp)>flxmax)flxmax=fabs(tmp);
        gpt_flx_old[m].put_data(g,tmp);
      };
    };

    //cout<<"# CALGPT : "<<k<<" "<<flxmax<<"\n";
    //cout<<"# CALGPT : "<<k<<" "<<errmax<<" (Maximum flux : "<<flxmax<<")\n";

    if(errmax<eps_f&&k!=0){
      break;
    };

    if(k==0){
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m].copy(gpt_flx_old[m]);
      };
    }else{
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m]=gpt_flx_store[m]+gpt_flx_old[m];
      };
    };

    if(k==itermax-1){
      cout<<"# Warning in GeneralSystem::CalGPT.\n";
      cout<<"# Iteration reaches to maximum ("<<itermax<<").\n";
      cout<<"# Iteration residual : "<<errmax<<"\n";
    };
  };

  for(int m=0;m<TotM;m++){
    mesh[m].GetFlux().copy(gpt_flx_store[m]-gpt_flx_old[m]*(iter-1));
  };

  //mesh[0].GetFlux().show_self();

};

// +++++++++++++++++++++
// for source extrapolation
// +++++++++++++++++++++

real GeneralSystem::SourceAndResidualRevision(vector<real> &fsold, vector<real> &res)
{
  real lamda=0.;
  int ind=0;
  for(int i=0;i<TotM;i++){
    if(mesh[i].GetFissionSrc()>0.){
      fsold[i]=mesh[i].GetFissionSrc();
      if(opt.Forward()){mesh[i].CalFissionSrc();}
      else{mesh[i].CalFissionSrcAdjoint();};
      real tmp=mesh[i].GetFissionSrc()-fsold[i];
      real tmp2=tmp/res[i];
      res[i]=tmp;
      if(tmp2>0.&&tmp2<1.){
        lamda+=tmp2;
	ind++;
      };
    };
  };
  lamda/=ind;
  return lamda;
};

void GeneralSystem::SourceExtrapolation(vector<real> &fsold,real omega)
{
  for(int i=0;i<TotM;i++){
    if(fsold[i]>0.){
      real tmp=fsold[i]+omega*(mesh[i].GetFissionSrc()-fsold[i]);
      if(tmp>0.)mesh[i].PutFissionSrc(tmp);
    };
  };
};

real GeneralSystem::GetSourceError(vector<real> &fsold)
{
  real errs=0.;
  for(int i=0;i<TotM;i++){
    real tmp=mesh[i].GetFissionSrc();
    if(tmp>0.){
      real err=fabs(tmp/fsold[i]-1.);
      if(err>errs)errs=err;
    };
  };
  return errs;
};

void GeneralSystem::WriteSourceExtrapolationInfo(real lambda,real omega)
{
  if(print){
    cout<<"#   (Extrapolation for fission source)\n";
    cout.setf(ios::showpoint);
    cout.precision(5);
    cout<<"#     dominance ratio="<<lambda<<" , omega = "<<omega<<"\n";
  };
};


// **********************************************
// * For perturbation calculation               *
// **********************************************

void GeneralSystem::CheckSameMesh(GeneralSystem *sec)
{
  bool error=false;
  if(Ndim!=sec->GetNdim())error=true;
  if(TotM!=sec->GetTotM())error=true;
  for(int i=0;i<Ndim;i++){
    if(name!="PJI"&&name!="MEC"){
      if(mi.GetFMesh(i)!=sec->GetMI().GetFMesh(i))error=true;
    };
  };
  if(error){
    cout<<"Error in GeneralSystem.\n";
    cout<<"There is inconsistency between two systems.\n";
    exit(0);
  };
};

void GeneralSystem::CheckAdjointForward(GeneralSystem *sec)
{
  if(GetGeneralOption().Forward()){
    cout<<"# Error in sensitivity calculation.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"# Error in sensitivity calculation.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };
};

real GeneralSystem::CalPerturbDenominator(GeneralSystem *sec)
{
  real denom=0.;
  for(int r=0;r<TotM;r++){
    mesh[r].CalFissionSrcAdjoint();
    denom+=sec->GetMesh(r).GetFlux(0)*mesh[r].GetFissionSrcSigf();
  };
  return denom;
};

real GeneralSystem::CalPerturbDenominatorWithFissionSpectrumMatrix(GeneralSystem *sec)
{
  real denom=0.;
  for(int r=0;r<TotM;r++){
    real vol=mesh[r].GetVolume();
    for(int g=0;g<grp;g++){
      real adj=mesh[r].GetFlux(0).get_dat(g);
      for(int g2=0;g2<grp;g2++){
	denom+=adj*sec->GetMesh(r).GetMed()->GetMacxs().GetData2d(chi).get_dat(g2,g)
	  *sec->GetMesh(r).GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g2)
          *sec->GetMesh(r).GetFlux(0).get_dat(g2)*vol;
      };
    };
  };
  return denom;
};

real GeneralSystem::CalPerturbDenominatorFromSource(GeneralSystem *sec)
// This is appropriate only for the FOP calculations.
{
  real denom=0.;
  for(int r=0;r<TotM;r++){
    real vol=mesh[r].GetVolume();
    denom+=mesh[r].GetFissionSrc()*sec->GetMesh(r).GetFissionSrc()/vol;
  };
  return denom;
};

real GeneralSystem::CalBetaEffective(GeneralSystem *sec,DelayedNeutronData &dnd,bool print,bool fismat)
{
  if(GetGeneralOption().Forward()){
    cout<<"# Error in GeneralSystem::CalBetaEffective.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"# Error in GeneralSystem::CalBetaEffective.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  CheckSameMesh(sec);
  real dp;
  if(fismat){
    dp=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
  }else{
    dp=CalPerturbDenominator(sec);
  };

  int famnum=dnd.GetFamnum();
  int nucnum=dnd.GetNucnum();

  vector< vector<real> > frac;
  frac.resize(nucnum);
  for(int i=0;i<nucnum;i++){
    frac[i].resize(famnum,0.);
  };

  for(int i=0;i<nucnum;i++){
    int nucid=dnd.GetNucid(i);
    for(int k=0;k<TotM;k++){
      if(mesh[k].GetMed()->ExistNuclide(nucid)){
        real tmp1=(mesh[k].GetMed()->GetNuclide(nucid).GetMicxs().GetData1d(sigf).mult(dnd.GetYield(i)))
	          *sec->GetMesh(k).GetFlux()*mesh[k].GetMed()->GetNuclide(nucid).GetDensity()
                  *mesh[k].GetVolume(); 
        for(int j=0;j<famnum;j++){
 	  frac[i][j]+=tmp1*(dnd.GetKai(i,j)*mesh[k].GetFlux())*dnd.GetFraction(i,j);
	};
      };
    };
  };

  real fsum=0;
  for(int i=0;i<nucnum;i++){
    for(int j=0;j<famnum;j++){
      frac[i][j]/=dp;
      fsum+=frac[i][j];
    };
  };

  if(print){
    cout<<"###############################################################\n";
    cout<<"# Delayed neutron parameter calculation\n#\n";

    cout.setf(ios::scientific);
    cout.precision(6);
    cout<<"# Beta_eff = "<<fsum<<"\n#\n";

    vector<bool> nonzero(nucnum,false);
    cout<<"# [Nuclide-wise (relative contribution)]\n#\n";
    for(int i=0;i<nucnum;i++){
      real sum=0.;
      for(int j=0;j<famnum;j++){
	sum+=frac[i][j];
      };
      int nucid=dnd.GetNucid(i);
      real dat=sum;
      if(dat!=0.){
	nonzero[i]=true;
        cout<<"#  "<<nucid<<" : "<<dat<<"  ("<<dat/fsum<<")\n";
      };
    };

    vector<real> frac_sum(famnum,0.);
    cout<<"#\n# [Nuclide&family-wise]\n#\n";
    cout<<"#           ";
    for(int i=0;i<famnum;i++){
      cout<<i+1<<"         ";
    };
    cout<<"\n";
    for(int i=0;i<nucnum;i++){
      if(nonzero[i]){
        int nucid=dnd.GetNucid(i);
        cout.setf(ios::scientific);
        cout.precision(3);
        cout<<"#  "<<nucid<<" ";
        for(int j=0;j<famnum;j++){
	  cout<<frac[i][j]<<" ";
          frac_sum[j]+=frac[i][j];
        };
        cout<<"\n";
      };
    };
    cout<<"#   (sum) ";
    for(int i=0;i<famnum;i++){
      cout<<frac_sum[i]<<" ";
    };
    cout<<"\n";

    // averaged decay constant calculation
    cout<<"#\n# [Averaged decay constant]\n";
    cout<<"#   (beta-weight)  (beta-weight-inverse)\n";
    for(int i=0;i<famnum;i++){
      real nume1=0.;
      real nume2=0.;
      real denom=0.;
      for(int j=0;j<nucnum;j++){
	nume1+=dnd.GetLambda(j,i)*frac[j][i];
	if(frac[j][i]!=0.)nume2+=frac[j][i]/dnd.GetLambda(j,i);
	denom+=frac[j][i];
      };
      real ave1=nume1/denom;
      real ave2=1./(nume2/denom);
      cout<<"# "<<i+1<<"   "<<ave1<<"        "<<ave2<<"\n";
    };

    cout<<"###############################################################\n";

  };

  return fsum;
};

real GeneralSystem::CalFundamentalBeta(DelayedNeutronData &dnd,bool print)
{
  if(!GetGeneralOption().Forward()){
    cout<<"# Error in GeneralSystem::CalFundamentalBeta.\n";
    cout<<"# Forward flux should be calculated in reference system.\n";
    exit(0);
  };

  GroupData1D yld;
  yld=GetIntegratedReactionRate(nusigf);
  real yldtot=yld.get_sum();

  int nucnum=dnd.GetNucnum();

  vector<real> frac;
  frac.resize(nucnum,0.);

  real tmp=0;
  for(int i=0;i<nucnum;i++){
    int nucid=dnd.GetNucid(i);
    for(int k=0;k<TotM;k++){
      if(mesh[k].GetMed()->ExistNuclide(nucid)){
	real den=mesh[k].GetMed()->GetNuclide(nucid).GetDensity();
        real tmp1=(mesh[k].GetMed()->GetNuclide(nucid).GetMicxs().GetData1d(sigf).mult(dnd.GetYield(i)))
	          *mesh[k].GetFlux();
        real vol=mesh[k].GetVolume();
	frac[i]+=vol*den*tmp1;
        tmp+=vol*den*tmp1;
      };
    };
  };

  if(print){
    cout<<"###############################################################\n";
    cout<<"# Delayed neutron parameter calculation\n#\n";

    cout.setf(ios::scientific);
    cout.precision(6);
    cout<<"# Fundamental Beta = "<<tmp/yldtot<<"\n#\n";

    cout<<"# [Nuclide-wise (relative contribution)]\n#\n";
    for(int i=0;i<nucnum;i++){
      int nucid=dnd.GetNucid(i);
      if(frac[i]!=0.){
        cout<<"#  "<<nucid<<" : "<<frac[i]/yldtot<<"  ("<<frac[i]/tmp<<")\n";
      };
    };
    cout<<"############################################################\n";
  };

  return tmp/yldtot;
};

real GeneralSystem::CalNeutronLifeTime(GeneralSystem *sec,bool fismat)
{
  if(GetGeneralOption().Forward()){
    cout<<"# Error in GeneralSystem::CalNeutronLifeTime.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"# Error in GeneralSystem::CalNeutronLifeTime.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  CheckSameMesh(sec);
  real dp;
  if(fismat){
    dp=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
  }else{
    dp=CalPerturbDenominator(sec);
  };

  real m_n=1.674927e-27; // [kg]
  real J_eV=1.602e-19; // [J/eV]
  real sum=0.;
  for(int k=0;k<TotM;k++){
    real vol=mesh[k].GetVolume();
    for(int g=0;g<grp-1;g++){ // contribution of the final energy group is igonred
      real e_ave=sqrt(med[0].GetEnband().get_dat(g)*med[0].GetEnband().get_dat(g+1));
      real v_ave=sqrt(2.*e_ave*J_eV/m_n);// [m/s]
      v_ave*=100.; // [cm/s]
      sum+=mesh[k].GetFlux().get_dat(g)*sec->GetMesh(k).GetFlux().get_dat(g)*vol/v_ave;
    };
  };

  real l=sum/dp;
  cout<<"###############################################################\n";
  cout<<"# Neutron generation time : "<<l<<" [s]\n";
  cout<<"###############################################################\n";

  return l;
};

void GeneralSystem::CalPerturbYieldTerm(GeneralSystem *sec,bool *flag,real *val,real kp)
{
  real inv_kp=1./kp;
  GroupData1D dnusigf(grp);
  GroupData1D dkai(grp);
  for(int i=0;i<grp;i++){val[i]=0.;};
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      // fission term
      dnusigf=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf()
  	     -mesh[m].GetMed()->GetMacxs().GetNusigf();
      real a1=mesh[m].GetMed()->GetMacxs().GetKai()*mesh[m].GetFlux();
      a1*=vol*inv_kp;
      GroupData1D tmp=dnusigf.mult(sec->GetMesh(m).GetFlux(0));
      for(int g=0;g<grp;g++){
        val[g]+=a1*tmp.get_dat(g);
      };
      // fission spectrum term
      dkai=sec->GetMesh(m).GetMed()->GetMacxs().GetKai()
  	          -mesh[m].GetMed()->GetMacxs().GetKai();
      real a2=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf()*sec->GetMesh(m).GetFlux();
      tmp=dkai.mult(mesh[m].GetFlux(0));
      a2*=vol*inv_kp;
      for(int g=0;g<grp;g++){
        val[g]+=a2*tmp.get_dat(g);
      };
    };
  };
};

void GeneralSystem::CalPerturbAbsorptionTerm(GeneralSystem *sec,bool *flag,real *val)
{
  GroupData1D dsiga(grp);
  for(int i=0;i<grp;i++){val[i]=0.;};
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      dsiga=sec->GetMesh(m).GetMed()->GetMacxs().GetSiga()
  	           -mesh[m].GetMed()->GetMacxs().GetSiga();
      GroupData1D bil_flx=mesh[m].GetFlux(0).mult(sec->GetMesh(m).GetFlux(0));
      real vol=mesh[m].GetVolume();
      for(int g=0;g<grp;g++){
        val[g]+=-1.*bil_flx.get_dat(g)*dsiga.get_dat(g)*vol;
      };
    };
  };
};

real GeneralSystem::CalPerturbAbsorptionTerm(GeneralSystem *sec,bool *flag,int g)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      ret+=CalPerturbAbsorptionTerm(sec,m,g);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbAbsorptionTerm(GeneralSystem *sec,int m,int g)
{
  real dsiga=sec->GetMesh(m).GetMed()->GetMacxs().GetSiga().get_dat(g)
                    -mesh[m].GetMed()->GetMacxs().GetSiga().get_dat(g);
  real vol=mesh[m].GetVolume();
  real ret=-1.*mesh[m].GetFlux(0).get_dat(g)*dsiga*sec->GetMesh(m).GetFlux(0).get_dat(g)*vol;
  return ret;
};


void GeneralSystem::CalPerturbTotalTerm(GeneralSystem *sec,bool *flag,real *val)
{
  GroupData1D dsiga(grp);
  for(int i=0;i<grp;i++){val[i]=0.;};
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      dsiga=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt()
  	           -mesh[m].GetMed()->GetMacxs().GetSigt();
      GroupData1D bil_flx=mesh[m].GetFlux(0).mult(sec->GetMesh(m).GetFlux(0));
      real vol=mesh[m].GetVolume();
      for(int g=0;g<grp;g++){
        val[g]+=-1.*bil_flx.get_dat(g)*dsiga.get_dat(g)*vol;
      };
    };
  };
};

real GeneralSystem::CalPerturbTotalTerm(GeneralSystem *sec,bool *flag,int g)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      ret+=CalPerturbTotalTerm(sec,m,g);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbTotalTerm(GeneralSystem *sec,int m,int g)
{
  real dsiga=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt().get_dat(g)
                    -mesh[m].GetMed()->GetMacxs().GetSigt().get_dat(g);
  real vol=mesh[m].GetVolume();
  real ret=-1.*mesh[m].GetFlux(0).get_dat(g)*dsiga*sec->GetMesh(m).GetFlux(0).get_dat(g)*vol;
  return ret;
};



real GeneralSystem::CalPerturbYieldTerm(GeneralSystem *sec,bool *flag,int g,real kp)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      ret+=CalPerturbYieldTerm(sec,m,g,kp);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbYieldTerm(GeneralSystem *sec,int m,int g,real kp)
{
  real ret=0.;
  real inv_kp=1./kp;

  real vol=mesh[m].GetVolume();
  // fission term
  real dnusigf=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf().get_dat(g)
	              -mesh[m].GetMed()->GetMacxs().GetNusigf().get_dat(g);
  real a1=mesh[m].GetMed()->GetMacxs().GetKai()*mesh[m].GetFlux();
  a1*=vol*inv_kp;
  ret+=a1*dnusigf*sec->GetMesh(m).GetFlux(0).get_dat(g);
  // fission spectrum term
  real dkai=sec->GetMesh(m).GetMed()->GetMacxs().GetKai().get_dat(g)
                   -mesh[m].GetMed()->GetMacxs().GetKai().get_dat(g);
  real a2=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf()*sec->GetMesh(m).GetFlux();
  a2*=vol*inv_kp;
  ret+=a2*dkai*mesh[m].GetFlux(0).get_dat(g);

  return ret;
};

real GeneralSystem::CalPerturbYieldTermWithFissionSpectrumMatrix(GeneralSystem *sec,bool *flag,int g,real kp)
{
  real ret=0.;
  real inv_kp=1./kp;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      // fission term
      real dnusigf=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf().get_dat(g)
   	                  -mesh[m].GetMed()->GetMacxs().GetNusigf().get_dat(g);
      real a1=0.;
      for(int g2=0;g2<grp;g2++){
	a1+=mesh[m].GetMed()->GetMacxs().GetData2d(chi).get_dat(g,g2)
	  *mesh[m].GetFlux().get_dat(g2);
      };
      a1*=vol*inv_kp;
      ret+=a1*dnusigf*sec->GetMesh(m).GetFlux(0).get_dat(g);
      // fission spectrum term
      real a2=0.;
      /*
      for(int g2=0;g2<grp;g2++){
        real dkai=sec->GetMesh(m).GetMed()->GetMacxs().GetData2d(chi).get_dat(g2,g)
	                 -mesh[m].GetMed()->GetMacxs().GetData2d(chi).get_dat(g2,g);
        a2+=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf().get_dat(g2)
	   *sec->GetMesh(m).GetFlux().get_dat(g2)*dkai;
      };
      a2*=vol*inv_kp;
      ret+=a2*mesh[m].GetFlux(0).get_dat(g);
      */
      real nsf=sec->GetMesh(m).GetMed()->GetMacxs().GetNusigf().get_dat(g);
      for(int g2=0;g2<grp;g2++){
        real dkai=sec->GetMesh(m).GetMed()->GetMacxs().GetData2d(chi).get_dat(g,g2)
	                 -mesh[m].GetMed()->GetMacxs().GetData2d(chi).get_dat(g,g2);
        a2+=nsf*dkai*mesh[m].GetFlux().get_dat(g2);
      };
      a2*=vol*inv_kp;
      ret+=a2*sec->GetMesh(m).GetFlux(0).get_dat(g);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbN2NTerm(GeneralSystem *sec,bool *flag,int g)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      ret+=CalPerturbN2NTerm(sec,m,g);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbN2NTerm(GeneralSystem *sec,int m,int g)
{
  real dsign2n=sec->GetMesh(m).GetMed()->GetMacxs().GetSign2n().get_dat(g)
                      -mesh[m].GetMed()->GetMacxs().GetSign2n().get_dat(g);
  real vol=mesh[m].GetVolume();
  real ret=mesh[m].GetFlux(0).get_dat(g)*dsign2n*sec->GetMesh(m).GetFlux(0).get_dat(g)*vol;
  return ret;
};

real GeneralSystem::CalPerturbScatteringTermDiffusion(GeneralSystem *sec,bool *flag,int g)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      ret+=CalPerturbScatteringTermDiffusion(sec,m,g);
    };
  };
  return ret;
};

real GeneralSystem::CalPerturbScatteringTermDiffusion(GeneralSystem *sec,int m,int g)
{
  real f2v=sec->GetMesh(m).GetFlux(0).get_dat(g)*mesh[m].GetVolume();
  real f1=mesh[m].GetFlux(0).get_dat(g);
  real tmp=0.;
  for(int j=0;j<grp;j++){ // w up scattering
    //for(int j=g+1;j<grp;j++){ // w/o up scattering
    tmp+=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g,j)
                 -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,j))*
         (f1-mesh[m].GetFlux(0).get_dat(j));
  };
  real ret=-tmp*f2v;
  return ret;
};

real GeneralSystem::CalPerturbScatteringTermDiffusion(GeneralSystem *sec,bool *flag,int g1,int g2)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real tmp=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g1,g2)
 	               -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g1,g2))*
	     (mesh[m].GetFlux(0).get_dat(g1)-mesh[m].GetFlux(0).get_dat(g2));
      ret-=tmp*sec->GetMesh(m).GetFlux(0).get_dat(g1)*mesh[m].GetVolume();
    };
  };
  return ret;
};

void GeneralSystem::PutPerturbationDenominator(real val)
{
  ip_input=true;
  ip_input_value=val;
};

SensitivityData GeneralSystem::CalDirectTermPointReactionRateSigf(int meshid, int matid)
{
  SensitivityData sens_dir;
  sens_dir.PutGroup(grp);
  GroupData1D snsdat(grp);
  snsdat=GetMesh(meshid).GetFlux().mult(GetMesh(meshid).GetMed()->GetNuclide(matid).GetMicxs().GetData1d(sigf));
  snsdat=snsdat*(1./snsdat.get_sum());
  sens_dir.PutSensitivity1D(matid,18,snsdat);
  return sens_dir;
};

GroupData1D GeneralSystem::GetPositionFlux(real x,real y,real z)
{
  // Linearly interpolation

  GroupData1D ret(grp);

  int xpos=mi.GetXMeshPosition(x);
  int ypos=mi.GetYMeshPosition(y);
  int zpos=mi.GetZMeshPosition(z);

  int mid=meshid[zpos][ypos][xpos];

  for(int i=0;i<grp;i++){
    real flx_cent=mesh[mid].GetFlux().get_dat(i);

    real xcor=1.;
    real xcent=mi.GetXMeshLocation(xpos);
    if(xcent>x&&xpos!=0){
      real flx2=mesh[meshid[zpos][ypos][xpos-1]].GetFlux().get_dat(i);
      real x2=mi.GetXMeshLocation(xpos-1);
      real xnew=flx_cent+(x-xcent)/(x2-xcent)*(flx2-flx_cent);
      xcor=xnew/flx_cent;
    };
    if(xcent<x&&xpos!=mi.GetXF()-1){
      real flx2=mesh[meshid[zpos][ypos][xpos+1]].GetFlux().get_dat(i);
      real x2=mi.GetXMeshLocation(xpos+1);
      real xnew=flx_cent+(x-xcent)/(x2-xcent)*(flx2-flx_cent);
      xcor=xnew/flx_cent;
    };

    real ycor=1.;
    real ycent=mi.GetYMeshLocation(ypos);
    if(ycent>y&&ypos!=0){
      real flx2=mesh[meshid[zpos][ypos-1][xpos]].GetFlux().get_dat(i);
      real x2=mi.GetYMeshLocation(ypos-1);
      real xnew=flx_cent+(y-ycent)/(x2-ycent)*(flx2-flx_cent);
      ycor=xnew/flx_cent;
    };
    if(ycent<y&&ypos!=mi.GetYF()-1){
      real flx2=mesh[meshid[zpos][ypos+1][xpos]].GetFlux().get_dat(i);
      real x2=mi.GetYMeshLocation(ypos+1);
      real xnew=flx_cent+(y-ycent)/(x2-ycent)*(flx2-flx_cent);
      ycor=xnew/flx_cent;
    };

    real zcor=1.;
    real zcent=mi.GetZMeshLocation(zpos);
    if(zcent>z&&zpos!=0){
      real flx2=mesh[meshid[zpos-1][ypos][xpos]].GetFlux().get_dat(i);
      real x2=mi.GetZMeshLocation(zpos-1);
      real xnew=flx_cent+(z-zcent)/(x2-zcent)*(flx2-flx_cent);
      zcor=xnew/flx_cent;
    };
    if(zcent<z&&zpos!=mi.GetZF()-1){
      real flx2=mesh[meshid[zpos+1][ypos][xpos]].GetFlux().get_dat(i);
      real x2=mi.GetZMeshLocation(zpos+1);
      real xnew=flx_cent+(z-zcent)/(x2-zcent)*(flx2-flx_cent);
      zcor=xnew/flx_cent;
    };

    ret.put_data(i,flx_cent*xcor*ycor*zcor);
  };

  //ret.copy(mesh[meshid[zpos][ypos][xpos]].GetFlux());

  return ret;
};

real GeneralSystem::CalPositionReactionRate(real x,real y,real z,int matnum,enum xstype react)
{
  int xpos=mi.GetXMeshPosition(x);
  int ypos=mi.GetYMeshPosition(y);
  int zpos=mi.GetZMeshPosition(z);
  int mid=meshid[zpos][ypos][xpos];

  GroupData1D flx=GetPositionFlux(x,y,z);

  return mesh[mid].GetMed()->GetNuclide(matnum).GetMicxs().GetData1d(react)*flx;
};

Medium GeneralSystem::Homogenize(int hreg, int *nreg)
{
  Medium ret(grp);
  int pl=med[0].GetPL();
  int pltot=med[0].GetPLT();
  ret.PutPL(pl,pltot);
  ret.GetEnband().copy(med[0].GetEnband());

  GroupData1D totf(grp),totc(grp); // total flux & current
  totf.set_zero();
  totc.set_zero();

  vector<GroupData1D> cur;
  cur.resize(TotM);
  for(int i=0;i<TotM;i++){
    cur[i].put_imax(grp);
    cur[i]=mesh[i].GetFlux(0).mult(mesh[i].GetMed()->GetMacxs().GetData1d(d));
    // Current is defined as D*Flux
  };

  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    totf=totf+mesh[r].GetFlux()*mesh[r].GetVolume();
    totc=totc+cur[r]*mesh[r].GetVolume();
  };

  // 1D flux-weighted cross section
  enum xstype ss[]={siga,nusigf,sigt,sign2n};
  for(int i=0;i<4;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetData1D(ss[i]).mult(mesh[r].GetFlux())
           *mesh[r].GetVolume();
    };
    rr=rr/totf;
    ret.GetData1D(ss[i]).copy(rr);
  };

  // 1D current-weighted cross section
  for(int i=0;i<1;i++){
    GroupData1D rr(grp);
    enum xstype ss;
    if(i==0)ss=sigtr;
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetData1D(ss).mult(cur[r])
           *mesh[r].GetVolume();
    };
    rr=rr/totc;
    ret.GetData1D(ss).copy(rr);
  };

  // total cross section
  for(int i=1;i<=pltot;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetSigt(i).mult(cur[r])
           *mesh[r].GetVolume();
    };
    rr=rr/totc;
    ret.GetSigt(i).copy(rr);
  };

  // Scattering cross section
  GroupData2D mat(grp,grp);
  for(int l=0;l<=pl;l++){
    mat.set_zero();
    for(int i=0;i<grp;i++){
      for(int j=0;j<grp;j++){
        real bs=0.0;
        for(int k=0;k<hreg;k++){
    	  int r=nreg[k];
          real tmp=mesh[r].GetMed()->GetDataSigs(l,i,j)*mesh[r].GetVolume();
	  if(l==0)tmp*=mesh[r].GetFlux().get_dat(i);
	  if(l!=0)tmp*=cur[r].get_dat(i);
	  bs+=tmp;
        };
        if(l==0)bs/=totf.get_dat(i);
        if(l!=0)bs/=totc.get_dat(i);
        mat.put_data(i,j,bs);
      };
    };
    ret.GetSigs(l).copy(mat);
  };

  // Fission spectrum
  real *fiss=new real[hreg];
  real totfis=0.0;
  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    fiss[i]=mesh[r].GetMed()->GetData1D(nusigf)*mesh[r].GetFlux()
           *mesh[r].GetVolume();
    totfis+=fiss[i];
  };
  if(totfis>0.){
    totfis=1./totfis;
    real kais;
    for(int i=0;i<grp;i++){
      kais=0.0;
      for(int j=0;j<hreg;j++){
        int r=nreg[j];
        kais+=mesh[r].GetMed()->GetData1D(chi).get_dat(i)*fiss[j];
      };
      kais*=totfis;
      ret.GetData1D(chi).put_data(i,kais);
    };
  };
  delete []fiss;

  // Fission spectrum vector
  GroupData2D inp_vec(grp,grp);
  inp_vec.set_zero();
  for(int i=0;i<grp;i++){
    for(int j=0;j<grp;j++){
      for(int k=0;k<hreg;k++){
	int r=nreg[k];
        real nsf=mesh[r].GetMed()->GetData1D(nusigf).get_dat(i);
	if(nsf!=0.){
	  inp_vec.add_data
          (i,j,nsf*mesh[r].GetVolume()*mesh[r].GetFlux().get_dat(i)*mesh[r].GetMed()->GetMacxs().GetData2d(chi).get_dat(i,j));
	};
      };
    };
  };
  for(int i=0;i<grp;i++){
    real tot=0.;
    for(int j=0;j<grp;j++){
      tot+=inp_vec.get_dat(i,j);
    };
    if(tot!=0.){
      tot=1./tot;
      for(int j=0;j<grp;j++){
	real org=inp_vec.get_dat(i,j);
	inp_vec.put_data(i,j,org*tot);
      };
    };
  };
  ret.GetMacxs().GetData2d(chi).copy(inp_vec);

  ret.CalSigtr();
  // Diffusion coefficient
  ret.CalDFromSigtr();
  ret.GetMacxs().GetData1d(dr).copy(ret.GetMacxs().GetData1d(d));
  ret.GetMacxs().GetData1d(dz).copy(ret.GetMacxs().GetData1d(d));

  return ret;
}

real GeneralSystem::GetNeutronMultiplicationInfo(real keff)
{
  GroupData1D abs;
  GroupData1D yld;
  GroupData1D n2n;
  abs=GetIntegratedReactionRate(siga);
  yld=GetIntegratedReactionRate(nusigf);
  n2n=GetIntegratedReactionRate(sign2n);
  real fisr=0.;
  for(int i=0;i<TotM;i++){
    real vol=GetMesh(i).GetVolume();
    fisr+=GetMesh(i).GetFlux()*GetMesh(i).GetMed()->GetMacroSigf()*vol;
  };
  real abstot=abs.get_sum();
  real yldtot=yld.get_sum();
  real n2ntot=n2n.get_sum()*0.5;
  abstot/=yldtot;
  n2ntot/=yldtot;
  fisr=yldtot/fisr;
  yldtot=1.;
  real leaktot=yldtot/keff+n2ntot-abstot;

  real ff=yldtot/keff;

  cout.setf(ios::showpoint);
  cout.precision(7);
  cout<<"#-----------------------------------------------------------\n";
  cout<<"# Keff contribution normalized by fission yield\n#\n";
  cout<<"#  Total absorption   : "<<abstot/ff<<"\n";
  cout<<"#  Total n2n yield    : "<<n2ntot/ff<<"\n";
  cout<<"#  Total leakage      : "<<leaktot/ff<<"\n";
  cout<<"#\n#     (Averaged nu-value : "<<fisr<<")\n";
  cout<<"#-----------------------------------------------------------\n";

  // Yield fraction calculation
  real yld_fast=0.;
  real yld_inter=0.;
  real yld_thermal=0.;
  real e_fi=1e5;
  real e_it=0.625;
  for(int i=0;i<grp;i++){
    real etop=med[0].GetEnband().get_dat(i);
    real ebot=med[0].GetEnband().get_dat(i+1);
    real letwid=log(etop/ebot);
    real yldg=yld.get_dat(i);
    if(etop>e_fi){
      if(ebot>e_fi){
        yld_fast+=yldg;
      }else{
	real letwid1=log(etop/e_fi);
	yld_fast+=yldg*letwid1/letwid;
        yld_inter+=yldg*(letwid-letwid1)/letwid;
      };
    }else if(etop<=e_fi&&etop>e_it){
      if(ebot>e_it){
	yld_inter+=yldg;
      }else{
	real letwid1=log(etop/e_it);
	yld_inter+=yldg*letwid1/letwid;
	yld_thermal+=yldg*(letwid-letwid1)/letwid;
      };
    }else{
      yld_thermal+=yldg;
    };
  };
  real yld_sum=yld.get_sum();
  cout<<"# Energy contribution of fission neutrons\n#\n";
  cout<<"#   Fast range         (>100keV)  : "<<yld_fast/yld_sum<<"\n";
  cout<<"#   Intermediate range (>0.625eV) : "<<yld_inter/yld_sum<<"\n";
  cout<<"#   Thermal range      (<0.625eV) : "<<yld_thermal/yld_sum<<"\n#\n";

  yld_thermal/=yld_sum;

  int hm_num=1+1+6+1+5+3+5;
  int hm_id[]={
    902320,
    912330,
    922330,922340,922350,922360,922370,922380,
    932370,
    942380,942390,942400,942410,942420,
    952410,952421,952430,
    962420,962430,962440,962450,962460
  };
  string hm_name[]={
    "Th232 ",
    "Pa233 ",
    "U233  ","U234  ","U235  ","U236  ","U237  ","U238  ",
    "Np237 ",
    "Pu238 ","Pu239 ","Pu240 ","Pu241 ","Pu242 ",
    "Am241 ","Am242m","Am243 ",
    "Cm242 ","Cm243 ","Cm244 ","Cm245 ","Cm246 "
  };
  real *hm_fis=new real[hm_num];
  real *hm_abs=new real[hm_num];
  for(int i=0;i<hm_num;i++){
    hm_fis[i]=0.;
    hm_abs[i]=0.;
  };

  GroupData1D fis(grp);
  GroupData1D abs2(grp);
  real nume=0.;
  real deno=0.;
  real nume2=0.;
  real deno2=0.;
  real labs=0.;
  real mabs=0.;
  real habs=0.;
  for(int i=0;i<grp;i++){
    real mede=sqrt(med[0].GetEnband().get_dat(i)*med[0].GetEnband().get_dat(i+1));
    real medl=log(2e7/mede);
    real sum_fis=0.;
    real sum_abs=0.;
    for(int j=0;j<TotM;j++){
      real vol=mesh[j].GetVolume();
      real volflx=vol*mesh[j].GetFlux().get_dat(i);
      nume2+=medl*volflx;
      deno2+=volflx;
      int nucnum=mesh[j].GetMed()->GetNucnum();
      for(int k=0;k<nucnum;k++){
	real sf=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigf).get_dat(i);
	real sc=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigc).get_dat(i);
        real den=mesh[j].GetMed()->GetNuclideInTurn(k).GetDensity();
        real tmp=den*sf*volflx;
        real tmp2=den*sc*volflx;
        sum_fis+=tmp;
        sum_abs+=tmp+tmp2;
        nume+=medl*tmp;
        deno+=tmp;
        int id=mesh[j].GetMed()->GetNuclideInTurn(k).GetMatnum();
	bool check=false;
	for(int l=0;l<hm_num;l++){
	  if(id==hm_id[l]){
	    check=true;
	    hm_fis[l]+=tmp;
	    hm_abs[l]+=tmp+tmp2;
	  };
	};
	if(!check){
	  if(id<100000){labs+=tmp+tmp2;}
	  else if(id<900000){mabs+=tmp+tmp2;}
	  else{habs+=tmp+tmp2;};
	};
      };
    };
    fis.put_data(i,sum_fis);
    abs2.put_data(i,sum_abs);
  };

  real totfis=fis.get_sum();
  real totabs=abs2.get_sum();
  real sum=0.;
  int ggg=-1;
  real ratio=1.;
  for(int i=0;i<grp;i++){
    real tmp=fis.get_dat(i);
    sum+=tmp;
    if(sum>totfis*0.5&&ggg==-1){
      ratio=(tmp-(sum-totfis*0.5))/tmp;
      ggg=i;
    };
  };

  real etop=med[0].GetEnband().get_dat(ggg);
  real ebtm=med[0].GetEnband().get_dat(ggg+1);
  real dlet=log(etop/ebtm);
  cout<<"#   AVG fission Eng. : "<<etop*exp(-dlet*ratio)<<" [eV]\n";

  real ealf=2e7*exp(-nume/deno);
  cout<<"#   Eng. Corr. to the AVG neutron lethargy causing fission : "<<ealf<<" [eV]\n";
  real eavefl=2e7*exp(-nume2/deno2);
  cout<<"#   Eng. Corr. to the AVG neutron lethargy : "<<eavefl<<" [eV]\n";

  cout<<"#-----------------------------------------------------------\n";
  cout<<"# Nuclide-wise fission contribution\n#\n";
  for(int i=0;i<hm_num;i++){
    if(hm_fis[i]>0.)cout<<"#   "<<hm_name[i]<<" "<<hm_fis[i]/totfis<<"\n";
  };

  cout<<"#-----------------------------------------------------------\n";
  cout<<"# Nuclide-wise absorption contribution\n#\n";
  for(int i=0;i<hm_num;i++){
    if(hm_abs[i]>0.)cout<<"#   "<<hm_name[i]<<" "<<hm_abs[i]/totabs<<"\n";
  };
  if(labs>0.)     cout<<"#   light elements  (Z<10)  "<<labs/totabs<<"\n";
  if(mabs>0.)     cout<<"#   medium elements (Z<90)  "<<mabs/totabs<<"\n";
  if(habs>0.)     cout<<"#   heavy elements  (Z>=90) "<<habs/totabs<<"\n";
  cout<<"#-----------------------------------------------------------\n";

  delete [] hm_fis;
  delete [] hm_abs;
  return yld_thermal;
};

void GeneralSystem::GetNeutronMultiplicationDetailedInfo()
{
  real fe=0.;
  real cr=0.;
  real ni=0.;
  real na=0.;
  real mn=0.;
  real u8=0.;
  real u8f=0.;
  real u5=0.;
  real u5f=0.;
  real pu=0.;
  real puf=0.;
  real h=0.;
  real o=0.;
  real n14=0.;
  real f19=0.;
  real zr=0.;

  GroupData1D tmp(grp);
  tmp.set_zero();

  real sum=0.;
  for(int i=0;i<grp;i++){
    for(int j=0;j<TotM;j++){
      real volflx=mesh[j].GetVolume()*mesh[j].GetFlux().get_dat(i);
      int nucnum=mesh[j].GetMed()->GetNucnum();
      for(int k=0;k<nucnum;k++){
        real den=mesh[j].GetMed()->GetNuclideInTurn(k).GetDensity();
	real cc=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigc).get_dat(i);
	real ff=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigf).get_dat(i);
        real rcc=den*cc*volflx;
	real rff=den*ff*volflx;
        int id=mesh[j].GetMed()->GetNuclideInTurn(k).GetMatnum();
	sum+=rcc+rff;
	if(id==922350){u5+=rcc+rff; u5f+=rff;};
	if(id==922380){u8+=rcc+rff; u8f+=rff;};
	if(id>=260000&&id<270000){
          fe+=rcc+rff;
	  tmp.add_data(i,rcc+rff);
        };
	if(id>=940000&&id<950000){pu+=rcc+rff; puf+=rff;};
	if(id>=240000&&id<250000){cr+=rcc+rff;};
	if(id>=280000&&id<290000){ni+=rcc+rff;};
	if(id>=400000&&id<410000){zr+=rcc+rff;};
	if(id==250550){mn+=rcc+rff;};
	if(id==110230){na+=rcc+rff;};
	if(id==10010){h+=rcc+rff;};
	if(id==80160){o+=rcc+rff;};
	if(id==90190){f19+=rcc+rff;};
	if(id==70140){n14+=rcc+rff;};
      };
    };
  };

  cout<<"#\n# ** Absorption fraction (fission) **\n";
  cout<<"#   Na : "<<na/sum<<"\n";
  cout<<"#   Fe : "<<fe/sum<<"\n";
  cout<<"#   Cr : "<<cr/sum<<"\n";
  cout<<"#   Ni : "<<ni/sum<<"\n";
  cout<<"#   Mn : "<<mn/sum<<"\n";
  cout<<"#   Zr : "<<zr/sum<<"\n";
  cout<<"#   U8 : "<<u8/sum<<" ("<<u8f/sum<<")\n";
  cout<<"#   U5 : "<<u5/sum<<" ("<<u5f/sum<<")\n";
  cout<<"#   Pu : "<<pu/sum<<" ("<<puf/sum<<")\n";
  cout<<"#   H  : "<<h/sum<<"\n";
  cout<<"#   O  : "<<o/sum<<"\n";
  cout<<"#   N14: "<<n14/sum<<"\n";
  cout<<"#   F19: "<<f19/sum<<"\n";

  /*
  for(int i=0;i<grp;i++){
    real e0=med[0].GetEnband().get_dat(i);
    real e1=med[0].GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<tmp.get_dat(i)/letwid<<"\n";
  };
  */
};

void GeneralSystem::GetNuclideWiseAbsorptionRate(int nuc,int *nucid)
{
  vector<real> arr(nuc,0.);
  GroupData1D tmp(grp);
  tmp.set_zero();

  real sum=0.;
  for(int i=0;i<grp;i++){
    for(int j=0;j<TotM;j++){
      real volflx=mesh[j].GetVolume()*mesh[j].GetFlux().get_dat(i);
      int nucnum=mesh[j].GetMed()->GetNucnum();
      for(int k=0;k<nucnum;k++){
        real den=mesh[j].GetMed()->GetNuclideInTurn(k).GetDensity();
	real cc=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigc).get_dat(i);
	real ff=mesh[j].GetMed()->GetNuclideInTurn(k).GetMicxs().GetData1d(sigf).get_dat(i);
        real rcc=den*cc*volflx;
	real rff=den*ff*volflx;
        int id=mesh[j].GetMed()->GetNuclideInTurn(k).GetMatnum();
        sum+=rcc+rff;
	for(int ii=0;ii<nuc;ii++){
  	  if(id==nucid[ii]){arr[ii]+=rcc+rff;};
          tmp.add_data(i,rcc+rff);
	};
      };
    };
  };

  cout<<"# ** Absorption fraction **\n";
  for(int i=0;i<nuc;i++){
    cout<<"#  "<<nucid[i]<<" : "<<arr[i]/sum<<"\n";
  };
  /*
  for(int i=0;i<grp;i++){
    real e0=med[0].GetEnband().get_dat(i);
    real e1=med[0].GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<tmp.get_dat(i)/letwid*0.25<<"\n";
  };
  */
};

void GeneralSystem::AllVectorClear()
{
  mesh.clear();
  med.clear();
  meshid.clear();
  cmeshid.clear();
  xedgel.clear();
  xedger.clear();
  yedgel.clear();
  yedger.clear();
  edge.clear();
};

void GeneralSystem::FluxNormalization(real factor)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<mesh[i].GetMaxPL();l++){
      mesh[i].GetFlux(l)=mesh[i].GetFlux(l)*factor;
    };
  };
};

void GeneralSystem::AddFlux(GeneralSystem &sec)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)+sec.GetMesh(i).GetFlux(l);
      mesh[i].GetFlux(l).copy(tmp);
    };
  };
};

void GeneralSystem::NegFlux(GeneralSystem &sec,real factor)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)-sec.GetMesh(i).GetFlux(l)*factor;
      mesh[i].GetFlux(l).copy(tmp);
    };
  };
};

void GeneralSystem::MedClear()
{
  //med.clear();
  /*
  for(int i=0;i<nmed;i++){
    med[i].MacxsVectorClear();
    med[i].NuclideClear();
  };
  */

  vector<Medium>().swap(med);

  nmed=0;
  med.resize(MAX_MED);
};

int GeneralSystem::GetTotalInnerIterationSum()
{
  int ret=0.;
  for(int i=0;i<grp;i++){
    ret+=total_inner_iteration[i];
  };
  return ret;
};
