#include <cstdlib>
#include <iostream>
#include "GeneralMesh.h"

using namespace std;

GeneralMesh::GeneralMesh()
{
  pl = -1;
  grp = -1;
  dim = 0;
  FSrc=0.;
};

GeneralMesh::~GeneralMesh()
{
};

void GeneralMesh::PutDim(int n)
{
  if(n<1||n>3){
    cout<<"Error in 'PutDim' of GeneralMesh.\n";
    exit(0);
  };
  dim=n;
}

void GeneralMesh::PutDim(int n,real *l)
{
  if(n<1||n>3){
    cout<<"Error in 'PutDim' of GeneralMesh.\n";
    exit(0);
  };
  dim=n;

  Volume=1.;
  for(int i=0;i<3;i++){
    if(i<dim){Len[i]=l[i];}
    else{Len[i]=1.;};
    Volume*=Len[i];
  };
  for(int i=0;i<3;i++){
    real ret=1.;
    for(int j=0;j<3;j++){
      if(j!=i)ret*=Len[j];
    };
    Sur[i*2]=ret;
    Sur[i*2+1]=ret;
  };
}

void GeneralMesh::PutDimSphere(real l,real ll,real lr)
{
  dim=1;
  Volume=0.33333333*PI4*(pow(lr,3)-pow(ll,3));
  Len[0]=l;
  Sur[1]=PI4*lr*lr;
  Sur[0]=PI4*ll*ll;
}

void GeneralMesh::PutDimCylinder(real l,real ll,real lr,real z)
{
  dim=2;
  Volume=PI*(lr*lr-ll*ll)*z;
  Len[0]=l;
  Len[1]=z;
  Sur[0]=PI2*ll*z;
  Sur[1]=PI2*lr*z;
  Sur[2]=PI*(lr*lr-ll*ll);
  Sur[3]=Sur[2];
}

void GeneralMesh::SetVectorSizeForFlux()
{
  if(pl!=-1&&grp!=-1){
    Flux.resize(maxpl);
    SSrc.resize(maxpl);
    srcin.resize(maxpl);
    for(int i=0;i<maxpl;i++){
      Flux[i].put_imax(grp);
      SSrc[i].put_imax(grp);
    };
  };
};

void GeneralMesh::PutGrp(int g)
{
  grp=g;
  SetVectorSizeForFlux();
};

void GeneralMesh::PutPL(int l)
{
  if(dim==0){
    cout<<"Do 'PutDim' before 'PutPL' for GeneralMesh.\n";
    exit(0);
  };
  pl=l;
  if(dim==1)maxpl=pl+1;
  if(dim==2||dim==3){
    int ii=0;
    for(int i=0;i<=pl;i++){
      if(dim==2)ii+=i+1; 
      if(dim==3)ii+=i*2+1;
    };
    maxpl=ii;
  };
  SetVectorSizeForFlux();
};

void GeneralMesh::SetZeroScatSrc(int g)
{
  for(int l=0;l<maxpl;l++){
    SSrc[l].put_data(g,0.);
  };
};

void GeneralMesh::SetZeroScalarFlux(int g)
{
  Flux[0].put_data(g,0.);
};

void GeneralMesh::SetZeroUpScatSrc(int g)
{
  for(int l=0;l<maxpl;l++){
    UpSSrc[l].put_data(g,0.);
  };
};

void GeneralMesh::AddDownScatSrc(int DepGrp,int pls,bool selfscat)
{
  int LowestScatGrp=Med->GetLowestDownScatGrp(DepGrp);
  int stt=DepGrp+1;
  if(selfscat)stt=DepGrp;

  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=stt;i<=LowestScatGrp;i++){
        SSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(DepGrp,i)*tmp);
      };
      ind++;
    };
  };
}

void GeneralMesh::AddDownScatSrcReactionWise(int DepGrp,enum xstype xst,int pls,bool selfscat)
{
  int LowestScatGrp=Med->GetLowestDownScatGrp(DepGrp);
  int stt=DepGrp+1;
  if(selfscat)stt=DepGrp;

  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=stt;i<=LowestScatGrp;i++){
        int nucn=Med->GetNucnum();
        real xs=0.;
        for(int n=0;n<nucn;n++){
          real den=Med->GetNuclideInTurn(n).GetDensity();
	  if(l==0||xst==sigel)xs+=den*Med->GetNuclideInTurn(n).GetMicxs().GetData2d(xst,l).get_dat(DepGrp,i);
	};
        SSrc[ind].add_data(i,xs*tmp);
      };
      ind++;
    };
  };
}

void GeneralMesh::AddUpScatSrc(int DepGrp,int pls)
{
  if(DepGrp==0)return;
  int hgrp=Med->GetHighestUpScatGrp();
  if(hgrp==-1)return;

  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=hgrp;i<DepGrp;i++){
        //SSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(DepGrp,i)*tmp);
        UpSSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(DepGrp,i)*tmp);
      };
      ind++;
    };
  };
}

void GeneralMesh::AddDownScatSrcAdjoint(int DepGrp,int pls)
{
  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=DepGrp-1;i>=0;i--){
        SSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(i,DepGrp)*tmp);
      };
      ind++;
    };
  };
}

void GeneralMesh::AddDownScatSrcAdjointArv(int DepGrp,int ArvMinGrp,int pls)
{
  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=DepGrp-1;i>=ArvMinGrp;i--){
        SSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(i,DepGrp)*tmp);
      };
      ind++;
    };
  };
}

void GeneralMesh::AddUpScatSrcAdjoint(int DepGrp,int pls)
{
  if(DepGrp==grp-1)return;

  real Volp=Volume*INV_PI4;
  int ind=0;
  for(int l=0;l<=pls;l++){
    int is=0;
    int ie=0;
    if(dim!=1)ie=l;
    if(dim==3)is=-l;
    for(int m=is;m<=ie;m++){
      real tmp=Flux[ind].get_dat(DepGrp)*Volp;
      for(int i=DepGrp+1;i<grp;i++){
        //SSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(i,DepGrp)*tmp);
        UpSSrc[ind].add_data(i,Med->GetMacxs().GetSigs(l).get_dat(i,DepGrp)*tmp);
      };
      ind++;
    };
  };
}

real GeneralMesh::GetScatSrc(int g,int ind)
{
  //if(ind>maxpl){
  //  cout<<"There is no moment in GeneralMesh.\n";
  //  cout<<"  The number of Pl array is "<<maxpl<<"\n";
  //  cout<<"  You requested "<<ind<<"\n";
  //  exit(0);
  //};
  return SSrc[ind].get_dat(g);
};

real GeneralMesh::GetUpScatSrc(int g,int ind)
{
  //if(ind>maxpl){
  //  cout<<"There is no moment in GeneralMesh.\n";
  //  cout<<"  The number of Pl array is "<<maxpl<<"\n";
  //  cout<<"  You requested "<<ind<<"\n";
  //  exit(0);
  //};
  return UpSSrc[ind].get_dat(g);
};

void GeneralMesh::CalFissionSrc()
{
  //FSrc=Med->GetMacxs().GetData1d("nusigf")*Flux[0]*Volume;
  FSrc=Med->GetMacxs().GetNusigf()*Flux[0]*Volume;
};

real GeneralMesh::GetFissionSrcChi(int g)
{
  return FSrc*Med->GetMacxs().GetKai().get_dat(g);
};

real GeneralMesh::GetReactionRate(enum xstype ss)
{
  return Flux[0]*Med->GetMacxs().GetData1d(ss)*Volume;
};

real GeneralMesh::GetReactionRate(int matnum,enum xstype ss)
{
  return Flux[0]*Med->GetNuclide(matnum).GetMicxs().GetData1d(ss);
};

real GeneralMesh::GetVolumeFlux()
{
  return Flux[0].get_sum()*Volume;
};

void GeneralMesh::CalFissionSrcAdjoint()
{
  FSrc=Med->GetMacxs().GetKai()*Flux[0]*Volume;
};

real GeneralMesh::GetFissionSrcSigf(int g)
{
  return FSrc*Med->GetMacxs().GetNusigf().get_dat(g);
};

GroupData1D GeneralMesh::GetFissionSrcSigf()
{
  return Med->GetMacxs().GetNusigf()*FSrc;
};

real GeneralMesh::GetDif(int g,int flag)
{
  switch(flag){
  case 0:
    //return Med->GetMacxs().GetData1d("d").get_dat(g);
    return Med->GetMacxs().GetD().get_dat(g);
  case 1:
    return Med->GetMacxs().GetDr().get_dat(g);
  case 2:
    return Med->GetMacxs().GetDz().get_dat(g);
  };
  cout<<"Error in GetDif in GeneralMesh.\n";
  cout<<"Requested flag is "<<flag<<"\n";
  exit(0);
};

void GeneralMesh::InitializeUpScatSrc(int wgrp)
{
  UpSSrc.resize(maxpl);
  for(int i=0;i<maxpl;i++){
    UpSSrc[i].put_imax(grp);
    for(int g=0;g<grp;g++){
      UpSSrc[i]=0.;
    };
  };
};

void GeneralMesh::SrcDataClear()
{
  SSrc.clear();
  UpSSrc.clear();
  srcin.clear();
};

real GeneralMesh::GetEnergyIntegratedFlux(real etop,real elow,int mom)
{
  bool in=false;
  bool out=false;

  real ret=0.;
  real e0=Med->GetEnband().get_dat(0);
  for(int i=0;i<grp;i++){
    real e1=Med->GetEnband().get_dat(i+1);
    real letwid=log(e0/e1);
    real eu=e0;
    real el=e1;
    if(e1<etop&&!in){
      in=true;
      eu=etop;
    };
    if(in&&!out){
      if(e1<elow){
        out=true;
        el=elow;
      };
      ret+=Flux[mom].get_dat(i)*log(eu/el)/letwid;
    };
  };

  return ret;
};
