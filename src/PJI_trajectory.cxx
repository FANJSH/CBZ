#include <cstdlib>
#include "PJI_trajectory.h"

using namespace std;

void Trajectory::Initialize(int inp)
{
  numreg = inp;
  reg.resize(numreg,0);
  opt.resize(numreg,0.);
  dist.resize(numreg,0.);
  xsinv.resize(numreg,0.);
};

void Trajectory::PutData(int *ireg, real *idist)
{
  vector<int> reg2(numreg);
  vector<real> dist2(numreg);
  for(int i=0;i<numreg;i++){
    reg2[i]=ireg[i];
    dist2[i]=idist[i];
  };
  PutData(reg2,dist2);
};

void Trajectory::PutData(vector<int> ireg, vector<real> idist)
{
  reg[0] =ireg[0];
  dist[0]=idist[0];
  for(int i=1;i<numreg;i++){
    reg[i] =ireg[i];
    real tmp=idist[i]-idist[i-1];
    if(tmp<=0.){
      cout<<"Error in instance of Trajectory class.\n";    
      cout<<"Error in order of distance.\n";
      for(int i=0;i<numreg;i++){
	cout<<i<<" "<<idist[i]<<"\n";
      };
      exit(0);
    };
    dist[i]=tmp;
  };
};

void Trajectory::PutXS(real *xs)
{
  TotalOpt=0.;
  for(int i=0;i<numreg;i++){
    real x=xs[reg[i]];
    real tmp=dist[i]*x;
    opt[i]=tmp;
    xsinv[i]=1./x;
    TotalOpt+=tmp;
  };
};

void Trajectory::GetRegionwiseLength(int rgt,real *length)
{
  for(int i=0;i<rgt;i++){length[i]=0;};
  for(int i=0;i<numreg;i++){
    length[reg[i]]+=dist[i];
  };
};

void Trajectory::CalPij(int rgt,fktab &fin,real *retpij,bool aniso)
{
  int val=3;
  if(aniso)val=5;

  int rgt1=rgt+1;
  for(int i=0;i<rgt1*rgt1;i++){retpij[i]=0;};
  real tau,taum,taun,taumold,taunold;
  real opti;
  int regi,regj;
  real fkin0=fin.fkin(val,0.0);

  // i=0
  regi=reg[0];
  opti=opt[0];
  tau =0.0;
  taum=fkin0;
  taun=fin.fkin(val,opti);

  int regi_rgt=regi*rgt;
  int rgt_rgt1=rgt*rgt1;
  int rgt_rgt=rgt*rgt;
  retpij[regi_rgt+regi]+=opti-taum+taun; // (i->i)
  retpij[rgt_rgt1+regi]+=taum-taun; // (s->i)
  for(int j=1;j<numreg;j++){
    regj=reg[j];
    taumold=taum;
    taunold=taun;
    tau+=opt[j];
    taum=fin.fkin(val,tau);
    taun=fin.fkin(val,tau+opti);
    retpij[regi_rgt+regj]+=taumold-taunold-taum+taun; // (i->j)
    retpij[rgt_rgt1+regj]+=taunold-taun; // (s->j)
  };
  retpij[rgt_rgt+regi]+=taum-taun;   // (i->s)
  retpij[rgt1*rgt1-1]+=taun; // (s->s)

  for(int i=1;i<numreg;i++){
    regi=reg[i];
    opti=opt[i];
    tau =0.0;
    taum=fkin0-fin.fkin(val,opti);
    retpij[regi*rgt+regi]+=opti-taum; // (i->i)
    int tmp=regi*rgt;
    for(int j=i+1;j<numreg;j++){
      taumold=taum;
      tau+=opt[j];
      taum=fin.fkin(val,tau)-fin.fkin(val,tau+opti);
      retpij[tmp+reg[j]]+=taumold-taum; // (i->j)
    };
    retpij[rgt_rgt+regi]+=taum;   // (i->s)
  };
};

void Trajectory::CalPijSphereWithExpTable(int rgt,exptab &etab,real *retpij)
{
  int rgt1=rgt+1;
  for(int i=0;i<rgt1*rgt1;i++){retpij[i]=0;};

  for(int i=0;i<numreg;i++){
    int reg1=reg[i];
    real opt1=opt[i];
    real val1=1.-etab.e(-opt1);
    retpij[reg1*rgt1+reg1]+=(opt1-val1); // i -> i
    if(i==0)retpij[rgt*rgt1+reg1]+=val1; // s -> i
    real optc=0.;
    for(int j=i+1;j<numreg;j++){
      int reg2=reg[j];
      real opt2=opt[j];
      real val2=1.-etab.e(-opt2); 
      retpij[reg1*rgt1+reg2]+=(val1*etab.e(-optc)*val2);    // i -> j
      if(i==0)retpij[rgt*rgt1+reg2]+=etab.e(-optc-opt1)*val2; // s -> j
      optc+=opt2;
    };
    retpij[reg1*rgt1+rgt]+=val1*etab.e(-optc); // i -> s
    if(i==0)retpij[rgt*rgt1+rgt]+=etab.e(-optc-opt1); // s -> s
  };
};

void Trajectory::CalPijSphere(int rgt,real *retpij)
{
  int rgt1=rgt+1;
  for(int i=0;i<rgt1*rgt1;i++){retpij[i]=0;};

  for(int i=0;i<numreg;i++){
    int reg1=reg[i];
    real opt1=opt[i];
    real val1=1.-exp(-opt1);
    retpij[reg1*rgt1+reg1]+=(opt1-val1); // i -> i
    if(i==0)retpij[rgt*rgt1+reg1]+=val1; // s -> i
    real optc=0.;
    for(int j=i+1;j<numreg;j++){
      int reg2=reg[j];
      real opt2=opt[j];
      real val2=1.-exp(-opt2); 
      retpij[reg1*rgt1+reg2]+=(val1*exp(-optc)*val2);    // i -> j
      if(i==0)retpij[rgt*rgt1+reg2]+=exp(-optc-opt1)*val2; // s -> j
      optc+=opt2;
    };
    retpij[reg1*rgt1+rgt]+=val1*exp(-optc); // i -> s
    if(i==0)retpij[rgt*rgt1+rgt]+=exp(-optc-opt1); // s -> s
  };
};

void Trajectory::CalPijAniso(int rgt,fktab &fin,real *retpij)
{
  CalPij(rgt,fin,retpij,true);
};

void Trajectory::CalPij
(int rgt,Trajectory &sec,real Optb,real *retpij,fktab &fin,bool aniso,bool oppo1,bool oppo2)
{
  int val=3;
  if(aniso)val=5;

  real tau2,taumn,taumnold;
  real tau=TotalOpt+Optb;
  int regi;

  for(int i=0;i<numreg;i++){
    int ii=i;
    if(oppo1)ii=numreg-1-i;
    regi=reg[ii];
    real opti=opt[ii];
    tau-=opti;
    tau2=tau;
    taumn=fin.fkin(val,tau2+opti)-fin.fkin(val,tau2);
    int tmp=regi*rgt;
    for(int j=0;j<sec.numreg;j++){
      int jj=j;
      if(oppo2)jj=sec.numreg-1-j;
      taumnold=taumn;
      tau2+=sec.opt[jj];
      taumn=fin.fkin(val,tau2+opti)-fin.fkin(val,tau2);
      retpij[tmp+sec.reg[jj]]+=taumn-taumnold;
    };
  };
};

void Trajectory::CalPijAniso(int rgt,Trajectory &sec,real Optb,real *retpij,fktab &fin,bool oppo1,bool oppo2)
{
  CalPij(rgt,sec,Optb,retpij,fin,true,oppo1,oppo2);
};

void Trajectory::CalPijIso(int rgt,Trajectory &sec,real Optb,real *retpij,fktab &fin,bool oppo1,bool oppo2)
{
  CalPij(rgt,sec,Optb,retpij,fin,false,oppo1,oppo2);
};

void Trajectory::show_self()
{
  if(numreg>0){
    for(int i=0;i<numreg;i++){
      cout<<dist[i]<<":"<<reg[i]<<"//";
    };
    cout<<"\n";
  };
};

void Trajectory::AdjustTrajectoryLength(Trajectory &sec,Trajectory &thi)
  // ***
  //  This method is to increase accuracy in macroband calculation.
  //  With neighboring two trajectories, length is fitted by 
  //  second order polynomials.
  // ***
{
  int numreg2=sec.GetNumreg();
  int numreg3=thi.GetNumreg();
  if(numreg2<1||numreg3<1)return;

  vector<int> reg2(numreg);
  vector<int> reg3(numreg);
  vector<real> dist2(numreg);
  vector<real> dist3(numreg);

  reg2[0]=999;
  reg3[0]=999;

  if(numreg==numreg2){
    for(int i=0;i<numreg;i++){
      reg2[i]=sec.GetReg(i);
      dist2[i]=sec.GetDist(i);
    };
  }else if(numreg2==numreg+1){
    int id=0;
    int is=0;
    for(int i=0;i<numreg;i++){
      if(reg[i]==sec.GetReg(id)){
        reg2[i]=reg[i];
	dist2[i]=sec.GetDist(id);
	id++;
      }
      else{
	is++;
	id++;
	reg2[i]=sec.GetReg(id);
	dist2[i]=sec.GetDist(id);
	id++;
      };
    };
  }else if(numreg2==numreg-1){
    int id=0;
    int is=0;
    for(int i=0;i<numreg;i++){
      if(reg[i]==sec.GetReg(id)||is==1){
        reg2[i]=reg[i];
	dist2[i]=sec.GetDist(id);
	id++;
      }
      else{
	is++;
	reg2[i]=reg[i];
	dist2[i]=dist[i];
      };
    };
  };

  if(numreg==numreg3){
    for(int i=0;i<numreg;i++){
      reg3[i]=thi.GetReg(i);
      dist3[i]=thi.GetDist(i);
    };
  }else if(numreg3==numreg+1){
    int id=0;
    int is=0;
    for(int i=0;i<numreg;i++){
      if(reg[i]==thi.GetReg(id)||is==1){
        reg3[i]=reg[i];
	dist3[i]=thi.GetDist(id);
	id++;
      }
      else{
	is++;
	id++;
	reg3[i]=thi.GetReg(id);
	dist3[i]=thi.GetDist(id);
	id++;
      };
    };
  }else if(numreg3==numreg-1){
    int id=0;
    int is=0;
    for(int i=0;i<numreg;i++){
      if(reg[i]==thi.GetReg(id)||is==1){
        reg3[i]=reg[i];
	dist3[i]=thi.GetDist(id);
	id++;
      }
      else{
	is++;
	reg3[i]=reg[i];
	dist3[i]=dist[i];
      };
    };
  };

  bool agree2=true;
  for(int i=0;i<numreg;i++){
    if(reg[i]!=reg2[i]){
      agree2=false;
      break;
    };
  };

  bool agree3=true;
  for(int i=0;i<numreg;i++){
    if(reg[i]!=reg3[i]){
      agree3=false;
      break;
    };
  };

  if(agree2&&agree3){
    for(int i=0;i<numreg;i++){
      real d1=dist2[i];
      real d3=dist3[i];
      real d2=dist[i];
      real a=(d3-2*d2+d1)*0.125;
      real b=(-d3+3*d2-2*d1)*0.5;
      real c=d1-a-b;
      dist[i]=12*a+3*b+c;
    };
  }else if(agree2&&!agree3){
    for(int i=0;i<numreg;i++){
      dist[i]=(2*dist[i]+dist2[i])*0.33333333;
    };
  }else if(!agree2&&agree3){
    for(int i=0;i<numreg;i++){
      dist[i]=(2*dist[i]+dist3[i])*0.33333333;
    };
  };
};

// +++ for CBG/MEC

real Trajectory::SolveCFormTransport(exptab &e,real flxin,real *src,real *aveflx)
{
  real flxout=0.;
  for(int i=0;i<numreg;i++){
    real optlen=opt[i];
    //real dem=exp(-optlen);
    real dem=e.e(-optlen);
    real optlen_inv=1./optlen;
    flxout=flxin*dem+src[i]*optlen_inv*dist[i]*(1.-dem);
    aveflx[i]=(flxin-flxout+src[i]*dist[i])*optlen_inv;
    flxin=flxout;
  };
  return flxout;
};

real Trajectory::SolveCFormTransport(exptab &e,real flxin,real *src,real *aveflx,real sin_inv)
{
  real flxout=0.;
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    //real dem=exp(-optlen);
    real dem=e.e(-optlen);
    real optlen_inv=1./optlen;
    real distlen=dist[i]*sin_inv;
    real srctmp=src[i];
    flxout=flxin*dem+srctmp*optlen_inv*distlen*(1.-dem);
    aveflx[i]=(flxin-flxout+srctmp*distlen)*optlen_inv;
    flxin=flxout;
  };
  return flxout;
};

void Trajectory::SolveCFormTransport(exptab &e,real flxin,real *src,real *aveflx,real *outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real dem=e.e(-optlen);
    //real dem=exp(-optlen);
    real optlen_inv=1./optlen;
    real distlen=dist[i]*sin_inv;
    real srctmp=src[i];
    real flxout=flxin*dem+srctmp*optlen_inv*distlen*(1.-dem);
    aveflx[i]=(flxin-flxout+srctmp*distlen)*optlen_inv;
    if(aveflx[i]>500.){
      cout<<flxin<<" "<<dem<<" "<<optlen<<" "<<srctmp<<" "<<distlen<<" "<<flxout<<" "<<aveflx[i]<<"\n";
    };
    flxin=flxout;
    outflx[i]=flxout; //
  };
  return;
};

void Trajectory::SolveCFormTransportOpposite(exptab &e,real flxin,real *src,real *aveflx,real *outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real optlen=opt[ii]*sin_inv;
    real dem=e.e(-optlen);
    //real dem=exp(-optlen);
    real optlen_inv=1./optlen;
    real distlen=dist[ii]*sin_inv;
    real srctmp=src[ii];
    real flxout=flxin*dem+srctmp*optlen_inv*distlen*(1.-dem);
    aveflx[ii]=(flxin-flxout+srctmp*distlen)*optlen_inv;
    flxin=flxout;
    outflx[ii]=flxout; //
  };
  return;
};

void Trajectory::SolveCFormTransport(funcmoc &fmoc,real flxin,real *src,real *aveflx,real *outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real tmp=distlen*src[i];
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      //real dem=(1.-exp(-optlen))/optlen;
      real flxout=flxin+dem*(tmp-optlen*flxin);
      aveflx[i]=(flxin-flxout+tmp)/optlen*dist[i];
      outflx[i]=flxout;
      flxin=flxout;
    }else{
      real flxout=flxin+tmp;
      outflx[i]=flxout;
      aveflx[i]=(flxin+outflx[i])*0.5*dist[i];
      flxin=flxout;
    };
  };

  return;
};

void Trajectory::SolveCFormTransportOpposite(funcmoc &fmoc,real flxin,real *src,real *aveflx,real *outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real optlen=opt[ii]*sin_inv;
    real distlen=dist[ii]*sin_inv;
    real tmp=distlen*src[ii];
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      //real dem=(1.-exp(-optlen))/optlen;
      real flxout=flxin+dem*(tmp-optlen*flxin);
      aveflx[ii]=(flxin-flxout+tmp)/optlen*dist[ii];
      outflx[ii]=flxout;
      flxin=flxout;
    }else{
      outflx[ii]=flxin+tmp;
      aveflx[ii]=(flxin+outflx[ii])*0.5*dist[ii];
    };
  };

  return;
};

void Trajectory::SolveCFormTransport(funcmoc &fmoc,real flxin,vector<real> &src,vector<real> &aveflx,real &outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real tmp=distlen*src[i];
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      real flxout=flxin+dem*(tmp-optlen*flxin);
      aveflx[i]=(flxin-flxout+tmp)/optlen*dist[i];
      flxin=flxout;
    }else{
      real flxout=flxin+tmp;
      aveflx[i]=(flxin+flxout)*0.5*dist[i];
      flxin=flxout;
    };
  };
  outflx=flxin;

  return;
};

void Trajectory::SolveCFormTransport(funcmoc &fmoc,real flxin,real *src,real *aveflx,real &outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real tmp=distlen*src[i];
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      real flxout=flxin+dem*(tmp-optlen*flxin);
      aveflx[i]=(flxin-flxout+tmp)/optlen*dist[i];
      flxin=flxout;
    }else{
      real flxout=flxin+tmp;
      aveflx[i]=(flxin+flxout)*0.5*dist[i];
      flxin=flxout;
    };
  };
  outflx=flxin;

  /*
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real srctmp=src[i];
    if(optlen!=0.){
      //real dem=fmoc.get(optlen);
      real dem=(1.-exp(-optlen))/optlen;
      real optlen_inv=1./optlen;
      real flxout=flxin+distlen*dem*srctmp-optlen*dem*flxin;
      //aveflx[i]=(flxin-flxout+srctmp*distlen)*optlen_inv;
      aveflx[i]=(flxin-flxout+srctmp*distlen)*optlen_inv*dist[i];
      outflx=flxout; //
    }else{
      outflx=flxin+srctmp*distlen;
      aveflx[i]=(flxin+outflx)*0.5*dist[i];
    };
    flxin=outflx;
  };
  */

  return;
};

void Trajectory::SolveCFormTransportOpposite(funcmoc &fmoc,real flxin,real *src,real *aveflx,real &outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real optlen=opt[ii]*sin_inv;
    real distlen=dist[ii]*sin_inv;
    real tmp=distlen*src[ii];
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      real flxout=flxin+dem*(tmp-optlen*flxin);
      aveflx[ii]=(flxin-flxout+tmp)/optlen*dist[ii];
      flxin=flxout;
    }else{
      real flxout=flxin+tmp;
      aveflx[ii]=(flxin+flxout)*0.5*dist[ii];
      flxin=flxout;
    };
  };
  outflx=flxin;

  /*
  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real optlen=opt[ii]*sin_inv;
    real srctmp=src[ii];
    real distlen=dist[ii]*sin_inv;
    if(optlen!=0.){
      real dem=fmoc.get(optlen);
      //real dem=(1.-exp(-optlen))/optlen;
      real optlen_inv=1./optlen;
      real flxout=flxin+distlen*dem*srctmp-optlen*dem*flxin;
      //aveflx[ii]=(flxin-flxout+srctmp*distlen)*optlen_inv;
      aveflx[ii]=(flxin-flxout+srctmp*distlen)*optlen_inv*dist[ii];
      outflx=flxout;
    }else{
      outflx=flxin+srctmp*distlen;
      aveflx[ii]=(flxin+outflx)*0.5*dist[ii];
    };
    flxin=outflx;
  };
  */
  return;
};

void Trajectory::SolveCFormTransportDoubleDirection(funcmoc &fmoc,real flxin1,real flxin2,real *src,real *aveflx1,real *aveflx2,real &outflx1,real &outflx2,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real tmp1=opt[i]*sin_inv;
    real tmp3=dist[i]*sin_inv*src[i];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout1=flxin1+dem*(tmp3-tmp1*flxin1);
      aveflx1[i]=(flxin1-flxout1+tmp3)/tmp1*dist[i];
      flxin1=flxout1;
    }else{
      real flxout1=flxin1+tmp3;
      aveflx1[i]=(flxin1+flxout1)*0.5*dist[i];
      flxin1=flxout1;
    };
  };
  outflx1=flxin1;

  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real tmp1=opt[ii]*sin_inv;
    real tmp3=dist[ii]*sin_inv*src[ii];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout2=flxin2+dem*(tmp3-tmp1*flxin2);
      aveflx2[ii]=(flxin2-flxout2+tmp3)/tmp1*dist[ii];
      flxin2=flxout2;
    }else{
      real flxout2=flxin2+tmp3;
      aveflx2[ii]=(flxin2+flxout2)*0.5*dist[ii];
      flxin2=flxout2;
    };
  };
  outflx2=flxin2;
};

void Trajectory::SolveCFormTransportDoubleDirection(funcmoc &fmoc,real flxin1,real flxin2,real *src,real *aveflx1,real *aveflx2,vector<real> &outflx1,vector<real> &outflx2,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real tmp1=opt[i]*sin_inv;
    real tmp3=dist[i]*sin_inv*src[i];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout1=flxin1+dem*(tmp3-tmp1*flxin1);
      aveflx1[i]=(flxin1-flxout1+tmp3)/tmp1*dist[i];
      flxin1=flxout1;
    }else{
      real flxout1=flxin1+tmp3;
      aveflx1[i]=(flxin1+flxout1)*0.5*dist[i];
      flxin1=flxout1;
    };
    outflx1[i]=flxin1;
  };
  //outflx1=flxin1;

  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real tmp1=opt[ii]*sin_inv;
    real tmp3=dist[ii]*sin_inv*src[ii];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout2=flxin2+dem*(tmp3-tmp1*flxin2);
      aveflx2[ii]=(flxin2-flxout2+tmp3)/tmp1*dist[ii];
      flxin2=flxout2;
    }else{
      real flxout2=flxin2+tmp3;
      aveflx2[ii]=(flxin2+flxout2)*0.5*dist[ii];
      flxin2=flxout2;
    };
    outflx2[i]=flxin2;
  };
  //outflx2=flxin2;
};


void Trajectory::SolveCFormTransportDoubleDirection(funcmoc &fmoc,real flxin1,real flxin2,real *src,real *aveflx1,real *aveflx2,real *outflx1,real *outflx2,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real tmp1=opt[i]*sin_inv;
    real tmp3=dist[i]*sin_inv*src[i];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout1=flxin1+dem*(tmp3-tmp1*flxin1);
      aveflx1[i]=(flxin1-flxout1+tmp3)/tmp1*dist[i];
      flxin1=flxout1;
    }else{
      real flxout1=flxin1+tmp3;
      aveflx1[i]=(flxin1+flxout1)*0.5*dist[i];
      flxin1=flxout1;
    };
    outflx1[i]=flxin1;
  };
  //outflx1=flxin1;

  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real tmp1=opt[ii]*sin_inv;
    real tmp3=dist[ii]*sin_inv*src[ii];
    if(tmp1!=0.){
      real dem=fmoc.get(tmp1);
      real flxout2=flxin2+dem*(tmp3-tmp1*flxin2);
      aveflx2[ii]=(flxin2-flxout2+tmp3)/tmp1*dist[ii];
      flxin2=flxout2;
    }else{
      real flxout2=flxin2+tmp3;
      aveflx2[ii]=(flxin2+flxout2)*0.5*dist[ii];
      flxin2=flxout2;
    };
    outflx2[i]=flxin2;
  };
  //outflx2=flxin2;
};


void Trajectory::SolveCFormTransport(funcmoc &fmoc,real flxin,real *src,real *aveflx,real *outflx,real *ep,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real srctmp=src[i];
    if(optlen!=0.){
    real flxout=flxin+ep[i]*(srctmp-optlen/distlen*flxin);
    aveflx[i]=(flxin-flxout+srctmp*distlen)/optlen;
    flxin=flxout;
    outflx[i]=flxout; //
      }else{
      outflx[i]=flxin+srctmp*distlen;
      aveflx[i]=(flxin+outflx[i])*0.5;
    };
  };
  return;
};

void Trajectory::SolveCFormTransportOpposite(funcmoc &fmoc,real flxin,real *src,real *aveflx,real *outflx,real *ep,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    int ii=numreg-1-i;
    real optlen=opt[ii]*sin_inv;
    real distlen=dist[ii]*sin_inv;
    real srctmp=src[ii];
    if(optlen!=0.){
      real flxout=flxin+ep[ii]*(srctmp-optlen/distlen*flxin);
      aveflx[ii]=(flxin-flxout+srctmp*distlen)/optlen;
      flxin=flxout;
      outflx[ii]=flxout; //
    }else{
      outflx[ii]=flxin+srctmp*distlen;
      aveflx[ii]=(flxin+outflx[ii])*0.5;
    };
  };
  return;
};

void Trajectory::SolveCFormTransportDD(real flxin,real *src,real *aveflx,real *outflx,real sin_inv)
{
  for(int i=0;i<numreg;i++){
    real optlen=opt[i]*sin_inv;
    real distlen=dist[i]*sin_inv;
    real srctmp=src[i];
    real coef=1./(2.+optlen);
    real flxout=(flxin*(2.-optlen)+srctmp*2.*distlen)*coef;
    if(flxout<0.){
      flxout=0.;
      aveflx[i]=(flxin+srctmp*distlen)/optlen;
    }else{
      aveflx[i]=(flxin*2.+srctmp*distlen)*coef;
    };
    flxin=flxout;
    outflx[i]=flxout;
  };
};

void Trajectory::LengthFactorize(real factor)
{
  for(int i=0;i<numreg;i++){
    dist[i]*=factor;
    opt[i]*=factor;
  };
  TotalOpt*=factor;
};
