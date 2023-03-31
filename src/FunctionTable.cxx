#include<iostream>
#include<cmath>
#include"FunctionTable.h"

using namespace std;

// +++ Bickley function
//
// fk(4) is already multiplied by 1.5 
// for anisotropic Pij calculation

void fktab::init(real err)
{
  real aatmp[]={
   0.96572105   , -0.28959120   , 0.030610238  , -0.0011370233 ,
   0.36668658e-2, -0.82618978e-3, 0.63930345e-4, -0.16954459e-5,
   0.43228799   , -0.12254097   , 0.012167113  , -0.42191691e-3,
   0.16388649e-2, -0.35676160e-3, 0.26622503e-4, -0.67972543e-6,
   0.084790066  , -0.022044221  , 0.0019904763 , -0.62257550e-4,
   0.32297971e-3, -0.65260999e-4, 0.45013738e-5, -0.10581271e-6};
  real bbtmp[]={
       1.0,  0.45118240  , 0.73261267   , -0.039791833  , 0.039573198  ,
       1.0, -0.78064825  , 0.24422565   , -0.035315751  , 0.0023728056 ,
       1.0, -0.33898117  , 0.67917163   , -0.097397491  , 0.029176632  ,
       1.0, -0.72078149  , 0.20596925   , -0.027412106  , 0.0016314771 ,
       1.0, -0.91375070  , 0.50917132   , -0.093411994  , 0.012892154  ,
       1.0, -0.62018356  , 0.15021226   , -0.016967420  , 0.80823019e-3};
  real cctmp[]={-0.625, 0.80368280, -1.125, 2.0285838, -1.625, 3.7538403};
  real pxtmp[]={5.9, 9.4 , 6.5, 9.8, 7.5, 11.0};

  for(int i=0;i<24;i++){aa[i]=aatmp[i];};
  for(int i=0;i<30;i++){bb[i]=bbtmp[i];};
  for(int i=0;i<6;i++){cc[i]=cctmp[i];};
  for(int i=0;i<6;i++){px[i]=pxtmp[i];};
  
  real fk_x[2202];
  real fk_y[2201*5];

  fk_x[0]=1.0;
  fk_y[0]=1.5707963267949;
  fk_y[1*2201]=1.0000000000000;
  fk_y[2*2201]=0.7853981633974;
  fk_y[3*2201]=0.6666666666667;
  //fk_y[4*2201]=0.58904862256;
  fk_y[4*2201]=0.58904862256*1.5;

  for(int i=1;i<2201;i++){
    fk_x[i]=real(i)/200.0;
    for(int n=0;n<5;n++){
      fk_y[n*2201+i]=fkin0(fk_x[i],n,err);
      if(n==4)fk_y[n*2201+i]*=1.5; // (Multiply by 1.5)
    };
  };

  int i=0;
  for(int j=0;j<1099;j++){
    i++;
    for(int n=0;n<5;n++){
      //fk_a[n*1100+j]=(4.0*fk_y[n*2201+i]-3.0*fk_y[n*2201+i-1]-fk_y[n*2201+i+1])/0.01;
      fk_a[n*1100+j]=(4.0*fk_y[n*2201+i]-3.0*fk_y[n*2201+i-1]-fk_y[n*2201+i+1])*100.;
      //fk_b[n*1100+j]=(fk_y[n*2201+i+1]-2.0*fk_y[n*2201+i]+fk_y[n*2201+i-1])/0.00005;
      fk_b[n*1100+j]=(fk_y[n*2201+i+1]-2.0*fk_y[n*2201+i]+fk_y[n*2201+i-1])*20000.;
      fk_c[n*1100+j]=fk_y[n*2201+i-1];
    };
    i++;
  };

  for(int n=0;n<5;n++){
    fk_a[n*1100+1099]=0.0;
    fk_b[n*1100+1099]=0.0;
    fk_c[n*1100+1099]=0.0;
  };

};

real fktab::xkini(real x,real y,int n)
{
  real z=cos(x);
  real ep,ret,arg,tmp;
  int ip;
  real zero=0.0;

  if(z<=0)return zero;

  if(y>0){
    ip=1;
    arg=y/z;
    for(int i=0;i<3;i++){
      if(arg>150.0){
	arg=arg/2.0;
	ip=ip*2;
      };
    };
    if(arg>150.0)return zero;

    tmp=exp(-arg);
    if(fabs(tmp)<=1e-32&&ip>=2){
      ep=0.0;
    }
    else{
      ep=pow(exp(-arg),ip);
    };
  }
  else{
    ep=1.0;
  };
  z=pow(z,n-1);
  ret=z*ep;
  return ret;
};

real fktab::simpsk(real err,real y,int n)
{
  real a=0.0;
  real b=1.5707963267949;

  real x,stw0,del,ans;

  int ndiv=1;
  real cur=1.0;
  real prev=0.0;
  real sone=(b-a)*(xkini(a,y,n)+xkini(b,y,n))/2.0;

  while(err*fabs(cur)<fabs(cur-prev)){
    ndiv*=2;
    stw0=0.0;
    del=(b-a)/real(ndiv);
    for(int i=1;i<ndiv+1;i+=2){
      x=a+del*real(i);
      stw0+=xkini(x,y,n);
    };
    prev=cur;
    cur=sone+4.0*del*stw0;
    sone=(sone+cur)/4.0;
  };
  ans=cur/3.0;
  return ans;
};

real fktab::fkin0(real y, int n, real err)
{
  real ret=simpsk(err,y,n+1);
  return ret;
};

real fktab::fkin(int n,real x)
{
  if(x>20.)return 0.;
  // This is to reduce calculation time for large value of x.
  // This function converge to 0 when the value of x increases. 
  // Modification is done in 2011/09/18.

  int nn=int(100.0*x)+1;
  real dx=x-real(nn-1)*0.01;
  int nm=n-1;

  if(nn<=1100){
    nn=nn+1100*nm-1;
    return (fk_b[nn]*dx+fk_a[nn])*dx+fk_c[nn];
  };

  real fkin;

  int j=-1;
  for(int i=0;i<2;i++){
    if(x<px[nm*2+i]&&j==-1)j=i;
  };

  if(j==1||j==0){
    real qq=aa[nm*8+j*4+3];
    int tmp=nm*8+j*4;
    for(int mm=2;mm>=0;mm--){
      qq=qq*x+aa[tmp+mm];
    };
    real pp=bb[nm*10+j*5+4];
    tmp=nm*10+j*5;
    for(int mm=3;mm>=0;mm--){
      pp=pp*x+bb[tmp+mm];
    };
    fkin=qq/pp;
  }
  else{
    fkin=1.2533141*exp(-x)/sqrt(x)*(1.0+(cc[nm*2]+cc[nm*2+1]/x)/x);
  };
  if(n==5)fkin*=1.5;
  return fkin;
};

// Exponential function

exptab::exptab()
{
  for(int i=0;i<1001;i++){
    real x=i*(-0.01);
    p0[i]=exp(x);
  };
  real dx=0.01;
  real dx_inv=1./dx;
  for(int i=0;i<1000;i++){
    p1[i]=(p0[i+1]-p0[i])*dx_inv;
  };
};

real exptab::e(real x)
{
  if(x<-10.){return 0.;};
  int i=int(x*(-100));
  float x2=i*-0.01;
  float dx=x2-x;
  return p0[i]+p1[i]*dx;
};

// FUNCMOC::E(x)=(1-exp(-x))/x

funcmoc::funcmoc()
{
  // The range of 0 to 100 is divided into [num] meshes.

  num=100000;

  num2=num/100;
  wid=100./num;

  p0.resize(num+1);
  p1.resize(num);

  p0[0]=1.;
  for(int i=1;i<num+1;i++){
    //real x=i*0.01;
    real x=i*wid;
    p0[i]=(1.-exp(-x))/x;
  };
  for(int i=0;i<num;i++){
    p1[i]=(p0[i+1]-p0[i])*num2;
  };
};

real funcmoc::get(real x)
{
  if(x>100.){
    return (1.-exp(-x))/x;
  }else{;
    int i=int(x*num2);
    float x2=i*wid;
    float dx=x2-x;
    return p0[i]+p1[i]*dx;
  };
};
