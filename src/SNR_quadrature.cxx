#include <cstdlib>
#include <iostream>
#include "SNR_quadrature.h"

using namespace std;

SNRQuadrature::SNRQuadrature(int p)
{
  exist=false;
  pl=p;
};

SNRQuadrature::~SNRQuadrature()
{
}

void SNRQuadrature::PutSN(int sninp)
{
  sn=sninp;
  omega.resize(sn,0.);
  mu.resize(sn,0.);
  mc.resize(sn,0.);
  val.resize(sn+1);
  xref.resize(sn,0);
  for(int i=0;i<sn+1;i++){
    val[i].resize(pl+1,0.);
  };
  exist=true;
};

real SNRQuadrature::GetMoment(int m, int is)
{
  return val[is][m];
};

real SNRQuadrature::GetMoment(int m, real i)
{
  return Legendre(m,i);
};

void SNRQuadrature::CalMC()
{
  real sum=0.;
  for(int i=0;i<sn;i++){
    mc[i]=sum;
    sum+=omega[i]*mu[i];
  };
};

void SNRQuadrature::CalWeightGaussian()
{
  for(int i=0;i<sn;i++){
    real md=mu[i];
    omega[i]=2.*(1-pow(mu[i],2))/(sn*sn*pow(GetMoment(sn-1,md),2));
  };
  WeightNormalize();
};

void SNRQuadrature::WeightNormalize()
{
  real tot=0.;
  for(int i=0;i<sn;i++){
    tot+=omega[i];
  };
  for(int i=0;i<sn;i++){
    omega[i]*=2./tot;
  };
};

void SNRQuadrature::SymmetricMU()
{
  for(int i=sn/2;i<sn;i++){
    mu[i]=-mu[sn-i-1];
  };
};

void SNRQuadrature::PutGaussian(int s,bool odd_permit)
{
  if(s%2!=0&&!odd_permit){
    cout<<"Error in SNRQuadrature::PutGaussian.\n";
    cout<<"Gaussian quadrature order should be even.\n";
    //exit(0);
  };

  PutSN(s);

  real *tmp=new real[s];
  RootOfLegendre(s,tmp);
  for(int i=0;i<s;i++){
    mu[i]=tmp[i];
  };
  CalWeightGaussian();
  CalMC();
  CalValue();
  CalXYZref();
};

void SNRQuadrature::PutDoubleGaussian(int s)
{
  if(s%2!=0){
    cout<<"Error in SNRQuadrature::PutDoubleGaussian.\n";
    cout<<"Quadrature order should be even.\n";
    exit(0);
  };

  PutSN(s);

  SNRQuadrature gs(0);
  gs.PutGaussian(s/2,true);

  for(int i=0;i<s/2;i++){
    mu[i]=(gs.GetMu(i)-1.)*0.5;
    omega[i]=gs.GetOmega(i)*0.5;
    omega[s-1-i]=omega[i];
  };
  SymmetricMU();
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYZref();
};

void SNRQuadrature::CalValue()
{
  for(int i=0;i<sn;i++){
    real mui=mu[i];
    for(int l=0;l<=pl;l++){
      val[i][l]=Legendre(l,mui);
    };
  };
  for(int l=0;l<=pl;l++){
    val[sn][l]=Legendre(l,-1.);
  };
};

void SNRQuadrature::show_self()
{
  if(exist){
    cout<<"sn : "<<sn<<"\n";
    for(int i=0;i<sn;i++){
      cout<<mu[i]<<"  "<<omega[i]<<"\n";
    };
  };
};

void SNRQuadrature::CalXYZref()
{
  for(int i=0;i<sn;i++){
    real m=mu[i];
    //real e=eata[i];
    //real x=xi[i];
    for(int j=0;j<sn;j++){
      real m2=mu[j];
      //real e2=eata[j];
      //real x2=xi[j];
      if(fabs(m+m2)<0.0001)xref[i]=j;
    };
  };
};
