#include <cstdlib>
#include <iostream>
#include "SNRZ_quadrature.h"

using namespace std;

SNRZQuadrature::SNRZQuadrature(int p)
{
  exist=false;
  pl=p;
  int i=0;
  for(int j=0;j<=p;j++){
    i+=(j+1); // for 2D
  };
  plnum=i;
};

SNRZQuadrature::~SNRZQuadrature()
{
}

void SNRZQuadrature::PutSnnum(int sninp)
{
  snnum=sninp;
  omega.resize(snnum,0.);
  mu.resize(snnum,0.);
  xi.resize(snnum,0.);
  xref.resize(snnum,0);
  yref.resize(snnum,0);
  xyref.resize(snnum,0);
  mc.resize(snnum,0.);
  val.resize(snnum);
  for(int i=0;i<snnum;i++){
    val[i].resize(plnum,0.);
  };
  exist=true;
};

void SNRZQuadrature::PutLevelSymmetric(int sn)
{
  if(sn==2){PutS2();}
  else if(sn==4){PutS4();}
  else if(sn==6){PutS6();}
  else if(sn==8){PutS8();}
  else if(sn==12){PutS12();}
  else if(sn==16){PutS16();}
  else if(sn==20){PutS20();}
  else{
    cout<<"Error in SNRZQuadrature.\n";
    cout<<"  You requested "<<sn<<"-order.\n";
    cout<<"  You cannot.\n";
    exit(0);
  };
};

void SNRZQuadrature::CheckOrthogonalityTo00(real eps)
{
  int index=0;
  for(int l=0;l<=pl;l++){
    for(int m=0;m<=l;m++){
      real sum=0.;
      for(int i=0;i<snnum;i++){
        sum+=GetOmega(i)*GetMoment(index,i);
      };
      index++;
      if(l!=0&&fabs(sum)>eps){
	cout<<"!! Caution !!\n";
	cout<<" This quadrature set does NOT satisfy the orthogonal relation to (0,0).\n";
	cout<<"   l="<<l<<", m="<<m<<"   ("<<sum<<")\n";
      };
    };
  };
};

void SNRZQuadrature::PutS2()
{
  real p[]={0.5773503};
  real w[]={0.33333333};
  int wm[]={0,0,0,0};
  sn=2;

  sizex.resize(sn);
  int index=0;
  for(int i=0;i<sn/2;i++){
    sizex[index]=(i+1)*2;
    index++;
  };
  for(int i=0;i<sn/2;i++){
    sizex[index]=(sn-i*2);
    index++;
  };

  PutSnnum(4);
  int id=0;
  for(int i=0;i<sn/2;i++){
    for(int j=0;j<i+1;j++){
      mu[id]=-1.*p[i-j];
      xi[id]=-1.*p[sn/2-1-i];
      id++;
    };
    for(int j=0;j<i+1;j++){
      mu[id]=p[j];
      xi[id]=-1.*p[sn/2-1-i];
      id++;
    };
  };
  for(int i=0;i<sn/2;i++){
    for(int j=0;j<sn/2-i;j++){
      mu[id]=-1.*p[sn/2-i-1-j];
      xi[id]=p[i];
      id++;
    };
    for(int j=0;j<sn/2-i;j++){
      mu[id]=p[j];
      xi[id]=p[i];
      id++;
    };
  };

  for(int i=0;i<snnum;i++){
    omega[i]=w[wm[i]];
  };

  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS4()
{
  real p[]={0.30163878, 0.90444905};
  real w[]={0.33333333};
  int wm[]={0,0,0, 0,0,0, 0,0,0, 0,0,0};
  sn=4;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS6()
{
  real p[]={0.23009194, 0.6881343, 0.9455768};
  real w[]={0.40971693, 0.42361639};
  int wm[]={0,0, 1,1,1,1, 0,1,0,0,1,0,
            0,1,0,0,1,0, 1,1,1,1, 0,0};
  sn=6;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS8()
{
  // (Kobayashi's text)
  real p[]={0.2182179, 0.5773503, 0.7867958, 0.9511897};
  real w[]={0.1209877, 0.0907407, 0.0925926};
  // (TWODANT)
  //real p[]={0.19232747, 0.57735026, 0.79352176, 0.96229947};
  //real w[]={0.029197117, 0.023313807, 0.0225258};
  int wm[]={0,0, 1,1,1,1, 1,2,1,1,2,1, 0,1,1,0,0,1,1,0,
            0,1,1,0,0,1,1,0, 1,2,1,1,2,1, 1,1,1,1, 0,0};
  sn=8;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS12()
{
  // From TWODANT
  //real p[]={0.15395746, 0.45769112, 0.6286966, 0.76225828,
  //           0.87568027, 0.97600932};
  //real w[]={0.07332178, 0.0526674, 0.04161495, 0.03895667,
  //           0.03249018 };
  // From ENDOSN
  real p[]={0.1672126528, 0.4595476346, 0.6280190966, 0.760021014834,
            0.8722705439, 0.97163771925};
  real w[]={0.008845323746, 0.00698513769, 0.0046672092199, 0.0062852376325,
             0.0032314114569688 };
  int wm[]={0,0, 1,1,1,1, 2,3,2,2,3,2, 2,4,4,2,2,4,4,2, 
            1,3,4,3,1,1,3,4,3,1, 0,1,2,2,1,0,0,1,2,2,1,0,
            0,1,2,2,1,0,0,1,2,2,1,0, 1,3,4,3,1,1,3,4,3,1,
            2,4,4,2,2,4,4,2, 2,3,2,2,3,2, 1,1,1,1, 0,0};
  sn=12;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS16()
{
  // From Dr. Endo
  real p[]={0.138956875, 0.3922892614, 0.537096561301, 0.65042645062,
            0.746750573615, 0.83199655691, 0.90928550094, 0.98050087901};
  real w[]={0.00612340489, 0.0051661997, 0.0028059497113, 0.0042023308617,
             0.00305070985, 0.0019592377, 0.004615716388, 0.00076102049208};
  // TWODANT
  //real p[]={0.13344572, 0.39119433, 0.53689687, 0.65075610,
  //           0.74746822, 0.83302700, 0.91058181, 0.98203079};
  //real w[]={0.05415425, 0.03679653, 0.02777273, 0.02494275,
  //           0.02580284, 0.01962325, 0.01879762, 0.01544801};
  int wm[]={0,0,  1,1,1,1, 2,3,2,2,3,2, 4,5,5,4,4,5,5,4,
            4,6,7,6,4,4,6,7,6,4,
            2,5,7,7,5,2,2,5,7,7,5,2,
            1,3,5,6,5,3,1,1,3,5,6,5,3,1,
            0,1,2,4,4,2,1,0,0,1,2,4,4,2,1,0,
            0,1,2,4,4,2,1,0,0,1,2,4,4,2,1,0,
            1,3,5,6,5,3,1,1,3,5,6,5,3,1,
            2,5,7,7,5,2,2,5,7,7,5,2,
            4,6,7,6,4,4,6,7,6,4,
            4,5,5,4,4,5,5,4, 2,3,2,2,3,2, 1,1,1,1, 0,0};
  sn=16;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

void SNRZQuadrature::PutS20()
{
  // From Dr. Endo
  real p[]={0.120603343039, 0.34757429231644, 0.476519266144, 0.5773502691896,
            0.66302040365313, 0.7388225619101, 0.8075404016608, 0.87085258376,
            0.92986393895477, 0.9853474855580};
  real w[]={0.004627631133, 0.0041605270671, 0.00174587686158, 0.00363564154006,
            0.000778991255756, 0.0032770837505, 0.000285942423601, 0.004549891127614,
            0.0011238245018, 0.0011238245018, 0.000372008640195, 0.0013696348440};
  int wm[]={0,0, 1,1,1,1, 2,5,2,2,5,2, 3,6,6,3,3,6,6,3, 
            4,7,9,7,4,4,7,9,7,4,
            4,8,10,10,8,4,4,8,10,10,8,4,
            3,7,10,11,10,7,3,3,7,10,11,10,7,3,
	    2,6,9,10,10,9,6,2,2,6,9,10,10,9,6,2,
            1,5,6,7,8,7,6,5,1,1,5,6,7,8,7,6,5,1,
            0,1,2,3,4,4,3,2,1,0,0,1,2,3,4,4,3,2,1,0,
            0,1,2,3,4,4,3,2,1,0,0,1,2,3,4,4,3,2,1,0,
            1,5,6,7,8,7,6,5,1,1,5,6,7,8,7,6,5,1,
	    2,6,9,10,10,9,6,2,2,6,9,10,10,9,6,2,
            3,7,10,11,10,7,3,3,7,10,11,10,7,3,
            4,8,10,10,8,4,4,8,10,10,8,4,
            4,7,9,7,4,4,7,9,7,4,
            3,6,6,3,3,6,6,3, 2,5,2,2,5,2, 1,1,1,1, 0,0
            };
  sn=20;
  ArrangeDiscretedPoint(p,w,wm);
  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();
};

// Rectangular-type Double-Gaussian Tchebyshev

void SNRZQuadrature::PutRectangularDPnTn(int dpn,int tn)
{
  if(dpn<=0){
    cout<<"Error in PutDPnEW.\n";
    cout<<"You should set dpn>0.\n";
    exit(0);
  };
  if(tn<1){
    cout<<"Error in PutDPnTn.\n";
    cout<<"You should put Tn>0.\n";
    exit(0);
  };

  sn=dpn*2;
  sizex.resize(sn);
  for(int i=0;i<sn;i++){
    sizex[i]=tn*2;
  };

  PutSnnum(sn*tn*2);

  real *tmp2=new real[dpn];
  real *ep=new real[dpn];
  real *w=new real[dpn];
  RootOfLegendre(dpn,tmp2);
  for(int i=0;i<dpn;i++){
    real temp=tmp2[dpn-1-i];
    ep[i]=(1.+temp)*0.5;
    w[i]=WeightOfGaussian(dpn,temp)/tn;
  };

  int index=0;
  for(int i=0;i<dpn;i++){
    real xitmp=-ep[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w[i];
      index++;
    };
  };
  for(int i=0;i<dpn;i++){
    real xitmp=ep[dpn-1-i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w[i];
      index++;
    };
  };

  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();

  delete [] tmp2;
  delete [] ep;
  delete [] w;
};

void SNRZQuadrature::PutBiasedRectangularDPnTn(int dpn1,int dpn2,int tn,real mu_b)
{
  if(dpn1<=0||dpn2<=0){
    cout<<"# Error in SNRZQuadrature::PutBiasedRectangularDPnTn.\n";
    cout<<"# You should set dpn1>0 or dpn2>0.\n";
    exit(0);
  };
  if(tn<1){
    cout<<"# Error in SNRZQuadrature::PutBiasedRectangularDPnTn.\n";
    cout<<"# You should set tn>0.\n";
    exit(0);
  };

  sn=dpn1+dpn2;
  sizex.resize(sn);
  for(int i=0;i<sn;i++){
    sizex[i]=tn*2;
  };

  PutSnnum(sn*tn*2);

  real *tmp=new real[dpn1];
  real *ep=new real[dpn1];
  real *w=new real[dpn1];
  RootOfLegendre(dpn1,tmp);
  for(int i=0;i<dpn1;i++){
    real temp=(1.+tmp[i])*0.5;
    ep[i]=-1.+(mu_b+1.)*temp;
    w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*(mu_b+1.);
    //w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*(1.-0.03971);
    //w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*(1.-0.0506145);
  };

  real *tmp2=new real[dpn2];
  real *ep2=new real[dpn2];
  real *w2=new real[dpn2];
  RootOfLegendre(dpn2,tmp2);
  for(int i=0;i<dpn2;i++){
    real temp=(1.+tmp2[i])*0.5;
    ep2[i]=mu_b+(1.-mu_b)*temp;
    w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*(1.-mu_b);
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*0.0506145;
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*0.03971;
  };

  int index=0;
  for(int i=0;i<dpn1;i++){
    real xitmp=ep[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w[i];
      index++;
    };
  };
  for(int i=0;i<dpn2;i++){
    real xitmp=ep2[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w2[i];
      index++;
    };
  };

  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();

  delete [] tmp;
  delete [] ep;
  delete [] w;
  delete [] tmp2;
  delete [] ep2;
  delete [] w2;
};

void SNRZQuadrature::PutBiasedRectangularDPnTn(int dpn1,int dpn2,int dpn3,int tn,real mu_b)
{
  // + S3 quadrature +
  // 
  //  mu[]   ={-0.774597, 0.. 0.774597};
  //  omega[]={0.555556, 0.888889, 0.555556};
  //
  // So the boundary should be mu=0.549194,
  // but this quadrature does not result in accurate integration on spherical harmonics.

  if(dpn1<=0||dpn2<=0||dpn3<=0){
    cout<<"# Error in SNRZQuadrature::PutBiasedRectangularDPnTn.\n";
    cout<<"# You should set dpn1>0, dpn2>0 and dpn3>0.\n";
    exit(0);
  };
  if(tn<1){
    cout<<"# Error in SNRZQuadrature::PutBiasedRectangularDPnTn.\n";
    cout<<"# You should set tn>0.\n";
    exit(0);
  };

  sn=dpn1+dpn2+dpn3;
  sizex.resize(sn);
  for(int i=0;i<sn;i++){
    sizex[i]=tn*2;
  };

  PutSnnum(sn*tn*2);

  // [-1,-mu_b]
  real *tmp=new real[dpn1];
  real *ep=new real[dpn1];
  real *w=new real[dpn1];
  RootOfLegendre(dpn1,tmp);
  for(int i=0;i<dpn1;i++){
    real tt=(1.+tmp[i])*0.5; // [0,1]
    ep[i]=-1.+(-mu_b+1.)*tt;
    w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*(-mu_b+1.);
    //w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*0.555556;
    //w[i]=WeightOfGaussian(dpn1,tmp[i])/tn*0.101229;
    //cout<<ep[i]<<" "<<w[i]<<"\n";
  };

  // [-mu_b,+mu_b]
  real *tmp2=new real[dpn2];
  real *ep2=new real[dpn2];
  real *w2=new real[dpn2];
  RootOfLegendre(dpn2,tmp2);
  for(int i=0;i<dpn2;i++){
    real tt=(1.+tmp2[i])*0.5;
    ep2[i]=-mu_b+(mu_b+mu_b)*tt;
    w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*(mu_b+mu_b);
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*0.888889;
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*(2.-0.101229*2.);
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*0.03971;
    //cout<<ep2[i]<<" "<<w2[i]<<"\n";
  };

  // [+mu_b,+1]
  real *tmp3=new real[dpn3];
  real *ep3=new real[dpn3];
  real *w3=new real[dpn3];
  RootOfLegendre(dpn3,tmp3);
  for(int i=0;i<dpn3;i++){
    real tt=(1.+tmp3[i])*0.5;
    ep3[i]=mu_b+(1.-mu_b)*tt;
    w3[i]=WeightOfGaussian(dpn3,tmp3[i])/tn*(1.-mu_b);
    //w3[i]=WeightOfGaussian(dpn3,tmp3[i])/tn*0.555556;
    //w3[i]=WeightOfGaussian(dpn3,tmp3[i])/tn*0.101229;
    //w2[i]=WeightOfGaussian(dpn2,tmp2[i])/tn*0.03971;
    //cout<<ep3[i]<<" "<<w3[i]<<"\n";
  };

  int index=0;
  for(int i=0;i<dpn1;i++){
    real xitmp=ep[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w[i];
      index++;
    };
  };
  for(int i=0;i<dpn2;i++){
    real xitmp=ep2[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w2[i];
      index++;
    };
  };
  for(int i=0;i<dpn3;i++){
    real xitmp=ep3[i];
    real fac=sqrt(1.-xitmp*xitmp);
    for(int j=0;j<tn*2;j++){
      mu[index]=fac*cos(PI-PI*0.5/(tn*2)*(2*j+1));
      xi[index]=xitmp;
      omega[index]=w3[i];
      index++;
    };
  };

  WeightNormalize();
  CalMC();
  CalValue();
  CalXYref();

  delete [] tmp;
  delete [] ep;
  delete [] w;
  delete [] tmp2;
  delete [] ep2;
  delete [] w2;
  delete [] tmp3;
  delete [] ep3;
  delete [] w3;
};

void SNRZQuadrature::ArrangeDiscretedPoint(real *p,real *w,int *wm)
{
  sizex.resize(sn);
  int index=0;
  for(int i=0;i<sn/2;i++){
    sizex[index]=(i+1)*2;
    index++;
  };
  for(int i=0;i<sn/2;i++){
    sizex[index]=(sn-i*2);
    index++;
  };

  PutSnnum((sn/2+1)*sn);

  int id=0;
  for(int i=0;i<sn/2;i++){
    for(int j=0;j<i+1;j++){
      mu[id]=-1.*p[i-j];
      xi[id]=-1.*p[sn/2-1-i];
      id++;
    };
    for(int j=0;j<i+1;j++){
      mu[id]=p[j];
      xi[id]=-1.*p[sn/2-1-i];
      id++;
    };
  };
  for(int i=0;i<sn/2;i++){
    for(int j=0;j<sn/2-i;j++){
      mu[id]=-1.*p[sn/2-i-1-j];
      xi[id]=p[i];
      id++;
    };
    for(int j=0;j<sn/2-i;j++){
      mu[id]=p[j];
      xi[id]=p[i];
      id++;
    };
  };

  for(int i=0;i<snnum;i++){
    omega[i]=w[wm[i]];
  };

};

real SNRZQuadrature::GetMoment(int m, int is)
{
  return val[is][m];
};

real SNRZQuadrature::GetMoment(int ll, int mm, int is)
{
  int ind=0;
  for(int l=0;l<=pl;l++){
    for(int m=0;m<=l;m++){
      if(l==ll&&m==mm)return val[is][ind];
      ind++;
    };
  };
  cout<<"# Error in SNRZQuadrature::GetMoment.\n";
  cout<<"# Required moment cannot be found.\n";
  exit(0);
};

real SNRZQuadrature::GetMoment(int m, real i)
{
  return Legendre(m,i);
};

real SNRZQuadrature::GetMoment(int l,int m,real mu,real xi)
{
  return Harmonics(l,m,mu,xi,0.);
};

void SNRZQuadrature::CalMC()
{
  int index=0;
  for(int i=0;i<sn;i++){
    real sum=0.;
    for(int j=0;j<sizex[i];j++){
      mc[index]=sum;
      sum-=omega[index]*mu[index];
      index++;
    };
  };
};

void SNRZQuadrature::CalValue()
{
  for(int i=0;i<snnum;i++){
    int ind=0;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=l;m++){ // specific for 2D
	real eata=sqrt(1.-mu[i]*mu[i]-xi[i]*xi[i]);
	val[i][ind]=Harmonics(l,m,mu[i],xi[i],eata);
	ind++;
      };
    };
  };
};

void SNRZQuadrature::WeightNormalize()
{
  real tot=0.;
  for(int i=0;i<snnum;i++){
    tot+=omega[i];
  };
  for(int i=0;i<snnum;i++){
    omega[i]*=2./tot;
  };
};

void SNRZQuadrature::show_self()
{
  if(exist){
    cout<<"#\n";
    cout<<"# Angular quadrature information\n";
    cout<<"#\n";
    cout<<"# Total number of angular directions : "<<snnum<<"\n";
    cout<<"#\n";
    cout<<"# ID    Omega_r        Omega_z        Omega_{r_perp}Weight\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    for(int i=0;i<snnum;i++){
      if(i<10){
	cout<<"   ";
      }else{
	cout<<"  ";
      };
      cout<<i<<"   ";
      if(mu[i]>0.)cout<<" ";
      cout<<mu[i]<<"   ";
      if(xi[i]>0.)cout<<" ";
      cout<<xi[i]<<"   ";
      real tmp=sqrt(1.-mu[i]*mu[i]-xi[i]*xi[i]);
      if(tmp>0.)cout<<" ";
      cout<<tmp<<"   ";
      cout<<omega[i];
      cout<<"\n";
    };
    /*
    cout<<"sn : "<<snnum<<"\n";
    for(int i=0;i<snnum;i++){
      cout<<i<<" : "<<mu[i]<<"  "<<xi[i]<<"  "<<omega[i]<<"   "<<mc[i]<<"\n";
    };
    */
    /*
    cout<<"\n";
    cout<<" *Value*\n";
    for(int i=0;i<plnum;i++){
      cout<<"plnum="<<i<<"\n";
      for(int j=0;j<snnum;j++){
	cout<<"   "<<j<<" : "<<val[j][i]<<"\n";
      };
    };
    */
  };

  //cout<<"** Reflection **\n";
  for(int i=0;i<snnum;i++){
    //cout<<i<<":  "<<xref[i]<<"&"<<yref[i]<<"\n";
  };
};

void SNRZQuadrature::CalXYref()
{
  for(int i=0;i<snnum;i++){
    real m=mu[i];
    real x=xi[i];
    for(int j=0;j<snnum;j++){
      real m2=mu[j];
      real x2=xi[j];
      if(fabs(m+m2)<0.0001&&
         fabs(x-x2)<0.0001)xref[i]=j;
      if(fabs(m-m2)<0.0001&&
         fabs(x+x2)<0.0001)yref[i]=j;
      if(fabs(m+m2)<0.0001&&
         fabs(x+x2)<0.0001)xyref[i]=j;
    };
  };
};
