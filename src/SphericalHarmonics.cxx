#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "SphericalHarmonics.h"

using namespace std;

real Kaijo(int m)
{
  real ret=1.;
  for(int i=m;i>1;i--){
    ret*=i;
  };
  return ret;
};

real AssociatedLegendre(int l,int m,real mu)
{
  if(m<0){
    int mm=-m;
    return pow(-1.,mm)*Kaijo(l-mm)/Kaijo(l+mm)*AssociatedLegendre(l,mm,mu);
  };

  real pmm=1.0;
  real x=mu;
  if(m>0){
    real somx2=sqrt((1.0-x)*(1.0+x));
    real fact=1.;
    for(int i=1;i<=m;++i){
      pmm*=(-fact)*somx2;
      fact+=2.0;
    };
  };
  if(l==m)return pmm;
  real pmmp1=x*(2.*m+1.)*pmm;
  if(l==m+1)return pmmp1;
  real pll=0.;
  for(int ll=m+2;ll<=l;++ll){
    pll=((2.0*ll-1.0)*x*pmmp1-(ll+m-1.0)*pmm)/(ll-m);
    pmm=pmmp1;
    pmmp1=pll;
  }
  return pll;
  /*

  if(l<abs(m)||l>=6){
    cout<<"Error in AssociatedLegendre : l="<<l<<",m="<<m<<"\n";
    exit(0);
  };


  switch(l){
  case 0:
    if(m==0) return 1.;
    break;
  case 1:
    if(m==0){return mu;}
    else{return sqrt(1-mu*mu);};
    break;
  case 2:
    if(m==0){return 0.5*(3*mu*mu-1);}
    else if(m==1){return 3.*sqrt(1-mu*mu)*mu;}
    else {return 3.*(1-mu*mu);};
    break;
  case 3:
    if(m==0){return 0.5*(5*mu*mu*mu-3*mu);}
    else if(m==1){return 0.5*sqrt(1-mu*mu)*(15.*mu*mu-3.);}
    else if(m==2){return (1.-mu*mu)*15.*mu;}
    else {return 15.*(1-mu*mu)*sqrt(1-mu*mu);};
    break;
  case 4:
    if(m==0){return 0.125*(35*mu*mu*mu*mu-30*mu*mu+3);}
    else if(m==1){return 2.5*sqrt(1-mu*mu)*mu*(7*mu*mu-3);}
    else if(m==2){return 7.5*(1-mu*mu)*(7*mu*mu-1);}
    else if(m==3){return 105*(1-mu*mu)*sqrt(1-mu*mu)*mu;}
    else {return 105*(1-mu*mu)*(1-mu*mu);};
    break;
  case 5:
    if(m==0){return 0.125*(63*mu*mu*mu*mu*mu-70*mu*mu*mu+15*mu);}
    else if(m==1){return 15./8.*sqrt(1-mu*mu)*(21*mu*mu*mu*mu-14*mu*mu+1);}
    else if(m==2){return 52.5*(1-mu*mu)*mu*(3*mu*mu-1);}
    else if(m==3){return 52.5*(1-mu*mu)*sqrt(1-mu*mu)*(9*mu*mu-1);}
    else if(m==4){return 945*(1-mu*mu)*(1-mu*mu)*mu;}
    else {return 945*(1-mu*mu)*(1-mu*mu)*sqrt(1-mu*mu);};
    break;
  };

  if(m==0)       return Legendre(l,mu);
  if(l==1&&m==1) return sqrt(1-mu*mu);
  if(l==2&&m==1) return 3.*sqrt(1-mu*mu)*mu;
  if(l==2&&m==2) return 3.*(1-mu*mu);
  if(l==3&&m==1) return 0.5*sqrt(1-mu*mu)*(15.*mu*mu-3.);
  if(l==3&&m==2) return (1.-mu*mu)*15.*mu;
  if(l==3&&m==3) return 15.*(1-mu*mu)*sqrt(1-mu*mu); 
  if(l==4){
    if(m==1){return 2.5*sqrt(1-mu*mu)*mu*(7*mu*mu-3);}
    else if(m==2){return 7.5*(1-mu*mu)*(7*mu*mu-1);}
    else if(m==3){return 105*(1-mu*mu)*sqrt(1-mu*mu)*mu;}
    else if(m==4){return 105*(1-mu*mu)*(1-mu*mu);};
  }else if(l==5){
    if(m==1){return 15./8.*sqrt(1-mu*mu)*(21*mu*mu*mu*mu-14*mu*mu+1);}
    else if(m==2){return 52.5*(1-mu*mu)*mu*(3*mu*mu-1);}
    else if(m==3){return 52.5*(1-mu*mu)*sqrt(1-mu*mu)*(9*mu*mu-1);}
    else if(m==4){return 945*(1-mu*mu)*(1-mu*mu)*mu;}
    else if(m==5){return 945*(1-mu*mu)*(1-mu*mu)*sqrt(1-mu*mu);};
  };

  cout<<"Error in AssociatedLegendre in SphericalHarmonics.\n";
  cout<<"    "<<l<<" "<<m<<"\n";
  exit(0);
  */
}

real Legendre(int l,real mui)
{
  if(l==0)return 1.;
  if(l==1)return mui;

  real p1=1.;
  real p2=mui;
  real p3=1.;
  for(int j=0;j<l-1;j++){
    int pl=j+2;
    p3=((2*pl-1)*mui*p2-(pl-1)*p1)/pl;
    p1=p2;
    p2=p3;
  };
  return p3;
};

real Ylm(int l,int m,real mu)
{
  real tmp=sqrt(Kaijo(l-m)/Kaijo(l+m));
  return tmp*pow(-1.,m)*AssociatedLegendre(l,m,mu);
};

real Harmonics(int l,int m,real mu,real et,real xi)
{
  real sint=sqrt(1-mu*mu);
  real cosf=0.;
  real sinf=0.;
  real fai =0.;
  if(sint>0.){
    cosf=et/sint;
    sinf=xi/sint;
    fai =fabs(acos(cosf));
  };
  if(cosf>=1.)fai=0.;
  if(cosf<=-1.)fai=PI;
  if(sinf<0.)fai=-fai;

  if(m>=0){
    real eps=2.;
    if(m==0)eps=1.;
    return sqrt(eps)*Ylm(l,m,mu)*cos(fai*m);
  }
  else{
    real eps=sqrt(2.);
    return eps*Ylm(l,m,mu)*sin(fai*m);
  };
};

void RootOfLegendre(int pl, real *root)
{
  /*
  real pold=0.;
  int iin=0;
  for(int i=-1000000;i<=1000000;i++){
    real vur=i*0.000001;
    real pnow=Legendre(pl,vur);
    if(pnow*pold<0.||pnow==0.){
      root[iin]=vur;
      iin++;
      if(iin>pl){
        cout<<"Error in RootOfLegendre in Spherical Harmonics.\n";
        cout<<"The number of roots is too large.\n";
	for(int j=0;j<iin;j++){
	  cout<<root[j]<<"\n";
	};
        exit(0);
      };
    };
    pold=pnow;
  };
  */

  real delta=1e-5;
  int iter_half_max=10000;
  real eps=1e-15;

  real pold=0.;
  int iin=0;
  int range=int(1/delta);
  for(int i=-range;i<=range;i++){
    real vur=i*delta;
    real pnow=Legendre(pl,vur);
    if(pnow==0.){
      root[iin]=vur;
      iin++;
      if(iin>pl){
        cout<<"# Error in RootOfLegendre in Spherical Harmonics.\n";
        cout<<"# The number of roots is too large.\n";
	for(int j=0;j<iin;j++){
	  cout<<root[j]<<"\n";
	};
        exit(0);
      };
      pnow=pold*-1.;
    }else if(pnow*pold<0.){
      real lp=vur-delta;
      real rp=vur;
      real lv=pold;
      real rv=pnow;
      // (2-bunshi method)
      real np=rp;
      for(int j=0;j<iter_half_max;j++){
        np=(lp+rp)*0.5;
        real nv=Legendre(pl,np);
	if(nv==0.){
	  j=iter_half_max;
	}else if(nv*lv>0.){
	  lp=np;
	  lv=nv;
	}else{
	  rp=np;
	  rv=nv;
	};
	if(fabs(lp/rp-1.)<eps){
	  //cout.setf(ios::scientific);
	  //cout.precision(15);
	  //cout<<j<<" "<<fabs(lp/rp-1.)<<" "<<lp<<"\n";
          j=iter_half_max;
	};
      };
      root[iin]=np;
      /*
      cout.setf(ios::showpoint);
      cout.precision(10);
      cout<<vur<<" "<<np<<"\n";
      */
      iin++;
      if(iin>pl){
        cout<<"# Error in RootOfLegendre in Spherical Harmonics.\n";
        cout<<"# The number of roots is too large.\n";
	for(int j=0;j<iin;j++){
	  cout<<root[j]<<"\n";
	};
        exit(0);
      };
    };
    pold=pnow;
  };


  if(iin<pl){
    cout<<"# Error in RootOfLegendre in Spherical Harmonics.\n";
    cout<<"# Root could not be found.\n";
    exit(0);
  };

};

real WeightOfGaussian(int pl, real mu)
// pl should be consistent 'pl' of 'RootOfLegendre'
{
  real plmu=Legendre(pl-1,mu);
  return 2.*(1.-mu*mu)/(pl*pl*plmu*plmu);
};

real Chebyshev(int l,real x)
{
  if(l==0){
    return 1.;
  }else if(l==1){
    return x;
  }else{
    return 2.*x*Chebyshev(l-1,x)-Chebyshev(l-2,x);
  };
};

real ChebyshevSecondKind(int l,real x)
{
  if(l==0){
    return 1.;
  }else if(l==1){
    return 2.*x;
  }else{
    return 2.*x*ChebyshevSecondKind(l-1,x)-ChebyshevSecondKind(l-2,x);
  };
};
