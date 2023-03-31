#include <cstdlib>
#include "Quadrature.h"

using namespace std;

void Quadrature::Initialize(int dim, int p)
{
  Ndim=dim;
  exist=false;
  pl=p;
  int i=0;
  for(int j=0;j<=p;j++){
    i+=j*2+1;
  };
  plnum=i;
};

Quadrature::~Quadrature()
{
}

void Quadrature::PutSN(int sninp)
{
  sn=sninp;
  omega.resize(sn,0.);
  mu.resize(sn,0.);
  eata.resize(sn,0.);
  xi.resize(sn,0.);
  xref.resize(sn,0);
  yref.resize(sn,0);
  zref.resize(sn,0);
  xyzref.resize(sn,0);
  val.resize(sn);
  for(int i=0;i<sn;i++){
    val[i].resize(plnum,0.);
  };
  exist=true;
};


void Quadrature::CalXYZref()
{
  for(int i=0;i<sn;i++){
    real m=mu[i];
    real e=eata[i];
    real x=xi[i];
    for(int j=0;j<sn;j++){
      real m2=mu[j];
      real e2=eata[j];
      real x2=xi[j];
      if(fabs(m+m2)<0.0001&&
         fabs(e-e2)<0.0001&&
         fabs(x-x2)<0.0001)xref[i]=j;
      if(fabs(m-m2)<0.0001&&
         fabs(e+e2)<0.0001&&
         fabs(x-x2)<0.0001)yref[i]=j;
      if(fabs(m-m2)<0.0001&&
         fabs(e-e2)<0.0001&&
         fabs(x+x2)<0.0001)zref[i]=j;
      if(Ndim==3){
        if(fabs(m+m2)<0.0001&&
           fabs(e+e2)<0.0001&&
           fabs(x+x2)<0.0001)xyzref[i]=j;
      };
      if(Ndim==2){
        if(fabs(m+m2)<0.0001&&
           fabs(e+e2)<0.0001)xyzref[i]=j;
      };
    };
  };
};

real Quadrature::GetMoment(int ll,int mm,int is)
{
  int ind=0;
  for(int l=0;l<=pl;l++){
    int lst=-l;
    if(Ndim==2)lst=0;
    for(int m=lst;m<=l;m++){
      if(l==ll&&m==mm)return val[is][ind];
      ind++;
    };
  };
  cout<<"# Error in Quadrature::GetMoment.\n";
  cout<<"# Required moment cannot be found.\n";
  exit(0);
};

void Quadrature::CalDiscretedSphericalHarmonics()
{
  for(int i=0;i<sn;i++){
    int ind=0;
    for(int l=0;l<=pl;l++){
      int lst=-l;
      if(Ndim==2)lst=0;
      for(int m=lst;m<=l;m++){
	val[i][ind]=Harmonics(l,m,mu[i],eata[i],xi[i]);
	if(l==4&&m==0){
	  cout.setf(ios::scientific);
	  cout.precision(10);
	  //	  cout<<mu[i]<<" "<<val[i][ind]<<"\n";
	};
	ind++;
      };
    };
  };
};

void Quadrature::WeightNormalize(real sum)
{
  real tot=0.;
  for(int i=0;i<sn;i++){
    tot+=omega[i];
  };
  real coef=sum/tot;
  for(int i=0;i<sn;i++){
    omega[i]*=coef;
  };
};

void Quadrature::show_self()
{
  if(exist){
    cout<<"sn : "<<sn<<"\n";
    for(int i=0;i<sn;i++){
      cout<<mu[i]<<"  "<<eata[i]<<" "<<xi[i]<<" :  "<<omega[i];
      cout<<"("<<xref[i]<<"&"<<yref[i]<<"&"<<zref[i]<<")\n";
    };
  };
};

void Quadrature::show_self_konno()
{
  if(exist){
    for(int i=0;i<sn/8;i++){
      cout.setf(ios::scientific);
      cout.precision(8);
      cout<<mu[i]<<"  "<<eata[i]<<" "<<xi[i]<<"   "<<omega[i]<<"\n";
    };
  };
};

void Quadrature::PutData(real *mui,real *eai,real *xii,real *w)
{
  for(int i=0;i<sn;i++){
    mu[i]=mui[i];
    eata[i]=eai[i];
    xi[i]=xii[i];
    omega[i]=w[i];
  };
};

void Quadrature::ArrangeDiscretedPoint(real pp,real *mui,real *eai,real *xii)
{
  int id=0;
  for(int z=0;z<(Ndim-1);z++){
    for(int y=0;y<2;y++){
      for(int x=0;x<2;x++){
        for(int i=0;i<pp;i++){
	  if(x==0){mu[id]=-mui[i];}
	  else{mu[id]=mui[i];};
	  if(y==0){eata[id]=-eai[i];}
	  else{eata[id]=eai[i];};
	  if(z==0){xi[id]=-xii[i];}
	  else{xi[id]=xii[i];};
	  omega[id]=omega[i];
	  id++;
	};
      };
    };
  };
};

void Quadrature::Rotation(real theta,string ax)
{
  real cs=cos(theta);
  real si=sin(theta);
  if(ax=="x"){
    for(int i=0;i<sn;i++){
      real new_mu=cs*eata[i]-si*xi[i];
      real new_et=si*eata[i]+cs*xi[i];
      if(new_mu*eata[i]<0.||new_et*xi[i]<0.){
	cout<<"Quadrature::Rotation ... Sign change !.\n";
      };
      eata[i]=new_mu;
      xi[i]=new_et;
    };
  };
  if(ax=="y"){
    for(int i=0;i<sn;i++){
      real new_mu=cs*xi[i]-si*eata[i];
      real new_et=si*xi[i]+cs*eata[i];
      if(new_mu*xi[i]<0.||new_et*eata[i]<0.){
	cout<<"Quadrature::Rotation ... Sign change !.\n";
      };
      xi[i]=new_mu;
      eata[i]=new_et;
    };
  };
  if(ax=="z"){
    for(int i=0;i<sn;i++){
      real new_mu=cs*mu[i]-si*eata[i];
      real new_et=si*mu[i]+cs*eata[i];
      if(new_mu*mu[i]<0.||new_et*eata[i]<0.){
	cout<<"Quadrature::Rotation ... Sign change !.\n";
      };
      mu[i]=new_mu;
      eata[i]=new_et;
    };
  };
  CalDiscretedSphericalHarmonics();
  CalXYZref();
};

void Quadrature::CheckOrthogonalityTo00(real eps)
{
  int index=0;
  for(int l=0;l<=pl;l++){
    int is=-l;
    if(Ndim==2)is=0;
    for(int m=is;m<=l;m++){
      real sum=0.;
      for(int i=0;i<sn;i++){
        sum+=GetOmega(i)*GetMoment(index,i);
      };
      index++;
      if(l!=0&&fabs(sum)>eps){
	cout<<"!! Caution !!\n";
	cout<<" This quadrature set does NOT satisfy the orthogonal relation to (0,0).\n";
	cout<<"   l="<<l<<", m="<<m<<"  ("<<fabs(sum)<<")\n";
      };
    };
  };
};

// +++ Various angular-quadrature set +++

// Level-symmetric quadrature set

void Quadrature::PutLevelSymmetric(int sn)
{
  // Conventional level symmetric quadrature sets to preserve even-moment.
  // Note that there are other "even-moment" quadratures; Lee's one for example.
  // Please see the CRC handbook written by Alcouffe, page 390.

  switch(sn){
  case 2:
    PutLSS2();
    break;
  case 4:
    PutLSS4();
    break;
  case 6:
    PutLSS6();
    break;
  case 8:
    PutLSS8();
    break;
  case 10:
    PutLSS10();
    break;
  case 12:
    PutLSS12();
    break;
  case 16:
    PutLSS16();
    break;
  case 20:
    PutLSS20();
    break;
  default:
    cout<<"Not coded : "<<sn<<"-order level-symmetric quadrature.\n";
    exit(0);
  };
};

void Quadrature::PutLebedev()
{
  // precision-3
  /*
  int pnt=6;
  real data[]={ 
    0.000000000000000,    90.000000000000000,     0.166666666666667,
    180.000000000000000,    90.000000000000000,     0.166666666666667,
    90.000000000000000,   90.000000000000000,     0.166666666666667,
    -90.000000000000000,    90.000000000000000,     0.166666666666667,
    90.000000000000000,     0.000000000000000,     0.166666666666667,
    90.000000000000000,   180.000000000000000,     0.166666666666667
  };
  */

  // precision-7
  int pnt=26;
  real data[]={
    0.000000000000000   , 90.000000000000000,     0.047619047619048,
    180.000000000000000    ,90.000000000000000,     0.047619047619048,
    90.000000000000000    ,90.000000000000000  ,   0.047619047619048,
    -90.000000000000000   , 90.000000000000000  ,   0.047619047619048,
    90.000000000000000   ,  0.000000000000000   ,  0.047619047619048,
    90.000000000000000   ,180.000000000000000   ,  0.047619047619048,
    90.000000000000000   , 45.000000000000000   ,  0.038095238095238,
    90.000000000000000   ,135.000000000000000   ,  0.038095238095238,
    -90.000000000000000    ,45.000000000000000  ,   0.038095238095238,
    -90.000000000000000   ,135.000000000000000  ,   0.038095238095238,
    0.000000000000000   , 45.000000000000000   ,  0.038095238095238,
    0.000000000000000   ,135.000000000000000   ,  0.038095238095238,
    180.000000000000000    ,45.000000000000000 ,    0.038095238095238,
    180.000000000000000   ,135.000000000000000 ,    0.038095238095238,
    45.000000000000000    ,90.000000000000000  ,   0.038095238095238,
    -45.000000000000000    ,90.000000000000000 ,    0.038095238095238,
    135.000000000000000    ,90.000000000000000 ,    0.038095238095238,
    -135.000000000000000   , 90.000000000000000 ,    0.038095238095238,
    45.000000000000000   , 54.735610317245346 ,    0.032142857142857,
    45.000000000000000   ,125.264389682754654  ,   0.032142857142857,
    -45.000000000000000  ,  54.735610317245346 ,    0.032142857142857,
    -45.000000000000000,   125.264389682754654 ,    0.032142857142857,
    135.000000000000000   , 54.735610317245346 ,    0.032142857142857,
    135.000000000000000  , 125.264389682754654 ,    0.032142857142857,
    -135.000000000000000  ,  54.735610317245346 ,    0.032142857142857,
    -135.000000000000000,   125.264389682754654  ,   0.032142857142857
  };

  /*
  // precision-15
  int pnt=86;
  real data[]={
    0.000000000000000,    90.000000000000000    , 0.011544011544012,
    180.000000000000000,    90.000000000000000  ,   0.011544011544012,
    90.000000000000000,    90.000000000000000   ,  0.011544011544012,
    -90.000000000000000,    90.000000000000000  ,   0.011544011544012,
    90.000000000000000  ,   0.000000000000000   ,  0.011544011544012,
    90.000000000000000  , 180.000000000000000   ,  0.011544011544012,
    45.000000000000000  ,  54.735610317245346   ,  0.011943909085856,
    45.000000000000000  , 125.264389682754654   ,  0.011943909085856,
    -45.000000000000000  ,  54.735610317245346  ,   0.011943909085856,
    -45.000000000000000  , 125.264389682754654  ,   0.011943909085856,
    135.000000000000000   , 54.735610317245346  ,   0.011943909085856,
    135.000000000000000  , 125.264389682754654  ,   0.011943909085856,
    -135.000000000000000   , 54.735610317245346 ,    0.011943909085856,
    -135.000000000000000  , 125.264389682754654 ,    0.011943909085856,
    45.000000000000000   , 31.513359490876244   ,  0.011110555710603,
    45.000000000000000  , 148.486640509123760   ,  0.011110555710603,
    -45.000000000000000   , 31.513359490876244  ,   0.011110555710603,
    -45.000000000000000   ,148.486640509123760  ,   0.011110555710603,
    135.000000000000000   , 31.513359490876244   ,  0.011110555710603,
    135.000000000000000   ,148.486640509123760  ,   0.011110555710603,
    -135.000000000000000    ,31.513359490876244 ,    0.011110555710603,
    -135.000000000000000  , 148.486640509123760 ,    0.011110555710603,
    66.561222392026522   , 68.308874114180426   ,  0.011110555710603,
    -66.561222392026522   , 68.308874114180426  ,   0.011110555710603,
    66.561222392026522   ,111.691125885819588   ,  0.011110555710603,
    -66.561222392026522   ,111.691125885819588  ,   0.011110555710603,
    113.438777607973492   , 68.308874114180426  ,   0.011110555710603,
    -113.438777607973492   , 68.308874114180426 ,    0.011110555710603,
    113.438777607973492   ,111.691125885819588  ,   0.011110555710603,
    -113.438777607973492  , 111.691125885819588 ,    0.011110555710603,
    23.438777607973478    ,68.308874114180426   ,  0.011110555710603,
    156.561222392026536  ,  68.308874114180426  ,   0.011110555710603,
    23.438777607973478  , 111.691125885819588   ,  0.011110555710603,
    156.561222392026536  , 111.691125885819588  ,   0.011110555710603,
    -23.438777607973478  ,  68.308874114180426  ,   0.011110555710603,
    -156.561222392026536  ,  68.308874114180426 ,    0.011110555710603,
    -23.438777607973478  , 111.691125885819588  ,   0.011110555710603,
    -156.561222392026536 ,  111.691125885819588 ,    0.011110555710603,
    45.000000000000000  ,  79.101860732695940   ,  0.011876501294537,
    45.000000000000000 ,  100.898139267304060   ,  0.011876501294537,
    -45.000000000000000 ,   79.101860732695940  ,   0.011876501294537,
    -45.000000000000000  , 100.898139267304060  ,   0.011876501294537,
    135.000000000000000  ,  79.101860732695940  ,   0.011876501294537,
    135.000000000000000  , 100.898139267304060  ,   0.011876501294537,
    -135.000000000000000  ,  79.101860732695940 ,    0.011876501294537,
    -135.000000000000000 ,  100.898139267304060 ,    0.011876501294537,
    15.231635441931418  ,  46.024237785326378   ,  0.011876501294537,
    -15.231635441931418  ,  46.024237785326378  ,   0.011876501294537,
    15.231635441931418 ,  133.975762214673637   ,  0.011876501294537,
    -15.231635441931418 ,  133.975762214673637  ,   0.011876501294537,
    164.768364558068583  ,  46.024237785326378  ,   0.011876501294537,
    -164.768364558068583  ,  46.024237785326378 ,    0.011876501294537,
    164.768364558068583  , 133.975762214673637  ,   0.011876501294537,
    -164.768364558068583 ,  133.975762214673637 ,    0.011876501294537,
    74.768364558068569  ,  46.024237785326378   ,  0.011876501294537,
    105.231635441931431  ,  46.024237785326378  ,   0.011876501294537,
    74.768364558068569  , 133.975762214673637   ,  0.011876501294537,
    105.231635441931431  , 133.975762214673637  ,   0.011876501294537,
    -74.768364558068569  ,  46.024237785326378  ,   0.011876501294537,
    -105.231635441931431  ,  46.024237785326378 ,    0.011876501294537,
    -74.768364558068569  , 133.975762214673637  ,   0.011876501294537,
    -105.231635441931431  , 133.975762214673637 ,    0.011876501294537,
    68.022464238570777   , 90.000000000000000   ,  0.011812303746904,
    -68.022464238570777   , 90.000000000000000  ,   0.011812303746904,
    111.977535761429237   , 90.000000000000000  ,   0.011812303746904,
    -111.977535761429237   , 90.000000000000000 ,    0.011812303746904,
    21.977535761429227  ,  90.000000000000000   ,  0.011812303746904,
    -21.977535761429227   , 90.000000000000000  ,   0.011812303746904,
    158.022464238570763   , 90.000000000000000  ,   0.011812303746904,
    -158.022464238570763   , 90.000000000000000 ,    0.011812303746904,
    0.000000000000000   , 21.977535761429227    , 0.011812303746904,
    0.000000000000000   ,158.022464238570763   ,  0.011812303746904,
    180.000000000000000  ,  21.977535761429227  ,   0.011812303746904,
    180.000000000000000  , 158.022464238570763   ,  0.011812303746904,
    0.000000000000000  ,  68.022464238570777    , 0.011812303746904,
    0.000000000000000 ,  111.977535761429237    , 0.011812303746904,
    180.000000000000000 ,   68.022464238570777  ,   0.011812303746904,
    180.000000000000000  , 111.977535761429237  ,   0.011812303746904,
    90.000000000000000  ,  21.977535761429227   ,  0.011812303746904,
    90.000000000000000  , 158.022464238570763   ,  0.011812303746904,
    -90.000000000000000  ,  21.977535761429227   ,  0.011812303746904,
    -90.000000000000000  , 158.022464238570763  ,   0.011812303746904,
    90.000000000000000  ,  68.022464238570777  ,   0.011812303746904,
    90.000000000000000  , 111.977535761429237  ,   0.011812303746904,
    -90.000000000000000  ,  68.022464238570777 ,    0.011812303746904,
    -90.000000000000000  , 111.977535761429237 ,    0.011812303746904
  };
  */
  PutSN(pnt);

  for(int i=0;i<pnt;i++){
    real theta=data[i*3]/180.*PI;
    real phi=data[i*3+1]/180.*PI;
    mu[i]=cos(theta)*sin(phi);
    eata[i]=sin(theta)*sin(phi);
    xi[i]=cos(phi);
    omega[i]=data[i*3+2];
  };
  
  WeightNormalize();
  CalValue();
};

void Quadrature::LevelSymmetric(int sn,real *p,real *w,int *wm)
{
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int id=0;
  for(int i=0;i<sn/2;i++){
    int is;
    if(i<sn/4){is=i;}else{is=sn/4-1-(i-sn/4);};
    for(int j=0;j<i+1;j++){
      mui[id]=p[i-j];
      eai[id]=p[j];
      xii[id]=p[sn/2-1-i];
      id++;
    };
  };
  for(int i=0;i<tmp;i++){
    omega[i]=w[wm[i]];
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutLSS2()
{
  //real p[]={0.5773503};
  real p[]={0.5773502691896257}; // From Endo
  real w[]={1.0};
  int wm[]={0};
  LevelSymmetric(2,p,w,wm);
};

void Quadrature::PutLSS4()
{
  // (TWODANT or TRITAC)
  //real p[]={0.30163878, 0.90444905};
  //real w[]={0.33333333};
  // (Kobayashi's text & EndoSN) 
  /*
  real p[]={0.35002117458, 0.86889030072};
  real w[]={0.33333333333};
  */
  // (Endo)
  real p[]={0.3500211745815407, 0.8688903007222012};
  real w[]={0.33333333333};

  int wm[]={0,0,0};
  LevelSymmetric(4,p,w,wm);

  // +++ Galerkin matrix
  if(pl>=5){

  gmatrix.put_yx(sn,sn);

  int l_list[]={0,1,1,3,4, 3,3,4, 5, 1,  3,3,2, 2, 4,  3, 4, 4, 4, 2, 5,2,2,3};
  int m_list[]={0,0,1,1,1, 2,3,3,-4,-1, -3,0,2,-1,-1, -2,-2,-4,-3,-2, -2,1,0,-1};
  for(int i=0;i<sn;i++){
    for(int j=0;j<sn;j++){
      int ll=l_list[j];
      int mm=m_list[j];
      gmatrix.put_data(i,j,GetMoment(ll,mm,i));
    };
  };

  gmatrix_inv=gmatrix;
  gmatrix_inv.DoInverse();


  l_array_gmatrix.resize(sn);
  for(int i=0;i<sn;i++){
    l_array_gmatrix[i]=l_list[i];
  };

  };

};

void Quadrature::PutLSS6()
{
  // (Endo)
  real p[]={0.2666354015167047, 0.6815077265365469, 0.9261809355174892};
  real w[]={0.02201576635792293, 0.01965090030874374};

  int wm[]={0, 1,1, 0,1,0};
  LevelSymmetric(6,p,w,wm);
};

void Quadrature::PutLSS8()
{
  // (TWODANT)
  //real p[]={0.19232747, 0.57735026, 0.79352176, 0.96229947};
  //real w[]={0.029197117, 0.023313807, 0.0225258};
  // (TRITAC)
  //real p[]={0.19232747, 0.57735027, 0.79352178, 0.96229948};
  //real w[]={0.11678847, 0.09325523, 0.0901032};
  // (Kobayashi's text & EndoSN)
  //real p[]={0.2182179, 0.5773503, 0.7867958, 0.9511897};
  //real w[]={0.1209877, 0.0907407, 0.0925926};
  // (Endo)
  real p[]={0.2182178902359924, 0.5773502691896257, 0.7867957924694432, 0.9511897312113419};
  real w[]={0.01512345679012345, 0.01134259259259259, 0.01157407407407407};
  int wm[]={0, 1,1, 1,2,1, 0,1,1,0};
  LevelSymmetric(8,p,w,wm);

  // +++ Galerkin matrix
  if(pl>=9){

  gmatrix.put_yx(sn,sn);

  int l_list[]={
  0,1,1,1,2,2,2,2,3,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,6,5,6,6,5,6,6,7,7,7,8,6,6,5,5,6,5,6,6,6,7,4,7,7,7,7,7,8,7,7,8,6,7,7,7,8,5,8,7,8,8,8,8,6,8,8,8,9,9,9,9,
  };
  int m_list[]={
  0,-1,0,1,-2,-1,0,1,-3,2,-2,-1,0,1,2,3,-3,1,-1,3,-4,-2,0,2,-5,-4,0,1,2,0,3,-6,-5,-2,1,3,-6,0,1,1,-2,-1,-1,4,2,5,4,5,6,-7,4,-5,-4,2,-2,-1,-7,4,6,-5,-4,5,3,7,-8,-3,-6,-3,-4,-3,-2,-1,-3,3,5,7,-8,-6,-4,-2,
  };
  for(int i=0;i<sn;i++){
    for(int j=0;j<sn;j++){
      int ll=l_list[j];
      int mm=m_list[j];
      gmatrix.put_data(i,j,GetMoment(ll,mm,i));
    };
  };

  gmatrix_inv=gmatrix;
  gmatrix_inv.DoInverse();


  l_array_gmatrix.resize(sn);
  for(int i=0;i<sn;i++){
    l_array_gmatrix[i]=l_list[i];
  };

  };
};

void Quadrature::PutLSS10()
{
  // (TWODANT or TRITAC)
  //real p[]={0.16962228,0.50714192,0.69686020,0.84500612,0.97080202};
  //real w[]={0.089842043, 0.067288705, 0.055780071, 0.053133809};
  //real p[]={0.189321326478, 0.5088817555826, 0.694318887594, 0.8397599622366, 0.9634909811105};
  //real w[]={0.011162893498, 0.0090661439640, 0.00563047092955, 0.00674101431098};
  // (Endo)
  real p[]={0.1893213264780105, 0.5088817555826189, 0.6943188875943843, 0.8397599622366848, 0.9634909811104685};
  real w[]={0.01116289349804459, 0.009066143964045693, 0.005630470929551081, 0.006741014310979617};
  int wm[]={0, 1,1, 2,3,2, 1,3,3,1, 0,1,2,1,0};
  LevelSymmetric(10,p,w,wm);
};

void Quadrature::PutLSS12()
{
  // From TWODANT
  //real p[]={0.15395746, 0.45769112, 0.6286966, 0.76225828,
  //           0.87568027, 0.97600932};
  //real w[]={0.07332178, 0.0526674, 0.04161495, 0.03895667,
  //           0.03249018 };
  // From ENDOSN
  /*
  real p[]={0.1672126528, 0.4595476346, 0.6280190966, 0.760021014834,
            0.8722705439, 0.97163771925};
  real w[]={0.008845323746, 0.00698513769, 0.0046672092199, 0.0062852376325,
             0.0032314114569688 };
  */
  // (Endo)
  real p[]={0.1672126528227133, 0.4595476346425947, 0.6280190966421309, 0.760021014833664,
            0.8722705430257215, 0.9716377192513583};
  real w[]={0.008845323746261383, 0.006985137695611102, 0.004667209219853574, 0.006285237632507141,
             0.0032314114569688 };
  int wm[]={0, 1,1, 2,3,2, 2,4,4,2, 1,3,4,3,1, 0,1,2,2,1,0};
  LevelSymmetric(12,p,w,wm);

  // +++ Galerkin matrix
  if(pl>=13){

  gmatrix.put_yx(sn,sn);

  int l_list[]={
0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,8,8,8,8,8,9,7,7,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,11,11,11,12,8,9,9,9,9,9,7,7,6,7,9,9,9,9,9,9,7,8,10,10,10,10,10,10,10,10,8,11,10,8,10,11,10,10,10,10,12,11,11,11,11,11,8,11,11,11,8,11,8,8,10,11,11,11,11,11,11,11,11,11,12,10,12,12,12,12,12,12,12,12,12,12,9,12,12,12,12,12,13,13,13,13,13,13,
  };
  int m_list[]={
0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-5,1,-3,3,-1,5,0,-6,2,-4,4,-2,-7,-6,-5,-4,0,1,2,3,4,-8,-7,0,1,2,0,7,-2,-6,-5,3,5,-9,-8,-6,1,2,3,-10,0,1,-9,3,-10,0,1,1,-1,-5,-4,-3,-2,-1,6,-1,6,5,4,5,6,7,8,9,-3,6,-8,-7,-6,-5,-4,-3,2,-1,4,-11,-2,7,4,2,7,6,9,10,-11,-8,-6,-4,-2,6,8,-5,-9,-3,-4,-1,-3,-2,5,3,4,5,-7,7,8,9,10,11,-12,8,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,-7,3,5,7,9,11,-12,-10,-8,-6,-4,-2,
  };

  for(int i=0;i<sn;i++){
    for(int j=0;j<sn;j++){
      int ll=l_list[j];
      int mm=m_list[j];
      gmatrix.put_data(i,j,GetMoment(ll,mm,i));
    };
  };

  gmatrix_inv=gmatrix;
  gmatrix_inv.DoInverse();

  l_array_gmatrix.resize(sn);
  for(int i=0;i<sn;i++){
    l_array_gmatrix[i]=l_list[i];
  };

  };
};

void Quadrature::PutLSS16()
{
 // From Dr. Endo
  real p[]={0.1389568750676416, 0.3922892614447836, 0.5370965613008739, 0.6504264506287802,
            0.7467505736146995, 0.8319965569100706, 0.9092855009437586, 0.9805008790117792};
  real w[]={0.006123404894745294, 0.00516619973373803, 0.002805949711275305, 0.004202330861723087,
             0.003050709853805049, 0.001959237716203034, 0.004615716388076537, 0.0007610204920789216};
  // TWODANT
  //real p[]={0.13344572, 0.39119433, 0.53689687, 0.65075610,
  //           0.74746822, 0.83302700, 0.91058181, 0.98203079};
  //real w[]={0.05415425, 0.03679653, 0.02777273, 0.02494275,
  //           0.02580284, 0.01962325, 0.01879762, 0.01544801};
   int wm[]={0, 1,1, 2,3,2, 4,5,5,4, 4,6,7,6,4, 2,5,7,7,5,2,
            1,3,5,6,5,3,1, 0,1,2,4,4,2,1,0};
  LevelSymmetric(16,p,w,wm);

  if(pl>=17){

  gmatrix.put_yx(sn,sn);


  int l_list[]={
0,1,1,3,2,3,4,5,4,3,6,6,7,5,7,6,7,7,7,7,6,7,7,9,8,9,8,7,4,5,7,9,9,8,9,5,1,4,6,6,4,7,8,7,6,8,9,8,8,8,9,9,9,8,8,9,8,9,10,8,9,10,10,10,5,4,5,10,10,11,10,11,3,8,3,8,6,7,11,7,4,6,11,10,10,6,6,11,11,9,4,12,6,10,10,10,10,11,9,10,11,8,11,11,11,12,12,12,2,12,3,5,12,12,12,12,13,13,12,13,11,13,13,13,14,13,13,14,14,14,14,14,15,10,11,12,12,11,11,12,13,15,15,16,9,9,7,5,12,6,12,5,9,8,12,12,4,11,8,10,12,10,12,12,12,12,12,10,12,10,13,10,13,13,11,13,13,13,13,11,2,10,11,13,11,11,14,13,13,13,14,13,13,13,13,15,14,14,14,14,14,15,14,14,14,14,14,13,16,15,14,9,14,14,15,15,11,9,14,2,14,14,14,14,14,15,13,14,3,15,13,15,12,15,15,15,15,15,15,15,2,11,15,15,15,15,15,15,15,15,15,16,15,15,15,5,15,16,16,16,16,14,16,16,16,16,16,16,16,16,16,16,5,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,
  };
  int m_list[]={
0,0,1,-1,1,-3,0,-4,-4,1,0,-2,-6,3,-7,-1,-4,-2,-3,6,-6,2,4,0,1,-9,-5,3,3,1,0,-6,-7,-3,-4,5,-1,-1,-4,1,-2,-5,0,-1,5,3,-5,2,-6,-1,1,2,-8,-4,6,3,-2,4,-9,5,5,-7,0,2,-2,-3,2,1,3,0,4,2,3,-8,0,4,4,7,-11,1,1,-3,-10,-10,-8,-5,6,-9,1,6,2,1,2,-2,-1,5,7,-8,7,6,-6,-7,3,4,5,-12,-11,-10,0,-9,2,-5,0,2,3,5,-13,-12,-5,-10,-2,0,1,2,-13,-9,3,-8,-14,0,1,3,-14,-6,-7,-7,-3,-1,-4,7,-4,0,1,1,-3,-1,5,-3,-8,3,-6,4,9,7,-2,-1,4,-3,8,9,4,-5,6,8,9,10,11,-4,12,-3,-11,8,-8,-7,6,-5,-2,-3,-1,8,-1,10,7,4,9,10,-11,5,6,7,2,9,10,11,12,-13,-10,-12,-9,-7,-6,2,5,-5,-4,-3,-1,8,-11,-12,7,-2,9,11,-10,-8,-5,8,8,2,10,6,12,13,14,-15,-6,4,-2,-11,13,-9,-4,-7,-6,-5,-4,-3,-2,-1,-2,11,3,15,4,5,6,7,8,9,10,-16,11,13,14,-1,12,-15,-14,-13,-12,-2,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,3,5,7,9,11,13,15,-14,-16,-12,-10,-8,-6,-4,-2,
  };

  for(int i=0;i<sn;i++){
    for(int j=0;j<sn;j++){
      int ll=l_list[j];
      int mm=m_list[j];
      gmatrix.put_data(i,j,GetMoment(ll,mm,i));
    };
  };

  gmatrix_inv=gmatrix;
  gmatrix_inv.DoInverse();

  l_array_gmatrix.resize(sn);
  for(int i=0;i<sn;i++){
    l_array_gmatrix[i]=l_list[i];
  };

  };

};

void Quadrature::PutLSS20()
{
 // From Dr. Endo
  real p[]={0.1206033430392688, 0.3475742923164429, 0.4765192661438829, 0.5773502691896257,
            0.6630204036531319, 0.7388225619100911, 0.8075404016607585, 0.8708525837599884,
            0.9298639389547678, 0.9853474855580162};
  real w[]={0.004627631132711769, 0.00416052706706861, 0.00174587686158126, 0.003635641540060942,
            0.0007789912557568425, 0.003277083750555231, 0.000285942423601845, 0.00454989112761428,
            0.001123824501815679, 0.001123824501815679, 0.0003720086401950338, 0.001369634844032048};
  int wm[]={0, 1,1, 2,5,2, 3,6,6,3, 4,7,9,7,4, 4,8,10,10,8,4, 3,7,10,11,10,7,3,
	    2,6,9,10,10,9,6,2, 1,5,6,7,8,7,6,5,1, 0,1,2,3,4,4,3,2,1,0};
  LevelSymmetric(20,p,w,wm);
};

// Triangular-type Double-Gaussian Tchebyshev

void Quadrature::PutTriangularDPnTn(int dpn, string dir)
{
  if(dpn<=0){
    cout<<"Error in PutDPnEW.\n";
    cout<<"You should set dpn>0.\n";
    exit(0);
  };
  real *tmp2=new real[dpn];
  real *ep=new real[dpn];
  real *w=new real[dpn];
  RootOfLegendre(dpn,tmp2);
  for(int i=0;i<dpn;i++){
    real temp=tmp2[dpn-1-i];
    ep[i]=(1.+temp)*0.5;
    w[i]=WeightOfGaussian(dpn,temp)/(i+1);
  };

  int tmp=(dpn+1)*dpn/2;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  for(int i=0;i<dpn;i++){
    real eat=ep[i];
    real fac=sqrt(1.-eat*eat);
    for(int j=0;j<=i;j++){
      real a=fac*cos(PI*0.5/((i+1)*2)*(2*j+1));
      real b=eat;
      real c=sqrt(1.-a*a-b*b);
      if(dir=="z"){
	mui[index]=a;
	eai[index]=c;
	xii[index]=b;
      } else if(dir=="x"){
	mui[index]=b;
	eai[index]=c;
	xii[index]=a;
      } else {
	mui[index]=a;
	eai[index]=b;
	xii[index]=c;
      };
      omega[index]=w[i];
      index++;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
  delete [] tmp2;
  delete [] ep;
  delete [] w;
};

void Quadrature::PutTriangularPnTn(int pn, string dir)
{
  if(pn<=0){
    cout<<"# Error in Quadrature::PutTriangularPnTn.\n";
    cout<<"# You should set pn>0.\n";
    exit(0);
  };
  real *tmp2=new real[pn*2];
  real *ep=new real[pn];
  real *w=new real[pn];
  RootOfLegendre(pn*2,tmp2);
  for(int i=0;i<pn;i++){
    real temp=tmp2[pn*2-1-i];
    ep[i]=temp;
    w[i]=WeightOfGaussian(pn*2,temp)/(i+1);
  };

  int tmp=(pn+1)*pn/2;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  for(int i=0;i<pn;i++){
    real eat=ep[i];
    real fac=sqrt(1.-eat*eat);
    for(int j=0;j<=i;j++){
      real a=fac*cos(PI*0.5/((i+1)*2)*(2*j+1));
      real b=eat;
      real c=sqrt(1.-a*a-b*b);
      if(dir=="z"){
	mui[index]=a;
	eai[index]=c;
	xii[index]=b;
      } else if(dir=="x"){
	mui[index]=b;
	eai[index]=c;
	xii[index]=a;
      } else {
	mui[index]=a;
	eai[index]=b;
	xii[index]=c;
      };
      omega[index]=w[i];
      index++;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
  delete [] tmp2;
  delete [] ep;
  delete [] w;
};

// Rectangular-type Double-Gaussian Tchebyshev

void Quadrature::PutRectangularDPnTn(int dpn,int tn,string dir)
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

  real *tmp2=new real[dpn];
  real *ep=new real[dpn];
  real *w=new real[dpn];
  RootOfLegendre(dpn,tmp2);
  for(int i=0;i<dpn;i++){
    real temp=tmp2[dpn-1-i];
    ep[i]=(1.+temp)*0.5;
    //cout<<i<<" "<<ep[i]<<"\n";
    w[i]=WeightOfGaussian(dpn,temp)/tn;
  };

  int tmp=dpn*tn;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  for(int i=0;i<dpn;i++){
    real eat=ep[i];
    real fac=sqrt(1.-eat*eat);
    for(int j=0;j<tn;j++){
      real a=fac*cos(PI*0.5/(tn*2)*(2*j+1));
      real b=eat;
      real c=sqrt(1.-a*a-b*b);
      if(dir=="z"){
	mui[index]=a;
	eai[index]=c;
	xii[index]=b;
      } else if(dir=="x"){
	mui[index]=b;
	eai[index]=c;
	xii[index]=a;
      } else {
	mui[index]=a;
	eai[index]=b;
	xii[index]=c;
      };
      omega[index]=w[i];
      index++;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
  delete [] tmp2;
  delete [] ep;
  delete [] w;
};

void Quadrature::PutRectangularPnTn(int pn,int tn,string dir)
{
  if(pn<=0){
    cout<<"# Error in PutRectangularPnTn.\n";
    cout<<"# You should set pn>0.\n";
    exit(0);
  };
  if(tn<1){
    cout<<"# Error in PutRectangularPnTn.\n";
    cout<<"# You should put Tn>0.\n";
    exit(0);
  };

  real *tmp2=new real[pn*2];
  real *ep=new real[pn];
  real *w=new real[pn];
  RootOfLegendre(pn*2,tmp2);
  for(int i=0;i<pn;i++){
    real temp=tmp2[pn*2-1-i];
    ep[i]=temp;
    //cout<<i<<" "<<ep[i]<<"\n";
    w[i]=WeightOfGaussian(pn*2,temp)/tn;
  };

  int tmp=pn*tn;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  for(int i=0;i<pn;i++){
    real eat=ep[i];
    real fac=sqrt(1.-eat*eat);
    for(int j=0;j<tn;j++){
      real a=fac*cos(PI*0.5/(tn*2)*(2*j+1));
      real b=eat;
      real c=sqrt(1.-a*a-b*b);
      if(dir=="z"){
	mui[index]=a;
	eai[index]=c;
	xii[index]=b;
      } else if(dir=="x"){
	mui[index]=b;
	eai[index]=c;
	xii[index]=a;
      } else {
	mui[index]=a;
	eai[index]=b;
	xii[index]=c;
      };
      omega[index]=w[i];
      index++;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
  delete [] tmp2;
  delete [] ep;
  delete [] w;
};

// Triangular-Type Tchebyshev-Tchebyshev

void Quadrature::PutTriangularTn(int x)
{
  int tmpsn=(x+2)*x/8;
  real *mui=new real[tmpsn];
  real *eai=new real[tmpsn];
  real *xii=new real[tmpsn];

  int totsn=tmpsn*(Ndim-1)*4;
  PutSN(totsn);

  real *mu_tmp=new real[x];
  real *ww_tmp=new real[x];
  for(int i=0;i<x;i++){
    mu_tmp[i]=cos(PI/(2*x)*(2*i+1));
    real mu0=cos(PI/x*i);
    real mu1=cos(PI/x*(i+1));
    ww_tmp[i]=mu0-mu1;
  };

  int index=0;
  for(int z=0;z<x/2;z++){
    real om=ww_tmp[z]/(z+1);
    real fac=sqrt(1.-mu_tmp[z]*mu_tmp[z]);
    for(int i=0;i<=z;i++){
      mui[index]=mu_tmp[z];
      eai[index]=fac*mu_tmp[i];
      xii[index]=fac*sqrt(1.-mu_tmp[i]*mu_tmp[i]);
      omega[index]=om;
      index++;
    };
  };

  delete [] mu_tmp;
  delete [] ww_tmp;

  ArrangeDiscretedPoint(tmpsn,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();
  delete [] mui;
  delete [] eai;
  delete [] xii;
};

// Rectangular-type Double-Gaussian Double-Gaussian

void Quadrature::PutRectangularDPnNxN(int x1,int x2,string dir)
{
  if(dir!="x"&&dir!="y"&&dir!="z"){
    cout<<"Error in PutDPnNxN.\n";
    cout<<"dir cannot be accepted.\n";
    exit(0);
  };

  if(x2==-1)x2=x1;

  int tmpsn=x1*x2;
  real *mui=new real[tmpsn];
  real *eai=new real[tmpsn];
  real *xii=new real[tmpsn];

  int totsn=tmpsn*(Ndim-1)*4;
  PutSN(totsn);

  real *tmp1=new real[x1];
  real *ep=new real[x1];
  real *w1=new real[x1];
  RootOfLegendre(x1,tmp1);
  for(int i=0;i<x1;i++){
    real temp=tmp1[x1-1-i];
    ep[i]=(1.+temp)*0.5;
    w1[i]=WeightOfGaussian(x1,temp);
  };

  real *tmp2=new real[x2];
  real *w2=new real[x2];
  RootOfLegendre(x2,tmp2);
  real *phi=new real[x2];
  for(int i=0;i<x2;i++){
    real temp=tmp2[x2-1-i];
    phi[i]=0.5*PI*(1.+temp)*0.5;
    w2[i]=WeightOfGaussian(x2,temp);
  };

  int index=0;
  for(int i=0;i<x1;i++){
    real xi=ep[i];
    for(int j=0;j<x2;j++){
      real mu=sqrt(1.-xi*xi)*cos(phi[j]);
      real et=sqrt(1.-xi*xi)*sin(phi[j]);
      if(dir=="z"){
        mui[index]=mu;
        eai[index]=et;
        xii[index]=xi;
      }else if(dir=="x"){
	mui[index]=xi;
	eai[index]=et;
	xii[index]=mu;
      }else{
	mui[index]=mu;
	eai[index]=xi;
	xii[index]=et;
      };
      omega[index]=w1[i]*w2[j];
      index++;
    };
  };

  delete [] tmp1;
  delete [] ep;
  delete [] w1;
  delete [] tmp2;
  delete [] w2;
  delete [] phi;

  ArrangeDiscretedPoint(tmpsn,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();
  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO4()
{
  int sn=4;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  for(int i=0;i<1;i++){
    if(i==0){
      pinp[0]=0.2958758547680685;
      pinp[1]=0.9082482904638630;
      winp=0.041666666666;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO6()
{
  int sn=6;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  for(int i=0;i<2;i++){
    if(i==0){
      pinp[0]=0.1436334558683300;
      pinp[1]=0.9791521131625265;
      winp=0.01244626137612538;
    };
    if(i==1){
      pinp[0]=0.6897995893236758;
      pinp[1]=0.2198932766998037;
      winp=0.02922040529054129;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO8()
{
  int sn=8;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  real tmp1=0.5773502691896258;
  mui[index]=tmp1;
  eai[index]=tmp1;
  xii[index]=tmp1;
  omega[index++]=0.01724403382698509;

  for(int i=0;i<3;i++){
    if(i==0){
      pinp[0]=0.1155882274614839;
      pinp[1]=0.9865488955670796;
      winp=0.007278598892030222;
    };
    if(i==1){
      pinp[0]=0.7059905246056226;
      pinp[1]=0.05616723541492438;
      winp=0.011303854062073894;
    };
    if(i==2){
      pinp[0]=0.3445413425107939;
      pinp[1]=0.8732597131447893;
      winp=0.01733620243623581;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO10()
{
  int sn=10;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  for(int i=0;i<3;i++){
    if(i==0){
      pinp[0]=0.08854339907163095;
      pinp[1]=0.9921290908756198;
      winp=0.004233324830259320;
    };
    if(i==1){
      pinp[0]=0.2850048087327168;
      pinp[1]=0.9151746833493821;
      winp=0.009890082229875418;
    };
    if(i==2){
      pinp[0]=0.6660539086474957;
      pinp[1]=0.3357743015044285;
      winp=0.01508167827282694;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  int x1[]={0,0,1,1,2,2};
  int x2[]={1,2,0,2,0,1};
  int x3[]={2,1,2,0,1,0};
  for(int i=0;i<1;i++){
    if(i==0){
      pinp[0]=0.05573276265112507;
      pinp[1]=0.5023501980082433;
      pinp[2]=0.8628662339716117;
      winp=0.006230790666852494;
    };
    for(int j=0;j<6;j++){
      mui[index]=pinp[x1[j]];
      eai[index]=pinp[x2[j]];
      xii[index]=pinp[x3[j]];
      omega[index++]=winp;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO12()
{
  int sn=12;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  for(int i=0;i<3;i++){
    if(i==0){
      pinp[0]=0.07668616425146965;
      pinp[1]=0.9941018380552333;
      winp=0.003064844706191271;
    };
    if(i==1){
      pinp[0]=0.2454508067179744;
      pinp[1]=0.9378207733692996;
      winp=0.008399197314589212;
    };
    if(i==2){
      pinp[0]=0.4628296650685621;
      pinp[1]=0.7560273819545462;
      winp=0.01147705941274838;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  int x1[]={0,0,1,1,2,2};
  int x2[]={1,2,0,2,0,1};
  int x3[]={2,1,2,0,1,0};
  for(int i=0;i<2;i++){
    if(i==0){
      pinp[0]=0.01982409694301226;
      pinp[1]=0.4311826481193306;
      pinp[2]=0.9020468552914508;
      winp=0.003262057706317985;
    };
    if(i==1){
      pinp[0]=0.1577019659607366;
      pinp[1]=0.6149489311638969;
      pinp[2]=0.7726369794363323;
      winp=0.006100724910250915;
    };
    for(int j=0;j<6;j++){
      mui[index]=pinp[x1[j]];
      eai[index]=pinp[x2[j]];
      xii[index]=pinp[x3[j]];
      omega[index++]=winp;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO14()
{
  int sn=14;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  real tmp1=0.5773502691896258;
  mui[index]=tmp1;
  eai[index]=tmp1;
  xii[index]=tmp1;
  omega[index++]=0.008080816270852758;

  for(int i=0;i<3;i++){
    if(i==0){
      pinp[0]=0.02805168459462732;
      pinp[1]=0.9992127931440865;
      winp=0.0005358799060520819;
    };
    if(i==1){
      pinp[0]=0.4243600729477484;
      pinp[1]=0.7998981541268634;
      winp=0.007965951316354112;
    };
    if(i==2){
      pinp[1]=0.2283518463067843;
      pinp[0]=0.6884240823388890;
      winp=0.008241164702920527;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  int x1[]={0,0,1,1,2,2};
  int x2[]={1,2,0,2,0,1};
  int x3[]={2,1,2,0,1,0};
  for(int i=0;i<3;i++){
    if(i==0){
      pinp[0]=0.01434849841296581;
      pinp[1]=0.6443413732198754;
      pinp[2]=0.7646033712654021;
      winp=0.001740154916756853;
    };
    if(i==1){
      pinp[0]=0.04513452112438377;
      pinp[1]=0.2015730694667628;
      pinp[2]=0.9784330189995737;
      winp=0.002561193800372772;
    };
    if(i==2){
      pinp[0]=0.1401725432202530;
      pinp[1]=0.3966561951922939;
      pinp[2]=0.9072020287360138;
      winp=0.006813683941731556;
    };
    for(int j=0;j<6;j++){
      mui[index]=pinp[x1[j]];
      eai[index]=pinp[x2[j]];
      xii[index]=pinp[x3[j]];
      omega[index++]=winp;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutEO16()
{
  int sn=16;
  int tmp=(sn+2)*sn/8;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=tmp*(Ndim-1)*4;
  PutSN(totsn);

  int index=0;
  real winp;
  real pinp[3];

  for(int i=0;i<4;i++){
    if(i==0){
      pinp[0]=0.0243277493095127;
      pinp[1]=0.9994079853728747;
      winp=0.0003786929095150918;
    };
    if(i==1){
      pinp[1]=0.3457912005135811;
      pinp[0]=0.6634864149503652;
      winp=0.004516919686491092;
    };
    if(i==2){
      pinp[0]=0.5008009847472008;
      pinp[1]=0.7059722001272202;
      winp=0.004690337133527791;
    };
    if(i==3){
      pinp[1]=0.1017609854155053;
      pinp[0]=0.7034361029429984;
      winp=0.005202608997169839;
    };
    mui[index]=pinp[1];
    eai[index]=pinp[0];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[1];
    xii[index]=pinp[0];
    omega[index++]=winp;
    mui[index]=pinp[0];
    eai[index]=pinp[0];
    xii[index]=pinp[1];
    omega[index++]=winp;
  };

  int x1[]={0,0,1,1,2,2};
  int x2[]={1,2,0,2,0,1};
  int x3[]={2,1,2,0,1,0};
  for(int i=0;i<4;i++){
    if(i==0){
      pinp[0]=0.03767117495065577;
      pinp[1]=0.1643090746848965;
      pinp[2]=0.9856893073144449;
      winp=0.001676516384733566;
    };
    if(i==1){
      pinp[0]=0.02052821406587529;
      pinp[1]=0.4941178958527522;
      pinp[2]=0.8691525167801761;
      winp=0.002061020005761503;
    };
    if(i==2){
      pinp[0]=0.1308815532023498;
      pinp[1]=0.3035009313489606;
      pinp[2]=0.9437993450419713;
      winp=0.003807703743521037;
    };
    if(i==3){
      pinp[0]=0.2544227835824535;
      pinp[1]=0.4710496539585243;
      pinp[2]=0.8446190091986508;
      winp=0.005893813835965320;
    };
    for(int j=0;j<6;j++){
      mui[index]=pinp[x1[j]];
      eai[index]=pinp[x2[j]];
      xii[index]=pinp[x3[j]];
      omega[index++]=winp;
    };
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

void Quadrature::PutAEO8()
{
  // Advanced EO8

  int tmp=31;
  real *mui=new real[tmp];
  real *eai=new real[tmp];
  real *xii=new real[tmp];

  int totsn=31*8;
  PutSN(totsn);

  real omegain[]={
    0.001137877548384876, 0.001137877548384876, 0.001137877548384876,
    0.001600637513808256, 0.001600637513808256, 0.001600637513808256,
    0.001600637513808256, 0.001600637513808256, 0.001600637513808256,
    0.003421036575518821, 0.003421036575518821, 0.003421036575518821,
    0.003883262983893857, 0.003883262983893857, 0.003883262983893857, 
    0.003883262983893857, 0.003883262983893857, 0.003883262983893857, 
    0.004845695722472774, 0.004845695722472774, 0.004845695722472774,
    0.004845695722472774, 0.004845695722472774, 0.004845695722472774,
    0.006502464080380583,
    0.006816100385266683, 0.006816100385266683, 0.006816100385266683,
    0.007464971690352985, 0.007464971690352985, 0.007464971690352985,
  };
  real con[]={
    0.04609696113773013, 0.9978728077003267, 0.01945840472801398,
    0.2760302416227376, 0.9609519635211398, 0.1569168423849870,
    0.9750662588521101, 0.05540059111863282, 0.5589236629456942,
    0.8273663719919212, 0.2026744515387446, 0.3895737576397144,
    0.8984182511792278, 0.5773502691896258, 0.4434446717974471,
    0.7789182538038950, 0.2474556441181923, 0.6851152108200673
  };
  int muid[]={
    0,0,1,
    2,3,2,4,3,4,
    5,5,6,
    7,8,7,9,8,9,
    10,11,10,12,11,12,
    13,
    14,14,15,
    16,17,17
  };
  int eaid[]={
    0,1,0,
    3,2,4,2,4,3,
    5,6,5,
    8,7,9,7,9,8,
    11,10,12,10,12,11,
    13,
    14,15,14,
    17,16,17
  };
  int xiid[]={
    1,0,0,
    4,4,3,3,2,2,
    6,5,5,
    9,9,8,8,7,7,
    12,12,11,11,10,10,
    13,
    15,14,14,
    17,17,16
  };
  for(int i=0;i<31;i++){
    mui[i]=con[muid[i]];
    eai[i]=con[eaid[i]];
    xii[i]=con[xiid[i]];
    omega[i]=omegain[i];
  };

  ArrangeDiscretedPoint(tmp,mui,eai,xii);
  WeightNormalize();
  CalXYZref();
  CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
};

// The following quadrature devide mu[-1,0] into two ranges
/*
void Quadrature::PutDPnDPnTn(int dpn1,int dpn2,int tn,real bmu,string dir)
{
  if(dpn1<=0||dpn2<=0){
    cout<<"Error in PutDPnDPnTn.\n";
    cout<<"You should set dpn1>0 and dpn2>0.\n";
    exit(0);
  };
  if(tn<1){
    cout<<"Error in PutDPnDPnTn.\n";
    cout<<"You should put Tn>0.\n";
    exit(0);
  };
  if(bmu<=0.||bmu>=1.){
    cout<<"Error in PutDPnDPnTn.\n";
    cout<<"You should set 0<bmu<1.\n";
    exit(0);
  };

  real bmu2=1.-bmu;
  real *mu1=new real[dpn1];
  real *mu2=new real[dpn2];
  RootOfLegendre(dpn1,mu1);
  RootOfLegendre(dpn2,mu2);
  real *w1=new real[dpn1];
  real *w2=new real[dpn2];
  for(int i=0;i<dpn1;i++){
    w1[i]=WeightOfGaussian(dpn1,mu1[i])*bmu2;
  };
  for(int i=0;i<dpn2;i++){
    w2[i]=WeightOfGaussian(dpn2,mu2[i])*bmu;
  };
  
  real *ep=new real[dpn1+dpn2];
  real *w=new real[dpn1+dpn2];

  int id=0;
  for(int i=0;i<dpn1;i++){
    real temp=mu1[dpn1-1-i];
    ep[id]=bmu+(1.+temp)*bmu2/2.;
    w[id]=w1[dpn1-1-i];
    id++;
  };
  for(int i=0;i<dpn2;i++){
    real temp=mu2[dpn2-1-i];
    ep[id]=(1.+temp)*bmu/2.;
    w[id]=w2[dpn2-1-i];
    id++;
  };
  DPnTn(dpn1+dpn2,tn,ep,w,dir);

  delete [] mu1;
  delete [] mu2;
  delete [] w1;
  delete [] w2;
  delete [] ep;
  delete [] w;
};
*/

