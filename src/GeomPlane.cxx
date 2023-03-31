#include <iostream>
#include <fstream>
#include <cstdio>
#include "GeomPlane.h"
#include "Numeric.h"

using namespace std;

void GeomPlane::PutVector(GeomVector i,GeomVector j)
{
  Exist=true;
  p1=i;
  p2=j;
  if(p1.BeIdentical(p2)){
    //cout<<"In this GeomPlane instance, p1 is same as p2.\n";
  }
  put_line();
}

void GeomPlane::PutVector(real x1,real y1,real x2,real y2)
{
  p1.put_data(x1,y1);
  p2.put_data(x2,y2);
  PutVector(p1,p2);
}

void GeomPlane::put_line()
{
  real xx=fabs(p2.getx()-p1.getx());
  real yy=fabs(p2.gety()-p1.gety());

  if(xx<0.00001){a=1.0;
    b=0.0;
    c=-p2.getx();
    a1=0.0;
    b1=1.0;
  }
  else if(yy<0.00001){a=0.0;
    b=1.0;
    c=-p2.gety();
    a1=1.0;
    b1=0.0;
  }
  else{
    a=(p1.gety()-p2.gety())/(p1.getx()-p2.getx());
    b=-1.0;
    c=p1.gety()-a*p1.getx();
    a1=1.0;
    b1=a;
  };

  GeomVector r;
  r=(p1*1.235+p2*2.765)*0.25;
  c1=-a1*r.getx()-b1*r.gety();

  pver.put_data(b1,-a1);
  pver.unity();
}

void GeomPlane::show_self()
{
  cout<<"*** Information of GeomPlane ***\n";
  if(Exist){
    cout <<"*p1*\n";
    p1.show_self();
    cout <<"*p2*\n";
    p2.show_self();
    cout<<"The line:"<<a<<"x+"<<b<<"y+"<<c<<"=0\n";
    cout<<"The vertical line:"<<a1<<"x+"<<b1<<"y+"<<c1<<"=0\n";
  }
  else{
    cout<<"This plane is not realized.\n";
  };
  cout<<"********************************\n";
}

real GeomPlane::CalcDistance(GeomVector &p)
{
  real tmp=1./sqrt(a*a+b*b);
  real ret=fabs(a*p.getx()+b*p.gety()+c)*tmp;
  return ret;
}

real GeomPlane::CalcDistance(real x,real y)
{
  real tmp=1./sqrt(a*a+b*b);
  real ret=fabs(a*x+b*y+c)*tmp;
  return ret;
}

bool GeomPlane::BeIdentical(GeomPlane sec,real judge)
{
  if(p1.BeIdentical(sec.get_p1(),judge)&&p2.BeIdentical(sec.get_p2(),judge))return true;
  if(p1.BeIdentical(sec.get_p2(),judge)&&p2.BeIdentical(sec.get_p1(),judge))return true;
  return false;
};

bool GeomPlane::BePartial(GeomPlane sec,real judge)
{
  if(sec.VectorOnPlane(p1,judge)&&sec.VectorOnPlane(p2,judge))return true;
  return false;
};

bool GeomPlane::BeCrossing(GeomPlane sec,real judge)
{
  if(a*sec.get_b()-b*sec.get_a()!=0){
    GeomVector pp;
    pp=CrossPlane(sec);
    if(VectorOnPlane(pp,1e-6)&&sec.VectorOnPlane(pp,1e-6))return true;
  };
  return false;
}

GeomVector GeomPlane::CrossPlane(GeomPlane &sec)
{
  /*
  real am[4];
  real bm[2];
  */
  vector<real> am(4);
  vector<real> bm(2);
  am[0]=a;
  am[1]=b;
  am[2]=sec.get_a();
  am[3]=sec.get_b();
  bm[0]=-c;
  bm[1]=-sec.get_c();
  gauss22(am,bm); 
  GeomVector pp(bm[0],bm[1]);
  return pp;
};

bool GeomPlane::VectorOnPlane(GeomVector p, real judge)
{
  int ret=0;
  if(a!=0){
    real tmp=1./a;
    real xx=-b*tmp*p.gety()-c*tmp;
    if(fabs(xx-p.getx())<judge)ret=1;//On line made by plane
  }
  if(a==0){
    real tmp=1./b;
    real yy=-c*tmp;
    if(fabs(yy-p.gety())<judge)ret=1;//On line made by plane
  }

  if(ret==1){
    real d1=p.CalcDistance(p1);
    real d2=p.CalcDistance(p2);
    real lo=p1.CalcDistance(p2);
    //if(fabs(d1-lo)<judge||fabs(d2-lo)<judge)return true; //On plane
    if(fabs(d1-lo)<judge&&d2<judge)return true; //On plane
    if(fabs(d2-lo)<judge&&d1<judge)return true; //On plane
    if(d1<=lo&&d2<=lo)return true; //On plane
  }

  return false;
}

real GeomPlane::get_xsi(GeomVector r)
{
  real x1=p1.getx();
  real x2=p2.getx();
  real x=r.getx();

  real a=fabs(x2-x1);
  if(a<0.0001){
    x1=p1.gety();
    x2=p2.gety();
    x=r.gety();
  };

  real xsi=2.0*(x-x1)/(x2-x1)-1.0;
  return xsi;
}

GeomVector GeomPlane::GetVecXsi(real xsi)
{
  GeomVector ret;
  ret=p1*(1-xsi)*0.5+p2*(1+xsi)*0.5;
  return ret;
};

GeomPlane GeomPlane::OperatingElementaryReflector(GeomVector u)
{
  GeomPlane ret;
  GeomVector v1=p1.OperatingElementaryReflector(u);
  GeomVector v2=p2.OperatingElementaryReflector(u);
  ret.PutVector(v1,v2);
  return ret;
};

GeomPlane GeomPlane::OperatingRotator(real rot)
{
  GeomPlane ret;
  GeomVector v1=p1.OperatingRotator(rot);
  GeomVector v2=p2.OperatingRotator(rot);
  ret.PutVector(v1,v2);
  return ret;
};

void GeomPlane::Shift(real dx, real dy)
{
  p1.Shift(dx,dy);
  p2.Shift(dx,dy);
  pver.Shift(dx,dy);
  put_line();
};
