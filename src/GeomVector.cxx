#include <iostream>
#include <fstream>
#include <cstdio>
#include "GeomVector.h"

using namespace std;

void GeomVector::Initialize(real a,real b)
{
  x=a; y=b;
}

GeomVector GeomVector::operator+(const GeomVector &sec)
{
  return GeomVector(x+sec.x,y+sec.y);
}

GeomVector GeomVector::operator-(const GeomVector &sec)
{
  return GeomVector(x-sec.x,y-sec.y);
}

GeomVector GeomVector::operator*(real a)
{
  return GeomVector(x*a,y*a);
}

real GeomVector::CalcDistance(GeomVector sec)
{
  real a=(x-sec.x)*(x-sec.x)+(y-sec.y)*(y-sec.y);
  a=sqrt(a);
  return a;
}

void GeomVector::show_self()
{
  cout.setf(ios::showpoint);
  cout.precision(5);
  cout << "** GeomVector **\n";
  cout << "x:"<<x<<"  y:"<<y<<"\n" ;
}

void GeomVector::NormalizedToUnity()
{
  real enorm=sqrt(x*x+y*y);
  if(enorm==0.){
    cout<<"# Error in GeomVector::NormalizedToUnity().\n";
    cout<<"# The norm of the vector is zero.\n";
    exit(0);
  };

  real ar=1./enorm;
  x*=ar;  y*=ar;
}

int GeomVector::SameOrNot(GeomVector sec,real Judge)
{
  int ret=1;
  if(fabs(x-sec.getx())>Judge)return 0;
  if(fabs(y-sec.gety())>Judge)return 0;
  return ret;
}

bool GeomVector::BeIdentical(GeomVector sec,real Judge)
{
  if(fabs(x)<0.00001){
    if(fabs(sec.getx()-x)>Judge)return false;
  }else{
    if(fabs((sec.getx()-x)/x)>Judge)return false;
  };

  if(fabs(y)<0.00001){
    if(fabs(sec.gety()-y)>Judge)return false;
  }else{
    if(fabs((sec.gety()-y)/y)>Judge)return false;
  };

  return true;
}

GeomVector GeomVector::OperatingElementaryReflector(GeomVector u)
{
  // Elementary reflector about u^{perp}
  //
  //   R = I - 2 uu^T/(u^T u)
  //

  u.NormalizedToUnity();
  real tmp=u.getx()*x+u.gety()*y;
  GeomVector ret(x-2*tmp*u.getx(),y-2*tmp*u.gety());
  return ret;
};

GeomVector GeomVector::OperatingRotator(real rot)
{
  // rot = 1 : 2PI

  real ang=PI2*rot;
  real c=cos(ang);
  real s=sin(ang);
  real xnew=c*x-s*y;
  real ynew=s*x+c*y;
  GeomVector ret(xnew,ynew);
  //cout<<"rot : "<<rot<<" "<<x<<" "<<y<<" "<<xnew<<" "<<ynew<<" "<<"\n";
  return ret;
};

void GeomVector::Shift(real dx,real dy)
{
  x+=dx;
  y+=dy;
};
