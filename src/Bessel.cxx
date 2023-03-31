#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Bessel.h"

using namespace std;

real bessi0(real x)
{
  real p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y;
  real ax,ret;

  p1=1.0;
  p2=3.5156229;
  p3=3.0899424;
  p4=1.2067492;
  p5=0.2659732;
  p6=0.360768e-1;
  p7=0.45813e-2;

  q1=0.39894228;
  q2=0.1328592e-1;
  q3=0.225319e-2;
  q4=-0.157565e-2;
  q5=0.916281e-2;
  q6=-0.2057706e-1;
  q7=0.2635537e-1;
  q8=-0.1647633e-1;
  q9=0.392377e-2;

  if (x>-3.75 && x<3.75){
    y=(x/3.75)*(x/3.75);
    ret = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
  }
  else{
    if(x>0.0){ax=x;}
    else{ax=-x;};
    y=3.75/ax;
    ret = (exp(ax)/sqrt(ax))*(q1+y*(q2+y
    * (q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
  };
  return ret;
}

real bessk0(real x)
{
  real p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y;
  real ret;

  p1=-0.57721566;
  p2=0.42278420;
  p3=0.23069756e0;
  p4=0.3488590e-1;
  p5=0.262698e-2;
  p6=0.10750e-3;
  p7=0.74e-5;

  q1=1.25331414;
  q2=-0.7832358e-1;
  q3=0.2189568e-1;
  q4=-0.1062446e-1;
  q5=0.587872e-2;
  q6=-0.251540e-2;
  q7=0.53208e-3;

  if (x<=2.0){
    y=x*x/4.0;
    ret = (-log(x/2.0)*bessi0(x))
      + (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
  }
  else{
    y=2.0/x;
     ret= (exp(-x)/sqrt(x))*(q1+y*(q2
      +y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
  };

  return ret;
}

real bessi1(real x)
{
  real p1,p2,p3,p4,p5,p6,p7;
  real q1,q2,q3,q4,q5,q6,q7,q8,q9;
  real y,ax;
  real ret;

  p1=0.5e0;
  p2=0.87890594e0;
  p3=0.51498869e0;
  p4=0.15084934;
  p5=0.2658733e-1;
  p6=0.301532e-2;
  p7=0.32411e-3;
  q1=0.39894228;
  q2=-0.3988024e-1;
  q3=-0.362018e-2;
  q4=0.163801e-2;
  q5=-0.1031555e-1;
  q6=0.2282967e-1;
  q7=-0.2895312e-1;
  q8=0.1787654e-1;
  q9=-0.420059e-2;

  if(x<3.75 && x>-3.75){
    y=(x/3.75)*(x/3.75);
    ret = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
  }
  else{
    ax=x;
    if(x<0.0)ax=-x;
    y=3.75/ax;
    ret = (exp(ax)/sqrt(ax))
         * (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y
		 * (q7+y*(q8+y*q9))))))));
    if(x<0.0) ret=-ret;
  };

  return ret;
}

real bessk1(real x)
{
  real p1,p2,p3,p4,p5,p6,p7;
  real q1,q2,q3,q4,q5,q6,q7;
  real y;
  real ret;

  p1=1.0;
  p2=0.15443144;
  p3=-0.67278579;
  p4=-0.18156897;
  p5=-0.1919402e-1;
  p6=-0.110404e-2;
  p7=-0.4686e-4;
  q1=1.25331414e0;
  q2=0.23498619e0;
  q3=-0.3655620e-1;
  q4=0.1504268e-1;
  q5=-0.780353e-2;
  q6=0.325614e-2;
  q7=-0.68245e-3;

  if(x<2.0){
    y=x*x/4.0;
    ret= (log(x/2.0)*bessi1(x)) + (1.0/x)
      *(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));}
  else{
    y=2.0/x;
    ret= (exp(-x)/sqrt(x))
      * (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));};

  return ret;
}

real gammln(real xx)
{
  real cof[]={76.18009173, -86.50532033,
               24.01409822, -1.231739516,
               0.120858003e-2, -0.536382e-5};
  real stp=2.50662827465;
  real x,tmp,ser;

  x=xx-1.0;
  tmp=x+5.5;
  tmp=(x+0.5)*log(tmp)-tmp;
  ser=1.0;
  for (int i=0;i<6;i++){
    x+=1.0;
    ser+=cof[i]/x;
  };
  return tmp+log(stp*ser);
}

real gamma(real x)
{
  real ret;
  if(x>0.0){
    if(x>=1.0){ret=exp(gammln(x));}
    else{ret=exp(gammln(x+1.0)-log(x));};
  }else{
    cout<<"Not coded gamma(negative value) in Bessel.\n";
    exit(0);
  };
  return ret;
}

real struveh(int iv,real z)
{
  real stv0,stv,cof1,cof2,dstv,ret,r,r2;
  real r3,r4;

  stv0=pow(z*0.5,iv+1);
  stv=0.0;
  for(int k=0;k<100;k++){
    r3=real(k)+1.5;
    r4=real(k+iv)+1.5;
    cof1=gamma(r3);
    cof2=gamma(r4);
    dstv=pow(-1.0,k)*pow(z*0.5,k*2)/cof1/cof2;
    r=stv*1e-14;
    if(r<0.0)r=-r;
    r2=dstv;
    if(r2<0.0)r2=-r2;
    if(r2<r)k=100;
    stv=stv+dstv;
  };
  ret=stv0*stv;
  return ret;
};

real struvel(int iv,real z)
{
  real stv0,stv,cof1,cof2,dstv,ret;
  real r,r2,r3,r4;

  stv0=pow(z*0.5,iv+1);
  stv=0.0;
  for (int i=0;i<100;i++){
    r3=real(i)+1.5;
    r4=real(i+iv)+1.5;
    cof1=gamma(r3);
    cof2=gamma(r4);
    dstv=pow(z*0.5,2*i)/cof1/cof2;
    r=dstv;
    if(r<0.0)r=-r;
    r2=stv;
    if(r2<0.0)r2=-r2;
    r2*=1e-14;
    if(r<r2)i=100;
    stv=stv+dstv;
  };
  ret =stv0*stv;
  return ret;
};

real ibessy0(real x)
{
  real ret;

#ifdef __linux__
  //real y0t=y0f(x); // "y0f" returns a float value.
  real y0t=y0(x);
  //ret=x*y0t+PI*0.5*x*(struveh(0,x)*y1f(x)-struveh(1,x)*y0t);
 // "y1f" returns a float value.
  ret=x*y0t+PI*0.5*x*(struveh(0,x)*y1(x)-struveh(1,x)*y0t);
#else
#ifdef WIN32
  real y0=_y0(x);
  real y1=_y1(x);
  ret=x*y0+PI*0.5*x*(struveh(0,x)*y1-struveh(1,x)*y0);
#else
#ifdef __APPLE__
  real y0x=y0(x);
  real y1x=y1(x);
  ret=x*y0x+PI*0.5*x*(struveh(0,x)*y1x-struveh(1,x)*y0x);
#else
  cerr << "ibessy0: Not Implemented Error" << endl;
  return(1);
#endif
#endif
#endif
  return ret;
};

real ibessj0(real x)
{
  real ret;

#ifdef __linux__
  real j0=j0f(x);
  ret=x*j0+PI*0.5*x*(struveh(0,x)*j1f(x)-struveh(1,x)*j0);
#else
#ifdef WIN32
  real j0=_j0(x);
  real j1=_j1(x);
  ret=x*j0+PI*0.5*x*(struveh(0,x)*j1-struveh(1,x)*j0);
#else
#ifdef __APPLE__
  real j0x=j0(x);
  real j1x=j1(x);
  ret=x*j0x+PI*0.5*x*(struveh(0,x)*j1x-struveh(1,x)*j0x);
#else
  cerr << "ibessj0: Not Implemented Error" << endl;
  return(1);
#endif
#endif
#endif
  return ret;
};

real ibessk0(real x)
{
  real ret;
  ret=x*bessk0(x)+PI*0.5*x*
    (struvel(0,x)*bessk1(x)+struvel(1,x)*bessk0(x));
  return ret;
};
