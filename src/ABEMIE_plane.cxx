#include <iostream>
#include <cstdlib>
#include "ABEMIE_plane.h"

using namespace std;

void plane::set_imax(int i)
{
  imax=i;
  for(int j=0;j<pnt;j++){
    flx[j].put_imax(i);
    cur[j].put_imax(i);
    dfdk[j].put_imax(i);
  };
}

real plane::CalcDistance2(GeomVector p)
{
  real ret=CalcDistance(p);
  real d=a*p.getx()+b*p.gety()+c;
  real e=a*(p1.getx()+pver.getx())+b*(p1.gety()+pver.gety())+c;
  if(d*e>0.0)ret=-ret;
  return ret;
}

void plane::cal_ghij(real b2, real bdk, GeomVector &sp,
                     real *gr, real *hr, real *dgr, real *dhr)
{
  real l=get_long();
  int igauss=24;
  if(l<=25)igauss=16;
  real *ww=new real[igauss];
  real *pp=new real[igauss];

  if(l>25){
    for(int i=0;i<igauss;i++){
      ww[i]=ww24[i];
      pp[i]=pp24[i];
    };
  }
  else{
    for(int i=0;i<igauss;i++){
      ww[i]=ww16[i];
      pp[i]=pp16[i];
    };
  };

  real *po=new real[igauss*pnt];

  switch(pnt){
  case 1:
    for(int i=0;i<igauss;i++){
      po[i]=ww[i];
    };
    break;
  case 2:
    for(int i=0;i<igauss;i++){
      real co=0.5*ww[i]*htp2_inv;
      po[i*2]=(htp2-pp[i])*co;
      po[i*2+1]=(htp2+pp[i])*co;
    };
    break;
  case 3:
    for(int i=0;i<igauss;i++){
      real mh=pp[i]-htp;
      real ph=pp[i]+htp;
      real co=ww[i]*htp_inv*htp_inv;
      po[i*3]  = 0.5*pp[i]*mh*co;
      po[i*3+1]=-mh*ph*co;
      po[i*3+2]= 0.5*pp[i]*ph*co;
    };
    break;
  case 4:
    real a2b2=ht2*ht2-htp4*htp4;
    for(int i=0;i<igauss;i++){
      real mh =pp[i]-htp4;
      real mh2=pp[i]-ht2;
      real ph =pp[i]+htp4;
      real ph2=pp[i]+ht2;
      real co =0.5/a2b2*ww[i]*mh*mh2*ph*ph2;
      po[i*4]=   co*htp4_inv/ph;
      po[i*4+1]=-co*ht2_inv/ph2;
      po[i*4+2]= co*ht2_inv/mh2;
      po[i*4+3]=-co*htp4_inv/mh;
    };
  }

  real b;
  if(b2<0.0){b=sqrt(-b2);}else{b=sqrt(b2);};

  for(int i=0;i<pnt;i++){
    gr[i]=0.0;  hr[i]=0.0;  dgr[i]=0.0; dhr[i]=0.0;
  };
  
  GeomVector pon;
  real rri,dst;
  if(b2<0.0){
    for(int i=0;i<igauss;i++){
      pon=p1*((1.0-pp[i])*0.5)+p2*((1.0+pp[i])*0.5);
      rri=sp.CalcDistance(pon);
      dst=CalcDistance2(sp);
      real rrb=rri*b;
      real bk0=bessk0(rrb);
      real bk1=bessk1(rrb);
      real dst_rri_inv=dst/rri;
      for(int j=0;j<pnt;j++){
        real poi=po[i*pnt+j];
	real tmp=poi*bk0;
	real tmp2=poi*bk1;
        gr[j] +=tmp;
        hr[j] -=tmp2*dst_rri_inv;
        dgr[j]-=tmp2*rri;
        dhr[j]+=tmp*dst;
      };
    };
  }
  else{
    for(int i=0;i<igauss;i++){
      pon=p1*((1.0-pp[i])*0.5)+p2*((1.0+pp[i])*0.5);
      rri=sp.CalcDistance(pon);
      dst=CalcDistance2(sp);
      real rrb=rri*b;
      //real by0=y0f(rrb); // y0f is for "float"
      //real by1=y1f(rrb); // y1f is for "float"
      real by0=y0(rrb);
      real by1=y1(rrb);
      real dst_rri_inv=dst/rri;
      for(int j=0;j<pnt;j++){
	real poi=po[i*pnt+j];
	real tmp=poi*by0;
	real tmp2=poi*by1;
        gr[j] -=tmp;
        hr[j] +=tmp2*dst_rri_inv;
        dgr[j]+=tmp2*rri;
        dhr[j]+=tmp*dst;
      };
    };
  };

  real gm1=p1.CalcDistance(p2)*0.5;
  if(b2<0.0){
    real pi2g=0.5*INV_PI*gm1;
    for(int i=0;i<pnt;i++){
      gr[i] *=pi2g;
      hr[i] *=pi2g*b;
      dgr[i]*=bdk*pi2g;
      dhr[i]*=bdk*pi2g*b;
    };
  }
  else{
    real g4=gm1*0.25;
    for(int i=0;i<pnt;i++){
      gr[i] *=g4;
      hr[i] *=g4*b;
      dgr[i]*=bdk*g4;
      dhr[i]*=bdk*g4*b;
    };
  };

  delete []ww;
  delete []pp;
  delete []po;
}

void plane::cal_ghii(real b2,real bdk, GeomVector &sp,
                     real *gr, real *hr, real *dgr, real *dhr)
{
  real *aa=new real[4*pnt];
  real *ap=new real[4*pnt*2];
  int pnt4=pnt*4;

  switch(pnt){
  case 1:
    aa[0]=0; aa[1]=0; aa[2]=0; aa[3]=1; 
    break;
  case 2:
    aa[0]=0; aa[1]=0; aa[2]=-0.5/htp2; aa[3]=0.5;
    aa[4]=0; aa[5]=0; aa[6]= 0.5/htp2; aa[7]=0.5;
    break;
  case 3:
    aa[0]=0; aa[1]=0.5*htp_inv*htp_inv; aa[2]=-0.5*htp_inv; aa[3]=0;
    aa[4]=0; aa[5]=-1*htp_inv*htp_inv;  aa[6]=0;        aa[7]=1;
    aa[8]=0; aa[9]=0.5*htp_inv*htp_inv; aa[10]=0.5*htp_inv;  aa[11]=0;
    break;
  case 4:
    real ab=ht2*ht2-htp4*htp4;
    real fow=0.5/htp4/ab;
    aa[0]=fow;  aa[1]=-fow*htp4; aa[2]=-fow*ht2*ht2; aa[3]=fow*ht2*ht2*htp4;
    fow=-0.5/ht2/ab;
    aa[4]=fow;  aa[5]=-fow*ht2; aa[6]=-fow*htp4*htp4; aa[7]=fow*ht2*htp4*htp4;
    fow=0.5/ht2/ab;
    aa[8]=fow;  aa[9]=fow*ht2; aa[10]=-fow*htp4*htp4; aa[11]=-fow*ht2*htp4*htp4;
    fow=-0.5/htp4/ab;
    aa[12]=fow; aa[13]=fow*htp4; aa[14]=-fow*ht2*ht2; aa[15]=-fow*ht2*ht2*htp4;
    break;
  }

  for(int i=0;i<pnt;i++){
    gr[i]=0.0;  hr[i]=0.0; dgr[i]=0.0; dhr[i]=0.0;
  };

  real b;
  if(b2<0.0){b=sqrt(-b2);}else{b=sqrt(b2);}

  real gm1=p1.CalcDistance(p2);
  real xsi=get_xsi(sp);
  real bg=b*gm1;
  real bg2=b*b*gm1*gm1;
  real tmp1=1./(bg2*bg);
  real xsi2=xsi*xsi;
  real bg_inv=1./bg;
  real bg2_inv=1./bg2;

  for(int i=0;i<pnt;i++){
    int i4=i*4;
    int i4p=i4+pnt4;
    ap[i4]   =8*aa[i4]*tmp1;
    ap[i4p]  =-ap[i4];
    ap[i4+1] =(12*aa[i4]*xsi+4*aa[i4+1])*bg2_inv;
    ap[i4p+1]=ap[i*4+1];
    ap[i4+2] =(6*aa[i4]*xsi2+4*aa[i4+1]*xsi+2*aa[i4+2])*bg_inv;
    ap[i4p+2]=-ap[i*4+2];
    ap[i4+3] =aa[i4]*xsi2*xsi+aa[i4+1]*xsi2+aa[i4+2]*xsi+aa[i4+3];
    ap[i4p+3]=ap[i4+3];
  };

  real *x=new real[2];
  x[0]=0.5*b*(1.0-xsi)*gm1;
  x[1]=0.5*b*(1.0+xsi)*gm1;
  real b_inv=1./b;
  real b2_inv=1./(b*b);
  if(b2<0.0){
    for(int i=0;i<pnt;i++){
      for(int j=0;j<2;j++){
	real x1=x[j];
	real x2=x1*x1;
	real x3=x2*x1;
	int indx=j*pnt4+i*4;
	real k0=bessk0(x1);
	real k1=bessk1(x1);
	real ik0=ibessk0(x1);
	real i1k1=-x1*k0+ik0;
	gr[i]+= ap[indx+3]*ik0;
	dgr[i]+=ap[indx+3]*i1k1;
	if(pnt>1){
	  real i1k0=-x1*k1;
          real i2k1=-x2*k0+2*i1k0;
          gr[i]+= ap[indx+2]*i1k0;
	  dgr[i]+=ap[indx+2]*i2k1;
  	  if(pnt>2){
	    real i2k0=-x2*k1+i1k1;
	    real i3k1=-x3*k0+3*i2k0;
            gr[i]+= ap[indx+1]*i2k0;
	    dgr[i]+=ap[indx+1]*i3k1;
  	    if(pnt>3){
	      real i2k1=-x2*k0+2*i1k0;
	      real i3k0=-x3*k1+2*i2k1;
	      real i4k1=-x3*x1*k0+4*i3k0;
              gr[i]+= ap[indx]*i3k0;
	      dgr[i]+=ap[indx]*i4k1;
	    };
	  };
	};
      };
      gr[i] *=0.5*INV_PI*b_inv;
      dgr[i]*=-0.5*INV_PI*bdk*b2_inv;
    };
  }
  else{
    for(int i=0;i<pnt;i++){
      for(int j=0;j<2;j++){
	int indx=j*pnt4+i*4;
	real x1=x[j];
	real x2=x1*x1;
	real x3=x2*x1;
	//real y0=y0f(x1);
	//real y1=y1f(x1);
	real by0=y0(x1);
	real by1=y1(x1);
	real iy0=ibessy0(x1);
	real i1y1=-x1*by0+iy0;
	gr[i]+=ap[indx+3]*iy0;
	dgr[i]+=ap[indx+3]*i1y1;
	if(pnt>1){
	  real i1y0=x1*by1;
	  real i2y1=-x2*by0+2*i1y0;
	  gr[i]+=ap[indx+2]*i1y0;
	  dgr[i]+=ap[indx+2]*i2y1;
  	  if(pnt>2){
	    real i2y0=x2*by1-i1y1;
	    real i3y1=-x3*by0+3*i2y0;
	    gr[i]+=ap[indx+1]*i2y0;
            dgr[i]+=ap[indx+1]*i3y1;
  	    if(pnt>3){
	      real i3y0=x3*by1-2*i2y1;
	      real i4y1=-x3*x1*by0+4*i3y0;
	      gr[i]+=ap[indx]*i3y0;
              dgr[i]+=ap[indx]*i4y1;
	    };
	  };
	};
      };
      gr[i] *=-0.25*b_inv;
      dgr[i]*=0.25*bdk*b2_inv;
    };
  };
  delete []x;
  delete []aa;
  delete []ap;
}

int plane::opposite_p(int i)
{
  if(kind==5)return pnt-i-1;
  if(kind==2)return i;
  cout<<"Error in ABEMIE_plane::opposite_p.\n";
  cout<<" kind is "<<kind<<"\n";
  show_self();
  exit(0);
}

GeomVector plane::get_pos(int ip)
{
  real xsi=0.;

  switch(pnt){
    case 1:
      xsi=0;
      break;
    case 2:
      if(ip==0)xsi=-htp2;
      if(ip==1)xsi=htp2;
      break;
    case 3:
      if(ip==0)xsi=-htp;
      if(ip==1)xsi=0;
      if(ip==2)xsi=htp;
      break;
    case 4:
      if(ip==0)xsi=-htp4;
      if(ip==1)xsi=-ht2;
      if(ip==2)xsi=+ht2;
      if(ip==3)xsi=+htp4;
      break;
  }

  GeomVector sp=GetVecXsi(xsi);
  return sp;
}

void plane::put_data1da(int i,char *ss,GroupData1D inp)
{
  if(strcmp(ss,"flx")==0)     flx[i].copy(inp);
  if(strcmp(ss,"cur")==0)     cur[i].copy(inp);
  if(strcmp(ss,"dfdk")==0)    dfdk[i].copy(inp);
}

void plane::change_pnt(int ii)
{
  real fl[imax],cu[imax];
  for(int i=0;i<imax;i++){
    fl[i]=flx[0].get_dat(i);
    cu[i]=cur[0].get_dat(i);
  };

  put_pnt(ii);
  set_imax(imax);
  for(int i=0;i<imax;i++){
    for(int j=0;j<pnt;j++){
      flx[j].put_data(i,fl[i]/pnt);
      cur[j].put_data(i,cu[i]/pnt);
    };
  };
}

void plane::put_pnt(int i)
{
  pnt=i;
  flx.resize(pnt);
  cur.resize(pnt);
  dfdk.resize(pnt);
};
