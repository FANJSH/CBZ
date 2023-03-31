#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "Numeric.h"

using namespace std;

inline void ChangeValue(real &a,real &b)
{
  real c=a;
  a=b;
  b=c;
}

inline void ChangeValue(int &a,int &b)
{
  int c=a;
  a=b;
  b=c;
}

real kaijo(int i1, int i2)
{
  real ret=1;
  for(int j=i1;j>=i2;j--){
    ret*=j;
  };
  return ret;
};

void gauss(real *a,real *b,int n,int m)
{
  cout<<"!!!\n";
  cout<<"!!!\n";
  cout<<"!!!\n";
  cout<<"!!! More efficient algorithm has been implemented into `gauss' with vector members\n";
  cout<<"!!!\n";

  real p;
  real det=1.0;

  for(int id=0;id<n-1;id++){
    det*=a[id*n+id];
    for(int i=id+1;i<n;i++){
      p=a[i*n+id]/a[id*n+id];
      for (int j=id;j<n;j++){
	a[i*n+j]-=p*a[id*n+j];
      };
      for (int j=0;j<m;j++){
	b[i*m+j]-=p*b[id*m+j];
      };
    };
  };
  det*=a[n*n-1];

  for(int j=0;j<m;j++){
    b[(n-1)*m+j]/=a[n*n-1];
  };

  for (int ii=2;ii<n+1;ii++){
    int i=n+1-ii-1;
    for (int j=0;j<m;j++){
      for (int k=i+1;k<n;k++){
        b[i*m+j]-=a[i*n+k]*b[k*m+j];
      };
      b[i*m+j]/=a[i*n+i];
    };
  };
}

void gauss(vector<real> &a,vector<real> &b,int n,int m)
{
  for(int id=0;id<n-1;id++){
    int tmp=id*n+id;
    //if(fabs(a[tmp])<1e-15)cout<<"# Pivot is too small : "<<id<<"\n";
    real a_inv=1./a[tmp];
    for(int i=id+1;i<n;i++){
      real val=a[i*n+id];
      if(val!=0.){
        real p=val*a_inv;
        int tmp2=i*n+id;
        int tmp3=id*n+id+1;
        a[tmp2++]=0.;
        for (int j=id+1;j<n;j++){
	  real vtmp=a[tmp3++];
	  if(vtmp!=0.)a[tmp2]-=p*vtmp;
	  tmp2++;
        };
        tmp2=i*m;
        tmp3=id*m;
        for (int j=0;j<m;j++){
	  real vtmp=b[tmp3++];
	  if(vtmp!=0.)b[tmp2]-=p*vtmp;
	  tmp2++;
        };
      };
    };
  };

  int tmp=(n-1)*m;
  real a_inv=1./a[n*n-1];
  for(int j=0;j<m;j++){
    b[tmp++]*=a_inv;
  };

  for (int ii=2;ii<n+1;ii++){
    int i=n-ii;
    for (int k=i+1;k<n;k++){
      real piv=a[i*n+k];
      if(piv!=0.){
        int tmp=i*m;
        int tmp2=k*m;
        for (int j=0;j<m;j++){
          real vtmp=b[tmp2++];
	  if(vtmp!=0.)b[tmp]-=piv*vtmp;
	  tmp++;
        };
      };
    };

    real a_inv=1./a[i*n+i];
    int tmp=i*m;
    for(int j=0;j<m;j++){
      b[tmp++]*=a_inv;
    };

  };
}

void gauss_inverse(vector<real> &a,int n)
{
  int n2=n*2;

  vector<real> b(n*n2,0.);
  int ii=0;
  int ia=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      b[ii++]=a[ia++];
    };
    ii+=i;
    b[ii]=1.;
    ii+=n-i;
  };

  for(int id=0;id<n-1;id++){
    int tmp=id*n2+id;

    /*
    if(fabs(b[tmp])<1e-10){
      cout<<"# Pivot is too small : "<<id<<"\n";
      for(int ii=id+1;ii<n-1;ii++){
      real ttt=fabs(b[id*n2+ii]);
    	if(ttt>1e-10){
    	  cout<<"#    ... you can change to the column : "<<ii<<"\n";
    	  ii=n-1;
    	};
    };
    };
    */

    real a_inv=1./b[tmp];
    for(int i=id+1;i<n;i++){
      int tmp2=i*n2+id;
      real val=b[tmp2];
      if(val!=0.){
        real p=val*a_inv;
        int tmp3=id*n2+id+1;
        b[tmp2++]=0.;
        for (int j=id+1;j<n+id+1;j++){
	  real tmp=b[tmp3++];
	  if(tmp!=0.)b[tmp2]-=p*tmp;
	  tmp2++;
        };
      };
    };
  };

  int tmp=n*n2-(n+1);
  real a_inv=1./b[tmp];
  tmp++;
  for(int j=0;j<n;j++){
    b[tmp++]*=a_inv;
  };

  for (int ii=2;ii<n+1;ii++){
    int i=n-ii;
    for (int k=i+1;k<n;k++){
      real piv=b[i*n2+k];
      if(piv!=0.){
        int tmp=i*n2+n;
        int tmp2=k*n2+n;
        for(int j=0;j<n;j++){
	  real vtmp=b[tmp2++];
          if(vtmp!=0.)b[tmp]-=piv*vtmp;
	  tmp++;
        };
      };
    };

    real a_inv=1./b[i*n2+i];
    int tmp=i*n2+n;
    for(int j=0;j<n;j++){
      b[tmp++]*=a_inv;
    };

  };

  int ii1=0;
  int ii2=0;
  for(int i=0;i<n;i++){
    ii2+=n;
    for(int j=0;j<n;j++){
      a[ii1++]=b[ii2++];
    };
  };


}

void gauss_pivot_check(vector<real> &a,vector<int> &pos, int n, real eps)
{
  pos.resize(n);
  for(int i=0;i<n;i++){
    pos[i]=i;
  };

  for(int id=0;id<n-1;id++){

    int tmp=id*n+id;

    if(fabs(a[tmp])<eps){
      for(int ii=id+1;ii<n;ii++){
	real ttt=fabs(a[id*n+ii]);
	if(ttt>eps){
	  /*
          pos1.push_back(id);
          pos2.push_back(ii);
	  */
	  ChangeValue(pos[id],pos[ii]);
	  //cout<<id<<" "<<ii<<"\n";
	  for(int j=0;j<n;j++){
	    real tmp=a[j*n+id];
	    a[j*n+id]=a[j*n+ii];
	    a[j*n+ii]=tmp;
	  };
	  ii=n-1;
	};
      };
    };
    
    real a_inv=1./a[tmp];
    for(int i=id+1;i<n;i++){
      int tmp2=i*n+id;
      real val=a[tmp2];
      if(val!=0.){
        real p=val*a_inv;
        int tmp3=id*n+id+1;
        a[tmp2++]=0.;
        for (int j=id+1;j<n;j++){
	  real tmp=a[tmp3++];
	  if(tmp!=0.)a[tmp2]-=p*tmp;
	  tmp2++;
        };
      };
    };
  };

  //int nnn=pos1.size();
  /*
  for(int i=0;i<n;i++){
    cout<<i<<" "<<pos[i]<<" ";
    if(i!=pos[i])cout<<"!";
    cout<<"\n";
  };
  */
};

void lu_decomposition(vector<real> &a, int n)
{
  for(int id=0;id<n-1;id++){
    int idn=id*n;
    int tmp=idn+id;
    real a_inv=1./a[tmp];
    for(int i=id+1;i<n;i++){
      int inid=i*n+id;
      real val=a[inid];
      if(val!=0.){
        real p=val*a_inv;
        int tmp2=inid;
        int tmp3=idn+id+1;
        a[tmp2++]=p; // ! specific for LU decomposition
        for (int j=id+1;j<n;j++){
	  real vtmp=a[tmp3++];
	  if(vtmp!=0.)a[tmp2]-=p*vtmp;
	  tmp2++;
        };
      };
    };
  };
};

void lu_decomposition_inverse(vector<real> &a, int n)
{
  lu_decomposition(a,n);

  vector<real> b(n*n,0.);
  for(int i=0;i<n;i++){b[i*n+i]=1.;};

  // Forward substitution
  for(int i=0;i<n;i++){
    int in=i*n;
    for(int j=0;j<i;j++){
      real tmp2=a[in+j];
      if(tmp2!=0.){
        int jn=j*n;
        for(int ii=0;ii<n;ii++){
          b[in+ii]-=tmp2*b[jn];
	  jn++;
	};
      };
    };
  };

  // Backward substitution
  for(int i=n-1;i>=0;i--){
    int in=i*n;
    for(int j=n-1;j>i;j--){
      real tmp2=a[in+j];
      if(tmp2!=0.){
	int jn=j*n;
        for(int ii=0;ii<n;ii++){
          b[in+ii]-=tmp2*b[jn];
	  jn++;
        };
      };
    };
    real tmp=1./a[in+i];
    for(int ii=0;ii<n;ii++){
      b[in+ii]*=tmp;
      a[in+ii]=b[in+ii];
    };
  };

};

void gauss_3p(real *a,real *b,int n,int m)
{
  real p;

  int tmp1;
  int tmp2=1;
  for(int id=0;id<n-1;id++){
    tmp1=tmp2+2;
    p=a[tmp1]/a[tmp2];
    a[tmp1+1]-=p*a[tmp2+1];
    b[id+1]-=p*b[id];
    tmp2+=3;
  };
  
  b[n-1]/=a[3*n-2];

  tmp1=3*n-4;
  for (int i=n-2;i>=0;i--){
    b[i]=(b[i]-b[i+1]*a[tmp1])/a[tmp1-1];
    tmp1-=3;
  };
}

void gauss22(real *a,real *b)
{
  if(a[0]==0&&a[2]==0){
    cout<<"gauss22 method:Cannot solve.\n";
  }

  int ij=0;
  if(a[0]==0&&a[2]!=0){
    ij=1;
    ChangeValue(a[0],a[1]);
    ChangeValue(a[2],a[3]);
  }

  gauss(a,b,2,1);
  if(ij==1)ChangeValue(b[0],b[1]);
}

void gauss22(vector<real> &a,vector<real> &b)
{
  if(a[0]==0&&a[2]==0){
    cout<<"gauss22 method:Cannot solve.\n";
  }

  int ij=0;
  if(a[0]==0&&a[2]!=0){
    ij=1;
    ChangeValue(a[0],a[1]);
    ChangeValue(a[2],a[3]);
  }

  gauss(a,b,2,1);
  if(ij==1)ChangeValue(b[0],b[1]);
}

void ChangeOrder(real *a,int ip)
{
  int *typ=new int[ip];
  real *aorg=new real[ip];
  for(int i=0;i<ip;i++){
    aorg[i]=a[i];
  };
  ChangeOrder(a,ip,typ);
  for(int i=0;i<ip;i++){
    a[i]=aorg[typ[i]];
  };
  delete [] typ;
  delete [] aorg;
};

void ChangeOrder(real *a, int ip, int *ret)
{
  vector<real> a2(ip);
  vector<int> ret2(ip);
  for(int i=0;i<ip;i++){
    a2[i]=a[i];
  };
  ChangeOrder(a2,ip,ret2);
  for(int i=0;i<ip;i++){
    ret[i]=ret2[i];
  };
};

void ChangeOrder(vector<real>& a,int ip,vector<int> &ret)
{
  real *tmp=new real[ip];
  int  *typ=new int[ip];
  real maxval=1e10;

  if(ip==1){
    ret[0]=0;
    return;
  };

  for(int i=0;i<ip;i++){
    tmp[i]=a[i];
    typ[i]=0;
    if(tmp[i]>maxval){
      cout<<"# Error in ChangeOrder of Numeric.\n";
      cout<<"# Please change a value of maxval.\n";
      exit(0);
    };
  };

  for(int i=0;i<ip;i++){
    real min=maxval;
    int ic=0;
    for(int j=0;j<ip;j++){
      if(tmp[j]<min&&typ[j]==0){
        min=tmp[j];
        ic=j;
      };
    };
    if(min>=maxval){
      cout<<"Error in Change Order.\n";
      exit(0);
    };
    ret[i]=ic;
    typ[ic]=1;
  };

  delete [] tmp;
  delete [] typ;
};

real GetEnx(int m,real x)
{
  real ex=exp(-x);
  int mm1=m-1;
  real xm=real(mm1);
  if(x<=0.)return 1./xm; // En(0)=1/(n-1) when n>1
  if(x<1.){
    real e1=-(log(x)+(((((( -2.8E-5*x+2.31E-4)*x-1.667E-3)*x+1.0417E-2)*x
			 -5.55556E-2)*x+0.25)*x-1.0)*x+0.5772157);
    // Euler's constant : 0.5772156649 
    if(m!=1){
      real sum=0.;
      real xkai=1./xm;
      for(int k=1;k<=mm1;k++){
	sum+=xkai;
	xm-=1.0;
	if(xm==0.)xm=1.;
	xkai=-xkai*x/xm;
      };
      return ex*sum+xkai*e1;
    }
    else{
      return e1;
    };
  }else{
    real ax=(((x+8.5733287)*x+18.059017)*x+8.6347609)*x+0.26777373;
    real bx=(((x+9.5733223)*x+25.632956)*x+21.099653)*x+3.9584969;
    if(m!=1){
      real sum=0.;
      real xkai=1./xm;
      for(int k=1;k<=mm1;k++){
        sum+=xkai;
        xm-=1.;
        if(k!=mm1)xkai=-xkai*x/xm;
      };
      return ex*(sum-xkai*ax/bx);
    }else{
      return ex*ax/bx/x;
    };
  };
  return 0.;
};

int TranslateNuclideIDFromJFS(int i)
{
  switch(i){
  case 1:
    return 125;
  case 2:
    return 128;
  case 307:
    return 328;
  case 4:
    return 425;
  case 115:
    return 528;
  case 6:
    return 600;
  case 147:
    return 725;
  case 8:
    return 825;
  case 9:
    return 925;
  case 11:
    return 1125;
  case 13:
    return 1325;
  case 14:
    return 1400;
  case 15:
    return 1525;
  case 16:
    return 1600;
  case 19:
    return 1900;
  case 20:
    return 2000;
  case 22:
    return 2200;
  case 23:
    return 2300;
  case 24:
    return 2400;
  case 242:
    return 2431;
  case 26:
    return 2600;
  case 266:
    return 2631;
  case 28:
    return 2800;
  case 288:
    return 2825;
  case 280:
    return 2831;
  case 29:
    return 2900;
  case 293:
    return 2925;
  case 295:
    return 2931;
  case 335:
    return 3325;
  case 40:
    return 4000;
  case 400:
    return 4025;
  case 401:
    return 4028;
  case 402:
    return 4031;
  case 404:
    return 4037;
  case 41:
    return 4125;
  case 413:
    return 4125;
  case 42:
    return 4200;
  case 47:
    return 4700;
  case 31:
    return 3100;
  case 25:
    return 2525;
  case 74:
    return 7400;
  case 742:
    return 7431;
  case 743:
    return 7434;
  case 744:
    return 7437;
  case 746:
    return 7443;
  case 902:
    return 9040;
  case 923:
    return 9222;
  case 924:
    return 9225;
  case 925:
    return 9228;
  case 926:
    return 9231;
  case 928:
    return 9237;
  case 937:
    return 9346;
  case 948:
    return 9434;
  case 949:
    return 9437;
  case 940:
    return 9440;
  case 941:
    return 9443;
  case 942:
    return 9446;
  case 951:
    return 9543;
  case 953:
    return 9549;
  case 962:
    return 9631;
  case 963:
    return 9634;
  case 964:
    return 9637;
  case 965:
    return 9640;
  case 966:
    return 9643;
  };
  cout<<"Error in TranslateNuclideIDFromJFS.\n";
  cout<<" JFS ID "<<i<<"\n";
  cout<<"  (See Numeric.cxx.)\n";
  exit(0);
};

void MakeFineMap(int x,int y,int *inp,int xst,int xed,int c1,int *cell1,int *cell2,int *cell3,int *cell4,int *mapnew)
{
  int index1=0;
  int index2=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      if(inp[index1]==-1){
	for(int k=0;k<c1;k++){
	  mapnew[index2]=cell1[k];
	  index2++;
	};
      }else if(inp[index1]==-2){
	for(int k=0;k<c1;k++){
	  mapnew[index2]=cell2[k];
	  index2++;
	};
      }else if(inp[index1]==-3){
	for(int k=0;k<c1;k++){
	  mapnew[index2]=cell3[k];
	  index2++;
	};
      }else if(inp[index1]==-4){
	for(int k=0;k<c1;k++){
	  mapnew[index2]=cell4[k];
	  index2++;
	};
      }else if(j>=xst&&j<=xed){
	for(int k=0;k<c1;k++){
  	  mapnew[index2]=inp[index1];
  	  index2++;
	};
      }else {
	mapnew[index2]=inp[index1];
	index2++;
      };
      index1++;
    };
  };
};

void MakeXMirrorMap(int x,int y,int *inp,int *mapnew)
{
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      mapnew[index++]=inp[i*x+(x-j-1)];
    };
  };
};

void MakeYMirrorMap(int x,int y,int *inp,int *mapnew)
{
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      mapnew[index++]=inp[(y-i-1)*x+j];
    };
  };
};

void WriteXMirrorMap(int x,int y,int *inp)
{
  int *mapnew=new int[x*y];
  MakeXMirrorMap(x,y,inp,mapnew);
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      cout<<""<<mapnew[index++]<<",";
    };
    cout<<"\n";
  };
};

void WriteYMirrorMap(int x,int y,int *inp)
{
  int *mapnew=new int[x*y];
  MakeYMirrorMap(x,y,inp,mapnew);
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      cout<<""<<mapnew[index++]<<",";
    };
    cout<<"\n";
  };
};

real CalAreaMaximizeYPoint(real x1,real y1)
{
  bool err=false;
  if(x1<0.||x1>1.)err=true;
  if(y1<0.||y1>1.)err=true;
  real ymax=sqrt(1.-x1*x1);
  if(ymax<=y1)err=true;
  if(err){
    cout<<"Error in CalAreaMaximizeYPoint.\n";
    exit(0);
  };
  real sold=0.;
  real ret=0.;
  for(int i=99999;i>0;i--){
    real yp=i*0.00001;
    if(yp<ymax&&yp>y1){
      real x2=sqrt(1.-yp*yp);
      real s=(x2-x1)*(yp-y1);
      if(s>sold){
        sold=s;
	ret=yp;
      };
    };
  };
  return ret;
};

void CircleToSquare(real r,int n)
{
  real *ypoint=new real[10000];
  real area=0.;
  int pnum=0;
  for(int i=0;i<n;i++){
    if(i==0){
      real t=CalAreaMaximizeYPoint(0.,0.);
      ypoint[pnum++]=t;
      area+=t*t*0.5;
    }else if(i==1){
      real t=CalAreaMaximizeYPoint(0.,ypoint[0]);
      real t2=sqrt(1-t*t);
      ypoint[pnum++]=t;
      area+=(t-ypoint[0])*t2;
    }else{
      int pnumtmp=pnum;
      real *tmp=new real[pnumtmp];
      for(int j=0;j<pnumtmp;j++){tmp[j]=ypoint[j];};
      real x1=0.;
      for(int j=0;j<pnumtmp;j++){
	real t=tmp[pnumtmp-j-1];
	ypoint[pnum++]=CalAreaMaximizeYPoint(x1,t);
	real t2=ypoint[pnum-1];
	area+=(t2-t)*(sqrt(1-t2*t2)-x1);
	x1=sqrt(1-t*t);
      };
      delete [] tmp;
    };
    ChangeOrder(ypoint,pnum);
  };

  real ratio=sqrt(3.14159265*0.125/area)*r;
  cout<<"    (x & y)\n";
  for(int i=0;i<pnum;i++){
    cout<<i<<" "<<sqrt(1-ypoint[i]*ypoint[i])*ratio<<" "<<ypoint[i]*ratio<<"\n";
  };

  delete [] ypoint;
};

void RedividingEquivolumeRing(real *rinp, int pos)
{
  real volume=PI*rinp[pos]*rinp[pos]/(pos+1);
  
  for(int i=0;i<pos;i++){
    rinp[i]=sqrt(volume*(i+1)/PI);
  };
};

bool CalCubicSplineCoefficient
(real *x,real *y,real *dy,int n,real *c,real *d,real *e)
{
  if(n<2)return false;

  int n1=n-1;
  int n2=n-2;

  real x_org=x[0];
  for(int i=1;i<n;i++){
    real x_new=x[i];
    if(x_org>=x_new)return false;
    x_org=x_new;
  };

  d[0]=dy[0];
  d[n-1]=dy[1];
  if(n!=2){
    int n3=n-3;
    real h1=x[1]-x[0];
    real h=x[2]-x[1];
    real yh1=(y[1]-y[0])/h1;
    real yh=(y[2]-y[1])/h;
    d[1]=0.5/(h1+h);
    c[1]=6.*(yh-yh1)-h1*d[0];

    for(int i=2;i<n1;i++){
      h1=h;
      yh1=yh;
      h=x[i+1]-x[i];
      yh=(y[i+1]-y[i])/h;
      real tm=h1*d[i-1];
      d[i]=1./(2.*(h1+h)-h1*tm);
      e[i]=tm;
      c[i]=6.*(yh-yh1)-tm*c[i-1];
    };

    d[n2]=(c[n2]-h*d[n1])*d[n2];
    for(int i=0;i<n3;i++){
      int i0=n1-i-2;
      int i1=n-i-2;
      d[i0]=c[i0]*d[i0]-d[i1]*e[i1];
    };
  };

  real h=0.;
  real h_inv=0.;
  for(int i=0;i<n1;i++){
    real dw=d[i];
    real dw1=d[i+1];
    h=x[i+1]-x[i];
    h_inv=1./h;
    c[i]=(y[i+1]-y[i])*h_inv-h*(dw1+2.*dw)*0.166666666;
    e[i]=(dw1-dw)*h_inv*0.166666666;
    d[i]=dw*0.5;
  };
  c[n1]=(y[n1]-y[n2])*h_inv+h*(d[n1]+d[n2])*0.333333333;
  d[n1]=d[n1]*0.5;
  e[n1]=0.;

  return true;
};

real CubicSplineFitting(real *x,real *y,int n,real xinp)
{
  real *b=new real[n];
  real *c=new real[n];
  real *d=new real[n];

  real dy[]={0.,0.};
  bool cal=CalCubicSplineCoefficient(x,y,dy,n,b,c,d);
  if(!cal){
    cout<<"# Error in CalCubicSplineCoefficient,\n";
    exit(0);
  };

  int pos=0;
  if(xinp>=x[n-1]){
    pos=n-1;
    c[n-1]=0.;
  }else if(xinp<x[0]){
    pos=0;
    c[0]=0.;
  }else{
    real x_org=x[0];
    for(int i=0;i<n-1;i++){
      real x_new=x[i+1];
      if(xinp>=x_org&&xinp<x_new){
        pos=i;
	break;
      };
      x_org=x_new;
    };
  };
  real dx=xinp-x[pos];
  real dpos=d[pos];
  real cpos=c[pos];
  real bpos=b[pos];

  delete [] b;
  delete [] c;
  delete [] d;

  return y[pos]+((dpos*dx+cpos)*dx+bpos)*dx;
};

real CubicSplineFittingNew(real *x,real *y,int n,real xinp)
{
  int pos=0;
  if(xinp>x[n-1]||xinp<x[0]){
    cout<<"# Error in CubicSplineFitting.\n";
    cout<<"# Input value is out of range.\n";
    exit(0);
  }else{
    real x_org=x[0];
    for(int i=0;i<n-1;i++){
      real x_new=x[i+1];
      if(xinp>=x_org&&xinp<=x_new){
        pos=i;
	break;
      };
      x_org=x_new;
    };
  };

  vector<real> c(n);
  c[0]=0.;
  c[n-1]=0.;
  int maxn=max(pos+1,n-1);

  for(int i=1;i<maxn;i++){
    c[i]=3.*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]));
  };

  vector<real> w(maxn+1);
  w[0]=0.;
  for(int i=1;i<maxn;i++){
    w[i]=2.*(x[i+1]-x[i-1]);
  };
  
  for(int i=1;i<maxn;i++){
    if(i!=n-2){
      real tmp=(x[i+1]-x[i])/w[i];
      w[i+1]-=(x[i+1]-x[i])*tmp;
      c[i+1]-=c[i]*tmp;
    };
    if(i!=1){
      real tmp=(x[i]-x[i-1])/w[i];
      c[i-1]-=c[i]*tmp;
    };
  };

  real tmp=x[pos+1]-x[pos];
  real dx=xinp-x[pos];
  real tmp_inv=1./tmp;
  real cpos=c[pos];
  real cpos1=c[pos+1];
  if(pos!=0)  cpos/=w[pos];
  if(pos!=n-2)cpos1/=w[pos+1];
  real dpos=(cpos1-cpos)*0.333333333333333*tmp_inv;
  real bpos=(y[pos+1]-y[pos])*tmp_inv-(cpos1+2.*cpos)*tmp*0.33333333333333;

  return y[pos]+((dpos*dx+cpos)*dx+bpos)*dx;
};

real LinearLinearInterpolation
(real x1,real x2,real y1,real y2,real x)
{
  real tmp=1./(x1-x2);
  return ((y1-y2)*x+(y2*x1-y1*x2))*tmp;
};

real GetNuclideWeight(int id)
{
  int total=38;
  int nucid[]={
    525,7328,9228,4200,2525, 4800,4125,2900,4000,7400,
    2600,2400,2800,2725,4234, 4240,4243,4249,4525,4634,
    5525,6028,6034,6240,6331, 4731,9237,9040,125,600,
    2200,8200,8325,1200,425, 2300,1400,1325
    // 105,731,925,42,25, 48,413,29,40,74,
    //26,24,28,279,425, 427,428,420,453,465,
    //553,603,605,629,633, 479,928,902,1,6,
    //22,82,839,12,4, 23,14,13

  };
  real weight[]={
    9.92692, 179.394, 233.025, 95.1072, 54.4661,
    111.446, 92.1082, 62.9997, 90.45, 182.269,
    55.367, 51.5492, 58.1838, 58.4269,94.0905,
    96.0735, 97.0643, 99.0492, 102.021,104.004,
    131.763, 141.682, 143.668, 147.638,151.608,
    107.969, 236.006, 230.045, 0.999167, 11.8969,
    47.4818, 205.435, 207.185, 24.0963, 8.93476,
    50.5063, 27.8443, 26.75, 
  };
  for(int i=0;i<total;i++){
    if(nucid[i]==id){return weight[i];};
  };
  cout<<"Error in Numeric::GetNuclideWeight.\n";
  cout<<"No data for ID="<<id<<"\n";
  exit(0);
};

bool SameReal(real a,real b,real f)
{
  real tmp=fabs(a/b-1.);
  if(tmp<f)return true;
  return false;
};

string IntToString(int i)
{
  char chr[256]="chardata";
  sprintf(chr,"%i",i);
  string ret=string(chr);
  return ret;
};

int StringToInt(string inp)
{
  int sz=inp.size();
  int ret=0;
  bool neg=false;
  for(int i=0;i<sz;i++){
    string sp=inp.substr(sz-1-i,1);
    int tmp=0;
    if(sp=="1"){tmp=1;}
    else if(sp=="2"){tmp=2;}
    else if(sp=="3"){tmp=3;}
    else if(sp=="4"){tmp=4;}
    else if(sp=="5"){tmp=5;}
    else if(sp=="6"){tmp=6;}
    else if(sp=="7"){tmp=7;}
    else if(sp=="8"){tmp=8;}
    else if(sp=="9"){tmp=9;}
    else if(sp=="-"){neg=true;};
    ret+=tmp*pow(10.,i);
  };

  if(neg)ret*=-1;
  return ret;
};

real StringToReal(string inp)
{
  if(inp=="0")return 0.;
  
  int sz=inp.size();
  real ret=0.;

  int pos;
  for(int i=0;i<sz;i++){
    string sp=inp.substr(i,1);
    if(sp=="."){pos=i;}
  };

  int epos=-1;
  for(int i=0;i<sz;i++){
    string sp=inp.substr(i,1);
    if(sp=="e"){epos=i;}
  };
  real factor=1.;
  if(epos!=-1){
    int tmp=StringToInt(inp.substr(epos+1,sz-epos-1));
    factor=pow(10,tmp);
  };

  int ppp=pos-1;
  int maxi=sz;
  if(epos>-1)maxi=epos;
  for(int i=0;i<maxi;i++){
    if(i!=pos){
      string sp=inp.substr(i,1);
      real tmp=0;
      if(sp=="1"){tmp=1;}
      else if(sp=="2"){tmp=2;}
      else if(sp=="3"){tmp=3;}
      else if(sp=="4"){tmp=4;}
      else if(sp=="5"){tmp=5;}
      else if(sp=="6"){tmp=6;}
      else if(sp=="7"){tmp=7;}
      else if(sp=="8"){tmp=8;}
      else if(sp=="9"){tmp=9;}
      ret+=tmp*pow(10.,ppp);
      ppp-=1;
    };
  };

  return ret*factor;
};

// Conjugate gradient

void cg(real *a,real *b,int n)
{
  real *x=new real[n];
  real *res=new real[n];
  real *p=new real[n];

  for(int i=0;i<n;i++){
    x[i]=b[i];
  };

  for(int i=0;i<n;i++){
    real tmp=0.;
    for(int j=0;j<n;j++){
      tmp+=a[i*n+j]*x[j];
    };
    res[i]=b[i]-tmp;
  };

  for(int i=0;i<n;i++){
    p[i]=res[i];
  };

  for(int k=0;k<n;k++){

    real tmp1=0.;
    for(int i=0;i<n;i++){
      tmp1+=p[i]*res[i];
    };
    real tmp2=0.;
    for(int i=0;i<n;i++){
      real tmp3=0.;
      for(int j=0;j<n;j++){
	tmp3+=a[i*n+j]*p[j];
      };
      tmp2+=p[i]*tmp3;
    };
    real alp=tmp1/tmp2;

    for(int i=0;i<n;i++){
      x[i]+=alp*p[i];
    };

    for(int i=0;i<n;i++){
      real tmp3=0.;
      for(int j=0;j<n;j++){
	tmp3+=a[i*n+j]*p[j];
      };
      res[i]-=alp*tmp3;
    };

    tmp1=0.;
    tmp2=0.;
    for(int i=0;i<n;i++){
      real tmp3=0.;
      for(int j=0;j<n;j++){
	tmp3+=a[i*n+j]*p[j];
      };
      tmp1+=res[i]*tmp3;
      tmp2+=p[i]*tmp3;
    };
    real bet=-tmp1/tmp2;

    for(int i=0;i<n;i++){
      p[i]=res[i]+bet*p[i];
    };

  };

  for(int i=0;i<n;i++){
    b[i]=x[i];
  };

  delete [] p;
  delete [] x;
  delete [] res;
};

void bicg(real *a,real *b,int n)
{
  real eps=1e-2;

  real *x=new real[n];
  real *res=new real[n];
  real *resp=new real[n];
  real *p=new real[n];
  real *pp=new real[n];
  real *ap=new real[n];

  for(int i=0;i<n;i++){
    //x[i]=b[i];
    x[i]=0.;
  };

  for(int i=0;i<n;i++){
    real tmp=0.;
    for(int j=0;j<n;j++){
      tmp+=a[i*n+j]*x[j];
    };
    res[i]=b[i]-tmp;
  };

  for(int i=0;i<n;i++){
    resp[i]=res[i];
  };

  real bet=0.;
  real betnum=0.;
  for(int k=0;k<n;k++){

    if(k!=0){
      for(int i=0;i<n;i++){
        p[i]=res[i]+bet*p[i];
	pp[i]=resp[i]+bet*pp[i];
      };
    }else{
      for(int i=0;i<n;i++){
	p[i]=res[i];
	pp[i]=resp[i];
	betnum+=res[i]*resp[i];
      };
    };

    real tmp2=0.;
    for(int i=0;i<n;i++){
      real tmp3=0.;
      for(int j=0;j<n;j++){
	tmp3+=a[i*n+j]*p[j];
      };
      ap[i]=tmp3;
      tmp2+=pp[i]*tmp3;
    };

    if(tmp2==0.){
      cout<<"alpha break-down in BiCG.\n";
      exit(0);
    };
    real alp=betnum/tmp2;

    real errmax=0.;
    for(int i=0;i<n;i++){
      real err=alp*p[i]/x[i];
      if(fabs(err)>errmax)errmax=fabs(err);
      x[i]+=alp*p[i];
    };
    if(errmax<eps){
      cout<<"BiCG is converged. ("<<k<<"/"<<n<<")\n";
      break;
    };

    for(int i=0;i<n;i++){
      real tmp4=0.;
      for(int j=0;j<n;j++){
	tmp4+=a[j*n+i]*pp[j];
      };
      res[i]-=alp*ap[i];
      resp[i]-=alp*tmp4;
    };

    real tmp1=0.;
    for(int i=0;i<n;i++){
      tmp1+=res[i]*resp[i];
    };
    if(betnum==0.){
      cout<<"beta break-down in BiCG.\n";
      exit(0);
    };
    bet=tmp1/betnum;
    betnum=tmp1;

  };

  for(int i=0;i<n;i++){
    b[i]=x[i];
  };

  delete [] p;
  delete [] pp;
  delete [] x;
  delete [] res;
  delete [] resp;
  delete [] ap;
};

void MatrixDiagonalization(real *a, real *u, int n)
{
  //
  // Jacobi diagonalization method
  //

  real max_val=0.;
  int count_non_zero=0;
  for(int i=0;i<n*n;i++){
    if(fabs(a[i])>max_val)max_val=fabs(a[i]);
    if(fabs(a[i])>0.)count_non_zero++;
  };

  if(max_val==0.){
    cout<<"# Error in MatrixDiagonalization.\n";
    cout<<"# Matrix is zero.\n";
    exit(0);
  };

  // Check for matrix symmetry
  bool err=false;
  for(int i=0;i<n;i++){
    for(int j=i;j<n;j++){
      real t1=a[i*n+j];
      real t2=a[j*n+i];
      if(t2==0.){
	if(t1!=0.)err=true;
      }else{
        real dif=fabs(t1-t2)/t2;
	if(dif>1e-8)err=true;
      };
    };
  };
  if(err){
    cout<<"\n#Error in MatrixDiagonalization.\n";
    cout<<"# Matrix is not symmetric.\n";
    exit(0);
  };

  real *uold=new real[n*n];
  real *utmp=new real[n*n];
  real *anew=new real[n*n];
  
  int itermax=100000000;
  real eps=max_val*1e-20;
  //real eps=max_val*1e-40;  

  //real eps=1e-10;
  //real eps=1e-20;
  //real eps=1e-25;
  //real eps=1e-40;

  // initialization for u_mat
  int id=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(i==j){
	u[id]=1.;
      }else{
        u[id]=0.;
      };
      id++;
    };
  };

  // Iteration
  for(int iter=0;iter<itermax;iter++){

    // store old u
    int id=0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	uold[id]=u[id];
	if(i==j){
	  utmp[id]=1.;
	}else{
	  utmp[id]=0.;
	};
	id++;
      };
    };

    // Search for maximum absolute value
    int maxi=-1;
    int maxj=-1;
    real maxv=0.;
    int ii=0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        if(i!=j){
	  real v=fabs(a[ii]);
	  if(v>maxv){
	    maxi=i;
	    maxj=j;
	    maxv=v;
	  };
        };
	ii++;
      };
    };
    if(maxi==-1||maxj==-1){
      //if(iter!=0&&n-1!=iter){
      if(iter!=0&&n-1!=iter&&iter<count_non_zero-n){	
        cout<<"# Warning in MatrixDiagonalization.\n";
        cout<<"# Matrix becomes diagonal at iteration "<<iter<<"\n";
	cout<<"# The size of matrix is "<<n<<"\n";
	cout<<"# The number of non-zero elements is "<<count_non_zero<<"\n";
      };
      break;
    };

    //cout<<iter<<" "<<maxv<<" for "<<eps<<"\n";
    if(maxv<eps){
      //cout<<maxv<<" "<<iter<<"\n";
      break;
    };

    real aii=a[maxi*n+maxi];
    real ajj=a[maxj*n+maxj];
    real theta;
    if(aii==ajj){
      theta=PI*0.25;
    }else{
      theta=atan(-2*a[maxi*n+maxj]/(aii-ajj))*0.5;
    };
    real cost=cos(theta);
    real sint=sin(theta);

    // Put matrix u
    int tmp1=maxi;
    int tmp2=maxj;
    for(int i=0;i<n;i++){
      u[tmp1]=uold[tmp1]*cost-uold[tmp2]*sint;
      u[tmp2]=uold[tmp1]*sint+uold[tmp2]*cost;
      tmp1+=n;
      tmp2+=n;
    };

    int tmpi=maxi*n+maxi;
    int tmpj=maxj*n+maxj;
    int tmpij=maxi*n+maxj;
    real cc=cost*cost;
    real ss=sint*sint;
    real cs=2.*cost*sint*a[tmpij];
    anew[tmpi]=a[tmpi]*cc+a[tmpj]*ss-cs;
    anew[tmpj]=a[tmpi]*ss+a[tmpj]*cc+cs;
    anew[maxi*n+maxj]=anew[maxj*n+maxi]=0.;
    for(int i=0;i<n;i++){
      if(i!=maxi&&i!=maxj){
	anew[maxi*n+i]=a[maxi*n+i]*cost-a[maxj*n+i]*sint;
	anew[i*n+maxi]=anew[maxi*n+i];
	anew[maxj*n+i]=a[maxi*n+i]*sint+a[maxj*n+i]*cost;
	anew[i*n+maxj]=anew[maxj*n+i];
      };
    };

    for(int i=0;i<n;i++){
      if(i!=maxi&&i!=maxj){
        for(int j=0;j<n;j++){
  	  if(j!=maxi&&j!=maxj){
	    anew[i*n+j]=a[i*n+j];
	  };
	};
      };
    };

    id=0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	a[id]=anew[id];
	id++;
      };
    };

    if(iter==itermax-1){
      cout<<"# The number of loops reach itermax.\n";
    };

  };

  delete [] uold;
  delete [] utmp;
  delete [] anew;
};

// Random number generator

double ransu()
{
  //return (drand48()-0.5)*2.;
  double i=double(random())*fact;
  return (i-0.5)*2.;
};

double ransu01()
{
  //return (drand48()-0.5)*2.;
  return double(random())*fact;
};

real ransu_gauss(real mu, real sigma)
{
  /*
  random();
  random();
  */

  real alpha=double(random())*fact;
  real beta=double(random())*fact*PI2;
  
  real BoxMuller1=sqrt(-2*log(alpha));
  real BoxMuller2=sin(beta);

  return sigma*(BoxMuller1*BoxMuller2)+mu;
};

void GetMeanAndVariance(vector<real> data, int num, real &mean, real &var)
{
  mean=0.;
  for(int i=0;i<num;i++){
    mean+=data[i];
  };
  mean/=num;

  var=0.;
  for(int i=0;i<num;i++){
    var+=pow(data[i]-mean,2);
  };
  var/=num-1;
};

real GetStatics(vector<real> data, int num, bool print, bool positivity, int print_div)
{
  real maxv=data[0];
  real minv=data[0];
  for(int i=0;i<num;i++){
    if(data[i]>maxv)maxv=data[i];
    if(data[i]<minv)minv=data[i];
  };

  real vnum=real(num);

  real avg=0.;
  for(int i=0;i<num;i++){
    avg+=data[i];
  };
  avg/=num;

  real var=0.;
  real skew=0.;
  real kurto=0.;
  for(int i=0;i<num;i++){
    var+=pow(data[i]-avg,2);
    skew+=pow(data[i]-avg,3);
    kurto+=pow(data[i]-avg,4);
  };
  real var2=var/num;
  var/=num-1;
  real stdev=sqrt(var);

  skew*=vnum;
  skew/=(vnum-1)*(vnum-2);  // unbiased
  skew/=pow(stdev,3);

  /*
  kurto/=vnum-1; // biased
  kurto/=pow(stdev,4);
  kurto-=3;
  */

  kurto*=vnum*(vnum+1);
  kurto/=(vnum-1)*(vnum-2)*(vnum-3);
  kurto/=pow(stdev,4);

  real tmp=3.*(vnum-1)*(vnum-1);
  tmp/=(vnum-2)*(vnum-3);  // unbiased
  kurto-=tmp;



  cout.setf(ios::scientific);
  cout.precision(5);
  //if(print){
    cout<<"# Number of samples    : "<<num<<"\n";
    cout<<"# Maximum              : "<<maxv<<"\n";
    cout<<"# Minimum              : "<<minv<<"\n";
    cout<<"# Expectation          : "<<avg<<"\n";
    cout<<"# Variance (unbiased)  : "<<var<<"\n";
    cout<<"# Std. Dev. (unbiased) : "<<sqrt(var)<<"\n"; 
    cout<<"# Skewness             : "<<skew<<"\n";
    cout<<"# Kurtosis             : "<<kurto<<"\n";
    cout<<"# (Std. Dev. of expectation      : "<<sqrt(var/num)<<")\n";
    cout<<"# (Rel. Std. Dev. of expectation : "<<sqrt(var/num)/avg<<")\n";
    //};

  /*
  cout<<"# Number of   Expec-      Vari-       Std.        Skew-   Kurto-\n";
  cout<<"# samples     tation      ance        Dev.        ness    sis\n";
  real sum1=0.;
  real sum2=0.;
  real sum3=0.;
  real sum4=0.;
  real sum5=0.;
  for(int i=0;i<num;i++){
    sum1+=data[i];
    sum2+=pow(data[i],2);
    sum3+=pow(data[i],3);
    sum4+=pow(data[i],4);
    if(i!=0){
      cout<<"     "<<i+1<<"    ";
      int pnt=i+1;
      real avg=sum1/pnt;
      real mom2=sum2/pnt-pow(avg,2);
      real mom3=sum3/pnt-3*avg*sum2/pnt+2*pow(avg,3);
      real mom4=sum4/pnt-4*avg*sum3/pnt+6*pow(avg,2)*sum2/pnt-3*pow(avg,4);
      mom2*=real(pnt)/(pnt-1);
      mom3*=real(pnt)/(pnt-1);
      mom4*=real(pnt)/(pnt-1);
      cout<<avg<<" ";
      cout<<mom2<<" ";
      real stdev=sqrt(mom2);
      cout<<stdev<<" ";
      cout<<mom3/pow(stdev,3)<<" ";
      cout<<mom4/pow(stdev,4)<<" ";
      cout<<data[i]<<" ";
      cout<<"\n";
    };
  };

  */

  if(print){
    real std=sqrt(var);
    cout<<"#\n# Frequency distribution\n#\n";
    cout<<"# Lower     Higher      Relative     Cumulative  Normal\n";
    cout<<"#                       frequency    frequency   Dist. \n";
    cout<<"#\n";
    int div=print_div;
    real vmin=avg-std*6.;
    if(positivity)vmin=0.;
    real vmax=avg+std*6.;
#if 0
    vmin=0.;
    vmax=1.;
    div=100;
#endif
    real leng=vmax-vmin;
    real unit=leng/div;
    real vcur=vmin;
    real cons1=1./sqrt(PI2)/std;
    real cons2=0.5/var;
    real cons3=1./num;
    real cum=0.;
    real cum_fx2=0.;
    real cum_fx2_fix=0.;    
    for(int i=0;i<=div;i++){
      real vcur2=vcur+unit;
      int count=0;
      for(int j=0;j<num;j++){
	if(data[j]>=vcur&&data[j]<vcur2)count++;
      };
      cum+=real(count)*cons3;
      cout<<vcur<<" "<<vcur2<<" "<<real(count)*cons3<<" ";
      cout<<cum<<" ";

      // ... normal dist.
      real vcent=(vcur+vcur2)*0.5;
      real fx=cons1*exp(-(vcent-avg)*(vcent-avg)*cons2)*unit;
      cout<<" "<<fx<<" ";

#if 0      
      // ... log-normal dist.
      //real vcent_geom=sqrt(vcur*vcur2);
      //if(vcur==0.)vcent_geom=sqrt(1e-15*vcur2);
      real vcent_geom=vcent;
      real cons5=log(1.+var/(avg*avg));
      real sigma_ln=sqrt(cons5);
      real mu_ln=log(avg)-0.5*cons5;
      real cons4=pow(log(vcent_geom)-mu_ln,2)*0.5/(sigma_ln*sigma_ln);
      real fx2=1./(sqrt(2*PI)*sigma_ln*vcent_geom)*exp(-cons4)*unit;
      //cout<<cons1<<" "<<vcent_geom<<" "<<cons4<<" "<<unit<<" ";
      cout<<" "<<fx2<<" ";
      cum_fx2+=fx2;

      // ... log-normal dist.
      //real var_fix=10.;
      //real avg_fix=-10.;
      //real cons5_fix=log(1.+var_fix/(avg_fix*avg_fix));
      //real sigma_ln_fix=sqrt(cons5_fix);
      //real mu_ln_fix=log(avg_fix)-0.5*cons5_fix;
      real sigma_ln_fix=sqrt(10.);
      real mu_ln_fix=-10.;
      real cons4_fix=pow(log(vcent_geom)-mu_ln_fix,2)*0.5/(sigma_ln_fix*sigma_ln_fix);
      real fx2_fix=1./(sqrt(2*PI)*sigma_ln_fix*vcent_geom)*exp(-cons4_fix)*unit;
      //cout<<cons1<<" "<<vcent_geom<<" "<<cons4<<" "<<unit<<" ";
      cout<<" "<<fx2_fix<<" ";
      cum_fx2_fix+=fx2_fix;

      cout<<sigma_ln*sigma_ln<<" "<<mu_ln<<"\n";
#endif      
      
      cout<<"\n";
      vcur=vcur2;
    };
#if 0    
    cout<<"# The cumulative probability for log-normal       : "<<cum_fx2<<"\n";
    cout<<"# The cumulative probability for fixed log-normal : "<<cum_fx2_fix<<"\n";
#endif    
  };

  return var;

};

void GetStatics(vector<real> &data, vector<real> &result)
{
  int sz=data.size();
  GetStatics(data,sz,result);
};

void GetStatics(vector<real> &data, int st, int ed, vector<real> &result)
{
  int sz=data.size();
  if(st<0||st>=sz){
    cout<<"# Error in GetStatics.\n";
    cout<<"# Starting position "<<st<<" is inappropriate.\n";
    exit(0);
  };
  if(ed<0||ed<=st||ed>=sz){
    cout<<"# Error in GetStatics.\n";
    cout<<"# Ending position "<<ed<<" is inappropriate.\n";
    exit(0);
  };

  vector<real> dummy;
  for(int i=st;i<=ed;i++){
    dummy.push_back(data[i]);
  };
  
  GetStatics(dummy,ed-st+1,result);
};

void GetStatics(vector<real> &data, int num, vector<real> &result)  
{
  real avg=0.;
  real max=data[0];
  real min=data[0];

  real vnum=real(num);

  for(int i=0;i<num;i++){
    avg+=data[i];
    if(data[i]>max)max=data[i];
    if(data[i]<min)min=data[i];
  };
  avg/=vnum;

  real var=0.;
  real skew=0.;
  real kurto=0.;
  for(int i=0;i<num;i++){
    var+=pow(data[i]-avg,2);
    skew+=pow(data[i]-avg,3);
    kurto+=pow(data[i]-avg,4);
  };
  real var2=var/vnum;
  var/=vnum-1;
  real stdev=sqrt(var);

  /*
  skew/=num-1;  // biased
  skew/=pow(stdev,3);
  */

  skew*=vnum;
  skew/=(vnum-1)*(vnum-2);  // unbiased
  skew/=pow(stdev,3);

  /*
  kurto/=vnum-1; // biased
  kurto/=pow(stdev,4);
  kurto-=3;
  */

  kurto*=vnum*(vnum+1);
  kurto/=(vnum-1)*(vnum-2)*(vnum-3);
  kurto/=pow(stdev,4);

  real tmp=3.*(vnum-1)*(vnum-1);
  tmp/=(vnum-2)*(vnum-3);  // unbiased
  kurto-=tmp;


  result.resize(6);
  result[0]=avg;
  result[1]=var;
  result[2]=skew;
  result[3]=kurto;
  result[4]=max;
  result[5]=min;
};

real GetStatics(vector<real> data1, vector<real> data2, int num)
{
  real avg1=0.;
  real avg2=0.;
  for(int i=0;i<num;i++){
    avg1+=data1[i];
    avg2+=data2[i];
  };
  avg1/=num;
  avg2/=num;

  real var1=0.;
  real var2=0.;
  real cov=0.;
  for(int i=0;i<num;i++){
    var1+=pow(data1[i]-avg1,2);
    var2+=pow(data2[i]-avg2,2);
    cov+=(data1[i]-avg1)*(data2[i]-avg2);
  };

  var1/=num-1;
  var2/=num-1;
  cov/=num-1;

  real stdev1=sqrt(var1);
  real stdev2=sqrt(var2);

  /*
  cout.setf(ios::scientific);
  cout.precision(5);

  cout<<"# Number of samples    : "<<num<<"\n";
  cout<<"# Variance1 (unbiased)  : "<<var1<<"\n";
  cout<<"# Variance2 (unbiased)  : "<<var2<<"\n";
  cout<<"# Covariance (unbiased) : "<<cov<<"\n";
  cout<<"# Correlation           : "<<cov/stdev1/stdev2<<"\n";
  */

  return cov/stdev1/stdev2;

};

void GetMeanAndVarianceForLognormal(real mu_org, real var_org, real &mu_log, real &var_log)
{
  var_log=log(1.+var_org/(pow(mu_org,2)));
  mu_log=log(mu_org)-0.5*var_log;             
};

// +++ Polynomial handler

real GetValue(vector<real> &z, real x)
{
  int order=z.size();
  real sol=0.;
  for(int i=0;i<order;i++){
    sol+=z[i]*pow(x,i);
  };
  return sol;
};

real SolveZeroPoint(vector<real> &z, real vl,real vh)
{
  int itermax=10000;
  real eps=1e-25;
  int order=z.size();

  real vc;
  for(int i=0;i<itermax;i++){
    real tmp=(vl+vh)*0.5;
    if(fabs(tmp/vc-1.)<eps)return vc;
    vc=tmp;
    real soll=GetValue(z,vl);
    real solh=GetValue(z,vh);
    real solc=GetValue(z,vc);
    if(solc*soll>0.){
      vl=vc;
    }else{
      vh=vc;
    };
  };

  return vc;
};

void PolynomialDivision(vector<real> &z,real &p,real &q,vector<real> &b,vector<real> &c)
{
  // 
  // (z[n]x^n + z[n-1]x^{n-1} + ... + z[1]x + z[0] )
  // = (x^2 + p x + q) ( x^{n-2} + b[n-3] x^{n-3] + ... + b[1] x + b[0] ) + ( c[1] x + c[0] )
  // = (x^2 + p x + q) ( x^{n-2} + b[n-3] x^{n-3] + ... + b[1] x + b[0] ) + ( r    x + s    )
  //
  // !! CAUTION !!
  //
  // z[n] should be 1.0

  int sz=z.size(); // (n+1)
  int order=sz;
  int max_order=sz-1; // (n)
  if(z[max_order]!=1.){
    cout<<"# Error in PolynomialDivision.\n";
    exit(0);
  };

  vector<real> btmp(sz-1);

  b.resize(sz-3,0.); 
  c.resize(2,0.);

  if(max_order==2){
    c[1]=z[sz-2]-p;
    c[0]=-q;
    return;
  };

  btmp[0]=z[order-2]-p;
  btmp[1]=z[order-3]-p*btmp[0]-q;
  for(int i=2;i<order-2;i++){
    btmp[i]=z[order-2-i]-p*btmp[i-1]-q*btmp[i-2];
  };
  btmp[order-2]=z[0]-q*btmp[order-4];

  for(int i=0;i<sz-3;i++){
    b[i]=btmp[order-4-i];
  };
  c[1]=btmp[order-3];
  c[0]=btmp[order-2];
};

void BairstowHitchcock(vector<real> &z,real &p,real &q)
{
  int order=z.size()-1; // original polynomial order

  //cout<<"# Bairstow Hitchcock\n";
  //cout<<"# (p,q) : ("<<p<<" , "<<q<<" ) -> ";
  for(int ii=0;ii<10;ii++){

  vector<real> b;
  vector<real> c;
  PolynomialDivision(z,p,q,b,c);

  real s=c[0];
  real r=c[1];

  vector<real> g(order); // (order-1)
  vector<real> bn;
  vector<real> cn;
  g[0]=0.;
  for(int i=1;i<order-1;i++){
    g[i]=b[i-1];
  };
  g[order-1]=1.;

  PolynomialDivision(g,p,q,bn,cn);
  real dsdp=-cn[0];
  real drdp=-cn[1];

  real dsdq=0.;
  real drdq=0.;
  if(order==3){
    dsdq=-b[0];
    drdq=-1.0;
  }else{
    vector<real> g2(order-1); // (order-2)
    for(int i=0;i<order-2;i++){
      g2[i]=b[i];
    };
    g2[order-2]=1.;
    PolynomialDivision(g2,p,q,bn,cn);
    dsdq=-cn[0]; 
    drdq=-cn[1];
  };

  real a00=drdp;
  real a01=drdq;
  real a10=dsdp;
  real a11=dsdq;
  real det_inv=1./(a00*a11-a01*a10);
  real am00=a11*det_inv;
  real am11=a00*det_inv;
  real am01=-a01*det_inv;
  real am10=-a10*det_inv;
  real dp=am00*-r+am01*-s;
  real dq=am10*-r+am11*-s;

  cout.setf(ios::scientific);
  cout.precision(5);
  //cout<<p<<" "<<q<<" "<<dp<<" "<<dq<<" "<<r<<" "<<s<<"\n";
  p+=dp;
  q+=dq;

  };

  //cout<<" ("<<p<<" , "<<q<<" )\n";
};

// Output format

void WriteOut(int a,int col)
{
  int cola=0;

  int ord=0;
  int val=1;
  while(cola==0){
    val*=10;
    ord++;
    if(a<val)cola=ord;
  };

  int spc=col-cola;
  if(spc<0){
    cout<<"# Error in WriteOut.\n";
    exit(0);
  };

  for(int i=0;i<spc;i++)cout<<" ";
  cout<<a;
};

void WriteOut(string a,int col)
{
  int cola=a.size();
  
  int spc=col-cola;
  if(spc<0){
    cout<<"# Error in WriteOut.\n";
    exit(0);
  };

  cout<<a;
  for(int i=0;i<spc;i++)cout<<" ";
};

void WriteOut(real a,string b)
{
  printf(b.data(),a);
};

void NoiseReductionBySVD_org(string filename, string outname, real dx, int m, int pnt)
{
  // m   : dimensions of delayed coordinate space
  // pnt : the number of data for averaging
  // dx

  // ++++++++++++++++++++++++++++++++++++++++

  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# Name : "<<filename<<"\n";
    exit(1);
  };

  ofstream fout;
  fout.open(outname.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# Name : "<<outname<<"\n";
    exit(1);
  };

  vector<real> a;

  int n=0;
  int n_org=0;
  int nn=0;
  real sum=0.;

  while(!fin.eof()){
    real tmp;
    fin>>tmp;
    sum+=tmp;
    nn++;
    n_org++;
    if(nn==pnt){
      sum/=pnt;
      a.push_back(sum);
      sum=0.;
      nn=0.;
      n++;    
    };
    //a.push_back(tmp);
  };

  /*
  fout<<"##########################################\n";
  fout<<"#\n";                                       
  fout<<"# Noise reduction by SVD\n";
  fout<<"#\n";
  fout<<"#   Number of original data points : "<<n_org<<"\n";
  fout<<"#   Number of reduced data points  : "<<n<<"\n";
  fout<<"#     (Number of averaging points  : "<<pnt<<" )\n";
  fout<<"#   Number of SVD order            : "<<m<<"\n";
  fout<<"#\n";
  fout<<"##########################################\n";
  */

  cout<<"##########################################\n";
  cout<<"#\n";                                       
  cout<<"# Noise reduction by SVD\n";
  cout<<"#\n";
  cout<<"#   Number of original data points : "<<n_org<<"\n";
  cout<<"#   Number of reduced data points  : "<<n<<"\n";
  cout<<"#     (Number of averaging points  : "<<pnt<<" )\n";
  cout<<"#   Number of SVD order            : "<<m<<"\n";
  cout<<"#\n";
  cout<<"##########################################\n";
  
  if(n<m){
    cout<<"# Error!\n";
    cout<<"# Space dimensions are smaller than the number of data points.\n";
    exit(0);
  };

  int num=n-m+1;
  real n_inv=1./num;

  real *cmat=new real[m*m];
  for(int i=0;i<m*m;i++){
    cmat[i]=0.;
  };

  for(int i=0;i<num;i++){
    int index=0;
    for(int j=0;j<m;j++){
      real tmp=a[i+j];
      for(int k=0;k<m;k++){
	cmat[index++]+=tmp*a[i+k];
      };
    };    
  };

  for(int i=0;i<m*m;i++){
    cmat[i]*=n_inv;
  };

  real *umat=new real[m*m];
  MatrixDiagonalization(cmat,umat,m);

  real *eigen=new real[m];
  for(int i=0;i<m;i++){
    eigen[i]=cmat[i*m+i];
  };

  int *order=new int[m];
  ChangeOrder(eigen,m,order);

  int mp=min(m,10);

  //
  /*
  real dxx=dx*pnt;
  real x=dxx*0.5;
  for(int i=0;i<num;i++){
    cout<<x<<" ";
    real sum=0.;
    cout<<a[i]<<" ";
    for(int jj=0;jj<mp;jj++){ // 
      int j=order[m-jj-1];
      real tmp=0.;
      for(int k=0;k<m;k++){
	tmp+=umat[j+k*m]*a[i+k];
      };
      tmp*=umat[j]; // the first entry of the vector u_m
      sum+=tmp;
      cout<<sum<<" ";
    };
    x+=dxx;
    cout<<"\n";
  };
  */

  /*
  for(int i=0;i<mp;i++){
    int j=order[m-i-1];
    cout<<umat[j]<<" "<<eigen[j]<<"\n";
  };
  exit(0);
  */

  real dxx=dx*pnt;
  real x=dxx*0.5;
  for(int i=0;i<n;i++){
    fout<<x<<" ";
    real sum=0.;
    fout<<a[i]<<" ";
    int ii=i;
    int j1=0;
    if(i>=num){
      ii=num-1;
      j1=i-num;
    };
    for(int jj=0;jj<mp;jj++){ // 
      int j=order[m-jj-1];
      real tmp=0.;
      for(int k=0;k<m;k++){
	tmp+=umat[j+k*m]*a[ii+k];
      };
      tmp*=umat[j1*m+j]; // the first entry of the vector u_m (i<num)
      sum+=tmp;
      fout<<sum<<" ";
      //fout<<tmp<<" ";
    };
    x+=dxx;
    fout<<"\n";
  };


  delete [] cmat;
  delete [] umat;
  delete [] eigen;
  delete [] order;

  return;
}


void NoiseReductionBySVD(string filename, string outname, real dx, int m, int pnt)
{
  // m   : dimensions of delayed coordinate space
  // pnt : the number of data for averaging
  // dx

  // ++++++++++++++++++++++++++++++++++++++++

  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# Name : "<<filename<<"\n";
    exit(1);
  };

  ofstream fout;
  fout.open(outname.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# Name : "<<outname<<"\n";
    exit(1);
  };
  
  vector<real> data_in;
  while(!fin.eof()){
    real tmp;
    fin>>tmp;
    data_in.push_back(tmp);
  };

  vector< vector<real> > data_out=NoiseReductionBySVD(data_in, dx, m, pnt);

  int n=data_out.size();
  {
  for(int i=0;i<n;i++){
    fout<<data_out[i][0]<<" ";
    fout<<data_out[i][1]<<" ";
    int mp=data_out[i].size()-2;
    for(int jj=0;jj<mp;jj++){ // 
      fout<<data_out[i][2+jj]<<" ";
    };
    fout<<"\n";
  };
  };

  
  return;
}


vector< vector<real> > NoiseReductionBySVD(vector<real> &data_in, int st, int ed, real dx, int m, int pnt)
{
  int sz=data_in.size();
  if(st<0||st>sz-1||st>ed){
    cout<<"# Error in NoiseReductionBySVD.\n";
    cout<<"# Starting data point "<<st<<" is inappropriate.\n";
    exit(0);
  };
  if(ed<0||ed>sz-1){
    cout<<"# Error in NoiseReductionBySVD.\n";
    cout<<"# End data point "<<ed<<" is inappropriate.\n";
    exit(0);
  };

  vector<real> data_in2;
  for(int i=st;i<=ed;i++){
    data_in2.push_back(data_in[i]);
  };

  return NoiseReductionBySVD(data_in2,dx,m,pnt);
};


vector< vector<real> > NoiseReductionBySVD(vector<real> &data_in, real dx, int m, int pnt)    
{
  // m   : dimensions of delayed coordinate space
  // pnt : the number of data for averaging
  // dx

  // ++++++++++++++++++++++++++++++++++++++++

  vector<real> a;

  int n_org=data_in.size();
  int n=0;
  int nn=0;
  real sum=0.;

  for(int i=0;i<n_org;i++){
    sum+=data_in[i];
    nn++;
    if(nn==pnt){
      sum/=pnt;
      a.push_back(sum);
      sum=0.;
      nn=0.;
      n++;    
    };
  };

  cout<<"##########################################\n";
  cout<<"#\n";                                       
  cout<<"# Noise reduction by SVD\n";
  cout<<"#\n";
  cout<<"#   Number of original data points : "<<n_org<<"\n";
  cout<<"#   Number of reduced data points  : "<<n<<"\n";
  cout<<"#     (Number of averaging points  : "<<pnt<<" )\n";
  cout<<"#   Number of SVD order            : "<<m<<"\n";
  cout<<"#\n";
  cout<<"##########################################\n";
  
  if(n<m){
    cout<<"# Error!\n";
    cout<<"# Space dimensions are smaller than the number of data points.\n";
    exit(0);
  };

  int num=n-m+1;
  real n_inv=1./num;

  real *cmat=new real[m*m];
  for(int i=0;i<m*m;i++){
    cmat[i]=0.;
  };

  for(int i=0;i<num;i++){
    int index=0;
    for(int j=0;j<m;j++){
      real tmp=a[i+j];
      for(int k=0;k<m;k++){
	cmat[index++]+=tmp*a[i+k];
      };
    };    
  };

  for(int i=0;i<m*m;i++){
    cmat[i]*=n_inv;
  };

  real *umat=new real[m*m];
  MatrixDiagonalization(cmat,umat,m);

  real *eigen=new real[m];
  for(int i=0;i<m;i++){
    eigen[i]=cmat[i*m+i];
  };

  int *order=new int[m];
  ChangeOrder(eigen,m,order);

  int mp=min(m,10);

  vector< vector<real> > data_out;

  {
  real dxx=dx*pnt;
  real x=dxx*0.5;
  for(int i=0;i<n;i++){
    vector<real> tmp_out;
    //fout<<x<<" ";
    tmp_out.push_back(x);
    real sum=0.;
    //fout<<a[i]<<" ";
    tmp_out.push_back(a[i]);    
    int ii=i;
    int j1=0;
    if(i>=num){
      ii=num-1;
      j1=i-num;
    };
    for(int jj=0;jj<mp;jj++){ // 
      int j=order[m-jj-1];
      real tmp=0.;
      for(int k=0;k<m;k++){
	tmp+=umat[j+k*m]*a[ii+k];
      };
      tmp*=umat[j1*m+j]; // the first entry of the vector u_m (i<num)
      sum+=tmp;
      //fout<<sum<<" ";
      tmp_out.push_back(sum);      
      //fout<<tmp<<" ";
    };
    x+=dxx;
    //fout<<"\n";
    data_out.push_back(tmp_out);
  };
  };

  delete [] cmat;
  delete [] umat;
  delete [] eigen;
  delete [] order;

  return data_out;
}


enum xstype XSType(int i)
{
  switch(i){
  case 1:
    return sigt;
  case 18:
    return sigf;
  case 102:
    return sigc;
  case 452:
    return nu;
  case 2:
    return sigel;
  case 251:
    return mu;
  case 181:
    return chi;
  case 4:
    return siginel;
  case 16:
    return sign2n;
  };

  cout<<"# Error in XSType in Numeric.\n";
  cout<<"# MT-index "<<i<<" is NOT assigned.\n";
  exit(0);

};

// The following are implemented by Yanagihara-kun for the NFI joint research

void LeastSquaresMethod(int data_num, double *x, double *y, double &a, double &b)
{
  int n=2;
  double sum_xy=0., sum_y=0., sum_x=0., sum_x2=0.; 

  for(int i=0;i<data_num;i++){
    sum_xy+=x[i]*y[i];
    sum_y+=y[i];
    sum_x+=x[i];
    sum_x2+=pow(x[i],2);
  };

  vector< vector<double> > A(n);
  for(int i=0;i<n;i++){
    A[i].resize(n);
  };
  A[0][0]=sum_x2;
  A[0][1]=sum_x;
  A[1][0]=sum_x;
  A[1][1]=data_num;

  //inverse of A
  vector<double> A_inv(n*n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A_inv[i*n+j]=A[i][j];
    };
  };
  gauss_inverse(A_inv, n); // [A_inv] is inversed.

  a=A_inv[0]*sum_xy+A_inv[1]*sum_y;
  b=A_inv[2]*sum_xy+A_inv[3]*sum_y;
};

void LeastSquaresMethodQuadr(int data_num, double *x, double *y, double &a, double &b, double &c)
{
  int n=3;
  double sum_xy=0., sum_x2y=0., sum_y=0., sum_x=0., sum_x2=0., sum_x3=0., sum_x4=0.;

  for(int i=0;i<data_num;i++){
    sum_xy+=x[i]*y[i];
    sum_x2y+=pow(x[i],2)*y[i];
    sum_y+=y[i];
    sum_x+=x[i];
    sum_x2+=pow(x[i],2);
    sum_x3+=pow(x[i],3);
    sum_x4+=pow(x[i],4);
  };

  vector< vector<double> > A(n);
  for(int i=0;i<n;i++){
    A[i].resize(n);
  };
  A[0][0]=sum_x4;
  A[0][1]=sum_x3;
  A[0][2]=sum_x2;
  A[1][0]=sum_x3;
  A[1][1]=sum_x2;
  A[1][2]=sum_x;
  A[2][0]=sum_x2;
  A[2][1]=sum_x;
  A[2][2]=data_num;

  //inverse of A
  vector<double> A_inv(n*n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A_inv[i*n+j]=A[i][j];
    };
  };
  gauss_inverse(A_inv, n); // [A_inv] is inversed.

  a=A_inv[0]*sum_x2y+A_inv[1]*sum_xy+A_inv[2]*sum_y;
  b=A_inv[3]*sum_x2y+A_inv[4]*sum_xy+A_inv[5]*sum_y;
  c=A_inv[6]*sum_x2y+A_inv[7]*sum_xy+A_inv[8]*sum_y;

};


void LeastSquaresMethodCubic(int data_num, double *x, double *y, double &a, double &b, double &c, double &d)
{
  int n=4;
  double sum_xy=0., sum_x2y=0., sum_x3y=0., sum_y=0., sum_x=0., sum_x2=0., sum_x3=0., sum_x4=0., sum_x5=0., sum_x6=0.;

  for(int i=0;i<data_num;i++){
    sum_xy+=x[i]*y[i];
    sum_x2y+=pow(x[i],2)*y[i];
    sum_x3y+=pow(x[i],3)*y[i];
    sum_y+=y[i];
    sum_x+=x[i];
    sum_x2+=pow(x[i],2);
    sum_x3+=pow(x[i],3);
    sum_x4+=pow(x[i],4);
    sum_x5+=pow(x[i],5);
    sum_x6+=pow(x[i],6);
  };

  vector< vector<double> > A(n);
  for(int i=0;i<n;i++){
    A[i].resize(n);
  };
  A[0][0]=sum_x6;
  A[0][1]=sum_x5;
  A[0][2]=sum_x4;
  A[0][3]=sum_x3;
  A[1][0]=sum_x5;
  A[1][1]=sum_x4;
  A[1][2]=sum_x3;  
  A[1][3]=sum_x2;
  A[2][0]=sum_x4;
  A[2][1]=sum_x3;
  A[2][2]=sum_x2;
  A[2][3]=sum_x;
  A[3][0]=sum_x3;
  A[3][1]=sum_x2;
  A[3][2]=sum_x;
  A[3][3]=data_num;

  //inverse of A
  vector<double> A_inv(n*n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A_inv[i*n+j]=A[i][j];
    };
  };
  gauss_inverse(A_inv, n); // [A_inv] is inversed.

  a=A_inv[0]*sum_x3y+A_inv[1]*sum_x2y+A_inv[2]*sum_xy+A_inv[3]*sum_y;
  b=A_inv[4]*sum_x3y+A_inv[5]*sum_x2y+A_inv[6]*sum_xy+A_inv[7]*sum_y;
  c=A_inv[8]*sum_x3y+A_inv[9]*sum_x2y+A_inv[10]*sum_xy+A_inv[11]*sum_y;
  d=A_inv[12]*sum_x3y+A_inv[13]*sum_x2y+A_inv[14]*sum_xy+A_inv[15]*sum_y;

};

void MovingAveraging(vector<real> &data, int points)
{
  int sz=data.size();
  if(points>sz){
    cout<<"# Error in MovingAveraging.\n";
    cout<<"# The number of averaging points is larger than the total number of data.\n";
    exit(0);
  };

  real sum=0.;
  for(int i=0;i<points;i++){
    sum+=data[i];
  };

  for(int i=0;i<sz-points;i++){
    real org=data[i];
    data[i]=sum/points;
    sum+=data[points+i]-org;
  };

  for(int i=0;i<points-1;i++){
    data.pop_back();
  };
};

void Normalize(vector<real> &data, real factor)
{
  int sz=data.size();
  real sum=0.;
  for(int i=0;i<sz;i++){
    sum+=data[i];
  };
  factor/=sum;
  for(int i=0;i<sz;i++){
    data[i]*=factor;
  };
};
