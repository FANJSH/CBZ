#include<cstdlib>
#include "ABEMIE_core.h"

using namespace std;

core::~core(){
};

void core::Init()
{
  numreg=0;
  numib=0;
  numibp=0;
  count=0;
  nummed=0;
  sp3=false;
};

void core::AddRegion(region i)
{
  reg[numreg]=i;
  reg[numreg].put_vert();
  numreg++;
  if(numreg==maxreg){
    cout<<"Please increase maxreg in ABEMIE_const.h.\n";
    exit(0);
  };
}

void core::AddMedium(Medium i)
{
  BemMedium inp;
  inp.PutMedium(i);
  med[nummed]=inp;
  nummed++;
};

void core::put_boundary_info()
{
  for(int i=0;i<numreg;i++){
    for(int j=0;j<reg[i].get_num();j++){
      int is=reg[i].get_plane_kind(j);
      if(is==0){
	reg[i].set_plane_kind(j,1); // treated as 'inner boundary'
	ibr[numib]=i; //region ID
	ibp[numib]=j; //plane ID
	for(int i2=0;i2<numreg;i2++){
	  if(i!=i2){
            for(int j2=0;j2<reg[i2].get_num();j2++){
 	      if(reg[i].get_be(j).BeIdentical(reg[i2].get_be(j2),0.1)){
                if(reg[i].get_be(j).get_p1().BeIdentical(reg[i2].get_be(j2).get_p1(),0.1))
                  {reg[i2].set_plane_kind(j2,2);} // treated as 'opposite side of Inner Boundary'
		else{reg[i2].set_plane_kind(j2,5);}; // treated as 'opposite direction & side of IB'
	        ibrr[numib]=i2; //region ID(opposite side)
	        ibrp[numib]=j2; //plane ID (opposite side)
	      };
	    };
	  };
	};
	numibp+=get_ibe(numib).get_pnt();
	numib++;
	if(numib==maxib){
	  cout<<"Please increase maxib in ABEMIE_const.h.\n";
	  exit(0);
	};
      };	
    };
  };
  set_array_uppercal();
}

void core::show_boundary_info()
{
  cout<<" (Region:Plane:Kind) ";
  for(int i=0;i<numreg;i++){
    for(int j=0;j<reg[i].get_num();j++){
      int is=reg[i].get_be(j).get_kind();
      cout <<i<<" : "<<j<<" : "<<is<<"\n";
    };
  };
}

void core::show_innerb_info()
{
  cout<<" Total Inner Boundary = "<<numib<<"\n"; 
  cout<<"i:Reg:Plane:(Reg):(Plane) \n";
  for(int i=0;i<numib;i++){
    cout<<i<<":"<<ibr[i]<<":"<<ibp[i]<<":"<<ibrr[i]<<":"<<
      ibrp[i]<<"\n";
  };
}

void core::put_flx_innerr()
{
  GroupData2D cc(imax,imax);
  for(int i=0;i<numib;i++){
    cc=get_cc(ibrr[i],ibr[i]);
    for(int j=0;j<get_pnt(i);j++){
      get_iber(i).GetFlx(j).copy
	(cc*get_ibe(i).GetFlx(get_iber(i).opposite_p(j)));
    };
  };
}

void core::put_outerb(real x1,real y1,real x2,real y2,
                      BCondition b,GroupData1D f)
{
  int kind=0;
  if(b==Zeroflux)   kind=3;
  if(b==Reflective) kind=4;
  if(b==Vacuum)     kind=6;
  if(kind==0){
    cout<<"Inappropriate boundary condition.\n";
    exit(0);
  };

  plane ob(x1,y1,x2,y2,1,1);

  for(int i=0;i<numreg;i++){
    for(int j=0;j<reg[i].get_num();j++){
      bool ir=reg[i].get_be(j).BePartial(ob,1e-4);
      int bpnt=reg[i].get_be(j).get_pnt();
      if(ir&&kind==3){
        reg[i].set_plane_kind(j,3); // treated as Outer Boundary(flux is fixed)
        for(int k=0;k<bpnt;k++){
          reg[i].get_be(j).GetFlx(k).copy(f);
        };
      };
      if(ir&&kind==4){
        reg[i].set_plane_kind(j,4); // treated as OB(current is fixed)
        for(int k=0;k<bpnt;k++){
          reg[i].get_be(j).GetCur(k).copy(f);
	};
      };
      if(ir&&kind==6){
        reg[i].set_plane_kind(j,6); // treated as OB(Extrapolation)
	for(int k=0;k<bpnt;k++){
	  reg[i].get_be(j).GetCur(k).copy(f);
	  reg[i].get_be(j).GetFlx(k).copy(f);
	};
      };
    };
  };
}

void core::cal_f()
{
  GroupData1D dcc(imax);
  real y=0.0;
  int index=0;
  for(int i=0;i<numib;i++){
    int ir=ibr[i];
    int irr=ibrr[i];
    for(int i2=0;i2<get_pnt(i);i2++){
      int i2r=get_iber(i).opposite_p(i2);
      dcc=reg[ir].get_med()->GetDcmat()*get_ibe(i).GetCur(i2)
         +reg[irr].get_med()->GetDcmat()*get_iber(i).GetCur(i2r);
      dcc=dcc*-1.0;
      f.put_data(index,dcc);
      index+=imax;
      y+=get_ibe(i).GetFlx(i2).get_sum();
    };
  };
  f.put_data(numibp*imax,-1.0+y);
}

void core::modf_flx(bool accel)
{
  GroupData1D dat(imax);
  real omega=1.;
  if(accel){ 
    omega=1./(1.-lamda);
    cout<<"  (Aitkin accelaration) omega= "<<omega<<"\n";
  };
  int icol=0;
  int jcol=0;
  real rrsum=0.;
  for(int i=0;i<numib;i++){
    for(int i2=0;i2<get_pnt(i);i2++){
      for(int j=0;j<imax;j++){
	real oflx=get_ibe(i).GetFlx(i2).get_dat(j);
	real dflx=df.get_dat(icol);
	real rr=dflx/res.get_dat(icol);
	if(rr>0.&&rr<1.){
          rrsum+=rr;
	  jcol++;
	};
	res.put_data(icol,dflx);
        dat.put_data(j,oflx+dflx*omega);
	icol++;
      };
      get_ibe(i).GetFlx(i2).copy(dat);
    };
  };

  if(jcol!=0){
    lamda=rrsum/jcol;
  }else{
    lamda=0.;
  };

  df.put_data(numibp*imax,df.get_dat(numibp*imax)*omega);
}
	  
void core::get_power()
{
  real *pow=new real[numreg];
  real *crt=new real[imax];
  real dst;
  real powtot=0;
  real voltot=0;

  //FILE *power=fopen("power","w");

  for(int i=0;i<numreg;i++){
    pow[i]=0.0;
    for(int ig=0;ig<imax;ig++){
      crt[ig]=0.0;
      for(int j=0;j<reg[i].get_num();j++){
	dst=reg[i].get_be(j).get_long()*0.5;
	int pnt=reg[i].get_be(j).get_pnt();
	switch(pnt){
	case 1:
	  crt[ig]+=dst*2*reg[i].get_be(j).GetCur(0).get_dat(ig);
	  break;
	case 2:
	  crt[ig]+=dst*reg[i].get_be(j).GetCur(0).get_dat(ig);
	  crt[ig]+=dst*reg[i].get_be(j).GetCur(1).get_dat(ig);
	  break;
	case 3:
          crt[ig]+=dst/(3*htp*htp)
	    *reg[i].get_be(j).GetCur(0).get_dat(ig);
          crt[ig]+=dst*(-2/(3*htp*htp)+2)
	    *reg[i].get_be(j).GetCur(1).get_dat(ig);
          crt[ig]+=dst/(3*htp*htp)
	    *reg[i].get_be(j).GetCur(2).get_dat(ig);
	  break;
	case 4:
	  real ab=dst*(ht2*ht2-htp4*htp4);
          crt[ig]+=ab/htp4*(-htp4/3+htp4*ht2*ht2)*
	    reg[i].get_be(j).GetCur(0).get_dat(ig);
	  crt[ig]-=ab/ht2*(-ht2/3+htp4*htp4*ht2)*
	    reg[i].get_be(j).GetCur(1).get_dat(ig);
	  crt[ig]+=ab/ht2*(ht2/3-htp4*htp4*ht2)*
	    reg[i].get_be(j).GetCur(2).get_dat(ig);
	  crt[ig]-=ab/htp4*(htp4/3-htp4*ht2*ht2)*
	    reg[i].get_be(j).GetCur(3).get_dat(ig);
	  break;
	}
      };
    };

    for(int j=0;j<imax;j++){
      for(int k=0;k<imax;k++){
	pow[i]-=reg[i].get_med()->GetMacxs().GetData1d(nusigf).get_dat(j)*
	  reg[i].get_med()->GetCmat().get_dat(j,k)/
	  reg[i].get_med()->GetBmat().get_dat(k)*crt[k];
      };
    };
    powtot+=pow[i];
    if(pow[i]>0)voltot+=reg[i].get_vol();
  };

  real sta=powtot/voltot;

  for(int i=0;i<numreg;i++){
    pow[i]/=sta;
    pow[i]/=reg[i].get_vol();
    if(pow[i]>0){
      cout<<i<<":"<<pow[i]<<"\n";
      //fprintf(power,"%6.4f \n",pow[i]);
    };
  };
  delete [] pow;
  delete [] crt;

  //fclose(power);
}

real core::get_integrated_flux(int r,int g)
{
  return reg[r].get_integrated_flux(g);
}

void core::cal_df_p0()
{
  df.set_zero();

  GroupData2D cc1(imax,imax);
  GroupData2D cc2(imax,imax);

  for(int i=0;i<numib;i++){
    int ir=ibr[i];
    int irr=ibrr[i];
    int pnt1=get_pnt(i);
    int ac1=reg[ir].get_address(ibp[i]);
    int ad1=reg[irr].get_address(ibrp[i]);
    for(int i2=0;i2<numib;i2++){
      int ir2=ibr[i2];
      int irr2=ibrr[i2];
      int pnt2=get_pnt(i2);
      int ac2=reg[ir2].get_address(ibp[i2]);
      int ad2=reg[irr2].get_address(ibrp[i2]);
      int index=0;
      if(ir==irr2) cc1=get_cc(ir,ir2);
      if(irr==irr2)cc2=get_cc(irr,ir2);
      real *x=new real[imax*pnt1*imax*pnt2];
      for(int k1=0;k1<pnt1;k1++){
	int is1=k1;
	int isr1=get_iber(i).opposite_p(k1);
        for(int j=0;j<imax;j++){
	  for(int k2=0;k2<pnt2;k2++){
	    int is2=k2;
	    int isr2=get_iber(i2).opposite_p(k2);
       	    for(int j2=0;j2<imax;j2++){
  	      real xx=0.0;
              if(ir==ir2){
		xx+=reg[ir].get_med()->GetDcmat().get_dat(j,j2)*
	            reg[ir].get_dfdf(j2).get_dat(ac1+is1,ac2+is2);
	      }else if(ir==irr2){
	        xx+=reg[ir].get_med()->GetDcmat().get_dat(j,j2)*
		    cc1.get_dat(j2,j2)*
	            reg[ir].get_dfdf(j2).get_dat(ac1+is1,ad2+isr2);
              };
              if(irr==ir2){
	        xx+=reg[irr].get_med()->GetDcmat().get_dat(j,j2)*
	            reg[irr].get_dfdf(j2).get_dat(ad1+isr1,ac2+is2);
	      }else if(irr==irr2){
	        xx+=reg[irr].get_med()->GetDcmat().get_dat(j,j2)*
		    cc2.get_dat(j2,j2)*
	            reg[irr].get_dfdf(j2).get_dat(ad1+isr1,ad2+isr2);
	      };
	      x[index]=xx;
	      index++;
	    };
	  };
	};
      };
      if(i==i2){
        jself[i].put_data(x);
	jself[i]=jself[i].inverse();
      };
      int index2=jinfo[i*numib+i2];
      if(index2!=0){
	jother[index2-1].put_data(x);
      };
      delete [] x;
    };
  };
};

void core::cal_df_direct_inversion()
{
  int mat_size=numibp*imax+1;
  real *a=new real[mat_size*mat_size];

  GroupData2D cc1(imax,imax);
  GroupData2D cc2(imax,imax);

  int index_y=0;
  for(int i=0;i<numib;i++){
    int ir=ibr[i];
    int irr=ibrr[i];
    int pnt1=get_pnt(i);
    int ac1=reg[ir].get_address(ibp[i]);
    int ad1=reg[irr].get_address(ibrp[i]);
    for(int k1=0;k1<pnt1;k1++){
      int is1=k1;
      int isr1=get_iber(i).opposite_p(k1);
      for(int j=0;j<imax;j++){
	int index_x=0;
        for(int i2=0;i2<numib;i2++){
          int ir2=ibr[i2];
          int irr2=ibrr[i2];
          int pnt2=get_pnt(i2);
          int ac2=reg[ir2].get_address(ibp[i2]);
          int ad2=reg[irr2].get_address(ibrp[i2]);
          if(ir==irr2) cc1=get_cc(ir,ir2);
          if(irr==irr2)cc2=get_cc(irr,ir2);
	  for(int k2=0;k2<pnt2;k2++){
	    int is2=k2;
	    int isr2=get_iber(i2).opposite_p(k2);
       	    for(int j2=0;j2<imax;j2++){
  	      real xx=0.0;
              if(ir==ir2){
		xx+=reg[ir].get_med()->GetDcmat().get_dat(j,j2)*
	            reg[ir].get_dfdf(j2).get_dat(ac1+is1,ac2+is2);
	      }else if(ir==irr2){
	        xx+=reg[ir].get_med()->GetDcmat().get_dat(j,j2)*
		    cc1.get_dat(j2,j2)*
	            reg[ir].get_dfdf(j2).get_dat(ac1+is1,ad2+isr2);
              };
              if(irr==ir2){
	        xx+=reg[irr].get_med()->GetDcmat().get_dat(j,j2)*
	            reg[irr].get_dfdf(j2).get_dat(ad1+isr1,ac2+is2);
	      }else if(irr==irr2){
	        xx+=reg[irr].get_med()->GetDcmat().get_dat(j,j2)*
		    cc2.get_dat(j2,j2)*
	            reg[irr].get_dfdf(j2).get_dat(ad1+isr1,ad2+isr2);
	      };
	      a[index_y*mat_size+index_x]=xx;
	      index_x++;
	    };
	  };
	};
        index_y++;
      };
    };
  };

  int index=0;
  int tmp1=numibp*imax;
  for(int i=0;i<numib;i++){
    int pnt1=get_pnt(i);
    for(int i2=0;i2<pnt1;i2++){
      GroupData1D tmp(imax);
      tmp=cal_jfk(i,i2,get_iber(i).opposite_p(i2));
      for(int j=0;j<imax;j++){
        a[index*mat_size+tmp1]=tmp.get_dat(j);
	a[tmp1*mat_size+index]=-1.;
	index++;
      };
    };
  };


  real *b=new real[mat_size];
  for(int i=0;i<mat_size;i++){
    b[i]=f.get_dat(i);
  };


  // Solving ax=b

  // +++ with fast algolythm
  int n=mat_size;
  real p;
  real det=1.0;
  for(int id=0;id<n-1;id++){
    det*=a[id*n+id];
    for(int i=id+1;i<n;i++){
      if(ledge[i]<=id){
        p=a[i*n+id]/a[id*n+id];
        for (int j=id;j<=redge[i];j++){
          a[i*n+j]-=p*a[id*n+j];
        };
        a[i*n+(n-1)]-=p*a[id*n+(n-1)];
        b[i]-=p*b[id];
      };
    };
  };
  det*=a[n*n-1];

  b[n-1]/=a[n*n-1];
  for (int ii=2;ii<n+1;ii++){
    int i=n+1-ii-1;
    for (int k=i+1;k<n+1;k++){
      b[i]-=a[i*n+k]*b[k];
    };
    b[i]/=a[i*n+i];
  };
  // +++ convenional gauss 
  //gauss(a,b,mat_size,1);
  // +++ conjugate gradient
  //bicg(a,b,mat_size);

  for(int i=0;i<mat_size;i++){
    df.put_data(i,b[i]);
  };

  delete [] a;
  delete [] b;
};

void core::set_edgedata_for_direct_inversion()
{
  int n=numibp*imax+1;
  ledge.resize(n,0);
  redge.resize(n,0);

  vector<bool> a(n*n);
  for(int i=0;i<n*n;i++){a[i]=false;};

  int index_y=0;
  for(int i=0;i<numib;i++){
    int ir=ibr[i];
    int irr=ibrr[i];
    int pnt1=get_pnt(i);
    for(int k1=0;k1<pnt1;k1++){
      for(int j=0;j<imax;j++){
	int index_x=0;
        for(int i2=0;i2<numib;i2++){
          int ir2=ibr[i2];
          int irr2=ibrr[i2];
          int pnt2=get_pnt(i2);
	  for(int k2=0;k2<pnt2;k2++){
       	    for(int j2=0;j2<imax;j2++){
              if(ir==ir2||ir==irr2||irr==ir2||irr==irr2){
  	        a[index_y*n+index_x]=true;
	      };
	      index_x++;
	    };
	  };
	};
        index_y++;
      };
    };
  };

  for(int i=0;i<n;i++){
    int pos=i*n;
    int ll=0;
    for(int j=0;j<i;j++){
      if(a[pos+j]){
	ll=j;
	break;
      };
    };
    ledge[i]=ll;

    int rr=n-2;
    for(int j=n-2;j>=i-1;j--){
      if(a[pos+j]){
	rr=j;
	break;
      };
    };
    redge[i]=rr;
  };

  for(int i=1;i<n;i++){
    if(redge[i]<redge[i-1])redge[i]=redge[i-1];
  };

};

void core::cal_df_p(real errf,real errk,int iter)
{
  GroupData1D *ff=new GroupData1D[numib];
  GroupData1D *dp1=new GroupData1D[numib];

  df.set_zero();

  int index2=0;
  for(int i=0;i<numib;i++){
    int pnt1=get_pnt(i);
    ff[i].put_imax(imax*pnt1);
    dp1[i].put_imax(imax*pnt1);
    int index=0;
    for(int i2=0;i2<pnt1;i2++){
      jk[i].put_data(index,cal_jfk(i,i2,get_iber(i).opposite_p(i2)));
      ff[i].put_data(index,f,index2,imax);
      index+=imax;
      index2+=imax;
    };
  };

  int ind1=0;
  for(int i=0;i<numib;i++){
    ff[i]=jself[i]*ff[i];
    jk[i]=jself[i]*jk[i];
    for(int j=0;j<numib;j++){
      int index=jinfo[ind1];
      ind1++;
      if(index!=0){
	jother[index-1]=jself[i]*jother[index-1];
      };
    };
  };

  GroupData1D temp;
  for(int i=0;i<numib;i++){
    dp1[i]=ff[i];
  };
  for(int it2=0;it2<10;it2++){
    int ind=0;
    for(int i=0;i<numib;i++){
      temp.put_imax(get_pnt(i)*imax);
      temp.set_zero();
      for(int j=0;j<numib;j++){
        int index=jinfo[ind];
        ind++;
        if(index!=0){
          temp=temp-jother[index-1]*dp1[j];
        };
      };
      dp1[i]=temp+ff[i];
    };
  };

  if(errk>1e-4||iter==0){
    for(int i=0;i<numib;i++){
      dp2[i]=jk[i];
    };
    for(int it2=0;it2<10;it2++){
      int ind=0;
      for(int i=0;i<numib;i++){
        temp.put_imax(get_pnt(i)*imax);
	temp.set_zero();
        for(int j=0;j<numib;j++){
          int index=jinfo[ind];
          ind++;
          if(index!=0){
            temp=temp-jother[index-1]*dp2[j];
          };
        };
        dp2[i]=temp+jk[i];
      };
    };
  };

  real psum=0.0;
  real psum0=0.0;
  real dk0=0.001;
  for(int i=0;i<numib;i++){
    psum0+=dp1[i].get_sum();
    psum-=dp2[i].get_sum();
  };
  psum*=dk0;
  psum+=psum0;
  real dk=dk0-dk0/(-psum+psum0)*(-psum-f.get_dat(numibp*imax));

  int index=0;
  for(int i=0;i<numib;i++){
    df.put_data(index,dp1[i]-dp2[i]*dk);
    index+=get_pnt(i)*imax;
  };

  psum0=psum;
  psum=0.;
  for(int i=0;i<numib;i++){
    psum+=dp1[i].get_sum();
    psum-=dp2[i].get_sum()*dk;
  };
  real dk00=dk0;
  dk0=dk;
  dk=dk0-(dk0-dk00)/(-psum+psum0)*(-psum-f.get_dat(numibp*imax));

  df.put_data(numibp*imax,dk);

  delete [] ff;
  delete [] dp1;
}

void core::put_jinfo()
{
  for(int i=0;i<numreg;i++){
    for(int j=0;j<reg[i].get_num();j++){
      if(reg[i].get_plane_kind(j)==0)cout<<"Error in put_jinfo in core.\n";
    };
  };

  jinfo.resize(numib*numib,0);

  for(int i=0;i<numib;i++){
    int ir=ibr[i];
    int irr=ibrr[i];
    for(int i2=0;i2<numib;i2++){
      int ir2=ibr[i2];
      int irr2=ibrr[i2];
      int judge=0;
      if(ir==ir2||ir==irr2||irr==ir2||irr==irr2)judge=1;
      if(i!=i2&&judge==1){
        jinfo[i*numib+i2]=count+1;
        count++;
      }else{
        jinfo[i*numib+i2]=0;
      };
    };
  };

  jother.resize(count+1);

  for(int i=0;i<numib;i++){
    int pp=get_pnt(i)*imax;
    jself[i].put_yx(pp,pp);
    jk[i].put_imax(pp);
    dp2[i].put_imax(pp);
    for(int i2=0;i2<numib;i2++){
      int in2=jinfo[i*numib+i2];
      if(in2!=0){
	jother[in2-1].put_yx(pp,imax*get_pnt(i2));
      };
    };
  };
}

GroupData2D core::get_cc(int i,int j)
{
  return reg[i].get_med()->GetCmati()*reg[j].get_med()->GetCmat();
}

GroupData1D core::cal_jfk(int i,int is1,int is2)
{
  int ir=ibr[i];
  int irr=ibrr[i];
  return reg[ir].get_med()->GetMacxs().GetData1d(d).mult
    (reg[ir].get_med()->GetCmatdk()*get_ibe(i).GetCur(is1))
    +reg[ir].get_med()->GetDcmat()*get_ibe(i).GetDfdk(is1)
          +reg[irr].get_med()->GetMacxs().GetData1d(d).mult
    (reg[irr].get_med()->GetCmatdk()*get_iber(i).GetCur(is2)
     +reg[irr].get_med()->GetDcmat()*get_iber(i).GetDfdk(is2));
}

real core::cal_real_flx()
{
  real errmax=0.0;
  for(int i=0;i<numreg;i++){
    real err=reg[i].cal_real_flx();
    if(err>errmax){
      errmax=err;
    };
  };
  return errmax;
}

void core::out_flx(int ig)
{
  GeomVector sp;
  //FILE *out2=fopen("out2","w");
  for(int i=0;i<numreg;i++){
    for(int j=0;j<reg[i].get_num();j++){
      for(int k=0;k<reg[i].get_be(j).get_pnt();k++){
        sp=reg[i].get_be(j).get_pos(k);
	//real flx=reg[i].get_be(j).get_data1da(k,ig,"real_flx");
	//flx=reg[i].get_be(j).get_data1da(k,ig,"flx");
	//fprintf(out2,"%12.5f %12.5f %10.3e \n",sp.getx(),sp.gety(),flx);
      };
    };
  };
  //fclose(out2);
}

void core::change_pnt(int ii)
{
  for(int i=0;i<numreg;i++){
    reg[i].change_pnt(ii);
  };
  numibp*=ii;
  set_array_uppercal();
  put_jinfo();
}

void core::put_rmed(int *inp)
{
  for(int i=0;i<numreg;i++){rmed[i]=inp[i];};
}

void core::put_eql(int *inp)
{
  for(int i=0;i<numreg;i++){eql[i]=inp[i];};
}

void core::initial_flx()
{
  real *inp=new real[imax];
  real tmp=1./(numibp*imax);
  for(int i=0;i<imax;i++){
    inp[i]=tmp;
  };
  for(int i=0;i<numib;i++){
    for(int ii=0;ii<get_pnt(i);ii++){
      get_ibe(i).GetFlx(ii).put_data(inp);
    };
  };
  delete [] inp;
}
  
void core::set_array_uppercal()
{
  f.put_imax(imax*numibp+1);
  df.put_imax(imax*numibp+1);
  res.put_imax(imax*numibp);

  //jother.resize(numib*12);
  jself.resize(numib);
  jk.resize(numib);
  dp2.resize(numib);

  for(int i=0;i<numreg;i++){
    reg[i].set_array();
  };
}

void core::run(real epsk,real epsf,bool accel,bool direct)
{
  real time_start=clock();
  real time_l=0;
  real time_u=0;

  //FILE *out1=fopen("out1","w");

  put_boundary_info();
  put_jinfo();

  //show_innerb_info();
  //show_boundary_info();

  initial_flx();
  put_keff(1.0);

  for(int i=0;i<numreg;i++){
    reg[i].put_BemMedium(&med[rmed[i]]);
  };

  real t1,t2,t3;

  //*** Iterate ***
  if(direct)set_edgedata_for_direct_inversion();

  keffold=keff;
  real errf=1.0;
  real errk=1.0;

  for(int iter=0;iter<10000;iter++){

    if(iter!=0)keffold=keff;
    t1=clock();

    real dk=keff*1e-3;
    for(int i=0;i<nummed;i++){
      med[i].put_abc(keff,dk,sp3);
    };

    put_flx_innerr();

    for(int i=0;i<numreg;i++){
      if(eql[i]==-1){reg[i].cal_gh();}
      else{reg[i].copy_gh(reg[eql[i]]);};
      //cout<<"+++ Region : "<<i<<"+++\n";
      reg[i].cal_matrix(sp3);
    };

    errf=cal_real_flx();
    if(iter==0)errf=1e+10;
 
    t2=clock();
    cal_f();
    if(direct){
      cal_df_direct_inversion();
    }else{
      cal_df_p0();
      cal_df_p(errf,errk,iter);
    };

    bool accel_inp=false;
    if(accel&&iter!=0&&iter%5==0&&iter>4)accel_inp=true;
    modf_flx(accel_inp);

    keffold=keff;
    keff+=df.get_dat(numibp*imax);

    t3=clock();

    errf=cal_real_flx();
    if(iter==0)errf=1e+10;
    errk=fabs((keff-keffold)/keffold);
    if(iter==0)errk=1e+10;
    time_l+=t2-t1;
    time_u+=t3-t2;

    cout.width(3);
    cout<<iter<<":   ";
    cout.setf(ios::showpoint);
    cout.precision(7);
    cout<<keff<<"  ";
    cout.unsetf(ios::showpoint);
    cout.setf(ios::scientific);
    cout.precision(3);
    cout<<errk<<"  "<<errf<<"\n";
    cout.unsetf(ios::scientific);
    //fprintf(out1,"%d %12.5f %10.3e %10.3e\n",iter,keff,errk,errf);

    if(errk<epsk&&errf<epsf)break;
  };

  //get_power();
  //test_sys.out_flx(0);

  real time_end=clock();
  real time_run=time_end-time_start;
  time_run/=1e+6;
  time_l/=1e+6;
  time_u/=1e+6;
  cout<<time_run<<":"<<time_l<<":"<<time_u<<"\n";

  //fclose(out1);
}

void core::PutMeshXYFromCartCore(CartCore &cinp,real z,int pnt)
{
  MakeRegion ms(imax,pnt);

  int xr=cinp.GetXr();
  int yr=cinp.GetYr();

  int *mat=new int[xr*yr];
  cinp.GetMaterialMapXY(z,mat);

  real xw,yw;
  real xp;
  real yp=0.;
  for(int i=0;i<yr;i++){
    yw=cinp.GetYwid(i);
    xp=0.;
    for(int j=0;j<xr;j++){
      xw=cinp.GetXwid(j);
      AddRegion(ms.square(xp,yp,xw,yw,1,1));
      xp+=xw;
    };
    yp+=yw;
  };

  put_rmed(mat);

  int *eql=new int[xr*yr];
  real *xdef=new real[xr*yr];
  real *ydef=new real[xr*yr];
  int *meddef=new int[xr*yr];
  int *regdef=new int[xr*yr];
  int defid=0;
  for(int i=0;i<yr;i++){
    for(int j=0;j<xr;j++){
      real xtmp=cinp.GetXwid(j);
      real ytmp=cinp.GetYwid(i);
      int medtmp=mat[i*xr+j];
      int sfl=-1;
      for(int k=0;k<defid;k++){
	if(sfl==-1&&xdef[k]==xtmp&&ydef[k]==ytmp&&meddef[k]==medtmp)sfl=k;
      };
      if(sfl==-1){
	xdef[defid]=xtmp;
	ydef[defid]=ytmp;
	meddef[defid]=medtmp;
	regdef[defid]=i*xr+j;
	eql[i*xr+j]=-1;
	defid++;
      }else{
	eql[i*xr+j]=regdef[sfl];
      };
    };
  };

  put_eql(eql);

  GroupData1D inp(imax);
  inp.set_zero();

  real xleng=0.;
  real yleng=0.;
  for(int i=0;i<xr;i++){xleng+=cinp.GetXwid(i);};
  for(int j=0;j<yr;j++){yleng+=cinp.GetYwid(j);};
  if(cinp.GetLeftBC()==0){put_outerb(0,0,0,yleng,Zeroflux,inp);};
  if(cinp.GetLeftBC()==1){put_outerb(0,0,0,yleng,Reflective,inp);};
  if(cinp.GetLeftBC()==2){put_outerb(0,0,0,yleng,Vacuum,inp);};
  if(cinp.GetRightBC()==0){put_outerb(xleng,0,xleng,yleng,Zeroflux,inp);};
  if(cinp.GetRightBC()==1){put_outerb(xleng,0,xleng,yleng,Reflective,inp);};
  if(cinp.GetRightBC()==2){put_outerb(xleng,0,xleng,yleng,Vacuum,inp);};
  if(cinp.GetBackBC()==0){put_outerb(0,0,xleng,0,Zeroflux,inp);};
  if(cinp.GetBackBC()==1){put_outerb(0,0,xleng,0,Reflective,inp);};
  if(cinp.GetBackBC()==2){put_outerb(0,0,xleng,0,Vacuum,inp);};
  if(cinp.GetFrontBC()==0){put_outerb(0,yleng,xleng,yleng,Zeroflux,inp);};
  if(cinp.GetFrontBC()==1){put_outerb(0,yleng,xleng,yleng,Reflective,inp);};
  if(cinp.GetFrontBC()==2){put_outerb(0,yleng,xleng,yleng,Vacuum,inp);};

  delete [] mat;
  delete [] eql;
  delete [] xdef;
  delete [] ydef;
  delete [] meddef;
  delete [] regdef;
}
  
void core::PutIrregularGeometryInformation(IrregularGeometryInformation &inp,int *reg_med,int pnt)
{
  int r=inp.GetIpol();
  vector<GeomPlane> ptmp;
  GeomPolygon pol;
  region reg;
  int *rmed=new int[r];

  for(int i=0;i<r;i++){
    pol=inp.GetPolygon(i);
    rmed[i]=reg_med[pol.GetRegionID()];
    int ipl=pol.GetNumpl();
    ptmp.resize(ipl);
    reg.Init();
    real x1=pol.GetPlane(0).get_p1().getx();
    real y1=pol.GetPlane(0).get_p1().gety();
    for(int i=0;i<ipl;i++){
      real x2=pol.GetPlane(i).get_p2().getx();
      real y2=pol.GetPlane(i).get_p2().gety();
      reg.add_plane(x1,y1,x2,y2,imax,pnt);
      x1=x2;
      y1=y2;
    };
    AddRegion(reg);
  };

  put_rmed(rmed);
  delete [] rmed;

  int *eql=new int[r];
  for(int i=0;i<r;i++){eql[i]=-1;};
  put_eql(eql);
  delete [] eql;
  /*
  GroupData1D ini(imax);
  ini.set_zero();
  put_outerb(  0,  0,100,  0,Zeroflux,ini);
  put_outerb(100,  0,100,100,Zeroflux,ini);
  put_outerb(100,100,  0,100,Zeroflux,ini);
  put_outerb(  0,100,  0,  0,Zeroflux,ini);
  */
}
