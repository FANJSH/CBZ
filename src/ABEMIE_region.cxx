#include<iostream>
#include "ABEMIE_region.h"

using namespace std;

region::region()
{
  num=0;
  nump=0;
}

void region::add_plane(real x1,real y1,
                       real x2,real y2,int i,int j)
{
  be[num].put_data(x1,y1,x2,y2,i,j);
  num++;
  nump+=j;
  if(num>max_be)cout<<"Cannot add plane into region class\n";
}

void region::cal_gh()
{
  int imax=med->GetImax();
  int icol,jcol;

  //GeomVector *sp =new GeomVector;
  GeomVector sp;
  real b2,bdk;

  for (int i=0;i<imax;i++){
    b2 =med->GetBmat().get_dat(i);
    bdk=med->GetBmatdk().get_dat(i);
    icol=0;
    for (int ib=0;ib<num;ib++){
      for(int ip=0;ip<be[ib].get_pnt();ip++){
	//*sp=be[ib].get_pos(ip);
	sp=be[ib].get_pos(ip);
	jcol=0;
        for (int ib2=0;ib2<num;ib2++){
	  int pnt=be[ib2].get_pnt();
	  real *gr=new real[pnt];
	  real *hr=new real[pnt];
	  real *dgr=new real[pnt];
	  real *dhr=new real[pnt];
    	  //if(ib==ib2){be[ib2].cal_ghii(b2,bdk,*sp,gr,hr,dgr,dhr);}
    	  if(ib==ib2){be[ib2].cal_ghii(b2,bdk,sp,gr,hr,dgr,dhr);}
	  //else{be[ib2].cal_ghij(b2,bdk,*sp,gr,hr,dgr,dhr);};
	  else{be[ib2].cal_ghij(b2,bdk,sp,gr,hr,dgr,dhr);};
	  for(int j=0;j<be[ib2].get_pnt();j++){
	    gmat[i].put_data(icol,jcol,gr[j]);
	    hmat[i].put_data(icol,jcol,hr[j]);
	    dgmat[i].put_data(icol,jcol,dgr[j]);
	    dhmat[i].put_data(icol,jcol,dhr[j]);
	    jcol+=1;
	  };
	  delete [] gr;
          delete [] hr;
          delete [] dgr;
          delete [] dhr;
	};
      icol+=1;
      };
    };
    for(int j=0;j<nump;j++){
      hmat[i].put_data(j,j,hmat[i].get_dat(j,j)+0.5);
    };
  };

  //delete sp;
}

void region::cal_matrix(bool sp3)
{
  int gr=med->GetImax();

  GroupData1D *res=new GroupData1D(nump);
  GroupData1D *fl2=new GroupData1D(nump);
  GroupData2D *gmat2=new GroupData2D(nump,nump);
  GroupData2D *hmat2=new GroupData2D(nump,nump);
  GroupData2D *ginv=new GroupData2D(nump,nump);

  for(int i=0;i<gr;i++){

    gmat2->copy(gmat[i]);
    hmat2->copy(hmat[i]);

    fl2->copy(GetRegionFlx(i));

    int icol=0;
    for(int j=0;j<num;j++){
      int bpnt=be[j].get_pnt();
      int ty=be[j].get_kind();
      if(ty==4){
        for(int l=0;l<bpnt;l++){
	  for(int k=0;k<nump;k++){
	    gmat2->put_data(k,icol+l,-hmat[i].get_dat(k,icol+l));
	    hmat2->put_data(k,icol+l,-gmat[i].get_dat(k,icol+l));
	  };
	  fl2->put_data(icol+l,GetRegionCur(i).get_dat(icol+l));
        };
      };
      if(ty==6){
	if(!sp3){
        // Original (diffusion)
	real inv_dif=1./med->GetMacxs().GetData1d(d).get_dat(i);
	for(int l=0;l<bpnt;l++){
	  real a1=-0.4692*inv_dif;
	  real b1=a1/med->GetCmat().get_dat(i,i)*med->GetCmat().get_dat(i,1-i)*be[j].GetFlx(l).get_dat(1-i)
           -1/med->GetCmat().get_dat(i,i)*med->GetCmat().get_dat(i,1-i)*be[j].GetCur(l).get_dat(1-i);
	  for(int k=0;k<nump;k++){
	    gmat2->put_data(k,icol+l,gmat2->get_dat(k,icol+l)-hmat[i].get_dat(k,icol+l)/a1);
	  };
	  fl2->put_data(icol+l,-b1/a1);
	};
	}else{
	// SP3
	real c11=med->GetCmat().get_dat(0,0);
	real c12=med->GetCmat().get_dat(0,1);
	real c21=med->GetCmat().get_dat(1,0);
	real c22=med->GetCmat().get_dat(1,1);
        for(int l=0;l<bpnt;l++){
	  if(i==0){
  	    // (group 0)
            real a1=0.5*c11-3./8.*c21;
	    real a2=-0.5*c12+3./8.*c22;
            real d0=med->GetMacxs().GetData1d(d).get_dat(0);
            real a3=-d0*c11/a1;
            for(int k=0;k<nump;k++){
	      gmat2->put_data(k,icol+l,gmat2->get_dat(k,icol+l)-hmat[i].get_dat(k,icol+l)*a3);
	    }; 
            real b1=-d0*c12*be[j].GetCur(l).get_dat(1)+a2*be[j].GetFlx(l).get_dat(1);
            fl2->put_data(icol+l,b1/a1);
	  }else{
	    real a1=21./40.*c22-3./40.*c12;
            real a2=-21./40.*c21+3./40.*c11;
            real d1=med->GetMacxs().GetData1d(d).get_dat(1);
            real a3=-d1*c22/a1;
            for(int k=0;k<nump;k++){
	      gmat2->put_data(k,icol+l,gmat2->get_dat(k,icol+l)-hmat[i].get_dat(k,icol+l)*a3);
	    }; 
            real b1=-d1*c21*be[j].GetCur(l).get_dat(0)+a2*be[j].GetFlx(l).get_dat(0);
            fl2->put_data(icol+l,b1/a1);
	  };
	};
      };
      };
      icol+=bpnt;
    };

    *ginv=gmat2->inverse();
    
    *res=*ginv*(*hmat2**fl2);
    //for(int ii=0;ii<nump;ii++){
      //for(int jj=0l;jj<nump;jj++){
      //cout<<ii<<" "<<jj<<" "<<ginv->get_dat(ii,jj)<<"\n";
      //cout<<i<<" "<<ii<<" "<<res->get_dat(ii)<<"\n";
	//};
    //};

    PutCurToPlane(i,*res,sp3);

    dfdf[i]=*ginv**hmat2;

    *res=*ginv*(dhmat[i]*GetRegionFlx(i)-dgmat[i]*GetRegionCur(i));
    PutDfdkToPlane(i,*res);
  };

  delete res;
  delete fl2;
  delete gmat2;
  delete hmat2;
  delete ginv;
}

GroupData1D region::GetRegionFlx(int ig)
{
  GroupData1D ret(nump);
  int icol=0;
  for(int i=0;i<num;i++){
    for(int j=0;j<be[i].get_pnt();j++){
      ret.put_data(icol,be[i].GetFlx(j).get_dat(ig));
      icol++;
    };
  };
  return ret;
};

GroupData1D region::GetRegionCur(int ig)
{
  GroupData1D ret(nump);
  int icol=0;
  for(int i=0;i<num;i++){
    for(int j=0;j<be[i].get_pnt();j++){
      ret.put_data(icol,be[i].GetCur(j).get_dat(ig));
      icol++;
    };
  };
  return ret;
};

void region::PutDfdkToPlane(int ig, GroupData1D data)
{
  int icol=0;
  for(int i=0;i<num;i++){
    int bpnt=be[i].get_pnt();
    for(int j=0;j<bpnt;j++){
      be[i].GetDfdk(j).put_data(ig,data.get_dat(icol));
      icol++;
    };
  };
}

void region::PutCurToPlane(int ig, GroupData1D cur,bool sp3)
{
  int icol=0;
  for(int i=0;i<num;i++){
    for(int j=0;j<be[i].get_pnt();j++){
      int ki=get_be(i).get_kind();
      if(ki!=4&&ki!=6){
        be[i].GetCur(j).put_data(ig,cur.get_dat(icol));
      }else if(ki==4){
        be[i].GetFlx(j).put_data(ig,cur.get_dat(icol));}
      else{
	be[i].GetCur(j).put_data(ig,cur.get_dat(icol));
	// Original (for Diffusion)
	if(!sp3){
        real dif=med->GetMacxs().GetData1d(d).get_dat(ig);
        real a1=-0.4692/dif;
	real b1=-a1/med->GetCmat().get_dat(ig,ig)*
	  med->GetCmat().get_dat(ig,1-ig)*be[i].GetFlx(j).get_dat(1-ig)
	  -1/med->GetCmat().get_dat(ig,ig)*med->GetCmat().get_dat(ig,1-ig)
	  *be[i].GetCur(j).get_dat(1-ig);
	be[i].GetFlx(j).put_data(ig,(cur.get_dat(icol)-b1)/a1);
	}else{
	// SP3
	real c11=med->GetCmat().get_dat(0,0);
	real c12=med->GetCmat().get_dat(0,1);
	real c21=med->GetCmat().get_dat(1,0);
	real c22=med->GetCmat().get_dat(1,1);
        if(ig==0){
	  // (group 1)
          real a1=0.5*c11-3./8.*c21;
          real a2=-0.5*c12+3./8.*c22;
          real d0=med->GetMacxs().GetData1d(d).get_dat(0);
          real b1=-d0*c11*cur.get_dat(icol)-d0*c12*be[i].GetCur(j).get_dat(1)+a2*be[i].GetFlx(j).get_dat(1);
          be[i].GetFlx(j).put_data(ig,b1/a1);
	}else{
          real a1=21./40.*c22-3./40.*c12;
          real a2=-21./40.*c21+3./40.*c11;
          real d1=med->GetMacxs().GetData1d(d).get_dat(1);
          real b1=-d1*c21*be[i].GetCur(j).get_dat(0)-d1*c22*cur.get_dat(icol)+a2*be[i].GetFlx(j).get_dat(0);
          be[i].GetFlx(j).put_data(ig,b1/a1);
	};
	};
      };
      icol++;
    };
  };
}

void region::put_vert()
{
  real retx=0.;
  real rety=0.;
  GeomVector cent,ret;

  for(int i=0;i<num;i++){
    int ic1=0;
    int ic2=0;
    cent=(be[i].get_p1()+be[i].get_p2())*0.5;
    real a1=be[i].get_a1();
    real b1=be[i].get_b1();
    real c1=be[i].get_c1();
    if(b1<0.0){a1=-a1; b1=-b1; c1=-c1;};
    if(b1==0.0&&a1<0.0){a1=-a1; c1=-c1;};
    for(int j=0;j<num;j++){
      int ju=0;
      if(i==j)ju=1;
      real a2=be[j].get_a();
      real b2=be[j].get_b();
      real c2=be[j].get_c();
      if(a1==0.0&&a2==0.0)ju=1;
      if(b1==0.0&&b2==0.0)ju=1;
      if(b1!=0.0&&b2!=0.0&&fabs(a1/b1-a2/b2)<1e-5)ju=1;
      if(ju==0){
	real x,y;
	if(a1==0.0){y=-c1/b1; x=-b2/a2*y-c2/a2;}
        else if(a2==0.0){y=-c2/b2; x=-b1/a1*y-c1/a1;}
        else{
          y=(c2/a2-c1/a1)/(b1/a1-b2/a2);
          x=-b1/a1*y-c1/a1;
	};
        real x1=be[j].get_p1().getx();
        real x2=be[j].get_p2().getx();
        real y1=be[j].get_p1().gety();
        real y2=be[j].get_p2().gety();
        if((x1-x)*(x2-x)<=0.0&&(y1-y)*(y2-y)<=0.0){
          if(b1!=0.0){
	    if(cent.getx()<x){ic1++;}
	    else{ic2++;};
          }
          else{
            if(cent.gety()<y){ic1++;}
            else{ic2++;};
          };
          if(b1!=0.0){
            if(ic1%2==0){retx=b1; rety=-a1;}
            else{retx=-b1; rety=a1;};
          }
          else{
            if(ic1%2==0){retx=0; rety=a1;}
            else{retx=0; rety=-a1;};
          };
	};
      };
    };
    ret.put_data(retx,rety);
    be[i].put_pver(ret);
  };
};

int region::get_address(int ip)
{
  int ret=0;
  if(ip!=0){
    for(int i=0;i<ip;i++){
      ret+=be[i].get_pnt();
    };
  };
  return ret;
};

real region::cal_real_flx()
{
  int imax=med->GetImax();
  GroupData1D tmp(imax);
  GroupData1D err(imax);
  real errmax=0;

  for(int i=0;i<num;i++){
    int kind=get_be(i).get_kind();
    for(int j=0;j<get_be(i).get_pnt();j++){
      tmp=get_med()->GetCmat()*get_be(i).GetFlx(j);
      err=(tmp-RealFlux[i][j])/RealFlux[i][j];
      RealFlux[i][j]=tmp;
      if(kind==1||kind==2||kind==5){
        for(int k=0;k<imax;k++){
	  if(fabs(err.get_dat(k))>errmax)errmax=fabs(err.get_dat(k));
        };
      };
    };
  };

  return errmax;
};

void region::change_pnt(int ii)
{
  for(int i=0;i<num;i++){
    be[i].change_pnt(ii);
  };
  nump*=ii;
}

void region::copy_gh(region sec)
{
  int imax=med->GetImax();
  for(int i=0;i<imax;i++){
    gmat[i].copy(sec.gmat[i]);
    hmat[i].copy(sec.hmat[i]);
    dgmat[i].copy(sec.dgmat[i]);
    dhmat[i].copy(sec.dhmat[i]);
  };
}

void region::set_array()
{
  int imax=be[0].get_imax();
  gmat.resize(imax);
  hmat.resize(imax);
  dgmat.resize(imax);
  dhmat.resize(imax);
  dfdf.resize(imax);

  for(int i=0;i<imax;i++){
    gmat[i].put_yx(nump,nump);
    hmat[i].put_yx(nump,nump);
    dgmat[i].put_yx(nump,nump);
    dhmat[i].put_yx(nump,nump);
    dfdf[i].put_yx(nump,nump);
  };

  RealFlux.resize(num);
  for(int i=0;i<num;i++){
    int pnt=get_be(i).get_pnt();
    RealFlux[i].resize(pnt);
    for(int j=0;j<pnt;j++){
      RealFlux[i][j].put_imax(imax);
    };
  };
}

real region::get_integrated_flux(int g)
{
  int imax=med->GetImax();
  vector<real> crt(imax);

  for(int ig=0;ig<imax;ig++){
    crt[ig]=0.0;
    for(int j=0;j<num;j++){
      real dst=be[j].get_long()*0.5;
      int pnt=be[j].get_pnt();
      switch(pnt){
      case 1:
        crt[ig]+=dst*2*be[j].GetCur(0).get_dat(ig);
	break;
      case 2:
	crt[ig]+=dst*be[j].GetCur(0).get_dat(ig);
	crt[ig]+=dst*be[j].GetCur(1).get_dat(ig);
	break;
      case 3:
        crt[ig]+=dst/(3*htp*htp)*be[j].GetCur(0).get_dat(ig);
        crt[ig]+=dst*(-2/(3*htp*htp)+2)*be[j].GetCur(1).get_dat(ig);
        crt[ig]+=dst/(3*htp*htp)*be[j].GetCur(2).get_dat(ig);
        break;
      case 4:
	real ab=dst*(ht2*ht2-htp4*htp4);
        crt[ig]+=ab/htp4*(-htp4/3+htp4*ht2*ht2)*be[j].GetCur(0).get_dat(ig);
	crt[ig]-=ab/ht2*(-ht2/3+htp4*htp4*ht2)*be[j].GetCur(1).get_dat(ig);
	crt[ig]+=ab/ht2*(ht2/3-htp4*htp4*ht2)*be[j].GetCur(2).get_dat(ig);
        crt[ig]-=ab/htp4*(htp4/3-htp4*ht2*ht2)*be[j].GetCur(3).get_dat(ig);
        break;
      };
    };
  };

  real ret=0.;
  for(int k=0;k<imax;k++){
    ret-=med->GetCmat().get_dat(g,k)/med->GetBmat().get_dat(k)*crt[k];
  };
  return ret;
};

// ****************
// * MakeRegion   *
// ****************

region MakeRegion::tri(real x,real y,real d,int m,int j)
{
  region reg;
  if(j==0)j=pnt;

  for(int i=0;i<m;i++){
    reg.add_plane(x,y,x+d/m,y,imax,j);
    x+=d/m;
  };
  for(int i=0;i<m;i++){
    reg.add_plane(x,y,x,y+d/m,imax,j);
    y+=d/m;
  };
  for(int i=0;i<m;i++){
    reg.add_plane(x,y,x-d/m,y-d/m,imax,j);
    y-=d/m;
    x-=d/m;
  };
  reg.put_vol(d*d*0.5);

  return reg;
}

region MakeRegion::square(real x,real y,real dx,real dy,
                          int mx, int my, int j)
{
  region reg;
  if(j==0)j=pnt;

  for(int i=0;i<mx;i++){
    reg.add_plane(x,y,x+dx/mx,y,imax,j);
    x+=dx/mx;
  };
  for(int i=0;i<my;i++){
    reg.add_plane(x,y,x,y+dy/my,imax,j);
    y+=dy/my;
  };
  for(int i=0;i<mx;i++){
    reg.add_plane(x,y,x-dx/mx,y,imax,j);
    x-=dx/mx;
  };
  for(int i=0;i<my;i++){
    reg.add_plane(x,y,x,y-dy/my,imax,j);
    y-=dy/my;
  };

  reg.put_vol(dx*dy);

  return reg;
}

region MakeRegion::drawing(real xs,real ys,int pv,vector< vector<real> > vec)
{
  region ret;
  real x=xs;
  real y=ys;
  for(int i=0;i<pv;i++){
    ret.add_plane(x,y,x+vec[i][0],y+vec[i][1],imax,pnt);
    x+=vec[i][0];
    y+=vec[i][1];
  };
  if(x!=xs||y!=ys)cout<<"Not closed in region in drawing of MakeRegion.cxx\n";
  return ret;
}
