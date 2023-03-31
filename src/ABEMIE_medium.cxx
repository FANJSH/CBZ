#include <cstdlib>
#include <iostream>
#include "ABEMIE_medium.h"

using namespace std;

void BemMedium::PutImaxBem(int i)
{
  bmat.put_imax(i);
  bmatdk.put_imax(i);
  amat.put_yx(i,i);
  cmat.put_yx(i,i);
  cmati.put_yx(i,i);
  cmatdk.put_yx(i,i);
  dcmat.put_yx(i,i);
}

void BemMedium::PutMedium(Medium inp)
{
  PutImaxBem(inp.GetImax());
  PutImax(inp.GetImax());
  PutPL(inp.GetPL());

  Macxs.GetData1d(siga).copy(inp.GetMacxs().GetData1d(siga));
  Macxs.GetData1d(nusigf).copy(inp.GetMacxs().GetData1d(nusigf));
  Macxs.GetData1d(chi).copy(inp.GetMacxs().GetData1d(chi));
  Macxs.GetData1d(d).copy(inp.GetMacxs().GetData1d(d));
  Macxs.GetData1d(sigt).copy(inp.GetMacxs().GetData1d(sigt));
  Macxs.GetData1d(sigt,1).copy(inp.GetMacxs().GetData1d(sigt,1));
  Macxs.GetData1d(sigtr).copy(inp.GetMacxs().GetData1d(sigtr));
  flux[0].copy(inp.GetFlux(0));
  flux[1].copy(inp.GetFlux(1));
  for(int i=0;i<pl+1;i++){
    Macxs.GetData2d(sigs,i).copy(inp.GetMacxs().GetData2d(sigs,i));
  };
}

GroupData2D BemMedium::cal_amat(real keff)
{
  GroupData1D a1(imax);
  GroupData1D a2(imax);
  GroupData1D a3(imax);
  GroupData2D ret(imax,imax);

  a1=Macxs.GetData1d(chi)/Macxs.GetData1d(d);
  a2=Macxs.GetData1d(nusigf)/keff;
  a3=(Macxs.GetData2d(sigs).get_sumx()
     -Macxs.GetData2d(sigs).get_diag()
     +Macxs.GetData1d(siga))/Macxs.GetData1d(d);

  for (int i=0;i<imax;i++){
    for(int j=0;j<imax;j++){
      real x=a1.get_dat(i)*a2.get_dat(j);
      if(i==j){
	x-=a3.get_dat(i);
      };
      //if(j<i){
      if(j!=i){
	x+=Macxs.GetData2d(sigs).get_dat(j,i)
          /Macxs.GetData1d(d).get_dat(i);
      };
      ret.put_data(i,j,x);
    };
  };

  return ret;
};

GroupData2D BemMedium::cal_amat_SP3_1grp(real keff)
{
  GroupData2D ret(imax,imax);

  real d1=Macxs.GetData1d(d).get_dat(0);
  real st=Macxs.GetData1d(sigt).get_dat(0);
  real st0=Macxs.GetData1d(sigt).get_dat(0)-Macxs.GetData2d(sigs,0).get_dat(0);
  real d3=0.33333333333333333/st;
  real nsf=Macxs.GetData1d(nusigf).get_dat(0);
  
  // Diffusion coefficient for phi_2
  Macxs.GetData1d(d).put_data(1,27./35.*d3);

  ret.put_data(0,0,-st0/d1+1./d1/keff*nsf);
  ret.put_data(0,1,-2./d1/keff*nsf+2./d1*st0);
  ret.put_data(1,0,-35./27./d3*(-0.4*st0+0.4/keff*nsf));
  ret.put_data(1,1,-35./27./d3*(st+0.8*st0-0.8/keff*nsf));

  return ret;
};

GroupData1D BemMedium::cal_bmat(GroupData2D am)
{
  GroupData1D ret(imax);
  am.CalEigenvaluesByLeverrieFaddeev(ret);
  return ret;
  /*
  GroupData2D mat(imax,imax);
  GroupData2D unit(imax,imax);
  GroupData1D ret(imax);
  real *pval=new real[imax];

  for(int i=0;i<imax;i++){
    unit.put_unit();
    if(i==0){mat=am;}
    else{mat=am*mat;}
    pval[i]=mat.get_diag().get_sum()/(i+1);
    unit=unit*pval[i];
    mat=mat-unit;
  };

  if(imax==2){
    real x=sqrt(pval[0]*pval[0]+4*pval[1]);
    ret.put_data(0,(pval[0]+x)*0.5);
    ret.put_data(1,(pval[0]-x)*0.5);
  }else if(imax==1){
    ret.put_data(0,pval[0]);
  }else{
    cout<<"More than 3 group is not tested.\n";
    cout<<"B value calculations\n";
    vector<real> sol;
    bool conv=false;

    real vinp=1e-9;
    while(vinp<100.){
      real v0=vinp;
      for(int it=0;it<100;it++){
        real f=0.;
        for(int i=0;i<=imax;i++){
          real coef=1.;
  	  if(i!=imax)coef=-pval[imax-1-i];
	  f+=pow(v0,i)*coef;
        };
        real df=0.;
        for(int i=0;i<imax;i++){
	  real coef=1.;
	  if(i!=imax-1)coef=-pval[imax-1-i-1];
	  df+=pow(v0,i)*coef*(i+1);
        };
        real dv=f/df;
        v0-=dv;
        if(fabs(f)<1e-10){
          bool eql=false;
          for(int i=0;i<int(sol.size());i++){
	    if(fabs(sol[i]/v0-1.)<1e-5)eql=true;
	  };
	  if(!eql){
	    cout<<v0<<"\n";
            sol.push_back(v0);
	    if(int(sol.size())==imax)conv=true;
	  };
          break;
        };
      };
      if(vinp>0){
	vinp*=-1;
      }else{
	vinp*=-1;
	vinp+=1e-5;
      };
    };
  };

  delete [] pval;
  return ret;
  */
}

GroupData2D BemMedium::cal_cmat(GroupData2D am,GroupData1D bm)
{
  GroupData2D ret(imax,imax);
  if(imax==2){
    real c21=1.;
    real c22=1.;
    ret.put_data(1,0,c21);
    ret.put_data(1,1,c22);
    ret.put_data(0,0,c21*(bm.get_dat(0)-am.get_dat(1,1))/am.get_dat(1,0));
    ret.put_data(0,1,c22*(bm.get_dat(1)-am.get_dat(1,1))/am.get_dat(1,0));
  }else if(imax==1){
    ret.put_data(0,0,1.);
    //}else if(imax==3){
  }else{
    cout<<"# Error in BemMedium::cal_cmat.\n";
    cout<<"# Not coded yet.\n";
    exit(0);
  };
  return ret;
}

void BemMedium::put_abc(real keff,real dk,bool sp3)
{
  if(!sp3){
    amat=cal_amat(keff);
  }else{
    amat=cal_amat_SP3_1grp(keff);
  };
  bmat=cal_bmat(amat);
  cmat=cal_cmat(amat,bmat);
  cmati=cmat.inverse();

  for(int i=0;i<imax;i++){
    for(int j=0;j<imax;j++){
      dcmat.put_data(i,j,Macxs.GetData1d(d).get_dat(i)*cmat.get_dat(i,j));
    };
  };

  if(dk>0){
    GroupData2D amat2(imax,imax);
    GroupData1D bmat2(imax);
    GroupData2D cmat2(imax,imax);
    if(!sp3){
      amat2=cal_amat(keff+dk);
    }else{
      amat2=cal_amat_SP3_1grp(keff+dk);
    };
    bmat2=cal_bmat(amat2);
    cmat2=cal_cmat(amat2,bmat2);

    bmatdk=bmat2.sqrt1d()-bmat.sqrt1d();
    bmatdk=bmatdk*(1.0/dk);

    cmatdk=cmat2-cmat;
    cmatdk=cmatdk*(1.0/dk);
  };
}
