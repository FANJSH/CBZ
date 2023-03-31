#include <cstdlib>
#include "UNC_Covariance.h"

using namespace std;

void Covariance::PutType(string type)
{
  if(type=="variance"){Variance=true;}
  else{Variance=false;};
  size1=-1;
  size2=-1;
  if(Variance){
    val.resize(1);
  }else{
    val.resize(2);
  };
};

void Covariance::PutSize(int ginp)
{
  size1=ginp;
  size2=ginp;
  val[0].put_imax(size1);
  if(!Variance){
    val[1].put_imax(size2);
  };
  cov.put_yx(size1,size2);
};

void Covariance::PutSize(int ginp1, int ginp2)
{
  if(Variance){
    cout<<"This instance of Covariance is variance.\n";
    cout<<"You cannot set two ginp.\n";
    exit(0);
  };
  size1=ginp1;
  size2=ginp2;
  val[0].put_imax(size1);
  val[1].put_imax(size2);
  cov.put_yx(size1,size2);
};

void Covariance::CheckPosition(int i,int j)
{
  if(j==-1)j=i;
  if(i<0||i>=size1){
    cout<<"You chose larger value than covariance size.\n";
    cout<<"The value is "<<i<<"\n";
    exit(0);
  };
  if(!Variance&&(j<0||j>=size2)){
    cout<<"You chose larger value than covariance size.\n";
    cout<<"The value is "<<j<<"\n";
    exit(0);
  };
};

void Covariance::PutVal(real *inp)
{
  if(!Variance){
    cout<<"This instance of Covariance is covariance.\n";
    cout<<"You should use PutVal(i,inp2) method.\n";
    exit(0);
  };
  if(size1==-1){
    cout<<"Error in PutVal of Covariance.\n";
    cout<<"You should do PutSize before doing PutVal.\n";
    exit(0);
  };
  val[0].put_data(inp);
};

void Covariance::PutVal(int i,real inp)
{
  CheckPosition(i);
  val[0].put_data(i,inp);
};

void Covariance::PutVal(int i,GroupData1D &inp)
{
  CheckPosition(i);
  int tmp=size1;
  if(i==1)tmp=size2;
  if(tmp==-1){
    cout<<"Error in PutVal of Covariance.\n";
    cout<<"You should do PutSize before doing PutVal.\n";
    exit(0);
  };
  for(int j=0;j<tmp;j++){  
    val[i].put_data(j,inp.get_dat(j));
  };
  /*
  for(int i=0;i<tmp;i++){  
    val[0].put_data(i,inp.get_dat(i));
  };
  */
};

void Covariance::PutVal(int i,real *inp)
{
  if(i<0||i>=2){
    cout<<"# Error in PutVal of Covariance.\n";
    cout<<"# i should be 0 or 1.\n";
    exit(0);
  };
  if(Variance){
    cout<<"# This instance of Covariance is variance.\n";
    cout<<"# You should use PutVal(inp) method.\n";
    exit(0);
  };
  int tmp;
  if(i==0){tmp=size1;}else{tmp=size2;};
  if(tmp==-1){
    cout<<"# Error in PutVal of Covariance.\n";
    cout<<"# You should do PutSize before doing PutVal.\n";
    exit(0);
  };
  val[i].put_data(inp);
};

void Covariance::PutCov(real *inp, string type)
{
  if(size1==-1||size2==-1){
    cout<<"# Error in PutCov of Covariance.\n";
    cout<<"# You should do PutSize before doing PutCov.\n";
    exit(0);
  };
  if(type=="Absolute"||type=="absolute"){
    cov.put_data(inp);
  }else if(type=="Relative"||type=="relative"){
    int ind=0;
    for(int i=0;i<size1;i++){
      real tmp1=val[0].get_dat(i);
      for(int j=0;j<size2;j++){
	real tmp2;
	if(Variance){tmp2=val[0].get_dat(j);}
	else{tmp2=val[1].get_dat(j);};
	inp[ind]*=(tmp1*tmp2);
	ind++;
      };
    };
    cov.put_data(inp);
  }else{
    cout<<"Error in PutCov in Covariance.\n";
    cout<<"You have to put covariance data as absolute or relative value.\n";
    cout<<"You requested "<<type<<"\n";
    exit(0);
  };

  /*
  if(Variance){
    for(int i=0;i<size1;i++){
      if(cov.get_dat(i,i)==0.){
	for(int j=0;j<size1;j++){
	  cov.put_data(i,j,0.);
	  cov.put_data(j,i,0.);
	};
      };
    };
  };
  */
};

void Covariance::PutCov(real *dev, real *corr, string type)
{
  if(!Variance){
    cout<<"Error in PutCov of Covariance.\n";
    cout<<"You cannot do PutCov because this instance is covariance.\n";
    exit(0);
  };
  if(size1==-1||size2==-1){
    cout<<"Error in PutCov of Covariance.\n";
    cout<<"You should do PutSize before doing PutCov.\n";
    exit(0);
  };
  real *inp=new real[size1*size2];
  int ind=0;
  for(int i=0;i<size1;i++){
    for(int j=0;j<size2;j++){
      if(corr[ind]<-1.||corr[ind]>1.){
	cout<<"Error in PutCov in Covariance.\n";
	cout<<"Correlation is not from -1 to 1.\n";
	exit(0);
      };
      inp[ind]=corr[ind]*dev[i]*dev[j];
      ind++;
    };
  };
  PutCov(inp,type);
  delete [] inp;
};

void Covariance::PutCov(real *dev, GroupData2D corr, string type)
{
  if(corr.get_x()!=size1||corr.get_y()!=size1){
    cout<<"Error in PutCov in Covariance.\n";
    cout<<"There is inconsistency between correlation matrix and covariance.\n";
    exit(0);
  };

  real *cor=new real[size1*size1];
  for(int i=0;i<size1*size1;i++){
    cor[i]=corr.get_dat(i);
  };

  PutCov(dev,cor,type);
  delete [] cor;
};

void Covariance::PutCov(GroupData2D covin, string type)
{
  int size2tmp=size1;
  if(!Variance)size2tmp=size2;

  if(covin.get_x()!=size1||covin.get_y()!=size2tmp){
    cout<<"Error in PutCov in Covariance.\n";
    cout<<"There is inconsistency in inputted covariance.\n";
    exit(0);
  };

  for(int i=0;i<size1;i++){
    for(int j=0;j<size2tmp;j++){
      real tmp=covin.get_dat(i,j);
      if(type=="Relative"||type=="relative"){
        tmp=covin.get_dat(i,j)*val[0].get_dat(i);
        if(Variance){
  	  tmp*=val[0].get_dat(j);
        }else{
	  tmp*=val[1].get_dat(j);
        };
      };
      cov.put_data(i,j,tmp);
    };
  };
};

void Covariance::PutStandardDeviation(int i,real inp,string type)
{
  CheckPosition(i);
  if(!Variance){
    cout<<"PutStandardDeviation(i,inp,type) cannot be used for not-variance.\n";
    exit(0);
  };
  if(type=="Absolute"){
    cov.put_data(i,i,inp*inp);
  }else{
    real tmp=inp*val[0].get_dat(i);
    cov.put_data(i,i,tmp*tmp);
  };
};

GroupData2D Covariance::GetCovariance(string type)
{
  GroupData2D ret(size1,size2);
  ret.copy(cov);
  if(type=="relative"||type=="Relative"){
    for(int i=0;i<size1;i++){
      real tmp1=val[0].get_dat(i);
      for(int j=0;j<size2;j++){
	real tmp2;
	if(Variance){tmp2=val[0].get_dat(j);}
	else{tmp2=val[1].get_dat(j);};
	real org=cov.get_dat(i,j);
	if(tmp1!=0.&&tmp2!=0.){
  	  ret.put_data(i,j,org/tmp1/tmp2);
	}else{
  	  ret.put_data(i,j,0.);
	};
      };
    };
  };
  return ret;
};

GroupData1D Covariance::GetStandardDeviation(string type)
{
  if(!Variance){
    cout<<"This instance of Covariance is covariance.\n";
    cout<<"You cannot define standard deviation.\n";
    exit(0);
  };

  GroupData1D ret(size1);
  ret=cov.get_diag();
  for(int i=0;i<size1;i++){
    real tmp=sqrt(ret.get_dat(i));
    if(type=="Relative"||type=="relative")tmp/=fabs(val[0].get_dat(i));
    ret.put_data(i,tmp);
  };
  return ret;
};

GroupData2D Covariance::GetCorrelationMatrix()
{
  if(!Variance){
    cout<<"This instance of Covariance is covariance.\n";
    cout<<"You cannot define correlation matrix.\n";
    exit(0);
  };

  GroupData1D stddev(size1);
  stddev=GetStandardDeviation("Absolute");

  GroupData2D corr(size1,size1);
  for(int i=0;i<size1;i++){
    real tmp1=stddev.get_dat(i);
    for(int j=0;j<size2;j++){
      real tmp2=stddev.get_dat(j);
      real tmp=cov.get_dat(i,j);
      if(tmp1!=0.&&tmp2!=0.){
        tmp/=(tmp1*tmp2);
      }else{
	tmp=0.;
      };
      corr.put_data(i,j,tmp);
    };
  };

  return corr;
};

void Covariance::Factorize(real v)
{
  cov=cov*v;
};

void Covariance::Factorize(vector<real>& v)
{
  if(!Variance){
    cout<<"# Error in Covariance::Factorize.\n";
    cout<<"# This method cannot be used for non-variance matrix.\n";
    exit(0);
  };

  int sz=v.size();
  if(sz!=size1){
    cout<<"# Error in Covariance::Factorize.\n";
    cout<<"# 1D vector size is inconsistent with the size of covariance matrix.\n";
    exit(0);
  };

  for(int i=0;i<size1;i++){
    for(int j=i;j<size1;j++){
      real org=cov.get_dat(i,j);
      real mod=org*v[i]*v[j];
      cov.put_data(i,j,mod);
      cov.put_data(j,i,mod);
    };
  };
};

void Covariance::Factorize(vector<real>& v1, vector<real>& v2)
{
  int sz1=v1.size();
  int sz2=v2.size();
  if(sz1!=size1||sz2!=size2){
    cout<<"# Error in Covariance::Factorize.\n";
    cout<<"# 1D vector size is inconsistent with the size of covariance matrix.\n";
    exit(0);
  };

  for(int i=0;i<size1;i++){
    for(int j=0;j<size2;j++){
      real org=cov.get_dat(i,j);
      real mod=org*v1[i]*v2[j];
      cov.put_data(i,j,mod);
    };
  };
};


void Covariance::SetNoCorrelation()
{
  if(!Variance){
    cout<<"# Warning in Covariance::SetNoCorrelation.\n";
    cout<<"# Not coded for covariance data.\n";
  };
  for(int i=0;i<size1;i++){
    for(int j=0;j<size1;j++){
      if(i!=j)cov.put_data(i,j,0.);
    };
  };
};

Covariance Covariance::Cond(GroupData1D &wgt,int bgrp,int *bndgrp)
{
  if(!Variance){
    cout<<"# Error in Covariance::Cond.\n";
    cout<<"# Not coded for covariance data.\n";
    exit(0);
  };

  int grp=size1;
  if(bndgrp[bgrp-1]!=grp-1){
    cout<<"# Error in Covariance::Cond.\n";
    cout<<"# The lowest energy group is INCONSISTENT.\n";
    exit(0);
  };

  // 
  GroupData1D bval(bgrp);
  bval=val[0].Cond(wgt,bgrp,bndgrp);
  
  GroupData1D bwgt(bgrp);
  bwgt=wgt.CondSum(bgrp,bndgrp);

  Covariance bcov("variance");
  bcov.PutSize(bgrp);
  bcov.PutVal(0,bval);

  GroupData2D covorg=cov;

  GroupData2D covdat(bgrp,bgrp);
  for(int i=0;i<bgrp;i++){
    real wgti_inv=1./bwgt.get_dat(i);
    for(int j=i;j<bgrp;j++){
      real sum=0.;
      real wgtj_inv=1./bwgt.get_dat(j);
      int is=0;
      if(i!=0)is=bndgrp[i-1]+1;
      int ie=bndgrp[i];
      for(int i2=is;i2<=ie;i2++){
	int js=0;
	if(j!=0)js=bndgrp[j-1]+1; 
	int je=bndgrp[j];
	for(int j2=js;j2<=je;j2++){
	  sum+=wgt.get_dat(i2)*wgti_inv*wgt.get_dat(j2)*wgtj_inv*covorg.get_dat(i2,j2);
	};
      };
      covdat.put_data(i,j,sum);
    };
  };

  for(int i=0;i<bgrp;i++){
    for(int j=i+1;j<bgrp;j++){
      covdat.put_data(j,i,covdat.get_dat(i,j));
    };
  };

  bcov.PutCov(covdat,"Absolute");

  return bcov;
};

void Covariance::NormalizeWeightFunction(GroupData1D &wgt,int bgrp,int *bndgrp)
{
  // In prior to do group collapsing for covariance matrix,
  // weight functions are normalized by this method
  // because usually weight function w_g/w_G is used in the group collapsing.

  GroupData1D bwgt(bgrp);
  bwgt=wgt.CondSum(bgrp,bndgrp);

  for(int i=0;i<bgrp;i++){
    real wgti_inv=1./bwgt.get_dat(i);
    int is=0;
    if(i!=0)is=bndgrp[i-1]+1;
    int ie=bndgrp[i];
    for(int i2=is;i2<=ie;i2++){
      real org=wgt.get_dat(i2);
      wgt.put_data(i2,org*wgti_inv);
    };
  };
};

Covariance Covariance::CondRelative(GroupData1D &wgt,int bgrp,int *bndgrp)
{
  //
  // Only relative covariance data has physical meaning
  // because "val" is not well collapised
  //

  if(!Variance){
    cout<<"# Error in Covariance::Cond.\n";
    cout<<"# Not coded for covariance data.\n";
    exit(0);
  };

  int grp=size1;
  if(bndgrp[bgrp-1]!=grp-1){
    cout<<"# Error in Covariance::Cond.\n";
    cout<<"# The lowest energy group is INCONSISTENT.\n";
    exit(0);
  };

  // 
  GroupData1D bval(bgrp);
  bval=val[0].Cond(wgt,bgrp,bndgrp);
  
  GroupData1D bwgt(bgrp);
  bwgt=wgt.CondSum(bgrp,bndgrp);

  Covariance bcov("variance");
  bcov.PutSize(bgrp);
  bcov.PutVal(0,bval);

  GroupData2D covorg=GetCovariance("Relative");

  GroupData2D covdat(bgrp,bgrp);
  for(int i=0;i<bgrp;i++){
    for(int j=i;j<bgrp;j++){
      real sum=0.;
      int is=0;
      if(i!=0)is=bndgrp[i-1]+1;
      int ie=bndgrp[i];
      for(int i2=is;i2<=ie;i2++){
	int js=0;
	if(j!=0)js=bndgrp[j-1]+1; 
	int je=bndgrp[j];
	for(int j2=js;j2<=je;j2++){
	  sum+=wgt.get_dat(i2)*wgt.get_dat(j2)*covorg.get_dat(i2,j2);
	};
      };
      covdat.put_data(i,j,sum);
    };
  };

  for(int i=0;i<bgrp;i++){
    for(int j=i+1;j<bgrp;j++){
      covdat.put_data(j,i,covdat.get_dat(i,j));
    };
  };

  bcov.PutCov(covdat,"Relative");

  return bcov;
};

Covariance Covariance::GetCovarianceLognormal(real factor)
{
  if(!Variance){
    cout<<"# Error in Covariance::GetCovarianceLognormal.\n";
    cout<<"# This method can be used for variance matrix.\n";
    exit(0);
  };

  Covariance cov_ln("variance");
  cov_ln.PutSize(size1);

  GroupData2D covinp(size1,size1);

  for(int i=0;i<size1;i++){
    real tmp1,tmp2;
    GetMeanAndVarianceForLognormal(val[0].get_dat(i), cov.get_dat(i,i), tmp1,tmp2);
    cov_ln.PutVal(i,tmp1);
    covinp.put_data(i,i,tmp2);
  };

  for(int i=0;i<size1;i++){
    for(int j=i+1;j<size1;j++){

      // ... Correction to covariance matrix in log(p)
      //real tmp=log(1.+cov.get_dat(i,j)/(val[0].get_dat(i)*val[0].get_dat(j)));
      //tmp*=factor;

      // ... Correction to covariance matrix in p
      real tmp=log(1.+cov.get_dat(i,j)*factor/(val[0].get_dat(i)*val[0].get_dat(j)));

      // range check [-1,+1]
      real corr=tmp/sqrt(covinp.get_dat(i,i)*covinp.get_dat(j,j));
      if(corr>1){
	cout<<"# Error in Covariance::GetCovarianceLognormal.\n";
	cout<<"# Correlation larger than +1 ("<<corr<<") is calculated.\n";
	exit(0);
      };
      if(corr<-1){
	cout<<"# Error in Covariance::GetCovarianceLognormal.\n";
	cout<<"# Correlation smaller than -1 ("<<corr<<") is calculated.\n";
	exit(0);
      };
      covinp.put_data(i,j,tmp);
      covinp.put_data(j,i,tmp);
    };
  };

  cov_ln.PutCov(covinp);

  return cov_ln;
};

void Covariance::DoRandomSampling(int num, vector<GroupData1D> &sample, string type)
{
  if(type!="Absolute"&&type!="Relative"){
    cout<<"# Error in Covariance::DoRandomSampling.\n";
    cout<<"# Keyword [type] should be [Absolute] or [Relative].\n";
    cout<<"# You choose "<<type<<"\n";
    exit(0);
  };

  if(!Variance){
    cout<<"# Warning in Covariance::DoRandomSampling.\n";
    cout<<"# This method does NOT work because this covariance is NOT variance.\n";
    return;
  };

  GroupData2D covdat=GetCovariance(type);
  int n=covdat.get_x();

  // +++ Matrix diagonalization
  GroupData2D diagmat;
  GroupData2D evecmat;
  covdat.Diagonalization(diagmat,evecmat);
  // C=U^{-1}PU
  //   P:covdat, C:diagmat, U:evecmat

  vector<real> eigen(n);
  //int icnt=0;
  for(int i=0;i<n;i++){
    real tmp=diagmat.get_dat(i*n+i);
    eigen[i]=tmp;
    if(eigen[i]<0.)eigen[i]=0.;
    //if(tmp>1e-5)icnt++;
  };
  
  // +++ Random sampling of cross section
  sample.resize(num);
  for(int ii=0;ii<num;ii++){  
    sample[ii].put_imax(n);
    vector<real> eigen_sample(n);
    for(int i=0;i<n;i++){
      eigen_sample[i]=ransu_gauss(0.,sqrt(eigen[i]));
    };
    vector<real> dxs(n);
    for(int i=0;i<n;i++){
      real tmp=0.;
      for(int j=0;j<n;j++){
        tmp+=evecmat.get_dat(i*n+j)*eigen_sample[j];
      };
      sample[ii].put_data(i,tmp);
    };
  };
};




