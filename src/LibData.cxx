#include <cstdlib>
#include "LibData.h"

// LibDataPTable

LibDataPTable::LibDataPTable()
{
};

void LibDataPTable::PutGroup(int g)
{
  step.resize(g,0);
  probability.resize(g);
  total.resize(g);
  elastic.resize(g);
  fission.resize(g);
  capture.resize(g);
  sigs_gg_p.resize(g);
  downscat_distribution.resize(g);
  maxerr_total.resize(g);
  maxerr_elastic.resize(g);
  maxerr_fission.resize(g);
  maxerr_capture.resize(g);
};

void LibDataPTable::PutProbability(int g,int st,real *p)
{
  step[g]=st;
  probability[g].resize(st,0.);
  total[g].resize(st,0.);
  elastic[g].resize(st,0.);
  fission[g].resize(st,0.);
  capture[g].resize(st,0.);
  for(int i=0;i<st;i++){
    probability[g][i]=p[i];
  };
  sigs_gg_p[g].resize(st);
  for(int i=0;i<st;i++){
    sigs_gg_p[g][i].resize(st+1);
  };
  downscat_distribution[g].resize(st);
};

void LibDataPTable::PutTotal(int g,real *xs)
{
  int st=step[g];
  if(st==0){
    cout<<"# Error in LibDataPTable::PutTotal.\n";
    cout<<"# You have to put probability before putting xs.\n";
    exit(0);
  };
  for(int i=0;i<st;i++){
    total[g][i]=xs[i];
  };
};

void LibDataPTable::PutElastic(int g,real *xs)
{
  int st=step[g];
  if(st==0){
    cout<<"# Error in LibDataPTable::PutElastic.\n";
    cout<<"# You have to put probability before putting xs.\n";
    exit(0);
  };
  for(int i=0;i<st;i++){
    elastic[g][i]=xs[i];
  };
};

void LibDataPTable::PutFission(int g,real *xs)
{
  int st=step[g];
  if(st==0){
    cout<<"# Error in LibDataPTable::PutFission.\n";
    cout<<"# You have to put probability before putting xs.\n";
    exit(0);
  };
  for(int i=0;i<st;i++){
    fission[g][i]=xs[i];
  };
};

void LibDataPTable::PutCapture(int g,real *xs)
{
  int st=step[g];
  if(st==0){
    cout<<"# Error in LibDataPTable::PutCapture.\n";
    cout<<"# You have to put probability before putting xs.\n";
    exit(0);
  };
  for(int i=0;i<st;i++){
    capture[g][i]=xs[i];
  };
};

void LibDataPTable::ReadFile(string mdir,string ss,real maxerr)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  real vtmp;
  int ingrp=1;

  real me1;
  real me2;
  real me3;
  real me4;
  while(ingrp!=-1){
    fin>>ingrp;
    if(ingrp!=-1){
      ingrp--;
      fin>>me1;
      fin>>me2;
      fin>>me3;
      fin>>me4;
      bool accept=false;
      if(me1<maxerr&&
         me2<maxerr&&
         me3<maxerr&&
         me4<maxerr)accept=true;
      if(accept){
	maxerr_total[ingrp]=me1;
	maxerr_elastic[ingrp]=me2;
	maxerr_capture[ingrp]=me3;
	maxerr_fission[ingrp]=me4;
      };
      int st;
      fin>>st;
      real *pp=new real[st];
      for(int i=0;i<st;i++){
	fin>>vtmp;
	pp[i]=vtmp;
      };
      if(accept)PutProbability(ingrp,st,pp);
      for(int i=0;i<st;i++){
	fin>>vtmp;
	pp[i]=vtmp;
      };
      if(accept)PutTotal(ingrp,pp);
      for(int i=0;i<st;i++){
	fin>>vtmp;
	pp[i]=vtmp;
      };
      if(accept)PutElastic(ingrp,pp);
      for(int i=0;i<st;i++){
	fin>>vtmp;
	pp[i]=vtmp;
      };
      if(accept)PutCapture(ingrp,pp);
      for(int i=0;i<st;i++){
	fin>>vtmp;
	pp[i]=vtmp;
      };
      if(accept)PutFission(ingrp,pp);
      delete [] pp;      
    };
  };

  fin.close();
};

void LibDataPTable::WriteFile(string mdir,string ss,int stgrp)
{
  mdir.append(ss);
  ofstream fout;
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  fout.setf(ios::scientific);
  fout.precision(6);
  int sz=step.size();
  for(int i=stgrp;i<sz;i++){
    if(step[i]!=0){
      fout<<i+1<<"\n";
      fout<<maxerr_total[i]<<"\n";
      fout<<maxerr_elastic[i]<<"\n";
      fout<<maxerr_capture[i]<<"\n";
      fout<<maxerr_fission[i]<<"\n";
      fout<<step[i]<<"\n";
      for(int j=0;j<step[i];j++){
	fout<<probability[i][j]<<"\n";
      };
      for(int j=0;j<step[i];j++){
	fout<<total[i][j]<<"\n";
      };
      for(int j=0;j<step[i];j++){
	fout<<elastic[i][j]<<"\n";
      };
      for(int j=0;j<step[i];j++){
	fout<<capture[i][j]<<"\n";
      };
      for(int j=0;j<step[i];j++){
	fout<<fission[i][j]<<"\n";
      };
    };
  };
  fout<<"   -1\n";

  fout.close();
};

void LibDataPTable::ReadMatrixFile(string mdir,string ss)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"(file name) "<<mdir<<"\n";
    exit(0);
  };

  real vtmp;
  int ingrp=1;

  while(ingrp!=0){
    fin>>ingrp;
    if(ingrp!=0){
      ingrp--;
      int st;
      fin>>st;
      if(st!=step[ingrp]){
	cout<<"Error in LibDataPTable::ReadMatrixFile\n";
        cout<<"There is inconsistency between PT and PTMAT in the number of steps.\n";
	exit(0);
      };
      for(int i=0;i<st;i++){
	for(int j=0;j<st+1;j++){
  	  fin>>vtmp;
	  sigs_gg_p[ingrp][i][j]=vtmp;
	};
      };
      for(int i=0;i<st;i++){
	fin>>vtmp;
	downscat_distribution[ingrp][i]=vtmp;
      };
    };
  };

  fin.close();
};

void LibDataPTable::show_self()
{
  int sz=step.size();
  if(sz==0)return;

  cout.setf(ios::scientific);
  cout.precision(3);
  for(int i=0;i<sz;i++){
    if(step[i]!=0){
      cout<<"###################################################\n";
      cout<<"# Group : "<<i<<"\n#\n";
      cout<<"#   Number of divisions : "<<step[i]<<"\n#\n";
      cout<<"#    P      Total     Elastic   Capture   Fission\n";
      for(int j=0;j<step[i];j++){
        cout<<"# "<<probability[i][j];
	cout<<" "<<total[i][j];
	cout<<" "<<elastic[i][j];
	cout<<" "<<capture[i][j];
	cout<<" "<<fission[i][j];
	cout<<"\n";
      };
      cout<<"# (Max.Err. in %)\n";
      cout<<"#           ";
      cout<<maxerr_total[i]<<" ";
      cout<<maxerr_elastic[i]<<" ";
      cout<<maxerr_capture[i]<<" ";
      cout<<maxerr_fission[i]<<" ";
      cout<<"\n";
    };
  };
};

// FTable

void FTable::Initialize(int r,int t,int s)
{
  data.resize(r);
  for(int i=0;i<r;i++){
    data[i].resize(t);
    for(int j=0;j<t;j++){
      data[i][j].resize(s,0.);
    };
  };
};

// LibDataFTable

LibDataFTable::LibDataFTable()
{
  nsig0=-1;
  maxnr=-1;
};

void LibDataFTable::PutNomtft(int i)
{
  nomtft=i;
  no_temp.resize(nomtft,0);
  startgrp_temp.resize(nomtft,0);
  no_rpara.resize(nomtft,0);
  startgrp_ftab.resize(nomtft,0);
  endgrp_ftab.resize(nomtft,0);
  mt_ftab.resize(nomtft,0);
};

void LibDataFTable::PutNsig0(int i)
{
  nsig0=i;
  val_sig0.resize(nsig0,0.);
  val_sig0_log.resize(nsig0,0.);
  ResizeVal_rpara();
};

void LibDataFTable::PutMaxtemp(int i)
{
  maxtemp=i;
  val_temp.resize(i,0.);
  val_temp_log.resize(i,0.);
};

void LibDataFTable::PutMaxnr(int i)
{
  maxnr=i;
  ResizeVal_rpara();
};

void LibDataFTable::ResizeVal_rpara()
{
  if(nsig0!=-1&&maxnr!=-1){
    val_rpara.resize(maxnr);
    for(int i=0;i<maxnr;i++){
      val_rpara[i].resize(nsig0,0.);
    };
  };
};

void LibDataFTable::PutValSig0(int i,real j)
{
  val_sig0[i]=j;
  val_sig0_log[i]=log(j);
};

void LibDataFTable::PutValTemp(int i,real j)
{
  val_temp[i]=j;
  val_temp_log[i]=log(j);
};

void LibDataFTable::Initialize()
{
  fdata.resize(nomtft);
  for(int mt=0;mt<nomtft;mt++){
    int startgrp=startgrp_ftab[mt];
    int endgrp=endgrp_ftab[mt];
    int grplength=endgrp-startgrp+1;
    int norpara=no_rpara[mt];
    int notemp=no_temp[mt];
    fdata[mt].resize(grplength);
    for(int i=0;i<grplength;i++){
      fdata[mt][i].resize(MAX_RID);
      fdata[mt][i][0].Initialize(norpara,notemp,nsig0);
    };
  };
  for(int i=0;i<group;i++){
    num_fdata[i]=1;
  };
};

void LibDataFTable::PutFData(int mt,int g,int rid,int r,int t,int sig0,real val)
{
  CheckInput(mt,g,rid);
  int sttgrp=startgrp_ftab[mt];
  fdata[mt][g-sttgrp][rid].PutData(r,t,sig0,val);
};

void LibDataFTable::PutFData(int mt,int g,int rid,FTable sec)
{
  CheckInput(mt,g,rid);
  int sttgrp=startgrp_ftab[mt];
  fdata[mt][g-sttgrp][rid]=sec;
};

real LibDataFTable::GetFData(int mt,int g,int rid,int r,int t,int sig0)
{
  CheckInput(mt,g,rid);
  int sttgrp=startgrp_ftab[mt];
  return fdata[mt][g-sttgrp][rid].GetData(r,t,sig0);
};

FTable& LibDataFTable::GetFData(int mt,int g,int rid)
{
  CheckInput(mt,g,rid);
  int sttgrp=startgrp_ftab[mt];
  return fdata[mt][g-sttgrp][rid];
};

bool LibDataFTable::CheckInput(int mt,int g,int rid)
{
  if(mt>=nomtft){
    cout<<"Error in CheckInput in LibDataFTable.\n";
    cout<<"You request f-table of mt "<<mt<<"\n";
    cout<<"But f-table exists below "<<nomtft<<"\n";
    exit(0);
  };
  int sttgrp=startgrp_ftab[mt];
  int endgrp=endgrp_ftab[mt];
  if(g<sttgrp||g>endgrp){
    cout<<"Error in CheckInput in LibDataFTable.\n";
    cout<<"You request f-table in group-"<<g<<"\n";
    cout<<"But f-table exists from "<<sttgrp<<" to "<<endgrp<<"\n";
    exit(0);
  };
  if(rid<0||rid>=num_fdata[g]){
    cout<<"Error in CheckInput in LibDataFTable.\n";
    cout<<"You request f-table of R-ID "<<rid<<"\n";
    cout<<"But f-table exists below "<<num_fdata[g]<<"\n";
    exit(0);
  };
  return true;
};

void LibDataFTable::WriteFtableWithSig0Temp(int mt,int g,int r)
{
  string reaction_name[]={
    "(n,f)","(n,g)","(n,n)","(n,t1)","(n,ne)","(n,n')"
  };
  CheckInput(mt,g,0);

  cout.setf(ios::fixed);
  cout.precision(1);
  cout<<"# -------------------------------------------------------------\n";
  cout<<"# Sig0 [b]   f-factor in group "<<g<<" of "<<reaction_name[mt]<<"\n";
  cout<<"#           (Temperature [K])\n";
  cout<<"#            ";
  for(int i=0;i<no_temp[mt];i++){
    cout<<val_temp[i]<<"        ";
  };
  cout<<"\n";
  cout<<"# -------------------------------------------------------------\n";
  cout.unsetf(ios::fixed);

  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<nsig0;i++){
    cout<<val_sig0[i]<<" ";
    for(int j=0;j<no_temp[mt];j++){
      cout<<GetFData(mt,g,0,r,j,i)<<"  ";
    };
    cout<<"\n";
  };
  cout<<"# -------------------------------------------------------------\n";
  cout.unsetf(ios::scientific);
};

void LibDataFTable::WriteFtableWithSig0Rpara(int mt,int g,int t)
{
  string reaction_name[]={
    "(n,f)","(n,g)","(n,n)","(n,t1)","(n,ne)","(n,n')"
  };
  cout.setf(ios::scientific);
  cout.precision(3);
  CheckInput(mt,g,0);
  cout<<"# -------------------------------------------------------------\n";
  cout<<"# Sig0    f-factor in group "<<g<<" of "<<reaction_name[mt]<<"\n";
  cout<<"# [b]    (R-parameter to "<<id_rpara[g][0]<<")\n";
  cout<<"#         ";
  for(int i=0;i<no_rpara[mt];i++){
    cout<<val_rpara[i][0]<<" ";
  };
  cout<<"\n";
  cout<<"# -------------------------------------------------------------\n";
  for(int i=0;i<nsig0;i++){
    cout<<val_sig0[i]<<" ";
    for(int j=0;j<no_rpara[mt];j++){
      cout<<GetFData(mt,g,0,j,t,i)<<" ";
    };
    cout<<"\n";
  };
  cout<<"# -------------------------------------------------------------\n";
};

void LibDataFTable::WriteFtableWithRparaSig0(int mt,int g,int t)
{
  string reaction_name[]={
    "(n,f)","(n,g)","(n,n)","(n,t1)","(n,ne)","(n,n')"
  };
  cout.setf(ios::scientific);
  cout.precision(3);
  CheckInput(mt,g,0);
  cout<<"# -------------------------------------------------------------\n";
  cout<<"# R-para    f-factor in group "<<g<<" of "<<reaction_name[mt]<<"\n";
  cout<<"#("<<id_rpara[g][0]<<")\n";
  cout<<"#         ";
  for(int i=0;i<nsig0;i++){
    cout<<val_sig0[i]<<" ";
  };
  cout<<"\n";
  cout<<"# -------------------------------------------------------------\n";
  for(int i=0;i<no_rpara[mt];i++){
    cout<<val_rpara[i][0]<<" ";
    for(int j=0;j<nsig0;j++){
      cout<<GetFData(mt,g,0,i,t,j)<<" ";
    };
    cout<<"\n";
  };
  cout<<"# -------------------------------------------------------------\n";
};

void LibDataFTable::PutGroup(int i)
{
  group=i;
  id_rpara.resize(group);
  num_fdata.resize(group,0);
  for(int i=0;i<group;i++){
    id_rpara[i].resize(MAX_RID,0);
  };
};

void LibDataFTable::PutIDRpara(int g,int nid,int rid)
{
  if(g<0||g>=group){
    cout<<"Error in LibDataFTable::PutIDRpara.\n";
    cout<<"Energy group "<<g<<" is not appropriated.\n";
    exit(0);
  };
  id_rpara[g][nid]=rid;
};

real LibDataFTable::GetF(int mt,int g,int rid,real rx,real r_sigt,real rtemp,real x)
{
  /*
  int sttgrp=startgrp_ftab[mt];
  int endgrp=endgrp_ftab[mt];
  if(mt>=nomtft||g<sttgrp||g>endgrp){
    return 1.;
  };
  */

  if(no_rpara[mt]==1){
    return GetF(mt,g,rid,0.,rtemp,x+r_sigt);
  }else{
    return GetF(mt,g,rid,rx,rtemp,x);
  };
};

real LibDataFTable::GetF(int mt,int g,int rid,real rx,real rtemp,real x)
{
  int sttgrp=startgrp_ftab[mt];
  int endgrp=endgrp_ftab[mt];
  if(mt>=nomtft||g<sttgrp||g>endgrp){
    return 1.;
  };

  int grp=g-sttgrp;
  int ntemp=no_temp[mt];
  int nxr=no_rpara[mt];

  if(nxr<=1&&x>=val_sig0[nsig0-1])return GetFLargeSig0(mt,g,rtemp);

  real *tabx=new real[nxr];
  real *taby=new real[nxr];

  vector< vector<real> > ftab(ntemp);
  for(int i=0;i<ntemp;i++){
    ftab[i].resize(nsig0,0.);
    for(int j=0;j<nsig0;j++){
      ftab[i][j]=log(fdata[mt][grp][rid].GetData(0,i,j));
      if(nxr!=1&&val_rpara[0][j]!=-1){
	//if(rx>val_rpara[0][j]&&rx>0.01){
	if(rx>val_rpara[0][j]&&rx>val_rpara[1][j]*0.1){	  
	  if(rx>=val_rpara[nxr-1][j]){
	    ftab[i][j]=log(fdata[mt][grp][rid].GetData(nxr-1,i,j));
	  }else if(nxr==2){
	    real r1=val_rpara[0][j];
	    real r2=val_rpara[1][j];
	    real z1=fdata[mt][grp][rid].GetData(0,i,j);
	    real z2=fdata[mt][grp][rid].GetData(1,i,j);
	    real ans=LinearLinearInterpolation(r1,r2,z1,z2,rx);
	    if(ans<=0.)ans=1.;
	    ftab[i][j]=log(ans);
	  }else{
	    real rxlog=log(rx);
	    for(int lr=0;lr<nxr;lr++){
	      real save=val_rpara[lr][j];
	      //if(save<=0.)save=0.01;
	      if(save<=0.)save=val_rpara[1][j]*0.1;	      
	      tabx[lr]=log(save);
	      taby[lr]=log(fdata[mt][grp][rid].GetData(lr,i,j));
	    };
	    real save=CubicSplineFitting(tabx,taby,nxr,rxlog);
	    ftab[i][j]=save;
	  };
	};
      };
    };
  };

  /*
  for(int i=0;i<ntemp;i++){
    for(int j=0;j<nsig0;j++){
      cout<<exp(ftab[i][j])<<" ";
    };
    cout<<"\n";
  };
  */

  delete [] tabx;
  delete [] taby;


  real *xx=new real[ntemp];
  if(ntemp>1){
    for(int i=0;i<ntemp;i++){
      xx[i]=val_temp_log[i];
    };
  };

  real *tabx2=new real[nsig0];
  for(int i=0;i<nsig0;i++){
    tabx2[i]=val_sig0_log[i];
  };

  real temp=rtemp;
  real *yy=new real[ntemp];
  real *taby2=new real[nsig0];

  if(x<=val_sig0[0]||nsig0==1){
    for(int i=0;i<ntemp;i++){
      yy[i]=ftab[i][0];
    };
  }else if(x>=val_sig0[nsig0-1]){
    for(int i=0;i<ntemp;i++){
      yy[i]=ftab[i][nsig0-1];
    };
  }else{
    real xxx=log(x);
    for(int i=0;i<ntemp;i++){
      for(int l=0;l<nsig0;l++){
	taby2[l]=ftab[i][l];
      };
      yy[i]=CubicSplineFitting(tabx2,taby2,nsig0,xxx);
    };
  };

  real ans;
  if(ntemp!=1&&temp>val_temp[0]){
    if(temp>=val_temp[ntemp-1]){
      ans=exp(yy[ntemp-1]);
    }else{
      real tt=log(temp);
      ans=exp(CubicSplineFitting(xx,yy,ntemp,tt));
    };
  }else{
    ans=exp(yy[0]);
  };

  delete [] tabx2;
  delete [] taby2;
  delete [] xx;
  delete [] yy;

  return ans;
};

real LibDataFTable::GetFLargeSig0(int mt,int g,real rtemp)
{
  int ntemp=no_temp[mt];
  int sttgrp=startgrp_ftab[mt];
  int grp=g-sttgrp;

  real *xx=new real[ntemp];
  if(ntemp>1){
    for(int i=0;i<ntemp;i++){
      //xx[i]=log(val_temp[i]);
      xx[i]=val_temp_log[i];
    };
  };

  real temp=rtemp;

  real *yy=new real[ntemp];
  for(int i=0;i<ntemp;i++){
    yy[i]=log(fdata[mt][grp][0].GetData(0,i,nsig0-1));
  };

  real ans=exp(yy[0]);
  if(ntemp!=1&&temp>val_temp[0]){
    if(temp>=val_temp[ntemp-1]){
      ans=exp(yy[ntemp-1]);
    }else{
      real tt=log(temp);
      ans=exp(CubicSplineFitting(xx,yy,ntemp,tt));
    };
  };

  delete [] xx;
  delete [] yy;

  return ans;
};

bool LibDataFTable::CheckSamePoint(LibDataFTable &sec,int g)
{
  if(nomtft!=sec.GetNomtft()){
    cout<<"In LibDataFTable::CheckSamePoint.\n";
    cout<<"NOMTFT is inconsistent.\n";
    cout<<"By which be replaced : "<<sec.GetNomtft()<<"\n";
    cout<<"Is replaced          : "<<nomtft<<"\n";
    return false;
  };
  if(nsig0!=sec.GetNsig0()){
    cout<<"In LibDataFTable::CheckSamePoint.\n";
    cout<<"NSIG0 is inconsistent.\n";
    cout<<"By which be replaced : "<<sec.GetNsig0()<<"\n";
    cout<<"Is replaced          : "<<nsig0<<"\n";
    return false;
  };

  for(int i=0;i<nomtft;i++){
    if(no_rpara[i]!=sec.GetNoRpara(i)){
      cout<<"In LibDataFTable::CheckSamePoint.\n";
      cout<<"NO_RPARA is inconsistent.\n";
      return false;
    };
    if(no_temp[i]!=sec.GetNoTemp(i)){
      cout<<"In LibDataFTable::CheckSamePoint.\n";
      cout<<"NO_TEMP is inconsistent.\n";
      return false;
    };
  };

  return true;
};

void LibDataFTable::CopyData(LibDataFTable &sec,int g)
{
  if(!CheckSamePoint(sec,g)){
    cout<<"LibDataFTable::CopyData is not done.\n";
    return;
  };

  for(int i=0;i<nomtft;i++){
    for(int j=0;j<no_rpara[i];j++){
      for(int k=0;k<no_temp[i];k++){
	for(int l=0;l<nsig0;l++){
	  real tmp=sec.GetFData(i,g,0,j,k,l);
	  PutFData(i,g,0,j,k,l,tmp);
	};
      };
    };
  };

  for(int i=0;i<num_fdata[g];i++){
    id_rpara[g][i]=sec.GetIDRpara(g,i);
  };
};

void LibDataFTable::AddData(LibDataFTable &sec,int g)
{
  if(!CheckSamePoint(sec,g)){
    cout<<"# LibDataFTable::AddData is not done.\n";
    return;
  };
  if(num_fdata[g]==MAX_RID){
    cout<<"# Please increase MAX_RID in LibData.h.\n";
    exit(0);
  };

  int add_rid=sec.GetIDRpara(g,0);
  for(int i=0;i<num_fdata[g];i++){
    int rid=GetIDRpara(g,i);
    if(rid==add_rid){
      cout<<"# f-table already exists in group "<<g<<".\n";
      cout<<"# ID : "<<rid<<"\n";
    };
  };
  id_rpara[g][num_fdata[g]]=add_rid;
  num_fdata[g]++;

  int ttt=num_fdata[g]-1;
  for(int i=0;i<nomtft;i++){
    int sttgrp=startgrp_ftab[i];
    if(sttgrp!=-1){
      fdata[i][g-sttgrp][ttt].Initialize(no_rpara[i],no_temp[i],nsig0);
      for(int j=0;j<no_rpara[i];j++){
        for(int k=0;k<no_temp[i];k++){
  	  for(int l=0;l<nsig0;l++){
	    real tmp=sec.GetFData(i,g,0,j,k,l);
	    PutFData(i,g,ttt,j,k,l,tmp);
	  };
	};
      };
    };
  };
};

void LibDataFTable::ShowSelf()
{
  for(int i=0;i<nomtft;i++){
    cout<<i<<" "<<mt_ftab[i]<<"\n";
  };
};

// LibDataChiVector

LibDataChiVector::LibDataChiVector()
{
  no_vector=-1;
  group=-1;
};

void LibDataChiVector::Initialize()
{
  if(no_vector!=-1&&group!=-1){
    chiv.resize(no_vector);
    for(int i=0;i<no_vector;i++){
      chiv[i].put_imax(group);
    };
  };
};

void LibDataChiVector::PutGrp(int g)
{
  group=g;
  vectorid.resize(group);
  Initialize();
};

void LibDataChiVector::PutNoVector(int i)
{
  no_vector=i;
  Initialize();
};

void LibDataChiVector::PutVectorID(int i,int j)
{
  if(i<0||i>=group){
    cout<<"# Error in PutVectorID in LibDataChiVector.\n";
    cout<<"# You request incorrect group : "<<i<<"\n";
    exit(0);
  };
  if(j<0||j>=no_vector){
    cout<<"# Error in PutVectorID in LibDataChiVector.\n";
    cout<<"# You request incorrect vector ID : "<<j<<"\n";
    exit(0);
  };
  vectorid[i]=j;
};

void LibDataChiVector::ShowVectorID()
{
  for(int i=0;i<group;i++){
    cout<<i<<" "<<vectorid[i]<<"\n";
  };
};

void LibDataChiVector::AssignChiVectorToEachGroup()
{
  vector<GroupData1D> chi_org(no_vector);
  for(int i=0;i<no_vector;i++){
    chi_org[i].copy(chiv[i]);
  };
  vector<int> vecid_org(group);
  for(int i=0;i<group;i++){
    vecid_org[i]=vectorid[i];
  };

  PutNoVector(group);
  for(int i=0;i<group;i++){
    int id=vecid_org[i];
    PutVectorID(i,i);
    chiv[i].copy(chi_org[id]);
  };
};

// ThermalScatteringData

ThermalScatteringData::ThermalScatteringData()
{
  exist=false;
  grp=-1;
  init_grp=-1;
  pl=-1;
  ntemp=-1;
};

void ThermalScatteringData::PutTemperature(int i,vector<real> t)
{
  ntemp=i;
  temp.resize(ntemp,0.);
  for(int j=0;j<ntemp;j++){
    temp[j]=t[j];
  };
  for(int j=1;j<ntemp;j++){
    if(temp[j]<=temp[j-1]){
      cout<<"Error in ThermalScatteringData::PutTemperature.\n";
      cout<<"Not-appropriated value "<<temp[j]<<"&"<<temp[j-1]<<"\n";
      cout<<"(Please set aseinding order.)\n";
      exit(0);
    };
  };
};

void ThermalScatteringData::Initialize(int ig,int igrp,int ipl)
{
  exist=true;
  grp=ig;

  if(ntemp==-1){
    cout<<"Error in ThermalScatteringData::Initialize.\n";
    cout<<"Please do [PutTemperature] before [Initialize].\n";
    exit(0);
  };
  if(igrp<0||igrp>=grp){
    cout<<"Error in ThermalScatteringData::Initialize.\n";
    cout<<"Not-appropriated initial_grp ("<<igrp<<")\n";
    exit(0);
  };
  init_grp=igrp;
  pl=ipl;

  int grpband=grp-init_grp;

  data.resize(ntemp);
  for(int i=0;i<ntemp;i++){
    data[i].resize(pl+1);
    for(int j=0;j<=pl;j++){
      data[i][j].put_yx(grpband,grpband);
      data[i][j].set_zero();
    };
  };
};

void ThermalScatteringData::ReadFile(string mdir,string ss)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  int tmp;
  real vtmp;

  int ig; // Total number of energy groups
  fin>>ig; 

  fin>>tmp;
  vector<real> tin(tmp);
  for(int i=0;i<tmp;i++){
    fin>>tin[i];
  };
  PutTemperature(tmp,tin);

  int ip;  // Legendre order
  fin>>ip;
  int ig2; // The highest energy group of thermal scattering
  fin>>ig2;
  Initialize(ig,ig2-1,ip); 

  int igband=ig-ig2+1;

  for(int tt=0;tt<ntemp;tt++){
  for(int i=0;i<=ip;i++){
    for(int j=0;j<igband;j++){
      for(int k=0;k<igband;k++){
	fin>>vtmp;
	if(vtmp!=0.)data[tt][i].put_data(j,k,vtmp);
      };
    };
  };
  };

  init_src_grp=ig2-1;
  bool nonzero_src=true;
  for(int i=0;i<igband;i++){
    bool nonzero=false;
    for(int j=0;j<igband;j++){
      if(data[0][0].get_dat(i,j)!=0.)nonzero=true;
    };
    if(!nonzero&&nonzero_src){
      init_src_grp++;
    };
    if(nonzero&&nonzero_src){
      nonzero_src=false;
    };
  };
};

void ThermalScatteringData::ReadFileSlide(string mdir,string ss,int dg)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  int tmp;
  real vtmp;

  int ig; // Total number of energy groups
  fin>>ig;

  ig+=dg; // SLIDE

  fin>>tmp;
  vector<real> tin(tmp);
  for(int i=0;i<tmp;i++){
    fin>>tin[i];
  };
  PutTemperature(tmp,tin);

  int ip;  // Legendre order
  fin>>ip;
  int ig2; // The highest energy group of thermal scattering
  fin>>ig2;

  ig2+=dg; // SLIDE

  Initialize(ig,ig2-1,ip); 

  int igband=ig-ig2+1;

  for(int tt=0;tt<ntemp;tt++){
  for(int i=0;i<=ip;i++){
    for(int j=0;j<igband;j++){
      for(int k=0;k<igband;k++){
	fin>>vtmp;
	if(vtmp!=0.)data[tt][i].put_data(j,k,vtmp);
      };
    };
  };
  };

  init_src_grp=ig2-1;
  bool nonzero_src=true;
  for(int i=0;i<igband;i++){
    bool nonzero=false;
    for(int j=0;j<igband;j++){
      if(data[0][0].get_dat(i,j)!=0.)nonzero=true;
    };
    if(!nonzero&&nonzero_src){
      init_src_grp++;
    };
    if(nonzero&&nonzero_src){
      nonzero_src=false;
    };
  };
};

real ThermalScatteringData::GetData(int itemp,int ipl,int srcg,int sinkg)
{
  if(srcg<init_grp||sinkg<init_grp)return 0.;
  return data[itemp][ipl].get_dat(srcg-init_grp,sinkg-init_grp);
};

real ThermalScatteringData::GetData(int ipl,int srcg,int sinkg,real tt)
{
  if(srcg<init_grp||sinkg<init_grp)return 0.;
  int pos_srcg=srcg-init_grp;
  int pos_sinkg=sinkg-init_grp;
  if(tt<temp[0])return data[0][ipl].get_dat(pos_srcg,pos_sinkg);
  if(tt>temp[ntemp-1])return data[ntemp-1][ipl].get_dat(pos_srcg,pos_sinkg);
  for(int i=0;i<ntemp;i++){
    real diff=abs(temp[i]/tt-1.);
    if(diff<1e-5)return data[i][ipl].get_dat(pos_srcg,pos_sinkg);
  };

  // (interpolation with temperature)
  real mindel=10000.;
  int idmax=-1;
  for(int t=0;t<ntemp;t++){
    real tmp=fabs(tt-temp[t]);
    if(mindel>tmp){
      idmax=t;
      mindel=tmp;
    };
  };

  if(tt<temp[1]){
    real sig1=data[0][ipl].get_dat(pos_srcg,pos_sinkg);
    real sig2=data[1][ipl].get_dat(pos_srcg,pos_sinkg);
    real t1=temp[0];
    real t2=temp[1];
    return sig1+(sig2-sig1)*(tt-t1)/(t2-t1);
  };
  if(tt>temp[ntemp-2]){
    real sig1=data[ntemp-2][ipl].get_dat(pos_srcg,pos_sinkg);
    real sig2=data[ntemp-1][ipl].get_dat(pos_srcg,pos_sinkg);
    real t1=temp[ntemp-2];
    real t2=temp[ntemp-1];
    return sig1+(sig2-sig1)*(tt-t1)/(t2-t1);
  };

  real sig1=data[idmax-1][ipl].get_dat(pos_srcg,pos_sinkg);
  real sig2=data[idmax][ipl].get_dat(pos_srcg,pos_sinkg);
  real sig3=data[idmax+1][ipl].get_dat(pos_srcg,pos_sinkg);
  real t1=temp[idmax-1];
  real t2=temp[idmax];
  real t3=temp[idmax+1];

  real tmp1=(tt-t2)*(tt-t3)*sig1/((t1-t2)*(t1-t3))
    +(tt-t1)*(tt-t3)*sig2/((t2-t1)*(t2-t3))
    +(tt-t1)*(tt-t2)*sig3/((t3-t1)*(t3-t2));

  if(ipl==0&&tmp1<-1e-4){
    cout<<"!! Warning !!\n";
    cout<<"Negative value ("<<tmp1<<") is calculated in temperature intrapolation in P0 thermal matrix.\n";
    //cout<<"   "<<t1<<" : "<<sig1<<"\n";
    //cout<<"   "<<t2<<" : "<<sig2<<"\n";
    //cout<<"   "<<t3<<" : "<<sig3<<"\n";
    //cout<<" (You requested)\n";
    //cout<<"   "<<tt<<" : "<<tmp1<<"\n";
  };

  if(ipl==0&&tmp1<0.)tmp1=0.; // Negative value is set to zero.
  return tmp1;
};

// LibData

LibData::LibData()
{
  xsdata.Init("MicroCrossSection");
  exist_chi_vector=false;
  exist_bell_factor=false;
  mat=-1;
};

void LibData::ReadFile(string mdir,string ss,MATIDTranslator &midt)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  int tmp;
  real vtmp;

  // control data
  //fin>>mat;
  string matnn;
  fin>>matnn;
  string matnn1=matnn.substr(0,1);
  if(matnn1=="1"||
     matnn1=="2"||
     matnn1=="3"||
     matnn1=="4"||
     matnn1=="5"||
     matnn1=="6"||
     matnn1=="7"||
     matnn1=="8"||
     matnn1=="9"){
    mat=StringToInt(matnn);
    mat=midt.GetMATIDFromENDFID(mat);    
  }else{
    mat=99999; // dummy
  };
  fin>>tmp;
  PutGroup(tmp);
  ftable.PutGroup(tmp);
  int ifiss;
  fin>>ifiss;
  exist_fiss=false;
  if(ifiss==1)exist_fiss=true;
  int ichvec;
  fin>>ichvec;
  fin>>tmp;
  int nomt2d=tmp;
  //maxpl1.resize(nomt2d,0);
  SetSizeMaxpl1(nomt2d);
  fin>>tmp;
  ftable.PutNomtft(tmp);
  int nomtft=tmp;
  fin>>tmp;
  ftable.PutNsig0(tmp);
  int nsig0=tmp;
  fin>>tmp;
  ftable.PutMaxtemp(tmp);
  int maxtemp=tmp;
  fin>>tmp;
  ftable.PutMaxnr(tmp);
  int maxnr=tmp;

  for(int i=0;i<nomt2d;i++){
    fin>>tmp;
    //maxpl1.at(i)=tmp+1;
    PutMaxpl1(i,tmp+1);
    xsdata.PutDim2d(i,tmp+1);
  };
  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    ftable.PutNoTemp(i,tmp);
  };
  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    ftable.PutStartGrpTemp(i,tmp);
  };
  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    ftable.PutNoRpara(i,tmp);
  };

  int stt=group-1;
  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    if(tmp<stt&&tmp!=0)stt=tmp;
    tmp--;
    ftable.PutStartGrpFtab(i,tmp);
  };
  stt--;
  ftable.PutStartGrpAllMF(stt);

  int sed=0;
  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    if(tmp>sed&&tmp!=0)sed=tmp;
    tmp--;
    ftable.PutEndGrpFtab(i,tmp);
  };
  sed--;
  ftable.PutEndGrpAllMF(sed);

  for(int i=0;i<nomtft;i++){
    fin>>tmp;
    ftable.PutMtFtab(i,tmp);
  };
  for(int i=0;i<nsig0;i++){
    fin>>vtmp;
    ftable.PutValSig0(i,vtmp);
  };
  for(int i=0;i<maxtemp;i++){
    fin>>vtmp;
    ftable.PutValTemp(i,vtmp);
  };

  for(int i=0;i<maxnr;i++){
    for(int j=0;j<nsig0;j++){
      fin>>vtmp;
      ftable.PutValRpara(i,j,vtmp);
    };
  };

  /*
  cout<<"# f-table information\n";
  cout<<"#\n";
  cout<<"# The number of reactions : "<<ftable.GetNomtft()<<"\n";
  for(int i=0;i<ftable.GetNomtft();i++){
    cout<<"# Mt : "<<ftable.GetMtFtab(i)<<"\n";
  };
  */

  // 1dxs 
  if(ifiss==1){
    for(int i=0;i<group;i++){
      fin>>vtmp;
      xsdata.GetData1d(sigf).put_data(i,vtmp);
    };
    for(int i=0;i<group;i++){
      fin>>vtmp;
      xsdata.GetData1d(nu).put_data(i,vtmp);
    };
  }else{
    xsdata.GetData1d(sigf).set_zero();
    xsdata.GetData1d(nu).set_zero();
  };

  for(int i=0;i<group;i++){
    fin>>vtmp;
    xsdata.GetData1d(sigc).put_data(i,vtmp);
  };
  for(int i=0;i<group;i++){
    fin>>vtmp;
    xsdata.GetData1d(siginel).put_data(i,vtmp);
  };
  for(int i=0;i<group;i++){
    fin>>vtmp;
    xsdata.GetData1d(sigel).put_data(i,vtmp);
  };
  for(int i=0;i<group;i++){
    fin>>vtmp;
    xsdata.GetData1d(mu).put_data(i,vtmp);
  };
  for(int i=0;i<group;i++){
    fin>>vtmp;
    xsdata.GetData1d(sign2n).put_data(i,vtmp);
  };
  if(ichvec!=0){
    for(int i=0;i<group;i++){
      fin>>vtmp;
      xsdata.GetData1d(chi).put_data(i,vtmp);
    };
  };
  for(int i=0;i<group;i++){
    real sum=0.;
    sum+=xsdata.GetData1d(sigc).get_dat(i);
    sum+=xsdata.GetData1d(sigf).get_dat(i);
    sum+=xsdata.GetData1d(sigel).get_dat(i);
    sum+=xsdata.GetData1d(siginel).get_dat(i);
    sum+=xsdata.GetData1d(sign2n).get_dat(i);
    xsdata.GetData1d(sigt).put_data(i,sum);
  };

  // 2dxs
  for(int mt=0;mt<nomt2d;mt++){
    int endgrp;
    fin>>endgrp;
    for(int p=0;p<maxpl1[mt];p++){
      xsdata.GetData2d(mt,p).set_zero();
      for(int k=0;k<endgrp;k++){
	int endgrp2;
	fin>>endgrp2;
	for(int j=0;j<endgrp2;j++){
	  fin>>vtmp;
	  xsdata.GetData2d(mt,p).put_data(k,k+j,vtmp);
	};
      };
      //if(mt==0)xsdata.GetData2d(mt,p).ReducedForm();
    };
  };

  ftable.Initialize();
  // f-table
  for(int mt=0;mt<nomtft;mt++){

#if 0
    cout<<"# MT-index : "<<mt<<"\n";
    cout<<"#   start-grp     : "<<ftable.GetStartGrpFtab(mt)<<"\n";
    cout<<"#   end-grp       : "<<ftable.GetEndGrpFtab(mt)<<"\n";
    cout<<"#   No. of R-para : "<<ftable.GetNoRpara(mt)<<"\n";
    cout<<"#   No. of Temp   : "<<ftable.GetNoTemp(mt)<<"\n";
    cout<<"#   No. of sig0   : "<<ftable.GetNsig0()<<"\n";
#endif
    
    for(int i=ftable.GetStartGrpFtab(mt);i<=ftable.GetEndGrpFtab(mt);i++){
      for(int n=0;n<ftable.GetNoRpara(mt);n++){
	for(int k=0;k<ftable.GetNoTemp(mt);k++){
	  for(int j=0;j<ftable.GetNsig0();j++){
	    fin>>vtmp;
	    //cout<<i<<" "<<n<<" "<<k<<" "<<j<<" "<<vtmp<<"\n";
	    if(i!=-1)ftable.PutFData(mt,i,0,n,k,j,vtmp);
	  };
	};
      };
    };
  };

  // chi-vector
  if(ichvec>1){
    exist_chi_vector=true;
    fin>>ichvec;
    chiv.PutNoVector(ichvec);
    fin>>tmp;
    int mxd=tmp;
    //cout<<"# ichvec / mxd : "<<ichvec<<" "<<mxd<<"\n";
    for(int i=0;i<group;i++){
      fin>>tmp;
      chiv.PutVectorID(i,tmp-1);
    };
    for(int k=0;k<ichvec;k++){
      for(int j=0;j<mxd;j++){
	fin>>vtmp;
	chiv.GetChiVector(k).put_data(j,vtmp);
      };
    };
  };

  // Nuclide ID for R-parameter
  if(maxnr>1){
    int tmp;
    fin>>tmp; // if -1, group-wise parameter is defined
    if(tmp>0){
      tmp=TranslateNuclideIDFromJFS(tmp);
      tmp=midt.GetMATIDFromENDFID(tmp);
      for(int i=0;i<group;i++){
        ftable.PutIDRpara(i,0,tmp);
      };
    }else if(tmp==-1){
      for(int i=0;i<group;i++){
        int tmp;
        fin>>tmp;
	if(tmp>0&&tmp<1000)tmp=TranslateNuclideIDFromJFS(tmp); // To ENDFID
        tmp=midt.GetMATIDFromENDFID(tmp);
        ftable.PutIDRpara(i,0,tmp);
      };
    };
  };

  fin.close();
};

void LibData::WriteFile(string mdir,string ss)
{
  mdir.append(ss);
  ofstream fout;
  fout.open(mdir.data(),ios::out);
  fout.setf(ios::scientific);
  fout.precision(8);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<ss<<"\n";
    exit(0);
  };

  // control data
  fout<<mat<<"\n";
  fout<<group<<"\n";

  int ifiss=0;
  if(exist_fiss)ifiss=1;
  fout<<ifiss<<"\n";

  int ichvec=0;
  if(exist_fiss)ichvec=1;
  if(exist_chi_vector){
    ichvec=chiv.GetNoVector();
  };
  fout<<ichvec<<"\n";

  int nomt2d=maxpl1.size();
  fout<<nomt2d<<"\n";

  int nomtft=ftable.GetNomtft();
  fout<<nomtft<<"\n";
  int nsig0=ftable.GetNsig0();
  fout<<nsig0<<"\n";
  int maxtemp=ftable.GetMaxtemp();
  fout<<maxtemp<<"\n";
  int maxnr=ftable.GetMaxnr();
  fout<<maxnr<<"\n";
  for(int i=0;i<nomt2d;i++){
    fout<<xsdata.GetDim2d(i)-1<<"\n";
  };
  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetNoTemp(i)<<"\n";
  };
  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetStartGrpTemp(i)<<"\n";
  };
  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetNoRpara(i)<<"\n";
  };

  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetStartGrpFtab(i)+1<<"\n";
  };

  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetEndGrpFtab(i)+1<<"\n";
  };

  for(int i=0;i<nomtft;i++){
    fout<<ftable.GetMtFtab(i)<<"\n";
  };
  for(int i=0;i<nsig0;i++){
    fout<<ftable.GetSig0(i)<<"\n";
  };
  for(int i=0;i<maxtemp;i++){
    fout<<ftable.GetTemp(i)<<"\n";
  };

  for(int i=0;i<maxnr;i++){
    for(int j=0;j<nsig0;j++){
      fout<<ftable.GetRpara(i,j)<<"\n";
    };
  };

  // 1dxs 
  if(ifiss==1){
    for(int i=0;i<group;i++){
      fout<<xsdata.GetData1d(sigf).get_dat(i)<<"\n";
    };
    for(int i=0;i<group;i++){
      fout<<xsdata.GetData1d(nu).get_dat(i)<<"\n";
    };
  };

  for(int i=0;i<group;i++){
    fout<<xsdata.GetData1d(sigc).get_dat(i)<<"\n";
  };
  for(int i=0;i<group;i++){
    fout<<xsdata.GetData1d(siginel).get_dat(i)<<"\n";
  };
  for(int i=0;i<group;i++){
    fout<<xsdata.GetData1d(sigel).get_dat(i)<<"\n";
  };
  for(int i=0;i<group;i++){
    fout<<xsdata.GetData1d(mu).get_dat(i)<<"\n";
  };
  for(int i=0;i<group;i++){
    fout<<xsdata.GetData1d(sign2n).get_dat(i)<<"\n";
  };

  if(ichvec!=0){
    for(int i=0;i<group;i++){
      fout<<xsdata.GetData1d(chi).get_dat(i)<<"\n";
    };
  };

  // 2dxs
  for(int mt=0;mt<nomt2d;mt++){
    int src_grp_end=0;
    for(int p=0;p<maxpl1[mt];p++){
      for(int k=0;k<group;k++){
	for(int j=k;j<group;j++){
	  real tmp=xsdata.GetData2d(mt,p).get_dat(k,j);
	  if(tmp!=0.&&src_grp_end<k)src_grp_end=k;
	};
      };
    };
    src_grp_end++;
    fout<<src_grp_end<<"\n";
    for(int p=0;p<maxpl1[mt];p++){
      for(int k=0;k<src_grp_end;k++){
        int sink_grp_end=k-1;      
	for(int j=k;j<group;j++){
          real tmp=xsdata.GetData2d(mt,p).get_dat(k,j);
	  if(tmp!=0.&&sink_grp_end<j)sink_grp_end=j;
	};
        sink_grp_end++;
	fout<<sink_grp_end-k<<"\n";
	for(int j=k;j<sink_grp_end;j++){
	  fout<<xsdata.GetData2d(mt,p).get_dat(k,j)<<"\n";
	};
      };
    };
  };

  // f-table
  for(int mt=0;mt<nomtft;mt++){
    for(int i=ftable.GetStartGrpFtab(mt);i<=ftable.GetEndGrpFtab(mt);i++){
      // ..... The following part is modified in 2020/3/10 .....
      /*
      if(i!=-1){
      for(int n=0;n<ftable.GetNoRpara(mt);n++){
	for(int k=0;k<ftable.GetNoTemp(mt);k++){
	  for(int j=0;j<ftable.GetNsig0();j++){
            fout<<ftable.GetFData(mt,i,0,n,k,j)<<"\n";
	  };
	};
      };
      };
      */
      for(int n=0;n<ftable.GetNoRpara(mt);n++){
	for(int k=0;k<ftable.GetNoTemp(mt);k++){
	  for(int j=0;j<ftable.GetNsig0();j++){
            if(i!=-1){
	      fout<<ftable.GetFData(mt,i,0,n,k,j)<<"\n";
	    }else{
	      fout<<"   0\n"; // to keep consistensy to [ReadFile].
	    };
	  };
	};
      };
      // ..........................................................
    };
  };

  // chi-vector
  if(ichvec>1){
    fout<<ichvec<<"\n";
    fout<<group<<"\n";
    //cout<<"# ichvec/group : "<<ichvec<<" "<<group<<"\n";
    for(int i=0;i<group;i++){
      fout<<chiv.GetVectorID(i)+1<<"\n";
      //cout<<chiv.GetVectorID(i)+1<<"\n";
    };
    for(int k=0;k<ichvec;k++){
      for(int j=0;j<group;j++){
	fout<<chiv.GetChiVector(k).get_dat(j)<<"\n";
      };
    };
  };

  // Nuclide ID for R-parameter
  if(maxnr>1){
    int tmp=-1;
    fout<<tmp<<"\n";
    for(int i=0;i<group;i++){
      int tmp=ftable.GetIDRpara(i,0);
      fout<<tmp<<"\n";
    };
  };

  fout.close();
};

void LibData::WriteFileUNCFormat(string mdir,string ss)
{
  mdir.append(ss);
  ofstream fout;
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"(file name) "<<ss<<"\n";
    exit(0);
  };

  fout<<group<<"\n";

  // 1D cross section
  int mt1d[]={18,102,452,251,2,181,4,16,1};
  for(int i=0;i<9;i++){
    fout<<" 1\n";
    fout<<" "<<mat<<"\n";
    fout<<" "<<mt1d[i]<<"\n";
    for(int j=0;j<group;j++){
      real f=1.;
      //if(i==1){f=GetF(1,j,0,0.,300.,0.);};
      //if(i==4){f=GetF(2,j,0,0.,300.,0.);};
      fout<<" "<<xsdata.GetData1d(i,0).get_dat(j)*f<<"\n";
    };
  };

  // 2D cross section
  int mt2d[]={2,4,16};
  for(int l=0;l<2;l++){
    for(int i=0;i<3;i++){
      fout<<" 2\n";
      fout<<" "<<mat<<"\n";
      int tmp=mt2d[i];
      if(l==1)tmp+=1000;
      fout<<" "<<tmp<<"\n";
      for(int j=0;j<group;j++){
        real f=1.;
        //if(i==0){f=GetF(2,j,0,0.,300.,0.);};
	for(int k=j;k<group;k++){
	  fout<<" "<<xsdata.GetData2d(i,l).get_dat(j,k)*f<<"\n";
	};
      };
    };
  };

  fout.close();
};

void LibData::PutGroup(int g)
{
  group=g;
  xsdata.PutGrp(g);
  chiv.PutGrp(g);
  ptable.PutGroup(g);
  bell_factor.resize(g,0.);
};

LibDataChiVector &LibData::GetChiVector()
{
  if(!exist_chi_vector){
    cout<<"# Chi vector does not exist.\n";
    exit(0);
  };
  return chiv;
};

bool LibData::ExistThermalData()
{
  if(!thscat.DoesExist())return false;

  int grp_thscat=thscat.GetGrp();
  if(grp_thscat!=group){
    cout<<"\n";
    cout<<"# Error in XSLibrary::ExistThermalData.\n";
    cout<<"# Thermal data exists in this library.\n";
    cout<<"# but the number of groups is inconsistent with the library.\n";
    cout<<"\n";    
    exit(0);
  };

  return true;
};


void LibData::PutBellFactor(real *inp)
{
  for(int i=0;i<group;i++){
    real tmp=inp[i];
    /*
    if(tmp<1.||tmp>2.){
      cout<<" !!! Warning !!! \n";
      cout<<"Error in LibData::PutBellFactor.\n";
      cout<<"Bell factor should be 1.0 - 2.0.\n";
      cout<<"You set "<<tmp<<"\n";
    };
    */
    bell_factor[i]=tmp;
  };
  exist_bell_factor=true;
};

void LibData::PutBellFactorAllGroup(real inp)
{
  for(int i=0;i<group;i++){
    bell_factor[i]=inp;
  };
  exist_bell_factor=true;
};

void LibData::PutBellFactor(int g, real inp)
{
  if(!exist_bell_factor){
    cout<<"# Error in LibData::PutBellFactor.\n";
    cout<<"# Group-wise Bell factor cannot be set\n";
    cout<<"# because initial Bell factors have not yet been given.\n";
    exit(0);
  };

  if(g<0||g>=group){
    cout<<"# Error in LibData::PutBellFactor.\n";
    cout<<"# Energy group "<<g<<" is inappropriate.\n";
    exit(0);
  };

  bell_factor[g]=inp;
};


void LibData::CalF(int g,int rid,real r,real r_sigt,real t,real sig0,vector<real> &f)
{
  f[0]=GetF(0,g,rid,r,r_sigt,t,sig0); // fis
  f[1]=GetF(1,g,rid,r,r_sigt,t,sig0); // cap
  f[2]=GetF(2,g,rid,r,r_sigt,t,sig0); // ela
  f[3]=GetF(3,g,rid,r,r_sigt,t,sig0); // tot
  f[4]=GetF(4,g,rid,r,r_sigt,t,sig0); // rem
  f[5]=GetF(5,g,rid,r,r_sigt,t,sig0); // ine
};

void LibData::Update1DChiDataForFastReactorCalculation()
{
  if(exist_chi_vector){
    // +++ chi matrix construction
    GroupData2D chimat(group,group);
    for(int g=0;g<group;g++){
      int vecid=GetChiVector().GetVectorID(g);
      for(int g2=0;g2<group;g2++){
	real val=GetChiVector(vecid).get_dat(g2);
	chimat.put_data(g,g2,val);
      };
    };
    // +++ 1D chi data is updated by assuming chi as neutron flux
    GroupData1D chi1d(group);
    for(int g=0;g<group;g++){
      real tmp1=0.;
      real tmp2=0.;
      for(int g2=0;g2<group;g2++){
        real tmp3=GetXSData().GetData1d(sigf).get_dat(g2)*GetXSData().GetData1d(chi).get_dat(g2);
	tmp1+=tmp3*chimat.get_dat(g2,g);
	tmp2+=tmp3;
      };
      chi1d.put_data(g,tmp1/tmp2);
    };
    real sum=chi1d.get_sum();
    sum=1./sum;
    chi1d.Factorize(sum);
    GetXSData().GetData1d(chi).copy(chi1d);
  };
};

void LibData::CrossSectionPerturbation(int g, enum xstype sigx, real r_delta)
{
  real org=GetXSData().GetData1d(sigx).get_dat(g);
  GetXSData().GetData1d(sigx).add_data(g,org*r_delta);

  if(sigx==chi){
    // Fission spectrum normalization for 1D vector
    real sum=GetXSData().GetData1d(sigx).get_sum();
    sum=1./sum;
    GetXSData().GetData1d(sigx).copy(GetXSData().GetData1d(sigx)*sum);
    // Fission spectrum matrix perturbation
    if(exist_chi_vector){
      int nv=GetChiVector().GetNoVector();
      for(int i=0;i<nv;i++){
	real org=GetChiVector().GetChiVector(i).get_dat(g);
	GetChiVector().GetChiVector(i).add_data(g,org*r_delta);
        real sum=GetChiVector().GetChiVector(i).get_sum();
	sum=1./sum;
	GetChiVector().GetChiVector(i).copy(GetChiVector().GetChiVector(i)*sum);
      };
    };
  };

  if(sigx==sigf||sigx==sigc||sigx==sigel||sigx==siginel||sigx==sign2n){
    GetXSData().GetData1d(sigt).add_data(g,org*r_delta);
  };

  if(sigx==sigel||sigx==siginel||sigx==sign2n){
    GetXSData().XSMatrixNormalizationTo1DXS();
  };

};


void LibData::ConstantCrossSectionData(int mat_inp, int g, real sigc_inp)
{
  mat=mat_inp;
  PutGroup(g);  
  exist_fiss=false;

  xsdata.GetData1d(sigf).set_zero();
  xsdata.GetData1d(nu).set_zero();
  xsdata.GetData1d(mu).set_zero();
  xsdata.GetData1d(sigel).set_zero();
  xsdata.GetData1d(chi).set_zero();
  xsdata.GetData1d(siginel).set_zero();
  xsdata.GetData1d(sign2n).set_zero();
  for(int i=0;i<group;i++){
    xsdata.GetData1d(sigc).put_data(i, sigc_inp);
    xsdata.GetData1d(sigt).put_data(i, sigc_inp);
  };

  xsdata.GetData2d(sigel).set_zero();
  xsdata.GetData2d(siginel).set_zero();
  xsdata.GetData2d(sign2n).set_zero();
};

void LibData::AddFtableAtNewTemperaturePoint(LibData &sec)
{
  // Detail can be found in my notebook in 2021/9/15.
  
  LibDataFTable secf=sec.GetFtable();
  real temp_sec=secf.GetTemp(0);

  
  // ... consistency check
  bool consistency=false;
  while(consistency==false){

    if(mat!=sec.GetMat())break;
    if(group!=sec.GetGroup())break;

    if(ftable.GetNsig0()!=secf.GetNsig0())break;
    for(int i=0;i<ftable.GetNsig0();i++){
      real tmp1=ftable.GetSig0(i);
      real tmp2=secf.GetSig0(i);
      if(fabs(tmp1/tmp2-1.)>1e-5)break;
    }

    if(ftable.GetNomtft()!=secf.GetNomtft())break;
    for(int i=0;i<ftable.GetNomtft();i++){
      if(ftable.GetMtFtab(i)!=secf.GetMtFtab(i))break;
    };

    for(int i=0;i<ftable.GetNomtft();i++){
      if(ftable.GetNoRpara(i)>1||secf.GetNoRpara(i)>1){
	cout<<"# R-parameter-implemented library is detected.\n";
	break;
      };
      if(ftable.GetStartGrpFtab(i)!=secf.GetStartGrpFtab(i)){
	cout<<"# Start-group of f-table is inconsistent.\n";
	break;
      };
      if(ftable.GetEndGrpFtab(i)!=secf.GetEndGrpFtab(i)){
	cout<<"# End-group of f-table is inconsistent.\n";
	break;
      };
    };

    for(int i=0;i<ftable.GetNomtft();i++){    
      if(secf.GetNoTemp(i)>1){
	cout<<"# The number of temperature points is larger than 1.\n";
	break;
      };
      int num_temp=ftable.GetNoTemp(i);
      for(int j=0;j<num_temp;j++){
	real temp=ftable.GetTemp(j);
        if(fabs(temp/temp_sec-1.)<1e-4){
          cout<<"# Error in LibData::AddFtableAtNewTemperaturePoint.\n";
	  cout<<"# The same temperature point is defined in the base library.\n";
	  exit(0);
	};
      };
    };
    
    consistency=true;
  };

  if(!consistency){
    cout<<"# Error in LibData::AddFtableAtNewTemperaturePoint.\n";
    cout<<"# Two library data are inconsistent with each other.\n";
    exit(0);
  };
    
  enum xstype sigx[]={sigf,sigc,sigel,sigt,sigt,siginel};
  for(int i=0;i<ftable.GetNomtft();i++){

    if(ftable.GetStartGrpFtab(i)!=-1){
      
    for(int g=ftable.GetStartGrpFtab(i);g<=ftable.GetEndGrpFtab(i);g++){

      real siginf=0.;
      if(i!=4){
	siginf=GetXSData().GetData1d(sigx[i]).get_dat(g);
      }else{
	for(int g2=g+1;g2<group;g2++){
	  siginf+=GetXSData().GetData2d(sigel).get_dat(g,g2);
	};
      };

      real siginf2=0.;
      if(i!=4){
	siginf2=sec.GetXSData().GetData1d(sigx[i]).get_dat(g);
      }else{
	for(int g2=g+1;g2<group;g2++){
	  siginf2+=sec.GetXSData().GetData2d(sigel).get_dat(g,g2);
	};
      };

      int num_temp=ftable.GetNoTemp(i);
      FTable ftab_new;
      
      ftab_new.Initialize(1,num_temp+1,ftable.GetNsig0());
      for(int j=0;j<num_temp;j++){
	for(int k=0;k<ftable.GetNsig0();k++){
	  ftab_new.PutData(0,j,k,ftable.GetFData(i,g,0).GetData(0,j,k));
	};
      };
      
      for(int j=0;j<ftable.GetNsig0();j++){
	real sigeff=siginf2*secf.GetF(i,g,0,0.,temp_sec,ftable.GetSig0(j));
	real f=sigeff/siginf;
	/*
	if(i==1){
	  cout<<g<<" "<<j<<" "<<siginf<<" "<<siginf2<<" "<<f<<"\n";
	};
	*/
        ftab_new.PutData(0,num_temp,j,f);
      };

      ftable.PutFData(i,g,0,ftab_new);
    };

    };
  };

  for(int i=0;i<ftable.GetNomtft();i++){    
    int org=ftable.GetNoTemp(i);
    ftable.PutNoTemp(i,org+1);
  };
  ftable.AddMaxtemp();
  ftable.AddValTemp(temp_sec);
  ftable.AddValTempLog(log(temp_sec));
  
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++ For FRENDY-MG

void LibData::ReadFileFrendyMG(string frendy_dir, string inpname, int matid, int group_inp)
{
  // ... Abondoned? ...
  
  mat=matid;

  PutGroup(group_inp);
  ftable.PutGroup(group_inp);

  int pl=5;
 
  exist_fiss=true;
  string filename_nuchi=frendy_dir+"FMNuChi_"+inpname+".txt";
  {
    ifstream fin;
    fin.open(filename_nuchi.data(),ios::in);
    if(fin.fail())exist_fiss=false;
  };
  

  int nomt2d=3;
  SetSizeMaxpl1(nomt2d);

  for(int i=0;i<nomt2d;i++){
    PutMaxpl1(i,pl+1);
    xsdata.PutDim2d(i,pl+1);
  };


  // ... f-table initialization

  int nomtft=5;
  ftable.PutNomtft(nomtft); 

  int num_sigb=9;
  real sigb[]={1e10,1e6,1e5,1e4,1e3,1e2,1e1,1,0.1};
  ftable.PutNsig0(num_sigb);
  for(int i=0;i<num_sigb;i++){
    //ftable.PutValSig0(i,sigb[i]);
    ftable.PutValSig0(i,sigb[num_sigb-i-1]);
  };

  int num_temp=1;
  real temp[]={600.};
  ftable.PutMaxtemp(num_temp);
  for(int i=0;i<num_temp;i++){
    ftable.PutValTemp(i,temp[i]);
  };

  ftable.PutMaxnr(1);
  for(int j=0;j<num_sigb;j++){
    ftable.PutValRpara(0,j,0.);
  };

  int mtlist[]={18,102,2,1,998};
  for(int i=0;i<nomtft;i++){
    ftable.PutNoTemp(i,num_temp);
    ftable.PutNoRpara(i,1);
    int grpstt=0;
    if(i==0&&!exist_fiss)grpstt=group;
    ftable.PutStartGrpFtab(i,grpstt);
    int grpend=group-1;
    if(i==4)grpend-=1;
    ftable.PutEndGrpFtab(i,grpend);
    ftable.PutMtFtab(i,mtlist[i]);
  };
  ftable.PutStartGrpAllMF(0);
  ftable.PutEndGrpAllMF(group-1);

  ftable.Initialize();

  // ... inelastic cheking

  int max_inelalevel=51;
  for(int i=51;i<=90;i++){
    string filename=frendy_dir+"FMMicXS2d_"+inpname+"_MT"+IntToString(i)+"_bg0.txt";
    ifstream fin;
    fin.open(filename.data(),ios::in);
    if(fin.fail()){
      i=100;
    }else{
      max_inelalevel=i;
    };
  };
  cout<<"# max_inelalevel : "<<max_inelalevel<<"\n";


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ... 1D-XS

  vector<real> data_chi;
  vector<real> data_chid;
  vector<real> nud;
  if(exist_fiss)ReadFileFrendyMG_NuChi(filename_nuchi, data_chi,data_chid,nud);

  int num_1dxs=8;
  if(!exist_fiss)num_1dxs=6;
  int list_1dxs[]={-1,1,2,102,4,16,18,452};
  vector<real> ebnd;
  vector< vector< vector<real> > > data1d(num_1dxs);  

#if 0
  for(int i=0;i<num_1dxs;i++){
    string filename=frendy_dir+"FMMicXS1d_"+inpname+"_MT"+IntToString(list_1dxs[i])+".txt";
    ReadFileFrendyMG_XS1D(filename, data1d[i], ebnd);
  };
#endif

#if 1
  for(int i=0;i<num_1dxs;i++){
    if(i!=4){ // (inela)
      string filename=frendy_dir+"FMMicXS1d_"+inpname+"_MT"+IntToString(list_1dxs[i])+".txt";
      ReadFileFrendyMG_XS1D(filename, data1d[i], ebnd);
    };
  };

  // (inela)
  vector<real> data1d_ine;
  {
  string filename1=frendy_dir+"FMMicXS1d_"+inpname+"_MT";
  string filename2=".txt";
  ReadFileFrendyMG_XS1D_MT4(filename1, filename2, data1d_ine, max_inelalevel);
  };
#endif

  // (siginf)
  xsdata.GetData1d(sigt).put_data(data1d[1][0]); 
  xsdata.GetData1d(sigel).put_data(data1d[2][0]); 
  xsdata.GetData1d(sigc).put_data(data1d[3][0]); 
  //xsdata.GetData1d(siginel).put_data(data1d[4][0]); 
  xsdata.GetData1d(siginel).put_data(data1d_ine); 
  xsdata.GetData1d(sign2n).put_data(data1d[5][0]); 
  if(exist_fiss){
    xsdata.GetData1d(sigf).put_data(data1d[6][0]); 
    xsdata.GetData1d(nu).put_data(data1d[7][0]); 
    xsdata.GetData1d(chi).put_data(data_chi); 
  };
    
  // (f-table w/o elastic removal)
  int list_ftable[]={6,3,2,0,-1};
  int istt=0;
  if(!exist_fiss)istt=1;
  for(int i=istt;i<nomtft;i++){
    int tt=list_ftable[i];
    if(tt!=-1){
      for(int g=0;g<group;g++){
	for(int n=0;n<num_sigb;n++){
	  real f_factor=1.;
	  if(data1d[tt][0][g]>0.)f_factor=data1d[tt][n][g]/data1d[tt][0][g];
	  //ftable.PutFData(i,g,0,0,0,n,f_factor);
	  ftable.PutFData(i,g,0,0,0,num_sigb-n-1,f_factor);
	};
      };
    };
  };
  

#if 0
  for(int i=0;i<group;i++){
    cout<<ebnd[i]<<" ";
    for(int j=0;j<num_1dxs;j++){
      cout<<data1d[j][0][i]<<" ";
    };
    cout<<"\n";
  };
#endif

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ... 2D-XS

  if(exist_fiss){

  // ... chi
  vector< vector< vector<real> > > chi2d;
  string filename_chi=frendy_dir+"FMMicXS2d_"+inpname+"_MT18_bg0.txt";
  ReadFileFrendyMG_XS2D(filename_chi, chi2d, 0);

  // (adding delayed neutron contribution)
  for(int g=0;g<group;g++){
    real nud_sigf=nud[g]*data1d[6][0][g];
    for(int g2=0;g2<group;g2++){
      chi2d[0][g2][g]+=data_chid[g2]*nud_sigf;
    };
  };

#if 0
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      sum+=chi2d[0][g2][g];
    };
    real org=data1d[7][0][g]*data1d[6][0][g];
    cout<<g<<" "<<sum<<" "<<org<<" "<<sum/org<<"\n";
  };
  exit(0);
#endif

  NormalizeMatrixData(chi2d[0]);
  exist_chi_vector=true;

  chiv.PutNoVector(group);
  for(int g=0;g<group;g++){
    chiv.PutVectorID(g,g);
    for(int g2=0;g2<group;g2++){
      chiv.GetChiVector(g).put_data(g2,chi2d[0][g2][g]);
    };
  };

# if 0
  for(int g=0;g<group;g++){
    cout<<ebnd[g]<<" ";
    real letwid=log(ebnd[g]/ebnd[g+1]);
    for(int g2=0;g2<group;g2++){
      cout<<chi2d[0][g][g2]/letwid<<" ";
    };
    cout<<"\n";
  };
#endif

  };

  int num_2dxs=3;
  int list_2dxs[]={2,16,4};
  vector< vector< vector< vector< vector<real> > > > > data2d(num_2dxs);
  for(int i=0;i<num_2dxs;i++){
    data2d[i].resize(num_sigb);
  };

  // ... (n,n) and (n,2n)
  for(int i=0;i<num_2dxs;i++){
    //for(int i=0;i<num_2dxs-1;i++){
    for(int j=0;j<num_sigb;j++){
      string filename=frendy_dir+"FMMicXS2d_"+inpname+"_MT"+IntToString(list_2dxs[i])+"_bg"+IntToString(j)+".txt";
      ReadFileFrendyMG_XS2D(filename, data2d[i][j], pl);
    };
  };

  // ... mu-bar calculation
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      sum+=data2d[0][0][1][g2][g];
    };
    real val=sum/data1d[2][0][g];
    xsdata.GetData1d(mu).put_data(g,val);
    //real ref=xslib.GetLibData(922350).GetXSData().GetData1d(mu).get_dat(g);
    //cout<<ebnd[g]<<" "<<val<<" "<<ref<<" "<<val/ref<<"\n";
    //cout<<ebnd[g]<<" "<<val<<"\n";
  };

  // ... f-table for elastic removal
  for(int g=0;g<group-1;g++){
    for(int n=0;n<num_sigb;n++){
      real sum0=0.;
      real sum1=0.;
      for(int g2=g+1;g2<group;g2++){
	sum0+=data2d[0][0][0][g2][g];
	sum1+=data2d[0][n][0][g2][g];
      };
      real f_factor=sum1/sum0;
      //ftable.PutFData(4,g,0,0,0,n,f_factor);
      ftable.PutFData(4,g,0,0,0,num_sigb-n-1,f_factor);
    };
  };

#if 0
  // ... (n,n')
  num_sigb=1;

  string filename1=frendy_dir+"FMMicXS2d_"+inpname+"_MT";
  for(int i=0;i<num_sigb;i++){
    string filename2="_bg"+IntToString(i)+".txt";
    ReadFileFrendyMG_XS2D_MT4(filename1, filename2, data2d[num_2dxs-1][i], pl, max_inelalevel);
  };
#endif

  int mt_order[]={0,2,1}; // (n,n), (n,n'), (n,2n)
  for(int i=0;i<num_2dxs;i++){
    for(int l=0;l<=pl;l++){
      real factor=2.*l+1;
      for(int g=0;g<group;g++){
	for(int g2=0;g2<group;g2++){
          xsdata.GetData2d(mt_order[i],l).put_data(g,g2,data2d[i][0][l][g2][g]*factor);
	};
      };
    };
  };


#if 0
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      //sum+=data2d[0][0][g2][g];
      //sum+=data2d[1][0][g2][g]; // (n,2n)
      sum+=data2d[2][0][0][g2][g];  // (n,n')
    };
    //cout<<ebnd[g]<<" "<<sum<<" "<<data1d[2][g][0]<<"\n";
    //cout<<ebnd[g]<<" "<<sum<<" "<<data1d[6][g][0]<<" "<<sum/data1d[6][g][0]<<"\n"; // (n,2n)
    cout<<ebnd[g]<<" "<<sum<<" "<<data1d[5][0][g]<<" "<<sum/data1d[5][0][g]<<"\n"; // (n,n')
  };
#endif
  

};

void LibData::ReadFileFrendyMG(int tempnum, int *temp_list, string frendy_dir, string inpname_i, int matid, int group_inp)
{
  MATIDTranslator midt;
  mat=matid;

  PutGroup(group_inp);
  ftable.PutGroup(group_inp);

  int pl=5;

  vector<string> inpname(tempnum);
  for(int i=0;i<tempnum;i++){
    string tempname=IntToString(temp_list[i]);
    if(tempname.size()==3)tempname="0"+tempname;
    //inpname[i]="TMP_"+inpname_i;
    inpname[i]=midt.Name(matid)+"_"+tempname+"K_"+inpname_i;    
  };
 
  exist_fiss=true;
  string filename_nuchi=frendy_dir+"FMNuChi_"+inpname[0]+".txt";
  {
    ifstream fin;
    fin.open(filename_nuchi.data(),ios::in);
    if(fin.fail())exist_fiss=false;
    fin.close();
  };

  bool exist_n2n=true;
  string filename_n2n=frendy_dir+"FMMicXS1d_"+inpname[0]+"_MT16.txt";
  {
    ifstream fin;
    fin.open(filename_n2n.data(),ios::in);
    if(fin.fail())exist_n2n=false;
    fin.close();
  };

  int nomt2d=3;
  SetSizeMaxpl1(nomt2d);

  for(int i=0;i<nomt2d;i++){
    PutMaxpl1(i,pl+1);
    xsdata.PutDim2d(i,pl+1);
  };


  // ... f-table initialization

  int nomtft=5;
  ftable.PutNomtft(nomtft); 


  int num_sigb=9;
  real sigb[]={1e10,1e6,1e5,1e4,1e3,1e2,1e1,1,0.1};

  /*
  int num_sigb=2;
  real sigb[]={1e10,1e-10};
  */


  ftable.PutNsig0(num_sigb);
  for(int i=0;i<num_sigb;i++){
    ftable.PutValSig0(i,sigb[num_sigb-i-1]);
  };

  int num_temp=tempnum;
  ftable.PutMaxtemp(num_temp);
  for(int i=0;i<num_temp;i++){
    ftable.PutValTemp(i,real(temp_list[i]));
  };

  ftable.PutMaxnr(1);
  for(int j=0;j<num_sigb;j++){
    ftable.PutValRpara(0,j,0.);
  };

  int mtlist[]={18,102,2,1,998};
  for(int i=0;i<nomtft;i++){
    ftable.PutNoTemp(i,num_temp);
    ftable.PutNoRpara(i,1);
    int grpstt=0;
    if(i==0&&!exist_fiss)grpstt=group;
    ftable.PutStartGrpFtab(i,grpstt);
    int grpend=group-1;
    if(i==4)grpend-=1;
    ftable.PutEndGrpFtab(i,grpend);
    ftable.PutMtFtab(i,mtlist[i]);
  };
  ftable.PutStartGrpAllMF(0);
  ftable.PutEndGrpAllMF(group-1);

  ftable.Initialize();

  // ... inelastic checking
  int max_inelalevel=50;
  for(int i=51;i<=90;i++){
    string filename=frendy_dir+"FMMicXS2d_"+inpname[0]+"_MT"+IntToString(i)+"_bg0.txt";
    ifstream fin;
    fin.open(filename.data(),ios::in);
    if(fin.fail()){
      i=100;
    }else{
      max_inelalevel=i;
    };
    fin.close();
  };
  if(max_inelalevel==50){
    cout<<"# No inelastic scattering data.\n";
  }else{
    cout<<"# max_inelalevel : "<<max_inelalevel<<"\n";
  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ... 1D-XS

  vector<real> data_chi;
  vector<real> data_chid;
  vector<real> nud;
  if(exist_fiss)ReadFileFrendyMG_NuChi(filename_nuchi, data_chi,data_chid,nud);

  int num_1dxs=8;
  if(!exist_fiss)num_1dxs=6;
  int list_1dxs[]={-1,1,2,102,4,16,18,452};
  vector<real> ebnd;
  vector< vector< vector< vector<real> > > > data1d(num_1dxs);  
  // [reaction][temperature][sig0][group]

  for(int i=0;i<num_1dxs;i++){
    data1d[i].resize(num_temp);
    if(i!=4&&!(i==5&&!exist_n2n)){ // (non-inela, n2n)
      for(int j=0;j<num_temp;j++){
        string filename=frendy_dir+"FMMicXS1d_"+inpname[j]+"_MT"+IntToString(list_1dxs[i])+".txt";
        ReadFileFrendyMG_XS1D(filename, data1d[i][j], ebnd);
      };
    };
  };

  // (inela)
  vector<real> data1d_ine;
  if(max_inelalevel>50){
    string filename1=frendy_dir+"FMMicXS1d_"+inpname[0]+"_MT";
    string filename2=".txt";
    ReadFileFrendyMG_XS1D_MT4(filename1, filename2, data1d_ine, max_inelalevel);
  };

  // (capture)
  for(int j=0;j<num_temp;j++){
    string filename1=frendy_dir+"FMMicXS1d_"+inpname[j]+"_MT";
    string filename2=".txt";
    ReadFileFrendyMG_XS1D_MT102(filename1, filename2, data1d[3][j]);
  };
	
  // (siginf)
  xsdata.GetData1d(sigt).put_data(data1d[1][0][0]); 
  xsdata.GetData1d(sigel).put_data(data1d[2][0][0]); 
  xsdata.GetData1d(sigc).put_data(data1d[3][0][0]); 
  if(max_inelalevel>50)xsdata.GetData1d(siginel).put_data(data1d_ine); 
  if(exist_n2n)xsdata.GetData1d(sign2n).put_data(data1d[5][0][0]); 
  if(exist_fiss){
    xsdata.GetData1d(sigf).put_data(data1d[6][0][0]); 
    xsdata.GetData1d(nu).put_data(data1d[7][0][0]); 
    xsdata.GetData1d(chi).put_data(data_chi); 
  };
    
  // (f-table w/o elastic removal)
  int list_ftable[]={6,3,2,0,-1}; // 18,102,2,-1
  int istt=0;
  if(!exist_fiss)istt=1;
  for(int i=istt;i<nomtft;i++){
    int tt=list_ftable[i];
    if(tt!=-1){
      for(int g=0;g<group;g++){
	for(int t=0;t<num_temp;t++){
  	  for(int n=0;n<num_sigb;n++){
	    real f_factor=1.;
	    if(data1d[tt][0][0][g]>0.)f_factor=data1d[tt][t][n][g]/data1d[tt][0][0][g];
	    //ftable.PutFData(i,g,0,0,0,n,f_factor);
	    ftable.PutFData(i,g,0,0,t,num_sigb-n-1,f_factor);
	  };
	};
      };
    };
  };
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ... 2D-XS

  if(exist_fiss){

  // ... chi
  vector< vector< vector<real> > > chi2d;
  string filename_chi=frendy_dir+"FMMicXS2d_"+inpname[0]+"_MT18_bg0.txt";
  ReadFileFrendyMG_XS2D(filename_chi, chi2d, 0);

  // (adding delayed neutron contribution)
  for(int g=0;g<group;g++){
    real nud_sigf=nud[g]*data1d[6][0][0][g];
    for(int g2=0;g2<group;g2++){
      chi2d[0][g2][g]+=data_chid[g2]*nud_sigf;
    };
  };

#if 0
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      sum+=chi2d[0][g2][g];
    };
    real org=data1d[7][0][0][g]*data1d[6][0][0][g];
    cout<<g<<" "<<sum<<" "<<org<<" "<<sum/org<<"\n";
  };
  exit(0);
#endif

  NormalizeMatrixData(chi2d[0]);
  exist_chi_vector=true;

  chiv.PutNoVector(group);
  for(int g=0;g<group;g++){
    chiv.PutVectorID(g,g);
    for(int g2=0;g2<group;g2++){
      chiv.GetChiVector(g).put_data(g2,chi2d[0][g2][g]);
    };
  };

# if 0
  for(int g=0;g<group;g++){
    cout<<ebnd[g]<<" ";
    real letwid=log(ebnd[g]/ebnd[g+1]);
    for(int g2=0;g2<group;g2++){
      cout<<chi2d[0][g][g2]/letwid<<" ";
    };
    cout<<"\n";
  };
#endif

  };

  int num_2dxs=3;
  int list_2dxs[]={2,16,4};
  vector< vector< vector< vector< vector< vector<real> > > > > > data2d(num_2dxs);
  // [reaction][temperature][sig0][pl][gin][gout]
  for(int i=0;i<num_2dxs;i++){
    data2d[i].resize(num_temp);
  };

  // ... (n,n), (n,2n) and (n,n')
  {
    int i=0;
    for(int t=0;t<num_temp;t++){
      data2d[i][t].resize(num_sigb);
      for(int j=0;j<num_sigb;j++){
        string filename=frendy_dir+"FMMicXS2d_"+inpname[t]+"_MT"+IntToString(list_2dxs[i])+"_bg"+IntToString(j)+".txt";
        ReadFileFrendyMG_XS2D(filename, data2d[i][t][j], pl);
      };
    };
  };
  if(exist_n2n){
    int i=1;
    for(int t=0;t<num_temp;t++){
      data2d[i][t].resize(num_sigb);
      for(int j=0;j<num_sigb;j++){
        string filename=frendy_dir+"FMMicXS2d_"+inpname[t]+"_MT"+IntToString(list_2dxs[i])+"_bg"+IntToString(j)+".txt";
        ReadFileFrendyMG_XS2D(filename, data2d[i][t][j], pl);
      };
    };
  };
  if(max_inelalevel>50){
    int i=2;
    for(int t=0;t<num_temp;t++){
      data2d[i][t].resize(num_sigb);
      for(int j=0;j<num_sigb;j++){
        string filename=frendy_dir+"FMMicXS2d_"+inpname[t]+"_MT"+IntToString(list_2dxs[i])+"_bg"+IntToString(j)+".txt";
        ReadFileFrendyMG_XS2D(filename, data2d[i][t][j], pl);
      };
    };
  };
    
  // ... mu-bar calculation
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      sum+=data2d[0][0][0][1][g2][g];
    };
    real val=sum/data1d[2][0][0][g];
    xsdata.GetData1d(mu).put_data(g,val);
    //real ref=xslib.GetLibData(922350).GetXSData().GetData1d(mu).get_dat(g);
    //cout<<ebnd[g]<<" "<<val<<" "<<ref<<" "<<val/ref<<"\n";
    //cout<<ebnd[g]<<" "<<val<<"\n";
  };

  // ... f-table for elastic removal
  for(int g=0;g<group-1;g++){
    for(int t=0;t<num_temp;t++){
      for(int n=0;n<num_sigb;n++){
        real sum0=0.;
        real sum1=0.;
        for(int g2=g+1;g2<group;g2++){
    	  sum0+=data2d[0][0][0][0][g2][g];
 	  sum1+=data2d[0][t][n][0][g2][g];
        };
        real f_factor=sum1/sum0;
        //ftable.PutFData(4,g,0,0,0,n,f_factor);
        ftable.PutFData(4,g,0,0,t,num_sigb-n-1,f_factor);
      };
    };
  };

  int mt_order[]={0,2,1}; // (n,n), (n,n'), (n,2n)
  {
    int i=0;
    for(int l=0;l<=pl;l++){
      real factor=2.*l+1;
      for(int g=0;g<group;g++){
	for(int g2=0;g2<group;g2++){
          xsdata.GetData2d(mt_order[i],l).put_data(g,g2,data2d[i][0][0][l][g2][g]*factor);
	};
      };
    };
  };
  if(exist_n2n){
    int i=1;
    for(int l=0;l<=pl;l++){
      real factor=2.*l+1;
      for(int g=0;g<group;g++){
	for(int g2=0;g2<group;g2++){
          xsdata.GetData2d(mt_order[i],l).put_data(g,g2,data2d[i][0][0][l][g2][g]*factor);
	};
      };
    };
  };
  if(max_inelalevel>50){
    int i=2;
    for(int l=0;l<=pl;l++){
      real factor=2.*l+1;
      for(int g=0;g<group;g++){
	for(int g2=0;g2<group;g2++){
          xsdata.GetData2d(mt_order[i],l).put_data(g,g2,data2d[i][0][0][l][g2][g]*factor);
	};
      };
    };
  };



};

void LibData::ReadFileFrendyMG_NuChi(string filename, vector<real> &data, vector<real> &data_chid, vector<real> &nud)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<filename<<"\n";
    exit(0);
  };

  data.resize(group);
  data_chid.resize(group);
  nud.resize(group);

  string dummy;
  /*
  for(int i=0;i<20;i++){
    fin>>dummy;
  };
  */
  for(int i=0;i<9;i++){
    fin>>dummy;
  };
  real dnf;
  fin>>dnf;
  for(int i=0;i<10;i++){
    fin>>dummy;
  };

  vector<real> nu_d_frac;

  int nd_group=6;
  
  if(dnf>0.){
    
  real sum=0.;
  for(int i=0;i<6;i++){
    real tmp;
    fin>>tmp;
    nu_d_frac.push_back(tmp);
    sum+=tmp;
  };

  if(sum<0.99){
    nd_group=8;
    real tmp;
    fin>>tmp;
    nu_d_frac.push_back(tmp);
    fin>>tmp;
    nu_d_frac.push_back(tmp);
  };

  }else{
    cout<<"#    ..... no delayed neutron data .....\n";
    nd_group=0;
  };

  //for(int i=0;i<27;i++){
  for(int i=0;i<21+nd_group;i++){    
    fin>>dummy;
  };

  vector< vector<real> > chi_d_frac(nd_group);
  for(int i=0;i<nd_group;i++){
    chi_d_frac[i].resize(group);
  };

  for(int g=0;g<group;g++){
    for(int i=0;i<6;i++){
      fin>>dummy;
    };
    fin>>nud[g];
    fin>>data[g];
    fin>>dummy;
    fin>>data_chid[g];
    for(int i=0;i<nd_group;i++){
      //fin>>dummy;
      fin>>chi_d_frac[i][g];
    };
  };

#if 0
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int i=0;i<6;i++){
      sum+=chi_d_frac[i][g]*nu_d_frac[i];
    };
    cout<<g<<" "<<sum<<" "<<data_chid[g]<<"\n";
  };
#endif

  fin.close();
};
 
void LibData::ReadFileFrendyMG_XS1D(string filename, vector< vector<real> > &data, vector<real> &ebnd)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<filename<<"\n";
    exit(0);
  };

  int nsig0=ftable.GetNsig0();
  int maxtemp=ftable.GetMaxtemp();

  string dummy;
  for(int i=0;i<7;i++){
    fin>>dummy;
  };

  // sig-0
  for(int i=0;i<nsig0;i++){
    fin>>dummy;
  };

  // temp
  fin>>dummy;
  for(int i=0;i<nsig0;i++){
    fin>>dummy;
  };

  data.resize(nsig0);
  for(int i=0;i<nsig0;i++){
    data[i].resize(group);
  };

  ebnd.resize(group+1);

  int dummy_i;
  real dummy_r;
  for(int g=0;g<group;g++){
    fin>>dummy_i;
    fin>>ebnd[g];
    fin>>dummy_r;
    fin>>dummy_r;
    if(g==group-1)ebnd[g+1]=dummy_r;
    for(int i=0;i<nsig0;i++){
      fin>>data[i][g];
      //fin>>dummy_r;
    };
  };

  fin.close();
};
 
void LibData::ReadFileFrendyMG_XS2D(string filename, vector< vector< vector<real> > > &data, int pl)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<filename<<"\n";
    exit(0);
  };

  string dummy;
  for(int i=0;i<7;i++){
    fin>>dummy;
  };

  data.resize(pl+1);
  for(int i=0;i<=pl;i++){
    data[i].resize(group);
    fin>>dummy; // PL
    fin>>dummy; // =
    fin>>dummy; // l
    for(int g1=0;g1<group;g1++){
      data[i][g1].resize(group);
      fin>>dummy;
      for(int g2=0;g2<group;g2++){
	fin>>data[i][g1][g2];
	//if(data[i][g1][g2]>0.&&g1<g2)cout<<filename<<" "<<i<<" "<<g1<<" "<<g2<<" "<<data[i][g1][g2]<<"\n";
      };
    };
  };

  fin.close();
};

void LibData::ReadFileFrendyMG_XS2D_MT4(string filename1, string filename2, vector< vector< vector<real> > > &data, int pl, int max_inelalevel)
{
  bool exist_mt91=true;
  {
    ifstream fin;
    string tmp=filename1+IntToString(91)+filename2;
    fin.open(tmp.data(),ios::in);
    if(fin.fail())exist_mt91=false;
    fin.close();
  };
  
  data.resize(pl+1);
  for(int i=0;i<=pl;i++){
    data[i].resize(group);
    for(int g=0;g<group;g++){
      data[i][g].resize(group);
      for(int g2=0;g2<group;g2++){
	data[i][g][g2]=0.;
      };
    };
  };

  for(int mt=51;mt<=91;mt++){
    if(mt<=max_inelalevel||(mt==91&&exist_mt91)){
    vector< vector< vector<real> > > tmp;
    ReadFileFrendyMG_XS2D(filename1+IntToString(mt)+filename2,tmp,pl);
    for(int i=0;i<=pl;i++){
      for(int g=0;g<group;g++){
	for(int g2=0;g2<group;g2++){
	  data[i][g][g2]+=tmp[i][g][g2];
	};
      };
    };
    };
  };
};

void LibData::ReadFileFrendyMG_XS1D_MT4(string filename1, string filename2, vector<real> &data, int max_inelalevel)
{
  bool exist_mt91=true;
  {
    ifstream fin;
    string tmp=filename1+IntToString(91)+filename2;
    fin.open(tmp.data(),ios::in);
    if(fin.fail())exist_mt91=false;
    fin.close();
  };
  
  
  data.resize(group);
  for(int g=0;g<group;g++){
    data[g]=0.;
  };

  for(int mt=51;mt<=91;mt++){
    if(mt<=max_inelalevel||(mt==91&&exist_mt91)){    
      vector< vector<real> >tmp;
      vector<real> tmp2;
      ReadFileFrendyMG_XS1D(filename1+IntToString(mt)+filename2,tmp,tmp2);
      for(int g=0;g<group;g++){
        data[g]+=tmp[0][g];
      };
    };
  };
};

void LibData::ReadFileFrendyMG_XS1D_MT102(string filename1, string filename2, vector< vector<real> > &data)
{
#if 0  

(capture-correction for MF=103-107 and 113)
   103  (n,p)    "np"
   104  (n,d)    "nd"
   105  (n,t)    "nt"
   106  (n,3He)  "nh"?
   107  (n,a)    "na"
   108  (n,2a)
   109  (n,3a)
   111  (n,2p)   "n2p"
   112  (n,pa)   
   113  (n,t2a)  "nt2a"
   114  (n,d2a)
   115  (n,pd)
   116  (n,pt)
   117  (n,da)

#endif
  
  for(int mt=103;mt<=113;mt++){
    if(mt<=107||mt==111||mt==113){
    ifstream fin;
    string tmp=filename1+IntToString(mt)+filename2;
    fin.open(tmp.data(),ios::in);
    if(!fin.fail()){
      vector< vector<real> >tmp;
      vector<real> tmp2;
      ReadFileFrendyMG_XS1D(filename1+IntToString(mt)+filename2,tmp,tmp2);
      for(int g=0;g<group;g++){
	for(int n=0;n<data.size();n++){	
          data[n][g]+=tmp[n][g];
	};
      };
    };
    };
  };
};

void LibData::NormalizeMatrixData(vector< vector<real> > &data)
{
  for(int g=0;g<group;g++){
    real sum=0.;
    for(int g2=0;g2<group;g2++){
      sum+=data[g2][g];
    };
    for(int g2=0;g2<group;g2++){
      data[g2][g]/=sum;
    };
  };
};

void LibData::TSL(string tslname, string aceid, int mtnum, int *mtlist, string libdir, string libname)
{
  PutGroup(211);
  
  string frendy_dir="/home/roko-de-go/CODE/frendy-MG/run_tsl/";

  int pl=5;
  vector< vector< vector< vector<real> > > > data(mtnum);

  for(int i=0;i<mtnum;i++){

    string filename=frendy_dir+"FMMicXS2d_"+tslname+"_300K_"+aceid+"_MT"+IntToString(mtlist[i])+"_bg0.txt";
    ReadFileFrendyMG_XS2D(filename, data[i], pl);

  };

  int grp_nonzero_xs=group;
  for(int j=0;j<mtnum;j++){
  for(int i=0;i<=pl;i++){
    for(int g=0;g<group;g++){
      for(int g2=0;g2<group;g2++){
	if(data[j][i][g2][g]!=0.){
	  if(g<grp_nonzero_xs)grp_nonzero_xs=g;
	  if(j!=0){
	    data[0][i][g2][g]+=data[j][i][g2][g];
	  };
	};
      };
    };
  };
  };

  //cout<<"# Beginning : "<<grp_nonzero_xs<<"\n";

  //string mdir="/home/roko-de-go/CBGLIB/j4.211g.20210301/Thermal/tsl_tmp";
  string mdir=libdir+libname;
  ofstream fout;
  
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };

  fout.setf(ios::scientific);
  fout.precision(8);

  fout<<group<<"\n";

  fout<<"   1\n"; // Number of temperature points
  fout<<"   300.\n"; // Temperature
  
  fout<<"   5\n"; // Legendre order
  fout<<"   "<<grp_nonzero_xs+1<<"\n";

  // For higher-order Pl components of thermal scattering cross section,
  // (2l+1) is NOT multiplied to the library,
  // and this factor is multiplied in the cross section generation phase.
  // Please see [OnePointCalculator::CalThermalScatteringMatrix].
  
  for(int i=0;i<=pl;i++){
    //real factor=2.*i+1;    
    for(int j=grp_nonzero_xs;j<group;j++){
      for(int k=grp_nonzero_xs;k<group;k++){
	//fout<<"  "<<data[0][i][k][j]*factor<<"\n";
	fout<<"  "<<data[0][i][k][j]<<"\n";	
      };
    };
  };

  fout.close();
};

// +++ XSLibrary

XSLibrary::XSLibrary()
{
  warning_print=false;
};

bool XSLibrary::ExistLibData(int mat)
{
  map<int,LibData>::iterator i=data.find(mat);
  if(i==this->data.end()){
    return false;
  }else{
    return true;
  };
};

LibData& XSLibrary::GetLibData(int mat)
{
  map<int,LibData>::iterator i=data.find(mat);
  if(i==this->data.end()){
    cout<<"# You cannot find LibData in XSLibrary.\n";
    cout<<"# MAT number is "<<mat<<"\n";
    exit(0);
  }else{
    return i->second;
  };
};

void XSLibrary::ReadFile(int nucnum,string mdir,string *filename,int *matno)
{
  //LibData tmp;
  for(int i=0;i<nucnum;i++){
    LibData tmp;
    tmp.ReadFile(mdir,filename[i],midt);
    AddLibData(matno[i],tmp);
  };
};

void XSLibrary::ReadFile(int nucnum,string mdir,string *filename)
{
  for(int i=0;i<nucnum;i++){
    LibData tmp;
    tmp.ReadFile(mdir,filename[i],midt);
    AddLibData(tmp.GetMat(),tmp);
  };
};

void XSLibrary::ReadPTFile(string mdir,string filename,int matno,real maxerr)
{
  GetLibData(matno).GetPtable().ReadFile(mdir,filename,maxerr);
};

void XSLibrary::ReadPTMatrixFile(string mdir,string filename,int matno)
{
  GetLibData(matno).GetPtable().ReadMatrixFile(mdir,filename);
};

void XSLibrary::ReadNEnergy(string mdir,string ss)
{
  cout<<"# Reading N-ENERGY data at the directory ["<<mdir<<"]\n";
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# (file name) "<<mdir<<"\n";
    exit(0);
  };
  int grp;
  fin>>grp;
  enband.put_imax(grp+1);
  wtflux.put_imax(grp);
  for(int i=0;i<grp+1;i++){
    real tmp;
    fin>>tmp;
    enband.put_data(i,tmp);
  };
  bool zero_flux=true;
  for(int i=0;i<grp;i++){
    real tmp;
    fin>>tmp;
    if(tmp>0.)zero_flux=false;
    wtflux.put_data(i,tmp);
  };
  fin.close();
  if(zero_flux){
    cout<<"Weight function is zero in the [N-ENERGY] file.\n";
    exit(0);
  };
};

void XSLibrary::PutBellFactor(int mat,real *inp)
{
  GetLibData(mat).PutBellFactor(inp);
};

void XSLibrary::ReadBellFactor(string mdir,string ss,int mat)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"(file name) "<<mdir<<"\n";
    exit(0);
  };
  int grp;
  fin>>grp;
  real *inp=new real[grp];
  for(int i=0;i<grp;i++){
    real tmp;
    fin>>tmp;
    inp[i]=tmp;
  };
  PutBellFactor(mat,inp);  
  fin.close();
  delete [] inp;
};

void XSLibrary::ShowSelf()
{
  map<int,LibData>::iterator it=data.begin();
  while(it!=data.end()){
    cout<<(*it).second.GetMat()<<"\n";
    it++;
  };
};

void XSLibrary::ShowCrossSectionData1D(int mat)
{
  if(ExistLibData(mat)){
    data[mat].GetXSData().ShowData1D(enband);
  };
};

void XSLibrary::Update1DChiDataForFastReactorCalculation()
{
  map<int,LibData>::iterator it=data.begin();
  while(it!=data.end()){
    (*it).second.Update1DChiDataForFastReactorCalculation();
    it++;
  };
};


LibData XSLibrary::LibraryMixing(int newnucid, int nucn, int *nucid, real *wgt)
{
  // non-zero check for weight
  for(int i=0;i<nucn;i++){
    if(wgt[i]<=0.)wgt[i]=1e-20;
  };

  int group=GetGroup();

  LibData isonuc;
  isonuc.PutGroup(group);
  isonuc.PutExistFiss(false);
  isonuc.PutMat(newnucid);

  // +++ Averaged "GroupDataSet" +++

  vector<GroupDataSet> cr(nucn);
  for(int i=0;i<nucn;i++){
    cr[i]=GetLibData(nucid[i]).GetXSData();
  };

  int maxid=-1;
  int maxnum=0;
  for(int i=0;i<nucn;i++){
    int tmp=0;
    for(int j=0;j<3;j++){
      tmp+=cr[i].GetDim2d(j);
    };
    if(tmp>maxnum){
      maxnum=tmp;
      maxid=i;
    };
  };
  if(maxid==-1){
    cout<<"# Error !\n";
    exit(0);
  };

  // mu averaging
  GroupData1D nume(group);
  GroupData1D denom(group);
  nume.set_zero();
  denom.set_zero();
  for(int i=0;i<nucn;i++){
    nume=nume+(cr[i].GetData1d(mu).mult(cr[i].GetData1d(sigel)))*wgt[i];
    denom=denom+(cr[i].GetData1d(sigel)*wgt[i]);
  };
  nume=nume/denom;

  // cross section averaging
  for(int i=0;i<nucn;i++){
    cr[i].Multiply(wgt[i]);
  };
  GroupDataSet gds=cr[maxid];
  for(int i=0;i<nucn;i++){
    if(i!=maxid)gds.Add(cr[i]);
  };

  // over-writing averaged mu
  gds.GetData1d(mu).copy(nume);

  // (original xsdata is recovered)
  for(int i=0;i<nucn;i++){
    cr[i].Multiply(1./wgt[i]);
  };

  isonuc.SetSizeMaxpl1(gds.GetNum2d());
  for(int i=0;i<cr[maxid].GetNum2d();i++){
    isonuc.PutMaxpl1(i,gds.GetDim2d(i));
  };

  //cout<<"# finish-1.\n";

  // +++ Averaged "LibDataFTable" +++

  LibDataFTable ftab;
  ftab.PutGroup(group);

  vector<LibDataFTable> crf(nucn);
  for(int i=0;i<nucn;i++){
    crf[i]=GetLibData(nucid[i]).GetFtable();
  };

  for(int i=1;i<nucn;i++){
    bool check=false;
    if(crf[i].GetNomtft()!=crf[0].GetNomtft())check=true;
    if(crf[i].GetNsig0()!=crf[0].GetNsig0())check=true;
    if(crf[i].GetMaxtemp()!=crf[0].GetMaxtemp())check=true;
    if(crf[i].GetMaxnr()!=crf[0].GetMaxnr())check=true;
    if(check){
      cout<<"# Nomtft/nsig0/maxtemp/maxnr is inconsistent.\n";
      exit(0);
    };
  };

  int nomtft=crf[0].GetNomtft();
  ftab.PutNomtft(nomtft);
  ftab.PutNsig0(crf[0].GetNsig0());
  ftab.PutMaxtemp(crf[0].GetMaxtemp());
  ftab.PutMaxnr(crf[0].GetMaxnr());

  for(int i=0;i<ftab.GetNsig0();i++){
    ftab.PutValSig0(i,crf[0].GetSig0(i));
  };
  for(int i=0;i<ftab.GetMaxtemp();i++){
    ftab.PutValTemp(i,crf[0].GetTemp(i));
  };

  for(int j=0;j<nomtft;j++){
    int nt=crf[0].GetNoTemp(j);
    int mt=crf[0].GetMtFtab(j);
    int nr=crf[0].GetNoRpara(j);
    for(int i=1;i<nucn;i++){
      bool check=false;
      if(nt!=crf[i].GetNoTemp(j))check=true;
      if(mt!=crf[i].GetMtFtab(j))check=true;
      if(nr!=crf[i].GetNoRpara(j))check=true;
      if(check){
	cout<<"# NoTemp/mtftab/NoRpara is inconsistent.\n";
	exit(0);
      };
    };
    ftab.PutNoTemp(j,nt);
    ftab.PutMtFtab(j,mt);
    ftab.PutNoRpara(j,nr);
  };

  for(int mt=0;mt<nomtft;mt++){
    int stt=crf[0].GetStartGrpFtab(mt);
    int edd=crf[0].GetEndGrpFtab(mt);
    if(stt==-1)stt=group-1;
    if(edd==-1)edd=0;
    for(int i=1;i<nucn;i++){    
      int st2=crf[i].GetStartGrpFtab(mt);
      int ed2=crf[i].GetEndGrpFtab(mt);
      if(st2<stt&&st2!=-1)stt=st2;
      if(ed2>edd&&ed2!=-1)edd=ed2;
    };
    if(stt>edd){
      stt=-1;
      edd=-1;
    };
    ftab.PutStartGrpFtab(mt,stt);
    ftab.PutEndGrpFtab(mt,edd);
  };

  ftab.Initialize();

  vector<real> sigtinf(nucn); 
  vector<real> sigt0(nucn);

  for(int mt=0;mt<nomtft;mt++){
    int stt=ftab.GetStartGrpFtab(mt);
    int edd=ftab.GetEndGrpFtab(mt);
    if(stt!=-1){
      int nt=ftab.GetNoTemp(mt);
      int ns0=ftab.GetNsig0();
      for(int g=stt;g<=edd;g++){
	for(int t=0;t<nt;t++){
          real temp=ftab.GetTemp(t);
	  for(int s0=0;s0<ns0;s0++){

	    //cout<<mt<<" "<<g<<" "<<t<<" "<<s0<<"\n";
	    real s0b=ftab.GetSig0(s0);

	    for(int i=0;i<nucn;i++){
	      sigtinf[i]=cr[i].GetData1d(sigt).get_dat(g);
              sigt0[i]=sigtinf[i]; // initial setting
	    };
	    // (total cross section iteration)
	    bool conv=false;
	    int iternum=0;
	    while(!conv){
	      conv=true;
  	      for(int i=0;i<nucn;i++){
		real s0v=0.;
		for(int j=0;j<nucn;j++){
		  if(j!=i)s0v+=wgt[j]*sigt0[j];
		};
		s0v+=s0b;
		s0v/=wgt[i];
                real fc=crf[i].GetF(1,g,0,0.,temp,s0v);
                real fe=crf[i].GetF(2,g,0,0.,temp,s0v);
                real fi=crf[i].GetF(5,g,0,0.,temp,s0v);
	        real xsc=cr[i].GetData1d(sigc).get_dat(g);
	        real xse=cr[i].GetData1d(sigel).get_dat(g);
	        real xsi=cr[i].GetData1d(siginel).get_dat(g);
                real ft=(fc*xsc+fe*xse+fi*xsi)/(xsc+xse+xsi);
                real xst=ft*sigtinf[i];
		//if(g==19&&t==0)cout<<s0b<<" "<<xst<<" "<<s0v<<"\n";
	        if(fabs(xst/sigt0[i]-1.)>1e-5)conv=false;
		sigt0[i]=xst;
		//cout<<"# Iteration : "<<iternum<<"\n";
		iternum++;
	      };
	    };

            real s0v=0.;
      	    for(int j=0;j<nucn;j++){
	      s0v+=wgt[j]*sigt0[j];
	    };

	    real nume=0.;
	    for(int i=0;i<nucn;i++){
	      real f=crf[i].GetF(mt,g,0,0.,temp,(s0b+s0v)/wgt[i]-sigt0[i]);
	      //if(g==19&&t==0&&mt==3)cout<<(s0b+s0v)/wgt[i]-sigt0[i]<<" "<<f<<"\n";
              real xs=0.;
	      if(mt==0)xs=cr[i].GetData1d(sigf).get_dat(g);
	      if(mt==1)xs=cr[i].GetData1d(sigc).get_dat(g);
	      if(mt==2)xs=cr[i].GetData1d(sigel).get_dat(g);
	      if(mt==3)xs=cr[i].GetData1d(sigt).get_dat(g);
	      if(mt==4){
                xs=cr[i].GetData2d(sigel).get_sumx().get_dat(g)
		  -cr[i].GetData2d(sigel).get_dat(g,g);
	      };
	      if(mt==5)xs=cr[i].GetData1d(siginel).get_dat(g);
	      //if(g==19&&t==0&&s0==1&&mt==3)cout<<xs<<" "<<f<<" "<<wgt[i]<<" "<<nume<<"\n";
	      nume+=f*xs*wgt[i];
	    };
            real xs=0.;
            if(mt==0)xs=gds.GetData1d(sigf).get_dat(g);
	    if(mt==1)xs=gds.GetData1d(sigc).get_dat(g);
	    if(mt==2)xs=gds.GetData1d(sigel).get_dat(g);
	    if(mt==3)xs=gds.GetData1d(sigt).get_dat(g);
	    if(mt==4){
              xs=gds.GetData2d(sigel).get_sumx().get_dat(g)
                -gds.GetData2d(sigel).get_dat(g,g);
	    };
	    if(mt==5)xs=gds.GetData1d(siginel).get_dat(g);
	    //if(g==19&&t==0&&s0==1&&mt==3)cout<<nume<<" "<<xs<<"\n";
	    if(xs==0.){
	      nume=1.;
	    }else{
   	      nume/=xs;
	    };
	    //if(g==19&&t==0&&s0==1&&mt==3)cout<<mt<<" "<<g<<" "<<t<<" "<<s0<<" "<<nume<<"\n";
	    ftab.PutFData(mt,g,0,0,t,s0,nume);
	  };
	};
      };
    };
  };

  //cout<<ftab.GetF(3,19,0,0.,300.,1.)<<"\n";
  isonuc.GetXSData()=gds;
  isonuc.GetFtable()=ftab;

  //isonuc.GetFtable().WriteFtableWithSig0Temp(1,20,0);
  //crf[0].WriteFtableWithSig0Temp(1,20,0);
  //crf[1].WriteFtableWithSig0Temp(1,20,0);
  //crf[2].WriteFtableWithSig0Temp(1,20,0);
  //crf[3].WriteFtableWithSig0Temp(1,20,0);
  //lib.GetLibData(newnucid).GetFtable().WriteFtableWithSig0Temp(1,20,0);
  /*
  cout<<cr[0].GetData1d(sigt).get_dat(20)<<"\n";
  cout<<cr[1].GetData1d(sigt).get_dat(20)<<"\n";
  cout<<cr[2].GetData1d(sigt).get_dat(20)<<"\n";
  cout<<cr[3].GetData1d(sigt).get_dat(20)<<"\n";
  */
  //cr0.GetData1d(siginel).show_self();
  //gds.GetData1d(siginel).show_self();

  //isonuc.WriteFile(libdir,outfile);

  //delete [] nucid;
  //delete [] wgt;

  return isonuc;

};


void XSLibrary::ConstantCrossSectionData(int mat_inp, int g, real sigc_inp)
{
  LibData newlib;
  newlib.ConstantCrossSectionData(mat_inp, g, sigc_inp);
  AddLibData(mat_inp, newlib);
};

