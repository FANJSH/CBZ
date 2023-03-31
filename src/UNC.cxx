#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "UNC.h"
#include "UNC_Sensitivity.h"

using namespace std;

UncertaintyCalculation::UncertaintyCalculation
(int g, real *bnd,LibraryCovariance *cinp,
 SensitivityContainer *sinp,LibraryContainer *linp,
 ParametersContainer *cvalinp, Parameters *evalinp,
 ParameterCovariance *ccovinp, ParameterCovariance *ecovinp)
{
  cov_con=cinp;
  sens_con=sinp;
  lib_con=linp;
  cval_con=cvalinp;
  cval_cov=ccovinp;
  eval=evalinp;
  eval_cov=ecovinp;

  PutEnergyGroupInfo(g,bnd);

  uq_on=false;

  count=0;  
};

UncertaintyCalculation::UncertaintyCalculation(int g, real *bnd)
{
  PutEnergyGroupInfo(g,bnd);

  uq_on=false;
  count=0;    
};

UncertaintyCalculation::UncertaintyCalculation(XSLibrary &xslib)
{
  int g=xslib.GetGroup();
  real *bnd=new real[g+1];
  for(int i=0;i<g+1;i++){
    bnd[i]=xslib.GetEnband().get_dat(i);
  };
  PutEnergyGroupInfo(g,bnd);
  delete [] bnd;

  uq_on=false;
  count=0;    
};

UncertaintyCalculation::UncertaintyCalculation()
{
  uq_on=false;
  count=0;    
};


// 

void UncertaintyCalculation::PutEnergyGroupInfo(int g, real *bnd)
{
  grp=g;
  ebnd.put_imax(grp+1);
  for(int i=0;i<grp+1;i++){
    ebnd.put_data(i,bnd[i]);
  };
};

void UncertaintyCalculation::GetCEValue(ParametersContainer &cval_con,Parameters &eval)
{
  int clib=cval_con.GetSize();
  if(clib==0){
    cout<<"There is no parameters.\n";
    exit(0);
  };

  int size=eval.GetParameterList()->GetSize();
  cout<<"*************************************************\n";
  cout<<"           ";
  for(int i=0;i<clib;i++){
    cout<<cval_con.GetName(i)<<" ";
  };
  cout<<"\n";
  cout<<"*************************************************\n";

  for(int i=0;i<size;i++){
    string core=eval.GetParameterList()->GetCoreTag(i);
    string chara=eval.GetParameterList()->GetCharaTag(i);
    int step=eval.GetParameterList()->GetStepTag(i);
    real e=eval.GetValue(i);
    cout<<core<<" "<<chara<<" "<<step<<" : ";
    for(int j=0;j<clib;j++){
      int pos=cval_con.GetParameters(j).GetParameterList()->FindData(core,chara,step);
      real c;
      if(pos!=-1){
        c=cval_con.GetParameters(j).GetValue(pos);
      }else{
	c=0.;
      };
      cout<<c/e<<"  ";
    };
    cout<<"\n";
  };
};

void UncertaintyCalculation::GetCEValue(ParametersContainer &cval_con,Parameters &eval,string libname)
{
  int tmp=cval_con.FindData(libname);
  if(tmp==-1){
    cout<<"There is no parameters "<<libname<<"\n";
    exit(0);
  };

  Parameters cval=cval_con.GetParameters(libname);

  int size=eval.GetParameterList()->GetSize();
  for(int i=0;i<size;i++){
    string core=eval.GetParameterList()->GetCoreTag(i);
    string chara=eval.GetParameterList()->GetCharaTag(i);
    int step=eval.GetParameterList()->GetStepTag(i);
    int pos=cval.GetParameterList()->FindData(core,chara,step);
    if(pos!=-1){
      real c=cval.GetValue(pos);
      real e=eval.GetValue(i);
      cout<<core<<" "<<chara<<" "<<step<<" : "<<c/e<<"\n";
    };
  };
};

real UncertaintyCalculation::CalCrossSectionUncertainty
(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,
 string core,string chara,int step,real mean,bool print)
{
  if(print){
  cout<<"#********************************************************\n";
  cout<<"#* Cross section-induced uncertainty\n";
  cout<<"#********************************************************\n";
  cout<<"#* Integral Data : core   "<<core<<"\n";
  cout<<"#*               : Chara. "<<chara<<"\n";
  cout<<"#*               : step   "<<step<<"\n";
  cout<<"#********************************************************\n";
  };

  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"# There is no sensitivity data ...\n";
    exit(0);
  };

  SensitivityData sens=sens_con.GetSensitivityData(core,chara,step);

  //
  vector<int> d_mat1;
  vector<int> d_mt1;
  vector<int> d_mat2;
  vector<int> d_mt2;
  vector<real> d_var;
  //

  real totsol=0;
  real maxcomp=0.;
  //int max1,max2;
  for(int i=0;i<lib_cov.GetSize();i++){
    CrossSectionCovariance xscov=lib_cov.GetCrossSectionCovariance(i);
    int mat1 =xscov.GetMat1();
    int mt1  =xscov.GetMt1();
    int mat2 =xscov.GetMat2();
    int mt2  =xscov.GetMt2();
    int tmp1=sens.FindData1D(mat1,mt1);
    int tmp2=sens.FindData1D(mat2,mt2);
    if(tmp1!=-1&&tmp2!=-1){
      GroupData1D s1=sens.GetSensitivity1D(mat1,mt1);
      GroupData1D s2=sens.GetSensitivity1D(mat2,mt2);
      GroupData2D cov=xscov.GetCovariance("relative");
      real sol=s1*(cov*s2);

      d_mat1.push_back(mat1);
      d_mt1.push_back(mt1);
      d_mat2.push_back(mat2);
      d_mt2.push_back(mt2);
      d_var.push_back(sol);

      if(mat1!=mat2||mt1!=mt2)sol*=2;
      totsol+=sol;
      if(print&&sqrt(fabs(sol))>mean){
	cout<<"#* ";
        WriteOut(mat1,7);
	cout<<":";
        WriteOut(mt1,3);
	cout<<" & ";
	WriteOut(mat2,7);
	cout<<":";
	WriteOut(mt2,3);
	cout<<"  ";
	if(sol>0.)cout<<" ";
        cout.setf(ios::scientific);
        cout.precision(4);
	cout<<sol;
        cout.unsetf(ios::scientific);
        cout.setf(ios::showpoint);
        cout.precision(4);
        cout<<"  (";
        if(sol>0.){cout<<" ";}else{cout<<"-";};
        cout<<sqrt(fabs(sol))<<")\n";
        cout.unsetf(ios::showpoint);
      };
      //
      if(uq_on){
        if(mt1!=181){
          for(int g=0;g<grp;g++){
            real s1d=s1.get_dat(g);
	    for(int g2=0;g2<grp;g2++){
	      real tmpv=s1d*s2.get_dat(g2)*cov.get_dat(g,g2);
	      if(tmpv>maxcomp){
                maxcomp=tmpv;
	    //max1=mat1;
	    //max2=mt1;
	      };
	    };
	  };
	};
      };
      //
    };
  };

  if(print){
    cout<<"#********************************************************\n";
    cout.setf(ios::scientific);
    cout.precision(4);
    cout<<"#* Total uncertainty     "<<totsol;
    cout.unsetf(ios::scientific);
    cout.setf(ios::showpoint);
    cout.precision(4);
    cout<<"  ("<<sqrt(fabs(totsol))<<")\n";
    cout.unsetf(ios::showpoint);
    cout<<"#********************************************************\n";
  };

  if(uq_on){
  maxcomp=1./maxcomp;
  for(int i=0;i<lib_cov.GetSize();i++){
    CrossSectionCovariance xscov=lib_cov.GetCrossSectionCovariance(i);
    int mat1 =xscov.GetMat1();
    int mt1  =xscov.GetMt1();
    int mat2 =xscov.GetMat2();
    int mt2  =xscov.GetMt2();
    int tmp1=sens.FindData1D(mat1,mt1);
    int tmp2=sens.FindData1D(mat2,mt2);
    if(tmp1!=-1&&tmp2!=-1){
      GroupData1D s1=sens.GetSensitivity1D(mat1,mt1);
      GroupData1D s2=sens.GetSensitivity1D(mat2,mt2);
      GroupData2D cov=xscov.GetCovariance("relative");
      string outfile=filename_uq+"_"+IntToString(mat1)+"_"+IntToString(mt1);
      if(mat2!=mat1||mt2!=mt1)outfile+="_"+IntToString(mat2)+"_"+IntToString(mt2);
      ShowCrossSectionUncertaintyComponentXYPlot(s1,s2,cov,maxcomp,outfile);
    };
  };
  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++
  {
  int sz=d_mat1.size(); 
  vector<int> matid;
  vector<int> mtid;
  for(int i=0;i<sz;i++){
    int tmp=d_mat1[i];
    bool exist=false;
    for(int j=0;j<matid.size();j++){
      if(tmp==matid[j])exist=true;
    };
    if(!exist)matid.push_back(tmp);
    tmp=d_mt1[i];
    exist=false;
    for(int j=0;j<mtid.size();j++){
      if(tmp==mtid[j])exist=true;
    };
    if(!exist)mtid.push_back(tmp);
  };

  int matsz=matid.size();
  int mtsz=mtid.size();

  vector< vector<real> > var(matsz);
  for(int i=0;i<matsz;i++){
    var[i].resize(mtsz+1,0.);
  };

  for(int i=0;i<sz;i++){
    if(d_mat1[i]==d_mat2[i]){
      int jj=-1;
      for(int j=0;j<matsz;j++){
       if(matid[j]==d_mat1[i])jj=j;
      };
      int kk=-1;
      for(int k=0;k<mtsz;k++){
        if(mtid[k]==d_mt1[i])kk=k;
      };
      var[jj][mtsz]+=d_var[i];
      if(d_mt1[i]==d_mt2[i]){
	var[jj][kk]+=d_var[i];
      }else{
	int j2=-1;
        for(int j=0;j<matsz;j++){
         if(matid[j]==d_mat2[i])j2=j;
        };
        int k2=-1;
        for(int k=0;k<mtsz;k++){
          if(mtid[k]==d_mt2[i])k2=k;
        };
	//var[jj][kk]+=d_var[i]*0.5;
	//var[j2][k2]+=d_var[i]*0.5;
      };
    };
  };

  /*
  // For plotting of cumulative data
  cout.setf(ios::scientific);
  cout.precision(3);
  cout<<"#       ";
  for(int j=0;j<mtid.size();j++){
    WriteOut(mtid[j],4);
    cout<<"      ";
  };
  cout<<"\n";
  for(int i=0;i<matid.size();i++){
    //WriteOut(matid[i],8);
    WriteOut(midt.Name(matid[i]),8);
    real sum=0.;
    for(int j=0;j<mtid.size();j++){
      cout<<" "<<sum+sqrt(var[i][j]);
      sum+=sqrt(var[i][j]);
    };
    cout<<" "<<sqrt(var[i][mtid.size()])<<"\n";
  };
  */

  };
  // +++++++++++++++++++++++++++++++++++++++++++++++++++


  return totsol;
};

void UncertaintyCalculation::CalCrossSectionUncertaintyDetail
(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,
 int mat,int mt,string core,string chara,int step)
{
  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"# There is no sensitivity data ...\n";
    cout<<"# "<<core<<" "<<chara<<" "<<step<<"\n";
    exit(0);
  };
  SensitivityData sens=sens_con.GetSensitivityData(core,chara,step);

  int tmp1=sens.FindData1D(mat,mt);
  if(tmp1==-1){
    cout<<"# There is no sensitivity data ...\n";
    cout<<"# (MAT/MT)=("<<mat<<"/"<<mt<<")\n";
    exit(0);
  };
  GroupData1D sen1d=sens.GetSensitivity1D(mat,mt);

  int tmp2=lib_cov.FindData(mat,mt,mat,mt);
  if(tmp2==-1){
    cout<<"# There is no covariance data ...\n";
    cout<<"# (MAT/MT)=("<<mat<<"/"<<mt<<")\n";
    exit(0);
  };
  CrossSectionCovariance xscov=lib_cov.GetCrossSectionCovariance(tmp2);

  cout<<"**************************************************\n";
  cout<<"* Cross section-induced uncertainty\n";
  cout<<"**************************************************\n";
  cout<<"* Integral Data : core   "<<core<<"\n";
  cout<<"*               : Chara. "<<chara<<"\n";
  cout<<"*               : step   "<<step<<"\n";
  cout<<"* Nuclear Data  : (MAT/MT)=("<<mat<<"/"<<mt<<")\n";
  cout<<"**************************************************\n";

  real e_fast=0.;
  real e_inter=0.;
  real e_thml=0.;
  int g_fast=0;
  int g_inter=0;
  int g_thml=0;
  for(int i=0;i<grp;i++){
    real e=ebnd.get_dat(i+1);
    if(e_fast==0.&&e<1e5){
      e_fast=e;
      g_fast=i+1;
    };
    if(e_inter==0.&&e<1e2){
      e_inter=e;
      g_inter=i+1;
    };
    if(e_thml==0.&&e<0.625){
      e_thml=e;
      g_thml=i+1;
    };
  };

  vector<int> g_st(4);
  vector<int> g_ed(4);
  g_st[0]=0;
  g_ed[0]=g_fast;
  g_st[1]=g_fast+1;
  g_ed[1]=g_inter;
  g_st[2]=g_inter+1;
  g_ed[2]=g_thml;
  g_st[3]=g_thml+1;
  g_ed[3]=grp-1;

  GroupData2D cov(4,4);

  real sum=0.;
  for(int i=0;i<4;i++){
    GroupData1D tmp1(grp);
    tmp1.set_zero();
    for(int g=g_st[i];g<=g_ed[i];g++){
      tmp1.put_data(g,sen1d.get_dat(g));
    };
    for(int j=0;j<4;j++){
      GroupData1D tmp2(grp);
      tmp2.set_zero();
      for(int g=g_st[j];g<=g_ed[j];g++){
        tmp2.put_data(g,sen1d.get_dat(g));
      };
      real unc=tmp1*(xscov.GetCovariance("Relative")*tmp2);
      cov.put_data(i,j,unc);
      sum+=unc;
    };
  };

  cout<<"# [Enregy range definition]\n#\n";
  cout<<"#  Fast range     : "<<e_fast<<" [eV]-\n";
  cout<<"#  Inter.-1 range : "<<e_inter<<" [eV]-\n";
  cout<<"#  Inter.-2 range : "<<e_thml<<" [eV]-\n";
  cout<<"#  Thermel range  : below "<<e_thml<<" [eV]\n";
  cout<<"#\n";

  string name[]={
    "Fast   ","Inter-1","Inter-2","Thermal"
  };
  cout<<"# Range-wise variance & Std.Dev.\n#\n";
  for(int i=0;i<4;i++){
    real val=cov.get_dat(i,i);
    cout<<"# "<<name[i]<<" : "<<val<<"   "<<sqrt(val)<<"\n";
  };
  cout<<"#\n";

  //cout<<"# Range-range covariance & Std.Dev.\n#\n";
  cout<<"# Range-range covariance & correlation \n#\n";
  for(int i=0;i<4;i++){
    for(int j=i+1;j<4;j++){
      real val=cov.get_dat(i,j);
      if(val!=0.){
        //val+=cov.get_dat(j,i);
        cout<<"# "<<name[i]<<"-"<<name[j]<<" : ";
        cout<<val<<"   ";
	cout<<val/sqrt(cov.get_dat(i,i)*cov.get_dat(j,j));
	/*
        if(val>=0.){
	  cout<<"+"<<sqrt(val);
        }else{
	  cout<<"-"<<sqrt(-val);
        };
	*/
	cout<<"\n";
      };
    };
  };

  cout<<"**************************************************\n";
}

real UncertaintyCalculation::CalG1MG2
(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,
 string core,string chara,int step,string core2,string chara2,int step2)
{
  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  int tmp2=sens_con.GetParameterList()->FindData(core2,chara2,step2);
  if(tmp==-1||tmp2==-1){
    cout<<"There is no sensitivity data ...\n";
    exit(0);
  };

  SensitivityData sens1=sens_con.GetSensitivityData(core,chara,step);
  SensitivityData sens2=sens_con.GetSensitivityData(core2,chara2,step2);

  return CalG1MG2(sens1,sens2,lib_cov);
}

real UncertaintyCalculation::CalG1MG2
(SensitivityData &sens1, SensitivityData &sens2, LibraryCovariance &lib_cov)
{
  real totsol=0;
  for(int i=0;i<lib_cov.GetSize();i++){
    CrossSectionCovariance xscov=lib_cov.GetCrossSectionCovariance(i);
    int mat1 =xscov.GetMat1();
    int mt1  =xscov.GetMt1();
    int mat2 =xscov.GetMat2();
    int mt2  =xscov.GetMt2();
    int tmp1=sens1.FindData1D(mat1,mt1);
    int tmp2=sens2.FindData1D(mat2,mt2);
    if(tmp1!=-1&&tmp2!=-1){
      GroupData1D dum(grp);
      dum=xscov.GetCovariance("relative")*sens2.GetSensitivity1D(mat2,mt2);
      real sol=sens1.GetSensitivity1D(mat1,mt1)*dum;
      if(mat1!=mat2||mt1!=mt2){
        int tmp3=sens1.FindData1D(mat2,mt2);
        int tmp4=sens2.FindData1D(mat1,mt1);
        if(tmp3!=-1&&tmp4!=-1){
          dum=xscov.GetCovariance("relative")*sens1.GetSensitivity1D(mat2,mt2);
          sol+=sens2.GetSensitivity1D(mat1,mt1)*dum;
	};
      };
      totsol+=sol;
    };
  };
  return totsol;
}


real UncertaintyCalculation::CalUncertaintyWithoutAnyInformation
(string core,string chara,int step)
{
  //real xs_cov=CalCrossSectionUncertainty(core,chara,step,false);
  //real method_cov=cval_cov->GetStandardDeviation(core,chara,step,"Relative");
  //method_cov*=method_cov;
  //cout<<"Cross section-induced component : "<<xs_cov<<" ("<<sqrt(xs_cov)<<")\n";
  //cout<<"Method component                : "<<method_cov<<" ("<<sqrt(method_cov)<<")\n";
  //real tot_cov=xs_cov+method_cov;
  //cout<<"Total uncertainty               : "<<tot_cov<<" ("<<sqrt(tot_cov)<<")\n";
  //return tot_cov;
  return 0.;
};

real UncertaintyCalculation::CalCrossSectionUncertaintyBiasMethod
(string core1,string chara1,string core2,string chara2,int step1,int step2,bool print)
{
  if(print){
  cout<<"****************************************\n";
  cout<<"* Cross section-induced uncertainty with Bias method\n";
  cout<<"****************************************\n";
  cout<<"* Mock-up Core   : "<<core1<<"\n";
  cout<<"*         Chara. : "<<chara1<<"\n";
  cout<<"*         Step   : "<<step1<<"\n";
  cout<<"* Target  Core   : "<<core2<<"\n";
  cout<<"*         Chara. : "<<chara2<<"\n";
  cout<<"*         Step   : "<<step2<<"\n";
  cout<<"****************************************\n";
  };

  int tmp1=sens_con->GetParameterList()->FindData(core1,chara1,step1);
  int tmp2=sens_con->GetParameterList()->FindData(core2,chara2,step2);
  if(tmp1==-1||tmp2==-1){
    cout<<"There is no sensitivity data ...\n";
    exit(0);
  };

  int numcov=cov_con->GetSize();

  real totsol=0;
  for(int i=0;i<numcov;i++){
    int mat1 =cov_con->GetCrossSectionCovariance(i).GetMat1();
    int mt1  =cov_con->GetCrossSectionCovariance(i).GetMt1();
    int mat2 =cov_con->GetCrossSectionCovariance(i).GetMat2();
    int mt2  =cov_con->GetCrossSectionCovariance(i).GetMt2();
    GroupData1D sens11(grp);
    GroupData1D sens12(grp);
    GroupData1D sens21(grp);
    GroupData1D sens22(grp);
    if(sens_con->GetSensitivityData(tmp1).FindData1D(mat1,mt1)!=-1){
      sens11.copy(sens_con->GetSensitivityData(tmp1).GetSensitivity1D(mat1,mt1));
    };
    if(sens_con->GetSensitivityData(tmp2).FindData1D(mat1,mt1)!=-1){
      sens21.copy(sens_con->GetSensitivityData(tmp2).GetSensitivity1D(mat1,mt1));
    };
    if(sens_con->GetSensitivityData(tmp1).FindData1D(mat2,mt2)!=-1){
      sens12.copy(sens_con->GetSensitivityData(tmp1).GetSensitivity1D(mat2,mt2));
    };
    if(sens_con->GetSensitivityData(tmp2).FindData1D(mat2,mt2)!=-1){
      sens22.copy(sens_con->GetSensitivityData(tmp2).GetSensitivity1D(mat2,mt2));
    };
    GroupData1D dum(grp);
    GroupData1D difsens1(grp);
    GroupData1D difsens2(grp);
    difsens1=sens11-sens21;
    difsens2=sens12-sens22;
    dum=cov_con->GetCrossSectionCovariance(i).GetCovariance("Relative")*difsens2;
    real sol=difsens1*dum;
    if(mat1!=mat2||mt1!=mt2)sol*=2;
    totsol+=sol;
    if(sol!=0.){
      if(print){
      cout<<"* "<<mat1<<":"<<mt1<<" & "<<mat2<<":"<<mt2<<"  "<<sol;
      cout.setf(ios::showpoint);
      cout.precision(6);
      cout<<"("<<sqrt(fabs(sol))<<")\n";
      cout.unsetf(ios::showpoint);
      };
    };
  };
  if(print){
  cout<<"****************************************\n";
  cout<<"* Total Unc.: "<<totsol;
  cout.setf(ios::showpoint);
  cout.precision(6);
  cout<<" ("<<sqrt(fabs(totsol))<<")\n";
  cout.unsetf(ios::showpoint);
  cout<<"****************************************\n";
  };

  return totsol;
}

real UncertaintyCalculation::CalUncertaintyBiasMethod
(string core1,string chara1,string core2,string chara2,int step1,int step2,bool print)
{
  real xs_cov=CalCrossSectionUncertaintyBiasMethod(core1,chara1,core2,chara2,step1,step2,false);
  real exp_cov=eval_cov->GetStandardDeviation(core1,chara1,step1,"Relative");
  real mtd1_cov=cval_cov->GetStandardDeviation(core1,chara1,step1,"Relative");
  real mtd2_cov=cval_cov->GetStandardDeviation(core2,chara2,step2,"Relative");
  real corr_mtd=cval_cov->GetCorrelation(core1,chara1,core2,chara2,step1,step2);
  real corr=-1.*(corr_mtd*mtd1_cov*mtd2_cov);
  exp_cov*=exp_cov;
  mtd1_cov*=mtd1_cov;
  mtd2_cov*=mtd2_cov;
  real tot_cov=xs_cov+exp_cov+mtd1_cov+mtd2_cov+corr;
  if(print){
  cout<<"Cross section-induced component : "<<xs_cov<<" ("<<sqrt(xs_cov)<<")\n";
  cout<<"Experimental component          : "<<exp_cov<<" ("<<sqrt(exp_cov)<<")\n";
  cout<<"Mockup method component         : "<<mtd1_cov<<" ("<<sqrt(mtd1_cov)<<")\n";
  cout<<"Target method component         : "<<mtd2_cov<<" ("<<sqrt(mtd2_cov)<<")\n";
  cout<<"Target correlation component         : "<<corr<<" ("<<sqrt(mtd2_cov)<<")\n";
  cout<<"Total uncertainty               : "<<tot_cov<<" ("<<sqrt(tot_cov)<<")\n";
  };
  return tot_cov;
};

void UncertaintyCalculation::SearchBestMockup(string core,string chara,int step)
{
  int size=eval->GetParameterList()->GetSize();
  real errmin=1000000.;
  int tag=-1;
  for(int i=0;i<size;i++){
    string core2=eval->GetParameterList()->GetCoreTag(i);
    string chara2=eval->GetParameterList()->GetCharaTag(i);
    int step2=eval->GetParameterList()->GetStepTag(i);
    if(core!=core2||chara!=chara2||step!=step2){
      real err=CalUncertaintyBiasMethod(core2,chara2,core,chara,step2,step,false);
      if(err<errmin){
        errmin=err;
	tag=i;
      };
    };
  };

  cout<<"The best mockup system for "<<core<<" "<<chara<<" "<<step<<" is...\n";
  string core2=eval->GetParameterList()->GetCoreTag(tag);
  string chara2=eval->GetParameterList()->GetCharaTag(tag);
  int step2=eval->GetParameterList()->GetStepTag(tag);
  cout<<"   "<<core2<<" "<<chara2<<" "<<step2<<"\n";
  cout<<"  Uncertainty is "<<errmin<<" ("<<sqrt(errmin)<<")\n";
};

real UncertaintyCalculation::CalLibraryEffect
(SensitivityContainer &sens_con,LibraryContainer &lib_con,
 string core,string chara,int step,string libname1,string libname2,real mean,bool en,bool print)
{
  if(print){
  cout<<"********************************************************************************\n";
  cout<<"* Library effect calculation            \n";
  cout<<"*       by UncertaintyCalculation       \n";
  cout<<"********************************************************************************\n";
  cout<<"* Integral Data : core   "<<core<<"\n";
  cout<<"*               : Chara. "<<chara<<"\n";
  cout<<"*               : step   "<<step<<"\n";
  cout<<"*   "<<libname1<<" -> "<<libname2<<"\n";
  cout<<"********************************************************************************\n";
  };


  
  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"There is no sensitivity data ...\n";
    exit(0);
  };

  Library lib1=lib_con.GetLibrary(libname1);
  Library lib2=lib_con.GetLibrary(libname2);

  real totval=0.;

  // for 1D xs
  int nxs1=lib1.GetSize();
  for(int i=0;i<nxs1;i++){
    int mat=lib1.GetMatList(i);
    int mt=lib1.GetMtList(i);
    int tmp=lib2.FindData(mat,mt);
    int tmpsens=sens_con.GetSensitivityData(core,chara,step).FindData1D(mat,mt);
    bool exist2d=false;
    if(mt==2||mt==4||mt==16)exist2d=true;
    if(tmp!=-1&&tmpsens!=-1&&!exist2d){
      GroupData1D diff(grp);
      // (fission spectrum normalization)
      real sum1=0.;
      real sum2=0.;
      if(mt==181){
	sum1=lib1.GetCrossSection(i).get_sum();
	sum2=lib2.GetCrossSection(tmp).get_sum();
	if(sum1!=0.)sum1=1./sum1;
	if(sum2!=0.)sum2=1./sum2;
      };
      for(int j=0;j<grp;j++){
	real tmp1=lib1.GetCrossSection(i).get_dat(j);
	real tmp2=lib2.GetCrossSection(tmp).get_dat(j);
	if(mt==181){
	  tmp1*=sum1;
	  tmp2*=sum2;
	};
	if(mt==251||mt==181){
	  diff.put_data(j,tmp2-tmp1);
	}else{
	  if(tmp1!=0.&&tmp2!=0.){
            diff.put_data(j,(tmp2-tmp1)/tmp1);
	  }
          else{
            diff.put_data(j,0.);
          };
	};
      };
      real pos_effect=0.;
      real neg_effect=0.;
      GroupData1D effect(grp);
      effect=diff.mult(sens_con.GetSensitivity1D(mat,mt,core,chara,step));
      for(int j=0;j<grp;j++){
	real tmp=effect.get_dat(j);
	if(tmp>0.){pos_effect+=tmp;}
	else{neg_effect+=tmp;};
      };
      real value=pos_effect+neg_effect;
      if(print){
      if(en&&(pos_effect>mean||neg_effect<mean*-1)){
	if(count!=0)cout<<"\n\n";	
        cout<<"# Energy group-contribution : "<<mat<<" / "<<mt<<"\n";
	cout<<"# Plot as [plot 'output' i "<<count<<":"<<count<<" w st] with Gnuplot.\n";
	cout<<"#              [/lethargy*0.25]\n#\n";
	cout.setf(ios::scientific);
	cout.precision(5);
        for(int j=0;j<grp;j++){
	  real e1=ebnd.get_dat(j);
	  real e2=ebnd.get_dat(j+1);
          real letwid=log(e1/e2);
	  letwid/=0.25;
          cout<<"    "<<e1<<" "<<effect.get_dat(j)/letwid<<"\n";
	};
	count++;
      };
      };
      totval+=value;
      if(print){
      if(fabs(value)>mean||fabs(pos_effect)>mean||fabs(neg_effect)>mean){
	cout.setf(ios::showpoint);
	cout.precision(5);
        cout<<"#* mat:";
        WriteOut(mat,6);
        cout<<" & mt:";
        WriteOut(mt,3);
        cout<<" : ";
        if(value>=0.)cout<<" ";
        cout<<value;
	cout<<"  ("<<pos_effect<<","<<neg_effect<<")\n";
      };
      };
    };
  };

  // for 2D xs
  nxs1=lib1.GetSize2D();
  for(int i=0;i<nxs1;i++){
    int mat=lib1.GetMatList2D(i);
    int mt=lib1.GetMtList2D(i);
    int tmp=lib2.FindData2D(mat,mt);
    int tmpsens=sens_con.GetSensitivityData(core,chara,step).FindData2D(mat,mt);
    if(tmp!=-1&&tmpsens!=-1){
      GroupData2D diff(grp,grp);
      for(int j=0;j<grp;j++){
	//for(int k=j;k<grp;k++){
	  for(int k=0;k<grp;k++){
	  real tmp1=0.;
	  real tmp2=0.;
	  if(lib1.GetCrossSection2D(i).get_y()!=0){
  	    tmp1=lib1.GetCrossSection2D(i).get_dat(j,k);
	  };
	  if(lib2.GetCrossSection2D(tmp).get_y()!=0){
  	    tmp2=lib2.GetCrossSection2D(tmp).get_dat(j,k);
	  };
  	  //real tmp1=lib1.GetCrossSection2D(i).get_dat(j,k);
	  //real tmp2=lib2.GetCrossSection2D(tmp).get_dat(j,k);
          if(mt!=2){
            diff.put_data(j,k,tmp2-tmp1);
	  }else{
	  if(tmp1!=0.&&tmp2!=0.){
            diff.put_data(j,k,(tmp2-tmp1)/tmp1);}
          else{
            diff.put_data(j,k,0.);
          };
	  };
	};
      };
      real pos_effect=0.;
      real neg_effect=0.;
      vector<real> effect(grp);
      GroupData2D senstmp=sens_con.GetSensitivityData(core,chara,step).GetSensitivity2D(mat,mt);
      for(int j=0;j<grp;j++){
        real tmp=0.;
	//for(int k=j;k<grp;k++){
	for(int k=0;k<grp;k++){
	  tmp+=diff.get_dat(j,k)*senstmp.get_dat(j,k);
	};
	if(tmp>0.){pos_effect+=tmp;}
	else{neg_effect+=tmp;};
	effect[j]=tmp;
      };
      real value=pos_effect+neg_effect;
      if(print){
      if(en&&(pos_effect>mean||neg_effect<mean*-1)){
	if(count!=0)cout<<"\n\n";
        cout<<"# Energy group-contribution : "<<mat<<" / "<<mt<<"\n";
	cout<<"# Plot as [plot 'output' i "<<count<<":"<<count<<" w st] with Gnuplot.\n";
	cout<<"#              [/lethargy*0.25]\n#\n";
	cout.setf(ios::scientific);
	cout.precision(5);
        for(int j=0;j<grp;j++){
	  real e1=ebnd.get_dat(j);
	  real e2=ebnd.get_dat(j+1);
          real letwid=log(e1/e2);
	  letwid/=0.25;
          if(print)cout<<"    "<<e1<<"   "<<effect[j]/letwid<<"\n";
	};
	count++;
      };
      };
      totval+=value;
      if(print){
      if(fabs(value)>mean||fabs(pos_effect)>mean||fabs(neg_effect)>mean){
	cout.setf(ios::showpoint);
	cout.precision(5);
        cout<<"#* mat:";
        WriteOut(mat,6);
        cout<<" & mt:";
        WriteOut(mt,3);
        cout<<" : ";
        if(value>=0.)cout<<" ";
	cout<<value;
	cout<<"  ("<<pos_effect<<","<<neg_effect<<")\n";
      };
      };
    };
  };

  if(print){
  cout<<"********************************************************************************\n";
  cout<<"* Total :"<<totval<<"\n";
  /*
  int tmp1=cval_con.FindData(libname1);
  int tmp2=cval_con.FindData(libname2);
  if(tmp1!=-1&&tmp2!=-1){
    int tmp1=cval_con.GetParameters(libname1).GetParameterList()->FindData(core,chara,step);
    int tmp2=cval_con.GetParameters(libname2).GetParameterList()->FindData(core,chara,step);
    if(tmp1!=-1&&tmp2!=-1){
      real c1=cval_con.GetParameters(libname1).GetValue(tmp1);
      real c2=cval_con.GetParameters(libname2).GetValue(tmp2);
      cout<<"* Direct Calc. : "<<(c2-c1)/c1<<"\n";
    };
  };
  */
  cout<<"********************************************************************************\n";
  };

  return totval;
}

void UncertaintyCalculation::CalLibraryEffect
(SensitivityContainer &sens_con,BurnupChain &chain1, BurnupChain &chain2, string core,string chara,int step,real mean)
{
  SensitivityData snsd=sens_con.GetSensitivityData(core,chara,step);
  int sz0=snsd.GetSize0D();
  real sum=0.;
  for(int i=0;i<sz0;i++){
    int mt=snsd.GetMtList0D(i);
    int mat=snsd.GetMatList0D(i);
    if(mt==8888){ // half-life sensitivity
      real dc1=chain1.GetDecayConstant(mat);
      real dc2=chain2.GetDecayConstant(mat);
      if(dc1!=0.&&dc2!=0.){
	real hl1=1./dc1;
	real hl2=1./dc2;
	real rdif=(hl2-hl1)/hl1;
	real eff=snsd.GetSensitivity0D(mat,mt)*rdif;
	sum+=eff;
	if(fabs(eff)>mean){
   	  cout<<"# Half-life of "<<midt.Name(mat)<<" : "<<eff<<"\n";
	};
      };
    }else if(mt>18000000){
      int matf=mt-18000000;
      real fy1=chain1.GetFissionYield(matf,mat);
      real fy2=chain2.GetFissionYield(matf,mat);
      if(fy1>0.&&fy2>0.){
	real rdif=(fy2-fy1)/fy1;
	real eff=snsd.GetSensitivity0D(mat,mt)*rdif;
	sum+=eff;
	if(fabs(eff)>mean){
   	  cout<<"# Yield of "<<midt.Name(mat)<<" from "<<midt.Name(matf)<<" : "<<eff<<"\n";
	};
      };
    };
  };
  cout<<"#\n# Total effect : "<<sum<<"\n#\n";
};

void UncertaintyCalculation::EstimateCEValueFromSensitivity
(SensitivityContainer &sens_con,LibraryContainer &lib_con,string libname1,string libname2,
 Parameters &cal_val, Parameters &exp_val, ParameterCovariance &exp_cov, string *name)
{
  int sz=sens_con.GetParameterList()->GetSize();
  GroupData1D exp_std=exp_cov.GetStandardDeviation("Relative");
  GroupData2D exp_covmat=exp_cov.GetCovariance("Relative");
  real chi2_1=0.;
  real chi2_2=0.;
  cout<<"#                       ";
  WriteOut(libname1,9);
  WriteOut(libname2,9);
  cout<<"Exp. Unc.\n";
  GroupData1D ce_pre(sz);
  GroupData1D ce_post(sz);
  for(int i=0;i<sz;i++){
    string core=sens_con.GetParameterList()->GetCoreTag(i);
    string chara=sens_con.GetParameterList()->GetCharaTag(i);
    int step=sens_con.GetParameterList()->GetStepTag(i);
    real dif=CalLibraryEffect(sens_con,lib_con,core,chara,step,libname1,libname2,0.0001,false,false);
    real ce=cal_val.GetValue(i)/exp_val.GetValue(i);
    real ce2=ce*(1.+dif);
    ce_pre.put_data(i,ce-1.);
    ce_post.put_data(i,ce2-1.);
    real err=exp_std.get_dat(i);
    chi2_1+=(ce-1)*(ce-1)/err/err;
    chi2_2+=(ce2-1)*(ce2-1)/err/err;
    WriteOut(i+1,2);
    cout<<" ";
    WriteOut(name[i],20);
    cout<<"  ";
    WriteOut(ce,"%6.4f");
    cout<<"   ";
    WriteOut(ce2,"%6.4f");
    cout<<"   ";
    WriteOut(err,"%6.4f");
    cout<<" 1.\n";
  };
  exp_covmat.DoInverse();
  chi2_1=(exp_covmat*ce_pre)*ce_pre;
  chi2_2=(exp_covmat*ce_post)*ce_post;
  chi2_1/=sz;
  chi2_2/=sz;
  cout<<"#\n";
  cout<<"#  Chi-square/DOF  Base library      : "<<chi2_1<<"\n";
  cout<<"#                  Perturbed library : "<<chi2_2<<"\n";
};


void UncertaintyCalculation::EstimateCEValueFromSensitivity
(SensitivityContainer &sens_con,LibraryContainer &lib_con,string libname1,string libname2,
 Parameters &cal_val, Parameters &exp_val, ParameterCovariance &exp_cov)
{
  int sz=sens_con.GetParameterList()->GetSize();
  string *name=new string[sz];

  for(int i=0;i<sz;i++){
    string core=sens_con.GetParameterList()->GetCoreTag(i);
    string chara=sens_con.GetParameterList()->GetCharaTag(i);
    int step=sens_con.GetParameterList()->GetStepTag(i);
    name[i]=core+"_"+chara+"_"+IntToString(step);
  };

  EstimateCEValueFromSensitivity(sens_con,lib_con,libname1,libname2,cal_val,exp_val,exp_cov,name);

  delete [] name;
};


ParameterCovariance UncertaintyCalculation::GetGMG
(Library &lib,LibraryCovariance &xscov_con,Parameters &cal_val,SensitivityContainer &sens_con)
{
  int cal_size=cal_val.GetParameterList()->GetSize();
  int cov_size=xscov_con.GetSize();

  ParameterCovariance gmg(cal_val);

  GroupData2D gmg_data(cal_size,cal_size);

  // GMG
  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=xscov_con.GetCrossSectionCovariance(i);
    GroupData2D covdata=cov.GetCovariance("Relative");
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    int mat2=cov.GetMat2();
    int mt2 =cov.GetMt2();
    for(int j1=0;j1<cal_size;j1++){
      for(int j2=0;j2<cal_size;j2++){
	int tmp1=sens_con.GetSensitivityData(j1).FindData1D(mat1,mt1);
	int tmp2=sens_con.GetSensitivityData(j2).FindData1D(mat2,mt2);
	if(tmp1!=-1&&tmp2!=-1){
	  gmg_data.add_data
          (j1,j2,sens_con.GetSensitivityData(j1).GetSensitivity1D(mat1,mt1)*
                 (covdata*sens_con.GetSensitivityData(j2).GetSensitivity1D(mat2,mt2)));
	};
	tmp1=sens_con.GetSensitivityData(j1).FindData1D(mat2,mt2);
	tmp2=sens_con.GetSensitivityData(j2).FindData1D(mat1,mt1);
        if(mat1!=mat2||mt1!=mt2){
    	  if(tmp1!=-1&&tmp2!=-1){
	    gmg_data.add_data(j1,j2,sens_con.GetSensitivityData(j1).GetSensitivity1D(mat2,mt2)
	      *(covdata.GetTransposedMatrix()*sens_con.GetSensitivityData(j2).GetSensitivity1D(mat1,mt1)));
	  };
	};
      };
    };
  };

  gmg.PutCov(gmg_data,"Relative");

  return gmg;
};

GroupData2D UncertaintyCalculation::GetGMG
(LibraryCovariance &xscov_con,Parameters &cal_val1,Parameters &cal_val2,SensitivityContainer &sens_con1,SensitivityContainer &sens_con2)
{
  int cal_size1=cal_val1.GetParameterList()->GetSize();
  int cal_size2=cal_val2.GetParameterList()->GetSize();
  int cov_size=xscov_con.GetSize();

  GroupData2D gmg_data(cal_size1,cal_size2);

  // GMG
  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=xscov_con.GetCrossSectionCovariance(i);
    GroupData2D covdata=cov.GetCovariance("Relative");
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    int mat2=cov.GetMat2();
    int mt2 =cov.GetMt2();
    
    for(int j1=0;j1<cal_size1;j1++){
      for(int j2=0;j2<cal_size2;j2++){
	int tmp1=sens_con1.GetSensitivityData(j1).FindData1D(mat1,mt1);
	int tmp2=sens_con2.GetSensitivityData(j2).FindData1D(mat2,mt2);
	if(tmp1!=-1&&tmp2!=-1){
	  gmg_data.add_data
          (j1,j2,sens_con1.GetSensitivityData(j1).GetSensitivity1D(mat1,mt1)*
                 (covdata*sens_con2.GetSensitivityData(j2).GetSensitivity1D(mat2,mt2)));
	};
	tmp1=sens_con1.GetSensitivityData(j1).FindData1D(mat2,mt2);
	tmp2=sens_con2.GetSensitivityData(j2).FindData1D(mat1,mt1);
        if(mat1!=mat2||mt1!=mt2){
    	  if(tmp1!=-1&&tmp2!=-1){
	    gmg_data.add_data(j1,j2,sens_con1.GetSensitivityData(j1).GetSensitivity1D(mat2,mt2)
	      *(covdata.GetTransposedMatrix()*sens_con2.GetSensitivityData(j2).GetSensitivity1D(mat1,mt1)));
	  };
	};
      };
    };
  };

  return gmg_data;
};

void UncertaintyCalculation::DoUncertaintyLibraryAdjustment
(SensitivityContainer &sens_con, ParameterCovariance &cal_cov,
 ParameterCovariance &exp_cov, Parameters &cal_val,
 Parameters &exp_val, Library &lib, LibraryCovariance &lib_cov)
{
  // --------------------------------------------------------------------------------
  //   Fundamental equations
  // --------------------------------------------------------------------------------
  //   T' = T+MG^t (GMG^t + Ve + Vm)^{-1} (Rexp - R(T))
  //
  //   M' = M-MG^t (GMG^t + Ve + Vm)^{-1} GM
  //
  //   R(T') = R(T) + G(T'-T)
  //         = R(T) + GMG^t (GMG^t + Ve + Vm)^{-1} (Rexp - R(T))
  //
  //   Rexp - R(T') = (Rexp - R(T)) - GMG^t (GMG^t + Ve + Vm)^{-1} (Rexp - R(T))
  //
  //   GM'G^t   = GMG^t + GMG^t (GMG^t + Ve + Vm)^{-1} GMG^t
  //   G'M'G'^t = G'MG^t + G'MG^t (GMG^t + Ve + Vm)^{-1} GMG'^t
  // --------------------------------------------------------------------------------

  int sens_size=sens_con.GetParameterList()->GetSize();

  Parameters ce(cal_val.GetParameterList()); // C/E
  for(int i=0;i<sens_size;i++){
    ce.PutValue(i,cal_val.GetValue(i)/exp_val.GetValue(i));
  };

  ParameterCovariance xs_cov=GetGMG(lib,lib_cov,cal_val,sens_con);
  // GMG^t : Covariance matrix for cross section induced component

  GroupData2D inv(sens_size,sens_size);
  inv=xs_cov.GetCovariance("Relative")+cal_cov.GetCovariance("Relative")+exp_cov.GetCovariance("Relative");

  ParameterCovariance total_cov(ce);
  total_cov.PutCov(inv,"Relative");

  inv=inv.inverse();

  GroupData1D ce_m1(sens_size); // (C/E-1.0)
  for(int i=0;i<sens_size;i++){
    ce_m1.put_data(i,ce.GetValue(i)-1.0);
  };

  GroupData1D kai(sens_size);
  kai=xs_cov.GetCovariance()*(inv*ce_m1); // GMG^t (GMG^t + Ve + Vm)^{-1} (Rexp-R(T)) 

  GroupData2D pgmg(sens_size,sens_size); // 
  pgmg=xs_cov.GetCovariance()-xs_cov.GetCovariance()*(inv*xs_cov.GetCovariance());
  // GM'G^t = GMG^t - GMG^t (GMG^t + Ve + Vm)^{-1} GMG^t

  GroupData1D ce_m1_post=ce_m1-kai; // (C/E - 1.0) after adjustment

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  GroupData1D covt=total_cov.GetStandardDeviation("Relative");
  GroupData1D cov1=xs_cov.GetStandardDeviation("Relative");
  GroupData1D cov2=cal_cov.GetStandardDeviation("Relative");
  GroupData1D cov3=exp_cov.GetStandardDeviation("Relative");

  cout<<"****************************************\n";
  cout<<"* Cross section adjustment              \n";
  cout<<"*       by UncertaintyCalculation       \n";
  cout<<"****************************************\n";
  cout<<"* Number of integral data    : "<<sens_size<<"\n";
  cout<<"****************************************\n";  

  cout<<"* Prior C/E values with total uncertainties\n";
  cout<<"*\n";
  cout<<"*                                 C/E             Uncertainty\n";
  cout<<"*                                       Total    ND       Cal.     Exp.\n";
  for(int i=0;i<sens_size;i++){
    cout<<"* ";
    WriteOut(i,4);
    cout<<" ";
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+"_"+
      sens_con.GetParameterList()->GetCharaTag(i)+"_"+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    WriteOut(tmp,25);
    real tmp1=ce.GetValue(i);
    WriteOut(tmp1,"%6.4f");
    WriteOut(covt.get_dat(i),"%9.5f");
    WriteOut(cov1.get_dat(i),"%9.5f");
    WriteOut(cov2.get_dat(i),"%9.5f");
    WriteOut(cov3.get_dat(i),"%9.5f");
    cout<<"\n";
  };
  cout<<"****************************************\n";  
  cout<<"* Correlation matrix of GMG+Ve+Vm\n";
  cout<<"****************************************\n";  
  total_cov.GetCorrelationMatrix().show_self(false);
  cout<<"****************************************\n";  

  real kais=ce_m1*(inv*ce_m1);
  cout<<"*    Kai square = "<<kais<<" ("<<kais/sens_size<<")\n";

  // Only diagonal element of covariance matirx
  GroupData2D tt(sens_size,sens_size);
  for(int i=0;i<sens_size;i++){
    tt.put_data(i,i,inv.get_dat(i,i));
  };
  real kais2=ce_m1*(tt*ce_m1);
  cout<<"*    Kai square (diag) = "<<kais2<<" ("<<kais2/sens_size<<")\n";
  cout<<"****************************************\n";  
  kais2/=sens_size;

  GroupData1D cov1r=pgmg.get_diag();
  pgmg=pgmg+cal_cov.GetCovariance("Relative")+exp_cov.GetCovariance("Relative");

  cout<<"****************************************\n";
  cout<<"* Adjustment\n";
  cout<<"****************************************\n";
  cout<<"* Posterior C/E values with total uncertainties\n";
  cout<<"*\n";
  cout<<"*                                 C/E             Uncertainty\n";
  cout<<"*                                       Total    ND       Cal.     Exp.\n";
  cout<<"*\n";
  for(int i=0;i<sens_size;i++){
    cout<<"* ";
    WriteOut(i,4);
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+"_"+
      sens_con.GetParameterList()->GetCharaTag(i)+"_"+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    cout<<" ";
    WriteOut(tmp,25);
    real tmp1=ce.GetValue(i)-kai.get_dat(i);
    WriteOut(tmp1,"%6.4f");
    WriteOut(sqrt(pgmg.get_dat(i,i)),"%9.5f");
    WriteOut(sqrt(cov1r.get_dat(i)),"%9.5f");
    WriteOut(cov2.get_dat(i),"%9.5f");
    WriteOut(cov3.get_dat(i),"%9.5f");
    cout<<"\n";
  };

  GroupData1D covtr=pgmg.get_diag();
  pgmg=pgmg.inverse();
  real kaisp=(ce_m1-kai)*(pgmg*(ce_m1-kai));
  cout<<"****************************************\n";
  cout<<"*    Kai square = "<<kaisp<<" ("<<kaisp/sens_size<<")\n";
  cout<<"****************************************\n";
  kaisp/=sens_size;

  cout<<"###########################################################\n";
  cout<<"#  Prior and posterior C/E values for data handling\n#\n";
  cout<<"#                                            [(C/E-1) normalized by STDEV]\n#\n";
  for(int i=0;i<sens_size;i++){
    WriteOut(i,4);
    //WriteOut(i%18,4);    
    cout<<" ";
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+"_"+
      sens_con.GetParameterList()->GetCharaTag(i)+"_"+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    WriteOut(tmp,25);
    real cem1_pre=ce_m1.get_dat(i);
    real cem1_post=ce_m1_post.get_dat(i);
    cout<<" ";
    WriteOut(cem1_pre+1.,"%6.4f");
    cout<<" ";
    WriteOut(cem1_post+1.,"%6.4f");
    cout<<"      ";
    if(cem1_pre>0)cout<<" ";
    WriteOut(cem1_pre/covt.get_dat(i),"%6.4f");
    cout<<"  ";
    if(cem1_post>0)cout<<" ";    
    WriteOut(cem1_post/sqrt(covtr.get_dat(i)),"%6.4f");    
    cout<<"\n";
    //if(i%18==17)cout<<"\n\n";
  };
  cout<<"###########################################################\n";
  //
  vector<int> mmat;
  vector<int> mmt;
  int mmnn=0;
  int cov_size=lib_cov.GetSize();
  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=lib_cov.GetCrossSectionCovariance(i);
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    if(mmnn==0){
      mmat.push_back(mat1);
      mmt.push_back(mt1);
      mmnn++;
    }else{
      bool exist=false;
      for(int j=0;j<mmnn;j++){
	if(mat1==mmat[j]&&mt1==mmt[j])exist=true;
      };
      if(!exist){
	mmat.push_back(mat1);
	mmt.push_back(mt1);
	mmnn++;
      };
    };
    int mat2=cov.GetMat2();
    int mt2=cov.GetMt2();
    bool exist=false;
    for(int j=0;j<mmnn;j++){
      if(mat2==mmat[j]&&mt2==mmt[j])exist=true;
    };
    if(!exist){
      mmat.push_back(mat2);
      mmt.push_back(mt2);
      mmnn++;
    };
  };

  vector< vector<GroupData1D> > tmpmat(mmnn);
  for(int i=0;i<mmnn;i++){
    tmpmat[i].resize(sens_size);
    for(int j=0;j<sens_size;j++){
      tmpmat[i][j].put_imax(grp);
      tmpmat[i][j].set_zero();
    };
  };

  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=lib_cov.GetCrossSectionCovariance(i);
    GroupData2D covdata=cov.GetCovariance("Relative");
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    int mat2=cov.GetMat2();
    int mt2 =cov.GetMt2();
    int index1=-1;
    int index2=-1;
    for(int j=0;j<mmnn;j++){
      if(mat1==mmat[j]&&mt1==mmt[j])index1=j;
      if(mat2==mmat[j]&&mt2==mmt[j])index2=j;
    };
    for(int j=0;j<sens_size;j++){
      int tmp1=sens_con.GetSensitivityData(j).FindData1D(mat1,mt1);
      int tmp2=sens_con.GetSensitivityData(j).FindData1D(mat2,mt2);
      if(index1!=-1&&tmp2!=-1){
        tmpmat[index1][j]=tmpmat[index1][j]+
	  covdata*sens_con.GetSensitivityData(j).GetSensitivity1D(tmp2);
      };
      if(index2!=-1&&tmp1!=-1&&index2!=index1){
        tmpmat[index2][j]=tmpmat[index2][j]+
	  covdata.GetTransposedMatrix()*sens_con.GetSensitivityData(j).GetSensitivity1D(tmp1);
      };
    };
  };



  cout<<"###########################################################\n";
  cout<<"# Change in cross section\n";
  cout<<"###########################################################\n";
  GroupData1D ce_n1(sens_size);
  for(int i=0;i<sens_size;i++){
    ce_n1.put_data(i,1.-ce.GetValue().get_dat(i));
  };
  int index=0;
  for(int i=0;i<mmnn;i++){ //
    int mat=mmat[i];
    int mt=mmt[i];
    int fd=lib_cov.FindData(mat,mt,mat,mt);
    int fd2=lib.FindData(mat,mt);
    if(fd!=-1&&fd2!=-1){
      //if(fd!=-1&&fd2!=-1&&mt==18){      
      cout<<"#  MAT: "<<mat<<"  MT: "<<mt<<"\n";
      cout<<"#  Index: "<<index<<"\n";
      index++;
      cout<<"#grp Upper Eng.  XS-Rel.Chg.  XS-Org.     XS-Mod.  Rel.Std.Dev.\n";
      GroupData1D xschange(grp);
      for(int j=0;j<grp;j++){
        real tmp=0.;
        for(int k=0;k<sens_size;k++){
  	  tmp+=tmpmat[i][k].get_dat(j)*(inv*ce_n1).get_dat(k);
        };
        xschange.put_data(j,tmp);
        real xs=lib.GetCrossSection(mat,mt).get_dat(j);
        real stdev=lib_cov.GetCrossSectionCovariance(mat,mt,mat,mt).GetCovariance("Relative").get_diag().sqrt1d().get_dat(j);
        printf("%3d %11.4e %10.5f %11.4e %11.4e %10.5f\n",j,ebnd.get_dat(j),tmp,xs,xs+xs*tmp,stdev);
        //lib.GetCrossSection(mat,mt).put_data(j,xs*(1.+tmp));
      };
#if 1
      cout<<"###########################################################\n";
      cout<<"# ( Impact on responce/integral parameter )\n";
      cout<<"#  MAT: "<<mat<<"  MT: "<<mt<<"\n";
      for(int k=0;k<sens_size;k++){
        int ttt=sens_con.GetSensitivityData(k).FindData1D(mat,mt);
        if(ttt!=-1){
	  real dce=sens_con.GetSensitivityData(k).GetSensitivity1D(ttt)*xschange;
	  if(dce>1e-3){
  	    cout<<"#   + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	    cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	    cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	    cout<<dce<<"\n";
	  };
	};
      };
      cout<<"###########################################################\n";
#endif
      cout<<"\n\n";
    };
  };

  
  int id=0;
  for(int i=0;i<mmnn;i++){
    int mat=mmat[i];
    int mt=mmt[i];
    int fd=lib_cov.FindData(mat,mt,mat,mt);
    int fd2=lib.FindData(mat,mt);
    if(fd!=-1&&fd2!=-1){
      cout<<"#  "<<id<<" "<<mmat[i]<<" "<<mmt[i]<<"\n";
      id++;
    };
  };
  /*
  //MG(GMG+VE+VM)^{-1}
  vector< vector<GroupData1D> > tmpmat2(ic);
  for(int i=0;i<ic;i++){
    tmpmat2[i].resize(is);
    for(int j=0;j<is;j++){
      tmpmat2[i][j].put_imax(grp);
      tmpmat2[i][j].set_zero();
      for(int k=0;k<grp;k++){
	real tmp=0.;
	for(int l=0;l<is;l++){
	  tmp+=tmpmat[i][l].get_dat(k)*inv.get_dat(j,l);
	};
	tmpmat2[i][j].put_data(k,tmp);
      }; 
    };
  };

  cout<<"****************************************\n";
  cout<<"* Reduction of standard deviation\n";
  cout<<"****************************************\n";
  cout<<"*                  (before  -> after)\n";
  for(int ii=0;ii<ic;ii++){
    cout<<"("<<mat[ii]<<":"<<mt[ii]<<")\n";
    for(int i=0;i<grp;i++){
      real tmp=0.;
      for(int j=0;j<is;j++){
        tmp+=tmpmat2[ii][j].get_dat(i)*tmpmat[ii][j].get_dat(i);
      };
      real var=GetXSCovMatrix(mat[ii],mt[ii],mat[ii],mt[ii],"relative").get_diag().get_dat(i);
      real newvar=var-tmp;
      printf("%3d %11.4e %10.5f %10.5f\n",i,ebnd.get_dat(i),sqrt(var),sqrt(newvar));
    };
  };
  */
};

void UncertaintyCalculation::DoUncertaintyLibraryAdjustmentWithNewIntegralDataPrediction
(SensitivityContainer &sens_con, ParameterCovariance &cal_cov,
 ParameterCovariance &exp_cov, Parameters &cal_val,
 Parameters &exp_val, Library &lib, LibraryCovariance &lib_cov,
 SensitivityContainer &sens_con_im, Parameters &cal_val_im)
{
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // The following is taken from the "DoUncertaintyLibraryAdjustment" method

  int sens_size=sens_con.GetParameterList()->GetSize();

  Parameters ce(cal_val.GetParameterList()); // C/E
  for(int i=0;i<sens_size;i++){
    ce.PutValue(i,cal_val.GetValue(i)/exp_val.GetValue(i));
  };

  ParameterCovariance xs_cov=GetGMG(lib,lib_cov,cal_val,sens_con);
  // GMG^t : Covariance matrix for cross section induced component

  GroupData2D inv(sens_size,sens_size);
  inv=xs_cov.GetCovariance("Relative")+cal_cov.GetCovariance("Relative")+exp_cov.GetCovariance("Relative");

  ParameterCovariance total_cov(ce);
  total_cov.PutCov(inv,"Relative");

  inv=inv.inverse();

  GroupData1D ce_m1(sens_size); // (C/E-1.0)
  for(int i=0;i<sens_size;i++){
    ce_m1.put_data(i,ce.GetValue(i)-1.0);
  };

  GroupData1D kai(sens_size);
  kai=xs_cov.GetCovariance()*(inv*ce_m1); // GMG^t (GMG^t + Ve + Vm)^{-1} (Rexp-R(T)) 

  GroupData2D pgmg(sens_size,sens_size); // 
  pgmg=xs_cov.GetCovariance()-xs_cov.GetCovariance()*(inv*xs_cov.GetCovariance());
  // GM'G^t = GMG^t - GMG^t (GMG^t + Ve + Vm)^{-1} GMG^t

  GroupData1D ce_m1_post=ce_m1-kai; // (C/E - 1.0) after adjustment
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int sens_size_im=sens_con_im.GetParameterList()->GetSize();
  /*
  ParameterList plist_im=sens_con_im.GetParameterList();
  Parameters cal_val_im(&plist_im);
  for(int i=0;i<sens_size_im;i++){
    real tmp=sens_con_im.GetSensitivityData(i).GetValue();
    cal_val_im.PutValue(i,tmp);
  };
  */

  // G'MG'
  GroupData2D gmg_im=GetGMG(lib,lib_cov,cal_val_im,sens_con_im).GetCovariance("Relative");
  // GMG'
  GroupData2D gmg_im1=GetGMG(lib_cov,cal_val,cal_val_im,sens_con,sens_con_im);
  // G'MG
  GroupData2D gmg_im2=GetGMG(lib_cov,cal_val_im,cal_val,sens_con_im,sens_con);

  GroupData1D kai_im(sens_size_im);
  kai_im=gmg_im2*(inv*ce_m1);

  GroupData2D pgmg_im(sens_size_im,sens_size_im);
  pgmg_im=gmg_im-gmg_im2*(inv*gmg_im1);

  cout<<"****************************************\n";
  cout<<"* Adjustment results in imaginary case\n";
  cout<<"****************************************\n";
  cout<<"* Number of imaginary case    : "<<sens_size_im<<"\n";
  cout<<"****************************************\n";
  cout<<"* Prior and posterior calculation values\n";
  cout<<"* with nuclear data induced relative uncertainties\n";
  cout<<"*\n";
  for(int i=0;i<sens_size_im;i++){
    cout<<"* "<<sens_con_im.GetParameterList()->GetCoreTag(i);
    cout<<","<<sens_con_im.GetParameterList()->GetCharaTag(i);
    cout<<","<<sens_con_im.GetParameterList()->GetStepTag(i);
    cout.setf(ios::showpoint);
    cout.precision(5);
    cout<<": "<<cal_val_im.GetValue(i);
    cout<<" (";
    cout<<sqrt(gmg_im.get_dat(i,i));
    cout<<") -> ";
    cout<<cal_val_im.GetValue(i)*(1.-kai_im.get_dat(i));
    cout<<" (";
    cout<<sqrt(pgmg_im.get_dat(i,i));
    cout<<")\n";
    cout.unsetf(ios::showpoint);
  };
  cout<<"\n";
  //+++++++++++++++++++++++++++++++++++++++++

  cout<<"# For editing.\n#\n";
  for(int i=0;i<sens_size_im;i++){
    cout<<i<<" ";
    cout<<i+1<<" ";
    cout<<cal_val_im.GetValue(i)<<" ";    
    cout<<sqrt(gmg_im.get_dat(i,i))<<" ";
    cout<<cal_val_im.GetValue(i)*(1.-kai_im.get_dat(i))<<" ";    
    cout<<sqrt(pgmg_im.get_dat(i,i))<<" ";
    cout<<"\n";
  };


};

void UncertaintyCalculation::CalCEwithUncertainty
  (SensitivityContainer &sens_con, ParameterCovariance &cal_cov, ParameterCovariance &exp_cov,
   Parameters &cal_val, Parameters &exp_val,
   IndependentYieldCovariance &icov, HalfLifeCovariance &hcov,
   DecayEnergyCovariance &dcov, BranchingRatioCovariance &bcov, bool adjustment)
{
  int sens_size=sens_con.GetParameterList()->GetSize();

  cout<<"*---------------------------------------------------------------------\n";
  cout<<"* Summary of C/E values with their uncertainties\n*\n";
  cout<<"*    Number of integral data : "<<sens_size<<"\n";
  cout<<"*---------------------------------------------------------------------\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* C and E values with thier ABSOLUTE uncertainties\n";
  cout<<"*\n*   (ND-induced uncertainty is NOT considered.)\n";
  cout<<"*---------------------------------------------------------------------\n";
  cout<<"*  Integral data             Cal.       Unc.       Exp.       Unc.\n";
  cout<<"*---------------------------------------------------------------------\n";
  for(int i=0;i<sens_size;i++){
    cout<<"* ";
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+" / "+
      sens_con.GetParameterList()->GetCharaTag(i)+" / "+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    WriteOut(tmp,25);
    cout<<cal_val.GetValue(i)<<"  ";
    cout<<sqrt(cal_cov.GetCovariance("Absolute").get_dat(i,i))<<"  ";
    cout<<exp_val.GetValue(i)<<"  ";
    cout<<sqrt(exp_cov.GetCovariance("Absolute").get_dat(i,i))<<"  ";
    cout<<"\n";
  };
  cout<<"*---------------------------------------------------------------------\n";

  //ParameterCovariance gmg_xs=GetGMG(lib,lib_cov,cal_val,sens_con);
  UncertaintyCalculationForYieldDecay ucyd;
  GroupData2D gmg_nd=
    //gmg_xs.GetCovariance("Relative")+
     ucyd.GetFissionYieldGMG(sens_con,icov)
    +ucyd.GetHalfLifeGMG(sens_con,hcov)
    +ucyd.GetDecayEnergyGMG(sens_con,dcov)
    +ucyd.GetBranchingRatioGMG(sens_con,bcov);

  GroupData2D inv(sens_size,sens_size);
  inv=gmg_nd+cal_cov.GetCovariance("Relative")
            +exp_cov.GetCovariance("Relative");

  /*
  cout<<"#\n# ND-induced\n#\n";
  gmg_nd.show_self();

  cout<<"#\n# Exp.\n#\n";
  exp_cov.GetCovariance("Relative").show_self();

  cout<<"#\n# All\n#\n";
  inv.show_self();
  */
  
  Parameters ce(cal_val.GetParameterList());
  for(int i=0;i<sens_size;i++){
    ce.PutValue(i,cal_val.GetValue(i)/exp_val.GetValue(i));
  };

  ParameterCovariance total_cov(ce);
  total_cov.PutCov(inv,"Relative");

  ParameterCovariance nd_cov(ce);
  nd_cov.PutCov(gmg_nd,"Relative");

  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* C/E values with RELATIVE uncertainties\n";
  cout<<"*---------------------------------------------------------------------\n";
  cout<<"*  Integral data               C/E       Uncertainty\n";
  cout<<"*                                      All    ND     Exp\n";
  cout<<"*---------------------------------------------------------------------\n";
  for(int i=0;i<sens_size;i++){
    cout<<"* ";
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+" / "+
      sens_con.GetParameterList()->GetCharaTag(i)+" / "+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    WriteOut(tmp,25);
    cout<<" : ";
    WriteOut(ce.GetValue(i),"%6.4f");
    cout<<"  ";
    WriteOut(sqrt(inv.get_dat(i,i)),"%5.3f");
    cout<<"  ";
    WriteOut(sqrt(gmg_nd.get_dat(i,i)),"%5.3f");
    cout<<"  ";
    WriteOut(sqrt(exp_cov.GetCovariance("Relative").get_dat(i,i)),"%5.3f");
    cout<<"\n";
  };
  cout<<"*---------------------------------------------------------------------\n";

  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* Correlation matrix of GMG+Ve+Vm\n";
  cout<<"*---------------------------------------------------------------------\n";
  total_cov.GetCorrelationMatrix().show_self(false);
  total_cov.GetCorrelationMatrix().show_plot("./","cor_matrix");
  nd_cov.GetCorrelationMatrix().show_plot("./","cor_matrix_nd");
  cout<<"*---------------------------------------------------------------------\n";

  GroupData2D inv2=inv;

  inv=inv.inverse();
  GroupData1D diff(sens_size);
  for(int i=0;i<sens_size;i++){
    diff.put_data(i,ce.GetValue(i)-1.0);
  };
  real kais=diff*(inv*diff);

  GroupData2D tt(sens_size,sens_size);
  for(int i=0;i<sens_size;i++){
    tt.put_data(i,i,1./inv2.get_dat(i,i));   // Only diagonal element of covariance matirx
  };
  real kais2=diff*(tt*diff);

  cout.setf(ios::showpoint);
  cout.precision(5);
  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* Chi-square over DOF\n*\n";
  cout<<"*    DOF : "<<sens_size<<"\n";
  cout<<"*                    "<<kais/sens_size<<"\n";
  cout<<"*  (w/o correlation) "<<kais2/sens_size<<"\n";
  cout<<"*---------------------------------------------------------------------\n";

  if(!adjustment)return;

  // +++ NUCLEAR DATA ADJUSTMENT PROCEDURE

  ofstream fout_xschg;
  fout_xschg.open("./output_xschg",ios::out);
  if(fout_xschg.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# [output_xschg]\n";
    exit(0);
  };

  GroupData1D kai(sens_size);
  kai=gmg_nd*(inv*diff);
  GroupData2D pgmg(sens_size,sens_size);
  pgmg=gmg_nd-gmg_nd*(inv*gmg_nd);

  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* Adjustment\n";
  cout<<"*---------------------------------------------------------------------\n";
  cout<<"* Posterior C/E values with nuclear data induced uncertainties\n";
  cout<<"*\n";
  for(int i=0;i<sens_size;i++){
    cout<<"* ";
    string tmp=sens_con.GetParameterList()->GetCoreTag(i)+" / "+
      sens_con.GetParameterList()->GetCharaTag(i)+" / "+
      IntToString(sens_con.GetParameterList()->GetStepTag(i));
    WriteOut(tmp,25);
    cout<<" : ";
    //real ce_post=ce.GetValue(i)-kai.get_dat(i);//Absolute?(Original)
    real ce_post=ce.GetValue(i)*(1.-kai.get_dat(i));//Relative?
    WriteOut(ce_post,"%6.4f");
    cout<<"  ";
    WriteOut(sqrt(pgmg.get_dat(i,i)),"%5.3f");
    cout<<"\n";
  };

  pgmg=pgmg+cal_cov.GetCovariance("Relative")+exp_cov.GetCovariance("Relative");
  pgmg=pgmg.inverse();

  real kaisp=(diff-kai)*(pgmg*(diff-kai));
  cout.setf(ios::showpoint);
  cout.precision(5);
  cout<<"\n*---------------------------------------------------------------------\n";
  cout<<"* Chi-square over DOF after adjustment\n*\n";
  cout<<"*    DOF : "<<sens_size<<"\n";
  cout<<"*                    "<<kaisp/sens_size<<"\n";
  cout<<"*---------------------------------------------------------------------\n";

  GroupData1D ce_n1(sens_size);
  for(int i=0;i<sens_size;i++){
    ce_n1.put_data(i,1.-ce.GetValue().get_dat(i));
  };
  GroupData1D tmp1d=inv*ce_n1;

  //++++++(Fission yield change)+++++++
  int fisnucnum=icov.GetFissionNuclideNumber();
  for(int i=0;i<fisnucnum;i++){
    int size=icov.GetSize(i);//size of covarinace matrix
    string fisnucname=icov.GetFissionNuclideName(i);
    int mt=18000000+midt.ID(fisnucname);
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=icov.GetID(i,k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));	
      };
    };
    
    GroupData1D change(size);//change of independent fission yield
    GroupData2D M=icov.GetCovarianceMatrix(i);
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in fission yield\n";
    cout<<"****************************************\n";
    cout<<"*  Fission nuclide: "<<fisnucname<<"\n";
    cout<<"*Nucl.      Yield.      Independent      Rel.\n";
    cout<<"*          Rel.Chg.       yield        Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(icov.GetID(i,ii));
      real chng=change.get_dat(ii);
      real yield=icov.GetIndependentYield(i,ii);
      real stdev=icov.GetIndependentYieldStandardDeviation(i,ii);
      if(yield>0.){
        WriteOut(name,6);
        WriteOut(chng,"%14.3e");
        WriteOut(yield,"%14.3e");
        WriteOut(stdev,"%14.3e");
        cout<<"\n";
	if(chng!=0.){
          fout_xschg<<icov.GetID(i,ii)<<"\n";
          fout_xschg<<18000000+midt.ID(fisnucname)<<"\n";
          fout_xschg<<chng<<"\n";
	};
      };
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";
    cout<<"*  Fission nuclide: "<<fisnucname<<"\n";

    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=icov.GetID(i,kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";
  };

  //++++++(Half life change)+++++++++++
    int size=hcov.GetSize();//size of covarinace matrix
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=hcov.GetID(k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,8888);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));
      };
    };
    
    GroupData1D change(size);//change of half life
    GroupData2D M=hcov.GetCovarianceMatrix();
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in half life\n";
    cout<<"****************************************\n";
    cout<<"*Nuclide   HalfLife-Rel.Chg. Half life[s]  Rel.Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(hcov.GetID(ii));
      real chng=change.get_dat(ii);
      real halflife=hcov.GetHalfLife(ii);
      real stdev=hcov.GetHalfLifeStandardDeviation(ii);
      if(halflife!=0.){
        printf("%6s %16.5e %16.5e %16.5e\n",name.c_str(),chng,halflife,stdev);
        if(chng!=0.){
          fout_xschg<<hcov.GetID(ii)<<"\n";
          fout_xschg<<8888<<"\n";
          fout_xschg<<chng<<"\n";
	};
      };
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";
    cout<<"*  Half life\n";
    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=hcov.GetID(kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,8888);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";

  //+++++++++++++++++++++++++++++++++++

  //++++++(Branching ratio change)+++++++
    int nucnum=bcov.GetNucNum();//number of nuclide
    int size_b=bcov.GetSize();//size of covarinace matrix
    GroupData1D ggmgr_b(size_b);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr_b.set_zero();
    for(int j=0;j<sens_size;j++){
      int counter=0;
      for(int k=0;k<nucnum;k++){
	int mat=bcov.GetID(k);
	int channel=bcov.GetChannel(k);
	for(int ii=0;ii<channel;ii++){
	  int mt=88880+ii;
	  int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	  real tmp;
	  if(data==-1){
	    tmp=0.;
	  }else{
	    tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	  };
	  ggmgr_b.add_data(counter,tmp*tmp1d.get_dat(j));
	  counter++;
	};
      };
    };

    GroupData1D change_b(size_b);//change of branching ratio
    GroupData2D M_b=bcov.GetCovarianceMatrix();
    M_b.ReducedForm();
    change_b=M_b*ggmgr_b;
    cout<<"****************************************\n";
    cout<<"* Change in branching ratio\n";
    cout<<"****************************************\n";
    cout<<"*Nuclide\n";
cout<<"*Channel-No. BranchingRatio-Rel.Chg. BranchingRatio  Rel.Std.Dev.\n";
    int counter=0;
    for(int ii=0;ii<nucnum;ii++){
      string name=midt.GetName(bcov.GetID(ii));
      int channel=bcov.GetChannel(ii);      
      if(channel>1)cout<<name<<":\n";
      for(int iii=0;iii<channel;iii++){
	if(channel>1){
	  real chng=change_b.get_dat(counter);
	  real ratio=bcov.GetRatio(ii,iii);
	  real stdev=bcov.GetBranchingRatioStandardDeviation(counter);
	  printf("%6d %23.5e %18.5e %14.5e\n",iii,chng,ratio,stdev);
          fout_xschg<<bcov.GetID(ii)<<"\n";
          fout_xschg<<88880+iii<<"\n";
          fout_xschg<<chng<<"\n";
	};
	counter++;
      };
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";
    cout<<"*  Branching ratio\n";
    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      int counter=0;
      for(int kk=0;kk<nucnum;kk++){
	int mat=bcov.GetID(kk);
	int channel=bcov.GetChannel(kk);
	for(int kkk=0;kkk<channel;kkk++){
	  int mt=88880+kkk;
	  int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	  if(data==-1){
	    tmp+=0.;
	  }else{
	    tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change_b.get_dat(counter);
	  };
	  counter++;
	};
      };      
      cout.setf(ios::scientific);
      cout.precision(4);
      cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
      cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
      cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
      cout<<tmp<<"\n";
    };
    cout<<"\n";
    
    //+++++++++++++++++++++++++++++++++++++

  //++++++(Decay energy change)++++++++
  for(int i=0;i<3;i++){
    int size=dcov.GetSize();//size of covarinace matrix
    int mt=99990+i;
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=dcov.GetID(k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));	
      };
    };
    
    GroupData1D change(size);//change of decay energy
    GroupData2D M=dcov.GetCovarianceMatrix(i);
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in decay energy\n";
    cout<<"****************************************\n";
    if(i==0)cout<<"* Beta decay energy \n";
    if(i==1)cout<<"* Gamma decay energy \n";
    if(i==2)cout<<"* Alpha decay energy \n";
    cout<<"*Nuclide   Energy-Rel.Chg. Decay energy[eV]  Rel.Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(dcov.GetID(ii));
      real chng=change.get_dat(ii);
      real energy=dcov.GetDecayEnergy(ii,i);
      real stdev=dcov.GetDecayEnergyStandardDeviation(ii,i);
      if(energy!=0.){
        printf("%6s %16.5e %16.5e %16.5e\n",name.c_str(),chng,energy,stdev);
        fout_xschg<<dcov.GetID(ii)<<"\n";
        fout_xschg<<99990+i<<"\n";
        fout_xschg<<chng<<"\n";
      };
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";    
    if(i==0)cout<<"* Beta decay energy \n";
    if(i==1)cout<<"* Gamma decay energy \n";
    if(i==2)cout<<"* Alpha decay energy \n";

    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=dcov.GetID(kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";
  };

  //+++++++++++++++++++++++++++++++++++

  fout_xschg<<" -1\n";
  fout_xschg.close();

};

void UncertaintyCalculation::DoUncertaintyLibraryAdjustmentForBurn
  (SensitivityContainer &sens_con, SensitivityContainer &sens_con_im,
   ParameterCovariance &cal_cov, ParameterCovariance &exp_cov,
   Parameters &cal_val, Parameters &cal_val_im,
   Parameters &exp_val, Library &lib,
   LibraryCovariance &lib_cov, IndependentYieldCovariance &icov,
   HalfLifeCovariance &hcov, DecayEnergyCovariance &dcov,
   BranchingRatioCovariance &bcov)
{
  int sens_size=sens_con.GetParameterList()->GetSize();

  cout<<"****************************************\n";
  cout<<"* Cross section adjustment              \n";
  cout<<"*       by UncertaintyCalculation       \n";
  cout<<"****************************************\n";
  cout<<"* Number of integral data    : "<<sens_size<<"\n";
  cout<<"****************************************\n";

  ParameterCovariance gmg_xs=GetGMG(lib,lib_cov,cal_val,sens_con);
  UncertaintyCalculationForYieldDecay ucyd;
  GroupData2D gmg_nd=gmg_xs.GetCovariance("Relative")
    +ucyd.GetFissionYieldGMG(sens_con,icov)
    +ucyd.GetHalfLifeGMG(sens_con,hcov)
    +ucyd.GetDecayEnergyGMG(sens_con,dcov)
    +ucyd.GetBranchingRatioGMG(sens_con,bcov);

  GroupData2D inv(sens_size,sens_size);
  inv=gmg_nd
    +cal_cov.GetCovariance("Relative")
    +exp_cov.GetCovariance("Relative");
  
  Parameters ce(cal_val.GetParameterList());
  for(int i=0;i<sens_size;i++){
    ce.PutValue(i,cal_val.GetValue(i)/exp_val.GetValue(i));
  };

  ParameterCovariance total_cov(ce);
  total_cov.PutCov(inv,"Relative");

  cout<<"* Prior C/E values with total uncertainties\n";
  cout<<"*\n";
  for(int i=0;i<sens_size;i++){
    cout<<"*  "<<sens_con.GetParameterList()->GetCoreTag(i);
    cout<<","<<sens_con.GetParameterList()->GetCharaTag(i);
    cout<<","<<sens_con.GetParameterList()->GetStepTag(i);
    cout<<" : "<<ce.GetValue(i)<<" ("<<sqrt(inv.get_dat(i,i))<<")\n";
  };
  cout<<"****************************************\n";  

  cout<<"* Correlation matrix of GMG+Ve+Vm\n";
  cout<<"****************************************\n";  
  total_cov.GetCorrelationMatrix().show_self(false);
  cout<<"****************************************\n";  

  inv=inv.inverse();

  GroupData1D diff(sens_size);
  for(int i=0;i<sens_size;i++){
    diff.put_data(i,ce.GetValue(i)-1.0);
  };
  real kais=diff*(inv*diff);
  cout<<"*    Kai square = "<<kais<<" ("<<kais/sens_size<<")\n";

  // Only diagonal element of covariance matirx

  GroupData2D tt(sens_size,sens_size);
  for(int i=0;i<sens_size;i++){
    tt.put_data(i,i,inv.get_dat(i,i));
  };
  real kais2=diff*(tt*diff);
  cout<<"*    Kai square (diag) = "<<kais2<<" ("<<kais2/sens_size<<")\n";
  cout<<"****************************************\n";  

  GroupData1D kai(sens_size);
  kai=gmg_nd*(inv*diff);
  GroupData2D pgmg(sens_size,sens_size);
  pgmg=gmg_nd-gmg_nd*(inv*gmg_nd);

  cout<<"****************************************\n";
  cout<<"* Adjustment\n";
  cout<<"****************************************\n";
  cout<<"* Posterior C/E values with nuclear data induced uncertainties\n";
  cout<<"*\n";
  for(int i=0;i<sens_size;i++){
    cout<<"*  "<<sens_con.GetParameterList()->GetCoreTag(i);
    cout<<","<<sens_con.GetParameterList()->GetCharaTag(i);
    cout<<","<<sens_con.GetParameterList()->GetStepTag(i);
    cout.setf(ios::showpoint);
    cout.precision(5);
    //cout<<" : "<<ce.GetValue(i)-kai.get_dat(i);//Absolute?(Original)
    cout<<" : "<<ce.GetValue(i)*(1.-kai.get_dat(i));//Relative?
    cout<<" ("<<sqrt(pgmg.get_dat(i,i))<<")\n";
    cout.unsetf(ios::showpoint);
  };

  pgmg=pgmg+cal_cov.GetCovariance("Relative")+exp_cov.GetCovariance("Relative");
  pgmg=pgmg.inverse();

  /*
  for(int i=0;i<sens_size;i++){
    tmp.put_data(i,ce.GetValue(i)-1.0);
  };
  */
  //real kaisp=(ce.GetValue()-kai)*(pgmg*(ce.GetValue()-kai));
  real kaisp=(diff-kai)*(pgmg*(diff-kai));
  cout<<"****************************************\n";
  cout<<"*    Kai square = "<<kaisp<<" ("<<kaisp/sens_size<<")\n";
  cout<<"****************************************\n";
  //

  //++++++(Cross section change)+++++++
  vector<int> mmat;
  vector<int> mmt;
  int mmnn=0;
  int cov_size=lib_cov.GetSize();
  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=lib_cov.GetCrossSectionCovariance(i);
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    if(mmnn==0){
      mmat.push_back(mat1);
      mmt.push_back(mt1);
      mmnn++;
    }else{
      bool exist=false;
      for(int j=0;j<mmnn;j++){
	if(mat1==mmat[j]&&mt1==mmt[j])exist=true;
      };
      if(!exist){
	mmat.push_back(mat1);
	mmt.push_back(mt1);
	mmnn++;
      };
    };
    int mat2=cov.GetMat2();
    int mt2=cov.GetMt2();
    bool exist=false;
    for(int j=0;j<mmnn;j++){
      if(mat2==mmat[j]&&mt2==mmt[j])exist=true;
    };
    if(!exist){
      mmat.push_back(mat2);
      mmt.push_back(mt2);
      mmnn++;
    };
  };

  vector< vector<GroupData1D> > tmpmat(mmnn);
  for(int i=0;i<mmnn;i++){
    tmpmat[i].resize(sens_size);
    for(int j=0;j<sens_size;j++){
      tmpmat[i][j].put_imax(grp);
      tmpmat[i][j].set_zero();
    };
  };

  for(int i=0;i<cov_size;i++){
    CrossSectionCovariance cov=lib_cov.GetCrossSectionCovariance(i);
    GroupData2D covdata=cov.GetCovariance("Relative");
    int mat1=cov.GetMat1();
    int mt1 =cov.GetMt1();
    int mat2=cov.GetMat2();
    int mt2 =cov.GetMt2();
    int index1=-1;
    int index2=-1;
    for(int j=0;j<mmnn;j++){
      if(mat1==mmat[j]&&mt1==mmt[j])index1=j;
      if(mat2==mmat[j]&&mt2==mmt[j])index2=j;
    };
    for(int j=0;j<sens_size;j++){
      int tmp1=sens_con.GetSensitivityData(j).FindData1D(mat1,mt1);
      int tmp2=sens_con.GetSensitivityData(j).FindData1D(mat2,mt2);
      if(index1!=-1&&tmp2!=-1){
        tmpmat[index1][j]=tmpmat[index1][j]+
	  covdata*sens_con.GetSensitivityData(j).GetSensitivity1D(tmp2);
      };
      if(index2!=-1&&tmp1!=-1&&index2!=index1){
        tmpmat[index2][j]=tmpmat[index2][j]+
	  covdata.GetTransposedMatrix()*sens_con.GetSensitivityData(j).GetSensitivity1D(tmp1);
      };
    };
  };

  cout<<"****************************************\n";
  cout<<"* Change in cross section\n";
  cout<<"****************************************\n";
  GroupData1D ce_n1(sens_size);
  for(int i=0;i<sens_size;i++){
    ce_n1.put_data(i,1.-ce.GetValue().get_dat(i));
  };
  for(int i=0;i<mmnn;i++){ //
    int mat=mmat[i];
    int mt=mmt[i];
    int fd=lib_cov.FindData(mat,mt,mat,mt);
    int fd2=lib.FindData(mat,mt);
    if(fd!=-1&&fd2!=-1){
      cout<<"*  MAT: "<<mat<<"  MT: "<<mt<<"\n";
      cout<<"*grp Upper Eng.  XS-Rel.Chg.  XS        Rel.Std.Dev.\n";
      GroupData1D xschange(grp);
      for(int j=0;j<grp;j++){
        real tmp=0.;
        for(int k=0;k<sens_size;k++){
  	  tmp+=tmpmat[i][k].get_dat(j)*(inv*ce_n1).get_dat(k);
        };
        xschange.put_data(j,tmp);
        real xs=lib.GetCrossSection(mat,mt).get_dat(j);
        real stdev=lib_cov.GetCrossSectionCovariance(mat,mt,mat,mt).GetCovariance("Relative").get_diag().sqrt1d().get_dat(j);
        printf("%3d %11.4e %10.5f %10.5f %10.5f\n",j,ebnd.get_dat(j),tmp,xs,stdev);
      };
      cout<<"****************************************\n";
      cout<<" ( Impact on responce/integral parameter )\n";          
      cout<<"*  MAT: "<<mat<<"  MT: "<<mt<<"\n";
      for(int k=0;k<sens_size;k++){
        int ttt=sens_con.GetSensitivityData(k).FindData1D(mat,mt);
        if(ttt!=-1){
	  real dce=sens_con.GetSensitivityData(k).GetSensitivity1D(ttt)*xschange;
	  cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	  cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	  cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	  cout<<dce<<"\n";
	};
      };
      cout<<"\n";
    };
  };
  //+++++++++++++++++++++++++++++++++++

  GroupData1D tmp1d=inv*ce_n1;

  //++++++(Fission yield change)+++++++
  int fisnucnum=icov.GetFissionNuclideNumber();
  for(int i=0;i<fisnucnum;i++){
    int size=icov.GetSize(i);//size of covarinace matrix
    string fisnucname=icov.GetFissionNuclideName(i);
    int mt=18000000+midt.ID(fisnucname);
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=icov.GetID(i,k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));	
      };
    };
    
    GroupData1D change(size);//change of independent fission yield
    GroupData2D M=icov.GetCovarianceMatrix(i);
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in fission yield\n";
    cout<<"****************************************\n";
    cout<<"*  Fission nuclide: "<<fisnucname<<"\n";
    cout<<"*Nuclide   Yield-Rel.Chg. Independent yield   Rel.Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(icov.GetID(i,ii));
      real chng=change.get_dat(ii);
      real yield=icov.GetIndependentYield(i,ii);
      real stdev=icov.GetIndependentYieldStandardDeviation(i,ii);
      printf("%6s %16.5e %16.5e %16.5e\n",name.c_str(),chng,yield,stdev);
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";              
    cout<<"*  Fission nuclide: "<<fisnucname<<"\n";

    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=icov.GetID(i,kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";
  };

  //++++++(Decay energy change)++++++++
  for(int i=0;i<3;i++){
    int size=dcov.GetSize();//size of covarinace matrix
    int mt=99990+i;
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=dcov.GetID(k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));	
      };
    };
    
    GroupData1D change(size);//change of decay energy
    GroupData2D M=dcov.GetCovarianceMatrix(i);
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in decay energy\n";
    cout<<"****************************************\n";
    if(i==0)cout<<"* Beta decay energy \n";
    if(i==1)cout<<"* Gamma decay energy \n";
    if(i==2)cout<<"* Alpha decay energy \n";
    cout<<"*Nuclide   Energy-Rel.Chg. Decay energy[eV]  Rel.Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(dcov.GetID(ii));
      real chng=change.get_dat(ii);
      real energy=dcov.GetDecayEnergy(ii,i);
      real stdev=dcov.GetDecayEnergyStandardDeviation(ii,i);
      if(energy!=0.)printf("%6s %16.5e %16.5e %16.5e\n",name.c_str(),chng,energy,stdev);
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";                  
    if(i==0)cout<<"* Beta decay energy \n";
    if(i==1)cout<<"* Gamma decay energy \n";
    if(i==2)cout<<"* Alpha decay energy \n";

    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=dcov.GetID(kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";
  };

  //+++++++++++++++++++++++++++++++++++

  //++++++(Half life change)+++++++++++
    int size=hcov.GetSize();//size of covarinace matrix
    GroupData1D ggmgr(size);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr.set_zero();
    for(int j=0;j<sens_size;j++){
      for(int k=0;k<size;k++){
	int mat=hcov.GetID(k);
	int data=sens_con.GetSensitivityData(j).FindData0D(mat,8888);
	real tmp;
	if(data==-1){
	  tmp=0.;
	}else{
	  tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	};
	ggmgr.add_data(k,tmp*tmp1d.get_dat(j));
      };
    };
    
    GroupData1D change(size);//change of half life
    GroupData2D M=hcov.GetCovarianceMatrix();
    M.ReducedForm();
    change=M*ggmgr;
    
    cout<<"****************************************\n";
    cout<<"* Change in half life\n";
    cout<<"****************************************\n";
    cout<<"*Nuclide   HalfLife-Rel.Chg. Half life[s]  Rel.Std.Dev.\n";
    for(int ii=0;ii<size;ii++){
      string name=midt.GetName(hcov.GetID(ii));
      real chng=change.get_dat(ii);
      real halflife=hcov.GetHalfLife(ii);
      real stdev=hcov.GetHalfLifeStandardDeviation(ii);
      if(halflife!=0.)printf("%6s %16.5e %16.5e %16.5e\n",name.c_str(),chng,halflife,stdev);
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";                      
    cout<<"*  Half life\n";
    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      for(int kk=0;kk<size;kk++){
	int mat=hcov.GetID(kk);
	int data=sens_con.GetSensitivityData(k).FindData0D(mat,8888);
	if(data==-1){
	  tmp+=0.;
	}else{
	  tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change.get_dat(kk);
	};
      };
      cout.setf(ios::scientific);
      cout.precision(4);
	cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
	cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
	cout<<tmp<<"\n";
    };
    cout<<"\n";

  //+++++++++++++++++++++++++++++++++++


  //++++++(Branching ratio change)+++++++
    int nucnum=bcov.GetNucNum();//number of nuclide
    int size_b=bcov.GetSize();//size of covarinace matrix
    GroupData1D ggmgr_b(size_b);//G^t(GMG^t+Ve+Vm)^-1(Re-Rc(T)) 
    ggmgr_b.set_zero();
    for(int j=0;j<sens_size;j++){
      int counter=0;
      for(int k=0;k<nucnum;k++){
	int mat=bcov.GetID(k);
	int channel=bcov.GetChannel(k);
	for(int ii=0;ii<channel;ii++){
	  int mt=88880+ii;
	  int data=sens_con.GetSensitivityData(j).FindData0D(mat,mt);
	  real tmp;
	  if(data==-1){
	    tmp=0.;
	  }else{
	    tmp=sens_con.GetSensitivityData(j).GetSensitivity0D(data);
	  };
	  ggmgr_b.add_data(counter,tmp*tmp1d.get_dat(j));
	  counter++;
	};
      };
    };

    GroupData1D change_b(size_b);//change of branching ratio
    GroupData2D M_b=bcov.GetCovarianceMatrix();
    M_b.ReducedForm();
    change_b=M_b*ggmgr_b;
    cout<<"****************************************\n";
    cout<<"* Change in branching ratio\n";
    cout<<"****************************************\n";
    cout<<"*Nuclide\n";
cout<<"*Channel-No. BranchingRatio-Rel.Chg. BranchingRatio  Rel.Std.Dev.\n";
    int counter=0;
    for(int ii=0;ii<nucnum;ii++){
      string name=midt.GetName(bcov.GetID(ii));
      int channel=bcov.GetChannel(ii);      
      if(channel>1)cout<<name<<":\n";
      for(int iii=0;iii<channel;iii++){
	if(channel>1){
	  real chng=change_b.get_dat(counter);
	  real ratio=bcov.GetRatio(ii,iii);
	  real stdev=bcov.GetBranchingRatioStandardDeviation(counter);
	  printf("%6d %23.5e %18.5e %14.5e\n",iii,chng,ratio,stdev);
	};
	counter++;
      };
    };
    cout<<"****************************************\n";
    cout<<" ( Impact on responce/integral parameter )\n";                          
    cout<<"*  Branching ratio\n";
    for(int k=0;k<sens_size;k++){
      real tmp=0.;
      int counter=0;
      for(int kk=0;kk<nucnum;kk++){
	int mat=bcov.GetID(kk);
	int channel=bcov.GetChannel(kk);
	for(int kkk=0;kkk<channel;kkk++){
	  int mt=88880+kkk;
	  int data=sens_con.GetSensitivityData(k).FindData0D(mat,mt);
	  if(data==-1){
	    tmp+=0.;
	  }else{
	    tmp+=sens_con.GetSensitivityData(k).GetSensitivity0D(data)*change_b.get_dat(counter);
	  };
	  counter++;
	};
      };      
      cout.setf(ios::scientific);
      cout.precision(4);
      cout<<"    + "<<sens_con.GetParameterList()->GetCoreTag(k)<<" ";
      cout<<sens_con.GetParameterList()->GetCharaTag(k)<<" ";
      cout<<sens_con.GetParameterList()->GetStepTag(k)<<" ";
      cout<<tmp<<"\n";
    };
    cout<<"\n";
    
    //+++++++++++++++++++++++++++++++++++++
    
    
    //++++++(Imaginary case calculation)+++++++
    //G'MG'(G':imaginary case G)
    GroupData2D gmg_xs_im=GetGMG(lib,lib_cov,cal_val_im,sens_con_im).GetCovariance("Relative");
    GroupData2D gmg_fy_im=ucyd.GetFissionYieldGMG(sens_con_im,icov);
    GroupData2D gmg_hl_im=ucyd.GetHalfLifeGMG(sens_con_im,hcov);
    GroupData2D gmg_de_im=ucyd.GetDecayEnergyGMG(sens_con_im,dcov);
    GroupData2D gmg_br_im=ucyd.GetBranchingRatioGMG(sens_con_im,bcov);
    GroupData2D gmg_nd_im=gmg_xs_im+gmg_fy_im+gmg_hl_im+gmg_de_im+gmg_br_im;
    //GMG'
    GroupData2D gmg_xs_im1=GetGMG(lib_cov,cal_val,cal_val_im,sens_con,sens_con_im);
    GroupData2D gmg_fy_im1=ucyd.GetFissionYieldGMG(sens_con,icov,sens_con_im);
    GroupData2D gmg_hl_im1=ucyd.GetHalfLifeGMG(sens_con,hcov,sens_con_im);
    GroupData2D gmg_de_im1=ucyd.GetDecayEnergyGMG(sens_con,dcov,sens_con_im);
    GroupData2D gmg_br_im1=ucyd.GetBranchingRatioGMG(sens_con,bcov,sens_con_im);
    GroupData2D gmg_nd_im1=gmg_xs_im1+gmg_fy_im1+gmg_hl_im1+gmg_de_im1+gmg_br_im1;
    //G'MG
    GroupData2D gmg_xs_im2=GetGMG(lib_cov,cal_val_im,cal_val,sens_con_im,sens_con);
    GroupData2D gmg_fy_im2=ucyd.GetFissionYieldGMG(sens_con_im,icov,sens_con);
    GroupData2D gmg_hl_im2=ucyd.GetHalfLifeGMG(sens_con_im,hcov,sens_con);
    GroupData2D gmg_de_im2=ucyd.GetDecayEnergyGMG(sens_con_im,dcov,sens_con);
    GroupData2D gmg_br_im2=ucyd.GetBranchingRatioGMG(sens_con_im,bcov,sens_con);
    GroupData2D gmg_nd_im2=gmg_xs_im2+gmg_fy_im2+gmg_hl_im2+gmg_de_im2+gmg_br_im2;

    int sens_size_im=sens_con_im.GetParameterList()->GetSize();  
    GroupData1D kai_im(sens_size_im);
    kai_im=gmg_nd_im2*(inv*diff);
    GroupData2D pgmg_im(sens_size_im,sens_size_im);
    pgmg_im=gmg_nd_im-gmg_nd_im2*(inv*gmg_nd_im1);

    cout<<"****************************************\n";
    cout<<"* Adjustment results in imaginary case\n";
    cout<<"****************************************\n";
    cout<<"* Number of imaginary case    : "<<sens_size_im<<"\n";
    cout<<"****************************************\n";
    cout<<"* Prior and posterior calculation values with nuclear data induced relative uncertainties\n";
    cout<<"*\n";
    for(int i=0;i<sens_size_im;i++){
      cout<<"* "<<sens_con_im.GetParameterList()->GetCoreTag(i);
      cout<<","<<sens_con_im.GetParameterList()->GetCharaTag(i);
      cout<<","<<sens_con_im.GetParameterList()->GetStepTag(i);
      cout.setf(ios::showpoint);
      cout.precision(5);
      cout<<": "<<cal_val_im.GetValue(i);
      cout<<" (";
      cout<<sqrt(gmg_nd_im.get_dat(i,i));
      cout<<") -> ";
      cout<<cal_val_im.GetValue(i)*(1.-kai_im.get_dat(i));
      cout<<" (";
      cout<<sqrt(pgmg_im.get_dat(i,i));
      cout<<")\n";
    
      cout.unsetf(ios::showpoint);
    };
    cout<<"\n";
    //+++++++++++++++++++++++++++++++++++++++++
  
};

/*
void UncertaintyCalculation::CondSensitivity(int cgrp, int *bgrp, string name)
{
  GroupData1D dummy(cgrp);
  cout<<cgrp<<"\n";
  for(int i=0;i<numsens;i++){
    if(sensname[i]==name){
      cout<<sensmat[i]<<"\n";
      cout<<sensmt[i]<<"\n";
      dummy.set_zero();
      dummy=GetSensitivity1D(sensmat[i],sensmt[i],name).CondSum(cgrp,bgrp);
      for(int j=0;j<cgrp;j++){
	cout<<dummy.get_dat(j)<<"\n";
      };
    };
  };
};

void UncertaintyCalculation::MixCrossSection
(int mixn,int *matlist,string libname,real *weight,int newmat)
{
  ofstream fout;
  fout.open("output",ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  fout<<grp<<"\n";

  int mtn=6;
  int mtlist[]={2,102,4,16,1,251};

  real tmp=0.;
  for(int i=0;i<mixn;i++){
    tmp+=weight[i];
  };
  for(int i=0;i<mixn;i++){
    weight[i]/=tmp;
  };

  for(int i=0;i<mtn;i++){
    bool exist=true;
    for(int j=0;j<mixn;j++){
      if(!ExistCrossSection(matlist[j],mtlist[i],libname))exist=false;
    };
    if(exist){
      fout<<newmat<<"\n";
      fout<<mtlist[i]<<"\n";
      GroupData1D tmp(grp);
      GroupData1D tmp3(grp);
      tmp.set_zero();
      tmp3.set_zero();
      for(int j=0;j<mixn;j++){
	if(mtlist[i]!=251){
	  tmp=tmp+GetCrossSection(matlist[j],mtlist[i],libname)*weight[j];
	}else{
	  GroupData1D tmp2(grp);
	  tmp2=GetCrossSection(matlist[j],2,libname)*weight[j];
	  tmp=tmp+GetCrossSection(matlist[j],mtlist[i],libname).mult(tmp2);
	  tmp3=tmp3+tmp2;
	};
      };
      if(mtlist[i]==251){
	tmp=tmp/tmp3;
      };
      for(int i=0;i<grp;i++){
	fout<<tmp.get_dat(i)<<"\n";
      };
    };
  };
  fout.close();
};

void UncertaintyCalculation::CovMultiply(real v,int mat,int mt)
{
  string mm;
  for(int i=0;i<numcov;i++){
    int mat1=XSCov[i].GetMat1();
    int mt1 =XSCov[i].GetMt1();
    int mat2=XSCov[i].GetMat2();
    int mt2 =XSCov[i].GetMt2();
    if(mat1==mat&&mt1==mt)XSCov[i].Multiply(v);
    if(mat2==mat&&mt2==mt)XSCov[i].Multiply(v);
  };
};

void UncertaintyCalculation::CovMultiply(real v,int mat,int mt,int g,string nn)
{
  int mats1,mts1,mats2,mts2;
  string mm;
  for(int i=0;i<icov;i++){
    mats1=CovDat[i].GetMat1();
    mts1 =CovDat[i].GetMt1();
    mats2=CovDat[i].GetMat2();
    mts2 =CovDat[i].GetMt2();
    mm   =CovDat[i].GetName();
    if(mats1==mat&&mts1==mt&&mm==nn){
      int ig=CovDat[i].GetGrp();
      for(int j=0;j<ig;j++){
        real org=CovDat[i].GetCov().get_dat(g,j);
	CovDat[i].GetCov().put_data(g,j,org*v);
      };
    };
    if(mats2==mat&&mts2==mt&&mm==nn){
      int ig=CovDat[i].GetGrp();
      for(int j=0;j<ig;j++){
	real org=CovDat[i].GetCov().get_dat(j,g);
	CovDat[i].GetCov().put_data(j,g,org*v);
      };
    };
  };
};

void DataSet::CondSensitivity(int cgrp, int *bgrp, char *name)
{
  cout<<cgrp<<"\n";
  for(int i=0;i<isens;i++){
    if(GetSensitivity1D(i).GetName()==name){
      cout<<GetSensitivity1D(i).GetMat()<<"\n";
      cout<<GetSensitivity1D(i).GetMt()<<"\n";
      GetSensitivity1D(i).Cond(cgrp,bgrp);
    };
  };
};

void DataSet::ShowDiffLibrary(int mat,int mt,char *para1,char *para2)
{
  cout<<"*** Relative difference in cross section ***\n";
  cout<<"   "<<para1<<" - "<<para2<<"\n";
  cout<<"   (mat:"<<mat<<", mt:"<<mt<<")\n";
  int ind1,ind2;
  ind1=GetCrossSectionIndex(mat,mt,para1);
  ind2=GetCrossSectionIndex(mat,mt,para2);
  if(ind1==-1||ind2==-1)return;

  int g1=GetCrossSection(ind1).GetGrp();
  int g2=GetCrossSection(ind2).GetGrp();
  if(g1!=g2)return;

  for(int i=0;i<g1;i++){
    real xs1=GetCrossSection(ind1).GetDAT().get_dat(i);
    real xs2=GetCrossSection(ind2).GetDAT().get_dat(i);
    printf("%11.4e  %11.5f\n",ebnd.get_dat(i),(xs1-xs2)/xs2);
  };
};
*/

void UncertaintyCalculation::ShowCrossSectionStandardDeviation(LibraryCovariance &libcov,int mat,int mt,bool excel)
{
  int id=libcov.FindData(mat,mt,mat,mt);
  if(id==-1){
    cout<<"# Warning!\n";
    cout<<"# No variance information for mat/mt="<<mat<<"/"<<mt<<"\n";
    exit(0);
  };

  cout<<"# Upper     Relative\n";
  cout<<"# energy    standard Deviation\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  GroupData1D stddev=libcov.GetCrossSectionCovariance(mat,mt,mat,mt).GetStandardDeviation();
  for(int i=0;i<grp;i++){
    real e0=ebnd.get_dat(i);
    real e1=ebnd.get_dat(i+1);
    cout<<e0<<" "<<stddev.get_dat(i)<<"\n";
    if(excel)cout<<e1<<" "<<stddev.get_dat(i)<<"\n";
  };
};

void UncertaintyCalculation::PutCvalueFromSensitivity(SensitivityContainer &sens_con,Parameters &cal_val)
{
  int sz=sens_con.GetSize();
  for(int i=0;i<sz;i++){
    string core=sens_con.GetParameterList()->GetCoreTag(i);
    string chara=sens_con.GetParameterList()->GetCharaTag(i);
    int step=sens_con.GetParameterList()->GetStepTag(i);
    real value=sens_con.GetSensitivityData(i).GetValue();
    cal_val.PutValue(core,chara,value,step);
  };
};

void UncertaintyCalculation::ShowCorrelationForXYPlot(LibraryCovariance &libcov,int mat,int mt)
{
  GroupData2D cor=libcov.GetCrossSectionCovariance(mat,mt,mat,mt).GetCorrelationMatrix();
  
  for(int i=0;i<grp;i++){
    for(int k=0;k<2;k++){
    for(int j=0;j<grp;j++){
      real dat=cor.get_dat(i,j);
      //if(i>=j)dat=0.;
      cout<<ebnd.get_dat(i+k)<<" "<<ebnd.get_dat(j)<<" "<<dat<<"\n";
      cout<<ebnd.get_dat(i+k)<<" "<<ebnd.get_dat(j+1)<<" "<<dat<<"\n";
    };
    cout<<"\n";
    };
  };
};

void UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXYPlot
(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,
 string core,string chara,int step,int mat1,int mt1,int mat2,int mt2, string filename)
{
  string outfile="./"+filename+"_"+IntToString(mat1)+"_"+IntToString(mt1);
  if(mat2!=mat1||mt2!=mt2)outfile+="_"+IntToString(mat2)+"_"+IntToString(mt2);

  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"# Error in UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXSPlot.\n";
    cout<<"# There is no sensitivity data ...\n";
    exit(0);
  };
  SensitivityData sens=sens_con.GetSensitivityData(core,chara,step);

  int ids1=sens.FindData1D(mat1,mt1);
  int ids2=sens.FindData1D(mat2,mt2);
  if(ids1==-1||ids2==-1){
    cout<<"# Error in UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXSPlot.\n";
    cout<<"# There is no sensitivity data ...\n";
    exit(0);
  };

  int idcov=lib_cov.FindData(mat1,mt1,mat2,mt2);
  if(idcov==-1){
    cout<<"# Error in UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXSPlot.\n";
    cout<<"# There is no covariance data ("<<mat1<<","<<mt1<<"/"<<mat2<<","<<mt2<<"\n";
    exit(0);
  };

  real totvar=CalCrossSectionUncertainty(sens_con,lib_cov,core,chara,step,1.,false);

  GroupData1D s1=sens.GetSensitivity1D(ids1);
  GroupData1D s2=sens.GetSensitivity1D(ids2);
  GroupData2D cov=lib_cov.GetCrossSectionCovariance(idcov).GetCovariance("relative");

  ShowCrossSectionUncertaintyComponentXYPlot(s1,s2,cov,totvar,outfile);
};

void UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXYPlot
(GroupData1D &s1, GroupData1D &s2, GroupData2D &covi, real factor, string outfile)
{
  GroupData2D cov(grp,grp);

  for(int g=0;g<grp;g++){
    for(int g2=0;g2<grp;g2++){
      real org=covi.get_dat(g,g2);
      cov.put_data(g,g2,org*s1.get_dat(g)*s2.get_dat(g2)*factor);
    };
  };

  ofstream fout;
  fout.open(outfile.data(),ios::out);
  if(fout.fail()){
    cout<<"# Error in UncertaintyCalculation::ShowCrossSectionUncertaintyComponentXSPlot.\n";
    cout<<"# The following file cannot be open.\n";
    cout<<"# file name : "<<outfile<<"\n";
    exit(0);
  };

  for(int i=0;i<grp;i++){
    for(int k=0;k<2;k++){
    for(int j=0;j<grp;j++){
      real dat=cov.get_dat(i,j);
      fout<<ebnd.get_dat(i+k)<<" "<<ebnd.get_dat(j)<<" "<<dat<<"\n";
      fout<<ebnd.get_dat(i+k)<<" "<<ebnd.get_dat(j+1)<<" "<<dat<<"\n";
    };
    fout<<"\n";
    };
  };
  fout.close();

};

void UncertaintyCalculation::CalVarianceReductionFactor
(SensitivityContainer &sens_con, LibraryCovariance &lib_cov,string core,string chara,int step,real mean)
{
  real para_var=CalCrossSectionUncertainty(sens_con,lib_cov,core,chara,step,1.,false);

  int tmp=sens_con.GetParameterList()->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"# Error in UncertaintyCalculation::MIMP.\n";
    cout<<"# There is no sensitivity data ...\n";
    exit(0);
  };
  SensitivityData sens=sens_con.GetSensitivityData(core,chara,step);

  vector<int> matid_store;
  vector<int> mtid_store;
  vector<real> maxvrf_store;

  int sz=sens.GetSize1D();
  for(int ii=0;ii<sz;ii++){
    int mat=sens.GetMatList1D(ii);
    int mt=sens.GetMtList1D(ii);
    int id=lib_cov.FindData(mat,mt,mat,mt);
    if(id!=-1){
      GroupData2D cov=lib_cov.GetCrossSectionCovariance(id).GetCovariance("Relative");
      vector<real> val(grp,0.);
      for(int i1=0;i1<sz;i1++){
        int mat1=sens.GetMatList1D(i1);
        int mt1=sens.GetMtList1D(i1);
        int id11=lib_cov.FindData(mat,mt,mat1,mt1);
        int id12=lib_cov.FindData(mat1,mt1,mat,mt);
        if(id11!=-1||id12!=-1){
          GroupData2D cov1;
          if(id11!=-1){
  	    cov1=lib_cov.GetCrossSectionCovariance(id11).GetCovariance("Relative");
	  }else{
	    cov1=lib_cov.GetCrossSectionCovariance(id12).GetCovariance("Relative").GetTransposedMatrix();
	  };
          GroupData1D s1=sens.GetSensitivity1D(i1);
          for(int i2=0;i2<sz;i2++){
            int mat2=sens.GetMatList1D(i2);
            int mt2=sens.GetMtList1D(i2);
            int id21=lib_cov.FindData(mat,mt,mat2,mt2);
            int id22=lib_cov.FindData(mat2,mt2,mat,mt);
	    if(id21!=-1||id22!=-1){
              GroupData2D cov2;
              if(id21!=-1){
  	        cov2=lib_cov.GetCrossSectionCovariance(id21).GetCovariance("Relative");
	      }else{
	        cov2=lib_cov.GetCrossSectionCovariance(id22).GetCovariance("Relative").GetTransposedMatrix();
  	      };
              GroupData1D s2=sens.GetSensitivity1D(i2);
              for(int g=0;g<grp;g++){
		real cv=cov.get_dat(g,g);
		if(cv>0.){
                  for(int g1=0;g1<grp;g1++){
                    real s1v=s1.get_dat(g1);
  	  	    for(int g2=0;g2<grp;g2++){
		      val[g]+=s1v*s2.get_dat(g2)*cov1.get_dat(g,g1)*cov2.get_dat(g,g2)/cv;
		    };
		  };
	        };
	      };
	    };
	  };
	};
      };

      real max=0.;
      for(int g=0;g<grp;g++){
        real tmp=val[g]/para_var;
	if(tmp>max)max=tmp;
      };
      if(max>mean){
	matid_store.push_back(mat);
	mtid_store.push_back(mt);
	maxvrf_store.push_back(max);
        cout<<"# "<<mat<<"/"<<mt<<"\n";
        for(int g=0;g<grp;g++){
  	  //cout<<g<<" "<<val[g]<<" "<<sqrt(val[g])<<"\n";
	  cout<<ebnd.get_dat(g)<<" "<<val[g]/para_var<<"\n";
        };
        cout<<"\n\n";
      };
    };
  };

  cout<<"#\n";
  cout<<"# Summary of variance reduction factor\n";
  cout<<"#\n";
  int sz2=matid_store.size();
  for(int i=0;i<sz2;i++){
    cout<<"#  ";
    WriteOut(i,2);
    cout<<" ";
    WriteOut(matid_store[i],7);
    cout<<" ";
    WriteOut(mtid_store[i],3);
    cout<<" "<<maxvrf_store[i]<<"\n";
  };


};

//+++(kawamoto)+++

void UncertaintyCalculationForYieldDecay::CalFissionYieldUncertainty
(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,string core,string chara,int step)
{
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  int fisnucnum=cov.GetFissionNuclideNumber();
  
  cout<<"********************************************************\n";
  cout<<"* Fission yield-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  real total=0;
  vector<string> name_store;
  for(int i=0;i<fisnucnum;i++){
    int size=cov.GetSize(i);
    string fisnucname=cov.GetFissionNuclideName(i);
    for(int j=0;j<name_store.size();j++){
      if(name_store[j]==fisnucname){
	cout<<"# Error in UncertaintyCalculationForYieldDecay::CalFissionYieldUncertainty.\n";
	cout<<"# Duplicate fission yield covariance are defined.\n";
	exit(0);
      };
    };
    name_store.push_back(fisnucname);
    GroupData1D G;
    G.put_imax(size);
    for(int j=0;j<size;j++){
      int mat=cov.GetID(i,j);
      int mt=18000000+midt.ID(fisnucname);
      int data=sens_data.FindData0D(mat,mt);
      if(data==-1){
	G.put_data(j,0.);
      }else{
	real tmp1=sens_data.GetSensitivity0D(data);
	real tmp2=tmp1;
	G.put_data(j,tmp2);
      };
    };

    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    real unc=G*(M*G);

    cout<<"*  Fission Yield from "<<fisnucname<<"    "<<unc<<" ( "<<sqrt(unc)<<" )\n";
    total=total+unc;
  };
  
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";
};

void UncertaintyCalculationForYieldDecay::CalFissionYieldUncertaintyDetail
(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,string core,string chara,int step,real mean)
{
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  int fisnucnum=cov.GetFissionNuclideNumber();
  
  cout<<"********************************************************\n";
  cout<<"* Fission yield-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"*                 value   "<<sens_data.GetValue()<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  real total=0;
  vector<string> name_store;
  for(int i=0;i<fisnucnum;i++){
    int size=cov.GetSize(i);
    GroupData2D M=cov.GetCovarianceMatrix(i);
    
    for(int j=0;j<size;j++){
      real tmp=M.get_dat(j,j);
      //if(tmp<0.)cout<<"# Negative variance value is detected : "<<tmp<<"\n";
    };

    M.ReducedForm();
    string fisnucname=cov.GetFissionNuclideName(i);
    for(int j=0;j<name_store.size();j++){
      if(name_store[j]==fisnucname){
	cout<<"# Error in UncertaintyCalculationForYieldDecay::CalFissionYieldUncertainty.\n";
	cout<<"# Duplicate fission yield covariance are defined.\n";
	exit(0);
      };
    };
    name_store.push_back(fisnucname);
    for(int ii=1;ii<300;ii++){
      int mass_counter=0;
      vector<real> mini_G_tmp;
      vector<int> mini_id;
      vector<int> mini_position;
      for(int j=0;j<size;j++){
	int mat=cov.GetID(i,j);
	int lz,la,li;
	midt.GetParameter(mat,lz,la,li);
	if(ii==la){
	  mass_counter++;
	  mini_id.push_back(mat);
	  mini_position.push_back(j);
	  int mt=18000000+midt.ID(fisnucname);
	  int data=sens_data.FindData0D(mat,mt);
	  if(data==-1){
	    mini_G_tmp.push_back(0.);
	  }else{
	    real tmp1=sens_data.GetSensitivity0D(data);
	    real tmp2=tmp1;
	    mini_G_tmp.push_back(tmp2);
	  };
	};
      };

      if(mass_counter!=0){
	GroupData1D mini_G;
	GroupData2D mini_cov;
	mini_cov.put_yx(mass_counter,mass_counter);
	mini_G.put_imax(mass_counter);
	for(int iii=0;iii<mass_counter;iii++){
	  mini_G.put_data(iii,mini_G_tmp[iii]);
	  for(int iiii=0;iiii<mass_counter;iiii++){
	    mini_cov.put_data(iii,iiii,M.get_dat(mini_position[iii],mini_position[iiii]));
	  };
	};
	real vali=mini_G*(mini_cov*mini_G);
	if(sqrt(vali)>=mean){
	  cout<<"* "<<fisnucname<<" A="<<ii<<" Mass Chain: "<<vali<<" ( "<<sqrt(vali)<<" ) \n";
	  //cout<<"   "<<ii<<"  "<<sqrt(vali)<<"\n";
	};
	total=total+vali;
      };
    };
  }; 
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";

};

void UncertaintyCalculationForYieldDecay::CalDecayEnergyUncertainty(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step){  
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Decay energy-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  real total=0;
  for(int i=0;i<3;i++){
    int size=cov.GetSize();
    GroupData1D G;
    G.put_imax(size);
    for(int j=0;j<size;j++){
      int mat=cov.GetID(j);
      int mt=99990+i;
      int data=sens_data.FindData0D(mat,mt);
      if(data==-1){
	G.put_data(j,0.);
      }else{
	real tmp1=sens_data.GetSensitivity0D(data);
	real tmp2=tmp1;
	G.put_data(j,tmp2);
      };
    };
    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    real unc=G*(M*G);

    string type;
    if(i==0)type=" Beta";
    if(i==1)type="Ganma";
    if(i==2)type="Alpha";
    
    cout<<"* "<<type<<" decay energy: "<<unc<<" ( "<<sqrt(unc)<<" )\n";
    total=total+unc;
};
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";
  
};

real UncertaintyCalculationForYieldDecay::CalDecayEnergyUncertaintyNuclideWise(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step,string nucname)
{  
  MATIDTranslator midt;
  int id=midt.ID(nucname);

  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);

  real total=0.;
  for(int i=0;i<3;i++){
    int size=cov.GetSize();
    for(int j=0;j<size;j++){
      int mat=cov.GetID(j);
      if(mat==id){
        int mt=99990+i;
        int data=sens_data.FindData0D(mat,mt);
        if(data!=-1){
          real ss=sens_data.GetSensitivity0D(data);
          real var=cov.GetCovarianceMatrix(i).get_dat(j,j);
          total+=ss*ss*var;
	};
      };
    };
  };
  return total; // variance  
};

void UncertaintyCalculationForYieldDecay::CalDecayEnergyUncertaintyDetail(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,string core,string chara,int step,real mean){

  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Decay energy-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  real total=0;
  for(int i=0;i<3;i++){
    int size=cov.GetSize();
    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    for(int j=0;j<size;j++){
      int mat=cov.GetID(j);
      int mt=99990+i;
      int data=sens_data.FindData0D(mat,mt);
      real tmp;
      if(data==-1){
	tmp=0.;
      }else{
	real tmp1=sens_data.GetSensitivity0D(data);
	tmp=tmp1;
      };
      real sd=sqrt(M.get_dat(j,j));
      real vali=tmp*tmp*M.get_dat(j,j);
      string type;
      if(i==0)type=" Beta";
      if(i==1)type="Gamma";
      if(i==2)type="Alpha";      
      if(sqrt(vali)>=mean){
	cout<<"* "<<midt.GetName(mat)<<" "<<type<<" decay energy (Error: "<<sd<<"): "<<vali<<" ( "<<sqrt(vali)<<" )\n";
      };
      total+=vali;
    };
};
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";

};

void UncertaintyCalculationForYieldDecay::CalHalfLifeUncertainty(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step){

  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Half life-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";
  
  cout.setf(ios::scientific);
  cout.precision(4);
  
  int size=cov.GetSize();
  GroupData1D G;
  G.put_imax(size);
  
  for(int j=0;j<size;j++){
    real sns_total;
    int mat=cov.GetID(j);
    int data=sens_data.FindData0D(mat,8888);
    if(data==-1){
      sns_total=0.;
    }else{
      real tmp1=sens_data.GetSensitivity0D(data);
      sns_total=tmp1;
    };
    G.put_data(j,sns_total);
  };
  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();
  real unc=G*(M*G);

  cout<<"* Total uncertainty from half life "<<unc<<" ( "<<sqrt(unc)<<" )\n";
  cout<<"********************************************************\n";

};

real UncertaintyCalculationForYieldDecay::CalHalfLifeUncertaintyNuclideWise(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step,string nucname)
{
  MATIDTranslator midt;
  int id=midt.ID(nucname);

  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  int size=cov.GetSize();
  for(int j=0;j<size;j++){
    real sns_total;
    int mat=cov.GetID(j);
    if(mat==id){
      int data=sens_data.FindData0D(mat,8888);
      if(data==-1){
        return 0.;
      }else{
        real tmp1=sens_data.GetSensitivity0D(data);
        return tmp1*tmp1*cov.GetCovarianceMatrix().get_dat(j,j);
      };
    };
  };

  return 0.;
};

void UncertaintyCalculationForYieldDecay::CalHalfLifeUncertaintyDetail(SensitivityContainer &sens_con,HalfLifeCovariance &cov,string core,string chara,int step,real mean){

  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Half life-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";
  
  cout.setf(ios::scientific);
  cout.precision(4);
  
  int size=cov.GetSize();
  real total=0.;
  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();
  for(int j=0;j<size;j++){
    real sns_total;
    int mat=cov.GetID(j);
    int data=sens_data.FindData0D(mat,8888);
    if(data==-1){
      sns_total=0.;
    }else{
      real tmp1=sens_data.GetSensitivity0D(data);
      sns_total=tmp1;
    };
    real sd=sqrt(M.get_dat(j,j));
    real vali=sns_total*sns_total*M.get_dat(j,j);
    if(sqrt(vali)>=mean){
      cout<<"* "<<midt.GetName(mat)<<" Half Life (Error: "<<sd<<"): "<<vali<<" ( "<<sqrt(vali)<<" ) \n";
    };
    total+=vali;
  };

  cout<<"********************************************************\n";
  cout<<"* Total uncertainty from half life "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";

};

void UncertaintyCalculationForYieldDecay::CalBranchingRatioUncertainty(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,string core,string chara,int step){
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Branching ratio-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";
  
  cout.setf(ios::scientific);
  cout.precision(4);
  
  int size=cov.GetSize();
  int nucnum=cov.GetNucNum();
  GroupData1D G;
  G.put_imax(size);
  
  int counter=0;
  for(int j=0;j<nucnum;j++){
    int mat=cov.GetID(j);    
    int channel=cov.GetChannel(j);
    for(int i=0;i<channel;i++){
      int mt=88880+i;
      int data=sens_data.FindData0D(mat,mt);
      real sns;
      if(data==-1){
	sns=0.;
      }else{
	real tmp=sens_data.GetSensitivity0D(data);
	sns=tmp;
      };
      G.put_data(counter,sns);
      counter++;
    };
  };
  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();
  real unc=G*(M*G);
  
  cout<<"* Total uncertainty from branching ratio "<<unc<<" ( "<<sqrt(unc)<<" )\n";
  cout<<"********************************************************\n";

};

void UncertaintyCalculationForYieldDecay::CalBranchingRatioUncertaintyDetail(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,string core,string chara,int step,real mean,bool decay)
{
  int mtid_base=88880;
  if(!decay)mtid_base=1020;

  string type="decay";
  if(!decay)type="reaction";

  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  
  cout<<"********************************************************\n";
  cout<<"* Branching ratio-induced uncertainty ("<<type<<")\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";
  
  cout.setf(ios::scientific);
  cout.precision(4);
  
  int nucnum=cov.GetNucNum();
  int counter=0;
  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();

  int size=cov.GetSize();
  GroupData1D G;
  G.put_imax(size);

  real total=0.;
  for(int j=0;j<nucnum;j++){
    int mat=cov.GetID(j);
    int channel=cov.GetChannel(j);
    real vali=0.;
    for(int i=0;i<channel;i++){
      int mt=mtid_base+i;
      int data=sens_data.FindData0D(mat,mt);
      real sns=0.;
      if(data!=-1)sns=sens_data.GetSensitivity0D(data);
      for(int j=0;j<channel;j++){
	int mt2=mtid_base+j;
	int data2=sens_data.FindData0D(mat,mt2);
	real sns2=0.;
	if(data2!=-1)sns2=sens_data.GetSensitivity0D(data2);
	vali+=sns*sns2*M.get_dat(counter+i,counter+j);
      };
    };
    if(sqrt(vali)>=mean){
      cout<<"* "<<midt.GetName(mat)<<" Branching Ratio : "<<vali<<" ( "<<sqrt(vali)<<" ) ";
      for(int j=0;j<channel;j++){
	int mt2=mtid_base+j;
	int data2=sens_data.FindData0D(mat,mt2);
	real sns2=0.;
	if(data2!=-1)sns2=sens_data.GetSensitivity0D(data2);
	cout<<sns2<<" ";
      };
      cout<<"\n";
    };
    counter+=channel;
    total+=vali;
  };
  
  /*
  // Kawamoto-original
  for(int j=0;j<nucnum;j++){
    int mat=cov.GetID(j);
    int channel=cov.GetChannel(j);
    for(int i=0;i<channel;i++){
      int mt=88880+i;
      int data=sens_data.FindData0D(mat,mt);
      real sns;
      if(data==-1){
	sns=0.;
      }else{
	real tmp=sens_data.GetSensitivity0D(data);
	sns=tmp;
      };
      G.put_data(counter,sns);
      real sd=sqrt(M.get_dat(counter,counter));
      real vali=sns*sns*M.get_dat(counter,counter);
      if(sqrt(vali)>=mean){
	cout<<"* "<<midt.GetName(mat)<<" Branching Ratio No."<<i<<" (Error: "<<sd<<"): "<<vali<<" ( "<<sqrt(vali)<<" ) \n";
      };
      total+=vali;
      counter++;
    };
  };
  */

  //total=G*(M*G);
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty from branching ratio "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";

};

GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldGMG(SensitivityContainer &sens_con,IndependentYieldCovariance &cov,int chara_num,string *core,string *chara,int *step){
  MATIDTranslator midt;

  GroupData2D gmg_total;
  gmg_total.put_yx(chara_num,chara_num);
  gmg_total.set_zero();

  int fisnucnum=cov.GetFissionNuclideNumber();

  //cout<<fisnucnum<<" "<<chara_num<<"\n";

  for(int i=0;i<fisnucnum;i++){
    int size=cov.GetSize(i);
    //cout<<i<<" "<<size<<"...\n";
    string fisnucname=cov.GetFissionNuclideName(i);
    GroupData2D G;
    G.put_yx(chara_num,size);
    for(int k=0;k<chara_num;k++){
      SensitivityData sens_data;
      sens_data=sens_con.GetSensitivityData(core[k],chara[k],step[k]);
      for(int j=0;j<size;j++){
	int mat=cov.GetID(i,j);
	int mt=18000000+midt.ID(fisnucname);
	int data=sens_data.FindData0D(mat,mt);
	//cout<<i<<" "<<k<<" "<<j<<" "<<mat<<" "<<mt<<" "<<data<<"\n";
	if(data==-1){
	  G.put_data(k,j,0.);
	}else{
	  real tmp1=sens_data.GetSensitivity0D(data);
	  real tmp2=tmp1;
	  //cout<<"     "<<k<<" "<<chara_num<<" "<<j<<" "<<size<<" "<<tmp2<<"\n";
	  G.put_data(k,j,tmp2);
	};
      };
    };

    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    GroupData2D gmg=G*(M*G.T());

    gmg_total=gmg_total+gmg;
  };

  return gmg_total;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldGMG(SensitivityContainer &sens_con,IndependentYieldCovariance &cov){

  int chara_num=sens_con.GetSize();
  string * core = new string[chara_num];
  string * chara = new string[chara_num];
  int * step = new int[chara_num];
  
  for(int i=0;i<chara_num;i++){
    core[i]=sens_con.GetParameterList()->GetCoreTag(i);
    chara[i]=sens_con.GetParameterList()->GetCharaTag(i);
    step[i]=sens_con.GetParameterList()->GetStepTag(i);
  };
  GroupData2D gmg=GetFissionYieldGMG(sens_con,cov,chara_num,core,chara,step);

  delete[] core;
  delete[] chara;
  delete[] step;

  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldGMG(SensitivityContainer &sens_con1,IndependentYieldCovariance &cov,SensitivityContainer &sens_con2){
  MATIDTranslator midt;

  GroupData2D gmg_total;

  int chara_num1=sens_con1.GetSize();
  int chara_num2=sens_con2.GetSize();
  gmg_total.put_yx(chara_num1,chara_num2);
  gmg_total.set_zero();
  vector<string> core1;
  vector<string> chara1;
  vector<int> step1;
  for(int i=0;i<chara_num1;i++){
    core1.push_back(sens_con1.GetParameterList()->GetCoreTag(i));
    chara1.push_back(sens_con1.GetParameterList()->GetCharaTag(i));
    step1.push_back(sens_con1.GetParameterList()->GetStepTag(i));
  };
  vector<string> core2;
  vector<string> chara2;
  vector<int> step2;
  for(int i=0;i<chara_num2;i++){
    core2.push_back(sens_con2.GetParameterList()->GetCoreTag(i));
    chara2.push_back(sens_con2.GetParameterList()->GetCharaTag(i));
    step2.push_back(sens_con2.GetParameterList()->GetStepTag(i));
  };
  
  int fisnucnum=cov.GetFissionNuclideNumber();

  for(int i=0;i<fisnucnum;i++){
    int size=cov.GetSize(i);
    string fisnucname=cov.GetFissionNuclideName(i);
    GroupData2D G1;
    GroupData2D G2;
    G1.put_yx(chara_num1,size);
    G2.put_yx(chara_num2,size);

    for(int k=0;k<chara_num1;k++){
      SensitivityData sens_data1;
      sens_data1=sens_con1.GetSensitivityData(core1[k],chara1[k],step1[k]);
      for(int j=0;j<size;j++){
	int mat1=cov.GetID(i,j);
	int mt1=18000000+midt.ID(fisnucname);
	int data1=sens_data1.FindData0D(mat1,mt1);
	if(data1==-1){
	  G1.put_data(k,j,0.);
	}else{
	  real tmp1=sens_data1.GetSensitivity0D(data1);
	  real tmp2=tmp1;
	  G1.put_data(k,j,tmp2);
	};
      };
    };
    for(int k=0;k<chara_num2;k++){
      SensitivityData sens_data2;
      sens_data2=sens_con2.GetSensitivityData(core2[k],chara2[k],step2[k]);
      for(int j=0;j<size;j++){
	int mat2=cov.GetID(i,j);
	int mt2=18000000+midt.ID(fisnucname);
	int data2=sens_data2.FindData0D(mat2,mt2);
	if(data2==-1){
	  G2.put_data(k,j,0.);
	}else{
	  real tmp3=sens_data2.GetSensitivity0D(data2);
	  real tmp4=tmp3;
	  G2.put_data(k,j,tmp4);
	};
      };
    };

    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    GroupData2D gmg=G1*(M*G2.T());

    gmg_total=gmg_total+gmg;
  };

  return gmg_total;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetDecayEnergyGMG(SensitivityContainer &sens_con,DecayEnergyCovariance &cov,int chara_num,string *core,string *chara,int *step){
  MATIDTranslator midt;
 
  GroupData2D gmg_total;
  gmg_total.put_yx(chara_num,chara_num);
  gmg_total.set_zero();

  for(int i=0;i<3;i++){
    int size=cov.GetSize();
    GroupData2D G;
    G.put_yx(chara_num,size);
    for(int k=0;k<chara_num;k++){
      SensitivityData sens_data;
      sens_data=sens_con.GetSensitivityData(core[k],chara[k],step[k]);
      for(int j=0;j<size;j++){
	int mat=cov.GetID(j);
	int mt=99990+i;
	int data=sens_data.FindData0D(mat,mt);
	if(data==-1){
	  G.put_data(k,j,0.);
	}else{
	  real tmp1=sens_data.GetSensitivity0D(data);
	  real tmp2=tmp1;
	  G.put_data(k,j,tmp2);
	  /*
	  real ttt=sens_data.GetSensitivity0D(data);
	  if(ttt>0.){
  	    cout<<i<<" "<<k<<" "<<j<<" "<<mat<<" "<<mt<<" "<<ttt<<" "<<cov.GetCovarianceMatrix(i).get_dat(j,j)<<"\n";
	  };
	  */
	};
      };
    };
    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    GroupData2D gmg=G*(M*G.T());
    gmg_total=gmg_total+gmg;
  };

  return gmg_total;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetDecayEnergyGMG(SensitivityContainer &sens_con,DecayEnergyCovariance &cov){

  int chara_num=sens_con.GetSize();
  string * core = new string[chara_num];
  string * chara = new string[chara_num];
  int * step = new int[chara_num];
  
  for(int i=0;i<chara_num;i++){
    core[i]=sens_con.GetParameterList()->GetCoreTag(i);
    chara[i]=sens_con.GetParameterList()->GetCharaTag(i);
    step[i]=sens_con.GetParameterList()->GetStepTag(i);
  };
  GroupData2D gmg=GetDecayEnergyGMG(sens_con,cov,chara_num,core,chara,step);

  delete[] core;
  delete[] chara;
  delete[] step;

  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetDecayEnergyGMG(SensitivityContainer &sens_con1,DecayEnergyCovariance &cov,SensitivityContainer &sens_con2){
  MATIDTranslator midt;

  GroupData2D gmg_total;

  int chara_num1=sens_con1.GetSize();
  int chara_num2=sens_con2.GetSize();
  gmg_total.put_yx(chara_num1,chara_num2);
  gmg_total.set_zero();
  vector<string> core1;
  vector<string> chara1;
  vector<int> step1;
  for(int i=0;i<chara_num1;i++){
    core1.push_back(sens_con1.GetParameterList()->GetCoreTag(i));
    chara1.push_back(sens_con1.GetParameterList()->GetCharaTag(i));
    step1.push_back(sens_con1.GetParameterList()->GetStepTag(i));
  };
  vector<string> core2;
  vector<string> chara2;
  vector<int> step2;
  for(int i=0;i<chara_num2;i++){
    core2.push_back(sens_con2.GetParameterList()->GetCoreTag(i));
    chara2.push_back(sens_con2.GetParameterList()->GetCharaTag(i));
    step2.push_back(sens_con2.GetParameterList()->GetStepTag(i));
  };

  for(int i=0;i<3;i++){
    int size=cov.GetSize();
    GroupData2D G1;
    GroupData2D G2;
    G1.put_yx(chara_num1,size);
    G2.put_yx(chara_num2,size);

    for(int k=0;k<chara_num1;k++){
      SensitivityData sens_data1;
      sens_data1=sens_con1.GetSensitivityData(core1[k],chara1[k],step1[k]);
      for(int j=0;j<size;j++){
	int mat1=cov.GetID(j);
	int mt1=99990+i;
	int data1=sens_data1.FindData0D(mat1,mt1);
	if(data1==-1){
	  G1.put_data(k,j,0.);
	}else{
	  real tmp1=sens_data1.GetSensitivity0D(data1);
	  real tmp2=tmp1;
	  G1.put_data(k,j,tmp2);
	};
      };
    };
    for(int k=0;k<chara_num2;k++){
      SensitivityData sens_data2;
      sens_data2=sens_con2.GetSensitivityData(core2[k],chara2[k],step2[k]);
      for(int j=0;j<size;j++){
	int mat2=cov.GetID(j);
	int mt2=99990+i;
	int data2=sens_data2.FindData0D(mat2,mt2);
	if(data2==-1){
	  G2.put_data(k,j,0.);
	}else{
	  real tmp3=sens_data2.GetSensitivity0D(data2);
	  real tmp4=tmp3;
	  G2.put_data(k,j,tmp4);
	};
      };
    };

    GroupData2D M=cov.GetCovarianceMatrix(i);
    M.ReducedForm();
    GroupData2D gmg=G1*(M*G2.T());
    gmg_total=gmg_total+gmg;
  };

  return gmg_total;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetHalfLifeGMG(SensitivityContainer &sens_con,HalfLifeCovariance &cov,int chara_num,string *core,string *chara,int *step){
  MATIDTranslator midt;
 
  int size=cov.GetSize();
  GroupData2D G;
  G.put_yx(chara_num,size);
  GroupData2D M=cov.GetCovarianceMatrix();

  for(int k=0;k<chara_num;k++){
    SensitivityData sens_data;
    sens_data=sens_con.GetSensitivityData(core[k],chara[k],step[k]);
    for(int j=0;j<size;j++){
      real sns_total;
      int mat=cov.GetID(j);
      int data=sens_data.FindData0D(mat,8888);
      if(data==-1){
	sns_total=0.;
      }else{
	real tmp1=sens_data.GetSensitivity0D(data);
	sns_total=tmp1;
	/*
	real tmp2=M.get_dat(j,j);
	if(fabs(tmp1*tmp1*tmp2)>1e-6){
	  tmp2=sqrt(tmp2);
	  cout<<k<<" "<<j<<" "<<mat<<" "<<tmp1<<" "<<tmp2<<" "<<tmp1*tmp2<<"\n";
	};
	*/
      };
      G.put_data(k,j,sns_total);
    };
  };

  M.ReducedForm();
  GroupData2D gmg=G*(M*G.T());

  return gmg;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetHalfLifeGMG(SensitivityContainer &sens_con,HalfLifeCovariance &cov){

  int chara_num=sens_con.GetSize();
  string * core = new string[chara_num];
  string * chara = new string[chara_num];
  int * step = new int[chara_num];
  
  for(int i=0;i<chara_num;i++){
    core[i]=sens_con.GetParameterList()->GetCoreTag(i);
    chara[i]=sens_con.GetParameterList()->GetCharaTag(i);
    step[i]=sens_con.GetParameterList()->GetStepTag(i);
  };
  GroupData2D gmg=GetHalfLifeGMG(sens_con,cov,chara_num,core,chara,step);

  delete[] core;
  delete[] chara;
  delete[] step;

  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetHalfLifeGMG(SensitivityContainer &sens_con1,HalfLifeCovariance &cov,SensitivityContainer &sens_con2){
  MATIDTranslator midt;

  int chara_num1=sens_con1.GetSize();
  int chara_num2=sens_con2.GetSize();
  vector<string> core1;
  vector<string> chara1;
  vector<int> step1;
  for(int i=0;i<chara_num1;i++){
    core1.push_back(sens_con1.GetParameterList()->GetCoreTag(i));
    chara1.push_back(sens_con1.GetParameterList()->GetCharaTag(i));
    step1.push_back(sens_con1.GetParameterList()->GetStepTag(i));
  };
  vector<string> core2;
  vector<string> chara2;
  vector<int> step2;
  for(int i=0;i<chara_num2;i++){
    core2.push_back(sens_con2.GetParameterList()->GetCoreTag(i));
    chara2.push_back(sens_con2.GetParameterList()->GetCharaTag(i));
    step2.push_back(sens_con2.GetParameterList()->GetStepTag(i));
  };
 
  int size=cov.GetSize();
  GroupData2D G1;
  GroupData2D G2;
  G1.put_yx(chara_num1,size);
  G2.put_yx(chara_num2,size);

  for(int k=0;k<chara_num1;k++){
    SensitivityData sens_data1;
    sens_data1=sens_con1.GetSensitivityData(core1[k],chara1[k],step1[k]);
    for(int j=0;j<size;j++){
      real sns_total1;
      int mat1=cov.GetID(j);
      int data1=sens_data1.FindData0D(mat1,8888);
      if(data1==-1){
	sns_total1=0.;
      }else{
	real tmp1=sens_data1.GetSensitivity0D(data1);
	sns_total1=tmp1;
      };
      G1.put_data(k,j,sns_total1);
    };
  };
  for(int k=0;k<chara_num2;k++){
    SensitivityData sens_data2;
    sens_data2=sens_con2.GetSensitivityData(core2[k],chara2[k],step2[k]);
    for(int j=0;j<size;j++){
      real sns_total2;
      int mat2=cov.GetID(j);
      int data2=sens_data2.FindData0D(mat2,8888);
      if(data2==-1){
	sns_total2=0.;
      }else{
	real tmp2=sens_data2.GetSensitivity0D(data2);
	sns_total2=tmp2;
      };
      G2.put_data(k,j,sns_total2);
    };
  };

  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();
  GroupData2D gmg=G1*(M*G2.T());

  return gmg;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetBranchingRatioGMG(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,int chara_num,string *core,string *chara,int *step,bool decay)
{
  int mtid_base=88880;
  if(!decay)mtid_base=1020;

  MATIDTranslator midt;
 
  int size=cov.GetSize();
  int nucnum=cov.GetNucNum();
  GroupData2D G;
  G.put_yx(chara_num,size);
  
  for(int k=0;k<chara_num;k++){
    SensitivityData sens_data;
    sens_data=sens_con.GetSensitivityData(core[k],chara[k],step[k]);
    int counter=0;
    for(int j=0;j<nucnum;j++){
      int mat=cov.GetID(j);
      int channel=cov.GetChannel(j);
      for(int i=0;i<channel;i++){
	//int mt=88880+i;
	int mt=mtid_base+i;
	int data=sens_data.FindData0D(mat,mt);
	real sns;
	if(data==-1){
	  sns=0.;
	}else{
	  real tmp=sens_data.GetSensitivity0D(data);
	  sns=tmp;
	};
	G.put_data(k,counter,sns);
	counter++;
      };
    };
  };
  GroupData2D M=cov.GetCovarianceMatrix();

  M.ReducedForm();
  GroupData2D gmg=G*(M*G.T());

  //gmg.show_self();

  return gmg;  
};

GroupData2D UncertaintyCalculationForYieldDecay::GetBranchingRatioGMG(SensitivityContainer &sens_con,BranchingRatioCovariance &cov,bool decay){

  int chara_num=sens_con.GetSize();
  string * core = new string[chara_num];
  string * chara = new string[chara_num];
  int * step = new int[chara_num];
  
  for(int i=0;i<chara_num;i++){
    core[i]=sens_con.GetParameterList()->GetCoreTag(i);
    chara[i]=sens_con.GetParameterList()->GetCharaTag(i);
    step[i]=sens_con.GetParameterList()->GetStepTag(i);
  };
  GroupData2D gmg=GetBranchingRatioGMG(sens_con,cov,chara_num,core,chara,step,decay);

  delete[] core;
  delete[] chara;
  delete[] step;

  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetBranchingRatioGMG(SensitivityContainer &sens_con1,BranchingRatioCovariance &cov,SensitivityContainer &sens_con2,bool decay)
{
  int mtid_base=88880;
  if(!decay)mtid_base=1020;

  MATIDTranslator midt;
 
  int chara_num1=sens_con1.GetSize();
  int chara_num2=sens_con2.GetSize();
  vector<string> core1;
  vector<string> chara1;
  vector<int> step1;
  for(int i=0;i<chara_num1;i++){
    core1.push_back(sens_con1.GetParameterList()->GetCoreTag(i));
    chara1.push_back(sens_con1.GetParameterList()->GetCharaTag(i));
    step1.push_back(sens_con1.GetParameterList()->GetStepTag(i));
  };
  vector<string> core2;
  vector<string> chara2;
  vector<int> step2;
  for(int i=0;i<chara_num2;i++){
    core2.push_back(sens_con2.GetParameterList()->GetCoreTag(i));
    chara2.push_back(sens_con2.GetParameterList()->GetCharaTag(i));
    step2.push_back(sens_con2.GetParameterList()->GetStepTag(i));
  };

  int size=cov.GetSize();
  int nucnum=cov.GetNucNum();
  GroupData2D G1;
  GroupData2D G2;
  G1.put_yx(chara_num1,size);
  G2.put_yx(chara_num2,size);
  
  for(int k=0;k<chara_num1;k++){
    SensitivityData sens_data1;
    sens_data1=sens_con1.GetSensitivityData(core1[k],chara1[k],step1[k]);
    int counter1=0;
    for(int j=0;j<nucnum;j++){
      int mat1=cov.GetID(j);
      int channel1=cov.GetChannel(j);
      for(int i=0;i<channel1;i++){
	//int mt1=88880+i;
	int mt1=mtid_base+i;
	int data1=sens_data1.FindData0D(mat1,mt1);
	real sns1;
	if(data1==-1){
	  sns1=0.;
	}else{
	  real tmp1=sens_data1.GetSensitivity0D(data1);
	  sns1=tmp1;
	};
	G1.put_data(k,counter1,sns1);
	counter1++;
      };
    };
  };
  for(int k=0;k<chara_num2;k++){
    SensitivityData sens_data2;
    sens_data2=sens_con2.GetSensitivityData(core2[k],chara2[k],step2[k]);
    int counter2=0;
    for(int j=0;j<nucnum;j++){
      int mat2=cov.GetID(j);
      int channel2=cov.GetChannel(j);
      for(int i=0;i<channel2;i++){
	//int mt2=88880+i;
	int mt2=mtid_base+i;
	int data2=sens_data2.FindData0D(mat2,mt2);
	real sns2;
	if(data2==-1){
	  sns2=0.;
	}else{
	  real tmp2=sens_data2.GetSensitivity0D(data2);
	  sns2=tmp2;
	};
	G2.put_data(k,counter2,sns2);
	counter2++;
      };
    };
  };

  GroupData2D M=cov.GetCovarianceMatrix();
  M.ReducedForm();
  GroupData2D gmg=G1*(M*G2.T());

  return gmg;  
};
GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldCorrGMG(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &fbcov,int chara_num,string *core,string *chara,int *step){
  MATIDTranslator midt;

  int fisnucnum=fbcov.GetFissionNuclideNumber();
  int size=fbcov.GetSize();
  vector<string> fisnucname(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    fisnucname[i]=fbcov.GetFissionNuclideName(i);
  };
  GroupData2D G;
  G.put_yx(chara_num,size);
  for(int k=0;k<chara_num;k++){
    SensitivityData sens_data;
    sens_data=sens_con.GetSensitivityData(core[k],chara[k],step[k]);
    for(int j=0;j<size;j++){
      int mat=fbcov.GetID(j);
      int mt=18000000+midt.ID(fisnucname[floor((j)/(size/fisnucnum))]);
      int data=sens_data.FindData0D(mat,mt);
      if(data==-1){
	G.put_data(k,j,0.);
      }else{
	real tmp=sens_data.GetSensitivity0D(data);
	G.put_data(k,j,tmp);
      };
    };
  };
  //G.show_self();
  GroupData2D M=fbcov.GetCovarianceMatrix();
  M.ReducedForm();  // ?
  GroupData2D gmg=G*(M*G.T());
  
  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldCorrGMG(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &fbcov){ // okumura

  int chara_num=sens_con.GetSize();
  string * core = new string[chara_num];
  string * chara = new string[chara_num];
  int * step = new int[chara_num];

  for(int i=0;i<chara_num;i++){
    core[i]=sens_con.GetParameterList()->GetCoreTag(i);
    chara[i]=sens_con.GetParameterList()->GetCharaTag(i);
    step[i]=sens_con.GetParameterList()->GetStepTag(i);
  };
  GroupData2D gmg=GetFissionYieldCorrGMG(sens_con,fbcov,chara_num,core,chara,step);

  delete[] core;
  delete[] chara;
  delete[] step;

  return gmg;
};

GroupData2D UncertaintyCalculationForYieldDecay::GetFissionYieldCorrGMG(SensitivityContainer &sens_con1,FissionYieldCovarianceCorr &fbcov,SensitivityContainer &sens_con2){
  MATIDTranslator midt;

  int chara_num1=sens_con1.GetSize();
  int chara_num2=sens_con2.GetSize();
  vector<string> core1;
  vector<string> chara1;
  vector<int> step1;
  for(int i=0;i<chara_num1;i++){
    core1.push_back(sens_con1.GetParameterList()->GetCoreTag(i));
    chara1.push_back(sens_con1.GetParameterList()->GetCharaTag(i));
    step1.push_back(sens_con1.GetParameterList()->GetStepTag(i));
  };
  vector<string> core2;
  vector<string> chara2;
  vector<int> step2;
  for(int i=0;i<chara_num2;i++){
    core2.push_back(sens_con2.GetParameterList()->GetCoreTag(i));
    chara2.push_back(sens_con2.GetParameterList()->GetCharaTag(i));
    step2.push_back(sens_con2.GetParameterList()->GetStepTag(i));
  };

  int fisnucnum=fbcov.GetFissionNuclideNumber();
  GroupData2D G1;
  GroupData2D G2;
  
  int size=fbcov.GetSize();  // fisnucnum x fpnucnum
  G1.put_yx(chara_num1,size);
  G2.put_yx(chara_num2,size);
  vector<string> fisnucname(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    fisnucname[i]=fbcov.GetFissionNuclideName(i);
  };
  
  for(int k=0;k<chara_num1;k++){
    SensitivityData sens_data1;
    sens_data1=sens_con1.GetSensitivityData(core1[k],chara1[k],step1[k]);
    for(int j=0;j<size;j++){
      int mat1=fbcov.GetID(j);
      int mt1=18000000+midt.ID(fisnucname[floor((j)/(size/fisnucnum))]);
      int data1=sens_data1.FindData0D(mat1,mt1);
      if(data1==-1){
	G1.put_data(k,j,0.);
      }else{
	real tmp1=sens_data1.GetSensitivity0D(data1);
	G1.put_data(k,j,tmp1);
      };
    };
  };

  for(int k=0;k<chara_num2;k++){
    SensitivityData sens_data2;
    sens_data2=sens_con2.GetSensitivityData(core2[k],chara2[k],step2[k]);
    for(int j=0;j<size;j++){
      int mat2=fbcov.GetID(j);
      int mt2=18000000+midt.ID(fisnucname[floor((j)/(size/fisnucnum))]);
      int data2=sens_data2.FindData0D(mat2,mt2);
      if(data2==-1){
	G2.put_data(k,j,0);
      }else{
	real tmp2=sens_data2.GetSensitivity0D(data2);
	G2.put_data(k,j,tmp2);
      };
    };
  };

  GroupData2D M=fbcov.GetCovarianceMatrix();
  //M.ReducedForm();
  GroupData2D gmg=G1*(M*G2.T());
  
  return gmg;
};

real UncertaintyCalculationForYieldDecay::CalFissionYieldCorrUncertainty(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &cov,string core,string chara,int step){ // okumura
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  int fisnucnum=cov.GetFissionNuclideNumber();
  
  cout<<"********************************************************\n";
  cout<<"* Fission yield-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  real total=0;
  vector<string> name_store;
  int size=cov.GetSize();
  vector<string> fisnucname(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    fisnucname[i]=cov.GetFissionNuclideName(i);
  };

  GroupData2D M=cov.GetCovarianceMatrix();

  //cout<<M.get_dat(734,739)<<"\n";

  GroupData1D G;
  G.put_imax(size);
  for(int j=0;j<size;j++){
    int mat=cov.GetID(j);
    int mt=18000000+midt.ID(fisnucname[floor((j)/(size/fisnucnum))]);
    int data=sens_data.FindData0D(mat,mt);
    if(data==-1){
      G.put_data(j,0.);
    }else{
      real tmp1=sens_data.GetSensitivity0D(data);
      real tmp2=tmp1;
      G.put_data(j,tmp2);

      //real tmp=tmp2*tmp2*M.get_dat(j,j);
      //if(fabs(tmp)>1e-10)cout<<j<<" "<<tmp<<" "<<mat<<" "<<tmp2<<" "<<M.get_dat(j,j)<<"\n";
      //if(j<800&&j>700)cout<<j<<" "<<tmp<<" "<<mat<<" "<<tmp2<<" "<<M.get_dat(j,j)<<"\n";
      //if(j==734||j==739)cout<<j<<" "<<tmp<<" "<<mat<<" "<<tmp2<<" "<<M.get_dat(j,j)<<"\n";                  

    };
  };
  

  M.ReducedForm();
  real unc=G*(M*G);
  
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<unc<<" ( "<<sqrt(unc)<<" )\n";
  cout<<"********************************************************\n";
  return unc; // okumura
};


void UncertaintyCalculationForYieldDecay::CalFissionYieldCorrUncertaintyDetail(SensitivityContainer &sens_con,FissionYieldCovarianceCorr &cov,string core,string chara, int step,real mean)
{ 
  MATIDTranslator midt;
  SensitivityData sens_data;
  sens_data=sens_con.GetSensitivityData(core,chara,step);
  int fisnucnum=cov.GetFissionNuclideNumber();

  cout<<"********************************************************\n";
  cout<<"* Fission yield-induced uncertainty\n";
  cout<<"********************************************************\n";
  cout<<"* Integral Data : core    "<<core<<"\n";
  cout<<"*                 Chara.  "<<chara<<"\n";
  cout<<"*                 step    "<<step<<"\n";
  cout<<"*                 value   "<<sens_data.GetValue()<<"\n";
  cout<<"********************************************************\n";

  cout.setf(ios::scientific);
  cout.precision(4);

  int size=cov.GetSize();
  int fpnum=size/fisnucnum;
  GroupData2D M=cov.GetCovarianceMatrix();

  for(int i=0;i<fisnucnum;i++){

    string fisnucname=cov.GetFissionNuclideName(i);
    int mt=18000000+midt.ID(fisnucname);
    GroupData1D mini_G(fpnum);
    GroupData2D mini_cov(fpnum,fpnum);

    for(int j=0;j<fpnum;j++){
      int mat=cov.GetID(j);
      int data=sens_data.FindData0D(mat,mt);
      real sss=0.;
      if(data!=-1)sss=sens_data.GetSensitivity0D(data);
      //if(mat!=521331&&mat!=531330)sss=0.;
      mini_G.put_data(j,sss);
      for(int k=0;k<fpnum;k++){
	mini_cov.put_data(j,k,M.get_dat(i*fpnum+j,i*fpnum+k));
      };
    };

    real vali=mini_G*(mini_cov*mini_G);
    if(sqrt(vali)>=mean){
      cout<<"* "<<fisnucname<<" : "<<vali<<" ( "<<sqrt(vali)<<" ) \n";
    };
  };
  /*
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";
  */


  /*
  // !! caution
  // This method does NOT work because correlation among different mass chains should exit
  // in reduced chain


 
  real total=0;
  vector<string> name_store;

    
  for(int i=0;i<fisnucnum;i++){
    string fisnucname=cov.GetFissionNuclideName(i);
    for(int ii=1;ii<300;ii++){
      int mass_counter=0;
      vector<real> mini_G_tmp;
      vector<int> mini_id;
      vector<int> mini_position;
      for(int j=0;j<fpnum;j++){
	int mat=cov.GetID(j);
	int lz,la,li;
	midt.GetParameter(mat,lz,la,li);
	if(ii==la){
	  mass_counter++;
	  mini_id.push_back(mat);
	  mini_position.push_back(j);
	  int mt=18000000+midt.ID(fisnucname);
	  int data=sens_data.FindData0D(mat,mt);
	  if(data==-1){
	    mini_G_tmp.push_back(0.);
	  }else{
	    real tmp1=sens_data.GetSensitivity0D(data);
	    real tmp2=tmp1;
	    mini_G_tmp.push_back(tmp2);
	  };
	};
      };
      
      if(mass_counter!=0){
	GroupData1D mini_G;
	GroupData2D mini_cov;
	mini_cov.put_yx(mass_counter,mass_counter);
	mini_G.put_imax(mass_counter);
	for(int iii=0;iii<mass_counter;iii++){
	  mini_G.put_data(iii,mini_G_tmp[iii]);
	  for(int iiii=0;iiii<mass_counter;iiii++){
	    mini_cov.put_data(iii,iiii,M.get_dat(mini_position[iii],mini_position[iiii]));
	  };
	};
	real vali=mini_G*(mini_cov*mini_G);
	if(sqrt(vali)>=mean){
	  cout<<"* "<<fisnucname<<" A="<<ii<<" Mass Chain: "<<vali<<" ( "<<sqrt(vali)<<" ) \n";
	  //cout<<"   "<<ii<<"  "<<sqrt(vali)<<"\n";
	};
	total=total+vali;
      };	
    };
  };
  cout<<"********************************************************\n";
  cout<<"* Total uncertainty "<<total<<" ( "<<sqrt(total)<<" )\n";
  cout<<"********************************************************\n";
  */
  
};

