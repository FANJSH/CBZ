#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "SUTool.h"

using namespace std;

int SUTool::GetSensitivityIndex(string name)
{
  int sz=sns_data.size();
  for(int i=0;i<sz;i++){
    if(name==sns_name[i])return i;
  };
  return -1;
};

void SUTool::ReadSensitivityDataFromFile(string mdir,string filename,string data_name)
{
  SensitivityData snstmp;

  sns_data.push_back(snstmp);
  sns_name.push_back(data_name);

  int sz=sns_data.size();
  sns_data[sz-1].ReadFile(mdir,filename);
};

SensitivityData &SUTool::GetSensitivityData(string name)
{
  int id=GetSensitivityIndex(name);
  if(id==-1){
    cout<<"# Error in SUTool::GetSensitivityData.\n";
    exit(0);
  };
  return sns_data[id];
};

enum xstype SUTool::GetXSType(int mt)
{
  enum xstype xst;
  if(mt==18){
    xst=sigf;
  }else if(mt==452){
    xst=nu;
  }else if(mt==181){
    xst=chi;
  }else if(mt==251){
    xst=mu;
  }else if(mt==102){
    xst=sigc;
  }else if(mt==2){
    xst=sigel;
  }else if(mt==4){
    xst=siginel;
  }else if(mt==16){
    xst=sign2n;
  }else{
    cout<<"# Error in SUTool:GetXSType.\n";
    cout<<"# Not coded for MT="<<mt<<"\n";
    exit(0);
  };
  return xst;
};

void SUTool::CalLibraryEffect(string data_name,XSLibrary &lib1,XSLibrary &lib2,bool grp_wise_prt)
{
  real threshold=0.0001;
  
  int id=GetSensitivityIndex(data_name);
  if(id==-1){
    cout<<"# Error in SUTool::CalLibraryEffect.\n";
    cout<<"# Sensitivity named as "<<data_name<<" does NOT exist.\n";
    exit(0);
  };

  int grp=sns_data[id].GetGroup();
  int grp_lib1=lib1.GetGroup();
  int grp_lib2=lib2.GetGroup();
  if(grp!=grp_lib1||grp!=grp_lib2){
    cout<<"# Error in SUTool::CalLibraryEffect.\n";
    cout<<"# Energy group is inconsistent.\n";
    cout<<"#   sensitivity data : "<<grp<<"\n";
    cout<<"#   library 1        : "<<grp_lib1<<"\n";
    cout<<"#   library 2        : "<<grp_lib2<<"\n";
    exit(0);
  };

  cout.setf(ios::showpoint);
  cout.precision(6);

  cout<<"#\n# Library effect calculation\n";
  cout<<"# data name is "<<data_name<<"\n#\n";
  real totval=0.;

  int size1d=sns_data[id].GetSize1D();
  int size2d=sns_data[id].GetSize2D();

  // +++ for 1D xs
  for(int i=0;i<size1d;i++){
    int mat=sns_data[id].GetMatList1D(i);
    int mt=sns_data[id].GetMtList1D(i);
    if(lib1.ExistLibData(mat)&&lib2.ExistLibData(mat)){
      enum xstype xst=GetXSType(mt);
      GroupData1D diff(grp);
      // (fission spectrum normalization)
      real sum1=0.;
      real sum2=0.;
      if(mt==181){
	sum1=lib1.GetLibData(mat).GetXSData().GetData1d(chi).get_sum();
	sum2=lib2.GetLibData(mat).GetXSData().GetData1d(chi).get_sum();
	sum1=1./sum1;
	sum2=1./sum2;
      };
      for(int j=0;j<grp;j++){
	real tmp1=lib1.GetLibData(mat).GetXSData().GetData1d(xst).get_dat(j);
	real tmp2=lib2.GetLibData(mat).GetXSData().GetData1d(xst).get_dat(j);
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
      effect=diff.mult(sns_data[id].GetSensitivity1D(mat,mt));
      for(int j=0;j<grp;j++){
	real tmp=effect.get_dat(j);
	if(tmp>0.){pos_effect+=tmp;}
	else{neg_effect+=tmp;};
      };
      real value=pos_effect+neg_effect;
      if(grp_wise_prt){
	cout<<"#\n#  mat/mt = "<<mat<<"/"<<mt<<"\n";
        for(int j=0;j<grp;j++){
	  real e1=lib1.GetEnband().get_dat(j);
	  real e2=lib1.GetEnband().get_dat(j+1);
          real letwid=log(e1/e2);
	  letwid/=0.25;
          cout<<"    "<<e1<<" "<<effect.get_dat(j)/letwid<<"\n";
	};
	cout<<"\n";
      };
      totval+=value;
      if(fabs(value)>threshold){
        cout<<"* mat:"<<mat<<" & mt:"<<mt<<" : "<<value;
        cout<<"   ("<<pos_effect<<", "<<neg_effect<<")\n";
      };
    };
  };

  // +++ for 2D xs
  for(int i=0;i<size2d;i++){
    int mat=sns_data[id].GetMatList2D(i);
    int mt=sns_data[id].GetMtList2D(i);
    if(lib1.ExistLibData(mat)&&lib2.ExistLibData(mat)){
      enum xstype xst=GetXSType(mt);
      GroupData2D diff(grp,grp);
      for(int j=0;j<grp;j++){
	for(int k=0;k<grp;k++){
  	  real tmp1=lib1.GetLibData(mat).GetXSData().GetData2d(xst).get_dat(j,k);
  	  real tmp2=lib2.GetLibData(mat).GetXSData().GetData2d(xst).get_dat(j,k);
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
      GroupData2D senstmp=sns_data[id].GetSensitivity2D(mat,mt);
      for(int j=0;j<grp;j++){
        real tmp=0.;
	for(int k=0;k<grp;k++){
	  tmp+=diff.get_dat(j,k)*senstmp.get_dat(j,k);
	};
	if(tmp>0.){pos_effect+=tmp;}
	else{neg_effect+=tmp;};
	effect[j]=tmp;
      };
      real value=pos_effect+neg_effect;
      if(grp_wise_prt){
	cout<<"#\n#  mat/mt = "<<mat<<"/"<<mt<<"\n";
        for(int j=0;j<grp;j++){
	  real e1=lib1.GetEnband().get_dat(j);
	  real e2=lib1.GetEnband().get_dat(j+1);
          real letwid=log(e1/e2);
	  letwid/=0.25;
          cout<<"    "<<e1<<" "<<effect[j]/letwid<<"\n";
	};
	cout<<"\n";
      };
      totval+=value;
      if(fabs(value)>threshold){      
        cout<<"* mat:"<<mat<<" & mt:"<<mt<<" : "<<value;
        cout<<"   ("<<pos_effect<<", "<<neg_effect<<")\n";
      };
    };
  };

  cout<<"\n* Total :"<<totval<<"\n\n";
};
