#include <cstdlib>
#include "GData.h"
#include<fstream>

using namespace std;

GData::GData(int num_y_inp, bool point_inp)
{
  PutNumY(num_y_inp);
  point=true;
};

void GData::PutNumY(int inp)
{
  num_y=inp;
  y.resize(inp);
  tag_y.resize(inp);
};

void GData::PutTagY(int i, string tag_inp)
{
  if(i<0||i>=num_y){
    cout<<"# Error in GData::PutTagY.\n";
    exit(0);
  };
  tag_y[i]=tag_inp;
};

void GData::push_back_y(int i, real yinp)
{
  if(i<0||i>=num_y){
    cout<<"# Error in GData::push_back_y.\n";
    exit(0);
  };
  y[i].push_back(yinp);
};

void GData::show_self()
{
  /*
  //cout<<"# The size of x : "<<x.size()<<"\n";
  cout<<"# Tag of x : "<<tag_x<<"\n";
  for(int i=0;i<num_y;i++){
    //cout<<"# The size of y"<<i<<" : "<<y[i].size()<<"\n";
    cout<<"# Tag of y"<<i<<" : "<<tag_y[i]<<"\n";
  };
  exit(0);    
  */

  cout<<"# "<<tag_x<<"    ";
  for(int j=0;j<num_y;j++){
    cout<<tag_y[j]<<"  ";
  };
  cout<<"\n";
  int num=x.size();
  for(int i=0;i<num;i++){
    cout<<"   "<<x[i]<<"    ";
    for(int j=0;j<num_y;j++){
      cout<<y[j][i]<<" ";
    };
    cout<<"\n";
  };
};

void GData::WriteFile(string mdir, string name)
{
  mdir.append(name);
  
  ofstream fout;
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Error in GData::WriteFile.\n";
    cout<<"# File ["<<mdir<<"] cannot be found.\n";
    exit(0);
  };

  fout.setf(ios::scientific);
  fout.precision(8);
  
  int sz=x.size();
  fout<<num_y<<"\n";   
  fout<<sz<<"\n";
  fout<<tag_x<<" "; 
  for(int i=0;i<num_y;i++){
    fout<<tag_y[i]<<" ";
  };
  fout<<"\n";

  for(int i=0;i<sz;i++){
    fout<<x[i]<<" ";
    for(int j=0;j<num_y;j++){
      fout<<y[j][i]<<" ";
   };
    fout<<"\n";
  };

  fout.close();
};

void GData::ReadFile(string mdir, string name)
{
  mdir.append(name);
  
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in GData::ReadFile.\n";
    cout<<"# File ["<<mdir<<"] cannot be found.\n";
    exit(0);
  };

  int tmp;
  fin>>tmp;
  PutNumY(tmp);

  int sz;
  fin>>sz;

  string dummy;
  fin>>dummy;
  PutTagX(dummy);

  for(int i=0;i<num_y;i++){
    fin>>dummy;
    PutTagY(i,dummy);
  };

  real vtmp;
  for(int i=0;i<sz;i++){
    fin>>vtmp;
    push_back_x(vtmp);
    for(int j=0;j<num_y;j++){
      fin>>vtmp;
      push_back_y(j,vtmp);
    };
  };

  fin.close();
};

void GData::SlideX(real inp)
{
  int sz=x.size();
  for(int i=0;i<sz;i++){
    x[i]+=inp;
  };
};

void GData::FactorizeX(real inp)
{
  int sz=x.size();
  for(int i=0;i<sz;i++){
    x[i]*=inp;
  };
};

real GData::GetY_Lin(real x_inp, int id_y)
{
  if(id_y<0||id_y>=num_y){
    cout<<"# Error in GData::GetY_Lin.\n";
    exit(0);
  };

  int sz=x.size();
  for(int j=0;j<sz-1;j++){
    if(x_inp>=x[j]&&x_inp<x[j+1]){
      real x2=x[j+1];
      real x1=x[j];
      real xx=x_inp;
      real y_ret=y[id_y][j]+(y[id_y][j+1]-y[id_y][j])/(x2-x1)*(xx-x1);
      return y_ret;
    };
  };

  cout<<"# Error in GData::GetY_Lin.\n";
  cout<<"# X "<<x_inp<<" is out of the data range.\n";
  exit(0);
};

real GData::GetY_Log(real x_inp, int id_y)
{
  if(id_y<0||id_y>=num_y){
    cout<<"# Error in GData::GetY_Log.\n";
    exit(0);
  };

  int sz=x.size();
  for(int j=0;j<sz-1;j++){
    if(x_inp>=x[j]&&x_inp<x[j+1]){
      real x2=log(x[j+1]);
      real x1=log(x[j]);
      real xx=log(x_inp);
      real y_ret=y[id_y][j]+(y[id_y][j+1]-y[id_y][j])/(x2-x1)*(xx-x1);
      return y_ret;
    };
  };

  cout<<"# Error in GData::GetY_Log.\n";
  cout<<"# X "<<x_inp<<" is out of the data range.\n";
  exit(0);
};

void GData::DeleteData(real x_val, real thresh)
{
  int sz=x.size();
  for(int i=0;i<sz;i++){
    if(fabs(x[i]-x_val)/x_val<thresh){
      x.erase(x.begin()+i);
      for(int j=0;j<num_y;j++){
	y[j].erase(y[j].begin()+i);
      };
      return;
    };
  };
};

vector<real> GData::GetY(int i)
{
  if(i<0||i>=num_y){
    cout<<"# Error in GData::GetY\n";
    cout<<"# Inappropriate index "<<i<<" is assgined.\n";
    exit(0);
  };
  return y[i];
};

void GData::TakeDifferenceFrom(GData &ref, bool absolute)
{
  bool consistency=true;
  if(num_y!=ref.GetNumY())consistency=false;
  int data_size=x.size();
  int data_size_ref=ref.GetDataNum();
  if(data_size!=data_size_ref)consistency=false;
  if(!consistency){
    cout<<"# Error in GData::TakeAbsoluteDifferenceFrom.\n";
    cout<<"# Data size is inconsistent with each other.\n";
    exit(0);
  };

  for(int j=0;j<num_y;j++){
    vector<real> y_ref=ref.GetY(j);
    for(int i=0;i<data_size;i++){
      if(absolute){
        y[j][i]-=y_ref[i];
      }else{
	if(y_ref[i]==0.){
	  y[j][i]=0.;
	}else{
          y[j][i]=(y[j][i]-y_ref[i])/y_ref[i];	  
	};
      };
    };
  };
};

vector<real> GData::GetAbsoluteMaximumValue()
{
  vector<real> absmax(num_y);
  for(int i=0;i<num_y;i++){
    real tmp=0.;
    for(int j=0;j<x.size();j++){
      if(abs(y[i][j])>tmp)tmp=abs(y[i][j]);
    };
    absmax[i]=tmp;
  };
  return absmax;
};

