#include <cstdlib>
#include <iostream>
#include "GroupDataSet.h"

using namespace std;

GroupDataSet::~GroupDataSet()
{
  AllVectorClear();
};

void GroupDataSet::Init(string naminp)
{
  name=naminp;

  if(name=="MacroCrossSection"){
    PutNum1d(12);
    map1d[sigt]=0;
    map1d[nusigf]=1;
    map1d[chi]=2;
    map1d[sigtr]=3;
    map1d[sigr]=4;
    map1d[siga]=5;
    map1d[d]=6;
    map1d[dr]=7; // perpendicular
    map1d[dz]=8; // parallel
    map1d[sign2n]=9;
    map1d[sigst]=10;
    map1d[kerma]=11; // Kerma factor
    for(int i=0;i<num1d;i++){
      PutDim1d(i,1);
    };
    PutNum2d(2);
    map2d[sigs]=0;
    map2d[chi]=1;
    PutDim2d(0,1);
    PutDim2d(1,1);
    return;
  };

  if(name=="MicroCrossSection"){
    PutNum1d(9);
    map1d[sigf]=0;
    map1d[sigc]=1;
    map1d[nu]=2;
    map1d[mu]=3;
    map1d[sigel]=4;
    map1d[chi]=5;
    map1d[siginel]=6;
    map1d[sign2n]=7;
    map1d[sigt]=8;
    for(int i=0;i<num1d;i++){
      PutDim1d(i,1);
    };
    PutDim1d(8,2);
    PutNum2d(3);
    map2d[sigel]=0;
    map2d[siginel]=1;
    map2d[sign2n]=2;
    PutDim2d(0,2);
    PutDim2d(1,2);
    PutDim2d(2,2);
    return;
  };

  if(name=="MicroCrossSectionSimple"){
    PutNum1d(3);
    map1d[sigf]=0;
    map1d[sigc]=1;
    map1d[nu]=2;
    for(int i=0;i<num1d;i++){
      PutDim1d(i,1);
    };
    return;
  };

  if(name=="GammaCrossSection"){
    PutNum1d(6);
    map1d[sigt]=0;
    map1d[coh_el]=1;
    map1d[incoh_el]=2;
    map1d[pair_pro]=3;
    map1d[siga]=4;
    map1d[kerma]=5;
    for(int i=0;i<num1d;i++){
      PutDim1d(i,1);
    };
    PutNum2d(3);
    map2d[coh_el]=0;
    map2d[incoh_el]=1;
    map2d[pair_pro]=2;
    PutDim2d(0,2);
    PutDim2d(1,2);
    PutDim2d(2,2);
    return;
  };

  cout<<"Error in Init of GroupDataSet!.\n";
  cout<<" (name ="<<name<<")\n";
  exit(0);
};

void GroupDataSet::PutGrp(int g)
{
  grp=g;
  for(int i=0;i<num1d;i++){
    for(int j=0;j<dim1d[i];j++){
      data1d[i][j].put_imax(g);
      data1d[i][j].set_zero();
    };
  };

  for(int i=0;i<num2d;i++){
    for(int j=0;j<dim2d[i];j++){
      data2d[i][j].put_yx(g,g);
      data2d[i][j].set_zero();
    };
  };
};

void GroupDataSet::PutNum1d(int i)
{
  num1d=i;
  data1d.resize(i);
  dim1d.resize(i,0);
};

void GroupDataSet::PutNum2d(int i)
{
  num2d=i;
  if(i>0){
    data2d.resize(i);
    dim2d.resize(i,0);
  };
};

void GroupDataSet::PutDim1d(int l,int m)
{
  dim1d[l]=m;
  data1d[l].resize(m);
  if(grp!=-1){
    for(int i=0;i<m;i++){
      data1d[l][i].put_imax(grp);
    };
  };
};

void GroupDataSet::PutDim2d(int l,int m)
{
  dim2d[l]=m;
  data2d[l].resize(m);
  if(grp!=-1){
    for(int i=0;i<m;i++){
      data2d[l][i].put_yx(grp,grp);
    };
  };
};

void GroupDataSet::PutDim1d(enum xstype name,int l)
{
  PutDim1d(GetPos1d(name),l);
};

void GroupDataSet::PutDim2d(enum xstype name,int l)
{
  PutDim2d(GetPos2d(name),l);
};

GroupData1D &GroupDataSet::GetData1d(int i,int l)
{
  /*
  if(dim1d[i]<=l){
    cout<<"Error in GetData1D of GroupDataSet.\n";
    cout<<"maximam order = "<<dim1d[i]-1<<"\n";
    cout<<"requested     = "<<l<<"\n";
    cout<<" (Reaction name is "<<name<<")\n";	
  };
  */
  return data1d[i][l];
};

GroupData2D &GroupDataSet::GetData2d(int i,int l)
{
  /*
  if(dim2d[i]<=l){
    cout<<"Error in GetData2D of GroupDataSet.\n";
    cout<<"maximum order = "<<dim2d[i]-1<<"\n";
    cout<<"requested     = "<<l<<"\n";
    cout<<" (Reaction name is "<<name<<")\n";	
  };
  */
  return data2d[i][l];
};

GroupData1D &GroupDataSet::GetData1d(enum xstype name,int l)
{
  return GetData1d(GetPos1d(name),l);
};

GroupData2D &GroupDataSet::GetData2d(enum xstype name,int l)
{
  return GetData2d(GetPos2d(name),l);
};

GroupDataSet GroupDataSet::Cond(int ngrp,int *bgrp,GroupData1D fl,GroupData1D cu)
{
  GroupDataSet ret(name);
  ret.PutGrp(ngrp);

  for(int i=1;i<ngrp;i++){
    if(bgrp[i]<=bgrp[i-1]){
      cout<<"# Error in Cond of GroupDataSet.cxx.\n";
      cout<<"# Error in group structure.\n";
      exit(0);
    };
  };
  if(bgrp[ngrp-1]!=grp-1){
    cout<<"# Error in Cond of GroupDataSet.cxx\n";
    cout<<"# The lowest energy boundary is inconsistent!\n";
    cout<<"#  Energy group before condensation ="<<grp<<"\n";
    cout<<"#  Incorrect lowest energy boundary ="<<bgrp[ngrp-1]<<"\n";
    exit(0);
  };

  // check for sign of current
  for(int i=0;i<ngrp;i++){
    int is;
    if(i==0){is=0;}
    else{is=bgrp[i-1]+1;};
    int ie=bgrp[i];
    bool pos=false;
    bool neg=false;
    for(int j=is;j<=ie;j++){
      if(cu.get_dat(j)>0.)pos=true;
      if(cu.get_dat(j)<0.)neg=true;
    };
    if(pos&&neg){
      //cout<<"Positive & Negative current in the same broad group "<<i<<"\n";
    };
  };


  for(int i=0;i<grp;i++){
    real cuabs=fabs(cu.get_dat(i));
    cu.put_data(i,cuabs);
  };


  map<enum xstype,int>::iterator itr;
  for(itr=map1d.begin();itr!=map1d.end();itr++){
    int i=itr->second;
    ret.PutDim1d(i,dim1d[i]);
    if(itr->first==chi){
      ret.GetData1d(i,0).copy(GetData1d(i,0).CondSum(ngrp,bgrp));
    }else if(itr->first==sigtr){
      ret.GetData1d(i,0).copy(GetData1d(i,0).Cond(cu,ngrp,bgrp));
    }else if(itr->first==nu){
      GroupData1D tmp(grp);
      tmp=fl.mult(GetData1d(sigf));
      ret.GetData1d(i,0).copy(GetData1d(i,0).Cond(tmp,ngrp,bgrp));
    }else{
      ret.GetData1d(i,0).copy(GetData1d(i,0).Cond(fl,ngrp,bgrp));
      for(int j=1;j<dim1d[i];j++){
        ret.GetData1d(i,j).copy(GetData1d(i,j).Cond(cu,ngrp,bgrp));
      };
    };
  };

  for(itr=map2d.begin();itr!=map2d.end();itr++){
    int i=itr->second;
    ret.PutDim2d(i,dim2d[i]);
    ret.GetData2d(i,0).copy(GetData2d(i,0).Cond(fl,ngrp,bgrp));
    for(int j=1;j<dim2d[i];j++){
      ret.GetData2d(i,j).copy(GetData2d(i,j).Cond(cu,ngrp,bgrp));
    };
  };

  return ret;
};

void GroupDataSet::ShowSelf()
{
  map<enum xstype,int>::iterator itr;
  cout<<"#\n";
  cout<<"# 1D cross section data\n";
  cout<<"#\n";
  for(itr=map1d.begin();itr!=map1d.end();itr++){
    int i=itr->second;
    //cout<<itr->first<<"\n";
    cout<<xstype_string[itr->first]<<"\n";    
    for(int j=0;j<dim1d[i];j++){
      cout<<"  order : "<<j<<"\n";
      //data1d[i][j].show_self();
    };
    cout<<"\n";
  };

  cout<<"#\n";
  cout<<"# 2D cross section data\n";
  cout<<"#\n";
  for(itr=map2d.begin();itr!=map2d.end();itr++){
    int i=itr->second;
    //cout<<itr->first<<"\n";
    cout<<xstype_string[itr->first]<<"\n";    
    for(int j=0;j<dim2d[i];j++){
      cout<<"  order : "<<j<<"\n";
      //data2d[i][j].show_self();
    };
    cout<<"\n";
  };
};

void GroupDataSet::ShowData1D(GroupData1D &enband)
{
  if(grp!=enband.get_imax()-1){
    cout<<"# Error in GroupDataSet::ShowData1D\n";
    cout<<"# The number of energy groups is inconsistent.\n";
    cout<<"#   GroupDataSet    : "<<grp<<"-group\n";
    cout<<"#   Energy boundary : "<<enband.get_imax()-1<<"\n";
    exit(0);
  };
  cout<<"#\n# 1D-data stored in GroupDataSet\n#\n";

  map<enum xstype,int>::iterator it;
  int num=0;
  for(it=map1d.begin();it!=map1d.end();it++){
    num++;
  };
  vector<bool> non_zero(num,false);

  int index=0;
  for(it=map1d.begin();it!=map1d.end();it++){
    int i=it->second;
    for(int g=0;g<grp;g++){
      if(data1d[i][0].get_dat(g)!=0.){
        non_zero[index]=true;
	break;
      };
    };
    index++;
  };

  cout<<"#E_up     ";
  index=0;
  for(it=map1d.begin();it!=map1d.end();it++){
    if(non_zero[index]){
      enum xstype xt=it->first;
      cout<<xstype_string[xt]<<" ";
    };
    index++;
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  
  for(int g=0;g<grp+1;g++){
    cout<<enband.get_dat(g)<<" ";
    int gg=g;
    if(g==grp)gg=g-1;
    index=0;
    for(it=map1d.begin();it!=map1d.end();it++){
      if(non_zero[index]){
        int i=it->second;
        cout<<data1d[i][0].get_dat(gg)<<" ";
      };
      index++;
    };
    cout<<"\n";
  };
};


int GroupDataSet::GetPos1d(enum xstype name)
{
  map<enum xstype,int>::iterator iter;
  iter=map1d.find(name);
  if(iter==map1d.end()){
    cout<<"# Error in GetData1d of GroupDataSet.\n";
    cout<<"# There is not a member you requested.\n";
    cout<<"# You requested "<<name<<"\n";
    exit(0);
  };
  return iter->second;

  //return map1d[name];
};

int GroupDataSet::GetPos2d(enum xstype name)
{

  map<enum xstype,int>::iterator iter;
  iter=map2d.find(name);
  if(iter==map1d.end()){
    cout<<"# Error in GetData2d of GroupDataSet.\n";
    cout<<"# There is not a member you requested.\n";
    cout<<"# You requested "<<name<<"\n";
    exit(0);
  };
  return iter->second;

  //return map2d[name];
};

bool GroupDataSet::CheckSameType(GroupDataSet &sec)
{
  if(name!=sec.GetName()){
    cout<<"Warning in GroupDataSet::CheckSameType.\n";
    cout<<"Difference type of GroupDataSet.\n";
    cout<<" One is "<<name<<".\n";
    cout<<" The other is "<<sec.GetName()<<".\n";
    return false;
  };
  return true;
};

void GroupDataSet::DataCopy(GroupDataSet &sec)
{
  if(!CheckSameType(sec))exit(0);

  for(int i=0;i<num1d;i++){
    int tmp=sec.GetDim1d(i);
    PutDim1d(i,tmp);
    for(int l=0;l<dim1d[i];l++){
      if(sec.GetDim1d(i)>l){
	data1d[i][l].copy(sec.GetData1d(i,l));
      };
    };
  };
  for(int i=0;i<num2d;i++){
    int tmp=sec.GetDim2d(i);
    if(tmp!=0){
      PutDim2d(i,tmp);
      for(int l=0;l<dim2d[i];l++){
        if(sec.GetDim2d(i)>l){
	  data2d[i][l].copy(sec.GetData2d(i,l));
        };
      };
    };
  };
};

void GroupDataSet::DataCopyPL(GroupDataSet &sec,int plmax_2d)
{
  if(!CheckSameType(sec))exit(0);

  for(int i=0;i<num1d;i++){
    int tmp=sec.GetDim1d(i);
    PutDim1d(i,tmp);
    for(int l=0;l<dim1d[i];l++){
      if(sec.GetDim1d(i)>l){
	data1d[i][l].copy(sec.GetData1d(i,l));
      };
    };
  };
  for(int i=0;i<num2d;i++){
    int tmp=sec.GetDim2d(i);
    if(tmp!=0){
      if(tmp>plmax_2d+1)tmp=plmax_2d+1;
      PutDim2d(i,tmp);
      for(int l=0;l<dim2d[i];l++){
        if(sec.GetDim2d(i)>l){
	  data2d[i][l].copy(sec.GetData2d(i,l));
        };
      };
    };
  };
};

void GroupDataSet::DataCopy(GroupDataSet &sec,int lwgrp)
{
  DataCopy(sec,0,lwgrp);
};

void GroupDataSet::DataCopy(GroupDataSet &sec,int upgrp,int lwgrp)
{
  if(!CheckSameType(sec))exit(0);

  if(lwgrp>=grp){
    cout<<"Error in GroupDataSet::DataCopy.\n";
    cout<<"Lower group is bigger than group.\n";
    cout<<"Requested lower group : "<<lwgrp<<"\n";
    cout<<"Group in GroupDataSet : "<<grp<<"\n";
    exit(0);
  };
  for(int i=0;i<num1d;i++){
    int tmp=sec.GetDim1d(i);
    PutDim1d(i,tmp);
    for(int l=0;l<dim1d[i];l++){
      if(sec.GetDim1d(i)>l){
	for(int k=upgrp;k<=lwgrp;k++){
  	  data1d[i][l].put_data(k,sec.GetData1d(i,l).get_dat(k));
	};
      };
    };
  };

  int grp2=sec.GetGrp();
  for(int i=0;i<num2d;i++){
    int tmp=sec.GetDim2d(i);
    PutDim2d(i,tmp);
    for(int l=0;l<dim2d[i];l++){
      if(sec.GetDim2d(i)>l){
	for(int k=upgrp;k<=lwgrp;k++){
	  //for(int kk=0;kk<=lwgrp;kk++){
	  for(int kk=0;kk<grp2;kk++){
   	    data2d[i][l].put_data(k,kk,sec.GetData2d(i,l).get_dat(k,kk));
	  };
	};
      };
    };
  };
};

void GroupDataSet::MultiplyExceptForChi(real factor,int g)
{
  for(int i=0;i<num1d;i++){
    if(name!="MacroCrossSection"||i!=2){
    for(int l=0;l<dim1d[i];l++){
      real org=data1d[i][l].get_dat(g);
      data1d[i][l].put_data(g,org*factor);
    };
    };
  };
  for(int i=0;i<num2d;i++){
    for(int l=0;l<dim2d[i];l++){
      for(int j=0;j<grp;j++){
	real org=data2d[i][l].get_dat(g,j);
        data2d[i][l].put_data(g,j,org*factor);
      };
    };
  };
};

void GroupDataSet::Multiply(real factor,int g)
{
  for(int i=0;i<num1d;i++){
    for(int l=0;l<dim1d[i];l++){
      real org=data1d[i][l].get_dat(g);
      data1d[i][l].put_data(g,org*factor);
    };
  };
  for(int i=0;i<num2d;i++){
    for(int l=0;l<dim2d[i];l++){
      for(int j=0;j<grp;j++){
	real org=data2d[i][l].get_dat(g,j);
        data2d[i][l].put_data(g,j,org*factor);
      };
    };
  };
};

void GroupDataSet::Multiply(real factor)
{
  for(int g=0;g<grp;g++){
    Multiply(factor,g);
  };
};

void GroupDataSet::Add(GroupDataSet &sec)
{
  if(!CheckSameType(sec))exit(0);

  for(int i=0;i<num1d;i++){
    if(dim1d[i]!=sec.GetDim1d(i)){
      cout<<"Warning in GroupDataSet::Add\n";
      cout<<"Legendre component is inconsistent in 1D XS.\n";
    };
    for(int l=0;l<dim1d[i];l++){
      data1d[i][l]=data1d[i][l]+sec.GetData1d(i,l);
    };
  };
  for(int i=0;i<num2d;i++){
    int tmp=dim2d[i];
    if(tmp<sec.GetDim2d(i)){
      cout<<"Warning in GroupDataSet::Add\n";
      cout<<"Higher Legendre component is NOT added.\n";
    };
    if(sec.GetDim2d(i)<tmp)tmp=sec.GetDim2d(i);
    for(int l=0;l<tmp;l++){
      data2d[i][l]=data2d[i][l]+sec.GetData2d(i,l);
    };
  };
};

void GroupDataSet::AllVectorClear()
{
  /*
  dim1d.clear();
  dim2d.clear();
  data1d.clear();
  data2d.clear();
  */
  //cout<<"# !!! "<<num1d<<" "<<num2d<<" "<<grp<<" "<<name<<"\n";
  for(int i=0;i<num1d;i++){
    vector<GroupData1D> ().swap(data1d[i]);
  };
  vector<int> ().swap(dim1d);
  num1d=0;

  AllVectorClear2DData();

  /*
  for(int i=0;i<num2d;i++){
    vector<GroupData2D> ().swap(data2d[i]);
  };
  vector<int> ().swap(dim2d);
  num2d=0;
  */
  //vector< vector<GroupData1D> > ().swap(data1d);
  //vector< vector<GroupData2D> > ().swap(data2d);

};

void GroupDataSet::AllVectorClear2DData()
{
  /*
  dim2d.clear();
  data2d.clear();
  */

  for(int i=0;i<num2d;i++){
    vector<GroupData2D> ().swap(data2d[i]);
  };
  vector<int> ().swap(dim2d);
  num2d=0;

  /*
  vector<int> ().swap(dim2d);
  vector< vector<GroupData2D> > ().swap(data2d);
  */

};

void GroupDataSet::Allzero()
{
  for(int g=0;g<grp;g++){
    for(int i=0;i<num1d;i++){
      for(int l=0;l<dim1d[i];l++){
	data1d[i][l].put_data(g,0.);
      };
    };
    for(int i=0;i<num2d;i++){
      for(int l=0;l<dim2d[i];l++){
	for(int j=0;j<grp;j++){
	  data2d[i][l].put_data(g,j,0.);
	};
      };
    };
  };
};


void GroupDataSet::XSMatrixNormalizationTo1DXS()
{
  // !! Caution !!
  //
  // The sum of (n,2n) secondary cross sections should be HALF
  // of the 1D (n,2n) reaction cross section.

  if(name!="MicroCrossSection"){
    cout<<"# Warning in GroupDataSet::2DXSNormalizationTo1DXS.\n";
    cout<<"# Data type should be micro cross section.\n";
    return;
  };

  enum xstype tt[]={sigel,siginel,sign2n};
  for(int i=0;i<3;i++){
    GroupData1D xs1d=GetData1d(tt[i]);
    for(int g=0;g<grp;g++){
      if(xs1d.get_dat(g)>0.){
        real sum=0.;
        for(int g2=0;g2<grp;g2++){
  	  sum+=GetData2d(tt[i]).get_dat(g,g2);
        };
	if(i==2)sum*=0.5; // for (n,2n)
	if(sum>0.){
          real ratio=xs1d.get_dat(g)/sum;
          for(int pl=0;pl<GetDim2d(tt[i]);pl++){
            for(int g2=0;g2<grp;g2++){
      	      real orgxs=GetData2d(tt[i],pl).get_dat(g,g2);
              GetData2d(tt[i],pl).put_data(g,g2,orgxs*ratio);
	    };
	  };
	};
      };
    };
  };

};

void GroupDataSet::P1MatrixNormalizationToMubar()
{
  if(name!="MicroCrossSection"){
    cout<<"# Warning in GroupDataSet::P1MatrixNormalizationToMubar.\n";
    cout<<"# Data type should be micro cross section.\n";
    return;
  };

  for(int g=0;g<grp;g++){
    real vmu=GetData1d(mu).get_dat(g);
    real sum0=0.;
    real sum1=0.;
    for(int g2=0;g2<grp;g2++){
      sum0+=GetData2d(sigel).get_dat(g,g2);
      sum1+=GetData2d(sigel,1).get_dat(g,g2);
    };
    real vmu2=sum1/sum0/3.;
    real ratio=vmu/vmu2;
    for(int g2=0;g2<grp;g2++){
      real org=GetData2d(sigel,1).get_dat(g,g2);
      GetData2d(sigel,1).put_data(g,g2,org*ratio);
    };
  };

};

// The following are implemented by Yanagihara-kun for the NFI joint research

GroupDataSet GroupDataSet::Caldsdx1(GroupDataSet &a, GroupDataSet &b, real c)
{  
  GroupDataSet delta=a;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    delta.PutDim1d(l,dim1d);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D delta_xs=(xs2-xs1)*(1./c);
      for(int i=0;i<grp;i++){
	inp=delta_xs.get_dat(i);
	delta.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    delta.PutDim2d(l,dim2d);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs=(xs2-xs1)*(1./c);
      for(int i=0;i<grp;i++){
	for(int j=0;j<grp;j++){
	  inp=delta_xs.get_dat(i,j);
	  delta.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  GroupDataSet ref=delta;
  return ref;
};

GroupDataSet GroupDataSet::Caldsdx2(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d, bool nu)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=a;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1;
      GroupData1D xs2;
      GroupData1D xs3;      
      if(l==0&&nu){
	xs1=a.GetData1d(l,m).mult(a.GetData1d(2,0)); //Read for base Macro XSfile1
	xs2=b.GetData1d(l,m).mult(b.GetData1d(2,0)); //Read for base Macro XSfile2
	xs3=c.GetData1d(l,m).mult(c.GetData1d(2,0)); //Read for base Macro XSfile3
      }else{
	xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
	xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
	xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      };
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
	inp=delta_xs1.get_dat(i);
	delta1.GetData1d(l,m).put_data(i,inp);
	inp=delta_xs2.get_dat(i);
	delta2.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
	  inp=delta_xs1.get_dat(i,j);
	  delta1.GetData2d(l,m).put_data(i,j,inp);
	  inp=delta_xs2.get_dat(i,j);
	  delta2.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  GroupDataSet ret=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      for(int j=0;j<grp;j++){
	double tmp=y.get_dat(j);
	double tmp1=y1.get_dat(j);
	real tmp2[]={tmp,0.,tmp1};
	double tmpa=0.;
	double tmpb=0.;
	LeastSquaresMethod(3,d,tmp2,tmpa,tmpb);
	ret.GetData1d(l,m).put_data(j,tmpa);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      for(int j=0;j<grp;j++){
	for(int k=0;k<grp;k++){
	  double tmp=y.get_dat(j,k);
	  double tmp1=y1.get_dat(j,k);
	  real tmp2[]={tmp,0.,tmp1};
	  double tmpa=0.;
	  double tmpb=0.;
	  LeastSquaresMethod(3,d,tmp2,tmpa,tmpb);
	  ret.GetData2d(l,m).put_data(j,k,tmpa);
	};
      };
    };
  };  
  GroupDataSet ref=ret;
  return ref;
};

GroupDataSet GroupDataSet::Caldsdx3(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
   
  GroupDataSet ret=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        real tmp3[]={tmp,0.,tmp1,tmp2};
        double tmpa=0.;
        double tmpb=0.;
        LeastSquaresMethod(4,e,tmp3,tmpa,tmpb);
        ret.GetData1d(l,m).put_data(j,tmpa);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  real tmp3[]={tmp,0.,tmp1,tmp2};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  LeastSquaresMethod(4,e,tmp3,tmpa,tmpb);
      	  ret.GetData2d(l,m).put_data(j,k,tmpa);
	};
      };
    };
  };  
  GroupDataSet ref=ret;
  return ret;
};

GroupDataSet GroupDataSet::Caldsdx4(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, real *f)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  GroupDataSet delta4=d;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs5=e.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      GroupData1D delta_xs4=xs5-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs4.get_dat(i);
        delta4.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs5=e.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      GroupData2D delta_xs4=xs5-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs4.get_dat(i,j);
      	  delta4.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  GroupDataSet ret=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      GroupData1D y3=delta4.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        double tmp3=y3.get_dat(j);
        real tmp4[]={tmp,0.,tmp1,tmp2,tmp3};
        double tmpa=0.;
        double tmpb=0.;
        LeastSquaresMethod(5,f,tmp4,tmpa,tmpb);
        ret.GetData1d(l,m).put_data(j,tmpa);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      GroupData2D y3=delta4.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  double tmp3=y3.get_dat(j,k);
      	  real tmp4[]={tmp,0.,tmp1,tmp2,tmp3};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  LeastSquaresMethod(5,f,tmp4,tmpa,tmpb);
      	  ret.GetData2d(l,m).put_data(j,k,tmpa);
	};
      };
    };
  };  
  GroupDataSet ref=ret;
  return ret;  
};

GroupDataSet GroupDataSet::Caldsdx5(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, real *g)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  GroupDataSet delta4=d;
  GroupDataSet delta5=e;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs5=e.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs6=f.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      GroupData1D delta_xs4=xs5-xs1;
      GroupData1D delta_xs5=xs6-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs4.get_dat(i);
        delta4.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs5.get_dat(i);
        delta5.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs5=e.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs6=f.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      GroupData2D delta_xs4=xs5-xs1;
      GroupData2D delta_xs5=xs6-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs4.get_dat(i,j);
      	  delta4.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs5.get_dat(i,j);
      	  delta5.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  GroupDataSet ret=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      GroupData1D y3=delta4.GetData1d(l,m);
      GroupData1D y4=delta5.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        double tmp3=y3.get_dat(j);
        double tmp4=y4.get_dat(j);
        real tmp5[]={tmp,0.,tmp1,tmp2,tmp3,tmp4};
        double tmpa=0.;
        double tmpb=0.;
        LeastSquaresMethod(6,g,tmp5,tmpa,tmpb);
        ret.GetData1d(l,m).put_data(j,tmpa);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      GroupData2D y3=delta4.GetData2d(l,m);
      GroupData2D y4=delta5.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  double tmp3=y3.get_dat(j,k);
      	  double tmp4=y4.get_dat(j,k);
      	  real tmp5[]={tmp,0.,tmp1,tmp2,tmp3,tmp4};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  LeastSquaresMethod(6,g,tmp5,tmpa,tmpb);
      	  ret.GetData2d(l,m).put_data(j,k,tmpa);
	};
      };
    };
  };  
  GroupDataSet ref=ret;
  return ret;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx1Quadr(GroupDataSet &a, GroupDataSet &b, real c)
{  
  GroupDataSet delta=a;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1 
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D delta_xs1=xs2-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  real cc[]={c,0.};
  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        real tmp2[]={tmp,0.};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        LeastSquaresMethodQuadr(2,cc,tmp2,tmpa,tmpb,tmpc);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  real tmp2[]={tmp,0.};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
      	  LeastSquaresMethodQuadr(2,cc,tmp2,tmpa,tmpb,tmpc);
      	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
      	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx2Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=a;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1 
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
	  inp=delta_xs1.get_dat(i,j);
	  delta1.GetData2d(l,m).put_data(i,j,inp);
	  inp=delta_xs2.get_dat(i,j);
	  delta2.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };  
  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        real tmp2[]={tmp,0.,tmp1};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
	//if(tmp!=0&&tmp1!=0){
	LeastSquaresMethodQuadr(3,d,tmp2,tmpa,tmpb,tmpc);
	ret[0].GetData1d(l,m).put_data(j,tmpa);
	ret[1].GetData1d(l,m).put_data(j,tmpb);
	//}else{
	//  ret[0].GetData1d(l,m).put_data(j,0);
	//  ret[1].GetData1d(l,m).put_data(j,0);
	//};
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  real tmp2[]={tmp,0.,tmp1};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
	  //if(tmp!=0&&tmp1!=0){
	  LeastSquaresMethodQuadr(3,d,tmp2,tmpa,tmpb,tmpc);
	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	  //}else{
	  //  ret[0].GetData2d(l,m).put_data(j,k,0);
	  //  ret[1].GetData2d(l,m).put_data(j,k,0);
	  //};
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx3Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };

  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        real tmp3[]={tmp,0.,tmp1,tmp2};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        LeastSquaresMethodQuadr(4,e,tmp3,tmpa,tmpb,tmpc);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  real tmp3[]={tmp,0.,tmp1,tmp2};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
      	  LeastSquaresMethodQuadr(4,e,tmp3,tmpa,tmpb,tmpc);
      	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
      	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx5Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, real *g)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  GroupDataSet delta4=d;
  GroupDataSet delta5=e;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs5=e.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs6=f.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      GroupData1D delta_xs4=xs5-xs1;
      GroupData1D delta_xs5=xs6-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs4.get_dat(i);
        delta4.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs5.get_dat(i);
        delta5.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs5=e.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs6=f.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      GroupData2D delta_xs4=xs5-xs1;
      GroupData2D delta_xs5=xs6-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs4.get_dat(i,j);
      	  delta4.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs5.get_dat(i,j);
      	  delta5.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      GroupData1D y3=delta4.GetData1d(l,m);
      GroupData1D y4=delta5.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        double tmp3=y3.get_dat(j);
        double tmp4=y4.get_dat(j);
        real tmp5[]={tmp,0.,tmp1,tmp2,tmp3,tmp4};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        LeastSquaresMethodQuadr(6,g,tmp5,tmpa,tmpb,tmpc);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      GroupData2D y3=delta4.GetData2d(l,m);
      GroupData2D y4=delta5.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  double tmp3=y3.get_dat(j,k);
      	  double tmp4=y4.get_dat(j,k);
      	  real tmp5[]={tmp,0.,tmp1,tmp2,tmp3,tmp4};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
      	  LeastSquaresMethodQuadr(6,g,tmp5,tmpa,tmpb,tmpc);
      	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
      	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx6Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, GroupDataSet &g, real *h)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  GroupDataSet delta4=d;
  GroupDataSet delta5=e;
  GroupDataSet delta6=g;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs5=e.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs6=f.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs7=g.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      GroupData1D delta_xs4=xs5-xs1;
      GroupData1D delta_xs5=xs6-xs1;
      GroupData1D delta_xs6=xs7-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs4.get_dat(i);
        delta4.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs5.get_dat(i);
        delta5.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs6.get_dat(i);
        delta6.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs5=e.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs6=f.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs7=g.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      GroupData2D delta_xs4=xs5-xs1;
      GroupData2D delta_xs5=xs6-xs1;
      GroupData2D delta_xs6=xs7-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
  	  inp=delta_xs1.get_dat(i,j);
  	  delta1.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs2.get_dat(i,j);
  	  delta2.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs3.get_dat(i,j);
  	  delta3.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs4.get_dat(i,j);
  	  delta4.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs5.get_dat(i,j);
  	  delta5.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs6.get_dat(i,j);
  	  delta6.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      GroupData1D y3=delta4.GetData1d(l,m);
      GroupData1D y4=delta5.GetData1d(l,m);
      GroupData1D y5=delta6.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        double tmp3=y3.get_dat(j);
        double tmp4=y4.get_dat(j);
        double tmp5=y5.get_dat(j);
        real tmp6[]={tmp,0.,tmp1,tmp2,tmp3,tmp4,tmp5};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        LeastSquaresMethodQuadr(7,h,tmp6,tmpa,tmpb,tmpc);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      GroupData2D y3=delta4.GetData2d(l,m);
      GroupData2D y4=delta5.GetData2d(l,m);
      GroupData2D y5=delta6.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
  	  double tmp=y.get_dat(j,k);
  	  double tmp1=y1.get_dat(j,k);
  	  double tmp2=y2.get_dat(j,k);
  	  double tmp3=y3.get_dat(j,k);
  	  double tmp4=y4.get_dat(j,k);
  	  double tmp5=y5.get_dat(j,k);
  	  real tmp6[]={tmp,0.,tmp1,tmp2,tmp3,tmp4,tmp5};
  	  double tmpa=0.;
  	  double tmpb=0.;
  	  double tmpc=0.;
  	  LeastSquaresMethodQuadr(7,h,tmp6,tmpa,tmpb,tmpc);
  	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
  	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx7Quadr(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, GroupDataSet &e, GroupDataSet &f, GroupDataSet &g, GroupDataSet &h, real *i)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  GroupDataSet delta4=d;
  GroupDataSet delta5=e;
  GroupDataSet delta6=f;
  GroupDataSet delta7=g;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs5=e.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs6=f.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs7=g.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D xs8=h.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      GroupData1D delta_xs4=xs5-xs1;
      GroupData1D delta_xs5=xs6-xs1;
      GroupData1D delta_xs6=xs7-xs1;
      GroupData1D delta_xs7=xs8-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs4.get_dat(i);
        delta4.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs5.get_dat(i);
        delta5.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs6.get_dat(i);
        delta6.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs7.get_dat(i);
        delta7.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs5=e.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs6=f.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs7=g.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs8=h.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      GroupData2D delta_xs4=xs5-xs1;
      GroupData2D delta_xs5=xs6-xs1;
      GroupData2D delta_xs6=xs7-xs1;
      GroupData2D delta_xs7=xs8-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
    	  inp=delta_xs1.get_dat(i,j);
    	  delta1.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs2.get_dat(i,j);
    	  delta2.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs3.get_dat(i,j);
    	  delta3.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs4.get_dat(i,j);
    	  delta4.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs5.get_dat(i,j);
    	  delta5.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs6.get_dat(i,j);
    	  delta6.GetData2d(l,m).put_data(i,j,inp);
    	  inp=delta_xs7.get_dat(i,j);
    	  delta7.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  vector<GroupDataSet> ret(2);
  ret[0]=a;
  ret[1]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      GroupData1D y3=delta4.GetData1d(l,m);
      GroupData1D y4=delta5.GetData1d(l,m);
      GroupData1D y5=delta6.GetData1d(l,m);
      GroupData1D y6=delta7.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        double tmp3=y3.get_dat(j);
        double tmp4=y4.get_dat(j);
        double tmp5=y5.get_dat(j);
        double tmp6=y6.get_dat(j);
        real tmp7[]={tmp,0.,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        LeastSquaresMethodQuadr(8,i,tmp7,tmpa,tmpb,tmpc);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      GroupData2D y3=delta4.GetData2d(l,m);
      GroupData2D y4=delta5.GetData2d(l,m);
      GroupData2D y5=delta6.GetData2d(l,m);
      GroupData2D y6=delta7.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  double tmp2=y2.get_dat(j,k);
      	  double tmp3=y3.get_dat(j,k);
      	  double tmp4=y4.get_dat(j,k);
      	  double tmp5=y5.get_dat(j,k);
      	  double tmp6=y6.get_dat(j,k);
      	  real tmp7[]={tmp,0.,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
      	  LeastSquaresMethodQuadr(8,i,tmp7,tmpa,tmpb,tmpc);
      	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
      	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};
/*
vector<GroupDataSet> GroupDataSet::CaldsdxQuadr(int num, vector<GroupDataSet> &a, real *b)
{
  vector<GroupDataSet> delta(num-1);
  vector<GroupData1D> xs1D(num);
  vector<GroupData2D> xs2D(num);
  vector<GroupData1D> delta_xs1D(num-1);
  vector<GroupData2D> delta_xs2D(num-1);
  for(int i=0;i<num;i++){
    delta[i]=a[i];
  };
  real inp;
  int num1d=a[0].GetNum1d();
  int num2d=a[0].GetNum2d();
  for(int l=0;l<num1d;l++){
    for(int i=0;i<num;i++){
      xs1D[i]=a[i].GetData1d(l,m); //Read for base Macro XSfile1
    };
    for(int i=1;i<num;i++){
      delta_xs1D[i]=xs1D[i]-xs1D[0];
    };
    for(int i=0;i<grp;i++){
      for(int j=0;j<num;j++){
	inp=delta_xs1D[j].get_dat(i);
	delta[j].GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    for(int i=0;i<num;i++){
      xs2D[i]=a[i].GetData2d(l,m); //Read for base Macro XSfile1
    };
    for(int i=1;i<num;i++){
      delta_xs2D[i]=xs2D[i]-xs2D[0];
    };
    for(int i=0;i<grp;i++){
      for(int j=0;j<grp;j++){
	for(int k=0;k<num;k++){
	  inp=delta_xs2D[k].get_dat(i,j);
	  delta[j].GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };  
  vector<GroupDataSet> ret(2);
  vector<GroupData1D> y1D(num-1);
  vector<GroupData2D> y2D(num-1);
  vector<real> y(num);
  ret[0]=a[0];
  ret[1]=a[0];
  for(int l=0;l<num1d;l++){
    for(int i=0;i<num-1;i++){
      y1D[i]=delta[i].GetData1d(l,m); //Read for base Macro XSfile1
    };
    for(int j=0;j<grp;j++){
      y[0]=0.;
      for(int i=0;i<num-1;i++){
	y[i+1]=y1D[i].get_dat(j);
      };
      double tmpa=0.;
      double tmpb=0.;
      double tmpc=0.;
      LeastSquaresMethodQuadr(3,b,y,tmpa,tmpb,tmpc);
      ret[0].GetData1d(l,m).put_data(j,tmpa);
      ret[1].GetData1d(l,m).put_data(j,tmpb);
    };
  };
  for(int l=0;l<num2d;l++){
    for(int i=0;i<num-1;i++){
      y2D[i]=delta[i].GetData2d(l,m); //Read for base Macro XSfile1
    };
    for(int j=0;j<grp;j++){
      for(int k=0;k<grp;k++){
	y[0]=0.;
	for(int i=0;i<num-1;i++){
	  y[i+1]=y2D[i].get_dat(j);
	};
	double tmpa=0.;
	double tmpb=0.;
	double tmpc=0.;
	LeastSquaresMethodQuadr(3,b,y,tmpa,tmpb,tmpc);
	ret[0].GetData2d(l,m).put_data(j,k,tmpa);
	ret[1].GetData2d(l,m).put_data(j,k,tmpb);
      };
    };
  };  
  vector<GroupDataSet> ref(2);
  ref=ret;
  return ref;  
};
*/
vector<GroupDataSet> GroupDataSet::Caldsdx2Cubic(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, real *d)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=a;
  real inp;
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1 
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
  	  inp=delta_xs1.get_dat(i,j);
  	  delta1.GetData2d(l,m).put_data(i,j,inp);
  	  inp=delta_xs2.get_dat(i,j);
  	  delta2.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };  
  vector<GroupDataSet> ret(3);
  ret[0]=a;
  ret[1]=a;
  ret[2]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        real tmp2[]={tmp,0.,tmp1};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        double tmpd=0.;
        LeastSquaresMethodCubic(3,d,tmp2,tmpa,tmpb,tmpc,tmpd);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
        ret[2].GetData1d(l,m).put_data(j,tmpc);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
      	  double tmp=y.get_dat(j,k);
      	  double tmp1=y1.get_dat(j,k);
      	  real tmp2[]={tmp,0.,tmp1};
      	  double tmpa=0.;
      	  double tmpb=0.;
      	  double tmpc=0.;
      	  double tmpd=0.;
      	  LeastSquaresMethodCubic(3,d,tmp2,tmpa,tmpb,tmpc,tmpd);
      	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
      	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
      	  ret[2].GetData2d(l,m).put_data(j,k,tmpc);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(3);
  ref=ret;
  return ref;  
};

vector<GroupDataSet> GroupDataSet::Caldsdx3Cubic(GroupDataSet &a, GroupDataSet &b, GroupDataSet &c, GroupDataSet &d, real *e)
{
  GroupDataSet delta1=a;
  GroupDataSet delta2=b;
  GroupDataSet delta3=c;
  real inp;  
  int num1d=a.GetNum1d();
  int num2d=a.GetNum2d();
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D xs1=a.GetData1d(l,m); //Read for base Macro XSfile1
      GroupData1D xs2=b.GetData1d(l,m); //Read for base Macro XSfile2
      GroupData1D xs3=c.GetData1d(l,m); //Read for base Macro XSfile3
      GroupData1D xs4=d.GetData1d(l,m); //Read for base Macro XSfile4
      GroupData1D delta_xs1=xs2-xs1;
      GroupData1D delta_xs2=xs3-xs1;
      GroupData1D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        inp=delta_xs1.get_dat(i);
        delta1.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs2.get_dat(i);
        delta2.GetData1d(l,m).put_data(i,inp);
        inp=delta_xs3.get_dat(i);
        delta3.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D xs1=a.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs2=b.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs3=c.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs4=d.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D delta_xs1=xs2-xs1;
      GroupData2D delta_xs2=xs3-xs1;
      GroupData2D delta_xs3=xs4-xs1;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
      	  inp=delta_xs1.get_dat(i,j);
      	  delta1.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs2.get_dat(i,j);
      	  delta2.GetData2d(l,m).put_data(i,j,inp);
      	  inp=delta_xs3.get_dat(i,j);
      	  delta3.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  vector<GroupDataSet> ret(3);
  ret[0]=a;
  ret[1]=a;
  ret[2]=a;
  for(int l=0;l<num1d;l++){
    int dim1d=a.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D y=delta1.GetData1d(l,m);
      GroupData1D y1=delta2.GetData1d(l,m);
      GroupData1D y2=delta3.GetData1d(l,m);
      for(int j=0;j<grp;j++){
        double tmp=y.get_dat(j);
        double tmp1=y1.get_dat(j);
        double tmp2=y2.get_dat(j);
        real tmp3[]={tmp,0.,tmp1,tmp2};
        double tmpa=0.;
        double tmpb=0.;
        double tmpc=0.;
        double tmpd=0.;
        LeastSquaresMethodCubic(4,e,tmp3,tmpa,tmpb,tmpc,tmpd);
        ret[0].GetData1d(l,m).put_data(j,tmpa);
        ret[1].GetData1d(l,m).put_data(j,tmpb);
        ret[2].GetData1d(l,m).put_data(j,tmpc);
      };
    };
  };
  for(int l=0;l<num2d;l++){
    int dim2d=a.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D y=delta1.GetData2d(l,m);
      GroupData2D y1=delta2.GetData2d(l,m);
      GroupData2D y2=delta3.GetData2d(l,m);
      for(int j=0;j<grp;j++){
        for(int k=0;k<grp;k++){
  	  double tmp=y.get_dat(j,k);
  	  double tmp1=y1.get_dat(j,k);
  	  double tmp2=y2.get_dat(j,k);
  	  real tmp3[]={tmp,0.,tmp1,tmp2};
  	  double tmpa=0.;
  	  double tmpb=0.;
  	  double tmpc=0.;
  	  double tmpd=0.;
  	  LeastSquaresMethodCubic(4,e,tmp3,tmpa,tmpb,tmpc,tmpd);
  	  ret[0].GetData2d(l,m).put_data(j,k,tmpa);
  	  ret[1].GetData2d(l,m).put_data(j,k,tmpb);
  	  ret[2].GetData2d(l,m).put_data(j,k,tmpc);
	};
      };
    };
  };  
  vector<GroupDataSet> ref(3);
  ref=ret;
  return ref;  
};

GroupDataSet GroupDataSet::CaldeltaXS(GroupDataSet &baseXS, GroupDataSet &coeff, real delta)
{
  real inp;
  GroupDataSet ret=baseXS;
  int d=baseXS.GetNum1d();
  int e=baseXS.GetNum2d();
  for(int l=0;l<d;l++){
    int dim1d=baseXS.GetDim1d(l);
    ret.PutDim1d(l,dim1d);
    for(int m=0;m<dim1d;m++){
      GroupData1D dsdx=coeff.GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D xs=baseXS.GetData1d(l,m);
      GroupData1D delta_xs=dsdx*delta;
      for(int i=0;i<grp;i++){
	inp=xs.get_dat(i)+delta_xs.get_dat(i);
	ret.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<e;l++){
    int dim2d=baseXS.GetDim2d(l);
    ret.PutDim2d(l,dim2d);
    for(int m=0;m<dim2d;m++){
      GroupData2D dsdx=coeff.GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs=baseXS.GetData2d(l,m);
      GroupData2D delta_xs=dsdx*delta;
      for(int i=0;i<grp;i++){
	for(int j=0;j<grp;j++){
	  inp=xs.get_dat(i,j)+delta_xs.get_dat(i,j);
	  ret.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  return ret;
};

GroupDataSet GroupDataSet::CaldeltaXSQuadr(GroupDataSet &baseXS, vector<GroupDataSet> &coeff, real delta)
{
  real inp;
  GroupDataSet ret=baseXS;
  int d=baseXS.GetNum1d();
  int e=baseXS.GetNum2d();
  for(int l=0;l<d;l++){
    int dim1d=baseXS.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D dsdx1=coeff[0].GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D dsdx2=coeff[1].GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D xs=baseXS.GetData1d(l,m);
      GroupData1D delta_xs=dsdx1*delta*delta+dsdx2*delta;
      for(int i=0;i<grp;i++){
        inp=xs.get_dat(i)+delta_xs.get_dat(i);
        ret.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<e;l++){
    int dim2d=baseXS.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D dsdx1=coeff[0].GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D dsdx2=coeff[1].GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs=baseXS.GetData2d(l,m);
      GroupData2D delta_xs=dsdx1*delta*delta+dsdx2*delta;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
  	  inp=xs.get_dat(i,j)+delta_xs.get_dat(i,j);
  	  ret.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  return ret;
};

GroupDataSet GroupDataSet::CaldeltaXSCubic(GroupDataSet &baseXS, vector<GroupDataSet> &coeff, real delta)
{
  real inp;
  GroupDataSet ret=baseXS;
  int d=baseXS.GetNum1d();
  int e=baseXS.GetNum2d();
  for(int l=0;l<d;l++){
    int dim1d=baseXS.GetDim1d(l);
    for(int m=0;m<dim1d;m++){
      GroupData1D dsdx1=coeff[0].GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D dsdx2=coeff[1].GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D dsdx3=coeff[2].GetData1d(l,m); //Read for base Macro XSfile1 of nusigf
      GroupData1D xs=baseXS.GetData1d(l,m);
      GroupData1D delta_xs=dsdx1*delta*delta*delta+dsdx2*delta*delta+dsdx3*delta;
      for(int i=0;i<grp;i++){
        inp=xs.get_dat(i)+delta_xs.get_dat(i);
        ret.GetData1d(l,m).put_data(i,inp);
      };
    };
  };
  for(int l=0;l<e;l++){
    int dim2d=baseXS.GetDim2d(l);
    for(int m=0;m<dim2d;m++){
      GroupData2D dsdx1=coeff[0].GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D dsdx2=coeff[1].GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D dsdx3=coeff[2].GetData2d(l,m); //Read for base Macro XSfile1 of sigs
      GroupData2D xs=baseXS.GetData2d(l,m);
      GroupData2D delta_xs=dsdx1*delta*delta+dsdx2*delta*delta+dsdx3*delta;
      for(int i=0;i<grp;i++){
        for(int j=0;j<grp;j++){
   	  inp=xs.get_dat(i,j)+delta_xs.get_dat(i,j);
   	  ret.GetData2d(l,m).put_data(i,j,inp);
	};
      };
    };
  };
  return ret;
};
