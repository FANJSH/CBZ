#include <cstdlib>
#include "GammaLibData.h"

// +++ GammaLibData

GammaLibData::GammaLibData()
{
  xsdata.Init("GammaCrossSection");
};

void GammaLibData::ReadFile(string mdir,string ss,MATIDTranslator &midt)
{
  mdir.append(ss);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"(file name) "<<mdir<<"\n";
    exit(0);
  };

  int tmp;
  real vtmp;

  // control data
  fin>>mat;
  fin>>tmp;
  PutNGroup(tmp);
  if(mat<10000&&tmp!=0)mat=midt.GetMATIDFromENDFID(mat);
  fin>>tmp;
  PutGGroup(tmp);

  int mtn=5;
  int mtg=7;
  int mtg1=3;
  /*
  for(int i=0;i<3;i++){
    fin>>tmp; // dummy for mtn, mtg, mtg1
  };
  */
  fin>>mtn;
  fin>>mtg;
  fin>>mtg1;

  yield.resize(5);
  for(int i=0;i<5;i++){
    yield[i].put_yx(ngroup,ggroup);
    yield[i].set_zero();
  };
  for(int j=0;j<mtn;j++){
    for(int n=0;n<ngroup;n++){
      int kss,kgv;
      fin>>kss;
      fin>>kgv;
      kss--;
      kgv--;
      for(int i=kss;i<=kgv;i++){
	fin>>vtmp;
	yield[j].put_data(n,i,vtmp);
      };
    };
  };

  enum xstype xst[]={sigt,coh_el,incoh_el,pair_pro,siga,kerma};
  for(int j=0;j<mtg;j++){
    for(int i=0;i<ggroup;i++){
      fin>>vtmp;
      if(j!=0){
        xsdata.GetData1d(xst[j-1]).put_data(i,vtmp);
      };
    };
  };

  enum xstype xst2[]={coh_el,incoh_el,pair_pro};
  for(int j=0;j<mtg1;j++){
    fin>>tmp;
    xsdata.PutDim2d(xst2[j],tmp);
    for(int m=0;m<tmp;m++){
      for(int i=0;i<ggroup;i++){
        int lss,lgv;
	fin>>lss;
	fin>>lgv;
	lss--;
	lgv--;
	for(int i1=lss;i1<=lgv;i1++){
	  fin>>vtmp;
	  xsdata.GetData2d(xst2[j],m).put_data(i,i1,vtmp);
	};
      };
    };
  };

  fin.close();
};

void GammaLibData::PutNGroup(int i)
{
  ngroup=i;
};

void GammaLibData::PutGGroup(int i)
{
  ggroup=i;
  xsdata.PutGrp(ggroup);
};

void GammaLibData::ShowSelf()
{
  cout<<"# Material ID              : "<<mat<<"\n";
  cout<<"# Numbers of neutron group : "<<ngroup<<"\n";
  cout<<"# Numbers of gamma group   : "<<ggroup<<"\n"; 
};


// +++ GammaXSLibrary

void GammaXSLibrary::ReadGEnergy(string mdir,string ss)
{
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
    cout<<"Weight function is zero in the [G-ENERGY] file.\n";
    exit(0);
  };
};

GammaLibData& GammaXSLibrary::GetGammaLibData(int mat)
{
  map<int,GammaLibData>::iterator i;
  i=data.find(mat);
  if(i==this->data.end()){
    cout<<"You cannot find GammaLibData in GammaXSLibrary.\n";
    cout<<"MAT number is "<<mat<<"\n";
    exit(0);
  }else{
    return i->second;
  };
};

bool GammaXSLibrary::ExistGammaLibData(int mat)
{
  map<int,GammaLibData>::iterator i;
  i=data.find(mat);
  if(i==this->data.end()){
    return false;
  }else{
    return true;
  };
};

void GammaXSLibrary::ReadFile(int nucnum,string mdir,string *filename)
{
  GammaLibData tmp;
  for(int i=0;i<nucnum;i++){
    tmp.ReadFile(mdir,filename[i],midt);
    AddGammaLibData(tmp.GetMat(),tmp);
  };
};

void GammaXSLibrary::ReadFile(int nucnum,string mdir,string *filename,int *matid)
{
  GammaLibData tmp;
  for(int i=0;i<nucnum;i++){
    tmp.ReadFile(mdir,filename[i],midt);
    AddGammaLibData(matid[i],tmp);
  };
};

void GammaXSLibrary::AddDelayedFissionGamma(string mdir,string fname,real fraction)
{
  mdir.append(fname);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"(file name) "<<fname<<"\n";
    exit(0);
  };

  int grp;
  fin>>grp;
  int nuc;
  fin>>nuc;

  for(int i=0;i<nuc;i++){
    int mat;
    fin>>mat;
    if(mat<10000)mat=midt.GetMATIDFromENDFID(mat);
    bool ext=ExistGammaLibData(mat);
    int ngroup=0;
    if(!ext){
      cout<<"#  No gamma data : "<<mat<<"\n"; 
    }else{
      ngroup=GetGammaLibData(mat).GetNGroup();
    };
    for(int j=0;j<grp;j++){
      real val;
      fin>>val;
      if(ext){
	for(int k=0;k<ngroup;k++){
	  GetGammaLibData(mat).GetYield(1).add_data(k,j,val*fraction);
	};
      };
    };
  };
};
