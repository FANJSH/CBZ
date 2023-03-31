#include <cstdlib>
#include "UNC_Parameters.h"

// +++++++++++++
// ParameterList
// +++++++++++++

ParameterList::ParameterList()
{
};

int ParameterList::GetSize()
{
  int size_list=core_tag.size();
  return size_list;
};

int ParameterList::FindData(string core,string chara,int step)
{
  int size_list=core_tag.size();
  for(int i=0;i<size_list;i++){
    if(core_tag[i]==core&&chara_tag[i]==chara&&step_tag[i]==step)return i;
  };
  return -1;
};

void ParameterList::ShowSelf()
{
  int size_list=core_tag.size();
  if(size_list==0){
    cout<<"There is no data in ParameterList.\n";
  };
  for(int i=0;i<size_list;i++){
    cout<<i<<" : "<<core_tag[i]<<" "<<chara_tag[i]<<" "<<step_tag[i]<<"\n";
  };
};

void ParameterList::AddNewData(string core,string chara,int step)
{
  if(FindData(core,chara,step)!=-1){
    cout<<"The following data already exists in ParameterList.\n";
    cout<<"  Core            : "<<core<<"\n";
    cout<<"  Characteristics : "<<chara<<"\n";
    cout<<"  Step            : "<<step<<"\n";
  };
  core_tag.push_back(core);
  chara_tag.push_back(chara);
  step_tag.push_back(step);
};

ParameterList ParameterList::operator+(ParameterList &other_list)
{
  ParameterList ret;
  int other_size=other_list.GetSize();
  int size=GetSize();

  for(int i=0;i<size;i++){
    ret.AddNewData(core_tag[i],chara_tag[i],step_tag[i]);
  };

  for(int i=0;i<other_size;i++){
    string other_core=other_list.GetCoreTag(i);
    string other_chara=other_list.GetCharaTag(i);
    int other_step=other_list.GetStepTag(i);
    if(FindData(other_core,other_chara,other_step)==-1){
      ret.AddNewData(other_core,other_chara,other_step);
    };
  };

  return ret;
};

// ++++++++++
// Parameters 
// ++++++++++

Parameters::Parameters(ParameterList *inp)
{
  plist=inp;
  int size=plist->GetSize();
  val.put_imax(size);
};

void Parameters::PutValue(string core,string chara,real value,int step)
{
  int pos=plist->FindData(core,chara,step);
  if(pos==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core<<"\n";
    cout<<"  Characteristics : "<<chara<<"\n";
    cout<<"  Step            : "<<step<<"\n";
    exit(0);
  };
  val.put_data(pos,value);
};

void Parameters::ShowSelf()
{
  int size=plist->GetSize();
  for(int i=0;i<size;i++){
    cout<<i<<": ["<<plist->GetCoreTag(i)<<"]";
    cout<<"["<<plist->GetCharaTag(i)<<"]";
    cout<<"["<<plist->GetStepTag(i)<<"] : ";
    cout<<val.get_dat(i)<<"\n";
  };
};

// +++++++++++++++++++
// ParameterCovariance
// +++++++++++++++++++

ParameterCovariance::ParameterCovariance(Parameters& inp):Covariance()
{
  plist=inp.GetParameterList();
  int size=plist->GetSize();
  PutType("variance");
  PutGrp(size);
  for(int i=0;i<size;i++){
    PutValue(plist->GetCoreTag(i),plist->GetCharaTag(i),inp.GetValue(i),plist->GetStepTag(i));
  };
};

void ParameterCovariance::PutValue(string core,string chara,real value,int step)
{
  int pos=plist->FindData(core,chara,step);
  if(pos==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core<<"\n";
    cout<<"  Characteristics : "<<chara<<"\n";
    cout<<"  Step            : "<<step<<"\n";
    exit(0);
  };
  PutVal(pos,value);
};

void ParameterCovariance::PutStandardDeviation(string core,string chara,real value,int step,string type)
{
  int pos=plist->FindData(core,chara,step);
  if(pos==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core<<"\n";
    cout<<"  Characteristics : "<<chara<<"\n";
    cout<<"  Step            : "<<step<<"\n";
    exit(0);
  };
  Covariance::PutStandardDeviation(pos,value,type);
};

real ParameterCovariance::GetStandardDeviation(string core,string chara,int step,string type)
{
  int pos=plist->FindData(core,chara,step);
  if(pos==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core<<"\n";
    cout<<"  Characteristics : "<<chara<<"\n";
    cout<<"  Step            : "<<step<<"\n";
    exit(0);
  };
  return Covariance::GetStandardDeviation(type).get_dat(pos);
};

void ParameterCovariance::PutCorrelation(real cor)
{
  int size=plist->GetSize();
  for(int i=0;i<size;i++){
    for(int j=i+1;j<size;j++){
      real dev1=sqrt(cov.get_dat(i,i));
      real dev2=sqrt(cov.get_dat(j,j));
      real tmp=dev1*dev2*cor;
      cov.put_data(i,j,tmp);
      cov.put_data(j,i,tmp);
    };
  };
};

void ParameterCovariance::PutCorrelation
(string core1,string chara1,string core2,string chara2,real cor)
{
  int size=plist->GetSize();

  if(size==0){
    cout<<"You cannot do PutCorrelation\n";
    cout<<"since there is no data.\n";
    exit(0);
  };

  for(int i=0;i<size;i++){
    if(plist->GetCoreTag(i)==core1||core1==""){
      if(chara1==""||chara1==plist->GetCharaTag(i)){
	for(int j=0;j<size;j++){
	  if(plist->GetCoreTag(j)==core2||core2==""){
	    if(chara2==""||chara2==plist->GetCharaTag(j)){
	      if(i!=j){
		real dev1=sqrt(cov.get_dat(i,i));
		real dev2=sqrt(cov.get_dat(j,j));
		real tmp=dev1*dev2*cor;
		cov.put_data(i,j,tmp);
		cov.put_data(j,i,tmp);
	      };
	    };
	  };
	};
      };
    };
  };
};

void ParameterCovariance::PutCorrelation(string core,string chara,real cor)
{
  int size=plist->GetSize();

  vector<bool> coind(size,false);
  bool whole_coind=false;
  for(int i=0;i<size;i++){
    if(plist->GetCoreTag(i)==core&&plist->GetCharaTag(i)==chara){
      coind[i]=true;
      whole_coind=true;
    };
  };

  if(!whole_coind){
    cout<<"# Error in ParameterCovariance::PutCorrelation.\n";
    cout<<"# No perameters : "<<core<<" / "<<chara<<"\n";
    exit(0);
  };

  for(int i=0;i<size;i++){
    for(int j=i+1;j<size;j++){
      if(coind[i]&&coind[j]){
        real dev1=sqrt(cov.get_dat(i,i));
	real dev2=sqrt(cov.get_dat(j,j));
	real tmp=dev1*dev2*cor;
	cov.put_data(i,j,tmp);
	cov.put_data(j,i,tmp);
      };
    };
  };
};

real ParameterCovariance::GetCorrelation
(string core1,string chara1,string core2,string chara2,int step1,int step2)
{
  int pos1=plist->FindData(core1,chara1,step1);
  int pos2=plist->FindData(core2,chara2,step2);
  if(pos1==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core1<<"\n";
    cout<<"  Characteristics : "<<chara1<<"\n";
    cout<<"  Step            : "<<step1<<"\n";
    exit(0);
  };
  if(pos2==-1){
    cout<<"We cannot find the following data in ParameterList.\n";
    cout<<"  Core            : "<<core2<<"\n";
    cout<<"  Characteristics : "<<chara2<<"\n";
    cout<<"  Step            : "<<step2<<"\n";
    exit(0);
  };
  return GetCorrelationMatrix().get_dat(pos1,pos2);
};

// +++++++++++++++++++
// ParametersContainer
// +++++++++++++++++++

int ParametersContainer::FindData(string nameinp)
{
  int size=name.size();
  int tmp=-1;
  for(int i=0;i<size;i++){
    if(nameinp==name[i])tmp=i;
  };
  return tmp;
};

void ParametersContainer::PutParameters(Parameters inp,string nameinp)
{
  int tmp=FindData(nameinp);
  if(tmp==-1){
    paras.push_back(inp);
    name.push_back(nameinp);
  }else{
    paras[tmp]=inp;
  };
};

Parameters &ParametersContainer::GetParameters(string nameinp)
{
  int tmp=FindData(nameinp);
  if(tmp==-1){
    cout<<"You cannot find Parameters "<<nameinp<<"\n";
    exit(0);
  };
  return paras[tmp];
};
