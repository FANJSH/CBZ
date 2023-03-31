#include <cstdlib>
#include "DelayedNeutronData.h"

using namespace std;

DelayedNeutronData::DelayedNeutronData()
{
  nucnum=0;
};

void DelayedNeutronData::PutFamilyNumber(int famin)
{
  famnum=famin;
};

int DelayedNeutronData::GetNumFromNucid(int ii)
{
  for(int i=0;i<nucnum;i++){
    if(nucid[i]==ii)return i;
  };
  //cout<<"# Warning! No nuclide in DelayedNeutronData.\n";
  //cout<<"#  ID is "<<ii<<"\n";
  return -1;
};

void DelayedNeutronData::PutDataFromFile(string mdir,string ss)
{
  mdir.append(ss);

  ifstream fin;
  ofstream fout;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<mdir.data()<<"\n";
    exit(1);
  };

  int ig;
  fin>>ig;

  int inucid;
  fin>>inucid;
  if(inucid<1000){
    inucid=TranslateNuclideIDFromJFS(inucid);
  };
  if(inucid<10000)inucid=midt.GetMATIDFromENDFID(inucid);

  int ii=GetNumFromNucid(inucid);
  if(ii!=-1){
    cout<<"# Error in DelayedNeytronData::PutDataFromFile.\n";
    cout<<"# Delayed neutron data of "<<inucid<<" already exist.\n";
    exit(0);
  };

  nucid.push_back(inucid);

  int famin;
  fin>>famin;

  if(nucnum==0)famnum=famin;

  if(famin!=famnum){
    cout<<"# Error in DelayedNeytronData::PutDataFromFile.\n";
    cout<<"# The number of families is inconsistent.\n";
    cout<<"# Please give up to use the nuclide "<<inucid<<".\n";
    exit(0);
  };

  vector<real> fraction_inp(famin);
  for(int j=0;j<famin;j++){
    fin>>fraction_inp[j];
  };
  fraction.push_back(fraction_inp);

  GroupData1D yield_inp(ig);
  for(int j=0;j<ig;j++){
    real tmp;
    fin>>tmp;
    yield_inp.put_data(j,tmp);
  };
  yield.push_back(yield_inp);

  vector<GroupData1D> chi_inp(famin);
  for(int j=0;j<famin;j++){
    chi_inp[j].put_imax(ig);
    for(int k=0;k<ig;k++){
      real tmp;
      fin>>tmp;
      chi_inp[j].put_data(k,tmp);
    };
  };
  chi.push_back(chi_inp);

  vector<real> lambda_inp(famin);
  for(int j=0;j<famin;j++){
    fin>>lambda_inp[j];
  };
  lambda.push_back(lambda_inp);

  nucnum++;
};

GroupData1D DelayedNeutronData::GetAveragedChi(int inucid)
{
  int ii=GetNumFromNucid(inucid);
  if(ii==-1){
    cout<<"Error in DelayedNeutronData::GetAveragedChi.\n";
    cout<<"Nuclide "<<ii<<" is not included.\n";
    exit(0);
  };
  int ig=chi[ii][0].get_imax();
  GroupData1D avechi(ig);
  avechi.set_zero();
  for(int i=0;i<famnum;i++){
    avechi=avechi+chi[ii][i]*fraction[ii][i];
  };
  return avechi;
};

// ***

void DelayedNeutronData::ShiftChi(GroupData1D &ebnd, real factor)
{
  int grp=ebnd.get_imax()-1;

  real etop_ln=log10(ebnd.get_dat(0));
  real elow_ln=log10(ebnd.get_dat(grp-1));

  for(int i=0;i<nucnum;i++){
    for(int j=0;j<famnum;j++){
      int sz=chi[i][0].get_imax();
      if(sz==grp){
      for(int g=0;g<grp;g++){
	real e_ln=log10(ebnd.get_dat(g));
	real f=(1.-factor)*(etop_ln-e_ln)+(1.+factor)*(e_ln-elow_ln);
	f/=(etop_ln-elow_ln);
	real org=chi[i][j].get_dat(g);
	chi[i][j].put_data(g,org*f);
      };
      real sum=chi[i][j].get_sum();
      sum=1./sum;
      chi[i][j]=chi[i][j]*sum;
    };
    };
  };
};

void DelayedNeutronData::SetConstantYield(int g_set)
{
  int group=yield[0].get_imax();  
  if(g_set<0||g_set>=group){
    cout<<"# Error in DelayedNeutronData::SetConstantYield.\n";
    cout<<"# The target energy group "<<g_set<<" is inappropriate.\n";
    exit(0);
  };
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      yield[i].put_data(g,yield[i].get_dat(g_set));
    };
  };
};

void DelayedNeutronData::ShowSelf()
{
  int group=yield[0].get_imax();
  cout.setf(ios::scientific);
  cout.precision(2);
  cout<<"# Nu_d\n";
  cout<<"# ";
  for(int j=0;j<nucnum;j++){
    cout<<nucid[j]<<" ";
  };
  cout<<"\n";
  for(int g=0;g<group;g++){
    cout<<g<<" ";
    for(int j=0;j<nucnum;j++){
      cout<<yield[j].get_dat(g)<<" ";
    };
    cout<<"\n";
  };

  cout.precision(3);
  cout<<"#\n# Family-wise abundance\n#\n";
  for(int i=0;i<nucnum;i++){
    cout<<nucid[i]<<" ";
    for(int j=0;j<famnum;j++){
      cout<<fraction[i][j]<<" ";
    };
    cout<<"\n";
  };

  cout<<"#\n# Family-wise decay constant\n#\n";
  for(int i=0;i<nucnum;i++){
    cout<<nucid[i]<<" ";
    for(int j=0;j<famnum;j++){
      cout<<lambda[i][j]<<" ";
    };
    cout<<"\n";
  };
};
