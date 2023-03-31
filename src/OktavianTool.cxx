#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "OktavianTool.h"

using namespace std;

OktavianTool::OktavianTool(int ginp)
{
  PutGroup(ginp);
};

void OktavianTool::PutGroup(int ginp)
{
  grp=ginp;
};

void OktavianTool::ReadSensFile(string mdir,string filename)
{
  ifstream fin;
  mdir.append(filename);
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the sensitivity file.\n";
    cout<<"File name is "<<mdir<<"\n";
    exit(0);
  };

  int tmp; 
  fin>>tmp; // number of groups

  while(!fin.eof()){

  fin>>tmp; // (dimension)
  if(tmp!=1&&tmp!=2){
    cout<<"Error in OktavianTool::ReadSensFile,\n";
    cout<<"Dimension "<<tmp<<" is inappropriate.\n";
    exit(0);
  };

  int mat,mt;
  fin>>mat; // (mat_id)
  fin>>mt; // (mt_id)

  if(tmp==1){ // 1D data
    sns1d_mat.push_back(mat);
    sns1d_mt.push_back(mt);
    GroupData1D tmp1d;
    tmp1d.put_imax(grp);
    for(int i=0;i<grp;i++){
      real tmp;
      fin>>tmp;
      tmp1d.put_data(i,tmp);
    };
    sns1d.push_back(tmp1d);
  }else{ // 2D data
    sns2d_mat.push_back(mat);
    sns2d_mt.push_back(mt);
    GroupData2D tmp2d(grp,grp);
    tmp2d.set_zero();
    if(mt==2){ // P0 elastic
      for(int i=0;i<grp;i++){
        int maxg=i+10;
        if(maxg>grp)maxg=grp;
        for(int j=i;j<maxg;j++){
          real tmp;
          fin>>tmp;
          tmp2d.put_data(i,j,tmp);
        };
      };
    }else if(mt==4){ // inelastic
      fin>>tmp; 
      for(int i=0;i<tmp;i++){
        for(int j=i;j<grp;j++){
          real tmp;
          fin>>tmp;
          tmp2d.put_data(i,j,tmp);
	};
      };
    }else if(mt==16){ // n2n
      fin>>tmp;
      for(int i=0;i<tmp;i++){
        for(int j=i;j<grp;j++){
          real tmp;
          fin>>tmp;
          tmp2d.put_data(i,j,tmp);
        };
      };
    }else{
      cout<<"Error in OktavianTool::ReadSensFile.\n";
      cout<<"Not coding for 2D sensitivity, mt="<<mt<<"\n";
      exit(0);
    };
    sns2d.push_back(tmp2d);
  };

  };

  fin.close();
};

void OktavianTool::ShowSensData()
{
  int size1d=sns1d.size();
  cout<<"+++ 1D sensitivity +++\n";
  for(int i=0;i<size1d;i++){
    cout<<"   "<<sns1d_mat[i]<<" "<<sns1d_mt[i]<<"\n";
  };
  int size2d=sns2d.size();
  cout<<"+++ 2D sensitivity +++\n";
  for(int i=0;i<size2d;i++){
    cout<<"   "<<sns2d_mat[i]<<" "<<sns2d_mt[i]<<"\n";
  };
};

GroupData1D& OktavianTool::GetSens1D(int mati,int mti)
{
  int size1d=sns1d.size();
  for(int i=0;i<size1d;i++){
    if(mati==sns1d_mat[i]&&mti==sns1d_mt[i])return sns1d[i];
  };
  cout<<"Error in OktavianTool::GetSens1D.\n";
  cout<<"No 1D sensitivity data (mat/mt="<<mati<<"/"<<mti<<")\n";
};

GroupData2D& OktavianTool::GetSens2D(int mati,int mti)
{
  int size2d=sns2d.size();
  for(int i=0;i<size2d;i++){
    if(mati==sns2d_mat[i]&&mti==sns2d_mt[i])return sns2d[i];
  };
  cout<<"Error in OktavianTool::GetSens2D.\n";
  cout<<"No 2D sensitivity data (mat/mt="<<mati<<"/"<<mti<<")\n";
};

void OktavianTool::LibraryEffectCal(XSLibrary &lib1,XSLibrary &lib2,int mati,bool print_eng)
{
  GroupData1D sns1d;
  GroupData2D sns2d;

  if(print_eng)cout<<"\n# Capture\n";
  sns1d=sens.GetSensitivity1D(mati,102);
  real sum=0.;
  for(int i=0;i<grp;i++){
    real xs1=lib1.GetLibData(mati).GetXSData().GetData1d(sigc).get_dat(i);
    real xs2=lib2.GetLibData(mati).GetXSData().GetData1d(sigc).get_dat(i);
    real tmp=(xs2-xs1)/xs1*sns1d.get_dat(i);
    sum+=tmp;
    if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmp<<"\n";
  };
  real capc=sum;

  if(print_eng)cout<<"\n# Mu-bar\n";
  sns1d=sens.GetSensitivity1D(mati,251);
  sum=0.;
  for(int i=0;i<grp;i++){
    real xs1=lib1.GetLibData(mati).GetXSData().GetData1d(mu).get_dat(i);
    real xs2=lib2.GetLibData(mati).GetXSData().GetData1d(mu).get_dat(i);
    real tmp=(xs2-xs1)*sns1d.get_dat(i);
    sum+=tmp;
    if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmp<<"\n";
  };
  real muc=sum;

  if(print_eng)cout<<"\n# Elastic\n";
  sns2d=sens.GetSensitivity2D(mati,2);
  sum=0.;
  for(int i=0;i<grp;i++){
    real tmpsum=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=lib1.GetLibData(mati).GetXSData().GetData2d(sigel).get_dat(i,j);
      real xs2=lib2.GetLibData(mati).GetXSData().GetData2d(sigel).get_dat(i,j);
      if(xs1>0.){
        real tmp=(xs2-xs1)/xs1*sns2d.get_dat(i,j);
        sum+=tmp;
        tmpsum+=tmp;
      };
    };
    if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmpsum<<"\n";
  };
  real el0c=sum;

  if(print_eng)cout<<"\n# Inelastic\n";
  sns2d=sens.GetSensitivity2D(mati,4);
  sum=0.;
  for(int i=0;i<grp;i++){
    real tmpsum=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=lib1.GetLibData(mati).GetXSData().GetData2d(siginel).get_dat(i,j);
      real xs2=lib2.GetLibData(mati).GetXSData().GetData2d(siginel).get_dat(i,j);
      real tmp=(xs2-xs1)*sns2d.get_dat(i,j);
      sum+=tmp;
      tmpsum+=tmp;
    };
    if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmpsum<<"\n";
  };
  real iec=sum;

  if(print_eng)cout<<"\n# (n,2n)\n";
  sns2d=sens.GetSensitivity2D(mati,16);
  sum=0.;
  for(int i=0;i<grp;i++){
    real tmpsum=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=lib1.GetLibData(mati).GetXSData().GetData2d(sign2n).get_dat(i,j);
      real xs2=lib2.GetLibData(mati).GetXSData().GetData2d(sign2n).get_dat(i,j);
      real tmp=(xs2-xs1)*sns2d.get_dat(i,j);
      sum+=tmp;
      tmpsum+=tmp;
    };
    if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmpsum<<"\n";
  };
  real n2nc=sum;


  cout<<"  (Total)   : "<<capc+el0c+iec+n2nc+muc<<"\n";
  cout<<"  Capture   : "<<capc<<"\n";
  cout<<"  Elastic   : "<<el0c<<"\n";
  cout<<"  Inelastic : "<<iec<<"\n";
  cout<<"  N2N       : "<<n2nc<<"\n";
  cout<<"  Mu_bar    : "<<muc<<"\n";
};

void OktavianTool::LibraryEffectCal(Medium &med1,Medium &med2,int mati,bool print_eng)
{
  GroupData1D sns1d;
  GroupData2D sns2d;

  vector<real> letwid(grp);
  for(int i=0;i<grp;i++){
    real e0=med1.GetEnband().get_dat(i);
    real e1=med1.GetEnband().get_dat(i+1);
    letwid[i]=log(e0/e1);
  };

  if(print_eng)cout<<"\n\n# Capture\n";
  sns1d=sens.GetSensitivity1D(mati,102);
  real sum=0.;
  for(int i=0;i<grp;i++){
    real xs1=med1.GetNuclide(mati).GetMicxs().GetData1d(sigc).get_dat(i);
    real xs2=med2.GetNuclide(mati).GetMicxs().GetData1d(sigc).get_dat(i);
    real tmp=(xs2-xs1)/xs1*sns1d.get_dat(i);
    sum+=tmp;
    if(print_eng)cout<<med1.GetEnband().get_dat(i)<<" "<<tmp/letwid[i]<<"\n";
  };
  real capc=sum;

  if(print_eng)cout<<"\n\n# Mu-bar\n";
  sns1d=sens.GetSensitivity1D(mati,251);
  sum=0.;
  for(int i=0;i<grp;i++){
    real xs1=med1.GetNuclide(mati).GetMicxs().GetData1d(mu).get_dat(i);
    real xs2=med2.GetNuclide(mati).GetMicxs().GetData1d(mu).get_dat(i);
    real tmp=(xs2-xs1)*sns1d.get_dat(i);
    sum+=tmp;
    if(print_eng)cout<<med1.GetEnband().get_dat(i)<<" "<<tmp/letwid[i]<<"\n";
  };
  real muc=sum;

  if(print_eng)cout<<"\n\n# Elastic\n";
  sns2d=sens.GetSensitivity2D(mati,2);
  sum=0.;
  for(int i=0;i<grp;i++){
    real sum2=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=med1.GetNuclide(mati).GetMicxs().GetData2d(sigel).get_dat(i,j);
      real xs2=med2.GetNuclide(mati).GetMicxs().GetData2d(sigel).get_dat(i,j);
      if(xs1>0.){
        real tmp=(xs2-xs1)/xs1*sns2d.get_dat(i,j);
        sum+=tmp;
        sum2+=tmp;
      };
    };
    if(print_eng)cout<<med1.GetEnband().get_dat(i)<<" "<<sum2/letwid[i]<<"\n";
  };
  real el0c=sum;

  if(print_eng)cout<<"\n\n# Inelastic\n";
  sns2d=sens.GetSensitivity2D(mati,4);
  sum=0.;
  for(int i=0;i<grp;i++){
    real sum2=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=med1.GetNuclide(mati).GetMicxs().GetData2d(siginel).get_dat(i,j);
      real xs2=med2.GetNuclide(mati).GetMicxs().GetData2d(siginel).get_dat(i,j);
      real tmp=(xs2-xs1)*sns2d.get_dat(i,j);
      sum+=tmp;
      sum2+=tmp;
    };
    if(print_eng)cout<<med1.GetEnband().get_dat(i)<<" "<<sum2/letwid[i]<<"\n";
  };
  real iec=sum;

  if(print_eng)cout<<"\n\n# (n,2n)\n";
  sns2d=sens.GetSensitivity2D(mati,16);
  sum=0.;
  for(int i=0;i<grp;i++){
    real sum2=0.;
    for(int j=i;j<grp;j++){ 
      real xs1=med1.GetNuclide(mati).GetMicxs().GetData2d(sign2n).get_dat(i,j);
      real xs2=med2.GetNuclide(mati).GetMicxs().GetData2d(sign2n).get_dat(i,j);
      real tmp=(xs2-xs1)*sns2d.get_dat(i,j);
      sum+=tmp;
      sum2+=tmp;
    };
    if(print_eng)cout<<med1.GetEnband().get_dat(i)<<" "<<sum2/letwid[i]<<"\n";
  };
  real n2nc=sum;

  cout.setf(ios::showpoint);
  cout.precision(4);
  cout<<"#  (Total)   : "<<capc+el0c+iec+n2nc+muc<<"\n";
  cout<<"#  Capture   : "<<capc<<"\n";
  cout<<"#  Elastic   : "<<el0c<<"\n";
  cout<<"#  Inelastic : "<<iec<<"\n";
  cout<<"#  N2N       : "<<n2nc<<"\n";
  cout<<"#  Mu_bar    : "<<muc<<"\n";
};

void OktavianTool::LibraryEffectCalInelastic
(string lib1,string lib2,string fname,int mati)
{
  string name_reaction[]={
   "n01","n02","n03","n04","n05", "n06","n07","n08","n09","n10",
   "n11","n12","n13","n14","n15", "n16","n17","n18","n19","n20",
   "n21","n22","n23","n24","n25", "n26","n27","n28","n29","n30",
   "n31","n32","n33","n34","n35", "n36","n37","n38","n39","n40",
   "n41","n42","n43","n44","n45", "n46","n47","n48","n49","n50",
   "n51","n52","n53","n54","n55", "n56","n57","n58","n59","ncn",
   "nna","nnp","nnd","nnt"
  };
  int inmie=64;

  lib1.append(fname);
  lib2.append(fname);
  lib1.append(".ine");
  lib2.append(".ine");
    
  vector<GroupData2D> ine1(inmie);
  vector<GroupData2D> ine2(inmie);
  for(int i=0;i<inmie;i++){
    ine1[i].put_yx(grp,grp);
    ine2[i].put_yx(grp,grp);
  };

  ifstream fin1;
  fin1.open(lib1.data(),ios::in);
  if(fin1.fail()){
    cout<<"Failed to open the inelastic matrix file.\n";
    cout<<"File name is "<<lib1.data()<<"\n";
    exit(0);
  };
  
  while(!fin1.eof()){
    int id;
    fin1>>id;
    if(id==-1)break;
    for(int g=0;g<grp;g++){
      int i1,i2;
      fin1>>i1;
      fin1>>i2;
      for(int g2=i1;g2>=i2;g2--){
        real tmp;
        fin1>>tmp;
	ine1[id-1].put_data(g2-1,g,tmp);
      };
    };
  };

  ifstream fin2;
  fin2.open(lib2.data(),ios::in);
  if(fin2.fail()){
    cout<<"Failed to open the inelastic matrix file.\n";
    cout<<"File name is "<<lib2.data()<<"\n";
    exit(0);
  };
  
  while(!fin2.eof()){
    int id;
    fin2>>id;
    if(id==-1)break;
    for(int g=0;g<grp;g++){
      int i1,i2;
      fin2>>i1;
      fin2>>i2;
      for(int g2=i1;g2>=i2;g2--){
        real tmp;
        fin2>>tmp;
	ine2[id-1].put_data(g2-1,g,tmp);
      };
    };
  };

  GroupData2D sns2d=GetSens2D(mati,4);

  real totsum=0.;
  for(int ii=0;ii<inmie;ii++){
    real sum=0.;
    for(int i=0;i<grp;i++){
      real tmpsum=0.;
      for(int j=i;j<grp;j++){ 
        real xs1=ine1[ii].get_dat(i,j);
        real xs2=ine2[ii].get_dat(i,j);
        real tmp=(xs2-xs1)*sns2d.get_dat(i,j);
        sum+=tmp;
        tmpsum+=tmp;
      };
      //if(print_eng)cout<<lib1.GetEnband().get_dat(i)<<" "<<tmpsum<<"\n";
    };
    if(sum!=0.)cout<<name_reaction[ii]<<" : "<<sum<<"\n";
    totsum+=sum;
  };
  cout<<"\nTotal : "<<totsum<<"\n";
};

void OktavianTool::ShowSecondaryDistribution(string lib1,string fname,int srcg,GroupData1D &ebnd)
{
  string name_reaction[]={
   "n01","n02","n03","n04","n05", "n06","n07","n08","n09","n10",
   "n11","n12","n13","n14","n15", "n16","n17","n18","n19","n20",
   "n21","n22","n23","n24","n25", "n26","n27","n28","n29","n30",
   "n31","n32","n33","n34","n35", "n36","n37","n38","n39","n40",
   "n41","n42","n43","n44","n45", "n46","n47","n48","n49","n50",
   "n51","n52","n53","n54","n55", "n56","n57","n58","n59","ncn",
   "nna","nnp","nnd","nnt"
  };
  int inmie=64;

  lib1.append(fname);
  lib1.append(".ine");
    
  vector<GroupData2D> ine1(inmie);
  for(int i=0;i<inmie;i++){
    ine1[i].put_yx(grp,grp);
  };

  ifstream fin1;
  fin1.open(lib1.data(),ios::in);
  if(fin1.fail()){
    cout<<"Failed to open the inelastic matrix file.\n";
    cout<<"File name is "<<lib1.data()<<"\n";
    exit(0);
  };
  
  int index=0;
  while(!fin1.eof()){
    int id;
    fin1>>id;
    if(id==-1)break;
    for(int g=0;g<grp;g++){
      int i1,i2;
      fin1>>i1;
      fin1>>i2;
      for(int g2=i1;g2>=i2;g2--){
        real tmp;
        fin1>>tmp;
	ine1[id-1].put_data(g2-1,g,tmp);
      };
    };
    cout<<"# "<<name_reaction[id-1]<<" ("<<index<<")\n";
    for(int g=0;g<grp;g++){
      cout<<ebnd.get_dat(g)<<" "<<ine1[id-1].get_dat(srcg,g)<<"\n";
    };
    cout<<"\n\n";
    index++;
  };
};

