#include <cstdlib>
#include "SensitivityData.h"

SensitivityData::SensitivityData()
{
  group=-1;
};

void SensitivityData::PutGroup(int i)
{
  group=i;
  enband.put_imax(group+1);
};

void SensitivityData::PutEnergyBoundaryData(int i,real j)
{
  if(group==-1){
    cout<<"# Error in SensitivityData::PutEnergyBoundaryData.\n";
    cout<<"# Energy group has NOT been assigned.\n";
    exit(0);
  };

  if(i<0||i>group){
    cout<<"# Error in SensitivityData::PutEnergyBoundaryData.\n";
    cout<<"# You choose inappropriate energy group number "<<i<<"\n";
    exit(0);
  };

  enband.put_data(i,j);
};

int SensitivityData::FindData0D(int mat,int mt)
{
  int size=sens0d.size();
  for(int i=0;i<size;i++){
    if(mat_list0d[i]==mat&&mt_list0d[i]==mt)return i;
  };
  return -1;
};

int SensitivityData::FindData1D(int mat,int mt)
{
  int size=sens1d.size();
  for(int i=0;i<size;i++){
    if(mat_list1d[i]==mat&&mt_list1d[i]==mt)return i;
  };
  return -1;
};

int SensitivityData::FindData2D(int mat,int mt)
{
  int size=sens2d.size();
  for(int i=0;i<size;i++){
    if(mat_list2d[i]==mat&&mt_list2d[i]==mt)return i;
  };
  return -1;
};

real &SensitivityData::GetSensitivity0D(int mat,int mt)
{
  int tmp=FindData0D(mat,mt);
  if(tmp==-1){
    cout<<"# Error in SensitivityData::GetSensitivity0D.\n";
    cout<<"# Sensitivity data does NOT exist for MAT/MT="<<mat<<"/"<<mt<<"\n";
    exit(0);
  };
  return sens0d[tmp];
};

GroupData1D &SensitivityData::GetSensitivity1D(int mat,int mt)
{
  int tmp=FindData1D(mat,mt);
  if(tmp==-1){
    cout<<"# Error in SensitivityData::GetSensitivity1D.\n";
    cout<<"# Sensitivity data does NOT exist for MAT/MT="<<mat<<"/"<<mt<<"\n";
    exit(0);
  };
  return sens1d[tmp];
};

GroupData2D &SensitivityData::GetSensitivity2D(int mat,int mt)
{
  int tmp=FindData2D(mat,mt);
  if(tmp==-1){
    cout<<"# Error in SensitivityData::GetSensitivity2D.\n";
    cout<<"# Sensitivity data does NOT exist for MAT/MT="<<mat<<"/"<<mt<<"\n";
    exit(0);
  };
  return sens2d[tmp];
};

void SensitivityData::PutSensitivity0D(int mat,int mt,real dat)
{
  int pos=FindData0D(mat,mt);
  if(pos==-1){
    sens0d.push_back(dat);
    mat_list0d.push_back(mat);
    mt_list0d.push_back(mt);
  }else{
    //cout<<"# Note in SensitivityData::PutSensitivity0D.\n";
    //cout<<"# Sensitivity is added : "<<mat<<" "<<mt<<"\n";
    sens0d[pos]+=dat;
  };
};

void SensitivityData::PutSensitivity1D(int mat,int mt,GroupData1D dat)
{
  int grp=dat.get_imax();
  if(grp!=group){
    cout<<"# Error in SensitivityData::PutSensitivity1D.\n";
    cout<<"# Energy group is inconsistent.\n";
    exit(0);
  };

  int pos=FindData1D(mat,mt);
  if(pos==-1){
    sens1d.push_back(dat);
    mat_list1d.push_back(mat);
    mt_list1d.push_back(mt);
  }else{
    /*
    bool non_zero=false;
    for(int g=0;g<group;g++){
      if(sens1d[pos].get_dat(g)>0.)non_zero=true;
    };
    if(non_zero){
    */
    //cout<<"# Note in SensitivityData::PutSensitivity1D.\n";
    //cout<<"# Sensitivity is added : "<<mat<<" "<<mt<<"\n";
      //};
    sens1d[pos]=sens1d[pos]+dat;
  };
};

void SensitivityData::PutSensitivity2D(int mat,int mt,GroupData2D dat)
{
  int grp=dat.get_x();
  if(grp!=group){
    cout<<"# Error in SensitivityData::PutSensitivity2D.\n";
    cout<<"# Energy group is inconsistent.\n";
    exit(0);
  };

  int pos=FindData2D(mat,mt);
  if(pos==-1){
    sens2d.push_back(dat);
    mat_list2d.push_back(mat);
    mt_list2d.push_back(mt);
  }else{
    //cout<<"# Note in SensitivityData::PutSensitivity2D.\n";
    //cout<<"# Sensitivity is added : "<<mat<<" "<<mt<<"\n";
    sens2d[pos]=sens2d[pos]+dat;
  };
};

void SensitivityData::AddSensitivityData(SensitivityData &sns2,bool replace)
{
  for(int i=0;i<sns2.GetSize0D();i++){
    int mat=sns2.GetMatList0D(i);
    int mt=sns2.GetMtList0D(i);
    int id=FindData0D(mat,mt);
    if(id==-1){
      PutSensitivity0D(mat,mt,sns2.GetSensitivity0D(mat,mt));
    }else{
      real org=GetSensitivity0D(mat,mt);
      if(replace){
        GetSensitivity0D(mat,mt)=sns2.GetSensitivity0D(mat,mt);
      }else{
        GetSensitivity0D(mat,mt)=org+sns2.GetSensitivity0D(mat,mt);
      };
    };
  };

  for(int i=0;i<sns2.GetSize1D();i++){
    int mat=sns2.GetMatList1D(i);
    int mt=sns2.GetMtList1D(i);
    int id=FindData1D(mat,mt);
    if(id==-1){
      PutSensitivity1D(mat,mt,sns2.GetSensitivity1D(mat,mt));
    }else{
      GroupData1D org=GetSensitivity1D(mat,mt);
      if(replace){
        GetSensitivity1D(mat,mt)=sns2.GetSensitivity1D(mat,mt);
      }else{
        GetSensitivity1D(mat,mt)=org+sns2.GetSensitivity1D(mat,mt);
      };
    };
  };

  for(int i=0;i<sns2.GetSize2D();i++){
    int mat=sns2.GetMatList2D(i);
    int mt=sns2.GetMtList2D(i);
    int id=FindData2D(mat,mt);
    if(id==-1){
      PutSensitivity2D(mat,mt,sns2.GetSensitivity2D(mat,mt));
    }else{
      GroupData2D org=GetSensitivity2D(mat,mt);
      if(replace){
        GetSensitivity2D(mat,mt)=sns2.GetSensitivity2D(mat,mt);
      }else{
        GetSensitivity2D(mat,mt)=org+sns2.GetSensitivity2D(mat,mt);
      };
    };
  };

};

void SensitivityData::AverageSensitivityData(vector<SensitivityData> &sns)//kawamoto
{
  int num=sns.size();
  vector<vector<real> > org0D;
  vector<vector<GroupData1D> > org1D;
  vector<vector<GroupData2D> >org2D;
  vector<real> val_array;
  org0D.resize(num);
  org1D.resize(num);
  org2D.resize(num);
  val_array.resize(num);
  real tmp0D;
  GroupData1D tmp1D;
  GroupData2D tmp2D;

  //+++Set parameter+++
  PutName(sns[0].GetName1(),sns[0].GetName2(),sns[0].GetName3());
  PutGroup(sns[0].GetGroup());
  GroupData1D tmpe=sns[0].GetEnband();
  enband=tmpe;
  int size0D=sns[0].GetSize0D();
  int size1D=sns[0].GetSize1D();
  int size2D=sns[0].GetSize2D();
  for(int i=0;i<size0D;i++){
    mat_list0d.push_back(sns[0].GetMatList0D(i));
    mt_list0d.push_back(sns[0].GetMtList0D(i));
  };
  for(int i=0;i<size1D;i++){
    mat_list1d.push_back(sns[0].GetMatList1D(i));
    mt_list1d.push_back(sns[0].GetMtList1D(i));
  };
  for(int i=0;i<size2D;i++){
    mat_list2d.push_back(sns[0].GetMatList2D(i));
    mt_list2d.push_back(sns[0].GetMtList2D(i));
  };
  sens0d.resize(size0D);
  sens1d.resize(size1D);
  sens2d.resize(size2D);
  //+++++++++++++++++++

  for(int j=0;j<num;j++){  
    val_array[j]=sns[j].GetValue();
    org0D[j].resize(sns[j].GetSize0D());
    for(int i=0;i<sns[j].GetSize0D();i++){
      tmp0D=sns[j].GetSensitivity0D(i);
      org0D[j][i]=tmp0D;
    };
    org1D[j].resize(sns[j].GetSize1D());
    for(int i=0;i<sns[j].GetSize1D();i++){
      tmp1D=sns[j].GetSensitivity1D(i);
      org1D[j][i]=tmp1D;
    };
    org2D[j].resize(sns[j].GetSize2D());
    for(int i=0;i<sns[j].GetSize2D();i++){
      tmp2D=sns[j].GetSensitivity2D(i);
      org2D[j][i]=tmp2D;
    };
  };

  real sum_val=0.;
  for(int i=0;i<num;i++){
    sum_val+=val_array[i];
  };
  val=sum_val/num;

  for(int j=0;j<size0D;j++){
    real ave0D=0.;
    for(int i=0;i<num;i++){
      ave0D=ave0D+org0D[i][j]*val_array[i];
    };
    sens0d[j]=ave0D/sum_val;
  };
  for(int j=0;j<size1D;j++){
    GroupData1D ave1D=org1D[0][j]*val_array[0];
    for(int i=1;i<num;i++){
      ave1D=ave1D+org1D[i][j]*val_array[i];
    };
    sens1d[j]=ave1D*(1./sum_val);
  };
  for(int j=0;j<size2D;j++){
    GroupData2D ave2D=org2D[0][j]*val_array[0];
    for(int i=1;i<num;i++){
      ave2D=ave2D+org2D[i][j]*val_array[i];
    };
    sens2d[j]=ave2D*(1./sum_val);
  };

};

void SensitivityData::WithdrawSensitivityData(SensitivityData &sns2)
{
  for(int i=0;i<GetSize0D();i++){
    int mat=GetMatList0D(i);
    int mt=GetMtList0D(i);
    int id=sns2.FindData0D(mat,mt);
    if(id==-1){
      //PutSensitivity0D(mat,mt,sns2.GetSensitivity0D(mat,mt));
    }else{
      real org=GetSensitivity0D(mat,mt);
      GetSensitivity0D(mat,mt)=org-sns2.GetSensitivity0D(mat,mt);
    };
  };

  for(int i=0;i<GetSize1D();i++){
    int mat=GetMatList1D(i);
    int mt=GetMtList1D(i);
    int id=sns2.FindData1D(mat,mt);
    if(id==-1){
      //PutSensitivity1D(mat,mt,sns2.GetSensitivity1D(mat,mt));
    }else{
      GroupData1D org=GetSensitivity1D(mat,mt);
      GetSensitivity1D(mat,mt)=org-sns2.GetSensitivity1D(mat,mt);
    };
  };

  for(int i=0;i<GetSize2D();i++){
    int mat=GetMatList2D(i);
    int mt=GetMtList2D(i);
    int id=sns2.FindData2D(mat,mt);
    if(id==-1){
      //PutSensitivity2D(mat,mt,sns2.GetSensitivity2D(mat,mt));
    }else{
      GroupData2D org=GetSensitivity2D(mat,mt);
      GetSensitivity2D(mat,mt)=org-sns2.GetSensitivity2D(mat,mt);
    };
  };

};

SensitivityData SensitivityData::CalSensitivityDataForReactivity(SensitivityData &sns2)
{
  // (dr/r)/(ds/s) = -(dk/k)/(ds/s) (1/k) + (dk'/k')/(ds/s)(1/k')
  // 
  // Actually, CBZ calculates k-reactivity as s=(dk/kk')/(ds/s), 
  // so (dr/r)/(ds/s) is approximated as -(dk/kk)/(ds/s)+(dk'/k'k')/(ds/s)

  real k1=GetValue();
  real k2=sns2.GetValue();
  real rho=1./k1-1./k2;
  real k1_inv=1./k1;
  real k2_inv=1./k2;
  real rho_inv=1./rho;

  int grp1=GetGroup();
  int grp2=sns2.GetGroup();
  if(grp1!=grp2){
    cout<<"# Error in SensitivityData::CalSensitivityDataForReactivity.\n";
    cout<<"# Inconsistency in energy group.\n";
    exit(0);
  };  

  SensitivityData snsr;
  snsr.PutName(GetName1(),"reactivity",GetName3());
  snsr.PutValue(rho);
  snsr.PutGroup(grp1);
  snsr.GetEnband()=GetEnband();

  int size0d=GetSize0D();
  int size0d_sec=sns2.GetSize0D();
  vector<bool> calc0d(size0d_sec,false);
  for(int i=0;i<size0d;i++){
    int mat=GetMatList0D(i);
    int mt=GetMtList0D(i);
    int tmp=sns2.FindData0D(mat,mt);
    if(tmp!=-1){
      //real snsnew=GetSensitivity0D(mat,mt)*(-k1_inv)+sns2.GetSensitivity0D(mat,mt)*k2_inv;
      real snsnew=GetSensitivity0D(mat,mt)*(-1)+sns2.GetSensitivity0D(mat,mt);
      snsnew*=rho_inv;
      snsr.PutSensitivity0D(mat,mt,snsnew);
      calc0d[tmp]=true;
    }else{
      real snsnew=GetSensitivity0D(mat,mt)*(-1);
      snsnew*=rho_inv;
      snsr.PutSensitivity0D(mat,mt,snsnew);
    };
  };
  for(int i=0;i<size0d_sec;i++){
    if(!calc0d[i]){
      int mat=sns2.GetMatList0D(i);
      int mt=sns2.GetMtList0D(i);
      real snsnew=sns2.GetSensitivity0D(mat,mt);
      snsnew*=rho_inv;
      snsr.PutSensitivity0D(mat,mt,snsnew);
    };
  };

  int size1d=GetSize1D();
  int size1d_sec=sns2.GetSize1D();
  vector<bool> calc1d(size1d_sec,false);
  for(int i=0;i<size1d;i++){
    int mat=GetMatList1D(i);
    int mt=GetMtList1D(i);
    int tmp=sns2.FindData1D(mat,mt);
    if(tmp!=-1){
      //GroupData1D snsnew=GetSensitivity1D(mat,mt)*(-k1_inv)+sns2.GetSensitivity1D(mat,mt)*k2_inv;
      GroupData1D snsnew=GetSensitivity1D(mat,mt)*(-1)+sns2.GetSensitivity1D(mat,mt);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity1D(mat,mt,snsnew);
      calc1d[tmp]=true;
    }else{
      GroupData1D snsnew=GetSensitivity1D(mat,mt)*(-1);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity1D(mat,mt,snsnew);
    };
  };
  for(int i=0;i<size1d_sec;i++){
    if(!calc1d[i]){
      int mat=sns2.GetMatList1D(i);
      int mt=sns2.GetMtList1D(i);
      GroupData1D snsnew=sns2.GetSensitivity1D(mat,mt);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity1D(mat,mt,snsnew);
    };
  };

  int size2d=GetSize2D();
  int size2d_sec=sns2.GetSize2D();
  vector<bool> calc2d(size2d_sec,false);
  for(int i=0;i<size2d;i++){
    int mat=GetMatList2D(i);
    int mt=GetMtList2D(i);
    int tmp=sns2.FindData2D(mat,mt);
    if(tmp!=-1){
      //GroupData2D snsnew=GetSensitivity2D(mat,mt)*(-k1_inv)+sns2.GetSensitivity2D(mat,mt)*k2_inv;
      GroupData2D snsnew=GetSensitivity2D(mat,mt)*(-1)+sns2.GetSensitivity2D(mat,mt);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity2D(mat,mt,snsnew);
      calc2d[tmp]=true;
    }else{
      GroupData2D snsnew=GetSensitivity2D(mat,mt)*(-1);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity2D(mat,mt,snsnew);
    };
  };
  for(int i=0;i<size2d_sec;i++){
    if(!calc2d[i]){
      int mat=sns2.GetMatList2D(i);
      int mt=sns2.GetMtList2D(i);
      GroupData2D snsnew=sns2.GetSensitivity2D(mat,mt);
      snsnew=snsnew*rho_inv;
      snsr.PutSensitivity2D(mat,mt,snsnew);
    };
  };


  return snsr;
};

SensitivityData SensitivityData::CalSensitivityDataForBeta(SensitivityData &sns2)
{
  real k1=GetValue();
  real k2=sns2.GetValue();
  real beta=(k2-k1)/k1;
  real factor=k2/k1/beta;

  int grp1=GetGroup();
  int grp2=sns2.GetGroup();
  if(grp1!=grp2){
    cout<<"# Error in SensitivityData::CalSensitivityDataForReactivity.\n";
    cout<<"# Inconsistency in energy group.\n";
    exit(0);
  };  

  SensitivityData snsr;
  snsr.PutName(GetName1(),"beta",GetName3());
  snsr.PutValue(beta);
  snsr.PutGroup(grp1);
  snsr.GetEnband()=GetEnband();

  int size0d=GetSize0D();
  for(int i=0;i<size0d;i++){
    int mat=GetMatList0D(i);
    int mt=GetMtList0D(i);
    int tmp=sns2.FindData0D(mat,mt);
    if(tmp!=-1){
      //real snsnew=GetSensitivity0D(mat,mt)*(-k1)+sns2.GetSensitivity0D(mat,mt)*k2;
      real snsnew=GetSensitivity0D(mat,mt)*(-1.)+sns2.GetSensitivity0D(mat,mt);
      snsnew*=factor;
      snsr.PutSensitivity0D(mat,mt,snsnew);
    };
  };

  int size1d=GetSize1D();
  for(int i=0;i<size1d;i++){
    int mat=GetMatList1D(i);
    int mt=GetMtList1D(i);
    int tmp=sns2.FindData1D(mat,mt);
    if(tmp!=-1){
      //GroupData1D snsnew=GetSensitivity1D(mat,mt)*(-k1)+sns2.GetSensitivity1D(mat,mt)*k2;
      GroupData1D snsnew=GetSensitivity1D(mat,mt)*(-1.)+sns2.GetSensitivity1D(mat,mt);
      snsnew=snsnew*factor;
      snsr.PutSensitivity1D(mat,mt,snsnew);
    };
  };

  int size2d=GetSize2D();
  for(int i=0;i<size2d;i++){
    int mat=GetMatList2D(i);
    int mt=GetMtList2D(i);
    int tmp=sns2.FindData2D(mat,mt);
    if(tmp!=-1){
      //GroupData2D snsnew=GetSensitivity2D(mat,mt)*(-k1)+sns2.GetSensitivity2D(mat,mt)*k2;
      GroupData2D snsnew=GetSensitivity2D(mat,mt)*(-1.)+sns2.GetSensitivity2D(mat,mt);
      snsnew=snsnew*factor;
      snsr.PutSensitivity2D(mat,mt,snsnew);
    };
  };

  return snsr;
};

void SensitivityData::CalSensitivity1DFrom2D(int mat,int mt)
{
  if(FindData2D(mat,mt)!=-1&&FindData1D(mat,mt)==-1){
    GroupData1D s1d(group);
    s1d.set_zero();
    for(int g=0;g<group;g++){
      real tmp=0.;
      for(int g2=0;g2<group;g2++){
	tmp+=GetSensitivity2D(mat,mt).get_dat(g,g2);
      };
      s1d.put_data(g,tmp);
    };
    PutSensitivity1D(mat,mt,s1d);
  }else{
    cout<<"# Warning in SensitivityData::CalSensitivity1DFrom2D.\n";
    cout<<"# Process is NOT done.\n";
  };
};

void SensitivityData::ShowSensitivity1D(int mat,int mt)
{
  int id=FindData1D(mat,mt);
  if(id!=-1){
    cout<<"# Sensitivity printing\n#\n";
    cout<<"#   Name1 : "<<name1<<"\n";
    cout<<"#   Name2 : "<<name2<<"\n";
    cout<<"#   Name3 : "<<name3<<"\n";
    cout<<"#   MAT   : "<<mat<<"\n";
    cout<<"#   MT    : "<<mt<<"\n#\n";
    cout<<"# Upper       Sensitivity/    Sensitivity   Sensitivity/\n";
    cout<<"# energy      lethargy*0.25                 lethargy\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    real sum=0.;
    real sum1=0.;
    real sum2=0.;
    real sum3=0.;
    int eb0,eb1;
    for(int g=0;g<group;g++){
      real e0=enband.get_dat(g);
      real e1=enband.get_dat(g+1);
      real letwid=log(e0/e1);
      real sns=GetSensitivity1D(mat,mt).get_dat(g);
      sum+=sns;
      if(e0>1e5){
	sum1+=sns;
	eb0=g;
      }else if(e0>4.){
	sum2+=sns;
	eb1=g;
      }else{
	sum3+=sns;
      };
      cout<<e0<<"   ";
      if(sns>0.)cout<<" ";
      cout<<sns/letwid*0.25<<"   ";
      if(sns>0.)cout<<" ";
      cout<<sns<<"   ";
      if(sns>0.)cout<<" ";
      cout<<sns/letwid<<"   ";
      cout<<"\n";
    };
    cout<<"# Total      : "<<sum<<"\n";
    cout<<"#  > 0.1 MeV : "<<sum1<<" (1-grp to "<<eb0+1<<"-grp)\n";
    cout<<"#  > 4 eV    : "<<sum2<<" ("<<eb0+2<<"-grp to "<<eb1+1<<"-grp)\n";
    cout<<"#  < 4 eV    : "<<sum3<<" ("<<eb1+2<<"-grp to "<<group+1<<"-grp)\n";
  };
};

void SensitivityData::ShowSensitivity1DExcel(int mat,int mt)
{
  int id=FindData1D(mat,mt);
  if(id!=-1){
    cout<<"# Sensitivity printing\n#\n";
    cout<<"#   Name1 : "<<name1<<"\n";
    cout<<"#   Name2 : "<<name2<<"\n";
    cout<<"#   Name3 : "<<name3<<"\n#\n";
    cout<<"# Energy      Sensitivity/\n";
    cout<<"#             lethary*0.25\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    for(int g=0;g<group;g++){
      real e0=enband.get_dat(g);
      real e1=enband.get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<e0<<"   "<<GetSensitivity1D(mat,mt).get_dat(g)/letwid*0.25<<"\n";
      cout<<e1<<"   "<<GetSensitivity1D(mat,mt).get_dat(g)/letwid*0.25<<"\n";
    };
  };
};

void SensitivityData::ShowSelf(real minimum)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  int s0d=GetSize0D();
  int s1d=GetSize1D();
  int s2d=GetSize2D();
  cout<<"# Sensitivity data\n";
  cout<<"#    Tag-1 : "<<name1<<"\n";
  cout<<"#    Tag-2 : "<<name2<<"\n";
  cout<<"#    Tag-3 : "<<name3<<"\n";
  cout<<"#    Value : "<<val<<"\n";
  cout<<"#\n";
  cout<<"#  0D data : "<<s0d<<"\n";
  for(int i=0;i<s0d;i++){
    if(fabs(sens0d[i])>minimum){
      cout<<"#       "<<mat_list0d[i]<<" "<<mt_list0d[i]<<" "<<sens0d[i]<<"\n";
    };
  };
  cout<<"#  1D data : "<<s1d<<"\n";
  for(int i=0;i<s1d;i++){
    real tmp=sens1d[i].get_sum();
    if(fabs(tmp)>minimum){
      cout<<"#       "<<mat_list1d[i]<<" "<<mt_list1d[i]<<" ";
      cout<<tmp;
      cout<<"\n";
    };
  };
  cout<<"#  2D data : "<<s2d<<"\n";
  for(int i=0;i<s2d;i++){
    real tmp=sens2d[i].get_sumx().get_sum();
    if(fabs(tmp)>minimum){
      cout<<"#       "<<mat_list2d[i]<<" "<<mt_list2d[i]<<" ";
      cout<<tmp;
      cout<<"\n";
    };
  };
};

void SensitivityData::ShowSelf(MATIDTranslator &midt,real minimum)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  int s0d=GetSize0D();
  int s1d=GetSize1D();
  int s2d=GetSize2D();
  cout<<"# Sensitivity data\n";
  cout<<"#    "<<name1<<"\n";
  cout<<"#    "<<name2<<"\n";
  cout<<"#    "<<name3<<"\n";
  cout<<"#  0D data : "<<s0d<<"\n";
  real yldsns_sum=0.;
  for(int i=0;i<s0d;i++){
    if(fabs(sens0d[i])>minimum){
      cout<<"#       "<<midt.Name(mat_list0d[i])<<" "<<mt_list0d[i]<<" ";
      if(mt_list0d[i]==8888){
	cout<<" (half-life) ";
      }else if(mt_list0d[i]>18000000){
	cout<<" (yield:"<<midt.Name(mt_list0d[i]-18000000)<<") ";
	yldsns_sum+=sens0d[i];
      };
      cout<<" "<<sens0d[i]<<"\n";
    };
  };
  if(yldsns_sum!=0.){
    cout<<"# (Sum of yield sensitivity : "<<yldsns_sum<<" )\n";
  };
  cout<<"#  1D data : "<<s1d<<"\n";
  for(int i=0;i<s1d;i++){
    real tmp=sens1d[i].get_sum();
    if(fabs(tmp)>minimum){
      cout<<"#       "<<midt.Name(mat_list1d[i])<<" "<<mt_list1d[i]<<" ";
      if(mt_list1d[i]==18){
	cout<<" (fission) ";
      }else if(mt_list1d[i]==16){
	cout<<" (n,2n)    ";
      }else if(mt_list1d[i]==102){
	cout<<" (capture) ";
      };
      cout<<tmp;
      cout<<"\n";
    };
  };
  cout<<"#  2D data : "<<s2d<<"\n";
  for(int i=0;i<s2d;i++){
    real tmp=sens2d[i].get_sumx().get_sum();
    if(fabs(tmp)>minimum){
      cout<<"#       "<<midt.Name(mat_list2d[i])<<" "<<mt_list2d[i]<<" ";
      cout<<tmp;
      cout<<"\n";
    };
  };
};

real SensitivityData::GetFPYieldSensitivityFissileWise(int fisid)
{
  int s0d=GetSize0D();
  real sum=0.;
  for(int i=0;i<s0d;i++){
    int tmp=18000000+fisid;
    if(mt_list0d[i]==tmp){
      sum+=sens0d[i];
    };
  };  
  return sum;
};

real SensitivityData::GetFPYieldSensitivityFPWise(int fpid)
{
  int s0d=GetSize0D();
  real sum=0.;
  for(int i=0;i<s0d;i++){
    if(mat_list0d[i]==fpid){
      if(mt_list0d[i]>18000000&&mt_list0d[i]<19000000){
        sum+=sens0d[i];
      };
    };
  };  
  return sum;
};

void SensitivityData::ShowFPYieldSummary(MATIDTranslator &midt,real minimum)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  int s0d=GetSize0D();

  vector<int> fisid_list;
  vector<real> fissum;

  vector<int> fpid_list;
  vector<real> fpsum;

  real sum=0.;
  for(int i=0;i<s0d;i++){
    //cout<<"\n"<<i<<" "<<mat_list0d[i]<<" "<<mt_list0d[i]<<" "<<sens0d[i]<<"\n";
    if(mt_list0d[i]>18000000){
      sum+=sens0d[i];

      {
      int tmp=mt_list0d[i]-18000000;
      int sz=fisid_list.size();
      //cout<<"   [mt] "<<tmp<<"\n";
      //cout<<"   [sz] "<<sz<<"\n";
      int index=-1;
      for(int j=0;j<sz;j++){
	if(tmp==fisid_list[j])index=j;
	//cout<<"       .... "<<j<<" "<<fisid_list[j]<<" index="<<index<<"\n";
      };
      if(index==-1){
	fisid_list.push_back(tmp);
	fissum.push_back(sens0d[i]);
      }else{
	fissum[index]+=sens0d[i];
      };
      };

      {
      int tmp=mat_list0d[i];
      int sz=fpid_list.size();
      int index=-1;
      for(int j=0;j<sz;j++){
	if(tmp==fpid_list[j])index=j;
      };
      if(index==-1){
	fpid_list.push_back(tmp);
	fpsum.push_back(sens0d[i]);
      }else{
	fpsum[index]+=sens0d[i];
      };
      };

    };

  };

  cout<<"\n# Sum of yield sensitivity : "<<sum<<"\n";

  {
    cout<<"\n\n# Fissile-wise yield sensitivity\n#\n";
  int sz=fisid_list.size();
  for(int i=0;i<sz;i++){
    if(fabs(fissum[i])>minimum){
      cout<<"# "<<midt.Name(fisid_list[i])<<" "<<fisid_list[i]<<" "<<fissum[i]<<"\n";
    };
  };
  };

  {
    cout<<"\n\n# FP-wise yield sensitivity\n#\n";
  int sz=fpid_list.size();
  for(int i=0;i<sz;i++){
    if(fabs(fpsum[i])>minimum){
      cout<<"# "<<midt.Name(fpid_list[i])<<" "<<fpid_list[i]<<" "<<fpsum[i]<<"\n";
    };
  };
  };

};


void SensitivityData::WriteFile(string mdir,string filename)
{
  mdir.append(filename);

  ofstream fout;
  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# "<<mdir<<"\n";
    exit(1);
  };

  fout.setf(ios::scientific);
  fout.precision(8);

  if(name1=="")name1="dummy";
  if(name2=="")name2="dummy";
  if(name3=="")name3="dummy";

  fout<<name1<<"\n";
  fout<<name2<<"\n";
  fout<<name3<<"\n";
  fout<<" "<<val<<"\n";
  fout<<" "<<group<<"\n";

  int size0d=sens0d.size();
  int size1d=sens1d.size();
  int size2d=sens2d.size();

  fout<<" "<<size0d<<"\n";
  fout<<" "<<size1d<<"\n";
  fout<<" "<<size2d<<"\n";

  for(int i=0;i<group+1;i++){
    fout<<"  "<<enband.get_dat(i)<<"\n";
  };

  for(int i=0;i<size0d;i++){
    fout<<" "<<mat_list0d[i]<<"\n";
    fout<<" "<<mt_list0d[i]<<"\n";
    fout<<"  "<<sens0d[i]<<"\n";
  };

  for(int i=0;i<size1d;i++){
    int stgrp=0;
    int edgrp=group-1;
    for(int g=0;g<group;g++){
      if(sens1d[i].get_dat(g)==0.&&stgrp==g)stgrp=g+1;
    };
    for(int g=group-1;g>=0;g--){
      if(sens1d[i].get_dat(g)==0.&&edgrp==g)edgrp=g-1;
    };
    fout<<" "<<mat_list1d[i]<<"\n";
    fout<<" "<<mt_list1d[i]<<"\n";
    fout<<"  "<<stgrp<<"\n";
    fout<<"  "<<edgrp<<"\n";
    for(int g=stgrp;g<=edgrp;g++){
      fout<<"   "<<sens1d[i].get_dat(g)<<"\n";
    };
  };


  for(int i=0;i<size2d;i++){
    fout<<" "<<mat_list2d[i]<<"\n";
    fout<<" "<<mt_list2d[i]<<"\n";
    for(int g=0;g<group;g++){
      int stgrp=0;
      int edgrp=group-1;
      for(int g2=0;g2<group;g2++){
	if(sens2d[i].get_dat(g,g2)==0.&&stgrp==g2)stgrp=g2+1;
      };
      for(int g2=group-1;g2>=0;g2--){
	if(sens2d[i].get_dat(g,g2)==0.&&edgrp==g2)edgrp=g2-1;
      };
      if(stgrp<=edgrp){
	fout<<"  "<<g<<"\n";
	fout<<"  "<<stgrp<<"\n";
	fout<<"  "<<edgrp<<"\n";
	for(int g2=stgrp;g2<=edgrp;g2++){
	  fout<<"   "<<sens2d[i].get_dat(g,g2)<<"\n";
	};
      };
    };
    fout<<" -1\n";
  };

  fout.close();
};

void SensitivityData::ReadFile(string mdir,string filename)
{
  mdir.append(filename);

  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# "<<mdir<<"\n";
    exit(1);
  };

  string tmp1,tmp2,tmp3;
  fin>>tmp1;
  fin>>tmp2;
  fin>>tmp3;
  PutName(tmp1,tmp2,tmp3);

  real rtmp;
  fin>>rtmp;
  PutValue(rtmp);

  int itmp;
  fin>>itmp;
  PutGroup(itmp);

  int size0d,size1d,size2d;

  fin>>size0d;
  fin>>size1d;
  fin>>size2d;

  for(int i=0;i<group+1;i++){
    fin>>rtmp;
    PutEnergyBoundaryData(i,rtmp);
  };

  GroupData1D sns1d(group);
  GroupData2D sns2d(group,group);

  for(int i=0;i<size0d;i++){
    int mat,mt;
    fin>>mat;
    fin>>mt;
    fin>>rtmp;
    PutSensitivity0D(mat,mt,rtmp);
  };

  for(int i=0;i<size1d;i++){
    int mat,mt;
    fin>>mat;
    fin>>mt;
    int stgrp,edgrp;
    fin>>stgrp;
    fin>>edgrp;
    sns1d.set_zero();
    for(int g=stgrp;g<=edgrp;g++){
      fin>>rtmp;
      sns1d.put_data(g,rtmp);
    };
    PutSensitivity1D(mat,mt,sns1d);
  };

  for(int i=0;i<size2d;i++){
    int mat,mt;
    fin>>mat;
    fin>>mt;
    sns2d.set_zero();
    bool cont=true;
    while(cont){
      fin>>itmp;
      if(itmp!=-1){
        int stgrp,edgrp;
        fin>>stgrp;
        fin>>edgrp;
        for(int g=stgrp;g<=edgrp;g++){
          fin>>rtmp;
          sns2d.put_data(itmp,g,rtmp);
        };
      }else{
	cont=false;
        PutSensitivity2D(mat,mt,sns2d);
      };
    };
  };

  fin.close();
};

void SensitivityData::ReadFileOldFormat(string mdir,string filename)
{
  MATIDTranslator midt;
  
  mdir.append(filename);

  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# "<<mdir<<"\n";
    exit(1);
  };

  PutName("dummy","dummy","dummy");

  real rtmp=0.; // dummy 
  PutValue(rtmp);

  int itmp;
  fin>>itmp;
  PutGroup(itmp);

  GroupData1D sns1d(group);
  GroupData2D sns2d(group,group);

  itmp=1;

  while(itmp==1||itmp==2){

    fin>>itmp;
    if(itmp==1){ // 1D data
      int mat,mt;
      fin>>mat;
      fin>>mt;
      mat=midt.GetMATIDFromENDFID(mat);
      for(int g=0;g<group;g++){
        fin>>rtmp;
        sns1d.put_data(g,rtmp);
      };
      PutSensitivity1D(mat,mt,sns1d);
    }else if(itmp==2){ // 2D data
      int mat,mt;
      fin>>mat;
      fin>>mt;
      mat=midt.GetMATIDFromENDFID(mat);
      sns2d.set_zero();
      for(int g=0;g<group;g++){
        for(int g2=g;g2<group;g2++){
          fin>>rtmp;
          sns2d.put_data(g,g2,rtmp);
        };
      };
      PutSensitivity2D(mat,mt,sns2d);      
    };

  };

  fin.close();
};


void SensitivityData::ChangeMatMt(int mato,int mto,int matn,int mtn)
{
  for(int i=0;i<GetSize0D();i++){
    int mat=GetMatList0D(i);
    int mt=GetMtList0D(i);
    if(mat==mato&&mt==mto){
      mat_list0d[i]=matn;
      mt_list0d[i]=mtn;
    };
  };

  for(int i=0;i<GetSize1D();i++){
    int mat=GetMatList1D(i);
    int mt=GetMtList1D(i);
    if(mat==mato&&mt==mto){
      mat_list1d[i]=matn;
      mt_list1d[i]=mtn;
    };
  };

  for(int i=0;i<GetSize2D();i++){
    int mat=GetMatList2D(i);
    int mt=GetMtList2D(i);
    if(mat==mato&&mt==mto){
      mat_list2d[i]=matn;
      mt_list2d[i]=mtn;
    };
  };
};

void SensitivityData::TransformToConstrainedSensitivity(int mat,int mt,GroupData1D &xs)
{
  GroupData1D org;
  org=GetSensitivity1D(mat,mt);
  real sum=org.get_sum();
  GroupData1D newd(group);
  newd=org-xs*sum;
  GetSensitivity1D(mat,mt)=newd;
};

void SensitivityData::Factorize(real factor)
{
  for(int i=0;i<GetSize0D();i++){
    sens0d[i]*=factor;
  };

  for(int i=0;i<GetSize1D();i++){
    sens1d[i]=sens1d[i]*factor;
  };

  for(int i=0;i<GetSize2D();i++){
    sens2d[i]=sens2d[i]*factor;
  };
};

void SensitivityData::Diff(SensitivityData& sns2,real threshold,real sns_threshold)
{
  cout<<"#\n";
  cout<<"# Difference to reference sensitivity data\n";
  cout<<"#\n";
  cout<<"# Threshold of sensitivity        : "<<sns_threshold<<"\n";
  cout<<"# Threshold of relative differnce : "<<threshold<<"\n";
  cout<<"#\n";

  // +++ 0D
  cout<<"##############################################################\n";
  cout<<"# 0D data\n";
  cout<<"##############################################################\n";
  for(int i=0;i<GetSize0D();i++){
    int mat=mat_list0d[i];
    int mt=mt_list0d[i];
    real dat1=GetSensitivity0D(mat,mt);
    if(fabs(dat1)>sns_threshold){
      int id=sns2.FindData0D(mat,mt);
      if(id==-1){
        cout<<"# "<<mat<<"/"<<mt<<" [no data in the reference]\n";
      }else{
        real dat2=sns2.GetSensitivity0D(mat,mt);
        real rdif=(dat2-dat1)/dat1;
        if(fabs(rdif)>threshold){
          cout<<"# "<<mat<<"/"<<mt<<" : "<<rdif<<"  (Ref.-value : "<<dat2<<")\n";
	};
      };
    };
  };

  // +++ 1D
  cout<<"##############################################################\n";
  cout<<"# 1D data\n";
  cout<<"##############################################################\n";
  for(int i=0;i<GetSize1D();i++){
    int mat=mat_list1d[i];
    int mt=mt_list1d[i];
    int id=sns2.FindData1D(mat,mt);
    if(id==-1){
      cout<<"# "<<mat<<"/"<<mt<<" [no data in the reference]\n";
    }else{
      GroupData1D dat1v=GetSensitivity1D(mat,mt);
      GroupData1D dat2v=sns2.GetSensitivity1D(mat,mt);
      int grp=dat1v.get_imax();
      for(int g=0;g<grp;g++){
        real dat1=dat1v.get_dat(g);
        if(fabs(dat1)>sns_threshold){
	  real dat2=dat2v.get_dat(g);
          real rdif=(dat2-dat1)/dat1;
          if(fabs(rdif)>threshold){
	    cout<<"# "<<mat<<"/"<<mt<<"/"<<g<<" : "<<rdif<<"  (Ref.-value : "<<dat2<<")\n";
	  };
	};
      };
    };
  };

};

SensitivityData SensitivityData::EnergyGroupTranslation(GroupData1D &ebnd)
{
  int grp2=ebnd.get_imax()-1;

  SensitivityData sens2;
  sens2.PutName(name1,name2,name3);
  sens2.PutValue(val);
  sens2.PutGroup(grp2);
  for(int i=0;i<=grp2;i++){
    sens2.PutEnergyBoundaryData(i,ebnd.get_dat(i));
  };

  for(int i=0;i<GetSize1D();i++){

    GroupData1D sorg=GetSensitivity1D(i);
    GroupData1D snew(grp2);
    snew.set_zero();

    for(int g=0;g<group;g++){

      real e=enband.get_dat(g);
      real e2=enband.get_dat(g+1);
      int gs=-1;
      int ge=-1;
      for(int j=0;j<grp2;j++){
	real es=ebnd.get_dat(j);
	real es2=ebnd.get_dat(j+1);
	if(es>=e&&es2<e)gs=j;
	if(es>e2&&es2<=e2)ge=j;
      };
      if(gs==-1&&ge==-1){
      }else if(gs==ge){
        snew.add_data(gs,sorg.get_dat(g));
      }else{
	for(int j=gs;j<=ge;j++){
	  real es=ebnd.get_dat(j);
	  real es2=ebnd.get_dat(j+1);
          real wgt=0.;
	  if(es>=e&&es2>e2){
	    wgt=log(e/es2)/log(e/e2);
	  }else if(es<e&&es2>e2){
	    wgt=log(es/es2)/log(e/e2);
	  }else if(es<e&&es2>e2){
            wgt=log(es/es2)/log(e/e2);
	  }else if(es<e&&es2<=e2){
            wgt=log(es/e2)/log(e/e2);
	  };
	  if(wgt==0){
	    cout<<"# Error!\n";
	    exit(0);
	  };
          snew.add_data(j,sorg.get_dat(g)*wgt);
	};
      };
    };
    /*
    int g=0;
    bool detect=false;
    real e,e2;
    while(!detect){
      e=enband.get_dat(g);
      e2=enband.get_dat(g+1);
      if(e2>ebnd.get_dat(0)){
        g++;
      }else{
        detect=true;
      };
    };

    int gs=0;
    real es=ebnd.get_dat(gs);
    real es2=ebnd.get_dat(gs+1);
    
    if(e2>es2){
      snew.add_data(gs,sorg.get_dat(g)*(log(es/e2)/log(e/e2)));
      g++;
    }else{
    */

    sens2.PutSensitivity1D(mat_list1d[i],mt_list1d[i],snew);

  };

  return sens2;
  
};

void SensitivityData::CompareAnotherSensitivityData(SensitivityData& sd2)
{
  if(group!=sd2.GetGroup()){
    cout<<"# Error in SensitivityData::CompareAnotherSensitivityData.\n";
    cout<<"# The number of energy groups is inconsistent.\n";
    exit(0);
  };

  cout.setf(ios::scientific);
  cout.precision(3);

  cout<<"# 1D data\n";
  for(int i=0;i<GetSize1D();i++){
    int mat=mat_list1d[i];
    int mt=mt_list1d[i];
    for(int i2=0;i2<sd2.GetSize1D();i2++){
      if(mat==sd2.GetMatList1D(i2)&&mt==sd2.GetMtList1D(i2)){
        real sum_sq=0.;
        real sum_sq_rel=0.;
        real max=0.;
        real max_rel=0.;
        int count=0;
	for(int g=0;g<group;g++){
	  real ss1=sens1d[i].get_dat(g);
	  real ss2=sd2.GetSensitivity1D(i2).get_dat(g);
	  real diff=ss2-ss1;
          real diff_rel=0.;
	  if(ss1!=0.)count++;
          if(ss1!=0.)diff_rel=(ss2-ss1)/ss1;
	  sum_sq+=diff*diff;
	  sum_sq_rel+=diff_rel*diff_rel;
	  if(max<fabs(diff))max=fabs(diff);
	  if(max_rel<fabs(diff_rel))max_rel=fabs(diff_rel);
	};
	sum_sq/=count;
	sum_sq=sqrt(sum_sq);
	sum_sq_rel/=count;
	sum_sq_rel=sqrt(sum_sq_rel);
	cout<<"# "<<mat<<" "<<mt<<" "<<sum_sq<<" "<<max<<" "<<sum_sq_rel<<" "<<max_rel<<"\n";
      };
    };
  };

  cout<<"# 2D data\n";
  for(int i=0;i<GetSize2D();i++){
    int mat=mat_list2d[i];
    int mt=mt_list2d[i];
    for(int i2=0;i2<sd2.GetSize2D();i2++){
      if(mat==sd2.GetMatList2D(i2)&&mt==sd2.GetMtList2D(i2)){
        real sum_sq=0.;
        real sum_sq_rel=0.;
        real max=0.;
        real max_rel=0.;
	int count=0;
	for(int g=0;g<group;g++){
	  for(int g2=0;g2<group;g2++){
  	    real ss1=sens2d[i].get_dat(g,g2);
	    real ss2=sd2.GetSensitivity2D(i2).get_dat(g,g2);
	    real diff=ss2-ss1;
            real diff_rel=0.;
	    if(ss1!=0.)count++;
            if(ss1!=0.)diff_rel=(ss2-ss1)/ss1;
	    sum_sq+=diff*diff;
	    sum_sq_rel+=diff_rel*diff_rel;
	    if(max<fabs(diff))max=fabs(diff);
	    if(max_rel<fabs(diff_rel))max_rel=fabs(diff_rel);
	  };
	};
	sum_sq/=count;
	sum_sq=sqrt(sum_sq);
	sum_sq_rel/=count;
	sum_sq_rel=sqrt(sum_sq_rel);
	cout<<"# "<<mat<<" "<<mt<<" "<<sum_sq<<" "<<max<<" "<<sum_sq_rel<<" "<<max_rel<<"\n";
      };
    };
  };

};

SensitivityData SensitivityData::Compressing(int num, string *nuclist)
{
  MATIDTranslator midt;

  SensitivityData ret;
  ret.PutName(name1,name2,name3);
  ret.PutValue(val);
  ret.PutGroup(group);
  ret.GetEnband().copy(enband);

  int sz0d=mat_list0d.size();
  for(int i=0;i<sz0d;i++){
    int mat=mat_list0d[i];
    int mt=mt_list0d[i];
    ret.PutSensitivity0D(mat,mt,sens0d[i]);
  };

  int sz1d=mat_list1d.size();
  for(int i=0;i<sz1d;i++){
    int mat=mat_list1d[i];
    bool exist=false;
    for(int j=0;j<num;j++){
      if(mat==midt.ID(nuclist[j]))exist=true;
    };
    if(exist){
      int mt=mt_list1d[i];
      ret.PutSensitivity1D(mat,mt,sens1d[i]);
    };
  };

  int sz2d=mat_list2d.size();
  for(int i=0;i<sz2d;i++){
    int mat=mat_list2d[i];
    bool exist=false;
    for(int j=0;j<num;j++){
      if(mat==midt.ID(nuclist[j]))exist=true;
    };
    if(exist){
      int mt=mt_list2d[i];
      ret.PutSensitivity2D(mat,mt,sens2d[i]);
    };
  };

  return ret;
};

real SensitivityData::CalEuclideanNormForAbsorptionXS()
{
  real norm=0.;
  
  int size1d=GetSize1D();
  for(int i=0;i<size1d;i++){
    int mat=GetMatList1D(i);
    int mt=GetMtList1D(i);
    if(mt==102){
      GroupData1D sns_abs=GetSensitivity1D(mat,mt);
      int tmp=FindData1D(mat,18);
      if(tmp!=-1){
	sns_abs=sns_abs+GetSensitivity1D(mat,18);
      };
      norm+=pow(sns_abs.GetEuclideanNorm(),2);
    };
  };

  return sqrt(norm);
};
