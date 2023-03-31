#include <cstdlib>
#include "FRDesign.h"

using namespace std;

// +++ FRDSubAssembly

FRDSubAssembly::FRDSubAssembly()
{
  n_ax_zone=0;
  area_xy=0.;
  irradiated_day=0.;
};

void FRDSubAssembly::PutNumberAxialZone(int i)
{
  n_ax_zone=i;
  matn_ax_zone.resize(i);
  tag_ax_zone.resize(i);
  medid_ax_zone.resize(i);
  mat_ax_zone.resize(i);
  den_ax_zone.resize(i);
  len_ax_zone.resize(i);
  heavy_metal_weight.resize(i,0.);
  generated_power.resize(i,0.);
  fast_flux_fluence.resize(i,0.);
};

void FRDSubAssembly::PutAxialLength(real *inp)
{
  if(n_ax_zone==0){
    cout<<"Error in FRDSubAssembly::PutLengthData.\n";
    cout<<"Number of axial zone is not yet put.\n";
    exit(0);
  };
  for(int i=0;i<n_ax_zone;i++){
    if(inp[i]<0.){
      cout<<"Error in FRDSubAssembly::PutLenghData.\n";
      cout<<"Zero value is put as "<<i<<"-th length.\n";
      exit(0);
    };
    len_ax_zone[i]=inp[i];
  };
};

void FRDSubAssembly::PutAxialTag(string *inp)
{
  if(n_ax_zone==0){
    cout<<"Error in FRDSubAssembly::PutLengthData.\n";
    cout<<"Number of axial zone is not yet put.\n";
    exit(0);
  };
  for(int i=0;i<n_ax_zone;i++){
    tag_ax_zone[i]=inp[i];
  };
};

void FRDSubAssembly::PutMatnum(int i1, int i2)
{
  if(i1<0||i1>n_ax_zone){
    cout<<"Error in FRDSubAssembly::PutMatnum.\n";
    cout<<"Zone ID "<<i1<<" is inappropriate.\n";
    exit(0);
  };
  matn_ax_zone[i1]=i2;
  mat_ax_zone[i1].resize(i2,0);
  den_ax_zone[i1].resize(i2,0.);
};

void FRDSubAssembly::PutMediumID(int i1, int i2)
{
  if(i1<0||i1>n_ax_zone){
    cout<<"Error in FRDSubAssembly::PutMatnum.\n";
    cout<<"Zone ID "<<i1<<" is inappropriate.\n";
    exit(0);
  };
  medid_ax_zone[i1]=i2;
};

void FRDSubAssembly::PutDensityData(int iax,int inum,int *id,real *den)
{
  PutMatnum(iax,inum);
  for(int i=0;i<inum;i++){
    mat_ax_zone[iax][i]=id[i];
    den_ax_zone[iax][i]=den[i];
  };
};

void FRDSubAssembly::PutDensityData(int iax, FRDTSubAssembly &frdt_sa)
{
  frdt_sa.CheckExistAllData();

  int fuel_num=frdt_sa.GetFuelComposition()->GetNucnum();
  int sus_num=frdt_sa.GetSUSComposition()->GetNucnum();
  int tot_num=fuel_num+sus_num+1;

  int *nucid=new int[tot_num];
  real *den=new real[tot_num];
  
  real sa_vol=frdt_sa.GetFuelSAGeometry()->GetAssemblyVolume();
  real pellet_vol=frdt_sa.GetFuelSAGeometry()->GetPelletVolume();
  real na_vol=frdt_sa.GetFuelSAGeometry()->GetCoolantVolume();
  real duct_vol=frdt_sa.GetFuelSAGeometry()->GetDuctVolume();
  real pin_vol=frdt_sa.GetFuelSAGeometry()->GetPinVolume();
  real wire_vol=frdt_sa.GetFuelSAGeometry()->GetSpacerwireVolume();

  int index=0;
  real factor=pellet_vol/sa_vol;

  for(int i=0;i<fuel_num;i++){
    int id=frdt_sa.GetFuelComposition()->GetNucid(i);
    real d=frdt_sa.GetFuelComposition()->GetDensity(i);
    nucid[index]=id;
    den[index++]=d*factor;
  };
  
  real na_den=frdt_sa.GetSodiumNumberDensity();
  nucid[index]=1125;
  den[index++]=na_den*na_vol/sa_vol;

  real factor2=(duct_vol+pin_vol+wire_vol)/sa_vol;
  for(int i=0;i<sus_num;i++){
    int id=frdt_sa.GetSUSComposition()->GetNucid(i);
    real d=frdt_sa.GetSUSComposition()->GetDensity(i);
    nucid[index]=id;
    den[index++]=d*factor2;
  };

  PutDensityData(iax,tot_num,nucid,den);

  delete [] nucid;
  delete [] den;
};

void FRDSubAssembly::ShowDensityData(int iax)
{
  for(int i=0;i<matn_ax_zone[iax];i++){
    cout<<mat_ax_zone[iax][i]<<" : "<<den_ax_zone[iax][i]<<"\n";
  };
};

void FRDSubAssembly::PutDataToAssembly(Assembly &asmi)
{
  int zdiv=GetNumberAxialZone();
  string *nn=new string[zdiv];
  real *ll=new real[zdiv];
  for(int i=0;i<zdiv;i++){
    nn[i]=GetAxialTag(i);
    ll[i]=GetAxialLength(i);
  };
  asmi.Init(zdiv,nn,ll,tag,"width");

  delete [] nn;
  delete [] ll;
};

void FRDSubAssembly::PutDensityDataToMedium(Medium &med,int z)
{
  med.SetZeroNumberDensity();
  for(int j=0;j<matn_ax_zone[z];j++){
    int matid=mat_ax_zone[z][j];
    real den=den_ax_zone[z][j];
    if(!med.ExistNuclide(matid)){
      cout<<"!! Warning !!\n";
      cout<<" FRDSubAssembly::PutDensityDataToMedium.\n";
      cout<<" Material "<<matid<<" is not included in medium.\n";
      return;
    };
    med.GetNuclide(matid).PutDensity(den);
  };
};

void FRDSubAssembly::GetDensityDataFromMedium(Medium &med,int z)
{
  for(int j=0;j<matn_ax_zone[z];j++){
    den_ax_zone[z][j]=0.;
  };

  for(int i=0;i<med.GetNucnum();i++){
    int matid=med.GetNuclideInTurn(i).GetMatnum();
    real den=med.GetNuclideInTurn(i).GetDensity();
    for(int j=0;j<matn_ax_zone[z];j++){
      if(mat_ax_zone[z][j]==matid){
	den_ax_zone[z][j]=den;
      };
    };
  };
};

void FRDSubAssembly::PutDensityDataToBurnup(Burnup &bu,int z)
{
  for(int n=0;n<matn_ax_zone[z];n++){
    int mid=mat_ax_zone[z][n];
    int id_bu=bu.SearchNuclide(mid);
    if(id_bu!=-1){
      real den=den_ax_zone[z][n];
      bu.PutDensity(id_bu,den);
    };
  };
};

void FRDSubAssembly::GetDensityDataFromBurnup(Burnup &bu,int z)
{
  for(int n=0;n<bu.GetNucnum();n++){
    int id_bu=bu.GetNuclideID(n);
    real den=bu.GetDensity(n);
    for(int nn=0;nn<matn_ax_zone[z];nn++){
      int mid=GetMatID(z,nn);
      if(id_bu==mid)PutDensity(z,nn,den);
    };
  };
};

void FRDSubAssembly::CalHeavyMetalWeight(Burnup &bu)
{
  real avo=0.60221367;

  for(int i=0;i<n_ax_zone;i++){
    real sum=0.;
    real vol=len_ax_zone[i]*area_xy;
    int matnum=matn_ax_zone[i];
    for(int j=0;j<matnum;j++){
      int matid=mat_ax_zone[i][j];
      if(matid>9000){
	real atomw=bu.GetAtomicWeight(matid);
	if(atomw>0.){
	  real den=den_ax_zone[i][j];
	  real mol=den/avo;
	  sum+=mol*atomw*vol;
	};
      };
    };
    heavy_metal_weight[i]=sum;
  };
};

real FRDSubAssembly::GetHeavyMetalWeight(string itag)
{
  real ret=0.;
  for(int i=0;i<n_ax_zone;i++){
    if(itag==""||itag==tag_ax_zone[i])ret+=heavy_metal_weight[i];
  };
  return ret;
};

real FRDSubAssembly::GetGeneratedPower(string itag)
{
  real ret=0.;
  for(int i=0;i<n_ax_zone;i++){
    if(itag==""||itag==tag_ax_zone[i])ret+=generated_power[i];
  };
  return ret;
};

real FRDSubAssembly::GetFastFluxFluence(string itag)
{
  real ret=0.;
  for(int i=0;i<n_ax_zone;i++){
    if(itag==""||itag==tag_ax_zone[i])ret+=fast_flux_fluence[i];
  };
  return ret;
};

real FRDSubAssembly::GetMaximumFastFluxFluence()
{
  real ret=0.;
  for(int i=0;i<n_ax_zone;i++){
    if(fast_flux_fluence[i]>ret)ret=fast_flux_fluence[i];
  };
  return ret;
};

void FRDSubAssembly::ReadFile(string mdir,string ss)
{
  mdir.append(ss);

  ifstream fin;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(0);
  };
  string stmp;
  int itmp;
  real tmp;
  fin>>stmp;
  PutTag(stmp);
  fin>>itmp;
  PutNumberAxialZone(itmp);
  for(int i=0;i<n_ax_zone;i++){
    fin>>itmp;
    PutMatnum(i,itmp);
    fin>>stmp;
    PutAxialTag(i,stmp);
    fin>>itmp;
    PutMediumID(i,itmp);
    fin>>tmp;
    PutAxialLength(i,tmp);
    for(int j=0;j<matn_ax_zone[i];j++){
      fin>>itmp;
      PutMatID(i,j,itmp);
      fin>>tmp;
      PutDensity(i,j,tmp);
    };
    fin>>tmp;
    PutHeavyMetalWeight(i,tmp);
    fin>>tmp;
    PutGeneratedPower(i,tmp);
    fin>>tmp;
    PutFastFluxFluence(i,tmp);
  };
  fin>>tmp;
  PutIrradiatedDay(tmp);

  fin.close();
};

void FRDSubAssembly::WriteFile(string mdir,string ss)
{
  mdir.append(ss);

  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(0);
  };
  fout<<tag<<"\n";
  fout<<n_ax_zone<<"\n";
  for(int i=0;i<n_ax_zone;i++){
    fout<<matn_ax_zone[i]<<"\n";
    fout<<tag_ax_zone[i]<<"\n";
    fout<<medid_ax_zone[i]<<"\n";
    fout<<len_ax_zone[i]<<"\n";
    for(int j=0;j<matn_ax_zone[i];j++){
      fout<<mat_ax_zone[i][j]<<"\n";
      fout.setf(ios::scientific);
      fout.precision(6);
      fout<<den_ax_zone[i][j]<<"\n";
    };
    fout<<heavy_metal_weight[i]<<"\n";
    fout<<generated_power[i]<<"\n";
    fout<<fast_flux_fluence[i]<<"\n";
  };
  fout<<irradiated_day<<"\n";

  fout.close();
};

// +++ FRDMediumSet

FRDMediumSet::FRDMediumSet(int igrp,int ipl,int iplt)
{
  mednum=0;
  grp=igrp;
  pl=ipl;
  plt=iplt;
  if(iplt==-1)plt=pl;
};

void FRDMediumSet::PutData(int n,int *id,real *den,string itag,real temp)
{
  mednum++;
  Medium min;
  min.PutImax(grp);
  min.PutPL(pl,plt);
  min.PutNuclide(n,id,den);
  min.PutTemperatureForAllNuclide(temp);
  med.push_back(min);
  medtag.push_back(itag);
};

void FRDMediumSet::PutData(int n,int *id,string itag)
{
  real *den=new real[n];
  for(int i=0;i<n;i++){den[i]=0.;};
  PutData(n,id,den,itag,0.);
  delete [] den;
};

void FRDMediumSet::PutDensityData(FRDTSubAssembly &frdt_sa,string itag)
{
  int ii=GetIDFromTag(itag);
  if(ii==-1){
    cout<<"!! Warning !!\n";
    cout<<"There is no data in FRDMediumSet : "<<itag<<"\n";
  };
  frdt_sa.PutHomogenizedNumberDensity(med[ii]);
  med[ii].PutTemperatureForAllNuclide(frdt_sa.GetTemperature());
  if(med[ii].ExistNuclide(1125)){
    med[ii].GetNuclide(1125).PutTemperature(frdt_sa.GetSodiumTemperature());
  };
};

int FRDMediumSet::GetIDFromTag(string itag)
{
  for(int i=0;i<mednum;i++){
    if(itag==medtag[i])return i;
  };
  return -1;
};

void FRDMediumSet::GroupCollapsing(int ngrp,int *bgrp,GroupData1D &wgt1,GroupData1D &wgt2)
{
  grp=ngrp;

  for(int i=0;i<mednum;i++){
    Medium bmed=med[i].Cond(ngrp,bgrp,wgt1,wgt2,true);
    med[i]=bmed;
  };
};

void FRDMediumSet::GroupCollapsing(string mtag,int ngrp,int *bgrp,GroupData1D &wgt1,GroupData1D &wgt2)
{
  int id=GetIDFromTag(mtag);
  if(id!=-1){
    grp=ngrp;
    Medium bmed=med[id].Cond(ngrp,bgrp,wgt1,wgt2,true);
    med[id]=bmed;
  };
};

// +++ FRDCoreLocationManager

FRDCoreLocationManager::FRDCoreLocationManager(int ix,int iy)
{
  mx=ix;
  my=iy;
  tag_map.resize(mx*my,"");
};

void FRDCoreLocationManager::SetCenterPosition(int ix,int iy)
{
  cx=ix;
  cy=iy;
};

void FRDCoreLocationManager::GetPositionFromAddress(string inp,int &ix,int &iy)
{
  int i1;
  string in;
  int i2;
  AddressStringTransformer(inp,i1,in,i2);
  GetPositionFromAddress(i1,in,i2,ix,iy);
};

void FRDCoreLocationManager::GetPositionFromAddress(int i1,string in,int i2,int &ix,int &iy)
{
  if(i1==0){
    ix=cx;
    iy=cy;
    return;
  };

  if(i1<i2){
    cout<<"Error in FRDCoreLocationManager::GetPositionFromAddress.\n";
    cout<<"Inappropriate address setting "<<i1<<" & "<<i2<<"\n";
    exit(0);
  };

  bool cmove=false;

  int xx=cx;
  int yy=cy;
  if(in=="A"||in=="F"){
    yy-=i1;
    if(i1%2==1)cmove=true;
    if(in=="A"){
      xx+=i1/2;
    }else{
      xx-=(i1+1)/2;
    };
  }else if(in=="B"){
    xx+=i1;
  }else if(in=="E"){
    xx-=i1;
  }else if(in=="C"||in=="D"){
    yy+=i1;
    if(i1%2==1)cmove=true;
    if(in=="C"){
      xx+=i1/2;
    }else{
      xx-=(i1+1)/2;
    };
  }else{
    cout<<"Error in FRDCoreLocationManager::GetPositionFromAddress.\n";
    cout<<"Inappropriate address : "<<i1<<in<<i2<<"\n";
    exit(0);
  };

  i2--;
  if(i2==0){
    ix=xx;
    iy=yy;
    return;
  };
  int tmpa=0;
  int tmpb=1;
  if(cmove){
    tmpa=1;
    tmpb=0;
  };
  if(in=="A"||in=="B"){
    yy+=i2;
    if(in=="A"){
      xx+=(i2+tmpa)/2;
    }else{
      xx-=(i2+tmpb)/2;
    };
  }else if(in=="C"){
    xx-=i2;
  }else if(in=="F"){
    xx+=i2;
  }else{
    yy-=i2;
    if(in=="E"){
      xx+=(tmpa+i2)/2;
    }else{
      xx-=(tmpb+i2)/2;
    };
  };

  ix=xx;
  iy=yy;
};


void FRDCoreLocationManager::SetAssembly(int i1,string in,int i2,FRDSubAssembly &frd_sa)
{
  if(in=="*"&&i1!=0){
    string tmp[]={"A","B","C","D","E","F"};
    for(int i=0;i<6;i++){
      SetAssembly(i1,tmp[i],i2,frd_sa);
    };
    return;
  };
  if(i2==0&&i1!=0){
    for(int i=0;i<i1;i++){
      SetAssembly(i1,in,i+1,frd_sa);
    };
    return;
  };

  int ix,iy;
  GetPositionFromAddress(i1,in,i2,ix,iy);
  tag_map[iy*mx+ix]=frd_sa.GetTag();
};

void FRDCoreLocationManager::SetAssembly(string inp,FRDSubAssembly &frd_sa)
{
  int i1;
  string in;
  int i2;
  AddressStringTransformer(inp,i1,in,i2);
  SetAssembly(i1,in,i2,frd_sa);
};

void FRDCoreLocationManager::AddressStringTransformer(string inp,int &i1,string &in,int &i2)
{
  int sz=inp.size();

  string name[]={"A","B","C","D","E","F","*"};
  int nam_pos=-1;
  for(int i=0;i<sz;i++){
    string tmp=inp.substr(i,1);
    for(int j=0;j<7;j++){
      if(tmp==name[j])nam_pos=i;
    };
  };
  if(nam_pos==-1){
    cout<<"Error in FRDCoreLocationManager::AddressStringTransformer.\n";
    cout<<"This address does not include A-F or *\n";
    exit(0);
  };

  i1=StringToInt(inp.substr(0,nam_pos));
  in=inp.substr(nam_pos,1);
  i2=StringToInt(inp.substr(nam_pos+1,sz-nam_pos-1));
};

void FRDCoreLocationManager::ShowMap()
{
  int *in=new int[mx*my];
  CreateIntmap(in);
  for(int y=0;y<my;y++){
    if(y%2==0)cout<<" ";
    for(int x=0;x<mx;x++){
      int ii=y*mx+x;
      cout<<in[ii]<<" ";
    };
    cout<<"\n";
  };
  delete [] in;
};

int FRDCoreLocationManager::GetTagNumber()
{
  vector<string> rtag;
  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      string ntag=GetTagmap(x,y);
      if(ntag!=""){
	int num=rtag.size();
	bool exist=false;
	for(int i=0;i<num;i++){
	  if(ntag==rtag[i])exist=true;
	};
	if(!exist){
	  rtag.push_back(ntag);
	};
      };
    };
  };
  return rtag.size();
};

bool FRDCoreLocationManager::RepresentativeTag(int ix,int iy)
{
  string rtag=GetTagmap(ix,iy);
  for(int i=0;i<iy*mx+ix;i++){
    if(tag_map[i]==rtag)return false;
  };
  return true;
};

string FRDCoreLocationManager::GetTagmap(int i1,string in,int i2)
{
  int ix,iy;
  GetPositionFromAddress(i1,in,i2,ix,iy);
  return tag_map[iy*mx+ix];
};

void FRDCoreLocationManager::CreateIntmap(int *in)
{
  int index=0;
  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      int tmp=y*mx+x;
      string ntag=tag_map[tmp];
      if(ntag==""){in[tmp]=-1;}
      else if(RepresentativeTag(x,y)){
	for(int i=0;i<my*mx;i++){
	  if(tag_map[i]==ntag)in[i]=index;
	};
	index++;
      };
    };
  };

};

bool FRDCoreLocationManager::LoadedAssembly(string taginp)
{
  for(int i=0;i<mx*my;i++){
    if(tag_map[i]==taginp)return true;
  };
  return false;
};

void FRDCoreLocationManager::LoadedAssembly(string taginp,int &x,int &y)
{
  for(int i=0;i<my;i++){
    for(int j=0;j<mx;j++){
      if(tag_map[i*mx+j]==taginp){
	x=j;
	y=i;
      };
    };
  };
};

void FRDCoreLocationManager::ShowLoadedAssembly()
{
  vector<string> loaded_tag;
  vector<int> count_sa;
  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      string tmp=tag_map[y*mx+x];
      int sz=tmp.size();
      int fpos=sz;
      for(int i=1;i<sz;i++){
	if(tmp.substr(i,1)=="*")fpos=i;
      };
      string tmp2=tmp.substr(0,fpos);
      bool exist=false;
      int csz=loaded_tag.size();
      for(int i=0;i<csz;i++){
	if(tmp2==loaded_tag[i]){
          count_sa[i]++;
	  exist=true;
	};
      };
      if(!exist){
	loaded_tag.push_back(tmp2);
	count_sa.push_back(1);
      };
    };
  };

  cout<<"# * Loaded sub-assembly *\n";
  int csz=loaded_tag.size();
  for(int i=1;i<csz;i++){
    cout<<"#    "<<loaded_tag[i]<<" : "<<count_sa[i]<<"\n";
  };
  cout<<"\n";
};

void FRDCoreLocationManager::ShowLoadedAssemblyMap()
{
  vector<string> loaded_tag;
  vector<int> map_as(mx*my);

  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      string tmp=tag_map[y*mx+x];
      int sz=tmp.size();
      int fpos=sz;
      for(int i=1;i<sz;i++){
	if(tmp.substr(i,1)=="*")fpos=i;
      };
      string tmp2=tmp.substr(0,fpos);
      int exist=-1;
      int csz=loaded_tag.size();
      for(int i=0;i<csz;i++){
	if(tmp2==loaded_tag[i]){
	  exist=i;
	};
      };
      if(exist==-1){
	map_as[y*mx+x]=csz;
	loaded_tag.push_back(tmp2);
      }else{
	map_as[y*mx+x]=exist;
      };
    };
  };

  string ppp[]={
    "A","B","C","D","E", "F","G","H","I","J",
    "K","L","M","N","O", "P","Q","R","S","T",
    "U","V","W","X","Y", "Z"
  };

  for(int y=0;y<my;y++){
    if(y%2==0)cout<<" ";
    for(int x=0;x<mx;x++){
      int ii=y*mx+x;
      if(map_as[ii]<10){
        cout<<map_as[ii]<<" ";
      }else{
        cout<<ppp[map_as[ii]-10]<<" ";
      };
    };
    cout<<"\n";
  };
};

void FRDCoreLocationManager::WriteFileTagMap(string mdir,string filename)
{
  mdir.append(filename);

  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(0);
  };

  int index=0;
  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      if(tag_map[index]==""){
        fout<<"**EMPTY**\n";
      }else{
        fout<<tag_map[index]<<"\n";
      };
      index++;
    };
  };

  fout.close();
};

void FRDCoreLocationManager::ReadFileTagMap(string mdir,string filename)
{
  mdir.append(filename);

  ifstream fin;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(0);
  };

  int index=0;
  for(int y=0;y<my;y++){
    for(int x=0;x<mx;x++){
      fin>>tag_map[index];
      if(tag_map[index]=="**EMPTY**")tag_map[index]="";
      index++;
    };
  };

  fin.close();
};

// +++ FRDFuelLoadingData +++

int FRDFuelLoadingData::GetIDFromTag(string in)
{
  int sz=tag.size();
  for(int i=0;i<sz;i++){
    if(in==tag[i])return i;
  };
  return -1;
};

void FRDFuelLoadingData::PutData(string itag,int ii,string *in)
{
  if(GetIDFromTag(itag)!=-1){
    cout<<"Error in FRDFuelLoadingData::PutData.\n";
    cout<<"Tag name "<<itag<<" already exists.\n";
    exit(0);
  };

  tag.push_back(itag);
  vector<string> inp(ii);
  for(int i=0;i<ii;i++){
    inp[i]=in[i];
    for(int j=i+1;j<ii;j++){
      if(in[i]==in[j]){
	cout<<"Error in FRDFuelLoadingData::PutData.\n";
	cout<<"Same address "<<in[i]<<" is used.\n";
	exit(0);
      };
    };
  };
  address.push_back(inp);
};

int FRDFuelLoadingData::GetNumber(string itag)
{
  int ii=GetIDFromTag(itag);
  if(ii==-1){
    cout<<"Error in FRDFuelLoadingData::GetNumber.\n";
    cout<<"Tag name "<<itag<<" does not exist.\n";
    exit(0);
  };
  return address[ii].size();
};

string FRDFuelLoadingData::GetAddress(string itag,int ij)
{
  int ii=GetIDFromTag(itag);
  if(ii==-1){
    cout<<"Error in FRDFuelLoadingData::GetAddress.\n";
    cout<<"Tag name "<<itag<<" does not exist.\n";
    exit(0);
  };
  return address[ii][ij];
};

bool FRDFuelLoadingData::CheckAddress()
{
  int sz=tag.size();
  for(int i=0;i<sz;i++){
    for(int j=i+1;j<sz;j++){
      int sz1=address[i].size();
      for(int k=0;k<sz1;k++){
	string nam1=address[i][k];
	int sz2=address[j].size();
	for(int l=0;l<sz2;l++){
	  if(address[j][l]==nam1)return false;
	};
      };
    };
  };
  return true;
};

void FRDFuelLoadingData::ShowLoadingMap(int naminp,string *name,FRDCoreLocationManager &clm)
{
  int mx=clm.GetMX();
  int my=clm.GetMY();
  vector<int> in(mx*my,0);

  for(int i=0;i<naminp;i++){
    int id=GetIDFromTag(name[i]);
    if(id==-1){
      cout<<"++ Warning in FRDFuelLoadingData::ShowLoadingMap ++\n";
      cout<<" The tag "<<name[i]<<" is not found in FRDFuelLoadingData.\n";
    }else{
      for(int j=0;j<GetNumber(id);j++){
	string ads=GetAddress(id,j);
        int xx,yy;
        clm.GetPositionFromAddress(ads,xx,yy);
        if(in[yy*mx+xx]!=0){
	  cout<<"Error in FRDFuelLoadingData::ShowLoadingMap.\n";
	  cout<<"The address "<<ads<<" is included in the two patterns.\n";
	  exit(0);
	};
	in[yy*mx+xx]=1;
      };
    };
  };

  for(int y=0;y<my;y++){
    if(y%2==0)cout<<" ";
    for(int x=0;x<mx;x++){
      int ii=y*mx+x;
      cout<<in[ii]<<" ";
    };
    cout<<"\n";
  };
};
