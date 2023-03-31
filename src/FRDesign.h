#ifndef FRDESIGN
#define FRDESIGN

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "Numeric.h"
#include "FRDesignTool.h"
#include "Medium.h"
#include "CartCore.h"
#include "Burnup.h"

using namespace std;

class FRDSubAssembly{
 protected:
  string tag;
  int n_ax_zone;
  vector<int> matn_ax_zone;
  vector<string> tag_ax_zone;
  vector<int> medid_ax_zone; // Medium ID in loaded core to know the burnup region from SYSTEM
  vector< vector<int> > mat_ax_zone;
  vector< vector<real> > den_ax_zone;
  vector<real> len_ax_zone;
  real area_xy; 
  vector<real> heavy_metal_weight; // unit : g
  vector<real> generated_power; // unit : Wd
  vector<real> fast_flux_fluence; // unit : nvt
  real irradiated_day;
 public:
  FRDSubAssembly();
  void PutNumberAxialZone(int i);
  void PutAxialLength(real *i);
  void PutAxialTag(string *i);
  void PutAxialLength(int i,real j){len_ax_zone[i]=j;};
  void PutAxialTag(int i,string j){tag_ax_zone[i]=j;};
  void PutTag(string i){tag=i;};
  void PutMatnum(int i,int j);
  void PutMediumID(int i,int j);
  void PutDensityData(int iax,int inum,int *id,real *den);
  void PutDensityData(int i, FRDTSubAssembly& frdt_sa);
  void PutMatID(int i,int j,int k){mat_ax_zone[i][j]=k;};
  void PutDensity(int i,int j,real k){den_ax_zone[i][j]=k;};
  void ShowDensityData(int i);
  void PutAreaXY(real i){area_xy=i;};
  void AddGeneratedPower(int z,real p){generated_power[z]+=p;};
  void AddFastFluxFluence(int z,real f){fast_flux_fluence[z]+=f;};
  void AddIrradiatedDay(real d){irradiated_day+=d;};
  
  string GetTag(){return tag;};
  int GetNumberAxialZone(){return n_ax_zone;};
  int GetMatnum(int i){return matn_ax_zone[i];};
  int GetMatID(int i,int j){return mat_ax_zone[i][j];};
  int GetMediumID(int i){return medid_ax_zone[i];};
  real GetDensity(int i,int j){return den_ax_zone[i][j];};
  string GetAxialTag(int i){return tag_ax_zone[i];};
  real GetAxialLength(int i){return len_ax_zone[i];};
  void PutDataToAssembly(Assembly &asmi);

  void PutDensityDataToMedium(Medium &med,int z);
  void GetDensityDataFromMedium(Medium &med,int z);
  void PutDensityDataToBurnup(Burnup &bu,int z);
  void GetDensityDataFromBurnup(Burnup &bu,int z);
  real GetIrradiatedDay(){return irradiated_day;};
  void PutIrradiatedDay(real i){irradiated_day=i;};

  real GetHeavyMetalWeight(int i){return heavy_metal_weight[i];}; // unit : g
  real GetHeavyMetalWeight(string itag="");
  void PutHeavyMetalWeight(int i,real j){heavy_metal_weight[i]=j;};
  real GetGeneratedPower(string itag="");
  real GetGeneratedPower(int i){return generated_power[i];};
  void PutGeneratedPower(int i,real j){generated_power[i]=j;};
  real GetFastFluxFluence(string itag="");
  void PutFastFluxFluence(int i,real j){fast_flux_fluence[i]=j;};
  real GetMaximumFastFluxFluence();

  void CalHeavyMetalWeight(Burnup &bu);

  void WriteFile(string mdir,string ss);
  void ReadFile(string mdir,string ss);
};

class FRDMediumSet{
 protected:
  int grp,pl,plt;
  int mednum;
  vector<Medium> med;
  vector<string> medtag;
 public:
  FRDMediumSet(int igrp,int ipl,int iplt=-1);
  void PutData(int n,int *id,real *den,string itag,real temp);
  void PutData(int n,int *id,string itag);
  void PutDensityData(FRDTSubAssembly &frdt_sa,string itag);
  int GetIDFromTag(string itag);
  int GetMednum(){return mednum;};
  int GetGroup(){return grp;};
  int GetPL(){return pl;};
  int GetPLT(){return plt;};
  Medium &GetMedium(int i){return med[i];};
  Medium &GetMedium(string itag){return med[GetIDFromTag(itag)];};
  void GroupCollapsing(int ngrp,int *bgrp,GroupData1D &wgt1,GroupData1D &wgt2);
  void GroupCollapsing(string mtag,int ngrp,int *bgrp,GroupData1D &wgt1,GroupData1D &wgt2);
  string GetMedtag(int i){return medtag[i];};
  void PutMedtag(int i,string inp){medtag[i]=inp;};
};

class FRDCoreLocationManager{
 protected:
  int mx,my;
  int cx,cy;
  vector<string> tag_map;
 public:
  FRDCoreLocationManager(int ix,int iy);
  void SetCenterPosition(int ix,int iy);
  void GetPositionFromAddress(int i1,string in,int i2,int& ix,int& iy);
  void GetPositionFromAddress(string inp,int& ix,int& iy);
  void SetAssembly(int i1,string in,int i2,FRDSubAssembly &frd_sa);
  void SetAssembly(string inp,FRDSubAssembly &frd_sa);
  void AddressStringTransformer(string inp,int &i1,string &in,int &i2);
  void ShowMap();
  int GetMX(){return mx;};
  int GetMY(){return my;};
  string GetTagmap(int x,int y){return tag_map[y*mx+x];};
  string GetTagmap(int i1,string in,int i2);
  void CreateIntmap(int *in);
  int GetTagNumber();
  bool RepresentativeTag(int x,int y);
  bool LoadedAssembly(string taginp);
  void LoadedAssembly(string taginp,int &x,int &y);
  void ShowLoadedAssembly();
  void ShowLoadedAssemblyMap();
  void WriteFileTagMap(string mdir,string filename);
  void ReadFileTagMap(string mdir,string filename);
};

class FRDFuelLoadingData{
 protected:
  vector<string> tag;
  vector< vector<string> > address;
 public:
  FRDFuelLoadingData(){};
  int GetIDFromTag(string in);
  void PutData(string itag,int i,string *in);
  int GetNumber(string itag);
  int GetNumber(int i){return address[i].size();};
  string GetAddress(string itag,int ii);
  string GetAddress(int i,int i1){return address[i][i1];};
  bool CheckAddress();
  void ShowLoadingMap(int naminp,string *name,FRDCoreLocationManager &clm);
};
#endif
