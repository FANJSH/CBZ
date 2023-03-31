#ifndef LIBDATA
#define LIBDATA

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <fstream>
#include <fstream>
#include <map>
#include <string>
#include "Numeric.h"
#include "GroupData.h"
#include "GroupDataSet.h"
#include "MATIDTranslator.h"

const int MAX_RID=5;

using namespace std;

class LibDataPTable{
  vector<int> step;
  vector< vector<real> > probability;
  vector< vector<real> > total;
  vector< vector<real> > elastic;
  vector< vector<real> > fission;
  vector< vector<real> > capture;
  vector< vector< vector<real> > > sigs_gg_p;
  vector<real> maxerr_total;
  vector<real> maxerr_elastic;
  vector<real> maxerr_fission;
  vector<real> maxerr_capture;
  vector< vector<real> > downscat_distribution;
 public:
  LibDataPTable();
  void PutGroup(int g);
  void PutProbability(int g,int st,real *p);
  void PutTotal(int g,real *xs);
  void PutElastic(int g,real *xs);
  void PutFission(int g,real *xs);
  void PutCapture(int g,real *xs);
  void ReadFile(string mdir,string ss,real maxerr);
  void WriteFile(string mdir,string ss,int stgrp=0);
  void ReadMatrixFile(string mdir,string ss);
  int GetStep(int g){return step[g];};
  real GetProbability(int g,int st){return probability[g][st];};
  real GetTotal(int g,int st){return total[g][st];};
  real GetElastic(int g,int st){return elastic[g][st];};
  real GetCapture(int g,int st){return capture[g][st];};
  real GetFission(int g,int st){return fission[g][st];};
  real GetSigsGGP(int g,int st1,int st2){return sigs_gg_p[g][st1][st2];};
  real GetDownscatDistribution(int g,int st){return downscat_distribution[g][st];};
  void show_self();
};

class FTable{
  vector< vector< vector<real> > > data;
  // data[r][t][sig0]
 public:
  FTable(){};
  void Initialize(int r,int t,int s);
  void PutData(int r,int t,int s,real i){data[r][t][s]=i;};
  real GetData(int r,int t,int s){return data[r][t][s];};
};

class LibDataFTable{
  int group;
  int nomtft,nsig0,maxtemp,maxnr;
  vector<real> val_sig0;
  vector<real> val_sig0_log;
  vector<real> val_temp;
  vector<real> val_temp_log;
  vector< vector<real> > val_rpara; // [nr][sig0]
  vector< vector<int> > id_rpara;  // [grp][rid];
  vector<int> no_temp; // temperature point per nomtft
  vector<int> no_rpara;
  vector<int> mt_ftab; // mt number per nomtft

  vector<int> startgrp_temp; // (not used?)
  vector<int> startgrp_ftab; 
  vector<int> endgrp_ftab;

  int startgrp_allmf,endgrp_allmf;

  vector<int> num_fdata;
  vector< vector< vector<FTable> > > fdata;
  // fdata[mt][grp][rid];
  // mt=0 : fission
  // mt=1 : capture
  // mt=2 : elastic
  // mt=3 : total (current-weighted)
  // mt=4 : elastic removal
  // mt=5 : inelastic
 public:
  LibDataFTable();
  void Initialize();
  void PutFData(int mt,int g,int rid,int r,int t,int sig0,real val);
  void PutFData(int mt,int g,int rid, FTable sec);
  real GetFData(int mt,int g,int rid,int r,int t,int sig0);
  FTable& GetFData(int i,int j,int k);
  bool CheckInput(int mt,int g,int rid);
  void WriteFtableWithSig0Temp(int mt,int g,int r);
  void WriteFtableWithSig0Rpara(int mt,int g,int t);
  void WriteFtableWithRparaSig0(int mt,int g,int t);
  real GetF(int mt,int g,int rid,real rval,real temp,real sig0);
  real GetFLargeSig0(int mt,int g,real temp);
  real GetF(int mt,int g,int rid,real rval,real r_sigt,real temp,real sig0);
  real GetF(int mt,int g,int rid,int r,int t,int s0){return fdata[mt][g][rid].GetData(r,t,s0);};
  void CopyData(LibDataFTable &sec,int g);
  void AddData(LibDataFTable &sec,int g);
  bool CheckSamePoint(LibDataFTable &sec,int g);
  // for input
  void PutNomtft(int i);
  void PutNsig0(int i);
  void PutMaxtemp(int i);
  void AddMaxtemp(){maxtemp++;};
  void PutMaxnr(int i);
  void ResizeVal_rpara();
  void PutNoTemp(int i,int j){no_temp[i]=j;};
  void PutNoRpara(int i,int j){no_rpara[i]=j;};
  void PutStartGrpTemp(int i,int j){startgrp_temp[i]=j;};
  void PutStartGrpFtab(int i,int j){startgrp_ftab[i]=j;};
  void PutEndGrpFtab(int i,int j){endgrp_ftab[i]=j;};
  void PutMtFtab(int i,int j){mt_ftab[i]=j;};
  void PutValSig0(int i,real j);
  void PutValTemp(int i,real j);
  void PutValRpara(int i,int j,real k){val_rpara[i][j]=k;};
  void PutStartGrpAllMF(int i){startgrp_allmf=i;};
  void PutEndGrpAllMF(int i){endgrp_allmf=i;};
  void PutGroup(int i);
  void PutIDRpara(int g,int nid,int id);
  // for output
  int GetNsig0(){return nsig0;};
  int GetNoTemp(int i){return no_temp[i];};
  int GetNoRpara(int i){return no_rpara[i];};
  int GetStartGrpTemp(int i){return startgrp_temp[i];};
  int GetStartGrpFtab(int i){return startgrp_ftab[i];};
  int GetEndGrpFtab(int i){return endgrp_ftab[i];};
  int GetStartGrpAllMF(){return startgrp_allmf;};
  int GetEndGrpAllMF(){return endgrp_allmf;};
  int GetNomtft(){return nomtft;};
  int GetMaxtemp(){return maxtemp;};
  int GetMaxnr(){return maxnr;};
  int GetMtFtab(int i){return mt_ftab[i];};
  int GetIDRpara(int g,int nid){return id_rpara[g][nid];};
  int GetNumFdata(int g){return num_fdata[g];};
  real GetSig0(int i){return val_sig0[i];};
  real GetTemp(int i){return val_temp[i];};
  real GetRpara(int i,int j){return val_rpara[i][j];};
  void AddValTemp(real i){val_temp.push_back(i);};
  void AddValTempLog(real i){val_temp_log.push_back(i);};
  
  //
  void ShowSelf();
};

class LibDataChiVector{
  int no_vector,group;
  vector<int> vectorid;
  vector<GroupData1D> chiv;
 public:
  LibDataChiVector();
  void PutGrp(int g);
  void PutNoVector(int i);
  void Initialize();
  void PutVectorID(int i,int j);
  void ShowVectorID();
  int GetVectorID(int i){return vectorid[i];};
  int GetNoVector(){return no_vector;};
  int GetGroup(){return group;};
  GroupData1D &GetChiVector(int i){return chiv[i];};
  GroupData1D &GetChiVectorGroup(int i){return chiv[vectorid[i]];};
  void AssignChiVectorToEachGroup();
};

class ThermalScatteringData{
  int grp,init_grp, pl, ntemp;
  int init_src_grp;
  bool exist;
  vector<real> temp;
  vector< vector<GroupData2D> > data;
 public:
  ThermalScatteringData();
  void PutTemperature(int i,vector<real> t);
  void Initialize(int igrp,int initgrp,int ipl);
  void ReadFile(string mdir,string ss);
  void ReadFileSlide(string mdir,string ss,int dg);  
  real GetData(int itemp,int ipl,int srcg,int sinkg);
  real GetData(int ipl,int srcg,int sinkg,real tt);
  int GetInitGrp(){return init_grp;};
  int GetInitSrcGrp(){return init_src_grp;};
  int GetPL(){return pl;};
  int GetGrp(){return grp;};
  bool DoesExist(){return exist;};
  GroupData2D &GetData(int itemp,int ipl){return data[itemp][ipl];};
};

class LibData{
  int mat,group;
  bool exist_fiss;
  bool exist_chi_vector;
  vector<int> maxpl1;
  vector<real> bell_factor;
  bool exist_bell_factor;
  GroupDataSet xsdata;
  LibDataFTable ftable;
  LibDataChiVector chiv;
  LibDataPTable ptable;
  ThermalScatteringData thscat;
 public:
  LibData();
  void PutGroup(int g);
  void ReadFile(string mdir,string ss,MATIDTranslator& midt);
  void WriteFile(string mdir,string ss);
  void WriteFileUNCFormat(string mdir,string ss);
  GroupDataSet &GetXSData(){return xsdata;};
  LibDataFTable &GetFtable(){return ftable;};
  LibDataPTable &GetPtable(){return ptable;};
  ThermalScatteringData &GetThScat(){return thscat;};
  LibDataChiVector &GetChiVector();
  void PutChiVector(LibDataChiVector tmp){chiv=tmp; exist_chi_vector=true;};
  GroupData1D &GetChiVector(int i){return GetChiVector().GetChiVector(i);};
  real GetF(int mt,int g,int rid,real rval,real temp,real sig0){return ftable.GetF(mt,g,rid,rval,temp,sig0);};
  real GetF(int mt,int g,int rid,real rval,real r_sigt,real temp,real sig0){return ftable.GetF(mt,g,rid,rval,r_sigt,temp,sig0);};
  void CalF(int g,int rid,real rval,real r_sigt,real temp,real sig0,vector<real> &f);
  int GetGroup(){return group;};
  void PutExistFiss(bool i){exist_fiss=i;};
  void PutMat(int i){mat=i;};
  void SetSizeMaxpl1(int i){maxpl1.resize(i,0);};
  void PutMaxpl1(int i,int j){maxpl1[i]=j;};
  bool ExistChiVector(){return exist_chi_vector;};
  bool fissile(){return exist_fiss;};
  bool ExistThermalData();
  void PutBellFactor(real *inp);
  void PutBellFactorAllGroup(real inp);
  void PutBellFactor(int g, real inp);    
  bool ExistBellFactor(){return exist_bell_factor;};
  real GetBellFactor(int i){return bell_factor[i];};
  int GetMat(){return mat;};
  void Update1DChiDataForFastReactorCalculation();
  void CrossSectionPerturbation(int g, enum xstype sigx, real r_delta);
  void NormalizeMatrixData(vector< vector<real> > &data);
  void ConstantCrossSectionData(int mat_inp, int g, real sigc_inp);
  void AddFtableAtNewTemperaturePoint(LibData &sec);
  //
  // FRENDY-MG
  void ReadFileFrendyMG(string frendy_dir, string inpname, int matid, int group_inp=175);
  void ReadFileFrendyMG(int tempnum, int* temp_list, string frendy_dir, string inpname_i, int matid, int group_inp=175);
  void ReadFileFrendyMG_NuChi(string filename, vector<real> &data_chi, vector<real> &data_chid, vector<real> &nud);
  void ReadFileFrendyMG_XS1D(string filename, vector< vector<real> > &data, vector<real> &ebnd);
  void ReadFileFrendyMG_XS2D(string filename, vector< vector< vector<real> > > &data, int pl);
  void ReadFileFrendyMG_XS2D_MT4(string filename1, string filename2, vector< vector< vector<real> > > &data, int pl, int max_inelalevel);
  void ReadFileFrendyMG_XS1D_MT4(string filename1, string filename2, vector<real> &data, int max_inelalevel);
  void ReadFileFrendyMG_XS1D_MT102(string filename1, string filename2, vector< vector<real> > &data);
  void TSL(string tslname, string aceid, int mtnum, int *mtlist, string libdir, string libname);  

};

class XSLibrary{
  map<int,LibData> data;
  GroupData1D enband;
  GroupData1D wtflux;
  MATIDTranslator midt;
  bool warning_print; // Warning message for no data in xslibrary for self-shielding calculation
  string filename;
 public:
  XSLibrary();
  XSLibrary(string dir,string ss){Initialize(dir,ss);};
  void Initialize(string dir,string ss){ReadNEnergy(dir,ss);};
  void AddLibData(int mat,LibData &libinp){data[mat]=libinp;};
  bool ExistLibData(int mat);
  LibData& GetLibData(int mat);
  LibDataPTable& GetPtable(int mat){return GetLibData(mat).GetPtable();};
  void ReadFile(int nucnum,string mdir,string *filename,int *matno);
  void ReadFile(int nucnum,string mdir,string *filename);
  void ReadPTFile(string mdir,string filename,int matno,real maxerr=0.1); // [%]
  void ReadPTMatrixFile(string mdir,string filename,int matno);
  void ReadNEnergy(string dir,string ss);
  GroupData1D &GetEnband(){return enband;};
  GroupData1D &GetWtflux(){return wtflux;};
  void PutBellFactor(int mat,real *inp);
  void ReadBellFactor(string mdir,string filename,int mat);
  int GetGroup(){return enband.get_imax()-1;};
  void ShowSelf();
  void ShowCrossSectionData1D(int mat);
  map<int,LibData>& GetData(){return data;};
  void Update1DChiDataForFastReactorCalculation();

  void ConstantCrossSectionData(int mat_inp, int g, real sigc_inp);

  LibData LibraryMixing(int newnucid, int nucn, int *nucid, real *wgt);

  bool WarningPrint(){return warning_print;};
  void WarningPrintTrue(){warning_print=true;};
};

#endif


