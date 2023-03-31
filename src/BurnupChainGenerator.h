#ifndef BURNUPCHAINGENERATOR
#define BURNUPCHAINGENERATOR

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "MATIDTranslator.h"

using namespace std;

typedef double real;

class BCGNuclide{
 private:
  MATIDTranslator midt;
  int atomic_number;
  int mass_number;
  int ex_level;
  int r_num; // The number of reactions
  real awr;
  real half_life;
  real delta_half_life; // (standard deviation of half life)
  vector<real> decay_energy; // (eV) [LP/EM/HP]
  vector<real> delta_decay_energy; // (standard deviation of decay energy)
  //
  vector<real> yield;
  vector<real> delta_yield; // (standard deviation of fission yield)
  vector<string> yield_tag;

  real n_per_nuc; // number of total emitted neutrons including daughter nuclides
  vector<real> n_per_nuc_ch; // channel-wise

  // (decay)
  int channel;
  vector<int> rr_num; // reaction先の通し番号、KM作
  vector<int> atn_next;
  vector<int> mas_next;
  vector<int> lev_next;
  vector<int> decay_type; // 0:Beta-, 1:EC, 2:IT, 3:Alpha
  vector<int> dn;
  vector<real> br;
  vector<real> delta_br;
  // (reaction) 0:ng / 1:n2n
  vector<int> r_channel;
  vector< vector<int> > r_atn_next;
  vector< vector<int> > r_mas_next;
  vector< vector<int> > r_lev_next;
  vector< vector<real> > r_br;
  //
  bool flag;

  real ten; // KM
 public:
  BCGNuclide(int at,int ms,int lv);
  //
  int GetAtomicNumber(){return atomic_number;};
  int GetMassNumber(){return mass_number;};
  int GetExLevel(){return ex_level;};
  int GetID(){return midt.GetMATID(atomic_number,mass_number,ex_level);};//kawamoto
  int GetChannel(){return channel;};
  real GetHalflife(){return half_life;};
  real GetDeltaHalflife(){return delta_half_life;};
  real GetAWR(){return awr;};
  //
  void PutChannel(int i);
  void PutNextNuclideData(int i,int at,int ms,int lv,real br,real dbr,int type,int dnin);
  void PutHalflife(real i,real j){half_life=i; delta_half_life=j;};
  //
  void PutTotalNumberEmittedNeutrons(real i){n_per_nuc=i;};
  real GetTotalNumberEmittedNeutrons(){return n_per_nuc;};
  void PutTotalNumberEmittedNeutronsCh(int i,vector<real> &inp);
  real GetTotalNumberEmittedNeutronsCh(int i){return n_per_nuc_ch[i];};
  //
  void PutDecayEnergy(int i,real j){decay_energy[i]=j;};
  void PutDeltaDecayEnergy(int i,real j){delta_decay_energy[i]=j;};//kawamoto
  void AddDecayEnergy(int i,real j){decay_energy[i]+=j;};
  real GetDecayEnergy(int i){return decay_energy[i];};
  real GetDeltaDecayEnergy(int i){return delta_decay_energy[i];};
  //
  void PutYield(string tagname,real i,real i2);
  void AddYield(string tagname,real j,real j2);
  real GetYield(string tagname);
  real GetDeltaYield(string tagname);
  int GetYieldID(string tagname);
  bool ZeroYield();
  string GetYieldTag(int i){return yield_tag[i];};
  int GetYieldNum(){return yield.size();};
  //
  void PutAWR(real i){awr=i;};
  //
  void MakeFlagTrue(){flag=true;};
  void MakeFlagFalse(){flag=false;};

  int GetAtomicNumberNext(int i){return atn_next[i];};
  int GetMassNumberNext(int i){return mas_next[i];};
  int GetExLevelNext(int i){return lev_next[i];};
  int GetDecayType(int i){return decay_type[i];};
  int GetEmittedNeutron(int i){return dn[i];};
  real GetBr(int i){return br[i];};
  real GetDeltaBr(int i){return delta_br[i];};
  void PutBr(int i,real j){br[i]=j;};
  void NormalizeBr();

  bool SameNuclide(int atm,int mas,int lev);
  bool Flag(){return flag;};
  //
  void PutDecayData(int i,vector<int> i1,vector<int> i2,vector<int>i3,vector<real>i4,vector<real>i5,vector<real>i6);

  void PutReactionChannel(int i,int j);
  void PutReactionData(int i,int ic,real *br);
  void PutReactionData(int i,real br1,real br2);
  void PutReactionData(int i,real br1,real br2,real br3);
  void PutReactionData(int i,int ic,vector<int>i1,vector<int>i2,vector<int>i3,vector<real>i4);
  int GetReactionChannel(int i){return r_channel[i];};
  int GetReactionAtomicNumberNext(int i,int j){return r_atn_next[i][j];};
  int GetReactionMassNumberNext(int i,int j){return r_mas_next[i][j];};
  int GetReactionExLevelNext(int i,int j){return r_lev_next[i][j];};
  real GetReactionBr(int i,int j){return r_br[i][j];};
  int GetRnum(){return r_num;};
  // Matsuura
  void PutRRnum(int i,int j){rr_num[i]=j;};
  int GetRRnum(int i){return rr_num[i];};
  void PutTEN(real i){ten=i;};
  real GetTEN(){return ten;};
};

class BCGManager{
 private:
  vector<BCGNuclide> nuc;
  vector<string> yield_tagname;
  MATIDTranslator midt;
  void ShowFlagedNuclide(bool flaged,bool name);
 public:
  BCGManager(){};
  void ReadDecayDataFromFile(string filename,bool stable_only=false);
  void WriteDecayDataToFile(string filename);
  void ReadFPYieldDataFromFile(string filename,int pnt,string tagname,bool print=true);
  void WriteFPYieldDataToFile(string filename,string tagname,real eng);
  void ReadNGBranchingRatioDataFromFile(string filename);
  void ReadBranchingRatioDataFromFile(string filename);
  void SetNGBranchingRatioDataForFR();
  int GetNuclideIndex(int atm,int mas,int lev);
  real GetYield(string tagname,int atm,int mas,int lev);
  real GetDeltaYield(string tagname,int atm,int mas,int lev);
  void MakeFlagFalseAllNuclide();
  void MakeFlagTrueAllNuclide();
  void MakeFlagTrue(int atm,int mas,int lev);
  void MakeFlagTrueAtom(int atm);
  void MakeFlagTrue(int nuc,string *nuc_nam);
  void MakeFlagFalse(int nuc,string *nuc_nam);
  void MakeFlagTrue(int nuc,int *matid);
  void MakeFlagFalse(int atm,int mas,int lev){GetNuclide(atm,mas,lev).MakeFlagFalse();};
  void MakeFlagReversed();
  void ShowYieldSum(string tagname);
  void ShowYieldSum();
  void ShowNumberOfTotalEmittedNeutrons();
  real GetYieldSum(string tagname);
  real GetYieldSumFlag(string tagname,bool flag);
  void ShowFlagedNuclideList(bool flaged=true);
  void ShowFlagedNuclideName(bool flaged=true){ShowFlagedNuclide(flaged,true);};
  void ShowFlagedNuclideMatID(bool flaged=true){ShowFlagedNuclide(flaged,false);};
  void ShowFlagedNuclideShortHalflife(real hl,bool flaged=true);
  //
  void CalCumulativeYield(real half_life_limit=0., bool print=true);
  void CalReactionBranch();
  void CalDecayBranch(real half_life_limit=0., bool print=true);
  void CalTotalEmittedNeutrons(bool repeating_cal=false);
  void CalCumulativeDecayEnergy();
  real CalTotalDelayedNeutron(string tagname);
  real CalTotalDelayedNeutronWithCumulativeYield(string tagname);
  void BranchingRatioSumCheck();
  void IgnoreNGReaction(int num,string *nuc_nam);
  BCGNuclide &GetNuclide(int atm,int mas,int lev=0);
  BCGNuclide &GetNuclideFromID(int id);
  BCGNuclide &GetNuclide(int i){return nuc[i];};
  // (kawamoto)
  int GetAtomicNumber(int i){return nuc[i].GetAtomicNumber();};
  int GetMassNumber(int i){return nuc[i].GetMassNumber();};
  int GetExLevel(int i){return nuc[i].GetExLevel();};
  int GetID(int i){return nuc[i].GetID();};
  int GetAtomicNumberNext(int i,int j){return nuc[i].GetAtomicNumberNext(j);};
  int GetMassNumberNext(int i,int j){return nuc[i].GetMassNumberNext(j);};
  int GetExLevelNext(int i,int j){return nuc[i].GetExLevelNext(j);};
  int GetChannel(int i){return nuc[i].GetChannel();};
  real GetHalflife(int i){return nuc[i].GetHalflife();};
  void MakeFlagTrue(int i){return nuc[i].MakeFlagTrue();};
  void MakeFlagFalse(int i){return nuc[i].MakeFlagFalse();};
  bool Flag(int i){return nuc[i].Flag();};
  int GetSize(){return nuc.size();};

  void ShowNGBranch();
  void ShowTotalNumberEmittedNeutrons(string tagname,real low=0.);
  int GetNuclideNumber(){return nuc.size();};

  // (file output)
  void WriteFileChainDataCBGFormat(string filename,bool ng_only=true);
  void WriteFileChainDataCBGFormat(int num,string *nuc_nam,string filename,bool ng_only=true);
  void WriteFileYieldDataCBGFormat(int num,string *nuc_nam,string filename);
  void WriteFileYieldDataCBGFormatCmSF(int num,string *nuc_nam,string filename);
  void WriteFileYieldDataCBGFormatDetail(int num,string *nuc_nam,string filename);
  void WriteFileYieldDataCBGFormatDetail2(int num,string *nuc_nam,string filename);
  void WriteFileYieldDataCBGFormat(string filename,bool detail=false);
  void WriteFileYieldDataCBGFormatStandard(string filename,int ynum,string *yname);
  // (MATIDTranslator function)
  string Name(int an,int ms,int ex){return midt.Name(an,ms,ex);};
  string SearchNameFromParameter(int an,int ms,int ex){return Name(an,ms,ex);};
  // (Method of creating specific burnup chain)
  void SearchShortHalflivedNuclide(real hl_upper,real hl_lower,int fpmat_upper=68,int fpmat_lower=32);
  void SearchShortHalflivedNuclideForSpecificNuclides(real hl_upper,real hl_lower,int num_include_SourceTerm,int *AtomicNumber_include_SourceTerm,int fpmat_upper=68,int fpmat_lower=32);
  // (Output utility)
  void ShowNuclideList(string yld_tag="");
  void ShowYieldData(string yield_tag,bool flaged=true);
  void ShowYieldDataForXYPlot(string yield_tag);
  void ShowDecayDataForXYPlot();
  void ShowDNDataForXYPlot();
  void AddNuclide(int iz,int ia,int il);
  // Matsuura
  void CalTotalNeutrons();//made by km
  void CalRRnums();//made by km
  void CalTEN();//made by km

  void DataPerturbation(string mdir, string filename);
};

#endif
