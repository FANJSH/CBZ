#ifndef MEDIUM
#define MEDIUM

#include "GroupData.h"
#include "GroupDataSet.h"
#include "Nuclide.h"
#include "MATIDTranslator.h"

#include<string>
#include<iostream>
#include<fstream>

using namespace std;

class Medium{
 protected:
  GroupData1D enband; // Energy boundary
  GroupDataSet Macxs; // Macroscopic cross section
  vector<GroupData1D> flux; // Angular flux (0 and 1)
  int imax,pl;
  int pltot;  // Pl order for total cross section
  int nucnum; // Number of nuclide
  vector<Nuclide> nuc;
  void CalHomoB1(real bsq,GroupData1D src);

  vector<int> LowestDownScatGrp;
  int HighestUpScatGrp;
  int fissionable_flag; // 2:yes, 1:no, 0:not-checked
 public:
  Medium(){Init();};
  ~Medium();
  Medium(int i){Init(); PutImax(i);};
  void Init();
  void InitSigST();
  void PutImax(int i);
  void PutPL(int i,int itot=-1);
  int GetImax(){return imax;};
  int GetPL(){return pl;};
  int GetPLT(){return pltot;};

  void PutMacxs(GroupDataSet macinp);
  GroupDataSet &GetMacxs(){return Macxs;};
  void PutData(enum xstype ss,real *i){Macxs.GetData1d(ss,0).put_data(i);};
  void PutDataSigt(int l,real *i){Macxs.GetData1d(sigt,l).put_data(i);};
  void PutDataSigs(int l,real *i){Macxs.GetData2d(sigs,l).put_data(i);};
  void PutData(enum xstype ss,vector<real> i);
  void PutDataSigt(int l,vector<real> i);
  void PutDataSigs(int l,vector<real> i);

  GroupData1D &GetData1D(enum xstype ss,int l=0){return Macxs.GetData1d(ss,l);};
  GroupData2D &GetSigs(int ipl){return Macxs.GetData2d(sigs,ipl);};
  GroupData1D &GetSigt(int ipl){return Macxs.GetData1d(sigt,ipl);};
  GroupData1D &GetEnband(){return enband;};
  GroupData1D &GetFlux(int i=0){return flux[i];};
  real GetData1D(int i,enum xstype ss){return GetData1D(ss).get_dat(i);};
  real GetDataSigs(int m,int i,int j){return Macxs.GetData2d(sigs,m).get_dat(i,j);};
  real GetDataSigt(int m,int i){return Macxs.GetData1d(sigt,m).get_dat(i);}
  void CalRemovalCrossSection();
  void SetZeroChi();

  void ReadPDSFile(string mdir,string ss,int plt=1,bool upscat=false);
  void ReadFile(string mdir,string ss,int plt=1);
  void ReadKramxsFile(string mdir,string ss,bool up_scat=false);  
  void ReadPDSFileUpScat(string mdir,string ss,int plt=1){ReadPDSFile(mdir,ss,plt,true);};
  void WritePDSFile(string mdir,string ss);
  void WriteFile(string mdir,string ss,bool mic=false);
  void WriteFileNumberDensity(string mdir,string ss);
  void ReadFileNumberDensity(string mdir,string ss,bool iso=false);

  int GetUpScatteringGrp(int plmax=0);
  int GetUpScatteringSinkGrp(int plmax=0);
  void PutUpScatteringToSelfScattering();
  void PutUpScatteringToSelfScattering(int g0, int g1);  

  void CalHomoB1(real bsq){CalHomoB1(bsq,GetData1D(chi));};
  void CalHomoB1WithFixedSource(real bsq,GroupData1D &src);
  void CalHomoB1Adjoint(real bsq){CalHomoB1Adjoint(bsq,GetData1D(nusigf));};
  void CalHomoB1Adjoint(real bsq,GroupData1D src);
  void CalHomoB1FixedSource(real bsq,GroupData1D src){CalHomoB1(bsq,src);};
  real CalKeff();
  real CalKeffAdjoint();
  real CalKinf();
  real BucklingSearch();
  void AddPseudoAbsorption(real bsq);

  void CalSigt();
  void CalSigtr(int sigt_l=1);
  void CalSigtrDiagonal(int sigt_l=1);
  void PutSigt0AsSigt1(int g0=0, int g1=-1);
  void CalDFromSigt(int i=0);
  void CalDFromSigtr();
  void CalSigtrFromDav();
  void CalSigtFromDav();
  void TransportApproximation();
  void ConsistentPApproximation();
  void AdjustSelfScattering();

  Medium Cond(int ngrp,vector<int> bgrp);
  Medium Cond(int ngrp,int *bgrp);
  Medium Cond(int ngrp,vector<int> bgrp,GroupData1D f,GroupData1D c,bool micro=true);
  Medium Cond(int ngrp,int *bgrp,GroupData1D f,GroupData1D c,bool micro=true);  

  Nuclide &GetNuclide(int matnum);
  real GetNuclideDensity(int matnum);
  real GetNuclideDensity(int num, vector<int> &matid);
  real GetNuclideDensity(int num, vector<string> &matname, MATIDTranslator &midt);
  Nuclide &GetNuclideInTurn(int i){return nuc[i];};
  int GetNuclideID(int i){return nuc[i].GetMatnum();};
  bool ExistNuclide(int matnum);
  int SearchNuclide(int matnum);
  void PutNucnum(int i);
  int GetNucnum(){return nucnum;};
  void PutNuclide(int i,Nuclide ninp){nuc[i]=ninp;};
  void PutNuclide(int i,int mat,real den);
  void AddNuclide(Nuclide ninp){nuc.push_back(ninp); nucnum++;};
  void PutNuclide(int nuc,int *matno);
  void PutNuclide(int nuc,vector<int> matno);
  void AddNuclide(int nuc,int *matno);
  void PutNuclide(int nuc,int *matno,real *density);
  void PutNuclide(int nuc,int *matno,real *density,MATIDTranslator &midt);
  void PutNuclide(string file);
  void PutNuclideNew(int nuc,int *matno,real *density,bool each_iso=false);
  void PutNuclideFe(int nuc,int *matno,real *density);  
  void PutNuclide(int nuc,vector<int> matno,vector<real> density);
  void PutTemperatureForAllNuclide(real i);

  void ShowNuclideList();
  void ShowNumberDensity(bool zero_ignore=true);
  void ShowNumberDensityCBZStyle(MATIDTranslator &midt);
  void ShowMacroSection(){Macxs.ShowSelf();};
  void ShowMacroXS1D(){ShowMacroCrossSection1D();};
  void ShowMacroCrossSection1D();
  void ShowMicroCrossSection1D(int mat);
  void ShowNeutronFlux();
  void PrintMacroSectionTable();

  void PutFictitiousVacuum(int grp,int ipl=1,real xs=1e-6);
  void PutLowestDownScatGrp();
  void PutHighestUpScatGrp();
  int GetLowestDownScatGrp(int i){return LowestDownScatGrp[i];};
  int GetHighestUpScatGrp(){return HighestUpScatGrp;};

  void FissionSpectrumVectorReconstruction();
  void CalMacroFromMicro(bool chi_calc=true);
  void AddMacroFromMicro(int id);  
  void CalMacroFromMicroSimple();
  GroupData1D GetMacroSigf();
  GroupData1D GetMacroSigc(); 
  GroupData2D GetMacro2D(enum xstype ss);
  void MacxsVectorClear(){Macxs.AllVectorClear();};
  void MicxsVectorClear2DData();
  void MicxsVectorClear();
  //void NuclideClear(){nuc.clear(); nucnum=0;};
  void NuclideClear();
  void SetZeroSelfScattering();
  bool IsFissionable();

  void AdjustNumberDensity();
  void SetZeroNumberDensity();
  void FactorizeNumberDensity(real factor);

  void CalMuFromSigel_p1();
  void CorrectionForLastGroupSigt1();

  void ExtendPLTOrder(int plin);
  void MATIDTranslationFromENDFID(MATIDTranslator &midt);

  GroupData2D CalAmatForAbemie(real keff);

  void CopyWithoutMicroscopicScatteringMatrices(Medium &sec);
  void FourierAnalysisForUpScattering();
};

#endif
