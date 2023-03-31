//
// BURNUP CHAIN revised in 2014/8/7 (See the notebook)
//
//     - Decay data for heavy nuclides are taken from ENDF/B-VII.1.
//     - Am-241(n,g) and Np-237(n,2n) branching ratios 
//       are taken from JENDL-4.0.
//
// ! NOTE !
//  
//   -  "SetDefault" chain is specific for fast reactor analyses.

#ifndef BURNUPCHAIN
#define BURNUPCHAIN

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "GroupData.h"
#include "GeneralSystem.h"
#include "MATIDTranslator.h"
#include "BurnupChainGenerator.h"

using namespace std;

class NuclideChainData{
 private:
  int id; // nuclide ID
  vector<int> ndiv; // number of yielded nuclides
  // fission, capture, n2n, decay
  vector< vector<real> > ratio; // ratio of each yielded nuclide
  vector< vector<int> > idnext; // ID of yielded nuclide
  vector<real> decay_energy; // (eV) [elp,eem,ehp]
 public:
  NuclideChainData();
  void PutID(int i);
  int GetID(){return id;};
  void PutNdiv(int i,int j);
  void PutNdiv(int *in);
  void PutNdiv(int i1,int i2,int i3,int i4);
  void PutNdivFission(int i){PutNdiv(0,i);};
  // fission,capture,n2n,decay
  void PutData(int in,real *rin,int *idin);
  void PutData(int in,real rin,int idin);
  void PutIDnextDecay(int id){PutData(3,1.,id);};
  void PutIDnextCapture(int id){PutData(1,1.,id);};
  void PutIDnextN2N(int id){PutData(2,1.,id);};
  void PutIDnextFission(int i,int id,real w){ratio[0][i]=w; idnext[0][i]=id;};
  void PutRatio(int i,int j,real val);
  int GetNdiv(int i){return ndiv[i];};
  real GetRatio(int i,int j){return ratio[i][j];};
  int GetIDnext(int i,int j){return idnext[i][j];};  
  void ShowFissionYieldData(MATIDTranslator &midt);
  void PutDecayEnergy(int i,real j){decay_energy[i]=j;};
  real GetDecayEnergy(int i){return decay_energy[i];};
  void ShowSelf(MATIDTranslator &midt);
  void SetPseudoFP(int fpid);
  real GetFissionYield(int fpid);
  void PutFissionYield(int fpid, real yld);
  void AddFissionYield(int fpid, real yld);
  void SetSpontaneousFission(real br);
  void SetSpontaneousFission(real br,NuclideChainData &ncd2);
};

class BurnupChain{
 private:
  map<int,NuclideChainData> data;
  map<int,real> decay_const;
  map<int,real> dosecoef_ingestion;
  map<int,real> dosecoef_inhalation;
  MATIDTranslator midt;
  bool read_decay_constant;
 public:
  BurnupChain();
  void SetDefault();
  void SetDefaultPlusNp239();
  void SetOKUMURAChainFine(bool fr=false);
  void Set21HeavyMetalChain(bool fr=false);
  void Set28HeavyMetalChain(bool fr=false);
  void Set37HeavyMetalChain(bool fr=false);
  void AddThoriumDecaySeries(); // (partial series from Th-228 to Pb-208)
  void AddNeptuniumDecaySeries(); // (partial series from Th-229 to Bi-209)
  void AddUraniumDecaySeries(); // (partial series from Th-230 to Rn-222)
  void AddActiniumDecaySeries(); // (partial series from Th-231 to Rn-219)
  void SetHeavyMetalChainForADS(); 
  void AddHigherCmPlusBk(); 
  void AddPseudoFP();
  void AddPseudoFPMolybdenum();
  void AddPseudoFP(int fissile_num, int *matid_fp, real *yield_fp);
  void AddCobalt59Activation();
  int GetNdiv(int id,int i){return data[id].GetNdiv(i);};
  int GetNextID(int id,int ic,int i=0){return data[id].GetIDnext(ic,i);};
  real GetRatio(int id,int ic,int i=0){return data[id].GetRatio(ic,i);};
  int GetNdivFission(int id){return data[id].GetNdiv(0);};
  int GetNdivCapture(int id){return data[id].GetNdiv(1);};
  int GetNdivN2N(int id){return data[id].GetNdiv(2);};
  int GetNdivDecay(int id){return data[id].GetNdiv(3);};
  real GetRatioFission(int id,int i=0){return data[id].GetRatio(0,i);};
  int GetNextIDFission(int id,int i=0){return data[id].GetIDnext(0,i);};
  real GetRatioCapture(int id,int i=0){return data[id].GetRatio(1,i);};
  int GetNextIDCapture(int id,int i=0){return data[id].GetIDnext(1,i);};
  real GetRatioN2N(int id,int i=0){return data[id].GetRatio(2,i);};
  int GetNextIDN2N(int id,int i=0){return data[id].GetIDnext(2,i);};
  real GetRatioDecay(int id,int i=0){return data[id].GetRatio(3,i);};
  int GetNextIDDecay(int id,int i=0){return data[id].GetIDnext(3,i);};
  //
  void SetZeroDecayConstant(string nucnam);
  void FactorizeHalfLife(string nucnam,real factor);
  void SetZeroFissionYield(string nuc_fis,string nuc_fp);
  void CutChain(string nucnam,int i);
  void CutChainForCapture(string nucnam){CutChain(nucnam,1);};
  void CutChainForN2N(string nucnam){CutChain(nucnam,2);};
  void CutChainForDecay(string nucnam){CutChain(nucnam,3);};
  // 
  real GetDecayConstant(int id){return decay_const[id];};
  void PutDecayConstant(int id, real dc){decay_const[id]=dc;};
  real GetDoseCoefficientIngestion(int id){return dosecoef_ingestion[id];};
  real GetDoseCoefficientInhalation(int id){return dosecoef_inhalation[id];};
  real GetDecayEnergy(int id,int j){return data[id].GetDecayEnergy(j);};
  real GetFissionYield(int fid,int fpid);
  // (direct access to NuclideChainData)
  void PutNdiv(string nuc,int i1,int i2,int i3,int i4);
  void PutIDnextDecay(string nuc, string nuc2);
  void PutIDnextCapture(string nuc,string nuc2);
  void PutIDnextN2N(string nuc,string nuc2);
  void PutData(string nuc,int in,real *rin,string *idin);
  void PutData(int id);
  // (branching ratio treatment)
  void PutNGBranchingRatio(string nuc,int ch,real val){GetNuclideChainData(nuc).PutRatio(1,ch,val);};
  real GetNGBranchingRatio(string nuc,int ch){return GetNuclideChainData(nuc).GetRatio(1,ch);};
  // (MATIDTranslator)
  int ID(string name){return midt.ID(name);};
  string Name(int id){return midt.Name(id);};
  void OverWritingChainData(string cbglibdir,string filename);
  void OverWritingChainData(string filename);
  void OverWritingChainData(BCGManager &bm);
  void ReadDoseCoefficientData(string cbglibdir);
  // (FP yield input for 193 chain)
  void ReadFPYieldDataFromFile(string cbglibdir,string filename);
  void ReadFPYieldDataFromFile(string filename);
  void RetreiveFPYieldData(BCGManager &bm);
  void ChangeFPYieldData(string nuc1,string nuc2); // copy nuc1 <- nuc2  
  void ShowFPYield(int id);
  void ShowFPYield(string name){ShowFPYield(midt.ID(name));};
  void ShowFPYieldForFP(string name);
  // (Decay constant)
  void WriteDecayConstant();
  void ReadDecayConstantFromFile(string cbglibdir,string filename);
  void ShowHalfLife();
  void ShowDoseCoefficient();
  bool DecayConstantIsRead(){return read_decay_constant;};
  NuclideChainData& GetNuclideChainData(int id){return data[id];};
  NuclideChainData& GetNuclideChainData(string name){return GetNuclideChainData(ID(name));};
  map<int,NuclideChainData>& GetNuclideChainData(){return data;};
  void ShowChainData(string name);
  void ShowChainData(int id){ShowChainData(Name(id));};
  void ShowChainData();
  void ShowChainDataToOrigin(string name);
  void ShowNGBranchingRatio();
  bool IsConnected(int mat1,int mat2); // mat1:parent, mat2:daughter
  void SetPseudoFP(int mat,int fpid);
  void ModifyDecayConstantForSpontaneousFission();
  void SetZeroDecayEnergy();
  void PutDecayEnergy(string name,real val);
  void PutDecayEnergy(int num,string *name,real *val);
  void DataPerturbation(string mdir, string filename);
  void AddPseudoFPDecay(int fissile_num, int *matid_fp, real *yield_fp, real dc);
  void AddPseudoFPDecay(int fissile_num, int *matid_fp, real *yield_fp, int *matid_fp2, real dc);
  void AddPseudoFPDecay2(int fissile_num, int *matid_fp, real *yield_fp, int *matid_fp2, real *yield_fp2, real dc);
  void AddPseudoFPDecay3(int fissile_num, int *matid_fp, real *yield_fp, int *matid_fp2, real *yield_fp2, int *matid_fp3, real *yield_fp3, real dc, real dc2);
  void AddPseudoFPDecay5(int fissile_num, int *matid_fp, real *yield_fp, int *matid_fp2, real *yield_fp2, int *matid_fp3, real *yield_fp3, int *matid_fp4, real *yield_fp4, int *matid_fp5, real *yield_fp5, real dc, real dc2, real dc3, real dc4);
  void AddPseudoFPDecay5(int matid_fp, real yield_fp, int matid_fp2, real yield_fp2, int matid_fp3, real yield_fp3, int matid_fp4, real yield_fp4, int matid_fp5, real yield_fp5, real dc, real dc2, real dc3, real dc4);
  void AddPseudoFPD5(int fissile_num, int *matid_fp, real *yield_fp, real *dc);
  
};

#endif


