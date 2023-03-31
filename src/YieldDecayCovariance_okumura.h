#ifndef YIELDDECAYCOVARIANCE
#define YIELDDECAYCOVARIANCE

#include <fstream>
#include <iostream>
#include <cstdio>

#include "GroupData.h"
#include "ENDFIDTranslator.h"
#include "BurnupChainGenerator.h"

class IndependentYieldCovariance{
 protected:
 public:
  IndependentYieldCovariance(){
  fisnucnum=0;
  theor_err=1.;
};
  ~IndependentYieldCovariance(){};

  vector<GroupData2D> matrix; // might be RELATIVE covariance // !!
  
  vector<string> fisnucname_array; // !!
  vector<string> filename_array;
  vector<int> energy_point_array;//0,1,2,...
  vector<real> energy_array;
  vector<int> size_array;
  int fisnucnum;
  real theor_err;

  vector<vector<real> > mass;
  vector<vector<real> > yield_idp;
  vector<vector<real> > yield_cum;
  vector<vector<real> > unc_idp;
  vector<vector<real> > unc_cum;
  vector<vector<int> > id;//nuclide ID !!


  ENDFIDTranslator eidt;
  MATIDTranslator midt;

  void MakeCovarianceMatrix(BCGManager &bm);
  void MakeCovarianceMatrixWithoutCorrelation();
  void ReadYieldDataFromFile(string fisnucname, string filename, int epnt);
  void ReadDataFromFile(int num);
  void CalCovariance(BCGManager &bm,int num);
  void CalCovarianceWithoutCorrelation(int num);
  
  int GetNumberFromFissionNuclideName(string fisnucname);
  int GetSize(int i){return size_array[i];};
  int GetSize(string fisnucname){return GetSize(GetNumberFromFissionNuclideName(fisnucname));};
  int GetEnergyPoint(int i){return energy_point_array[i];};
  int GetEnergyPoint(string fisnucname){return GetEnergyPoint(GetNumberFromFissionNuclideName(fisnucname));};
  real GetEnergy(int i){return energy_array[i];};
  real GetEnergy(string fisnucname){return GetEnergy(GetNumberFromFissionNuclideName(fisnucname));};
  int GetID(int i,int j){return id[i][j];};
  int GetID(string fisnucname,int j){return GetID(GetNumberFromFissionNuclideName(fisnucname),j);};
  GroupData2D GetCovarianceMatrix(int i){return matrix[i];};
  GroupData2D GetCovarianceMatrix(string fisnucname){return GetCovarianceMatrix(GetNumberFromFissionNuclideName(fisnucname));};
  int GetFissionNuclideNumber(){return fisnucnum;};
  string GetFissionNuclideName(int i){return fisnucname_array[i];};
  int GetFissionNuclideID(int i){return midt.GetMATID(fisnucname_array[i]);};
  void PlotCovarianceMatrix(int i);
  void PlotCovarianceMatrix(string fisnucname){PlotCovarianceMatrix(GetNumberFromFissionNuclideName(fisnucname));};
  void DeleteCovarianceExceptAssignedMassChainFissionNuclide(int mass,int mat);
  real GetIndependentYield(int i,int j){return yield_idp[i][j];};
  real GetIndependentYieldStandardDeviation(int i,int j){return sqrt(matrix[i].get_dat(j,j));};//relative value
  void SetTheoreticalValueError(real e){theor_err=e;};

  void ShowVarianceData(string fisnucname);
  GroupData2D GetPartialMatrix(string fisnucname, int num, int *idlist);

  real GetYieldIDP(string fisnucname, int mat1, int mat2);
  void PutReducedChainNumber(int num); // okumura
  void PutReducedChainData(string fisnucname, int fnuc, int x); // okumura
  //void PutReducedChainData(string fisnucname, int fnuc, int x, real *average, real *meansq);//, real figure[][shortennucnum]); // okumura
  void ReadReducedChainData(const char *wfilen);
};

class HalfLifeCovariance{
 protected:
 public:
  HalfLifeCovariance(){
    theor_err=1.;
  };
  ~HalfLifeCovariance(){};

  GroupData2D matrix;
  int size;
  real theor_err;
  vector<real> halflife;
  vector<real> unc;
  vector<int> id;//nuclide ID
  
  MATIDTranslator midt;
  
  void MakeCovarianceMatrix(BCGManager &bm);
  
  int GetSize(){return size;};
  int GetID(int i){return id[i];};
  GroupData2D GetCovarianceMatrix(){return matrix;};
  void DeleteCovarianceExceptAssignedNuclide(int mat);

  real GetHalfLife(int i){return halflife[i];};
  real GetHalfLifeStandardDeviation(int i){return sqrt(matrix.get_dat(i,i));};//relative value
  void SetTheoreticalValueError(real e){theor_err=e;};

};

class BranchingRatioCovariance{
 protected:
 public:
  BranchingRatioCovariance(){
    theor_err=1.;
  };
  ~BranchingRatioCovariance(){};

  GroupData2D matrix; // RELATIVE covariance matrix
  int size;
  int nucnum;
  real theor_err;
  vector<vector<real> > ratio; // branching ratio [nuc][ch]
  vector<vector<real> > unc;   // DELTA-branching ratio [nuc][ch]
  vector<int> id;//nuclide ID
  vector<int> channel_num;
  vector<int> list_id;      // ID list corresponding to matrix size
  vector<int> list_channel; // channel list corresponding to matrix size

  MATIDTranslator midt;
  
  void MakeCovarianceMatrix(BCGManager &bm);
  
  int GetNucNum(){return nucnum;};
  int GetSize(){return size;};
  int GetSizeIDarray(){return id.size();};
  int GetID(int i){return id[i];};
  int GetIDFromList(int i){return list_id[i];};//corresponding matrix size
  int GetChannel(int i){return ratio[i].size();};
  int GetChannelFromList(int i){return list_channel[i];};
  real GetRatio(int i,int j){return ratio[i][j];};
  real GetUnc(int i,int j){return unc[i][j];};
  GroupData2D GetCovarianceMatrix(){return matrix;};
  GroupData2D GetCovarianceMatrix(int id);
  void DeleteCovarianceExceptAssignedNuclide(int mat);
  real GetBranchingRatioStandardDeviation(int i){return sqrt(matrix.get_dat(i,i));};//relative value
  void SetTheoreticalValueError(real e){theor_err=e;};

};

class DecayEnergyCovariance{
 protected:
 public:
  DecayEnergyCovariance(){
    theor_err=1.;
  };
  ~DecayEnergyCovariance(){};

  vector<GroupData2D> matrix;
  int size;
  real theor_err;
  vector<vector<real> > decay_energy;
  vector<vector<real> > unc;
  vector<int> id;//nuclide ID
  
  MATIDTranslator midt;
  
  void MakeCovarianceMatrix(BCGManager &bm);
  
  int GetSize(){return size;};
  int GetID(int i){return id[i];};
  GroupData2D GetCovarianceMatrix(int i){return matrix[i];};
  void DeleteCovarianceExceptAssignedNuclideDecayType(int mat,int type);
  void DeleteCovarianceExceptAssignedNuclideDecayType(int mat,string type);

  real GetDecayEnergy(int i,int type){return decay_energy[i][type];};
  real GetDecayEnergyStandardDeviation(int i,int type){return sqrt(matrix[type].get_dat(i,i));};//relative value
  void SetTheoreticalValueError(real e){theor_err=e;};

};

#endif
