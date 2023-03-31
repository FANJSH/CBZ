#ifndef UNC_CROSSSECTION
#define UNC_CROSSSECTION

#include <iostream>
#include <time.h>
#include <string>
#include <vector>
#include <math.h>

#include "GroupData.h"
#include "UNC_Covariance.h"
#include "Numeric.h"
#include "LibData.h"
#include "MATIDTranslator.h"

using namespace std;

class Library{
  vector<GroupData1D> xs;
  vector<GroupData2D> xs2d;
  vector<int> mat_list;
  vector<int> mt_list;
  vector<int> mat_list2d;
  vector<int> mt_list2d;
 public:
  Library(){};
  int FindData(int mat,int mt);
  int FindData2D(int mat,int mt);
  void ReadCrossSectionFromFiles(string dirname,int nuc,string *filename);
  void ReadCrossSectionFromFile(string dirname,string filename);
  GroupData1D &GetCrossSection(int mat,int mt);
  GroupData1D &GetCrossSection(int i){return xs[i];};
  GroupData2D &GetCrossSection2D(int mat,int mt);
  GroupData2D &GetCrossSection2D(int i){return xs2d[i];};
  int GetMatList(int i){return mat_list[i];};
  int GetMatList2D(int i){return mat_list2d[i];};
  int GetMtList(int i){return mt_list[i];};
  int GetMtList2D(int i){return mt_list2d[i];};
  int GetSize(){int tmp=mat_list.size(); return tmp;};
  int GetSize2D(){int tmp=mat_list2d.size(); return tmp;};

  void GetInformationFromLibData(LibData &ldat,bool thermal_chi=false);
  void GetInformationFromXSLibrary(XSLibrary &xslib,bool thermal_chi=false);
};

class LibraryContainer{
  vector<Library> lib;
  vector<string> libname;
 public:
  LibraryContainer(){};
  void PutLibrary(Library inp,string name);
  void PutXSLibrary(XSLibrary &xslib,string name,bool thermal_chi=false);
  int FindData(string name);
  Library &GetLibrary(string name);
};

class CrossSectionCovariance:public Covariance{
  int mat1,mt1,mat2,mt2;
 public:
  CrossSectionCovariance();
  CrossSectionCovariance(string type);
  ~CrossSectionCovariance(){};
  void PutMatMtInformation(int imat,int imt);
  void PutMatMtInformation(int imat1,int imt1,int imat2,int imt2);  
  void ShowMatMt();
  int GetMat1(){return mat1;};
  int GetMat2(){return mat2;};
  int GetMt1(){return mt1;};
  int GetMt2(){return mt2;};
  void WriteFile(string dirname,string filename);
};

class LibraryCovariance
{
  vector<CrossSectionCovariance> xscov;
  MATIDTranslator midt;
 public:
  LibraryCovariance(){};
  int FindData(int mat1,int mt1,int mat2,int mt2);
  void ReadCovarianceFromFile(string dirname, string filename);
  void ReadCovarianceFromFileOverWriting(string dirname, string filename, real emax);
  void ReadCovarianceFromFile(string dirname, string filename,MATIDTranslator &midt);
  CrossSectionCovariance &GetCrossSectionCovariance(int mat1,int mt1,int mat2,int mt2);
  CrossSectionCovariance &GetCrossSectionCovariance(int mat1,int mt1)  
    {return GetCrossSectionCovariance(mat1,mt1,mat1,mt1);};
  CrossSectionCovariance &GetCrossSectionCovariance(int i){return xscov[i];};
  void CalCovarianceForXSDifference(int mat,int mt1,int mt2);
  void ChangeMTNumber(int mat,int mt1,int mt2);
  void AddCovarianceData(CrossSectionCovariance inp);
  void AddConstantRelativeCovariance(int mat, int mt, int grp, real rel_stdev, bool full_correlation=true);
  void ResetDiagonalElement(int mat, int mt, int grp, real rel_stdev);
  void ReplaceNutCovByNupCov(int mat);
  int GetSize();
  void ShowSelf();
  void SetNoCorrelation();
  void CovarianceForSummedQuantity(int matid, int mtnum, int *mtid, int newmt);
  GroupData2D GetCovarianceMatrixMultipleParameters(int num, int *matid, int *mtid,bool corr=false);
  // ... Relative covariance data is returned. if [corr] is true, correlation matrix is returned.
  void CheckCorrelationMatrix(int mat1, int mt1, int mat2, int mt2);
  void ModifyingCorrelation(int mat1, int mt1);  
};

#endif
