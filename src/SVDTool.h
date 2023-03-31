#ifndef SVDTOOL
#define SVDTOOL

#include "GroupData.h"
#include "Numeric.h"

class SVDTool{
 protected:
  int x,y;
  GroupData2D mtx; // row:y, column:x
  vector<real> eigen_m; // eigenvalues for A A^T (<y)
  vector<real> eigen_n; // eigenvalues for A^T A (<x)
  vector<int> order_m;
  vector<int> order_n;
  GroupData2D umat; // size:y*y, eigenvectors for A A^T 
  GroupData2D vmat; // size:x*x, eigenvectors for A^T A
 public:
  SVDTool(){};
  void SVD(GroupData2D &mtxin);
  real GetUdata(int order,int pos);
  real GetVdata(int order,int pos);
  GroupData1D GetUVector(int order);
  GroupData1D GetVVector(int order);
  real GetEigenM(int order);
  real GetEigenN(int order);
  void ShowEigenM();
  void ShowEigenN();
  vector<real> GetSingularValue(real eps=1e-10, bool print=true);
  void MatrixReconstructionCheck(int order);
  GroupData2D GetUmat(real eps=1e-10);
  GroupData2D GetVmat(real eps=1e-10);
  GroupData2D CalPseudoInverse(real eps=1e-20); // It is strongly recommended to read the source code where an important note exists.
};

#endif
