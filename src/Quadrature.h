#ifndef QUADRATURE
#define QUADRATURE

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "SphericalHarmonics.h"
#include "Numeric.h"
#include "GroupData.h"

using namespace std;

class Quadrature{
 protected:
  int Ndim, sn, pl, plnum;
  vector<real> omega, mu, eata, xi;
  vector< vector<real> > val;
  vector<int> xref, yref, zref, xyzref;
  GroupData2D gmatrix, gmatrix_inv; // Galerkin matrix
  vector<int> l_array_gmatrix;// Legendre-moment array
  bool exist;
 public:
  Quadrature(){};
  Quadrature(int Ndim, int pl=0){Initialize(Ndim,pl);};
  ~Quadrature();
  void Initialize(int Ndim,int pl);
  bool Exist(){return exist;};
  void PutSN(int sninp);
  void CalXYZref();
  void CalDiscretedSphericalHarmonics();
  void WeightNormalize(real sum=2.); // Normalized to 2.
  real GetMoment(int i,int is){return val[is][i];};
  real GetMoment(int l,int m,int is);
  void PutData(real *mui,real *eai,real *xii,real *w);
  void ArrangeDiscretedPoint(real pp,real *mui,real *eai,real *xii);
  void Rotation(real theta,string ax);
  void CheckOrthogonalityTo00(real eps=1e-4);
  void show_self();
  void show_self_konno();
  // Access function
  int GetSN(){return sn;};
  real GetMu(int i){return mu[i];};
  real GetEata(int i){return eata[i];};
  real GetXi(int i){return xi[i];};
  int GetXref(int i){return xref[i];};
  int GetYref(int i){return yref[i];};
  int GetZref(int i){return zref[i];};
  int GetXYZref(int i){return xyzref[i];};
  real GetOM(int i){return omega[i]*mu[i];};
  real GetOmega(int i){return omega[i];};
  int GetPlnum(){return plnum;};

  // ++ Quadrature set ++
  // Level-symmetric quadrature set
  void PutLevelSymmetric(int sn);
  void LevelSymmetric(int sn,real *p,real *w,int *wm);
  void PutLSS2();
  void PutLSS4();
  void PutLSS6();
  void PutLSS8();
  void PutLSS10();
  void PutLSS12();
  void PutLSS16();
  void PutLSS20();
  // Triangular-type Double-Gaussian Tchebyshev
  void PutTriangularDPnTn(int dpn, string dir="z");
  void PutTriangularPnTn(int pn, string dir="z");
  // Rectangular-type Double-Gaussian Tchebyshev
  void PutRectangularDPnTn(int dpn, int tn, string dir="z");
  void PutRectangularPnTn(int pn, int tn, string dir="z");
  // Rectangular-type Double-Gaussian Double-Gaussian (x1:dir(cos), x2:others(phi))
  void PutRectangularDPnNxN(int x1,int x2=-1,string dir="z");
  // Triangular-type TChebyshev-TChebyshev
  void PutTriangularTn(int sn);
  // EOn
  void PutEO4();
  void PutEO6();
  void PutEO8();
  void PutEO10();
  void PutEO12();
  void PutEO14();
  void PutEO16();
  void PutAEO8();
  // Lebedev
  void PutLebedev();

  // Galerkin matrix
  GroupData2D GetGalerkinMatrix(){return gmatrix;};
  GroupData2D GetGalerkinMatrixInverse(){return gmatrix_inv;};
  int GetLArrayGMatrix(int i){return l_array_gmatrix[i];};

  // The following functions should be removed in future.
  void CalValue(){CalDiscretedSphericalHarmonics();};
  // DPnDPnTn
  //void PutDPnDPnTn(int dpn1,int dpn2,int tn,real bmu,string dir="z");
  // * dpn1 is applied into bmu~1, dpn2 is applied into 0~bmu *  
};


#endif
