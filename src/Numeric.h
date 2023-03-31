#ifndef PIDEF
#define PIDEF

#include <string>
#include <vector>
#include <math.h>

typedef double real;
typedef float real4;
typedef double real8;

const real PI =3.1415926535897932384626433;
const real PI2=PI*2.;
const real PI4=PI*4.;
const real INV_PI4=1./PI4;
const real INV_PI=1./PI;

const real mev_to_j=1.60219e-19*1e6;

const real avo=0.60221367;
const real ev_to_j=1.60219e-19;

// For random number generator
//const double fact=1./RAND_MAX;
const double fact=1./2147483647;

enum xstype
  {sigt,nusigf,sigtr,sigr,siga,d,dr,dz,sign2n,sigs,sigf,sigc,nu,mu,sigel,chi,siginel,sigst,kerma,
   coh_el,incoh_el,pair_pro};


const std::string xstype_string[]={
  "total    ","product. ","transport","removal  ","absorp.  ",
  "diffusion","diff_r   ","diff_z   ","(n,2n)   ","scatter. ",
  "fission  ","capture  ","nu       ","mu       ","elastic  ",
  "chi      ","inelastic","tot-scat ","kerma    ","coh-elast",
  "incoh-ela","pair-prod"
};

using namespace std;

void ChangeValue(real &a,real &b);
void ChangeValue(int &a,int &b);
real kaijo(int i,int j);

void gauss(real *a,real *b,int n,int m);
void gauss(vector<real> &a,vector<real> &b,int n,int m);
void lu_decomposition(vector<real> &a,int n);
void lu_decomposition_inverse(vector<real> &a,int n);
void gauss_inverse(vector<real> &a, int n);
void gauss_pivot_check(vector<real> &a, vector<int> &pos, int n, real eps=1e-10);
void gauss_3p(real *a,real *b,int n,int m); // for diagonal matrix(width=3)
void gauss22(real *a,real *b); // for
void gauss22(vector<real> &a,vector<real> &b); // for
void ChangeOrder(vector<real>& a,int ip,vector<int>& ret);
void ChangeOrder(real *a,int ip,int *ret);
void ChangeOrder(real *a,int ip);
real GetEnx(int n,real x);  
int TranslateNuclideIDFromJFS(int i);
real GetNuclideWeight(int id);

void MakeFineMap(int x,int y,int *inp,int xst,int xed,int c1,int *cell1,int *cell2,int *cell3,int *cell4,int *mapnew);
void MakeXMirrorMap(int x,int y,int *inp,int *mapnew);
void MakeYMirrorMap(int x,int y,int *inp,int *mapnew);
void WriteXMirrorMap(int x,int y,int *inp);
void WriteYMirrorMap(int x,int y,int *inp);

real CalAreaMaximizeYPoint(real x1,real y1);
void CircleToSquare(real r,int n);
void RedividingEquivolumeRing(real *rinp, int pos);

bool CalCubicSplineCoefficient
(real *x, real *y, real *dy,int n,real *c,real *d,real *e);
real CubicSplineFittingNew(real *x,real *y,int n,real xinp);
real CubicSplineFitting(real *x,real *y,int n,real xinp);
real LinearLinearInterpolation(real x1,real x2,real y1,real y2,real x);

bool SameReal(real a,real b,real f);
string IntToString(int i);
int StringToInt(string i);
real StringToReal(string i);

// Conjugate Gradient method

void cg(real *a,real *b,int n);
void bicg(real *a,real *b,int n);

void MatrixDiagonalization(real *a, real *u, int n);

// Vector manipulation
void MovingAveraging(vector<real> &data, int points);
void Normalize(vector<real> &data, real factor=1.);

// Random number generator
double ransu(); // [-1,1]
double ransu01(); // [0,1]
real ransu_gauss(real mu, real sigma);
real GetStatics(vector<real> data, int num, bool print=false, bool positivity=false, int print_div=100);
void GetStatics(vector<real> &data, int num, vector<real>& result);
void GetStatics(vector<real> &data, vector<real>& result);
void GetStatics(vector<real> &data, int st, int ed, vector<real>& result);
void GetMeanAndVariance(vector<real> data, int num, real &mean, real &var);
real GetStatics(vector<real> data1, vector<real> data2, int num);
void GetMeanAndVarianceForLognormal(real mu_org, real var_org, real &mu_log, real &var_log);

// Polynomial handler
real GetValue(vector<real> &z, real x);
real SolveZeroPoint(vector<real> &z, real vl,real vh);
void PolynomialDivision(vector<real> &z,real &p,real &q,vector<real> &b,vector<real> &c);
void BairstowHitchcock(vector<real> &z,real &p,real &q);

// Noise reduction
void NoiseReductionBySVD_org(string filename, string outname, real dx, int svd_order, int pnt=1);
void NoiseReductionBySVD(string filename, string outname, real dx, int svd_order, int pnt=1);
vector< vector<real> > NoiseReductionBySVD(vector<real> &data_in, real dx, int svd_order, int pnt=1);
vector< vector<real> > NoiseReductionBySVD(vector<real> &data_in, int st, int ed, real dx, int svd_order, int pnt=1);

// Output format
void WriteOut(int a,int col);
void WriteOut(string a,int col);
void WriteOut(real a,string b);

// Cross section type
enum xstype XSType(int i);

// Prepared by Yanagihara-kun
void LeastSquaresMethod(int data_num, double *x, double *y, double &a, double &b);
void LeastSquaresMethodQuadr(int data_num, double *x, double *y ,double &a, double &b, double &c);
void LeastSquaresMethodCubic(int data_num, double *x, double *y ,double &a, double &b, double &c, double &d);

#endif

