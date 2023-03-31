#ifndef GROUPDATA
#define GROUPDATA

#include <cstdlib>
#include <cmath> // add by hazama
#include <iostream>
#include "Numeric.h"
#include <vector>
#include <string>

using namespace std;

class GroupData{
 protected:
  bool reduced_form;
  int imax;
  vector<real> dat;
 public:
  GroupData(){reduced_form=false;};
  GroupData(int i){put_imax(i); reduced_form=false;};
  ~GroupData();

  void put_data(int i,real f){};
  real get_sum();
  real get_dat(int i)const{return dat[i];};
  void show_self(){};
  void put_imax(int i){imax=i; dat.resize(imax,0.);};
  void set_zero(){for(int i=0;i<imax;i++){dat[i]=0;};};
  int get_imax(){return imax;};
  real GetEuclideanNorm(); // Frobenius norm for matrix case
  vector<real> &get_dat(){return dat;};
  void Factorize(real c){for(int i=0;i<imax;i++){dat[i]*=c;};};
  //void DataClear(){dat.clear();};
  void DataClear();
  bool IsReducedForm() const{return reduced_form;};
};


class GroupData1D:public GroupData{
 public:
  GroupData1D():GroupData(){};
  GroupData1D(int i):GroupData(i){};
  void put_data(real *inp);
  void put_data(vector<real> &inp);
  void put_data(int i,real f){dat[i]=f;};
  void put_data(int i,GroupData1D inp,int j=0,int k=0);
  void add_data(int i,real f){dat[i]+=f;};
  void add_data(GroupData1D &sec);
  real get_dat(int i)const{return dat[i];};
  vector<real> &get_dat(){return dat;};
  GroupData1D get_dat(int i,int j);
  GroupData1D operator+(const GroupData1D &second);
  GroupData1D operator-(const GroupData1D &second);
  GroupData1D operator*(real a);
  GroupData1D mult(GroupData1D sec);
  void Multiply(GroupData1D &sec);
  void Multiply(real a);
  real operator*(const GroupData1D &second);
  GroupData1D operator/(real a);
  GroupData1D operator/(const GroupData1D &second);
  GroupData1D sqrt1d();
  GroupData1D GetInverse();
  void show_self();
  void copy(GroupData1D s);
  void CheckSize(GroupData1D s);
  real Cond(GroupData1D &w);
  GroupData1D Cond(GroupData1D &w,int grp,vector<int> bgrp);
  GroupData1D Cond(GroupData1D &w,int grp,int *bgrp);
  GroupData1D CondSum(int grp,vector<int> bgrp);
  GroupData1D CondSum(int grp,int *bgrp);
};

class GroupData2D:public GroupData{
 private:
  int x,y;
  vector< vector<int> > xpos; // [y][pos]
  vector< vector<real> > val;
 public:
  GroupData2D():GroupData(){y=0; x=0;};
  GroupData2D(int i,int j):GroupData(i*j){y=i; x=j;};
  real get_sum();
  real GetEuclideanNorm();
  void put_data(real *inp);
  void put_data(int i,int j,real inp){dat[i*x+j]=inp;};
  void put_data(int i,real inp){dat[i]=inp;};
  void put_data(int x1,int x2,int y1,int y2,real inp);
  void add_data(int i,int j,real f);
  void put_yx(int i,int j);
  real get_dat(int i);
  real get_dat(int i,int j)const;
  vector<real> &get_dat(){return dat;};
  real InfiniteNorm();
  void show_self(bool fine="true");
  void show_plot(real shift=0.);
  void show_plot(string dir,string filename);
  void show_plot(GroupData1D &val);
  int get_x(){return x;};
  int get_y(){return y;};
  GroupData1D get_diag();
  GroupData1D get_sumx();
  GroupData2D Slicing(int y0,int y1,int x0,int x1);
  GroupData1D Slicing(int y0,int y1,int x0);
  GroupData2D operator*(const GroupData2D &second);
  GroupData2D operator*(real a);
  GroupData2D operator+(const GroupData2D &second);
  GroupData2D operator-(const GroupData2D &second);
  //GroupData1D operator*(GroupData1D g1d);
  GroupData1D operator*(const GroupData1D &g1d);
  //GroupData1D operator*(GroupData1D &g1d);
  //GroupData2D operator/(GroupData1D g1d);
  GroupData2D inverse();
  GroupData2D inverseKatagiri();
  void DoInverse(bool pivot_move=false);
  vector<int> PivotCheck(real eps=1e-10);
  //GroupData2D SquaringForSparceMatrix();
  GroupData2D GetTransposedMatrix();
  GroupData2D T(){return GetTransposedMatrix();};
  void Transposition();

  GroupData2D solveaxb(GroupData2D &sec);
  void solveaxb_mod(GroupData2D &sec);
  GroupData1D solveaxb(GroupData1D &sec);
  void solveaxb_mod(GroupData1D &sec);
  void LUdecomposition();
  void solveaxb_LUdecomposition(GroupData1D &sec);

  void Diagonalization(GroupData2D &diagmat,GroupData2D &evecmat); // M=UDU^{T}
  void Diagonalization(GroupData2D &evecmat); // M=AA^T
  void CalEigenvaluesByLeverrieFaddeev(GroupData1D &eigen);
  void CalEigenvaluesByLF(GroupData1D &eigen){CalEigenvaluesByLeverrieFaddeev(eigen);};
  void Multiply(GroupData1D &sec);
  void PasteGroupData2D(int yp, int xp, GroupData2D &inp);
  void PasteGroupData1D(int yp, int xp, GroupData1D &inp);
  void put_unit();
  void PutUnitMatrix(){put_unit();};
  void copy(GroupData2D sec);
  void CheckSizeX(GroupData1D sec);
  void CheckSizeY(GroupData1D sec);
  void CheckSizeXY(GroupData2D sec);
  void CheckSquare();
  GroupData2D Cond(GroupData1D &w,int grp,vector<int> bgrp);
  GroupData2D Cond(GroupData1D &w,int grp,int *bgrp);
  // (Matrix exponential calculation)
  GroupData2D CalMatrixExponentialByPade(real factor);
  GroupData2D CalMatrixExponentialByPadeKatagiri(real factor);
  GroupData2D CalMatrixExponentialByTaylor(real factor);
  GroupData1D CalMatrixExponentialByKrylov(GroupData1D &inp,real factor,int dim);
  int CalOrderForMatrixExponential(real factor);
  //GroupData2D CalMatrixExponentialByChebyshev(real factor,int order);
  GroupData2D CalMatrixExponentialByChebyshev14(real factor);
  GroupData1D CalMatrixExponentialByChebyshev14(GroupData1D &inp,real factor);
  GroupData2D CalMatrixExponentialByChebyshev16(real factor);
  GroupData1D CalMatrixExponentialByChebyshev16(GroupData1D &inp,real factor);
  // KYMT
  GroupData2D CalMatrixExponentialByChebyshevNew14(real factor);
  GroupData1D CalMatrixExponentialByChebyshevNew14(GroupData1D &inp,real factor);
  GroupData1D CalMatrixExponentialByMMPA18(GroupData1D &inp,real factor);//kawamoto
  GroupData2D CalMatrixExponentialByMMPA32(real factor);//kawamoto
  GroupData1D CalMatrixExponentialByMMPA32(GroupData1D &inp,real factor);//kawamoto
  void LUDecompositionForMatrixExponentialByMMPA32(real factor);
  GroupData1D CalMatrixExponentialByLUDMMPA32(GroupData1D &inp);
  void MultiStepCalc(GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);//kawamoto
  //
  void CalSimultaneousDifferentialEquationByDecomposedMatrixExponential
    (int i1,int i2,real dt,real theta,GroupData2D &mmat, GroupData2D &emat, GroupData2D &xmat,
     GroupData2D &mat0, GroupData2D &mat1, GroupData2D &mat2); 
  // REDUCED FORM
  void ReducedForm();
  void NormalForm();
  int GetXpos(int i,int j) const{return xpos[i][j];};
  real GetVal(int i,int j) const{return val[i][j];};
  void PutVal(int i,int j,real k){val[i][j]=k;};
  void PutVal(int i,real k){val[i].push_back(k);};
  void PutXpos(int i,int j){xpos[i].push_back(j);};
  int GetXposSize(int i) const{int sz=xpos[i].size(); return sz;};

  //Covariance Matrix
  GroupData2D MakeCorrelationMatrix();

  GroupData2D GetPartialCorrelationMatrix();
};


#endif


