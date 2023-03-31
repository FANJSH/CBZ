#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <vector>
#include "SVDTool.h"

using namespace std;

void SVDTool::SVD(GroupData2D &mtxin)
{
  mtx.copy(mtxin);

  GroupData2D mtxt=mtx.GetTransposedMatrix();

  GroupData2D mmat=mtx*mtxt; // A A^T
  GroupData2D nmat=mtxt*mtx; // A^T A

  x=mtx.get_x();
  y=mtx.get_y();

  // +++ Matrix diagonalization
  real *a_m=new real[y*y];
  real *b_m=new real[y*y];
  real *a_n=new real[x*x];
  real *b_n=new real[x*x];

  int id=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<y;j++){
      a_m[id++]=mmat.get_dat(i,j);
    };
  };

  id=0;
  for(int i=0;i<x;i++){
    for(int j=0;j<x;j++){
      a_n[id++]=nmat.get_dat(i,j);
    };
  };

  MatrixDiagonalization(a_m,b_m,y);
  MatrixDiagonalization(a_n,b_n,x);
  // C=U^{-1}PU
  // (input)  a:P 
  // (output) a:C, b:U

  eigen_m.resize(y);
  eigen_n.resize(x);
  for(int i=0;i<y;i++){
    eigen_m[i]=a_m[i*y+i];
  };
  for(int i=0;i<x;i++){
    eigen_n[i]=a_n[i*x+i];
  };

  order_m.resize(y);
  order_n.resize(x);
  ChangeOrder(eigen_m,y,order_m);
  ChangeOrder(eigen_n,x,order_n);

#if 0  
  cout<<"# A set of eigenvalues of A A^T\n";
  for(int i=0;i<y;i++){
    cout<<i<<" "<<eigen_m[i]<<" "<<order_m[y-1-i]<<"\n";
  };
  cout<<"\n";
  cout<<"# A set of eigenvalues of A^T A\n";  
  for(int i=0;i<x;i++){
    cout<<i<<" "<<eigen_n[i]<<" "<<order_n[x-1-i]<<"\n";
  };
  cout<<"\n";
#endif  

  // --- Sign checking ---
  int xymin=min(x,y);
  for(int i=0;i<xymin;i++){
    real tmp=0.;
    for(int j=0;j<x;j++){
      tmp+=mtx.get_dat(0,j)*b_n[j*x+order_n[x-1-i]];
    };
    if(tmp*b_m[order_m[y-1-i]]<0.){
      for(int j=0;j<x;j++){
	b_n[j*x+order_n[x-1-i]]*=-1.;
      };
    };
  };

  umat.put_yx(y,y);
  vmat.put_yx(x,x);
  umat.put_data(b_m);
  vmat.put_data(b_n);

  delete [] a_m;
  delete [] b_m;
  delete [] a_n;
  delete [] b_n;
};

GroupData1D SVDTool::GetUVector(int order)
{
  GroupData1D ret(y);
  for(int i=0;i<y;i++){
    ret.put_data(i,umat.get_dat(i,order_m[y-1-order]));
  };
  return ret;
};

real SVDTool::GetUdata(int order,int pos)
{
  return umat.get_dat(pos,order_m[y-1-order]);
};

GroupData1D SVDTool::GetVVector(int order)
{
  GroupData1D ret(x);
  for(int i=0;i<x;i++){
    ret.put_data(i,vmat.get_dat(i,order_n[x-1-order]));
  };
  return ret;
};

real SVDTool::GetVdata(int order,int pos)
{
  return vmat.get_dat(pos,order_n[x-1-order]);
};

real SVDTool::GetEigenM(int order)
{
  return eigen_m[order_m[y-1-order]];
};

real SVDTool::GetEigenN(int order)
{
  return eigen_n[order_n[x-1-order]];
};

void SVDTool::ShowEigenM()
{
  cout<<"#\n# Eigenvalue of M (AA^T) matrix\n#\n";
  for(int i=0;i<y;i++){
    cout<<i<<" "<<GetEigenM(i)<<"\n";
  };
};

void SVDTool::ShowEigenN()
{
  cout<<"#\n# Eigenvalue of N (A^TA) matrix\n#\n";
  for(int i=0;i<x;i++){
    cout<<i<<" "<<GetEigenN(i)<<"\n";
  };
};

vector<real> SVDTool::GetSingularValue(real eps, bool print)
{
  cout.setf(ios::scientific);
  cout.precision(10);
  if(print){
    cout<<"#\n";
    cout<<"# Singular value from M (AA^T) and N (A^T A)\n";
    cout<<"#   (Imaginary values are regarded as zero.)\n";
    cout<<"#   (If two values are different from each other,\n";
    cout<<"#    you need to be careful about the results.)\n#\n";
  };
  int minxy=min(x,y);
  vector<real> sv;
  for(int i=0;i<minxy;i++){
    real valm=GetEigenM(i);
    real valn=GetEigenN(i);
    if(valm<0.)valm=0.;
    if(valn<0.)valn=0.;
    if(valm>eps&&valn>eps){
      if(print)cout<<i<<"    "<<sqrt(valm)<<" "<<sqrt(valn)<<"\n";
      real diff=(valm-valn)/valn;
      if(diff<1e-5){
        sv.push_back(sqrt(valm));
      }else{
        cout<<"# Error in SVDTool::GetSingularValue.\n";
        cout<<"# Two eigenvalues for A^T A and A A^T are inconsistent.\n";
	cout<<"# "<<sqrt(valm)<<" "<<sqrt(valn)<<"\n";
        exit(0);
      };
    };
  };
  return sv;
};

void SVDTool::MatrixReconstructionCheck(int order)
{
  cout<<"#\n# Matrix reconstruction check\n#\n";
  cout<<"# The number of basis vectors (dimension) : "<<order<<"\n#\n";
  for(int g=0;g<x;g++){
    for(int g2=0;g2<y;g2++){
      real tmp=0.;
      for(int l=0;l<order;l++){
	tmp+=GetUdata(l,g2)*sqrt(GetEigenN(l))*GetVdata(l,g);
      };
      real adiff=tmp-mtx.get_dat(g2,g);      
      real ratio=tmp/mtx.get_dat(g2,g);
      cout<<"# y="<<g2<<", x="<<g<<" : ";
      cout<<"Org. data = "<<mtx.get_dat(g2,g)<<" / ";
      cout<<"Abs. Diff. to Ref. = "<<adiff<<" / ";
      if(mtx.get_dat(g2,g)!=0.)cout<<"Ratio to Ref. = "<<ratio;
      cout<<"\n";
    };
  };
};

GroupData2D SVDTool::GetUmat(real eps)
{
  vector<real> dummy=GetSingularValue(eps, false);
  int sz=dummy.size();
  GroupData2D ret(y,sz);
  for(int i=0;i<sz;i++){
    GroupData1D tmp=GetUVector(i);
    ret.PasteGroupData1D(0,i,tmp);
  };
  return ret;
};

GroupData2D SVDTool::GetVmat(real eps)
{
  vector<real> dummy=GetSingularValue(eps, false);
  int sz=dummy.size();
  GroupData2D ret(x,sz);
  for(int i=0;i<sz;i++){
    GroupData1D tmp=GetVVector(i);
    ret.PasteGroupData1D(0,i,tmp);
  };
  return ret;
};

GroupData2D SVDTool::CalPseudoInverse(real eps)
{
  // This method should be used carefully.
  //
  // To calculate the pseudo inverse, we have to use the inverse of the singular values
  // as [dmat.inverse].
  // If the singular values which are significantly small,
  // these contributions should become large in the final pseudo inverse matrix.
  // It means that the choice of the "meaningful" singular values is quite important.
  // Relevant explanations are given in the master thesis of Watanabe-kun of Nagoya Univ.
  // which can be retreived from the reactor physics division website of AESJ.
  //
  // The default value for [eps] was set at 1e-10,
  // but this was changed to 1e-20 in 2022/5/22.
  
  vector<real> dummy=GetSingularValue(eps, false);
  int sz=dummy.size();

  GroupData2D umat=GetUmat(eps);
  GroupData2D vmat=GetVmat(eps);

  GroupData2D dmat(sz,sz);
  dmat.set_zero();
  for(int i=0;i<sz;i++){
    dmat.put_data(i,i,dummy[i]);
  };

  //dmat.show_self();

  // ... Original matrix calculation for testing
  //GroupData2D ret=umat*dmat*vmat.T();

  GroupData2D ret=vmat*dmat.inverse()*umat.T();
  
  return ret;
};


