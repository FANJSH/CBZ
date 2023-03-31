#ifndef CARTMESHINFO
#define CARTMESHINFO

#include <time.h>
#include <string>
#include <math.h>
#include "CartCore.h"
#include <vector>

using namespace std;

class CartMeshInfo{
  int FMesh[3],CMesh[3];
  int BC[6];
  vector< vector<real> > FMeshL;
  vector< vector<int> >  CMeshF;
  vector< vector<real> > CMeshL;
  vector<int> FMat; // including "-1" medium
  vector<int> FMat_par_mesh; // NOT including "-1" medium
  vector<int> CMat;
 public:
  CartMeshInfo(){};
  void PutMeshXY(int *xm,int *ym,real z,CartCore &core);
  void PutMeshXYZ(int *xm,int *ym,real zwidf,CartCore &core,real zwidc=20.);
  void PutMeshXYZ(int *xm,int *ym,int *zm,CartCore &core,real zwidc=20.);
  void PutMeshXY(vector<int> xm,vector<int> ym,real z,CartCore &core);
  void PutMeshXYZ(vector<int> xm,vector<int> ym,real zwidf,CartCore &core,real zwidc=20.);
  void PutMeshInfo(int xr,int yr,int zr,int *fmx,int *fmy,int *fmz,
               real *xl,real *yl,real *zl,int *mat,string type="width");
  void PutMeshInfo(int xr,int yr,int *fmx,int *fmy,real *xl,real *yl,int *mat,string type="width");
  void PutMeshInfo(int xr,int *fmx,real *xl,int *mat,string type="width");
  void PutMeshInfo(int xr,int yr,int zr,vector<int> fmx,vector<int> fmy,vector<int> fmz,
               vector<real> xl,vector<real> yl,vector<real> zl,vector<int> mat,string type="width");
  void PutMeshInfo(int xr,int yr,vector<int> fmx,vector<int> fmy,vector<real> xl,vector<real> yl,vector<int> mat,string type="width");
  void PutMeshInfo(int xr,vector<int>fmx,vector<real>xl,vector<int>mat,string type="width");
  int GetFMesh(int i){return FMesh[i];};
  int GetXF(){return FMesh[0];};
  int GetYF(){return FMesh[1];};
  int GetZF(){return FMesh[2];};
  int GetXC(){return CMesh[0];};
  int GetYC(){return CMesh[1];};
  int GetZC(){return CMesh[2];};
  int GetCMesh(int i){return CMesh[i];};
  real GetFMeshL(int i,int j){return FMeshL[i][j];};
  int GetCMeshF(int i,int j){return CMeshF[i][j];};
  real GetCMeshL(int i,int j){return CMeshL[i][j];};
  int GetFMat(int i){return FMat[i];};
  int GetFMatParMesh(int i){return FMat_par_mesh[i];};
  int GetCMat(int i){return CMat[i];};
  void PutBoundaryCondition(int *binp);
  void PutBoundaryCondition(vector<int> binp);
  void PutBoundaryCondition
    (string xl="Vacuum",string xr="Vacuum",
     string yl="Vacuum",string yr="Vacuum",
     string zl="Vacuum",string zr="Vacuum");
  int GetBC(int i){return BC[i];};

  void ReconstructCoarseMesh(real maxxl,real maxyl,real maxzl);
  bool CheckCoarseMeshMixedVacuum();

  int GetMeshPosition(int dir,real x);
  int GetXMeshPosition(real x){return GetMeshPosition(0,x);};
  int GetYMeshPosition(real x){return GetMeshPosition(1,x);};
  int GetZMeshPosition(real x){return GetMeshPosition(2,x);};
  real GetMeshLocation(int dir,int x);
  real GetXMeshLocation(int x){return GetMeshLocation(0,x);};
  real GetYMeshLocation(int x){return GetMeshLocation(1,x);};
  real GetZMeshLocation(int x){return GetMeshLocation(2,x);};

  void show_self();
  void GetPositionFromMeshID(int id);
};

#endif
