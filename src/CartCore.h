#ifndef CARTCORE
#define CARTCORE

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "Numeric.h"

using namespace std;

const int MAX_ASM=600;
const int MAX_MAT=9000;

class Assembly{
 protected:
  vector<int> MaterialID;  
  vector<string> MaterialName;
  vector<real> zdist; // Boundary point between different material along to z-axis
  int    zdiv;     //Number of division along to z-axis
  string AsmName;  //Assembly name
  bool exist;
 public:
  Assembly(){exist=false;};
  Assembly(int zdivi,string *zmapi,real *zdisti,string inp,string type="cumulative")
    {Init(zdivi,zmapi,zdisti,inp,type);};
  Assembly(int zdivi,vector<string> zmapi,vector<real> zdisti,string inp,string type="cumulative")
    {Init(zdivi,zmapi,zdisti,inp,type);};
  //Assembly &operator=(Assembly &sec);
  ~Assembly(){};
  void Init(int zdivi,string *zmapi,real *zdisti,string inp,string type="cumulative");
  void Init(int zdivi,vector<string> zmapi,vector<real> zdisti,string inp,string type="cumulative");
  void show_self();
  int GetZdiv(){return zdiv;};
  void PutMaterialID(int i,int j){MaterialID[i]=j;};
  int GetMaterialID(int i){return MaterialID[i];};
  real GetZdist(int i){return zdist[i];};
  void PutMaterialName(int i,string inp){MaterialName[i]=inp;};
  string GetMaterialName(int i){return MaterialName[i];};
  string GetName(){return AsmName;};
  bool Exist(){return exist;};
};

class AssemblySet{
 protected:
  int AsmNum;                  //Number of included 'Assembly' instances
  int MatNum;                  //Number of included material
  Assembly Asm[MAX_ASM];       //Included 'Assembly' instances
  string MaterialName[MAX_MAT]; //Name of included material
  bool exist;
 public:
  AssemblySet(int matinp);
  ~AssemblySet(){};
  void PutAsm(Assembly &inp);  //Put instance of 'Assembly' into this class
  void PutMaterialName(vector<string> inp);  
  void PutMaterialName(string *inp);  
  void show_self();
  int GetAsmNum(){return AsmNum;};
  int GetMatNum(){return MatNum;};
  string GetMaterialName(int i){return MaterialName[i];};
  Assembly &GetAsm(int i){return Asm[i];};
  Assembly GetJointAssembly(int n, string *r);
};

class XYLattice{
 protected:
  int xm,ym;
  vector<real> xl,yl;
  vector<int> map;
 public:
  XYLattice(){xm=0; ym=0;};
  void PutXY(int x,int y);
  void PutXL(real *inp);
  void PutYL(real *inp);
  void PutMap(int *mapinp);
  void PutMap(int x,int y,int *mapinp);
  void PutData(int x,int y,real *xl,real *yl,int *mapinp);
  void PutSingleData(int mat);
  //
  int GetXm(){return xm;};
  int GetYm(){return ym;};
  real GetXl(int i){return xl[i];};
  real GetYl(int i){return yl[i];};
  real GetXlSum();
  real GetYlSum();
  int GetMap(int x,int y){return map[y*xm+x];};
};

class XYLatticeSet{
 protected:
  int latnum;
  vector<XYLattice> xylat;
 public:
  XYLatticeSet(){latnum=0;};
  void PutXYLattice(XYLattice &latinp);
  int GetLatnum(){return latnum;};
  XYLattice &GetXYLattice(int i){return xylat[i];};
};

class CartCore{
 protected:
  int xr,yr;          //Region number along to x- and y-axis
  vector<int> AsmMap; //Assembly map on xy-plane
  vector<real> xwid;
  vector<real> ywid;
  int LeftBC,RightBC,BackBC,FrontBC,UpperBC,BottomBC; 
        // Boundary condition (0/1:vacuum/reflective)
  AssemblySet *AsmSet; //'AssemblySet' instance
 public:
  /** AssemblyMap(AssemblyID is described in XY plane)   **/
  /** MaterialMap(MaterialID is described in XYZ system) **/
  CartCore(){};
  CartCore(int x,int y){Initialize(x,y);};
  ~CartCore();
  void Initialize(int x,int y);
  void PutAsmMap(vector<int> inp);
  void PutAsmMap(int *inp);
  void PutAsmMap(int x1,int x2,vector<int> inp);
  void PutAsmMap(int x1,int x2,int *inp);
  void PutAsmMap_Hex_to_XY(int *inp);
  void PutAsmSet(AssemblySet *inp){AsmSet=inp;};
  void PutWidthXY(vector<real> xinp,vector<real> yinp);
  void PutWidthXY(real *xinp,real *yinp);
  void PutWidthFromPitch(real pitch);
  void PutBC(vector<int> inp);
  void PutBC(int *inp);
  void PutBoundaryCondition
    (string xl="Vacuum",string xr="Vacuum",
     string yl="Vacuum",string yr="Vacuum",
     string zl="Vacuum",string zr="Vacuum");
  void ShowAssemblyMapNew();
  void ShowAssemblyMap();
  void ShowAssemblyMap(int ii, int ij=-1);
  void ShowMaterialMap();
  void MakeMaterialMap(int &zr,int *MaterialMap,real *zwid);
  void GetMaterialMapXY(real zinp,int *map); /** To get M.Map at z=zinp **/
  int GetUnifiedZMesh(); /** To get number of unified Z mesh **/
  // for XYLattice
  void PutLatticeMap(int x,int y,int *latmap,real *xwid,real *ywid,XYLatticeSet &latset,bool print=false);
  // 
  void ChangeAssembly(int x,int y,int asmid);
  void ChangeAssembly(int p,int* xp,int *yp,int asmid);
  void ChangeAssembly(int asmid1,int asmid2);
  //Output of member function
  int GetLeftBC(){return LeftBC;};
  int GetRightBC(){return RightBC;};
  int GetBackBC(){return BackBC;};
  int GetFrontBC(){return FrontBC;};
  int GetUpperBC(){return UpperBC;};
  int GetBottomBC(){return BottomBC;}; 
  int GetXr(){return xr;};
  int GetYr(){return yr;};
  real GetXwid(int i){return xwid[i];};
  real GetYwid(int i){return ywid[i];};
  int GetAsmMap(int i){return AsmMap[i];};
  int UsedOrNotAssembly(int i); // 0/1:No/Yes
  AssemblySet *GetAset(){return AsmSet;};
  void CountAssembly(int i);
};

#endif
