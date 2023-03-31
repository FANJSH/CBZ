#ifndef REGION
#define REGION

#include <math.h>
#include <vector>
#include "GroupData.h"
#include "ABEMIE_plane.h"
#include "ABEMIE_const.h"
#include "ABEMIE_medium.h"

using namespace std;

const int max_be=40;

class region{
  plane be[max_be];
  BemMedium *med;
  int num, nump;
  vector<GroupData2D> gmat,hmat,dgmat,dhmat;
  vector<GroupData2D> dfdf;
  real volume;
  vector< vector<GroupData1D> > RealFlux;
 public:
  //Constructor
  region();
  void Init(){num=0; nump=0;};
  void add_plane(real x1,real y1,real x2,real y2,int i,int j);
  void cal_gh();
  void copy_gh(region sec);
  plane& get_be(int i){return be[i];};
  BemMedium *get_med(){return med;};
  void set_plane_kind(int i,int j){be[i].set_kind(j);};
  int  get_plane_kind(int i){return be[i].get_kind();};
  void put_BemMedium(BemMedium *i){med=i;};
  int get_num(){return num;};
  int get_nump(){return nump;};
  void cal_matrix(bool sp3=false);
  GroupData1D GetRegionFlx(int grp);
  GroupData1D GetRegionCur(int grp);
  void PutDfdkToPlane(int grp,GroupData1D inp);
  void PutCurToPlane(int grp, GroupData1D cur,bool sp3);
  GroupData2D& get_dfdf(int i){return dfdf[i];};
  void  put_vol(real i){volume=i;};
  real get_vol(){return volume;};
  void put_vert();
  int get_address(int ip);
  real cal_real_flx();
  void change_pnt(int i);
  void set_array();
  real get_integrated_flux(int g);
};

class MakeRegion{
  int imax,pnt;
 public:
  MakeRegion(int i,int ip){imax=i; pnt=ip;};
  region square(real x,real y,real xx,real yy,int i,int j,int k=0);
  region tri(real x,real y,real d,int i,int j=0);
  region drawing(real xs,real ys,int pv,vector <vector<real> > vec);
};

#endif
