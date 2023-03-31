#ifndef CORE
#define CORE

#include <time.h>
#include <string.h>
#include <math.h>
#include "GroupData.h"
#include "CartCore.h"
#include "IrregularGeometryInformation.h"
#include "ABEMIE_const.h"
#include "ABEMIE_region.h"
#include "ABEMIE_medium.h"
#include <vector>
#include<iostream>

using namespace std;

enum BCondition{Zeroflux,Reflective,Vacuum};

class core{
  region reg[maxreg];
  BemMedium med[maxmed];
  int numreg,nummed,numib,numibp,imax,count;
  int ibr[maxib],ibp[maxib],ibrr[maxib],ibrp[maxib];
  int rmed[maxreg];
  int eql[maxreg];
  real lamda;
  GroupData1D f, df;
  GroupData1D res;
  vector<int> jinfo;
  real keff,keffold;
  vector<GroupData2D> jother;
  vector<GroupData2D> jself;
  vector<GroupData1D> jk;
  vector<GroupData1D> dp2;
  vector<int> ledge;
  vector<int> redge;
  bool sp3;
 public:
  core(){Init();};
  core(int i){Init(); put_imax(i);};
  ~core();
  void Init();
  void put_imax(int i){imax=i;};
  void AddRegion(region i);
  void AddMedium(Medium i);
  void put_rmed(int *i);
  void put_eql(int *i);
  void put_boundary_info();
  void show_boundary_info();
  void put_outerb(real x1,real y1,real x2,real y2,BCondition b,GroupData1D inp);
  void show_innerb_info();
  void modf_flx(bool accel);
  void put_flx_innerr();
  int get_numib(){return numib;};
  int get_numreg(){return numreg;};
  void set_array_uppercal();
  void cal_f();
  void cal_df_p(real errf,real errk,int iter);
  void cal_df_p0();

  void cal_df_direct_inversion();
  void set_edgedata_for_direct_inversion();

  void get_power();
  real get_integrated_flux(int r,int g);
  void put_jinfo();
  plane& get_ibe(int i){return reg[ibr[i]].get_be(ibp[i]);};
  plane& get_iber(int i){return reg[ibrr[i]].get_be(ibrp[i]);};
  GroupData2D get_cc(int i,int j);
  int get_pnt(int i){return reg[ibr[i]].get_be(ibp[i]).get_pnt();};
  GroupData1D cal_jfk(int ib,int is1,int is2);
  real cal_real_flx();
  void out_flx(int ig);
  void change_pnt(int i);
  void put_keff(real k){keff=k;};
  real get_keff(){return keff;};
  void initial_flx();
  void run(real epsk=1e-5,real epsf=1e-4,bool accel=true,bool direct=false);
  void PutMeshXYFromCartCore(CartCore &cinp,real z,int pnt);
  void PutIrregularGeometryInformation(IrregularGeometryInformation &inp,int *reg_med,int pnt);

  void SP3_Calculation(){sp3=true;};
  void put_volume(int ir,real vinp){reg[ir].put_vol(vinp);};
  real get_volume(int ir){return reg[ir].get_vol();};
};
#endif
