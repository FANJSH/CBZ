#ifndef PLANE
#define PLANE

#include <string.h>
#include <math.h>
#include <vector>

#include "GroupData.h"
#include "GeomVector.h"
#include "GeomPlane.h"
#include "Numeric.h"
#include "Bessel.h"
#include "ABEMIE_const.h"

using namespace std;

class plane:public GeomPlane{
  int imax,kind,pnt;
  vector<GroupData1D> flx;
  vector<GroupData1D> cur;
  vector<GroupData1D> dfdk;
 public:
  //Constructor
  plane():GeomPlane(){};
  plane(GeomVector i1,GeomVector i2):GeomPlane(i1,i2){};
  plane(GeomVector i1,GeomVector i2,int i):GeomPlane(i1,i2){set_imax(i);};
  plane(real x1,real y1,real x2,real y2,int i,int j):GeomPlane(x1,y1,x2,y2){put_pnt(j); set_imax(i);};

  void put_data(real x1,real y1,real x2,real y2,int i,int j){PutVector(x1,y1,x2,y2); put_pnt(j); set_imax(i);};

  void cal_ghij(real b,real bdk,GeomVector &pp,real *gr,real *hr,real *dgr,real *dhr);
  void cal_ghii(real b,real bdk,GeomVector &pp,real *gr,real *hr,real *dgr,real *dhr);
  real CalcDistance2(GeomVector p);
  int opposite_p(int i);
  GeomVector get_pos(int i);
  void change_pnt(int i);

  //Access function for input
  void set_imax(int i);
  void set_kind(int i){kind=i;};
  void put_data1da(int i,char *ss,GroupData1D inp);
  void put_pnt(int i);

  //Access funtion for output
  int get_pnt(){return pnt;};
  int get_imax(){return imax;};
  int get_kind(){return kind;};

  GroupData1D &GetFlx(int i){return flx[i];};
  GroupData1D &GetCur(int i){return cur[i];};
  GroupData1D &GetDfdk(int i){return dfdk[i];};
};

#endif
