#ifndef IRREGULAR_GEOMETRY_INFORMATION
#define IRREGULAR_GEOMETRY_INFORMATION

#include <fstream>
#include <vector>
#include "GeomPlane.h"
#include "GeomVector.h"
#include "Numeric.h"

using namespace std;

// Geometry is expressed by combinations of the following FUNDAMENTAL geometries
//
//  - GeomPolygon
//  - GeomCircle
//  - GeomDividedCircle

enum GeomKind{Circle,Polygon,QCircle};

class GeomGeneral{
 protected:
  GeomVector center;
  int RegionID;
  real outer_circle_r;
public:
  GeomGeneral();
  void PutRegionID(int inp){RegionID=inp;};
  int GetRegionID(){return RegionID;};
  GeomVector GetCenter(){return center;};
  void PutCenter(real x,real y);
  void PutCenter(GeomVector &vinp);
  void PutOuterCircleR(real r){outer_circle_r=r;};
};

class GeomCircle:public GeomGeneral{
 protected:
  real r;
 public:
  GeomCircle():GeomGeneral(){};
  GeomCircle(real xinp,real yinp,real rinp):GeomGeneral(){Init(xinp,yinp,rinp);};
  void Init(real xinp,real yinp,real rinp);
  void PutR(real i);
  real GetR(){return r;};
  GeomPlane CrossPlane(GeomPlane &pl);
  GeomCircle GetReflectedCopy(real xvec,real yvec,int d_regid=0);
  GeomCircle GetRotatedCopy(real rot, int d_regid=0);
  void Shift(real dx, real dy);
};

class GeomPolygon:public GeomGeneral{
 protected:
  int numpl;
  vector<GeomPlane> plane;
  int kind; // 1:rectangular 2:hexagon 3:specific-triangle
 public:
  GeomPolygon():GeomGeneral(){kind=0; numpl=0;};
  void PutRectangular(real x,real y,real ai,real bi);
  void PutHexagon(real x,real y,real ai); // ai is TEIHEN NO NAGASA
  void PutTriangle(real x,real y,real len,int quad);
  void PutTriangleHalf(real x,real y,real len);  
  void PutTriangle(real x,real y,real lenx,real leny,int quad);
  void PutEquilateralTriangle(real x,real y,real len,real angle=0.);
  GeomPlane CrossPlane(GeomPlane &pl);
  GeomVector SymmVec(GeomVector &inp);
  GeomPlane GetPlane(int i){return plane[i];};
  real GetSurface();
  int GetKind(){return kind;};
  int GetNumpl(){return numpl;};
  void PutKind(int i){kind=i;};
  void PutNumpl(int i);
  real GetLen(int i){return plane[i].get_long();};
  void PutPlane(real x1,real y1,real x2,real y2);
  void PutPlane(GeomVector& v1,GeomVector& v2);
  void PutPlane(GeomPlane pl);
  GeomPolygon GetReflectedCopy(real xvec,real yvec,int d_regid=0);
  GeomPolygon GetRotatedCopy(real rot,int d_regid=0);
  void Shift(real dx, real dy);
};

class GeomDividedCircle:public GeomGeneral{
 protected:
  real r;
  int nline;
  vector<GeomPlane> line;
  bool quadrant[8]; // For 1/8 circle
  bool octant;
  void PutHalfCircle(real x,real y,real r,real x1,real y1,real x2,real y2);
 public:
  GeomDividedCircle();
  void PutNline(int i);
  void PutR(real i);
  void PutQuarterCircle(real x,real y,real r,int quad);
  void PutHalfQuarterCircle(real x,real y,real r,int quad);
  void PutTopHalfCircle(real x,real y,real r);
  void PutBottomHalfCircle(real x,real y,real r);
  void PutRightHalfCircle(real x,real y,real r);
  void PutRightBottomHalfCircle(real x,real y,real r);
  void PutLeftBottomHalfCircle(real x,real y,real r);  
  void PutSpecificHalfCircle(real x,real y,real r);
  void PutLeftHalfCircle(real x,real y,real r);
  void PutPlane(GeomPlane pl);
  void OctantTrue(){octant=true;};
  void QuadrantTrue(int i){quadrant[i]=true;};
  void QuadrantFalse(int i){quadrant[i]=false;};
  GeomPlane CrossPlane(GeomPlane &pl);
  GeomDividedCircle GetReflectedCopy(real xvec,real yvec,int d_regid=0);
  GeomDividedCircle GetRotatedCopy(real rot,int d_regid=0);
  void Shift(real dx, real dy);
};

class IrregularGeometryInformation{
 protected:
  vector<GeomPolygon> pol;
  vector<GeomCircle> cir;
  vector<GeomDividedCircle> qcir;
  vector<GeomKind> gek;
  vector<unsigned> geid;
  vector<bool> macbnd; // flag of macro-cell boundary or not (for C5G7)
  int icir,ipol,iqcir,geomn,Region;
public:
  IrregularGeometryInformation();
  void AddGeom(GeomCircle &inp,bool macb=true);
  void AddGeom(GeomPolygon &inp,bool macb=true);
  void AddGeom(GeomDividedCircle &inp,bool macb=true);
  void AddHexagonLayer(int n,real *a,int *rid,real x=0.0,real y=0.0);
  void AddHexagonRing(real a, int npin, int pinl, real *r, int *rid, real x=0.0, real y=0.0);
  void AddHexagonRingFine(real a, int npin, int pinl, real *r, int *rid, real x=0.0, real y=0.0);
  void AddCircleRing(int n,real *r,int *rid,real x=0.0,real y=0.0);
  void AddCircleRing(real x, real y,int n,real *r,int &meshid);
  void AddHalfCircleRingTop(real x,real y,int n,real *r,int &meshid);
  void AddHalfCircleRingRight(real x,real y,int n,real *r,int &meshid);  
  void AddHalfCircleRingRightBottom(real x,real y,int n,real *r,int &meshid);
  void AddHalfCircleRingLeftBottom(real x,real y,int n,real *r,int &meshid);  
  void AddCircleRing(int n,vector<real> r,int *rid,real x=0.0,real y=0.0);
  void AddHalfQuarterCircleRing(real x,real y,int ring, real *rr,int quad,int& meshid);
  void AddHalfQuarterCircleRingAll(real x,real y,int ring, real *rr,int& meshid);
  void AddHalfQuarterCircleRingTopHalf(real x,real y,int ring, real *rr,int& meshid);
  void AddHalfQuarterCircleRingRightHalf(real x,real y,int ring, real *rr,int& meshid);  
  void AddHalfQuarterCircleRingRightBottomHalf(real x,real y,int ring, real *rr,int& meshid);
  void AddHalfQuarterCircleRingLeftBottomHalf(real x,real y,int ring, real *rr,int& meshid);
  real GetEdge();
  GeomKind GetGek(int i){return gek[i];};
  int GetGeid(int i){return geid[i];};
  int GetGeomn(){return geomn;};
  int GetRegion(){return Region;};
  void PutRegion(int inp){Region=inp;};
  real GetSurface();
  int GetIcir(){return icir;};
  int GetIpol(){return ipol;};
  int GetRegionID(int i);
  GeomCircle &GetCircle(int i){return cir[i];};
  GeomPolygon &GetPolygon(int i){return pol[i];};
  GeomVector SymmVec(GeomVector &inp);
  void DrawTrajectory
  (GeomPlane &pl1, real &xsi_stt, real &xsi_end,
   int &no_segment, vector<real> &length_segment, vector<int> &regid_segment);
  void WriteGnuplotFile(real deltax, int reg_st=0, int reg_ed=9999);
  // +++ Default geometry
  void ExHexagon();
  void ExMonjuHex();
  void ExMonjuSqr();
  void ExCir();
  void ExCircle(int div,real *r,int *med);
  void PutLWRCell(real x,real y,real pitch,real rpel,real rcld,int dpel,int dmod,int &regid,bool triangular=false);
  void PutLWRCellNoDivision(real x,real y,real pitch,real rpel,real rcld,int &regid);
  void PutLWRCellNoDivisionHalfQuarter(real x,real y,real pitch,real rpel,real rcld,int &regid);
  void PutSquareCell(real x,real y,real pitch,int &regid);
  void PutSquareCellHalfQuarter(real x,real y,real pitch,int &regid);
  void PutSquareDividedCell(real x,real y,real pitch,int &regid,bool triangular=false);
  void PutTriangleAssembly(real x,real y,real len,int layer,int &regid);
  void PutTriangleAssemblySquare(real x,real y,real lenx,real leny,int layerx,int layery,int &regid,bool quadadd=false);
  void PutTriangleAssemblyTriangle(real x,real y,real lenx,real leny,int layerx,int layery,int &regid,bool quadadd=false);  
  //
  void ExpandingHalfQuarterGeometry(bool new_regid=false);
  void ExpandingQuarterGeometry(bool new_regid=false);  
  void ExpandingHalfQuarterGeometryNew(bool new_regid=false);
  void ShiftGeometry(real dx,real dy);
  void RotatingGeometry(real rot); // [rot] in unit of 2PI radian: in the case of [0.5], 180oC rotation 
};

#endif
