#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "IrregularGeometryInformation.h"

using namespace std;

/** IrregularGeometryInformation **/

IrregularGeometryInformation::IrregularGeometryInformation()
{
  geomn=0;
  icir=0;
  ipol=0;
  iqcir=0;

  Region=0;
};

int IrregularGeometryInformation::GetRegionID(int i)
{
  if(gek[i]==Circle){
    return cir[geid[i]].GetRegionID();
  }else if(gek[i]==Polygon){
    return pol[geid[i]].GetRegionID();
  }else if(gek[i]==QCircle){
    return qcir[geid[i]].GetRegionID();
  }else{
    cout<<"Error in GetRegionID.\n";
    exit(0);
  };
};

void IrregularGeometryInformation::AddGeom(GeomCircle &inp,bool macb)
{
  //cir[icir]=inp;
  //gek[geomn]=Circle;
  //geid[geomn]=icir;
  cir.push_back(inp);
  gek.push_back(Circle);
  geid.push_back(icir);
  macbnd.push_back(macb);
  icir++;
  //if(icir==MAX_GEOM)cout<<"Too many circle!\n";
  geomn++;
  //if(geomn==MAX_GEOM)cout<<"Too many Geometry!\n";
  int tmp=inp.GetRegionID();
  if(tmp<0){
    cout<<"You should put region ID for geometry.\n";
    exit(0);
  };
  if(tmp+1>Region)PutRegion(tmp+1);
};

void IrregularGeometryInformation::AddGeom(GeomPolygon &inp,bool macb)
{
  // Check for background polygon
  if(geomn==0){
    real x=inp.GetCenter().getx();
    real y=inp.GetCenter().gety();
    if(x!=0.||y!=0.){
      cout<<"# +++ Warning !!! +++++++++++++++++++++++++\n";
      cout<<"# In IrregularGeometryInformation::AddGeom.\n";
      cout<<"# Background polygon should be put (0,0).  \n";
      cout<<"# +++++++++++++++++++++++++++++++++++++++++\n";
    };
  };
  //pol[ipol]=inp;
  //gek[geomn]=Polygon;
  //geid[geomn]=ipol;
  pol.push_back(inp);
  gek.push_back(Polygon);
  geid.push_back(ipol);
  macbnd.push_back(macb);
  ipol++;
  //if(ipol==MAX_GEOM)cout<<"Too many Polygon!\n";
  geomn++;
  //if(geomn==MAX_GEOM)cout<<"Too many Geometry!\n";
  int tmp=inp.GetRegionID();
  if(tmp<0){
    cout<<"You should put region ID for geometry.\n";
    exit(0);
  };
  if(tmp+1>Region)PutRegion(tmp+1);
};

void IrregularGeometryInformation::AddGeom(GeomDividedCircle &inp,bool macb)
{
  //qcir[iqcir]=inp;
  //gek[geomn]=QCircle;
  //geid[geomn]=iqcir;
  qcir.push_back(inp);
  gek.push_back(QCircle);
  geid.push_back(iqcir);
  macbnd.push_back(macb);
  iqcir++;
  //if(iqcir==MAX_GEOM)cout<<"Too many QuarterCircle!\n";
  geomn++;
  //if(geomn==MAX_GEOM)cout<<"Too many Geometry!\n";
  int tmp=inp.GetRegionID();
  if(tmp<0){
    cout<<"You should put region ID for geometry.\n";
    exit(0);
  };
  if(tmp+1>Region)PutRegion(tmp+1);
};

void IrregularGeometryInformation::AddHexagonLayer(int n,real *a,int *med,real x,real y)
{
  GeomPolygon h;
  for(int i=0;i<n;i++){
    h.PutHexagon(x,y,a[i]);
    h.PutRegionID(med[i]);
    AddGeom(h);
  };
};

void IrregularGeometryInformation::AddHexagonRing(real a,int npin,int pinl,real *r,int *med, real x,real y)
{
  real h=a*0.5*sqrt(3.0);
  real pp=1.0/npin;

  GeomVector p(x-a*0.5,y+h);
  GeomVector v[6];
  v[0].put_data(a,0);
  v[1].put_data(a*0.5,-h);
  v[2].put_data(-a*0.5,-h);
  v[3].put_data(-a,0);
  v[4].put_data(-a*0.5,h);
  v[5].put_data(a*0.5,h);
  for(int i=0;i<6;i++){
    for(int j=0;j<npin;j++){
      AddCircleRing(pinl, r, med, p.getx(), p.gety());
      p=p+v[i]*pp;
    };
  };
};

void IrregularGeometryInformation::AddHexagonRingFine(real a,int npin,int pinl,real *r,int *med, real x,real y)
{
  real h=a*0.5*sqrt(3.0);
  real pp=1.0/npin;

  int *med2=new int[pinl];

  GeomVector p(x-a*0.5,y+h);
  GeomVector v[6];
  v[0].put_data(a,0);
  v[1].put_data(a*0.5,-h);
  v[2].put_data(-a*0.5,-h);
  v[3].put_data(-a,0);
  v[4].put_data(-a*0.5,h);
  v[5].put_data(a*0.5,h);
  for(int i=0;i<6;i++){
    int id=0;
    for(int j=0;j<npin;j++){
      for(int l=0;l<pinl;l++)med2[l]=med[l]+id*pinl;
      AddCircleRing(pinl, r, med2, p.getx(), p.gety());
      p=p+v[i]*pp;
      if(npin%2==0){
        if((j+2)*2<=npin+2){
          id++;
	}else{
	  id--;
	};
      }else{
	if((j+2)*2<=npin+1)id++;
        if((j+1)*2>npin+1)id--;
      };
    };
  };

  delete [] med2;
};

void IrregularGeometryInformation::AddCircleRing(int n,real *r,int *med,real x,real y)
{
  GeomCircle c;
  for(int i=0;i<n;i++){
    c.Init(x,y,r[i]);
    c.PutRegionID(med[i]);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddCircleRing(real x,real y,int ring,real *r,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomCircle c;
    c.Init(x,y,r[ring-1-i]);
    c.PutRegionID(meshid++);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddHalfCircleRingTop(real x,real y,int ring,real *r,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomDividedCircle c;
    c.PutTopHalfCircle(x,y,r[ring-1-i]);
    c.PutRegionID(meshid++);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddHalfCircleRingRight(real x,real y,int ring,real *r,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomDividedCircle c;
    c.PutRightHalfCircle(x,y,r[ring-1-i]);
    c.PutRegionID(meshid++);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddHalfCircleRingRightBottom(real x,real y,int ring,real *r,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomDividedCircle c;
    c.PutRightBottomHalfCircle(x,y,r[ring-1-i]);
    c.PutRegionID(meshid++);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddHalfCircleRingLeftBottom(real x,real y,int ring,real *r,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomDividedCircle c;
    c.PutLeftBottomHalfCircle(x,y,r[ring-1-i]);
    c.PutRegionID(meshid++);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddCircleRing(int n,vector<real> r,int *med,real x,real y)
{
  GeomCircle c;
  for(int i=0;i<n;i++){
    c.Init(x,y,r[i]);
    c.PutRegionID(med[i]);
    AddGeom(c);
  };
};

void IrregularGeometryInformation::AddHalfQuarterCircleRing(real x,real y,int ring,real *rr,int quad,int& meshid)
{
  for(int i=0;i<ring;i++){
    GeomDividedCircle cir;
    cir.PutHalfQuarterCircle(x,y,rr[ring-1-i],quad);
    cir.PutRegionID(meshid++);
    AddGeom(cir);
  };
};

void IrregularGeometryInformation::AddHalfQuarterCircleRingAll(real x,real y,int ring,real *rr,int& meshid)
{
  for(int i=0;i<8;i++){
    AddHalfQuarterCircleRing(x,y,ring,rr,i+1,meshid);
  };
};

void IrregularGeometryInformation::AddHalfQuarterCircleRingTopHalf(real x,real y,int ring,real *rr,int& meshid)
{
  for(int i=0;i<4;i++){
    AddHalfQuarterCircleRing(x,y,ring,rr,4-i,meshid);
  };
};

void IrregularGeometryInformation::AddHalfQuarterCircleRingRightHalf(real x,real y,int ring,real *rr,int& meshid)
{
  AddHalfQuarterCircleRing(x,y,ring,rr,7,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,8,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,1,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,2,meshid);  
};

void IrregularGeometryInformation::AddHalfQuarterCircleRingRightBottomHalf(real x,real y,int ring,real *rr,int& meshid)
{
  AddHalfQuarterCircleRing(x,y,ring,rr,6,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,7,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,8,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,1,meshid);
};

void IrregularGeometryInformation::AddHalfQuarterCircleRingLeftBottomHalf(real x,real y,int ring,real *rr,int& meshid)
{
  AddHalfQuarterCircleRing(x,y,ring,rr,4,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,5,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,6,meshid);
  AddHalfQuarterCircleRing(x,y,ring,rr,7,meshid);
};

void IrregularGeometryInformation::ExHexagon()
{
  Region=23;
 
  real scon=2.0/sqrt(3.0);
  int medi[]={22,21, 20,19,18,17,16,15,14};
  real hring[]={scon*4.075,scon*3.925,scon*3.735,scon*3.36191,
                 scon*2.80159,scon*2.24127,scon*1.68095,scon*1.12064,scon*0.56032};

  AddHexagonLayer(9,hring,medi);

  real cring[]={0.275,0.2315};
  int medj[]={1,0};
  AddCircleRing(2,cring,medj);
  int medj2[]={3,2};
  AddHexagonRing(scon*0.56032,1,2,cring,medj2);
  int medj3[]={5,4};
  AddHexagonRing(scon*1.12064,2,2,cring,medj3);
  int medj4[]={7,6};
  AddHexagonRing(scon*1.68095,3,2,cring,medj4);
  int medj5[]={9,8};
  AddHexagonRing(scon*2.24127,4,2,cring,medj5);
  int medj6[]={11,10};
  AddHexagonRing(scon*2.80159,5,2,cring,medj6);
  int medj7[]={13,12};
  AddHexagonRing(scon*3.36191,6,2,cring,medj7);
}

void IrregularGeometryInformation::ExMonjuHex()
{
  Region=3;
  GeomPolygon h;
  h.PutHexagon(0,0,0.66532);
  h.PutRegionID(2);
  AddGeom(h);

  real rinp[]={0.325,0.278};
  int medi[]={1,0};
  AddCircleRing(2,rinp,medi);
}

void IrregularGeometryInformation::ExCir()
{
  Region=4;

  GeomCircle c[6];
  c[0].Init(0,0,2);
  c[0].PutRegionID(0);
  c[1].Init(0,0,1);
  c[1].PutRegionID(1);
  c[2].Init(0,1.5,0.4);
  c[2].PutRegionID(2);
  c[3].Init(0,-1.5,0.4);
  c[3].PutRegionID(2);
  c[4].Init(1.5,0,0.4);
  c[4].PutRegionID(3);
  c[5].Init(-1.5,0,0.4);
  c[5].PutRegionID(3);

  for(int i=0;i<6;i++){
    AddGeom(c[i]);
  };
};

void IrregularGeometryInformation::ExMonjuSqr()
{
  Region=3;
  GeomPolygon s;
  s.PutRectangular(0,0,0.5362,0.5362);
  s.PutRegionID(2);
  AddGeom(s);

  real rinp[]={0.325,0.278};
  int medi[]={1,0};
  AddCircleRing(2,rinp,medi);
};

void IrregularGeometryInformation::ExCircle(int div,real *r,int *med)
{
  Region=div;
  AddCircleRing(div,r,med);
};

GeomVector IrregularGeometryInformation::SymmVec(GeomVector &inp)
{
  GeomVector ret;
  int ig=gek[0];
  if(ig==Polygon)ret=pol[0].SymmVec(inp);
  return ret;
};

real IrregularGeometryInformation::GetEdge()
{
  //real factor=1.01;
  real factor=1.0;
  if(gek[0]==Circle)return cir[0].GetR()*factor;
  if(gek[0]==Polygon){
    int numpl=pol[0].GetNumpl();
    real max_dist_x=0.;
    real max_dist_y=0.;    
    for(int j=0;j<numpl;j++){
      real x=fabs(pol[0].GetPlane(j).get_p1().getx());
      if(x>max_dist_x)max_dist_x=x;
      real y=fabs(pol[0].GetPlane(j).get_p1().gety());
      if(y>max_dist_y)max_dist_y=y;
      x=fabs(pol[0].GetPlane(j).get_p2().getx());
      if(x>max_dist_x)max_dist_x=x;
      y=fabs(pol[0].GetPlane(j).get_p2().gety());
      if(y>max_dist_y)max_dist_y=y;
    };
    if(max_dist_y>max_dist_x)max_dist_x=max_dist_y;
    return max_dist_x*factor;
    /*
    if(pol[0].GetKind()==1){
      real a=pol[0].GetPlane(0).get_long()*0.5;
      real b=pol[0].GetPlane(1).get_long()*0.5;
      if(a>b){return a*factor;}else{return b*factor;};
    }else if(pol[0].GetKind()==2){ // hexagon
      return pol[0].GetPlane(0).get_long()*factor;
    }else if(pol[0].GetKind()==3){ // specific triangle
      return pol[0].GetPlane(0).get_long()*0.5*factor;
    };
    */
    
  };
  return 0.0;
};

real IrregularGeometryInformation::GetSurface()
{
  if(gek[0]==Circle){
    return cir[0].GetR()*PI2;
  };
  if(gek[0]==Polygon){
    return pol[0].GetSurface();
  };
  return 0.0;
};		 

void IrregularGeometryInformation::DrawTrajectory
(GeomPlane &pl1, real &xsi_stt, real &xsi_end,
 int &no_segment, vector<real> &length_segment, vector<int>& regid_segment)
{
  vector<real> xsi;
  vector<int> genm;

  real planelong=pl1.get_long();

  vector<int> medk(geomn);

  int ipoint=0;
  GeomPlane ret;

  bool macskip=false;
  for(int j=0;j<geomn;j++){

      if(macbnd[j]||!macskip){

        if(gek[j]==Circle){
          medk[j]=cir[geid[j]].GetRegionID();
          ret=cir[geid[j]].CrossPlane(pl1);
        }else if(gek[j]==Polygon){
          medk[j]=pol[geid[j]].GetRegionID();
          ret=pol[geid[j]].CrossPlane(pl1);
        }else if(gek[j]==QCircle){
          medk[j]=qcir[geid[j]].GetRegionID();
          ret=qcir[geid[j]].CrossPlane(pl1);
        };

        if(ret.GetExist()){
          macskip=false;
          xsi.push_back(pl1.get_xsi(ret.get_p1()));
          xsi.push_back(pl1.get_xsi(ret.get_p2()));
	  if(j==0){
            real factor=1.;
	    if(xsi[ipoint]>xsi[ipoint+1])factor=-1;
	    xsi[ipoint]-=0.00000001*factor;
	    xsi[ipoint+1]+=0.00000001*factor;
	  };
          genm.push_back(j);
          genm.push_back(j);
          ipoint+=2;
        }else{
	  if(macbnd[j])macskip=true;
	};

      };

  };

  if(ipoint==0){
    no_segment=0;
    return;
  };

  vector<int> re(ipoint);

  ChangeOrder(xsi,ipoint,re);

  bool open_background_geom=false;

  xsi_stt=xsi[re[0]];
  //xsi_end=xsi[re[ipoint-1]];

  vector<bool> Open(geomn,false); // initial:close
  int LayTop=genm[re[0]]; // Top-plane-ID
  Open[LayTop]=true; // open

  if(LayTop==0){
    open_background_geom=true;
  };

  no_segment=0;
  length_segment.clear();
  regid_segment.clear();

  bool end_loop=false;
  for(int i=1;i<ipoint;i++){
    int rr=genm[re[i]];
    if(rr==0){
      if(!open_background_geom){
        xsi_stt=xsi[re[i]];
        open_background_geom=true;
      }else{
	end_loop=true;
      };
    };
    if(!Open[rr]){ // close
      if(rr>=LayTop){
	if(open_background_geom){
          regid_segment.push_back(medk[LayTop]);
          length_segment.push_back((xsi[re[i]]-xsi_stt)*0.5*planelong);
	  no_segment++;
	  xsi_end=xsi[re[i]];
	};
	LayTop=rr;
      };
      Open[rr]=true;
    }else{
      if(rr==LayTop||rr==0){
	int maxlay=0;
	for(int j=0;j<geomn;j++){
	  if(rr!=j){
	    if(Open[j]&&j>maxlay)maxlay=j;
	  };
	};
	LayTop=maxlay; 
        if(open_background_geom){
          regid_segment.push_back(medk[rr]);
          length_segment.push_back((xsi[re[i]]-xsi_stt)*0.5*planelong);
          no_segment++;
	  xsi_end=xsi[re[i]];
	};
      };
      Open[rr]=false;
    };
    if(end_loop)i=ipoint;
  };

  real Judge=0.000001;
  //If segment is shorter that this value, the segment is ignored. 
  for(int i=-1;i<no_segment-1;i++){
    real diff;
    if(i==-1){diff=length_segment[0];}
    else{diff=length_segment[i+1]-length_segment[i];};
    if(diff<Judge){
      for(int j=i+1;j<no_segment-1;j++){
        length_segment[j]=length_segment[j+1];
	regid_segment[j]=regid_segment[j+1];
      };
      i--;
      no_segment--;
    };
  };

  for(int i=0;i<no_segment-1;i++){
    if(regid_segment[i+1]==regid_segment[i]){
      for(int j=i;j<no_segment-1;j++){
	length_segment[j]=length_segment[j+1];
	regid_segment[j]=regid_segment[j+1];
      };
      i--;
      no_segment--;
    };
  };

};

void IrregularGeometryInformation::WriteGnuplotFile(real deltax, int reg_st, int reg_ed)
{
  ofstream fout;
  string filename="plotdata";

  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  real sizex=GetEdge()*1.1;

  real max=sizex;
  vector<real> l;
  vector<int> id;

  real y=-max;
  while(y<max){
    GeomPlane pl(-max,y,max,y);
    real xsi1,xsi2;
    int n;
    DrawTrajectory(pl,xsi1,xsi2,n,l,id);
    if(n!=0){
      real x1=(-max*(1-xsi1)+max*(1+xsi1))*0.5;
      int reg0=id[0];
      if(reg0>=reg_st&&reg0<=reg_ed)fout<<" "<<x1<<" "<<y<<"\n";
      for(int i=0;i<n;i++){
        bool prt=false;
        int reg1=id[i];
	if(reg1>=reg_st&&reg1<=reg_ed)prt=true;
	if(i!=n-1){
	  int reg2=id[i+1];
  	  if(reg2>=reg_st&&reg2<=reg_ed)prt=true;
	};
        if(prt)fout<<" "<<x1+l[i]<<" "<<y<<"\n";
      };
    };
    y+=deltax;
  };

  real x=-max;
  while(x<max){
    GeomPlane pl(x,-max,x,max);
    real xsi1,xsi2;
    int n;
    DrawTrajectory(pl,xsi1,xsi2,n,l,id);
    if(n!=0){
      real y1=(-max*(1-xsi1)+max*(1+xsi1))*0.5;
      int reg0=id[0];
      if(reg0>=reg_st&&reg0<=reg_ed)fout<<" "<<x<<" "<<y1<<"\n";
      for(int i=0;i<n;i++){
        bool prt=false;
        int reg1=id[i];
	if(reg1>=reg_st&&reg1<=reg_ed)prt=true;
	if(i!=n-1){
	  int reg2=id[i+1];
  	  if(reg2>=reg_st&&reg2<=reg_ed)prt=true;
	};
        if(prt)fout<<" "<<x<<" "<<y1+l[i]<<"\n";
      };
    };
    x+=deltax;
  };

  fout.close();
};

void IrregularGeometryInformation::ExpandingHalfQuarterGeometry(bool new_regid)
{
  real xvecd[]={1.,0.,1.};
  real yvecd[]={1.,1.,0.};

  for(int k=0;k<3;k++){

    int d_regid=0; // Delta region-ID
    if(new_regid)d_regid=Region;

    real xvec=xvecd[k];
    real yvec=yvecd[k];

    int max=geomn;
    for(int i=1;i<max;i++){
      int id=geid[i];
      if(gek[i]==Circle){
        GeomCircle newcir=cir[id].GetReflectedCopy(xvec,yvec,d_regid);
        AddGeom(newcir);
      }else if(gek[i]==Polygon){
      GeomPolygon newpol=pol[id].GetReflectedCopy(xvec,yvec,d_regid);
      AddGeom(newpol);
      }else{
      GeomDividedCircle newcir=qcir[id].GetReflectedCopy(xvec,yvec,d_regid);
      AddGeom(newcir);
      };
    };

  };
};

void IrregularGeometryInformation::ExpandingQuarterGeometry(bool new_regid)
{
  real xvecd[]={1.,0.,1.};
  real yvecd[]={1.,1.,0.};

  for(int k=1;k<3;k++){

    int d_regid=0; // Delta region-ID
    if(new_regid)d_regid=Region;

    real xvec=xvecd[k];
    real yvec=yvecd[k];

    int max=geomn;
    for(int i=1;i<max;i++){
      int id=geid[i];
      if(gek[i]==Circle){
        GeomCircle newcir=cir[id].GetReflectedCopy(xvec,yvec,d_regid);
        AddGeom(newcir);
      }else if(gek[i]==Polygon){
      GeomPolygon newpol=pol[id].GetReflectedCopy(xvec,yvec,d_regid);
      AddGeom(newpol);
      }else{
      GeomDividedCircle newcir=qcir[id].GetReflectedCopy(xvec,yvec,d_regid);
      AddGeom(newcir);
      };
    };

  };
};

void IrregularGeometryInformation::ShiftGeometry(real dx,real dy)
{
  for(int i=0;i<icir;i++){
    cir[i].Shift(dx,dy);
  };
  for(int i=0;i<ipol;i++){
    pol[i].Shift(dx,dy);
  };
  for(int i=0;i<iqcir;i++){
    qcir[i].Shift(dx,dy);
  };
};

void IrregularGeometryInformation::RotatingGeometry(real rot)
{
  for(int i=0;i<icir;i++){
    cir[i]=cir[i].GetRotatedCopy(rot);
  };
  for(int i=0;i<ipol;i++){
    pol[i]=pol[i].GetRotatedCopy(rot);
  };
  for(int i=0;i<iqcir;i++){
    qcir[i]=qcir[i].GetRotatedCopy(rot);
  };
};

void IrregularGeometryInformation::ExpandingHalfQuarterGeometryNew(bool new_regid)
{
  int d_regid=0; // Delta region-ID
  if(new_regid)d_regid=Region;

  // first 45-degree symmetry
  int max=geomn;
  for(int i=1;i<max;i++){
    int id=geid[i];
    if(gek[i]==Circle){
      GeomCircle newcir=cir[id].GetReflectedCopy(1.,1.,d_regid);
      AddGeom(newcir);
    }else if(gek[i]==Polygon){
      GeomPolygon newpol=pol[id].GetReflectedCopy(1.,1.,d_regid);
      AddGeom(newpol);
    }else{
      GeomDividedCircle newcir=qcir[id].GetReflectedCopy(1.,1.,d_regid);
      AddGeom(newcir);
    };
  };

  real rot[]={0.25, 0.5};
  for(int ii=0;ii<2;ii++){

    int d_regid=0; // Delta region-ID
    if(new_regid)d_regid=Region;

    int max=geomn;
    for(int i=1;i<max;i++){
      int id=geid[i];
      if(gek[i]==Circle){
        GeomCircle newcir=cir[id].GetRotatedCopy(rot[ii],d_regid);
        AddGeom(newcir);
      }else if(gek[i]==Polygon){
        GeomPolygon newpol=pol[id].GetRotatedCopy(rot[ii],d_regid);
        AddGeom(newpol);
      }else{
        GeomDividedCircle newcir=qcir[id].GetRotatedCopy(rot[ii],d_regid);
        AddGeom(newcir);
      };
    };
  };

};

// +++++++++++++++++++
// +++ GeomGeneral +++
// +++++++++++++++++++

GeomGeneral::GeomGeneral()
{
  RegionID=-1;
  outer_circle_r=0.;
};

void GeomGeneral::PutCenter(real x,real y)
{
  center.Initialize(x,y);
};

void GeomGeneral::PutCenter(GeomVector &vinp)
{
  center=vinp;
};

// ++++++++++++++++++
// +++ GeomCircle +++
// ++++++++++++++++++

void GeomCircle::Init(real xinp,real yinp,real rinp)
{
  center.put_data(xinp,yinp);
  PutR(rinp);
};

void GeomCircle::PutR(real i)
{
  r=i;
  outer_circle_r=i;
};

GeomPlane GeomCircle::CrossPlane(GeomPlane &pl)
{
  GeomPlane ret;

  real dist=pl.CalcDistance(center);
  if(dist>outer_circle_r)return ret;

  real xx=center.getx();
  real yy=center.gety();
  real rr=r;
  real r2=rr*rr;

  real x1,x2,y1,y2,dd;
  real a=pl.get_a();
  real b=pl.get_b();
  real c=pl.get_c();
  if(a!=0.0){
    /*
    real inv_a=1./a;
    real inv_a2=inv_a*inv_a;
    real alp=c*inv_a+xx;
    real aa=b*b*inv_a2+1.0;
    real bb=2*b*alp*inv_a-2*yy;
    real cc=yy*yy+alp*alp-r2;
    real abc=bb*bb-4*aa*cc;
    if(abc<=0){dd=0;}
    else{dd=sqrt(abc);};
    real half_aa_inv=0.5/aa;
    y1=(-bb+dd)*half_aa_inv;
    y2=(-bb-dd)*half_aa_inv;
    x1=-b*inv_a*y1-c*inv_a;
    x2=-b*inv_a*y2-c*inv_a;
    */
    real alpha=-a*xx-b*yy-c;
    real aa=a*a+b*b;
    real bb=-alpha*b;
    real cc=alpha*alpha-a*a*r2;
    real abc=bb*bb-aa*cc;
    if(abc<=0){dd=0.;}
    else{dd=sqrt(abc);};
    real aainv=1./aa;
    y1=(-bb+dd)*aainv+yy;
    y2=(-bb-dd)*aainv+yy;
    real inv_a=1./a;
    x1=-b*inv_a*y1-c*inv_a;
    x2=-b*inv_a*y2-c*inv_a;
  }
  else{
    real tmp=1./b;
    real c_bar_b=c*tmp;
    real ee=sqrt(r2-(c_bar_b+yy)*(c_bar_b+yy));
    x1=xx+ee;
    x2=xx-ee;
    y1=-c_bar_b;
    y2=y1;
  };
  ret.PutVector(x1,y1,x2,y2);
  return ret;
};

GeomCircle GeomCircle::GetReflectedCopy(real xvec,real yvec,int d_regid)
{
  GeomVector u(-yvec,xvec);
  GeomVector center2=center.OperatingElementaryReflector(u);

  GeomCircle ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);
  ret.PutR(r);

  return ret;
};

GeomCircle GeomCircle::GetRotatedCopy(real rot,int d_regid)
{
  GeomVector center2=center.OperatingRotator(rot);
  
  GeomCircle ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);
  ret.PutR(r);

  return ret;
};

void GeomCircle::Shift(real dx, real dy)
{
  center.Shift(dx,dy);
};

// +++++++++++++++++++
// +++ GeomPolygon +++
// +++++++++++++++++++

void GeomPolygon::PutNumpl(int i)
{
  numpl=i;
  plane.resize(numpl);
};

void GeomPolygon::PutRectangular(real x,real y,real a,real b)
{
  center.put_data(x,y);

  PutNumpl(4);
  PutKind(1);

  plane[0].PutVector(x-a,y-b,x+a,y-b);
  plane[1].PutVector(x+a,y-b,x+a,y+b);
  plane[2].PutVector(x+a,y+b,x-a,y+b);
  plane[3].PutVector(x-a,y+b,x-a,y-b);

  real rinp=sqrt(a*a+b*b);
  outer_circle_r=rinp;
};

void GeomPolygon::PutHexagon(real x,real y,real a)
{
  center.put_data(x,y);

  PutNumpl(6);
  PutKind(2);

  real r3a=a*sqrt(3.0)*0.5;
  real a05=a*0.5;
  plane[0].PutVector(x+a05,y+r3a,x-a05,y+r3a);
  plane[1].PutVector(x-a05,y+r3a,x-a,y);
  plane[2].PutVector(x-a,y,x-a05,y-r3a);
  plane[3].PutVector(x-a05,y-r3a,x+a05,y-r3a);
  plane[4].PutVector(x+a05,y-r3a,x+a,y);
  plane[5].PutVector(x+a,y,x+a05,y+r3a);

  outer_circle_r=a;
};

void GeomPolygon::PutTriangleHalf(real x,real y,real len)
{
  center.put_data(x,y);

  PutNumpl(3);
  PutKind(4); // dummy

  outer_circle_r=len*0.5*sqrt(2.);

  real hl=len*0.5;  
  plane[0].PutVector(x-hl,y-hl,x+hl,y-hl);
  plane[1].PutVector(x+hl,y-hl,x,y);
  plane[2].PutVector(x,y,x-hl,y-hl);    
};

void GeomPolygon::PutTriangle(real x,real y,real len,int quad)
{
  if(quad<1||quad>4){
    cout<<"# Error in GeomPolygon::PutTriangle.\n";
    cout<<"# You have to set [quad] from 1 to 4.\n";
    cout<<"# You set [quad] as "<<quad<<"\n";
    exit(0);
  };
  center.put_data(x,y);

  PutNumpl(3);
  PutKind(4);
  if(quad==4)PutKind(3);

  outer_circle_r=len*0.5*sqrt(2.);

  real hl=len*0.5;
  vector<real> xx(5);
  vector<real> yy(5);
  xx[0]=x+hl;
  xx[1]=x-hl;
  xx[2]=x-hl;
  xx[3]=x+hl;
  xx[4]=xx[0];
  yy[0]=y+hl;
  yy[1]=y+hl;
  yy[2]=y-hl;
  yy[3]=y-hl;
  yy[4]=yy[0];
  
  int ind=quad-1;
  int ind2=ind+1;
  for(int i=0;i<3;i++){
    if(ind>=4)ind-=4;
    if(ind2>=4)ind2-=4;
    plane[i].PutVector(xx[ind],yy[ind],xx[ind2],yy[ind2]);
    ind++;
    ind2++;
    if(i==0){
      ind2++;
    };
    if(i==1){
      ind++;
    };
  };
};

void GeomPolygon::PutTriangle(real x,real y,real lenx,real leny,int quad)
{
  if(quad<1||quad>4){
    cout<<"# Error in GeomPolygon::PutTriangle.\n";
    cout<<"# You have to set [quad] from 1 to 4.\n";
    cout<<"# You set [quad] as "<<quad<<"\n";
    exit(0);
  };
  center.put_data(x,y);

  PutNumpl(3);
  PutKind(4);
  if(quad==4)PutKind(3);

  real hlx=lenx*0.5;
  real hly=leny*0.5;

  outer_circle_r=sqrt(hlx*hlx+hly*hly);

  vector<real> xx(5);
  vector<real> yy(5);
  xx[0]=x+hlx;
  xx[1]=x-hlx;
  xx[2]=x-hlx;
  xx[3]=x+hlx;
  xx[4]=xx[0];
  yy[0]=y+hly;
  yy[1]=y+hly;
  yy[2]=y-hly;
  yy[3]=y-hly;
  yy[4]=yy[0];
  
  int ind=quad-1;
  int ind2=ind+1;
  for(int i=0;i<3;i++){
    if(ind>=4)ind-=4;
    if(ind2>=4)ind2-=4;
    plane[i].PutVector(xx[ind],yy[ind],xx[ind2],yy[ind2]);
    ind++;
    ind2++;
    if(i==0){
      ind2++;
    };
    if(i==1){
      ind++;
    };
  };
};

void GeomPolygon::PutEquilateralTriangle(real x,real y,real len,real angle)
{
  center.put_data(x,y);
  real cirr=sqrt(3.)*0.333333333*len;
  outer_circle_r=cirr;

  PutNumpl(3);
  PutKind(3); // No-use

  real theta=-0.166666666666*PI+angle/180.*PI;
  real x1=cirr*cos(theta);
  real y1=cirr*sin(theta);
  for(int i=0;i<3;i++){
    theta+=0.6666666666666*PI;
    real x2=cirr*cos(theta);
    real y2=cirr*sin(theta);
    plane[i].PutVector(x+x1,y+y1,x+x2,y+y2);
    x1=x2;
    y1=y2;
  };
};

GeomPlane GeomPolygon::CrossPlane(GeomPlane &pl)
{
  GeomPlane ret;

  if(outer_circle_r>0.){
    real dist=pl.CalcDistance(center);
    if(dist>outer_circle_r)return ret;
  };

  int ip=0;
  GeomVector pp[2];

  for(int i=0;i<numpl;i++){
    if(pl.BeCrossing(plane[i])){
      if(ip==0){
        pp[ip]=pl.CrossPlane(plane[i]);
        ip++;
      }else{
	int judge=0;
	GeomVector tmp=pl.CrossPlane(plane[i]);
        for(int j=0;j<ip;j++){
	  if(pp[j].SameOrNot(tmp)==1)judge=1;
	};
	if(judge==0){
          if(ip==2){
            cout<<"Error in CrossPlane in GeomSquare!\n";
	    cout<<"There are three crossing points.. ?\n";
	    cout<<" +++ Crossing points +++ \n";
            tmp.show_self();
	    pp[0].show_self();
	    pp[1].show_self();
	    cout<<" +++ GeomPlane +++\n";
	    for(int i=0;i<numpl;i++){
	      plane[i].show_self();
	    };
	    exit(0);
	  };
          pp[ip]=tmp;
          ip++;
	};
      };
    };
  };

  if(ip==2)ret.PutVector(pp[0],pp[1]);
  return ret;
};

GeomVector GeomPolygon::SymmVec(GeomVector &inp)
{
  GeomVector ret(0.,0.);
  real xsi;
  int iopp;

  if(kind==1){ // rectangular
    for(int i=0;i<4;i++){
      if(plane[i].VectorOnPlane(inp,0.00001)){
        xsi=plane[i].get_xsi(inp);
        if(fabs(xsi)<0.99999)xsi=-xsi;
        iopp=i+2;
        if(iopp>3)iopp-=4;
        ret=plane[iopp].GetVecXsi(xsi);
        return ret;
      };
    };
  }else if(kind==2){
    int ij=-1;
    for(int i=0;i<6;i++){
      if(plane[i].VectorOnPlane(inp,0.00001))ij=i;
      if(ij!=-1)break;
    };
    if(ij==-1)cout<<"# Error in SymmVec of GeomSquare!\n";
    xsi=plane[ij].get_xsi(inp);

    iopp=ij+3;
    if(iopp==6)iopp=0;
    if(iopp==7)iopp=1;
    if(iopp==8)iopp=2;
    GeomVector ret=plane[iopp].GetVecXsi(-xsi);
    return ret;
  };

  cout<<"# Error in SymmVec of GeomSquare!\n";
  cout<<"# kind : "<<kind<<"\n";
  inp.show_self();

  return ret;
};

real GeomPolygon::GetSurface()
{
  real ret=0.0;
  for(int i=0;i<numpl;i++){
    ret+=plane[i].get_long();
  };
  return ret;
};

void GeomPolygon::PutPlane(real x1,real y1,real x2,real y2)
{
  GeomPlane inp(x1,y1,x2,y2);
  plane.push_back(inp);
  numpl++;
};

void GeomPolygon::PutPlane(GeomVector &v1,GeomVector &v2)
{
  GeomPlane inp(v1,v2);
  plane.push_back(inp);
  numpl++;
};

void GeomPolygon::PutPlane(GeomPlane pl)
{
  plane.push_back(pl);
  numpl++;
};

GeomPolygon GeomPolygon::GetReflectedCopy(real xvec,real yvec,int d_regid)
{
  GeomVector u(-yvec,xvec);
  GeomVector center2=center.OperatingElementaryReflector(u);

  GeomPolygon ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);

  ret.PutKind(kind);
  for(int i=0;i<numpl;i++){
    ret.PutPlane(plane[i].OperatingElementaryReflector(u));
  };

  return ret;
};

GeomPolygon GeomPolygon::GetRotatedCopy(real rtot,int d_regid)
{
  GeomVector center2=center.OperatingRotator(rtot);

  GeomPolygon ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);

  ret.PutKind(kind);
  for(int i=0;i<numpl;i++){
    ret.PutPlane(plane[i].OperatingRotator(rtot));
  };

  return ret;
};

void GeomPolygon::Shift(real dx, real dy)
{
  center.Shift(dx,dy);
  for(int i=0;i<numpl;i++){
    plane[i].Shift(dx,dy);
  };
};

// +++++++++++++++++++++++++
// +++ GeomDividedCircle +++
// +++++++++++++++++++++++++

GeomDividedCircle::GeomDividedCircle():GeomGeneral()
{
  nline=0;
  octant=false;
};

void GeomDividedCircle::PutNline(int i){
  nline=i;
  line.resize(nline);
};

void GeomDividedCircle::PutR(real i)
{
  r=i;
  outer_circle_r=i;
};

void GeomDividedCircle::PutQuarterCircle(real xin,real yin,real rin,int quad)
{
  if(quad<1||quad>4){
    cout<<"Error in GeomDividedCircle::PutQuarterCircle.\n";
    cout<<"You have to put quadrant from 1 to 4.\n";
    cout<<"You set quadrant as "<<quad<<"\n";
    exit(0);
  };

  for(int i=0;i<8;i++){
    quadrant[i]=false;
  };
  quadrant[(quad-1)*2]=true;
  quadrant[(quad-1)*2+1]=true;

  center.put_data(xin,yin);
  PutR(rin);
  real x1=xin+r;
  real y2=yin+r;
  if(quad==2||quad==3){
    x1=xin-r;
  };
  if(quad==3||quad==4){
    y2=yin-r;
  };

  PutNline(2);
  line[0].PutVector(xin,yin,x1,yin);
  line[1].PutVector(xin,yin,xin,y2);
};

void GeomDividedCircle::PutHalfQuarterCircle(real xin,real yin,real rin,int quad)
{
  OctantTrue();

  if(quad<1||quad>8){
    cout<<"Error in GeomDividedCircle::PutHalfQuarterCircle.\n";
    cout<<"You have to put quadrant from 1 to 8.\n";
    cout<<"You set quadrant as "<<quad<<"\n";
    exit(0);
  };

  for(int i=0;i<8;i++){
    quadrant[i]=false;
  };
  quadrant[quad-1]=true;

  center.put_data(xin,yin);
  PutR(rin);
  real x1,y1,x2,y2;
  real fact=rin*sqrt(2.)*0.5;
  if(quad==1){
    x1=xin+r; y1=yin; x2=xin+fact; y2=yin+fact;
  }else if(quad==2){
    x1=xin+fact; y1=yin+fact; x2=xin; y2=yin+r;
  }else if(quad==3){
    x1=xin; y1=yin+r; x2=xin-fact; y2=yin+fact;
  }else if(quad==4){
    x1=xin-fact; y1=yin+fact; x2=xin-r; y2=yin;
  }else if(quad==5){
    x1=xin-r; y1=yin; x2=xin-fact; y2=yin-fact;
  }else if(quad==6){
    x1=xin-fact; y1=yin-fact; x2=xin; y2=yin-r;
  }else if(quad==7){
    x1=xin; y1=yin-r; x2=xin+fact; y2=yin-fact;
  }else{
    x1=xin+fact; y1=yin-fact; x2=xin+r; y2=yin;
  }

  PutNline(2);
  line[0].PutVector(xin,yin,x1,y1);
  line[1].PutVector(xin,yin,x2,y2);
};

void GeomDividedCircle::PutHalfCircle(real xin,real yin,real rin,real x1,real y1,real x2,real y2)
{
  center.put_data(xin,yin);
  PutR(rin);
  PutNline(1);
  line[0].PutVector(x1,y1,x2,y2);
};

void GeomDividedCircle::PutTopHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<4;i++){
    quadrant[i]=true;
    quadrant[4+i]=false;
  };
  PutHalfCircle(xin,yin,rin,xin+rin,yin,xin-rin,yin);
};

void GeomDividedCircle::PutBottomHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<4;i++){
    quadrant[i]=false;
    quadrant[4+i]=true;
  };
  PutHalfCircle(xin,yin,rin,xin+rin,yin,xin-rin,yin);
};

void GeomDividedCircle::PutSpecificHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<8;i++){
    quadrant[i]=true;
  };
  for(int i=1;i<5;i++){
    quadrant[i]=false;
  };
  real tmp=rin/sqrt(2.);
  PutHalfCircle(xin,yin,rin,xin+tmp,yin+tmp,xin-tmp,yin-tmp);
};

void GeomDividedCircle::PutRightHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<4;i++){
    quadrant[2+i]=false;
  };
  quadrant[0]=true;
  quadrant[1]=true;
  quadrant[6]=true;
  quadrant[7]=true;
  PutHalfCircle(xin,yin,rin,xin,yin+rin,xin,yin-rin);
};

void GeomDividedCircle::PutRightBottomHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<8;i++){
    quadrant[i]=false;
  };
  quadrant[0]=true;
  quadrant[5]=true;
  quadrant[6]=true;
  quadrant[7]=true;
  real tmp=rin/sqrt(2.);
  PutHalfCircle(xin,yin,rin,xin+tmp,yin+tmp,xin-tmp,yin-tmp);
};

void GeomDividedCircle::PutLeftBottomHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<8;i++){
    quadrant[i]=false;
  };
  quadrant[3]=true;
  quadrant[4]=true;
  quadrant[5]=true;
  quadrant[6]=true;
  real tmp=rin/sqrt(2.);
  PutHalfCircle(xin,yin,rin,xin+tmp,yin-tmp,xin-tmp,yin+tmp);
};

void GeomDividedCircle::PutLeftHalfCircle(real xin,real yin,real rin)
{
  for(int i=0;i<4;i++){
    quadrant[2+i]=true;
  };
  quadrant[0]=false;
  quadrant[1]=false;
  quadrant[6]=false;
  quadrant[7]=false;
  PutHalfCircle(xin,yin,rin,xin,yin+rin,xin,yin-rin);
};

GeomPlane GeomDividedCircle::CrossPlane(GeomPlane &pl)
{
  GeomPlane ret;

  real dist=pl.CalcDistance(center);
  if(dist>outer_circle_r)return ret;

  real xx=center.getx();
  real yy=center.gety();
  real rr=r;
  real r2=rr*rr;

  int ip=0;
  GeomVector pp[2];

  // The following is from GeomPolygon::GetPlane
  for(int i=0;i<nline;i++){
    if(pl.BeCrossing(line[i])){
      if(ip==0){
        pp[ip]=pl.CrossPlane(line[i]);
        ip++;
      }
      else{
	bool judge=false;
	GeomVector tmp=pl.CrossPlane(line[i]);
        for(int j=0;j<ip;j++){
	  if(pp[j].BeIdentical(tmp))judge=true;
	};
	if(!judge){
          if(ip==2){
            cout<<"Error in GeomDividedCircle::CrossPlane.\n";
	    exit(0);
	  };
          pp[ip]=tmp;
          ip++;
	};
      };
    };
  };

  if(ip==2){
    ret.PutVector(pp[0],pp[1]);
    return ret;
  };

  // The following is from GeomCircle::GetPlane
  real x1,x2,y1,y2,dd;
  real a=pl.get_a();
  real b=pl.get_b();
  real c=pl.get_c();
  if(a!=0.0){
    real inv_a=1./a;
    real inv_a2=inv_a*inv_a;
    real alp=c*inv_a+xx;
    real aa=b*b*inv_a2+1.0;
    real bb=2*b*alp*inv_a-2*yy;
    real cc=yy*yy+alp*alp-r2;
    real abc=bb*bb-4*aa*cc;
    if(abc<=0){dd=0;}
    else{dd=sqrt(abc);};
    y1=(-bb+dd)*0.5/aa;
    y2=(-bb-dd)*0.5/aa;
    x1=-b*inv_a*y1-c*inv_a;
    x2=-b*inv_a*y2-c*inv_a;
  }
  else{
    real ee=sqrt(r2-(c/b+yy)*(c/b+yy));
    x1=xx+ee;
    x2=xx-ee;
    real tmp=1./b;
    y1=-c*tmp;
    y2=y1;
  };

  GeomVector ptmp[2];
  ptmp[0].put_data(x1,y1);
  ptmp[1].put_data(x2,y2);
  bool ident=ptmp[0].BeIdentical(ptmp[1]);

  bool point1=false;
  bool point2=false;
  real relx1=x1-xx;
  real rely1=y1-yy;
  real relx2=x2-xx;
  real rely2=y2-yy;
  if(quadrant[0]){
    if(relx1>=0.&&rely1>=0.&&relx1>=rely1)point1=true;
    if(relx2>=0.&&rely2>=0.&&relx2>=rely2)point2=true;
  };
  if(quadrant[1]){
    if(relx1>=0.&&rely1>=0.&&relx1<=rely1)point1=true;
    if(relx2>=0.&&rely2>=0.&&relx2<=rely2)point2=true;
  };
  if(quadrant[2]){
    if(relx1<=0.&&rely1>=0.&&-relx1<=rely1)point1=true;
    if(relx2<=0.&&rely2>=0.&&-relx2<=rely2)point2=true;
  };
  if(quadrant[3]){
    if(relx1<=0.&&rely1>=0.&&-relx1>=rely1)point1=true;
    if(relx2<=0.&&rely2>=0.&&-relx2>=rely2)point2=true;
  };
  if(quadrant[4]){
    if(relx1<=0.&&rely1<=0.&&-relx1>=-rely1)point1=true;
    if(relx2<=0.&&rely2<=0.&&-relx2>=-rely2)point2=true;
  };
  if(quadrant[5]){
    if(relx1<=0.&&rely1<=0.&&-relx1<=-rely1)point1=true;
    if(relx2<=0.&&rely2<=0.&&-relx2<=-rely2)point2=true;
  };
  if(quadrant[6]){
    if(relx1>=0.&&rely1<=0.&&relx1<=-rely1)point1=true;
    if(relx2>=0.&&rely2<=0.&&relx2<=-rely2)point2=true;
  };
  if(quadrant[7]){
    if(relx1>=0.&&rely1<=0.&&relx1>=-rely1)point1=true;
    if(relx2>=0.&&rely2<=0.&&relx2>=-rely2)point2=true;
  };

  if(point1&&point2){
    if(ip==0){
      if(!ident){
	ret.PutVector(ptmp[0],ptmp[1]);
      };
      return ret;
    }else{
      if(ident){
	ret.PutVector(pp[0],ptmp[0]);
      };
      return ret;
    };
  }else if(point1&&ip==1){
    if(!ptmp[0].BeIdentical(pp[0])){
      ret.PutVector(pp[0],ptmp[0]);
      return ret;
    };
  }else if(point2&&ip==1){
    if(!ptmp[1].BeIdentical(pp[0])){
      ret.PutVector(pp[0],ptmp[1]);
      return ret;
    };
  };

  return ret;
};

void GeomDividedCircle::PutPlane(GeomPlane pl)
{
  line.push_back(pl);
  nline++;
};

GeomDividedCircle GeomDividedCircle::GetReflectedCopy(real xvec,real yvec,int d_regid)
{
  GeomVector u(-yvec,xvec);
  GeomVector center2=center.OperatingElementaryReflector(u);

  GeomDividedCircle ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);

  ret.PutR(r);
  for(int i=0;i<nline;i++){
    ret.PutPlane(line[i].OperatingElementaryReflector(u));
  };
 
  if(octant)ret.OctantTrue();

  if(xvec==1.&&yvec==1.){
    int quad2[]={1,0,7,6,5,4,3,2}; 
    for(int i=0;i<8;i++){
      if(quadrant[i]){
	ret.QuadrantTrue(quad2[i]);
      }else{
	ret.QuadrantFalse(quad2[i]);
      };
    };
  }else if(xvec==1.&&yvec==0.){
    int quad2[]={7,6,5,4,3,2,1,0};
    for(int i=0;i<8;i++){
      if(quadrant[i]){
	ret.QuadrantTrue(quad2[i]);
      }else{
	ret.QuadrantFalse(quad2[i]);
      };
    };
  }else  if(xvec==0.&&yvec==1.){
    int quad2[]={3,2,1,0,7,6,5,4};
    for(int i=0;i<8;i++){
      if(quadrant[i]){
	ret.QuadrantTrue(quad2[i]);
      }else{
	ret.QuadrantFalse(quad2[i]);
      };
    };
  }else{
    cout<<"# Error in GeomDividedCircle::GetReflectedCopy.\n";
    cout<<"# Not coded for symmetric line.\n";
    exit(0);
  };  

  return ret;
};

GeomDividedCircle GeomDividedCircle::GetRotatedCopy(real rot,int d_regid)
{
  GeomVector center2=center.OperatingRotator(rot);

  GeomDividedCircle ret;
  ret.PutCenter(center2);
  ret.PutRegionID(RegionID+d_regid);
  ret.PutOuterCircleR(outer_circle_r);

  ret.PutR(r);
  for(int i=0;i<nline;i++){
    ret.PutPlane(line[i].OperatingRotator(rot));
  };
 
  if(octant)ret.OctantTrue();

  int delta=0;
  if(rot>0.1245&&rot<0.1255){
    delta=1;
  }else if(rot>0.249&&rot<0.251){    
    delta=2;
  }else if(rot>0.499&rot<0.501){
    delta=4;
  }else if(rot>0.6245&rot<0.6255){
    delta=5;
  }else{
    cout<<"# Error in GeomDividedCircle::GetRotatedCopy.\n";
    exit(0);
  };

  for(int i=0;i<8;i++){

    int ii=i+delta;
    if(ii>7)ii-=8;
    if(quadrant[i]){
      ret.QuadrantTrue(ii);
    }else{
      ret.QuadrantFalse(ii);
    };

  };
  return ret;
};

void GeomDividedCircle::Shift(real dx, real dy)
{
  center.Shift(dx,dy);
  for(int i=0;i<nline;i++){
    line[i].Shift(dx,dy);
  };
};

// +++++++++++++++++++++++++++++++++++++++++++++++++

void IrregularGeometryInformation::PutLWRCell
(real x,real y,real pitch,real rpel,real rcld,int dpel,int dmod,int &regid,bool triangular)
{
  real hpitch=pitch*0.5;

  // (background triangular for moderator)
  int quadaa[]={3,1,2,4,2,4,3,1};
  bool tri_on[]={false,false,false,true,false,true,true,true};
  real xtmp=x-hpitch*0.5;
  real ytmp=y+hpitch*0.5;
  for(int i=0;i<8;i++){
    if(!triangular||(triangular&&tri_on[i])){
      GeomPolygon pol;
      pol.PutTriangle(xtmp,ytmp,hpitch,quadaa[i]);
      pol.PutRegionID(regid++);
      AddGeom(pol);
    };
    if(i==1||i==5)xtmp+=hpitch;
    if(i==3){
      xtmp-=hpitch;
      ytmp-=hpitch;
    };
  };

  // (moderator)
  real vol_pin=rcld*rcld*PI;
  real vol_pel=rpel*rpel*PI;
  real vol_mod=pitch*pitch-vol_pin;
  vector<real> ring(dpel+1+dmod-1);
  for(int i=0;i<dmod-1;i++){
    real vol_tmp=vol_pin+vol_mod/dmod*(dmod-i-1);
    real r_tmp=sqrt(vol_tmp/PI);
    if(r_tmp>hpitch)r_tmp=hpitch;
    ring[i]=r_tmp;
  };
  ring[dmod-1]=rcld;
  ring[dmod]=rpel;
  for(int i=0;i<dpel-1;i++){
    real vol_tmp=vol_pel/dpel*(dpel-i-1);
    ring[dmod+i+1]=sqrt(vol_tmp/PI);
  };
    
  int quaddd[]={4,3,2,1,5,6,7,8};
  for(int i=0;i<dpel+dmod;i++){
    for(int j=0;j<8;j++){
      if(!triangular||(triangular&&tri_on[j])){
        GeomDividedCircle cir;
        cir.PutHalfQuarterCircle(x,y,ring[i],quaddd[j]);
        cir.PutRegionID(regid++);
        AddGeom(cir);
      };
    };
  };

};

void IrregularGeometryInformation::PutLWRCellNoDivision
(real x,real y,real pitch,real rpel,real rcld,int &regid)
{
  real hpitch=pitch*0.5;

  // (background triangular for moderator)
  GeomPolygon pol;
  pol.PutRectangular(x,y,hpitch,hpitch);
  pol.PutRegionID(regid++);
  AddGeom(pol);

  for(int i=0;i<2;i++){
    real r=rcld;
    if(i==1)r=rpel;
    GeomCircle cir(x,y,r);
    cir.PutRegionID(regid++);
    AddGeom(cir);
  };

};

void IrregularGeometryInformation::PutLWRCellNoDivisionHalfQuarter
(real x,real y,real pitch,real rpel,real rcld,int &regid)
{
  // (background triangular for moderator)
  GeomPolygon pol;
  pol.PutTriangle(x,y,pitch,4);
  pol.PutRegionID(regid++);
  AddGeom(pol);

  for(int i=0;i<2;i++){
    real r=rcld;
    if(i==1)r=rpel;
    GeomDividedCircle cir;
    cir.PutR(r);
    cir.PutCenter(x,y);
    cir.PutRegionID(regid++);
    cir.QuadrantTrue(0);
    cir.QuadrantFalse(1);
    cir.QuadrantFalse(2);
    cir.QuadrantFalse(3);
    cir.QuadrantFalse(4);
    cir.QuadrantTrue(5);
    cir.QuadrantTrue(6);
    cir.QuadrantTrue(7);
    cir.PutNline(1);
    GeomPlane pl;
    real par=r/sqrt(2.);
    pl.PutVector(x+par,y+par,x-par,y-par);
    cir.PutPlane(pl);
    AddGeom(cir);
  };

};

void IrregularGeometryInformation::PutSquareDividedCell(real x,real y,real pitch,int &regid,bool triangular)
{
  real hpitch=pitch*0.5;

  // (background triangular for moderator)
  int quadaa[]={3,1,2,4,2,4,3,1};
  bool tri_on[]={false,false,false,true,false,true,true,true};
  real xtmp=x-hpitch*0.5;
  real ytmp=y+hpitch*0.5;
  for(int i=0;i<8;i++){
    if(!triangular||(triangular&&tri_on[i])){
      GeomPolygon pol;
      pol.PutTriangle(xtmp,ytmp,hpitch,quadaa[i]);
      pol.PutRegionID(regid++);
      AddGeom(pol);
    };
    if(i==1||i==5)xtmp+=hpitch;
    if(i==3){
      xtmp-=hpitch;
      ytmp-=hpitch;
    };
  };
};

void IrregularGeometryInformation::PutSquareCell(real x,real y,real pitch,int &regid)
{
  real hpitch=pitch*0.5;
  GeomPolygon pol;
  pol.PutRectangular(x,y,hpitch,hpitch);
  pol.PutRegionID(regid++);
  AddGeom(pol);
};

void IrregularGeometryInformation::PutSquareCellHalfQuarter(real x,real y,real pitch,int &regid)
{
  GeomPolygon pol;
  pol.PutTriangle(x,y,pitch,4);
  pol.PutRegionID(regid++);
  AddGeom(pol);
};

void IrregularGeometryInformation::PutTriangleAssembly(real x,real y,real len,int layer,int &regid)
{
  int quadorg[]={4,3,1,2};
  real unit_size=len*0.5;

  int num=0;
  real yu=1.;
  for(int i=layer-1;i>=0;i--){
    int part=i*2+1;
    real xu=1.+(layer-1-i)*2.;
    for(int j=0;j<part;j++){
      GeomPolygon pol;
      pol.PutTriangle(x+xu*unit_size,y+yu*unit_size,len,quadorg[j%4]);
      pol.PutRegionID(regid++);
      AddGeom(pol);
      if(j%2==0)xu+=2.;      
    };
    yu+=2.;
  };
};

void IrregularGeometryInformation::PutTriangleAssemblySquare(real x,real y,real lenx,real leny,int layerx,int layery,int &regid,bool quadadd)
{
  int quadorg[]={2,4,3,1};
  real unit_size_x=lenx*0.5;
  real unit_size_y=leny*0.5;

  int num=0;
  real yu=1.;
  for(int i=0;i<layery;i++){
    real xu=1.;
    for(int j=0;j<layerx*2;j++){
      GeomPolygon pol;
      int quadpos=j;
      if(quadadd)quadpos+=2;
      if(i%2==1)quadpos+=2;
      int quad=quadorg[quadpos%4];
      pol.PutTriangle(x+xu*unit_size_x,y+yu*unit_size_y,lenx,leny,quad);
      pol.PutRegionID(regid++);
      AddGeom(pol);
      if(j%2==1)xu+=2.;      
    };
    yu+=2.;
  };
};

void IrregularGeometryInformation::PutTriangleAssemblyTriangle(real x,real y,real lenx,real leny,int layerx,int layery,int &regid,bool quadadd)
{
  int quadorg[]={2,4,3,1};
  real unit_size_x=lenx*0.5;
  real unit_size_y=leny*0.5;

  int num=0;
  real yu=1.;
  for(int i=0;i<layery;i++){

    int quadpos=1;    
    real xu=2*i+1;
    {
    GeomPolygon pol;
    int quad=quadorg[quadpos%4];
    pol.PutTriangle(x+xu*unit_size_x,y+yu*unit_size_y,lenx,leny,quad);
    pol.PutRegionID(regid++);
    AddGeom(pol);
    };

    xu+=2;
    int boxnum=layerx-2-i*2;
    for(int j=0;j<boxnum*2;j++){
      GeomPolygon pol;
      quadpos++;
      int quad=quadorg[quadpos%4];
      pol.PutTriangle(x+xu*unit_size_x,y+yu*unit_size_y,lenx,leny,quad);
      pol.PutRegionID(regid++);
      AddGeom(pol);
      if(j%2==1)xu+=2.;      
    };

    {
      quadpos++;
      GeomPolygon pol;
      int quad=quadorg[quadpos%4];
      pol.PutTriangle(x+xu*unit_size_x,y+yu*unit_size_y,lenx,leny,quad);
      pol.PutRegionID(regid++);
      AddGeom(pol);
    };
    
    yu+=2.;
  };
};

