#include<iostream>
#include<cstdlib>
#include "CartMeshInfo.h"


using namespace std;


void CartMeshInfo::PutMeshXY(int *fmx,int *fmy,real z,CartCore &core)
{
  int xr=core.GetXr();
  int yr=core.GetYr();

  int *mat=new int[xr*yr];
  core.GetMaterialMapXY(z,mat);

  real *xwid = new real[xr];
  real *ywid = new real[yr];
  real *zwid = new real[1];
  for(int x=0;x<xr;x++){
    xwid[x]=core.GetXwid(x);
  };
  for(int y=0;y<yr;y++){
    ywid[y]=core.GetYwid(y);
  };
  PutMeshInfo(xr,yr,fmx,fmy,xwid,ywid,mat);

  int bc[6];
  bc[0]=core.GetLeftBC();
  bc[1]=core.GetRightBC();
  bc[2]=core.GetBackBC();
  bc[3]=core.GetFrontBC();
  bc[4]=0;
  bc[5]=0;
  PutBoundaryCondition(bc);

  delete [] mat;
  delete [] xwid;
  delete [] ywid;
  delete [] zwid;
}

void CartMeshInfo::PutMeshXYZ(int *fmx,int *fmy,real zd,CartCore &core,real zwidc)
{
  int xr=core.GetXr();
  int yr=core.GetYr();
  int zr=core.GetUnifiedZMesh();

  int *mat=new int[xr*yr*zr];
  real *zwid=new real[zr];

  core.MakeMaterialMap(zr,mat,zwid);

  cout<<"# Z-mesh division : ";
  int *fmz = new int[zr];
  for(int z=0;z<zr;z++){
    fmz[z]=int(ceil(zwid[z]/zd));
    //if(z==7)fmz[z]=22;
    //if(z==8)fmz[z]=36;
    cout<<fmz[z]<<" ";
    if(fmz[z]<0){
      cout<<"\n# Error in CartMeshInfo::PutMeshXYZ.\n";
      cout<<"# Negative mesh is detected along to the Z-mesh.\n";
      exit(0);
    };
  }
  cout<<"\n";

  int zr2=0;
  int fmzdm[999];
  real zwiddm[999];
  int *matdm = new int[xr*yr*999];
  for(int i=0;i<zr;i++){
    real cmz=0.;
    int ind=0;
    for(int j=0;j<fmz[i];j++){
      ind++;
      cmz+=zwid[i]/fmz[i];
      if(cmz>=zwidc||j==fmz[i]-1){
	cmz=0.;
	fmzdm[zr2]=ind;
	zwiddm[zr2]=zwid[i]/fmz[i]*ind;
	for(int k=0;k<xr*yr;k++){matdm[xr*yr*zr2+k]=mat[xr*yr*i+k];};
	ind=0;
	zr2++;
      };
    };
  };
  int *fmz2 = new int[zr2];
  real *zwid2 = new real[zr2];
  int *mat2 = new int[xr*yr*zr2];
  for(int i=0;i<zr2;i++){
    fmz2[i]=fmzdm[i];
    zwid2[i]=zwiddm[i];
    for(int j=0;j<xr*yr;j++){
      mat2[i*xr*yr+j]=matdm[i*xr*yr+j];
    };
  };

  real *xwid = new real[xr];
  real *ywid = new real[yr];
  for(int x=0;x<xr;x++){xwid[x]=core.GetXwid(x);};
  for(int y=0;y<yr;y++){ywid[y]=core.GetYwid(y);};

  PutMeshInfo(xr,yr,zr2,fmx,fmy,fmz2,xwid,ywid,zwid2,mat2);

  int bc[6];
  bc[0]=core.GetLeftBC();
  bc[1]=core.GetRightBC();
  bc[2]=core.GetBackBC();
  bc[3]=core.GetFrontBC();
  bc[4]=core.GetUpperBC();
  bc[5]=core.GetBottomBC();
  PutBoundaryCondition(bc);

  delete [] mat;
  delete [] zwid;
  delete [] fmz;
  delete [] matdm;
  delete [] fmz2;
  delete [] zwid2;
  delete [] mat2;
  delete [] xwid;
  delete [] ywid;
}

void CartMeshInfo::PutMeshXYZ(int *fmx,int *fmy,int *fmz,CartCore &core,real zwidc)
{
  int xr=core.GetXr();
  int yr=core.GetYr();
  int zr=core.GetUnifiedZMesh();

  int *mat=new int[xr*yr*zr];
  real *zwid=new real[zr];

  core.MakeMaterialMap(zr,mat,zwid);

  cout<<"# Z-mesh division : ";
  for(int z=0;z<zr;z++){
    cout<<fmz[z]<<" ";
  }
  cout<<"\n";

  int zr2=0;
  int fmzdm[999];
  real zwiddm[999];
  int *matdm = new int[xr*yr*999];
  for(int i=0;i<zr;i++){
    real cmz=0.;
    int ind=0;
    for(int j=0;j<fmz[i];j++){
      ind++;
      cmz+=zwid[i]/fmz[i];
      if(cmz>=zwidc||j==fmz[i]-1){
	cmz=0.;
	fmzdm[zr2]=ind;
	zwiddm[zr2]=zwid[i]/fmz[i]*ind;
	for(int k=0;k<xr*yr;k++){matdm[xr*yr*zr2+k]=mat[xr*yr*i+k];};
	ind=0;
	zr2++;
      };
    };
  };
  int *fmz2 = new int[zr2];
  real *zwid2 = new real[zr2];
  int *mat2 = new int[xr*yr*zr2];
  for(int i=0;i<zr2;i++){
    fmz2[i]=fmzdm[i];
    zwid2[i]=zwiddm[i];
    for(int j=0;j<xr*yr;j++){
      mat2[i*xr*yr+j]=matdm[i*xr*yr+j];
    };
  };

  real *xwid = new real[xr];
  real *ywid = new real[yr];
  for(int x=0;x<xr;x++){xwid[x]=core.GetXwid(x);};
  for(int y=0;y<yr;y++){ywid[y]=core.GetYwid(y);};

  PutMeshInfo(xr,yr,zr2,fmx,fmy,fmz2,xwid,ywid,zwid2,mat2);

  int bc[6];
  bc[0]=core.GetLeftBC();
  bc[1]=core.GetRightBC();
  bc[2]=core.GetBackBC();
  bc[3]=core.GetFrontBC();
  bc[4]=core.GetUpperBC();
  bc[5]=core.GetBottomBC();
  PutBoundaryCondition(bc);

  delete [] mat;
  delete [] zwid;
  delete [] matdm;
  delete [] fmz2;
  delete [] zwid2;
  delete [] mat2;
  delete [] xwid;
  delete [] ywid;
}

void CartMeshInfo::PutMeshInfo(int xr,int *fmx,real *xl,int *mat,string type)
{
  int fmz[]={1};
  real zl[]={1.};
  PutMeshInfo(xr,1,1,fmx,fmz,fmz,xl,zl,zl,mat,type);
};

void CartMeshInfo::PutMeshInfo(int xr,int yr,int *fmx,int *fmy,
                               real *xl,real *yl,int *mat,string type)
{
  int fmz[]={1};
  real zl[]={1.};
  PutMeshInfo(xr,yr,1,fmx,fmy,fmz,xl,yl,zl,mat,type);
};

void CartMeshInfo::PutMeshInfo(int xr,int yr,int zr,int *fmx,int *fmy,int *fmz,
                               real *xl,real *yl,real *zl,int *mat,string type)
{
  if(type!="cumulative"&&type!="width"&&
     type!="Cumulative"&&type!="Width"){
    cout<<"# Error in PutMeshInfo for CartMeshInfo.\n";
    cout<<"# Invalid type : "<<type<<"\n";
    exit(0);
  };

  CMesh[0]=xr;
  CMesh[1]=yr;
  CMesh[2]=zr;

  CMeshF.resize(3);
  CMeshL.resize(3);
  FMeshL.resize(3);

  vector<real> xll(xr);
  vector<real> yll(yr);
  vector<real> zll(zr);

  for(int i=0;i<3;i++){
    CMeshF[i].resize(CMesh[i]);
    CMeshL[i].resize(CMesh[i]);
  };

  for(int i=0;i<xr;i++){xll[i]=xl[i];};
  for(int i=0;i<yr;i++){yll[i]=yl[i];};
  for(int i=0;i<zr;i++){zll[i]=zl[i];};


  if(type=="cumulative"||type=="Cumulative"){
    for(int i=CMesh[0]-1;i>0;i--){xll[i]=xl[i]-xl[i-1];};
    for(int i=CMesh[1]-1;i>0;i--){yll[i]=yl[i]-yl[i-1];};
    for(int i=CMesh[2]-1;i>0;i--){zll[i]=zl[i]-zl[i-1];};
  };


  bool zero_xl=false;
  for(int i=0;i<CMesh[0];i++){
    if(fmx[i]<=0){
      cout<<"# Error in PutMeshInfo.\n";
      cout<<"# The number of x-meshes is less than 0.\n";
      exit(0);
    };
    if(xl[i]<=0.){
      zero_xl=true;
    };
    CMeshF[0][i]=fmx[i];
    CMeshL[0][i]=xll[i];
  };
  if(zero_xl){
    //cout<<"++ Warning in CartMeshInfo::PutMeshInfo.\n ++";
    //cout<<"   The length of x-mesh is less than 0.\n";
  };

  bool zero_yl=false;
  for(int i=0;i<CMesh[1];i++){
    if(fmy[i]<=0){
      cout<<"# Error in PutMeshInfo.\n";
      cout<<"# The number of y-meshes is less than 0.\n";
      exit(0);
    };
    if(yl[i]<=0.){
      zero_yl=true;
    };
    CMeshF[1][i]=fmy[i];
    CMeshL[1][i]=yll[i];
  };
  if(zero_yl){
    //cout<<"++ Warning in CartMeshInfo::PutMeshInfo.++\n";
    //cout<<"   The length of y-mesh is less than 0.\n";
  };


  for(int i=0;i<CMesh[2];i++){
    if(fmz[i]<=0){
      cout<<"# Error in PutMeshInfo.\n";
      cout<<"# The number of z-meshes is less than 0.\n";
      exit(0);
    };
    if(zl[i]<=0.){
      cout<<"# Error in PutMeshInfo.\n";
      cout<<"# The length of z-mesh is less than 0.\n";
      exit(0);
    };
    CMeshF[2][i]=fmz[i];
    CMeshL[2][i]=zll[i];
  };

  for(int i=0;i<3;i++){
    int k=0;
    for(int j=0;j<CMesh[i];j++){
      k+=CMeshF[i][j];
    };
    FMesh[i]=k;
  };

  for(int i=0;i<3;i++){
    FMeshL[i].resize(FMesh[i]);
    int l=0;
    for(int j=0;j<CMesh[i];j++){
      for(int k=0;k<CMeshF[i][j];k++){
	FMeshL[i][l]=CMeshL[i][j]/CMeshF[i][j];
	l++;
      };
    };
  };

  CMat.resize(CMesh[0]*CMesh[1]*CMesh[2]);
  for(int i=0;i<zr*yr*xr;i++){
    CMat[i]=mat[i];
  };

  FMat.resize(FMesh[0]*FMesh[1]*FMesh[2]);
  //FMat_par_mesh.resize(FMesh[0]*FMesh[1]*FMesh[2]);
  FMat_par_mesh.clear();

  int index=0;
  for(int z=0;z<zr;z++){
    for(int z2=0;z2<CMeshF[2][z];z2++){
      for(int y=0;y<yr;y++){
        for(int y2=0;y2<CMeshF[1][y];y2++){
          for(int x=0;x<xr;x++){
	    for(int x2=0;x2<CMeshF[0][x];x2++){
	      int tmp=mat[z*(xr*yr)+y*xr+x];
	      FMat[index]=tmp;
	      index++;
	      if(tmp!=-1){
  	        FMat_par_mesh.push_back(tmp);
	      };
	    };
          };
	};
      };
    };
  };

}


void CartMeshInfo::PutBoundaryCondition(string xl,string xr,string yl,string yr,string zl,string zr)
{
  string bcname[]={"Zeroflux","Reflective","Vacuum","Periodic"};
  for(int i=0;i<4;i++){
    if(xl==bcname[i])BC[0]=i;
    if(xr==bcname[i])BC[1]=i;
    if(yl==bcname[i])BC[2]=i;
    if(yr==bcname[i])BC[3]=i;
    if(zl==bcname[i])BC[4]=i;
    if(zr==bcname[i])BC[5]=i;
  };
  for(int i=0;i<6;i++){
    if(BC[i]<0||BC[i]>3){
      cout<<"Incorrect boundary condition in CartMeshInfo.\n";
      exit(0);
    };
  };
};

void CartMeshInfo::PutBoundaryCondition(vector<int> binp)
{
  int inp[6];
  for(int i=0;i<6;i++){
    inp[i]=binp[i];
  };
  PutBoundaryCondition(inp);
};

void CartMeshInfo::PutBoundaryCondition(int *binp)
{
  string bcname[]={"Zeroflux","Reflective","Vacuum","Periodic"};
  for(int i=0;i<6;i++){
    if(binp[i]<0||binp[i]>3){
      cout<<"Error in PutBoundaryCondition in CartMeshInfo.\n";
      exit(0);
    };
  };
  PutBoundaryCondition(bcname[binp[0]],bcname[binp[1]],
                       bcname[binp[2]],bcname[binp[3]],
                       bcname[binp[4]],bcname[binp[5]]);
};

// Converter for vector input into pointer input

void CartMeshInfo::PutMeshXY(vector<int> fmx,vector<int> fmy,real z,CartCore &core)
{
  int xr=core.GetXr();
  int yr=core.GetYr();
  int *fmxi=new int[xr];
  int *fmyi=new int[yr];
  for(int i=0;i<xr;i++){
    fmxi[i]=fmx[i];
    fmyi[i]=fmy[i];
  };
  PutMeshXY(fmxi,fmyi,z,core);
  delete [] fmxi;
  delete [] fmyi;
};

void CartMeshInfo::PutMeshXYZ(vector<int> fmx,vector<int> fmy,real zd,CartCore &core,real zwidc)
{
  int xr=core.GetXr();
  int yr=core.GetYr();
  int *fmxi=new int[xr];
  int *fmyi=new int[yr];
  for(int i=0;i<xr;i++){
    fmxi[i]=fmx[i];
    fmyi[i]=fmy[i];
  };
  PutMeshXYZ(fmxi,fmyi,zd,core,zwidc);
  delete [] fmxi;
  delete [] fmyi;
}
void CartMeshInfo::PutMeshInfo(int xr,int yr,int zr,vector<int> fmx,vector<int> fmy,vector<int> fmz,
                               vector<real> xl,vector<real> yl,vector<real> zl,vector<int> mat,string type)
{
  int *tmp1=new int[xr];
  int *tmp2=new int[yr];
  int *tmp3=new int[zr];
  real *tmp4=new real[xr];
  real *tmp5=new real[yr];
  real *tmp6=new real[zr];
  int *tmp7=new int[xr*yr*zr];
  for(int i=0;i<xr;i++){
    tmp1[i]=fmx[i];
    tmp4[i]=xl[i];
  };
  for(int i=0;i<yr;i++){
    tmp2[i]=fmy[i];
    tmp5[i]=yl[i];
  };
  for(int i=0;i<zr;i++){
    tmp3[i]=fmz[i];
    tmp6[i]=zl[i];
  };
  for(int i=0;i<xr*yr*zr;i++){
    tmp7[i]=mat[i];
  };
  PutMeshInfo(xr,yr,zr,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,type);
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3;
  delete [] tmp4;
  delete [] tmp5;
  delete [] tmp6;
  delete [] tmp7;
};

void CartMeshInfo::PutMeshInfo(int xr,int yr,vector<int> fmx,vector<int> fmy,
                               vector<real> xl,vector<real> yl,vector<int> mat,string type)
{
  int *tmp1=new int[xr];
  int *tmp2=new int[yr];
  real *tmp4=new real[xr];
  real *tmp5=new real[yr];
  int *tmp7=new int[xr*yr];
  for(int i=0;i<xr;i++){
    tmp1[i]=fmx[i];
    tmp4[i]=xl[i];
  };
  for(int i=0;i<yr;i++){
    tmp2[i]=fmy[i];
    tmp5[i]=yl[i];
  };
  for(int i=0;i<xr*yr;i++){
    tmp7[i]=mat[i];
  };
  PutMeshInfo(xr,yr,tmp1,tmp2,tmp4,tmp5,tmp7,type);
  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp4;
  delete [] tmp5;
  delete [] tmp7;
};

void CartMeshInfo::PutMeshInfo(int xr,vector<int> fmx,
                               vector<real> xl,vector<int> mat,string type)
{
  int *tmp1=new int[xr];
  real *tmp4=new real[xr];
  int *tmp7=new int[xr];
  for(int i=0;i<xr;i++){
    tmp1[i]=fmx[i];
    tmp4[i]=xl[i];
  };
  for(int i=0;i<xr;i++){
    tmp7[i]=mat[i];
  };
  PutMeshInfo(xr,tmp1,tmp4,tmp7,type);
  delete [] tmp1;
  delete [] tmp4;
  delete [] tmp7;
};

void CartMeshInfo::ReconstructCoarseMesh(real maxxl,real maxyl,real maxzl)
{
  int xf=FMesh[0];
  int yf=FMesh[1];
  int zf=FMesh[2];

  // x-direction
  if(maxxl>0.){
  vector<bool> boundary_x(xf-1);
  for(int y=0;y<yf;y++){
    for(int z=0;z<zf;z++){
      for(int x=0;x<xf-1;x++){
	int mat_left=FMat[z*xf*yf+y*xf+x];
	int mat_right=FMat[z*xf*yf+y*xf+x+1];
	if((mat_left==-1&&mat_right!=-1)||
           (mat_left!=-1&&mat_right==-1)){
          boundary_x[x]=true;
	};
      };
    };
  };
  real x=0.;
  for(int i=0;i<xf-1;i++){
    if(boundary_x[i]){
      x=0.;
    }else{
      x+=FMeshL[0][i];
      if(x+FMeshL[0][i+1]>maxxl){
	boundary_x[i]=true;
	x=0.;
      };
    };
  };
  int num=0;
  for(int i=0;i<xf-1;i++){
    if(boundary_x[i]){
      num++;
    };
  };
  num++;
  CMesh[0]=num;
  CMeshF[0].resize(num,0);
  CMeshL[0].resize(num,0.);
  x=0.;
  num=0;
  int cc=0;
  for(int i=0;i<xf-1;i++){
    x+=FMeshL[0][i];
    num++;
    if(boundary_x[i]){
      CMeshF[0][cc]=num;
      CMeshL[0][cc]=x;
      num=0;
      cc++;
      x=0.;
    };
  };
  x+=FMeshL[0][xf-1];
  num++;
  CMeshF[0][cc]=num;
  CMeshL[0][cc]=x;
  };

  // y-direction
  if(maxyl>0.){
  vector<bool> boundary_y(yf-1);
  for(int x=0;x<xf;x++){
    for(int z=0;z<zf;z++){
      for(int y=0;y<yf-1;y++){
	int mat_left=FMat[z*xf*yf+y*xf+x];
	int mat_right=FMat[z*xf*yf+(y+1)*xf+x];
	if((mat_left==-1&&mat_right!=-1)||
           (mat_left!=-1&&mat_right==-1)){
          boundary_y[y]=true;
	};
      };
    };
  };
  real x=0.;
  for(int i=0;i<yf-1;i++){
    if(boundary_y[i]){
      x=0.;
    }else{
      x+=FMeshL[1][i];
      if(x+FMeshL[1][i+1]>maxyl){
	boundary_y[i]=true;
	x=0.;
      };
    };
  };
  int num=0;
  for(int i=0;i<yf-1;i++){
    if(boundary_y[i]){
      num++;
    };
  };
  num++;
  CMesh[1]=num;
  CMeshF[1].resize(num,0);
  CMeshL[1].resize(num,0.);
  x=0.;
  num=0;
  int cc=0;
  for(int i=0;i<yf-1;i++){
    x+=FMeshL[1][i];
    num++;
    if(boundary_y[i]){
      CMeshF[1][cc]=num;
      CMeshL[1][cc]=x;
      num=0;
      cc++;
      x=0.;
    };
  };
  x+=FMeshL[1][yf-1];
  num++;
  CMeshF[1][cc]=num;
  CMeshL[1][cc]=x;
  };

  // z-direction
  if(maxzl>0.){
  vector<bool> boundary_z(zf-1);
  for(int x=0;x<xf;x++){
    for(int y=0;y<yf;y++){
      for(int z=0;z<zf-1;z++){
	int mat_left=FMat[z*xf*yf+y*xf+x];
	int mat_right=FMat[(z+1)*xf*yf+y*xf+x];
	if((mat_left==-1&&mat_right!=-1)||
           (mat_left!=-1&&mat_right==-1)){
          boundary_z[z]=true;
	};
      };
    };
  };
  real x=0.;
  for(int i=0;i<zf-1;i++){
    if(boundary_z[i]){
      x=0.;
    }else{
      x+=FMeshL[2][i];
      if(x+FMeshL[2][i+1]>maxzl){
	boundary_z[i]=true;
	x=0.;
      };
    };
  };
  int num=0;
  for(int i=0;i<zf-1;i++){
    if(boundary_z[i]){
      num++;
    };
  };
  num++;
  CMesh[2]=num;
  CMeshF[2].resize(num,0);
  CMeshL[2].resize(num,0.);
  x=0.;
  num=0;
  real cc=0;
  for(int i=0;i<zf-1;i++){
    x+=FMeshL[2][i];
    num++;
    if(boundary_z[i]){
      CMeshF[2][cc]=num;
      CMeshL[2][cc]=x;
      num=0;
      cc++;
      x=0.;
    };
  };
  x+=FMeshL[2][zf-1];
  num++;
  CMeshF[2][cc]=num;
  CMeshL[2][cc]=x;
  };

  CMat.resize(CMesh[0]*CMesh[1]*CMesh[2],0);

  int xr=CMesh[0];
  int yr=CMesh[1];
  int zr=CMesh[2];

  for(int i=0;i<xr*yr*zr;i++){CMat[i]=-1;};

  int index=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<CMeshF[2][z1];z2++){
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<CMeshF[1][y1];y2++){
          for(int x1=0;x1<xr;x1++){
	    for(int x2=0;x2<CMeshF[0][x1];x2++){
	      int tmp=z1*(xr*yr)+y1*xr+x1;
	      if(FMat[index]!=-1)CMat[tmp]=1;
	      index++;
	    };
	  };
	};
      };
    };
  };

  if(!CheckCoarseMeshMixedVacuum()){
    cout<<"Error in ReconstructCoarseMesh.\n";
    cout<<"Vacuum and real regions are mixed.\n";
    exit(0);
  };
};

bool CartMeshInfo::CheckCoarseMeshMixedVacuum()
{
  int xr=CMesh[0];
  int yr=CMesh[1];
  int zr=CMesh[2];

  int index=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<CMeshF[2][z1];z2++){
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<CMeshF[1][y1];y2++){
          for(int x1=0;x1<xr;x1++){
	    for(int x2=0;x2<CMeshF[0][x1];x2++){
	      int tmp=z1*(xr*yr)+y1*xr+x1;
	      if((FMat[index]==-1&&CMat[tmp]==1)||
                 (FMat[index]!=-1&&CMat[tmp]==-1)){
		return false;
	      };
	      index++;
	    };
	  };
	};
      };
    };
  };
  return true;
};

void CartMeshInfo::show_self()
{
  int dim=3;
  if(CMesh[2]==1)dim=2;

  if(dim==2){
    cout<<"#\n";
    cout<<"# Material map in coarse mesh division\n";
    cout<<"\n";
    cout<<"     ";
    for(int j=0;j<CMesh[0];j++){
      if(j<10){
	cout<<"  ";
      }else{
	cout<<" ";
      };
      cout<<j;
    };
    cout<<"\n";
    cout<<"#-----";
    for(int j=0;j<CMesh[0];j++){
      cout<<"---";
    };
    cout<<"\n";
    int index=0;
    for(int i=0;i<CMesh[1];i++){
      if(i<10){
	cout<<"  ";
      }else{
	cout<<" ";
      };
      cout<<i<<" |";
      for(int j=0;j<CMesh[0];j++){
	if(CMat[index]<10){
	  cout<<"  ";
	}else{
	  cout<<" ";
	};
        cout<<CMat[index];
	index++;
      };
      cout<<"\n";
    };
    cout<<"#-----";
    for(int j=0;j<CMesh[0];j++){
      cout<<"---";
    };
    cout<<"\n";
  };

  cout<<"#\n";
  cout<<"# Mesh division information\n";
  cout<<"#\n";
  
  for(int i=0;i<dim;i++){
    if(i==0){
      if(dim==3){
	cout<<"# X direction\n#\n";
      }else{
	cout<<"# X or R direction\n#\n";	
      };
    };
    if(i==1){
      if(dim==3){
        cout<<"# Y direction\n#\n";
      }else{
        cout<<"# X or Z direction\n#\n";	
      };
    };
    if(i==2)cout<<"# Z direction\n#\n";
    cout<<"#    Total number of fine meshes   : "<<FMesh[i]<<"\n";
    cout<<"#    Total number of corase meshes : "<<CMesh[i]<<"\n#\n";
    real tot=0.;
    cout<<"#  Coarse    Width      Number of     Mesh outer-edge \n";
    cout<<"#  mesh ID              fine meshes   position        \n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    for(int j=0;j<CMesh[i];j++){
      tot+=CMeshL[i][j];
      cout<<"#    ";
      WriteOut(j,4);
      cout<<"  "<<CMeshL[i][j]<<"   ";
      WriteOut(CMeshF[i][j],5);
      cout<<"         "<<tot<<"\n";
    };
    cout<<"\n";
  };
};

int CartMeshInfo::GetMeshPosition(int dir,real x)
{
  real xpos=0.;
  for(int i=0;i<FMesh[dir];i++){
    xpos+=FMeshL[dir][i];
    if(xpos>=x)return i;
  };
  cout<<"Error in CartMeshInfo::GetMeshPosition.\n";
  exit(0);
};

real CartMeshInfo::GetMeshLocation(int dir,int x)
{
  real ret=0.;
  for(int i=0;i<x;i++){
    ret+=FMeshL[dir][i];
  };
  ret+=FMeshL[dir][x]*0.5;
  return ret;
};

void CartMeshInfo::GetPositionFromMeshID(int id)
{
  int index=0;
  int index2=0;
  for(int z=0;z<FMesh[2];z++){
    for(int y=0;y<FMesh[1];y++){
      for(int x=0;x<FMesh[0];x++){
	if(FMat[index2]!=-1){
	  if(index==id){
	    cout<<" (The position of meshid "<<id<<")\n";
	    cout<<"    X : "<<x<<"\n";
	    cout<<"    Y : "<<y<<"\n";
	    cout<<"    Z : "<<z<<"\n";
	    return;
	  };
          index++;
	};
	index2++;
      };
    };
  };
};
