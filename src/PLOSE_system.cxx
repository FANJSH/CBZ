#include<cstdlib>
#include<iostream>
#include "PLOSE_system.h"

using namespace std;

void PLOSESystem::Init(int n,int g,int i)
{
  GeneralSystem::Init(n,g,i);
  Cylinder = false; 
  Sphere   = false; 
  difc[0]  = 0;
  difc[1]  = 0;
  difc[2]  = 0;
  name="PLOSE";
  PutPL(0);
  cmfdimp = false;  // CMFD is not implemented
};

void PLOSESystem::PutCartMeshInfo(CartMeshInfo cm, string geom)
{
  bool cart=false; 
  if(geom=="Cartesian")        cart    =true;
  if(geom=="Sphere"&&Ndim==1)  Sphere  =true;
  if(geom=="Cylinder"&&Ndim==2)Cylinder=true;

  if(!cart&&!Sphere&&!Cylinder){
    cout<<"Error in PLOSESystem::PutCartMeshInfo.\n";
    cout<<"...Cartesian? Sphere? Cylinder?\n";
    exit(0);
  };

  mi=cm;
  if(Ndim==2&&mi.GetZF()!=1){
    cout<<"This is not two-dimensional mesh.\n";
    cout<<"You must Ndim in system to be 3.\n";
    exit(0);
  };
  if(Ndim==1&&(mi.GetYF()!=1||mi.GetZF()!=1)){
    cout<<"This is not one-dimensional mesh.\n";
    cout<<"You must Ndim in system to be 2 or 3.\n";
    exit(0);
  };

  int xf=mi.GetXF();
  int yf=mi.GetYF();
  int zf=mi.GetZF();
  meshid.resize(zf);
  for(int i=0;i<zf;i++){
    meshid[i].resize(yf);
    for(int j=0;j<yf;j++){
      meshid[i][j].resize(xf);
      for(int k=0;k<xf;k++){
	meshid[i][j][k]=-1;
      };
    };
  };

  int tmp=0;
  int tmp2=0;
  for(int z=0;z<zf;z++){
    for(int y=0;y<yf;y++){
      for(int x=0;x<xf;x++){
	int tm=mi.GetFMat(tmp2);
	if(tm>=nmed){
	  cout<<"Error in PLOSESystem::PutCartMeshInfo.\n";
	  cout<<"You requested not-existing medium ID.\n";
	  cout<<"  (NMED : "<<nmed<<" and TM : "<<tm<<")\n";
          cout<<"Please check Medium ID.\n";
	  exit(0);
	};
	tmp2++;
	if(tm<-1||tm>=nmed){
	  cout<<"Invalid material number for CartMeshInfo (Number:"<<tm<<")\n";
	  exit(0);
	};
	if(tm!=-1)tmp++;
      };
    };
  };
  TotM=tmp;
  if(print)cout<<"#*** Total Mesh : "<<TotM<<"\n";
  mesh.resize(TotM);

  real *di=new real[Ndim];
  int index=0;
  int ind2=0;
  for(int z=0;z<zf;z++){
    for(int y=0;y<yf;y++){
      real tmp=0.;
      for(int x=0;x<xf;x++){
	real tmp2=tmp;
        tmp+=mi.GetFMeshL(0,x);
	int tm=mi.GetFMat(ind2);
	ind2++;
	if(tm!=-1){
	  if(Sphere){
	    mesh[index].PutDimSphere(mi.GetFMeshL(0,x),tmp2,tmp);
  	  }
          else if(Cylinder){
	    mesh[index].PutDimCylinder(mi.GetFMeshL(0,x),tmp2,tmp,mi.GetFMeshL(1,y));
  	  }
	  else{	  // Cartesian
  	    di[0]=mi.GetFMeshL(0,x);
            if(Ndim>1)di[1]=mi.GetFMeshL(1,y);
            if(Ndim>2)di[2]=mi.GetFMeshL(2,z);
            mesh[index].PutDim(Ndim,di);
	  };
          mesh[index].PutPL(0);
          mesh[index].PutGrp(grp);
          mesh[index].PutMedium(&med[tm]);
	  meshid[z][y][x]=index;
	  index++;
	}else{
          meshid[z][y][x]=-1;
        };
      };
    };
  };
  delete [] di;

  EdgeCalculation();

  for(int i=0;i<Ndim*2;i++){
    if(mi.GetBC(i)==2)BC[i]=2; // Vacuum
    if(mi.GetBC(i)==1)BC[i]=1; // Reflective
    if(mi.GetBC(i)==0)BC[i]=0; // Zero flux
    if(BC[i]<0||BC[i]>2){
      cout<<"Boundary Condition Error in PLOSE.\n";
      exit(0);
    };
  };
  if(Ndim!=3){BC[4]=0; BC[5]=0;};
  if(Ndim==1){BC[2]=0; BC[3]=0;};

  // ... plose

  // Edge for corner mesh

  edge_e.resize(3);

  vector<int> num_fmesh(3,0);
  for(int i=0;i<3;i++){
    int tmp=1;
    if(i<=Ndim){
      tmp=mi.GetFMesh(i);
    };
    num_fmesh[i]=tmp;
  };

  for(int i=0;i<3;i++){
    int axis1=1;
    int axis2=2;
    if(i==1)axis1=0; 
    if(i==2){ 
      axis1=0;
      axis2=1;
    };
    int fic_x=num_fmesh[axis1];
    int fic_y=num_fmesh[axis2];

    edge_e[i].resize(fic_x+1);
    for(int j=0;j<fic_x+1;j++){
      edge_e[i][j].resize(fic_y+1);
      for(int k=0;k<fic_y+1;k++){
        edge_e[i][j][k].resize(2,-1);
      };
    };
    for(int j=0;j<fic_x+1;j++){
      for(int k=0;k<fic_y+1;k++){

	if(i>Ndim-1){
	  edge_e[i][j][k][0]=0;
	  edge_e[i][j][k][1]=0;
	}else{

	int min=100000;
	int max=-1;
	if(j!=0&&k!=0){
	  int tmp=edge[i][j-1][k-1][0];
	  int tmp2=edge[i][j-1][k-1][1];
	  if(tmp<min&&tmp!=-1)min=tmp;
	  if(tmp2>max)max=tmp2;
	};
	if(j!=0&&k!=fic_y){
	  int tmp=edge[i][j-1][k][0];
	  int tmp2=edge[i][j-1][k][1];
	  if(tmp<min&&tmp!=-1)min=tmp;
	  if(tmp2>max)max=tmp2;
	};
	if(j!=fic_x&&k!=0){
	  int tmp=edge[i][j][k-1][0];
	  int tmp2=edge[i][j][k-1][1];
	  if(tmp<min&&tmp!=-1)min=tmp;
	  if(tmp2>max)max=tmp2;
	};
	if(j!=fic_x&&k!=fic_y){
	  int tmp=edge[i][j][k][0];
	  int tmp2=edge[i][j][k][1];
	  if(tmp<min&&tmp!=-1)min=tmp;
	  if(tmp2>max)max=tmp2;
	};
	edge_e[i][j][k][0]=min;
	edge_e[i][j][k][1]=max+1;

	if(edge_e[i][j][k][0]>10000){
	  edge_e[i][j][k][0]=-1;
	};
	if(edge_e[i][j][k][0]<0){
	  edge_e[i][j][k][1]=edge_e[i][j][k][0]-1;
	};

	};
      };
    };
  };

  int zf1=1;  
  if(Ndim==3)zf1=zf+1;

  if(Ndim!=1){
    TotME=(xf+1)*(yf+1)*zf1;
  }else{
    TotME=xf+1;
  };
};

void PLOSESystem::CalCoef()
{
  coefs.resize(grp);
  coefx.resize(grp);
  coefy.resize(grp);
  if(Ndim==3)coefz.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotME,0.); 
    coefx[i].resize(TotME,0.); 
    coefy[i].resize(TotME,0.); 
    if(Ndim==3)coefz[i].resize(TotME,0.);  
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };

  if(Sphere){
    CalCoefSphere();
    return;
  };
  if(Cylinder){
    CalCoefCylinder();
    return;
  };

  int x1=mi.GetXF()+1;
  int y1=2;
  if(Ndim>1)y1=mi.GetYF()+1;
  int z1=2;
  if(Ndim>2)z1=mi.GetZF()+1;
  int xy1=x1*y1;

  GroupData1D sigr(grp);
  real volfra=0.25;
  real suffra=0.5;
  real fact=0.4692;
  if(Ndim==1){
    volfra=0.5;
    suffra=1.0;
  };
  if(Ndim==3){
    volfra=0.125;
    suffra=0.25;
  };

  for(int z=0;z<z1-1;z++){
    for(int y=0;y<y1-1;y++){
      for(int x=0;x<x1-1;x++){
	int i=meshid[z][y][x];
	if(i!=-1){
  	  sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
	  real v8=mesh[i].GetVolume()*volfra;
	  for(int g=0;g<grp;g++){
	    real tmp=sigr.get_dat(g)*v8;
	    coefs[g][z    *xy1+y    *x1+x    ]+=tmp;     
	    coefs[g][z    *xy1+y    *x1+(x+1)]+=tmp;
            if(Ndim>1){
  	      coefs[g][z    *xy1+(y+1)*x1+x    ]+=tmp;
	      coefs[g][z    *xy1+(y+1)*x1+(x+1)]+=tmp;
	      if(Ndim>2){
  	        coefs[g][(z+1)*xy1+y    *x1+x    ]+=tmp;
	        coefs[g][(z+1)*xy1+(y+1)*x1+x    ]+=tmp;
	        coefs[g][(z+1)*xy1+y    *x1+(x+1)]+=tmp;
  	        coefs[g][(z+1)*xy1+(y+1)*x1+(x+1)]+=tmp;
	      };
	    };
	  };
	  real xl=mesh[i].GetLen(0);
          real yl=1.;
          if(Ndim>1)yl=mesh[i].GetLen(1);
	  real zl=1.;
	  if(Ndim>2)zl=mesh[i].GetLen(2);
          for(int g=0;g<grp;g++){
	    real tmp1=mesh[i].GetDif(g,difc[0])/xl*(yl*zl*suffra); // x-
            real tmp2=0.;
            if(Ndim>1)tmp2=mesh[i].GetDif(g,difc[1])/yl*(xl*zl*suffra); // y-
	    real tmp3=0.;
	    if(Ndim>2)tmp3=mesh[i].GetDif(g,difc[2])/zl*(xl*yl*suffra); // z-
	    real tsum=tmp1+tmp2+tmp3;
	    coefs[g][z    *xy1+y    *x1+x    ]+=tsum;     
	    coefs[g][z    *xy1+y    *x1+(x+1)]+=tsum;
  	    coefx[g][z    *xy1+y    *x1+x    ]-=tmp1;  
            if(Ndim>1){
  	      coefs[g][z    *xy1+(y+1)*x1+x    ]+=tsum;
	      coefs[g][z    *xy1+(y+1)*x1+(x+1)]+=tsum;
	      coefx[g][z    *xy1+(y+1)*x1+x    ]-=tmp1;
	      coefy[g][z    *xy1+y    *x1+x    ]-=tmp2;     
	      coefy[g][z    *xy1+y    *x1+(x+1)]-=tmp2;
	      if(Ndim>2){
  	        coefs[g][(z+1)*xy1+y    *x1+x    ]+=tsum;
  	        coefs[g][(z+1)*xy1+(y+1)*x1+x    ]+=tsum;
	        coefs[g][(z+1)*xy1+y    *x1+(x+1)]+=tsum;
  	        coefs[g][(z+1)*xy1+(y+1)*x1+(x+1)]+=tsum;
  	        coefx[g][(z+1)*xy1+y    *x1+x    ]-=tmp1;
  	        coefx[g][(z+1)*xy1+(y+1)*x1+x    ]-=tmp1;
  	        coefy[g][(z+1)*xy1+y    *x1+x    ]-=tmp2;
	        coefy[g][(z+1)*xy1+y    *x1+(x+1)]-=tmp2;
  	        coefz[g][z    *xy1+y    *x1+x    ]-=tmp3;     
	        coefz[g][z    *xy1+(y+1)*x1+x    ]-=tmp3;
	        coefz[g][z    *xy1+y    *x1+(x+1)]-=tmp3;
	        coefz[g][z    *xy1+(y+1)*x1+(x+1)]-=tmp3;
	      };
	    };
	  };

	  // Vacuum boundary for Plane 1 (X-)
	  int pos_ref=edge[0][y][z][0];
	  if(BC[0]==2&&x==pos_ref){
  	    real tmp=fact*(yl*zl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref<=edge_e[0][y][z][0])judge1=true;
	    if(Ndim>1&&pos_ref<=edge_e[0][y+1][z][0])judge2=true;
	    if(Ndim==3&&pos_ref<=edge_e[0][y][z+1][0])judge3=true;
	    if(Ndim==3&&pos_ref<=edge_e[0][y+1][z+1][0])judge4=true;
	    int tmp1=z*xy1+y*x1+x;
	    int tmp2=z*xy1+(y+1)*x1+x;
	    int tmp3=(z+1)*xy1+y*x1+x;
	    int tmp4=(z+1)*xy1+(y+1)*x1+x;
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	  // Vacuum boundary for Plane 2 (X+)
	  pos_ref=edge[0][y][z][1];
	  if(BC[1]==2&&x==pos_ref){
  	    real tmp=fact*(yl*zl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref+1>=edge_e[0][y][z][1])judge1=true;
	    if(Ndim>1&&pos_ref+1>=edge_e[0][y+1][z][1])judge2=true;
	    if(Ndim==3&&pos_ref+1>=edge_e[0][y][z+1][1])judge3=true;
	    if(Ndim==3&&pos_ref+1>=edge_e[0][y+1][z+1][1])judge4=true;
	    int tmp1=z*xy1+y*x1+x+1;
	    int tmp2=z*xy1+(y+1)*x1+x+1;
	    int tmp3=(z+1)*xy1+y*x1+x+1;
	    int tmp4=(z+1)*xy1+(y+1)*x1+x+1;
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	  // Vacuum boundary for Plane 3 (Y-)
	  pos_ref=edge[1][x][z][0];
	  if(Ndim>1&&BC[2]==2&&y==pos_ref){
  	    real tmp=fact*(xl*zl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref<=edge_e[1][x][z][0])judge1=true;
	    if(pos_ref<=edge_e[1][x+1][z][0])judge2=true;
	    if(Ndim==3&&pos_ref<=edge_e[1][x][z+1][0])judge3=true;
	    if(Ndim==3&&pos_ref<=edge_e[1][x+1][z+1][0])judge4=true;
	    int tmp1=z*xy1+y*x1+x;
	    int tmp2=z*xy1+y*x1+(x+1);
	    int tmp3=(z+1)*xy1+y*x1+x;
	    int tmp4=(z+1)*xy1+y*x1+(x+1);
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	  // Vacuum boundary for Plane 4 (Y+)
	  pos_ref=edge[1][x][z][1];
	  if(Ndim>1&&BC[3]==2&&y==pos_ref){
  	    real tmp=fact*(xl*zl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref+1>=edge_e[1][x][z][1])judge1=true;
	    if(pos_ref+1>=edge_e[1][x+1][z][1])judge2=true;
	    if(Ndim==3&&pos_ref+1>=edge_e[1][x][z+1][1])judge3=true;
	    if(Ndim==3&&pos_ref+1>=edge_e[1][x+1][z+1][1])judge4=true;
	    int tmp1=z*xy1+(y+1)*x1+x;
	    int tmp2=z*xy1+(y+1)*x1+x+1;
	    int tmp3=(z+1)*xy1+(y+1)*x1+x;
	    int tmp4=(z+1)*xy1+(y+1)*x1+x+1;
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	  // Vacuum boundary for Plane 5 (Z-)
	  pos_ref=edge[2][x][y][0];
	  if(Ndim>2&&BC[4]==2&&z==pos_ref){
  	    real tmp=fact*(xl*yl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref<=edge_e[2][x][y][0])judge1=true;
	    if(pos_ref<=edge_e[2][x+1][y][0])judge2=true;
	    if(pos_ref<=edge_e[2][x][y+1][0])judge3=true;
	    if(pos_ref<=edge_e[2][x+1][y+1][0])judge4=true;
	    int tmp1=z*xy1+y*x1+x;
	    int tmp2=z*xy1+y*x1+(x+1);
	    int tmp3=z*xy1+(y+1)*x1+x;
	    int tmp4=z*xy1+(y+1)*x1+(x+1);
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	  // Vacuum boundary for Plane 6 (Z+)
	  pos_ref=edge[2][x][y][1];
	  if(Ndim>2&&BC[5]==2&&z==pos_ref){
  	    real tmp=fact*(xl*yl*suffra);
	    bool judge1=false;
	    bool judge2=false;
	    bool judge3=false;
	    bool judge4=false;
  	    if(pos_ref+1>=edge_e[2][x][y][1])judge1=true;
	    if(pos_ref+1>=edge_e[2][x+1][y][1])judge2=true;
	    if(pos_ref+1>=edge_e[2][x][y+1][1])judge3=true;
	    if(pos_ref+1>=edge_e[2][x+1][y+1][1])judge4=true;
	    int tmp1=(z+1)*xy1+y*x1+x;
	    int tmp2=(z+1)*xy1+y*x1+x+1;
	    int tmp3=(z+1)*xy1+(y+1)*x1+x;
	    int tmp4=(z+1)*xy1+(y+1)*x1+x+1;
	    for(int g=0;g<grp;g++){
	      if(judge1)coefs[g][tmp1]+=tmp;
	      if(judge2)coefs[g][tmp2]+=tmp;
	      if(judge3)coefs[g][tmp3]+=tmp;
	      if(judge4)coefs[g][tmp4]+=tmp;
	    };
	  };
	};
      };
    };
  };

  DetermineSizeCoef();
};

void PLOSESystem::CalOmegaInFixedSource()
{
  for(int i=0;i<TotM;i++){
    mesh[i].PutFissionSrc(0.);
  };

  for(int g=0;g<grp;g++){
    CalSrcMultiplySystem(g,1.,0);
    CalFluxOmega(g,1e-1); // Initial guess for omega
    CalFluxOmega(g,1e-3); // The second (and final) guess for omega
    AddDownScatSrc(g,0);
  };
};

void PLOSESystem::CalCoefSphere()
{
  int x1=mi.GetXF()+1;

  real fact=4./3.*PI;

  GroupData1D vsigr(grp);
  real xlt=0.;
  for(int x=0;x<x1-1;x++){
    vsigr=mesh[x].GetMed()->GetMacxs().GetData1d(sigr);
    real xl=mesh[x].GetLen(0);
    real xrt=xlt+xl;
    real xmid=xlt+xl*0.5;
    real volr=fact*(xrt*xrt*xrt-xmid*xmid*xmid);
    real vol=mesh[x].GetVolume();
    real voll=vol-volr;
    for(int g=0;g<grp;g++){
      real tmp=vsigr.get_dat(g);
      coefs[g][x  ]+=tmp*voll;
      coefs[g][x+1]+=tmp*volr;
    };
    real xsur=4.*PI*xmid*xmid;
    for(int g=0;g<grp;g++){
      real tmp1=mesh[x].GetDif(g,difc[0])/xl*xsur; // x-
      coefs[g][x  ]+=tmp1;
      coefs[g][x+1]+=tmp1;
      coefx[g][x]-=tmp1;     
      if(x==x1-2&&BC[1]==2){
        real xsur=PI4*xrt*xrt;
        real tmp=0.4692*xsur;
	for(int g=0;g<grp;g++){
	  coefs[g][x+1]+=tmp;
	};
      };
      xlt=xrt;
    };
  };
  DetermineSizeCoef();
};

void PLOSESystem::CalCoefCylinder()
{
  int x1=mi.GetXF()+1;
  int y1=mi.GetYF()+1;

  GroupData1D vsigr(grp);
  //int i=0;
  for(int y=0;y<y1-1;y++){
    real xlt=0.;
    for(int x=0;x<x1-1;x++){
      int i=meshid[0][y][x];
      if(i!=-1){
	vsigr=mesh[i].GetMed()->GetMacxs().GetData1d(sigr);
        real xl=mesh[i].GetLen(0);
        real yl=mesh[i].GetLen(1);
        real xrt=xlt+xl;
        real xmid=xlt+xl*0.5;
        real volr=PI*(xrt*xrt-xmid*xmid)*yl;
        real vol=mesh[i].GetVolume();
        real voll=vol-volr;
        volr/=2.;
        voll/=2.;
        for(int g=0;g<grp;g++){
	  real tmp=vsigr.get_dat(g);
          coefs[g][y    *x1+x    ]+=tmp*voll;
          coefs[g][(y+1)*x1+x    ]+=tmp*voll;
          coefs[g][y    *x1+(x+1)]+=tmp*volr;
          coefs[g][(y+1)*x1+(x+1)]+=tmp*volr;
        };
        real xsur=PI2*xmid*yl*0.5;
        real ysurr=PI*(xrt*xrt-xmid*xmid);
        real ysurl=PI*(xmid*xmid-xlt*xlt);
        for(int g=0;g<grp;g++){
  	  real tmp1=mesh[i].GetDif(g,difc[0])/xl*xsur; // x-
	  real tmp2=mesh[i].GetDif(g,difc[1])/yl*ysurr; // y-
	  real tmp3=mesh[i].GetDif(g,difc[1])/yl*ysurl; // y+
          coefs[g][y    *x1+x    ]+=tmp1+tmp3;     
	  coefs[g][(y+1)*x1+x    ]+=tmp1+tmp3;
	  coefs[g][y    *x1+(x+1)]+=tmp1+tmp2;
	  coefs[g][(y+1)*x1+(x+1)]+=tmp1+tmp2;
	  coefx[g][y    *x1+x    ]-=tmp1;     
	  coefx[g][(y+1)*x1+x    ]-=tmp1;
	  coefy[g][y    *x1+x    ]-=tmp3;
	  coefy[g][y    *x1+(x+1)]-=tmp2;
        };

	// Vacuum boundary for Plane 2 (X+)
        int pos_ref=edge[0][y][0][1];	
	//if(x==x1-2&&BC[1]==2){
        if(x==pos_ref&&BC[1]==2){	
          real xsur=PI2*xrt*yl*0.5;
          real tmp=0.4692*xsur;
	  for(int g=0;g<grp;g++){
	    coefs[g][y    *x1+(x+1)]+=tmp;
	    coefs[g][(y+1)*x1+(x+1)]+=tmp;
	  };
        };

	// Vacuum boundary for Plane 3 (Y-)
        pos_ref=edge[1][x][0][0];
        //if(y==0&&BC[2]==2){
        if(y==pos_ref&&BC[2]==2){	
	  for(int g=0;g<grp;g++){
	    coefs[g][y*x1+x]+=0.4692*ysurl;
	    coefs[g][y*x1+(x+1)]+=0.4692*ysurr;
          };
        };

	// Vaccum boundary for Plane 4 (Y+)
        pos_ref=edge[1][x][0][1];      
        //if(y==y1-2&&BC[3]==2){
        if(y==pos_ref&&BC[3]==2){	
	  for(int g=0;g<grp;g++){
            coefs[g][(y+1)*x1+x    ]+=0.4692*ysurl;
	    coefs[g][(y+1)*x1+(x+1)]+=0.4692*ysurr;
	  };
        };
        xlt=xrt;
      };
    };
  };
  DetermineSizeCoef();
};

real PLOSESystem::CalFlux(int g,int itout,real epsif)
{
  int itmax=999;

  real *SrcE=new real[TotME];
  real *FlE =new real[TotME];
  for(int i=0;i<TotME;i++){SrcE[i]=0.; FlE[i]=0.;};

  SetEdgeSource(SrcE);
  SetInitialFluxInnerIteration(g,FlE);
  SweepInnerIteration(itmax,epsif,g,FlE,SrcE);
  real errf=RenewFluxInnerIteration(g,FlE);

  delete [] FlE;
  delete [] SrcE;
  return errf;
};

real PLOSESystem::CalFluxOmega(int g,real epsif)
{
  int itmax=299;

  // *** for FluxOmega ***
  int xm=mi.GetXF();
  int ym=mi.GetYF();
  int zm=mi.GetZF();
  if(Ndim==2)zm=0;
  if(Ndim==1){
    zm=0;
    ym=0;
  };
  int xm1=xm+1;
  int ym1=ym+1;
  int zm1=zm+1;
  int xym1=xm1*ym1;
  real *Res= new real[TotME]; 
  real tt1,tt2; 
  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  //real omega=1.; 
  // ********************

  real *SrcE=new real[TotME];
  real *FlE =new real[TotME];
  for(int i=0;i<TotME;i++){SrcE[i]=0.; FlE[i]=0.;};
  SetEdgeSource(SrcE);
  SetInitialFluxInnerIteration(g,FlE);

  real *tmparray=new real[size_coef];
  CalInnerCoef(g,tmparray);

  bool convergence=SweepInnerIteration(itmax,epsif,g,FlE,SrcE);

  if(convergence){

  // ***  Calculation for omega
  itmax=2;
  for(int iter=0;iter<itmax;iter++){
    int id=0;
    tt1=0.;    //fluxomega
    int cnt=0; //fluxomega
    int idt=0;
    for(int z=0;z<zm1;z++){
      for(int y=0;y<ym1;y++){
	int xel=edge_e[0][y][z][0];
	int xer=edge_e[0][y][z][1];
	if(xel>=0){
	//int xel=xedgel_e[y];
	//int xer=xedger_e[y];
	int xnum=xer-xel+1;
	real *a=new real[xnum*3];
	real *b=new real[xnum];
	id+=xel;
	int tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flt=SrcE[id];
          if(y!=edge_e[1][x][z][0]){flt-=tmparray[idt++]*FlE[id-xm1];};
          if(y!=edge_e[1][x][z][1]){flt-=tmparray[idt++]*FlE[id+xm1];};
          if(z!=edge_e[2][x][y][0]){flt-=tmparray[idt++]*FlE[id-xym1];};
          if(z!=edge_e[2][x][y][1]){flt-=tmparray[idt++]*FlE[id+xym1];};
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =tmparray[idt++];};
          if(x!=xer){a[tmp*3+2]=tmparray[idt++];};
	  a[tmp*3+1]=tmparray[idt++];
	  tmp++;
	  id++;
	};
        gauss_3p(a,b,xnum,1);
	id-=xnum;
	tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flo=FlE[id];
	  real fln=b[tmp];
	  // ** for FluxOmega **
	  real reso=Res[id]; 
	  real resn=pow(fln-flo,2); 
	  Res[id]=resn; 
	  tt2=sqrt(resn/reso);
	  if(tt2>0&&tt2<1){
	    tt1+=tt2;
	    cnt++;
	  };
	  // *******************
	  FlE[id]=flo+omega*(fln-flo);
	  id++;
	  tmp++;
	};
	id+=xm-xer;
	delete [] a;
	delete [] b;
	}else{
	  id+=xm1;
	};
      };
    };
    tt1/=cnt; //fluxomega
    if(cnt==0)tt1=0.; // fluxomega
  };
  if(tt1!=0.){
    real rnew=(omega-1+tt1)/omega/sqrt(tt1);
    rnew=rnew*rnew;
    if(rnew>1.){
      omega=1.;
    }else{
      omega=2.0/(1+sqrt(1-rnew)); 
    //if(omega>1.3)omega=1.3;
    //cout<<"  omega in group "<<g<<" : "<<omega<<"\n";
    };
  }else{
    omega=1.;
  };
  opt.PutOmegai(g,omega);
  //cout<<"   ** omega for inner iteration in group "<<g<<" : "<<omega<<" "<<tt1<<"\n";
  if(opt.Print())cout<<"   ** omega for inner iteration : "<<omega<<"\n";
  // *** end

  };

  real errf=RenewFluxInnerIteration(g,FlE);

  delete [] Res; // fluxomega
  delete [] FlE;
  delete [] SrcE;
  delete [] tmparray;
  return errf;
};

void PLOSESystem::CalFluxAutoConv(int g,real epsif)
{
  int xm  =mi.GetXF();
  int ym  =mi.GetYF();
  int zm  =mi.GetZF();
  if(Ndim==2)zm=0;
  if(Ndim==1){
    zm=0;
    ym=0;
  };
  int xm1 =xm+1;
  int ym1 =ym+1;
  int zm1 =zm+1;
  int xym1=xm1*ym1;
  real *SrcE=new real[TotME];
  real *FlE =new real[TotME];

  real *tmparray=new real[size_coef];
  CalInnerCoef(g,tmparray);

  real omega=opt.GetOmegai(g);
  int itmax  =int(log10(epsif)/log10(omega-1.))+1;

  for(int i=0;i<TotME;i++){SrcE[i]=0.; FlE[i]=0.;};
  SetEdgeSource(SrcE);

  for(int iter=0;iter<itmax;iter++){
    int id=0;
    int idt=0;
    for(int z=0;z<zm1;z++){
      for(int y=0;y<ym1;y++){
	//int xel=xedgel_e[y];
	//int xer=xedger_e[y];
	int xel=edge_e[0][y][z][0];
	int xer=edge_e[0][y][z][1];
	if(xel>=0){
	int xnum=xer-xel+1;
	real *a=new real[xnum*3];
	real *b=new real[xnum];
	id+=xel;
	int tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flt=SrcE[id];
          //if(y!=yedgel_e[x]){flt-=tmparray[idt++]*FlE[id-xm1];};
          //if(y!=yedger_e[x]){flt-=tmparray[idt++]*FlE[id+xm1];};
          if(y!=edge_e[1][x][z][0]){flt-=tmparray[idt++]*FlE[id-xm1];};
          if(y!=edge_e[1][x][z][1]){flt-=tmparray[idt++]*FlE[id+xm1];};
          //if(z!=0) {flt-=tmparray[idt++]*FlE[id-xym1];};
          //if(z!=zm){flt-=tmparray[idt++]*FlE[id+xym1];};
          if(z!=edge_e[2][x][y][0]) {flt-=tmparray[idt++]*FlE[id-xym1];};
          if(z!=edge_e[2][x][y][1]){flt-=tmparray[idt++]*FlE[id+xym1];};
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =tmparray[idt++];};
          if(x!=xer){a[tmp*3+2]=tmparray[idt++];};
	  a[tmp*3+1]=tmparray[idt++];
	  tmp++;
	  id++;
	};
        gauss_3p(a,b,xnum,1);
	id-=xnum;
	tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flo=FlE[id];
	  FlE[id]=flo+omega*(b[tmp]-flo);
	  id++;
	  tmp++;
	};
        id+=xm-xer;
	delete [] a;
	delete [] b;
	}else{
	  id+=xm1;
	};
      };
    };
  };

  RenewFluxInnerIteration(g,FlE);

  delete [] FlE;
  delete [] SrcE;
  delete [] tmparray;
};

bool PLOSESystem::SweepInnerIteration(int itmax,real epsif,int g,real *FlE,real *SrcE)
{
  bool convergence=false;

  int xm  =mi.GetXF();
  int ym  =mi.GetYF();
  int zm  =mi.GetZF();
  if(Ndim==2)zm=0;
  if(Ndim==1){
    zm=0;
    ym=0;
  };
  int xm1 =xm+1;
  int ym1 =ym+1;
  int zm1 =zm+1;
  int xym1=xm1*ym1;

  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };

  real *tmparray=new real[size_coef];

  CalInnerCoef(g,tmparray);
  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    int id=0;
    int idt=0;
    for(int z=0;z<zm1;z++){
      for(int y=0;y<ym1;y++){
	int xel=edge_e[0][y][z][0];
	int xer=edge_e[0][y][z][1];
	if(xel>=0){
	int xnum=xer-xel+1;
	real *a=new real[xnum*3];
	real *b=new real[xnum];
	id+=xel;
	int tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flt=SrcE[id];
          if(y!=edge_e[1][x][z][0]){flt-=tmparray[idt++]*FlE[id-xm1];};
          if(y!=edge_e[1][x][z][1]){flt-=tmparray[idt++]*FlE[id+xm1];};
          if(z!=edge_e[2][x][y][0]){flt-=tmparray[idt++]*FlE[id-xym1];};
          if(z!=edge_e[2][x][y][1]){flt-=tmparray[idt++]*FlE[id+xym1];};
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =tmparray[idt++];};
          if(x!=xer){a[tmp*3+2]=tmparray[idt++];};
	  a[tmp*3+1]=tmparray[idt++];
	  tmp++;
	  id++;
	};
        gauss_3p(a,b,xnum,1);
	id-=xnum;
	tmp=0;
	for(int x=xel;x<=xer;x++){
	  real flo=FlE[id];
	  real fln=b[tmp];
	  if(errmax<epsif){
 	    real err=fabs(fln/flo-1.);
	    if(err>errmax)errmax=err;
	  };
	  FlE[id]=flo+omega*(fln-flo);
	  id++;
	  tmp++;
	};
	id+=xm-xer;
	delete [] a;
	delete [] b;
	}else{
	  id+=xm1;
	};
      };
    };
    if(errmax<epsif&&iter>0){
      //cout<<"  (group:"<<g<<") "<<iter<<" : "<<errmax<<"\n";
      convergence=true;
      break;
    };
    if(iter==itmax-1&&itmax!=1){
      cout<<"#Inner iteration is exceeded for "<<itmax;
      cout<<"#  (group:"<<g<<")(errmax:"<<errmax<<")\n";
    };
  };

  delete [] tmparray;

  return convergence;
};

void PLOSESystem::DetermineSizeCoef()
{
  int ym  =mi.GetYF();
  int zm  =mi.GetZF();
  if(Ndim==2)zm=0;
  if(Ndim==1){
    zm=0;
    ym=0;
  };
  int ym1 =ym+1;
  int zm1 =zm+1;
  int idt=0;
  for(int z=0;z<zm1;z++){
    for(int y=0;y<ym1;y++){
      int xel=edge_e[0][y][z][0];
      int xer=edge_e[0][y][z][1];
      if(xel>=0){
      for(int x=xel;x<=xer;x++){
        if(y!=edge_e[1][x][z][0])idt++;
        if(y!=edge_e[1][x][z][1])idt++;
        if(z!=edge_e[2][x][y][0])idt++;
        if(z!=edge_e[2][x][y][1])idt++;
        if(x!=xel) idt++;
	if(x!=xer) idt++;
	idt++;
      };
      };
    };
  };
  size_coef=idt;
};

void PLOSESystem::CalInnerCoef(int g,real *tmparray)
{
  int xm  =mi.GetXF();
  int ym  =mi.GetYF();
  int zm  =mi.GetZF();
  if(Ndim==2)zm=0;
  if(Ndim==1){
    zm=0;
    ym=0;
  };
  int xm1 =xm+1;
  int ym1 =ym+1;
  int zm1 =zm+1;
  int xym1=xm1*ym1;

  int idt=0;
  int id=0;
  for(int z=0;z<zm1;z++){
    for(int y=0;y<ym1;y++){
      int xel=edge_e[0][y][z][0];
      int xer=edge_e[0][y][z][1];
      if(xel>=0){
      id+=xel;
      for(int x=xel;x<=xer;x++){
        if(y!=edge_e[1][x][z][0])tmparray[idt++]=coefy[g][id-xm1];
        if(y!=edge_e[1][x][z][1])tmparray[idt++]=coefy[g][id];
        if(z!=edge_e[2][x][y][0])tmparray[idt++]=coefz[g][id-xym1];
        if(z!=edge_e[2][x][y][1])tmparray[idt++]=coefz[g][id];
        if(x!=xel) tmparray[idt++]=coefx[g][id-1];
        if(x!=xer) tmparray[idt++]=coefx[g][id];
	tmparray[idt++]=coefs[g][id];
        id++;
      };
      id+=xm-xer;
      }else{
	id+=xm1;
      };
    };
  };
};

real PLOSESystem::RenewFluxInnerIteration(int g,real *FlE)
{
  if(Cylinder)return RenewFluxInnerIterationCylinder(g,FlE);
  if(Sphere)  return RenewFluxInnerIterationSphere(g,FlE);

  int xm  =mi.GetXF();
  int xm1 =xm+1;
  int ym=1;
  if(Ndim>1)ym=mi.GetYF();
  int xym1=xm1*(ym+1);
  int zm=1;
  if(Ndim==3)zm=mi.GetZF();

  real volfra=0.25;
  if(Ndim==3)volfra=0.125;
  if(Ndim==1)volfra=0.5;

  real errf=0.;
  for(int i=0;i<zm;i++){
    for(int j=0;j<ym;j++){
      for(int k=0;k<xm;k++){
	int id=meshid[i][j][k];
	if(id!=-1){
  	  int tag=i*xym1+j*xm1+k;
	  real tmp=FlE[tag]+FlE[tag+1];
          if(Ndim>1){
	    tmp+=FlE[tag+xm1]+FlE[tag+xm1+1];
	    if(Ndim>2){
              tmp+=FlE[tag+xym1]+FlE[tag+xym1+xm1]
 	          +FlE[tag+xym1+1]+FlE[tag+xym1+xm1+1];
	    };
	  };
          tmp*=volfra;
	  real err=fabs(tmp/mesh[id].GetFlux().get_dat(g)-1.0);
	  if(err>errf)errf=err;
          mesh[id].GetFlux().put_data(g,tmp);
	};
      };
    };
  };
  return errf;
};

real PLOSESystem::RenewFluxInnerIterationSphere(int g,real *FlE)
{
  int xm  =mi.GetXF();

  real fact=4./3.*PI;
  real errf=0.;
  real xlt=0.;
  for(int k=0;k<xm;k++){
    real vol=mesh[k].GetVolume();
    real xl=mesh[k].GetLen(0);
    real xmid=xlt+xl*0.5;
    real xrt=xlt+xl;
    real voll=fact*(xmid*xmid*xmid-xlt*xlt*xlt);
    real volr=vol-voll;
    real tmp=(FlE[k]*voll+FlE[k+1]*volr)/vol;
    real err=fabs(tmp/mesh[k].GetFlux().get_dat(g)-1.0);
    if(err>errf)errf=err;
    mesh[k].GetFlux().put_data(g,tmp);
    xlt=xrt;
  };
  return errf;
};

real PLOSESystem::RenewFluxInnerIterationCylinder(int g,real *FlE)
{
  int xm  =mi.GetXF();
  int ym  =mi.GetYF();
  int xm1 =mi.GetXF()+1;

  real errf=0.;
  int id=0;
  for(int j=0;j<ym;j++){
    real xlt=0.;
    for(int k=0;k<xm;k++){
      int id=meshid[0][j][k];
      if(id!=-1){
        real vol=mesh[id].GetVolume();
        real xl=mesh[id].GetLen(0);
        real yl=mesh[id].GetLen(1);
        real xmid=xlt+xl*0.5;
        real xrt=xlt+xl;
        real voll=PI*(xmid*xmid-xlt*xlt)*yl*0.5;
        real volr=PI*(xrt*xrt-xmid*xmid)*yl*0.5;
        int tag=j*xm1+k;
        real tmp=( FlE[tag]*voll+
                   FlE[tag+xm1]*voll+
                   FlE[tag+1]*volr+
                   FlE[tag+xm1+1]*volr)/vol;
        real err=fabs(tmp/mesh[id].GetFlux().get_dat(g)-1.0);
        if(err>errf)errf=err;
        mesh[id].GetFlux().put_data(g,tmp);
        xlt=xrt;
      };
    };
  };
  return errf;
};

void PLOSESystem::SetEdgeSource(real *SrcE)
{
  if(Sphere){
    SetEdgeSourceSphere(SrcE);
    return;
  };
  if(Cylinder){
    SetEdgeSourceCylinder(SrcE);
    return;
  };

  real *Src=new real[TotM];
  for(int i=0;i<TotM;i++){
    Src[i]=mesh[i].GetSrcin()*PI4;
  };

  int xm =mi.GetXF();
  int xm1=xm+1;
  int ym=1;
  if(Ndim>1)ym=mi.GetYF();
  int xym1=xm1*(ym+1);
  int zm=1;
  if(Ndim>2)zm=mi.GetZF();

  real volfra=0.25;
  if(Ndim==3)volfra=0.125;
  if(Ndim==1)volfra=0.5;

  for(int i=0;i<zm;i++){
    for(int j=0;j<ym;j++){
      for(int k=0;k<xm;k++){
	int id=meshid[i][j][k];
	if(id!=-1){
  	  real tmp=Src[id]*volfra;
	  int tag=i*xym1+j*xm1+k;
	  SrcE[tag]+=  tmp;
	  SrcE[tag+1]+=tmp;
          if(Ndim>1){
	    SrcE[tag+xm1]+= tmp;
  	    SrcE[tag+xm1+1]+=tmp;
  	    if(Ndim>2){
  	      SrcE[tag+xym1]+=tmp;
	      SrcE[tag+xym1+xm1]+=tmp;
	      SrcE[tag+xym1+1]+=  tmp;
	      SrcE[tag+xym1+xm1+1]+=tmp;
	    };
	  };
	};
      };
    };
  };

  delete [] Src;
};

void PLOSESystem::SetEdgeSourceSphere(real *SrcE)
{
  real *Src=new real[TotM];
  for(int i=0;i<TotM;i++){
    Src[i]=mesh[i].GetSrcin()*PI4;
  };

  real fact=4./3.*PI;
  int xm  =mi.GetXF();
  real xlt=0.;
  for(int k=0;k<xm;k++){
    real vol=mesh[k].GetVolume();
    real xl=mesh[k].GetLen(0);
    real xmid=xlt+xl*0.5;
    real xrt=xlt+xl;
    real voll=fact*(xmid*xmid*xmid-xlt*xlt*xlt);
    real volr=vol-voll;
    real tmp=Src[k]/vol;
    SrcE[k]+=tmp*voll;
    SrcE[k+1]+=tmp*volr;
    xlt=xrt;
  };
  delete [] Src;
};

void PLOSESystem::SetEdgeSourceCylinder(real *SrcE)
{
  real *Src=new real[TotM];
  for(int i=0;i<TotM;i++){
    Src[i]=mesh[i].GetSrcin()*PI4;
  };

  int xm  =mi.GetXF();
  int ym  =mi.GetYF();
  int xm1 =mi.GetXF()+1;

  //int id=0;
  for(int j=0;j<ym;j++){
    real xlt=0.;
    for(int k=0;k<xm;k++){
      int id=meshid[0][j][k];
      if(id!=-1){
        real vol=mesh[id].GetVolume();
        real xl=mesh[id].GetLen(0);
        real yl=mesh[id].GetLen(1);
        real xmid=xlt+xl*0.5;
        real xrt=xlt+xl;
        real voll=PI*(xmid*xmid-xlt*xlt)*yl*0.5;
        real volr=PI*(xrt*xrt-xmid*xmid)*yl*0.5;
        real tmp=Src[id]/vol;
        int tag=j*xm1+k;
        SrcE[tag]+=tmp*voll;
        SrcE[tag+xm1]+=tmp*voll;
        SrcE[tag+1]+=tmp*volr;
        SrcE[tag+xm1+1]+=tmp*volr;
        xlt=xrt;
      };
    };
  };
  delete [] Src;
};

void PLOSESystem::SetInitialFluxInnerIteration(int g,real *FlE)
{
  if(Cylinder){
    SetInitialFluxInnerIterationCylinder(g,FlE);
    return;
  };
  if(Sphere){
    SetInitialFluxInnerIterationSphere(g,FlE);
    return;
  };

  int xm=mi.GetXF();
  int ym=1;
  if(Ndim>1)ym=mi.GetYF();
  int zm=1;
  if(Ndim>2)zm=mi.GetZF();
  int xm1=xm+1;
  int ym1=ym+1;
  int xym1=xm1*ym1;

  int *num=new int[TotME];
  for(int i=0;i<TotME;i++){num[i]=0;};
  for(int i=0;i<zm;i++){
    for(int j=0;j<ym;j++){
      for(int k=0;k<xm;k++){
	int id=meshid[i][j][k];
	if(id!=-1){
          int tag=i*xym1+j*xm1+k;
      	  real tmp=mesh[id].GetFlux().get_dat(g);
	  FlE[tag]+=  tmp;
	  FlE[tag+1]+=tmp;
	  num[tag]++;
	  num[tag+1]++;
          if(Ndim>1){
  	    FlE[tag+xm1]+=tmp;
	    FlE[tag+xm1+1]+=tmp;
	    num[tag+xm1]++;
  	    num[tag+xm1+1]++;
	    if(Ndim>2){
	      FlE[tag+xym1]+=tmp;
	      FlE[tag+xym1+xm1]+=tmp;
	      FlE[tag+xym1+1]+=  tmp;
	      FlE[tag+xym1+xm1+1]+=tmp;
	      num[tag+xym1]++;
	      num[tag+xym1+xm1]++;
	      num[tag+xym1+1]++;
	      num[tag+xym1+xm1+1]++;
	    };
	  };
	};
      };
    };
  };
  for(int i=0;i<TotME;i++){FlE[i]/=num[i];};
  delete [] num;
};

void PLOSESystem::SetInitialFluxInnerIterationSphere(int g,real *FlE)
{
  int xm=mi.GetXF();

  real *num=new real[TotME];
  for(int i=0;i<TotME;i++){num[i]=0.;};

  real fact=4./3.*PI;

  real xlt=0.;
  for(int k=0;k<xm;k++){
    real xl=mesh[k].GetLen(0);
    real xmid=xlt+xl*0.5;
    real xrt=xlt+xl;
    real voll=fact*(xmid*xmid*xmid-xlt*xlt*xlt);
    real volr=mesh[k].GetVolume()-voll;
    real tmp=mesh[k].GetFlux().get_dat(g);
    FlE[k]+=tmp*voll;
    FlE[k+1]+=tmp*volr;
    num[k]+=voll;
    num[k+1]+=volr;
    xlt=xrt;
  };
  for(int i=0;i<TotME;i++){FlE[i]/=num[i];};
  delete [] num;
};

void PLOSESystem::SetInitialFluxInnerIterationCylinder(int g,real *FlE)
{
  int xm=mi.GetXF();
  int ym=mi.GetYF();
  int xm1=mi.GetXF()+1;

  real *num=new real[TotME];
  for(int i=0;i<TotME;i++){num[i]=0.;};
  //int id=0;
  for(int j=0;j<ym;j++){
    real xlt=0.;
    for(int k=0;k<xm;k++){
      int id=meshid[0][j][k];
      if(id!=-1){		       
        real xl=mesh[id].GetLen(0);
        real yl=mesh[id].GetLen(1);
        real xmid=xlt+xl*0.5;
        real xrt=xlt+xl;
        real voll=PI*(xmid*xmid-xlt*xlt)*yl*0.5;
        real volr=PI*(xrt*xrt-xmid*xmid)*yl*0.5;
        int tag=j*xm1+k;
        real tmp=mesh[id].GetFlux().get_dat(g);
        FlE[tag]+=tmp*voll;
        FlE[tag+xm1]+=tmp*voll;
        FlE[tag+1]+=tmp*volr;
        FlE[tag+xm1+1]+=tmp*volr;
        num[tag]+=voll;
        num[tag+xm1]+=voll;
        num[tag+1]+=volr;
        num[tag+xm1+1]+=volr;
      };
    };
  };
  for(int i=0;i<TotME;i++){FlE[i]/=num[i];};
  delete [] num;
};

void PLOSESystem::SetInitialFlux()
{
  SetInitialFlatFlux();
};

real PLOSESystem::CalFluxGeneral(int grp,real epsif,int iter)
{
  real err;
  if(iter==0&&Ndim!=1){
    err=CalFluxOmega(grp,0.1);
    err=CalFluxOmega(grp,1e-3);
  }
  else{
    err=CalFlux(grp,iter,epsif);
  };
  return err;
};

void PLOSESystem::SetDifc(string d1, string d2, string d3)
{
  difc[0]=-1;
  if(d1=="d") difc[0]=0;
  if(d1=="dr"||d1=="dperp")difc[0]=1;
  if(d1=="dz"||d1=="dpara")difc[0]=2;
  if(Ndim>1){
    difc[1]=-1;
    if(d2=="d") difc[1]=0;
    if(d2=="dr"||d2=="dperp")difc[1]=1;
    if(d2=="dz"||d2=="dpara")difc[1]=2;
  };
  if(Ndim>2){
    difc[2]=-1;
    if(d3=="d") difc[2]=0;
    if(d3=="dr"||d3=="dperp")difc[2]=1;
    if(d3=="dz"||d3=="dpara")difc[2]=2;
  };
  for(int i=0;i<Ndim;i++){
    if(difc[i]<0||difc[i]>2){
      cout<<"Error in Diffusion coefficient setting\n";
      cout<<" Direction(1:x 2:y 3:z) : "<<i<<"\n";
      if(i==0)cout<<" You requested : "<<d1<<"\n";
      if(i==1)cout<<" You requested : "<<d2<<"\n";
      if(i==2)cout<<" You requested : "<<d3<<"\n";
      exit(0);
    };
  };
};

void PLOSESystem::AllVectorClear()
{
  GeneralSystem::AllVectorClear();

  coefs.clear();
  coefx.clear();
  coefy.clear();
  coefz.clear();
  xedgel_e.clear();
  xedger_e.clear();
  yedgel_e.clear();
  yedger_e.clear();
  edge_e.clear();
};
