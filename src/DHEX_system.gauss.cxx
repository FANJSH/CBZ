#include "DHEX_system.h"

using namespace std;

void DHEXSystem::Init(int n,int g,int i)
{
  GeneralSystem::Init(n,g,i);
  difc[0]  = 0;
  difc[1]  = 0;
  difc[2]  = 0;
  name="DHEX";
  PutPL(0);
  cmfdimp = true; // CMFD is implemented in PLOS
  cmfd_factor=1.0;
  cal_mcor = false;
  cal_tri  = false;

  if(Ndim!=3){
    cout<<"DHEX can be applied only into 3D calculations.\n";
    exit(0);
  };
};

void DHEXSystem::EdgeCalculationForHex()
{
  int xr=mi->GetXF();
  int yr=mi->GetYF();

  ledge_yl.resize(xr);
  ledge_yr.resize(xr);
  redge_yl.resize(xr);
  redge_yr.resize(xr);
  for(int i=0;i<xr;i++){
    ledge_yl[i].resize(yr,false);
    ledge_yr[i].resize(yr,false);
    redge_yl[i].resize(yr,false);
    redge_yr[i].resize(yr,false);
  };

  for(int y=0;y<yr;y++){
    for(int x=0;x<xr;x++){
      if(meshid[0][y][x]!=-1){
      // YL direction
      bool ledge=false;
      if(y==0||(y%2==1&&x==0)){
	ledge=true;
      }else{
	int tmp=x;
	if(y%2==1)tmp--;
	if(meshid[0][y-1][tmp]==-1)ledge=true;
      };
      if(ledge)ledge_yl[x][y]=true;
      bool redge=false;
      if(y==yr-1||(y%2==0&&x==xr-1)){
	redge=true;
      }else{
	int tmp=x;
        if(y%2==0)tmp++;
	if(meshid[0][y+1][tmp]==-1)redge=true;
      };
      if(redge)redge_yl[x][y]=true;
      // YR direction
      ledge=false;
      if(y==0||(y%2==0&&x==xr-1)){
	ledge=true;
      }else{
	int tmp=x;
	if(y%2==0)tmp++;
	if(meshid[0][y-1][tmp]==-1)ledge=true;
      };
      if(ledge)ledge_yr[x][y]=true;
      redge=false;
      if(y==yr-1||(y%2==1&&x==0)){
	redge=true;
      }else{
	int tmp=x;
        if(y%2==1)tmp--;
	if(meshid[0][y+1][tmp]==-1)redge=true;
      };
      if(redge)redge_yr[x][y]=true;
      };
    };
  };
};

void DHEXSystem::PutCartMeshInfo(CartMeshInfo *cm, real pitchinp)
{
  pitch=pitchinp;
  len_plane=pitch*1.732050808*0.33333333;
  half_pitch=pitch*0.5;
  inv_pitch=1./pitch;
  hex_area=pitch*pitch*1.732050808*0.5;

  mi=cm;

  int xf=mi->GetXF();
  int yf=mi->GetYF();
  int zf=mi->GetZF();
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
	int tm=mi->GetFMat(tmp2++);
	if(tm<-1||tm>=nmed){
	  cout<<"Invalid material number for CartMeshInfo (Number:"<<tm<<")\n";
	  exit(0);
	};
	if(tm!=-1)tmp++;
      };
    };
  };
  TotM=tmp;
  if(print)cout<<"*** Total Mesh : "<<TotM<<"\n";
  mesh.resize(TotM);

  int index=0;
  int ind2=0;
  for(int z=0;z<zf;z++){
    for(int y=0;y<yf;y++){
      for(int x=0;x<xf;x++){
	int tm=mi->GetFMat(ind2);
	ind2++;
	if(tm!=-1){
	  real zz=mi->GetFMeshL(2,z);
	  mesh[index].PutDim(3);
	  mesh[index].PutLen(0,pitch);
	  mesh[index].PutLen(1,pitch);
	  mesh[index].PutLen(2,zz);
	  mesh[index].PutVolume(hex_area*zz);
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

  EdgeCalculation();
  EdgeCalculationForHex();

  SetCoarseMeshInfo();

  for(int i=0;i<Ndim*2;i++){
    if(mi->GetBC(i)==2)BC[i]=2; // Vacuum
    if(mi->GetBC(i)==1)BC[i]=1; // Reflective
    if(mi->GetBC(i)==0)BC[i]=0; // Zero flux
    if(BC[i]<0||BC[i]>2){
      cout<<"Boundary Condition Error in PLOS.\n";
      exit(0);
    };
  };
};

void DHEXSystem::CalCoef(bool mcorrect)
{
  coefs.resize(grp);
  coefx.resize(grp);
  coefyl.resize(grp);
  coefyr.resize(grp);
  coefz.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotM,0.);
    coefx[i].resize(TotM,0.);
    coefyl[i].resize(TotM,0.);
    coefyr[i].resize(TotM,0.);
    coefz[i].resize(TotM,0.);
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };

  GroupData1D sigr(grp);
  real d1[3];
  real d2,dl;
  int id;
  int xyz[3];
  bool ledge,redge;

  for(int z=0;z<mi->GetZF();z++){
    real zlen=mi->GetFMeshL(2,z);
    real half_zlen=zlen*0.5;
    real vol=zlen*hex_area;
    real surface_xy=zlen*len_plane;
    for(int y=0;y<mi->GetYF();y++){
      for(int x=0;x<mi->GetXF();x++){
	int i=meshid[z][y][x];
	if(i!=-1){
  	  xyz[0]=x;
	  xyz[1]=y;
	  xyz[2]=z;
	  sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
          for(int g=0;g<grp;g++){
            real tmp=sigr.get_dat(g)*vol;
	    if(cal_tri)tmp*=0.16666666666;
	    coefs[g][i]=tmp;
	    // mcorrect
	    if(mcorrect){
	      coefs[g][i]*=(1.-mcor[i][g]);
	    };
	    //
            for(int nn=0;nn<3;nn++){
  	      d1[nn]=mesh[i].GetDif(g,difc[nn]);
	      // mcorrect
	      if(mcorrect&&nn!=2){
		d1[nn]*=(1.-0.66666666*mcor[i][g]);
	      };
	      //
	    };
	    // *** X-direction
	    ledge=false;
	    if(x==edge[0][y][z][0])ledge=true;
            if(!ledge){
              d2=mesh[meshid[z][y][x-1]].GetDif(g,difc[0]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[meshid[z][y][x-1]][g]);
	      };
	      //
	      dl=(d1[0]+d2)*half_pitch;
	      real tmp=d1[0]*d2/dl*surface_xy;
              if(!cal_tri)coefs[g][i]+=tmp;
	      if(cal_tri)tmp*=3.;
              coefx[g][i]=-tmp;
            }else{
              real tmp=CalEdgeCurVacuum(i,0,g,difc[0],mcorrect)*surface_xy;
              if(!cal_tri)coefs[g][i]-=tmp;
	    };
	    //
            redge=false;
	    if(x==mi->GetXF()-1){
              redge=true;
            }else if(meshid[z][y][x+1]==-1){
	      redge=true;
	    };
	    if(!redge){
              d2=mesh[meshid[z][y][x+1]].GetDif(g,difc[0]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[meshid[z][y][x+1]][g]);
	      };
	      //
              dl=(d1[0]+d2)*half_pitch;
              if(!cal_tri)coefs[g][i]+=d1[0]*d2/dl*surface_xy;
            }else{
              if(!cal_tri)coefs[g][i]-=CalEdgeCurVacuum(i,0,g,difc[0],mcorrect)*surface_xy;
	    };
	    // *** Z-direction
            if(z!=0){
	      int id=meshid[z-1][y][x];
              d2=mesh[id].GetDif(g,difc[2]);
              dl=d1[2]*mesh[id].GetLen(2)*0.5+d2*half_zlen;
              real tmp=d1[2]*d2/dl*hex_area;
	      if(cal_tri)tmp*=0.166666666;
              coefz[g][i]=-tmp;
	      if(mcorrect){
		tmp*=(1.-mcor[i][g]);
	      };
              coefs[g][i]+=tmp;
	    }else{
              switch(BC[4]){
	        case 2:
                   real tmp=CalEdgeCurVacuum(i,2,g,difc[2])*hex_area;
		   if(mcorrect){
		     tmp*=(1.-mcor[i][g]);
		   };
		   if(cal_tri)tmp*=0.16666666666;
                   coefs[g][i]-=tmp;
		  break;
              };
	    };
            if(z!=mi->GetZF()-1){
              int id=meshid[z+1][y][x];
              d2=mesh[id].GetDif(g,difc[2]);
              dl=d1[2]*mesh[id].GetLen(2)*0.5+d2*half_zlen;
              real tmp=d1[2]*d2/dl*hex_area;
	      if(mcorrect){
		tmp*=(1.-mcor[i][g]);
	      };
	      if(cal_tri)tmp*=0.166666666666;
              coefs[g][i]+=tmp;
	    }else{
	      switch(BC[5]){
	        case 2:
                  real tmp=CalEdgeCurVacuum(i,2,g,difc[2]);
                  if(mcorrect){
		    tmp*=(1.-mcor[i][g]);
		  };
		  if(cal_tri)tmp*=0.1666666666;
                  coefs[g][i]-=tmp*hex_area;
		  break;
	      };
	    };
	    // *** Y-direction (Left)
            if(!ledge_yl[x][y]){
	      if(y%2==0){id=meshid[z][y-1][x];}
	      else if(x!=0){id=meshid[z][y-1][x-1];};
              d2=mesh[id].GetDif(g,difc[1]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[id][g]);
	      };
	      //
              dl=(d1[1]+d2)*half_pitch;
              real tmp=d1[1]*d2/dl*surface_xy;
              if(!cal_tri)coefs[g][i]+=tmp;
	      if(cal_tri)tmp*=3.;
              coefyl[g][i]=-tmp;
	    }else{
              if(!cal_tri)coefs[g][i]-=CalEdgeCurVacuum(i,1,g,difc[1],mcorrect)*surface_xy;
            };
	    //
            if(!redge_yl[x][y]){
	      if(y%2==1){id=meshid[z][y+1][x];}
	      else if(x!=mi->GetXF()-1){id=meshid[z][y+1][x+1];};
              d2=mesh[id].GetDif(g,difc[1]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[id][g]);
	      };
	      //
              dl=(d1[1]+d2)*half_pitch;
              if(!cal_tri)coefs[g][i]+=d1[1]*d2/dl*surface_xy;
            }else{
              if(!cal_tri)coefs[g][i]-=CalEdgeCurVacuum(i,1,g,difc[1],mcorrect)*surface_xy;
            };
	    // *** Y-direction (Right)
            if(!ledge_yr[x][y]){
	      if(y%2==1){id=meshid[z][y-1][x];}
	      else if(x!=mi->GetXF()-1){id=meshid[z][y-1][x+1];};
              d2=mesh[id].GetDif(g,difc[1]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[id][g]);
	      };
	      //

              dl=(d1[1]+d2)*half_pitch;
              real tmp=d1[1]*d2/dl*surface_xy;
              if(!cal_tri)coefs[g][i]+=tmp;
	      if(cal_tri)tmp*=3.;
              coefyr[g][i]=-tmp;
	    }else{
              if(!cal_tri)coefs[g][i]-=CalEdgeCurVacuum(i,1,g,difc[1],mcorrect)*surface_xy;
            };
	    //
            if(!redge_yr[x][y]){
	      if(y%2==0){id=meshid[z][y+1][x];}
	      else if(x!=0){id=meshid[z][y+1][x-1];};
              d2=mesh[id].GetDif(g,difc[1]);
	      // mcorrect
	      if(mcorrect){
		d2*=(1.-0.66666666*mcor[id][g]);
	      };
	      //
              dl=(d1[1]+d2)*half_pitch;
              if(!cal_tri)coefs[g][i]+=d1[1]*d2/dl*surface_xy;
            }else{
              if(!cal_tri)coefs[g][i]-=CalEdgeCurVacuum(i,1,g,difc[1],mcorrect)*surface_xy;
            };
	  };
	};
      };
    };
  };

};

real DHEXSystem::CalEdgeCurVacuum(int m,int dir,int g,int flag,bool mcorrect)
{
  real d=mesh[m].GetDif(g,flag);
  if(mcorrect){
    d*=(1.-0.666666666*mcor[m][g]);
  };
  return d/(2.1312*d+(mesh[m].GetLen(dir)*0.5))*-1.; 
};

void DHEXSystem::SetInitialFlux()
{
  SetInitialFlatFlux();
};

void DHEXSystem::SetDifc(string d1, string d2, string d3)
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

real DHEXSystem::CalFluxGeneral(int grp,real epsif,int iter)
{
  real err;

  if(!cal_tri){
    if(iter==0){
      err=CalFluxOmega(grp,0.1);
      //err=CalFlux(grp,iter,epsif);
      err=CalFluxOmega(grp,1e-3);
    }
    else{err=CalFlux(grp,iter,epsif);};
  }else{
    if(iter==0){
      err=CalFluxOmegaTriangle(grp,0.1);
      //err=CalFlux(grp,iter,epsif);
      err=CalFluxOmegaTriangle(grp,1e-3);
    }
    else{err=CalFlux(grp,iter,epsif);};
  };

  return err;
};

real DHEXSystem::CalFlux(int g,int iter,real epsif)
{
  int itmax=99;
  real *Src=new real[TotM];
  real *Fl =new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  bool InConv=false;
  if(cal_tri){
    InConv=SweepInnerIterationTriangle(itmax,epsif,g,omega,Fl,Src);
  }else{
    InConv=SweepInnerIteration(itmax,epsif,g,omega,Fl,Src);
  };

  real errf=RenewFluxInnerIteration(g,Fl);
  if(!InConv){
    bool InConv2=false;
    if(cal_tri){
      InConv2=SweepInnerIterationTriangle(itmax,epsif,g,1.0,Fl,Src);
    }else{
      InConv2=SweepInnerIteration(itmax,epsif,g,1.0,Fl,Src);
    };
    errf=RenewFluxInnerIteration(g,Fl);
    if(!InConv2)cout<<" Inner iteration is not converged in group "<<g<<"\n";
  };
  delete [] Fl;
  delete [] Src;
  return errf;
};

void DHEXSystem::PutTriangularMeshCorrection()
{
  cal_tri=true;

  real linv=1./(pitch*0.3333333333);
  real surf_tri=pitch*0.5/1.7320508*2.;

  tricor1.resize(TotM);
  tricor2.resize(TotM);
  for(int i=0;i<TotM;i++){
    tricor1[i].resize(grp);
    tricor2[i].resize(grp);
    real zlen=mesh[i].GetVolume()/hex_area;
    for(int g=0;g<grp;g++){
      real difc=mesh[i].GetMed()->GetMacxs().GetData1d(d).get_dat(g);
      tricor2[i][g]=difc*linv*surf_tri*zlen;
      real tmp=CalEdgeCurVacuum(i,0,g,0,false)*surf_tri*zlen;
      tricor1[i][g]=(2.1312*difc+pitch*0.5)/(2.1312*difc+pitch*0.1666666666)*tmp;

    };
  };

};

bool DHEXSystem::SweepInnerIteration(int itmax,real epsif,int g,real omega,real *Fl,real *Src)
{
  int ym=mi->GetYF();
  int zm=mi->GetZF();

  bool conv=true;
  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){

	bool ymod2eq0=true;
	if(y%2==1)ymod2eq0=false;
	int xleft_corr=-1;
	if(ymod2eq0)xleft_corr=0;

	int xl=edge[0][y][z][0];
	int xr=edge[0][y][z][1];
	for(int x=xl;x<=xr;x++){
	  int id=meshid[z][y][x];
	  if(id!=-1){
	    real flt=Src[id];
	    if(x!=xl)flt-=coefx[g][id]*Fl[id-1];
	    if(x!=xr)flt-=coefx[g][id+1]*Fl[id+1];
            if(!ledge_yl[x][y]){
	      int tmp=meshid[z][y-1][x+xleft_corr];
              flt-=coefyl[g][id]*Fl[tmp];
	    };
            if(!ledge_yr[x][y]){
              int tmp=meshid[z][y-1][x+xleft_corr+1];
	      flt-=coefyr[g][id]*Fl[tmp];
	    };
            if(!redge_yr[x][y]){
      	      int tmp=meshid[z][y+1][x+xleft_corr];
              flt-=coefyr[g][tmp]*Fl[tmp];
            };
	    if(!redge_yl[x][y]){
     	      int tmp=meshid[z][y+1][x+xleft_corr+1]; 
	      flt-=coefyl[g][tmp]*Fl[tmp];
	    };
	    if(z!=0){
              int tmp=meshid[z-1][y][x];
	      real fl=Fl[tmp];
	      if(cal_mcor){
	        fl*=(1.-mcor[tmp][g]);
	      };
              flt-=coefz[g][id]*fl;
	    };
	    if(z!=zm-1){
	      int tmp=meshid[z+1][y][x];
              real fl=Fl[tmp];
	      if(cal_mcor){
	      	fl*=(1.-mcor[tmp][g]);
	      };
              flt-=coefz[g][tmp]*fl;
	    };
	    real flo=Fl[id];
 	    real fln=flt/coefs[g][id];
  	    if(errmax<epsif){
 	      real err=fabs(fln/flo-1.);
	      if(err>errmax)errmax=err;
	    };
	    Fl[id]=flo+omega*(fln-flo);
	  };
	};
      };
    };
    if(errmax<epsif&&iter>0){
      //cout<<"  (group:"<<g<<") "<<iter<<"  ("<<errmax<<")\n";
      break;
    };
    if(iter==itmax-1&&itmax!=1){
      //cout<<"Inner iteration is exceeded for "<<itmax;
      //cout<<"  (group:"<<g<<")(errmax:"<<errmax<<")\n";
      conv=false;
    };
  };

  return conv;
};

bool DHEXSystem::SweepInnerIterationTriangle(int itmax,real epsif,int g,real omega,real *Fl,real *Src)
{
  vector< vector<real> > flx_tri(TotM);
  for(int i=0;i<TotM;i++){
    flx_tri[i].resize(6);
    real tmp=Fl[i];
    for(int j=0;j<6;j++){
      flx_tri[i][j]=tmp;
    };
  };

  for(int i=0;i<TotM;i++){
    Src[i]*=0.1666666666666;
  };

  int ym=mi->GetYF();
  int zm=mi->GetZF();

  real *a=new real[6*6];
  real *b=new real[6];

  bool conv=true;
  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){

	bool ymod2eq0=true;
	if(y%2==1)ymod2eq0=false;
	int xleft_corr=-1;
	if(ymod2eq0)xleft_corr=0;

	int xl=edge[0][y][z][0];
	int xr=edge[0][y][z][1];
	for(int x=xl;x<=xr;x++){
	  int id=meshid[z][y][x];
	  if(id!=-1){

	    for(int i=0;i<6*6;i++){a[i]=0.;};
            real tmp=coefs[g][id];
            for(int i=0;i<6;i++){
	      a[i*6+i]=tmp;
	    };
	    real flt=Src[id];
	    for(int i=0;i<6;i++){b[i]=flt;};

	    if(x!=xl){
              real tmp=coefx[g][id];
              b[4]-=tmp*flx_tri[id-1][1];
              a[4*6+4]-=tmp;
	    }else{
	      a[4*6+4]-=tricor1[id][g];
	    };
	    if(x!=xr){
              real tmp=coefx[g][id+1];
              b[1]-=tmp*flx_tri[id+1][4];
              a[1*6+1]-=tmp;
	    }else{
	      a[1*6+1]-=tricor1[id][g];
	    };
            if(!ledge_yl[x][y]){
	      int tmp=meshid[z][y-1][x+xleft_corr];
              real tmp2=coefyl[g][id];
              b[5]-=tmp2*flx_tri[tmp][2];
              a[5*6+5]-=tmp2;
	    }else{
	      a[5*6+5]-=tricor1[id][g];
	    };
            if(!ledge_yr[x][y]){
              int tmp=meshid[z][y-1][x+xleft_corr+1];
              real tmp2=coefyr[g][id];
	      b[0]-=tmp2*flx_tri[tmp][3];
              a[0]-=tmp2;
	    }else{
	      a[0*6+0]-=tricor1[id][g];
	    };
            if(!redge_yr[x][y]){
      	      int tmp=meshid[z][y+1][x+xleft_corr];
              real tmp2=coefyr[g][tmp];
              b[3]-=tmp2*flx_tri[tmp][0];
              a[3*6+3]-=tmp2;
	    }else{
	      a[3*6+3]-=tricor1[id][g];
            };
	    if(!redge_yl[x][y]){
     	      int tmp=meshid[z][y+1][x+xleft_corr+1]; 
              real tmp2=coefyl[g][tmp];
	      b[2]-=tmp2*flx_tri[tmp][5];
              a[2*6+2]-=tmp2;
	    }else{
	      a[2*6+2]-=tricor1[id][g];
	    };
	    if(z!=0){
              int tmp=meshid[z-1][y][x];
              real tmp2=coefz[g][id];
	      for(int s=0;s<6;s++){
                b[s]-=tmp2*flx_tri[tmp][s];
	      };
	    };
	    if(z!=zm-1){
	      int tmp=meshid[z+1][y][x];
              real tmp2=coefz[g][tmp];
	      for(int s=0;s<6;s++){
		b[s]-=tmp2*flx_tri[tmp][s];
	      };
	    };

            real coef=tricor2[id][g];
            for(int s=0;s<6;s++){
	      a[s*6+s]+=coef*2.;
	    };
	    a[0*6+1]-=coef;
	    a[0*6+5]-=coef;
	    a[1*6+0]-=coef;
	    a[1*6+2]-=coef;
	    a[2*6+1]-=coef;
	    a[2*6+3]-=coef;
	    a[3*6+2]-=coef;
	    a[3*6+4]-=coef;
	    a[4*6+3]-=coef;
	    a[4*6+5]-=coef;
	    a[5*6+0]-=coef;
	    a[5*6+4]-=coef;
            gauss(a,b,6,1);

            for(int s=0;s<6;s++){
	      real flo=flx_tri[id][s];
	      real fln=b[s];
	      if(errmax<epsif){
		real err=fabs(fln/flo-1.);
		if(err>errmax)errmax=err;
	      };
	      flx_tri[id][s]=flo+omega*(fln-flo);
	    };

	  };
	};
      };
    };
    if(errmax<epsif&&iter>0){
      //cout<<"  (group:"<<g<<") "<<iter<<"  ("<<errmax<<")\n";
      break;
    };
    if(iter==itmax-1&&itmax!=1){
      //cout<<"Inner iteration is exceeded for "<<itmax;
      //cout<<"  (group:"<<g<<")(errmax:"<<errmax<<")\n";
      conv=false;
    };
  };

  delete [] a;
  delete [] b;

  for(int i=0;i<TotM;i++){
    real sum=0.;
    for(int j=0;j<6;j++){
      sum+=flx_tri[i][j];
    };
    Fl[i]=sum*0.16666666666666;
  };

  return conv;
};

real DHEXSystem::CalFluxOmega(int g,real epsif)
{
  int itmax=999;
  real *Src= new real[TotM];
  real *Fl = new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  SweepInnerIteration(itmax,epsif,g,omega,Fl,Src);
  // *** for FluxOmega ***
  real *Res= new real[TotM]; 
  real tt1, tt2; 
  omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  //real omega=1.;
  // ********************

  int ym=mi->GetYF();
  int zm=mi->GetZF();

  itmax=2;

  for(int iter=0;iter<itmax;iter++){

    real errmax=0.;
    tt1=0.;    //fluxomega
    int cnt=0; //fluxomega

    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){

	bool ymod2eq0=true;
	if(y%2==1)ymod2eq0=false;
	int xleft_corr=-1;
	if(ymod2eq0)xleft_corr=0;

	int xl=edge[0][y][z][0];
	int xr=edge[0][y][z][1];
	for(int x=xl;x<=xr;x++){
	  int id=meshid[z][y][x];
	  if(id!=-1){
	    real flt=Src[id];
	    if(x!=xl)flt-=coefx[g][id]*Fl[id-1];
	    if(x!=xr)flt-=coefx[g][id+1]*Fl[id+1];
            if(!ledge_yl[x][y]){
	      int tmp=meshid[z][y-1][x+xleft_corr];
              flt-=coefyl[g][id]*Fl[tmp];
	    };
            if(!ledge_yr[x][y]){
              int tmp=meshid[z][y-1][x+xleft_corr+1];
	      flt-=coefyr[g][id]*Fl[tmp];
	    };
            if(!redge_yr[x][y]){
      	      int tmp=meshid[z][y+1][x+xleft_corr];
              flt-=coefyr[g][tmp]*Fl[tmp];
            };
	    if(!redge_yl[x][y]){
     	      int tmp=meshid[z][y+1][x+xleft_corr+1]; 
	      flt-=coefyl[g][tmp]*Fl[tmp];
	    };
	    if(z!=0){
              int tmp=meshid[z-1][y][x];
              real fl=Fl[tmp];
	      if(cal_mcor){
	        fl*=(1.-mcor[tmp][g]);
	      };
              flt-=coefz[g][id]*fl;
	    };
	    if(z!=zm-1){
	      int tmp=meshid[z+1][y][x];
              real fl=Fl[tmp];
	      if(cal_mcor){
	      	fl*=(1.-mcor[tmp][g]);
	      };
              flt-=coefz[g][tmp]*fl;
	    };
	    real flo=Fl[id];
 	    real fln=flt/coefs[g][id];
  	    if(errmax<epsif){
 	      real err=fabs(fln/flo-1.);
	      if(err>errmax)errmax=err;
	    };
	    Fl[id]=flo+omega*(fln-flo);
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
	  };
	};
      };
    };
    tt1/=cnt; //fluxomega
    if(cnt==0)tt1=0.; // fluxomega
  };

  real errf=RenewFluxInnerIteration(g,Fl);

  // ** Fluxomega **
  if(tt1!=0.){
    real rnew=(omega-1+tt1)/omega/sqrt(tt1);
    rnew=rnew*rnew;
    if(rnew>1.){
      omega=1.0;
    }else{
      omega=2.0/(1+sqrt(1-rnew)); 
    };
  }else{
    omega=1.;
  };
  opt.PutOmegai(g,omega);
  cout<<"   ** omega for inner iteration : "<<omega<<" (group:"<<g<<")\n";
  // ***************

  delete [] Fl;
  delete [] Src;
  delete [] Res; 
  return errf;
};

real DHEXSystem::CalFluxOmegaTriangle(int g,real epsif)
{
  int itmax=999;
  real *Src= new real[TotM];
  real *Fl = new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  SweepInnerIterationTriangle(itmax,epsif,g,omega,Fl,Src);
  // *** for FluxOmega ***
  real *Res= new real[TotM*6]; 
  real tt1, tt2; 
  omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  //real omega=1.;
  // ********************

  vector< vector<real> > flx_tri(TotM);
  for(int i=0;i<TotM;i++){
    flx_tri[i].resize(6);
    real tmp=Fl[i];
    for(int j=0;j<6;j++){
      flx_tri[i][j]=tmp;
    };
  };

  int ym=mi->GetYF();
  int zm=mi->GetZF();

  real *a=new real[6*6];
  real *b=new real[6];

  itmax=2;

  for(int iter=0;iter<itmax;iter++){

    real errmax=0.;
    tt1=0.;
    int cnt=0;

    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){

	bool ymod2eq0=true;
	if(y%2==1)ymod2eq0=false;
	int xleft_corr=-1;
	if(ymod2eq0)xleft_corr=0;

	int xl=edge[0][y][z][0];
	int xr=edge[0][y][z][1];
	for(int x=xl;x<=xr;x++){
	  int id=meshid[z][y][x];
	  if(id!=-1){

	    for(int i=0;i<6*6;i++){a[i]=0.;};
            real tmp=coefs[g][id];
            for(int i=0;i<6;i++){
	      a[i*6+i]=tmp;
	    };
	    real flt=Src[id];
	    for(int i=0;i<6;i++){b[i]=flt;};

	    if(x!=xl){
              real tmp=coefx[g][id];
              b[4]-=tmp*flx_tri[id-1][1];
              a[4*6+4]-=tmp;
	    }else{
	      a[4*6+4]-=tricor1[id][g];
	    };
	    if(x!=xr){
              real tmp=coefx[g][id+1];
              b[1]-=tmp*flx_tri[id+1][4];
              a[1*6+1]-=tmp;
	    }else{
	      a[1*6+1]-=tricor1[id][g];
	    };
            if(!ledge_yl[x][y]){
	      int tmp=meshid[z][y-1][x+xleft_corr];
              real tmp2=coefyl[g][id];
              b[5]-=tmp2*flx_tri[tmp][2];
              a[5*6+5]-=tmp2;
	    }else{
	      a[5*6+5]-=tricor1[id][g];
	    };
            if(!ledge_yr[x][y]){
              int tmp=meshid[z][y-1][x+xleft_corr+1];
              real tmp2=coefyr[g][id];
	      b[0]-=tmp2*flx_tri[tmp][3];
              a[0]-=tmp2;
	    }else{
	      a[0*6+0]-=tricor1[id][g];
	    };
            if(!redge_yr[x][y]){
      	      int tmp=meshid[z][y+1][x+xleft_corr];
              real tmp2=coefyr[g][tmp];
              b[3]-=tmp2*flx_tri[tmp][0];
              a[3*6+3]-=tmp2;
	    }else{
	      a[3*6+3]-=tricor1[id][g];
            };
	    if(!redge_yl[x][y]){
     	      int tmp=meshid[z][y+1][x+xleft_corr+1]; 
              real tmp2=coefyl[g][tmp];
	      b[2]-=tmp2*flx_tri[tmp][5];
              a[2*6+2]-=tmp2;
	    }else{
	      a[2*6+2]-=tricor1[id][g];
	    };
	    if(z!=0){
              int tmp=meshid[z-1][y][x];
              real tmp2=coefz[g][id];
	      for(int s=0;s<6;s++){
                b[s]-=tmp2*flx_tri[tmp][s];
	      };
	    };
	    if(z!=zm-1){
	      int tmp=meshid[z+1][y][x];
              real tmp2=coefz[g][tmp];
	      for(int s=0;s<6;s++){
		b[s]-=tmp2*flx_tri[tmp][s];
	      };
	    };

            real coef=tricor2[id][g];
            for(int s=0;s<6;s++){
	      a[s*6+s]+=coef*2.;
	    };
	    a[0*6+1]-=coef;
	    a[0*6+5]-=coef;
	    a[1*6+0]-=coef;
	    a[1*6+2]-=coef;
	    a[2*6+1]-=coef;
	    a[2*6+3]-=coef;
	    a[3*6+2]-=coef;
	    a[3*6+4]-=coef;
	    a[4*6+3]-=coef;
	    a[4*6+5]-=coef;
	    a[5*6+0]-=coef;
	    a[5*6+4]-=coef;
            gauss(a,b,6,1);

            for(int s=0;s<6;s++){
	      real flo=flx_tri[id][s];
	      real fln=b[s];
	      if(errmax<epsif){
		real err=fabs(fln/flo-1.);
		if(err>errmax)errmax=err;
	      };
	      flx_tri[id][s]=flo+omega*(fln-flo);
	      // +++
	      real reso=Res[id*6+s];
	      real resn=pow(fln-flo,2);
	      Res[id*6+s]=resn;
	      tt2=sqrt(resn/reso);
	      if(tt2>0&&tt2<1){
		tt1+=tt2;
		cnt++;
	      };
	      // +++
	    };

	  };
	};
      };
    };
    tt1/=cnt;
    if(cnt==0)tt1=0.;
  };

  delete [] a;
  delete [] b;

  for(int i=0;i<TotM;i++){
    real sum=0.;
    for(int j=0;j<6;j++){
      sum+=flx_tri[i][j];
    };
    Fl[i]=sum*0.16666666666666;
  };

  real errf=RenewFluxInnerIteration(g,Fl);

  // ** Fluxomega **
  if(tt1!=0.){
    real rnew=(omega-1+tt1)/omega/sqrt(tt1);
    rnew=rnew*rnew;
    if(rnew>1.){
      omega=1.0;
    }else{
      omega=2.0/(1+sqrt(1-rnew)); 
    };
  }else{
    omega=1.;
  };
  opt.PutOmegai(g,omega);
  cout<<"   ** omega for inner iteration : "<<omega<<" (group:"<<g<<")\n";
  // ***************

  delete [] Fl;
  delete [] Src;
  delete [] Res; 

  return errf;
};

real DHEXSystem::RenewFluxInnerIteration(int g,real *Fl)
{
  real errf=0.;
  for(int i=0;i<TotM;i++){
    real tmp=Fl[i];
    real err=fabs(tmp/mesh[i].GetFlux().get_dat(g)-1.0);
    if(err>errf)errf=err;
    mesh[i].GetFlux().put_data(g,tmp);
  };
  return errf;
};

void DHEXSystem::GetFluxAndSrcInnerIteration(int g,real *Fl,real *Src)
{
  for(int i=0;i<TotM;i++){
    Fl[i] =mesh[i].GetFlux().get_dat(g);
    Src[i]=mesh[i].GetSrcin()*PI4;
  };
};

real DHEXSystem::CalCurrentX(int x,int y,int z,int g)
{
  bool redge=false;
  if(x-1==edge[0][y][z][1])redge=true;

  if(!redge){
    // calculation for left-side current
    bool ledge=false;
    if(x==edge[0][y][z][0])ledge=true;
    int id=meshid[z][y][x];
    real fl=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      fl/=(1.-mcor[id][g]);
    };
    real coefl=coefx[g][id];
    if(!ledge){
      real tmp=mesh[id-1].GetFlux().get_dat(g);
      if(cal_mcor){
	tmp/=(1.-mcor[id-1][g]);
      };
      return coefl*(fl-tmp);
    }else{
      real zlen=mi->GetFMeshL(2,z);
      real tmp=CalEdgeCurVacuum(id,0,g,difc[0],cal_mcor);
      return tmp*fl*zlen*len_plane;
    };
  }else{
    // calculation for right-edge current
    int id=meshid[z][y][x-1];
    real zlen=mi->GetFMeshL(2,z);
    real flx=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      flx/=(1.-mcor[id][g]);
    };
    real tmp=flx*zlen*len_plane;
    return -1.*CalEdgeCurVacuum(id,0,g,difc[0],cal_mcor)*tmp;
  };
};

real DHEXSystem::CalCurrentZ(int x,int y,int z,int g)
{
  bool redge=false;
  if(z==mi->GetZF())redge=true;

  if(!redge){
    // calculation for left-side current
    int id=meshid[z][y][x];
    real fl=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      //fl/=(1.-mcor[id][g]);
    };
    real coefl=coefz[g][id];
    if(z!=0){
      int id2=meshid[z-1][y][x];
      real tmp=mesh[id2].GetFlux().get_dat(g);
      if(cal_mcor){
	//tmp/=(1.-mcor[id2][g]);
      };
      return coefl*(fl-tmp);
    }else{
      if(BC[4]==1){
	return 0.;
      }else{
        //return CalEdgeCurVacuum(id,2,g,difc[2],cal_mcor)*fl*hex_area;
        return CalEdgeCurVacuum(id,2,g,difc[2])*fl*hex_area;
      };
    };
  }else{
    // calculation for right-edge current
    int id=meshid[z-1][y][x];
    real tmp=mesh[id].GetFlux().get_dat(g)*hex_area;
    if(cal_mcor){
      //tmp/=(1.-mcor[id][g]);
    };
    if(BC[5]==1){
      return 0.;
    }else{
      //return -1.*CalEdgeCurVacuum(id,2,g,difc[2],cal_mcor)*tmp;
      return -1.*CalEdgeCurVacuum(id,2,g,difc[2])*tmp;
    };
  };
};

real DHEXSystem::CalCurrentYL(int x,int y,int z,int g)
{
  bool redge=false;
  if(y==mi->GetYF()){
    redge=true;
  }else if(y%2==1&&x==mi->GetXF()){
    redge=true;
  }else{
    if(meshid[z][y][x]==-1)redge=true;
  };

  if(!redge){
    // calculation for left-side current
    int id=meshid[z][y][x];
    real fl=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      fl/=(1.-mcor[id][g]);
    };
    real coefl=coefyl[g][id];
    if(ledge_yl[x][y]){
      real zlen=mi->GetFMeshL(2,z);
      return CalEdgeCurVacuum(id,1,g,difc[1],cal_mcor)*fl*zlen*len_plane;
    }else{
      int xpos;
      if(y%2==0){
	xpos=x;
      }else{
	xpos=x-1;
      };
      int id2=meshid[z][y-1][xpos];
      real tmp=mesh[id2].GetFlux().get_dat(g);
      if(cal_mcor){
	tmp/=(1.-mcor[id2][g]);
      };
      return coefl*(fl-tmp);
    };
  }else{
    // calculation for right-edge current
    int xnew=x;
    if(y%2==1)xnew--;
    int id=meshid[z][y-1][xnew];
    real zlen=mi->GetFMeshL(2,z);
    real flx=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      flx/=(1.-mcor[id][g]);
    };
    real tmp=flx*zlen*len_plane;
    return -1.*CalEdgeCurVacuum(id,1,g,difc[1],cal_mcor)*tmp;
  };
};

real DHEXSystem::CalCurrentYR(int x,int y,int z,int g)
{
  bool redge=false;
  if(y==mi->GetYF()){
    redge=true;
  }else if(y%2==0&&x==-1){
    redge=true;
  }else{
    if(meshid[z][y][x]==-1)redge=true;
  };

  if(!redge){
    // calculation for left-side current
    int id=meshid[z][y][x];
    real fl=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      fl/=(1.-mcor[id][g]);
    };
    real coefl=coefyr[g][id];
    if(ledge_yr[x][y]){
      real zlen=mi->GetFMeshL(2,z);
      return CalEdgeCurVacuum(id,1,g,difc[1],cal_mcor)*fl*zlen*len_plane;
    }else{
      int xpos;
      if(y%2==0){
	xpos=x+1;
      }else{
	xpos=x;
      };
      int id2=meshid[z][y-1][xpos];
      real tmp=mesh[id2].GetFlux().get_dat(g);
      if(cal_mcor){
	tmp/=(1.-mcor[id2][g]);
      };
      return coefl*(fl-tmp);
    };
  }else{
    // calculation for right-edge current
    int xnew=x;
    if(y%2==0)xnew++;
    int id=meshid[z][y-1][xnew];
    real zlen=mi->GetFMeshL(2,z);
    real flx=mesh[id].GetFlux().get_dat(g);
    if(cal_mcor){
      flx/=(1.-mcor[id][g]);
    };
    real tmp=flx*zlen*len_plane;
    return -1.*CalEdgeCurVacuum(id,1,g,difc[1],cal_mcor)*tmp;
  };
};

// **********************************************
// * For CMFD acceleration                      *
// **********************************************

void DHEXSystem::DoAcceleration(int iter,real errs,real fiss)
{
  if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+0.5);
  //if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+1e9);
};

void DHEXSystem::DoCMFDAcceleration(real delk)
{
  int ngrp=1;

  int xr=mi->GetXC();
  int yr=mi->GetYC();
  int zr=mi->GetZC();

  // +++ Calculation for fine-mesh current

  vector< vector< vector<real> > > fcur_x; // [z+1][y+1][x+1]
  vector< vector< vector<real> > > fcur_yl;
  vector< vector< vector<real> > > fcur_yr;
  vector< vector< vector<real> > > fcur_z;
  vector< vector< vector<real> > > moddif_x;
  vector< vector< vector<real> > > moddif_yl;
  vector< vector< vector<real> > > moddif_yr;
  vector< vector< vector<real> > > moddif_z;
  fcur_x.resize(zr+1);
  fcur_yl.resize(zr+1);
  fcur_yr.resize(zr+1);
  fcur_z.resize(zr+1);
  moddif_x.resize(zr+1);
  moddif_yl.resize(zr+1);
  moddif_yr.resize(zr+1);
  moddif_z.resize(zr+1);
  for(int i=0;i<zr+1;i++){
    fcur_x[i].resize(yr+1);
    fcur_yl[i].resize(yr+1);
    fcur_yr[i].resize(yr+1);
    fcur_z[i].resize(yr+1);
    moddif_x[i].resize(yr+1);
    moddif_yl[i].resize(yr+1);
    moddif_yr[i].resize(yr+1);
    moddif_z[i].resize(yr+1);
    for(int j=0;j<yr+1;j++){
      fcur_x[i][j].resize(xr+1,0.);
      fcur_yl[i][j].resize(xr+1,0.);
      fcur_yr[i][j].resize(xr+1,0.);
      fcur_z[i][j].resize(xr+1,0.);
      moddif_x[i][j].resize(xr+1,0.);
      moddif_yl[i][j].resize(xr+1,0.);
      moddif_yr[i][j].resize(xr+1,0.);
      moddif_z[i][j].resize(xr+1,0.);
    };
  };

  for(int i=0;i<grp;i++){
    CalCoarseCur(i,fcur_x,fcur_yl,fcur_yr,fcur_z,true);
  };

  // +++ Calculation for coarse-group cross sections

  int xyzr=xr*yr*zr;
  DHEXSystem cm(Ndim,ngrp,xyzr);
  if(!print)cm.NoPrint();

  vector<real> cc_nusigf(CTotM);
  vector<real> cc_sigt(CTotM);
  vector<real> cc_d(CTotM);
  vector<real> cc_flx(CTotM);
  vector<real> cc_sigs(CTotM);
  vector<real> cc_chi(CTotM);
  for(int i=0;i<CTotM;i++){
    cc_nusigf[i]=0.;
    cc_sigt[i]=0.;
    cc_d[i]=0.;
    cc_flx[i]=0.;
    cc_sigs[i]=0.;
    cc_chi[i]=0.;
  };

  GroupData1D volflx(grp);
  int iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi->GetCMeshF(2,z1);z2++){
      for(int y1=0;y1<yr;y1++){
        for(int x1=0;x1<xr;x1++){
          int id=meshid[iz][y1][x1];
          if(id!=-1){
            int index=cmeshid[x1][y1][z1];
	    real vol=mesh[id].GetVolume();
   	    volflx=mesh[id].GetFlux()*vol;
   	    cc_d[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetD();
   	    cc_sigt[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetSigt();
   	    if(opt.Forward()){
              //cc_sigs[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetSigs().get_sumx();
              cc_sigs[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetData1d(sigst);
              cc_nusigf[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetNusigf();
	    }else{
   	      real flx_chi=volflx*mesh[id].GetMed()->GetMacxs().GetData1d(chi);
   	      cc_chi[index]+=flx_chi;
   	      cc_nusigf[index]+=flx_chi*mesh[id].GetMed()->GetMacxs().GetData1d(nusigf).get_sum();
  	      for(int i=0;i<grp;i++){
  	        real tmp=volflx.get_dat(i);
   	        for(int j=0;j<grp;j++){
	  	  cc_sigs[index]+=tmp*mesh[id].GetMed()->GetMacxs().GetSigs().get_dat(j,i);
		};
    	      };
            };
   	    //cc_flx[index]+=mesh[id].GetVolumeFlux();
   	    cc_flx[index]+=volflx.get_sum();
          };
        };
      };
      iz++;
    };
  };

  Medium minp(ngrp);
  minp.PutPL(0);
  real inv_delk=1./delk;
  for(int i=0;i<CTotM;i++){
    real tmp=1./cc_flx[i];
    real fis=cc_nusigf[i]*tmp;
    if(opt.Forward()){
      minp.GetMacxs().GetNusigf().put_data(0,fis);
      minp.GetMacxs().GetKai().put_data(0,1.);
    }else{
      real ttt=0.;
      if(cc_chi[i]!=0.)ttt=cc_nusigf[i]/cc_chi[i];
      minp.GetMacxs().GetNusigf().put_data(0,ttt);
      minp.GetMacxs().GetKai().put_data(0,cc_chi[i]*tmp);
    };
    real sigtfic=cc_sigt[i]*tmp-fis*inv_delk;
    minp.GetMacxs().GetSigt().put_data(0,sigtfic);
    minp.GetMacxs().GetSigs().put_data(0,0,cc_sigs[i]*tmp);
    minp.GetMacxs().GetD().put_data(0,cc_d[i]*tmp*cmfd_factor);
    cm.AddMedium(minp);
  };

  // ***

  vector<real> xwid(xr);
  vector<real> ywid(yr);
  vector<real> zwid(zr);
  vector<int> fmx(xr);
  vector<int> fmy(yr);
  vector<int> fmz(zr);
  for(int i=0;i<xr;i++){
    xwid[i]=mi->GetCMeshL(0,i);
    fmx[i]=1;
  };
  for(int i=0;i<yr;i++){
    ywid[i]=mi->GetCMeshL(1,i);
    fmy[i]=1;
  };
  for(int i=0;i<zr;i++){
    zwid[i]=mi->GetCMeshL(2,i);
    fmz[i]=1;
  };
  vector<int> asmmap(xyzr);
  int index=0;
  for(int z1=0;z1<zr;z1++){
    for(int y1=0;y1<yr;y1++){
      for(int x1=0;x1<xr;x1++){
	asmmap[index]=cmeshid[x1][y1][z1];
	index++;
      };
    };
  };
  CartMeshInfo cmi;
  cmi.PutMeshInfo(xr,yr,zr,fmx,fmy,fmz,xwid,ywid,zwid,asmmap);
  int bci[6];
  for(int i=0;i<6;i++){
    bci[i]=mi->GetBC(i);
  };
  cmi.PutBoundaryCondition(bci);

  cm.PutCartMeshInfo(&cmi,pitch);

  // Calculation of initial flux for CMFD 
  vector< vector<real> > FlxF(ngrp);  //[g][m]
  for(int i=0;i<ngrp;i++){
    FlxF[i].resize(CTotM,0.);
  };
  for(int i=0;i<CTotM;i++){
    real inv_vol=1./cm.GetMesh(i).GetVolume();
    for(int g=0;g<ngrp;g++){
      FlxF[g][i]=cc_flx[i]*inv_vol;
    };
  };

  GeneralOption option;
  real cm_epsf=opt.GetEpsf()*0.5;
  if(cm_epsf>1e-4)cm_epsf=1e-4;
  if(!opt.Forward())option.PutAdjointCal();
  option.PutEpsf(1e-5);
  option.PutEpsk(1e-6);
  option.PutEpss(1e-6);

  cm.PutGeneralOption(option);

  cm.CalCoef();

  cm.PreIgenCoarse(fcur_x,fcur_yl,fcur_yr,fcur_z,
                   moddif_x,moddif_yl,moddif_yr,moddif_z,FlxF);
  // Iteration in CMFD
  real k2=cm.CalIgenCoarse(moddif_x,moddif_yl,moddif_yr,moddif_z);
  //real k2=cm.CalIgen();

  real ke=1./k2+1/delk;
  ke=1./ke;
  ke=ke/k2;

  for(int i=0;i<CTotM;i++){
    real vol=cm.GetMesh(i).GetVolume();
    real tmp=vol*ke;
    for(int g=0;g<ngrp;g++){ // ngrp=1
      cc_flx[i]=cm.GetMesh(i).GetFlux().get_dat(g)*tmp/cc_flx[i];
    };
  };

  iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi->GetCMeshF(2,z1);z2++){
      int iy=0;
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<mi->GetCMeshF(1,y1);y2++){
          int ix=0;
          for(int x1=0;x1<xr;x1++){
  	    int index=cmeshid[x1][y1][z1];
	    for(int x2=0;x2<mi->GetCMeshF(0,x1);x2++){
  	      int id=meshid[iz][iy][ix];
	      if(id!=-1){
		real tmp=cc_flx[index];
		for(int g=0;g<grp;g++){
    	          mesh[id].GetFlux().put_data(g,mesh[id].GetFlux().get_dat(g)*tmp);
		};
	      };
	      ix++;
	    };
          };
	  iy++;
        };
      };
      iz++;
    };
  };
};

real DHEXSystem::CalIgenCoarse(vector< vector< vector<real> > > &moddif_x,
                     vector< vector< vector<real> > > &moddif_yl,
                     vector< vector< vector<real> > > &moddif_yr,
                     vector< vector< vector<real> > > &moddif_z)
{
  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  //SetInitialFlux();  // o
  NormalizeFissionSrc();

  real cin=0.1;
  fiss=1.;  

  bool conv=false;
  for(int iter=0;iter<999;iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    real tmp=5e-7;
    if(opt.GetEpsf()<5e-7)tmp=opt.GetEpsf()*0.5;
    if(cin<tmp)cin=tmp;

    errf=0.;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,fiss,0);

      if(iter==0&&g==0){
	CalFluxOmega(ginp,0.1);
	CalFluxOmega(ginp,1e-3);
      };
      real err=CalFluxModifiedLeakage(ginp,cin,moddif_x,moddif_yl,moddif_yr,moddif_z); // o
      //real err=CalFlux(ginp,0,cin); // o
      if(err>errf)errf=err;
      SetZeroScatSrc(ginp);
      AddDownScatSrc(ginp,0);
    };
	
    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real tmp2=mesh[i].GetFissionSrc();
	real err=fabs(tmp2/tmp-1.0);
	if(err>errs)errs=err;
      };
    };

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    //WriteIterationInfo(iter,fiss,errk,errf,errs); // o

    if(conv){
      if(opt.Converged(errf,errk,errs))break;
    };
    if(opt.Converged(errf,errk,errs))conv=true;
  };

  if(!conv){
    cout<<"Coarse group eigenvalue calulation is not converged.\n";
  };

  return fiss;
}

real DHEXSystem::CalFluxModifiedLeakage(int g,real epsif,
    vector< vector< vector<real> > > &moddif_x,
    vector< vector< vector<real> > > &moddif_yl,
    vector< vector< vector<real> > > &moddif_yr,
    vector< vector< vector<real> > > &moddif_z)
{
  //real omega=1.;  // fluxmodifiedleakage
  real omega=opt.GetOmegai(g);

  int itmax=999;
  real *Src=new real[TotM];
  real *Fl =new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  int ym=mi->GetYF();
  int zm=mi->GetZF();

  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){

	bool ymod2eq0=true;
	if(y%2==1)ymod2eq0=false;
	int xleft_corr=-1;
	if(ymod2eq0)xleft_corr=0;

	int xl=edge[0][y][z][0];
	int xr=edge[0][y][z][1];
	for(int x=xl;x<=xr;x++){
	  int id=meshid[z][y][x];
	  if(id!=-1){
	    real flt=Src[id];
	    //if(x!=xl)flt-=(coefx[g][id]-moddif_x[z][y][x])*Fl[id-1];
	    if(x!=xl)flt-=coefx[g][id]*Fl[id-1];
	    if(x!=xr)flt-=(coefx[g][id+1]+moddif_x[z][y][x+1])*Fl[id+1];
            if(!ledge_yl[x][y]){
	      int tmp=meshid[z][y-1][x+xleft_corr];
              //flt-=(coefyl[g][id]-moddif_yl[z][y][x])*Fl[tmp];
              flt-=coefyl[g][id]*Fl[tmp];
	    };
            if(!ledge_yr[x][y]){
              int tmp=meshid[z][y-1][x+xleft_corr+1];
	      //flt-=(coefyr[g][id]-moddif_yr[z][y][x])*Fl[tmp];
	      flt-=coefyr[g][id]*Fl[tmp];
	    };
            if(!redge_yr[x][y]){
      	      int tmp=meshid[z][y+1][x+xleft_corr];
              flt-=(coefyr[g][tmp]+moddif_yr[z][y+1][x+xleft_corr])*Fl[tmp];
            };
	    if(!redge_yl[x][y]){
     	      int tmp=meshid[z][y+1][x+xleft_corr+1]; 
	      flt-=(coefyl[g][tmp]+moddif_yl[z][y+1][x+xleft_corr+1])*Fl[tmp];
	    };
            //if(z!=0)flt-=(coefz[g][id]-moddif_z[z][y][x])*Fl[meshid[z-1][y][x]];
            if(z!=0)flt-=coefz[g][id]*Fl[meshid[z-1][y][x]];
	    if(z!=zm-1){
	      int tmp=meshid[z+1][y][x];
              flt-=(coefz[g][tmp]+moddif_z[z+1][y][x])*Fl[tmp];
	    };
	    real flo=Fl[id];
 	    real fln=flt/coefs[g][id];
  	    if(errmax<epsif){
 	      real err=fabs(fln/flo-1.);
	      if(err>errmax)errmax=err;
	    };
	    Fl[id]=flo+omega*(fln-flo);
	  };
	};
      };
    };
    if(errmax<epsif&&iter>0){
      //cout<<"  (group:"<<g<<") "<<iter<<"  ("<<errmax<<")\n";
      break;
    };
    if(iter==itmax-1&&itmax!=1){
      //cout<<"Inner iteration is exceeded for "<<itmax;
      //cout<<"  (group:"<<g<<")(errmax:"<<errmax<<")\n";
    };
  };

  real errf=RenewFluxInnerIteration(g,Fl);

  delete [] Fl;
  delete [] Src;
  return errf;
};



void DHEXSystem::PreIgenCoarse
(vector< vector< vector<real> > > &fcur_x,
 vector< vector< vector<real> > > &fcur_yl,
 vector< vector< vector<real> > > &fcur_yr,
 vector< vector< vector<real> > > &fcur_z,
 vector< vector< vector<real> > > &moddif_x,
 vector< vector< vector<real> > > &moddif_yl,
 vector< vector< vector<real> > > &moddif_yr,
 vector< vector< vector<real> > > &moddif_z,
 vector< vector<real> > &flxf)
{
  for(int i=0;i<TotM;i++){
    for(int g=0;g<grp;g++){
      mesh[i].GetFlux().put_data(g,flxf[g][i]);
    };
  };
  for(int g=0;g<grp;g++){
    CalCoarseCur(g,moddif_x,moddif_yl,moddif_yr,moddif_z,false);
    // ModD <- Current in coarse mesh
  };

  CalCorDif(moddif_x,moddif_yl,moddif_yr,moddif_z,
            fcur_x,fcur_yl,fcur_yr,fcur_z); // ModD <- Diffusion coefficient 
  ModCoef(moddif_x,moddif_yl,moddif_yr,moddif_z);
};

void DHEXSystem::CalCorDif
(vector< vector< vector<real> > > &ccur_x,
 vector< vector< vector<real> > > &ccur_yl,
 vector< vector< vector<real> > > &ccur_yr,
 vector< vector< vector<real> > > &ccur_z,
 vector< vector< vector<real> > > &fcur_x,
 vector< vector< vector<real> > > &fcur_yl,
 vector< vector< vector<real> > > &fcur_yr,
 vector< vector< vector<real> > > &fcur_z)
{
  // !! caution !!
  // This method should be only used for coarse-mesh calculation
  // (Coarse-mesh = Fine-Mesh)

  int xr=mi->GetXC();
  int yr=mi->GetYC();
  int zr=mi->GetZC();

  for(int z=0;z<zr;z++){
    for(int y=0;y<yr;y++){
      int corr=y%2;
      for(int x=0;x<xr;x++){
        int id=meshid[z][y][x];
	if(id!=-1){
          real fl=mesh[id].GetFlux().get_dat(0);
	  real fl_inv=1./fl;
          // X-direction
          real tmp=0.;
	  if(x!=edge[0][y][z][0]){
	    tmp=mesh[id-1].GetFlux().get_dat(0);
	  };
          ccur_x[z][y][x]=(fcur_x[z][y][x]-ccur_x[z][y][x])/(fl+tmp);
	  if(x==edge[0][y][z][1]){
	    ccur_x[z][y][x+1]=(fcur_x[z][y][x+1]-ccur_x[z][y][x+1])*fl_inv;
	  };
	  // YL-direction
	  tmp=0.;
	  if(!ledge_yl[x][y]){
	    tmp=mesh[meshid[z][y-1][x-corr]].GetFlux().get_dat(0);
	  };
          ccur_yl[z][y][x]=(fcur_yl[z][y][x]-ccur_yl[z][y][x])/(fl+tmp);
	  //
	  if(redge_yl[x][y]){
	    ccur_yl[z][y+1][x+(1-corr)]=(fcur_yl[z][y+1][x+1-corr]-ccur_yl[z][y+1][x+1-corr])*fl_inv;
	  };
	  // YR-direction
	  tmp=0.;
	  if(!ledge_yr[x][y]){
	    tmp=mesh[meshid[z][y-1][x+(1-corr)]].GetFlux().get_dat(0);
	  };
          ccur_yr[z][y][x]=(fcur_yr[z][y][x]-ccur_yr[z][y][x])/(fl+tmp);
	  if(redge_yr[x][y]){
	    int ttt=x-corr;
	    if(ttt==-1)ttt=xr;
	    ccur_yr[z][y+1][ttt]=(fcur_yr[z][y+1][ttt]-ccur_yr[z][y+1][ttt])*fl_inv;
	  };
	  // Z-direction
	  tmp=0.;
	  if(z!=0){
	    tmp=mesh[meshid[z-1][y][x]].GetFlux().get_dat(0);
	  };
          ccur_z[z][y][x]=(fcur_z[z][y][x]-ccur_z[z][y][x])/(fl+tmp);
	  if(z==zr-1){
	    ccur_z[z+1][y][x]=(fcur_z[z+1][y][x]-ccur_z[z+1][y][x])*fl_inv;
	  };
	};
      };
    };
  };
};

void DHEXSystem::ModCoef(vector< vector< vector<real> > > &moddif_x,
               vector< vector< vector<real> > > &moddif_yl,
               vector< vector< vector<real> > > &moddif_yr,
               vector< vector< vector<real> > > &moddif_z)
{
  int xr=mi->GetXC();
  int yr=mi->GetYC();
  int zr=mi->GetZC();

  for(int z=0;z<zr;z++){
    for(int y=0;y<yr;y++){
      int corr=y%2;
      for(int x=0;x<xr;x++){
	int id=meshid[z][y][x];
	if(id!=-1){
	  // X
          real tmp=-moddif_x[z][y][x]+moddif_x[z][y][x+1];
	  coefx[0][id]-=moddif_x[z][y][x];
	  moddif_x[z][y][x]*=2.;
	  // YL
          tmp+=-moddif_yl[z][y][x]+moddif_yl[z][y+1][x+1-corr];
	  coefyl[0][id]-=moddif_yl[z][y][x];
	  moddif_yl[z][y][x]*=2.;
	  // YR
	  int ttt=x-corr;
	  if(ttt==-1)ttt=xr;
	  tmp+=-moddif_yr[z][y][x]+moddif_yr[z][y+1][ttt];
	  coefyr[0][id]-=moddif_yr[z][y][x];
	  moddif_yr[z][y][x]*=2.;
	  // Z
	  tmp+=-moddif_z[z][y][x]+moddif_z[z+1][y][x];
	  coefz[0][id]-=moddif_z[z][y][x];
	  moddif_z[z][y][x]*=2.;
          coefs[0][id]+=tmp;
	};
      };
    };
  };
};

void DHEXSystem::CalCoarseCur(int g, 
                              vector< vector< vector<real> > > &curx,
			      vector< vector< vector<real> > > &curyl,
			      vector< vector< vector<real> > > &curyr,
			      vector< vector< vector<real> > > &curz, bool cumulative)
{
  int xr=mi->GetXC();
  int yr=mi->GetYC();
  int zr=mi->GetZC();

  if(!cumulative){
    for(int i=0;i<zr+1;i++){
      for(int j=0;j<yr+1;j++){
        for(int k=0;k<xr+1;k++){
	  curx[i][j][k]=0.;
	  curyl[i][j][k]=0.;
	  curyr[i][j][k]=0.;
	  curz[i][j][k]=0.;
	};
      };
    };
  };

  int iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi->GetCMeshF(2,z1);z2++){
      for(int y1=0;y1<yr;y1++){
        for(int x1=0;x1<xr;x1++){
          int index=meshid[iz][y1][x1];
          if(index!=-1){
            // X-direction
            curx[z1][y1][x1]+=CalCurrentX(x1,y1,iz,g);
            if(x1==edge[0][y1][iz][1])curx[z1][y1][x1+1]+=CalCurrentX(x1+1,y1,iz,g);
            // YL-direction
 	    curyl[z1][y1][x1]+=CalCurrentYL(x1,y1,iz,g);
	    int tmpx=x1;
	    if(y1%2==0)tmpx++;
	    if(redge_yl[x1][y1])curyl[z1][y1+1][tmpx]+=CalCurrentYL(tmpx,y1+1,iz,g);
	    //                        ~~~~~~~~~~~~~~ caution !!
	    // YR-direction
	    curyr[z1][y1][x1]+=CalCurrentYR(x1,y1,iz,g);
	    tmpx=x1;
	    if(y1%2==1)tmpx--;
	    if(redge_yr[x1][y1]){
	      int ttt=tmpx;
	      if(tmpx==-1)ttt=xr;
              curyr[z1][y1+1][ttt]+=CalCurrentYR(tmpx,y1+1,iz,g);
	    }; //  ~~~~~~~~~~~~~~~  caution !!
            // Z-direction
	    if(z2==0)curz[z1][y1][x1]+=CalCurrentZ(x1,y1,iz,g);
            if(iz==mi->GetZF()-1)curz[zr][y1][x1]+=CalCurrentZ(x1,y1,iz+1,g);
	  };
	};
      };
      iz++;
    };
  };
};



// **********************************************
// * For perturbation calculation               *
// **********************************************

real DHEXSystem::CalPerturbScatteringTerm(DHEXSystem *sec,int g,bool *flag)
// The same as PLOS
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real f2v=sec->GetMesh(m).GetFlux(0).get_dat(g)*mesh[m].GetVolume();
      real f1=mesh[m].GetFlux(0).get_dat(g);
      real tmp=0.;
      //for(int j=0;j<grp;j++){ // w up scattering
      for(int j=g+1;j<grp;j++){ // w up scattering
        tmp+=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g,j)
	             -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,j))*
             (f1-mesh[m].GetFlux(0).get_dat(j));
      };
      ret-=tmp*f2v;
    };
  };
  return ret;
};

real DHEXSystem::CalPerturbScatteringTerm(DHEXSystem *sec,int g1,int g2,bool *flag)
// The same as PLOS
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real tmp=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g1,g2)
 	               -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g1,g2))*
	     (mesh[m].GetFlux(0).get_dat(g1)-mesh[m].GetFlux(0).get_dat(g2));
      ret-=tmp*sec->GetMesh(m).GetFlux(0).get_dat(g1)*mesh[m].GetVolume();
    };
  };
  return ret;
};

real DHEXSystem::CalPerturbN2NTerm(DHEXSystem *sec,bool *flag,int g)
// The same as PLOS
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real dsign2n=sec->GetMesh(m).GetMed()->GetMacxs().GetSign2n().get_dat(g)
       	                  -mesh[m].GetMed()->GetMacxs().GetSign2n().get_dat(g);
      real vol=mesh[m].GetVolume();
      ret+=mesh[m].GetFlux(0).get_dat(g)*dsign2n*sec->GetMesh(m).GetFlux(0).get_dat(g)*vol;
    };
  };
  return ret;
};

real DHEXSystem::CalPerturbLeakageTerm(int dir,DHEXSystem *sec,int i,bool *flag)
{
  int xm=mi->GetXF();
  int ym=mi->GetYF();
  int zm=mi->GetZF();

  real ret=0.;
  for(int z=0;z<zm;z++){
    real zlen=mi->GetFMeshL(2,z);
    for(int y=0;y<ym;y++){
      for(int x=0;x<xm;x++){
	int m=meshid[z][y][x];
	if(m!=-1){
	  if(flag[m]){
	    real surface=hex_area;
	    if(dir!=2)surface=zlen*len_plane;
            real vr=mesh[m].GetVolume()*0.5;
            real vl=vr;
	    if(dir!=2){
	      vr=pitch*len_plane*zlen*0.5;
	      vl=vr;
	    };

            int difc1=GetDifc(dir);
            int difc2=sec->GetDifc(dir);

            real d1=mesh[m].GetDif(i,difc1);
            real d2=sec->GetMesh(m).GetDif(i,difc2);
            real cl1=0.;
            real cl2=0.;
	    real cr1=0.;
	    real cr2=0.;
            if(dir==0){
              cl1=CalCurrentX(x,y,z,i);
              cl2=sec->CalCurrentX(x,y,z,i);
              cr1=CalCurrentX(x+1,y,z,i);
              cr2=sec->CalCurrentX(x+1,y,z,i);
	    };
	    if(dir==1){
	      cl1=CalCurrentYL(x,y,z,i);
	      cl2=sec->CalCurrentYL(x,y,z,i);
	      int xpos=x;
	      if(y%2==0)xpos++;
	      cr1=CalCurrentYL(xpos,y+1,z,i);
	      cr2=sec->CalCurrentYL(xpos,y+1,z,i);
	    };
	    if(dir==2){
	      cl1=CalCurrentZ(x,y,z,i);
	      cl2=sec->CalCurrentZ(x,y,z,i);
	      cr1=CalCurrentZ(x,y,z+1,i);
	      cr2=sec->CalCurrentZ(x,y,z+1,i);
	    };

            cl1/=(d1*surface); 
            cl2/=(d2*surface);
            cr1/=(d1*surface);
            cr2/=(d2*surface);
            ret+=-(d2-d1)*(cl1*cl2*vl+cr1*cr2*vr);

	    if(dir==1){
	      cl1=CalCurrentYR(x,y,z,i);
	      cl2=sec->CalCurrentYR(x,y,z,i);
	      int xpos=x;
	      if(y%2==1)xpos--;
	      cr1=CalCurrentYR(xpos,y+1,z,i);
	      cr2=sec->CalCurrentYR(xpos,y+1,z,i);
              cl1/=(d1*surface); 
              cl2/=(d2*surface);
              cr1/=(d1*surface);
              cr2/=(d2*surface);
              ret+=-(d2-d1)*(cl1*cl2*vl+cr1*cr2*vr);
	    };

	  };
	};
      };
    };
  };
  return ret;
};


real DHEXSystem::CalReactivity(DHEXSystem *sec,real kunp,real kp,bool pr)
{
  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };
  real ret=CalReactivity(sec,kunp,kp,flag,pr);
  delete [] flag;
  return ret;
};

real DHEXSystem::CalReactivity(DHEXSystem *sec,real kunp,real kp,bool* flag,bool pr)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in DHEXSystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in DHEXSystem::CalReactivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  real *yld=new real[grp];
  real *abs=new real[grp];
  real *sct=new real[grp];
  real *n2n=new real[grp];
  real *lx =new real[grp];
  real *ly =new real[grp];
  real *lz =new real[grp];

  real ip=CalPerturbDenominator(sec);
  for(int i=0;i<grp;i++){
    sct[i]=CalPerturbScatteringTerm(sec,i,flag);
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    abs[i]=CalPerturbAbsorptionTerm(sec,flag,i);
    n2n[i]=CalPerturbN2NTerm(sec,flag,i);
    lx[i]=CalPerturbLeakageTerm(0,sec,i,flag);
    ly[i]=CalPerturbLeakageTerm(1,sec,i,flag);
    lz[i]=CalPerturbLeakageTerm(2,sec,i,flag);
  };

  real yldsum=0.;
  real abssum=0.;
  real sctsum=0.;
  real n2nsum=0.;
  real lxsum=0.;
  real lysum=0.;
  real lzsum=0.;
  real inv_ip=1./ip;
  if(pr){
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scattering  ";
    cout<<"Leakage     ";
    cout<<"(n,2n)      ";
    cout<<"Total\n";
  };

  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    sct[i]*=inv_ip;
    n2n[i]*=inv_ip;
    lx[i]*=inv_ip;
    ly[i]*=inv_ip;
    lz[i]*=inv_ip;
    yldsum+=yld[i];
    abssum+=abs[i];
    sctsum+=sct[i];
    n2nsum+=n2n[i];
    lxsum+=lx[i];
    if(Ndim>1)lysum+=ly[i];
    if(Ndim>2)lzsum+=lz[i];
    if(pr){
      cout.width(3);
      cout<<i<<"   ";
      cout.setf(ios::scientific);
      cout.precision(4);
      real ltmp=lx[i];
      if(Ndim>1)ltmp+=ly[i];
      if(Ndim>2)ltmp+=lz[i];
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real tot=yld[i]+abs[i]+sct[i]+ltmp+n2n[i];
      //real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      //real leth=(log(en/en_next)/0.25);
      real leth=1.;
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sct[i]/leth;
      cout<<"  "<<ltmp/leth<<" "<<n2n[i]/leth<<" "<<tot/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  real tot=yldsum+abssum+sctsum+n2nsum+lxsum+lysum+lzsum;

  if(pr){
    cout<<"#\n";
    cout<<"# +++ Results +++\n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# Yield      : "<<yldsum<<"\n";
    cout<<"# Absorption : "<<abssum<<"\n";
    cout<<"# Scattering : "<<sctsum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";
    cout<<"# (Non-Leak) : "<<yldsum+abssum+sctsum+n2nsum<<"\n";
    cout<<"# Leakage-x  : "<<lxsum<<"\n";
    cout<<"# Leakage-y  : "<<lysum<<"\n";
    cout<<"# Leakage-z  : "<<lzsum<<"\n";
    cout<<"# (Leakage)  : "<<lxsum+lysum+lzsum<<"\n";
    cout.unsetf(ios::scientific);
    cout<<"# \n";
    cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
    cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
  };

  delete [] yld;
  delete [] abs;
  delete [] sct;
  delete [] n2n;
  delete [] lx;
  delete [] ly;
  delete [] lz;

  return tot;
};


void DHEXSystem::CalSensitivity(DHEXSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  for(int i=0;i<Ndim;i++){
    if(GetDifc(i)!=0){
      cout<<"Anisotropic diffusion coefficient ";
      cout<<"cannot be used for sensitivity calculation.\n";
      exit(0);
    };
  };

  if(GetGeneralOption().Forward()){
    cout<<"Error in DHEXSystem::CalSensitivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in DHEXSystem::CalSensitivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  real *nsforg=new real[nmed];
  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  real *tot1org=new real[nmed];
  real *sigsorg=new real[nmed];
  real *sigs1org=new real[nmed];

  real *fiss_frac=new real[nmed];

  real delta=1.;

  cout.setf(ios::scientific);
  cout.precision(5);

  CheckSameMesh(sec);

  real ip=CalPerturbDenominator(sec);
  real inv_ip=1./ip;

  bool *flag=new bool[TotM];

  cout<<grp<<"\n";

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    //int nidendf=TranslateNuclideIDFromJFS(nid);
    int nidendf=nid;

    vector<real> den(nmed);
    for(int i=0;i<nmed;i++){
      if(med[i].ExistNuclide(nid)){
	den[i]=med[i].GetNuclide(nid).GetDensity();
      }else{
	den[i]=0.;
      };
    };

    for(int i=0;i<TotM;i++){
      if(mesh[i].GetMed()->ExistNuclide(nid)){
        flag[i]=true;
      }else{
	flag[i]=false;
      };
    };

    bool fissile=false;
    for(int i=0;i<nmed;i++){
      if(med[i].ExistNuclide(nid)){
        for(int j=0;j<grp;j++){
	  if(med[i].GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(j)>0.0){
	    fissile=true;
	    break;
	  };
        };
      };
      if(fissile)break;
    };

    if(fissile){
      // Fission
      cout<<1<<"\n";
      cout<<nidendf<<"\n";
      cout<<"  18\n";
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
            absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
            totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
  	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
  	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den[j]*micnu*micsigf);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den[j]*micsigf);
  	    real tmp=0.33333333/totorg[j];
	    sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den[j]*micsigf));
  	  };
        };
        real re=0.;
        re+=CalPerturbLeakageTerm(0,sec,i,flag);
        re+=CalPerturbYieldTerm(sec,flag,i,k2);
        re+=CalPerturbAbsorptionTerm(sec,flag,i);
        if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
        if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
        re*=inv_ip;
        cout<<"  "<<re<<"\n";
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
            sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
  	  };
        };
      };

      // Nu
      cout<<1<<"\n";
      cout<<nidendf<<"\n";
      cout<<"  452\n";
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
  	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den[j]*micsigf*micnu);
  	  };
        };
        real re=0.;
        re+=CalPerturbYieldTerm(sec,flag,i,k2);
        re*=inv_ip;
        cout<<"  "<<re<<"\n";
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
  	  };
        };
      };

      // Chi
      cout<<1<<"\n";
      cout<<nidendf<<"\n";
      cout<<"  181\n";
      for(int i=0;i<nmed;i++){
	if(med[i].ExistNuclide(nid)){
 	  real total_fiss=sec->GetIntegratedFlux(i)*sec->GetMed(i).GetMacxs().GetNusigf();
          real part_fiss=den[i]*(sec->GetIntegratedFlux(i)*sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicNusigf());
          if(total_fiss!=0.){
  	    fiss_frac[i]=part_fiss/total_fiss;
	  }else{
	    fiss_frac[i]=0.;
	  };
	}else{
	  fiss_frac[i]=0.;
	};
      };
      for(int i=0;i<grp;i++){
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    totorg[j]=sec->GetMed(j).GetMacxs().GetKai().get_dat(i);
	    sec->GetMed(j).GetMacxs().GetKai().add_data(i,fiss_frac[j]*1.);
	  };
	};
	real re=0.;
        re+=CalPerturbYieldTerm(sec,flag,i,k2);
	re*=inv_ip;
	cout<<"  "<<re<<"\n";
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetKai().put_data(i,totorg[j]);
	  };
	};
      };

    };

    // Capture
    cout<<1<<"\n";
    cout<<nidendf<<"\n";
    cout<<"  102\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den[j]*micsigc);
	  real tmp=0.33333333/totorg[j];
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den[j]*micsigc));
        };
      };
      real re=0.;
      re+=CalPerturbLeakageTerm(0,sec,i,flag);
      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
      if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
      re*=inv_ip;
      cout<<"  "<<re<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	};
      };
    };

    // Elastic P0
    cout<<2<<"\n";
    cout<<nidendf<<"\n";
    cout<<2<<"\n";
    for(int i=0;i<grp;i++){
      for(int k=i;k<grp;k++){
	if(k<=i+2){
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
              totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
              sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
	      real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
              sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den[j]*micsigs);
 	      real tmp=0.33333333/totorg[j];
              real vmu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(mu).get_dat(i);
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+(1-vmu)*den[j]*micsigs));
	    };
	  };
          real re=0.;
          re+=CalPerturbLeakageTerm(0,sec,i,flag);
          if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
          if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
          re+=CalPerturbScatteringTerm(sec,i,k,flag);
          re*=inv_ip;
          cout<<"  "<<re<<"\n";
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	      sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	    };
	  };
        }else{
          cout<<" 0.\n";
	};
      };
    };

    // (scattering term)
    vector< vector<real> > scat_term(grp);
    for(int i=0;i<grp;i++){
      scat_term[i].resize(grp,0.);
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den[j]*delta);
	  };
	};
        scat_term[i][k]=CalPerturbScatteringTerm(sec,i,k,flag);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	  };
	};
      };
    };
    // Inelastic & n2n P0
    for(int ii=0;ii<2;ii++){
      cout<<2<<"\n";
      cout<<nidendf<<"\n";
      if(ii==0){cout<<4<<"\n";} // inelastic
      else {cout<<16<<"\n";}; // n2n
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	    sigsorg[j]=sec->GetMed(j).GetMacxs().GetSign2n().get_dat(i);
            real dtot=delta;
	    if(ii==1)dtot=delta*0.5;
            real tmp=0.33333333/totorg[j];
            sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den[j]*dtot));
	    if(ii==1)sec->GetMed(j).GetMacxs().GetSign2n().add_data(i,den[j]*dtot);
	  };
	};
        real leak=CalPerturbLeakageTerm(0,sec,i,flag);
        if(Ndim>1)leak+=CalPerturbLeakageTerm(1,sec,i,flag);
        if(Ndim>2)leak+=CalPerturbLeakageTerm(2,sec,i,flag);
	if(ii==1)leak+=CalPerturbN2NTerm(sec,flag,i);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	    if(ii==1)sec->GetMed(j).GetMacxs().GetSign2n().put_data(i,sigsorg[j]);
          };
	};
        for(int k=i;k<grp;k++){
          real re=leak+scat_term[i][k];
          re*=inv_ip;
          cout<<"  "<<re<<"\n";
	};
      };
    };

    // Elastic-p1(mu)
    cout<<1<<"\n";
    cout<<nidendf<<"\n";
    cout<<"  251\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	  real vsigel=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigel).get_dat(i);
	  real tmp=0.33333333/totorg[j];
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp-den[j]*delta*vsigel));
        };
      };
      real re=0.;
      re+=CalPerturbLeakageTerm(0,sec,i,flag);
      if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
      if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
      re*=inv_ip;
      cout<<"  "<<re<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	};
      };
    };
  };

  delete [] fiss_frac;
  delete [] flag;
  delete [] sigsorg;
  delete [] sigs1org;
  delete [] absorg;
  delete [] nsforg;
  delete [] totorg;
  delete [] tot1org;
};

void DHEXSystem::PrintFluxDistribution(int zz)
{
  int xx=mi->GetFMesh(0);
  int yy=mi->GetFMesh(1);
  for(int y=0;y<yy;y++){
    if(y%2==0)cout<<" ";
    for(int x=0;x<xx;x++){
      int id=meshid[zz][y][x];
      if(id==-1){cout<<" ";}
      else{
        real tmp=GetMesh(id).GetFlux().get_dat(0);
        if(tmp>0){
	  cout<<"O";
        }else{
	  cout<<"x";
        };
      };
      cout<<" ";
    };
    cout<<"\n";
  };
};

void DHEXSystem::MemoryReductionForHighMomentCal()
{
  for(int i=0;i<nmed;i++){
    med[i].MacxsVectorClear();
  };
  coefs.clear();
  coefx.clear();
  coefyl.clear();
  coefyr.clear();
  coefz.clear();
  for(int i=0;i<TotM;i++){
    mesh[i].SrcDataClear();
  };
};

void DHEXSystem::CopyCoef(DHEXSystem &sec)
{
  coefs.resize(grp);
  coefx.resize(grp);
  coefyl.resize(grp);
  coefyr.resize(grp);
  coefz.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotM,0.);
    coefx[i].resize(TotM,0.);
    coefyl[i].resize(TotM,0.);
    coefyr[i].resize(TotM,0.);
    coefz[i].resize(TotM,0.);
    for(int j=0;j<TotM;j++){
      coefs[i][j]=sec.GetCoefs(i,j);
      coefx[i][j]=sec.GetCoefx(i,j);
      coefyl[i][j]=sec.GetCoefyl(i,j);
      coefyr[i][j]=sec.GetCoefyr(i,j);
      coefz[i][j]=sec.GetCoefz(i,j);
    };
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };
};

void DHEXSystem::ModifiedCoarseMeshCorrectionCal(real keff)
{
  // 
  // correction factor is [1/9 h^2 b^2]
  //

  real pp=pitch*pitch/9.;

  cal_mcor=true;
  mcor.resize(TotM);
  for(int i=0;i<TotM;i++){
    mcor[i].resize(grp);
  };

  for(int i=0;i<TotM;i++){
    real fsrc=mesh[i].GetFissionSrc()/mesh[i].GetVolume()/keff;
    for(int g=0;g<grp;g++){
      real sss=0.;
      real ch=0.;
      if(opt.Forward()){
        for(int g2=0;g2<g;g2++){
	  sss+=mesh[i].GetFlux().get_dat(g2)*
	       mesh[i].GetMed()->GetMacxs().GetData2d(sigs).get_dat(g2,g);
        };
        ch=mesh[i].GetMed()->GetMacxs().GetData1d(chi).get_dat(g);
      }else{
        for(int g2=g+1;g2<grp;g2++){
	  sss+=mesh[i].GetFlux().get_dat(g2)*
  	       mesh[i].GetMed()->GetMacxs().GetData2d(sigs).get_dat(g,g2);
        };
        ch=mesh[i].GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g);
      };
      real flx_inv=1./mesh[i].GetFlux().get_dat(g);
      real sr=mesh[i].GetMed()->GetMacxs().GetSigr().get_dat(g);
      real tmp=ch*fsrc*flx_inv+sss*flx_inv-sr;
      mcor[i][g]=tmp/mesh[i].GetMed()->GetMacxs().GetData1d(d).get_dat(g)*pp;
    };
  };
};

real DHEXSystem::CalIgenMCMCor(int pl)
{
  //real ktmp=CalIgen();
  //ModifiedCoarseMeshCorrectionCal(ktmp);
  //CalCoef(true);

  real fiss,fissold,errk,errf,errs;
  errf=1.0;
  errs=1.0;

  //SetInitialFlux();
  NormalizeFissionSrc();
  
  real cin=0.1;
  fiss=1.;

  //cmfd=false;
  string acname="none"; // o
  if(cmfd)acname="CMFD"; // o
  WriteSolverName(acname);

  for(int iter=0;iter<opt.GetOutitermax();iter++){

    if(errf*0.02<cin)cin=errf*0.02;
    real tmp=5e-5;
    if(opt.GetEpsf()<5e-5)tmp=opt.GetEpsf()*0.5;
    if(cin<tmp)cin=tmp;

    errf=0.;
    bool InnerConvergence=true;

    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,fiss,pl);
      real err=CalFluxGeneral(ginp,cin,iter);
      if(err<0.){
	InnerConvergence=false;
	err=-err;
      };
      if(err>errf)errf=err;
      // +++++
      for(int m=0;m<TotM;m++){
	real org=mesh[m].GetFlux().get_dat(ginp);
	real oth=mcor[m][ginp];
	mesh[m].GetFlux().put_data(ginp,org*(1.-oth));
      };
      // +++++
      AddDownScatSrc(ginp,pl);
      SetZeroScatSrc(ginp);
    };


    if(cmfd){
      DoAcceleration(iter,errs,fiss); // o
      // Correction for up-scattering
      if(IsUpScatterSystem){
        for(int g=UpScatSinkGrp;g<grp;g++){
          int ginp=g;
          if(!opt.Forward())ginp=grp-1-g;
          if(IsUpScatterSystem)SetZeroUpScatSrc(ginp);
          AddUpScatSrc(ginp,pl);
	};
      };
    };


    errs=0.;
    for(int i=0;i<TotM;i++){
      real tmp=mesh[i].GetFissionSrc();
      if(tmp>0.){
        if(opt.Forward()){mesh[i].CalFissionSrc();}
	else{mesh[i].CalFissionSrcAdjoint();};
	real tmp2=mesh[i].GetFissionSrc();
	real err=fabs(tmp2/tmp-1.0);
	if(err>errs)errs=err;
      };
    };

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    errf=1e-10;
    WriteIterationInfo(iter,fiss,errk,errf,errs);

    if(opt.Converged(errf,errk,errs)){
      if(InnerConvergence){break;}
      else{cout<<"  ... Inner iteration is not converged ...\n";
      };
    };
  };

  return fiss;
}


/*


real PLOSSystem::CalEdgeCurZeroFlux(int m,int dir,int g,int flag)
{
  return mesh[m].GetDif(g,flag)/(mesh[m].GetLen(dir)*0.5)*-1.;
};

*/
