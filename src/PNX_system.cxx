#include <cstdlib>
#include "PNX_system.h"

using namespace std;

PNXSystem::PNXSystem(int n,int g,int i):GeneralSystem(n,g,i)
{
  name="PNX";
  cmfdimp=false; // CMFD is not implemented
  plexp=-1;
};

PNXSystem::~PNXSystem()
{
};

void PNXSystem::PutPLexp(int i)
{
  plexp=i;
  if(plexp%2==0){
    cout<<"Error in PutPLexp.\n";
    cout<<"Odd number should be used.\n";
    exit(0);
  };

  plexp1=plexp+1;

  pl_pl.resize(plexp1);
  for(int i=0;i<plexp1;i++){
    pl_pl[i].resize(plexp1);
  };

  int m_u=1000;
  real w_u=1./m_u;
  for(int i=0;i<plexp1;i++){
    for(int j=i;j<plexp1;j++){
      real sum=0.;
      for(int u=0;u<m_u;u++){
	real vu=w_u*u+w_u*0.5;
	//real vu=-1.+w_u*2.*u+w_u;
	sum+=w_u*Legendre(i,vu)*Legendre(j,vu);
      };
      pl_pl[i][j]=sum;
      if(i!=j)pl_pl[j][i]=sum;
    };
  };

  //pl_pl[1][1]*=2.1312/2.*1.065;
  //pl_pl[1][1]*=2.1312/2.;
  //pl_pl[1][1]*=1.731/2.;
};

void PNXSystem::PutCartMeshInfo(CartMeshInfo cm)
{
  mi=cm;

  TotM=mi.GetXF();
  mesh.resize(TotM);

  meshid.resize(1);
  meshid[0].resize(1);
  meshid[0][0].resize(TotM);
  for(int i=0;i<TotM;i++){
    meshid[0][0][i]=i;
  };

  real l[3];
  for(int i=0;i<TotM;i++){
    int tm=mi.GetFMat(i);
    if(tm>=nmed){
      cout<<"Error in PNXSystem::PutCartMeshInfo.\n";
      cout<<"You requested not-existing medium ID.\n";
      cout<<"Please check Medium ID.\n";
      exit(0);
    };
    mesh[i].PutMedium(&med[mi.GetFMat(i)]);
    l[0]=mi.GetFMeshL(0,i);
    mesh[i].PutDim(1,l);
  };

  for(int i=0;i<2;i++){
    BC[i]=-1;
    if(mi.GetBC(i)==1)BC[i]=1; // Reflective
    if(mi.GetBC(i)==2)BC[i]=2; // Vacuum
    if(BC[i]==-1){
      cout<<"Incorrect B. C. in PNX_system.\n";
      exit(0);
    };
  };
};

void PNXSystem::SetArray()
{
  if(pl==-1){
    cout<<"Do 'PutPL'.\n";
    exit(0);
  };
  if(plexp==-1){
    cout<<"Do PutPLexp.\n";
    exit(0);
  };

  if(plexp<pl){
    cout<<"Error in PNX::SetArray.\n";
    cout<<"Pl order for flux should be larger than that of scattering XS.\n";
    exit(0);
  };

  for(int i=0;i<TotM;i++){
    mesh[i].PutPL(plexp);
    mesh[i].PutGrp(grp);
  };

  pos.resize(plexp1);

  sz=0;
  if(BC[0]==1&&BC[1]==1){
    for(int l=0;l<plexp1;l++){
      pos[l]=sz;
      if(l%2==0){
        sz+=TotM+1;
      }else{
        sz+=TotM-1; // ref/ref
      };
    };
  }else if(BC[0]==1&&BC[1]==2){
    for(int l=0;l<plexp1;l++){
      pos[l]=sz;
      if(l%2==0){
        sz+=TotM+1;
      }else{
        sz+=TotM; // ref/vac
      };
    };
  }else{
    cout<<"Not coded.\n";
    exit(0);
  };

  CalAmat();

  // Current-weighted total cross section check
  for(int i=0;i<nmed;i++){
    if(med[i].GetPLT()>=1){
      bool zero_xs=false;
      for(int g=0;g<grp;g++){
	if(med[i].GetMacxs().GetData1d(sigt,1).get_dat(g)==0.){
	  zero_xs=true;
	};
      };
      if(zero_xs){
	cout<<"!! Warning !!\n";
	cout<<" Zero current-weighted total cross section is detected.\n";
      };
    };
  };
};

void PNXSystem::SetInitialFlux()
{
  SetArray();
  SetInitialFlatFlux();
};

real PNXSystem::CalFluxGeneral(int g,real cin,int iter)
{
  if(BC[0]==1&&BC[1]==1){
    return CalFluxReflectiveReflective(g);
  }else if(BC[0]==1&&BC[1]==2){
    return CalFluxReflectiveVacuum(g);
    //return CalFluxReflectiveVacuumCullen(g);
  }else{
    cout<<"Not coded.\n";
    exit(0);
  };
};

real PNXSystem::CalFluxReflectiveReflective(int g)
{
  // Only for reflective boundary conditions
  //
  // Boundary condition is implicitely treated

  vector<real> amat(sz*sz,0.);
  vector<real> bmat(sz,0.);

  // spatial integration each mesh
  int index=0;
  for(int l=0;l<plexp1;l++){
    real fact=1./(2.*l+1.)*PI4;
    for(int m=0;m<TotM;m++){
      if(l<plnum)bmat[index]=mesh[m].GetSrcin(l)*fact; // out-group source
      index++;
    };
  };

  for(int i=0;i<sz*sz;i++){amat[i]=amatg[g][i];};
  gauss(amat,bmat,sz,1);
  /*
  vector<real> sumt(sz);
  for(int i=0;i<sz;i++){
    real sum=0.;
    for(int j=0;j<sz;j++){
      sum+=amatg[g][i*sz+j]*bmat[j];
    };
    sumt[i]=sum;
  };
  for(int i=0;i<sz;i++){
    bmat[i]=sumt[i];
  };
  */


  real errf=0.;
  // Put new flux into GeneralMesh
  index=0;
  for(int l=0;l<plexp1;l++){
    for(int i=0;i<TotM;i++){
      real err=0.;
      real newflx;
      if(l%2==0){
        newflx=(bmat[index]+bmat[index+1])*0.5;
	index++;
	if(i==TotM-1)index++;
      }else{
        if(i==0){
  	  newflx=0.5*bmat[index];
	}else if(i==TotM-1){
  	  newflx=0.5*bmat[index];
          index++;
        }else{
          newflx=(bmat[index]+bmat[index+1])*0.5;
  	  index++;
	};
      };
      if(l==0)err=fabs(newflx/mesh[i].GetFlux(l).get_dat(g)-1.0);
      mesh[i].GetFlux(l).put_data(g,newflx);
      if(err>errf)errf=err;
    };
  };

  return errf;
};

real PNXSystem::CalFluxReflectiveVacuum(int g)
{
  // Boundary condition is implicitely treated

  vector<real> amat(sz*sz,0.);
  vector<real> bmat(sz,0.);

  // spatial integration each mesh
  int index=0;
  for(int l=0;l<plexp1;l++){
    real fact=1./(2.*l+1.)*PI4;
    for(int m=0;m<TotM;m++){
      if(l<plnum)bmat[index]=mesh[m].GetSrcin(l)*fact; // out-group source
      index++;
    };
  };

  for(int i=0;i<sz*sz;i++){amat[i]=amatg[g][i];};

  gauss(amat,bmat,sz,1);

  real errf=0.;
  // Put new flux into GeneralMesh
  index=0;
  for(int l=0;l<plexp1;l++){
    for(int i=0;i<TotM;i++){
      real err=0.;
      real newflx;
      if(l%2==0){
        newflx=(bmat[index]+bmat[index+1])*0.5;
	index++;
	if(i==TotM-1)index++;
      }else{
        if(i==0){
  	  newflx=0.5*bmat[index];
        }else{
          newflx=(bmat[index]+bmat[index+1])*0.5;
	  index++;
	  if(i==TotM-1)index++;
	};
      };
      if(l==0)err=fabs(newflx/mesh[i].GetFlux(l).get_dat(g)-1.0);
      mesh[i].GetFlux(l).put_data(g,newflx);
      if(err>errf)errf=err;
    };
  };

  return errf;
};

real PNXSystem::CalFluxReflectiveVacuumCullen(int g)
{
  int plmed=3;

  // Boundary condition is implicitely treated

  vector<real> amat(sz*sz,0.);
  vector<real> bmat(sz,0.);

  // spatial integration each mesh
  int index=0;
  for(int l=0;l<plexp1;l++){
    real fact=1./(2.*l+1.)*PI4;
    for(int m=0;m<TotM;m++){
      if(l<plnum)bmat[index]=mesh[m].GetSrcin(l)*fact; // out-group source
      index++;
    };
  };

  //for(int i=0;i<sz*sz;i++){amat[i]=amatg[g][i];};
  int sz2=0;
  for(int l=0;l<=plmed;l++){
    if(l%2==0){
      sz2+=TotM+1;
    }else{
      sz2+=TotM;
    };
  };

  index=0;
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      if(i<sz2&&j<sz2){
	amat[index++]=amatg[g][i*sz+j];
      };
    };
  };

  gauss(amat,bmat,sz2,1);

  real errf=0.;
  // Put new flux into GeneralMesh
  index=0;
  for(int l=0;l<=plmed;l++){
    for(int i=0;i<TotM;i++){
      real err=0.;
      real newflx;
      if(l%2==0){
        newflx=(bmat[index]+bmat[index+1])*0.5;
	index++;
	if(i==TotM-1)index++;
      }else{
        if(i==0){
  	  newflx=0.5*bmat[index];
        }else{
          newflx=(bmat[index]+bmat[index+1])*0.5;
	  index++;
	  if(i==TotM-1)index++;
	};
      };
      if(l==0)err=fabs(newflx/mesh[i].GetFlux(l).get_dat(g)-1.0);
      mesh[i].GetFlux(l).put_data(g,newflx);
      if(err>errf)errf=err;
    };
  };

  // (l=plmed+2)
  int st=sz2-TotM-1;
  real fll=0.;
  real fln;
  int plmed1=plmed+1;
  for(int i=0;i<TotM;i++){
    mesh[i].GetFlux(plmed1).put_data(g,0.);
    if(i==0){
      fln=fll-real(plmed1)/(plmed1+1)*bmat[st+i];
    }else{
      fln=fll-real(plmed1)/(plmed1+1)*(bmat[st+i]-bmat[st+i-1]);
    };
    mesh[i].GetFlux(plmed+2).put_data(g,(fll+fln)*0.5);
    fll=fln;
  };


  sz2+=TotM+1;

  int sz3=0;
  for(int l=plmed+2;l<plexp1;l++){
    if(l%2==0){
      sz3+=TotM+1;
    }else{
      sz3+=TotM;
    };
  };

  for(int i=0;i<sz2;i++){bmat[i]=0.;};

  index=0;
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      if(i>=sz2&&j>=sz2){
	amat[index++]=amatg[g][i*sz+j];
      };
    };
  };

  gauss(amat,bmat,sz3,1);
  /*
  // Put new flux into GeneralMesh
  index=0;
  for(int l=plmed+2;l<=plexp;l++){
    for(int i=0;i<TotM;i++){
      real newflx;
      if(l%2==0){
        newflx=(bmat[index]+bmat[index+1])*0.5;
	index++;
	if(i==TotM-1)index++;
      }else{
        if(i==0){
  	  newflx=0.5*bmat[index];
        }else{
          newflx=(bmat[index]+bmat[index+1])*0.5;
	  index++;
	  if(i==TotM-1)index++;
	};
      };
      mesh[i].GetFlux(l).put_data(g,newflx);
    };
  };
  */
  return errf;
};

real PNXSystem::CalFluxGeneral2(int g,real epsif,int iter)
{
  int totm1=TotM+1;
  int sz=totm1*plnum;

  real *amat=new real[sz*sz];
  real *bmat=new real[sz];
  for(int i=0;i<sz*sz;i++){amat[i]=0.;};
  for(int i=0;i<sz;i++){bmat[i]=0.;};

  // spatial integration each mesh
  int index=0;
  for(int l=0;l<plnum;l++){
    for(int m=0;m<TotM;m++){
      real mlen=mesh[m].GetVolume();
      int isz=index*sz;
      if(l!=0){
	amat[isz+(l-1)*totm1+(m+1)]=l/(2.*l+1.);
	amat[isz+(l-1)*totm1+m]=-l/(2.*l+1.);
      };
      if(l!=pl){
	amat[isz+(l+1)*totm1+(m+1)]=(l+1.)/(2.*l+1.);
	amat[isz+(l+1)*totm1+m]=-(l+1.)/(2.*l+1.);
      };
      real s=mesh[m].GetMed()->GetMacxs().GetData1d(sigt,l).get_dat(g);
      s-=mesh[m].GetMed()->GetMacxs().GetData2d(sigs,l).get_dat(g,g)/(2.*l+1.);
      amat[isz+l*totm1+m]=s*0.5*mlen;
      amat[isz+l*totm1+(m+1)]=s*0.5*mlen;
      bmat[index]=mesh[m].GetSrcin(l)*PI4/(2.*l+1.); // out-group source
      index++;
    };
  };
  // boundary condition
  for(int l=0;l<plnum;l++){
    // (vacuum)
    /*
    if(l%2==1){
      for(int m=0;m<plnum;m++){
	amat[index*sz+m*totm1+0]=(2.*m+1.)*pl_pl[l][m]; // left
      };
      index++;
      for(int m=0;m<plnum;m++){
        real tmp=0.;
	if(l==m)tmp=2./(2.*l+1);
	amat[index*sz+m*totm1+TotM]=(2.*m+1.)*(tmp-pl_pl[l][m]); // right
      };
      index++;
    };
    */
    // (reflective)

    if(l%2==1){
      amat[index*sz+l*totm1+0]=1.;   // left 
      index++;
      amat[index*sz+l*totm1+TotM]=1.; // right
      index++;
    };

  };

  gauss(amat,bmat,sz,1);

  real errf=0.;
  // Put new flux into GeneralMesh
  for(int i=0;i<TotM;i++){
    real err=0.;
    for(int l=0;l<plnum;l++){
      real newflx=(bmat[l*totm1+i]+bmat[l*totm1+(i+1)])*0.5;
      if(l==0)err=fabs(newflx/mesh[i].GetFlux(l).get_dat(g)-1.0);
      mesh[i].GetFlux(l).put_data(g,newflx);
    };
    if(err>errf)errf=err;
  };

  delete [] amat;
  delete [] bmat;

  return errf;
};

real PNXSystem::GetAngularFlux(int m,int g,real mu,real pli)
{
  if(pli==-1)pli=plexp;

  real sum=0.;
  for(int l=0;l<=pli;l++){
    sum+=(2*l+1)*mesh[m].GetFlux(l).get_dat(g)*Legendre(l,mu);
  };
  return sum;
};

void PNXSystem::CalAmat()
{
  //vector<real> cmat(sz*sz,0.);

  amatg.resize(grp);
  for(int g=0;g<grp;g++){
    amatg[g].resize(sz*sz,0.);
    // (reflective/reflective)
    if(BC[0]==1&&BC[1]==1){
      int index=0;
      for(int l=0;l<plexp1;l++){
        real fact=1./(2.*l+1.);
        for(int m=0;m<TotM;m++){
          real mlen=mesh[m].GetVolume();
          int isz=index*sz;
          int mpos=0;
          if(l%2==1)mpos=1;
          if(l!=0){
  	    if(m!=TotM-1||mpos==1)amatg[g][isz+pos[l-1]+(m+1)-(1-mpos)]=l*fact;
	    if(m!=0||mpos==1)     amatg[g][isz+pos[l-1]+m-(1-mpos)]=-l*fact;
          };
          if(l!=plexp){
	    if(m!=TotM-1||mpos==1)amatg[g][isz+pos[l+1]+(m+1)-(1-mpos)]=(l+1.)*fact;
	    if(m!=0||mpos==1)     amatg[g][isz+pos[l+1]+m-(1-mpos)]=-(l+1.)*fact;
          };
          int lt=l;
          int plt=mesh[m].GetMed()->GetPLT();
          if(lt>plt)lt=plt;
          real s=mesh[m].GetMed()->GetMacxs().GetData1d(sigt,lt).get_dat(g);
          if(l<plnum)s-=mesh[m].GetMed()->GetMacxs().GetData2d(sigs,l).get_dat(g,g)*fact;
          real tmp=s*0.5*mlen;
          if(mpos==0){
            amatg[g][isz+pos[l]+m]=tmp;
            amatg[g][isz+pos[l]+(m+1)]=tmp;
          }else{
            if(m!=0)amatg[g][isz+pos[l]+m-1]=tmp;
            if(m!=TotM-1)amatg[g][isz+pos[l]+(m+1)-1]=tmp;
          };
          index++;
        };
      };
    }else{
      // (reflective/vacuum)
      int index=0;
      for(int l=0;l<plexp1;l++){
        real fact=1./(2.*l+1.);
        int mpos=0;
        if(l%2==1)mpos=1;
        for(int m=0;m<TotM;m++){
          real mlen=mesh[m].GetVolume();
          int isz=index*sz;
          if(l!=0){
	    amatg[g][isz+pos[l-1]+(m+1)-(1-mpos)]=l*fact;
   	    if(m!=0||mpos==1)amatg[g][isz+pos[l-1]+m-(1-mpos)]=-l*fact;
          };
          if(l!=plexp){
	    amatg[g][isz+pos[l+1]+(m+1)-(1-mpos)]=(l+1.)*fact;
    	    if(m!=0||mpos==1)amatg[g][isz+pos[l+1]+m-(1-mpos)]=-(l+1.)*fact;
          };
          int lt=0;
          if(l!=0&&mesh[m].GetMed()->GetPLT()>=1)lt=1; 
          real s=mesh[m].GetMed()->GetMacxs().GetData1d(sigt,lt).get_dat(g);
          if(l<plnum)s-=mesh[m].GetMed()->GetMacxs().GetData2d(sigs,l).get_dat(g,g)*fact;
          real tmp=s*0.5*mlen;
          if(mpos==0){
            amatg[g][isz+pos[l]+m]=tmp;
            amatg[g][isz+pos[l]+(m+1)]=tmp;
          }else{
            if(m!=0)amatg[g][isz+pos[l]+m-1]=tmp;
            amatg[g][isz+pos[l]+(m+1)-1]=tmp;
          };
          index++;
        };
        // Marshak boundary condition at right surface
        if(mpos==1){
          for(int m=0;m<plexp1;m++){
            real tmp=0.;
	    if(l==m)tmp=2.*fact;
            int mpos=0;
	    if(m%2==1)mpos=1;
	    amatg[g][index*sz+pos[m]+TotM-mpos]=(2.*m+1.)*(tmp-pl_pl[l][m]); // right
	    /*
	    if(l==1&&m==1){
	      amatg[g][index*sz+pos[m]+TotM-mpos]=(2.*m+1.)*(tmp-pl_pl[l][m])*1.731/2.; // S2-consistent (Mark BC)
	    };
	    */
          };
          index++;
        };
      };
    };
    /*
    for(int i=0;i<sz*sz;i++){cmat[i]=0.;};
    for(int i=0;i<sz;i++){cmat[i*sz+i]=1.;};
    gauss(amatg[g],cmat,sz,sz);
    for(int i=0;i<sz*sz;i++){amatg[g][i]=cmat[i];};
    */
  };
};

void PNXSystem::ShowAngularFlux(int idm,int g)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int im=-100;im<=100;im++){
    real mu=real(im*0.01);
    real sum=0.;
    for(int l=0;l<=plexp;l++){
      sum+=(2*l+1)*GetMesh(idm).GetFlux(l).get_dat(g)*Legendre(l,mu);
    };
    cout<<mu<<" "<<sum<<"\n";
  };

};
