#include <cstdlib>
#include "PLOS_system.h"

using namespace std;

#define heat_conduction_calculation 0
#define rwmc_calculation 0
  
void PLOSSystem::Init(int n,int g,int i)
{
  GeneralSystem::Init(n,g,i);
  Cylinder = false;
  Sphere   = false;
  difc[0]  = 0;
  difc[1]  = 0;
  difc[2]  = 0;
  name="PLOS";
  PutPL(0);
  cmfdimp = true; // CMFD is implemented in PLOS
  cmfd_factor=1.0;
  
  ip_input=false;

  xr_bc=-1.; //
}

void PLOSSystem::CalCoef()
{
  coefs.resize(grp);
  coefx.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotM,0.);
    coefx[i].resize(TotM,0.);
  };
  if(Ndim>1){
    coefy.resize(grp);
    for(int i=0;i<grp;i++){
      coefy[i].resize(TotM,0.);
    };
  };
  if(Ndim>2){
    coefz.resize(grp);
    for(int i=0;i<grp;i++){
      coefz[i].resize(TotM,0.);
    };
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };

  GroupData1D sigr(grp);
  real d1[3];
  real d2,dl,sfl[3],sfr[3],lenh[3];
  int id;
  int xyz[3];
  for(int z=0;z<mi.GetZF();z++){
    for(int y=0;y<mi.GetYF();y++){
      for(int x=0;x<mi.GetXF();x++){
	int i=meshid[z][y][x];
	if(i!=-1){
  	  xyz[0]=x;
	  xyz[1]=y;
	  xyz[2]=z;
	  sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
	  real vol=mesh[i].GetVolume();
	  for(int nd=0;nd<Ndim;nd++){
	    sfl[nd]=mesh[i].GetSurL(nd);
	    sfr[nd]=mesh[i].GetSurR(nd);
	    lenh[nd]=mesh[i].GetLen(nd)*0.5;
	  };
          for(int g=0;g<grp;g++){
	    coefs[g][i]=sigr.get_dat(g)*vol;
            for(int nn=0;nn<Ndim;nn++){
  	      d1[nn]=mesh[i].GetDif(g,difc[nn]);
	    };
  	    for(int nd=0;nd<Ndim;nd++){
	      bool ledge=false;
	      bool redge=false;
              if(nd==0){
	        if(x==edge[0][y][z][0])ledge=true;
	        if(x==edge[0][y][z][1])redge=true;
	      }else if(nd==1){
                if(y==edge[1][x][z][0])ledge=true;
	        if(y==edge[1][x][z][1])redge=true;
 	      }else{
                if(z==edge[2][x][y][0])ledge=true;
	        if(z==edge[2][x][y][1])redge=true;
	      };

              if(!ledge){
	        xyz[nd]-=1;
	        id=meshid[xyz[2]][xyz[1]][xyz[0]];
	        xyz[nd]+=1;
                d2=mesh[id].GetDif(g,difc[nd]);
	        dl=d1[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	        real tmp=d1[nd]*d2/dl*sfl[nd];
	        coefs[g][i]+=tmp;
	        switch(nd){
	        case 0:
                  coefx[g][i]=-tmp;
		  break;
	        case 1:
                  coefy[g][i]=-tmp;
		  break;
	        case 2:
	          coefz[g][i]=-tmp;
		  break;
	        };
              }else{
	        switch(BC[nd*2]){
	        case 0:
                  coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfl[nd];
                  break;
	        case 2:
                  coefs[g][i]-=CalEdgeCurVacuum(i,nd,g,difc[nd])*sfl[nd];
		  break;
	        };
              };
	      if(!redge){
	        xyz[nd]+=1;
	        id=meshid[xyz[2]][xyz[1]][xyz[0]];
	        xyz[nd]-=1;
                d2=mesh[id].GetDif(g,difc[nd]);
                dl=d1[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	        coefs[g][i]+=d1[nd]*d2/dl*sfr[nd];
              }else{
	        switch(BC[nd*2+1]){
	        case 0:
                  coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfr[nd];
	  	  break;
	        case 2:
                  coefs[g][i]-=CalEdgeCurVacuum(i,nd,g,difc[nd])*sfr[nd];
		  break;
		};
	      };
            };
	  };
	};
      };
    };
  };

  DetermineSizeCoef();
};


void PLOSSystem::CalCoef1DWithDF()
{
  coefs.resize(grp);
  coefx.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotM,0.);
    coefx[i].resize(TotM,0.);
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };

  GroupData1D sigr(grp);
  for(int i=0;i<mi.GetXF();i++){
    sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
    real vol=mesh[i].GetVolume();
    real lenh=vol*0.5;
    for(int g=0;g<grp;g++){
      // (Df is temporary stored in Medium::Flux)
      real df_l=mesh[i].GetMed()->GetFlux(0).get_dat(g);
      real df_r=mesh[i].GetMed()->GetFlux(1).get_dat(g);
      //
      coefs[g][i]=sigr.get_dat(g)*vol;
      real d1=mesh[i].GetDif(g,0);
      if(i!=0){
        real d2=mesh[i-1].GetDif(g,0);
	real dfl_r=mesh[i-1].GetMed()->GetFlux(1).get_dat(g);
	d2/=dfl_r;
        real dl=(d1/df_l)*mesh[i-1].GetLen(0)*0.5+d2*lenh;
        real tmp=(d1/df_l)*d2/dl;
        coefs[g][i]+=tmp*df_l;
        coefx[g][i]=-tmp;
      }else{
        switch(BC[0]){
        case 0:
          coefs[g][i]-=CalEdgeCurZeroFlux(i,0,g,0);
          break;
        case 2:
          coefs[g][i]-=CalEdgeCurVacuum(i,0,g,0);
          break;
	};
      };
      if(i!=mi.GetXF()-1){
        real d2=mesh[i+1].GetDif(g,0);
	real dfr_l=mesh[i+1].GetMed()->GetFlux(0).get_dat(g);
	d2/=dfr_l;
        real dl=(d1/df_r)*mesh[i+1].GetLen(0)*0.5+d2*lenh;
        coefs[g][i]+=(d1/df_r)*d2/dl*df_r;
      }else{
        switch(BC[1]){
        case 0:
          coefs[g][i]-=CalEdgeCurZeroFlux(i,0,g,0);
  	  break;
        case 2:
          coefs[g][i]-=CalEdgeCurVacuum(i,0,g,0);
	  break;
	};
      };
    };
  };

  DetermineSizeCoef();
};

void PLOSSystem::CopyCoef(PLOSSystem &sec)
{
  coefs.resize(grp);
  coefx.resize(grp);
  for(int i=0;i<grp;i++){
    coefs[i].resize(TotM,0.);
    coefx[i].resize(TotM,0.);
    for(int j=0;j<TotM;j++){
      coefs[i][j]=sec.GetCoefs(i,j);
      coefx[i][j]=sec.GetCoefx(i,j);
    };
  };
  if(Ndim>1){
    coefy.resize(grp);
    for(int i=0;i<grp;i++){
      coefy[i].resize(TotM,0.);
      for(int j=0;j<TotM;j++){
	coefy[i][j]=sec.GetCoefy(i,j);
      };
    };
  };
  if(Ndim>2){
    coefz.resize(grp);
    for(int i=0;i<grp;i++){
      coefz[i].resize(TotM,0.);
      for(int j=0;j<TotM;j++){
	coefz[i][j]=sec.GetCoefz(i,j);
      };
    };
  };

  for(int i=0;i<nmed;i++){
    med[i].CalRemovalCrossSection();
  };


  DetermineSizeCoef();
};

real PLOSSystem::CalFlux(int g,int iter,real epsif)
{
  //int itmax=100000;
  //int itmax=10000;
  int itmax=1000; // Default
  //int itmax=100;

#if rwmc_calculation == 1
  itmax=100000;
#endif
  
  real *Src=new real[TotM];
  real *Fl =new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  real omega=opt.GetOmegai(g);
  if(Ndim==1){
    itmax=1;
    omega=1.;
  };
  bool InConv=SweepInnerIteration(itmax,epsif,g,omega,Fl,Src);
  real errf=RenewFluxInnerIteration(g,Fl);
  if(!InConv){
    bool InConv2=SweepInnerIteration(itmax,epsif,g,1.0,Fl,Src);
    errf=RenewFluxInnerIteration(g,Fl);
    if(!InConv2)cout<<"# Inner iteration is not converged in group "<<g<<"\n";
  };
  delete [] Fl;
  delete [] Src;
  return errf;
};

real PLOSSystem::CalFluxOmega(int g,real epsif)
{
  int itmax=299;
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

  int ym=mi.GetYF();
  int zm=mi.GetZF();
  itmax=2;

  real *tmparray=new real[size_coef];
  CalInnerCoef(g,tmparray);

  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    tt1=0.;    //fluxomega
    int cnt=0; //fluxomega
    int idt=0;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){
        int xer=edge[0][y][z][1];
	int xel=edge[0][y][z][0];
	if(xel!=-1){

	int xmeshnum=xer-xel+1;
	real *a=new real[xmeshnum*3];
	real *b=new real[xmeshnum];
	int tmp=0;
	int id=meshid[z][y][xel];
	for(int x=xel;x<=xer;x++){
	  real flt=Src[id];
          if(y!=edge[1][x][z][0]){flt-=tmparray[idt++]*Fl[meshid[z][y-1][x]];};
          if(y!=edge[1][x][z][1]){flt-=tmparray[idt++]*Fl[meshid[z][y+1][x]];};
          if(z!=edge[2][x][y][0]){flt-=tmparray[idt++]*Fl[meshid[z-1][y][x]];};
          if(z!=edge[2][x][y][1]){flt-=tmparray[idt++]*Fl[meshid[z+1][y][x]];};
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =tmparray[idt++];};
          if(x!=xer){a[tmp*3+2]=tmparray[idt++];};
	  // (for discontinuity factor consideration)
          //if(x!=xel){a[tmp*3]  =tmparray[idt++]*mesh[x-1].GetMed()->GetFlux(1).get_dat(g);};
          //if(x!=xer){a[tmp*3+2]=tmparray[idt++]*mesh[x+1].GetMed()->GetFlux(0).get_dat(g);};
	  a[tmp*3+1]=tmparray[idt++];
	  tmp++;
	  id++;
	};
        gauss_3p(a,b,xmeshnum,1);
	tmp=0;
	id-=xmeshnum;
	for(int x=xel;x<=xer;x++){
	  real flo =Fl[id];
	  real fln =b[tmp];
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
	  tmp++;
	  id++;
	};
	delete [] a;
	delete [] b;

	};
      };
    };
    tt1/=cnt; //fluxomega
    if(cnt==0)tt1=0.; // fluxomega
  };

  delete [] tmparray;

  real errf=RenewFluxInnerIteration(g,Fl);

  // ** Fluxomega **
  if(tt1!=0.){
    real rnew=(omega-1+tt1)/omega/sqrt(tt1);
    rnew=rnew*rnew;
    if(rnew>1.){
      omega=1.;
    }else{
      omega=2.0/(1+sqrt(1-rnew)); 
    };
  }else{
    omega=1.;
  };
  opt.PutOmegai(g,omega);
  //cout<<"   ** omega for inner iteration : "<<omega<<"\n";
  //if(opt.Print())cout<<"   ** omega for inner iteration : "<<omega<<"\n";
  // ***************

  delete [] Fl;
  delete [] Src;
  delete [] Res; 
  return errf;
};

real PLOSSystem::CalFluxModifiedLeakage(int g,real epsif,
		 vector< vector< vector< vector<real> > > > &d1)
{
  real omega=1.;  // fluxmodifiedleakage
  //real omega=1.05;  // fluxmodifiedleakage
  int itmax=999;
  if(Ndim==1)itmax=1;
  real *Src=new real[TotM];
  real *Fl =new real[TotM];
  GetFluxAndSrcInnerIteration(g,Fl,Src);

  int ym=mi.GetYF();
  int zm=mi.GetZF();

  if(Ndim==1)itmax=1;

  bool inner_conv=true;
  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){
	int xer=edge[0][y][z][1];
	int xel=edge[0][y][z][0];
	if(xel!=-1){

	int xmeshnum=xer-xel+1;
	real *a=new real[xmeshnum*3];
	real *b=new real[xmeshnum];
	int tmp=0;
	for(int x=xel;x<=xer;x++){
	  int id=meshid[z][y][x];
	  real flt=Src[id];
          if(y!=edge[1][x][z][0]){flt-=(coefy[g][id]-d1[z][y][x][1])*Fl[meshid[z][y-1][x]];}; 
          if(y!=edge[1][x][z][1]){
            int tmpid=meshid[z][y+1][x];
            flt-=(coefy[g][tmpid]+d1[z][y+1][x][1])*Fl[tmpid];
          };
          if(z!=edge[2][x][y][0]){flt-=(coefz[g][id]-d1[z][y][x][2])*Fl[meshid[z-1][y][x]];};
          if(z!=edge[2][x][y][1]){
            int tmpid=meshid[z+1][y][x];
            flt-=(coefz[g][tmpid]+d1[z+1][y][x][2])*Fl[tmpid];
          };
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =coefx[g][id]-d1[z][y][x][0];};
          if(x!=xer){a[tmp*3+2]=coefx[g][id+1]+d1[z][y][x+1][0];};
	  a[tmp*3+1]=coefs[g][id];
	  tmp++;
	};
	gauss_3p(a,b,xmeshnum,1);
	tmp=0;
	for(int x=xel;x<=xer;x++){
	  int id=meshid[z][y][x];
	  real flo =Fl[id];
	  real fln =b[tmp];
	  if(errmax<epsif){
 	    real err=fabs(fln/flo-1.);
	    if(err>errmax)errmax=err;
	  };
	  Fl[id]=flo+omega*(fln-flo);
	  tmp++;
        };
	delete [] a;
	delete [] b;

	};
      };
    };
    //if(iter>0){cout<<errmax<<" ";};
    if(errmax<epsif&&iter>0){
      //cout<<"Inner iteration : "<<iter<<" "<<errmax<<" "<<epsif<<"\n";
      break;
    };
    if(iter==itmax-1&&itmax!=1){
      //cout<<"Number of inner iteration is exceeded in group "<<g;
      //cout<<" ("<<errmax<<"/"<<epsif<<")\n";
      inner_conv=false;
    };
  };

  real errf=RenewFluxInnerIteration(g,Fl);
  if(!inner_conv)errf*=-1.;

  delete [] Fl;
  delete [] Src;
  return errf;
};

void PLOSSystem::DetermineSizeCoef()
{
  int ym=mi.GetYF();
  int zm=mi.GetZF();
  int idt=0;
  for(int z=0;z<zm;z++){
    for(int y=0;y<ym;y++){
      int xer=edge[0][y][z][1];
      int xel=edge[0][y][z][0];
      if(xel!=-1){
        for(int x=xel;x<=xer;x++){
          if(y!=edge[1][x][z][0])idt++;
          if(y!=edge[1][x][z][1])idt++;
          if(z!=edge[2][x][y][0])idt++;
          if(z!=edge[2][x][y][1])idt++;
          if(x!=xel)      idt++;
          if(x!=xer)      idt++;
          idt++;
	};
      };
    };
  };
  size_coef=idt;
};

void PLOSSystem::CalInnerCoef(int g, real *tmparray)
{
  int ym=mi.GetYF();
  int zm=mi.GetZF();
  int idt=0;
  for(int z=0;z<zm;z++){
    for(int y=0;y<ym;y++){
      int xer=edge[0][y][z][1];
      int xel=edge[0][y][z][0];
      if(xel!=-1){
      int id=meshid[z][y][xel];
      for(int x=xel;x<=xer;x++){

	/*
        real factor=1.;
        real dif_s=mesh[id].GetDif(g,difc[0]);
        real len_s=mesh[id].GetLen(0);
        real dif_l=1.;
        real len_l=0.;
	if(x!=xel){
	  int idl=id-1;
          dif_l=mesh[idl].GetDif(g,difc[0]);
          len_l=mesh[idl].GetLen(0);
	};
        real dif_r=1.;
        real len_r=0.;
	if(x!=xer){
          int idr=id+1;
          dif_r=mesh[idr].GetDif(g,difc[0]);
          len_r=mesh[idr].GetLen(0);
	};
        factor=(4.*len_s/dif_s)/(len_l/dif_l+2*len_s/dif_s+len_r/dif_r);
	*/

        if(y!=edge[1][x][z][0]){tmparray[idt++]=coefy[g][id];};
        if(y!=edge[1][x][z][1]){tmparray[idt++]=coefy[g][meshid[z][y+1][x]];};
        if(z!=edge[2][x][y][0]){tmparray[idt++]=coefz[g][id];};
        if(z!=edge[2][x][y][1]){tmparray[idt++]=coefz[g][meshid[z+1][y][x]];};

        if(x!=xel)      {tmparray[idt++]=coefx[g][id];};
        if(x!=xer)      {tmparray[idt++]=coefx[g][id+1];};
	tmparray[idt++]=coefs[g][id];
	/*
        if(x!=xel)      {tmparray[idt++]=factor*coefx[g][id];};
        if(x!=xer)      {tmparray[idt++]=factor*coefx[g][id+1];};
	tmparray[idt++]=factor*coefs[g][id];
	*/
        id++;
      };
      };
    };
  };
};

bool PLOSSystem::SweepInnerIteration(int itmax,real epsif,int g,real omega,real *Fl,real *Src)
{
  int ym=mi.GetYF();
  int zm=mi.GetZF();

  real *tmparray=new real[size_coef];
  CalInnerCoef(g,tmparray);

  bool conv=true;

  for(int iter=0;iter<itmax;iter++){
    real errmax=0.;
    int idt=0;
    for(int z=0;z<zm;z++){
      for(int y=0;y<ym;y++){
	int xer=edge[0][y][z][1];
	int xel=edge[0][y][z][0];
	if(xel!=-1){
	int xmeshnum=xer-xel+1;
	real *a=new real[xmeshnum*3];
	real *b=new real[xmeshnum];
	int tmp=0;
	int id=meshid[z][y][xel];
	for(int x=xel;x<=xer;x++){
	  real flt=Src[id];
          if(y!=edge[1][x][z][0]){flt-=tmparray[idt++]*Fl[meshid[z][y-1][x]];};
          if(y!=edge[1][x][z][1]){flt-=tmparray[idt++]*Fl[meshid[z][y+1][x]];};
          if(z!=edge[2][x][y][0]){flt-=tmparray[idt++]*Fl[meshid[z-1][y][x]];};
          if(z!=edge[2][x][y][1]){flt-=tmparray[idt++]*Fl[meshid[z+1][y][x]];};
	  b[tmp]=flt;
          if(x!=xel){a[tmp*3]  =tmparray[idt++];};
          if(x!=xer){a[tmp*3+2]=tmparray[idt++];};


#if rwmc_calculation == 1
	  // RWMC-2020&2021 study	  
	  if(y==edge[1][x][z][0]){
	    //real boundary_value=15.; // TRU-repository
	    real boundary_value=36;  // HLW-repository	    
	    real dc=mesh[id].GetMed()->GetMacxs().GetData1d(d).get_dat(g);
	    b[tmp]+=dc/(mi.GetFMeshL(1,0)*0.5)*boundary_value*mesh[id].GetSurL(1);
	  };
	  if(y==edge[1][x][z][1]){
	    //real boundary_value=75.; // TRU-repository
	    real boundary_value=51.; // HLW-repository	    
	    real dc=mesh[id].GetMed()->GetMacxs().GetData1d(d).get_dat(g);	    
	    b[tmp]+=dc/(mi.GetFMeshL(1,mi.GetYF()-1)*0.5)*boundary_value*mesh[id].GetSurR(1);
	  };
#endif

#if heat_conduction_calculation == 1
	  // Heat conduction calculation
	  if(xr_bc>0.&&x==xer){
	    real dc=mesh[id].GetMed()->GetMacxs().GetData1d(d).get_dat(g);
	    b[tmp]+=dc/(mi.GetFMeshL(0,mi.GetXF()-1)*0.5)*xr_bc*mesh[id].GetSurR(0);
	  };
#endif
	  
	  // (for discontinuity factor consideration)
          //if(x!=xel){a[tmp*3]  =tmparray[idt++]*mesh[x-1].GetMed()->GetFlux(1).get_dat(g);};
          //if(x!=xer){a[tmp*3+2]=tmparray[idt++]*mesh[x+1].GetMed()->GetFlux(0).get_dat(g);};
	  a[tmp*3+1]=tmparray[idt++];
	  tmp++;
	  id++;
	};
        gauss_3p(a,b,xmeshnum,1);
	tmp=0;
	id-=xmeshnum;
	for(int x=xel;x<=xer;x++){
	  real flo=Fl[id];
	  real fln=b[tmp];
	  if(errmax<epsif){
 	    real err=fabs(fln/flo-1.);
	    if(err>errmax)errmax=err;
	  };
	  Fl[id]=flo+omega*(fln-flo);
	  tmp++;
	  id++;
	};
	delete [] a;
	delete [] b;
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

  delete [] tmparray;
  return conv;
};

real PLOSSystem::RenewFluxInnerIteration(int g,real *Fl)
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

void PLOSSystem::GetFluxAndSrcInnerIteration(int g,real *Fl,real *Src)
{
  for(int i=0;i<TotM;i++){
    Fl[i] =mesh[i].GetFlux().get_dat(g);
    Src[i]=mesh[i].GetSrcin()*PI4;
  };
};

void PLOSSystem::SetInitialFlux()
{
  SetInitialFlatFlux();
};

real PLOSSystem::CalFluxGeneral(int grp,real epsif,int iter)
{
  real err;
  if(iter==0&&Ndim!=1){
    // This part is modified in 2019/2/26
    // for kinetics calculation 
    // to achieve null-transient condition
    err=CalFluxOmega(grp,0.1);   // Initial guess for omega
    err=CalFluxOmega(grp,0.001); // The second guess for omega
    err=CalFlux(grp,iter,epsif); // Iteration to specific criteria
  }
  else{err=CalFlux(grp,iter,epsif);};

  //real err=CalFlux(grp,iter,epsif);
  return err;
};

void PLOSSystem::DoAcceleration(int iter,real errs,real fiss)
{
  if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+0.5);
  //if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+1.0);
  //if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(1e9);
};

void PLOSSystem::PutCartMeshInfo(CartMeshInfo cm, string geom)
{
  bool cart=false;
  if(geom=="Cartesian")        cart    =true;
  if(geom=="Sphere"&&Ndim==1)  Sphere  =true;
  if(geom=="Cylinder"&&Ndim==2)Cylinder=true;

  if(!cart&&!Sphere&&!Cylinder){
    cout<<"Error in PutCartMeshInfo.\n";
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

  xedgel.resize(yf);
  xedger.resize(yf);
  yedgel.resize(xf);
  yedger.resize(xf);

  int tmp=0;
  int tmp2=0;
  for(int z=0;z<zf;z++){
    for(int y=0;y<yf;y++){
      for(int x=0;x<xf;x++){
	int tm=mi.GetFMat(tmp2);
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
	if(tm>=nmed){
	  cout<<"Error in PLOSSystem::PutCartMeshInfo.\n";
	  cout<<"You requested not-existing medium ID.\n";
          cout<<"Please check Medium ID.\n";
	  exit(0);
	};
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

  SetCoarseMeshInfo();

  for(int i=0;i<Ndim*2;i++){
    if(mi.GetBC(i)==2)BC[i]=2; // Vacuum
    if(mi.GetBC(i)==1)BC[i]=1; // Reflective
    if(mi.GetBC(i)==0)BC[i]=0; // Zero flux
    if(BC[i]<0||BC[i]>2){
      cout<<"Boundary Condition Error in PLOS.\n";
      exit(0);
    };
  };
  if(Ndim!=3){BC[4]=0; BC[5]=0;};
  if(Ndim==1){BC[2]=0; BC[3]=0;};
};

real PLOSSystem::CalCurrent(int x,int y,int z,int g,int dir)
{
  if(dir>=Ndim){
    cout<<"Error in CalCurrent of PLOSSystem.\n";
    cout<<"  Requested direction(0:x,1:y,2:z) is "<<dir<<"\n";
    cout<<"  Calculated dimension is "<<Ndim<<"\n";
    exit(0);
  };

  int xyz[3];
  xyz[0]=x;
  xyz[1]=y;
  xyz[2]=z;

  bool redge=false;
  //if(dir==0&&x==xedger[y]+1)redge=true;
  if(dir==0&&x==edge[0][y][z][1]+1)redge=true;
  //if(dir==1&&y==yedger[x]+1)redge=true;
  if(dir==1&&y==edge[1][x][z][1]+1)redge=true;
  //if(dir==2&&z==mi.GetFMesh(2))redge=true;
  if(dir==2&&z==edge[2][x][y][1]+1)redge=true;

  if(!redge){
    // calculation for left-side current
    int id=meshid[z][y][x];
    real fl=mesh[id].GetFlux().get_dat(g);
    real coefl;
    switch(dir){
      case 0:
        coefl=coefx[g][id];
	break;
      case 1:
        coefl=coefy[g][id];
	break;
      case 2:
        coefl=coefz[g][id];
	break;
      default:
	cout<<"Error in CalCurrent.\n";
	exit(0);
    };
    bool ledge=false;
    //if(dir==0&&x==xedgel[y])ledge=true;
    //if(dir==1&&y==yedgel[x])ledge=true;
    //if(dir==2&&z==0)ledge=true;
    if(dir==0&&x==edge[0][y][z][0])ledge=true;
    if(dir==1&&y==edge[1][x][z][0])ledge=true;
    if(dir==2&&z==edge[2][x][y][0])ledge=true;
    if(!ledge){
      xyz[dir]--;
      return coefl*(fl-mesh[meshid[xyz[2]][xyz[1]][xyz[0]]].GetFlux().get_dat(g));
    };
    switch(BC[dir*2]){
    case 0:
      return CalEdgeCurZeroFlux(id,dir,g,difc[dir])*fl*mesh[id].GetSurL(dir);
    case 1:
      return 0.;
    case 2:
      return CalEdgeCurVacuum(id,dir,g,difc[dir])*fl*mesh[id].GetSurL(dir);
    };
  };

  // calculation for right-edge current
  xyz[dir]--;
  int id=meshid[xyz[2]][xyz[1]][xyz[0]];
  real tmp=mesh[id].GetFlux().get_dat(g)*mesh[id].GetSurR(dir);
  switch(BC[dir*2+1]){
  case 0:
    return -1.*CalEdgeCurZeroFlux(id,dir,g,difc[dir])*tmp;
  case 1:
    return 0.;
  case 2:
    return -1.*CalEdgeCurVacuum(id,dir,g,difc[dir])*tmp;
  };

  cout<<"Error in CalCurrent.\n";
  exit(0);
};

void PLOSSystem::CalCoarseCur(int g, vector< vector< vector< vector<real> > > > &ret, bool cumulative)
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  if(!cumulative){
    for(int i=0;i<zr+1;i++){
      for(int j=0;j<yr+1;j++){
        for(int k=0;k<xr+1;k++){
	  for(int l=0;l<Ndim;l++){
	    ret[i][j][k][l]=0.;
	  };
	};
      };
    };
  };

  int iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi.GetCMeshF(2,z1);z2++){
      int iy=0;
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<mi.GetCMeshF(1,y1);y2++){
          int ix=0;
          for(int x1=0;x1<xr;x1++){
	    for(int x2=0;x2<mi.GetCMeshF(0,x1);x2++){
	      int index=meshid[iz][iy][ix];
	      if(index!=-1){
  	        // X-direction
	        if(x2==0){ret[z1][y1][x1][0]+=CalCurrent(ix,iy,iz,g,0);};
	        //if(ix==xedger[iy]){ret[z1][y1][xr][0]+=CalCurrent(ix+1,iy,iz,g,0);};
	        if(ix==edge[0][iy][iz][1]){ret[z1][y1][xr][0]+=CalCurrent(ix+1,iy,iz,g,0);};
	        // Y-direction
	        if(Ndim>1){
  	          if(y2==0){ret[z1][y1][x1][1]+=CalCurrent(ix,iy,iz,g,1);};
	          //if(iy==yedger[ix]){ret[z1][yr][x1][1]+=CalCurrent(ix,iy+1,iz,g,1);};
	          if(iy==edge[1][ix][iz][1]){ret[z1][yr][x1][1]+=CalCurrent(ix,iy+1,iz,g,1);};
	        };
	        // Z-direction
	        if(Ndim>2){
  	          if(z2==0){ret[z1][y1][x1][2]+=CalCurrent(ix,iy,iz,g,2);};
	          //if(iz==mi.GetZF()-1){ret[zr][y1][x1][2]+=CalCurrent(ix,iy,iz+1,g,2);};
	          if(iz==edge[2][ix][iy][1]){ret[zr][y1][x1][2]+=CalCurrent(ix,iy,iz+1,g,2);};
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

void PLOSSystem::DoCMFDAcceleration(real delk)
{
  int ngrp=1;
  vector<int> bgrp(ngrp);
  /*
  bgrp[0]=0;
  bgrp[1]=1;
  bgrp[2]=18;
  bgrp[3]=27;
  bgrp[4]=36;
  bgrp[5]=45;
  bgrp[6]=69;
  */
  bgrp[0]=grp-1;

  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  // Calculation for fine-mesh current
  vector< vector< vector< vector< vector<real> > > > > CurFF(ngrp); // [g][z+1][y+1][x+1][dim]
  vector< vector< vector< vector< vector<real> > > > > ModD(ngrp);  // [g][z+1][y+1][x+1][dim]

  for(int g=0;g<ngrp;g++){
    CurFF[g].resize(zr+1);
    ModD[g].resize(zr+1);
    for(int i=0;i<zr+1;i++){
      CurFF[g][i].resize(yr+1);
      ModD[g][i].resize(yr+1);
      for(int j=0;j<yr+1;j++){
        CurFF[g][i][j].resize(xr+1);
        ModD[g][i][j].resize(xr+1);
        for(int k=0;k<xr+1;k++){
  	  CurFF[g][i][j][k].resize(Ndim,0.);
	  ModD[g][i][j][k].resize(Ndim,0.);
	};
      };
    };
  };

  int g=0;
  for(int i=0;i<grp;i++){
    if(i>bgrp[g])g++;
    CalCoarseCur(i,CurFF[g],true);
  };

  // Calculation for coarse-group cross sections

  int xyzr=xr*yr*zr;
  PLOSSystem cm(Ndim,ngrp,xyzr);
  if(!print)cm.NoPrint();

  // For multi-group CMFD treatment
  /*
  vector<GroupData1D> fc_siga(TotM); // fine mesh & coarse group
  vector<GroupData1D> fc_nusigf(TotM);
  vector<GroupData1D> fc_kai(TotM);
  vector<GroupData1D> fc_sigt(TotM);
  vector<GroupData1D> fc_d(TotM);
  vector<GroupData1D> fc_flx(TotM);
  vector<GroupData2D> fc_sigs(TotM);
  for(int i=0;i<TotM;i++){
    fc_siga[i].put_imax(ngrp);
    fc_nusigf[i].put_imax(ngrp);
    fc_kai[i].put_imax(ngrp);
    fc_sigt[i].put_imax(ngrp);
    fc_d[i].put_imax(ngrp);
    fc_flx[i].put_imax(ngrp);
    fc_sigs[i].put_yx(ngrp,ngrp);
    GroupData1D flx;
    flx.copy(mesh[i].GetFlux());
    fc_siga[i]=mesh[i].GetMed()->GetMacxs().GetSiga().Cond(flx,ngrp,bgrp);
    fc_nusigf[i]=mesh[i].GetMed()->GetMacxs().GetNusigf().Cond(flx,ngrp,bgrp);
    fc_kai[i]=mesh[i].GetMed()->GetMacxs().GetKai().CondSum(ngrp,bgrp);
    fc_sigt[i]=mesh[i].GetMed()->GetMacxs().GetSigt().Cond(flx,ngrp,bgrp);
    fc_d[i]=mesh[i].GetMed()->GetMacxs().GetD().Cond(flx,ngrp,bgrp);
    fc_flx[i]=flx.CondSum(ngrp,bgrp);
    fc_sigs[i]=mesh[i].GetMed()->GetMacxs().GetSigs().Cond(flx,ngrp,bgrp);
  };

  vector< vector<real> > cc_siga(CTotM); // coarse mesh & coarse mesh
  vector< vector<real> > cc_nusigf(CTotM);
  vector< vector<real> > cc_kai(CTotM);
  vector< vector<real> > cc_sigt(CTotM);
  vector< vector<real> > cc_d(CTotM);
  vector< vector<real> > cc_flx(CTotM);
  vector< vector<real> > cc_sigs(CTotM);
  for(int i=0;i<CTotM;i++){
    cc_siga[i].resize(ngrp,0.);
    cc_nusigf[i].resize(ngrp,0.);
    cc_kai[i].resize(ngrp,0.);
    cc_sigt[i].resize(ngrp,0.);
    cc_d[i].resize(ngrp,0.);
    cc_flx[i].resize(ngrp,0.);
    cc_sigs[i].resize(ngrp*ngrp,0.);
  };

  int iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi.GetCMeshF(2,z1);z2++){
      int iy=0;
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<mi.GetCMeshF(1,y1);y2++){
          int ix=0;
          for(int x1=0;x1<xr;x1++){
	    for(int x2=0;x2<mi.GetCMeshF(0,x1);x2++){
  	      int id=meshid[iz][iy][ix];
	      if(id!=-1){
    	        int index=cmeshid[x1][y1][z1];
		real vol=mesh[id].GetVolume();
		for(int g=0;g<ngrp;g++){
	          real volflx=vol*fc_flx[id].get_dat(g);
	          cc_siga[index][g]+=volflx*fc_siga[id].get_dat(g);
	          real tmp=volflx*fc_nusigf[id].get_dat(g);
	          cc_nusigf[index][g]+=tmp;
	          cc_kai[index][g]+=volflx*fc_kai[id].get_dat(g);
	          cc_d[index][g]+=volflx*fc_d[id].get_dat(g);
  	          cc_sigt[index][g]+=volflx*fc_sigt[id].get_dat(g);
  	          cc_flx[index][g]+=volflx;
		  for(int g2=g;g2<ngrp;g2++){
  		    cc_sigs[index][g*ngrp+g2]+=volflx*fc_sigs[id].get_dat(g,g2);
		  };
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

  Medium minp(ngrp);
  minp.PutPL(0);
  for(int i=0;i<CTotM;i++){
    for(int g=0;g<ngrp;g++){
      real tmp=cc_flx[i][g];
      real fis=cc_nusigf[i][g]/tmp;
      minp.GetMacxs().GetNusigf().put_data(g,fis);
      real sigtfic=cc_sigt[i][g]/tmp-fis/delk;
      minp.GetMacxs().GetSigt().put_data(g,sigtfic);
      minp.GetMacxs().GetSiga().put_data(g,cc_siga[i][g]/tmp-fis/delk);
      real kai=cc_kai[i][g]/tmp;
      minp.GetMacxs().GetKai().put_data(g,kai);
      minp.GetMacxs().GetD().put_data(g,cc_d[i][g]/tmp);
      for(int g2=g;g2<ngrp;g2++){
        minp.GetMacxs().GetSigs().put_data(g,g2,cc_sigs[i][g*ngrp+g2]/tmp);
      };
    };
    cm.AddMedium(minp);
  };
  */

  // For 1-group CMFD treatment
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
    for(int z2=0;z2<mi.GetCMeshF(2,z1);z2++){
      int iy=0;
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<mi.GetCMeshF(1,y1);y2++){
          int ix=0;
          for(int x1=0;x1<xr;x1++){
	    for(int x2=0;x2<mi.GetCMeshF(0,x1);x2++){
  	      int id=meshid[iz][iy][ix];
	      if(id!=-1){
    	        int index=cmeshid[x1][y1][z1];
		real vol=mesh[id].GetVolume();
		volflx=mesh[id].GetFlux()*vol;
		cc_d[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetD();
		cc_sigt[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetSigt();
		if(opt.Forward()){
                  cc_sigs[index]+=volflx*mesh[id].GetMed()->GetMacxs().GetSigs().get_sumx();
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
		cc_flx[index]+=mesh[id].GetVolumeFlux();
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
    xwid[i]=mi.GetCMeshL(0,i);
    fmx[i]=1;
  };
  for(int i=0;i<yr;i++){
    ywid[i]=mi.GetCMeshL(1,i);
    fmy[i]=1;
  };
  for(int i=0;i<zr;i++){
    zwid[i]=mi.GetCMeshL(2,i);
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
    bci[i]=mi.GetBC(i);
  };
  cmi.PutBoundaryCondition(bci);

  string ss="Cartesian";
  if(Cylinder)ss="Cylinder";
  cm.PutCartMeshInfo(cmi,ss);

  // Calculation of initial flux for CMFD 
  vector< vector<real> > FlxF(ngrp);  //[g][m]
  for(int i=0;i<ngrp;i++){
    FlxF[i].resize(CTotM,0.);
  };
  for(int i=0;i<CTotM;i++){
    real inv_vol=1./cm.GetMesh(i).GetVolume();
    for(int g=0;g<ngrp;g++){
      FlxF[g][i]=cc_flx[i]*inv_vol;
      //FlxF[g][i]=cc_flx[i][g]*inv_vol;
    };
  };

  GeneralOption option;
  /*
  real cm_epsf=opt.GetEpsf()*0.1;
  real cm_epsk=opt.GetEpsk()*0.1;
  real cm_epss=opt.GetEpss()*0.1;
  if(cm_epsf>1e-7)cm_epsf=1e-7;
  if(cm_epsk>1e-7)cm_epsk=1e-7;
  if(cm_epss>1e-7)cm_epss=1e-7;
  option.PutEpsf(cm_epsf);
  option.PutEpsk(cm_epsk);
  option.PutEpss(cm_epss);
  */
  real cm_epsf=opt.GetEpsf()*0.1;
  //if(cm_epsf>1e-5)cm_epsf=1e-5;
  if(cm_epsf>1e-7)cm_epsf=1e-7;

  if(!opt.Forward())option.PutAdjointCal();

  option.PutEpsf(cm_epsf);
  /*
  option.PutEpsk(opt.GetEpsk()*0.5);
  option.PutEpss(opt.GetEpss()*0.5);
  */
  //option.PutEpsf(1e-5);
  option.PutEpsk(1e-6);
  option.PutEpss(1e-6);

  cm.PutGeneralOption(option);

  cm.CalCoef();
  cm.PreIgenCoarse(CurFF,FlxF,ModD);

  // Iteration in CMFD
  real k2=cm.CalIgenCoarse(ModD);

  real ke=1./k2+1/delk;
  ke=1./ke;
  ke=ke/k2;

  for(int i=0;i<CTotM;i++){
    real vol=cm.GetMesh(i).GetVolume();
    real tmp=vol*ke;
    for(int g=0;g<ngrp;g++){ // ngrp=1
      cc_flx[i]=cm.GetMesh(i).GetFlux().get_dat(g)*tmp/cc_flx[i];
      //cc_flx[i][g]=cm.GetMesh(i).GetFlux().get_dat(g)*tmp/cc_flx[i][g];
    };
  };

  iz=0;
  for(int z1=0;z1<zr;z1++){
    for(int z2=0;z2<mi.GetCMeshF(2,z1);z2++){
      int iy=0;
      for(int y1=0;y1<yr;y1++){
        for(int y2=0;y2<mi.GetCMeshF(1,y1);y2++){
          int ix=0;
          for(int x1=0;x1<xr;x1++){
  	    int index=cmeshid[x1][y1][z1];
	    for(int x2=0;x2<mi.GetCMeshF(0,x1);x2++){
  	      int id=meshid[iz][iy][ix];
	      if(id!=-1){
		//int g2=0;
		//real tmp=cc_flx[index][g2];
		real tmp=cc_flx[index];
		for(int g=0;g<grp;g++){
		  /*
		  if(g>bgrp[g2]){ // multi-group CMFD
                    g2++;
		    //tmp=cc_flx[index][g2];
		    tmp=cc_flx[index];
		  };
		  */
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

void PLOSSystem::PreIgenCoarse
 (vector< vector< vector< vector < vector<real> > > > > &CurFF,
  vector< vector<real> > &flxf,
  vector< vector< vector< vector < vector<real> > > > > &ModD)
{
  for(int i=0;i<TotM;i++){
    for(int g=0;g<grp;g++){
      mesh[i].GetFlux().put_data(g,flxf[g][i]);
    };
  };
  for(int g=0;g<grp;g++){
    CalCoarseCur(g,ModD[g],false); // ModD <- Current in coarse mesh
  };

  CalCorDif(ModD,CurFF); // ModD <- Diffusion coefficient 
  ModCoef(ModD);
};

void PLOSSystem::CalCorDif
(vector< vector< vector< vector< vector<real> > > > > &curc,
 vector< vector< vector< vector< vector<real> > > > > &curf)
{
  // !! caution !!
  // This method should be only used for coarse-mesh calculation
  // (Coarse-mesh = Fine-Mesh)

  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  for(int z=0;z<zr;z++){
    for(int y=0;y<yr;y++){
      for(int x=0;x<xr;x++){
        int id=meshid[z][y][x];
	if(id!=-1){
	  for(int g=0;g<grp;g++){
            real fl=mesh[id].GetFlux().get_dat(g);
	    for(int a=0;a<Ndim;a++){
	      real tmp=0.;
  	      //if(a==0&&x!=xedgel[y])tmp=mesh[id-1].GetFlux().get_dat(g);
	      //if(a==1&&y!=yedgel[x])tmp=mesh[meshid[z][y-1][x]].GetFlux().get_dat(g);
	      //if(a==2&&z!=0)        tmp=mesh[meshid[z-1][y][x]].GetFlux().get_dat(g);
  	      if(a==0&&x!=edge[0][y][z][0])tmp=mesh[id-1].GetFlux().get_dat(g);
	      if(a==1&&y!=edge[1][x][z][0])tmp=mesh[meshid[z][y-1][x]].GetFlux().get_dat(g);
	      if(a==2&&z!=edge[2][x][y][0])tmp=mesh[meshid[z-1][y][x]].GetFlux().get_dat(g);
  	      curc[g][z][y][x][a]=(curf[g][z][y][x][a]-curc[g][z][y][x][a])/(fl+tmp);
  	      //if(a==0&&x==xedger[y])curc[g][z][y][xr][0]=(curf[g][z][y][xr][0]-curc[g][z][y][xr][0])/fl;
  	      //if(a==1&&y==yedger[x])curc[g][z][yr][x][1]=(curf[g][z][yr][x][1]-curc[g][z][yr][x][1])/fl;
  	      //if(a==2&&z==zr-1)     curc[g][zr][y][x][2]=(curf[g][zr][y][x][2]-curc[g][zr][y][x][2])/fl;
  	      if(a==0&&x==edge[0][y][z][1])curc[g][z][y][xr][0]=(curf[g][z][y][xr][0]-curc[g][z][y][xr][0])/fl;
  	      if(a==1&&y==edge[1][x][z][1])curc[g][z][yr][x][1]=(curf[g][z][yr][x][1]-curc[g][z][yr][x][1])/fl;
  	      if(a==2&&z==edge[2][x][y][1])curc[g][zr][y][x][2]=(curf[g][zr][y][x][2]-curc[g][zr][y][x][2])/fl;
	    };
	  };
	};
      };
    };
  };
};

void PLOSSystem::ModCoef(vector< vector< vector< vector< vector<real> > > > > &moddif)
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  for(int z=0;z<zr;z++){
    for(int y=0;y<yr;y++){
      for(int x=0;x<xr;x++){
	int id=meshid[z][y][x];
	if(id!=-1){
	  for(int g=0;g<grp;g++){
  	    real tmp=-moddif[g][z][y][x][0];
	    //if(x!=xedger[y]){
	    if(x!=edge[0][y][z][1]){
              tmp+=moddif[g][z][y][x+1][0];}
	    else{
	      tmp+=moddif[g][z][y][xr][0];
            };
	    if(Ndim>1){
  	      tmp+=-moddif[g][z][y][x][1];
	      //if(y!=yedger[x]){
	      if(y!=edge[1][x][z][1]){
	        tmp+=moddif[g][z][y+1][x][1];
              }else{
	        tmp+=moddif[g][z][yr][x][1];
	      };
	      if(Ndim>2){
	        tmp+=-moddif[g][z][y][x][2];
	        //if(z!=zr-1){
	        if(z!=edge[2][x][y][1]){
		  tmp+=moddif[g][z+1][y][x][2];}
	        else{
	  	  tmp+=moddif[g][zr][y][x][2];
	        };
	      };
	    };
	    coefs[g][id]+=tmp;
	  };
	};
      };
    };
  };
};

void PLOSSystem::SetDifc(string d1, string d2, string d3)
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

real PLOSSystem::CalEdgeCurZeroFlux(int m,int dir,int g,int flag)
{
  return mesh[m].GetDif(g,flag)/(mesh[m].GetLen(dir)*0.5)*-1.;
};

real PLOSSystem::CalEdgeCurVacuum(int m,int dir,int g,int flag)
{
  real d=mesh[m].GetDif(g,flag);
  return d/(2.1312*d+(mesh[m].GetLen(dir)*0.5))*-1.;
};

real PLOSSystem::CalIgenCoarse(
	vector< vector< vector< vector< vector<real> > > > > &moddif)
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
    bool inner_conv=true;
    for(int g=0;g<grp;g++){
      int ginp=g;
      if(!opt.Forward())ginp=grp-1-g;
      CalSrcMultiplySystem(ginp,fiss,0);
      real err=CalFluxModifiedLeakage(ginp,cin,moddif[ginp]); // o
      //real err=CalFlux(ginp,0,cin); // o
      if(err>errf)errf=err;
      if(err<0.)inner_conv=false;
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
      if(opt.Converged(errf,errk,errs)){
        if(!inner_conv){
	  cout<<"   (Inner iteration is not converged in CMFD calculation.)\n";
	};
        break;
      };
    };
    if(opt.Converged(errf,errk,errs))conv=true;
  };

  return fiss;
}

// **********************************************
// * For perturbation calculation               *
// **********************************************

real PLOSSystem::CalPerturbLeakageTerm(int dir,PLOSSystem *sec,int i,bool *flag)
{
  int xm=mi.GetXF();
  int ym=mi.GetYF();
  int zm=mi.GetZF();

  real ret=0.;
  for(int z=0;z<zm;z++){
    for(int y=0;y<ym;y++){
      for(int x=0;x<xm;x++){
	int m=meshid[z][y][x];
        if(flag[m]&&m!=-1){
          real sfl=mesh[m].GetSurL(dir);
          real sfr=mesh[m].GetSurR(dir);
          real vl,vr;
          if(dir==0&&sphere()){
            real rl=sqrt(sfl/PI4);
            real rr=sqrt(sfr/PI4);
            real rc=(rl+rr)*0.5;
            vr=PI4*0.33333333*pow(rr,3);
            vl=PI4*0.33333333*pow(rc,3);
            real tmp=PI4*0.33333333*pow(rl,3);
            vr-=vl;
            vl-=tmp;
          }else if(dir==0&&cylinder()){
            real z=sec->GetMesh(m).GetLen(1);
            real rl=sfl/(PI2*z);
            real rr=sfr/(PI2*z);
            real rc=(rl+rr)*0.5;
            vr=PI*pow(rr,2)*z;
            vl=PI*pow(rc,2)*z;
            real tmp=PI*pow(rl,2)*z;
            vr-=vl;
            vl-=tmp;
          }else{
            vr=mesh[m].GetVolume()*0.5;
            vl=vr;
          };

          int x1=x;
          int y1=y;
          int z1=z;
          if(dir==0)x1++;
          if(dir==1)y1++;
          if(dir==2)z1++;

          int difc1=GetDifc(dir);
          int difc2=sec->GetDifc(dir);

          real d1=mesh[m].GetDif(i,difc1);
          real d2=sec->GetMesh(m).GetDif(i,difc2);
          real cl1=CalCurrent(x,y,z,i,dir);
          real cl2=sec->CalCurrent(x,y,z,i,dir);
          if(sfl>0.00001){
            cl1/=(d1*sfl); 
            cl2/=(d2*sfl);
          }else{
            cl1=0.;
            cl2=0.;
          };
          real cr1=CalCurrent(x1,y1,z1,i,dir);
          real cr2=sec->CalCurrent(x1,y1,z1,i,dir);
          cr1/=(d1*sfr);
          cr2/=(d2*sfr);
          ret+=-(d2-d1)*(cl1*cl2*vl+cr1*cr2*vr);
	};
      };
    };
  };
  return ret;
};

real PLOSSystem::CalReactivity(PLOSSystem *sec,real kunp,real kp,bool pr)
{
  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };
  real ret=CalReactivity(sec,kunp,kp,flag,pr);
  delete [] flag;
  return ret;
};

real PLOSSystem::CalReactivity(PLOSSystem *sec,real kunp,real kp,bool* flag,bool pr)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in PLOSSystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in PLOSSystem::CalReactivity.\n";
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
  //real ip=1.;

  for(int i=0;i<grp;i++){
    sct[i]=CalPerturbScatteringTermDiffusion(sec,flag,i);
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    abs[i]=CalPerturbAbsorptionTerm(sec,flag,i);
    n2n[i]=CalPerturbN2NTerm(sec,flag,i);
    lx[i]=CalPerturbLeakageTerm(0,sec,i,flag);
    if(Ndim>1)ly[i]=CalPerturbLeakageTerm(1,sec,i,flag);
    if(Ndim>2)lz[i]=CalPerturbLeakageTerm(2,sec,i,flag);
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
    cout<<"#   Component- and group-wise reactivity per unit lathergy x 0.25\n";
    cout<<"#\n";
    cout<<"#   Energy [eV]   ";
    cout<<"Yield        ";
    cout<<"Absorption   ";
    cout<<"Scattering   ";
    cout<<"Leakage      ";
    cout<<"(n,2n)       ";
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
      cout<<i<<" ";
      cout.setf(ios::scientific);
      cout.precision(5);
      real ltmp=lx[i];
      if(Ndim>1)ltmp+=ly[i];
      if(Ndim>2)ltmp+=lz[i];
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real tot=yld[i]+abs[i]+sct[i]+ltmp+n2n[i];
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sct[i]/leth;
      cout<<"  "<<ltmp/leth<<"  "<<n2n[i]/leth<<"  "<<tot/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  real tot=yldsum+abssum+sctsum+n2nsum+lxsum+lysum+lzsum;

  if(pr){
    cout<<"#\n";
    cout<<"# +++ Summary (energy group-integrated value) +++\n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# Yield      : "<<yldsum<<"\n";
    cout<<"# Absorption : "<<abssum<<"\n";
    cout<<"# Scattering : "<<sctsum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";
    cout<<"# (Non-Leak) : "<<yldsum+abssum+sctsum+n2nsum<<"\n";
    cout<<"#\n";
    if(cylinder()||sphere()){
      cout<<"# Leakage-r  : "<<lxsum<<"\n";}
    else{
      cout<<"# Leakage-x  : "<<lxsum<<"\n";
    };
    if(Ndim>1&&cylinder()){
      cout<<"# Leakage-z  : "<<lysum<<"\n";}
    else{
      cout<<"# Leakage-y  : "<<lysum<<"\n";
    };
    if(Ndim>2)cout<<"# Leakage-z  : "<<lzsum<<"\n";
    cout<<"# (Leakage)  : "<<lxsum+lysum+lzsum<<"\n";
    //cout.unsetf(ios::scientific);
    cout<<"# \n";
    cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
    cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
    cout<<"#\n";
    cout<<"#  ( Perturbation denominator : "<<ip<<" )\n";
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

real PLOSSystem::CalReactivitySP3(PLOSSystem *adj_p2,PLOSSystem *fwd_p0,PLOSSystem *fwd_p2,real kunp,real kp,bool pr)
{
  CheckSameMesh(adj_p2);
  CheckSameMesh(fwd_p0);
  CheckSameMesh(fwd_p2);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()||adj_p2->GetGeneralOption().Forward()){
    cout<<"# Error in PLOSSystem::CalReactivitySP3.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!fwd_p0->GetGeneralOption().Forward()||!fwd_p2->GetGeneralOption().Forward()){
    cout<<"# Error in PLOSSystem::CalReactivitySP3.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real *yld=new real[grp];
  real *abs=new real[grp];
  real *sct=new real[grp];
  real *sct_p2=new real[grp];
  real *n2n=new real[grp];
  real *lx =new real[grp];
  real *ly =new real[grp];
  real *lz =new real[grp];

  real ip00=CalPerturbDenominator(fwd_p0);
  real ip20=adj_p2->CalPerturbDenominator(fwd_p0);
  real ip=ip00+ip20;
  //cout<<ip00<<" "<<ip20<<" "<<ip<<"\n";

  for(int i=0;i<grp;i++){
    sct[i]=CalPerturbScatteringTermDiffusion(fwd_p0,flag,i);
    sct_p2[i]=adj_p2->CalPerturbScatteringTermDiffusion(fwd_p0,flag,i);

    yld[i]=CalPerturbYieldTerm(fwd_p0,flag,i,kp);
    yld[i]+=adj_p2->CalPerturbYieldTerm(fwd_p0,flag,i,kp);

    abs[i]=CalPerturbAbsorptionTerm(fwd_p0,flag,i);
    abs[i]+=adj_p2->CalPerturbAbsorptionTerm(fwd_p0,flag,i);
    //abs[i]+=adj_p2->CalPerturbTotalTerm(fwd_p2,flag,i)*-2.5;
    // ... This term is now added to Leakage-x component below

    n2n[i]=CalPerturbN2NTerm(fwd_p0,flag,i);
 
    lx[i]=CalPerturbLeakageTerm(0,fwd_p0,i,flag);
    lx[i]+=CalPerturbLeakageTerm(0,fwd_p2,i,flag)*2.0;
    lx[i]+=adj_p2->CalPerturbLeakageTerm(0,fwd_p2,i,flag)*-27./14.;
    lx[i]+=adj_p2->CalPerturbTotalTerm(fwd_p2,flag,i)*-2.5;
    if(Ndim>1){
      ly[i]=CalPerturbLeakageTerm(1,fwd_p0,i,flag);
      ly[i]+=CalPerturbLeakageTerm(1,fwd_p2,i,flag)*2.0;
      ly[i]+=adj_p2->CalPerturbLeakageTerm(1,fwd_p2,i,flag)*-27./14.;
    };
    if(Ndim>2){
      lz[i]=CalPerturbLeakageTerm(2,fwd_p0,i,flag);
      lz[i]+=CalPerturbLeakageTerm(2,fwd_p2,i,flag)*2.0;
      lz[i]+=adj_p2->CalPerturbLeakageTerm(2,fwd_p2,i,flag)*-27./14.;
    };
  };


  real yldsum=0.;
  real abssum=0.;
  real sctsum=0.;
  real sctp2sum=0.;
  real n2nsum=0.;
  real lxsum=0.;
  real lysum=0.;
  real lzsum=0.;
  real inv_ip=1./ip;
  if(pr){
    cout<<"#   Component- and group-wise reactivity per unit lathergy x 0.25\n";
    cout<<"#\n";
    cout<<"#     - P2-flux components are considered as the leakage component.\n";
    cout<<"#       (Total component is counted as a r- or x-direction leakage.)\n";
    cout<<"#     - P2-adjoint flux component is considered as the scattering component.\n";
    cout<<"#\n";
    cout<<"#   Energy [eV]   ";
    cout<<"Yield        ";
    cout<<"Absorption   ";
    cout<<"Scat(P0-adj) ";
    cout<<"Scat(P2-adj) ";
    cout<<"Leakage      ";
    cout<<"(n,2n)      ";
    cout<<"Total\n";
  };

  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    sct[i]*=inv_ip;
    sct_p2[i]*=inv_ip;
    n2n[i]*=inv_ip;
    lx[i]*=inv_ip;
    ly[i]*=inv_ip;
    lz[i]*=inv_ip;
    yldsum+=yld[i];
    abssum+=abs[i];
    sctsum+=sct[i];
    sctp2sum+=sct_p2[i];
    n2nsum+=n2n[i];
    lxsum+=lx[i];
    if(Ndim>1)lysum+=ly[i];
    if(Ndim>2)lzsum+=lz[i];
    if(pr){
      cout.width(3);
      cout<<i<<" ";
      cout.setf(ios::scientific);
      cout.precision(5);
      real ltmp=lx[i];
      if(Ndim>1)ltmp+=ly[i];
      if(Ndim>2)ltmp+=lz[i];
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real tot=yld[i]+abs[i]+sct[i]+ltmp+n2n[i];
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sct[i]/leth<<"  "<<sct_p2[i]/leth;
      cout<<"  "<<ltmp/leth<<" "<<n2n[i]/leth<<" "<<tot/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  real tot=yldsum+abssum+sctsum+sctp2sum+n2nsum+lxsum+lysum+lzsum;

  if(pr){
    cout<<"#\n";
    cout<<"# +++ Summary (energy group-integrated value) +++\n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# Yield      : "<<yldsum<<"\n";
    cout<<"# Absorption : "<<abssum<<"\n";
    cout<<"# Scattering : "<<sctsum+sctp2sum<<"\n";
    cout<<"#  (P0-adj)  : "<<sctsum<<"\n";
    cout<<"#  (P2-adj)  : "<<sctp2sum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";
    cout<<"# (Non-Leak) : "<<yldsum+abssum+sctsum+sctp2sum+n2nsum<<"\n";
    cout<<"#\n";
    if(cylinder()||sphere()){
      cout<<"# Leakage-r  : "<<lxsum<<"\n";}
    else{
      cout<<"# Leakage-x  : "<<lxsum<<"\n";
    };
    if(Ndim>1&&cylinder()){
      cout<<"# Leakage-z  : "<<lysum<<"\n";}
    else{
      cout<<"# Leakage-y  : "<<lysum<<"\n";
    };
    if(Ndim>2)cout<<"# Leakage-z  : "<<lzsum<<"\n";
    cout<<"# (Leakage)  : "<<lxsum+lysum+lzsum<<"\n";
    //cout.unsetf(ios::scientific);
    cout<<"# \n";
    cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
    cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
    cout<<"#\n";
    cout<<"#  ( Perturbation denominator : "<<ip<<" )\n";    
  };

  delete [] yld;
  delete [] abs;
  delete [] sct;
  delete [] sct_p2;
  delete [] n2n;
  delete [] lx;
  delete [] ly;
  delete [] lz;


  delete [] flag;

  return tot;
};

void PLOSSystem::CalSensitivity(PLOSSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  bool nonl=true;
  bool leak=true;

  for(int i=0;i<Ndim;i++){
    if(GetDifc(i)!=0){
      cout<<"Anisotropic diffusion coefficient ";
      cout<<"cannot be used for sensitivity calculation.\n";
      exit(0);
    };
  };

  if(GetGeneralOption().Forward()){
    cout<<"Error in PLOSSystem::CalSensitivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in PLOSSystem::CalSensitivity.\n";
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
  cout.precision(7);

  CheckSameMesh(sec);

  real ip=CalPerturbDenominator(sec);
  if(ip==0.){
    cout<<"Error in PLOSsystem::CalSensitivity.\n";
    cout<<"Perturbation denominator is zero.\n";
    exit(0);
  };
  // for reaction rate ratio sensitivity
  //ip=1;
  ip=3.67241e+6;
  real inv_ip=1./ip;

  bool *flag=new bool[TotM];

  cout<<grp<<"\n";

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    int nidendf=nid;
    //int nidendf=TranslateNuclideIDFromJFS(nid);

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
        if(nonl){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
          re+=CalPerturbAbsorptionTerm(sec,flag,i);
	};
	if(leak){
          re+=CalPerturbLeakageTerm(0,sec,i,flag);
          if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
          if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
	};
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
        if(nonl){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
          re*=inv_ip;
	};
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
        if(nonl){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
	  re*=inv_ip;
	};
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
      if(nonl)re+=CalPerturbAbsorptionTerm(sec,flag,i);
      if(leak){
        re+=CalPerturbLeakageTerm(0,sec,i,flag);
        if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
        if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
      };
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
    int ndd=2;
    if(nidendf==125)ndd=3;
    cout<<ndd<<"\n";  // `3' means matrix with upscattering
    cout<<nidendf<<"\n";
    cout<<2<<"\n";
    for(int i=0;i<grp;i++){
      cout<<" "<<i<<"\n";
      int st=i;
      if(ndd==3)st=0;
      cout<<" "<<st<<"\n";
      cout<<" "<<grp-1<<"\n";
      for(int k=st;k<grp;k++){
	if(nidendf<1000||k<=i+2){
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
              totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
              sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
	      real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
              sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den[j]*micsigs);
 	      real tmp=0.33333333/totorg[j];
              real s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData2d(sigel,0).get_sumx().get_dat(i);
              real s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData2d(sigel,1).get_sumx().get_dat(i);
              real vmu=s1*0.33333333/s0;
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+(1-vmu)*den[j]*micsigs));
	    };
	  };
          real re=0.;
          if(leak){
            re+=CalPerturbLeakageTerm(0,sec,i,flag);
            if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
            if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
	  };
          if(nonl)re+=CalPerturbScatteringTermDiffusion(sec,flag,i,k);
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
        if(nonl)scat_term[i][k]=CalPerturbScatteringTermDiffusion(sec,flag,i,k);
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
	real lt=0.;
        if(leak){
          lt+=CalPerturbLeakageTerm(0,sec,i,flag);
          if(Ndim>1)lt+=CalPerturbLeakageTerm(1,sec,i,flag);
          if(Ndim>2)lt+=CalPerturbLeakageTerm(2,sec,i,flag);
	};
	if(nonl&&ii==1)lt+=CalPerturbN2NTerm(sec,flag,i);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	    if(ii==1)sec->GetMed(j).GetMacxs().GetSign2n().put_data(i,sigsorg[j]);
          };
	};
        for(int k=i;k<grp;k++){
          real re=lt+scat_term[i][k];
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
      if(leak){
        re+=CalPerturbLeakageTerm(0,sec,i,flag);
        if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
        if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
        re*=inv_ip;
      };
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

SensitivityData PLOSSystem::CalSensitivityNew(PLOSSystem *sec,real keff,int nucnum,int *nucid,bool fiss_matrix)
{
  bool fission_spectrum_matrix=fiss_matrix; // default:false
  if(fission_spectrum_matrix){
    cout<<"#\n#\n# !! Warning : Sensitivity calculation is carried out with fission spectrum matrix.\n#\n#\n";
  };

  for(int i=0;i<Ndim;i++){
    if(GetDifc(i)!=0){
      cout<<"# Anisotropic diffusion coefficient ";
      cout<<"# cannot be used for sensitivity calculation.\n";
      exit(0);
    };
  };

  CheckAdjointForward(sec);
  CheckSameMesh(sec);

  SensitivityData sens;
  sens.PutValue(keff);
  sens.PutGroup(grp);
  sens.GetEnband().copy(med[0].GetEnband());

  GroupData1D sns1d(grp);
  GroupData2D sns2d(grp,grp);

  real *nsforg=new real[nmed];
  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  real *tot1org=new real[nmed];
  real *sigsorg=new real[nmed];
  real *sigs1org=new real[nmed];

  real *fiss_frac=new real[nmed];

  real delta=1.;

  bool *flag=new bool[TotM];

  real ip=0.;
  if(ip_input){
    ip=ip_input_value;
  }else{
    ip=CalPerturbDenominator(sec);
    if(fission_spectrum_matrix){
      ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
    };
    if(ip==0.){
      cout<<"# Error in PLOSsystem::CalSensitivity.\n";
      cout<<"# Perturbation denominator is zero.\n";
      exit(0);
    };
  };

  real factor=keff/ip;

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    cout<<"# Sensitivity calculation for nuclide : "<<nid<<" ("<<nc<<"/"<<nucnum<<")\n";

    bool nuclide_in_system=false;
    for(int i=0;i<TotM;i++){
      if(mesh[i].GetMed()->ExistNuclide(nid)){
        flag[i]=true;
	nuclide_in_system=true;
      }else{
	flag[i]=false;
      };
    };

    if(nuclide_in_system){

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
        for(int i=0;i<grp;i++){
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
              nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
              absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
              totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
  	      real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
  	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
              sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
  	      real tmp=0.33333333/totorg[j];
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den*micsigf));
  	    };
          };
          real re=0.;
          //re+=CalPerturbYieldTerm(sec,flag,i,keff);
          if(!fission_spectrum_matrix){
            re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  }else{
            re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,keff);
	  };
          re+=CalPerturbAbsorptionTerm(sec,flag,i);
          re+=CalPerturbLeakageTerm(0,sec,i,flag);
          if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
          if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
          re*=factor;
          sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
  	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	      sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
              sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
  	    };
          };
        };
        sens.PutSensitivity1D(nid,18,sns1d);

        // Nu
        for(int i=0;i<grp;i++){
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
              nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
  	      real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
  	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micsigf*micnu);
  	    };
          };
          real re=0.;
          re+=CalPerturbYieldTerm(sec,flag,i,keff);
          re*=factor;
	  sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
  	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
  	    };
          };
        };
        sens.PutSensitivity1D(nid,452,sns1d);

        // Chi
        for(int i=0;i<nmed;i++){
	  if(med[i].ExistNuclide(nid)){
  	    real den=med[i].GetNuclide(nid).GetDensity();
  	    real total_fiss=sec->GetIntegratedFlux(i)*sec->GetMed(i).GetMacxs().GetNusigf();
            real part_fiss=den*(sec->GetIntegratedFlux(i)*sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicNusigf());
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
	  if(!fission_spectrum_matrix){
            re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  }else{
            re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,keff);
	  };
          //re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  re*=factor;
	  sns1d.put_data(i,re);
  	  for(int j=0;j<nmed;j++){
	    if(med[j].ExistNuclide(nid)){
	      sec->GetMed(j).GetMacxs().GetKai().put_data(i,totorg[j]);
	    };
	  };
        };
      };
      sens.PutSensitivity1D(nid,181,sns1d);
    }; //(end for fissile nuclide)

    // Capture
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
	  real tmp=0.33333333/totorg[j];
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den*micsigc));
        };
      };
      real re=0.;
      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      re+=CalPerturbLeakageTerm(0,sec,i,flag);
      if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
      if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
      re*=factor;
      sns1d.put_data(i,re);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	};
      };
    };
    sens.PutSensitivity1D(nid,102,sns1d);
    
    // Elastic P0
    sns2d.set_zero();
    for(int i=0;i<grp;i++){
      int st=i;
      for(int k=st;k<grp;k++){
	if(nid<100000||k<=i+2){
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
              totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
              sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
	      real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
              sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
 	      real tmp=0.33333333/totorg[j];
              real s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData2d(sigel,0).get_sumx().get_dat(i);
              real s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData2d(sigel,1).get_sumx().get_dat(i);
              real vmu=s1*0.33333333/s0;
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+(1-vmu)*den*micsigs));
	    };
	  };
          real re=0.;
          re+=CalPerturbLeakageTerm(0,sec,i,flag);
          if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
          if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
          re+=CalPerturbScatteringTermDiffusion(sec,flag,i,k);
          re*=factor;
          sns2d.put_data(i,k,re);
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
	      sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	      sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	    };
	  };
	};
      };
    };
    sens.PutSensitivity2D(nid,2,sns2d);

    // (scattering term)
    vector< vector<real> > scat_term(grp);
    for(int i=0;i<grp;i++){
      scat_term[i].resize(grp,0.);
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
	  };
	};
        scat_term[i][k]=CalPerturbScatteringTermDiffusion(sec,flag,i,k);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	  };
	};
      };
    };

    // Inelastic & n2n P0
    for(int ii=0;ii<2;ii++){
      sns2d.set_zero();
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	    sigsorg[j]=sec->GetMed(j).GetMacxs().GetSign2n().get_dat(i);
            real dtot=delta;
	    if(ii==1)dtot=delta*0.5;
            real tmp=0.33333333/totorg[j];
            sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp+den*dtot));
	    if(ii==1)sec->GetMed(j).GetMacxs().GetSign2n().add_data(i,den*dtot);
	  };
	};
	real lt=0.;
        lt+=CalPerturbLeakageTerm(0,sec,i,flag);
        if(Ndim>1)lt+=CalPerturbLeakageTerm(1,sec,i,flag);
        if(Ndim>2)lt+=CalPerturbLeakageTerm(2,sec,i,flag);
	if(ii==1)lt+=CalPerturbN2NTerm(sec,flag,i);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	    if(ii==1)sec->GetMed(j).GetMacxs().GetSign2n().put_data(i,sigsorg[j]);
          };
	};
        for(int k=i;k<grp;k++){
          real re=lt+scat_term[i][k];
          re*=factor;
          sns2d.put_data(i,k,re);
	};
      };
      int mt=4;
      if(ii==1)mt=16;
      sens.PutSensitivity2D(nid,mt,sns2d);
    };

    // Elastic-p1(mu)
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
	  totorg[j]=sec->GetMed(j).GetMacxs().GetD().get_dat(i);
	  real vsigel=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigel).get_dat(i);
	  real tmp=0.33333333/totorg[j];
	  sec->GetMed(j).GetMacxs().GetD().put_data(i,0.33333333/(tmp-den*delta*vsigel));
        };
      };
      real re=0.;
      re+=CalPerturbLeakageTerm(0,sec,i,flag);
      if(Ndim>1)re+=CalPerturbLeakageTerm(1,sec,i,flag);
      if(Ndim>2)re+=CalPerturbLeakageTerm(2,sec,i,flag);
      re*=factor;
      sns1d.put_data(i,re);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetD().put_data(i,totorg[j]);
	};
      };
    };
    sens.PutSensitivity1D(nid,251,sns1d);

  };

  delete [] fiss_frac;
  delete [] flag;
  delete [] sigsorg;
  delete [] sigs1org;
  delete [] absorg;
  delete [] nsforg;
  delete [] totorg;
  delete [] tot1org;

  return sens;
};

void PLOSSystem::MemoryReductionForPerturbation()
{
  coefs.clear();
  for(int i=0;i<TotM;i++){
    mesh[i].SrcDataClear();
  };
};

void PLOSSystem::PrintFluxDistributionForCylinder()
{
  int xx=mi.GetFMesh(0);
  int yy=mi.GetFMesh(1);
  for(int y=0;y<yy;y++){
    for(int x=0;x<xx;x++){
      int id=meshid[0][y][x];
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

void PLOSSystem::MemoryReductionForHighMomentCal()
{
  for(int i=0;i<nmed;i++){
    med[i].MacxsVectorClear();
  };
  coefs.clear();
  coefx.clear();
  coefy.clear();
  coefz.clear();
  for(int i=0;i<TotM;i++){
    mesh[i].SrcDataClear();
  };
};

void PLOSSystem::PrintFluxDistribution(int zz)
{
  int xx=mi.GetFMesh(0);
  int yy=mi.GetFMesh(1);
  for(int y=0;y<yy;y++){
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

real PLOSSystem::CalSP3(bool rigorous_bc,real alpha)
{
  PLOSSystem sys2;
  return CalSP3(sys2,rigorous_bc,alpha);
};

real PLOSSystem::CalSP3(PLOSSystem &sys2, bool rigorous_bc,real alpha)
{
  //[alpha] in the argument
  //
  // - 0 : diffusion
  // - 2 : SP3 and default

  if(!opt.Forward()){
    return CalSP3Adjoint(sys2,rigorous_bc,alpha);
  };

 
  int oit_max=100; // maximum number of outer iteration
  real epsf_p02=1e-5; // convergence criteria for inner iteration
  real epsf_p2=1e-5; // convergence criteria for inner iteration
  int iit_max=1000; // maximum number of P02 and P2 iteration

  //PLOSSystem sys2(Ndim,grp,nmed); 
  sys2.Init(Ndim,grp,nmed); 
  if(!IsPrintTrue())sys2.NoPrint();

  for(int i=0;i<nmed;i++){
    sys2.AddMedium(med[i]);
    sys2.GetMed(i).GetMacxs().GetData2d(sigs).set_zero();
    // +++ modification for diffusion coefficient
    for(int g=0;g<grp;g++){
      real org=med[i].GetMacxs().GetData1d(d).get_dat(g);
      sys2.GetMed(i).GetMacxs().GetData1d(d).put_data(g,org*27./70.*alpha);
    };

    // +++ special treatment to use (phi0+2phi2) in RHS of the phi2 equation
    for(int g=0;g<grp;g++){
      real sigr0=GetMed(i).GetMacxs().GetData1d(sigr).get_dat(g);
      sys2.GetMed(i).GetMacxs().GetData1d(sigt).add_data(g,sigr0*0.8);
    };

  };

  if(Cylinder){
    sys2.PutCartMeshInfo(mi,"Cylinder");
  }else if(Sphere){
    sys2.PutCartMeshInfo(mi,"Sphere");
  }else{
    sys2.PutCartMeshInfo(mi);
  };

  sys2.PutGeneralOption(opt);
  sys2.CalCoef();

  int totm=TotM;

  SetInitialFlux();
  NormalizeFissionSrc();

  sys2.SetZeroScalarFlux();

  SetZeroScatSrc();
  sys2.SetZeroScatSrc();

  // *** for use of extrapolation
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold;
  int expit=0;
  bool exp=false;
  // ***

  // outer iteration
  vector<real> o_src(totm);
  vector<real> o_src_2(totm);
  vector<real> sr_v_ipi4(totm);

  real fiss=1.;
  real fissold,errk;
  for(int iter=0;iter<oit_max;iter++){

    for(int g=0;g<grp;g++){

      // (fission source)+(out-scattering source)
      CalSrcMultiplySystem(g,fiss,pl);
      for(int m=0;m<totm;m++){
	o_src[m]=mesh[m].GetSrcin(); // (fission source + scattering source) of phi0
        sr_v_ipi4[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigr).get_dat(g)*mesh[m].GetVolume()*INV_PI4;
      };

      for(int iit=0;iit<iit_max;iit++){
        int iter_2=iter;
	if(iit!=0)iter_2=1;
        // +++ phi0+phi2 calculation
        for(int m=0;m<totm;m++){
          real org=o_src[m]; // (fission source + scattering source) of phi0
          real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
          real tmp=alpha*fl2*sr_v_ipi4[m];               // 2*sigr*phi2
	  mesh[m].PutSrcin(org+tmp); // source term for (phi0+phi2)
          sys2.GetMesh(m).PutSrcin(-0.4*org); 
        };

        real err=CalFluxGeneral(g,epsf_p02,iter_2); // (phi0+2phi2) is obtained

        // +++ phi2 calculation
        for(int m=0;m<totm;m++){
          real fl02=mesh[m].GetFlux().get_dat(g); // (phi0+2phi2)
          sys2.GetMesh(m).AddSrcin(0.4*fl02*sr_v_ipi4[m]); // 2/5 x sigr x (phi0+2phi2)
         };
        sys2.CalFluxGeneral(g,epsf_p2,iter_2); // phi2 is obtained.
	if(err<1e-5)break;
      };

      // +++ P2 component is extracted from the first [sys]
      for(int m=0;m<totm;m++){
        real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
        mesh[m].GetFlux().add_data(g,-alpha*fl2);      // [-2phi2] is added to [phi0+2phi2]
      };

      AddDownScatSrc(g,pl);
      SetZeroScatSrc(g);
    };

    // +++ For RIGOROUS vacuum boundary treatment
    //
    // See the notebook in 2010/7/20

    if(rigorous_bc){
      real v21_40=21./40.;
      real v3_40=3./40.;
      real v3_8=3./8.;
    real sfl[3],sfr[3], lenh[3];
    int xyz[3];
    for(int z=0;z<mi.GetZF();z++){
      for(int y=0;y<mi.GetYF();y++){
        for(int x=0;x<mi.GetXF();x++){
	  int i=meshid[z][y][x];
	  if(i!=-1){
  	    xyz[0]=x;
	    xyz[1]=y;
	    xyz[2]=z;
	    GroupData1D sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
            GroupData1D sigr2=sys2.GetMesh(i).GetMed()->GetMacxs().GetSigr();
	    real vol=mesh[i].GetVolume();
	    for(int nd=0;nd<Ndim;nd++){
	      sfl[nd]=mesh[i].GetSurL(nd);
	      sfr[nd]=mesh[i].GetSurR(nd);
  	      lenh[nd]=mesh[i].GetLen(nd)*0.5;
	    };
            for(int g=0;g<grp;g++){
  	      coefs[g][i]=sigr.get_dat(g)*vol;
	      real coefs2=sigr2.get_dat(g)*vol;
              real fl2=sys2.GetMesh(i).GetFlux().get_dat(g);
              real fl02=mesh[i].GetFlux().get_dat(g)+alpha*fl2;
       	      for(int nd=0;nd<Ndim;nd++){
		real dc0=mesh[i].GetDif(g,difc[nd]);
                real dc2=sys2.GetMesh(i).GetDif(g,difc[nd]);
                real length=lenh[nd]; // half_length
                real const1=v21_40*length+dc2;
                real const2=0.5*length+dc0;
	        bool ledge=false;
	        bool redge=false;
	        if(nd==0){
		  if(x==edge[0][y][z][0])ledge=true;
	          if(x==edge[0][y][z][1])redge=true;
		}else if(nd==1){
                  if(y==edge[1][x][z][0])ledge=true;
	          if(y==edge[1][x][z][1])redge=true;
		}else{
                  if(z==edge[2][x][y][0])ledge=true;
	          if(z==edge[2][x][y][1])redge=true;
		};
                if(!ledge){
	          xyz[nd]-=1;
	          int id=meshid[xyz[2]][xyz[1]][xyz[0]];
	          xyz[nd]+=1;
                  real d2=mesh[id].GetDif(g,difc[nd]);
	          real dl=dc0*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          real tmp=dc0*d2/dl*sfl[nd];
	          coefs[g][i]+=tmp;
                  d2=sys2.GetMesh(id).GetDif(g,difc[nd]);
                  dl=dc2*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
                  tmp=dc2*d2/dl*sfl[nd];
                  coefs2+=tmp;
                }else{
	          switch(BC[nd*2]){
	          case 0:
                    coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfl[nd];
                    coefs2-=sys2.CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfl[nd];
                    break;
	          case 2:
		    real coef0=0.5*length+dc0-v3_8*length*v3_40*length/const1;
                    real coef00=v3_8*length/const1*(v3_40)-0.5;
                    real coef02=v3_8*length/const1*(-v21_40)+v3_8;
                    mesh[i].PutScatSrc(g,coef02/coef0*fl2*dc0*sfl[nd]*INV_PI4);
                    real coef2=v21_40*length+dc2-v3_40*length*v3_8*length/const2;
                    real coef20=v3_40*length/const2*(-0.5)+v3_40;
                    real coef22=v3_40*length/const2*(v3_8)-v21_40;
                    sys2.GetMesh(i).PutScatSrc(g,coef20/coef2*fl02*dc2*sfl[nd]*INV_PI4);
                    coefs[g][i]-=dc0/coef0*coef00*sfl[nd];
                    coefs2-=dc2/coef2*coef22*sfl[nd];
                    //coefs[g][i]-=CalEdgeCurVacuum(i,nd,g,difc[nd])*sfl[nd];
                    //coefs2-=sys2.CalEdgeCurVacuum(i,nd,g,difc[nd])*sfl[nd];
		    break;
	          };
                };
	        if(!redge){
	          xyz[nd]+=1;
	          int id=meshid[xyz[2]][xyz[1]][xyz[0]];
	          xyz[nd]-=1;
                  real d2=mesh[id].GetDif(g,difc[nd]);
                  real dl=dc0*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          coefs[g][i]+=dc0*d2/dl*sfr[nd];
                  d2=sys2.GetMesh(id).GetDif(g,difc[nd]);
                  dl=dc2*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          coefs2+=dc2*d2/dl*sfr[nd];
                }else{
	          switch(BC[nd*2+1]){
	          case 0:
                    coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfr[nd];
                    coefs2-=sys2.CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfr[nd];
	  	    break;
	          case 2:
		    real coef0=0.5*length+dc0-v3_8*length*v3_40*length/const1;
                    real coef00=v3_8*length/const1*(v3_40)-0.5;
                    real coef02=v3_8*length/const1*(-v21_40)+v3_8;
                    mesh[i].PutScatSrc(g,coef02/coef0*fl2*dc0*sfr[nd]*INV_PI4);
                    real coef2=v21_40*length+dc2-v3_40*length*v3_8*length/const2;
                    real coef20=v3_40*length/const2*(-0.5)+v3_40;
                    real coef22=v3_40*length/const2*(v3_8)-v21_40;
                    sys2.GetMesh(i).PutScatSrc(g,coef20/coef2*fl02*dc2*sfr[nd]*INV_PI4);
                    coefs[g][i]-=dc0/coef0*coef00*sfr[nd];
                    coefs2-=dc2/coef2*coef22*sfr[nd];

                    //coefs[g][i]-=CalEdgeCurVacuum(i,nd,g,difc[nd])*sfr[nd];
                    //coefs2-=sys2.CalEdgeCurVacuum(i,nd,g,difc[nd])*sfr[nd];
		    break;
		  };
	        };
              };
              sys2.PutCoefs(g,i,coefs2);
	    };
	  };
	};
      };
    };
    };


    // +++ Coefficient modification for vacuum boundary
    //     
    //  ... FAILED
    /*
    real d1[3], sfl[3],sfr[3], lenh[3];
    real d1_2[3];
    int xyz[3];
    for(int z=0;z<mi.GetZF();z++){
      for(int y=0;y<mi.GetYF();y++){
        for(int x=0;x<mi.GetXF();x++){
	  int i=meshid[z][y][x];
	  if(i!=-1){
  	    xyz[0]=x;
	    xyz[1]=y;
	    xyz[2]=z;
	    GroupData1D sigr=mesh[i].GetMed()->GetMacxs().GetSigr();
            GroupData1D sigr2=sys2.GetMesh(i).GetMed()->GetMacxs().GetSigr();
	    real vol=mesh[i].GetVolume();
	    for(int nd=0;nd<Ndim;nd++){
	      sfl[nd]=mesh[i].GetSurL(nd);
	      sfr[nd]=mesh[i].GetSurR(nd);
  	      lenh[nd]=mesh[i].GetLen(nd)*0.5;
	    };
            for(int g=0;g<grp;g++){
  	      coefs[g][i]=sigr.get_dat(g)*vol;
	      real coefs2=sigr2.get_dat(g)*vol;
              for(int nn=0;nn<Ndim;nn++){
  	        d1[nn]=mesh[i].GetDif(g,difc[nn]);
		d1_2[nn]=sys2.GetMesh(i).GetDif(g,difc[nn]);
	      };
       	      for(int nd=0;nd<Ndim;nd++){
	        bool ledge=false;
	        if(nd==0&&x==edge[0][y][z][0])ledge=true;
	        if(nd==1&&y==edge[1][x][z][0])ledge=true;
	        if(nd==2&&z==edge[2][x][y][0])ledge=true;
                if(!ledge){
	          xyz[nd]-=1;
	          int id=meshid[xyz[2]][xyz[1]][xyz[0]];
	          xyz[nd]+=1;
                  real d2=mesh[id].GetDif(g,difc[nd]);
	          real dl=d1[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          real tmp=d1[nd]*d2/dl*sfl[nd];
	          coefs[g][i]+=tmp;
                  d2=sys2.GetMesh(id).GetDif(g,difc[nd]);
                  dl=d1_2[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
                  tmp=d1_2[nd]*d2/dl*sfl[nd];
                  coefs2+=tmp;
                }else{
	          switch(BC[nd*2]){
	          case 0:
                    coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfl[nd];
                    coefs2-=sys2.CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfl[nd];
                    break;
	          case 2:
                    real f0=mesh[i].GetFlux().get_dat(g);
                    real f2=sys2.GetMesh(i).GetFlux().get_dat(g);
                    real rhs=f0*0.5+f2*5./8.;
		    real lhs=f0+f2*2.;
		    real alpha=rhs/lhs;
                    alpha=1./alpha;
                    real d=mesh[i].GetDif(g,difc[nd]);
                    coefs[g][i]-=d/(alpha*d+(mesh[i].GetLen(nd)*0.5))*-1*sfl[nd];
                    rhs=-f0*3./40.+3./8.*f2;
		    lhs=f2;
		    alpha=rhs/lhs;
		    alpha=1./alpha;
                    d=sys2.GetMesh(i).GetDif(g,difc[nd]);
                    coefs2-=d/(alpha*d+(mesh[i].GetLen(nd)*0.5))*-1*sfl[nd];
		    break;
	          };
                };
	        bool redge=false;
	        if(nd==0&&x==edge[0][y][z][1])redge=true;
	        if(nd==1&&y==edge[1][x][z][1])redge=true;
	        if(nd==2&&z==edge[2][x][y][1])redge=true;
	        if(!redge){
	          xyz[nd]+=1;
	          int id=meshid[xyz[2]][xyz[1]][xyz[0]];
	          xyz[nd]-=1;
                  real d2=mesh[id].GetDif(g,difc[nd]);
                  real dl=d1[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          coefs[g][i]+=d1[nd]*d2/dl*sfr[nd];
                  d2=sys2.GetMesh(id).GetDif(g,difc[nd]);
                  dl=d1_2[nd]*mesh[id].GetLen(nd)*0.5+d2*lenh[nd];
	          coefs2+=d1_2[nd]*d2/dl*sfr[nd];
                }else{
	          switch(BC[nd*2+1]){
	          case 0:
                    coefs[g][i]-=CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfr[nd];
                    coefs2-=sys2.CalEdgeCurZeroFlux(i,nd,g,difc[nd])*sfr[nd];
	  	    break;
	          case 2:
                    real f0=mesh[i].GetFlux().get_dat(g);
                    real f2=sys2.GetMesh(i).GetFlux().get_dat(g);
                    real rhs=f0*0.5+f2*5./8.;
		    real lhs=f0+f2*2.;
		    real alpha=rhs/lhs;
                    alpha=1./alpha;
                    real d=mesh[i].GetDif(g,difc[nd]);
                    coefs[g][i]-=d/(alpha*d+(mesh[i].GetLen(nd)*0.5))*-1.*sfl[nd];
                    rhs=-f0*3./40.+3./8.*f2;
		    lhs=f2;
		    alpha=rhs/lhs;
		    alpha=1./alpha;
                    d=sys2.GetMesh(i).GetDif(g,difc[nd]);
                    coefs2-=d/(alpha*d+(mesh[i].GetLen(nd)*0.5))*-1.*sfl[nd];
		    break;
		  };
	        };
              };
              sys2.PutCoefs(g,i,coefs2);
	    };
	  };
	};
      };
    };
    */

    // ** for use of extrapolation
    lamdaold=lamda;
    lamda=SourceAndResidualRevision(fsold,res);
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>5){
      real omega=1./(1.-lamda);
      exp=true;
      WriteSourceExtrapolationInfo(lamda,omega);
      SourceExtrapolation(fsold,omega);
    }else if(exp){
      expit=0;
      exp=false;
    };
    real errs=GetSourceError(fsold);
    // ****

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,0.,errs);

    if(errk<1e-5&&errs<1e-5)break;
  };
 

  for(int i=0;i<nmed;i++){
    sys2.GetMed(i)=med[i];
  };
  sys2.CalCoef();

  for(int m=0;m<totm;m++){
    mesh[m].CalFissionSrc(); // re-calculate fission source
  };

  return fiss;
};

real PLOSSystem::CalSP3Adjoint(bool rigorous_bc,real alpha)
{
  PLOSSystem sys2;
  return CalSP3Adjoint(sys2,rigorous_bc,alpha);
};

real PLOSSystem::CalSP3Adjoint(PLOSSystem &sys2,bool rigorous_bc,real alpha)
{
  //real alpha=2.; // 2 for SP3, 0 for Diffusion

  int oit_max=100; // maximum number of outer iteration
  real epsf_p02=1e-5; // convergence criteria for inner iteration
  real epsf_p2=1e-5; // convergence criteria for inner iteration
  int iit_max=1000; // maximum number of P02 and P2 iteration

  //PLOSSystem sys2(Ndim,grp,nmed); 
  sys2.Init(Ndim,grp,nmed); 
  if(!IsPrintTrue())sys2.NoPrint();

  for(int i=0;i<nmed;i++){
    sys2.AddMedium(med[i]);
    sys2.GetMed(i).GetMacxs().GetData2d(sigs).set_zero();
    // +++ modification for diffusion coefficient
    for(int g=0;g<grp;g++){
      real org=med[i].GetMacxs().GetData1d(d).get_dat(g);
      sys2.GetMed(i).GetMacxs().GetData1d(d).put_data(g,org*(27./35.));
    };
  };

  if(Cylinder){
    sys2.PutCartMeshInfo(mi,"Cylinder");
  }else if(Sphere){
    sys2.PutCartMeshInfo(mi,"Sphere");
  }else{
    sys2.PutCartMeshInfo(mi);
  };

  sys2.PutGeneralOption(opt);
  sys2.CalCoef();

  int totm=TotM;

  SetInitialFlux();
  NormalizeFissionSrc();

  sys2.SetZeroScalarFlux();

  SetZeroScatSrc();
  sys2.SetZeroScatSrc();


  // *** for use of extrapolation
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold;
  int expit=0;
  bool exp=false;
  // ***


  // outer iteration
  vector<real> o_src(totm);
  vector<real> o_src_2(totm);
  vector<real> src_tmp1(totm);
  vector<real> src_tmp2(totm);

  real fiss=1.;
  real fissold,errk;
  for(int iter=0;iter<oit_max;iter++){

    for(int g=grp-1;g>=0;g--){

      // (fission source)+(out-scattering source)
      CalSrcMultiplySystem(g,fiss,pl);
      for(int m=0;m<totm;m++){
	o_src[m]=mesh[m].GetSrcin(); // (fission source + scattering source) of (phi0+phi2)
        src_tmp1[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigr).get_dat(g)*mesh[m].GetVolume()*INV_PI4;
        src_tmp2[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)*mesh[m].GetVolume()*INV_PI4*28./27.;
      };

      for(int iit=0;iit<iit_max;iit++){
        int iter_2=iter;
	if(iit!=0)iter_2=1;

        // +++ phi0+phi2 calculation
        for(int m=0;m<totm;m++){
          real org=o_src[m]; // (fission source + scattering source) of (phi0+phi2)
          real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
          real tmp1=fl2*src_tmp1[m];  // sigr*phi2
	  mesh[m].PutSrcin(org-tmp1); // source term for phi0 equation
        };

        real err=CalFluxGeneral(g,epsf_p02,iter_2); // phi0 is obtained
	//cout<<g<<" "<<iit<<" "<<err<<"\n";

        // +++ phi2 calculation
        for(int m=0;m<totm;m++){	  
          real fl0=mesh[m].GetFlux().get_dat(g); // phi0
          sys2.GetMesh(m).PutSrcin(fl0*src_tmp2[m]); // 28/27 x sigt x phi0
        };
        sys2.CalFluxGeneral(g,epsf_p2,iter_2); // -(phi2 - 28/27 phi0) is calculated.

	//for(int m=0;m<totm;m++){
	//  cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetSrcin()<<" "<<sys2.GetMesh(m).GetSrcin()<<"\n";
	//  cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetFlux().get_dat(g)<<" "<<sys2.GetMesh(m).GetFlux().get_dat(g)<<"\n";
	//};

	for(int m=0;m<totm;m++){
          real fl02=sys2.GetMesh(m).GetFlux().get_dat(g); // -(phi2 - 28/27 phi0)
          real fl0=mesh[m].GetFlux().get_dat(g);          // phi0
          sys2.GetMesh(m).GetFlux().put_data(g,-fl02+28./27.*fl0); // phi2 is overwritten
	};

	//for(int m=0;m<totm;m++){
	for(int m=0;m<1;m++){
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetSrcin()<<" "<<sys2.GetMesh(m).GetSrcin()<<"\n";
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetFlux().get_dat(g)<<" "<<sys2.GetMesh(m).GetFlux().get_dat(g)<<"\n";
	};

	//if(err<1e-7)break;
	if(err<1e-5)break;	
      };
      //exit(0);

      // +++ In the first instance, [P0] is replaced by [P0+P2] to calculate fission and scattering sources
      for(int m=0;m<totm;m++){
	real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
	mesh[m].GetFlux().add_data(g,fl2);             // phi2 is added to phi0 to calculate fission and scattering source of (phi0 + phi2)
      };
      AddDownScatSrc(g,pl);
      SetZeroScatSrc(g);
    };


    // ** for use of extrapolation

    lamdaold=lamda;
    lamda=SourceAndResidualRevision(fsold,res);
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>5){
      real omega=1./(1.-lamda);
      exp=true;
      WriteSourceExtrapolationInfo(lamda,omega);
      SourceExtrapolation(fsold,omega);
    }else if(exp){
      expit=0;
      exp=false;
    };
    real errs=GetSourceError(fsold);

    // ****

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,0.,errs);

    if(errk<1e-5&&errs<1e-5)break;
  };

  // P2 component is extracted from (P0+P2) 

  for(int m=0;m<totm;m++){
    for(int g=0;g<grp;g++){
      real p2=sys2.GetMesh(m).GetFlux().get_dat(g);
      mesh[m].GetFlux().add_data(g,-p2);
    };
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };


  for(int i=0;i<nmed;i++){
    sys2.GetMed(i)=med[i];
  };
  sys2.CalCoef();

  for(int m=0;m<totm;m++){
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };


  return fiss;
};

real PLOSSystem::CalSP3AdjointOld(PLOSSystem &sys2,bool rigorous_bc,real alpha)
{
  //real alpha=2.; // 2 for SP3, 0 for Diffusion
 
  int oit_max=100; // maximum number of outer iteration
  real epsf_p02=1e-5; // convergence criteria for inner iteration
  real epsf_p2=1e-5; // convergence criteria for inner iteration
  int iit_max=1000; // maximum number of P02 and P2 iteration

  //PLOSSystem sys2(Ndim,grp,nmed); 
  sys2.Init(Ndim,grp,nmed); 
  if(!IsPrintTrue())sys2.NoPrint();

  for(int i=0;i<nmed;i++){
    sys2.AddMedium(med[i]);
    sys2.GetMed(i).GetMacxs().GetData2d(sigs).set_zero();
    for(int g=0;g<grp;g++){
      real sigt_org=GetMed(i).GetMacxs().GetData1d(sigt).get_dat(g);
      sys2.GetMed(i).GetMacxs().GetData1d(sigt).put_data(g,sigt_org*35./27.);
    };
  };

  if(Cylinder){
    sys2.PutCartMeshInfo(mi,"Cylinder");
  }else if(Sphere){
    sys2.PutCartMeshInfo(mi,"Sphere");
  }else{
    sys2.PutCartMeshInfo(mi);
  };

  sys2.PutGeneralOption(opt);
  sys2.CalCoef();

  int totm=TotM;

  SetInitialFlux();
  NormalizeFissionSrc();

  sys2.SetZeroScalarFlux();

  SetZeroScatSrc();
  sys2.SetZeroScatSrc();


  // *** for use of extrapolation
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold;
  int expit=0;
  bool exp=false;
  // ***


  // outer iteration
  vector<real> o_src(totm);
  vector<real> o_src_2(totm);
  vector<real> src_tmp1(totm);
  vector<real> src_tmp2(totm);

  real fiss=1.;
  real fissold,errk;
  for(int iter=0;iter<oit_max;iter++){

    for(int g=grp-1;g>=0;g--){

      // (fission source)+(out-scattering source)
      CalSrcMultiplySystem(g,fiss,pl);
      for(int m=0;m<totm;m++){
	o_src[m]=mesh[m].GetSrcin(); // (fission source + scattering source) of phi0
        src_tmp1[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigr).get_dat(g)*mesh[m].GetVolume()*INV_PI4;
        src_tmp2[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)*mesh[m].GetVolume()*INV_PI4*70./27.;
      };

      for(int iit=0;iit<iit_max;iit++){
        int iter_2=iter;
	if(iit!=0)iter_2=1;

        // +++ phi0+phi2 calculation
        for(int m=0;m<totm;m++){
          real org=o_src[m]; // (fission source + scattering source) of phi0
          real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
          real tmp1=fl2*src_tmp1[m];  // sigr*phi2
	  mesh[m].PutSrcin(org-tmp1); // source term for (phi0+phi2)
        };

        real err=CalFluxGeneral(g,epsf_p02,iter_2); // (phi0+2phi2) is obtained
	//cout<<g<<" "<<iit<<" "<<err<<"\n";

        // +++ phi2 calculation
        for(int m=0;m<totm;m++){	  
          real fl0=mesh[m].GetFlux().get_dat(g); // (phi0)
          sys2.GetMesh(m).PutSrcin(fl0*src_tmp2[m]); // 70/27 x sigt x phi0
        };
        sys2.CalFluxGeneral(g,epsf_p2,iter_2); // 2phi0 - 27/14 phi2


	//for(int m=0;m<totm;m++){
	for(int m=0;m<1;m++){
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetSrcin()<<" "<<sys2.GetMesh(m).GetSrcin()<<"\n";
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetFlux().get_dat(g)<<" "<<sys2.GetMesh(m).GetFlux().get_dat(g)<<"\n";
	};


	for(int m=0;m<totm;m++){
          real fl02=sys2.GetMesh(m).GetFlux().get_dat(g); // 2phi0 - 27/14 phi2
          real fl0=mesh[m].GetFlux().get_dat(g); // phi0
          sys2.GetMesh(m).GetFlux().put_data(g,(fl0*2.-fl02)*14./27.); // phi2
	};

	//for(int m=0;m<totm;m++){
	for(int m=0;m<1;m++){
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetSrcin()<<" "<<sys2.GetMesh(m).GetSrcin()<<"\n";
	  //cout<<g<<" "<<iit<<" "<<m<<" "<<mesh[m].GetFlux().get_dat(g)<<" "<<sys2.GetMesh(m).GetFlux().get_dat(g)<<"\n";
	};

	if(err<1e-5)break;
      };
      //exit(0);

      // +++ P0 component is calculated from P0+2 and P2
      for(int m=0;m<totm;m++){
	real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi2
	mesh[m].GetFlux().add_data(g,fl2);      // [phi2] is added to [phi0]
      };
      AddDownScatSrc(g,pl);
      SetZeroScatSrc(g);
    };


    // ** for use of extrapolation

    lamdaold=lamda;
    lamda=SourceAndResidualRevision(fsold,res);
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>5){
      real omega=1./(1.-lamda);
      exp=true;
      WriteSourceExtrapolationInfo(lamda,omega);
      SourceExtrapolation(fsold,omega);
    }else if(exp){
      expit=0;
      exp=false;
    };
    real errs=GetSourceError(fsold);

    // ****

    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,0.,errs);

    if(errk<1e-5&&errs<1e-5)break;
  };

  // P2 component is extracted from (P0+P2) 

  for(int m=0;m<totm;m++){
    for(int g=0;g<grp;g++){
      real p2=sys2.GetMesh(m).GetFlux().get_dat(g);
      mesh[m].GetFlux().add_data(g,-p2);
    };
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };

  /*
  for(int i=0;i<nmed;i++){
    sys2.GetMed(i)=med[i];
  };
  sys2.CalCoef();
  for(int m=0;m<totm;m++){
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };
  */

  return fiss;
};

/*
void PLOSSystem::PutPerturbationDenominator(real val)
{
  ip_input=true;
  ip_input_value=val;
};
*/


//The above SP3 function calculation was developed based on A. Yamamoto's work: a SP3 equations set.
//Now, a new SP3 equations set which was derived from P3 equations set is implemented below. 
real PLOSSystem::CalOSP3Adjoint(bool rigorous_bc)
{
  PLOSSystem sys2;
  return CalOSP3Adjoint(sys2,rigorous_bc);
};
real PLOSSystem::CalOSP3Adjoint(PLOSSystem &sys2,bool rigorous_bc)
{
  int oit_max=100; // maximum number of outer iteration
  real epsf_eq1=1e-5; // convergence criteria for inner iteration
  real epsf_eq2=1e-5; // convergence criteria for inner iteration
  int iit_max=1000; // maximum number of P02 and P2 iteration

  sys2.Init(Ndim,grp,nmed); 
  if(!IsPrintTrue())sys2.NoPrint();

  for(int i=0;i<nmed;i++){
    sys2.AddMedium(med[i]);
    sys2.GetMed(i).GetMacxs().GetData2d(sigs).set_zero();
    // +++ modification for diffusion coefficient
    for(int g=0;g<grp;g++){
      real org=med[i].GetMacxs().GetData1d(d).get_dat(g);
      sys2.GetMed(i).GetMacxs().GetData1d(d).put_data(g,org*(27./35.));
    };
  };
  if(Cylinder){
    sys2.PutCartMeshInfo(mi,"Cylinder");
  }else if(Sphere){
    sys2.PutCartMeshInfo(mi,"Sphere");
  }else{
    sys2.PutCartMeshInfo(mi);
  };

  sys2.PutGeneralOption(opt);
  sys2.CalCoef();

  int totm=TotM;

  SetInitialFlux();
  NormalizeFissionSrc();

  sys2.SetZeroScalarFlux();

  SetZeroScatSrc();
  sys2.SetZeroScatSrc();


  // *** for use of extrapolation
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold;

  int expit=0;
  bool exp=false;
  // ***
  // outer iteration
  vector<real> o_src(totm);
  vector<real> o_src_2(totm);
  vector<real> src_tmp1(totm);
  vector<real> src_tmp2(totm);
  real fiss=1.;
  real fissold,errk;
  for(int iter=0;iter<oit_max;iter++){
    for(int g=grp-1;g>=0;g--){//From high energy gorup to low energy group.
      // (fission source)+(out-scattering source)
      CalSrcMultiplySystem(g,fiss,pl);
      for(int m=0;m<totm;m++){
	    o_src[m]=mesh[m].GetSrcin(); //(fission source + scattering source) at RHS of Eq(1). 
      src_tmp1[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigr).get_dat(g)*mesh[m].GetVolume()*INV_PI4;//sigma_r for eq-1.
      src_tmp2[m]=mesh[m].GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)*mesh[m].GetVolume()*INV_PI4;//sigma_t for eq-2.
      };
      
      for(int iit=0;iit<iit_max;iit++){
        int iter_2=iter;
	      if(iit!=0)iter_2=1;
        // +++ (phi0 + 2/5 phi2) preparation, this is for Eq(1).
        for(int m=0;m<totm;m++){ //this loop prepares the so-called source term for 1st equation.
          real org=o_src[m]; // (fission source + scattering source) 
          real fl2=sys2.GetMesh(m).GetFlux().get_dat(g); // get phi2
          real tmp1=0.4*fl2*src_tmp1[m];  // (2/5) * sigr * phi2
	        mesh[m].PutSrcin(org+tmp1); // source term for (phi0+2/5phi2) equation, which is the eq-1.
        };
        
        real err=CalFluxGeneral(g,epsf_eq1,iter_2); // (phi0+2/5*phi2) is obtained: 1st equation is solved in this line. 
	      //cout<<g<<" "<<iit<<" "<<err<<"\n";

        // +++ Eq(2) calculation
        for(int m=0;m<totm;m++){//this loop prepares the so-called source term for 2nd equation. 
          real fl02_1=mesh[m].GetFlux().get_dat(g);// get the solved (phi0+2/5*phi2).
          sys2.GetMesh(m).PutSrcin(2.*fl02_1*src_tmp2[m]); // 2*(phi0+2/5*phi2)*sigma_t, which is the source term for Eq(2).
        };
        sys2.CalFluxGeneral(g,epsf_eq2,iter_2); // (2phi_0 + 11/7*phi_2) is calculated, and was updated into sys2. 
	      for(int m=0;m<totm;m++){
          real fl02_2=sys2.GetMesh(m).GetFlux().get_dat(g); // (2phi_0 + 11/7*phi_2) in Eq(2).
          real fl02_1=mesh[m].GetFlux().get_dat(g); // read the solved (phi0+2/5*phi2) in the Eq(1).
          sys2.GetMesh(m).GetFlux().put_data(g,(-70./27.)*fl02_1+(35./27.)*fl02_2); // phi2 is saved in instance "sys2".
	      };
	      if(err<1e-5)break;
      };
      
      // +++ phi0 preparing, this should be done after the inner iteration finished.
      for(int m=0;m<totm;m++){
        real fl02_1=mesh[m].GetFlux().get_dat(g);//get the solved (phi0+2/5*phi2).
        real fl02_2=sys2.GetMesh(m).GetFlux().get_dat(g); // phi_2
        mesh[m].GetFlux().put_data(g,fl02_1-(2./5.)*fl02_2);//phi0 is updated, saved in instance "sys".
      };
      AddDownScatSrc(g,pl);
      SetZeroScatSrc(g);
    };

    // ** for use of extrapolation
    lamdaold=lamda;
    lamda=SourceAndResidualRevision(fsold,res);
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>5){
      real omega=1./(1.-lamda);
      exp=true;
      WriteSourceExtrapolationInfo(lamda,omega);
      SourceExtrapolation(fsold,omega);
    }else if(exp){
      expit=0;
      exp=false;
    };
    real errs=GetSourceError(fsold);

    // ****
    fissold=fiss;
    fiss=GetFissionSum();
    errk=fabs(fiss-fissold)/fissold;

    WriteIterationInfo(iter,fiss,errk,0.,errs);

    if(errk<1e-5&&errs<1e-5)break;
  };
  /*
  // P2 component is extracted from (P0+P2) 
  for(int m=0;m<totm;m++){
    for(int g=0;g<grp;g++){
      real fl02_2=sys2.GetMesh(m).GetFlux().get_dat(g);
      mesh[m].GetFlux().add_data(g,-fl02_2);
    };
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };*/
  
  for(int i=0;i<nmed;i++){
    sys2.GetMed(i)=med[i];
  };
  sys2.CalCoef();

  for(int m=0;m<totm;m++){
    mesh[m].CalFissionSrcAdjoint(); // re-calculate adjoint fission source
  };

  return fiss;
};

real PLOSSystem::CalReactivityOSP3(PLOSSystem *adj_p2,PLOSSystem *fwd_p0,PLOSSystem *fwd_p2,real kunp,real kp,bool pr)//*adj_p2 是一个指针，指向包含adjoint phi 2实例/对象的地址
{
  CheckSameMesh(adj_p2);
  CheckSameMesh(fwd_p0);
  CheckSameMesh(fwd_p2);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()||adj_p2->GetGeneralOption().Forward()){
    cout<<"# Error in PLOSSystem::CalReactivityOSP3.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!fwd_p0->GetGeneralOption().Forward()||!fwd_p2->GetGeneralOption().Forward()){
    cout<<"# Error in PLOSSystem::CalReactivityOSP3.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real *yld=new real[grp];
  real *abs=new real[grp];
  real *sct=new real[grp];
  real *sct_p2=new real[grp];
  real *n2n=new real[grp];
  real *lx =new real[grp];
  real *ly =new real[grp];
  real *lz =new real[grp];

  real ip=CalPerturbDenominator(fwd_p0); //<phi+', F'phi'>

  for(int i=0;i<grp;i++){
    //Scattering XS changes relating parts: removal XS + scattering parts. 
    sct[i]=CalPerturbScatteringTermDiffusion(fwd_p0,flag,i);
    //sct_p2[i]=adj_p2->CalPerturbScatteringOnly(fwd_p2,flag,i);//This can be removed directly if the 2-nd order scattering XS is set as zero. 
    
    //Yield component 
    yld[i]=CalPerturbYieldTerm(fwd_p0,flag,i,kp);
  
    //Absorption component
    abs[i]=CalPerturbAbsorptionTerm(fwd_p0,flag,i);//reactivity contributes by Sigma_r. 

    n2n[i]=CalPerturbN2NTerm(fwd_p0,flag,i);
 
    lx[i]=CalPerturbLeakageTerm(0,fwd_p0,i,flag);
    lx[i]+=CalPerturbLeakageTerm(0,fwd_p2,i,flag)*2.0;
    lx[i]+=adj_p2->CalPerturbLeakageTerm(0,fwd_p0,i,flag)*0.4;
    lx[i]+=adj_p2->CalPerturbLeakageTerm(0,fwd_p2,i,flag)*11./7.;

    //Total XS perturbation part
    lx[i]+=adj_p2->CalPerturbTotalTerm(fwd_p2,flag,i);//Due to the total XS part relates to phi-2, it was reagared as leakage component. 
    
    if(Ndim>1){
    ly[i]=CalPerturbLeakageTerm(1,fwd_p0,i,flag);
    ly[i]+=CalPerturbLeakageTerm(1,fwd_p2,i,flag)*2.0;
    ly[i]+=adj_p2->CalPerturbLeakageTerm(1,fwd_p0,i,flag)*0.4;
    ly[i]+=adj_p2->CalPerturbLeakageTerm(1,fwd_p2,i,flag)*11./7.;
    };
    if(Ndim>2){
    lz[i]=CalPerturbLeakageTerm(2,fwd_p0,i,flag);
    lz[i]+=CalPerturbLeakageTerm(2,fwd_p2,i,flag)*2.0;
    lz[i]+=adj_p2->CalPerturbLeakageTerm(2,fwd_p0,i,flag)*0.4;
    lz[i]+=adj_p2->CalPerturbLeakageTerm(2,fwd_p2,i,flag)*11./7.;
    };
  };

  real yldsum=0.;
  real abssum=0.;
  real sctsum=0.;
  real sctp2sum=0.;
  real n2nsum=0.;
  real lxsum=0.;
  real lysum=0.;
  real lzsum=0.;
  real inv_ip=1./ip;
  if(pr){
    cout<<"#   Component- and group-wise reactivity per unit lathergy x 0.25\n";
    cout<<"#\n";
    cout<<"#     - P2-flux components are considered as the leakage component.\n";
    cout<<"#       (Total component is counted as a r- or x-direction leakage.)\n";
    cout<<"#     - P2-adjoint flux component is considered as the scattering component.\n";
    cout<<"#\n";
    cout<<"#   Energy [eV]   ";
    cout<<"Yield        ";
    cout<<"Absorption   ";
    cout<<"Scat(P0-adj) ";
    cout<<"Scat(P2-adj) ";
    cout<<"Leakage      ";
    cout<<"(n,2n)      ";
    cout<<"Total\n";
  };

  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    sct[i]*=inv_ip;
    sct_p2[i]*=inv_ip;
    n2n[i]*=inv_ip;
    lx[i]*=inv_ip;
    ly[i]*=inv_ip;
    lz[i]*=inv_ip;
    yldsum+=yld[i];
    abssum+=abs[i];
    sctsum+=sct[i];
    sctp2sum+=sct_p2[i];
    n2nsum+=n2n[i];
    lxsum+=lx[i];
    if(Ndim>1)lysum+=ly[i];
    if(Ndim>2)lzsum+=lz[i];
    if(pr){
      cout.width(3);
      cout<<i<<" ";
      cout.setf(ios::scientific);
      cout.precision(5);
      real ltmp=lx[i];
      if(Ndim>1)ltmp+=ly[i];
      if(Ndim>2)ltmp+=lz[i];
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real tot=yld[i]+abs[i]+sct[i]+ltmp+n2n[i];
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sct[i]/leth<<"  "<<sct_p2[i]/leth;
      cout<<"  "<<ltmp/leth<<" "<<n2n[i]/leth<<" "<<tot/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  real tot=yldsum+abssum+sctsum+sctp2sum+n2nsum+lxsum+lysum+lzsum;

  if(pr){
    cout<<"#\n";
    cout<<"# +++ Summary (energy group-integrated value) +++\n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# Yield      : "<<yldsum<<"\n";
    cout<<"# Absorption : "<<abssum<<"\n";
    cout<<"# Scattering : "<<sctsum+sctp2sum<<"\n";
    cout<<"#  (P0-adj)  : "<<sctsum<<"\n";
    cout<<"#  (P2-adj)  : "<<sctp2sum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";
    cout<<"# (Non-Leak) : "<<yldsum+abssum+sctsum+sctp2sum+n2nsum<<"\n";
    cout<<"#\n";
    if(cylinder()||sphere()){
      cout<<"# Leakage-r  : "<<lxsum<<"\n";}
    else{
      cout<<"# Leakage-x  : "<<lxsum<<"\n";
    };
    if(Ndim>1&&cylinder()){
      cout<<"# Leakage-z  : "<<lysum<<"\n";}
    else{
      cout<<"# Leakage-y  : "<<lysum<<"\n";
    };
    if(Ndim>2)cout<<"# Leakage-z  : "<<lzsum<<"\n";
    cout<<"# (Leakage)  : "<<lxsum+lysum+lzsum<<"\n";
    //cout.unsetf(ios::scientific);
    cout<<"# \n";
    cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
    cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
    cout<<"#\n";
    cout<<"#  ( Perturbation denominator : "<<ip<<" )\n";    
  };

  delete [] yld;
  delete [] abs;
  delete [] sct;
  delete [] sct_p2;
  delete [] n2n;
  delete [] lx;
  delete [] ly;
  delete [] lz;


  delete [] flag;

  return tot;
};

