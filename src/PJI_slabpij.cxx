#include "PJI_slabpij.h"

void PJISlabPij::PutMesh(int i)
{
  mesh=i;
  regid.resize(i,0);
  width.resize(i,0.);
};

void PJISlabPij::PutRegionID(int *inp)
{
  for(int i=0;i<mesh;i++){regid[i]=inp[i];};
};

void PJISlabPij::PutWidth(real *inp)
{
  for(int i=0;i<mesh;i++){
    width[i]=inp[i];
  };
};

void PJISlabPij::CalculationPij(real *xsinp,real *pij,bool aniso,bool black)
{
  vector<real> lmd(mesh);
  vector<real> inv_lmd(mesh);

  real epsf=3.0*2.3025;

  for(int i=0;i<mesh;i++){
    lmd[i]=xsinp[regid[i]]*width[i];
    //lmd[i]=xsinp[i]*width[i];
    inv_lmd[i]=1./lmd[i];
  };

  for(int i=0;i<mesh*mesh;i++){pij[i]=0.;};

  if(!aniso){
    for(int i=0;i<mesh;i++){
      real dist=0.;
      pij[i*mesh+i]+=1.-inv_lmd[i]*(0.5-GetEnx(3,lmd[i]));
      for(int j=i+1;j<mesh;j++){
        real tmp =0.5*(GetEnx(3,dist)-GetEnx(3,dist+lmd[i])-GetEnx(3,dist+lmd[j])+GetEnx(3,dist+lmd[i]+lmd[j]));
        pij[i*mesh+j]+=tmp*inv_lmd[i];
        pij[j*mesh+i]+=tmp*inv_lmd[j];
	dist+=lmd[j];
      };
      if(!black){
        while(dist<epsf){
	  for(int j=0;j<mesh;j++){
            real tmp =0.5*(GetEnx(3,dist)-GetEnx(3,dist+lmd[i])-GetEnx(3,dist+lmd[j])+GetEnx(3,dist+lmd[i]+lmd[j]));
            pij[i*mesh+j]+=tmp*inv_lmd[i];
            pij[j*mesh+i]+=tmp*inv_lmd[j];
	    dist+=lmd[j];
	  };
	};
      };
    };

  }else{
    for(int i=0;i<mesh;i++){
      real dist=0.;
      pij[i*mesh+i]+=1.-3*inv_lmd[i]*(GetEnx(5,0.)-GetEnx(5,lmd[i]));
      for(int j=i+1;j<mesh;j++){
        real tmp2=1.5*(GetEnx(5,dist)-GetEnx(5,dist+lmd[i])-GetEnx(5,dist+lmd[j])+GetEnx(5,dist+lmd[i]+lmd[j]));
        pij[i*mesh+j]+=tmp2*inv_lmd[i];
        pij[j*mesh+i]+=tmp2*inv_lmd[j];
	dist+=lmd[j];
      };
      if(!black){
        while(dist<epsf){
  	  for(int j=0;j<mesh;j++){
            real tmp2=1.5*(GetEnx(5,dist)-GetEnx(5,dist+lmd[i])-GetEnx(5,dist+lmd[j])+GetEnx(5,dist+lmd[i]+lmd[j]));
            pij[i*mesh+j]+=tmp2*inv_lmd[i];
            pij[j*mesh+i]+=tmp2*inv_lmd[j];
	    dist+=lmd[j];
	  };
	};
      };
    };
  };

  if(!black){
    for(int i=0;i<mesh;i++){
      real tmp=0.;
      for(int j=0;j<mesh;j++){
        tmp+=pij[i*mesh+j];
      };
      tmp=1./tmp;
      for(int j=0;j<mesh;j++){
        pij[i*mesh+j]*=tmp;
      };
    };
  };
};
