#include <cstdlib>
#include "CrossSectionEdit.h"

Medium CrossSectionEdit::HomogenizeAll(GeneralSystem &sys,int cur_id)
{
  int hreg=sys.GetTotM();
  int *nreg=new int[hreg];
  for(int i=0;i<hreg;i++){
    nreg[i]=i;
  };
  Medium ret=Homogenize(sys,hreg,nreg,cur_id);
  delete [] nreg;
  return ret;
};

Medium CrossSectionEdit::Homogenize(GeneralSystem &sys,int hreg, int *nreg,int cur_id)
{
  int grp=sys.GetGrp();
  int pl=sys.GetMed(0).GetPL();
  int pltot=sys.GetMed(0).GetPLT();

  Medium ret(grp);
  ret.PutPL(pl,pltot);
  ret.GetEnband().copy(sys.GetMed(0).GetEnband());

  GroupData1D totf(grp); // total flux
  GroupData1D totc(grp); // total current
  totf.set_zero();
  totc.set_zero();

  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    totf=totf+sys.GetMesh(r).GetFlux()*sys.GetMesh(r).GetVolume();
    for(int g=0;g<grp;g++){
      real tmp=sys.GetMesh(r).GetFlux(cur_id).get_dat(g);
      if(tmp<0.){
	sys.GetMesh(r).GetFlux(cur_id).put_data(g,tmp*-1);
      };
    };
    totc=totc+sys.GetMesh(r).GetFlux(cur_id)*sys.GetMesh(r).GetVolume();
  };

  // 1D flux-weighted cross section
  enum xstype ss[]={siga,nusigf,sigt,sign2n,d,dr,dz,sigtr};
  for(int i=0;i<8;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+sys.GetMesh(r).GetMed()->GetData1D(ss[i]).mult(sys.GetMesh(r).GetFlux())
      	   *sys.GetMesh(r).GetVolume();
    };
    rr=rr/totf;
    ret.GetData1D(ss[i]).copy(rr);
  };

  // 1D current-weighted cross section
  //for(int i=0;i<1;i++){
  //  GroupData1D rr(grp);
  //  enum xstype ss;
  //  if(i==0)ss=sigtr;
  //  rr.set_zero();
  //  for(int j=0;j<hreg;j++){
  //    int r=nreg[j];
  //    rr=rr+mesh[r].GetMed()->GetData1D(ss).mult(cur[r])
  //         *mesh[r].GetVolume();
  //  };
  //  rr=rr/totc;
  //  ret.GetData1D(ss).copy(rr);
  //};

  // total cross section
  for(int i=1;i<=pltot;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+sys.GetMesh(r).GetMed()->GetSigt(i).mult(sys.GetMesh(r).GetFlux(cur_id))
   	   *sys.GetMesh(r).GetVolume();
    };
    rr=rr/totc;
    ret.GetSigt(i).copy(rr);
  };

  // Scattering cross section
  GroupData2D mat(grp,grp);
  for(int l=0;l<=pl;l++){
    mat.set_zero();
    for(int i=0;i<grp;i++){
      for(int j=0;j<grp;j++){
        real bs=0.0;
        for(int k=0;k<hreg;k++){
    	  int r=nreg[k];
          real tmp=sys.GetMesh(r).GetMed()->GetDataSigs(l,i,j)
 	          *sys.GetMesh(r).GetVolume()
	          *sys.GetMesh(r).GetFlux().get_dat(i);
	  bs+=tmp;
        };
        bs/=totf.get_dat(i);
        mat.put_data(i,j,bs);
      };
    };
    ret.GetSigs(l).copy(mat);
  };

  // Fission spectrum
  real *fiss=new real[hreg];
  real totfis=0.0;
  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    fiss[i]=sys.GetMesh(r).GetMed()->GetData1D(nusigf)
           *sys.GetMesh(r).GetFlux()
           *sys.GetMesh(r).GetVolume();
    totfis+=fiss[i];
  };
  if(totfis>0.){
    totfis=1./totfis;
    real kais;
    for(int i=0;i<grp;i++){
      kais=0.0;
      for(int j=0;j<hreg;j++){
        int r=nreg[j];
        kais+=sys.GetMesh(r).GetMed()->GetData1D(chi).get_dat(i)*fiss[j];
      };
      kais*=totfis;
      ret.GetData1D(chi).put_data(i,kais);
    };
  };
  delete []fiss;

  // Fission spectrum vector
  GroupData2D inp_vec(grp,grp);
  inp_vec.set_zero();
  for(int i=0;i<grp;i++){
    for(int j=0;j<grp;j++){
      for(int k=0;k<hreg;k++){
        int r=nreg[k];
        real nsf=sys.GetMesh(r).GetMed()->GetData1D(nusigf).get_dat(i);
  	if(nsf!=0.){
          inp_vec.add_data
	    (i,j,nsf*sys.GetMesh(r).GetVolume()*sys.GetMesh(r).GetFlux().get_dat(i)
	     *sys.GetMesh(r).GetMed()->GetMacxs().GetData2d(chi).get_dat(i,j));
        };
      };
    };
  };
  for(int i=0;i<grp;i++){
    real tot=0.;
    for(int j=0;j<grp;j++){
      tot+=inp_vec.get_dat(i,j);
    };
    if(tot!=0.){
      tot=1./tot;
      for(int j=0;j<grp;j++){
	real org=inp_vec.get_dat(i,j);
	inp_vec.put_data(i,j,org*tot);
      };
    };
  };
  ret.GetMacxs().GetData2d(chi).copy(inp_vec);

  return ret;
}

Medium CrossSectionEdit::HomogenizeCartesian2D(GeneralSystem &sys,int x1,int x2,int y1,int y2)
{
  // Current is assumed to be the same as flux 

  Medium ret;

  int hreg=(x2-x1+1)*(y2-y1+1);
  if(hreg<=0){
    cout<<"Error in CrossSectionEdit::HomogenizeCartesian2D.\n";
    exit(0);
  };
  int *nreg=new int[hreg];

  int tmp=0;
  for(int i=y1;i<=y2;i++){
    for(int j=x1;j<=x2;j++){
      nreg[tmp++]=sys.GetMeshID(j,i,0);
    };
  };

  ret=Homogenize(sys,hreg,nreg,0);
  delete [] nreg;
  return ret;
};
