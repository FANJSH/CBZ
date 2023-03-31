#include <cstdlib>
#include "PJI_system.h"

using namespace std;

PJISystem::PJISystem(int i,int j):GeneralSystem(2,i,j)
{
  pij.resize(grp);
  pijr.resize(grp);
  pijz.resize(grp);
  name="PJI";
  slab=false;
  write_proc=true;
  PutLeakageTreatment("PseudoAbsorption");
  PutSigmaCol(sigt);
  b2=0.;
  b2r=0.;
  b2z=0.;
  put_pijr=false;
  pij_transform.resize(grp);
  for(int k=0;k<grp;k++){
    pij_transform[k]=false;
  };
};

void PJISystem::PutLeakageTreatment(string inp)
{
  if(inp=="PseudoAbsorption"||
     inp=="ModifiedSource"||
     inp=="ModifiedSelfScattering"||
     inp=="Tibere"){LeakTreat=inp;}
  else{
    cout<<"Error in PutLeakageTreatment.\n";
    cout<<"Inappropriate option.\n";
    cout<<"Your option is "<<inp<<"\n";
  };
};

void PJISystem::PutSigmaCol(enum xstype inp)
{
  if(inp==sigt||inp==sigtr){SigmaCol=inp;}
  else{
    cout<<"Error in PutSigmaCol.\n";
    cout<<"Inappropriate option.\n";
    cout<<"Your option is "<<inp<<"\n";
  };
};

void PJISystem::PutTrajectorySet(TrajectorySet *cinp)
{
  cpij=cinp;
  TotM=cpij->GetRegnum();
  if(print)cout<<"#*** Total Mesh : "<<TotM<<"\n";
  medium_id_par_region.resize(TotM);
  mesh.resize(TotM);
  for(int i=0;i<TotM;i++){
    mesh[i].PutVolume(cpij->GetVol(i));
    mesh[i].PutDim(2);
    mesh[i].PutPL(1);
    mesh[i].PutGrp(grp);
  };
};

void PJISystem::PutCartMeshInfo(CartMeshInfo cm)
{
  mi=cm;
  TotM=cm.GetFMesh(0);
  if(print)cout<<"#*** Total Mesh : "<<TotM<<"\n";
  mesh.resize(TotM);
  medium_id_par_region.resize(TotM);
  for(int i=0;i<TotM;i++){
    mesh[i].PutVolume(cm.GetFMeshL(0,i));
    mesh[i].PutDim(2);
    mesh[i].PutPL(1);
    mesh[i].PutGrp(grp);
    mesh[i].PutMedium(&med[cm.GetFMat(i)]);
    medium_id_par_region[i]=cm.GetFMat(i);
  };
  slab=true;

  slabpij.PutMesh(TotM);
  int *regid=new int[TotM];
  real *wid=new real[TotM];
  for(int i=0;i<TotM;i++){
    regid[i]=mi.GetFMat(i);
    wid[i]=mi.GetFMeshL(0,i);
  };
  slabpij.PutRegionID(regid);
  slabpij.PutWidth(wid);
  delete [] regid;
  delete [] wid;

};

real PJISystem::CalHomoSigt(int g,enum xstype ss)
{
  real hsigtr=0.;
  real base=0.;
  for(int i=0;i<TotM;i++){
    hsigtr+=mesh[i].GetMed()->GetData1D(g,ss)
           *mesh[i].GetFlux().get_dat(g)*mesh[i].GetVolume();
    base+=mesh[i].GetFlux().get_dat(g)*mesh[i].GetVolume();
  };
  if(base==0.){
    cout<<"# Error in PJISystem::CalHomoSigt.\n";
    cout<<"# Neutron flux is not yet calculated.\n";
    exit(0);
  };
  return hsigtr/base;
};

void PJISystem::PutPij(enum xstype ss,real b2)
{
  real *xs=new real[TotM];
  real *pijinp=new real[TotM*TotM];
  for(int ig=0;ig<grp;ig++){
    //cout<<"Group : "<<ig<<"\n"; // cbg
    pij[ig].put_yx(TotM,TotM);
    real hsigtr=0.;
    if(b2!=0.)hsigtr=b2*0.33333333/CalHomoSigt(ig,ss);
    //cout<<hsigtr<<"\n";
    for(int i=0;i<TotM;i++){
      xs[i]=mesh[i].GetMed()->GetData1D(ig,ss)+hsigtr;
      //cout<<xs[i]<<" "; // cbg
      if(xs[i]<0.){
	cout<<"Error in PJISystem::PutPij.\n";
	cout<<"Negative cross section is detected.\n";
	exit(0);
      };
    };
    //cout<<"\n"; // cbg
    if(slab){
      slabpij.CalculationPij(xs,pijinp,false);
    }else{
      cpij->CalculationPij(xs,pijinp,false);
    };
    pij[ig].put_data(pijinp);
    //if(write_proc)cout<<"# Terminated Pij Cal. in group "<<ig<<"\n";
    if(print)cout<<"# Terminated Pij Cal. in group "<<ig<<"\n";
  };

  PijTransformFalse();
  delete [] xs;
  delete [] pijinp;
}

void PJISystem::PutPijr(enum xstype ss,real b2)
{
  put_pijr=true;

  real *xs=new real[TotM];
  real *pijinp=new real[TotM*TotM];
  for(int ig=0;ig<grp;ig++){
    pijr[ig].put_yx(TotM,TotM);
    pijz[ig].put_yx(TotM,TotM);
    real hsigtr=0.;
    if(b2!=0.)hsigtr=b2*0.33333333/CalHomoSigt(ig,ss);
    for(int i=0;i<TotM;i++){
      xs[i]=mesh[i].GetMed()->GetData1D(ig,ss)+hsigtr;
    };
    if(slab){
      slabpij.CalculationPij(xs,pijinp,true);
    }else{
      cpij->CalculationPij(xs,pijinp,true);
    };
    pijr[ig].put_data(pijinp); // perpendicular
    if(slab){
      pijz[ig].copy((pij[ig]*3.-pijr[ig])*0.5); // parallel
    }else{
      pijz[ig].copy(pij[ig]*3.-pijr[ig]*2.); // parallel
    };

    //if(write_proc)cout<<"# Terminated Pijk Cal. in group "<<ig<<"\n";
    if(print)cout<<"# Terminated Pijk Cal. in group "<<ig<<"\n";
  };
  delete [] xs;
  delete [] pijinp;
}

void PJISystem::PutPijSphere()
{
  real *xs=new real[TotM];
  real *pijinp=new real[TotM*TotM];
  for(int ig=0;ig<grp;ig++){
    pij[ig].put_yx(TotM,TotM);
    for(int i=0;i<TotM;i++){
      xs[i]=mesh[i].GetMed()->GetData1D(ig,SigmaCol);
      if(xs[i]<0.){
	cout<<"Error in PJISystem::PutPij.\n";
	cout<<"Negative cross section is detected.\n";
	exit(0);
      };
    };
    cpij->CalculationPijSphere(xs,pijinp);
    pij[ig].put_data(pijinp);
    //if(write_proc)cout<<"# Terminated Pij Cal. in group "<<ig<<"\n";
    if(print)cout<<"# Terminated Pij Cal. in group "<<ig<<"\n";
  };

  delete [] xs;
  delete [] pijinp;
}


real PJISystem::CalFluxGeneral(int ng,real cin,int iter)
{
  real *flx=new real[TotM];
  for(int i=0;i<TotM;i++){
    flx[i]=mesh[i].GetFlux().get_dat(ng);
  };

  if(LeakTreat=="Tibere"){CalFluxTibere(ng);}
  else if(LeakTreat=="PseudoAbsorption"){CalFluxPijNew(ng,b2,iter);}
  //else if(LeakTreat=="PseudoAbsorption"){CalFluxPij(ng,b2);}
  else if(LeakTreat=="ModifiedSelfScattering"){CalFluxPij(ng,0.);}
  else if(LeakTreat=="ModifiedSource"){CalFluxPijModifiedSource(ng);}
  else{
    cout<<"Error in CalFlux.\n";
    cout<<LeakTreat<<"\n";
    exit(0);
  };

  real errmax=0.;
  for(int i=0;i<TotM;i++){
    real flnew=mesh[i].GetFlux().get_dat(ng);
    real err=fabs(flnew-flx[i])/flx[i];
    if(err>errmax)errmax=err;
  };
  delete [] flx;

  return errmax;
};

void PJISystem::CalFluxTibere(int ng)
{
  GroupData2D Atemp(TotM*3,TotM*3);
  GroupData1D Btemp(TotM*3);
  Atemp.set_zero();
  Btemp.set_zero();

  for(int i1=0;i1<TotM;i1++){
    real a=0.;
    real b=0.;
    real c=0.;
    for(int i2=0;i2<TotM;i2++){
      a+=pij[ng].get_dat(i2,i1)*mesh[i2].GetSrcin(0)*PI4;
      b+=pijr[ng].get_dat(i2,i1)*mesh[i2].GetSrcin(1)*PI4;
      c+=pijz[ng].get_dat(i2,i1)*mesh[i2].GetSrcin(2)*PI4;
    };
    Btemp.put_data(i1,      a);
    Btemp.put_data(TotM+i1,  b);
    Btemp.put_data(TotM*2+i1,c);

    for(int i2=0;i2<TotM;i2++){
      real diagp0=mesh[i2].GetMed()->GetDataSigs(0,ng,ng);
      real diagp1=mesh[i2].GetMed()->GetDataSigs(1,ng,ng);
      Atemp.put_data(i1,      i2,      -pij[ng].get_dat(i2,i1)*diagp0);
      Atemp.put_data(i1,      TotM+i2,   b2r*pijr[ng].get_dat(i2,i1));
      Atemp.put_data(i1,      TotM*2+i2, b2z*pijz[ng].get_dat(i2,i1));
      Atemp.put_data(TotM+i1,  i2,      -pijr[ng].get_dat(i2,i1));
      Atemp.put_data(TotM*2+i1,i2,      -pijz[ng].get_dat(i2,i1));
      Atemp.put_data(TotM+i1,  TotM+i2,  -pijr[ng].get_dat(i2,i1)*diagp1);
      Atemp.put_data(TotM*2+i1,TotM*2+i2,-pijz[ng].get_dat(i2,i1)*diagp1);
    };
    real tot=mesh[i1].GetMed()->GetDataSigt(0,ng);
    Atemp.add_data(i1,      i1,      tot);
    Atemp.add_data(TotM+i1,  TotM+i1,  3.0*tot);
    Atemp.add_data(TotM*2+i1,TotM*2+i1,3.0*tot);
  };
  Atemp.solveaxb_mod(Btemp);

  for(int i=0;i<TotM;i++){
    real vv=mesh[i].GetVolume();
    real flnew=Btemp.get_dat(i)/vv;
    mesh[i].GetFlux(0).put_data(ng,flnew);
    mesh[i].GetFlux(1).put_data(ng,Btemp.get_dat(TotM+i)/vv); // perpendicular
    mesh[i].GetFlux(2).put_data(ng,Btemp.get_dat(TotM*2+i)/vv); // parallel
  };
};

void PJISystem::PijTransformTrue()
{
  for(int i=0;i<grp;i++){
    pij_transform[i]=true;
  };
};

void PJISystem::PijTransformFalse()
{
  for(int i=0;i<grp;i++){
    pij_transform[i]=false;
  };
};

void PJISystem::CalFluxPijNew(int i,real b2,int iter)
{
  real hsigtr=0.;
  if(b2!=0.)hsigtr=b2*0.33333333/CalHomoSigt(i,SigmaCol);

  //if(iter==0&&!pij_transform[i]){
  if(!pij_transform[i]){
    pij_transform[i]=true;
    GroupData2D Amat(TotM,TotM);
    GroupData2D Smat(TotM,TotM);
    for(int j=0;j<TotM;j++){
      for(int k=0;k<TotM;k++){
        real sscat=mesh[k].GetMed()->GetDataSigs(0,i,i)
  	           -mesh[k].GetMed()->GetDataSigt(0,i)
                   +mesh[k].GetMed()->GetData1D(i,SigmaCol);
        real tmp=-pij[i].get_dat(k,j)*sscat;
        if(j==k)tmp+=(mesh[k].GetMed()->GetData1D(i,SigmaCol)+hsigtr);
        Amat.put_data(j,k,tmp);
      };
      real tmp=0.;
      for(int k=0;k<TotM;k++){
        Smat.put_data(j,k,pij[i].get_dat(k,j));
      };
    };
    //Amat.show_self();
    //Smat.show_self();
    /*
    Amat.DoInverse();
    Amat=Amat*Smat;
    */
    Amat.solveaxb_mod(Smat);
    for(int j=0;j<TotM;j++){
      for(int k=0;k<TotM;k++){
	pij[i].put_data(j,k,Smat.get_dat(j,k));
      };
    };
  };

  GroupData1D Bmat(TotM);
  bool zero_src=true;
  for(int j=0;j<TotM;j++){
    real tmp=mesh[j].GetSrcin(0)*PI4;
    Bmat.put_data(j,tmp);
    if(tmp!=0.)zero_src=false;
  };

  if(zero_src){
    for(int j=0;j<TotM;j++){
      mesh[j].GetFlux().put_data(i,0.);
    };
  }else{
    GroupData1D sol=pij[i]*Bmat;
    for(int j=0;j<TotM;j++){
      real flnew=sol.get_dat(j)/mesh[j].GetVolume();
      if(flnew<0.){
	//if(flnew>-1e-9){
	if(flnew>-1e-20){
	  flnew=0.;
	}else{
          cout<<"# Warning in PJISystem::CalFluxPij.\n";
	  /*
          cout<<"# Error in PJISystem::CalFluxPij.\n";
          cout<<"# Negative flux ["<<flnew<<"] is detected in group "<<i<<" and mesh "<<j<<"\n";
          exit(0);
	  */
          cout<<"# Warning in PJISystem::CalFluxPij.\n";
          cout<<"# Negative flux ["<<flnew<<"] is detected in group "<<i<<" and mesh "<<j<<"\n";
	  cout<<"# Reset to zero.\n";
	  flnew=0.;
	/*
	cout<<"#\n# Collision probability matrix is shown below:\n";
	pij[i].show_self();
	for(int ii=0;ii<TotM;ii++){
	  for(int jj=0;jj<TotM;jj++){
	    if(pij[i].get_dat(ii,jj)<0.){
	      cout<<"#   Negative Pij is detected from "<<ii<<" to "<<jj<<" : "<<pij[i].get_dat(ii,jj)<<"\n";
	    };
	  };
	};
	cout<<"#\n# Collision cross section is shown below:\n";
	for(int k=0;k<TotM;k++){
	  cout<<k<<" "<<mesh[k].GetMed()->GetData1D(i,SigmaCol)<<"\n";
	};
	*/

	};
      };
      //if(flnew<0.)flnew=1e-20;
      mesh[j].GetFlux().put_data(i,flnew);
    };
  };

};

void PJISystem::CalFluxPij(int i,real b2)
{
  real hsigtr=0.;
  if(b2!=0.)hsigtr=b2*0.33333333/CalHomoSigt(i,SigmaCol);

  vector<real> Amat(TotM*TotM);
  vector<real> Bmat(TotM);

  bool zero_src=true;
  for(int j=0;j<TotM;j++){
    for(int k=0;k<TotM;k++){
      real sscat=mesh[k].GetMed()->GetDataSigs(0,i,i)
	         -mesh[k].GetMed()->GetDataSigt(0,i)
                 +mesh[k].GetMed()->GetData1D(i,SigmaCol);
      real tmp=-pij[i].get_dat(k,j)*sscat;
      if(j==k)tmp+=(mesh[k].GetMed()->GetData1D(i,SigmaCol)+hsigtr);
      Amat[j*TotM+k]=tmp;
    };
    real tmp=0.;
    for(int k=0;k<TotM;k++){
      tmp+=pij[i].get_dat(k,j)*mesh[k].GetSrcin(0)*PI4;
    };
    Bmat[j]=tmp;
    if(tmp!=0.)zero_src=false;
  };

  if(zero_src){
    for(int j=0;j<TotM;j++){
      mesh[j].GetFlux().put_data(i,0.);
    };
  }else{

    gauss(Amat,Bmat,TotM,1);

    for(int j=0;j<TotM;j++){
      real flnew=Bmat[j]/mesh[j].GetVolume();
      if(flnew<0.){
        cout<<"# Error in PJISystem::CalFluxPij.\n";
        cout<<"# Negative flux is detected in group "<<i<<"\n";
	cout<<"#\n# Collision probability matrix is shown below:\n";
	pij[i].show_self();
	for(int ii=0;ii<TotM;ii++){
	  for(int jj=0;jj<TotM;jj++){
	    if(pij[i].get_dat(ii,jj)<0.){
	      cout<<"#   Negative Pij is detected from "<<ii<<" to "<<jj<<" : "<<pij[i].get_dat(ii,jj)<<"\n";
	    };
	  };
	};
	cout<<"#\n# Collision cross section is shown below:\n";
	for(int k=0;k<TotM;k++){
	  cout<<k<<" "<<mesh[k].GetMed()->GetData1D(i,SigmaCol)<<"\n";
	};
        exit(0);
      };
      //if(flnew<0.)flnew=1e-20;
      mesh[j].GetFlux().put_data(i,flnew);
    };
  };

  //delete [] Amat;
  //delete [] Bmat;
};

void PJISystem::CalFluxPijModifiedSource(int i)
{
  real sig=CalHomoSigt(i,sigtr);
  real sig1=CalHomoSigt(i,SigmaCol);
  //real mod=sig1/(sig1+0.33333333/sig*b2);
  real mod=sig1/(sig1+0.33333333/sig1*b2);

  /*
  GroupData2D Amat(TotM,TotM);
  GroupData1D Bmat(TotM);
  GroupData1D Cmat(TotM);
  */
  vector<real> Amat(TotM*TotM);
  vector<real> Bmat(TotM);

  for(int j=0;j<TotM;j++){
    for(int k=0;k<TotM;k++){
      real sscat=mesh[k].GetMed()->GetDataSigs(0,i,i)
                 -mesh[k].GetMed()->GetDataSigt(0,i)
                 +mesh[k].GetMed()->GetData1D(i,SigmaCol);
      sscat*=mod;
      real tmp=-pij[i].get_dat(k,j)*sscat;
      if(j==k)tmp+=mesh[k].GetMed()->GetData1D(i,SigmaCol);
      //Amat.put_data(j,k,tmp);
      Amat[j*TotM+k]=tmp;
    };
    real tmp=0.;
    for(int k=0;k<TotM;k++){
      tmp+=pij[i].get_dat(k,j)*mesh[k].GetSrcin(0)*mod*PI4;
    };
    //Bmat.put_data(j,tmp);
    Bmat[j]=tmp;
  };

  //Cmat=Amat.solveaxb(Bmat);
  gauss(Amat,Bmat,TotM,1);

  for(int j=0;j<TotM;j++){
    //real flnew=Cmat.get_dat(j)/mesh[j].GetVolume();
    //real flnew=Bmat.get_dat(j)/mesh[j].GetVolume();
    real flnew=Bmat[j]/mesh[j].GetVolume();
    mesh[j].GetFlux().put_data(i,flnew);
  };
};

void PJISystem::SetInitialFlux()
{
  SetInitialFlatFlux();
};

real PJISystem::BucklingSearch(string inp,enum xstype stinp)
{
  PutLeakageTreatment(inp);
  PutSigmaCol(stinp);

  cout<<"*************************************************\n";
  cout<<"* Buckling search with "<<LeakTreat<<" model    *\n";

  real critb;
  if(inp=="Tibere"){critb=BucklingSearchTibere();}
  else if(inp=="ModifiedSelfScattering"){critb=BucklingSearchMSS();}
  else{critb=BucklingSearchPij();};

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"*************************************************\n";
  cout<<"* Buckling search end\n";
  cout<<"*    Critical square buckling : "<<critb<<"\n";
  cout<<"*************************************************\n";

  if(LeakTreat!="Tibere"){
    if(LeakTreat=="PseudoAbsorption"){
      PutApproximatedCurrent(b2);
    }else{
      PutApproximatedCurrent(0.);
    };
  };

  return critb;
};

real PJISystem::BucklingSearchTibere()
{
  real rtoz=1.0; // o
  b2r=-0.0001;   // o

  real keff=1.0; 
  real keffold;
  real b2old=0.0;
  real b2new=0.0;
  b2=0.;
  real b2next,db;

  PutPij();
  PutPijr();

  for(int it=0;it<20;it++){

    cout<<"*************************************************\n";
    cout<<"* Iteration : "<<it<<"\n";
    cout<<"* (Square buckling) r-Dir(perp) : "<<b2r<<"\n";
    cout<<"*                   z-Dir(para) : "<<b2r*rtoz<<"\n";

    b2z=b2r*(rtoz*rtoz);
    keffold=keff;
    keff=CalIgenTibere();
    if(keff>0.9999&&keff<1.0001)break;

    if(it==0){
      b2old=b2r+b2z;
      if(b2old>0.0){
        b2new=b2old*keff;
      }
      else{
        b2new=b2old/keff;
      };
    }
    else{
      b2next=(b2new*(1-keffold)-b2old*(1-keff))/(keff-keffold);
      b2old=b2new;
      if(b2next<0.0){
	db=b2next-b2new;
	b2next=b2new+db*0.3;
      };
      b2new=b2next;
    };
    b2r=b2new/(1+rtoz*rtoz);

  };

  return b2new;
};


real PJISystem::BucklingSearchPij()
{
  bool modsrc=true; // 
  if(LeakTreat=="PseudoAbsorption")modsrc=false; // o

  real keff=1.0;
  real keffold;
  real b2old=+0.0001;
  real b2new=b2old;
  real b2next,db;

  // Pre-calculation
  PutPij(SigmaCol);
  CalIgenPij();

  for(int it=0;it<20;it++){

    cout<<"*************************************************\n";
    cout<<"* Iteration : "<<it<<"\n";
    cout<<"* (Square buckling) "<<b2new<<"\n";

    b2=b2new;

    if(!modsrc){
      PutPij(SigmaCol,b2);  
    };

    keffold=keff;

    keff=CalIgenPij();

    if(keff>0.9999&&keff<1.0001)break;

    if(it==0){
      b2old=b2;
      if(b2old>0.0){
        b2new=b2old*keff;
      }
      else{
        b2new=b2old/keff;
      };
    }
    else{
      b2next=(b2new*(1-keffold)-b2old*(1-keff))/(keff-keffold);
      b2old=b2new;
      if(b2next<0.0){
	db=b2next-b2new;
	b2next=b2new+db*0.3;
      };
      b2new=b2next;
    };
  };

  return b2new;
};


real PJISystem::BucklingSearchMSS()
{
  real keff=1.0;
  real keffold;
  real b2old=+0.0001;
  real b2new=b2old;
  real b2next,db;

  // Pre-calculation
  PutPij(SigmaCol);
  CalIgenPij();

  for(int it=0;it<20;it++){

    cout<<"*************************************************\n";
    cout<<"* Iteration : "<<it<<"\n";
    cout<<"* (Square buckling) "<<b2new<<"\n";

    b2=b2new;
    keffold=keff;

    for(int g=0;g<grp;g++){
      real cor=-b2*0.33333333/CalHomoSigt(g,SigmaCol);
      for(int i=0;i<nmed;i++){
	med[i].GetMacxs().GetData2d(sigs).add_data(g,g,cor);
      };
    };

    keff=CalIgenPij();

    for(int g=0;g<grp;g++){
      real cor=b2*0.33333333/CalHomoSigt(g,SigmaCol);
      for(int i=0;i<nmed;i++){
	med[i].GetMacxs().GetData2d(sigs).add_data(g,g,cor);
      };
    };


    if(keff>0.9999&&keff<1.0001)break;

    if(it==0){
      b2old=b2;
      if(b2old>0.0){
        b2new=b2old*keff;
      }
      else{
        b2new=b2old/keff;
      };
    }
    else{
      b2next=(b2new*(1-keffold)-b2old*(1-keff))/(keff-keffold);
      b2old=b2new;
      if(b2next<0.0){
	db=b2next-b2new;
	b2next=b2new+db*0.3;
      };
      b2new=b2next;
    };
  };

  return b2new;
};

real PJISystem::CalIgenTibere()
{
  name="PJI (Tibere)";
  PutPL(1);
  real keff=CalIgen();
  name="PJI";
  return keff;
};

real PJISystem::CalIgenPij(bool chi_mat)
{
  PutPL(0);

  if(chi_mat){
    return CalIgenWithFissionSpectrumMatrix();
  }else{
    return CalIgen();
  }
};

void PJISystem::PutApproximatedCurrent(real b2)
{
  PutPijr(SigmaCol,b2);
  // Calculation for current
  for(int g=0;g<grp;g++){
    real hsigtr=0.;
    if(b2!=0.)hsigtr=b2*0.33333333/CalHomoSigt(g,SigmaCol);
    for(int r=0;r<TotM;r++){
      real p1=0.;
      real p2=0.;
      for(int i=0;i<TotM;i++){
	p1+=mesh[i].GetFlux().get_dat(g)*pijr[g].get_dat(r*TotM+i)
          /(hsigtr+mesh[i].GetMed()->GetData1D(g,SigmaCol));
	p2+=mesh[i].GetFlux().get_dat(g)*pijz[g].get_dat(r*TotM+i)
          /(hsigtr+mesh[i].GetMed()->GetData1D(g,SigmaCol));
      };
      mesh[r].GetFlux(1).put_data(g,p1); // perp
      mesh[r].GetFlux(2).put_data(g,p2); // para
    };
  };
};

void PJISystem::PutFluxAsCurrent()
{
  for(int r=0;r<TotM;r++){
    mesh[r].GetFlux(1).copy(mesh[r].GetFlux(0));
    mesh[r].GetFlux(2).copy(mesh[r].GetFlux(0));
  };
};

Medium PJISystem::HomogenizeAll(bool benoist_d, bool micxs)
{
  Medium ret(grp);
  int nr[TotM];
  for(int i=0;i<TotM;i++){
    nr[i]=i;
  };
  ret=Homogenize(TotM,nr,benoist_d,micxs);
  return ret;
};

Medium PJISystem::Homogenize(int hreg, int *nreg,bool benoist_d,bool micxs)
{
  // -- neutron current checking --
  bool nonzero_cur=false;
  for(int i=0;i<TotM;i++){
    for(int g=0;g<grp;g++){
      if(mesh[i].GetFlux(1).get_dat(g)!=0.)nonzero_cur=true;
    };
  };
  if(!nonzero_cur){
    cout<<"\n# Error in PJISystem::Homogenize\n";
    cout<<"# Neutron current is NOT defined.\n";
    cout<<"# It is recommended to use PJISystem::PutFluxAsCurrent().\n";
    exit(0);
  };

  // -- Anisotropic Pij checking for Benoist's D --
  if(benoist_d&&!put_pijr){
    cout<<"\n# Error in PJISystem::Homogenize.\n";
    cout<<"# Anisotropic Pij is not calculated.\n";
    cout<<"# It is necessary for Benoist's D calculation.\n";
    exit(0);
  };

  Medium ret(grp);
  int pl=med[0].GetPL();
  int pltot=med[0].GetPLT();
  ret.PutPL(pl,pltot);
  ret.GetEnband().copy(med[0].GetEnband());

  GroupData1D totf(grp),totc(grp); // total flux & current
  totf.set_zero();
  totc.set_zero();

  vector<GroupData1D> cur;
  cur.resize(TotM);
  for(int i=0;i<TotM;i++){
    cur[i].put_imax(grp);
    //cur[i]=mesh[i].GetFlux();
    if(!slab)cur[i]=(mesh[i].GetFlux(1)*2.+mesh[i].GetFlux(2))*0.3333333;
    if(slab) cur[i]=(mesh[i].GetFlux(1)+mesh[i].GetFlux(2)*2.)*0.3333333;
  };

  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    totf=totf+mesh[r].GetFlux()*mesh[r].GetVolume();
    totc=totc+cur[r]*mesh[r].GetVolume();
  };

  ret.GetFlux().copy(totf);
  ret.GetFlux(1).copy(totc);

  // 1D flux-weighted cross section
  enum xstype ss[]={siga,nusigf,sigt,sign2n,d,dr,dz};
  for(int i=0;i<7;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetData1D(ss[i]).mult(mesh[r].GetFlux())
           *mesh[r].GetVolume();
    };
    rr=rr/totf;
    ret.GetData1D(ss[i]).copy(rr);
  };

  // 1D current-weighted cross section
  for(int i=0;i<1;i++){
    GroupData1D rr(grp);
    enum xstype ss;
    if(i==0)ss=sigtr;
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetData1D(ss).mult(cur[r])
           *mesh[r].GetVolume();
    };
    rr=rr/totc;
    ret.GetData1D(ss).copy(rr);
  };

  // total cross section
  for(int i=1;i<=pltot;i++){
    GroupData1D rr(grp);
    rr.set_zero();
    for(int j=0;j<hreg;j++){
      int r=nreg[j];
      rr=rr+mesh[r].GetMed()->GetSigt(i).mult(cur[r])
           *mesh[r].GetVolume();
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
          real tmp=mesh[r].GetMed()->GetDataSigs(l,i,j)*mesh[r].GetVolume();
	  if(l==0)tmp*=mesh[r].GetFlux().get_dat(i);
	  if(l!=0)tmp*=cur[r].get_dat(i);
	  bs+=tmp;
        };
        if(l==0)bs/=totf.get_dat(i);
        if(l!=0)bs/=totc.get_dat(i);
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
    fiss[i]=mesh[r].GetMed()->GetData1D(nusigf)*mesh[r].GetFlux()*mesh[r].GetVolume();
    totfis+=fiss[i];
  };
  if(totfis>0.){
    totfis=1./totfis;
    real kais;
    for(int i=0;i<grp;i++){
      kais=0.0;
      for(int j=0;j<hreg;j++){
        int r=nreg[j];
        kais+=mesh[r].GetMed()->GetData1D(chi).get_dat(i)*fiss[j];
      };
      kais*=totfis;
      ret.GetData1D(chi).put_data(i,kais);
    };
  };
  delete []fiss;


  // Fission spectrum (matrix)
  GroupData2D inp_vec(grp,grp);
  inp_vec.set_zero();
  for(int i=0;i<grp;i++){
    for(int j=0;j<grp;j++){
      for(int k=0;k<hreg;k++){
	int r=nreg[k];
        real nsf=mesh[r].GetMed()->GetData1D(nusigf).get_dat(i);
	if(nsf!=0.){
	  inp_vec.add_data
          (i,j,nsf*mesh[r].GetVolume()*mesh[r].GetFlux().get_dat(i)*mesh[r].GetMed()->GetMacxs().GetData2d(chi).get_dat(i,j));
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


  // ---------------------------------------------------------------------------------------------
  // [Fission spectrum vector is reconstructed from matrix]

  // ... matrix check
  bool non_zero=false;
  for(int g=0;g<grp;g++){
    for(int g2=0;g2<grp;g2++){
      if(inp_vec.get_dat(g,g2)>0.)non_zero=true;
    };
  };

  // ... reconstruction
  if(non_zero){
    for(int g=0;g<grp;g++){
      real sum1=0.;
      real sum2=0.;      
      for(int g2=0;g2<grp;g2++){
	real tmp=ret.GetData1D(nusigf).get_dat(g2)*ret.GetFlux().get_dat(g2);
	sum1+=inp_vec.get_dat(g2,g)*tmp;
	sum2+=tmp;
      };
      ret.GetData1D(chi).put_data(g,sum1/sum2);
    };
  };
  // ---------------------------------------------------------------------------------------------  
  

  if(benoist_d){
    // Benoist's classical definition
    for(int g=0;g<grp;g++){
      real denom=0.;
      real nume_r=0.;
      real nume_z=0.;
      for(int i=0;i<hreg;i++){
	int r=nreg[i];
	real vol=mesh[r].GetVolume();
	real flx=mesh[r].GetFlux().get_dat(g);
	denom+=vol*flx;
	for(int j=0;j<hreg;j++){
	  int rr=nreg[j];
  	  real str=mesh[rr].GetMed()->GetData1D(g,SigmaCol);
	  str+=b2*0.33333333/CalHomoSigt(g,SigmaCol);
  	  nume_r+=vol*flx/str*pijr[g].get_dat(r,rr);
  	  nume_z+=vol*flx/str*pijz[g].get_dat(r,rr);
	};
      };
      real vdr=0.33333333*nume_r/denom;
      real vdz=0.33333333*nume_z/denom;
      real tmp;
      if(slab){
	tmp=0.333333333*(vdr+vdz+vdz);
      }else{
	tmp=0.333333333*(vdr+vdr+vdz);
      };
      ret.GetData1D(d).put_data(g,tmp);
      ret.GetData1D(dr).put_data(g,vdr);
      ret.GetData1D(dz).put_data(g,vdz);
    };
  }else{
    // transport cross section homogenization with directional current
    for(int g=0;g<grp;g++){
      real tmp =0.0;
      real tmpr=0.0;
      real tmpz=0.0;
      real tmp2 =0.0;
      real tmp2r=0.0;
      real tmp2z=0.0;
      for(int i=0;i<hreg;i++){
        int r=nreg[i];
        real vol=mesh[r].GetVolume();
        real sigtrv=mesh[r].GetMed()->GetData1D(g,sigtr)*vol;
        tmp+=cur[r].get_dat(g)*sigtrv;
        tmpr+=mesh[r].GetFlux(1).get_dat(g)*sigtrv;
        tmpz+=mesh[r].GetFlux(2).get_dat(g)*sigtrv;
        tmp2+=cur[r].get_dat(g)*vol;
        tmp2r+=mesh[r].GetFlux(1).get_dat(g)*vol;
        tmp2z+=mesh[r].GetFlux(2).get_dat(g)*vol;
      };
      tmp=0.33333333*tmp2/tmp;
      tmpr=0.33333333*tmp2r/tmpr;
      tmpz=0.33333333*tmp2z/tmpz;
      ret.GetData1D(d).put_data(g,tmp);
      ret.GetData1D(dr).put_data(g,tmpr);
      ret.GetData1D(dz).put_data(g,tmpz);
    };
  };


  // +++ For microscopic cross section
  // (homogenized number density calculation)
  // (from SelfShieldingCalculator::CalMixture)
  //
  // ihnuc : number of nuclides in homogenized medium
  // ihnid[ihnuc] : nuclide id of homogenized medium
  // denh[ihnuc]  : number densities of homogenized medium
  // temp[ihnuc]  : volume-number density averated temperature

  if(micxs){
    
  vector<int> ihnid;
  int ihnuc=0;
  for(int i=0;i<hreg;i++){
    int r=nreg[i];
    for(int j=0;j<mesh[r].GetMed()->GetNucnum();j++){
      if(ihnuc==0){
	ihnuc++;
	ihnid.push_back(mesh[r].GetMed()->GetNuclideID(j));
      }else{
	int itmp=0;
	for(int k=0;k<ihnuc;k++){
	  if(ihnid[k]==mesh[r].GetMed()->GetNuclideID(j))itmp=1;
	};
	if(itmp==0){
	  ihnuc++;
	  ihnid.push_back(mesh[r].GetMed()->GetNuclideID(j));
	};
      };
    };
  };

  if(ihnuc!=0){

    vector<real> denh(ihnuc,0.);
    vector<real> temp(ihnuc,0.);
    vector<int> group_nuc(ihnuc);
    for(int i=0;i<hreg;i++){
      int r=nreg[i];
      real vol=mesh[r].GetVolume();
      for(int j=0;j<mesh[r].GetMed()->GetNucnum();j++){
        int id=mesh[r].GetMed()->GetNuclideID(j);
        real den=mesh[r].GetMed()->GetNuclideInTurn(j).GetDensity();
	if(den==0.)den=1e-20;
        // If homogenized number density is zero, spatially-averaged cross section
        // cannot be calculated.  To avoid this, extremely-small number density
	// is assumed.
        real t=mesh[r].GetMed()->GetNuclideInTurn(j).GetTemperature();
        int grp=mesh[r].GetMed()->GetNuclideInTurn(j).GetGrp();
        for(int k=0;k<ihnuc;k++){
	  if(id==ihnid[k]){
	    denh[k]+=den*vol;
            temp[k]+=t*den*vol;
	    group_nuc[k]=grp;
	  };
        };
      };
    };



    real total_vol=0.;
    for(int i=0;i<hreg;i++){
      int r=nreg[i];
      total_vol+=mesh[r].GetVolume();
    };
    total_vol=1./total_vol;
    for(int i=0;i<ihnuc;i++){
      if(denh[i]>0.)temp[i]/=denh[i];
      denh[i]*=total_vol;
    };

    ret.PutNucnum(ihnuc);
    for(int nn=0;nn<ihnuc;nn++){
      int id=ihnid[nn];
      Nuclide inpnuc;

      inpnuc.PutMatnum(id);
      inpnuc.PutDensity(denh[nn]);
      inpnuc.PutTemperature(temp[nn]);

      //cout.setf(ios::scientific);
      //cout.precision(5);
      //cout<<nn<<" "<<id<<" "<<denh[nn]<<" "<<temp[nn]<<"\n";

      if(group_nuc[nn]>0&&denh[nn]>0.){

        bool siginel_p1=true;
	for(int i=0;i<hreg;i++){
          int r=nreg[i];  
          if(mesh[r].GetMed()->ExistNuclide(id)){
       	    if(mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetDim2d(siginel)==1)siginel_p1=false;
	  };
	};

        inpnuc.PutGrp(group_nuc[nn]);
        for(int g=0;g<grp;g++){

          real totflx=0.;
          real totcur=0.;
          real cap=0.;
          real fis=0.;
          real nufis=0.;
          real n2n=0.;
          real tot0=0.;
          real tot1=0.;
          vector<real> el0(grp);
          vector<real> el1(grp);
          vector<real> ie0(grp);
          vector<real> ie1(grp);
          vector<real> n2n0(grp);
          vector<real> n2n1(grp);
          for(int i=0;i<grp;i++){
            el0[i]=0.; el1[i]=0.;
            ie0[i]=0.; ie1[i]=0.;
            n2n0[i]=0.; n2n1[i]=0.;
          };

          for(int i=0;i<hreg;i++){
            int r=nreg[i];  
     	    real volflx=mesh[r].GetFlux(0).get_dat(g)*mesh[r].GetVolume();
	    real volcur=fabs(mesh[r].GetFlux(1).get_dat(g))*mesh[r].GetVolume();
	    totflx+=volflx;
            totcur+=volcur;
            if(mesh[r].GetMed()->ExistNuclide(id)){
              real den=mesh[r].GetMed()->GetNuclide(id).GetDensity();
	      cap+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sigc).get_dat(g)*volflx;
	      fis+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sigf).get_dat(g)*volflx;
	      nufis+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(nu).get_dat(g)*
	                 mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sigf).get_dat(g)*volflx;
	      n2n+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sign2n).get_dat(g)*volflx;
	      tot0+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sigt).get_dat(g)*volflx;
	      tot1+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(sigt,1).get_dat(g)*volcur;
              for(int g2=0;g2<grp;g2++){
	        el0[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(sigel).get_dat(g,g2)*volflx;
	        el1[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(sigel,1).get_dat(g,g2)*volflx;
	        ie0[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(siginel).get_dat(g,g2)*volflx;
	        if(siginel_p1)ie1[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(siginel,1).get_dat(g,g2)*volflx;
	        n2n0[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(sign2n).get_dat(g,g2)*volflx;
	    //n2n1[g2]+=den*mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData2d(sign2n,1).get_dat(g,g2)*volflx;
	      };
	    };
          };

          inpnuc.GetMicxs().GetData1d(sigc).put_data(g,cap/totflx/denh[nn]);
          inpnuc.GetMicxs().GetData1d(sign2n).put_data(g,n2n/totflx/denh[nn]);
          inpnuc.GetMicxs().GetData1d(sigt).put_data(g,tot0/totflx/denh[nn]);
          inpnuc.GetMicxs().GetData1d(sigt,1).put_data(g,tot1/totflx/denh[nn]);
          real tmp=fis/totflx/denh[nn];
          inpnuc.GetMicxs().GetData1d(sigf).put_data(g,tmp);
          if(tmp!=0.){
            inpnuc.GetMicxs().GetData1d(nu).put_data(g,nufis/fis);
          }else{
            inpnuc.GetMicxs().GetData1d(nu).put_data(g,0.);
          };
          for(int g2=0;g2<grp;g2++){
            inpnuc.GetMicxs().GetData2d(sigel).put_data(g,g2,el0[g2]/totflx/denh[nn]);
            if(siginel_p1)inpnuc.GetMicxs().GetData2d(sigel,1).put_data(g,g2,el1[g2]/totcur/denh[nn]);
            inpnuc.GetMicxs().GetData2d(siginel).put_data(g,g2,ie0[g2]/totflx/denh[nn]);
            inpnuc.GetMicxs().GetData2d(siginel,1).put_data(g,g2,ie1[g2]/totcur/denh[nn]);
            inpnuc.GetMicxs().GetData2d(sign2n).put_data(g,g2,n2n0[g2]/totflx/denh[nn]);
            //inpnuc.GetMicxs().GetData2d(sign2n,1).put_data(g,g2,n2n1[g2]/totcur/denh[nn]);
          };
        };

        // fission spectrum (to be modified?)
        if(inpnuc.GetMicxs().GetData1d(nu).get_dat(0)>0.){
          for(int i=0;i<hreg;i++){
            int r=nreg[i];  
            if(mesh[r].GetMed()->ExistNuclide(id)){
	      inpnuc.GetMicxs().GetData1d(chi).copy(mesh[r].GetMed()->GetNuclide(id).GetMicxs().GetData1d(chi));
  	      i=hreg;
	    };
	  };
        };

      };

      ret.PutNuclide(nn,inpnuc);

    };


  };

  };

  return ret;
}

void PJISystem::PutRegMed(vector<int> inp)
{
  int sz=inp.size();
  int *in=new int[sz];
  for(int i=0;i<sz;i++){in[i]=inp[i];};
  PutRegMed(in);
  delete []in;
};

void PJISystem::PutRegMed(int *inp)
{
  for(int i=0;i<TotM;i++){
    if(inp[i]<0||inp[i]>=nmed){
      cout<<"Medium ID is not appropriated in PJISystem!\n";
      cout<<"Medium ID : "<<inp[i]<<"\n";
      exit(0);
    };
    mesh[i].PutMedium(&med[inp[i]]);
    medium_id_par_region[i]=inp[i];
  };


  int *fmx=new int[TotM];
  real *xl=new real[TotM];
  int fmy[]={1};
  real yl[]={1.};
  int *mat=new int[TotM];
  for(int i=0;i<TotM;i++){
    fmx[i]=1;
    xl[i]=mesh[i].GetVolume();
    mat[i]=medium_id_par_region[i];
  };
  mi.PutMeshInfo(TotM,1,fmx,fmy,xl,yl,mat,"width");
  
  delete [] fmx;
  delete [] xl;
  delete [] mat;

};


void PJISystem::TransformForAdjointFlux()
{
  vector<real> src(TotM);
  for(int i=0;i<TotM;i++){
    src[i]=mesh[i].GetFissionSrc()/mesh[i].GetVolume();
  };
  for(int g=0;g<grp;g++){
    vector<real> newflx(TotM);
    for(int i=0;i<TotM;i++){
      real tmp=0.;
      for(int g2=0;g2<grp;g2++){
	tmp+=mesh[i].GetMed()->GetMacxs().GetData2d(sigs).get_dat(g,g2)
	  *mesh[i].GetFlux().get_dat(g2);
      };
      newflx[i]=tmp+mesh[i].GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g)
	*src[i];
    };
    for(int i=0;i<TotM;i++){ 
      mesh[i].GetFlux().put_data(g,newflx[i]);
    };
  };
};

// **********************************************
// * For perturbation calculation               *
// **********************************************

real PJISystem::CalReactivity(PJISystem *sec,real kunp,real kp,bool pr,bool ipcal)
{
  vector<real> tmp;
  tmp.resize(grp,0.);
  return CalReactivity(sec,kunp,kp,tmp,pr,ipcal);
};


real PJISystem::CalReactivity(PJISystem *sec,real kunp,real kp,vector<real> &rho,bool pr,bool ipcal)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  real *yld=new real[grp];
  real *abs=new real[grp];
  real *n2n=new real[grp];    
  vector< vector<real> > sct(plnum);
  vector< vector<real> > leak(plnum);
  for(int i=0;i<plnum;i++){
    sct[i].resize(grp,0.);
    leak[i].resize(grp,0.);
  };

  real ip=1.;
  if(ipcal){
    ip=CalPerturbDenominator(sec);
  //cout<<"# Perturb denominator is : "<<ip<<"\n";
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real *scttmp=new real[plnum];
  real *leaktmp=new real[plnum];
  for(int i=0;i<grp;i++){
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    abs[i]=CalPerturbAbsorptionTerm(sec,flag,i);
    n2n[i]=CalPerturbN2NTerm(sec,flag,i);            
    sct[0][i]=CalPerturbScatteringTerm(sec,i,flag);
    //CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
    for(int j=0;j<plnum;j++){
      //sct[j][i]=scttmp[j];
      //leak[j][i]=leaktmp[j];
    };
  };
  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;

  real yldsum=0.;
  real abssum=0.;
  real n2nsum=0.;
  real tsctsum=0.;
  real tleaksum=0.;
  real *sctsum=new real[plnum];
  real *leaksum=new real[plnum];
  for(int i=0;i<plnum;i++){
    sctsum[i]=0.; leaksum[i]=0.;
  };

  if(pr){
    cout<<"# (Unit of Reactivity : (Reactivity)/(Unit_lethargy:0.25)\n";
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scattering  ";
    cout<<"Total\n";
  };

  real inv_ip=1./ip;
  for(int i=0;i<grp;i++){
    //          for(int i=27;i<91;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    n2n[i]*=inv_ip;    
    real sl=0.;
    real le=0.;
    for(int j=0;j<plnum;j++){
      sct[j][i]*=inv_ip; sctsum[j]+=sct[j][i]; sl+=sct[j][i];
      leak[j][i]*=inv_ip; leaksum[j]+=leak[j][i]; le+=leak[j][i];
    };
    yldsum+=yld[i];
    abssum+=abs[i];
    n2nsum+=n2n[i];    
    tsctsum+=sl;
    tleaksum+=le;
    real grpsum=yld[i]+abs[i]+sl+le;
    rho[i]=grpsum;
    if(pr){
      cout.width(3);
      cout<<i<<"   ";
      cout.setf(ios::scientific);
      cout.precision(4);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<grpsum/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# Yield       :"<<yldsum<<"\n";
    cout<<"# Absorption  :"<<abssum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";        
    int id=0;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th scattering  :"<<sctsum[id]<<"\n";
        id++;
      };
    };
    cout<<"# Scattering   : "<<tsctsum<<"\n";
    id=0;
    for(int l=1;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th leakage     :"<<leaksum[id]<<"\n";
        id++;
      };
    };
    cout<<"#     Higher leakage :"<<leaksum[plnum-1]<<"\n";
    cout<<"# Leakage      :"<<tleaksum<<"\n";
  };

  cout.setf(ios::scientific);
  cout.precision(4);

  real tot=yldsum+abssum+tsctsum+tleaksum;
  if(pr)cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
  if(pr)cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
  cout.unsetf(ios::scientific);

  delete [] sctsum;
  delete [] leaksum;
  delete [] yld;
  delete [] abs;
  delete [] n2n;
  
  return tot;
};

real PJISystem::CalReactivity(PJISystem *sec,real kunp,real kp,bool pr,int meshid,bool ipcal)
{
  vector<real> rho;
  rho.resize(grp,0.);

  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };


  real *yld=new real[grp];
  real *abs=new real[grp];
  real *n2n=new real[grp];  
  vector< vector<real> > sct(plnum);
  vector< vector<real> > leak(plnum);
  for(int i=0;i<plnum;i++){
    sct[i].resize(grp,0.);
    leak[i].resize(grp,0.);
  };

  real ip=1.;
  if(ipcal){
    ip=CalPerturbDenominator(sec);
  //cout<<"# Perturb denominator is : "<<ip<<"\n";
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=false;
  };
  flag[meshid]=true;

  real *scttmp=new real[plnum];
  real *leaktmp=new real[plnum];
  for(int i=0;i<grp;i++){
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    abs[i]=CalPerturbAbsorptionTerm(sec,flag,i);
    n2n[i]=CalPerturbN2NTerm(sec,flag,i);    
    sct[0][i]=CalPerturbScatteringTerm(sec,i,flag);
    //CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
    for(int j=0;j<plnum;j++){
      //sct[j][i]=scttmp[j];
      //leak[j][i]=leaktmp[j];
    };
  };
  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;

  real yldsum=0.;
  real abssum=0.;
  real n2nsum=0.;
  real tsctsum=0.;
  real tleaksum=0.;
  real *sctsum=new real[plnum];
  real *leaksum=new real[plnum];
  for(int i=0;i<plnum;i++){
    sctsum[i]=0.; leaksum[i]=0.;
  };

  if(pr){
    cout<<"# (Unit of Reactivity : (Reactivity)/(Unit_lethargy:0.25)\n";
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scattering  ";
    cout<<"Total\n";
  };

  real inv_ip=1./ip;
  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    n2n[i]*=inv_ip;
    real sl=0.;
    real le=0.;
    for(int j=0;j<plnum;j++){
      sct[j][i]*=inv_ip; sctsum[j]+=sct[j][i]; sl+=sct[j][i];
      leak[j][i]*=inv_ip; leaksum[j]+=leak[j][i]; le+=leak[j][i];
    };
    yldsum+=yld[i];
    abssum+=abs[i];
    n2nsum+=n2n[i];
    tsctsum+=sl;
    tleaksum+=le;
    real grpsum=yld[i]+abs[i]+sl+le;
    rho[i]=grpsum;
    if(pr){
      cout.width(3);
      cout<<i<<"   ";
      cout.setf(ios::scientific);
      cout.precision(4);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<grpsum/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# Yield       :"<<yldsum<<"\n";
    cout<<"# Absorption  :"<<abssum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";    
    int id=0;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th scattering  :"<<sctsum[id]<<"\n";
        id++;
      };
    };
    cout<<"# Scattering   : "<<tsctsum<<"\n";
    id=0;
    for(int l=1;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th leakage     :"<<leaksum[id]<<"\n";
        id++;
      };
    };
    cout<<"#     Higher leakage :"<<leaksum[plnum-1]<<"\n";
    cout<<"# Leakage      :"<<tleaksum<<"\n";
  };

  cout.setf(ios::scientific);
  cout.precision(4);

  real tot=yldsum+abssum+tsctsum+tleaksum;
  if(pr)cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
  if(pr)cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
  cout.unsetf(ios::scientific);

  delete [] sctsum;
  delete [] leaksum;
  delete [] yld;
  delete [] abs;
  delete [] n2n;  

  return tot;
};

GData PJISystem::CalReactivityGData(PJISystem *sec,real kunp,real kp,bool pr,bool ipcal)
{
  GData ret(3);
  ret.PutTagX("eV");
  ret.PutTagY(0,"yield");
  ret.PutTagY(1,"absorption");
  ret.PutTagY(2,"scattering");    
  for(int g=0;g<grp;g++){
    ret.push_back_x(mesh[0].GetMed()->GetEnband().get_dat(g));
  };

  
  vector<real> rho;
  rho.resize(grp,0.);

  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in PJISystem::CalReactivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };


  real *yld=new real[grp];
  real *abs=new real[grp];
  real *n2n=new real[grp];  
  vector< vector<real> > sct(plnum);
  vector< vector<real> > leak(plnum);
  for(int i=0;i<plnum;i++){
    sct[i].resize(grp,0.);
    leak[i].resize(grp,0.);
  };

  real ip=1.;
  if(ipcal){
    ip=CalPerturbDenominator(sec);
  //cout<<"# Perturb denominator is : "<<ip<<"\n";
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real *scttmp=new real[plnum];
  real *leaktmp=new real[plnum];
  for(int i=0;i<grp;i++){
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    abs[i]=CalPerturbAbsorptionTerm(sec,flag,i);
    n2n[i]=CalPerturbN2NTerm(sec,flag,i);    
    sct[0][i]=CalPerturbScatteringTerm(sec,i,flag);
    //CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
    for(int j=0;j<plnum;j++){
      //sct[j][i]=scttmp[j];
      //leak[j][i]=leaktmp[j];
    };
  };
  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;

  real yldsum=0.;
  real abssum=0.;
  real n2nsum=0.;
  real tsctsum=0.;
  real tleaksum=0.;
  real *sctsum=new real[plnum];
  real *leaksum=new real[plnum];
  for(int i=0;i<plnum;i++){
    sctsum[i]=0.; leaksum[i]=0.;
  };

  if(pr){
    cout<<"# (Unit of Reactivity : (Reactivity)/(Unit_lethargy:0.25)\n";
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scattering  ";
    cout<<"Total\n";
  };

  real inv_ip=1./ip;
  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    n2n[i]*=inv_ip;
    real sl=0.;
    real le=0.;
    for(int j=0;j<plnum;j++){
      sct[j][i]*=inv_ip; sctsum[j]+=sct[j][i]; sl+=sct[j][i];
      leak[j][i]*=inv_ip; leaksum[j]+=leak[j][i]; le+=leak[j][i];
    };
    yldsum+=yld[i];
    abssum+=abs[i];
    n2nsum+=n2n[i];
    tsctsum+=sl;
    tleaksum+=le;
    real grpsum=yld[i]+abs[i]+sl+le;
    rho[i]=grpsum;

    ret.push_back_y(0,yld[i]);
    ret.push_back_y(1,abs[i]);
    ret.push_back_y(2,sl);

    if(pr){
      cout.width(3);
      cout<<i<<"   ";
      cout.setf(ios::scientific);
      cout.precision(4);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<grpsum/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# Yield       :"<<yldsum<<"\n";
    cout<<"# Absorption  :"<<abssum<<"\n";
    cout<<"# N2N        : "<<n2nsum<<"\n";    
    int id=0;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th scattering  :"<<sctsum[id]<<"\n";
        id++;
      };
    };
    cout<<"# Scattering   : "<<tsctsum<<"\n";
    id=0;
    for(int l=1;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th leakage     :"<<leaksum[id]<<"\n";
        id++;
      };
    };
    cout<<"#     Higher leakage :"<<leaksum[plnum-1]<<"\n";
    cout<<"# Leakage      :"<<tleaksum<<"\n";
  };

  cout.setf(ios::scientific);
  cout.precision(4);

  real tot=yldsum+abssum+tsctsum+tleaksum;
  if(pr)cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
  if(pr)cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
  cout.unsetf(ios::scientific);

  delete [] sctsum;
  delete [] leaksum;
  delete [] yld;
  delete [] abs;
  delete [] n2n;  

  return ret;
};

real PJISystem::CalPerturbScatteringTerm(PJISystem *sec,int g,bool *flag)
{
  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real f2v=sec->GetMesh(m).GetFlux(0).get_dat(g)*mesh[m].GetVolume();
      // (forward-flux)*(volume)
      real tmp=0.;
      //real tmp2=0.;
      for(int j=0;j<grp;j++){ // w up scattering
	//for(int j=g;j<grp;j++){ // w/o up scattering
        tmp+=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g,j)
	             -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,j))*
	  mesh[m].GetFlux(0).get_dat(j);
	//tmp2+= sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g,j)
  	//              -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,j);

      };
      //tmp+=-tmp2*mesh[m].GetFlux(0).get_dat(g);
      tmp+=
	-( (sec->GetMesh(m).GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)-
            sec->GetMesh(m).GetMed()->GetMacxs().GetData1d(siga).get_dat(g) )
        -  (mesh[m].GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)-
            mesh[m].GetMed()->GetMacxs().GetData1d(siga).get_dat(g))
          )*mesh[m].GetFlux(0).get_dat(g);
      ret+=tmp*f2v;
    };
  };
  return ret;
};

void PJISystem::CalSensitivity(PJISystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  CheckAdjointForward(sec);
  CheckSameMesh(sec); // o

  real *nsforg=new real[nmed];
  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  real *sigsorg=new real[nmed*grp];
  real *fiss_frac=new real[nmed];

  //real delta=1.;

  cout.setf(ios::scientific);
  cout.precision(5);

  real ip=CalPerturbDenominator(sec);

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  cout<<grp<<"\n";
  for(int nc=0;nc<nucnum;nc++){
    int nid=nucid[nc];
    int nidendf=nid;

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
            real den=med[j].GetNuclide(nid).GetDensity();
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
            absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
            real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
	  };
        };
        real re=0.;
        re+=CalPerturbYieldTerm(sec,flag,i,k2);
        re+=CalPerturbAbsorptionTerm(sec,flag,i);
        re/=ip;
	cout<<"  "<<re<<"\n";
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
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
            real den=med[j].GetNuclide(nid).GetDensity();
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micsigf*micnu);
	  };
        };
        real re=0.;
        re+=CalPerturbYieldTerm(sec,flag,i,k2);
        re/=ip;
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
  	  real den=med[i].GetNuclide(nid).GetDensity();
 	  real total_fiss=sec->GetIntegratedFlux(i)*sec->GetMed(i).GetMacxs().GetNusigf();
          real part_fiss=den*(sec->GetIntegratedFlux(i)*sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicNusigf());
	  fiss_frac[i]=part_fiss/total_fiss;
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
	re/=ip;
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
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
        };
      };
      real re=CalPerturbAbsorptionTerm(sec,flag,i);
      re/=ip;
      cout<<"  "<<re<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	};
      };
    };

    // Scattering P0
    for(int ii=0;ii<3;ii++){
      cout<<2<<"\n";
      cout<<nidendf<<"\n";
      //
      if(ii==0){cout<<2<<"\n";}
      else if(ii==1){cout<<4<<"\n";}
      else {cout<<16<<"\n";};
      //
      for(int i=0;i<grp;i++){
        for(int k=i;k<grp;k++){
	  if((ii!=0)||(nid<100000)||(ii==0&&k<=i+2)){
            for(int j=0;j<nmed;j++){
              if(med[j].ExistNuclide(nid)){
                real den=med[j].GetNuclide(nid).GetDensity();
                sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
                totorg[j]=sec->GetMed(j).GetMacxs().GetSigt().get_dat(i);
	        real micsigs=0.;
	        switch(ii){
	        case 0:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	          break;
	        case 1:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_dat(i,k);
		  break;
	        case 2:
		  micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_dat(i,k);
	          break;
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
	        sec->GetMed(j).GetMacxs().GetSigt().add_data(i,den*micsigs);
              };
            };
            real re=CalPerturbScatteringTerm(sec,i,flag);
            re/=ip;
            cout<<"  "<<re<<"\n";
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
	        sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	        sec->GetMed(j).GetMacxs().GetSigt().put_data(i,totorg[j]);
	      };
	    };
	  }else{
	    cout<<"  0.0\n";
	  };
	};
      };
    };
  };

  delete [] totorg;
  delete [] fiss_frac;
  delete [] flag;
  delete [] sigsorg;
  delete [] absorg;
  delete [] nsforg;
};

SensitivityData PJISystem::CalSensitivityNew(PJISystem *sec,real keff,int nucnum,int *nucid,bool ipcal)
{
  // When sensitivity of k_inf is calculated, 
  // calculation results should be converted from reactivity (dk/kk) to sensitivity (dk/k),
  // so [keff] is multiplyed to calculation results.
  // On the other hand, when sensitivity of non-k_inf quantity is calculated,
  // conversions of calculation results are NOT required.
  // In such cases, [keff] should be unity;
  // but [keff] is also used in fission yield term calculations.
  // Thus, the following modification is made on the algorithm.
  //  - [keff] should be correctly input.
  //  - if ipcal is false, the above conversion is NOT made.

  CheckAdjointForward(sec);
  CheckSameMesh(sec); // o

  SensitivityData sens;
  sens.PutValue(keff);
  sens.PutGroup(grp);
  sens.GetEnband().copy(med[0].GetEnband());

  GroupData1D sns1d(grp);
  GroupData2D sns2d(grp,grp);

  real *nsforg=new real[nmed];
  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  real *sigsorg=new real[nmed*grp];
  real *fiss_frac=new real[nmed];

  real delta=1.;

  bool *flag=new bool[TotM];

  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real ip=1.;
  if(ipcal){
    ip=CalPerturbDenominator(sec);
    if(ip==0.){
      cout<<"# Error in PJIsystem::CalSensitivityNew.\n";
      cout<<"# Perturbation denominator is zero.\n";
      exit(0);
    };
  };

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
        // +++ Fission
        for(int i=0;i<grp;i++){
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
              nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
              absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
              real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
              sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
	    };
          };
          real re=0.;
          re+=CalPerturbYieldTerm(sec,flag,i,keff);
          re+=CalPerturbAbsorptionTerm(sec,flag,i);
          re/=ip;
          if(ipcal)re*=keff;
          sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
  	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	      sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	    };
	  };
        };
        sens.PutSensitivity1D(nid,18,sns1d);

        // +++ Nu
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
          re/=ip;
          if(ipcal)re*=keff;
	  sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    };
          };
        };
        sens.PutSensitivity1D(nid,452,sns1d);

        // +++ Chi
        for(int i=0;i<nmed;i++){
	  if(med[i].ExistNuclide(nid)){
  	    real den=med[i].GetNuclide(nid).GetDensity();
 	    real total_fiss=sec->GetIntegratedFlux(i)*sec->GetMed(i).GetMacxs().GetNusigf();
            real part_fiss=den*(sec->GetIntegratedFlux(i)*sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicNusigf());
	    fiss_frac[i]=part_fiss/total_fiss;
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
          re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  re/=ip;
          if(ipcal)re*=keff;
	  sns1d.put_data(i,re);
  	  for(int j=0;j<nmed;j++){
	    if(med[j].ExistNuclide(nid)){
	      sec->GetMed(j).GetMacxs().GetKai().put_data(i,totorg[j]);
	    };
	  };
        };
        sens.PutSensitivity1D(nid,181,sns1d);
      }; //(end for fissile nuclide)


      // +++ Capture
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
    	  if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
	    real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
          };
        };
        real re=CalPerturbAbsorptionTerm(sec,flag,i);
        re/=ip;
        if(ipcal)re*=keff;
        sns1d.put_data(i,re);   
        for(int j=0;j<nmed;j++){
  	  if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  };
        };
      };
      sens.PutSensitivity1D(nid,102,sns1d);

      // +++ Scattering P0
      for(int ii=0;ii<3;ii++){
        int mt=2;
	if(ii==1)mt=4;
	if(ii==2)mt=16;
	sns2d.set_zero();
        for(int i=0;i<grp;i++){
          //for(int k=0;k<grp;k++){
	  for(int k=i;k<grp;k++){
  	    if((ii!=0)||(nid<100000)||(ii==0&&k<=i+2)){
              for(int j=0;j<nmed;j++){
                if(med[j].ExistNuclide(nid)){
                  real den=med[j].GetNuclide(nid).GetDensity();
                  sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
                  totorg[j]=sec->GetMed(j).GetMacxs().GetSigt().get_dat(i);
	          real micsigs=0.;
	          switch(ii){
	          case 0:
	            micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	            break;
	          case 1:
	            micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_dat(i,k);
		    break;
	          case 2:
		    micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_dat(i,k);
                    absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);		    
	            break;
	          };
		  if(ii!=0)micsigs=delta;
		  // 'absolute' sensitivity for inelastic & (n,2n) 
	          sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
		  real factor=1.;
		  if(ii==2)factor=0.5;
	          sec->GetMed(j).GetMacxs().GetSigt().add_data(i,den*micsigs*factor);
		  if(ii==2)sec->GetMed(j).GetMacxs().GetSiga().add_data(i,-den*micsigs*factor);
                };
              };
              real re=CalPerturbScatteringTerm(sec,i,flag);
              re+=CalPerturbAbsorptionTerm(sec,flag,i);	      
              re/=ip;
              if(ipcal)re*=keff;
              sns2d.put_data(i,k,re);
              for(int j=0;j<nmed;j++){
  	        if(med[j].ExistNuclide(nid)){
	          sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	          sec->GetMed(j).GetMacxs().GetSigt().put_data(i,totorg[j]);
		  if(ii==2)sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	        };
	      };
	    };
	  };
	};
        sens.PutSensitivity2D(nid,mt,sns2d);
      };

    };
  };

  delete [] totorg;
  delete [] fiss_frac;
  delete [] flag;
  delete [] sigsorg;
  delete [] absorg;
  delete [] nsforg;

  return sens;

};

SensitivityData PJISystem::CalSensitivityRRR(PJISystem &lat,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh,int nucnum,int *nucid,real keff)
{
  real nume=lat.CalMacroscopicReactionRate(nume_nuc,nume_id,nume_xs,on_mesh);
  real denom=lat.CalMacroscopicReactionRate(denom_nuc,denom_id,denom_xs,on_mesh);
  real rrr=nume/denom;

  vector<GroupData1D> gpt_flx(TotM);
  for(int i=0;i<TotM;i++){
    gpt_flx[i].put_imax(grp);
  };

  CalGPT(keff,nume_nuc,nume_id,nume_xs,on_mesh,nume);
  for(int i=0;i<TotM;i++){
    gpt_flx[i]=GetMesh(i).GetFlux();
  };

  CalGPT(keff,denom_nuc,denom_id,denom_xs,on_mesh,denom);
  for(int i=0;i<TotM;i++){
    gpt_flx[i]=gpt_flx[i]-mesh[i].GetFlux();
    mesh[i].GetFlux().copy(gpt_flx[i]);
  };

  SensitivityData sns=CalSensitivityNew(&lat,keff,nucnum,nucid,false);// `false' means that perturbation denominator calculation is neglected.
  //SensitivityData sns=CalSensitivityNew(&lat,1.,nucnum,nucid,false);// `false' means that perturbation denominator calculation is neglected.
  sns.PutValue(rrr);
  sns.PutName("pin","RRR","unknown");
  sns.WriteFile("./","sns.flx");
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // +++ Sensitivity calculation for direct term
  SensitivityData sns_dir;
  sns_dir.PutName("pin","RRR","unknown");
  sns_dir.PutValue(rrr);
  sns_dir.PutGroup(grp);
  sns_dir.GetEnband().copy(lat.GetMesh(0).GetMed()->GetEnband());
  GroupData1D snsdat(grp);
  // (numerator)
  for(int i=0;i<nume_nuc;i++){
    if(lat.GetMesh(0).GetMed()->ExistNuclide(nume_id[i])){
      real den=lat.GetMesh(0).GetMed()->GetNuclide(nume_id[i]).GetDensity();
      for(int g=0;g<grp;g++){
        real xs=lat.GetMesh(0).GetMed()->GetNuclide(nume_id[i]).GetMicxs().GetData1d(nume_xs[i]).get_dat(g);
        real tmp=0.;
        for(int m=0;m<TotM;m++){
  	  if(on_mesh[m]){
            real volflx=lat.GetMesh(m).GetFlux().get_dat(g)*lat.GetMesh(m).GetVolume();
   	    tmp+=volflx*xs*den;
	  };
        };
        snsdat.put_data(g,tmp/nume);
      };
      int mt=102;
      if(nume_xs[i]==sigf)mt=18;
      sns_dir.PutSensitivity1D(nume_id[i],mt,snsdat);
    };
  };
  // (denomenator)
  for(int i=0;i<denom_nuc;i++){
    if(lat.GetMesh(0).GetMed()->ExistNuclide(denom_id[i])){
      real den=lat.GetMesh(0).GetMed()->GetNuclide(denom_id[i]).GetDensity();
      for(int g=0;g<grp;g++){
        real xs=lat.GetMesh(0).GetMed()->GetNuclide(denom_id[i]).GetMicxs().GetData1d(denom_xs[i]).get_dat(g);
        real tmp=0.;
        for(int m=0;m<TotM;m++){
	  if(on_mesh[m]){
            real volflx=lat.GetMesh(m).GetFlux().get_dat(g)*lat.GetMesh(m).GetVolume();
  	    tmp+=volflx*xs*den;
	  };
        };
        snsdat.put_data(g,-tmp/denom);
      };
      int mt=102;
      if(denom_xs[i]==sigf)mt=18;
      sns_dir.PutSensitivity1D(denom_id[i],mt,snsdat);
    };
  };
  sns_dir.WriteFile("./","sns.dir");

  sns.AddSensitivityData(sns_dir);
  return sns;
};

GroupData1D PJISystem::GetIntegratedFlux(int medid,int mom)
{
  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  for(int i=0;i<TotM;i++){
    if(medium_id_par_region[i]==medid){
      real vol=mesh[i].GetVolume();
      for(int g=0;g<grp;g++){
	sum[g]+=fabs(mesh[i].GetFluxData(mom,g))*vol;
      };
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};

real PJISystem::GetVolumePerMedium(int medid)
{
  real ret=0.;
  for(int i=0;i<TotM;i++){
    if(medium_id_par_region[i]==medid){
      real vol=mesh[i].GetVolume();
      ret+=vol;
    };
  };
  return ret;
};

void PJISystem::CopyPij(PJISystem &sys2)
{

  real *pijinp=new real[TotM*TotM];
  for(int ig=0;ig<grp;ig++){
    pij[ig].put_yx(TotM,TotM);
    GroupData2D tmp=sys2.GetPij(ig);
    for(int i=0;i<TotM*TotM;i++){
      pijinp[i]=tmp.get_dat(i);
    };
    pij[ig].put_data(pijinp);
  };

  delete [] pijinp;
}

/*
void PJISystem::AddFlux(PJISystem &sec)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)+sec.GetMesh(i).GetFlux(l);
      mesh[i].GetFlux(l).copy(tmp);
    };
  };
};

void PJISystem::NegFlux(PJISystem &sec)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)-sec.GetMesh(i).GetFlux(l);
      mesh[i].GetFlux(l).copy(tmp);
    };
  };
};

void PJISystem::NegFlux(PJISystem &sec,real fact)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)-sec.GetMesh(i).GetFlux(l)*fact;
      mesh[i].GetFlux(l).copy(tmp);
    };
  };
};
*/

void PJISystem::ClearFlux()
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      mesh[i].GetFlux(l).set_zero();
    };
  };
};

