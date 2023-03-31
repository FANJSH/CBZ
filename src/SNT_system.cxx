#include <cstdlib>
#include "SNT_system.h"

using namespace std;

void SNTSystem::Init(int n,int g,int i)
{
  GeneralSystem::Init(n,g,i);
  psys.Init(n,g,i);
  pl = -1;
  nquad = 1;
  quad.resize(nquad);
  quadid.resize(grp,0);
  for(int i=0;i<grp;i++){quadid[i]=0;};
  transport=true;
  name="SNT";
  cmfdimp = true;  // CMFD IS implemented
  wrtflx_endgrp=grp-1;
  dsa=true;
  itmax_inner=30;
};

SNTSystem::~SNTSystem()
{
};

void SNTSystem::SetQuadratureNum(int num)
{
  nquad=num;
  quad.resize(nquad);
};

void SNTSystem::SetQuadratureID(int *gbnd, int *id)
{
  int now=0;
  for(int i=0;i<grp;i++){
    if(id[now]<0||id[now]>=nquad){
      cout<<"Error in SetQuadratureID.\n";
      cout<<" Max quadrature number  = "<<nquad<<"\n";
      cout<<" your ID number         = "<<id[now]<<"\n";
      exit(0);
    };
    quadid[i]=id[now];
    if(gbnd[now]==i)now++;
  };
};

void SNTSystem::SetQuadrature(Quadrature *qinp,int id)
{
  if(id<0||id>=nquad){
    cout<<"Error in SetQuadratureID.\n";
    cout<<"Chack your ID number.\n";
    exit(0);
  };
  quad[id]=qinp;
};

void SNTSystem::SetArray()
{
  if(pl==-1){
    cout<<"Do 'PutPL' for system class.\n";
    exit(0);
  };

  for(int i=0;i<grp;i++){
    if(!quad[quadid[i]]->Exist()){
      cout<<"Quadrature is not defined yet.\n";
      exit(0);
    };
  };

  for(int i=0;i<TotM;i++){
    mesh[i].PutPL(pl);
    mesh[i].PutGrp(grp);
    if(wrtflx[i]){
      aflux[i].resize(wrtflx_endgrp+1);
      for(int j=0;j<=wrtflx_endgrp;j++){
        aflux[i][j].put_imax(quad[quadid[j]]->GetSN());
      };
    };
  };

  if(pl>0&&transport){
    for(int i=0;i<nmed;i++){
      int plt=med[i].GetPLT();
      for(int j=1;j<pl+1;j++){
        int jj=j;
        if(j>plt)jj=plt;
        for(int k=0;k<grp;k++){
	  real tmp=med[i].GetDataSigs(j,k,k);
	  real tmp2=med[i].GetDataSigt(0,k)-med[i].GetDataSigt(jj,k);
	  real tmp3=tmp+(j*2+1.)*tmp2;
	  real tmp4=med[i].GetSigs(0).get_dat(k,k);
	  if(j==1&&fabs(tmp3*0.33333333)>tmp4){
	    cout<<"#   ... warning... mu (sp1/sp0) is out of range : "<<tmp3*0.33333333/tmp4<<"\n";
	    cout<<"#         medium "<<i<<"  group "<<k<<"\n";
	  }
          else{
            med[i].GetSigs(j).put_data(k,k,tmp3);
	  };
	};
      };
    };
  };

  if(pl==0&&transport){
    for(int i=0;i<nmed;i++){
      //med[i].CalSigtrFromDav();
      med[i].TransportApproximation();
    };
  };

  for(int i=0;i<nquad;i++){
    quad[i]->CheckOrthogonalityTo00();
  };

  // For CMFD
  curxp.resize(TotM,0.);
  curxn.resize(TotM,0.);
  curyp.resize(TotM,0.);
  curyn.resize(TotM,0.);
  curzp.resize(TotM,0.);
  curzn.resize(TotM,0.);
};

void SNTSystem::PutWriteFlux(int med1,int med2)
{
  for(int i=med1;i<=med2;i++){
    PutWriteFlux(i);
  };
};

void SNTSystem::PutWriteFlux(int med)
{
  int index=0;
  for(int k=0;k<mi.GetZF();k++){
    for(int i=0;i<mi.GetYF();i++){
      for(int j=0;j<mi.GetXF();j++){
	int id=meshid[k][i][j];
	if(id!=-1){
 	  int tmp=mi.GetFMat(index);
	  if(tmp==med||med==-2)wrtflx[id]=true;
	};
	index++;
      };
    };
  };
};

void SNTSystem::PutWriteFluxEndGrp(int i)
{
  if(i<0||i>=grp){
    cout<<"Error in PutWriteFluxEndGrp.\n";
    cout<<"Max number of defined group is "<<grp+1<<"\n";
    cout<<"You shose group "<<i+1<<"("<<i<<")\n";
    exit(0);
  };
  wrtflx_endgrp=i;
};

void SNTSystem::PutCartMeshInfo(CartMeshInfo cm, string geom)
{
  bool cart=false;
  if(geom=="Cartesian")        cart    =true;

  if(!cart){
    cout<<"Error in PutCartMeshInfo.\n";
    cout<<"...Cartesian? \n";
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
  real *di=new real[Ndim];

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
	if(tm>=nmed){
	  cout<<"Error in SNTSystem::PutCartMeshInfo.\n";
	  cout<<"You requested not-existing medium ID.\n";
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
  wrtflx.resize(TotM,false);
  aflux.resize(TotM);

  int index=0;
  int ind2=0;
  for(int z=0;z<zf;z++){
    for(int y=0;y<yf;y++){
      real tmp=0.;
      for(int x=0;x<xf;x++){
        tmp+=mi.GetFMeshL(0,x);
        di[0]=mi.GetFMeshL(0,x);
        if(Ndim>1)di[1]=mi.GetFMeshL(1,y);
        if(Ndim>2)di[2]=mi.GetFMeshL(2,z);
	int tm=mi.GetFMat(ind2);
	ind2++;
	if(tm!=-1){
          mesh[index].PutDim(Ndim,di);
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

  int ind=0;
  for(int i=0;i<Ndim;i++){
    for(int j=0;j<2;j++){
      if(mi.GetBC(ind)==2)BC[ind]=2; // Vacuum
      if(mi.GetBC(ind)==1)BC[ind]=1; // Reflective
      if(mi.GetBC(ind)==3)BC[ind]=3; // Periodic
      if(BC[ind]<1||BC[ind]>3){
        cout<<"Incorrect B. C. in SNT_system.\n";
        exit(0);
      };
      if(BC[ind]==3&&Ndim==3){
	cout<<"Periodic boundary condition cannot be adopted to 3D Calc.\n";
	exit(0);
      };
      if(Ndim==3&&j==1&&BC[ind]==1){
	cout<<"Not coded of reflective B.C. in the right boundary.\n";
	exit(0);
      };
      ind++;
    };
  };

  if(dsa)psys.PutCartMeshInfo(cm,geom);

  // for CMFD
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();
  CurFF.resize(1);
  for(int g=0;g<1;g++){
    CurFF[g].resize(zr+1); 
    for(int i=0;i<zr+1;i++){
      CurFF[g][i].resize(yr+1);
      for(int j=0;j<yr+1;j++){
        CurFF[g][i][j].resize(xr+1);
        for(int k=0;k<xr+1;k++){
          CurFF[g][i][j][k].resize(Ndim,0.);
	};
      };
    };
  };
};

void SNTSystem::CalSelfScatteringSource(real *Src,real *SSrc,real *flxmom, real *sigsself,int qid,int iss)
{
  vector<real> moment(plnum);
  for(int i=0;i<plnum;i++){
    moment[i]=quad[qid]->GetMoment(i,iss);
  };

  int ind3=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    real tmp=0.;
    int index=0;
    for(int l=0;l<=pl;l++){
      real sigs=sigsself[ind3++]; // sigs*vol/pi4
      int m1=-l;
      if(Ndim==2)m1=0;
      for(int m=m1;m<=l;m++){
        //tmp+=(sigs*flxmom[ind2]+Src[ind2])*quad[qid]->GetMoment(index++,iss);
        tmp+=(sigs*flxmom[ind2]+Src[ind2])*moment[index++];
	ind2++;
      };
    };
    SSrc[i]=tmp;
  };
};

void SNTSystem::PutSigmaForInnerIteration(int g,real *sigtvol, real *sigsself)
{
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    sigtvol[i]=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g)*vol;
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      sigsself[ind]=mesh[i].GetMed()->GetMacxs().GetSigs(l).get_dat(g,g)*volpi4;
      ind++;
    };
  };
};

void SNTSystem::InitialFluxGuessInnerIteration(int g,real *flxmom)
{
  int ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };
};

void SNTSystem::PutSourceInnerIteration(real *Src)
{
  int ind=0;
  for(int j=0;j<TotM;j++){
    for(int l=0;l<plnum;l++){
      Src[ind++]=mesh[j].GetSrcin(l);
    };
  };
};

void SNTSystem::RenewFluxMomentInnerIteration(int sn,int qid,real *flxmom,real *anflx)
{
  //vector<real> moment(sn*plnum);
  //int index=0;
  //for(int j=0;j<sn;j++){
  //  for(int l=0;l<plnum;l++){
  //    moment[index++]=quad[qid]->GetMoment(l,j);
  //  };
  //};

  real *finp=new real[plnum];
  int ind=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){finp[j]=0.;};
    //int index=0;
    for(int j=0;j<sn;j++){
      real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
      ind++;
      finp[0]+=tmp;
      for(int l=1;l<plnum;l++){
	//for(int l=0;l<plnum;l++){
	  finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	  //finp[l]+=tmp*moment[index++];
      };
    };
    for(int l=0;l<plnum;l++){
      flxmom[ind2]=real(finp[l]);
      ind2++;
    };
  };
  delete [] finp;
};

void SNTSystem::WriteAngularFlux(int sn,int g,real *anflx)
{
  if(g<=wrtflx_endgrp){
    int ind=0;
    for(int m=0;m<TotM;m++){
      if(wrtflx[m]){
        for(int is=0;is<sn;is++){
          aflux[m][g].put_data(is,anflx[ind]);
          ind++;
        };
      }else{
        ind+=sn;
      };
    };
  };
};

real SNTSystem::PutFluxAfterInnerIteration(int g,real *flxmom)
{
  int ind=0;
  real errf=0.;

  for(int i=0;i<TotM;i++){
    real tmp=flxmom[ind];
    real err=fabs(tmp/mesh[i].GetFlux(0).get_dat(g)-1.0);
    mesh[i].GetFlux(0).put_data(g,tmp);
    ind++;
    for(int j=1;j<plnum;j++){
      mesh[i].GetFlux(j).put_data(g,flxmom[ind]);
      ind++;
    };
    if(err>errf)errf=err;
  };
  return errf;
};

real SNTSystem::CalFluxXY(int g,int oiter,real epsif,bool angular_dependent_source)
{
  bool cmfd_on=false;
  if(cmfd&&oiter%opt.GetItcmfd()==0)cmfd_on=true;

  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src     =new real[TotM*plnum]; // External source (fission and down-scattering)
  real *SSrc    =new real[TotM]; // Source (self-scattering)
  real *anflx   =new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol =new real[TotM]; // Sigma_t * Volume
  real *sigsself=new real[TotM*(pl+1)]; // Self-scattering cross section * Volume / PI4
  real *flxd    =new real[TotM]; // Scalar flux in old inner iteration

  vector<real> siga(TotM);

  // Put sigma_t and sigma_s to temporal array
  PutSigmaForInnerIteration(g,sigtvol,sigsself);

  for(int i=0;i<TotM;i++){
    siga[i]=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
  };

  int xm =mi.GetXF();
  int ym =mi.GetYF();

  real flx;
  vector<real> fly(xm); // Initial flux in X axis
  vector<real> bflx(ym*sn,0.);  
  vector<real> bfly(xm*sn,0.);

  // Initial Flux guess
  InitialFluxGuessInnerIteration(g,flxmom);

  // Put fission & slowing down neutron source 
  PutSourceInnerIteration(Src);

  real outx,outy,cfl;
  real wx, wy, wxi, wyi;
  
  vector<real> mu_str(sn);
  vector<real> et_str(sn);
  vector<real> ww_str(sn);
  vector<real> muu_str(sn);
  vector<real> ett_str(sn);
  vector<real> muw_str(sn);
  vector<real> etw_str(sn);
  vector<int> xrsn_str(sn);
  vector<int> yrsn_str(sn);
  vector<int> xnref_str(sn,0);
  vector<int> ynref_str(sn,0);
  vector<int> xpref_str(sn,0);
  vector<int> ypref_str(sn,0);
  for(int is=0;is<sn;is++){
    real mu=quad[qid]->GetMu(is);
    real et=quad[qid]->GetEata(is);
    real ww=quad[qid]->GetOmega(is);
    mu_str[is]=mu;
    et_str[is]=et;
    ww_str[is]=ww;
    muu_str[is]=fabs(mu);
    ett_str[is]=fabs(et);
    muw_str[is]=mu*ww*PI2;
    etw_str[is]=et*ww*PI2;
    xrsn_str[is]=quad[qid]->GetXref(is);
    yrsn_str[is]=quad[qid]->GetYref(is);

    if((mu<0.&&BC[0]==1)||(mu>0.&&BC[1]==1)){xnref_str[is]=1;} // (reflective)
    else if((mu<0.&&BC[0]==3)||(mu>0.&&BC[1]==3)){xnref_str[is]=2;}; // (periodic)
    if((mu>0.&&BC[0]==1)||(mu<0.&&BC[1]==1)){xpref_str[is]=1;}
    else if((mu>0.&&BC[0]==3)||(mu<0.&&BC[1]==3)){xpref_str[is]=2;};

    if((et<0.&&BC[2]==1)||(et>0.&&BC[3]==1)){ynref_str[is]=1;}
    else if((et<0.&&BC[2]==3)||(et>0.&&BC[3]==3)){ynref_str[is]=2;};
    if((et>0.&&BC[2]==1)||(et<0.&&BC[3]==1)){ypref_str[is]=1;}
    else if((et>0.&&BC[2]==3)||(et<0.&&BC[3]==3)){ypref_str[is]=2;};
  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Check for sweeping from reflective (or periodic) boundary
  //
  //   The flux moment is NOT updated in the first inner
  //   iteretion since initial in-going angular flux from
  //   reflective boundary is zero and the obtained flux moment
  //   is far from the real solution.
  //
  // Check for reflective-reflective boundary
  //
  //   Convergence criteria for inner iteration is re-set to
  //   be severer.
  //
  bool spec_bc=false;
  if(BC[1]!=2||BC[3]!=2)spec_bc=true;
  if((BC[0]!=2&&BC[1]!=2)||(BC[2]!=2&&BC[3]!=2)){
    epsif*=0.1;
  };
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int itmax=itmax_inner;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    // for CMFD
    if(cmfd_on){
      for(int i=0;i<TotM;i++){
        curxn[i]=0.;
        curxp[i]=0.;
        curyn[i]=0.;
        curyp[i]=0.;
      };
    };

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    for(int is=0;is<sn;is++){

      real mu =mu_str[is];
      real et =et_str[is];
      real muu=muu_str[is];
      real ett=ett_str[is];
      real muw=muw_str[is];
      real etw=etw_str[is];
      int xrsn=xrsn_str[is];
      int yrsn=yrsn_str[is];
      int xnref=xnref_str[is];
      int ynref=ynref_str[is];
      int xpref=xpref_str[is];
      if(ypref_str[is]==1){
	int tmp=yrsn*xm;
	for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
      }else if(ypref_str[is]==2){
        int tmp=is*xm;
	for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
      }else{
	for(int i=0;i<xm;i++){fly[i]=0.;}; // vacuum
      };
      int iss=is;
      if(!opt.Forward())iss=quad[qid]->GetXYZref(is);
      int v1=xrsn*ym;
      int v2=is*ym;
      CalSelfScatteringSource(Src,SSrc,flxmom,sigsself,qid,iss);

      if(angular_dependent_source){
	for(int i=0;i<TotM;i++){
          SSrc[i]+=aflux[i][g].get_dat(is);
	};
      };

      for(int yy=0;yy<ym;yy++){
        int y=yy;
        if(et<0.)y=ym-yy-1;
        real yl=mi.GetFMeshL(1,y);

        flx=0.;
        if(xpref==1){flx=bflx[v1+y];}
	else if(xpref==2){flx=bflx[v2+y];};
        int tmp=y*xm;
        for(int xx=0;xx<xm;xx++){
	  int x=xx;
	  if(mu<0.)x=xm-xx-1;
 	  int m    =tmp+x;
  	  int index=m*sn+iss;
  	  real ss =SSrc[m];
	  real sigtv=sigtvol[m];
	  real xl=mi.GetFMeshL(0,x);
	  // AWDD
	  /*
	    real taux=siga[m]*xl/muu;
	    real tauy=siga[m]*yl/ett;
	    if(taux<0.1){
	      wx=0.5+taux/12.;
	    }else{
  	      wx=1./(1.-exp(-taux))-1./taux;
	    };
	    if(tauy<0.1){
	      wy=0.5+tauy/12.;
	    }else{
  	      wy=1./(1.-exp(-tauy))-1./tauy;
	    };
	  */
	  // DD
	  wx=0.5;
	  wy=0.5;
	  // 
          // Constant
	  //wx=1.0;
	  //wy=1.0;

	  wxi=1./wx;
	  wyi=1./wy;
          real tx=muu*yl*wxi;
          real ty=ett*xl*wyi;
	  real inx=flx;
	  real iny=fly[x];
	  // Sweep
          cfl=(tx*inx+ty*iny+ss)/(tx+ty+sigtv);
          outx=cfl*wxi-(1-wx)*wxi*inx;
          outy=cfl*wyi-(1-wy)*wyi*iny;

          while(outx<0.||outy<0.){
	    real tmp1=1.;
	    real tmp2=1.;
	    real tmp3=1.;
	    real tmp4=1.;
            if(outx<=0.){
              tmp1=wx;
              tmp3=0.;
	    };
            if(outy<=0.){
              tmp2=wy;
              tmp4=0.;
	    };
            cfl=(tx*inx*tmp1+ty*iny*tmp2+ss)/
                (tx*tmp3+ty*tmp4+sigtv);
            outx=tmp3*(cfl*wxi-(1-wx)*wxi*inx);
            outy=tmp4*(cfl*wyi-(1-wy)*wyi*iny);
	  };
  	  anflx[index]=cfl;
	  flx         =outx;
 	  fly[x]      =outy;

          // for CMFD
	  if(cmfd_on){
	    if(mu<0.){curxn[m]+=outx*muw;}
	    else{curxp[m]+=outx*muw;};
	    if(et<0.){curyn[m]+=outy*etw;}
	    else{curyp[m]+=outy*etw;};
          };

	};
	if(xnref>0){bflx[v2+y]=flx;}
      };
      if(ynref>0){
        int tmp=is*xm;
	for(int i=0;i<xm;i++){bfly[tmp+i]=fly[i];};
      };
    };

    // Renew flux moment
    if(it!=0||!spec_bc)RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);

    if(dsa)AccelerationByDSA(g,flxmom,flxd,sigsself,anflx,false);

    real errmax=0.;
    int ind=0;
    real flxmax=0.;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      //if(g==0)cout<<flxmom[ind]<<" ";
      if(err>errmax){
        errmax=err;
        flxmax=flxmom[ind];
      };
      ind+=plnum;
    };
    //cout<<g<<" "<<it<<" "<<errmax<<" "<<flxmax<<"\n";
    if(it!=0||!spec_bc){
      if(errmax<epsif){
        //cout<<"# ... group "<<g<<" ... "<<errmax<<" ("<<it<<")\n";
        break;
      };
    };
    if(it==itmax-1){
      cout<<"#   ... Not converged in group "<<g<<"\n";
      cout<<"#       maximum residual error : "<<errmax<<"\n";
      InnerConvergence=false;
    };
  };

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  WriteAngularFlux(sn,g,anflx);

  // for CMFD
  if(cmfd_on){
    if((g==0&&opt.Forward())||(g==grp-1&&!opt.Forward()))SetZeroCurFF();
    CalCoarseCur();
  };

  delete [] SSrc;
  delete [] Src;
  delete [] anflx; 
  delete [] flxmom;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigsself;

  return errf;
};

real SNTSystem::CalFluxXYZAllAbsorption(int g,int oiter,real epsif)
{
  bool cmfd_on=false;
  if(cmfd&&oiter%opt.GetItcmfd()==0)cmfd_on=true;

  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src     =new real[TotM*plnum]; // External source in plm
  real *SSrc    =new real[TotM]; // Source (self-scattering)
  real *anflx   =new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol =new real[TotM]; // Sigma_t * Volume
  real *sigsself=new real[TotM*(pl+1)]; // Self-scattering cross section * Volume / PI4
  real *flxd    =new real[TotM]; // Scalar flux in old inner iteration
  real *st      =new real[TotM]; // Sigma_t

  // Put sigma_t and sigma_s to temporal array
  PutSigmaForInnerIteration(g,sigtvol,sigsself);
  for(int i=0;i<TotM;i++){
    st[i]=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
  };

  int xm =mi.GetXF();
  int ym =mi.GetYF();
  int zm =mi.GetZF();
  int xym=xm*ym;

  real flx;
  vector<real> fly(xm); // Initial flux in X axis
  vector<real> flz(xm*ym); // Initial flux on XY plane
  vector<real> bflx(ym*zm*sn);  
  vector<real> bfly(xm*zm*sn);
  vector<real> bflz(xm*ym*sn);

  // Initial Flux guess
  InitialFluxGuessInnerIteration(g,flxmom);

  // Put fission & slowing down source in temporal array
  PutSourceInnerIteration(Src);

  real outx,outy,outz,cfl;

  real f1[]={1., 1., 1., 0.5, 0.5, 0.5};
  real f2[]={1., 1., 1., 0., 0., 0.};

  int itmax=1;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    // for CMFD
    if(cmfd_on){
      for(int i=0;i<TotM;i++){
        curxn[i]=0.;
        curxp[i]=0.;
        curyn[i]=0.;
        curyp[i]=0.;
        curzn[i]=0.;
        curzp[i]=0.;
      };
    };

    /*
    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };
    */

    for(int is=0;is<sn;is++){

      real mu =quad[qid]->GetMu(is);
      real et =quad[qid]->GetEata(is);
      real xi =quad[qid]->GetXi(is);
      real ww =quad[qid]->GetOmega(is);
      real muw=mu*ww*PI2;
      real etw=et*ww*PI2;
      real xiw=xi*ww*PI2;
      real muu=fabs(mu)*2.;
      real ett=fabs(et)*2.;
      real xii=fabs(xi)*2.;
      int xrsn =quad[qid]->GetXref(is);
      int yrsn =quad[qid]->GetYref(is);
      int zrsn =quad[qid]->GetZref(is);
      bool xnref = false;
      if(mu<0.&&BC[0]==1)xnref=true;
      bool ynref = false;
      if(et<0.&&BC[2]==1)ynref=true;
      bool xpref = false;
      if(mu>0.&&BC[0]==1)xpref=true;
      if(xi>0.&&BC[4]==1){
	int tmp=zrsn*xym;
	for(int i=0;i<xym;i++){flz[i]=bflz[tmp+i];};
      }else{
        for(int i=0;i<xym;i++){flz[i]=0.;};
      };
      int iss=is;
      if(!opt.Forward())iss=quad[qid]->GetXYZref(is);
      CalSelfScatteringSource(Src,SSrc,flxmom,sigsself,qid,iss);
      for(int zz=0;zz<zm;zz++){
	int z=zz;
	if(xi<0.)z=zm-zz-1;
	real zl=mi.GetFMeshL(2,z);
	if(et>0.&&BC[2]==1){
	  int tmp=yrsn*xm*zm+z*xm;
	  for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	}else{
	  for(int i=0;i<xm;i++){fly[i]=0.;};
	};
	int v1=xrsn*(ym*zm)+z*ym;
	int v2=is*(ym*zm)+z*ym;
        for(int yy=0;yy<ym;yy++){
	  int y=yy;
	  if(et<0.)y=ym-yy-1;
	  real yl=mi.GetFMeshL(1,y);
	  flx=0.;
	  if(xpref){flx=bflx[v1+y];};
	  //int xer=xedger[y];
	  //int xel=xedgel[y];
	  int xer=edge[0][y][z][1];
	  int xel=edge[0][y][z][0];
	  int xnum=xer-xel+1;
          for(int xx=0;xx<xnum;xx++){
	    int x=xx+xel;
	    if(mu<0.)x=xer-xx;
	    real xl=mi.GetFMeshL(0,x);
	    int m=meshid[z][y][x];
  	    int index=m*sn+iss;
  	    real ss=SSrc[m];
	    // Sweep
	    real inx=flx;
	    real iny=fly[x];
	    real inz=flz[y*xm+x];

	    // Step Characteristic
	    /*
	    real taux=st[m]*xl/mu;
	    real tauy=st[m]*yl/et;
	    real tauz=st[m]*zl/xi;
	    real wx=1./(1.-exp(-taux))-1./taux;
	    real wy=1./(1.-exp(-tauy))-1./tauy;
	    real wz=1./(1.-exp(-tauz))-1./tauz;
	    //real wx=1./(1.-etab.ee(-taux))-1./taux;
	    //real wy=1./(1.-etab.ee(-tauy))-1./tauy;
	    //real wz=1./(1.-etab.ee(-tauz))-1./tauz;
	    real wxi=0.5/wx;
	    // (Since `muu' is multiplied by 2.0.)
	    real wyi=0.5/wy;
	    real wzi=0.5/wz;
            real t1=muu*yl*zl*wxi;
            real t2=ett*xl*zl*wyi;
            real t3=xii*xl*yl*wzi;
	    */
	    // ++++++++++++++++++++++++++++++++++
	    // +++ Step differencing
	    /*
            real stv=sigtvol[m];
            real t1=muu*yl*zl*0.5; // since muu is multilplied by 2.0
            real t2=ett*xl*zl*0.5;
            real t3=xii*xl*yl*0.5;
	    cfl=(t1*inx+t2*iny+t3*inz+ss)/(t1+t2+t3+stv);
            outx=cfl;
	    outy=cfl;
	    outz=cfl;
	    */
	    // ++++++++++++++++++++++++++++++++++
	    // Diamond differencing

            real t1=muu*yl*zl;
            real t2=ett*xl*zl;
            real t3=xii*xl*yl;
	    real stv=sigtvol[m];
            cfl=(t1*inx+t2*iny+t3*inz+ss)/(t1+t2+t3+stv);
            real cfl2=cfl*2.;
            outx=cfl2-inx;
            outy=cfl2-iny;
            outz=cfl2-inz;
            while(outx<0.||outy<0.||outz<0.){
              int ifi=0;
              int ifj=0;
              int ifk=0;
              if(outx<=0.)ifi=1;
              if(outy<=0.)ifj=1;
              if(outz<=0.)ifk=1;
              cfl=(t1*inx*f1[ifi*3]+t2*iny*f1[1+ifj*3]+t3*inz*f1[2+ifk*3]+ss)/
                  (t1*f2[ifi*3]+t2*f2[1+ifj*3]+t3*f2[2+ifk*3]+stv);
              cfl2=cfl*2.;
              outx=f2[ifi*3]*(cfl2-inx);
              outy=f2[1+ifj*3]*(cfl2-iny);
              outz=f2[2+ifk*3]*(cfl2-inz);
	    };

	    // +++++++++++++++++++++++++++++++
  	    anflx[index]=cfl;
	    flx         =outx;
 	    fly[x]      =outy;
	    flz[y*xm+x] =outz;

	    // for CMFD
	    if(cmfd_on){
	      if(mu<0.){curxn[m]+=outx*muw;}
	      else{curxp[m]+=outx*muw;};
	      if(et<0.){curyn[m]+=outy*etw;}
	      else{curyp[m]+=outy*etw;};
	      if(xi<0.){curzn[m]+=outz*xiw;}
	      else{curzp[m]+=outz*xiw;};
	    };
	  };
	  if(xnref){
            bflx[v2+y]=flx;
          };
	};
	if(ynref){
          int tmp=is*(xm*zm)+z*xm;
	  for(int i=0;i<xm;i++){
            bfly[tmp+i]=fly[i];
          };
	};
      };
      if(xi<0.&&BC[4]==1){
        int tmp=is*xym;
	for(int i=0;i<xym;i++){
          bflz[tmp+i]=flz[i];
        };
      };
    };

    // Renew flux moment
    RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);

    if(dsa)AccelerationByDSA(g,flxmom,flxd,sigsself,anflx,cmfd_on);

    /*
    real errmax=0.;
    int ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    //cout<<it<<" "<<errmax<<"\n";
    if(errmax<epsif){
      //cout<<" ... group "<<g<<" ... "<<errmax<<"\n";
      break;
    };
    if(it==itmax-1){
      //cout<<"#   ... Not converged in group "<<g<<" ("<<errmax<<"/"<<epsif<<")\n";
      InnerConvergence=false;
    };
    */
  };

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  WriteAngularFlux(sn,g,anflx);

  // for CMFD
  if(cmfd_on){
    if((g==0&&opt.Forward())||(g==grp-1&&!opt.Forward()))SetZeroCurFF();
    CalCoarseCur();
  };

  delete [] SSrc;
  delete [] Src;
  delete [] anflx; 
  delete [] flxmom;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigsself;
  delete [] st;

  return errf;
};

real SNTSystem::CalFluxXYZ(int g,int oiter,real epsif)
{
  bool cmfd_on=false;
  if(cmfd&&oiter%opt.GetItcmfd()==0)cmfd_on=true;

  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  // +++++ Galerkin quadrature test ++++++
  /*
  vector<real> a_l(pl+1);
  for(int l=0;l<=pl;l++){
    //a_l[l]=1.; // Delta scattering
    real xs0=med[0].GetMacxs().GetData2d(sigs,0).get_dat(0,0);
    real xs1=med[0].GetMacxs().GetData2d(sigs,l).get_dat(0,0);
    a_l[l]=xs1/(2.*l+1.)/xs0;
  };
  GroupData2D mmat=quad[qid]->GetGalerkinMatrix();
  GroupData2D sigmat(sn,sn);
  sigmat.set_zero();
  for(int i=0;i<sn;i++){
    sigmat.put_data(i,i,med[0].GetMacxs().GetData2d(sigs,0).get_dat(0,0)*a_l[quad[qid]->GetLArrayGMatrix(i)]);
  };
  GroupData2D dmat=quad[qid]->GetGalerkinMatrixInverse();
  GroupData2D pmat=mmat*sigmat*dmat;
  */
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++

  real *Src     =new real[TotM*plnum]; // External source in plm
  real *SSrc    =new real[TotM]; // Source (self-scattering)
  real *anflx   =new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol =new real[TotM]; // Sigma_t * Volume
  real *sigsself=new real[TotM*(pl+1)]; // Self-scattering cross section * Volume / PI4
  real *flxd    =new real[TotM]; // Scalar flux in old inner iteration
  real *st      =new real[TotM]; // Sigma_t

  // Put sigma_t and sigma_s to temporal array
  PutSigmaForInnerIteration(g,sigtvol,sigsself);
  for(int i=0;i<TotM;i++){
    st[i]=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
    //st[i]=mesh[i].GetMed()->GetMacxs().GetSiga().get_dat(g); // AWDD
  };

  int xm =mi.GetXF();
  int ym =mi.GetYF();
  int zm =mi.GetZF();
  int xym=xm*ym;

  real flx;
  vector<real> fly(xm); // Initial flux in X axis
  vector<real> flz(xm*ym); // Initial flux on XY plane
  vector<real> bflx(ym*zm*sn);  
  vector<real> bfly(xm*zm*sn);
  vector<real> bflz(xm*ym*sn);

  // Initial Flux guess
  InitialFluxGuessInnerIteration(g,flxmom);

  // Put fission & slowing down source in temporal array
  PutSourceInnerIteration(Src);

  real outx,outy,outz,cfl;

  real f1[]={1., 1., 1., 0.5, 0.5, 0.5};
  real f2[]={1., 1., 1., 0., 0., 0.};

  int itmax=itmax_inner;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    // for CMFD
    if(cmfd_on){
      for(int i=0;i<TotM;i++){
        curxn[i]=0.;
        curxp[i]=0.;
        curyn[i]=0.;
        curyp[i]=0.;
        curzn[i]=0.;
        curzp[i]=0.;
      };
    };

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    for(int is=0;is<sn;is++){

      real mu =quad[qid]->GetMu(is);
      real et =quad[qid]->GetEata(is);
      real xi =quad[qid]->GetXi(is);
      real ww =quad[qid]->GetOmega(is);
      real muw=mu*ww*PI2;
      real etw=et*ww*PI2;
      real xiw=xi*ww*PI2;
      real muu=fabs(mu)*2.;
      real ett=fabs(et)*2.;
      real xii=fabs(xi)*2.;
      int xrsn =quad[qid]->GetXref(is);
      int yrsn =quad[qid]->GetYref(is);
      int zrsn =quad[qid]->GetZref(is);
      bool xnref = false;
      if(mu<0.&&BC[0]==1)xnref=true;
      bool ynref = false;
      if(et<0.&&BC[2]==1)ynref=true;
      bool xpref = false;
      if(mu>0.&&BC[0]==1)xpref=true;
      if(xi>0.&&BC[4]==1){
	int tmp=zrsn*xym;
	for(int i=0;i<xym;i++){flz[i]=bflz[tmp+i];};
      }else{
        for(int i=0;i<xym;i++){flz[i]=0.;};
      };
      int iss=is;
      if(!opt.Forward())iss=quad[qid]->GetXYZref(is);


      // USUAL SNT
      CalSelfScatteringSource(Src,SSrc,flxmom,sigsself,qid,iss);
      // GALERKIN OPTION

      /*
      int ind2=0;
      for(int i=0;i<TotM;i++){
	SSrc[i]=Src[ind2];
	real tmp=0.;
	for(int l=0;l<sn;l++){
	  tmp+=pmat.get_dat(is,l)*anflx[i*sn+l];
	};
	SSrc[i]+=tmp*GetMesh(i).GetVolume();
        for(int l=0;l<=pl;l++){
          for(int m=-l;m<=l;m++){
    	    ind2++;
	  };
        };
      };
      */




      for(int zz=0;zz<zm;zz++){
	int z=zz;
	if(xi<0.)z=zm-zz-1;
	real zl=mi.GetFMeshL(2,z);
	if(et>0.&&BC[2]==1){
	  int tmp=yrsn*xm*zm+z*xm;
	  for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	}else{
	  for(int i=0;i<xm;i++){fly[i]=0.;};
	};
	int v1=xrsn*(ym*zm)+z*ym;
	int v2=is*(ym*zm)+z*ym;
        for(int yy=0;yy<ym;yy++){
	  int y=yy;
	  if(et<0.)y=ym-yy-1;
	  real yl=mi.GetFMeshL(1,y);
	  flx=0.;
	  if(xpref){flx=bflx[v1+y];};
	  //int xer=xedger[y];
	  //int xel=xedgel[y];
	  int xer=edge[0][y][z][1];
	  int xel=edge[0][y][z][0];
	  int xnum=xer-xel+1;
          for(int xx=0;xx<xnum;xx++){
	    int x=xx+xel;
	    if(mu<0.)x=xer-xx;
	    real xl=mi.GetFMeshL(0,x);
	    int m=meshid[z][y][x];
  	    int index=m*sn+iss;
  	    real ss=SSrc[m];
	    // Sweep
	    real inx=flx;
	    real iny=fly[x];
	    real inz=flz[y*xm+x];

	    // +++ Step Characteristic +++
	    /*
	    real taux=st[m]*xl/(muu*0.5);
	    real tauy=st[m]*yl/(ett*0.5);
	    real tauz=st[m]*zl/(xii*0.5);
	    real wx=1./(1.-exp(-taux))-1./taux;
	    real wy=1./(1.-exp(-tauy))-1./tauy;
	    real wz=1./(1.-exp(-tauz))-1./tauz;
	    //real wx=1./(1.-etab.ee(-taux))-1./taux;
	    //real wy=1./(1.-etab.ee(-tauy))-1./tauy;
	    //real wz=1./(1.-etab.ee(-tauz))-1./tauz;
	    real wxi=1./wx; 
	    real wyi=1./wy;
	    real wzi=1./wz;
            real t1=0.5*muu*yl*zl*wxi; // (Since `muu' is multiplied by 2.0.)
            real t2=0.5*ett*xl*zl*wyi;
            real t3=0.5*xii*xl*yl*wzi;
	    real stv=sigtvol[m];
            cfl=(t1*inx+t2*iny+t3*inz+ss)/(t1+t2+t3+stv);
            outx=cfl*wxi-(1-wx)*wxi*inx;
            outy=cfl*wyi-(1-wy)*wyi*iny;
            outz=cfl*wzi-(1-wz)*wzi*inz;
            while(outx<0.||outy<0.||outz<0.){
              int ifi=0;
              int ifj=0;
              int ifk=0;
              if(outx<=0.)ifi=1;
              if(outy<=0.)ifj=1;
              if(outz<=0.)ifk=1;
              cfl=(t1*inx*f1[ifi*3]+t2*iny*f1[1+ifj*3]+t3*inz*f1[2+ifk*3]+ss)/
                  (t1*f2[ifi*3]+t2*f2[1+ifj*3]+t3*f2[2+ifk*3]+stv);
              outx=f2[ifi*3]*(cfl*wxi-(1-wx)*wxi*inx);
              outy=f2[1+ifj*3]*(cfl*wyi-(1-wy)*wyi*iny);
              outz=f2[2+ifk*3]*(cfl*wzi-(1-wz)*wzi*inz);
	    };
	    */
	    // ++++++++++++++++++++++++++++
	    // +++ Step differencing +++
	    /*
            real stv=sigtvol[m];
            real t1=muu*yl*zl*0.5; // since muu is multilplied by 2.0
            real t2=ett*xl*zl*0.5;
            real t3=xii*xl*yl*0.5;
	    cfl=(t1*inx+t2*iny+t3*inz+ss)/(t1+t2+t3+stv);
            outx=cfl;
	    outy=cfl;
	    outz=cfl;
	    */
	    // +++++++++++++++++++++++++
	    // +++ Diamond differencing ++++
            real t1=muu*yl*zl;
            real t2=ett*xl*zl;
            real t3=xii*xl*yl;
	    real stv=sigtvol[m];
            cfl=(t1*inx+t2*iny+t3*inz+ss)/(t1+t2+t3+stv);
            real cfl2=cfl*2.;
            outx=cfl2-inx;
            outy=cfl2-iny;
            outz=cfl2-inz;
            while(outx<0.||outy<0.||outz<0.){
              int ifi=0;
              int ifj=0;
              int ifk=0;
              if(outx<=0.)ifi=1;
              if(outy<=0.)ifj=1;
              if(outz<=0.)ifk=1;
	      cfl=(t1*inx*f1[ifi*3]+t2*iny*f1[1+ifj*3]+t3*inz*f1[2+ifk*3]+ss)/
                  (t1*f2[ifi*3]+t2*f2[1+ifj*3]+t3*f2[2+ifk*3]+stv);
              cfl2=cfl*2.;
              outx=f2[ifi*3]*(cfl2-inx);
              outy=f2[1+ifj*3]*(cfl2-iny);
              outz=f2[2+ifk*3]*(cfl2-inz);
	    };
	    
	    // +++++++++++++++++++++++++++++
  	    anflx[index]=cfl;
	    flx         =outx;
 	    fly[x]      =outy;
	    flz[y*xm+x] =outz;

	    // for CMFD
	    if(cmfd_on){
	      if(mu<0.){curxn[m]+=outx*muw;}
	      else{curxp[m]+=outx*muw;};
	      if(et<0.){curyn[m]+=outy*etw;}
	      else{curyp[m]+=outy*etw;};
	      if(xi<0.){curzn[m]+=outz*xiw;}
	      else{curzp[m]+=outz*xiw;};
	    };
	  };
	  if(xnref){
            bflx[v2+y]=flx;
          };
	};
	if(ynref){
          int tmp=is*(xm*zm)+z*xm;
	  for(int i=0;i<xm;i++){
            bfly[tmp+i]=fly[i];
          };
	};
      };
      if(xi<0.&&BC[4]==1){
        int tmp=is*xym;
	for(int i=0;i<xym;i++){
          bflz[tmp+i]=flz[i];
        };
      };
    };

    // Renew flux moment
    RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);

    if(dsa)AccelerationByDSA(g,flxmom,flxd,sigsself,anflx,cmfd_on);

    real errmax=0.;
    int ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    //cout<<it<<" "<<errmax<<"\n";
    if(errmax<epsif){
      //cout<<" ... group "<<g<<" ... "<<it<<"  ("<<errmax<<")\n";
      break;
    };
    if(it==itmax-1){
      cout<<"#   ... Not converged in group "<<g<<" ("<<errmax<<"/"<<epsif<<")\n";
      InnerConvergence=false;
    };
  };

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  WriteAngularFlux(sn,g,anflx);

  // for CMFD
  if(cmfd_on){
    if((g==0&&opt.Forward())||(g==grp-1&&!opt.Forward()))SetZeroCurFF();
    CalCoarseCur();
  };

  delete [] SSrc;
  delete [] Src;
  delete [] anflx; 
  delete [] flxmom;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigsself;
  delete [] st;

  return errf;
};

void SNTSystem::AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself,real *anflx,bool cmfd_on)
{
  vector<real> fltt(TotM);
  // *** Source calculation
  int ind=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    fltt[i]=(flxmom[ind]-flxd[i])*sigsself[ind2]; 
    ind2+=pl+1;
    ind+=plnum;
  };

  // *** positive source term

  for(int i=0;i<TotM;i++){
    real tmp=fltt[i];
    //if(tmp<0.)tmp=0.;
    psys.GetMesh(i).PutSrcin(tmp);  
  };
  psys.CalFluxAutoConv(g,0.01);
  //psys.CalFlux(g,0,0.01);
  ind=0;
  for(int i=0;i<TotM;i++){
    real tmp1=flxmom[ind];
    real tmp2=psys.GetMesh(i).GetFlux(0).get_dat(g);
    real tmp3=tmp1+tmp2;
    if(tmp3>0.){
      flxmom[ind]=tmp3;
      // for CMFD
      if(cmfd_on){
        tmp3/=tmp1;
        curxp[i]*=tmp3;
        curxn[i]*=tmp3;
        curyp[i]*=tmp3;
        curyn[i]*=tmp3;
        curzp[i]*=tmp3;
        curzn[i]*=tmp3;
      };
    };
    ind+=plnum;
  };
  // *** negative source term
  /*
  for(int i=0;i<TotM;i++){
    real tmp=-fltt[i];
    if(tmp<0)tmp=0.;
    psys.GetMesh(i).PutSrcin(tmp);  
  };
  psys.CalFluxAutoConv(g,0.01);
  ind=0;
  for(int i=0;i<TotM;i++){
    flxmom[ind]-=psys.GetMesh(i).GetFlux(0).get_dat(g);
    ind+=plnum;
  };
  */
};

void SNTSystem::SetInitialFlux()
{
  SetArray();

  if(dsa){
  GeneralOption optdif;
  optdif.PutEpsf(5e-2);
  optdif.PutEpsk(5e-4);
  optdif.PutEpss(5e-3);
  optdif.PutOutitermax(30);
  if(!opt.Forward())optdif.PutAdjointCal();
  psys.PutGeneralOption(optdif);
  psys.CalCoef();
  psys.CalIgen("extrapolation");
  for(int i=0;i<TotM;i++){
    mesh[i].GetFlux().copy(psys.GetMesh(i).GetFlux());
  };
  }else{
    SetInitialFlatFlux();
  };
};

real SNTSystem::CalFluxGeneral(int ginp,real cin,int iter)
{
  if(Ndim==3)return CalFluxXYZ(ginp,iter,cin);
  if(Ndim==2)return CalFluxXY(ginp,iter,cin);
  return 0.;
};

void SNTSystem::AddMedium(Medium inp)
{
  GeneralSystem::AddMedium(inp);
  Medium inp2=inp;
  /*
  for(int g=0;g<grp;g++){
    //real org=inp2.GetMacxs().GetData1d(sigtr).get_dat(g);
    //inp2.GetMacxs().GetData1d(d).put_data(g,0.281/org);
    real org=inp2.GetMacxs().GetData1d(d).get_dat(g);
    inp2.GetMacxs().GetData1d(d).put_data(g,0.281*3.*org);
  };
  */
  psys.AddMedium(inp2); // specific for snt
  //psys.AddMedium(inp); // specific for snt
}

// **********************************************
// * For perturbation calculation               *
// **********************************************

real SNTSystem::CalReactivity(SNTSystem *sec,real kunp,real kp,bool pr)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in SNTSystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in SNTSystem::CalReactivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };


  real *yld=new real[grp];
  real *abs=new real[grp];
  vector< vector<real> > sct(plnum);
  vector< vector<real> > leak(plnum);
  for(int i=0;i<plnum;i++){
    sct[i].resize(grp,0.);
    leak[i].resize(grp,0.);
  };

  real ip=CalPerturbDenominator(sec);
  if(ip==0.){
    cout<<"Error in SNTsystem::CalSensitivity.\n";
    cout<<"Perturbation denominator is zero.\n";
    exit(0);
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
    CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
    for(int j=0;j<plnum;j++){
      sct[j][i]=scttmp[j];
      leak[j][i]=leaktmp[j];
    };
  };
  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;

  real yldsum=0.;
  real abssum=0.;
  real tsctsum=0.;
  real tleaksum=0.;
  real *sctsum=new real[plnum];
  real *leaksum=new real[plnum];
  for(int i=0;i<plnum;i++){
    sctsum[i]=0.; leaksum[i]=0.;
  };

  if(pr){
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scattering  ";
    cout<<"Leakage     ";
    cout<<"Total\n";
  };
  for(int i=0;i<grp;i++){
    yld[i]/=ip;
    abs[i]/=ip;
    real sl=0.;
    real le=0.;
    for(int j=0;j<plnum;j++){
      sct[j][i]/=ip; sctsum[j]+=sct[j][i]; sl+=sct[j][i];
      leak[j][i]/=ip; leaksum[j]+=leak[j][i]; le+=leak[j][i];
    };
    yldsum+=yld[i];
    abssum+=abs[i];
    tsctsum+=sl;
    tleaksum+=le;
    if(pr){
      cout.width(3);
      cout<<i<<"   ";
      cout.setf(ios::scientific);
      cout.precision(4);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<le/leth;
      cout<<" "<<(yld[i]+abs[i]+sl+le)/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  real tot=yldsum+abssum+tsctsum+tleaksum;
  if(pr){
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"#\n";
    cout<<"# Yield       : "<<yldsum<<"\n";
    cout<<"# Absorption  : "<<abssum<<"\n";
    real tmp=yldsum+abssum;
    int id=0;
    for(int l=0;l<=pl;l++){
      int mst=-l; // snt
      if(Ndim==2)mst=0; // snt
      for(int m=mst;m<=l;m++){ // snt
        cout<<"#  ("<<l<<","<<m<<")th scattering  : "<<sctsum[id]<<"\n";
	if(l==0){tmp+=sctsum[id];}else{tleaksum+=sctsum[id];};
        id++;
      };
    };
    cout<<"# Scattering    : "<<tsctsum<<"\n";
    cout<<"# (Non-Leakage) : "<<tmp<<"\n";
    id=0;
    for(int l=1;l<=pl;l++){
      int mst=-1; // snt
      if(Ndim==2)mst=0; // snt
      for(int m=mst;m<=l;m++){ // snt
        cout<<"#  ("<<l<<","<<m<<")th leakage     : "<<leaksum[id]<<"\n";
        id++;
      };
    };
    cout<<"#     Higher leakage : "<<leaksum[plnum-1]<<"\n";
    cout<<"# Leakage      : "<<tleaksum<<"\n";
    cout.unsetf(ios::scientific);
  };

  cout.setf(ios::scientific);
  cout.precision(4);

  if(pr)cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
  if(pr)cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
  cout.unsetf(ios::scientific);

  delete [] sctsum;
  delete [] leaksum;
  delete [] yld;
  delete [] abs;

  return tot;
};

void SNTSystem::CalPerturbLeakScat(SNTSystem *sec,bool *flag,int g,real *sct,real *leak)
{
  // The absorption component is calculated by the method [CalPerturbAbsorptionTerm]
  // implemented in [GeneralSystem].
  
  for(int i=0;i<plnum;i++){
    sct[i]=0.; leak[i]=0.;
  };

  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      real dtot=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt().get_dat(g)
                        -mesh[m].GetMed()->GetMacxs().GetSigt().get_dat(g);
      real dabs=sec->GetMesh(m).GetMed()->GetMacxs().GetSiga().get_dat(g)
                        -mesh[m].GetMed()->GetMacxs().GetSiga().get_dat(g);
      real dsct=dtot-dabs; // The 0th-order scattering cross section perturbation

      // ----------------------------------------------------------------
      // (Perturbation in total cross section)
      //
      
      // ... effect of the 0th-order scattering on the 0th-order angular moment
      //     categorized as the scattering component
      sct[0]-=dsct*sec->GetMesh(m).GetFlux(0).get_dat(g)*
                           mesh[m].GetFlux(0).get_dat(g)*vol;

      // ... effect of the 0th-order scattering on the higher-order angular moment
      //     categorized as the leakage component
      int id=0;
      for(int l=1;l<=pl;l++){
	int mst=-l; // snt
	if(Ndim==2)mst=0; // snt 
 	for(int mm=mst;mm<=l;mm++){ // snt
	  leak[id]-=(2.*l+1.)*dtot
                   *sec->GetMesh(m).GetFlux(id+1).get_dat(g)
                           *mesh[m].GetFlux(id+1).get_dat(g)*vol;
	  id++;
	};
      };

      // ... effect of the total cross section perturbation
      //     on the higher-order angular moment
      //
      // flp : [Angular flux] - [Component expanded by the [pl]th-order Legendre moment]
      //
      if(g<=wrtflx_endgrp&&GetWrtflx(m)&&sec->GetWrtflx(m)){
        real tmp2=0.;
        int sn=GetQuadrature(g)->GetSN(); // snt
        for(int i=0;i<sn;i++){
	  real om =GetQuadrature(g)->GetOmega(i);
          real flp=sec->GetAFlux(m,g).get_dat(i);
          int id=0;
	  for(int l=0;l<=pl;l++){
            int mst=-l; // snt
	    if(Ndim==2)mst=0; // snt
	    for(int mm=mst;mm<=l;mm++){ // snt
      	      real mom=GetQuadrature(g)->GetMoment(id,i);
              flp-=(2.*l+1.)/PI4*mom*sec->GetMesh(m).GetFlux(id).get_dat(g);
              id++;
            };
          };
	  tmp2+=om*flp*aflux[m][g].get_dat(i);
        };
	leak[plnum-1]-=PI2*PI4*dtot*tmp2*vol;
      };
      // ------------------------------------------------------------------------

      
      // ------------------------------------------------------------------------      
      // (Perturbation in scattering cross section in the scattering source term)
      //
      //   Results are added to the scattering component

      id=0;
      for(int l=0;l<=pl;l++){
	int mst=-l; // snt
	if(Ndim==2)mst=0; // snt
	for(int mm=mst;mm<=l;mm++){ // snt
  	  real tmp=0.;
	  for(int gg=g;gg<grp;gg++){
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    tmp+=dsigs*mesh[m].GetFlux(id).get_dat(gg);
	  };
  	  sct[id]+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	  id++;
	};
      };
      // ------------------------------------------------------------------------            

    };
  };
};

void SNTSystem::CalSensitivity(SNTSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  if(GetGeneralOption().Forward()){
    cout<<"Error in SNTSystem::CalSensitivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in SNTSystem::CalSensitivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  cout<<"******\n";
  cout<<" Caution !!!\n";
  cout<<" Sensitivity calculation with SNT has not been well validated.\n";

  real *nsforg=new real[nmed];
  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  real *tot1org=new real[nmed];
  real *sigsorg=new real[nmed*grp];
  real *sigs1org=new real[nmed*grp];

  real *fiss_frac=new real[nmed];

  real delta=1.;

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o

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
            totorg[j]=sec->GetMed(j).GetSigt(0).get_dat(i);
            tot1org[j]=sec->GetMed(j).GetSigt(1).get_dat(i);
   	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigf).get_dat(i);
	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den[j]*micnu*micsigf);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den[j]*micsigf);
	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetSigt(l).add_data(i,den[j]*micsigf);
	    };
          };
        };
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int k=0;k<plnum;k++){
          re+=scttmp[k]+leaktmp[k];
        };
	//re-=scttmp[0]; // L
        //re+=scttmp[0]; // NL

        re+=CalPerturbYieldTerm(sec,flag,i,k2);
        re+=CalPerturbAbsorptionTerm(sec,flag,i);
        re*=inv_ip;
        cout<<"  "<<re<<"\n";
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    for(int l=1;l<=pl;l++){
              sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	    };
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
  	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigf).get_dat(i);
	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den[j]*micnu*micsigf);
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
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
          tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigc).get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den[j]*micsigc);
	  for(int l=0;l<=pl;l++){
	    sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den[j]*micsigc);
	  };
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      //re+=scttmp[0]; // NL
      //re-=scttmp[0]; // L

      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      re*=inv_ip;
      cout<<"  "<<re<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  for(int l=1;l<=pl;l++){
            sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	  };
        };
      };
    };

    // Scattering P0
    for(int ii=0;ii<3;ii++){
      int ndd=2;
      if(nidendf==125&&ii==0)ndd=3;
      cout<<ndd<<"\n";
      cout<<nidendf<<"\n";
      //
      if(ii==0){cout<<2<<"\n";}
      else if(ii==1){cout<<4<<"\n";}
      else {cout<<16<<"\n";};
      //
      for(int i=0;i<grp;i++){
        int st=i;
        if(ndd==3)st=0;
	for(int k=st;k<grp;k++){
	  if(ii!=0||nidendf<1000||(ii==0&&k<=i+2)){
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
                totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
                tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
  	        sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
	        sigs1org[j]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
                real micsigs=0.;
                real micsigs1=0.;
                real s0=0.;
                real s1=0.;
		switch(ii){
		case 0:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	          micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_dat(i,k);
		  s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
		  s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_sumx().get_dat(i);
	  	  break;
	        case 1:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_dat(i,k);
	          micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(1).get_dat(i,k);
		  s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_sumx().get_dat(i);
		  s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(1).get_sumx().get_dat(i);
		  break;
	        case 2:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_dat(i,k);
	          //micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(1).get_dat(i,k);
		  break;
	        };
	        real a1=0.;
		//real s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
		//eal s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_sumx().get_dat(i);
		if(s0!=0.)a1=s1/s0;
	        //real a1;
	        //if(micsigs!=0.){a1=micsigs1/micsigs;}
	        //else{a1=0.;};
                if(ii!=0.)micsigs=delta;
	        real dtot;
	        if(ii==2){
  	  	  dtot=micsigs*0.5;
	        }else{
   	          dtot=micsigs;
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den[j]*micsigs);
	        sec->GetMed(j).GetMacxs().GetSigs(1).add_data(i,k,den[j]*a1*micsigs);
	        for(int l=0;l<=pl;l++){
	          sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den[j]*dtot);
	        };
		//if(pl==0&&ii==0){
		//  real vmu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(mu).get_dat(i);
		//  real tmp=-den[j]*vmu*micsigs;
		//  sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,tmp);
    	        //  sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,i,tmp);
		//};
              };
            };
            real re=0.;
            CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
            for(int l=0;l<plnum;l++){
  	      re+=scttmp[l]+leaktmp[l];
            };
	    //re+=scttmp[0]; // NL
	    //re-=scttmp[0]; // L

            re*=inv_ip;
            cout<<"  "<<re<<"\n";
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
                sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	        for(int l=1;l<=pl;l++){
                  sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	        sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigs1org[j]);
	      };
	    };
	  }else{
	    cout<<"  0.0\n";
	  };
	};
      };
    };

    // Elastic-p1(mu)
    cout<<1<<"\n";
    cout<<nidendf<<"\n";
    cout<<"  251\n";
    if(pl!=0){
     for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  int maxg=grp;
	  if(nidendf>1000&&i+3<maxg)maxg=i+3;
          int st=i;
          if(nidendf==125)st=0;
	  for(int k=st;k<maxg;k++){
	    sigs1org[j*grp+k]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
	    real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	    sec->GetMed(j).GetSigs(1).add_data(i,k,den[j]*delta*micsigs);
	  };
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      //re+=scttmp[0]; // NL
      //re-=scttmp[0]; // L

      re*=inv_ip;
      cout<<"  "<<re*3.<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          int maxg=grp;
	  if(nidendf>1000&&i+3<maxg)maxg=i+3;
	  for(int k=i;k<maxg;k++){
	    sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigs1org[j*grp+k]);
	  };
	};
      };
     };
    }else{
      for(int i=0;i<grp;i++){
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	    sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,i);
	    real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigel).get_dat(i);
	    sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,-den[j]*delta*micsigs);
	    sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,i,-den[j]*delta*micsigs);
	  };
	};
        real re=0.;
	CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
	re+=scttmp[0]+leaktmp[0];
	//re+=scttmp[0]; // NL
	//re-=scttmp[0]; // L
	re*=inv_ip;
	cout<<" "<<re<<"\n";
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,i,sigsorg[j]);
	  };
	};
      };
    };
  };

  delete [] scttmp;
  delete [] leaktmp;

  delete [] fiss_frac;
  delete [] flag;
  delete [] sigsorg;
  delete [] sigs1org;
  delete [] absorg;
  delete [] nsforg;
  delete [] totorg;
  delete [] tot1org;
};

SensitivityData SNTSystem::CalSensitivityNew(SNTSystem *sec,real keff,int nucnum,int *nucid)
{
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
  real *sigsorg=new real[nmed*grp];
  real *sigs1org=new real[nmed*grp];

  real *fiss_frac=new real[nmed];

  real delta=1.;

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o

  bool *flag=new bool[TotM];

  real ip=0.;
  if(ip_input){
    ip=ip_input_value;
  }else{
    ip=CalPerturbDenominator(sec);
    //if(fission_spectrum_matrix){
    //  ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
    //};
    if(ip==0.){
      cout<<"# Error in SNTsystem::CalSensitivity.\n";
      cout<<"# Perturbation denominator is zero.\n";
      exit(0);
    };
  };

/*
  real ip=CalPerturbDenominator(sec);
  if(ip==0.){
    cout<<"# Error in SNTsystem::CalSensitivityNew.\n";
    cout<<"# Perturbation denominator is zero.\n";
    exit(0);
  };
*/
  //real inv_ip=1./ip;

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
              totorg[j]=sec->GetMed(j).GetSigt(0).get_dat(i);
              tot1org[j]=sec->GetMed(j).GetSigt(1).get_dat(i);
   	      real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigf).get_dat(i);
	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
              sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
	      for(int l=0;l<=pl;l++){
	        sec->GetMed(j).GetSigt(l).add_data(i,den*micsigf);
	      };
            };
          };
          real re=0.;
          CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
          for(int k=0;k<plnum;k++){
            re+=scttmp[k]+leaktmp[k];
          };
          re+=CalPerturbYieldTerm(sec,flag,i,keff);
          re+=CalPerturbAbsorptionTerm(sec,flag,i);
          re*=factor;
          sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
  	    if(med[j].ExistNuclide(nid)){
    	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	      sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
              sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	      for(int l=1;l<=pl;l++){
                sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	      };
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
  	      real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigf).get_dat(i);
	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
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
          re+=CalPerturbYieldTerm(sec,flag,i,keff);
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
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
          tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigc).get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
	  for(int l=0;l<=pl;l++){
	    sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*micsigc);
	  };
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      re*=factor;
      sns1d.put_data(i,re);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  for(int l=1;l<=pl;l++){
            sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	  };
        };
      };
    };
    sens.PutSensitivity1D(nid,102,sns1d);

    // Scattering P0
    for(int ii=0;ii<3;ii++){
      int mt=2;
      if(ii==1)mt=4;
      if(ii==2)mt=16;
      sns2d.set_zero();
      for(int i=0;i<grp;i++){
	for(int k=i;k<grp;k++){
	  if(ii!=0||nid<100000||(ii==0&&k<=i+2)){
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
                real den=med[j].GetNuclide(nid).GetDensity();
                totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
                tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
  	        sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
	        sigs1org[j]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
                real micsigs=0.;
                real micsigs1=0.;
                real s0=0.;
                real s1=0.;
		switch(ii){
		case 0:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	          micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_dat(i,k);
		  s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
		  s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_sumx().get_dat(i);
	  	  break;
	        case 1:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_dat(i,k);
	          micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(1).get_dat(i,k);
		  s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_sumx().get_dat(i);
		  s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(1).get_sumx().get_dat(i);
		  break;
	        case 2:
	          micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_dat(i,k);
                  absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);		    
		  break;
	        };
                if(ii!=0.)micsigs=delta; 
                // 'absolute' sensitivity for inelastic & (n,2n) 
	        sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
	        real factor=1.;
	        if(ii==2)factor=0.5;
	        for(int l=0;l<=1;l++){
	          sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*micsigs*factor);
	        };
	        real a1=0.;
		if(s0!=0.)a1=s1/s0;
	        sec->GetMed(j).GetMacxs().GetSigs(1).add_data(i,k,den*a1*micsigs);
                if(ii==2)sec->GetMed(j).GetMacxs().GetSiga().add_data(i,-den*micsigs*factor);
              };
            };

            real re=0.;
            re+=CalPerturbAbsorptionTerm(sec,flag,i);	      
            CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
            for(int l=0;l<plnum;l++){
  	      re+=scttmp[l]+leaktmp[l];
            };
            re*=factor;
            sns2d.put_data(i,k,re);
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
                sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	        for(int l=1;l<=pl;l++){
                  sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[j]);
	        sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigs1org[j]);
      	        if(ii==2)sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	      };
	    };
	  };
	};
      };
      sens.PutSensitivity2D(nid,mt,sns2d);
    };

    // Elastic-p1(mu)
    if(pl!=0){
     for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
	  int maxg=grp;
	  if(nid>1000&&i+3<maxg)maxg=i+3;
	  for(int k=i;k<maxg;k++){
	    sigs1org[j*grp+k]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
	    real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	    sec->GetMed(j).GetSigs(1).add_data(i,k,den*delta*micsigs);
	  };
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      re*=factor;
      sns1d.put_data(i,re*3.);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          int maxg=grp;
	  if(nid>1000&&i+3<maxg)maxg=i+3;
	  for(int k=i;k<maxg;k++){
	    sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigs1org[j*grp+k]);
	  };
	};
      };
     };
    }else{
      for(int i=0;i<grp;i++){
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
	    totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	    sigsorg[j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,i);
	    real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigel).get_dat(i);
	    sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,-den*delta*micsigs);
	    sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,i,-den*delta*micsigs);
	  };
	};
        real re=0.;
	CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
	re+=scttmp[0]+leaktmp[0];
	re*=factor;
        sns1d.put_data(i,re*3.);
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,i,sigsorg[j]);
	  };
	};
      };
    };
    sens.PutSensitivity1D(nid,251,sns1d);
  };

  delete [] scttmp;
  delete [] leaktmp;

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

// For CMFD acceleration

void SNTSystem::DoAcceleration(int iter,real errs,real fiss)
{
  if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+0.5);
  //if(iter==50)DoCMFDAcceleration(fiss+0.5);
};


void SNTSystem::DoCMFDAcceleration(real delk)
{
  int ngrp=1;
  vector<int> bgrp(ngrp);
  bgrp[0]=grp-1;

  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();

  // Calculation for fine-mesh current
  vector< vector< vector< vector< vector<real> > > > > ModD(ngrp);  // [g][z+1][y+1][x+1][dim]
  for(int g=0;g<ngrp;g++){
    ModD[g].resize(zr+1);
    for(int i=0;i<zr+1;i++){
      ModD[g][i].resize(yr+1);
      for(int j=0;j<yr+1;j++){
        ModD[g][i][j].resize(xr+1);
        for(int k=0;k<xr+1;k++){
	  ModD[g][i][j][k].resize(Ndim,0.);
	};
      };
    };
  };

  // Calculation for coarse-group cross sections

  int xyzr=xr*yr*zr;
  PLOSSystem cm(Ndim,ngrp,xyzr);
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
		/*
		cout.setf(ios::scientific);
		cout.precision(5);
		cout<<id<<" "<<volflx.get_dat(0)/vol<<"\n";
		*/
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
      real dinp=cc_d[i]*tmp*cmfd_factor;
      if(dinp>50.)dinp=50.;
      // +++ For voided region +++
      //minp.GetMacxs().GetD().put_data(0,cc_d[i]*tmp*cmfd_factor);
      minp.GetMacxs().GetD().put_data(0,dinp);
      cm.AddMedium(minp);

      /*
    cout.setf(ios::showpoint);
    cout.precision(10);
    cout<<i<<" "<<fis<<"\n";
      */
  };
  //cout<<"\n\n";

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
  //if(Cylinder)ss="Cylinder";
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
    };
  };

  GeneralOption option;
  real cm_epsf=opt.GetEpsf()*0.1;
  real cm_epsk=opt.GetEpsk()*0.1;
  real cm_epss=opt.GetEpss()*0.1;
  real low_eps=1e-7;
  if(cm_epsf>low_eps)cm_epsf=low_eps;
  if(cm_epsk>low_eps)cm_epsk=low_eps;
  if(cm_epss>low_eps)cm_epss=low_eps;
  option.PutEpsf(cm_epsf);
  option.PutEpsk(cm_epsk);
  option.PutEpss(cm_epss);
  if(!opt.Forward())option.PutAdjointCal();
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
    cc_flx[i]=cm.GetMesh(i).GetFlux().get_dat(0)*tmp/cc_flx[i];
    /*
    cout.setf(ios::showpoint);
    cout.precision(10);
    cout<<i<<" "<<cc_flx[i]<<"\n";
    */
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
		int g2=0;
		real tmp=cc_flx[index];
		for(int g=0;g<grp;g++){
		  if(g>bgrp[g2]){
                    g2++;
		    tmp=cc_flx[index];
		  };
		  /*
    cout.setf(ios::showpoint);
    cout.precision(10);
    real tt=mesh[id].GetFlux().get_dat(g);
    cout<<ix<<" "<<(tmp-1)*tt<<"\n";
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
  //cout<<"\n\n";
};

void SNTSystem::SetZeroCurFF()
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();
  for(int i=0;i<zr+1;i++){
    for(int j=0;j<yr+1;j++){
      for(int k=0;k<xr+1;k++){
        for(int l=0;l<Ndim;l++){
    	  CurFF[0][i][j][k][l]=0.;
	};
      };
    };
  };
};

void SNTSystem::CalCoarseCur()
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();
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
  	        real sx=mesh[index].GetSurL(0);
	        real sy=mesh[index].GetSurL(1);
	        real sz=mesh[index].GetSurL(2);
  	        // X-direction
		if(x2==0){
                  //if(ix!=xedgel[iy]){
		  if(ix!=edge[0][iy][iz][0]){
		    CurFF[0][z1][y1][x1][0]+=(curxn[index]+curxp[index-1])*sx;
		  }else if(BC[0]!=1){
		    CurFF[0][z1][y1][x1][0]+=curxn[index]*sx;
		  };
		};
                //if(ix==xedger[iy]&&BC[1]!=1){
                if(ix==edge[0][iy][iz][1]&&BC[1]!=1){
		  CurFF[0][z1][y1][xr][0]+=curxp[index]*sx;
		};
	        // Y-direction
	        if(Ndim>1){
	  	  if(y2==0){
		    //if(iy!=yedgel[ix]){
		    if(iy!=edge[1][ix][iz][0]){
		      CurFF[0][z1][y1][x1][1]+=(curyn[index]+curyp[meshid[iz][iy-1][ix]])*sy;
		    }else if(BC[2]!=1){
		      CurFF[0][z1][y1][x1][1]+=curyn[index]*sy;
		    };
		  };
                  //if(iy==yedger[ix]&&BC[3]!=1){
                  if(iy==edge[1][ix][iz][1]&&BC[3]!=1){
		    CurFF[0][z1][yr][x1][1]+=curyp[index]*sy;
		  };
		};
	        // Z-direction
	        if(Ndim>2){
	  	  if(z2==0){
		    //if(iz!=0){
		    if(iz!=edge[2][ix][iy][0]){
		      CurFF[0][z1][y1][x1][2]+=(curzn[index]+curzp[meshid[iz-1][iy][ix]])*sz;
		    }else if(BC[4]!=1){
		      CurFF[0][z1][y1][x1][2]+=curzn[index]*sz;
		    };
		  };
                  //if(iz==mi.GetZF()-1&&BC[5]!=1){
                  if(iz==edge[2][ix][iy][1]&&BC[5]!=1){
		    CurFF[0][zr][y1][x1][2]+=curzp[index]*sz;
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
};

void SNTSystem::MemoryReductionForPerturbation()
{
  CurFF.clear();
  curxp.clear();
  curxn.clear();
  curyp.clear();
  curyn.clear();
  curzp.clear();
  curzn.clear();

  GetPLOSESystem().AllVectorClear();

  for(int i=0;i<TotM;i++){
    mesh[i].SrcDataClear();
  };
};

void SNTSystem::AddFlux(SNTSystem &sec)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)+sec.GetMesh(i).GetFlux(l);
      mesh[i].GetFlux(l).copy(tmp);
    };
    for(int g=0;g<grp;g++){
      GroupData1D tmp=aflux[i][g]+sec.GetAFlux(i,g);
      aflux[i][g].copy(tmp);
    };
  };
};

void SNTSystem::NegFlux(SNTSystem &sec)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)-sec.GetMesh(i).GetFlux(l);
      mesh[i].GetFlux(l).copy(tmp);
    };
    for(int g=0;g<grp;g++){
      GroupData1D tmp=aflux[i][g]-sec.GetAFlux(i,g);
      aflux[i][g].copy(tmp);
    };
  };
};

void SNTSystem::NegFlux(SNTSystem &sec,real fact)
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      GroupData1D tmp=mesh[i].GetFlux(l)-sec.GetMesh(i).GetFlux(l)*fact;
      mesh[i].GetFlux(l).copy(tmp);
    };
    for(int g=0;g<grp;g++){
      GroupData1D tmp=aflux[i][g]-sec.GetAFlux(i,g)*fact;
      aflux[i][g].copy(tmp);
    };
  };
};

