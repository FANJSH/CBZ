#include <cstdlib>
#include "SNR_system.h"

using namespace std;

SNRSystem::SNRSystem(int n,int g,int i):GeneralSystem(n,g,i)
{
  //cout<<"# !! Caution !!\n";
  //cout<<"# DSA Acceleration is applied as default.\n";

  psys.Init(n,g,i);
  pl = -1;
  nquad = 1;
  quad.resize(nquad);
  quadid.resize(grp,0);
  for(int i=0;i<grp;i++){quadid[i]=0;};
  transport=true;
  etransport=false;
  name="SNR";
  cmfdimp=false; // CMFD is not implemented
  dsa=true;
  step_d=false;
  surface_aflx_store=false;
  leakage.put_imax(grp);
  // for leakage sensitivity calculation
  leak_sens=false;
  leak_mesh=26;
  //top_grp=0;
  top_grp=20;
  //btm_grp=g-1;
  btm_grp=47;
  // +++ for white boundary condition treatment
  inflx_white.resize(grp,0.);
};

SNRSystem::~SNRSystem()
{
};

void SNRSystem::SetQuadratureNum(int num)
{
  nquad=num;
  quad.resize(nquad);
};

void SNRSystem::SetQuadratureID(int *gbnd, int *id)
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

void SNRSystem::SetQuadrature(SNRQuadrature *qinp,int id)
{
  if(id<0||id>=nquad){
    cout<<"Error in SetQuadratureID.\n";
    cout<<"Chack your ID number.\n";
    exit(0);
  };
  quad[id]=qinp;
};

void SNRSystem::SetArray()
{
  if(pl==-1){
    cout<<"Do 'PutPL'.\n";
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
      aflux[i].resize(grp);
      for(int j=0;j<grp;j++){
	aflux[i][j].put_imax(quad[quadid[j]]->GetSN());
	aflux[i][j].set_zero();
      };
    };
  };

  if(pl>0){
    if(!etransport&&transport){
    for(int i=0;i<nmed;i++){
      for(int j=1;j<pl+1;j++){
        for(int k=0;k<grp;k++){
	  real tmp=med[i].GetDataSigs(j,k,k);
	  //int pltt=j;
	  //if(med[i].GetPLT()<j)pltt=med[i].GetPLT();
	  //real tmp2=med[i].GetDataSigt(0,k)-med[i].GetDataSigt(pltt,k);
	  real tmp2=med[i].GetDataSigt(0,k)-med[i].GetDataSigt(1,k);
	  real tmp3=tmp+(j*2+1.)*tmp2;
	  real tmp4=med[i].GetSigs(0).get_dat(k,k);

	  if(j==1&&fabs(tmp3*0.33333333)>tmp4){
	    cout<<"#   ... warning... mu (sp1/sp0) is out of range : "<<tmp3*0.33333333/tmp4<<"\n";
	    cout<<"#         medium "<<i<<"  group "<<k<<"\n";
	    cout<<"#         this value IS USED in the present calculation.\n";
            //med[i].GetSigs(j).put_data(k,k,tmp3);	  	    
	  }else{
            //med[i].GetSigs(j).put_data(k,k,tmp3);
	  };
          med[i].GetSigs(j).put_data(k,k,tmp3);	  
	};
      };
    };
    };
    if(etransport){
      // (extended transport approximation)
      for(int i=0;i<nmed;i++){
        for(int k=0;k<grp;k++){
	  real a1=0.;
	  real a2=0.;
	  for(int j=pl+1;j<=5;j++){
	    a1+=med[i].GetMacxs().GetData2d(sigs,j).get_sumx().get_dat(k)/(2*j+1)/(2*j+1);
	    a2+=1./(2*j+1);
	  };
	  //real totsk=a1/a2;
	  real totsk=med[i].GetMacxs().GetData2d(sigs,pl+1).get_sumx().get_dat(k)/(2.*pl+3);

	  real total=med[i].GetMacxs().GetData1d(sigt,1).get_dat(k)-totsk;

 	  real totorg=med[i].GetMacxs().GetData1d(sigt).get_dat(k);
 	  real totorg1=med[i].GetMacxs().GetData1d(sigt,1).get_dat(k);
	  med[i].GetMacxs().GetData1d(sigt).put_data(k,total);
	  for(int j=0;j<=pl;j++){
	    real tmp=total-totorg;
	    if(j!=0)tmp=total-totorg1;
	    med[i].GetMacxs().GetData2d(sigs,j).add_data(k,k,tmp*(2*j+1));
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
};

void SNRSystem::PutWriteFlux(int med)
{
  int ind=0;
  for(int j=0;j<mi.GetXF();j++){
    int tmp=mi.GetFMat(ind);
    if(tmp==med||med==-1)wrtflx[ind]=true;
    ind++;
  };
};

void SNRSystem::PutCartMeshInfo(CartMeshInfo cm, string geom)
{
  if(geom!="Sphere"&&geom!="Cartesian"){
    cout<<"Error in SNRSystem::PutCartMeshInfo.\n";
    exit(0);
  };
  sphere=false;
  if(geom=="Sphere")sphere=true;

  mi=cm;

  TotM=mi.GetXF();
  mesh.resize(TotM);
  wrtflx.resize(TotM,false);
  aflux.resize(TotM);

  meshid.resize(1);
  meshid[0].resize(1);
  meshid[0][0].resize(TotM);
  for(int i=0;i<TotM;i++){
    meshid[0][0][i]=i;
  };

  real r1=0.;
  real r2;
  real l[3];
  for(int i=0;i<TotM;i++){
    int tm=mi.GetFMat(i);
    if(tm>=nmed){
      cout<<"Error in SNRSystem::PutCartMeshInfo.\n";
      cout<<"You requested not-existing medium ID.\n";
      cout<<"Please check Medium ID.\n";
      exit(0);
    };
    mesh[i].PutMedium(&med[mi.GetFMat(i)]);
    if(sphere){
      r2=r1+mi.GetFMeshL(0,i);
      mesh[i].PutDimSphere(mi.GetFMeshL(0,i),r1,r2);
      r1=r2;
    }else{
      l[0]=mi.GetFMeshL(0,i);
      mesh[i].PutDim(1,l);
    };
  };

  for(int i=0;i<6;i++){
    BC[i]=-1;
    if(mi.GetBC(i)==0)BC[i]=0; // Zeroflux
    if(mi.GetBC(i)==1)BC[i]=1; // Reflective
    if(mi.GetBC(i)==2)BC[i]=2; // Vacuum
    if(mi.GetBC(i)==3)BC[i]=3; // Periodic
    if(BC[i]<=0||BC[i]>3){
      cout<<"Incorrect B. C. in SNX_system.\n";
      exit(0);
    };
  };

  if(sphere&&BC[1]==1){
    cout<<"#\n";
    cout<<"# White boundary is assigned at the right edge.\n";
    cout<<"# DSA acceleration is NOT adopted.\n";
  };

  if(dsa){
    if(!print)psys.NoPrint();
    psys.PutCartMeshInfo(cm,geom);
  };
};

void SNRSystem::CalFluxInSource(real *SSrc,real *flxmom,real *selfsc,int qid,int sn)
{
  vector<real> moment(sn*plnum);
  int index=0;
  for(int l=0;l<=pl;l++){
    for(int k=0;k<sn;k++){
      moment[index++]=quad[qid]->GetMoment(l,k);
    };
  };

  real *tmp=new real[sn];
  real ssc;
  int ind=0;
  int ind2=0;
  int ind3=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<sn;j++){tmp[j]=0.;};
    int index=0;
    for(int l=0;l<=pl;l++){
      real sigs=selfsc[ind3]; // self-scattering xs * volume * 1/4PI
      ind3++; 	 
      for(int m=-0;m<=0;m++){ // for 1D
        ssc=sigs*flxmom[ind2];
	ind2++;
        for(int k=0;k<sn;k++){
	  //tmp[k]+=ssc*quad[qid]->GetMoment(index,k);  
	  tmp[k]+=ssc*moment[index++];
        };
	//index++;
      };
    };
    for(int j=0;j<sn;j++){
      //real tt=tmp[j];
      //if(tt<0.)tt=0;
      //SSrc[ind]=tt;
      SSrc[ind]=tmp[j];  // original
      ind++;
    };

  };
  delete [] tmp;
};

void SNRSystem::CalFluxInSourceSP(real *SSrc,real *flxmom,real *selfsc,int qid,int sn)
{
  real mu=-1.;
  //if(!opt.Forward())mu=1.;

  vector<real> moment(plnum);
  for(int i=0;i<plnum;i++){
    moment[i]=quad[qid]->GetMoment(i,mu);
  };
  int ind2=0;
  int ind3=0;
  for(int i=0;i<TotM;i++){
    real tmp=0.;
    for(int l=0;l<=pl;l++){
      real sigs=selfsc[ind3];
      ind3++;
      for(int m=-0;m<=0;m++){ // for 1D
        //tmp+=sigs*flxmom[ind2]*quad[qid]->GetMoment(l,mu);  
        tmp+=sigs*flxmom[ind2]*moment[l];
	ind2++;
      };
    };
    //real tt=tmp;
    //if(tt<0.)tt=0.;
    //SSrc[i]=tt;
    SSrc[i]=tmp;
  };
};

real SNRSystem::CalFluxGeneral(int ginp,real cin,int iter)
{
  if(sphere){
    if(cin>1e-5)cin=1e-5;
    return CalFlux(ginp,cin);
  }else{
    return CalFluxSlab(ginp,cin);
  };
};

real SNRSystem::CalFlux(int g,real epsif)
{
  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src    = new real[TotM*sn]; // External source
  real *SSrc   = new real[TotM*sn]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM];
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  // Put sigma_t and sigma_s to temporal array
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    sigtvol[i]=mesh[i].GetMed()->GetData1D(g,sigt)*vol;
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      sigsself[ind]=mesh[i].GetMed()->GetDataSigs(l,g,g)*volpi4;
      ind++;
    };
  };

  // Initial Flux guess
  ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };

  // Put fission & slowing down neutron source
  real *srcin=new real[sn];
  ind=0;
  for(int j=0;j<TotM;j++){
    real tmp=mesh[j].GetSrcin(0);
    for(int i=0;i<sn;i++){srcin[i]=tmp;};
    for(int l=1;l<plnum;l++){
      real tmp=mesh[j].GetSrcin(l);
      for(int i=0;i<sn;i++){
        srcin[i]+=tmp*quad[qid]->GetMoment(l,i);
      };
    };
    for(int i=0;i<sn;i++){
      Src[ind]=srcin[i];
      //if(g==43||g==44)cout<<j<<" "<<i<<" "<<Src[ind]<<"\n";
      ind++;
    };
  };
  //if(g==44)exit(0);
  delete [] srcin;

  // Put fission & slowing down source for mu=-1.
  real *SrcSP  =new real[TotM];
  real mu=-1.;
  //if(!opt.Forward())mu=1.;
  for(int i=0;i<TotM;i++){
    real tmp=0.;
    for(int l=0;l<=pl;l++){
      tmp+=mesh[i].GetSrcin(l)*quad[qid]->GetMoment(l,mu);
    };
    SrcSP[i]=tmp;
  };

  real *flangle=new real[TotM];
  real *SSrcSP =new real[TotM];

  // (for white boundary treatment)
  real flx_wht,flx_wht_denom;

  real *inifl=new real[sn/2];
  int itmax=9999;
  bool inconv=true;
  for(int it=0;it<itmax;it++){

    // for leakage calculation
    //leakage.put_data(g,0.);

    CalFluxInSourceSP(SSrcSP,flxmom,sigsself,qid,sn);
    CalFluxInSource(SSrc,flxmom,sigsself,qid,sn);

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    // ******************************
    // Starting direction method
    real afl=0.;
    if(BC[1]==1)afl=inflx_white[g]; // (white boundary)
    for(int m=TotM-1;m>=0;m--){
      real ss =SrcSP[m]+SSrcSP[m];
      // *** for sns
      if(leak_sens&&m==leak_mesh&&!opt.Forward()&&g>=top_grp&&g<=btm_grp){ss=1.;};

      real len=mesh[m].GetLen(0);
      real st=mesh[m].GetMed()->GetDataSigt(0,g);
      real ss_len=ss/mesh[m].GetVolume()*len;
      real bfl=(2.*afl+ss_len)/(2.+st*len);
      real out=2.*bfl-afl;
      if(out<0.){
	out=0.;
	bfl=(ss_len+afl)/(st*len);
      };
      flangle[m]=bfl;
      afl=out;
    };

    flx_wht=0.;
    flx_wht_denom=0.;
    for(int iss=0;iss<sn;iss++){
      afl=0.;
      if(iss<sn/2&&BC[1]==1){afl=inflx_white[g];}; // (white boundary)
      if(iss>=sn/2&&BC[0]==1){afl=inifl[sn-iss-1];};
      int is=iss;
      //if(!opt.Forward())is=sn-iss-1;
      real mu =quad[qid]->GetMu(is);
      real muu=fabs(mu);
      real om =quad[qid]->GetOmega(is);
      real mneg=quad[qid]->GetMC(is); // alpha(m-1/2)
      real mpos=mneg+mu*om; // alpha(m+1/2)
      for(int mm=0;mm<TotM;mm++){
	int m=mm;
	if(iss<sn/2){m=TotM-mm-1;};
	int index2=m*sn+is;
	real ss =Src[index2]+SSrc[index2];
	// *** for sns
	if(leak_sens&&m==leak_mesh&&!opt.Forward()&&mu<0.&&g>=top_grp&&g<=btm_grp){ss-=mu;};
	real bfl;
	real sl=mesh[m].GetSurL(0);
	real sr=mesh[m].GetSurR(0);
	real sarea;
	if(is<sn/2){sarea=sl;}else{sarea=sr;};
        real om_inv=1./om;
	real lhs=2.*muu*sarea-2.*om_inv*(sr-sl)*mpos+sigtvol[m];
	real rhs=muu*(sr+sl)*afl-(sr-sl)*(mneg+mpos)*flangle[m]*om_inv+ss;
	bfl=rhs/lhs;
	real outx=2.*bfl-afl;
	real outm=2.*bfl-flangle[m];
	real t1=1.;
	real t2=1.;
	while(outx<0.||outm<0.){
	  if(outx<0.)t1=0.;
	  if(outm<0.)t2=0.;
	  lhs=2.*muu*sarea*t1-2.*om_inv*(sr-sl)*mpos*t2+sigtvol[m];
	  rhs=muu*(sr+sl-sarea*(1.-t1))*afl-(sr-sl)*(mneg+mpos*t2)*flangle[m]*om_inv+ss;
	  bfl=rhs/lhs;
	  outx=(2.*bfl-afl)*t1;
	  outm=(2.*bfl-flangle[m])*t2;
	};
	anflx[index2]=bfl;
	afl=outx;
	flangle[m]=outm;
	// for leakage calculation
        //if(mu>0.&&mm==TotM-1){
	//  leakage.add_data(g,outx*mu*om*PI2*sr);
	//};
      };
      if(iss<sn/2){
        inifl[iss]=afl;
      }else if(BC[1]==1){
        flx_wht+=afl*om*mu;
        flx_wht_denom+=om*mu;
      };
    };

    // Renew flux moment
    real *finp=new real[plnum];
    int ind=0;
    int ind2=0;
    for(int i=0;i<TotM;i++){
      for(int j=0;j<plnum;j++){finp[j]=0.;};
      for(int j=0;j<sn;j++){
	real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
	ind++;
        finp[0]+=tmp;
	for(int l=1;l<plnum;l++){
	  finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	};
      };
      for(int l=0;l<plnum;l++){
	flxmom[ind2]=finp[l];
	ind2++;
      };
    };
    delete [] finp;

    if(dsa)AccelerationByDSA(g,flxmom,flxd,sigsself);

    real errmax=0.;
    ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      //if(i==0&&g==69)cout<<it<<" : "<<flxmom[ind]<<" "<<err<<"\n";
      //if(i==0&&g>=44)cout<<g<<" "<<i<<" "<<it<<" : "<<flxmom[ind]<<" "<<err<<"\n";
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    if(errmax<epsif){
      //cout<<g<<" "<<it<<" ("<<errmax<<")\n";
      break;
    };
    if(it==itmax-1){
      cout<<"   ... Not converged in group "<<g<<"\n";
      inconv=false;
    };
  };

  real errf=0.;
  // Put new flux into GeneralMesh
  ind=0;
  for(int i=0;i<TotM;i++){
    real err=fabs(flxmom[ind]/mesh[i].GetFlux(0).get_dat(g)-1.0);
    for(int j=0;j<plnum;j++){
      mesh[i].GetFlux(j).put_data(g,flxmom[ind]);
      ind++;
    };
    if(err>errf)errf=err;
  };

  // Write angular flux for perturbation calculation **
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]);
      ind++;
    };
  };

  // (for white boundary treatment)
  if(BC[1]==1)inflx_white[g]=flx_wht/flx_wht_denom;

  delete [] flangle;
  delete [] SrcSP;
  delete [] SSrcSP;
  delete [] sigtvol;
  delete [] sigsself;
  delete [] inifl;
  delete [] flxmom;
  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxd;

  if(!inconv)errf*=-1;
  return errf;
};

//real SNRSystem::CalFluxSlabDF(int g,real epsif,vector<real> inflx_df,vector< vector<real> >outflx_df_l,vector< vector<real> > outflx_df_r)
real SNRSystem::CalFluxSlabDF(int g,real epsif,vector<real> inflx_df, vector<real> outflx_df_l, vector<real> outflx_df_r)
{
  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src    = new real[TotM*sn]; // External source
  real *SSrc   = new real[TotM*sn]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM];
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  // Put sigma_t and sigma_s to temporal array
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    sigtvol[i]=mesh[i].GetMed()->GetData1D(g,sigt)*vol;
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      sigsself[ind]=mesh[i].GetMed()->GetDataSigs(l,g,g)*volpi4;
      ind++;
    };
  };

  // Initial Flux guess
  ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };

  // Put fission & slowing down neutron source
  real *srcin=new real[sn];
  ind=0;
  for(int j=0;j<TotM;j++){
    real tmp=mesh[j].GetSrcin(0);
    for(int i=0;i<sn;i++){srcin[i]=tmp;};
    for(int l=1;l<plnum;l++){
      real tmp=mesh[j].GetSrcin(l);
      for(int i=0;i<sn;i++){
        srcin[i]+=tmp*quad[qid]->GetMoment(l,i);
      };
    };
    for(int i=0;i<sn;i++){
      Src[ind]=srcin[i];
      ind++;
    };
  };
  delete [] srcin;

  vector<real> inifl(sn,0.);

  int itmax=20000;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    CalFluxInSource(SSrc,flxmom,sigsself,qid,sn);

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    for(int iss=0;iss<sn;iss++){
      int is=iss;
      real afl=0.;
      int xrsn=quad[qid]->GetXref(is);
      if(iss>=sn/2){ // ====>
	if(BC[0]==1){afl=inifl[xrsn];}
	else if(BC[0]==3){afl=inifl[is];};
      }else{ // <====
	if(BC[1]==1){afl=inifl[xrsn];}
	else if(BC[1]==3){afl=inifl[is];};
      };
      real mu =fabs(quad[qid]->GetMu(is));
      for(int mm=0;mm<TotM;mm++){
	int m=mm;
	if(iss<sn/2){m=TotM-mm-1;};
	// +++++ incoming flx DF
	if(iss<sn/2){
          if((m+1)%70==0){ // From right side
	    afl*=inflx_df[(m+1)/70-1];
	  };
	}else{ 
	  if(m%70==0){ // From left side
	    afl*=inflx_df[m/70];
	  };
	};
	// +++++++++++++++++++++
	int index2=m*sn+is;
	real ss =Src[index2]+SSrc[index2];
        real rhs,bfl,cfl;
	rhs=(mu-0.5*sigtvol[m])*afl+ss;
	bfl=rhs/(mu+0.5*sigtvol[m]);
	cfl=(afl+bfl)*0.5;
	if(bfl<0.){
	  bfl=0.;
	  cfl=(ss+mu*afl)/sigtvol[m];
	};
	anflx[index2]=cfl;
	afl=bfl;

	// +++++ outgoing flx DF +++
	if(iss<sn/2){
	  if(m%70==0){
	    afl*=outflx_df_l[m/70];
	  };
	}else{
	  if((m+1)%70==0){
	    afl*=outflx_df_r[(m+1)/70-1];
	  };
	}; 
	// ++++++++++++++++++++++++++
	if(surface_aflx_store)aflux[m][g].put_data(is,afl);// temporary for current Calc.
      };
      inifl[is]=afl;
    };

    // Renew flux moment
    real *finp=new real[plnum];
    int ind=0;
    int ind2=0;
    for(int i=0;i<TotM;i++){
      for(int j=0;j<plnum;j++){finp[j]=0.;};
      for(int j=0;j<sn;j++){
	real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
	ind++;
        finp[0]+=tmp;
	for(int l=1;l<plnum;l++){
	  finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	};
      };
      for(int l=0;l<plnum;l++){
	flxmom[ind2]=finp[l];
	ind2++;
      };
    };
    delete [] finp;

    real errmax=0.;
    ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    if(errmax<epsif){
      //cout<<" Inner iteration in group "<<g<<" : "<<it<<" ("<<errmax<<")\n";
      break;
    };
    if(it==itmax-1){
      //cout<<"   ... Not converged in group "<<g<<" ("<<errmax<<"/"<<epsif<<")\n";
      InnerConvergence=false;
    };
  };

  real errf=0.;
  // Put new flux into GeneralMesh
  ind=0;
  for(int i=0;i<TotM;i++){
    real err=fabs(flxmom[ind]/mesh[i].GetFlux(0).get_dat(g)-1.0);
    for(int j=0;j<plnum;j++){
      mesh[i].GetFlux(j).put_data(g,flxmom[ind]);
      ind++;
    };
    if(err>errf)errf=err;
  };
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  if(!surface_aflx_store){
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]); 
      ind++;
    };
  };
  };

  delete [] sigtvol;
  delete [] sigsself;
  delete [] flxmom;
  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxd;

  return errf;
};

real SNRSystem::CalFluxSlabDF(int g,real epsif, vector<real> outflx_df_l, vector<real> outflx_df_r)
{
  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src    = new real[TotM*sn]; // External source
  real *SSrc   = new real[TotM*sn]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM];
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  // Put sigma_t and sigma_s to temporal array
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    sigtvol[i]=mesh[i].GetMed()->GetData1D(g,sigt)*vol;
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      sigsself[ind]=mesh[i].GetMed()->GetDataSigs(l,g,g)*volpi4;
      ind++;
    };
  };

  // Initial Flux guess
  ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };

  // Put fission & slowing down neutron source
  real *srcin=new real[sn];
  ind=0;
  for(int j=0;j<TotM;j++){
    real tmp=mesh[j].GetSrcin(0);
    for(int i=0;i<sn;i++){srcin[i]=tmp;};
    for(int l=1;l<plnum;l++){
      real tmp=mesh[j].GetSrcin(l);
      for(int i=0;i<sn;i++){
        srcin[i]+=tmp*quad[qid]->GetMoment(l,i);
      };
    };
    for(int i=0;i<sn;i++){
      Src[ind]=srcin[i];
      ind++;
    };
  };
  delete [] srcin;

  vector<real> inifl(sn,0.);

  int itmax=20000;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    CalFluxInSource(SSrc,flxmom,sigsself,qid,sn);

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    for(int iss=0;iss<sn;iss++){
      int is=iss;
      real afl=0.;
      int xrsn=quad[qid]->GetXref(is);
      if(iss>=sn/2){ // ====>
	if(BC[0]==1){afl=inifl[xrsn];}
	else if(BC[0]==3){afl=inifl[is];};
      }else{ // <====
	if(BC[1]==1){afl=inifl[xrsn];}
	else if(BC[1]==3){afl=inifl[is];};
      };
      real mu =fabs(quad[qid]->GetMu(is));
      for(int mm=0;mm<TotM;mm++){
	int m=mm;
	if(iss<sn/2){m=TotM-mm-1;};
	int index2=m*sn+is;
	real ss =Src[index2]+SSrc[index2];
        real rhs,bfl,cfl;
	rhs=(mu-0.5*sigtvol[m])*afl+ss;
	bfl=rhs/(mu+0.5*sigtvol[m]);
	cfl=(afl+bfl)*0.5;
	if(bfl<0.){
	  bfl=0.;
	  cfl=(ss+mu*afl)/sigtvol[m];
	};
	anflx[index2]=cfl;
	afl=bfl;

	// +++++ outgoing flx DF +++
	if(iss<sn/2){
	  if(m%70==0){
	    afl*=outflx_df_l[m/70];
	  };
	}else{
	  if((m+1)%70==0){
	    afl*=outflx_df_r[(m+1)/70-1];
	  };
	}; 
	// ++++++++++++++++++++++++++
	if(surface_aflx_store)aflux[m][g].put_data(is,afl);// temporary for current Calc.
      };
      inifl[is]=afl;
    };

    // Renew flux moment
    real *finp=new real[plnum];
    int ind=0;
    int ind2=0;
    for(int i=0;i<TotM;i++){
      for(int j=0;j<plnum;j++){finp[j]=0.;};
      for(int j=0;j<sn;j++){
	real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
	ind++;
        finp[0]+=tmp;
	for(int l=1;l<plnum;l++){
	  finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	};
      };
      for(int l=0;l<plnum;l++){
	flxmom[ind2]=finp[l];
	ind2++;
      };
    };
    delete [] finp;

    real errmax=0.;
    ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    if(errmax<epsif){
      //cout<<" Inner iteration in group "<<g<<" : "<<it<<" ("<<errmax<<")\n";
      break;
    };
    if(it==itmax-1){
      //cout<<"   ... Not converged in group "<<g<<" ("<<errmax<<"/"<<epsif<<")\n";
      InnerConvergence=false;
    };
  };

  real errf=0.;
  // Put new flux into GeneralMesh
  ind=0;
  for(int i=0;i<TotM;i++){
    real err=fabs(flxmom[ind]/mesh[i].GetFlux(0).get_dat(g)-1.0);
    for(int j=0;j<plnum;j++){
      mesh[i].GetFlux(j).put_data(g,flxmom[ind]);
      ind++;
    };
    if(err>errf)errf=err;
  };
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  if(!surface_aflx_store){
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]); 
      ind++;
    };
  };
  };

  delete [] sigtvol;
  delete [] sigsself;
  delete [] flxmom;
  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxd;

  return errf;
};

real SNRSystem::CalFluxSlab(int g,real epsif,bool angular_dependent_source)
{
  if(epsif>1e-6)epsif=1e-6;

  int qid = quadid[g];
  int sn  = quad[qid]->GetSN();

  real *Src    = new real[TotM*sn]; // External source
  real *SSrc   = new real[TotM*sn]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  =new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM];
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  // Put sigma_t and sigma_s to temporal array
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    sigtvol[i]=mesh[i].GetMed()->GetData1D(g,sigt)*vol;
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      sigsself[ind]=mesh[i].GetMed()->GetDataSigs(l,g,g)*volpi4;
      //if(angular_dependent_source)sigsself[ind]=0.;
      //sigsself[ind]=0.;
      ind++;
    };
  };

  // Initial Flux guess
  ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };

  // Put fission & slowing down neutron source
  real *srcin=new real[sn];
  ind=0;
  for(int j=0;j<TotM;j++){
    real tmp=mesh[j].GetSrcin(0);
    for(int i=0;i<sn;i++){srcin[i]=tmp;};
    for(int l=1;l<plnum;l++){
      real tmp=mesh[j].GetSrcin(l);
      for(int i=0;i<sn;i++){
        srcin[i]+=tmp*quad[qid]->GetMoment(l,i);
      };
    };
    for(int i=0;i<sn;i++){
      Src[ind]=srcin[i];
      if(angular_dependent_source){
	Src[ind]+=aflux[j][g].get_dat(i);
      };
      ind++;
    };
  };
  delete [] srcin;

  vector<real> inifl(sn,0.);

  int itmax=20000;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    CalFluxInSource(SSrc,flxmom,sigsself,qid,sn);

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    for(int iss=0;iss<sn;iss++){
      int is=iss;
      //if(!opt.Forward())is=sn-iss-1;
      real afl=0.;
      int xrsn=quad[qid]->GetXref(is);
      if(iss>=sn/2){ // ====>
	if(BC[0]==1){afl=inifl[xrsn];}
	else if(BC[0]==3){afl=inifl[is];};
      }else{ // <====
	if(BC[1]==1){afl=inifl[xrsn];}
	else if(BC[1]==3){afl=inifl[is];};
      };
      //if(quad[qid]->GetMu(is)>0.){
      //if(BC[0]==1)afl=inifl[xrsn];
      //if(BC[0]==3)afl=inifl[is];
      //}else{
      //if(BC[1]==1)afl=inifl[xrsn];
      //if(BC[1]==3)afl=inifl[is];
      //};
      real mu =fabs(quad[qid]->GetMu(is));
      for(int mm=0;mm<TotM;mm++){
	int m=mm;
	if(iss<sn/2){m=TotM-mm-1;};
	int index2=m*sn+is;
	real ss =Src[index2]+SSrc[index2];
        real rhs,bfl,cfl;
	// Step differencing
	if(step_d){
	rhs=(ss+fabs(mu)*afl);
        bfl=rhs/(fabs(mu)+sigtvol[m]);
	cfl=bfl; 
	}else{
	// Diamond differencing

	  // (For direction-dependent)
	  /*
	  if(g==0&&quad[qid]->GetMu(is)>0.){
	    sigtvol[m]=mesh[m].GetMed()->GetMacxs().GetData1d(dr).get_dat(g)*mesh[m].GetVolume();
	  }else if(g==0&&quad[qid]->GetMu(is)<0.){
	    sigtvol[m]=mesh[m].GetMed()->GetMacxs().GetData1d(dz).get_dat(g)*mesh[m].GetVolume();
	  };
	  */

	rhs=(mu-0.5*sigtvol[m])*afl+ss;
	bfl=rhs/(mu+0.5*sigtvol[m]);
	cfl=(afl+bfl)*0.5;
	//if(angular_dependent_source&&g==0&&m==0&&is==0&&ss>0)cout<<afl<<" -> "<<bfl<<" "<<it<<"\n";
	if(bfl<0.){
	  bfl=0.;
	  cfl=(ss+mu*afl)/sigtvol[m];
	};
	};
	anflx[index2]=cfl;
	afl=bfl;

	// +++ For outgoing current discontinuity factor
	/*
	if(iss<sn/2){
	  if(m%7==0){
	    //if(g==0)afl*=0.74097;
	    if(g==0)afl*=1.+(0.74097-1.)*1;
            if(g==1)afl*=1.7513;
	    //if(g==0)afl*=0.71651;
            //if(g==1)afl*=1.8718;
	  };
	}else{
	  if((m+1)%7==0){
	    //if(g==0)afl*=1.1549;
	    if(g==0)afl*=1.+(1.1549-1.)*1;
	    if(g==1)afl*=0.095664;
	    //if(g==0)afl*=1.1789;
	    //if(g==1)afl*=0.088786;
	  };
	};
	*/
	// ++++++++++++++++++++++++++

	if(surface_aflx_store)aflux[m][g].put_data(is,afl);// temporary for current Calc.

      };
      inifl[is]=afl;
    };

    // Renew flux moment
    real *finp=new real[plnum];
    int ind=0;
    int ind2=0;
    for(int i=0;i<TotM;i++){
      for(int j=0;j<plnum;j++){finp[j]=0.;};
      for(int j=0;j<sn;j++){
	real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
	ind++;
        finp[0]+=tmp;
	for(int l=1;l<plnum;l++){
	  finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	};
      };
      for(int l=0;l<plnum;l++){
	flxmom[ind2]=finp[l];
	ind2++;
      };
    };
    delete [] finp;

    real errmax=0.;
    ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    if(errmax<epsif){
      //cout<<" Inner iteration in group "<<g<<" : "<<it<<" ("<<errmax<<")\n";
      break;
    };
    if(it==itmax-1){
      //cout<<"   ... Not converged in group "<<g<<" ("<<errmax<<"/"<<epsif<<")\n";
      InnerConvergence=false;
    };
  };

  real errf=0.;
  // Put new flux into GeneralMesh
  ind=0;
  for(int i=0;i<TotM;i++){
    real err=fabs(flxmom[ind]/mesh[i].GetFlux(0).get_dat(g)-1.0);
    for(int j=0;j<plnum;j++){
      mesh[i].GetFlux(j).put_data(g,flxmom[ind]);
      ind++;
    };
    if(err>errf)errf=err;
  };
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  if(!surface_aflx_store){
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]); 
      ind++;
    };
  };
  };

  delete [] sigtvol;
  delete [] sigsself;
  delete [] flxmom;
  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxd;

  return errf;
};

void SNRSystem::AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself)
{
  //real *fltt=new real[TotM];
  // *** Source calculation
  int ind=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    real tmp=(flxmom[ind]-flxd[i])*sigsself[ind2]; 
    ind2+=pl+1;
    ind+=plnum;
    //fltt[i]=tmp;   
    psys.GetMesh(i).PutSrcin(tmp);
  };

  psys.CalFluxAutoConv(g,0.01);
  ind=0;
  for(int i=0;i<TotM;i++){
    real tmp1=flxmom[ind];
    real tmp2=psys.GetMesh(i).GetFlux(0).get_dat(g);
    real tmp3=tmp1+tmp2;
    if(tmp3>0.){
      flxmom[ind]=tmp3;
    };
    ind+=plnum;
  };
};


void SNRSystem::SetInitialFlux()
{
  if(dsa){
    GeneralOption optdif;
    optdif.PutEpsf(5e-2);
    optdif.PutEpsk(5e-4);
    optdif.PutEpss(5e-3);
    optdif.PutOutitermax(30);
    //optdif.PutOutitermax(3);
    if(!opt.Forward())optdif.PutAdjointCal();
    psys.PutGeneralOption(optdif);
    psys.CalCoef();
    psys.CalIgen("extrapolation");
    SetArray();
    for(int i=0;i<TotM;i++){
      mesh[i].GetFlux().copy(psys.GetMesh(i).GetFlux());
    };
  }else{
    SetArray();
    SetInitialFlatFlux();
  };
};

void SNRSystem::AddMedium(Medium inp)
{
  GeneralSystem::AddMedium(inp);
  Medium inp2=inp;
  for(int g=0;g<grp;g++){
    real org=inp2.GetMacxs().GetData1d(sigtr).get_dat(g);
    inp2.GetMacxs().GetData1d(d).put_data(g,0.281/org);
  };
  psys.AddMedium(inp2); // specific for snt
  //psys.AddMedium(inp); // specific for snt
}


void SNRSystem::ClearFlux()
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      mesh[i].GetFlux(l).set_zero();
    };
    for(int g=0;g<grp;g++){
      aflux[i][g].set_zero();
    };
  };
};

void SNRSystem::AddFlux(SNRSystem &sec)
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

void SNRSystem::NegFlux(SNRSystem &sec)
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

void SNRSystem::NegFlux(SNRSystem &sec,real fact)
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

void SNRSystem::FluxAngleReverse()
{
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      if(l%2==1){
        mesh[i].GetFlux(l)=mesh[i].GetFlux(l)*-1.;
      };
    };
  };
  for(int g=0;g<grp;g++){
    int qid = quadid[g];
    int sn  = quad[qid]->GetSN();
    for(int i=0;i<TotM;i++){
      if(wrtflx[i]){
        vector<real> tmpsn(sn);
        for(int j=0;j<sn;j++){
          tmpsn[j]=aflux[i][g].get_dat(j);
	};
        for(int j=0;j<sn;j++){
          int jj=sn-j-1;
	  aflux[i][g].put_data(j,tmpsn[jj]);
	};
      };
    };
  };
};

// **********************************************
// * For perturbation calculation               *
// **********************************************

real SNRSystem::CalReactivity(SNRSystem *sec,real kunp,real kp,int g)
{
  CheckSameMesh(sec);
  real ip=CalPerturbDenominator(sec);
  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };
  real *scttmp=new real[plnum];
  real *leaktmp=new real[plnum];
  real ret=0.;
  ret+=CalPerturbYieldTerm(sec,flag,g,kp);
  ret+=CalPerturbAbsorptionTerm(sec,flag,g);
  CalPerturbLeakScat(sec,flag,g,scttmp,leaktmp);
  for(int i=0;i<plnum;i++){
    ret+=scttmp[i]+leaktmp[i];
  };
  ret/=ip;

  delete [] flag;
  delete [] scttmp;
  delete [] leaktmp;

  return ret;
};

real SNRSystem::CalReactivity(SNRSystem *sec,real kunp,real kp,bool pr)
{
  vector<real> tmp;
  tmp.resize(grp,0.);
  return CalReactivity(sec,kunp,kp,tmp,pr);
};


real SNRSystem::CalReactivity(SNRSystem *sec,real kunp,real kp,vector<real> &rho,bool pr)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in SNRSystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in SNRSystem::CalReactivity.\n";
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
  //real ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  real *scttmp=new real[plnum];
  real *leaktmp=new real[plnum];
  for(int i=0;i<grp;i++){
    yld[i]=CalPerturbYieldTerm(sec,flag,i,kp);
    //yld[i]=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,kp);
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

  if(pr)cout<<"# Grp : Yield  : Absorption : Scattering : Leakage : P0 scattering\n";
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
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<le/leth<<"  "<<grpsum/leth<<" "<<sct[0][i]/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# Yield       :"<<yldsum<<"\n";
    cout<<"# Absorption  :"<<abssum<<"\n";
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

  return tot;
};

void SNRSystem::CalPerturbLeakScat(SNRSystem *sec,bool *flag,int g,real *sct,real *leak)
{
  real PI4_inv=1./PI4;

  for(int i=0;i<plnum;i++){
    sct[i]=0.; leak[i]=0.;
  };

  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      real dtot=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt(0).get_dat(g)
	               -mesh[m].GetMed()->GetMacxs().GetSigt(0).get_dat(g);
      real dabs=sec->GetMesh(m).GetMed()->GetMacxs().GetSiga().get_dat(g)
                       -mesh[m].GetMed()->GetMacxs().GetSiga().get_dat(g);
      real dsct=dtot-dabs;
      // total 
      sct[0]-=dsct*sec->GetMesh(m).GetFlux(0).get_dat(g)*
                           mesh[m].GetFlux(0).get_dat(g)*vol;
      int id=0;
      for(int l=1;l<=pl;l++){
	for(int mm=0;mm<=0;mm++){ // snr
	  leak[id]-=(2.*l+1.)*dtot
                   *sec->GetMesh(m).GetFlux(id+1).get_dat(g)
                           *mesh[m].GetFlux(id+1).get_dat(g)*vol;
	  id++;
	};
      };
      // (Higher leakage term)
      if(GetWrtflx(m)&&sec->GetWrtflx(m)){
        real tmp2=0.;
        int sn=GetQuadrature(g)->GetSN(); // snr
        for(int i=0;i<sn;i++){
	  real om =GetQuadrature(g)->GetOmega(i);
          real flp=sec->GetAFlux(m,g).get_dat(i);
          int id=0;
	  for(int l=0;l<=pl;l++){
	    for(int mm=0;mm<=0;mm++){ // snr
      	      real mom=GetQuadrature(g)->GetMoment(id,i);
              flp-=(2.*l+1.)*PI4_inv*mom*sec->GetMesh(m).GetFlux(id).get_dat(g);
              id++;
            };
          };
	  tmp2+=om*flp*aflux[m][g].get_dat(i);
        };
	leak[plnum-1]-=PI2*PI4*dtot*tmp2*vol;
      };
      // scattering term
      id=0;
      for(int l=0;l<=pl;l++){
	for(int mm=0;mm<=0;mm++){ // snr
  	  real tmp=0.;
	  //  for(int gg=g;gg<grp;gg++){
	  for(int gg=0;gg<grp;gg++){
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    tmp+=dsigs*mesh[m].GetFlux(id).get_dat(gg);
	  };
  	  sct[id]+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	  id++;
	};
      };
    };
  };
};

real SNRSystem::CalPerturbLeakScatSimplifiedNew(SNRSystem *sec,bool *flag,int g,int gg)
{
  real invpi4=1./PI4;

  real ret=0.;

  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      real dtot=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt().get_dat(g)
                       -mesh[m].GetMed()->GetMacxs().GetSigt().get_dat(g);
      real fact=PI2*PI4*dtot*vol;
      // total 
      if(dtot!=0.&&GetWrtflx(m)&&sec->GetWrtflx(m)){
        real tmp2=0.;
        //int sn=GetQuadrature(g)->GetSnnum(); // snrz
        int sn=GetQuadrature(g)->GetSN(); // snr
        for(int i=0;i<sn;i++){
	  real om =GetQuadrature(g)->GetOmega(i);
          real flp=sec->GetAFlux(m,g).get_dat(i);
	  tmp2+=om*flp*aflux[m][g].get_dat(i);
        };
	ret-=fact*tmp2;
      };
      // scattering term
      if(gg!=-1){
        int id=0;
        for(int l=0;l<=pl;l++){
  	  //for(int mm=0;mm<=l;mm++){ // snrz
  	  for(int mm=0;mm<=0;mm++){ // snr
  	    real tmp=0.;
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    if(dsigs!=0.)tmp=dsigs*mesh[m].GetFlux(id).get_dat(gg);
    	    ret+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	    id++;
	  };
	};
      };

    };
  };

  return ret;
};

void SNRSystem::CalSensitivity(SNRSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  // !!!
  // This method should NOT be delted,
  // since PKR(prompt k-ratio method module) uses this.
  // !!!

  CheckAdjointForward(sec);
  CheckSameMesh(sec);

  //bool fission_spectrum_matrix=false; // default:false
  bool fission_spectrum_matrix=true; // default:false
  if(fission_spectrum_matrix){
    cout<<"#\n#\n# !! Warning : Sensitivity calculation is carried out with fission spectrum matrix.\n#\n#\n";
  };

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

  real ip=CalPerturbDenominator(sec);
  // ++++++++++++++++++++++++++++++++
  if(fission_spectrum_matrix){
    ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
  };
  // ++++++++++++++++++++++++++++++++
  if(ip==0.){
    cout<<"Error in SNRsystem::CalSensitivity.\n";
    cout<<"Perturbation denominator is zero.\n";
    exit(0);
  };

  cout.setf(ios::scientific);
  cout.precision(5);

  // for reaction rate ratio sensitivity
  //ip=PI4;
  //ip=1.;

  cout<<grp<<"\n";

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];

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
      cout<<1<<"\n";
      cout<<nid<<"\n";
      cout<<"  18\n";
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
            absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
            real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
	    for(int l=0;l<=1;l++){
	      sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*micsigf);
	    };
	  };
        };
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int k=0;k<plnum;k++){
          re+=scttmp[k]+leaktmp[k];
        };
	// +++++++++++++++++++++++++++++++++
        if(!fission_spectrum_matrix){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
	}else{
          re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,k2);
	};
	// ++++++++++++++++++++++++++++++++
        re+=CalPerturbAbsorptionTerm(sec,flag,i);
        re/=ip;
	cout<<"  "<<re<<"\n";
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    for(int l=1;l<=1;l++){
              sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	    };
	  };
	};
      };

      // Nu
      cout<<1<<"\n";
      cout<<nid<<"\n";
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
	// +++++++++++++++++++++++++
	if(!fission_spectrum_matrix){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
	}else{
          re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,k2);
	};
	// ++++++++++++++++++++++++++++
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
      cout<<nid<<"\n";
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
    cout<<nid<<"\n";
    cout<<"  102\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
          tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
	  for(int l=0;l<=1;l++){
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
      re/=ip;
      cout<<"  "<<re<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	  for(int l=1;l<=1;l++){
            sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	  };
	};
      };
    };

    // Scattering P0
    for(int ii=0;ii<3;ii++){
      cout<<2<<"\n";
      cout<<nid<<"\n";
      //
      if(ii==0){cout<<2<<"\n";}
      else if(ii==1){cout<<4<<"\n";}
      else {cout<<16<<"\n";};
      //
      for(int i=0;i<grp;i++){
        for(int k=i;k<grp;k++){
	  if((ii!=0)||(nid<1000)||(ii==0&&k<=i+2)){
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
	          //micsigs1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(1).get_dat(i,k);
	          break;
	        };
	        real a1=0.;
		//real s0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
		//real s1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(1).get_sumx().get_dat(i);
		if(s0!=0.)a1=s1/s0;
	        //if(micsigs!=0.){a1=micsigs1/micsigs;}
	        //else{a1=0.;};
		// !
                if(ii!=0.)micsigs=delta;
	        real dtot;
	        if(ii==2){
  	  	  dtot=micsigs*0.5;
	        }else{
   	          dtot=micsigs;
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
	        sec->GetMed(j).GetMacxs().GetSigs(1).add_data(i,k,den*a1*micsigs);
	        for(int l=0;l<=1;l++){
	          sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*dtot);
	        };
              };
            };
            real re=0.;
            CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
            for(int l=0;l<plnum;l++){
  	      re+=scttmp[l]+leaktmp[l];
            };
            re/=ip;
            cout<<"  "<<re<<"\n";
            for(int j=0;j<nmed;j++){
	      if(med[j].ExistNuclide(nid)){
                sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	        for(int l=1;l<=1;l++){
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

    /*
    // +++ Elastic matrix : P_1-P_L
    for(int l=1;l<=pl;l++){
      real delta=1.0;
      cout<<"  "<<2<<"\n";
      cout<<"  "<<nid<<"\n";
      cout<<"  "<<20000+l<<"\n";
      for(int i=0;i<grp;i++){
        int kmax=i+2;
	if(kmax>grp)kmax=grp;
        for(int k=i;k<kmax;k++){
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
	      sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*delta);
            };
	  };
          real re=0.;
          CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
          for(int ll=0;ll<plnum;ll++){
  	    re+=scttmp[ll]+leaktmp[ll];
          };
          re/=ip;
          cout<<"  "<<re<<"\n";
          for(int j=0;j<nmed;j++){
	    if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
  	      sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,-den*delta);
	    };
	  };
	};
      };
    };
    */

    // Elastic-p1
    cout<<1<<"\n";
    cout<<nid<<"\n";
    cout<<"  251\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
	  int maxg=grp;
	  if(nid>1000&&i+3<maxg)maxg=i+3;
	  //for(int k=i;(k<=i+2&&k<grp);k++){
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
      re/=ip;
      cout<<"  "<<re*3.<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          int maxg=grp;
	  if(nid>1000&&i+3<maxg)maxg=i+3;
	  //for(int k=i;(k<=i+2&&k<grp);k++){
	  for(int k=i;k<maxg;k++){
	    sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigs1org[j*grp+k]);
	  };
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

SensitivityData SNRSystem::CalSensitivityNew(SNRSystem *sec,real keff,int nucnum,int *nucid,bool fiss_matrix,bool ipcal)
{
  bool fission_spectrum_matrix=fiss_matrix; // default:false
  if(fission_spectrum_matrix){
    cout<<"#\n#\n# !! Warning : Sensitivity calculation is carried out with fission spectrum matrix.\n#\n#\n";
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
  real *sigsorg=new real[nmed*grp];
  real *sigs1org=new real[nmed*grp];

  real *fiss_frac=new real[nmed];

  real delta=1.;

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o

  bool *flag=new bool[TotM];
  
  real ip=1.;
  if(ipcal){
    ip=CalPerturbDenominator(sec);
  //cout<<"# Perturb denominator is : "<<ip<<"\n";
  };

  // ++++++++++++++++++++++++++++++++
  if(fission_spectrum_matrix){
    ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
  };
  // ++++++++++++++++++++++++++++++++
  if(ip==0.){
    cout<<"# Error in SNRsystem::CalSensitivityNew.\n";
    cout<<"# Perturbation denominator is zero.\n";
    exit(0);
  };

  // for reaction rate ratio sensitivity
  //ip=PI4;
  //ip=1.;

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
              totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
              tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
              real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	      real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(nu).get_dat(i);
              sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
              sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
	      for(int l=0;l<=1;l++){
	        sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*micsigf);
	      };
	    };
          };
          real re=0.;
          CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
          for(int k=0;k<plnum;k++){
            re+=scttmp[k]+leaktmp[k];
          };
  	  // +++++++++++++++++++++++++++++++++
          //re+=CalPerturbYieldTerm(sec,flag,i,keff);
          if(!fission_spectrum_matrix){
            re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  }else{
            re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,keff);
	  };
	  // ++++++++++++++++++++++++++++++++
          re+=CalPerturbAbsorptionTerm(sec,flag,i);
          re/=ip;
          sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
	    if(med[j].ExistNuclide(nid)){
  	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	      sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
              sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	      for(int l=1;l<=1;l++){
                sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	      };
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
          //re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  // +++++++++++++++++++++++++
	  if(!fission_spectrum_matrix){
            re+=CalPerturbYieldTerm(sec,flag,i,keff);
	  }else{
            re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,keff);
	  };
  	  // ++++++++++++++++++++++++++++
	  re/=ip;
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
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
	    real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
	    for(int l=0;l<=1;l++){
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
        re/=ip;
        sns1d.put_data(i,re);
        for(int j=0;j<nmed;j++){
  	  if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	    for(int l=1;l<=1;l++){
              sec->GetMed(j).GetMacxs().GetSigt(l).put_data(i,tot1org[j]);
	    };
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
          for(int k=i;k<grp;k++){
	    if((ii!=0)||(nid<100000)||(ii==0&&k<=i+2)){
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
              //real re=CalPerturbLeakScatSimplifiedNew(sec,flag,i,k);
              re/=ip;
              sns2d.put_data(i,k,re);
              for(int j=0;j<nmed;j++){
	        if(med[j].ExistNuclide(nid)){
                  sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	          for(int l=1;l<=1;l++){
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
      
      // +++ Elastic-p1
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
        re/=ip;
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
      sens.PutSensitivity1D(nid,251,sns1d);

      // +++ Elastic-p2
      if(plnum>1){
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
    	  if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
	    int maxg=grp;
	    if(nid>1000&&i+3<maxg)maxg=i+3;
  	    for(int k=i;k<maxg;k++){
	      sigs1org[j*grp+k]=sec->GetMed(j).GetMacxs().GetSigs(2).get_dat(i,k);
	      real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	      sec->GetMed(j).GetSigs(2).add_data(i,k,den*delta*micsigs);
	    };
          };
        };
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int k=0;k<plnum;k++){
          re+=scttmp[k]+leaktmp[k];
        };
        re/=ip;
        sns1d.put_data(i,re*5.);
        for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
            int maxg=grp;
	    if(nid>1000&&i+3<maxg)maxg=i+3;
	    for(int k=i;k<maxg;k++){
	      sec->GetMed(j).GetMacxs().GetSigs(2).put_data(i,k,sigs1org[j*grp+k]);
	    };
 	  };
        };
      };
      sens.PutSensitivity1D(nid,252,sns1d);
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

  return sens;
};


SensitivityData SNRSystem::CalSensitivityRRR(SNRSystem &snr,int nume_nuc,enum xstype nume_xs,int denom_nuc,enum xstype denom_xs,int pos,int nucnum,int *nucid,real keff)
{
  int itermax=20;

  //
  real vol=snr.GetMesh(pos).GetVolume();
  real rr0=snr.GetMesh(pos).GetFlux()*med[0].GetNuclide(denom_nuc).GetMicxs().GetData1d(denom_xs)*vol;
  real rr1=snr.GetMesh(pos).GetFlux()*med[0].GetNuclide(nume_nuc).GetMicxs().GetData1d(nume_xs)*vol;
  cout.setf(ios::showpoint);
  cout.precision(8);
  real rr=rr1/rr0;
  cout<<"# Reation rate : "<<rr<<"\n"; 
  rr0=1./rr0;
  rr1=1./rr1;
  GroupData1D src0,src1;
  src0=med[0].GetNuclide(denom_nuc).GetMicxs().GetData1d(denom_xs)*rr0;
  src1=med[0].GetNuclide(nume_nuc).GetMicxs().GetData1d(nume_xs)*rr1;

  // Generalized adjoint flux calculation
  SNRSystem snr_g1=*this;
  SNRSystem snr_g2=*this;
  SNRSystem snr_st=*this;

  snr_g1.SetArray();
  snr_g2.SetArray();
  snr_st.SetArray();

  snr_g1.PutIsotropicSourceParVolume(pos,src0);
  snr_g2.PutIsotropicSourceParVolume(pos,src1);

  for(int iter=0;iter<itermax;iter++){

    snr_g1.CalFixedSource(1e-4,1,false);
    snr_g2.CalFixedSource(1e-4,1,false);

    snr_g2.NegFlux(snr_g1);
    snr_g2.CalFissionSrcAllMesh();

    real sum1=CalPerturbDenominatorFromSource(&snr);
    real sum2=snr_g2.CalPerturbDenominatorFromSource(&snr);
    real fact=sum2/sum1;

    snr_g2.NegFlux(*this,fact); // Fundamental mode extraction
    snr_g2.CalFissionSrcAllMesh();
    snr_st.AddFlux(snr_g2); 

    snr_g1.SetZeroScatSrc(); // for negative source 
    snr_g2.SetZeroScatSrc(); // for positive source
    for(int i=0;i<snr_g1.GetTotM();i++){
      GroupData1D tmp1(grp);
      GroupData1D tmp2(grp);
      tmp1.set_zero();
      tmp2.set_zero();
      real fsrc=snr_g2.GetMesh(i).GetFissionSrc();
      if(fsrc<0.){
        for(int j=0;j<grp;j++){
          tmp1.put_data(j,-fsrc*med[0].GetMacxs().GetData1d(nusigf).get_dat(j));
        };
      }else{
        for(int j=0;j<grp;j++){
          tmp2.put_data(j,fsrc*med[0].GetMacxs().GetData1d(nusigf).get_dat(j));
        };
      };
      real vol_inv=1./snr_g1.GetMesh(i).GetVolume();
      tmp1=tmp1*vol_inv/keff;
      tmp2=tmp2*vol_inv/keff;
      snr_g1.PutIsotropicSourceParVolume(i,tmp1);
      snr_g2.PutIsotropicSourceParVolume(i,tmp2);
    };
  };

  snr_st.CalFissionSrcAllMesh();
  SensitivityData sens=snr_st.CalSensitivityNew(&snr,1.,nucnum,nucid,false,false);
  // The first [false] is for fission spectrum matrix calculation and 
  // the second [false] is for perturbation denominator calculation
  sens.PutName("dummy","rrr","dummy");
  sens.PutValue(rr);

  // +++ Sensitivity calculation for direct term
  SensitivityData sens_dir;
  sens_dir.PutGroup(grp);
  for(int i=0;i<nucnum;i++){
    if(nucid[i]==nume_nuc){
      GroupData1D snsdat(grp);
      snsdat=snr.GetMesh(pos).GetFlux().mult(med[0].GetNuclide(nume_nuc).GetMicxs().GetData1d(nume_xs));
      snsdat=snsdat*(1./snsdat.get_sum());
      int mtid=18;
      if(nume_xs==sigc)mtid=102;
      sens_dir.PutSensitivity1D(nume_nuc,mtid,snsdat);
    };
    if(nucid[i]==denom_nuc){
      GroupData1D snsdat(grp);
      snsdat=snr.GetMesh(pos).GetFlux().mult(med[0].GetNuclide(denom_nuc).GetMicxs().GetData1d(denom_xs));
      snsdat=snsdat*(-1./snsdat.get_sum());
      int mtid=18;
      if(denom_xs==sigc)mtid=102;
      sens_dir.PutSensitivity1D(denom_nuc,mtid,snsdat);
    };
  };
  sens.AddSensitivityData(sens_dir);

  return sens;
};

/*
void SNRSystem::CalGPT_SNR(real keff,int nuc,int *nuc_id,enum xstype *xs,bool *on_mesh,real rr,int itermax)
{
  GroupData1D src(grp);

  SetZeroScatSrc();
  for(int i=0;i<TotM;i++){
    if(on_mesh[i]){
      src.set_zero();
      for(int j=0;j<nuc;j++){
	if(mesh[i].GetMed()->ExistNuclide(nuc_id[j])){
          src=src+mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetMicxs().GetData1d(xs[j])*mesh[i].GetMed()->GetNuclide(nuc_id[j]).GetDensity()/rr;
        };
      };
      PutIsotropicSourcePerUnitVolume(i,src);
    };
  };
  CalGPT_SNR(keff,1e-4,itermax);
};

void SNRSystem::CalGPT_SNR(real keff,real eps_f,int itermax)
{
  CheckUpScattering();

  vector<GroupData1D> gpt_flx_store(TotM);
  vector<GroupData1D> gpt_flx_old(TotM);
  for(int i=0;i<TotM;i++){
    gpt_flx_store[i].put_imax(grp);
    gpt_flx_old[i].put_imax(grp);
  };

  vector< vector<GroupData1D> > gpt_aflx_store;
  gpt_aflx_store.resize(TotM);
  for(int i=0;i<TotM;i++){
    if(wrtflx[i]){
      gpt_aflx_store[i].resize(grp);
      for(int g=0;g<grp;g++){
        int totsn=GetQuadrature(g)->GetSN();
	gpt_aflx_store[i][g].put_imax(totsn);
      };
    };
  };

  int iter=0;
  real sumold=0.;
  real sum=0.;
  for(int k=0;k<itermax;k++){

    sumold=sum;
    sum=0.;
    for(int m=0;m<TotM;m++){
      if(k!=0)mesh[m].CalFissionSrcAdjoint();
      sum+=mesh[m].GetFissionSrc();
    };
    //if(k!=0)cout<<"     "<<sumold<<" -> "<<sum<<" : "<<sum/sumold<<"\n";
    
    if(k==0){
      if(IsUpScatterSystem){
        CalFixedSourceUpScat();
      }else{
        CalFixedSource(1e-4,1,false);
      };
    }else{
      SetZeroScatSrc();
      if(IsUpScatterSystem){
        CalFixedFissionSourceUpScat(keff);
      }else{
        CalFixedFissionSource(keff,k);
      };
    };
    iter++;

    real errmax=0.;
    real flxmax=0.;
    for(int m=0;m<TotM;m++){
      for(int g=0;g<grp;g++){
        real tmp=mesh[m].GetFlux().get_dat(g);
        real flxold=gpt_flx_old[m].get_dat(g);
	//cout<<m<<" "<<g<<" "<<flxold<<" -> "<<tmp<<" "<<tmp/flxold<<"\n";
        if(flxold!=0.&&k!=0){
	  //if(flxold>0.&&k!=0){
          real err=fabs(tmp/flxold-1.);
          if(err>errmax)errmax=err;
        };
	if(fabs(tmp)>flxmax)flxmax=fabs(tmp);
        gpt_flx_old[m].put_data(g,tmp);
      };
    };

    //cout<<"# CALGPT : "<<k<<" "<<flxmax<<"\n";
    //cout<<"# CALGPT : "<<k<<" "<<errmax<<" (Maximum flux : "<<flxmax<<")\n";

    if(errmax<eps_f&&k!=0){
      break;
    };

    if(k==0){
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m].copy(gpt_flx_old[m]);
        if(wrtflx[m])gpt_aflx_store[m]=aflux[m];
      };
    }else{
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m]=gpt_flx_store[m]+gpt_flx_old[m];
	if(wrtflx[m]){
	  for(int g=0;g<grp;g++){
	    gpt_aflx_store[m][g]=gpt_aflx_store[m][g]+aflux[m][g];
	  };
	};
      };
    };

    if(k==itermax-1){
      cout<<"# Warning in GeneralSystem::CalGPT.\n";
      cout<<"# Iteration reaches to maximum ("<<itermax<<").\n";
      cout<<"# Iteration residual : "<<errmax<<"\n";
    };
  };

  for(int m=0;m<TotM;m++){
    mesh[m].GetFlux().copy(gpt_flx_store[m]-gpt_flx_old[m]*(iter-1));
    if(wrtflx[m]){
      for(int g=0;g<grp;g++){
	aflux[m][g]=gpt_aflx_store[m][g]-aflux[m][g]*(iter-1);
      };
    };
  };

  //mesh[0].GetFlux().show_self();

};
*/


void SNRSystem::CalSensitivityScatteringMatrix(SNRSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  if(GetGeneralOption().Forward()){
    cout<<"Error in SNRSystem::CalSensitivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in SNRSystem::CalSensitivity.\n";
    cout<<"Forward flux should be calculated in perturbated system.\n";
    exit(0);
  };

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o
  
  cout.setf(ios::scientific);
  cout.precision(5);

  CheckSameMesh(sec); // o

  real ip=CalPerturbDenominator(sec);
  if(ip==0.){
    cout<<"Error in SNRsystem::CalSensitivity.\n";
    cout<<"Perturbation denominator is zero.\n";
    exit(0);
  };

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  cout<<grp<<"\n";
  for(int nc=0;nc<nucnum;nc++){
    int nid=nucid[nc];
    //int nidendf=TranslateNuclideIDFromJFS(nid);
    int nidendf=nid;
    //if(jfsid){nidendf=TranslateNuclideIDFromJFS(nid);};

    for(int i=0;i<TotM;i++){
      if(mesh[i].GetMed()->ExistNuclide(nid)){
        flag[i]=true;
      }else{
	flag[i]=false;
      };
    };

    // +++ Elastic matrix : P_1-P_L
    for(int l=1;l<=pl;l++){
      real delta=1.0;
      cout<<"  "<<2<<"\n";
      cout<<"  "<<nidendf<<"\n";
      cout<<"  "<<20000+l<<"\n";
      for(int i=0;i<grp;i++){
        int kmax=i+2;
	if(kmax>grp)kmax=grp;
        for(int k=i;k<kmax;k++){
          for(int j=0;j<nmed;j++){
            if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
	      sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*delta);
            };
	  };
          real re=0.;
          CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
          for(int ll=0;ll<plnum;ll++){
  	    re+=scttmp[ll]+leaktmp[ll];
          };
          re/=ip;
          cout<<"  "<<re<<"\n";
          for(int j=0;j<nmed;j++){
	    if(med[j].ExistNuclide(nid)){
              real den=med[j].GetNuclide(nid).GetDensity();
  	      sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,-den*delta);
	    };
	  };
	};
      };
    };
  };

  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;
};

void SNRSystem::CalSensitivityFixedSource(SNRSystem *sec,int nucnum,int *nucid)
{
  CheckAdjointForward(sec);

  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  vector< vector<real> > sigsorg((pl+1),vector<real>(nmed*grp));

  real delta=1.;
  real PI4_INV=1./PI4;
  // Because transport theory-based perturbation results are multiplied by 4\PI in CBG.
  // This is due to consistency in perturbation denominator between transport and diffusion theories.

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o
  
  cout.setf(ios::scientific);
  cout.precision(5);

  CheckSameMesh(sec); // o

  bool *flag=new bool[TotM];

  cout<<grp<<"\n";
  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];

    for(int i=0;i<TotM;i++){
      if(mesh[i].GetMed()->ExistNuclide(nid)){
        flag[i]=true;
      }else{
	flag[i]=false;
      };
    };

    // Capture
    cout<<"  1\n";
    cout<<"  "<<nid<<"\n";
    cout<<"  102\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
          sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigc);
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      cout<<"  "<<re*PI4_INV<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	};
      };
    };

    // P0 Elastic scattering
    cout<<"  2\n";
    cout<<"  "<<nid<<"\n";
    cout<<"  2\n";
    for(int i=0;i<grp;i++){
      int maxg=i+10;
      if(maxg>grp)maxg=grp;
      for(int k=i;k<maxg;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	    real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
	    for(int l=0;l<=pl;l++){
              sigsorg[l][j]=sec->GetMed(j).GetMacxs().GetSigs(l).get_dat(i,k);
	      real a1=1.;
	      if(l!=0){
	        real st1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(l).get_sumx().get_dat(i);
                a1=st1/st0;
	      };
              sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*a1*micsigs);
	    };
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigs);
          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        cout<<"  "<<re*PI4_INV<<"\n";
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetMacxs().GetSigs(l).put_data(i,k,sigsorg[l][j]);
	    };
	  };
	};
      };
    };

    // Inelastic scattering
    cout<<"  2\n";
    cout<<"  "<<nid<<"\n";
    cout<<"  4\n";
    int maxinel=119;
    cout<<"  "<<maxinel<<"\n";
    for(int i=0;i<maxinel;i++){
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	    /*
            real micsigs=delta;
	    real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_sumx().get_dat(i);
	    for(int l=0;l<=pl;l++){
              sigsorg[l][j]=sec->GetMed(j).GetMacxs().GetSigs(l).get_dat(i,k);
	      real a1=1.;
	      if(l!=0){
	        real st1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(l).get_sumx().get_dat(i);
                a1=st1/st0;
	      };
              sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*a1*micsigs);
	    };
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigs);
	    */

            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta);

          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        cout<<"  "<<re*PI4_INV<<"\n";
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    /*
	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetMacxs().GetSigs(l).put_data(i,k,sigsorg[l][j]);
	    };
	    */
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
	};
      };
    };

    // n2n
    cout<<"  2\n";
    cout<<"  "<<nid<<"\n";
    cout<<"  16\n";
    int maxn2n=119;
    cout<<"  "<<maxn2n<<"\n";
    for(int i=0;i<maxn2n;i++){
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta*0.5);
          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        cout<<"  "<<re*PI4_INV<<"\n";
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
	};
      };
    };

    // Elastic P1
    cout<<"  1\n";
    cout<<"  "<<nid<<"\n";
    cout<<"  251\n";
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  real den=med[j].GetNuclide(nid).GetDensity();
	  for(int k=i;k<grp;k++){
	    sigsorg[1][j*grp+k]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
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
      cout<<"  "<<re*3.*PI4_INV<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  for(int k=i;k<grp;k++){
	    sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigsorg[1][j*grp+k]);
	  };
	};
      };
    };

  };

  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;
  delete [] absorg;
  delete [] totorg;
};

SensitivityData SNRSystem::CalSensitivityNewFixedSource(SNRSystem *sec,real val,int nucnum,int *nucid)
{
  CheckAdjointForward(sec);

  SensitivityData sens;
  sens.PutValue(val);
  sens.PutGroup(grp);
  sens.GetEnband().copy(med[0].GetEnband());

  GroupData1D sns1d(grp);
  GroupData2D sns2d(grp,grp);

  real *absorg=new real[nmed];
  real *totorg=new real[nmed];
  vector< vector<real> > sigsorg((pl+1),vector<real>(nmed*grp));

  real delta=1.;
  real PI4_INV=1./PI4;
  // Because transport theory-based perturbation results are multiplied by 4\PI in CBG.
  // This is due to consistency in perturbation denominator between transport and diffusion theories.

  real *scttmp=new real[plnum]; // o
  real *leaktmp=new real[plnum]; // o
  
  cout.setf(ios::scientific);
  cout.precision(5);

  CheckSameMesh(sec); // o

  bool *flag=new bool[TotM];

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    cout<<"# Sensitivity calculation for nuclide : "<<nid<<" ("<<nc<<"/"<<nucnum<<")\n";

    for(int i=0;i<TotM;i++){
      if(mesh[i].GetMed()->ExistNuclide(nid)){
        flag[i]=true;
      }else{
	flag[i]=false;
      };
    };

    // Capture
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigc);
          sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigc);
        };
      };
      real re=0.;
      CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
      for(int k=0;k<plnum;k++){
        re+=scttmp[k]+leaktmp[k];
      };
      re+=CalPerturbAbsorptionTerm(sec,flag,i);
      sns1d.put_data(i,re*PI4_INV/val);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	  sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	};
      };
    };
    sens.PutSensitivity1D(nid,102,sns1d);

    // P0 Elastic scattering
    sns2d.set_zero();
    for(int i=0;i<grp;i++){
      for(int k=i;k<grp;k++){

	bool zeroxs=true;
	for(int j=0;j<nmed;j++){
	  if(med[j].ExistNuclide(nid)){
	    real den=med[j].GetNuclide(nid).GetDensity();
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	    if(den>0.&&micsigs!=0.)zeroxs=false;
	  };
	};

	if(!zeroxs){

        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	    real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
	    for(int l=0;l<=pl;l++){
              sigsorg[l][j]=sec->GetMed(j).GetMacxs().GetSigs(l).get_dat(i,k);
	      real a1=1.;
	      if(l!=0){
	        real st1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(l).get_sumx().get_dat(i);
                a1=st1/st0;
	      };
              sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*a1*micsigs);
	    };
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigs);
          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        sns2d.put_data(i,k,re*PI4_INV/val);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetMacxs().GetSigs(l).put_data(i,k,sigsorg[l][j]);
	    };
	  };
	};

	};

      };
    };
    sens.PutSensitivity2D(nid,2,sns2d);

    // Inelastic scattering
    int egrp=0;
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
        if(med[j].ExistNuclide(nid)){
          real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_sumx().get_dat(i);
          if(micsigs!=0.)egrp=i;
	};
      };
    };
    egrp+=5; // 
    if(egrp>grp)egrp=grp;

    sns2d.set_zero();
    for(int i=0;i<egrp;i++){
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	    /*
            real micsigs=delta;
	    real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_sumx().get_dat(i);
            int ineldim=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetDim2d(siginel);
            for(int l=0;l<=pl;l++){
              sigsorg[l][j]=sec->GetMed(j).GetMacxs().GetSigs(l).get_dat(i,k);
	      if(l<=ineldim){
	        real a1=1.;
	        if(l!=0){
	          real st1=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(l).get_sumx().get_dat(i);
                  a1=st1/st0;
	        };
                sec->GetMed(j).GetMacxs().GetSigs(l).add_data(i,k,den*a1*micsigs);
	      };
	    };
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*micsigs);
	    */
	    
            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta);

          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        sns2d.put_data(i,k,re*PI4_INV/val);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    /*
	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetMacxs().GetSigs(l).put_data(i,k,sigsorg[l][j]);
	    };
	    */
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
	};
      };
    };
    sens.PutSensitivity2D(nid,4,sns2d);

    // n2n
    egrp=0;
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
        if(med[j].ExistNuclide(nid)){
          real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_sumx().get_dat(i);
          if(micsigs!=0.)egrp=i;
	};
      };
    };
    egrp+=5; // 
    if(egrp>grp)egrp=grp;

    sns2d.set_zero();
    for(int i=0;i<egrp;i++){
      for(int k=i;k<grp;k++){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta*0.5);
          };
	};
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int l=0;l<plnum;l++){
          re+=scttmp[l]+leaktmp[l];
        };
        sns2d.put_data(i,k,re*PI4_INV/val);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
	};
      };
    };
    sens.PutSensitivity2D(nid,16,sns2d);

    // Elastic P1
    for(int i=0;i<grp;i++){
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  real den=med[j].GetNuclide(nid).GetDensity();
	  for(int k=i;k<grp;k++){
	    sigsorg[1][j*grp+k]=sec->GetMed(j).GetMacxs().GetSigs(1).get_dat(i,k);
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
      sns1d.put_data(i,re*3.*PI4_INV/val);
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
	  for(int k=i;k<grp;k++){
	    sec->GetMed(j).GetMacxs().GetSigs(1).put_data(i,k,sigsorg[1][j*grp+k]);
	  };
	};
      };
    };
    sens.PutSensitivity1D(nid,251,sns1d);

  };

  delete [] scttmp;
  delete [] leaktmp;
  delete [] flag;
  delete [] absorg;
  delete [] totorg;

  return sens;
};
