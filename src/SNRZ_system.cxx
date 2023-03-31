#include <cstdlib>
#include "SNRZ_system.h"

using namespace std;

void SNRZSystem::Init(int n,int g,int i)
{
  psys.Init(n,g,i);
  GeneralSystem::Init(n,g,i);
  pl = -1;
  nquad = 1;
  quad.resize(nquad);
  quadid.resize(grp,0);
  for(int i=0;i<grp;i++){quadid[i]=0;};
  transport=true;
  etransport=false;
  name="SNRZ";
  cmfdimp=true; // CMFD is not implemented
  dsa=true;
  step_d=false;
  array_set=false;
  leakage.put_imax(grp);
  itmax_inner=99;
};

SNRZSystem::~SNRZSystem()
{
};

void SNRZSystem::SetQuadratureNum(int num)
{
  nquad=num;
  quad.resize(nquad);
};

void SNRZSystem::SetQuadratureID(int *gbnd, int *id)
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

void SNRZSystem::SetQuadrature(SNRZQuadrature *qinp,int id)
{
  if(id<0||id>=nquad){
    cout<<"Error in SetQuadratureID.\n";
    cout<<"Chack your ID number.\n";
    exit(0);
  };
  quad[id]=qinp;
};

void SNRZSystem::SetArray()
{
  if(array_set)return;
  array_set=true;

  if(pl==-1){
    cout<<"# Do 'PutPL for system class'.\n";
    exit(0);
  };

  for(int i=0;i<grp;i++){
    if(!quad[quadid[i]]->Exist()){
      cout<<"# Quadrature is not defined yet.\n";
      exit(0);
    };
  };

  for(int i=0;i<TotM;i++){
    mesh[i].PutPL(pl);
    mesh[i].PutGrp(grp);
    if(wrtflx[i]){
      aflux[i].resize(grp);
      for(int j=0;j<grp;j++){
        aflux[i][j].put_imax(quad[quadid[j]]->GetSnnum()); // snrz
      };
    };
  };

  if(pl>0){
    if(!etransport&&transport){
    for(int i=0;i<nmed;i++){
      for(int j=1;j<pl+1;j++){
        for(int k=0;k<grp;k++){
	  real tmp=med[i].GetMacxs().GetData2d(sigs,j).get_dat(k,k);
	  real tmp2=med[i].GetDataSigt(0,k)-med[i].GetDataSigt(1,k);
	  real tmp3=tmp+(j*2+1.)*tmp2;
	  real tmp4=med[i].GetSigs(0).get_dat(k,k);
	  
	  if(j==1&&fabs(tmp3*0.33333333)>tmp4){
	    if(print)cout<<"#   ... warning... mu (sp1/sp0) is out of range : "<<tmp3*0.33333333/tmp4<<"\n";
	    if(print)cout<<"#         medium "<<i<<"  group "<<k<<"\n";
	  }
          else{
            med[i].GetSigs(j).put_data(k,k,tmp3);
	  };

          //med[i].GetMacxs().GetData2d(sigs,j).put_data(k,k,tmp3);
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

  for(int i=0;i<nquad;i++){
    quad[i]->CheckOrthogonalityTo00();
  };

  // For CMFD
  curxp.resize(TotM,0.);
  curxn.resize(TotM,0.);
  curyp.resize(TotM,0.);
  curyn.resize(TotM,0.);
};

void SNRZSystem::PutWriteFlux(int med)
{
  int ind=0;
  int indall=0;
  for(int i=0;i<mi.GetYF();i++){
    for(int j=0;j<mi.GetXF();j++){
      int tmp=mi.GetFMat(indall);
      if(tmp!=-1){
	if(tmp==med||med==-1)wrtflx[ind]=true;
	ind++;
      };
      indall++;
    };
  };
};


void SNRZSystem::PutCartMeshInfo(CartMeshInfo cm, string geom)
{
  if(geom!="Cylinder")exit(0);

  mi=cm;

  int xf=mi.GetXF();
  int yf=mi.GetYF();
  int zf=1; // o

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
  //int index=0;
  for(int y=0;y<yf;y++){
    //real tmp=0.;
    for(int x=0;x<xf;x++){
      //real tmp2=tmp;
      int tm=mi.GetFMat(tmp2);
      //tmp+=mi.GetFMeshL(0,x);
      //mesh[index].PutDimCylinder(mi.GetFMeshL(0,x),tmp2,tmp,mi.GetFMeshL(1,y));
      //int tm=mi.GetFMat(index);
      if(tm>=nmed){
	cout<<"Error in SNRZSystem::PutCartMeshInfo.\n";
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
      //mesh[index].PutMedium(&med[tm]);
      //meshid[0][y][x]=index;
      //index++;
    };
  };
  //TotM=xf*yf;
  TotM=tmp;
  if(print)cout<<"#*** Total Mesh : "<<TotM<<"\n";
  mesh.resize(TotM);
  wrtflx.resize(TotM,false);
  aflux.resize(TotM);

  int index=0;
  int ind2=0;
  for(int y=0;y<yf;y++){
    real tmp=0.;
    for(int x=0;x<xf;x++){
      real tmpold=tmp;
      tmp+=mi.GetFMeshL(0,x);
      int tm=mi.GetFMat(ind2);
      ind2++;
      if(tm!=-1){
        mesh[index].PutDimCylinder(mi.GetFMeshL(0,x),tmpold,tmp,mi.GetFMeshL(1,y));	
	mesh[index].PutMedium(&med[tm]);
	meshid[0][y][x]=index;
	index++;
      }else{
	meshid[0][y][x]=-1;
      };
    };
  };

  /*
  for(int y=0;y<yf;y++){
    for(int x=0;x<xf;x++){
      WriteOut(meshid[0][y][x],4);
      cout<<" ";
    };
    cout<<"\n";
  };
  exit(0);
  */ 
 
  EdgeCalculation();
 
  SetCoarseMeshInfo();

  for(int i=0;i<4;i++){
    if(mi.GetBC(i)==2)BC[i]=2; // Vacuum
    if(mi.GetBC(i)==1)BC[i]=1; // Reflective
    if(BC[i]<0||BC[i]>2){
      cout<<"Boundary Condition Error in PLOS.\n";
      exit(0);
    };
  };
  BC[4]=0;
  BC[5]=0;

  if(BC[0]==2){
    cout<<"# Boundary condition at left side should be 'Reflective'\n";
    cout<<"# This condition is modified.\n";
    BC[0]=1;
  };
  if(BC[1]==1){
    cout<<"# Boundary condition at right side should be 'Vacuum'\n";
    cout<<"# This condition is modified.\n";
    BC[1]=0;
  };

  if(dsa)psys.PutCartMeshInfo(cm,geom);

  // for CMFD
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  int zr=mi.GetZC();
  //int zr=0; // o

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

void SNRZSystem::CalFluxInSource(real *Src,real *SSrc,real *flxmom,real *selfsc,int qid,int iss)
{
  vector<real> moment(plnum);
  for(int i=0;i<plnum;i++){
    moment[i]=quad[qid]->GetMoment(i,iss);
  };

  real ssc;
  int ind=0;
  int ind2=0;
  int ind3=0;
  for(int i=0;i<TotM;i++){
    real tmp=0.;
    int index=0;
    for(int l=0;l<=pl;l++){
      // chibatmp
      real sigs=selfsc[ind3++];
      for(int m=0;m<=l;m++){ // for 2D
	// chibatmp
        ssc=sigs*flxmom[ind2]+Src[ind2];
        //ssc=Src[ind2];
	ind2++;
        //tmp+=ssc*quad[qid]->GetMoment(index,iss);  
        tmp+=ssc*moment[index];
	index++;
      };
    };
    //if(tmp<0.)tmp=0.;
    SSrc[i]=tmp;
    ind++;
  };
};

void SNRZSystem::CalFluxInSourceSP(real *Src,real *SSrc,real *flxmom,real *selfsc,int qid,real xi,real mu)
{
  vector<real> moment(plnum);
  int index=0;
  for(int l=0;l<=pl;l++){
    for(int m=0;m<=l;m++){
      moment[index++]=quad[qid]->GetMoment(l,m,mu,xi);
    };
  };

  int ind=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    real srcin=0.;
    int index=0;
    for(int l=0;l<=pl;l++){
      // chibatmp
      real sigs=selfsc[ind++];
      for(int m=0;m<=l;m++){ // for 2D
	// chibatmp
        srcin+=(sigs*flxmom[ind2]+Src[ind2])*moment[index++];
        //srcin+=(sigs*flxmom[ind2]+Src[ind2])*quad[qid]->GetMoment(l,m,mu,xi);
        //srcin+=(Src[ind2])*quad[qid]->GetMoment(l,m,mu,xi);
	ind2++;
      };
    };
    if(srcin<0.)srcin=0.;
    SSrc[i]=srcin; 
  };
};

void SNRZSystem::PutSigmaForInnerIteration(int g,real *sigtvol,real *sigsself,real *sigt)
{
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    real sigttmp=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
    sigt[i]=sigttmp;
    sigtvol[i]=sigttmp*vol;
    // chibatmp
    real volpi4=vol*INV_PI4;
    for(int l=0;l<=pl;l++){
      // chibatmp
      sigsself[ind++]=mesh[i].GetMed()->GetMacxs().GetSigs(l).get_dat(g,g)*volpi4;
    };
  };
};

void SNRZSystem::InitialFluxGuessInnerIteration(int g,real *flxmom)
{
  int ind=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){
      flxmom[ind]=mesh[i].GetFlux(j).get_dat(g);
      ind++;
    };
  };
};

void SNRZSystem::PutSourceInnerIteration(real *Src)
{
  int ind=0;
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      Src[ind++]=mesh[i].GetSrcin(l);
    };
  };
};

void SNRZSystem::RenewFluxMomentInnerIteration(int sn,int qid,real *flxmom,real *anflx)
{
  vector<real> moment(sn*plnum);
  int index=0;
  for(int j=0;j<sn;j++){
    for(int l=0;l<plnum;l++){
      moment[index++]=quad[qid]->GetMoment(l,j);
    };
  };

  real *finp=new real[plnum];
  int ind=0;
  int ind2=0;
  for(int i=0;i<TotM;i++){
    for(int j=0;j<plnum;j++){finp[j]=0.;};
    int index=0;
    for(int j=0;j<sn;j++){
      real tmp=anflx[ind]*quad[qid]->GetOmega(j)*PI2;
      ind++;
      //finp[0]+=tmp;
      //for(int l=1;l<plnum;l++){
      for(int l=0;l<plnum;l++){
	//finp[l]+=tmp*quad[qid]->GetMoment(l,j);
	finp[l]+=tmp*moment[index++];
      };
    };
    for(int l=0;l<plnum;l++){
      flxmom[ind2]=real(finp[l]);
      ind2++;
    };
  };
  delete [] finp;
};

real SNRZSystem::PutFluxAfterInnerIteration(int g,real *flxmom)
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
    if(err>errf){
      errf=err;
    };
  };

  return errf;
};

real SNRZSystem::CalFluxGeneral(int ginp,real cin,int iter)
{
  return CalFlux(ginp,iter,cin);
  //return CalFluxFirstCollision(ginp,iter,cin);
};

real SNRZSystem::CalFlux(int g,int oiter,real epsif)
{
  bool cmfd_on=false;
  if(cmfd&&oiter%opt.GetItcmfd()==0)cmfd_on=true;

  int qid = quadid[g];
  int sn  = quad[qid]->GetSnnum();
  int sn_z = quad[qid]->GetSn(); // o
  int sn_z_div2=sn_z/2;
 
  real *Src    = new real[TotM*plnum]; // External source in pl
  real *SSrc   = new real[TotM]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  = new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM]; // Sigma_t * Volume
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  real *sigt = new real[TotM]; // o
  vector<real> vol_inv(TotM);

  // o
  vector<real> srr(TotM);
  vector<real> sll(TotM);
  vector<real> sbb(TotM);
  for(int i=0;i<TotM;i++){
    srr[i]=mesh[i].GetSurR(0);
    sll[i]=mesh[i].GetSurL(0);
    sbb[i]=mesh[i].GetSurR(1);
  };
  // o

  // Put sigma_t and sigma_s to temporal array
  PutSigmaForInnerIteration(g,sigtvol,sigsself,sigt);
  for(int i=0;i<TotM;i++){
    vol_inv[i]=sigt[i]/sigtvol[i];
  };

  int xm=mi.GetXF();
  int ym=mi.GetYF();

  real flx;
  vector<real> fly(xm);
  vector<real> bflx(ym*sn);
  vector<real> bfly(xm*sn);
  vector<real> bfly3(xm*sn); 
  for(int i=0;i<ym*sn;i++){bflx[i]=0.;};
  for(int i=0;i<xm*sn;i++){bfly[i]=0.;};
  for(int i=0;i<xm*sn;i++){bfly3[i]=0.;};

  // Initial Flux guess
  //chibatmp
  InitialFluxGuessInnerIteration(g,flxmom);

  // Put fission & slowing down source in temporal array
  PutSourceInnerIteration(Src);

  real *flangle=new real[TotM];
  real *SSrcSP =new real[TotM];
  real *bflys=new real[TotM*sn_z_div2];
  real *bflys3=new real[TotM*sn_z_div2];

  real tt1[]={1.,0.};
  real tt2[]={1.,0.5};
  real tt3[]={1.,0.};
  real tt4[]={0.,1.};
  real tt5[]={0.,1.};
  real tt6[]={1.,0.};

  int itmax=itmax_inner;
  bool InnerConvergence=true;

  for(int it=0;it<itmax;it++){
    //cout<<"ITERATION : "<<it<<"\n";
    // for leakage calculation
    //leakage.put_data(g,0.);
    // for CMFD
    if(cmfd_on){
      for(int i=0;i<TotM;i++){
        curxn[i]=0.;
        curxp[i]=0.;
        curyn[i]=0.;
        curyp[i]=0.;
        //curzn[i]=0.;
        //curzp[i]=0.;
      };
    };

    // chibatmp
    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    int is=0;
    for(int n=0;n<sn_z;n++){
      int sizex=quad[qid]->GetSizex(n); 
      // ********************************************
      // Starting direction method
      real xi=quad[qid]->GetXi(is);
      real mu=-sqrt(1.-xi*xi);
      if(!opt.Forward()){
        xi*=-1.;
	      mu*=-1.;
      };
      CalFluxInSourceSP(Src,SSrcSP,flxmom,sigsself,qid,xi,mu);
      real muu=fabs(mu)*2.;
      real xii=fabs(xi)*2.;
      if(n>=sn_z_div2&&BC[2]==1){
	      int tmp=(sn_z-n-1)*xm;
	      for(int i=0;i<xm;i++){
	        fly[i]=bflys[tmp+i];
	        };
      }else if(n<sn_z_div2&&BC[3]==1){
	      int tmp=n*xm;
	      for(int i=0;i<xm;i++){
	        fly[i]=bflys3[tmp+i];
	      };
      }else{
        for(int i=0;i<xm;i++){fly[i]=0.;};
      };
      for(int yy=0;yy<ym;yy++){
	      int y=yy;
	      if(n<sn_z_div2)y=ym-yy-1;
	      flx=0.;
	      real yl=1./mi.GetFMeshL(1,y);
	      //cout<<n<<" "<<y<<"/"<<ym<<"\n";
	      //cout<<edge[0][y][0][1]<<"\n";
	      int xer=edge[0][y][0][1];
        int xnum=xer+1;
	      //cout<<n<<" "<<y<<"/"<<ym<<" "<<xer<<" "<<xnum<<"\n";
	      //for(int x=xm-1;x>=0;x--){
	      for(int x=xnum-1;x>=0;x--){	  
	        real xl=1./mi.GetFMeshL(0,x);
	        //int m=y*xm+x;
	        int m=meshid[0][y][x];
	        //cout<<x<<" "<<m<<"\n";
	        real ss=SSrcSP[m]*vol_inv[m];
	        //sweep
	        real inx=flx;
	        real iny=fly[x];
	        real cfl,outx,outy;
	        if(step_d){
            // Step differencing
	          cfl=(fabs(mu)*xl*inx+fabs(xi)*yl*iny*ss)/(fabs(mu)*xl+fabs(xi)*yl+sigt[m]);
	          outx=cfl;
	          outy=cfl;
	        }else{
	          // Diamond differencing
	          cfl=(muu*xl*inx+xii*yl*iny+ss)/(muu*xl+xii*yl+sigt[m]);
	          real cfl2=cfl*2.;
	          outx=cfl2-inx;
	          outy=cfl2-iny;
	          while(outx<0.||outy<0.){
	            real t1=1.;
	            real t2=1.;
	            real t3=1.;
	            real t4=1.;
	            if(outx<=0.){t1=0.5; t3=0.;};
	            if(outy<=0.){t2=0.5; t4=0.;};
	            cfl=(t1*muu*xl*inx+t2*xii*yl*iny+ss)/(t3*muu*xl+t4*xii*yl+sigt[m]);
	            cfl2=cfl*2.; 
	            if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	            if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	          };
	        };
	        flangle[m]=cfl;
	        flx=outx;
	        fly[x]=outy;
	      };
	      if(BC[2]==1&&n<sn_z_div2){
	        int tmp=n*xm;
	        for(int i=0;i<xm;i++){
	          bflys[tmp+i]=fly[i];
	        };
	      }else if(BC[3]==1&&n>=sn_z_div2){
	        int tmp=(sn_z-n-1)*xm;
	        for(int i=0;i<xm;i++){
	          bflys3[tmp+i]=fly[i];
	        };
	      };
      };
      // end of starting direction method
      // *********************************************
      for(int nn=0;nn<sizex;nn++){
	      int iss=is;
	      if(!opt.Forward())iss=quad[qid]->GetXYref(is);
        CalFluxInSource(Src,SSrc,flxmom,sigsself,qid,iss);
        real mu=quad[qid]->GetMu(is);
        real xi=quad[qid]->GetXi(is);
	      real om=quad[qid]->GetOmega(is);
        real muw=mu*om*PI2;
        real xiw=xi*om*PI2;
	      real muu=fabs(mu);
	      real xii=fabs(xi)*2.;
        real om_inv=1./om;
	      real aneg=quad[qid]->GetMC(is);
	      real apos=aneg-mu*om;
	      int xrsn=quad[qid]->GetXref(is);
	      int yrsn=quad[qid]->GetYref(is);
	      if(xi>0.&&BC[2]==1){
          int tmp=yrsn*xm;
	        for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	      }else if(xi<0.&&BC[3]==1){
	        int tmp=yrsn*xm;
	        for(int i=0;i<xm;i++){fly[i]=bfly3[tmp+i];};
	      }else{
	        for(int i=0;i<xm;i++){fly[i]=0.;};
	      };
	      for(int yy=0;yy<ym;yy++){
	        int y=yy;
	        if(xi<0.)y=ym-yy-1;
	        if(mu>0.){flx=bflx[xrsn*ym+y];}else{flx=0.;};  	  
	        int yxm=y*xm;
	        int xer=edge[0][y][0][1];
          int xnum=xer+1;
  	      //for(int xx=0;xx<xm;xx++){
  	      for(int xx=0;xx<xnum;xx++){	    
	          int x=xx;
	          //if(mu<0.)x=xm-xx-1;
	          if(mu<0.)x=xer-xx;	    
	          //int m=yxm+x;
	          int m=meshid[0][y][x];
	          int index=m*sn+iss;
	          real ss=SSrc[m];
	          real sr=srr[m];
	          real sl=sll[m];
	          real sb=sbb[m];
  	        //sweep
	          real inx=flx;
	          real iny=fly[x];
	          real inm=flangle[m];
	          real stv=sigtvol[m];
	          // (for angular-dependent total)
	          //if(!dsa)stv=aflux[m][g].get_dat(is)*mesh[m].GetVolume();
	          real cfl,outx,outy,outm;
            // The following is modified for the BSA calculation in 2021/2/16.
	          // In the calculations for a system including void region,
	          // negative flux fix-up results in zero angular neutron flux.
	          // To avoid this problem, the step differencing scheme is adopted
	          // if a product of cross section and volume is small.
	          // The threshold is tentatively determined as [1e-6] 
	          //if(step_d){
	          if(step_d||stv<1e-6){
	            // ... Step differencing
	            real rhs=fabs(xi)*sb*iny+aneg*om_inv*(sr-sl)*inm+ss;
  	          if(mu>0.){
	              rhs+=mu*sl*inx;
	            }else{
	            rhs+=-mu*sr*inx;
	            };
	            real lhs=fabs(xi)*sb+apos*om_inv*(sr-sl)+stv;
	            if(mu>0.){
	              lhs+=mu*sr;
	            }else{
	              lhs+=-mu*sl;
	            };
	            cfl=rhs/lhs;
	            outx=cfl;
	            outy=cfl;
	            outm=cfl;
	          }else{
	            // ... Diamond differencing
	            real c1=muu*(sr+sl);
	            real c2=xii*sb;
	            real coef=om_inv*(sr-sl);
	            real c3=coef*(apos+aneg);
	            real rhs=c1*inx+c2*iny+c3*inm+ss;
	            real lhs=c1+c2+c3+stv;
	            cfl=rhs/lhs;
	            real cfl2=cfl*2.;
	            outx=cfl2-inx;
	            outy=cfl2-iny;
	            outm=cfl2-inm;
	            while(outx<0.||outy<0.||outm<0.){
	              int tx=0;
	              int ty=0;
	              int tm=0;
	              if(outx<=0.)tx=1;
	              if(outy<=0.)ty=1;
	              if(outm<=0.)tm=1;
	              real coef2;
	              if(mu<0.){coef2=sl*muu;}else{coef2=sr*muu;};
	              real cc1=c1-tt5[tx]*coef2;
	              real cc3=c3-tt4[tm]*coef*apos;
	              real rhs=cc1*inx+tt2[ty]*c2*iny+cc3*inm+ss;
	              real c4=tt1[tx]*2.*coef2;
	              real c6=tt6[tm]*2.*coef*apos;
	              real lhs=c4+tt3[ty]*c2+c6+stv;
	              cfl=rhs/lhs;
	              cfl2=cfl*2.;
                if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	              if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	              if(outm>0.){outm=cfl2-inm;}else{outm=0.;};
	            };
	          };
	          anflx[index]=cfl;
	          flx=outx;
	          fly[x]=outy;
	          flangle[m]=outm;
	          // for leakage calculation
	          /*
	          if(mu>0.&&xx==xm-1&&BC[1]!=1){
	            leakage.add_data(g,outx*muw*sr);
	          };
	          if(xi>0.&&yy==ym-1&&BC[3]!=1){
	            leakage.add_data(g,outy*xiw*sb);
	          };
	          if(xi<0.&&yy==0&&BC[2]!=1){
	            leakage.add_data(g,-outy*xiw*sb);
	          };
	          */
	          // for CMFD
	          if(cmfd_on){
	            if(mu<0.){curxn[m]+=outx*muw;}
	            else{curxp[m]+=outx*muw;};
	            //if(et<0.){curyn[m]+=outy*etw;}
	            //else{curyp[m]+=outy*etw;};
	            if(xi<0.){curyn[m]+=outy*xiw;}
	            else{curyp[m]+=outy*xiw;};
	          };
      	  };
	        if(mu<0.&&BC[0]==1){bflx[is*ym+y]=flx;};
	      };
	      if(xi<0.&&BC[2]==1){
	        int tmp=is*xm;
	        for(int i=0;i<xm;i++){
	          bfly[tmp+i]=fly[i];
	        };
        }else if(xi>0.&&BC[3]==1){
	        int tmp=is*xm;
	        for(int i=0;i<xm;i++){
	          bfly3[tmp+i]=fly[i];
	        };
	      };
	      is++;
      };
    };

    // Renew flux moment
    RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);
    if(it<30&&dsa&&itmax!=1)AccelerationByDSA(g,flxmom,flxd,sigsself,cmfd_on);
    real errmax=0.;
    int ind=0;
    // chibatmp
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    //cout<<" "<<errmax;
    //if(g==69)cout<<" . "<<errmax;
    if(errmax<epsif){
      //cout<<"   group "<<g<<" : "<<it<<" ("<<epsif<<")\n";
      break;
    };
    if(it==itmax-1){
      cout<<"#   ... Not converged in group "<<g<<" within "<<itmax<<" inner iterations.\n";
      cout<<"#       (maximum error is : "<<errmax<<")\n";
      InnerConvergence=false;
    };
  };

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  int ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]);
      // (for angular dependent total)
      //if(dsa&&wrtflx[m])aflux[m][g].put_data(is,anflx[ind]);
      //if(wrtflx[m]&&oiter==20)aflux[m][g].put_data(is,anflx[ind]);
      ind++;
    };
  };

  // for CMFD
  if(cmfd_on){
    if((g==0&&opt.Forward())||(g==grp-1&&!opt.Forward()))SetZeroCurFF();
    CalCoarseCur();
  };

  delete [] flangle;
  delete [] SSrcSP;
  delete [] bflys;
  delete [] bflys3;

  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxmom;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigt;
  delete [] sigsself;

  return errf;
};

real SNRZSystem::CalFluxAllAbsorption(int g,int oiter,real epsif)
// neutron source is assumed to be isotropic
{
  int qid = quadid[g];
  int sn  = quad[qid]->GetSnnum();
  int sn_z = quad[qid]->GetSn(); // o
  int sn_z_div2=sn_z/2;
 
  real *Src    = new real[TotM*plnum]; // External source in pl
  real *SSrc   = new real[TotM]; // Self-scattering source
  real *anflx  = new real[TotM*sn]; // Angular flux
  real *flxmom  = new real[TotM*plnum]; // Flux moment
  real *sigtvol = new real[TotM]; // Sigma_t * Volume
  real *flxd   = new real[TotM];

  real *sigt = new real[TotM]; // o
  vector<real> vol_inv(TotM);

  // o
  vector<real> srr(TotM);
  vector<real> sll(TotM);
  vector<real> sbb(TotM);
  for(int i=0;i<TotM;i++){
    srr[i]=mesh[i].GetSurR(0);
    sll[i]=mesh[i].GetSurL(0);
    sbb[i]=mesh[i].GetSurR(1);
  };
  // o

  // Put sigma_t and sigma_s to temporal array
  int ind=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    real sigttmp=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
    sigt[i]=sigttmp;
    sigtvol[i]=sigttmp*vol;
    vol_inv[i]=1./vol;
  };

  int xm=mi.GetXF();
  int ym=mi.GetYF();

  real flx;
  vector<real> fly(xm);
  vector<real> bflx(ym*sn);
  vector<real> bfly(xm*sn);
  for(int i=0;i<ym*sn;i++){bflx[i]=0.;};
  for(int i=0;i<xm*sn;i++){bfly[i]=0.;};

  // Put fission & slowing down source in temporal array
  PutSourceInnerIteration(Src);

  real *flangle=new real[TotM];
  real *SSrcSP =new real[TotM];
  real *bflys=new real[TotM*sn_z_div2];

  real tt1[]={1.,0.};
  real tt2[]={1.,0.5};
  real tt3[]={1.,0.};
  real tt4[]={0.,1.};
  real tt5[]={0.,1.};
  real tt6[]={1.,0.};

    int is=0;
    for(int n=0;n<sn_z;n++){

      int sizex=quad[qid]->GetSizex(n); 
      // ********************************************
      // Starting direction method

      real xi=quad[qid]->GetXi(is);
      real mu=-sqrt(1.-xi*xi);
      if(!opt.Forward()){
        xi*=-1.;
	mu*=-1.;
      };

      //CalFluxInSourceSP(Src,SSrcSP,flxmom,sigsself,qid,xi,mu);
      int ind2=0;
      for(int i=0;i<TotM;i++){
        //real srcin=0.;
        //for(int l=0;l<=pl;l++){
        //  for(int m=0;m<=l;m++){ // for 2D
        //    srcin+=(Src[ind2])*quad[qid]->GetMoment(l,m,mu,xi);
	//    ind2++;
        //  };
        //};
        //if(srcin<0.)srcin=0.;
        //SSrc[i]=srcin; 
        SSrc[i]=Src[ind2];
        ind2+=plnum;
      };

      real muu=fabs(mu)*2.;
      real xii=fabs(xi)*2.;
      if(n>=sn_z_div2&&BC[2]==1){
	int tmp=(sn_z-n-1)*xm;
	for(int i=0;i<xm;i++){
	  fly[i]=bflys[tmp+i];
	};
      }else{
        for(int i=0;i<xm;i++){fly[i]=0.;};
      };
      for(int yy=0;yy<ym;yy++){
	int y=yy;
	if(n<sn_z_div2)y=ym-yy-1;
	flx=0.;
	real yl=1./mi.GetFMeshL(1,y);
	for(int x=xm-1;x>=0;x--){
	  real xl=1./mi.GetFMeshL(0,x);
	  int m=y*xm+x;
	  real ss=SSrcSP[m]*vol_inv[m];
	  //sweep
	  real inx=flx;
	  real iny=fly[x];
	  real cfl,outx,outy;
	  if(step_d){
          // Step differencing
	  cfl=(fabs(mu)*xl*inx+fabs(xi)*yl*iny+ss)/(fabs(mu)*xl+fabs(xi)*yl+sigt[m]);
	  outx=cfl;
	  outy=cfl;
	  }else{
	  // Diamond differencing
	  cfl=(muu*xl*inx+xii*yl*iny+ss)/(muu*xl+xii*yl+sigt[m]);
	  real cfl2=cfl*2.;
	  outx=cfl2-inx;
	  outy=cfl2-iny;
	  while(outx<0.||outy<0.){
	    real t1=1.;
	    real t2=1.;
	    real t3=1.;
	    real t4=1.;
	    if(outx<=0.){t1=0.5; t3=0.;};
	    if(outy<=0.){t2=0.5; t4=0.;};
	    cfl=(t1*muu*xl*inx+t2*xii*yl*iny+ss)/(t3*muu*xl+t4*xii*yl+sigt[m]);
	    cfl2=cfl*2.; 
	    if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	    if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	  };
	  };

	  flangle[m]=cfl;
	  flx=outx;
	  fly[x]=outy;
	};
	if(n<sn_z_div2){
	  int tmp=n*xm;
	  for(int i=0;i<xm;i++){
	    bflys[tmp+i]=fly[i];
	  };
	};
      };
      // end of starting direction method
      // *********************************************
 
      for(int nn=0;nn<sizex;nn++){

	int iss=is;
	if(!opt.Forward())iss=quad[qid]->GetXYref(is);

        int ind2=0;
        for(int i=0;i<TotM;i++){
          //real tmp=0.;
          //int index=0;
          //for(int l=0;l<=pl;l++){
          //  for(int m=0;m<=l;m++){ // for 2D
          //    tmp+=Src[ind2++]*quad[qid]->GetMoment(index,iss);  
	  //    index++;
          //  };
          //};
          //SSrc[i]=tmp;
	  SSrc[i]=Src[ind2];
	  ind2+=plnum;
        };

        real mu=quad[qid]->GetMu(is);
        real xi=quad[qid]->GetXi(is);
	real muu=fabs(mu);
	real xii=fabs(xi)*2.;
	real om=quad[qid]->GetOmega(is);
        real om_inv=1./om;
	real aneg=quad[qid]->GetMC(is);
	real apos=aneg-mu*om;
	int xrsn=quad[qid]->GetXref(is);
	int yrsn=quad[qid]->GetYref(is);
	if(xi>0.&&BC[2]==1){
          int tmp=yrsn*xm;
	  for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	}else{
	  for(int i=0;i<xm;i++){fly[i]=0.;};
	};
	for(int yy=0;yy<ym;yy++){
	  int y=yy;
	  if(xi<0.)y=ym-yy-1;
	  if(mu>0.){flx=bflx[xrsn*ym+y];}else{flx=0.;};  	  
	  int yxm=y*xm;
  	  for(int xx=0;xx<xm;xx++){
	    int x=xx;
	    if(mu<0.)x=xm-xx-1;
	    int m=yxm+x;
	    int index=m*sn+iss;
	    real ss=SSrc[m];
	    real sr=srr[m];
	    real sl=sll[m];
	    real sb=sbb[m];
  	    //sweep
	    real inx=flx;
	    real iny=fly[x];
	    real inm=flangle[m];
	    real stv=sigtvol[m];
	    real cfl,outx,outy,outm;
	    if(step_d){
	    // Step differencing
            real rhs=fabs(xi)*sb*iny+aneg*om_inv*(sr-sl)*inm+ss;
	    if(mu>0.){
	      rhs+=mu*sl*inx;
	    }else{
	      rhs+=-mu*sr*inx;
	    };
	    real lhs=fabs(xi)*sb+apos*om_inv*(sr-sl)+stv;
	    if(mu>0.){
	      lhs+=mu*sr;
	    }else{
	      lhs+=-mu*sl;
	    };
	    cfl=rhs/lhs;
	    outx=cfl;
	    outy=cfl;
	    outm=cfl;
	    }else{
	    // Diamond differencing
	    real c1=muu*(sr+sl);
	    real c2=xii*sb;
	    real coef=om_inv*(sr-sl);
	    real c3=coef*(apos+aneg);
	    real rhs=c1*inx+c2*iny+c3*inm+ss;
	    real lhs=c1+c2+c3+stv;

	    cfl=rhs/lhs;
	    real cfl2=cfl*2.;
	    outx=cfl2-inx;
	    outy=cfl2-iny;
	    outm=cfl2-inm;

	    while(outx<0.||outy<0.||outm<0.){
	      int tx=0;
	      int ty=0;
	      int tm=0;
	      if(outx<=0.)tx=1;
	      if(outy<=0.)ty=1;
	      if(outm<=0.)tm=1;
	      real coef2;
	      if(mu<0.){coef2=sl*muu;}else{coef2=sr*muu;};
	      real cc1=c1-tt5[tx]*coef2;
	      real cc3=c3-tt4[tm]*coef*apos;
	      real rhs=cc1*inx+tt2[ty]*c2*iny+cc3*inm+ss;
	      real c4=tt1[tx]*2.*coef2;
	      real c6=tt6[tm]*2.*coef*apos;
	      real lhs=c4+tt3[ty]*c2+c6+stv;
	      cfl=rhs/lhs;
	      cfl2=cfl*2.;
              if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	      if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	      if(outm>0.){outm=cfl2-inm;}else{outm=0.;};
	    };
	    };

	    //

	    anflx[index]=cfl;
	    flx=outx;
	    fly[x]=outy;
	    flangle[m]=outm;

	  };
	  if(mu<0.&&BC[0]==1){bflx[is*ym+y]=flx;};
	};
	if(xi<0.&&BC[2]==1){
	  int tmp=is*xm;
	  for(int i=0;i<xm;i++){
	    bfly[tmp+i]=fly[i];
	  };
	};
	is++;
      };
    };

    // Renew flux moment
    RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);

  // Write angular flux for perturbation calculation **
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m])aflux[m][g].put_data(is,anflx[ind]);
      ind++;
    };
  };

  delete [] flangle;
  delete [] SSrcSP;
  delete [] bflys;

  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] flxmom;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigt;

  return errf;
};

real SNRZSystem::CalFluxFirstCollision(int g,int oiter,real epsif)
{
  int qid = quadid[g];
  int sn  = quad[qid]->GetSnnum();
  int sn_z = quad[qid]->GetSn(); // o
  int sn_z_div2=sn_z/2;
 
  real *Src    = new real[TotM*plnum];      // External source in pl
  real *SSrc   = new real[TotM];            // Self-scattering source
  real *anflx  = new real[TotM*sn];         // Angular flux
  real *anflx_unc = new real[TotM*sn];      // Angular flux (Uncollided)
  real *flxmom  = new real[TotM*plnum];     // Flux moment
  real *flxmom_unc  = new real[TotM*plnum]; // Flux moment (Uncollided)
  for(int i=0;i<TotM*plnum;i++){
    flxmom_unc[i]=0.;
  };
  real *sigtvol = new real[TotM]; // Sigma_t * Volume
  real *sigsself = new real[TotM*(pl+1)]; 
  real *flxd   = new real[TotM];

  real *sigt = new real[TotM]; // o
  vector<real> vol_inv(TotM);

  // o
  vector<real> srr(TotM);
  vector<real> sll(TotM);
  vector<real> sbb(TotM);
  for(int i=0;i<TotM;i++){
    srr[i]=mesh[i].GetSurR(0);
    sll[i]=mesh[i].GetSurL(0);
    sbb[i]=mesh[i].GetSurR(1);
  };
  // o

  // Put sigma_t and sigma_s to temporal array
  PutSigmaForInnerIteration(g,sigtvol,sigsself,sigt);
  for(int i=0;i<TotM;i++){
    vol_inv[i]=sigt[i]/sigtvol[i];
  };

  int xm=mi.GetXF();
  int ym=mi.GetYF();

  real flx;
  vector<real> fly(xm);
  vector<real> bflx(ym*sn);
  vector<real> bfly(xm*sn);
  for(int i=0;i<ym*sn;i++){bflx[i]=0.;};
  for(int i=0;i<xm*sn;i++){bfly[i]=0.;};

  // ++++++++++++++++++++++++++++
  // +++ First collision part +++
  // ++++++++++++++++++++++++++++
  //
  // Put fission & slowing down source in temporal array
  PutSourceInnerIteration(Src);

  real *flangle=new real[TotM];
  real *SSrcSP =new real[TotM];
  real *bflys=new real[TotM*sn_z_div2];

  real tt1[]={1.,0.};
  real tt2[]={1.,0.5};
  real tt3[]={1.,0.};
  real tt4[]={0.,1.};
  real tt5[]={0.,1.};
  real tt6[]={1.,0.};

    int is=0;
    for(int n=0;n<sn_z;n++){

      int sizex=quad[qid]->GetSizex(n); 
      // ********************************************
      // Starting direction method

      real xi=quad[qid]->GetXi(is);
      real mu=-sqrt(1.-xi*xi);
      if(!opt.Forward()){
        xi*=-1.;
	mu*=-1.;
      };
      CalFluxInSourceSP(Src,SSrcSP,flxmom_unc,sigsself,qid,xi,mu);
      real muu=fabs(mu)*2.;
      real xii=fabs(xi)*2.;
      if(n>=sn_z_div2&&BC[2]==1){
	int tmp=(sn_z-n-1)*xm;
	for(int i=0;i<xm;i++){
	  fly[i]=bflys[tmp+i];
	};
      }else{
        for(int i=0;i<xm;i++){fly[i]=0.;};
      };
      for(int yy=0;yy<ym;yy++){
	int y=yy;
	if(n<sn_z_div2)y=ym-yy-1;
	flx=0.;
	real yl=1./mi.GetFMeshL(1,y);
	for(int x=xm-1;x>=0;x--){
	  real xl=1./mi.GetFMeshL(0,x);
	  int m=y*xm+x;
	  real ss=SSrcSP[m]*vol_inv[m];
	  //sweep
	  real inx=flx;
	  real iny=fly[x];
	  real cfl,outx,outy;
	  if(step_d){
          // Step differencing
	  cfl=(fabs(mu)*xl*inx+fabs(xi)*yl*iny*ss)/(fabs(mu)*xl+fabs(xi)*yl+sigt[m]);
	  outx=cfl;
	  outy=cfl;
	  }else{
	  // Diamond differencing
	  cfl=(muu*xl*inx+xii*yl*iny+ss)/(muu*xl+xii*yl+sigt[m]);
	  real cfl2=cfl*2.;
	  outx=cfl2-inx;
	  outy=cfl2-iny;
	  while(outx<0.||outy<0.){
	    real t1=1.;
	    real t2=1.;
	    real t3=1.;
	    real t4=1.;
	    if(outx<=0.){t1=0.5; t3=0.;};
	    if(outy<=0.){t2=0.5; t4=0.;};
	    cfl=(t1*muu*xl*inx+t2*xii*yl*iny+ss)/(t3*muu*xl+t4*xii*yl+sigt[m]);
	    cfl2=cfl*2.; 
	    if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	    if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	  };
	  };
	  flangle[m]=cfl;
	  flx=outx;
	  fly[x]=outy;
	};
	if(n<sn_z_div2){
	  int tmp=n*xm;
	  for(int i=0;i<xm;i++){
	    bflys[tmp+i]=fly[i];
	  };
	};
      };
      // end of starting direction method
      // *********************************************
 
      for(int nn=0;nn<sizex;nn++){

	int iss=is;
	if(!opt.Forward())iss=quad[qid]->GetXYref(is);

        CalFluxInSource(Src,SSrc,flxmom_unc,sigsself,qid,iss);

        real mu=quad[qid]->GetMu(is);
        real xi=quad[qid]->GetXi(is);
	real muu=fabs(mu);
	real xii=fabs(xi)*2.;
	real om=quad[qid]->GetOmega(is);
        real om_inv=1./om;
	real aneg=quad[qid]->GetMC(is);
	real apos=aneg-mu*om;
	int xrsn=quad[qid]->GetXref(is);
	int yrsn=quad[qid]->GetYref(is);
	if(xi>0.&&BC[2]==1){
          int tmp=yrsn*xm;
	  for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	}else{
	  for(int i=0;i<xm;i++){fly[i]=0.;};
	};
	for(int yy=0;yy<ym;yy++){
	  int y=yy;
	  if(xi<0.)y=ym-yy-1;
	  if(mu>0.){flx=bflx[xrsn*ym+y];}else{flx=0.;};  	  
	  int yxm=y*xm;
  	  for(int xx=0;xx<xm;xx++){
	    int x=xx;
	    if(mu<0.)x=xm-xx-1;
	    int m=yxm+x;
	    int index=m*sn+iss;
	    real ss=SSrc[m];
	    real sr=srr[m];
	    real sl=sll[m];
	    real sb=sbb[m];
  	    //sweep
	    real inx=flx;
	    real iny=fly[x];
	    real inm=flangle[m];
	    real stv=sigtvol[m];

	    real c1=muu*(sr+sl);
	    real c2=xii*sb;
	    real coef=om_inv*(sr-sl);
	    real c3=coef*(apos+aneg);
	    real rhs=c1*inx+c2*iny+c3*inm+ss;
	    real lhs=c1+c2+c3+stv;

	    real cfl=rhs/lhs;
	    real cfl2=cfl*2.;
	    real outx=cfl2-inx;
	    real outy=cfl2-iny;
	    real outm=cfl2-inm;

	    while(outx<0.||outy<0.||outm<0.){
	      int tx=0;
	      int ty=0;
	      int tm=0;
	      if(outx<=0.)tx=1;
	      if(outy<=0.)ty=1;
	      if(outm<=0.)tm=1;
	      real coef2;
	      if(mu<0.){coef2=sl*muu;}else{coef2=sr*muu;};
	      real cc1=c1-tt5[tx]*coef2;
	      real cc3=c3-tt4[tm]*coef*apos;
	      real rhs=cc1*inx+tt2[ty]*c2*iny+cc3*inm+ss;
	      real c4=tt1[tx]*2.*coef2;
	      real c6=tt6[tm]*2.*coef*apos;
	      real lhs=c4+tt3[ty]*c2+c6+stv;
	      cfl=rhs/lhs;
	      cfl2=cfl*2.;
              if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	      if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	      if(outm>0.){outm=cfl2-inm;}else{outm=0.;};
	    };
	    anflx_unc[index]=cfl;
	    flx=outx;
	    fly[x]=outy;
	    flangle[m]=outm;
	  };
	  if(mu<0.&&BC[0]==1){bflx[is*ym+y]=flx;};
	};
	if(xi<0.&&BC[2]==1){
	  int tmp=is*xm;
	  for(int i=0;i<xm;i++){
	    bfly[tmp+i]=fly[i];
	  };
	};
	is++;
      };
    };

    // Uncollided angular flux moment -> flxmom
    RenewFluxMomentInnerIteration(sn,qid,flxmom_unc,anflx_unc);
    int ind=0;
    int ind2=0;
    for(int i=0;i<TotM;i++){
      for(int l=0;l<=pl;l++){
	real sigs=sigsself[ind++];
	for(int m=0;m<=l;m++){
	  Src[ind2]=flxmom_unc[ind2]*sigs;
	  ind2++;
	};
      };
    };

    // +++++++++++++++++++++++++++++++++++++++++++
    // + Calculation with first-collision source +
    // +++++++++++++++++++++++++++++++++++++++++++
  
  for(int i=0;i<ym*sn;i++){bflx[i]=0.;};
  for(int i=0;i<xm*sn;i++){bfly[i]=0.;};

  // Initial Flux guess
  InitialFluxGuessInnerIteration(g,flxmom);

  int itmax=itmax_inner;
  bool InnerConvergence=true;
  for(int it=0;it<itmax;it++){

    int da=0;
    for(int i=0;i<TotM;i++){
      flxd[i]=flxmom[da];
      da+=plnum;
    };

    int is=0;
    for(int n=0;n<sn_z;n++){

      int sizex=quad[qid]->GetSizex(n); 
      // ********************************************
      // Starting direction method

      real xi=quad[qid]->GetXi(is);
      real mu=-sqrt(1.-xi*xi);
      if(!opt.Forward()){
        xi*=-1.;
	mu*=-1.;
      };
      CalFluxInSourceSP(Src,SSrcSP,flxmom,sigsself,qid,xi,mu);
      real muu=fabs(mu)*2.;
      real xii=fabs(xi)*2.;
      if(n>=sn_z_div2&&BC[2]==1){
	int tmp=(sn_z-n-1)*xm;
	for(int i=0;i<xm;i++){
	  fly[i]=bflys[tmp+i];
	};
      }else{
        for(int i=0;i<xm;i++){fly[i]=0.;};
      };
      for(int yy=0;yy<ym;yy++){
	int y=yy;
	if(n<sn_z_div2)y=ym-yy-1;
	flx=0.;
	real yl=1./mi.GetFMeshL(1,y);
	for(int x=xm-1;x>=0;x--){
	  real xl=1./mi.GetFMeshL(0,x);
	  int m=y*xm+x;
	  real ss=SSrcSP[m]*vol_inv[m];
	  //sweep
	  real inx=flx;
	  real iny=fly[x];
	  real cfl=(muu*xl*inx+xii*yl*iny+ss)/(muu*xl+xii*yl+sigt[m]);
	  real cfl2=cfl*2.;
	  real outx=cfl2-inx;
	  real outy=cfl2-iny;
	  while(outx<0.||outy<0.){
	    real t1=1.;
	    real t2=1.;
	    real t3=1.;
	    real t4=1.;
	    if(outx<=0.){t1=0.5; t3=0.;};
	    if(outy<=0.){t2=0.5; t4=0.;};
	    cfl=(t1*muu*xl*inx+t2*xii*yl*iny+ss)/(t3*muu*xl+t4*xii*yl+sigt[m]);
	    cfl2=cfl*2.; 
	    if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	    if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	  };
	  flangle[m]=cfl;
	  flx=outx;
	  fly[x]=outy;
	};
	if(n<sn_z_div2){
	  int tmp=n*xm;
	  for(int i=0;i<xm;i++){
	    bflys[tmp+i]=fly[i];
	  };
	};
      };
      // end of starting direction method
      // *********************************************
 
      for(int nn=0;nn<sizex;nn++){

	int iss=is;
	if(!opt.Forward())iss=quad[qid]->GetXYref(is);

        CalFluxInSource(Src,SSrc,flxmom,sigsself,qid,iss);

        real mu=quad[qid]->GetMu(is);
        real xi=quad[qid]->GetXi(is);
	real muu=fabs(mu);
	real xii=fabs(xi)*2.;
	real om=quad[qid]->GetOmega(is);
        real om_inv=1./om;
	real aneg=quad[qid]->GetMC(is);
	real apos=aneg-mu*om;
	int xrsn=quad[qid]->GetXref(is);
	int yrsn=quad[qid]->GetYref(is);
	if(xi>0.&&BC[2]==1){
          int tmp=yrsn*xm;
	  for(int i=0;i<xm;i++){fly[i]=bfly[tmp+i];};
	}else{
	  for(int i=0;i<xm;i++){fly[i]=0.;};
	};
	for(int yy=0;yy<ym;yy++){
	  int y=yy;
	  if(xi<0.)y=ym-yy-1;
	  if(mu>0.){flx=bflx[xrsn*ym+y];}else{flx=0.;};  	  
	  int yxm=y*xm;
  	  for(int xx=0;xx<xm;xx++){
	    int x=xx;
	    if(mu<0.)x=xm-xx-1;
	    int m=yxm+x;
	    int index=m*sn+iss;
	    real ss=SSrc[m];
	    real sr=srr[m];
	    real sl=sll[m];
	    real sb=sbb[m];
  	    //sweep
	    real inx=flx;
	    real iny=fly[x];
	    real inm=flangle[m];
	    real stv=sigtvol[m];

	    real c1=muu*(sr+sl);
	    real c2=xii*sb;
	    real coef=om_inv*(sr-sl);
	    real c3=coef*(apos+aneg);
	    real rhs=c1*inx+c2*iny+c3*inm+ss;
	    real lhs=c1+c2+c3+stv;

	    real cfl=rhs/lhs;
	    real cfl2=cfl*2.;
	    real outx=cfl2-inx;
	    real outy=cfl2-iny;
	    real outm=cfl2-inm;

	    while(outx<0.||outy<0.||outm<0.){
	      int tx=0;
	      int ty=0;
	      int tm=0;
	      if(outx<=0.)tx=1;
	      if(outy<=0.)ty=1;
	      if(outm<=0.)tm=1;
	      real coef2;
	      if(mu<0.){coef2=sl*muu;}else{coef2=sr*muu;};
	      real cc1=c1-tt5[tx]*coef2;
	      real cc3=c3-tt4[tm]*coef*apos;
	      real rhs=cc1*inx+tt2[ty]*c2*iny+cc3*inm+ss;
	      real c4=tt1[tx]*2.*coef2;
	      real c6=tt6[tm]*2.*coef*apos;
	      real lhs=c4+tt3[ty]*c2+c6+stv;
	      cfl=rhs/lhs;
	      cfl2=cfl*2.;
              if(outx>0.){outx=cfl2-inx;}else{outx=0.;};
	      if(outy>0.){outy=cfl2-iny;}else{outy=0.;};
	      if(outm>0.){outm=cfl2-inm;}else{outm=0.;};
	    };
	    anflx[index]=cfl;
	    flx=outx;
	    fly[x]=outy;
	    flangle[m]=outm;
	  };
	  if(mu<0.&&BC[0]==1){bflx[is*ym+y]=flx;};
	};
	if(xi<0.&&BC[2]==1){
	  int tmp=is*xm;
	  for(int i=0;i<xm;i++){
	    bfly[tmp+i]=fly[i];
	  };
	};
	is++;
      };
    };

    // Renew flux moment
    RenewFluxMomentInnerIteration(sn,qid,flxmom,anflx);

    if(dsa)AccelerationByDSA(g,flxmom,flxd,sigsself,false);

    real errmax=0.;
    int ind=0;
    for(int i=0;i<TotM;i++){
      real err=fabs(flxmom[ind]/flxd[i]-1.);
      ind+=plnum;
      if(err>errmax)errmax=err;
    };
    if(errmax<epsif){
      //cout<<"   group "<<g<<" : "<<it<<" ("<<epsif<<")\n";
      break;
    };
    if(it==itmax-1){
      cout<<"   ... Not converged in group "<<g<<"\n";
      InnerConvergence=false;
    };
  };

  for(int i=0;i<TotM*plnum;i++){
    flxmom[i]+=flxmom_unc[i];
  };

  // Put new flux into GeneralMesh
  real errf=PutFluxAfterInnerIteration(g,flxmom);
  if(!InnerConvergence)errf=-errf;

  // Write angular flux for perturbation calculation **
  ind=0;
  for(int m=0;m<TotM;m++){
    for(int is=0;is<sn;is++){
      if(wrtflx[m]){
        aflux[m][g].put_data(is,anflx[ind]+anflx_unc[ind]);
	//cout<<m<<" ";
      };
      ind++;
    };
  };

  delete [] flangle;
  delete [] SSrcSP;
  delete [] bflys;

  delete [] SSrc;
  delete [] Src;
  delete [] anflx;
  delete [] anflx_unc;
  delete [] flxmom;
  delete [] flxmom_unc;
  delete [] flxd;
  delete [] sigtvol;
  delete [] sigt;
  delete [] sigsself;

  return errf;
};

void SNRZSystem::AccelerationByDSA(int g,real *flxmom,real *flxd,real *sigsself,bool cmfd_on)
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
  /*
  // *** positive source term
  for(int i=0;i<TotM;i++){
    real tmp=fltt[i];
    if(tmp<0.)tmp=0.;
    psys.GetMesh(i).PutSrcin(tmp);  
  };
  psys.CalFluxAutoConv(g,0.01);
  ind=0;
  for(int i=0;i<TotM;i++){
    flxmom[ind]+=psys.GetMesh(i).GetFlux(0).get_dat(g);
    ind+=plnum;
  };
  */
  // *** negative source term
  //for(int i=0;i<TotM;i++){
    //real tmp=-fltt[i];
    //if(tmp<0)tmp=0.;
    //psys.GetMesh(i).PutSrcin(fltt[i]);  
  //};
  psys.CalFluxAutoConv(g,0.01);
  ind=0;
  for(int i=0;i<TotM;i++){
    //real fold=flxmom[ind];
    //real fnew=fold+psys.GetMesh(i).GetFlux(0).get_dat(g);
    //for(int j=0;j<plnum;j++){
    //  flxmom[ind]*=fnew/fold;
    //  ind++;
    //};
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
        //curzp[i]*=tmp3;
        //curzn[i]*=tmp3;
      };
    };
    ind+=plnum;
  };
  //delete[] fltt;
};

void SNRZSystem::SetInitialFlux()
{
  if(dsa){
  GeneralOption optdif;
  optdif.PutEpsf(5e-2);
  optdif.PutEpsk(5e-4);
  optdif.PutEpss(5e-3);
  optdif.PutOutitermax(30);
  //optdif.PutOutitermax(1);
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

void SNRZSystem::AddMedium(Medium inp)
{
  GeneralSystem::AddMedium(inp);
  Medium inp2=inp;
  /*
  for(int g=0;g<grp;g++){
    real org=inp2.GetMacxs().GetData1d(sigtr).get_dat(g);
    inp2.GetMacxs().GetData1d(d).put_data(g,0.281/org);
  };
  */
  psys.AddMedium(inp2); // specific for snt
  //psys.AddMedium(inp); // specific for snt
}

void SNRZSystem::AddMedium(string mdir,string ss,int plt,bool upscat)
{
  GeneralSystem::AddMedium(mdir,ss,plt,upscat);
  if(dsa)psys.AddMedium(mdir,ss,plt,upscat); // specific for snt
}

void SNRZSystem::AddFlux(SNRZSystem &sec)
{
  for(int i=0;i<TotM;i++){
    GroupData1D tmp=mesh[i].GetFlux()+sec.GetMesh(i).GetFlux();
    mesh[i].GetFlux().copy(tmp);
    for(int g=0;g<grp;g++){
      GroupData1D tmp=aflux[i][g]+sec.GetAFlux(i,g);
      aflux[i][g].copy(tmp);
    };
  };
};


// **********************************************
// * For perturbation calculation               *
// **********************************************

real SNRZSystem::CalReactivity(SNRZSystem *sec,real kunp,real kp,bool pr)
{
  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };
  real ret=CalReactivity(sec,kunp,kp,flag,pr);
  delete [] flag;
  return ret;
};

real SNRZSystem::CalReactivity(SNRZSystem *sec,real kunp,real kp,bool *flag,bool pr)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"Error in SNRZSystem::CalReactivity.\n";
    cout<<"Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"Error in SNRZSystem::CalReactivity.\n";
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
    cout<<"Error in SNRZsystem::CalSensitivity.\n";
    cout<<"Perturbation denominator is zero.\n";
    exit(0);
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
    cout<<"#   Component- and group-wise reactivity per unit lathergy x 0.25\n";
    cout<<"#\n";
    cout<<"#      - [HO-scat] means high-order scattering components.\n";
    cout<<"#\n";
    cout<<"#   Energy [eV]   ";
    cout<<"Yield        ";
    cout<<"Absorption   ";
    cout<<"Scattering   ";
    cout<<"(P0-Scat.)   ";
    cout<<"Leakage      ";
    cout<<"(+HO-scat*)  ";
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
      cout<<i<<" ";
      cout.setf(ios::scientific);
      cout.precision(5);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  ";
      cout<<sct[0][i]/leth<<"  "; // P0 scattering
      cout<<le/leth<<"  "; // leakage
      cout<<(le+(sl-sct[0][i]))/leth<<"  "; // leakage + higher-order scattering
      cout<<(yld[i]+abs[i]+sl+le)/leth<<"  "; // total
      cout<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# +++ Summary (energy group-integrated value) +++\n";
    cout<<"#\n";
    cout.setf(ios::scientific);
    cout.precision(5);
    cout<<"# Yield       : "<<yldsum<<"\n";
    cout<<"# Absorption  : "<<abssum<<"\n";
    cout<<"# Scattering  : "<<tsctsum<<"\n";    
    real tmp=yldsum+abssum;
    int id=0;
    real sctpnsum=0.;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=l;m++){ // snrz
        cout<<"#    ("<<l<<","<<m<<")th Comp. :"<<sctsum[id]<<"\n";
	//	if(l==0){tmp+=sctsum[id];}else{tleaksum+=sctsum[id];};
	tmp+=sctsum[id];
	if(l!=0)sctpnsum+=sctsum[id];
        id++;
      };
    };
    cout<<"# (Non-Leak)  : "<<tmp<<"\n";
    cout<<"#\n";

    id=0;
    for(int l=1;l<=pl;l++){
      for(int m=0;m<=l;m++){ // snrz
        cout<<"# ("<<l<<","<<m<<")th leakage   :"<<leaksum[id]<<"\n";
        id++;
      };
    };
    cout<<"# Higher-order   : "<<leaksum[plnum-1]<<"\n";    
    cout<<"# (Leakage)      : "<<tleaksum<<"\n";    

    cout<<"#\n#\n";

    cout<<"# (For alternative categorization)\n";
    cout<<"#\n";
    cout<<"#   [Non-Leakage] - [high-order scattering] : "<<tmp-sctpnsum<<"\n";
    cout<<"#   [Leakage] + [high-order scattering]     : "<<tleaksum+sctpnsum<<"\n";    
  };

  cout.setf(ios::scientific);
  cout.precision(5);

  real tot=yldsum+abssum+tsctsum+tleaksum;
  if(pr){
    cout<<"# \n";    
    cout<<"# ** Perturbation Cal.  : "<<tot<<"\n";
    cout<<"# ** Direct Cal.        : "<<1/kunp-1/kp<<"\n";
    cout<<"#\n";
    cout<<"#  ( Perturbation denominator : "<<ip<<" )\n";
  };
  cout.unsetf(ios::scientific);

  delete [] sctsum;
  delete [] leaksum;
  delete [] yld;
  delete [] abs;

  return tot;
};

void SNRZSystem::CalPerturbLeakScat(SNRZSystem *sec,bool *flag,int g,real *sct,real *leak)
{
  real invpi4=1./PI4;

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
      real dsct=dtot-dabs;
      real fact=PI2*PI4*dtot*vol;
      // total 
      if(dsct!=0.){
        sct[0]-=dsct*sec->GetMesh(m).GetFlux(0).get_dat(g)*
                             mesh[m].GetFlux(0).get_dat(g)*vol;
      };
      if(dtot!=0.){
      int id=0;
      for(int l=1;l<=pl;l++){
	for(int mm=0;mm<=l;mm++){ // snrz
	  leak[id]-=(2.*l+1.)*dtot
                   *sec->GetMesh(m).GetFlux(id+1).get_dat(g)
                           *mesh[m].GetFlux(id+1).get_dat(g)*vol;
	  id++;
	};
      };
      // (Higher leakage term)
      if(GetWrtflx(m)&&sec->GetWrtflx(m)){
        real tmp2=0.;
        int sn=GetQuadrature(g)->GetSnnum(); // snrz
        for(int i=0;i<sn;i++){
	  real om =GetQuadrature(g)->GetOmega(i);
          real flp=sec->GetAFlux(m,g).get_dat(i);
          int id=0;
	  for(int l=0;l<=pl;l++){
            real factor=(2.*l+1.)*invpi4;
	    for(int mm=0;mm<=l;mm++){ // snrz
      	      real mom=GetQuadrature(g)->GetMoment(id,i);
              flp-=factor*mom*sec->GetMesh(m).GetFlux(id).get_dat(g);
              id++;
            };
          };
	  tmp2+=om*flp*aflux[m][g].get_dat(i);
        };
	leak[plnum-1]-=fact*tmp2;
      };
      };
      // scattering term
      int id=0;
      for(int l=0;l<=pl;l++){
	for(int mm=0;mm<=l;mm++){ // snrz
  	  real tmp=0.;
          int st=g;
	  if(IsUpScatterSystem)st=0;
	  for(int gg=st;gg<grp;gg++){
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    if(dsigs!=0.)tmp+=dsigs*mesh[m].GetFlux(id).get_dat(gg);
	  };
  	  sct[id]+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	  id++;
	};
      };

    };
  };
};

real SNRZSystem::CalPerturbLeakScatSimplified(SNRZSystem *sec,bool *flag,int g)
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
        int sn=GetQuadrature(g)->GetSnnum(); // snrz
        for(int i=0;i<sn;i++){
	  real om =GetQuadrature(g)->GetOmega(i);
          real flp=sec->GetAFlux(m,g).get_dat(i);
	  tmp2+=om*flp*aflux[m][g].get_dat(i);
        };
	ret-=fact*tmp2;
      };
      // scattering term
      int id=0;
      for(int l=0;l<=pl;l++){
	for(int mm=0;mm<=l;mm++){ // snrz
  	  real tmp=0.;
          int st=g;
	  if(IsUpScatterSystem)st=0;
	  for(int gg=st;gg<grp;gg++){
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    if(dsigs!=0.)tmp+=dsigs*mesh[m].GetFlux(id).get_dat(gg);
	  };
  	  ret+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	  id++;
	};
      };

    };
  };

  return ret;
};

real SNRZSystem::CalPerturbLeakScatSimplifiedNew(SNRZSystem *sec,bool *flag,int g,int gg)
{
  real invpi4=1./PI4;

  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      real dtot=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt().get_dat(g)
                       -mesh[m].GetMed()->GetMacxs().GetSigt().get_dat(g);
      // total 
      if(dtot!=0.&&GetWrtflx(m)&&sec->GetWrtflx(m)){
        real fact=PI2*PI4*dtot*vol;
        real tmp2=0.;
        int sn=GetQuadrature(g)->GetSnnum(); // snrz
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
  	  for(int mm=0;mm<=l;mm++){ // snrz
  	    real tmp=0.;
	    real sigsp=sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real sigsu=mesh[m].GetMed()->GetMacxs().GetSigs(l).get_dat(g,gg);
  	    real dsigs=sigsp-sigsu;
	    if(dsigs!=0.){
              tmp=dsigs*mesh[m].GetFlux(id).get_dat(gg);
    	      ret+=tmp*sec->GetMesh(m).GetFlux(id).get_dat(g)*vol;
	    };
	    id++;
	  };
	};
      };

    };
  };

  return ret;
};

void SNRZSystem::CalSensitivity(SNRZSystem *sec,real k1,real k2,int nucnum,int *nucid)
{
  CheckAdjointForward(sec);

  bool fission_spectrum_matrix=false; // default:false
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

  cout.setf(ios::scientific);
  cout.precision(7);

  CheckSameMesh(sec); // o

  real ip=CalPerturbDenominator(sec);
  if(fission_spectrum_matrix){
    ip=CalPerturbDenominatorWithFissionSpectrumMatrix(sec);
  };

  if(ip==0.){
    cout<<"Error in SNRZsystem::CalSensitivity.\n";
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
    int nidendf=nid;
    //int nidendf=TranslateNuclideIDFromJFS(nid);

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
      cout<<nidendf<<"\n";
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
  	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micnu*micsigf);
            sec->GetMed(j).GetMacxs().GetSiga().add_data(i,den*micsigf);
  	    for(int l=0;l<=pl;l++){
	      sec->GetMed(j).GetMacxs().GetSigt(l).add_data(i,den*micsigf);
	    };
          };
        };
        real re=0.;
        CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
        for(int k=0;k<plnum;k++){
          re+=scttmp[k]+leaktmp[k];
        };
        if(!fission_spectrum_matrix){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
	}else{
          re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,k2);
	};
        re+=CalPerturbAbsorptionTerm(sec,flag,i);
        re/=ip;
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
            real den=med[j].GetNuclide(nid).GetDensity();
            nsforg[j]=sec->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
	    real micsigf=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
  	    real micnu=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
            sec->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*micsigf*micnu);
	  };
        };
        real re=0.;
	if(!fission_spectrum_matrix){
          re+=CalPerturbYieldTerm(sec,flag,i,k2);
	}else{
          re+=CalPerturbYieldTermWithFissionSpectrumMatrix(sec,flag,i,k2);
	};
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
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
          tot1org[j]=sec->GetMed(j).GetMacxs().GetSigt(1).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
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
      re/=ip;
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
      //for(int ii=0;ii<2;ii++){
      cout<<2<<"\n";
      cout<<nidendf<<"\n";
      //
      if(ii==0){cout<<2<<"\n";}
      else if(ii==1){cout<<4<<"\n";}
      else {cout<<16<<"\n";};
      //
      for(int i=0;i<grp;i++){
        for(int k=i;k<grp;k++){
	  if(ii!=0||nid<1000||(ii==0&&k<=i+2)){
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
                if(ii!=0.)micsigs=delta;
	        real dtot;
	        if(ii==2){
  	  	  dtot=micsigs*0.5;
	        }else{
   	          dtot=micsigs;
	        };
	        sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*micsigs);
	        sec->GetMed(j).GetMacxs().GetSigs(1).add_data(i,k,den*a1*micsigs);
	        for(int l=0;l<=pl;l++){
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

  delete [] fiss_frac;

  delete [] flag;

  delete [] scttmp;
  delete [] leaktmp;

  delete [] sigsorg;
  delete [] sigs1org;
  delete [] absorg;
  delete [] nsforg;
  delete [] totorg;
  delete [] tot1org;
};

SensitivityData SNRZSystem::CalSensitivityNew(SNRZSystem *sec,real keff,int nucnum,int *nucid,bool fiss_matrix)
{
  //cout<<"# Please modify (n,2n) sensitivity.\n"; exit(0);
  // Modification has been done in 2020/2/26
  // Same as SNRSystem (20131218)

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
  
  real ip=CalPerturbDenominator(sec);
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
          //re+=CalPerturbYieldTerm(sec,flag,i,keff);
  	  // +++++++++++++++++++++++++++++++++
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

real SNRZSystem::CalFluxMomentFromAngularFlux(int g,int mom,int meshid)
{
  if(!wrtflx[meshid]){
    cout<<"Error in CalFluxMomentFromAngularFlux of SNRZSystem\n";
    cout<<"You have to use PutWriteFlux method.\n";
    exit(0);
  };

  int qid=quadid[g];
  int snnum=quad[qid]->GetSnnum();
  real ret=0.;
  for(int j=0;j<snnum;j++){
    real tmp=aflux[meshid][g].get_dat(j)*quad[qid]->GetOmega(j)*PI2;
    ret+=tmp*quad[qid]->GetMoment(mom,j);
  };
  return ret;
};

GroupData1D SNRZSystem::GetIntegratedFlux(int medid,int mom)
{
  for(int i=0;i<nquad;i++){
    if(mom>quad[i]->GetPlnum()-1){
      cout<<"Pl number in quadrature should be increased\n";
      exit(0);
    };
  };

  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  for(int i=0;i<TotM;i++){
    if(mi.GetFMatParMesh(i)==medid){
      real vol=mesh[i].GetVolume();
      for(int g=0;g<grp;g++){
	if(mom<=plnum-1){
  	  sum[g]+=fabs(mesh[i].GetFluxData(mom,g))*vol;
	}else{
	  sum[g]+=fabs(CalFluxMomentFromAngularFlux(g,mom,i))*vol;
	};
      };
    };
  };
  GroupData1D ret(grp);
  ret.put_data(sum);
  delete []sum;
  return ret;
};


// For CMFD acceleration

void SNRZSystem::DoAcceleration(int iter,real errs,real fiss)
{
  if(iter%opt.GetItcmfd()==0&&iter>=opt.GetItmin_cmfd())DoCMFDAcceleration(fiss+0.5);
};


void SNRZSystem::DoCMFDAcceleration(real delk)
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
  };

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

  //string ss="Cartesian";
  //if(Cylinder)ss="Cylinder";
  string ss="Cylinder";
  cm.PutCartMeshInfo(cmi,ss);

  // Calculation of initial flux for CMFD 
  vector< vector<real> > FlxF(ngrp);  //[g][m]
  for(int i=0;i<ngrp;i++){
    FlxF[i].resize(CTotM,0.);
  };
  for(int i=0;i<CTotM;i++){
    real vol_inv=1./cm.GetMesh(i).GetVolume();
    for(int g=0;g<ngrp;g++){
      FlxF[g][i]=cc_flx[i]*vol_inv;
    };
  };

  GeneralOption option;
  real cm_epsf=opt.GetEpsf()*0.1;
  real cm_epsk=opt.GetEpsk()*0.1;
  real cm_epss=opt.GetEpss()*0.1;
  if(cm_epsf>1e-7)cm_epsf=1e-7;
  if(cm_epsk>1e-7)cm_epsk=1e-7;
  if(cm_epss>1e-7)cm_epss=1e-7;
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
    for(int g=0;g<ngrp;g++){
      cc_flx[i]=cm.GetMesh(i).GetFlux().get_dat(g)*tmp/cc_flx[i];
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
		int g2=0;
		real tmp=cc_flx[index];
		for(int g=0;g<grp;g++){
		  if(g>bgrp[g2]){
                    g2++;
		    tmp=cc_flx[index];
		  };
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

void SNRZSystem::SetZeroCurFF()
{
  int xr=mi.GetXC();
  int yr=mi.GetYC();
  //int zr=mi.GetZC(); // o
  int zr=0; // o

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

void SNRZSystem::CalCoarseCur()
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
	        //real sz=mesh[index].GetSurL(2);
  	        real sx2=mesh[index].GetSurR(0);
	        real sy2=mesh[index].GetSurR(1);
	        //real sz2=mesh[index].GetSurR(2);
  	        // X-direction
		if(x2==0){
                  //if(ix!=xedgel[iy]){
                  if(ix!=0){
		    CurFF[0][z1][y1][x1][0]+=(curxn[index]+curxp[index-1])*sx;
		  }else if(BC[0]!=1){
		    CurFF[0][z1][y1][x1][0]+=curxn[index]*sx;
		  };
		};
                //if(ix==xedger[iy]&&BC[1]!=1){
                if(ix==mi.GetXF()-1&&BC[1]!=1){
		  CurFF[0][z1][y1][xr][0]+=curxp[index]*sx2;
		};
	        // Y-direction
	        if(Ndim>1){
	  	  if(y2==0){
		    //if(iy!=yedgel[ix]){
		    if(iy!=0){
		      CurFF[0][z1][y1][x1][1]+=(curyn[index]+curyp[meshid[iz][iy-1][ix]])*sy;
		    }else if(BC[2]!=1){
		      CurFF[0][z1][y1][x1][1]+=curyn[index]*sy;
		    };
		  };
                  //if(iy==yedger[ix]&&BC[3]!=1){
                  if(iy==mi.GetYF()-1&&BC[3]!=1){
		    CurFF[0][z1][yr][x1][1]+=curyp[index]*sy2;
		  };
		};
	        // Z-direction
		/*
	        if(Ndim>2){
	  	  if(z2==0){
		    if(iz!=0){
		      CurFF[0][z1][y1][x1][2]+=(curzn[index]+curzp[meshid[iz-1][iy][ix]])*sz;
		    }else if(BC[4]!=1){
		      CurFF[0][z1][y1][x1][2]+=curzn[index]*sz;
		    };
		  };
                  if(iz==mi.GetZF()-1&&BC[5]!=1){
		    CurFF[0][zr][y1][x1][2]+=curzp[index]*sz;
		  };
		};
		*/
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

void SNRZSystem::ShowKeffContribution()
{
  GroupData1D abs;
  GroupData1D yld;
  GroupData1D n2n;
  GroupData1D leak;
  abs=GetIntegratedReactionRate(siga);
  yld=GetIntegratedReactionRate(nusigf);
  n2n=GetIntegratedReactionRate(sign2n);
  leak=GetLeakage();
  real abstot=abs.get_sum();
  real yldtot=yld.get_sum();
  real leaktot=leak.get_sum();
  real n2ntot=n2n.get_sum();
  real fact=1./yldtot;
  abs=abs*fact;
  yld=yld*fact;
  leak=leak*fact;
  n2n=n2n*fact;
  abstot*=fact;
  yldtot*=fact;
  leaktot*=fact;
  n2ntot*=fact;
  cout.setf(ios::showpoint);
  cout.precision(6);
  cout<<"# Total yield      : "<<yldtot<<"\n";
  for(int i=0;i<nmed;i++){
    real tmp=GetIntegratedReactionRate(nusigf,i).get_sum()*fact;
    cout<<"#  (Medium "<<i<<" : "<<tmp<<"\n";
  }; 
  cout<<"# Total n2n        : "<<n2ntot<<"\n";
  cout<<"# Total absorption : "<<abstot<<"\n";
  for(int i=0;i<nmed;i++){
    real tmp=GetIntegratedReactionRate(siga,i).get_sum()*fact;
    cout<<"#  (Medium "<<i<<" : "<<tmp<<"\n";
  }; 
  cout<<"# Total leakage    : "<<leaktot<<"\n";
  cout<<"# keff             : "<<(yldtot+n2ntot)/(abstot+leaktot)<<"\n";

  cout<<"\n";
  cout<<"#  yield, absorption, leakage, n2n\n";
  for(int i=0;i<grp;i++){
    real en=med[0].GetEnband().get_dat(i);
    cout<<en<<" "<<yld.get_dat(i)<<" "<<abs.get_dat(i)<<" "<<leak.get_dat(i)<<" "<<n2n.get_dat(i)<<"\n";
  };
};

real SNRZSystem::CalZCurrent(int r,int z,int g,bool pos)
{
  int m=meshid[0][z][r];

  int qid=quadid[g];
  int sn=quad[qid]->GetSnnum();

  real ret=0.;
  for(int i=0;i<sn;i++){
    real aflx=aflux[m][g].get_dat(i);
    real xi=quad[qid]->GetXi(i);
    if(pos&&xi>0.)ret+=aflx*xi;
    if(!pos&&xi<0.)ret+=-aflx*xi;
  };

  return ret;
};

SensitivityData SNRZSystem::CalSensitivityNewFixedSource(SNRZSystem *sec,real val,int nucnum,int *nucid)
{
  CheckAdjointForward(sec);
  CheckSameMesh(sec); 

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
      //cout<<"# capture : "<<i<<"\n";
      for(int j=0;j<nmed;j++){
	if(med[j].ExistNuclide(nid)){
          real den=med[j].GetNuclide(nid).GetDensity();
          absorg[j]=sec->GetMed(j).GetMacxs().GetSiga().get_dat(i);
          totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
	  real micsigc=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigc().get_dat(i);
          real tmp=den*micsigc;
          sec->GetMed(j).GetMacxs().GetSiga().add_data(i,tmp);
          sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,tmp);
        };
      };
      real re=CalPerturbLeakScatSimplifiedNew(sec,flag,i);
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

	  //cout<<"# elastic-p0 : "<<i<<" "<<k<<"\n";
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_dat(i,k);
	    if(micsigs!=0.){
              real den=med[j].GetNuclide(nid).GetDensity();
              totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
  	      real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigel(0).get_sumx().get_dat(i);
	      //real st0=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetData1d(sigel).get_dat(i);
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
	};
        
        real re=CalPerturbLeakScatSimplifiedNew(sec,flag,i,k);
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
    //egrp+=5; // 
    //if(egrp>grp)egrp=grp;

    sns2d.set_zero();
    for(int i=0;i<egrp;i++){
      //for(int i=0;i<grp;i++){
      for(int k=i;k<grp;k++){

        //cout<<"inelastic : "<<i<<" "<<k<<"\n";

        bool do_calc=false;
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSiginel(0).get_dat(i,k);
            if(micsigs!=0.)do_calc=true;
	  };
	};

	if(do_calc){
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta);
          };
	};
        real re=CalPerturbLeakScatSimplifiedNew(sec,flag,i,k);
        sns2d.put_data(i,k,re*PI4_INV/val);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
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
    //egrp+=5; // 
    //if(egrp>grp)egrp=grp;

    sns2d.set_zero();
    for(int i=0;i<egrp;i++){
      //for(int i=0;i<grp;i++){
      for(int k=i;k<grp;k++){

	//cout<<"n2n : "<<i<<" "<<k<<"\n";

        bool do_calc=false;
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real micsigs=sec->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSign2n(0).get_dat(i,k);
            if(micsigs!=0.)do_calc=true;
	  };
	};
	if(do_calc){

        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            real den=med[j].GetNuclide(nid).GetDensity();
            totorg[j]=sec->GetMed(j).GetMacxs().GetSigt(0).get_dat(i);
            sigsorg[0][j]=sec->GetMed(j).GetMacxs().GetSigs(0).get_dat(i,k);
            sec->GetMed(j).GetMacxs().GetSigs(0).add_data(i,k,den*delta);
            sec->GetMed(j).GetMacxs().GetSigt(0).add_data(i,den*delta*0.5);
          };
	};
        real re=CalPerturbLeakScatSimplifiedNew(sec,flag,i,k);
        sns2d.put_data(i,k,re*PI4_INV/val);
        for(int j=0;j<nmed;j++){
          if(med[j].ExistNuclide(nid)){
            sec->GetMed(j).GetMacxs().GetSigt(0).put_data(i,totorg[j]);
	    sec->GetMed(j).GetMacxs().GetSigs(0).put_data(i,k,sigsorg[0][j]);
	  };
	};

	};
      };
    };
    sens.PutSensitivity2D(nid,16,sns2d);

    // Elastic P1
    for(int i=0;i<grp;i++){

      //cout<<"# elastic-p1 : "<<i<<"\n";
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
      real re=CalPerturbLeakScatSimplified(sec,flag,i);
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

  delete [] flag;
  delete [] absorg;
  delete [] totorg;

  return sens;
};



void SNRZSystem::CalFixedSourceUpScat(real epsf, int oiter, bool print)
{
  /*
  cout<<edge.size()<<"\n";
  cout<<edge[0].size()<<"\n";
  for(int y=0;y<mi.GetYF();y++){
    cout<<y<<" "<<edge[0][y].size()<<"\n";
    //cout<<y<<" "<<edge[0][y][0][0]<<" "<<edge[0][y][0][1]<<"\n";
  };
  exit(0);
  */
 
  // - External source is stored in the scattering source array.
  // - Fission source is set to zero.

  CheckUpScattering();
  SetZeroFissionSrc();
  
  if(!IsUpScatterSystem){
    cout<<"#\n";
    cout<<"# Warning in GeneralSystem::CalFixedSourceUpScat.\n";
    cout<<"# Since this system is NOT up-scattering-included system,\n";
    cout<<"# [CalFixedSource] is runned instead.\n";
    cout<<"#\n";
    CalFixedSource(epsf,oiter,print);
    return;
  };


  // ... Calculation of error mode for two-grid acceleration
  for(int i=0;i<nmed;i++){
    med[i].FourierAnalysisForUpScattering();
  };

  //PLOSSystem init(2,grp,nmed);

  //PLOSSystem cm_pos(2,1,nmed);
  PLOSESystem cme_pos(2,1,nmed);    
  //PLOSSystem cm_neg(2,1,nmed);    
  //if(!print)cm_pos.NoPrint();
  if(!print)cme_pos.NoPrint();  
  //if(!print)cm_neg.NoPrint();

  int max_eigen=-1;
  real eigen_max=0.;
  for(int i=0;i<nmed;i++){
    real tmp=med[i].GetMacxs().GetData1d(dz).get_dat(0);
    if(tmp>eigen_max){      
      eigen_max=tmp;
      max_eigen=i;
    };
  };
  cout<<"# Medium having the maximum eigenvalue : "<<max_eigen<<"\n";
  // Medium-including check
  /*
  vector<bool> med_include(nmed,false);
  for(int i=0;i<TotM;i++){
    med_include[mi.GetFMatParMesh(i)]=true;
  };
  */

  cout<<"# 1-group cross section for two-grid acceleration.\n";
  for(int i=0;i<nmed;i++){
    Medium minp(1);
    minp.PutPL(0);
    real sigma_a=0.;
    real dc=0.;
    for(int g=0;g<grp;g++){
      real errormode=med[i].GetMacxs().GetData1d(dr).get_dat(g);
      sigma_a+=med[i].GetMacxs().GetData1d(siga).get_dat(g)*errormode;
      dc+=med[i].GetMacxs().GetData1d(d).get_dat(g)*errormode;
    };
    //if(sigma_a==0.||i!=max_eigen){
    //if(sigma_a==0.||(i!=6&&i!=8)){
    if(sigma_a==0.){            
      for(int g=0;g<grp;g++){
        sigma_a+=med[i].GetMacxs().GetData1d(siga).get_dat(g)/grp;
        dc+=med[i].GetMacxs().GetData1d(d).get_dat(g)/grp;
      };
    };
    minp.GetMacxs().GetData1d(d).put_data(0,dc);
    minp.GetMacxs().GetData1d(siga).put_data(0,sigma_a);
    minp.GetMacxs().GetData1d(sigt).put_data(0,sigma_a);
    minp.GetMacxs().GetData1d(nusigf).put_data(0,0.);
    minp.GetMacxs().GetData2d(sigs).put_data(0,0,0.);
    //cm_pos.AddMedium(minp);
    cme_pos.AddMedium(minp);    
    //cm_neg.AddMedium(minp);
    //init.AddMedium(med[i]);
    real eigen=med[i].GetMacxs().GetData1d(dz).get_dat(0);
    cout<<"# "<<i<<" "<<eigen<<" "<<dc<<" "<<sigma_a<<"\n";
  };
  //exit(0);


  //cm_pos.PutCartMeshInfo(mi,"Cylinder");
  cme_pos.PutCartMeshInfo(mi,"Cylinder");  
  //cm_neg.PutCartMeshInfo(mi,"Cylinder");  
  
  GeneralOption option;
  //cm_pos.PutGeneralOption(option);
  cme_pos.PutGeneralOption(option);  
  //cm_neg.PutGeneralOption(option);  

  //cm_pos.CalCoef();
  cme_pos.CalCoef();  
  //cm_neg.CalCoef();  

  //int itermax=10;
  int itermax=5000;  
  real err;

  vector<real> flxold(TotM);
  vector<real> R(TotM);

    // (initial guess of up-scattering source)
    if(oiter!=0){
      for(int g=0;g<grp;g++){
        if(g>=UpScatSrcGrp)AddUpScatSrc(g,pl);
        SetZeroUpScatSrc(g);
      };
    };

    // +++ Fixed source term calculation for thermal iteration
    for(int g=0;g<UpScatSinkGrp;g++){
      if(print)cout<<"# Fixed source calculation in group "<<g<<"\n";      
      CalSrcMultiplySystem(g,1.,pl);
      err=CalFluxGeneral(g,epsf,0);
      //cout<<"            "<<g<<" "<<err<<"\n";
      AddDownScatSrc(g,pl);
    };

    vector< vector< vector<real> > > src_fixed(TotM);
    for(int i=0;i<TotM;i++){
      src_fixed[i].resize(plnum);
      for(int l=0;l<plnum;l++){
        src_fixed[i][l].resize(grp);
        for(int g=UpScatSinkGrp;g<grp;g++){
          src_fixed[i][l][g]=mesh[i].GetScatSrc(g,l);
	};
      };
    };

    // +++ Thermal iteration
    vector<bool> convergence(grp,false);

    for(int i=0;i<itermax;i++){

      for(int j=0;j<TotM;j++){
	R[j]=0.;
      };
      
      bool conv=true;
      if(i==0)conv=false;
      real errmax=0.;      
      for(int g=UpScatSinkGrp;g<grp;g++){
	if(i!=0){
          for(int j=0;j<TotM;j++){
	    for(int l=0;l<plnum;l++){
      	      mesh[j].AddScatSrc(g,src_fixed[j][l][g],l);
	    };
	  };
	};	    
	if(!convergence[g]){

          //if(print)cout<<"# Fixed source calculation in group "<<g<<"\n";      	  

	  for(int j=0;j<TotM;j++){
	    flxold[j]=GetMesh(j).GetFlux().get_dat(g);
	  };
	  
  	  CalSrcMultiplySystem(g,1.,pl);
	  err=CalFluxGeneral(g,epsf,i);
	  err=fabs(err);
	  //if(g==100)cout<<i<<"         "<<g<<" : "<<err<<"\n";
	  cout<<"# "<<i<<"         "<<g<<" : "<<err<<"\n";
	  //if(err<epsf*0.1){
	  if(err<epsf){	    
	    convergence[g]=true;
	  };
	  if(err>errmax)errmax=err;

	  for(int j=0;j<TotM;j++){
	    real diff=GetMesh(j).GetFlux().get_dat(g)-flxold[j];
	    real tmp=0.;	    
	    for(int g2=UpScatSinkGrp;g2<g;g2++){
	      tmp+=GetMesh(j).GetMed()->GetMacxs().GetData2d(sigs).get_dat(g,g2)*diff;
	    };
  	    R[j]+=tmp;
	  };

	  
	};
        AddDownScatSrc(g,pl);
        AddUpScatSrc(g,pl);
        SetZeroScatSrc(g);
        SetZeroUpScatSrc(g);
        if(err>epsf)conv=false;
      };
      if(conv){
	cout<<"#    Number of iterations : "<<i<<"\n";
        break;
      };
      //cout<<i<<"    "<<errmax<<"\n";

      // Two-Grid acceleration
#if 1      
      if(i!=0){

      cout<<"#   ... Two-Grid acceleration ...\n";

      //cm_pos.SetZeroScatSrc();
      cme_pos.SetZeroScatSrc();      
      //cm_neg.SetZeroScatSrc();
      GroupData1D inp(1);
      //inp.set_zero();

      for(int j=0;j<TotM;j++){
        inp.put_data(0,R[j]); 	  
	//cm_pos.PutIsotropicSourcePerUnitVolume(j,inp);
	cme_pos.PutIsotropicSourcePerUnitVolume(j,inp);	
	
	/*
	if(R[j]>0.){
          inp.put_data(0,R[j]); 	  
  	  cm_pos.PutIsotropicSourcePerUnitVolume(j,inp);
        }else{
          inp.put_data(0,-R[j]); 	  	  
  	  cm_neg.PutIsotropicSourcePerUnitVolume(j,inp);
	};
	*/
      };
      
      //cm_pos.CalFixedSource();
      //cm_pos.CalFixedSource(1e-4,0);
      cme_pos.CalFixedSource(1e-4,0); // The second argument is [0] since if it is set as [i-1], the convergence becomes slow.      
      //cm_pos.CalFixedSource(1e-8,0);            
      //cm_neg.CalFixedSource(1e-4,0);      

      for(int g=UpScatSinkGrp;g<grp;g++){
	for(int j=0;j<TotM;j++){
	  real errmode=mesh[j].GetMed()->GetMacxs().GetData1d(dr).get_dat(g);
	  if(errmode>0.){
	    //cout<<j<<" "<<g<<" "<<mesh[j].GetFlux().get_dat(g)<<" "<<errmode*cm_pos.GetMesh(j).GetFlux().get_dat(0)<<"\n";
	    //mesh[j].GetFlux().add_data(g,errmode*cm_pos.GetMesh(j).GetFlux().get_dat(0));
	    mesh[j].GetFlux().add_data(g,errmode*cme_pos.GetMesh(j).GetFlux().get_dat(0));	    
	    //mesh[j].GetFlux().add_data(g,-errmode*cm_neg.GetMesh(j).GetFlux().get_dat(0));	    
	  };
	};
      };

      for(int g=UpScatSinkGrp;g<grp;g++){
	SetZeroUpScatSrc(g);
	AddUpScatSrc(g,pl);
      };

      

      };

#endif      
      // ...
      
      if(i==itermax-1){
        cout<<"# !! Caution !!\n";
        cout<<"# Fixed source calculation is not converged.\n";
        cout<<"# Maximum residual is "<<err<<"\n";
        cout<<"# Convergence criteria "<<epsf<<"\n";
      };
    };


};
