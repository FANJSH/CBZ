#include<cstdlib>
#include "MEC_system.h"

using namespace std;

MECSystem::MECSystem(int i,int j):GeneralSystem(2,i,j)
{
  pl=-1;
  name="MEC";
  cmfdimp=true; // CMFD IS implemented
  Vacuum=false;
  cmesh=-1;
  cmr=true;
  transport=true;
  degree360=false;
  wrtflx=false;
  set_array=false;
  already_run=false;
};

void MECSystem::SetInitialFlux()
{
  SetArray();
  if(!already_run)SetInitialFlatFlux();
  already_run=true;
};

void MECSystem::SetArray()
{
  if(set_array){
    //cout<<"# Error in MECSytem::SetArray.\n";
    //cout<<"# Array has been already set.\n";
    //exit(0);
    cout<<"# Warning in MECSytem::SetArray.\n";
    cout<<"# Array has been already set.\n";
    //return;
  };

  set_array=true;

  if(pl==-1){
    cout<<"# Do 'PutPL' for system class.\n";
    exit(0);
  };

  if(pl>0&&transport){
    for(int i=0;i<nmed;i++){
      for(int j=1;j<pl+1;j++){
        for(int k=0;k<grp;k++){
	  real tmp=med[i].GetDataSigs(j,k,k);
	  real tmp2=med[i].GetDataSigt(0,k)-med[i].GetDataSigt(j,k);
	  real tmp3=tmp+(j*2+1.)*tmp2;
	  real tmp4=med[i].GetSigs(0).get_dat(k,k);
	  if(j==1&&fabs(tmp3*0.33333333)>tmp4){
	    cout<<"   ... warning...sigs_p1 > sigs_p0 \n";
	    cout<<"     sigs_p1/sigs_p0 : "<<tmp3*0.33333333/tmp4<<"\n";
	    cout<<"   medium "<<i<<"  group "<<k<<"\n";
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

  if(wrtflx){
    int totsn=quad.GetSN();
    aflux.resize(TotM);
    for(int i=0;i<TotM;i++){
      aflux[i].resize(grp);
      for(int j=0;j<grp;j++){
	aflux[i][j].put_imax(totsn);
      };
    };
  };
};

void MECSystem::PutTrajectorySet(TrajectorySet *cinp,bool deg720)
{
  if(cinp->GetBoundaryCondition()==Black){
    PutVacuum();
  }else if(cinp->GetBoundaryCondition()!=Periodic&&cinp->GetBoundaryCondition()!=Reflective){
    cout<<"# Error in MECSystem::PutTrajectorySet.\n";
    cout<<"# Black or Periodic B.C. can be chosen.\n";
    exit(0);
  };

  tset=cinp;
  TotM=tset->GetRegnum();
  if(print)cout<<"*** Total Mesh : "<<TotM<<"\n";
  mesh.resize(TotM);
  medid_fmesh.resize(TotM);
  for(int i=0;i<TotM;i++){
    mesh[i].PutVolume(tset->GetVol(i));
    mesh[i].PutDim(2);
    mesh[i].PutPL(1);
    mesh[i].PutGrp(grp);
  };
  degree360=tset->IsAngle360();

  max_segments=cinp->GetMaximumSegments();

  // Leonard
  /*
  PutPolarAngleDivision(2);
  polar_angle_sin_inv[0]=1./0.273658;
  polar_angle_sin_inv[1]=1./0.865714;
  polar_angle_weight[0]=0.139473;
  polar_angle_weight[1]=0.860527;
  */
  // TY 3-point
  /*
  PutPolarAngleDivision(3);
  polar_angle_sin_inv[0]=1./0.166648;
  polar_angle_sin_inv[1]=1./0.537707;
  polar_angle_sin_inv[2]=1./0.932954;
  polar_angle_weight[0]=0.046233;
  polar_angle_weight[1]=0.283619;
  polar_angle_weight[2]=0.670148;
  */
  // TY 2-point

  PutPolarAngleDivision(2);
  polar_angle_sin_inv[0]=1./0.3639;
  polar_angle_sin_inv[1]=1./0.8999;
  polar_angle_weight[0]=0.212854;
  polar_angle_weight[1]=0.787146;

  // Gaussian
  /*
  int pnt=3;
  real cosin[]={0.238619, 0.661209, 0.932469};
  real wetin[]={0.467909, 0.36076, 0.171331};
  */
  /*
  int pnt=4;
  real cosin[]={0.183435, 0.525532, 0.796666, 0.96029};
  real wetin[]={0.362684, 0.313707, 0.222381, 0.101228};
  */
  /*
  int pnt=6;
  real cosin[]={0.125233, 0.367831, 0.587318, 0.769903, 0.904117, 0.981561};
  real wetin[]={0.249147, 0.233493, 0.203167, 0.160078, 0.106939, 0.0471752};
  */
  /*
  int pnt=4; // DG
  real cosin[]={0.0694318, 0.330009, 0.669991, 0.930568};
  real wetin[]={0.173927, 0.326073, 0.326073, 0.173927};
  */

  // (for Gaussian)
  /*
  PutPolarAngleDivision(pnt);
  for(int i=0;i<pnt;i++){
    real sinin=sqrt(1.-cosin[i]*cosin[i]);
    polar_angle_sin_inv[i]=1./sinin;
    polar_angle_weight[i]=wetin[i];
  };
  */

  PutAzimuthalAngleDivision(tset->GetAngularDivision());

  int nt=tset->GetNumTrajectory();
  if(degree360)nt*=2;
  int pad=polar_angle_division;
  if(deg720)pad*=2;
  if(deg720)influx.resize(pad);
  for(int i=0;i<pad;i++){
    influx[i].resize(nt);
    for(int j=0;j<nt;j++){
      influx[i][j].resize(grp,0.);
    };
  };

  quad.Initialize(2,1); // 2D, P1

  // +++ The following is under construction
  // Only for 90 degree
  int totsn=polar_angle_division*azimuthal_angle_division;
  if(degree360)totsn*=2;
  quad.PutSN(totsn);
  real *mui=new real[totsn];
  real *eai=new real[totsn];
  real *xii=new real[totsn];
  real *w=new real[totsn];

  int sn=0;
  for(int i=0;i<polar_angle_division;i++){
    real s=1./polar_angle_sin_inv[i];
    real z=sqrt(1-s*s);
    real xy=sqrt(1-z*z);
    int aza=azimuthal_angle_division;
    for(int j=0;j<aza;j++){
      real theta=tset->GetDiscretedAngle(j)/180.*PI;
      mui[sn]=cos(theta)*xy;
      eai[sn]=sin(theta)*xy;
      xii[sn]=z;
      w[sn]=polar_angle_weight[i]/azimuthal_angle_division;
      sn++;
    };
    if(degree360){
      for(int j=0;j<aza;j++){
        real theta=tset->GetDiscretedAngle(j)/180.*PI;
        mui[sn]=-cos(theta)*xy;
        eai[sn]=-sin(theta)*xy;
        xii[sn]=z;
        w[sn]=polar_angle_weight[i]/azimuthal_angle_division;
        sn++;
      };
    };
  };
  quad.PutData(mui,eai,xii,w);
  real sum=PI4;
  quad.WeightNormalize(sum); // This is different from quadratures for SN codes
  quad.CalValue();

  delete [] mui;
  delete [] eai;
  delete [] xii;
  delete [] w;
};

void MECSystem::PutPolarAngleDivision(int i)
{
  polar_angle_division=i;
  polar_angle_sin_inv.resize(i);
  polar_angle_weight.resize(i);
  influx.resize(i);
};

void MECSystem::PutAzimuthalAngleDivision(int i)
{
  azimuthal_angle_division=i;
  azimuthal_angle_weight.resize(i);
};

void MECSystem::PutRelationRegionMedium(vector<int> &inp)
{
  for(int i=0;i<TotM;i++){
    if(inp[i]<0||inp[i]>=nmed){
      cout<<"# Medium ID is not appropriated in MECSystem!\n";
      cout<<"# Medium ID : "<<inp[i]<<"\n";
      cout<<"# Medium ID should range from 0 to "<<nmed-1<<"\n";
      exit(0);
    };
    mesh[i].PutMedium(&med[inp[i]]);
    medid_fmesh[i]=inp[i];
  };
};

void MECSystem::PutRelationRegionMedium(int *inp)
{
  vector<int> inp2(TotM);
  for(int i=0;i<TotM;i++){
    inp2[i]=inp[i];
  };
  PutRelationRegionMedium(inp2);
};

real MECSystem::CalFluxGeneral(int ng,real cin,int iter)
{
  // Original version of this method is renamed as "CalFluxGeneralDetail" in 2012/1/24.
  //
  // CMR/CMFD/P1 options are discarded in the present version.

  if(degree360){
    if(Vacuum){
      cout<<"# Error in MECSystem::CalFluxGeneral.\n";
      cout<<"# Whole angle calculation (360-degree) cannot be done under vacuum BCs.\n";
      exit(0);
    };
    return CalFluxAngle360(ng,cin,iter);
    //return CalFluxAngle360Old(ng,cin,iter,false); // to calculate neutron current
  };

  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM*(pl+1)];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src_out=new real[TotM*plnum];  // Out-group source
  real *flux=new real[TotM*plnum]; // Flux moment

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source to `src_out' (angluar moment-wise)
  PutSourceInnerIteration(src_out);

  int aad=azimuthal_angle_division;
  int pad=polar_angle_division;
  int tot_ang=pad*aad;
  vector< vector<real> > aveflx(tot_ang,(vector<real>)TotM);
  // volume-integrated mesh-wise angular flux

  int innermax=1000;
  //real convergence_condition_influx=cin*0.05;

  // At the initial outer iteration, incoming angular flux in the g-th group
  // is assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&!Vacuum){
    int ng2=ng-1;
    if(!opt.Forward())ng2=ng+1;
    if(ng2>=0&&ng2<ng){
      for(int i=0;i<polar_angle_division;i++){
        for(int j=0;j<numt;j++){
	  influx[i][j][ng]=influx[i][j][ng-1];
        };
      };
    };
  };

  real *putsrc=new real[max_segments]; // Source term
  real *ave=new real[max_segments];    // Mesh-averaged flux
  vector<int> reg_id(max_segments);

  real consttmp=aad*INV_PI;
  for(int inner=0;inner<innermax;inner++){

    for(int j=0;j<tot_ang;j++){
      for(int i=0;i<TotM;i++){
        aveflx[j][i]=0.;
      };
    };

    int azm_angle_id=0;
    for(int tj=0;tj<numt;tj++){

      int r=tset->GetTrajectory(tj).GetNumreg();
      for(int i=0;i<r;i++){
	reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
      };
      real width=tset->GetTrajectory(tj).GetWeight()*consttmp;

      for(int pp=0;pp<pad;pp++){ // polar-angle division

        int sn=pp*aad+azm_angle_id;
        //real omega=quad.GetOmega(sn);

        real sin_inv=polar_angle_sin_inv[pp];

	// +++ source calculation
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  putsrc[i]=src_out[id]+selfscat[id]*flux[id];
	  /*
	  real sum=0.;
          int pos=id*plnum;
	  for(int l=0;l<plnum;l++){
  	    sum+=(src_out[pos]+selfscat[pos]*flux[pos])*quad.GetMoment(l,sn);
	    pos++;
	  };
	  putsrc[i]=sum;
	  */
        };

	// (Step characteristics)
	//real ooo;
	/*
        if(Vacuum){
	// +++ vacuum boundary
          //tset->GetTrajectory(tj).SolveCFormTransport(fmoc,0.,putsrc,ave,esc,polar_angle_sin_inv[pp]);
	}else{
	*/
     	  real outflx;
          tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,ave,outflx,sin_inv);
  	  int next_id=tset->GetNextTrajectory(tj);
          influx[pp][next_id][ng]=outflx;
	  //};

        for(int i=0;i<r;i++){
	  aveflx[sn][reg_id[i]]+=ave[i]*width;
        };

      };

      if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id))azm_angle_id++;

    };

    // +++ Renew of angular moment of mesh-wise flux
    real errmax=0.;
    {
      int index=0;
      for(int i=0;i<TotM;i++){
        real tmp=0.;
        for(int j=0;j<tot_ang;j++){
          tmp+=aveflx[j][i]*quad.GetOmega(j);
        };
        tmp/=mesh[i].GetVolume();
	real err=fabs(tmp/flux[index]-1.);
	if(err>errmax)errmax=err;
        flux[index]=tmp;
	index++;
      };
    };

    bool loop_fin=false;
    if(inner==innermax-1){
      cout<<"# Inner iteration does not converge in group "<<ng<<".\n";
      cout<<"# Iteration residual is "<<errmax<<"\n";
      loop_fin=true;
    };

    if(errmax<cin){
      loop_fin=true;
      //cout<<"   (grp:"<<ng<<") "<<inner<<" "<<errmax<<"\n";
      //  break;
    };

    // ---- system-rebalancing -----------------------------------
    //     (valid only for zero-leakage B.C.)
    //if(!cmr&&!Vacuum){
    if(!cmr&&!Vacuum&&!loop_fin){
      real srcsum=0.;
      real colsum=0.;
      for(int i=0;i<TotM;i++){
	real vol=mesh[i].GetVolume();
	srcsum+=src_out[i*plnum]*vol*PI4;
	colsum+=flux[i*plnum]*vol*(xs[i]-selfscat[i*(pl+1)]*PI4);
      };
      if(colsum>0.){
        real factor=srcsum/colsum;
        for(int i=0;i<TotM;i++){
          flux[i]*=factor; // rebalance factor
        };
        for(int i=0;i<polar_angle_division;i++){
	  for(int j=0;j<numt;j++){
            influx[i][j][ng]*=factor;
	  };
        };
      };
    };

    if(loop_fin)break;

  };// --- inner iteration loop-end

  // +++ flux-moment renew & error calculation of scalar flux
  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i*plnum]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    for(int l=0;l<plnum;l++){
      mesh[i].GetFlux(l).put_data(ng,flux[i*plnum+l]);
    };
    //if(!opt.Forward())cout<<ng<<" "<<i<<" "<<flux[i]<<"\n";
  };

  if(wrtflx){
    for(int i=0;i<TotM;i++){
      real vol_inv=1./mesh[i].GetVolume();
      for(int j=0;j<tot_ang;j++){
	aflux[i][ng].put_data(j,aveflx[j][i]*vol_inv);
      };
    };
  };

  delete [] xs;
  delete [] selfscat;
  delete [] src_out;
  delete [] flux;
  delete [] putsrc;
  delete [] ave;

  return maxerr;
};

real MECSystem::CalFluxGeneralDetail(int ng,real cin,int iter)
{
  // This method was seperated and stored in 2012/1/24.
  //
  // CMR/CMFD/P1 options were tried, (but failed).
  // Detail can be found in the notebook in 2008/10/29, 11/10, 12/16


  if(degree360)return CalFluxAngle360(ng,cin,iter);

  // +++ 180-degree symmetry +++

  bool cmfd_on=false;
  if(cmfd&&iter%opt.GetItcmfd()==0)cmfd_on=true;

  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM*(pl+1)];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src=new real[TotM*plnum];  // External source
  real *flux=new real[TotM*plnum]; // Flux moment
  real *ssc=new real[TotM*plnum];  // Self-scattering source
  real *flxd=new real[TotM]; // Old scalar-flux

  // ++ current data storing : crt[cmesh][0-4] (cmesh -> 0-4)
  //
  //           1
  //           ^
  //           |
  //   2 <- (Reg(4)) -> 0
  //           |
  //           Y
  //           3
  //
  //    (4 means ''a current to outer boundary'')
  vector< vector<real> > crt;
  if(cmr||cmfd_on){
    crt.resize(cmesh);
    for(int i=0;i<cmesh;i++){
      crt[i].resize(5);
    };
  };

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source to `src' (angluar moment-wise)
  PutSourceInnerIteration(src);

  int aad=azimuthal_angle_division;
  int pad=polar_angle_division;
  int tot_ang=pad*aad;
  vector< vector<real> > aveflx(tot_ang); // volume-integrated mesh-wise angular flux
  for(int i=0;i<tot_ang;i++){
    aveflx[i].resize(TotM,0.);
  };

  int innermax=1000;
  real convergence_condition_influx=cin*0.05;

  // At the initial outer iteration,
  // incoming angular flux in g-th group is 
  // assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&ng!=0&&!Vacuum){
    for(int i=0;i<polar_angle_division;i++){
      for(int j=0;j<numt;j++){
	influx[i][j][ng]=influx[i][j][ng-1];
      };
    };
  };

  real consttmp=aad*INV_PI;
  for(int inner=0;inner<innermax;inner++){

    for(int i=0;i<TotM;i++){flxd[i]=flux[i*plnum];};

    // +++ self-scattering source 
    int index=0;
    int index2=0;
    for(int i=0;i<TotM;i++){
      for(int l=0;l<=pl;l++){
	real tmp=selfscat[index2++];
	for(int m=0;m<=l;m++){
          ssc[index]=flux[index]*tmp;
	  index++;
	};
      };
    };

    // +++ initialize for current-matrix
    if(cmr||cmfd_on){
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<5;j++){
          crt[i][j]=0.;
        };
      };
    };

    for(int j=0;j<tot_ang;j++){
      for(int i=0;i<TotM;i++){
        aveflx[j][i]=0.;
      };
    };

    real maxerr_bflx=0.;
    int sn=0;
    real omega=quad.GetOmega(sn);
    real mu=quad.GetMu(sn);
    real et=quad.GetEata(sn);
    real sqq=sqrt(mu*mu+et*et);

    real sin_inv=1.;
    //real sin_inv_old=1.;
    for(int pp=0;pp<pad;pp++){

      //sin_inv_old=sin_inv;
      sin_inv=polar_angle_sin_inv[pp];

      int azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){

        int r=tset->GetTrajectory(tj).GetNumreg();
	vector<int> reg_id(r);
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real width=tset->GetTrajectory(tj).GetWeight()*consttmp;

	real *putsrc=new real[r]; // Source term
	real *ave=new real[r];    // Mesh-averaged flux
        real *esc=new real[r];    // Mesh-escaping flux

	// +++ source calculation
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=0.;
	  for(int l=0;l<plnum;l++){
	    real tmp=src[id*plnum+l]+ssc[id*plnum+l];
  	    sum+=tmp*quad.GetMoment(l,sn);
	  };
	  putsrc[i]=sum;
        };

	// (Step characteristics)
        if(Vacuum){
	// +++ vacuum boundary
          //tset->GetTrajectory(tj).SolveCFormTransport(etab,0.,putsrc,ave,esc,polar_angle_sin_inv[pp]);
          tset->GetTrajectory(tj).SolveCFormTransport(fmoc,0.,putsrc,ave,esc,polar_angle_sin_inv[pp]);
	}else{
          //tset->GetTrajectory(tj).SolveCFormTransport(etab,influx[pp][tj][ng],putsrc,ave,esc,polar_angle_sin_inv[pp]);
          tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,ave,esc,sin_inv);
	};
	// (Diamond differencing)
	/*
        if(!Vacuum){
          tset->GetTrajectory(tj).SolveCFormTransportDD(influx[pp][tj][ng],putsrc,ave,esc,polar_angle_sin_inv[pp]);
	}else{
	// +++ vacuum boundary
          tset->GetTrajectory(tj).SolveCFormTransportDD(0.,putsrc,ave,esc,polar_angle_sin_inv[pp]);
	};
	*/
	real outflx=esc[r-1];

	// -- For CMR --------------------------------------------------------------
	if(cmr||cmfd_on){
          //real factor=omega*width*sqq*0.125; // for 45degree symmetric
          real factor=omega*width*sqq*0.5; // for 180degree symmetric
	  for(int i=0;i<r-1;i++){
	    int id1=reg_id[i];
	    int id2=reg_id[i+1];
            int cd1=cmeshid[id1];
	    int cd2=cmeshid[id2];
	    if(cd1!=cd2){
	      if(cd1+1==cd2){
		crt[cd1][0]+=esc[i]*factor;
	      }else if(cd1-1==cd2){
		crt[cd1][2]+=esc[i]*factor;
	      }else if(cd1+cmeshx==cd2){
		crt[cd1][1]+=esc[i]*factor;
	      }else{
  	        //if(inner==0)cout<<pp<<" "<<tj<<" : "<<cd1<<" "<<cd2<<"\n";
		crt[cd1][0]+=esc[i]*factor*0.5;
		crt[cd1][1]+=esc[i]*factor*0.5;
	      };
	    };
	  };
          if(!Vacuum){
	    int next_id=tset->GetNextTrajectory(tj);
            int next_reg=tset->GetTrajectory(next_id).GetReg(0);
            int cd1=cmeshid[reg_id[r-1]];
            int cid=cmeshid[next_reg];
	    if(cd1-cmeshx+1==cid){
	      crt[cd1][0]+=esc[r-1]*factor;
	    }else if(cd1+cmeshx-1==cid){
	      crt[cd1][2]+=esc[r-1]*factor;
	    }else if(cid+cmeshx==cd1){
	      crt[cd1][1]+=esc[r-1]*factor;
	    };
	  }else{
  	    int id1=reg_id[r-1];
	    int cd1=cmeshid[id1];
	    crt[cd1][4]+=esc[r-1]*factor; // Outward current in outer boundary
	  };
	};
	// ------------------------------------------------------------------------
	// 
	if(!Vacuum){
  	  int next_id=tset->GetNextTrajectory(tj);
  	  real err=fabs(influx[pp][next_id][ng]/outflx-1.);
	  if(err>maxerr_bflx)maxerr_bflx=err;
          influx[pp][next_id][ng]=outflx;
        };
        for(int i=0;i<r;i++){
	  aveflx[sn][reg_id[i]]+=ave[i]*width;
        };
	delete [] putsrc;
	delete [] ave;
        delete [] esc;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Shift to the next azimuthal-direction trajectories
	  sn++;
          azm_angle_id++;
	  if(tj==numt-1)sn+=aad;
	  if(sn!=tot_ang){ // Set-up for next azimuthal angle calculations
            omega=quad.GetOmega(sn);
            real mu=quad.GetMu(sn);
            real et=quad.GetEata(sn);
            sqq=sqrt(mu*mu+et*et);
	  };
	};

      };
    };

    // +++ Renew of angular moment of mesh-wise flux
    {
      int index=0;
      for(int i=0;i<TotM;i++){
        real tmp=0.;
        real tmp_x=0.;
        real tmp_y=0.;
        for(int j=0;j<tot_ang;j++){
          real flom=aveflx[j][i]*quad.GetOmega(j);
          tmp+=flom;
          if(pl>0){
  	    tmp_x+=flom*quad.GetMu(j);
	    tmp_y+=flom*quad.GetEata(j);
	  };
        };
        real vinv=1./mesh[i].GetVolume();
	tmp*=vinv;
        flux[index++]=tmp;
        if(pl>0){
	  tmp_x*=vinv;
	  tmp_y*=vinv;
          flux[index++]=tmp_x;
          flux[index++]=tmp_y;
        };
      };
    };

    // +++ Symmetric correction for flux
    /*
    int ttt[]={4,5,6,7,0,1,2,3};
    regid_180.resize(8);
    for(int i=0;i<8;i++){
      regid_180[i]=ttt[i];
    };

    vector<real> flux_org(TotM*plnum);
    */
    /*
    // [45 degree]
    for(int i=0;i<TotM*plnum;i++){flux_org[i]=flux[i];};
    for(int i=0;i<TotM;i++){
      int tmp=regid_45[i];
      flux[tmp*plnum]+=flux_org[i*plnum];
      if(pl>0){
	flux[tmp*plnum+1]+=flux_org[i*plnum+2]; // cur_x
	flux[tmp*plnum+2]+=flux_org[i*plnum+1]; // cur_y
      };
    };
    // [90 degree]
    for(int i=0;i<TotM*plnum;i++){flux_org[i]=flux[i];};
    for(int i=0;i<TotM;i++){
      int tmp=regid_90[i];
      flux[tmp*plnum]+=flux_org[i*plnum];
      if(pl>0){
	flux[tmp*plnum+1]-=flux_org[i*plnum+1]; // cur_x
	flux[tmp*plnum+2]+=flux_org[i*plnum+2]; // cur_y
      };
    };
    */
    // [180 degree]
    /*
    for(int i=0;i<TotM*plnum;i++){flux_org[i]=flux[i];};
    for(int i=0;i<TotM;i++){
      int tmp=regid_180[i];
      flux[tmp*plnum]+=flux_org[i*plnum];
      if(pl>0){
	flux[tmp*plnum+1]+=flux_org[i*plnum+1]; // cur_x
	flux[tmp*plnum+2]-=flux_org[i*plnum+2]; // cur_y
      };
    };
    for(int i=0;i<TotM*plnum;i++){
      //flux[i]*=0.125; // 45 degree symmetry
      flux[i]*=0.5; // 180 degree symmetry
    };
    */
    // ---- system-rebalancing -----------------------------------
    //     (valid only for relfective B.C.)
    if(!cmr&&!Vacuum){
      real srcsum=0.;
      real sscsum=0.;
      real stsum=0.;
      for(int i=0;i<TotM;i++){
        real vol=mesh[i].GetVolume();
        srcsum+=src[i*plnum]*vol;
        sscsum+=flux[i*plnum]*selfscat[i*(pl+1)]*vol;
        stsum+=xs[i]*flux[i*plnum]*vol;
      };
      srcsum*=PI4;
      sscsum*=PI4;
      real factor=srcsum/(stsum-sscsum);
      for(int i=0;i<TotM;i++){
        flux[i*plnum]*=factor; // rebalance factor
      };
    };
    // ------------------------------------------------------------

    real errmax=0.;
    for(int i=0;i<TotM;i++){
      real err=fabs(flux[i*plnum]/flxd[i]-1.);
      if(err>errmax)errmax=err;
    };
    if(errmax<cin){
      if((!Vacuum&&maxerr_bflx<convergence_condition_influx)||Vacuum){
        //cout<<"   (grp:"<<ng<<") "<<inner<<"\n";
        break;
      };
    };

    // +++ Symmetric correction for current ---------------------------
    if(cmr||cmfd_on){
      // symmetric effect for current-matrix 
      vector< vector<real> > crt_org(cmesh);
      for(int i=0;i<cmesh;i++){
        crt_org[i].resize(5);
	for(int j=0;j<5;j++){
	  crt_org[i][j]=crt[i][j];
	};
      };
      /*
      //  - 45 degree symmetric
      for(int y=0;y<cmeshy;y++){
        for(int x=0;x<cmeshx;x++){
	  int id1=y*cmeshx+x;
          int nd1=x*cmeshx+y; 
	  crt[nd1][0]+=crt_org[id1][1];
	  crt[nd1][1]+=crt_org[id1][0];
	  crt[nd1][2]+=crt_org[id1][3];
	  crt[nd1][3]+=crt_org[id1][2];
	  crt[nd1][4]+=crt_org[id1][4];
        };
      };
      // - 90 degree symmtric
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<5;j++){crt_org[i][j]=crt[i][j];};
      };
      for(int y=0;y<cmeshy;y++){
        for(int x=0;x<cmeshx;x++){
  	  int id1=y*cmeshx+x;
          int nd1=y*cmeshx+(cmeshx-1-x);
	  crt[nd1][0]+=crt_org[id1][2];
  	  crt[nd1][1]+=crt_org[id1][1];
	  crt[nd1][2]+=crt_org[id1][0];
	  crt[nd1][3]+=crt_org[id1][3];
	  crt[nd1][4]+=crt_org[id1][4];
        };1
      };
      */
      // - 180 degree symmtric
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<5;j++){crt_org[i][j]=crt[i][j];};
      };
      for(int y=0;y<cmeshy;y++){
        for(int x=0;x<cmeshx;x++){
 	  int id1=y*cmeshx+x;
          int nd1=(cmeshy-1-y)*cmeshx+x;
	  crt[nd1][0]+=crt_org[id1][0];
	  crt[nd1][1]+=crt_org[id1][3];
	  crt[nd1][2]+=crt_org[id1][2];
	  crt[nd1][3]+=crt_org[id1][1];
	  crt[nd1][4]+=crt_org[id1][4];
        };
      };
    };
    // -------------------------------------------------

    // --- Coarse Mesh Rebalance Acceleration -------------------
    if(cmr){

      /*
      for(int i=0;i<4;i++){
	cout<<i<<" ";
	real sum=0.;
	for(int j=0;j<5;j++){
	  cout<<crt[i][j]<<" ";
	  sum+=crt[i][j];
	};
	cout<<" : "<<sum<<"\n";
      };
      */

      vector<real> lhs_diag(cmesh);
      vector<real> rhs(cmesh);
      for(int i=0;i<cmesh;i++){
	lhs_diag[i]=0.;
	rhs[i]=0.;
      };
      for(int i=0;i<cmesh;i++){
	int stt=0;
	if(i!=0)stt=cmesh_bound[i-1]+1;
        real srcsum=0.;
        real sscsum=0.;
        real stsum=0.;
        for(int j=stt;j<=cmesh_bound[i];j++){
  	  real vol=mesh[j].GetVolume();
          real tmp=flux[j*plnum]*vol;
	  srcsum+=src[j*plnum]*vol;
	  sscsum+=tmp*selfscat[j*(pl+1)];
	  stsum+=xs[j]*tmp;
        };
        srcsum*=PI4;
        sscsum*=PI4;
	real tmpdiag=stsum-sscsum;
	for(int j=0;j<5;j++){
  	  tmpdiag+=crt[i][j];
	};
	cout<<i<<" "<<srcsum<<" "<<sscsum<<" "<<stsum<<" "<<tmpdiag<<"\n";
	if(i==3)exit(0);
        lhs_diag[i]=1./tmpdiag;
        rhs[i]=srcsum;
      };

      vector<real> sol(cmesh);
      for(int i=0;i<cmesh;i++){sol[i]=rhs[i];};

      int itermax=1000000;
      real eps=1e-9;

      for(int iter=0;iter<itermax;iter++){
	int id=0;
	real errmax=0.;
	for(int y=0;y<cmeshy;y++){
	  for(int x=0;x<cmeshx;x++){
	    real tmp=rhs[id];
	    if(x!=0)        tmp+=crt[id-1][0]*sol[id-1];
	    if(x!=cmeshx-1) tmp+=crt[id+1][2]*sol[id+1];
	    if(y!=0)        tmp+=crt[id-cmeshx][1]*sol[id-cmeshx];
	    if(y!=cmeshy-1) tmp+=crt[id+cmeshx][3]*sol[id+cmeshx];
	    tmp*=lhs_diag[id];
	    real err=fabs(tmp/sol[id]-1.);
	    if(err>errmax)errmax=err;
	    sol[id]=tmp;
	    id++;
	  };
	};
	if(errmax<eps)break;
      };

      for(int i=0;i<cmesh;i++){
        int stt=0;
        if(i!=0)stt=cmesh_bound[i-1]+1;
        for(int j=stt;j<=cmesh_bound[i];j++){
	  for(int l=0;l<plnum;l++){
  	    flux[j*plnum+l]*=sol[i];
	  };
        };
      };
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<5;j++){
	  crt[i][j]*=sol[i];
        };
      };
    };
    // ---------------------------------------------
    // +++ end of CMR Acceleration

    if(inner==innermax-1){
      cout<<"# !! Caution !!";
      cout<<" Inner iteration is not converged. ("<<ng<<")\n";
    };
  };// --- inner iteration loop-end

  // +++ flux-moment renew & error calculation of scalar flux
  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i*plnum]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    for(int l=0;l<plnum;l++){
      mesh[i].GetFlux(l).put_data(ng,flux[i*plnum+l]);
      //cout<<i<<" "<<flux[i*plnum+l]<<"\n";
    };
  };

  // --- for CMFD acceleration ----------------------------
  if(cmfd_on){
    if(ng==0)SetZeroCurFF();
    // 
    for(int y=0;y<cmeshy;y++){
      for(int x=0;x<cmeshx;x++){
        int nd1=y*cmeshx+x;
        // X- plane
	if(x!=0){
          CurFF[0][0][y][x][0]-=crt[nd1][2];
	};
        // X+ plane
	if(x!=cmeshx-1){
          CurFF[0][0][y][x+1][0]+=crt[nd1][0];
        };
	// Y- plane
	if(y!=0){
          CurFF[0][0][y][x][1]-=crt[nd1][3];
        };
	// Y+ plane
	if(y!=cmeshy-1){
	  CurFF[0][0][y+1][x][1]+=crt[nd1][1];
        };
	if(Vacuum){
          real tmp=crt[nd1][4];
  	  // X- boundary
	  if(x==0){
            int tt=1;
	    if(y==0)tt++;
	    if(y==cmeshy-1)tt++;
	    if(x==cmeshx-1)tt++;
            real ratio=1./real(tt);
     	    CurFF[0][0][y][x][0]-=tmp*ratio;
	  };
	  // X+ boundary
	  if(x==cmeshx-1){
            int tt=1;
	    if(y==0)tt++;
	    if(y==cmeshy-1)tt++;
	    if(x==0)tt++;
            real ratio=1./real(tt);
	    CurFF[0][0][y][x+1][0]+=tmp*ratio;
	  };
	  // Y- boundary
	  if(y==0){
            int tt=1;
	    if(x==0)tt++;
	    if(x==cmeshx-1)tt++;
	    if(y==cmeshy-1)tt++;
            real ratio=1./real(tt);
	    CurFF[0][0][y][x][1]-=tmp*ratio;
	  };
          // Y+ boundary
	  if(y==cmeshy-1){
            int tt=1;
	    if(x==0)tt++;
	    if(x==cmeshx-1)tt++;
	    if(y==0)tt++;
            real ratio=1./real(tt);
	    CurFF[0][0][y+1][x][1]+=tmp*ratio;
	  };
	};

      };
    };
    //CalCoarseCur();
  };
  // ----------------------------------------------------------

  delete [] xs;
  delete [] selfscat;
  delete [] src;
  delete [] flux;
  delete [] ssc;
  delete [] flxd;


  return maxerr;
};

real MECSystem::CalFluxAngle360(int ng,real cin,int iter)
{
  bool cmr=false;

  // +++ CMR setting
  vector< vector<real> > crt_x;
  vector< vector<real> > crt_y;
  vector<real> srcsum;
  vector<real> sscsum;
  vector<real> stsum;
  GroupData2D mat1;
  GroupData1D vec1;

  if(cmr){
  cmeshid.resize(TotM);
  /*
  cmesh=3;
  cmeshid[0]=0;
  for(int i=0;i<4;i++){cmeshid[1+i]=1;};
  for(int i=0;i<4;i++){cmeshid[5+i]=2;};
  for(int i=0;i<12;i++){cmeshid[9+i]=0;};
  for(int i=0;i<28;i++){cmeshid[21+i]=1;};
  for(int i=0;i<28;i++){cmeshid[49+i]=2;};
  */

  /*
  cmesh=TotM;
  for(int i=0;i<TotM;i++){cmeshid[i]=i;};
  */
  /*
  cmesh=1;
  for(int i=0;i<TotM;i++){cmeshid[i]=0;};
  */
  /*
  cmesh=2;
  cmeshid[0]=0;
  cmeshid[1]=1;
  cmeshid[2]=1;
  cmeshid[3]=1;
  cmeshid[4]=1;
  cmeshid[5]=1;
  cmeshid[6]=1;
  cmeshid[7]=1;
  */
  /*
  cmesh=9;
  for(int i=0;i<9;i++){cmeshid[i]=i;};
  for(int i=0;i<12;i++){cmeshid[9+i]=0;};
  for(int i=0;i<8;i++){
    for(int j=0;j<7;j++){
      cmeshid[21+i*7+j]=i+1;
    };
  };
  */

  crt_x.resize(cmesh);
  crt_y.resize(cmesh);
  for(int i=0;i<cmesh;i++){
    crt_x[i].resize(cmesh,0.);
    crt_y[i].resize(cmesh,0.);
  };

  srcsum.resize(cmesh,0.);
  sscsum.resize(cmesh,0.);
  stsum.resize(cmesh,0.);

  mat1.put_yx(cmesh,cmesh);
  vec1.put_imax(cmesh);

  };


  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM*(pl+1)];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src_out=new real[TotM*plnum];  // Out-group source
  real *flux=new real[TotM*plnum]; // Flux moment
  real *fluxold=new real[TotM];

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source to `src_out' (angluar moment-wise)
  PutSourceInnerIteration(src_out);

  // +++ (sero-source check for GPT calculations) ++++++++++

  // In some specific cases, no scattering sources are given
  // in specific energy groups.
  // This is implemented in 2019/11/21.
  bool zero_src=true;
  for(int i=0;i<TotM;i++){
    if(src_out[i]>0){
      zero_src=false;
      break;
    };
  };
  if(zero_src)return 0.;
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /*
  for(int i=0;i<TotM;i++){
    cout<<i<<" "<<xs[i]<<" "<<selfscat[i]<<" "<<flux[i]<<" "<<src_out[i]<<"\n";
  };
  exit(0);
  */

  int aad=azimuthal_angle_division;
  int pad=polar_angle_division;
  int tot_ang=pad*aad*2;
  vector< vector<real> > aveflx(tot_ang,(vector<real>)TotM);
  // volume-integrated mesh-wise angular flux

  int innermax=1000;
  //real convergence_condition_influx=cin*0.05;

  // At the initial outer iteration, incoming angular flux in g-th group
  // is assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&!Vacuum){
    int ng2=ng-1;
    if(!opt.Forward())ng2=ng+1;
    if(ng2>=0&&ng2<ng){
      for(int i=0;i<polar_angle_division;i++){
        for(int j=0;j<numt*2;j++){
	  influx[i][j][ng]=influx[i][j][ng-1];
        };
      };
    };
  };

  real *putsrc=new real[max_segments]; // Source term
  real *ave=new real[max_segments];    // Mesh-averaged flux
  real *ave2=new real[max_segments];    // Mesh-averaged flux
  int *reg_id=new int[max_segments];

  // +++ CMR
  /*
  real *outflx1=new real[max_segments];
  real *outflx2=new real[max_segments];
  */

  real consttmp=aad*INV_PI;
  for(int inner=0;inner<innermax;inner++){

    total_inner_iteration[ng]++;

    // +++ CMR
    /*
    if(cmr){
    for(int i=0;i<cmesh;i++){
      for(int j=0;j<cmesh;j++){
	crt_x[i][j]=0.;
	crt_y[i][j]=0.;
      };
    };
    };
    */

    for(int j=0;j<tot_ang;j++){
      for(int i=0;i<TotM;i++){
        aveflx[j][i]=0.;
      };
    };


    int azm_angle_id=0;
    for(int tj=0;tj<numt;tj++){

      int r=tset->GetTrajectory(tj).GetNumreg();
      real width=tset->GetTrajectory(tj).GetWeight()*consttmp; // ray-width
      int next_id=tset->GetNextTrajectory(tj);
      int pre_id=tset->GetPreTrajectory(tj);

      for(int i=0;i<r;i++){
        reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
        // +++ source calculation 
	int id=reg_id[i];
	putsrc[i]=src_out[id]+selfscat[id]*flux[id];
      };

      for(int pp=0;pp<pad;pp++){ // polar-angle division

        int sn=pp*aad*2+azm_angle_id;
        int sn2=sn+aad;

        real sin_inv=polar_angle_sin_inv[pp];


        real outflx,outflx2;
        tset->GetTrajectory(tj).SolveCFormTransportDoubleDirection(fmoc,influx[pp][tj][ng],influx[pp][numt+tj][ng],putsrc,ave,ave2,outflx,outflx2,sin_inv);
        influx[pp][next_id][ng]=outflx;
	if(pre_id>numt)pre_id-=numt*2;
        influx[pp][numt+pre_id][ng]=outflx2;

	// +++ CMR
	/*
        tset->GetTrajectory(tj).SolveCFormTransportDoubleDirection(fmoc,influx[pp][tj][ng],influx[pp][numt+tj][ng],putsrc,ave,ave2,outflx1,outflx2,sin_inv);
        influx[pp][next_id][ng]=outflx1[r-1];
	if(pre_id>numt)pre_id-=numt*2;
        influx[pp][numt+pre_id][ng]=outflx2[r-1];
	*/

        for(int i=0;i<r;i++){
          int id=reg_id[i];
	  aveflx[sn][id]+=ave[i]*width;
	  aveflx[sn2][id]+=ave2[i]*width;
        };

	// +++ CMR
	/*
	if(cmr){
        real omega=quad.GetOmega(sn);
        real mu=quad.GetMu(sn);
        real et=quad.GetEata(sn);
        real sqq=sqrt(mu*mu+et*et);
        real cosine_inv=sqq/mu;
        real sine_inv=sqq/et;
        real factor_x=omega*width*cosine_inv*mu;
        real factor_y=omega*width*sine_inv*et;
        int id1=reg_id[0];
	int id2;
	for(int i=0;i<r-1;i++){
	  id2=reg_id[i+1];
	  int cd1=cmeshid[id1];
	  int cd2=cmeshid[id2];
	  if(cd1!=cd2){
	    crt_x[cd1][cd2]+=outflx1[i]*factor_x;
	    crt_y[cd1][cd2]+=outflx1[i]*factor_y;
	    crt_x[cd2][cd1]+=outflx2[r-2-i]*factor_x;
	    crt_y[cd2][cd1]+=outflx2[r-2-i]*factor_y;
	  };
	  id1=id2;
	};
	};
	*/

      };

      if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id))azm_angle_id++;

    };

    // +++ Renew of angular moment of mesh-wise flux
    for(int i=0;i<TotM;i++){
      real tmp=0.;
      for(int j=0;j<tot_ang;j++){
        tmp+=aveflx[j][i]*quad.GetOmega(j);
      };
      tmp/=mesh[i].GetVolume();
      fluxold[i]=flux[i];
      flux[i]=tmp;
    };

    // +++ CMR acceleration

    /*
    for(int i=0;i<cmesh;i++){
      for(int j=0;j<cmesh;j++){
	cout<<i<<" "<<j<<" "<<crt_x[i][j]<<" "<<crt_y[i][j]<<"\n";
      };
    }; exit(0);
    */

    if(cmr){

      real gcmr_factor=0.5;
      //vector<real> gcmr_factor(cmesh);
      //gcmr_factor[0]=0.5;
      //gcmr_factor[1]=0.5;
      //gcmr_factor[2]=0.5;
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<cmesh;j++){
	  if(i!=j){
	    crt_x[i][j]=0.5*((1+2*gcmr_factor)*crt_y[i][j]-(1-2*gcmr_factor)*crt_y[j][i]);
	    //real f=gcmr_factor[j]/(gcmr_factor[i]+gcmr_factor[j]);
	    //crt_x[i][j]=f*((1+2*gcmr_factor[i])*crt_y[i][j]-(1-2*gcmr_factor[i])*crt_y[j][i]);
	  };
        };
      };

    for(int i=0;i<cmesh;i++){
      srcsum[i]=0.;
      sscsum[i]=0.;
      stsum[i]=0.;
    };

    for(int i=0;i<TotM;i++){
      int cid=cmeshid[i];
      real vol=mesh[i].GetVolume();
      real tmp=flux[i]*vol;
      srcsum[cid]+=src_out[i]*vol;
      sscsum[cid]+=tmp*selfscat[i];
      stsum[cid]+=xs[i]*tmp;
    };

    mat1.set_zero();
    for(int i=0;i<cmesh;i++){
      real tmp=stsum[i]-sscsum[i]*PI4;
      for(int j=0;j<cmesh;j++){
	//cout<<i<<" "<<j<<" "<<crt_x[i][j]<<" "<<crt_y[i][j]<<"\n";
	tmp+=crt_x[i][j];
	mat1.put_data(i,j,-crt_x[j][i]);
      };
      mat1.put_data(i,i,tmp);
      vec1.put_data(i,srcsum[i]*PI4);
    };
    //mat1.show_self(); exit(0);
    mat1.solveaxb_mod(vec1);

    for(int i=0;i<TotM;i++){
      flux[i]*=vec1.get_dat(cmeshid[i]); // rebalance factor
    };
    for(int i=0;i<polar_angle_division;i++){
      for(int j=0;j<numt;j++){
        int r=tset->GetTrajectory(j).GetNumreg();
        int regids=tset->GetTrajectory(j).GetReg(0);
        int regide=tset->GetTrajectory(j).GetReg(r-1);
        //influx[i][j][ng]*=vec1.get_dat(cmeshid[regids]);
        //influx[i][j+numt][ng]*=vec1.get_dat(cmeshid[regide]);
	int next_id=tset->GetNextTrajectory(j);
	int pre_id=tset->GetPreTrajectory(j);
        influx[i][next_id][ng]*=vec1.get_dat(cmeshid[regide]);
	if(pre_id>numt)pre_id-=numt*2; // (for reflective boundary)
        influx[i][numt+pre_id][ng]*=vec1.get_dat(cmeshid[regids]);

      };
    };


    };

    /*
    cout.setf(ios::showpoint);
    cout.precision(10);
    for(int i=0;i<cmesh;i++){
      cout<<vec1.get_dat(i)<<" ";
    };
    cout<<"\n";
    */
    /*
    for(int i=0;i<cmesh;i++){
      if(vec1.get_dat(i)<0.){
	cout<<"# Rebalance factor is negative!\n";
      };
    };
    */

    // ---- system-rebalancing -----------------------------------

    if(!cmr&&!Vacuum){
      real srcsum=0.;
      real colsum=0.;
      for(int i=0;i<TotM;i++){
        real vol=mesh[i].GetVolume();
        srcsum+=src_out[i]*vol*PI4;
	colsum+=flux[i]*vol*(xs[i]-selfscat[i]*PI4);
      };
      real factor=srcsum/colsum;
      for(int i=0;i<TotM;i++){
        flux[i]*=factor; // rebalance factor
      };
      for(int i=0;i<polar_angle_division;i++){
        for(int j=0;j<numt;j++){
          influx[i][j][ng]*=factor;
        };
      };
    };

    // ------------------------------------------------------------

    real errmax=0.;
    for(int i=0;i<TotM;i++){
      real err=fabs(flux[i]/fluxold[i]-1.);
      if(err>errmax)errmax=err;
    };
    //if(!opt.Forward())cout<<"   (grp:"<<ng<<") "<<inner<<" "<<errmax<<"\n";

    bool loop_fin=false;
    if(inner==innermax-1){
      cout<<"# Inner iteration does not converge in group "<<ng<<".\n";
      cout<<"# Iteration residual is "<<errmax<<"\n";
      loop_fin=true;
    };
    if(errmax<cin){
      loop_fin=true;
      //cout<<"   (grp:"<<ng<<") "<<inner<<" "<<errmax<<"\n";
    };
    if(loop_fin)break;

  };// --- inner iteration loop-end

  // +++ flux-moment renew & error calculation of scalar flux
  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    mesh[i].GetFlux().put_data(ng,flux[i]);
  };

  if(wrtflx){
    for(int i=0;i<TotM;i++){
      real vol_inv=1./mesh[i].GetVolume();
      for(int j=0;j<tot_ang;j++){
	aflux[i][ng].put_data(j,aveflx[j][i]*vol_inv);
      };
    };
  };

  delete [] reg_id;
  delete [] xs;
  delete [] selfscat;
  delete [] src_out;
  delete [] flux;
  delete [] fluxold;
  delete [] putsrc;
  delete [] ave;
  delete [] ave2;

  // +++ CMR
  /*
  delete [] outflx1;
  delete [] outflx2;
  */

  return maxerr;
};

real MECSystem::CalFluxAngle360Old(int ng,real cin,int iter,bool angular_dependent_source)
{
  if(angular_dependent_source)return CalFluxAngle360AngularDependentSource(ng,cin,iter);

  // +++ 360-degree calculation +++
  //
  //  - Limited for P0 scattering only
  //  - Currents along x- and y-directions are stored.
  //
  // CMR/CMFD acceleration is being developed (not successfully worked).
  // Detail can be found in the notebook
  //   2008/10/29, 11/10, 12/16

  //bool cmfd_on=false;
  //if(cmfd&&iter%opt.GetItcmfd()==0)cmfd_on=true;

  if(Vacuum){
    cout<<"# Error in MECSystem::CalFluxAngle360.\n";
    cout<<"# This routine can be applied only to periodic BCs.\n";
    exit(0);
  };

  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM*(pl+1)];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src_out=new real[TotM*plnum];  // Out-group source
  real *flux=new real[TotM*3]; // Flux moment (P1 components are stored)
  //real *flux=new real[TotM*plnum]; // Flux moment (P1 components are stored)

  // +++ Leakage data storing : leak[cmesh1][cmesh2] (cmesh1->cmesh2)
  /*
  vector< vector<real> > leak;
  if(cmr||cmfd_on){
    leak.resize(cmesh);
    for(int i=0;i<cmesh;i++){
      leak[i].resize(cmesh);
    };
  };
  */

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source to `src_out' (angluar moment-wise)
  PutSourceInnerIteration(src_out);

  int aad=azimuthal_angle_division;
  int pad=polar_angle_division;
  int tot_ang=pad*aad*2;
  vector< vector<real> > aveflx(tot_ang,(vector<real>)TotM);
  // volume-integrated mesh-wise angular flux

  int innermax=1000;
  //real convergence_condition_influx=cin*0.05;

  // At the initial outer iteration, incoming angular flux in g-th group
  // is assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&!Vacuum){
    int ng2=ng-1;
    if(!opt.Forward())ng2=ng+1;
    if(ng2>=0&&ng2<ng){
      for(int i=0;i<polar_angle_division;i++){
        for(int j=0;j<numt;j++){
	  influx[i][j][ng]=influx[i][j][ng-1];
        };
      };
    };
  };

  /*
  vector< vector<real> > escape_probability(polar_angle_division*numt);
  for(int i=0;i<numt;i++){
    int r=tset->GetTrajectory(i).GetNumreg();
    for(int j=0;j<polar_angle_division;j++){
      escape_probability[j*numt+i].resize(r,0.);
    };
    for(int j=0;j<r;j++){
      real opt=tset->GetTrajectory(i).GetOpt(j);
      real len=tset->GetTrajectory(i).GetDist(j);
      for(int k=0;k<polar_angle_division;k++){
	real sin_inv=polar_angle_sin_inv[k];
	escape_probability[k*numt+i][j]=len*sin_inv*fmoc.get(opt*sin_inv);
      };
    };
  };
  */

  real *putsrc=new real[max_segments]; // Source term
  real *ave=new real[max_segments];    // Mesh-averaged flux
  vector<int> reg_id(max_segments);

  real consttmp=aad*INV_PI;
  for(int inner=0;inner<innermax;inner++){

    for(int j=0;j<tot_ang;j++){
      for(int i=0;i<TotM;i++){
        aveflx[j][i]=0.;
      };
    };

    // +++ initialize leakage-matrix
    /*
    if(cmr||cmfd_on){
      for(int i=0;i<cmesh;i++){
	for(int j=0;j<cmesh;j++){
	  leak[i][j]=0.;
	};
      };
    };
    */

    //real maxerr_bflx=0.;
    int sn=0;
    //real omega=quad.GetOmega(sn);
    //real mu=quad.GetMu(sn);
    //real et=quad.GetEata(sn);
    //real sqq=sqrt(mu*mu+et*et); 

    real sin_inv=1.;
    //real sin_inv_old=1.;
    for(int pp=0;pp<pad;pp++){ // polar-angle division
 
      //sin_inv_old=sin_inv;
      sin_inv=polar_angle_sin_inv[pp];
      
      int azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){

        int sn2=sn+aad;
        int r=tset->GetTrajectory(tj).GetNumreg();
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real width=tset->GetTrajectory(tj).GetWeight()*consttmp; // ray-width

	real *esc=new real[r];    // Mesh-escaping flux (0-180 degree)
	real *esc2=new real[r];   // Mesh-escaping flux (180-360 degree)
	//real outflx,outflx2;
        //real *ep=new real[r];     // Escape-probability

	// +++ source calculation 
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=src_out[id]+selfscat[id]*flux[id];
	  putsrc[i]=sum;
	  //ep[i]=escape_probability[pp*numt+tj][i];
        };

	// +++ 0-180 
	/*
        tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,ave,esc,sin_inv);
        real outflx=esc[r-1];
	*/
        real outflx;
        tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,ave,outflx,sin_inv);

	int next_id=tset->GetNextTrajectory(tj);
  	//real err=fabs(influx[pp][next_id][ng]/outflx-1.);
	//if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][next_id][ng]=outflx;
        for(int i=0;i<r;i++){
	  aveflx[sn][reg_id[i]]+=ave[i]*width;
        };
	// +++ 180-360
	/*
        tset->GetTrajectory(tj).SolveCFormTransportOpposite(fmoc,influx[pp][numt+tj][ng],putsrc,ave,esc2,sin_inv);
        real outflx2=esc2[0];
	*/
        real outflx2;
        tset->GetTrajectory(tj).SolveCFormTransportOpposite(fmoc,influx[pp][numt+tj][ng],putsrc,ave,outflx2,sin_inv);
	int next_id2=tset->GetPreTrajectory(tj);
        int pre_id=next_id2;
	if(next_id2>numt)pre_id-=numt*2; // (for reflective boundary)
  	//err=fabs(influx[pp][numt+pre_id][ng]/outflx2-1.);
	//if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][numt+pre_id][ng]=outflx2;
        for(int i=0;i<r;i++){
	  aveflx[sn2][reg_id[i]]+=ave[i]*width;
        };

	// ------ leakage calculation --------------------------------------
	/*
	if(cmr||cmfd_on){
	  real factor=omega*width*sqq;
	  for(int i=0;i<r-1;i++){
	    int id1=reg_id[i];
	    int id2=reg_id[i+1];
	    int cd1=cmeshid[id1];
	    int cd2=cmeshid[id2];
	    if(cd1!=cd2){
	      leak[cd1][cd2]+=esc[i]*factor;
	      leak[cd2][cd1]+=esc2[i+1]*factor;
	    };
	  };
	  //
          int next_reg;
	  if(next_id<numt){
	    next_reg=tset->GetTrajectory(next_id).GetReg(0);
	  }else{
	    next_reg=tset->GetTrajectory(next_id-numt).GetLastReg();
	  };
	  int cd1=cmeshid[reg_id[r-1]];
	  int cd2=cmeshid[next_reg];
	  if(cd1!=cd2)leak[cd1][cd2]+=esc[r-1]*factor;
	  //
	  int pre_reg;
	  if(next_id2<numt){
	    pre_reg=tset->GetTrajectory(next_id2).GetLastReg();
	  }else{
	    pre_reg=tset->GetTrajectory(next_id2-numt).GetReg(0);
	  };
	  cd1=cmeshid[reg_id[0]];
	  cd2=cmeshid[pre_reg];
	  if(cd1!=cd2)leak[cd1][cd2]+=esc2[0]*factor;
	};
	*/
	// -----------------------------------------------------------------

	delete [] esc;
	delete [] esc2;
        //delete [] ep;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Sweep-end for azimuthal angle
	  sn++;
          azm_angle_id++;
	  if(tj==numt-1)sn+=aad; // only for 360-degree treatment
	  if(sn!=tot_ang){ // Set-up for next azimuthal angle calculations
            //omega=quad.GetOmega(sn);
            //real mu=quad.GetMu(sn);
            //real et=quad.GetEata(sn);
            //sqq=sqrt(mu*mu+et*et);
	  };
	};

      };

    };

    // +++ Renew of angular moment of mesh-wise flux
    real errmax=0.;
    {
      int index=0;
      for(int i=0;i<TotM;i++){
        real tmp=0.;
        for(int j=0;j<tot_ang;j++){
          tmp+=aveflx[j][i]*quad.GetOmega(j);
        };
        tmp/=mesh[i].GetVolume();
	real err=fabs(tmp/flux[index]-1.);
	if(err>errmax)errmax=err;
        flux[index]=tmp;
	index++;
      };
    };

    bool loop_fin=false;
    if(inner==innermax-1){
      cout<<"# Inner iteration does not converge in group "<<ng<<".\n";
      cout<<"# Iteration residual is "<<errmax<<"\n";
      loop_fin=true;
    };

    if(errmax<cin){
      loop_fin=true;
      //cout<<"   (grp:"<<ng<<") "<<inner<<" "<<errmax<<"\n";
      //  break;
    };

    // ---- system-rebalancing -----------------------------------
    //     (valid only for zero-leakage B.C.)

    if(!cmr&&!Vacuum&&!loop_fin){
      //if(!cmr&&!Vacuum){
      real srcsum=0.;
      real sscsum=0.;
      real stsum=0.;
      for(int i=0;i<TotM;i++){
        real vol=mesh[i].GetVolume();
        srcsum+=src_out[i]*vol;
        sscsum+=flux[i]*selfscat[i]*vol;
        stsum+=xs[i]*flux[i]*vol;
      };
      real factor=(srcsum*PI4)/(stsum-sscsum*PI4);
      cout.setf(ios::showpoint);
      cout.precision(6);
      //cout<<"Rebalance factor  : "<<factor<<"\n";
      for(int i=0;i<TotM;i++){
        flux[i]*=factor; // rebalance factor
	/*
        flux[TotM+i]*=factor; // rebalance factor
        flux[TotM*2+i]*=factor; // rebalance factor
	*/
      };
      for(int i=0;i<polar_angle_division;i++){
        for(int j=0;j<numt;j++){
          influx[i][j][ng]*=factor;
        };
      };
    };


    // ------------------------------------------------------------

    // --- Coarse Mesh Rebalance Acceleration -------------------

    if(cmr){
      /*
      vector<real> lhs_diag(cmesh);
      vector<real> rhs(cmesh);

      int stt=0;
      for(int i=0;i<cmesh;i++){
        real srcsum=0.;
        real sscsum=0.;
        real stsum=0.;
        for(int j=stt;j<=cmesh_bound[i];j++){
  	  real vol=mesh[j].GetVolume();
          real tmp=flux[j*plnum]*vol;
	  srcsum+=src[j*plnum]*vol;
	  sscsum+=tmp*selfscat[j*(pl+1)];
	  stsum+=xs[j]*tmp;
        };
        srcsum*=PI4;
        sscsum*=PI4;
	real tmpdiag=stsum-sscsum;
	for(int j=0;j<cmesh;j++){
	  tmpdiag+=leak[i][j];
	};
        lhs_diag[i]=1./tmpdiag;
        rhs[i]=srcsum;
	stt=cmesh_bound[i]+1;
      };

      vector<real> sol(cmesh);
      for(int i=0;i<cmesh;i++){sol[i]=rhs[i];};

      int itermax=1000000;
      real eps=cin*0.1;
      //real eps=1e-9;

      for(int iter=0;iter<itermax;iter++){
	real errmax=0.;
	for(int i=0;i<cmesh;i++){
	  real tmp=rhs[i];
	  for(int j=0;j<cmesh;j++){
	    tmp+=leak[j][i]*sol[j];
	  };
	  tmp*=lhs_diag[i];
	  real err=fabs(tmp/sol[i]-1.);
	  if(err>errmax)errmax=err;
	  sol[i]=tmp;
	};
	if(errmax<eps){
          break;
	};
      };

      stt=0;
      for(int i=0;i<cmesh;i++){
	//cout.setf(ios::showpoint);
	//cout.precision(6);
	//cout<<"Rebalance factor "<<i<<" : "<<sol[i]<<"\n";
        for(int j=stt;j<=cmesh_bound[i];j++){
	  for(int l=0;l<plnum;l++){
  	    flux[j*plnum+l]*=sol[i];
	  };
        };
        stt=cmesh_bound[i]+1;
      };
      //cout<<"\n";


      for(int j=0;j<numt;j++){
	int id=tset->GetPreTrajectory(j);
        int idd;
	if(id<numt){
    	  idd=tset->GetTrajectory(id).GetLastReg();
	}else{
	  idd=tset->GetTrajectory(id-numt).GetReg(0);
	};
	int id2=tset->GetNextTrajectory(j);
	int idd2;
	if(id2<numt){
	  idd2=tset->GetTrajectory(id2).GetReg(0);
	}else{
	  idd2=tset->GetTrajectory(id2-numt).GetLastReg();
	};
        for(int i=0;i<pad;i++){
	  influx[i][j][ng]*=sol[cmeshid[idd]];
	  influx[i][j+numt][ng]*=sol[cmeshid[idd2]];
	};
      };
      */

      /*
      for(int i=0;i<cmesh;i++){
        for(int j=0;j<cmesh;j++){
	  leak[i][j]*=sol[i];
        };
      };
      */

    };
    // ---------------------------------------------
    // +++ end of CMR Acceleration


    if(loop_fin){

      // (current calculation)
      for(int i=0;i<TotM;i++){
        real tmp_x=0.;
        real tmp_y=0.;
        for(int j=0;j<tot_ang;j++){
          real flom=aveflx[j][i]*quad.GetOmega(j);
	  tmp_x+=flom*quad.GetMu(j);
	  tmp_y+=flom*quad.GetEata(j);
	  //tmp_x+=flom*quad.GetMoment(1,j);
	  //tmp_y+=flom*quad.GetMoment(2,j);
        };
        real vinv=1./mesh[i].GetVolume();
	flux[TotM+i]=tmp_x*vinv;
	flux[TotM*2+i]=tmp_y*vinv;
      };

      break;
    };


  };// --- inner iteration loop-end

  // +++ flux-moment renew & error calculation of scalar flux
  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    mesh[i].GetFlux().put_data(ng,flux[i]);

    mesh[i].GetFlux(1).put_data(ng,flux[TotM+i]); // x-current
    mesh[i].GetFlux(2).put_data(ng,flux[TotM*2+i]); // y-current

  };

  if(wrtflx){
    for(int i=0;i<TotM;i++){
      real vol_inv=1./mesh[i].GetVolume();
      for(int j=0;j<tot_ang;j++){
	aflux[i][ng].put_data(j,aveflx[j][i]*vol_inv);
      };
    };
  };

  delete [] xs;
  delete [] selfscat;
  delete [] src_out;
  delete [] flux;
  delete [] putsrc;
  delete [] ave;

  return maxerr;
};

real MECSystem::CalFluxAngle360AngularDependentSource(int ng,real cin,int iter,bool benoist_consistent)
{
  // Only in-coming angular flux residual is watched.
  // (Residual of mesh-averaged flux is ignored)

  // (Limited for P0 scattering)

  if(Vacuum){
    cout<<"Error in MECSystem::CalFluxAngle360.\n";
    cout<<"It is not yet coded for vacuum boundary conditions.\n";
    exit(0);
  };

  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM*(pl+1)];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src=new real[TotM*plnum];  // External source
  real *flux=new real[TotM*plnum]; // Flux moment
  real *ssc=new real[TotM*plnum];  // Self-scattering source

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source
  PutSourceInnerIteration(src);

  int tot_ang=polar_angle_division*azimuthal_angle_division*2;
  vector< vector<real> > aveflx(tot_ang); // mesh-averaged angular flux
  for(int i=0;i<tot_ang;i++){
    aveflx[i].resize(TotM,0.);
  };

  real *nume=new real[TotM];

  int innermax=1000;
  real convergence_condition=cin*0.05;
  //real convergence_condition=1e-7;

  // At the initial outer iteration,
  // incoming angular flux in g-th group is 
  // assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&ng!=0){
    int ttr=numt*2;
    for(int i=0;i<polar_angle_division;i++){
      for(int j=0;j<ttr;j++){
	influx[i][j][ng]=influx[i][j][ng-1];
      };
    };
  };

  real consttmp2=PI/azimuthal_angle_division;
  for(int inner=0;inner<innermax;inner++){

    // +++ self-scattering source
    for(int i=0;i<TotM;i++){
      ssc[i]=selfscat[i]*flux[i];
      if(benoist_consistent)ssc[i]=0.; // For Benoist's DC-consistent calculation
    };

    real maxerr_bflx=0.;
    int sn=0;

    for(int pp=0;pp<polar_angle_division;pp++){

      real sin_inv=polar_angle_sin_inv[pp];

      // +++ 0 - 180 degree
      for(int i=0;i<TotM;i++){
	nume[i]=0.;
      };
      int azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){
        int r=tset->GetTrajectory(tj).GetNumreg();
	vector<int> reg_id(r);
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real weight=tset->GetTrajectory(tj).GetWeight();

	real *putsrc=new real[r]; // Source term
	real *tmp=new real[r];    // Mesh-averaged flux
	real outflx;

	// +++ source calculation 
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=src[id]+ssc[id];
	  sum+=aflux[id][ng].get_dat(sn); // angular_dependence_source
	  putsrc[i]=sum;
        };
        tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,tmp,outflx,sin_inv);
	int next_id=tset->GetNextTrajectory(tj);
  	real err=fabs(influx[pp][next_id][ng]/outflx-1.);
	if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][next_id][ng]=outflx;
        for(int i=0;i<r;i++){
  	  int id=reg_id[i];
	  nume[id]+=tmp[i]*weight;
        };
	delete [] putsrc;
	delete [] tmp;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Sweep-end for azimuthal angle
          for(int i=0;i<TotM;i++){
	    aveflx[sn][i]=nume[i]/(consttmp2*mesh[i].GetVolume());
          };
	  sn++;
	  if(tj!=numt-1){ // Initialization for next azimuthal angle calculations
            azm_angle_id++;
            for(int i=0;i<TotM;i++){
	      nume[i]=0.;
	    };
          };
	};
      };

      // +++ 180 - 360 degree (opposite direction)
      for(int i=0;i<TotM;i++){
	nume[i]=0.;
      };
      azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){
        int r=tset->GetTrajectory(tj).GetNumreg();
	vector<int> reg_id(r);
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real weight=tset->GetTrajectory(tj).GetWeight();

	real *putsrc=new real[r]; // Source term
	real *tmp=new real[r];    // Mesh-averaged flux

	real outflx;

	// +++ source calculation
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=src[id]+ssc[id];
	  sum+=aflux[id][ng].get_dat(sn);
	  putsrc[i]=sum;
        };
        tset->GetTrajectory(tj).SolveCFormTransportOpposite(fmoc,influx[pp][numt+tj][ng],putsrc,tmp,outflx,sin_inv);
        int next_id=tset->GetPreTrajectory(tj);
  	real err=fabs(influx[pp][numt+next_id][ng]/outflx-1.);
	if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][numt+next_id][ng]=outflx;
        for(int i=0;i<r;i++){
  	  int id=reg_id[i];
	  nume[id]+=tmp[i]*weight;
        };
	delete [] putsrc;
	delete [] tmp;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Sweep-end for azimuthal angle
          for(int i=0;i<TotM;i++){
	    aveflx[sn][i]=nume[i]/(consttmp2*mesh[i].GetVolume());
          };
	  sn++;
	  if(tj!=numt-1){ // Initialization for next azimuthal angle calculations
            azm_angle_id++;
            for(int i=0;i<TotM;i++){
	      nume[i]=0.;
	    };
          };
	};

      };
    };

    // +++ Renew of angular moment of mesh-wise flux
    {
      int index=0;
      for(int i=0;i<TotM;i++){
        real tmp=0.;
        for(int j=0;j<tot_ang;j++){
          tmp+=aveflx[j][i]*quad.GetOmega(j);
        };
        flux[index++]=tmp;
      };
    };

    if(inner==innermax-1){
      cout<<" Inner iteration does not converge in group "<<ng<<".\n";
    };
    if(maxerr_bflx<convergence_condition){
      //cout<<"   (grp:"<<ng<<") "<<inner<<"\n";
      break;
    };
  }; 

  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i*plnum]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    mesh[i].GetFlux().put_data(ng,flux[i]);
  };

  if(wrtflx){
    for(int i=0;i<TotM;i++){
      for(int j=0;j<tot_ang;j++){
	aflux[i][ng].put_data(j,aveflx[j][i]);
      };
    };
  };

  delete [] xs;
  delete [] selfscat;
  delete [] src;
  delete [] flux;
  delete [] ssc;
  delete [] nume;

  return maxerr;
};

real MECSystem::CalFluxDegree720(int ng,real cin,int iter)
{
  // Only in-coming angular flux residual is watched.
  // (Residual of mesh-averaged flux is ignored)

  if(Vacuum){
    cout<<"Error in MECSystem::CalFluxDegree720.\n";
    cout<<"It is not yet coded for vacuum boundary conditions.\n";
    exit(0);
  };

  // +++ Put sigma_t and sigma_s to temporal arrays
  real *xs=new real[TotM];
  real *selfscat=new real[TotM];
  PutSigmaForInnerIteration(ng,xs,selfscat);
  tset->PutXS(xs);

  int numt=tset->GetNumTrajectory();

  real *src=new real[TotM];  // External source
  real *flux=new real[TotM]; // Flux moment
  real *ssc=new real[TotM];  // Self-scattering source

  // +++ Initial Flux guess
  InitialFluxGuessInnerIteration(ng,flux);

  // +++ Put fission & slowing down source
  PutSourceInnerIteration(src);

  int tot_ang=polar_angle_division*2*azimuthal_angle_division*2;
  vector< vector<real> > aveflx(tot_ang); // mesh-averaged angular flux
  for(int i=0;i<tot_ang;i++){
    aveflx[i].resize(TotM,0.);
  };

  real *nume=new real[TotM];

  // At the initial outer iteration,
  // incoming angular flux in g-th group is 
  // assumed to be same as incoming flux in (g-1)-th group
  if(iter==0&&ng!=0&&!Vacuum){
    int ttr=numt*2;
    for(int i=0;i<polar_angle_division*2;i++){
      for(int j=0;j<ttr;j++){
	influx[i][j][ng]=influx[i][j][ng-1];
      };
    };
  };

  int innermax=2000;
  real convergence_condition=cin*0.05;
  //real convergence_condition=1e-7;

  real consttmp2=PI/azimuthal_angle_division;
  for(int inner=0;inner<innermax;inner++){

    // +++ self-scattering source
    for(int i=0;i<TotM;i++){
      ssc[i]=selfscat[i]*flux[i];
    };

    real maxerr_bflx=0.;
    int sn=0;

    for(int pp=0;pp<polar_angle_division*2;pp++){
      int ppp=pp;
      if(pp>=polar_angle_division)ppp-=polar_angle_division;
      real sin_inv=polar_angle_sin_inv[ppp];

      // +++ 0 - 180 degree
      for(int i=0;i<TotM;i++){
	nume[i]=0.;
      };
      int azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){
        int r=tset->GetTrajectory(tj).GetNumreg();
	vector<int> reg_id(r);
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real weight=tset->GetTrajectory(tj).GetWeight();

	real *putsrc=new real[r]; // Source term
	real *tmp=new real[r];    // Mesh-averaged flux
        real outflx;

	// +++ source calculation 
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=src[id]+ssc[id];
	  if(pp<polar_angle_division)sum+=aflux[id][ng].get_dat(sn);
	  putsrc[i]=sum;
        };
        //tset->GetTrajectory(tj).SolveCFormTransport(etab,influx[pp][tj][ng],putsrc,tmp,tmp2,sin_inv);
        //tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,tmp,tmp2,sin_inv);
        tset->GetTrajectory(tj).SolveCFormTransport(fmoc,influx[pp][tj][ng],putsrc,tmp,outflx,sin_inv);

	int next_id=tset->GetNextTrajectory(tj);
  	real err=fabs(influx[pp][next_id][ng]/outflx-1.);
	if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][next_id][ng]=outflx;

        for(int i=0;i<r;i++){
  	  int id=reg_id[i];
	  nume[id]+=tmp[i]*weight;
        };
	delete [] putsrc;
	delete [] tmp;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Sweep-end for azimuthal angle
          for(int i=0;i<TotM;i++){
	    aveflx[sn][i]=nume[i]/(consttmp2*mesh[i].GetVolume());
          };
	  sn++;
	  if(tj!=numt-1){ // Initialization for next azimuthal angle calculations
            azm_angle_id++;
            for(int i=0;i<TotM;i++){
	      nume[i]=0.;
	    };
          };
	};
      };

      // +++ 180 - 360 degree (opposite direction)
      for(int i=0;i<TotM;i++){
	nume[i]=0.;
      };
      azm_angle_id=0;
      for(int tj=0;tj<numt;tj++){
        int r=tset->GetTrajectory(tj).GetNumreg();
	vector<int> reg_id(r);
	for(int i=0;i<r;i++){
	  reg_id[i]=tset->GetTrajectory(tj).GetReg(i);
	};
        real weight=tset->GetTrajectory(tj).GetWeight();

	real *putsrc=new real[r]; // Source term
	real *tmp=new real[r];    // Mesh-averaged flux
        real outflx;

	// +++ source calculation
        for(int i=0;i<r;i++){
	  int id=reg_id[i];
	  real sum=src[id]+ssc[id];
	  if(pp<polar_angle_division)sum+=aflux[id][ng].get_dat(sn);
	  putsrc[i]=sum;
        };

        //tset->GetTrajectory(tj).SolveCFormTransportOpposite(etab,influx[pp][numt+tj][ng],putsrc,tmp,tmp2,sin_inv);
        //tset->GetTrajectory(tj).SolveCFormTransportOpposite(fmoc,influx[pp][numt+tj][ng],putsrc,tmp,tmp2,sin_inv);
        tset->GetTrajectory(tj).SolveCFormTransportOpposite(fmoc,influx[pp][numt+tj][ng],putsrc,tmp,outflx,sin_inv);

        int next_id=tset->GetPreTrajectory(tj);
  	real err=fabs(influx[pp][numt+next_id][ng]/outflx-1.);
	if(err>maxerr_bflx)maxerr_bflx=err;
        influx[pp][numt+next_id][ng]=outflx;

        for(int i=0;i<r;i++){
  	  int id=reg_id[i];
	  nume[id]+=tmp[i]*weight;
        };
	delete [] putsrc;
	delete [] tmp;

	if(tj==tset->GetTrajectoryOnAngularBoundary(azm_angle_id)){
	  // Sweep-end for azimuthal angle
          for(int i=0;i<TotM;i++){
	    aveflx[sn][i]=nume[i]/(consttmp2*mesh[i].GetVolume());
          };
	  sn++;
	  if(tj!=numt-1){ // Initialization for next azimuthal angle calculations
            azm_angle_id++;
            for(int i=0;i<TotM;i++){
	      nume[i]=0.;
	    };
          };
	};

      };
    };

    // +++ Renew of angular moment of mesh-wise flux
    {
      for(int i=0;i<TotM;i++){
        real tmp=0.;
        for(int j=0;j<tot_ang/2;j++){
          tmp+=(aveflx[j][i]+aveflx[tot_ang/2+j][i])*quad.GetOmega(j)*0.5;
        };
        flux[i]=tmp;
      };
    };

    if(inner==innermax-1){
      cout<<" Inner iteration does not converge in group "<<ng<<".\n";
    };
    if(maxerr_bflx<convergence_condition){
      //cout<<"   (grp:"<<ng<<") "<<inner<<"\n";
      break;
    };
  }; 

  real maxerr=0.;
  for(int i=0;i<TotM;i++){
    real err=fabs(flux[i]/mesh[i].GetFlux().get_dat(ng)-1);
    if(err>maxerr)maxerr=err;
    mesh[i].GetFlux().put_data(ng,flux[i]);
    //cout<<i<<" "<<flux[i*plnum+l]<<"\n";
  };

  if(wrtflx){
    for(int i=0;i<TotM;i++){
      for(int j=0;j<tot_ang/2;j++){
	aflux[i][ng].put_data(j,aveflx[j][i]-aveflx[tot_ang/2+j][i]);
      };
    };
  };

  delete [] xs;
  delete [] selfscat;
  delete [] src;
  delete [] flux;
  delete [] ssc;
  delete [] nume;

  return maxerr;
};

void MECSystem::PutSigmaForInnerIteration(int g,real *sigt,real *sigsself)
{
  int index=0;
  for(int i=0;i<TotM;i++){
    sigt[i]=mesh[i].GetMed()->GetMacxs().GetSigt().get_dat(g);
    for(int l=0;l<=pl;l++){
      sigsself[index++]=mesh[i].GetMed()->GetMacxs().GetSigs(l).get_dat(g,g)*INV_PI4;
    };
  };
};

void MECSystem::InitialFluxGuessInnerIteration(int g,real *flx)
{
  int index=0;
  for(int i=0;i<TotM;i++){
    for(int l=0;l<plnum;l++){
      flx[index++]=mesh[i].GetFlux(l).get_dat(g);
    };
  };
};

void MECSystem::PutSourceInnerIteration(real *src)
{
  int index=0;
  for(int i=0;i<TotM;i++){
    real vol=mesh[i].GetVolume();
    real volinv=0.;
    if(vol>0)volinv=1./vol;
    for(int l=0;l<plnum;l++){
      src[index++]=mesh[i].GetSrcin(l)*volinv;
    };
  };
};

void MECSystem::PutCMesh(int i,int j)
{
  cmeshx=i;
  cmeshy=j;
  cmesh=i*j;
  CTotM=cmesh;
  cmesh_bound.resize(cmesh);
  cmeshid.resize(TotM,-1);
  fmesh_par_cmesh.resize(cmesh);
  cmxl.resize(cmeshx);
  cmyl.resize(cmeshy);
  // for CMFD
  CurFF.resize(1);
  for(int g=0;g<1;g++){
    CurFF[g].resize(1); 
    for(int i=0;i<1;i++){
      CurFF[g][i].resize(cmeshy+1);
      for(int j=0;j<cmeshy+1;j++){
        CurFF[g][i][j].resize(cmeshx+1);
        for(int k=0;k<cmeshx+1;k++){
          CurFF[g][i][j][k].resize(Ndim,0.);
	};
      };
    };
  };

};

void MECSystem::PutCMeshLength(real *xin,real *yin)
{
  for(int i=0;i<cmeshx;i++){cmxl[i]=xin[i];};
  for(int i=0;i<cmeshy;i++){cmyl[i]=yin[i];};
};

void MECSystem::PutCMeshBound(int *inp)
{
  if(cmesh<1){
    cout<<"Error in MECSystem::PutCMeshBound.\n";
    cout<<"CMesh is less than 1. ("<<cmesh<<")\n";
    cout<<"Please do MECSystem::PutCMesh.\n";
    exit(0);
  };

  for(int i=0;i<cmesh;i++){
    cmesh_bound[i]=inp[i];
    if(i!=0){
      if(inp[i]<=inp[i-1]){
	cout<<"Error in MECSystem::PutCMeshBound.\n";
	cout<<"Please check your input data.\n";
	exit(0);
      };
    };
    if(i==0){
      fmesh_par_cmesh[i]=cmesh_bound[i]+1;
    }else{
      fmesh_par_cmesh[i]=cmesh_bound[i]-cmesh_bound[i-1];
    };
  };

  if(cmesh_bound[cmesh-1]!=TotM-1){
    cout<<"Error in MECSystem::PutCMeshBound.\n";
    cout<<"Final Mesh ID for CMesh is inconsistent.\n";
    exit(0);
  };

  int id1=0;
  for(int i=0;i<cmesh;i++){
    int id2=cmesh_bound[i];
    for(int j=id1;j<=id2;j++){
      cmeshid[j]=i;
    };
    id1=id2+1;
  };
};

void MECSystem::PutCellSymmetricData(int inp,int *i45,int *i90,int *i180)
{
  int reg_cell=inp;

  vector< vector< vector<int> > > cm_fm(cmeshy);
  int index=0;
  for(int i=0;i<cmeshy;i++){
    cm_fm[i].resize(cmeshx);
    for(int j=0;j<cmeshx;j++){
      cm_fm[i][j].resize(reg_cell,0);
      int tmp=reg_cell;
      if(fmesh_par_cmesh[i*cmeshx+j]==1)tmp=1;
      for(int k=0;k<tmp;k++){
	cm_fm[i][j][k]=index++;
      };
    };
  };

  regid_45.resize(TotM);
  regid_90.resize(TotM);
  regid_180.resize(TotM);
  for(int y=0;y<cmeshy;y++){
    for(int x=0;x<cmeshx;x++){
      if(fmesh_par_cmesh[y*cmeshx+x]==1){
	// 45-degree symmetric
        int id0=cm_fm[y][x][0];
	int id1=cm_fm[x][y][0];
	regid_45[id0]=id1;
	// 90-degree symmetric
	int id2=cm_fm[y][cmeshx-1-x][0];
	regid_90[id0]=id2;
	// 180-degree symmetric
	int id3=cm_fm[cmeshy-1-y][x][0];
	regid_180[id0]=id3;
      }else{
        for(int j=0;j<reg_cell;j++){
	  // 45-degree symmetric
	  int id0=cm_fm[y][x][j];
	  int id1=cm_fm[x][y][i45[j]];
          regid_45[id0]=id1;
	  // 90-degree symmetric
	  int id2=cm_fm[y][cmeshx-1-x][i90[j]];
	  regid_90[id0]=id2;
	  //180-degree symmetric
	  int id3=cm_fm[cmeshy-1-y][x][i180[j]];
	  regid_180[id0]=id3;
	};
      };
    };
  };
};

void MECSystem::PutCellSymmetricDataForAllOctant(int sqr)
{
  if(sqr<1||sqr>7||sqr==4||sqr==6){
    cout<<"Not coded::MECSystem::PutCellSymmetricDataForAllOctant.\n";
    exit(0);
  };

  vector< vector< vector<int> > > cm_fm(cmeshy);
  int index=0;
  for(int i=0;i<cmeshy;i++){
    cm_fm[i].resize(cmeshx);
    for(int j=0;j<cmeshx;j++){
      int tmp=fmesh_par_cmesh[i*cmeshx+j];
      cm_fm[i][j].resize(tmp,0);
      for(int k=0;k<tmp;k++){
	cm_fm[i][j][k]=index++;
      };
    };
  };

  int sym45[]={6,7,3,2,5,4,0,1};
  int sym90[]={3,2,1,0,7,6,5,4};
  int sym180[]={4,5,6,7,0,1,2,3};
  //   1 2 
  // 0     3
  // 4     7
  //   5 6 

  int sym45_sq2[]={3,1,2,0};
  int sym90_sq2[]={1,0,3,2};
  int sym180_sq2[]={2,3,0,1};
  // 0 1
  // 2 3

  int sym45_sq3[]={8,5,2,7,4,1,6,3,0};
  int sym90_sq3[]={2,1,0,5,4,3,8,7,6};
  int sym180_sq3[]={6,7,8,3,4,5,0,1,2};
  // 0 1 2
  // 3 4 5
  // 6 7 8

  int sym45_sq5[]={24,19,14,9,4,
                   23,18,13,8,3,
		   22,17,12,7,2,
		   21,16,11,6,1,
		   20,15,10,5,0};
  int sym90_sq5[]={4,3,2,1,0,
		   9,8,7,6,5,
		   14,13,12,11,10,
		   19,18,17,16,15,
		   24,23,22,21,20};
  int sym180_sq5[]={20,21,22,23,24,
		    15,16,17,18,19,
		    10,11,12,13,14,
		    5,6,7,8,9,
		    0,1,2,3,4};
  // 0 1 2 3 4
  // 5 6 7 8 9
  //1011121314
  //1516171819
  //2021222324

  int sym45_sq7[]={
    48,41,34,27,20,13, 6,
    47,40,33,26,19,12, 5,
    46,39,32,25,18,11, 4,
    45,38,31,24,17,10, 3,
    44,37,30,23,16, 9, 2,
    43,36,29,22,15, 8, 1,
    42,35,28,21,14, 7, 0
  };
  int sym90_sq7[]={
    6,5,4,3,2,1,0,
    13,12,11,10,9,8,7,
    20,19,18,17,16,15,14,
    27,26,25,24,23,22,21,
    34,33,32,31,30,29,28,
    41,40,39,38,37,36,35,
    48,47,46,45,44,43,42
  };
  int sym180_sq7[]={
    42,43,44,45,46,47,48,
    35,36,37,38,39,40,41,
    28,29,30,31,32,33,34,
    21,22,23,24,25,26,27,
    14,15,16,17,18,19,20,
    7,8,9,10,11,12,13,
    0,1,2,3,4,5,6};
  // 0 1 2 3 4 5 6
  // 7 8 910111213
  //14151617181920
  //21222324252627
  //28293031323334
  //35363738394041
  //42434445464748

  regid_45.resize(TotM);
  regid_90.resize(TotM);
  regid_180.resize(TotM);
  for(int y=0;y<cmeshy;y++){
    for(int x=0;x<cmeshx;x++){
      int tmp=fmesh_par_cmesh[y*cmeshx+x];
      int rrr=tmp%8;
      if(rrr!=0){
	if(sqr==1){
  	  // 45-degree symmetric
          int id0=cm_fm[y][x][0];
	  int id1=cm_fm[x][y][0];
	  regid_45[id0]=id1;
	  // 90-degree symmetric
	  int id2=cm_fm[y][cmeshx-1-x][0];
	  regid_90[id0]=id2;
	  // 180-degree symmetric
	  int id3=cm_fm[cmeshy-1-y][x][0];
	  regid_180[id0]=id3;
	}else if(sqr==2){
	  for(int j=0;j<4;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][j];
	    int id1=cm_fm[x][y][sym45_sq2[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][sym90_sq2[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][sym180_sq2[j]];
	    regid_180[id0]=id3;
	  };
	}else if(sqr==3){
	  for(int j=0;j<9;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][j];
	    int id1=cm_fm[x][y][sym45_sq3[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][sym90_sq3[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][sym180_sq3[j]];
	    regid_180[id0]=id3;
	  };
	}else if(sqr==5){
	  for(int j=0;j<5*5;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][j];
	    int id1=cm_fm[x][y][sym45_sq5[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][sym90_sq5[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][sym180_sq5[j]];
	    regid_180[id0]=id3;
	  };
	}else{
	  for(int j=0;j<7*7;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][j];
	    int id1=cm_fm[x][y][sym45_sq7[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][sym90_sq7[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][sym180_sq7[j]];
	    regid_180[id0]=id3;
	  };
	};
      }else{
        int maxx=tmp/8;
        for(int k=0;k<maxx;k++){
          for(int j=0;j<8;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][k*8+j];
	    int id1=cm_fm[x][y][k*8+sym45[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][k*8+sym90[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][k*8+sym180[j]];
	    regid_180[id0]=id3;
	  };
	};
      };
    };
  };
};

void MECSystem::PutCellSymmetricDataForAllQuarter()
{
  vector< vector< vector<int> > > cm_fm(cmeshy);
  int index=0;
  for(int i=0;i<cmeshy;i++){
    cm_fm[i].resize(cmeshx);
    for(int j=0;j<cmeshx;j++){
      int tmp=fmesh_par_cmesh[i*cmeshx+j];
      cm_fm[i][j].resize(tmp,0);
      for(int k=0;k<tmp;k++){
	cm_fm[i][j][k]=index++;
      };
    };
  };

  int sym45[]={3,1,2,0};
  int sym90[]={1,0,3,2};
  int sym180[]={2,3,0,1};
  // 0 1
  // 2 3


  regid_45.resize(TotM);
  regid_90.resize(TotM);
  regid_180.resize(TotM);
  for(int y=0;y<cmeshy;y++){
    for(int x=0;x<cmeshx;x++){
      int tmp=fmesh_par_cmesh[y*cmeshx+x];
      int rrr=tmp%4;
      if(rrr==1){
	// 45-degree symmetric
        int id0=cm_fm[y][x][0];
	int id1=cm_fm[x][y][0];
	regid_45[id0]=id1;
	// 90-degree symmetric
	int id2=cm_fm[y][cmeshx-1-x][0];
	regid_90[id0]=id2;
	// 180-degree symmetric
	int id3=cm_fm[cmeshy-1-y][x][0];
	regid_180[id0]=id3;
      }else{
        int maxx=tmp/4;
        for(int k=0;k<maxx;k++){
          for(int j=0;j<4;j++){
	    // 45-degree symmetric
	    int id0=cm_fm[y][x][k*4+j];
	    int id1=cm_fm[x][y][k*4+sym45[j]];
            regid_45[id0]=id1;
	    // 90-degree symmetric
	    int id2=cm_fm[y][cmeshx-1-x][k*4+sym90[j]];
	    regid_90[id0]=id2;
	    //180-degree symmetric
	    int id3=cm_fm[cmeshy-1-y][x][k*4+sym180[j]];
	    regid_180[id0]=id3;
	  };
	};
      };
    };
  };
};

void MECSystem::PutCellSymmetricDataForRectangular(int m)
{
  // 0 1 2 
  // 3 4 5
  // 6 7 8
  vector< vector< vector<int> > > cm_fm(cmeshy);
  int index=0;
  for(int i=0;i<cmeshy;i++){
    cm_fm[i].resize(cmeshx);
    for(int j=0;j<cmeshx;j++){
      int tmp=fmesh_par_cmesh[i*cmeshx+j];
      cm_fm[i][j].resize(tmp,0);
      for(int k=0;k<tmp;k++){
	cm_fm[i][j][k]=index++;
      };
    };
  };

  vector<int> sym45(m*m);
  vector<int> sym90(m*m);
  vector<int> sym180(m*m);
  int ii=0;
  for(int y=0;y<m;y++){
    for(int x=0;x<m;x++){
      sym45[ii]=(m-1-x)*m+(m-1-y);
      sym90[ii]=y*m+(m-1-x);
      sym180[ii]=(m-1-y)*m+x;
      ii++;
    };
  };

  regid_45.resize(TotM);
  regid_90.resize(TotM);
  regid_180.resize(TotM);
  for(int y=0;y<cmeshy;y++){
    for(int x=0;x<cmeshx;x++){
      int tmp=fmesh_par_cmesh[y*cmeshx+x];
      if(tmp==1){
	// 45-degree symmetric
        int id0=cm_fm[y][x][0];
	int id1=cm_fm[x][y][0];
	regid_45[id0]=id1;
	// 90-degree symmetric
	int id2=cm_fm[y][cmeshx-1-x][0];
	regid_90[id0]=id2;
	// 180-degree symmetric
	int id3=cm_fm[cmeshy-1-y][x][0];
	regid_180[id0]=id3;
      }else{
        for(int k=0;k<tmp;k++){
          // 45-degree symmetric
	  int id0=cm_fm[y][x][k];
	  int id1=cm_fm[x][y][sym45[k]];
          regid_45[id0]=id1;
          // 90-degree symmetric
	  int id2=cm_fm[y][cmeshx-1-x][sym90[k]];
	  regid_90[id0]=id2;
          //180-degree symmetric
	  int id3=cm_fm[cmeshy-1-y][x][sym180[k]];
	  regid_180[id0]=id3;
	};
      };
    };
  };
};

void MECSystem::PutCellSymmetricDataNoSymmetry()
{
  vector< vector< vector<int> > > cm_fm(cmeshy);
  int index=0;
  for(int i=0;i<cmeshy;i++){
    cm_fm[i].resize(cmeshx);
    for(int j=0;j<cmeshx;j++){
      int tmp=fmesh_par_cmesh[i*cmeshx+j];
      cm_fm[i][j].resize(tmp,0);
      for(int k=0;k<tmp;k++){
	cm_fm[i][j][k]=index++;
      };
    };
  };

  regid_45.resize(TotM);
  regid_90.resize(TotM);
  regid_180.resize(TotM);
  for(int y=0;y<cmeshy;y++){
    for(int x=0;x<cmeshx;x++){
      int tmp=fmesh_par_cmesh[y*cmeshx+x];
      for(int k=0;k<tmp;k++){
	int id0=cm_fm[y][x][k];
 	// 45-degree symmetric
  	int id1=cm_fm[x][y][k];
	regid_45[id0]=id1;
	// 90-degree symmetric
	int id2=cm_fm[y][cmeshx-1-x][k];
	regid_90[id0]=id2;
	// 180-degree symmetric
	int id3=cm_fm[cmeshy-1-y][x][k];
	regid_180[id0]=id3;
      };
    };
  };
};


// +++ For CMFD acceleration

void MECSystem::DoAcceleration(int iter,real errs,real fiss)
{
  if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+0.5);
  //if(iter%opt.GetItcmfd()==0)DoCMFDAcceleration(fiss+1e9);
};


void MECSystem::DoCMFDAcceleration(real delk)
{
  int ngrp=1;
  vector<int> bgrp(ngrp);
  bgrp[0]=grp-1;

  int xr=cmeshx;
  int yr=cmeshy;
  int zr=1;

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
  for(int i=0;i<CTotM;i++){
    cc_nusigf[i]=0.;
    cc_sigt[i]=0.;
    cc_d[i]=0.;
    cc_flx[i]=0.;
    cc_sigs[i]=0.;
  };

  GroupData1D volflx(grp);
  for(int i=0;i<cmesh;i++){
    int stt=0;
    if(i!=0)stt=cmesh_bound[i-1]+1;
    for(int j=stt;j<=cmesh_bound[i];j++){
      real vol=mesh[j].GetVolume();
      volflx=mesh[j].GetFlux()*vol;
      cc_d[i]+=volflx*mesh[j].GetMed()->GetMacxs().GetD();
      cc_sigt[i]+=volflx*mesh[j].GetMed()->GetMacxs().GetSigt();
      cc_sigs[i]+=volflx*mesh[j].GetMed()->GetMacxs().GetSigs().get_sumx();
      cc_nusigf[i]+=volflx*mesh[j].GetMed()->GetMacxs().GetNusigf();
      cc_flx[i]+=mesh[j].GetVolumeFlux();
    };
  };
  Medium minp(ngrp);
  minp.PutPL(0);
  real inv_delk=1./delk;
  for(int i=0;i<CTotM;i++){
    real tmp=1./cc_flx[i];
    real fis=cc_nusigf[i]*tmp;
    minp.GetMacxs().GetNusigf().put_data(0,fis);
    minp.GetMacxs().GetKai().put_data(0,1.);
    real sigtfic=cc_sigt[i]*tmp-fis*inv_delk;
    minp.GetMacxs().GetSigt().put_data(0,sigtfic);
    minp.GetMacxs().GetSigs().put_data(0,0,cc_sigs[i]*tmp);
    real dinp=cc_d[i]*tmp*cmfd_factor;
    if(dinp>50.)dinp=50.;
    // +++ For voided region +++
    minp.GetMacxs().GetD().put_data(0,dinp);
    cm.AddMedium(minp);
  };

  vector<real> xwid(xr);
  vector<real> ywid(yr);
  vector<int> fmx(xr);
  vector<int> fmy(yr);
  for(int i=0;i<xr;i++){
    xwid[i]=cmxl[i];
    fmx[i]=1;
  };
  for(int i=0;i<yr;i++){
    ywid[i]=cmyl[i];
    fmy[i]=1;
  };
  vector<int> asmmap(xyzr);
  int index=0;
  for(int z1=0;z1<zr;z1++){
    for(int y1=0;y1<yr;y1++){
      for(int x1=0;x1<xr;x1++){
	asmmap[index]=index;
	index++;
      };
    };
  };
  CartMeshInfo cmi;
  cmi.PutMeshInfo(xr,yr,fmx,fmy,xwid,ywid,asmmap);
  int bci[6];
  int org=1; // reflective
  if(Vacuum)org=2; // Vacuum
  for(int i=0;i<6;i++){
    bci[i]=org;
  };
  cmi.PutBoundaryCondition(bci);

  string ss="Cartesian";
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
    cc_flx[i]=cm.GetMesh(i).GetFlux().get_dat(0)*tmp/cc_flx[i];
  };

  for(int i=0;i<cmesh;i++){
    int stt=0;
    if(i!=0)stt=cmesh_bound[i-1]+1;
    real tmp=cc_flx[i];
    for(int j=stt;j<=cmesh_bound[i];j++){
      for(int g=0;g<grp;g++){
        mesh[j].GetFlux().put_data(g,mesh[j].GetFlux().get_dat(g)*tmp);
      };
    };
  };
};

void MECSystem::SetZeroCurFF()
{
  for(int i=0;i<1;i++){
    for(int j=0;j<cmeshy+1;j++){
      for(int k=0;k<cmeshx+1;k++){
        for(int l=0;l<Ndim;l++){
    	  CurFF[0][i][j][k][l]=0.;
	};
      };
    };
  };
};

GroupData1D MECSystem::GetIntegratedFlux(int medid,int mom)
{
  real *sum=new real[grp];
  for(int i=0;i<grp;i++){sum[i]=0.;};
  for(int i=0;i<TotM;i++){
    if(medid_fmesh[i]==medid){
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

GroupData1D MECSystem::GetIntegratedReactionRate(enum xstype react, int medid)
{
  GroupData1D ret(grp);
  ret.set_zero();
  for(int i=0;i<nmed;i++){
    if(medid==-1||medid==i){
      ret=ret+med[i].GetMacxs().GetData1d(react).mult(GetIntegratedFlux(i));
    };
  };
  return ret;
};

real MECSystem::GetVolumeParMedium(int medid)
{
  real ret=0.;
  for(int i=0;i<TotM;i++){
    if(medid_fmesh[i]==medid)ret+=mesh[i].GetVolume();
  };
  return ret;
};

real MECSystem::CalReactivity(MECSystem *sec,real kunp,real kp,bool pr,bool ipcal)
{
  vector<real> tmp;
  tmp.resize(grp,0.);
  return CalReactivity(sec,kunp,kp,tmp,pr,ipcal);
};


real MECSystem::CalReactivity(MECSystem *sec,real kunp,real kp,vector<real> &rho,bool pr,bool ipcal)
{
  CheckSameMesh(sec);
  if(pr)WritePerturbName();

  if(GetGeneralOption().Forward()){
    cout<<"# Error in MECSystem::CalReactivity.\n";
    cout<<"# Adjoint flux should be calculated in reference system.\n";
    exit(0);
  };
  if(!sec->GetGeneralOption().Forward()){
    cout<<"# Error in MECSystem::CalReactivity.\n";
    cout<<"# Forward flux should be calculated in perturbated system.\n";
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
    sct[0][i]=CalPerturbScatteringTerm(sec,i,flag);
    leak[0][i]=0.;
    //leak[0][i]=CalPerturbP1LeakageTerm(sec,i,flag);
    //CalPerturbLeakScat(sec,flag,i,scttmp,leaktmp);
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
    cout<<"# (Unit of Reactivity : (Reactivity)/(Unit_lethargy:0.25)\n";
    cout<<"#     Energy      ";
    cout<<"Yield       ";
    cout<<"Absorption  ";
    cout<<"Scat+Leak   ";
    //cout<<"Leakage     ";
    cout<<"Total\n";
  };

  real inv_ip=1./ip;
  for(int i=0;i<grp;i++){
    yld[i]*=inv_ip;
    abs[i]*=inv_ip;
    real sl=0.;
    real le=0.;
    for(int j=0;j<plnum;j++){
      sct[j][i]*=inv_ip; sctsum[j]+=sct[j][i]; sl+=sct[j][i];
      leak[j][i]*=inv_ip; leaksum[j]+=leak[j][i]; le+=leak[j][i];
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
      //cout.precision(10);
      cout.precision(4);
      real en=mesh[0].GetMed()->GetEnband().get_dat(i);
      real en_next=mesh[0].GetMed()->GetEnband().get_dat(i+1);
      real leth=(log(en/en_next)/0.25);
      cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth;
      //cout<<en<<"  "<<yld[i]/leth<<"  "<<abs[i]/leth<<"  "<<sl/leth<<"  "<<le/leth;
      //cout<<en<<"  "<<yld[i]<<"  "<<abs[i]<<"  "<<sl<<"  "<<le;
      cout<<" "<<grpsum/leth<<"\n";
      cout.unsetf(ios::scientific);
    };
  };

  if(pr){
    cout<<"#\n";
    cout<<"# Yield       : "<<yldsum<<"\n";
    cout<<"# Absorption  : "<<abssum<<"\n";
    /*
    int id=0;
    for(int l=0;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th scattering  :"<<sctsum[id]<<"\n";
        id++;
      };
    };
    */
    cout<<"# Scat+Leak   : "<<tsctsum<<"\n";
    /*
    id=0;
    for(int l=1;l<=pl;l++){
      for(int m=0;m<=0;m++){ // snr
        cout<<"#  ("<<l<<","<<m<<")th leakage     :"<<leaksum[id]<<"\n";
        id++;
      };
    };
    */
    //cout<<"#     Higher leakage :"<<leaksum[plnum-1]<<"\n";
    //cout<<"# Leakage      :"<<tleaksum<<"\n";
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

real MECSystem::CalPerturbScatteringTerm(MECSystem *sec,int g,bool *flag)
{
  // Weights of angular quadrature set is normalized to 4PI.

  int sn=GetQuad().GetSN();
  int pdiv=polar_angle_division;
  int sn_per_pang=sn/pdiv;

  /*
  real sum=0.;
  for(int i=0;i<sn;i++){
    sum+=GetQuad().GetOmega(i);
  };
  cout<<sum<<"\n"; exit(0);
  */

  real ret=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){

      /*
      real tmp1=0.;
      real tmp2=0.;
      for(int j=0;j<sn;j++){
	real tmp=sec->GetAFlux(m,g).get_dat(j)*GetQuad().GetOmega(j);
	tmp1+=tmp*GetQuad().GetMoment(1,j);
	tmp2+=tmp*GetQuad().GetMoment(2,j);
      };
      cout<<tmp1<<" "<<sec->GetMesh(m).GetFlux(1).get_dat(g)<<" ";
      cout<<tmp2<<" "<<sec->GetMesh(m).GetFlux(2).get_dat(g)<<" ";
      cout<<"\n";
      */

      real vol=GetMesh(m).GetVolume();
      real f2v=sec->GetMesh(m).GetFlux(0).get_dat(g)*vol; // (forward-flux)*(volume)
      // (scattering term)
      real tmp=0.;
      for(int j=0;j<grp;j++){ // w up scattering
        tmp+=(sec->GetMesh(m).GetMed()->GetMacxs().GetSigs(0).get_dat(g,j)
	             -mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,j))*
	  mesh[m].GetFlux(0).get_dat(j);
      };
      // (total term)
      real sum=0.;
      real dsig=sec->GetMesh(m).GetMed()->GetMacxs().GetData1d(sigt).get_dat(g)
                    -GetMesh(m).GetMed()->GetMacxs().GetData1d(sigt).get_dat(g);
      // ... isotropic approximation
      //sum+=dsig*vol*sec->GetMesh(m).GetFlux(0).get_dat(g)*mesh[m].GetFlux(0).get_dat(g);     
      // ... based on angular flux

      for(int i=0;i<sn;i++){
        int ii=i;
        int mm=m;
	if(degree360){
          int tmp=i%sn_per_pang;
          if(tmp<sn_per_pang/2){
	    ii=i+sn_per_pang/2;
	  }else{
	    ii=i-sn_per_pang/2;
	  };
	}else{
	  // The following are commented out in 2022/2/26
	  // during the resonance scattering treatment study with Mosteller's benchmark
	  // to use the MOC-perturbation capability.
	  // If the following works, the results are very different from the expected one.
	  /*
          int hm=TotM/2;
	  if(m>=hm){
	    mm=m-hm;
	  }else{
	    mm=m+hm;
	  };
	  */
	};
        real omega=GetQuad().GetOmega(i);
        sum+=dsig*omega*sec->GetAFlux(m,g).get_dat(i)*GetAFlux(mm,g).get_dat(ii)*vol*PI4;
	// PI4 is multiplied because of consistency
      };

      // ... based on angular moment
      /*
      sum+=dsig*vol*sec->GetMesh(m).GetFlux(0).get_dat(g)*mesh[m].GetFlux(0).get_dat(g);     
      sum-=3.*dsig*vol*sec->GetMesh(m).GetFlux(1).get_dat(g)*mesh[m].GetFlux(1).get_dat(g);     
      sum-=3.*dsig*vol*sec->GetMesh(m).GetFlux(2).get_dat(g)*mesh[m].GetFlux(2).get_dat(g);     
      */
      // (absorption term)
      real dsiga=sec->GetMesh(m).GetMed()->GetMacxs().GetData1d(siga).get_dat(g)
                     -GetMesh(m).GetMed()->GetMacxs().GetData1d(siga).get_dat(g);
      tmp+=dsiga*mesh[m].GetFlux(0).get_dat(g);
      ret+=tmp*f2v-sum;
    };
  };
  return ret;

};

real MECSystem::CalPerturbP1LeakageTerm(MECSystem *sec,int g,bool *flag)
{
  real tmp=0.;
  for(int m=0;m<TotM;m++){
    if(flag[m]){
      real vol=mesh[m].GetVolume();
      real dtot=sec->GetMesh(m).GetMed()->GetMacxs().GetSigt(0).get_dat(g)
	               -mesh[m].GetMed()->GetMacxs().GetSigt(0).get_dat(g);
      //tmp-=3.*dtot*sec->GetMesh(m).GetFlux(1).get_dat(g)*mesh[m].GetFlux(1).get_dat(g)*vol;
      //tmp-=3.*dtot*sec->GetMesh(m).GetFlux(2).get_dat(g)*mesh[m].GetFlux(2).get_dat(g)*vol;
      tmp+=3.*dtot*sec->GetMesh(m).GetFlux(1).get_dat(g)*mesh[m].GetFlux(1).get_dat(g)*vol;
      tmp+=3.*dtot*sec->GetMesh(m).GetFlux(2).get_dat(g)*mesh[m].GetFlux(2).get_dat(g)*vol;
      // +++ signs are opposite because currents in adjoint calculation should be inversed. 
    };
  };
  return tmp;
};

SensitivityData MECSystem::CalSensitivityNew(MECSystem *sec,real keff,int nucnum,int *nucid,bool ipcal)
{
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
      cout<<"# Error in MECsystem::CalSensitivityNew.\n";
      cout<<"# Perturbation denominator is zero.\n";
      exit(0);
    };
  };

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    cout<<"# Sensitivity calculation for nuclide : "<<nid<<" ("<<nc<<"/"<<nucnum<<")\n";

    bool nuclide_in_system=false;
    for(int i=0;i<TotM;i++){
      if(sec->GetMesh(i).GetMed()->ExistNuclide(nid)){
        flag[i]=true;
	nuclide_in_system=true;
      }else{
	flag[i]=false;
      };
    };

    if(nuclide_in_system){

      bool fissile=false;
      for(int i=0;i<nmed;i++){
        if(sec->GetMedium(i).ExistNuclide(nid)){
          for(int j=0;j<grp;j++){
	    //if(med[i].GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(j)>0.0){
	    if(sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(j)>0.0){
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
  	    if(sec->GetMedium(j).ExistNuclide(nid)){
              real den=sec->GetMedium(j).GetNuclide(nid).GetDensity();
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
  	    if(sec->GetMedium(j).ExistNuclide(nid)){
  	      sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	      sec->GetMed(j).GetMacxs().GetSiga().put_data(i,absorg[j]);
	    };
	  };
        };
        sens.PutSensitivity1D(nid,18,sns1d);

        // +++ Nu
        for(int i=0;i<grp;i++){
          for(int j=0;j<nmed;j++){
    	    if(sec->GetMedium(j).ExistNuclide(nid)){
              real den=sec->GetMedium(j).GetNuclide(nid).GetDensity();
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
            if(sec->GetMedium(j).ExistNuclide(nid)){
  	    sec->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
	    };
          };
        };
        sens.PutSensitivity1D(nid,452,sns1d);

        // +++ Chi
        for(int i=0;i<nmed;i++){
	  if(sec->GetMedium(i).ExistNuclide(nid)){
  	    real den=sec->GetMedium(i).GetNuclide(nid).GetDensity();
 	    real total_fiss=sec->GetIntegratedFlux(i)*sec->GetMed(i).GetMacxs().GetNusigf();
            real part_fiss=den*(sec->GetIntegratedFlux(i)*sec->GetMed(i).GetNuclide(nid).GetMicxs().GetMicNusigf());
	    fiss_frac[i]=part_fiss/total_fiss;
	  }else{
	    fiss_frac[i]=0.;
	  };
        };
        for(int i=0;i<grp;i++){
	  for(int j=0;j<nmed;j++){
	    if(sec->GetMedium(j).ExistNuclide(nid)){
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
	    if(sec->GetMedium(j).ExistNuclide(nid)){
	      sec->GetMed(j).GetMacxs().GetKai().put_data(i,totorg[j]);
	    };
	  };
        };
        sens.PutSensitivity1D(nid,181,sns1d);
      }; //(end for fissile nuclide)


      // +++ Capture
      for(int i=0;i<grp;i++){
        for(int j=0;j<nmed;j++){
    	  if(sec->GetMedium(j).ExistNuclide(nid)){
            real den=sec->GetMedium(j).GetNuclide(nid).GetDensity();
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
  	  if(sec->GetMedium(j).ExistNuclide(nid)){
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
          for(int k=i;k<grp;k++){
  	    if((ii!=0)||(nid<100000)||(ii==0&&k<=i+2)){
              for(int j=0;j<nmed;j++){
                if(sec->GetMedium(j).ExistNuclide(nid)){
                  real den=sec->GetMedium(j).GetNuclide(nid).GetDensity();
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
  	        if(sec->GetMedium(j).ExistNuclide(nid)){
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


SensitivityData MECSystem::CalSensitivityRRR(MECSystem &lat,int nume_nuc,int *nume_id,enum xstype *nume_xs,int denom_nuc,int *denom_id,enum xstype *denom_xs,bool *on_mesh,int nucnum,int *nucid,real keff)
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

  //SensitivityData sns;
  return sns;
};

void MECSystem::CalGPT_MEC(real keff,real eps_f,int itermax)
{
  //eps_f=1e-20;
  //itermax=20;
  
  CheckUpScattering();

  /*
  // .....
  vector<GroupData1D> extsrc(TotM);
  for(int m=0;m<TotM;m++){
    extsrc[m].put_imax(grp);
    for(int g=0;g<grp;g++){
      extsrc[m].put_data(g,mesh[m].GetScatSrc(g,0));
    };
  };
  // .....
  */
  
  vector<GroupData1D> gpt_flx_store(TotM);
  vector<GroupData1D> gpt_flx_old(TotM);
  for(int i=0;i<TotM;i++){
    gpt_flx_store[i].put_imax(grp);
    gpt_flx_old[i].put_imax(grp);
  };

  int totsn=GetQuad().GetSN();
  vector< vector<GroupData1D> > gpt_aflx_store;
  if(wrtflx){
    gpt_aflx_store.resize(TotM);
    for(int i=0;i<TotM;i++){
      gpt_aflx_store[i].resize(grp);
      for(int g=0;g<grp;g++){
	gpt_aflx_store[i][g].put_imax(totsn);
      };
    };
  };

  int iter=0;
  //real sumold=0.;
  real sum=0.;
  for(int k=0;k<itermax;k++){

    //sumold=sum;
    sum=0.;
    for(int m=0;m<TotM;m++){
      if(k!=0)mesh[m].CalFissionSrcAdjoint();
      sum+=mesh[m].GetFissionSrc();
    };
    //cout.setf(ios::showpoint);
    //cout.precision(7);
    //cout<<k<<"  :  "<<sum<<"\n";    
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

    int grp_errmax=-1;
    int mesh_errmax=-1;
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
          if(err>errmax){
	    errmax=err;
	    grp_errmax=g;
	    mesh_errmax=m;
	  };
        };
	if(fabs(tmp)>flxmax)flxmax=fabs(tmp);
        gpt_flx_old[m].put_data(g,tmp);
      };
    };
      

    //cout<<"# CALGPT : "<<k<<" "<<flxmax<<"\n";
    cout<<"# CALGPT : "<<k<<" "<<errmax<<" in mesh "<<mesh_errmax<<", group "<<grp_errmax<<" (Maximum flux : "<<flxmax<<")\n";

    if(k==0){
      gpt_flx_store=gpt_flx_old;
      if(wrtflx)gpt_aflx_store=aflux;
      /*
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m].copy(gpt_flx_old[m]);
	if(wrtflx){
	  gpt_aflx_store=aflux;
	};
      };
      */
    }else{
      for(int m=0;m<TotM;m++){
	gpt_flx_store[m]=gpt_flx_store[m]+gpt_flx_old[m];
	if(wrtflx){
	  for(int g=0;g<grp;g++){
	    gpt_aflx_store[m][g]=gpt_aflx_store[m][g]+aflux[m][g];
	  };
	};
      };
    };

    if(errmax<eps_f&&k!=0){

      /*
      cout<<"#\n# Source distribution with volume\n#\n";
      for(int m=0;m<TotM;m++){
	for(int g=0;g<grp;g++){
  	  cout<<m<<"   "<<g<<"   "<<extsrc[m].get_dat(g)<<"\n";
	};
      };
      cout<<"\n\n\n";
      */
       
      break;
    };

    if(k==itermax-1){
      cout<<"# Warning in GeneralSystem::CalGPT.\n";
      cout<<"# Iteration reaches to maximum ("<<itermax<<").\n";
      cout<<"# Iteration residual : "<<errmax<<"\n";

      /*
      cout<<"#\n# Source distribution with volume\n#\n";
      for(int m=0;m<TotM;m++){
	for(int g=0;g<grp;g++){
  	  cout<<m<<"   "<<g<<"   "<<extsrc[m].get_dat(g)<<"\n";
	};
      };
      cout<<"\n\n\n";

      cout<<"#\n# Flux and total/self-scattering xs distribution\n#\n";
      for(int m=0;m<TotM;m++){
	cout<<"#\n# mesh : "<<m<<"\n";
	for(int g=0;g<grp;g++){
	  cout<<g<<" "<<gpt_flx_store[m].get_dat(g)<<" "<<mesh[m].GetMed()->GetMacxs().GetSigt().get_dat(g)<<" "<<mesh[m].GetMed()->GetMacxs().GetSigs(0).get_dat(g,g)<<"\n";
	};
      };
      */
      exit(0);
      
    };
  };

  for(int m=0;m<TotM;m++){
    mesh[m].GetFlux().copy(gpt_flx_store[m]-gpt_flx_old[m]*iter);
    if(wrtflx){
      for(int g=0;g<grp;g++){
	aflux[m][g]=gpt_aflx_store[m][g]-aflux[m][g]*iter;
      };
    };
  };

  //mesh[0].GetFlux().show_self();

};

