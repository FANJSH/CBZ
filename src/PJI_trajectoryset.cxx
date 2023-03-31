#include <cstdlib>
#include "PJI_trajectoryset.h"

TrajectorySet::TrajectorySet()
{
  fkk.init(0.00000001);
  cal_vol_sph=false;
  angle360=false;
};

void TrajectorySet::PutAngularDivision(int i)
{
  angular_division=i;
  trajectory_on_angular_boundary.resize(angular_division,0);
  discreted_angle.resize(angular_division);
  vdir.resize(angular_division);
};

void TrajectorySet::CalculationPij(real *xs,real *pijinp,bool aniso)
{
  switch (BoundaryCondition) {
  case White:
    Cal_White(xs,pijinp,0,aniso);
    break;
  case Periodic:
    Cal_Periodic(xs,pijinp,aniso);
    break;
  case Reflective:
    Cal_Periodic(xs,pijinp,aniso);
    break;
  case Black:
    Cal_White(xs,pijinp,1,aniso);
    break;
  default:
    cout<<"Error in TrajectorySet::CalculationPij.\n";
    cout<<"Boundary condition is not yet fixed.\n";
    exit(0);
  };
};

void TrajectorySet::PutXS(real *xs)
{
  for(int i=0;i<NumTrajectory;i++){
    line[i].PutXS(xs);
  };
};

void TrajectorySet::CalVolume()
{
  int mk=regnum;
  real inv_pi=1./PI;

  //volume.resize(mk,0.);
  volume.resize(mk);
  for(int i=0;i<mk;i++){
    volume[i]=0.;
  };

  real *length=new real[mk];
  for(int i=0;i<NumTrajectory;i++){
    real width_frac=line[i].GetWeight()*inv_pi; 
    line[i].GetRegionwiseLength(mk,length);
    for(int j=0;j<mk;j++){
      volume[j]+=length[j]*width_frac;
    };
  };

  for(int i=0;i<mk;i++){
    if(volume[i]==0.){
      cout<<"# !! Warning !!\n";
      cout<<"# Zero volume is detected in TrajectorySet::CalVol. (region : "<<i<<")\n";
      cout<<"# Since zero volume is not allowed, so volume is reset to 1e-20.\n";
      cout<<"# This treatment is temporal, and shoud be revised in future.\n";
      volume[i]=1e-20;
    };
  };

  // ++ Adjustment of segment-length 
  // ++ to preserve region-wise volume for each azimuthal angle

  vector<real> voltmp(mk);
  int id1=0;
  int id2;
  for(int i=0;i<angular_division;i++){
    id2=trajectory_on_angular_boundary[i];
    for(int j=0;j<mk;j++){voltmp[j]=0.;};
    for(int j=id1;j<=id2;j++){
      real width=line[j].GetWeight()*(angular_division*inv_pi);
      line[j].GetRegionwiseLength(mk,length);
      for(int k=0;k<mk;k++){
	voltmp[k]+=length[k]*width;
      };
    };
    for(int j=0;j<mk;j++){
      voltmp[j]=volume[j]/voltmp[j]; // volume adjustment factor
    };
    for(int j=id1;j<=id2;j++){
      int rr=line[j].GetNumreg();
      for(int k=0;k<rr;k++){
	int irg=line[j].GetReg(k);
	line[j].PutDist(k,line[j].GetDist(k)*voltmp[irg]);
      };
    };
    id1=id2+1;
  };


  delete [] length;
}

void TrajectorySet::CalVolumeSphere(int r,int *rid,real *rinp)
{
  cal_vol_sph=true;
  for(int i=0;i<regnum;i++){
    volume[i]=0.;
  };
  real volt=0.;
  for(int i=0;i<r;i++){
    real tmp=rinp[r-i-1];
    real v=PI*4.*0.3333333333333*tmp*tmp*tmp;
    volume[rid[r-i-1]]+=v-volt;
    volt=v;
  };
  /*
  for(int i=0;i<regnum;i++){
    cout<<i<<" "<<volume[i]<<"\n";
  };
  */

};

void TrajectorySet::CalculationPijSphere(real *xs,real *pij)
{
  if(!cal_vol_sph){
    cout<<"\n# Error in TrajectorySet::CalculationPijSphere.\n";
    cout<<"# The method [CalVolumeSphere] is not yet done.\n";
    exit(0);
  };

  if(BoundaryCondition!=White&&BoundaryCondition!=Black){
    cout<<"\n# Error in TrajectorySet::CalculationPijSphere.\n";
    cout<<"# White or Black boundary conditions are only acceptable.\n";
    exit(0);
  };

  vector<real> pis(regnum);
  vector<real> psi(regnum);
  real pss=0.;

  PutXS(xs);

  real tsum=0.;
  for(int i=0;i<regnum;i++){
    tsum+=volume[i];
  };
  real r=tsum*INV_PI*3.*0.25;
  r=pow(r,0.333333333333);

  int nr=regnum;
  int nr1=nr+1;

  for(int i=0;i<nr*nr;i++){pij[i]=0.;};

  real *pijtmp=new real[nr1*nr1];
  int numt=NumTrajectory;
  for(int tr=0;tr<numt;tr++){
    int numr=line[tr].GetNumreg();
    real width=line[tr].GetWeight()*INV_PI;

    real tl=0.;
    for(int j=0;j<numr;j++){
      tl+=line[tr].GetDist(j);
    };
    tl*=0.5;
    real rho=0.;
    real tmp=r*r-tl*tl;
    if(tmp>0.)rho=sqrt(tmp);

    real weight=width*rho;
    line[tr].CalPijSphere(nr,pijtmp);
    //line[tr].CalPijSphereWithExpTable(nr,etab,pijtmp);

    for(int i=0;i<nr;i++){
      for(int j=0;j<nr;j++){
	pij[i*nr+j]+=pijtmp[i*nr1+j]*weight;
      };
      pis[i]+=pijtmp[i*nr1+nr]*weight;
      psi[i]+=pijtmp[nr*nr1+i]*weight;
    };
    pss+=pijtmp[nr*nr1+nr]*weight;
  };

  delete [] pijtmp;


  for(int i=0;i<nr;i++){
    real sigt=xs[i];
    real factor=2.*PI/volume[i]/sigt*0.5;
    for(int j=0;j<nr;j++){
      pij[i*nr+j]*=factor;
    };
    pis[i]*=factor;
  };

  real surface=4*PI*r*r;
  real factor=4.*PI/surface;
  for(int i=0;i<nr;i++){
    psi[i]*=factor;
  };
  pss*=factor;

  // (for white BC)
  if(BoundaryCondition==White){
  real pss_inv=1./(1.-pss);
  for(int i=0;i<nr;i++){
    real sum=0.;
    for(int j=0;j<nr;j++){
      pij[i*nr+j]+=pis[i]*psi[j]*pss_inv;
      sum+=pij[i*nr+j];
    };
    sum=1./sum;
    for(int j=0;j<nr;j++){
      pij[i*nr+j]*=sum;
    };
  };
  };
  
};

void TrajectorySet::Cal_White(real *xs,real *pij,bool black,bool aniso)
{
  int mk=regnum;
  int mkmk=mk*mk;
  int mk1mk1=(mk+1)*(mk+1);

  real *pijdummy=new real[mk1mk1];

  vector<real> pis(mk);
  vector<real> psi(mk);
  real pss=0.0;

  PutXS(xs);


  for(int i=0;i<mkmk;i++){pij[i]=0.;};
  for(int i=0;i<mk;i++){pis[i]=0.; psi[i]=0.;};

  for(int i=0;i<NumTrajectory;i++){

    if(!aniso){
      line[i].CalPij(mk,fkk,pijdummy);
    }else{
      line[i].CalPijAniso(mk,fkk,pijdummy);
    };


    real wgt=line[i].GetWeight();
    for(int j=0;j<mk1mk1;j++){
      pijdummy[j]*=wgt;
    };

    int tmp=0;
    for(int j=0;j<mk;j++){
      for(int k=0;k<mk;k++){
	pij[tmp]+=pijdummy[k*mk+j]+pijdummy[tmp];
	tmp++;
      };
    };

    //cout<<pijdummy[0]<<"\n";

    for(int j=0;j<mk;j++){
      pis[j]+=pijdummy[mkmk+j]+pijdummy[mk*(mk+1)+j];
    };

    pss+=pijdummy[mk1mk1-1]*2;

  };

  real fact=1.0;
  if(aniso)fact=1.125;

  for(int i=0;i<mk;i++){psi[i]=pis[i];};

  real tmp=4./(surface*PI2*fact);
  for(int i=0;i<mk;i++){
    real xsvol=1./(xs[i]*volume[i]*PI2);
    for(int j=0;j<mk;j++){
      pij[i*mk+j]*=xsvol;
    };
    pis[i]*=xsvol;
    psi[i]*=tmp;
  };
  pss*=tmp;

  /*
  cout.setf(ios::showpoint);
  cout.precision(8);
  cout<<"# pii : "<<pij[0]<<" "<<pij[1]<<" "<<pij[2]<<" "<<pij[3]<<"\n";
  cout<<"# pis : "<<pis[0]<<" "<<pis[1]<<"\n";
  cout<<"# psi : "<<psi[0]<<" "<<psi[1]<<"\n";
  cout<<"# pss : "<<pss<<"\n";
  */

  if(!black){
    tmp=1./(1.-pss);
    for(int i=0;i<mk;i++){
      for(int j=0;j<mk;j++){
        pij[i*mk+j]+=pis[i]*psi[j]*tmp;
      };
    };

    for(int i=0;i<mk;i++){
      real sum=0;
      for(int j=0;j<mk;j++){
        sum+=pij[i*mk+j];
      };
      sum=1./sum;
      for(int j=0;j<mk;j++){
        pij[i*mk+j]*=sum;
      };

    };
  };

  delete [] pijdummy;
}

void TrajectorySet::Cal_Periodic(real *xs,real *pij,bool aniso)
{
  real epsf= 3.0*2.3025; // exp(epsf)=10^{-3}
  epsf     = 6.; // Recommended value in the SRAC-2K6 manual (p.50)

  int max_trace_cell=100;

  bool opposite_direction=angle360; // angle=360

  int mk=regnum;
  int mkmk=mk*mk;
  int mk1mk1=(mk+1)*(mk+1);
  int numt=NumTrajectory;

  real *pijtemp=new real[mk1mk1];

  PutXS(xs);
  for(int i=0;i<mkmk;i++){pij[i]=0.;};

  for(int i=0;i<numt;i++){

    if(!aniso){
      line[i].CalPij(mk,fkk,pijtemp);
    }else{
      line[i].CalPijAniso(mk,fkk,pijtemp);
    };
    if(opposite_direction){
      // for opposite direction (angle=360)
      real *pijdummy=new real[mk1mk1];
      for(int j=0;j<mkmk;j++){pijdummy[j]=pijtemp[j];};
      for(int j=0;j<mk;j++){
        int tmp=j*mk;
        for(int k=0;k<mk;k++){
	  pijtemp[tmp++]+=pijdummy[k*mk+j];
        };
      };
      delete [] pijdummy;
    };

    if(BoundaryCondition==Periodic){
      // (Periodic)
      real Opttemp=0.;
      int NextT=NextTrajectory[i];
      int sumss=0;
      while(Opttemp<epsf&&sumss<max_trace_cell){
        if(!aniso){
          line[i].CalPijIso(mk,line[NextT],Opttemp,pijtemp,fkk,false,false);
        }
        else{
          line[i].CalPijAniso(mk,line[NextT],Opttemp,pijtemp,fkk,false,false);
        };
        Opttemp+=line[NextT].GetTotalOpt();
        NextT=NextTrajectory[NextT];
        sumss++;
      };

      // ** for opposite direction (angle=360)
      if(opposite_direction){
        Opttemp=0.;
        int PreT=PreTrajectory[i];
        int sumss=0;
        while(Opttemp<epsf&&sumss<max_trace_cell){
          if(!aniso){
            line[i].CalPijIso(mk,line[PreT],Opttemp,pijtemp,fkk,true,true);
          }
          else{
            line[i].CalPijAniso(mk,line[PreT],Opttemp,pijtemp,fkk,true,true);
          };
          Opttemp+=line[PreT].GetTotalOpt();
          PreT=PreTrajectory[PreT];
	  sumss++;
        };
      };

    }else{
      // (Reflective)
      real Opttemp=0.;
      bool reverse=false;
      int NextT=NextTrajectory[i];
      while(Opttemp<epsf){
        if(NextT>=numt){
          reverse=!reverse;
          NextT-=numt;
        };
        if(!aniso){
          line[i].CalPijIso(mk,line[NextT],Opttemp,pijtemp,fkk,false,reverse);
        }else{
          line[i].CalPijAniso(mk,line[NextT],Opttemp,pijtemp,fkk,false,reverse);
        };
        Opttemp+=line[NextT].GetTotalOpt();
        if(reverse){
	  NextT=PreTrajectory[NextT];
        }else{
          NextT=NextTrajectory[NextT];
        };
      };

      // ** for opposite direction (angle=360)
      Opttemp=0.;
      reverse=true;
      int PreT=PreTrajectory[i];
      while(Opttemp<epsf){
        if(PreT>=numt){
	  reverse=!reverse;
	  PreT-=numt;
        };
        if(!aniso){
          line[i].CalPijIso(mk,line[PreT],Opttemp,pijtemp,fkk,true,reverse);
        }else{
          line[i].CalPijAniso(mk,line[PreT],Opttemp,pijtemp,fkk,true,reverse);
        };
        Opttemp+=line[PreT].GetTotalOpt();
        if(reverse){
	  PreT=PreTrajectory[PreT];
        }else{
          PreT=NextTrajectory[PreT];
        };
      };
    };

    real wgt=line[i].GetWeight();
    for(int j=0;j<mkmk;j++){
      pij[j]+=pijtemp[j]*wgt;
    };

  };

  int ii=0;
  for(int i=0;i<mk;i++){
    real xsvol=xs[i]*volume[i]*PI2;
    xsvol=1./xsvol;
    for(int j=0;j<mk;j++){
      pij[ii++]*=xsvol;
    };
  };

  for(int i=0;i<mk;i++){
    real sum=0;
    int imk=i*mk;
    for(int j=0;j<mk;j++){
      sum+=pij[imk++];
    };
    sum=1./sum;
    imk-=mk;
    for(int j=0;j<mk;j++){
      pij[imk++]*=sum;
    };
  };

  delete [] pijtemp;
}

void TrajectorySet::CalTrajectory(IrregularGeometryInformation &ginp,int AngleDiv,real dlinp,real angle)
{
  line.clear();
  NextTrajectory.clear();
  PreTrajectory.clear();
  volume.clear();

  if(angle==360.){
    angle360=true;
    angle=180.;
    AngleDiv/=2;
  };

  if(BoundaryCondition==Reflective){
    if(!angle360){
      cout<<"#\n# Error in TrajectorySet::CalTrajectory.\n";
      cout<<"# When reflective BC is assigned,\n";
      cout<<"# angle value should be 360.\n";
      exit(0);
    };
  };

  PutAngularDivision(AngleDiv);

  cout<<"# Performing trajectory calculation\n";
  cout<<"#     (Theta / dl / TJs / accumulated-TJs)\n";

  regnum=ginp.GetRegion();

  real edge=ginp.GetEdge();
  surface=ginp.GetSurface();

  NumTrajectory=0;

  real dl;

  vector<GeomVector> pst;
  vector<GeomVector> ped;

  GeomVector vec1,vec2;
  Trajectory tratemp;
  GeomPlane pl1;
  int index,ist,is;
  real ll;

  bool rectangle=false;
  bool hexagon=false;
  bool s_triangle=false;

  real len1=0.;
  real len2=0.;
  if(ginp.GetGek(0)==Polygon){
    len1=ginp.GetPolygon(0).GetLen(0);
    len2=ginp.GetPolygon(0).GetLen(1);
    int tmp=ginp.GetPolygon(0).GetKind();
    if(tmp==1){ // square
      rectangle=true;
      if(angle==180.&&AngleDiv%2==1){
	cout<<"# Error in TrajectorySet::CalTrajectory.\n";
	cout<<"# You cannot choose angle=180 and anglediv=2n+1.\n";
	exit(0);
      };
    }else if(tmp==2){ // hexagon
      bool same=true;
      for(int i=1;i<6;i++){
  	real len=ginp.GetPolygon(0).GetLen(i);
	if(fabs(len1-len)/len>1e-5)same=false;
      };
      if(same)hexagon=true;
    }else if(tmp==3){ // specific-triangle
      s_triangle=true;
      len1=ginp.GetPolygon(0).GetLen(2);
      len2=ginp.GetPolygon(0).GetLen(0);
    }else if(tmp==4){ // defined by [PutTriangleHalf]
      s_triangle=true;
      len1=ginp.GetPolygon(0).GetLen(1);
      len2=ginp.GetPolygon(0).GetLen(2);
    };
  };

  /*
  if(hexagon&&BoundaryCondition==Periodic&&angle!=60.){
    cout<<"Error in TrajectorySet::CalTrajectory.\n";
    cout<<"Please set other boundary condition\n";
    cout<<" or set angle = 60.\n";
    exit(0);
  };
  */

  real UnitAngle=(angle/180.*PI)/AngleDiv; // Degree -> radian
  real weight_adjust=1.;
  if(!angle360)weight_adjust*=2.; // for opposite direction (y-reflection)
  real ang=PI/AngleDiv;

  for(int ann=0;ann<AngleDiv;ann++){

    cout<<"#  "<<ann<<" : ";
    real cosv,sinv;

    real target=UnitAngle*0.5+ann*UnitAngle;
    cosv=cos(target);
    sinv=sqrt(1.0-cosv*cosv);
    dl=dlinp;

    if(BoundaryCondition==Periodic||BoundaryCondition==Reflective){
      if(rectangle||s_triangle){
        // +++ Rectangle/Specific triangle ++++++++++++++++++++++++++++++++++++++;
        real tanv=abs(sinv/cosv);
        real diff_max=1.;
	int min_m=0;
        int min_n=0;
        for(int m=1;m<=100;m++){
	  //for(int m=1;m<=100;m++){
	  //for(int m=1;m<=1000;m++){
	  for(int n=1;n<=m;n++){
	    //if(n%2==0&&m%2==0){
            real tmp1=len2/len1;
            real tmp2=real(n)/m;
	    real ptan=tmp1*tmp2;
            real angtmp=atan(ptan);
	    if(target>PI*0.5)angtmp=PI-angtmp;
	    //real diff=fabs(ptan/tanv-1.0);
	    real diff=fabs(angtmp/target-1.0);
	    if(diff<diff_max*0.9999){
	      diff_max=diff;
	      min_m=m;
	      min_n=n;
            };
	    ptan=tmp1/tmp2;
            angtmp=atan(ptan);
	    if(target>PI*0.5)angtmp=PI-angtmp;
	    //diff=fabs(ptan/tanv-1.0);
	    diff=fabs(angtmp/target-1.0);
	    if(diff<diff_max*0.9999){
	      diff_max=diff;
	      min_m=n;
	      min_n=m;
            };
          };
	  //};
        };
        if(min_m==0||min_n==0){
	  cout<<"\n# Error in PJI_trajectoryset::CalTrejectory.\n";
	  exit(0);
	};
	real ptan=len2/len1*min_n/min_m;
	//cout<<min_m<<" "<<min_n<<" "<<diff_max<<"\n";
	cosv=sqrt(1./(1.+ptan*ptan));
	if(target>PI*0.5)cosv*=-1.;
	sinv=sqrt(1.-cosv*cosv);
	dl=len2*fabs(cosv)/min_m;
	if(dl>dlinp){
	  int itmp=int(dl/dlinp);
	  dl/=(itmp+1);
	};
	//if(dl*2.<dlinp)dl*=2.;
      };
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(hexagon){
        // +++ Hexagon
	bool over_90=false;

        // +++ for symmetrical angular setting
	//bool symmetric_angular_setting=false;
	bool symmetric_angular_setting=true;

	if(symmetric_angular_setting){
  	  if(target>PI*0.5){
	    over_90=true;
	    target=PI-target;
	  };
	};

        int iun=int(target/(PI*0.33333333333333));
        target=target-iun*PI*0.33333333333333;
	// The following part is modified in 2010/02/08
	/*
	real target_tan=sin(target)/cos(target);
	real diff_max=100.;
	int min_m=0;
        int min_n=0;
	for(int n=1;n<=100;n++){
	  for(int m=-n/2;m<=n;m++){
	    real present_tan=(2.*m+n)/(sqrt(3.0)*n);
	    real diff=fabs(present_tan/target_tan-1.0);
	    if(diff<diff_max){
	      diff_max=diff;
	      min_m=m;
	      min_n=n;
	    };
	  };
	};
        if(min_n==0){
	  cout<<"Error in PJI_trajectoryset::CalTrejectory.\n";
	  exit(0);
	};
	real present_tan=(2.*min_m+min_n)/(sqrt(3.0)*min_n);
	cosv=sqrt(1./(1.+present_tan*present_tan));
	sinv=sqrt(1.-cosv*cosv);
	dl=(len1*sqrt(3.0))*cosv/min_n;
	if(dl>dlinp){
  	  int itmp=int(dl/dlinp);
	  dl/=itmp;
	};
	*/
        real eps_theta=0.05;
        real diff_theta=1.;
        int n=1;
        real factor=0.8;
        while(diff_theta>eps_theta||dl>dlinp||dl<dlinp*factor){
	  for(int m=-n/2;m<=n;m++){
	    real present_tan=(2.*m+n)/(sqrt(3.0)*n);
            real theta1=atan(present_tan);
 	    cosv=sqrt(1./(1.+present_tan*present_tan));
	    sinv=sqrt(1.-cosv*cosv);
	    dl=(len1*sqrt(3.0))*cosv/n;
  	    if(dl>dlinp){
  	      int itmp=int(dl/dlinp);
	      if(itmp!=0)dl/=itmp;
	    };
            diff_theta=fabs(theta1/target-1.0);
            if(diff_theta<eps_theta&&dl<dlinp&&dl>dlinp*factor)break;
	  };
	  n++;
	  if(n==100)factor=0.7;
	  if(n==500)factor=0.6;
	  if(n==1000)factor=0.5;
	  if(n==2000)factor=0.4;
	  if(n==3000)factor=0.3;
	  if(n==5000)factor=0.;
	  if(n==10000){
	    cout<<"\n# Appropriate angle cannot be found.\n";
	    cout<<"# Please tune PJI_trajectoryset::CalTrajectory.\n";
	    exit(0);
	  };
	};

        real cosa=1.;
	real sina=0.;
	for(int i=0;i<iun;i++){
	  real cosanew=cosa*0.5-sina*sqrt(3.)*0.5;
	  real sinanew=sina*0.5+cosa*sqrt(3.)*0.5;
	  cosa=cosanew;
	  sina=sinanew;
	};
	real cosnewv=cosv*cosa-sinv*sina;
	real sinnewv=sinv*cosa+cosv*sina;
        cosv=cosnewv;
	sinv=sinnewv;
        target+=iun*PI*0.3333333333333;

	if(symmetric_angular_setting){
	if(over_90){
	  cosv=-cosv;
	  target=PI-target;
	};
	};

      };
    };

    // The modified following algorithm is implemented in 2009/05/28
    real xin=-edge;
    real yin=edge;
    if(target>PI*0.5)xin*=-1.;
    if(hexagon){
      yin=edge*0.5*sqrt(3.);
      xin=-edge*0.5;
      if(target>PI*0.33333333333333){
	xin=-edge;
	yin=0.;
      };
      if(target>PI*0.33333333333333*2.){
	xin=-edge*0.5;
	yin=-edge*0.5*sqrt(3.);
      };
    };

    GeomVector start_point(xin,yin);

    vec1.put_data(cosv,sinv);
    vec2.put_data(sinv,-cosv);
    if(!hexagon&&target>PI*0.5)vec2.Initialize(-sinv,cosv);

    real max_len=2.*edge*sqrt(2.)*1.01;
    if(hexagon)max_len=edge*2.*1.01;

    index=0;
    ist=NumTrajectory;
    is=0;
    ll=dl*0.5;

    while(ll<max_len){
      pl1.PutVector(start_point+vec2*ll-vec1*max_len,start_point+vec2*ll+vec1*max_len);
      //pl1.show_self();
      GeomVector vec1,vec2;
      tratemp=PlaneToTrajectory(pl1,vec1,vec2,ginp);
      //cout<<vec1.getx()<<" "<<vec1.gety()<<" "<<vec2.getx()<<" "<<vec2.gety()<<"\n";
      if(tratemp.GetExist()){
	//tratemp.PutWeight(dl*UnitAngle*weight_adjust);
	tratemp.PutWeight(dl*ang);
        line.push_back(tratemp);
	//line[NumTrajectory].show_self();
        pst.push_back(vec1);
        ped.push_back(vec2);
        NumTrajectory++;
	index++;
      };
      is++;
      ll=dl*is+dl*0.5;
    };

    real theta=atan(sinv/cosv)*180/PI;
    if(theta<0.)theta+=180.;
    cout.setf(ios::showpoint);
    cout.precision(4);
    cout<<theta<<" / ";
    cout.precision(3);
    cout<<dl<<" / "<<NumTrajectory-ist<<" / "<<NumTrajectory<<" \n";

    trajectory_on_angular_boundary[ann]=NumTrajectory-1;
    discreted_angle[ann]=theta;
    vdir[ann]=vec1;

  };

  cout<<"# Total Trajectory="<<NumTrajectory<<"\n";


  // +++ Trajectory-connection calculation on outer boundary

  NextTrajectory.resize(NumTrajectory);
  PreTrajectory.resize(NumTrajectory);

  int numt=NumTrajectory;
  if(BoundaryCondition==Periodic){
    // +++++++++ periodic boundary condition +++++++++++++++++++++
    if(ginp.GetGek(0)==Circle){
      for(int i=0;i<numt;i++){
	NextTrajectory[i]=i;
	PreTrajectory[i]=i;
      };
    }else{
      int angid=0;
      int tst=0;
      int ted=trajectory_on_angular_boundary[0];
      for(int i=0;i<numt;i++){
        int iidx=-1;
        GeomVector aop=ginp.SymmVec(ped[i]);
	real tmpx=aop.getx();
	real tmpy=aop.gety();
        for(int j=tst;j<=ted;j++){
          if(fabs(pst[j].getx()-tmpx)<0.0001&&
             fabs(pst[j].gety()-tmpy)<0.0001){
            iidx=j;
	    break;
  	  };
        };
        if(iidx<0){
          cout<<"\n# Error in CalPij!\n";
	  cout<<"# No consistency between pst and ped.\n";
	  cout<<"# The following pst cannot be found.\n";
	  aop.show_self();
	  exit(0);
        };
        NextTrajectory[i]=iidx;
        PreTrajectory[iidx]=i;
        if(i==trajectory_on_angular_boundary[angid]){
	  angid++;
	  tst=ted+1;
	  ted=trajectory_on_angular_boundary[angid];
	};
      };
    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }else if(BoundaryCondition==Reflective){
    // +++++++++ Reflective boundary condition ++++++++++++++++++
    if(ginp.GetGek(0)==Circle){
      cout<<"# Error in TrajectorySet::CalTrajectory.\n";
      cout<<"# Not coded for cylinder geometry with reflective boundary condition.\n";
      exit(0);
    }else{
      /*
      int angid=0;
      int tst=-1;
      int ted,tst2,ted2;
      for(int i=0;i<numt;i++){

	if(tst==-1){
          int angid2=AngleDiv-1-angid;
          tst=trajectory_on_angular_boundary[angid2];
	  if(angid2==0){
	    ted=0;
	  }else{
      	    ted=trajectory_on_angular_boundary[angid2-1]+1;
	  };
          tst2=tst;
          ted2=ted;
	};

        int iidx=-1;
        for(int j=tst;j>=ted;j--){
          if(ped[i].BeIdentical(pst[j],0.0001)){
            iidx=j;
            break;
	  }else if(ped[i].BeIdentical(ped[j],0.0001)){
	    iidx=numt+j;
	    break;
          };
	};
        if(iidx==-1){
          cout<<"\n# Error in CalPij!\n";
	  cout<<"# No consistency between pst and ped.\n";
	  cout<<"# The following point cannot be found.\n";
	  ped[i].show_self();
	  exit(0);
        };
        NextTrajectory[i]=iidx;

        iidx=-1;
        for(int j=tst;j>=ted;j--){
          if(pst[i].BeIdentical(pst[j],0.0001)){
            iidx=numt+j;
	    break;
	  }else if(pst[i].BeIdentical(ped[j],0.0001)){
	    iidx=j;
	    break;
	  };
        };
        if(iidx==-1){
          cout<<"\n# Error in CalPij!\n";
	  cout<<"# No consistency between pst and ped.\n";
	  cout<<"# The following point cannot be found.\n";
	  pst[i].show_self();
	  exit(0);
        };
        PreTrajectory[i]=iidx;

        if(i==trajectory_on_angular_boundary[angid]){
	  angid++;
	  tst=-1;
 	};

      };
      */

      int angid=0;
      for(int i=0;i<numt;i++){

	// Next trajectory calculation
	{
        int angid2=AngleDiv-1-angid;
	if(fabs(ped[i].getx()-ped[i].gety())<1e-6){
	  if(angid>=AngleDiv/2){
            angid2=AngleDiv/4*3+(AngleDiv/4*3-angid)-1;
	  }else{
            angid2=AngleDiv/4+(AngleDiv/4-angid)-1;
	  };
	};

        int tst=trajectory_on_angular_boundary[angid2];
        int ted=0;
	if(angid2!=0)ted=trajectory_on_angular_boundary[angid2-1]+1;

	//cout<<i<<"/"<<numt<<"  "<<angid2<<" "<<tst<<" "<<ted<<"\n";
        // 
        int iidx=-1;
        for(int j=tst;j>=ted;j--){
          if(ped[i].BeIdentical(pst[j],0.0001)){
            iidx=j;
            break;
	  }else if(ped[i].BeIdentical(ped[j],0.0001)){
	    iidx=numt+j;
	    break;
          };
	};
        if(iidx==-1){
          cout<<"\n# Error in TrajectorySet::CalTrajectory!\n";
	  cout<<"# Cannot find the next trajectory at the following point.\n";
	  ped[i].show_self();
	  cout<<"# angid : "<<angid<<"  angid2 : "<<angid2<<"\n";
	  exit(0);
        };
        NextTrajectory[i]=iidx;
	};

	// Pre trajectory calculation
	{
        int angid2=AngleDiv-1-angid;
	if(fabs(pst[i].getx()-pst[i].gety())<1e-6){
	  if(angid>=AngleDiv/2){
            angid2=AngleDiv/4*3+(AngleDiv/4*3-angid)-1;
	  }else{
            angid2=AngleDiv/4+(AngleDiv/4-angid)-1;
	  };
	};

        int tst=trajectory_on_angular_boundary[angid2];
        int ted=0;
	if(angid2!=0)ted=trajectory_on_angular_boundary[angid2-1]+1;

        int iidx=-1;
        for(int j=tst;j>=ted;j--){

	  /*
	  if(fabs(pst[i].getx()-pst[i].gety())<0.0001){
	    cout.precision(6);
	    cout<<pst[i].getx()<<" "<<ped[j].getx()<<" "<<ped[j].gety()<<"\n";
	  };
	  */

          if(pst[i].BeIdentical(pst[j],0.0001)){
            iidx=numt+j;
	    break;
	  }else if(pst[i].BeIdentical(ped[j],0.0001)){
	    iidx=j;
	    break;
	  };
        };
        if(iidx==-1){
          cout<<"\n# Error in TrajectorySet::CalTrajectory!\n";
	  cout<<"# Cannot find the pre trajectory at the following point.\n";
	  pst[i].show_self();
	  cout<<"# angid : "<<angid<<"  angid2 : "<<angid2<<"\n";
	  exit(0);
        };
        PreTrajectory[i]=iidx;
	};

        if(i==trajectory_on_angular_boundary[angid])angid++;

      };

    };
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  };

  CalVolume();
  
}


void TrajectorySet::CalTrajectoryMB(IrregularGeometryInformation &ginp,int AngleDiv,real dlinp,real angle)
{
  PutAngularDivision(AngleDiv);
  //angular_division=AngleDiv;
  //trajectory_on_angular_boundary.resize(angular_division,0);
  //discreted_angle.resize(angular_division);

  if(BoundaryCondition==Periodic){
    cout<<"Error in TrajectorySet::CalTrajectoryMB.\n";
    cout<<"MacroBand method cannot be applied to the periodic B.C.\n";
    exit(0);
  };

  cout<<"Trajectory calculation is being performed...\n";

  regnum=ginp.GetRegion();

  real edge=ginp.GetEdge();
  surface=ginp.GetSurface();

  NumTrajectory=0;

  real dl=0.;
  real llim,rlim;

  real Wangle,cosv,sinv;
  GeomVector dummy1,dummy2;
  GeomVector p1(-edge,-edge);
  GeomVector vec1,vec2,v1,v2;
  Trajectory tratemp,tt2,tt3;
  GeomPlane pl1,pl2,pl3;

  real UnitAngle=(angle/90.*0.5*PI)/AngleDiv;

  int max_divxsi=1000;
  real *divxsi=new real[max_divxsi];

  for(int ann=0;ann<AngleDiv;ann++){

    cout<<"*Angle "<<ann<<"\n";

    cosv=cos(UnitAngle*0.5+UnitAngle*ann);
    sinv=sqrt(1.0-cosv*cosv);
    Wangle=0.5*PI/AngleDiv;

    vec1.put_data(cosv,sinv);
    vec2.put_data(-sinv,cosv);
    llim=2*edge*cosv+dlinp*0.5;
    rlim=2*edge*sinv+dlinp*0.5;
    v1=p1+vec2*llim;
    v2=p1-vec2*rlim;
    GeomPlane lin(v1,v2);

    int divmb=2;
    divxsi[0]=-1.0;
    divxsi[1]=1.0;

    GeomCircle *tmpc=new GeomCircle;
    GeomVector c1,c3;
    GeomPlane cc;

    for(int i=0;i<ginp.GetIcir();i++){
      *tmpc=ginp.GetCircle(i);
      c1=tmpc->GetCenter();
      cc.PutVector(c1,c1+vec1*10);
      c3=cc.CrossPlane(lin);
      divxsi[divmb]=lin.get_xsi(c3);
      divxsi[divmb+1]=lin.get_xsi(c3+vec2*tmpc->GetR());
      divxsi[divmb+2]=lin.get_xsi(c3-vec2*tmpc->GetR());
      divmb+=3;
      if(divmb>max_divxsi)cout<<"Increase max_divxsi in CalPij.\n";
    };
    delete tmpc;

    GeomPolygon *tmph=new GeomPolygon;
    for(int i=0;i<ginp.GetIpol();i++){
      *tmph=ginp.GetPolygon(i);
      for(int j=0;j<tmph->GetNumpl();j++){
	c1=tmph->GetPlane(j).get_p1();
	cc.PutVector(c1,c1+vec1*10);
	c3=cc.CrossPlane(lin);
	divxsi[divmb]=lin.get_xsi(c3);
        divmb++;
        if(divmb>max_divxsi)cout<<"Increase max_divxsi in CalPij.\n";
      };
    };
    delete tmph;

    int *divo=new int[divmb];
    ChangeOrder(divxsi,divmb,divo);

    GeomVector tmp,tmp2,tmp3;
    for(int i=0;i<divmb-1;i++){
      real divz=divxsi[divo[i+1]]-divxsi[divo[i]];
      if(fabs(divz)>0.0000001){
        real ll=divz*lin.get_long()*0.5;
        int div=int(ll/dlinp)+1;
        dl=ll/div;
        for(int j=0;j<div;j++){
	  real xsi1=divxsi[divo[i]]+divz/div*0.5*(j*2+1);
	  real xsi2=xsi1-divz/div*0.33333333;
	  real xsi3=xsi1+divz/div*0.33333333;
	  tmp=lin.GetVecXsi(xsi1);
	  tmp2=lin.GetVecXsi(xsi2);
	  tmp3=lin.GetVecXsi(xsi3);
	  pl1.PutVector(tmp,tmp+vec1*3*edge);
	  pl2.PutVector(tmp2,tmp2+vec1*3*edge);
	  pl3.PutVector(tmp3,tmp3+vec1*3*edge);
          tratemp=PlaneToTrajectory(pl1,dummy1,dummy2,ginp);
          tt2=PlaneToTrajectory(pl2,dummy1,dummy2,ginp);
          tt3=PlaneToTrajectory(pl3,dummy1,dummy2,ginp);
          tratemp.AdjustTrajectoryLength(tt2,tt3);
          if(tratemp.GetExist()){
	    tratemp.PutWeight(dl*Wangle*2);  //"2" means (pi/2 ~ pi)
            line.push_back(tratemp);
            NumTrajectory++;
	  };
	};
      };
    };

    real theta=atan(sinv/cosv)*180/PI;
    if(theta<0.)theta+=180.;
    cout<<"    ("<<theta<<" & "<<dl<<") "<<NumTrajectory<<"\n";

    trajectory_on_angular_boundary[ann]=NumTrajectory-1;
    discreted_angle[ann]=theta;

    delete [] divo;
  };

  cout<<"Total Trajectory="<<NumTrajectory<<"\n";

  CalVolume();

  delete [] divxsi;
}

Trajectory TrajectorySet::PlaneToTrajectory(GeomPlane &pl1,GeomVector &pst, GeomVector &ped, IrregularGeometryInformation &gset)
{
  int no_segment;
  vector<real> length_segment;
  vector<int> regid_segment;
  real xsi_stt,xsi_end;

  gset.DrawTrajectory(pl1,xsi_stt,xsi_end,no_segment,length_segment,regid_segment);

  if(no_segment==0){
    Trajectory ret;
    return ret;
  };

  pst=(pl1.get_p1()*(1-xsi_stt)+pl1.get_p2()*(1+xsi_stt))*0.5;
  ped=(pl1.get_p1()*(1-xsi_end)+pl1.get_p2()*(1+xsi_end))*0.5;

  Trajectory rettra;
  rettra.Initialize(no_segment);
  rettra.PutData(regid_segment,length_segment);
  return rettra;
};

void TrajectorySet::PutBoundaryCondition(BCondition i)
{
  if(i!=Black&&i!=White&&i!=Periodic&&i!=Reflective){
    cout<<"Boundary Condition is not appropriated!\n";
    exit(0);
  };
  BoundaryCondition=i;
};

bool TrajectorySet::IsFixedBoundaryCondition()
{
  if(BoundaryCondition==Black)return true;
  if(BoundaryCondition==White)return true;
  if(BoundaryCondition==Periodic)return true;
  return false;
};

void TrajectorySet::WriteFile(string mdir,string ss)
{
  cout<<"# Writing the TrajectorySet data on file ["<<ss<<"]\n";
  mdir.append(ss);

  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  fout.setf(ios::scientific);
  fout.precision(12);
  fout<<NumTrajectory<<"\n";
  fout<<regnum<<"\n";

  fout<<angular_division<<"\n";
  for(int i=0;i<angular_division;i++){
    fout<<trajectory_on_angular_boundary[i]<<"\n";
    fout<<discreted_angle[i]<<"\n";
  };

  fout<<surface<<"\n";

  for(int i=0;i<NumTrajectory;i++){
    // trajectorydata
    int tmp=line[i].GetNumreg();
    fout<<tmp<<"\n";
    fout<<line[i].GetWeight()<<"\n";
    for(int j=0;j<tmp;j++){
      fout<<line[i].GetReg(j)<<"\n";
      fout<<line[i].GetDist(j)<<"\n";
    };
    fout<<NextTrajectory[i]<<"\n";
    fout<<PreTrajectory[i]<<"\n";
  };

  for(int i=0;i<regnum;i++){
    fout<<volume[i]<<"\n";
  };

  fout.close();
};

void TrajectorySet::ReadFile(string mdir,string ss)
{
  cout<<"# Reading the TrajectorySet data from file ["<<ss<<"]\n";
  mdir.append(ss);

  ifstream fin;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  int tmp;
  real rtmp;

  fin>>tmp;
  NumTrajectory=tmp;
  line.resize(NumTrajectory);
  NextTrajectory.resize(NumTrajectory);
  PreTrajectory.resize(NumTrajectory);

  fin>>tmp;
  regnum=tmp;
  volume.resize(regnum);

  fin>>tmp;
  PutAngularDivision(tmp);

  for(int i=0;i<angular_division;i++){
    fin>>tmp;
    trajectory_on_angular_boundary[i]=tmp;
    fin>>rtmp;
    discreted_angle[i]=rtmp;
  };

  fin>>rtmp;
  surface=rtmp;

  for(int i=0;i<NumTrajectory;i++){
    // trajectorydata
    fin>>tmp;
    int nn=tmp;
    line[i].Initialize(tmp);
    fin>>rtmp;
    line[i].PutWeight(rtmp);
    vector<int> ireg(tmp);
    vector<real> idist(tmp);
    real sum=0.;
    for(int j=0;j<nn;j++){
      fin>>tmp;
      ireg[j]=tmp;
      fin>>rtmp;
      sum+=rtmp;
      idist[j]=sum;
    };
    line[i].PutData(ireg,idist);
    fin>>tmp;
    NextTrajectory[i]=tmp;
    fin>>tmp;
    PreTrajectory[i]=tmp;
  };

  for(int i=0;i<regnum;i++){
    fin>>rtmp;
    volume[i]=rtmp;
  };


  fin.close();
};

void TrajectorySet::LengthFactorize(real factor)
{
  for(int i=0;i<NumTrajectory;i++){
    line[i].LengthFactorize(factor);
  };
};

int TrajectorySet::GetMaximumSegments()
{
  int max_segments=0;
  for(int i=0;i<NumTrajectory;i++){
    if(line[i].GetNumreg()>max_segments)max_segments=line[i].GetNumreg();
  };
  return max_segments;
};
