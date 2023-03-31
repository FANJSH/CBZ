#include <cstdlib>
#include "HigherModeCalculator.h"

HigherModeCalculator::HigherModeCalculator()
{
  iter_max = 500;
  epsk = 1e-6;
  xsss=15;
  ysss=15;
  zsss=0;
  accel=true;
};

void HigherModeCalculator::PutInitialSourcePosition(int x,int y,int z)
{
  xsss=x;
  ysss=y;
  zsss=z;
};

void HigherModeCalculator::RunDHEX
(vector<DHEXSystem> &sysf, vector<DHEXSystem> &sysa,
 vector<real> &kf, int max_order, bool cor_tri)
{
  int group=sysf[0].GetGrp();

  vector<real> ipnn(max_order);

  int totm=sysf[0].GetTotM();
  vector<bool> fissile_mesh(totm,false);
  for(int i=0;i<totm;i++){
    for(int g=0;g<group;g++){
      if(sysf[0].GetMesh(i).GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g)>0.)fissile_mesh[i]=true;
    };
  };

  // +++ forward fundamental mode
  if(cor_tri)sysf[0].PutTriangularMeshCorrection();
  sysf[0].CalCoef();
  if(cor_tri){
    kf[0]=sysf[0].CalIgen();
  }else{
    kf[0]=sysf[0].CalIgen("cmfd");
  };
  //sysf[0].MemoryReductionForHighMomentCal();
  // +++ adjoint fundamental mode
  if(cor_tri)sysa[0].PutTriangularMeshCorrection();
  sysa[0].CopyCoef(sysf[0]);
  sysa[0].CalIgen();
  sysa[0].MemoryReductionForHighMomentCal();

  ipnn[0]=sysa[0].CalPerturbDenominatorFromSource(&sysf[0]);

  // (For source extrapolation)
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold=0.;
  int expit=0;
  bool exp=false;
  //**

  for(int od=1;od<max_order;od++){

  // forward  
  {
  cout<<"\n +++ Forward calculation ("<<od<<"-th mode +++)\n\n";
  if(cor_tri)sysf[od].PutTriangularMeshCorrection();
  sysf[od].CopyCoef(sysf[0]);
  for(int g=0;g<group;g++){
    sysf[od].GetMesh(xsss,ysss,zsss).GetFlux().put_data(g,1.);
  };
  sysf[od].NormalizeFissionSrc();

  real fiss=1.;
  real fissold;
  real errk;
  lamda=0.;
  expit=0;
  exp=false;

  for(int iter=0;iter<iter_max;iter++){

    for(int g=0;g<group;g++){
      sysf[od].CalSrcMultiplySystem(g,fiss,0);
      sysf[od].CalFluxGeneral(g,1e-2,iter);
      if(g!=group-1)sysf[od].AddDownScatSrc(g,0);
      sysf[od].SetZeroScatSrc(g);
    };

    vector<real> ip2(od);
    for(int od2=0;od2<od;od2++){
      ip2[od2]=sysa[od2].CalPerturbDenominatorFromSource(&sysf[od])/ipnn[od2];
    };

    lamdaold=lamda;
    lamda=0.;
    int ind=0;

    real sum=0.;
    for(int i=0;i<totm;i++){
      if(fissile_mesh[i]){
        real tmp=sysf[od].GetMesh(i).GetFissionSrc();
        fsold[i]=tmp;
        sysf[od].GetMesh(i).CalFissionSrc();
        real tmp2=sysf[od].GetMesh(i).GetFissionSrc();
        for(int od2=0;od2<od;od2++){
          real tmp3=sysf[od2].GetMesh(i).GetFissionSrc();
          tmp2-=tmp3*ip2[od2];
        };
        sysf[od].GetMesh(i).PutFissionSrc(tmp2);
        real tmp4=tmp2-tmp;
	real tmp5=tmp4/res[i];
	res[i]=tmp4;
	if(tmp5>0.&&tmp5<1.){
	  lamda+=tmp5;
	  ind++;
	};
        sum+=fabs(tmp2);
      };
    };

    lamda/=ind;
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>0){
      real omega=1./(1.-lamda);
      exp=true;
      cout<<"   (Extrapolation for fission source)\n";
      cout<<"     dominance ratio="<<lamda<<" , omega = "<<omega<<"\n";
      sum=0.;
      for(int i=0;i<totm;i++){
	if(fissile_mesh[i]){
	  real tmp=fsold[i]+omega*(sysf[od].GetMesh(i).GetFissionSrc()-fsold[i]);
	  sysf[od].GetMesh(i).PutFissionSrc(tmp);
	  sum+=fabs(sysf[od].GetMesh(i).GetFissionSrc());
	};
      };
    }else if(exp){
      expit=0;
      exp=false;
    };

    fissold=fiss;
    fiss=sum;
    errk=fabs(fiss-fissold)/fissold;
    sysf[od].WriteIterationInfo(iter,fiss,errk,0.0,0.0);
    if(errk<epsk)break;
  };
  kf[od]=fiss;
  sysf[od].MemoryReductionForHighMomentCal();
  };

  // adjoint
  if(od!=max_order-1)
  {
  cout<<"\n +++ Adjoint calculation ("<<od<<"-th mode) +++\n\n";
  if(cor_tri)sysa[od].PutTriangularMeshCorrection();
  sysa[od].CopyCoef(sysf[0]);
  for(int g=0;g<group;g++){
    sysa[od].GetMesh(xsss,ysss,zsss).GetFlux().put_data(g,1.);
  };
  sysa[od].NormalizeFissionSrc();

  real fiss=1.;
  real fissold;
  lamda=0.;
  expit=0;
  exp=false;
  for(int iter=0;iter<iter_max;iter++){

    for(int g=0;g<group;g++){
      int ginp=group-1-g;      
      sysa[od].CalSrcMultiplySystem(ginp,fiss,0);
      sysa[od].CalFluxGeneral(ginp,1e-2,iter);
      if(g!=group-1)sysa[od].AddDownScatSrc(ginp,0);
      sysa[od].SetZeroScatSrc(ginp);
    };

    vector<real> ip2(od);
    for(int od2=0;od2<od;od2++){
      ip2[od2]=sysa[od].CalPerturbDenominatorFromSource(&sysf[od2])/ipnn[od2];
    };

    lamdaold=lamda;
    lamda=0.;
    int ind=0;

    real sum=0.;
    for(int i=0;i<totm;i++){
      if(fissile_mesh[i]){
        real tmp=sysa[od].GetMesh(i).GetFissionSrc();
        fsold[i]=tmp;
        sysa[od].GetMesh(i).CalFissionSrcAdjoint();
        real tmp2=sysa[od].GetMesh(i).GetFissionSrc();
        for(int od2=0;od2<od;od2++){
	  real tmp3=sysa[od2].GetMesh(i).GetFissionSrc();
	  tmp2-=tmp3*ip2[od2];
        };
        sysa[od].GetMesh(i).PutFissionSrc(tmp2);
        real tmp4=tmp2-tmp;
        real tmp5=tmp4/res[i];
	res[i]=tmp4;
	if(tmp5>0.&&tmp5<1.){
	  lamda+=tmp5;
	  ind++;
	};
        sum+=fabs(tmp2);
      };
    };
    lamda/=ind;
    expit++;
    if(fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>0){
      real omega=1./(1.-lamda);
      exp=true;
      cout<<"   (Extrapolation for fission source)\n";
      cout<<"     dominance ratio="<<lamda<<" , omega = "<<omega<<"\n";
      sum=0.;
      for(int i=0;i<totm;i++){
	if(fissile_mesh[i]){
	  real tmp=fsold[i]+omega*(sysa[od].GetMesh(i).GetFissionSrc()-fsold[i]);
	  sysa[od].GetMesh(i).PutFissionSrc(tmp);
	  sum+=fabs(sysa[od].GetMesh(i).GetFissionSrc());
	};
      };
    }else if(exp){
      expit=0;
      exp=false;
    };

    fissold=fiss;
    fiss=sum;
    real errk=fabs(fiss-fissold)/fissold;
    sysa[od].WriteIterationInfo(iter,fiss,errk,0.0,0.0);
    if(errk<epsk)break;
  };
  real dif_k=fabs(fiss/kf[od]-1.);
  if(dif_k>0.005){
    cout<<"Error !\n";
    cout<<"Large difference is observed between forward k and adjoint k.\n";
    exit(0);
  };
  };
  sysa[od].MemoryReductionForHighMomentCal();
  ipnn[od]=sysa[od].CalPerturbDenominatorFromSource(&sysf[od]);
  };

};

void HigherModeCalculator::RunPLOS
(vector<PLOSSystem> &sysf, vector<PLOSSystem> &sysa,
 vector<real> &kf, int max_order)
{
  bool memory_reduction=true;

  // The following is the same as RunDHEX.

  int group=sysf[0].GetGrp();

  vector<real> ipnn(max_order);

  int totm=sysf[0].GetTotM();
  vector<bool> fissile_mesh(totm,false);
  for(int i=0;i<totm;i++){
    for(int g=0;g<group;g++){
      if(sysf[0].GetMesh(i).GetMed()->GetMacxs().GetData1d(nusigf).get_dat(g)>0.)fissile_mesh[i]=true;
    };
  };

  // +++ forward fundamental mode
  sysf[0].CalCoef();
  kf[0]=sysf[0].CalIgen("cmfd");
  //sysf[0].MemoryReductionForHighMomentCal();
  // +++ adjoint fundamental mode
  sysa[0].CopyCoef(sysf[0]);
  sysa[0].CalIgen();
  if(memory_reduction)sysa[0].MemoryReductionForHighMomentCal();

  ipnn[0]=sysa[0].CalPerturbDenominatorFromSource(&sysf[0]);

  // (For source extrapolation)
  vector<real> fsold(totm);
  vector<real> res(totm);
  real lamda=0.;
  real lamdaold=0.;
  int expit=0;
  bool exp=false;
  //**

  for(int od=1;od<max_order;od++){

  // forward  
  {
  sysf[od].CopyCoef(sysf[0]);
  for(int g=0;g<group;g++){
    sysf[od].GetMesh(xsss,ysss,zsss).GetFlux().put_data(g,1.);
  };
  sysf[od].NormalizeFissionSrc();

  real fiss=1.;
  real fissold;
  real errk;
  cout<<"\n# +++ Forward calculation ("<<od<<"-th mode +++)\n#\n";
  lamda=0.;
  expit=0;
  exp=false;

  for(int iter=0;iter<iter_max;iter++){

    for(int g=0;g<group;g++){
      sysf[od].CalSrcMultiplySystem(g,fiss,0);
      sysf[od].CalFluxGeneral(g,1e-2,iter);
      if(g!=group-1)sysf[od].AddDownScatSrc(g,0);
      sysf[od].SetZeroScatSrc(g);
    };

    vector<real> ip2(od);
    for(int od2=0;od2<od;od2++){
      ip2[od2]=sysa[od2].CalPerturbDenominatorFromSource(&sysf[od])/ipnn[od2];
    };

    lamdaold=lamda;
    lamda=0.;
    int ind=0;

    real sum=0.;
    for(int i=0;i<totm;i++){
      if(fissile_mesh[i]){
        real tmp=sysf[od].GetMesh(i).GetFissionSrc();
        fsold[i]=tmp;
        sysf[od].GetMesh(i).CalFissionSrc();
        real tmp2=sysf[od].GetMesh(i).GetFissionSrc();
        for(int od2=0;od2<od;od2++){
          real tmp3=sysf[od2].GetMesh(i).GetFissionSrc();
          tmp2-=tmp3*ip2[od2];
        };
        sysf[od].GetMesh(i).PutFissionSrc(tmp2);
        real tmp4=tmp2-tmp;
	real tmp5=tmp4/res[i];
	res[i]=tmp4;
	if(tmp5>0.&&tmp5<1.){
	  lamda+=tmp5;
	  ind++;
	};
        sum+=fabs(tmp2);
      };
    };

    lamda/=ind;
    expit++;
    if(accel&&fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>0){
      real omega=1./(1.-lamda);
      exp=true;
      cout<<"#   (Extrapolation for fission source)\n";
      cout<<"#     dominance ratio="<<lamda<<" , omega = "<<omega<<"\n";
      sum=0.;
      for(int i=0;i<totm;i++){
	if(fissile_mesh[i]){
	  real tmp=fsold[i]+omega*(sysf[od].GetMesh(i).GetFissionSrc()-fsold[i]);
	  sysf[od].GetMesh(i).PutFissionSrc(tmp);
	  sum+=fabs(sysf[od].GetMesh(i).GetFissionSrc());
	};
      };
    }else if(exp){
      expit=0;
      exp=false;
    };

    fissold=fiss;
    fiss=sum;
    errk=fabs(fiss-fissold)/fissold;
    sysf[od].WriteIterationInfo(iter,fiss,errk,0.0,0.0);
    if(errk<epsk)break;
  };
  kf[od]=fiss;
  if(memory_reduction)sysf[od].MemoryReductionForHighMomentCal();
  };

  // adjoint
  {
  sysa[od].CopyCoef(sysf[0]);
  for(int g=0;g<group;g++){
    sysa[od].GetMesh(xsss,ysss,zsss).GetFlux().put_data(g,1.);
  };
  sysa[od].NormalizeFissionSrc();

  real fiss=1.;
  real fissold;
  cout<<"\n# +++ Adjoint calculation ("<<od<<"-th mode) +++\n#\n";
  lamda=0.;
  expit=0;
  exp=false;
  for(int iter=0;iter<iter_max;iter++){

    for(int g=0;g<group;g++){
      int ginp=group-1-g;      
      sysa[od].CalSrcMultiplySystem(ginp,fiss,0);
      sysa[od].CalFluxGeneral(ginp,1e-2,iter);
      if(g!=group-1)sysa[od].AddDownScatSrc(ginp,0);
      sysa[od].SetZeroScatSrc(ginp);
    };

    vector<real> ip2(od);
    for(int od2=0;od2<od;od2++){
      ip2[od2]=sysa[od].CalPerturbDenominatorFromSource(&sysf[od2])/ipnn[od2];
    };

    lamdaold=lamda;
    lamda=0.;
    int ind=0;

    real sum=0.;
    for(int i=0;i<totm;i++){
      if(fissile_mesh[i]){
        real tmp=sysa[od].GetMesh(i).GetFissionSrc();
        fsold[i]=tmp;
        sysa[od].GetMesh(i).CalFissionSrcAdjoint();
        real tmp2=sysa[od].GetMesh(i).GetFissionSrc();
        for(int od2=0;od2<od;od2++){
	  real tmp3=sysa[od2].GetMesh(i).GetFissionSrc();
	  tmp2-=tmp3*ip2[od2];
        };
        sysa[od].GetMesh(i).PutFissionSrc(tmp2);
        real tmp4=tmp2-tmp;
        real tmp5=tmp4/res[i];
	res[i]=tmp4;
	if(tmp5>0.&&tmp5<1.){
	  lamda+=tmp5;
	  ind++;
	};
        sum+=fabs(tmp2);
      };
    };

    lamda/=ind;
    expit++;
    if(accel&&fabs(lamda-lamdaold)<0.01&&expit>2&&lamda<1.0&&!exp&&iter>0){
      real omega=1./(1.-lamda);
      exp=true;
      cout<<"#   (Extrapolation for fission source)\n";
      cout<<"#     dominance ratio="<<lamda<<" , omega = "<<omega<<"\n";
      sum=0.;
      for(int i=0;i<totm;i++){
	if(fissile_mesh[i]){
	  real tmp=fsold[i]+omega*(sysa[od].GetMesh(i).GetFissionSrc()-fsold[i]);
	  sysa[od].GetMesh(i).PutFissionSrc(tmp);
	  sum+=fabs(sysa[od].GetMesh(i).GetFissionSrc());
	};
      };
    }else if(exp){
      expit=0;
      exp=false;
    };

    fissold=fiss;
    fiss=sum;
    real errk=fabs(fiss-fissold)/fissold;
    sysa[od].WriteIterationInfo(iter,fiss,errk,0.0,0.0);
    if(errk<epsk)break;
  };
  };
  if(memory_reduction)sysa[od].MemoryReductionForHighMomentCal();
  ipnn[od]=sysa[od].CalPerturbDenominatorFromSource(&sysf[od]);
  };

};

/*
real HigherModeCalculator::CalPerturbDenominatorFromSource(GeneralSystem *sec)
{
  real denom=0.;
  for(int r=0;r<TotM;r++){
    real vol=mesh[r].GetVolume();
    denom+=mesh[r].GetFissionSrc()*sec->GetMesh(r).GetFissionSrc()/vol;
  };
  return denom;
};
*/
