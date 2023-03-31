#include<cstdlib>
#include "KineticUtility.h"

using namespace std;

OnePointSolverExplicitFPModel::OnePointSolverExplicitFPModel(bool cinp)
{
  chain_effect=cinp;
};

GroupData1D OnePointSolverExplicitFPModel::GetNfwd(int i)
{
  int time_step=second.size();
  if(i<time_step)return n_fwd[i];
};

void OnePointSolverExplicitFPModel::DataInputFromChainData(BCGManager &man,string tagname)  
{
  int sz=man.GetNuclideNumber();
  pos.resize(sz);
  for(int i=0;i<sz;i++){
    int an=man.GetNuclide(i).GetAtomicNumber();
    int ms=man.GetNuclide(i).GetMassNumber();
    int ex=man.GetNuclide(i).GetExLevel();
    int ch=man.GetNuclide(i).GetChannel();
    real hl=man.GetNuclide(i).GetHalflife();
    real lambda=0.;
    if(hl!=0.)lambda=log(2.)/hl;
    real yld=man.GetNuclide(i).GetYield(tagname);
    real dntt=man.GetNuclide(i).GetTotalNumberEmittedNeutrons();
    real dnbr=0.; // the number of neutrons emitted per one decay
    bool parent=false;
    for(int j=0;j<ch;j++){
      real dn=man.GetNuclide(i).GetEmittedNeutron(j);
      real br=man.GetNuclide(i).GetBr(j);
      dnbr+=dn*br;  
      int an2=man.GetNuclide(i).GetAtomicNumberNext(j);
      int ms2=man.GetNuclide(i).GetMassNumberNext(j);
      int ex2=man.GetNuclide(i).GetExLevelNext(j);
      int id=man.GetNuclideIndex(an2,ms2,ex2);
      if(id!=-1){
        int ch2=man.GetNuclide(id).GetChannel();
        for(int k=0;k<ch2;k++){
  	  if(man.GetNuclide(id).GetEmittedNeutron(k)>0.)parent=true;
        };
      };
    };
    if((yld*dnbr>0||parent)&&an>3){ // H, He, Li are abondaned
      //if(yld*dnbr>0){
      string namtmp=man.SearchNameFromParameter(an,ms,ex);
      name.push_back(namtmp);
      matid.push_back(an*10000+ms*10+ex);
      atomn.push_back(an);
      massn.push_back(ms);
      level.push_back(ex);
      lmd.push_back(lambda);
      dnt.push_back(dnbr);
      dnt_dtr.push_back(dntt);
      pos[i]=name.size()-1;
      pos2.push_back(i);
    };
  };
};

void OnePointSolverExplicitFPModel::DataInputFromChainDataAll(BCGManager &man,string tagname)
{
  int sz=man.GetNuclideNumber();
  pos.resize(sz);
  for(int i=0;i<sz;i++){
    int an=man.GetNuclide(i).GetAtomicNumber();
    int ms=man.GetNuclide(i).GetMassNumber();
    int ex=man.GetNuclide(i).GetExLevel();
    int ch=man.GetNuclide(i).GetChannel();
    real hl=man.GetNuclide(i).GetHalflife();
    real lambda=0.;
    if(hl!=0.)lambda=log(2.)/hl;
    real yld=man.GetNuclide(i).GetYield(tagname);
    real dntt=man.GetNuclide(i).GetTotalNumberEmittedNeutrons();
    real dnbr=0.; // the number of neutrons emitted per one decay
    for(int j=0;j<ch;j++){
      real dn=man.GetNuclide(i).GetEmittedNeutron(j);
      real br=man.GetNuclide(i).GetBr(j);
      dnbr+=dn*br;  
    };
    if(an>3){ // H, He, Li are abondaned
      //if(lambda>0.&&an>3){ // H, He, Li are abondaned
      //if(yld*dnbr>0){
      string namtmp=man.SearchNameFromParameter(an,ms,ex);
      name.push_back(namtmp);
      matid.push_back(an*10000+ms*10+ex);
      atomn.push_back(an);
      massn.push_back(ms);
      level.push_back(ex);
      lmd.push_back(lambda);
      dnt.push_back(dnbr);
      dnt_dtr.push_back(dntt);
      pos[i]=name.size()-1;
    };
  };
};

void OnePointSolverExplicitFPModel::ChainModification(BCGManager &man)
{
  man.MakeFlagFalseAllNuclide();
  int num=atomn.size();
  for(int i=0;i<num;i++){
    man.MakeFlagTrue(atomn[i],massn[i],level[i]);
  };

  if(chain_effect)man.CalCumulativeYield(); 
  // +++
  //  When chain effect is ignored, cumulative yield data is originally used.
  //  So `CalCumulativeYield' method is unnecessary.
  man.CalDecayBranch();
};

void OnePointSolverExplicitFPModel::YieldDataSetting(BCGManager &man, string tagname)
{
  //int ii=0;
  int sz=man.GetNuclideNumber();
  for(int i=0;i<sz;i++){
    int an=man.GetNuclide(i).GetAtomicNumber();
    int ms=man.GetNuclide(i).GetMassNumber();
    int ex=man.GetNuclide(i).GetExLevel();
    if(man.GetNuclide(i).Flag()){
      yield.push_back(man.GetNuclide(i).GetYield(tagname));
      delta_yield.push_back(man.GetNuclide(i).GetDeltaYield(tagname));
      /*
      int ch=man.GetNuclide(i).GetChannel();
      real dnbr=0.;
      for(int j=0;j<ch;j++){
	dnbr+=man.GetNuclide(i).GetBr(j)*man.GetNuclide(i).GetEmittedNeutron(j);
      };
      dnt[ii]=dnbr;
      ii++;
      */
    };
  };

  if(!chain_effect){
    real tot=0.;
    int sz2=name.size();
    for(int i=0;i<sz2;i++){
      tot+=yield[i]*dnt[i];
    };
    cout.setf(ios::showpoint);
    cout.precision(8);
    cout<<"# Total number of delayed neutrons : "<<tot<<"\n";
  };

  nu_d=man.CalTotalDelayedNeutron(tagname);  
};

GroupData2D OnePointSolverExplicitFPModel::CalDecayMatrix(BCGManager &man)
{
  int num=name.size();
  int num1=num+1;

  GroupData2D tmat(num1,num1);
  tmat.set_zero();

  tmat.put_data(0,0,(rho-beta)/ge);
  for(int i=0;i<num;i++){
    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    tmat.put_data(0,i+1,lmd[i]*dnt[i]);
    tmat.put_data(i+1,0,yield[i]/(nu*ge));
    tmat.put_data(i+1,i+1,-lmd[i]);
    int ch=man.GetNuclide(at,ms,lv).GetChannel();
    for(int j=0;j<ch;j++){
      int at2=man.GetNuclide(at,ms,lv).GetAtomicNumberNext(j);
      int ms2=man.GetNuclide(at,ms,lv).GetMassNumberNext(j);
      int lv2=man.GetNuclide(at,ms,lv).GetExLevelNext(j);
      real br=man.GetNuclide(at,ms,lv).GetBr(j);
      int id=man.GetNuclideIndex(at2,ms2,lv2);
      int idi=pos[id];
      if(chain_effect)tmat.put_data(idi+1,i+1,lmd[i]*br);
    };
  };

  return tmat;
};

GroupData2D OnePointSolverExplicitFPModel::CalDecayMatrixKp(BCGManager &man)
{
  int num=name.size();
  int num1=num+1;

  GroupData2D tmat(num1,num1);
  tmat.set_zero();

  tmat.put_data(0,0,(k_p-1.)/l);
  for(int i=0;i<num;i++){
    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    tmat.put_data(0,i+1,lmd[i]*dnt[i]);
    tmat.put_data(i+1,0,k_p/l*yield[i]/nu_p);
    tmat.put_data(i+1,i+1,-lmd[i]);
    int ch=man.GetNuclide(at,ms,lv).GetChannel();
    for(int j=0;j<ch;j++){
      int at2=man.GetNuclide(at,ms,lv).GetAtomicNumberNext(j);
      int ms2=man.GetNuclide(at,ms,lv).GetMassNumberNext(j);
      int lv2=man.GetNuclide(at,ms,lv).GetExLevelNext(j);
      real br=man.GetNuclide(at,ms,lv).GetBr(j);
      int id=man.GetNuclideIndex(at2,ms2,lv2);
      if(id!=-1){
        int idi=pos[id];
        if(chain_effect)tmat.put_data(idi+1,i+1,lmd[i]*br);
      };
    };
  };

  return tmat;
};

void OnePointSolverExplicitFPModel::YieldFactorize(real fact)
{
  int sz=yield.size();
  for(int i=0;i<sz;i++){
    yield[i]*=fact;
    delta_yield[i]*=fact;
  };

  nu_d*=fact;
};

void OnePointSolverExplicitFPModel::SetDefaultInitialCondition()
{
  second.clear();
  n_fwd.clear();
  k_p_hist.clear();

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val(num1);
  val.set_zero();

  real n=1.; // initial value for neutron
  val.put_data(0,n);
  for(int i=0;i<num;i++){
    if(lmd[i]!=0.){
      val.put_data(i+1,yield[i]/nu_p/lmd[i]/l*k_p);
    }else{
      val.put_data(i+1,0.);
    };
  };

  n_fwd.push_back(val);
  second.push_back(0.);
  k_p_hist.push_back(1.); // dummy
};

void OnePointSolverExplicitFPModel::SetStationaryStateAsInitialCondition(BCGManager &man, real k_p_inp)
{
  real k_p_org=k_p;
  PutK_p(k_p_inp);
  
  real delta_t=1.; // [sec]
  int stepmax=1000;
  real period_condition=1e5;

  SetDefaultInitialCondition();
  
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating
  tmat.LUDecompositionForMatrixExponentialByMMPA32(delta_t);
 
  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val=n_fwd[0];

  Burnup bu;

  real t=0;
  real n_old;
  real n=val.get_dat(0);
  for(int i=0;i<stepmax;i++){
    n_old=n;
    val=tmat.CalMatrixExponentialByLUDMMPA32(val);
    n=val.get_dat(0);
    real lmd=log(n/n_old)/delta_t;
    real period=1./lmd;
    t+=delta_t;
    if(fabs(period)>period_condition)break;
    if(i==stepmax-1){
      cout<<"# Error in OnePointSolverExplicitFPModel::SetStationaryStateAsInitialCondition.\n";
      cout<<"# Iteration reaches maximum.\n";
      cout<<"# The current period is "<<period<<"\n";
      exit(0);
    };
  };

  n_fwd[0]=val;
  PutK_p(k_p_org);
  k_p_hist[0]=k_p_org;
};


void OnePointSolverExplicitFPModel::ShowNumberOfNeutrons()
{


  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Time[s]  Neutrons   Period\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  
  int time_step=second.size();
  for(int i=0;i<time_step;i++){
    cout<<second[i]<<" "<<n_fwd[i].get_dat(0);

    // ... reactor power
    //cout<<" "<<kappa*k_p_hist[i]/nu_p/l*n_fwd[i].get_dat(0);
    
    if(i!=0){
      real delta_t=second[i]-second[i-1];
      real lmd=log(n_fwd[i].get_dat(0)/n_fwd[i-1].get_dat(0))/delta_t;
      real period=1./lmd;
      cout<<" "<<period;
    };
    cout<<"\n";
  };
};
  
void OnePointSolverExplicitFPModel::ShowDecayHeat(BCGManager &man)
{
  real kappa=180e6; // eV
  
  int num=GetTreatedNuclideNumber();
  vector<int> tmpid(num);
  for(int i=0;i<num;i++){
    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    tmpid[i]=man.GetNuclideIndex(at,ms,lv);
  };
  
  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Time      Heat energy [eV]\n";
  cout<<"#  [s]      Prompt     Delayed    (beta)     (gamma)    (alpha)\n";
  cout.setf(ios::scientific);
  cout.precision(4);


  int time_step=second.size();
  for(int i=0;i<time_step;i++){
    real dh_beta=0.;
    real dh_gamma=0.;
    real dh_alpha=0.;
    for(int j=0;j<num;j++){
      int id=tmpid[j];
      real bq=lmd[j]*n_fwd[i].get_dat(1+j);
      if(bq>0.){
	dh_beta+=man.GetNuclide(id).GetDecayEnergy(0)*bq;
	dh_gamma+=man.GetNuclide(id).GetDecayEnergy(1)*bq;
	dh_alpha+=man.GetNuclide(id).GetDecayEnergy(2)*bq;	
      };
    };
    real dh_sum=dh_beta+dh_gamma+dh_alpha;
    cout<<second[i]<<" ";
    cout<<" "<<kappa*k_p_hist[i]/nu_p/l*n_fwd[i].get_dat(0); // ... reactor power
    cout<<" "<<dh_sum;
    cout<<" "<<dh_beta<<" "<<dh_gamma<<" "<<dh_alpha;
    cout<<"\n";
  };
};

  
void OnePointSolverExplicitFPModel::CalTransientNew(BCGManager &man,int step,real delta_t)  
{
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating
  tmat.LUDecompositionForMatrixExponentialByMMPA32(delta_t);

  Burnup bu;

  int time_step=second.size();
  real t=second[time_step-1];
  for(int i=0;i<step;i++){
    GroupData1D val=tmat.CalMatrixExponentialByLUDMMPA32(n_fwd[time_step-1+i]);
    t+=delta_t;
    second.push_back(t);
    n_fwd.push_back(val);
    k_p_hist.push_back(k_p);
  };
};

void OnePointSolverExplicitFPModel::Restarting()
{
  int time_step=second.size();
  GroupData1D n_fwd_tmp=n_fwd[time_step-1];
  real k_p_hist_tmp=k_p_hist[time_step-1];

  second.clear();
  n_fwd.clear();
  k_p_hist.clear();

  second.push_back(0.);
  n_fwd.push_back(n_fwd_tmp);
  k_p_hist.push_back(k_p_hist_tmp);
};

void OnePointSolverExplicitFPModel::CalTransient(BCGManager &man,int step,real delta_t)
{
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val(num1);
  val.set_zero();

  real n=1.; // initial value for neutron
  val.put_data(0,n);
  for(int i=0;i<num;i++){
    if(lmd[i]!=0.){
      val.put_data(i+1,yield[i]/nu_p/lmd[i]/l*k_p);
    }else{
      val.put_data(i+1,0.);
    };
  };

  Burnup bu;
  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Time[s]  Neutrons   Period\n";
  cout.setf(ios::scientific);
  cout.precision(4);

  real t=0;
  real n_old;

  //tmat=tmat.CalMatrixExponentialByChebyshev14(delta_t);
  //GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);
  //GroupData2D tmp=tmat.CalMatrixExponentialByPade(delta_t);
  //GroupData2D tmp=tmat.CalMatrixExponentialByMMPA32(delta_t);
  tmat.LUDecompositionForMatrixExponentialByMMPA32(delta_t);

  for(int i=0;i<step;i++){
    n_old=n;
    //val=tmat*val;
    //val=tmp*val;
    //val=tmat.CalMatrixExponentialByMMPA32(val,delta_t);
    val=tmat.CalMatrixExponentialByLUDMMPA32(val);
    n=val.get_dat(0);
    real lmd=log(n/n_old)/delta_t;
    real period=1./lmd;
    t+=delta_t;
    cout<<t<<" "<<n<<" "<<period<<"\n";
  };

  cout<<"# Time[s]  Neutrons   Period\n";

};

void OnePointSolverExplicitFPModel::CalTransientFromBurstFission(BCGManager &man,int step,real delta_t)
{
  k_p=0.; // prompt k is set to zero (no fission multiplication)
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val(num1);
  val.set_zero();

  real n=0.; // initial value for neutron
  val.put_data(0,n);
  for(int i=0;i<num;i++){
    val.put_data(i+1,yield[i]);
  };

  Burnup bu;
  cout<<"# Time[s]    t*Neutrons\n";

  // +++ CRAM
  GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);

  real t=0;
  cout.setf(ios::scientific);
  cout.precision(6);
  for(int i=0;i<step;i++){
    val=tmp*val;
    n=val.get_dat(0);
    t+=delta_t;
    real sum=0.;
    for(int j=0;j<num;j++){
      sum+=val.get_dat(j+1)*lmd[j]*dnt[j];
    };
    cout<<t<<" "<<sum*t<<"\n";
  };

  cout<<"# Time[s]    t*Neutrons\n";

};

void OnePointSolverExplicitFPModel::CalTransientNeutronEmission(string srcnuc,BCGManager &man,int step,real delta_t)
{
  // At initial condition, one specific precursor is given to system

  k_p=0.;
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating

  int num=GetTreatedNuclideNumber();
  GroupData1D val(num+1);
  val.set_zero();

  for(int i=0;i<num;i++){
    if(name[i]==srcnuc)val.put_data(i+1,1.);
  };

  Burnup bu;
  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Time[s]  Neutrons\n";

  // +++ CRAM
  GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);

  real t=0;  
  real totn_init;
  real half=0.;
  real halfhalf=0.;
  for(int i=0;i<step;i++){
    real totn=0.;
    for(int j=0;j<num;j++){
      real n=val.get_dat(j+1);
      real lambda=lmd[j];
      real dn=dnt[j];
      totn+=n*lambda*dn;
    };
    if(i==0)totn_init=totn;
    if(half==0.&&totn<=totn_init*0.5)half=t;
    if(halfhalf==0.&&totn<=totn_init*0.25){
      halfhalf=t;
      break;
    };
    val=tmp*val;
    t+=delta_t;
    cout<<t<<" "<<totn<<"\n";
  };

  cout<<"# half     : "<<half<<"\n";
  cout<<"# halfhalf : "<<halfhalf<<"\n";


};


void OnePointSolverExplicitFPModel::CalTransientWithPlotXY(BCGManager &man,int step,real delta_t)
{
  int zmax=60;
  int nmax=100;
  int zmin=30;
  int nmin=50;

  k_p=0.993495; // to make critical state
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val(num1);
  val.set_zero();
  val.put_data(0,1.);
  for(int i=0;i<num;i++){
    if(lmd[i]!=0.){
      val.put_data(i+1,yield[i]/nu_p/lmd[i]/l*k_p);
    }else{
      val.put_data(i+1,0.);
    };
  };


  Burnup bu;
  cout.setf(ios::scientific);
  cout.precision(5);

  GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);
  real t=0;
  real n_old;
  real n=1.;
  for(int i=0;i<step;i++){

    if(i>=50){
    vector< vector<real> > dat(zmax+1,vector<real>(nmax+1,0.));
    real sum=0.;
    for(int j=1;j<num;j++){
      int at=atomn[j-1];
      int ms=massn[j-1];
      real v=val.get_dat(j);
      real lm=lmd[j-1];
      real dn=dnt[j-1];  
      int nn=ms-at;
      real tmp=v*lm*dn;
      if(at<=zmax&&nn<=nmax)dat[at][nn]+=tmp;
      sum+=tmp;
    };

    sum=1./sum;
    for(int ii=nmin;ii<nmax;ii++){
      for(int jj=0;jj<2;jj++){
        for(int j=zmin;j<zmax;j++){
          cout<<ii-0.5+jj<<" "<<j-0.5<<" "<<dat[j][ii]*sum<<"\n";
          cout<<ii-0.5+jj<<" "<<j+0.5<<" "<<dat[j][ii]*sum<<"\n";
        };
        cout<<"\n";
      };
    };
    cout<<"\n";
    };

    n_old=n;
    val=tmp*val;
    n=val.get_dat(0);
    real lmd=log(n/n_old)/delta_t;
    real period=1./lmd;
    t+=delta_t;

    //cout<<t<<" "<<n<<" "<<period<<"\n";
  };

  k_p=0.996;
  tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating
  tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);
  t=0;
  for(int i=0;i<step;i++){

    vector< vector<real> > dat(zmax+1,vector<real>(nmax+1,0.));
    real sum=0.;
    for(int j=1;j<num;j++){
      int at=atomn[j-1];
      int ms=massn[j-1];
      real v=val.get_dat(j);
      real lm=lmd[j-1];
      real dn=dnt[j-1];  
      int nn=ms-at;
      real tmp=v*lm*dn;
      if(at<=zmax&&nn<=nmax)dat[at][nn]+=tmp;
      sum+=tmp;
    };

    sum=1./sum;
    for(int ii=nmin;ii<nmax;ii++){
      for(int jj=0;jj<2;jj++){
        for(int j=zmin;j<zmax;j++){
          cout<<ii-0.5+jj<<" "<<j-0.5<<" "<<dat[j][ii]*sum<<"\n";
          cout<<ii-0.5+jj<<" "<<j+0.5<<" "<<dat[j][ii]*sum<<"\n";
        };
        cout<<"\n";
      };
    };
    cout<<"\n";

    n_old=n;
    val=tmp*val;
    n=val.get_dat(0);
    real lmd=log(n/n_old)/delta_t;
    real period=1./lmd;
    t+=delta_t;

    //cout<<t<<" "<<n<<" "<<period<<"\n";
  };

};

void OnePointSolverExplicitFPModel::CalTransientWithAdjoint(BCGManager &man,int istep,real idelta_t)
{
  cout<<"##########################################################\n";
  cout<<"# Transient calculation\n";

  step=istep;
  delta_t=idelta_t;

  n_fwd.resize(step+1);
  n_adj.resize(step+1);

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData2D tmat=CalDecayMatrixKp(man);

  // +++ Burnup calculation
  Burnup bu;

  // (forward calculation)
  GroupData1D val(num1);
  val.set_zero();
  val.put_data(0,1.);
  for(int i=0;i<num;i++){
    if(lmd[i]!=0.){
      val.put_data(i+1,yield[i]/nu_p/lmd[i]/l*k_p);
    }else{
      val.put_data(i+1,0.);
    };
  };

  // +++ By CRAM
  GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);
  for(int i=0;i<step;i++){
    cout<<"#  (fwd) "<<i<<"/"<<step<<"\n";
    n_fwd[i].copy(val);
    val=tmp*val;
  };

  n_fwd[step].copy(val);

  /*
  real dt=(step-step1)*delta_t;
  real dn=n[step].get_dat(0)/n[step1].get_dat(0);
  real period=dt/(log(dn));
  real end_nuc=n[step].get_dat(0);

  //cout<<k_p<<" "<<period<<"\n";

  cout<<"# k_p    : "<<k_p<<"\n";
  cout<<"# Period : "<<period<<"[s]\n";
  //exit(0);
  cout<<"# Neutron number at final time step : "<<end_nuc<<"\n";
  */

  cout<<"#  (adjoint)\n";
  // (adjoint calculation)
  val.set_zero();
  val.put_data(0,1.);
  GroupData2D tmatt=tmat.GetTransposedMatrix();


  // [CRAM]
  GroupData2D tmpt=tmatt.CalMatrixExponentialByChebyshev14(delta_t);
  for(int i=0;i<step;i++){
    cout<<"#  (adj) "<<i<<"/"<<step<<"\n";
    n_adj[i].copy(val);
    val=tmpt*val;
  };

  n_adj[step].copy(val);

  cout<<"##########################################################\n";
};

void OnePointSolverExplicitFPModel::CalTransientWithAdjointForDNEmission
(BCGManager &man,int istep,real idelta_t)
{
  cout<<"##########################################################\n";
  cout<<"# Transient calculation\n";

  step=istep;
  delta_t=idelta_t;

  n_fwd.resize(step+1);
  n_adj.resize(step+1);

  k_p=0.;
  GroupData2D tmat=CalDecayMatrixKp(man);

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  // (forward calculation)

  GroupData1D val(num1);
  val.set_zero();

  real n=0.; // initial value for neutron
  val.put_data(0,n);
  for(int i=0;i<num;i++){
    val.put_data(i+1,yield[i]);
  };

  Burnup bu;

  // +++ CRAM
  GroupData2D tmp=tmat.CalMatrixExponentialByChebyshev14(delta_t);

  for(int i=0;i<step;i++){
    cout<<"#  (fwd) "<<i<<"/"<<step<<"\n";
    n_fwd[i].copy(val);
    val=tmp*val;
  };
  n_fwd[step].copy(val);

  // (adjoint calculation)

  val.set_zero();
  val.put_data(0,0.); // final condition for number of neutrons
  for(int i=0;i<num;i++){
    val.put_data(i+1,lmd[i]*dnt[i]);
  };

  GroupData2D tmatt=tmat.GetTransposedMatrix();

  // +++ CRAM
  tmp=tmatt.CalMatrixExponentialByChebyshev14(delta_t);
  for(int i=0;i<step;i++){
    cout<<"#  (adj) "<<i<<"/"<<step<<"\n";
    n_adj[i].copy(val);
    val=tmp*val;
  };

  n_adj[step].copy(val);

  cout<<"##########################################################\n";
};

SensitivityData OnePointSolverExplicitFPModel::SensitivityCalculation
(BCGManager &man,int step1,real nu_d,string tagname,int fisid)
{
  man.CalTotalEmittedNeutrons();
  real nu_d_org=man.CalTotalDelayedNeutron(tagname);

  real end_t=step*delta_t;

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  sens1.resize(num);
  sens2.resize(num);
  sens3.resize(num);

  vector<int> tmpid(num);
  for(int i=0;i<num;i++){
    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    tmpid[i]=man.GetNuclideIndex(at,ms,lv);
  };

  real dt=(step-step1)*delta_t;
  real dn=n_fwd[step].get_dat(0)/n_fwd[step1].get_dat(0);
  real period=dt/(log(dn));
  real end_nuc=n_fwd[step].get_dat(0);
  cout<<"##########################################################\n";
  cout<<"# Sensitivity calculation for reactor period\n";
  cout<<"#\n";
  cout<<"#   k_p    : "<<k_p<<"\n";
  cout<<"#   Period : "<<period<<"[s]\n";
  cout<<"#   Neutron number at final time step : "<<end_nuc<<"\n";
  cout<<"##########################################################\n";

  SensitivityData sns;
  sns.PutName("one_point","period","unknown");
  sns.PutValue(period);
  sns.PutGroup(0);

  // +++ Anscestor contribution calculation
  /*
  int step_i=1500;
  real sum=n_adj[step-step_i].get_dat(0)*n[step_i].get_dat(0)/end_nuc;
  for(int i=0;i<num;i++){
    real ww=n_adj[step-step_i].get_dat(i+1)*n[step_i].get_dat(i+1);
    ww/=end_nuc;
    sum+=ww;
    if(ww>0.01)cout<<name[i]<<" "<<ww<<"\n";
  };
  cout<<sum<<"\n";
  */

  // +++ neutron number to yield sensitivity +++
  /*
  for(int i=0;i<num;i++){
    real tmp=0.;
    for(int j=0;j<step;j++){
      real avg=(n[j].get_dat(0)*n_adj[step-j].get_dat(i+1)+
                n[j+1].get_dat(0)*n_adj[step-j-1].get_dat(i+1))*0.5;
      tmp+=delta_t*avg/(nu*ge);
    };
    cout<<name[i]<<" "<<tmp*yld[i]/n[step].get_dat(0)<<"\n";
  };
  */

  real aa=(period*period)/dt/n_fwd[step1].get_dat(0);
  real bb=-(period*period)/dt/n_fwd[step].get_dat(0);

  // +++ period to yield sensitivity +++
  vector<real> sens1_un(num);
  for(int i=0;i<num;i++){
    real tmp=0.;
    for(int j=0;j<step;j++){
      real avg=(n_fwd[j].get_dat(0)*n_adj[step-j].get_dat(i+1)+
                n_fwd[j+1].get_dat(0)*n_adj[step-j-1].get_dat(i+1))*0.5;
      tmp+=delta_t*avg;
    };
    tmp*=k_p/(l*nu_p);
    real tmp2=0.;
    for(int j=0;j<step1;j++){
      real avg=(n_fwd[j].get_dat(0)*n_adj[step1-j].get_dat(i+1)+
                n_fwd[j+1].get_dat(0)*n_adj[step1-j-1].get_dat(i+1))*0.5;
      tmp2+=delta_t*avg;
    };
    tmp2*=k_p/(l*nu_p);
    sens1_un[i]=(aa*tmp2+bb*tmp)/period*yield[i];
  };

  // +++ Normalization condition consideration for yield sensitivity
  for(int i=0;i<num;i++){
    sens1[i]=sens1_un[i];
    real coef=dnt_dtr[i]/nu_d*yield[i];
    for(int j=0;j<num;j++){
       sens1[i]-=coef*sens1_un[j];
    };
  };

  for(int i=0;i<num;i++){
    sns.PutSensitivity0D(matid[i],18000000+fisid,sens1[i]);
  };

  // +++ decay constant sensitivity +++
  for(int i=0;i<num;i++){

    real n_adj_dtdl_l=CalDecayConstantSensitivityPart(man,n_adj[step],i);
    real tmp=0.;
    for(int j=0;j<step;j++){
      real n_adj_dtdl_r=CalDecayConstantSensitivityPart(man,n_adj[step-1-j],i);
      real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
      tmp+=delta_t*avg;
      n_adj_dtdl_l=n_adj_dtdl_r;
    };

    n_adj_dtdl_l=CalDecayConstantSensitivityPart(man,n_adj[step1],i);
    real tmp2=0.;
    for(int j=0;j<step1;j++){
      real n_adj_dtdl_r=CalDecayConstantSensitivityPart(man,n_adj[step1-1-j],i);
      real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
      tmp2+=delta_t*avg;
      n_adj_dtdl_l=n_adj_dtdl_r;
    };

    sens2[i]=(aa*tmp2+bb*tmp)/period*lmd[i];
  };

  for(int i=0;i<num;i++){
    //sns.PutSensitivity0D(matid[i],8887,sens2[i]);
    sns.PutSensitivity0D(matid[i],8888,-sens2[i]);
  };

  // +++ branching ratio sensitivity
  vector< vector<real> > sens3_un(num);
  for(int i=0;i<num;i++){
    int id=tmpid[i];
    int ch=man.GetNuclide(id).GetChannel();
    sens3[i].resize(ch);
    sens3_un[i].resize(ch);
    for(int k=0;k<ch;k++){

      real n_adj_dtdl_l=CalDecayConstantSensitivityPart2(man,n_adj[step],i,id,k);
      real tmp=0.;
      for(int j=0;j<step;j++){
        real n_adj_dtdl_r=CalDecayConstantSensitivityPart2(man,n_adj[step-1-j],i,id,k);
        real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                  n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
        tmp+=delta_t*avg;
        n_adj_dtdl_l=n_adj_dtdl_r;
      };

      n_adj_dtdl_l=CalDecayConstantSensitivityPart2(man,n_adj[step1],i,id,k);
      real tmp2=0.;
      for(int j=0;j<step1;j++){
        real n_adj_dtdl_r=CalDecayConstantSensitivityPart2(man,n_adj[step1-1-j],i,id,k);
        real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                  n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
        tmp2+=delta_t*avg;
        n_adj_dtdl_l=n_adj_dtdl_r;
      };

      sens3_un[i][k]=(aa*tmp2+bb*tmp)/period*man.GetNuclide(id).GetBr(k);
    };
  };

  // +++ normalization
  for(int i=0;i<num;i++){
    int id=tmpid[i];
    int ch=man.GetNuclide(id).GetChannel();
    vector<real> br_org(ch);
    for(int j=0;j<ch;j++){
      br_org[j]=man.GetNuclide(id).GetBr(j);
    };
    for(int j=0;j<ch;j++){
      real br=br_org[j];
      sens3[i][j]=sens3_un[i][j];
      for(int k=0;k<ch;k++){
	sens3[i][j]-=sens3_un[i][k]*br;
      };
      // (yield normalization)
      real delta=br*0.01;
      real factor=1./(1.+delta);
      for(int k=0;k<ch;k++){
	if(j==k){
          man.GetNuclide(id).PutBr(k,(br_org[k]+delta)*factor);
	}else{
	  man.GetNuclide(id).PutBr(k,br_org[k]*factor);
	};
      };
      man.CalTotalEmittedNeutrons();
      real nu_d_delta=man.CalTotalDelayedNeutron(tagname);
      real coef=(nu_d_delta-nu_d_org)/nu_d_org*100.;
      for(int k=0;k<ch;k++){
	man.GetNuclide(id).PutBr(k,br_org[k]);
      };
      for(int k=0;k<num;k++){
        sens3[i][j]-=sens1_un[k]*coef;
      };
    };
  };
  for(int i=0;i<num;i++){
    int ch=man.GetNuclide(tmpid[i]).GetChannel();
    for(int j=0;j<ch;j++){
      sns.PutSensitivity0D(matid[i],88880+j,sens3[i][j]);
    };
  };

  return sns;
};

SensitivityData OnePointSolverExplicitFPModel::SensitivityCalculationForDNEmission
(BCGManager &man,string tagname,int fisid)
{
  man.CalTotalEmittedNeutrons();
  real nu_d=man.CalTotalDelayedNeutron(tagname);

  real end_t=step*delta_t;

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  sens1.resize(num);
  sens2.resize(num);
  sens3.resize(num);

  vector<int> tmpid(num);
  for(int i=0;i<num;i++){
    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    tmpid[i]=man.GetNuclideIndex(at,ms,lv);
  };

  real parameter=0.;
  for(int i=0;i<num;i++){
    parameter+=n_fwd[step].get_dat(i+1)*lmd[i]*dnt[i];
  };

  cout<<"##########################################################\n";
  cout<<"# Sensitivity calculation for DN emission\n#\n";
  cout<<"#   DN emission rate at t= "<<end_t<<" : "<<parameter<<"\n";
  cout<<"##########################################################\n";

  SensitivityData sns;
  sns.PutName("one_point","dn_emission_rate","unknown");
  sns.PutValue(parameter);
  sns.PutGroup(0);

  // +++ Sensitivity to yield +++
  vector<real> sens1_un(num);
  for(int i=0;i<num;i++){
    sens1_un[i]=(n_adj[step].get_dat(i+1)*n_fwd[0].get_dat(i+1))/parameter; 
  };

  // +++ Normalization condition consideration for yield sensitivity
  for(int i=0;i<num;i++){
    sens1[i]=sens1_un[i];
    /*
    real coef=dnt_dtr[i]/nu_d*yield[i];
    for(int j=0;j<num;j++){
      sens1[i]-=coef*sens1_un[j];
    };
    */
  };

  for(int i=0;i<num;i++){
    sns.PutSensitivity0D(matid[i],18000000+fisid,sens1[i]);
  };

  // +++ Decay constant sensitivity +++
  for(int i=0;i<num;i++){
    real sum=0.;
    real n_adj_dtdl_l=CalDecayConstantSensitivityPart3(man,n_adj[step],i);
    for(int j=0;j<step;j++){
      real n_adj_dtdl_r=CalDecayConstantSensitivityPart3(man,n_adj[step-1-j],i);
      real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
      sum+=delta_t*avg;
      n_adj_dtdl_l=n_adj_dtdl_r;
    };
    sens2[i]=(sum+dnt[i]*n_fwd[step].get_dat(i+1))*lmd[i]/parameter; 
    // (direct term is also considered.)
  };

  for(int i=0;i<num;i++){
    sns.PutSensitivity0D(matid[i],8888,-sens2[i]);
  };

  // +++ branching ratio sensitivity
  vector< vector<real> > sens3_un(num);
  for(int i=0;i<num;i++){
    int id=tmpid[i];
    int ch=man.GetNuclide(id).GetChannel();
    sens3[i].resize(ch);
    sens3_un[i].resize(ch);
    for(int k=0;k<ch;k++){
      real br=man.GetNuclide(id).GetBr(k);
      real n_adj_dtdl_l=CalDecayConstantSensitivityPart4(man,n_adj[step],i,id,k);
      real sum=0.;
      for(int j=0;j<step;j++){
        real n_adj_dtdl_r=CalDecayConstantSensitivityPart4(man,n_adj[step-1-j],i,id,k);
        real avg=(n_adj_dtdl_l*n_fwd[j].get_dat(i+1)+
                  n_adj_dtdl_r*n_fwd[j+1].get_dat(i+1))*0.5;
        sum+=delta_t*avg;
        n_adj_dtdl_l=n_adj_dtdl_r;
      };
      real dnt_ch=man.GetNuclide(id).GetEmittedNeutron(k);
      sens3_un[i][k]=(sum+dnt_ch*n_fwd[step].get_dat(i+1)*br*lmd[i])/parameter;
    };
  };

  // +++ BR normalization
  for(int i=0;i<num;i++){
    int id=tmpid[i];
    int ch=man.GetNuclide(id).GetChannel();
    vector<real> br_org(ch);
    for(int j=0;j<ch;j++){
      br_org[j]=man.GetNuclide(id).GetBr(j);
    };
    for(int j=0;j<ch;j++){
      real br=br_org[j];
      sens3[i][j]=sens3_un[i][j];
      for(int k=0;k<ch;k++){
	sens3[i][j]-=sens3_un[i][k]*br;
      };
      /*
      // (yield normalization)
      real delta=br*0.01;
      real factor=1./(1.+delta);
      for(int k=0;k<ch;k++){
	if(j==k){
          man.GetNuclide(id).PutBr(k,(br_org[k]+delta)*factor);
	}else{
	  man.GetNuclide(id).PutBr(k,br_org[k]*factor);
	};
      };
      man.CalTotalEmittedNeutrons();
      real nu_d_delta=man.CalTotalDelayedNeutron(tagname);
      real coef=(nu_d_delta-nu_d_org)/nu_d_org*100.;
      for(int k=0;k<ch;k++){
	man.GetNuclide(id).PutBr(k,br_org[k]);
      };
      for(int k=0;k<num;k++){
        sens3[i][j]-=sens1_un[k]*coef;
      };
      */
    };
  };

  for(int i=0;i<num;i++){
    int ch=man.GetNuclide(tmpid[i]).GetChannel();
    if(ch>1){
      for(int j=0;j<ch;j++){
        sns.PutSensitivity0D(matid[i],88880+j,sens3[i][j]);
      };
    };
  };

  return sns;
};

void OnePointSolverExplicitFPModel::BranchingRatioCheck(BCGManager &man)
{
  int num=GetTreatedNuclideNumber();
  for(int i=0;i<num;i++){

    int at=atomn[i];
    int ms=massn[i];
    int lv=level[i];
    int id=man.GetNuclideIndex(at,ms,lv);
    int ch=man.GetNuclide(id).GetChannel();
    real sum_n=0.;
    real sum_br=0.;
    for(int k=0;k<ch;k++){
      sum_n+=man.GetNuclide(id).GetBr(k)*man.GetNuclide(id).GetEmittedNeutron(k);
      sum_br+=man.GetNuclide(id).GetBr(k);
    };
    if(sum_n!=dnt[i]){
      cout<<"# OnePointSolverExplicitFPModel::BranchingRatioCheck.\n";
      cout<<"# Total emitted neutron is inconsistent.\n";
      cout<<"# at/ms/lv = "<<at<<"/"<<ms<<"/"<<lv<<"\n";
      cout<<"# Summation  : "<<sum_n<<"\n";
      cout<<"# dnt        : "<<dnt[i]<<"\n"; 
      cout<<"# sum of branching ratio : "<<sum_br<<"\n";
      exit(0);
    };
  };
};

real OnePointSolverExplicitFPModel::CalDecayConstantSensitivityPart
(BCGManager &man,GroupData1D &n_adj,int ii)
{
  int num=GetTreatedNuclideNumber();

  real ret=n_adj.get_dat(0)*dnt[ii]-n_adj.get_dat(ii+1);

  int at=atomn[ii];
  int ms=massn[ii];
  int lv=level[ii];
  int id=man.GetNuclideIndex(at,ms,lv);
  int ch=man.GetNuclide(id).GetChannel();
  for(int j=0;j<ch;j++){
    int at2=man.GetNuclide(id).GetAtomicNumberNext(j);
    int ms2=man.GetNuclide(id).GetMassNumberNext(j);
    int lv2=man.GetNuclide(id).GetExLevelNext(j);
    real br=man.GetNuclide(id).GetBr(j);
    int id2=man.GetNuclideIndex(at2,ms2,lv2);
    if(id2!=-1){
      int idi=pos[id2];
      if(chain_effect)ret+=n_adj.get_dat(idi+1)*br;
    };
  };

  return ret;
};

real OnePointSolverExplicitFPModel::CalDecayConstantSensitivityPart2
(BCGManager &man,GroupData1D &n_adj,int ii,int id,int jj)
{
  int at2=man.GetNuclide(id).GetAtomicNumberNext(jj);
  int ms2=man.GetNuclide(id).GetMassNumberNext(jj);
  int lv2=man.GetNuclide(id).GetExLevelNext(jj);
  real br=man.GetNuclide(id).GetBr(jj);

  real ret=n_adj.get_dat(0)*(lmd[ii]*man.GetNuclide(id).GetEmittedNeutron(jj));

  int id2=man.GetNuclideIndex(at2,ms2,lv2);
  if(id2!=-1){
    int idi=pos[id2];
    if(chain_effect)ret+=n_adj.get_dat(idi+1)*lmd[ii];
  };

  return ret;
};

real OnePointSolverExplicitFPModel::CalDecayConstantSensitivityPart3
(BCGManager &man,GroupData1D &n_adj,int ii)
{
  real ret=-n_adj.get_dat(ii+1);

  int at=atomn[ii];
  int ms=massn[ii];
  int lv=level[ii];
  int id=man.GetNuclideIndex(at,ms,lv);
  int ch=man.GetNuclide(id).GetChannel();
  for(int j=0;j<ch;j++){
    int at2=man.GetNuclide(id).GetAtomicNumberNext(j);
    int ms2=man.GetNuclide(id).GetMassNumberNext(j);
    int lv2=man.GetNuclide(id).GetExLevelNext(j);
    real br=man.GetNuclide(id).GetBr(j);
    int id2=man.GetNuclideIndex(at2,ms2,lv2);
    if(id2!=-1){
      int idi=pos[id2];
      ret+=n_adj.get_dat(idi+1)*br;
    };
  };

  return ret;
};


real OnePointSolverExplicitFPModel::CalDecayConstantSensitivityPart4
(BCGManager &man,GroupData1D &n_adj,int ii,int id,int jj)
{
  int at2=man.GetNuclide(id).GetAtomicNumberNext(jj);
  int ms2=man.GetNuclide(id).GetMassNumberNext(jj);
  int lv2=man.GetNuclide(id).GetExLevelNext(jj);
  real br=man.GetNuclide(id).GetBr(jj);
  int id2=man.GetNuclideIndex(at2,ms2,lv2);
  if(id2!=-1){
    int idi=pos[id2];
    return n_adj.get_dat(idi+1)*br*lmd[ii];
  };
  return 0.;
};


void OnePointSolverExplicitFPModel::DensityDataPrinting(int nucp,string *nucinp,string type)
{
  real end_t=step*delta_t;
  int num=GetTreatedNuclideNumber();
  int num1=num+1;
  real end_nuc=n_fwd[step].get_dat(0);

  vector<int> nucid(nucp);
  for(int i=0;i<nucp;i++){
    for(int j=0;j<num;j++){
      if(nucinp[i]==name[j])nucid[i]=j;
    };
  };

  if(type=="fwd"){ 
    cout<<"# Forward number density\n";
  }else if(type=="adj"){
    cout<<"# Adjoint number density\n";
  }else if(type=="fwdadj"||type=="adjfwd"){
    cout<<"# Contribution function (fwd*adj/n_e)\n";
  }else{
    cout<<"# Error in OnePointSolverExplcitFPModel::DensityDataPrinting.\n";
    exit(0);
  };

  cout<<"# t[s]     t-t_end[s]  ";
  for(int i=0;i<nucp;i++){
    cout<<nucinp[i]<<"     ";
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(3);
  real t=0.;
  for(int i=0;i<step+1;i++){
    cout<<t<<" "<<t-end_t<<" ";
    real sum=0.;
    for(int j=0;j<nucp;j++){
      real ww=0.;
      if(type=="fwd"){
	ww=n_fwd[i].get_dat(nucid[j]+1);
      }else if(type=="adj"){
	ww=n_adj[step-i].get_dat(nucid[j]+1);
      }else{
        ww=n_adj[step-i].get_dat(nucid[j]+1)*n_fwd[i].get_dat(nucid[j]+1);
        ww/=end_nuc;
      };
      sum+=ww;
      cout<<ww<<" ";
      //cout<<n[i].get_dat(nucid[j]+1)<<" ";
      //cout<<n[i].get_dat(nucid[j]+1)*lmd[nucid[j]]*dnt[nucid[j]]<<" ";
    };
    cout<<sum<<"\n";
    //cout<<"\n";
    t+=delta_t;
  };
};

void OnePointSolverExplicitFPModel::DensityDataPrintingForXYPlot()
{
  int zmax=60;
  int nmax=100;
  int zmin=30;
  int nmin=50;

  real end_t=step*delta_t;
  int num=GetTreatedNuclideNumber();
  int num1=num+1;
  real end_nuc=n_fwd[step].get_dat(0);

  cout.setf(ios::scientific);
  cout.precision(5);
  real t=0.;
  for(int i=step;i>=0;i--){

    vector< vector<real> > dat(zmax+1,vector<real>(nmax+1,0.));
    for(int j=1;j<num;j++){
      real ww=n_adj[step-i].get_dat(j)*n_fwd[i].get_dat(j);
      ww/=end_nuc;
      int at=atomn[j-1];
      int ms=massn[j-1];
      int nn=ms-at;
      if(at<=zmax&&nn<=nmax)dat[at][nn]+=ww;
    };
    for(int ii=nmin;ii<nmax;ii++){
      for(int jj=0;jj<2;jj++){
        for(int j=zmin;j<zmax;j++){
          cout<<ii-0.5+jj<<" "<<j-0.5<<" "<<dat[j][ii]<<"\n";
          cout<<ii-0.5+jj<<" "<<j+0.5<<" "<<dat[j][ii]<<"\n";
        };
        cout<<"\n";
      };
    };
    cout<<"\n";

  };
};


void OnePointSolverExplicitFPModel::UncertaintyCalculation(BCGManager &man)
{
  real var1=0.;
  real var2=0.;
  real var3=0.;

  int sz=man.GetNuclideNumber();

  // +++ decay constant uncertainty
  cout<<"# Decay constant uncertainty\n";
  int ii=0;
  for(int i=0;i<sz;i++){
    if(man.GetNuclide(i).Flag()){
      real hl=man.GetNuclide(i).GetHalflife();
      real dhl=man.GetNuclide(i).GetDeltaHalflife();
      if(hl!=0.){
        real dlambda=dhl/hl;
        real tmp=dlambda*dlambda*sens2[ii]*sens2[ii];
        var1+=tmp;
	//if(tmp>0.001)cout<<"#  "<<at<<" "<<ms<<" "<<lv<<" : "<<sqrt(tmp)<<" "<<dlambda<<"\n";
	if(tmp>1e-6)cout<<"#  "<<name[ii]<<" : "<<sqrt(tmp)<<" "<<dlambda<<"\n";
      };
      ii++;
    };
  };

  cout<<"#\n#\n# Fission yield uncertainty\n";
  ii=0;
  for(int i=0;i<sz;i++){
    if(man.GetNuclide(i).Flag()){
      if(yield[ii]!=0.){
        real dy=delta_yield[ii]/yield[ii];
        real tmp=dy*dy*sens1[ii]*sens1[ii];
        var2+=tmp;
	if(tmp>1e-6)cout<<"#  "<<name[ii]<<" : "<<sqrt(tmp)<<" "<<dy<<"\n";
      };
      ii++;
    };
  };

  cout<<"#\n#\n# Branching ratio uncertainty\n";
  ii=0;
  for(int i=0;i<sz;i++){
    if(man.GetNuclide(i).Flag()){
      int ch=man.GetNuclide(i).GetChannel();
      real tmpsum=0.;
      for(int j=0;j<ch;j++){
	if(man.GetNuclide(i).GetBr(j)!=0.){
	  real dbr=man.GetNuclide(i).GetDeltaBr(j)/man.GetNuclide(i).GetBr(j);
	  real tmp=dbr*dbr*sens3[ii][j]*sens3[ii][j];
	  var3+=tmp;
	  tmpsum+=tmp;
	};
      };
      if(tmpsum>1e-6)cout<<"#  "<<name[ii]<<" : "<<sqrt(tmpsum)<<"\n";
      ii++;
    };
  };

  real var=var1+var2+var3;
  cout<<"#\n#\n SUMMARY\n#\n";
  cout<<"# Decay constant uncertainty  : "<<sqrt(var1)<<"\n";
  cout<<"# Fission yield uncertainty   : "<<sqrt(var2)<<"\n";
  cout<<"# Branching ratio uncertainty : "<<sqrt(var3)<<"\n";
  cout<<"# Total uncertainty           : "<<sqrt(var)<<"\n";
};



void OnePointSolverExplicitFPModel::MakeStationaryState(BCGManager &man)
{
  real delta_t=1.; // [sec]
  int stepmax=100000; 
  
  GroupData2D tmat=CalDecayMatrixKp(man);  // +++ Decay matrix creating

  int num=GetTreatedNuclideNumber();
  int num1=num+1;

  GroupData1D val(num1);
  val.set_zero();

  real n=1.; // initial value for neutron
  val.put_data(0,n);
  for(int i=0;i<num;i++){
    if(lmd[i]!=0.){
      val.put_data(i+1,yield[i]/nu_p/lmd[i]/l*k_p);
    }else{
      val.put_data(i+1,0.);
    };
  };

  Burnup bu;

  real t=0;
  real n_old;

  tmat.LUDecompositionForMatrixExponentialByMMPA32(delta_t);

  for(int i=0;i<stepmax;i++){
    n_old=n;
    val=tmat.CalMatrixExponentialByLUDMMPA32(val);
    n=val.get_dat(0);
    real lmd=log(n/n_old)/delta_t;
    real period=1./lmd;
    t+=delta_t;
    if(fabs(period)>1e5)break;
    if(i==stepmax-1){
      cout<<"# Error in OnePointSolverExplicitFPModel::MakeStationaryState.\n";
      cout<<"# Iteration reaches maximum.\n";
      exit(0);
    };
  };

};

