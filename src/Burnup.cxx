#include <cstdlib>
#include "Burnup.h"

AtomicMassData::AtomicMassData()
{
  // Atomic weight (atomic mass unit)
  real anmu=1.008665; // From SRAC burnup-chain data
  // (from neutron mass unit to atomic mass unit)

  // Th
  atomic_weight[902270]=225.077*anmu;
  atomic_weight[902280]=226.070*anmu;
  atomic_weight[902290]=227.064*anmu;
  atomic_weight[902300]=228.057*anmu;
  atomic_weight[902310]=229.052*anmu;
  atomic_weight[902320]=230.045*anmu;
  // Pa
  atomic_weight[912310]=229.051*anmu;
  atomic_weight[912330]=229.051*anmu;
  // U
  atomic_weight[922320]=230.044*anmu;
  atomic_weight[922330]=231.038*anmu;
  atomic_weight[922340]=232.03*anmu;
  atomic_weight[922350]=233.025*anmu;
  atomic_weight[922360]=234.018*anmu;
  atomic_weight[922370]=235.012*anmu;
  atomic_weight[922380]=236.006*anmu;
  // Np
  atomic_weight[932360]=233.973*anmu;
  atomic_weight[932370]=235.012*anmu;
  atomic_weight[932390]=236.999*anmu;
  // Pu
  atomic_weight[942360]=234.018*anmu;
  atomic_weight[942380]=236.005*anmu;
  atomic_weight[942390]=236.999*anmu;
  atomic_weight[942400]=237.992*anmu;
  atomic_weight[942410]=238.986*anmu;
  atomic_weight[942420]=239.979*anmu;
  // Am
  atomic_weight[952410]=238.986*anmu;
  atomic_weight[952420]=239.98*anmu;
  atomic_weight[952421]=239.98*anmu;
  atomic_weight[952430]=240.973*anmu;
  // Cm
  atomic_weight[962420]=239.980*anmu;
  atomic_weight[962430]=240.972*anmu;
  atomic_weight[962440]=241.966*anmu;
  atomic_weight[962450]=242.961*anmu;
  atomic_weight[962460]=243.953*anmu;
  atomic_weight[962470]=244.948*anmu;
  atomic_weight[962480]=245.941*anmu;
  atomic_weight[962490]=246.936*anmu;
  atomic_weight[962500]=247.930*anmu;
  // Bk
  atomic_weight[972490]=246.935*anmu;
  // PseudoFP
  atomic_weight[9950000]=atomic_weight[922350];
  atomic_weight[9980000]=atomic_weight[922380];
  atomic_weight[9990000]=atomic_weight[942390];
  atomic_weight[9910000]=atomic_weight[942410];

  atomic_weight[40090]=9.012182;
  atomic_weight[70150]=14.8713*anmu;
  atomic_weight[80160]=15.9949146;
  atomic_weight[110230]=22.98977;
  //atomic_weight[240000]=51.9661;
  atomic_weight[240000]=51.99561;
  //atomic_weight[250550]=55.938;
  atomic_weight[250550]=54.938;
  //atomic_weight[260000]=55.845;
  atomic_weight[260000]=55.847;
  atomic_weight[270590]=58.4269;
  atomic_weight[270600]=50.41896;
  //atomic_weight[280000]=58.6934;
  atomic_weight[280000]=58.69;
  //atomic_weight[420000]=95.96;
  atomic_weight[420000]=95.94;
  atomic_weight[740000]=183.84;
};

void AtomicMassData::ReadFile(string cbglibdir,string fname)
{
  string mdir=cbglibdir;

  mdir=mdir+"CBGLIB_BURN/atomic_mass/";

  mdir.append(fname);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"File name is "<<mdir<<"\n";
    exit(0);
  };

  int nuc;
  real par;

  bool ed=false;
  while(!ed){
    fin>>nuc;
    if(nuc==-1){
      ed=true;
    }else{
      if(nuc<10000)nuc=midt.GetMATIDFromENDFID(nuc);
      fin>>par;
      atomic_weight[nuc]=par;
    };    
  };

  fin.close();
};

void AtomicMassData::ShowSelf()
{
  map<int,real>::iterator it;
  for(it=atomic_weight.begin();it!=atomic_weight.end();it++){
    int i=it->first;
    cout<<i<<" "<<atomic_weight[i]<<"\n";
  };
};

// +++++++++++++++++++++++++++++++++
//    Reaction energy data
// +++++++++++++++++++++++++++++++++

ReactionEnergyData::ReactionEnergyData()
{
};

void ReactionEnergyData::PutGammaEnergySeparatedData()
{
  // Fission : Beta energy and FP energy are only considered.
  // Capture : Decay heat of daughter nuclide is considered.
  // Other heat source is considered through gamma-ray.
  //
  // Unit [J/fission]

  fis.clear();
  cap.clear();

  // +++ Fission
  //
  // Beta-ray energy and FP energy are only considered.

  // U
  fis[922340]=3.038e-11;
  fis[922350]=2.8138e-11; //
  fis[922360]=3.075e-11;
  fis[922380]=2.8490e-11; //
  // Np
  fis[932370]=3.106e-11;
  fis[932390]=3.137e-11;
  // Pu
  fis[942380]=3.162e-11;
  fis[942390]=2.9014e-11; //
  fis[942400]=2.8862e-11; //
  fis[942410]=2.9150e-11; //
  fis[942420]=3.190e-11;
  // Am
  fis[952410]=3.236e-11;
  fis[ID("Am242m")]=3.227e-11;
  fis[ID("Am243")]=3.220e-11;
  // Cm
  fis[ID("Cm242")]=3.267e-11;
  fis[ID("Cm243")]=3.271e-11;
  fis[ID("Cm244")]=3.253e-11;
  fis[ID("Cm245")]=3.279e-11;
  fis[ID("Cm246")]=3.283e-11;
  fis[ID("Cm247")]=3.287e-11;

  // +++ Capture
  // 
  // Decay heat of daughter nuclide is considered here.

  // U
  cap[922360]=0.3192*mev_to_j; // U-237
  cap[922380]=0.8619*mev_to_j; // U-239 (0.46MeV)+Np-239(0.41MeV)
  // Np
  cap[932370]=0.8080*mev_to_j; // Np-238 
  cap[932390]=1.7880*mev_to_j; // Np-240
  // Pu
  cap[942420]=0.1947*mev_to_j; // Pu-243
  // Am
  cap[952410]=4.5250*mev_to_j; // ?
  cap[ID("Am243")]=0.8840*mev_to_j; // ?

  //cap[525]=2.78868*mev_to_j;
  cap[50100]=2.36*mev_to_j; // gamma energy is removed
  cap[110230]=4.68*mev_to_j; 
  cap[250550]=2.52*mev_to_j;
  //
  cap[260580]=1.31*mev_to_j;
  //
  cap[280640]=1.18*mev_to_j;
};


void ReactionEnergyData::ReadFile(string mdir,string filename)
{
  fis.clear();
  cap.clear();

  ifstream fin;
  string fname=mdir+filename;
  fin.open(fname.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<fname<<"\n";
    exit(0);
  };

  int unit;
  fin>>unit; // 0:J, 1:MeV

  // fission
  int mat=0;
  while(mat!=-1){
    fin>>mat;
    if(mat!=-1){
      real tmp;
      fin>>tmp;
      if(unit==1)tmp*=mev_to_j;
      fis[mat]=tmp;
    };
  };

  // capture
  mat=0;
  while(mat!=-1){
    fin>>mat;
    if(mat!=-1){
      real tmp;
      fin>>tmp;
      if(unit==1)tmp*=mev_to_j;
      cap[mat]=tmp;
    };
  };

};

// ++++++++++++++++++++++++++++++++++++++++++++++
//   Burnup class
// ++++++++++++++++++++++++++++++++++++++++++++++

Burnup::Burnup()
{
  order_krylov=20;
};

void Burnup::PutNucnum(int i)
{
  nucnum=i;
  nucid.resize(i);
  sgf.resize(i);
  sgc.resize(i);
  sgn2n.resize(i);
  dens.resize(i);
  trmat_flxdep.put_yx(i,i);
  trmat_flxindep.put_yx(i,i);
};

void Burnup::PutNuclideData(int i,int id,real den,real sf,real sc,real sn2n)
{
  if(i>=nucnum){
    cout<<"Error in Burnup::PutNuclideData.\n";
    exit(0);
  };
  nucid[i]=id;
  sgf[i]=sf;
  sgc[i]=sc;
  sgn2n[i]=sn2n;
  dens[i]=den;
};

int Burnup::SearchNuclide(int id)
{
  for(int i=0;i<nucnum;i++){
    if(nucid[i]==id)return i;
  };
  return -1;
};

GroupData2D Burnup::CalTransitionMatrix(real flx,bool decay)
{
  flx*=1e-24; // since unit of XS is barn
  GroupData2D tmat(nucnum,nucnum);
  tmat.set_zero();

  if(decay){
    tmat=trmat_flxdep*flx+trmat_flxindep;
  }else{
    tmat=trmat_flxdep*flx;
  };
  return tmat;
};

GroupData1D Burnup::CalMatrixExponentialByKrylov(GroupData2D &tmat,GroupData1D &inp,real delt)
{
  return tmat.CalMatrixExponentialByKrylov(inp,delt,order_krylov);
};


void Burnup::BurnupCalculationNew(real flx,real delt,bool print)
{
  // Matrix exponential is calculated by the Taylor expansion method.
  // The efficient algorithm of SWAT is adopted.
  // Only the non-zero elements are considered in the matrix calculation
  //
  // (See the notebook : 2008/10/10)

  int order_limit=50;

  //GroupData2D tmat=CalTransitionMatrix(flx);
  flx*=1e-24;
  GroupData2D tmat=trmat_flxdep*flx+trmat_flxindep;

  vector<real> non0num(nucnum);
  vector< vector<int> > xpos(nucnum);
  vector< vector<real> > val(nucnum);

  for(int i=0;i<nucnum;i++){
    int tmp=0;
    for(int j=0;j<nucnum;j++){
      real tval=tmat.get_dat(i,j);
      if(tval!=0.){
	xpos[i].push_back(j);
	val[i].push_back(tval);
	tmp++;
      };
    };
    non0num[i]=tmp;
  };

  int order=tmat.CalOrderForMatrixExponential(delt);
  if(order>order_limit){
    cout<<"# Error in Burnup::BurnupCalulationNew.\n";
    cout<<"# You have to set smaller delta_t.\n";
    cout<<"# Necessary order is "<<order<<"\n";
    exit(0);
  };

  vector<real> den_new(nucnum);
  vector<real> den_org(nucnum);
  vector<real> den_tmp(nucnum);
  for(int n=0;n<nucnum;n++){
    den_org[n]=dens[n];
    den_new[n]=dens[n];
  };

  for(int i=0;i<=order;i++){
    real fact=delt/(i+1.);
    for(int n=0;n<nucnum;n++){
      real tmp=0.;
      for(int x=0;x<non0num[n];x++){
        int pos=xpos[n][x];
	tmp+=tmat.get_dat(n,pos)*dens[pos];
      };
      den_tmp[n]=tmp*fact;
      den_new[n]+=tmp*fact;
    };
    for(int n=0;n<nucnum;n++){
      dens[n]=den_tmp[n];
    };
  };

    for(int n=0;n<nucnum;n++){
      dens[n]=den_new[n];
    };

  if(print){
    for(int i=0;i<nucnum;i++){
      cout<<nucid[i]<<" : "<<den_org[i]<<" -> "<<dens[i]<<"\n";
    };
  };
};

void Burnup::BurnupCalculation(real flx,real delt,string opt,bool print)
{
  // Matrix exponential is calculated at the first step with the arbitrary scheme,
  // then the number density vector of the next time step is calculated
  // by multiplying the matrix exponential by the number density vector at the present time step.

  GroupData2D tmat=trmat_flxindep+trmat_flxdep*flx*1e-24;

  /*
  // burnup matrix printing
  //
  cout.setf(ios::scientific);
  cout.precision(5);
  int sz=tmat.get_y();
  for(int i=0;i<sz;i++){
    for(int j=0;j<sz;j++){
      real dat=tmat.get_dat(i,j);
      if(dat==0.){
	cout<<"0 ";
      }else{
        cout<<tmat.get_dat(i,j)<<" ";
      };
    };
    cout<<"\n";
  };
  exit(0);
  */

  GroupData2D last;
  if(opt=="taylor"||opt=="Taylor"){
    last=tmat.CalMatrixExponentialByTaylor(delt);
  }else if(opt=="pade"||opt=="Pade"){
    //last=CalMatrixExponentialByPade(tmat,delt);
    last=tmat.CalMatrixExponentialByPade(delt);
  }else if(opt=="chebyshev"||opt=="Chebyshev"){
    //last=tmat.CalMatrixExponentialByChebyshev(delt,16);
    last=tmat.CalMatrixExponentialByChebyshev14(delt);
  }else if(opt=="mmpa"||opt=="MMPA"){
    last=tmat.CalMatrixExponentialByMMPA32(delt);
  }else{
    cout<<"# Error in Burnup::BurnupCalculation.\n";
    cout<<"# Inappropriate option : "<<opt<<"\n";
    exit(0);
  };

  GroupData1D den_org(nucnum);
  GroupData1D den_new(nucnum);
  for(int i=0;i<nucnum;i++){
    den_org.put_data(i,dens[i]);
    //cout<<i<<" "<<last.get_dat(i,i)<<"\n";
  };

  den_new=last*den_org;

  if(print){
    for(int i=0;i<nucnum;i++){
      real tmp1=den_org.get_dat(i);
      real tmp2=den_new.get_dat(i);
      if(tmp1!=0.||tmp2!=0.){
        cout<<midt.Name(nucid[i])<<" "<<nucid[i]<<" : "<<tmp1<<" -> "<<tmp2<<"\n";
      };
    };
  };

  for(int i=0;i<nucnum;i++){
    //dens[i]=den_new.get_dat(i);
    real tmp=den_new.get_dat(i);
    if(tmp<0.)tmp=0;
    dens[i]=tmp;
  };
};

void Burnup::BurnupCalculation(Medium &med,real flx_level,real delt,bool putmed)
{
  if(putmed)PutMediumData(med);

  //BurnupCalculationNew(flx_level,delt,false);
  //BurnupCalculation(flx_level,delt,"pade",false);
  BurnupCalculation(flx_level,delt,"chebyshev",false);

  for(int i=0;i<nucnum;i++){
    real den=GetDensity(i);
    med.GetNuclideInTurn(i).PutDensity(den);
  };
};

void Burnup::BurnupCalculationByKrylov(real flx,real dt)
{
  GroupData2D tmat=trmat_flxindep+trmat_flxdep*(flx*1e-24);

  GroupData1D nuc(nucnum);
  for(int i=0;i<nucnum;i++){
    nuc.put_data(i,dens[i]); 
  };

  // 
  int ord=tmat.CalOrderForMatrixExponential(dt);
  real ordl=log(ord);
  int div=int(ordl);

  real dtdiv=dt/div;
  for(int ii=0;ii<div;ii++){
    nuc=tmat.CalMatrixExponentialByKrylov(nuc,dtdiv,order_krylov);
  };

  for(int i=0;i<nucnum;i++){
    real tmp=nuc.get_dat(i);
    if(tmp<0.)tmp=0.;
    dens[i]=tmp;
  };

};

void Burnup::BurnupCalculationByKrylov(Medium &med,real flx,real dt,bool putmed)
{
  if(putmed)PutMediumData(med);
  BurnupCalculationByKrylov(flx,dt);

  for(int i=0;i<nucnum;i++){
    med.GetNuclideInTurn(i).PutDensity(dens[i]);
  };
};

void Burnup::BurnupCalculationByPade(Medium &med,real flx_level,real delt,bool putmed)
{
  if(putmed)PutMediumData(med);
  BurnupCalculation(flx_level,delt,"chebyshev",false);

  for(int i=0;i<nucnum;i++){
    real den=GetDensity(i);
    med.GetNuclideInTurn(i).PutDensity(den);
  };
};

void Burnup::BurnupCalculationByChebyshev(Medium &med,real flx_level,real delt,bool putmed)
{
  /*
  if(putmed)PutMediumData(med);
  BurnupCalculation(flx_level,delt,"chebyshev",false);
  for(int i=0;i<nucnum;i++){
    med.GetNuclideInTurn(i).PutDensity(dens[i]);
  };
  */
  if(putmed)PutMediumData(med);

  GroupData2D tmat=trmat_flxindep+trmat_flxdep*(flx_level*1e-24);
  GroupData1D nuc(nucnum);
  for(int i=0;i<nucnum;i++){
    nuc.put_data(i,dens[i]); 
  };

  nuc=tmat.CalMatrixExponentialByChebyshev14(nuc,delt);
  //nuc=tmat.CalMatrixExponentialByMMPA32(nuc,delt);

  for(int i=0;i<nucnum;i++){
    real tmp=nuc.get_dat(i);
    //if(tmp<0.)tmp=0.; // Negative number density is reset to zero.
    med.GetNuclideInTurn(i).PutDensity(tmp);
    dens[i]=tmp;
  };

};

void Burnup::BurnupCalculationByMMPA(Medium &med,real flx_level,real delt,bool putmed)
{
  if(putmed)PutMediumData(med);

  GroupData2D tmat=trmat_flxindep+trmat_flxdep*(flx_level*1e-24);
  GroupData1D nuc(nucnum);
  for(int i=0;i<nucnum;i++){
    nuc.put_data(i,dens[i]); 
  };

  //nuc=tmat.CalMatrixExponentialByMMPA18(nuc,delt);
  nuc=tmat.CalMatrixExponentialByMMPA32(nuc,delt);

  for(int i=0;i<nucnum;i++){
    real tmp=nuc.get_dat(i);
    if(tmp<0.)tmp=0.;  // Negative number density is reset to zero.
    med.GetNuclideInTurn(i).PutDensity(tmp);
    dens[i]=tmp;
  };

};

void Burnup::BurnupCalculationByMultiStepCalc(Medium &med,GroupData1D &nuc_int,real flx_level,real delt,bool putmed)
{
  if(putmed)PutMediumData(med);

  GroupData2D tmat=trmat_flxindep+trmat_flxdep*(flx_level*1e-24);
  GroupData1D nuc(nucnum);
  for(int i=0;i<nucnum;i++){
    nuc.put_data(i,dens[i]); 
  };

  vector<GroupData1D> ans(2);
  tmat.MultiStepCalc(nuc,ans,delt,2);

  for(int i=0;i<nucnum;i++){
    real tmp1=ans[1].get_dat(i);
    if(tmp1<0.)tmp1=0.;
    med.GetNuclideInTurn(i).PutDensity(tmp1);
    dens[i]=tmp1;
    real tmp2=ans[0].get_dat(i);
    if(tmp2<0.)tmp2=0.;
    nuc_int.put_data(i,tmp2);    
  };
};

real Burnup::CalPowerNormalizationFactor(GeneralSystem &test,real power)
{
  real sum=0.;
  for(int m=0;m<test.GetTotM();m++){
    sum+=GetIntegratedPowerParMesh(test,m);
  };
  return power/sum;
};

real Burnup::GetIntegratedPowerParMedium(GeneralSystem &test,int m)
{
  real sum=0.;
  for(int i=0;i<test.GetTotM();i++){
    if(test.GetMI().GetFMat(i)==m){
      sum+=GetIntegratedPowerParMesh(test,i);
    };
  };
  return sum;
};

real Burnup::GetIntegratedPowerParMesh(GeneralSystem &sys,int x,int y,int z,bool cap)
{
  return GetIntegratedPowerParMesh(sys,sys.GetMeshID(x,y,z),cap);
};

real Burnup::GetIntegratedPowerParMesh(GeneralSystem &test,int m,bool cap)
{
  real sum=0.;
  int num=test.GetMesh(m).GetMed()->GetNucnum();
  GroupData1D flx=test.GetMesh(m).GetFlux();
  real volume=test.GetMesh(m).GetVolume();
  for(int i=0;i<num;i++){
    int ng=test.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      int id=test.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetMatnum();
      real den=test.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetDensity();
      real fission=test.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetMicxs().GetData1d(sigf)*flx;
      real capture=test.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetMicxs().GetData1d(sigc)*flx;
      if(red.GetFissionEnergy(id)>0.)sum+=den*fission*volume*red.GetFissionEnergy(id);
      if(cap&&red.GetCaptureEnergy(id)>0.)sum+=den*capture*volume*red.GetCaptureEnergy(id);
    };
  }; 
 return sum;
};

real Burnup::GetIntegratedPower(Medium &med,GroupData1D &flx,bool cap)
{
  real sum=0.;
  int num=med.GetNucnum();
  //GroupData1D flx=med.GetFlux();
  for(int i=0;i<num;i++){
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      int id=med.GetNuclideInTurn(i).GetMatnum();
      real den=med.GetNuclideInTurn(i).GetDensity();
      real fission=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf)*flx;
      real capture=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc)*flx;
      if(red.GetFissionEnergy(id)>0.)sum+=den*fission*red.GetFissionEnergy(id);
      if(cap&&red.GetCaptureEnergy(id)>0.)sum+=den*capture*red.GetCaptureEnergy(id);
      //cout<<id<<" "<<red.GetFissionEnergy(id)<<" "<<red.GetCaptureEnergy(id)<<"\n";
    };
  };
  //exit(0);
  return sum;
};

real Burnup::GetIntegratedPower(Medium &med,bool cap)
{
  /*
  real sum=0.;
  int num=med.GetNucnum();
  GroupData1D flx=med.GetFlux();
  for(int i=0;i<num;i++){
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      int id=med.GetNuclideInTurn(i).GetMatnum();
      real den=med.GetNuclideInTurn(i).GetDensity();
      real fission=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf)*flx;
      real capture=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc)*flx;
      if(red.GetFissionEnergy(id)>0.)sum+=den*fission*red.GetFissionEnergy(id);
      if(cap&&red.GetCaptureEnergy(id)>0.)sum+=den*capture*red.GetCaptureEnergy(id);
      //cout<<id<<" "<<red.GetFissionEnergy(id)<<" "<<red.GetCaptureEnergy(id)<<"\n";
    };
  };
  //exit(0);
  return sum;
  */
  return GetIntegratedPower(med,med.GetFlux(),cap);
};

real Burnup::GetIntegratedPower(real total_flux,bool cap)
{
  real sum=0.;
  int num=nucnum;
  for(int i=0;i<num;i++){
    int id=nucid[i];
    real den=dens[i];
    real fission=sgf[i]*total_flux;
    real capture=sgc[i]*total_flux;
    if(red.GetFissionEnergy(id)>0.)sum+=den*fission*red.GetFissionEnergy(id);
    if(cap&&red.GetCaptureEnergy(id)>0.)sum+=den*capture*red.GetCaptureEnergy(id);
      //cout<<id<<" "<<red.GetFissionEnergy(id)<<" "<<red.GetCaptureEnergy(id)<<"\n";
  };
  //exit(0);
  return sum;
};

real Burnup::GetDecayPower(Medium &med,int type)
{
  real sum=0.;
  int num=med.GetNucnum();
  for(int i=0;i<num;i++){
    int id=med.GetNuclideInTurn(i).GetMatnum();
    real den=med.GetNuclideInTurn(i).GetDensity();
    real decay_const=bc.GetDecayConstant(id);
    real decay_energy=bc.GetDecayEnergy(id,type);
    sum+=den*1e24*decay_const*decay_energy;
  };
  return sum;
};

real Burnup::CalTotalPower(GeneralSystem &test,bool cap)
{
  real powsum=0.;
  for(int i=0;i<test.GetTotM();i++){
    real pow=GetIntegratedPowerParMesh(test,i,cap);
    powsum+=pow;
  };
  return powsum;
};

real Burnup::CalMaximumFastNeutronFlux(GeneralSystem &test,real e)
{
  int bnd=-1;
  int grp=test.GetGrp();
  for(int i=0;i<grp;i++){
    real en=test.GetMed(0).GetEnband().get_dat(i+1);
    if(en<e&&bnd==-1)bnd=i;
  };

  real en1=test.GetMed(0).GetEnband().get_dat(bnd);
  real en2=test.GetMed(0).GetEnband().get_dat(bnd+1);
  real letup=log(en1/e);
  real letdn=log(e/en2);
  real weight=letup/(letup+letdn);

  real max_fast_flx=0.;
  for(int i=0;i<test.GetTotM();i++){
     real sumflx=0.;
     for(int j=0;j<bnd;j++){
       sumflx+=test.GetMesh(i).GetFlux().get_dat(j);
     };
     sumflx+=test.GetMesh(i).GetFlux().get_dat(bnd)*weight;
     if(sumflx>max_fast_flx)max_fast_flx=sumflx;
  };

  return max_fast_flx;
};

real Burnup::CalConversionRatio(GeneralSystem &test)
{
  real nume=0.;
  real denom=0.;
  for(int i=0;i<test.GetTotM();i++){
    real vol=test.GetMesh(i).GetVolume();
    if(test.GetMesh(i).GetMed()->ExistNuclide(942400)){
      real den_p0=test.GetMesh(i).GetMed()->GetNuclide(942400).GetDensity();
      nume+=test.GetMesh(i).GetMed()->GetNuclide(942400).GetMicxs().GetData1d(sigc)*
            test.GetMesh(i).GetFlux()*vol*den_p0;
    };
    if(test.GetMesh(i).GetMed()->ExistNuclide(922380)){
      real den_u8=test.GetMesh(i).GetMed()->GetNuclide(922380).GetDensity();
      nume+=test.GetMesh(i).GetMed()->GetNuclide(922380).GetMicxs().GetData1d(sigc)*
            test.GetMesh(i).GetFlux()*vol*den_u8;
    };
    if(test.GetMesh(i).GetMed()->ExistNuclide(942390)){
      real den_p9=test.GetMesh(i).GetMed()->GetNuclide(942390).GetDensity();
      denom+=(test.GetMesh(i).GetMed()->GetNuclide(942390).GetMicxs().GetData1d(sigc)+
              test.GetMesh(i).GetMed()->GetNuclide(942390).GetMicxs().GetData1d(sigf))*
              test.GetMesh(i).GetFlux()*vol*den_p9;
    };
    if(test.GetMesh(i).GetMed()->ExistNuclide(942410)){
      real den_p1=test.GetMesh(i).GetMed()->GetNuclide(942410).GetDensity();
      denom+=(test.GetMesh(i).GetMed()->GetNuclide(942410).GetMicxs().GetData1d(sigc)+
              test.GetMesh(i).GetMed()->GetNuclide(942410).GetMicxs().GetData1d(sigf))*
              test.GetMesh(i).GetFlux()*vol*den_p1;
    };
  };

  return nume/denom;
};

real Burnup::CalWeightOfHeavyNuclideParUnitVolume(Medium &med, string keyword)
{
  int matid_min=900000;
  int matid_max=990000;
  if(keyword=="u"||keyword=="U"){
    matid_min=920000;
    matid_max=930000;
  }else if(keyword=="tru"||keyword=="TRU"){
    matid_min=930000;
    matid_max=990000;
  };
  
  real avo=0.60221367;

  real weight=0.;
  int nucnum=med.GetNucnum();
  for(int i=0;i<nucnum;i++){
    int matid=med.GetNuclideInTurn(i).GetMatnum();
    //if(matid>900000){
    if(matid>=matid_min&&matid<matid_max){      
      if(amd.GetData(matid)>0.){
	//if(atomic_weight[matid]>0.){
        real den=med.GetNuclideInTurn(i).GetDensity();
        real mol=den/avo;
        //weight+=mol*atomic_weight[matid];
        weight+=mol*amd.GetData(matid);
      };
    };
  };

  return weight;
};

real Burnup::CalWeightOfHeavyNuclide(GeneralSystem &sys,int medid1,int medid2)
{
  real weight=0.;
  for(int i=medid1;i<=medid2;i++){
    real vol=sys.GetVolumeParMedium(i);
    real w=CalWeightOfHeavyNuclideParUnitVolume(sys.GetMed(i));
    weight+=w*vol;
  };
  return weight;
};

real Burnup::CalWeightOfAllNuclideParUnitVolume(Medium &med)
{
  real avo=0.60221367;

  real weight=0.;
  int nucnum=med.GetNucnum();
  for(int i=0;i<nucnum;i++){
    int matid=med.GetNuclideInTurn(i).GetMatnum();
    //if(matid>900000){
    //if(atomic_weight[matid]>0.){
    if(amd.GetData(matid)>0.){
        real den=med.GetNuclideInTurn(i).GetDensity();
        real mol=den/avo;
        //weight+=mol*atomic_weight[matid];
        weight+=mol*amd.GetData(matid);
      };
      //};
  };

  return weight;
};

real Burnup::CalWeightOfAllNuclide(GeneralSystem &sys,int medid1,int medid2)
{
  real weight=0.;
  for(int i=medid1;i<=medid2;i++){
    real vol=sys.GetVolumeParMedium(i);
    //real w=CalWeightOfHeavyNuclideParUnitVolume(sys.GetMed(i));
    real w=CalWeightOfAllNuclideParUnitVolume(sys.GetMed(i));
    weight+=w*vol;
  };
  return weight;
};

GroupData2D Burnup::CalTransitionMatrix(Medium &med,real flx_level)
{
  int bgrp=1;
  int bgrp2[]={0};
  bgrp2[0]=med.GetImax()-1;

  int nucnum=med.GetNucnum();
  PutNucnum(nucnum);

  GroupData1D flx;
  flx=med.GetFlux();
  for(int i=0;i<nucnum;i++){
    real sig_f=0.;
    real sig_c=0.;
    real sig_n=0.;
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      GroupData1D fis=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).Cond(flx,bgrp,bgrp2);
      GroupData1D cap=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).Cond(flx,bgrp,bgrp2);
      GroupData1D n2n=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).Cond(flx,bgrp,bgrp2);
      sig_f=fis.get_dat(0);
      sig_c=cap.get_dat(0);
      sig_n=n2n.get_dat(0);
    };
    real den=med.GetNuclideInTurn(i).GetDensity();
    int id=med.GetNuclideInTurn(i).GetMatnum();
    PutNuclideData(i,id,den,sig_f,sig_c,sig_n);
  };

  GroupData2D ret;
  CalTransitionMatrixFluxDependentPart();
  ret=CalTransitionMatrix(flx_level);
  return ret;
};

GroupData2D Burnup::CalTransitionMatrixByGroup(Medium &med,int grp,real flx_level)
{
  int nucnum=med.GetNucnum();
  PutNucnum(nucnum);

  GroupData1D flx;
  flx=med.GetFlux();
  real flxsum=flx.get_sum();
  for(int i=0;i<nucnum;i++){
    real fis=0.;
    real cap=0.;
    real n2n=0.;
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      fis=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(grp)/flxsum;
      cap=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).get_dat(grp)/flxsum;
      n2n=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).get_dat(grp)/flxsum;
    };
    real den=med.GetNuclideInTurn(i).GetDensity();
    int id=med.GetNuclideInTurn(i).GetMatnum();
    PutNuclideData(i,id,den,fis,cap,n2n);
  };

  GroupData2D ret;

  CalTransitionMatrixFluxDependentPart();
  ret=CalTransitionMatrix(flx_level,false);
  //ret.ReducedForm();
  return ret;
};

GroupData2D Burnup::CaldTMatdFlux(Medium &med,int grp,int nn)
{
  int nucnum=nn;
  if(nn==-1)nucnum=med.GetNucnum();

  PutNucnum(nucnum);

  for(int i=0;i<nucnum;i++){
    real fis=0.;
    real cap=0.;
    real n2n=0.;
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      fis=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(grp);
      cap=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).get_dat(grp);
      n2n=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).get_dat(grp);
    };
    //real den=med.GetNuclideInTurn(i).GetDensity();
    int id=med.GetNuclideInTurn(i).GetMatnum();
    //PutNuclideData(i,id,den,fis,cap,n2n);
    PutNuclideData(i,id,0.,fis,cap,n2n);
  };

  GroupData2D ret;

  CalTransitionMatrixFluxDependentPart();
  ret=CalTransitionMatrix(1.,false); // false:no-decay
  return ret;
};


void Burnup::PutMediumData(Medium &med)
{
  int bgrp=1;
  int bgrp2[]={0};
  bgrp2[0]=med.GetImax()-1;

  int nucnum=med.GetNucnum();
  PutNucnum(nucnum);

  GroupData1D flx;
  flx=med.GetFlux();
  for(int i=0;i<nucnum;i++){
    real sig_f=0.;
    real sig_c=0.;
    real sig_n=0.;
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      GroupData1D fis=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).Cond(flx,bgrp,bgrp2);
      GroupData1D cap=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).Cond(flx,bgrp,bgrp2);
      GroupData1D n2n=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).Cond(flx,bgrp,bgrp2);
      sig_f=fis.get_dat(0);
      sig_c=cap.get_dat(0);
      sig_n=n2n.get_dat(0);
    };
    /*
    if(ng==-1&&med.GetNuclideInTurn(i).GetMatnum()==6446){
      sig_c=1000.;
    };
    */
    real den=med.GetNuclideInTurn(i).GetDensity();
    int id=med.GetNuclideInTurn(i).GetMatnum();
    PutNuclideData(i,id,den,sig_f,sig_c,sig_n);
    //cout<<id<<" : "<<fis.get_dat(0)<<" "<<cap.get_dat(0)<<" "<<n2n.get_dat(0)<<"\n";
    //cout<<"Fission xs : "<<id<<" "<<fis.get_dat(0)<<"\n";
    //cout<<"Capture xs : "<<id<<" "<<cap.get_dat(0)<<"\n";
  };

  CalTransitionMatrixFluxDependentPart();
  CalTransitionMatrixFluxInDependentPart();
};

void Burnup::PutDensity(vector<real> &j)
{
  int sz=j.size();
  if(sz!=nucnum){
    cout<<"# Error in Burnup::PutDensity.\n";
    cout<<"# Vector size is inconsistent.\n";
    exit(0);
  };
  for(int i=0;i<nucnum;i++){
    dens[i]=j[i];
  };
};

void Burnup::SetZeroDensity()
{
  for(int i=0;i<nucnum;i++){
    dens[i]=0.;
  };
};

void Burnup::ShowDensity()
{
  int num=0;
  for(int i=0;i<nucnum;i++){
    if(dens[i]>0.){
      cout<<nucid[i]<<" "<<midt.Name(nucid[i])<<" "<<dens[i]<<"\n";
      num++;
    };
  };
  cout<<"# Total number of non-zero-density nuclides : "<<num<<"\n";
};

void Burnup::CalOneGroupCrossSection(Medium &med,GroupData1D &flx)
{
  for(int i=0;i<med.GetNucnum();i++){
    int id=med.GetNuclideInTurn(i).GetMatnum();
    int bu_id=SearchNuclide(id);
    if(bu_id!=-1){
      real flxsum_inv=1./flx.get_sum();
      real tmp1=0.;
      real tmp2=0.;
      real tmp3=0.;
      if(med.GetNuclideInTurn(i).GetGrp()!=-1){
        for(int g=0;g<med.GetImax();g++){
  	  real flxg=flx.get_dat(g);
	  tmp1+=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(g)*flxg;
	  tmp2+=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).get_dat(g)*flxg;
	  tmp3+=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).get_dat(g)*flxg;
        };
        tmp1*=flxsum_inv;
        tmp2*=flxsum_inv;
        tmp3*=flxsum_inv;
      };
      sgf[bu_id]=tmp1;
      sgc[bu_id]=tmp2;
      sgn2n[bu_id]=tmp3;
    };
  };
};


void Burnup::PutMediumDataWithZeroCrossSection(Medium &med)
{
  int nucnum=med.GetNucnum();
  PutNucnum(nucnum);

  for(int i=0;i<nucnum;i++){
    real den=med.GetNuclideInTurn(i).GetDensity();
    int id=med.GetNuclideInTurn(i).GetMatnum();
    PutNuclideData(i,id,den,0.,0.,0.);
  };

  trmat_flxdep.set_zero();
  CalTransitionMatrixFluxInDependentPart();
};

void Burnup::CalTransitionMatrixFluxDependentPart()
{
  trmat_flxdep.set_zero();

  for(int i=0;i<nucnum;i++){
    int id=nucid[i];
    // absorption
    trmat_flxdep.add_data(i,i,-(sgf[i]+sgc[i]+sgn2n[i]));
    int tmp=bc.GetNdivCapture(id);
    for(int j=0;j<tmp;j++){
      int id2=bc.GetNextIDCapture(id,j);
      int pos=SearchNuclide(id2);
      if(pos!=-1){
        real rat=bc.GetRatioCapture(id,j);
	trmat_flxdep.add_data(pos,i,sgc[i]*rat);
      };
    };
    // fission
    int tmp2=bc.GetNdivFission(id);
    for(int j=0;j<tmp2;j++){
      int id2=bc.GetNextIDFission(id,j);
      int pos=SearchNuclide(id2);
      if(pos!=-1){
        real rat=bc.GetRatioFission(id,j);
	trmat_flxdep.add_data(pos,i,sgf[i]*rat);
      };
    };
    // n2n
    int tmp3=bc.GetNdivN2N(id);
    for(int j=0;j<tmp3;j++){
      int id2=bc.GetNextIDN2N(id,j);
      int pos=SearchNuclide(id2);
      if(pos!=-1){
        real rat=bc.GetRatioN2N(id,j);
	trmat_flxdep.add_data(pos,i,sgn2n[i]*rat);
      };
    };
  };
};

void Burnup::CalTransitionMatrixFluxInDependentPart()
{
  if(!bc.DecayConstantIsRead()){
    cout<<"# Error in Burnup::CalTransitionMatrixFluxInDependentPart.\n";
    cout<<"# Decay constant data has not yet read.\n";
    exit(0);
  };

  trmat_flxindep.set_zero();

  for(int i=0;i<nucnum;i++){
    int id=nucid[i];
    // Decay
    real dcon=bc.GetDecayConstant(id);
    if(dcon!=0.){
      trmat_flxindep.add_data(i,i,-dcon);
      int tmp4=bc.GetNdivDecay(id);
      for(int j=0;j<tmp4;j++){
        int id2=bc.GetNextIDDecay(id,j);
        int pos=SearchNuclide(id2);
        if(pos!=-1){
          real rat=bc.GetRatioDecay(id,j);
  	  trmat_flxindep.add_data(pos,i,dcon*rat);
        };
      };
    };
  };
};

void Burnup::UpdateTransitionMatrixFluxInDependentPartForSpontaneousFission(int id,real half_life)
{
  for(int i=0;i<nucnum;i++){
    if(id==nucid[i]){
      real dcon=0.6931471806/half_life;
      trmat_flxindep.add_data(i,i,-dcon);
      // fission
      int tmp2=bc.GetNdivFission(id);
      for(int j=0;j<tmp2;j++){
        int id2=bc.GetNextIDFission(id,j);
        int pos=SearchNuclide(id2);
        if(pos!=-1){
          real rat=bc.GetRatioFission(id,j);
	  trmat_flxindep.add_data(pos,i,dcon*rat);
        };
      };
    };
  };
};

void Burnup::CalKermaFactor(Medium &med,bool cap)
{
  int grp=med.GetImax();
  int num=med.GetNucnum();
  
  for(int g=0;g<grp;g++){
    real sum=0.;
    for(int i=0;i<num;i++){
      int ng=med.GetNuclideInTurn(i).GetGrp();
      if(ng!=-1){
        int id=med.GetNuclideInTurn(i).GetMatnum();
        real den=med.GetNuclideInTurn(i).GetDensity();
        real fission=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(g);
        real capture=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).get_dat(g);
        if(red.GetFissionEnergy(id)>0.)sum+=den*fission*red.GetFissionEnergy(id);
        if(cap&&red.GetCaptureEnergy(id)>0.)sum+=den*capture*red.GetCaptureEnergy(id);
      };
    };
    med.GetMacxs().GetData1d(kerma).put_data(g,sum);
  };
};

// +++++

//void Burnup::SetFissionYieldAsInitialNumberDensity(string fisnuc)
void Burnup::SetFissionYieldAsInitialNumberDensity(int matid)
{
  SetZeroDensity();
  int iz,ia,il;
  midt.GetParameter(matid,iz,ia,il);
  //midt.GetParameter(fisnuc,iz,ia,il);
  NuclideChainData ncd=bc.GetNuclideChainData(midt.ID(iz,ia,il));

  int ch=ncd.GetNdiv(0);
  for(int i=0;i<ch;i++){
    int id=ncd.GetIDnext(0,i);
    int jj=SearchNuclide(id);
    if(jj!=-1){
      real yld=ncd.GetRatio(0,i);
      dens[jj]=yld;
    };
  };
};

void Burnup::PutNuclideDataFromBurnupChain()
{
  int num=0;
  vector<int> nid;

  map<int,NuclideChainData>::iterator it=bc.GetNuclideChainData().begin();
  while(it!=bc.GetNuclideChainData().end()){
    int tmp=(*it).first;
    num++;
    nid.push_back(tmp);
    it++;
  };

  PutNucnum(num);
  for(int i=0;i<num;i++){
    PutNuclideData(i,nid[i],0.,0.,0.,0.);
  };
};

void Burnup::CoolingCalculation(real init_dt,int step,int sub_step)
{
  CalTransitionMatrixFluxInDependentPart();
  //cout<<"# ... matrix exponential calculation ...\n";
  GroupData2D matexp=trmat_flxindep.CalMatrixExponentialByChebyshev14(init_dt);
  //GroupData2D matexp=trmat_flxindep.CalMatrixExponentialByChebyshev(init_dt,14);
  //GroupData2D matexp=trmat_flxindep.CalMatrixExponentialByPade(init_dt/sub_step);

  /*
  cout.setf(ios::scientific);
  cout.precision(10);
  for(int i=0;i<nucnum;i++){
    for(int j=0;j<nucnum;j++){
      cout<<matexp.get_dat(i,j)<<"\n";
    };
  };
  exit(0);
  cout<<"... end.\n";
  */

  vector<GroupData1D> density(step*sub_step+1);
  vector<real> t(step*sub_step+1);
  density[0].put_imax(nucnum);
  for(int i=0;i<nucnum;i++){
    density[0].put_data(i,dens[i]);
  };

  t[0]=0.;
  real dt=init_dt/sub_step;
  int cstep=0;
  for(int i=0;i<step;i++){
    cout<<"time step : "<<i<<"\n";
    //GroupData2D matexp=trmat_flxindep.CalMatrixExponentialByChebyshev16(dt);
    //GroupData2D matexp=trmat_flxindep.CalMatrixExponentialByPade(init_dt);
    /*
  cout.setf(ios::scientific);
  cout.precision(10);
  for(int ii=0;ii<nucnum;ii++){
    for(int j=0;j<nucnum;j++){
      if(matexp.get_dat(ii,j)!=0.)cout<<ii<<" "<<j<<" "<<matexp.get_dat(ii,j)<<"\n";
    };
  };
    */

  //exit(0);
    for(int j=0;j<sub_step;j++){
      t[cstep+1]=t[cstep]+dt;
      density[cstep+1]=matexp*density[cstep];
      cstep++;
    };
    dt*=2;
    matexp=matexp*matexp;
  };

  // +++ Decay heat calculation
  cout<<"#\n";
  cout<<"# Decay heat f(t)*t [MeV/fission]\n";
  cout<<"# t[s]       Beta         Gamma        Heavy-part.  Total\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<step*sub_step+1;i++){
    real en=0.;
    real en_b=0.;
    real en_g=0.;
    real en_h=0.;
    for(int j=0;j<nucnum;j++){
      int id=nucid[j];
      real lambda=GetDecayConstant(id);
      real factor=lambda*density[i].get_dat(j)*1e-6; // [MeV]
      for(int k=0;k<3;k++){
	en+=GetBurnupChain().GetDecayEnergy(id,k)*factor;
      };
      en_b+=GetBurnupChain().GetDecayEnergy(id,0)*factor;
      en_g+=GetBurnupChain().GetDecayEnergy(id,1)*factor;
      en_h+=GetBurnupChain().GetDecayEnergy(id,2)*factor;
    }; 
    cout<<t[i]<<"  "<<en_b*t[i]<<"  "<<en_g*t[i]<<"  "<<en_h*t[i]<<"  "<<en*t[i]<<"\n";
  };
  cout<<"\n";

  /*
  for(int i=0;i<nucnum;i++){
    cout<<i<<" "<<density[step].get_dat(i)<<"\n";
  };
  */

  /*
  for(int i=0;i<step+1;i++){
    cout<<i<<" "<<t[i]<<" "<<density[i].get_dat(10)<<" "<<density[i].get_dat(20)<<"\n";
  };
  */

};

void Burnup::CoolingCalculationByKrylov(real init_dt,int step)
{
  CalTransitionMatrixFluxInDependentPart();
  cout<<"... matrix exponential calculation ...\n";

  vector<GroupData1D> density(step+1);
  vector<real> t(step+1);
  density[0].put_imax(nucnum);
  for(int i=0;i<nucnum;i++){
    density[0].put_data(i,dens[i]);
  };

  t[0]=0.;
  real dt=init_dt;
  for(int i=0;i<step;i++){
    cout<<"time step : "<<i<<"\n";
    t[i+1]=t[i]+dt;
    density[i+1]=trmat_flxindep.CalMatrixExponentialByKrylov(density[i],dt,10);
    dt*=2;
  };

  for(int i=0;i<nucnum;i++){
    cout<<i<<" "<<density[step].get_dat(i)<<"\n";
  };

  /*
  for(int i=0;i<step+1;i++){
    cout<<i<<" "<<t[i]<<" "<<density[i].get_dat(10)<<" "<<density[i].get_dat(20)<<"\n";
  };
  */

};

void Burnup::AddNuclideToMediumFromBurnupChain(Medium &med)
{
  real temp=300.;
  int nucnum_org=med.GetNucnum();
  vector<int> nucid_org(nucnum_org);
  vector<real> nucden_org(nucnum_org);
  for(int i=0;i<nucnum_org;i++){
    nucid_org[i]=med.GetNuclideInTurn(i).GetMatnum();
    nucden_org[i]=med.GetNuclideInTurn(i).GetDensity();
  };
  if(nucnum_org!=0.)temp=med.GetNuclideInTurn(0).GetTemperature();
  med.NuclideClear();

  for(int iz=1;iz<100;iz++){
    for(int ia=iz;ia<300;ia++){
      for(int il=0;il<3;il++){
        int id=iz*10000+ia*10+il;
        real den=0.;
        for(int j=0;j<nucnum_org;j++){
          if(id==nucid_org[j]){
            den=nucden_org[j];
          };
        };
	/*
	bool inp=false;
	if(den!=0.)inp=true;
	if(bc.GetNuclideChainData().find(id)!=bc.GetNuclideChainData().end()){
	*/
	//if(den!=0.||bc.GetNuclideChainData().find(id)!=bc.GetNuclideChainData().end()){
	if(den!=0.||bc.GetNuclideChainData(id).GetID()==id){
          Nuclide nucinp;
          nucinp.PutMatnum(id);
          nucinp.PutDensity(den);
          nucinp.PutTemperature(temp);
          med.AddNuclide(nucinp);
	};
      };
    };
  };

  // pseudo FP
  int num=4;
  int fpid[]={9910000,9950000,9980000,9990000};
  for(int i=0;i<num;i++){
    int id=fpid[i];
    real den=0.;
    for(int j=0;j<nucnum_org;j++){
      if(id==nucid_org[j]){
        den=nucden_org[j];
      };
    };
    if(den!=0.||bc.GetNuclideChainData().find(id)!=bc.GetNuclideChainData().end()){
      Nuclide nucinp;
      nucinp.PutMatnum(id);
      nucinp.PutDensity(den);
      nucinp.PutTemperature(temp);
      med.AddNuclide(nucinp);
    };
  };


  for(int i=9990000;i<9997000;i++){
    int id=i*10;
    if(bc.GetNuclideChainData().find(id)!=bc.GetNuclideChainData().end()){
      Nuclide nucinp;
      nucinp.PutMatnum(id);
      nucinp.PutDensity(0.);
      nucinp.PutTemperature(temp);
      med.AddNuclide(nucinp);
    };
  };

};
