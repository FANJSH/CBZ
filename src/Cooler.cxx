#include <cstdlib>
#include "Cooler.h"
#include "SensitivityData.h"

Cooler::Cooler(Burnup &bui)
{
  bu=bui;
  beta=true;
  alpha=true;
  gamma=true;
  spontaneous_fission=false;
};

void Cooler::PutTotalStep(int i) //total step is defined as num of |
{
  total_step=i+1;
  density.resize(i+1);
  adj_density.resize(i+1);
  int_density.resize(i);  //for integral cal
  int_adj_density.resize(i);
  ave_density.resize(i);
  ave_adj_density.resize(i);
  time.resize(i+1);
};

void Cooler::SpontaneousFissionConsideration(int matid,real hl)
{
  spontaneous_fission=true;
  id_sf=matid;
  hl_sf=hl;
};

GroupData2D Cooler::GetBurnupMatrix()
{
  PutNucnum(bu.GetNucnum());
  bu.CalTransitionMatrixFluxInDependentPart();
  GroupData2D tmp=bu.GetTrmatFlxInDep();
  return tmp;
};

void Cooler::CoolingCalculation(real init_dt,int step,int sub_step)
{
  //int max_mat_mult=1;
  int max_mat_mult=20;

  // Numerical round-off errors appear in the CRAM calculation
  // when a large number of matrix multiplications are done,
  // so the number of matrix multiplications in CRAM is limited 
  // below the "max_mat_mult".

  PutNucnum(bu.GetNucnum());
  bu.CalTransitionMatrixFluxInDependentPart();
  // A transition (burnup) matrix is calculated and stored in the burnup class.
  // Note that the neutron-reaction is NOT considered
  // because it is just for a cooling calculation

  if(spontaneous_fission){
    int pos=bu.SearchNuclide(id_sf);
    if(pos==-1){
      cout<<"# Error in Cooler::CoolingCalculation.\n";
      cout<<"# Nuclide "<<id_sf<<" cannot be found in the burnup chain.\n";
      exit(0);
    };
    bu.UpdateTransitionMatrixFluxInDependentPartForSpontaneousFission(id_sf,hl_sf);
    bu.PutDensity(pos,1.);
  };

  PutTotalStep(step*sub_step);

  density[0].put_imax(nucnum);
  // The size of nuclide vector is set.
  for(int i=0;i<nucnum;i++){
    density[0].put_data(i,bu.GetDensity(i));
    //cout<<i<<" "<<density[0].get_dat(i)<<"\n";
    // Initial number density is set from the instance of "burnup" class.
  };

  time[0]=0.;
  real dt=init_dt/sub_step; // Delta_t
  int cstep=0;

  /*
  for(int i=0;i<1400;i++){
    real tmp=bu.GetTrmatFlxInDep().get_dat(i,1400);
    if(tmp!=0.)cout<<i<<" "<<tmp<<"\n";
  };
  */

  // +++ Explcit matrix exponential calculation
  /*
  //GroupData2D matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByChebyshevNew14(dt);
  //GroupData2D matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByChebyshev14(dt);
  GroupData2D matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByMMPA32(dt); 

  // Matrix Exponential of burnup matrix
  for(int i=0;i<step;i++){
    //cout<<"# time step : "<<i<<"\n";
    for(int j=0;j<sub_step;j++){
      time[cstep+1]=time[cstep]+dt;
      density[cstep+1]=matexp*density[cstep];
      //cout<<i<<" "<<j<<" "<<density[cstep+1].get_dat(400)<<"\n";
      cstep++;
    };
    dt*=2.; // !! this should be 2 because matexp is multiplied twice !!
    if((i+1)%max_mat_mult==0){
      //matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByChebyshevNew14(dt);
      //matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByChebyshev14(dt);
      matexp=bu.GetTrmatFlxInDep().CalMatrixExponentialByMMPA32(dt);
    }else{
      matexp=matexp*matexp;
    };
  };
  */

  // +++ Implicit matrix exponential calculation
  /*
  for(int i=0;i<step;i++){
    for(int j=0;j<sub_step;j++){
      time[cstep+1]=time[cstep]+dt;
      density[cstep+1]=bu.GetTrmatFlxInDep().CalMatrixExponentialByMMPA32(density[cstep],dt);
      //density[cstep+1]=bu.GetTrmatFlxInDep().CalMatrixExponentialByChebyshev14(density[cstep],dt);
      cstep++;
    };
    dt*=2.;
  };
  */

  // +++ Implicit matrix exponential calculation with LU-decompositioned MMPA
  //  
  // Very fast !

  for(int i=0;i<step;i++){
    //cout<<"# step : "<<i<<"\n";
    GroupData2D tmp=bu.GetTrmatFlxInDep();
    tmp.LUDecompositionForMatrixExponentialByMMPA32(dt);
    for(int j=0;j<sub_step;j++){
      //cout<<"# ...   "<<j<<"\n";
      time[cstep+1]=time[cstep]+dt;
      density[cstep+1]=tmp.CalMatrixExponentialByLUDMMPA32(density[cstep]);
      cstep++;
    };
    dt*=2.; 
  };

  //density[cstep-1].show_self();

};

GroupData1D Cooler::GetDensity(int st)
{
  if(st<0||st>total_step){
    cout<<"# Error in Cooler::GetDensity.\n";
    cout<<"# Step ID "<<st<<" is inappropriate.\n";
    exit(0);
  };
  return density[st];
};

real Cooler::GetTime(int st)
{
  if(st<0||st>total_step){
    cout<<"# Error in Cooler::GetTime.\n";
    cout<<"# Step ID "<<st<<" is inappropriate.\n";
    exit(0);
  };
  return time[st];
};

void Cooler::CoolingCalculationWithReaction(real flux_level,real init_dt,int step,int sub_step)
{
  bu.CalTransitionMatrixFluxInDependentPart();
  bu.CalTransitionMatrixFluxDependentPart();    
  
  //int max_mat_mult=1;
  int max_mat_mult=20;

  // Numerical round-off errors appear in the CRAM calculation
  // when a large number of matrix multiplications are done,
  // so the number of matrix multiplications in CRAM is limited 
  // below the "max_mat_mult".

  PutNucnum(bu.GetNucnum());

  PutTotalStep(step*sub_step);

  density[0].put_imax(nucnum);
  // The size of nuclide vector is set.
  for(int i=0;i<nucnum;i++){
    density[0].put_data(i,bu.GetDensity(i));
  };

  time[0]=0.;
  real dt=init_dt/sub_step; // Delta_t
  int cstep=0;

  // +++ Implicit matrix exponential calculation with LU-decompositioned MMPA
  //  

  for(int i=0;i<step;i++){
    GroupData2D tmp=bu.CalTransitionMatrix(flux_level);    
    //GroupData2D tmp=bu.GetTrmatFlxInDep();    
    tmp.LUDecompositionForMatrixExponentialByMMPA32(dt);
    for(int j=0;j<sub_step;j++){
      time[cstep+1]=time[cstep]+dt;
      density[cstep+1]=tmp.CalMatrixExponentialByLUDMMPA32(density[cstep]);
      cstep++;
    };
    dt*=2.; 
  };

};

void Cooler::CoolingCalculationArbitralTimeStep(int step_num,real *time_step)
{
  PutNucnum(bu.GetNucnum());
  bu.CalTransitionMatrixFluxInDependentPart();
  
  PutTotalStep(step_num);
  
  density[0].put_imax(nucnum);
  
  for(int i=0;i<nucnum;i++){
    density[0].put_data(i,bu.GetDensity(i));
  };
  
  time[0]=0.;
  for(int i=1;i<step_num+1;i++){
    time[i]=time_step[i-1];
  };
  
  /*
    for(int i=0;i<step_num;i++){
    real dt=time[i+1]-time[i];
    density[i+1]=bu.GetTrmatFlxInDep().CalMatrixExponentialByMMPA32(density[i],dt);
    };
  */
  for(int i=0;i<step_num;i++){
    real dt=time[i+1]-time[i];
    vector<GroupData1D> ans(2);
    bu.GetTrmatFlxInDep().MultiStepCalc(density[i],ans,dt,2);
    density[i+1]=ans[1];
    int_density[i]=ans[0];
  };
};

void Cooler::Adjoint(int step_num,real *time_step,real *power)
{
  PutNucnum(bu.GetNucnum());
  bu.CalTransitionMatrixFluxInDependentPart();
  bu.CalTransitionMatrixFluxDependentPart();
  
  //    GroupData2D tmat=mat;
  
  adj_density[step_num].put_imax(nucnum);
  vector<bool> decay;
  decay.push_back(beta);
  decay.push_back(gamma);
  decay.push_back(alpha);
  
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    real decay_const=bu.GetBurnupChain().GetDecayConstant(id);
    real energy=0.;
    for(int j=0;j<3;j++){
      if(decay[j]==true){
	energy+=bu.GetBurnupChain().GetDecayEnergy(id,j);
      };
    };
    adj_density[step_num].put_data(i,energy*decay_const);
  };
  
  time[0]=0.;
  for(int i=1;i<step_num+1;i++){
    time[i]=time_step[i-1];
  };
  
  //  adj_density[step_num].show_self();
  

  //normal calc
  /*
    for(int i=step_num;i>0;i--){
    real dt=time[i]-time[i-1];
    adj_density[i-1]=tmat.CalMatrixExponentialByMMPA32(adj_density[i],dt);
    };
  */   

  //sub step calculation
  int ssv=40;
  for(int i=step_num;i>0;i--){
    GroupData2D mat=bu.CalTransitionMatrix(1e24*power[i-1],true); // burnup matrix with flux level of unity
    GroupData2D tmat=mat.GetTransposedMatrix();
    real dt=time[i]-time[i-1];
    vector<GroupData1D> ans(ssv);
    tmat.MultiStepCalc(adj_density[i],ans,dt,ssv);
    adj_density[i-1]=ans[ssv-1];
    //+++averaging+++
    int_adj_density[i-1]=adj_density[i]*(0.5/ssv);
    for(int j=0;j<ssv-1;j++){
    int_adj_density[i-1]=int_adj_density[i-1]+ans[j]/ssv;
    };
    int_adj_density[i-1]=int_adj_density[i-1]+ans[ssv-1]*(0.5/ssv);
    //+++++++++++++++
  };
    
    //  adj_density[299].show_self();
    //  adj_density[200].show_self();
    //  cout<<adj_density[299].get_dat(979)<<"\n";
    //  cout<<adj_density[200].get_dat(979)<<"\n";
  
  
  int show_nuc;
  for(int i=1;i<4;i++){
    switch(i){
    case 1:
      show_nuc=441031;
      //show_nuc=541360;
      //show_nuc=551380;
      break;
    case 2:
      //show_nuc=531360;
      show_nuc=541380;
      break;
    case 3:
      //show_nuc=521360;
      show_nuc=531380; 
      break;
    }
    //int show_nuc=631550;
    //cout<<"Adjoint Number Density of "<<midt.GetName(show_nuc)<<"\n";
    //cout<<"# time[s]   AdjointDen  Density     CF      \n";
    for(int i=0;i<step_num;i++){
      int pos=-1;
      for(int j=0;j<nucnum;j++){
	if(show_nuc==bu.GetNuclideID(j))pos=j;
      };
      cout.setf(ios::scientific);
      cout.precision(5);
      real cf = adj_density[i].get_dat(pos)*density[i+1].get_dat(pos);
      //cout<<pos<<"/n";
      //cout<<time[i+1]<<" "<<adj_density[i].get_dat(pos)<<" "<<density[i+1].get_dat(pos)<<" "<<cf<<"\n";
    };  
  };
};

////

void Cooler::ABGSwitch(bool a, bool b, bool g)
{
  alpha=a;
  beta=b;
  gamma=g;
};

void Cooler::IntegrationForward()
{
  for(int i=0;i<total_step-1;i++){
    ave_density[i].put_imax(nucnum);
    real dt=time[i+1]-time[i];
    for(int k=0;k<nucnum;k++){
      real n1=density[i].get_dat(k);
      real n2=int_density[i].get_dat(k);
      real n3=density[i+1].get_dat(k);
      real tmp;
      if(n1==n3){
	tmp=n1;
      }else if(2*n2==(n1+n3)){
	tmp=n2;
      }else if(n1>n3&&n2>=n1){
	tmp=n1;
      }else if(n1>n3&&n2<=n3){
	tmp=n3;
      }else if(n1<n3&&n2>=n3){
	tmp=n3;
      }else if(n1<n3&&n2<=n1){
          tmp=n1;
      }else {
	real a=-(n1-n2)*(n1-n2)/(2*n2-n1-n3);
	real b=log((n2-n3)/(n1-n2))/(dt*0.5);
	real c=(n2*n2-n1*n3)/(2*n2-n1-n3);
	tmp=(a*exp(b*dt)/b+c*dt-a/b)/dt;
      };
      ave_density[i].put_data(k,tmp);
    };
  };
};

void Cooler::IntegrationAdjoint()
{
  for(int i=0;i<total_step-1;i++){
    ave_adj_density[i].put_imax(nucnum); 
    real dt=time[i+1]-time[i];
    for(int k=0;k<nucnum;k++){
      real n1=adj_density[i].get_dat(k);
      real n2=int_adj_density[i].get_dat(k);
      real n3=adj_density[i+1].get_dat(k);
      real tmp;
      if(n1==n3){
	tmp=n1;
      }else if(2*n2==(n1+n3)){
	tmp=n2;
      }else if(n1>n3&&n2>=n1){
	tmp=n1;
      }else if(n1>n3&&n2<=n3){
	tmp=n3;
      }else if(n1<n3&&n2>=n3){
	tmp=n3;
      }else if(n1<n3&&n2<=n1){
          tmp=n1;
      }else {
	real a=-(n1-n2)*(n1-n2)/(2*n2-n1-n3);
	real b=log((n2-n3)/(n1-n2))/(dt*0.5);
	real c=(n2*n2-n1*n3)/(2*n2-n1-n3);
	tmp=(a*exp(b*dt)/b+c*dt-a/b)/dt;
      };
      ave_adj_density[i].put_data(k,tmp);
    };
  };
};

void Cooler::SensitivityCalculation(int mat,int step_num,real *time_step,real *power,string addname)
{
  string sys_name="fission_pul_fin";
  string para_name="decayheat";
  string lib_name="unknown";
  //string addname="decayheat_pf";

  SensitivityData sns;
  sns.PutName(sys_name,para_name,lib_name);
  sns.PutGroup(1);
  GroupData1D band;
  band.put_imax(1);
  band.put_data(0,0.);
  sns.GetEnband().copy(band);
  PutTotalStep(step_num);

  FiniteIrradiationCalculation(mat,step_num,time_step,power);
  //CoolingCalculationArbitralTimeStep(step_num,time_step); // # of step, time
  Adjoint(step_num,time_step,power);

  /////step_modification/////

  for(int ii=0;ii<step_num;ii++){
    //ave_density[ii].put_imax(nucnum);
    //ave_adj_density[ii].put_imax(nucnum);
    for(int i=0;i<nucnum;i++){
      //arithmetic average
      //ave_density[ii].put_data(i,(density[ii].get_dat(i)+density[ii+1].get_dat(i))*0.5);
      //ave_adj_density[ii].put_data(i,(adj_density[ii].get_dat(i)+adj_density[ii+1].get_dat(i))*0.5);
      //gemetric average
      //ave_density[ii].put_data(i,sqrt(density[ii].get_dat(i)*density[ii+1].get_dat(i)));
      //ave_adj_density[ii].put_data(i,sqrt(adj_density[ii].get_dat(i)*adj_density[ii+1].get_dat(i)));
    };
  };

  //Integral Caintegratlculation
  //IntegrationForward();
  for(int ii=0;ii<step_num;ii++){
    ave_density[ii]=int_density[ii];
  };
  //IntegrationAdjoint();
  for(int ii=0;ii<step_num;ii++){
    ave_adj_density[ii]=int_adj_density[ii];
  };

  // Energy //
  //cout<<"Energy"<<"\n";  
  //real ev_to_j=1.60219e-19;
  real tmp=0.;
  real tmp1;
  real end_decayheat=0.;//[W/cm3]
  vector<vector<real> > energy_sens;
  energy_sens.resize(nucnum);
  vector<bool> decay;
  decay.push_back(beta);
  decay.push_back(gamma);
  decay.push_back(alpha);
  
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    tmp1=0.;
    real decay_const=bu.GetBurnupChain().GetDecayConstant(id);
    real dens=density[step_num].get_dat(i);
    energy_sens[i].resize(3);
    for(int j=0;j<3;j++){
      if(decay[j]==true){
	real energy=bu.GetBurnupChain().GetDecayEnergy(id,j);
	tmp1=tmp1+energy;
	energy_sens[i][j]=energy*decay_const*dens;
	//cout<<"ID,bgaenergy,decay,dens"<<" "<<id<<" "<<j<<" "<<energy<<" "<<decay_const<<" "<<dens<<"\n";
      }else if(decay[j]==false){
	energy_sens[i][j]=0.;
      };
    };
    tmp=tmp1*decay_const;
    end_decayheat=end_decayheat+tmp*dens;
  };
  
  real end_decayheat_inverse=1./end_decayheat;
  real response = end_decayheat;
  sns.PutValue(response);
  
  for(int i=0;i<nucnum;i++){
   int id=bu.GetNuclideID(i);
   real sume=0.;
   for(int j=0;j<3;j++){
     real tmp=energy_sens[i][j];
     energy_sens[i][j]=tmp/end_decayheat;
     sume+=energy_sens[i][j];
     int mt=99990+j;
     sns.PutSensitivity0D(id,mt,energy_sens[i][j]);
   };
   //cout.setf(ios::scientific);
   //cout.precision(5);
   //if(sume!=0.)cout<<id<<" "<<sume<<"\n";
  };
  
  //half life//
  //cout<<"half life"<<"\n";  
  //cout<<"ID"<<"    "<<"sensitivity"<<"\n";
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    real decay_c=bu.GetDecayConstant(id);
    if(decay_c!=0.){
      real ddecay_c=-decay_c;
      real sum=0.;
      for(int ii=0;ii<step_num;ii++){
	real dt=time[ii+1]-time[ii];	
	real den=ave_density[ii].get_dat(i);
	real adj=ave_adj_density[ii].get_dat(i);
	
	//if(id == 441031){
	//cout<<"id,den,adj,dt,decay_c"<<" "<<id<<" "<<density[ii].get_dat(i)<<" "<<adj_density[ii].get_dat(i)<<" "<<dt<<" "<<ddecay_c<<"\n";
	//};
	
	real value=-ddecay_c*den*adj;
	//cout<<value<<" "<<ddecay_c<<" "<<den<<" "<<adj<<"\n";
	int channel=bu.GetBC().GetNdivDecay(id);
	for(int j=0;j<channel;j++){
	  int id2=bu.GetBC().GetNextIDDecay(id,j);
	  int pos=bu.SearchNuclide(id2);
	  if(pos!=-1){
	    real bratio=bu.GetBC().GetRatioDecay(id,j);
	    value+=bratio*ddecay_c*den*ave_adj_density[ii].get_dat(pos);
	    //cout<<"bratio,ddecay_c"<<" "<<bratio<<" "<<ddecay_c<<"\n";
	  };
	  //cout<<value<<"\n";
	};
	value*=dt;
	sum+=value;
      };
      sum*=end_decayheat_inverse;
      //sum=0.;
      
      //+++direct term+++//
      real sns_direct=0.;
      //     for(int j=0;j<nucnum;j++){
      real tmpsum=0.;
      for(int jj=0;jj<3;jj++){
	//	  tmpsum+=energy_sens[j][jj];
	tmpsum+=energy_sens[i][jj];
      };
      sns_direct=tmpsum;
      //};
      sum-=sns_direct;
      //++++++++++++++++//     
      
      sns.PutSensitivity0D(id,8888,sum);
      
      //cout.setf(ios::scientific);
      //cout.precision(5);
      //cout<<id<<" "<<sum<<"\n";
      
    };  
  };

  
  //yield//
  int idfisnum=21;
  int idfisorg[]={
    922340,922350,922360,922370,922380,
    932370,932390,
    942380,942390,942400,942410,942420,
    952410,952420,952421,952430,
    962420,962430,962440,962450,962460,
  };
  
  //cout<<"yield"<<"\n";
  //cout<<"ID"<<"    "<<"sensitivity"<<"\n";
  
  for(int i=0;i<idfisnum;i++){
    int idfis=idfisorg[i];
    int posf=bu.SearchNuclide(idfis);
    int nuct=bu.GetBC().GetNdivFission(idfis);
    for(int ii=0;ii<nuct;ii++){
      int id=bu.GetBC().GetNextIDFission(idfis,ii);
      real ratio=bu.GetBC().GetRatioFission(idfis,ii);
      int pos=bu.SearchNuclide(id);
      if(posf!=-1){
	real val=0.;
	for(int j=0;j<step_num;j++){
	  real dt=time[j+1]-time[j];
	  val+=(ave_adj_density[j].get_dat(pos)*power[j]*ave_density[j].get_dat(posf)*ratio)*dt;
	};
	val*=end_decayheat_inverse;
	int mt=18000000+idfis;
	sns.PutSensitivity0D(id,mt,val);
	/*
	if(val!=0.){
	  cout.setf(ios::scientific);
	  cout.precision(5);
	  cout<<idfis<<" "<<id<<" "<<val<<"\n";
	};  
	*/
      };
    };
  };
  
  // branching ratio //
  //cout<<"branching ratio"<<"/n";
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);  
    real decay_c=bu.GetDecayConstant(id);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(id);
      if(channel>1){
	vector<real> sns_tmp1;
	vector<real> ratio;
	sns_tmp1.resize(channel);
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(id,j);
	  ratio[j]=rat;

	  for(int k=0;k<step_num;k++){
	    real dt=time[k+1]-time[k];
	    real den=ave_density[k].get_dat(i);
	    int id2=bu.GetBC().GetNextIDDecay(id,j);
	    int pos=bu.SearchNuclide(id2);
	    real val=0.;
	    if(pos!=-1){
	      val+=ave_adj_density[k].get_dat(pos)*decay_c*rat*den;
	    };
	    val*=dt;
	    sum+=val;
	  };
	  sum=sum*end_decayheat_inverse;
	  sns_tmp1[j]=sum;
	};
	
	vector<real> sns_tmp2;
	sns_tmp2.resize(channel);
	for(int j=0;j<channel;j++){
	  real tmps=0.;
	  for(int jj=0;jj<channel;jj++){
	    tmps+=sns_tmp1[jj];
	  };
	  sns_tmp2[j]=sns_tmp1[j]-tmps*ratio[j];
	};
	
	for(int j=0;j<channel;j++){
	  int mt=88880+j;
	  sns.PutSensitivity0D(id,mt,sns_tmp2[j]);
	};
      };
    };
  	      
    };
  string filename="sns."+addname;
  sns.WriteFile("./",filename);
  cout<<"# Sensitivity data is stored in a file ["<<filename<<"]\n";
 
};
 
void Cooler::FiniteIrradiationCalculation(int mat,real dt)
{
  // mat : ID of fissile nuclide 
  
  PutNucnum(bu.GetNucnum()); // The number of nuclides is taken from "Burnup" class.
  
  int mat_turn=-1;
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    bu.PutNuclideData(i,id,0.,0.,0.,0.); // initial ND/sf/sc/sn2n
    if(mat==bu.GetNuclideID(i)){
      mat_turn=i;
      bu.PutNuclideData(i,id,1.,1.,0.,0.); // initial ND/sf/sc/sn2n
    };
  };
  if(mat_turn==-1){
    cout<<"# Error in Cooler::FiniteIrradiationCalculation.\n";
    cout<<"# Nuclide "<<mat<<" is NOT registered to Cooler.\n";
    exit(0);
  };

  GroupData1D den(nucnum);
  den.set_zero();
  den.put_data(mat_turn,1.); // number density is given only to the fissile nuclide

  bu.CalTransitionMatrixFluxInDependentPart();
  bu.CalTransitionMatrixFluxDependentPart();

  GroupData2D trmat=bu.CalTransitionMatrix(1.*1e24,true); // burnup matrix with flux level of unity
  trmat.put_data(mat_turn,mat_turn,0.); // fissile element is NOT depleted

  den=trmat.CalMatrixExponentialByMMPA18(den,dt);

  for(int i=0;i<nucnum;i++){
    bu.PutDensity(i,den.get_dat(i));
  };
};

void Cooler::FiniteIrradiationCalculation(int mat,int step_num, real *time_step,real *power)
{
  // mat : ID of fissile nuclide 
  time[0]=0.;
  for(int i=1;i<step_num+1;i++){
    time[i]=time_step[i-1];
  };

  PutNucnum(bu.GetNucnum()); // The number of nuclides is taken from "Burnup" class.
  
  int mat_turn=-1;
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    bu.PutNuclideData(i,id,0.,0.,0.,0.); // initial ND/sf/sc/sn2n
    if(mat==bu.GetNuclideID(i)){
      mat_turn=i;
      bu.PutNuclideData(i,id,1.,1.,0.,0.); // initial ND/sf/sc/sn2n
    };
  };
  if(mat_turn==-1){
    cout<<"# Error in Cooler::FiniteIrradiationCalculation.\n";
    cout<<"# Nuclide "<<mat<<" is NOT registered to Cooler.\n";
    exit(0);
  };

  density[0].put_imax(nucnum);
  density[0].set_zero();
  density[0].put_data(mat_turn,1.); // number density is given only to the fissile nuclide

  bu.CalTransitionMatrixFluxInDependentPart();
  bu.CalTransitionMatrixFluxDependentPart();

  int ssv=40;
  for(int i=0;i<step_num;i++){

    GroupData2D trmat=bu.CalTransitionMatrix(1e24*power[i],true); // burnup matrix with flux level of unity
    trmat.put_data(mat_turn,mat_turn,0.); // 

    real dt=time[i+1]-time[i];
    /*
    GroupData2D matexp=trmat.CalMatrixExponentialByChebyshev14(dt);
    density[i+1]=matexp*density[i];
    */

    //Multistep(Averaging)    
    vector<GroupData1D> ans(ssv);
    trmat.MultiStepCalc(density[i],ans,dt,ssv);
    density[i+1].put_imax(nucnum);
    int_density[i].put_imax(nucnum);
    density[i+1]=ans[ssv-1];

    int_density[i]=density[i]*(0.5/ssv);
    for(int j=0;j<ssv-1;j++){
      int_density[i]=int_density[i]+ans[j]/ssv;
    };
    int_density[i]=int_density[i]+ans[ssv-1]*(0.5/ssv);

  };

  /*
  for(int i=0;i<nucnum;i++){
    int idid=bu.GetNuclideID(i);
    cout<<"INVESTIGATION" <<"idid"<<" "<<idid<<" "<<density[1].get_dat(i)<<"\n";
  };
  */
};

void Cooler::ShowRadioactivity(int prt_nuc,string *prt_nuc_nam)
{
  vector<int> prt_nuc_turn(prt_nuc);

  for(int i=0;i<prt_nuc;i++){
    int id=midt.ID(prt_nuc_nam[i]);
    for(int j=0;j<nucnum;j++){
      if(bu.GetNuclideID(j)==id){
	prt_nuc_turn[i]=j;
	break;
      };
    };
    if(prt_nuc_turn[i]==-1){
      cout<<"# Error in Cooler::ShowRadioactivity.\n";
      cout<<"# Nuclide "<<prt_nuc_nam[i]<<" is NOT included in the burnup chain.\n";
      exit(0);
    };
  };

  cout.setf(ios::scientific);
  cout.precision(3);

  cout<<"#\n";
  cout<<"# Radioactivity [Bq]\n";
  cout<<"# Time[s] Time[y]   ";
  for(int j=0;j<prt_nuc;j++){
    WriteOut(prt_nuc_nam[j],10);
  };
  cout<<" All";
  cout<<"\n#\n";

  for(int i=0;i<total_step;i++){
    cout<<time[i]<<" "<<time[i]/(60*60*24*365)<<" ";
    real sum=0.;
    for(int j=0;j<prt_nuc;j++){

      real tmp=0.;
      for(int k=0;k<nucnum;k++){
	string matname=midt.Name(bu.GetNuclideID(k));
	int sz=prt_nuc_nam[j].size();
	//if(i==0)cout<<sz<<" "<<prt_nuc_nam[j]<<" "<<matname.substr(0,sz)<<"\n";
	if(prt_nuc_nam[j]==matname.substr(0,sz)){
          int id=bu.GetNuclideID(k);
          real lambda=bu.GetDecayConstant(id);
          tmp+=lambda*density[i].get_dat(k)*1e24;
	  //cout<<id<<" "<<lambda<<" "<<density[i].get_dat(k)<<"\n";
	};
      };
      cout<<tmp<<" ";
      sum+=tmp;

      /*
      int id=midt.ID(prt_nuc_nam[j]);
      real lambda=bu.GetDecayConstant(id);
      cout<<lambda*density[i].get_dat(prt_nuc_turn[j])<<" "; // [Bq]
      */
    };
    cout<<" "<<sum<<"\n";
  };

};
 
void Cooler::ShowDecayHeat(real tcor, int precision)
{
  cout<<"#\n";
  cout<<"# Decay heat f(t)*t [MeV/fission]\n";
  cout<<"# time[s]    Beta         Gamma        Heavy-Part.  Total\n";
  cout.setf(ios::scientific);
  cout.precision(precision);
  real ent_b=0.;
  real ent_g=0.;
  real ent_h=0.;
  for(int i=0;i<total_step;i++){
    time[i]-=tcor;
    real en=0.;
    real en_b=0.;
    real en_g=0.;
    real en_h=0.;
    for(int j=0;j<nucnum;j++){
      int id=bu.GetNuclideID(j);
      real lambda=bu.GetDecayConstant(id);
      real factor=lambda*density[i].get_dat(j)*1e-6; // [MeV]
      for(int k=0;k<3;k++){
	en+=bu.GetBurnupChain().GetDecayEnergy(id,k)*factor;
      };
      en_b+=bu.GetBurnupChain().GetDecayEnergy(id,0)*factor;
      // for beta component
      en_g+=bu.GetBurnupChain().GetDecayEnergy(id,1)*factor;
      // for gamma component
      en_h+=bu.GetBurnupChain().GetDecayEnergy(id,2)*factor;
      // for heavy particle component
    }; 
    cout<<time[i]<<"  "<<en_b*time[i]<<"  "<<en_g*time[i]<<"  "<<en_h*time[i]<<"  "<<en*time[i]<<"\n";
    if(i>1){
      real th=log(time[i]);
      real tl=log(time[i-1]);
      real dw=th-tl;
      ent_b+=dw*en_b*time[i];
      ent_g+=dw*en_g*time[i];
      ent_h+=dw*en_h*time[i];
    };
  };
  cout<<"\n";

  cout<<"# Time-integrated decay heat\n";
  cout<<"#\n";
  cout<<"#   Beta energy  : "<<ent_b<<" [MeV]\n";
  cout<<"#   Gamma energy : "<<ent_g<<" [MeV]\n";
  cout<<"#   H.P. energy  : "<<ent_h<<" [MeV]\n";
  /*
  for(int i=0;i<nucnum;i++){
    cout<<i<<" "<<density[step].get_dat(i)<<"\n";
  };
  */

  /*
  for(int i=0;i<step+1;i++){
    cout<<i<<" "<<time[i]<<" "<<density[i].get_dat(10)<<" "<<density[i].get_dat(20)<<"\n";
  };
  */
};

void Cooler::ShowRadiationToxity(real tcor)
{
  cout<<"#\n";
  cout<<"# Radiation toxity [Sv]\n";
  cout<<"# time[s]    time[y]      Ingestion    Inhalation\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<total_step;i++){
    time[i]-=tcor;
    real sum_ing=0.;
    real sum_inh=0.;
    for(int j=0;j<nucnum;j++){
      int id=bu.GetNuclideID(j);
      real lambda=bu.GetDecayConstant(id);
      real factor=lambda*density[i].get_dat(j)*1e24;
      sum_ing+=factor*bu.GetBurnupChain().GetDoseCoefficientIngestion(id);
      sum_inh+=factor*bu.GetBurnupChain().GetDoseCoefficientInhalation(id);
    }; 
    if(sum_ing<=1e-100)sum_ing=0.;
    if(sum_inh<=1e-100)sum_inh=0.;
    cout<<time[i]<<"  ";
    cout<<time[i]/(60*60*24*365)<<"  "<<sum_ing<<"  "<<sum_inh<<"\n";
  };
  cout<<"\n";
};

void Cooler::ShowDecayHeatJoule()
{
  cout<<"#\n";
  cout<<"# ABSOLUTE VALUE of decay heat from ABSOLUTE number density data\n";
  cout<<"#\n";
  cout<<"# time[s]   time[d]     time[y]     DecayHeat[W] [Actinide]  [F.P.]\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<total_step;i++){
    real en=0.;
    real en_ac=0.;
    real en_fp=0.;
    for(int j=0;j<nucnum;j++){
      int id=bu.GetNuclideID(j);
      real lambda=bu.GetDecayConstant(id);
      real factor=lambda*density[i].get_dat(j)*1e24*1e-6*mev_to_j;
      for(int k=0;k<3;k++){
	real tmp=bu.GetBurnupChain().GetDecayEnergy(id,k)*factor;
	en+=tmp;
	if(id>=900000){
	  en_ac+=tmp;
	}else{
	  en_fp+=tmp;
	};
      };
    }; 
    cout<<time[i]<<"  "<<time[i]/(60*60*24)<<"  "<<time[i]/(60*60*24*365);
    cout<<"  "<<en<<"  "<<en_ac<<"  "<<en_fp<<"\n";
  };
  cout<<"\n";
};

void Cooler::ShowDecayHeatJoule(int num,string *name)
{
  vector<int> id(num);
  for(int i=0;i<num;i++){
    id[i]=midt.ID(name[i]);
  };

  cout<<"#\n";
  cout<<"# time[s]   time[d]     time[y]";
  for(int j=0;j<num;j++){
    cout<<name[j]<<" ";
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(4);
  for(int i=0;i<total_step;i++){
    cout<<time[i]<<"  "<<time[i]/(60*60*24)<<"  "<<time[i]/(60*60*24*365)<<" ";
    real sum=0.;
    for(int j=0;j<num;j++){
      real lambda=bu.GetDecayConstant(id[j]);
      int pos=GetNucpos(id[j]);
      real factor=lambda*density[i].get_dat(pos)*1e24*1e-6*mev_to_j;
      real en=0.;
      for(int k=0;k<3;k++){
	real tmp=bu.GetBurnupChain().GetDecayEnergy(id[j],k)*factor;
	en+=tmp;
      };
      sum+=en;
      cout<<en<<" ";
    };
    cout<<sum<<"\n";
  };
};

int Cooler::GetNucpos(int id)
{
  for(int i=0;i<nucnum;i++){
    if(bu.GetNuclideID(i)==id)return i;
  };
  cout<<"# Error in Cooler::GetNucpos.\n";
  exit(0);
};

void Cooler::ShowRadioactivityRatio(string nuc1,string nuc2)
{
  int id1=midt.ID(nuc1);
  int id2=midt.ID(nuc2);
  int idt1=-1;
  int idt2=-1;
  for(int i=0;i<nucnum;i++){
    if(bu.GetNuclideID(i)==id1)idt1=i;
    if(bu.GetNuclideID(i)==id2)idt2=i;
  };
  if(idt1==-1||idt2==-1){
    cout<<"# Error in Cooler::ShowRadioactivityRatio.\n";
    cout<<"# Nuclide "<<nuc1<<" or "<<nuc2<<" cannot be found.\n";
    exit(0);
  };

  real lambda1=bu.GetDecayConstant(id1);
  real lambda2=bu.GetDecayConstant(id2);

  cout<<"#\n";
  cout<<"# time[s]  [day]   Radioactivity ratio "<<nuc1<<"/"<<nuc2<<"\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<total_step;i++){
    real bq1=lambda1*density[i].get_dat(idt1);
    real bq2=lambda2*density[i].get_dat(idt2);
    cout<<time[i]<<"  "<<time[i]/(24*60*60)<<" "<<bq1/bq2<<"\n";
  };
  cout<<"\n";
};

void Cooler::ShowNumberDensityForXYPlot(int st)
{
  int zmin=0;
  int nmin=0;
  int zmax=130;
  int nmax=180;
  vector< vector<real> > dat(zmax,vector<real>(nmax,0.));

  int stt=0;
  int edd=total_step-1;
  if(st!=-1){
    stt=st;
    edd=st;
  };

  for(int step=stt;step<=edd;step++){

  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    int at,ms,lv;
    eidt.GetParameterNew(id,at,ms,lv);
    if(lv==0){
      real val=density[step].get_dat(i);
      bool cont=true;
      if(i==nucnum-1)cont=false;
      while(cont){
	int id2=bu.GetNuclideID(i+1);
	int at2,ms2,lv2;
	eidt.GetParameterNew(id2,at2,ms2,lv2);
	if(lv2!=0){
	  val+=density[step].get_dat(i+1);
	  i++;
	}else{
	  cont=false;
	};
      };
      int nn=ms-at;
      dat[at][nn]=density[step].get_dat(i);
    };
  };

  for(int i=nmin;i<nmax;i++){
    for(int jj=0;jj<2;jj++){
      for(int j=zmin;j<zmax;j++){
        cout<<i-0.5+jj<<" "<<j-0.5<<" "<<dat[j][i]<<"\n";
        cout<<i-0.5+jj<<" "<<j+0.5<<" "<<dat[j][i]<<"\n";
      };
      cout<<"\n";
    };
  };
  cout<<"\n\n";

  };

};

void Cooler::ShowDecayEnergyForXYPlot(int st)
{
  int zmin=0;
  int nmin=0;
  int zmax=130;
  int nmax=180;
  vector< vector<real> > dat(zmax,vector<real>(nmax,0.));

  int stt=0;
  int edd=total_step-1;
  if(st!=-1){
    stt=st;
    edd=st;
  };

  for(int step=stt;step<=edd;step++){

  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    int at,ms,lv;
    eidt.GetParameterNew(id,at,ms,lv);
    if(lv==0){
      real en=0.;
      real lambda=bu.GetDecayConstant(id);
      real val=density[step].get_dat(i);
      real factor=lambda*val*1e-6; // [MeV]
      for(int k=0;k<3;k++){
	en+=bu.GetBurnupChain().GetDecayEnergy(id,k)*factor;
      };
      bool cont=true;
      if(i==nucnum-1)cont=false;
      while(cont){
	int id2=bu.GetNuclideID(i+1);
	int at2,ms2,lv2;
	eidt.GetParameterNew(id2,at2,ms2,lv2);
	if(lv2!=0){
	  real lambda=bu.GetDecayConstant(id2);
          real val=density[step].get_dat(i+1);
          real factor=lambda*val*1e-6;
	  for(int k=0;k<3;k++){
	    en+=bu.GetBurnupChain().GetDecayEnergy(id2,k)*factor;
	  };
	  i++;
	}else{
	  cont=false;
	};
      };
      int nn=ms-at;
      dat[at][nn]=en;
    };
  };

  for(int i=nmin;i<nmax;i++){
    for(int jj=0;jj<2;jj++){
      for(int j=zmin;j<zmax;j++){
        cout<<i-0.5+jj<<" "<<j-0.5<<" "<<dat[j][i]<<"\n";
        cout<<i-0.5+jj<<" "<<j+0.5<<" "<<dat[j][i]<<"\n";
      };
      cout<<"\n";
    };
  };
  cout<<"\n\n";

  };

};

void Cooler::ShowDecayGammaSpectrum(DecayGammaSpectrum &dgs)
{
  int ggroup=dgs.GetGgroup();
  GroupData1D gspec(ggroup);

  for(int i=0;i<total_step;i++){
    cout<<"# Time step    : "<<i<<"/"<<total_step<<"\n";
    cout<<"# Elapsed time\n";
    cout<<"#     "<<time[i]<<" [Sec.]\n";
    cout<<"#     "<<time[i]/60.<<" [Min.]\n";
    cout<<"#     "<<time[i]/3600.<<" [Hours]\n";
    cout<<"#     "<<time[i]/(3600.*24.)<<" [Days]\n";
    cout<<"#     "<<time[i]/(3600.*24.*365.)<<" [Years]\n";
    gspec.set_zero();
    for(int j=0;j<nucnum;j++){
      int id=bu.GetNuclideID(j);
      real lambda=bu.GetDecayConstant(id);
      real factor=lambda*density[i].get_dat(j); // [Bq]
      if(dgs.ExistData(id))gspec=gspec+dgs.GetData(id)*factor;
    };
    dgs.ShowGammaSpectrum(gspec);
    cout<<"\n\n";
  };
};

void Cooler::Calculation(int mat,int step_num,real *time_step,real *power){

  PutTotalStep(step_num);
  time[0]=0.;
  for(int i=1;i<step_num+1;i++){
    time[i]=time_step[i-1];
  };

  FiniteIrradiationCalculation(mat,step_num,time_step,power);
};

void Cooler::ShowNumberDensity(int step_inp)
{
  if(step_inp==-1)step_inp=total_step-1;
  
  if(step_inp<0||step_inp>=total_step){
    cout<<"# Error in Cooler::ShowNumberDensity.\n";
    cout<<"# Step number "<<step_inp<<" is inappropriate.\n";
    exit(0);
  };

  cout<<"#\n";
  cout<<"# Number density list at time : "<<time[step_inp]<<" Sec.\n";
  cout<<"#\n";
  cout.setf(ios::scientific);
  cout.precision(20);
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    //cout<<id<<" "<<density[total_step-1].get_dat(i)<<"\n";
    cout<<id<<" "<<density[step_inp].get_dat(i)<<"\n";    
  };
};

void Cooler::ShowNuclideInfo()
{
  cout<<"#\n";
  cout<<"# Nuclide information\n";
  cout<<"\n";
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    int at2,ms2,lv2;
    eidt.GetParameterNew(id,at2,ms2,lv2);
    WriteOut(i,4);
    cout<<" ";
    WriteOut(id,9);
    cout<<"  ";
    WriteOut(midt.Name(id),6);
    cout<<"  ";
    WriteOut(at2,3);
    cout<<"  "; 
    WriteOut(ms2,3);
    cout<<"  ";
    WriteOut(lv2,1);
    cout<<"\n";
  };
};

void Cooler::ShowDecayEnergyInfo()
{
  cout<<"#\n";
  cout<<"# Decay energy information\n";
  cout<<"\n";
  for(int i=0;i<nucnum;i++){
    int id=bu.GetNuclideID(i);
    int at2,ms2,lv2;
    eidt.GetParameterNew(id,at2,ms2,lv2);
    WriteOut(i,4);
    cout<<" ";
    WriteOut(midt.Name(id),6);
    cout<<"  ";
    cout<<bu.GetBurnupChain().GetNuclideChainData(id).GetDecayEnergy(0)<<" ";
    cout<<bu.GetBurnupChain().GetNuclideChainData(id).GetDecayEnergy(1)<<" ";
    cout<<bu.GetBurnupChain().GetNuclideChainData(id).GetDecayEnergy(2)<<" ";    
    cout<<"\n";
  };
};

GData Cooler::GetTimeDependentNumberDensity()
{
  GData gdata(nucnum);
  gdata.PutTagX("second");
  for(int i=0;i<nucnum;i++){
    gdata.PutTagY(i,IntToString(bu.GetNuclideID(i)));
  };

  for(int i=0;i<total_step;i++){
    gdata.push_back_x(time[i]);
    for(int j=0;j<nucnum;j++){
      gdata.push_back_y(j,density[i].get_dat(j));
    };
  };

  return gdata;
};


