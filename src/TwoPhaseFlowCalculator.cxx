#include <cstdlib>
#include "TwoPhaseFlowCalculator.h"

TwoPhaseFlowCalculator::TwoPhaseFlowCalculator(real lin, real din)
{
  length=lin;
  diameter=din;
  area=(diameter*diameter)/4.*PI; 

  //real c_f=0.002; // Friction coefficient
  //c_f=0.8; // Friction coefficient
  c_f=0.01; // Friction coefficient

  g=9.80665; // [m/s2]
  //real g=0.; // [m/s2]

  itermax_rho=10;
  eps_rho=1e-7;
};

void TwoPhaseFlowCalculator::PutNumNode(int i)
{
  num_node=i;

  p.resize(num_node);
  rho.resize(num_node-1);
  h.resize(num_node-1);
  x.resize(num_node-1);
  vr.resize(num_node-1);
  temp.resize(num_node-1);
  v.resize(num_node);
  q_ex.resize(num_node-1);

  dz=length/(num_node-1); 
};

void TwoPhaseFlowCalculator::PutConstantPowerDensity(real q_in)
{
  q_ex.resize(num_node-1);
  for(int i=0;i<num_node-1;i++){
    q_ex[i]=q_in;
  };
};

void TwoPhaseFlowCalculator::PutPowerDensity(int i, real q_in)
{
  if(i<0||i>num_node){
    cout<<"# Error in TwoPhaseFlowCalculator::PutPowerDensity.\n";
    cout<<"# The position assignment "<<i<<" is incorrect.\n";
    exit(0);    
  };
  q_ex[i]=q_in;
};

void TwoPhaseFlowCalculator::TransientFixedOutletPressure()
{
  //int n=51; // Number of boundaries;

  // PBL condition
  double flow_rate=0.1736;     // [kg/s]
  double pressure_in=1.4e7;  // [Pa]
  double rho_in=762.47;     // [kg/m3]
  double t_in=280+273.15;    // [K], calculated from pressure and density

  double factor=1.00002; // Perturbation factor of exit pressure

  int itermax_vel=1000;
  double eps_vel=1e-10;
  double adj_factor=1.; // velocity adjusting factor

  int step=1000;
  double dt=0.1; // [sec]

  bool print=false;

  // Node-variables
  vector< vector<double> > pp(step);   // pressure
  /*
  vector< vector<double> > rhorho(step); // density
  vector< vector<double> > hh(step);   // enthalpy
  vector< vector<double> > xx(step);   // quality
  vector< vector<double> > vrvr(step);   // void rate
  */
  vector< vector<double> > vv(step); // velocity [m/s]
  for(int i=0;i<step;i++){
    pp[i].resize(num_node);
    /*
    rhorho[i].resize(num_node-1);
    hh[i].resize(num_node-1);
    xx[i].resize(num_node-1);
    vrvr[i].resize(num_node-1);
    */
    vv[i].resize(num_node);
  };

  double velocity_in=flow_rate/rho_in/area; // [m/s]
  double h_in=htp2_(&t_in,&pressure_in);
  
  vv[0][0]=velocity_in;

  // ... Initial stationary condition
  Run(flow_rate,pressure_in,t_in);


  pp[0]=GetP();
  /*
  rhorho[0]=GetRho();
  hh[0]=GetH();
  xx[0]=GetX();
  vrvr[0]=GetVr();
  */
  vv[0]=GetV();
  
  cout<<"# Initial exit pressure : "<<pp[0][num_node-1]<<" [Pa]\n";  

  // ... Transient 

  // Perturbation in pressure
  double pressure_out=pp[0][num_node-1]; // [Pa]  
  pressure_out*=factor;


  cout<<0.<<" "<<vv[0][0]*area*rho_in<<"\n";

  for(int ii=1;ii<step;ii++){

    if(print)cout<<"# Time step : "<<ii<<"\n";
    vv[ii][0]=vv[ii-1][0]; // initial entrance velocity is assumed one at the previous time step

    double v_in;
    double p_out=0.;
    double v_in_old, p_out_old;

    v_in=vv[ii][0];

    // Loop of entrance velocity
    for(int jj=0;jj<itermax_vel;jj++){

      if(print)cout<<"#    Velocity-loop : "<<jj<<"\n";
      p_out_old=p_out;      
      p_out=RunTransient(vv[ii][0]*area*rho_in, pressure_in, t_in, dt, false);      
      double new_p=p_out;

      if(jj==0){
        vv[ii][0]*=1.0+adj_factor*(new_p/pressure_out-1.0);
        v_in_old=v_in;
        v_in=vv[ii][0];
      }else{
        double tmp;
        if((pressure_out-p_out)*(pressure_out-p_out_old)<0.){
          tmp=(v_in+v_in_old)*0.5;
        }else{
          tmp=v_in+(v_in-v_in_old)*(pressure_out-p_out)/(p_out-p_out_old);
        };
        v_in_old=v_in;
        v_in=tmp;
        vv[ii][0]=tmp;
      };
      //cout<<"         "<<v_in_old<<" "<<new_p<<" "<<pressure_out<<"\n";
      double res2=fabs(1.-new_p/pressure_out);
      if(res2<eps_vel){
        break;
      };
      if(jj==itermax_vel-1){
        cout<<"# Warning : velocity iteration is NOT converged : "<<res2<<"\n";
      //exit(0);
      };

    }; // entrance velocity-loop end

    // ... after convergence
    RunTransient(vv[ii][0]*area*rho_in, pressure_in, t_in, dt, true);

    
    pp[ii]=GetP();
    /*
    rhorho[ii]=GetRho();
    hh[ii]=GetH();
    xx[ii]=GetX();
    vrvr[ii]=GetVr();
    */
    vv[ii]=GetV();

    //cout<<dt*ii<<" "<<pp[ii][n-2]<<" "<<rhorho[ii][n-2]<<" "<<hh[ii][n-2]<<" "<<vv[ii][n-1]<<" "<<v[ii][0]<<"\n";
    cout<<dt*ii<<" "<<vv[ii][0]*area*rho_in<<"\n";

  }; // end of time-step loop
  
};

real TwoPhaseFlowCalculator::RunTransient(real flow_rate_in, real pressure_in, real t_in, real dt, bool conv)
{
  bool print=false;
  double theta=0.;
  
  // ... Quantities at the next time point
  vector<real> pp;        // pressure [Pa] at mesh-center + at exit
  vector<real> rhorho;    // density at mesh-center
  vector<real> hh;        // enthalpy at mesh-center
  vector<real> xx;        // quality at mesh-center
  vector<real> vrvr;      // void rate at mesh-center
  vector<real> temptemp;  // temperature at mesh-center  
  vector<real> vv;        // velocity [m/s] at node point

  pp.resize(num_node);
  rhorho.resize(num_node-1);
  hh.resize(num_node-1);
  xx.resize(num_node-1);
  vrvr.resize(num_node-1);
  temptemp.resize(num_node-1);
  vv.resize(num_node);

  real rho_in=1./vtp2_(&t_in,&pressure_in);
  double h_in=htp2_(&t_in,&pressure_in);
  
  double velocity_in=flow_rate_in/rho_in/area; // [m/s]
  vv[0]=velocity_in;

  // Loop of node
  double v1=vv[0];
  double v1p=v[0];
  double rho1=rho_in;
  double rho1p=rho_in;
  double p1=pressure_in; 
  double h1=h_in;     
  for(int i=0;i<num_node-1;i++){

    double dzz=dz;
    if(i==0)dzz*=0.5;
	
    // density loop
    rhorho[i]=rho[i];
    for(int j=0;j<itermax_rho;j++){
 
      if(print)cout<<"#        density-loop : "<<j<<"\n";

      // mass equation
      vv[i+1]=1./rhorho[i]*(rho1*vv[i]-dz/dt*(rhorho[i]-rho[i]));

      // momentum equation
      double rho_h0=(rhorho[i]+rho1)*0.5; // rho at n+1 and i-1/2 
      double rho_h1=(rho[i]+rho1p)*0.5; // rho at n and i-1/2
      double term1=dzz/dt*(rho_h0*vv[i]-rho_h1*v[i]);
      double term2=rhorho[i]*vv[i+1]*(vv[i+1]+vv[i])*0.5-rho1*vv[i]*(vv[i]+v1)*0.5+dzz*c_f*0.5*rho_h0*vv[i]*vv[i]+dzz*rho_h0*g;
      double term3=rho[i]*v[i+1]*(v[i+1]+v[i])*0.5-rho1p*v[i]*(v[i]+v1p)*0.5+dzz*c_f*0.5*rho_h1*v[i]*v[i]+dzz*rho_h1*g;
      double dp=-term1-theta*term2-(1.-theta)*term3;
      pp[i]=p1+dp;

      // energy equation
      double ave_v1=(vv[i]+v1)*0.5;
      double ave_v1p=(v[i]+v1p)*0.5; // v at n and i+1/2
      double ave_v=(vv[i]+v1)*0.5;      // v at n+1 and i+1/2
      double t1=rhorho[i]/dt+vv[i+1]*rhorho[i]/dz;
      double t2=(pp[i]-p[i])/dt+q_ex[i]-rhorho[i]*vv[i+1]*g+rho[i]*(h[i]+0.5*ave_v1p*ave_v1p)/dt+vv[i]*rho1*(h1+0.5*ave_v1*ave_v1)/dz;
      hh[i]=t2/t1-0.5*ave_v*ave_v;

      xx[i]=quality(pp[i],hh[i]);
      double new_rho;
      vrvr[i]=void_rate(pp[i],xx[i],hh[i],new_rho);
      double res=fabs(new_rho-rhorho[i])/rhorho[i];
      rhorho[i]=new_rho;

      temptemp[i]=get_temp(pp[i], hh[i]);      

      if(res<eps_rho)break;

      if(j==itermax_rho-1){
        cout<<"# Warning : density iteration is NOT converged at node : "<<i<<"\n";
        cout<<"# Residual is "<<res<<"\n";
        //exit(0);
      };
    }; // density-loop end

    rho1=rhorho[i];
    rho1p=rho[i];
    v1=vv[i];
    v1p=v[i];
    p1=pp[i];       
    h1=hh[i];

  }; // node-loop end

  // momentum equation at exit boundary
  double dzz=dz*0.5;
  double dp=-rhorho[num_node-2]*vv[num_node-1]*(vv[num_node-1]+vv[num_node-2])*0.5+rhorho[num_node-2]*vv[num_node-2]*(vv[num_node-2]+vv[num_node-3])*0.5-dzz*c_f*0.5*rhorho[num_node-2]*vv[num_node-1]*vv[num_node-1]-dzz*rhorho[num_node-2]*g;
  pp[num_node-1]=pp[num_node-2]+dp;

  // ... Quantities are overwritten
  if(conv){
    rho=rhorho;
    p=pp;
    h=hh;
    x=xx;
    vr=vrvr;
    v=vv;
    temp=temptemp;
  };
  
  return pp[num_node-1];

};

void TwoPhaseFlowCalculator::Run(real flow_rate_in, real pressure_in, real t_in)
{
  real rho_in=1./vtp2_(&t_in,&pressure_in);
  real h_in=htp2_(&t_in,&pressure_in);

  v[0]=flow_rate_in/rho_in/area; // [m/s]
  rho[0]=rho_in;

  real v1=v[0];
  real rho1=rho_in;  
  real p1=pressure_in;
  real h1=h_in;  
  for(int i=0;i<num_node-1;i++){

    real dzz=dz;
    if(i==0)dzz*=0.5;

    // density loop
    rho[i]=rho1; // initial assumption using the up-stream mesh data
    for(int j=0;j<itermax_rho;j++){

      // mass equation
      v[i+1]=v[i]*rho1/rho[i];

      // momentum equation
      real dp=-rho[i]*v[i+1]*(v[i+1]+v[i])*0.5+rho1*v[i]*(v[i]+v1)*0.5-dzz*c_f*0.5*(rho[i]+rho1)*0.5*v[i]*v[i]-dzz*(rho[i]+rho1)*0.5*g;
      p[i]=p1+dp;

      // energy equation
      real dh_total=dz*(q_ex[i]/(rho1*v[i])-g);      
      real ave_v=(v[i]+v1)*0.5;
      real h_total=(h1+ave_v*ave_v*0.5)+dh_total;
      real ave_v2=(v[i]+v[i+1])*0.5;
      h[i]=h_total-ave_v2*ave_v2*0.5;

      x[i]=quality(p[i],h[i]);
      real new_rho;
      vr[i]=void_rate(p[i],x[i],h[i],new_rho);
      real res=fabs(new_rho-rho[i])/rho[i];
      rho[i]=new_rho;

      temp[i]=get_temp(p[i], h[i]);

      if(res<eps_rho)break;

      if(j==itermax_rho-1)cout<<"# Warning : density iteration is NOT converged : "<<res<<"\n";
    }; // density-loop end

    v1=v[i];
    rho1=rho[i];
    p1=p[i];
    h1=h[i];    
  };

  // momentum equation at exit boundary
  real dzz=dz*0.5;
  real dp=-rho[num_node-2]*v[num_node-1]*(v[num_node-1]+v[num_node-2])*0.5+rho[num_node-2]*v[num_node-2]*(v[num_node-2]+v[num_node-3])*0.5-dzz*c_f*0.5*rho[num_node-2]*v[num_node-1]*v[num_node-1]-dzz*rho[num_node-2]*g;
  p[num_node-1]=p[num_node-2]+dp; // [Pa]

#if 0  
  cout<<"# Initial exit pressure : "<<p[num_node-1]<<" [Pa]\n";  
  cout<<"# IN   "<<pressure_in<<" "<<rho_in<<" "<<h_in<<"\n";
  Printing();
  cout<<"# Exit pressure : "<<p[num_node-1]<<"\n";
#endif  
};

void TwoPhaseFlowCalculator::Printing(int digit)
{
  cout.setf(ios::scientific);
  cout.precision(digit);
  
  real zpos=0.;
  cout<<"# Z        P          Rho        Specific   quality    velocity   Temp      void-ratio\n";
  cout<<"#                                enthalpy\n";
  cout<<"# [m]      [Pa]       [kg/m3]    [J/kg]                [m/s]      [K]\n";
  for(int i=0;i<num_node-1;i++){
    zpos+=dz*0.5;
    cout<<zpos<<" "<<p[i]<<" "<<rho[i]<<" "<<h[i]<<" "<<x[i]<<" "<<v[i]<<" "<<temp[i]<<" "<<vr[i]<<"\n";
    zpos+=dz*0.5;
  };
};

vector<real> TwoPhaseFlowCalculator::GetHeatTransferCoefficient()
{
  real d_e=0.011; // [m]
  real lambda_c=0.545; // [W/m/K]
  real nu=1.21e-7; // [m2/s]
  real a=1.33e-7;  // [m2/s]

  real Pr=nu/a;
  
  vector<real> htc(num_node);
  for(int i=0;i<num_node;i++){
    real Re=d_e*v[i]/nu;
    real nu=0.023*pow(Re,0.8)*pow(Pr,0.4);
    htc[i]=lambda_c/d_e*nu;
  };

  return htc;
};

real TwoPhaseFlowCalculator::get_temp(real p, real h)
{
  real sat=stemp2_(&p);

  real sat_n=sat-0.001;
  real sat_p=sat+0.001;

  real h_l=htp2_(&sat_n,&p);
  real h_g=htp2_(&sat_p,&p);

  /*
  real h_l=htp2_(&sat,&p)*0.99999999999999999999;
  real h_g=htp2_(&sat,&p)*1.00000000000000000001;
  */

  //cout<<h_l<<" "<<h_g<<" "<<h<<"\n";

  real h_l_org=h_l;
  real h_g_org=h_g;

  if(h>=h_l&&h<=h_g)return sat;

  //real delta_t=500.;
  int itermax=1000;
  //real eps=1e-7;
  real eps=1e-6;

  // liquid phase
  real t_h,t_l;
  if(h<h_l){
    t_h=sat;
    t_l=273.2;
  }else{
    //t_h=sat+delta_t;
    t_h=273.15+799;
    t_l=sat;
  };
  
  {
  real h_h=htp2_(&t_h,&p);
  real h_l=htp2_(&t_l,&p);

  if(h>h_h||h<h_l){
    cout<<"# Error \n";
    cout<<"# Temperature cannot be calculated from enthalpy "<<h<<"\n";
    exit(0);
  };

  //cout<<"# Input "<<t_h<<" "<<t_l<<"\n";
  for(int i=0;i<itermax;i++){
    real t_c=(t_h+t_l)*0.5;
    //real t_c=t_l+(h-h_l)/(h_h-h_l)*(t_h-t_l);
    real h_c=htp2_(&t_c,&p);
    real rat;
    cout.setf(ios::scientific);
    cout.precision(10);
    if(h>h_c){
      //cout<<t_l<<" -> "<<t_c<<"\n";
      //cout<<"!! "<<h_l<<" -> "<<h_c<<"\n";
      t_l=t_c;
      h_l=h_c;
      rat=h_l/h;
      //rat=t_l/t_c;
    }else{
      //cout<<t_h<<" -> "<<t_c<<"\n";
      //cout<<h_h<<" -> "<<h_c<<"\n";
      t_h=t_c;
      h_h=h_c;
      rat=h_h/h;
      //rat=t_h/t_c;

    };
    if(fabs(rat-1.)<eps){
      return t_c;
    };
    if(i==itermax-1){
      cout<<"# Temperature calculation from pressure and enthaly is NOT converged.\n";
      //return t_c;
    };
  };

  cout.setf(ios::scientific);
  cout.precision(15);
  cout<<"#Error!\n";
  cout<<"Pressure : "<<p<<"\n";
  cout<<"Enthalpy : "<<h<<"\n";
  cout<<"Enthalpy_top : "<<h_h<<"\n";
  cout<<"Enthalpy_low : "<<h_l<<"\n";
  cout<<"Temperature top : "<<t_h<<" "<<htp2_(&t_h,&p)<<"\n";
  cout<<"Temperature low : "<<t_l<<" "<<htp2_(&t_l,&p)<<"\n";
  cout<<"h_g_org : "<<h_g_org<<"\n";
  cout<<"h_l_org : "<<h_l_org<<"\n";
  exit(0);

  };

  return 0.;
};

real TwoPhaseFlowCalculator::quality(real p, real h) 
{
  real sat=stemp2_(&p);

  real sat_n=sat-0.001;
  real sat_p=sat+0.001;

  real h_l=htp2_(&sat_n,&p);
  real h_g=htp2_(&sat_p,&p);

  if(h<h_l){
    return 0.; // pure-liquid
  }else if(h>h_g){
    return 1.; // pure-gas
  };

  return (h-h_l)/(h_g-h_l);
};

real TwoPhaseFlowCalculator::void_rate(real p, real x, real h, real &density)
{

  if(x<0.00000001||x>0.999999999){
    real t=get_temp(p,h);
    density=1./vtp2_(&t,&p);
    if(x<0.1){
      return 0.;
    }else{
      return 1.;
    };
  };

  real sat=stemp2_(&p);

  real sat_n=sat-0.001;
  real sat_p=sat+0.001;

  real rho_l=1./vtp2_(&sat_n,&p);
  real rho_g=1./vtp2_(&sat_p,&p);

  real vr=x*rho_l/(x*(rho_l-rho_g)+rho_g);
  //cout<<p<<" "<<sat<<" "<<rho_l<<" "<<rho_g<<"\n";
  density=rho_l*(1.-vr)+rho_g*vr;

  return vr;
};
real TwoPhaseFlowCalculator::get_specific_volume(real t, real p)
{
  return vtp2_(&t,&p);
};

real TwoPhaseFlowCalculator::get_specific_enthalpy(real t, real p)
{
  return htp2_(&t,&p);
};

real TwoPhaseFlowCalculator::get_saturation_temperature(real p)
{
  return stemp2_(&p);
};

real TwoPhaseFlowCalculator::get_dynamic_viscocity(real t, real p)
{
  return dvtp2_(&t,&p);
};

/*
real TwoPhaseFlowCalculator::get_kinetic_viscocity(real t, real p)
{
  return kvtp2_(&t,&p);
};
*/
