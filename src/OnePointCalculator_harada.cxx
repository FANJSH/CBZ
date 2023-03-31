#include "OnePointCalculator.h"

void OnePointCalculator::B1SpectrumCalculationWithFixedSource(Medium &med,GroupData1D &src,real bsq)
{
  if(med.GetPL()==0){
    cout<<"Caution !!!\n";
    cout<<"B1 spectrum calculation cannot be performed\n";
    cout<<"since P1 scattering XS is not given.\n";
    return;
  };

  if(bsq==0.0)bsq=1e-10;

  int imax=med.GetImax();

  real *s0=new real[imax];
  real *s1=new real[imax];
  real *f0=new real[imax];
  real *f1=new real[imax];

  GroupData1D diagp0,diagp1;
  GroupData1D coef0,coef1;

  for(int i=0;i<imax;i++){
    s0[i]=0.; s1[i]=0.;
  };

  diagp0=med.GetMacxs().GetData2d(sigs).get_diag();
  diagp1=med.GetMacxs().GetData2d(sigs,1).get_diag();

  coef0=med.GetMacxs().GetData1d(sigt)-diagp0;
  coef1=med.GetMacxs().GetData1d(sigt,1)*3.0-diagp1;

  for(int i=0;i<imax;i++){
    f0[i]=((s0[i]+src.get_dat(i))*coef1.get_dat(i)-s1[i])
      /(coef0.get_dat(i)*coef1.get_dat(i)+bsq);
    f1[i]=(s1[i]+bsq*f0[i])/coef1.get_dat(i);
    for(int k=0;k<imax;k++){
      s0[k]+=med.GetMacxs().GetData2d(sigs).get_dat(i,k)*f0[i];
      s1[k]+=med.GetMacxs().GetData2d(sigs,1).get_dat(i,k)*f1[i];
    };
  };

  med.GetFlux(0).put_data(f0);
  med.GetFlux(1).put_data(f1);

  delete [] f0;
  delete [] f1;
  delete [] s0;
  delete [] s1;
}

void OnePointCalculator::B1SpectrumCalculation(Medium &med,real bsq)
{
  return B1SpectrumCalculationWithFixedSource(med,med.GetMacxs().GetData1d(chi),bsq);
};

real OnePointCalculator::BucklingSearchB1Method(Medium &med)
{
  real bsq=1e-10;
  real keff,kinf;
  for(int i=0;i<50;i++){
    B1SpectrumCalculation(med,bsq);
    keff=med.CalKeff();
    kinf=med.CalKinf();
    real z=keff/kinf-1.0;
    if(z<1e-10&&z>-1e-10){
      bsq=bsq*1000;
    }
    else{
      real amigra=(kinf/keff-1.0)/bsq;
      bsq=(kinf-1.0)/amigra;
    };
    if(keff<1.0001&&keff>0.9999)break;
  };
  return bsq;
};

void OnePointCalculator::GiveInfiniteDillutionCrossSection(Medium &med,XSLibrary &xslib)
{
  CheckEnergyGroupConsistency(med,xslib);
  int pl=med.GetPL();
  int nucnum=med.GetNucnum();

  // Number density check
  bool nd_non_zero=false;
  for(int i=0;i<nucnum;i++){
    if(med.GetNuclideInTurn(i).GetDensity()>0.)nd_non_zero=true;
  };
  if(!nd_non_zero){
    cout<<"# Error in OnePointCalculator::GiveInfiniteDilutionCrossSection.\n";
    cout<<"# Concerned medium does not contain any nuclides whose ND is larger than 0.0\n";
    cout<<"# (Number densities of all the included nuclides are zero.)\n";
    exit(0);
  };
  
  vector<int> nuc_id_not_include;
  for(int i=0;i<nucnum;i++){
    int id=med.GetNuclideID(i);
    if(xslib.ExistLibData(id)){
      int group=xslib.GetLibData(id).GetGroup();
      med.GetNuclideInTurn(i).PutGrp(group);
      med.GetNuclideInTurn(i).GetMicxs().DataCopyPL(xslib.GetLibData(id).GetXSData(),pl);
      med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,1).copy(med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,0));
    }else{
      nuc_id_not_include.push_back(id);
      //cout<<"# Warning !!! Nuclide "<<id<<" is NOT included in the library.\n";
    };
  };

  int sz=nuc_id_not_include.size();
  if(sz>0&&!xslib.WarningPrint()){
    xslib.WarningPrintTrue();
    cout<<"# Warning !  The following nuclides are NOT included in the library.\n";
    cout<<"#   ";
    for(int i=0;i<sz;i++){
      WriteOut(nuc_id_not_include[i],7);
      cout<<" ";
      if((i+1)%8==0&&i!=sz-1)cout<<"\n#   ";
    };
    cout<<"\n";
  };
  med.GetFlux().copy(xslib.GetWtflux());
  med.CalMacroFromMicro();
  med.GetEnband().copy(xslib.GetEnband());
};

void OnePointCalculator::CalSelfShieldingWithHeterogeneousCorrection
(Medium &med,XSLibrary &xslib,real average_chord,real bell_factor)
{
  GroupData1D tmp(med.GetImax());
  for(int i=0;i<med.GetImax();i++){
    tmp.put_data(i,0.);
  };
  CalSelfShieldingWithDancoffCorrection(med,xslib,average_chord,bell_factor,tmp);
};

void OnePointCalculator::CalSelfShieldingInfiniteSystemWithSig0Correction
(Medium &med,XSLibrary &xslib,real cor)
{
  CheckEnergyGroupConsistency(med,xslib);

  // (No-R-parameter)
  // (Not NRWR)

  GiveInfiniteDillutionCrossSection(med,xslib);

  int nucnum=med.GetNucnum();
  int group=med.GetImax();

  int iter_max=2;
  for(int i=0;i<group;i++){

    for(int iter=0;iter<=iter_max;iter++){

      real sum=0.;
      for(int j2=0;j2<nucnum;j2++){
        real den=med.GetNuclideInTurn(j2).GetDensity();
	if(med.GetNuclideInTurn(j2).GetGrp()!=-1){
	  sum+=den*med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigt).get_dat(i);
	};
      };
      med.GetMacxs().GetData1d(sigt).put_data(i,sum);

      for(int j=0;j<nucnum;j++){
        real den=med.GetNuclideInTurn(j).GetDensity();
        real ng=med.GetNuclideInTurn(j).GetGrp();
	if(den!=0.&&ng!=-1){
          int id=med.GetNuclideID(j);
          real den_inv=1./den;
  	  int sttgrp=xslib.GetLibData(id).GetFtable().GetStartGrpAllMF();
	  int endgrp=xslib.GetLibData(id).GetFtable().GetEndGrpAllMF();
          if(i>=sttgrp&&i<=endgrp){
	    real sig0;
            real micsigt=med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigt).get_dat(i);
            real macsigt=med.GetMacxs().GetData1d(sigt).get_dat(i);
            sig0=(macsigt-den*micsigt+cor)*den_inv;
            real temp=med.GetNuclideInTurn(j).GetTemperature();
  	    bool matrix_cal=false;
	    if(iter==iter_max)matrix_cal=true;
	    vector<real> f(6);
 	    xslib.GetLibData(id).CalF(i,0,0.,0.,temp,sig0,f);
	    med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f,matrix_cal);
	  };
        };
      };
    };
  };
  med.CalMacroFromMicro();
};

void OnePointCalculator::CalSelfShieldingInfiniteSystem(Medium &med,XSLibrary &xslib)
{
  GiveInfiniteDillutionCrossSection(med,xslib);

  GroupData1D tmp(med.GetImax());
  for(int i=0;i<med.GetImax();i++){
    tmp.put_data(i,0.);
  };
  CalSelfShieldingWithDancoffCorrection(med,xslib,1e10,1.0,tmp);
};

void OnePointCalculator::CalSelfShieldingWithDancoffCorrection
(Medium &med,XSLibrary &xslib,real average_chord,real bell,GroupData1D dancoff_correction)
{
  int grp=med.GetImax();
  GroupData1D bell_factor(grp);
  for(int i=0;i<grp;i++){
    bell_factor.put_data(i,bell);
  };
  CalSelfShieldingWithDancoffCorrection(med,xslib,average_chord,bell_factor,dancoff_correction);
};

void OnePointCalculator::CalSelfShieldingWithDancoffCorrection
(Medium &med,XSLibrary &xslib,real average_chord,GroupData1D bell,GroupData1D dancoff_correction)
// Energy-dependent Dancoff
{
  int nucnum=med.GetNucnum();
  int group=med.GetImax();

  // +++ Advanced Bondarenko model option
  //bool cw_correction=false;
  bool cw_correction=true;
  if(average_chord>1e5)cw_correction=false; 
  if(group!=107&&group!=172)cw_correction=false; 
  // if medium is homogeneous or the number of groups is NOT 107, this option is off.
  // Note that most of the CW-total is calculated as FW-total.

  // (P1-total/P0-total)
  // (for 107-group library)
  // the 38th group to the 58th group
  real cw_factor[]={
    2.36741,     2.51232 ,    2.47352 ,    1.38654 ,    4.9612 ,
    3.10858 ,    10.551 ,    1.36165 ,    5.19444 ,    3.51566 ,
    5.88284 ,    0.996285 ,    1.72766 ,    7.18388 ,    1.04683 ,
    4.13905 ,    1.00297 ,    1.00693 ,    1.01314 ,    5.56379 ,
    2.00438 ,
  };  
  real cw_factor_p2=3.15; // pu-242 for 61grp 

  // (for 172-group library)
  // the 58th group to the 90th group
  real cw_factor_172[]={
     2.7387, 2.66749, 1.65458,
     2.09778, 4.11369, 2.90873, 7.26437, 1.50617,  3.828, 2.75492, 0.970958, 6.03587, 1.00189,
     0.995975, 0.990646, 1.00794, 1.46461, 5.15195,  0.980899, 0.998701, 1.00154, 1.03306, 3.53786,
     0.989975, 1.00297, 1.00675, 1.00295, 1.00293, 0.996875, 1.02908, 3.21599, 1.05408, 1.00076, 
  };

  /*
  // (P1-total/total-inf)
  real cw_factor[]={
    2.25412, 2.37487, 2.31219, 1.36464, 4.03517,
    2.6881, 5.73137, 1.33266, 1.87684, 2.95021,
    1.92763, 0.996336, 1.56665, 0.526194, 1.04328,
    0.391435, 1.00196, 1.00531, 1.01177, 0.428908,
    1.58807,
  };
  real cw_factor_p2=2.15; // pu-242 for 61grp
  */


  CheckEnergyGroupConsistency(med,xslib);

  //real black_limit=0.;
  real black_limit=1e-2;
  int wr_nucid=400000; // WR is applied for A>40 nuclides.

  real inv_chord=1./average_chord;

  //GiveInfiniteDillutionCrossSection(med,xslib);
  // ^-- This should be done in prior (20110128)


  vector<bool> exist_rpara(nucnum,false);
  for(int i=0;i<nucnum;i++){
    int id=med.GetNuclideID(i);
    if(xslib.ExistLibData(id)){
      if(xslib.GetLibData(id).GetFtable().GetMaxnr()>1)exist_rpara[i]=true;
    };
  };

  vector<real> f_base(6);
  vector<real> f_calc(6);
  vector<real> f_tmp(6);
  vector<real> f(6);
  vector<int> num_rpara(nucnum);

  
  vector< vector<real> > bgxs;
  if(background_xs_printing){
    bgxs.resize(nucnum);
    for(int i=0;i<nucnum;i++){
      bgxs[i].resize(group,0.);
    };
  };
  
  
  //int iter_max=0;
  int iter_max=2;
  for(int i=0;i<group;i++){

    bool nrwr=true;
    real edwn=xslib.GetEnband().get_dat(i);
    if(edwn>40.)nrwr=false;

    for(int n=0;n<nucnum;n++){
      if(exist_rpara[n]){
	num_rpara[n]=xslib.GetLibData(med.GetNuclideID(n)).GetFtable().GetNumFdata(i);
      };
    };

    for(int iter=0;iter<=iter_max;iter++){

      real macsigt=0.;
      for(int j2=0;j2<nucnum;j2++){
        real ng=med.GetNuclideInTurn(j2).GetGrp();
	if(ng!=-1){
          real den=med.GetNuclideInTurn(j2).GetDensity();
    	  int id=med.GetNuclideID(j2);
	  if(nrwr&&id>=wr_nucid){
	    macsigt+=den*(med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigc).get_dat(i)
	  	         +med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigf).get_dat(i));
	  }else{
            macsigt+=den*med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigt).get_dat(i);
	  };
        };
      };

      for(int j=0;j<nucnum;j++){
        real ng=med.GetNuclideInTurn(j).GetGrp();
        real den=med.GetNuclideInTurn(j).GetDensity();
	  if(((exist_rpara[j]&&den!=0.)||den>1e-6)&&ng!=-1){
	    //if(((exist_rpara[j]&&den!=0.)||den>0.)&&ng!=-1){
          int id=med.GetNuclideID(j);
          real den_inv=1./den;
  	  int sttgrp=xslib.GetLibData(id).GetFtable().GetStartGrpAllMF();
	  int endgrp=xslib.GetLibData(id).GetFtable().GetEndGrpAllMF();
          real bell_factor=bell.get_dat(i);
	  if(xslib.GetLibData(id).ExistBellFactor()){
	    bell_factor=xslib.GetLibData(id).GetBellFactor(i);
	  };
	  if(i>=sttgrp&&i<=endgrp){
  	    real sig0;
            real micsigt=med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigt).get_dat(i);
	    // +++ for NR-WR
	    if(nrwr&&id>=wr_nucid){
	      sig0=(macsigt-den*(med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigc).get_dat(i)
	     	                +med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigf).get_dat(i)))*den_inv;
            }else{
              sig0=(macsigt-den*micsigt)*den_inv;
	    };
	    /*
	    if(id==942420&&iter==iter_max){
	      cout<<"    "<<i<<"\n";
              cout<<" "<<sig0<<"\n";
	    };
	    */
            real temp=med.GetNuclideInTurn(j).GetTemperature();
            real rval=0.;
	    real r_sigt=0.;
            real c=0.;
            if(den*average_chord*micsigt>black_limit){
              c=dancoff_correction.get_dat(i);
	    };
            real sig0c=inv_chord*bell_factor*(1.-c)/(1.+(bell_factor-1.)*c);
	    if(num_rpara[j]<=1){
	      // No R-parameter, or resonance interference with single nuclide
	      //           if(exist_rpara[j]&&num_rpara[j]>0){
              if(exist_rpara[j]){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,0);
		//if(id==922380)cout<<"GRP : "<<i<<" id_ri : "<<id_ri<<"\n";
	        if(med.ExistNuclide(id_ri)){
  	          real den_ri=med.GetNuclide(id_ri).GetDensity();
  	          rval=den_ri*den_inv;
	          if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
  	          sig0-=r_sigt;
	        };
	      };
	      sig0+=sig0c*den_inv;
	      if(iter==iter_max){
  	        xslib.GetLibData(id).CalF(i,0,rval,r_sigt,temp,sig0,f);

		if(background_xs_printing){
		  bgxs[j][i]=sig0+r_sigt;
		};

		// P1/P0 factor
		if(cw_correction){
 	  	  f[3]=0.; // CW-total is same as FW-total in dault.
                  // for uranium-238
		  if(group==107&&id==922380&&i>=37&&i<=57){
		    f[3]=-cw_factor[i-37]; // negative means use of different equation 
		  };
		  if(group==172&&id==922380&&i>=57&&i<=89){
		    f[3]=-cw_factor_172[i-57]; // negative means use of different equation 
		  };
                  // for plutonium-242 (2.38-3.06eV)
		  if(group==107&&id==942420&&i==60){
		    f[3]=-cw_factor_p2;
		  };
		};
	        med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f,true);
	      }else{
	        real fc,ff,fe,fi;
	        fc=xslib.GetLibData(id).GetF(1,i,0,rval,r_sigt,temp,sig0);
	        ff=xslib.GetLibData(id).GetF(0,i,0,rval,r_sigt,temp,sig0);
	        fe=xslib.GetLibData(id).GetF(2,i,0,rval,r_sigt,temp,sig0);
	        fi=xslib.GetLibData(id).GetF(5,i,0,rval,r_sigt,temp,sig0);
	        real xc,xf,xe,xi;
	        xc=xslib.GetLibData(id).GetXSData().GetData1d(sigc).get_dat(i)*fc;
	        xf=xslib.GetLibData(id).GetXSData().GetData1d(sigf).get_dat(i)*ff;
	        xe=xslib.GetLibData(id).GetXSData().GetData1d(sigel).get_dat(i)*fe;
	        xi=xslib.GetLibData(id).GetXSData().GetData1d(siginel).get_dat(i)*fi;
	        med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigt).put_data(i,xc+xf+xe+xi);
	        if(nrwr&&id>=wr_nucid){
		  med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigc).put_data(i,xc);
		  med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigf).put_data(i,xf);
	        };
	      };
	    }else{
	      // +++ Multiple resonance interference treatment +++
	      // base (no-interference)
	      real sig0tmp=sig0+sig0c*den_inv;
  	      xslib.GetLibData(id).CalF(i,0,0.,0.,temp,sig0tmp,f_base);
	      for(int k=0;k<6;k++){f_calc[k]=f_base[k];};
	      // with interference
	      for(int k=0;k<num_rpara[j];k++){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,k);
	        rval=0.;
	        r_sigt=0.;
	        if(med.ExistNuclide(id_ri)){
		  real den_ri=med.GetNuclide(id_ri).GetDensity();
		  rval=den_ri*den_inv;
	          if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
		  sig0tmp=sig0-r_sigt;
  	          sig0tmp+=sig0c*den_inv;
	          xslib.GetLibData(id).CalF(i,k,rval,r_sigt,temp,sig0tmp,f_tmp);
	          for(int l=0;l<6;l++){f_calc[l]*=f_tmp[l]/f_base[l];};
	        };
	      };
	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;

	      // P1/P0 factor
	      if(cw_correction){
 	  	f_calc[3]=0.; // CW-total is same as FW-total in dault.
                // for uranium-238
		if(group==107&&id==922380&&i>=37&&i<=57){
		  f_calc[3]=-cw_factor[i-37];
		};
                // for plutonium-242 (2.38-3.06eV)
		if(group==107&&id==942420&&i==60){
		  f_calc[3]=-cw_factor_p2;
		};
	      };

	      med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f_calc,matrix_cal);

	      if(matrix_cal&&background_xs_printing){
		bgxs[j][i]=sig0+sig0c*den_inv;
	      };

	    };
	  };
        };
      };
    };
  };

  if(background_xs_printing){
    cout.setf(ios::scientific);
    cout.precision(6);
    cout<<" Background cross section (group[0-]/upper energy/BG XS/f-factor(n,g)(n,f)(n,n)/+10% BG XS/+10%ff/-10% BG XS/-10%ff)\n#\n";
    for(int i=0;i<nucnum;i++){
      if(med.GetNuclideInTurn(i).GetGrp()!=-1){
	int s=med.GetNuclideInTurn(i).GetMatnum();
      cout<<"#\n#    Nuclide : "<<s<<"\n#\n";
      cout<<"#group[0-] upperEng    BG XS        f-factor(n,g)  (n,f)       (n,e)        +10% BG XS   ";
      cout<<"+10%ff(n,g)  +10%ff(n,f)  +10%ff(n,e)   -10% BG XS   -10%ff(n,g)  -10%ff(n,f)  -10%ff(n,e)  ";
      cout<<"変動量(n,g)  変動量(n,f)   変動量(n,e) 変動量はabs(+10%ff-(-10%ff)\n";
      real temp=300.;
      for(int g=0;g<group;g++){
	
	//	if(bgxs[i][g]>0.){
	  real e_top=xslib.GetEnband().get_dat(g);
          real f_sigc=xslib.GetLibData(s).GetF(1,g,0,0.,temp,bgxs[i][g]);  // (n,g)
          real f_sigf=xslib.GetLibData(s).GetF(0,g,0,0.,temp,bgxs[i][g]);  // (n,f)
          real f_sigel=xslib.GetLibData(s).GetF(2,g,0,0.,temp,bgxs[i][g]); // (n,e)
	  real f_sigcP10=xslib.GetLibData(s).GetF(1,g,0,0.,temp,bgxs[i][g]*1.1);  // (n,g)
          real f_sigfP10=xslib.GetLibData(s).GetF(0,g,0,0.,temp,bgxs[i][g]*1.1);  // (n,f)
          real f_sigelP10=xslib.GetLibData(s).GetF(2,g,0,0.,temp,bgxs[i][g]*1.1); // (n,e)
	  real f_sigcM10=xslib.GetLibData(s).GetF(1,g,0,0.,temp,bgxs[i][g]*0.9);  // (n,g)
          real f_sigfM10=xslib.GetLibData(s).GetF(0,g,0,0.,temp,bgxs[i][g]*0.9);  // (n,f)
          real f_sigelM10=xslib.GetLibData(s).GetF(2,g,0,0.,temp,bgxs[i][g]*0.9); // (n,e)
	  
	  cout<<"   ";
	  WriteOut(g,4);
          cout<<"  "<<med.GetEnband().get_dat(g)<<"  "<<bgxs[i][g]<<" "<<f_sigc<<" "<<f_sigf<<" "<<f_sigel<<"  ";
	  cout<<bgxs[i][g]*1.1<<" "<<f_sigcP10<<" "<<f_sigfP10<<" "<<f_sigelP10<<"  ";
	  cout<<bgxs[i][g]*0.9<<" "<<f_sigcM10<<" "<<f_sigfM10<<" "<<f_sigelM10<<"　";
	  cout<<abs(f_sigcP10-f_sigcM10)<<" "<<abs(f_sigfP10-f_sigfM10)<<" "<<abs(f_sigelP10-f_sigelM10)<<"\n";
	  //	};
	
      };
      cout<<"\n\n";
      };
    };
  };

  med.CalMacroFromMicro();

  
};

void OnePointCalculator::CalSelfShieldingWithDancoffCorrection
(Medium &med,XSLibrary &xslib,real average_chord,GroupData1D bell,vector<real> dancoff_correction)
// Nuclide-wise Dancoff
{
  CheckEnergyGroupConsistency(med,xslib);

  real black_limit=1e-2;
  int wr_nucid=4000; // WR is applied for A>40 nuclides.

  real inv_chord=1./average_chord;

  GiveInfiniteDillutionCrossSection(med,xslib);

  int nucnum=med.GetNucnum();
  int group=med.GetImax();

  vector<bool> exist_rpara(nucnum,false);
  for(int i=0;i<nucnum;i++){
    int id=med.GetNuclideID(i);
    if(xslib.ExistLibData(id)){
      if(xslib.GetLibData(id).GetFtable().GetMaxnr()>1)exist_rpara[i]=true;
    };
  };

  int iter_max=2;
  for(int i=0;i<group;i++){

    bool nrwr=true;
    real edwn=xslib.GetEnband().get_dat(i+1);
    if(edwn>40.)nrwr=false;
    //if(etop>30.)nrwr=false;

    vector<int> num_rpara(nucnum,0);
    for(int n=0;n<nucnum;n++){
      if(exist_rpara[n]){
	num_rpara[n]=xslib.GetLibData(med.GetNuclideID(n)).GetFtable().GetNumFdata(i);
      };
    };

    for(int iter=0;iter<=iter_max;iter++){

      real sum=0.;

      for(int j2=0;j2<nucnum;j2++){
        real ng=med.GetNuclideInTurn(j2).GetGrp();
	if(ng!=-1){
          real den=med.GetNuclideInTurn(j2).GetDensity();
  	  real micsigt=med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigt).get_dat(i);
	  sum+=den*micsigt;
	};
      };
      med.GetMacxs().GetData1d(sigt).put_data(i,sum);

      for(int j=0;j<nucnum;j++){
        real den=med.GetNuclideInTurn(j).GetDensity();
        real ng=med.GetNuclideInTurn(j).GetGrp();
	if(den>0.&&ng!=-1){
          int id=med.GetNuclideID(j);
          real den_inv=1./den;
  	  int sttgrp=xslib.GetLibData(id).GetFtable().GetStartGrpAllMF();
	  int endgrp=xslib.GetLibData(id).GetFtable().GetEndGrpAllMF();
          real bell_factor=bell.get_dat(i);
	  if(xslib.GetLibData(id).ExistBellFactor()){
	    bell_factor=xslib.GetLibData(id).GetBellFactor(i);
	  };
          if(i>=sttgrp&&i<=endgrp){

	    real sig0;
            real micsigt=med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigt).get_dat(i);
	    // +++ for NR-WR
	    if(nrwr){
	      sig0=0.;
	      for(int k=0;k<nucnum;k++){
                real ng=med.GetNuclideInTurn(k).GetGrp();
	        if(j!=k&&ng!=-1){
	          real den=med.GetNuclideInTurn(k).GetDensity();
	          real micsigt;
		  if(med.GetNuclideID(k)>=wr_nucid){
	  	    micsigt=med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc).get_dat(i)
	  	           +med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf).get_dat(i);
	          }else{
	    	    micsigt=med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigt).get_dat(i);
	          };
	          sig0+=den*micsigt;
	        };
	      };
              sig0*=den_inv;
              }else{
                real macsigt=med.GetMacxs().GetData1d(sigt).get_dat(i);
                sig0=(macsigt-den*micsigt)*den_inv;
	    };
            real temp=med.GetNuclideInTurn(j).GetTemperature();
            real rval=0.;
	    real r_sigt=0.;
            real c=0.;
            if(den*average_chord*micsigt>black_limit){
              c=dancoff_correction[j];
              //c=dancoff_correction.get_dat(i);
	    };
            real sig0c=inv_chord*bell_factor*(1.-c)/(1.+(bell_factor-1.)*c);
	    if(num_rpara[j]<=1){
	      // No-R-parameter, or resonance interference with single nuclide
	      //           if(exist_rpara[j]&&num_rpara[j]>0){
              if(exist_rpara[j]){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,0);
	        if(med.ExistNuclide(id_ri)){
  	          real den_ri=med.GetNuclide(id_ri).GetDensity();
  	          rval=den_ri*den_inv;
	          if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
  	          sig0-=r_sigt;
	        };
	      };
	      sig0+=sig0c*den_inv;

  	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;
	      vector<real> f(6);
	      xslib.GetLibData(id).CalF(i,0,rval,r_sigt,temp,sig0,f);
	      med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f,matrix_cal);
	    }else{
	      // +++ Multiple resonance interference treatment +++
	      vector<real> f_base(6);
	      vector<real> f_calc(6);
	      // base (no-interference)
	      real sig0tmp=sig0+sig0c*den_inv;
  	      xslib.GetLibData(id).CalF(i,0,0.,0.,temp,sig0tmp,f_base);
	      for(int k=0;k<6;k++){f_calc[k]=f_base[k];};
	      // with interference
	      for(int k=0;k<num_rpara[j];k++){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,k);
	        rval=0.;
	        r_sigt=0.;
	        if(med.ExistNuclide(id_ri)){
		  real den_ri=med.GetNuclide(id_ri).GetDensity();
	  	  rval=den_ri*den_inv;
		  if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
		  sig0tmp=sig0-r_sigt;
  	          sig0tmp+=sig0c*den_inv;
	          vector<real> f_tmp(6);
  	          xslib.GetLibData(id).CalF(i,k,rval,r_sigt,temp,sig0tmp,f_tmp);
	          for(int l=0;l<6;l++){f_calc[l]*=f_tmp[l]/f_base[l];};
	        };
	      };
	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;
	      med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f_calc,matrix_cal);
	    };
	  };
        };
      };
    };
  };

  med.CalMacroFromMicro();
};

void OnePointCalculator::CalSelfShieldingWithDancoffCorrection
(Medium &med,XSLibrary &xslib,real average_chord,GroupData1D bell,vector<real> dancoff_correction,int gg)
// Nuclide-wise Dancoff
{
  CheckEnergyGroupConsistency(med,xslib);

  real black_limit=1e-2;
  int wr_nucid=4000; // WR is applied for A>40 nuclides.

  real inv_chord=1./average_chord;

  //GiveInfiniteDillutionCrossSection(med,xslib);

  int nucnum=med.GetNucnum();

  vector<bool> exist_rpara(nucnum,false);
  for(int i=0;i<nucnum;i++){
    int id=med.GetNuclideID(i);
    if(xslib.ExistLibData(id)){
      if(xslib.GetLibData(id).GetFtable().GetMaxnr()>1)exist_rpara[i]=true;
    };
  };

  int iter_max=2;

  int i=gg;

    bool nrwr=true;
    real etop=xslib.GetEnband().get_dat(i);
    if(etop>30.)nrwr=false;

    vector<int> num_rpara(nucnum,0);
    for(int n=0;n<nucnum;n++){
      if(exist_rpara[n]){
	num_rpara[n]=xslib.GetLibData(med.GetNuclideID(n)).GetFtable().GetNumFdata(i);
      };
    };

    for(int iter=0;iter<=iter_max;iter++){

      real sum=0.;

      for(int j2=0;j2<nucnum;j2++){
        real ng=med.GetNuclideInTurn(j2).GetGrp();
	if(ng!=-1){
          real den=med.GetNuclideInTurn(j2).GetDensity();
    	  real micsigt=med.GetNuclideInTurn(j2).GetMicxs().GetData1d(sigt).get_dat(i);
	  sum+=den*micsigt;
	};
      };
      med.GetMacxs().GetData1d(sigt).put_data(i,sum);

      for(int j=0;j<nucnum;j++){
        real den=med.GetNuclideInTurn(j).GetDensity();
        real ng=med.GetNuclideInTurn(j).GetGrp();
	if(den>0.&&ng!=-1){
          int id=med.GetNuclideID(j);
          real den_inv=1./den;
      	  int sttgrp=xslib.GetLibData(id).GetFtable().GetStartGrpAllMF();
	  int endgrp=xslib.GetLibData(id).GetFtable().GetEndGrpAllMF();
          real bell_factor=bell.get_dat(i);
	  if(xslib.GetLibData(id).ExistBellFactor()){
	    bell_factor=xslib.GetLibData(id).GetBellFactor(i);
	  };
          if(i>=sttgrp&&i<=endgrp){

	    real sig0;
            real micsigt=med.GetNuclideInTurn(j).GetMicxs().GetData1d(sigt).get_dat(i);
	    // +++ for NR-WR
	    if(nrwr){
	      sig0=0.;
	      for(int k=0;k<nucnum;k++){
                real ng=med.GetNuclideInTurn(k).GetGrp();
	        if(j!=k&&ng!=-1){
	          real den=med.GetNuclideInTurn(k).GetDensity();
	          real micsigt;
		  if(med.GetNuclideID(k)>=wr_nucid){
	  	    micsigt=med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigc).get_dat(i)
	  	           +med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigf).get_dat(i);
	          }else{
	    	    micsigt=med.GetNuclideInTurn(k).GetMicxs().GetData1d(sigt).get_dat(i);
	          };
	          sig0+=den*micsigt;
	        };
	      };
              sig0*=den_inv;
            }else{
              real macsigt=med.GetMacxs().GetData1d(sigt).get_dat(i);
              sig0=(macsigt-den*micsigt)*den_inv;
  	    };
            real temp=med.GetNuclideInTurn(j).GetTemperature();
            real rval=0.;
	    real r_sigt=0.;
            real c=0.;
            if(den*average_chord*micsigt>black_limit){
              c=dancoff_correction[j];
              //c=dancoff_correction.get_dat(i);
	    };
            real sig0c=inv_chord*bell_factor*(1.-c)/(1.+(bell_factor-1.)*c);
	    if(num_rpara[j]<=1){
	      // No-R-parameter, or resonance interference with single nuclide
	      //           if(exist_rpara[j]&&num_rpara[j]>0){
              if(exist_rpara[j]){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,0);
	        if(med.ExistNuclide(id_ri)){
  	          real den_ri=med.GetNuclide(id_ri).GetDensity();
  	          rval=den_ri*den_inv;
	          if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
  	          sig0-=r_sigt;
	        };
	      };
	      sig0+=sig0c*den_inv;

  	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;
	      vector<real> f(6);
  	      xslib.GetLibData(id).CalF(i,0,rval,r_sigt,temp,sig0,f);
	      med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f,matrix_cal);
	    }else{
	      // +++ Multiple resonance interference treatment +++
	      vector<real> f_base(6);
	      vector<real> f_calc(6);
	      // base (no-interference)
	      real sig0tmp=sig0+sig0c*den_inv;
  	      xslib.GetLibData(id).CalF(i,0,0.,0.,temp,sig0tmp,f_base);
	      for(int k=0;k<6;k++){f_calc[k]=f_base[k];};
	      // with interference
	      for(int k=0;k<num_rpara[j];k++){
	        int id_ri=xslib.GetLibData(id).GetFtable().GetIDRpara(i,k);
	        rval=0.;
	        r_sigt=0.;
	        if(med.ExistNuclide(id_ri)){
		  real den_ri=med.GetNuclide(id_ri).GetDensity();
		  rval=den_ri*den_inv;
		  if(!nrwr||id_ri<wr_nucid){
  	            r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigt).get_dat(i)*rval;
	          }else{
                    r_sigt=med.GetNuclide(id_ri).GetMicxs().GetData1d(sigc).get_dat(i)*rval
                          +med.GetNuclide(id_ri).GetMicxs().GetData1d(sigf).get_dat(i)*rval;
	          };
		  sig0tmp=sig0-r_sigt;
  	          sig0tmp+=sig0c*den_inv;
	          vector<real> f_tmp(6);
  	          xslib.GetLibData(id).CalF(i,k,rval,r_sigt,temp,sig0tmp,f_tmp);
	          for(int l=0;l<6;l++){f_calc[l]*=f_tmp[l]/f_base[l];};
	        };
	      };
	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;
	      med.GetNuclideInTurn(j).CalMicroFromF(i,xslib.GetLibData(id),f_calc,matrix_cal);
   	    };
	  };
        };
      };
    };
};

bool OnePointCalculator::IsFissionMaterial(Medium &med,XSLibrary &xslib)
{
  CheckEnergyGroupConsistency(med,xslib);

  bool fission_material=false;
  for(int i=0;i<med.GetNucnum();i++){
    int ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      int id=med.GetNuclideID(i);
      if(xslib.GetLibData(id).fissile())fission_material=true;
    };
  };
  return fission_material;
};

void OnePointCalculator::CalIncidentEnergyDependentChi(Medium &med,XSLibrary &xslib)
{
  CheckEnergyGroupConsistency(med,xslib);

  bool fission_material=IsFissionMaterial(med,xslib);
  int imax=med.GetImax();
  if(!fission_material){
    for(int i=0;i<imax;i++){
      med.GetMacxs().GetData1d(chi).put_data(i,0.);
    };
    return;
  };
  BucklingSearchB1Method(med);
  GroupData1D flux(imax);
  flux.copy(med.GetFlux());
  CalIncidentEnergyDependentChi(med,xslib,flux);
};

void OnePointCalculator::CalIncidentEnergyDependentChi(Medium &med,XSLibrary &xslib,GroupData1D &flux)
{
  CheckEnergyGroupConsistency(med,xslib);
  bool fission_material=IsFissionMaterial(med,xslib);

  int imax=med.GetImax();
  if(!fission_material){
    for(int i=0;i<imax;i++){
      med.GetMacxs().GetData1d(chi).put_data(i,0.);
    };
    return;
  };

  CalFissionSpectrumMatrix(med,xslib);
  CalEquivalentFissionSpectrumVector(med,flux);
};

void OnePointCalculator::CalEquivalentFissionSpectrumVector(Medium &med,GroupData1D &flux)
{
  int imax=med.GetImax();
  real *sum=new real[imax];
  for(int i=0;i<imax;i++){
    sum[i]=0.;
    for(int j=0;j<imax;j++){
      real nsf_flx=med.GetMacxs().GetData1d(nusigf).get_dat(j)*flux.get_dat(j);
      real chimat=med.GetMacxs().GetData2d(chi).get_dat(j,i);
      sum[i]+=chimat*nsf_flx;
    };
  };

  real tot=0.;
  for(int i=0;i<imax;i++){
    tot+=sum[i];
  };
  tot=1./tot;
  for(int i=0;i<imax;i++){
    sum[i]*=tot;
    med.GetMacxs().GetData1d(chi).put_data(i,sum[i]);
  };

  delete [] sum;
};

void OnePointCalculator::CalFissionSpectrumMatrix(Medium &med,XSLibrary &xslib)
{
  CheckEnergyGroupConsistency(med,xslib);
  bool fission_material=IsFissionMaterial(med,xslib);

  int imax=med.GetImax();
  if(!fission_material){
    for(int i=0;i<imax;i++){
      for(int j=0;j<imax;j++){
        med.GetMacxs().GetData2d(chi).put_data(i,j,0.);
      };
    };
    return;
  };

  GroupData2D inpvec(imax,imax);
  inpvec.set_zero();

  for(int i=0;i<med.GetNucnum();i++){
    real ng=med.GetNuclideInTurn(i).GetGrp();
    if(ng!=-1){
      real den=med.GetNuclideInTurn(i).GetDensity();
      int id=med.GetNuclideID(i);
      if(xslib.GetLibData(id).ExistChiVector()){
        for(int j=0;j<imax;j++){
          real micsigf=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(j);
          real micnu=med.GetNuclideInTurn(i).GetMicxs().GetData1d(nu).get_dat(j);
          real factor=den*micnu*micsigf;
          int vecid=xslib.GetLibData(id).GetChiVector().GetVectorID(j);
          for(int k=0;k<imax;k++){
            real micchivec=xslib.GetLibData(id).GetChiVector(vecid).get_dat(k);
	    inpvec.add_data(j,k,micchivec*factor);
          };
        };
      };
    };
  };

  for(int i=0;i<imax;i++){
    real tot=0.;
    for(int j=0;j<imax;j++){
      tot+=inpvec.get_dat(i,j);
    };
    if(tot!=0.){
      tot=1./tot;
      for(int j=0;j<imax;j++){
	real org=inpvec.get_dat(i,j);
	inpvec.put_data(i,j,org*tot);
      };
    };
  };

  med.GetMacxs().GetData2d(chi).copy(inpvec);
};


void OnePointCalculator::CalThermalScatteringMatrix(Medium &med,XSLibrary &xslib,real cutoff)
{
  CheckEnergyGroupConsistency(med,xslib);
  int grp=med.GetImax();
  int init_grp=-1;
  for(int i=0;i<med.GetImax();i++){
    if(cutoff>med.GetEnband().get_dat(i)&&init_grp==-1){
      init_grp=i;
    };
  };

  int pl_med=med.GetPL();

  for(int i=0;i<med.GetNucnum();i++){
    int grp=med.GetNuclideInTurn(i).GetGrp();
    if(grp!=-1){
      int id=med.GetNuclideID(i);
      if(xslib.GetLibData(id).ExistThermalData()){

        int pl_max=xslib.GetLibData(id).GetThScat().GetPL();
        if(pl_max>pl_med)pl_max=pl_med;

        real tt=med.GetNuclideInTurn(i).GetTemperature();
        int xslib_init_grp=xslib.GetLibData(id).GetThScat().GetInitGrp();
        int xslib_init_src_grp=xslib.GetLibData(id).GetThScat().GetInitSrcGrp();
        int ing=init_grp;
        int ing1=init_grp;
        if(xslib_init_src_grp>init_grp)ing1=xslib_init_src_grp;
        if(xslib_init_grp>init_grp)ing=xslib_init_grp;
        for(int l=0;l<=pl_max;l++){
          for(int g=ing1;g<grp;g++){
	    // +++ mu-bar calculation and overwritten +++++++++++++
            if(xslib.GetLibData(id).GetThScat().GetPL()>0&&pl_max==0){
              real sige0_sum=0.;
	      real sige1_sum=0.;
	      for(int g2=ing;g2<grp;g2++){
		sige0_sum+=xslib.GetLibData(id).GetThScat().GetData(0,g,g2,tt);
		sige1_sum+=xslib.GetLibData(id).GetThScat().GetData(1,g,g2,tt);
	      };
	      med.GetNuclideInTurn(i).GetMicxs().GetData1d(mu).put_data(g,sige1_sum/sige0_sum);
	      // (2l+1) is NOT multiplied in cross sections stored in the library.
	    }; // +++++++++++++++++++++++++++++++++++++++++++++++++
            real sum=0.;
            real sigt0_org=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,0).get_dat(g);
            real sigt1_org=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,1).get_dat(g);
            real sigt0_new=sigt0_org;
            for(int g2=ing;g2<grp;g2++){
	      real tmp=xslib.GetLibData(id).GetThScat().GetData(l,g,g2,tt);
	      tmp*=(2*l+1);
              sum+=tmp;
	      real org=med.GetNuclideInTurn(i).GetMicxs().GetData2d(sigel,l).get_dat(g,g2);
	      med.GetNuclideInTurn(i).GetMicxs().GetData2d(sigel,l).put_data(g,g2,tmp);
	      if(l==0){
	        real delta=tmp-org;
	        sigt0_new+=tmp-org;
  	        //med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,1).add_data(g,delta);
	      };
	    };
	    if(l==0){
              med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigel).put_data(g,sum);
	      med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt,1).put_data(g,sigt1_org*sigt0_new/sigt0_org);
	    };

	  };
        };
      };
    };
  };

  // +++ Macroscopic cross section is revised
  for(int g=init_grp;g<grp;g++){
    // (scattering matrix)
    for(int g2=init_grp;g2<grp;g2++){
      for(int l=0;l<=pl_med;l++){
        real tmp=0.;
        for(int i=0;i<med.GetNucnum();i++){
          int ng=med.GetNuclideInTurn(i).GetGrp();
	  if(ng!=-1){
  	    real den=med.GetNuclideInTurn(i).GetDensity();
	    for(int j=0;j<3;j++){
	      if(l<med.GetNuclideInTurn(i).GetMicxs().GetDim2d(j)){
 	        tmp+=den*med.GetNuclideInTurn(i).GetMicxs().GetData2d(j,l).get_dat(g,g2);
	      };
	    };
	  };
	};
	med.GetMacxs().GetData2d(sigs,l).put_data(g,g2,tmp);
      };
    };

    real sa=med.GetMacxs().GetData1d(siga).get_dat(g);
    real ss=0.;
    real ss1=0.;
    for(int g2=0;g2<grp;g2++){
      ss+=med.GetMacxs().GetData2d(sigs).get_dat(g,g2);
      if(pl_med>0)ss1+=med.GetMacxs().GetData2d(sigs,1).get_dat(g,g2);
    };
    if(pl_med==0){
      for(int jj=0;jj<med.GetNucnum();jj++){
        int ng=med.GetNuclideInTurn(jj).GetGrp();
	if(ng!=-1){
        ss1+=med.GetNuclideInTurn(jj).GetMicxs().GetData1d(mu).get_dat(g)
           *med.GetNuclideInTurn(jj).GetMicxs().GetData1d(sigel).get_dat(g)
           *med.GetNuclideInTurn(jj).GetDensity();
	};
      };
    };
    real s2n=med.GetMacxs().GetData1d(sign2n).get_dat(g);
    real st=sa+ss-s2n;
    med.GetMacxs().GetData1d(sigt).put_data(g,st);
    if(pl_med>0){
      med.GetMacxs().GetData1d(sigtr).put_data(g,st-ss1*0.33333333333333333);
    }else{
      med.GetMacxs().GetData1d(sigtr).put_data(g,st-ss1);
      med.CalDFromSigtr();
    };
    // with transport XS in low energy range is calculated by flux-weighted total
  };

  // Pn-total cross section is replaced by P0-total
  int plt=med.GetPLT();
  for(int i=init_grp;i<med.GetImax();i++){
    real tmp=med.GetMacxs().GetData1d(sigt).get_dat(i);
    for(int l=1;l<=plt;l++){
      med.GetMacxs().GetData1d(sigt,l).put_data(i,tmp);
    };
  };

};

void OnePointCalculator::CoarseGroupCorrectionForCurrentWeightCrossSection(Medium &inp)
{
  real den=inp.GetNuclide(928).GetDensity();
  real fact[]={3.74,3.14,5.68};
  int gg[]={50,52,56};
  for(int i=0;i<3;i++){
    real micst0=inp.GetNuclide(928).GetMicxs().GetData1d(sigt).get_dat(gg[i]);
    inp.GetMacxs().GetData1d(sigt).add_data(gg[i],den*(fact[i]-1)*micst0);
    inp.GetMacxs().GetData2d(sigs).add_data(gg[i],gg[i],den*(fact[i]-1)*micst0);
  };
};

void OnePointCalculator::CheckEnergyGroupConsistency(Medium &inp,XSLibrary &xslib)
{
  int grp_med=inp.GetImax();
  int grp_lib=xslib.GetGroup();
  if(grp_med!=grp_lib){
    cout<<"# Error in OnePointCalculator::CheckEnergyGroupConsistency.\n";
    cout<<"# Energy group is inconsistent between Medium and XSLibrary.\n";
    cout<<"#\n";
    cout<<"#   The number of energy groups of medium  : "<<grp_med<<"\n";
    cout<<"#   The number of energy groups of library : "<<grp_lib<<"\n";
    if(grp_med==0){
      cout<<"#\n";
      cout<<"# Please do the [PutImax] method for this medium before using.\n";
    };
    exit(0);
  };
};
