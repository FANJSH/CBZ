#include <cstdlib>
#include "SelfShieldingCalculator.h"

SelfShieldingCalculator::SelfShieldingCalculator()
{
  background_xs_printing=false;  
  cw_correction=true;
};

void SelfShieldingCalculator::PutPijCalculator(TrajectorySet *tset)
{
  pijcalculator_general=tset;
  pij_c=general;
  region=GetRegion();
  medium_id.resize(region,0);
};

void SelfShieldingCalculator::PutPijCalculator(PJISlabPij *pspij)
{
  pijcalculator_slab=pspij;
  pij_c=slab;
  region=GetRegion();
  medium_id.resize(region,0);
};

int SelfShieldingCalculator::GetRegion()
{
  if(pij_c==general){
    return pijcalculator_general->GetRegnum();
  }else{
    return pijcalculator_slab->GetMesh();
  };
  return 1;
};

real SelfShieldingCalculator::GetVolume(int i)
{
  if(pij_c==general){
    return pijcalculator_general->GetVol(i);
  }else{
    return pijcalculator_slab->GetWidth(i);
  };
  return 0.;
};

void SelfShieldingCalculator::CalculationPij(real *xs,real *pij)
{
  if(pij_c==general){
    pijcalculator_general->CalculationPij(xs,pij,false);
    //pijcalculator_general->CalculationPijSphere(xs,pij);
  }else{
    pijcalculator_slab->CalculationPij(xs,pij);
  };
};

void SelfShieldingCalculator::SetNumberOfMedium(int i)
{
  num_medium=i;
  medium.resize(num_medium);
};

void SelfShieldingCalculator::PutMedium(int number,Medium &medinp)
{
  if(number<0||number>=num_medium){
    cout<<"Error in SelfShieldingCalculator::PutMedium.\n";
    cout<<"Medium ID should be from 0 to "<<num_medium-1<<"\n";
    cout<<"You set medium ID as "<<number<<"\n";
    exit(0);
  };
  medium[number]=medinp;
};

void SelfShieldingCalculator::PutMediumID(int *reg_med)
{
  for(int i=0;i<region;i++){
    medium_id[i]=reg_med[i];
  };
};

void SelfShieldingCalculator::GiveInfiniteDilutionCrossSection(XSLibrary &xslib)
{
  OnePointCalculator opc;
  for(int i=0;i<num_medium;i++){
    medium[i].GetEnband().copy(xslib.GetEnband());
    opc.GiveInfiniteDillutionCrossSection(medium[i],xslib);
  };
};

void SelfShieldingCalculator::WithToneMethod(XSLibrary &xslib, bool print)
{
  GiveInfiniteDilutionCrossSection(xslib);

  int group=medium[0].GetImax();
  int reg=GetRegion();

  vector<real> volume(reg);
  for(int i=0;i<reg;i++){
    volume[i]=GetVolume(i);
  };

  int iter_max=0;
  for(int iter=0;iter<=iter_max;iter++){

    for(int i=0;i<group;i++){
      if(print)cout<<"#   Tone method in group "<<i<<"\n";
      real *xs=new real[reg];
      for(int j=0;j<reg;j++){
        xs[j]=medium[j].GetMacxs().GetData1d(sigt).get_dat(i);
      };
      real *pij=new real[reg*reg];
      CalculationPij(xs,pij);
      for(int j=0;j<reg;j++){
        int nuc=medium[j].GetNucnum();
        for(int k=0;k<nuc;k++){
	  real den=medium[j].GetNuclideInTurn(k).GetDensity();
	  int ng=medium[j].GetNuclideInTurn(k).GetGrp();
	  if(den>0.&&ng!=-1){
            int id=medium[j].GetNuclideID(k);
    	    int sttgrp=xslib.GetLibData(id).GetFtable().GetStartGrpAllMF();
	    int endgrp=xslib.GetLibData(id).GetFtable().GetEndGrpAllMF();
            if(i>=sttgrp&&i<=endgrp){
    	      real denom=0.;
	      real nume=0.;
	      real r_nume=0.;
	      real r_sigt=0.;
	      for(int l=0;l<reg;l++){
	        int nuc2=medium[l].GetNucnum();
  	        real vol2=volume[l];
	        for(int m=0;m<nuc2;m++){
	          int id2=medium[l].GetNuclideID(m);
        	  int ng2=medium[l].GetNuclideInTurn(m).GetGrp();
		  if(ng2!=-1){
  	            real den2=medium[l].GetNuclideInTurn(m).GetDensity();
                    if(id==id2){
                      denom+=vol2*den2*pij[l*reg+j];
	            }else if((id==922380&&id2==942390)||id2==922380){
	              r_nume+=vol2*den2*pij[l*reg+j];
  	              real micsigt=medium[l].GetNuclideInTurn(m).GetMicxs().GetData1d(sigt).get_dat(i);
  	              r_sigt+=vol2*den2*micsigt*pij[l*reg+j];
	            }else{
  	              real micsigt=medium[l].GetNuclideInTurn(m).GetMicxs().GetData1d(sigt).get_dat(i);
  	              nume+=vol2*den2*micsigt*pij[l*reg+j];
		    };
	          };
	        };
	      };
	      real denom_inv=1./denom;
              real sig0=nume*denom_inv;
	      real rval=r_nume*denom_inv;
  	      real rsigt=r_sigt*denom_inv;
	      real temp=medium[j].GetNuclideInTurn(k).GetTemperature();
	      bool matrix_cal=false;
	      if(iter==iter_max)matrix_cal=true;
	      vector<real> f(6);
 	      xslib.GetLibData(id).CalF(i,0,rval,rsigt,temp,sig0,f);

#if 0
	      cout.setf(ios::scientific);
	      cout.precision(6);
              if(i==279){
		cout<<id<<" "<<j<<" "<<" "<<sig0<<" "<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<f[4]<<" "<<f[5]<<"\n";
	      };
#endif
	      
#if 0
	      if(j==0){
		real siginf=xslib.GetLibData(id).GetXSData().GetData1d(sigf).get_dat(i);
		if(siginf!=0.){
		  cout.setf(ios::scientific);
		  cout.precision(6);
  	  	  cout<<iter<<" "<<i<<" "<<id<<" "<<siginf<<" "<<f[0]<<" "<<f[0]*siginf<<"\n";
		};
	      };
#endif	      

              medium[j].GetNuclideInTurn(k).CalMicroFromF(i,xslib.GetLibData(id),f,matrix_cal);
	      

	    };

	  };

        };
      };
      delete [] xs;
      delete [] pij;
    };

    for(int j=0;j<reg;j++){
      medium[j].CalMacroFromMicro();
    };
    
  };

  OnePointCalculator opc;
  for(int i=0;i<reg;i++){
    if(medium[i].GetMacxs().GetData1d(nusigf).get_dat(0)>0.){
      opc.CalIncidentEnergyDependentChi(medium[i],xslib);
    };
  };

};

void SelfShieldingCalculator::CalMixture()
{
  vector<int> ihnid;
  int ihnuc=0;
  for(int i=0;i<num_medium;i++){
    for(int j=0;j<medium[i].GetNucnum();j++){
      if(ihnuc==0){
	ihnuc++;
	ihnid.push_back(medium[i].GetNuclideID(j));
      }else{
	int itmp=0;
	for(int k=0;k<ihnuc;k++){
	  if(ihnid[k]==medium[i].GetNuclideID(j))itmp=1;
	};
	if(itmp==0){
	  ihnuc++;
	  ihnid.push_back(medium[i].GetNuclideID(j));
	};
      };
    };
  };

  vector<real> denh;
  denh.resize(ihnuc,0.);
  for(int i=0;i<num_medium;i++){
    for(int j=0;j<medium[i].GetNucnum();j++){
      int id=medium[i].GetNuclideID(j);
      real den=medium[i].GetNuclideInTurn(j).GetDensity();
      real vol=GetVolume(i);
      for(int k=0;k<ihnuc;k++){
	if(id==ihnid[k]){
	  denh[k]+=den*vol;
	};
      };
    };
  };

  real total_vol=0.;
  for(int i=0;i<num_medium;i++){
    total_vol+=GetVolume(i);
  };
  total_vol=1./total_vol;

  for(int i=0;i<ihnuc;i++){
    denh[i]*=total_vol;
    cout<<ihnid[i]<<" : "<<denh[i]<<"\n";
  };
};

void SelfShieldingCalculator::DancoffMethod
(XSLibrary &xslib, TrajectorySet &tset, vector<Medium> &med, real eng)
{
  int mednum=med.size();
  if(mednum==3){
    ThreeRegionDancoffMethod(xslib,tset,med[0],med[1],med[2],true,eng);
  }else if(mednum==4){
    FourRegionDancoffMethod(xslib,tset,med[0],med[1],med[2],med[3],true,eng);
  }else{
    cout<<"# Error in SelfShieldingCalculator::DancoffMethod.\n";
    cout<<"# Three- or four-medium system only can be treated.\n";
    exit(0);
  };
  return;
};

void SelfShieldingCalculator::ThreeRegionDancoffMethod
(XSLibrary &xslib, TrajectorySet &sys,
 Medium &med0, Medium &med1, Medium &med2, bool clad,real eng)
{
  OnePointCalculator opc;
  if(background_xs_printing)opc.BackgroundXSPrinting();
  if(!cw_correction)opc.CWCorrectionOff();
  
  int group=med0.GetImax();

  int gg=0;
  if(eng>0.){
    real et=xslib.GetEnband().get_dat(0);
    for(int i=0;i<group;i++){
      real el=xslib.GetEnband().get_dat(i+1);
      if(eng<=et&&eng>el)gg=i;
      et=el;
    };
    cout<<"\n# Dancoff factor is calculated in group "<<gg;
    cout<<" (energy : "<<eng<<" eV)\n\n";
  };

  real vol1=sys.GetVolume(0);
  real vol2=sys.GetVolume(1);
  vol2+=vol1;

  real r=sqrt(vol1/PI);
  real r2=sqrt(vol2/PI);

  GroupData1D c(group);
  GroupData1D b(group);
  real xs[3];
  real pij[3*3];

  int matmax=1;
  if(clad)matmax=2;
  dancoff.resize(matmax);
  for(int mat=0;mat<matmax;mat++){
    if(mat==1)r=r2-r;
    if(eng<0.){
      // (energy group-dependent Dancoff)
      for(int i=0;i<group;i++){
        xs[2]=med2.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[1]=med1.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[0]=med0.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[mat]=1e5;
        sys.CalculationPij(xs,pij,false);
        real dancoff_corr=1.-xs[mat]*2.*r*(1.-pij[mat*3+mat]);
        c.put_data(i,dancoff_corr);
        b.put_data(i,1.2); // bell factor
      };
    }else{
      // (energy group-independent Dancoff)
      xs[2]=med2.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[1]=med1.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[0]=med0.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[mat]=1e5;
      sys.CalculationPij(xs,pij,false);
      real dancoff_corr=1.-xs[mat]*2.*r*(1.-pij[mat*3+mat]);
      for(int i=0;i<group;i++){
        c.put_data(i,dancoff_corr);
        b.put_data(i,1.2); // bell factor
      };
    };
    if(mat==0)opc.CalSelfShieldingWithDancoffCorrection(med0,xslib,r*2.,b,c);
    if(mat==1)opc.CalSelfShieldingWithDancoffCorrection(med1,xslib,r*2.,b,c);
    dancoff[mat]=c;
  };

  /*
  cout.setf(ios::scientific);
  cout.precision(6);
  for(int g=0;g<group;g++){
    real etop=med1.GetEnband().get_dat(g);
    real sig=med1.GetNuclide(400000).GetMicxs().GetData1d(sigc).get_dat(g);
    cout<<etop<<" "<<sig<<"\n";
  };
  exit(0);
  */
};

void SelfShieldingCalculator::FourRegionDancoffMethod
(XSLibrary &xslib, TrajectorySet &sys,
 Medium &med0, Medium &med1, Medium &med2, Medium &med3, bool clad,real eng)
{
  OnePointCalculator opc;
  if(background_xs_printing)opc.BackgroundXSPrinting();
  if(!cw_correction)opc.CWCorrectionOff();

  int group=med0.GetImax();

  int gg=0;
  if(eng>0.){
    real et=xslib.GetEnband().get_dat(0);
    for(int i=0;i<group;i++){
      real el=xslib.GetEnband().get_dat(i+1);
      if(eng<=et&&eng>el)gg=i;
      et=el;
    };
    cout<<"\n# Dancoff factor is calculated in group "<<gg;
    cout<<" (energy : "<<eng<<" eV)\n\n";
  };

  real vol1=sys.GetVolume(0);
  real vol2=sys.GetVolume(1);
  vol2+=vol1;

  real r=sqrt(vol1/PI);
  real r2=sqrt(vol2/PI);

  GroupData1D c(group);
  GroupData1D b(group);
  real xs[4];
  real pij[4*4];

  int matmax=1;
  if(clad)matmax=3;
  dancoff.resize(matmax);
  for(int mat=0;mat<matmax;mat++){
    if(mat==1)r=r2-r;
    if(eng<0.){
      // (energy group-dependent Dancoff)
      for(int i=0;i<group;i++){
        xs[3]=med3.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[2]=med2.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[1]=med1.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[0]=med0.GetMacxs().GetData1d(sigt).get_dat(i);
        xs[mat]=1e5;
        sys.CalculationPij(xs,pij,false);
        real dancoff_corr=1.-xs[mat]*2.*r*(1.-pij[mat*4+mat]);
        c.put_data(i,dancoff_corr);
        b.put_data(i,1.2); // bell factor
      };
    }else{
      // (energy group-independent Dancoff)
      xs[3]=med3.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[2]=med2.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[1]=med1.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[0]=med0.GetMacxs().GetData1d(sigt).get_dat(gg);
      xs[mat]=1e5;
      sys.CalculationPij(xs,pij,false);
      real dancoff_corr=1.-xs[mat]*2.*r*(1.-pij[mat*4+mat]);
      for(int i=0;i<group;i++){
        c.put_data(i,dancoff_corr);
        b.put_data(i,1.2); // bell factor
      };
    };
    if(mat==0)opc.CalSelfShieldingWithDancoffCorrection(med0,xslib,r*2.,b,c);
    if(mat==1)opc.CalSelfShieldingWithDancoffCorrection(med1,xslib,r*2.,b,c);
    if(mat==2)opc.CalSelfShieldingWithDancoffCorrection(med2,xslib,r*2.,b,c);    
    dancoff[mat]=c;
  };
};

/*
void SelfShieldingCalculator::WithDancoffMethod(XSLibrary &xslib,TrajectorySet &sys)
{
  OnePointCalculator opc;

  int number_medium=medium.size();
  int group=medium[0].GetImax();
  int region=sys.GetRegnum();

  for(int i=0;i<number_medium;i++){
    opc.GiveInfiniteDillutionCrossSection(medium[i],xslib);
  };

  GroupData1D c(group);
  GroupData1D b(group);

  real *xs=new real[region];
  for(int i=0;i<group;i++){
    for(int j=0;j<region;j++){
      xs[j]=medium[medium_id[j]].GetMacxs().GetData1d(sigt).get_dat(i);
      if(medium_id[j]==0)xs[j]=30000.;
    };
    real *pij=new real[region*region];
    sys.CalculationPij(xs,pij,false);
    real pff=0.;
    for(int j=0;j<region;j++){
      if(medium_id[j]==0)pff+=pij[0*region+j];
    };
    real pesc=1.-pff;
    real dancoff_corr=1.-xs[0]*2.*r*pesc;
    c.put_data(i,dancoff_corr);
    if(jfs3j32r.GetEnband().get_dat(i)<100.){
      b.put_data(i,1.1); // bell factor
    }else{
      b.put_data(i,1.2);
    };
  };
  opc.CalSelfShieldingWithDancoffCorrection(medium[0],jfs3j32r,r*2.,b,c);
  delete [] xs;
};
*/



/*
void SelfShieldingCalculator::CalDancoffCorrection(real average_chord)
{
  GiveInfiniteDilutionCrossSection(xslib);

  real black_judge=0.01;

  real e=10.;
  int group=0;
  for(int i=0;i<medium[0].GetImax();i++){
    real ee=medium[0].GetEnband().get_dat(i+1);
    if(ee<e){
      group=i;
      break;
    };
  };

  real *xs=new real[region];
  real *pij=new real[region*region];
  vector< vector<real> > dancoff(region);

  for(int i=0;i<region;i++){
    int mid=medium_id[i];
    int nuc=medium[mid].GetNucnum();
    dancoff[i].resize(nuc,0.);
    for(int j=0;j<nuc;j++){
      int nucid=medium[mid].GetNuclideID(j);
      for(int k=0;k<region;k++){
	xs[k]=medium[medium_id[k]].GetMacxs().GetData1d(sigt).get_dat(group);
	if(medium[medium_id[k]].ExistNuclide(nucid)){
	  if(average_chord*medium[medium_id[k]].GetNuclide(nucid).GetMicxs().GetData1d(sigt).get_dat(group)>black_judge){
	    xs[k]=30000.;
	  };
	};
      };
      CalculationPij(xs,pij);
      real pff=0.;
      for(int k=0;k<region;k++){
        if(xs[k]==30000.)pff+=pij[i*region+k];
      };
      real pesc=1.-pff;
      real dancoff_corr=1.-xs[i]*average_chord*pesc;
      dancoff[i][j]=dancoff_corr;
      cout<<i<<" "<<nucid<<" "<<dancoff_corr<<"\n";
    };
  };

  delete [] xs;
  delete [] pij;
};

void SelfShieldingCalculator::CalDancoffCorrection(TrajectorySet &tset,int group,int nucid,real average_chord)
{
  //int number_medium=medium.size();
  int region=tset.GetRegnum();
  real black_judge=0.01;
  dancoff.resize(region,0.);

  real *xs=new real[region];
  vector<bool> black_region(region);

  for(int i=0;i<region;i++){
    xs[i]=medium[i].GetMacxs().GetData1d(sigt).get_dat(group);
    black_region[i]=false;
    if(medium[i].ExistNuclide(nucid)){
      real den=medium[i].GetNuclide(nucid).GetDensity();
      real micsigt=medium[i].GetNuclide(nucid).GetMicxs().GetData1d(sigt).get_dat(group);
      if(den*micsigt*average_chord>black_judge){
        xs[i]=30000.;
	black_region[i]=true;
      };
    };
  };

  real *pij=new real[region*region];
  tset.CalculationPij(xs,pij,false);

  for(int i=0;i<region;i++){
    if(black_region[i]){
      real pff=0.;
      for(int j=0;j<region;j++){
        if(black_region[j])pff+=pij[i*region+j];
      };
      real pesc=1.-pff;
      real dancoff_corr=1.-xs[i]*average_chord*pesc;
      cout<<dancoff_corr<<"\n";
      dancoff[i]=dancoff_corr;
    };
  };

  delete [] xs;
  delete [] pij;
};

void SelfShieldingCalculator::WriteDancoffCorrection()
{
  int reg=dancoff.size();
  for(int i=0;i<reg;i++){
    cout<<i<<" "<<dancoff[i]<<"\n";
  };
};

void SelfShieldingCalculator::WithDancoffMethod()
{
  if(solver=="tone"){


  //pin[1].CalMacroFromMicro();
  //pin[2].CalMacroFromMicro();

  };
  };


  //PJISystem lat(70,16);
  //lat.PutTrajectorySet(&sys);
  //for(int i=0;i<region;i++){
  //  lat.AddMedium(pin[reg_med[i]]);
  //};
  //int region_medium[]={0,1,2};
  //lat.PutRegMed(region_medium);
  //GeneralOption opt;
  //lat.PutGeneralOption(opt);
  //lat.PutPij(sigt);
  //real k1=lat.CalIgenPij();
*/
