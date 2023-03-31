#include <cstdlib>
#include "Medium.h"

Medium::~Medium()
{
  MacxsVectorClear();
  NuclideClear();
};

void Medium::Init()
{
  pl   =-1;
  pltot=-1;
  imax = 0;
  nucnum=0;
  flux.resize(2);
  fissionable_flag=0;
};

void Medium::PutPL(int i,int itot)
{
  pl=i;
  if(itot==-1){itot=pl;};
  if(itot==0)itot=1;
  pltot=itot;
  InitSigST();
};

void Medium::PutImax(int i)
{
  imax=i;
  flux[0].put_imax(i);
  flux[1].put_imax(i);
  enband.put_imax(i+1);
  LowestDownScatGrp.resize(i,0);
  InitSigST();
}

void Medium::InitSigST()
{
  if(imax!=0&&pl!=-1&&pltot!=-1){
    Macxs.Init("MacroCrossSection");
    Macxs.PutDim1d(sigt,pltot+1);
    Macxs.PutDim2d(sigs,pl+1);
    Macxs.PutGrp(imax);
  };
}

void Medium::ReadPDSFile(string mdir,string ss,int plt,bool upscat)
{
  //cout<<"Reading PDS File:"<<ss<<"\n";
  mdir.append(ss);

  if(plt==0){
    cout<<"Pl order for sigt is reset to 1.\n";
    plt=1;
  };

  ifstream fin;
  ofstream fout;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  int ig;
  fin>>ig;
  int ipl;
  fin>>ipl;
  PutPL(ipl,plt);
  PutImax(ig);

  real inp;
  for(int i=0;i<ig+1;i++){
    fin>>inp;
    enband.put_data(i,inp);
  };
  
  enum xstype inpname[]={nusigf,siga,chi,sign2n,sigt};
  for(int in=0;in<5;in++){
    //if(in!=3){
      for(int i=0;i<ig;i++){
        fin>>inp;
        Macxs.GetData1d(inpname[in]).put_data(i,inp);
      };
      //}else{
      //for(int i=0;i<ig;i++){
      //fin>>inp;
      //};
      //};
  };

  for(int i=0;i<ig;i++){
    fin>>inp;
    for(int l=1;l<=pltot;l++){
      Macxs.GetData1d(sigt,l).put_data(i,inp);
    };
  };

  enum xstype inpname2[]={d,dr,dz}; // dr:perp, dz:para
  for(int in=0;in<3;in++){
    for(int i=0;i<ig;i++){
      fin>>inp;
      Macxs.GetData1d(inpname2[in]).put_data(i,inp);
    };
  };

  for(int l=0;l<=pl;l++){
    int ind=0;
    for(int i=0;i<ig;i++){
      for(int j=0;j<ig;j++){
	if(!upscat&&i>j){inp=0.;}
	else{fin>>inp;}
	Macxs.GetData2d(sigs,l).put_data(ind,inp);
	ind++;
      };
    };
  };

  // Read microscopic cross section
  int ninp;
  fin>>ninp;
  PutNucnum(ninp);
  if(nucnum!=0){
    int imtx;
    fin>>imtx; // 0:no matrix   1:matrix
    for(int i=0;i<nucnum;i++){
      nuc[i].PutGrp(ig);
      int nucid;
      fin>>nucid;
      nuc[i].PutMatnum(nucid);
      real den;
      fin>>den;
      nuc[i].PutDensity(den);
      real tmp;
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigf).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigc).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
	real tmpinp=0.;
	real sf=nuc[i].GetMicxs().GetData1d(sigf).get_dat(g);
	if(sf>0.)tmpinp=tmp/sf;
        nuc[i].GetMicxs().GetData1d(nu).put_data(g,tmpinp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(chi).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(mu).put_data(g,tmp);
      };
      if(imtx==1){
	enum xstype scname[]={sigel,siginel,sign2n};
	for(int l=0;l<2;l++){
  	  for(int j=0;j<3;j++){
            int ind=0;
            for(int g=0;g<ig;g++){
	      for(int gg=0;gg<ig;gg++){
	        if(gg>=g){
  	          fin>>tmp;
	        }else{
                  tmp=0.;
	        };
                nuc[i].GetMicxs().GetData2d(scname[j],l).put_data(ind,tmp);
	        ind++;
	      };
	    };
	  };
	};
	for(int g=0;g<ig;g++){
	  real tmp=0.;
	  for(int gg=g;gg<ig;gg++){
	    tmp+=nuc[i].GetMicxs().GetData2d(sigel).get_dat(g,gg);
	  };
	  nuc[i].GetMicxs().GetData1d(sigel).put_data(g,tmp);
	};
      };
    };
  };

  fin.close();
  CalSigtr();

  //cout<<"  ...Finished.\n";
}

void Medium::ReadFile(string mdir,string ss,int plt)
{
  // [plt]
  //
  //  - Maximum Legendre order for energy-averaged total cross section definition.
  //  - Generally this is set to 1: flux- and current-weighted total cross sections are used.

  /*
  cout<<"#################################################\n";
  cout<<"#\n";                                     
  cout<<"# WARNING IN Medium::ReadFile.\n";
  cout<<"#\n";
  cout<<"# Flux/current data are added in 2021/6/27,\n";
  cout<<"# so previously-created files cannot be properly read now.\n";
  cout<<"#\n";
  cout<<"#################################################\n";  
  */
  
  //cout<<"# Reading Medium data from File : "<<ss<<"\n";
  mdir.append(ss);

  if(plt==0){
    cout<<"#   Pl order for sigt is reset to 1.\n";
    plt=1;
  };

  ifstream fin;
  ofstream fout;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# "<<mdir.data()<<"\n";
    exit(1);
  };

  int ig;  
  fin>>ig;  // The number of energy groups
  int ipl;
  fin>>ipl; // Legendre order to represent scattering cross section
  PutPL(ipl,plt);
  PutImax(ig);

  real inp;
  for(int i=0;i<ig+1;i++){
    fin>>inp;
    enband.put_data(i,inp); // energy boundaries in this group structure
  };

  // 1D macroscopic cross section: Nu-Sigma_f, Sigma_a, Fission spectra, (n,2n), Sigma_t  
  enum xstype inpname[]={nusigf,siga,chi,sign2n,sigt};
  for(int in=0;in<5;in++){
    for(int i=0;i<ig;i++){
      fin>>inp;
      Macxs.GetData1d(inpname[in]).put_data(i,inp);
    };
  };

  // 1D macroscopic total cross section weighted by high-order angular moment of neutron flux
  for(int i=0;i<ig;i++){
    fin>>inp;
    for(int l=1;l<=pltot;l++){
      Macxs.GetData1d(sigt,l).put_data(i,inp);
    };
  };

  // Average and anisotropic diffusion coefficient
  enum xstype inpname2[]={d,dr,dz}; // dr:perp, dz:para
  for(int in=0;in<3;in++){
    for(int i=0;i<ig;i++){
      fin>>inp;
      Macxs.GetData1d(inpname2[in]).put_data(i,inp);
    };
  };

  // Flux & Current
  for(int m=0;m<2;m++){
    for(int i=0;i<ig;i++){
      fin>>inp; 
      flux[m].put_data(i,inp);
    };
  };
  
  // Macroscopic scattering cross section in a matrix form
  for(int l=0;l<=pl;l++){
    for(int i=0;i<ig;i++){
      int tmp;
      fin>>tmp;
      int tmp2;
      fin>>tmp2;
      if(tmp!=-1){
        for(int j=tmp;j<=tmp2;j++){
	  fin>>inp;
	  Macxs.GetData2d(sigs,l).put_data(i,j,inp);
        };
      };
    };
  };

  // Read microscopic cross section
  int ninp;
  fin>>ninp;
  PutNucnum(ninp);
  if(nucnum!=0){
    int imtx;
    fin>>imtx; // 0:no matrix   1:matrix
    for(int i=0;i<nucnum;i++){

      int nucid;
      fin>>nucid;
      nuc[i].PutMatnum(nucid);
      real den;
      fin>>den;
      nuc[i].PutDensity(den);
      int grpinp;
      fin>>grpinp;

      if(grpinp!=-1){

      nuc[i].PutGrp(grpinp);
      real tmp;
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigf).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigc).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
	real tmpinp=0.;
	real sf=nuc[i].GetMicxs().GetData1d(sigf).get_dat(g);
	if(sf>0.)tmpinp=tmp/sf;
        nuc[i].GetMicxs().GetData1d(nu).put_data(g,tmpinp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sign2n).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigel).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(siginel).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(sigt).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(chi).put_data(g,tmp);
      };
      for(int g=0;g<ig;g++){
        fin>>tmp;
        nuc[i].GetMicxs().GetData1d(mu).put_data(g,tmp);
      };
      if(imtx==1){
	enum xstype scname[]={sigel,siginel,sign2n};
	for(int l=0;l<2;l++){
  	  for(int j=0;j<3;j++){
	    if(l==1&&j==2){
	      nuc[i].GetMicxs().GetData2d(scname[j],l).set_zero();
	    }else{
	      for(int g=0;g<imax;g++){
		int stgrp,edgrp;
		fin>>stgrp;
		fin>>edgrp;
		for(int gg=stgrp;gg<=edgrp;gg++){
		  real tmp;
		  fin>>tmp;
                  nuc[i].GetMicxs().GetData2d(scname[j],l).put_data(g,gg,tmp);
		};
	      };
	    };
	  };
	};
      };
      };

    };
  };

  fin.close();
  CalSigtr();

  //cout<<"#  ...Finished.\n";
}

void Medium::ReadKramxsFile(string mdir,string ss,bool up_scat)
{
  // [plt]
  //
  //  - Maximum Legendre order for energy-averaged total cross section definition.
  //  - Generally this is set to 1: flux- and current-weighted total cross sections are used.

  if(imax==0){
    cout<<"# Error in Medium::ReadKramxsFile.\n";
    cout<<"# The number of enregy groups is not yet set.\n";
    cout<<"# Please do [PutImax] before this.\n";
    exit(0);
  };
  if(pl==-1){
    cout<<"# Error in Medium::ReadKramxsFile.\n";
    cout<<"# The number of Pl truncation order is not yet set.\n";
    cout<<"# Please do [PutPL] before this.\n";    
    exit(0);
  };
  
  cout<<"# Reading Medium data from file in the KRAMXS format named as "<<ss<<"\n";
  mdir.append(ss);

  ifstream fin;
  ofstream fout;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    exit(1);
  };
  
  string dummy;
  fin>>dummy;
  fin>>dummy;

  for(int i=0;i<imax;i++){
    for(int j=0;j<imax;j++){
      real tmp;
      fin>>tmp;
      if(up_scat){
        Macxs.GetData2d(sigs).put_data(j,i,tmp);
      }else{
      if(i<j){
        if(tmp>0.)cout<<"# Warning : up-scattering cross section "<<tmp<<" from "<<j<<" to "<<i<<" is NOT considered.\n";
      }else{
        Macxs.GetData2d(sigs).put_data(j,i,tmp);	
      };
      };
      //Macxs.GetData2d(sigs).put_data(j,i,tmp);	     
    };
  };

  enum xstype inpname[]={nusigf,siga,sigt,d};
  for(int in=0;in<4;in++){
    int l=0;
    if(in==2)l=1;
    for(int i=0;i<imax;i++){
      real tmp;
      fin>>tmp;
      Macxs.GetData1d(inpname[in],l).put_data(i,tmp);      
    };
  };

  Macxs.GetData1d(sigt)=Macxs.GetData1d(siga)+Macxs.GetData2d(sigs).get_sumx();  
  
  fin>>dummy;
  fin>>dummy;  
  for(int i=0;i<imax;i++){
    real tmp;
    fin>>tmp;
    Macxs.GetData1d(chi).put_data(i,tmp);
  };

  fin>>dummy;
  fin>>dummy;

  int order;
  fin>>order; // P_l order
  for(int l=0;l<order;l++){
    for(int i=0;i<imax;i++){
      for(int j=0;j<imax;j++){
        real tmp;
        fin>>tmp;
        if(l+1<=pl){
	  if(up_scat){
  	    Macxs.GetData2d(sigs,l+1).put_data(j,i,tmp*(2*l+3));
	  }else{
	  if(i<j){
	    //cout<<"# Warning : up-scattering cross section "<<tmp<<" is NOT considered.\n";
	  }else{
  	    Macxs.GetData2d(sigs,l+1).put_data(j,i,tmp*(2*l+3));
	  };
	  };
	};
      };
    };
  };

  //Macxs.GetData1d(sigt,0)=Macxs.GetData1d(sigt,1);
  
  fin.close();
  CalSigtr();
  CalDFromSigtr();

  cout<<"#  ...Finished.\n";
}

void Medium::WritePDSFile(string mdir,string ss)
{
  cout<<"# Writing Medium data on file : "<<ss<<"\n";
  mdir.append(ss);

  ifstream fin;
  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };
  fout<<imax<<"\n";
  fout<<pl<<"\n";
  int ig=imax;

  for(int i=0;i<ig+1;i++){
    fout<<enband.get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(nusigf).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(siga).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(chi).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sign2n).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sigt).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sigt,1).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(d).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(dr).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(dz).get_dat(i)<<"\n";
  };
  for(int j=0;j<=pl;j++){
    for(int i=0;i<ig;i++){
      for(int k=i;k<ig;k++){
	fout<<Macxs.GetData2d(sigs,j).get_dat(i,k)<<"\n";
      };
    };
  };

  // nucnum=0; // temp
  fout<<0;

  fout.close();
}

void Medium::WriteFile(string mdir,string ss,bool mic)
{
  cout<<"# Writing Medium data on file : "<<ss<<"\n";
  mdir.append(ss);

  ifstream fin;
  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };
  fout<<imax<<"\n";
  fout<<pl<<"\n";
  int ig=imax;

  fout.setf(ios::scientific);
  fout.precision(6);
  for(int i=0;i<ig+1;i++){
    fout<<enband.get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(nusigf).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(siga).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(chi).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sign2n).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sigt).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(sigt,1).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(d).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(dr).get_dat(i)<<"\n";
  };
  for(int i=0;i<ig;i++){
    fout<<Macxs.GetData1d(dz).get_dat(i)<<"\n";
  };
  for(int m=0;m<2;m++){
    for(int i=0;i<ig;i++){
      fout<<flux[m].get_dat(i)<<"\n";
    };
  };

  
  for(int j=0;j<=pl;j++){
    for(int i=0;i<ig;i++){
      int ist=-1;
      int ied=-1;
      for(int k=0;k<ig;k++){
        real tmp=Macxs.GetData2d(sigs,j).get_dat(i,k);
        if(tmp!=0.){
          if(ist==-1)ist=k;
	  if(ied!=-1)ied=-1;
	}else{
  	  if(ied==-1)ied=k;
	};
      };
      if(ied==-1)ied=ig;
      ied--;
      fout<<ist<<"\n";
      fout<<ied<<"\n";
      if(ist!=-1){
        for(int k=ist;k<=ied;k++){
	  fout<<Macxs.GetData2d(sigs,j).get_dat(i,k)<<"\n";
	};
      };
    };
  };

  /*
  int numtmp=0;
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetMatnum()>290000&&nuc[i].GetMatnum()<300000)numtmp++;
  };
  */
      

  // nucnum=0; // temp
  if(!mic){
    fout<<0;
  }else{
    //fout<<"  "<<numtmp<<"  \n";
    fout<<"  "<<nucnum<<"  \n";
    fout<<" 1\n"; // with matrix
    for(int i=0;i<nucnum;i++){

      //if(nuc[i].GetMatnum()>290000&&nuc[i].GetMatnum()<300000){

      fout<<"  "<<nuc[i].GetMatnum()<<"  \n";
      fout<<"  "<<nuc[i].GetDensity()<<"\n";
      fout<<"  "<<nuc[i].GetGrp()<<"\n";
      if(nuc[i].GetGrp()!=-1){

      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(sigf).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(sigc).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
        real sf=nuc[i].GetMicxs().GetData1d(sigf).get_dat(g);
	if(sf>0.){
    	  fout<<sf*nuc[i].GetMicxs().GetData1d(nu).get_dat(g)<<"\n";
	}else{
	  fout<<" 0.\n";
	};
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(sign2n).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(sigel).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(siginel).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(sigt).get_dat(g)<<"\n";
      };
      for(int g=0;g<imax;g++){
	fout<<nuc[i].GetMicxs().GetData1d(chi).get_dat(g)<<"\n";
	//fout<<" 0\n"; // dummy for chi
      };
      for(int g=0;g<imax;g++){
	fout<<" 0\n"; // dummy for mu
      };
      enum xstype scname[]={sigel,siginel,sign2n};
      for(int l=0;l<2;l++){
	for(int j=0;j<3;j++){
	  if(l!=1||j!=2){ // P1 of (n,2n) is ignored.
  	    for(int g=0;g<imax;g++){
              int stgrp=0;
              int edgrp=imax-1;
	      for(int gg=0;gg<imax;gg++){
                real xs=nuc[i].GetMicxs().GetData2d(scname[j],l).get_dat(g,gg);              
  	        if(stgrp==gg&&fabs(xs)<1e-20)stgrp=gg+1;
                int gg2=imax-1-gg;
                real xs2=nuc[i].GetMicxs().GetData2d(scname[j],l).get_dat(g,gg2);              
		if(edgrp==gg2&&fabs(xs2)<1e-20)edgrp=gg2-1;
	      };
	      fout<<"   "<<stgrp<<"\n";
	      fout<<"   "<<edgrp<<"\n";
  	      for(int gg=stgrp;gg<=edgrp;gg++){
  	        fout<<nuc[i].GetMicxs().GetData2d(scname[j],l).get_dat(g,gg)<<"\n";              
	      };
	    };
	  };
	};
      };

      };      

    };
  };

  fout.close();
}

void Medium::ReadFileNumberDensity(string mdir,string ss,bool iso)
{
  cout<<"# Reading Medium number density data from File : "<<ss<<"\n";
  mdir.append(ss);

  ifstream fin;

  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  int tmp;
  fin>>tmp;
  PutImax(tmp);
  int ninp;
  fin>>ninp;

  if(!iso){

  PutNucnum(ninp);
  for(int i=0;i<ninp;i++){
    int nucid;
    fin>>nucid;
    real den;
    fin>>den;
    real temp;
    fin>>temp;
    int tmp;
    fin>>tmp;
    Nuclide nuc;
    if(tmp!=-1)nuc.PutGrp(tmp);
    nuc.PutMatnum(nucid);
    nuc.PutDensity(den);
    nuc.PutTemperature(temp);
    PutNuclide(i,nuc);
  };

  }else{

    int *matno_in=new int[ninp];
    real *den_in=new real[ninp];

    real temp=0.;
    for(int i=0;i<ninp;i++){
      fin>>matno_in[i];
      fin>>den_in[i];
      real temp_tmp;
      fin>>temp_tmp;
      if(i==0){
        temp=temp_tmp;
      }else{
	if(temp_tmp!=temp){
	  cout<<"# Error in Medium::ReadFileNumberDensity.\n";
	  cout<<"# Nuclide-wise temperature cannot be assigned.\n";
	  exit(0);
	};
      };
      int tmp;
      fin>>tmp;
    };

    PutNuclideNew(ninp,matno_in,den_in,true);
    PutTemperatureForAllNuclide(temp);
    delete [] matno_in;
    delete [] den_in;

  };

  fin.close();
}

void Medium::WriteFileNumberDensity(string mdir,string ss)
{
  cout<<"# Writing Medium number density data on file : "<<ss<<"\n";
  mdir.append(ss);

  ofstream fout;

  fout.open(mdir.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  fout<<"  "<<imax<<"\n";     // The number of energy groups of this medium
  fout<<"  "<<nucnum<<"  \n"; // The number of nuclides contained in this medium
  fout.setf(ios::scientific);
  fout.precision(9);
  for(int i=0;i<nucnum;i++){
    fout<<"  "<<nuc[i].GetMatnum()<<"  \n"; // Material ID of the i-th nuclide
    fout<<"  "<<nuc[i].GetDensity()<<"\n";  // Number density of the i-th nuclide
    fout<<"  "<<nuc[i].GetTemperature()<<"\n";  // Temperature of the i-th nuclide
    fout<<"  "<<nuc[i].GetGrp()<<"\n"; 
    // Number of energy groups of the i-th nuclide:
    // If this is zero, this nuclide has no cross section data.
  };
  fout.close();
}

void Medium::PutNucnum(int i)
{
  nucnum=i;
  if(nucnum!=0){
    nuc.resize(nucnum);
  };
};

void Medium::PutNuclide(string file)
{
  ifstream fin;
  fin.open(file.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  int nucnum;
  fin>>nucnum;
  int *mat=new int[nucnum];
  real *den=new real[nucnum];
  for(int j=0;j<nucnum;j++){
    fin>>mat[j];
    fin>>den[j];
  };
  PutNuclide(nucnum,mat,den);
  
  delete [] mat;
  delete [] den;
};

void Medium::PutNuclide(int nuc,int *matno)
{
  PutNucnum(nuc);
  for(int i=0;i<nucnum;i++){
    Nuclide tmp;
    tmp.PutMatnum(matno[i]);
    tmp.PutDensity(0.);
    PutNuclide(i,tmp);
  };
};

void Medium::PutNuclide(int nuc,vector<int> matno)
{
  PutNucnum(nuc);
  for(int i=0;i<nucnum;i++){
    Nuclide tmp;
    tmp.PutMatnum(matno[i]);
    tmp.PutDensity(0.);
    PutNuclide(i,tmp);
  };
};

void Medium::AddNuclide(int nuc,int *matno)
{
  nucnum+=nuc;
  for(int i=0;i<nuc;i++){
    Nuclide tmp;
    tmp.PutMatnum(matno[i]);
    tmp.PutDensity(0.);
    AddNuclide(tmp);
  };
};

void Medium::PutNuclide(int i,int mat,real den)
{
  nuc[i].PutMatnum(mat);
  nuc[i].PutDensity(den);
};

void Medium::PutNuclide(int nuc,int *matno,real *density,MATIDTranslator &midt)
{
  int *matno2=new int[nuc];
  for(int i=0;i<nuc;i++){
    matno2[i]=midt.GetMATIDFromENDFID(matno[i]);
  };
  PutNuclide(nuc,matno2,density);
  delete [] matno2;
};

void Medium::PutNuclide(int nuc,int *matno,real *density)
{
  PutNucnum(nuc);
  for(int i=0;i<nucnum;i++){
    Nuclide tmp;
    tmp.PutMatnum(matno[i]);
    if(density[i]>0.2){
      cout<<"# Error(?) in Medium::PutNuclide.\n";
      cout<<"# Too large values for number density : "<<density[i]<<"\n";
      cout<<"# Material ID is "<<matno[i]<<"\n";
      exit(0);
    };
    tmp.PutDensity(density[i]);
    PutNuclide(i,tmp);
  };
};

void Medium::PutNuclide(int nuc,vector<int> matno,vector<real> density)
{
  int *matno_in=new int[nuc];
  real *density_in=new real[nuc];

  for(int i=0;i<nuc;i++){
    matno_in[i]=matno[i];
    density_in[i]=density[i];
  };

  PutNuclide(nuc,matno_in,density_in);

  delete [] matno_in;
  delete [] density_in;
};

void Medium::PutNuclideNew(int nuc,int *matno,real *density,bool each_iso)
{
  if(!each_iso){
    PutNuclide(nuc,matno,density);
    return;
  };

  int isonuc=16;
  int isonucmat[]={
    120000,140000,190000,200000,220000, 240000,260000,280000,290000,310000,
    640000,400000,420000,480000,740000, 820000
  };
  int isonucnum[]={
    3,3,3,6,5,4, 4,5,2,2,7,
    5,7,8,4,4
  };

  vector<int> isonucpos(isonuc);
  int nn=0;
  for(int i=0;i<isonuc;i++){
    isonucpos[i]=nn;
    nn+=isonucnum[i];
  };

  int matin[]={
    120240,120250,120260, // Mg
    140280,140290,140300, // Si
    190390,190400,190410, // K
    200400,200420,200430,200440,200460,200480, // Ca
    220460,220470,220480,220490,220500, // Ti
    240500,240520,240530,240540, // Cr
    260540,260560,260570,260580, // Fe
    280580,280600,280610,280620,280640, // Ni
    290630,290650, // Cu
    310690,310710, // Ga
    641520,641540,641550,641560,641570,641580,641600, // Gd    
    400900,400910,400920,400940,400960, // Zr
    420920,420940,420950,420960,420970,420980,421000, // Mo
    481060,481080,481100,481110,481120,481130,481140,481160, // Cd
    741820,741830,741840,741860, // W
    822040,822060,822070,822080, // Pb
  };
  real ratio[]={
    0.7899, 0.1000, 0.1101, // Mg
    0.922296, 0.046832, 0.030872, // Si
    0.9326, 0.0001, 0.0673, // K
    0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187, // Ca
    0.0825, 0.0744, 0.7372, 0.0541, 0.0518, // Ti
    0.04345, 0.83789, 0.09501, 0.02365,
    0.05845, 0.91752, 0.02119, 0.00282,
    0.68077, 0.26223, 0.01140, 0.03634, 0.00926,
    0.6917, 0.3083, // Cu
    0.60108, 0.39892, // Ga
    0.002, 0.0218, 0.148, 0.2047, 0.1565, 0.2484, 0.2186, //Gd
    0.5145, 0.1122, 0.1715, 0.1738, 0.028, // Zr 
    0.1484, 0.0925, 0.1592, 0.1668, 0.0955, 0.2413, 0.0963, // Mo
    0.0125, 0.0089, 0.1249, 0.1280, 0.2413, 0.1222, 0.2873, 0.0749,// Cd
    0.265318, 0.143272, 0.306768, 0.284642, // W
    0.014, 0.241, 0.221, 0.524, // Pb
  };

  vector<int> matno_in;
  vector<real> density_in;
  for(int i=0;i<nuc;i++){
    int jpos=-1;
    for(int j=0;j<isonuc;j++){
      if(matno[i]==isonucmat[j])jpos=j;
    };
    if(jpos==-1){
      matno_in.push_back(matno[i]);
      density_in.push_back(density[i]);
    }else{
      int nnn=isonucnum[jpos];
      int nn2=isonucpos[jpos];
      for(int j=0;j<nnn;j++){
	matno_in.push_back(matin[nn2+j]);
	density_in.push_back(density[i]*ratio[nn2+j]);
      };
    };
  };

  int nuc_in=matno_in.size();
  PutNuclide(nuc_in,matno_in,density_in);
};

void Medium::PutNuclideFe(int nuc,int *matno,real *density)
{
  int isonuc=1;
  int isonucmat[]={260000};
  int isonucnum[]={4};

  vector<int> isonucpos(isonuc);
  int nn=0;
  for(int i=0;i<isonuc;i++){
    isonucpos[i]=nn;
    nn+=isonucnum[i];
  };

  int matin[]={
    260540,260560,260570,260580, // Fe
  };
  real ratio[]={
    0.05845, 0.91752, 0.02119, 0.00282,
  };

  vector<int> matno_in;
  vector<real> density_in;
  for(int i=0;i<nuc;i++){
    int jpos=-1;
    for(int j=0;j<isonuc;j++){
      if(matno[i]==isonucmat[j])jpos=j;
    };
    if(jpos==-1){
      matno_in.push_back(matno[i]);
      density_in.push_back(density[i]);
    }else{
      int nnn=isonucnum[jpos];
      int nn2=isonucpos[jpos];
      for(int j=0;j<nnn;j++){
	matno_in.push_back(matin[nn2+j]);
	density_in.push_back(density[i]*ratio[nn2+j]);
      };
    };
  };

  int nuc_in=matno_in.size();
  PutNuclide(nuc_in,matno_in,density_in);
};

void Medium::PutTemperatureForAllNuclide(real intmp)
{
  for(int i=0;i<nucnum;i++){
    nuc[i].PutTemperature(intmp);
  };
};

void Medium::ShowNuclideList()
{
  for(int i=0;i<nucnum;i++){
    cout<<" "<<i<<"    "<<nuc[i].GetMatnum()<<"\n";
  };
};

void Medium::ShowNumberDensity(bool zero_ignore)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  real sum=0.;
  for(int i=0;i<nucnum;i++){
    real d=nuc[i].GetDensity();
    if(!zero_ignore||(zero_ignore&&d>0.))cout<<"   "<<nuc[i].GetMatnum()<<"  "<<d<<"\n";
    sum+=d;
  };
  cout<<"    TOTAL   "<<sum<<"\n";
};

void Medium::ShowNumberDensityCBZStyle(MATIDTranslator &midt)
{
  for(int i=0;i<nucnum;i++){
    //cout<<nuc[i].GetMatnum()<<", ";
    cout<<"\""<<midt.Name(nuc[i].GetMatnum())<<"\", ";    
    if(i%5==4)cout<<"\n";
  };
  cout<<"\n";

  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<nucnum;i++){
    cout<<nuc[i].GetDensity()<<", ";
    if(i%5==4)cout<<"\n";
  };
  cout<<"\n";
};

void Medium::ShowMacroCrossSection1D()
{
  Macxs.ShowData1D(enband);
};

void Medium::ShowMicroCrossSection1D(int mat)
{
  if(ExistNuclide(mat)){
    GetNuclide(mat).GetMicxs().ShowData1D(enband);
  };
};

void Medium::ShowNeutronFlux()
{
  cout<<"# E_up     Flux       Flux/lethargy [normalized]\n";
  cout.setf(ios::scientific);
  cout.precision(4);
  real sum=flux[0].get_sum();
  for(int g=0;g<imax;g++){
    real e0=enband.get_dat(g);
    real e1=enband.get_dat(g+1);
    real letwid=log(e0/e1);
    cout<<e0<<" "<<flux[0].get_dat(g)<<" "<<flux[0].get_dat(g)/letwid<<"    "<<flux[0].get_dat(g)/letwid/sum<<"\n";
    if(g==imax-1)cout<<e0<<" "<<flux[0].get_dat(g)<<" "<<flux[0].get_dat(g)/letwid<<"    "<<flux[0].get_dat(g)/letwid/sum<<"\n";
  };

};


void Medium::CalHomoB1(real bsq,GroupData1D src)
{
  // 
  // Neutron current is treated as (iB)J
  //

  if(pl==0){
    cout<<"Caution !!!\n";
    cout<<"B1 spectrum calculation cannot be performed\n";
    cout<<"since P1 scattering XS is not given.\n";
    return;
  };

  if(bsq==0.0)bsq=1e-20;

  real *s0=new real[imax];
  real *s1=new real[imax];
  real *f0=new real[imax];
  real *f1=new real[imax];

  GroupData1D diagp0,diagp1;
  GroupData1D coef0,coef1;

  for(int i=0;i<imax;i++){
    s0[i]=0.; s1[i]=0.;
  };

  diagp0=Macxs.GetData2d(sigs).get_diag();
  diagp1=Macxs.GetData2d(sigs,1).get_diag();

  coef0=Macxs.GetData1d(sigt)-diagp0;
  if(pltot==0){
    coef1=Macxs.GetData1d(sigt)*3.0-diagp1;
  }else{
    coef1=Macxs.GetData1d(sigt,1)*3.0-diagp1;
  };

  int upg=GetUpScatteringGrp();
  int itermax=1;
  if(upg<imax-1)itermax=1000;

  vector<real> f0old(imax);
  for(int iter=0;iter<itermax;iter++){

    for(int i=0;i<imax;i++){
      f0old[i]=f0[i];
    };

    for(int i=0;i<imax;i++){
      f0[i]=((s0[i]+src.get_dat(i))*coef1.get_dat(i)-s1[i])
        /(coef0.get_dat(i)*coef1.get_dat(i)+bsq);
      f1[i]=(s1[i]+bsq*f0[i])/coef1.get_dat(i);
      for(int k=0;k<imax;k++){
        s0[k]+=Macxs.GetData2d(sigs).get_dat(i,k)*f0[i];
        s1[k]+=Macxs.GetData2d(sigs,1).get_dat(i,k)*f1[i];
      };
    };

    real errmax=0.;
    if(iter!=0){
      for(int i=0;i<imax;i++){
	real err=fabs(f0[i]/f0old[i]-1.);
	if(err>errmax)errmax=err;
      };
    };
    if(errmax<1e-5&&iter!=0)break;
    if(iter==itermax-1)break;

    for(int i=0;i<imax;i++){
      s0[i]=0.; s1[i]=0.;
    };
    for(int i=0;i<imax;i++){
      for(int k=0;k<i;k++){
        s0[k]+=Macxs.GetData2d(sigs).get_dat(i,k)*f0[i];
        s1[k]+=Macxs.GetData2d(sigs,1).get_dat(i,k)*f1[i];
      };
    };
    
  };

  flux[0].put_data(f0);
  flux[1].put_data(f1);

  delete [] f0;
  delete [] f1;
  delete [] s0;
  delete [] s1;
}

void Medium::CalHomoB1WithFixedSource(real bsq,GroupData1D &src)
{
  if(pl==0){
    cout<<"# !!! Caution !!!\n";
    cout<<"#   B1 spectrum calculation cannot be performed\n";
    cout<<"#   since P1 scattering XS is not given.\n";
    return;
  };

  if(bsq==0.0)bsq=1e-20;

  vector<real> s0(imax,0.);
  vector<real> s1(imax,0.);
  vector<real> f0(imax,0.);
  vector<real> f1(imax,0.);

  GroupData1D diagp0,diagp1;
  GroupData1D coef0,coef1;

  diagp0=Macxs.GetData2d(sigs).get_diag();
  diagp1=Macxs.GetData2d(sigs,1).get_diag();

  coef0=Macxs.GetData1d(sigt)-diagp0;
  coef1=Macxs.GetData1d(sigt,1)*3.0-diagp1;

  vector<real> f0old(imax);
  int itermax=10000;
  for(int iter=0;iter<itermax;iter++){

    real sum=0.;
    for(int i=0;i<imax;i++){
      sum+=f0[i]*Macxs.GetData1d(nusigf).get_dat(i);
    };

    for(int i=0;i<imax;i++){
      f0[i]=((s0[i]+Macxs.GetData1d(chi).get_dat(i)*sum+src.get_dat(i))*coef1.get_dat(i)-s1[i])
        /(coef0.get_dat(i)*coef1.get_dat(i)+bsq);
      f1[i]=(s1[i]+bsq*f0[i])/coef1.get_dat(i);
      for(int k=0;k<imax;k++){
        s0[k]+=Macxs.GetData2d(sigs).get_dat(i,k)*f0[i];
        s1[k]+=Macxs.GetData2d(sigs,1).get_dat(i,k)*f1[i];
      };
      s0[i]=0.;
      s1[i]=0.;
    };

    real errmax=0.;
    if(iter!=0){
      for(int i=0;i<imax;i++){
	real err=fabs(f0[i]/f0old[i]-1.);
	if(err>errmax)errmax=err;
      };
    };
    if(errmax<1e-5&&iter!=0)break;
    if(iter==itermax-1){
      cout<<"# Warning !!!\n";
      cout<<"# Fixed source calculation is NOT converged.\n";
      break;
    };

    for(int i=0;i<imax;i++){
      f0old[i]=f0[i];
    };
  };

  flux[0].put_data(f0);
  flux[1].put_data(f1);
}

void Medium::CalHomoB1Adjoint(real bsq,GroupData1D src)
{
  if(pl==0){
    cout<<"Caution !!!\n";
    cout<<"B1 spectrum calculation cannot be performed\n";
    cout<<"since P1 scattering XS is not given.\n";
    return;
  };

  if(bsq==0.0)bsq=1e-20;

  real *s0=new real[imax];
  real *s1=new real[imax];
  real *f0=new real[imax];
  real *f1=new real[imax];

  GroupData1D diagp0,diagp1;
  GroupData1D coef0,coef1;

  for(int i=0;i<imax;i++){
    s0[i]=0.; s1[i]=0.;
  };

  diagp0=Macxs.GetData2d(sigs).get_diag();
  diagp1=Macxs.GetData2d(sigs,1).get_diag();

  coef0=Macxs.GetData1d(sigt)-diagp0;
  coef1=Macxs.GetData1d(sigt,1)*3.0-diagp1;

  int upg=GetUpScatteringGrp();
  int itermax=1;
  if(upg<imax-1)itermax=1000;

  vector<real> f0old(imax);

  for(int iter=0;iter<itermax;iter++){

    for(int i=0;i<imax;i++){
      f0old[i]=f0[i];
    };

    for(int i=imax-1;i>=0;i--){
      f1[i]=(-s0[i]-src.get_dat(i)+coef0.get_dat(i)*s1[i])
        /(coef0.get_dat(i)*coef1.get_dat(i)+bsq);
      f0[i]=(s1[i]-coef1.get_dat(i)*f1[i]);
      for(int k=0;k<imax;k++){
        s0[k]+=Macxs.GetData2d(sigs).get_dat(k,i)*f0[i];
        s1[k]+=Macxs.GetData2d(sigs,1).get_dat(k,i)*f1[i];
      };
    };

    real errmax=0.;
    if(iter!=0){
      for(int i=0;i<imax;i++){
	real err=fabs(f0[i]/f0old[i]-1.);
	if(err>errmax)errmax=err;
      };
    };
    if(errmax<1e-5&&iter!=0)break;
    if(iter==itermax-1)break;

    for(int i=0;i<imax;i++){
      s0[i]=0.; s1[i]=0.;
    };
    for(int i=0;i<imax;i++){
      for(int k=0;k<i;k++){
        s0[i]+=Macxs.GetData2d(sigs).get_dat(i,k)*f0[k];
        s1[i]+=Macxs.GetData2d(sigs,1).get_dat(i,k)*f1[k];
      };
    };
  };

  flux[0].put_data(f0);
  flux[1].put_data(f1);

  delete [] f0;
  delete [] f1;
  delete [] s0;
  delete [] s1;
}

real Medium::CalKeff()
{
  //return (Macxs.GetData1d(nusigf)*flux[0]+Macxs.GetData1d(sign2n)*flux[0])/
  //       (Macxs.GetData1d(siga)*flux[0]+flux[1].get_sum());
  return Macxs.GetData1d(nusigf)*flux[0];
};

real Medium::CalKeffAdjoint()
{
  return Macxs.GetData1d(chi)*flux[0];
};

real Medium::CalKinf()
{
  return (Macxs.GetData1d(nusigf)*flux[0])/(Macxs.GetData1d(siga)*flux[0]-Macxs.GetData1d(sign2n)*flux[0]);
};

real Medium::BucklingSearch()
{
  real bsq=1e-10;
  real keff,kinf;
  for(int i=0;i<50;i++){
    CalHomoB1(bsq);
    keff=CalKeff();
    kinf=CalKinf();
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

void Medium::AddPseudoAbsorption(real bsq)
{
  for(int l=0;l<=pltot;l++){
    GroupData1D tmp=Macxs.GetData1d(sigt,l)+Macxs.GetData1d(d)*bsq;
    Macxs.GetData1d(sigt,l).copy(tmp);
  };
  GroupData1D tmp=Macxs.GetData1d(sigtr)+Macxs.GetData1d(d)*bsq;
  Macxs.GetData1d(sigtr).copy(tmp);
  tmp=Macxs.GetData1d(siga)+Macxs.GetData1d(d)*bsq;
  Macxs.GetData1d(siga).copy(tmp);
};

Medium Medium::Cond(int ngrp,vector<int> bgrp)
{
  int *bgrp2=new int[ngrp];
  for(int i=0;i<ngrp;i++){
    bgrp2[i]=bgrp[i];
  };
  return Cond(ngrp,bgrp2);
};

Medium Medium::Cond(int ngrp,int *bgrp)
{
  Medium ret(ngrp);
  ret=Cond(ngrp,bgrp,flux[0],flux[1]);
  return ret;
};

Medium Medium::Cond(int ngrp,vector<int> bgrp,GroupData1D fl,GroupData1D cu,bool micro)
{
  int *bgrp2=new int[ngrp];
  for(int i=0;i<ngrp;i++){
    bgrp2[i]=bgrp[i];
  };
  return Cond(ngrp,bgrp2,fl,cu,micro);
};

Medium Medium::Cond(int ngrp,int *bgrp,GroupData1D fl,GroupData1D cu,bool micro)  
{
  Medium ret(ngrp);
  ret.PutPL(pl,pltot);

  for(int i=1;i<ngrp;i++){
    if(bgrp[i]<=bgrp[i-1]){
      cout<<"# Error in Cond of Medium.cxx.\n";
      cout<<"$ Error in group structure.\n";
      exit(0);
    };
  };
  if(bgrp[ngrp-1]!=imax-1){
    cout<<"# Error in Cond of Medium.cxx\n";
    cout<<"# Lower boundary is inconsistent!\n";
    cout<<"  (lower  boundary) : "<<bgrp[ngrp-1]<<"\n";
    exit(0);
  };
  if(fl.get_imax()!=imax||cu.get_imax()!=imax){
    cout<<"# Error in Medium::Cond.\n";
    cout<<"# Energy group of weight function is inconsistent with cross sections.\n";
    exit(0);
  };

  ret.GetEnband().put_data(0,enband.get_dat(0));
  for(int i=0;i<ngrp;i++){
    ret.GetEnband().put_data(i+1,enband.get_dat(bgrp[i]+1));
  };

  ret.PutMacxs(Macxs.Cond(ngrp,bgrp,fl,cu));

  ret.GetFlux(0).copy(fl.CondSum(ngrp,bgrp));
  ret.GetFlux(1).copy(cu.CondSum(ngrp,bgrp));

  if(micro){
    ret.PutNucnum(nucnum);
    for(int i=0;i<nucnum;i++){
      if(nuc[i].GetGrp()!=-1){
        ret.PutNuclide(i,nuc[i].Cond(ngrp,bgrp,fl,cu));
      }else{
	ret.PutNuclide(i,nuc[i].GetMatnum(),nuc[i].GetDensity());
      };
    };
  }else{
    ret.PutNucnum(0);
  };

  return ret;
};

void Medium::CalSigtr(int sigt_l)
{
  //cout<<"  (Transport XS is calculated from ";
  if(pl>0){
    //cout<<"p1-total and p1-scattering)\n";
    Macxs.GetData1d(sigtr)=Macxs.GetData1d(sigt,sigt_l)-Macxs.GetData2d(sigs,1).get_sumx()*0.33333333;
  }else{
    //cout<<"averaged diffusion coefficient)\n";
    CalSigtrFromDav();
  };
};

void Medium::PutSigt0AsSigt1(int g0,int g1)
{
  if(g1==-1)g1=imax-1;
  for(int g=g0;g<=g1;g++){
    real st0=GetMacxs().GetData1d(sigt,0).get_dat(g);
    GetMacxs().GetData1d(sigt,1).put_data(g,st0);
  };
};

void Medium::CalSigtrDiagonal(int sigt_l)
{
  if(pl>0){
    //cout<<"p1-total and p1-scattering)\n";
    Macxs.GetData1d(sigtr)=Macxs.GetData1d(sigt,sigt_l)-Macxs.GetData2d(sigs,1).get_diag()*0.33333333;
  }else{
    cout<<"Error in Medium::CalSigtrDiagonal.\n";
    exit(0);
  };
};

void Medium::CalSigt()
{
  Macxs.GetData1d(sigt)=Macxs.GetData1d(siga)+Macxs.GetData2d(sigs).get_sumx()-Macxs.GetData1d(sign2n);
};

void Medium::CalSigtrFromDav()
{
  for(int i=0;i<imax;i++){
    Macxs.GetData1d(sigtr).put_data(i,0.333333333/Macxs.GetData1d(d).get_dat(i));
  };
};

void Medium::AdjustSelfScattering()
{
  for(int g=0;g<imax;g++){
    real st=Macxs.GetData1d(sigt).get_dat(g);
    real sa=Macxs.GetData1d(siga).get_dat(g);
    real ss=Macxs.GetData2d(sigs).get_sumx().get_dat(g);
    ss-=Macxs.GetData2d(sigs).get_dat(g,g);
    Macxs.GetData2d(sigs).put_data(g,g,st-sa-ss);
  };
};

void Medium::CalSigtFromDav()
{
  for(int i=0;i<imax;i++){
    Macxs.GetData1d(sigt).put_data(i,0.333333333/Macxs.GetData1d(d).get_dat(i));
  };
};

void Medium::CalDFromSigt(int n)
{
  for(int i=0;i<imax;i++){
    Macxs.GetData1d(d).put_data(i,0.33333333/Macxs.GetData1d(sigt,n).get_dat(i));
  };
};

void Medium::CalDFromSigtr()
{
  for(int i=0;i<imax;i++){
    Macxs.GetData1d(d).put_data(i,0.33333333/Macxs.GetData1d(sigtr).get_dat(i));
  };
};

void Medium::TransportApproximation()
{
  vector<bool> neg_ss(imax,false);
  bool neg_ss_all=false;
  for(int i=0;i<imax;i++){
    real tr=Macxs.GetData1d(sigtr).get_dat(i);
    real tt=Macxs.GetData1d(sigt).get_dat(i);
    real ss=Macxs.GetData2d(sigs).get_dat(i,i);
    real tmp=ss+tr-tt;
    //if(tmp>0.){
      Macxs.GetData1d(sigt).put_data(i,tr);
      Macxs.GetData1d(sigtr).put_data(i,tt);
      Macxs.GetData2d(sigs).put_data(i,i,tmp);
      //  };

    if(tmp<0.){
      //cout<<"!! Warning !! Negative self-scattering cross section in group "<<i<<"\n";
      neg_ss[i]=true;
      neg_ss_all=true;

      // self-scattering is set to be 0.0
      /*
      //if(i==179){
      Macxs.GetData1d(sigt).put_data(i,tt-ss);
      Macxs.GetData1d(sigtr).put_data(i,tt);
      Macxs.GetData2d(sigs).put_data(i,i,0.);
      //};
      */
    };
  };

  if(neg_ss_all){
    cout<<"# !! Warning !! Negative self-scattering XSs are found in groups ";
    for(int i=0;i<imax;i++){
      if(neg_ss[i])cout<<i<<" ";
    };
    cout<<"\n";
  };
}

void Medium::ConsistentPApproximation()
{
  for(int j=1;j<=pl;j++){
    for(int k=0;k<imax;k++){
      real tmp=GetDataSigs(j,k,k);
      real tmp2=GetDataSigt(0,k)-GetDataSigt(1,k);
      real tmp3=tmp+(j*2+1.)*tmp2;
      GetSigs(j).put_data(k,k,tmp3);
    };
  };
};

void Medium::CalRemovalCrossSection()
{
  GroupData1D ret(imax);
  ret=Macxs.GetSigt()-Macxs.GetSigs().get_diag();
  //ret=Macxs.GetSiga()+Macxs.GetSigs().get_sumx()
  //  -Macxs.GetSigs().get_diag()-Macxs.GetSign2n();
  // The above assumes that n2n does not have self-scattering
  Macxs.GetData1d(sigr).copy(ret);
};

Nuclide &Medium::GetNuclide(int matnum)
{
  for(int i=0;i<nucnum;i++){
    if(matnum==nuc[i].GetMatnum())return nuc[i];
  };
  cout<<"There is not the nuclide in the medium class.\n";
  cout<<"You requested "<<matnum<<"\n";
  for(int i=0;i<nucnum;i++){
    cout<<nuc[i].GetMatnum()<<" ";
  };
  cout<<"\n";
  exit(0);
};

real Medium::GetNuclideDensity(int matnum)
{
  if(!ExistNuclide(matnum)){
    return 0.;
  }else{
    return GetNuclide(matnum).GetDensity();
  };
};

real Medium::GetNuclideDensity(int num, vector<int> &matid)
{
  real sum=0.;
  for(int i=0;i<num;i++){
    sum+=GetNuclideDensity(matid[i]);
  };
  return sum;
};

real Medium::GetNuclideDensity(int num, vector<string> &matname, MATIDTranslator &midt)
{
  real sum=0.;
  for(int i=0;i<num;i++){
    sum+=GetNuclideDensity(midt.ID(matname[i]));
  };
  return sum;
};

bool Medium::ExistNuclide(int matnum)
{
  for(int i=0;i<nucnum;i++){
    if(matnum==nuc[i].GetMatnum())return true;
  };
  return false;
};

int Medium::SearchNuclide(int matnum)
{
  for(int i=0;i<nucnum;i++){
    if(matnum==nuc[i].GetMatnum())return i;
  };
  return -1;
};

void Medium::PutMacxs(GroupDataSet macinp)
{
  if(macinp.GetName()!="MacroCrossSection"){
    cout<<"Error in Medium class.\n";
    cout<<"You are trying to put GroupDataSet into Medium class.\n";
    cout<<"The name of GroupDataSet should be MacroCrossSection.\n";
    cout<<"Your case is "<<macinp.GetName()<<"\n";
    exit(0);
  };
  Macxs=macinp;
};

void Medium::PutData(enum xstype ss,vector<real> in)
{
  real *inp=new real[imax];
  for(int i=0;i<imax;i++){
    inp[i]=in[i];
  };
  PutData(ss,inp);
  delete [] inp;
};

void Medium::PutDataSigt(int l,vector<real> in)
{
  real *inp=new real[imax];
  for(int i=0;i<imax;i++){
    inp[i]=in[i];
  };
  PutDataSigt(l,inp);
  delete [] inp;
};

void Medium::PutDataSigs(int l,vector<real> in)
{
  real *inp=new real[imax*imax];
  for(int i=0;i<imax*imax;i++){
    inp[i]=in[i];
  };
  PutDataSigs(l,inp);
  delete [] inp;
};

int Medium::GetUpScatteringGrp(int plmax)
{
  if(imax==1)return -1;
  //for(int i=1;i<imax-1;i++){
  for(int i=1;i<imax;i++){
    for(int pl=0;pl<=plmax;pl++){
      for(int j=0;j<i;j++){
	if(Macxs.GetSigs(pl).get_dat(i,j)!=0.)return i;
      };
    };
  };
  return -1;
};

int Medium::GetUpScatteringSinkGrp(int plmax)
{
  if(imax==1)return -1;
  for(int i=0;i<imax-1;i++){
    for(int pl=0;pl<=plmax;pl++){
      //for(int j=i+1;j<imax-1;j++){
      for(int j=i+1;j<imax;j++){
	if(Macxs.GetSigs(pl).get_dat(j,i)!=0.)return i;
      };
    };
  };
  return -1;
};

void Medium::PutUpScatteringToSelfScattering()
{
  for(int l=0;l<=pl;l++){
    for(int j=1;j<imax;j++){
      real tmp=0.;
      for(int k=0;k<j;k++){
	tmp+=Macxs.GetSigs(l).get_dat(j,k);
	Macxs.GetSigs(l).put_data(j,k,0.);
      };
      Macxs.GetSigs(l).add_data(j,j,tmp);
    };
  };
};

void Medium::PutUpScatteringToSelfScattering(int g0, int g1)
{
  for(int l=0;l<=pl;l++){
    for(int j=g0;j<=g1;j++){
      real tmp=0.;
      for(int k=0;k<j;k++){
	tmp+=Macxs.GetSigs(l).get_dat(j,k);
	Macxs.GetSigs(l).put_data(j,k,0.);
      };
      Macxs.GetSigs(l).add_data(j,j,tmp);
    };
  };
};

void Medium::CalMacroFromMicro(bool chi_calc)
{
  // +++ IMPORTANT NOTICE FOR FUTURE IMPLEMENTATION +++
  //
  // Transport cross section is calculated with NOT Sig_{t,1}-Sig_{s,1} BUT Sig_{t,1}-Sig_{e,1} here,
  // and diffusion coefficients are calculated from this tranrport cross section.
  //
  // In the [CalSigtr] method of the [Medium] class, however,
  // transport cross section is redefined from Sig_{t,1}-Sig_{s,1},
  // so inconsistency should occur between transport cross section and diffusion coefficient.


  GroupData1D tmp(imax);
  GroupData2D tmp2(imax,imax);

  // nusigf
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
        tmp=tmp+nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*nuc[i].GetDensity();
      };
    };
  };
  Macxs.GetData1d(nusigf).copy(tmp);
  //tmp.show_self();

  // siga
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    //cout<<i<<"\n";
    if(nuc[i].GetGrp()!=-1){
      tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigc)+nuc[i].GetMicxs().GetData1d(sigf))*nuc[i].GetDensity();
      //nuc[i].GetMicxs().GetData1d(sigc).show_self();
      //nuc[i].GetMicxs().GetData1d(sigf).show_self();
    };


  };
  Macxs.GetData1d(siga).copy(tmp);
  //tmp.show_self();

  // sign2n
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      tmp=tmp+nuc[i].GetMicxs().GetData1d(sign2n)*nuc[i].GetDensity();
    };
  };
  Macxs.GetData1d(sign2n).copy(tmp);

  // sigt_p1
  if(pltot>0){
    tmp.set_zero();
    for(int i=0;i<nucnum;i++){
      if(nuc[i].GetGrp()!=-1){
        tmp=tmp+nuc[i].GetMicxs().GetData1d(sigt,1)*nuc[i].GetDensity();
      };
    };
    Macxs.GetData1d(sigt,1).copy(tmp);
  };

  // sigtr

  //pl=0;  pltot=0;
  
  int plt=0;
  if(pltot>0)plt=1;
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      real den=nuc[i].GetDensity();
      if(pl>0){
        tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigt,plt)-nuc[i].GetMicxs().GetData2d(sigel,1).get_sumx()*0.33333333)*den;
      }else{
	tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigt,plt)-nuc[i].GetMicxs().GetData1d(mu).mult(nuc[i].GetMicxs().GetData2d(sigel,0).get_sumx()))*den;
      };
    };
    Macxs.GetData1d(sigtr).copy(tmp);
  };

  /*
  if(pl>0&&pltot>0){
    tmp.set_zero();
    for(int i=0;i<nucnum;i++){
      if(nuc[i].GetGrp()!=-1){
        tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigt,1)-nuc[i].GetMicxs().GetData2d(sigel,1).get_sumx()*0.33333333)*nuc[i].GetDensity();
      };
    };
    Macxs.GetData1d(sigtr).copy(tmp);
  };
  */

  // chi
  if(chi_calc){
    
  tmp.set_zero();
  vector<real> fission_rate(nucnum,0.);
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
        fission_rate[i]=
          nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*flux[0]*nuc[i].GetDensity();
      }; 
    };
  };
  for(int i=0;i<nucnum;i++){
    if(fission_rate[i]>0.){
      tmp=tmp+nuc[i].GetMicxs().GetData1d(chi)*fission_rate[i];
    };
  };
  real sum=tmp.get_sum();

  if(sum!=0.){
    sum=1./sum;
    tmp=tmp*sum;
  }else{
    real tt=1./imax;
    for(int i=0;i<imax;i++){
      tmp.put_data(i,tt);
    };     
  };

  Macxs.GetData1d(chi).copy(tmp);

  };

  // Scattering matrix
  for(int l=0;l<=pl;l++){
    tmp2.set_zero();
    for(int i=0;i<nucnum;i++){
      if(nuc[i].GetGrp()!=-1){
        real den=nuc[i].GetDensity();
        for(int j=0;j<3;j++){
  	  if(l<nuc[i].GetMicxs().GetDim2d(j)){
	    tmp2=tmp2+nuc[i].GetMicxs().GetData2d(j,l)*den;
	  };
	};
      };
    };
    Macxs.GetData2d(sigs,l).copy(tmp2); 
  };

  CalSigt();
  CalDFromSigtr();
};

void Medium::AddMacroFromMicro(int id)
{
  // Fission spectrum and diffusion coefficient are NOT yet implemented.

  
  int pos=-1;
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetMatnum()==id)pos=i;
  };
  if(pos==-1){
    cout<<"# Error in Medium::AddMacroFromMicro.\n";
    cout<<"# The nuclide "<<id<<" cannot be found in this medium.\n";
    exit(0);
  };

  if(nuc[pos].GetGrp()==-1){
    cout<<"# Error in Medium::AddMacroFromMicro.\n";
    cout<<"# The nuclide "<<id<<" does not have any microscopic cross section data.\n";
    exit(0);
  };

  GroupData1D tmp(imax);
  GroupData2D tmp2(imax,imax);

  real den=nuc[pos].GetDensity();  
  
  // nusigf
  if(nuc[pos].GetMicxs().GetData1d(nu).get_dat(0)>0.){
    tmp=Macxs.GetData1d(nusigf);
    tmp=tmp+nuc[pos].GetMicxs().GetData1d(nu).mult(nuc[pos].GetMicxs().GetData1d(sigf))*den;
    Macxs.GetData1d(nusigf).copy(tmp);
  };

  // siga
  tmp=Macxs.GetData1d(siga);
  tmp=tmp+(nuc[pos].GetMicxs().GetData1d(sigc)+nuc[pos].GetMicxs().GetData1d(sigf))*den;
  Macxs.GetData1d(siga).copy(tmp);

  // sign2n
  tmp=Macxs.GetData1d(sign2n);
  tmp=tmp+nuc[pos].GetMicxs().GetData1d(sign2n)*den;
  Macxs.GetData1d(sign2n).copy(tmp);

  // sigt_p0
  tmp=Macxs.GetData1d(sigt);
  //tmp=tmp+nuc[pos].GetMicxs().GetData1d(sigt)*den;
  tmp=tmp+nuc[pos].GetMicxs().GetData1d(sigc)*den;
  tmp=tmp+nuc[pos].GetMicxs().GetData1d(sigf)*den;
  tmp=tmp+nuc[pos].GetMicxs().GetData2d(sigel).get_sumx()*den;
  tmp=tmp+nuc[pos].GetMicxs().GetData2d(siginel).get_sumx()*den;
  tmp=tmp+nuc[pos].GetMicxs().GetData2d(sign2n).get_sumx()*den*0.5;          
  Macxs.GetData1d(sigt).copy(tmp);

  // sigt_p1
  if(pltot>0){
    tmp=Macxs.GetData1d(sigt,1);
    tmp=tmp+nuc[pos].GetMicxs().GetData1d(sigt,1)*den;
    Macxs.GetData1d(sigt,1).copy(tmp);
  };

  // sigtr

  //pl=0;  pltot=0;
  
  int plt=0;
  if(pltot>0)plt=1;
  tmp=Macxs.GetData1d(sigtr);
  if(pl>0){
    tmp=tmp+(nuc[pos].GetMicxs().GetData1d(sigt,plt)-nuc[pos].GetMicxs().GetData2d(sigel,1).get_sumx()*0.33333333)*den;
  }else{
    tmp=tmp+(nuc[pos].GetMicxs().GetData1d(sigt,plt)-nuc[pos].GetMicxs().GetData1d(mu).mult(nuc[pos].GetMicxs().GetData2d(sigel,0).get_sumx()))*den;
  };
  Macxs.GetData1d(sigtr).copy(tmp);

  // chi ... skip at present
  /*
  tmp.set_zero();
  vector<real> fission_rate(nucnum,0.);
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
        fission_rate[i]=
          nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*flux[0]*nuc[i].GetDensity();
      };
    };
  };
  for(int i=0;i<nucnum;i++){
    if(fission_rate[i]>0.){
      tmp=tmp+nuc[i].GetMicxs().GetData1d(chi)*fission_rate[i];
    };
  };
  real sum=tmp.get_sum();
  if(sum!=0.){
    sum=1./sum;
    tmp=tmp*sum;
  }else{
    real tt=1./imax;
    for(int i=0;i<imax;i++){
      tmp.put_data(i,tt);
    };     
  };
  Macxs.GetData1d(chi).copy(tmp);
  */

  // Scattering matrix
  for(int l=0;l<=pl;l++){
    tmp2=Macxs.GetData2d(sigs,l);
    for(int j=0;j<3;j++){
      if(l<nuc[pos].GetMicxs().GetDim2d(j)){
	tmp2=tmp2+nuc[pos].GetMicxs().GetData2d(j,l)*den;
      };
    };
    Macxs.GetData2d(sigs,l).copy(tmp2); 
  };

  // Diffusion coefficient ... skip at present
  //CalDFromSigtr();
};

void Medium::FissionSpectrumVectorReconstruction()
{
  GroupData1D tmp(imax);

  tmp.set_zero();
  vector<real> fission_rate(nucnum,0.);
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
        fission_rate[i]=
          nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*flux[0]*nuc[i].GetDensity();
      };
    };
  };
  for(int i=0;i<nucnum;i++){
    if(fission_rate[i]>0.){
      tmp=tmp+nuc[i].GetMicxs().GetData1d(chi)*fission_rate[i];
    };
  };
  real sum=tmp.get_sum();
  if(sum!=0.){
    sum=1./sum;
    tmp=tmp*sum;
  }else{
    real tt=1./imax;
    for(int i=0;i<imax;i++){
      tmp.put_data(i,tt);
    };     
  };
  Macxs.GetData1d(chi).copy(tmp);
};

void Medium::CalMacroFromMicroSimple()
{
  GroupData1D tmp(imax);

  // nusigf
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
        tmp=tmp+nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*nuc[i].GetDensity();
      };
    };
  };
  Macxs.GetData1d(nusigf).copy(tmp);

  // siga
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigc)+nuc[i].GetMicxs().GetData1d(sigf))*nuc[i].GetDensity();
    };
  };
  Macxs.GetData1d(siga).copy(tmp);

  // sign2n
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      tmp=tmp+nuc[i].GetMicxs().GetData1d(sign2n)*nuc[i].GetDensity();
    };
  };
  Macxs.GetData1d(sign2n).copy(tmp);

  // sigtr
  tmp.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      tmp=tmp+(nuc[i].GetMicxs().GetData1d(sigt,1)-nuc[i].GetMicxs().GetData1d(mu))*nuc[i].GetDensity();
    };
  };
  Macxs.GetData1d(sigtr).copy(tmp);

  // chi
  tmp.set_zero();
  /*
  vector<real> fission_rate(nucnum,0.);
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetMicxs().GetData1d(nu).get_dat(0)>0.){
      fission_rate[i]=
        nuc[i].GetMicxs().GetData1d(nu).mult(nuc[i].GetMicxs().GetData1d(sigf))*flux[0]*nuc[i].GetDensity();
    };
  };
  */
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetMatnum()==942390)tmp=tmp+nuc[i].GetMicxs().GetData1d(chi)*0.9;
    if(nuc[i].GetMatnum()==922380)tmp=tmp+nuc[i].GetMicxs().GetData1d(chi)*0.1;
  };
  real sum=tmp.get_sum();
  if(sum!=0.){
    sum=1./sum;
    tmp=tmp*sum;
  }else{
    real tt=1./imax;
    for(int i=0;i<imax;i++){
      tmp.put_data(i,tt);
    };     
  };

  Macxs.GetData1d(chi).copy(tmp);

  // Only-P0 scattering matrix without up-scattering
  GroupData2D tmp2(imax,imax);
  tmp2.set_zero();
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      real den=nuc[i].GetDensity();
      for(int j=0;j<3;j++){
        for(int g=0;g<imax;g++){
  	  for(int g2=g;g2<imax;g2++){
  	    tmp2.add_data(g,g2,nuc[i].GetMicxs().GetData2d(j,0).get_dat(g,g2)*den);
	  };
	};
      };
    };
  };
  Macxs.GetData2d(sigs).copy(tmp2); 

  CalSigt();
  CalDFromSigtr();
};

void Medium::PutFictitiousVacuum(int grp,int ipl,real fic_sig)
{
  PutImax(grp);
  PutPL(ipl);
  for(int l=0;l<=ipl;l++){
    Macxs.GetData2d(sigs,l).set_zero();
  };

  for(int i=0;i<imax;i++){
    Macxs.GetData1d(siga).put_data(i,fic_sig);
    Macxs.GetData1d(sigt).put_data(i,fic_sig);
    Macxs.GetData1d(sigtr).put_data(i,fic_sig);
    Macxs.GetData1d(sigt,1).put_data(i,fic_sig);
  };
  CalDFromSigt();
  Macxs.GetData1d(dr).copy(Macxs.GetData1d(d));
  Macxs.GetData1d(dz).copy(Macxs.GetData1d(d));
};

void Medium::PutLowestDownScatGrp()
{
  for(int i=0;i<imax;i++){
    for(int j=imax-1;j>i;j--){
      real tmp=Macxs.GetData2d(sigs).get_dat(i,j);
      if(tmp>0.){
        LowestDownScatGrp[i]=j;
	break;
      };
    };
  };
};

void Medium::PutHighestUpScatGrp()
{
  if(imax!=1){
    for(int i=0;i<imax;i++){
      //for(int j=1;j<i;j++){ // Bug found in 2020/1/15
      for(int j=i+1;j<imax;j++){ // fixed in 2020/1/15
        if(Macxs.GetSigs().get_dat(j,i)!=0.){
	  HighestUpScatGrp=i;
	  return;
        };
      };
    };
  };
  HighestUpScatGrp=-1;
  return;
};

void Medium::SetZeroChi()
{
  for(int i=0;i<imax;i++){
    if(Macxs.GetData1d(nusigf).get_dat(i)>0.)return;
  };
  Macxs.GetData1d(chi).set_zero();
};

GroupData1D Medium::GetMacroSigf()
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    real sum=0.;
    for(int k=0;k<nucnum;k++){
      if(nuc[i].GetGrp()!=-1){
        real sf=nuc[k].GetMicxs().GetData1d(sigf).get_dat(i);
        if(sf>0.){
          real den=nuc[k].GetDensity();
          sum+=den*sf;
	};
      };
    };
    ret.put_data(i,sum);
  };
  return ret;
};

GroupData1D Medium::GetMacroSigc()
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    real sum=0.;
    for(int k=0;k<nucnum;k++){
      if(nuc[i].GetGrp()!=-1){
        real sf=nuc[k].GetMicxs().GetData1d(sigc).get_dat(i);
        if(sf>0.){
          real den=nuc[k].GetDensity();
          sum+=den*sf;
        };
      };
    };
    ret.put_data(i,sum);
  };
  return ret;
};

GroupData2D Medium::GetMacro2D(enum xstype ss)
{
  int upscg=GetUpScatteringGrp(0);

  real *res=new real[imax*imax];
  for(int i=0;i<imax*imax;i++){res[i]=0.;};
  for(int k=0;k<nucnum;k++){
    if(nuc[k].GetGrp()!=-1){
      real den=nuc[k].GetDensity();
      if(den>0.){
        for(int i=0;i<imax;i++){
          int ist=0;
  	  if(upscg==-1)ist=i;
          for(int j=ist;j<imax;j++){
            real sf=nuc[k].GetMicxs().GetData2d(ss).get_dat(i,j);
            if(sf>0.){
              res[i*imax+j]+=den*sf;
	    };
	  };
        };
      };
    };
  };
  GroupData2D ret(imax,imax);
  ret.put_data(res);
  delete [] res;
  return ret;
};

void Medium::SetZeroSelfScattering()
{
  for(int g=0;g<imax;g++){
    real tmp=GetMacxs().GetData2d(sigs).get_dat(g,g);
    real tmp1=GetMacxs().GetData2d(sigs,1).get_dat(g,g)*0.33333333;
    real vmu=tmp1/tmp;
    for(int l=0;l<=pltot;l++){
      GetMacxs().GetData1d(sigt,l).add_data(g,-tmp);
    };
    for(int l=0;l<=pl;l++){
      GetMacxs().GetData2d(sigs,l).put_data(g,g,0.);
    };
    if(g!=imax-1){
      real tmp2=tmp*(1.-vmu);
      GetMacxs().GetData2d(sigs).add_data(g,g+1,tmp2);
      for(int l=0;l<=pltot;l++){
	GetMacxs().GetData1d(sigt,l).add_data(g,tmp2);
      };
    };
  };
};

void Medium::MicxsVectorClear2DData()
{
  for(int i=0;i<nucnum;i++){
    nuc[i].GetMicxs().AllVectorClear2DData();
  };
};

void Medium::MicxsVectorClear()
{
  for(int i=0;i<nucnum;i++){
    nuc[i].GetMicxs().AllVectorClear();
  };
};

void Medium::SetZeroNumberDensity()
{
  for(int i=0;i<nucnum;i++){
    nuc[i].PutDensity(0.);
  };
};

void Medium::AdjustNumberDensity()
{
  int nn=3+3+4+2+3+6+5+2+4+4+5+2+2+2+7+5+7+8+4+4;
  int matin[]={
    120240,120250,120260, // Mg
    140280,140290,140300, // Si
    160320,160330,160340,160360, // S
    170350,170370, // Cl
    190390,190400,190410, // K
    200400,200420,200430,200440,200460,200480, // Ca
    220460,220470,220480,220490,220500, // Ti
    230500,230510, // V
    240500,240520,240530,240540, // Cr
    260540,260560,260570,260580, // Fe
    280580,280600,280610,280620,280640, // Ni
    290630,290650, // Cu
    310690,310710, // Ga
    471070,471090, // Ag
    641520,641540,641550,641560,641570,641580,641600, // Gd    
    400900,400910,400920,400940,400960, // Zr
    420920,420940,420950,420960,420970,420980,421000, // Mo
    481060,481080,481100,481110,481120,481130,481140,481160, // Cd
    741820,741830,741840,741860, // W
    822040,822060,822070,822080, // Pb
  };
  real ratio[]={
    0.7899, 0.1000, 0.1101, // Mg
    0.922296, 0.046832, 0.030872, // Si 
    0.9499, 0.0075, 0.0425, 0.0001, // S
    0.7578, 0.2422, // Cl
    0.9326, 0.0001, 0.0673, // K
    0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187, // Ca
    0.0825, 0.0744, 0.7372, 0.0541, 0.0518, // Ti
    0.0025, 0.9975, // V
    0.04345, 0.83789, 0.09501, 0.02365,
    0.05845, 0.91752, 0.02119, 0.00282,
    0.68077, 0.26223, 0.01140, 0.03634, 0.00926,
    0.6917, 0.3083, // Cu
    0.60108, 0.39892, // Ga
    0.51839, 0.48161, // Ag
    0.002, 0.0218, 0.148, 0.2047, 0.1565, 0.2484, 0.2186, //Gd
    0.5145, 0.1122, 0.1715, 0.1738, 0.028, // Zr
    0.1484, 0.0925, 0.1592, 0.1668, 0.0955, 0.2413, 0.0963, // Mo
    0.0125, 0.0089, 0.1249, 0.1280, 0.2413, 0.1222, 0.2873, 0.0749,// Cd
    0.265318, 0.143272, 0.306768, 0.284642, // W
    0.014, 0.241, 0.221, 0.524, // Pb
  };

  for(int i=0;i<nucnum;i++){
    int id=GetNuclideID(i);
    if(id<10000){
      cout<<"# Warning in Medium::AdjustNumberDensity.\n";
      cout<<"# This method is prepared for CBZ nuclide ID,\n";
      cout<<"# NOT for ENDF nuclide ID, which you are intending.\n";
      cout<<"# ID is "<<id<<"\n";
    };
    for(int j=0;j<nn;j++){
      if(id==matin[j]){
	real org=GetNuclideInTurn(i).GetDensity();
        GetNuclideInTurn(i).PutDensity(org*ratio[j]);
      };
    };
  };
};

void Medium::FactorizeNumberDensity(real factor)
{
  for(int i=0;i<nucnum;i++){
    real org=GetNuclideInTurn(i).GetDensity();
    GetNuclideInTurn(i).PutDensity(org*factor);
  };
};

void Medium::CalMuFromSigel_p1()
{
  for(int i=0;i<nucnum;i++){
    if(nuc[i].GetGrp()!=-1){
      //nuc[i].GetMicxs().GetData1d(mu).copy(nuc[i].GetMicxs().GetData2d(sigel,1).get_sumx()*0.33333333);
      GroupData1D s1=nuc[i].GetMicxs().GetData2d(sigel,1).get_sumx()*0.3333333333;
      GroupData1D s0=nuc[i].GetMicxs().GetData2d(sigel,0).get_sumx();
      GroupData1D inp(imax);
      for(int g=0;g<imax;g++){
        inp.put_data(g,s1.get_dat(g)/s0.get_dat(g));
      };
      nuc[i].GetMicxs().GetData1d(mu).copy(inp);
    };
  };
};

void Medium::CorrectionForLastGroupSigt1()
{
  real st0=GetMacxs().GetData1d(sigt).get_dat(imax-1);
  GetMacxs().GetData1d(sigt,1).put_data(imax-1,st0);
};

void Medium::PrintMacroSectionTable()
{
  cout.setf(ios::scientific);
  cout.precision(3);

  cout<<"# Upper   Total     Produc-   Absorp-   (n,2n)    Chi       Diffusion Transport\n";
  cout<<"# energy            tion      tion                          Coef.\n";
  for(int g=0;g<imax;g++){
    cout<<enband.get_dat(g)<<" ";
    cout<<Macxs.GetData1d(sigt).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(nusigf).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(siga).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(sign2n).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(chi).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(d).get_dat(g)<<" ";
    cout<<Macxs.GetData1d(sigtr).get_dat(g)<<" ";
    cout<<"\n";
  };
};

void Medium::ExtendPLTOrder(int plin)
{
  if(plin<pltot){
    cout<<"+ Medium::ExtendPLTOrder.\n";
    cout<<"+   input   : "<<plin<<"\n";
    cout<<"+   existed : "<<pltot<<"\n";
    return;
  };

  vector<GroupData1D> storg(pltot+1);
  for(int i=0;i<=pltot;i++){
    storg[i].copy(Macxs.GetData1d(sigt,i));
  };

  int pltorg=pltot;
  pltot=plin;
  Macxs.PutDim1d(sigt,pltot+1);
  for(int i=0;i<=pltorg;i++){
    Macxs.GetData1d(sigt,i).copy(storg[i]);
  };
  for(int i=pltorg+1;i<=pltot;i++){
    Macxs.GetData1d(sigt,i).copy(storg[pltorg]);
  };
};

void Medium::MATIDTranslationFromENDFID(MATIDTranslator &midt)
{
  for(int i=0;i<nucnum;i++){
    int id_old=nuc[i].GetMatnum();
    if(id_old<10000){
      id_old=midt.GetMATIDFromENDFID(id_old);
      nuc[i].PutMatnum(id_old);
    };
  };
};

GroupData2D Medium::CalAmatForAbemie(real keff)
{
  GroupData1D a1(imax);
  GroupData1D a2(imax);
  GroupData1D a3(imax);
  GroupData2D ret(imax,imax);

  a1=Macxs.GetData1d(chi)/Macxs.GetData1d(d);
  a2=Macxs.GetData1d(nusigf)/keff;
  a3=(Macxs.GetData2d(sigs).get_sumx()
     -Macxs.GetData2d(sigs).get_diag()
     +Macxs.GetData1d(siga))/Macxs.GetData1d(d);

  for (int i=0;i<imax;i++){
    for(int j=0;j<imax;j++){
      real x=a1.get_dat(i)*a2.get_dat(j);
      if(i==j){
	x-=a3.get_dat(i);
      };
      if(j<i){
	x+=Macxs.GetData2d(sigs).get_dat(j,i)
          /Macxs.GetData1d(d).get_dat(i);
      };
      ret.put_data(i,j,x);
    };
  };

  return ret;
};

bool Medium::IsFissionable()
{
  if(fissionable_flag==2){
    return true;
  }else if(fissionable_flag==1){
    return false;
  };

  for(int i=0;i<nucnum;i++){
    if(nuc[i].IsFissionable()){
      fissionable_flag=2;
      return true;
    };
  };
  fissionable_flag=1;
  return false;
};

void Medium::NuclideClear()
{
  /*
  for(int i=0;i<nucnum;i++){
    nuc[i].PutGrp(1);
    nuc[i].GetMicxs().Init("MicroCrossSection");
    //nuc[i].GetMicxs().AllVectorClear();
  };
  */

  vector<Nuclide>().swap(nuc);
  nucnum=0;
};

void Medium::CopyWithoutMicroscopicScatteringMatrices(Medium &sec)
{
  PutImax(sec.GetImax());
  PutPL(sec.GetPL(), sec.GetPLT());

  enband=sec.GetEnband();
  Macxs=sec.GetMacxs();
  flux[0]=sec.GetFlux(0);
  flux[1]=sec.GetFlux(0);

  for(int i=0;i<sec.GetNucnum();i++){
    Nuclide nuctmp(true);
    nuctmp.PutMatnum(sec.GetNuclideInTurn(i).GetMatnum());
    nuctmp.PutDensity(sec.GetNuclideInTurn(i).GetDensity());
    nuctmp.PutTemperature(sec.GetNuclideInTurn(i).GetTemperature());
    if(sec.GetNuclideInTurn(i).GetGrp()>0){
      nuctmp.PutGrp(sec.GetNuclideInTurn(i).GetGrp());
      nuctmp.GetMicxs().GetData1d(sigf)=sec.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf);
      nuctmp.GetMicxs().GetData1d(sigc)=sec.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc);
      nuctmp.GetMicxs().GetData1d(nu)=sec.GetNuclideInTurn(i).GetMicxs().GetData1d(nu);
    };
    AddNuclide(nuctmp);
  };
};

void Medium::FourierAnalysisForUpScattering()
{
  GroupData2D t(imax,imax);
  GroupData2D sdo(imax,imax);
  GroupData2D suo(imax,imax);

  t.set_zero();
  sdo.set_zero();
  suo.set_zero();  
  for(int g=0;g<imax;g++){
    t.put_data(g,g,Macxs.GetData1d(sigt).get_dat(g));
    for(int g2=0;g2<imax;g2++){
      real tmp=Macxs.GetData2d(sigs).get_dat(g2,g);
      if(g2<=g){
        sdo.put_data(g,g2,tmp);
      }else{
        suo.put_data(g,g2,tmp);	
      };
    };
  };

  t=t-sdo;
  t.DoInverse();
  t=t*suo;

  // Initial vector
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,1.);
  };

  int itermax=100;

  GroupData1D ret2;
  real eigen=0;
  for(int i=0;i<itermax;i++){
    ret2=t*ret;
    real norm=ret.GetEuclideanNorm();
    real norm2=ret2.GetEuclideanNorm();
    real errk=fabs((norm2/norm)/eigen-1.);
    eigen=norm2/norm;
    /*
    cout.setf(ios::scientific);
    cout.precision(10);
    cout<<eigen<<"\n";
    */
    if(eigen==0.||errk<1e-6)break;
    ret=ret2*(1./norm2);    
  };

  cout<<"# Eigenvalue for two-grid matrix : "<<eigen<<"\n";
  //real norm_inv=1./ret2.GetEuclideanNorm();
  real norm_inv=1./ret2.get_sum();
  ret2=ret2*norm_inv;
  //if(eigen<1e-6)ret2.set_zero();
  if(eigen<0.5)ret2.set_zero();
  //if(eigen<0.9)ret2.set_zero();  
  Macxs.GetData1d(dr).copy(ret2);
  Macxs.GetData1d(dz).put_data(0,eigen);  
  //ret2.show_self();

  //t.CalEigenvaluesByLeverrieFaddeev(ret);
  //ret.show_self();
};
