#include <cstdlib>
#include "GammaTool.h"

void DecayGammaSpectrum::PutGgroup(int i)
{
  ggroup=i;
  ebnd.put_imax(ggroup+1);
};

void DecayGammaSpectrum::ReadFile(string cbglibdir, string filename)
{
  ifstream fin;
  string mdir=cbglibdir+"CBGLIB_BURN/DecayGammaSpectra/"+filename;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in DecayGammaSpectrum::ReadFile.\n";
    cout<<"# Failed to open the file ["<<mdir<<"].\n";
    exit(1);
  };
  int tmp;
  fin>>tmp;
  PutGgroup(tmp);
  for(int i=0;i<ggroup+1;i++){
    real tmp;
    fin>>tmp;
    ebnd.put_data(i,tmp);
  };

  int mat=0;
  while(mat!=-1){
    fin>>mat;
    if(mat!=-1){
      GroupData1D tmp(ggroup);
      for(int i=0;i<ggroup;i++){
	real val;
	fin>>val;
	tmp.put_data(i,val);
      };
      data[mat]=tmp;
    };
  };
};

bool DecayGammaSpectrum::ExistData(int mat)
{
  if(data.find(mat)!=data.end())return true;
  return false;
};

GroupData1D &DecayGammaSpectrum::GetData(int mat)
{
  if(ExistData(mat))return data[mat];

  cout<<"# Error in DecayGammaSpectrum::GetData.\n";
  cout<<"# No data for material ID ["<<mat<<"].\n";
  exit(0);
};

void DecayGammaSpectrum::ShowGammaSpectrum(GroupData1D &spec)
{
  cout<<"# Total gamma-ray emission spectrum\n";
  cout<<"#    [Up-Energy] [Spec]      [Spec/lethergy]\n";

  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<ggroup;i++){
    real e0=ebnd.get_dat(i);
    real e1=ebnd.get_dat(i+1);
    real letwid=log(e0/e1);
    WriteOut(i,4);
    cout<<" "<<e0<<" ";
    cout<<spec.get_dat(i)<<" ";
    cout<<spec.get_dat(i)/letwid;
    cout<<"\n";
  };
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GammaTool::GammaXSOUT(GammaXSLibrary &gxslib,Medium &med,string mdir,string ss,bool print)
{
  int ng=0;
  int gg=0;

  int nucnum=med.GetNucnum();

  for(int i=0;i<nucnum;i++){
    int matid=med.GetNuclideInTurn(i).GetMatnum();
    if(gxslib.ExistGammaLibData(matid)){
      ng=gxslib.GetGammaLibData(matid).GetNGroup();
      gg=gxslib.GetGammaLibData(matid).GetGGroup();
      break;
    };
  };

  // +++
  enum xstype g1dname[]={sigt,coh_el,incoh_el,pair_pro,siga,kerma};
  enum xstype g2dname[]={coh_el,incoh_el,pair_pro};

  vector<GroupData2D> yield(5);
  vector<GroupData1D> gamma1d(6);
  vector< vector<GroupData2D> > gamma2d(3);
  for(int i=0;i<5;i++){
    yield[i].put_yx(ng,gg);
    yield[i].set_zero();
  };
  for(int i=0;i<6;i++){
    gamma1d[i].put_imax(gg);
    gamma1d[i].set_zero();
  };
  for(int i=0;i<3;i++){
    gamma2d[i].resize(7);
    for(int j=0;j<7;j++){
      gamma2d[i][j].put_yx(gg,gg);
      gamma2d[i][j].set_zero();
    };
  };

  ofstream fout;
  
  vector< vector<real> > enuc(nucnum);
  for(int i=0;i<nucnum;i++){
    enuc[i].resize(5,0.);
  };

  if(print&&ng!=0){
    cout<<"#\n#   GAMMA ENERGY PER REACTION (MeV)\n";
    cout<<"#  MAT  DENSITY     (n,2n)    (n,f)     (n,c)     (n,n')    (non-e)\n#\n";
  };

  for(int i=0;i<nucnum;i++){
    int mat=med.GetNuclideInTurn(i).GetMatnum();
    real den=med.GetNuclideInTurn(i).GetDensity();
    int g_check=med.GetNuclideInTurn(i).GetGrp();
    if(g_check!=-1){
    cout.setf(ios::scientific);
    cout.precision(3);
    if(print&&ng!=0){
      WriteOut(mat,7);
      cout<<" "<<den<<" : ";
    };
    bool ext=gxslib.ExistGammaLibData(mat);
    if(!ext){
      if(print)cout<<" (no gamma data)\n";
    }else{
      vector<real> denom(5);
      vector<real> nume(5);
      for(int j=0;j<5;j++){
	denom[j]=0.;
	nume[j]=0.;
      };
      for(int g=0;g<ng;g++){
        vector<real> mic(5);
        mic[0]=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sign2n).get_dat(g);
        mic[1]=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigf).get_dat(g);
        mic[2]=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigc).get_dat(g);
        mic[3]=med.GetNuclideInTurn(i).GetMicxs().GetData1d(siginel).get_dat(g);
        real mice=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigel).get_dat(g);
        real mict=med.GetNuclideInTurn(i).GetMicxs().GetData1d(sigt).get_dat(g);
	mic[4]=mict-mice;
        real nflx=med.GetFlux().get_dat(g);
	for(int j=0;j<5;j++){
	  for(int k=0;k<gg;k++){
	    real tmp=gxslib.GetGammaLibData(mat).GetYield(j).get_dat(g,k);
	    yield[j].add_data(g,k,mic[j]*tmp*den);
            real aveng=gxslib.GetEnband().get_dat(k)*gxslib.GetEnband().get_dat(k+1);
	    aveng=sqrt(aveng);
	    nume[j]+=aveng*tmp*mic[j]*nflx;
	    if(k==0)denom[j]+=mic[j]*nflx;
	  };
	};
      };
      //
      for(int j=0;j<6;j++){
        for(int k=0;k<gg;k++){
          real xs=gxslib.GetGammaLibData(mat).GetXSData().GetData1d(g1dname[j]).get_dat(k);
	  gamma1d[j].add_data(k,xs*den);
        };
      };
      for(int j=0;j<3;j++){
        int mm=gxslib.GetGammaLibData(mat).GetXSData().GetDim2d(g2dname[j]);
        for(int m=0;m<mm;m++){
	  for(int k=0;k<gg;k++){
	    for(int l=0;l<gg;l++){
   	      real xs=gxslib.GetGammaLibData(mat).GetXSData().GetData2d(g2dname[j],m).get_dat(k,l);
	      gamma2d[j][m].add_data(k,l,xs*den);
	    };
	  };
        };
      };
      for(int j=0;j<5;j++){
        if(denom[j]!=0.){
  	  if(print&&ng!=0)cout<<nume[j]/denom[j]*1e-6<<" ";
	  enuc[i][j]=nume[j]*den*1e-6;
	}else{
	  if(print&&ng!=0)cout<<0.<<" ";
	};
      };
      if(print&&ng!=0)cout<<"\n";
      };
    };
  };

  real sum=0.;
  for(int i=0;i<nucnum;i++){
    for(int j=0;j<5;j++){
      sum+=enuc[i][j];
    };
  };
  if(print)cout<<"#\n# Ratio of gamma energy origen\n#\n";
  for(int i=0;i<nucnum;i++){
    real sump=0.;
    for(int j=0;j<5;j++){
      sump+=enuc[i][j];
    };
    if(sump>0.&&print){
      WriteOut(med.GetNuclideInTurn(i).GetMatnum(),7);
      cout<<" : ";
      cout<<sump/sum<<" ( ";
      for(int j=0;j<5;j++){
        cout<<enuc[i][j]/sump<<" ";
      };
      cout<<")\n";
    };
  };
  if(print)cout<<"#\n";

  // +++ Writing gamma yield data
  string yldname=mdir;
  yldname.append(ss);
  yldname.append(".yield");

  fout.open(yldname.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };
  fout.setf(ios::scientific);
  fout.precision(8);
  fout<<ng<<"\n";
  fout<<gg<<"\n";
  for(int i=0;i<ng;i++){
    for(int j=0;j<gg;j++){
      real sum=0.;
      for(int k=0;k<5;k++){
	sum+=yield[k].get_dat(i,j);
      };
      fout<<sum<<"\n";
    };
  };
  fout.close();

  // +++ Writing gamma cross section data as CBGXS format
  string gxsname=mdir;
  gxsname.append(ss);
  gxsname.append(".gxs");
  fout.open(gxsname.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };
  fout<<gg<<"\n";
  fout<<" 6\n";
  // (energy boundary)
  for(int i=0;i<gg+1;i++){
    fout<<gxslib.GetEnband().get_dat(i)<<"\n";
  };
  // (Kerma factor as fission yield data)
  for(int i=0;i<gg;i++){
    fout<<gamma1d[5].get_dat(i)<<"\n";
  };
  // (absorption)
  for(int i=0;i<gg;i++){
    fout<<gamma1d[4].get_dat(i)<<"\n";
  };
  // (dummy for chi)
  for(int i=0;i<gg;i++){
    fout<<" 0.\n";
  };
  // (pair-production as n2n)
  for(int i=0;i<gg;i++){
    fout<<gamma1d[3].get_dat(i)<<"\n";
  };
  // (total)
  for(int i=0;i<2;i++){
    for(int j=0;j<gg;j++){
      fout<<gamma1d[0].get_dat(j)<<"\n";
    };
  };
  // (dummy for D and flux/current)
  for(int i=0;i<3+2;i++){
    for(int j=0;j<gg;j++){
      fout<<" 0.\n";
    };
  };
  //
  for(int i=0;i<7;i++){
    for(int j=0;j<gg;j++){
      fout<<j<<"\n";
      fout<<gg-1<<"\n";
      for(int k=j;k<gg;k++){
	real sum=0.;
	for(int l=0;l<3;l++){
	  sum+=gamma2d[l][i].get_dat(j,k);
	};
	fout<<sum<<"\n";
      };
    };
  };

  fout<<" 0\n";

  fout.close();

  //
  cout.setf(ios::scientific);
  cout.precision(4);
  cout<<"# Macroscopic cross section for photon.\n";
  cout<<"#\n";
  cout<<"# Eng.     Total      Coh.       Inc.       P.P.       Abs.\n";
  for(int i=0;i<gg;i++){
    cout<<gxslib.GetEnband().get_dat(i)<<" ";
    cout<<gamma1d[0].get_dat(i)<<" "; // Total
    cout<<gamma1d[1].get_dat(i)<<" "; // Coherent-scattering
    cout<<gamma1d[2].get_dat(i)<<" "; // Incoherent-scattering
    cout<<gamma1d[3].get_dat(i)<<" "; // Pair-production
    cout<<gamma1d[4].get_dat(i)<<" "; // Absorption
    cout<<"\n";
  };
  cout<<"\n\n";
};

GroupData2D GammaTool::ReadYieldDataFromFile(string mdir,string filename)
{
  GroupData2D yld;

  ifstream fin;
  string mtmp=mdir;
  mtmp.append(filename);
  mtmp.append(".yield");
  fin.open(mtmp.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  int group,ggroup;
  fin>>group;
  fin>>ggroup;
  yld.put_yx(ggroup,group);

  for(int j=0;j<group;j++){
    for(int k=0;k<ggroup;k++){
      real tmp;
      fin>>tmp;
      yld.put_data(k,j,tmp);
    };
  };

  return yld;
};

void GammaTool::ReadDataFromFile(string mdir,string filename,Medium &med,GroupData2D &yld)
{
  // (Gamma transport data)
  string tmp=filename;
  tmp.append(".gxs");
  med.ReadFile(mdir,tmp);

  // (Gamma yield data)
  yld=ReadYieldDataFromFile(mdir,filename);
};

void GammaTool::HeatingCalculation(GeneralSystem &sys,GeneralSystem &gam,Burnup &bu,int *reg_med,bool decay)
{
  real ev_to_j=1.60219e-19;

  int mednum=sys.GetNmed();
  int group=sys.GetGrp();
  int totm=gam.GetTotM();

  GroupData1D enband=sys.GetMedium(0).GetEnband();

  for(int i=0;i<totm;i++){
    gam.GetMesh(i).CalFissionSrc();
  };

  // (Heating distribution)
  vector<real> med_fis(mednum,0);
  vector<real> med_gamma(mednum,0);
  vector<real> med_sct(mednum,0);
  vector<real> med_decay(mednum,0);
  cout<<"#\n";
  cout<<"#  Region-wise volume-integrated heat [J]\n";
  cout<<"#    (fission/gamma/elastic/volume/decay/total)\n#\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int m=0;m<totm;m++){
    real vol=gam.GetMesh(m).GetVolume();
    real fis=bu.GetIntegratedPowerParMesh(sys,m);
    real gamma=gam.GetMesh(m).GetFissionSrc()*ev_to_j;
    real sct=0.;
    GroupData2D macsigel=sys.GetMesh(m).GetMed()->GetMacro2D(sigel);
    for(int g=0;g<group;g++){
      real flx=sys.GetMesh(m).GetFlux().get_dat(g);
      real emed1=sqrt(enband.get_dat(g)*enband.get_dat(g+1));
      for(int g2=0;g2<group;g2++){
  	real emed2=sqrt(enband.get_dat(g2)*enband.get_dat(g2+1));
	sct+=flx*vol*(emed1-emed2)*macsigel.get_dat(g,g2)*ev_to_j;
      };
    };
    real dec=0.;
    if(decay){
    int nn=sys.GetMesh(m).GetMed()->GetNucnum();
    for(int i=0;i<nn;i++){
      int mat=sys.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetMatnum();
      real dc=bu.GetBurnupChain().GetDecayConstant(mat);
      if(dc!=0.){
        real den=sys.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetDensity();
        if(den>0.){
	  real coef=dc*den*1e24*vol*ev_to_j;
	  dec+=coef*(bu.GetBurnupChain().GetDecayEnergy(mat,0)+bu.GetBurnupChain().GetDecayEnergy(mat,2));
	};
      };
    };
    };
    WriteOut(m,3);
    cout<<" "<<fis<<" "<<gamma<<" "<<sct<<" "<<dec<<" "<<fis+gamma+sct+dec<<"\n";
    int mid=reg_med[m];
    med_fis[mid]+=fis;
    med_gamma[mid]+=gamma;
    med_sct[mid]+=sct;
    med_decay[mid]+=dec;
  };
  cout<<"#\n";


  cout<<"# Medium-wise volume-integrated heat [J]\n";
  cout<<"#   (fission/gamma/elastic/decay/total)\n#\n";
  real sum_fis=0.;
  real sum_gamma=0.;
  real sum_sct=0.;
  real sum_decay=0.;
  for(int m=0;m<mednum;m++){
    WriteOut(m,3);
    cout<<" "<<med_fis[m]<<" "<<med_gamma[m]<<" "<<med_sct[m]<<" "<<med_decay[m]<<" ";
    cout<<med_fis[m]+med_gamma[m]+med_sct[m]+med_decay[m]<<"\n";
    sum_sct+=med_sct[m];
    sum_fis+=med_fis[m];
    sum_gamma+=med_gamma[m];
    sum_decay+=med_decay[m];
  };
  cout<<"\n";


  cout.setf(ios::scientific);
  cout.precision(4);
  real sum=sum_fis+sum_gamma+sum_sct+sum_decay;
  cout<<"# Heat by fission (beta & FP) "<<sum_fis<<" ("<<sum_fis/sum<<")\n";
  cout<<"# Heat by gamma               "<<sum_gamma<<" ("<<sum_gamma/sum<<")\n";
  cout<<"# Heat by elastic scattering  "<<sum_sct<<" ("<<sum_sct/sum<<")\n";
  cout<<"# Heat by delayed decay       "<<sum_decay<<" ("<<sum_decay/sum<<")\n";
  cout<<"# Total power                 "<<sum<<"\n";

  int macro_med[]={
    0,1,0,1,0,1,0,1,0,1,
    0,1,0,1,0,1,0,1,0,1,
    0,1,0,1,0,1,0,1,0,1,
    0,1,0,1,0,1,0,1,0,1,
    2,2,2,2,2,2,2,2,3,2,
  };
  vector<real> macro_eng(4,0.);
  for(int m=0;m<mednum;m++){
    real tote=med_sct[m]+med_fis[m]+med_gamma[m]+med_decay[m];
    macro_eng[macro_med[m]]+=tote;
  };

  for(int i=0;i<4;i++){
    cout<<"# Heat of macro region "<<i<<" : "<<macro_eng[i]<<"\n";
  };


};

void GammaTool::GammaSourceCalculation(GeneralSystem &sys,GeneralSystem &gam,vector<GroupData2D> &yld,int *reg_med)
{
  int totm=gam.GetTotM();
  int ggroup=gam.GetGrp();
  int group=sys.GetGrp();

  for(int m=0;m<totm;m++){
    for(int g=0;g<ggroup;g++){
      real sum=0.;
      for(int ng=0;ng<group;ng++){
        real flx=sys.GetMesh(m).GetFlux().get_dat(ng);
        sum+=flx*yld[reg_med[m]].get_dat(g,ng);
      };
      gam.GetMesh(m).GetFlux().put_data(g,sum); // per unit volume
    };
  };
};

void GammaTool::AddGammaSourceFromDecay
(GeneralSystem &sys,GeneralSystem &gam,map<int,GroupData1D> &spec,Burnup &bu)
{
  int totm=gam.GetTotM();
  int ggroup=gam.GetGrp();
 
  for(int m=0;m<totm;m++){
    int nn=sys.GetMesh(m).GetMed()->GetNucnum();
    for(int i=0;i<nn;i++){
      int mat=sys.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetMatnum();
      real dc=bu.GetBurnupChain().GetDecayConstant(mat);
      if(dc!=0.){
        real den=sys.GetMesh(m).GetMed()->GetNuclideInTurn(i).GetDensity();
	if(den>0.){
          real coef=dc*den*1e24;
	  if(spec.find(mat)!=spec.end()){
	    for(int g=0;g<ggroup;g++){
              gam.GetMesh(m).GetFlux().add_data(g,coef*spec[mat].get_dat(g)); // per unit volume
	    };
	  };
	};
      };
    };
  };
};

void GammaTool::GammaTransportCalculation(GeneralSystem &gam)
{
  int totm=gam.GetTotM();
  int ggroup=gam.GetGrp();

  gam.SetZeroScatSrc();
  for(int i=0;i<totm;i++){
    GroupData1D src(ggroup);
    for(int g=0;g<ggroup;g++){
      real tmp=gam.GetMesh(i).GetFlux().get_dat(g);
      src.put_data(g,tmp);
    };
    gam.PutIsotropicSourceParVolume(i,src);
  };
  gam.CalFixedSource(1e-4);
};

void GammaTool::SetGammaEnergyAsProductionXS(GeneralSystem &gam)
{
  // gamma energy deposit without transport
  int ggroup=gam.GetGrp();
  int mednum=gam.GetNmed();
  for(int i=0;i<mednum;i++){
    real e=gam.GetMedium(i).GetEnband().get_dat(0);
    for(int g=0;g<ggroup;g++){
      real e1=gam.GetMedium(i).GetEnband().get_dat(g+1);
      e=sqrt(e*e1);
      gam.GetMed(i).GetMacxs().GetData1d(nusigf).put_data(g,e);
      e=e1;
    };
  };
};
