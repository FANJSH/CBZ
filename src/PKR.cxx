#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "PKR.h"

using namespace std;

PromptKratio::PromptKratio
()
{
};

void PromptKratio::DelayedNeutronDataExtractionFromChi
(XSLibrary& xslib, DelayedNeutronData& dnd,int matid)
{
  if(!xslib.ExistLibData(matid)){
    cout<<"# Warning in PromptKratio::DelayedNeutronDataExtractionFromChi.\n";
    cout<<"# XSLibrary does NOT contain nuclide "<<matid<<"\n";
  };

  int group=xslib.GetLibData(matid).GetGroup();
  xslib.GetLibData(matid).GetChiVector().AssignChiVectorToEachGroup();

  int dn=dnd.GetNumFromNucid(matid);
  if(dn==-1){
    cout<<"# Error in PromptKratio::DelayedNeutronDataExtractionFromChi.\n";
    cout<<"# DelayedNeutronData does NOT include the material ID : "<<matid<<"\n";
    exit(0);
  };

  GroupData1D chid=dnd.GetAveragedChi(matid);
  for(int i=0;i<group;i++){
    real totnu=xslib.GetLibData(matid).GetXSData().GetData1d(nu).get_dat(i);
    real delay=dnd.GetYield(dn).get_dat(i);
    real propt=totnu-delay;
    xslib.GetLibData(matid).GetXSData().GetData1d(nu).put_data(i,propt);
    GroupData1D tmp;
    tmp.put_imax(group);
    tmp=xslib.GetLibData(matid).GetChiVector().GetChiVectorGroup(i);
    tmp=tmp*totnu-chid*delay;
    real ttt=1./propt;
    tmp=tmp*ttt;
    xslib.GetLibData(matid).GetChiVector().GetChiVectorGroup(i).copy(tmp);
  };
};

void PromptKratio::DelayedNeutronDataAdditionToChi
(XSLibrary& xslib, DelayedNeutronData& dnd,int matid, real factor)
{
  if(!xslib.ExistLibData(matid)){
    cout<<"#\n# Warning in PromptKratio::DelayedNeutronDataAdditionToChi.\n";
    cout<<"# XSLibrary does NOT contain nuclide "<<matid<<"\n#\n";
    return;
  };

  int group=xslib.GetLibData(matid).GetGroup();
  xslib.GetLibData(matid).GetChiVector().AssignChiVectorToEachGroup();

  int dn=dnd.GetNumFromNucid(matid);
  if(dn==-1){
    cout<<"# Error in PromptKratio::DelayedNeutronDataAdditionToChi.\n";
    cout<<"# DelayedNeutronData does NOT include the material ID : "<<matid<<"\n";
    exit(0);
  };

  GroupData1D chid=dnd.GetAveragedChi(matid);
  for(int i=0;i<group;i++){
    real totnu=xslib.GetLibData(matid).GetXSData().GetData1d(nu).get_dat(i);
    real delay=dnd.GetYield(dn).get_dat(i);
    real newnu=totnu+delay*factor;
    xslib.GetLibData(matid).GetXSData().GetData1d(nu).put_data(i,newnu);
    GroupData1D tmp;
    tmp.put_imax(group);
    tmp=xslib.GetLibData(matid).GetChiVector().GetChiVectorGroup(i);
    tmp=tmp*totnu+chid*delay*factor;
    real ttt=1./newnu;
    tmp=tmp*ttt;
    xslib.GetLibData(matid).GetChiVector().GetChiVectorGroup(i).copy(tmp);
  };
};

SensitivityData PromptKratio::CalSensitivity
//    (SNRSystem* fwd, SNRSystem* adj, XSLibrary &xslib, DelayedNeutronData& dnd, 
    (GeneralSystem* fwd, GeneralSystem* adj, XSLibrary &xslib, DelayedNeutronData& dnd, 
     real keff, int nucnum, int* nucid, real factor)
{
  // - Fission spectrum matrix is EXPLICITLY considered.
  // - Only fissile nuclides are treated.
  //
  // (factor) delayed neutron data addition 

  /*
  // +++ Fissile nuclide checking
  for(int i=0;i<nucnum;i++){
    int id=nucid[i];
    if(!xslib.GetLibData(id).fissile()){
      cout<<"# Error in PromptKratio::CalSensitivity.\n";
      cout<<"# Nuclide ID "<<id<<" is NOT fissile.\n";
      exit(0);
    };
  };
  */

  // +++ Forward/Adjoint checking
  if(!fwd->GetGeneralOption().Forward()||
      adj->GetGeneralOption().Forward()){
    cout<<"# Error in PromptKratio::CalSensitivity.\n";
    cout<<"# Please check forward/adjoint calcution.\n";
    exit(0);
  };

  // +++ Mesh consistency checking
  fwd->CheckSameMesh(adj);
  
  // +++ Perturbation denominator calculation
  real ip=adj->CalPerturbDenominatorWithFissionSpectrumMatrix(fwd);
  if(ip==0.){
    cout<<"# Error in PromptKratio::CalSensitivity.\n";
    cout<<"# Perturbation denominator is zero.\n";
    exit(0);
  };

  int nmed=fwd->GetNmed();
  int grp=fwd->GetGrp();
  int TotM=fwd->GetTotM();

  real *nsforg=new real[nmed];
  vector<GroupData1D> chimat_org(nmed);
  for(int i=0;i<nmed;i++){
    chimat_org[i].put_imax(grp);
  };

  cout.setf(ios::scientific);
  cout.precision(5);

  bool *flag=new bool[TotM];
  for(int i=0;i<TotM;i++){
    flag[i]=true;
  };

  SensitivityData sens;
  sens.PutValue(keff);
  sens.PutGroup(grp);
  sens.GetEnband().copy(fwd->GetMed(0).GetEnband());

  GroupData1D sns1d(grp);

  for(int nc=0;nc<nucnum;nc++){

    int nid=nucid[nc];
    bool nuclide_in_system=false;
    for(int i=0;i<TotM;i++){
      if(fwd->GetMesh(i).GetMed()->ExistNuclide(nid)){
        flag[i]=true;
	nuclide_in_system=true;
      }else{
	flag[i]=false;
      };
    };

    // +++ Delayed neutron data
    int dn=dnd.GetNumFromNucid(nid);
    /*
    if(dn==-1){
      cout<<"# Error in PromptKratio::CalSensitivity.\n";
      cout<<"# DelayedNeutronData does NOT include the material ID : "<<nid<<"\n";
      exit(0);
    };
    */

    if(nuclide_in_system&&dn!=-1){

      GroupData1D chid=dnd.GetAveragedChi(nid);

      for(int it=0;it<2;it++){
	// it=0:Nu_p(mt=452) / it=1:Nu_d(mt=455)
        int mtnum=452;
	if(it==1)mtnum=455;
        for(int i=0;i<grp;i++){
          real micnu_d=dnd.GetYield(dn).get_dat(i)*(1.+factor);
          for(int j=0;j<nmed;j++){
  	    if(fwd->GetMedium(j).ExistNuclide(nid)){
              nsforg[j]=fwd->GetMed(j).GetMacxs().GetNusigf().get_dat(i);
              real den=fwd->GetMedium(j).GetNuclide(nid).GetDensity();
	      real micsigf=fwd->GetMed(j).GetNuclide(nid).GetMicxs().GetMicSigf().get_dat(i);
	      real micnu_t=fwd->GetMed(j).GetNuclide(nid).GetMicxs().GetMicNu().get_dat(i);
              real micnu_p=micnu_t-micnu_d;
              real del_nu=micnu_p;
              if(it==1)del_nu=micnu_d;
              fwd->GetMed(j).GetMacxs().GetNusigf().add_data(i,den*del_nu*micsigf);
              // (re-calculation of chi matrix)
	      for(int ii=0;ii<grp;ii++){
	        real org_mac_chi=fwd->GetMed(j).GetMacxs().GetData2d(chi).get_dat(i,ii);
                chimat_org[j].put_data(ii,org_mac_chi);
                real micchi=0.;
                if(it==0){
                  int vecid=xslib.GetLibData(nid).GetChiVector().GetVectorID(i);
                  micchi=xslib.GetLibData(nid).GetChiVector(vecid).get_dat(ii);
  	  	  // --- prompt chi is assumed to be same as total chi.
		}else{
		  micchi=chid.get_dat(ii);
		};
                real newchi=(org_mac_chi*nsforg[j]+den*del_nu*micsigf*micchi);
	        newchi/=nsforg[j]+den*del_nu*micsigf;
                fwd->GetMed(j).GetMacxs().GetData2d(chi).put_data(i,ii,newchi);
	      };
	    };
          };
          real re=adj->CalPerturbYieldTermWithFissionSpectrumMatrix(fwd,flag,i,keff)/ip;
          //cout<<"  "<<re<<"\n";
	  sns1d.put_data(i,re);
          for(int j=0;j<nmed;j++){
            if(fwd->GetMed(j).ExistNuclide(nid)){
	      fwd->GetMed(j).GetMacxs().GetNusigf().put_data(i,nsforg[j]);
              for(int ii=0;ii<grp;ii++){
                fwd->GetMed(j).GetMacxs().GetData2d(chi).put_data(i,ii,chimat_org[j].get_dat(ii));
	      };
	    };
          };
        };
        sens.PutSensitivity1D(nid,mtnum,sns1d);
      };
    };
  };

  delete [] flag;
  delete [] nsforg;

  return sens;
};
