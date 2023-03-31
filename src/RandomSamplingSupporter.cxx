#include<cstdlib>
#include "RandomSamplingSupporter.h"

using namespace std;

void RandomSamplingSupporter::ReadXSLibrary(int nucnum,string libdir,string *fname)
{
  xslib.Initialize(libdir,"N-ENERGY");
  xslib.ReadFile(nucnum,libdir,fname);
  grp=xslib.GetGroup();
};

void RandomSamplingSupporter::ReadCovarianceFromFile(int nucnum,string covdir,string *fname)
{
  for(int i=0;i<nucnum;i++){
    xscov.ReadCovarianceFromFile(covdir,fname[i]);
  };
};

void RandomSamplingSupporter::PutSampledNuclearDataInfo(int i1,int *i2, int *i3)
{
  rs_nuc_data=i1;

  rs_matid.resize(rs_nuc_data);
  rs_mtid.resize(rs_nuc_data);
  for(int i=0;i<rs_nuc_data;i++){
    rs_matid[i]=i2[i];
    rs_mtid[i]=i3[i];
  };
};

void RandomSamplingSupporter::PutSampleNum(int i)
{
  sample=i;

  if(rs_nuc_data==0.){
    cout<<"# Error in RandomSamplingSupporter::PutSampleNum.\n";
    cout<<"# Please do [PutSampledNuclearDataInfo] in prior.\n";
    exit(0);
  };

  sampled_data.resize(rs_nuc_data, (vector<GroupData1D>) sample);
};

void RandomSamplingSupporter::PrepareSamplingData(int i1,real factor)
{
  PutSampleNum(i1);
  for(int i=0;i<rs_nuc_data;i++){
    int mat=rs_matid[i];
    int mt=rs_mtid[i];
    CrossSectionCovariance csv=xscov.GetCrossSectionCovariance(mat,mt,mat,mt);
    csv.DoRandomSampling(sample,sampled_data[i],"Relative");
  };

  // Calculation sample average
  SampleAverageCalculation(factor);
};

void RandomSamplingSupporter::SampleAverageCalculation(real factor)
{
  sampled_avg.resize(rs_nuc_data); // relative perturbation
  for(int i=0;i<rs_nuc_data;i++){
    GroupData1D sum(grp);
    sum.set_zero();
    for(int j=0;j<sample;j++){
      sum=sum+sampled_data[i][j]*factor;
    };
    real tmp=1./sample;
    sum=sum*tmp;
    sampled_avg[i]=sum;
  };

};

void RandomSamplingSupporter::PrepareOriginalData()
{
  orgxs.resize(rs_nuc_data, (vector<real>) grp);
  for(int i=0;i<rs_nuc_data;i++){
    int mat=rs_matid[i];
    int mt=rs_mtid[i];
    for(int g=0;g<grp;g++){
      orgxs[i][g]=xslib.GetLibData(mat).GetXSData().GetData1d(XSType(mt)).get_dat(g);
    };
  };
};

void RandomSamplingSupporter::PerturbXSLib(int id,real factor)
{
  for(int i=0;i<rs_nuc_data;i++){
    int mat=rs_matid[i];
    int mt=rs_mtid[i];
    for(int g=0;g<grp;g++){
      real tmp=sampled_data[i][id].get_dat(g)*factor*orgxs[i][g];
      xslib.GetLibData(mat).GetXSData().GetData1d(XSType(mt)).add_data(g,tmp);
      if(mt==2||mt==4||mt==18||mt==102){
        xslib.GetLibData(mat).GetXSData().GetData1d(sigt).add_data(g,tmp);
      };
      //dxs[i][g]=tmp;
    };
    if(mt==2||mt==4)xslib.GetLibData(mat).GetXSData().XSMatrixNormalizationTo1DXS();
    if(mt==251)xslib.GetLibData(mat).GetXSData().P1MatrixNormalizationToMubar();
  };
};

void RandomSamplingSupporter::PutOutputNum(int i)
{
  outp_data.clear();
  output_avg.clear();
  outp_num=i;
  outp_data.resize(i);
  output_avg.resize(i);
  //outp_data.resize(i,(vector<real>) sample);
  //output_avg.resize(i);
};

void RandomSamplingSupporter::PutOutputData(int i,int j,real k)
{
  outp_data[i][j]=k;
};

void RandomSamplingSupporter::PutOutputData(int i,real k)
{
  outp_data[i].push_back(k);
};

void RandomSamplingSupporter::ShowOutputStatics(int i)
{
  int sz=outp_data[i].size();
  GetStatics(outp_data[i],sz,true);
};

void RandomSamplingSupporter::CalOutputAverage()
{
  for(int i=0;i<outp_num;i++){
    real sum=0.;
    int sz=outp_data[i].size();
    for(int j=0;j<sz;j++){
      sum+=outp_data[i][j];
    };
    output_avg[i]=sum/sz;
  };
};

SensitivityData RandomSamplingSupporter::SensitivityCalculation(real red_factor)
{
  CalOutputAverage();

  int outp_sz=outp_data[0].size();
  if(outp_sz!=sample){
    cout<<"# Error in RandomSamplingSupporter::SensitivityCalculation.\n";
    cout<<"# Inconsistency between the number of XS samples and that of output samples.\n";
    cout<<"# The number of XS samples     : "<<sample<<"\n";
    cout<<"# The number of output samples : "<<outp_sz<<"\n";
    exit(0);
  };

  // +++ Sample relative covariance
  int size1=rs_nuc_data*grp;
  GroupData2D covdat(size1,size1); // sigma-sigma
  GroupData1D cov_oi(size1);       // sigma-output

  int ii=0;
  for(int i1=0;i1<rs_nuc_data;i1++){
    for(int i2=0;i2<grp;i2++){

      if(orgxs[i1][i2]!=0.){

      real avgxs1_sample=orgxs[i1][i2]*(1.+sampled_avg[i1].get_dat(i2));
      real alpha1=orgxs[i1][i2]/avgxs1_sample;
      real beta1=sampled_avg[i1].get_dat(i2);
      real sum=0.;
      for(int k=0;k<sample;k++){
	sum+=(outp_data[0][k]-output_avg[0])/output_avg[0]
	  *(sampled_data[i1][k].get_dat(i2)*red_factor-beta1)*alpha1;
      };
      sum/=(sample-1);
      cov_oi.put_data(ii,sum);
      int jj=0;
      for(int j1=0;j1<rs_nuc_data;j1++){
	for(int j2=0;j2<grp;j2++){

          if(orgxs[j1][j2]!=0.){

          real avgxs2_sample=orgxs[j1][j2]*(1.+sampled_avg[j1].get_dat(j2));
          real alpha2=orgxs[j1][j2]/avgxs2_sample;
          real beta2=sampled_avg[j1].get_dat(j2);
	  real sum=0.;
	  for(int k=0;k<sample;k++){
	    sum+=(sampled_data[i1][k].get_dat(i2)*red_factor-beta1)*alpha1
  	        *(sampled_data[j1][k].get_dat(j2)*red_factor-beta2)*alpha2;
	  };
	  sum/=(sample-1);
	  covdat.put_data(ii,jj,sum);

	  };

	  jj++;
	};
      };

      };

      ii++;
    };
  };

  // True covariance is used for cross section
  /*
  covdat.set_zero();
  for(int i=0;i<rs_nuc_data;i++){
    int mat=rs_matid[i];
    int mt=rs_mtid[i];
    //CrossSectionCovariance csv=xscov.GetCrossSectionCovariance(mat,mt,mat,mt);
    GroupData2D csvdat=xscov.GetCrossSectionCovariance(mat,mt,mat,mt).GetCovariance("Relative");
    int tmp=i*grp;
    for(int g1=0;g1<grp;g1++){
      for(int g2=0;g2<grp;g2++){
        covdat.put_data(tmp+g1,tmp+g2,csvdat.get_dat(g1,g2));
      };
    };
  };
  */

  cout<<"# Matrix diagonalization ...\n";
  // +++ matrix diagonalization
  GroupData2D diagmat2(size1,size1);
  GroupData2D evecmat2(size1,size1);
  covdat.Diagonalization(diagmat2,evecmat2);
  for(int g=0;g<size1;g++){
    real org=diagmat2.get_dat(g,g);
    if(org<1e-10)org=0.;
    if(org!=0.)org=1./org;
    diagmat2.put_data(g,g,org);
  };
  GroupData2D tmpmat=evecmat2.GetTransposedMatrix();
  cov_oi=evecmat2*diagmat2*tmpmat*cov_oi;

  SensitivityData snsd;
  snsd.PutGroup(grp);
  snsd.GetEnband()=xslib.GetEnband();
  for(int i1=0;i1<rs_nuc_data;i1++){
    int mat=rs_matid[i1];
    int mt=rs_mtid[i1];
    GroupData1D dat(grp);
    for(int i2=0;i2<grp;i2++){
      dat.put_data(i2,cov_oi.get_dat(i1*grp+i2));
    };
    snsd.PutSensitivity1D(mat,mt,dat);
  };

  return snsd;

  /*
  for(int g=0;g<grp;g++){
    real e0=xslib.GetEnband().get_dat(g);
    real e1=xslib.GetEnband().get_dat(g+1);
    real letwid=log(e0/e1);
    cout<<e0<<" ";
    for(int i=0;i<rs_nuc_data;i++){
      cout<<cov_oi.get_dat(i*grp+g)/letwid*0.25<<" ";
    };
    cout<<"\n";
  };
  */
  

  /*
  for(int i=0;i<rs_nuc_data;i++){

    GroupData2D covp(grp,grp);
    for(int g=0;g<grp;g++){
      for(int g2=0;g2<grp;g2++){
	covp.put_data(g,g2,covdat.get_dat(i*grp+g,i*grp+g2));
      };
    };

    GroupData1D covp_oi(grp);
    for(int g=0;g<grp;g++){
      covp_oi.put_data(g,cov_oi.get_dat(i*grp+g));
    };

    // +++ matrix diagonalization
    GroupData2D diagmat2(grp,grp);
    GroupData2D evecmat2(grp,grp);
    covp.Diagonalization(diagmat2,evecmat2);
    for(int g=0;g<grp;g++){
      real org=diagmat2.get_dat(g,g);
      if(org<1e-10)org=0.;
      if(org!=0.)org=1./org;
      diagmat2.put_data(g,g,org);
    };
    GroupData2D tmpmat=evecmat2.GetTransposedMatrix();
    covp_oi=evecmat2*diagmat2*tmpmat*covp_oi;

    for(int g=0;g<grp;g++){
      real e0=xslib.GetEnband().get_dat(g);
      real e1=xslib.GetEnband().get_dat(g+1);
      real letwid=log(e0/e1);
      cout<<e0<<" ";
      cout<<covp_oi.get_dat(g)/letwid*0.25<<" ";
      cout<<"\n";
    };

    cout<<"\n\n";
  };

  */
};



void RandomSamplingSupporter::OutputPartialUncertaintyCalculation()
{
  CalOutputAverage();

  /*
  int outp_sz=outp_data[0].size();
  if(outp_sz!=sample){
    cout<<"# Error in RandomSamplingSupporter::SensitivityCalculation.\n";
    cout<<"# Inconsistency between the number of XS samples and that of output samples.\n";
    cout<<"# The number of XS samples     : "<<sample<<"\n";
    cout<<"# The number of output samples : "<<outp_sz<<"\n";
    exit(0);
  };
  */

  // +++ Sample relative covariance
  GroupData1D cov_oi(grp);       // sigma-output
  for(int i1=0;i1<rs_nuc_data;i1++){
    int mat=rs_matid[i1];
    int mt=rs_mtid[i1];

    for(int i2=0;i2<grp;i2++){
      if(orgxs[i1][i2]!=0.){
        real avgxs1_sample=orgxs[i1][i2]*(1.+sampled_avg[i1].get_dat(i2));
        real alpha1=orgxs[i1][i2]/avgxs1_sample;
        real beta1=sampled_avg[i1].get_dat(i2);
        real sum=0.;
        for(int k=0;k<sample;k++){
	  sum+=(outp_data[0][k]-output_avg[0])/output_avg[0]
	    *(sampled_data[i1][k].get_dat(i2)-beta1)*alpha1;
        };
        sum/=(sample-1);
        cov_oi.put_data(i2,sum);
      };
    };

    GroupData2D csvdat=xscov.GetCrossSectionCovariance(mat,mt,mat,mt).GetCovariance("Relative");
    // +++ matrix diagonalization
    GroupData2D diagmat2(grp,grp);
    GroupData2D evecmat2(grp,grp);
    csvdat.Diagonalization(diagmat2,evecmat2);
    for(int g=0;g<grp;g++){
      real org=diagmat2.get_dat(g,g);
      if(org<1e-10)org=0.;
      if(org!=0.)org=1./org;
      diagmat2.put_data(g,g,org);
    };
    GroupData2D tmpmat=evecmat2.GetTransposedMatrix();
    tmpmat=evecmat2*diagmat2*tmpmat;
    real val=cov_oi*(tmpmat*cov_oi);
    cout<<mat<<" "<<mt<<" "<<val<<" "<<sqrt(val)<<"\n";

  };
};



