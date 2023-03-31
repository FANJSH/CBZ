#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "YieldDecayCovariance.h"

using namespace std;

//++++++(IndependentYieldCovariance)++++++

void IndependentYieldCovariance::SetDefaultCovariance(bool corr_mass_chain)
{
  BCGManager ibm;

  int fnuc=23;
  string nucname_fnuc[]={
    "Th232","Pa231","U232","U233","U234",
    "U235","U236","U237","U238","Np237",
    "Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242m","Am243","Cm242","Cm243",
    "Cm244","Cm245","Cm246"
  };
  int einc[]={0,0,0,0,0, 0,0,0,0,1, 0,0,1,0,1, 0,0,0,0,0, 0,0,0}; // thermal 

  for(int i=0;i<fnuc;i++){
    ReadYieldDataFromFile(nucname_fnuc[i],"JENDL-FPY-2011/"+nucname_fnuc[i],einc[i]);
  };

  ibm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.endf71");
  ibm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.jendl2011");

  if(corr_mass_chain){
    MakeCovarianceMatrix(ibm);
  }else{
    MakeCovarianceMatrixWithoutCorrelation();
  };
};

void IndependentYieldCovariance::SetCovariance(int fnuc, string *nucname_fnuc, int *einc,bool corr_mass_chain)
{
  BCGManager ibm;

  /*
  int fnuc=23;
  string nucname_fnuc[]={
    "Th232","Pa231","U232","U233","U234",
    "U235","U236","U237","U238","Np237",
    "Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242m","Am243","Cm242","Cm243",
    "Cm244","Cm245","Cm246"
  };
  int einc[]={0,0,0,0,0, 0,0,0,0,1, 0,0,1,0,1, 0,0,0,0,0, 0,0,0}; // thermal 
  */

  for(int i=0;i<fnuc;i++){
    ReadYieldDataFromFile(nucname_fnuc[i],"JENDL-FPY-2011/"+nucname_fnuc[i],einc[i]);
  };

  ibm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.endf71");
  ibm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.jendl2011");

  if(corr_mass_chain){
    MakeCovarianceMatrix(ibm);
  }else{
    MakeCovarianceMatrixWithoutCorrelation();
  };
};

void IndependentYieldCovariance::ReadYieldDataFromFile(string fisnucname, string filename, int epnt)
{
  int sz=fisnucname_array.size();
  for(int i=0;i<sz;i++){
    if(fisnucname==fisnucname_array[i]){
      cout<<"# Error in IndependentYieldCovariance::ReadYieldDataFromFile.\n";
      cout<<"# Fission yield covariance already exists for the name "<<fisnucname<<"\n";
      exit(0);
    };
  };

  fisnucnum++;
  fisnucname_array.push_back(fisnucname);
  filename_array.push_back(filename);
  energy_point_array.push_back(epnt);
};

void IndependentYieldCovariance::MakeCovarianceMatrix(BCGManager &bm)
{
  id.resize(fisnucnum);
  yield_idp.resize(fisnucnum);
  yield_cum.resize(fisnucnum);
  unc_idp.resize(fisnucnum);
  unc_cum.resize(fisnucnum);
  energy_array.resize(fisnucnum);
  size_array.resize(fisnucnum);
  matrix.resize(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    ReadDataFromFile(i);
    CalCovariance(bm,i);
  };
};

void IndependentYieldCovariance::MakeCovarianceMatrixGLSUpdating(BCGManager &bm)
{
  bool condition[]={true,true,true,true,true};
  MakeCovarianceMatrixGLSUpdating(bm,condition);
};

void IndependentYieldCovariance::MakeCovarianceMatrixGLSUpdating(BCGManager &bm, bool *condition)
{
  id.resize(fisnucnum);
  yield_idp.resize(fisnucnum);
  yield_cum.resize(fisnucnum);
  unc_idp.resize(fisnucnum);
  unc_cum.resize(fisnucnum);
  energy_array.resize(fisnucnum);
  size_array.resize(fisnucnum);
  matrix.resize(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    ReadDataFromFile(i);
    CalCovarianceGLSUpdating(bm,i,condition);
  };
};

void IndependentYieldCovariance::MakeCovarianceMatrixWithoutCorrelation(){
  id.resize(fisnucnum);
  yield_idp.resize(fisnucnum);
  yield_cum.resize(fisnucnum);
  unc_idp.resize(fisnucnum);
  unc_cum.resize(fisnucnum);
  energy_array.resize(fisnucnum);
  size_array.resize(fisnucnum);
  matrix.resize(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    ReadDataFromFile(i);
    CalCovarianceWithoutCorrelation(i);
  };
};

void IndependentYieldCovariance::ReadDataFromFile(int num){
  vector<int> id_tmp;
  vector<real> yield_idp_tmp;
  vector<real> yield_cum_tmp;
  vector<real> unc_idp_tmp;
  vector<real> unc_cum_tmp;

  string filename=filename_array[num];
  int energy_point=energy_point_array[num];
  string file_idp="../../CBGLIB_BURN/FPYfile/"+filename;
  string file_cum="../../CBGLIB_BURN/FPYfile_cum/"+filename;
  ifstream fin_idp;
  ifstream fin_cum;
  fin_idp.open(file_idp.data(),ios::in);
  fin_cum.open(file_cum.data(),ios::in);
  if(fin_idp.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<file_idp<<"\n";
    exit(0);
  };
  if(fin_cum.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<file_cum<<"\n";
    exit(0);
  };

  int epnt_num;
  fin_idp>>epnt_num;
  fin_cum>>epnt_num;

  for(int i=0;i<epnt_num;i++){
    real energy_tmp;
    int size_tmp;

    fin_idp>>energy_tmp;
    fin_cum>>energy_tmp;
    fin_idp>>size_tmp;
    fin_cum>>size_tmp;

    id[num].resize(size_tmp);
    yield_idp[num].resize(size_tmp);
    yield_cum[num].resize(size_tmp);
    unc_idp[num].resize(size_tmp);
    unc_cum[num].resize(size_tmp);

    id_tmp.resize(size_tmp);
    yield_idp_tmp.resize(size_tmp);
    yield_cum_tmp.resize(size_tmp);
    unc_idp_tmp.resize(size_tmp);
    unc_cum_tmp.resize(size_tmp);

    energy_array[num]=energy_tmp;
    size_array[num]=size_tmp;
    matrix[num].put_yx(size_tmp,size_tmp);

    for(int j=0;j<size_tmp;j++){
      int iz,ia,il;
      fin_idp>>iz;
      fin_cum>>iz;
      fin_idp>>ia;
      fin_cum>>ia;
      fin_idp>>il;
      fin_cum>>il;
      real yc_idp,yc_cum,dyc_idp,dyc_cum;
      fin_idp>>yc_idp;
      fin_cum>>yc_cum;
      fin_idp>>dyc_idp;
      fin_cum>>dyc_cum;
      if(i==energy_point){
	size_array[num]=size_tmp;
	id_tmp[j]=eidt.NewID(iz,ia,il);
	yield_idp_tmp[j]=yc_idp;
	yield_cum_tmp[j]=yc_cum;
	unc_idp_tmp[j]=dyc_idp;
	unc_cum_tmp[j]=dyc_cum;
      };
    };
  };

  int counter=0;
  for(int i=0;i<999992;i++){
    for(int j=0;j<size_array[num];j++){
      if(i==id_tmp[j]){
	id[num][counter]=id_tmp[j];
	yield_idp[num][counter]=yield_idp_tmp[j];
	yield_cum[num][counter]=yield_cum_tmp[j];
	//If the error is zero even though the yield is not zero, we use assumed value (dafault; 100%).
	if(unc_idp_tmp[j]==0.&&yield_idp_tmp[j]!=0.){
	  unc_idp[num][counter]=yield_idp_tmp[j]*theor_err;//assumed
	}else{
	  unc_idp[num][counter]=unc_idp_tmp[j];
	};
	if(unc_cum_tmp[j]==0.&&yield_cum_tmp[j]!=0.){
	  unc_cum[num][counter]=yield_cum_tmp[j]*theor_err;//assumed
	}else{
	  unc_cum[num][counter]=unc_cum_tmp[j];
	};
	//+++++++++++++++
	counter++;
      };
    };
  };
};
/*
void IndependentYieldCovariance::CalCovariance(BCGManager &bm,int num){

  int z;
  int a;
  int l;
  matrix[num].set_zero();
  
  for(int i=1;i<999;i++){

    //+++Searching Mass Yield+++    
    vector<real> mass_unc;
    vector<int> mass_nuc;
    vector<int> mass_z;
    vector<real> sum_unc_square;
    int mass_nuc_num=0;
    vector<int> mass_chain_id;//1,2,3,...
    mass_chain_id.resize(size_array[num]);
    for(int j=0;j<size_array[num];j++){
      mass_chain_id[j]=0.;
    };

    for(int j=0;j<size_array[num];j++){
      eidt.GetParameterNew(id[num][j],z,a,l);
      if(a==i){
	int channel=bm.GetNuclide(z,a,l).GetChannel();
	real hl=bm.GetNuclide(z,a,l).GetHalflife();
	if(hl==0.){
	  mass_unc.push_back(unc_cum[num][j]);
	  mass_nuc.push_back(id[num][j]);
	  mass_z.push_back(z);
	  mass_nuc_num++;
	  mass_chain_id[j]=mass_nuc_num;
	}else if(channel==1){
	  int decay_type=bm.GetNuclide(z,a,l).GetDecayType(0);
	  if(decay_type==3){
	    mass_unc.push_back(unc_cum[num][j]);
	    mass_nuc.push_back(id[num][j]);
	    mass_z.push_back(z);
	    mass_nuc_num++;
	    mass_chain_id[j]=mass_nuc_num;
	  };
	};
      };
    };

    //mass_nuc_num=1;

    real tmp2;
    for(int jj=0;jj<mass_nuc_num-1;jj++){
      tmp2=0.;
      for(int j=0;j<size_array[num];j++){
	eidt.GetParameterNew(id[num][j],z,a,l);
	if(a==i&&mass_z[jj]>z&&mass_chain_id[j]==0){
	  mass_chain_id[j]=jj+1;
	  tmp2+=pow(unc_idp[num][j],2.);
	};	
      };
      sum_unc_square.push_back(tmp2);
      };

      tmp2=0.;
      for(int j=0;j<size_array[num];j++){
      eidt.GetParameterNew(id[num][j],z,a,l);
      if(a==i&&mass_chain_id[j]==0){
      mass_chain_id[j]=mass_nuc_num;
      tmp2+=pow(unc_idp[num][j],2.);
      };
      };
      sum_unc_square.push_back(tmp2);

      //++++++++++++++++++++++++++
  
      for(int ii=0;ii<mass_nuc_num;ii++){
    
      real tmp;
      real tmp1;
     
      for(int k=0;k<size_array[num];k++){
      eidt.GetParameterNew(id[num][k],z,a,l);
      if(a==i){
      for(int kk=0;kk<size_array[num];kk++){
      eidt.GetParameterNew(id[num][kk],z,a,l);
      if(a==i&&kk==k&&mass_chain_id[k]==ii+1){
      tmp1=yield_idp[num][k]*yield_idp[num][kk];
      if(tmp1==0.){
      tmp=0.;
      }else{
      tmp=pow(unc_idp[num][k],2.)*(1.-(pow(unc_idp[num][k],2.))/(pow(mass_unc[ii],2.)+sum_unc_square[ii]))/tmp1;
      if(tmp<0.)cout<<"!!!:"<<pow(unc_idp[num][k],2.)<<" "<<sum_unc_square[ii]<<"\n";
      };
      matrix[num].put_data(k,kk,tmp);
      }else if(a==i&&kk!=k&&mass_chain_id[k]==ii+1&&mass_chain_id[kk]==ii+1){
      tmp1=yield_idp[num][k]*yield_idp[num][kk];
      if(tmp1==0.){
      tmp=0.;
      }else{
      tmp=-(pow(unc_idp[num][k],2.)*pow(unc_idp[num][kk],2.))/(pow(mass_unc[ii],2.)+sum_unc_square[ii])/tmp1;
      };
      matrix[num].put_data(k,kk,tmp);
      };
      };
      };
      };
      };
      };
      };
*/

void IndependentYieldCovariance::CalCovariance(BCGManager &bm,int num)
{
  // Condition of the mass yield search
  //
  //  - Its half-life is zero.  (=stable nuclide)
  //  - It has only one decay channel, and it is the alpha-decay.
  //  - Its half-life is longer than 1e20 sec = 3.1e12 year.
  //
  //  The third condition is added in 2022/9/13
  //  since non-zero half-life is given to Nd-150 in JENDL-5 and ENDF/B-VIII.0.
  //  Its half-life is around 2.5e26 sec, so it is regarded as unstable nuclide.
  //  In JENDL/FPD-2011, the zero half-life is given to this nuclide
  
  int sz=size_array[num]; // # of FP nuclides

  int z;
  int a;
  int l;
  matrix[num].set_zero();

  vector<real> unc_idp_sqr(sz);
  for(int i=0;i<sz;i++){
    unc_idp_sqr[i]=pow(unc_idp[num][i],2.);
  };

  vector<int> z_array;
  vector<int> a_array;
  vector<int> l_array;
  int a_max=0;
  for(int i=0;i<sz;i++){
    eidt.GetParameterNew(id[num][i],z,a,l);
    z_array.push_back(z);
    a_array.push_back(a);
    l_array.push_back(l);
    if(a>a_max)a_max=a;
  };

  vector<bool> exist(a_max,false);
  for(int i=0;i<sz;i++){
    exist[a_array[i]]=true;
  };
  
  cout<<"# Covariance matrix of fission yield is constructed for "<<fisnucname_array[num]<<"\n";

  for(int i=1;i<=a_max;i++){
    if(exist[i]){

    //+++Searching Mass Yield+++    
    vector<real> mass_unc;     // absolute variance
    vector<real> mass_unc_idp; // absolute variance
    vector<int> mass_nuc;
    vector<int> mass_z;
    int mass_nuc_num=0;
    vector<int> mass_chain_id;//1,2,3,...
    mass_chain_id.resize(sz,0);

    for(int j=0;j<sz;j++){

      z=z_array[j];
      a=a_array[j];
      l=l_array[j];

      if(bm.GetNuclideIndex(z,a,l)==-1){
	bm.AddNuclide(z,a,l);
      };

      if(a==i){
	int channel=bm.GetNuclide(z,a,l).GetChannel();
	real hl=bm.GetNuclide(z,a,l).GetHalflife();
        bool end_nuclide=false; // To detect the end nuclide in the mass chain
	//if(hl==0.){
	if(hl==0.||hl>1e20){ // modified in 2022/9/13.	  
          end_nuclide=true;
	}else if(channel==1){
	  int decay_type=bm.GetNuclide(z,a,l).GetDecayType(0);
	  if(decay_type==3)end_nuclide=true; // (Alpha-decay)
	};
	if(end_nuclide){
	  /*
	  cout<<"# The final nuclide in the mass chain "<<a<<" is detected.\n";
	  cout<<"#    (Z,A,L) of this nuclide : "<<z<<","<<a<<","<<l<<"\n";
	  cout<<"#    Cumulative yield uncertainty : "<<unc_cum[num][j]<<"\n";
	  */
	  /*
	  mass_unc.push_back(unc_cum[num][j]);
	  mass_unc_idp.push_back(unc_idp[num][j]); // absolute standard deviation
	  */
	  mass_unc.push_back(pow(unc_cum[num][j],2));
	  //mass_unc_idp.push_back(pow(unc_idp[num][j],2)); // absolute standard deviation
	  mass_unc_idp.push_back(unc_idp_sqr[j]); // absolute standard deviation
	  mass_nuc.push_back(id[num][j]);
	  mass_z.push_back(z);
	  mass_nuc_num++;
	  mass_chain_id[j]=mass_nuc_num;
	};
      };

    };

    vector<real> sum_unc_square;
    sum_unc_square.resize(mass_nuc_num);

    if(mass_nuc_num>=1){

      //cout<<"# Mass number : "<<i<<" / Number of nuclides : "<<mass_nuc_num<<"\n";
      /*
      cout<<"#\n# Mass number is : "<<i<<"\n#\n";
      for(int jj=0;jj<mass_nuc_num;jj++){
	cout<<"    "<<jj<<" "<<mass_nuc[jj]<<" "<<mass_z[jj]<<" "<<mass_unc[jj]<<" "<<mass_unc_idp[jj]<<"\n";
      };
      */
      real tmp2;
      for(int jj=0;jj<mass_nuc_num-1;jj++){
	tmp2=0.;
	for(int j=0;j<sz;j++){
          z=z_array[j];
          a=a_array[j];
	  if(a==i&&mass_z[jj]>z&&mass_chain_id[j]==0){ // Beta-minus decay series
	    mass_chain_id[j]=jj+1;
	    //tmp2+=pow(unc_idp[num][j],2.);
	    tmp2+=unc_idp_sqr[j];
	  };
	}; 
	//sum_unc_square[jj]=tmp2+mass_unc_idp[jj]; // Kawamoto-Original
	//sum_unc_square[jj]=tmp2+pow(mass_unc_idp[jj],2); // Chiba-Modified
	sum_unc_square[jj]=tmp2+mass_unc_idp[jj]; // Chiba-Modified
      };

      tmp2=0.; // Beta-plus decay series
      for(int j=0;j<sz;j++){
        a=a_array[j];
	if(a==i&&mass_chain_id[j]==0){
	  mass_chain_id[j]=mass_nuc_num;
	  //tmp2+=pow(unc_idp[num][j],2.);
	  tmp2+=unc_idp_sqr[j];
	};
      };
      //sum_unc_square[mass_nuc_num-1]=tmp2+mass_unc_idp[mass_nuc_num-1]; // Kawamoto-Original
      //sum_unc_square[mass_nuc_num-1]=tmp2+pow(mass_unc_idp[mass_nuc_num-1],2); // Chiba-modified
      sum_unc_square[mass_nuc_num-1]=tmp2+mass_unc_idp[mass_nuc_num-1]; // Chiba-modified

      //++++++++++++++++++++++++++
  
      for(int ii=0;ii<mass_nuc_num;ii++){
    
	for(int k=0;k<sz;k++){
	  if(mass_chain_id[k]==ii+1){
	    for(int kk=0;kk<sz;kk++){
	      if(kk==k){
  	        real tmp1=yield_idp[num][k]*yield_idp[num][kk];
	        if(tmp1!=0.){
	          //tmp=pow(unc_idp[num][k],2.)*(1.-(pow(unc_idp[num][k],2.))/(pow(mass_unc[ii],2.)+sum_unc_square[ii]))/tmp1;
                  //real tmp2=pow(unc_idp[num][k],2);
                  real tmp2=unc_idp_sqr[k];
	          real tmp=tmp2*(1.-tmp2/(mass_unc[ii]+sum_unc_square[ii]))/tmp1;
  	          matrix[num].put_data(k,kk,tmp);
	        };
  	      }else if(mass_chain_id[kk]==ii+1){
  	        real tmp1=yield_idp[num][k]*yield_idp[num][kk];
	        if(tmp1!=0.){
	          //tmp=-(pow(unc_idp[num][k],2.)*pow(unc_idp[num][kk],2.))/(pow(mass_unc[ii],2.)+sum_unc_square[ii])/tmp1;
	          //real tmp=-(pow(unc_idp[num][k],2.)*pow(unc_idp[num][kk],2.))/(mass_unc[ii]+sum_unc_square[ii])/tmp1;
	          real tmp=-unc_idp_sqr[k]*unc_idp_sqr[kk]/(mass_unc[ii]+sum_unc_square[ii])/tmp1;
  	          matrix[num].put_data(k,kk,tmp);
	        };
	      };
	    };
	  };
	};

      };
      // ++++++++++++++++++++++++++++++++

    };
  };

  };
};

void IndependentYieldCovariance::CalCovarianceGLSUpdating(BCGManager &bm,int num, bool *condition)
{
  int z_lcp=14;    // if Z<=[z_lcp], this nuclide is regarded light nuclide
  real std_ch=0.1; // Relative standard deviation of chain yield
  real std_zy=0.0001;
  real std_ay=0.0001;
  real std_fy=0.0001;
  real std_fyupper=0.0001;

  int fpnum=GetSize(num); // The number of fission products 
  cout<<"# The number of FP nuclides : "<<fpnum<<"\n";

  GroupData2D dmat(fpnum,fpnum);
  dmat.set_zero();

  vector<real> halflife(fpnum);

  for(int i=0;i<fpnum;i++){
    int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
    int z,a,l;
    midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
    halflife[i]=bm.GetNuclide(z,a,l).GetHalflife();

    int channel=bm.GetNuclide(z,a,l).GetChannel();
    bool end_nuclide=false; // To detect the end nuclide in the mass chain
    if(halflife[i]==0.){
      end_nuclide=true;}
    else if(channel==1){
      int decay_type=bm.GetNuclide(z,a,l).GetDecayType(0);
      if(decay_type==3)end_nuclide=true; // (Alpha-decay)
    };
    if(end_nuclide){
      //if(halflife[i]==0.){
      dmat.put_data(i,i,1.);
    }else{
      int chnum=bm.GetNuclide(z,a,l).GetChannel();
      for(int j=0;j<chnum;j++){
        int z2=bm.GetNuclide(z,a,l).GetAtomicNumberNext(j);
        int a2=bm.GetNuclide(z,a,l).GetMassNumberNext(j);
        int l2=bm.GetNuclide(z,a,l).GetExLevelNext(j);
	real br=bm.GetNuclide(z,a,l).GetBr(j);
	int fpid2=midt.GetMATID(z2,a2,l2);
	bool datainp=false;
	for(int k=0;k<fpnum;k++){
	  if(GetID(num,k)==fpid2){
	    dmat.put_data(k,i,br);
	    datainp=true;
	  };
	};
      };

    };
  };

  GroupData2D emat=dmat;
  int itermax=100;
  for(int iter=0;iter<itermax;iter++){
    GroupData2D tmp=emat;
    emat=dmat*emat;
    tmp=tmp-emat;
    bool diff=false;
    for(int j=0;j<fpnum*fpnum;j++){
      if(fabs(tmp.get_dat(j))>1e-20)diff=true;
    };
    if(!diff){
      cout<<"# Matrix is unchangend at iteration "<<iter<<"\n";
      break;
    };
  };

  int cy_num =0.;
  for(int i=0;i<fpnum;i++){
    bool non_zero=false;
    for(int j=0;j<fpnum;j++){
      if(emat.get_dat(i,j)>0.)non_zero=true;
    };
    int fpid=GetID(num,i); 
    if(non_zero){
      cy_num=cy_num+1;
      //cout<<fpid<<" "<<halflife[i]<<"\n";
    };
  };
  cout<<"#\n";
  //////////////////////////////////  ////  ////  /////  ///////////////////////////////////////////////

  real mass_sum=0.;
  real p_sum=0.;
  real fy_sum=0.;
  for(int i=0;i<fpnum;i++){
    int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
    int z,a,l;
    midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
    real fy =GetYieldIDP(num,i); // Independent fission yield
    real dfy=GetDeltaYieldIDP(num,i); // Absolute uncertainty of independent FY
    //cout<<i<<" "<<icov.GetID(tag,i)<<" "<<z<<" "<<a<<" "<<l<<" "<<fy<<" "<<dfy<<"\n";
    mass_sum+=a*fy;
    p_sum+=z*fy;
    fy_sum+=fy;
  };
  cout<<"# Sum of mass   : "<<mass_sum<<"\n";
  cout<<"# Sum of proton : "<<p_sum<<"\n";
  cout<<"# Sum of yield  : "<<fy_sum<<"\n";
  cout<<"# Number of chain yield : "<<cy_num<<"\n";

  real vp = 2.4169; //for  U235

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  GroupData1D y(fpnum);          // Fission yield mean vector
  GroupData2D cov(fpnum,fpnum);    // Fission yield ABSOLUTE covariance matrix
 
  for(int i=0;i<fpnum;i++){      // Fission yield vector
    y.put_data(i,GetYieldIDP(num,i));
  };

  cov.set_zero();
  for(int i=0;i<fpnum;i++){      // Covariance matrix
    cov.put_data(i,i,pow(GetDeltaYieldIDP(num,i),2));
  };
  real nd_sum_v=cov.get_sum();

  GroupData2D tmp11,tmp10;
  GroupData1D Cy,Cy_const;
  GroupData2D gx; // sensitivity matrix
  GroupData2D Vem; 

  ////////////////////////// 1      Chain yields             //////////////////////
  ////////////////////////// 2 Sum of ZY is equal to Zcn-Zlcp //////////////////////
  //////////////////////////// 3 Mass number AY          ////////////////////////////
  //////////////////////////// 4 Sum of fy is equal to 2 //////////////////////////////
  //////////////////////////// 5 Sum of fy is equal to 1. //////////////////////////////

  int num_output_list[]={fpnum, 1, 1, 1, 1};
  //bool condition[]={true,true,true,true,true};
  for(int ii=0;ii<5;ii++){

    if(condition[ii]){

    int num_output=num_output_list[ii];
    Cy.put_imax(num_output);
    Vem.put_yx(num_output,num_output); // Covariance of constrained output data
    gx.put_yx(num_output,fpnum);       // Sensitivity matrix
    gx.set_zero();

    if(ii==0){
  ////////////////////////// 1      Chain yields             //////////////////////
      Cy=emat*y;       // Output data calculated from input data
      Cy_const=Cy;  
      gx=emat;         // Sensitivity matrix
      for(int i=0;i<fpnum;i++){     
        real tmp=pow(GetDeltaYieldCUM(num,i),2);
        //real tmp=pow(Cy_const.get_dat(i)*std_ch,2);
        Vem.put_data(i,i,tmp);
      };
    }else if(ii==1){
  ////////////////////////// 2 Sum of ZY is equal to Zcn-Zlcp //////////////////////
      real Zlcp = 0.;
      real Zall = 0.;
      for(int i=0;i<fpnum;i++){
        int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
        int z,a,l;
        midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
        real fy =GetYieldIDP(num,i); // Independent fission yield
        Zall+=z*fy;
        if(z<=z_lcp){
          Zlcp+=z*fy;
        }else{
          gx.put_data(0,i,z); //Sensitivity vector ZY
        };
      };
      cout<<"#      Z_all : "<<Zall<<"\n";
      cout<<"#      Z_lcp : "<<Zlcp<<"\n";
      Cy.put_data(0,Zall-Zlcp);
      Cy_const=Cy;

      for(int i=0;i<num_output;i++){     
        real tmp=pow(Cy_const.get_dat(i)*std_zy,2);
        Vem.put_data(i,i,tmp);
      };
    }else if(ii==2){
  //////////////////////////// 3 Mass number AY          ////////////////////////////

      real Aall  = 0.;
      real Alcp = 0.;
      for(int i=0;i<fpnum;i++){
        int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
        int z,a,l;
        midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
        real fy =GetYieldIDP(num,i); // Independent fission yield
        Aall+=a*fy;
        if(z<=z_lcp){
          Alcp+=a*fy;
        }else{
          gx.put_data(0,i,a); //Sensitivity vector ZY
        };
      };

      Cy.put_data(0,Aall-Alcp-vp);
      Cy_const=Cy;
  
      for(int i=0;i<num_output;i++){     
        real tmp=pow(Cy_const.get_dat(i)*std_ay,2);
        Vem.put_data(i,i,tmp);
      };
    }else if(ii==3){
  //////////////////////////// 4 Sum of fy is equal to 2 //////////////////////////////
      real FYall = 0.;
      real FYlcp = 0.;
      for(int i=0;i<fpnum;i++){
        int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
        int z,a,l;
        midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
        real fy =GetYieldIDP(num,i); // Independent fission yield
        FYall+=fy;
        if(z<=z_lcp){
          FYlcp+=fy;
        }else{
          gx.put_data(0,i,1); //Sensitivity vector ZY
        };
      };
      Cy.put_data(0,FYall-FYlcp);
      Cy_const=Cy;

      for(int i=0;i<num_output;i++){     
        real tmp=pow(Cy_const.get_dat(i)*std_fy,2);
        Vem.put_data(i,i,tmp);
      };
    }else if(ii==4){
  //////////////////////////// 5 Sum of fy is equal to 1. //////////////////////////////
      real FYall5   = 0.;
      real FYlower  = 0.;
      real FYupper  = 0.;
      real Ai = mass_sum/2 ;
      for(int i=0;i<fpnum;i++){
        int fpid=GetID(num,i); // FP index of the i-th FP is retrieved.
        int z,a,l;
        midt.GetParameter(fpid,z,a,l); // A pair (Z,A,L) is obtained from nuclide index [fpid].
        real fy =GetYieldIDP(num,i); // Independent fission yield
        FYall5+=fy;
        if(a<=Ai){
          FYlower+=fy;
        }else{
          gx.put_data(0,i,1); //Sensitivity vector ZY
        };
      };
      FYupper = FYall5 - FYlower;
      Cy.put_data(0,FYupper);
      Cy_const=Cy;

      for(int i=0;i<num_output;i++){     
        real tmp=pow(Cy_const.get_dat(i)*std_fyupper,2);
        Vem.put_data(i,i,tmp);
      };
    };

    tmp10=cov*gx.T();
    tmp11=tmp10* (gx*tmp10 + Vem).inverse();
    y=y+tmp11*(Cy_const-Cy);  // Modified FY
    cov = cov-tmp11*gx*cov;   // Modified covariance matrix

    cout<<"########  Condition "<<ii+1<<" ########" <<"\n";
    cout<<"# Sum of modified yield  : "    << y.get_sum() <<"\n";
    cout<<"# Variance of modified yield : "<< cov.get_sum() << "\n";

    };
  };

  int sz=size_array[num];
  for(int i=0;i<sz;i++){
    real idp1=yield_idp[num][i];
    for(int j=i;j<sz;j++){
      real idp2=yield_idp[num][j];
      if(idp1>0.&&idp2>0.){
        real tmp=cov.get_dat(i,j)/(idp1*idp2);
        matrix[num].put_data(i,j,tmp);
        if(i!=j)matrix[num].put_data(j,i,tmp);
      };
    };
  };
  
};

void IndependentYieldCovariance::CalCovarianceWithoutCorrelation(int num){
  real tmp;
  real tmp1;
  matrix[num].set_zero();

  for(int k=0;k<size_array[num];k++){
    tmp1=yield_idp[num][k]*yield_idp[num][k];
    if(tmp1==0.){
      tmp=0.;
    }else{
      tmp=unc_idp[num][k]*unc_idp[num][k]/tmp1;
    };
    matrix[num].put_data(k,k,tmp);
  };
};

void IndependentYieldCovariance::PlotCovarianceMatrix(int num){
  cout<<"#size:"<<GetSize(num)<<"\n";
  for(int i=0;i<size_array[num];i++){
    for(int j=0;j<size_array[num];j++){
      if(matrix[num].get_dat(i,j)!=0.){
	cout<<i+1<<" "<<j+1<<"\n";
      };
    };
  };
};

int IndependentYieldCovariance::GetNumberFromFissionNuclideName(string fisnucname){
  int tmp=-1;
  for(int i=0;i<fisnucnum;i++){
    if(fisnucname_array[i]==fisnucname){
      tmp=i;
    };
  };
  return tmp;
};

void IndependentYieldCovariance::DeleteCovarianceExceptAssignedMassChainFissionNuclide(int mass,int mat){
  for(int i=0;i<fisnucnum;i++){
    if(mat==GetFissionNuclideID(i)){
      for(int j=0;j<size_array[i];j++){
	int z;
	int a;
	int l;
	midt.GetParameter(id[i][j],z,a,l);
	if(mass!=a){
	  for(int k=0;k<size_array[i];k++){
	    matrix[i].put_data(j,k,0.);
	  };
	  for(int k=0;k<size_array[i];k++){
	    matrix[i].put_data(k,j,0.);
	  };
	};
      };
    }else{
      matrix[i].set_zero();
    };    
  };
};

void IndependentYieldCovariance::ShowVarianceData(string fisnucname)
{
  int ind=GetNumberFromFissionNuclideName(fisnucname);
  if(ind!=-1){
    int sz=size_array[ind];
    for(int i=0;i<sz;i++){
      cout<<id[ind][i]<<" "<<yield_idp[ind][i]<<" "<<yield_cum[ind][i]<<"\n";
    };
  };
};

GroupData2D IndependentYieldCovariance::GetPartialMatrix(string fisnucname, int num, int *idlist)
{
  int ind=GetNumberFromFissionNuclideName(fisnucname);
  if(ind==-1){
    cout<<"# Warning in IndependentYieldCovariance::GetPatrialMatrix.\n";
    cout<<"# No fissile data named as"<<fisnucname<<"\n";
    exit(0);
  };

  int sz=size_array[ind];

  vector<int> pos(num,-1);
  for(int j=0;j<num;j++){
    for(int i=0;i<sz;i++){
      if(idlist[j]==id[ind][i])pos[j]=i;
    };
    if(pos[j]==-1){
      cout<<"# Error in IndependentYieldCovariance::GetPatrialMatrix.\n";
      cout<<"# No data whose ID is "<<idlist[j]<<"\n";
      exit(0);
    };
  };

  GroupData2D pmat(num,num);
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      int ii=pos[i];
      int jj=pos[j];
      real dat=matrix[ind].get_dat(ii,jj);
      pmat.put_data(i,j,dat);
      if(i!=j)pmat.put_data(j,i,dat);
    };
  };

  return pmat;
};

real IndependentYieldCovariance::GetYieldIDP(string fisnucname, int i)
{
  int tmp=GetNumberFromFissionNuclideName(fisnucname);
  if(tmp!=-1)return yield_idp[tmp][i];
};

real IndependentYieldCovariance::GetDeltaYieldIDP(string fisnucname, int i)
{
  int tmp=GetNumberFromFissionNuclideName(fisnucname);
  if(tmp!=-1)return unc_idp[tmp][i];
};

real IndependentYieldCovariance::GetDeltaYieldCUM(string fisnucname, int i)
{
  int tmp=GetNumberFromFissionNuclideName(fisnucname);
  if(tmp!=-1)return unc_cum[tmp][i];
};

real IndependentYieldCovariance::CalVarianceForMixedYield(string fisnucname, int mat1, int mat2)
{
  GroupData2D mat=GetCovarianceMatrix(fisnucname);
  real relcov1=mat.get_dat(mat1,mat1);
  real relcov2=mat.get_dat(mat2,mat2);
  real relcorr=mat.get_dat(mat1,mat2);
  int tmp=GetNumberFromFissionNuclideName(fisnucname);
  real abscov1=relcov1*yield_idp[tmp][mat1]*yield_idp[tmp][mat1];
  real abscov2=relcov2*yield_idp[tmp][mat2]*yield_idp[tmp][mat2];
  real abscorr=relcorr*yield_idp[tmp][mat1]*yield_idp[tmp][mat2];
  //cout<<relcov1<<" "<<relcov2<<" "<<relcorr<<" "<<"\n";
  //cout<<abscov1<<" "<<abscov2<<" "<<abscorr<<" "<<"\n";
  //cout<<yield_idp[tmp][mat1]<<" "<<yield_idp[tmp][mat2]<<"\n";
  //cout<<abscov1<<"*+"<<abscov2<<"+2x"<<abscorr<<"\n";
  cout<<abscov1+abscov2+2*abscorr<<"\n";
  //if(tmp!=-1)return abscov1+abscov2+2*abscorr;
  //if(tmp!=-1)return yield_idp[tmp][i];
};

// Prepared by Okumura

void IndependentYieldCovariance::ReadReducedChainData(const char *wfilen)
{
  int num;
  ifstream fin;
  //cout<<"# Reading data from file :" <<wfilen<<"\n";
  fin.open(wfilen,ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file : "<<wfilen<<"\n";
    exit(0);
  };

  fin>>fisnucnum;
  size_array.resize(fisnucnum);
  fin>>num;
  //vector< vector< vector<real> > > cov_sh(num);
  //cov_sh.resize(fisnucnum);
  fisnucname_array.resize(fisnucnum);
  matrix.resize(fisnucnum);
  id.resize(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    //cov_sh[i].resize(num);
    id[i].resize(num);
    size_array[i]=num;
    for(int j=0;j<num;j++){
      //cov_sh[i][j].resize(num);
    };
  };

  //cout<<"# ShortenBurnupChain..."<<num<<"nuclides  "<<fisnucnum<<" Fissile nuclide "<<"\n";
  cout<<"# Covariance data for FP yield of reduced chain are read from file.\n#\n";
  cout<<"#   The number of fissile nuclides : "<<fisnucnum<<"\n";
  cout<<"#   The number of FP nuclides of the reduced chain : "<<num<<"\n";

  for(int i=0;i<fisnucnum;i++){
    matrix[i].put_yx(num,num);
  };

  for(int j=0;j<num;j++){
    int tmp;
    fin>>tmp;
    for(int i=0;i<fisnucnum;i++){
      id[i][j]=tmp;
    };
  };
    
  for(int i=0;i<fisnucnum;i++){
    int tmp1;
    real tmp2;
    fin>>tmp1;
    fisnucname_array[i]=midt.Name(tmp1);
    //fin>>fisnucname_array[i];
    for(int j=0;j<num;j++){
      for(int k=0;k<num;k++){
	fin>>tmp2;
	matrix[i].put_data(j,k,tmp2);
	//fin>>cov_sh[i][j][k];
      };
    };
    for(int j=0;j<num;j++){
      for(int k=0;k<num;k++){
	//matrix[i].put_data(j,k,cov_sh[i][j][k]);
      };
    };
  }; // i-loop

  //test//
  /*
  cout<<"check1 "<<fisnucnum<<"\n";
  for(int i=0;i<fisnucnum;i++){
    cout<<"check2 "<<size_array[i]<<"\n";
  };
  for(int i=0;i<fisnucnum;i++){
    cout<<"check3 "<<fisnucname_array[i]<<"\n";
  };
  for(int i=0;i<fisnucnum;i++){
    for(int j=0;j<size_array[i];j++){
      cout<<i<<" "<<j<<" check4 "<<id[i][j]<<"\n";
    };
  };
  for(int i=0;i<fisnucnum;i++){
    matrix[i].show_self();
  };
  */  
  fin.close();
};

void IndependentYieldCovariance::PutCovarianceMatrix(string fisnucname, GroupData2D &matinp)
{
  int id=GetNumberFromFissionNuclideName(fisnucname);
  if(id==-1){
    cout<<"# Error in IndependentYieldCovariance::PutCovarianceMatrix.\n";
    cout<<"# No tag ["<<fisnucname<<"] is found.\n";
    exit(0);
  };
  matrix[id]=matinp;
};

void IndependentYieldCovariance::ConversionFromAbsoluteToRelative(string fisnucname)
{
  int id=GetNumberFromFissionNuclideName(fisnucname);
  if(id==-1){
    cout<<"# Error in IndependentYieldCovariance::ConversionFromAbsoluteToRelative.\n";
    cout<<"# No tag ["<<fisnucname<<"] is found.\n";
    exit(0);
  };

  int sz=size_array[id];
  for(int i=0;i<sz;i++){
    real idp1=yield_idp[id][i];
    for(int j=i;j<sz;j++){
      real idp2=yield_idp[id][j];
      if(idp1>0.&&idp2>0.){
        real tmp=matrix[id].get_dat(i,j)/(idp1*idp2);
        matrix[id].put_data(i,j,tmp);
        if(i!=j)matrix[id].put_data(j,i,tmp);
      };
    };
  };
};


//++++++(HalfLifeCovariance)++++++


void HalfLifeCovariance::SetDefaultCovariance()
{
  BCGManager bm;
  //bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.endf71");
  // If the above line is NOT commended, uncertainty data
  // of non-FP nuclides are taken into account.
  bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.jendl2011");
  MakeCovarianceMatrix(bm);
};

void HalfLifeCovariance::MakeCovarianceMatrix(BCGManager &bm)
{
  size=bm.GetSize();
  halflife.resize(size);
  unc.resize(size);
  id.resize(size);
  matrix.put_yx(size,size);
  matrix.set_zero();

  for(int i=0;i<size;i++){
    id[i]=bm.GetNuclide(i).GetID();    
    real tmph=bm.GetNuclide(i).GetHalflife();
    real tmpu=bm.GetNuclide(i).GetDeltaHalflife();
    halflife[i]=tmph;
    //If the error is zero even though the half life is not zero, we use assumed value (dafault; 100%).
    if(tmph!=0.&&tmpu==0.){
      unc[i]=tmph*theor_err;//assumed
    }else{
      unc[i]=tmpu;
    };
  };

  real tmp1;
  real tmp2;
  for(int i=0;i<size;i++){
    tmp1=halflife[i]*halflife[i];
    if(tmp1==0.){
      tmp2=0.;
    }else{
      tmp2=unc[i]*unc[i]/tmp1;
    };
    matrix.put_data(i,i,tmp2);
  };
  
};

void HalfLifeCovariance::DeleteCovarianceExceptAssignedNuclide(int mat){
  for(int i=0;i<size;i++){
    if(mat!=id[i]){
      matrix.put_data(i,i,0.);
    };
  };
};

void HalfLifeCovariance::ShowSelf(){
  for(int i=0;i<size;i++){
    real var=matrix.get_dat(i,i);
    cout<<"# "<<id[i]<<" "<<halflife[i]<<" "<<var<<" "<<sqrt(var)<<"\n";
  };
};

//++++++(BranchingRatioCovariance)++++++

void BranchingRatioCovariance::SetDefaultCovariance()
{
  BCGManager bm;
  //bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.endf71");
  bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.jendl2011");
  MakeCovarianceMatrix(bm);
};

void BranchingRatioCovariance::MakeCovarianceMatrix(BCGManager &bm){
  
  nucnum=bm.GetSize();
  
  ratio.resize(nucnum);
  unc.resize(nucnum);
  id.resize(nucnum);
  channel_num.resize(nucnum);
  
  int counter=0;
  for(int i=0;i<nucnum;i++){
    id[i]=bm.GetNuclide(i).GetID();    
    int channel=bm.GetNuclide(i).GetChannel();
    channel_num[i]=channel;
    ratio[i].resize(channel);
    unc[i].resize(channel);
    for(int j=0;j<channel;j++){
      counter++;
      real br=bm.GetNuclide(i).GetBr(j);
      real dbr=bm.GetNuclide(i).GetDeltaBr(j);
      ratio[i][j]=br;
      //If the error is zero even though the half life is not zero, we use assumed value (dafault; 100%).
      if(dbr==0.&&channel!=1&&br!=0.){
	unc[i][j]=br*theor_err;//assumed
      }else{
	unc[i][j]=dbr;
      };
    };
  };
  
  size=counter;
  matrix.put_yx(size,size);
  matrix.set_zero();
  list_id.resize(size);
  list_channel.resize(size);

  /*
  // +++ relative covariance matrix construction
  counter=0;
  for(int i=0;i<nucnum;i++){
    int channel=bm.GetNuclide(i).GetChannel();
    for(int j=0;j<channel;j++){
      list_id[counter+j]=bm.GetNuclide(i).GetID();
      list_channel[counter+j]=j;
      real tmp;
      real tmp1=ratio[i][j]*ratio[i][j];
      if(tmp1==0.||channel==1){
	tmp=0.;
      }else{
	tmp=unc[i][j]*unc[i][j]/tmp1;
      };
      matrix.put_data(counter+j,counter+j,tmp);
    };
    counter+=channel;
  };
  */

  // +++ ABSOLUTE covariance matrix construction
  counter=0;
  for(int i=0;i<nucnum;i++){
    int channel=bm.GetNuclide(i).GetChannel();
    for(int j=0;j<channel;j++){
      list_id[counter+j]=bm.GetNuclide(i).GetID();
      list_channel[counter+j]=j;
      real tmp=0.;
      if(channel!=1&&ratio[i][j]!=0.)tmp=unc[i][j]*unc[i][j];
      matrix.put_data(counter+j,counter+j,tmp);
    };
    counter+=channel;
  };

  // +++ zero-sum rule correction
  GroupData2D matrix_new;
  matrix_new.put_yx(size,size);
  matrix_new.set_zero();

  counter=0;
  for(int i=0;i<nucnum;i++){
    int channel=bm.GetNuclide(i).GetChannel();
    for(int j1=0;j1<channel;j1++){
      for(int j2=0;j2<channel;j2++){
        real org=matrix.get_dat(counter+j1,counter+j2);
        for(int k1=0;k1<channel;k1++){
	  org-=ratio[i][j2]*matrix.get_dat(counter+j1,counter+k1);
	  org-=ratio[i][j1]*matrix.get_dat(counter+k1,counter+j2);
	  for(int k2=0;k2<channel;k2++){
	    org+=ratio[i][j1]*ratio[i][j2]*matrix.get_dat(counter+k1,counter+k2);
	  };
	};
	matrix_new.put_data(counter+j1,counter+j2,org);
      };
    };
    counter+=channel;
  };

  matrix.copy(matrix_new);

  // +++ absolute -> RELATIVE covariance
  counter=0;
  for(int i=0;i<nucnum;i++){
    int channel=bm.GetNuclide(i).GetChannel();
    for(int j1=0;j1<channel;j1++){
      real v1=ratio[i][j1];
      for(int j2=0;j2<channel;j2++){
	real v2=ratio[i][j2];
        real org=matrix.get_dat(counter+j1,counter+j2);
	if(org!=0.)matrix.put_data(counter+j1,counter+j2,org/v1/v2);
      };
    };
    counter+=channel;
  };

};

void BranchingRatioCovariance::DeleteCovarianceExceptAssignedNuclide(int mat){
  int counter=0;
  for(int i=0;i<nucnum;i++){
    for(int j=0;channel_num[i];j++){      
      if(mat!=id[i]){
	matrix.put_data(counter,counter,0.);
      };
      counter++;
    };
  };
};

GroupData2D BranchingRatioCovariance::GetCovarianceMatrix(int mat)
{
  int sz=id.size();
  int chnum=0;
  for(int i=0;i<sz;i++){
    if(id[i]==mat)chnum=channel_num[i];
  };
  
  if(chnum==0){
    cout<<"# Error in BranchingRatioCovariance::GetCovarianceMatrix.\n";
    cout<<"# No decay channel for MATID : "<<mat<<"\n";
    exit(0);
  };

  GroupData2D ret(chnum,chnum);
  for(int i=0;i<size;i++){
    if(list_id[i]==mat&&list_channel[i]==0){
      for(int j=0;j<chnum;j++){
	for(int k=0;k<chnum;k++){
	  ret.put_data(j,k,matrix.get_dat(i+j,i+k));
	};
      };
      return ret;
    };
  };

  cout<<"# Error in BranchingRatioCovariance::GetCovarianceMatrix.\n";
  cout<<"# Please check the source code.\n";
  exit(0);
};
void BranchingRatioCovariance::ReadReducedChainData(const char *wfilen)
{
  ifstream fin;
  fin.open(wfilen,ios::in);
  if(fin.fail()){
    cout<<"# Fail to open the file : "<<wfilen<<"\n";
    exit(0);
  };
  fin>>size;
  fin>>nucnum;
  id.resize(nucnum);
  ratio.resize(nucnum);
  channel_num.resize(nucnum);
  for(int i=0;i<nucnum;i++){
   int tmp;
    fin>>tmp;
    id[i]=tmp;
  };
  for(int i=0;i<nucnum;i++){
    int tmp_channel;
    fin>>tmp_channel;
    channel_num[i]=tmp_channel;
    ratio[i].resize(tmp_channel);
  };
  matrix.put_yx(size,size);
  for(int j=0;j<size;j++){
    for(int k=0;k<size;k++){
      real tmp;
      fin>>tmp;
      matrix.put_data(j,k,tmp);
    };
  };
  fin.close();
};

void BranchingRatioCovariance::ResetCovarianceMatrix()
{
  matrix.set_zero();
};

void BranchingRatioCovariance::DirectInputCovarianceData(int target_id, real var1, real var2, real cov)
{
  int index=0;
  for(int i=0;i<nucnum;i++){
    if(id[i]==target_id){
      if(channel_num[i]!=2){
	cout<<"# Error in BranchingRatioCovariance::DirectInputCovarianceData.\n";
	cout<<"# The number of branches should be 2.\n";
	exit(0);
      };
      matrix.put_data(index,index,var1);
      matrix.put_data(index+1,index+1,var2);
      matrix.put_data(index,index+1,cov);
      matrix.put_data(index+1,index,cov);                  
    };
    index+=channel_num[i];
  };
};
  
void BranchingRatioCovariance::ShowSelf()    
{
  cout<<"#\n";
  cout<<"# The number of nuclides : "<<nucnum<<"\n";
  cout<<"# The size of the matrix : "<<matrix.get_x()<<"\n";
  cout<<"#\n";

  int index=0;
  for(int i=0;i<nucnum;i++){
    cout<<id[i]<<" : "<<channel_num[i]<<"\n";
    if(channel_num[i]>0){
      for(int j=0;j<channel_num[i];j++){
	cout<<"           "<<matrix.get_dat(index,index)<<" "<<sqrt(matrix.get_dat(index,index))<<" "<<index<<"\n";
	index++;
      };
    };
  };
};

//++++++(DecayEnergyCovariance)++++++

void DecayEnergyCovariance::SetDefaultCovariance()
{
  BCGManager bm;
  //bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.endf71");
  // If the above line is NOT commended, uncertainty data
  // of non-FP nuclides are taken into account.
  bm.ReadDecayDataFromFile("../../CBGLIB_BURN/DDfile/dd.jendl2011");
  MakeCovarianceMatrix(bm);
};


void DecayEnergyCovariance::MakeCovarianceMatrix(BCGManager &bm){
  size=bm.GetSize();
  decay_energy.resize(size);
  unc.resize(size);
  id.resize(size);
  matrix.resize(3);
  matrix[0].put_yx(size,size);
  matrix[0].set_zero();
  matrix[1].put_yx(size,size);
  matrix[1].set_zero();
  matrix[2].put_yx(size,size);
  matrix[2].set_zero();

  for(int i=0;i<size;i++){
    id[i]=bm.GetNuclide(i).GetID();
    decay_energy[i].resize(3);
    unc[i].resize(3);
    for(int j=0;j<3;j++){
      real tmpe=bm.GetNuclide(i).GetDecayEnergy(j);
      real tmpu=bm.GetNuclide(i).GetDeltaDecayEnergy(j);
      decay_energy[i][j]=tmpe;
      //If the error is zero even though the energy is not zero, we use assumed value (dafault; 100%).
      if(tmpe!=0.&&tmpu==0.){
	unc[i][j]=tmpe*theor_err;//assumed
      }else{
	unc[i][j]=tmpu;
      };
      //+++++++++++++++
    };
  };
  
  real tmp1;
  real tmp2;
  for(int j=0;j<3;j++){
    for(int i=0;i<size;i++){
      tmp1=decay_energy[i][j]*decay_energy[i][j];
      if(tmp1==0.){
	tmp2=0.;
      }else{
	tmp2=unc[i][j]*unc[i][j]/tmp1;
      };
      matrix[j].put_data(i,i,tmp2);
    };
  };

};

void DecayEnergyCovariance::DeleteCovarianceExceptAssignedNuclideDecayType(int mat,int type){
  for(int j=0;j<3;j++){
    if(j==type){
      for(int i=0;i<size;i++){
	if(mat!=id[i]){
	  matrix[j].put_data(i,i,0.);
	};
      };
    }else{
      matrix[j].set_zero();
    };
  };
};

void DecayEnergyCovariance::DeleteCovarianceExceptAssignedNuclideDecayType(int mat,string type){
  if(type=="beta"||type=="Beta"||type=="BETA"){
    DeleteCovarianceExceptAssignedNuclideDecayType(mat,0);
  }else if(type=="gamma"||type=="Gamma"||type=="GAMMA"){
    DeleteCovarianceExceptAssignedNuclideDecayType(mat,1);
  }else if(type=="alpha"||type=="Alpha"||type=="ALPHA"){
    DeleteCovarianceExceptAssignedNuclideDecayType(mat,2);
  };
};

void FissionYieldCovarianceCorr::ReadReducedChainData(const char *wfilen){
  int num;  // Number of FP
  ifstream fin;
  fin.open(wfilen,ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file : "<<wfilen<<"\n";
    exit(0);
  };
  fin>>fisnucnum;
  fin>>num;
  fisnucname_array.resize(fisnucnum);
  size=fisnucnum*num;
  matrix.put_yx(size,size);
  idold.resize(fisnucnum);
  for(int i=0;i<fisnucnum;i++){
    idold[i].resize(num);
  };
  id.resize(size);
  cout<<"# Covariance data for FP yield of reduced chain are read from file.\n#\n";
  cout<<"#   The number of fissile nuclides : "<<fisnucnum<<"\n";
  cout<<"#   The number of FP nuclides of the reduced chain : "<<num<<"\n";
  for(int i=0;i<fisnucnum;i++){
    int tmp;
    fin>>tmp;
    fisnucname_array[i]=midt.GetName(tmp);
  };
  for(int j=0;j<num;j++){
    int tmp;
    fin>>tmp;
    for(int i=0;i<fisnucnum;i++){
      id[num*i+j]=tmp;
    };
  };

  for(int i=0;i<fisnucnum*num;i++){
    for(int j=0;j<fisnucnum*num;j++){
      real tmp2;
      fin>>tmp2;
      matrix.put_data(i,j,tmp2);
    };
  };
  fin.close();
};

GroupData2D FissionYieldCovarianceCorr::GetCovarianceMatrix(string fisnuc,int fpnum,int *fpid)
{
  bool detect=false;
  GroupData2D ret(fpnum,fpnum);

  int fisnucid=-1;
  for(int i=0;i<fisnucnum;i++){
    if(fisnucname_array[i]==fisnuc)fisnucid=i;
  };
  if(fisnucid==-1){
    cout<<"# Error in FissionYieldCovarianceCorr::GetCovarianceMatrix.\n";
    cout<<"# Cannot find fissile : "<<fisnuc<<"\n";
    exit(0);
  };

  int num=size/fisnucnum;
  vector<int> fporder(fpnum);
  for(int i=0;i<fpnum;i++){
    for(int j=0;j<num;j++){
      if(fpid[i]==id[j])fporder[i]=j;
    };
  };

  int tmp=num*fisnucid;
  for(int i=0;i<fpnum;i++){
    for(int j=0;j<fpnum;j++){
      ret.put_data(i,j,matrix.get_dat(tmp+fporder[i],tmp+fporder[j]));
    };
  };

  return ret;
};

void FissionYieldCovarianceCorr::ShowSelf()
{
  int fpnuc=size/fisnucnum;

  int index=0;
  for(int i=0;i<fisnucnum;i++){

    cout<<"#\n";
    cout<<"# Fissile : "<<fisnucname_array[i]<<"\n";
    cout<<"#\n";
    
    for(int j=0;j<fpnuc;j++){
      cout<<j<<" "<<id[index]<<" "<<sqrt(matrix.get_dat(index,index))<<"\n";
      index++;
    };
    
  };
};
