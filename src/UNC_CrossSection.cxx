#include <cstdlib>
#include "UNC_CrossSection.h"


// +++++++
// Library
// +++++++

int Library::FindData(int mat,int mt)
{
  int size=xs.size();
  for(int i=0;i<size;i++){
    if(mat_list[i]==mat&&mt_list[i]==mt)return i;
  };
  return -1;
};

int Library::FindData2D(int mat,int mt)
{
  int size=xs2d.size();
  for(int i=0;i<size;i++){
    if(mat_list2d[i]==mat&&mt_list2d[i]==mt)return i;
  };
  return -1;
};

void Library::ReadCrossSectionFromFiles(string dirname,int nuc,string *filename)
{
  for(int i=0;i<nuc;i++){
    ReadCrossSectionFromFile(dirname,filename[i]);
  };
};


void Library::ReadCrossSectionFromFile(string dirname,string filename)
{
  dirname.append(filename);

  ifstream fin;
  ofstream fout;
  fin.open(dirname.data(),ios::in);
  if(fin.fail()){
    cout<<"Open failure...in ReadCrossSectionFromFile\n";
    cout<<dirname<<"\n";
    return;
  }

  int grp,mat,mt,dim;
  fin>>grp;

  while(fin.eof()!=true){
    fin>>dim;
    fin>>mat;
    //if(mat==600)mat=625;
    fin>>mt;
    if(fin.eof()){break;};
    int pos;
    if(dim==1){pos=FindData(mat,mt);}
    else{pos=FindData2D(mat,mt);};
    if(pos!=-1){
      cout<<"Caution !!\n";
      cout<<"The following cross section already exists.\n";
      cout<<"  mat : "<<mat<<"   mt : "<<mt<<"\n";
    };
    if(dim==1){
      GroupData1D input_data(grp);
      for(int j=0;j<grp;j++){
        real tmp;
        fin>>tmp;
        input_data.put_data(j,tmp);
      };
      if(pos==-1){
        xs.push_back(input_data);
        mat_list.push_back(mat);
        mt_list.push_back(mt);
      }else{
        xs[pos].put_imax(grp);
        xs[pos].copy(input_data);
      };
    }else {
      GroupData2D input_data(grp,grp);
      for(int j=0;j<grp;j++){
	for(int k=j;k<grp;k++){
          real tmp;
          fin>>tmp;
          input_data.put_data(j,k,tmp);
	};
      };
      if(pos==-1){
        xs2d.push_back(input_data);
        mat_list2d.push_back(mat);
        mt_list2d.push_back(mt);
      }else{
        xs2d[pos].put_yx(grp,grp);
        xs2d[pos].copy(input_data);
      };
    };
  };
};

GroupData1D &Library::GetCrossSection(int mat,int mt)
{
  int tmp=FindData(mat,mt);
  if(tmp==-1){
    cout<<"# Error in Library::GetCrossSection.\n";
    cout<<"# Data of Mat "<<mat<<" and Mt "<<mt<<" cannot be found.\n";
    exit(0);
  };
  return xs[tmp];
};

GroupData2D &Library::GetCrossSection2D(int mat,int mt)
{
  int tmp=FindData2D(mat,mt);
  if(tmp==-1){
    cout<<"Error in Library::GetCrossSection2D.\n";
    exit(0);
  };
  return xs2d[tmp];
};

void Library::GetInformationFromXSLibrary(XSLibrary &xslib,bool thermal_chi)
{
  map<int,LibData>::iterator it=xslib.GetData().begin();
  while(it!=xslib.GetData().end()){
    int ii=(*it).first;
    GetInformationFromLibData(xslib.GetLibData(ii),thermal_chi);
    it++;
  };
};

void Library::GetInformationFromLibData(LibData &ldat,bool thermal_chi)
{
  int mat=ldat.GetMat();
  //int mat=TranslateNuclideIDFromJFS(mat_jfs);
  // for 1D
  int mt1d[]={18,102,452,251,181,2,4,16};
  enum xstype type1d[]={sigf,sigc,nu,mu,chi,sigel,siginel,sign2n};
  for(int i=0;i<8;i++){
    mat_list.push_back(mat);
    mt_list.push_back(mt1d[i]);
    GroupData1D input;
    if(i==4&&thermal_chi&&ldat.ExistChiVector()){
      int grp=ldat.GetChiVector().GetGroup();
      int vid=ldat.GetChiVector().GetVectorID(grp-1);
      input.copy(ldat.GetChiVector().GetChiVector(vid));
    }else{
      input.copy(ldat.GetXSData().GetData1d(type1d[i]));
    }
    xs.push_back(input);
  };

  {
  // P2 component of elastic as 1D
  mat_list.push_back(mat);
  mt_list.push_back(252);

  int grp=ldat.GetGroup();  
  GroupData1D coefp2(grp);
  for(int g=0;g<grp;g++){
    real sige0=0.;
    real sige1=0.;    
    real sige2=0.;
    for(int g2=0;g2<grp;g2++){
      sige0+=ldat.GetXSData().GetData2d(sigel,0).get_dat(g,g2);
      sige1+=ldat.GetXSData().GetData2d(sigel,1).get_dat(g,g2);
      sige2+=ldat.GetXSData().GetData2d(sigel,2).get_dat(g,g2);      
    };
    coefp2.put_data(g,sige2/sige0/5.);
    //cout<<g<<" "<<sige1/sige0/3.<<" "<<sige2/sige0/5.<<"\n";
  };
  xs.push_back(coefp2);
  };

  // for 2D
  int mt2d[]={2,4,16};
  enum xstype type2d[]={sigel,siginel,sign2n};
  for(int i=0;i<3;i++){
    mat_list2d.push_back(mat);
    mt_list2d.push_back(mt2d[i]);
    GroupData2D input;
    input.copy(ldat.GetXSData().GetData2d(type2d[i]));
    if(i==0&&ldat.ExistThermalData()){
      int grp=ldat.GetThScat().GetGrp();
      int srcg=ldat.GetThScat().GetInitSrcGrp();
      int snkg=ldat.GetThScat().GetInitGrp();
      for(int j=srcg;j<grp;j++){
        real s1=0.;
        real s0=0.;
	for(int k=snkg;k<grp;k++){
	  //input.put_data(j,k,ldat.GetThScat().GetData(0,0,j,k));
	  input.put_data(j,k,ldat.GetThScat().GetData(0,j,k,300.));
          s1+=ldat.GetThScat().GetData(1,j,k,300.);
          s0+=ldat.GetThScat().GetData(0,j,k,300.);
	};
        real mu=s1/s0;
	GetCrossSection(mat,251).put_data(j,mu);
      };
    };
    xs2d.push_back(input);
  };
};

// ++++++++++++++++
// LibraryContainer
// ++++++++++++++++

int LibraryContainer::FindData(string name)
{
  int size=libname.size();
  int tmp=-1;
  for(int i=0;i<size;i++){
    if(name==libname[i])tmp=i;
  };
  return tmp;
};

void LibraryContainer::PutLibrary(Library inp,string name)
{
  int tmp=FindData(name);
  if(tmp==-1){
    lib.push_back(inp);
    libname.push_back(name);
  }else{
    lib[tmp]=inp;
  };
};

void LibraryContainer::PutXSLibrary(XSLibrary &xslib,string name,bool thermal_chi)
{
  Library inp;
  inp.GetInformationFromXSLibrary(xslib,thermal_chi);
  PutLibrary(inp,name);
};

Library &LibraryContainer::GetLibrary(string name)
{
  int tmp=FindData(name);
  if(tmp==-1){
    cout<<"You cannot find library "<<name<<"\n";
    exit(0);
  };
  return lib[tmp];
};

// ++++++++++++++++++++++
// CrossSectionCovariance
// ++++++++++++++++++++++

CrossSectionCovariance::CrossSectionCovariance():Covariance()
{
};

CrossSectionCovariance::CrossSectionCovariance(string type):Covariance(type)
{
};

void CrossSectionCovariance::PutMatMtInformation(int imat,int imt)
{
  if(!Variance){
    cout<<"# Error in PutMatMtInformation of CrossSectionCovariance.\n";
    cout<<"# This method is for variance instance.\n";
    exit(0);
  };
  PutMatMtInformation(imat,imt,imat,imt);
};

void CrossSectionCovariance::PutMatMtInformation(int imat1,int imt1,int imat2,int imt2)
{
  if(Variance){
    if(imat1!=imat2||imt1!=imt2){
      cout<<"# Error in PutMatMtInformation of CrossSectionCovariance.\n";
      cout<<"# MAT or MT number is inconsistent.\n";
      cout<<"# This instance is variance.\n";
      exit(0);
    };
  };
  mat1=imat1;
  mt1=imt1;
  mat2=imat2;
  mt2=imt2;
};

void CrossSectionCovariance::ShowMatMt()
{
  cout<<"  "<<mat1<<" : "<<mt1<<"  &  "<<mat2<<" : "<<mt2<<"\n";
};

void CrossSectionCovariance::WriteFile(string dirname,string filename)
{
  dirname.append(filename);

  ofstream fout;
  fout.open(dirname.data(),ios::out);
  if(fout.fail()){
    cout<<"Failed to open the file.\n";
    exit(1);
  };

  fout<<"  1\n";
  fout<<mat1<<"\n";
  fout<<mt1<<"\n";
  fout<<mat2<<"\n";
  fout<<mt2<<"\n";
  fout<<size1<<"\n";
  for(int i=0;i<size1+1;i++){
    fout<<"  0.\n"; // dummy for energy boundary
  };
  for(int i=0;i<2;i++){
    fout<<size1<<"\n";
    for(int i=0;i<size1;i++){
      fout<<val[0].get_dat(i)<<"\n";
    };
  };
  for(int i=0;i<size1*size1;i++){
    fout<<GetCovariance("Relative").get_dat(i)<<"\n";
  };
  fout.close();
};

// ++++++++++++++++++++++++++++++++
//        LibraryCovariance
// ++++++++++++++++++++++++++++++++

int LibraryCovariance::FindData(int mat1,int mt1,int mat2,int mt2)
{
  int size=xscov.size();
  for(int i=0;i<size;i++){
    int mats1=xscov[i].GetMat1();
    int mts1=xscov[i].GetMt1();
    int mats2=xscov[i].GetMat2();
    int mts2=xscov[i].GetMt2();
    if(mat1==mats1&&mt1==mts1&&mat2==mats2&&mt2==mts2)return i;
  };
  return -1;
};

void LibraryCovariance::ReadCovarianceFromFile(string dirname, string filename)
{
  dirname.append(filename);

  ifstream fin;
  ofstream fout;
  fin.open(dirname.data(),ios::in);
  if(fin.fail()){
    cout<<"# Open failure ...\n";
    cout<<"# "<<dirname<<"\n";
    exit(1);
  }

  CrossSectionCovariance inp;

  int number,mat1,mt1,mat2,mt2,g;
  fin>>number;
  for(int i=0;i<number;i++){

    fin>>mat1;
    fin>>mt1;
    fin>>mat2;
    fin>>mt2;
    if(mat1<10000)mat1=midt.GetMATIDFromENDFID(mat1);
    if(mat2<10000)mat2=midt.GetMATIDFromENDFID(mat2);

    if(mt1==251251)mt1=252;
    if(mt2==251251)mt2=252;    

    bool exist=false;
    int idp=FindData(mat1,mt1,mat2,mt2);
    if(idp!=-1){
      cout<<"# CovarianceData already exists ... updated\n";
      cout<<"#  Mat1:"<<mat1<<" Mt1:"<<mt1<<" Mat2:"<<mat2<<" Mt2:"<<mt2<<"\n";
      exist=true;
    };

    fin>>g;
    real *en=new real[g+1];
    for(int i=0;i<g+1;i++){
      fin>>en[i];
    };
 
    fin>>g;
    real *xs1=new real[g];
    for(int i=0;i<g;i++){
      fin>>xs1[i];
    };

    fin>>g;
    real *xs2=new real[g];
    for(int i=0;i<g;i++){
      fin>>xs2[i];
    };
    real *cov=new real[g*g];
    for(int i=0;i<g*g;i++){
      fin>>cov[i];
    };

    if(mat1==mat2&&mt1==mt2){
      inp.PutType("variance");
      inp.PutSize(g);
      inp.PutVal(xs1);
      inp.PutCov(cov,"relative");
      inp.PutMatMtInformation(mat1,mt1);
    }else{
      inp.PutType("covariance");
      inp.PutSize(g,g);
      inp.PutVal(0,xs1);
      inp.PutVal(1,xs2);
      inp.PutCov(cov,"relative");
      inp.PutMatMtInformation(mat1,mt1,mat2,mt2);
    };

    if(!exist){
      xscov.push_back(inp);
    }else{
      xscov[idp]=inp;
    };

    delete [] en;
    delete [] xs1;
    delete [] xs2;
    delete [] cov;
  };
};

void LibraryCovariance::ReadCovarianceFromFileOverWriting(string dirname, string filename,real emax)
{
  dirname.append(filename);

  ifstream fin;
  ofstream fout;
  fin.open(dirname.data(),ios::in);
  if(fin.fail()){
    cout<<"# Open failure ...\n";
    cout<<"# "<<dirname<<"\n";
    exit(1);
  }

  CrossSectionCovariance inp;

  int number,mat1,mt1,mat2,mt2,g;
  real tmp;

  fin>>number;
  for(int i=0;i<number;i++){

    fin>>mat1;
    fin>>mt1;
    fin>>mat2;
    fin>>mt2;
    if(mat1<10000)mat1=midt.GetMATIDFromENDFID(mat1);
    if(mat2<10000)mat2=midt.GetMATIDFromENDFID(mat2);

    int idp=FindData(mat1,mt1,mat2,mt2);
    if(idp==-1){
      cout<<"# Error in LibraryCovariance::ReadCovarianceFromFileOverWriting.\n";
      cout<<"# CovarianceData does not exist.\n";
      cout<<"#  Mat1:"<<mat1<<" Mt1:"<<mt1<<" Mat2:"<<mat2<<" Mt2:"<<mt2<<"\n";
      exit(0);
    };

    fin>>g;
    int gtop=-2;

    real *xs1=new real[g];
    real *xs2=new real[g];
    real *cov=new real[g*g];
 
    for(int i=0;i<g+1;i++){
      fin>>tmp;
      if(gtop==-2&&tmp<emax)gtop=i;
    };
    if(gtop==-2){
      cout<<"# Error in LibraryCovariance::ReadCovarianceFromFileOverWriting.\n";
      cout<<"# emax is inappropriate.\n";
      exit(0);
    };
    if(gtop==-1)gtop=0;

    fin>>g;
    for(int i=0;i<g;i++){
      xs1[i]=xscov[idp].GetVal(0).get_dat(i);
      fin>>tmp;
      if(gtop<=i)xs1[i]=tmp;
    };

    fin>>g;
    for(int i=0;i<g;i++){
      fin>>tmp;
      if(mt1!=mt2){
        xs2[i]=xscov[idp].GetVal(1).get_dat(i);
        if(gtop<=i)xs2[i]=tmp;
      };
    };

    GroupData2D covorg=xscov[idp].GetCovariance();

    int ii=0;
    for(int i=0;i<g;i++){
      for(int j=0;j<g;j++){
        cov[ii]=covorg.get_dat(i,j);
        fin>>tmp;
	if(gtop<=i&&gtop<=j){
          cov[ii]=tmp;
	};
	ii++;
      };
    };

    if(mt1==mt2){
      xscov[idp].PutVal(xs1);
    }else{
      xscov[idp].PutVal(0,xs1);
      xscov[idp].PutVal(1,xs2);
    };
    xscov[idp].PutCov(cov,"Relative");

    delete [] xs1;
    delete [] xs2;
    delete [] cov;

  };
};

void LibraryCovariance::ReadCovarianceFromFile(string dirname, string filename,MATIDTranslator &midt)
{
  dirname.append(filename);

  ifstream fin;
  ofstream fout;
  fin.open(dirname.data(),ios::in);
  if(fin.fail()){
    cout<<"# Open failure ...\n";
    cout<<"# "<<dirname<<"\n";
    exit(1);
  }

  CrossSectionCovariance inp;

  int number,mat1,mt1,mat2,mt2,g;
  fin>>number;
  for(int i=0;i<number;i++){

    fin>>mat1;
    fin>>mt1;
    fin>>mat2;
    fin>>mt2;

    mat1=midt.GetMATIDFromENDFID(mat1);
    mat2=midt.GetMATIDFromENDFID(mat2);

    bool exist=false;
    int idp=FindData(mat1,mt1,mat2,mt2);
    if(idp!=-1){
      cout<<"CovarianceData already exists ... updated\n";
      cout<<"  Mat1:"<<mat1<<" Mt1:"<<mt1<<" Mat2:"<<mat2<<" Mt2:"<<mt2<<"\n";
      exist=true;
    };

    fin>>g;
    real *en=new real[g+1];
    for(int i=0;i<g+1;i++){
      fin>>en[i];
    };
 
    fin>>g;
    real *xs1=new real[g];
    for(int i=0;i<g;i++){
      fin>>xs1[i];
    };

    fin>>g;
    real *xs2=new real[g];
    for(int i=0;i<g;i++){
      fin>>xs2[i];
    };
    real *cov=new real[g*g];
    for(int i=0;i<g*g;i++){
      fin>>cov[i];
    };

    if(mat1==mat2&&mt1==mt2){
      inp.PutType("variance");
      inp.PutSize(g);
      inp.PutVal(xs1);
      inp.PutCov(cov,"relative");
      inp.PutMatMtInformation(mat1,mt1);
    }else{
      inp.PutType("covariance");
      inp.PutSize(g,g);
      inp.PutVal(0,xs1);
      inp.PutVal(1,xs2);
      inp.PutCov(cov,"relative");
      inp.PutMatMtInformation(mat1,mt1,mat2,mt2);
    };

    if(!exist){
      xscov.push_back(inp);
    }else{
      xscov[idp]=inp;
    };

    delete [] en;
    delete [] xs1;
    delete [] xs2;
    delete [] cov;
  };
};

CrossSectionCovariance &LibraryCovariance::GetCrossSectionCovariance(int mat1,int mt1,int mat2,int mt2)
{
  int pos=FindData(mat1,mt1,mat2,mt2);
  if(pos==-1){
    cout<<"# Error in LibraryCovariance::GetCrossSectionCovariance.\n";
    cout<<"# Covariance data between ("<<mat1<<","<<mt1<<") and ";
    cout<<"("<<mat2<<","<<mt2<<") cannot be found.\n";
    exit(0);
  };
  return xscov[pos];
};

int LibraryCovariance::GetSize()
{
  int ret=xscov.size();
  return ret;
};

void LibraryCovariance::ShowSelf()
{
  int size=xscov.size();
  for(int i=0;i<size;i++){
    int mats1=xscov[i].GetMat1();
    int mts1=xscov[i].GetMt1();
    int mats2=xscov[i].GetMat2();
    int mts2=xscov[i].GetMt2();
    cout<<i<<"/"<<size<<" : "<<mats1<<" "<<mts1<<"  :  "<<mats2<<" "<<mts2<<"\n";
  };
};

void LibraryCovariance::CalCovarianceForXSDifference(int mat,int mt1,int mt2)
{
  if(FindData(mat,mt1,mat,mt1)==-1||FindData(mat,mt2,mat,mt2)==-1){
    cout<<"# Covariance matrix cannot be founded in LibraryCovariance.\n";
    return;
  };
  
  GroupData1D xs1=GetCrossSectionCovariance(mat,mt1,mat,mt1).GetVal();
  GroupData1D xs2=GetCrossSectionCovariance(mat,mt2,mat,mt2).GetVal();

  GroupData2D acov1=GetCrossSectionCovariance(mat,mt1,mat,mt1).GetCovariance("Absolute");
  GroupData2D acov2=GetCrossSectionCovariance(mat,mt2,mat,mt2).GetCovariance("Absolute");

  GroupData2D acov12;
  bool mt1mt2=false;
  if(FindData(mat,mt1,mat,mt2)!=-1){
    mt1mt2=true;
    acov12=GetCrossSectionCovariance(mat,mt1,mat,mt2).GetCovariance("Absolute");
  };
  GroupData2D acov21;
  bool mt2mt1=false;
  if(FindData(mat,mt2,mat,mt1)!=-1){
    mt2mt1=true;
    acov21=GetCrossSectionCovariance(mat,mt2,mat,mt1).GetCovariance("Absolute");
  };
  if(mt1mt2&&mt2mt1){
    cout<<"#???\n";
    exit(0);
  };

  int grp=xs1.get_imax();

  GroupData2D covdat(grp,grp);
  covdat.set_zero();
  for(int g=0;g<grp;g++){
    for(int g2=0;g2<grp;g2++){
      covdat.add_data(g,g2,acov1.get_dat(g,g2));
      covdat.add_data(g,g2,acov2.get_dat(g,g2));
      if(mt1mt2){
	covdat.add_data(g,g2,-acov12.get_dat(g,g2));
	covdat.add_data(g2,g,-acov12.get_dat(g2,g));
      };
      if(mt2mt1){
	covdat.add_data(g,g2,-acov21.get_dat(g,g2));
	covdat.add_data(g2,g,-acov21.get_dat(g2,g));
      };
    };
  };

  Covariance cov("variance");
  cov.PutSize(grp);
  xs1=xs1-xs2;
  cov.PutVal(0,xs1);
  cov.PutCov(covdat,"Absolute");
  cov.GetStandardDeviation("Relative").show_self();
};

void LibraryCovariance::ChangeMTNumber(int mat,int mt1,int mt2)
{
  int sz=xscov.size();
  for(int i=0;i<sz;i++){
    int mat1=xscov[i].GetMat1();
    int mta=xscov[i].GetMt1();
    int mat2=xscov[i].GetMat2();
    int mtb=xscov[i].GetMt2();
    if(midt.GetMATIDFromENDFID(mat1)==mat){
      if(mta==mt1){
        mta=mt2;
      }else if(mta==mt2){
	mta=mt1;
      };
    };
    if(midt.GetMATIDFromENDFID(mat2)==mat){
      if(mtb==mt1){
        mtb=mt2;
      }else if(mtb==mt2){
	mtb=mt1;
      };
    };
    xscov[i].PutMatMtInformation(mat1,mta,mat2,mtb);
  };
};

void LibraryCovariance::SetNoCorrelation()
{  
  int sz=xscov.size();
  for(int i=0;i<sz;i++){
    xscov[i].SetNoCorrelation();
  };
};

void LibraryCovariance::AddCovarianceData(CrossSectionCovariance inp)
{
  int mat1=inp.GetMat1();
  int mat2=inp.GetMat2();
  int mt1=inp.GetMt1();
  int mt2=inp.GetMt2();
  if(FindData(mat1,mt1,mat2,mt2)!=-1){
    cout<<"# Error in LibraryCovariance::AddCovarianceData.\n";
    cout<<"# Covariance data between ("<<mat1<<","<<mt1<<") and ";
    cout<<"("<<mat2<<","<<mt2<<") already exists.\n";
    exit(0);
  };
  xscov.push_back(inp);
};

void LibraryCovariance::AddConstantRelativeCovariance(int mat, int mt, int grp, real rel_stdev, bool full_correlation)
{
  // *** CAUTION !! : Temporal implementation ***
  //
  // Covariance data should be registered with their cross section
  // to calculate both absolute and relative covariance data,
  // but this method does NOT use proper cross section data.
  // Instead of the proper cross section data, cross section
  // values are assumed to unity, so absolute and relative
  // covariance data are identical to each other.
  
  GroupData2D covdat(grp,grp);
  GroupData1D xsdummy(grp); // dummy cross section
  real tmp=rel_stdev*rel_stdev;
  for(int g=0;g<grp;g++){
    xsdummy.put_data(g,1.);
    for(int g2=0;g2<grp;g2++){
      covdat.put_data(g,g2,tmp);
    };
  };

  if(!full_correlation){
    for(int g=0;g<grp;g++){
      for(int g2=0;g2<grp;g2++){
	if(g!=g2)covdat.put_data(g,g2,0.);
      };
    };
  };

      

  CrossSectionCovariance cov("variance");
  cov.PutSize(grp);
  cov.PutVal(0,xsdummy);
  cov.PutCov(covdat,"Relative");
  cov.PutMatMtInformation(mat,mt);
  AddCovarianceData(cov);
  
};

void LibraryCovariance::ResetDiagonalElement(int mat, int mt, int grp, real rel_stdev)
{
  GetCrossSectionCovariance(mat,mt).GetCov().put_data(grp,grp,rel_stdev*rel_stdev);
};

void LibraryCovariance::ReplaceNutCovByNupCov(int mat)
{
  int pos=FindData(mat,456,mat,456);
  if(pos==-1){
    cout<<"# Error in LibraryCovariance::ReplaceNutCovByNupCov.\n";
    cout<<"# Nu_p covariance of materal "<<mat<<" cannot be found.\n";
    exit(0);
  };

  int pos2=FindData(mat,452,mat,452);
  if(pos2!=-1){
    cout<<"# Error in LibraryCovariance::ReplaceNutCovByNupCov.\n";
    cout<<"# Nu_t covariance of materal "<<mat<<" is already defined.\n";
    exit(0);
  };

  CrossSectionCovariance tmp=GetCrossSectionCovariance(mat,456);
  tmp.PutMatMtInformation(mat,452,mat,452);
  xscov.push_back(tmp);
};

void LibraryCovariance::CovarianceForSummedQuantity(int matid, int mtnum, int *mtid, int newmt)
{
  int grp;
  GroupData1D val_new;
  for(int i=0;i<mtnum;i++){
    int tmp=FindData(matid,mtid[i],matid,mtid[i]);
    if(tmp==-1){
      cout<<"# Error in LibraryCovariance::CovarianceForSummedQuantity.\n";
      cout<<"# The covariance data "<<matid<<"/"<<mtid[i]<<" does not exist.\n";
      exit(0);
    };
    GroupData1D val=GetCrossSectionCovariance(tmp).GetVal();
    //val.show_self();
    if(i==0){
      grp=val.get_imax();
      val_new=val;
    }else{
      val_new=val_new+val;
    };
  };
  //val_new.show_self();

  GroupData2D cov_new(grp,grp);
  cov_new.set_zero();
  int sz=xscov.size();
  for(int i=0;i<sz;i++){
    int mat1=xscov[i].GetMat1();
    int mat2=xscov[i].GetMat2();
    if(mat1==matid&&mat2==matid){
      int mt1=xscov[i].GetMt1();
      int mt2=xscov[i].GetMt2();
      int mt1pos=-1;
      int mt2pos=-1;
      for(int j=0;j<mtnum;j++){
        if(mtid[j]==mt1)mt1pos=j;
	if(mtid[j]==mt2)mt2pos=j;
      };
      if(mt1pos!=-1&&mt2pos!=-1){
        cov_new=cov_new+xscov[i].GetCovariance("Absolute");
	if(mt1pos!=mt2pos){
          GroupData2D tmp=xscov[i].GetCovariance("Absolute").GetTransposedMatrix();
	  cov_new=cov_new+tmp;
	};
      };
      /*
      if(mt1pos==mt2pos&&mt1pos!=-1){
        cov_new=cov_new+xscov[i].GetCovariance("Absolute");
      };
      */
    };
  };

  CrossSectionCovariance inp;
  inp.PutType("variance");
  inp.PutSize(grp);
  inp.PutVal(0,val_new);
  inp.PutCov(cov_new,"absolute");
  inp.PutMatMtInformation(matid,newmt);

  int tmp=FindData(matid,newmt,matid,newmt);
  if(tmp==-1){
    xscov.push_back(inp);
  }else{
    xscov[tmp]=inp;
  };

  //

  vector<int> pos1;
  vector<int> pos2;
  for(int i=0;i<sz;i++){
    int mat1=xscov[i].GetMat1();
    int mat2=xscov[i].GetMat2();
    if(mat1==matid&&mat2==matid){
      int mt1=xscov[i].GetMt1();
      int mt2=xscov[i].GetMt2();
      int mt1pos=-1;
      int mt2pos=-1;
      for(int j=0;j<mtnum;j++){
        if(mtid[j]==mt1)mt1pos=j;
	if(mtid[j]==mt2)mt2pos=j;
      };
      if(mt1==newmt&&mt2pos==-1)pos1.push_back(i);
      if(mt2==newmt&&mt1pos==-1)pos2.push_back(i);
    };
  };

  int sz2=pos1.size();
  for(int i=0;i<sz2;i++){
    int p=pos1[i];
    int mt_pair=xscov[p].GetMt2();
    for(int j=0;j<sz;j++){
      int mat1=xscov[j].GetMat1();
      int mat2=xscov[j].GetMat2();
      if(p!=j&&mat1==matid&&mat2==matid){
        int mt1=xscov[j].GetMt1();
	int mt2=xscov[j].GetMt2();
        int mt1pos=-1;
        int mt2pos=-1;
        for(int k=0;k<mtnum;k++){
          if(mtid[k]==mt1)mt1pos=k;
	  if(mtid[k]==mt2)mt2pos=k;
        };
	if(mt1==mt_pair&&mt2pos!=-1){
          GroupData2D tmp=xscov[p].GetCovariance("Absolute");
          GroupData2D tmp2=xscov[j].GetCovariance("Absolute").GetTransposedMatrix();
	  tmp=tmp+tmp2;
	  xscov[p].PutCov(tmp,"Absolute");
	};
	if(mt2==mt_pair&&mt1pos!=-1){
          GroupData2D tmp=xscov[p].GetCovariance("Absolute");
          GroupData2D tmp2=xscov[j].GetCovariance("Absolute");
	  tmp=tmp+tmp2;
	  xscov[p].PutCov(tmp,"Absolute");
	};
      };
    };
    GroupData1D tmp=GetCrossSectionCovariance(matid,newmt).GetVal();
    xscov[p].PutVal(0,tmp);
  };

  sz2=pos2.size();
  for(int i=0;i<sz2;i++){
    int p=pos2[i];
    int mt_pair=xscov[p].GetMt1();
    //cout<<i<<" "<<mt_pair<<"\n";
    for(int j=0;j<sz;j++){
      int mat1=xscov[j].GetMat1();
      int mat2=xscov[j].GetMat2();
      if(p!=j&&mat1==matid&&mat2==matid){
        int mt1=xscov[j].GetMt1();
	int mt2=xscov[j].GetMt2();
        int mt1pos=-1;
        int mt2pos=-1;
        for(int k=0;k<mtnum;k++){
          if(mtid[k]==mt1)mt1pos=k;
	  if(mtid[k]==mt2)mt2pos=k;
        };
	if(mt2==mt_pair&&mt1pos!=-1){
          GroupData2D tmp=xscov[p].GetCovariance("Absolute");
          GroupData2D tmp2=xscov[j].GetCovariance("Absolute").GetTransposedMatrix();
	  tmp=tmp+tmp2;
	  xscov[p].PutCov(tmp,"Absolute");
	};
	if(mt1==mt_pair&&mt2pos!=-1){
	  //cout<<"   "<<j<<" "<<mt2pos<<" "<<mtid[mt2pos]<<"\n";
          GroupData2D tmp=xscov[p].GetCovariance("Absolute");
          GroupData2D tmp2=xscov[j].GetCovariance("Absolute");
	  //if(i==0){
	  //cout<<tmp.get_dat(69,69)<<" "<<tmp2.get_dat(69,69)<<"\n";
	  //};
	  tmp=tmp+tmp2;
	  xscov[p].PutCov(tmp,"Absolute");
	};
      };
    };
    GroupData1D tmp=GetCrossSectionCovariance(matid,newmt).GetVal();
    xscov[p].PutVal(1,tmp);
  };

};

GroupData2D LibraryCovariance::GetCovarianceMatrixMultipleParameters(int num, int *matid, int *mtid, bool corr)
{
  vector<int> check;
  for(int i=0;i<num;i++){
    for(int j=i;j<num;j++){
      int tmp=FindData(matid[i],mtid[i],matid[j],mtid[j]);
      check.push_back(tmp);
      if(i==j&&tmp==-1){
        cout<<"# Error in LibraryCovariance::GetCovarianceMatrixMultipleParameters.\n";
        cout<<"# No covariance data for "<<matid[i]<<" "<<mtid[i]<<"\n";
        exit(0);
      };
    };
  };

  int group=GetCrossSectionCovariance(check[0]).GetSize();

  int tot=check.size();
  for(int i=1;i<tot;i++){
    if(check[i]!=-1){
      bool var=GetCrossSectionCovariance(check[i]).GetVariance();
      if(var){
	int tmp=GetCrossSectionCovariance(check[i]).GetSize();
	if(tmp!=group){
          cout<<"# Error in LibraryCovariance::GetCovarianceMatrixMultipleParameters.\n";
          cout<<"# Inconsistent-size data.\n";
          exit(0);
	};
      }else{
	int tmp1=GetCrossSectionCovariance(check[i]).GetSize1();
	int tmp2=GetCrossSectionCovariance(check[i]).GetSize2();
	if(tmp1!=group||tmp2!=group){
          cout<<"# Error in LibraryCovariance::GetCovarianceMatrixMultipleParameters.\n";
          cout<<"# Inconsistent-size data.\n";
          exit(0);
	};
      };
    };
  };

  GroupData2D COV(group*num,group*num);

  int index=0;
  for(int j=0;j<num;j++){
    for(int k=j;k<num;k++){
      GroupData2D covmat=GetCrossSectionCovariance(check[index]).GetCovariance("Relative"); 
      COV.PasteGroupData2D(group*j,group*k,covmat);
      if(k!=j){
        covmat.Transposition();
        COV.PasteGroupData2D(group*k,group*j,covmat);
      };
      index++;
    };
  };

  if(corr){
    vector<real> stdev_inv(group*num);
    for(int i=0;i<group*num;i++){
      stdev_inv[i]=1./(sqrt(COV.get_dat(i,i)));
    };
    for(int i=0;i<group*num;i++){
      for(int j=0;j<group*num;j++){
	real tmp=COV.get_dat(i,j);
	COV.put_data(i,j,tmp*stdev_inv[i]*stdev_inv[j]);
      };
    };
  };
  
  return COV;
};

void LibraryCovariance::CheckCorrelationMatrix(int mat1, int mt1, int mat2, int mt2)
{
  int pos1=FindData(mat1,mt1,mat1,mt1);
  int pos2=FindData(mat2,mt2,mat2,mt2);
  int pos3=FindData(mat1,mt1,mat2,mt2);
  int pos4=FindData(mat2,mt2,mat1,mt1);

  if(pos1==-1){
    cout<<"# Error in LibraryCovariance::CheckCorrelationMatrix.\n";
    cout<<"# No variance data is defined for (MAT/MT)="<<mat1<<"/"<<mt1<<"\n";
    exit(0);
  };
  if(pos2==-1){
    cout<<"# Error in LibraryCovariance::CheckCorrelationMatrix.\n";
    cout<<"# No variance data is defined for (MAT/MT)="<<mat2<<"/"<<mt2<<"\n";
    exit(0);    
  };
  if(pos3==-1&&pos4==-1){
    cout<<"# Error in LibraryCovariance::CheckCorrelationMatrix.\n";
    cout<<"# No covariance data is defined between ("<<mat1<<"/"<<mt1<<") and ("<<mat2<<"/"<<mt2<<")\n";
    exit(0);    
  };

  GroupData2D covmat1=GetCrossSectionCovariance(mat1,mt1).GetCovariance();
  GroupData2D covmat2=GetCrossSectionCovariance(mat2,mt2).GetCovariance();  
  GroupData2D corrmat;
  if(pos3!=-1){
    corrmat=GetCrossSectionCovariance(mat1,mt1,mat2,mt2).GetCovariance();
  }else{
    corrmat=GetCrossSectionCovariance(mat2,mt2,mat1,mt1).GetCovariance();    
  };

  int ig=corrmat.get_y();
  for(int g=0;g<ig;g++){
    real tmp1=sqrt(covmat1.get_dat(g,g));
    for(int g2=0;g2<ig;g2++){
      real tmp2=sqrt(covmat2.get_dat(g2,g2));
      if(tmp1>0.&&tmp2>0.){
	real tmp=corrmat.get_dat(g,g2);
	if(pos3==-1)tmp=corrmat.get_dat(g2,g);
        real cor=tmp/tmp1/tmp2;
	if(fabs(cor)>1.00000001){
	  cout<<g<<" "<<g2<<" "<<cor<<"\n";
	};
      };
    };
  };
};

void LibraryCovariance::ModifyingCorrelation(int mat1, int mt1)
{
  int pos1=FindData(mat1,mt1,mat1,mt1);

  if(pos1==-1){
    cout<<"# Error in LibraryCovariance::ModifyingCorrelation.\n";
    cout<<"# No variance data is defined for (MAT/MT)="<<mat1<<"/"<<mt1<<"\n";
    exit(0);
  };

  GroupData2D covmat1=GetCrossSectionCovariance(mat1,mt1).GetCovariance();

  int ig=covmat1.get_y();
  for(int g=0;g<ig;g++){
    real tmp1=sqrt(covmat1.get_dat(g,g));
    for(int g2=0;g2<ig;g2++){
      real tmp2=sqrt(covmat1.get_dat(g2,g2));
      if(tmp1>0.&&tmp2>0.){
	real corr=covmat1.get_dat(g,g2)/tmp1/tmp2;
        real org=covmat1.get_dat(g,g2);

	//GetCrossSectionCovariance(mat1,mt1).GetCov().put_data(g,g2,0.);

	if(fabs(corr)>1){
	  cout<<"# Absolute value of correlation is greater than unity.\n";
	  //GetCrossSectionCovariance(mat1,mt1).GetCov().put_data(g,g2,org/corr);
	  GetCrossSectionCovariance(mat1,mt1).GetCov().put_data(g,g2,0.);
	};

	/*
	if(corr>1.){
	  GetCrossSectionCovariance(mat1,mt1).GetCov().put_data(g,g2,tmp1*tmp2);
	};
	if(corr<-1.){
	  GetCrossSectionCovariance(mat1,mt1).GetCov().put_data(g,g2,-tmp1*tmp2);
	};
	*/
	
      };
    };
  };
};
