#include <cstdlib>
#include "UNC_Sensitivity.h"

SensitivityContainer::SensitivityContainer(ParameterList* inp)
{
  plist=inp;
  int size=plist->GetSize();
  sens.resize(size);
};

void SensitivityContainer::ReadSensitivityDataFromFile
(string dirname,string filename,string core,string chara,int step)
{
  int tmp=plist->FindData(core,chara,step);
  if(tmp==-1){
    cout<<"Error in SensitivityContainer::ReadSensitivityFile.\n";
    cout<<"There is no parameter in parameter list.\n";
    cout<<"  "<<core<<" "<<chara<<" "<<step<<"\n";
    exit(0);
  };
  sens[tmp].ReadFile(dirname,filename);
};

real &SensitivityContainer::GetSensitivity0D
(int mat,int mt,string core,string chara,int step)
{
  return GetSensitivityData(core,chara,step).GetSensitivity0D(mat,mt);
};

GroupData1D &SensitivityContainer::GetSensitivity1D
(int mat,int mt,string core,string chara,int step)
{
  return GetSensitivityData(core,chara,step).GetSensitivity1D(mat,mt);
};

GroupData2D &SensitivityContainer::GetSensitivity2D
(int mat,int mt,string core,string chara,int step)
{
  return GetSensitivityData(core,chara,step).GetSensitivity2D(mat,mt);
};

SensitivityData &SensitivityContainer::GetSensitivityData(string core,string chara,int step)
{
  int pos=plist->FindData(core,chara,step);
  if(pos==-1){
    cout<<"Error in SensitivityContainer::GetSensitivityData.\n";
    cout<<"  core :"<<core<<"\n";
    cout<<"  chara:"<<chara<<"\n";
    cout<<"  step : "<<step<<"\n";
    exit(0);
  };
  return sens[pos];
};

void SensitivityContainer::TransformToConstrainedSensitivityForChi(Library &lib)
{
  int size=sens.size();
  for(int i=0;i<size;i++){
    for(int j=0;j<sens[i].GetSize1D();j++){  
      int mat=sens[i].GetMatList1D(j);
      int mt=sens[i].GetMtList1D(j);
      if(mt==181){
	int tmp=lib.FindData(mat,mt);
	if(tmp!=-1){
          sens[i].TransformToConstrainedSensitivity(mat,mt,lib.GetCrossSection(tmp));
	};
      };
    };
  };
};


void SensitivityContainer::SetRelativeSensitivity(Library &lib)
{
  int size=sens.size();
  for(int i=0;i<size;i++){
    // for 1D-xs
    for(int j=0;j<sens[i].GetSize1D();j++){
      int mat=sens[i].GetMatList1D(j);
      int mt=sens[i].GetMtList1D(j);
      if(mt==251||mt==181||mt==252){ // [252] is for P1 coef of (n,n)
	//if(mt==251||mt==181){	
	int tmp=lib.FindData(mat,mt);
	if(tmp!=-1){
	  sens[i].GetSensitivity1D(j).copy
	    (sens[i].GetSensitivity1D(j).mult(lib.GetCrossSection(tmp)));
	  /*
	}else{
	  cout<<"# Error in SensitivityContainer::SetRelativeSensitivity.\n";
	  cout<<"# There is no 1D cross section data for (mat/mt)= "<<mat<<"/"<<mt<<"\n";
	  exit(0);
	  */
	};
	/*
      }else if(mt==252){
	int tmp=lib.FindData(mat,2);
	if(tmp!=-1){
	  
	  sens[i].GetSensitivity1D(j).copy
	    (sens[i].GetSensitivity1D(j).mult(lib.GetCrossSection(tmp)));
	};
	*/
      };
    };
    // for 2D-xs
    for(int j=0;j<sens[i].GetSize2D();j++){
      int mat=sens[i].GetMatList2D(j);
      int mt=sens[i].GetMtList2D(j);
      if(mt==4||mt==16){
	int tmp=lib.FindData2D(mat,mt);
	if(tmp!=-1){
	  int imax=sens[i].GetSensitivity2D(j).get_y();
	  for(int k1=0;k1<imax;k1++){
	    for(int k2=0;k2<imax;k2++){
  	      real xs=lib.GetCrossSection2D(tmp).get_dat(k1,k2);
	      real ss=sens[i].GetSensitivity2D(j).get_dat(k1,k2);
	      sens[i].GetSensitivity2D(j).put_data(k1,k2,ss*xs);
	    };
	  };
	  /*
	}else{
	  cout<<"# Error in SensitivityContainer::SetRelativeSensitivity.\n";
	  cout<<"# There is no 2D cross section data for (mat/mt)= "<<mat<<"/"<<mt<<"\n";
	  exit(0);
	  */
	};
      };
    };
  };

};

void SensitivityContainer::Produce1DSensitivity()
{
  int size=sens.size();
  for(int i=0;i<size;i++){
    for(int j=0;j<sens[i].GetSize2D();j++){
      int mat=sens[i].GetMatList2D(j);
      int mt=sens[i].GetMtList2D(j);
      if(mt==2||mt==4||mt==16){
	int imax=sens[i].GetSensitivity1D(j).get_imax();
	GroupData1D input_data(imax);
	for(int k1=0;k1<imax;k1++){
	  real sum=0.;
	  for(int k2=0;k2<imax;k2++){
	    sum+=sens[i].GetSensitivity2D(j).get_dat(k1,k2);
	  };
	  input_data.put_data(k1,sum);
	};
	sens[i].PutSensitivity1D(mat,mt,input_data);
      };
    };
  };
};

void SensitivityContainer::GetDiscreteInelasticSensitivity(string mdir, string filename, int matid)
{
  int grp=sens[0].GetGroup();
  int size=sens.size();

  int inmie=64;
  string name_reaction[]={
   "n01","n02","n03","n04","n05", "n06","n07","n08","n09","n10",
   "n11","n12","n13","n14","n15", "n16","n17","n18","n19","n20",
   "n21","n22","n23","n24","n25", "n26","n27","n28","n29","n30",
   "n31","n32","n33","n34","n35", "n36","n37","n38","n39","n40",
   "n41","n42","n43","n44","n45", "n46","n47","n48","n49","n50",
   "n51","n52","n53","n54","n55", "n56","n57","n58","n59","ncn",
   "nna","nnp","nnd","nnt"
  };

  mdir=mdir+filename;
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the inelastic matrix file.\n";
    cout<<"# File name is "<<mdir<<"\n";
    exit(0);
  };

  vector<GroupData2D> ine(inmie);
  for(int i=0;i<inmie;i++){
    ine[i].put_yx(grp,grp);
    ine[i].set_zero();
  };

  while(!fin.eof()){
    int id;
    fin>>id;
    if(id==-1)break;
    for(int g=0;g<grp;g++){
      int i1,i2;
      fin>>i1;
      fin>>i2;
      for(int g2=i1;g2>=i2;g2--){
        real tmp;
        fin>>tmp;
	ine[id-1].put_data(g2-1,g,tmp);
      };
    };

    if((id>=1&&id<=59)||id==60){
      int mtn=50+id;
      if(id==60)mtn=91;
      for(int i=0;i<size;i++){
        for(int j=0;j<sens[i].GetSize2D();j++){
          int mat=sens[i].GetMatList2D(j);
          int mt=sens[i].GetMtList2D(j);
          if(mat==matid&&mt==4){
	    GroupData1D sns_vec(grp);
            for(int g=0;g<grp;g++){
  	      real tmp=0.;
	      for(int g2=g;g2<grp;g2++){
		real xs=ine[id-1].get_dat(g,g2);
		if(xs!=0.){
  	  	  real sv=sens[i].GetSensitivity2D(j).get_dat(g,g2);
		  tmp+=sv*xs;
		};
	      };
	      sns_vec.put_data(g,tmp);
	    };
	    sens[i].PutSensitivity1D(mat,mtn,sns_vec);
	  };
	};
      };
    };
  };

};


void SensitivityContainer::ChangeMTNumber(int mat,int mt1,int mt2)
{
  int sz=sens.size();
  for(int i=0;i<sz;i++){
    sens[i].ChangeMatMt(mat,mt1,mat,mt2);
  };
};
