#include <cstdlib>
#include "CartCore.h"

using namespace std;

// *********************************************************
//    Assembly
// *********************************************************

/*
Assembly &Assembly::operator=(Assembly &sec){
  if(sec.exist){
    this->exist=true;
    this->zdiv=sec.zdiv;
    this->MaterialID.resize(this->zdiv);
    this->zdist.resize(this->zdiv);
    this->MaterialName.resize(this->zdiv);
    for(int i=0;i<this->zdiv;i++){
      this->MaterialID[i]=sec.MaterialID[i];
      this->zdist[i]=sec.zdist[i];
      this->MaterialName[i]=sec.MaterialName[i];
    };
    this->AsmName=sec.AsmName;
  };
  return *this;
};
*/

void Assembly::Init(int zdivi,string *zmapi,real *zdisti,string inp,string type)
{
  vector<string> zmapt;
  vector<real> zdistt;
  zmapt.resize(zdivi);
  zdistt.resize(zdivi);
  for(int i=0;i<zdivi;i++){
    zmapt[i]=zmapi[i];
    zdistt[i]=zdisti[i];
  };
  Init(zdivi,zmapt,zdistt,inp,type);
}

void Assembly::Init(int zdivi,vector<string> zmapi, vector<real> zdisti,string inp,string type)
{
  if(type!="cumulative"&&type!="width"){
    cout<<"Error in Init method of Assembly class.\n";
    cout<<"[type] should be cumulative or width.\n";
    cout<<"You requested "<<type<<"\n";
    exit(0);
  };

  exist=true;
  AsmName=inp;
  zdiv=zdivi;
  MaterialID.resize(zdiv);
  zdist.resize(zdiv);
  MaterialName.resize(zdiv);

  real zedge=0.;
  for(int i=0;i<zdiv;i++){
    MaterialName[i]=zmapi[i];
    if(type=="cumulative"){
      zdist[i]=zdisti[i];
    }else{
      zdist[i]=zedge+zdisti[i];
      zedge+=zdisti[i];
    };
  };

  for(int i=1;i<zdiv;i++){
    if(zdist[i-1]>=zdist[i]){
      cout<<"Error in Assembly!\n";
      cout<<i-1<<"-th edge : "<<zdist[i-1]<<"\n";
      cout<<i<<"-th edge : "<<zdist[i]<<"\n";
      exit(0);
    };
  };
};

void Assembly::show_self()
{
  cout<<"*Assembly*\n";
  cout<<" name:"<<AsmName<<"\n";
  cout<<" zdiv= "<<zdiv<<"\n";
  for(int i=0;i<zdiv;i++){
    cout<<"  "<<i<<":  "<< MaterialID[i]<<" ("<<MaterialName[i];
    cout<<") "<<zdist[i]<<"\n";
  };
}

// *********************************************************
//    AssemblySet 
// *********************************************************

AssemblySet::AssemblySet(int matinp)
{
  AsmNum=0;
  if(matinp>MAX_MAT){
    cout<<"You should set bigger MAX_MAT in CartCore.h.\n";
    exit(0);
  };
  MatNum=matinp;
  exist=false;
};

void AssemblySet::show_self()
{
  cout<<"*AssemblySet*";
  cout<<"\n  Number of Material : "<<MatNum<<"\n";
  for(int i=0;i<MatNum;i++){
    cout<<"  "<<i<<":"<<MaterialName[i]<<"\n";
  };

  cout<<"\n  Number of Assembly : "<<AsmNum<<"\n";
  for(int i=0;i<AsmNum;i++){
    cout<<"  "<<i<<":"<<Asm[i].GetName()<<"\n";
  };

  cout<<"\n ** Assembly Information **\n";
  for(int i=0;i<AsmNum;i++){
    Asm[i].show_self();
  };
};

void AssemblySet::PutMaterialName(string *inp)
{
  vector<string> tmp;
  tmp.resize(MatNum);
  for(int i=0;i<MatNum;i++){
    tmp[i]=inp[i];
  };
  PutMaterialName(tmp);
};

void AssemblySet::PutMaterialName(vector<string> inp)
{
  exist=true;
  for(int i=0;i<MatNum;i++){
    MaterialName[i]=inp[i];
  };
};

void AssemblySet::PutAsm(Assembly &inp)
{
  for(int i=0;i<inp.GetZdiv();i++){
    if(inp.GetMaterialID(i)>MatNum)cout<<"There is not Material No. "<<inp.GetMaterialID(i)<<".\n";
  };

  if(AsmNum==MAX_ASM){
    cout<<"You have to increase MAX_ASM in AssemblySet.\n";
    exit(0);
  };

  Asm[AsmNum]=inp;
  for(int i=0;i<Asm[AsmNum].GetZdiv();i++){
    if(Asm[AsmNum].GetMaterialName(i)=="vacuum"||
       Asm[AsmNum].GetMaterialName(i)=="Vacuum"){
      Asm[AsmNum].PutMaterialID(i,-1);
    }else{
      int tmp=-1;
      for(int j=0;j<MatNum;j++){
        if(MaterialName[j]==Asm[AsmNum].GetMaterialName(i))tmp=j;
      };
      if(tmp==-1){
        cout<<"# Error in PutAsm in AssemblySet.\n";
        cout<<"# Material name ["<<Asm[AsmNum].GetMaterialName(i)<<"] does not exist in AssemblySet.\n";
	cout<<"# AssemblySet includes the following materials:\n";
	cout<<"#   ";
	for(int j=0;j<MatNum;j++){
	  cout<<MaterialName[j]<<", ";
	};
	cout<<"\n";
	exit(0);
       };
      Asm[AsmNum].PutMaterialID(i,tmp);
    };
  };

  AsmNum++;
};

Assembly AssemblySet::GetJointAssembly(int n,string *r)
{
  Assembly ret;

  int zdiv=0;
  for(int i=0;i<n;i++){
    int id=-1;
    for(int j=0;j<AsmNum;j++){
      if(Asm[j].GetName()==r[i])id=j;
    };
    if(id==-1){
      cout<<"Assembly "<<r[i]<<" does not exist in this assemblyset.\n";
      exit(0);
    };
    zdiv+=Asm[id].GetZdiv();
  };

  vector<string> zmap(zdiv);
  vector<real> zdist(zdiv);

  int index=0;
  for(int i=0;i<n;i++){
    int id=-1;
    for(int j=0;j<AsmNum;j++){
      if(Asm[j].GetName()==r[i])id=j;
    };
    float zz=0;
    for(int j=0;j<Asm[id].GetZdiv();j++){
      zmap[index]=Asm[id].GetMaterialName(j);
      zdist[index]=Asm[id].GetZdist(j)-zz;
      zz=Asm[id].GetZdist(j);
      index++;
    };
  };

  ret.Init(zdiv,zmap,zdist,"","width");
  return ret;
};

// *********************************************************
//  XYLattice
// *********************************************************

void XYLattice::PutXY(int x,int y)
{
  xm=x;
  ym=y;
  map.resize(x*y,-1);
  xl.resize(x,0.);
  yl.resize(y,0.);
};

void XYLattice::PutMap(int *mapinp)
{
  if(xm==0||ym==0){
    cout<<"Error in XYLattice::PutMap.\n";
    cout<<"X mesh or Y mesh is zero.\n";
    exit(0);
  };

  for(int i=0;i<xm*ym;i++){
    map[i]=mapinp[i];
  };
};

void XYLattice::PutXL(real *inp)
{
  if(xm==0){
    cout<<"Error in XYLattice::PutXL.\n";
    cout<<"X mesh is zero.\n";
    exit(0);
  };

  for(int i=0;i<xm;i++){
    xl[i]=inp[i];
  };
};

void XYLattice::PutYL(real *inp)
{
  if(ym==0){
    cout<<"Error in XYLattice::PutYL.\n";
    cout<<"Y mesh is zero.\n";
    exit(0);
  };

  for(int i=0;i<ym;i++){
    yl[i]=inp[i];
  };
};

void XYLattice::PutMap(int x,int y,int *mapinp)
{
  PutXY(x,y);
  PutMap(mapinp);
};

void XYLattice::PutData(int x,int y,real *xli,real *yli,int *mapinp)
{
  PutXY(x,y);
  PutXL(xli);
  PutYL(yli);
  PutMap(mapinp);
};

void XYLattice::PutSingleData(int mat)
{
  int *tmp=new int[1];
  tmp[0]=mat;
  PutMap(1,1,tmp);
  delete [] tmp;
};

real XYLattice::GetXlSum()
{
  real sum=0.;
  for(int i=0;i<xm;i++){
    sum+=xl[i];
  };
  return sum;
};

real XYLattice::GetYlSum()
{
  real sum=0.;
  for(int i=0;i<ym;i++){
    sum+=yl[i];
  };
  return sum;
};

void XYLatticeSet::PutXYLattice(XYLattice &latinp)
{
  latnum++;
  xylat.push_back(latinp);
};

// *********************************************************
//  CartCore 
// *********************************************************

CartCore::~CartCore()
{
};

void CartCore::Initialize(int x,int y)
{
  xr=x;
  yr=y;
  AsmMap.resize(xr*yr);
  xwid.resize(xr);
  ywid.resize(yr);
}

void CartCore::PutAsmMap(int *inp)
{
  vector<int> inp2(yr*xr);
  for(int i=0;i<xr*yr;i++){
    inp2[i]=inp[i];
  };
  PutAsmMap(inp2);
}

void CartCore::PutAsmMap(vector<int> inp)
{
  for(int i=0;i<xr*yr;i++){
    AsmMap[i]=inp[i];
  };
};

void CartCore::PutAsmMap(int x1,int x2,vector<int> inp)
{
  int index=0;
  int index2=0;
  for(int j=0;j<yr;j++){
    for(int k=0;k<xr;k++){
      if(k>=x1&&k<=x2){
        AsmMap[index]=inp[index2];
	index2++;
      };
      index++;
    };
  };
}

void CartCore::PutAsmMap(int x1,int x2,int *inp)
{
  int tmp=yr*(x2-x1+1);
  vector<int> inp2(tmp);
  int index2=0;
  for(int j=0;j<yr;j++){
    for(int k=0;k<xr;k++){
      if(k>=x1&&k<=x2){
        inp2[index2]=inp[index2];
	index2++;
      };
    };
  };
  PutAsmMap(x1,x2,inp2);
}

void CartCore::PutAsmMap_Hex_to_XY(int *inp)
{
  int xrold=xr;

  Initialize(xr*2+1,yr);

  cout<<"+++ Assembly map in hexagonal description is transformed to XY system.+++\n\n";
  cout<<"   New map size   x : "<<xr<<"\n";
  cout<<"                  y : "<<yr<<"\n\n";

  int *AsmInp2=new int[xr*yr];
  int index=0;
  for(int i=0;i<yr;i++){
    if(i%2==0)AsmInp2[index++]=-1;
    for(int j=0;j<xrold;j++){
      AsmInp2[index++]=inp[i*xrold+j];
      AsmInp2[index++]=inp[i*xrold+j];
    };
    if(i%2==1)AsmInp2[index++]=-1;
  };
  PutAsmMap(AsmInp2);

  ShowAssemblyMapNew();

  delete [] AsmInp2;
};

void CartCore::ShowAssemblyMapNew()
{
  int index=0;
  string cc[]={
    "A","B","C","D","E", "F","G","H","I","J",
    "K","L","M","N","O", "P","Q","R","S","T",
    "U","V","W","X","Y"
  };
  for(int i=0;i<yr;i++){
    for(int j=0;j<xr;j++){
      int tmp=AsmMap[index];
      if(tmp>=0&&tmp<=9){
	cout<<tmp;
      }else if(tmp==-1){
	cout<<"-";
      }else{
	int tt=tmp-10;
	tt=tt%25;
	cout<<cc[tt];
      };
      index++;
    };
    cout<<"\n";
  };
  cout<<"\n";
};

void CartCore::ShowAssemblyMap(int ii,int ij)
{
  if(ij==-1)ij=ii;
  int index=0;
  for(int j=0;j<yr;j++){
    if(j%2==0)cout<<" ";
    for(int k=0;k<xr;k++){
      if(AsmMap[index]>=ii&&AsmMap[index]<=ij){cout<<"X ";}
      else{cout<<"  ";};
      index++;
    };
    cout<<"\n";
  };
};

void CartCore::CountAssembly(int ii)
{
  int sum=0;
  int index=0;
  for(int i=0;i<yr;i++){
    for(int j=0;j<xr;j++){
      if(AsmMap[index]==ii)sum++;
      index++;
    };
  };
  cout<<sum<<"\n";
};

void CartCore::ShowAssemblyMap()
{
  int index=0;
  for(int j=0;j<yr;j++){
    for(int k=0;k<xr;k++){
      cout.width(3);
      cout<<AsmMap[index];
      index++;
    };
    cout<<"\n";
  };
};

void CartCore::GetMaterialMapXY(real z,int *map)
{
  int zr=GetUnifiedZMesh();
  int *MaterialMap=new int[xr*yr*zr];
  real *zwid=new real[zr];
  MakeMaterialMap(zr,MaterialMap,zwid);

  int index=-1;
  real za=0.;
  for(int i=0;i<zr;i++){
    za+=zwid[i];
    if(za>z){
      index=i;
      break;
    };
  }

  for(int i=0;i<yr;i++){
    for(int j=0;j<xr;j++){
      map[i*xr+j]=MaterialMap[index*(xr*yr)+i*xr+j];
    };
  };

  delete [] MaterialMap;
  delete [] zwid;
}

void CartCore::ShowMaterialMap()
{
  int zr=GetUnifiedZMesh();
  int *MaterialMap=new int[xr*yr*zr];
  real *zwid=new real[zr];

  MakeMaterialMap(zr,MaterialMap,zwid);

  int index=0;
  for(int i=0;i<zr;i++){
    cout<<"plane:"<<i<<"  width="<<zwid[i]<<"\n";
    for(int j=0;j<yr;j++){
      for(int k=0;k<xr;k++){
	cout.width(3);
	cout<<MaterialMap[index];
	index++;
      };
      cout<<"\n";
    };
  };

  delete [] MaterialMap;
  delete [] zwid;
};

int CartCore::GetUnifiedZMesh()
{
  int zdiv=0;
  real *zdist=new real[999];
  int AsmType=AsmSet->GetAsmNum();

  for(int i=0;i<AsmType;i++){
    //    if(UsedOrNotAssembly(i)==1){
      for(int j=0;j<AsmSet->GetAsm(i).GetZdiv();j++){
        real obj=AsmSet->GetAsm(i).GetZdist(j);
        int judge=0;
        for(int k=0;k<zdiv;k++){
	  if(fabs(obj-zdist[k])<0.001)judge=1;
        };
        if(judge==0){
	  zdist[zdiv]=obj;
	  zdiv++;
        };
      };
      //    };
  };

  int *tmp=new int[zdiv];
  ChangeOrder(zdist,zdiv,tmp);
  delete [] zdist;
  delete [] tmp;

  return zdiv;
};

void CartCore::PutBC(int *inp)
{
  LeftBC  =inp[0];
  RightBC =inp[1];
  BackBC  =inp[2];
  FrontBC =inp[3];
  UpperBC =inp[4];
  BottomBC=inp[5];
}

void CartCore::PutBC(vector<int> inp)
{
  LeftBC  =inp[0];
  RightBC =inp[1];
  BackBC  =inp[2];
  FrontBC =inp[3];
  UpperBC =inp[4];
  BottomBC=inp[5];
};

void CartCore::PutBoundaryCondition(string xl,string xr,string yl,string yr,string zl,string zr)
{
  string bcname[]={"Zeroflux","Reflective","Vacuum"};
  for(int i=0;i<3;i++){
    if(xl==bcname[i])LeftBC=i;
    if(xr==bcname[i])RightBC=i;
    if(yl==bcname[i])BackBC=i;
    if(yr==bcname[i])FrontBC=i;
    if(zl==bcname[i])UpperBC=i;
    if(zr==bcname[i])BottomBC=i;
  };

  int errcode=-1;
  if(LeftBC<0||LeftBC>2)    errcode=0;
  if(RightBC<0||RightBC>2)  errcode=1;
  if(BackBC<0||BackBC>2)    errcode=2;
  if(FrontBC<0||FrontBC>2)  errcode=3;
  if(UpperBC<0||UpperBC>2)  errcode=4;
  if(BottomBC<0||BottomBC>2)errcode=5;
  if(errcode!=-1){
    cout<<"Error in Boundary condition in CartCore.\n";
    cout<<"Boundary No."<<errcode<<"\n";
    cout<<" (0:x-,1:x+,2:y-,3:y+,4:z-,5:z+)\n";
    exit(0);
  };
};

void CartCore::PutWidthXY(real *xinp,real *yinp)
{
  vector<real> xinp2(xr);
  vector<real> yinp2(yr);
  for(int i=0;i<xr;i++){
    xinp2[i]=xinp[i];
  };
  for(int i=0;i<yr;i++){
    yinp2[i]=yinp[i];
  };
  PutWidthXY(xinp2,yinp2);
}

void CartCore::PutWidthXY(vector<real> xinp, vector<real> yinp)
{
  for(int i=0;i<xr;i++){
    xwid[i]=xinp[i];
    if(xwid[i]<=0.){
      cout<<"Error in PutWidthXY.\n";
      cout<<"X-width is zero or negative.\n";
      exit(0);
    };
  };
  for(int i=0;i<yr;i++){
    ywid[i]=yinp[i];
    if(ywid[i]<=0.){
      cout<<"Error in PutWidthXY.\n";
      cout<<"Y-width is zero or negative.\n";
      exit(0);
    };
  };
};

void CartCore::PutWidthFromPitch(real pitch)
{
  real vol=0.5*sqrt(3.)*pitch*pitch;
  real xlen=pitch*0.5;
  real ylen=vol/(xlen*2.);
  vector<real> xinp(xr,xlen);
  vector<real> yinp(yr,ylen);
  PutWidthXY(xinp,yinp);
};

void CartCore::MakeMaterialMap(int &zr,int *MaterialMap,real *zwid)
{
  for(int y=0;y<yr;y++){
    for(int x=0;x<xr;x++){
      if(!AsmSet->GetAsm(AsmMap[y*xr+x]).Exist()){
	cout<<"Error in MakeMaterialMap of CartCore class.\n";
	cout<<"Assembly "<<AsmMap[y*xr+x]<<" does not exist.\n";
	cout<<"  x="<<x<<", y="<<y<<"\n";
	exit(0);
      };
    };
  };

  int zdiv=0;
  real *zdist=new real[999];

  int AsmType=AsmSet->GetAsmNum();

  for(int i=0;i<AsmType;i++){
    //    if(UsedOrNotAssembly(i)==1){
      for(int j=0;j<AsmSet->GetAsm(i).GetZdiv();j++){
        real obj=AsmSet->GetAsm(i).GetZdist(j);
        int judge=0;
        for(int k=0;k<zdiv;k++){
  	  if(fabs(obj-zdist[k])<0.001)judge=1;
        };
        if(judge==0){
	  zdist[zdiv]=obj;
	  zdiv++;
        };
      };
      //    };
  };

  int *tmp=new int[zdiv];
  ChangeOrder(zdist,zdiv,tmp);

  for(int i=0;i<zdiv;i++){
    zwid[i]=zdist[tmp[i]];
  };

  zr=zdiv;

  int tmp3=0;
  int index=0;
  for(int i=0;i<zdiv;i++){
    for(int j=0;j<yr;j++){
      for(int k=0;k<xr;k++){
	if(AsmMap[j*xr+k]==-1){
	  MaterialMap[index]=-1;
	}else{
  	  for(int l=0;l<AsmSet->GetAsm(AsmMap[j*xr+k]).GetZdiv();l++){
	    real bt=AsmSet->GetAsm(AsmMap[j*xr+k]).GetZdist(l);
	    if(zwid[i]-bt<0.001){
	      tmp3=AsmSet->GetAsm(AsmMap[j*xr+k]).GetMaterialID(l);
	      break;
	    };
	  };
	  MaterialMap[index]=tmp3;
	};
	index++;
      };
    };
  };

  for(int i=zr-1;i>0;i--){
    zwid[i]=zwid[i]-zwid[i-1];
  };

  delete [] zdist;
  delete [] tmp;
};

int CartCore::UsedOrNotAssembly(int inp)
{
  for(int i=0;i<xr*yr;i++){
    if(AsmMap[i]==inp)return 1;
  };
  return 0;
};

void CartCore::ChangeAssembly(int x,int y,int asmid)
{
  AsmMap[y*xr+x]=asmid;
};

void CartCore::ChangeAssembly(int p,int* xp,int* yp,int asmid)
{
  for(int i=0;i<p;i++){
    ChangeAssembly(xp[i],yp[i],asmid);
  };
};

void CartCore::ChangeAssembly(int asmid1,int asmid2)
{
  for(int i=0;i<yr;i++){
    for(int j=0;j<xr;j++){
      if(AsmMap[i*xr+j]==asmid1)AsmMap[i*xr+j]=asmid2;
    };
  };
};

void CartCore::PutLatticeMap(int xi,int yi,int *latmap,real *xwid,real *ywid,XYLatticeSet &latset,bool print)
{
  int latnum=latset.GetLatnum();
  int retx=0;
  int rety=0;
  vector<real> retxl;
  vector<real> retyl;
  vector<int> divx;
  vector<int> divy;

  // X-direction sort
  for(int x=0;x<xi;x++){

    int ll=0;
    vector<int> llid;
    for(int y=0;y<yi;y++){
      int latno=latmap[y*xi+x];
      if(latno>=latnum){
	cout<<"# Error in CartCore::PutLatticeMap.\n";
	cout<<"#  Lattice ID         : "<<latno<<"\n";
	cout<<"#  Number of lattices : "<<latnum<<"\n";
      };
      if(ll==0){
	ll++;
	llid.push_back(latno);
      }else{
	bool same=false;
	for(int i=0;i<ll;i++){
	  if(latno==llid[i])same=true;
	};
	if(!same){
	  ll++;
	  llid.push_back(latno);
	};
      };
    };

    int div=0;
    vector<real> dl;
    real edge=xwid[x];
    for(int i=0;i<ll;i++){
      int id=llid[i];
      real x0=0.;
      // Total length check
      if(latset.GetXYLattice(id).GetXm()!=1){
	if(!SameReal(latset.GetXYLattice(id).GetXlSum(),edge,1e-6)){
	  cout<<"Error in CartCore::PutLatticeMap.\n";
	  cout<<"Lattice size is inconsistent in X.\n";
	  cout<<"   "<<latset.GetXYLattice(id).GetXlSum()<<"\n";
	  cout<<"   "<<edge<<"\n";
	  exit(0);
	};
      };
      for(int j=0;j<latset.GetXYLattice(id).GetXm()-1;j++){
	x0+=latset.GetXYLattice(id).GetXl(j);
	if(div==0){
	  div++;
	  dl.push_back(x0);
	}else{
	  bool same=false;
	  for(int k=0;k<div;k++){
	    if(SameReal(dl[k],x0,1e-5))same=true;
	  };
	  if(!same){
	    div++;
	    dl.push_back(x0);
	  };
	};
      };
    };

    real *dl2=new real[div+1];
    divx.push_back(div+1);
    for(int i=0;i<div;i++){
      dl2[i]=dl[i];
    };
    dl2[div]=edge;
    ChangeOrder(dl2,div);

    retx+=div+1;
    for(int i=0;i<div+1;i++){
      real tmp=dl2[i];
      if(i!=0)tmp-=dl2[i-1];
      retxl.push_back(tmp);
    };
  };
  // X-direction end

  // Y-direction sort
  for(int y=0;y<yi;y++){

    int ll=0;
    vector<int> llid;
    for(int x=0;x<xi;x++){
      int latno=latmap[y*xi+x];
      if(ll==0){
	ll++;
	llid.push_back(latno);
      }else{
	bool same=false;
	for(int i=0;i<ll;i++){
	  if(latno==llid[i])same=true;
	};
	if(!same){
	  ll++;
	  llid.push_back(latno);
	};
      };
    };

    int div=0;
    vector<real> dl;
    real edge=ywid[y];
    for(int i=0;i<ll;i++){
      int id=llid[i];
      real y0=0.;
      // Total length check
      if(latset.GetXYLattice(id).GetYm()!=1){
	if(!SameReal(latset.GetXYLattice(id).GetYlSum(),edge,1e-6)){
	  cout<<"Error in CartCore::PutLatticeMap.\n";
	  cout<<"Lattice size is inconsistent in Y.\n";
	  cout<<"   "<<latset.GetXYLattice(id).GetYlSum()<<"\n";
	  cout<<"   "<<edge<<"\n";
	  exit(0);
	};
      };
      //
      for(int j=0;j<latset.GetXYLattice(id).GetYm()-1;j++){
	y0+=latset.GetXYLattice(id).GetYl(j);
	if(div==0){
	  div++;
	  dl.push_back(y0);
	}else{
	  bool same=false;
	  for(int k=0;k<div;k++){
	    if(SameReal(dl[k],y0,1e-5))same=true;
	  };
	  if(!same){
	    div++;
	    dl.push_back(y0);
	  };
	};
      };
    };

    real *dl2=new real[div+1];
    divy.push_back(div+1);
    for(int i=0;i<div;i++){
      dl2[i]=dl[i];
    };
    dl2[div]=edge;
    ChangeOrder(dl2,div);

    rety+=div+1;
    for(int i=0;i<div+1;i++){
      real tmp=dl2[i];
      if(i!=0)tmp-=dl2[i-1];
      retyl.push_back(tmp);
    };
  };
  // Y-direction end

  // Create map
  Initialize(retx,rety);
  PutWidthXY(retxl,retyl);

  int index=0;
  int indy=0;
  for(int i=0;i<yi;i++){
    real y0=0.;
    for(int i2=0;i2<divy[i];i2++){
      y0+=retyl[indy++];
      int indx=0;
        for(int j=0;j<xi;j++){ 
        int latno=latmap[i*xi+j];
	XYLattice lat=latset.GetXYLattice(latno);
        real x0=0.;
	for(int j2=0;j2<divx[j];j2++){
          x0+=retxl[indx++];

	  int xpos=-1;
	  real x00=0.;
	  for(int k=0;k<lat.GetXm();k++){
	    if(lat.GetXm()==1){
	      x00+=xwid[j];
	    }else{
  	      x00+=lat.GetXl(k);
	    };
	    if((x00>=x0||SameReal(x00,x0,1e-5))&&xpos==-1){
	      xpos=k;
	    };
	  };
          int ypos=-1;
          real y00=0.;
	  for(int k=0;k<lat.GetYm();k++){
	    if(lat.GetYm()==1){
	      y00+=ywid[i];
	    }else{
  	      y00+=lat.GetYl(k);
	    };
	    if((y00>=y0||SameReal(y00,y0,1e-5))&&ypos==-1){
	      ypos=k;
	    };
	  };
	  AsmMap[index++]=lat.GetMap(xpos,ypos);

	};
      };
    };
  };

  if(print){
  cout<<"** CartCore::PutLatticeMap **\n\n";
  cout<<"[X-direction]\n";
  int tmp=0;
  for(int i=0;i<xi;i++){
    cout<<i<<" : "<<divx[i]<<"\n";
    for(int j=0;j<divx[i];j++){
      cout<<"     "<<tmp<<" "<<retxl[tmp]<<"\n";    
      tmp++;
    };
  };
  cout<<"\n[Y-direction]\n";
  tmp=0;
  for(int i=0;i<yi;i++){
    cout<<i<<" : "<<divy[i]<<"\n";
    for(int j=0;j<divy[i];j++){
      cout<<"     "<<tmp<<" "<<retyl[tmp]<<"\n";    
      tmp++;
    };
  };
  };

};
