#include "MATIDTranslator.h"

MATIDTranslator::MATIDTranslator()
{
  nuclnum=112;
  string nucln_inp[]={
    "H","He","Li","Be","B", "C","N","O","F","Ne",
    "Na","Mg","Al","Si","P", "S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn", "Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br", "Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh", "Pd","Ag","Cd","In","Sn",
    "Sb","Te","I","Xe","Cs", "Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb", "Dy","Ho","Er","Tm","Yb",
    "Lu","Hf","Ta","W","Re", "Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At", "Rn","Fr","Ra","Ac","Th",
    "Pa","U","Np","Pu","Am", "Cm","Bk","Cf","Es","Fm",
    "Md","No","Lr","Rf","Db", "Sg","Bh","Hs","Mt","Ds",
    "Rg","Cn"
  };
  int mass_b_inp[]={
    1,3,6,9,10, 12,14,16,19,20,
    23,24,27,28,31, 32,35,36,39,40,
    45,46,50,50,55, 54,59,58,63,64,
    69,70,75,74,79, 78,85,84,89,90,
    93,92,97,96,103, 102,107,106,113,112,
    121,120,127,124,133, 130,138,136,141,142,
    139,144,151,152,159, 156,165,162,169,168,
    175,174,180,180,185, 184,191,190,197,196,
    203,204,209,206,203, 211,212,223,225,227,
    229,234,230,235,235, 240,240,240,240,240,
    240,240,240,240,240, 240,240,240,240,240,
    240,240,240,
  };  
  nucln.resize(nuclnum);
  mass_b.resize(nuclnum);
  for(int i=0;i<nuclnum;i++){
    nucln[i]=nucln_inp[i];
    mass_b[i]=mass_b_inp[i];
  };
};

string MATIDTranslator::ExtractIsotopeName(string name)
{
  string tmp=name.substr(1,1);
  if(tmp=="0"||tmp=="1"||tmp=="2"){
    return name.substr(0,1);
  }else{
    return name.substr(0,2);
  };
};

string MATIDTranslator::NameFromAtomicNumber(int an)
{
  if(an==0)return "n";
  if(an-1<nucln.size())return nucln[an-1];
  cout<<"# Error in MATIDTranslator::NameFromAtomicNumber.\n";
  cout<<"# The requested atomic number : "<<an<<"\n";
  exit(0);
};

int MATIDTranslator::AtomicNumberFromName(string name)
{
  string isoname=ExtractIsotopeName(name);
  if(isoname=="n"){
    //cout<<name<<"\n";
    return 0;
  };
  for(int i=0;i<nuclnum;i++){
    if(isoname==nucln[i])return i+1;
  };
  cout<<"# Error in MATIDTranslator::AtomicNumberFromName.\n";
  cout<<"# Nuclide "<<name<<" cannot be found.\n";
  exit(0);
};

void MATIDTranslator::GetParameter(string name,int &iz,int &ia,int &il)
{
  iz=AtomicNumberFromName(name);
  int sz_isoname=1;
  if(iz!=0)sz_isoname=nucln[iz-1].size();
  int sz_all=name.size();
  ia=StringToInt(name.substr(sz_isoname,3));
  il=0;
  if(sz_all!=sz_isoname+3){
    if(name.substr(sz_all-1,1)=="m"||
       name.substr(sz_all-1,1)=="M")il=1;
    if(name.substr(sz_all-1,1)=="n"||
       name.substr(sz_all-1,1)=="N")il=2;
    if(name.substr(sz_all-1,1)=="o"||
       name.substr(sz_all-1,1)=="O")il=3;
  };
};

void MATIDTranslator::GetParameter(int id,int &iz,int &ia,int &il)
{
  il=id%10;
  int tmp=(id-il)/10;
  ia=tmp%1000;
  iz=(tmp-ia)/1000;
};

int MATIDTranslator::GetMATID(int iz,int ia,int il)
{
  return iz*10000+ia*10+il;
};

int MATIDTranslator::GetMATID(string name)
{
  int iz,ia,il;
  GetParameter(name,iz,ia,il);
  return GetMATID(iz,ia,il);
};

string MATIDTranslator::GetName(int iz,int ia,int il)
{
  string ret=NameFromAtomicNumber(iz);
  string mass_st=IntToString(ia);
  if(mass_st.size()==1){
    ret=ret+"00";
  }else if(mass_st.size()==2){
    ret=ret+"0";
  };
  ret=ret+mass_st;
  if(il==1)ret=ret+"m";
  if(il==2)ret=ret+"n";
  if(il==3)ret=ret+"o";
  return ret;
};

string MATIDTranslator::GetName(int id)
{
  if(id>9900000){
    if(id==9910000){
      return "FP.Pu241";
    }else if(id==9950000){
      return "FP.U235";
    }else if(id==9980000){
      return "FP.U238";
    }else if(id==9990000){
      return "FP.Pu239";
    };
    //cout<<"# Error in MATIDTranslator::GetName\n";
    //exit(0);
    return "NO_NAME";
  };
  int iz,ia,il;
  GetParameter(id,iz,ia,il);
  return GetName(iz,ia,il); 
};

int MATIDTranslator::GetMATIDFromENDFID(int mat)
{
  if(mat>9999)return mat;
  if(mat>=9990)return 9900000+(mat%10)*10000; // (pseudo FP)

  int no=mat%100;
  int atm=(mat-no)/100;
  int mass=0;
  int meta=0;

  if(no>=25){
    int dif=no-25;
    meta=dif%3;
    dif-=meta;
    mass=mass_b[atm-1]+dif/3;
  }else if(no!=0){
    int dif=25-no;
    meta=dif%3;
    if(meta!=0){
      dif+=(3-meta);
      meta=3-meta;
    };
    mass=mass_b[atm-1]-dif/3;
  };

  return GetMATID(atm,mass,meta);
};


