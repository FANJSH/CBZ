#include "ENDFIDTranslator.h"

ENDFIDTranslator::ENDFIDTranslator()
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
  string nucln_srac_inp[]={
    "H-","HE","LI","BE","B-", "C-","N-","O-","F-","NE",
    "NA","MG","AL","SI","P-", "S-","CL","AR","K-","CA",
    "SC","TI","V-","CR","MN", "FE","CO","NI","CU","ZN",
    "GA","GE","AS","SE","BR", "KR","RB","SR","Y-","ZR",
    "NB","MO","TC","RU","RH", "PD","AG","CD","IN","SN",
    "SB","TE","I-","XE","CS", "BA","LA","CE","PR","ND",
    "PM","SM","EU","GD","TB", "DY","HO","ER","TM","YB",
    "LU","HF","TA","W-","RE", "OS","IR","PT","AU","HG",
    "TL","PB","BI","PO","AT", "RN","FR","RA","AC","TH",
    "PA","U-","NP","PU","AM", "CM","BK","CF","ES","FM",
    "Md","No","Lr","Rf","Db", "Sg","Bh","Hs","Mt","Ds",
    "Rg","Cn"
  };
  int mass_b_inp[]={
    1,3,6,9,10, 12,14,16,19,20,
    23,24,27,28,31, 32,35,40,39,40,
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
  nucln_srac.resize(nuclnum);
  mass_b.resize(nuclnum);
  for(int i=0;i<nuclnum;i++){
    nucln[i]=nucln_inp[i];
    nucln_srac[i]=nucln_srac_inp[i];
    mass_b[i]=mass_b_inp[i];
  };
};

void ENDFIDTranslator::GetParameter(string name,int &iz,int &ia,int &il)
{
  iz=AtomicNumberFromName(name);
  int sz_isoname=nucln[iz-1].size();
  int sz_all=name.size();
  ia=StringToInt(name.substr(sz_isoname,3));
  il=0;
  if(sz_all!=sz_isoname+3){
    if(name.substr(sz_all-1,1)=="m"||
       name.substr(sz_all-1,1)=="M")il=1;
    if(name.substr(sz_all-1,1)=="n"||
       name.substr(sz_all-1,1)=="N")il=2;
  };
  
  /*
  int id=ID(name);
  
  int i1=id%100;
  il=(i1-1)%3;

  iz=(id-i1)/100;

  int size=name.size();
  int lpos=size-3;
  if(name.substr(size-1,1)=="m"||name.substr(size-1,1)=="n")lpos--;
  ia=StringToInt(name.substr(lpos,3));
  */
};

void ENDFIDTranslator::GetParameterNew(int id,int &iz,int &ia,int &il)
{
  il=id%10;
  int tmp=(id-il)/10;
  ia=tmp%1000;
  iz=(tmp-ia)/1000;
};

string ENDFIDTranslator::NewName(int id)
{
  int iz,ia,il;
  GetParameterNew(id,iz,ia,il);
  return Name(iz,ia,il);
};

string ENDFIDTranslator::SearchNameFromParameter(int iz,int ia,int il,bool srac)
{
  string ret=nucln[iz-1];
  if(srac)ret=nucln_srac[iz-1];
  
  string mass_st=IntToString(ia);
  if(mass_st.size()==1){
    ret=ret+"00";
  }else if(mass_st.size()==2){
    ret=ret+"0";
  };
  ret=ret+mass_st;
  if(il==1){
    if(srac){
      ret=ret+"M";
    }else{
      ret=ret+"m";
    };
  };
  if(il==2){
    if(srac){
      ret=ret+"N";
    }else{
      ret=ret+"n";
    };
  };

  return ret;
};

string ENDFIDTranslator::NameFromAtomicNumber(int an)
{
  if(an==0)return "n";
  if(an-1<nucln.size())return nucln[an-1];
  cout<<"# Error in ENDFIDTranslator::NameFromAtomicNumber.\n";
  cout<<"# The requested atomic number : "<<an<<"\n";
  exit(0);
};

int ENDFIDTranslator::AtomicNumberFromName(string name)
{
  string isoname=ExtractIsotopeName(name);
  for(int i=0;i<nuclnum;i++){
    if(isoname==nucln[i])return i+1;
  };
  cout<<"# Error in ENDFIDTranslator::AtomicNumberFromName.\n";
  cout<<"# Nuclide "<<name<<" cannot be found.\n";
  exit(0);
};

string ENDFIDTranslator::ExtractIsotopeName(string name)
{
  int sz=name.size();

  int meta=0;
  if(name.substr(sz-1,1)=="m"||
     name.substr(sz-1,1)=="M")meta=1;
  if(name.substr(sz-1,1)=="n"||
     name.substr(sz-1,1)=="N")meta=2;

  int name_len=sz-3;
  if(meta!=0)name_len--;
  string nucname=name.substr(0,name_len);
  return nucname;
};

int ENDFIDTranslator::NewID(int ia,int iz,int il)
{
  return ia*10000+iz*10+il;
};

/*
int ENDFIDTranslator::NewID(string name)
{
  int ia,iz,il;
  GetParameterNew(NewID(
};
*/

int ENDFIDTranslator::GetENDFID(string name)
{
  if(name=="C000")return 600;
  
  int sz=name.size();

  int meta=0;
  if(name.substr(sz-1,1)=="m"||
     name.substr(sz-1,1)=="M")meta=1;
  if(name.substr(sz-1,1)=="n"||
     name.substr(sz-1,1)=="N")meta=2;

  int name_len=sz-3;
  if(meta!=0)name_len--;
  string nucname=name.substr(0,name_len);

  int mass_pos=1;
  if(name_len!=1)mass_pos=2;  
  int mass=StringToInt(name.substr(mass_pos,3));

  int atm=-1;
  for(int i=0;i<nuclnum;i++){
    if(nucname==nucln[i])atm=i+1;
  };
  if(atm==-1){
    cout<<"# Error in ENDFIDTranslator::GetENDFID.\n";
    cout<<"# Nuclide "<<nucname<<" cannot be found.\n";
    exit(0);
  };

  int d_mass=mass-mass_b[atm-1];

  int matno=atm*100+25;
  matno+=d_mass*3+meta;

  int tmp=matno%100;
  tmp=(matno-tmp)/100;
  if(tmp!=atm){
    cout<<"\n# Error in ENDFIDTranslator::GetENDFID\n";
    cout<<"# ENDF-ID cannot be determined.\n";
    cout<<"# Nuclide name is "<<name<<"\n";
    cout<<"# Atomic number is "<<AtomicNumberFromName(nucname)<<"\n";
    exit(0);
  };

  if(mass==0)matno=atm*100;

  return matno;
};

string ENDFIDTranslator::GetNuclideName(int mat)
{
  if(mat>9990){
    if(mat==9991)return "FP.P1";
    if(mat==9995)return "FP.U5";
    if(mat==9998)return "FP.U8";
    if(mat==9999)return "FP.P9";
  };

  int no=mat%100;
  int atm=(mat-no)/100;
  int mass,meta;

  if(no>=25){
    int dif=no-25;
    meta=dif%3;
    dif-=meta;
    mass=mass_b[atm-1]+dif/3;
  }else{
    int dif=25-no;
    meta=dif%3;
    if(meta!=0)dif+=(3-meta);
    meta=3-meta;
    mass=mass_b[atm-1]-dif/3;
  };

  string ret=nucln[atm-1];

  if(no==0){
    ret=ret+"000";
    return ret;
  };

  string mass_st=IntToString(mass);
  if(mass_st.size()==1){
    ret=ret+"00";
  }else if(mass_st.size()==2){
    ret=ret+"0";
  };
  ret=ret+mass_st;
  if(meta==1)ret=ret+"m";
  if(meta==2)ret=ret+"n";
  return ret;
};

void ENDFIDTranslator::check()
{
};

