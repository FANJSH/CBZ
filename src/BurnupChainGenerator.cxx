#include "BurnupChainGenerator.h"

BCGNuclide::BCGNuclide(int at,int ms,int lv)
{
  atomic_number=at;
  mass_number=ms;
  ex_level=lv;
  flag=false;
  channel=0;
  decay_energy.resize(3);
  delta_decay_energy.resize(3);
  //
  r_num=2;
  r_channel.resize(r_num);
  r_atn_next.resize(r_num);
  r_mas_next.resize(r_num);
  r_lev_next.resize(r_num);
  r_br.resize(r_num);
  PutReactionChannel(0,1); // (n,g)
  PutReactionChannel(1,1); // (n,2n)
  // (default setting)
  r_atn_next[0][0]=atomic_number;
  r_mas_next[0][0]=mass_number+1;
  r_lev_next[0][0]=0;
  r_br[0][0]=1.;
  r_atn_next[1][0]=atomic_number;
  r_mas_next[1][0]=mass_number-1;
  r_lev_next[1][0]=0;
  r_br[1][0]=1.;
};

// (reaction data)

void BCGNuclide::PutReactionChannel(int i,int j)
{
  if(i<0||i>=r_num){
    cout<<"# Error in BCGNuclide::PutReactionChannel.\n";
    cout<<"# Channel ID "<<i<<" is not permitted.\n";
    exit(0);
  };
  r_channel[i]=j;
  r_atn_next[i].resize(j);
  r_mas_next[i].resize(j);
  r_lev_next[i].resize(j);
  r_br[i].resize(j);
};

void BCGNuclide::PutReactionData(int ii,int ic,real *br)
{
  PutReactionChannel(ii,ic);
  int dm=1;
  if(ii==1)dm=-1;
  for(int i=0;i<ic;i++){
    r_atn_next[ii][i]=atomic_number;
    r_mas_next[ii][i]=mass_number+dm;
    r_lev_next[ii][i]=i;
    r_br[ii][i]=br[i];
  };
};

void BCGNuclide::PutReactionData(int ii,real br1,real br2)
{
  real *br=new real[2];
  br[0]=br1;
  br[1]=br2;
  PutReactionData(ii,2,br);
  delete [] br;
};

void BCGNuclide::PutReactionData(int ii,real br1,real br2,real br3)
{
  real *br=new real[3];
  br[0]=br1;
  br[1]=br2;
  br[2]=br3;
  PutReactionData(ii,3,br);
  delete [] br;
};

void BCGNuclide::PutReactionData(int ii,int ic,vector<int>i1,vector<int>i2,vector<int>i3,vector<real>i4)
{
  PutReactionChannel(ii,ic);
  for(int i=0;i<ic;i++){
    r_atn_next[ii][i]=i1[i];
    r_mas_next[ii][i]=i2[i];
    r_lev_next[ii][i]=i3[i];
    r_br[ii][i]=i4[i];
  };
};

// (decay data)

void BCGNuclide::PutChannel(int i)
{
  channel=i;
  atn_next.resize(channel);
  mas_next.resize(channel);
  lev_next.resize(channel);
  decay_type.resize(channel);
  br.resize(channel);
  delta_br.resize(channel);
  dn.resize(channel);
};

void BCGNuclide::PutNextNuclideData(int i,int at,int ms,int lv,real b,real db,int type,int dnin)
{
  atn_next[i]=at;
  mas_next[i]=ms;
  lev_next[i]=lv;
  decay_type[i]=type;
  br[i]=b;
  delta_br[i]=db;
  dn[i]=dnin;
};

void BCGNuclide::PutDecayData(int ii,vector<int>i1,vector<int>i2,vector<int>i3,vector<real>i4,vector<real>i5,vector<real>i6)
{
  PutChannel(ii);
  for(int i=0;i<channel;i++){
    atn_next[i]=i1[i];
    mas_next[i]=i2[i];
    lev_next[i]=i3[i];
    br[i]=i4[i];
    delta_br[i]=i5[i];
    dn[i]=i6[i];
  };
};

void BCGNuclide::NormalizeBr()
{
  real sum=0.;
  for(int i=0;i<channel;i++){
    sum+=br[i];
  };
  sum=1./sum;
  for(int i=0;i<channel;i++){
    br[i]*=sum;
  };
};

bool BCGNuclide::SameNuclide(int atm,int mas,int lev)
{
  if(atm!=atomic_number)return false;
  if(mas!=mass_number)return false;
  if(lev!=ex_level)return false;
  return true;
};

int BCGNuclide::GetYieldID(string tagname)
{
  int sz=yield.size();
  for(int i=0;i<sz;i++){
    if(yield_tag[i]==tagname){
      return i;
    };
  };
  return -1;
};

void BCGNuclide::PutYield(string tagname,real inp,real inp2)
{
  int ii=GetYieldID(tagname);
  if(ii!=-1){
    yield[ii]=inp;
    delta_yield[ii]=inp2;
  }else{
    yield.push_back(inp);
    delta_yield.push_back(inp2);
    yield_tag.push_back(tagname);
  };
};

void BCGNuclide::AddYield(string tagname,real inp,real inp2)
{
  int ii=GetYieldID(tagname);
  if(ii!=-1){
    yield[ii]+=inp;
    real org=delta_yield[ii];
    delta_yield[ii]=sqrt(org*org+inp2*inp2);  // (no-correlation is assumed)
  }else{
    PutYield(tagname,inp,inp2); 
  };
};

real BCGNuclide::GetYield(string tagname)
{
  int ii=GetYieldID(tagname);
  if(ii==-1){
    return 0.;
  }else{
    return yield[ii];
  };
};

real BCGNuclide::GetDeltaYield(string tagname)
{
  int ii=GetYieldID(tagname);
  if(ii==-1){
    return 0.;
  }else{
    return delta_yield[ii];
  };
};

bool BCGNuclide::ZeroYield()
{
  int sz=yield.size();
  for(int i=0;i<sz;i++){
    if(yield[i]!=0.)return false;
  };
  return true;
};

void BCGNuclide::PutTotalNumberEmittedNeutronsCh(int ii, vector<real> &n_emit_ch)
{
  n_per_nuc_ch.resize(ii);
  for(int i=0;i<ii;i++){
    n_per_nuc_ch[i]=n_emit_ch[i];
  };
};


// -----------------

int BCGManager::GetNuclideIndex(int atm,int mat,int lev)
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].SameNuclide(atm,mat,lev))return i;
  };
  return -1;
};

void BCGManager::ReadFPYieldDataFromFile(string filename,int pnt,string tagname,bool print)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  yield_tagname.push_back(tagname);

  int epnt;
  fin>>epnt;

  if(pnt<0||pnt>=epnt){
    cout<<"Error in BCGManager::ReadFPYieldDataFromFile\n";
    exit(0);
  };

  for(int i=0;i<epnt;i++){

    real e;
    fin>>e;
    int nfp;
    fin>>nfp;

    if(print&&i==pnt)cout<<"# Fission yield data ["<<tagname<<"] is stored for fissions at "<<e<<" (eV).\n";

    //real sum=0.;
    for(int j=0;j<nfp;j++){
      int iz,ia,il;
      fin>>iz;
      fin>>ia;
      fin>>il;
      real yc,dyc;
      fin>>yc;
      fin>>dyc;

      if(i==pnt){
	int id=GetNuclideIndex(iz,ia,il);
	if(id==-1){
	  //cout<<"No decay data : "<<iz<<" "<<ia<<" "<<il<<" "<<yc<<"\n";
	}else{
	  nuc[id].PutYield(tagname,yc,dyc);
	  //sum+=yc;
	};
      };

    };
    //cout<<"# Total : "<<sum<<"\n";

  };

  // (Put zero-yield)
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].GetYieldID(tagname)==-1){
      nuc[i].PutYield(tagname,0.,0.);
    };
  };

};

void BCGManager::WriteFPYieldDataToFile(string filename,string tagname,real eng)
{
  int ii=yield_tagname.size();
  bool exist=false;
  for(int i=0;i<ii;i++){
    if(yield_tagname[i]==tagname)exist=true;
  };
  if(!exist){
    cout<<"# Error in BCGManager::WriteFPYieldDataToFile.\n";
    cout<<"# Tagname for fissile "<<tagname<<" does NOT exist.\n";
    exit(0);
  };
   

  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  fout.setf(ios::scientific);
  fout.precision(10);

  fout<<"  1\n";  // dummy for number of energy point
  fout<<"  "<<eng<<"\n";  // incident energy

  int nfp=nuc.size();
  fout<<nfp<<"\n";

  for(int j=0;j<nfp;j++){
    fout<<nuc[j].GetAtomicNumber()<<"\n";
    fout<<nuc[j].GetMassNumber()<<"\n";
    fout<<nuc[j].GetExLevel()<<"\n";
    fout<<nuc[j].GetYield(tagname)<<"\n";
    fout<<nuc[j].GetDeltaYield(tagname)<<"\n";
  };

  fout.close();
};

void BCGManager::AddNuclide(int iz,int ia,int il)
{
  BCGNuclide nucinp(iz,ia,il);
  nucinp.PutChannel(0);
  int id=GetNuclideIndex(iz,ia,il);
  if(id==-1)nuc.push_back(nucinp);
};

void BCGManager::ReadDecayDataFromFile(string filename,bool stable_only)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  int iz,ia,il;
  iz=0;

  while(iz!=-1){

    fin>>iz;
    if(iz==-1)return;
    fin>>ia;
    fin>>il;

    BCGNuclide nucinp(iz,ia,il);

    real awr;
    fin>>awr;
    nucinp.PutAWR(awr);
    real t,dt;
    fin>>t;
    fin>>dt;
    nucinp.PutHalflife(t,dt);

    real elp,delp,eem,deem,ehp,dehp;
    fin>>elp;
    fin>>delp;
    fin>>eem;
    fin>>deem;
    fin>>ehp;
    fin>>dehp;
    nucinp.PutDecayEnergy(0,elp);
    nucinp.PutDecayEnergy(1,eem);
    nucinp.PutDecayEnergy(2,ehp);
    nucinp.PutDeltaDecayEnergy(0,delp);
    nucinp.PutDeltaDecayEnergy(1,deem);
    nucinp.PutDeltaDecayEnergy(2,dehp);

    int lis;
    fin>>lis;
   
    int ndk;
    fin>>ndk;

    nucinp.PutChannel(ndk);

    real sum=0.;
    for(int i=0;i<ndk;i++){
      real br,dbr;
      fin>>br;
      fin>>dbr;
      sum+=br;
      int iz2,ia2,il2,dnin;
      fin>>iz2;
      fin>>ia2;
      fin>>il2;
      fin>>dnin;
      int tp=-1;
      if(ia==ia2&&iz==iz2){
	tp=2;
      }else if(ia==ia2&&iz+1==iz2){
	tp=0;
      }else if(ia==ia2&&iz-1==iz2){
	tp=1;
      }else if(ia-4==ia2&&iz-2==iz2){
	tp=3;
      }else{
	tp=-1;
      };
      nucinp.PutNextNuclideData(i,iz2,ia2,il2,br,dbr,tp,dnin);
    };
    /*
    cout.setf(ios::showpoint);
    cout.precision(5);
    if(sum<0.999||sum>1.001){
      cout<<iz<<" "<<ia<<" "<<il<<"\n";
      cout<<sum<<"\n";
      exit(0);
    };
    */
    
    int id=GetNuclideIndex(iz,ia,il);
    if(!stable_only||(stable_only&&ndk==0)){
      if(id==-1){
        nuc.push_back(nucinp);
      }else{
        nuc[id]=nucinp;
        //cout<<"# Decay data of nuclide ("<<iz<<","<<ia<<","<<il<<") is overwritten.\n";
      };
    };

  };
    
};

void BCGManager::WriteDecayDataToFile(string filename)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  fout.setf(ios::scientific);
  fout.precision(10);

  int sz=nuc.size();
  for(int ii=0;ii<sz;ii++){

    fout<<nuc[ii].GetAtomicNumber()<<"\n";
    fout<<nuc[ii].GetMassNumber()<<"\n";
    fout<<nuc[ii].GetExLevel()<<"\n";

    fout<<nuc[ii].GetAWR()<<"\n";
    fout<<nuc[ii].GetHalflife()<<"\n";
    fout<<nuc[ii].GetDeltaHalflife()<<"\n";

    fout<<nuc[ii].GetDecayEnergy(0)<<"\n";
    fout<<nuc[ii].GetDeltaDecayEnergy(0)<<"\n";
    fout<<nuc[ii].GetDecayEnergy(1)<<"\n";
    fout<<nuc[ii].GetDeltaDecayEnergy(1)<<"\n";
    fout<<nuc[ii].GetDecayEnergy(2)<<"\n";
    fout<<nuc[ii].GetDeltaDecayEnergy(2)<<"\n";

    fout<<"  0\n"; // dummy for [lis]

    int ndk=nuc[ii].GetChannel();
    fout<<ndk<<"\n";

    for(int i=0;i<ndk;i++){
      fout<<nuc[ii].GetBr(i)<<"\n";
      fout<<nuc[ii].GetDeltaBr(i)<<"\n";
      fout<<nuc[ii].GetAtomicNumberNext(i)<<"\n";
      fout<<nuc[ii].GetMassNumberNext(i)<<"\n";
      fout<<nuc[ii].GetExLevelNext(i)<<"\n";
      fout<<nuc[ii].GetEmittedNeutron(i)<<"\n";
    };

  };
    
  fout<<" -1\n";
  fout.close();
};

void BCGManager::SetNGBranchingRatioDataForFR()
{
  // Based on ORIGEN library for fast reactor core region.
  // Values are taken from JAERI-Data/Code 2004-015, p.32-33.

  int nucn=110;
  string nucnam[]={
    "Na023","Cl037","Sc045","Co059","Zn068", "Zn070","Ga071","Ge070","Ge074","Ge076",
    "Se076","Se078","Se080","Se082","Br079", "Br081","Kr078","Kr080","Kr082","Kr084",
    "Rb085","Sr084","Sr086","Y089","Y090", "Nb093","Mo092","Rh103","Rh105","Pd106",
    "Pd108","Pd110","Ag107","Ag109","Ag110m", "Cd110","Cd112","Cd114","Cd116","Cd118",
    "In113","In115","In117","In117m","In119", "In119m","Sn112","Sn116","Sn118","Sn120",
    "Sn122","Sn124","Sn126","Sb121","Sb123", "Sb125","Te120","Te122","Te124","Te126",
    "Te128","Te130","I129","Xe124","Xe126", "Xe128","Xe130","Xe132","Xe134","Cs133",
    "Cs134","Cs134m","Ba130","Ba132","Ba134", "Ba135","Ba136","Ce136","Ce138","Pr141",
    "Pm147","Eu151","Gd154","Dy164","Ho165", "Er166","Tm169","Yb174","Lu175","Lu176",
    "Hf177","Hf178","Hf179","Ta181","W182", "W184","Re187","Os190","Os192","Ir191",
    "Ir193","Pt192","Pt194","Pt196","Pt198", "Hg196","Hg198","Bi209","Po210","Pa233",
  };
  real br[]={ // Branching ratio to META-STABLE state
    2.033e-1, 1.290e-2, 3.622e-1, 5.352e-1, 6.721e-2,
    9.489e-2, 3.090e-2, 8.159e-2, 8.178e-3, 9.624e-2,
    9.970e-2, 9.014e-2, 6.372e-3, 4.395e-4, 7.932e-2,
    1.751e-1, 4.460e-2, 1.323e-1, 4.955e-1, 5.449e-2,
    6.187e-3, 6.789e-1, 8.301e-2, 5.448e-5, 8.690e-3,
    3.774e-3, 5.165e-3, 1.606e-1, 9.307e-1, 2.459e-3,
    3.833e-2, 1.677e-2, 1.753e-2, 1.904e-1, 5.000e-1,
    2.231e-3, 8.258e-2, 1.183e-2, 1.070e-3, 5.000e-1,
    3.341e-1, 9.037e-1, 5.000e-1, 5.000e-1, 6.000e-1,
    6.000e-1, 3.044e-1, 2.279e-1, 8.989e-2, 3.921e-4,
    2.300e-4, 2.709e-1, 1.273e-2, 6.603e-3, 2.149e-3,
    1.754e-3, 1.452e-1, 1.301e-1, 2.772e-4, 1.984e-2,
    1.761e-3, 4.451e-3, 9.419e-2, 1.719e-1, 6.501e-2,
    7.590e-3, 1.324e-2, 1.246e-3, 1.845e-4, 1.075e-1,
    9.240e-2, 5.000e-1, 1.853e-1, 7.409e-2, 2.598e-2,
    1.148e-3, 1.127e-3, 1.310e-1, 1.349e-2, 3.778e-2,
    5.701e-1, 3.474e-1, 2.055e-4, 6.025e-1, 3.775e-2,
    1.336e-1, 6.640e-2, 7.076e-1, 7.009e-1, 3.301e-3,
    2.999e-3, 6.164e-1, 7.600e-3, 4.899e-4, 3.787e-4,
    1.101e-3, 9.786e-1, 7.000e-1, 7.407e-3, 4.099e-4,
    5.008e-2, 1.572e-1, 7.500e-2, 6.762e-2, 7.300e-3,
    1.248e-1, 9.495e-3, 4.243e-1, 1.640e-2, 5.000e-1
  };

  real bri[]={0.,0.};
  int iz=0;
  int ia=0;
  int il=0;
  for(int i=0;i<nucn;i++){
    midt.GetParameter(nucnam[i],iz,ia,il);
    if(GetNuclideIndex(iz,ia,il)!=-1){
      if(GetNuclideIndex(iz,ia+1,0)!=-1&&GetNuclideIndex(iz,ia+1,1)!=-1){
	bri[0]=1.-br[i];
	bri[1]=br[i];
	GetNuclide(iz,ia).PutReactionData(0,2,bri);
      };
    };
  };

};

void BCGManager::ReadNGBranchingRatioDataFromFile(string filename)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  int tmp1,tmp2,tmp3;
  real vtmp;

  int iz=0;
  int ia=0;

  vector<real> bri;
  bri.clear();

  tmp1=0;

  while(tmp1!=-1){
    fin>>tmp1;
    if(tmp1==-1){
      int sz=bri.size();
      real *bri2=new real[sz];
      for(int i=0;i<sz;i++){
        bri2[i]=bri[i];
      };
      GetNuclide(iz,ia-1).PutReactionData(0,sz,bri2);
      delete [] bri2;
      return;
    };
    fin>>tmp2;
    fin>>tmp3;
    fin>>vtmp;
    if(tmp3==0){
      int sz=bri.size();
      if(sz!=0){
        real *bri2=new real[sz];
        for(int i=0;i<sz;i++){
          bri2[i]=bri[i];
        };
        GetNuclide(iz,ia-1).PutReactionData(0,sz,bri2);
        delete [] bri2;
        bri.clear();
      };
      bri.push_back(vtmp);
      iz=tmp1;
      ia=tmp2;
    }else{
      bri.push_back(vtmp);
    };
  };
};

void BCGManager::ReadBranchingRatioDataFromFile(string filename)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(0);
  };

  int iz,im,il;
  int mat=0;
  while(mat!=-1){
    fin>>mat;
    if(mat==-1)return;
    int flag;
    fin>>flag; // 0:capture, 1/(n2n)
    int ich;
    fin>>ich;
    real *bri=new real[ich];    
    for(int i=0;i<ich;i++){
      real tmp;
      fin>>tmp;
      bri[i]=tmp;
    };
    int matid=midt.GetMATIDFromENDFID(mat);
    midt.GetParameter(matid,iz,im,il);
    int id=GetNuclideIndex(iz,im,il);
    if(id!=-1){
      nuc[id].PutReactionData(flag,ich,bri);
    };
    delete [] bri;
  };


};

void BCGManager::ShowNuclideList(string yld_tag)
{ 
  int flag_count=0;
  int sz=nuc.size();
  cout<<"#\n# Nuclide List in BCGManager\n#\n";
  cout<<"#    Ch : Number of decay channels\n#\n";
  cout<<"#     No.     Z      M  Lv  Ch  Half-life[s]\n";  
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()){
      cout<<" o ";
    }else{
      cout<<" x ";
    };
    int at=nuc[i].GetAtomicNumber();
    int ms=nuc[i].GetMassNumber();
    int lv=nuc[i].GetExLevel();
    int ch=nuc[i].GetChannel();
    real halflife=nuc[i].GetHalflife();
    string name=midt.NameFromAtomicNumber(at);
    WriteOut(i,5);
    WriteOut(at,5);
    cout<<"[";
    WriteOut(name,2);
    cout<<"]";
    WriteOut(ms,5);
    WriteOut(lv,4);
    WriteOut(ch,4);
    cout<<"   "<<halflife<<" ";
    /*
    int id=nuc[i].GetYieldID(yld_tag);
    if(id!=-1){
      cout<<nuc[i].GetYield(yld_tag)<<" ";
    };
    */
    //cout<<nuc[i].Flag()<<" ";
    /*
    for(int j=0;j<nuc[i].GetYieldNum();j++){
      cout<<nuc[i].GetYieldTag(j)<<" ";
    };
    */
    if(nuc[i].Flag())flag_count++;
    cout<<"\n";
  };
  cout<<"#\n# Flaged nuclide   : "<<flag_count<<"\n";
  cout<<"# Unflaged nuclide : "<<sz-flag_count<<"\n";
};

real BCGManager::GetYieldSumFlag(string tagname,bool flag)
{
  real sum=0.;
  for(int i=0;i<GetNuclideNumber();i++){
    if(flag&&nuc[i].Flag())sum+=nuc[i].GetYield(tagname);
    if(!flag&&!nuc[i].Flag())sum+=nuc[i].GetYield(tagname);
  };
  return sum;
};

real BCGManager::GetYieldSum(string tagname)
{
  real sum=0.;
  for(int i=0;i<GetNuclideNumber();i++){
    sum+=nuc[i].GetYield(tagname);
    /*
    if(nuc[i].GetYield(tagname)!=0.){
      cout<<nuc[i].GetID()<<" "<<nuc[i].GetYield(tagname)<<"\n";
    };
    */
  };
  return sum;
};

void BCGManager::ShowYieldSum()
{
  int sz=yield_tagname.size();
  cout.setf(ios::showpoint);
  cout.precision(7);
  cout<<"#############################################################\n";
  cout<<"# Summation of yield (included/excluded)\n#\n";
  for(int i=0;i<sz;i++){
    cout<<"# "<<yield_tagname[i]<<" : ";
    cout<<GetYieldSum(yield_tagname[i])<<" ";
    cout<<"( "<<GetYieldSumFlag(yield_tagname[i],true)<<" / ";
    cout<<GetYieldSumFlag(yield_tagname[i],false)<<" )\n";
  };
  cout<<"#############################################################\n\n";
};

void BCGManager::ShowYieldSum(string tagname)
{
  cout.setf(ios::showpoint);
  cout.precision(10);
  cout<<"Sum of yields : "<<GetYieldSum(tagname)<<"\n";
};

void BCGManager::ShowNumberOfTotalEmittedNeutrons()
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    real tmp=nuc[i].GetTotalNumberEmittedNeutrons();
    if(tmp>0.){
      cout<<nuc[i].GetID()<<" "<<tmp<<"\n";
    };
  };
};

void BCGManager::MakeFlagFalseAllNuclide()
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    nuc[i].MakeFlagFalse();
  };
};

void BCGManager::MakeFlagTrueAllNuclide()
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    nuc[i].MakeFlagTrue();
  };
};

void BCGManager::MakeFlagTrue(int atm,int mas,int lev)
{
  int id=GetNuclideIndex(atm,mas,lev);
  if(id==-1){
    cout<<"Warning in BCGManager::MakeFlagTrue\n";
    cout<<"No nuclide data: "<<atm<<" "<<mas<<" "<<lev<<"\n";
  }else{
    nuc[id].MakeFlagTrue();
  };
};

void BCGManager::MakeFlagTrueAtom(int atm)
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].GetAtomicNumber()==atm){
      nuc[i].MakeFlagTrue();
    };
  };
};

void BCGManager::MakeFlagTrue(int num,string *nuc_nam)
{
  //MakeFlagFalseAllNuclide();

  for(int i=0;i<num;i++){
    int i1,i2,i3;
    midt.GetParameter(nuc_nam[i],i1,i2,i3);
    MakeFlagTrue(i1,i2,i3);
  };
};

void BCGManager::MakeFlagFalse(int num,string *nuc_nam)
{
  //MakeFlagFalseAllNuclide();

  for(int i=0;i<num;i++){
    int i1,i2,i3;
    midt.GetParameter(nuc_nam[i],i1,i2,i3);
    MakeFlagFalse(i1,i2,i3);
  };
};

void BCGManager::MakeFlagTrue(int num,int *matid)
{
  MakeFlagFalseAllNuclide();

  for(int i=0;i<num;i++){
    int i1,i2,i3;
    midt.GetParameter(matid[i],i1,i2,i3);
    MakeFlagTrue(i1,i2,i3);
  };
};

void BCGManager::MakeFlagReversed()
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()){
      nuc[i].MakeFlagFalse();
    }else{
      nuc[i].MakeFlagTrue();
    };
  };
};

void BCGManager::CalCumulativeDecayEnergy()
{
  int sz=nuc.size();
  vector<real> den(sz);
  vector<real> eng(3);

  for(int i=0;i<sz;i++){

    for(int k=0;k<sz;k++){
      den[k]=0.;
    };

    eng[0]=nuc[i].GetDecayEnergy(0);
    eng[1]=nuc[i].GetDecayEnergy(1);
    eng[2]=nuc[i].GetDecayEnergy(2);

    int div=nuc[i].GetChannel();
    for(int k=0;k<div;k++){
      int am=nuc[i].GetAtomicNumberNext(k);
      int ms=nuc[i].GetMassNumberNext(k);
      int lv=nuc[i].GetExLevelNext(k);
      real br=nuc[i].GetBr(k);
      int id=GetNuclideIndex(am,ms,lv);
      if(id!=-1&&!nuc[id].Flag()){
        den[id]=br;
      };
    };

    bool end=false;
    while(!end){
      end=true;
      for(int j=0;j<sz;j++){
        if(den[j]!=0.&&!nuc[j].Flag()){
          eng[0]+=nuc[j].GetDecayEnergy(0)*den[j];          
          eng[1]+=nuc[j].GetDecayEnergy(1)*den[j];          
          eng[2]+=nuc[j].GetDecayEnergy(2)*den[j];          
          int div=nuc[j].GetChannel();
          for(int k=0;k<div;k++){
            int am=nuc[j].GetAtomicNumberNext(k);
            int ms=nuc[j].GetMassNumberNext(k);
            int lv=nuc[j].GetExLevelNext(k);
            real br=nuc[j].GetBr(k);
            int id=GetNuclideIndex(am,ms,lv);
            if(id!=-1&&!nuc[id].Flag()){
              end=false;
              den[id]=den[j]*br;
            };
	  };
          den[j]=0.;
	};
      };
    };

    nuc[i].PutDecayEnergy(0,eng[0]);
    nuc[i].PutDecayEnergy(1,eng[1]);
    nuc[i].PutDecayEnergy(2,eng[2]);
  };

  
};

void BCGManager::CalCumulativeYield(real half_life_limit, bool print)
{
  int sz=nuc.size();
  int szy=yield_tagname.size();
  bool end=false;

  int nuc0=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag())nuc0++;
  };
  if(print){
  cout<<"#############################################################\n";
  cout<<"# Cumulative yield calculation by BurnupChainManager.\n#\n";
  cout<<"# The number of original nuclides : "<<sz<<"\n";
  cout<<"# The reduced number              : "<<nuc0<<"\n";
  };

  int iter=0;
  real half_life_max=0.;
  real half_life_min=1e10;
  int am1=0;
  int ms1=0;
  int lv1=0;
  while(!end){
    iter++;
    end=true;
    for(int i=0;i<sz;i++){
      bool flag=nuc[i].Flag();
      int div=nuc[i].GetChannel();
      int amb=nuc[i].GetAtomicNumber();
      int msb=nuc[i].GetMassNumber();
      int lvb=nuc[i].GetExLevel();
      real half_life=nuc[i].GetHalflife();
      if(!nuc[i].ZeroYield()&&!flag&&div!=0&&(half_life_limit==0.||half_life<half_life_limit)){
        end=false;
        for(int j=0;j<div;j++){
          int am=nuc[i].GetAtomicNumberNext(j);
          int ms=nuc[i].GetMassNumberNext(j);
          int lv=nuc[i].GetExLevelNext(j);
          real br=nuc[i].GetBr(j);
          int id=GetNuclideIndex(am,ms,lv);
	  if(id==-1){
	    if(print){
	      cout<<"# (No nuclide data in data base : "<<am<<" "<<ms<<" "<<lv<<")\n";
	    };
	  }else{
	    if(half_life>half_life_max){
              half_life_max=half_life;
	      am1=amb;
	      ms1=msb;
	      lv1=lvb;
	    };
	    if(half_life<half_life_min&&half_life!=0.){
	      half_life_min=half_life;
	    };
	    for(int k=0;k<szy;k++){
              string tag=nuc[i].GetYieldTag(k);
	      real yld=nuc[i].GetYield(tag);
              real d_yld=nuc[i].GetDeltaYield(tag);
  	      nuc[id].AddYield(tag,br*yld,br*d_yld);
	      if(j==div-1){
		nuc[i].PutYield(tag,0.,0.);
	      };
	    };
	  };
        };
      };
    };
  };

  if(print){
  cout<<"# Maximum halflife is "<<half_life_max<<" [s]\n";
  cout<<"# Minimum halflife is "<<half_life_min<<" [s]\n";
  //cout<<am1<<" "<<ms1<<" "<<lv1<<"\n";
  cout<<"#############################################################\n\n";
  };
};

void BCGManager::CalReactionBranch()
{
  int sz=nuc.size();
  int r_num=nuc[0].GetRnum();

  for(int iii=0;iii<r_num;iii++){

  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()&&nuc[i].GetReactionChannel(iii)!=0){
      bool end=false;
      while(!end){
        end=true;
        int ngc=nuc[i].GetReactionChannel(iii);        
        int ngc2=0;
        vector<int> atn;
        vector<int> mas;
        vector<int> lev;
        vector<real> bri;
	for(int j=0;j<ngc;j++){
          int am=nuc[i].GetReactionAtomicNumberNext(iii,j);
          int ms=nuc[i].GetReactionMassNumberNext(iii,j);
          int lv=nuc[i].GetReactionExLevelNext(iii,j);
          real br=nuc[i].GetReactionBr(iii,j);
          int id=GetNuclideIndex(am,ms,lv);
	  if(id!=-1){       
          bool fl=GetNuclide(am,ms,lv).Flag();
	  if(fl){
            bool exist=false;
            for(int l=0;l<ngc2;l++){
              if(am==atn[l]&&ms==mas[l]&&lv==lev[l]){
	        exist=true;
                bri[l]+=br;
	      };
	    };
            if(!exist){
  	      atn.push_back(am);
              mas.push_back(ms);
              lev.push_back(lv);
              bri.push_back(br);
              ngc2++;
	    };
	  }else{
            end=false;
            int id2=GetNuclideIndex(am,ms,lv);
            int ch=nuc[id2].GetChannel();
            for(int k=0;k<ch;k++){
              int am2=nuc[id2].GetAtomicNumberNext(k);
              int ms2=nuc[id2].GetMassNumberNext(k);
              int lv2=nuc[id2].GetExLevelNext(k);
              real br2=nuc[id2].GetBr(k);       
              bool exist=false;
              for(int l=0;l<ngc2;l++){
                if(am2==atn[l]&&ms2==mas[l]&&lv2==lev[l]){
		  exist=true;
                  bri[l]+=br*br2;
		};
	      };
              if(!exist){
                atn.push_back(am2);
	        mas.push_back(ms2);
	        lev.push_back(lv2);
                bri.push_back(br*br2);
	        ngc2++;
	      };
	    };
	  };
	  };
	};
        nuc[i].PutReactionData(iii,ngc2,atn,mas,lev,bri);
      };
    };
  };

  }; // loop for reaction type
};

void BCGManager::CalDecayBranch(real half_life_limit, bool print)
{
  // ! Notice !
  //
  // The number of emitted delayed neutrons, dn[i], is NOT well calculated
  // when daughter nuclides are NOT considered in the burnup chain.
  // It is because such decay path is ignored in this method.

  if(print){
  cout<<"#############################################################\n";
  cout<<"# Decay branching calculation\n#\n";
  };

  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()&&nuc[i].GetChannel()!=0){

      bool end=false;
      int iter=0;

      while(!end){
        end=true;
        int ngc=nuc[i].GetChannel();        
        int ngc2=0;
        vector<int> atn;
        vector<int> mas;
        vector<int> lev;
        vector<real> bri;
        vector<real> dbri;
        vector<real> dni;
	for(int j=0;j<ngc;j++){
          int am=nuc[i].GetAtomicNumberNext(j);
          int ms=nuc[i].GetMassNumberNext(j);
          int lv=nuc[i].GetExLevelNext(j);
          real br=nuc[i].GetBr(j); 
          real dbr=nuc[i].GetDeltaBr(j); 
          real dn=nuc[i].GetEmittedNeutron(j);   

	  if(GetNuclideIndex(am,ms,lv)!=-1){
          bool fl=GetNuclide(am,ms,lv).Flag();
          if(fl){
            bool exist=false;
            for(int l=0;l<ngc2;l++){
              if(am==atn[l]&&ms==mas[l]&&lv==lev[l]){
                exist=true;
                dni[l]=(dni[l]*bri[l]+dn*br)/(bri[l]+br);
                bri[l]+=br;
              };
            };
            if(!exist){
              atn.push_back(am);
              mas.push_back(ms);
              lev.push_back(lv);
              bri.push_back(br);
              dbri.push_back(dbr);
              dni.push_back(dn);
              ngc2++;
	    };
	  }else{
            end=false;
            int id2=GetNuclideIndex(am,ms,lv);
            int ch=nuc[id2].GetChannel();
            real half_life=nuc[id2].GetHalflife();
	    if(half_life>half_life_limit&&half_life_limit!=0.){
	      if(print){
	      cout<<"# Half-life limit is over : "<<midt.Name(am,ms,lv)<<" ";
	      cout<<"("<<half_life<<"[s])\n";
	      };
	    };
	    if(half_life_limit==0||half_life<half_life_limit){
              for(int k=0;k<ch;k++){
                int am2=nuc[id2].GetAtomicNumberNext(k);
                int ms2=nuc[id2].GetMassNumberNext(k);
                int lv2=nuc[id2].GetExLevelNext(k);
                real br2=nuc[id2].GetBr(k);       
                real dn2=nuc[id2].GetEmittedNeutron(k);       
                bool exist=false;
                for(int l=0;l<ngc2;l++){
                  if(am2==atn[l]&&ms2==mas[l]&&lv2==lev[l]){
		    exist=true;
                    dni[l]=(dni[l]*bri[l]+(dn+dn2)*(br*br2))/(bri[l]+br*br2);
                    bri[l]+=br*br2;
		  };
	        };
                if(!exist){
                  atn.push_back(am2);
	          mas.push_back(ms2);
	          lev.push_back(lv2);
                  bri.push_back(br*br2);
                  dbri.push_back(0.);
                  dni.push_back(dn+dn2);
	          ngc2++;
		};
	      };
	    };
	  };
	};

	};
        nuc[i].PutDecayData(ngc2,atn,mas,lev,bri,dbri,dni);
        iter++;
      };
      //cout<<i<<" "<<iter<<"\n";
    };

  };

  if(print){
  cout<<"#############################################################\n\n";
  };
};

void BCGManager::CalTotalEmittedNeutrons(bool repeating_cal)
{
  // Total number of emitted neutrons is calculated.
  // Neutron emittions of its daughter nuclides are also considered.
  //
  // repeating_cal : this is for sensitivity calculation.


  //cout<<"#############################################################\n";
  //cout<<"# Total emitted neutron calculation\n#\n";
  //int yld_num=yield_tagname.size();
  //vector<real> sumy(yld_num,0.);

  int sz=nuc.size();
  for(int i=0;i<sz;i++){

    if(repeating_cal&&nuc[i].GetTotalNumberEmittedNeutrons()==0.){
    }else{

    int at_org=nuc[i].GetAtomicNumber();
    int ms_org=nuc[i].GetMassNumber();
    int lv_org=nuc[i].GetExLevel();

    bool end=false;
    int iter=0;

    real n_emit=0.;

    // +++ initialize
    int ngc=nuc[i].GetChannel(); 
    int ngc_org=ngc;       
    vector<real> n_emit_ch(ngc,0.);
    vector<int> atn(ngc);
    vector<int> mas(ngc);
    vector<int> lev(ngc);
    vector<real> bri(ngc);
    vector<int> ch_id(ngc);
    for(int j=0;j<ngc;j++){
      atn[j]=nuc[i].GetAtomicNumberNext(j);
      mas[j]=nuc[i].GetMassNumberNext(j);
      lev[j]=nuc[i].GetExLevelNext(j);
      bri[j]=nuc[i].GetBr(j); 
      ch_id[j]=j;
      n_emit+=bri[j]*nuc[i].GetEmittedNeutron(j);
      n_emit_ch[j]+=bri[j]*nuc[i].GetEmittedNeutron(j);
    };

    while(!end){
      int ngc2=0;
      vector<int> atn2;
      vector<int> mas2;
      vector<int> lev2;
      vector<real> bri2;
      vector<int> ch_id2;

      end=true;
      for(int j=0;j<ngc;j++){
        int id=GetNuclideIndex(atn[j],mas[j],lev[j]);
        int ch_org=ch_id[j];
	if(id!=-1){        
          int ngc3=nuc[id].GetChannel();
          for(int l=0;l<ngc3;l++){
            int am=nuc[id].GetAtomicNumberNext(l);
  	    int ms=nuc[id].GetMassNumberNext(l);
	    int lv=nuc[id].GetExLevelNext(l);
	    if(am!=atn[j]||ms!=mas[j]||lv!=lev[j]){
	      end=false;
              real br=nuc[id].GetBr(l);
              real dn=nuc[id].GetEmittedNeutron(l);   
  	      n_emit+=(bri[j]*br)*dn;          
	      n_emit_ch[ch_org]+=(bri[j]*br)*dn;
              bool exist=false;
	      for(int m=0;m<ngc2;m++){
                if(am==atn2[m]&&ms==mas2[m]&&lv==lev2[m]){
                  exist=true;
                  bri2[m]+=bri[j]*br;
                };
              };
              if(!exist){
                atn2.push_back(am);
                mas2.push_back(ms);
                lev2.push_back(lv);
                bri2.push_back(bri[j]*br);
                ch_id2.push_back(ch_org);
                ngc2++;
	      };
	    };
	  };
	};
      };
      ngc=ngc2;
      atn.resize(ngc);
      mas.resize(ngc);
      lev.resize(ngc);
      bri.resize(ngc);
      ch_id.resize(ngc);
      for(int j=0;j<ngc;j++){
	atn[j]=atn2[j];
	mas[j]=mas2[j];
	lev[j]=lev2[j];
	bri[j]=bri2[j];
	ch_id[j]=ch_id2[j];
      };
      iter++;
    };

    nuc[i].PutTotalNumberEmittedNeutrons(n_emit);
    real sum=0.;
    for(int j=0;j<ngc_org;j++){
      sum+=n_emit_ch[j];
    };
    if(fabs(sum-n_emit)>1e-5){
      cout<<sum<<" "<<n_emit<<" "<<ngc<<"\n";
      exit(0);
    };
    nuc[i].PutTotalNumberEmittedNeutronsCh(ngc_org,n_emit_ch);

    };
    //cout<<i<<" "<<at_org<<" "<<ms_org<<" "<<lv_org<<" "<<n_emit<<"\n";
    /*
    for(int j=0;j<yld_num;j++){
      sumy[j]+=n_emit*nuc[i].GetYield("tmp");
    };
    */
  };
  //cout<<"# Total number of delayed neutron : "<<sumy[0]<<"\n";
  //cout<<"#############################################################\n\n";
};

void BCGManager::ShowYieldData(string yield_tag,bool flaged)
{
  cout<<"#\n# Yield data for ";
  if(flaged){
    cout<<"flaged ";
  }else{
    cout<<"not-flaged ";
  };
  cout<<"nuclides.\n";
  cout<<"# Fissile nuclide : "<<yield_tag<<"\n#\n";

  int sz=nuc.size();
  int nn=0;
  real sum=0.;
  real sum2=0.;
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<sz;i++){
    real yld=nuc[i].GetYield(yield_tag);
    bool flag=nuc[i].Flag();
    sum2+=yld;
    if((flaged&&flag&&yld>0.)||(!flaged&&!flag&&yld>0.)){
      //if((flaged&&flag)||(!flaged&&!flag)){
      int at=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      string name=midt.NameFromAtomicNumber(at);
      WriteOut(nn,5);
      WriteOut(at,5);
      cout<<"[";
      WriteOut(name,2);
      cout<<"]";
      WriteOut(ms,5);
      WriteOut(lv,4);
      cout<<"   "<<yld<<" ";
      cout<<"\n";
      nn++;
      sum+=yld;
    };
  };
  cout.setf(ios::showpoint);
  cout.precision(10);
  cout<<"# Sum of target nuclides : "<<sum<<"\n";
  cout<<"# Sum of all nuclides    : "<<sum2<<"\n";

  // +++ temporaty treatment 
  vector<string> name_store;
  vector<real> yld_store;
  int count=0;

  for(int i=0;i<sz;i++){
    real yld=nuc[i].GetYield(yield_tag);
    bool flag=nuc[i].Flag();
    if((flaged&&flag&&yld>0.)||(!flaged&&!flag&&yld>0.)){
      //if((flaged&&flag)||(!flaged&&!flag)){
      int at=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      string name=midt.Name(at,ms,lv);
      name_store.push_back(name);
      yld_store.push_back(yld);
      count++;
    };
  };

  cout<<"# Number of nuclides : "<<count<<"\n";
  for(int i=0;i<count;i++){
    cout<<"\""<<name_store[i]<<"\", ";
    if((i+1)%5==0)cout<<"\n";
  };
  cout<<"\n";
  cout.setf(ios::scientific);
  cout.precision(6);
  for(int i=0;i<count;i++){
    cout<<yld_store[i]<<", ";
    if((i+1)%5==0)cout<<"\n";
  };
  cout<<"\n";
 
};



void BCGManager::BranchingRatioSumCheck()
{
  cout.setf(ios::showpoint);
  cout.precision(10);

  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    int div=nuc[i].GetChannel();
    if(div!=0){
      real sum=0.;
      for(int j=0;j<div;j++){
	sum+=nuc[i].GetBr(j);
      };
      if(sum<0.9999999||sum>1.0000001){
	cout<<"# Sum of BR is not 1. ("<<sum<<")  ";
	cout<<nuc[i].GetAtomicNumber()<<" ";
	cout<<nuc[i].GetMassNumber()<<" ";
	cout<<nuc[i].GetExLevel()<<"\n";
      };
    };
  };
};

void BCGManager::IgnoreNGReaction(int num,string *nuc_nam)
{
  for(int i=0;i<num;i++){
    int atm,mas,lev;
    midt.GetParameter(nuc_nam[i],atm,mas,lev);
    GetNuclide(atm,mas,lev).PutReactionChannel(0,0);
  };
};

real BCGManager::GetYield(string tagname,int atm,int mas,int lev)
{
  int id=GetNuclideIndex(atm,mas,lev);
  if(id==-1){
    cout<<"# Error in BCGManager::GetYield.\n";
    cout<<"# No data in BCGManager : (IA,IM,IZ)=("<<atm<<", "<<mas<<", "<<lev<<")\n";
    exit(0);
  };
  return nuc[id].GetYield(tagname);
};

real BCGManager::GetDeltaYield(string tagname,int atm,int mas,int lev)
{
  int id=GetNuclideIndex(atm,mas,lev);
  if(id==-1){
    cout<<"# Error in BCGManager::GetDeltaYield.\n";
    cout<<"# No data in BCGManager : (IA,IM,IZ)=("<<atm<<", "<<mas<<", "<<lev<<")\n";
    exit(0);
  };
  return nuc[id].GetDeltaYield(tagname);
};

void BCGManager::ShowNGBranch()
{
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    int ch=nuc[i].GetReactionChannel(0);
    if(ch>0){
      cout<<nuc[i].GetAtomicNumber()<<" ";
      cout<<nuc[i].GetMassNumber()<<" ";
      cout<<nuc[i].GetExLevel()<<"\n";
      for(int j=0;j<ch;j++){
	cout<<"    ";
        cout<<nuc[i].GetReactionAtomicNumberNext(0,j)<<" ";
        cout<<nuc[i].GetReactionMassNumberNext(0,j)<<" ";
        cout<<nuc[i].GetReactionExLevelNext(0,j)<<" : ";
        cout<<nuc[i].GetReactionBr(0,j)<<"\n";
      };
    };
  };
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     FILE OUTPUT
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void BCGManager::WriteFileChainDataCBGFormat(string filename,bool ng_only)
{
  int num=0;
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag())num++;
  };

  string *nuc_num=new string[num];
  
  int cnt=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()){
      int atm=nuc[i].GetAtomicNumber();
      int mas=nuc[i].GetMassNumber();
      int lev=nuc[i].GetExLevel();
      nuc_num[cnt++]=midt.GetName(atm,mas,lev);
    };
  };

  WriteFileChainDataCBGFormat(num,nuc_num,filename,ng_only);

  delete [] nuc_num;
};

void BCGManager::WriteFileChainDataCBGFormat(int num,string *nuc_nam,string filename,bool ng_only)
{
  // if [ng_only] is true, (n,2n) reaction is ignored.

  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  int r_num=nuc[0].GetRnum();
  vector<int> div_r(r_num);
  vector<int> div_r_2(r_num);
  int div_d_2;

  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<num;i++){
    int atm1,mas1,lev1;
    midt.GetParameter(nuc_nam[i],atm1,mas1,lev1);
    int ii=GetNuclideIndex(atm1,mas1,lev1);
    fout<<midt.ID(atm1,mas1,lev1)<<"\n";
    fout<<GetNuclide(ii).GetHalflife()<<"\n";
    fout<<GetNuclide(ii).GetDecayEnergy(0)<<"\n";
    fout<<GetNuclide(ii).GetDecayEnergy(1)<<"\n";
    fout<<GetNuclide(ii).GetDecayEnergy(2)<<"\n";
    int div_d=GetNuclide(ii).GetChannel(); // The number of decay channels
    for(int r=0;r<r_num;r++){
      div_r[r]=GetNuclide(ii).GetReactionChannel(r); // The number of reaction channels
      if(ng_only&&r!=0)div_r[r]=0;
    };
    fout<<" 0\n"; // for fission

    // ????????
    for(int r=0;r<r_num;r++){
      div_r_2[r]=div_r[r];
      if(div_r[r]==1){
        real br=GetNuclide(ii).GetReactionBr(0,0);
        if(br<0.999999)div_r_2[r]*=-1;
      };
    };

    fout<<" "<<div_r_2[0]<<"\n"; // for (n,g)
    fout<<" "<<div_r_2[1]<<"\n"; // for (n,2n)
    div_d_2=div_d;
    if(div_d==1){
      real br=GetNuclide(ii).GetBr(0);
      if(br<0.999999)div_d_2*=-1;
    };
    fout<<" "<<div_d_2<<"\n";

    for(int r=0;r<r_num;r++){
      for(int k=0;k<div_r[r];k++){
        int atm2=GetNuclide(ii).GetReactionAtomicNumberNext(r,k);
        int mas2=GetNuclide(ii).GetReactionMassNumberNext(r,k);
        int lev2=GetNuclide(ii).GetReactionExLevelNext(r,k);
        real br=GetNuclide(ii).GetReactionBr(r,k);
        fout<<midt.ID(atm2,mas2,lev2)<<"\n";
        if(div_r[r]>1||div_r_2[r]==-1)fout<<br<<"\n";
      };
    };
    for(int j=0;j<div_d;j++){
      int atm2=GetNuclide(ii).GetAtomicNumberNext(j);
      int mas2=GetNuclide(ii).GetMassNumberNext(j);
      int lev2=GetNuclide(ii).GetExLevelNext(j);
      real br=GetNuclide(ii).GetBr(j);
      fout<<midt.ID(atm2,mas2,lev2)<<"\n";
      if(div_d>1||div_d_2==-1)fout<<br<<"\n";
    };
  };
  fout<<"  -1\n";

};

void BCGManager::WriteFileYieldDataCBGFormat(string filename,bool detail)
{
  int num=0;
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag())num++;
  };

  string *nuc_num=new string[num];
  
  int cnt=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()){
      int atm=nuc[i].GetAtomicNumber();
      int mas=nuc[i].GetMassNumber();
      int lev=nuc[i].GetExLevel();
      nuc_num[cnt++]=midt.Name(atm,mas,lev);
    };
  };

  if(detail){
    WriteFileYieldDataCBGFormatDetail(num,nuc_num,filename);
    //WriteFileYieldDataCBGFormatCmSF(num,nuc_num,filename);
  }else{
    WriteFileYieldDataCBGFormat(num,nuc_num,filename);
  };

  delete [] nuc_num;
};



void BCGManager::WriteFileYieldDataCBGFormat(int num,string *nuc_nam,string filename)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  int ynum=12;
  string ytag[]={
    "Th232","U233","U235","U236","U238",
    "Np237","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am243"
  };
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int j=0;j<ynum;j++){
    if(ytag[j]=="U235"){
      fout<<" 2\n 922350\n 932360\n"; // for Np-236
    }else if(ytag[j]=="U236"){
      fout<<" 3\n 922340\n 922370\n 942360\n"; // for U-234, U-237, Pu-236
    }else if(ytag[j]=="U238"){
      fout<<" 2\n 922380\n 942380\n"; // for Pu-238
    }else if(ytag[j]=="Np237"){
      fout<<" 2\n 932370\n 932390\n"; // for Np-239
    }else if(ytag[j]=="Pu241"){
      fout<<" 5\n 942410\n 952420\n 952421\n 962430\n 962450\n"; // for Am-242, -242m, Cm-243, Cm-245
    }else if(ytag[j]=="Pu242"){
      fout<<" 4\n 942420\n 962420\n 962440\n 962460\n"; // for Cm-242, -244, -246
    }else if(ytag[j]=="Th232"){
      fout<<" 4\n 902320\n 912310\n 912330\n 922320\n"; // for Pa-231, -233, U-232
    }else{
      fout<<" 1\n "<<midt.ID(ytag[j])<<"\n";
    };
    fout<<" "<<num<<"\n";
    for(int k=0;k<num;k++){
      int atm,mas,lev;
      midt.GetParameter(nuc_nam[k],atm,mas,lev);
      real dat=GetYield(ytag[j],atm,mas,lev);
      real dat2=GetDeltaYield(ytag[j],atm,mas,lev);
      fout<<midt.ID(atm,mas,lev)<<"\n";
      fout<<dat<<"\n";
    };
  };
  fout<<"-1\n";
};

void BCGManager::WriteFileYieldDataCBGFormatStandard(string filename,int ynum,string *ytag)
{
  int num=0;
  int sz=nuc.size();
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag())num++;
  };

  string *nuc_nam=new string[num];
  
  int cnt=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()){
      int atm=nuc[i].GetAtomicNumber();
      int mas=nuc[i].GetMassNumber();
      int lev=nuc[i].GetExLevel();
      nuc_nam[cnt++]=midt.Name(atm,mas,lev);
    };
  };

  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  fout.setf(ios::scientific);
  fout.precision(5);
  for(int j=0;j<ynum;j++){
    fout<<" 1\n "<<midt.ID(ytag[j])<<"\n";
    fout<<" "<<num<<"\n";
    for(int k=0;k<num;k++){
      int atm,mas,lev;
      midt.GetParameter(nuc_nam[k],atm,mas,lev);
      real dat=GetYield(ytag[j],atm,mas,lev);
      real dat2=GetDeltaYield(ytag[j],atm,mas,lev);
      fout<<midt.ID(atm,mas,lev)<<"\n";
      fout<<dat<<"\n";
    };
  };
  fout<<"-1\n";

  delete [] nuc_nam;
};

void BCGManager::WriteFileYieldDataCBGFormatCmSF(int num,string *nuc_nam,string filename)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  int ynum=2;
  string ytag[]={"Cm242.SF","Cm244.SF"};
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int j=0;j<ynum;j++){
    if(j==0){
      fout<<" 1\n 962420\n"; // u233
    }else if(j==1){
      fout<<" 1\n 962440\n";
    };
    fout<<" "<<num<<"\n";
    for(int k=0;k<num;k++){
      int atm,mas,lev;
      midt.GetParameter(nuc_nam[k],atm,mas,lev);
      real dat=GetYield(ytag[j],atm,mas,lev);
      real dat2=GetDeltaYield(ytag[j],atm,mas,lev);
      fout<<midt.ID(atm,mas,lev)<<"\n";
      fout<<dat<<"\n";
    };
  };
  fout<<"-1\n";
};

void BCGManager::WriteFileYieldDataCBGFormatDetail(int num,string *nuc_nam,string filename)
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    exit(1);
  };

  int ynum=23;
  string ytag[]={
    "Th232","Pa231","U232","U233","U234",
    "U235","U236","U237","U238","Np237",
    "Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242m","Am243","Cm242","Cm243",
    "Cm244","Cm245","Cm246"
  };
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int j=0;j<ynum;j++){
    if(ytag[j]=="Th232"){
      fout<<" 2\n 902320\n 912330\n"; // for Pa-233
    }else if(ytag[j]=="U235"){
      fout<<" 2\n 922350\n 932360\n"; // for Np-236
    }else if(ytag[j]=="U236"){
      fout<<" 2\n 922360\n 942360\n"; // for Pu-236
    }else if(ytag[j]=="Np237"){
      fout<<" 2\n 932370\n 932390\n"; // for Np-239
    }else if(ytag[j]=="Pu241"){
      fout<<" 2\n 942410\n 952420\n"; // fpr Am-242
    }else{
      fout<<" 1\n "<<midt.ID(ytag[j])<<"\n";
    };
    fout<<" "<<num<<"\n";
    for(int k=0;k<num;k++){
      int atm,mas,lev;
      midt.GetParameter(nuc_nam[k],atm,mas,lev);
      real dat=GetYield(ytag[j],atm,mas,lev);
      real dat2=GetDeltaYield(ytag[j],atm,mas,lev);
      fout<<midt.ID(atm,mas,lev)<<"\n";
      fout<<dat<<"\n";
    };
  };
  fout<<"-1\n";
};

void BCGManager::WriteFileYieldDataCBGFormatDetail2(int num,string *nuc_nam,string filename)
// For JEFF-3.1.1
{
  ofstream fout;
  fout.open(filename.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file\n";
    exit(1);
  };

  int ynum=18;
  string ytag[]={
    "Th232","U233","U234","U235","U236",
    "U238","Np237","Pu238","Pu239","Pu240",
    "Pu241","Pu242","Am241","Am242m","Am243",
    "Cm243","Cm244","Cm245"
  };
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int j=0;j<ynum;j++){
    if(ytag[j]=="Th232"){
      fout<<" 4\n 902320\n 912310\n 912330\n 922320\n"; // for Pa-231, -233, U-232
    }else if(ytag[j]=="U235"){
      fout<<" 2\n 922350\n 932360\n"; // for Np-236
    }else if(ytag[j]=="U236"){
      fout<<" 3\n 922360\n 922370\n 942360\n"; // for U-237, Pu-236
    }else if(ytag[j]=="Np237"){
      fout<<" 2\n 932370\n 932390\n"; // for Np-239
    }else if(ytag[j]=="Pu241"){
      fout<<" 2\n 942410\n 952420\n"; // fpr Am-242
    }else if(ytag[j]=="Pu242"){
      fout<<" 3\n 942420\n 962420\n 962460\n"; // for Cm-242, -246
    }else{
      fout<<" 1\n "<<midt.ID(ytag[j])<<"\n";
    };
    fout<<" "<<num<<"\n";
    for(int k=0;k<num;k++){
      int atm,mas,lev;
      midt.GetParameter(nuc_nam[k],atm,mas,lev);
      real dat=GetYield(ytag[j],atm,mas,lev);
      real dat2=GetDeltaYield(ytag[j],atm,mas,lev);
      fout<<midt.ID(atm,mas,lev)<<"\n";
      fout<<dat<<"\n";
    };
  };
  fout<<"-1\n";
};

void BCGManager::ShowFlagedNuclideList(bool flaged)
{
  int sz=nuc.size();
  int cnt=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()==flaged){
      int am=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      real half_life=nuc[i].GetHalflife();
      string name=Name(am,ms,lv);
      cout<<name<<" "<<half_life<<"\n";
      cnt++;
    };
  };
  cout<<"# Total number of nuclides : "<<cnt<<"\n";
};

void BCGManager::ShowFlagedNuclide(bool flaged,bool name)
{
  int sz=nuc.size();
  int cnt=0;
  cout<<"\n";
  cout<<"    ";
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()==flaged){
      int am=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      if(am>0&&am<100){
	if(name){
  	  cout<<"\""<<midt.Name(am,ms,lv)<<"\",";
	}else{
          int matid=midt.ID(am,ms,lv);
          cout<<matid<<", ";
	};
        cnt++;
        if(cnt%10==0)cout<<"\n    ";
      };
    };
  };
  cout<<"\n# Total number of nuclides : "<<cnt<<"\n\n";
};

void BCGManager::ShowFlagedNuclideShortHalflife(real hl,bool flaged)
{
  cout<<"\n";
  cout<<"# Flaged nuclide with short half-life\n";
  cout<<"# Half life is less than "<<hl<<" [Sec.]\n";
  int sz=nuc.size();
  int cnt=0;
  for(int i=0;i<sz;i++){
    if(nuc[i].Flag()==flaged){
      int am=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      real halflife=nuc[i].GetHalflife();
      if(halflife!=0.&&halflife<hl){
        //int matid=midt.SearchMatIDFromParameter(am,ms,lv);
	//cout<<"#    "<<midt.Name(matid)<<" : "<<halflife<<"\n";
        cout<<"\""<<midt.Name(am,ms,lv)<<" : "<<halflife<<"\n";
      };
    };
  };
};

BCGNuclide& BCGManager::GetNuclide(int atm,int mas,int lev)
{
  int id=GetNuclideIndex(atm,mas,lev);

  if(id==-1){
    cout<<"# Error in BCGManager::GetNuclide\n";
    cout<<"# (atm/mas/lev) = "<<atm<<"/"<<mas<<"/"<<lev<<"\n";
    exit(0);
  };

  return nuc[id];
};

BCGNuclide& BCGManager::GetNuclideFromID(int id)
{
  int lev=id%10;
  id=(id-lev)/10;
  int mas=id%1000;
  id=(id-mas)/1000;
  int atm=id;
  return GetNuclide(atm,mas,lev);
};

real BCGManager::CalTotalDelayedNeutron(string tagname)
{
  int sz=nuc.size();
  real sum=0.;
  for(int i=0;i<sz;i++){
    sum+=nuc[i].GetTotalNumberEmittedNeutrons()*nuc[i].GetYield(tagname);
  };
  //cout<<"#\n# Sum of delayed neutrons : "<<sum<<"\n#\n";
  return sum;
};

real BCGManager::CalTotalDelayedNeutronWithCumulativeYield(string tagname)
{
  int sz=nuc.size();
  real sum=0.;
  for(int i=0;i<sz;i++){
    real ch=nuc[i].GetChannel();
    real tmp=0.;
    for(int j=0;j<ch;j++){
      tmp+=nuc[i].GetBr(j)*nuc[i].GetEmittedNeutron(j);
    };
    sum+=tmp*nuc[i].GetYield(tagname);
  };
  //cout<<"#\n# Sum of delayed neutrons : "<<sum<<"\n#\n";
  return sum;
};

void BCGManager::ShowTotalNumberEmittedNeutrons(string tagname,real low)
{
  int sz=nuc.size();
  cout<<"# Total number of emitted neutrons multiplied by yield\n";
  cout<<"#  (including daughter nuclides-neutron emissions)\n";
  cout<<"#         At  Ms Lv   Yield*N     Yield       N\n";
  cout.setf(ios::showpoint);
  cout<<fixed;
  cout.precision(8);
  real sum=0.;
  for(int i=0;i<sz;i++){
    real n=nuc[i].GetTotalNumberEmittedNeutrons();
    real y=nuc[i].GetYield(tagname);
    if(n*y>low){
      int at=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      string name=midt.Name(at,ms,lv);
      cout<<"# ";
      WriteOut(name,6);
      cout<<" ";
      WriteOut(at,3);
      cout<<" ";
      WriteOut(ms,3);
      cout<<"  ";
      WriteOut(lv,1);
      cout<<" : ";
      cout<<y*n<<"  "<<y<<"  "<<n<<"\n";
    };
    sum+=y*n;
  };
  cout<<"#\n# Sum : "<<sum<<"\n#\n";
};

// +++ Method of creating specific burnup chain

void BCGManager::SearchShortHalflivedNuclide(real hl_upper,real hl_lower,int fpmat_upper,int fpmat_lower)
{
  // In this method, FP nuclides which should be discarded from burnup chain are identified
  // and tagged as 'TRUE'.
  // Long-lived FP, whose half lives are longer than `hl_upper', are NOT excluded from burnup chain.
  // Half-lived FP, whose half lives are shorter than `hl_lower', are excluded from burnup chain
  // since numerical calculation becomes unstable when such FP nuclides are included in a burnup chain.
  // This treatment, therefore, results in inappropriate treatment for
  // decay heat related to the excluded FP nuclides.
  //
  // fpmat_upper=68; // Upper limit of atomic number treated as FP nuclide
  // fpmat_lower=32; // Lower limit of atomic number treated as FP nuclide

  MakeFlagFalseAllNuclide(); 

  int sz=nuc.size();

  bool end_flag=false;
  while(!end_flag){
    end_flag=true;
    for(int i=0;i<sz;i++){
      int am=nuc[i].GetAtomicNumber();
      int ms=nuc[i].GetMassNumber();
      int lv=nuc[i].GetExLevel();
      real half_life=nuc[i].GetHalflife();
      // judgement to take this nuclide from chain
      // ++ half-life check
      bool discard=true;
      if(half_life>hl_upper||half_life==0.)discard=false;
      // ++ decay check
      for(int j=0;j<sz;j++){
  	int dch=nuc[j].GetChannel();
	for(int k=0;k<dch;k++){
          int am2=nuc[j].GetAtomicNumberNext(k);
          int ms2=nuc[j].GetMassNumberNext(k);
          int lv2=nuc[j].GetExLevelNext(k);
	  if(am==am2&&ms==ms2&&lv==lv2&&!nuc[j].Flag())discard=false;  
	};
      };
      if(half_life<hl_lower&&half_life!=0.)discard=true;
      if(am<fpmat_lower||am>fpmat_upper)discard=true;
      /*
      // ++ (n,g) check
      for(int j=0;j<sz;j++){
	int dch=nuc[j].GetReactionChannel(0);
	for(int k=0;k<dch;k++){
          int am2=nuc[j].GetReactionAtomicNumberNext(0,k);
          int ms2=nuc[j].GetReactionMassNumberNext(0,k);
          int lv2=nuc[j].GetReactionExLevelNext(0,k);
	  if(am==am2&&ms==ms2&&lv==lv2&&!nuc[j].Flag())discard=false;  
	};
      };
      */
      if(discard&&!nuc[i].Flag()){
	nuc[i].MakeFlagTrue();
	end_flag=false;
	//cout<<am<<" "<<ms<<" "<<lv<<"\n";
      };
    };      
  };
};

void BCGManager::SearchShortHalflivedNuclideForSpecificNuclides(real hl_upper,real hl_lower,int num_include_SourceTerm,int *AtomicNumber_include_SourceTerm,int fpmat_upper,int fpmat_lower)
{
  int sz=nuc.size();
  MakeFlagTrueAllNuclide();

  //+++(必要な元素のすべての同位体を入れる)+++
  for(int n=0;n<sz;n++){
    int sa=nuc[n].GetAtomicNumber();
    real shl=nuc[n].GetHalflife();
    for(int p=0;p<num_include_SourceTerm;p++){
      if(sa==AtomicNumber_include_SourceTerm[p]&&(shl>hl_lower||shl==0.)){
	nuc[n].MakeFlagFalse();
      };
    };
  };
  
  //+++(上で指定した核種の上流の核種をすべて入れる)+++
  bool end_flag=false;
  while(!end_flag){
    end_flag=true;
    for(int i=0;i<sz;i++){
      int am_p=nuc[i].GetAtomicNumber();
      int ms_p=nuc[i].GetMassNumber();
      int lv_p=nuc[i].GetExLevel();
      int dch_p=nuc[i].GetChannel();
      
      bool include=false;
      
      for(int k=0;k<dch_p;k++){
	int am_d1=nuc[i].GetAtomicNumberNext(k);
	int ms_d1=nuc[i].GetMassNumberNext(k);
	int lv_d1=nuc[i].GetExLevelNext(k);
	for(int j=0;j<sz;j++){
	  int am_d2=nuc[j].GetAtomicNumber();
	  int ms_d2=nuc[j].GetMassNumber();
	  int lv_d2=nuc[j].GetExLevel();
	  if(am_d1==am_d2&&ms_d1==ms_d2&&lv_d1==lv_d2&&!nuc[j].Flag())include=true;  
	};
      };
      if(am_p<fpmat_lower||am_p>fpmat_upper)include=false;
      if(include&&nuc[i].Flag()){
	nuc[i].MakeFlagFalse();
	end_flag=false;
      };
    };
  };

  //+++(上で入れた核種のうち半減期が上方下限以下のものを除外)+++
  for(int ii=0;ii<sz;ii++){
    real hl=nuc[ii].GetHalflife();
    if(hl<hl_upper&&hl!=0.)nuc[ii].MakeFlagTrue();
  };
   
  //+++(上流に除外してない核種がある場合は除外したものを復活させる)+++  
  bool end_flag2=false;
  while(!end_flag2){
    end_flag2=true;
    for(int jj=0;jj<sz;jj++){
      int am3=nuc[jj].GetAtomicNumber();
      int ms3=nuc[jj].GetMassNumber();
      int lv3=nuc[jj].GetExLevel();

      bool include2=false;

      for(int kk=0;kk<sz;kk++){
  	int dch2=nuc[kk].GetChannel();
	for(int nn=0;nn<dch2;nn++){
          int am4=nuc[kk].GetAtomicNumberNext(nn);
          int ms4=nuc[kk].GetMassNumberNext(nn);
          int lv4=nuc[kk].GetExLevelNext(nn);
	  if(am3==am4&&ms3==ms4&&lv3==lv4&&!nuc[kk].Flag())include2=true;
	  for(int mm=0;mm<num_include_SourceTerm;mm++){
	    if(am4==AtomicNumber_include_SourceTerm[mm])include2=false;
	  };
	};
      };
      if(include2&&nuc[jj].Flag()){
	nuc[jj].MakeFlagFalse();
	end_flag2=false;
      };
    };      
  };

  //+++(上方下限以下でソースターム核種でないものは除外する)+++

  for(int jjj=0;jjj<sz;jjj++){
    bool discard=false;
    int am5=nuc[jjj].GetAtomicNumber();
    real hl5=nuc[jjj].GetHalflife();
    if(hl5<hl_upper&&hl5!=0.)nuc[jjj].MakeFlagTrue();
    for(int pp=0;pp<num_include_SourceTerm;pp++){
      if(am5==AtomicNumber_include_SourceTerm[pp])nuc[jjj].MakeFlagFalse();
    };
  };
  
  //+++(下方下限以下のものを除外)+++
  for(int iii=0;iii<sz;iii++){
    real hl_lower=nuc[iii].GetHalflife();
    if(hl_lower<hl_lower&&hl_lower!=0.)nuc[iii].MakeFlagTrue();
  };

};


void BCGManager::ShowYieldDataForXYPlot(string yield_tag)
{
  int zmax=80;
  int nmax=110;
  int zmin=20;
  int nmin=30;

  vector< vector<real> > yld_dat(zmax);
  for(int i=0;i<zmax;i++){
    yld_dat[i].resize(nmax,0.);
  };
  cout<<"# Fissile nuclide : "<<yield_tag<<"\n#\n";

  int sz=nuc.size();
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<sz;i++){
    real yld=nuc[i].GetYield(yield_tag);
    int at=nuc[i].GetAtomicNumber();
    int ms=nuc[i].GetMassNumber();
    int lv=nuc[i].GetExLevel();
    if(lv==0){
      bool cont=true;
      while(cont){
        if(nuc[i+1].GetExLevel()!=0){
          yld+=nuc[i+1].GetYield(yield_tag);
          i++;
        }else{
          cont=false;
        };
      };
    };
    int nn=ms-at;
    if(yld!=0.)yld_dat[at][nn]=yld;
  };

  for(int i=nmin;i<nmax;i++){
    for(int jj=0;jj<2;jj++){
      for(int j=zmin;j<zmax;j++){
        cout<<i-0.5+jj<<" "<<j-0.5<<" "<<yld_dat[j][i]<<"\n";
        cout<<i-0.5+jj<<" "<<j+0.5<<" "<<yld_dat[j][i]<<"\n";
      };
      cout<<"\n";
    };
  };
};

void BCGManager::ShowDecayDataForXYPlot()
{
  int zmax=150;
  int nmax=300;
  int zmin=0;
  int nmin=0;

  vector< vector<real> > dat(150,vector<real>(300,0.));

  int sz=nuc.size();
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<sz;i++){
    real hl=nuc[i].GetHalflife();
    int at=nuc[i].GetAtomicNumber();
    int ms=nuc[i].GetMassNumber();
    int lv=nuc[i].GetExLevel();
    int nn=ms-at;
    real dc=0.;
    if(hl>0.)dc=0.693147/hl;
    if(lv==0)dat[at][nn]=dc;
  };

  for(int i=nmin;i<nmax;i++){
    for(int jj=0;jj<2;jj++){
      for(int j=zmin;j<zmax;j++){
        cout<<i-0.5+jj<<" "<<j-0.5<<" "<<dat[j][i]<<"\n";
        cout<<i-0.5+jj<<" "<<j+0.5<<" "<<dat[j][i]<<"\n";
      };
      cout<<"\n";
    };
  };

};

void BCGManager::ShowDNDataForXYPlot()
{
  CalTotalEmittedNeutrons();

  int zmax=130;
  int nmax=180;
  int zmin=0;
  int nmin=0;

  vector< vector<real> > dat(zmax);
  vector< vector<real> > dat2(zmax); // half-life
  for(int i=0;i<zmax;i++){
    dat[i].resize(nmax,0.);
    dat2[i].resize(nmax,0.);
  };

  int sz=nuc.size();
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<sz;i++){
    int at=nuc[i].GetAtomicNumber();
    int ms=nuc[i].GetMassNumber();
    int lv=nuc[i].GetExLevel();
    real hl=nuc[i].GetHalflife();
    real dnavg=0.;
    int inuc=1;
    for(int j=0;j<nuc[i].GetChannel();j++){
      dnavg+=nuc[i].GetEmittedNeutron(j)*nuc[i].GetBr(j);
    };

    if(lv==0&&i!=sz-1){ // meta-stable nuclides are added to stable nuclides
      bool cont=true;
      while(cont){
        if(nuc[i+1].GetExLevel()!=0){
          for(int j=0;j<nuc[i+1].GetChannel();j++){
            dnavg+=nuc[i+1].GetEmittedNeutron(j)*nuc[i+1].GetBr(j);
          };
          inuc++;
          i++;
        }else{
          cont=false;
        };
      };
    };

    int nn=ms-at;
    real val=dnavg/inuc;
    if(val!=0.){
      if(at<zmin||at>=zmax)cout<<"# Problem in atomic number : "<<at<<"\n";
      if(nn<nmin||nn>=nmax)cout<<"# Problem in neutron number : "<<nn<<"\n";
      dat[at][nn]=val;
      dat2[at][nn]=hl;
    };

  };

  cout<<"#           DN   Half-life [s]\n";
  for(int i=nmin;i<nmax;i++){
    for(int jj=0;jj<2;jj++){
      for(int j=zmin;j<zmax;j++){
        cout<<i-0.5+jj<<" "<<j-0.5<<" "<<dat[j][i]<<" "<<dat2[j][i]<<"\n";
        cout<<i-0.5+jj<<" "<<j+0.5<<" "<<dat[j][i]<<" "<<dat2[j][i]<<"\n";
      };
      cout<<"\n";
    };
  };

  /*
  for(int i=nmin;i<nmax;i++){
    for(int j=zmin;j<zmax;j++){
      if(dat[j][i]!=0.)cout<<j<<" "<<i<<" "<<dat[j][i]<<"\n";
    };
  };
  */

};

// Matsuura
void BCGManager::CalTotalNeutrons()//(cum収率用
{//応急措置的に、とりあえず従来のを上書きしますです。
    // Total number of emitted neutrons is calculated.
    // Neutron emittions of its daughter nuclides are also considered.
    
    //cout<<"#############################################################\n";
    //cout<<"# Total emitted neutron calculation\n#\n";
    
    //int yld_num=yield_tagname.size();
    //vector<real> sumy(yld_num,0.);
    
    int sz=nuc.size();
    for(int i=0;i<sz;i++){
        
        int at_org=nuc[i].GetAtomicNumber();
        int ms_org=nuc[i].GetMassNumber();
        int lv_org=nuc[i].GetExLevel();
        
        bool end=false;
        int iter=0;
        
        real n_emit=0.;
        
        // +++ initialize
        int ngc=nuc[i].GetChannel();
        vector<real> bri(ngc);
        for(int j=0;j<ngc;j++){
            bri[j]=nuc[i].GetBr(j);
            n_emit+=bri[j]*nuc[i].GetEmittedNeutron(j);
        };
        nuc[i].PutTotalNumberEmittedNeutrons(n_emit);
        //cout<<i<<" "<<at_org<<" "<<ms_org<<" "<<lv_org<<" "<<n_emit<<"\n";
        /*
         for(int j=0;j<yld_num;j++){
         sumy[j]+=n_emit*nuc[i].GetYield("tmp");
         };
         */
    };
    //cout<<"# Total number of delayed neutron : "<<sumy[0]<<"\n";
    //cout<<"#############################################################\n\n";
};

void BCGManager::CalTEN()//(cum収率用
{//応急措置的に、とりあえず従来のを上書きしますです。
    // Total number of emitted neutrons is calculated.
    // Neutron emittions of its daughter nuclides are also considered.
    
    //cout<<"#############################################################\n";
    //cout<<"# Total emitted neutron calculation\n#\n";
    
    //int yld_num=yield_tagname.size();
    //vector<real> sumy(yld_num,0.);
    
    int sz=nuc.size();
    for(int i=0;i<sz;i++){
        
        int at_org=nuc[i].GetAtomicNumber();
        int ms_org=nuc[i].GetMassNumber();
        int lv_org=nuc[i].GetExLevel();
        
        bool end=false;
        int iter=0;
        
        real n_emit=0.;
        
        // +++ initialize
        int ngc=nuc[i].GetChannel();
        vector<real> bri(ngc);
        for(int j=0;j<ngc;j++){
            bri[j]=nuc[i].GetBr(j);
            n_emit+=bri[j]*nuc[i].GetEmittedNeutron(j);
        };
        nuc[i].PutTEN(n_emit);
        //cout<<i<<" "<<at_org<<" "<<ms_org<<" "<<lv_org<<" "<<n_emit<<"\n";
        /*
         for(int j=0;j<yld_num;j++){
         sumy[j]+=n_emit*nuc[i].GetYield("tmp");
         };
         */
    };
    //cout<<"# Total number of delayed neutron : "<<sumy[0]<<"\n";
    //cout<<"#############################################################\n\n";
};

void BCGManager::CalRRnums()//(独立収率用
{
    int sz=nuc.size();
    for (int i=0; i<sz; i++) {
        // +++ initialize
        int ngc=nuc[i].GetChannel();
        for (int j=0; j<ngc; j++) {
            int ann = nuc[i].GetAtomicNumberNext(j);
            int mnn = nuc[i].GetMassNumberNext(j);
            int eln = nuc[i].GetExLevelNext(j);
            //cout<<"ann : "<<ann<<" mnn : "<<mnn<<" eln : "<<eln<<endl;
            for (int k=0; k<sz; k++) {
                int atm  = nuc[k].GetAtomicNumber();
                int mass = nuc[k].GetMassNumber();
                int Exl  = nuc[k].GetExLevel();
                //cout<<"atm : "<<atm<<" mass : "<<mass<<" Exl : "<<Exl<<endl;
                if (ann==atm && mnn==mass && eln==Exl) {
                    nuc[i].PutRRnum(j,k);
                    //cout<<k<<endl;
                    break;//
                }
            }
        }
    }//*/
};

void BCGManager::DataPerturbation(string mdir, string filename)
{
  mdir=mdir+filename;
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in BCGManager::DataPerturbation.\n";
    cout<<"# Failed to open the file "<<mdir<<"\n";
    exit(0);
  };

  while(1){
    int mat,mt;
    fin>>mat;
    if(mat==-1)break;
    fin>>mt;
    real chg;
    fin>>chg;
    int iz,ia,il;
    midt.GetParameter(mat,iz,ia,il);
    int ttt=GetNuclideIndex(iz,ia,il);
    if(ttt!=-1){
    if(mt>18000000){
      int fisid=mt-18000000;
      string nucname=midt.Name(fisid);
      real yldorg=GetNuclide(iz,ia,il).GetYield(nucname);
      if(yldorg>=0.){
        real yldnew=yldorg+yldorg*chg;
        GetNuclide(iz,ia,il).PutYield(nucname,yldnew,0.); // 0 is uncertainty
      };
    }else if(mt==8888){
      real hlorg=GetNuclide(iz,ia,il).GetHalflife();
      real hlnew=hlorg+hlorg*chg;
      GetNuclide(iz,ia,il).PutHalflife(hlnew,0.); // 0 is uncertainty
    }else if(mt>=88880&&mt<=88889){
      int ch=mt-88880;
      int chnum=GetNuclide(iz,ia,il).GetChannel();
      if(ch>=chnum){
	cout<<"# Warning : inconsistent decay chennel : "<<mat<<"\n";
	cout<<"#    Channel : "<<chnum<<" / "<<ch<<"\n";
      }else{
        real brorg=GetNuclide(iz,ia,il).GetBr(ch);
        real brnew=brorg+brorg*chg;
        GetNuclide(iz,ia,il).PutBr(ch,brnew);
      };
    }else if(mt>=99990&&mt<=99999){
      int type=mt-99990;
      real deorg=GetNuclide(iz,ia,il).GetDecayEnergy(type);
      real denew=deorg+deorg*chg;
      GetNuclide(iz,ia,il).PutDecayEnergy(type,denew);
    }else{
      cout<<"# Error in BCGManager::DataPerturbation.\n";
      cout<<"# Not coded for MT="<<mt<<"\n";
      exit(0);
    };
    };
  };
  fin.close();
};


