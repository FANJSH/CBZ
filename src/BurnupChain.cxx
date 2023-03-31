#include "BurnupChain.h"

// ++++++++++++++++++++++++++++++
//   NuclideChainData
// ++++++++++++++++++++++++++++++

NuclideChainData::NuclideChainData()
{
  id=-1;
  // +++ fission, capture, n2n, decay
  ndiv.resize(4);
  ratio.resize(4);
  idnext.resize(4);
  decay_energy.resize(3);
};

void NuclideChainData::PutID(int i)
{
  id=i;
};

void NuclideChainData::PutNdiv(int i,int j)
{
  if(i<0||i>3){
    cout<<"# Error in NuclideChainData::PutNdiv.\n";
    exit(0);
  };
  ndiv[i]=j;
  ratio[i].resize(j);
  idnext[i].resize(j);
  if(j==1){
    ratio[i][0]=1.;
    if(i==1)idnext[i][0]=id+10; // (n,g)
    if(i==2)idnext[i][0]=id-10; // (n,2n)
    if(i==3)idnext[i][0]=id-20040; // (alpha-decay)
  };
};

void NuclideChainData::PutNdiv(int *in)
{
  for(int i=0;i<4;i++){
    PutNdiv(i,in[i]);
  };
};

void NuclideChainData::PutNdiv(int i1,int i2,int i3,int i4)
{
  int *in=new int[4];
  in[0]=i1;
  in[1]=i2;
  in[2]=i3;
  in[3]=i4;
  PutNdiv(in);
  delete [] in;
};

void NuclideChainData::PutData(int in,real *rin,int *idin)
{
  for(int i=0;i<ndiv[in];i++){
    ratio[in][i]=rin[i];
    idnext[in][i]=idin[i];
  };
};

void NuclideChainData::PutData(int in,real rin,int idin)
{
  ratio[in][0]=rin;
  idnext[in][0]=idin;
};

void NuclideChainData::PutRatio(int i,int j,real val)
{
  if(j>=ndiv[i]){
    cout<<"# Error in NuclideChainData::PutRatio.\n";
    cout<<"# The number of channels is "<<ndiv[i]<<"\n";
    cout<<"# You are choosing a channel "<<j<<"\n";
    exit(0);
  };
  ratio[i][j]=val;
};

void NuclideChainData::ShowFissionYieldData(MATIDTranslator &midt)
{
  cout.setf(ios::scientific);
  cout.precision(5);
  for(int i=0;i<ndiv[0];i++){
    cout<<midt.Name(idnext[0][i])<<" "<<idnext[0][i]<<" "<<ratio[0][i]<<"\n";
  };
};

void NuclideChainData::ShowSelf(MATIDTranslator &midt)
{
  cout.setf(ios::showpoint);
  cout.precision(10);
  cout<<midt.Name(id)<<"\n";
  for(int i=0;i<4;i++){
    if(i==0)cout<<" Fission \n";
    if(i==1)cout<<" Capture\n";
    if(i==2)cout<<" (n,2n)\n";
    if(i==3)cout<<" Decay\n";
    cout<<"\n";
    cout<<"    Number of channels : "<<ndiv[i]<<"\n";
    real sum=0.;
    for(int j=0;j<ndiv[i];j++){
      cout<<"         "<<j<<" : "<<midt.Name(idnext[i][j])<<" "<<ratio[i][j]<<"\n";
      sum+=ratio[i][j];
    };
    //cout<<"# sum : "<<sum<<"\n";
    cout<<"\n";
  };
};

void NuclideChainData::SetPseudoFP(int fpid)
{
  real sum=0.;
  for(int i=0;i<ndiv[0];i++){
    sum+=ratio[0][i];
  };

  ndiv[0]+=1;
  ratio[0].push_back(2.-sum);
  idnext[0].push_back(fpid);
};

void NuclideChainData::SetSpontaneousFission(real br)
{
  int fp_num=GetNdiv(0);
  
  int decay_ndiv=GetNdiv(3);
  real factor=1.-br;
  for(int i=0;i<decay_ndiv;i++){
    ratio[3][i]*=factor;
  };

  ndiv[3]+=fp_num;
  for(int i=0;i<fp_num;i++){
    int fpid=idnext[0][i];
    real fpy=ratio[0][i];
    idnext[3].push_back(fpid);
    ratio[3].push_back(fpy*br);
  };
};

void NuclideChainData::SetSpontaneousFission(real br, NuclideChainData &ncd2)
{
  int fp_num=ncd2.GetNdiv(0);
  
  int decay_ndiv=GetNdiv(3);
  real factor=1.-br;
  for(int i=0;i<decay_ndiv;i++){
    ratio[3][i]*=factor;
  };

  ndiv[3]+=fp_num;
  for(int i=0;i<fp_num;i++){
    int fpid=ncd2.GetIDnext(0,i);
    real fpy=ncd2.GetRatio(0,i);
    idnext[3].push_back(fpid);
    ratio[3].push_back(fpy*br);
  };

};

// ++++++++++++++++++++++++++
//    BurnupChain
// ++++++++++++++++++++++++++

BurnupChain::BurnupChain()
{
  read_decay_constant=false;
};

void BurnupChain::PutNdiv(string nuc,int i1,int i2,int i3,int i4)
{
  int id=ID(nuc);
  data[id].PutNdiv(i1,i2,i3,i4);
};

void BurnupChain::PutIDnextDecay(string nuc,string nuc2)
{
  int id=ID(nuc);
  int id2=ID(nuc2);
  data[id].PutIDnextDecay(id2);
};

void BurnupChain::PutIDnextCapture(string nuc,string nuc2)
{
  int id=ID(nuc);
  int id2=ID(nuc2);
  data[id].PutIDnextCapture(id2);
};

void BurnupChain::PutIDnextN2N(string nuc,string nuc2)
{
  int id=ID(nuc);
  int id2=ID(nuc2);
  data[id].PutIDnextN2N(id2);
};

void BurnupChain::PutData(string nuc,int in,real *rin,string *nuc2)
{
  int id=ID(nuc);
  int ndiv=data[id].GetNdiv(in);
  int *idin=new int[ndiv];
  for(int i=0;i<ndiv;i++){
    idin[i]=ID(nuc2[i]);
  };
  data[id].PutData(in,rin,idin);
  delete [] idin;
};

void BurnupChain::PutData(int id)
{
  NuclideChainData dummy;
  data.insert(map<int,NuclideChainData>::value_type(id,dummy));
};


void BurnupChain::AddPseudoFP()
{
  int fpid_u5=9950000;
  int fpid_u8=9980000;
  int fpid_p1=9910000;
  int fpid_p9=9990000;

  int nucnum=1+2+7+3+6+4+9+1; // 33 Heavy nuclides
  string nucname[]={
    "Th232",
    "Pa231","Pa233",
    "U232","U233","U234","U235","U236","U237","U238",
    "Np236","Np237","Np239",
    "Pu236","Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242","Am242m","Am243",
    "Cm242","Cm243","Cm244","Cm245","Cm246","Cm247","Cm248","Cm249","Cm250",
    "Bk249",
  };
  int pfp[]={
    fpid_u5,
    fpid_u5,fpid_u5,
    fpid_u5,fpid_u5,fpid_u5,fpid_u5,fpid_u8,fpid_u8,fpid_u8,
    fpid_u5,fpid_u8,fpid_u8,
    fpid_u5,fpid_u8,fpid_p9,fpid_p9,fpid_p1,fpid_p1,
    fpid_p1,fpid_p1,fpid_p1,fpid_p1,
    fpid_p1,fpid_p1,fpid_p1,fpid_p1,fpid_p1,fpid_p1,fpid_p1,fpid_p1,fpid_p1,
    fpid_p1,
  };

  for(int i=0;i<nucnum;i++){
    int id=ID(nucname[i]);
    if(data[id].GetID()!=0){
      data[id].PutNdivFission(1);
      data[id].PutIDnextFission(0,pfp[i],1.);
    };
  };

  // fp-u5
  data[fpid_u5].PutID(fpid_u5);
  data[fpid_u5].PutNdiv(0,1,1,1); // f,c,2n,d
  data[fpid_u5].PutData(1,1.,fpid_u5);
  data[fpid_u5].PutData(2,1.,fpid_u5);
  data[fpid_u5].PutData(3,1.,fpid_u5);
  // fp-u8
  data[fpid_u8].PutID(fpid_u8);
  data[fpid_u8].PutNdiv(0,1,1,1); // f,c,2n,d
  data[fpid_u8].PutData(1,1.,fpid_u8);
  data[fpid_u8].PutData(2,1.,fpid_u8);
  data[fpid_u8].PutData(3,1.,fpid_u8);
  // fp-p9
  data[fpid_p9].PutID(fpid_p9);
  data[fpid_p9].PutNdiv(0,1,1,1); // f,c,2n,d
  data[fpid_p9].PutData(1,1.,fpid_p9);
  data[fpid_p9].PutData(2,1.,fpid_p9);
  data[fpid_p9].PutData(3,1.,fpid_p9);
  // fp-p1
  data[fpid_p1].PutID(fpid_p1);
  data[fpid_p1].PutNdiv(0,1,1,1); // f,c,2n,d
  data[fpid_p1].PutData(1,1.,fpid_p1);
  data[fpid_p1].PutData(2,1.,fpid_p1);
  data[fpid_p1].PutData(3,1.,fpid_p1);

};

void BurnupChain::AddPseudoFPMolybdenum()
{
  int mo_index=420000;

  int nucnum=1+2+7+3+6+4+9+1; // 33 Heavy nuclides
  string nucname[]={
    "Th232",
    "Pa231","Pa233",
    "U232","U233","U234","U235","U236","U237","U238",
    "Np236","Np237","Np239",
    "Pu236","Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242","Am242m","Am243",
    "Cm242","Cm243","Cm244","Cm245","Cm246","Cm247","Cm248","Cm249","Cm250",
    "Bk249",
  };

  for(int i=0;i<nucnum;i++){
    int id=ID(nucname[i]);
    if(data[id].GetID()!=0){
      data[id].PutNdivFission(1); // The number of paths of nuclide change through fission reaction
      data[id].PutIDnextFission(0,mo_index,2.); // Two natural molybdenum nuclides are generated by one fission reaction
    };
  };

  // Mo-100
  data[mo_index].PutID(mo_index);
  data[mo_index].PutNdiv(0,1,1,0); // f,c,2n,d
  data[mo_index].PutData(1,1.,mo_index); // Molybdenum is generated by capture reaction of Mo.
  data[mo_index].PutData(2,1.,mo_index); // Molybdenum is generated by (n,2n) reaction of Mo.
};

void BurnupChain::AddPseudoFP(int fissile_num, int *matid_fp, real *yield_fp)
{
  for(int i=0;i<fissile_num;i++){
    int id=matid_fp[i]-99000000;
    if(data[id].GetID()!=0){
      data[id].AddFissionYield(matid_fp[i], yield_fp[i]);

      data[matid_fp[i]].PutID(matid_fp[i]);
      data[matid_fp[i]].PutNdiv(0,1,1,0); // f,c,2n,d
      data[matid_fp[i]].PutData(1,1.,matid_fp[i]); // The same nuclide is generated by capture reaction. 
      data[matid_fp[i]].PutData(2,1.,matid_fp[i]); // The same nuclide is generated by (n,2n) reaction.
    
    };
  };
};

void BurnupChain::SetDefault()
{
  //
  // SPECIFIC FOR FAST REACTOR ANALYSES
  //

  // u-235
  data[ID("U235")].PutID(ID("U235"));
  data[ID("U235")].PutNdiv(0,1,0,0); // f,c,2n,d
  // u-236
  data[ID("U236")].PutID(ID("U236"));
  data[ID("U236")].PutNdiv(0,1,1,0);
  data[ID("U236")].PutData(1,1.,ID("U238"));
  // u-238
  data[ID("U238")].PutID(ID("U238"));
  data[ID("U238")].PutNdiv(0,1,1,0);
  data[ID("U238")].PutData(1,1.,ID("Pu239"));
  data[ID("U238")].PutData(2,1.,ID("Np237"));
  // np-237
  data[ID("Np237")].PutID(ID("Np237"));
  data[ID("Np237")].PutNdiv(0,1,0,0);
  data[ID("Np237")].PutData(1,1.,ID("Pu238"));
  // pu-238
  data[ID("Pu238")].PutID(ID("Pu238"));
  data[ID("Pu238")].PutNdiv(0,1,0,0);
  // pu-239
  data[ID("Pu239")].PutID(ID("Pu239"));
  data[ID("Pu239")].PutNdiv(0,1,1,1);
  // pu-240
  data[ID("Pu240")].PutID(ID("Pu240"));
  data[ID("Pu240")].PutNdiv(0,1,1,1);
  // pu-241
  data[ID("Pu241")].PutID(ID("Pu241"));
  data[ID("Pu241")].PutNdiv(0,1,1,1);
  data[ID("Pu241")].PutData(3,1.,ID("Am241")); // beta-decay
  // pu-242
  data[ID("Pu242")].PutID(ID("Pu242"));
  data[ID("Pu242")].PutNdiv(0,1,1,1);
  data[ID("Pu242")].PutData(1,1.,ID("Am243"));
  // am-241
  data[ID("Am241")].PutID(ID("Am241"));
  data[ID("Am241")].PutNdiv(0,3,0,1);
  real rin0[]={0.1494, 0.8506*0.827, 0.8506*0.173};
  int din0[]={ID("Am242m"), ID("Cm242"), ID("Pu242")};
  data[ID("Am241")].PutData(1,rin0,din0);
  data[ID("Am241")].PutData(3,1.,ID("Np237"));
  // am-242m
  data[ID("Am242m")].PutID(ID("Am242m"));
  data[ID("Am242m")].PutNdiv(0,1,1,1);
  data[ID("Am242m")].PutData(1,1.,ID("Am243"));
  data[ID("Am242m")].PutData(2,1.,ID("Am241"));
  data[ID("Am242m")].PutData(3,1.,ID("Cm242"));
  // am-243
  data[ID("Am243")].PutID(ID("Am243"));
  data[ID("Am243")].PutNdiv(0,1,1,1);
  data[ID("Am243")].PutData(1,1.,ID("Cm244"));
  data[ID("Am243")].PutData(2,1.,ID("Am242m"));
  data[ID("Am243")].PutData(3,1.,ID("Pu239"));
  // cm-242
  data[ID("Cm242")].PutID(ID("Cm242"));
  data[ID("Cm242")].PutNdiv(0,1,0,1);
  data[ID("Cm242")].PutData(3,1.,ID("Pu238"));
  // cm-243
  data[ID("Cm243")].PutID(ID("Cm243"));
  data[ID("Cm243")].PutNdiv(0,1,1,1);
  data[ID("Cm243")].PutData(3,1.,ID("Pu239"));
  // cm-244
  data[ID("Cm244")].PutID(ID("Cm244"));
  data[ID("Cm244")].PutNdiv(0,1,1,1);
  data[ID("Cm244")].PutData(3,1.,ID("Pu240"));
  // cm-245
  data[ID("Cm245")].PutID(ID("Cm245"));
  data[ID("Cm245")].PutNdiv(0,1,1,1);
  data[ID("Cm245")].PutData(3,1.,ID("Pu241"));
  // cm-246
  data[ID("Cm246")].PutID(ID("Cm246"));
  data[ID("Cm246")].PutNdiv(0,0,1,1);
  data[ID("Cm246")].PutData(3,1.,ID("Pu242"));

  // Fission reaction chain is added.
  AddPseudoFP();
};


void BurnupChain::Set28HeavyMetalChain(bool fr)
{
  // You can get detailed information from my notebook in 2012/9/13 and 2014/8/7
  //
  // - (n,g) and (n,2n) branching ratio is calculated from LWR burnup calculation results and JEFF-3.1A data in 2012/9/13.
  //
  // - Am-241(n,g) and Np-237(n,2n) branching ratio is re-calculated using JENDL-4.0 in 2014/8/7.
  //
  // - Decay branching ratio of Np-236m is based on ENDF/B-VII.1 evaluation.
  //   JENDL/DD-2015 gives the same evaluation (0.5 & 0.5)
  //
  // - (n,2n) branching ratio of Np-237 is re-calculated with 361-group library by [Burner/main.br.lwr.cxx].
  //   [0.3673/0.6327] is changed to [0.3351/0.6649] in 2022/8/12]
  //   Fast reactor case is unchanged.
  
  Set21HeavyMetalChain(fr);

  int add_num=7;
  string add_name[]={"Th232","Pa231","Pa233","U232","U233","Np236","Pu236"};

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  // 
  data[ID("Th232")].PutNdiv(0,1,1,0); // f,c,2n,d
  PutIDnextCapture("Th232","Pa233");
  PutIDnextN2N("Th232","Pa231");
  //
  data[ID("Pa231")].PutNdiv(0,1,0,0); // f,c,2n,d
  PutIDnextCapture("Pa231","U232");
  //
  data[ID("Pa233")].PutNdiv(0,1,1,1); // f,c,2n,d
  PutIDnextCapture("Pa233","U234");
  PutIDnextN2N("Pa233","U232");
  PutIDnextDecay("Pa233","U233");
  //
  data[ID("U232")].PutNdiv(0,1,0,0); // f,c,2n,d
  //
  data[ID("U233")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Np236")].PutNdiv(0,1,0,0); // f,c,2n,d
  //
  data[ID("Pu236")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Pu236","U232");
  //
  // (modification for already existing 21 HM)
  data[ID("U234")].PutNdiv(2,1);
  //
  data[ID("U236")].PutNdiv(3,1); // decay to Th-232
  //
  data[ID("Np237")].PutNdiv(2,3);
  int id_n7_n2n[]={ID("Np236"),ID("U236"),ID("Pu236")};
  //real r_n7_n2n[]={0.3673, 0.6327*0.5, 0.6327*0.5};
  real r_n7_n2n[]={0.3351, 0.6649*0.5, 0.6649*0.5};  // switched in 2022/8/12.
 
  if(fr){ // For fast reactor
    r_n7_n2n[0]=0.3631;
    r_n7_n2n[1]=0.6369*0.5;
    r_n7_n2n[2]=0.6369*0.5;
  };
  data[ID("Np237")].PutData(2,r_n7_n2n,id_n7_n2n);
};

void BurnupChain::Set37HeavyMetalChain(bool fr)
{
  Set28HeavyMetalChain(fr);

  int add_num=9;
  string add_name[]={"Th228","Ra224","Rn220","Po216","Pb212","Bi212","Tl208","Po212","Pb208"};

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  //
  data[ID("U232")].PutNdiv(0,1,0,1); // f,c,2n,d
  PutIDnextDecay("U232","Th228");
  //
  data[ID("Th228")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Th228","Ra224");
  //
  data[ID("Ra224")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ra224","Rn220");
  //
  data[ID("Rn220")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Rn220","Po216");
  //
  data[ID("Po216")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Po216","Pb212");
  //
  data[ID("Pb212")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Pb212","Bi212");
  //
  data[ID("Bi212")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d[]={0.6406, 0.3594}; // from ENDF/B-VII.1
  int id_d[]={ID("Po212"),ID("Tl208")};
  data[ID("Bi212")].PutData(3,r_d,id_d);
  //
  data[ID("Po212")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Po212","Pb208");
  //
  data[ID("Tl208")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Tl208","Pb208");
};

void BurnupChain::AddThoriumDecaySeries()
{
  if(data[ID("U232")].GetID()!=922320){
    cout<<"# Error in BurnupChain::AddThoriumDecaySeries.\n";
    cout<<"# Original chain does NOT contain uranium-232,\n";
    cout<<"# which is source nuclide of this (partial) decay series.\n";
    exit(0);
  };

  int add_num=9;
  string add_name[]={"Th228","Ra224","Rn220","Po216","Pb212","Bi212","Tl208","Po212","Pb208"};

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  //
  data[ID("U232")].PutNdiv(0,1,0,1); // f,c,2n,d
  PutIDnextDecay("U232","Th228");
  //
  data[ID("Th228")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Th228","Ra224");
  //
  data[ID("Ra224")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ra224","Rn220");
  //
  data[ID("Rn220")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Rn220","Po216");
  //
  data[ID("Po216")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Po216","Pb212");
  //
  data[ID("Pb212")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Pb212","Bi212");
  //
  data[ID("Bi212")].PutNdiv(0,0,0,3); // f,c,2n,d
  real r_d[]={0.64056, 0.35930, 0.00014}; // from JEFF-3.1.1
  int id_d[]={ID("Po212"),ID("Tl208"),ID("Pb208")};
  data[ID("Bi212")].PutData(3,r_d,id_d);
  //
  data[ID("Po212")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Po212","Pb208");
  //
  data[ID("Tl208")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Tl208","Pb208");
};

void BurnupChain::AddNeptuniumDecaySeries()
{
  if(data[ID("U233")].GetID()!=922330){
    cout<<"# Error in BurnupChain::AddNeptuniumDecaySeries.\n";
    cout<<"# Original chain does NOT contain uranium-233,\n";
    cout<<"# which is source nuclide of this (partial) decay series.\n";
    exit(0);
  };

  int add_num=11;
  string add_name[]={
   "Th229","Ra225","Ac225","Fr221","At217",
   "Rn217","Bi213","Po213","Tl209","Pb209",
   "Bi209"
  };

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  //
  data[ID("U233")].PutNdiv(0,1,1,1); // f,c,2n,d
  PutIDnextDecay("U233","Th229");
  //
  data[ID("Th229")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Th229","Ra225");
  //
  data[ID("Ra225")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ra225","Ac225");
  //
  data[ID("Ac225")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ac225","Fr221");
  //
  data[ID("Fr221")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Fr221","At217");
  //
  data[ID("At217")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d[]={0.99988, 0.00012}; // from JEFF-3.1.1
  int id_d[]={ID("Bi213"),ID("Rn217")};
  data[ID("At217")].PutData(3,r_d,id_d);
  //
  data[ID("Bi213")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d2[]={0.9784, 0.0216}; // from JEFF-3.1.1
  int id_d2[]={ID("Po213"),ID("Tl209")};
  data[ID("Bi213")].PutData(3,r_d2,id_d2);
  //
  data[ID("Rn217")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Rn217","Po213");
  //
  data[ID("Po213")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Po213","Pb209");
  //
  data[ID("Tl209")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Tl209","Pb209");
  //
  data[ID("Pb209")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Pb209","Bi209");
};

void BurnupChain::AddUraniumDecaySeries()
{
  if(data[ID("U234")].GetID()!=922340){
    cout<<"# Error in BurnupChain::AddUraniumDecaySeries.\n";
    cout<<"# Original chain does NOT contain uranium-234,\n";
    cout<<"# which is source nuclide of this (partial) decay series.\n";
    exit(0);
  };

  int add_num=3;
  string add_name[]={
   "Th230","Ra226","Rn222"
  };

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  //
  data[ID("U234")].PutNdiv(0,1,1,1); // f,c,2n,d
  PutIDnextDecay("U234","Th230");
  //
  data[ID("Th230")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Th230","Ra226");
  //
  data[ID("Ra226")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ra226","Rn222");
};

void BurnupChain::AddActiniumDecaySeries()
{
  if(data[ID("Pa231")].GetID()!=912310){
    cout<<"# Error in BurnupChain::AddActiniumDecaySeries.\n";
    cout<<"# Original chain does NOT contain protoactinium-231,\n";
    cout<<"# which is source nuclide of this (partial) decay series.\n";
    exit(0);
  };

  int add_num=8;
  string add_name[]={
    "Th231","Ac227","Th227","Ra223","Fr223",
    "At219","Rn219","Bi215"
  };

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  //
  data[ID("U235")].PutNdiv(0,1,1,1); // f,c,2n,d
  PutIDnextDecay("U235","Th231");
  //
  data[ID("Th231")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Th231","Pa231");
  //
  data[ID("Pa231")].PutNdiv(0,1,0,1); // f,c,2n,d
  PutIDnextCapture("Pa231","U232");
  PutIDnextDecay("Pa231","Ac227");
  //
  data[ID("Ac227")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d[]={0.9862, 0.0138}; // from JEFF-3.1.1
  int id_d[]={ID("Th227"),ID("Fr223")}; 
  data[ID("Ac227")].PutData(3,r_d,id_d);
  //
  data[ID("Th227")].PutNdiv(0,0,0,2); // f,c,2n,d
  PutIDnextDecay("Th227","Ra223");
  //
  data[ID("Fr223")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d2[]={0.99994, 0.00006}; // from JEFF-3.1.1
  int id_d2[]={ID("Ra223"),ID("At219")};
  data[ID("Fr223")].PutData(3,r_d2,id_d2);
  //
  data[ID("Ra223")].PutNdiv(0,0,0,1); // f,c,2n,d
  PutIDnextDecay("Ra223","Rn219");
  //
  data[ID("At219")].PutNdiv(0,0,0,2); // f,c,2n,d
  real r_d3[]={0.03, 0.97}; // from JEFF-3.1.1
  int id_d3[]={ID("Rn219"),ID("Bi215")};
  data[ID("At219")].PutData(3,r_d3,id_d3);
};


void BurnupChain::AddCobalt59Activation()
{
  data[270590].PutID(270590);
  data[270590].PutNdiv(0,1,0,0); // f,c,2n,d

  data[270600].PutID(270600);
};


void BurnupChain::SetHeavyMetalChainForADS()
{
  Set21HeavyMetalChain(true);

  int add_num=5;
  string add_name[]={"Cm247","Cm248","Cm249","Cm250","Bk249"};

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  // 
  data[ID("Cm247")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Cm248")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Cm249")].PutNdiv(0,1,1,1); // f,c,2n,d
  data[ID("Cm249")].PutIDnextDecay(ID("Bk249"));
  //
  data[ID("Cm250")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Bk249")].PutNdiv(0,1,1,0); // f,c,2n,d
  
  // (modification for already existing 21 HM)
  data[ID("Cm246")].PutNdiv(0,1,1,1);
  data[ID("Cm246")].PutIDnextDecay(ID("Pu242"));
};

void BurnupChain::AddHigherCmPlusBk()
{
  int add_num=5;
  string add_name[]={"Cm247","Cm248","Cm249","Cm250","Bk249"};

  for(int i=0;i<add_num;i++){
    int id=ID(add_name[i]);
    data[id].PutID(id);
  };

  // 
  data[ID("Cm247")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Cm248")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Cm249")].PutNdiv(0,1,1,1); // f,c,2n,d
  data[ID("Cm249")].PutIDnextDecay(ID("Bk249"));
  //
  data[ID("Cm250")].PutNdiv(0,1,1,0); // f,c,2n,d
  //
  data[ID("Bk249")].PutNdiv(0,1,1,0); // f,c,2n,d
  
  // (modification for already existing HM chain)
  data[ID("Cm246")].PutNdiv(0,1,1,1);
  data[ID("Cm246")].PutIDnextDecay(ID("Pu242")); // re-definition
};

void BurnupChain::Set21HeavyMetalChain(bool fr)
{
  int hm_num=5+2+5+4+5; 
  string hm_name[]={
   "U234","U235","U236","U237","U238",
   "Np237","Np239",
   "Pu238","Pu239","Pu240","Pu241","Pu242",
   "Am241","Am242","Am242m","Am243",
   "Cm242","Cm243","Cm244","Cm245","Cm246",
  };
  for(int i=0;i<hm_num;i++){
    int id=ID(hm_name[i]);
    data[id].PutID(id);
  };

  // u-234
  data[ID("U234")].PutNdiv(0,1,0,0);
  // u-235
  data[ID("U235")].PutNdiv(0,1,1,0); // f,c,2n,d
  // u-236
  data[ID("U236")].PutNdiv(0,1,1,0);
  // u-237
  data[ID("U237")].PutNdiv(0,1,1,1);
  data[ID("U237")].PutIDnextDecay(ID("Np237")); // beta-
  // u-238
  data[ID("U238")].PutNdiv(0,1,1,0);
  data[ID("U238")].PutIDnextCapture(ID("Np239"));
  // np-237
  data[ID("Np237")].PutNdiv(0,1,1,0);
  data[ID("Np237")].PutIDnextCapture(ID("Pu238"));
  data[ID("Np237")].PutIDnextN2N(ID("U237"));
  // np-239
  data[ID("Np239")].PutNdiv(0,1,0,1);
  data[ID("Np239")].PutIDnextCapture(ID("Pu240"));
  data[ID("Np239")].PutIDnextDecay(ID("Pu239"));
  // pu-238
  data[ID("Pu238")].PutNdiv(0,1,0,1);
  data[ID("Pu238")].PutIDnextDecay(ID("U234"));
  // pu-239
  data[ID("Pu239")].PutNdiv(0,1,1,1);
  data[ID("Pu239")].PutIDnextDecay(ID("U235"));
  // pu-240
  data[ID("Pu240")].PutNdiv(0,1,1,1);
  data[ID("Pu240")].PutIDnextDecay(ID("U236"));
  // pu-241
  data[ID("Pu241")].PutNdiv(0,1,1,2);
  real r_p1_d[]={0.9999755, 0.0000245}; // ENDF/B-VII.1 decay data
  int id_p1_d[]={ID("Am241"),ID("U237")};
  data[ID("Pu241")].PutData(3,r_p1_d,id_p1_d);
  // pu-242
  data[ID("Pu242")].PutNdiv(0,1,1,1);
  data[ID("Pu242")].PutIDnextCapture(ID("Am243"));
  data[ID("Pu242")].PutIDnextDecay(ID("U238"));
  // am-241
  data[ID("Am241")].PutNdiv(0,2,0,1);
  //real r_a1_c[]={0.88148,0.11852}; // from SRAC chain (based on ENDF/B-VI)
  real r_a1_c[]={0.877983, 0.122017}; // from ChainJ40. Replaced in 2021/7/16
  int id_a1_c[]={ID("Am242"),ID("Am242m") };
  if(fr){ // For fast reactor
    r_a1_c[0]=0.8506;
    r_a1_c[1]=0.1494;
  };
  data[ID("Am241")].PutData(1,r_a1_c,id_a1_c);
  data[ID("Am241")].PutIDnextDecay(ID("Np237"));
  // am-242
  data[ID("Am242")].PutNdiv(0,1,0,2);
  real r_a2_d[]={0.173, 0.827}; // ENDF/B-VII.1 decay data 
  int id_a2_d[]={ID("Pu242"),ID("Cm242")};
  data[ID("Am242")].PutData(3,r_a2_d,id_a2_d);
  // am-242m
  data[ID("Am242m")].PutNdiv(0,1,1,1);
  data[ID("Am242m")].PutIDnextDecay(ID("Am242"));
  data[ID("Am242m")].PutIDnextCapture(ID("Am243"));
  data[ID("Am242m")].PutIDnextN2N(ID("Am241"));
  // am-243
  data[ID("Am243")].PutNdiv(0,1,1,1);
  data[ID("Am243")].PutIDnextDecay(ID("Np239"));
  data[ID("Am243")].PutIDnextCapture(ID("Cm244"));
  // cm-242
  data[ID("Cm242")].PutNdiv(0,1,0,1);
  data[ID("Cm242")].PutIDnextDecay(ID("Pu238"));
  // cm-243
  data[ID("Cm243")].PutNdiv(0,1,1,1);
  data[ID("Cm243")].PutIDnextDecay(ID("Pu239"));
  // cm-244
  data[ID("Cm244")].PutNdiv(0,1,1,1);
  data[ID("Cm244")].PutIDnextDecay(ID("Pu240"));
  // cm-245
  data[ID("Cm245")].PutNdiv(0,1,1,1);
  data[ID("Cm245")].PutIDnextDecay(ID("Pu241"));
  // cm-246
  data[ID("Cm246")].PutNdiv(0,1,1,1);
  data[ID("Cm246")].PutIDnextDecay(ID("Pu242"));
};

void BurnupChain::SetOKUMURAChainFine(bool fr)
{
  // SRAC : fp193bp6 chain
  string fp_name[]={
   "Ge073","Ge074","Ge076",
   "As075",
   "Se076","Se077","Se078","Se079","Se080","Se082",
   "Br081",
   "Kr082","Kr083","Kr084","Kr085","Kr086",
   "Rb085","Rb086","Rb087",
   "Sr086","Sr087","Sr088","Sr089","Sr090",
   "Y089","Y090","Y091",
   "Zr090","Zr091","Zr092","Zr093","Zr094","Zr095","Zr096",
   "Nb093","Nb093m","Nb094","Nb095",
   "Mo092","Mo094","Mo095","Mo096","Mo097","Mo098","Mo099","Mo100",
   "Tc099",
   "Ru100","Ru101","Ru102","Ru103","Ru104","Ru105","Ru106",
   "Rh103","Rh105","Rh106",
   "Pd104","Pd105","Pd106","Pd107","Pd108","Pd110",
   "Ag107","Ag109","Ag110m",
   "Cd110","Cd111","Cd112","Cd113","Cd113m","Cd114","Cd116",
   "In113","In115",
   "Sn116","Sn117","Sn118","Sn119","Sn119m","Sn120","Sn121","Sn121m","Sn122","Sn123","Sn124","Sn126",
   "Sb121","Sb123","Sb124","Sb125","Sb126","Sb126m",
   "Te122","Te123","Te123m","Te124","Te125","Te125m","Te126","Te127m","Te128","Te129m","Te130","Te132",
   "I127","I129","I130","I131","I135",
   "Xe126","Xe128","Xe129","Xe130","Xe131","Xe132","Xe133","Xe134","Xe135","Xe136",
   "Cs133","Cs134","Cs135","Cs136","Cs137",
   "Ba134","Ba135","Ba136","Ba137","Ba137m","Ba138","Ba140",
   "La139","La140",
   "Ce140","Ce141","Ce142","Ce143","Ce144",
   "Pr141","Pr143","Pr144",
   "Nd142","Nd143","Nd144","Nd145","Nd146","Nd147","Nd148","Nd150",
   "Pm147","Pm148","Pm148m","Pm149","Pm151", 
   "Sm147","Sm148","Sm149","Sm150","Sm151","Sm152","Sm153","Sm154",
   "Eu151","Eu152","Eu153","Eu154","Eu155","Eu156","Eu157",
   "Gd152","Gd154","Gd155","Gd156","Gd157","Gd158","Gd160",
   "Tb159","Tb160",
   "Dy160","Dy161","Dy162","Dy163","Dy164",
   "Ho163","Ho165","Ho166m",
   "Er162","Er164","Er166","Er167","Er168","Er170",
   "Hf176","Hf177","Hf178","Hf179","Hf180",
  };

  int ii=0;
  bool flagend=false;
  while(!flagend){
    if(fp_name[ii]=="Hf180"){
      flagend=true;
    };
    ii++;
  };
  int fp_num=ii;

  for(int i=0;i<fp_num;i++){
    int id=ID(fp_name[i]);
    data[id].PutID(id);
  };

  Set21HeavyMetalChain(fr);

  // +++ FP
  // Ge73
  PutNdiv("Ge073",0,1,0,0); // f,c,2n,d
  // Ge74
  PutNdiv("Ge074",0,1,0,0);
  PutIDnextCapture("Ge074","As075");
  // As75
  // (In SRAC, path to Ge076 is neglected.)
  PutNdiv("As075",0,2,0,0);
  real r_as5_c[]={0.9998,0.0002};
  string nuc_as5_c[]={"Se076","Ge076"};
  PutData("As075",1,r_as5_c,nuc_as5_c);
  // Se76
  PutNdiv("Se076",0,1,0,0);
  // Se77
  PutNdiv("Se077",0,1,0,0);
  // Se78
  PutNdiv("Se078",0,1,0,0);
  // Se79
  PutNdiv("Se079",0,1,0,0);
  // Se80
  PutNdiv("Se080",0,1,0,0);
  PutIDnextCapture("Se080","Br081");
  // Br81
  data[ID("Br081")].PutNdiv(0,1,0,0);
  data[ID("Br081")].PutIDnextCapture(ID("Kr082"));
  // Kr82
  data[ID("Kr082")].PutNdiv(0,1,0,0);
  // Kr83
  data[ID("Kr083")].PutNdiv(0,1,0,0);
  // Kr84
  data[ID("Kr084")].PutNdiv(0,2,0,0);
  real r_kr84_c[]={ 0.46348, 0.53652};
  //real r_kr84_c[]={0.9455, 0.043}; // FAST
  string nuc_kr84_c[]={"Kr085","Rb085"};
  PutData("Kr084",1,r_kr84_c,nuc_kr84_c);
  // Kr85 
  data[ID("Kr085")].PutNdiv(0,1,0,1);
  data[ID("Kr085")].PutIDnextDecay(ID("Rb085"));
  // Kr86
  data[ID("Kr086")].PutNdiv(0,1,0,0);
  data[ID("Kr086")].PutIDnextCapture(ID("Rb087"));
  // Rb85
  data[ID("Rb085")].PutNdiv(0,1,0,0);
  // Rb86
  data[ID("Rb086")].PutNdiv(0,1,0,1);
  data[ID("Rb086")].PutIDnextDecay(ID("Sr086"));
  // Rb87
  data[ID("Rb087")].PutNdiv(0,1,0,1);
  data[ID("Rb087")].PutIDnextCapture(ID("Sr088"));
  data[ID("Rb087")].PutIDnextDecay(ID("Sr087"));
  // Sr86
  data[ID("Sr086")].PutNdiv(0,2,0,0);
  real r_sr6_c[]={0.997576,0.002424};
  int id_sr6_c[]={ID("Sr087"),ID("Rb087")};
  data[ID("Sr086")].PutData(1,r_sr6_c,id_sr6_c);
  // Sr87
  data[ID("Sr087")].PutNdiv(0,1,0,0);
  // Sr88
  data[ID("Sr088")].PutNdiv(0,1,0,0);
  // Sr89
  data[ID("Sr089")].PutNdiv(0,1,0,1);
  data[ID("Sr089")].PutIDnextDecay(ID("Y089"));
  // Sr90
  data[ID("Sr090")].PutNdiv(0,1,0,1);
  data[ID("Sr090")].PutIDnextCapture(ID("Y091"));
  data[ID("Sr090")].PutIDnextDecay(ID("Y090"));
  // Y89
  data[ID("Y089")].PutNdiv(0,2,0,0);
  real r_y9_c[]={0.99999812,0.00000188};
  int id_y9_c[]={ID("Y090"),ID("Zr090")};
  data[ID("Y089")].PutData(1,r_y9_c,id_y9_c);
  // Y90
  data[ID("Y090")].PutNdiv(0,1,0,1);
  data[ID("Y090")].PutIDnextDecay(ID("Zr090"));
  // Y91
  data[ID("Y091")].PutNdiv(0,1,0,1);
  data[ID("Y091")].PutIDnextCapture(ID("Zr092"));
  data[ID("Y091")].PutIDnextDecay(ID("Zr091"));
  // Zr90
  data[ID("Zr090")].PutNdiv(0,1,0,0);
  // Zr91
  data[ID("Zr091")].PutNdiv(0,1,0,0);
  // Zr92
  data[ID("Zr092")].PutNdiv(0,1,0,0);
  // Zr93
  data[ID("Zr093")].PutNdiv(0,1,0,2);
  real r_zr3_d[]={0.05,0.95};
  int id_zr3_d[]={ID("Nb093"),ID("Nb093m")};
  data[ID("Zr093")].PutData(3,r_zr3_d,id_zr3_d);
  // Zr94
  data[ID("Zr094")].PutNdiv(0,1,0,0);
  // Zr95
  data[ID("Zr095")].PutNdiv(0,1,0,2);
  real r_zr5_d[]={0.999775,0.000225};
  int id_zr5_d[]={ID("Nb095"),ID("Mo095")};
  data[ID("Zr095")].PutData(3,r_zr5_d,id_zr5_d);
  // Zr96
  // Nb93
  data[ID("Nb093")].PutNdiv(0,1,0,0);
  // Nb93m
  data[ID("Nb093m")].PutNdiv(0,0,0,1);
  data[ID("Nb093m")].PutIDnextDecay(ID("Nb093"));
  // Nb94
  data[ID("Nb094")].PutNdiv(0,2,0,1);
  real r_nb4_c[]={0.999025, 0.000975};
  string nuc_nb4_c[]={"Nb095","Mo095"};
  PutData("Nb094",1,r_nb4_c,nuc_nb4_c);
  data[ID("Nb094")].PutIDnextDecay(ID("Mo094"));
  // Nb95
  data[ID("Nb095")].PutNdiv(0,1,0,1);
  data[ID("Nb095")].PutIDnextCapture(ID("Mo096"));
  data[ID("Nb095")].PutIDnextDecay(ID("Mo095"));
  // Mo92
  data[ID("Mo092")].PutNdiv(0,2,0,0);
  real r_mo2_c[]={0.1, 0.9};
  string nuc_mo2_c[]={"Nb093","Nb093m"};
  PutData("Mo092",1,r_mo2_c,nuc_mo2_c);
  // Mo94
  data[ID("Mo094")].PutNdiv(0,1,0,0);
  // Mo95
  data[ID("Mo095")].PutNdiv(0,1,0,0);
  // Mo96
  data[ID("Mo096")].PutNdiv(0,1,0,0);
  // Mo97
  data[ID("Mo097")].PutNdiv(0,1,0,0);
  // Mo98
  data[ID("Mo098")].PutNdiv(0,1,0,0);
  // Mo99
  data[ID("Mo099")].PutNdiv(0,1,0,1);
  data[ID("Mo099")].PutIDnextDecay(ID("Tc099"));
  // Mo100
  // Tc99
  data[ID("Tc099")].PutNdiv(0,1,0,0);
  data[ID("Tc099")].PutIDnextCapture(ID("Ru100"));
  // Ru100
  data[ID("Ru100")].PutNdiv(0,1,0,0);
  // Ru101
  data[ID("Ru101")].PutNdiv(0,1,0,0);
  // Ru102
  data[ID("Ru102")].PutNdiv(0,1,0,0);
  // Ru103
  data[ID("Ru103")].PutNdiv(0,1,0,1);
  data[ID("Ru103")].PutIDnextDecay(ID("Rh103"));
  // Ru104
  data[ID("Ru104")].PutNdiv(0,1,0,0);
  // Ru105
  data[ID("Ru105")].PutNdiv(0,1,0,1);
  data[ID("Ru105")].PutIDnextDecay(ID("Rh105"));
  // Ru106
  data[ID("Ru106")].PutNdiv(0,0,0,1);
  data[ID("Ru106")].PutIDnextDecay(ID("Rh106"));
  // Rh103
  data[ID("Rh103")].PutNdiv(0,2,0,0);
  real r_rh3_c[]={0.99600056, 3.99944E-3};
  string nuc_rh3_c[]={"Pd104","Ru104"};
  PutData("Rh103",1,r_rh3_c,nuc_rh3_c);
  // Rh105
  data[ID("Rh105")].PutNdiv(0,2,0,1);
  data[ID("Rh105")].PutIDnextDecay(ID("Pd105"));
  real r_rh5_c[]={0.687, 0.313};
  int id_rh5_c[]={ID("Rh106"),ID("Pd106")};
  data[ID("Rh105")].PutData(1,r_rh5_c,id_rh5_c);
  // Rh106
  data[ID("Rh106")].PutNdiv(0,0,0,1);
  data[ID("Rh106")].PutIDnextDecay(ID("Pd106"));
  // Pd104
  data[ID("Pd104")].PutNdiv(0,1,0,0);
  // Pd105
  data[ID("Pd105")].PutNdiv(0,1,0,0);
  // Pd106
  data[ID("Pd106")].PutNdiv(0,1,0,0);
  // Pd107
  data[ID("Pd107")].PutNdiv(0,1,0,1);
  data[ID("Pd107")].PutIDnextDecay(ID("Ag107"));
  // Pd108
  data[ID("Pd108")].PutNdiv(0,1,0,0);
  data[ID("Pd108")].PutIDnextCapture(ID("Ag109"));
  // Pd110
  // Ag107
  // Ag109
  data[ID("Ag109")].PutNdiv(0,3,0,0);
  real r_ag9_c[]={0.002844,0.052,0.945156};
  //real r_ag9_c[]={0.0024288,0.1904,0.807171}; // FAST
  int id_ag9_c[]={ID("Pd110"),ID("Ag110m"),ID("Cd110")};
  data[ID("Ag109")].PutData(1,r_ag9_c,id_ag9_c);
  // Ag110m
  data[ID("Ag110m")].PutNdiv(0,0,0,2);
  real r_ag0m_d[]={0.000042, 0.999958};
  int id_ag0m_d[]={ID("Pd110"),ID("Cd110")};
  data[ID("Ag110m")].PutData(3,r_ag0m_d,id_ag0m_d);
  // Cd110
  data[ID("Cd110")].PutNdiv(0,1,0,0);
  // Cd111
  data[ID("Cd111")].PutNdiv(0,1,0,0);
  // Cd112
  data[ID("Cd112")].PutNdiv(0,2,0,0);
  real r_cd2_c[]={0.6628, 0.3372};
  int id_cd2_c[]={ID("Cd113m"),ID("Cd113")};
  data[ID("Cd112")].PutData(1,r_cd2_c,id_cd2_c);
  // Cd113
  data[ID("Cd113")].PutNdiv(0,1,0,1);
  data[ID("Cd113")].PutIDnextDecay(ID("In113"));
  // Cd113m
  data[ID("Cd113m")].PutNdiv(0,1,0,2);
  data[ID("Cd113m")].PutIDnextCapture(ID("Cd114"));
  real r_cd3m_d[]={0.001, 0.999};
  int id_cd3m_d[]={ID("Cd113"), ID("In113")};
  data[ID("Cd113m")].PutData(3,r_cd3m_d,id_cd3m_d);
  // Cd114
  data[ID("Cd114")].PutNdiv(0,1,0,0);
  data[ID("Cd114")].PutIDnextCapture(ID("In115"));
  // Cd116
  // In113
  // In115
  data[ID("In115")].PutNdiv(0,1,0,0);
  data[ID("In115")].PutIDnextCapture(ID("Sn116"));
  // Sn116
  data[ID("Sn116")].PutNdiv(0,1,0,0);
  // Sn117
  data[ID("Sn117")].PutNdiv(0,1,0,0);
  // Sn118
  data[ID("Sn118")].PutNdiv(0,2,0,0);
  real r_sn8_c[]={0.955, 0.045 };
  int id_sn8_c[]={ID("Sn119"), ID("Sn119m")};
  data[ID("Sn118")].PutData(1,r_sn8_c,id_sn8_c);
  // Sn119
  data[ID("Sn119")].PutNdiv(0,1,0,0);
  // Sn119m
  data[ID("Sn119m")].PutNdiv(0,0,0,1);
  data[ID("Sn119m")].PutIDnextDecay(ID("Sn119"));
  // Sn120
  data[ID("Sn120")].PutNdiv(0,2,0,0);
  real r_sn0_c[]={0.993, 0.007};
  int id_sn0_c[]={ID("Sn121"), ID("Sn121m")};
  data[ID("Sn120")].PutData(1,r_sn0_c,id_sn0_c);
  // Sn121
  data[ID("Sn121")].PutNdiv(0,0,0,1);
  data[ID("Sn121")].PutIDnextDecay(ID("Sb121"));
  // Sn121m
  data[ID("Sn121m")].PutNdiv(0,0,0,2);
  real r_sn1m_d[]={0.776, 0.224};
  int id_sn1m_d[]={ID("Sn121"), ID("Sb121")};
  data[ID("Sn121m")].PutData(3,r_sn1m_d,id_sn1m_d);
  // Sn122
  data[ID("Sn122")].PutNdiv(0,2,0,0);
  real r_sn2_c[]={0.005,0.995};
  int id_sn2_c[]={ID("Sn123"),ID("Sb123")};
  data[ID("Sn122")].PutData(1,r_sn2_c,id_sn2_c);
  // Sn123
  data[ID("Sn123")].PutNdiv(0,1,0,1);
  data[ID("Sn123")].PutIDnextDecay(ID("Sb123"));
  // Sn124
  data[ID("Sn124")].PutNdiv(0,1,0,0);
  data[ID("Sn124")].PutIDnextCapture(ID("Sb125"));
  // Sn126
  data[ID("Sn126")].PutNdiv(0,2,0,1);
  data[ID("Sn126")].PutIDnextDecay(ID("Sb126m"));
  real r_sn6_c[]={0.861,0.139};
  int id_sn6_c[]={ID("I127"),ID("Te127m")};
  data[ID("Sn126")].PutData(1,r_sn6_c,id_sn6_c);
  // Sb121
  data[ID("Sb121")].PutNdiv(0,2,0,0);
  real r_sb1_c[]={0.024,0.976};
  int id_sb1_c[]={ID("Sn122"),ID("Te122")};
  data[ID("Sb121")].PutData(1,r_sb1_c,id_sb1_c);
  // Sb123
  data[ID("Sb123")].PutNdiv(0,2,0,0);
  //real r_sb3_c[]={0.0026,0.9982}; // The sum is not 1
  real r_sb3_c[]={0.0027,0.9973}; // Re-calculated
  int id_sb3_c[]={ID("Te124"),ID("Sb124")};
  data[ID("Sb123")].PutData(1,r_sb3_c,id_sb3_c);
  // Sb124
  data[ID("Sb124")].PutNdiv(0,1,0,1);
  data[ID("Sb124")].PutIDnextDecay(ID("Te124"));
  // Sb125
  data[ID("Sb125")].PutNdiv(0,1,0,2);
  real r_sb5_d[]={0.782,0.218};
  int id_sb5_d[]={ID("Te125"), ID("Te125m")};
  data[ID("Sb125")].PutData(3,r_sb5_d,id_sb5_d);
  // Sb126
  data[ID("Sb126")].PutNdiv(0,2,0,1); 
  data[ID("Sb126")].PutIDnextDecay(ID("Te126"));
  real r_sb6_c[]={0.861,0.139};
  int id_sb6_c[]={ID("I127"),ID("Te127m")};
  data[ID("Sb126")].PutData(1,r_sb6_c,id_sb6_c);
  // Sb126m
  data[ID("Sb126m")].PutNdiv(0,2,0,2);
  real r_sb6m_d[]={0.14, 0.86};
  int id_sb6m_d[]={ID("Sb126"), ID("Te126")};
  real r_sb6m_c[]={0.861,0.139};
  int id_sb6m_c[]={ID("I127"),ID("Te127m")};
  data[ID("Sb126m")].PutData(3,r_sb6m_d,id_sb6m_d);
  data[ID("Sb126m")].PutData(1,r_sb6m_c,id_sb6m_c);
  // Te122
  data[ID("Te122")].PutNdiv(0,2,0,0);
  real r_te2_c[]={0.676,0.324};
  int id_te2_c[]={ID("Te123"), ID("Te123m")};
  data[ID("Te122")].PutData(2,r_te2_c,id_te2_c);
  // Te123
  data[ID("Te123")].PutNdiv(0,1,0,1);
  data[ID("Te123")].PutIDnextDecay(ID("Sb123"));
  // Te123m
  data[ID("Te123m")].PutNdiv(0,0,0,1);
  data[ID("Te123m")].PutIDnextDecay(ID("Te123"));
  // Te124
  data[ID("Te124")].PutNdiv(0,2,0,0);
  real r_te4_c[]={0.994,0.006};
  int id_te4_c[]={ID("Te125"), ID("Te125m")};
  data[ID("Te124")].PutData(2,r_te4_c,id_te4_c);
  // Te125
  data[ID("Te125")].PutNdiv(0,1,0,0);
  // Te125m
  data[ID("Te125m")].PutNdiv(0,0,0,1);
  data[ID("Te125m")].PutIDnextDecay(ID("Te125"));
  // Te126
  data[ID("Te126")].PutNdiv(0,2,0,0);
  real r_te6_c[]={0.13,0.87};
  int id_te6_c[]={ID("Te127m"),ID("I127")};
  data[ID("Te126")].PutData(1,r_te6_c,id_te6_c);
  // Te127m
  data[ID("Te127m")].PutNdiv(0,1,0,1);
  data[ID("Te127m")].PutIDnextCapture(ID("Te128"));
  data[ID("Te127m")].PutIDnextDecay(ID("I127"));
  // Te128
  data[ID("Te128")].PutNdiv(0,2,0,0);
  real r_te8_c[]={0.071,0.929};
  int id_te8_c[]={ID("Te129m"),ID("I129")};
  data[ID("Te128")].PutData(2,r_te8_c,id_te8_c);
  // Te129m
  data[ID("Te129m")].PutNdiv(0,1,0,1);
  data[ID("Te129m")].PutIDnextCapture(ID("Te130"));
  data[ID("Te129m")].PutIDnextDecay(ID("I129"));
  // Te130
  data[ID("Te130")].PutNdiv(0,1,0,0);
  data[ID("Te130")].PutIDnextCapture(ID("I131"));
  // Te132
  data[ID("Te132")].PutNdiv(0,0,0,1);
  data[ID("Te132")].PutIDnextDecay(ID("Xe132"));
  // I127
  data[ID("I127")].PutNdiv(0,2,0,0);
  real r_i7_c[]={0.939,0.061};
  string nuc_i7_c[]={"Xe128","Te128"};
  PutData("I127",1,r_i7_c,nuc_i7_c);
  // I129
  data[ID("I129")].PutNdiv(0,2,0,1);
  real r_i9_c[]={0.88811,0.11189};
  int id_i9_c[]={ID("I130"),ID("Xe130")};
  data[ID("I129")].PutData(1,r_i9_c,id_i9_c);
  data[ID("I129")].PutIDnextDecay(ID("Xe129"));
  // I130
  data[ID("I130")].PutNdiv(0,1,0,1);
  data[ID("I130")].PutIDnextDecay(ID("Xe130"));
  // I131
  data[ID("I131")].PutNdiv(0,1,0,1);
  data[ID("I131")].PutIDnextDecay(ID("Xe131"));
  data[ID("I131")].PutIDnextCapture(ID("Xe132"));
  // I135
  data[ID("I135")].PutNdiv(0,1,0,1);
  data[ID("I135")].PutIDnextDecay(ID("Xe135"));
  data[ID("I135")].PutIDnextCapture(ID("Xe136"));
  // Xe126
  // Xe128
  data[ID("Xe128")].PutNdiv(0,1,0,0);
  // Xe129
  data[ID("Xe129")].PutNdiv(0,1,0,0);
  // Xe130
  data[ID("Xe130")].PutNdiv(0,1,0,0);
  // Xe131
  data[ID("Xe131")].PutNdiv(0,1,0,0);
  // Xe132
  data[ID("Xe132")].PutNdiv(0,1,0,0);
  // Xe133
  data[ID("Xe133")].PutNdiv(0,1,0,1);
  data[ID("Xe133")].PutIDnextDecay(ID("Cs133"));
  // Xe134
  data[ID("Xe134")].PutNdiv(0,2,0,0);
  real r_xe4_c[]={0.99999952,0.00000048};
  int id_xe4_c[]={ID("Xe135"),ID("Cs135")};
  data[ID("Xe134")].PutData(1,r_xe4_c,id_xe4_c);
  // Xe135
  data[ID("Xe135")].PutNdiv(0,1,0,1);
  data[ID("Xe135")].PutIDnextDecay(ID("Cs135"));
  // Xe136
  data[ID("Xe136")].PutNdiv(0,1,0,0);
  data[ID("Xe136")].PutIDnextCapture(ID("Cs137"));
  // Cs133
  data[ID("Cs133")].PutNdiv(0,1,0,0);
  // Cs134
  data[ID("Cs134")].PutNdiv(0,1,0,1);
  data[ID("Cs134")].PutIDnextDecay(ID("Ba134"));
  // Cs135
  data[ID("Cs135")].PutNdiv(0,1,0,1);
  data[ID("Cs135")].PutIDnextDecay(ID("Ba135"));
  // Cs136
  data[ID("Cs136")].PutNdiv(0,1,0,1);
  data[ID("Cs136")].PutIDnextDecay(ID("Ba136"));
  // Cs137
  data[ID("Cs137")].PutNdiv(0,1,0,2);
  data[ID("Cs137")].PutIDnextCapture(ID("Ba138"));
  real r_cs7_d[]={0.053, 0.947};
  int id_cs7_d[]={ID("Ba137"),ID("Ba137m")};
  data[ID("Cs137")].PutData(3,r_cs7_d,id_cs7_d);
  // Ba134
  data[ID("Ba134")].PutNdiv(0,1,0,0);
  // Ba135
  data[ID("Ba135")].PutNdiv(0,1,0,0);
  // Ba136
  data[ID("Ba136")].PutNdiv(0,2,0,0);
  real r_ba6_c[]={0.975, 0.025};
  int id_ba6_c[]={ID("Ba137"), ID("Ba137m")};
  data[ID("Ba136")].PutData(1,r_ba6_c,id_ba6_c);
  // Ba137
  data[ID("Ba137")].PutNdiv(0,1,0,0);
  // Ba137m
  data[ID("Ba137m")].PutNdiv(0,0,0,1);
  data[ID("Ba137m")].PutIDnextDecay(ID("Ba137"));
  // Ba138
  data[ID("Ba138")].PutNdiv(0,1,0,0);
  data[ID("Ba138")].PutIDnextCapture(ID("La139"));
  // Ba140
  data[ID("Ba140")].PutNdiv(0,0,0,1);
  data[ID("Ba140")].PutIDnextDecay(ID("La140"));
  // La139
  data[ID("La139")].PutID(ID("La139"));
  data[ID("La139")].PutNdiv(0,1,0,0);
  // La140
  data[ID("La140")].PutNdiv(0,1,0,1);
  data[ID("La140")].PutIDnextCapture(ID("Ce141"));
  data[ID("La140")].PutIDnextDecay(ID("Ce140"));
  // Ce140
  data[ID("Ce140")].PutNdiv(0,1,0,0);
  // Ce141
  data[ID("Ce141")].PutNdiv(0,1,0,1);
  data[ID("Ce141")].PutIDnextDecay(ID("Pr141"));
  // Ce142
  data[ID("Ce142")].PutNdiv(0,1,0,0);
  // Ce143
  data[ID("Ce143")].PutNdiv(0,1,0,1);
  data[ID("Ce143")].PutIDnextDecay(ID("Pr143"));
  // Ce144
  data[ID("Ce144")].PutNdiv(0,1,0,1);
  data[ID("Ce144")].PutIDnextCapture(ID("Nd145"));
  data[ID("Ce144")].PutIDnextDecay(ID("Pr144"));
  // Pr141
  data[ID("Pr141")].PutNdiv(0,2,0,0);
  real r_pr1_c[]={0.999836,0.000164};
  string nuc_pr1_c[]={"Nd142","Ce142"};
  PutData("Pr141",1,r_pr1_c,nuc_pr1_c);
  // Pr143
  data[ID("Pr143")].PutNdiv(0,1,0,1);
  data[ID("Pr143")].PutIDnextCapture(ID("Nd144"));
  data[ID("Pr143")].PutIDnextDecay(ID("Nd143"));
  // Pr144
  data[ID("Pr144")].PutNdiv(0,0,0,1);
  data[ID("Pr144")].PutIDnextDecay(ID("Nd144"));
  // Nd142
  data[ID("Nd142")].PutNdiv(0,1,0,0);
  // Nd143
  data[ID("Nd143")].PutNdiv(0,1,0,0);
  // Nd144
  data[ID("Nd144")].PutNdiv(0,1,0,0);
  // Nd145
  data[ID("Nd145")].PutNdiv(0,1,0,0);
  // Nd146
  data[ID("Nd146")].PutNdiv(0,1,0,0);
  // Nd147
  data[ID("Nd147")].PutNdiv(0,1,0,1);
  data[ID("Nd147")].PutIDnextDecay(ID("Pm147"));
  // Nd148
  data[ID("Nd148")].PutNdiv(0,1,0,0);
  data[ID("Nd148")].PutIDnextCapture(ID("Pm149"));
  // Nd150
  data[ID("Nd150")].PutNdiv(0,1,0,0);
  data[ID("Nd150")].PutIDnextCapture(ID("Pm151"));
  // Pm147
  data[ID("Pm147")].PutNdiv(0,2,0,1);
  real r_pm7_c[]={0.53,0.47};
  int id_pm7_c[]={ID("Pm148"),ID("Pm148m")};
  data[ID("Pm147")].PutData(1,r_pm7_c,id_pm7_c);
  data[ID("Pm147")].PutIDnextDecay(ID("Sm147"));
  // Pm148
  data[ID("Pm148")].PutNdiv(0,1,0,1);
  data[ID("Pm148")].PutIDnextDecay(ID("Sm148"));
  // Pm148m
  data[ID("Pm148m")].PutNdiv(0,1,0,2);
  data[ID("Pm148m")].PutIDnextCapture(ID("Pm149"));
  real r_pm8m_d[]={0.046, 0.954};
  string nuc_pm8m_d[]={"Pm148","Sm148"};
  PutData("Pm148m",3,r_pm8m_d,nuc_pm8m_d);
  // Pm149
  data[ID("Pm149")].PutNdiv(0,1,0,1);
  data[ID("Pm149")].PutIDnextCapture(ID("Sm150"));
  data[ID("Pm149")].PutIDnextDecay(ID("Sm149"));
  // Pm151
  data[ID("Pm151")].PutNdiv(0,1,0,1);
  data[ID("Pm151")].PutIDnextCapture(ID("Sm152"));
  data[ID("Pm151")].PutIDnextDecay(ID("Sm151"));
  // Sm147
  data[ID("Sm147")].PutNdiv(0,1,0,0);
  // Sm148
  data[ID("Sm148")].PutNdiv(0,1,0,0);
  // Sm149
  data[ID("Sm149")].PutNdiv(0,1,0,0);
  // Sm150
  data[ID("Sm150")].PutNdiv(0,1,0,0);
  // Sm151
  data[ID("Sm151")].PutNdiv(0,1,0,1);
  data[ID("Sm151")].PutIDnextDecay(ID("Eu151"));
  // Sm152
  data[ID("Sm152")].PutNdiv(0,1,0,0);
  //data[ID("Sm152")].PutIDnextCapture(ID("Eu153")); // ??????
  // Sm153
  data[ID("Sm153")].PutNdiv(0,1,0,1);
  data[ID("Sm153")].PutIDnextDecay(ID("Eu153"));
  // Sm154
  data[ID("Sm154")].PutNdiv(0,1,0,0);
  data[ID("Sm154")].PutIDnextCapture(ID("Eu155"));
  // Eu151
  data[ID("Eu151")].PutNdiv(0,3,0,0);
  //real r_eu1_c[]={0.645, 0.2698}; // The sum is not 1
  // (SRAC-ChainJ33, Eu-151=>Sm-152 is neglected)
  real r_eu1_c[]={0.644, 0.27056, 0.08544};
  string id_eu1_c[]={"Eu152","Gd152","Sm152"};
  PutData("Eu151",1,r_eu1_c,id_eu1_c);
  // Eu152 
  // (SRAC-ChainJ33, Eu-152=>Sm-152 is neglected)
  PutNdiv("Eu152",0,1,0,2);
  real r_eu2_d[]={0.27,0.73}; 
  string id_eu2_d[]={"Gd152","Sm152"};
  PutData("Eu152",3,r_eu2_d,id_eu2_d);
  // Eu153
  data[ID("Eu153")].PutNdiv(0,1,0,0);
  // Eu154
  //  In SRAC, path to Sm-154 is neglected.
  data[ID("Eu154")].PutNdiv(0,1,0,2);
  real r_eu4_d[]={0.9998, 0.0002};
  string id_eu4_d[]={"Gd154","Sm154"};
  PutData("Eu154",3,r_eu4_d,id_eu4_d);
  // Eu155
  data[ID("Eu155")].PutNdiv(0,1,0,1);
  data[ID("Eu155")].PutIDnextDecay(ID("Gd155"));
  // Eu156
  data[ID("Eu156")].PutNdiv(0,1,0,1);
  data[ID("Eu156")].PutIDnextDecay(ID("Gd156"));
  // Eu157
  data[ID("Eu157")].PutNdiv(0,1,0,1);
  data[ID("Eu157")].PutIDnextDecay(ID("Gd157"));
  data[ID("Eu157")].PutIDnextCapture(ID("Gd158"));
  // Gd152
  // Gd154
  data[ID("Gd154")].PutNdiv(0,1,0,0);
  // Gd155
  data[ID("Gd155")].PutNdiv(0,1,0,0);
  // Gd156
  data[ID("Gd156")].PutNdiv(0,1,0,0);
  // Gd157
  data[ID("Gd157")].PutNdiv(0,1,0,0);
  // Gd158
  data[ID("Gd158")].PutNdiv(0,1,0,0);
  data[ID("Gd158")].PutIDnextCapture(ID("Tb159"));
  // Gd160
  // Tb159
  data[ID("Tb159")].PutNdiv(0,1,0,0);
  // Tb160
  data[ID("Tb160")].PutNdiv(0,1,0,1);
  data[ID("Tb160")].PutIDnextCapture(ID("Dy161"));
  data[ID("Tb160")].PutIDnextDecay(ID("Dy160"));
  // Dy160
  data[ID("Dy160")].PutNdiv(0,1,0,0);
  // Dy161
  data[ID("Dy161")].PutNdiv(0,1,0,0); 
  // Dy162
  data[ID("Dy162")].PutNdiv(0,1,0,0);
  // Dy163
  data[ID("Dy163")].PutNdiv(0,1,0,0);
  // Dy164
  data[ID("Dy164")].PutNdiv(0,1,0,0);
  data[ID("Dy164")].PutIDnextCapture(ID("Ho165"));
  // Ho163
  data[ID("Ho163")].PutNdiv(0,0,0,1);
  data[ID("Ho163")].PutIDnextDecay(ID("Dy163"));
  // Ho165
  data[ID("Ho165")].PutNdiv(0,2,0,0);
  real r_ho5_c[]={0.9555, 0.0445};
  int id_ho5_c[]={ID("Er166"), ID("Ho166m")};
  data[ID("Ho165")].PutData(1,r_ho5_c,id_ho5_c);
  // Er162
  data[ID("Er162")].PutNdiv(0,1,0,0);
  data[ID("Er162")].PutIDnextCapture(ID("Ho163"));
  // Er164
  data[ID("Er164")].PutNdiv(0,1,0,0);
  data[ID("Er164")].PutIDnextCapture(ID("Ho165"));
  // Er166
  data[ID("Er166")].PutNdiv(0,1,0,0);
  // Er167
  data[ID("Er167")].PutNdiv(0,1,0,0);
  // Er168
  // Er170
  // +++ Burnable poison
  // Hf176
  data[ID("Hf176")].PutNdiv(0,1,0,0);
  // Hf177
  data[ID("Hf177")].PutNdiv(0,1,0,0);
  // Hf178
  data[ID("Hf178")].PutNdiv(0,1,0,0);
  // Hf179
  data[ID("Hf179")].PutNdiv(0,1,0,0);
  // Hf180
};

void BurnupChain::OverWritingChainData(string cbglibdir,string filename)
{
  string mdir=cbglibdir+"CBGLIB_BURN/CBG_Chain/"+filename;
  OverWritingChainData(mdir);
};

void BurnupChain::ReadDoseCoefficientData(string cbglibdir)
{
  string filename[]={
    cbglibdir+"CBGLIB_BURN/DoseCoef/"+"lib_dosecoef_ingestion",
    cbglibdir+"CBGLIB_BURN/DoseCoef/"+"lib_dosecoef_inhalation",
  };

  for(int i=0;i<2;i++){
    ifstream fin;
    fin.open(filename[i].data(),ios::in);
    if(fin.fail()){
      cout<<"# Failed to open the file.\n";
      cout<<"# File name is "<<filename[i]<<".\n";
      exit(0);
    };
    int num;
    fin>>num;
    for(int j=0;j<num;j++){
      int mat;
      fin>>mat;
      real data;
      fin>>data;
      if(i==0){
	dosecoef_ingestion[mat]=data;
      }else{
	dosecoef_inhalation[mat]=data;
      };
    };
    fin.close();
  };

};

void BurnupChain::OverWritingChainData(string filename)
{
  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"File name is "<<filename<<"\n";
    exit(0);
  };

  vector<int> div(4);
  vector<int> adiv(4);
  int fl=0;
  while(fl!=-1){
    int matno;
    fin>>matno;
    if(matno==-1)break;
    if(matno<10000)matno=midt.GetMATIDFromENDFID(matno);
    data[matno].PutID(matno);
    real halflife;
    fin>>halflife;
    //if(halflife>1e-10)decay_const[matno]=0.693/halflife;
    if(fabs(halflife)>1e-10)decay_const[matno]=0.693/halflife;
    real elp,eem,ehp;
    fin>>elp;
    fin>>eem;
    fin>>ehp;
    data[matno].PutDecayEnergy(0,elp);
    data[matno].PutDecayEnergy(1,eem);
    data[matno].PutDecayEnergy(2,ehp);
    for(int i=0;i<4;i++){
      fin>>div[i];
      adiv[i]=abs(div[i]);
    };
    data[matno].PutNdiv(adiv[0],adiv[1],adiv[2],adiv[3]);
    for(int i=0;i<4;i++){
      if(adiv[i]>0){
        real *tmp1=new real[adiv[i]];
        int *tmp2=new int[adiv[i]];
        tmp1[0]=1.;
        for(int j=0;j<adiv[i];j++){
	  fin>>tmp2[j];
	  if(tmp2[j]<10000)tmp2[j]=midt.GetMATIDFromENDFID(tmp2[j]);
  	  if(adiv[i]!=1||div[i]==-1){
   	    fin>>tmp1[j];
  	  };
        };
        data[matno].PutData(i,tmp1,tmp2);
        delete [] tmp1;
        delete [] tmp2;
      };
    };
  };
  fin.close();
};

void BurnupChain::OverWritingChainData(BCGManager &bm)
{
  // (n,2n) reaction is considered as default

  int r_num=bm.GetNuclide(0).GetRnum(); // 

  vector<int> div(4);
  vector<int> adiv(4);
  
  int sz=bm.GetSize();
  for(int i=0;i<sz;i++){

    if(bm.GetNuclide(i).Flag()){

      int atm=bm.GetNuclide(i).GetAtomicNumber();
      int mas=bm.GetNuclide(i).GetMassNumber();
      int lev=bm.GetNuclide(i).GetExLevel();
      int matno=midt.ID(atm,mas,lev);

      data[matno].PutID(matno);

      int div_decay=bm.GetNuclide(i).GetChannel();

      if(div_decay>0){

        real halflife=bm.GetNuclide(i).GetHalflife();
        //if(halflife>1e-10)decay_const[matno]=0.693/halflife;
        if(fabs(halflife)>1e-10)decay_const[matno]=0.693/halflife;
        //if(halflife!=0.)decay_const[matno]=0.693/halflife;

        real elp=bm.GetNuclide(i).GetDecayEnergy(0);
        real eem=bm.GetNuclide(i).GetDecayEnergy(1);
        real ehp=bm.GetNuclide(i).GetDecayEnergy(2);

        data[matno].PutDecayEnergy(0,elp);
        data[matno].PutDecayEnergy(1,eem);
        data[matno].PutDecayEnergy(2,ehp);

      };

      int div_c=bm.GetNuclide(i).GetReactionChannel(0);
      int div_n2n=0;
      if(r_num>1)div_n2n=bm.GetNuclide(i).GetReactionChannel(1);

      data[matno].PutNdiv(0,div_c,div_n2n,div_decay);

      // (n,g)
      {
      real *tmp1=new real[div_c];
      int *tmp2=new int[div_c];
      for(int k=0;k<div_c;k++){
        int atm2=bm.GetNuclide(i).GetReactionAtomicNumberNext(0,k);
        int mas2=bm.GetNuclide(i).GetReactionMassNumberNext(0,k);
        int lev2=bm.GetNuclide(i).GetReactionExLevelNext(0,k);
        real br=bm.GetNuclide(i).GetReactionBr(0,k);
        int id2=midt.ID(atm2,mas2,lev2);
	tmp2[k]=id2;
	tmp1[k]=br;
      };
      data[matno].PutData(1,tmp1,tmp2);
      delete [] tmp1;
      delete [] tmp2;
      };

      // (n,2n)
      if(r_num>1){
      real *tmp1=new real[div_n2n];
      int *tmp2=new int[div_n2n];
      for(int k=0;k<div_n2n;k++){
        int atm2=bm.GetNuclide(i).GetReactionAtomicNumberNext(1,k);
        int mas2=bm.GetNuclide(i).GetReactionMassNumberNext(1,k);
        int lev2=bm.GetNuclide(i).GetReactionExLevelNext(1,k);
        real br=bm.GetNuclide(i).GetReactionBr(1,k);
        int id2=midt.ID(atm2,mas2,lev2);
	tmp2[k]=id2;
	tmp1[k]=br;
      };
      data[matno].PutData(2,tmp1,tmp2);
      delete [] tmp1;
      delete [] tmp2;
      };

      // (decay)
      {
      real *tmp1=new real[div_decay];
      int *tmp2=new int[div_decay];
      for(int k=0;k<div_decay;k++){
        int atm2=bm.GetNuclide(i).GetAtomicNumberNext(k);
        int mas2=bm.GetNuclide(i).GetMassNumberNext(k);
        int lev2=bm.GetNuclide(i).GetExLevelNext(k);
        real br=bm.GetNuclide(i).GetBr(k);
        int id2=midt.ID(atm2,mas2,lev2);
	tmp2[k]=id2;
	tmp1[k]=br;
      };
      data[matno].PutData(3,tmp1,tmp2);
      delete [] tmp1;
      delete [] tmp2;
      };
    };

  };
};

void BurnupChain::ReadFPYieldDataFromFile(string cbglibdir,string fname)
{
  string mdir=cbglibdir+"CBGLIB_BURN/CBG_FY/"+fname;
  ReadFPYieldDataFromFile(mdir);
};

void BurnupChain::ReadFPYieldDataFromFile(string fname)
{
  ifstream fin;

  fin.open(fname.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<fname<<"\n";
    exit(0);
  };

  bool ed=false;

  while(!ed){

  int fisnuc;
  fin>>fisnuc;

  if(fisnuc==-1){
    ed=true;
  }else{

  vector<int> fisid(fisnuc);
  vector<bool> fis_exist(fisnuc,false);
  for(int i=0;i<fisnuc;i++){
    fin>>fisid[i];
    if(fisid[i]<10000)fisid[i]=midt.GetMATIDFromENDFID(fisid[i]);
    if(data.find(fisid[i])!=data.end())fis_exist[i]=true;
  };

  int fpnuc;
  fin>>fpnuc;
  for(int i=0;i<fisnuc;i++){
    if(fis_exist[i])data[fisid[i]].PutNdivFission(fpnuc);
  };

  for(int i=0;i<fpnuc;i++){
    int fpid;
    real fpy;
    fin>>fpid;
    fin>>fpy;
    if(fpid<10000)fpid=midt.GetMATIDFromENDFID(fpid);
    for(int ii=0;ii<fisnuc;ii++){
      if(fis_exist[ii])data[fisid[ii]].PutIDnextFission(i,fpid,fpy);
    };
  };

  };

  };

  fin.close();

};

void BurnupChain::RetreiveFPYieldData(BCGManager &bm)
{
  // Counting the numgber of treated FPs
  int fpnuc=0;
  vector<int> fpid;
  int sz=bm.GetSize();
  for(int i=0;i<sz;i++){
    if(bm.GetNuclide(i).Flag()){
      int atm=bm.GetNuclide(i).GetAtomicNumber();
      int mas=bm.GetNuclide(i).GetMassNumber();
      int lev=bm.GetNuclide(i).GetExLevel();
      fpid.push_back(midt.ID(atm,mas,lev));
      fpnuc++;
    };
  };

  int ynum=28;
  string yname[]={
    "Th232",
    "Pa231","Pa233",
    "U232","U233","U234","U235","U236","U237","U238",
    "Np236","Np237","Np239",
    "Pu236","Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Am242","Am242m","Am243",
    "Cm242","Cm243","Cm244","Cm245","Cm246"
  };

  // (FPY data of the following nuclide data are used)
  string ytag[]={
    "Th232",
    "Pa231","Th232",
    "U232","U233","U234","U235","U236","U237","U238",
    "U235","Np237","Np237",
    "U236","Pu238","Pu239","Pu240","Pu241","Pu242",
    "Am241","Pu241","Am242m","Am243",
    "Cm242","Cm243","Cm244","Cm245","Cm246"
  };

  for(int i=0;i<ynum;i++){
    int fisid=midt.ID(yname[i]);
    if(data.find(fisid)!=data.end()){
      data[fisid].PutNdivFission(fpnuc);
      for(int j=0;j<fpnuc;j++){
        int atm,mas,lev;
        midt.GetParameter(fpid[j],atm,mas,lev);
	real fpy=bm.GetYield(ytag[i],atm,mas,lev);
        data[fisid].PutIDnextFission(j,fpid[j],fpy);
      };
    };
  };

};

void BurnupChain::ChangeFPYieldData(string nuc1, string nuc2)
{
  int id1=ID(nuc1);
  int id2=ID(nuc2);

  int ndiv=data[id1].GetNdiv(0);
  if(ndiv>0){
    data[id2].PutNdiv(0,ndiv);
    for(int i=0;i<ndiv;i++){
      data[id2].PutIDnextFission(i,data[id1].GetIDnext(0,i),data[id1].GetRatio(0,i));
    };
  };
};

void BurnupChain::ShowFPYield(int id)
{
  cout<<"# Fission yield data for "<<midt.Name(id)<<"\n";
  data[id].ShowFissionYieldData(midt);
};

// +++ Decay constant manupulation

void BurnupChain::WriteDecayConstant()
{
  cout.setf(ios::showpoint);
  cout.precision(6);
  map<int,real>::iterator it=decay_const.begin();
  while(it!=decay_const.end()){
    if(it->second!=0.){
      cout<<"  "<<it->first<<"\n";
      cout<<"  "<<it->second<<"\n";
    };
    it++;
  };
  cout<<"-1\n";
};

void BurnupChain::FactorizeDecayConstant(real factor)
{
  map<int,real>::iterator it=decay_const.begin();
  while(it!=decay_const.end()){
    if(it->second!=0.){
      int i=it->first;
      real tmp=it->second*factor;
      decay_const[i]=tmp;
    };
    it++;
  };
};

void BurnupChain::ReadDecayConstantFromFile(string cbglibdir,string fname)
{
  string mdir=cbglibdir;
  mdir=mdir+"CBGLIB_BURN/decay_constant/";

  mdir.append(fname);
  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"Failed to open the file.\n";
    cout<<"File name is "<<mdir<<"\n";
    exit(0);
  };

  read_decay_constant=true;

  int nuc;
  real par;

  bool ed=false;
  while(!ed){
    fin>>nuc;
    if(nuc==-1){
      ed=true;
    }else{
      if(nuc<10000)nuc=midt.GetMATIDFromENDFID(nuc);
      fin>>par;
      decay_const[nuc]=par;
    };    
  };

  fin.close();
};

void BurnupChain::ShowHalfLife()
{
  real ln2=0.69314718;

  cout<<"#\n# List of decay constant [1/s] & half-life\n#\n";
  map<int,real>::iterator it=decay_const.begin();
  while(it!=decay_const.end()){
    if(it->second!=0.){
      int i=it->first;
      string name=Name(i);
      real dc=it->second;
      real hl=ln2/dc;
      cout<<"# "<<i<<" ";
      WriteOut(name,7);
      cout<<" : "<<dc<<" / ";
      if(hl<60.){
	cout<<hl<<" [Sec.]\n";
      }else if(hl<60*60){
	cout<<hl/60.<<" [Min.]\n";
      }else if(hl<60*60*24){
	cout<<hl/(60*60)<<" [Hours]\n";
      }else if(hl<60*60*24*365){
	cout<<hl/(60*60*24)<<" [Days]\n";
      }else{
	cout<<hl/(60*60*24*365)<<" [Years]\n";
      };
    };
    it++;
  };
};

void BurnupChain::ShowDoseCoefficient()
{
  cout<<"#\n# List of dose coefficient [Sv/Bq]\n#\n";

  {
  cout<<"# Ingestion dose [KEIKOH]\n#\n";
  map<int,real>::iterator it=dosecoef_ingestion.begin();
  while(it!=dosecoef_ingestion.end()){
    if(it->second!=0.){
      int i=it->first;
      string name=Name(i);
      real dc=it->second;
      cout<<"# "<<i<<" ";
      WriteOut(name,6);
      cout.setf(ios::scientific);
      cout.precision(5);
      cout<<" : "<<dc<<"\n";
    };
    it++;
  };
  };

  {
  cout<<"#\n# Inhalation dose [KYU-NYU]\n#\n";
  map<int,real>::iterator it=dosecoef_inhalation.begin();
  while(it!=dosecoef_inhalation.end()){
    if(it->second!=0.){
      int i=it->first;
      string name=Name(i);
      real dc=it->second;
      cout<<"# "<<i<<" ";
      WriteOut(name,6);
      cout.setf(ios::scientific);
      cout.precision(5);
      cout<<" : "<<dc<<"\n";
    };
    it++;
  };
  };

};

void BurnupChain::ShowChainData()
{
  map<int,NuclideChainData>::iterator it=data.begin();
  while(it!=data.end()){
    ShowChainData((*it).first);
    cout<<"#\n";
    it++;
  };
};

void BurnupChain::ShowChainData(string name)
{
  int id=ID(name);
  cout<<"######################################\n";
  cout<<"# Chain data for "<<name<<" ("<<id<<")\n";
  cout<<"######################################\n#\n";
  int ncap=GetNdivCapture(id);
  int nn2n=GetNdivN2N(id);
  int ndec=GetNdivDecay(id);
  if(ncap!=0){
    cout<<"#   [Capture]\n";
    for(int i=0;i<ncap;i++){
      int idn=GetNextIDCapture(id,i);
      string name2=Name(idn);
      real frac=GetRatioCapture(id,i);
      cout<<"#       "<<name2<<" ("<<idn<<") : "<<frac<<"\n";
    };
  };
  if(nn2n!=0){
    cout<<"#   [n,2n]\n";
    for(int i=0;i<nn2n;i++){
      int idn=GetNextIDN2N(id,i);
      string name2=Name(idn);
      real frac=GetRatioN2N(id,i);
      cout<<"#       "<<name2<<" ("<<idn<<") : "<<frac<<"\n";
    };
  };
  if(ndec!=0){
    cout<<"#   [Decay]\n";
    for(int i=0;i<ndec;i++){
      int idn=GetNextIDDecay(id,i);
      string name2=Name(idn);
      real frac=GetRatioDecay(id,i);
      cout<<"#       "<<name2<<" ("<<idn<<") : "<<frac<<"\n";
    };
  };
};

void BurnupChain::ShowChainDataToOrigin(string name)
{
  int id=ID(name);
  cout<<"# Origin of "<<name<<" ("<<id<<")\n";

  map<int,NuclideChainData>::iterator it=data.begin();
  while(it!=data.end()){

    int i=(*it).first;

    int ncap=GetNdivCapture(i);
    int nn2n=GetNdivN2N(i);
    int ndec=GetNdivDecay(i);
    if(ncap!=0){
      for(int j=0;j<ncap;j++){
        int id2=GetNextIDCapture(i,j);
        if(id==id2){
	  cout<<"#   "<<Name(i)<<" ("<<i<<") : ";
	  cout<<"Capture   ";
	  cout<<GetRatioCapture(i,j)<<"\n";
	};
      };
    };
    if(nn2n!=0){
      for(int j=0;j<nn2n;j++){
        int id2=GetNextIDN2N(i,j);
        if(id==id2){
          cout<<"#   "<<Name(i)<<" ("<<i<<") : ";
          cout<<"(n,2n)    ";
          cout<<GetRatioN2N(i,j)<<"\n";
	};
      };
    };
    if(ndec!=0){
      for(int j=0;j<ndec;j++){
	int id2=GetNextIDDecay(i,j);
	if(id==id2){
	  cout<<"#   "<<Name(i)<<" ("<<i<<") : ";
	  cout<<"Decay     ";
	  cout<<GetRatioDecay(i,j)<<"\n";
	};
      };
    };

    it++;
  };

};

void BurnupChain::ShowFPYieldForFP(string name)
{
  cout<<"# Fission yield : "<<name<<"\n";
  cout.setf(ios::scientific);
  cout.precision(5);
  int idt=ID(name);

  map<int,NuclideChainData>::iterator it=data.begin();
  while(it!=data.end()){
    int nd=(*it).second.GetNdiv(0);
    if(nd!=0){
      for(int j=0;j<nd;j++){
	int id=(*it).second.GetIDnext(0,j);
	if(id==idt){
	  cout<<"# "<<Name((*it).first)<<" : "<<(*it).second.GetRatio(0,j)<<"\n";
	};
      };
    };
    it++;
  };

};

bool BurnupChain::IsConnected(int mat1,int mat2)
{
  for(int i=0;i<4;i++){
    int nd=data[mat1].GetNdiv(i);
    for(int j=0;j<nd;j++){
      if(data[mat1].GetIDnext(i,j)==mat2)return true;
    };
  };
  return false;
};

void BurnupChain::CutChain(string nucnam,int ii)
{
  int id=midt.ID(nucnam);
  if(data.find(id)==data.end()){
    cout<<"# Error in BurnupChain::CutChain.\n";
    cout<<"# No nuclide name : "<<nucnam<<"\n";
    exit(0);
  };
  data[id].PutNdiv(ii,0);
};

void BurnupChain::SetZeroDecayConstant(string nucnam)
{
  int id=midt.ID(nucnam);
  if(data.find(id)==data.end()){
    cout<<"# Error in BurnupChain::SetZeroDecayConstant.\n";
    cout<<"# No nuclide name : "<<nucnam<<"\n";
    exit(0);
  };
  //decay_const[id]=0.;
  decay_const[id]=1e-20;
};

void BurnupChain::FactorizeHalfLife(string nucnam,real factor)
{
  int id=midt.ID(nucnam);
  if(data.find(id)==data.end()){
    cout<<"# Error in BurnupChain::DoubleHalfLife.\n";
    cout<<"# No nuclide name : "<<nucnam<<"\n";
    exit(0);
  };
  decay_const[id]*=1/factor;
};

void BurnupChain::SetZeroFissionYield(string nuc_fis, string nuc_fp)
{
  int id_fis=midt.ID(nuc_fis);
  int id_fp=midt.ID(nuc_fp);

  if(data.find(id_fis)==data.end()||data.find(id_fp)==data.end()){
    cout<<"# Error in BurnupChain::SetZeroFissionYield.\n";
    cout<<"# No nuclide name : "<<nuc_fis<<" or "<<nuc_fp<<"\n";
    exit(0);
  };

  int fp_num=data[id_fis].GetNdiv(0);
  for(int i=0;i<fp_num;i++){
    if(data[id_fis].GetIDnext(0,i)==id_fp)data[id_fis].PutIDnextFission(i,id_fp,0.);
  };
};

void BurnupChain::ShowNGBranchingRatio()
{
  cout<<"#\n# (n,g) branching ratio data\n#\n";
  cout.setf(ios::scientific);
  cout.precision(5);

  map<int,NuclideChainData>::iterator it=data.begin();
  while(it!=data.end()){
    int ndiv=(*it).second.GetNdiv(1);
    if(ndiv>1){
      cout<<"# "<<midt.Name((*it).first)<<" ";
      for(int j=0;j<ndiv;j++){
        cout<<" "<<(*it).second.GetRatio(1,j);
      };
      cout<<"\n";
    };
    it++;
  };
};

void BurnupChain::SetPseudoFP(int mat,int fpid)
{
  data[mat].SetPseudoFP(fpid);
  data[fpid].PutID(fpid);
  data[fpid].PutNdiv(0,0,0,0);
  //data[fpid].PutNdiv(0,1,0,0);
  //data[fpid].PutData(1,1.,fpid);
};

real BurnupChain::GetFissionYield(int fid,int fpid)
{
  if(data.find(fid)==data.end())return -1;
  return data[fid].GetFissionYield(fpid);  
};

void BurnupChain::ModifyDecayConstantForSpontaneousFission()
{
  decay_const[922340]*=1.70e-11;
  decay_const[922350]*=7.20e-11;
  decay_const[922360]*=9.00e-10;
  decay_const[922380]*=5.46e-7;
  decay_const[942380]*=1.86e-9;
  decay_const[942390]*=3.10e-12;
  decay_const[942400]*=5.70e-8;
  decay_const[942420]*=5.50e-6;
  decay_const[952410]*=4.30e-12;
  decay_const[952421]*=1.80e-10;
  decay_const[952430]*=3.70e-11;
  decay_const[962420]*=6.10e-8;
  decay_const[962440]*=1.38e-6;
  decay_const[962460]*=2.61e-4;
};

void BurnupChain::SetZeroDecayEnergy()
{
  map<int,NuclideChainData>::iterator it=data.begin();
  while(it!=data.end()){
    (*it).second.PutDecayEnergy(0,0.);
    (*it).second.PutDecayEnergy(1,0.);
    (*it).second.PutDecayEnergy(2,0.);
    it++;
  };
};

void BurnupChain::PutDecayEnergy(string name,real val)
{
  GetNuclideChainData(name).PutDecayEnergy(0,val);
};

void BurnupChain::PutDecayEnergy(int num,string *name,real *val)
{
  for(int i=0;i<num;i++){
    PutDecayEnergy(name[i],val[i]);
  };
};

real NuclideChainData::GetFissionYield(int fpid)
{
  for(int i=0;i<ndiv[0];i++){
    if(idnext[0][i]==fpid)return ratio[0][i];
  };
  return -1.;
};

void NuclideChainData::PutFissionYield(int fpid, real yld)
{
  for(int i=0;i<ndiv[0];i++){
    if(idnext[0][i]==fpid){
      ratio[0][i]=yld;
      return;
    };
  };
  cout<<"# Error in NuclideChainData::PutFissionYield.\n";
  cout<<"# FP nuclide (index : "<<fpid<<") cannot be found.\n";
  exit(0);
};

void NuclideChainData::AddFissionYield(int fpid, real yld)
{
  ndiv[0]+=1;
  idnext[0].push_back(fpid);
  ratio[0].push_back(yld);
};

void BurnupChain::DataPerturbation(string mdir, string filename)
{
  mdir=mdir+filename;

  ifstream fin;
  fin.open(mdir.data(),ios::in);
  if(fin.fail()){
    cout<<"# Error in NuclideChainData::DataPerturbation.\n";
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
    if(mt>18000000){
      int fisid=mat-18000000;
      real yldorg=GetNuclideChainData(fisid).GetFissionYield(mat);
      if(yldorg>=0.){
        real yldnew=yldorg+yldorg*chg;
        GetNuclideChainData(fisid).PutFissionYield(mat,yldnew);
      };
    }else if(mt==8888){
      real dcorg=GetDecayConstant(mat);
      real dcnew=dcorg-dcorg*chg;
      PutDecayConstant(mat,dcnew);      
    }else if(mt>=88880&&mt<=88889){
      int ch=mt-88880;
      int chnum=GetNuclideChainData(mat).GetNdiv(3);
      if(ch>=chnum){
	cout<<"# Warning : inconsistent decay chennel : "<<mat<<"\n";
	cout<<"#    Channel : "<<chnum<<" / "<<ch<<"\n";
      }else{
        real brorg=GetNuclideChainData(mat).GetRatio(3,ch);
        real brnew=brorg+brorg*chg;
        GetNuclideChainData(mat).PutRatio(3,ch,brnew);
      };
    }else if(mt>=99990&&mt<=99999){
      int type=mt-99990;
      real deorg=GetNuclideChainData(mat).GetDecayEnergy(type);
      real denew=deorg+deorg*chg;
      GetNuclideChainData(mat).PutDecayEnergy(type,denew);
    }else{
      cout<<"# Error in NuclideChainData::DataPerturbation.\n";
      cout<<"# Not coded for MT="<<mt<<"\n";
      exit(0);
    };
  };
  fin.close();

};

// (For advanced pseudo FP model)

void BurnupChain::AddPseudoFPSeriesData(int fpnum, int *matid_fp, real *dc)
{
  for(int i=0;i<fpnum;i++){
    data[matid_fp[i]].PutID(matid_fp[i]);
    int decay_channel=1;
    if(i==fpnum-1)decay_channel=0;
    data[matid_fp[i]].PutNdiv(0,1,1,decay_channel); // f,c,2n,d
    data[matid_fp[i]].PutData(1,1.,matid_fp[i]); // The same nuclide is generated by capture reaction. 
    data[matid_fp[i]].PutData(2,1.,matid_fp[i]); // The same nuclide is generated by (n,2n) reaction.
    if(i!=fpnum-1)data[matid_fp[i]].PutData(3,1.,matid_fp[i+1]); // Another pseudo FP is generated by decay.
    PutDecayConstant(matid_fp[i], dc[i]);
  };
};

void BurnupChain::AddPseudoFPSeriesData1(int fpnum, int fissilenum, int *matid_fp, real *dc)
{
  for(int i=0;i<fissilenum;i++){
    for(int j=0;j<fpnum;j++){
      data[matid_fp[i*fpnum+j]].PutID(matid_fp[i*fpnum+j]);
      data[matid_fp[i*fpnum+j]].PutNdiv(0,1,1,1); // f,c,2n,d
      data[matid_fp[i*fpnum+j]].PutData(1,1.,matid_fp[i*fpnum+j]); // The same nuclide is generated by capture reaction.
      data[matid_fp[i*fpnum+j]].PutData(2,1.,matid_fp[i*fpnum+j]); // The same nuclide is generated by (n,2n) reaction.
      if(j!=fpnum-1){
	data[matid_fp[i*fpnum+j]].PutData(3,1.,matid_fp[i*fpnum+j+1]); // Another pseudo FP is generated by decay.
      }else{
	data[matid_fp[i*fpnum+j]].PutData(3,1.,909990); // Another pseudo FP is generated by decay.
      };
      PutDecayConstant(matid_fp[i*fpnum+j], dc[i*fpnum+j]);
    };
  };
};

void BurnupChain::AddYieldDataToPseudoFP(int fissile_id, int fpnum, int *matid_fp, real *yield_fp)
{
  if(data[fissile_id].GetID()==0){
    cout<<"# Warning in BurnupChain::AddYieldDataToPseudoFP.\n";
    cout<<"# No fission nuclide whose ID is "<<fissile_id<<" does NOT exist in the burnup chain.\n";
  }else{
    for(int i=0;i<fpnum;i++){
      data[fissile_id].AddFissionYield(matid_fp[i], yield_fp[i]);
      //cout<<i<<" "<<fissile_id<<" "<<matid_fp[i]<<" "<<yield_fp[i]<<"\n";
    };
  };
};
void BurnupChain::AddYieldDataToPseudoFP1(int *fissile_id, int fissilenum, int fpnum, int *matid_fp, real *yield_fp)
{
  for (int i=0;i<fissilenum;i++){
    if(data[fissile_id[i]].GetID()==0){
	cout<<"# Warning in BurnupChain::AddYieldDataToPseudoFP.\n";
	cout<<"# No fission nuclide whose ID is "<<fissile_id<<" does NOT exist in the burnup chain.\n";
      }else{
	for(int j=0;j<fpnum;j++){
	  data[fissile_id[i]].AddFissionYield(matid_fp[i*fpnum+j], yield_fp[j]);
	};
      };
  };
};
#if 0
void BurnupChain::AddPseudoFPD5(int fissile_num, int *matid_fp, real *yield_fp, real *dc)
{
   for(int i=0;i<fissile_num-1;i++){
    int id=matid_fp[i]-99000000;
    if(data[id].GetID()!=0){
      data[id].AddFissionYield(matid_fp[i], yield_fp[i]);
      data[matid_fp[i]].PutID(matid_fp[i]);
      data[matid_fp[i]].PutNdiv(0,1,1,1); // f,c,2n,d
      data[matid_fp[i]].PutData(1,1.,matid_fp[i]); // The same nuclide is generated by capture reaction. 
      data[matid_fp[i]].PutData(2,1.,matid_fp[i]); // The same nuclide is generated by (n,2n) reaction.
      data[matid_fp[i]].PutData(3,1.,matid_fp[i+1]); // Another pseudo FP is generated by decay.
      PutDecayConstant(matid_fp[i], dc[i]);
    };
   };
   int id=matid_fp[4]-99000000;
    if(data[id].GetID()!=0){
      data[id].AddFissionYield(matid_fp[4], yield_fp[4]);
      data[matid_fp[4]].PutID(matid_fp[4]);
      data[matid_fp[4]].PutNdiv(0,1,1,1); // f,c,2n,d
      data[matid_fp[4]].PutData(1,1.,matid_fp[4]); // The same nuclide is generated by capture reaction. 
      data[matid_fp[4]].PutData(2,1.,matid_fp[4]); // The same nuclide is generated by (n,2n) reaction.
      data[matid_fp[4]].PutData(3,1.,matid_fp[4]); // Another pseudo FP is generated by decay.
    };
};
#endif
