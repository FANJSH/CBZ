#include <cstdlib>
#include "FMG_Tool.h"

FMG_Tool::FMG_Tool()
{
  file_grp="";
  dir_acef="";
  slowing_down_calc=false;
  pnt=-1;

  wgt_inverse_e=false;
  
  int numbg=9;
  real bgxsinp[]={1e+10, 1e+6, 1e+5, 1e+4, 1e+3, 1e+2, 1e+1, 1.0, 0.1};
  SetBackgroundCrossSection(numbg, bgxsinp);
};

void FMG_Tool::SetBackgroundCrossSection(int numbg, real *xsinp)
{
  num_bgxs=numbg;
  bgxs.resize(num_bgxs);
  for(int i=0;i<num_bgxs;i++){
    bgxs[i]=xsinp[i];
  };
};

void FMG_Tool::SetGroupStructureFile(string inp)
{
  file_grp=inp;
};

void FMG_Tool::SetAceFileDirectory(string inp)
{
  dir_acef=inp;
};

void FMG_Tool::SlowingDownCalculationOn(real wgt, real fehi_inp)
{
  slowing_down_calc=true;
  sd_wgt=wgt;
  fehi=fehi_inp;
};

void FMG_Tool::ReadEnergyGroupStructure()
{
  ebnd.clear();
  
  ifstream fin2;
  fin2.open(file_grp.data(),ios::in);
  if(fin2.fail()){
    cout<<"# Failed to open the group structure file.\n";
  };
  fin2>>grp;
  for(int g=0;g<grp+1;g++){
    float tmp;
    fin2>>tmp;
    ebnd.push_back(tmp);
  };
};

void FMG_Tool::ReadWeightFunction(string inp)
{
  ifstream fin2;
  fin2.open(inp.data(),ios::in);
  if(fin2.fail()){
    cout<<"# Failed to open the weight function data file.\n";
  };

  fin2>>pnt;
  epnt.resize(pnt);
  wgt.resize(pnt);
  for(int g=0;g<pnt;g++){
    fin2>>epnt[g];
    fin2>>wgt[g];
  };
};

void FMG_Tool::PutMaterialInfo(int i, string *sinp, real *dinp, string *nameinp)
{
  nucnum=i;
  aceid.resize(nucnum);
  den.resize(nucnum);
  nucname.resize(nucnum);
  for(int i=0;i<nucnum;i++){
    aceid[i]=sinp[i];
    den[i]=dinp[i];
    nucname[i]=nameinp[i];
  };
};

void FMG_Tool::input_generator_lanlmix(string filenameadd)
{
  ofstream fout;
  fout.open("./input",ios::out);
  if(fout.fail()){
    cout<<"Failed to open the OUTPUT file.\n";
    exit(1);
  };

  ReadEnergyGroupStructure();
  
  string tmp1;

  // (ultra_fine_group)
  fout<<"ultra_fine_group\n";
  fout<<" 2.0e+07            10000  EL \n";
  fout<<" 52475.0            56000  EL \n";
  fout<<" 9118.8             12000  EL \n";
  fout<<" 4307.4             12000  EL \n";
  fout<<" 961.12             40000  EL \n";
  fout<<" 130.07             60000  EL \n";
  // ... lanl & Bigten
  //fout<<" 0.32242            50000  EL \n";
  //fout<<" 1.0e-5 \n";
  // ... FNS/Cu
  fout<<" 0.32242\n";
  fout<<"\n";

  // (multi_group)
  fout<<"multi_group\n";
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<grp+1;i++){
    fout<<" "<<ebnd[i];
    if(i%6==5)fout<<"\n";
  };
  fout<<"\n";
  fout<<"\n";

  // (default_spectrum)
  fout<<"default_spectrum\n";
  // ... FNS/Cu
  fout<<"     1/E\n";
  // ... lanl & Bigten
  /*
  fout<<"     input\n";
  for(int i=0;i<pnt;i++){
    fout<<"  "<<epnt[i]<<"   "<<wgt[i]<<"\n";
  };
  */
  fout<<"\n\n";

  // (material)
  fout<<"material\n";
  for(int i=0;i<nucnum;i++){
    fout<<"        "<<aceid[i]<<"     "<<den[i]<<"\n";
  };
  fout<<"\n";

  // (ace_file)
  fout<<"ace_file  \n";
  for(int i=0;i<nucnum;i++){
    fout<<"  "+dir_acef+"/"<<nucname[i]<<".t0300.dat\n";
  };
  fout<<"\n";

  fout<<"thermal_cutoff\n";
  fout<<"   1e-5\n";
  fout<<"\n";

  // (cal_cond)
  fout<<"calc_cond\n";
  if(slowing_down_calc){
    fout<<"   SLD    ";
    fout<<sd_wgt<<"   ";
    fout<<fehi<<"   1.0  \n";    
  }else{
    fout<<"   NR	0.9996167 \n";
    fout<<"    NoThermalXS\n";
  };
  fout<<"\n";
  fout<<"\n";

  // (bg_xs_set)
  fout<<"bg_xs_set\n";
  //fout<<"   1.0e+10  1e+6  1e+5  1e+4  1e+3  1e+2  1e+1  1.0  0.1\n";
  fout<<"   1.0e+10  0.\n";  
  fout<<"\n";
  fout<<"\n";

  fout<<"max_pl\n";
  fout<<"	5\n";
  fout<<"\n";

  fout<<"edit_xs\n";
  //fout<<"          KRAMXS  Macro\n";
  fout<<"          KRAMXS  Macro CurrentWeightTotalXS\n";    
  //fout<<"          Default\n";
  //fout<<"          -KRAMXS\n";
  //fout<<"          -MATXS\n";
  //fout<<"          -UFG\n";  
  fout<<"\n";

  fout<<"file_add_name\n";
  fout<<"  "<<filenameadd<<"\n";
  fout<<"\n";

  fout<<"end_mg\n";

  
  fout.close();
};

void FMG_Tool::input_generator(string aceid, string nucname, string temp)
{
  ofstream fout;
  fout.open("./input",ios::out);
  if(fout.fail()){
    cout<<"Failed to open the OUTPUT file.\n";
    exit(1);
  };

  ReadEnergyGroupStructure();
  /*
  ifstream fin2;
  fin2.open(file_grp.data(),ios::in);
  if(fin2.fail()){
    cout<<"Failed to open the group structure file.\n";
  };
  int grp;
  fin2>>grp;
  vector<float> ebnd;
  for(int g=0;g<grp+1;g++){
    float tmp;
    fin2>>tmp;
    ebnd.push_back(tmp);
  };
  */

  string tmp1;

  // (ultra_fine_group)
  fout<<"ultra_fine_group\n";
  fout<<" 1.0e+07            10000  EL \n";
  fout<<" 52475.0            56000  EL \n";
  fout<<" 9118.8             12000  EL \n";
  fout<<" 4307.4             12000  EL \n";
  fout<<" 961.12             40000  EL \n";
  fout<<" 130.07             60000  EL \n";
  fout<<" 0.32242            50000  EL \n";
  fout<<" 1.0e-5 \n";
  fout<<"\n";


  // (multi_group)
  fout<<"multi_group\n";
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<grp+1;i++){
    fout<<" "<<ebnd[i];
    if(i%6==5)fout<<"\n";
  };
  fout<<"\n";
  fout<<"\n";

  // (default_spectrum)
  fout<<"default_spectrum\n";
  if(pnt==-1){
    if(wgt_inverse_e){
      fout<<"     1/E\n";
    }else{
      fout<<"     Fission+1/E+Maxwell\n";
      fout<<"	2.0e+7, 1.0e-5, 1.6e+6, 600., 1.0e+6, 0.625\n";
    };
  }else{
    fout<<"     input\n";
    for(int i=0;i<pnt;i++){
      fout<<"  "<<epnt[i]<<"   "<<wgt[i]<<"\n";
    };
  };
  
  fout<<"\n\n";

  // (material)
  fout<<"material\n";
  fout<<"        ";
  fout<<aceid<<"     "<<1.0<<"\n";
  fout<<"\n";

  // (ace_file)
  fout<<"ace_file  \n";
  //fout<<"  /home/roko-de-go/XSDATA/ACEF/JENDL-4/"<<nucname<<".t"<<temp<<".dat\n";
  fout<<"  "+dir_acef+"/"<<nucname<<".t"<<temp<<".dat\n";    
  fout<<"\n";

  fout<<"thermal_cutoff\n";
  fout<<"   1e-5\n";
  fout<<"\n";

  // (cal_cond)
  fout<<"calc_cond\n";
  if(slowing_down_calc){
    fout<<"   SLD    ";
    fout<<sd_wgt<<"   ";
    fout<<fehi<<"   1.0  \n";    
  }else{
    fout<<"   NR	0.9996167 \n";
    fout<<"    NoThermalXS\n";
  };
  fout<<"\n";
  fout<<"\n";

  // (bg_xs_set)
  fout<<"bg_xs_set\n";
  fout<<"   ";
  for(int i=0;i<num_bgxs;i++){
    fout<<bgxs[i]<<" ";
  };
  fout<<"\n";
  fout<<"\n";
  fout<<"\n";

  fout<<"max_pl\n";
  fout<<"	5\n";
  fout<<"\n";

  fout<<"edit_xs\n";
  fout<<"          Default\n";
  fout<<"          -KRAMXS\n";
  fout<<"          -MATXS\n";
  fout<<"          -UFG\n";  
  fout<<"\n";

  fout<<"file_add_name\n"; 
  fout<<"  "<<nucname<<"_"<<temp<<"K\n"; 


  fout<<"\n";

  fout<<"end_mg\n";
  
  fout.close();
};

void FMG_Tool::input_generator_tsl(string aceid1, string aceid2, string nucname, string acefname, string sab_type)
{
  ofstream fout;
  fout.open("./input",ios::out);
  if(fout.fail()){
    cout<<"Failed to open the OUTPUT file.\n";
    exit(1);
  };

  ifstream fin2;
  fin2.open(file_grp.data(),ios::in);
  if(fin2.fail()){
    cout<<"Failed to open the group structure file.\n";
  };
  int grp;
  fin2>>grp;
  vector<float> ebnd;
  for(int g=0;g<grp+1;g++){
    float tmp;
    fin2>>tmp;
    ebnd.push_back(tmp);
  };

  string tmp1;

  // (ultra_fine_group)
  fout<<"ultra_fine_group\n";
  fout<<" 1.0e+07            10000  EL \n";
  fout<<" 52475.0            56000  EL \n";
  fout<<" 9118.8             12000  EL \n";
  fout<<" 4307.4             12000  EL \n";
  fout<<" 961.12             40000  EL \n";
  fout<<" 130.07             60000  EL \n";
  fout<<" 0.32242            50000  EL \n";
  fout<<" 1.0e-5 \n";
  fout<<"\n";

  // (multi_group)
  fout<<"multi_group\n";
  fout.setf(ios::scientific);
  fout.precision(5);
  for(int i=0;i<grp+1;i++){
    fout<<" "<<ebnd[i];
    if(i%6==5)fout<<"\n";
  };
  fout<<"\n";
  fout<<"\n";

  // (default_spectrum)
  fout<<"default_spectrum\n";
  /*
  fout<<"     Fission+1/E+Maxwell\n";
  fout<<"	1.0e+7, 1.0e-5, 1.6e+6, 350, 1.0e+6, 0.1\n";
  */
  fout<<"     1/E\n";  
  fout<<"\n\n";

  // (material)
  fout<<"material\n";
  fout<<"        ";
  fout<<aceid1<<"     "<<aceid2<<"     "<<sab_type<<"    1.0\n";
  fout<<"\n";

  // (ace_file)
  fout<<"ace_file  \n";
  fout<<"  "+dir_acef+"/"<<nucname<<".t0300.dat\n";
  fout<<"  "+dir_acef+"/"<<acefname<<"\n";    
  fout<<"\n";

  fout<<"thermal_cutoff\n";
  fout<<"   4.\n";
  fout<<"\n";

  // (cal_cond)
  fout<<"calc_cond\n";
  fout<<"   NR	0.9996167 \n";
  fout<<"\n";
  fout<<"\n";

  // (bg_xs_set)
  fout<<"bg_xs_set\n";
  fout<<"   1.0e+10\n";
  fout<<"\n";
  fout<<"\n";

  fout<<"max_pl\n";
  fout<<"	5\n";
  fout<<"\n";

  fout<<"edit_xs\n";
  fout<<"          Default\n";
  fout<<"          -KRAMXS\n";
  fout<<"          -MATXS\n";
  fout<<"          -UFG\n";  
  fout<<"\n";

  fout<<"file_add_name\n";
  fout<<"  "<<acefname<<"\n";
  fout<<"\n";
  
  fout<<"end_mg\n";
  
  fout.close();
};

void FMG_Tool::delayed_neutron_data_reader(int famnum, int group, string nucname, string frendy_dir, string cbglib_dir)
{
  MATIDTranslator midt;  
  int matid=midt.ID(nucname);
  string matidname=IntToString(matid);
  int leng=matidname.size();
  matidname=matidname.substr(0,leng-1);

  if(nucname=="Am242m")matidname="95642";
  string filename=frendy_dir+"FMNuChi_"+nucname+"_0300K_"+matidname+".30c.txt";
  string filename_out=cbglib_dir+nucname;

  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  string dummy;  
  for(int i=0;i<9;i++){fin>>dummy;};

  real totnud;
  fin>>totnud;

  if(totnud>0.){

  for(int i=0;i<10;i++){fin>>dummy;};

  vector<real> ai(famnum);
  for(int i=0;i<famnum;i++){
    fin>>ai[i];
  };
  
  for(int i=0;i<10;i++){fin>>dummy;};

  vector<real> lambda(famnum);
  for(int i=0;i<famnum;i++){
    fin>>lambda[i];
  };

  for(int i=0;i<11;i++){fin>>dummy;};

  vector<real> nu_d(group);
  vector< vector<real> > chi_d(famnum);
  for(int g=0;g<group;g++){

    for(int i=0;i<6;i++){fin>>dummy;};
    fin>>nu_d[g];
    for(int i=0;i<3;i++){fin>>dummy;};
    for(int i=0;i<famnum;i++){
      real tmp;
      fin>>tmp;
      chi_d[i].push_back(tmp);
    };

  };

  ofstream fout;
  fout.open(filename_out.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename_out<<"\n";
    exit(0);
  };

  fout<<"  "<<group<<"\n";
  fout<<"  "<<matid<<"\n";
  fout<<"  "<<famnum<<"\n";

  fout.setf(ios::scientific);
  fout.precision(8);
  for(int i=0;i<famnum;i++){
    fout<<ai[i]<<"\n";
  };
  for(int g=0;g<group;g++){
    fout<<" "<<nu_d[g]<<"\n";
  };
  for(int i=0;i<famnum;i++){
    fout<<"   ";
    for(int g=0;g<group;g++){
      fout<<"  "<<chi_d[i][g]<<"\n";
    };
  };
  for(int i=0;i<famnum;i++){
    fout<<" "<<lambda[i]<<"\n";
  };
  fout.close();  

  };

  fin.close();
};

void FMG_Tool::nenergy_generator(int group, string frendy_dir, string cbglib_dir, string nucname)  
{
  MATIDTranslator midt;  
  int matid=midt.ID(nucname);
  string matidname=IntToString(matid);
  int leng=matidname.size();
  matidname=matidname.substr(0,leng-1);
  
  string filename=frendy_dir+"FMMGFlux_"+nucname+"_0300K.txt";

  ifstream fin;
  fin.open(filename.data(),ios::in);
  if(fin.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename<<"\n";
    exit(1);
  };

  ofstream fout;
  string filename_out=cbglib_dir+"/N-ENERGY";
  fout.open(filename_out.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file.\n";
    cout<<"# File name is "<<filename_out<<"\n";
    exit(0);
  };

  string dummy;  
  for(int i=0;i<16;i++){fin>>dummy;};

  vector<real> ebnd(group+1);
  vector<real> flx(group);

  for(int g=0;g<group;g++){
    fin>>dummy;
    fin>>ebnd[g];
    fin>>ebnd[g+1];
    fin>>dummy;
    fin>>flx[g];
    for(int j=0;j<8;j++){fin>>dummy;};
  };

  fout.setf(ios::scientific);
  fout.precision(9);
  fout<<"    "<<group<<"\n";
  for(int g=0;g<group+1;g++){
    fout<<"  "<<ebnd[g]<<"\n";
  };
  for(int g=0;g<group;g++){
    fout<<"    "<<flx[g]<<"\n";
  };
  
  fin.close();
  fout.close();
};
