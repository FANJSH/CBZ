#include "ResonanceCalculationFan.h"
using namespace std;

void ResonanceCalculationInHomoSystemFan::WithSelfShieldingCalculation(int arg_grp, int arg_pl, int arg_nucnum, string arg_libdir, string *arg_filename, int *arg_nuclideID, real *arg_numberdensity, real arg_temperature, bool sig0state){
  if(!sig0state)
  {
    group = arg_grp;
    pl = arg_pl;
    nucnum = arg_nucnum;
    libdir = arg_libdir;
    temperature = arg_temperature;
    NuclideNameSaving.resize(nucnum);
    NuclideIDSaving.resize(nucnum);  
    NumberDensitySaving.resize(nucnum);
    for(int i=0;i<nucnum;i++){
      NuclideNameSaving[i] = arg_filename[i];
      NumberDensitySaving[i] = arg_numberdensity[i];
      NuclideIDSaving[i] = arg_nuclideID[i];  
    }
    BackgroundXS.resize(nucnum,GroupData1D(group));
    EffectiveXS1D.resize(nucnum,vector<GroupData1D>(6,GroupData1D(group)));
    EffectiveTotalXSSaving.resize(nucnum,vector<real>(group));
    FluxWeightedTotalXSSaving.resize(nucnum,GroupData1D(group));
  }

  #if 0
  // Background xs updated or not checking. 
  cout<<"#\t Background xs updating situation check\n";
  cout.precision(7);
  for(int g=0;g<group;g++){
    cout<<setw(8)<<this_medium.GetEnband().get_dat(g)<<"  ";
    for(int i=0;i<nucnum;i++){
      cout<<setw(8)<<BackgroundXS[i].get_dat(g)<<"  ";
    }
    cout<<"\n";
  }
  #endif 

/*
  cout<<"# Background xs check before resonance calculation\n";
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      if(arg_nuclideID[i]==922350){
        cout<<BackgroundXS[i].get_dat(g)<<"\n";
      }
    }
  }
*/

  //XSLibrary xslib(libdir,"N-ENERGY");
  xslib.Initialize(libdir,"N-ENERGY");
  xslib.ReadFile(nucnum,libdir,arg_filename);

  this_medium.PutImax(group);
  this_medium.PutPL(pl,1);
  this_medium.PutNuclide(nucnum,arg_nuclideID,arg_numberdensity);
  this_medium.PutTemperatureForAllNuclide(temperature);
  this_medium.GetEnband().copy(xslib.GetEnband());

  // Retrieve all xs info from CBZLIB
  vector<LibData> NuclideInfo(nucnum);
  for(int i=0;i<nucnum;i++){NuclideInfo[i] = xslib.GetLibData(arg_nuclideID[i]);}

  #if 0
  int index = 0;
  cout<<"#\t Infinite dilution xs check: \n";
  cout.precision(5);
  for(int i=0;i<nucnum;i++){
    cout<<"#\t Nuclide: "<<NuclideIDSaving[i]<<", index:"<<index<<"\n";
    for(int g=0;g<group;g++){
      cout<<setw(8)<<this_medium.GetEnband().get_dat(g)<<"\t";
      cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(sigf).get_dat(g)<<"\t";
      cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(sigc).get_dat(g)<<"\t";
      cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(sigel).get_dat(g)<<"\t";
      cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(siginel).get_dat(g)<<"\t";
      cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g)<<"\n";
    }
    cout<<"\n\n";
    index++;
  }
  #endif 
  printf("%s","# 1. xslib data retrieve: okay.\n");
  this_medium.GetFlux().copy(xslib.GetWtflux());
  
  real tmp1=0.;
  // Use infinite dilution cross-section to define background xs directly at first. 
  // Needs a GroupData1D array to save background xs tempoarily.
  if(!sig0state){ //对某介质首次执行背景截面计算时，需要通过无限稀释截面进行。其它情况跳过此步骤
    printf("%s","#\t Use infinite dilution xs to give background xs at first\n");
    for(int i=0.;i<nucnum;i++){
      for(int g=0.;g<group;g++){
        for(int j=0.;j<nucnum;j++){
          if(j!=i){
            tmp1 += NuclideInfo[j].GetXSData().GetData1d(sigt).get_dat(g) * arg_numberdensity[j];
          }
        }
        tmp1 = tmp1 / arg_numberdensity[i];
        BackgroundXS[i].put_data(g,tmp1);
        tmp1 = 0.;
      }
    }
    //cout<<"# U235 infinite dilution background xs at group 0 is : "<<BackgroundXS[4].get_dat(0)<<"\n";
  }else{
    printf("%s","#\t Skip infinite dilution xs step since already has background xs\n");
  };

  printf("%s","# 2. background xs cal: okay.\n");
  // After obtaining background xs, using self-shielding factor to give effective xs for various reaction. 
  // Therefore, the ralationship between the self-shielding factor and background xs is necessary. 
  // Such info has been stored in the CBZLIB. So retrieve it now. 
  // Example to retreive self-shielding factor from CBZLIB
  // mt=0 : fission
  // mt=1 : capture
  // mt=2 : elastic
  // mt=3 : total (current-weighted)
  // mt=4 : elastic removal
  // mt=5 : inelastic

  //vector<vector<GroupData1D>> EffectiveXS1D(nucnum,vector<GroupData1D>(6,GroupData1D(group))); // 6 is indexs number, corresponding to six kinds of 1D xs. 
  vector<vector<GroupData2D>> EffectiveElasticScatteringMatrixPL(pl+1,vector<GroupData2D>(nucnum, GroupData2D(group,group)));
  vector<vector<GroupData2D>> EffectiveInelasticScatteringMatrixPL(pl+1,vector<GroupData2D>(nucnum, GroupData2D(group,group)));
  vector<GroupData1D> XSError(nucnum,GroupData1D(group));
  GroupData2D EffectiveTotalXS(nucnum,group);  // For updating total xs. 

  int rid=0.;      // fixed
  real rval=0.;   // fixed
  real temp; // temperature [K]
  real sig0=0.; // backgroud XS
  real f_factor=0.; // self-shielding factor
  int iter = 0.;
  real deltasig, deltasig_old = 0;

  // Then, starts calculate effective xs according to background xs and self-sheidling factor. 
  for(deltasig=1;abs(deltasig-deltasig_old)>1e-10;){
    deltasig_old = deltasig;
    for(int i=0;i<nucnum;i++){
      cout<<"#\t Nuclide: "<<NuclideIDSaving[i]<<"\n";
      temp = this_medium.GetNuclideInTurn(i).GetTemperature();
      cout<<"#\t Temp check: "<<temp<<"\n";
      for(int g=0;g<group;g++){
        // It is necessary to check the consistency between input deck pl number and CBZLIB file pl number,
        // Omit this step at here.
        sig0 = BackgroundXS[i].get_dat(g);

        f_factor = NuclideInfo[i].GetFtable().GetF(0,g,rid,rval,temp,sig0);
        EffectiveXS1D[i][0].put_data(g,NuclideInfo[i].GetXSData().GetData1d(sigf).get_dat(g) * f_factor); // sigf
        //cout<<setw(8)<<this_medium.GetEnband().get_dat(g)<<"  "<<NuclideIDSaving[i]<<"  "<<f_factor<<"\t";

        f_factor = NuclideInfo[i].GetFtable().GetF(1,g,rid,rval,temp,sig0);
        EffectiveXS1D[i][1].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigc).get_dat(g) * f_factor)); // sigc
        //cout<<f_factor<<"\t";
        /*
        if(g==56&&NuclideIDSaving[i]==922380){
          cout<<this_medium.GetEnband().get_dat(g)<<"  "<<NuclideInfo[i].GetXSData().GetData1d(sigc).get_dat(g)<<"  "<<f_factor<<"  "<<EffectiveXS1D[i][1].get_dat(g)<<"\n";
        }*/

        f_factor = NuclideInfo[i].GetFtable().GetF(2,g,rid,rval,temp,sig0);
        EffectiveXS1D[i][2].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigel).get_dat(g) * f_factor)); // sigel (1d)
        //cout<<f_factor<<"\t";

        f_factor = NuclideInfo[i].GetFtable().GetF(3,g,rid,rval,temp,sig0);
        EffectiveXS1D[i][3].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigt,0).get_dat(g) * f_factor)); // (sigt,1) current-weighted.
        //cout<<f_factor<<"\n";

        EffectiveXS1D[i][5].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(siginel).get_dat(g))); // siginel  (1d)

        // Next, elastic scattering xs matrix
        for(int pl_iter=0;pl_iter<=pl;pl_iter++){
          for(int g2=0;g2<group;g2++){
            if(g2==g){
              if(g<(group-1)){ // confirm that the group g is not the final group.
                f_factor = NuclideInfo[i].GetFtable().GetF(4,g,rid,rval,temp,sig0); // Elastic removal self-shielding factor
                EffectiveElasticScatteringMatrixPL[pl_iter][i].put_data(g,g2+1,f_factor*NuclideInfo[i].GetXSData().GetData2d(sigel,pl_iter).get_dat(g,g2+1));
              }else{  // If not diagnoal position, set as zero.
              };
              //self-scattering xs.
              if(pl_iter==0){
                tmp1 = EffectiveXS1D[i][2].get_dat(g) - EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2+1); // el-removal - (g->g+1) = el-self-secttering.
                EffectiveElasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,tmp1);
              }else{
                // when pl>0, using 0-th order data to cal f-factor. o-th order effective elastic self-scattering xs has been caled above. 
                f_factor=EffectiveElasticScatteringMatrixPL[0][i].get_dat(g,g2)/NuclideInfo[i].GetXSData().GetData2d(sigel,0).get_dat(g,g2);
                EffectiveElasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,NuclideInfo[i].GetXSData().GetData2d(sigel,pl_iter).get_dat(g,g2)*f_factor); 
              }
            }else if(g2!=(g+1)){ // For all postion that not at diagnoal, set zero. 
              EffectiveElasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,0);
            }
          }
        }

        // Next, inelastic scattering xs matrix (this can be removed actually, just for test)
        for(int pl_iter=0;pl_iter<=pl;pl_iter++){
          for(int g2=0;g2<group;g2++){
            f_factor = NuclideInfo[i].GetFtable().GetF(5,g,rid,rval,temp,sig0);
            EffectiveInelasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,f_factor * NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2));
          }
        }
      }
      cout<<"\n\n";
    }
    #if 0
    // prepare total xs
    tmp1 = 0;
    for(int i=0;i<nucnum;i++){
      for(int g=0;g<group;g++){
        tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g); // fission + capture
        for(int pl_iter=0;pl_iter<=pl;pl_iter++){
          if(pl_iter==0){
          for(int g2=0;g2<group;g2++){
            tmp1 = tmp1 + EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2); // scattering
            //tmp1 = tmp1 + NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2);  // inelastic 
            tmp1 = tmp1 + EffectiveInelasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2);  // inelastic 
            if(i==3&&pl_iter>=2){
              tmp1 = tmp1 + 0;  // Na has no n2n info when pl>=2. 
            }else{
              tmp1 = tmp1 - NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2); // n2n
            }
          }
          }
        }
        EffectiveTotalXS.put_data(i,g,tmp1);
        EffectiveTotalXSSaving[i][g] = tmp1;
      }
      tmp1 = 0.;
    }
    #endif

    #if 0
    // ------------ Another way to get total xs ----------------
      // prepare total xs
      tmp1 = 0;
      for(int i=0;i<nucnum;i++){
        for(int g=0;g<group;g++){
          tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g) + EffectiveXS1D[i][2].get_dat(g) + EffectiveXS1D[i][5].get_dat(g); // fission + capture + elastic(1d) + inelast(1d)
          for(int pl_iter=0;pl_iter<=pl;pl_iter++){
            for(int g2=0;g2<group;g2++){
              //if(i==3&&pl_iter>=2){
              if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||    (arg_nuclideID[i]==280000&&pl_iter>=5)){  // Na, O, Mn, no such data
                tmp1 = tmp1 - 0;  // Na, O, Mn have no n2n info when pl>=2. 
              }else{
                //cout<<i<<" "<<g<<" "<<g2<<" "<<pl_iter<<"\n";
                tmp1 = tmp1 - NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2); // n2n
              }
            } 
          }
          EffectiveTotalXS.put_data(i,g,tmp1);
          EffectiveTotalXSSaving[i][g]=tmp1;
          FluxWeightedTotalXSSaving[i].put_data(g,tmp1);
        }
        tmp1 = 0.;
      }
    #endif

    #if 1
    // ------------ a third way to get total xs ----------------
    // prepare total xs
    tmp1 = 0;
    real tmp12 = 0.;
    for(int i=0;i<nucnum;i++){
      for(int g=0;g<group;g++){
        tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g) + EffectiveXS1D[i][2].get_dat(g) + EffectiveXS1D[i][5].get_dat(g) + NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g); // fission + capture + elastic(1d) + inelast(1d) + n2n(1d)
        EffectiveTotalXS.put_data(i,g,tmp1);
        EffectiveTotalXSSaving[i][g]=tmp1;
        FluxWeightedTotalXSSaving[i].put_data(g,tmp1);
      }
      tmp1 = 0.;
    }
    #endif
    
    if(sig0state){
      printf("%s","#\t （跳过Bondarenko循环）Generate effective crosse-sction with new background cross-section (by Tone's Method) directly\n");
      break;
    }

    // ---------------------------------------------------------
    // get the difference beteen current and previous background xs, and Update background xs
    tmp1 = 0.;
    deltasig = 0.;
    for(int i=0.;i<nucnum;i++){
      for(int g=0.;g<group;g++){
        for(int j=0.;j<nucnum;j++){
          if(j!=i){
            tmp1 += EffectiveTotalXS.get_dat(j,g) * arg_numberdensity[j];
          }
        }
        tmp1 = tmp1/arg_numberdensity[i]; // new background xs of nuclide i in group g. 
        XSError[i].put_data(g,abs(tmp1-BackgroundXS[i].get_dat(g))); // Save the differences on background xs. 
        BackgroundXS[i].put_data(g,tmp1);
        tmp1 = 0.;
      }
      deltasig += XSError[i].get_sum(); // Sum the largest differences of each group together.
    }

    iter++;
    if(iter>10000){
      //cout<<"WARNING: CANNOT GET CONVERGENCE WITHIN 10000 LOOPS"<<"\n";
      break;
    }
  }
  printf("%s" "%d\n","#    Total Iteration times: ",iter);
  printf("%s","# 3. effective xs cal: okay \n");

  #if 0
  cout<<"#\t Micro XS check with in class:\n";
  int index=0;
  cout.precision(6);
  for(int i=0;i<nucnum;i++){
    cout<<"#\t Nuclide: "<<NuclideIDSaving[i]<<"\n";
    for(int g=0;g<group;g++){
      if(g==56&&NuclideIDSaving[i]==922380){
        cout<<setw(8)<<this_medium.GetEnband().get_dat(g)<<"  ";
        cout<<setw(8)<<EffectiveXS1D[i][0].get_dat(g)<<"  ";// sigf
        cout<<setw(8)<<EffectiveXS1D[i][1].get_dat(g)<<"  ";// sigc
        cout<<setw(8)<<EffectiveXS1D[i][2].get_dat(g)<<"  ";// sigel
        cout<<setw(8)<<EffectiveXS1D[i][5].get_dat(g)<<"  ";// siginel
        cout<<setw(8)<<NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g)<<"  "; 
        cout<<setw(8)<<EffectiveXS1D[i][0].get_dat(g)+EffectiveXS1D[i][1].get_dat(g)+EffectiveXS1D[i][2].get_dat(g)+EffectiveXS1D[i][5].get_dat(g)+NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g)<<"  ";
        cout<<setw(8)<<EffectiveTotalXSSaving[i][g]<<"\n";
      }

    }
    cout<<"\n\n";
    index++;
  }
  #endif 

  // Next, Macroscopic data is prepared and save in medium_cbz
  printf("%s","# 4. Macro XS preparation.\n");
  tmp1 = 0.;
  real tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
  for(int g=0;g<group;g++){
    tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = tmp9 = tmp10 = tmp11 = 0.;
    for(int i=0;i<nucnum;i++){
      tmp1 = tmp1 + (EffectiveXS1D[i][3].get_dat(g) * arg_numberdensity[i]);  // current-weighted total
      tmp2 = tmp2 + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g)*EffectiveXS1D[i][0].get_dat(g)*arg_numberdensity[i]); // nuSigf
      tmp3 = tmp3 + (EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g))* arg_numberdensity[i];  // Siga = sigf + sigc
      tmp4 = tmp4 + NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g) * arg_numberdensity[i]; // sign2n
      for(int g2=0;g2<group;g2++){
        tmp6 = tmp6 + (EffectiveElasticScatteringMatrixPL[0][i].get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(siginel,0).get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(sign2n,0).get_dat(g,g2)) * arg_numberdensity[i]; // 0th order sigs
        tmp7 = tmp7 + (EffectiveElasticScatteringMatrixPL[1][i].get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(siginel,1).get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(sign2n,1).get_dat(g,g2)) * arg_numberdensity[i]; // 1th order sigs
        for(int pl_iter=0;pl_iter<=pl;pl_iter++){
          //if(i==3&&pl_iter>=2){
          if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||(arg_nuclideID[i]==280000&&pl_iter>=5)){
            tmp8 = tmp8 + 0;
          }else{
            tmp8 = tmp8 + (NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2) * arg_numberdensity[i]);
          }
        }
      }
    }
    
    // Macro scattering needs to be handled specifically. 
    for(int g2=0;g2<group;g2++){
      for(int pl_iter=0.;pl_iter<=pl;pl_iter++){
        for(int i=0.;i<nucnum;i++){
          if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||(arg_nuclideID[i]==280000&&pl_iter>=5)){
            tmp5 = tmp5 + 
            (EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2) + 
              NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2)) * arg_numberdensity[i]; // Socium has no n2n info when pl>1.
            //medium_cbz.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);
          }else{
            tmp5 = tmp5 + 
            (EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2) + 
             NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2) + 
             NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2)) * arg_numberdensity[i];
            //medium_cbz.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);  // Macro total scattering (el + inel + n2n) 
          }
        }
        this_medium.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);  // Macro total scattering (el + inel + n2n)
        tmp5 = 0.;
      }
      //tmp10 += this_medium.GetSigs(0).get_dat(g,g2);
      tmp10 += this_medium.GetMacxs().GetData2d(sigs,0).get_dat(g,g2); // 0-order Macro total scattering (el + inel + n2n)
    }

    // Transport xs = Siga  + \sum{Sigs(0)} - Sign2n -1/3\sum{sigs(1)}
    tmp9 = tmp1 - 0.3333333333*tmp7;
    //tmp9 = tmp3 + tmp6 - tmp4 - 0.33333333333*tmp7;

    this_medium.GetMacxs().GetData1d(sigt,1).put_data(g,tmp1); // current-weighted
    this_medium.GetMacxs().GetData1d(nusigf).put_data(g,tmp2); // nusigf
    this_medium.GetMacxs().GetData1d(siga).put_data(g,tmp3); // siga
    //medium_cbz.GetMacxs().GetData1d(sigc).put_data(g,tmp12); // sigc
    this_medium.GetMacxs().GetData1d(sign2n).put_data(g,tmp4); // sign2n
    this_medium.GetMacxs().GetData1d(sigtr).put_data(g,tmp9);  // sigtr
    this_medium.GetMacxs().GetData1d(d).put_data(g,(1/(3.*tmp9))); // diffusion coeffi
    this_medium.GetMacxs().GetSigt(0).put_data(g,tmp3+tmp10); // Macro flux-weighted xs. checked, correct. 
  }
  printf("%s","# 4. macroscopic cal and save: [okay].\n");

  // Group data: chi
  // chi could be defined in both 1d and 2d
  // fission spectrum of mixed medium is calculated from the fission production rate
  vector<real> MicroFissionProduction(nucnum);
  real denominator=0.;
  real numerator = 0.;
  tmp1 = 0.;
  tmp2 = 0.;
  for(int g=0;g<group;g++){
    for(int i=0;i<nucnum;i++){ // nu sigf phi N
      denominator = denominator + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g) * EffectiveXS1D[i][0].get_dat(g) * this_medium.GetFlux().get_dat(g) * arg_numberdensity[i]);
    }
  }
  //cout<<"Chi calcul, the deniminator is : "<<denominator<<"\n";
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      tmp1 = tmp1 + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g) * EffectiveXS1D[i][0].get_dat(g) * this_medium.GetFlux().get_dat(g) * arg_numberdensity[i]);
    }
    MicroFissionProduction[i] = tmp1;
    tmp1 = 0.;
  }

  for(int g=0;g<group;g++){
    for(int i=0;i<nucnum;i++){
      numerator += NuclideInfo[i].GetXSData().GetData1d(chi).get_dat(g) * MicroFissionProduction[i];
    }
    this_medium.GetMacxs().GetData1d(chi).put_data(g,numerator/denominator);
    numerator = 0.;
  }
  printf("%s","# 6. chi calculation: okay.\n");
  /*
  // Swap the vector info. 
  vector<LibData>().swap(NuclideInfo);
  for(int i=0;i<nucnum;i++){vector<GroupData1D>().swap(EffectiveXS1D[i]);}
  for(int i=0;i<=pl;i++){vector<GroupData2D>().swap(EffectiveElasticScatteringMatrixPL[i]);}
  vector<real>().swap(MicroFissionProduction);
  printf("%s","# 7. vectors is released.\n");
  */
};


void ResonanceCalculationInHomoSystemFan::AssignInfiniteDilutionCrossSection(int arg_grp, int arg_pl, int arg_nucnum, string arg_libdir, string *arg_filename, int *arg_nuclideID, real *arg_numberdensity, real arg_temperature)
{
  group = arg_grp;
  pl = arg_pl;
  nucnum = arg_nucnum;
  libdir = arg_libdir;
  temperature = arg_temperature;
  NuclideNameSaving.resize(nucnum);
  NuclideIDSaving.resize(nucnum);  
  NumberDensitySaving.resize(nucnum);
  for(int i=0;i<nucnum;i++){
    NuclideNameSaving[i] = arg_filename[i];
    NumberDensitySaving[i] = arg_numberdensity[i];
    NuclideIDSaving[i] = arg_nuclideID[i];  
  }
  BackgroundXS.resize(nucnum,GroupData1D(group));
  EffectiveXS1D.resize(nucnum,vector<GroupData1D>(6,GroupData1D(group)));
  EffectiveTotalXSSaving.resize(nucnum,vector<real>(group));
  FluxWeightedTotalXSSaving.resize(nucnum,GroupData1D(group));
  EffectiveN2NMatrixPL.resize(pl+1,vector<GroupData2D>(nucnum,GroupData2D(group,group)));

  XSLibrary xslib(libdir,"N-ENERGY");
  xslib.ReadFile(nucnum,libdir,arg_filename);

  this_medium.PutImax(group);
  this_medium.PutPL(pl,1);
  this_medium.PutNuclide(nucnum,arg_nuclideID,arg_numberdensity);
  this_medium.PutTemperatureForAllNuclide(temperature);
  this_medium.GetEnband().copy(xslib.GetEnband());

  // Retrieve all xs info from CBZLIB
  vector<LibData> NuclideInfo(nucnum);
  for(int i=0;i<nucnum;i++){NuclideInfo[i] = xslib.GetLibData(arg_nuclideID[i]);}
  printf("%s","# 1. xslib data retrieve: okay.\n");
  this_medium.GetFlux().copy(xslib.GetWtflux());

  real tmp1=0.;

  // After obtaining background xs, using self-shielding factor to give effective xs for various reaction. 
  // Therefore, the ralationship between the self-shielding factor and background xs is necessary. 
  // Such info has been stored in the CBZLIB. So retrieve it now. 
  // Example to retreive self-shielding factor from CBZLIB
  // mt=0 : fission
  // mt=1 : capture
  // mt=2 : elastic
  // mt=3 : total (current-weighted)
  // mt=4 : elastic removal
  // mt=5 : inelastic
  vector<vector<GroupData2D>> EffectiveElasticScatteringMatrixPL(pl+1,vector<GroupData2D>(nucnum, GroupData2D(group,group)));
  vector<vector<GroupData2D>> EffectiveInelasticScatteringMatrixPL(pl+1,vector<GroupData2D>(nucnum, GroupData2D(group,group)));
  vector<GroupData1D> XSError(nucnum,GroupData1D(group));
  GroupData2D EffectiveTotalXS(nucnum,group);  // For updating total xs. 

  // Then, starts calculate effective xs according to background xs and self-sheidling factor. 
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      // It is necessary to check the consistency between input deck pl number and CBZLIB file pl number,
      // Omit this step at here.
      EffectiveXS1D[i][0].put_data(g,NuclideInfo[i].GetXSData().GetData1d(sigf).get_dat(g)); // sigf
      EffectiveXS1D[i][1].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigc).get_dat(g))); // sigc
      EffectiveXS1D[i][2].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigel).get_dat(g))); // sigel (1d)
      EffectiveXS1D[i][3].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(sigt,0).get_dat(g))); // (sigt,1) current-weighted.
      EffectiveXS1D[i][5].put_data(g,(NuclideInfo[i].GetXSData().GetData1d(siginel).get_dat(g))); // siginel  (1d)

      // Next, elastic scattering xs matrix
      for(int pl_iter=0;pl_iter<=pl;pl_iter++){
        for(int g2=0;g2<group;g2++){
          EffectiveElasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,NuclideInfo[i].GetXSData().GetData2d(sigel,pl_iter).get_dat(g,g2));
          EffectiveInelasticScatteringMatrixPL[pl_iter][i].put_data(g,g2,NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2));
        }
      }
    }
  }

  #if 0
  // prepare total xs
  tmp1 = 0;
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g); // fission + capture
      for(int pl_iter=0;pl_iter<=pl;pl_iter++){
        if(pl_iter==0){
        for(int g2=0;g2<group;g2++){
          tmp1 = tmp1 + EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2); // scattering
          //tmp1 = tmp1 + NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2);  // inelastic 
          tmp1 = tmp1 + EffectiveInelasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2);  // inelastic 
          if(i==3&&pl_iter>=2){
            tmp1 = tmp1 + 0;  // Na has no n2n info when pl>=2. 
          }else{
            tmp1 = tmp1 - NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2); // n2n
          }
        }
        }
      }
      EffectiveTotalXS.put_data(i,g,tmp1);
      EffectiveTotalXSSaving[i][g] = tmp1;
    }
    tmp1 = 0.;
  }
  #endif

  #if 0
  // ------------ Another way to get total xs ----------------
  // prepare total xs
  tmp1 = 0;
  real tmp12 = 0.;
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g) + EffectiveXS1D[i][2].get_dat(g) + EffectiveXS1D[i][5].get_dat(g); // fission + capture + elastic(1d) + inelast(1d)
      for(int pl_iter=0;pl_iter<=pl;pl_iter++){
        for(int g2=0;g2<group;g2++){
          if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||(arg_nuclideID[i]==280000&&pl_iter>=5)){  // Na, O, Mn, no such data
            tmp1 = tmp1 - 0;  // Na, O, Mn have no n2n info when pl>=2.
          }else{
            tmp1 = tmp1 + NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2); // n2n
            tmp12 += NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2); // n2n
          }
        } 
        EffectiveN2NMatrixPL[pl_iter][i].put_data(g,tmp12);
        tmp12 = 0;
      }
        EffectiveTotalXS.put_data(i,g,tmp1);
        EffectiveTotalXSSaving[i][g]=tmp1;
        FluxWeightedTotalXSSaving[i].put_data(g,tmp1);
    }
    tmp1 = 0.;
  }
  #endif

  #if 1
  // ------------ a third way to get total xs ----------------
  // prepare total xs
  tmp1 = 0;
  real tmp12 = 0.;
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      tmp1 = EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g) + EffectiveXS1D[i][2].get_dat(g) + EffectiveXS1D[i][5].get_dat(g) + NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g); // fission + capture + elastic(1d) + inelast(1d) + n2n(1d)
      EffectiveTotalXS.put_data(i,g,tmp1);
      EffectiveTotalXSSaving[i][g]=tmp1;
      FluxWeightedTotalXSSaving[i].put_data(g,tmp1);
    }
    tmp1 = 0.;
  }
  #endif
  
  tmp1 = 0.;
  for(int i=0.;i<nucnum;i++){
    for(int g=0.;g<group;g++){
      for(int j=0.;j<nucnum;j++){
        if(j!=i){
          tmp1 += EffectiveTotalXS.get_dat(j,g) * arg_numberdensity[j];
        }
      }
      tmp1 = tmp1/arg_numberdensity[i]; // new background xs of nuclide i in group g. 
      BackgroundXS[i].put_data(g,tmp1); // New background xs is updated. Although infinite dilution xs is assigned initially. 
      tmp1 = 0.;
    }
  }

  // Next, Macroscopic data is prepared and save in medium_cbz
  // NOTE: Both 1D and 2D data are needed. So that 2D data: elastic, inelastic and n2n should be prepared. 
  printf("%s","# 4. Macro XS preparation.\n");
  tmp1 = 0.;
  real tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
  for(int g=0;g<group;g++){
    tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = tmp9 = tmp10 = tmp11 = 0.;
    for(int i=0;i<nucnum;i++){
      tmp1 = tmp1 + (EffectiveXS1D[i][3].get_dat(g) * arg_numberdensity[i]);  // current-weighted total
      tmp2 = tmp2 + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g)*EffectiveXS1D[i][0].get_dat(g)*arg_numberdensity[i]); // nuSigf
      tmp3 = tmp3 + (EffectiveXS1D[i][0].get_dat(g) + EffectiveXS1D[i][1].get_dat(g))* arg_numberdensity[i];  // Siga = sigf + sigc
      tmp4 = tmp4 + NuclideInfo[i].GetXSData().GetData1d(sign2n).get_dat(g) * arg_numberdensity[i]; // sign2n
      for(int g2=0;g2<group;g2++){
        tmp6 = tmp6 + (EffectiveElasticScatteringMatrixPL[0][i].get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(siginel,0).get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(sign2n,0).get_dat(g,g2)) * arg_numberdensity[i]; // 0th order sigs
        tmp7 = tmp7 + (EffectiveElasticScatteringMatrixPL[1][i].get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(siginel,1).get_dat(g,g2) + NuclideInfo[i].GetXSData().GetData2d(sign2n,1).get_dat(g,g2)) * arg_numberdensity[i]; // 1th order sigs
        for(int pl_iter=0;pl_iter<=pl;pl_iter++){
          if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||(arg_nuclideID[i]==280000&&pl_iter>=5)){
            tmp8 = tmp8 + 0;
          }else{
            tmp8 = tmp8 + (NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2) * arg_numberdensity[i]);
          }
        }
      }
    }
    
    // Macro scattering needs to be handled specifically. 
    for(int g2=0;g2<group;g2++){
      for(int pl_iter=0.;pl_iter<=pl;pl_iter++){
        for(int i=0.;i<nucnum;i++){
          if(((arg_nuclideID[i]==110230||arg_nuclideID[i]==80160||arg_nuclideID[i]==250550)&&pl_iter>=2)||(arg_nuclideID[i]==280000&&pl_iter>=5)){
            tmp5 = tmp5 + 
            (EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2) + 
              NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2)) * arg_numberdensity[i]; // Socium has no n2n info when pl>1.
            //medium_cbz.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);
          }else{
            tmp5 = tmp5 + 
            (EffectiveElasticScatteringMatrixPL[pl_iter][i].get_dat(g,g2) + 
             NuclideInfo[i].GetXSData().GetData2d(siginel,pl_iter).get_dat(g,g2) + 
             NuclideInfo[i].GetXSData().GetData2d(sign2n,pl_iter).get_dat(g,g2)) * arg_numberdensity[i];
            //medium_cbz.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);  // Macro total scattering (el + inel + n2n) 
          }
        }
        this_medium.GetMacxs().GetData2d(sigs,pl_iter).put_data(g,g2,tmp5);  // Macro total scattering (el + inel + n2n)
        tmp5 = 0.;
      }
      tmp10 += this_medium.GetSigs(0).get_dat(g,g2);
    }

    // Transport xs = Siga  + \sum{Sigs(0)} - Sign2n -1/3\sum{sigs(1)}
    tmp9 = tmp1 - 0.3333333333*tmp7;
    //tmp9 = tmp3 + tmp6 - tmp4 - 0.33333333333*tmp7;

    this_medium.GetMacxs().GetData1d(sigt,1).put_data(g,tmp1); // current-weighted
    this_medium.GetMacxs().GetData1d(nusigf).put_data(g,tmp2); // nusigf
    this_medium.GetMacxs().GetData1d(siga).put_data(g,tmp3); // siga
    //medium_cbz.GetMacxs().GetData1d(sigc).put_data(g,tmp12); // sigc
    this_medium.GetMacxs().GetData1d(sign2n).put_data(g,tmp4); // sign2n
    this_medium.GetMacxs().GetData1d(sigtr).put_data(g,tmp9);  // sigtr
    this_medium.GetMacxs().GetData1d(d).put_data(g,(1/(3.*tmp9))); // diffusion coeffi
    this_medium.GetMacxs().GetSigt(0).put_data(g,tmp3+tmp10-tmp4); // Macro flux-weighted xs. checked, correct. 
  }
  printf("%s","# 4. macroscopic cal and save: [okay].\n");

  // Group data: chi
  // chi could be defined in both 1d and 2d
  // fission spectrum of mixed medium is calculated from the fission production rate
  vector<real> MicroFissionProduction(nucnum);
  real denominator=0.;
  real numerator = 0.;
  tmp1 = 0.;
  tmp2 = 0.;
  for(int g=0;g<group;g++){
    for(int i=0;i<nucnum;i++){ // nu sigf phi N
      denominator = denominator + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g) * EffectiveXS1D[i][0].get_dat(g) * this_medium.GetFlux().get_dat(g) * arg_numberdensity[i]);
    }
  }
  //cout<<"Chi calcul, the deniminator is : "<<denominator<<"\n";
  for(int i=0;i<nucnum;i++){
    for(int g=0;g<group;g++){
      tmp1 = tmp1 + (NuclideInfo[i].GetXSData().GetData1d(nu).get_dat(g) * EffectiveXS1D[i][0].get_dat(g) * this_medium.GetFlux().get_dat(g) * arg_numberdensity[i]);
    }
    MicroFissionProduction[i] = tmp1;
    tmp1 = 0.;
  }

  for(int g=0;g<group;g++){
    for(int i=0;i<nucnum;i++){
      numerator += NuclideInfo[i].GetXSData().GetData1d(chi).get_dat(g) * MicroFissionProduction[i];
    }
    this_medium.GetMacxs().GetData1d(chi).put_data(g,numerator/denominator);
    numerator = 0.;
  }
  printf("%s","# 6. chi calculation: okay.\n");
  /*
  // Swap the vector info. 
  vector<LibData>().swap(NuclideInfo);
  for(int i=0;i<nucnum;i++){vector<GroupData1D>().swap(EffectiveXS1D[i]);}
  for(int i=0;i<=pl;i++){vector<GroupData2D>().swap(EffectiveElasticScatteringMatrixPL[i]);}
  vector<real>().swap(MicroFissionProduction);
  printf("%s","# 7. vectors is released.\n");
  */
};