void MulticellBurner::SensitivityCalculationKeffEOC(Burnup &bu,int med_normalize)
{
  if(!corrector_calc){
    cout<<"# Error in MulticellBurner::SensitivityCalculationKeffEOC.\n";
    cout<<"# Predictor-corrector method should be used.\n";
    exit(0);
  };

  bool isotropic_approx=false; // Isotropic approximation in perturbation integral

  int ssv=40;  // (sub-sub-time step division)

  wc_gpt_wpc=0.6; // !!!

  PreCalculation(bu);

  // +++ Array setting for adjoint calculation +++
  pow_adj.resize(burn_step);
  adj_nuc.resize(burn_step);
  macxs.resize(burn_step);

  pow_adj_c.resize(burn_step); // !!!
  adj_nuc_c.resize(burn_step); // !!!
  macxs_p.resize(burn_step); // !!!

  bilinear_flx.resize(burn_step);
  bilinear_flx_c.resize(burn_step); // !!!
  bilinear_aflx.resize(burn_step);
  bilinear_aflx_c.resize(burn_step); // !!!
  chi_gpt_fwdflx.resize(burn_step);
  chi_gpt_fwdflx_c.resize(burn_step); // !!!
  keff_p.resize(burn_step+1);

  for(int i=0;i<burn_step;i++){
    int sub_step=sub_step_list[i];

    pow_adj[i].resize(sub_step);
    pow_adj_c[i].resize(sub_step); // !!!
    adj_nuc[i].resize(sub_step);
    adj_nuc_c[i].resize(sub_step); // !!!
    for(int j=0;j<sub_step;j++){
      adj_nuc[i][j].resize(mednum_fuel);
      adj_nuc_c[i][j].resize(mednum_fuel); // !!!
      for(int k=0;k<mednum_fuel;k++){
        adj_nuc[i][j][k].put_imax(nucn); 
        adj_nuc_c[i][j][k].put_imax(nucn); // !!! 
      };
    };

    bilinear_flx[i].resize(mednum_fuel);
    bilinear_flx_c[i].resize(mednum_fuel); // !!!
    bilinear_aflx[i].resize(mednum_fuel);
    bilinear_aflx_c[i].resize(mednum_fuel); // !!!
    chi_gpt_fwdflx[i].resize(mednum_fuel);
    chi_gpt_fwdflx_c[i].resize(mednum_fuel); // !!!
    for(int k=0;k<mednum_fuel;k++){
      bilinear_flx[i][k].put_imax(group);
      bilinear_flx_c[i][k].put_imax(group); // !!!
      bilinear_aflx[i][k].put_imax(group);
      bilinear_aflx_c[i][k].put_imax(group); // !!!
      chi_gpt_fwdflx[i][k].put_imax(group);
      chi_gpt_fwdflx_c[i][k].put_imax(group); // !!!
    };

    macxs[i].resize(mednum_fuel);
    macxs_p[i].resize(mednum_fuel); // !!!
    for(int j=0;j<mednum_fuel;j++){
      macxs[i][j].Init("MacroCrossSection");
      macxs_p[i][j].Init("MacroCrossSection"); // !!!
      macxs[i][j].PutGrp(107);
      macxs_p[i][j].PutGrp(107);
    };
  };

  // +++ Forward burnup calculation +++

  ForwardCalculation(bu,med_normalize,true);

  // +++ Integrating forward number density

  IntegratingForwardNumberDensity(fwd_nuc,fwd_nuc_int);
  IntegratingForwardNumberDensity(fwd_nuc_p,fwd_nuc_p_int);

  // +++ Adjoint burnup calculation +++

  GeneralOption opta;
  opta.PutAdjointCal();

  // !!!

  // +++ EOC calculation
  GeneralOption opt;

  // +++ Pre-calculation of Dancoff factor +++
  SelfShieldingCalculator ssc;  
  ssc.ThreeRegionDancoffMethod(xslib,sys,med[0],med[med_clad],med[med_water],true);
  GroupData1D dancoff=ssc.GetDancoff(0);
  GroupData1D bell(group);
  for(int i=0;i<group;i++){bell.put_data(i,1.2);};

  vector<Medium> med_tmp(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    PutNuclideDataToMedium(fwd_nuc[burn_step][0][i],0);
    opc.GiveInfiniteDillutionCrossSection(med[0],xslib);
    opc.CalSelfShieldingWithDancoffCorrection(med[0],xslib,fuel_r*2.,bell,dancoff);
    opc.CalThermalScatteringMatrix(med[0],xslib,3.93);
    med[0].CalMacroFromMicro();
    med_tmp[i]=med[0];
  };


  MECSystem lata(group,mednum);
  //lata.NoPrint();
  lata.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum_fuel;i++){
    lata.AddMedium(med_tmp[i]);
  };
  lata.AddMedium(med[med_clad]);
  lata.AddMedium(med[med_water]);
  lata.PutRegMed(region_medium);
  lata.PutGeneralOption(opta);
  lata.PutPL(0);
  lata.NoCMRAcceleration();
  lata.PutWriteFlux();
  real k_adj=lata.CalIgen();

  MECSystem lat(group,mednum);
  //lata.NoPrint();
  lat.PutTrajectorySet(&sys_f);
  for(int i=0;i<mednum_fuel;i++){
    lat.AddMedium(med_tmp[i]);
  };
  lat.AddMedium(med[med_clad]);
  lat.AddMedium(med[med_water]);
  lat.PutRegMed(region_medium);
  lat.PutGeneralOption(opt);
  lat.PutPL(0);
  lat.NoCMRAcceleration();
  lat.PutWriteFlux();
  real k_fwd=lat.CalIgen();

  // (Direct term)
  int nucnum=0;
  int nnmax=med[0].GetNucnum();
  int *nucid=new int[nnmax];
  for(int i=0;i<med[0].GetNucnum();i++){
    if(med[0].GetNuclideInTurn(i).GetGrp()!=-1){
      int id=med[0].GetNuclideInTurn(i).GetMatnum();
      nucid[nucnum++]=id;
    };
  };
  SensitivityData sns_dir=lata.CalSensitivityNew(&lat,k_fwd,nucnum,nucid);
  delete [] nucid;

  sns_dir.PutName("mburner","k_EOC","unknown");
  sns_dir.WriteFile("./","sns.k_EOC_dir");

  real response=k_fwd;

  cout.setf(ios::scientific);
  cout.precision(5);
  cout<<"# Response (k-inf) : "<<response<<"\n";

  // +++ Adjoint calculation
  // (MEC)
  MECSystem lat_gpt(group,mednum);

  lat_gpt.NoPrint();
  lat_gpt.PutTrajectorySet(&sys_f);
  int totsn=lat_gpt.GetQuad().GetSN();

  for(int i=0;i<mednum_fuel;i++){
    lat_gpt.AddMedium(med[0]);
    lat_gpt.GetMedium(i).NuclideClear();
  };
  lat_gpt.AddMedium(med[med_clad]);
  lat_gpt.AddMedium(med[med_water]);
  lat_gpt.PutRegMed(region_medium);
  lat_gpt.PutGeneralOption(opta);
  lat_gpt.PutPL(0);

  lat_gpt.NoCMRAcceleration();
  //lat_gpt.NoTransportApprox();
  lat_gpt.PutWriteFlux();
  lat_gpt.SetArray();

  
  vector<GroupData1D> gpt_flx(totm);
  vector< vector<GroupData1D> > gpt_aflx(totm);
  for(int i=0;i<totm;i++){
    gpt_flx[i].put_imax(group);
    gpt_aflx[i].resize(group);
    for(int g=0;g<group;g++){
      gpt_aflx[i][g].put_imax(totsn);
    };
  };

  vector<GroupData1D> gpt_src(mednum_fuel);
  vector<GroupData1D> gpt_src2(mednum_fuel);
  for(int i=0;i<mednum_fuel;i++){
    gpt_src[i].put_imax(group);
    gpt_src2[i].put_imax(group);
  };

  vector<GroupData1D> adj_nuc_e(mednum_fuel);
  vector<GroupData1D> adj_nuc_e_c(mednum_fuel); // for corrector calc

  // +++ Final condition for adjoint number density
  for(int i=0;i<mednum_fuel;i++){
    adj_nuc_e[i].put_imax(nucn);
    adj_nuc_e[i].set_zero();
    for(int j=0;j<nucn;j++){
      if(med_tmp[i].GetNuclideInTurn(j).GetGrp()!=-1){
        real org=med_tmp[i].GetNuclideInTurn(j).GetDensity();
        if(org>1e-20){
          lat.GetMedium(i).CalMacroFromMicro();
          lata.GetMedium(i).CalMacroFromMicro();
          lat.GetMedium(i).GetNuclideInTurn(j).PutDensity(org*1.01);
          lat.GetMedium(i).CalMacroFromMicro();
          real dk=lata.CalReactivity(&lat,k_adj,k_fwd,false)*k_adj*k_fwd;
          lat.GetMedium(i)=med_tmp[i];
          //int mmt=med[0].GetNuclideInTurn(j).GetMatnum();
          //cout<<mmt<<" "<<midt.Name(mmt)<<" "<<dk*100.<<"\n";
          real src=0.;
          if(org!=0.)src=dk/(org*0.01);
          adj_nuc_e[i].put_data(j,src);
	};
      };
    };
    med_tmp[i].NuclideClear();
  };

  GroupData2D trmat_flxindep=bu.GetTrmatFlxInDep();
  for(int st=burn_step-1;st>=0;st--){

    real power_density=power_density_list[st];
    int sub_step=sub_step_list[st];

    cout<<"#   Adjoint calculation step : "<<st<<"\n";

    adj_nuc_e_c=adj_nuc_e; // to store initial adjoint for following predictor calculation

    // +++ CORRECTOR PART ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real keff_step=keff[st];
    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross section data
      // is updated in the preceding lines in MEC
    };

    // (PJI)
    /*
    if(pij_storing){
      for(int g=0;g<group;g++){
        lat_gpt.GetPij(g).copy(pij_store[st][1][g]);
      };
    }else{
      lat_gpt.PutPij();
    };
    */
    // (MEC)
    // no processing    

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g[st][m][i],xsc_1g[st][m][i],xsn2n_1g[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc[st][j][m]+fwd_nuc[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_normalize];
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);
  
        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor[st][j];
		};
	      };
	    };
	  };
        }else{
	  if(m==med_normalize){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc[st][j][m][g])*delt[st][j]*power_factor[st][j]/vol_med[m];
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc[st][j][m]))*delt[st][j]*power_factor[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    // ++++
    int sntot=lat_gpt.GetQuad().GetSN();
    int pdiv=lat_gpt.GetPolarAngleDivision();
    int sn_per_pang=sntot/pdiv;

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh[st][mm]));
          real tmp=gpt_flx[mm]*macxs[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
              real tmp=omega*volaflx_mesh[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // +++ CORRECTOR CALCULATION END ++++++++++++++++++++++++++++++++++++++++

    swap(adj_nuc_e, adj_nuc_e_c);
    swap(pow_adj[st], pow_adj_c[st]);
    swap(adj_nuc[st], adj_nuc_c[st]);
    swap(bilinear_flx[st], bilinear_flx_c[st]);
    swap(bilinear_aflx[st], bilinear_aflx_c[st]);
    swap(chi_gpt_fwdflx[st], chi_gpt_fwdflx_c[st]);
    
    // +++ PREDICTOR CALCULATION +++++++++++++++++++++++++++++++++++++++++++++

    for(int m=0;m<mednum_fuel;m++){
      lat_gpt.GetMed(m).GetMacxs().DataCopyPL(macxs_p[st][m],0);
      lat_gpt.GetMed(m).TransportApproximation();
      // (transport approximation) 
      // This procedure is necessary because macro cross sectio data
      // is updated in the preceding lines in MEC
    };

    // (PJI)
    /*
    if(pij_storing){
      for(int g=0;g<group;g++){
        lat_gpt.GetPij(g).copy(pij_store[st][0][g]);
      };
    }else{
      lat_gpt.PutPij();
    };
    */
    // (MEC)
    // no processing
    
    // jump condition for generalized adjoint at the end of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
              tmp2+=xsa*bilinear_flx_c[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx_c[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step; // keff should be one at corrector step
          };
          adj_nuc_e[m].add_data(k,(tmp2-tmp3)*wc_gpt_wpc/(1.-wc_gpt_wpc));
	};
      };
    };

    keff_step=keff_p[st];

    for(int j=sub_step-1;j>=0;j--){

      // (adjoint number density check)
      bool zero_adj=true;
      for(int m=0;m<mednum_fuel;m++){
        for(int n=0;n<nucn;n++){
          if(fabs(adj_nuc_e[m].get_dat(n))>1e-20){
	    zero_adj=false;
	    break;
	  };
        };
      };
      if(zero_adj){
        cout<<"# Error in MulticellBurner::SensitivityCalculationPC.\n";
        cout<<"# All adjoint number density is zero.\n";
        exit(0);
      };

      pow_adj[st][j]=0.;
      for(int m=0;m<mednum_fuel;m++){
        for(int i=0;i<nucn;i++){
	  int id=med[0].GetNuclideInTurn(i).GetID();
	  bu.PutNuclideData(i,id,0.,xsf_1g_p[st][m][i],xsc_1g_p[st][m][i],xsn2n_1g_p[st][m][i]);
	};
	bu.CalTransitionMatrixFluxDependentPart();
        GroupData2D mmat1=bu.GetTrmatFlxDep()*(total_flux_p[st][j][m]*1e-24);
	GroupData2D mmat2=trmat_flxindep+mmat1;
	mmat2.Transposition();

	GroupData1D ttt2=adj_nuc_e[m];
        adj_nuc[st][j][m]=ttt2*(0.5/ssv);
        vector<GroupData1D> ans(ssv);
        mmat2.MultiStepCalc(ttt2,ans,delt[st][j],ssv);
        for(int k=0;k<ssv-1;k++){
          adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[k]/ssv;
        };
        adj_nuc[st][j][m]=adj_nuc[st][j][m]+ans[ssv-1]*(0.5/ssv);
        adj_nuc_e[m]=ans[ssv-1];
        pow_adj[st][j]+=adj_nuc[st][j][m]*(mmat1*(fwd_nuc_p[st][j][m]+fwd_nuc_p[st][j+1][m]))*0.5*delt[st][j];
      };
      if(!input_flux_level){
	real factor=1./power_density;
	if(power_density==0.)factor=1e10;
        pow_adj[st][j]*=factor;
      }else{
        pow_adj[st][j]/=flux_level_list[st]*vol_med[med_normalize];
      };

      // (Jump condition by adjoint power)
      if(!input_flux_level){
        for(int m=0;m<mednum_fuel;m++){
	  if(med_normalize==-1||med_normalize==m){
	    for(int k=0;k<nucn;k++){
	      real xsf1g=xsf_1g_p[st][m][k];
	      if(xsf1g>0.){
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
                real tmp1=xsf1g*total_flux_p[st][j][m]*vol_med[m]
                   *bu.GetReactionEnergyData().GetFissionEnergy(id)*pow_adj[st][j];
                adj_nuc_e[m].add_data(k,-tmp1);
	      };
	    };
	  };
	};
      };

    }; // sub-step loop end

    // Generalized adjoint flux calculation
    // (source calculation)
    for(int m=0;m<mednum_fuel;m++){
      for(int g=0;g<group;g++){

        for(int k=0;k<nucn;k++){
          if(nuclide_info[k]!=0){
            PutMicroXSDataToMedium(mic_sigf[st][m][k],mic_sigc[st][m][k],mic_sign2n[st][m][k],k,0,nuclide_info[k]);
          };
        };
        GroupData2D dmdf=bu.CaldTMatdFlux(med[0],g);

        real tmp=0.;
        if(!input_flux_level){
          // (power term)
	  if(med_normalize==-1||med_normalize==m){
            for(int k=0;k<nucn;k++){
  	      if(nuclide_info[k]==1){ // fissile
	        int id=med[0].GetNuclideInTurn(k).GetMatnum();
	        for(int j=sub_step-1;j>=0;j--){
	          tmp+=mic_sigf[st][m][k].get_dat(g)*fwd_nuc_p[st][j][m].get_dat(k)
		    *bu.GetReactionEnergyData().GetFissionEnergy(id)
                    *pow_adj[st][j]*power_factor_p[st][j];
		};
	      };
	    };
	  };
        }else{
	  if(m==med_normalize){
	    for(int j=sub_step-1;j>=0;j--){
  	      tmp+=pow_adj[st][j]*power_factor_p[st][j];
	    };
	  };
	};
        // (number density term)
	real tmp2=0.;
	for(int j=sub_step-1;j>=0;j--){
          //tmp2+=(adj_nuc[st][j][m]*dmdf_nuc_p[st][j][m][g])*delt[st][j]*power_factor_p[st][j]/vol_med[m];
          tmp2+=(adj_nuc[st][j][m]*(dmdf*fwd_nuc_p[st][j][m]))*delt[st][j]*power_factor_p[st][j]/vol_med[m];
	};
        tmp-=tmp2;

        if(tmp>0.){
          gpt_src[m].put_data(g,tmp);
          gpt_src2[m].put_data(g,0.);
        }else{
	  gpt_src[m].put_data(g,0.);
	  gpt_src2[m].put_data(g,-tmp);
	};
      };
    };

    // (calculation with positive source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=lat_gpt.GetAFlux(rr,g);
      };
    };
    // (calculation with negative source)
    lat_gpt.SetZeroScalarFlux();
    lat_gpt.SetZeroScatSrc();
    for(int m=0;m<totm;m++){
      if(region_medium[m]<mednum_fuel){
        lat_gpt.PutIsotropicSourcePerUnitVolume(m,gpt_src2[region_medium[m]]);
      };
    };
    lat_gpt.CalGPT_MEC(keff_step,1e-4,100);
    for(int rr=0;rr<totm;rr++){
      gpt_flx[rr]=gpt_flx[rr]-lat_gpt.GetMesh(rr).GetFlux();
      for(int g=0;g<group;g++){
        gpt_aflx[rr][g]=gpt_aflx[rr][g]-lat_gpt.GetAFlux(rr,g);
      };
    };

    for(int m=0;m<mednum_fuel;m++){
      bilinear_flx[st][m].set_zero();
      bilinear_aflx[st][m].set_zero();
      chi_gpt_fwdflx[st][m].set_zero();
      for(int mm=0;mm<totm;mm++){
        if(region_medium[mm]==m){
          bilinear_flx[st][m]=bilinear_flx[st][m]+(gpt_flx[mm].mult(volflx_mesh_p[st][mm]));
          real tmp=gpt_flx[mm]*macxs_p[st][m].GetData1d(chi);
          chi_gpt_fwdflx[st][m]=chi_gpt_fwdflx[st][m]+(volflx_mesh_p[st][mm]*tmp);

          for(int sn=0;sn<totsn;sn++){
            real omega=lat_gpt.GetQuad().GetOmega(sn); 
            int sn2=sn;
	    int tmp=sn%sn_per_pang;
	    if(tmp<sn_per_pang/2){
              sn2=sn+sn_per_pang/2;
	    }else{
	      sn2=sn-sn_per_pang/2;
	    };
	    for(int g=0;g<group;g++){
              real tmp=omega*volaflx_mesh_p[st][mm][g].get_dat(sn)*gpt_aflx[mm][g].get_dat(sn2);
              bilinear_aflx[st][m].add_data(g,tmp);
	    };
	  };

        };
      };
      bilinear_aflx[st][m]=bilinear_aflx[st][m]*PI4;
    };

    // jump condition for generalized adjoint at the beginning of step
    for(int m=0;m<mednum_fuel;m++){
      for(int k=0;k<nucn;k++){
        real xsf1g=xsf_1g_p[st][m][k];

        if(nuclide_info[k]!=0){
          // (absorption or total term)
          real tmp2=0.;
          for(int g=0;g<group;g++){
            real xsa=mic_sigc[st][m][k].get_dat(g);
            if(xsf1g>0.)xsa+=mic_sigf[st][m][k].get_dat(g);
	    if(isotropic_approx){
	      tmp2+=xsa*bilinear_flx[st][m].get_dat(g); // absorption
	    }else{
  	      tmp2+=xsa*bilinear_aflx[st][m].get_dat(g); // total
	    };
	  };
	  // (yield term)
	  real tmp3=0.;
	  if(xsf1g>0.){
	    for(int mm=0;mm<totm;mm++){
	      if(region_medium[mm]==m){
	        real fsrc=0.;
	        for(int g=0;g<group;g++){
	          fsrc+=volflx_mesh_p[st][mm].get_dat(g)*mic_sigf[st][m][k].get_dat(g)
	          *xslib.GetLibData(med[0].GetNuclideInTurn(k).GetMatnum()).GetXSData().GetData1d(nu).get_dat(g);     
     	        };
	        for(int g=0;g<group;g++){
	          tmp3+=gpt_flx[mm].get_dat(g)*fsrc*macxs_p[st][m].GetData1d(chi).get_dat(g);
	        };
	      };
	    };
  	    tmp3/=keff_step;
          };
          adj_nuc_e[m].add_data(k,tmp2-tmp3);
	};

      };
    };

    // +++ PREDICTOR CALCULATION END +++++++++++++++++++++++++++++++++++++++++

    // Averaging of adjoint number density
    for(int m=0;m<mednum_fuel;m++){
      adj_nuc_e[m]=adj_nuc_e[m]*(1.-wc_gpt_wpc)+adj_nuc_e_c[m]*wc_gpt_wpc;
    };

  }; // the end ot step

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Sensitivity Printing

  vector< vector< vector<real> > > nadj_dm_nfwd;

  SensitivityData sns;
  sns.PutName("dummy","dummy","dummy");
  sns.PutValue(response);
  sns.PutGroup(group);
  sns.GetEnband().copy(med[0].GetEnband());

  GroupData1D sns1d(group);
  for(int i=0;i<nucn;i++){

    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    int rmax=3;
    if(matnum<900000)rmax=2;
    if(nuclide_info[i]==0)rmax=0; // no cross section data
    for(int r=0;r<rmax;r++){

      enum xstype sigxx=sigc;
      int mt=102;
      int bc_channel=1;
      if(r==1){
	sigxx=sign2n;
	mt=16;
	bc_channel=2;
      }else if(r==2){
        sigxx=sigf;
        mt=18;
        bc_channel=0;
      };

      // +++ predictor
      GroupData1D sns1d_p(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc_p,adj_nuc,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh_p[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc_p[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor_p[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor_p[k][l]*vol_med[m];
                sum-=pow_adj[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][0][m].get_dat(i);
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff_p[k];
	      sum-=tmp*chi_gpt_fwdflx[k][m].get_dat(j);
	    };
	  };
        };
        sns1d_p.put_data(j,sum/response);
      };

      // +++ corrector
      GroupData1D sns1d_c(group);

      CalNadjDMNfwd(i,matnum,bc_channel,bu,fwd_nuc,adj_nuc_c,nadj_dm_nfwd);

      for(int j=0;j<group;j++){
        real sum=0.;
        for(int k=0;k<burn_step;k++){
          int sub_step=sub_step_list[k];
	  for(int m=0;m<mednum_fuel;m++){
            real xs=0.;
            if(r==0)xs=mic_sigc[k][m][i].get_dat(j);
            if(r==1)xs=xslib.GetLibData(matnum).GetXSData().GetData1d(sign2n).get_dat(j);
            if(r==2)xs=mic_sigf[k][m][i].get_dat(j);
            real flx=0.;
	    for(int l=0;l<totm;l++){
	      if(region_medium[l]==m){
		flx+=volflx_mesh[k][l].get_dat(j);
	      };
	    };
	    flx/=vol_med[m];
            for(int l=0;l<sub_step;l++){
              real den=fwd_nuc[k][l][m].get_dat(i);
              // --- Number density term
              real dsig=xs*(flx*power_factor[k][l])*1e-24;
              sum+=dsig*nadj_dm_nfwd[k][l][m];
              // --- Power normalization term (fission case)
              if(sigxx==sigf&&!input_flux_level&&(med_normalize==-1||med_normalize==m)){
                int iid=med[0].GetNuclideInTurn(i).GetMatnum();
                real tmp=flx*power_factor[k][l]*vol_med[m];
                sum-=pow_adj_c[k][l]*tmp*xs*den*bu.GetReactionEnergyData().GetFissionEnergy(iid);
              };
	    };
            // --- flux term [(n,2n) reaction is not well treated yet.]
            real den0=fwd_nuc_p[k][sub_step][m].get_dat(i); // !!! CAUTION
            real nu_value=0.;
  	    if(sigxx==sigf)nu_value=xslib.GetLibData(med[0].GetNuclideInTurn(i).GetMatnum()).GetXSData().GetData1d(nu).get_dat(j);
	    if(isotropic_approx){
              sum+=den0*xs*bilinear_flx_c[k][m].get_dat(j); // (absorption term)
	    }else{
              sum+=den0*xs*bilinear_aflx_c[k][m].get_dat(j); // (total term)
	    };
            // (yield term)
	    if(sigxx==sigf){
	      real tmp=den0*xs*nu_value/keff[k];
	      sum-=tmp*chi_gpt_fwdflx_c[k][m].get_dat(j);
	    };
	  };
	};
        sns1d_c.put_data(j,sum/response);
      };

      sns1d=sns1d_p*(1.-wc_gpt_wpc)+sns1d_c*wc_gpt_wpc;
      sns.PutSensitivity1D(matnum,mt,sns1d);

    };
  }; // end of nuclide loop for cross section sensitivities

  // +++ For fission yield +++++++++++++++++++++++++++++++++++++++++++++++++++++
  int idfisn=21;//kawamoto
  int idfisorg[]={
    922340,922350,922360,922370,922380,
    932370,932390,
    942380,942390,942400,942410,942420,
    952410,952420,952421,952430,
    962420,962430,962440,962450,962460,
  };

  for(int ii=0;ii<idfisn;ii++){
    int idfis=idfisorg[ii];
    int pos0=bu.SearchNuclide(idfis);
    int nuct=bu.GetBC().GetNdivFission(idfis);
    for(int i=0;i<nuct;i++){
      int id=bu.GetBC().GetNextIDFission(idfis,i);
      real rat=bu.GetBC().GetRatioFission(idfis,i);
      int pos=bu.SearchNuclide(id);
      if(pos!=-1){

        vector<real> snsval(2,0.);
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
        for(int jj=0;jj<2;jj++){ // predictor & corrector

          real val=0.;
          for(int k=0;k<burn_step;k++){
            int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
              real rr=xsf_1g[k][m][pos0]*rat;
              for(int l=0;l<sub_step;l++){
                val+=(adj_nuc[k][l][m].get_dat(pos)*rr*total_flux[k][l][m]*1e-24*fwd_nuc[k][l][m].get_dat(pos0))*delt[k][l];
	      };
	    };
	  };
  	  snsval[jj]=val/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

        }; // end of loop-jj
        swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

        real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
        sns.PutSensitivity0D(id,18000000+idfis,snstot);
        //sns.PutSensitivity0D(id,18000000+idfis,val/response);

      };
    };
  };

  // +++ For Half-life +++++++++++++++++++++++++++++++++++++++++++
  // (For decay heat sensitivity, direct term should be taken into account)
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      decay_c*=-0.01; // dT=0.01T -> dlamba=-0.01 lambda

      vector<real> snsval(2,0.);
      swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
      for(int jj=0;jj<2;jj++){ // predictor & corrector

      real sum=0.;
      for(int k=0;k<burn_step;k++){
        int sub_step=sub_step_list[k];
        for(int m=0;m<mednum_fuel;m++){
          for(int l=0;l<sub_step;l++){
            real den=fwd_nuc[k][l][m].get_dat(i);
            real val=-decay_c*den*adj_nuc[k][l][m].get_dat(i);
            int tmp=bu.GetBC().GetNdivDecay(matnum);
            for(int j=0;j<tmp;j++){
              int id2=bu.GetBC().GetNextIDDecay(matnum,j);
              int pos=bu.SearchNuclide(id2);
              if(pos!=-1){
   	        real rat=bu.GetBC().GetRatioDecay(matnum,j);
	        val+=rat*decay_c*den*adj_nuc[k][l][m].get_dat(pos);
	      };
	    };
	    val*=delt[k][l];
	    sum+=val;
	  };
	};
      };
      sum*=100.;// because dT=0.01T

      snsval[jj]=sum/response;

      if(jj==0){
        swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
        swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
      };

      }; // end of loop-jj
      swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c

      real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
      sns.PutSensitivity0D(matnum,8888,snstot);
      //sns.PutSensitivity0D(matnum,8888,sum/response);
    };
  };

  // +++ For Branching Ratio +++++++++++++++++++++++++++++++++++++++++++
  for(int i=0;i<nucn;i++){
    int matnum=med[0].GetNuclideInTurn(i).GetMatnum();
    real decay_c=bu.GetDecayConstant(matnum);
    if(decay_c!=0.){
      int channel=bu.GetBC().GetNdivDecay(matnum);
      if(channel>1){

	vector<real> ratio;
	ratio.resize(channel);
	for(int j=0;j<channel;j++){
	  ratio[j]=bu.GetBC().GetRatioDecay(matnum,j);
	};
	vector<real> sns_tmp1;
	sns_tmp1.resize(channel);

	for(int j=0;j<channel;j++){

          vector<real> snsval(2,0.);
          swap(fwd_nuc_p, fwd_nuc); // fwd_nuc_p -> fwd_nuc
          for(int jj=0;jj<2;jj++){ // predictor & corrector

	  real sum=0.;
	  real rat=bu.GetBC().GetRatioDecay(matnum,j);
	  rat*=0.01;
	  real time=0.;
	  for(int k=0;k<burn_step;k++){
	    int sub_step=sub_step_list[k];
            for(int m=0;m<mednum_fuel;m++){
	      for(int l=0;l<sub_step;l++){
	        real den=fwd_nuc[k][l][m].get_dat(i);
	        real val=0.;
	        int id2=bu.GetBC().GetNextIDDecay(matnum,j);
	        int pos=bu.SearchNuclide(id2);
	        if(pos!=-1){
	  	  val+=adj_nuc[k][l][m].get_dat(pos)*decay_c*rat*den;
	        };
	        time+=delt[k][l];
	        val*=delt[k][l];
	        sum+=val;
	      };
	    };
	  };
	  sum*=100.;// because dr=0.01r
          snsval[jj]=sum/response;

          if(jj==0){
            swap(fwd_nuc_p, fwd_nuc); // fwd_nuc -> fwd_nuc_p
            swap(adj_nuc_c, adj_nuc); // adj_nuc_c -> adj_nuc
          };

          }; // end of loop-jj
          swap(adj_nuc_c, adj_nuc); // adj_nuc -> adj_nuc_c
          real snstot=snsval[0]*(1.-wc_gpt_wpc)+snsval[1]*wc_gpt_wpc; 
	  sns_tmp1[j]=snstot;

	};

	// (Make constrained sensitivity)
  	vector<real> sns_tmp2;
  	sns_tmp2.resize(channel);
  	for(int j=0;j<channel;j++){
	  real tmps=0.;
	  for(int jj=0;jj<channel;jj++){
	    tmps+=sns_tmp1[jj];
	  };
	  sns_tmp2[j]=sns_tmp1[j]-tmps*ratio[j];
	};

	for(int j=0;j<channel;j++){
	  int mt=88880+j;
	  sns.PutSensitivity0D(matnum,mt,sns_tmp2[j]);
	};
      };
    };
  };

  //

  sns.WriteFile("./","sns.k_EOC_indir");
  sns.AddSensitivityData(sns_dir);
  sns.WriteFile("./","sns.k_EOC");
  
};

