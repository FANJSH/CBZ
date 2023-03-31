#include <cstdlib>
#include "Pij_dyn.h"
ofstream fs1("phi_f.txt");
ofstream fs2("phi_t.txt");
ofstream fs3("DNP.txt");
ofstream fs4("DNPcal.txt");
ofstream fs5("an.txt");
ofstream fs7("look.txt");
ofstream fs9("output.txt");



using namespace std;

void Pij_dyn::Initialize(int ttotM, int allmat,int ene,int dp,int ffreg,bool doku,int brr)
{
    chainef = doku;
    hosei = 1;
    //br_seed=brr;
    br_seed = 2;
    use_mmpa = true;
    if (dp == 6) {
        chainef = false;
    }
    //とりあえずここでBCGManegerに取り込むっぽい？
    // +++ Decay Data input
    man.ReadDecayDataFromFile("../../CBGLIB/DDfile/dd.jendl2011");
    //man.ReadDecayDataFromFile("../../CBGLIB/DDfile/dd.endf71");
    // +++ brancing ratio input
    man.ReadNGBranchingRatioDataFromFile("../../CBGLIB/BR_DATA/ng_br_jndc_v2");
    // JNDC Ver.2 for thermal reactor
    // +++ Fission yield data input
    //man.ReadFPYieldDataFromFile("../../CBGLIB/FPYfile_cum/JENDL-FPY-2011/U235",0,"u235");
    //以上、取り込み完了*/
    //man.ShowNuclideList();
    if (ene==2) {
    fueld.ReadFile("../../CBGLIB/kmcal/","pinc2g.0");//2群
     cladd.ReadFile("../../CBGLIB/kmcal/","pinc2g.1");
     modd.ReadFile("../../CBGLIB/kmcal/","pinc2g.2");//*/
    }
    if (ene==4) {
    fueld.ReadFile("../../CBGLIB/kmcal/","pinc0_4");//4群
    cladd.ReadFile("../../CBGLIB/kmcal/","pinc1_4");
    modd.ReadFile("../../CBGLIB/kmcal/","pinc2_4");//*/
    }
    if (ene==107) {
    fueld.ReadFile("../../CBGLIB/kmcal/","pin0");//107群用
    cladd.ReadFile("../../CBGLIB/kmcal/","pin1");
    modd.ReadFile("../../CBGLIB/kmcal/","pin2");//*/
    }
    if (fueld.GetImax()!=ene||cladd.GetImax()!=ene||modd.GetImax()!=ene) {
        cout<<"data is not correct"<<endl;
    }
    int sz;
    if (chainef == true) {
        man.ReadFPYieldDataFromFile("../../CBGLIB/FPYfile/JENDL-FPY-2011/U235",0,"u235");
        man.ReadFPYieldDataFromFile("../../CBGLIB/FPYfile/JENDL-FPY-2011/U238",1,"u238");
        //man.ShowYieldDataForXYPlot("u235");
        sz = man.GetSize();
        fissN.resize(sz);
        fpid.resize(sz);
        fpid_nspect.resize(sz);
        nspect.resize(sz);
        for (int i=0; i<sz; i++) {
            nspect[i].resize(ene);
        }
        man.CalTotalEmittedNeutrons();
        //man.CalTotalNeutrons();
        man.CalTEN();
        man.CalRRnums();
        for (int i=0; i<sz; i++) {
            int ngc=man.GetNuclide(i).GetChannel();
            //cout<<"No."<<i<<" TEN:"<<man.GetNuclide(i).GetTEN()<<endl;
            if (ngc != 0) {
            for (int j=0; j<ngc;j++) {
                int nonb=man.GetNuclide(i).GetRRnum(j);
                //cout<<"No."<<i<<" to No."<<nonb<<" by reactNo."<<j<<endl;

            }}
        }
        cout<<"now ceswitch on"<<endl;
        if (ene==107) {
        ifstream fin("../../CBGLIB/kmcal/endffp_e71.107g");//endffp_e71.107g
        int iz,mmm;
        mmm=0;
        fin>>mmm;
        //cout<<"yomikomitest : "<<mmm<<endl;
        iz=0;
        int count=0;
        while(iz!=-1){
            fin>>iz;
            if(iz==-1)break;
            fpid_nspect[count]=iz;
            //cout<<iz<<endl;
            for (int i=0; i<107; i++) {//ここ、柔軟性無い書き方してるんで適当にお願いします。
                real buf;
                fin>>buf;
                nspect[count][i]=buf;
            }
            count++;
        };
        for (int i=0; i<count; i++) {
            real buf = 0;
            for (int j=0; j<ene; j++) {
                buf+=nspect[i][j];
            }
            buf=1/buf;
            for (int j=0; j<ene; j++) {
                nspect[i][j]*=buf;
            }
        }
        cout<<"count total : "<<count<<endl;
        fp_spect_non=count;
        }
    }else{
        man.ReadFPYieldDataFromFile("../../CBGLIB/FPYfile_cum/JENDL-FPY-2011/U235",0,"u235");
        man.ReadFPYieldDataFromFile("../../CBGLIB/FPYfile_cum/JENDL-FPY-2011/U238",1,"u238");
        sz = man.GetSize();
        fissN.resize(sz);
        fpid.resize(sz);
        fpid_nspect.resize(sz);
        nspect.resize(sz);
        for (int i=0; i<sz; i++) {
            nspect[i].resize(ene);
        }
        man.CalTotalNeutrons();
        cout<<"now ceswitch off"<<endl;
        if (ene==107) {
        ifstream fin("../../CBGLIB/kmcal/endffp_e71.107g");//endffp_e71.107g

        int iz,mmm;
        mmm=0;
        fin>>mmm;
        cout<<mmm<<endl;
        iz=0;
        int count=0;
        while(iz!=-1){
            fin>>iz;
            if(iz==-1)break;
            fpid_nspect[count]=iz;
            for (int i=0; i<107; i++) {//ここ、柔軟性無い書き方してるんで適当にお願いします。
                real buf;
                fin>>buf;
                nspect[count][i]=buf;
            }
            count++;
        };
        for (int i=0; i<count; i++) {
            real buf = 0;
            for (int j=0; j<ene; j++) {
                buf+=nspect[i][j];
            }
            buf=1/buf;
            for (int j=0; j<ene; j++) {
                nspect[i][j]*=buf;
            }
        }
        fp_spect_non=count;
        }
    }
    cout<<sz<<endl;
    toriaezu = 0;
    for (int i=0; i<sz; i++) {
        //int mass = man.GetNuclide(i).GetMassNumber();
        //int atm  = man.GetNuclide(i).GetAtomicNumber();
        //int Exl  = man.GetNuclide(i).GetExLevel();
        //real hl  = man.GetNuclide(i).GetHalflife();
        real Yel = man.GetNuclide(i).GetYield("u235");
        toriaezu += man.GetNuclide(i).GetTotalNumberEmittedNeutrons()*Yel;
        //cout<<hl<<endl;
        //cout<<"atm: "<<atm<<" mass: "<<mass<<" Exl: "<<Exl<<" Yield: "<<Yel<<" hl: "<<hl<<" TNE: "<<TNE<<endl;
    }//*/
    cout<<"yield sum : "<<man.GetYieldSum("u235")<<endl;;//
    cout<<"total Nutrons : "<<toriaezu<<endl;
    freg = ffreg;
    mednum = allmat;
    group = ene;
    totM = ttotM;
    time = 0;
    if (dp == 6) {
        totDNP = dp;//6群用？
    }else{
        totDNP = sz;//全核種網羅用
    }
    keff = 1;
    DNPlam.resize(totDNP);//totDNP*mednumなら複数種の核分裂にも対応可
    DNPbeta.resize(totDNP);//ただ、現状メモリの無駄として切りました
    u8y.resize(totDNP);
    //size=totM*(group+totDNP);
    size = totM*group+freg*totDNP;
    sizes=mednum*group*(2+group);
    cout<<"vector size : "<<size<<endl;
    cal_rad.resize(totM+1);
    vec_a.resize(size);
    vec_b.resize(size);
    mat_D.resize(mednum*group);//energy+gruop*matidで？
    mmpa_a.resize(size);
    mmpa_an.resize(size);
    mmpa_b.resize(size);
    mmpa_bg.resize(size);
    mmpa_mat.resize(size);
    mmpa_ans.resize(size);
    vec_bcopy.resize(size);
    vec_bgaus.resize(size);
    phi.resize(size);
    phi_n.resize(size);
    Pij_c.resize(totM);
    cs_c.resize(totM);
    med.resize(totM);
    chiC.resize(group*mednum);
    chiD.resize(group*totDNP);
    Pij.resize(totM*totM*group);
    Pij_n.resize(totM*totM*group);
    vol.resize(totM);
    cs.resize(sizes);
    cs_n.resize(sizes);
    vel.resize(group);
    temperature.resize(totM);
    heattrans.resize(totM);
    heatoutput.resize(totM);
    specific_heat.resize(mednum);
    rho.resize(mednum);
    radius.resize(totM);
    pijregid.resize(totM);
    fission_S.resize(totM);
    scat_S.resize(totM*group);
    br_cal.resize(group);
    for (int i=0; i<size; i++) {
        vec_b[i].resize(size);
        mmpa_a[i].resize(size);
        mmpa_b[i].resize(size);
        mmpa_bg[i].resize(size);
        mmpa_mat[i].resize(size);
        vec_bgaus[i].resize(size);
        vec_bcopy[i].resize(size);
        phi[i]=1.;//初期値、これも本来は関数で決めたいとか
        mmpa_ans[i]=0;
        for (int j=0; j<size; j++) {
            mmpa_a[i][j]=0;
            mmpa_b[i][j]=0;
            mmpa_bg[i][j]=0;
            mmpa_mat[i][j]=0;
        }
    }
//pijcal用
    real radius[]={0.714,0.664,0.611,0.553,0.488,0.412,0.291};//将来的にはこれも外部から入れたいよねー
    //real radius[]={0.714,0.664,0.611,0.553,0.488,0.460,0.412};//将来的にはこれも外部から入れたいよねー
    //real radius[]={0.714,0.412};//将来的にはこれも外部から入れたいよねー
    int regid[]={6,5,4,3,2,1,0};//無論これも
    //int regid[]={1,0};//無論これも
    igi.AddCircleRing(totM,radius,regid);

    vol[0]=radius[6]*radius[6]*3.14;
    for (int i=1; i<totM; i++) {
        vol[i]=(radius[6-i]*radius[6-i]-radius[7-i]*radius[7-i])*3.14;
    }
    /*
    vol[0]=radius[1]*radius[1]*3.14;
    vol[1]=radius[0]*radius[0]*3.14-vol[0];
    */
    igi.WriteGnuplotFile(0.01);
    sys.PutBoundaryCondition(White);
    sys.CalTrajectory(igi,1,0.001,45.); // anguler division/width/angle

};

void Pij_dyn::cal_Pij(int energyy)
{
    real *xs=new real[totM];
    for (int i=0; i<totM; i++) {
        xs[i]=cs[0*mednum*group+med[i]*group+energyy];
    }
     real *pij=new real[totM*totM];
     sys.CalculationPij(xs,pij);
    //Pij[i*totM+j+totM*totM*ene] iからjへの移行率。エネルギーはene
    for (int i=0; i<totM; i++) {
        for (int j=0; j<totM; j++) {
            Pij[i*totM+j+totM*totM*energyy]=pij[i*totM+j];
        }
    }
    /*for(int i=0;i<totM;i++){
     for(int j=0;j<totM;j++){
     cout<<i<<" "<<j<<" "<<pij[i*totM+j]<<"\n";
     };
     };//*/
    delete[] xs;
    delete[] pij;//
    //cout<<"-----------"<<endl;
}
void Pij_dyn::cal_xs_f(){
    enrich_f=3.4;//3.4%濃縮率
    real f_density = 10.0952388812265;//とりあえず燃料ペレットの密度を確定
    real u8_density = 0;
    real u5_density = 0;
    real o6_density = 0;//各密度をここで
    real u5_per=238*(enrich_f/100)/(235+238*(enrich_f/100)-235*(enrich_f/100));
    real mol = f_density/(u5_per*235+(1-u5_per)*238+2*16);
    real allperticle = mol*6.0221413E+23*1.0E-24;
    u5_density = allperticle*u5_per;
    u8_density = allperticle*(1-u5_per);
    o6_density = allperticle*2;
    //cout<<u5_density<<" , "<<u8_density<<" , "<<o6_density<<endl;
    //ここからした、２群である。
    /*real u5_xs[]={3.1369e+01,3.0235e+02,2.4375e+00*1.3227e+01,2.4363e+00*2.4455e+02,3.1369e+01-7.7263e+00- 6.9448e-03-1.32E+01,6.9448e-03,3.0235e+02-4.3648e+01-2.45E+02};///total1-2,nyu_sig_f1-2,scat1-1,1-2,2-2
    real u8_xs[]={1.4059e+01,1.0516e+01,2.8170e+00*4.7589e-02,2.3277e+00*7.9817e-06,1.4059e+01-1.9380e+00-5.1334e-03-4.76E-02,5.1334e-03,1.0516e+01-1.3002e+00-7.98E-06};
    real o6_xs[]={3.6869e+00,3.9971e+00,0.0000e+00,0.0000e+00,3.6869e+00-1.4090e-03-2.9904e-02,2.9904e-02,3.9971e+00-8.8648e-05};
    br_cal[0]=1.3227e+01*u5_density/(1.3227e+01*u5_density+4.7589e-02*u8_density);//ここも自動化したい。したい。
    br_cal[1]=2.4455e+02*u5_density/(2.4455e+02*u5_density+7.9817e-06*u8_density);
    cs[0*mednum*group+0*group+0]=u5_density*u5_xs[0]+u8_density*u8_xs[0]+o6_density*o6_xs[0];//fuel ver.micro
    cs[0*mednum*group+0*group+1]=u5_density*u5_xs[1]+u8_density*u8_xs[1]+o6_density*o6_xs[1];
    cs[1*mednum*group+0*group+0]=u5_density*u5_xs[2]+u8_density*u8_xs[2]+o6_density*o6_xs[2];
    cs[1*mednum*group+0*group+1]=u5_density*u5_xs[3]+u8_density*u8_xs[3]+o6_density*o6_xs[3];
    cs[2*mednum*group+0*group*group+0+group*0]=u5_density*u5_xs[4]+u8_density*u8_xs[4]+o6_density*o6_xs[4];
    cs[2*mednum*group+0*group*group+0+group*1]=u5_density*u5_xs[5]+u8_density*u8_xs[5]+o6_density*o6_xs[5];
    cs[2*mednum*group+0*group*group+1+group*1]=u5_density*u5_xs[6]+u8_density*u8_xs[6]+o6_density*o6_xs[6];
    //cout<<cs[0*mednum*group+0*group+0]<<" , "<<cs[0*mednum*group+0*group+1]<<" , "<<cs[1*mednum*group+0*group+0]<<" , "<<cs[1*mednum*group+0*group+1]<<" , "<<cs[2*mednum*group+0*group*group+0+group*0]<<" , "<<cs[2*mednum*group+0*group*group+0+group*1]<<" , "<<cs[2*mednum*group+0*group*group+1+group*1]<<endl;//*/
    //更にここから下が４群です。求む、他クラスからの引用。
    //ccc
    /*real u5_tot[]={9.91352e+00,4.49299e+01,1.05350e+02,5.06947e+02};//全断面積
    real u8_tot[]={9.89441e+00,1.66917e+01,9.70495e+00,1.13593e+01};
    real o6_tot[]={3.44917e+00,3.83721e+00,3.85941e+00,4.14020e+00};
    real u5_nf[]={2.52268e+00*1.43025e+00,2.43381e+00*2.06830e+01,2.43633e+00*7.61167e+01,2.43633e+00*4.19493e+02};//ν*Σf
    real u8_nf[]={2.81824e+00*1.22557e-01,2.32796e+00*2.04440e-04,2.32774e+00*3.44759e-06,2.32774e+00*1.26909e-05};
    real u5_fis[]={1.43025e+00,2.06830e+01,7.61167e+01,4.19493e+02};//Σf
    real u8_fis[]={1.22557e-01,2.04440e-04,3.44759e-06,1.26909e-05};
    real u5_scat[]={8.16375e+00,2.76414e-02,6.48382e-11,0.00000e+00,0.00000e+00,1.18093e+01,1.13344e-02,0.00000e+00,0.00000e+00,0.00000e+00,1.33212e+01,4.63637e-02,0.00000e+00,0.00000e+00,0.00000e+00,1.49027e+01};
    real u8_scat[]={9.58040e+00,2.73156e-02,1.45724e-10,0.00000e+00,0.00000e+00,1.35700e+01,8.37809e-03,0.00000e+00,0.00000e+00,0.00000e+00,9.07833e+00,2.95480e-02,0.00000e+00,0.00000e+00,0.00000e+00,9.30310e+00};
    real o6_scat[]={3.36667e+00,8.36729e-02,0.00000e+00,0.00000e+00,0.00000e+00,3.78792e+00,4.88056e-02,1.33335e-08,0.00000e+00,5.00883e-03,3.66905e+00,2.02802e-01,0.00000e+00,0.00000e+00,1.74197e-01,3.99182e+00};//散乱*/
    //922350 : u-235
    //922380 : u-238
    //80160  : o-16
    for (int i=0; i<group; i++) {
        u5_density=fueld.GetNuclide(922350).GetDensity();
        u8_density=fueld.GetNuclide(922380).GetDensity();
        o6_density=fueld.GetNuclide(80160).GetDensity();//*/
        br_cal[i]=fueld.GetNuclide(922350).GetMicxs().GetData1d(sigf).get_dat(i)*u5_density/(fueld.GetNuclide(922350).GetMicxs().GetData1d(sigf).get_dat(i)*u5_density+fueld.GetNuclide(922380).GetMicxs().GetData1d(sigf).get_dat(i)*u8_density);
        cs[0*mednum*group+0*group+i]=u5_density*fueld.GetNuclide(922350).GetMicxs().GetData1d(sigt).get_dat(i)+u8_density*fueld.GetNuclide(922380).GetMicxs().GetData1d(sigt).get_dat(i)+o6_density*fueld.GetNuclide(80160).GetMicxs().GetData1d(sigt).get_dat(i);
        cs[1*mednum*group+0*group+i]=u5_density*fueld.GetNuclide(922350).GetMicxs().GetData1d(sigf).get_dat(i)*fueld.GetNuclide(922350).GetMicxs().GetData1d(nu).get_dat(i)+u8_density*fueld.GetNuclide(922380).GetMicxs().GetData1d(sigf).get_dat(i)*fueld.GetNuclide(922380).GetMicxs().GetData1d(nu).get_dat(i);
        mat_D[i+group*0]=1/(3*cs[0*mednum*group+0*group+i]);
        for (int j=0; j<group; j++) {
            cs[2*mednum*group+0*group*group+i+group*j]=u5_density*(fueld.GetNuclide(922350).GetMicxs().GetData2d(sigel).get_dat(i,j)+fueld.GetNuclide(922350).GetMicxs().GetData2d(siginel).get_dat(i,j))+u8_density*(fueld.GetNuclide(922380).GetMicxs().GetData2d(sigel).get_dat(i,j)+fueld.GetNuclide(922380).GetMicxs().GetData2d(siginel).get_dat(i,j))+o6_density*(fueld.GetNuclide(80160).GetMicxs().GetData2d(sigel).get_dat(i,j)+fueld.GetNuclide(80160).GetMicxs().GetData2d(siginel).get_dat(i,j));
        }//*/
    }
    /*for (int i=0; i<4; i++) {
        br_cal[i]=u5_fis[i]*u5_density/(u5_fis[i]*u5_density+u8_fis[i]*u8_density);
        cs[0*mednum*group+0*group+i]=u5_density*u5_tot[i]+u8_density*u8_tot[i]+o6_density*o6_tot[i];//
        cs[1*mednum*group+0*group+i]=u5_density*u5_nf[i]+u8_density*u8_nf[i];
        cs[2*mednum*group+0*group*group+i+group*0]=u5_density*u5_scat[i*4+0]+u8_density*u8_scat[i*4+0]+o6_density*o6_scat[i*4+0];
        cs[2*mednum*group+0*group*group+i+group*1]=u5_density*u5_scat[i*4+1]+u8_density*u8_scat[i*4+1]+o6_density*o6_scat[i*4+1];
        cs[2*mednum*group+0*group*group+i+group*2]=u5_density*u5_scat[i*4+2]+u8_density*u8_scat[i*4+2]+o6_density*o6_scat[i*4+2];
        cs[2*mednum*group+0*group*group+i+group*3]=u5_density*u5_scat[i*4+3]+u8_density*u8_scat[i*4+3]+o6_density*o6_scat[i*4+3];

    }
    cout<<cs[0*mednum*group+0*group+0]<<" , "<<cs[0*mednum*group+0*group+1]<<" , "<<cs[1*mednum*group+0*group+0]<<" , "<<cs[1*mednum*group+0*group+1]<<" , "<<cs[2*mednum*group+0*group*group+0+group*0]<<" , "<<cs[2*mednum*group+0*group*group+0+group*1]<<" , "<<cs[2*mednum*group+0*group*group+1+group*1]<<endl;//*/

}
void Pij_dyn::fact_sigt_mod(real factor)
{
    for (int i=0; i<group; i++) {
        cs[0*mednum*group+2*group+i]=factor*modd.GetMacxs().GetData1d(sigt).get_dat(i);
    }
}
void Pij_dyn::fact_sig_mod(real factor)
{
    for (int i=0; i<group; i++) {
        cs[0*mednum*group+2*group+i]=factor*modd.GetMacxs().GetData1d(sigt).get_dat(i);
        cs[1*mednum*group+2*group+i]=factor*modd.GetMacxs().GetData1d(nusigf).get_dat(i);
        for (int j=0; j<group; j++) {
            cs[2*mednum*group+2*group*group+i+group*j]=factor*modd.GetMacxs().GetData2d(sigs).get_dat(i,j);
        }
    }//*/
}
void Pij_dyn::Adjust_bq(real factor)
{
    for (int k=0;k<totM ; k++) {
    for (int i=0; i<group; i++) {
        cs[0*mednum*group+k*group+i]+=mat_D[i+group*k]*factor;
    }//*/
    }
}


void Pij_dyn::UseData()
{
    theta = 1;
    Totaltime = 10;
    timestep = 1.0E-5;
    step = Totaltime/timestep;
    //step = 1000000;//メモ：即発臨界時、中性子寿命よりタイムステップを長くすると計算が収束しない可能性
    //timestep = Totaltime/step;
    //2gun
    if (group==2) {
    vel[0] = 1.0E+07;
    vel[1] = 2.0E+05;//*/
    }
    //ccc
    if (group==4) {
    vel[0] = 1.0E+07;
    vel[1] = 1.0E+06;
    vel[2] = 5.0E+05;
    vel[3] = 2.0E+05;//*/
    }
    if (group==107) {//nh
        real eneene[]={1.0000E+07,7.7880E+06,6.0653E+06,4.7237e+06,3.6788e+06,2.8650e+06,2.2313e+06,1.7377e+06,1.3534e+06,1.0540e+06,8.2085e+05,6.3928e+05,4.9787e+05,3.8774e+05,3.0197e+05,2.3518e+05,1.8316e+05,1.4264e+05,1.1109e+05,8.6517e+04,6.7379e+04,5.2475e+04,4.0868e+04,3.1828e+04,2.4788e+04,1.9305e+04,1.5034e+04,1.1709e+04,9.1188e+03,7.1017e+03,5.5308e+03,4.3074e+03,3.3546e+03,2.6126e+03,2.0347e+03,1.5846e+03,1.2341e+03,9.6112e+02,7.4852e+02,5.8295e+02,4.5400e+02,3.5358e+02,2.7536e+02,2.1445e+02,1.6702e+02,1.3007e+02,1.0130e+02,7.8893e+01,6.1442e+01,4.7851e+01,3.7267e+01,2.9023e+01,2.2603e+01,1.7603e+01,1.3710e+01,1.0677e+01,8.3153e+00,6.4760e+00,5.0435e+00,3.9279e+00,3.0590e+00,2.3824e+00,1.8554e+00,1.6374e+00,1.4450e+00,1.2752e+00,1.1254e+00,9.9312e-01,8.7642e-01,7.7344e-01,6.8256e-01,6.0236e-01,5.3158e-01,4.6912e-01,4.1399e-01,3.8926e-01,3.6528e-01,3.4205e-01,3.1959e-01,2.9790e-01,2.7689e-01,2.5681e-01,2.3740e-01,2.1875e-01,2.0087e-01,1.8375e-01,1.6739e-01,1.5180e-01,1.3697e-01,1.2290e-01,1.0960e-01,9.7052e-02,8.5375e-02,7.4258e-02,6.4004e-02,5.4508e-02,4.5776e-02,3.7805e-02,3.0595e-02,2.4149e-02,1.8463e-02,1.1354e-02,9.3788e-03,5.9796e-03,3.3419e-03,1.4662e-03,3.5234e-04};//eV
        real nmass = 1.67e-27;//kg
        real evtoj = 1.6e-19;//1eV=1.6e-19J
        for (int i=0; i<group; i++) {
            vel[i]=sqrt(eneene[i]*evtoj*2/nmass);
            //cout<<"vel["<<i<<"] = "<<vel[i]<<endl;
        }
    }//*/
    firstoutput = 160;//W/cm
    hanteiki=1.0e-12;
    nowtime=0;
    //keff=1;
    hantei=1;
    nyt=1.0;
    //改修後は、DNPbetaは独立収率、DNPlamは崩壊定数となります。
    //それと、とりあえず崩壊後の中性子は高速群に投げ入れてみるテスト
    if (totDNP == 6) {
        //beta = 0;
        beta = 0.0065;
        nyd = beta;
        hosei = 1.0;
        nyt = 1.0;
        DNPbeta[0] = beta*0.033;
        DNPbeta[1] = beta*0.219;
        DNPbeta[2] = beta*0.196;
        DNPbeta[3] = beta*0.395;
        DNPbeta[4] = beta*0.115;
        DNPbeta[5] = beta*0.042;
        DNPlam[0] = 0.01244;
        DNPlam[1] = 0.03054;
        DNPlam[2] = 0.1114;
        DNPlam[3] = 0.3014;
        DNPlam[4] = 1.136;
        DNPlam[5] = 3.014;//*/
        /*DNPlam[0] = 0.001244;//ダメな奴。
         DNPlam[1] = 0.003054;
         DNPlam[2] = 0.01114;
         DNPlam[3] = 0.03014;
         DNPlam[4] = 0.1136;
         DNPlam[5] = 0.3014;//*/

        for (int i=0; i<totDNP; i++) {
            fissN[i]=1.0;
        }
        for (int i=0; i<group; i++) {
            for (int j=0; j<totDNP; j++) {
                chiD[i+group*j]=fueld.GetMacxs().GetData1d(chi).get_dat(i);//これが基礎。
            }
        }
    }else{
        beta = 0.0065;
        hosei = beta/toriaezu;
        //beta = toriaezu;
        //beta=0;
        nyd = 0;
        for (int i = 0; i<totDNP; i++) {
            //ここで諸入力（引っ張ってきて。
            real a     = man.GetNuclide(i).GetYield("u235");//cum;
            real b     = 0;
            if (chainef == true) {
                b     = man.GetNuclide(i).GetTEN();
            }else{
                b     = man.GetNuclide(i).GetTotalNumberEmittedNeutrons();
            }
            DNPbeta[i] = a;
            u8y[i] = man.GetNuclide(i).GetYield("u238");
            DNPlam[i]  = log(2)/man.GetNuclide(i).GetHalflife();
            fissN[i]   = b;
            nyd +=a*b;
            fpid[i] = man.GetNuclide(i).GetID();
            //cout<<fpid[i]<<endl;
            if (man.GetNuclide(i).GetHalflife() == 0) {
                DNPlam[i] = 0;
            }
        }
        nyt = nyd/beta;
        cout<<"nd = "<<nyd<<endl;
        cout<<"nt = "<<nyt<<endl;
        for (int i=0; i<group; i++) {
            for (int j=0; j<totDNP; j++) {
                chiD[i+group*j]=fueld.GetMacxs().GetData1d(chi).get_dat(i);//これが基礎。
            }
        }
        if (group==107) {//nh
        for (int i=0;i<totDNP; i++) {
            int pid = fpid[i];
            for (int j=0; j<fp_spect_non; j++) {
                int bbuf=fpid_nspect[j];
                if (pid==bbuf) {
                    for (int k=0; k<group; k++) {
                        chiD[k+group*i]=nspect[j][k];
                    }
                    //cout<<"No : "<<bbuf<<" OK"<<endl;
                }
            }
        }
        }//*/
    }
    //2gun
    /*for (int i = 0; i<totDNP*2; i++) {
        if (i%2==0) {
            chiD[i]=1;
        }else{
            chiD[i]=0;
        }
    }//*/
    //4群
    //ccc
    /*for (int i = 0; i<totDNP*group; i++) {
     if (i%group==0) {
     chiD[i]=1;
     }else{
     chiD[i]=0;
     }
     }//*/

    /*chiD[0] = 1;
    chiD[1] = 0;
    chiD[2] = 1;
    chiD[3] = 0;
    chiD[4] = 1;
    chiD[5] = 0;
    chiD[6] = 1;
    chiD[7] = 0;
    chiD[8] = 1;
    chiD[9] = 0;
    chiD[10] = 1;
    chiD[11] = 0;//ここは、六群用ですので*/
    for (int i=0; i<totM; i++) {
        med[i] = 2;
    }
    for (int i=0; i<freg; i++) {
        med[i] = 0;
        heattrans[i] = 0.033472;//J/s/K/cm
    }//全体均質の燃料として扱うはずもなく、中央部ノミにしてみよう！
    //2群
    //ccc
    /*for (int i=0; i<mednum*group; i++) {
        if (i%group==0) {
            chiC[i]=1;
        }else{
            chiC[i]=0;
        }
    }//*/
    /*chiC[0] = 1;
    chiC[1] = 0;
    chiC[2] = 1;
    chiC[3] = 0;
    chiC[4] = 1;
    chiC[5] = 0 ;//mednum=3での*/
    for (int i=0; i<group; i++) {
        chiC[i]=fueld.GetMacxs().GetData1d(chi).get_dat(i);
        chiC[i+group*1]=cladd.GetMacxs().GetData1d(chi).get_dat(i);
        chiC[i+group*2]=modd.GetMacxs().GetData1d(chi).get_dat(i);
    }
//*/
    //cout<<"cp1"<<endl;
    
    //100度っぽい。
    /*cs[0*mednum*group+0*group+0]=4.95E-01;//fuel ver.yoko
    cs[0*mednum*group+0*group+1]=6.44E-01;
    cs[1*mednum*group+0*group+0]=2.79E-02;
    cs[1*mednum*group+0*group+1]=4.62E-01;
    cs[2*mednum*group+0*group*group+0+group*0]=4.34E-01;
    cs[2*mednum*group+0*group*group+0+group*1]=1.47E-03;
    cs[2*mednum*group+0*group*group+1+group*1]=3.92E-01;*/
    
    /*cs[0*mednum*group+0*group+0]=4.9550e-01;//fuel ver.micro
    cs[0*mednum*group+0*group+1]=6.4390e-01;
    cs[1*mednum*group+0*group+0]=2.7912e-02;
    cs[1*mednum*group+0*group+1]=4.6193e-01;
    cs[2*mednum*group+0*group*group+0+group*0]=4.345418E-01;
    cs[2*mednum*group+0*group*group+0+group*1]=1.4642e-03;
    cs[2*mednum*group+0*group*group+1+group*1]=3.92170E-01;//*/

    //ここのデータ類も２群ベースです。４群用の作成を。
    /*cs[0*mednum*group+0*group+0]=0.49619992027778359045;//fuel ver.micro -mk2
    cs[0*mednum*group+0*group+1]=0.64320587525823391672;
    cs[1*mednum*group+0*group+0]=0.027912044072392287047;
    cs[1*mednum*group+0*group+1]=0.46192158937837340948;
    cs[2*mednum*group+0*group*group+0+group*0]=0.43526133906053399159;
    cs[2*mednum*group+0*group*group+0+group*1]=0.0014642265852958034726;
    cs[2*mednum*group+0*group*group+1+group*1]=0.39113378194423586987;//
    
    cs[0*mednum*group+1*group+0]=2.82E-01;//cladding
    cs[0*mednum*group+1*group+1]=2.49E-01;
    cs[1*mednum*group+1*group+0]=0;
    cs[1*mednum*group+1*group+1]=0;
    cs[2*mednum*group+1*group*group+0+group*0]=2.80E-01;
    cs[2*mednum*group+1*group*group+0+group*1]=3.78E-04;
    cs[2*mednum*group+1*group*group+1+group*1]=2.45E-01;//
    cs[0*mednum*group+2*group+0]=9.91E-01;//mod
    cs[0*mednum*group+2*group+1]=2.067E+00;
    cs[1*mednum*group+2*group+0]=0;
    cs[1*mednum*group+2*group+1]=0;
    cs[2*mednum*group+2*group*group+0+group*0]=9.10E-01;
    cs[2*mednum*group+2*group*group+0+group*1]=8.01E-02;
    cs[2*mednum*group+2*group*group+1+group*1]=2.01E+00;//*/
    for (int i=0; i<group; i++) {
        cs[0*mednum*group+0*group+i]=fueld.GetMacxs().GetData1d(sigt).get_dat(i);
        cs[0*mednum*group+1*group+i]=cladd.GetMacxs().GetData1d(sigt).get_dat(i);
        cs[0*mednum*group+2*group+i]=modd.GetMacxs().GetData1d(sigt).get_dat(i);
        cs[1*mednum*group+0*group+i]=fueld.GetMacxs().GetData1d(nusigf).get_dat(i);
        cs[1*mednum*group+1*group+i]=cladd.GetMacxs().GetData1d(nusigf).get_dat(i);
        cs[1*mednum*group+2*group+i]=modd.GetMacxs().GetData1d(nusigf).get_dat(i);
        mat_D[i+group*0]=fueld.GetMacxs().GetData1d(d).get_dat(i);
        mat_D[i+group*1]=cladd.GetMacxs().GetData1d(d).get_dat(i);
        mat_D[i+group*2]=modd.GetMacxs().GetData1d(d).get_dat(i);
        for (int j=0; j<group; j++) {
            cs[2*mednum*group+0*group*group+i+group*j]=fueld.GetMacxs().GetData2d(sigs).get_dat(i,j);
            cs[2*mednum*group+1*group*group+i+group*j]=cladd.GetMacxs().GetData2d(sigs).get_dat(i,j);
            cs[2*mednum*group+2*group*group+i+group*j]=modd.GetMacxs().GetData2d(sigs).get_dat(i,j);
        }
    }//*/
    cout<<"modd:2 , "<<cs[0*mednum*group+2*group+1]<<endl;
    /*cout<<cs[0*mednum*group+0*group+0]<<endl;
    cout<<cs[0*mednum*group+0*group+1]<<endl;
    cout<<cs[1*mednum*group+0*group+0]<<endl;
    cout<<cs[1*mednum*group+0*group+1]<<endl;
    cout<<cs[2*mednum*group+0*group*group+0+group*0]<<endl;
    cout<<cs[2*mednum*group+0*group*group+0+group*1]<<endl;
    cout<<cs[2*mednum*group+0*group*group+1+group*1]<<endl;//
    
    cout<<cs[0*mednum*group+1*group+0]<<endl;//cladding
    cout<<cs[0*mednum*group+1*group+1]<<endl;
    cout<<cs[1*mednum*group+1*group+0]<<endl;
    cout<<cs[1*mednum*group+1*group+1]<<endl;
    cout<<cs[2*mednum*group+1*group*group+0+group*0]<<endl;
    cout<<cs[2*mednum*group+1*group*group+0+group*1]<<endl;
    cout<<cs[2*mednum*group+1*group*group+1+group*1]<<endl;//
    cout<<cs[0*mednum*group+2*group+0]<<endl;//mod
    cout<<cs[0*mednum*group+2*group+1]<<endl;
    cout<<cs[1*mednum*group+2*group+0]<<endl;
    cout<<cs[1*mednum*group+2*group+1]<<endl;
    cout<<cs[2*mednum*group+2*group*group+0+group*0]<<endl;
    cout<<cs[2*mednum*group+2*group*group+0+group*1]<<endl;
    cout<<cs[2*mednum*group+2*group*group+1+group*1]<<endl;*/
    //matID:: fuel:0 clad:1 mod:2
    //以下４群用。
    //ccc
     /*cs[0*mednum*group+0*group+0]=3.78450e-01;//fuel
     cs[0*mednum*group+0*group+1]=5.69480e-01;
     cs[0*mednum*group+0*group+2]=4.67405e-01;
     cs[0*mednum*group+0*group+3]=8.27200e-01;
     cs[1*mednum*group+0*group+0]=1.03097e-02;
     cs[1*mednum*group+0*group+1]=3.90378e-02;
     cs[1*mednum*group+0*group+2]=1.43776e-01;
     cs[1*mednum*group+0*group+3]=7.92375e-01;
     cs[2*mednum*group+0*group*group+0+group*0]=3.66372e-01;
     cs[2*mednum*group+0*group*group+0+group*1]=4.38501e-03;
     cs[2*mednum*group+0*group*group+0+group*2]=3.21977e-12;
     cs[2*mednum*group+0*group*group+0+group*3]=0.00000e+00;
    
     cs[2*mednum*group+0*group*group+1+group*0]=0.00000e+00;
     cs[2*mednum*group+0*group*group+1+group*1]=4.74949e-01;
     cs[2*mednum*group+0*group*group+1+group*2]=2.38970e-03;
     cs[2*mednum*group+0*group*group+1+group*3]=6.00676e-10;

     cs[2*mednum*group+0*group*group+2+group*0]=0.00000e+00;
     cs[2*mednum*group+0*group*group+2+group*1]=2.25648e-04;
     cs[2*mednum*group+0*group*group+2+group*2]=3.73072e-01;
     cs[2*mednum*group+0*group*group+2+group*3]=9.81483e-03;
    
     cs[2*mednum*group+0*group*group+3+group*0]=0.00000e+00;
     cs[2*mednum*group+0*group*group+3+group*1]=0.00000e+00;
     cs[2*mednum*group+0*group*group+3+group*2]=7.84759e-03;
     cs[2*mednum*group+0*group*group+3+group*3]=3.93728e-01;;//
    
    cs[0*mednum*group+1*group+0]=2.90507e-01;//cladding
    cs[0*mednum*group+1*group+1]=2.77505e-01;
    cs[0*mednum*group+1*group+2]=2.46184e-01;
    cs[0*mednum*group+1*group+3]=2.52170e-01;
    cs[1*mednum*group+1*group+0]=0;
    cs[1*mednum*group+1*group+1]=0;
    cs[1*mednum*group+1*group+2]=0;
    cs[1*mednum*group+1*group+3]=0;
    cs[2*mednum*group+1*group*group+0+group*0]=2.87880e-01;
    cs[2*mednum*group+1*group*group+0+group*1]=2.16513e-03;
    cs[2*mednum*group+1*group*group+0+group*2]=1.00130e-13;
    cs[2*mednum*group+1*group*group+0+group*3]=0.00000e+00;
    cs[2*mednum*group+1*group*group+1+group*0]=0.00000e+00;
    cs[2*mednum*group+1*group*group+1+group*1]=2.74091e-01;
    cs[2*mednum*group+1*group*group+1+group*2]=6.17082e-04;
    cs[2*mednum*group+1*group*group+1+group*3]=3.84268e-15;
    cs[2*mednum*group+1*group*group+2+group*0]=0.00000e+00;
    cs[2*mednum*group+1*group*group+2+group*1]=0.00000e+00;
    cs[2*mednum*group+1*group*group+2+group*2]=2.42573e-01;
    cs[2*mednum*group+1*group*group+2+group*3]=2.05703e-03;
    cs[2*mednum*group+1*group*group+3+group*0]=0.00000e+00;
    cs[2*mednum*group+1*group*group+3+group*1]=0.00000e+00;
    cs[2*mednum*group+1*group*group+3+group*2]=0.00000e+00;
    cs[2*mednum*group+1*group*group+3+group*3]=2.45965e-01;
    
    cs[0*mednum*group+2*group+0]=5.94040e-01;//mod
    cs[0*mednum*group+2*group+1]=1.24177e+00;
    cs[0*mednum*group+2*group+2]=1.40132e+00;
    cs[0*mednum*group+2*group+3]=2.67383e+00;
    cs[1*mednum*group+2*group+0]=0;
    cs[1*mednum*group+2*group+1]=0;
    cs[1*mednum*group+2*group+2]=0;
    cs[1*mednum*group+2*group+3]=0;
    cs[2*mednum*group+2*group*group+0+group*0]=4.39701e-01;
    cs[2*mednum*group+2*group*group+0+group*1]=1.54181e-01;
    cs[2*mednum*group+2*group*group+0+group*2]=2.25012e-05;
    cs[2*mednum*group+2*group*group+0+group*3]=1.41901e-06;
    cs[2*mednum*group+2*group*group+1+group*0]=0.00000e+00;
    cs[2*mednum*group+2*group*group+1+group*1]=1.11018e+00;
    cs[2*mednum*group+2*group*group+1+group*2]=1.22066e-01;
    cs[2*mednum*group+2*group*group+1+group*3]=8.72810e-03;
    cs[2*mednum*group+2*group*group+2+group*0]=0.00000e+00;
    cs[2*mednum*group+2*group*group+2+group*1]=2.01093e-04;
    cs[2*mednum*group+2*group*group+2+group*2]=1.06763e+00;
    cs[2*mednum*group+2*group*group+2+group*3]=3.26197e-01;
    cs[2*mednum*group+2*group*group+3+group*0]=0.00000e+00;
    cs[2*mednum*group+2*group*group+3+group*1]=0.00000e+00;
    cs[2*mednum*group+2*group*group+3+group*2]=6.10214e-02;
    cs[2*mednum*group+2*group*group+3+group*3]=2.58367e+00;//*/
    
    //入力データは最終的には関数で外部から入力するものとする
    specific_heat[0] = 14*4.184/270.03;//J/g/K
    rho[0] = 10.97*0.95;//g/cm3
    for (int i=0; i<totM*totM*group; i++) {
        Pij_n[i]=Pij[i];
    }
    //*/
    for (int i=0; i<sizes; i++) {
        cs_n[i]=cs[i];
    }
};
void Pij_dyn::c_pijcs()
{
    for (int i=0; i<totM*totM*group; i++) {
        Pij_n[i]=Pij[i];
    }
    for (int i=0; i<sizes; i++) {
        cs_n[i]=cs[i];
    }
}

void Pij_dyn::cal_dyn()
{
    if (chainef == true) {
        cout<<"cef cal"<<endl;
        c_pijcs();
        cal_vec_b_cef();
        gaus_method_vec_b();
        for (int xx=0; xx<step; xx++) {
            cal_vec_a_cef();//
            cal_nxfai_betu();//*/
            time += timestep;
            plot_phi();
        }
    }else{
        cout<<"not cef cal"<<endl;
        c_pijcs();
        cal_vec_b();
        gaus_method_vec_b();
        for (int xx=0; xx<step; xx++) {
            cal_vec_a();//
            cal_nxfai_betu();//*/
            time += timestep;
            plot_phi();
        }
    }
};
//チェーンエフェクトについて
//とりあえず、チェーンエフェクト用のベクトル作成関数作成だ！（長えとか言わない）
//共通のvec_aとかを利用、cal_vec_a、cal_vec_bを切り替え？（逆行列作成とか中性子束計算の共通化）
//cal_vec_a_cefとかを追加だ！

void Pij_dyn::cal_nxfai_betu()
{
    for (int i=0; i<size; i++) {
        phi[i]=0;
        for (int j=0; j<size; j++) {
            /*real buf = vec_bgaus[j][i];
            if (buf != 0) {
                phi[i]+=buf*vec_a[j];
            }//*/
            phi[i]+=vec_bgaus[j][i]*vec_a[j];
        }
    }
};


void Pij_dyn::reset_vec_a()
{
    for (int i=0; i<size; i++) {
        vec_a[i]=0;
    }
};
void Pij_dyn::reset_vec_b()
{
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            vec_b[i][j]=0;
        }
    }
};


void Pij_dyn::cal_vec_a()
{
    reset_vec_a();
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {//アクセス速度の問題的に、グループを外にしてメッシュを
            real no1=cs[mednum*0*group+group*med[i]+j]*vol[i];//cs[matid[i]][j]*v[i];
            vec_a[i+j*totM]-=(1-theta)*no1*phi[i+j*totM];
        }
    }
    for (int j=0; j<totM; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                int mat=med[j];
                for (int l=0; l<group; l++) {
                    real no2=Pij[j*totM+i+totM*totM*l]*vol[j]*(1-beta)*chiC[mat*group+l]*cs[mednum*1*group+mat*group+k];//csf[mat][k]*nyu[mat];
                    vec_a[i+totM*l]+=no2*(1-theta)*phi[j+totM*k];
                    real no3=Pij[j*totM+i+totM*totM*l]*vol[j]*cs[2*mednum*group+mat*group*group+k+l*group];//css[mat][k][l];
                    vec_a[i+totM*l]+=no3*(1-theta)*phi[j+totM*k];
                }
                real no4=Pij[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_a[i+totM*k]+=(1-theta)*no4*phi[j+totM*k];
                real no4_5=Pij_n[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_a[i+totM*k]+=theta*no4_5*phi[j+totM*k];
            }
        }
    }
    for (int j=0; j<freg; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                for (int m=0; m<totDNP; m++) {
                    real no5=Pij[j*totM+i+totM*totM*k]*vol[j]*DNPlam[m]*chiD[m*group+k]*fissN[m]*hosei;
                    vec_a[i+totM*k]+=no5*(1-theta)*phi[totM*group+j+freg*m];
                }
            }
        }
    }
    for (int k=0; k<group; k++) {
        for (int i=0; i<freg; i++) {
            for (int m=0; m<totDNP; m++) {
                real a=med[i];
                real No6=DNPbeta[m]*cs[mednum*1*group+a*group+k];//csf[a][k]*nyu[a];
                vec_a[totM*group+i+freg*m]+=No6*(1-theta)*phi[i+totM*k];
            }
        }
    }
    for (int i=0; i<freg; i++) {
        for (int m=0; m<totDNP; m++) {
            vec_a[totM*group+i+freg*m]+=1/timestep*phi[totM*group+i+freg*m];
            vec_a[totM*group+i+freg*m]-=(1-theta)*DNPlam[m]*phi[totM*group+i+freg*m];
        }
    }
    
};

void Pij_dyn::cal_vec_b()
{
    reset_vec_b();
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            real no1=cs_n[0*mednum*group+med[i]*group+j]*vol[i];//cs[matid[i]][j]*v[i];
            vec_b[i+j*totM][i+j*totM]+=no1*theta;
        }
    }
    for (int j=0; j<totM; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                int mat=med[j];
                for (int l=0; l<group; l++) {
                    real no2=Pij_n[j*totM+i+totM*totM*l]*vol[j]*(1-beta)*chiC[mat*group+l]*cs_n[mednum*1*group+mat*group+k];//csf[mat][k]*nyu[mat];
                    vec_b[j+totM*k][i+totM*l]-=no2*theta;
                    real no3=Pij_n[j*totM+i+totM*totM*l]*vol[j]*cs_n[2*group*mednum+mat*group*group+k+l*group];//css[mat][k][l];
                    vec_b[j+totM*k][i+totM*l]-=no3*theta;
                }
                real no4=Pij_n[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_b[j+totM*k][i+totM*k]+=no4*theta;
                real no4_5=Pij[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_b[j+totM*k][i+totM*k]+=no4_5*(1-theta);//新設
            }
        }
    }
    for (int j=0; j<freg; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                for (int m=0; m<totDNP; m++) {
                    real no5=Pij_n[j*totM+i+totM*totM*k]*vol[j]*DNPlam[m]*chiD[m*group+k]*fissN[m]*hosei;
                    vec_b[totM*group+j+freg*m][i+totM*k]-=no5*theta;
                }
            }
        }
    }
    for (int m=0; m<totDNP; m++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<freg; i++) {
                
                real a=med[i];
                real No6=DNPbeta[m]*cs_n[mednum*1*group+a*group+k];//csf[a][k]*nyu[a];
                vec_b[i+totM*k][totM*group+i+freg*m]-=No6*theta;
            }
        }
    }
    for (int i=0; i<freg; i++) {
        for (int m=0; m<totDNP; m++) {
            vec_b[totM*group+i+freg*m][totM*group+i+freg*m]+=1/timestep;
            vec_b[totM*group+i+freg*m][totM*group+i+freg*m]+=DNPlam[m]*theta;
        }
    }
    
}
//チェーンエフェクト考慮版
void Pij_dyn::cal_vec_a_cef()
{
    reset_vec_a();
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {//アクセス速度の問題的に、グループを外にしてメッシュを
            real no1=cs[mednum*0*group+group*med[i]+j]*vol[i];//cs[matid[i]][j]*v[i];
            vec_a[i+j*totM]-=(1-theta)*no1*phi[i+j*totM];
        }
    }
    for (int j=0; j<totM; j++) {//このループは、各領域での核分裂や散乱なんかを考えてる項だったりします。それにPijかけたり
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                int mat=med[j];
                for (int l=0; l<group; l++) {
                    real no2=Pij[j*totM+i+totM*totM*l]*vol[j]*(1-beta)*chiC[mat*group+l]*cs[mednum*1*group+mat*group+k];//csf[mat][k]*nyu[mat];
                    vec_a[i+totM*l]+=no2*(1-theta)*phi[j+totM*k];
                    real no3=Pij[j*totM+i+totM*totM*l]*vol[j]*cs[2*mednum*group+mat*group*group+k+l*group];//css[mat][k][l];
                    vec_a[i+totM*l]+=no3*(1-theta)*phi[j+totM*k];
                }
                real no4=Pij[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_a[i+totM*k]+=(1-theta)*no4*phi[j+totM*k];
                real no4_5=Pij_n[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_a[i+totM*k]+=theta*no4_5*phi[j+totM*k];
            }
        }
    }
    //結局、チェーンエフェクトの追加にはここから下の改修が必要よねー
    //下のさぁ、結局チェーンエフェクト考えるなら経路を考えなきゃなので複雑化の要アリ
    for (int j=0; j<freg; j++) {//ここは、各領域への遅発中性子先行核の崩壊による系への寄与っぽい。fissN=dnだ。
        //ここは結構そのままでも？各崩壊様式による中性子発生を考慮すればそれ意外（寄与する場所とか）は変更の余地が無さそう。
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                for (int m=0; m<totDNP; m++) {
                    int ngc = man.GetNuclide(m).GetChannel();
                    //real no5 = DNPlam[m]*chiD[m*group+k]*man.GetNuclide(m).GetTEN();/*
                    real no5 = 0;//Pij[j*totM+i+totM*totM*k]*vol[j]*DNPlam[m]*chiD[m*group+k]*fissN[m];
                    for (int n=0; n<ngc; n++) {
                        no5 += DNPlam[m]*chiD[m*group+k]*man.GetNuclide(m).GetBr(n)*man.GetNuclide(m).GetEmittedNeutron(n)*hosei;
                    }//*///あれ、ここが原因だったの……？いや、今度は少し過剰気味。発生中性子数が少ないっぽい。
                    //つまりは、崩壊部が問題なのか……！
                    no5 *= Pij[j*totM+i+totM*totM*k]*vol[j];
                    vec_a[i+totM*k]+=no5*(1-theta)*phi[totM*group+j+freg*m];
                }
            }
        }
    }
    //下の部分は、そのままでもいいだろうが”崩壊”による移行も考える必要があるためその項目の作成が必要
    for (int k=0; k<group; k++) {//遅発中性子先行核の生成部。とっても大事。なお、核分裂からのもの”のみ”
        //更に崩壊による寄与も必要っぽい。ある核種へ他の核種からの寄与、数密度×崩壊分岐比×崩壊定数、か？
        for (int i=0; i<freg; i++) {
            for (int m=0; m<totDNP; m++) {
                real a=med[i];
                real No6=DNPbeta[m]*cs[mednum*1*group+a*group+k];//csf[a][k]*nyu[a];
                vec_a[totM*group+i+freg*m]+=No6*(1-theta)*phi[i+totM*k];
            }
        }
    }
    for (int i=0; i<freg; i++) {//ということで、ここが崩壊による他の核種への移行部です。
        for (int m=0; m<totDNP; m++) {
            int ngc = man.GetNuclide(m).GetChannel();
            for (int n=0; n<ngc; n++) {
                real No7=DNPlam[m]*man.GetNuclide(m).GetBr(n);
                int nonb=man.GetNuclide(m).GetRRnum(n);
                vec_a[totM*group+i+freg*nonb]+=No7*(1-theta)*phi[totM*group+i+freg*m];
            }
        }
    }
    
    
    for (int i=0; i<freg; i++) {//自壊による減衰項。ここは手を入れる必要は無いっぽい？
        for (int m=0; m<totDNP; m++) {
            vec_a[totM*group+i+freg*m]+=1/timestep*phi[totM*group+i+freg*m];
            vec_a[totM*group+i+freg*m]-=(1-theta)*DNPlam[m]*phi[totM*group+i+freg*m];
        }
    }
    
};

void Pij_dyn::cal_vec_b_cef()
{
    reset_vec_b();
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            real no1=cs_n[0*mednum*group+med[i]*group+j]*vol[i];//cs[matid[i]][j]*v[i];
            vec_b[i+j*totM][i+j*totM]+=no1*theta;
        }
    }
    for (int j=0; j<totM; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                int mat=med[j];
                for (int l=0; l<group; l++) {
                    real no2=Pij_n[j*totM+i+totM*totM*l]*vol[j]*(1-beta)*chiC[mat*group+l]*cs_n[mednum*1*group+mat*group+k];//csf[mat][k]*nyu[mat];
                    vec_b[j+totM*k][i+totM*l]-=no2*theta;
                    real no3=Pij_n[j*totM+i+totM*totM*l]*vol[j]*cs_n[2*group*mednum+mat*group*group+k+l*group];//css[mat][k][l];
                    vec_b[j+totM*k][i+totM*l]-=no3*theta;
                }
                real no4=Pij_n[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_b[j+totM*k][i+totM*k]+=no4*theta;
                real no4_5=Pij[j*totM+i+totM*totM*k]*vol[j]/(vel[k]*timestep);
                vec_b[j+totM*k][i+totM*k]+=no4_5*(1-theta);//新設
            }
        }
    }
    for (int j=0; j<freg; j++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totM; i++) {
                for (int m=0; m<totDNP; m++) {
                    //real no5=Pij_n[j*totM+i+totM*totM*k]*vol[j]*DNPlam[m]*chiD[m*group+k]*fissN[m];
                    //vec_b[totM*group+j+freg*m][i+totM*k]-=no5*theta;
                    int ngc = man.GetNuclide(m).GetChannel();
                    //real no5 = DNPlam[m]*chiD[m*group+k]*man.GetNuclide(m).GetTEN();/*
                    real no5 = 0;//Pij[j*totM+i+totM*totM*k]*vol[j]*DNPlam[m]*chiD[m*group+k]*fissN[m];
                    for (int n=0; n<ngc; n++) {
                        no5 += DNPlam[m]*chiD[m*group+k]*man.GetNuclide(m).GetBr(n)*man.GetNuclide(m).GetEmittedNeutron(n)*hosei;
                    }//*/
                    no5 *= Pij[j*totM+i+totM*totM*k]*vol[j];
                    vec_b[totM*group+j+freg*m][i+totM*k]-=no5*theta;
                }
            }
        }
    }
    for (int m=0; m<totDNP; m++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<freg; i++) {
                real a=med[i];
                real No6=DNPbeta[m]*cs_n[mednum*1*group+a*group+k];//csf[a][k]*nyu[a];
                vec_b[i+totM*k][totM*group+i+freg*m]-=No6*theta;
            }
        }
    }
    for (int i=0; i<freg; i++) {//ということで、ここが崩壊による他の核種への移行部です。
        for (int m=0; m<totDNP; m++) {
            int ngc = man.GetNuclide(m).GetChannel();
            for (int n=0; n<ngc; n++) {
                real No7=DNPlam[m]*man.GetNuclide(m).GetBr(n);
                int nonb=man.GetNuclide(m).GetRRnum(n);
                //vec_a[totM*group+i+freg*nonb]+=No7*(1-theta)*phi[totM*group+i+freg*m];
                vec_b[totM*group+i+freg*m][totM*group+i+freg*nonb]-=No7*theta;

            }
        }
    }

    
    for (int i=0; i<freg; i++) {
        for (int m=0; m<totDNP; m++) {
            vec_b[totM*group+i+freg*m][totM*group+i+freg*m]+=1/timestep;
            vec_b[totM*group+i+freg*m][totM*group+i+freg*m]+=DNPlam[m]*theta;
        }
    }
    
}
//numeric
void Pij_dyn::gaus_method_vec_b()
{
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            vec_bcopy[i][j]=vec_b[i][j];
        }
    }
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            vec_bgaus[i][j]=(i==j)?1.0:0.0;
        }
    }
    real buf;
    for(int i=0;i<size;i++){
        buf=1/vec_bcopy[i][i];
        for(int j=0;j<size;j++){
            vec_bcopy[i][j]*=buf;
            vec_bgaus[i][j]*=buf;
        }
        for(int j=0;j<size;j++){
            if(i!=j){
                buf=vec_bcopy[j][i];
                for(int k=0;k<size;k++){
                    vec_bcopy[j][k]-=vec_bcopy[i][k]*buf;
                    vec_bgaus[j][k]-=vec_bgaus[i][k]*buf;
                }
            }
        }
    }}
void Pij_dyn::non_heatcal()
{
    for (int i=0; i<=sizes; i++) {
        cs_n[i]=cs[i];
    }
    for (int i=0; i<totM*totM*group; i++) {
        Pij_n[i]=Pij[i];
    }
}



//void Pij_dyn::input_cst(int matid,int ener,real data){cs[0*mednum*group+matid*group+ener]=data;}
//void Pij_dyn::input_csf(int matid,int ener,real data){cs[1*mednum*group+matid*group+ener]=data;}
//void Pij_dyn::input_css(int matid,int i,int j,real data){cs[2*mednum*group+matid*group*group+i+group*j]=data;}
//void Pij_dyn::input_pij(int i,int j,int ene,real data){Pij[i*totM+j+totM*totM*ene]=data;}
//void Pij_dyn::input_vol(int i,real data){vol[i]=data;}
//void Pij_dyn::input_firstoutput(real data){firstoutput = data;}
//void Pij_dyn::input_specific_heat(int i,real data){specific_heat[i] = data;}
//void Pij_dyn::input_rho(int i,real data){rho[i] = data;}
//void Pij_dyn::input_heattrans(int i,real data){heattrans[i] = data;}
//void Pij_dyn::input_vel(int i,real data){vel[i] = data;}

void Pij_dyn::input_beta(real data)
{
    beta = data;
    DNPbeta[0] = beta*0.033;
    DNPbeta[1] = beta*0.219;
    DNPbeta[2] = beta*0.196;
    DNPbeta[3] = beta*0.395;
    DNPbeta[4] = beta*0.115;
    DNPbeta[5] = beta*0.042;
    //でもこれって、実質的にはU-235燃料のっすよね。
}

void Pij_dyn::input_steps(real time,real tstep)
{
    Totaltime = time;
    timestep = tstep;
    step = Totaltime/timestep;
    //step = steps;
    //timestep = Totaltime/step;
}




void Pij_dyn::setup_Pij()
{
    //本来、ここで入力はなんか違わね？って思ってます。
    //ひとまず、下記の機能をここに移す模様。ひとまずradiusやregidは用意済み
    /*int ring=4;
     real radius[]={4.,3.,2.,1.};
     int regid[]={3,2,1,0};
     IrregularGeometryInformation igi;
     igi.AddCircleRing(ring,radius,regid);
     Pij_dyn test(ring,3,2);
     
     //igi.WriteGnuplotFile(0.01);
     int energy = 2;
     TrajectorySet sys;
     sys.PutBoundaryCondition(White);
     sys.CalTrajectory(igi,1,0.001,45.); // anguler division/width/angle*/
    
}
//定常計算用ツール、ベタ移植。これらの
void Pij_dyn::cssuesankaku(){
    for (int k=0; k<mednum; k++) {
        for (int i=1; i<group; i++) {
            for (int j=0; j<i; j++) {
                //css[k][i][j]=0;
                cs[2*mednum*group+k*group*group+i+group*j]=0;
            }
        }
    }
}
void Pij_dyn::setup_start(){
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            phi[i+totM*j]=1.;
        }
    }
    cal_fiss();
    real imoyoukan=1/get_totalfiss();
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            phi[i+totM*j]=imoyoukan;
        }
    }
    cal_fiss();
}
void Pij_dyn::setup_nxPhi(){
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            phi[i+totM*j]=phi_n[i+totM*j];
            phi_n[i+totM*j]=0;
        }
    }
}
void Pij_dyn::cal_fiss(){
    for (int i=0; i<totM; i++) {
        fission_S[i]=0;
    }
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            //fission_S[i]+=fai[i][j]*calcsf[i][j];//csf[matid[i]][j]*nyu[matid[i]];
            fission_S[i]+=phi[i+totM*j]*cs[1*mednum*group+med[i]*group+j];
        }
    }
}
void Pij_dyn::cal_scat(){
    for (int i=0; i<totM; i++) {
        for (int j=0; j<group; j++) {
            //scat_S[i][j]=0;
            scat_S[i+totM*j]=0;
        }
    }
    for (int i=0; i<totM; i++) {//計算対象の位置
        for (int j=0; j<group; j++) {//計算対象のエネルギー
            for (int l=0; l<group; l++) {//飛んでくるもと
                //scat_S[i][j]+=calcss[i][l][j]*fai[i][l];//css[matid[i]][l][j]
                scat_S[i+totM*j]+=cs[2*mednum*group+med[i]*group*group+l+group*j]*phi[i+totM*l];//css[matid[i]][l][j]
            }
        }
    }
}
double Pij_dyn::get_totalfiss(){
    real youkan=0;
    for (int i=0; i<totM; i++) {
        youkan+=fission_S[i];
    }
    return youkan;
}

void Pij_dyn::solver(){
    //real hanteiki=1.0e-10;
    //keff=1;//初回keff
    int count = 0;
    real check_totalfis;
    real kefftemp = 0;
    setup_start();//
    do{
        real hantemp=0;
        //fissionの処理を先にね？
        kefftemp=1/keff;
        check_totalfis=get_totalfiss();
        cal_scat();
        for (int j=0; j<group; j++) {
            for (int i=0; i<totM; i++) {
                real vtemp=vol[i];
                //real keitemp1=1/(cs[matid[i]][j]*vtemp);
                real keitemp1=1/(cs[0*mednum*group+med[i]*group+j]*vtemp);
                real keitemp2=0;
                for (int k=0; k<totM; k++) {
                    real zuzui=0;
                    zuzui=fission_S[k]*chiC[med[k]*group+j]*kefftemp+scat_S[k+totM*j];
                    //keitemp2+=P[k][i][j]*zuzui*vol[k];
                    keitemp2+=Pij[k*totM+i+totM*totM*j]*zuzui*vol[k];
                }
                //nxfai[i][j]=keitemp2*keitemp1;
                phi_n[i+totM*j]=keitemp2*keitemp1;
		//cout<<phi_n[i+totM*j]<<" "<<j<<" "<<i<<"\n";
            }
        }
        setup_nxPhi();
        cal_fiss();
        hantemp=get_totalfiss();
        //nextKeff
        hantei=fabs((hantemp-keff)*kefftemp);
        keff=hantemp;
        count++;
    }while(hanteiki<hantei);
    //cout.precision(10);
    cout<<"keff : "<<keff<<endl;
    cout<<"ene : "<<group<<endl;
    cout<<"DNP : "<<totDNP<<endl;
    //cout<<"roop : "<<count<<endl;//*/
}//*/

void Pij_dyn::plot_phi()
{
    fs1<<time;
    fs2<<time;
    fs3<<time;
    //cout<<time;
    for (int i=0; i<1; i++) {
        fs1<<" , "<<phi[i+totM*0];
        fs2<<" , "<<phi[i+totM*1];
        //cout<<" , "<<phi[i+totM*1];
    }
    for (int i=0; i<1; i++) {
        fs3<<" , "<<phi[totM*group+1+freg*i];
    }
    fs1<<endl;
    fs2<<endl;
    fs3<<endl;
    real output_en = 0;
    for (int i=0; i<freg; i++) {
        for (int j=0; j<group; j++) {
            output_en+=phi[i+totM*j]*cs[mednum*1*group+med[i]*group+j];
        }
    }
    fs9<<time<<" "<<output_en<<endl;

    //cout<<endl;
}

void Pij_dyn::test()
{
    /*for (int i=0; i<totM; i++) {
     cout<<"vol["<<i<<"] = "<<vol[i]<<endl;
     }*/
    /*for (int i=0; i<totDNP; i++) {
        cout<<phi[totM*group+0+totM*i]<<endl;
    }
    /*for (int i=0; i<totDNP-1; i++) {
     cout<<phi[totM*group+0+totM*i]<<" , ";
     }
     cout<<phi[totM*group+0+totM*(totDNP-1)]<<endl;//*/
    /*for (int i=0; i<totM; i++) {
     cout<<Pij[i]<<endl;
     }
     cout<<"-------"<<endl;
     for (int i=0; i<totM; i++) {
     cout<<Pij_n[i]<<endl;
     }//*/
    cout<<"beta = "<<beta<<endl;
    for (int i=0;i<size; i++) {
        //fs7<<i<<" "<<mmpa_ans[i]/phi[i]<<endl;
        fs7<<mmpa_ans[i]<<" "<<phi[i]<<endl;
    }
    /*for (int i=0; i<totDNP; i++) {
        cout<<"DNPbeta["<<i<<"] : "<<DNPbeta[i]<<endl;
    }//*/
    //cout<<phi[0]<<" , "<<phi[1]<<" , "<<phi[2]<<" , "<<phi[3]<<endl;
    //cout<<han_sta<<endl;
    //cout<<"keff = "<<keff<<endl;
}


void Pij_dyn::setup_dp()
{
    if (chainef == true) {//さて、ここにループ構造を作りFPの初期値を決めよう
        dim1d fpcal;//前のFP量
        dim1d fps;//総発生量
        fpcal.resize(totDNP);
        fps.resize(totDNP);
        int miru = 651;
        //real forcheck;
        for (int i=0; i<freg; i++) {
            for (int i=0; i<totDNP; i++) {
                fpcal[i]=0.;
                fps[i]=0.;
            }
            for (int xx=0; xx<20; xx++) {//とりあえず20回だ
                real bbuf = 0;
                real now  = 0;
                for (int j=0; j<group; j++) {//bbufで核分裂ソースの作成
                    bbuf+=cs[1*mednum*group+med[i]*group+j]*phi[i+totM*j];
                }
                //forcheck=bbuf*beta;
                //cout<<fps[miru]<<endl;
                for (int m=0; m<totDNP; m++) {//他の核種からの寄与項、何度も繰り返してたら定常になるでしょう。
                    int ngc = man.GetNuclide(m).GetChannel();
                    for (int n=0; n<ngc; n++) {
                        int nonb=man.GetNuclide(m).GetRRnum(n);
                        real ttt = DNPlam[m]*man.GetNuclide(m).GetBr(n)*fpcal[m];
                        fps[nonb]+=ttt;
                        /*if (nonb == miru &&xx==99) {
                            cout<<n<<endl;
                            cout<<m<<" , "<<ttt<<" , "<<fps[nonb]<<endl;
                        }//*/
                    }
                }
                //mmpa_a[totM*group+m+freg*i][totM*group+m+freg*nonb]=DNPlam[i]*man.GetNuclide(i).GetBr(n);
                //番号,分岐によるアレ、生成全体、崩壊定数、数密度、計算結果数密度
                fs4<<xx<<" , "<<fps[miru]<<" ,"<<fps[miru]+bbuf*DNPbeta[miru]<<","<<bbuf*DNPbeta[miru]<<","<<DNPlam[miru]<<" , "<<(fps[miru]+bbuf*DNPbeta[miru])/DNPlam[miru]<<" , ";
                //cout<<bbuf<<endl;;
                for (int m=0; m<totDNP; m++) {//核分裂からの寄与を計算
                    real a = DNPbeta[m];
                    real b = DNPlam[m];
                    real check = bbuf*a;
                    //cout<<"---"<<endl;
                    //cout<<fps[m]<<endl;
                    fps[m]+=check;
                    //cout<<fps[m]<<endl;
                    if (b!=0) {
                        fps[m]/=b;
                    }else{//崩壊しなきゃ、溜まるだけだしね。
                        fps[m]=0.;
                    }
                    //cout<<fps[m]<<endl;
                    phi[totM*group+i+freg*m]=fps[m];
                    fpcal[m]=fps[m];
                }
                for (int m=0; m<totDNP; m++) {
                    fps[m]=0.;
                }//*/
                fs4<<fpcal[miru]<<endl;
                //fs4<<xx<<" , "<<fpcal[570]<<endl;
            }
            /*real cyousei=0.;
            for (int m=0; m<totDNP; m++) {
                cyousei+=phi[totM*group+i+freg*m]*fissN[m]*DNPlam[m];
            }
            cout<<"emit neutron :"<<cyousei<<endl;
            cout<<"emit  fiss   :"<<forcheck<<endl;//*/

        }
        /*real buf = beta;//これが原因でしたね！滅びろ！そして、残りは高速中性子と651番遅発中性子先行核の異常な減少？
        for (int i=0; i<freg; i++) {
            real bbuf = 0;
            real now  = 0;
            for (int j=0; j<group; j++) {
                bbuf+=cs[1*mednum*group+med[i]*group+j]*phi[i+totM*j];
            }
            for (int m=0; m<totDNP; m++) {
                real a = DNPbeta[m];
                real b = DNPlam[m];
                real check;
                if (a==0||b==0) {
                    check = 0;
                }else{
                    check = bbuf*a/b;
                }
                phi[totM*group+i+freg*m]=check;
                now +=phi[totM*group+i+freg*m]*DNPlam[m]*fissN[m];
            }
        }//*/

    }else{
        for (int i=0; i<freg; i++) {
            real bbuf = 0;
            real now  = 0;
            for (int j=0; j<group; j++) {
                bbuf+=cs[1*mednum*group+med[i]*group+j]*phi[i+totM*j];
            }
            cout<<"freg["<<i<<"],total fiss = "<<bbuf<<endl;
            cout<<"freg["<<i<<"],delay(in) = "<<bbuf*beta*nyt<<endl;
            for (int m=0; m<totDNP; m++) {
                real a = DNPbeta[m];
                real b = DNPlam[m];
                real check;
                if (a==0||b==0) {
                    check = 0;
                }else{
                    check = bbuf*a/b;
                }
                phi[totM*group+i+freg*m]=check;
                now +=phi[totM*group+i+freg*m]*DNPlam[m]*fissN[m];
            }
            cout<<"freg["<<i<<"],delay(out) = "<<now<<endl;
            cout<<"freg["<<i<<"],delay(out/in) = "<<now/(bbuf*beta*nyt)<<endl;
            cout<<"freg["<<i<<"],beta(calculated) = "<<now/bbuf<<endl;
        }
    }//*/
}

//色々あるが、何も言うまい
void Pij_dyn::make_mmpa_a()
{
    for (int k=0; k<group; k++) {
        for (int i=0; i<totM; i++) {
            for (int j=0; j<totM; j++) {
                for (int l=0; l<group; l++) {
                    real No1 = 0;
                    int mat=med[j];
                    No1+= cs[2*mednum*group+mat*group*group+l+group*k];//l to kの弾性散乱
                    No1+= (1-beta)*chiC[mat*group+k]*cs[mednum*1*group+med[j]*group+l];//核分裂項
                    //csf[] En:l to k,area: j to i
                    No1*= Pij[j*totM+i+totM*totM*k]*vol[j];
                    mmpa_a[j+totM*l][i+totM*k]+= No1;
                }
            }
            mmpa_a[i+totM*k][i+totM*k]-= cs[0*mednum*group+med[i]*group+k]*vol[i];
            //流石にここまでは問題ないっぽい？
            //遅発中性子先行核からの寄与をここに
        }
    }
    for (int m=0; m<freg; m++) {
        for (int i=0; i<totM; i++) {
            for (int k=0; k<group; k++) {
                real pijnow = Pij[m*totM+i+totM*totM*k]*vol[m];//m to i energy:k
                for (int n=0; n<totDNP; n++) {
                    mmpa_a[totM*group+m+freg*n][i+totM*k] += pijnow*chiD[n*group+k]*fissN[n]*DNPlam[n]*hosei;
                }
            }
        }
    }//*/

    //作成するもの
    //遅発中性子先行核への寄与、バグ取り完了
    for (int m=0; m<freg; m++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totDNP; i++) {
                    mmpa_a[m+totM*k][totM*group+m+freg*i]+=DNPbeta[i]*cs[mednum*1*group+med[m]*group+k];
            }
        }
    }
    //phi[j+totM*k],j領域k群の中性子束
    //phi[totM*group+j+freg*m],j領域のm番DNP
    //mmpa_a[j][i] jからiへの寄与
    //mmpa_a[j+totM*k][totM*group+i+freg*m]j領域k群の中性子束からi領域のm番DNPへ
    //Pij[i*totM+j+totM*totM*ene] iからjへの移行率。エネルギーはene
    //chaineffectによる色々
    //cs[1*mednum*group+med[i]*group+j]*phi[i+totM*j]*DNPbeta[miru]
    if (chainef == true) {
        for (int m=0; m<freg; m++) {
            for (int i=0; i<totDNP; i++) {//ここで、chainエフェクトを
                int ngc=man.GetNuclide(i).GetChannel();
                for (int n=0; n<ngc; n++) {
                    int nonb=man.GetNuclide(i).GetRRnum(n);
                    //cout<<i<<" , "<<nonb<<endl;
                    mmpa_a[totM*group+m+freg*i][totM*group+m+freg*nonb]+=DNPlam[i]*man.GetNuclide(i).GetBr(n);
                }
            }
        }
    }//*/
    
    //遅発中性子先行核内の色々
    for (int i=0; i<freg; i++) {
        for (int m=0; m<totDNP;m++) {
            real buf = DNPlam[m];
            mmpa_a[totM*group+i+freg*m][totM*group+i+freg*m] = -buf;
            //mmpa_a[totM*group+i+freg*m][totM*group+i+freg*m]-=1/timestep;
            /*if (buf == 0) {
                for (int j=0; j<size; j++) {
                    mmpa_a[j][totM*group+i+freg*m] = 0;
                    cout<<"zap!zap!"<<endl;
                }
            }else{
                mmpa_a[totM*group+i+freg*m][totM*group+i+freg*m] = -buf;
            }//*/
            //mmpa_a[totM*group+i+freg*m][totM*group+i+freg*m] -= DNPlam[m];
        }
    }
    
    //ここに、ちょっと崩壊に寄る分の確認を。
    for (int j=0; j<freg; j++) {
    real tetetetete=0;
        real fff=0;
        for (int k=0; k<group; k++) {
            fff+= beta*cs[mednum*1*group+med[j]*group+k]*phi[j+totM*k]*vol[j];

        }
    for (int i=0; i<totDNP; i++) {
        for (int k = 0; k<totM; k++) {
        tetetetete+=fissN[i]*DNPlam[i]*phi[totM*group+j+freg*i]*Pij[j*totM+k+totM*totM*0]*vol[j];
        }
    }
        cout<<"tot fiss n :"<<fff<<endl;
    cout<<"tot emit n :"<<tetetetete<<endl;
    }
    /*for (int i=0; i<group*totM; i++) {
        for (int j=0; j<group*totM; j++) {
            mmpa_a[i][j]*=1;
        }
    }//*/
}
void Pij_dyn::make_mmpa_b()
{
    for (int k=0; k<group; k++) {
        for (int i=0; i<totM; i++) {
            for (int j=0; j<totM; j++) {
                mmpa_b[j+totM*k][i+totM*k]+=Pij[j*totM+i+totM*totM*k]*vol[j]/(vel[k]);//jtoi
            }
        }
    }
    for (int i=0; i<freg; i++) {
        for (int j=0; j<totDNP; j++) {
            mmpa_b[totM*group+i+freg*j][totM*group+i+freg*j]+=1.;
        }
    }
    
}
void Pij_dyn::cal_mmpa_a()
{
    make_mmpa_a();
    make_mmpa_b();
    real anooo=0;
    for (int i=0; i<size; i++) {
        mmpa_an[i]=0;
        for (int j=0; j<size; j++) {
            real buf = mmpa_a[j][i]*phi[j];
            if (i <=totM) {
                if (group*totM<=j) {
                    anooo+=buf;
                }
            }
            /*if (i == 1316) {
                if (buf !=0) {
                    fs7<<j<<" "<<buf<<endl;
                }
            }//*/
            mmpa_an[i]+=buf;//*/
            //mmpa_an[i]+=mmpa_a[j][i]*phi[j];
        }
        //fs5<<mmpa_an[i]<<" , "<<phi[i]<<endl;
        fs5<<i<<" "<<mmpa_an[i]<<endl;
        //cout<<mmpa_an[i]<<endl;/*
        //cout<<phi[i]<<endl;//*/
    }
    cout<<"calculated :"<<anooo<<endl;
}
void Pij_dyn::mmpa_clean()//
{
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            mmpa_a[i][j]=0;
            mmpa_b[i][j]=0;
            mmpa_bg[i][j]=0;
            mmpa_mat[i][j]=0.;
        }
    }
}
void Pij_dyn::cal_mmpa()//
{//MMPA法の計算。
    mmpa_clean();
    make_mmpa_a();
    make_mmpa_b();
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            mmpa_bg[i][j]=(i==j)?1.0:0.0;
        }
    }
    real buf;
    int ssize = totM*group;/*
    int ssize = size;//*/
    for(int i=0;i<ssize;i++){
        buf=1/mmpa_b[i][i];
        for(int j=0;j<ssize;j++){
            mmpa_b[i][j]*=buf;
            mmpa_bg[i][j]*=buf;
        }
        for(int j=0;j<ssize;j++){
            if(i!=j){
                buf=mmpa_b[j][i];
                for(int k=0;k<ssize;k++){
                    mmpa_b[j][k]-=mmpa_b[i][k]*buf;
                    mmpa_bg[j][k]-=mmpa_bg[i][k]*buf;
                }
            }
        }
    }
    cout<<"-----------"<<endl;
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            mmpa_mat[i][j]=0.;
        }
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            buf=mmpa_bg[i][j];
            if (buf!=0) {
                for (int k=0; k<size; k++) {
                    mmpa_mat[k][j]+=mmpa_bg[i][j]*mmpa_a[k][j];
                }
            }
        }
    }

            //mmpa_mat[i][j]=mmpa_a[i][j];
    //ちょっとしたテスト
    //これさあ……即発跳躍が出来てないんじゃない？
    /*for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            mmpa_mat[i][j]*=timestep;
        }
    }
    dim1d calcu;
    calcu.resize(size);
    for (int i=0; i<size; i++) {
        mmpa_mat[i][i]+=1.;
        calcu[i]=phi[i];
    }
    for (int j=0; j<step; j++) {
        for (int i=0; i<size;i++) {
            mmpa_ans[i]=0;
            for (int k=0; k<size; k++) {
                mmpa_ans[i]+=mmpa_mat[k][i]*calcu[k];
            }
        }
        for (int i=0; i<size; i++) {
            calcu[i]=mmpa_ans[i];
        }
        fs9<<mmpa_ans[0]<<endl;
    }//*/
    
    //MMPA計算場
    GroupData2D kawamoto(size,size); // burnup matrix
    GroupData1D kajihara(size); // initial values 
    GroupData1D GC(size); // 
    for (int i=0; i<size; i++) {
        kajihara.put_data(i,phi[i]);
        for (int j=0; j<size; j++) {
            kawamoto.put_data(i,i,j,j,mmpa_mat[j][i]);
        }
    }

    for(int i=0;i<26;i++){
      cout<<i<<" "<<kawamoto.get_dat(i,i)*timestep<<"\n";
    };
    //exit(0);

    /*
    GroupData1D eigen;
    kawamoto.CalEigenvaluesByLF(eigen);
    eigen.show_self();
    */
    //exit(0);


    GroupData2D tmp=kawamoto.CalMatrixExponentialByMMPA32(timestep);
    //GroupData2D tmp=kawamoto.CalMatrixExponentialByChebyshev16(timestep);
    //tmp.show_self(false); exit(0);

    for (int xx=0; xx<=step; xx++) {
        real output_en = 0;
        if (use_mmpa == true) {
	  GC = kawamoto.CalMatrixExponentialByMMPA32(kajihara, xx*timestep);//mmpa
          // GC = kawamoto.CalMatrixExponentialByMMPA18(kajihara, xx*timestep);//mmpa
        }else{
            GC = kawamoto.CalMatrixExponentialByChebyshev16(kajihara, xx*timestep);
        }
        fs1<<nowtime<<" ";
        fs2<<nowtime<<" ";
        fs3<<nowtime<<" ";
    for (int i=0;i<totM-1; i++) {
        //cout<<GC.get_dat(i)<<endl;
        fs1<<GC.get_dat(i)<<" ";
        fs2<<GC.get_dat(totM+i)<<" ";
        //fs7<<phi[i]<<" , "<<GC.get_dat(i)<<endl;
    }//*/
        fs3<<GC.get_dat(totM*group+0+freg*0)<<" ";
        fs3<<GC.get_dat(totM*group+1+freg*0)<<endl;
        fs1<<GC.get_dat(totM)<<endl;
        fs2<<GC.get_dat(2*totM)<<endl;

        for (int i=0; i<freg; i++) {
            for (int j=0; j<group; j++) {
                output_en+=GC.get_dat(i+totM*j)*cs[mednum*1*group+med[i]*group+j];
            }
        }
        fs9<<nowtime<<" "<<output_en<<endl;
        cout<<nowtime<<" "<<output_en<<endl;

        nowtime+=timestep;
    }

    for (int i=0; i<size; i++) {
        phi[i]=GC.get_dat(i);
    }
}

void Pij_dyn::cal_mmpa_n()//
{//MMPA法の計算。
    mmpa_clean();
    make_mmpa_a_n();
    make_mmpa_b();
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            mmpa_bg[i][j]=(i==j)?1.0:0.0;
        }
    }
    real buf;
    int ssize = totM*group;/*
                            int ssize = size;//*/
    for(int i=0;i<ssize;i++){
        buf=1/mmpa_b[i][i];
        for(int j=0;j<ssize;j++){
            mmpa_b[i][j]*=buf;
            mmpa_bg[i][j]*=buf;
        }
        for(int j=0;j<ssize;j++){
            if(i!=j){
                buf=mmpa_b[j][i];
                for(int k=0;k<ssize;k++){
                    mmpa_b[j][k]-=mmpa_b[i][k]*buf;
                    mmpa_bg[j][k]-=mmpa_bg[i][k]*buf;
                }
            }
        }
    }
    cout<<"-----------"<<endl;
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            mmpa_mat[i][j]=0.;
        }
    }
    cout.precision(5);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            buf=mmpa_bg[i][j];
            if (buf!=0) {
                for (int k=0; k<size; k++) {
                    mmpa_mat[k][j]+=mmpa_bg[i][j]*mmpa_a[k][j];
                }
            }
        }
    }
    cout.precision(15);
    //MMPA計算場
    GroupData2D kawamoto(size,size);
    GroupData1D kajihara(size);
    GroupData1D GC(size);
    for (int i=0; i<size; i++) {
        kajihara.put_data(i,phi[i]);
        for (int j=0; j<size; j++) {
            kawamoto.put_data(i,i,j,j,mmpa_mat[j][i]);
        }
    }
    for (int xx=0; xx<step; xx++) {
        real output_en = 0;
        if (use_mmpa == true) {
            GC = kawamoto.CalMatrixExponentialByMMPA32(kajihara, xx*timestep);//mmpa
        }else{
            GC = kawamoto.CalMatrixExponentialByChebyshev16(kajihara, xx*timestep);
        }

        fs1<<nowtime<<" ";
        fs2<<nowtime<<" ";
        fs3<<nowtime<<" ";
        for (int i=0;i<totM-1; i++) {
            //cout<<GC.get_dat(i)<<endl;
            fs1<<GC.get_dat(i)<<" ";
            fs2<<GC.get_dat(totM+i)<<" ";
            //fs7<<phi[i]<<" , "<<GC.get_dat(i)<<endl;
        }//*/
        fs1<<GC.get_dat(totM)<<endl;
        fs2<<GC.get_dat(2*totM)<<endl;
        fs3<<GC.get_dat(totM*group+0+freg*0)<<" ";
        fs3<<GC.get_dat(totM*group+1+freg*0)<<" ";//koko
        fs3<<GC.get_dat(totM*group+2+freg*0)<<" ";
        fs3<<GC.get_dat(totM*group+3+freg*0)<<" ";//koko
        fs3<<GC.get_dat(totM*group+4+freg*0)<<" ";
        fs3<<GC.get_dat(totM*group+5+freg*0)<<endl;//koko
        for (int i=0; i<freg; i++) {
            for (int j=0; j<group; j++) {
                output_en+=GC.get_dat(i+totM*j)*cs[mednum*1*group+med[i]*group+j];
            }
        }
        fs9<<nowtime<<" "<<output_en<<endl;

        nowtime+=timestep;
    }
    for (int i=0; i<size; i++) {
        phi[i]=GC.get_dat(i);
    }
}

void Pij_dyn::make_mmpa_a_n()
{
    for (int k=0; k<group; k++) {
        for (int i=0; i<totM; i++) {
            for (int j=0; j<totM; j++) {
                for (int l=0; l<group; l++) {
                    real No1 = 0;
                    int mat=med[j];
                    No1+= cs[2*mednum*group+mat*group*group+l+group*k];//l to kの弾性散乱
                    No1+= (1-beta)*chiC[mat*group+k]*cs[mednum*1*group+med[j]*group+l];//核分裂項
                    //csf[] En:l to k,area: j to i
                    No1*= Pij[j*totM+i+totM*totM*k]*vol[j];
                    mmpa_a[j+totM*l][i+totM*k]+= No1;
                }
            }
            mmpa_a[i+totM*k][i+totM*k]-= cs[0*mednum*group+med[i]*group+k]*vol[i];
            //流石にここまでは問題ないっぽい？
            //遅発中性子先行核からの寄与をここに
        }
    }
    for (int m=0; m<freg; m++) {
        for (int i=0; i<totM; i++) {
            for (int k=0; k<group; k++) {
                real pijnow = Pij[m*totM+i+totM*totM*k]*vol[m];//m to i energy:k
                for (int n=0; n<totDNP; n++) {
                    mmpa_a[totM*group+m+freg*n][i+totM*k] += pijnow*chiD[n*group+k]*fissN[n]*DNPlam[n]*hosei;
                }
            }
        }
    }//*/
    
    //作成するもの
    //遅発中性子先行核への寄与、バグ取り完了
    for (int m=0; m<freg; m++) {
        for (int k=0; k<group; k++) {
            for (int i=0; i<totDNP; i++) {
                mmpa_a[m+totM*k][totM*group+m+freg*i]+=(DNPbeta[i]*br_cal[k]+u8y[i]*(1-br_cal[k]))*cs[mednum*1*group+med[m]*group+k];
            }
        }
    }
    //phi[j+totM*k],j領域k群の中性子束
    //phi[totM*group+j+freg*m],j領域のm番DNP
    //mmpa_a[j][i] jからiへの寄与
    //mmpa_a[j+totM*k][totM*group+i+freg*m]j領域k群の中性子束からi領域のm番DNPへ
    //Pij[i*totM+j+totM*totM*ene] iからjへの移行率。エネルギーはene
    //chaineffectによる色々
    //cs[1*mednum*group+med[i]*group+j]*phi[i+totM*j]*DNPbeta[miru]
    if (chainef == true) {
        for (int m=0; m<freg; m++) {
            for (int i=0; i<totDNP; i++) {//ここで、chainエフェクトを
                int ngc=man.GetNuclide(i).GetChannel();
                for (int n=0; n<ngc; n++) {
                    int nonb=man.GetNuclide(i).GetRRnum(n);
                    //cout<<i<<" , "<<nonb<<endl;
                    mmpa_a[totM*group+m+freg*i][totM*group+m+freg*nonb]+=DNPlam[i]*man.GetNuclide(i).GetBr(n);
                }
            }
        }
    }//*/
    
    //遅発中性子先行核内の色々
    for (int i=0; i<freg; i++) {
        for (int m=0; m<totDNP;m++) {
            real buf = DNPlam[m];
            mmpa_a[totM*group+i+freg*m][totM*group+i+freg*m] = -buf;
        }
    }
    
    //ここに、ちょっと崩壊に寄る分の確認を。
    for (int j=0; j<freg; j++) {
        real tetetetete=0;
        real fff=0;
        for (int k=0; k<group; k++) {
            fff+= beta*cs[mednum*1*group+med[j]*group+k]*phi[j+totM*k]*vol[j];
            
        }
        for (int i=0; i<totDNP; i++) {
            for (int k = 0; k<totM; k++) {
                tetetetete+=fissN[i]*DNPlam[i]*phi[totM*group+j+freg*i]*Pij[j*totM+k+totM*totM*0]*vol[j];
            }
        }
        cout<<"tot fiss n :"<<fff<<endl;
        cout<<"tot emit n :"<<tetetetete<<endl;
    }
}
void Pij_dyn::mmpa_switch(bool mmpa)
{
    use_mmpa = mmpa;
}

