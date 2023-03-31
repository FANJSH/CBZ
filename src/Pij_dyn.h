#ifndef PIJDYN
#define PIJDYN

#include <fstream>
#include <vector>
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include"FunctionTable.h"
#include"GeomVector.h"
#include"Numeric.h"
#include"GeneralSystem.h"
#include "BurnupChainGenerator.h"
#include "GroupData.h"
#include "IrregularGeometryInformation.h"
#include <math.h>
#include "LibData.h"
#include "Medium.h"
#include "OnePointCalculator.h"
#include "PJI_trajectoryset.h"
#include "PJI_system.h"
#include "Burnup.h"
#include "SelfShieldingCalculator.h"
#include "GeneralSystem.h"
#include "BurnupChainGenerator.h"


using namespace std;
typedef vector< vector<real> > dim2d;
typedef vector< vector< vector<real> > > dim3d;
typedef vector< real > dim1d;
typedef vector< vector<int> > dim2di;
typedef vector< vector< vector<int> > > dim3di;
typedef vector<int> dim1di;

class Pij_dyn{
 private:
    int mednum;//Matidの数
    int group;//群数
    int totM;//全領域数
    int totDNP;//全遅発中性子先行核数
    int allmat;//現在Initializeで定義中、気をつけてー
    int size;//resize用の数
    int sizes;//同上
    dim1d DNPlam;//先行核の崩壊定数
    dim1d DNPbeta;//先行核への分岐比
    dim1d vec_a;//ここから計算用
    dim2d vec_b;//
    dim2d vec_bcopy;//
    dim2d vec_bgaus;//
    dim1d med;//領域におけるMatIDを示す
    dim1d phi;//Φ、中性子束ですね
    dim1d phi_n;
    dim1d chiC;//useDataで定義中
    dim1d chiD;//同上
    dim1d vol;//未定義
    dim1d Pij;//未定義
    dim1d Pij_n;//next
    real beta;
    dim1d cs;//未定義
    dim1d cs_n;//next
    dim1d cs_c;//Pijcal用
    dim1d Pij_c;//同上
    real theta;//現在cxx内UseDataにて定義中、後で適当にね？
    dim1d vel;//同上、二群なので後で色々と考えようねー
    real timestep;//時間幅、同上
    real Totaltime;//全時間、同上
    int step;//ステップ数、ここまで同上
    real time;//現時刻
    int dynstep;//動特性計算におけるステップ数（熱計算等の回数、基本はTotaltime/timestep）、ここから下は現在未実装
    int fuelreg;//ここと下二つは温度計算用、各領域の領域数
    int gapreg;//
    int cladreg;//
    dim1d temperature;//温度
    dim1d heattrans;//''各領域熱伝導率
    dim1d heatoutput;//''各領域熱出力
    dim1d specific_heat;//比熱
    dim1d rho;//密度
    dim1d cal_rad;//計算用の云々
    real firstoutput;//初期出力
    real hanteiki;//主に、定常計算用
    real keff;//同上
    dim1d radius;
    dim1d pijregid;
    dim1d fission_S;//定常計算フィッションソース、位置
    dim1d scat_S;//定常計算散乱項、位置＋totM*群
    real hantei;
    dim1d fissN;//
    real hosei;
    real nyd;
    real nyt;
    real toriaezu;
    BCGManager man;
    real nowtime;
    bool chainef;//独立収率か否か
    int freg;//燃料領域数
    dim2d mmpa_a;//MMPA用？
    dim1d mmpa_an;
    dim2d mmpa_b;
    dim2d mmpa_bg;
    dim2d mmpa_mat;
    dim1d mmpa_ans;
    IrregularGeometryInformation igi;
    TrajectorySet sys;
    real enrich_f;
    dim1d br_cal;//sigf_total内のu235の割合。3種以上の収率を用いる場合？拡張が必要ですね。
    int br_seed;
    dim1d u8y;
    Medium fueld;//燃料データ保持
    Medium cladd;//被覆材データ
    Medium modd;//減速材データ
    dim1di fpid;
    dim1di fpid_nspect;
    dim2d nspect;
    dim1d mat_D;//:D
    int fp_spect_non;
    bool use_mmpa;//mmpa?
 public:
  //Trajectory(){numreg=-1;};
  Pij_dyn(int ttotM,int allmat,int ene,int dp,int ffreg,bool chain,int brr){Initialize(ttotM,allmat,ene,dp,ffreg,chain,brr);};
  void Initialize(int ttotM, int allmat,int ene,int dp,int ffreg,bool chain,int brr);
    void UseData();
    void cal_nxfai_betu();
    void cal_dyn();
    void reset_vec_a();
    void reset_vec_b();
    void cal_vec_a();
    void cal_vec_b();
    void cal_vec_a_cef();
    void cal_vec_b_cef();
    void gaus_method_vec_b();
    void input_cst(int matid,int ener,real data){cs[0*mednum*group+matid*group+ener]=data;}
    void input_csf(int matid,int ener,real data){cs[1*mednum*group+matid*group+ener]=data;}
    void input_css(int matid,int i,int j,real data){cs[2*mednum*group+matid*group*group+i+group*j]=data;}//i to j (ener)
    void input_pij(int i,int j,int ene,real data){Pij[i*totM+j+totM*totM*ene]=data;}//i to j (ene)
    void input_vol(int i,real data){vol[i]=data;}
    void input_firstoutput(real data){firstoutput = data;}
    void input_specific_heat(int i,real data){specific_heat[i] = data;}
    void input_rho(int i,real data){rho[i] = data;}
    void input_heattrans(int i,real data){heattrans[i] = data;}
    void input_steps(real time,real tstep);
    void input_beta(real data);
    void input_vel(int i,real data){vel[i] = data;}
    void non_heatcal();
    void cal_Pij(int energy);
    void setup_Pij();
    void c_pijcs();
    void setup_dp();//遅発中性子先行核セット
    //以下定常用
    void test();
    //ベタ移植が下になります。
    double get_totalfiss();
    void cal_scat();
    void cssuesankaku();
    void setup_start();
    void solver();
    void setup_nxPhi();
    void cal_fiss();
    void plot_phi();
    void make_mmpa_a();
    void cal_mmpa_a();//A*nを行うだけのモジュール
    void make_mmpa_b();
    void cal_mmpa();
    void mmpa_clean();
    void cal_xs_f();
    void make_mmpa_a_n();//以下、複数種による分岐比考慮版。全体への適応？待たれよ。
    void cal_mmpa_n();
    void fact_sigt_mod(real factor);
    void fact_sig_mod(real factor);
    void Adjust_bq(real factor);
    void mmpa_switch(bool mmpa);

};
#endif
