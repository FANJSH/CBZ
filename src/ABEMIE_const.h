#ifndef CONST
#define CONST

const int maxreg=300;
const int maxmed=10;
const int maxib=900;
//const real htp=0.7746;
const real htp=0.85; // quadratic
const real htp2=0.65;  // linear
const real htp4=0.85; // cubic (change from 0.88 in 2010/7/5)
//const real htp4=0.88; // cubic
const real htp_inv=1./htp;
const real htp2_inv=1./htp2;
const real htp4_inv=1./htp4;
const real ht2=0.3;   // cubic
const real ht2_inv=1./ht2;

const real pp16[]={0.95012509387e-1, 0.28160355077,
                0.45801677766,    0.61787624440,
                0.75540440836,    0.86563120239,
                0.94457502307,    0.98940093499,
               -0.95012509387e-1,-0.28160355077,
               -0.45801677766,   -0.61787624440,
               -0.75540440836,   -0.86563120239,
               -0.94457502307,   -0.98940093499};
const real ww16[]={ 0.1894506105,   0.18260341504,
                 0.1691565194,   0.14959598882,
                 0.1246289713,   0.95158511682e-1,
                 0.6225352394e-1,0.27152459412e-1,
                 0.1894506105,   0.18260341504,
                 0.1691565194,   0.14959598882,
                 0.1246289713,   0.95158511682e-1,
                 0.6225352394e-1,0.27152459412e-1};
const real pp24[]={ 0.6405689286e-1, 0.19111886747,
                 0.3150426797   , 0.43379350763,
                 0.5454214714   , 0.64809365194,
                 0.7401241916   , 0.82000198597,
                 0.8864155270   , 0.93827455200,
                 0.9747285560   , 0.99518722000,
                -0.6405689286e-1,-0.19111886747,
                -0.3150426797   ,-0.43379350763,
                -0.5454214714   ,-0.64809365194,
                -0.7401241916   ,-0.82000198597,
                -0.8864155270   ,-0.93827455200,
                -0.9747285560   ,-0.99518722000};
const real ww24[]={ 0.1279381953   , 0.12583745635 ,
                 0.1216704729   , 0.11550566805 , 
                 0.1074442701   , 0.97618652104E-1 ,
                 0.8619016153E-1, 0.73346481411E-1 , 
                 0.5929858492E-1, 0.44277438817E-1 ,
                 0.2853138863E-1, 0.12341229800E-1 , 
                 0.1279381953   , 0.12583745635   ,
                 0.1216704729   , 0.11550566805   , 
                 0.1074442701   , 0.97618652104E-1 ,
                 0.8619016153E-1, 0.73346481411E-1 , 
                 0.5929858492E-1, 0.44277438817E-1 ,
                 0.2853138863E-1, 0.12341229800E-1};

#endif
