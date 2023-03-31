#ifndef HIGHERMODECALCULATOR
#define HIGHERMODECALCULATOR

#include "Numeric.h"
#include "DHEX_system.h"
#include "PLOS_system.h"

#include<string>
#include<iostream>
#include<fstream>
#include<vector>

using namespace std;

class HigherModeCalculator{
 private:
  int iter_max;
  real epsk;
  int xsss,ysss,zsss;
  bool accel;
 public:
  HigherModeCalculator();
  void PutIterMax(int i){iter_max=i;};
  void PutEpsk(real i){epsk=i;};
  void PutInitialSourcePosition(int x,int y,int z);
  void RunDHEX(vector<DHEXSystem> &sysf, vector<DHEXSystem> &sysa,
               vector<real> &kf, int max_order, bool cor_tri=false);
  void RunPLOS(vector<PLOSSystem> &sysf, vector<PLOSSystem> &sysa,
               vector<real> &kf, int max_order);
  void PutNoAcceleration(){accel=false;};

};

#endif
