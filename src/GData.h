#ifndef GDATA
#define GDATA

#include <cstdlib>
#include <cmath> // add by hazama
#include <iostream>
#include "Numeric.h"
#include <vector>
#include <string>

using namespace std;

class GData{
 protected:
  int num_y;
  vector<real> x;
  vector< vector<real> > y;
  string tag_x;
  vector<string> tag_y;
  bool point;
 public:
  GData(bool point_inp=true){point=true;};
  GData(int num_y_inp, bool point_inp=true);
  ~GData(){};
  void PutNumY(int inp);
  int GetNumY(){return num_y;};
  int GetDataNum(){return x.size();};
  string GetTagY(int i){return tag_y[i];};
  vector<real> GetX(){return x;};
  vector<real> GetY(int i);
  void PutTagX(string tag_inp){tag_x=tag_inp;};
  void PutTagY(int i, string tag_inp);
  void push_back_x(real xinp){x.push_back(xinp);};
  void push_back_y(int i, real yinp);
  void show_self();
  void WriteFile(string mdir, string name);
  void ReadFile(string mdir, string name);
  void SlideX(real inp);
  void FactorizeX(real inp);
  real GetY_Lin(real x_inp, int id_y=0);
  real GetY_Log(real x_inp, int id_y=0);  
  void DeleteData(real x_val, real thresh=1e-8);
  void TakeDifferenceFrom(GData &ref, bool absolute=true);
  void TakeAbsoluteDifferenceFrom(GData &ref){TakeDifferenceFrom(ref, true);};
  void TakeRelativeDifferenceFrom(GData &ref){TakeDifferenceFrom(ref, false);};  
  vector<real> GetAbsoluteMaximumValue();
};

#endif


