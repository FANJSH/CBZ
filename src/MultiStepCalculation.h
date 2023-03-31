#ifndef MULTISTEPCALCULATION
#define MULTISTEPCALCULATION

#include <cstdlib>
#include <cmath> // add by hazama
#include <iostream>
#include "Numeric.h"
#include "GroupData.h"
#include "MultiStepCalculation.h"
#include <vector>
#include <string>

using namespace std;

class MultiStepCalculation{
 private:
 public:
  MultiStepCalculation(){};
  ~MultiStepCalculation(){};

  void MultiStepCalc48th6stp(GroupData2D &mat,GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);
  void MultiStepCalc48th43stp(GroupData2D &mat,GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);
  void MultiStepCalc48th275stp(GroupData2D &mat,GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);
  void MultiStepCalc48th1725stp(GroupData2D &mat,GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);
  void MultiStepCalc48th10787stp(GroupData2D &mat,GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep);

};

#endif


