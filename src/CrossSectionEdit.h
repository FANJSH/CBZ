#ifndef CROSSSECTIONEDIT
#define CROSSSECTIONEDIT

#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "Numeric.h"
#include "GeneralSystem.h"
#include "Medium.h"

using namespace std;

class CrossSectionEdit{
 private:
 public:
  CrossSectionEdit(){};
  Medium Homogenize(GeneralSystem &sys,int hreg, int*nreg,int cur_id=1); // cur_id:flux moment treated as current
  Medium HomogenizeAll(GeneralSystem &sys,int cur_id=1);  
  Medium HomogenizeCartesian2D(GeneralSystem &sys,int x1,int x2,int y1,int y2);
};

#endif


