#ifndef BEMMEDIUM
#define BEMMEDIUM

#include "GroupData.h"
#include "Numeric.h"
#include "Medium.h"

class BemMedium:public Medium{
  GroupData1D bmat,bmatdk;
  GroupData2D amat,cmat,cmatdk,dcmat,cmati;
 public:
  BemMedium():Medium(){};
  BemMedium(int i):Medium(i){PutImaxBem(i);};
  void PutImaxBem(int i);
  void PutMedium(Medium inp);
  GroupData2D cal_amat(real keff);
  GroupData2D cal_amat_SP3_1grp(real keff);
  GroupData1D cal_bmat(GroupData2D am);
  GroupData2D cal_cmat(GroupData2D am,GroupData1D bm);
  void put_abc(real keff,real dk=0,bool sp3=false);
  GroupData1D& GetBmat(){return bmat;};
  GroupData1D& GetBmatdk(){return bmatdk;};
  GroupData2D& GetCmat(){return cmat;};
  GroupData2D& GetCmati(){return cmati;};
  GroupData2D& GetCmatdk(){return cmatdk;};
  GroupData2D& GetDcmat(){return dcmat;};
};

#endif
