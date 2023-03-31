#include<iostream>
#include "GeneralOption.h"

using namespace std;

GeneralOption::GeneralOption()
{
  forward=true;
  omegao=1.0;
  epsf=1e-4;
  epss=1e-4;
  epsk=1e-5;
  print=false;
  //outitermax=999;
  outitermax=49;
  cout<<"# Maximum OUTER ITERATION times is : "<<outitermax+1<<"\n";
  itcmfd=2;
  itmin_cmfd=0;
}

void GeneralOption::PutGrp(int i)
{
  grp=i;
  omegai.resize(i,1.0);
}

bool GeneralOption::Converged(real errf,real errk,real errs)
{
  if(errf<epsf&&errk<epsk&&errs<epss)return true;
  return false;
};
