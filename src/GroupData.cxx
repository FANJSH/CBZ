#include <cstdlib>
#include "GroupData.h"
#include "MultiStepCalculation.h"
#include<fstream>

using namespace std;

// GroupData //

real GroupData::get_sum()
{
  real ret=0.0;
  for(int i=0;i<imax;i++){
    ret+=dat[i];
  }
  return ret;
}

GroupData::~GroupData()
{
  //dat.clear();
}

real GroupData::GetEuclideanNorm()
{
  real sum=0.;
  for(int i=0;i<imax;i++){
    sum+=pow(dat[i],2);
  };
  sum=sqrt(sum);
  return sum;
};

void GroupData::DataClear()
{
  //vector<real> ().swap(dat);
  dat.clear();
};

// GroupData1D //

void GroupData1D::show_self()
{
  cout << " Group :" << imax << "\n";
  cout.setf(ios::scientific);
  cout.precision(6);
  for(int i=0;i<imax;i++){
    cout << i << "  " << dat[i] <<"\n";
  };
  cout.unsetf(ios::scientific);
}

GroupData1D GroupData1D::operator+(const GroupData1D &second)
{
  CheckSize(second);

  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,dat[i]+second.dat[i]);
  }
  return ret;
}

GroupData1D GroupData1D::operator-(const GroupData1D &second)
{
  CheckSize(second);

  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,dat[i]-second.dat[i]);
  }
  return ret;
}

GroupData1D GroupData1D::operator*(real a)
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,dat[i]*a);
  }
  return ret;
}

GroupData1D GroupData1D::mult(GroupData1D sec)
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,dat[i]*sec.get_dat(i));
  }
  return ret;
}

void GroupData1D::Multiply(GroupData1D &sec)
{
  for(int i=0;i<imax;i++){
    dat[i]*=sec.get_dat(i);
  };
}

void GroupData1D::Multiply(real a)
{
  for(int i=0;i<imax;i++){
    dat[i]*=a;
  };
}

real GroupData1D::operator*(const GroupData1D &second)
{
  real ret=0;
  for(int i=0;i<imax;i++){
    ret+=dat[i]*second.dat[i];
  };
  return ret;
}


GroupData1D GroupData1D::operator/(real a)
{
  GroupData1D ret(imax);
  real a_inv=1./a;
  for(int i=0;i<imax;i++){
    ret.put_data(i,dat[i]*a_inv);
  }
  return ret;
}

GroupData1D GroupData1D::operator/(const GroupData1D &second)
{
  CheckSize(second);

  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    if(second.dat[i]!=0.0){
      ret.put_data(i,dat[i]/second.dat[i]);}
    else{
      ret.put_data(i,0.0);};
  }
  return ret;
}

void GroupData1D::put_data(real *inp)
{
  for(int i=0;i<imax;i++){
    dat[i]=inp[i];
  };
}

void GroupData1D::put_data(vector<real> &inp)
{
  for(int i=0;i<imax;i++){
    dat[i]=inp[i];
  };
}

void GroupData1D::put_data(int i1,GroupData1D inp,int j,int k)
{
  if(k==0)k=inp.get_imax();
  for(int i=0;i<k;i++){
    dat[i1+i]=inp.get_dat(j+i);
  };
}

void GroupData1D::copy(GroupData1D s)
{
  put_imax(s.get_imax());
  for(int i=0;i<imax;i++){
    dat[i]=s.get_dat(i);
  };
}

GroupData1D GroupData1D::sqrt1d()
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    real a=dat[i];
    if(a>=0.0){
      ret.put_data(i,sqrt(a));
    }else{
      ret.put_data(i,sqrt(-a));
    };
  };
  return ret;
}

GroupData1D GroupData1D::get_dat(int i,int j)
{
  GroupData1D ret(j);
  for(int k=0;k<j;k++){
    ret.put_data(k,dat[i+k]);
  };
  return ret;
}

void GroupData1D::add_data(GroupData1D &sec)
{
  for(int i=0;i<imax;i++){
    dat[i]+=sec.get_dat(i);
  };
};

GroupData1D GroupData1D::GetInverse()
{
  GroupData1D ret(imax);
  for(int i=0;i<imax;i++){
    ret.put_data(i,1./get_dat(i));
  };
  return ret;
};

void GroupData1D::CheckSize(GroupData1D sec)
{
  if(imax!=sec.get_imax()){
    cout<<"There is inconsistency between two instances of GroupData1D.\n";
  };
} 

real GroupData1D::Cond(GroupData1D &w)
{
  int ngrp=1;
  vector<int> bgrp(1);
  bgrp[0]=imax-1;
  return Cond(w,ngrp,bgrp).get_dat(0);
};

GroupData1D GroupData1D::Cond(GroupData1D &w,int ngrp,vector<int> bgrp)
{
  if(imax!=w.get_imax()){
    cout<<"# Error in GroupData1D::Cond.\n";
    cout<<"# There is inconsistency in group\n";
  };

  GroupData1D ret(ngrp);

  int is=0;
  for(int i=0;i<ngrp;i++){
    real sum1=0;
    real sum2=0;
    for(int j=is;j<=bgrp[i];j++){
      sum1+=dat[j]*w.get_dat(j);
      sum2+=w.get_dat(j);
    };
    is=bgrp[i]+1;
    if(sum2!=0.){
      ret.put_data(i,sum1/sum2);
    }else{
      ret.put_data(i,0.);
    };
  };
  return ret;
};

GroupData1D GroupData1D::Cond(GroupData1D &w,int ngrp,int* bgrp)
{
  vector<int> bgrp2(ngrp);
  for(int i=0;i<ngrp;i++){
    bgrp2[i]=bgrp[i];
  };
  return Cond(w,ngrp,bgrp2);
};

GroupData1D GroupData1D::CondSum(int ngrp,vector<int> bgrp)
{
  GroupData1D ret(ngrp);
  ret.set_zero();

  int is=0;
  for(int i=0;i<ngrp;i++){
    real sum=0;
    for(int j=is;j<=bgrp[i];j++){
      sum+=dat[j];
    };
    is=bgrp[i]+1;
    ret.put_data(i,sum);
  };
  return ret;
};

GroupData1D GroupData1D::CondSum(int ngrp,int *bgrp)
{
  vector<int> bgrp2(ngrp);
  for(int i=0;i<ngrp;i++){
    bgrp2[i]=bgrp[i];
  };
  return CondSum(ngrp,bgrp2);
}

GroupData2D GroupData2D::GetPartialCorrelationMatrix()
{
  if(x!=y){
    cout<<"# Error in GroupData2D::GetPartialCorrelationMatrix.\n";
    cout<<"# Partial correlation matrix cannot be calculated from non-square matrix.\n";
    exit(0);
  };

  if(reduced_form){
    cout<<"# Error in GroupData2D::GetPartialCorrelationMatrix.\n";
    cout<<"# Reduced-form instance cannot be treated.\n";
    exit(0);
  };
  
  GroupData2D ret(x,y);
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      ret.put_data(i,j,dat[index]);
      index++;
    };
  };

  GroupData2D ret_inv=ret.inverse();
  for(int i=0;i<y;i++){
    if(ret_inv.get_dat(i,i)<0.){
      cout<<"# Error in GroupData2D::GetPartialCorrelationMatrix.\n";      
      cout<<"# Diagonal element in the inverse of correlation matrix becomes negative.\n";
      exit(0);
    };
  };
  
  for(int i=0;i<y;i++){
    real tmp1=sqrt(ret_inv.get_dat(i,i));
    for(int j=0;j<x;j++){
      real tmp2=sqrt(ret_inv.get_dat(j,j));      
      ret.put_data(i,j,-ret_inv.get_dat(i,j)/tmp1/tmp2);
    };
  };
  
  return ret;
};

/////////////////
// GroupData2D //
/////////////////

real GroupData2D::get_dat(int i)
{
  if(reduced_form)NormalForm();
  return GroupData::get_dat(i);
};

real GroupData2D::get_dat(int i,int j)const
{
  if(reduced_form){
    int sz=xpos[i].size();
    for(int k=0;k<sz;k++){
      if(xpos[i][k]==j)return val[i][k];
    };
    return 0.;
  }else{
    return dat[i*x+j];
  };
};

real GroupData2D::get_sum()
{
  if(!reduced_form){
    return GroupData::get_sum();
  }else{
    real sum=0.;
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	sum+=val[i][j];
      };
    };
    return sum;
  };
};

real GroupData2D::GetEuclideanNorm()
{
  if(!reduced_form){
    return GroupData::GetEuclideanNorm();
  }else{
    real sum=0.;
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	sum+=pow(val[i][j],2);
      };
    };
    return sqrt(sum);
  };
};

void GroupData2D::put_data(real *inp)
{
  if(reduced_form){
    cout<<"# Error in GroupData2D::put_data.\n";
    cout<<"# This instance is in a reduced form.\n";
    exit(0);
  };

  for(int i=0;i<x*y;i++){
    dat[i]=inp[i];
  }
}

void GroupData2D::put_data(int x1,int x2,int y1,int y2,real inp)
{
  if(reduced_form){
    cout<<"# Error in GroupData2D::put_data.\n";
    cout<<"# This instance is in a reduced form.\n";
    exit(0);
  };

  if(x1<0||x1>=x||x2<0||x2>=x||y1<0||y1>=y||y2<0||y2>=y){
    cout<<"# Error in put_data in GroupData2D.\n";
    cout<<"# It is impossible to put the value.\n";
    exit(0);
  };
  for(int j=y1;j<=y2;j++){
    for(int i=x1;i<=x2;i++){
      put_data(i,j,inp);
    };
  };
};

void GroupData2D::add_data(int i,int j,real f)
{
  if(reduced_form){
    int sz=xpos[i].size();
    for(int k=0;k<sz;k++){
      if(xpos[i][k]==j){
	val[i][k]+=f;
	return;
      };
    };
    PutXpos(i,j);
    PutVal(i,f);
  }else{
    dat[i*x+j]+=f;
  };
};

void GroupData2D::put_yx(int i,int j)
{
  /*
  if(reduced_form){
    cout<<"# Error in GroupData2D::put_yx.\n";
    cout<<"# This instance is in a reduced form.\n";
    exit(0);
  };
  */
  y=i; x=j; 
  if(reduced_form){
    xpos.resize(y);
    val.resize(y);
  }else{
    put_imax(i*j);
  };
}

void GroupData2D::show_self(bool fine)
{
  bool trans=false;
  if(reduced_form){
    trans=true;
    NormalForm();
  };

  cout <<"x:"<<x<<" & y:"<<y<<"\n";
  if(fine){
    for (int i=0;i<y;i++){
      cout<<"*"<<i<<"*";
      for(int j=0;j<x;j++){
        cout.setf(ios::scientific);
        cout.precision(6);
        cout << dat[i*x+j] <<":";
        cout.unsetf(ios::scientific);
      };
      cout << "\n";
    };
  }else{
    for (int i=0;i<y;i++){
      cout<<"*"<<i<<"*";
      for(int j=0;j<x;j++){
        cout.setf(ios::showpoint);
        cout.precision(2);
        cout << dat[i*x+j] <<":";
        cout.precision(6);
        cout.unsetf(ios::showpoint);
      };
      cout << "\n";
    };
    /*
    for (int i=0;i<y;i++){
      for(int j=0;j<x;j++){
	if(dat[i*x+j]==0.){
	  cout<<i<<" "<<j<<" 0\n";
	}else{
	  cout<<i<<" "<<j<<" 1\n";
	};
      };
      cout << "\n";
    };
    */
  };

  if(trans)ReducedForm();
}

GroupData1D GroupData2D::get_diag()
{
  CheckSquare();
  GroupData1D ret(y);

  if(reduced_form){
    ret.set_zero();
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	if(i==xpos[i][j]){
	  ret.put_data(i,val[i][j]);
	};
      };
    };
  }else{
    for(int i=0;i<y;i++){
      ret.put_data(i,dat[i*x+i]);
    }
  };

  return ret;
}

GroupData1D GroupData2D::get_sumx()
{
  GroupData1D ret(y);
  
  if(reduced_form){
    for(int i=0;i<y;i++){
      real sum=0.;
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	sum+=val[i][j];
      };
      ret.put_data(i,sum);
    };
  }else{
    for(int i=0;i<y;i++){
      real sum=0.;
      int ix=i*x;
      for(int j=0;j<x;j++){
        sum+=dat[ix+j];
      };
      ret.put_data(i,sum);
    };
  };

  return ret;
}

GroupData2D GroupData2D::Slicing(int y0, int y1, int x0, int x1)
{
  bool err_flag=false;
  if(y0<0||y0>=y)err_flag=true;
  if(y1<0||y1>=y)err_flag=true;
  if(x0<0||x0>=x)err_flag=true;
  if(x1<0||x1>=x)err_flag=true;
  if(y0>y1)err_flag=true;
  if(x0>x1)err_flag=true;

  if(err_flag){
    cout<<"# Error in GroupData2D::Slicing.\n";
    exit(0);
  };

  GroupData2D ret(y1-y0+1, x1-x0+1);
  for(int i=0;i<y1-y0+1;i++){
    for(int j=0;j<x1-x0+1;j++){
      ret.put_data(i,j,get_dat(y0+i,x0+j));
    };
  };

  return ret;
};

GroupData1D GroupData2D::Slicing(int y0, int y1, int x0)
{
  bool err_flag=false;
  if(y0<0||y0>=y)err_flag=true;
  if(y1<0||y1>=y)err_flag=true;
  if(x0<0||x0>=x)err_flag=true;
  if(y0>y1)err_flag=true;

  if(err_flag){
    cout<<"# Error in GroupData2D::Slicing.\n";
    exit(0);
  };

  GroupData1D ret(y1-y0+1);
  for(int i=0;i<y1-y0+1;i++){
    ret.put_data(i,get_dat(y0+i,x0));
  };

  return ret;
};

void GroupData2D::put_unit()
{
  if(reduced_form){
    cout<<"# Error in GroupData2D::put_unit.\n";
    cout<<"# This instance is in a reduced form.\n";
  };

  CheckSquare();
  set_zero();
  for(int i=0;i<y;i++){
    dat[i*x+i]=1.0;
  };
  //ReducedForm();
}

GroupData2D GroupData2D::operator+(const GroupData2D &second)
{
  CheckSizeXY(second);
  GroupData2D ret(y,x);

  if(reduced_form){
    ret.set_zero();
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	ret.put_data(i,xpos[i][j],val[i][j]);
      };
    };
  }else{
    for(int i=0;i<y*x;i++){
      ret.put_data(i,dat[i]);
    };
  };

  if(second.IsReducedForm()){
    for(int i=0;i<y;i++){
      int sz=second.GetXposSize(i);
      for(int j=0;j<sz;j++){
	ret.add_data(i,second.GetXpos(i,j),second.GetVal(i,j));
      };
    };
  }else{
    for (int i=0;i<y;i++){
      for(int j=0;j<x;j++){
        ret.add_data(i,j,second.get_dat(i,j));
      };
    };
  };

  return ret;
}

GroupData2D GroupData2D::GetTransposedMatrix()
{

  if(reduced_form){
    GroupData2D ret;
    ret.ReducedForm();
    ret.put_yx(x,y);
    for(int i=0;i<x;i++){
      for(int j=0;j<y;j++){
	int sz=xpos[j].size();
	for(int k=0;k<sz;k++){
	  if(xpos[j][k]==i){
	    ret.PutXpos(i,j);
	    ret.PutVal(i,val[j][k]);
	  };
	};
      };
    };
    return ret;
  }else{
    GroupData2D ret(x,y);
    for(int i=0;i<ret.get_y();i++){
      for(int j=0;j<ret.get_x();j++){
        ret.put_data(i,j,get_dat(j,i));
      };
    };
    return ret;
  };
};

void GroupData2D::Transposition()
{
  if(reduced_form){
    GroupData2D newmat=GetTransposedMatrix();
    this->copy(newmat);
  }else{
  for(int i=0;i<y;i++){
    int id=i*x+(i+1);
    for(int j=i+1;j<x;j++){
      real org=dat[id];
      dat[id]=dat[j*x+i];
      dat[j*x+i]=org;
      id++;
    };
  };
  };
};

GroupData2D GroupData2D::operator-(const GroupData2D &second)
{
  CheckSizeXY(second);
  GroupData2D ret(y,x);

  if(reduced_form){
    ret.set_zero();
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	ret.put_data(i,xpos[i][j],val[i][j]);
      };
    };
  }else{
    for(int i=0;i<y*x;i++){
      ret.put_data(i,dat[i]);
    };
  };

  if(second.IsReducedForm()){
    for(int i=0;i<y;i++){
      int sz=second.GetXposSize(i);
      for(int j=0;j<sz;j++){
	ret.add_data(i,second.GetXpos(i,j),-second.GetVal(i,j));
      };
    };
  }else{
    for (int i=0;i<y;i++){
      for(int j=0;j<x;j++){
        ret.add_data(i,j,-second.get_dat(i,j));
      };
    };
  };

  return ret;
}

GroupData2D GroupData2D::operator*(const GroupData2D &second)
{
  if(x!=second.y){
    cout<<"# Error in GroupData2D::*(G2D)\n";
    cout<<"# There is inconsistency\n";
    cout<<"   "<<x<<" & "<<second.y<<"\n";
    exit(0);
  };

  if(reduced_form&&second.IsReducedForm()){
    cout<<"# Error in GroupData2D::operator*.\n";
    cout<<"# Not yet coded for reduced-formed original and secondary matrices.\n";
    exit(0);
  };

  GroupData2D ret(y,second.x);
  real *ret2=new real[y*second.x];
  for(int i=0;i<y*second.x;i++){ret2[i]=0.;};

  if(reduced_form){
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      real tmp=0.;
      for(int j=0;j<sz;j++){
	real vv=val[i][j];
	int pos=xpos[i][j];
	int tmp2=pos*second.x;
	for(int k=0;k<second.x;k++){
	  ret2[i*second.x+k]+=vv*second.dat[tmp2++];
	};
      };
    };
  }else{
    if(second.IsReducedForm()){
      for(int i=0;i<second.y;i++){
	int sz=second.GetXposSize(i);
	for(int j=0;j<sz;j++){
	  int pos=second.GetXpos(i,j);
	  real vv=second.GetVal(i,j);
	  for(int k=0;k<y;k++){
	    ret2[k*second.x+pos]+=vv*dat[k*x+i];
	  };
	};
      };
    }else{
      for(int i=0;i<y;i++){
        int tmp1=i*x;
        for(int j=0;j<x;j++){
          real tmp=dat[tmp1++];
          if(tmp!=0.){
            int tmp2=j*second.x;
            for(int k=0;k<second.x;k++){
  	      ret2[i*second.x+k]+=tmp*second.dat[tmp2++];
            };
          };
        };
      };
    };

  };

  ret.put_data(ret2);
  delete [] ret2;

  return ret;
}

/*
GroupData2D GroupData2D::SquaringForSparceMatrix()
{
  if(x!=y){
    cout<<"# Error in GroupData2D::SquaringForSparceMatrix\n";
    cout<<"# This method works only for square matrix.\n";
    exit(0);
  };

  GroupData2D ret(x,x);
  ret.set_zero();

  vector<int> num_non_zero(y);
  vector< vector<int> > pos_non_zero(y);
  vector< vector<real> > val_non_zero(y);

  int index=0;
  for(int i=0;i<y;i++){
    int cnt=0;
    for(int j=0;j<x;j++){
      real tmp=dat[index++];
      if(tmp!=0.){
	cnt++;
	pos_non_zero[i].push_back(j);
	val_non_zero[i].push_back(tmp);
      };
    };
    num_non_zero[i]=cnt;
  };

  for(int i=0;i<y;i++){
    for(int j=0;j<num_non_zero[i];j++){
      real tmp=val_non_zero[i][j];
      int pos2=pos_non_zero[i][j];
      for(int k=0;k<num_non_zero[pos2];k++){
        int tmp2=pos_non_zero[pos2][k];
        ret.add_data(i,tmp2,tmp*val_non_zero[pos2][k]);
      };
    };
  };

  return ret;
}
*/

GroupData1D GroupData2D::operator*(const GroupData1D &g1d)
{
  CheckSizeX(g1d);
  GroupData1D ret(y);

  if(reduced_form){
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      real tmp=0.;
      for(int j=0;j<sz;j++){
        tmp+=val[i][j]*g1d.get_dat(xpos[i][j]);
      };
      ret.put_data(i,tmp);
    };
    return ret;
  };

  ret.set_zero();
  for (int i=0;i<y;i++){
    int pos=i*x;
    for (int j=0;j<x;j++){
      real tmp=dat[pos++];
      if(tmp!=0.)ret.add_data(i,tmp*g1d.get_dat(j));
      //if(tmp!=0.)ret.add_data(i,tmp*g1d.dat[j]);
    };
  };

  return ret;
}

/*
GroupData2D GroupData2D::operator/(GroupData1D g1d)
{
  CheckSizeY(g1d);

  GroupData2D ret(y,x);
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      ret.put_data(i,j,dat[i*x+j]/g1d.get_dat(i));
    };
  };
  return ret;
}
*/

GroupData2D GroupData2D::operator*(real a)
{

  if(reduced_form){
    GroupData2D ret;
    ret.copy(*this);
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      for(int j=0;j<sz;j++){
	ret.PutVal(i,j,val[i][j]*a);
      };
    };
    return ret;
  }else{
    GroupData2D ret(y,x);
    for(int i=0;i<y;i++){
      int tmp=i*x;
      for(int j=0;j<x;j++){
        real val=dat[tmp++];
        if(val!=0.)ret.put_data(i,j,val*a);
      };
    };
    return ret;
  };
}

void GroupData2D::Multiply(GroupData1D &sec)
{
  if(reduced_form){
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      real tmp=sec.get_dat(i);
      for(int j=0;j<sz;j++){
	val[i][j]*=tmp;
      };
    };
  }else{
  int id=0;
  for(int i=0;i<y;i++){
    real tmp=sec.get_dat(i);
    for(int j=0;j<x;j++){
      dat[id]*=tmp;
      id++;
    };
  };
  };
}

void GroupData2D::PasteGroupData2D(int yp, int xp, GroupData2D &inp)
{
  if(reduced_form){
    cout<<"# Error in GroupData2D::PasteGroupData2D.\n";
    cout<<"# Not yet coded for Reduced-formed instance.\n";
    exit(0);
  };

  int inp_y=inp.get_y();
  int inp_x=inp.get_x();

  if(yp<0||xp<0){
    cout<<"# Error in GroupData2D::PasteGroupData2D.\n";
    cout<<"# Positions to be pasted should be positive integer.\n";
    exit(0);
  };

  if(yp+inp_y>y||xp+inp_x>x){
    cout<<"# Error in GroupData2D::PasteGroupData2D.\n";
    cout<<"# The size of original GroupData2D instance is inappropriate to paste.\n";
    exit(0);
  };

  int index=0;
  for(int yy=0;yy<inp_y;yy++){
    for(int xx=0;xx<inp_x;xx++){
      real tmp=inp.get_dat(index);
      index++;
      dat[(yp+yy)*x+xp+xx]=tmp;
    };
  };

};

void GroupData2D::PasteGroupData1D(int yp, int xp, GroupData1D &inp)
{
  if(reduced_form){
    cout<<"# Error in GroupData2D::PasteGroupData1D.\n";
    cout<<"# Not yet coded for Reduced-formed instance.\n";
    exit(0);
  };

  int inp_y=inp.get_imax();

  if(yp<0||xp<0){
    cout<<"# Error in GroupData2D::PasteGroupData1D.\n";
    cout<<"# Positions to be pasted should be positive integer.\n";
    exit(0);
  };

  if(yp+inp_y>y){
    cout<<"# Error in GroupData2D::PasteGroupData1D.\n";
    cout<<"# The size of original GroupData2D instance is inappropriate to paste.\n";
    exit(0);
  };

  int index=0;
  for(int yy=0;yy<inp_y;yy++){
    real tmp=inp.get_dat(index);
    index++;
    dat[(yp+yy)*x+xp]=tmp;
  };

};

GroupData2D GroupData2D::inverseKatagiri()
{
  int N=y;
  GroupData2D ret(N,N);

  vector< vector<real> > B1(N);
  vector< vector<real> > B2(N);
  for(int i=0;i<N;i++){
    B1[i].resize(N);
    B2[i].resize(N);
    for(int j=0;j<N;j++){
      B1[i][j]=dat[i*N+j];
    };
  };

  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      if(i==j){
	B2[i][j]=1.;
      }
    }
  }
  double w;//記憶用変数
  for(int i=0;i<N;i++){
    w=1/B1[i][i];
    for(int j=0;j<N;j++){
      B1[i][j]*=w;
      B2[i][j]*=w;
    }
    for(int j=0;j<N;j++){
      if(i!=j){
	w=B1[j][i];
	for(int k=0;k<N;k++){
	  B1[j][k]-=B1[i][k]*w;
	  B2[j][k]-=B2[i][k]*w;
	}
      }
    }
  }


  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret.put_data(i,j,B2[i][j]);
    };
  };

  return ret;

};

GroupData2D GroupData2D::inverse()
{
  if(reduced_form){
    cout<<"# Not coded in G2D inverse.\n";
    exit(0);
  };

  CheckSquare();

  GroupData2D amat(x,x);
  int index=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      amat.put_data(index,dat[index]);
      index++;
    };
  };

  GroupData2D bmat(x,x);

  // (LU decomposition)

  amat.LUdecomposition();
  GroupData1D tmp(x);

  for(int i=0;i<x;i++){
    tmp.set_zero();
    tmp.put_data(i,1.);
    amat.solveaxb_LUdecomposition(tmp);
    bmat.PasteGroupData1D(0,i,tmp);
  };

  // (direct calculation)
  /*
  bmat.put_unit();
  amat.solveaxb_mod(bmat);
  */

  return bmat;
}

void GroupData2D::DoInverse(bool pivot_move)
{
  CheckSquare();

  if(pivot_move){
    gauss_inverse(dat,y); // NOT-CODED
  }else{
    gauss_inverse(dat,y);
    //lu_decomposition_inverse(dat,y);
  };

}

void GroupData2D::LUdecomposition()
{
  if(y!=x){
    cout<<"# Error in GroupData2D::LUdecomposition.\n";
    cout<<"# This matrix is not square.\n";
    exit(0);
  };
  lu_decomposition(dat,y);
};

void GroupData2D::solveaxb_LUdecomposition(GroupData1D &sec)
{
  int n=y;
  GroupData1D ret(y);

  if(reduced_form){

  // Forward substitution
  for(int i=0;i<n;i++){
    real tmp=sec.get_dat(i);
    int sz=xpos[i].size();
    for(int j=0;j<sz;j++){
      int pos=xpos[i][j];
      if(pos<i){
	tmp-=val[i][j]*ret.get_dat(pos);
      }else{
	j=sz;
      };
    };
    ret.put_data(i,tmp);
  };

  // Backward substitution
  for(int i=n-1;i>=0;i--){
    real tmp=ret.get_dat(i);
    int sz=xpos[i].size();
    for(int j=sz-1;j>=0;j--){
      int pos=xpos[i][j];
      if(pos>i){
	tmp-=val[i][j]*sec.get_dat(pos);
      }else if(pos==i){
	tmp/=val[i][j];
      }else{
	j=-1;
      };
    };	  
    sec.put_data(i,tmp);
  };

  return;

  };

  // Forward substitution
  for(int i=0;i<n;i++){
    real tmp=sec.get_dat(i);
    for(int j=0;j<i;j++){
      real tmp2=dat[i*n+j];
      if(tmp2!=0.){tmp-=tmp2*ret.get_dat(j);};
    };
    ret.put_data(i,tmp);
  };

  // Backward substitution
  for(int i=n-1;i>=0;i--){
    real tmp=ret.get_dat(i);
    for(int j=n-1;j>i;j--){
      real tmp2=dat[i*n+j];
      if(tmp2!=0.){tmp-=tmp2*sec.get_dat(j);};
    };
    sec.put_data(i,tmp/dat[i*n+i]);
  };

};

vector<int> GroupData2D::PivotCheck(real eps)
{
  vector<int> order;
  gauss_pivot_check(dat,order,y,eps);
  return order;
};

GroupData2D GroupData2D::solveaxb(GroupData2D &g2d)
{
  if(reduced_form){
    cout<<"# Not coded in G2D solveaxb.\n";
    exit(0);
  };

  CheckSquare();
  CheckSizeXY(g2d);

  GroupData2D ret(y,x);
  real *amat=new real[y*x];
  real *bmat=new real[y*x];

  for(int i=0;i<y;i++){
    int tmp=i*x;
    for(int j=0;j<x;j++){
      amat[tmp]=get_dat(tmp);
      bmat[tmp]=g2d.get_dat(tmp);
      tmp++;
    };
  };

  gauss(amat,bmat,y,y);

  ret.put_data(bmat);

  delete[] amat;
  delete[] bmat;
  return ret;
}

void GroupData2D::solveaxb_mod(GroupData2D &g2d)
{
  if(reduced_form){
    cout<<"# Not coded in G2D solveaxb?mod.\n";
    exit(0);
  };

  CheckSquare();
  CheckSizeXY(g2d);
  gauss(dat,g2d.get_dat(),y,y);
}

void GroupData2D::solveaxb_mod(GroupData1D &g1d)
{
  if(reduced_form){
    cout<<"# Not coded in G2D solveaxb?mod.\n";
    exit(0);
  };
  CheckSquare();
  CheckSizeY(g1d);
  gauss(dat,g1d.get_dat(),y,1);
}

GroupData1D GroupData2D::solveaxb(GroupData1D &g1d)
{
  if(reduced_form){
    cout<<"# Not coded in G2D solveaxb.\n";
    exit(0);
  };

  CheckSquare();
  CheckSizeY(g1d);

  real *amat=new real[y*x];
  real *bmat=new real[y];

  GroupData1D ret(y);

  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      amat[i*x+j]=get_dat(i*x+j);
    };
    bmat[i]=g1d.get_dat(i);
  };

  gauss(amat,bmat,y,1);

  ret.put_data(bmat);

  delete[] amat;
  delete[] bmat;
  return ret;
}


void GroupData2D::copy(GroupData2D s)
{
  if(s.IsReducedForm()){
    reduced_form=true;
    put_yx(s.get_y(),s.get_x());
    for(int i=0;i<y;i++){
      int sz=s.GetXposSize(i);
      for(int j=0;j<sz;j++){
	PutVal(i,s.GetVal(i,j));
	PutXpos(i,s.GetXpos(i,j));
      };
    };
  }else{
    put_yx(s.get_y(),s.get_x());
    for(int i=0;i<x*y;i++){
      dat[i]=s.get_dat(i);
    };
  };
}

void GroupData2D::CheckSizeX(GroupData1D sec)
{
  if(x!=sec.get_imax()){
    cout<<"# There is inconsistency between instances\n";
    cout<<"# of GroupData1D and GroupData2D\n";
    cout<<"# "<<x<<" & "<<sec.get_imax()<<"\n";
  };
}

void GroupData2D::CheckSizeY(GroupData1D sec)
{
  if(y!=sec.get_imax()){
    cout<<"# There is inconsistency between instances\n";
    cout<<"# of GroupData1D and GroupData2D\n";
    cout<<"# "<<y<<" & "<<sec.get_imax()<<"\n";
  };
}

void GroupData2D::CheckSizeXY(GroupData2D sec)
{
  if(x!=sec.get_x()||y!=sec.get_y()){
    cout<<"# There is inconsistency between instances\n";
    cout<<"# of GroupData2D and GroupData2D\n";
    cout<<"# "<<x<<" & "<<sec.get_imax()<<"\n";
    cout<<"# "<<y<<" & "<<sec.get_imax()<<"\n";
  };
}

void GroupData2D::CheckSquare()
{
  if(x!=y){
    cout<<"# The GroupData2D is not Square matrix\n";
    cout<<"#   (x,y) = "<<x<<" "<<y<<"\n";
  };
}


GroupData2D GroupData2D::Cond(GroupData1D &w,int ngrp,vector<int> bgrp)
{
  GroupData2D ret(ngrp,ngrp);

  bool trans=false;
  if(reduced_form){
    NormalForm();
    trans=true;
  };

  int is=0;
  for(int i=0;i<ngrp;i++){
    int is2=0;
    for(int i2=0;i2<ngrp;i2++){
      //if(i2>=i){
        real sum1=0;
        real sum2=0;
        for(int j=is;j<=bgrp[i];j++){
	  sum2+=w.get_dat(j);
          for(int j2=is2;j2<=bgrp[i2];j2++){
	    sum1+=dat[j*x+j2]*w.get_dat(j);
          };
        };
        ret.put_data(i,i2,sum1/sum2);
	//};
      is2=bgrp[i2]+1;
    };
    is=bgrp[i]+1;
  };

  if(trans)ReducedForm();

  return ret;
};

GroupData2D GroupData2D::Cond(GroupData1D &w,int ngrp,int *bgrp)
{
  vector<int> bgrp2(ngrp);
  for(int i=0;i<ngrp;i++){
    bgrp2[i]=bgrp[i];
  };
  return Cond(w,ngrp,bgrp2);
};

real GroupData2D::InfiniteNorm()
{
  real maxv=0.;

  if(reduced_form){
    for(int i=0;i<y;i++){
      int sz=xpos[i].size();
      real sum=0.;
      for(int j=0;j<sz;j++){
	sum+=fabs(val[i][j]);
      };
      if(sum>maxv)maxv=sum;
    };
  }else{
  int index=0;
  for(int i=0;i<y;i++){
    real sum=0.;
    for(int j=0;j<x;j++){
      sum+=fabs(dat[index++]);
    };
    if(sum>maxv)maxv=sum;
  };
  };

  return maxv;
};

void GroupData2D::Diagonalization(GroupData2D &diagmat,GroupData2D &evecmat)
{
  if(reduced_form)NormalForm();

  if(x!=y){
    cout<<"# Error in GroupData2D::Diagonalization.\n";
    cout<<"# Matrix should be square.\n";
    exit(0);
  };

  int n=x;
  real *matorg=new real[n*n];
  for(int i=0;i<n*n;i++){
    matorg[i]=dat[i];
  };

  real *mat2=new real[n*n];  
  MatrixDiagonalization(matorg,mat2,n);
  // C=U^{-1}PU
  // (input)  matorg:P 
  // (output) matorg:C, mat2:U

  diagmat.put_yx(n,n);
  diagmat.put_data(matorg);

  evecmat.put_yx(n,n);
  evecmat.put_data(mat2);

  delete [] matorg;
  delete [] mat2;
};

void GroupData2D::Diagonalization(GroupData2D &evecmat)
{
  if(reduced_form)NormalForm();

  if(x!=y){
    cout<<"# Error in GroupData2D::Diagonalization.\n";
    cout<<"# Matrix should be square.\n";
    exit(0);
  };

  int n=x;
  real *matorg=new real[n*n];
  for(int i=0;i<n*n;i++){
    matorg[i]=dat[i];
  };

  real *mat2=new real[n*n];  
  MatrixDiagonalization(matorg,mat2,n);
  // C=U^{-1}PU
  // (input)  matorg:P 
  // (output) matorg:C, mat2:U

  for(int i=0;i<n;i++){
    if(matorg[i*n+i]<0){
      matorg[i*n+i]=0.;
    }else{
      matorg[i*n+i]=sqrt(matorg[i*n+i]);
    };
  };

  evecmat.put_yx(n,n);

  int index=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      evecmat.put_data(i,j,mat2[index]*matorg[j*n+j]);
      index++;
    };
  };

  delete [] matorg;
  delete [] mat2;
};

void GroupData2D::CalEigenvaluesByLeverrieFaddeev(GroupData1D &ret)
{
  if(reduced_form)NormalForm();

  if(x!=y){
    cout<<"# Error in GroupData2D::CalEigenvaluesByLeverrieFaddeev.\n";
    cout<<"# Matrix should be square.\n";
    exit(0);
  };
 
  ret.put_imax(x);

  GroupData2D mat(x,x);
  GroupData2D unit(x,x);
  real *pval=new real[x];

  for(int i=0;i<x;i++){
    unit.put_unit();
    if(i==0){mat=(*this);}
    else{mat=(*this)*mat;}
    pval[i]=mat.get_diag().get_sum()/(i+1);
    unit=unit*pval[i];
    mat=mat-unit;
  };


  if(x==1){
    ret.put_data(0,pval[0]);
  }else if(x==2){
    real x=sqrt(pval[0]*pval[0]+4*pval[1]);
    ret.put_data(0,(pval[0]+x)*0.5);
    ret.put_data(1,(pval[0]-x)*0.5);
  }else if(x==3){
    cout<<"# (B2)^3 -"<<pval[0]<<"*(B2)^2 -"<<pval[1]<<"*(B2) -"<<pval[2]<<"\n";
    vector<real> sol;

    bool conv=false;
    real vinp=1e-9;
    while(vinp<100.&&!conv){
      real v0=vinp;
      for(int it=0;it<100;it++){
        real f=0.;
        for(int i=0;i<=x;i++){
          real coef=1.;
  	  if(i!=x)coef=-pval[x-1-i];
	  f+=pow(v0,i)*coef;
        };
        real df=0.;
        for(int i=0;i<x;i++){
	  real coef=1.;
	  if(i!=x-1)coef=-pval[x-1-i-1];
	  df+=pow(v0,i)*coef*(i+1);
        };
        real dv=f/df;
        v0-=dv;
        if(fabs(f)<1e-10){
	  conv=true;
	  it=100;
          sol.push_back(v0);
          cout<<"# First solution is detected ! "<<v0<<"\n";
          real b1=v0-pval[0];
          real b2=pval[2]/v0;
	  real discriminant=b1*b1-4*b2;
	  if(discriminant>0.){
	    //cout<<"# Other two solutions are real.\n";
	    sol.push_back(0.5*(-b1+sqrt(discriminant)));
	    sol.push_back(0.5*(-b1-sqrt(discriminant)));
	  }else{
	    cout<<"# Other two solutions are complex conjugate.\n";
	    cout<<"#   Real part: "<<0.5*(-b1)<<" / Imaginary part(+/-): "<<0.5*sqrt(-discriminant)<<"\n";
	  };
        };
      };
      if(vinp>0){
	vinp*=-1;
      }else{
	vinp*=-1;
	vinp+=1e-5;
      };
      //if(conv)vinp=100.;
    };

    ret.put_data(sol);
  }else{
    /*
    for(int x=-1000;x<1000;x++){
      real xp=x*0.001;
      cout<<xp<<" "<<pow(xp,3)-pval[0]*pow(xp,2)-pval[1]*xp-pval[2]<<"\n";
    };
    exit(0);
    */
    vector<real> sol;
    bool conv=false;

    real vinp=1e-9;
    while(vinp<100.){
      real v0=vinp;
      for(int it=0;it<100;it++){
        real f=0.;
        for(int i=0;i<=x;i++){
          real coef=1.;
  	  if(i!=x)coef=-pval[x-1-i];
	  f+=pow(v0,i)*coef;
        };
        real df=0.;
        for(int i=0;i<x;i++){
	  real coef=1.;
	  if(i!=x-1)coef=-pval[x-1-i-1];
	  df+=pow(v0,i)*coef*(i+1);
        };
        real dv=f/df;
        v0-=dv;
        if(fabs(f)<1e-10){
          bool eql=false;
          for(int i=0;i<int(sol.size());i++){
	    if(fabs(sol[i]/v0-1.)<1e-5)eql=true;
	  };
	  if(!eql){
            sol.push_back(v0);
	    cout<<"# Detected ! "<<v0<<"\n";
	    if(int(sol.size())==x)conv=true;
	  };
          break;
        };
      };
      if(vinp>0){
	vinp*=-1;
      }else{
	vinp*=-1;
	vinp+=1e-5;
      };
      if(conv)vinp=100.;
    };
    ret.put_data(sol);
  };


  delete [] pval;
};

// (Matrix exponential calculation)

GroupData2D GroupData2D::CalMatrixExponentialByPade(real factor)
{
  // This method calculates the matrix exponential for (self*factor)
  //
  // Numerical procedure is based on the paper on Expokit by Sidje
  // (ACM Transactions on Methematical Software, Vol.24, No.1, March 1998)
  //
  // Order is set to be 6 as suggested in the above reference.

  if(reduced_form){
    cout<<"# Not coded in G2D::CalMatrixExponentialByPade.\n";
    exit(0);
  };

  int order=6;

  if(x!=y){
    cout<<"# Error in GroupData2D::CalMatrixExponentialByPade.\n";
    cout<<"# GroupData2D should be square matrix.\n";
    exit(0);
  };

  real norm=InfiniteNorm()*factor;
  real m_norm=norm;

  // scaling the matrix to satisfy that the maximum norm is less than 0.5
  int counter=0;
  while(m_norm>0.5){
    m_norm*=0.5;
    factor*=0.5;
    counter++;
  };
  //cout<<" Counter : "<<counter<<"\n";

  int order2=order*2;

  GroupData2D tmat2=(*this)*(*this)*(factor*factor);

  vector<real> c(order+1);
  c[0]=1.;
  for(int i=1;i<=order;i++){
    c[i]=c[i-1]*(order+1-i)/((order2+1-i)*i);
  };

  GroupData2D mat1(x,x);
  GroupData2D mat2(x,x);
  mat1.set_zero();
  mat2.set_zero();

  for(int i=0;i<x;i++){
    mat1.put_data(i,i,c[1]);
    mat2.put_data(i,i,c[0]);
  };

  GroupData2D multmat=tmat2;
  for(int i=1;i<=order/2;i++){
    mat2=mat2+multmat*c[2*i];
    if(i!=order/2){
      mat1=mat1+multmat*c[2*i+1];
      multmat=multmat*tmat2;
    };
  };

  mat1=((*this)*factor)*mat1; // numerator
  mat2=mat2-mat1;          // denominator

  mat2.solveaxb_mod(mat1);

  mat1=mat1*2.;
  for(int i=0;i<x;i++){
    mat1.add_data(i,i,1.);
  };

  // squaring the scaled matrix
  for(int i=0;i<counter;i++){
    mat1=mat1*mat1;
  };

  return mat1;
};

GroupData2D GroupData2D::CalMatrixExponentialByPadeKatagiri(real factor)
{
  int N=y;
  int p=6;
  int q=6;

  real dt=factor;

  GroupData2D A(N,N);
  for(int i=0;i<N*N;i++){
    A.put_data(i,dat[i]);
  };


  GroupData2D Npq(N,N);
  GroupData2D Dpq(N,N);
  /*
  vector< vector<double> >Npq(N);
  vector< vector<double> >Dpq(N);
  for(int i=0;i<N;i++){
    Npq[i].resize(N,0.);
    Dpq[i].resize(N,0.);
  
  }
  */
  //累乗計算用行列B,C (B=A^k)
  GroupData2D B(N,N);

  /*
  vector< vector<double> >B(N);
  for(int i=0;i<N;i++){
    B[i].resize(N,0.);
  }
  */
  GroupData2D C(N,N);
  /*
  vector< vector<double> >C(N);
  for(int i=0;i<N;i++){
    C[i].resize(N);
  }
  */
  B.put_unit();
  /*
  for(int i=0;i<N;i++){
    B[i][i]=1.;
  };
  */
  //fa1=p!,fa2=q!,fa3=(p+q)!
  /*
  double fa1=mk_fact(p);
  double fa2=mk_fact(q);
  double fa3=mk_fact(p+q);
  */
  double fa1=kaijo(p,1);
  double fa2=kaijo(q,1);
  double fa3=kaijo(p+q,1);

  //Calculate Npq and Dpq
  for(int k=0;k<=p;k++){
    //B=A^k
    if(k==0){
      //unit(C,N);
      C.put_unit();
    }else{
      //multi(B,A,C,N);
      C=B*A;
    }
    B=C;
    //fa4=k!,fa5=(p+q-k)!,fa6=(p-k)! or (q-k)!
    /*
    double fa4=mk_fact(k);
    double fa5=mk_fact(p+q-k);
    double fa6=mk_fact(p-k);
    */
    double fa4=kaijo(k,1);
    double fa5=kaijo(p+q-k,1);
    double fa6=kaijo(p-k,1);

    //Npqの計算
    double temp=fa5*fa1/(fa3*fa4*fa6);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	//Npq.add_data(i,j,[i][j]+=temp*pow(dt,k)*B[i][j];
	Npq.add_data(i,j,temp*pow(dt,k)*B.get_dat(i,j));
	//Dpq[i][j]+=temp*pow((-1.*dt),k)*B[i][j];
	Dpq.add_data(i,j,temp*pow((-1.*dt),k)*B.get_dat(i,j));
      }
    }
  }
  //Dpqの逆行列xを求める
  GroupData2D x(N,N);
  /*
  vector< vector<double> >x(N);
  for(int i=0;i<N;i++){
    x[i].resize(N);
  }
  */
  //x=Dpq.inverse();
  x=Dpq.inverseKatagiri();
  //inverse(Dpq,x,N); 
  //multi(Npq,x,Mtex,N);
  GroupData2D Mtex=Npq*x;

  return Mtex;
};

GroupData2D GroupData2D::CalMatrixExponentialByTaylor(real factor)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalMatrixExponentialByTalyor.\n";
    exit(0);
  };

  int order_limit=60;

  GroupData2D tmat(x,x);
  for(int i=0;i<x*x;i++){
    tmat.put_data(i,dat[i]);
  };

  int order=tmat.CalOrderForMatrixExponential(factor);

  if(order>order_limit){
    cout<<"# Error in GroupData2D::CalMatrixExponentialByTaylor.\n";
    cout<<"# You have to set smaller factor.\n";
    cout<<"# Necessary order is "<<order<<"\n";
    exit(0);
  };

  if(x!=y){
    cout<<"# Error in GroupData2D::CalMatrixExponentialByTaylor.\n";
    cout<<"# GroupData2D should be square matrix.\n";
    exit(0);
  };

  tmat=tmat*factor;

  GroupData2D last(x,x);
  last.put_unit();

  GroupData2D mul(x,x);

  mul=last*tmat;
  last=last+mul;

  for(int i=0;i<=order;i++){
    mul=mul*tmat;
    real tmp=1./(i+2);
    mul=mul*tmp;
    last=last+mul;
  };

  return last;
};

int GroupData2D::CalOrderForMatrixExponential(real factor)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixExponential.\n";
    exit(0);
  };

  real maxv=0.;
  real maxv1=0.;
  for(int i=0;i<x;i++){
    real sum=0.;
    real sum1=0.;
    for(int j=0;j<x;j++){
      sum+=fabs(get_dat(i,j));
      sum1+=fabs(get_dat(j,i));
    };
    if(sum>maxv)maxv=sum;
    if(sum1>maxv1)maxv1=sum1;
  };
  if(maxv1<maxv)maxv=maxv1;

  // +++ temporal treatment
  real tmp=3.5*maxv*factor+5.;
  if(tmp>1e9)return 1000000000;

  int order=int(3.5*maxv*factor+5.);
  if(order<0){
    cout<<"#\n";
    cout<<"# Error in GroupData2D::CalOrderForMatrixExponential.\n";
    cout<<"# Negative order ("<<order<<") is detected.\n";
    cout<<"# Matrix norm defined by Lapidus&Luus : "<<maxv<<"\n";
    cout<<"# Factor : "<<factor<<"\n#\n";
    exit(0);
  };
  return order;
};

GroupData1D GroupData2D::CalMatrixExponentialByKrylov(GroupData1D &inp,real factor,int dim)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByKrylov.\n";
    exit(0);
  };

  // Matrix A          : self
  // Vector n          : inp
  // Scalar variable t : factor
  // 
  // This method calculate a vector n'=exp(At)n by the Krylov sub-space method
  // with the dimension `dim'

  if(x!=y){
    cout<<"# Error in GroupData2D::CalMatrixExponentialByKrylov.\n";
    cout<<"# GroupData2D should be square matrix.\n";
    exit(0);
  };

  int sz=inp.get_imax();
  if(sz!=x){
    cout<<"# Error in GroupData2D::CalMatrixExponentialByKrylov.\n";
    cout<<"# The sizes of matrix and vector are inconsistent.\n";
    exit(0);
  };

  int order=dim;
  int order1=order+1;

  vector<GroupData1D> v(order1);
  GroupData2D hmat(order1,order1);
  GroupData1D ret(sz);

  hmat.set_zero();
  v[0]=inp;
  real norm=v[0].GetEuclideanNorm();
  v[0].Multiply(1./norm);

  for(int j=0;j<order;j++){
    GroupData1D p=(*this)*v[j];
    for(int k=0;k<=j;k++){
      real tmp=v[k]*p;
      hmat.put_data(k,j,tmp);
      p=p-v[k]*tmp;
    };
    real tmp=p.GetEuclideanNorm();
    hmat.put_data(j+1,j,tmp);
    v[j+1]=p/tmp;
  };

  //GroupData2D tmp=hmat.CalMatrixExponentialByPade(factor);
  GroupData2D tmp=hmat.CalMatrixExponentialByChebyshev14(factor);

  ret.set_zero();
  for(int j=0;j<order1;j++){
    real ss=tmp.get_dat(j,0);
    if(ss!=0.){
      for(int k=0;k<sz;k++){
        real tmp=v[j].get_dat(k)*ss;
        ret.add_data(k,tmp);
      };
    };
  };

  ret.Multiply(norm);
  return ret;
};

GroupData2D GroupData2D::CalMatrixExponentialByChebyshev14(real factor)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByChebyshev14.\n";
    exit(0);
  };

  int num=x;

  real alfa0=1.83217437825404E-14;
  real alfa_Re[]={
    -7.15428806358907E-05,     9.43902531073617E-03,    -3.76360038782270E-01,    -2.34982320910827E+01,
     4.69332744888313E+01,    -2.78751619401456E+01,     4.80711209883251E+00};
  real alfa_Im[]={
     1.43610433495413E-04, -1.71847919584830E-02, 3.35183470294501E-01, -5.80835912971421E+00,
     4.56436497688278E+01, -1.02147339990565E+02, -1.32097938374287E+00};
  real theta_Re[]={
    -8.89777318646889E+00, -3.70327504942345E+00, -2.08758638250130E-01, 3.99336971057857E+00,
     5.08934506058063E+00, 5.62314257274598E+00, 2.26978382923111E+00};
  real theta_Im[]={
     1.66309826199021E+01, 1.36563718714833E+01, 1.09912605619013E+01, 6.00483164223504E+00,
     3.58882402902701E+00, 1.19406904634397E+00, 8.46173797304022E+00};
  real c1[]={
    3.55759950581330429031e+02,
    2.00210738783922550965e+02,
    1.20851388908650918097e+02,
    5.20050046969535344488e+01,
    3.87810910569781199797e+01,
    3.30455332808650794618e+01,
    7.67529281558298492882e+01,
  };
  real c2[]={
    -3.02495494850359234601e-03,
    2.69637216443245342212e-01,
    -3.76265726723711768642e+00,
    1.28715346976247644761e+02,
    -4.02666655759470813791e+02,
    2.78716986676928399902e+02,
    2.66676105789059536555e-01,
  };

  GroupData2D tmatdt2=(*this)*(*this)*(factor*factor);
  //tmatdt2.ReducedForm();
  GroupData2D tmp1(num,num);
  GroupData2D tmp2(num,num);
  GroupData2D tmp(num,num);
  tmp.set_zero();

  //cout.setf(ios::scientific);
  //cout.precision(20);
  for(int k=0;k<7;k++){
    tmp1=tmatdt2-(*this)*(2*theta_Re[k]*factor);
    //if(k==0)tmp1.ReducedForm();
    tmp2=(*this)*(alfa_Re[k]*factor);
    real coef1=theta_Re[k]*theta_Re[k]+theta_Im[k]*theta_Im[k];
    real coef2=-alfa_Re[k]*theta_Re[k]-alfa_Im[k]*theta_Im[k];
    for(int jj=0;jj<num;jj++){
      //tmp1.add_data(jj,jj,coef1);
      //tmp2.add_data(jj,jj,coef2);
      tmp1.add_data(jj,jj,c1[k]);
      tmp2.add_data(jj,jj,c2[k]);
    };
    tmp1.solveaxb_mod(tmp2);
    tmp=tmp+tmp2;
  };

  tmp=tmp*2;
  for(int i=0;i<num;i++){
    tmp.add_data(i,i,alfa0);
  };

  return tmp;
};

GroupData2D GroupData2D::CalMatrixExponentialByChebyshev16(real factor)
{ //kawamoto
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByChebyshev16.\n";
    exit(0);
  };

  int num=x;

  real alfa0=2.1248537104952237488e-16;
  real alfa_Re[]={-5.0901521865224915650e-7, 2.1151742182466030907e-4, 1.1339775178483930527e2, 1.5059585270023467528e1, -6.4500878025539646595e1, -1.4793007113557999718e0, -6.2518392463207918892e1, 4.1023136835410021273e-2 };

  real alfa_Im[]={-2.4220017652852287970e-5, 4.3892969647380673918e-3, 1.0194721704215856450e2, -5.7514052776421819979e0, -2.2459440762652096056e2, 1.7686588323782937906e0, -1.1190391094283228480e1, -1.5743466173455468191e-1 };

  real theta_Re[]={-1.0843917078696988026e1, -5.2649713434426468895e0, 5.9481522689511774808e0, 3.5091036084149180974e0, 6.4161776990994341923e0, 1.4193758971856659786e0, 4.9931747377179963991e0, -1.4139284624888862114e0 };
		   
  real theta_Im[]={1.9277446167181652284e1, 1.6220221473167927305e1, 3.5874573620183222829e0, 8.4361989858843750826e0, 1.1941223933701386874e0, 1.0925363484496722585e1, 5.9968817136039422260e0, 1.3497725698892745389e1 };

  GroupData2D tmatdt2=(*this)*(*this)*(factor*factor);
  //tmatdt2.ReducedForm();
  GroupData2D tmp1(num,num);
  GroupData2D tmp2(num,num);
  GroupData2D tmp(num,num);
  tmp.set_zero();

  //cout.setf(ios::scientific);
  //cout.precision(20);
  for(int k=0;k<8;k++){
    tmp1=tmatdt2-(*this)*(2*theta_Re[k]*factor);
    //if(k==0)tmp1.ReducedForm();
    tmp2=(*this)*(alfa_Re[k]*factor);
    real coef1=theta_Re[k]*theta_Re[k]+theta_Im[k]*theta_Im[k];
    real coef2=-alfa_Re[k]*theta_Re[k]-alfa_Im[k]*theta_Im[k];
    for(int jj=0;jj<num;jj++){
      tmp1.add_data(jj,jj,coef1);
      tmp2.add_data(jj,jj,coef2);
    };
    tmp1.solveaxb_mod(tmp2);
    tmp=tmp+tmp2;
  };

  tmp=tmp*2;
  for(int i=0;i<num;i++){
    tmp.add_data(i,i,alfa0);
  };

  return tmp;
};

GroupData1D GroupData2D::CalMatrixExponentialByChebyshev14(GroupData1D &inp, real factor)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByChebyshev14.\n";
    exit(0);
  };

  int num=x;

  real alfa0=1.83217437825404E-14;
  real alfa_Re[]={-7.15428806358907E-05, 9.43902531073617E-03, -3.76360038782270E-01, -2.34982320910827E+01, 4.69332744888313E+01,
		     -2.78751619401456E+01, 4.80711209883251E+00};
  real alfa_Im[]={1.43610433495413E-04, -1.71847919584830E-02, 3.35183470294501E-01, -5.80835912971421E+00, 4.56436497688278E+01,
		     -1.02147339990565E+02, -1.32097938374287E+00};
  real theta_Re[]={-8.89777318646889E+00, -3.70327504942345E+00, -2.08758638250130E-01, 3.99336971057857E+00, 5.08934506058063E+00,
		      5.62314257274598E+00, 2.26978382923111E+00};
  real theta_Im[]={1.66309826199021E+01, 1.36563718714833E+01, 1.09912605619013E+01, 6.00483164223504E+00, 3.58882402902701E+00,
		      1.19406904634397E+00, 8.46173797304022E+00};

  GroupData2D tmatdt2=(*this)*(*this)*(factor*factor);

  //tmatdt2.ReducedForm();

  GroupData1D tmat_inp=(*this)*inp*factor;

  GroupData2D tmp1(num,num);
  GroupData1D tmp2(num);
  GroupData1D tmp(num);
  tmp.set_zero();

  for(int k=0;k<7;k++){

    tmp1=tmatdt2-(*this)*(factor*2*theta_Re[k]);
    //if(k==0)tmp1.ReducedForm();
    tmp2=tmat_inp*alfa_Re[k];
    real coef1=theta_Re[k]*theta_Re[k]+theta_Im[k]*theta_Im[k];
    real coef2=-alfa_Re[k]*theta_Re[k]-alfa_Im[k]*theta_Im[k];
    for(int jj=0;jj<num;jj++){
      tmp1.add_data(jj,jj,coef1);
    };
    tmp2=tmp2+inp*coef2;
    tmp1.solveaxb_mod(tmp2);
    tmp=tmp+tmp2;

  };

  tmp=tmp*2+inp*alfa0;

  return tmp;
};

GroupData1D GroupData2D::CalMatrixExponentialByChebyshev16(GroupData1D &inp, real factor)
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByChebyshev16.\n";
    exit(0);
  };

  int num=x;

  real alfa0=2.1248537104952237488e-16;
  real alfa_Re[]={-5.0901521865224915650e-7, 2.1151742182466030907e-4, 1.1339775178483930527e2, 1.5059585270023467528e1, -6.4500878025539646595e1, -1.4793007113557999718e0, -6.2518392463207918892e1, 4.1023136835410021273e-2 };

  real alfa_Im[]={-2.4220017652852287970e-5, 4.3892969647380673918e-3, 1.0194721704215856450e2, -5.7514052776421819979e0, -2.2459440762652096056e2, 1.7686588323782937906e0, -1.1190391094283228480e1, -1.5743466173455468191e-1 };

  real theta_Re[]={-1.0843917078696988026e1, -5.2649713434426468895e0, 5.9481522689511774808e0, 3.5091036084149180974e0, 6.4161776990994341923e0, 1.4193758971856659786e0, 4.9931747377179963991e0, -1.4139284624888862114e0 };
		   
  real theta_Im[]={1.9277446167181652284e1, 1.6220221473167927305e1, 3.5874573620183222829e0, 8.4361989858843750826e0, 1.1941223933701386874e0, 1.0925363484496722585e1, 5.9968817136039422260e0, 1.3497725698892745389e1 };
		 
  GroupData2D tmatdt2=(*this)*(*this)*(factor*factor);

  //tmatdt2.ReducedForm();

  GroupData1D tmat_inp=(*this)*inp*factor;

  GroupData2D tmp1(num,num);
  GroupData1D tmp2(num);
  GroupData1D tmp(num);
  tmp.set_zero();

  for(int k=0;k<8;k++){

    tmp1=tmatdt2-(*this)*(factor*2*theta_Re[k]);
    //if(k==0)tmp1.ReducedForm();
    tmp2=tmat_inp*alfa_Re[k];
    real coef1=theta_Re[k]*theta_Re[k]+theta_Im[k]*theta_Im[k];
    real coef2=-alfa_Re[k]*theta_Re[k]-alfa_Im[k]*theta_Im[k];
    for(int jj=0;jj<num;jj++){
      tmp1.add_data(jj,jj,coef1);
    };
    tmp2=tmp2+inp*coef2;
    tmp1.solveaxb_mod(tmp2);
    tmp=tmp+tmp2;

  };

  tmp=tmp*2+inp*alfa0;

  return tmp;
};


/*
GroupData2D GroupData2D::CalMatrixExponentialByChebyshev(real factor,int order)
{
  GroupData2D tmat(x,x);
  for(int i=0;i<x*x;i++){
    tmat.put_data(i,dat[i]);
  };

  int num=this->get_x();
  GroupData2D nume(num,num);

  real alfa0;
  real alfa_Re[order/2];
  real alfa_Im[order/2];
  real theta_Re[order/2];
  real theta_Im[order/2];

  //16
  real alfa0_16=2.1248537104952237488E-16;
  real alfa_Re_16[]={-5.0901521865224915650E-7, 2.1151742182466030907E-4, 1.1339775178483930527E+2, 1.5059585270023467528E+1, -6.4500878025539646595E+1,
		     -1.4793007113557999718E+0, -6.2518392463207918892E+1, 4.1023136835410021273E-2};
  real alfa_Im_16[]={ -2.4220017652852287970E-5, 4.3892969647380673918E-3, 1.0194721704215856450E+2, -5.7514052776421819979E+0, -2.2459440762652096056E+2,
		      1.7686588323782937906E+0, -1.1190391094283228480E+1, -1.5743466173455468191E-1};
  real theta_Re_16[]={-1.0843917078696988026E+1, -5.2649713434426468895E+0, 5.9481522689511774808E+0, 3.5091036084149180974E+0, 6.4161776990994341923E+0,
		      1.4193758971856659786E+0, 4.9931747377179963991E+0, -1.4139284624888862114E+0};
  real theta_Im_16[]={1.9277446167181652284E+1, 1.6220221473167927305E+1, 3.5874573620183222829E+0, 8.4361989858843750826E+0, 1.1941223933701386874E+0,
		      1.0925363484496722585E+1, 5.9968817136039422260E+0, 1.3497725698892745389E+1};
  
  //14
  real alfa0_14=1.83217437825404E-14;
  real alfa_Re_14[]={-7.15428806358907E-05, 9.43902531073617E-03, -3.76360038782270E-01, -2.34982320910827E+01, 4.69332744888313E+01,
		     -2.78751619401456E+01, 4.80711209883251E+00};
  real alfa_Im_14[]={1.43610433495413E-04, -1.71847919584830E-02, 3.35183470294501E-01, -5.80835912971421E+00, 4.56436497688278E+01,
		     -1.02147339990565E+02, -1.32097938374287E+00};
  real theta_Re_14[]={-8.89777318646889E+00, -3.70327504942345E+00, -2.08758638250130E-01, 3.99336971057857E+00, 5.08934506058063E+00,
		      5.62314257274598E+00, 2.26978382923111E+00};
  real theta_Im_14[]={1.66309826199021E+01, 1.36563718714833E+01, 1.09912605619013E+01, 6.00483164223504E+00, 3.58882402902701E+00,
		      1.19406904634397E+00, 8.46173797304022E+00};
  nume.set_zero();
  
  if(order==16){
    alfa0=alfa0_16;
    for(int ii=0;ii<order/2;ii++){
      alfa_Re[ii]=alfa_Re_16[ii];
      alfa_Im[ii]=alfa_Im_16[ii];
      theta_Re[ii]=theta_Re_16[ii];
      theta_Im[ii]=theta_Im_16[ii];
    };
  }else if(order==14){
    alfa0=alfa0_14;
    for(int jj=0;jj<order/2;jj++){
      alfa_Re[jj]=alfa_Re_14[jj];
      alfa_Im[jj]=alfa_Im_14[jj];
      theta_Re[jj]=theta_Re_14[jj];
      theta_Im[jj]=theta_Im_14[jj];
    };
  }else{
    cout<<"error: dimension of CRA is invalid \n";
    exit(0);
  };

  GroupData2D tmatdt(num,num);
  tmatdt=tmat*factor;
  GroupData2D tmatdt2(num,num);
  tmatdt2=tmatdt*tmatdt;

  GroupData2D tmp1(num,num);
  GroupData2D tmp2(num,num);
  GroupData2D tmp(num,num);

  for(int k=0;k<order/2;k++){

    tmp1=tmatdt2-tmatdt*(2*theta_Re[k]);
    tmp2=tmatdt*alfa_Re[k];
    real coef1=theta_Re[k]*theta_Re[k]+theta_Im[k]*theta_Im[k];
    real coef2=-alfa_Re[k]*theta_Re[k]-alfa_Im[k]*theta_Im[k];
    for(int jj=0;jj<num;jj++){
      tmp1.add_data(jj,jj,coef1);
      tmp2.add_data(jj,jj,coef2);
    };

    tmp1.solveaxb_mod(tmp2);
    tmp=tmp+tmp2;
  };

  nume=tmp*2;
  for(int i=0;i<num;i++){
    nume.add_data(i,i,alfa0);
  };

  return nume;
};
*/

GroupData1D GroupData2D::CalMatrixExponentialByChebyshevNew14(GroupData1D &inp, real factor)//kawamoto
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalOrderForMatrixByChebyshevNew14.\n";
    exit(0);
  };

  int num=x;
  int order=14;                 

  real alfa=1.832174378254041562E-14;
  
  real a[]={
    1.779554637293771435E+01,
    7.406550098846913954E+00,
    4.175172765007157949E-01,
    -4.539567658462848598E+00,
    -7.986739421157043495E+00,
    -1.017869012116137561E+01,
    -1.124628514549169722E+01,
  };
  
  real b[]={ 
    3.557599505813295764E+02,
    2.002107387839255352E+02,
    1.208513889086467117E+02,
    7.675292815583104300E+01,
    5.200500469695239047E+01,
    3.878109105697956238E+01,
    3.304553328086434050E+01,
  };

  real c[]={ 
    -1.430857612707635935E-04,
    1.887805062146882124E-02,
    -7.527200775643080322E-01,
    9.614224197666130678E+00,
    -4.699646418217149346E+01,
    9.386654897766216266E+01,
    -5.575032388028618158E+01,
  };

  real d[]={ 
    -6.049909897011568338E-03,
    5.392744328863667835E-01,
    -7.525314534448264148E+00,
    5.333522115626030402E-01,
    2.574306939525277471E+02,
    -8.053333115210730284E+02,
    5.574339733538109840E+02,
  };  

  GroupData2D At2=(*this)*(*this)*(factor*factor);
  GroupData1D AtN=(*this)*inp*factor;
  GroupData1D nuc(num);
  GroupData1D tmp1(num);
  GroupData2D q(num,num);

    for(int k=0;k<num;k++){
      nuc.put_data(k,alfa*inp.get_dat(k));
    };

  for(int i=0;i<order/2;i++){
    q=At2+(*this)*(a[i]*factor);
    for(int j=0;j<num;j++){
      tmp1.put_data(j,inp.get_dat(j)*d[i]+AtN.get_dat(j)*c[i]);
      q.add_data(j,j,b[i]);
    };
    q.solveaxb_mod(tmp1);
    nuc=nuc+tmp1;
  };
  return nuc;
};

GroupData2D GroupData2D::CalMatrixExponentialByChebyshevNew14(real factor)//kawamoto
{
  if(reduced_form){
    cout<<"# Not coded in G2D::CalMatrixExpoinentialByChebysjevNwe14.\n";
    exit(0);
  };

  int num=x;
  int order=14;                 

  real alfa=1.832174378254041562E-14;
  
  real a[]={
    1.779554637293771435E+01,
    7.406550098846913954E+00,
    4.175172765007157949E-01,
    -4.539567658462848598E+00,
    -7.986739421157043495E+00,
    -1.017869012116137561E+01,
    -1.124628514549169722E+01,
  };
  
  real b[]={ 
    // (Taken from original routine)
    /*
    3.55759950581330429031e+02,
    2.00210738783922550965e+02,
    1.20851388908650918097e+02,
    7.67529281558298492882e+01,
    5.20050046969535344488e+01,
    3.87810910569781199797e+01,
    3.30455332808650794618e+01,
    */
    // (Kawamoto-kun original)
    3.557599505813295764E+02,
    2.002107387839255352E+02,
    1.208513889086467117E+02,
    7.675292815583104300E+01,
    5.200500469695239047E+01,
    3.878109105697956238E+01,
    3.304553328086434050E+01,

  };

  real c[]={ 
    -1.430857612707635935E-04,
    1.887805062146882124E-02,
    -7.527200775643080322E-01,
    9.614224197666130678E+00,
    -4.699646418217149346E+01,
    9.386654897766216266E+01,
    -5.575032388028618158E+01,
  };

  real d[]={ 
    // (Taken from original routine)
    /*
    -3.02495494850359234601e-03*2,
    2.69637216443245342212e-01*2,
    -3.76265726723711768642e+00*2,
    2.66676105789059536555e-01*2,
    1.28715346976247644761e+02*2,
    -4.02666655759470813791e+02*2,
    2.78716986676928399902e+02*2,
    */
    // (Kawamoto-kun original)
    -6.049909897011568338E-03,
    5.392744328863667835E-01,
    -7.525314534448264148E+00,
    5.333522115626030402E-01,
    2.574306939525277471E+02,
    -8.053333115210730284E+02,
    5.574339733538109840E+02,

  };  

  GroupData2D At2=(*this)*(*this)*(factor*factor);

  GroupData2D ret(num,num);
  ret.set_zero();
  GroupData2D q(num,num);
  GroupData2D q2(num,num);

  for(int i=0;i<order/2;i++){
    q=At2+(*this)*(a[i]*factor);
    q2=(*this)*(c[i]*factor);
    for(int j=0;j<num;j++){
      q.add_data(j,j,b[i]);
      q2.add_data(j,j,d[i]);
    };
    q.solveaxb_mod(q2);
    ret=ret+q2;
  };

  for(int i=0;i<num;i++){
    ret.add_data(i,i,alfa);
  };

  return ret;
};

void GroupData2D::show_plot(real shift)
{
  if(reduced_form)NormalForm();

  int ind1=0;
  int ind2=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      cout<<j+0.5+shift<<" "<<i+0.5+shift<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      cout<<j+1.5+shift<<" "<<i+0.5+shift<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      ind1++;
    };
    cout<<"\n";
    for(int j=0;j<x;j++){
      cout<<j+0.5+shift<<" "<<i+1.5+shift<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      cout<<j+1.5+shift<<" "<<i+1.5+shift<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      ind2++;
    };
    cout<<"\n";
  };
};

void GroupData2D::show_plot(string dir,string filename)
{
  if(reduced_form)NormalForm();

  string output=dir+filename;
  ofstream fout;
  fout.open(output.data(),ios::out);
  if(fout.fail()){
    cout<<"# Failed to open the file for plotting\n";
    cout<<"# File name is "<<output<<"\n";
    exit(1);
  };
  fout.setf(ios::scientific);                                                   
  fout.precision(6);  

  int ind1=0;
  int ind2=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      fout<<j+0.5<<" "<<i+0.5<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      fout<<j+1.5<<" "<<i+0.5<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      ind1++;
    };
    fout<<"\n";
    for(int j=0;j<x;j++){
      fout<<j+0.5<<" "<<i+1.5<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      fout<<j+1.5<<" "<<i+1.5<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      ind2++;
    };
    fout<<"\n";
  };
  fout.close();
};

void GroupData2D::show_plot(GroupData1D &val)
{
  int imax=val.get_imax();
  if(imax!=x+1||imax!=y+1){
    cout<<"# Error in GroupData2D::show_plot.\n";
    exit(0);
  };
  
  if(reduced_form)NormalForm();

  int ind1=0;
  int ind2=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      cout<<val.get_dat(j)<<" "<<val.get_dat(i)<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      cout<<val.get_dat(j+1)<<" "<<val.get_dat(i)<<" "<<dat[ind1]<<" "<<fabs(dat[ind1])<<"\n";
      ind1++;
    };
    cout<<"\n";
    for(int j=0;j<x;j++){
      cout<<val.get_dat(j)<<" "<<val.get_dat(i+1)<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      cout<<val.get_dat(j+1)<<" "<<val.get_dat(i+1)<<" "<<dat[ind2]<<" "<<fabs(dat[ind2])<<"\n";
      ind2++;
    };
    cout<<"\n";
  };
};


GroupData1D GroupData2D::CalMatrixExponentialByMMPA18(GroupData1D &inp, real factor)//kawamoto
{
  if(reduced_form){
    cout<<"# [Reduced_form] Not coded in MMPA18.\n";
    exit(0);
  };
  int num=x;
  int order=18;
  real c=11.6;
  real a[]={
    9.15932623819800e-6,    -2.12690490494072e-4,
    2.25558826845951e-3,    -1.43527517394836e-2,
    6.05493236636062e-2,    -1.75678306731683e-1,
    3.46682453500716e-1,    -4.27402177134336e-1,
    2.27616348870777e-1,    1.42342618503802e-1,
    -2.40743301027114e-1,    -2.57330703975956e-2,
    1.66542198829366e-1,    -2.96247724265369e-4,
    -8.91182385704820e-2,    1.67265642326255e-3,
    3.14353670403391e-2,    -3.40028968333583e-4,
    -5.22890515538606e-3,
  };


  GroupData1D nuc=inp*a[0];
  GroupData1D tmp1=inp;

  /*
  GroupData2D S=(*this)*factor;
  for(int k=0;k<num;k++){
    S.add_data(k,k,-c);
  };

  S.DoInverse();

  S=S*(c*2.);
  for(int k=0;k<num;k++){
    S.add_data(k,k,1.);
  };

  S.ReducedForm();

  for(int i=1;i<order+1;i++){
    tmp1=S*tmp1;
    nuc=nuc+tmp1*a[i];
  };
  */

  GroupData2D S=(*this)*factor;
  for(int k=0;k<num;k++){
    S.add_data(k,k,-c);
  };

  S.LUdecomposition();
  S.ReducedForm();

  GroupData1D tmp2;
  for(int i=1;i<=order;i++){
    tmp2=tmp1;
    S.solveaxb_LUdecomposition(tmp1);
    tmp1=tmp1*(c*2.)+tmp2;
    nuc=nuc+tmp1*a[i];
  };


  return nuc;
};

GroupData1D GroupData2D::CalMatrixExponentialByMMPA32(GroupData1D &inp, real factor) //kawamoto
{
  if(reduced_form){
    cout<<"# [Reduced_form] Not coded in MMPA32.\n";
    exit(0);
  };
  int num=x;

  /*
  int order=10;
  real c=12.504; // case 1
  real a[]={
-9.21877742692655e-6, -2.92352400447805e-5, 2.08755397019021e-3, -9.16869195193606e-3, 
2.17756357601593e-2, -9.99888786617314e-2, 3.25186434321652e-1, -4.59443792425971e-1, 
2.13857621214331e-1, 6.86306899307289e-2, -6.29118198603515e-2, 
  };
  */
  /*
  int order=10;
  real c=12.463; // case 2
  real a[]={
-9.72202507517625e-6, -4.02029454492449e-5, 2.16797450524908e-3, -9.19815043716278e-3, 
2.22288912532799e-2, -1.03839351834854e-1, 3.30717930187500e-1, -4.56867790548517e-1, 
2.06135055386484e-1, 6.99526521406809e-2, -6.12472842485857e-2, 
  };
  */
  /*
  int order=10;
  real c=12.559; // case 3
  real a[]={
-9.82261787454618e-6, -1.20955603325446e-5, 2.06199002444275e-3, -9.18159043214877e-3, 
2.04359625825662e-2, -9.47098208045513e-2, 3.20245985245503e-1, -4.63006636438534e-1, 
2.20885997171252e-1, 6.69028624680187e-2, -6.36274904109506e-2, 
  };
  */
  /*
  int order=12;
  real c=11.8; // original
  real a[]={
    5.75675518822064e-6,    -1.74303098865960e-4,    2.10039799323683e-3,    -1.25401951057288e-2,
    5.04542890245522e-2,    -1.59170661868637e-1,    3.50070832130723e-1,    -4.33870454060926e-1,
    1.90537907279739e-1,    1.21316736259960e-1,    -1.18763288443338e-1,    -1.55621553712551e-2,
    2.55948302640294e-2,
  };
  */
  /*
  int order=12;
  real c=11.842; // case 1
  real a[]={
5.52393846542740e-6, -1.66286523474290e-4, 2.02859693896894e-3, -1.22176080888358e-2, 
4.92067179119885e-2, -1.55583035297711e-1, 3.45648187460189e-1, -4.35751950239325e-1, 
1.98758334856865e-1, 1.18207950705269e-1, -1.21858627448792e-1, -1.44890684240145e-2, 
2.62129597060598e-2, 
  };
  */
  /*
  int order=12;
  real c=11.808; // case 2
  real a[]={
5.71124053372124e-6, -1.72749057027109e-4, 2.08657801668722e-3, -1.24781332160876e-2, 
5.02137647493151e-2, -1.58482426805889e-1, 3.49233116408566e-1, -4.34244343855078e-1, 
1.92103250520163e-1, 1.20738803153321e-1, -1.19353351005513e-1, -1.53620240410216e-2, 
2.57118023356741e-2, 
  };
  */
  /*
  int order=12;
  real c=11.881; // case 3
  real a[]={
5.19835983571307e-6, -1.58932287192293e-4, 1.97440037640449e-3, -1.19306326423362e-2, 
4.79315685945573e-2, -1.52270388828235e-1, 3.42228226079096e-1, -4.37423479701040e-1, 
2.04774844519067e-1, 1.15279004316642e-1, -1.23114376772990e-1, -1.34946936997478e-2, 
2.62010246122738e-2, 
  };
  */
  /*
  int order=12;
  real c=11.88124; // case 4
  real a[]={
5.19639720874541e-6, -1.58887455974832e-4, 1.97407148544813e-3, -1.19289002583794e-2, 
4.79238202896580e-2, -1.52250075898093e-1, 3.42207028551248e-1, -4.37433515001117e-1, 
2.04811856320297e-1, 1.15260842884047e-1, -1.23122052981765e-1, -1.34885826765661e-2, 
2.62009616934767e-2, 
  };
  */
  /*
  int order=14;
  real c=11.580; // case 1
  real a[]={
9.10351898599396e-6, -2.16567497056083e-4, 2.32624126404127e-3, -1.45616399828950e-2, 
6.05149543724722e-2, -1.77230576579051e-1, 3.55085119424283e-1, -4.27383649368767e-1, 
1.95498365124026e-1, 1.46648357413824e-1, -1.71139053827321e-1, -3.06159546259200e-2, 
7.13743645942142e-2, 3.36002848425238e-3, -1.36693400515054e-2, 
  };
  */
  /*
  int order=14;
  real c=11.609; // case 3
  real a[]={
8.82918100919415e-6, -2.10724798492612e-4, 2.27325537234040e-3, -1.42714200182412e-2, 
5.94529602592659e-2, -1.74824175340244e-1, 3.52544816049316e-1, -4.28338244296778e-1, 
2.00542153602190e-1, 1.43289776925379e-1, -1.72638060687039e-1, -2.85171976598181e-2, 
7.14032612247240e-2, 2.87185889188562e-3, -1.35873437314569e-2, 
  };
  */
  /*
  int order=16;
  real c=11.5; // original
  real a[]={
    1.00893480690695e-5,    -2.33129670152107e-4,    2.45359574039542e-3,    -1.54079862484901e-2,
    6.40874197350830e-2,    -1.83794037592770e-1,    3.57381552529419e-1,    -4.25312835414091e-1,
    2.00116659272177e-1,    1.58124659927456e-1,    -2.07785492608093e-1,    -3.99592478630238e-2,
    1.20978880327727e-1,    7.34703261312408e-3,    -4.52608511768821e-2,    -7.64479984552456e-4,
    8.01816399829783e-3,
  };
  */
  /*
  int order=16;
  real c=11.530; // case 1
  real a[]={
9.79098498941607e-6, -2.26796650297272e-4, 2.39413793031682e-3, -1.50843371358783e-2, 
6.29622439288415e-2, -1.81325866194398e-1, 3.54649968211641e-1, -4.26096864592222e-1, 
2.05996041185721e-1, 1.53720285121783e-1, -2.10557170632784e-1, -3.62441946317224e-2, 
1.22147408195900e-1, 5.68791917807874e-3, -4.57238129607344e-2, -4.30145411948643e-4, 
8.12143291695198e-3, 
  };
  */
  /*
  int order=16;
  real c=11.554; // case 3
  real a[]={
9.55673616071427e-6, -2.21853699988005e-4, 2.34787382195520e-3, -1.48300407009397e-2, 
6.20676940819214e-2, -1.79367744030110e-1, 3.52517505308707e-1, -4.26660649630603e-1, 
2.10376440686283e-1, 1.50129630178696e-1, -2.12038277997235e-1, -3.32350699958465e-2, 
1.22178229111490e-1, 4.34479734899626e-3, -4.55075681146358e-2, -1.59048856859936e-4, 
8.04856677593516e-3, 
  };
  */
  /*
  int order=16;
  real c=11.55385; // case 4
  real a[]={
9.55818297829682e-6, -2.21884266813879e-4, 2.34816021434733e-3, -1.48316170495578e-2, 
6.20732505511125e-2, -1.79379941734168e-1, 3.52530860768316e-1, -4.26657273240572e-1, 
2.10349191740789e-1, 1.50152179516416e-1, -2.12029252834189e-1, -3.32538652145767e-2, 
1.22178202678754e-1, 4.35315041682305e-3, -4.55089821578535e-2, -1.60727919677468e-4, 
8.04903136564369e-3, 
  };
  */
  /*
  int order=18;
  real c=11.6; // original
  real a[]={
    9.15932623819800e-6,    -2.12690490494072e-4,    2.25558826845951e-3,    -1.43527517394836e-2,
    6.05493236636062e-2,    -1.75678306731683e-1,    3.46682453500716e-1,    -4.27402177134336e-1,
    2.27616348870777e-1,    1.42342618503802e-1,    -2.40743301027114e-1,    -2.57330703975956e-2,
    1.66542198829366e-1,    -2.96247724265369e-4,    -8.91182385704820e-2,    1.67265642326255e-3,
    3.14353670403391e-2,    -3.40028968333583e-4,    -5.22890515538606e-3,
  };
  */
  /*
  int order=18;
  real c=14.765; // case 3
  real a[]={
3.93648518289373e-7, -1.14346418748559e-5, 1.55838542734623e-4, -1.33408701176120e-3, 
7.82402013592279e-3, -3.26792965130909e-2, 1.00433361252508e-1, -2.29881874452124e-1, 
3.75775864495606e-1, -3.78619734678191e-1, 1.17380630920616e-1, 1.88825526846999e-1, 
-1.73268794001677e-1, -5.62431384343664e-2, 1.01127470786705e-1, 1.10286742262245e-2, 
-3.50949715203808e-2, -1.08463195091157e-3, 5.66618907058401e-3, 
  };
  */
  /*
  int order=20;
  real c=14.7; // original
  real a[]={
    4.14011494076567e-7,    -1.21389050176560e-5,    1.66036007429611e-4,    -1.40421255624477e-3,
    8.14300954554650e-3,    -3.39547425075120e-2,    1.04171615631136e-1,    -2.35238700580452e-1,
    3.76260500179568e-1,    -3.74210537031247e-1,    1.18096804136174e-1,    1.95228904688230e-1,
    -1.96413784679446e-1,    -6.29078572280224e-2,    1.40871764292171e-1,    1.44334184733382e-2,
    -7.04623701868630e-2,    -2.04779614779486e-3,    2.25739021964080e-2,    1.13661125146827e-4,
    -3.40789155164420e-3,
  };
  */
  /*
  int order=20;
  real c=14.674; // case 1
  real a[]={
4.24857702041211e-7, -1.24356984528263e-5, 1.69798689702040e-4, -1.43304533353300e-3, 
8.29043424421308e-3, -3.44814141423497e-2, 1.05487738941498e-1, -2.37362338631114e-1, 
3.77746119561873e-1, -3.72499976486242e-1, 1.13307097076996e-1, 1.98092208989805e-1, 
-1.93999562069480e-1, -6.61993177080266e-2, 1.40290899941064e-1, 1.64391183765991e-2, 
-7.05536484728613e-2, -2.77902209642216e-3, 2.27007917542749e-2, 2.36222727592072e-4, 
-3.44009558252788e-3, 
  };
  */
  /*
  int order=20;
  real c=14.694; // case 3
  real a[]={
4.16489503039799e-7, -1.22067562043404e-5, 1.66896965309797e-4, -1.41081632404875e-3, 
8.17680801326161e-3, -3.40756211264286e-2, 1.04474162428242e-1, -2.35728167632626e-1, 
3.76605537674522e-1, -3.73820531958901e-1, 1.16993194423470e-1, 1.95894618054336e-1, 
-1.95861943840487e-1, -6.36687874035241e-2, 1.40741897132615e-1, 1.48955635489777e-2, 
-7.04850635639192e-2, -2.21579173173015e-3, 2.26033711398655e-2, 1.41740782915132e-4, 
-3.41527739624383e-3, 
  };
  */
  /*
  int order=20;
  real c=14.69369; // case 4
  real a[]={
4.16617937693700e-7, -1.22102722125430e-5, 1.66941567355808e-4, -1.41115832902076e-3, 
8.17855788134105e-3, -3.40818772644543e-2, 1.04489813131612e-1, -2.35753466679296e-1, 
3.76623328898622e-1, -3.73800304490707e-1, 1.16936146488922e-1, 1.95928932846937e-1, 
-1.95833345799492e-1, -6.37080802397047e-2, 1.40735119353696e-1, 1.49194527081002e-2, 
-7.04862096272195e-2, -2.22448368710266e-3, 2.26048906890670e-2, 1.43194867024605e-4, 
-3.41565974216751e-3, 
  };
  */
  /*
  int order=24;
  real c=17.8; // original
  real a[]={
    1.85727706362754e-8,-6.62293286498090e-7,    1.11359131965495e-5,    -1.16959490572316e-4,
    8.59372414386465e-4,    -4.67728379562120e-3,    1.93964760373493e-2,    -6.19986270438348e-2,
    1.52693867084004e-1,    -2.84754109473843e-1,    3.80107454619660e-1,    -3.03861908580259e-1,
    2.26808786137166e-2,    2.29611295173640e-1,    -1.62026209814790e-1,    -1.00944100073890e-1,
    1.46610569213925e-1,    3.43074053045852e-2,    -9.05555737864560e-2,    -9.08731145408663e-3,
    3.99341579793824e-2,    1.68400175634436e-3,    -1.11859466876513e-2,    -1.61740046601756e-4,
    1.47379985234632e-3,
  };
  */
  /*
  int order=24;
  real c=17.822; // case 1
  real a[]={
1.81685608919666e-8, -6.48658453502610e-7, 1.09213669963751e-5, -1.14869101156187e-4, 
8.45261222168201e-4, -4.60793141616206e-3, 1.91440448596518e-2, -6.13196553942114e-2, 
1.51384269730683e-1, -2.83152238655078e-1, 3.79600733352948e-1, -3.06063169215033e-1, 
2.66415235806343e-2, 2.28147794823381e-1, -1.64916841456689e-1, -9.83839608190727e-2, 
1.48101848030460e-1, 3.23270464274958e-2, -9.12476718026753e-2, -8.09016644857577e-3, 
4.02237479421622e-2, 1.37520075988067e-3, -1.12754760015959e-2, -1.17402302536026e-4, 
1.48762103485865e-3, 
  };
  */
  /*
  int order=24;
  real c=17.839; // case 3
  real a[]={
1.78613053016960e-8, -6.38314083252574e-7, 1.07586966655067e-5, -1.13279054736428e-4, 
8.34495784240667e-4, -4.55501307755088e-3, 1.89513844164972e-2, -6.07993443474382e-2, 
1.50374037631736e-1, -2.81913568778984e-1, 3.79214622032284e-1, -3.07739997475025e-1, 
2.96184750642083e-2, 2.26981197943787e-1, -1.66915512334619e-1, -9.63816016197855e-2, 
1.48872467262053e-1, 3.07849858129218e-2, -9.13724459440875e-2, -7.31511902587161e-3, 
4.01657364905896e-2, 1.13534658215509e-3, -1.12329332686242e-2, -8.29686305210943e-5, 
1.47889632201372e-3, 
  };
  */
  /*
  int order=24;
  real c=17.83855; // case 4
  real a[]={
1.78693714038827e-8, -6.38585790751995e-7, 1.07629714909868e-5, -1.13320862128423e-4, 
8.34779025956038e-4, -4.55640640259608e-3, 1.89564611943183e-2, -6.08130681266726e-2, 
1.50400720179418e-1, -2.81946369641136e-1, 3.79225022145313e-1, -3.07695865568209e-1, 
2.95396776760794e-2, 2.27012429235457e-1, -1.66862871510024e-1, -9.64347808786503e-2, 
1.48852351679059e-1, 3.08258225095540e-2, -9.13693065612238e-2, -7.33559838313701e-3, 
4.01673232082259e-2, 1.14167163667979e-3, -1.12340640574192e-2, -8.38749188089628e-5, 
1.47912619399230e-3, 
  };
  */
  /*
  int order=28;
  real c=21.0; // original
  real a[]={
    7.59063464648194e-10,    -3.18471237173039e-8,    6.36552759591466e-7,    -8.05717433289618e-6,
    7.22266586138094e-5,    -4.86104590408253e-4,    2.54347940392935e-3,    -1.05573917534379e-2,
    3.50580072584662e-2,    -9.30230994107252e-2,    1.95068488371746e-1,    -3.14335794211098e-1,
    3.62206739642819e-1,    -2.34106034575124e-1,    -4.66000739481099e-2,    2.41768334946413e-1,
    -1.29556144006507e-1,    -1.26611252942053e-1,    1.47958394691771e-1,    4.95442624916426e-2,
    -1.07042173758764e-1,    -1.51304261858324e-2,    5.78917380013658e-2,    3.40138693075018e-3,
    -2.25100872099550e-2,    -4.87631594671312e-4,    5.55539967403326e-3,    3.18399153687889e-5,
    -6.46632091406698e-4,
  };
  */
  /*
  int order=28;
  real c=20.971; // case 1
  real a[]={
7.81352276011718e-10, -3.27378898270847e-8, 6.53449142261731e-7, -8.25829803387529e-6, 
7.39064253697677e-5, -4.96546485152077e-4, 2.59330905409602e-3, -1.07420736427465e-2, 
3.55882962046706e-2, -9.41810049642365e-2, 1.96884129467129e-1, -3.15990664771596e-1, 
3.61885744329124e-1, -2.30596344656369e-1, -5.11499884459546e-2, 2.42586763807499e-1, 
-1.25663668163141e-1, -1.29562089271278e-1, 1.46006732247775e-1, 5.23545883356463e-2, 
-1.06422717602898e-1, -1.68902310399588e-2, 5.78444652020299e-2, 4.15971011701345e-3, 
-2.25789554632392e-2, -6.91168568972892e-4, 5.59078855955660e-3, 5.73521760897817e-5, 
-6.52696045783164e-4, 
  };
  */
  /*
  int order=28;
  real c=20.986; // case 3
  real a[]={
7.69742556743028e-10, -3.22740750970350e-8, 6.44654897485491e-7, -8.15365795857101e-6, 
7.30328297366917e-5, -4.91118463575971e-4, 2.56741974867830e-3, -1.06461771098932e-2, 
3.53131264203556e-2, -9.35806819686082e-2, 1.95944048276983e-1, -3.15136350798709e-1, 
3.62057190688311e-1, -2.32417296370276e-1, -4.87981590610022e-2, 2.42173666284266e-1, 
-1.27684330099748e-1, -1.28042713853379e-1, 1.47027610334535e-1, 5.09033020455038e-2, 
-1.06753539064323e-1, -1.59797495424943e-2, 5.78756378197467e-2, 3.76675743306127e-3, 
-2.25463187423362e-2, -5.85546447906277e-4, 5.57329622803850e-3, 4.40947236556885e-5, 
-6.49660804022675e-4, 
  };
  */
  /*
  int order=28;
  real c=20.98637; // case 4
  real a[]={
7.69458396243680e-10, -3.22627179248301e-8, 6.44439457424640e-7, -8.15109343402427e-6, 
7.30114104026411e-5, -4.90985306962261e-4, 2.56678427090052e-3, -1.06438217432729e-2, 
3.53063629378636e-2, -9.35659121096769e-2, 1.95920884545500e-1, -3.15115232698687e-1, 
3.62061277744097e-1, -2.32462061889856e-1, -4.87401222134409e-2, 2.42163199304149e-1, 
-1.27733928410003e-1, -1.28005047860850e-1, 1.47052397448686e-1, 5.08674398882920e-2, 
-1.06761307507097e-1, -1.59572970902370e-2, 5.78761288109387e-2, 3.75708391755245e-3, 
-2.25453782741934e-2, -5.82950429866301e-4, 5.57282450020625e-3, 4.37693751702772e-5, 
-6.49580473175483e-4, 
  };
  */
  // The following are for sensitivity study in the order of 20
  /*
  real c=14.6;
  real a[]={
    4.57476599790406e-7,    -1.33205505181424e-5,    1.80936767730830e-4,    -1.51829338925222e-3,
    8.72538032480633e-3,    -3.60210555660917e-2,    1.09287112854811e-1,    -2.43448654769803e-1,
    3.81972697610476e-1,    -3.67306118016262e-1,    9.90145965908058e-2,    2.05837230467397e-1,
    -1.85424584789771e-1,    -7.53140213305319e-2,    1.36287321330401e-1,    2.20207851763837e-2,
    -6.88056300898309e-2,    -4.81240228544392e-3,    2.20927649845548e-2,    5.75851403598567e-4,
    -3.33105306596497e-3,
  };
  */
  /*
  real c=14.8;
  real a[]={
    3.74824492754109e-7,    -1.10624247968531e-5,    1.52314424904242e-4,    -1.29842906339757e-3,
    7.59885530874149e-3,    -3.19973711892996e-2,    9.92322385743807e-2,    -2.27137125431971e-1,
    3.70321756867182e-1,    -3.80293176126087e-1,    1.36330035306322e-1,    1.83704576193222e-1,
    -2.05140441183036e-1,    -5.01119048728325e-2,    1.42669176042810e-1,    6.79261402572323e-3,
    -6.99413478145378e-2,    6.89565845493742e-4,    2.20652279894855e-2,    -3.37688157839414e-4,
    -3.28819034074573e-3,
  };
  */
  /*
  real c=10.4;
  real a[]={
    3.04338005097814e-5,    -6.33064149188569e-4,    5.94984322658747e-3,    -3.31044924719991e-2,
    1.19544442585883e-1,    -2.86259877046089e-1,    4.35141479868454e-1,    -3.41591925772565e-1,
    -3.31524380616092e-2,    2.77638423658223e-1,    -6.20128778826213e-2,    -2.02665696669958e-1,
    5.61704658364902e-2,    1.42069570226044e-1,    -3.02929895180450e-2,    -7.90646556185228e-2,
    1.04916457243738e-2,    2.82961538730699e-2,    -2.00598629289046e-3,    -4.68443943477518e-3,
    1.35980712866816e-4,
  };
  */
  /*
  real c=20.4;
  real a[]={
    2.66060365106564e-9,    2.27011581281219e-8,    7.93657079543564e-7,    -1.95321920094646e-5,
    1.27828187247342e-4,    -6.12334591021742e-4,    3.60175627502864e-3,    -1.65322176953263e-2,
    4.90987632935842e-2,    -1.10598950272923e-1,    2.26394622716139e-1,    -3.76707227572211e-1,
    3.70635324242444e-1,    -7.93359962589959e-2,    -1.90357999229662e-1,    1.25650944311614e-1,
    4.50884940930692e-2,    -4.98353179553064e-2,    -4.45205137195922e-3,    7.99061340063306e-3,
    -1.37534523573306e-4,
  };
  */
  /*
  real c=25.2;
  real a[]={
    -4.02346413240874e-9,    -1.49458957107863e-7,    9.43676494780661e-7,    1.11118063331337e-5,
    -3.36103713055002e-5,    -2.72998544578544e-4,    6.85675272238871e-4,    1.87246691615582e-3,
    -1.29145801530128e-3,    -2.60360267159299e-2,    5.61768270205278e-2,    -3.76069884306840e-2,
    1.09090894908767e-1,    -3.86914343399243e-1,    4.86320942942664e-1,    -1.35960627367562e-1,
    -1.71892756095758e-1,    1.03862647357972e-1,    2.07843594505945e-2,    -1.89551001985154e-2,
    1.58185234541995e-4,
  };
  */
  /*
  int order=32;
  real c=24.1; // original
  real a[]={
    3.41366034346810763441e-11,    -1.64648752802981232517e-9,    3.80465674795813706161e-8,    -5.59799382133376361074e-7,
    5.88558416643299794299e-6,    -4.69954781431498324372e-5,    2.95331252981155031018e-4,    -1.49342612561658940208e-3,
    6.16185418583114456893e-3,    -2.08884457395397287047e-2,    5.81666843649626684195e-2,    -1.31988170767468710376e-1,
    2.39594549098354306472e-1,    -3.34580286145275369457e-1,    3.25793940419810518947e-1,    -1.47667478685995070947e-1,
    -1.19600677614316509554e-1,    2.40742697969870662426e-1,    -7.60967222595908615166e-2,    -1.58220954609283804027e-1,
    1.32429922343368914663e-1,    7.69052976527847896194e-2,    -1.14189797653561609782e-1,    -3.06105871275143780464e-2,
    7.30162061991511174207e-2,    9.95181782801371819122e-3,    -3.55155454832496982106e-2,    -2.48870477891655009776e-3,
    1.23929263457145990096e-2,    4.21590406862371374909e-4,    -2.74111882594745761616e-3,    -3.57929539250313916141e-5,
    2.86523961626939505591e-4,
    };
  */
  /*
  int order=32; 
  real c=24.123; // case 1
  real a[]={
3.33603968111001e-11, -1.61056134726998e-9, 3.72540365236593e-8, -5.48714797987259e-7, 
5.77531862575176e-6, -4.61682401180472e-5, 2.90493057070386e-4, -1.47091812917168e-3, 
6.07769798994104e-3, -2.06357394925998e-2, 5.75660058883860e-2, -1.30897449958902e-1, 
2.38213589142823e-1, -3.33764661350338e-1, 3.26791085398036e-1, -1.50766912654514e-1, 
-1.16553371834021e-1, 2.41041598193768e-1, -7.95456409890457e-2, -1.56441472549635e-1, 
1.34749964700017e-1, 7.47237083326763e-2, -1.15466350483208e-1, -2.89668365069765e-2, 
7.36505995771813e-2, 9.06218569629274e-3, -3.57981513903182e-2, -2.15138123223478e-3, 
1.24947496569467e-2, 3.41455915497725e-4, -2.76594186109408e-3, -2.68576983853753e-5, 
2.89458541285386e-4, 
  };
  */
  /*
  int order=32;
  real c=24.110; // case 2
  real a[]={
3.37969273558912e-11, -1.63077020491425e-9, 3.76999469900845e-8, -5.54952949342634e-7, 
5.83738925910984e-6, -4.66340354871015e-5, 2.93218124722499e-4, -1.48359965929711e-3, 
6.12513063991533e-3, -2.07782308703719e-2, 5.79048789258495e-2, -1.31513200250985e-1, 
2.38994055340841e-1, -3.34227344177324e-1, 3.26230930858747e-1, -1.49017527766451e-1, 
-1.18278390534823e-1, 2.40879376679756e-1, -7.75984057689985e-2, -1.57453802886706e-1, 
1.33443639954646e-1, 7.59614388488051e-2, -1.14749873957365e-1, -2.98987859200041e-2, 
7.32955343648534e-2, 9.56655151750474e-3, -3.56401583905796e-2, -2.34267585242734e-3, 
1.24377856154004e-2, 3.86921139517713e-4, -2.75203091931670e-3, -3.19301828218256e-5, 
2.87810623115800e-4, 
  };
  */
  /*
  int order=32;
  real c=24.137; // case 3
  real a[]={
3.28959971440906e-11, -1.58907698043670e-9, 3.67800269049377e-8, -5.42074692460092e-7, 
5.70917242091975e-6, -4.56717209417432e-5, 2.87587045862597e-4, -1.45737775602171e-3, 
6.02698283862533e-3, -2.04832767565142e-2, 5.72031584837072e-2, -1.30236514350340e-1, 
2.37371585981818e-1, -3.33261463433834e-1, 3.27393244144877e-1, -1.52643992641465e-1, 
-1.14720590978893e-1, 2.41197856248353e-1, -8.15505849368947e-2, -1.55334715619717e-1, 
1.35976586525059e-1, 7.33810031429538e-2, -1.15989191830025e-1, -2.79590084943435e-2, 
7.37836559030159e-2, 8.51778211730999e-3, -3.57914009391911e-2, -1.94517252560520e-3, 
1.24720529223180e-2, 2.92495127184344e-4, -2.75697100292969e-3, -2.13996732373782e-5, 
2.88139857318398e-4, 
  }; 
  */

  int order=32;
  real c=24.13694; // case 4
  real a[]={
3.28979735872738e-11, -1.58916844491099e-9, 3.67820454991818e-8, -5.42102978872146e-7, 
5.70945430126668e-6, -4.56738376635576e-5, 2.87599439571872e-4, -1.45743552929181e-3, 
6.02719933663635e-3, -2.04839279890385e-2, 5.72047094494963e-2, -1.30239342084535e-1, 
2.37375193939405e-1, -3.33263630944967e-1, 3.27390685838193e-1, -1.52635962819998e-1, 
-1.14728463172040e-1, 2.41197225880137e-1, -8.15420042770902e-2, -1.55339492153136e-1, 
1.35971357413276e-1, 7.33867740449378e-2, -1.15986977045075e-1, -2.79633313486442e-2, 
7.37831013682762e-2, 8.52011363362550e-3, -3.57914361107692e-2, -1.94605446533440e-3, 
1.24721515959643e-2, 2.92704265189227e-4, -2.75700953444178e-3, -2.14229591230169e-5, 
2.88145489362067e-4, 
  };


  GroupData1D nuc=inp*a[0];
  GroupData1D tmp1=inp;

  GroupData2D S=(*this)*factor;
  for(int k=0;k<num;k++){
    S.add_data(k,k,-c);
  };
  S.LUdecomposition();
  S.ReducedForm();

  GroupData1D tmp2;
  for(int i=1;i<=order;i++){
    tmp2=tmp1;
    S.solveaxb_LUdecomposition(tmp1);
    tmp1=tmp1*(c*2.)+tmp2;
    nuc=nuc+tmp1*a[i];
  };

  /*
  GroupData2D S=(*this)*factor;
  for(int k=0;k<num;k++){
    S.add_data(k,k,-c);
  };
  S.DoInverse();
  S=S*(c*2.);
  for(int k=0;k<num;k++){
    S.add_data(k,k,1.);
  };
  S.ReducedForm();
  for(int i=1;i<order+1;i++){
    tmp1=S*tmp1;
    nuc=nuc+tmp1*a[i];
  };
  */

  return nuc;
};

void GroupData2D::LUDecompositionForMatrixExponentialByMMPA32(real factor)
{
  if(reduced_form){
    cout<<"# [Reduced_form] Not coded in MMPA32.\n";
    exit(0);
  };

  int num=x;
  //real c=24.1;
  //real c=24.13694; // n=32 and case 4
  real c=11.88124; // n=12 and case 4
  Factorize(factor);
  for(int k=0;k<num;k++){
    add_data(k,k,-c);
  };
  LUdecomposition();
  ReducedForm();
};

GroupData1D GroupData2D::CalMatrixExponentialByLUDMMPA32(GroupData1D &inp)
{
  int num=x;
  /*
  int order=32;
  real c=24.1;
  real a[]={
    3.41366034346810763441e-11,    -1.64648752802981232517e-9,
    3.80465674795813706161e-8,    -5.59799382133376361074e-7,
    5.88558416643299794299e-6,    -4.69954781431498324372e-5,
    2.95331252981155031018e-4,    -1.49342612561658940208e-3,
    6.16185418583114456893e-3,    -2.08884457395397287047e-2,
    5.81666843649626684195e-2,    -1.31988170767468710376e-1,
    2.39594549098354306472e-1,    -3.34580286145275369457e-1,
    3.25793940419810518947e-1,    -1.47667478685995070947e-1,
    -1.19600677614316509554e-1,    2.40742697969870662426e-1,
    -7.60967222595908615166e-2,    -1.58220954609283804027e-1,
    1.32429922343368914663e-1,    7.69052976527847896194e-2,
    -1.14189797653561609782e-1,    -3.06105871275143780464e-2,
    7.30162061991511174207e-2,    9.95181782801371819122e-3,
    -3.55155454832496982106e-2,    -2.48870477891655009776e-3,
    1.23929263457145990096e-2,    4.21590406862371374909e-4,
    -2.74111882594745761616e-3,    -3.57929539250313916141e-5,
    2.86523961626939505591e-4,
  };
  */
  /*
  int order=32;
  real c=24.13694; // case 4
  real a[]={
3.28979735872738e-11, -1.58916844491099e-9, 3.67820454991818e-8, -5.42102978872146e-7, 
5.70945430126668e-6, -4.56738376635576e-5, 2.87599439571872e-4, -1.45743552929181e-3, 
6.02719933663635e-3, -2.04839279890385e-2, 5.72047094494963e-2, -1.30239342084535e-1, 
2.37375193939405e-1, -3.33263630944967e-1, 3.27390685838193e-1, -1.52635962819998e-1, 
-1.14728463172040e-1, 2.41197225880137e-1, -8.15420042770902e-2, -1.55339492153136e-1, 
1.35971357413276e-1, 7.33867740449378e-2, -1.15986977045075e-1, -2.79633313486442e-2, 
7.37831013682762e-2, 8.52011363362550e-3, -3.57914361107692e-2, -1.94605446533440e-3, 
1.24721515959643e-2, 2.92704265189227e-4, -2.75700953444178e-3, -2.14229591230169e-5, 
2.88145489362067e-4, 
  };
  */
  int order=12;
  real c=11.88124; // case 4
  real a[]={    
5.19639720874541e-6, -1.58887455974832e-4, 1.97407148544813e-3, -1.19289002583794e-2, 
4.79238202896580e-2, -1.52250075898093e-1, 3.42207028551248e-1, -4.37433515001117e-1, 
2.04811856320297e-1, 1.15260842884047e-1, -1.23122052981765e-1, -1.34885826765661e-2, 
2.62009616934767e-2, 
  };


  GroupData1D nuc=inp*a[0];
  GroupData1D tmp1=inp;
  GroupData1D tmp2;
  for(int i=1;i<=order;i++){
    tmp2=tmp1;
    solveaxb_LUdecomposition(tmp1);
    tmp1=tmp1*(c*2.)+tmp2;
    nuc=nuc+tmp1*a[i];
  };

  return nuc;
};


GroupData2D GroupData2D::CalMatrixExponentialByMMPA32(real factor)//kawamoto
{
  if(reduced_form){
    cout<<"# Not coded in MMPA32.\n";
    exit(0);
  };

  int num=x;
  /*
  int order=32;
  real c=24.1;
  real a[]={
    3.41366034346810763441e-11,    -1.64648752802981232517e-9,
    3.80465674795813706161e-8,    -5.59799382133376361074e-7,
    5.88558416643299794299e-6,    -4.69954781431498324372e-5,
    2.95331252981155031018e-4,    -1.49342612561658940208e-3,
    6.16185418583114456893e-3,    -2.08884457395397287047e-2,
    5.81666843649626684195e-2,    -1.31988170767468710376e-1,
    2.39594549098354306472e-1,    -3.34580286145275369457e-1,
    3.25793940419810518947e-1,    -1.47667478685995070947e-1,
    -1.19600677614316509554e-1,    2.40742697969870662426e-1,
    -7.60967222595908615166e-2,    -1.58220954609283804027e-1,
    1.32429922343368914663e-1,    7.69052976527847896194e-2,
    -1.14189797653561609782e-1,    -3.06105871275143780464e-2,
    7.30162061991511174207e-2,    9.95181782801371819122e-3,
    -3.55155454832496982106e-2,    -2.48870477891655009776e-3,
    1.23929263457145990096e-2,    4.21590406862371374909e-4,
    -2.74111882594745761616e-3,    -3.57929539250313916141e-5,
    2.86523961626939505591e-4,
  };
  */
  int order=32;
  real c=24.13694; // case 4
  real a[]={
3.28979735872738e-11, -1.58916844491099e-9, 3.67820454991818e-8, -5.42102978872146e-7, 
5.70945430126668e-6, -4.56738376635576e-5, 2.87599439571872e-4, -1.45743552929181e-3, 
6.02719933663635e-3, -2.04839279890385e-2, 5.72047094494963e-2, -1.30239342084535e-1, 
2.37375193939405e-1, -3.33263630944967e-1, 3.27390685838193e-1, -1.52635962819998e-1, 
-1.14728463172040e-1, 2.41197225880137e-1, -8.15420042770902e-2, -1.55339492153136e-1, 
1.35971357413276e-1, 7.33867740449378e-2, -1.15986977045075e-1, -2.79633313486442e-2, 
7.37831013682762e-2, 8.52011363362550e-3, -3.57914361107692e-2, -1.94605446533440e-3, 
1.24721515959643e-2, 2.92704265189227e-4, -2.75700953444178e-3, -2.14229591230169e-5, 
2.88145489362067e-4, 
  };

  GroupData2D S=(*this)*factor;
  for(int k=0;k<num;k++){
    S.add_data(k,k,-c);
  };

  S.DoInverse();
  S.Factorize(c*2.);
  for(int k=0;k<num;k++){
    S.add_data(k,k,1.);
  };

  GroupData2D E(num,num);
  E.set_zero();
  for(int k=0;k<num;k++){
    E.put_data(k,k,a[0]);
  };

  E=E+S*a[1];
  GroupData2D M=S;
  //S.ReducedForm();
  for(int i=2;i<order+1;i++){
    M=S*M;
    E=E+M*a[i];
  };

  return E;
};

void GroupData2D::CalSimultaneousDifferentialEquationByDecomposedMatrixExponential
    (int i1,int i2,real dt,real theta,GroupData2D &mmat, GroupData2D &emat, GroupData2D &xmat,
     GroupData2D &mat0, GroupData2D &mat1, GroupData2D &mat2)
{
  int ii=i1+i2;

  GroupData2D a11=Slicing(0,i1-1,0,i1-1);
  GroupData2D a12=Slicing(0,i1-1,i1,ii-1);
  GroupData2D a21=Slicing(i1,ii-1,0,i1-1);
  GroupData2D a22=Slicing(i1,ii-1,i1,ii-1);
  GroupData2D a22inv=a22.inverse();

  GroupData2D imat2(i2,i2);
  imat2.put_unit();

  real dtinv=1./dt;

  mmat=a22.CalMatrixExponentialByChebyshev14(dt);
  emat=a22inv*(mmat+a22inv*dtinv*(imat2-mmat))*a21;
  xmat=a22inv*((imat2*-1.)-a22inv*dtinv*(imat2-mmat))*a21;

  /*
  GroupData2D tmp=a22inv*mmat;
  a21.show_self();
  tmp.show_self(); exit(0);

  a22inv.show_self();
  mmat.show_self();
  a21.show_self();
  emat.show_self(); exit(0);
  */

  GroupData2D imat1(i1,i1);
  imat1.put_unit();

  mat0=imat1-(a11+a12*xmat)*theta*dt;
  mat1=imat1+a12*emat*theta*dt+a11*(1.-theta)*dt;
  mat2=a12*mmat*theta*dt+a12*(1.-theta)*dt;
}; 


void GroupData2D::ReducedForm()
{
  if(reduced_form)return;

  reduced_form=true;
  if(x==0&&y==0)return;

  xpos.resize(y);
  val.resize(y);

  int ii=0;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      real tmp=dat[ii++];
      if(tmp!=0.){
	xpos[i].push_back(j);
	val[i].push_back(tmp);
      };
    };
  };

  vector<real> ().swap(dat);
  //dat.clear();
  //dat.resize(0);

};

void GroupData2D::NormalForm()
{
  if(!reduced_form)return;

  reduced_form=false;
  put_yx(y,x);

  int y_base=0;
  for(int i=0;i<y;i++){
    int sz=xpos[i].size();
    for(int j=0;j<sz;j++){
      dat[y_base+xpos[i][j]]=val[i][j];
    };
    y_base+=x;
  };

  vector< vector<int> > ().swap(xpos);
  vector< vector<real> > ().swap(val);

};

void GroupData2D::MultiStepCalc(GroupData1D &inp,vector<GroupData1D> &nuc,real factor,int substep)//kawamoto
{
  if(reduced_form){
    cout<<"Not coded in MSC.\n";
    exit(0);
  };

  MultiStepCalculation mult;
  if(substep<=6){
    mult.MultiStepCalc48th6stp((*this),inp,nuc,factor,substep);
  }else if(substep>6&&substep<=43){
    mult.MultiStepCalc48th43stp((*this),inp,nuc,factor,substep);
  }else if(substep>43&&substep<=275){
    mult.MultiStepCalc48th275stp((*this),inp,nuc,factor,substep);
  }else if(substep>275&&substep<=1725){
    mult.MultiStepCalc48th1725stp((*this),inp,nuc,factor,substep);
  }else if(substep>1725&&substep<=10787){
    mult.MultiStepCalc48th10787stp((*this),inp,nuc,factor,substep);
  }else{
    cout<<"In MultiStepCalc: substep over 10787\nInstead of that, subsubstep matrix exponential is calculated by CRAM 14th\n";
    
    real dt=factor/substep;
    GroupData2D mexpt=(*this).CalMatrixExponentialByChebyshev14(dt);
    mexpt.ReducedForm();
    nuc[0]=mexpt*inp;
    for(int k=1;k<substep;k++){
      nuc[k]=mexpt*nuc[k-1];
    };
  };
};

GroupData2D GroupData2D::MakeCorrelationMatrix()
{
  if(reduced_form){
    cout<<"Not coded in MakeCorrelationMatrix.\n";
    exit(0);
  };

  if(x!=y){
    cout<<"This matrix is not square\n";
    exit(0);
  };
  GroupData2D corr;
  corr.put_yx(y,x);
  corr.set_zero();

  real tmp;
  real sdi;
  real sdj;
  for(int i=0;i<y;i++){
    sdi=sqrt((*this).get_dat(i,i));
    for(int j=0;j<x;j++){
      sdj=sqrt((*this).get_dat(j,j));
      tmp=((*this).get_dat(i,j))/sdi/sdj;
      corr.put_data(i,j,tmp);
    };
  };
  
  return corr;
};
