void LeastSquaresMethod(int data_num, double *x, double *y, double a, double b)
{
  int n=2;
  double sum_xy=0., sum_y=0., sum_x=0., sum_x2=0.;

  for(int i=0;i<data_num;i++){
    sum_xy+=x[i]*y[i];
    sum_y+=y[i];
    sum_x+=x[i];
    sum_x2+=pow(x[i],2);
  };

  vector<vector<double>> A(n);
  for(int i=0;i<n;i++){
    A[i].resize(n);
  };
  A[0][0]=sum_x2;
  A[0][1]=sum_x;
  A[1][0]=sum_x;
  A[1][1]=data_num;

  //inverse of A
  vector<double> A_inv(n*n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A_inv[i*n+j]=A[i][j];
    };
  };
  gauss_inverse(A_inv, n); // [A_inv] is inversed.

  a=A_inv[0]*sum_xy+A_inv[1]*sum_y;
  b=A_inv[2]*sum_xy+A_inv[3]*sum_y;
  //cout<<"a = "<<a<<", b = "<<b<<"\n";
};
