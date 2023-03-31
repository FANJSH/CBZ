void mk_mat_ex(int p,int q,int N,double dt,vector< vector<double> >&Mtex,vector<vector <double> >A){
  vector< vector<double> >Npq(N);
  vector< vector<double> >Dpq(N);
  for(int i=0;i<N;i++){
    Npq[i].resize(N,0.);
    Dpq[i].resize(N,0.);
  }
  //累乗計算用行列B,C (B=A^k)
  vector< vector<double> >B(N);
  for(int i=0;i<N;i++){
    B[i].resize(N);
  }
  vector< vector<double> >C(N);
  for(int i=0;i<N;i++){
    C[i].resize(N);
  }
  unit(B,N);
  //fa1=p!,fa2=q!,fa3=(p+q)!
  double fa1=mk_fact(p);
  double fa2=mk_fact(q);
  double fa3=mk_fact(p+q);
  //Calculate Npq and Dpq
  for(int k=0;k<=p;k++){
    //B=A^k
    if(k==0){
      unit(C,N);
    }else{
      multi(B,A,C,N);
    }
    B=C;
    //fa4=k!,fa5=(p+q-k)!,fa6=(p-k)! or (q-k)!
    double fa4=mk_fact(k);
    double fa5=mk_fact(p+q-k);
    double fa6=mk_fact(p-k);
    //Npqの計算
    double temp=fa5*fa1/(fa3*fa4*fa6);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	Npq[i][j]+=temp*pow(dt,k)*B[i][j];
	Dpq[i][j]+=temp*pow((-1.*dt),k)*B[i][j];
      }
    }
  }
  //Dpqの逆行列xを求める
  vector< vector<double> >x(N);
  for(int i=0;i<N;i++){
    x[i].resize(N);
  }
  inverse(Dpq,x,N); 
  multi(Npq,x,Mtex,N);
}
