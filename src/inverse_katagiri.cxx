void inverse(vector< vector<double> >&B1,vector< vector<double> >&B2,int N){
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
}
