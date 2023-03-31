void multi(vector< vector<double> >&A1,vector< vector<double> >&A2,vector< vector<double> >&A3,int N){
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      A3[i][j]=0;
    }
  }
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      if(A1[i][j]!=0){
	double temp=A1[i][j];
	for(int k=0;k<N;k++){
	  A3[i][k]+=temp*A2[j][k];
	}
      }
    }
  }  
}
