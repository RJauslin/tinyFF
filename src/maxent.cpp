#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix qfromw(NumericVector& w, int& n){
  
  int N = w.size();
  NumericMatrix expa(N,n);
  double tmp = 0.0;
  
  for(int i = 0; i < N; i++){
    tmp = 0.0;
    for(int j = i; j < N; j++){
      tmp = tmp +  w[j];
    }
    expa(i,0) = tmp;
  }
  for(int i = N - n  ; i < N; i++){
    tmp = 0.0;
    for(int j = i-1; j < N; j++){
      tmp = tmp +  log(w[j]);
    }
    expa(i-1,N-i) = exp(tmp);
  }
  
  for(int i = N-2;i >= 0; i--){
    int m = std::min(N-i,n);
    for(int z = 1; z < m; z++){
      expa(i,z) = w[i]*expa(i+1,z-1) + expa(i+1,z);
    }
  }
  
  NumericMatrix q(N,n);
  
  for(int i = 0; i < N; i++){
    q(i,0) = w[i]/expa(i,0);
  }
  for(int i = N - n  ; i < N; i++){
    q(i-1,N-i) = 1.0;
  }
  for(int i = N-2;i >= 0; i--){
    int m = std::min(N-i,n);
    for(int z = 1; z < m; z++){
      q(i,z) = w[i]*expa(i+1,z-1)/expa(i,z);
    }
  }
  
  return(q);
}




/*** R
pik=c(0.07,0.17,0.41,0.61,0.83,0.91)
# pik <- sampling::inclusionprobabilities(runif(20),5)

w = (pik)/(1 - pik)
qfromw(w,sum(pik))
UPMEqfromw(w, sum(pik))
*/
