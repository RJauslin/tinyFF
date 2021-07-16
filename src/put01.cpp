#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param x x
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//'
//'
//' @export
// [[Rcpp::export]]
void onestep3(arma::vec& pik,arma::vec u,double EPS=0.0000001){
  
  arma::uword N = pik.size(); 
  double l1 = 1e+200;
  double l2 = 1e+200;
  double l = 1e-9;
  
  for(arma::uword k = 0; k < N; k++){
    if(u[k]> 0){
      l1 = std::min(l1,(1.0 - pik[k])/u[k]);
      l2 = std::min(l2,pik[k]/u[k]);
    }
    if(u[k]< 0){
      l1 = std::min(l1,-pik[k]/u[k]);
      l2 = std::min(l2,(pik[k]-1.0)/u[k]);
    }
  }
  if(Rcpp::runif(1)[0]<l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  for(arma::uword k = 0; k < N; k++){
    pik[k] = pik[k] + l*u[k];
    if(pik[k] < EPS){
      pik[k] = 0;
    }
    if(pik[k] > (1-EPS)){
      pik[k] = 1;
    }
  }
  // return(pik);
}


//' @title blab
//'
//' @description skjdfng
//'
//' @return A vector
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' @export
// [[Rcpp::export]]
void put01(arma::uword& done,
           arma::uword& size,
           arma::uvec& index,
           arma::vec& pik){
  // put finished units at beginning of list
  double eps = 1e-8;
  arma::uword tempInt(0);
  // int done = 0;
  for(arma::uword i = done;i < size; i++){
    if(pik[index[i]]<eps || pik[index[i]]>1-eps){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  }
} 

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param x x
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' func
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector flightphase_arma4(Rcpp::NumericMatrix Xr,
                                      Rcpp::NumericVector pikr,
                                      double EPS=0.0000001){
// Rcpp::NumericVector flightphase_arma4(const arma::mat& X,arma::vec pik,double EPS=0.0000001){


  arma::mat X(Xr.begin(),Xr.nrow(),Xr.ncol(),false);
  arma::vec pik(pikr.begin(),pikr.size(),false);
  
  // ncol and nrow
  arma::uword N = X.n_rows;
  arma::uword J = X.n_cols;
  
  // calculate A
  // arma::mat D = arma::diagmat(1/pik);
  // arma::mat A = (D*X).t();
  arma::vec pikInit(N);
  
  
  // initiallizing 
  arma::uvec index(N);
  arma::uword done = 0;
  for(arma::uword i = 0; i < N; i++ ){
    index[i] = i;
    pikInit[i] = pik[i];
  }
  
  // change index at first to put done at the beginning and initializing end
  put01(done,N,index,pik);
  arma::uword end = done + (J+1);
  arma::uword size_sub = (J+1);
  // initializing B and kern pikstar and index_pikstar
  arma::mat B(J,J+1); 
  arma::mat kern;
  arma::vec pikstar(J+1);
  arma::uvec index_pikstar(J+1);
  
  
  while(done < N){
    
    
    // fill B, pikstar and index_pikstar
    for(arma::uword j = 0; j < size_sub; j++){
      index_pikstar[j] = index[done + j];
      pikstar[j] = pik[index[done + j]];
      for(arma::uword i = 0; i < J; i++){
        B(i,j) = X(index_pikstar[j],i)/pikInit[index_pikstar[j]];
        // B(i,j) = A(i,index_pikstar[j]);
      }
    }
    
    
    //find u
    kern = null(B);
    if(kern.empty()){
      break;
    };
    arma::vec u = kern.col(0); // vector in nullspace
    
    
    //modify pikstar
    onestep3(pikstar,u); // modify pikstar
    
    // update pik in the correct indices
    for(arma::uword j = 0;j < size_sub;j++){
      pik[index_pikstar[j]] = pikstar[j];
    }  
    
    // update done and index
    put01(done,end,index,pik);
    end = done + (J+1);
    if(end > N){
      size_sub = N-done;
      end = N;
      B.resize(J,size_sub);
      pikstar.resize(size_sub);
      index_pikstar.resize(size_sub);
    }
    
    // cout << size_sub << endl;
    
    
  }
  
  return(Rcpp::wrap(pik));
}







/*** R

rm(list = ls())
N = 10000
n = 100
p = 50
pik= sampling::inclusionprobabilities(runif(N),n)
# pik[3] <- 1
# pik[59] <- 0
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
system.time(test <- flightphase_arma4(X,pik))
system.time(test <- ffphase(X,pik))


system.time(test <- BalancedSampling::flightphase(pik,X))
system.time(test <- samptools::flightphase(pik,X))
system.time(test <- flightphase_arma4(X,pik))
sum(test)


# system.time(test <- ffphase(X,pik))
# system.time(test <- parallelffphase(X,pik)) 



*/