#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppThread.h>


using namespace Rcpp;
using namespace RcppParallel;
using namespace std;



double lambda(arma::vec& pik,arma::vec& u,double EPS=0.0000001){
  // arma::mat kern = arma::null(B);
  // arma::uword N = pik.size();
  // arma::vec u(N);
  // u = kern.col(0);
  // int ncol = kern.n_cols;
  double l1 = 1e+200;
  double l2 = 1e+200;
  double l = 1e-9;
  arma::uword N = pik.size();

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

  return(l);
}




arma::vec unull(const arma::mat& A,arma::vec& pik){

  double EPS = 0.00001;
  size_t J = A.n_rows;
  size_t N = A.n_cols;

  // find non integervalue
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
  arma::mat B = A.cols(i);

  // define u
  arma::vec u(N);
  u.fill(0.0);

  // calculates nullspace of B and check if empty
  arma::mat kern = null(B);
  if(kern.empty()){
    return u;
  }

  u.elem(i) = kern.col(0);
  return(u);
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
//'
//' @export
// [[Rcpp::export]]
void tinyFF(const arma::mat& X,
                       arma::vec& pik,
                       double EPS=0.0000001){

  size_t N = X.n_rows;

  const arma::mat D = arma::diagmat(1/pik);
  const arma::mat A = (D*X).t();


  arma::vec tmp(N);


  for(std::size_t begin = 0; begin < N; begin++ ){
    arma::mat A_tmp = A.cols(0,N-1);
    // std::cout << A << endl;
    // std::cout << N << endl;

    arma::vec u = unull(A,pik);
    double l = lambda(pik,u);

    tmp.fill(l);

    pik += tmp%u;
    pik.clamp(0.0,1.0);


  }
}




/*** R


rm(list = ls())
N = 100
n = 10
p = 5
pik= sampling::inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik

tinyFF(X,pik)
test <- parallelffphase(X,pik)
test
t(A)%*%test
t(A)%*%pik
*/


