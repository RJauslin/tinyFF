#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <mutex>
#include <thread>

#include <stdlib.h> 

// std::mutex m2;//you can use std::lock_guard if you want to be exception safe

class Pikstar {
private:
  
  arma::vec pikstar;
  std::mutex m3;
  
public:
  
  // CONSTRUCTOR 
  Pikstar() {}
  Pikstar(arma::vec pikstar) : pikstar(pikstar) {}
  
  //GET and SET
  arma::vec getpik(){
    return(pikstar);
  }
  void setpik(arma::vec newpik){
    pikstar = newpik;
  }
  
  
  //METHODS
  void onestepff(arma::mat B,arma::uvec i_tmp){
    
    
    // m3.lock();
    std::lock_guard<std::mutex> lockGuard(m3);
    
    double EPS=0.0000001;
    arma::mat kern = arma::null(B);
    arma::uword N = i_tmp.size();
    arma::vec u(N);
    u = kern.col(0);
    // int ncol = kern.n_cols;
    double l1 = 1e+200;
    double l2 = 1e+200;
    double l = 1e-9;
    
    
    for(arma::uword k = 0; k < N; k++){
      if(u[k]> 0){
        l1 = std::min(l1,(1.0 - pikstar[i_tmp[k]])/u[k]);
        l2 = std::min(l2,pikstar[i_tmp[k]]/u[k]);
      }
      if(u[k]< 0){
        l1 = std::min(l1,-pikstar[i_tmp[k]]/u[k]);
        l2 = std::min(l2,(pikstar[i_tmp[k]]-1.0)/u[k]);
      }
    }
    if(Rcpp::runif(1)[0]<l2/(l1+l2)){
      l = l1;
    }else{
      l = -l2;
    }
    for(arma::uword k = 0; k < N; k++){
      pikstar[i_tmp[k]] = pikstar[i_tmp[k]] + l*u[k];
      if(pikstar[i_tmp[k]] < EPS){
        pikstar[i_tmp[k]] = 0;
      }
      if(pikstar[i_tmp[k]] > (1-EPS)){
        pikstar[i_tmp[k]] = 1;
      }
    }
    
    // m3.unlock();
    
    
    return;
  }
  
  
};




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
//' @examples
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector parma(const arma::mat& X,
                arma::vec pik,
                int nthreads,
                double EPS=0.0000001){
  
  arma::vec pikInit = pik;
  unsigned int J = X.n_cols; 
  
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), nthreads*(J+1), "first"); // find first index of B
  
  
  Pikstar pikstar(pik);
 
  // int step = 0;
  // while(step < 30){
  while(i.size() > 0){

    arma::uvec p = arma::find(pikstar.getpik() > EPS && pikstar.getpik() < (1-EPS)); // find first index of B
    Rcpp::Rcout << p.size() << "\n";

    arma::mat A_tmp(X.n_cols,J+1);

    if(i.size() <  nthreads*(J+1) && nthreads > 1){
      nthreads--;
    }
    

    if(i.size() < nthreads*(J+1)){

      // std::cout << "end" << std::endl;

      arma::uvec i_tmp;

      A_tmp.resize(X.n_cols, i.size());
      i_tmp = i;


      for(int r = 0;r < X.n_cols; r++){
        for(int s = 0; s < i_tmp.size(); s++){
          A_tmp(r,s) = X(i_tmp[s],r)/pikInit[i_tmp[s]];
        }
      }


      arma::mat kern = arma::null(A_tmp);
      if(kern.empty()){
        break;
      }

      std::thread th = std::thread(&Pikstar::onestepff, &pikstar, A_tmp, i_tmp);
      th.join();


    }else{

      
      std::vector<std::thread> threads;
      threads.reserve(nthreads);
      std::vector<arma::uvec> is;
      
      
      for(int k = 0; k < nthreads; ++k){
        is.push_back(i.subvec(k*(J+1),(k+1)*J + k));
      }
  

      for(int k = 0; k < nthreads; ++k){
        arma::mat A_tmp(X.n_cols,J+1);
        for(int r = 0;r < X.n_cols; r++){
          for(int s = 0; s < is[k].size(); s++){
              A_tmp(r,s) = X(is[k][s],r)/pikInit[is[k][s]];
          }
        }
       threads.push_back(std::thread(&Pikstar::onestepff,&pikstar,A_tmp,is[k]));
        
      }
      for (std::thread & th : threads){
        if (th.joinable()){
          th.join();
        }
      }
      
      


    }


    i = arma::find(pikstar.getpik() > EPS && pikstar.getpik() < (1-EPS), nthreads*(J+1), "first");
    // step = step + 1;
  }
  
  return(Rcpp::wrap(pikstar.getpik()));
}


/*** R
rm(list = ls())
N = 5000
n = 100
p = 50
pik= sampling::inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
system.time(s <- parma(X,pik,4))
system.time(s <- flightphase_arma4(X,pik))
system.time(s <- BalancedSampling::flightphase(pik,X))
system.time(s <- tinyFF::flightphase(pik,X))
sum(s)
*/