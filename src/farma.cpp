#include <RcppArmadillo.h>
#include <mutex>
#include <thread>
#include <stdlib.h> 

#include <RcppThread.h>
// std::mutex m2;//you can use std::lock_guard if you want to be exception safe


class Block {
private:
  arma::mat B;
  arma::vec pik;
  std::mutex m2;
public:
  
  // CONSTRUCTOR 
  Block() {}
  Block(arma::mat B,arma::vec pik) : B(B), pik(pik) {}
  Block(const Block & obj) {
    B = obj.B;
    pik = obj.pik;
  }
  
  //GET and SET
  arma::mat getB(){
    return(B); 
  }
  arma::vec getpik(){
    return(pik);
  }
  
  void setB(arma::mat newB){
    B = newB;
  }
  void setpik(arma::vec newpik){
    pik = newpik;
  }
  
  
  //METHODS
  void onestepf(){
    
    std::lock_guard<std::mutex> lockGuard(m2);
    
    double EPS=0.0000001;
    arma::mat kern = arma::null(B);
    arma::uword N = pik.size();
    arma::vec u(N);
    u = kern.col(0);
    // int ncol = kern.n_cols;
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
    // m2.unlock();
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
arma::vec farma(const arma::mat& X,
                arma::vec pik,
                unsigned int nthreads,
                double EPS=0.0000001){
  
  // arma::mat D = arma::diagmat(1/pik);
  // arma::mat A = D*X;
  
  arma::vec pikInit = pik;
  unsigned int J = X.n_cols; // initial size of B
  // unsigned int nthreads = 3;
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), nthreads*(J+1), "first"); // find first index of B
  // arma::mat B = (A.rows(i)).t(); // extract B of A
  
  
  while(i.size() > 0){
    arma::uvec p = arma::find(pik > EPS && pik < (1-EPS)); // find first index of B
    RcppThread::Rcout << p.size() << std::endl;
    // while(step < 1){
    
    arma::mat A_tmp(X.n_cols,J+1);
    
    if(i.size() <  nthreads*(J+1) && nthreads > 1){
      nthreads--;
    }
    
    if(i.size() < (J+1)){
      
      // if(i.size() < nthreads*(J+1) ){
      //   
      arma::uvec i_tmp;
      //   if(i.size() < (J+1)){
      A_tmp.resize(X.n_cols, i.size());
      i_tmp = i;
      // }else{
      //   
      //   // std::cout << J << std::endl;
      //   // std::cout << i.size() << std::endl;
      //   
      //   i_tmp = i.subvec(0,J);
      // }
      // 
      
      
      for(int r = 0;r < X.n_cols; r++){
        for(int s = 0; s < i_tmp.size(); s++){
          A_tmp(r,s) = X(i_tmp[s],r)/pikInit[i_tmp[s]];
        }
      }
      
      
      arma::mat kern = arma::null(A_tmp);
      if(kern.empty()){
        break;
      }
      
      arma::vec pik_tmp = pik.elem(i_tmp);
      Block block(A_tmp,pik_tmp);
      
      // std::cout << A_tmp << std::endl;
      
      // arma::vec pik_tmp = pik.elem(i.subvec(0,J));
      // Block block((A.rows(i.subvec(0,J))).t(),pik_tmp);
      
      std::thread threadObj(&Block::onestepf,&block);
      threadObj.join();
      
      // block.onestepf();
      // pik.elem(i.subvec(0,J)) = block.getpik();
      
      pik.elem(i_tmp) = block.getpik();
      // onestepf(B,pik_tmp);
      // onestepf((A.rows(i.subvec(0,J))).t(),pik_tmp);
      // pik.elem(i.subvec(0,J)) = pik_tmp;
    }else{
      
      
      std::vector<arma::uvec> i_tmp;
      std::vector<std::thread> threads;
      std::vector<Block> B;
      
      B.reserve(nthreads);
      i_tmp.reserve(nthreads);
      threads.reserve(nthreads);
      
      for(int k = 0; k < nthreads; ++k){
        i_tmp.push_back(i.subvec(k*(J+1),(k+1)*J + k));
      }
      for(int k = 0; k < nthreads; ++k){
        
        // for(int s = 0; s < i_tmp[k].size(); s++){
        // std::cout << i_tmp[k][s] << std::endl;
        // }
        
        // std::cout << i_tmp[k].size() << std::endl;
        // std::cout << X.n_cols <<std::endl;
        
        // arma::mat A_tmp(X.n_cols,i_tmp[k].size());
        for(int r = 0;r < X.n_cols; r++){
          for(int s = 0; s < i_tmp[k].size(); s++){
            A_tmp(r,s) = X(i_tmp[k][s],r)/pikInit[i_tmp[k][s]];
          }
        }
        
        // std::cout << A_tmp << std::endl;
        // std::cout << (A.rows(i_tmp[k])).t() <<std::endl;
        
        // B.push_back(Block((A.rows(i_tmp[k])).t(), pik.elem(i_tmp[k])));
        B.push_back(Block(A_tmp, pik.elem(i_tmp[k])));
        threads.push_back(std::thread(&Block::onestepf,&B[k]));
      }
      
      
      // std::vector<Block>::iterator it;
      // for (it = B.begin(); it != B.end(); ++it) {
      //   std::cout << (*it).getpik();
      // }
      
      
      for (std::thread & th : threads){
        // If thread Object is Joinable then Join that thread.
        if (th.joinable()){
          th.join();
        }
      }
      for(int k = 0; k < nthreads; ++k){
        pik.elem(i_tmp[k]) = B[k].getpik();
      }
      
    }
    
    
    
    i = arma::find(pik > EPS && pik < (1-EPS), nthreads*(J+1), "first");
    
    
    
  }
  return(pik);
}


/*** R
rm(list = ls())
N = 5000
n = 100
p = 50
pik= sampling::inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
system.time(s <- farma(X,pik,25))
system.time(flightphase_arma4(X,pik))
system.time(s <- BalancedSampling::flightphase(pik,X))
system.time(s <- tinyFF::flightphase(pik,X))
sum(s)
*/