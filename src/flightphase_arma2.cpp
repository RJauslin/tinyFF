#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;


double lambda_2(arma::vec& pik,arma::vec& u,double EPS=0.0000001){
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




arma::vec unull_2(const arma::mat& A,arma::vec& pik){

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
arma::vec flightphase_arma_2(const arma::mat& X,arma::vec pik,double EPS=0.0000001){

  size_t N = X.n_rows;

  const arma::mat D = arma::diagmat(1/pik);
  const arma::mat A = (D*X).t();

  arma::vec tmp(N);
  for(std::size_t begin = 0; begin < N; begin++ ){
    arma::mat A_tmp = A.cols(0,N-1);
    // std::cout << A << endl;
    // std::cout << N << endl;

    arma::vec u = unull_2(A,pik);
    double l = lambda_2(pik,u);

    tmp.fill(l);

    pik += tmp%u;
    pik.clamp(0.0,1.0);


  }

  return(pik);
}



struct Task : public Worker
{
  // source
  const arma::mat A;
  const arma::vec pik;
  size_t N;
  
  // output
  arma::vec out;
  
  
  Task(const arma::mat A,const arma::vec pik,size_t N,arma::vec out) : A(A), pik(pik), N(N), out(out) {}
  Task(const Task& task, Split) : A(task.A), pik(task.pik), N(task.N), out(task.out) {}
  
  
  void operator()(std::size_t begin, std::size_t end) {
    
    arma::mat A_tmp = A.cols(begin,end);
    arma::vec pik_tmp = pik.subvec(begin,end);
    // arma::vec pik_tmp = pik.subvec(begin,end);
    
    
    // arma::mat A_tmp = A;
    // arma::vec pik_tmp = pik;
    
    arma::vec u = unull_2(A_tmp,pik_tmp);
    double l = lambda_2(pik_tmp,u);
    arma::vec tmp(N);
    tmp.fill(l);
    out += tmp%u;
    // out.clamp(0.0,1.0);
  }
  
  // join my value with that of another Sum
  void join(const Task& rhs) {
    // value += rhs.value;
    out += rhs.out;
    out.clamp(0.0,1.0);
  }
  
};



//
//' @title blab
//'
//' @description skjdfng
//'
//' @return A vector
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' @export
// [[Rcpp::export]]
arma::vec parallelffphase_2(const arma::mat X,const arma::vec pik) {
  // NumericVector parallelffphase(const arma::mat& X,const arma::vec& pik) {
  
  const arma::mat D = arma::diagmat(1/pik);
  const arma::mat A = (D*X).t();
  size_t N = X.n_rows;
  arma::vec s(pik);

  Task task(A,pik,N,s);
  
  parallelReduce(0, 10, task);
  
  return task.out;
}




// struct FP : public Worker
// {
//   // source
//   const arma::mat A;
//   size_t N;
// 
//   //output
//   arma::vec out;
// 
//   FP(const arma::mat A, size_t N,arma::vec out) : A(A), N(N), out(out) {}
//   // FP(const FP& fp, Split) : A(fp.A), N(fp.N), out(arma::vec(N)) { out.fill(0);}
// 
// 
//   void operator()(std::size_t begin, std::size_t end) {
//     arma::mat A_tmp = A.cols(begin,end);
//     arma::vec out_tmp = out.subvec(begin,end);
//     for(std::size_t begin = begin; begin < end; begin++ ){
// 
//       arma::vec u = unull_2(A_tmp,out_tmp);
//       double l = lambda_2(out_tmp,u);
//       arma::vec tmp(N);
//       tmp.fill(l);
// 
//       out += tmp%u;
//       out.clamp(0.0,1.0);
// 
//     }
//   }
// 
//   // join my value with that of another Sum
//   // void join(const FP& fp) {
//   //   out += fp.out;
//   // }
// 
// };
// 
// 
// //
// //' @title blab
// //'
// //' @description skjdfng
// //'
// //' @return A vector
// //'
// //' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
// //' @export
// // [[Rcpp::export]]
// arma::vec parallelffphase_2(const arma::mat X,const arma::vec pik) {
//   // NumericVector parallelffphase(const arma::mat& X,const arma::vec& pik) {
// 
//     const arma::mat D = arma::diagmat(1/pik);
//     const arma::mat A = (D*X).t();
//     size_t N = X.n_rows;
// 
// 
//     arma::vec s(pik);
//     // arma::vec s(pik);
//     FP fp(A,N,s);
// 
// 
//     parallelFor(0, N, fp);
// 
//     return s;
//   }




/*** R


rm(list = ls())
N = 100
n = 10
p = 5
pik= sampling::inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik

test1 <- flightphase_arma_2(X,pik)
sum(test1)
test <- parallelffphase_2(X,pik)
test
sum(test)
t(A)%*%test
t(A)%*%pik
*/


