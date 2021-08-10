#' @title Fast flight phase of the cube method
#'
#' @description 
#' 
#' This function computes the flight phase of the cube method proposed by Chauvet and Tillé (2006).
#'
#' @param X A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
#' @param pik A vector of inclusion probabilities.
#'
#' @details 
#' This function implements the method proposed by (Chauvet and Tillé 2006). It recursively transforms the vector of inclusion probabilities \code{pik} into a
#' sample that respects the balancing equations. The algorithm stops when the null space of the sub-matrix \eqn{B} is empty.
#' For more information see (Chauvet and Tillé 2006).
#' 
#' The function uses the function \code{\link[MASS:Null]{Null}} to find the null space of the sub-matrix \eqn{B}.
#'
#' @return Updated vector of \code{pik} that contains 0 and 1 for unit that are rejected or selected.
#'
#' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. \emph{Computational Statistics}, 21/1:53-62
#'
#'
#' @seealso \code{\link[sampling:samplecube]{fastflightphase}}, \code{\link[BalancedSampling:flightphase]{flightphase}}. 
#'
#' 
#' @export
#' @examples 
#' rm(list = ls())
#' N = 500
#' n = 10
#' p = 5
#' pik= sampling::inclusionprobabilities(runif(N),n)
#' X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
#' A=X/pik
#' 
ffphase2 <- function(X, pik){
  
  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  EPS = 1e-8
  A <- X/pik
  N <- length(pik)  
  p = ncol(X)
  
  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------
  
  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  
  ##----------------------------------------------------------------
  ##                          flight phase                         -
  ##----------------------------------------------------------------
  
  task(pik,A)
 

  simple_reduce(N-p,task,pik,A)
  
  return(pik)
}


simple_reduce <- function(n, f, pik, A) {
  out <- pik
  for (i in seq(2, n)) {
    out <- f(out, A)
  }
  out
}

task <- function(pik,A){
  
  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  
  ##------ Find B
  i_tmp <- i[1:(p+1)]
  pik_tmp <- pik[i_tmp]

  B <- as.matrix(A[i_tmp,])
  
  ##------ onestep and check if null
  tmp <-  onestep(B,pik_tmp,EPS)
  pik[i_tmp] <- tmp
    
  return(pik)
}



#' Internal function of ffphase
#' @noRd
onestep <- function(B,pik,EPS){
  
  kern <- MASS::Null(B)
  if(length(kern) == 0){
    return(NULL)
  }
  N <- length(pik)
  u = kern[,1]
  
  l1=min(pmax((1-pik)/u,-pik/u))
  l2=min(pmax((pik-1)/u,pik/u))
  
  if(stats::runif(1) < l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  pik = pik + l*u
  
  return(pik);
}


