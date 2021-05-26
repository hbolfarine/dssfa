#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Generate Random Multivariate Normal Distribution
//' @description Creates a random Multivariate Normal when given number of obs, mean, and sigma. 
//' @param n An \code{int}, which gives the number of observations.  (> 0)
//' @param mu A \code{vector} length m that represents the means of the normals.
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is a Multivariate Normal distribution
//' @seealso \code{\link{TwoPLChoicemcmc}} and \code{\link{probitHLM}}
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' rmvnorm(2, c(0,0), diag(2))
//' 
// [[Rcpp::export]]
arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S) {
  unsigned int ncols = S.n_cols;
  arma::mat Y(n, ncols);
  Y.imbue( norm_rand ) ;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(S);
}





// [[Rcpp::export]]
double rTruncNorm_lb(double mean,double sd, double b_lb){
  double p0 = R::pnorm(b_lb,mean,sd, 1, 0);
  double p1 = 1-p0;
  double uZ = R::runif(0,1);
  double Z = R::qnorm(p0 + uZ*p1, mean, sd, 1, 0);
  return(Z);
}


// [[Rcpp::export]]
void update_invClam(const arma::mat& Lambda,
                    arma::mat& invClam) {
  
  // No spike-and-slab, only column-wise variance
  
  unsigned int J = Lambda.n_rows;
  int M = Lambda.n_cols;
  double u = R::runif(0.0,1.0);
  double d = accu(pow(Lambda, 2));
  invClam.fill(R::qgamma(u, 
                         double(J*M-((M-1)*M/2) + 2)/2.0, 
                         2.0/(2.0 + d),1,0));
}


// [[Rcpp::export]]
void update_EFA_parms(unsigned int N,unsigned int J,unsigned int K,const arma::mat& Y,
                      const arma::mat& I_K,arma::mat& F,arma::mat& Lambda,arma::vec& psis_inv,
                      const arma::mat invClam){
  //update F
  arma::mat Psi_inv = arma::diagmat(psis_inv);
  arma::mat Psi_inv_Lam = Psi_inv*Lambda;
  arma::mat VF = inv(I_K + Lambda.t()*Psi_inv_Lam);
  arma::mat mF = Y*Psi_inv_Lam*VF;
  F = mF + arma::randn<arma::mat>(N,K)*arma::chol(VF); 
  
  //update Lambda + psis_inv
  arma::mat FpF = F.t()*F;
  arma::mat FpY = F.t()*Y;
  arma::vec Y_m_Flam(K);
  arma::rowvec lambdaj(K);  
  for(unsigned int j=0;j<J;++j){
    
    arma::mat Cj = inv(psis_inv(j)*FpF+invClam);
    arma::vec mj = psis_inv(j)*Cj*FpY.col(j);
    lambdaj = rmvnorm(1,mj,Cj);
    Lambda.row(j) = lambdaj;
      
    arma::rowvec lambdaj = Lambda.row(j);
    Y_m_Flam = Y.col(j)-F*lambdaj.t();
    double u = R::runif(0.0,1.0);
    
    //qgamma is parameterized as shape, scale rather than shape and rate
    psis_inv(j) = R::qgamma(u,double(N+2.)/2.0,
             2.0/(2.+arma::dot(Y_m_Flam,Y_m_Flam)),1,0);
  }
}


// [[Rcpp::export]]
Rcpp::List EFA_unconstrained(const arma::mat& Y,unsigned int K,double l0,double l1,
                       unsigned int burnin,unsigned int chain_length=10000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  unsigned int chain_m_burn = chain_length-burnin;
  unsigned int tmburn;
  arma::mat I_K = arma::eye(K,K);
  
  //Saving output
  arma::cube LAMBDA(J,K,chain_m_burn);
  arma::mat PSIs(J,chain_m_burn);
  arma::vec DELTA=arma::zeros<arma::vec>(K);
  
  // Initialize intial values
  arma::mat F = arma::randn<arma::mat>(N,K);
  arma::mat Cstart=inv(I_K+F.t()*F);
  arma::mat Lambda = (arma::randn<arma::mat>(J,K))*arma::chol(Cstart);
  for(unsigned int j=0;j<K;++j){
    Lambda(j,j) = rTruncNorm_lb(0.,1.,0.);
    for(unsigned int k=j+1;k<K;++k){
      Lambda(j,k) = 0.;
    }
  }
  arma::vec psis_inv(J);
  arma::rowvec varY=arma::var(Y);
  for(unsigned int j=0;j<J;++j){
    arma::rowvec lambdaj = Lambda.row(j);
    arma::vec Y_m_Flam = Y.col(j)-F*lambdaj.t();
    double u = R::runif(0.0,1.0);
    //qgamma is parameterized as shape, scale rather than shape and rate
    psis_inv(j) = R::qgamma(u,double(N+2.)/2.0,
             2.0/(2.+arma::dot(Y_m_Flam,Y_m_Flam)),1,0);
  }
  double omega=.5;
  arma::vec delta=arma::ones<arma::vec>(K);
  arma::mat invClam=l1*I_K;

  //Start Markov chain
  for(unsigned int t = 0; t < chain_length; ++t){
    
    update_invClam(Lambda, invClam);
    
    update_EFA_parms(N,J,K,Y,I_K,F,Lambda,psis_inv,invClam);
    if(t>burnin-1){
      tmburn = t-burnin;
      LAMBDA.slice(tmburn) = Lambda;
      PSIs.col(tmburn) = 1./psis_inv;
      DELTA+=delta;
    }
  }
  return Rcpp::List::create(Rcpp::Named("LAMBDA",LAMBDA),
                            Rcpp::Named("PSIs",PSIs),
                            Rcpp::Named("DS",DELTA/chain_m_burn)
  );
}







