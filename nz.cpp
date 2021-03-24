#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List calcGamma(arma::vec nu, arma::mat z, arma::mat b, arma::mat c,
                     arma::mat J, arma::mat J0, arma::mat K, arma::mat L, double gamma0, 
                     double DET_QVP, double DET_QVP0, arma::mat a, arma::vec mu, int d, arma::vec N,
                     arma::mat bJL, arma::mat bJL0, arma::mat cKL, arma::mat G, double dt,
                     double half_sum_elem_prod_P_a, double r0, double half_sigmaE2, 
                     double sigmaE, double rho, int niter) {
  
  int ncols = z.n_cols;
  arma::mat gamma(ncols,ncols);
  arma::mat del_m_i(d, ncols);
  arma::mat dz(d, ncols);
  arma::vec r_z(ncols);
  arma::vec m_z(ncols);
  
  for (int it=0; it<niter; it++){
    
    for (int i=0; i<ncols; i++){
      
      arma::colvec delta_i = z.col(i)-nu;
      arma::mat del_gamma_N(d,ncols);
      del_gamma_N.zeros();
      for (int j=0; j<ncols; j++){
        
        arma::colvec delta_ij = z.col(i)-z.col(j);
        if (j != i){
          gamma(i,j) = arma::as_scalar(
            gamma0 * DET_QVP *
              exp(
                (-0.5*(delta_i.t() * b * delta_i + delta_ij.t() * c * delta_ij)) +
                  (0.5*(delta_i.t() * J * delta_i + 0.5*(delta_ij.t() * K * delta_ij) + 
                  delta_i.t() * L * delta_ij))
              )
          );
          del_gamma_N.col(i)  +=  (-gamma(i,j)*(bJL*delta_i + cKL*delta_ij))*N(j);
          
        }else{
          gamma(i,j) = arma::as_scalar(
            gamma0 * DET_QVP0 *
              exp(
                (-0.5*(delta_i.t() * b * delta_i )) +
                  (0.5*(delta_i.t() * J0 * delta_i))
              )
          );
          del_gamma_N.col(i)  +=  (-gamma(i,j)*(bJL0*delta_i))*N(j);
        }
        
      }//j
      
      del_m_i.col(i) = -a*(z.col(i)-mu) - del_gamma_N.col(i);
      dz.col(i) = G*del_m_i.col(i)*dt;
      
      r_z(i) = as_scalar(r0 - half_sum_elem_prod_P_a - 0.5*(z.col(i)-mu).t()*a*((z.col(i)-mu)));
      m_z(i) = r_z(i) - half_sigmaE2 - accu(gamma.row(i)*N); 
      
    }//i
    
    z += dz;
    arma::vec dB = rep(rnorm(1, 0, sqrt(dt)),ncols);
    arma::vec dBi = rnorm(ncols, 0, sqrt(dt));
    
    N = exp(log(N) + (m_z*dt + sigmaE*sqrt(rho)*dB + sigmaE*sqrt(1-rho)*dBi));
  }
  
  return Rcpp::List::create(Rcpp::Named("N") = N,
                            Rcpp::Named("z") = z
  );
  
}
