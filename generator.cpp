// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat stick_breaking_cpp( double alpha=1, double mu_0=2, double K=0.01, double v_0=2, double u_0=1, double tol=0.01){
  arma::vec pis = arma::zeros(1);  arma::vec betas = arma::zeros(1);
  double new_beta = R::rbeta(1,alpha);
  betas.tail(1) = new_beta;  pis.tail(1) = new_beta;
  double s=0;
  while ( sum(pis) < (1.0-tol) ){
    s = prod(1-betas);
    new_beta = R::rbeta(1,alpha);
    betas.resize( betas.n_elem+1 );
    betas.tail(1) += new_beta;
    pis.resize( pis.n_elem+1 );
    pis.tail(1) += new_beta*s;
  }
  arma::mat theta = arma::zeros(pis.n_elem, 3);
  double var =0; double beta = u_0*v_0/2; // beta is rate rate
  for (int i = 0; i < pis.n_elem ; ++i ){
    theta(i,2) = R::rgamma(v_0/2, 1/beta ); //tau_i, parameter: alpha, scale.
    var = 1/(K*theta(i,2)); //inverse of precision
    theta(i,1) = R::rnorm(mu_0, pow((var), 0.5) ); //mu_i, parameter: mu, sd
    theta(i,0) = pis[i]; //weight p_i
  }
  return theta;
}

// [[Rcpp::export]]
Rcpp::List generator_cpp(int num_samples=5, double alpha=1, double mu_0=2, double K=0.01, double v_0=2, double u_0=1, double tol=0.01){
  arma::vec ys = arma::zeros(num_samples);
  arma::vec index_class = arma::zeros(num_samples);
  arma::mat theta = stick_breaking_cpp(alpha, mu_0, K, v_0, u_0, tol); 
  arma::vec pis = theta.col(0);
  arma::vec cum_pis = cumsum(pis);
  double rnd;  int ind; double var;
  for(int j=0; j<num_samples; ++j){
    ind = 0;    rnd = Rcpp::runif(1)[0];
    for(int i=0; i < pis.n_elem; ++i){
      if( rnd < cum_pis[i]){
        break;
      }else{
        ind += 1;
      }
    }
    index_class[j] = ind;
    var = 1/theta(ind,2);
    ys[j] = R::rnorm(theta(ind,1), pow((var), 0.5));
  }
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("ys") = ys, Rcpp::Named("class") = index_class);
}