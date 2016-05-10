// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>    
#include <vector>     

using namespace std;

// [[Rcpp::export]]
arma::mat algo2_cpp(arma::vec ys, arma::vec C,
                    int max_iter = 115000, int skip = 500,
                    double alpha=1, double mu_0=2, double K=100, double v_0=1, double u_0=1){
  
  int count = 1; int n = ys.n_elem;
  int max_count = max_iter/skip;
  arma::mat posterior = arma::zeros(max_count, 2);
  arma::mat phis = arma::zeros(n, 2);
    
  // initialize phi's from base measure. the number of elements in state cannot be larger than n.
  for( int i=0 ; i < n ; ++i ){ 
    phis(i,1) = R::rgamma(v_0/2, (v_0*u_0)/2 ); //tau_i
    phis(i,0) = R::rnorm(mu_0, (1/(K*phis(i,1)))); //mu_i 
  }

  arma::uvec ni, class_index;
  arma::vec uc, class_ys;
  arma::mat prob = arma::zeros(n,3); // likelihood, joint, posterior
  arma::vec prob_class;
  cout<< "Iteration start!!" <<endl;
  
  double joint, post, rnd;
  int class_ind, posterior_index;
  // Iteration
  while(max_iter >= count){
    // sample phi
    uc = unique(C); ni = hist(C, unique(C));
    for(int j=0; j<n; ++j){
      if( j < uc.n_elem){
        // sampling from conditional distribution
        class_index = find(C == (j+1)); class_ys = ys(class_index);
        phis(j,1) = R::rgamma( (v_0/2+ni[j]/2), (v_0*u_0/2 + (sum(pow(class_ys - mean(class_ys),2.0)))/2 + (ni[j]*K/(ni[j]+K))*pow((mean(class_ys)-mu_0),2.0)/2) ); //tau_i
        phis(j,0) = R::rnorm( (mean(class_ys)*(ni[j]/(ni[j]+K)) + K*mu_0/(ni[j]+K)), (1/(ni[j]*phis(j,1)+K*phis(j,1))) );  //mu_i 
      }else{
        // sampling from base measure, so we can pass this part.
      }
    }
    
    // sample C
    for(int j=0; j < n; ++j){
      // setting n_i 
      uc = unique(C); ni =  hist(C, uc);
      ni[(C[j]-1)] -= 1;
      
      // calculate probability 
      prob_class = arma::zeros(uc.n_elem+1);
      for( int k =0; k < uc.n_elem; ++k ){
        prob_class(k) = ni(k)*R::dnorm(ys[j], phis(k,0), (1/phis(k,1)), false); //likelihood  
      }
      joint = R::dgamma(phis(j,1), v_0/2, (v_0*u_0)/2, false )*R::dnorm(phis(j,0), mu_0, (1/(K*phis(j,1))), false)*prob(j,0); // joint
      post = R::dgamma(phis(j,1), (v_0/2+1/2), (v_0*u_0/2 + (K/(1+K))*pow((ys[j]-mu_0),2.0)/2), false )*R::dnorm(phis(j,0), (ys[j]/(1+K)+K*mu_0/(1+K)), (1/(phis(j,1)+K*phis(j,1))), false); // posterior
      prob_class[uc.n_elem] = joint/post;
      prob_class = prob_class/sum(prob_class);
      
      rnd = Rcpp::runif(1)[0];
      arma::vec cum_class = cumsum(prob_class);
      class_ind = 1;
      for( int k=0; k < uc.n_elem; ++k ){
        if( rnd < cum_class[k]){
          break;
        }else{
          class_ind += 1;
        }
      }  
      if ((ni[(C[j]-1)] == 0) | (class_ind == uc.n_elem)) {
        //removce phi_ci
        phis((C[j]-1),1) = R::rgamma( (v_0/2+1/2), (v_0*u_0/2 + (K/(1+K))*pow((ys[(C[j]-1)]-mu_0),2.0)/2) ); //tau_i
        phis((C[j]-1),0) = R::rnorm((ys[(C[j]-1)]/(1+K)+K*mu_0/(1+K)), (1/(phis((C[j]-1),1)+K*phis((C[j]-1),1)))); //mu_i 
      }
    }
    if( count % skip == 0 ){
      //save
      cout << count << endl;
      posterior_index = count/skip-1;
      posterior(posterior_index, 0)= phis(0,0);
      posterior(posterior_index, 1)= phis(0,1);
    }
    count += 1;
  }
  
  return posterior;
}
