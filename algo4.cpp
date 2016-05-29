// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>    
#include <vector>     

using namespace std;

// [[Rcpp::export]]
arma::mat algo4_cpp(arma::vec ys, arma::vec C,
                    int max_iter = 11500, int thining = 100, int burn_in = 1500,
                    double alpha=1.0, double mu_0=2.0, double K=2.0, double v_0=0.2, double u_0=1.0){
  int count = 1; int n = ys.n_elem;
  int max_count = (max_iter-burn_in)/thining;
  double beta = (v_0*u_0)/2;
  arma::mat posterior = arma::zeros(max_count, 2*n); 
  arma::mat phis = arma::zeros(n, 2); // state matrix
  
  // initialize phi's from base measure. the number of elements in state cannot be larger than n.
  double var=0;
  for( int i=0 ; i < n ; ++i ){ 
    beta = (v_0*u_0)/2;
    phis(i,1) = R::rgamma(v_0/2, 1/beta ); //tau_i
    var = 1/(K*phis(i,1));
    phis(i,0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i 
  }
  
  arma::uvec ni, class_index;  arma::vec uc, class_ys;
  arma::mat prob = arma::zeros(n,3); // likelihood, joint, posterior
  arma::vec prob_class;
  cout<< "Iteration start!!" <<endl;
  
  double rnd;
  int class_ind, posterior_index; int zero_class = 0;
  bool is_update=true; bool is_zero=false; double thre=0;
  // Iteration
  while(true){
    if(max_iter < count){
      break;
    }
    
    // update C
    for(int j=0; j < n; ++j){
      is_update=true; is_zero=false;
      // setting n_i 
      uc = unique(C); ni =  hist(C, uc); // order by class
      ni[(C[j]-1)] -= 1;
      
      // if c_j is associated with no other obsevation then initialize
      if(ni[(C[j]-1)] == 0){
        is_zero=true;
        zero_class = C[j];
        rnd = Rcpp::runif(1)[0];
        thre = uc.n_elem / (uc.n_elem+1);
        if ( rnd < thre){
          is_update = false;
        }
      }  
      
      if(is_update ==true){
        if(is_zero==true){
          phis((uc.n_elem),1)=phis((C[j]-1),1);
          phis((uc.n_elem),0)=phis((C[j]-1),0);
        }else{
          beta = (v_0*u_0)/2;
          phis((uc.n_elem),1) = R::rgamma(v_0/2, 1/beta ); //tau_i
          var = 1/(K*phis((uc.n_elem),1));
          phis((uc.n_elem),0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i 
        }
        
        // sampling probability by equation (3.6) in NEAL's paper
        // calculate probability for existing class
        prob_class = arma::zeros(uc.n_elem+1);
        for( int k = 0; k < uc.n_elem; ++k ){
          var = 1/phis(k,1);
          prob_class(k) = ni(k)*R::dnorm(ys[j], phis(k,0), pow((var),0.5), false); //likelihood  
        }
        //for new class
        var = 1/phis((uc.n_elem),1);
        prob_class[uc.n_elem] = alpha*R::dnorm(ys[j], phis((uc.n_elem),0), pow((var),0.5), false)/(uc.n_elem+1);
        prob_class = prob_class/sum(prob_class); 
        
        // class sampling
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
        // class update
        C[j] = class_ind;
        
        // standardize C
        if(is_zero==true){
          class_index = find(C > zero_class); C(class_index) -= 1;
          //remove
          phis.shed_row((zero_class-1));
          //insert
          phis.insert_rows((n-1),1);
          beta = (v_0*u_0)/2;
          phis((n-1),1) = R::rgamma(v_0/2, 1/beta ); //tau_i
          var = 1/(K*phis((n-1),1));
          phis((n-1),0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i
        }
      }
    }
    
    // update phi
    uc = unique(C); ni = hist(C, unique(C));
    for(int j=0; j<n; ++j){
      if( j < uc.n_elem){
        // sampling from conditional distribution
        class_index = find(C == (j+1)); class_ys = ys(class_index);
        beta = (v_0*u_0/2 + (sum(pow(class_ys - mean(class_ys),2.0)))/2 + (ni[j]*K/(ni[j]+K))*pow((mean(class_ys)-mu_0),2.0)/2);
        phis(j,1) = R::rgamma( (v_0/2+ni[j]/2), 1/beta ); //tau_i
        var =  (1/(ni[j]*phis(j,1)+K*phis(j,1)));
        phis(j,0) = R::rnorm( (mean(class_ys)*(ni[j]/(ni[j]+K)) + K*mu_0/(ni[j]+K)), pow((var),0.5));  //mu_i 
      }else{
        // sampling from base measure
        beta = (v_0*u_0)/2;
        phis(j,1) = R::rgamma(v_0/2, 1/beta ); //tau_i
        var = 1/(K*phis(j,1));
        phis(j,0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i 
      }
    }
    
    if((count % thining == 0) & (count > burn_in)){
      //save
      cout << count << endl;
      posterior_index = (count-burn_in)/thining-1;
      for( int k=0; k < n; ++k ){
        posterior(posterior_index, k)= phis((C[k]-1),0); // mu_i
        posterior(posterior_index, n+k)= phis((C[k]-1),1); // tau_i
      }
    }
    count += 1;
  }
  
  return posterior;
}












