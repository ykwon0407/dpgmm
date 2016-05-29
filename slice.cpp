// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>    
#include <vector>     

using namespace std;

// [[Rcpp::export]]
arma::mat slice_cpp(arma::vec ys, arma::vec C,
                    int max_iter = 115000, int thining = 1000, int burn_in = 15000,
                    double alpha=0.6, double mu_0=2.0, double K=1.0, double v_0=1.0, double u_0=1.0, double tol = 0.005){
  int count = 1; int n = ys.n_elem;
  int max_count = (max_iter-burn_in)/thining;
  double beta = (v_0*u_0)/2;
  arma::mat posterior = arma::zeros(max_count, 2*n); // result matrix
  arma::mat phis = arma::zeros(n, 2); // state matrix
  
  // initialize phi's from base measure. the number of elements in state cannot be larger than n.
  // initializer v and u
  arma::vec v_list = arma::zeros(n);
  arma::vec w_list = arma::zeros(n);
  double var=0; arma::vec u_list = arma::zeros(n);
  for( int i=0 ; i < n ; ++i ){ 
    beta = (v_0*u_0)/2;
    phis(i,1) = R::rgamma(v_0/2, 1/beta ); //tau_i
    var = 1/(K*phis(i,1));
    phis(i,0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i 
    
    v_list(i) = R::rbeta(1, alpha); // v_i
    if(C[i]==1){
      u_list[i] = Rcpp::runif(1,0, v_list[(C[i]-1)])[0]; // u_i
    }else{
      u_list[i] = Rcpp::runif(1,0, (prod(1-v_list(arma::span(0,(C[i]-2))))*v_list[(C[i]-1)]))[0]; // u_i
    }
  }
  
  arma::uvec ni, class_index;  arma::vec uc, class_ys;
  arma::vec prob_class; arma::vec cw;
  cout<< "Iteration start!!" <<endl;
  
  double rnd, v_lower, v_upper, temp, new_v;
  int class_ind, posterior_index, M; double w, thre;
  bool is_zero=false; int zero_class = 0;
  // Iteration
  while(true){
    // iteration criteria
    if(max_iter < count){
      break;
    }
    
    uc = unique(C); M=uc.n_elem;
    // update V
    for(int j=0; j < n ; ++j){
      if( j == M ){
        break;
      }
      class_index = find(C == (j+1));
      if( j == 0){
        v_lower = max(u_list(class_index)); 
      }else{
        v_lower = max(u_list(class_index))/prod(1-v_list(arma::span(0,(j-1))));
      }
      if(v_lower >= 1){
        v_lower = 1;
      }
      
      class_index = find(C > (j+1)); v_upper = 1;
      for(int k=0; k < class_index.n_elem ; ++k){
        w = prod(1-v_list(arma::span(0,(C[class_index(k)]-2))))*v_list[(C[class_index(k)]-1)]/(1-v_list[j]) ;  
        temp = 1-u_list[class_index(k)]/w;
        if((v_upper >= temp) & (temp>=v_lower)){
          v_upper = temp;
        }
      }
      while(true){
        new_v = R::rbeta(1, alpha);
        if((new_v>=v_lower)&(new_v<=v_upper)){
          v_list[j] = new_v;
          if(j==0){
            w_list[j] = new_v;
          }else{
            w_list[j] = new_v*prod(1-v_list(arma::span(0,(j-1))));  
          }
          break;
        }
      }
    }
    
    // update U
    for(int j=0; j < n; ++j){
      u_list[j] = Rcpp::runif(1, 0, w_list[(C[j]-1)])[0]; // u_i
    }
    
    // update C
    cw = cumsum(w_list);
    thre = 1+max(-u_list);
    for(int j=0; j < n; ++j){
      is_zero=false;
      // setting n_i 
      uc = unique(C); ni =  hist(C, uc); // order by class
      ni[(C[j]-1)] -= 1;
      
      // if c_j is associated with no other obsevation then initialize
      if(ni[(C[j]-1)] == 0){
        is_zero=true;
        zero_class = C[j];
      }
      
      M = 1;
      for( int k=0; k < n; ++k ){
        if((k>uc.n_elem)|(thre < cw[k])){
          break;
        }else{
          M += 1;
        }
      }
      if( M > 5){
        cout << M << endl;
      }
      prob_class = arma::zeros(M);
      
      // calculate probability for existing class
      for( int k = 0; k < M; ++k ){
        var = 1/phis(k, 1);
        prob_class(k) = R::dnorm(ys[j], phis(k, 0), pow((var),0.5), false); //likelihood  
      }
      class_index = find(u_list(j) > w_list(arma::span(0,(M-1))));
      prob_class(class_index) = arma::zeros(class_index.n_elem);
      prob_class = prob_class/sum(prob_class); 
      
      // class sampling
      rnd = Rcpp::runif(1)[0];
      arma::vec cum_class = cumsum(prob_class); class_ind = 1;
      for( int k=0; k < M; ++k ){
        if( rnd < cum_class[k]){
          break;
        }else{
          class_ind += 1;
        }
      }
      
      // class update
      C[j] = class_ind;
      
      // standardize C
      if( is_zero & (class_ind != zero_class)){
        
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
    
    // save result
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











