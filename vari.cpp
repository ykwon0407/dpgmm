// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>    
#include <vector>     

using namespace std;

// [[Rcpp::export]]
arma::mat vari_cpp(arma::vec ys, arma::vec C,
                    int max_iter = 11500, int thining = 100, int burn_in = 1500,
                    double alpha=1.0, double mu_0=2.0, double K=2.0, double v_0=0.2, double u_0=1.0, int T = 8){
  int count = 1; int n = ys.n_elem;
  int max_count = (max_iter-burn_in)/thining;
  double beta = (v_0*u_0)/2;
  arma::mat posterior = arma::zeros(max_count, 2*n); // result matrix
  arma::mat phis = arma::zeros(n, 2); // state matrix
  
  // initialize phi's from base measure. the number of elements in state cannot be larger than n.
  // initializer v and u
  double var=0;
  for( int i=0 ; i < n ; ++i ){ 
    beta = (v_0*u_0)/2;
    phis(i,1) = R::rgamma(v_0/2, 1/beta ); //tau_i
    var = 1/(K*phis(i,1));
    phis(i,0) = R::rnorm(mu_0, pow((var),0.5)); //mu_i 
  }
  arma::vec v_list = arma::zeros(T);
  arma::vec coef = arma::zeros(T);
  for( int i=0 ; i < T ; ++i ){ 
    v_list(i) = R::rbeta(1, alpha); // v_i
  }
  
  arma::uvec ni, class_index;  arma::vec uc, class_ys;
  arma::vec prob_class; arma::vec cw;
  cout<< "Iteration start!!" <<endl;
  
  double rnd;
  int class_ind, posterior_index, M;
  bool is_zero=false; int zero_class = 0; double a, b;
  // Iteration
  while(true){
    // iteration criteria
    if(max_iter < count){
      break;
    }
    
    uc = unique(C); ni = hist(C, unique(C)); M=uc.n_elem;
    coef = arma::zeros(T);
    // update V
    for(int j=0; j < T ; ++j){
      if(j < M){
        a = ni(j)+1; 
      }else{
        a = 1;
      }
      
      if( j == 0 ){
        b = n-a + alpha;  
      }else if( j < M){
        b = n-sum(ni(arma::span(0,j))) + alpha;
      }else{
        b = alpha;
      }
      
      v_list(j) = R::rbeta(a, b); // v_i
      if(j==0){
        coef(j) = 1;
      }else{
        coef(j) = coef((j-1))*b/(a+b);
      }
    }
    
    // update C
    for(int j=0; j < n; ++j){
      uc = unique(C); ni = hist(C, unique(C)); M=uc.n_elem;
      is_zero=false;
      // if c_j is associated with no other obsevation
      if(ni[(C[j]-1)] == 1){
        is_zero=true;
        zero_class = C[j];
      }
      
      // calculate probability for existing class
      prob_class = arma::zeros(T);
      for( int k = 0; k < T; ++k ){
        if(k < M){
          a = ni(k)+1; 
        }else{
          a = 1;
        }
        
        if( k == 0 ){
          b = n-a + alpha;  
        }else if(k < M){
          b = n-sum(ni(arma::span(0,k))) + alpha;
        }else{
          b = alpha;
        }
        
        var = 1/phis(k, 1);
        prob_class(k) = (a/(a+b))*coef[k]*R::dnorm(ys[j], phis(k, 0), pow((var),0.5), false); //likelihood  
      }
      prob_class = prob_class/sum(prob_class); 
      
      // class sampling
      rnd = Rcpp::runif(1)[0];
      arma::vec cum_class = cumsum(prob_class); class_ind = 1;
      for( int k=0; k < T; ++k ){
        if( rnd < cum_class[k]){
          break;
        }else{
          class_ind += 1;
        }
      }
      
      // class update
      C[j] = class_ind;
      
      if(class_ind>uc.n_elem){
        C[j] = uc.n_elem+1;
        phis((uc.n_elem),1)=phis((class_ind-1),1);
        phis((uc.n_elem),0)=phis((class_ind-1),0);
        // new phi's
        beta = (v_0*u_0)/2;
        phis((class_ind-1),1)=R::rgamma(v_0/2, 1/beta ); //tau_i
        var = 1/(K*phis((n-1),1));
        phis((class_ind-1),0)=R::rnorm(mu_0, pow((var),0.5)); //mu_i
      }
      
      // standardize C
      if(is_zero & (class_ind != zero_class)){
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
  
  cout << "Done!" << endl;
  return posterior;
}











