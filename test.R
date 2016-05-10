#--------------------------------------------------------------------
# Source cpp files
#--------------------------------------------------------------------
set.seed(100)
Rcpp::sourceCpp('~/project/bayesian/algorithm/generator.cpp')
Rcpp::sourceCpp('~/project/bayesian/algorithm/algo2.cpp')

#--------------------------------------------------------------------
# Generate data
#--------------------------------------------------------------------
list_data = generator_cpp(5, alpha = 0.6, tol=0.001)
theta=list_data[[1]]
data=list_data[[2]]

# Setting intial value for 'c' using k-means algorithm.
# Let the number of cluster be logarithm of the number of data.
require(cluster)
C <- kmeans(ys, round(log(length(ys)))+1)$cluster #initial cluster

#--------------------------------------------------------------------
# Algorithm 2
#--------------------------------------------------------------------
algorithm2 <- function(ys){
  print('Start: algorithm 2')
  post_algo2 <- algo2_cpp(ys, C)
  print('Done: algorithm 2')
}






# posterior 계산해보고
# algorithm별로 차이점 생각한뒤
# code 짜기
# posteriror distribution 구해서 histogram 그리기