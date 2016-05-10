#--------------------------------------------------------------------
# Advanced Bayesian - Dirichelet Process Mixture Model
# 2013-22897 Kwon, Yongchan
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Default setting
#--------------------------------------------------------------------
n = 5
alpha = 0.6
tol = 0.001
set.seed(100)

#--------------------------------------------------------------------
# Source cpp files
#--------------------------------------------------------------------
Rcpp::sourceCpp('~/project/bayesian/algorithm/generator.cpp')
Rcpp::sourceCpp('~/project/bayesian/algorithm/algo2.cpp')

#--------------------------------------------------------------------
# Generate data
#--------------------------------------------------------------------
list_data = generator_cpp(n, alpha = alpha, tol=tol)
theta=list_data[[1]]
data=as.vector(list_data[[2]])

# Setting intial value for 'c' using k-means algorithm.
# Let the number of cluster be logarithm of the number of data.
require(cluster)
C <- kmeans(data, round(log(length(data)))+1)$cluster #initial cluster

#--------------------------------------------------------------------
# Algorithm 2
#--------------------------------------------------------------------
algorithm2 <- function(data){
  print('Start: algorithm 2')
  post_algo2 <- algo2_cpp(data, C, skip=1000)
  print('Done: algorithm 2')
  return(post_algo2)
}


res_algo2 <- algorithm2(data)


# posterior 계산해보고
# algorithm별로 차이점 생각한뒤
# code 짜기
# posteriror distribution 구해서 histogram 그리기