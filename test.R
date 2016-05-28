#--------------------------------------------------------------------
# Advanced Bayesian Midterm - Dirichelet Process Mixture Model
# 2013-22897 Kwon, Yongchan
# 
# Short descriptions:
# The following script compares the Dirichlet Process Mixture Model
# algorithms discussed in the Advanced Bayesian lectures. I followed 
# the papers Neal(2000), , and Biel et al.(2006).
# 
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Default setting
#--------------------------------------------------------------------
n = 200
set.seed(2357)

#--------------------------------------------------------------------
# Source cpp files
#--------------------------------------------------------------------
setwd('~/Dropbox/algorithm/')
require(cluster)
require(data.table)
source('data_model.R')
Rcpp::sourceCpp('~/Dropbox/algorithm/algo2.cpp')
Rcpp::sourceCpp('~/Dropbox/algorithm/algo4.cpp')
Rcpp::sourceCpp('~/Dropbox/algorithm/algo8.cpp')

#--------------------------------------------------------------------
# Generate data and exploratory data analysis for hyperparameter
#--------------------------------------------------------------------
data1 = data_model1(n) # Model 1
data2 = data_model2(n) # Model 2

plot(data1)
plot(data2)

hist(data1, prob=TRUE, col='green')
lines(density(data1), col='red', lwd=3)

hist(data2, prob=TRUE, col='green')
lines(density(data2), col='red', lwd=3)

# Setting intial value for 'c' using k-means algorithm.
# Let the number of cluster be logarithm of the number of data.
C1 <- kmeans(data1, 2)$cluster #initial cluster
C2 <- kmeans(data2, 3)$cluster #initial cluster

plot(data1, col=C1) 
plot(data2, col=C2) 

#--------------------------------------------------------------------
# Algorithm 2
#--------------------------------------------------------------------
mat_algo2_model1 <- algo2_cpp(data1, C1, alpha=0.1, K=1.0, mu_0= mean(data1))
mat_algo2_model2 <- algo2_cpp(data2, C2, alpha=0.1, K=3.0, mu_0= mean(data2), u_0 = 0.7)

plot(data1, col=as.factor(mat_algo2_model1[1,1:200]))
plot(data2, col=as.factor(mat_algo2_model2[1,1:200]))
as.factor(mat_algo2_model2[1,1:200])

#--------------------------------------------------------------------
# Algorithm 4
#--------------------------------------------------------------------
mat_algo4_model1 <- algo4_cpp(data1, C1, alpha=0.1, K=1.0, mu_0= mean(data1))
mat_algo4_model2 <- algo4_cpp(data2, C2, alpha=0.1, K=1.5, mu_0= mean(data2), u_0 = 1.2)

plot(data1, col=as.factor(mat_algo4_model1[1,1:200]))
plot(data2, col=as.factor(mat_algo4_model2[1,1:200]))

#--------------------------------------------------------------------
# Algorithm 8
#--------------------------------------------------------------------
mat_algo8_model1 <- algo8_cpp(data1, C1, alpha=0.1, K=1.0, mu_0= mean(data1))
mat_algo8_model2 <- algo8_cpp(data2, C2, alpha=0.1, K=1.5, mu_0= mean(data2), u_0 = 1.2)

plot(data1, col=as.factor(mat_algo8_model1[1,1:200]))
plot(data2, col=as.factor(mat_algo8_model2[1,1:200]))

#--------------------------------------------------------------------
# Sliced algorithm
#--------------------------------------------------------------------
mat_slice_model1 <- slice_cpp(data1, C1, alpha=0.1, K=1.0, mu_0= mean(data1))
mat_slice_model2 <- slice_cpp(data2, C2, alpha=0.1, K=1.5, mu_0= mean(data2), u_0 = 1.2)

plot(data1, col=as.factor(mat_slice_model1[1,1:200]))
plot(data2, col=as.factor(mat_slice_model2[1,1:200]))

#--------------------------------------------------------------------
# Variational algorithm method
#--------------------------------------------------------------------
mat_vari_model1 <- algo_var_cpp(data1, C1, alpha=0.1, K=1.0, mu_0= mean(data1))
mat_vari_model2 <- algo_var_cpp(data2, C2, alpha=0.1, K=1.5, mu_0= mean(data2), u_0 = 1.2)

plot(data1, col=as.factor(mat_vari_model1[1,1:200]))
plot(data2, col=as.factor(mat_vari_model2[1,1:200]))












