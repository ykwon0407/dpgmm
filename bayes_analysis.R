#--------------------------------------------------------------------
# Advanced Bayesian Midterm - Dirichelet Process Mixture Model
# 2013-22897 Kwon, Yongchan
# 
# Short descriptions:
# The following script compares the Dirichlet Process Mixture Model
# algorithms discussed in the Advanced Bayesian lectures. I followed 
# the papers Neal(2000) and lecture notes.
# 
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Default setting
#--------------------------------------------------------------------
set.seed(2357)
n = 200

#--------------------------------------------------------------------
# Source cpp files
#--------------------------------------------------------------------
setwd('~/Dropbox/Adv_Bayes_algorithm/algorithm/')
require(cluster)
require(data.table)
require(rbenchmark)
source('data_model.R')
Rcpp::sourceCpp('algo2.cpp')
Rcpp::sourceCpp('algo4.cpp')
Rcpp::sourceCpp('algo8.cpp')
Rcpp::sourceCpp('slice.cpp')
Rcpp::sourceCpp('vari.cpp')

#--------------------------------------------------------------------
# Generate data and exploratory data analysis for hyperparameter
#--------------------------------------------------------------------
data1 = data_model1(n) # Model 1
data2 = data_model2(n) # Model 2

par(mfrow=c(1,2))
hist(data1, prob=TRUE, col='green', breaks=40)
lines(density(data1), col='red', lwd=3)

hist(data2, prob=TRUE, col='green', ylim =c(0,0.05), breaks=20)
lines(density(data2, bw=1.5), col='red', lwd=3)

# Setting intial value for cluster 'C' using k-means algorithm.
C1 <- kmeans(data1, 2)$cluster #initial cluster
C2 <- kmeans(data2, 5)$cluster #initial cluster

plot(data1, col=C1) 
plot(data2, col=C2) 

#--------------------------------------------------------------------
# Algorithm 2
#--------------------------------------------------------------------

algo2 <- function(){
  mat_algo2_model1 <- algo2_cpp(data1, C1, mu_0= mean(data1))
  mat_algo2_model2 <- algo2_cpp(data2, C2, mu_0= mean(data2), K=1.0, v_0= 100.0)
  return(list(mat_algo2_model1, mat_algo2_model2))
}
#--------------------------------------------------------------------
# Algorithm 4
#--------------------------------------------------------------------
algo4 <- function(){
  mat_algo4_model1 <- algo4_cpp(data1, C1, mu_0= mean(data1))
  mat_algo4_model2 <- algo4_cpp(data2, C2, mu_0= mean(data2), K=0.15, v_0= 100.0)
  return(list(mat_algo4_model1, mat_algo4_model2))
}

#--------------------------------------------------------------------
# Algorithm 8
#--------------------------------------------------------------------
algo8 <- function(){
  mat_algo8_model1 <- algo8_cpp(data1, C1, mu_0= mean(data1))
  mat_algo8_model2 <- algo8_cpp(data2, C2, mu_0= mean(data2), K=0.15, v_0= 100.0)
  return(list(mat_algo8_model1, mat_algo8_model2))
}
#--------------------------------------------------------------------
# Sliced algorithm
#--------------------------------------------------------------------
slice <- function(){
  mat_slice_model1 <- slice_cpp(data1, C1, mu_0= mean(data1))
  mat_slice_model2 <- slice_cpp(data2, C2, mu_0= mean(data2), K=0.15, v_0= 100.0)
  return(list(mat_slice_model1, mat_slice_model2))
}

#--------------------------------------------------------------------
# Variational inference algorithm
#--------------------------------------------------------------------
vari <- function(){
  mat_vari_model1 <- vari_cpp(data1, C1, mu_0= mean(data1))
  mat_vari_model2 <- vari_cpp(data2, C2, mu_0= mean(data2), K=0.15, v_0= 100.0)
  return(list(mat_vari_model1, mat_vari_model2))
}

#--------------------------------------------------------------------
# Comparing by visualization
#--------------------------------------------------------------------
res_algo2 <- algo2();
res_algo4 <- algo4();
res_algo8 <- algo8();
res_slice <- slice();
res_vari <- vari();

par(mfrow=c(2,5))
plot(data1, col=as.factor(res_algo2[[1]][100,1:200]))
plot(data1, col=as.factor(res_algo4[[1]][100,1:200]))
plot(data1, col=as.factor(res_algo8[[1]][100,1:200]))
plot(data1, col=as.factor(res_slice[[1]][100,1:200]))
plot(data1, col=as.factor(res_vari[[1]][100,1:200]))

plot(data2, col=as.factor(res_algo2[[2]][100,1:200]))
plot(data2, col=as.factor(res_algo4[[2]][100,1:200]))
plot(data2, col=as.factor(res_algo8[[2]][100,1:200]))
plot(data2, col=as.factor(res_slice[[2]][100,1:200]))
plot(data2, col=as.factor(res_vari[[2]][100,1:200]))

#--------------------------------------------------------------------
# Comparing by speed
#--------------------------------------------------------------------
benchmark( algo2(), algo4(), algo8(), slice(), vari(), replications=1)
save.image("./bayes1.RData")



