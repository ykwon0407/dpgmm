data_model1 <- function(n=5){
  rnd = runif(n) 
  data = (rnd <= 0.2)*rnorm(n, mean =0, sd=1) + (rnd > 0.2)*rnorm(n, mean =10, sd=2)
  return(data)
}

data_model2 <- function(n=5){
  proba_list= c(1,2,3,2,1)/9
  sd_list = c(1,2,3,2,1)
  mean_list = c(0,1,2,3,4)*10
  rnd = rmultinom(n, size=1, prob=proba_list)
  
  norm_mat = matrix(0, nr=5, nc=n)
  for( i in 1:5){
    norm_mat[i,] = rnorm(n, mean_list[i], sd_list[i])
  }
  data = colSums(rnd*norm_mat)
  return(data)
}

