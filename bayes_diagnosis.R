setwd('~/Dropbox/Adv_Bayes_algorithm/algorithm/')
load("./bayes1.RData")
ls()
require(data.table)

#--------------------------------------
# Trace plot
#--------------------------------------
klist = c(30, 60, 90, 120, 150)
par(mfrow=c(1,5))

for( k in klist){
  plot(res_algo2[[1]][,k], ylab="", xlab="")
  plot(res_algo4[[1]][,k], ylab="", xlab="")
  plot(res_algo8[[1]][,k], ylab="", xlab="")
  plot(res_slice[[1]][,k], ylab="", xlab="")
  plot(res_vari[[1]][,k], ylab="", xlab="")
}

par(mfrow=c(1,5))
for( k in klist){
  plot(res_algo2[[2]][,k], ylab="", xlab="")
  plot(res_algo4[[2]][,k], ylab="", xlab="")
  plot(res_algo8[[2]][,k], ylab="", xlab="")
  plot(res_slice[[2]][,k], ylab="", xlab="")
  plot(res_vari[[2]][,k], ylab="", xlab="")
}

#--------------------------------------
# Posterior plot
#--------------------------------------
res_list <- list(res_algo2, res_algo4, res_algo8, res_slice, res_vari)

# Model 1
par(mfrow=c(1,1))
x=-5000:18000/1000
density_x = rep(0, length(x))
density_x = 0.2*dnorm(x, mean =0, sd=1) + 0.8*dnorm(x, mean =10, sd=2)
plot( x, density_x, lwd=4.0, col=6, main="Model 1", cex=0.3, ylim=c(0,0.2))

count = 1
for( res in res_list){
  raw <- res[[1]][100,] 
  DT <- data.table(mu=raw[1:200], pre=raw[201:400], class= as.factor(raw[1:200]))
  uni_DT <- unique(DT)
  pred_p = pred_mu = pred_sd = rep(0, nrow(uni_DT))
  for( i in 1:nrow(uni_DT)){
    pred_p[i] <- mean(DT$class == (uni_DT$class[i]))
    pred_mu[i] <- uni_DT$mu[i]
    pred_sd[i] <- (1/uni_DT$pre[i])^(1/2)
  }
  density_x = rep(0, length(x))
  for( i in 1:nrow(uni_DT)){
    density_x = density_x + pred_p[i]*dnorm(x, mean = uni_DT$mu[i], sd = pred_sd[i]) 
  }
  points( x, density_x, lwd=0.3, cex=0.3, col=count)
  count = count +1;
}
legend("topleft", col=1:5, lty=1, lwd=2, legend=c("algo2","no gaps","algo8","slice","vari"))



# Model 2
par(mfrow=c(1,1))
x=-2000:10000/200
proba_list= c(1,2,3,2,1)/9
sd_list = c(1,2,3,2,1)
mean_list = c(0,1,2,3,4)*10
density_x=rep(0, length(x))

for( i in 1:5){
  density_x = density_x + proba_list[i]*dnorm(x, mean_list[i], sd_list[i])  
}
plot( x, density_x, col=6, lwd=3.0, main="Model 2", cex=0.3, ylim=c(0, 0.12))

count = 1
for( res in res_list){
  raw <- res[[2]][100,] 
  DT <- data.table(mu=raw[1:200], pre=raw[201:400], class= as.factor(raw[1:200]))
  uni_DT <- unique(DT)
  pred_p = pred_mu = pred_sd = rep(0, nrow(uni_DT))
  for( i in 1:nrow(uni_DT)){
    pred_p[i] <- mean(DT$class == (uni_DT$class[i]))
    pred_mu[i] <- uni_DT$mu[i]
    pred_sd[i] <- (1/uni_DT$pre[i])^(1/2)
  }
  density_x = rep(0, length(x))
  for( i in 1:nrow(uni_DT)){
    density_x = density_x + pred_p[i]*dnorm(x, mean = uni_DT$mu[i], sd = pred_sd[i]) 
  }
  points( x, density_x, lwd=0.3, cex=0.3, col=count)
  count = count +1;
}
legend("topleft", col=1:5, lty=1, lwd=2, legend=c("algo2","no gaps","algo8","slice","vari"))

plot(density(data2, bw=1.0))
