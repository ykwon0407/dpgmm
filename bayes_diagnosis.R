setwd('~/Dropbox/algorithm/')
load("./bayes.RData")
ls()
install.packages("traceplot")

require(traceplot)
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
