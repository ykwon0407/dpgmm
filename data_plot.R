#-------------------------------------------------
# Data visualization
#-------------------------------------------------
par(mfrow=c(1,2))
# Model 1
x=-5000:18000/1000
density_x= 0.2*dnorm(x, mean =0, sd=1) + 0.8*dnorm(x, mean =10, sd=2)
plot( x, density_x, lwd=0.3, main="Model 1", cex=0.3)

# Model 2
x=-2000:10000/200
proba_list= c(1,2,3,2,1)/9
sd_list = c(1,2,3,2,1)
mean_list = c(0,1,2,3,4)*10
density_x=rep(0, length(x))
for( i in 1:5){
  density_x = density_x + proba_list[i]*dnorm(x, mean_list[i], sd_list[i])  
}

plot( x, density_x, lwd=0.3, main="Model 2", cex=0.3)



