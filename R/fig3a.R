source("R/plotting.R")
source("R/functions.R")

par(pty="m")
d1 = read.table("results/alpha_beta_log10.csv",header=T,stringsAsFactors=F,sep=",")
d1 = add_random_switching(d1)

b0 = 2
n0 = 0
l0 = 1

d = subset(d1,d1$b == b0 & d1$alpha > 0 & d1$n == n0 & d1$lambda == l0)


strategies = c("naive_flexible","always_active","always_inactive","positive_feedback_active")
d$winner = strategies[apply(d[,paste(strategies,"fitness",sep="_")],1,which.max)]
d$random_switching_winner = d$winner
d$random_switching_winner[d$winner_fitness < d$maximum_naive_fitness_analytical] = "random_switching"
d$fitness_increase = apply(d[,paste(strategies,"fitness",sep="_")],1,max) - apply(d[,paste(strategies[1:3],"fitness",sep="_")],1,max)
d$random_switching_fitness_increase = apply(d[,paste(strategies,"fitness",sep="_")],1,max) - d$maximum_naive_fitness_analytical

pdf("plots/fig3a.pdf",width=4.5,height=4.5)
par(mgp=c(2,0.75,0))
for(i in 1:5){
  
  plot(0,0,type="n",xlim=c(-2,2),ylim=c(0,12),bty="n",xlab=expression("Observation rate (" * alpha * ")"),
       ylab=expression("Optimal interaction rate (" * beta * "*)"),bty="l",xaxt="n",yaxt="n")
  axis(side=1, at=seq(-2,2,1), labels=c(0.01,0.1,1,10,100), tck=-0.03, lwd=0, lwd.ticks = 1)
  axis(side=2, at=seq(0,12,3), tck=-0.03, lwd=0, lwd.ticks = 1, las = 2)
  
  d$x = log10(d$alpha)
  d$y = d$beta
  dx = diff(sort(unique(d$x)))[1]
  dy = diff(sort(unique(d$y)))[1]
  
  if(i==1){rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$winner])}
  if(i==2){rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$random_switching_winner])}
  if(i==3){rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=sapply(d$fitness_increase, function(x){fitness_increase_color(x, max_value = 0.08)}))}
  if(i==4){rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=sapply(d$random_switching_fitness_increase, function(x){fitness_increase_color(x, max_value = 0.08)}))}
  
  
  alphas = unique(d$alpha)
  optimal_betas = sapply(alphas, function(a){
    d2 = subset(d,d$alpha == a)
    d2 = d2[order(d2$winner_fitness,decreasing=T),]
    return(d2$beta[1])
  })
  
  if(i %in% c(1,2,5)){ # add analytical prediction
    # analytical prediction for low alpha
    predicted_values = l0*(b0 - 1 + sqrt(b0**2 - 1)) - alphas
    w = which(predicted_values > 0 & predicted_values < 12)
    points(log10(alphas)[w], predicted_values[w] , type="l",lwd=3,lty=1,col="coral")
    
    # analytical prediction for high alpha
    predicted_slope = uniroot(function(z){b0*z/(z+1) - b0*log(z+1) + z/(1-z) - log(1/(1-z))}, c(1e-6,1-1e-6))$root
    predicted_values = sapply(alphas, function(alpha){alpha * predicted_slope})
    w = which(predicted_values > 0 & predicted_values < 12)
    points(log10(alphas)[w], predicted_values[w], type="l",lwd=3,lty=1,col="cornflowerblue")
  }
  
  points(log10(alphas)[optimal_betas < 12],optimal_betas[optimal_betas < 12],type="l",lwd=3)
}
dev.off()

