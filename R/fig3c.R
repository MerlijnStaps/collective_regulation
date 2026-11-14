source("R/plotting.R")
source("R/functions.R")

a0 = 1

par(pty="m")
d1 = read.table("results/alpha_beta_tradeoff.csv",header=T,stringsAsFactors=F,sep=",")
d1 = add_random_switching(d1)
d1 = d1[d1$alpha0 == a0,]

strategies = c("naive_flexible","always_active","always_inactive","positive_feedback_active")
d1$winner = strategies[apply(d1[,paste(strategies,"fitness",sep="_")],1,which.max)]
d1$random_switching_winner = d1$winner
d1$random_switching_winner[d1$winner_fitness < d1$maximum_naive_fitness_analytical] = "random_switching"


d2 = read.table("results/lambda_beta_log10.csv",header=T,stringsAsFactors=F,sep=",")
d2 = add_random_switching(d2)
d2 = subset(d2,d2$b == 2 & d2$alpha == a0 & d2$n == 0)
strategies = c("naive_flexible","always_active","always_inactive","positive_feedback_active")
d2$winner = strategies[apply(d2[,paste(strategies,"fitness",sep="_")],1,which.max)]
d2$random_switching_winner = d2$winner
d2$random_switching_winner[d2$winner_fitness < d2$maximum_naive_fitness_analytical] = "random_switching"
d2$fitness_increase = apply(d2[,paste(strategies,"fitness",sep="_")],1,max) - apply(d2[,paste(strategies[1:3],"fitness",sep="_")],1,max)
d2$random_switching_fitness_increase = apply(d2[,paste(strategies,"fitness",sep="_")],1,max) - d2$maximum_naive_fitness_analytical


lambdas = sort(unique(d2$lambda))
optimal_betas = sapply(lambdas, function(l){
  d3 = subset(d2,d2$lambda == l)
  d3 = d3[order(d3$winner_fitness,decreasing=T),]
  return(d3$beta[1])
})



#plot(d$alpha, d$positive_feedback_active_fitness)

lambdas = sort(unique(d1$lambda))

pdf("plots/fig3c.pdf",width=4.5,height=4.5)
par(mgp=c(2,0.75,0))
for(i in 1:2){
  for(k0 in c(4)){
    plot(0,0,type="n",xlim=c(-2,2),ylim=c(0,12),bty="l",
         xlab=expression("Rate of environmental change (" * lambda * ")"),
         ylab=expression("Optimal interaction rate (" * beta * "*)"),
         xaxt="n",yaxt="n")
    axis(side=1, at=seq(-2,2,1), labels = c(0.01,0.1,1,10,100), tck=-0.03, lwd=0, lwd.ticks = 1)
    axis(side=2, at=seq(0,12,3), tck=-0.03, lwd=0, lwd.ticks = 1, las = 2)
    
    b0 = 2
    
    d = d1[d1$b == b0 & d1$k == k0,]
    
    
    d$x = log10(d$lambda)
    d$y = d$beta
    dx = diff(sort(unique(d$x)))[1]
    dy = diff(sort(unique(d$y)))[1]
    
    if(i==1){
      rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$random_switching_winner])
    }
    
    w = which(optimal_betas > 0.05 & optimal_betas < 12)
    points(log10(lambdas)[w], optimal_betas[w], col="gray", lwd=3, type="l")
    
    
    optimal_alphas = sapply(lambdas, function(lambda0){
      dd = d[d$lambda == lambda0,]
      j = which.max(dd$positive_feedback_active_fitness)
      #if(dd$random_switching_winner[j] == "random_switching" && dd$beta[j] > 0){
      #  return(a0)
      #}else{
        return(dd$alpha[which.max(dd$positive_feedback_active_fitness)])
      #}
    })
    x = log10(lambdas)
    y = k0*(a0-optimal_alphas)
    points(x,y,type="l",lwd=3)
    
  }
}
dev.off()

