source("R/plotting.R")
source("R/functions.R")

par(pty="m")
d1 = read.table("results/alpha_bc_log10.csv",header=T,stringsAsFactors=F,sep=",")
colnames(d1) = sapply(colnames(d1),function(u){gsub("\\.","_",u)})

d1 = add_random_switching(d1, add_random_inactivation = F)

beta0 = 4
d = subset(d1,d1$beta == beta0 & d1$n == 0)

strategies = c("naive_flexible","always_active","always_inactive","positive_feedback_active") 

d$winner = strategies[apply(d[,paste(strategies,"fitness",sep="_")],1,which.max)]
d$individual_winner = strategies[apply(d[,paste(strategies[1:3],"fitness",sep="_")],1,which.max)]
d$winner_fitness = apply(d[,paste(strategies,"fitness",sep="_")],1,max)
d$fitness_increase = (apply(d[,paste(strategies,"fitness",sep="_")],1,max) - apply(d[,paste(strategies[1:3],"fitness",sep="_")],1,max))/(d$b+d$c)
d$random_switching_winner = d$winner
d$random_switching_winner[d$winner_fitness < d$maximum_naive_fitness_analytical] = "random_switching"
d$random_switching_fitness_increase = 2*(apply(d[,paste(strategies,"fitness",sep="_")],1,max) - d$maximum_naive_fitness_analytical)/(d$b+d$c)

d$x = log2(d$b)
d$y = log10(d$alpha)

pdf("plots/optimal_strategies_log10_beta4_n0.pdf",width=4,height=4)
par(mgp=c(2,0.75,0))
for(i in 0:4){
  
  plot(0,0,type="n",xlim=c(-4,4),ylim=c(-1,2),bty="l",xlab="b/c",ylab=expression("Observation rate (" * alpha * ")"),main="",xaxt="n",yaxt="n")
  axis(side=1, at=c(-4,-2,0,2,4), labels=c(0.0625, 0.25, 1, 4, 16), tcl=-0.375, lwd=0, lwd.ticks = 1)
  axis(side=2, at=c(-1,0,1,2), labels=c(0.1,1,10,100), tcl=-0.375, las = 2, lwd =0, lwd.ticks=1)
  dx = diff(sort(unique(d$x)))[1]
  dy = diff(sort(unique(d$y)))[1]
  if(i==0){
    rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$individual_winner])
  }
  if(i==1){ # optimal strategies
    rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$winner])
  }
  if(i==2){ # fitness increase 
    rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=sapply(d$fitness_increase, fitness_increase_color))
  }
  if(i==3){ # fitness increase 
    rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=sapply(d$random_switching_fitness_increase, fitness_increase_color))
  }
  if(i==4){ # optimal strategies (with random switching)
    rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=strategy_colors[d$random_switching_winner])
  }
  if(i %in% 1:2){ # plot flexibility threshold
    y = seq(-1,2,0.01)
    x = log2(10** y + 1)
    y = y[abs(x) <= 4]
    x = x[abs(x) <= 4]
    points(x,y,type="l",lwd=2,lty=2)
    points(-x,y,type="l",lwd=2,lty=2)}
  if(i==4){ # plot contour of collective strategy
    for(j in 4:5){
      d$z = d$winner == strategies[j]
      M = xtabs(z ~ x + y, data = d)
      l = contourLines(x = sort(unique(d$x)), y = sort(unique(d$y)), z = M, levels = 0.5)
      points(l[[1]]$x, l[[1]]$y, type = "l", lwd =2, lty= 2)
      points(l[[2]]$x, l[[2]]$y, type ="l", lwd = 2, lty = 2)
    }
  }
}
dev.off()

