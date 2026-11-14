source("R/plotting.R")
source("R/functions.R")

par(pty="m")
d1 = read.table("results/lambda_bc_log10.csv",header=T,stringsAsFactors=F,sep=",")
d1 = add_random_switching(d1, add_random_inactivation = F)


pdf("plots/figS2.pdf",width=4,height=4)

for(beta0 in c(0.4,1,4,10)){
  
  d = d1[d1$n == 0 & d1$beta == beta0,]
  
  d = calculate_winners(d, "positive_feedback_active")
  d$random_switching_fitness_increase[d$b <= d$c] = 0
  
  d$x = log2(d$b)
  d$y = log10(d$lambda)
  d = d[abs(d$y) <= 1,]
  d = d[abs(d$x) <= 2, ]
  
  add_flexibility_threshold <- function(){
    y = seq(-1,1,0.01)
    x = log2(1/10**y + 1)
    w = which(abs(x) <= 2)
    y = y[w]
    x = x[w]
    points(x,y,type="l",lwd=2,lty=2)
    points(-x,y,type="l",lwd=2,lty=2)
  }
  
  add_collective_regime_contour <- function(){
    d$z = d$winner == "positive_feedback_active"
    M = xtabs(z ~ x + y, data = d)
    l = contourLines(x = sort(unique(d$x)), y = sort(unique(d$y)), z = M, levels = 0.5)
    points(l[[1]]$x, l[[1]]$y, type = "l", lwd =2, lty= 2)
    points(l[[2]]$x, l[[2]]$y, type ="l", lwd = 2, lty = 2)
  }
  
  par(mgp=c(2,0.75,0))
  lambda_bc_plot(d$x, d$y, strategy_colors[d$random_switching_winner])
  add_collective_regime_contour()
}
dev.off()

