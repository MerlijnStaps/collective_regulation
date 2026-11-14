source("R/plotting.R")
source("R/functions.R")

d1 = read.table("results/E1_processed.csv",header=T,stringsAsFactors=F,sep=",")
d1$strategy_pi = sapply(1:nrow(d1),function(i){
  call_strategy_pi(d1$piAA[i], d1$piAB[i], d1$piBA[i], d1$piBB[i])
})
d1$strategy_theta = sapply(1:nrow(d1),function(i){
  call_strategy_theta_rho(d1$thetaAA[i], d1$thetaAB[i], d1$thetaBA[i], d1$thetaBB[i], d1$rhoA[i], d1$rhoB[i])
})
d1$strategy_theta[d1$strategy_pi == "always active"] = "always active"
d1$strategy_theta[d1$strategy_pi == "always inactive"] = "always inactive"


d1$naive_fitness_analytical = d1$b*(d1$alpha+d1$lambda)/(2*(d1$alpha+2*d1$lambda)) - d1$c*(d1$lambda)/(2*(d1$alpha+2*d1$lambda))
d1$optimal_gamma_theoretical = d1$lambda*(pmax(d1$b/d1$c,d1$c/d1$b)-1+sqrt((pmax(d1$b/d1$c,d1$c/d1$b))^2-1))-d1$alpha
d1$maximum_gamma = 1e8
d1$optimal_gamma = pmax(pmin(d1$optimal_gamma_theoretical,d1$maximum_gamma),0)
d1$p_active_good = d1$optimal_gamma/(d1$alpha+d1$optimal_gamma) + d1$alpha/(d1$alpha+d1$optimal_gamma)*((d1$alpha+d1$optimal_gamma+d1$lambda)/(d1$alpha+d1$optimal_gamma+2*d1$lambda))
d1$p_active_bad = d1$optimal_gamma/(d1$alpha+d1$optimal_gamma) + d1$alpha/(d1$alpha+d1$optimal_gamma)*((d1$lambda)/(d1$alpha+d1$optimal_gamma+2*d1$lambda))
d1$maximum_naive_fitness_analytical[d1$b>d1$c] = 0.5*d1$p_active_good[d1$b>d1$c]*d1$b[d1$b>d1$c]  -0.5*d1$p_active_bad[d1$b>d1$c]*d1$c[d1$b>d1$c]
d1$maximum_naive_fitness_analytical[d1$b<=d1$c] = 0.5*(1-d1$p_active_bad[d1$b<=d1$c])*d1$b[d1$b<=d1$c]  -0.5*(1-d1$p_active_good[d1$b<=d1$c])*d1$c[d1$b<=d1$c]

d1$fitness_gain = d1$fitness - pmax(pmax(d1$naive_fitness_analytical,0),0.5*(d1$b-d1$c))
d1$fitness_gain_relative_to_random_switching = d1$fitness - pmax(pmax(pmax(d1$naive_fitness_analytical,0),0.5*(d1$b-d1$c)),d1$maximum_naive_fitness_analytical)

d1$unique_id <- match(do.call(paste, d1[c("alpha", "b")]), 
                      unique(do.call(paste, d1[c("alpha", "b")])))
d1 = d1[order(d1$unique_id, -d1$fitness),]
d1 = d1[!duplicated(d1$unique_id),]
d = d1[d1$b != 1,]
#d = d1

n0 = 5
beta0 = 4

examples = list(c(1,-0.6),c(1,-0.2),c(1,0.5))

pdf("plots/evolutionary_outcomes_n0_5_beta0_4.pdf",width=5.5,height=5)
par(mgp=c(2,0.75,0))
alpha_bc_plot(x = log2(d$b), y = log10(d$alpha), z = strategy_colors[d$strategy_pi], main = paste(c("n = ",n0,", beta = ",beta0),collapse=""))
alpha_bc_plot(x = log2(d$b), y = log10(d$alpha), z = strategy_colors[d$strategy_theta], main = paste(c("n = ",n0,", beta = ",beta0),collapse=""))
points(sapply(examples,function(x){x[[1]]}),sapply(examples,function(x){x[[2]]}),pch=21,lwd=2,bg="white")
gains = pmax(0, d$fitness_gain/(d$b+d$c))
alpha_bc_plot(x = log2(d$b), y = log10(d$alpha), z = intensityColor(gains,min=0,max=max(gains)), main = paste(c("n = ",n0,", beta = ",beta0),collapse=""))

dev.off()

pdf("plots/evolutionary_outcomes_n0_5_beta0_4_alleles.pdf",width=3,height=3)
par(mgp=c(2,0.75,0))
for(par in c("piAA","piAB","piBA","piBB","thetaAA","thetaAB","thetaBA","thetaBB","rhoA","rhoB","meanA","meanB")){
  alpha_bc_plot_small(x = log2(d$b), y = log10(d$alpha), z = intensityColor(d[,par],resolution=20), main = "")
}
dev.off()

# 
# isapprox <- function(x,y){abs(x-y) < 1e-6}
# # 
# pdf("plots/evolutionary_outcomes_examples.pdf",width=3,height=3)
# par(mgp=c(2,0.75,0))
# for(x in examples){
#   i = which(isapprox(log2(d$b), x[[1]]) & isapprox(log10(d$alpha),x[[2]]))
#   h = 0.4
#   plot(0,0,type="n",xlim=c(-0.5,5.5),ylim=c(0,1),bty="n",xlab="Active individuals",ylab="Frequency",yaxt="n")
#   #axis(side=1, at=0:5, tcl=-0.3, lwd =0, lwd.ticks=1)
#   axis(side=2, c(0,0.5,1), tcl=-0.3, las = 2, lwd =1, lwd.ticks=1)
#   rect(0:5, 0, 0:5+h, d[i,c("uA_0","uA_1","uA_2","uA_3","uA_4","uA_5")], border=NA, col=environment_colors["good"])
#   rect(0:5, 0, 0:5-h, d[i,c("uB_0","uB_1","uB_2","uB_3","uB_4","uB_5")], border=NA, col=environment_colors["bad"])
# }
# dev.off()
# 
pdf("plots/evolved_thetas.pdf",width=3.6,height=3.2)
par(mgp=c(2,0.75,0))

for(b0 in c(2)){
  d0 = d[d$b == b0,]
  d0 = d0[order(d0$alpha),]
  if(b0 == 2){
    alpha2 = d0[2,"alpha"]
    d0[2,] = d0[1,]
    d0[2,]$alpha = alpha2
  }
  for(p in c("thetaAA","thetaAB","thetaBA","thetaBB","rhoA","rhoB")){
    plot(0,0,type="n",xlim=c(-1,2),ylim=c(0,1),bty="l",xaxt="n",yaxt="n",xlab=expression("Observation rate (" * alpha * ")"),ylab="Evolved value",main=p)
    axis(side=1, at=c(-1,0,1,2), labels=c(0.1,1,10,100), tcl=-0.3, lwd =0, lwd.ticks=1)
    axis(side=2, at=c(0,.5,1), tcl=-0.3, las = 2, lwd =0, lwd.ticks=1)
    for(i in 1:nrow(d0)){
      rect(log10(d0$alpha)[i]-0.5*0.025, 0, log10(d0$alpha)[i] +0.5*0.025, 1, border=NA, col=add.alpha(strategy_colors[d0$strategy_theta[i]],1))
    }
    points(log10(d0$alpha), d0[,p], ylim=c(0,1),pch=19,lwd=3,cex=1,type="l")
    for(x in examples){
      i = which(isapprox(log2(d0$b), x[[1]]) & isapprox(log10(d0$alpha),x[[2]]))
      points(log10(d0$alpha)[i], d0[i,p], pch=21, cex=1.5, bg="white", lwd=2,xpd=T)
    }

  }
}

dev.off()
