source("R/plotting.R")
source("R/functions.R")

for(beta0 in c(4,40)){
  d1 = read.table("results/L3_processed.csv",header=T,stringsAsFactors=F,sep=",")
  d1 = d1[d1$beta == beta0,]
  
  d1$strategy_pi = sapply(1:nrow(d1),function(i){
    call_strategy_pi(d1$piAA[i], d1$piAB[i], d1$piBA[i], d1$piBB[i])
  })
  d1$strategy_theta = sapply(1:nrow(d1),function(i){
    call_strategy_theta_rho(d1$thetaAA[i], d1$thetaAB[i], d1$thetaBA[i], d1$thetaBB[i], d1$rhoA[i], d1$rhoB[i])
  })
  d1$strategy_theta[d1$strategy_pi == "always active"] = "always active"
  d1$strategy_theta[d1$strategy_pi == "always inactive"] = "always inactive"
  
  d1 = add_random_switching(d1)
  d1$individual_fitness = pmax(pmax(d1$naive_fitness_analytical,0),d1$b-d1$c)
  
  d1$fitness_gain = (d1$fitness - d1$individual_fitness)/(d1$b+d1$c)
  d1$fitness_gain_relative_to_random_switching = (d1$fitness - pmax(d1$individual_fitness,d1$maximum_naive_fitness_analytical))/(d1$b+d1$c)
  
  d1$unique_id <- match(do.call(paste, d1[c("lambda", "b")]), 
                        unique(do.call(paste, d1[c("lambda", "b")])))
  d1 = d1[order(d1$unique_id, -d1$fitness),]
  d1 = d1[!duplicated(d1$unique_id),]
  d = d1
  d$x = log2(d$b)
  d$y = log10(d$lambda)
  d = d[abs(d$x) <= 2,]
  d = d[d$x != 0,]
  
  pdf(paste("plots/fig4_",beta0,".pdf",sep=""),width=4,height=4)
  par(mgp=c(2,0.75,0))
  lambda_bc_plot(x = d$x, y = d$y, z = strategy_colors[d$strategy_pi])
  lambda_bc_plot(x = d$x, y = d$y, z = strategy_colors[d$strategy_theta])
  lambda_bc_plot(x = d$x, y = d$y, z = intensityColor(d$individual_fitness, max=2))
  lambda_bc_plot(x = d$x, y = d$y, z = sapply(d$fitness_gain, fitness_increase_color))
  dev.off()
  
  pdf(paste("plots/figS4_",beta0,"_alleles.pdf",sep=""),width=3,height=3)
  par(mgp=c(2,0.75,0))
  for(par in c("piAA","piAB","piBA","piBB","thetaAA","thetaAB","thetaBA","thetaBB","rhoA","rhoB","meanA","meanB")){
    lambda_bc_plot_small(x = d$x, y = d$y, z = sapply(d[,par],trait_color), bty="n", main=par)
  }
  plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  for(i in 0:10){
    rect(i/11,0,(i+1)/11,1,border=NA,col=trait_color(i/10))
  }
  rect(0,0,1,1)
  dev.off()
}


d$theta_max = apply(d[,c("thetaAA","thetaAB","thetaBA","thetaBB")],1,max)
lambda_bc_plot(x = d$x, y = d$y, z = sapply(d$theta_max, trait_color))
