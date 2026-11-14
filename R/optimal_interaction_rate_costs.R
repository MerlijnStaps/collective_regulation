source("R/plotting.R")
source("R/functions.R")

a0 = 1
b0 = 2
n0 = 0
k0 = 4
strategies = c("individual_regulation","always_active","always_inactive","collective_regulation_I")

ds = lapply(1:2, function(x){
  if(x==1){filename = "results/optimal_interaction_rates_costs.csv"}
  if(x==2){filename = "results/optimal_interaction_rates_lambda.csv"}
  d1 = read.table(filename,header=T,stringsAsFactors=F,sep=",")
  if(x==1){d = subset(d1,d1$b == b0 & d1$n == n0)}
  if(x==2){d = subset(d1,d1$b == b0 & d1$alpha == a0 & d1$n == n0)}
  d = add_random_switching(d, add_random_inactivation = FALSE)
  d = calculate_winners(d, "collective_regulation_I")
  return(d)
})
d_costly = ds[[1]]
d_costfree = ds[[2]]

lambdas = sort(unique(d_costfree$lambda))
optimal_betas_costfree = sapply(lambdas, function(l){
  d3 = subset(d_costfree,d_costfree$lambda == l)
  d3 = d3[order(d3$winner_fitness,decreasing=T),]
  return(d3$beta[1])
})
optimal_alphas_costly = sapply(lambdas, function(lambda0){
  dd = d_costly[d_costly$lambda == lambda0,]
  j = which.max(dd$collective_regulation_I_fitness)
  return(dd$alpha[which.max(dd$collective_regulation_I_fitness)])
})
optimal_betas_costly = k0*(a0-optimal_alphas_costly)

pdf("plots/optimal_interaction_rate_costs.pdf",width=4.5,height=4.5)
plot(0,0,type="n",xlim=c(-2,2),ylim=c(0,12),bty="l",
     xlab=expression("Rate of environmental change (" * lambda * ")"),
     ylab=expression("Optimal interaction rate (" * beta * "*)"),
     xaxt="n",yaxt="n")
axis(side=1, at=seq(-2,2,1), labels = c(0.01,0.1,1,10,100), tck=-0.03, lwd=0, lwd.ticks = 1)
axis(side=2, at=seq(0,12,3), tck=-0.03, lwd=0, lwd.ticks = 1, las = 2)
x = log10(lambdas)
w = which(optimal_betas_costfree > 0.05 & optimal_betas_costfree < 12)
points(x[w], optimal_betas_costfree[w], col="gray", lwd=3, type="l")
points(x, optimal_betas_costly, col="black", lwd=3, type="l")
dev.off()

