source("R/plotting.R")
source("R/functions.R")


add_flexibility_threshold <- function(){
  y = seq(-1,1,0.01)
  x = log2(1/10**y + 1)
  w = which(abs(x) <= 2)
  y = y[w]
  x = x[w]
  points(x,y,type="l",lwd=2,lty=2)
  points(-x,y,type="l",lwd=2,lty=2)
}

par(pty="m")
d1 = read.table("results/optimal_strategies.csv",header=T,stringsAsFactors=F,sep=",")

add_collective_regime_contour <- function(d, winner) {
  d$z <- d$winner == winner
  M <- xtabs(z ~ x + y, data = d)
  l <- contourLines(
    x = sort(unique(d$x)),
    y = sort(unique(d$y)),
    z = M,
    levels = 0.5
  )
  
  points(l[[1]]$x, l[[1]]$y, type = "l", lwd = 2, lty = 2)
  points(l[[2]]$x, l[[2]]$y, type = "l", lwd = 2, lty = 2)
}


for(motif in c("I", "II")) {
  strategy = paste("collective_regulation", motif, sep="_")
  
  if (motif == "I") {
    d1 <- add_random_switching(d1, add_random_inactivation = FALSE)
    d <- subset(d1, beta == 4 & n == 0)
    d <- calculate_winners(d, strategy)
    d$random_switching_fitness_increase[d$b <= d$c] <- 0
    
  } else if (motif == "II") {
    d1 <- add_random_switching(d1,add_random_inactivation = TRUE,add_random_activation = FALSE)
    d <- subset(d1, beta == 4 & n == 0)
    d <- calculate_winners(d, strategy)
    d$random_switching_fitness_increase[d$b >= d$c] <- 0
  }
  
  d$x <- log2(d$b)
  d$y <- log10(d$lambda)
  d <- d[abs(d$y) <= 1, ]
  d <- d[abs(d$x) <= 2, ]
  
  pdf(paste0("plots/optimal_strategies_motif_", motif, ".pdf"), width = 4, height = 4)
  
  par(mgp = c(2, 0.75, 0))
  
  lambda_bc_plot(
    d$x, d$y,
    strategy_colors[d$individual_winner],
    main = "Optimal strategy (no interactions)"
  )
  
  lambda_bc_plot(
    d$x, d$y,
    strategy_colors[d$winner],
    main = "Optimal strategy"
  )
  add_flexibility_threshold()
  
  lambda_bc_plot(
    d$x, d$y,
    sapply(d$fitness_increase, fitness_increase_color),
    main = "Fitness increase due to interactions"
  )
  add_flexibility_threshold()
  
  lambda_bc_plot(
    d$x, d$y,
    strategy_colors[d$random_switching_winner],
    main = "Optimal strategy\n(spontaneous switching included)"
  )
  add_collective_regime_contour(d, strategy)
  
  lambda_bc_plot(
    d$x, d$y,
    sapply(d$random_switching_fitness_increase, fitness_increase_color),
    main = "Fitness increase due to interactions\n(spontaneous switching included)"
  )
  dev.off()
}