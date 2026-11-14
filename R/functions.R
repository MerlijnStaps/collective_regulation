call_strategy_pi <- function(pAA, pAB, pBA, pBB){
  if(pAA == 0 & pBA == 0){return("always inactive")}
  if(pBB == 0 & pAB == 0){return("always active")}
  if(pAA == 1 & pAB == 0 & pBA == 0 & pBB == 1){return("naive flexible")}
  if(pAA == 1 & pAB > 0 & pBA == 0 & pBB == 1){return("anticipate inactivation")}
  if(pAA == 1 & pAB == 0 & pBA > 0 & pBB == 1){return("anticipate activation")}
  if(pAA == 1 & pAB == 0 & pBA == 0 & pBB < 1){return("delay inactivation")}
  if(pAA < 1 & pAB == 0 & pBA == 0 & pBB == 1){return("delay activation")}
  return("other")
}

call_strategy_theta <- function(tAA, tAB, tBA, tBB){
  if(tAA == 0 & tAB == 0 & tBA == 0 & tBB == 0){return("none")}
  if(tAA > 0 & tAB == 0 & tBA == 0 & tBB == 0){return("positive_feedback_active")}
  if(tAA == 0 & tAB > 0 & tBA == 0 & tBB == 0){return("negative_feedback_active")}
  if(tAA == 0 & tAB == 0 & tBA > 0 & tBB == 0){return("negative_feedback_inactive")}
  if(tAA == 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_inactive")}
  if(tAA > 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_both")}
  if(tAA > 0 & tAB == 0 & tBA > 0 & tBB == 0){return("activation")}
  if(tAA == 0 & tAB > 0 & tBA == 0 & tBB > 0){return("inactivation")}
  if(tAA > 0 & tAB > 0 & tBA == 0 & tBB == 0){return("positive_negative_activity")}
  if(tAA == 0 & tAB == 0 & tBA > 0 & tBB > 0){return("positive_negative_inactivity")}
  if(tAA > 0 & tAB == 0 & tBA > 0 & tBB > 0){return("positive_feedback_both_extra_activation")}
  if(tAA > 0 & tAB > 0 & tBA == 0 & tBB > 0){return("positive_feedback_both_extra_inactivation")}
  return("other")
}


call_strategy_theta_rho <- function(tAA, tAB, tBA, tBB,rA,rB){
  if(rA == 0 & rB == 0){
    if(tAA == 0 & tAB == 0 & tBA == 0 & tBB == 0){return("none")}
    if(tAA > 0 & tAB == 0 & tBA == 0 & tBB == 0){return("positive_feedback_active")}
    if(tAA == 0 & tAB > 0 & tBA == 0 & tBB == 0){return("negative_feedback_active")}
    if(tAA == 0 & tAB == 0 & tBA > 0 & tBB == 0){return("negative_feedback_inactive")}
    if(tAA == 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_inactive")}
    if(tAA > 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_both")}
    if(tAA > 0 & tAB == 0 & tBA > 0 & tBB > 0){return("positive_feedback_both_extra_activation")}
    if(tAA > 0 & tAB > 0 & tBA == 0 & tBB > 0){return("positive_feedback_both_extra_inactivation")}
  }
  if(rA > 0){
    if(tAA > 0 & tAB == 0 & tBA == 0 & tBB == 0){return("activation")}
    if(tAA == 0 & tAB == 0 & tBA > 0 & tBB > 0){return("positive_negative_inactivity")}
    if(tAA > 0 & tAB == 0 & tBA > 0 & tBB > 0){return("positive_feedback_both_extra_activation_and_spontaneous_activation")}
    if(tAA > 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_both_spontaneous_activation")}
  }
  if(rB > 0){
    if(tAA == 0 & tAB == 0 & tBA == 0 & tBB > 0){return("inactivation")}
    if(tAA > 0 & tAB > 0 & tBA == 0 & tBB == 0){return("positive_negative_activity")}
    if(tAA > 0 & tAB > 0 & tBA == 0 & tBB > 0){return("positive_feedback_both_extra_inactivation_and_spontaneous_inactivation")}
    if(tAA > 0 & tAB == 0 & tBA == 0 & tBB > 0){return("positive_feedback_both_spontaneous_inactivation")}
  }
  
  return("other")
}

clamp <- function(x,a,b){
  if(x < a){return(a)}
  if(x > b){return(b)}
  return(x)
}

add_random_switching <- function(d1, add_random_inactivation = T, add_random_activation = T){
  d1$naive_fitness_analytical = d1$b*(d1$alpha+d1$lambda)/(d1$alpha+2*d1$lambda) - d1$c*d1$lambda/(d1$alpha+2*d1$lambda)
  d1$optimal_gamma_theoretical = d1$lambda*(pmax(d1$b/d1$c,d1$c/d1$b)-1+sqrt((pmax(d1$b/d1$c,d1$c/d1$b))^2-1))-d1$alpha
  d1$optimal_gamma = sapply(d1$optimal_gamma_theoretical, function(x){clamp(x,0,1e8)})
  d1$p_active_good = d1$optimal_gamma/(d1$alpha+d1$optimal_gamma) + d1$alpha/(d1$alpha+d1$optimal_gamma)*((d1$alpha+d1$optimal_gamma+d1$lambda)/(d1$alpha+d1$optimal_gamma+2*d1$lambda))
  d1$p_active_bad = d1$optimal_gamma/(d1$alpha+d1$optimal_gamma) + d1$alpha/(d1$alpha+d1$optimal_gamma)*((d1$lambda)/(d1$alpha+d1$optimal_gamma+2*d1$lambda))
  d1$maximum_naive_fitness_analytical = d1$p_active_good*d1$b  -d1$p_active_bad*d1$c
  if(add_random_inactivation){
    d1$maximum_naive_fitness_analytical[d1$b<=d1$c] = (1-d1$p_active_bad[d1$b<=d1$c])*d1$b[d1$b<=d1$c]  -(1-d1$p_active_good[d1$b<=d1$c])*d1$c[d1$b<=d1$c]
  }else{
    d1$maximum_naive_fitness_analytical[d1$b <= d1$c] = d1$naive_fitness_analytical[d1$b <= d1$c]
  }
  if(!add_random_activation){
    d1$maximum_naive_fitness_analytical[d1$b>=d1$c] = d1$naive_fitness_analytical[d1$b >=d1$c]
  }
  return(d1)
}

calculate_winners <- function(d, strategy){
  strategies = c("naive_flexible","always_active","always_inactive",strategy) 
  d$winner = strategies[apply(d[,paste(strategies,"fitness",sep="_")],1,which.max)]
  d$winner_fitness = apply(d[,paste(strategies,"fitness",sep="_")],1,max)
  d$individual_winner = strategies[apply(d[,paste(strategies[1:3],"fitness",sep="_")],1,which.max)]
  d$individual_fitness = apply(d[,paste(strategies[1:3],"fitness",sep="_")],1,max)
  d$fitness_increase = (d$winner_fitness - d$individual_fitness)/(d$b+d$c)
  d$random_switching_winner = d$winner
  d$random_switching_winner[d$winner_fitness < d$maximum_naive_fitness_analytical] = "random_switching"
  d$random_switching_winner[d$random_switching_winner == "naive_flexible"] = "random_switching"
  d$random_switching_fitness_increase = (d$winner_fitness - d$maximum_naive_fitness_analytical)/(d$b+d$c)
  return(d)
}
