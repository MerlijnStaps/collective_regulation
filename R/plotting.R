library(RColorBrewer)


# Functions

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

intensityColor <- function(values,min=0,max=1,resolution=1000,colorList=brewer.pal(5,"PuRd")){
  colorsToUse = colorRampPalette(colorList)(resolution)
  normalizedValues = ((values-min)/(max-min))
  bins = seq(0,1,length.out=resolution+1)
  roundedNormalizedValues = sapply(normalizedValues, function(u){sum(bins <= u)})
  return(colorsToUse[pmin(roundedNormalizedValues,resolution)])
}

# Colors

environment_colors = setNames(c("#C6E9AF","#447821"),c("B","A"))
environment_colors["good"] = environment_colors["A"]
environment_colors["bad"] = environment_colors["B"]

strategy_colors = setNames(NA, "other")
strategy_colors["always active"] = "#666666"
strategy_colors["always inactive"] = "#cccccc"
strategy_colors["naive flexible"] = "#aaeeff"
strategy_colors["always_active"] = "#666666"
strategy_colors["always_inactive"] = "#cccccc"
strategy_colors["naive_flexible"] = "#aaeeff"
strategy_colors["random_switching"] = "#B8CBFC"
strategy_colors["individual_regulation"] = "#aaeeff"
strategy_colors["individual_regulation_spontaneous"] = "#B8CBFC"
strategy_colors["none"] = "#eeeeee"
strategy_colors["collective_regulation_I"] = "#ffa500"
strategy_colors["collective_regulation_II"] = "coral"
strategy_colors["collective_regulation_IV"] = "cornflowerblue"
strategy_colors["collective_regulation_III"] = "purple"
strategy_colors["positive_feedback_active"] = "#ffa500"
strategy_colors["negative_feedback_active"] = "coral"
strategy_colors["positive_feedback_inactive"] = "cornflowerblue"
strategy_colors["negative_feedback_inactive"] = "purple"
strategy_colors["positive_feedback_both"] = "#91ee92"
strategy_colors["positive_negative_inactivity"] = "#eeee00"
strategy_colors["positive_negative_activity"] = "violet"
strategy_colors["positive_feedback_both_extra_activation"] = "#CAE83D"
strategy_colors["positive_feedback_both_extra_inactivation"] = "#A98EF0"
strategy_colors["positive_feedback_both_extra_inactivation_and_spontaneous_inactivation"] = "#B864FA"
strategy_colors["positive_feedback_both_spontaneous_activation"] = "#44ce82" 
strategy_colors["positive_feedback_both_extra_activation_and_spontaneous_activation"] = "#708F0D"
strategy_colors["positive_feedback_both_spontaneous_inactivation"] = "#90C7BF"
strategy_colors["activation"] = "#8c5f09"
strategy_colors["inactivation"] = "#551a8b"

fitness_increase_color <- function(x, max_value = 0.05, resolution = 20){
  if(x <= 0){return("white")}
  else{
    x = min(x, max_value)
    return(intensityColor(x,  min = 0, max = max_value, resolution = resolution))
    
  }
}

trait_color <- function(x){
  return(intensityColor(x, min=0, max = 1, resolution = 20, colorList= brewer.pal(9,"Purples")))
}

# Plots


lambda_bc_plot <- function(x, y, z, main = ""){
  plot(0,0,type="n",xlim=c(-2,2),ylim=c(-1,1),bty="l",xlab="b/c",ylab=expression("Environmental variability (" * lambda * ")"),
       main=main, xaxt="n", yaxt="n")
  axis(side=1, at=c(-2,-1,0,1,2), labels = c(0.25, 0.5, 1, 2, 4), tcl=-0.3, lwd =0, lwd.ticks=1)
  axis(side=2, c(-1,0,1), labels=c(0.1,1,10), tcl=-0.3, las = 2, lwd =0, lwd.ticks=1)
  dx = diff(sort(unique(x)))[1]
  dy = diff(sort(unique(y)))[1]
  rect(x - 0.5*dx, y - 0.5*dy, x + 0.5*dx, y + 0.5*dy, border=NA, col=z)
}

lambda_bc_plot_small <- function(x, y, z, main = "", bty="o"){
  plot(0,0,type="n",xlim=c(-2,2),ylim=c(-1,1),bty=bty,xlab="",ylab="",xaxt="n",yaxt="n",main=main)
  dx = diff(sort(unique(x)))[1]
  dy = diff(sort(unique(y)))[1]
  rect(x - 0.5*dx, y - 0.5*dy, x + 0.5*dx, y + 0.5*dy, border=NA, col=z)
}



interaction_rate_plot_alpha <- function(d, background_colors, analytical_predictions = T, main = ""){
  
  plot(0,0,type="n",xlim=c(-2,2),ylim=c(0,12),bty="n",xlab=expression("Observation rate (" * alpha * ")"),
       ylab=expression("Optimal interaction rate (" * beta * "*)"),bty="l",xaxt="n",yaxt="n",main=main)
  axis(side=1, at=seq(-2,2,1), labels=c(0.01,0.1,1,10,100), tck=-0.03, lwd=0, lwd.ticks = 1)
  axis(side=2, at=seq(0,12,3), tck=-0.03, lwd=0, lwd.ticks = 1, las = 2)
  
  d$x = log10(d$alpha)
  d$y = d$beta
  dx = diff(sort(unique(d$x)))[1]
  dy = diff(sort(unique(d$y)))[1]
  
  rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=background_colors)
  
  alphas = unique(d$alpha)
  optimal_betas = sapply(alphas, function(a){
    d2 = subset(d,d$alpha == a)
    d2 = d2[order(d2$winner_fitness,decreasing=T),]
    return(d2$beta[1])
  })
  
  if(analytical_predictions){ # add analytical prediction
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




interaction_rate_plot_lambda <- function(d, background_colors, analytical_predictions = T, main = ""){
  
  plot(0,0,type="n",xlim=c(-2,2),ylim=c(0,12),bty="n",xlab=expression("Rate of environmental change (" * lambda * ")"),
       ylab=expression("Optimal interaction rate (" * beta * "*)"),bty="l",xaxt="n",yaxt="n")
  axis(side=1, at=seq(-2,2,1), labels=c(0.01,0.1,1,10,100), tck=-0.03, lwd=0, lwd.ticks = 1)
  axis(side=2, at=seq(0,12,3), tck=-0.03, lwd=0, lwd.ticks = 1, las = 2)
  
  d$x = log10(d$lambda)
  d$y = d$beta
  dx = diff(sort(unique(d$x)))[1]
  dy = diff(sort(unique(d$y)))[1]
  
  rect(d$x - 0.5*dx, d$y - 0.5*dy, d$x + 0.5*dx, d$y + 0.5*dy, border=NA, col=background_colors)
  
  lambdas = sort(unique(d$lambda))
  optimal_betas = sapply(lambdas, function(l){
    d2 = subset(d,d$lambda == l)
    d2 = d2[order(d2$winner_fitness,decreasing=T),]
    return(d2$beta[1])
  })
  
  if(analytical_predictions){ # add analytical prediction
    # analytical prediction for low lambda
    predicted_slope = uniroot(function(z){b0*z/(z+1) - b0*log(z+1) + z/(1-z) - log(1/(1-z))}, c(1e-6,1-1e-6))$root
    predicted_values = rep(predicted_slope * a0, length(lambdas))
    w = which(predicted_values > 0 & predicted_values < 12)
    points(log10(lambdas)[w], predicted_values[w], type="l",lwd=3,lty=1,col="cornflowerblue")
    
    # analytical prediction for high lambda
    predicted_values = sapply(lambdas, function(lambda){lambda*(b0 - 1 + sqrt(b0**2 - 1)) - a0 })
    w = which(predicted_values > 0 & predicted_values < 12)
    points(log10(lambdas)[w],  predicted_values[w], type="l",lwd=3,lty=1,col="coral")
  }
  
  
  points(log10(lambdas)[optimal_betas < 12 & optimal_betas > 0.05],optimal_betas[optimal_betas < 12 & optimal_betas > 0.05],
         type="l",lwd=3)
}



