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
strategy_colors["none"] = "#eeeeee"
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

alpha_bc_plot <- function(x, y, z, main = ""){
  plot(0,0,type="n",xlim=c(-3,3),ylim=c(-1,2),bty="l",xlab="b/c",ylab=expression("Observation rate (" * alpha * ")"),
       main=main, xaxt="n", yaxt="n")
  axis(side=1, at=c(-3,-2,-1,0,1,2,3), labels = c(0.125, 0.25, 0.5, 1, 2, 4, 8), tcl=-0.3, lwd =0, lwd.ticks=1)
  axis(side=2, c(-1,0,1,2), labels=c(0.1,1,10,100), tcl=-0.3, las = 2, lwd =0, lwd.ticks=1)
  dx = diff(sort(unique(x)))[1]
  dy = diff(sort(unique(y)))[1]
  rect(x - 0.5*dx, y - 0.5*dy, x + 0.5*dx, y + 0.5*dy, border=NA, col=z)
}


alpha_bc_plot_small <- function(x, y, z, main = "", bty = "o"){
  plot(0,0,type="n",xlim=c(-3,3),ylim=c(-1,2),xlab="b/c",ylab=expression("Observation rate (" * alpha * ")"),
       main=main, xaxt="n", yaxt="n",bty=bty)
  #axis(side=1, at=c(-3,-2,-1,0,1,2,3), labels = c(0.125, 0.25, 0.5, 1, 2, 4, 8), tcl=-0.3, lwd =0, lwd.ticks=1)
  #axis(side=2, c(-1,0,1,2), labels=c(0.1,1,10,100), tcl=-0.3, las = 2, lwd =0, lwd.ticks=1)
  dx = diff(sort(unique(x)))[1]
  dy = diff(sort(unique(y)))[1]
  rect(x - 0.5*dx, y - 0.5*dy, x + 0.5*dx, y + 0.5*dy, border=NA, col=z)
}


alpha_bc_plot_old <- function(x, y, z, main = ""){
  plot(0,0,type="n",xlim=c(-3,3),ylim=c(0,6),bty="l",xlab="b/c",ylab=expression("Observation rate (" * alpha * ")"),
       main=main, xaxt="n", yaxt="n")
  axis(side=1, at=c(-3,-2,-1,0,1,2,3), labels = c(0.125, 0.25, 0.5, 1, 2, 4, 8), tcl=-0.3, lwd =0, lwd.ticks=1)
  axis(side=2, c(0,2,4,6), tcl=-0.3, las = 2, lwd =0, lwd.ticks=1)
  dx = diff(sort(unique(x)))[1]
  dy = diff(sort(unique(y)))[1]
  rect(x - 0.5*dx, y - 0.5*dy, x + 0.5*dx, y + 0.5*dy, border=NA, col=z)
}

