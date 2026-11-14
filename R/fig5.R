source("R/plotting.R")
source("R/functions.R")

d1 = read.table(gzfile("results/M2_Strategy3_results.csv.gz"),header=T,stringsAsFactors=F,sep=",")
pi_pars = grep("pi", colnames(d1), value = T)
theta_pars = grep("theta", colnames(d1), value = T)
rho_pars = grep("rho", colnames(d1), value = T)

d1 = add_random_switching(d1)
d1$individual_fitness = pmax(pmax(d1$naive_fitness_analytical,0),d1$b-d1$c)
d1$fitness_gain = (d1$fitness - d1$individual_fitness)/(d1$b+d1$c)
d1$fitness_gain_relative_to_random_switching = (d1$fitness - pmax(d1$individual_fitness,d1$maximum_naive_fitness_analytical))/(d1$b+d1$c)

bs = sort(unique(d1$b))
lambdas = sort(unique(d1$alpha))
beta0 = 4


d2 = read.table("results/L3_processed.csv",header=T,stringsAsFactors=F,sep=",")

d_memory = d1[d1$beta == beta0,]
d_baseline = d2[d2$beta == beta0,]

ds = list(d_memory, d_baseline)

for(i in 1:2){
  ds[[i]]$unique_id <- match(do.call(paste, ds[[i]][c("lambda", "b")]), 
                             unique(do.call(paste, ds[[i]][c("lambda", "b")])))
  ds[[i]] = ds[[i]][order(ds[[i]]$unique_id, -ds[[i]]$fitness),]
  ds[[i]] = ds[[i]][!duplicated(ds[[i]]$unique_id),]
  ds[[i]] = ds[[i]][order(ds[[i]]$lambda,ds[[i]]$b),]
}

d = ds[[1]]
d0 = ds[[2]]


if(all(d$b == d0$b) && all(d$lambda == d0$lambda)){
  d$memory_fitness_gain = pmax(0,d$fitness - d0$fitness)/(d$b+d$c)
  d$interactions_fitness_gain = pmax(0,d0$fitness - d$maximum_naive_fitness_analytical)/(d$b+d$c)
}

d$x = log2(d$b)
d$y = log10(d$lambda)
d = d[abs(d$x) <= 2,]
d = d[abs(d$y) <= 1,]





# Swap Y and Z if necessary, to ensure that piAXY > piAXZ. 
swap_YZ <- function(x) {
  x <- gsub("Y", "__TEMP__", x)
  x <- gsub("Z", "Y", x)
  x <- gsub("__TEMP__", "Z", x)
  x
}
# Swap X and Z if necessary
swap_XZ <- function(x) {
  x <- gsub("X", "__TEMP__", x)
  x <- gsub("Z", "X", x)
  x <- gsub("__TEMP__", "Z", x)
  x
}
for(i in 1:nrow(d)){
  if(d$activityZ[i] == 1){
    if(d$piAXY[i] < d$piAXZ[i]){
      w = d[i, sapply(c(pi_pars,theta_pars,rho_pars),swap_YZ)]
      d[i,c(pi_pars,theta_pars,rho_pars)] = w
    }
  }
  if(d$activityZ[i] == 0){
    if(d$piBYX[i] < d$piBYZ[i]){
      w = d[i, sapply(c(pi_pars,theta_pars,rho_pars),swap_XZ)]
      d[i,c(pi_pars,theta_pars,rho_pars)] = w
    }
  }
}


pdf("plots/fig5_fitness_gains.pdf",width=4,height=4)
lambda_bc_plot(d$x, d$y, sapply(d$interactions_fitness_gain, function(x){fitness_increase_color(x)}))
lambda_bc_plot(d$x, d$y, sapply(d$memory_fitness_gain, function(x){fitness_increase_color(x)}))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
rect(0, 1/20, 0, 1, border=NA, col=fitness_increase_color(0))
for(i in 1:19){
  rect(i/20,0,(i+1)/20,1,border=NA,col=fitness_increase_color((i/20)*0.05+0.0001,max_value=0.05))
}
rect(0,0,1,1)
dev.off()

d = d[d$x != 0, ]

pdf("plots/figS5_memory_alleles.pdf",width=3,height=3)
for(p in c(pi_pars, theta_pars, rho_pars,"activityZ")){
  lambda_bc_plot_small(d$x, d$y, sapply(d[,p],trait_color),main=p,bty="n")
}
dev.off()




for(i in 1:nrow(d)){
  activity = c()
  activity["X"] = 0
  activity["Y"] = 1
  activity["Z"] = d$activityZ[i]
  d$thetaAA[i] = 0
  d$thetaAB[i] = 0
  d$thetaBA[i] = 0
  d$thetaBB[i] = 0
  d$rhoA[i] = 0
  d$rhoB[i] = 0
  d$Xinactivates[i] = 0
  d$Yinactivates[i] = 0
  d$Zinactivates[i] = 0
  d$Xactivates[i] = 0
  d$Yactivates[i] = 0
  d$Zactivates[i] = 0
  # calculate theta parameters
  for(p1 in c("X","Y","Z")){
    for(p2 in c("X","Y","Z")){
      for(p3 in c("X","Y","Z")){
        if(activity[p2] != activity[p3]){
          type = paste("theta",ifelse(activity[p1]==1,"A","B"),ifelse(activity[p3]==1,"A","B"),sep="")
          d[i, type] = d[i, type] + d[i, paste("theta",p1,p2,p3,sep="")]
          if(activity[p2] == 1 & d[i, paste("theta",p1,p2,p3,sep="")] > 0){ #inactivation
            d[i,paste(p1,"inactivates",sep="")] = 1
          }
          if(activity[p2] == 0 & d[i, paste("theta",p1,p2,p3,sep="")] > 0){ #activation
            d[i,paste(p1,"activates",sep="")] = 1
          }
        }
      }
    }
  }
  
  # calculate rho parameters
  for(p1 in c("X","Y","Z")){
    for(p2 in c("X","Y","Z")){
      if(activity[p1] != activity[p2]){
        type = paste("rho",ifelse(activity[p2]==1,"A","B"),sep="")
        d[i, type] = d[i, type] + d[i, paste("rho",p1,p2,sep="")]
      }
    }
  }
}

d$strategy_theta = sapply(1:nrow(d),function(i){
  call_strategy_theta_rho(d$thetaAA[i], d$thetaAB[i], d$thetaBA[i], d$thetaBB[i], d$rhoA[i], d$rhoB[i])
})



for(p in c("thetaAA","thetaAB","thetaBA","thetaBB","rhoA","rhoB")){
  lambda_bc_plot(d$x, d$y, intensityColor(d[,p],resolution=20,min=0,max=1),main=p)
}

pdf(width=4,height=4,"plots/figS6_memory_strategies.pdf")
lambda_bc_plot(d$x, d$y, strategy_colors[d$strategy_theta])
points(rep(1,2),c(log10(0.4),log10(4)),bg="white",pch=21,cex=2,lwd=2)
dev.off()





lambda1 = unique(d$lambda)[which.min(abs(unique(d$lambda)-0.4))]
lambda2 = unique(d$lambda)[which.min(abs(unique(d$lambda)-4))]

d[d$lambda == lambda1 & d$b == 2,]
d0[d0$lambda == lambda1 & d0$b == 2,]
