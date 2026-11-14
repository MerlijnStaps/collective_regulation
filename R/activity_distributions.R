source("R/plotting.R")

par(pty="m")

pdf("plots/activity_distributions.pdf",width=4,height=3.5)

for(motif in c("I","II","spontaneous")){
  d1 = read.table(paste("results/motif_",motif,"_examples.csv",sep=""),header=T,stringsAsFactors=F,sep=",")
  dm1 = read.table(paste("results/motif_",motif,"_examples_means.csv",sep=""),header=T,stringsAsFactors=F,sep=",")
  
  
  par(mgp=c(2,0.5,0))
  for(b0 in c(0,4)){
    d = d1[d1$beta == b0,]
    dm = dm1[dm1$beta == b0,]
    
    plot(0,0,type="n",xlim=c(0,1),ylim=c(0,10),bty="l",xlab="Activity",ylab="Density",xaxt="n",yaxt="n",
         main=bquote("Motif" ~ .(motif)*"," ~ lambda == .(l0) * "," ~ beta == .(b0)))
    axis(side=1, at=c(0,.2,.4,.6,.8,1), tck=-0.03, lwd =0, lwd.ticks=1)
    axis(side=2, at=c(0,5,10), tck=-0.03, las = 2, lwd =0, lwd.ticks=1)
    
    for(e in c("bad","good")){
      x = d$x[d$env == e]
      y = d$y[d$env == e]
      polygon(c(x,rev(x)),c(y,rep(0,length(y))),border=NA,col=add.alpha(environment_colors[e],0.3))
      points(x,y, type="l", lwd =3, col=environment_colors[e])
      segments(x0 = ifelse(e == "good", dm$meanA, dm$meanB), y0 = 0, y1 = 10, lwd=3, lty = 1, col =environment_colors[e])
    }
  }
  
  
  
}

dev.off()