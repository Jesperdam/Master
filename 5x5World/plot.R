#Plot for monthly
setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/5x5World")

x1<- read.table("oosy2005_y_bm_w.txt",T)
x2 <- read.table("oosy2005_y_bm_exus.txt",T)

fin <- as.numeric(x1[nrow(x1),1])
fin <- as.numeric(paste(format(fin/1000000,digits=4)))*10000+
  ((as.numeric(paste(format(fin/1000000,digits=6)))*10000-
      as.numeric(paste(format(fin/1000000,digits=4)))*10000)*100+1)/12

plot(y=x1[,2],x=seq(fin-(nrow(x1)*1/12)+1/12,to=fin,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="5x5 PFs, Yearly Returns",ylim=c(-8,8),xaxt="n",lwd=2,lty=1)
abline(0,0)
lines(y=x2[,2],x=seq(fin-(nrow(x2)*1/12)+1/12,to=fin,by=1/12),lty=2,col=2,lwd=2)
axis(1,xaxp=c(1997,2009,6),las=2)
legend("topleft",legend= c("World PF Predictor","Ex US Predictor"),col=c(1,2),lty=c(1,2))
