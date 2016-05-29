#Plot for monthly
setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x3")

x1<- read.table("frame23IS_1m.txt",T)
x2 <- read.table("frame2x3IS_1m_40y.txt",T)

fin <- as.numeric(x1[nrow(x1),4])
fin <- as.numeric(paste(format(fin/1000000,digits=4)))*10000+
  ((as.numeric(paste(format(fin/1000000,digits=6)))*10000-
      as.numeric(paste(format(fin/1000000,digits=4)))*10000)*100+1)/12

plot(y=x1[,5],x=seq(fin-(nrow(x1)*1/12)+1/12,to=fin,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="2x3 PFs, Monthly Returns",ylim=c(-30,30),xaxt="n",lwd=2,lty=1)
abline(0,0)
lines(y=x2[,5],x=seq(fin-(nrow(x2)*1/12)+1/12,to=fin,by=1/12),lty=2,col=2,lwd=2)
axis(1,xaxp=c(1940,2010,7),las=2)
legend("topleft",legend= c("1930-2015","1975-2015"),col=c(1,2),lty=c(1,2))

#Plot for yearly
setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x3")

x1<- read.table("frame23IS_1y.txt",T)
x2 <- read.table("frame2x3IS_1y_40y.txt",T)

fin <- as.numeric(x1[nrow(x1),4])
fin <- as.numeric(paste(format(fin/1000000,digits=4)))*10000+
  ((as.numeric(paste(format(fin/1000000,digits=6)))*10000-
      as.numeric(paste(format(fin/1000000,digits=4)))*10000)*100+1)/12

plot(y=x1[,5],x=seq(fin-(nrow(x1)*1/12)+1/12,to=fin,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="2x3 PFs, Yearly Returns",ylim=c(-30,30),xaxt="n",lwd=2,lty=1)
abline(0,0)
lines(y=x2[,5],x=seq(fin-(nrow(x2)*1/12)+1/12,to=fin,by=1/12),lty=2,col=2,lwd=2)
axis(1,xaxp=c(1940,2010,7),las=2)
legend("topleft",legend= c("1930-2015","1975-2015"),col=c(1,2),lty=c(1,2))

