setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/Multi")
library(ggplot2)
  x<-read.table("europe_oos.txt",header=T)
 #x<-read.table("world_oos.txt",header=T)
 #x<-read.table("asia_oos.txt",header=T)
# x<-read.table("scand_oos.txt",header=T)
#x<-read.table("us_oos.txt",header=T)



fin <- as.numeric(x[nrow(x),1])
fin <- as.numeric(paste(format(fin/1000000,digits=4)))*10000+
  ((as.numeric(paste(format(fin/1000000,digits=6)))*10000-
      as.numeric(paste(format(fin/1000000,digits=4)))*10000)*100-1)/12
date<-seq(fin-(nrow(x)*1/12)+1/12,to=fin,by=1/12)
dummy <- ifelse(x[,3]<0.05 & x[,2]>0,1,0)
x<-cbind(date,dummy,x)
subs<-subset(x,subset=dummy==1)

means<- x[,1]
sds<-x[,1]
for (i in 1:nrow(x)){
  means[i] <- apply(x[i,6:105],1,mean)
  sds[i] <- apply(x[i,6:105],1,sd)
}
upperone <- x[,1]
upperfive <- x[,1]
upperten <- x[,1]

for (i in 1:nrow(x)){
  upperone[i] <- means[i]+qt(0.995,99)*(sds[i]/sqrt(100))
  upperfive[i] <- means[i]+qt(0.975,99)*(sds[i]/sqrt(100))
  upperten[i] <- means[i]+qt(0.95,99)*(sds[i]/sqrt(100))
}

cols <- c("ENCNEW"="#f04546","LINE1"="#62c76b","99"="#042F6B","95"="#023F76","90"="#09689E")
ggplot(x,aes(x=x$date))+
  geom_rect(aes(xmin = subs[1,1], xmax = subs[nrow(subs),1],   ymin = 0, ymax = 10, fill ="ENCNEW"))+
  geom_ribbon(aes(ymax=upperone,ymin= -80, fill="99"))+
  geom_ribbon(aes(ymax=upperfive,ymin= -80, fill="95" ))+
  geom_ribbon(aes(ymax=upperten,ymin= -80, fill="90" ))+
  geom_line(aes(y=x$V2,color="LINE1"),size=2)+ylim(c(-80,10))+ 
  geom_line(y=0, size=1.5)+
  scale_colour_manual(name="OOSR2",values=cols)+
  scale_fill_manual(name=c("Areas"),values=cols) +
  labs(title = "Out-of-sample Predictive Power (Europe)",y="OOSR2",x="Time")+
  scale_x_continuous(expand = c(0,0))

