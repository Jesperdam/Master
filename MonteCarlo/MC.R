#Copyright Jesper Dam
#Supported by Pruitt, Kelly's Market Return paper

#Section one, Investigate US returns on American Factors
#Load Y's Market return 1926-2016
set.seed(99)
library(sandwich)
library(lmtest)
library(xtable)

setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/MonteCarlo")
Factmont<- read.table("Y.txt",
                      skip=586,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,2]-Factmont[,5])*0.01+1)

y <- matrix(Mktrt,nrow=length(Mktrt),ncol=1)

bm <- read.table("X.txt",
                na.strings=c(-99.99,99.99,NA),header=FALSE)
bm<- bm[1:nrow(y),2:length(bm)]
bm <- bm*0.01
#Remove obs with too many NA's

sims<-100
means <- apply(bm,2,FUN=mean,na.rm=T)
sds <- apply(bm,2,FUN=sd,na.rm=T)
ars <- rep(1,ncol(bm))
oos <- matrix(NA,ncol=sims,nrow=nrow(y))
pen <- matrix(NA,ncol=sims,nrow=nrow(y))

for (i in 1:ncol(bm)){
  ars[i] <- as.numeric(coef(arima(bm[i],order=c(1,0,0)))[1])
}

df <- bm
for (k in 1:sims){
  for (m in 1:ncol(bm)){
    df[,m] <- arima.sim(n=nrow(df),model=list(ar=ars[m]),
                        rand.gen=rnorm,mean=means[m],sd=sds[m])
  }

#Remove obs with too many NA's
df <- log(1+df)

#Recursive OOS
framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y))
for (l in 1:1){#ncol(y)){
prf <- rep(NA,nrow(y))
prfm <- rep(NA,nrow(y))

for (i in 60:(nrow(y)-1)) {
  tempr <- rep(NA,nrow(y))
  xoos <- matrix(NA,nrow=nrow(df),ncol=ncol(df))
  for (j in 1:ncol(xoos))
  {
    xoos[1:i,j] = df[1:i,j]
  }
  tempr[1:i] <- y[1:i,l]
  xoos <- xoos[ ,colSums(!is.na(xoos))>0.6*i]
  
  #stdize
  xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
  xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
  
  if (sum(!is.na(tempr))>60+1) {
    
    phi1 <- rep(NA,ncol(xoos))
    for (j in 1:ncol(xoos))
    {
      phi1[j] <- coef(lm(xoos[1:i,j]~tempr[2:(i+1)]))[2]
    }
    
    #Cholesky decomp, for normalization?
    phi1<-phi1/chol(cov(phi1,phi1)*((length(phi1)-1)/(length(phi1))))
    
    #second stage
    f1<- rep(NA,nrow(xoos))
    
    for (j in 1:nrow(xoos))
    {
      if (sum(!is.na(xoos[j,]))>5)
      {
        f1[j] <- coef(lm(as.numeric(xoos[j,])~phi1))[2]
      }
    }
    
    #final stage 1m
    coefs <- coef(lm(tempr[2:(nrow(xoos))]~f1[1:(nrow(xoos)-1)]))
    
    #prediktor
    prf[i+1] <- coefs[1]+coefs[2]*f1[i]
    prfm[i+1] <- mean(tempr,na.rm=T)
    print(i)
  }
}

#OOSR2
 frame <- na.omit(cbind(y[,l],prf,prfm))
 oosr2<- rep(NA,(nrow(frame)-60))
 encnew <- rep(NA,(nrow(frame)-60))
for (i in 1:(nrow(frame)-60))
{
  oosr2[i] = 1- sum( (frame[i:nrow(frame),1]-frame[i:nrow(frame),2])^2)/
    sum( ( frame[i:nrow(frame),1]-frame[i:nrow(frame),3])^2)
  fe1 <- frame[i:nrow(frame),1]-frame[i:nrow(frame),3]
  fe2 <- frame[i:nrow(frame),1]-frame[i:nrow(frame),2]
  encnew[i] <- sum(!is.na(fe1+fe2))*(sum(fe1^2-fe1*fe2))/sum(fe2^2)
}
oos[(61:(60+length(oosr2))),k] <-oosr2*100
pen[(61:(60+length(oosr2))),k] <- pt(-abs(encnew),100)*2
print(k)
}
}


#Find true OOS
framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y))
toos <- matrix(NA,ncol=1,nrow=nrow(y))
tpen <- matrix(NA,ncol=1,nrow=nrow(y))
for (l in 1:1){#ncol(y)){
  tprf <- rep(NA,nrow(y))
  tprfm <- rep(NA,nrow(y))
  
  for (i in 60:(nrow(y)-1)) {
    tempr <- rep(NA,nrow(y))
    xoos <- matrix(NA,nrow=nrow(df),ncol=ncol(df))
    for (j in 1:ncol(xoos))
    {
      xoos[1:i,j] = log(bm[1:i,j]+1)
    }
    tempr[1:i] <- y[1:i,l]
    xoos <- xoos[ ,colSums(!is.na(xoos))>0.6*i]
    
    #stdize
    xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
    xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
    
    if (sum(!is.na(tempr))>60+1) {
      
      phi1 <- rep(NA,ncol(xoos))
      for (j in 1:ncol(xoos))
      {
        phi1[j] <- coef(lm(xoos[1:i,j]~tempr[2:(i+1)]))[2]
      }
      
      #Cholesky decomp, for normalization?
      phi1<-phi1/chol(cov(phi1,phi1)*((length(phi1)-1)/(length(phi1))))
      
      #second stage
      f1<- rep(NA,nrow(xoos))
      
      for (j in 1:nrow(xoos))
      {
        if (sum(!is.na(xoos[j,]))>5)
        {
          f1[j] <- coef(lm(as.numeric(xoos[j,])~phi1))[2]
        }
      }
      
      #final stage 1m
      coefs <- coef(lm(tempr[2:(nrow(xoos))]~f1[1:(nrow(xoos)-1)]))
      
      #prediktor
      tprf[i+1] <- coefs[1]+coefs[2]*f1[i]
      tprfm[i+1] <- mean(tempr,na.rm=T)
      print(i)
    }
  }
  
  #OOSR2
  frame <- na.omit(cbind(y[,l],tprf,tprfm))
  oosr2<- rep(NA,(nrow(frame)-60))
  encnew <- rep(NA,(nrow(frame)-60))
  for (i in 1:(nrow(frame)-60))
  {
    oosr2[i] = 1- sum( (frame[i:nrow(frame),1]-frame[i:nrow(frame),2])^2)/
      sum( ( frame[i:nrow(frame),1]-frame[i:nrow(frame),3])^2)
    fe1 <- frame[i:nrow(frame),1]-frame[i:nrow(frame),3]
    fe2 <- frame[i:nrow(frame),1]-frame[i:nrow(frame),2]
    encnew[i] <- sum(!is.na(fe1+fe2))*(sum(fe1^2-fe1*fe2))/sum(fe2^2)
  }
  toos[(61:(60+length(oosr2)))] <-oosr2*100
  tpen[(61:(60+length(oosr2)))] <- pt(-abs(encnew),100)*2
}

test <- na.omit(cbind(toos,oos))
write.table(test,"montecarlores2.txt",row.names=F,col.names=F)



data <- read.table("montecarlores.txt")
plot(y=data[,1],x=seq(2012-(nrow(data)*1/12)+1/12,to=2012,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="2x3 Pfs",ylim=c(-30,4))

upper <-  rep(1,370)
for (i in 1:length(upper)){
  upper[i] <- mean(as.numeric(data[i,2:101]),na.rm=T)+1.96*(sd(as.numeric(data[i,2:101]),na.rm=T)/sqrt(100))
 }

abline(0,0)
lines(y=upper,x=seq(2012-(nrow(data)*1/12)+1/12,to=2012,by=1/12),col="red")

