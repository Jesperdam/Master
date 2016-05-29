setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/Multi")
library(xtable)
library(sandwich)
library(lmtest)

Factmont <- read.table("F-F_Research_Data_Factors.txt",skip=586,header=F,nrow=492)
y<-read.table("y.txt",T)
x<-read.table("x.txt",T,na.strings=c(NA))
y <- matrix(y[,8],ncol=1)
inf <- x$infl
x2 <- x[,2:101]
x <- x[,142:221]

framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y)*2)
for (l in 1:ncol(y)){
  prf <- rep(NA,nrow(y))
  prfm <- rep(NA,nrow(y))
  
  for (i in 60:(nrow(y)-12)) {
    tempr <- rep(NA,nrow(y))
    xoos <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    xoos2 <- matrix(NA,nrow=nrow(x2),ncol=ncol(x2))
    
    for (j in 1:ncol(xoos))
    {
      xoos[1:i,j] = x[1:i,j]
    }
    for (j in 1:ncol(xoos2))
    {
      xoos2[1:i,j] = x2[1:i,j]
    }
    tempr[1:(i-11)] <- y[1:(i-11),l]
    xoos <- xoos[ ,colSums(!is.na(xoos))>0.8*i]
    xoos2 <- xoos2[ ,colSums(!is.na(xoos2))>0.8*i]
    
    #stdize
    xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
    xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
    
    xoos_sd2 <- apply(xoos2, 2, function(col) sd(col, na.rm=TRUE))
    xoos2 <- t(apply(xoos2, 1, function(row) row/xoos_sd2))
    
    if (sum(!is.na(tempr))>60) {
      
      phi1 <- rep(NA,ncol(xoos))
      for (j in 1:ncol(xoos))
      {
        phi1[j] <- coef(lm(xoos[1:i,j]~tempr[2:(i+1)]))[2]
      }
      phi2 <- rep(NA,ncol(xoos2))
      for (j in 1:ncol(xoos2))
      {
        phi2[j] <- coef(lm(xoos2[1:i,j]~tempr[2:(i+1)]))[2]
      }
      
      #Cholesky decomp, for normalization?
      phi1<-phi1/chol(cov(phi1,phi1)*((length(phi1)-1)/(length(phi1))))
      phi2<-phi2/chol(cov(phi2,phi2)*((length(phi2)-1)/(length(phi2))))
      
      
      #second stage
      f1<- rep(NA,nrow(xoos))
      f2<- rep(NA,nrow(xoos2))
      
      
      for (j in 1:nrow(xoos))
      {
        if (sum(!is.na(xoos[j,]))>5)
        {
          f1[j] <- coef(lm(as.numeric(xoos[j,])~phi1))[2]
        }
      }
      
      for (j in 1:nrow(xoos2))
      {
        if (sum(!is.na(xoos2[j,]))>5)
        {
          f2[j] <- coef(lm(as.numeric(xoos2[j,])~phi2))[2]
        }
      }
      
      #final stage 1m
      coefs <- coef(lm(tempr[2:(nrow(xoos))]~f1[1:(nrow(xoos)-1)]+f2[1:(nrow(xoos2)-1)]+inf[1:(nrow(xoos)-1)]))
      
      #prediktor
      prf[i+1] <- coefs[1]+coefs[2]*f1[i]+coefs[3]*f2[i]+coefs[4]*inf[i]
      prfm[i+1] <- mean(tempr,na.rm=T)
      print(i)
    }
  }
  
  #OOSR2
  frame <- na.omit(cbind(y[,l],prf,prfm,Factmont[1:nrow(y),1]))
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
  framed[1:length(oosr2),(l*2-1):(l*2)] <- cbind(oosr2*100,pt(-abs(encnew),(length(oosr2)-1))*2)
}
oosres<- cbind(Factmont[73:(nrow(na.omit(framed))+72),1],na.omit(framed))

sims <- 100
#simulating inflation rate
arinf1 <- summary(lm(inf[12:492]~inf[11:491]+inf[1:481]))$coefficients[2,1]
arinf12 <- summary(lm(inf[12:492]~inf[11:491]+inf[1:481]))$coefficients[3,1]

means <- apply(x,2,FUN=mean,na.rm=T)
sds <- apply(x,2,FUN=sd,na.rm=T)
ars <- rep(1,ncol(x))
means2 <- apply(x2,2,FUN=mean,na.rm=T)
sds2 <- apply(x2,2,FUN=sd,na.rm=T)
ars2 <- rep(1,ncol(x2))
oos <- matrix(NA,ncol=sims,nrow=nrow(y))
pen <- matrix(NA,ncol=sims,nrow=nrow(y))

for (i in 1:ncol(x)){
  ars[i] <- as.numeric(coef(arima(x[i],order=c(1,0,0)))[1])
}

for (i in 1:ncol(x2)){
  ars2[i] <- as.numeric(coef(arima(x2[i],order=c(1,0,0)))[1])
}

df <- x
df2 <- x2

for (k in 1:sims){
  infsim <- arima.sim(model=list(ar=c(arinf1,0,0,0,0,0,0,0,0,0,0,arinf12)), n=492,mean=inf,sd=sd(inf))
    for (m in 1:ncol(x)){
        df[,m] <- arima.sim(n=nrow(df),model=list(ar=ars[m]),
              mean=means[m],sd=sds[m])
    }
  for (m in 1:ncol(x2)){
    df2[,m] <- arima.sim(n=nrow(df2),model=list(ar=ars2[m]),
                        mean=means2[m],sd=sds2[m])
  }
  framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y))
    for (l in 1:1){#ncol(y)){
      prf <- rep(NA,nrow(y))
      prfm <- rep(NA,nrow(y))
  
      for (i in 60:(nrow(y)-12)) {
      tempr <- rep(NA,nrow(y))
      xoos <- matrix(NA,nrow=nrow(df),ncol=ncol(df))
      for (j in 1:ncol(xoos))
      {
        xoos[1:i,j] = df[1:i,j]
      }
      xoos2 <- matrix(NA,nrow=nrow(df2),ncol=ncol(df2))
      for (j in 1:ncol(xoos2))
      {
        xoos2[1:i,j] = df2[1:i,j]
      }
      
      tempr[1:(i-11)] <- y[1:(i-11),l]
      xoos <- xoos[ ,colSums(!is.na(xoos))>0.6*i]
      xoos2 <- xoos2[ ,colSums(!is.na(xoos2))>0.6*i]
      
    
      #stdize
      xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
      xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
      xoos_sd2 <- apply(xoos2, 2, function(col) sd(col, na.rm=TRUE))
      xoos2 <- t(apply(xoos2, 1, function(row) row/xoos_sd2))
      
      if (sum(!is.na(tempr))>60+1) {
      
      phi1 <- rep(NA,ncol(xoos))
      for (j in 1:ncol(xoos))
      {
        phi1[j] <- coef(lm(xoos[1:i,j]~tempr[2:(i+1)]))[2]
      }
      
      phi2 <- rep(NA,ncol(xoos2))
      for (j in 1:ncol(xoos2))
      {
        phi2[j] <- coef(lm(xoos2[1:i,j]~tempr[2:(i+1)]))[2]
      }
      
      #Cholesky decomp, for normalization?
      phi1<-phi1/chol(cov(phi1,phi1)*((length(phi1)-1)/(length(phi1))))
      phi2<-phi2/chol(cov(phi2,phi2)*((length(phi2)-1)/(length(phi2))))
      
      #second stage
      f1<- rep(NA,nrow(xoos))
      f2<- rep(NA,nrow(xoos2))
      
      for (j in 1:nrow(xoos))
      {
        if (sum(!is.na(xoos[j,]))>5)
        {
          f1[j] <- coef(lm(as.numeric(xoos[j,])~phi1))[2]
        }
      }
      
      for (j in 1:nrow(xoos2))
      {
        if (sum(!is.na(xoos2[j,]))>5)
        {
          f2[j] <- coef(lm(as.numeric(xoos2[j,])~phi2))[2]
        }
      }
      
      
      #final stage 1m
      coefs <- coef(lm(tempr[2:(nrow(xoos))]~f1[1:(nrow(xoos)-1)]+f2[1:(nrow(xoos2)-1)]+infsim[1:(nrow(xoos)-1)]))
      
      #prediktor
      prf[i+1] <- coefs[1]+coefs[2]*f1[i]+coefs[3]*f2[i]+coefs[4]*infsim[i]
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
  oos[(72:(71+length(oosr2))),k] <-oosr2*100
  pen[(72:(71+length(oosr2))),k] <- pt(-abs(encnew),200)*2
  print(k)
}
}

write.table(cbind(oosres[1:nrow(na.omit(oos)),1:3],na.omit(oos)[,]),"asia_oos.txt",row.names=F)




