#Copyright Jesper Dam
#Supported by Pruitt, Kelly's Market Return paper

#Section one, Investigate US returns on American Factors
#Load Y's Market return 1926-2016
library(sandwich)
library(lmtest)
library(xtable)

setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x3 World, with and ex US")
Factmont<- read.table("Asia_Pacific_ex_Japan_6_Portfolios_ME_BE-ME.txt",
                      skip=22,
                      nrows=306,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,2:7])*0.01+1)
CSRets <- matrix(1,nrow=(nrow(Factmont)-11),ncol=ncol(Mktrt))
for(l in 1:ncol(CSRets)) {
  for (i in 1:(nrow(Factmont)-11))
  {
    CSRets[i,l] <- sum(Mktrt[i:(11+i),l])
  }
}
y<-  CSRets
bm <- read.table("Yrets.txt",
                 skip=186,nrows=306,
                 na.strings=c(-99.99,99.99,NA),header=FALSE)
bm<- bm[1:nrow(y),2:length(bm)]
x <- log(bm/100+1)

x <- x[ , colSums(is.na(x))/nrow(x) < 0.2]

#standardize x
xstd <- matrix(0,nrow=nrow(x),ncol=ncol(x))
for (i in 1:ncol(x))
{
  xstd[,i] <- x[,i]/sd(x[,i],na.rm=T)
}
results <- matrix(0,ncol=3,nrow=ncol(y))

for (l in 1:ncol(y)) {
  phi1 <- 0
  
  for (i in 1:ncol(xstd))
  {
    model <- summary(lm(xstd[1:(nrow(y)-1),i]~y[2:nrow(y),l]))
    phi1[i] <- model$coefficients[2,1]
  }
  
  #second stage
  f1<-0
  
  for (i in 1:nrow(y))
  {
    model <- summary(lm(as.numeric(xstd[i,])~phi1))
    f1[i] <- model$coefficients[2,1]
  }
  
  #final stage 1m
  fit3 <- lm(y[2:nrow(y),l]~f1[1:(nrow(y)-1)])
  beta <- coef(fit3)[2]
  res <- c(NA,resid(fit3))
  
  #t-KP
  #KP
  tkp <- 0
  tmp <- 0
  tct<-0
  fct<- NULL
  for (j in 1:nrow(y))
  {
    tmp1 <- res[j]^2*t(f1[j]-mean(f1))%*%(f1[j]-mean(f1))
    if (is.na(tmp1)) {
      tmp=tmp+0 
      tct=tct+0
      fct=fct
    }
    else {
      tmp=tmp+tmp1
      tct=tct+1
      fct= c(fct,f1[j])}
  } 
  tmp=tmp/tct
  avarbeta = ( tmp/(solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)) / 
    (solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)
  tstatkp <- beta/sqrt(diag(avarbeta)/tct)
  pkp <- pt(-abs(tstatkp),summary(fit3)$df[2])*2
  
  
  #newey west
  fit4 <- lm(y[2:(nrow(y)),l]~f1[1:(nrow(y)-1)])
  nyw <- coeftest(fit4, vcov.=NeweyWest(fit4, prewhite=F,lag=0))
  nywtstat <- nyw[2:dim(nyw)[1],3]
  pnyw <- pt(-abs(nywtstat),summary(fit3)$df[2])*2
  
  results[l,] <- cbind(summary(fit3)$r.squared*100,pkp,pnyw)
  
}
xtable(results,round=3)

#Recursive OOS
oos <- rep(1,ncol(y))
pen <- rep(1,ncol(y))
framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y)*2)
for (l in 1:ncol(y)){
  prf <- rep(NA,nrow(y))
  prfm <- rep(NA,nrow(y))
  
  for (i in 60:(nrow(y)-11)) {
    tempr <- rep(NA,nrow(y))
    xoos <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    for (j in 1:ncol(xoos))
    {
      xoos[1:i,j] = x[1:i,j]
    }
    tempr[1:(i-11)] <- y[1:(i-11),l]
    xoos <- xoos[ ,colSums(!is.na(xoos))>0.6*i]
    
    #stdize
    xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
    xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
    
    if (sum(!is.na(tempr))>60) {
      
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
  oos[l] <-oosr2[114]*100
  pen[l] <- pt(-abs(encnew[114]),(114-1))*2
}
write.table(cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed)),"oos_asia_y_220.txt",
            row.names=F)
xtable(rbind(cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,2:3],
             cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,4:5],
             cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,6:7],
             cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,8:9],
             cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,10:11],
             cbind(Factmont[72:(nrow(na.omit(framed))+71),1],na.omit(framed))[115,12:13]),digits=3)


