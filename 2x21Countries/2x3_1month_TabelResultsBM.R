#Copyright Jesper Dam
#Supported by Pruitt, Kelly's Market Return paper

#Section one, Investigate US returns on American Factors
#Load Y's Market return 1926-2016
library(sandwich)
library(lmtest)
library(xtable)

setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x21Countries")
Factmont<- read.table("Portfolios_Formed_on_BE-ME.txt",
                      skip=606,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,6:10])*0.01+1)

y <- Mktrt

bm <- read.table("Yrets.txt")
bm<- bm[1:nrow(y),2:length(bm)]
x <- log(bm/100+1)


#Remove obs with too many NA's
x <- x[ , colSums(is.na(x))/nrow(x) < 1]


#Recursive OOS
oos <- rep(1,ncol(y))
pen <- rep(1,ncol(y))
framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y)*2)
for (l in 1:ncol(y)){
prf <- rep(NA,nrow(y))
prfm <- rep(NA,nrow(y))

for (i in 60:(nrow(y)-1)) {
  tempr <- rep(NA,nrow(y))
  xoos <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  for (j in 1:ncol(xoos))
  {
    xoos[1:i,j] = x[1:i,j]
  }
  tempr[1:i] <- y[1:i,l]
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
framed[1:length(oosr2),(l*2-1):(l*2)] <- cbind(oosr2,pt(-abs(encnew),200)*2)
}
write.table(cbind(Factmont[61:(nrow(na.omit(framed))+60),1],na.omit(framed)),"oosy2005_1m.txt",
            row.names=F)

#yearly
#Copyright Jesper Dam


setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x21Countries")
Factmont<- read.table("Portfolios_Formed_on_BE-ME.txt",
                      skip=606,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,6:10])*0.01+1)


CSRets <- matrix(1,nrow=(nrow(Factmont)-11),ncol=ncol(Mktrt))
for(l in 1:ncol(CSRets)) {
for (i in 1:(nrow(Factmont)-11))
{
  CSRets[i,l] <- sum(Mktrt[i:(11+i),l])
}
}
y<-  CSRets

bm <- read.table("Yrets.txt")
bm<- bm[1:nrow(y),2:length(bm)]
x <- log(bm/100+1)


#Remove obs with too many NA's
x <- x[ , colSums(is.na(x))/nrow(x) < 1]


#Recursive OOS
oos <- rep(1,ncol(y))
pen <- rep(1,ncol(y))
framed <- matrix(NA,nrow=nrow(y),ncol=ncol(y)*2)
for (l in 1:ncol(y)){
  prf <- rep(NA,nrow(y))
  prfm <- rep(NA,nrow(y))
  
  for (i in 60:(nrow(y)-1)) {
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
  framed[1:length(oosr2),(l*2-1):(l*2)] <- cbind(oosr2,pt(-abs(encnew),200)*2)
}
write.table(cbind(Factmont[73:(nrow(na.omit(framed))+72),1],na.omit(framed)),"oosy2005_1y.txt",
            row.names=F)
res2 <- cbind(Factmont[73:(nrow(na.omit(framed))+72),1],na.omit(framed))







