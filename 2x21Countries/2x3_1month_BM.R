#Copyright Jesper Dam
#Supported by Pruitt, Kelly's Market Return paper

#Section one, Investigate US returns on American Factors
#Load Y's Market return 1926-2016
library(sandwich)
library(lmtest)
library(xtable)

setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/2x21Countries")
Factmont<- read.table("F-F_Research_Data_Factors.txt",
                      skip=586,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,2]-Factmont[,5])*0.01+1)

CSRets <- rep(1,(nrow(Factmont)-11))
for (i in 1:(nrow(Factmont)-11))
{
  CSRets[i] <- sum(Mktrt[i:(11+i)])
}
y <- CSRets[1:(length(CSRets)-1)]
y <- Mktrt 

bm <- read.table("BMrets.txt")
bm<- bm[1:length(y),2:length(bm)]
x <- log(bm/100+1)


#Remove obs with too many NA's
x <- x[ , colSums(is.na(x))/nrow(x) < 1]

#standardize x
xstd <- matrix(0,nrow=nrow(x),ncol=ncol(x))
for (i in 1:ncol(x))
{
  xstd[,i] <- x[,i]/sd(x[,i],na.rm=T)
}


#first stage
phi1 <- 0

for (i in 1:ncol(xstd))
{
  model <- summary(lm(xstd[1:(length(y)-1),i]~y[2:length(y)]))
  phi1[i] <- model$coefficients[2,1]
}

#second stage
f1<-0

for (i in 1:length(y))
{
  model <- summary(lm(as.numeric(xstd[i,])~phi1))
  f1[i] <- model$coefficients[2,1]
}

#final stage 1m
fit3 <- lm(y[2:length(y)]~f1[1:(length(y)-1)])
beta <- coef(fit3)[2]
res <- c(NA,resid(fit3))

#t-KP

#t-KP
#KP
tkp <- 0
tmp <- 0
tct<-0
fct<- NULL
for (j in 1:length(y))
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
fit4 <- lm(y[2:(length(y))]~f1[1:(length(y)-1)])
nyw <- coeftest(fit4, vcov.=NeweyWest(fit4, prewhite=F,lag=0))
nywtstat <- nyw[2:dim(nyw)[1],3]
pnyw <- pt(-abs(nywtstat),summary(fit3)$df[2])*2

results <- cbind(summary(fit3)$r.squared*100,pkp,pnyw)
xtable(results, digits=3)

#Recursive OOS
prf <- rep(NA,length(y))
prfm <- rep(NA,length(y))

for (i in 60:(length(y)-1)) {
  tempr <- rep(NA,length(y))
  xoos <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  for (j in 1:ncol(xoos))
  {
    xoos[1:i,j] = x[1:i,j]
  }
  tempr[1:i] <- y[1:i]
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
frame <- na.omit(cbind(y,prf,prfm,Factmont[1:length(y),1]))
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
oosr2<-oosr2*100
fin <- frame[length(oosr2),4]
fin<-as.numeric(paste(format(fin/1000000,digits=4)))*10000+
  ((as.numeric(paste(format(fin/1000000,digits=6)))*10000-
      as.numeric(paste(format(fin/1000000,digits=4)))*10000)*100+1)/12
plot(y=oosr2,x=seq(fin-(length(oosr2)*1/12)+1/12,to=fin,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="2x20 Pfs")
abline(0,0)
pencnew <- (1-pt(encnew,(length(oosr2)-1)))*2
write.table(cbind(frame,c(oosr2,rep(NA,60))),"frame1010BMrets_1m.txt",
            row.names=F)



