setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/TPRF-rep")
library(sandwich)
library(lmtest)
x <- read.table("6bm.txt",sep=";")
y <- read.table("y1y.txt",sep=";")

#standardize x
xstd <- matrix(0,nrow=nrow(x),ncol=ncol(x))
for (i in 1:ncol(x))
{
  xstd[,i] <- x[,i]/sd(x[,i])
}


#first stage
phi1 <- 0

for (i in 1:ncol(x))
{
  model <- summary(lm(xstd[1:(nrow(x)-1),i]~y[2:nrow(x),1]))
  phi1[i] <- model$coefficients[2,1]
}

#second stage
f1<-0

for (i in 1:nrow(x))
{
  model <- summary(lm(as.numeric(xstd[i,])~phi1))
  f1[i] <- model$coefficients[2,1]
}

#final stage 1m
fit3 <- lm(y[2:nrow(x),1]~f1[1:(nrow(x)-1)])
summary(fit3)
beta <- coef(fit3)[2]
res <- c(NA,resid(fit3))

#t-KP
#KP
tkp <- rep(1,12)

for (i in 1:12)
{
  tmp <- 0
  tct<-0
  fct<- NULL
  for (j in seq(i,nrow(x),by=12))
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
  tkp[i] <- beta/sqrt(diag(avarbeta)/tct)
}
tstatkp <- median(tkp)
pkp <- pt(-abs(tstatkp),summary(fit3)$df[2])*2


#newey west
fit4 <- lm(y[2:(nrow(y)-11),1]~f1[1:(nrow(x)-12)])
nyw <- coeftest(fit4, vcov.=NeweyWest(fit4, prewhite=F,lag=11))
nywtstat <- nyw[2:dim(nyw)[1],3]
pnyw <- (1-pt(nywtstat,summary(fit3)$df[2]))*2

results <- cbind(summary(fit3)$r.squared,pkp,pnyw)

#Recursive OOS
prf <- rep(NA,nrow(x))
prfm <- rep(NA,nrow(x))

for (i in 60:(nrow(x)-11)) {
  tempr <- rep(NA,nrow(y))
  xoos <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  for (j in 1:ncol(xoos))
  {
  xoos[1:i,j] = x[1:i,j]
  }
  tempr[1:(i-11)] <- y[1:(i-11),1]
  xoos <- xoos[ ,colSums(!is.na(xoos))>0.8*i]
  
  #stdize
  xoos_sd <- apply(xoos, 2, function(col) sd(col, na.rm=TRUE))
  xoos <- t(apply(xoos, 1, function(row) row/xoos_sd))
  
  if (sum(!is.na(tempr))>60) {
  
  phi1 <- rep(NA,ncol(xoos))
  for (j in 1:ncol(xoos))
  {
    model <- summary(lm(xoos[1:i,j]~tempr[2:(i+1)]))
    phi1[j] <- model$coefficients[2,1]
  }
  
  #Cholesky decomp, for normalization?
  phi1<-phi1/chol(cov(phi1,phi1)*((length(phi1)-1)/(length(phi1))))
  
  #second stage
  f1<- rep(NA,nrow(xoos))
  
  for (j in 1:nrow(xoos))
  {
    if (sum(!is.na(xoos[j,]))>5)
    {
    model <- summary(lm(as.numeric(xoos[j,])~phi1))
    f1[j] <- model$coefficients[2,1]
    }
  }
  
  #final stage 1m
  fit3 <- lm(tempr[2:(nrow(xoos))]~f1[1:(nrow(xoos)-1)])

  #prediktor
  prf[i+1] <- coef(fit3)[1]+coef(fit3)[2]*f1[i]
  prfm[i+1] <- mean(tempr,na.rm=T)
  print(i)
  }
}

#OOSR2
frame <- na.omit(cbind(y[,],prf,prfm))
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
oosr2*100
fin <- 2010
plot(y=oosr2,x=seq(fin-(length(oosr2)*1/12)+1/12,to=fin,by=1/12),type="l",
     ylab="Out-of-sample R.squared", xlab="Time",main="10x10 Pfs")
pencnew <- (1-pt(encnew,(length(oosr2)-1)))*2



