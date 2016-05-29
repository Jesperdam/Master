#Copyright Jesper Dam
#Supported by Pruitt, Kelly's Market Return paper

#Section one, Investigate US returns on American Factors
#Load Y's Market return 1930-2015
library(sandwich)
library(lmtest)
library(xtable)

setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/5x5")
Factmont<- read.table("F-F_Research_Data_Factors.txt",
                      skip=4+36,
                      nrows=1076-36,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,2]-Factmont[,5])*0.01+1)

y <- Mktrt

#Load X's 10x10 size/BM 1930-2015
Firms <- read.table("25_Portfolios_5x5.txt",
                    skip=2371+36,
                    nrows=1076-36,
                    na.strings=c(-99.99,99.99,0),
                    header=FALSE)

Avgme<- read.table("25_Portfolios_5x5.txt",
                      skip=3452+36,
                      nrows = 1076-36,
                      na.strings=c(-99.99,99.99,0),
                      header=FALSE)

Sumbe <- read.table("25_Portfolios_5x5.txt",
                skip=4723+3,
                nrows = 90-3,
                na.strings=c(-99.99,99.99,0),
                header=FALSE)
Sumbe[Sumbe==0] <-NA
Sumbe <- Sumbe[,2:ncol(Sumbe)]

Summe <-  cbind(Firms[,2:ncol(Firms)]*Avgme[,2:ncol(Avgme)])
Summe[Summe==0] <-NA


#Construct monthly panel of log(Be/ME)
bm <- matrix(1,nrow=((nrow(Sumbe))*12),ncol=ncol(Summe))

for(i in 1:(nrow(Sumbe)-1))
{
  for (j in 1:12)
  {
  bm[(i*12-12+j),] <- as.numeric(Sumbe[i,1:ncol(Sumbe)])/as.numeric(Summe[(i*12-12+j),])
  }
}
bm<- bm[1:length(y),]
x <- log(bm)


#Remove obs with too many NA's
x <- x[ , colSums(is.na(x))/nrow(x) < 0.2]

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

results <- cbind(100*summary(fit3)$r.squared,pkp,pnyw)
xtable(results,digits=3)

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
  xoos <- xoos[ ,colSums(!is.na(xoos))>0.8*i]
  
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
     ylab="Out-of-sample R.squared", xlab="Time",main="5x5 Pfs",ylim=c(-30,30),xaxt="n")
abline(0,0)
axis(1,xaxp=c(1940,2010,7),las=2)
pencnew <- pt(-abs(encnew),(length(oosr2)-1))*2
write.table(cbind(frame,c(oosr2,rep(NA,60)),c(pencnew,rep(NA,60))),"frame5x5IS_1m.txt",
            row.names=F,col.names=c("Ret","Predict","HistMean","Year","OOSR2","pENCNEW"))



