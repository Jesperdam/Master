library(sandwich)
Factmont<- read.table("F-F_Research_Data_Factors.txt",
                      skip=586,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed-Rf, divide by 100 and use logreturns: log(1+ret)
Mktrt <- log((Factmont[,2]-Factmont[,5])*0.01+1)

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
pkp <- (1-pt(tstatkp,summary(fit3)$df[2]))*2


#newey west
fit4 <- lm(y[2:(length(y))]~f1[1:(length(y)-1)])
nyw <- coeftest(fit4, vcov.=NeweyWest(fit4, prewhite=F,lag=0))
nywtstat <- nyw[2:dim(nyw)[1],3]
pnyw <- pt(-abs(nywtstat),summary(fit3)$df[2])*2

results <- cbind(summary(fit3)$r.squared,pkp,pnyw)
results



