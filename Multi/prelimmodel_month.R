setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/Multi")
library(xtable)
library(sandwich)
library(lmtest)
library(stargazer)

y<-read.table("y.txt",T)
x<-read.table("x.txt",T,na.strings=c(NA))


#Insample analysis
xstd <- matrix(0,nrow=nrow(x),ncol=(ncol(x)-1))
for (i in 2:ncol(x))
{
  xstd[,i-1] <- x[,i]/sd(x[,i],na.rm=T)
}

#finding phis
fac <- matrix(1,ncol=15,nrow=(nrow(x)-1))
for (l in 1:5){
  tempy<- y[,1+l]
for (k in 1:3){
  if (k==1){
    tempx=xstd[,2:101]
    }
  else{
  if (k==2){
      tempx=xstd[,102:141]
  }
  else{
  if (k==3){
    tempx=xstd[,142:221]
  }
  }
  }
  tempx <- tempx[ , colSums(is.na(tempx))/nrow(tempx) < 0.4]
  
  phi <-0
  for (i in 1:(ncol(tempx))){
    phi[i] <- coef(lm(tempx[1:(nrow(tempx)-1),i]~tempy[2:nrow(tempx)]))[2]
  }
  phi<-phi/chol(cov(phi,phi)*((length(phi)-1)/(length(phi))))
  f1 <- 0
  
  for (i in 1:(nrow(tempx)-1)){
    fac[i,k*5-5+l]  <- coef(lm(tempx[i,]~phi))[2]
  }

}
}
beta1f<-matrix(0,nrow=1,ncol=15)
beta2f<-matrix(0,nrow=2,ncol=15)
beta3f<-matrix(0,nrow=3,ncol=5)
resmatrix <- matrix(0,ncol=35,nrow=(nrow(x)-1))
rsq <- rep(1,35)
beta <- cbind(rbind(beta1f,NA,NA),rbind(beta2f,NA),beta3f)
pkp <- beta
pkp[,]=0
pnw<-beta
pnw[,]<-0

for (l in 1:5)
{
for(k in 1:3)
{
  model <- lm(y[2:nrow(y),l+1]~fac[,l*3-3+k])
  beta1f[,l*3-3+k] <- coef(model)[2]
  resmatrix[,l*3-3+k] <- residuals(model)
  rsq[l*3-3+k] <- summary(model)$r.squared
  nyw <- coeftest(model, vcov.=NeweyWest(model, prewhite=F,lag=0))
  pnw[1,l*3-3+k] <- pt(-abs(nyw[2:dim(nyw)[1],3]),(nrow(y)-1))*2
  }
}
for (l in 1:5)
{
    model <- lm(y[2:nrow(y),l+1]~fac[,l*3-2]+fac[,l*3-1])
    beta2f[,l] <- coef(model)[2:3]
    resmatrix[,l+15] <- residuals(model)
    rsq[l+15] <- summary(model)$r.squared
    nyw <- coeftest(model, vcov.=NeweyWest(model, prewhite=F,lag=0))
    pnw[1:2,15+l] <- pt(-abs(nyw[2:dim(nyw)[1],3]),(nrow(y)-1))*2
    model <- lm(y[2:nrow(y),l+1]~fac[,l*3-2]+fac[,l*3])
    beta2f[,l+5] <- coef(model)[2:3]
    resmatrix[,l+20] <- residuals(model)
    rsq[l+20] <- summary(model)$r.squared
    nyw <- coeftest(model, vcov.=NeweyWest(model, prewhite=F,lag=0))
    pnw[1:2,l+20] <- pt(-abs(nyw[2:dim(nyw)[1],3]),(nrow(y)-1))*2
    model <- lm(y[2:nrow(y),l+1]~fac[,l*3]+fac[,l*3-1])
    beta2f[,l+10] <- coef(model)[2:3]
    resmatrix[,l+25] <- residuals(model)
    rsq[l+25] <- summary(model)$r.squared
    nyw <- coeftest(model, vcov.=NeweyWest(model, prewhite=F,lag=0))
    pnw[1:2,l+25] <- pt(-abs(nyw[2:dim(nyw)[1],3]),(nrow(y)-1))*2
    model <- lm(y[2:nrow(y),l+1]~fac[,l*3-2]+fac[,l*3-1]+fac[,l*3])
    beta3f[,l] <- coef(model)[2:4]
    resmatrix[,l+30] <- residuals(model)
    rsq[l+30] <- summary(model)$r.squared
    nyw <- coeftest(model, vcov.=NeweyWest(model, prewhite=F,lag=0))
    pnw[1:3,l+30] <- pt(-abs(nyw[2:dim(nyw)[1],3]),(nrow(y)-1))*2
}
beta <- cbind(rbind(beta1f,NA,NA),rbind(beta2f,NA),beta3f)
rsq<-rsq*100

resufunc <- function(x) rbind(cbind(beta[1,(x*3-2)],pnw[1,(x*3-2)],NA,NA,NA,NA,rsq[(x*3-2)]),
      cbind(NA,NA,beta[1,(x*3-1)],pnw[1,(x*3-1)],NA,NA,rsq[(x*3-1)]),
      cbind(NA,NA,NA,NA,beta[1,(x*3)],pnw[1,(x*3)],rsq[(x*3)]),
      cbind(beta[1,(x+15)],pnw[1,(x+15)],beta[2,(x+15)],pnw[2,(x+15)],NA,NA,rsq[(x+15)]),
      cbind(beta[1,(20+x)],pnw[1,(20+x)],NA,NA,beta[2,(20+x)],pnw[2,(20+x)],rsq[(20+x)]),
      cbind(NA,NA,beta[2,(25+x)],pnw[2,(25+x)],beta[1,(25+x)],pnw[1,(25+x)],rsq[(25+x)]),
      cbind(beta[1,(30+x)],pnw[1,(30+x)],beta[2,(30+x)],pnw[2,(30+x)],beta[3,(30+x)],pnw[3,(30+x)],rsq[(30+x)]))

results <- resufunc(5)
colnames(results) <- c("$\beta_{10x10}$","p(NW)","$\beta_{2x20}$","p(NW)","$\beta_{10x10}$","p(NW)","$IS-R^2$")
rownames(results) <- c("1","2","3","4","5","6","7")
xtable(results,digits=3)

#extra model investigation
relfac <- cbind(fac[,1],rep(0,nrow(fac)),fac[,4],fac[,6],fac[,8],fac[,9],fac[,10],fac[,11],fac[,15],rep(0,nrow(fac)))
tx1 <- relfac[1:(nrow(relfac)-1),9]
tx2 <- relfac[1:(nrow(relfac)-1),10]
ty <- y[2:nrow(fac),6]
fit1 <- lm(ty~tx1+x$caymon[1:(nrow(x)-2)])#lettauludvigson
nyw1 <- coeftest(fit1, vcov.=NeweyWest(fit1, prewhite=F,lag=0))[,2]
fit2 <- lm(ty~tx1+x$infl[1:(nrow(x)-2)])#goyalwelch
nyw2 <- coeftest(fit2, vcov.=NeweyWest(fit2, prewhite=F,lag=0))[,2]
fit3 <- lm(ty~tx1+y[1:(nrow(x)-2),6])#rapach
nyw3 <- coeftest(fit3, vcov.=NeweyWest(fit3, prewhite=F,lag=0))[,2]
fit4 <- lm(ty~tx1+y[1:(nrow(x)-2),4])#rapach
nyw4 <- coeftest(fit4, vcov.=NeweyWest(fit4, prewhite=F,lag=0))[,2]
fit5 <- lm(ty~tx1+x$corpr[1:(nrow(x)-2)])    
nyw5 <- coeftest(fit5, vcov.=NeweyWest(fit5, prewhite=F,lag=0))[,2]
fit6 <- lm(ty~tx1+x$caymon[1:(nrow(x)-2)]+x$infl[1:(nrow(x)-2)]+
             y[1:(nrow(x)-2),6]+y[1:(nrow(x)-2),4]+x$corpr[1:(nrow(x)-2)]+x$corpr[1:(nrow(x)-2)])

nyw6 <- coeftest(fit6, vcov.=NeweyWest(fit6, prewhite=F,lag=0))[,2]
  
stargazer(fit1,fit2,fit3,fit4,fit5,fit6,se=list(nyw1,nyw2,nyw3,nyw4,nyw5,nyw6),dep.var.labels = "Return on the World portfolio",
          covariate.labels = c("2x40 Factor","cay","US inflation","Return on US","Return on Europe","Corp rate US","Constant"),
          omit.stat=c("LL","ser","f","n","rsq","adj.rsq"),no.space=T,omit.table.layout = "n",
          table.placement = "H",add.lines=list(c("Adj. $R^2(\\%)$",round(summary(fit1)$adj.r.squared*100,3),
          round(summary(fit2)$adj.r.squared*100,3),round(summary(fit3)$adj.r.squared*100,3),
          round(summary(fit4)$adj.r.squared*100,3),round(summary(fit5)$adj.r.squared*100,3),
          round(summary(fit6)$adj.r.squared*100,3))))



  



