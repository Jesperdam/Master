for (i in 1:15){
  tmp <- 0
  tct<-0
  fct<- NULL
  for (j in 1:nrow(resmatrix))
  {
    tmp1 <- resmatrix[j,i]^2*t(fac[j,i]-mean(fac[,i]))%*%(fac[j,i]-mean(fac[,i]))
    if (is.na(tmp1)) {
      tmp=tmp+0 
      tct=tct+0
      fct=fct
    }
    else {
      tmp=tmp+tmp1
      tct=tct+1
      fct= c(fct,fac[j,i])}
  } 
  tmp=tmp/tct
  avarbeta = ( tmp/(solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)) / 
    (solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)
  tstatkp <- beta[,i]/sqrt(diag(avarbeta)/tct)
  pkp[,i] <- pt(-abs(tstatkp),nrow(resmatrix)-1)*2
}
for (i in 1:5){
  for (l in 1:2){
    tmp <- 0
    tct<-0
    fct<- NULL
    for (j in 1:nrow(resmatrix))
    {
      tmp1 <- resmatrix[j,i+15]^2*t(fac[j,i*3-3+l]-mean(fac[,i*3-3-l]))%*%(fac[j,i*3-3+l]-mean(fac[,i*3-3+l]))
      if (is.na(tmp1)) {
        tmp=tmp+0 
        tct=tct+0
        fct=fct
      }
      else {
        tmp=tmp+tmp1
        tct=tct+1
        fct= c(fct,fac[j,i*3-3+l])}
    } 
    tmp=tmp/tct
    avarbeta = ( tmp/(solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)) / 
      (solve(tct)%*%t(fct)%*%(diag(tct)-(1/tct)*rep(1,tct))%*%fct)
    tstatkp <- beta[l,i+15]/sqrt(diag(avarbeta)/tct)
    pkp[l,i+15] <- pt(-abs(tstatkp),nrow(resmatrix)-2)*2
  }
}