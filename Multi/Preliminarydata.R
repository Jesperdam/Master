setwd("C:/Users/jesper/OneDrive/Master Thesis/Data testing/rstuff/Multi")

#load data
indall <- read.table("Ind_all.Dat",header=F,
                     skip=4,nrow=492)[2]
indasia <- read.table("Ind_Asia_Pacific.Dat",header=F,
                     skip=4,nrow=492)[2]
indeuro <- read.table("Ind_Eur_With_UK.Dat",header=F,
                     skip=4,nrow=492)[2]
induk <- read.table("Ind_UK.Dat",header=F,
                     skip=4,nrow=492)[2]
indscan <- read.table("Ind_Scandanavia.Dat",header=F,
                     skip=4,nrow=492)[2]
Factmont<- read.table("F-F_Research_Data_Factors.txt",
                      skip=586,
                      nrows=492,
                      na.strings=c(-99.99,99.99),
                      header=FALSE)

#Construct market returns Exceed+Rf, divide by 100 and use logreturns: log(1+ret)
indus <- (Factmont[,2]+Factmont[,5])
markets <- cbind(indall,indasia,indeuro,indscan,indus)
markets <- log(markets*0.01+1)

#cumulative returns
cummarkets <- matrix(0,ncol=ncol(markets),nrow=(nrow(markets)-11))

for (l in 1:ncol(markets)){
for (i in 1:(nrow(markets)-11))
{
  cummarkets[i,l] <- sum(markets[i:(11+i),l])
}
}

#100 us book to market
#Load X's 10x10 size/BM 1975-2015
Firms <- read.table("100_Portfolios_10x10.txt",
                    skip=2947,
                    nrows=498,
                    na.strings=c(-99.99,99.99,0),
                    header=FALSE)

Avgme<- read.table("100_Portfolios_10x10.txt",
                   skip=4028,
                   nrows = 498,
                   na.strings=c(-99.99,99.99,0),
                   header=FALSE)

Sumbe <- read.table("100_Portfolios_10x10.txt",
                    skip=4771,
                    nrows = 42,
                    na.strings=c(-99.99,99.99,0),
                    header=FALSE)
Sumbe[Sumbe==0] <-NA
Sumbe <- Sumbe[,2:ncol(Sumbe)]

Summe <-  cbind(Firms[,2:ncol(Firms)]*Avgme[,2:ncol(Avgme)])
Summe[Summe==0] <-NA


#Construct monthly panel of log(Be/ME)
bm <- matrix(0,nrow=((nrow(Sumbe))*12),ncol=ncol(Summe))

for(i in 1:nrow(Sumbe))
{
  for (j in 1:12)
  {
    bm[(i*12-12+j),] <- as.numeric(Sumbe[i,1:ncol(Sumbe)])/as.numeric(Summe[(i*12-12+j),])
  }
}
bm <- log(bm[7:498,])

#load 2x20 DP sort
dpsort <- read.table("Yrets.txt",na.strings="NA")

#create 4x20 thingie
bmsort <- read.table("BMrets.txt",na.strings="NA")
epsort <- read.table("EPrets.txt",na.strings="NA")
cepsort <- read.table("CEPrets.txt",na.strings="NA")

splits <- log(cbind(dpsort[,2:41],bmsort[,2:41],epsort[,2:41],cepsort[,2:41])*0.01+1)

splitdata <- matrix(0,ncol=80,nrow=nrow(splits))
for (i in 1:80)
{
  splitdata[,i] <- splits[,i*2-1]-splits[,i*2]
}

#goyal data
gwmon<- read.csv("gwmonth.csv",na.strings="NA",header=T,sep=";",dec=".")[1249:1740,3:18]

#construct monthly cay
gwqua<- read.csv("gwquart.csv",na.strings="NA",header=T,sep=";",dec=".")[417:580,]
caymon <- rep(0,492)
for (i in 1:nrow(gwqua))
{
  caymon[i*3-2] <- gwqua[i,10]
}
for (i in 1:nrow(gwqua))
{
  caymon[i*3-1] <- caymon[i*3-2]+(caymon[i*3-2]-caymon[(i+1)*3-2])*(1/3)
  caymon[i*3] <- caymon[i*3-1]+(caymon[i*3-2]-caymon[(i+1)*3-2])*(1/3)
}
caymon[(length(caymon)-1):length(caymon)]<-caymon[(length(caymon)-2)]

#extracting data
y <- cbind(Factmont[,1],markets,rbind(cummarkets,matrix(NA,ncol=5,nrow=11)))
colnames(y) <- c("date","world","asia","euro","scand","us","cumworld","cumasia","cumeuro","cumscand","cumus")
write.table(y,"y.txt",row.names=F)

x<- cbind(Factmont[,1],bm,dpsort,splitdata,gwmon,caymon)
write.table(x,"x.txt",row.names=F)
