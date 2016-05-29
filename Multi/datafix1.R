data<-read.table("rec.csv",sep=";",skip=11,F)
write.table(data,"rec.txt",col.names=F,row.names=F)

