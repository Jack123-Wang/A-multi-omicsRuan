data<-read.table("all.id.txt",header=T)
rownames(data)<-data$Geneid
data<-data[7:14]

#calculate RPKM
genelength<-read.table("all.id.txt",header = T)
genelength<-genelength[c(1,6)]

neededGeneLength=genelength[match(rownames(data), genelength$Geneid),2] 

for(i in 1:1){
counts<-data[i]
effLen<-neededGeneLength
N<-sum(counts)
b0<-exp( log(counts) + log(1e9) - log(effLen) - log(N))
names(b0)[1]<-"id"
}

for(i in 2:ncol(data)){
  counts<-data[i]
  effLen<-neededGeneLength
  N<-sum(counts)
  b<-exp( log(counts) + log(1e9) - log(effLen) - log(N))
  names(b)[1]<-"id"
  b0<-cbind(b0,b$id)
}
colnames(b0)<-colnames(data)
rpkm<-b0
trpkm<-t(rpkm)
write.table(rpkm,"RPKM.txt")
