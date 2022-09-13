data<-read.table('GBMLGG_RNAseq.txt',header = T,stringsAsFactors = F)

gene<-read.csv("gene1.csv")
genetotal<-as.data.frame(gene$Total)
colnames(genetotal)[1]<-"sample"
library(dplyr)
data<-inner_join(data,genetotal,by="sample")
colnames(data)[1]<-"Entrez.ID"

# Perform horizontal gene comparison processing on the expression profile
#step1: Restore count
data[,c(2:ncol(data))]<-2^data[,c(2:ncol(data))]

#step2:CPM
for(i in 2:ncol(data)){
  data[,i]<-((data[,i])*1000000)/sum(data[,i]) 
}
#step3
for (i in 2:nrow(data)){
  data[i,c(2:ncol(data))]<-data[i,c(2:ncol(data))]/(sum(data[i,c(2:ncol(data))])/(ncol(data)-1)) #centralized
  data[i,c(2:ncol(data))]<-log2(data[i,c(2:ncol(data))])
  data[i,c(2:ncol(data))]<-data[i,c(2:ncol(data))]-(sum(data[i,c(2:ncol(data))])/(ncol(data)-1)) #centralized
}

data_anno<-read.table("GBMLGG_annno.txt",sep="\t",header = T)
colnames(data_anno)[1]<-"sample"
data_anno<-data_anno[1]
library(dplyr)

#step2: Import the gene set
gene1<-read.table("ES exp1.txt")
gene1<-as.integer(gene1$V1)
gene1<-intersect(data$Entrez.ID,gene1)

gene<-read.csv("gene1.csv")
colnames(genetotal)[1]<-"sample"

gene2<-na.omit(gene$PAF) #PAF
gene2<-intersect(data$Entrez.ID,gene2)

gene3<-na.omit(gene$Myc)#Myc
gene3<-intersect(data$Entrez.ID,gene3)

gene4<-na.omit(gene$Core)#Core
gene4<-intersect(data$Entrez.ID,gene4)

gene5<-na.omit(gene$PCGF) #PCGF
gene5<-intersect(data$Entrez.ID,gene5)

gene6<-na.omit(gene$TBX) #TBX
gene6<-intersect(data$Entrez.ID,gene6)

gene7<-na.omit(gene$PRC) #PRC
gene7<-intersect(data$Entrez.ID,gene7)

#step3, 4: Extract genes with FC greater than or less than 1, and do hypergeometric distribution test.
N=20000
res<-as.data.frame(matrix(nrow=6,ncol=ncol(data)-1))
colnames(res)<-colnames(data)[2:703]
rownames(res)<-c("PAF","Myc","Core","PCGF","TBX","PRC")

genex<-gene7
j<-6
for (i in 2:ncol(data)){
  a<-data[,c(1,i)]
  colnames(a)<-c("V1","V2")
  b1<-filter(a,V2 > 1)
  c1<-b1$V1
  k1<-as.numeric(length(intersect(c1,genex)))
  if(length(c1)>length(genex)){M=as.numeric(length(c1));s=as.numeric(length(genex))}else{M=as.numeric(length(genex));s=as.numeric(length(c1))}
  
  L=N-M
  ke=s*M/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M, L, s))}    
  if(k1>ke){p1=sum(dhyper(k1:s, M, L, s))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  b2<-filter(a,V2< -1)
  c2<-b2$V1
  k2<-as.numeric(length(intersect(c2,genex)))
  if(length(c2)>length(genex)){M=as.numeric(length(c2));s=as.numeric(length(genex))}else{M=as.numeric(length(genex));s=as.numeric(length(c2))}
  
  L=N-M
  ke=s*M/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M, L, s))}    
  if(k2>ke){p2=sum(dhyper(k2:s, M, L, s))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  if(p1>0.05 &p2>0.05){res[j,i-1]=0}
  
  if(p1<0.05 &r1>0 & p2>0.05){res[j,i-1]<-  1}
  if(p1<0.05 &r1<0 & p2>0.05){res[j,i-1]=0}
  
  if(p2<0.05 &r2>0 & p1>0.05){res[j,i-1]<- -1}
  if(p2<0.05 &r2<0 & p1>0.05){res[j,i-1]=0}
  
  if(p1<0.05 & p2<0.05 & r1>0 & r1>r2){res[j,i-1]=1}
  if(p1<0.05 & p2<0.05 & r2>0 & r1<r2){res[j,i-1]=-1}
  
  if(p1<0.05 & p2<0.05 & r1<0 & r1>r2){res[j,i-1]=0}
  if(p1<0.05 & p2<0.05 & r2<0 & r1<r2){res[j,i-1]=0}  
}

tres<-as.data.frame(t(res))
tres$sample<-rownames(tres)

tres1<-inner_join(data_anno,tres,by="sample")
rownames(tres1)<-tres1$sample
res1<-t(tres1)
res1<-res1[-1,]
res1<-as.data.frame(res1)
write.csv(res1,"tmp.csv")
res1<-read.csv("tmp.csv")
rownames(res1)<-res1$X
res1<-res1[-1]

rownames(tres1)<-c(1:702)
library(pheatmap)
pheatmap(res1,
         scale = "none",
         cluster_row = FALSE,
         cluster_cols = FALSE,
         color=colorRampPalette(c("green","black","red"))(1000),
         border=T,
         border_color = "whtie",
         show_rownames=T,show_colnames=F,
         gaps_col = c(167,697))

library(reshape2)
res11<-res1
res11<-as.matrix(res11)
res12<-melt(res11)
res12$Var1<-factor(res12$Var1,levels = c("CTCF","PRC","TBX","PCGF","Core","Myc","PAF"),ordered=T)

res12$Var1<-factor(res12$Var1,levels = c("PCGF","TBX","PRC","PAF","Myc","Core"),ordered=T)

library(ggplot2)

red<-rgb(255,10,23,maxColorValue=255)
blue<-rgb(79,129,189,maxColorValue=255)
p<-ggplot(res12, aes(Var2, Var1))
p+  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = blue, mid = "white", high = red)+
  theme(     axis.ticks.x=element_line(color="black",size=0,lineend = 22), 
             axis.ticks.y=element_line(color="black",size=0,lineend = 22),
             axis.title.x = element_text(size = 0, face = "bold"), 
             axis.title.y = element_text(size = 0, face = "bold"),
             axis.text.x = element_text(size = 0, face = "bold" ),
             axis.text.y = element_text(size = 10, face = "bold" )
  ) 



#second picture
library(dplyr)

#step1: Prepare the result table
Res1<-as.data.frame(matrix(ncol=3,nrow=6))
Res2<-as.data.frame(matrix(ncol=3,nrow=6))
colnames(Res1)<-c("GBM","LGG","Normal")
colnames(Res2)<-c("GBM","LGG","Normal")
rownames(Res1)<-c("PAF","Myc","Core","PCGF","TBX","PRC")
rownames(Res2)<-c("PAF","Myc","Core","PCGF","TBX","PRC")

#step2: Prepare the overlap and calculate the p value
#1、Normal
tres1_normal<-tres1[c(1:167),]

#1-2、Normal-PAF
{
  a1<-as.numeric(nrow(filter(tres1,PAF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PAF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,PAF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,PAF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=1;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#1-3、Normal-Myc
{
  a1<-as.numeric(nrow(filter(tres1,Myc== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Myc==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,Myc== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,Myc==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=2;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#1-4、Normal-Core
{
  a1<-as.numeric(nrow(filter(tres1,Core== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Core==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,Core== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,Core==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=3;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  Res1[j,i]=1;Res2[j,i]=-log2(p2)
}
#1-5、Normal-PCGF
{
  a1<-as.numeric(nrow(filter(tres1,PCGF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PCGF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,PCGF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,PCGF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=4;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
  }
#1-6、Normal-TBX
{
  a1<-as.numeric(nrow(filter(tres1,TBX== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,TBX==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,TBX== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,TBX==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=5;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
  }
#1-7、Normal-PRC
{
  a1<-as.numeric(nrow(filter(tres1,PRC== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PRC==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_normal,PRC== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_normal,PRC==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_normal)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=6;i=1 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}


#2、ODG2
tres1_ODG2<-tres1[c(168:697),]

#2-2、Normal-PAF
{
  a1<-as.numeric(nrow(filter(tres1,PAF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PAF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,PAF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,PAF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=1;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#2-3、Normal-Myc
{
  a1<-as.numeric(nrow(filter(tres1,Myc== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Myc==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,Myc== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,Myc==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=2;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#2-4、Normal-Core
{
  a1<-as.numeric(nrow(filter(tres1,Core== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Core==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,Core== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,Core==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=3;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#2-5、Normal-PCGF
{
  a1<-as.numeric(nrow(filter(tres1,PCGF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PCGF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,PCGF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,PCGF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=4;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#2-6、Normal-TBX
{
  a1<-as.numeric(nrow(filter(tres1,TBX== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,TBX==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,TBX== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,TBX==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=5;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#2-7、Normal-PRC
{
  a1<-as.numeric(nrow(filter(tres1,PRC== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PRC==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_ODG2,PRC== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_ODG2,PRC==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_ODG2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=6;i=2 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}



#3、AC2
tres1_AC2<-tres1[c(698:703),]

#3-2、Normal-PAF
{
  a1<-as.numeric(nrow(filter(tres1,PAF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PAF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,PAF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,PAF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=1;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#3-3、Normal-Myc
{
  a1<-as.numeric(nrow(filter(tres1,Myc== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Myc==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,Myc== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,Myc==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=2;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#3-4、Normal-Core
{
  a1<-as.numeric(nrow(filter(tres1,Core== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,Core==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,Core== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,Core==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=3;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#3-5、Normal-PCGF
{
  a1<-as.numeric(nrow(filter(tres1,PCGF== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PCGF==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,PCGF== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,PCGF==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=4;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#3-6、Normal-TBX
{
  a1<-as.numeric(nrow(filter(tres1,TBX== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,TBX==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,TBX== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,TBX==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=5;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
}
#3-7、Normal-PRC
{
  a1<-as.numeric(nrow(filter(tres1,PRC== 1))) #over-total-ES sample
  a2<-as.numeric(nrow(filter(tres1,PRC==-1))) #under-total-ES sample
  
  b1<-as.numeric(nrow(filter(tres1_AC2,PRC== 1))) #over-normal-ES sample
  b2<-as.numeric(nrow(filter(tres1_AC2,PRC==-1))) #under-normal-ES sample
  
  c=as.numeric(nrow(tres1_AC2)) # number sample of group
  
  N=as.numeric(nrow(tres1)) #sample number
  

  
  if(a1>c){M1=a1;s1=c}else{M1=c;s1=a1} #over-a1和c
  if(a2>c){M2=a2;s2=c}else{M2=c;s2=a2} #under-a2和c
  k1= b1 #b intersect
  k2= b2 #b intersect
  
  L=N-M1
  ke=s1*M1/N
  if(k1<=ke){p1=sum(dhyper(0:k1, M1, L, s1))}    
  if(k1>ke){p1=sum(dhyper(k1:s1, M1, L, s1))}     
  if(k1<=ke){r1=-abs(log10(p1))}
  if(k1>ke){r1=abs(log10(p1))}
  
  L=N-M2
  ke=s2*M2/N
  if(k2<=ke){p2=sum(dhyper(0:k2, M2, L, s2))}    
  if(k2>ke){p2=sum(dhyper(k2:s2, M2, L, s2))}     
  if(k2<=ke){r2=-abs(log10(p2))}
  if(k2>ke){r2=abs(log10(p2))}
  
  j=6;i=3 
  if(p1<p2 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p2<p1 & r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1<p2 & r2<0 &r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1>p2 & r2<0 & r1<0 & r1<r2){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1<0 & r1>r2){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==1){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p2==1){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1>p2 & r2<0 & r1>0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1<p2 & r2>0 & r1<0){Res1[j,i]=-1;Res2[j,i]=log2(p2)}
  if(p1==p2 & r1>0 &r2<0){Res1[j,i]=1;Res2[j,i]=-log2(p1)}
  if(p1==p2 & r1<0 &r2>0){Res1[j,i]=-1;Res2[j,i]=log2(p1)}
  
}


library(pheatmap)
colnames(Res2)[1]<-c("GBM")
colnames(Res2)[2]<-c("LGG")
rownames(Res2)[1]<-c("PAF")
rownames(Res2)[2]<-c("MYC")

bk <- c(seq(-100,0,by=1),seq(0,100,by=1))

pheatmap(Res2,
         scale = "none",
         cluster_row = FALSE,
         cluster_cols = FALSE,
         color = c(colorRampPalette(colors = c("#4F81BD","#4F81BD","#E1E9F2"))(length(bk)/2+10),colorRampPalette(colors = c("#FCEFD3","#FFC14F","red","red"))(length(bk)/2)),
         border_color = "grey",
         border=T,
         show_rownames=T,show_colnames=T)

