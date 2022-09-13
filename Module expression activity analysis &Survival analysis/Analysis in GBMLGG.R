exprset<-read.table('GBMLGG_RNAseq.txt',header = T,stringsAsFactors = F)
sample<-read.table('GBMLGG_clinicalMatrix.txt',sep='\t',header = T,stringsAsFactors = F)
gene<-read.csv('human gene.csv',stringsAsFactors=F)

# Perform horizontal gene comparison processing on the expression profile
#step1:return count
exprset[,c(2:ncol(exprset))]<-2^exprset[,c(2:ncol(exprset))]

#step2:CPM
for(i in 2:ncol(exprset)){
  exprset[,i]<-((exprset[,i])*1000000)/sum(exprset[,i]) 
}

#step3: Select a specific row (the gene rows with 7 classes are extracted to reduce the computational load below)
totalgene<-unique(gene$human.gene)
totalgene<-as.data.frame(totalgene)
colnames(totalgene)<-"sample"
match(totalgene$sample,exprset$sample)
exprset<-exprset[match(totalgene$sample,exprset$sample),]
exprset<-na.omit(exprset)


##The number of samples in the expression matrix table is inconsistent with the number of samples in the sample table and needs to be processed

ref<-colnames(exprset)
ref<-as.data.frame(colnames(exprset)[-1])
colnames(ref)<-'sampleID'
sample<-sample[match(ref$sampleID,sample$sampleID),]

sample<-na.omit(sample)

unique(gene$type)
gene1<-gene[gene$type=='PCGF',]
gene2<-gene[gene$type=='PAF',]
gene3<-gene[gene$type=='TBX',]
gene4<-gene[gene$type=='PRC',]
gene5<-gene[gene$type=='Myc',]
gene6<-gene[gene$type=='Core',]
##Filter table by gene
exprset1<-exprset[match(gene1$human.gene,exprset$sample),] #PCGF
exprset1<-na.omit(exprset1)
exprset2<-exprset[match(gene2$human.gene,exprset$sample),] #PAF
exprset2<-na.omit(exprset2)
exprset3<-exprset[match(gene3$human.gene,exprset$sample),] #TBX
exprset3<-na.omit(exprset3)
exprset4<-exprset[match(gene4$human.gene,exprset$sample),] #PRC
exprset4<-na.omit(exprset4)
exprset5<-exprset[match(gene5$human.gene,exprset$sample),] #MYC
exprset5<-na.omit(exprset5)
exprset6<-exprset[match(gene6$human.gene,exprset$sample),] #Core
exprset6<-na.omit(exprset6)


unique(sample$type)
sample1<-sample[sample$type=='Normal',]
sample2<-sample[sample$type=='LGG',]
sample3<-sample[sample$type=='GBM',]

group1<-c('sample',sample1$sampleID) #normal
group2<-c('sample',sample2$sampleID) #LGG
group3<-c('sample',sample3$sampleID) #GBM


##Filter table according to gene and sample at the same time (3*7)
exprset11<-subset(exprset1, select=group1)
index11<-duplicated(exprset11$sample)
exprset11<-exprset11[!index11,]
rownames(exprset11)<-exprset11[ ,1]
exprset11<-exprset11[ ,-1]
mean(apply(exprset11,1,mean))
sd(apply(exprset11,1,mean))

exprset12<-subset(exprset1, select=group2)
index12<-duplicated(exprset12$sample)
exprset12<-exprset12[!index12,]
rownames(exprset12)<-exprset12[ ,1]
exprset12<-exprset12[ ,-1]
mean(apply(exprset12,1,mean))
sd(apply(exprset12,1,mean))

exprset13<-subset(exprset1, select=group3)
index13<-duplicated(exprset13$sample)
exprset13<-exprset13[!index13,]
rownames(exprset13)<-exprset13[ ,1]
exprset13<-exprset13[ ,-1]
mean(apply(exprset13,1,mean))
sd(apply(exprset13,1,mean))

exprset21<-subset(exprset2, select=group1)
index21<-duplicated(exprset21$sample)
exprset21<-exprset21[!index21,]
rownames(exprset21)<-exprset21[ ,1]
exprset21<-exprset21[ ,-1]
mean(apply(exprset21,1,mean))
sd(apply(exprset21,1,mean))

exprset22<-subset(exprset2, select=group2)
index22<-duplicated(exprset22$sample)
exprset22<-exprset22[!index22,]
rownames(exprset22)<-exprset22[ ,1]
exprset22<-exprset22[ ,-1]
mean(apply(exprset22,1,mean))
sd(apply(exprset22,1,mean))

exprset23<-subset(exprset2, select=group3)
index23<-duplicated(exprset23$sample)
exprset23<-exprset23[!index23,]
rownames(exprset23)<-exprset23[ ,1]
exprset23<-exprset23[ ,-1]
mean(apply(exprset23,1,mean))
sd(apply(exprset23,1,mean))

exprset31<-subset(exprset3, select=group1)
index31<-duplicated(exprset31$sample)
exprset31<-exprset31[!index31,]
rownames(exprset31)<-exprset31[ ,1]
exprset31<-exprset31[ ,-1]
mean(apply(exprset31,1,mean))
sd(apply(exprset31,1,mean))

exprset32<-subset(exprset3, select=group2)
index32<-duplicated(exprset32$sample)
exprset32<-exprset32[!index32,]
rownames(exprset32)<-exprset32[ ,1]
exprset32<-exprset32[ ,-1]
mean(apply(exprset32,1,mean))
sd(apply(exprset32,1,mean))

exprset33<-subset(exprset3, select=group3)
index33<-duplicated(exprset33$sample)
exprset33<-exprset33[!index33,]
rownames(exprset33)<-exprset33[ ,1]
exprset33<-exprset33[ ,-1]
mean(apply(exprset33,1,mean))
sd(apply(exprset33,1,mean))

exprset41<-subset(exprset4, select=group1)
index41<-duplicated(exprset41$sample)
exprset41<-exprset41[!index41,]
rownames(exprset41)<-exprset41[ ,1]
exprset41<-exprset41[ ,-1]
mean(apply(exprset41,1,mean))
sd(apply(exprset41,1,mean))

exprset42<-subset(exprset4, select=group2)
index42<-duplicated(exprset42$sample)
exprset42<-exprset42[!index42,]
rownames(exprset42)<-exprset42[ ,1]
exprset42<-exprset42[ ,-1]
mean(apply(exprset42,1,mean))
sd(apply(exprset42,1,mean))

exprset43<-subset(exprset4, select=group3)
index43<-duplicated(exprset43$sample)
exprset43<-exprset43[!index43,]
rownames(exprset43)<-exprset43[ ,1]
exprset43<-exprset43[ ,-1]
mean(apply(exprset43,1,mean))
sd(apply(exprset43,1,mean))

exprset51<-subset(exprset5, select=group1)
index51<-duplicated(exprset51$sample)
exprset51<-exprset51[!index51,]
rownames(exprset51)<-exprset51[ ,1]
exprset51<-exprset51[ ,-1]
mean(apply(exprset51,1,mean))
sd(apply(exprset51,1,mean))

exprset52<-subset(exprset5, select=group2)
index52<-duplicated(exprset52$sample)
exprset52<-exprset52[!index52,]
rownames(exprset52)<-exprset52[ ,1]
exprset52<-exprset52[ ,-1]
mean(apply(exprset52,1,mean))
sd(apply(exprset52,1,mean))

exprset53<-subset(exprset5, select=group3)
index53<-duplicated(exprset53$sample)
exprset53<-exprset53[!index53,]
rownames(exprset53)<-exprset53[ ,1]
exprset53<-exprset53[ ,-1]
mean(apply(exprset53,1,mean))
sd(apply(exprset53,1,mean))

exprset61<-subset(exprset6, select=group1)
index61<-duplicated(exprset61$sample)
exprset61<-exprset61[!index61,]
rownames(exprset61)<-exprset61[ ,1]
exprset61<-exprset61[ ,-1]
mean(apply(exprset61,1,mean))
sd(apply(exprset61,1,mean))

exprset62<-subset(exprset6, select=group2)
index62<-duplicated(exprset62$sample)
exprset62<-exprset62[!index62,]
rownames(exprset62)<-exprset62[ ,1]
exprset62<-exprset62[ ,-1]
mean(apply(exprset62,1,mean))
sd(apply(exprset62,1,mean))

exprset63<-subset(exprset6, select=group3)
index63<-duplicated(exprset63$sample)
exprset63<-exprset63[!index63,]
rownames(exprset63)<-exprset63[ ,1]
exprset63<-exprset63[ ,-1]
mean(apply(exprset63,1,mean))
sd(apply(exprset63,1,mean))

##Create the form required for the drawing
data_mean<-c(mean(apply(exprset11,2,mean)),mean(apply(exprset12,2,mean)),mean(apply(exprset13,2,mean)),mean(apply(exprset21,2,mean)),mean(apply(exprset22,2,mean)),mean(apply(exprset23,2,mean)),mean(apply(exprset31,2,mean)),mean(apply(exprset32,2,mean)),mean(apply(exprset33,2,mean)),mean(apply(exprset41,2,mean)),mean(apply(exprset42,2,mean)),mean(apply(exprset43,2,mean)),mean(apply(exprset51,2,mean)),mean(apply(exprset52,2,mean)),mean(apply(exprset53,2,mean)),mean(apply(exprset61,2,mean)),mean(apply(exprset62,2,mean)),mean(apply(exprset63,2,mean)))
data_gene<-c("PCGF","PCGF","PCGF","PAF","PAF","PAF","TBX","TBX","TBX","PRC","PRC","PRC","MYC","MYC","MYC","CORE","CORE","CORE")
data_sd<-c(sd(apply(exprset11,2,mean)),sd(apply(exprset12,2,mean)),sd(apply(exprset13,2,mean)),sd(apply(exprset21,2,mean)),sd(apply(exprset22,2,mean)),sd(apply(exprset23,2,mean)),sd(apply(exprset31,2,mean)),sd(apply(exprset32,2,mean)),sd(apply(exprset33,2,mean)),sd(apply(exprset41,2,mean)),sd(apply(exprset42,2,mean)),sd(apply(exprset43,2,mean)),sd(apply(exprset51,2,mean)),sd(apply(exprset52,2,mean)),sd(apply(exprset53,2,mean)),sd(apply(exprset61,2,mean)),sd(apply(exprset62,2,mean)),sd(apply(exprset63,2,mean)))
data<-data.frame('mean'=data_mean,'gene'=data_gene,'sd'=data_sd)
data$mean<-as.numeric(data$mean)
data$sd<-as.numeric(data$sd)
data$N<-c(rep(nrow(exprset1),3),rep(nrow(exprset2),3),rep(nrow(exprset3),3),rep(nrow(exprset4),3),rep(nrow(exprset5),3),rep(nrow(exprset6),3))
data$SEM<-data$sd/(data$N)^(1/2)
data$gene<-factor(data$gene,levels = c("PCGF", "PAF", "TBX", "CORE", "MYC", "PRC"))
data$group<-c(rep(c("Normal","LGG","GBM"),6))


library(ggplot2)
data$gene<-as.factor(data$gene)
data$gene<-factor(data$gene,levels =c("CORE","MYC","PAF","PRC","TBX","PCGF"))

#SEM
data$cut<-data$mean-data$SEM
data$add<-data$mean+data$SEM
data$group<-as.factor(data$group)
data$group<-factor(data$group,levels = c("Normal","LGG","GBM"),ordered = T)

red<-rgb(255,0,0,maxColorValue=255)
blue<-rgb(79,129,189,maxColorValue=255)
p<-ggplot(data,aes(x=gene,y=mean,fill=group))
p+geom_bar(color='black',stat = "identity", width=0.7, position="dodge")+ 
  scale_fill_manual(values=c(blue,"orange",red))+ 
  labs(x='',y='Average module expression')+
  geom_errorbar(data,mapping=aes(ymin=cut, ymax=add), width=0.3, color="black",position = position_dodge(0.7),size=0.2)+
  geom_hline(aes(yintercept=0))+ 
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#######Survival#########

mean(apply(exprset52,2,mean))
mean(apply(exprset52,1,mean))
library(dplyr)
library(survival)
library(survminer)

#MYC
  data5_mean<-as.data.frame(apply(exprset5[2:703],2,mean))
  colnames(data5_mean)<-"FC"
  data5_mean<-arrange(data5_mean,desc(FC))
  data5_mean$sample<-rownames(data5_mean)
  top<-data5_mean[c(1:(0.5*nrow(data5_mean))),2]
  down<-data5_mean[c((0.5*nrow(data5_mean)):nrow(data5_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"myc.csv")
  
  
#PRC
  data4_mean<-as.data.frame(apply(exprset4[2:703],2,mean))
  colnames(data4_mean)<-"FC"
  data4_mean<-arrange(data4_mean,desc(FC))
  data4_mean$sample<-rownames(data4_mean)
  top<-data4_mean[c(1:(0.5*nrow(data4_mean))),2]
  down<-data4_mean[c((0.5*nrow(data4_mean)):nrow(data4_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"prc.csv")

#CORE
  data6_mean<-as.data.frame(apply(exprset6[2:703],2,mean))
  colnames(data6_mean)<-"FC"
  data6_mean<-arrange(data6_mean,desc(FC))
  data6_mean$sample<-rownames(data6_mean)
  top<-data6_mean[c(1:(0.5*nrow(data6_mean))),2]
  down<-data6_mean[c((0.5*nrow(data6_mean)):nrow(data6_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"core.csv")
  
#PCGF
  data1_mean<-as.data.frame(apply(exprset1[2:703],2,mean))
  colnames(data1_mean)<-"FC"
  data1_mean<-arrange(data1_mean,desc(FC))
  data1_mean$sample<-rownames(data1_mean)
  top<-data1_mean[c(1:(0.5*nrow(data1_mean))),2]
  down<-data1_mean[c((0.5*nrow(data1_mean)):nrow(data1_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"n1.csv")
  
  
#PAF
  data2_mean<-as.data.frame(apply(exprset2[2:703],2,mean))
  colnames(data2_mean)<-"FC"
  data2_mean<-arrange(data2_mean,desc(FC))
  data2_mean$sample<-rownames(data2_mean)
  top<-data2_mean[c(1:(0.5*nrow(data2_mean))),2]
  down<-data2_mean[c((0.5*nrow(data2_mean)):nrow(data2_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"n2.csv")
  
  
#TBX
  data3_mean<-as.data.frame(apply(exprset3[2:703],2,mean))
  colnames(data3_mean)<-"FC"
  data3_mean<-arrange(data3_mean,desc(FC))
  data3_mean$sample<-rownames(data3_mean)
  top<-data3_mean[c(1:(0.5*nrow(data3_mean))),2]
  down<-data3_mean[c((0.5*nrow(data3_mean)):nrow(data3_mean)),2]
  write.csv(rbind(cbind(top,"top"),cbind(down,"down")),"n3.csv")
  



library(dplyr)
{
  top<-as.data.frame(top)
  top$type<-"top"
  colnames(top)[1]<-"sample"
  down<-as.data.frame(down)
  colnames(down)[1]<-"sample"
  down$type<-"down"
  total<-rbind(top,down)
  survival<-read.table("GBMLGG_survival.txt",sep="\t",header = T)
  total<-inner_join(total,survival,by="sample")
}


fit<-survfit(Surv(PFI.time, PFI) ~ type, data = total)

print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#3D6FB5", "#EB3324")
)



