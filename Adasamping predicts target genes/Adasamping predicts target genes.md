## Adasamping predicts target genes

###  1. Determine the proximal gene of the class and the distal gene corresponding to Hic

Randomly sample the same number of genes as the target gene set in R, as the non-target gene set, assign 1 and 0 respectively.

### 2. Prepare a feature matrix and remove genes without these features, which are generally predicted genes, and try not to lose target genes.

①EB0d, 4d, 6d, 9d, 12d, 20d RNA-seq data (unit Z-score, horizontal comparison, Z-score = single time point of the same gene minus the average of all time points of the same gene/all time points of the same gene standard deviation)

```
setwd("/Users/apple/Desktop/Adasampling/feature总表")
data<-read.csv("GSE120224_TPM.csv")
data<-filter(data,EB0d!=0 |EB4d!=0|EB6d!=0|EB9d!=0|EB12d!=0|EB20d!=0)
data<-data[c(1:7)]

for (i in 1:nrow(data)){
  data[i,c(2:ncol(data))]<-data[i,c(2:ncol(data))]/(sum(data[i,c(2:ncol(data))])/(ncol(data)-1)) #标准化
  data[i,c(2:ncol(data))]<-log2(data[i,c(2:ncol(data))]+0.001)
  data[i,c(2:ncol(data))]<-data[i,c(2:ncol(data))]-(sum(data[i,c(2:ncol(data))])/(ncol(data)-1)) #中心化
}

library(reshape2)
boxplot<-melt(data[,2:7])
p<-ggplot()
p+geom_boxplot(boxplot,mapping=aes(x=variable,y=value))

write.csv(data,"GSE120224_TPM_R.csv")
```



② TPM data of 10 histone proteins before and after the gene promoter of 5kb (unit log2TPM)

```bash
#step1：bam、bai、bed files preparation
GSE11724_H3K79me2.merged.positionsort.markdup.bam      GSE29218_H3K4me3.merged.positionsort.markdup.bam.bai
GSE11724_H3K79me2.merged.positionsort.markdup.bam.bai  GSE29413_H3K9me3.merged.positionsort.markdup.bam
GSE12241_H4K20me3.merged.positionsort.markdup.bam      GSE29413_H3K9me3.merged.positionsort.markdup.bam.bai
GSE12241_H4K20me3.merged.positionsort.markdup.bam.bai  GSE30203_H3K4me1.merged.positionsort.markdup.bam
GSE24164_H3K27ac.merged.positionsort.markdup.bam       GSE30203_H3K4me1.merged.positionsort.markdup.bam.bai
GSE24164_H3K27ac.merged.positionsort.markdup.bam.bai   GSE31284_H3K9ac.merged.positionsort.markdup.bam
GSE25532_H3K27me3.merged.positionsort.markdup.bam      GSE31284_H3K9ac.merged.positionsort.markdup.bam.bai
GSE25532_H3K27me3.merged.positionsort.markdup.bam.bai  GSE41589_H3K36me3.merged.positionsort.markdup.bam
GSE27827_H3K4me2.merged.positionsort.markdup.bam       GSE41589_H3K36me3.merged.positionsort.markdup.bam.bai
GSE27827_H3K4me2.merged.positionsort.markdup.bam.bai   mm10_tss10kb.bed
GSE29218_H3K4me3.merged.positionsort.markdup.bam

#step2：bedtools multicov
bedtools multicov -bams GSE25532_H3K27me3.merged.positionsort.markdup.bam GSE29413_H3K9me3.merged.positionsort.markdup.bam  GSE12241_H4K20me3.merged.positionsort.markdup.bam GSE29218_H3K4me3.merged.positionsort.markdup.bam GSE27827_H3K4me2.merged.positionsort.markdup.bam GSE31284_H3K9ac.merged.positionsort.markdup.bam GSE41589_H3K36me3.merged.positionsort.markdup.bam GSE11724_H3K79me2.merged.positionsort.markdup.bam GSE30203_H3K4me1.merged.positionsort.markdup.bam GSE24164_H3K27ac.merged.positionsort.markdup.bam -bed mm10_tss10kb.bed



```



Convert to TPM

```R
setwd("~/Data_HD_3/ry_chipseq/Hic/test/Hic-test/step5")
data<-read.csv("result.txt",header = F,sep="\t")
colnames(data)<-c("chr","start","end","genename","H3K27me3","H3K9me3","H4K20me3","H3K4me3","H3K4me2","H3K9ac","H3K36me3","H3K79me2","H3K4me1","H3K27ac")

#Step1: Prepare expression matrix and gene length list
#Import gene length list
gene<-data[c(4,3)]
gene$end<-10000
names(gene)<-c("ID","length")

#Import the expression matrix (counts) (the final form is: rownames is the gene name, colnames is the sample name, starting from the first column is the value, that is, the expression matrix)
exprSet<-data[c(5:14)]
rownames(exprSet)<-data$genename
exprSet<-exprSet[ rownames(exprSet) %in% gene$ID ,]#得出与gene长度list对应的exprSet
#exprSet<-exprSet[,-1]

total_count<- colSums(exprSet)

#Get the gene length corresponding to exprSet
neededGeneLength=gene[match(rownames(exprSet), gene$ID),2] 


#Step2:formula
countToTpm <- function(counts, effLen)
{ rate <- log(counts) - log(effLen)
denom <- log(sum(exp(rate)))
exp(rate - denom + log(1e6))}

############   Calculate count to TPM   ############
rm(b)
for (i in 1:1){
  i=1
  counts<-exprSet[i]
  effLen<-neededGeneLength
  N<-sum(counts)
  b<-countToTpm(counts, effLen)
}

for (i in 2:length(exprSet)){
  counts<-exprSet[i]
  effLen<-neededGeneLength
  N<-sum(counts)
  c<-countToTpm(counts, effLen)
  names(c)[1]<-"id"
  b<-cbind(b,c$id)
}
names(b)<-names(exprSet)
write.csv(b,"tss_histone10_TPM.csv")
```



### 3. Target gene prediction of each sub-classes

Determine the truncated p value of the initial training result of the transcriptome: the proximal target gene and the distal target gene are regarded as positive, and the random gene set with the same number as the training data set is found as negative, and the FPR corresponding to different truncated p values is calculated, and finally p =1 screening threshold.

①Core's predictive analysis (other sub-classes analysis is similar, no longer presented)

step1: Get Core_inter.txt

```bash
cd /disk_HD_3/wangjiaqi/ry_chipseq/Hic/test

bedtools intersect -a ESC-Hic-promoter.bed -b core.2939.bed |cut -f 4|grep EN|sort -u >core_inter_pre1.txt

cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $1}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $2}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $3}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $4}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $5}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $6}'|grep -v "^$" >>core_inter_pre2.txt
cat core_inter_pre1.txt |sed 's/|/\t/g'|awk '{print $7}'|grep -v "^$" >>core_inter_pre2.txt

cat core_inter_pre2.txt|sort -u|wc #922
cat core_inter_pre2.txt|sort -u >core_inter.txt
rm core_inter_pre2.txt core_inter_pre1.txt
```

Step2: Merge proximal and distal target genes

```bash
cat core_inter.txt core_prox.txt|sort -u >core_targetgene.csv
```

Step3：Adasampling

```R
setwd("/Users/apple/Desktop/Adasampling/MYC")
#step1Importing data
#1. Prepare feature feature data
RNAseq<-read.csv("feature1Z：EBRNA-seq.csv")
colnames(RNAseq)[1]<-"ID"

Histone<-read.csv("feature2logTPM：tss_histone10_TPM.csv")
colnames(Histone)[1]<-"ID"

#2. Import target genes and random genes
library(dplyr)
data.target<-read.csv("Core_targetgene.csv")
colnames(data.target)[1]<-"ID"
data.target<-inner_join(data.target,RNAseq,by="ID")
data.target<-inner_join(data.target,Histone,by="ID")

totalgene<-read.csv("totalgene.csv",header = F)
data.random <-sample(setdiff(totalgene$V1,data.target$Name),nrow(data.target))
data.random<-as.data.frame(data.random)
colnames(data.random)<-"ID"
data.random<-inner_join(data.random,RNAseq,by="ID")
data.random<-inner_join(data.random,Histone,by="ID")

#step2: Prepare the matrix
data.mat <- apply(X = data.target[-c(1)], MARGIN = 2, FUN = as.numeric)
random.mat<-apply(X = data.random[-c(1)], MARGIN = 2, FUN = as.numeric)
rownames(data.mat)<- data.target$ID
rownames(random.mat)<- data.random$ID

# step3:Identify positive and negative examples from the noisy dataset
Ps <- rownames(data.mat)
Ns <- rownames(random.mat)

total.mat<-rbind(data.mat,random.mat)
#step4:Adasampling
library(AdaSampling)
data.preds <- adaSample(Ps, Ns, train.mat=total.mat, test.mat=total.mat, classifier = "knn")

library(dplyr)
#FPR calculation
a<-sort(unique(data.preds[,1]))
for (i in 1:length(a)){
  i=1
  b<-a[i]
  data.p<-filter(as.data.frame(data.preds),P>=b)
  FP<-length(intersect(rownames(data.p),Ns))
  FPR<-FP/length(Ns)
  c<-cbind(b,FPR)
}
for (i in 2:length(a)){
  b<-a[i]
  data.p<-filter(as.data.frame(data.preds),P>=b)
  FP<-length(intersect(rownames(data.p),Ns))
  FPR<-FP/length(Ns)
  c1<-cbind(b,FPR)
  c<-rbind(c,c1)
}
colnames(c)<-c("Probablility","FPR");c<-as.data.frame(c)

#FPR-Probability Scatter Plot Determining P Values
library(ggplot2)
p<-ggplot()
p+theme( axis.ticks.x=element_line(color="black",size=0.7,lineend = 22),  
         axis.ticks.y=element_line(color="black",size=0.7,lineend = 22),
         axis.title.x = element_text(size = 15, face = "bold"), 
         axis.title.y = element_text(size = 15, face = "bold"),
         axis.text.x = element_text(size = 10, face = "bold" ),
         axis.text.y = element_text(size = 10, face = "bold" ))+
  geom_point(c,mapping=aes(x=Probablility,y=FPR)) + xlim(1,0)+ylim(0,1)+
  labs( x="AdaEnsemble prediction threshold 
        (probability positive)", y="Estimated false positive rate(FDR)")+
  theme(panel.grid.minor=element_blank())

#According to the above, to ensure that the FPR is less than 10%, select the p value, so choose p=1
data.p1<-filter(as.data.frame(data.preds),P>=1)
length(intersect(rownames(data.p1),Ps))
TP<-intersect(rownames(data.p1),Ps)

#Draw the expression profile to verify the trend of expression changes
targetgene<-as.data.frame(intersect(rownames(data.p1),Ps))
colnames(targetgene)<-"ID"
targetgene<-inner_join(targetgene,RNAseq)

targetgene<-targetgene[-c(1,8,9)]
library(reshape2)
boxplot<-melt(targetgene)
p<-ggplot()
p+geom_boxplot(boxplot,mapping=aes(x=variable,y=value))

#step5: Prepare to calculate AUC
#1. Prepare input data
P<-as.data.frame(Ps);colnames(P)[1]<-"ID"
P$PN<-1
N<-as.data.frame(Ns);colnames(N)[1]<-"ID"
N$PN<-0
PN<-rbind(P,N)
#import
data.pred.input<-data.preds
#
data.pred.input<-as.data.frame(data.pred.input)
data.pred.input$ID<-rownames(data.pred.input)
data.pred.input<-inner_join(data.pred.input,PN,by="ID")

AUC.input<-as.data.frame(data.pred.input$PN);colnames(AUC.input)[1]<-"PN"
AUC.input$P<-data.pred.input$P

library(ROCR)
pred<-prediction(AUC.input$P,AUC.input$PN);perf<-performance(pred,"tpr","fpr");performance(pred,'auc');auc <- performance(pred,'auc')
auc=unlist(slot(auc,"y.values"))
auc

#step6: output data
write.csv(TP,"Core-ada-predict.csv")
write.csv(c,"Core-FPR.csv")
write.csv(data.random,"Core-random.csv") 
write.csv(data.preds,"Core-data.preds.csv")
# Also output the FPR graph
```



