 

### **Calculate Z-Score**

#step1: keep the first three columns

```bash
ls *.narrowPeak|while read id;do awk -v OFS='\t' '{print $1,$2,$3}' $id > $(basename $id _peaks.narrowPeak).sed.peak;done 
```

#step2: Create peak random distribution file

```bash
ls *.sed.peak|while read id;do ls *.sed.peak|while read ID;do echo $id $ID;done;done> config
ls *.sed.peak|while read id;do for i in {1..2000};do bedtools shuffle -i $id -g chr.txt > $(basename $id .sed.peak).shuffle.$i;done;done
```

#step3-4: Record the length occupied by different peak files
#step3: Different

#Length of real merge peak (1 line)

```bash
cat config2 |while read id;do arr=($id); b=${arr[1]}; a=${arr[0]};bedtools intersect -a $a -b $b -wo|cut -f 7|awk '{a+=$1}END{print a}' >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig ;done
```

#Length distribution of random merge peaks (2-1001 lines); number of peak intersections of random merge peaks (1002-2001)

```bash
cat config2 |while read id;do arr=($id); b=${arr[1]}; a=${arr[0]};for i in {1..1000};do bedtools intersect -a $(basename $a .sed.peak).shuffle.$i -b $(basename $b .sed.peak).shuffle.$i -wo|cut -f 7 >$(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i;cat $(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i | awk '{a+=$1}END{print a}' >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig;done;for i in {1..1000};do cat $(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i | wc|sed 's/^ *//g'|cut -d ' ' -f 1 >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig ;rm $(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i ;done;done
```

#The number of peak intersections of the real merge peak (line 2002) (different)

```bash
cat config2 |while read id;do arr=($id); b=${arr[1]}; a=${arr[0]};bedtools intersect -a $a -b $b |wc|sed 's/^ *//g'|cut -d ' ' -f 1 >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig ;mv $(basename $a .sed.peak)_$(basename $b .sed.peak).fig $(basename $a .sed.peak)_$(basename $b .sed.peak).fig_ed ;done
```

#step4: the same
#The number of peak intersections of the real merge peak (1 line)

```bash
cat config1 |while read id;do arr=($id); b=${arr[1]}; a=${arr[0]};bedtools intersect -a $a -b $b -wo|cut -f 7|awk '{a+=$1}END{print a}' >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig ;done
```

#number of peak intersections of random merge peaks (2-1001)

```bash
cat config1 |while read id;do arr=($id); b=${arr[1]}; a=${arr[0]};for i in {1..1000};do cat $(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i | wc|sed 's/^ *//g'|cut -d ' ' -f 1 >>$(basename $a .sed.peak)_$(basename $b .sed.peak).fig ;rm $(basename $a .sed.peak)_$(basename $b .sed.peak).tmp.$i;done;done
```

#Calculate Z-score score

```R
setwd("Data_HD_2/ry_chipseq-task/Zscore")
temp<-list.files(pattern ="*.fig_ed")
options(stringsAsFactors = F)
library(data.table)
for (i in 1:1){
  assign("a",read.csv(temp[i],,sep = "",header=F,blank.lines.skip=F))
  a[is.na(a)] <- 0
  a1000<-a[c(2:1001),1]
  k2<-as.numeric((a[1,1]-mean(a1000))/sd(a1000)) #Z-score:peak distribution
  k0<-cbind(k2,temp[i])
}
for (i in 2:length(temp)){
  a<-fread(temp[i],blank.lines.skip=F)
  a<-as.data.frame(a)
  a[is.na(a)] <- 0
  a1000<-a[c(2:1001),1]
  k2<-as.numeric((a[1,1]-mean(a1000))/sd(a1000)) #Z-score:peak distribution
  k01<-cbind(k2,temp[i])
  k0<-rbind(k0,k01)
}
k1<-k0
k1<-as.data.frame(k1)
k1$k2<-as.numeric(k1$k2)

library(dplyr)
library(stringr)
str1<-str_replace(k1[,3],".fig_ed","")
str1<-str_replace(str1,".peaks.narrowPeak","")
str1
str2<-str_replace(str1,"_GSE",".GSE")
str2
k1$V4<-word(str2,1,sep=fixed('.'))
k1$V5<-word(str2,2,sep=fixed('.'))

k1<-k1[-3]
tk1<-k1
colnames(tk1)[3]<-"V5"
colnames(tk1)[4]<-"V4"
ktotal<-rbind(k1,tk1)
ktotal<-unique(ktotal)

r2<-matrix(rep(NA,171),171,171)#create empty matrix
r2<-data.frame(r2)
colnames(r2)<-unique(ktotal$V4)
rownames(r2)<-unique(ktotal$V4)

ktotal<-arrange(ktotal,V4)


for (i in 1:nrow(r2)){
  for(j in 1:nrow(r2)){
    r2[ktotal$V4[220*(i-1)+1],ktotal$V5[j]]<-filter(ktotal,V4==ktotal$V4[220*(i-1)+1] & V5==ktotal$V5[j])[2]
  }
}
write.csv(r2,"r2.csv")
col<- colorRampPalette(c("black","yellow","yellow", "yellow"))(20)
library(heatmap3)
library(gplots)

heatmap3<-heatmap3(r2,
                   margins=c(3,3),
                   scale="none",    
                   col = col, 
                   Colv = "Rowv", 
                   symm = TRUE,
                   showColDendro = T, showRowDendro = F, 
                   cexRow=0.1, 
                   cexCol=0.5
)

```

