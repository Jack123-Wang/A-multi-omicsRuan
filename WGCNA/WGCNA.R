#step1:Import datExpr和datTraits
library(dplyr)
##1.1 Organize data0 into an expression matrix, the sample name and gene name are both colnames and rownames
data0<-read.csv("expreSet.csv")
rownames(data0)<-data0$X
data0<-distinct(data0,X,.keep_all = T)
data0<-data0[,-1]

##1.2 After data0 is inverted, import datExpr
data01 = t(data0)
datExpr<-data01
#write.csv(datExpr,"datExpr.csv")
##1.3 Import datTraits
datTraits<-read.csv("datTraits.csv",header = T)

#step2:Determining the optimal beta value
library(WGCNA)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#step3:One-step construction of co-expression matrix
#After setting the table property of datExpr.csv to "number", read it again
datExpr<-read.csv("datExpr.csv")
rownames(datExpr)<-datExpr$X
datExpr<-datExpr[,-1]

net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors) 

#step4: Module Visualization
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.

#step5: Relationship between modules and traits
#the number of samples and the number of genes
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
names(datTraits)[2]<-"subtype"
design=model.matrix(~0+ datTraits$subtype)
colnames(design)=levels(datTraits$subtype)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); 
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#step6: Export the gene set corresponding to each module
library(stringr)
a<-rownames(moduleTraitCor) # extract module color
b<-str_replace(a,"ME","") # remove prefix "ME"

for (i in 1:length(b)){
  module = b[i] # Select module
  # Select module probes
  probes = colnames(datExpr) 
  inModule = (moduleColors==module);
  c = probes[inModule];  #c is the gene for each module
  write.csv(c,paste(b[i],".csv",sep=""))
}

