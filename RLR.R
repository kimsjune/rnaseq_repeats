if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Rtools")
install.packages("edgeR")
library("edgeR")
setwd("~/Lab/Spleen_RNA_seq/repenrich/counts/fraction/")
SD1<-read.delim("RLR_B_D1_fraction_counts.txt",header=F)
SD2<-read.delim("RLR_B_D2_fraction_counts.txt",header=F)
SD3<-read.delim("RLR_B_D3_fraction_counts.txt",header=F)
SG1<-read.delim("RLR_B_G1_fraction_counts.txt",header=F)
SG2<-read.delim("RLR_B_G2_fraction_counts.txt",header=F)
SG3<-read.delim("RLR_B_G3_fraction_counts.txt",header=F)

counts<- data.frame(
  row.names=SD1[,1],
  SD1=SD1[,4],SD2=SD2[,4],SD3=SD3[,4],SG1=SG1[,4],SG2=SG2[,4],SG3=SG3[,4]
)
meta <-data.frame(
  row.names=colnames(counts),
  group=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  libsize=c(26713498,27068674,21813303,27818411,27275411,28664320),
  replicate=c("1","2","3","1","2","3")
)
libsize=meta$libsize
group<-factor(meta$group)
replicate<-factor(meta$replicate)
design<-model.matrix(~replicate+group)
colnames(design)<-levels(meta$group,meta$replicate)



y<-DGEList(counts=counts,lib.size=libsize, group=group)
y<- calcNormFactors(y)
y$samples
plotMDS(y)
y<-estimateGLMCommonDisp(y,design)
y<-estimateGLMTrendedDisp(y,design)
y<-estimateGLMTagwiseDisp(y,design)

et<-exactTest(y,pair=c("DMSO","GSK"))
topTags(et)
plotBCV(y)
logcpm<-cpm(y,log=TRUE,lib.size=libsize)
logcpm<-as.data.frame(logcpm)
colnames(logcpm)<-factor(meta$group)
cpm<-cpm(y,lib.size=libsize)
cpm<-as.data.frame(cpm     )
colnames(cpm)<-factor(meta$group)
write.table(cpm,file="cpm.txt",quote=F,sep="\t")

o <- order(et$table$PValue)
write.table(cpm(y)[o[1:300],],file="test.txt",sep="\t")
yfit<-glmQLFit(y,design)
qlf<-glmQLFTest(yfit)
topTags(qlf)
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

my.contrasts <-makeContrasts(
  GSK_DMSO= GSK - DMSO,
  levels=design
)
