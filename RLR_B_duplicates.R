if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Rtools")
install.packages("edgeR")
install.packages("extrafont")
library("extrafont")
font_import()
library("edgeR")
setwd("~/Lab/Spleen_RNA_seq/repenrich/counts/fraction/")
# read count table
RD1<-read.delim("RLR_B_D1_fraction_counts.txt",header=F)
RD2<-read.delim("RLR_B_D2_fraction_counts.txt",header=F)
RD3<-read.delim("RLR_B_D3_fraction_counts.txt",header=F)
RG1<-read.delim("RLR_B_G1_fraction_counts.txt",header=F)
RG2<-read.delim("RLR_B_G2_fraction_counts.txt",header=F)
RG3<-read.delim("RLR_B_G3_fraction_counts.txt",header=F)

# cbind all counts
x.counts<- data.frame(
  row.names=RD1[,1],
  RD1=RD1[,4],RD2=RD2[,4],RD3=RD3[,4],RG1=RG1[,4],RG2=RG2[,4],RG3=RG3[,4]
)
# information about columns
meta <-data.frame(
  row.names=colnames(x.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  libsize=c(27780203,30044784,26603562,27200176,26838724,23460174),
  replicate=c("1","2","3","1","2","3")
  
)
libsize=meta$libsize
treatment<-factor(meta$treatment)
replicate<-factor(meta$replicate)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~replicate+treatment)
rownames(design)<-colnames(x)
#colnames(design)<-levels(meta$group,meta$replicate)



x<-DGEList(counts=x.counts,lib.size=libsize, group=treatment)
#keep.x<- filterByExpr(x, min.count.total=10000)
#x <- x[keep.x, , keep.lib.sizes=F]
x<- calcNormFactors(x, method="TMMwsp")
x$samples
# PCA plot equivalent
plotMDS(x)
# general linear model
x.glm<-estimateDisp(x,design)
# how well does the model fit?
plotBCV(x.glm)
# determine differentially expressed genes
x.glm.fit<-glmFit(x.glm,design, robust=F)
x.glm.fit.ql<-glmQLFit(x.glm,design, robust=F)
# testing the last coefficient in the linear model (GSK vs DMSO)
x.glm.lrt<-glmLRT(x.glm.fit)
topTags(x.glm.lrt)
summary(decideTests(x.glm.lrt))
x.glm.qlt<-glmQLFTest(x.glm.fit.ql)
topTags(x.glm.qlt)
summary(decideTests(x.glm.qlt))
# subset DGEs
x.dge.qlt<- topTags(x.glm.qlt,n=55, sort.by="PValue")
x.dge.lrt <- topTags(x.glm.lrt,n=86, sort.by="PValue")


tiff("RLR_B_MD.tiff",units='in',width=5,height=5,res=1200)
plotMD(x.glm.qlt, xlim=c(-5,15),main="GSK343 treated RLR B cells")
dev.off()
library(Glimma)


######################################################################################
######################################################################################
# heatmap
install.packages("gplots")
install.packages("RColorBrewer")
library(gplots)
library(RColorBrewer)
colfunc<-colorRampPalette(c("darkblue","white","red"))
colfunc2<-colorRampPalette(c("blue","black","red"))

colfunc(8)
# sort DGEs (by P value) now by log FC
x.dge.qlt.fcsort <- x.dge.qlt[order(x.dge.qlt$table$logFC),]
x.dge.lrt.fcsort <- x.dge.lrt[order(x.dge.lrt$table$logFC),]

# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-RD1$V2[match(rownames(x.dge.qlt.fcsort),rownames(x.counts))]
types.lrt<-RD1$V2[match(rownames(x.dge.lrt.fcsort),rownames(x.counts))]

#install.packages("stringr")
library(stringr)

types<-str_replace(types,"tRNA","pink2")
types<-str_replace(types,"LINE","green3")
types<-str_replace(types,"LINE\\?","green3")

types<-str_replace(types,"DNA","black")
types<-str_replace(types,"DNA?","black")
types<-str_replace(types,"RC\\?","black")
types<-str_replace(types,"RC","black")



types<-str_replace(types,"LTR","slateblue3")
types<-str_replace(types,"LTR\\?","slateblue3")

types<-str_replace(types,"SINE","green3")
types<-str_replace(types,"Satellite","yellow3")
types<-str_replace(types,"snRNA","blue")
types<-str_replace(types,"srpRNA","blue")
types<-str_replace(types,"rRNA","blue")
types<-str_replace(types,"RNA","blue")
types<-str_replace(types,"Other","blue")
types<-str_replace(types,"Unknown","blue")












# sort DGEs (by P value) now by log FC

tiff(filename="RLR_heatmap_qlt.tiff",units='in',width=7,height=12,res=600)
heatmap.2(cpm(x)[rownames(x.dge.qlt.fcsort),],
          col=colfunc2, scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
           lwid=c(3,0.3, 5),margins=c(5,20),cexRow=1)
dev.off()

tiff(filename="RLR_heatmap_lrt.tiff",units='in',width=7,height=12,res=600)
heatmap.2(cpm(x)[rownames(x.dge.lrt.fcsort),],
          col=colfunc, scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lwid=c(3,0.3, 5),margins=c(5,20),cexRow=1)
dev.off()


# up and downregulated repeat list
deg.rlr<- topTags(x.glm.qlt,n=Inf, p=0.05)$table
up.rlr<- row.names(deg.rlr[deg.rlr$logFC>0,])
down.rlr <- row.names(deg.rlr[deg.rlr$logFC<0,])




#####################################################################


o<-order(y.glm.lrt$table$PValue, y.glm.lrt$table$logFC)

degenes<-cpm(y)[o[1:82],]
write.table(degenes,file="WT_degenes.txt",quote=F,sep="\t")


logcpm<-cpm(y,log=TRUE,lib.size=libsize)
logcpm<-as.data.frame(logcpm)
colnames(logcpm)<-factor(meta$group)
cpm<-cpm(y,lib.size=libsize)
cpm<-as.data.frame(cpm     )
colnames(cpm)<-factor(meta$group)
write.table(cpm,file="cpm.txt",quote=F,sep="\t")


yfit<-glmQLFit(y,design)
qlf<-glmQLFTest(yfit)
topTags(qlf)
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

et<-exactTest(x,pair=c("DMSO","GSK"))
