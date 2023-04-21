# edgeR Repenrich2 analysis
# Author: Seung June Kim
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Rtools")
install.packages("edgeR")

library(edgeR)
setwd("~/Lab/Spleen_RNA_seq/repenrich/counts/fraction/")
# read count table
SD1<-read.delim("SPL_D1_fraction_counts.txt",header=F)
SD2<-read.delim("SPL_D2_fraction_counts.txt",header=F)
SD3<-read.delim("SPL_D3_fraction_counts.txt",header=F)
SG1<-read.delim("SPL_G1_fraction_counts.txt",header=F)
SG2<-read.delim("SPL_G2_fraction_counts.txt",header=F)
SG3<-read.delim("SPL_G3_fraction_counts.txt",header=F)

# cbind all counts
z.counts<- data.frame(
  row.names=SD1[,1],
  SD1=SD1[,4],SD2=SD2[,4],SD3=SD3[,4],SG1=SG1[,4],SG2=SG2[,4],SG3=SG3[,4]
)
# information about columns
z.meta <-data.frame(
  row.names=colnames(z.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  libsize=c(28758501,28990003,23278579,29988755,29379446,30543215),
  replicate=c("1","2","3","1","2","3")
  
)
z.libsize=z.meta$libsize
z.treatment<-factor(z.meta$treatment)
z.replicate<-factor(z.meta$replicate)

# paired t-test model: common effect of treatment considering paired samples
z.design<-model.matrix(~z.replicate+z.treatment)
#rownames(design)<-colnames(z.counts)


z<-DGEList(counts=z.counts,lib.size=z.libsize, group=treatment)
keep.z<- filterByExpr(z, min.count=2000)
z <- z[keep.z, , keep.lib.sizes=F]
z<- calcNormFactors(z)

# PCA plot equivalent
plotMDS(z)

# general linear model
z.glm<-estimateDisp(z,z.design)
# how well does the model fit?
plotBCV(z.glm)

# determine differentially expressed genes
z.glm.fit<-glmFit(z.glm,design)
# Likelihood ratio test
z.glm.lrt<-glmLRT(z.glm.fit)
topTags(z.glm.lrt)
summary(decideTests(z.glm.lrt))

## QLF Test -- obsolete
#####
z.glm.fit.ql<-glmQLFit(z.glm,design)
z.glm.qlt<-glmQLFTest(z.glm.fit.ql)
topTags(z.glm.qlt)
summary(decideTests(z.glm.qlt))
#####

# subset DGEs
z.dge.lrt<- topTags(z.glm.lrt,n=7, sort.by="PVal")


tiff("SPL_MD.tiff",units='in',width=5,height=5,res=1200)
z.md<-plotMD(z.glm.lrt, xlim=c(-2,20),main="GSK343 treated WT splenocytes")
dev.off()

# make RowSideColor vector
install.packages("gplots")
install.packages("RColorBrewer")
library(gplots)
library(RColorBrewer)
colfunc<-colorRampPalette(c("blue","black","red"))
colfunc(10)
# sort DGEs (by P value) now by log FC
z.dge.lrt.fcsort <- z.dge.lrt[order(z.dge.lrt$table$logFC),]
# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-SD1$V2[match(rownames(z.dge.lrt.fcsort),rownames(z.counts))]
types <- gsub("\\?.*","",types)
#install.packages("stringr")
types<-gsub("tRNA","pink2",types)
types<-gsub("LINE","green3",types)
types<-gsub("DNA","black",types)


types<-gsub("LTR","slateblue3",types)
types<-gsub("SINE","green3",types)
types<-gsub("Satellite","yellow3",types)
types<-gsub("snRNA","blue",types)
types<-gsub("srpRNA","blue",types)
types<-gsub("rRNA","blue",types)
types<-gsub("RNA","blue",types)
types<-gsub("Other","blue",types)
types<-gsub("Unknown","blue",types)



# heatmaps
svg("SPL_heatmap.svg", family='sans', width=6, height=7)
heatmap.2(cpm(z)[rownames(z.dge.lrt.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)), lhei=c(2,9),
          lwid=c(3,0.5, 5),cexRow=1.5,margins=c(1,10)
)
dev.off()
pdf(file="SPL_heatmap.pdf")
par(font=2, font.axis=2, font.lab=2,cex.lab=1.5, cex.axis=1.5,
    lty=1)
heatmap.2(cpm(z)[rownames(z.dge.lrt.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA,
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,4),
          lwid=c(1.75,0.4, 3.5),margins=c(5,12),cexRow=1.2)
dev.off()


#EnhancedVolcano
library(EnhancedVolcano)
z.res <- topTags(z.glm.lrt, sort.by="logFC", n="Inf")
svg("SPL_volcano.svg", family='sans', width=7, height=5)
EnhancedVolcano(data.frame(z.res),
                lab=rownames(z.res),
                selectLab=c("RLTR6_Mm","RLTR1B-int","MuRRS4-int",
                            "RLTR6-int","LSU-rRNA_Hsa"),
                x='logFC',
                y='FDR',
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                labFace='bold',
                title=NULL,
                subtitle=NULL,
                labSize=9,
                xlim=c(-2,2),
                #ylim=c(0,6),
                axisLabSize=24,
                drawConnectors=T,
                widthConnectors=0.5,
                boxedLabels=F,
                col=c('grey','grey','grey','red'),
                colAlpha=1,
                xlab =" ",
                ylab= " ",
                caption=NULL
)
dev.off()











tiff(filename="dummy_heatmap.tiff",units='in',width=7,height=6,res=1200)
heatmap.2(cpm(z)[rownames(z.dge.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = dummy,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1.5,6), lwid=c(2.5,0.3, 5),cexRow=1.4,margins=c(5,20))
dev.off()

