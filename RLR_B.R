# edgeR Repenrich2 analysis
# Author: Seung June Kim
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
packages <- c("edgeR","extrafont","ggplot2","gplots","Rtools","RColorBrewer","ggrepel","EnhancedVolcano")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
source("C:/users/jk/Documents/Lab/EnhancedVolcano.R")


font_import()

setwd("~/Lab/Spleen_RNA_seq/repenrich/counts/fraction/")
set.seed(0)

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
r.meta <-data.frame(
  row.names=colnames(x.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  libsize=c(29722003,32067122,28784483,29232676,28815362,25791057),
  replicate=c("1","2","3","1","2","3")
  
)

r.treatment<-factor(r.meta$treatment)
r.replicate<-factor(r.meta$replicate)
r.libsize <- r.meta$libsize

# paired t-test model: common effect of treatment considering paired samples
r.design<-model.matrix(~replicate+r.treatment)
rownames(design)<-colnames(x)

# histogram of rowSums
r.rowsum <- rowSums(x.counts)
plot(r.rowsum,log='y')


x<-DGEList(counts=x.counts,lib.size=colSums(x.counts), group=r.treatment)
keep.x<- filterByExpr(x, min.count=2000)
x <- x[keep.x, , keep.lib.sizes=F]
x<- calcNormFactors(x)
dim(x)

# PCA plot equivalent
plotMDS(x)

# general linear model
x.glm<-estimateDisp(x,r.design)
# how well does the model fit?
plotBCV(x.glm)

# determine differentially expressed genes
x.glm.fit<-glmFit(x.glm,r.design)
topTags(test)
x.glm.fit.ql<-glmQLFit(x.glm,r.design)
# testing the last coefficient in the linear model (GSK vs DMSO)

# LR test
#####
 x.glm.lrt<-glmLRT(x.glm.fit)
 topTags(x.glm.lrt)
 summary(decideTests(x.glm.lrt))
 x.dge.lrt <- topTags(x.glm.lrt,n=57, sort.by="PValue")
 x.dge.all<- topTags(x.glm.lrt,n="Inf", sort.by="logFC")
 x.dge.all.sort <- x.dge.all[order(x.dge.all$table$logFC),]
#####



tiff("RLR_B_MD.tiff",units='in',width=5,height=5,res=1200)
plotMD(x.glm.lrt, xlim=c(-5,15),main="GSK343 treated RLR B cells")
dev.off()



# RowSideColor vector

colfunc<-colorRampPalette(c("blue","black","red"))

# sort DGEs (by P value) now by log FC
x.dge.lrt.fcsort <- x.dge.lrt[order(x.dge.lrt$table$logFC),]

# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-RD1$V2[match(rownames(x.dge.lrt.fcsort),rownames(x.counts))]
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
types<-gsub("Other","grey",types)
types<-gsub("Unknown","grey",types)
types<-gsub("RC","black",types)

# heatmap
library(viridis)
colfunc<-colorRampPalette(c("blue","black","red"))
#svg("RLR_heatmap.svg", family='sans', width=6, height=7)
pdf(file="RLR_heatmap.pdf")
par(font=2, font.axis=2, font.lab=2,cex.lab=1.5, cex.axis=1.5,
    lty=1)
heatmap.2(cpm(x)[rownames(x.dge.lrt.fcsort),],
          col=colfunc2, scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,4),
          lwid=c(1.75,0.4, 3.5),margins=c(5,12),cexRow=1.2)
dev.off()
svg("RLR_allrepeats.svg", family='sans', width=11, height=7)

heatmap.2(cpm(x)[rownames(x.dge.all.sort),],
          col=colfunc2, scale="row", Rowv=NA, Colv=NA,
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,8), lwid=c(3,0.3, 5),margins=c(5,20),cexRow=1)
dev.off()


library(EnhancedVolcano)
r.res <- topTags(x.glm.lrt, sort.by="logFC", n="Inf")
svg("RLR_volcano.svg", family='sans', width=7, height=5)
EnhancedVolcano(data.frame(r.res),
                lab=rownames(r.res),
                selectLab=c("Lx3B","7SLRNA","FordPrefect",
                          "SSU-rRNA_Hsa",
                            "MMVL30-int",
                          "MMSAT4","L1Md_A",
                            "MMETn-int","MMERGLN-int"),
                #selectLab=c("Lx3B"),
                x='logFC',
                y='FDR',
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                title=NULL,
                subtitle=NULL,
                labSize=9,
                xlim=c(-0.4,0.9),
                ylim=c(0,20),
                axisLabSize=24,
                labFace='bold',
                drawConnectors=T,
                widthConnectors=0.5,
                boxedLabels=F,
                col=c('grey','grey','grey','red'),
                colAlpha=1,
                shadeAlpha=1,
                xlab =" ",
                ylab= " ",
                caption=NULL,
                max.overlaps=10)

dev.off()




#####
tiff(filename="RLR_heatmap_lrt.tiff",units='in',width=7,height=12,res=600)
heatmap.2(cpm(x)[rownames(x.dge.lrt.fcsort),],
          col=colfunc2, scale="row", Rowv=NA, Colv=NA, 
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lwid=c(3,0.3, 5),margins=c(5,20),cexRow=0.5)
dev.off()
#####

