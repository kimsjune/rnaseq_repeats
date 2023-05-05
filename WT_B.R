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
WD1<-read.delim("WT_B_D1_fraction_counts.txt",header=F)
WD2<-read.delim("WT_B_D2_fraction_counts.txt",header=F)
WD3<-read.delim("WT_B_D3_fraction_counts.txt",header=F)
WG1<-read.delim("WT_B_G1_fraction_counts.txt",header=F)
WG2<-read.delim("WT_B_G2_fraction_counts.txt",header=F)
WG3<-read.delim("WT_B_G3_fraction_counts.txt",header=F)

# cbind all counts
counts<- data.frame(
  row.names=WD1[,1],
  WD1=WD1[,4],WD2=WD2[,4],WD3=WD3[,4],WG1=WG1[,4],WG2=WG2[,4],WG3=WG3[,4]
)
# histogram of rowSums
w.rowsum <- rowSums(counts)
plot(w.rowsum,log='y')
# information about columns
meta <-data.frame(
  row.names=colnames(counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK"),
  libsize=c(28386658,29155288,26340151,26579955,25341353,30930295),
  replicate=c("1","2","3","1","2","3")
  
)
libsize=meta$libsize
treatment<-factor(meta$treatment)
replicate<-factor(meta$replicate)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~replicate+treatment)
rownames(design)<-colnames(y)


y<-DGEList(counts=counts,lib.size=libsize, group=treatment)
keep.y<- filterByExpr(y, min.count=2000)
y <- y[keep.y, , keep.lib.sizes=F]
y<- calcNormFactors(y)
dim(y)
y$samples
# PCA plot equivalent
plotMDS(y)
# general linear model
y.glm<-estimateDisp(y,design)
# how well does the model fit?
plotBCV(y.glm)
# determine differentially expressed genes
y.glm.fit<-glmFit(y.glm,design)

# testing the last coefficient in the linear model (GSK vs DMSO)

#####
 y.glm.lrt<-glmLRT(y.glm.fit)
 topTags(y.glm.lrt)
 summary(decideTests(y.glm.lrt))
 y.dge.lrt<- topTags(y.glm.lrt,n=27, sort.by="PValue")
 y.dge.all<- topTags(y.glm.lrt,n="Inf", sort.by="logFC")
 y.dge.all.sort <- y.dge.all[order(y.dge.all$table$logFC),]
#####
y.glm.fit.ql<-glmQLFit(y.glm,design)
y.glm.qlt<-glmQLFTest(y.glm.fit.ql)
topTags(y.glm.qlt)
summary(decideTests(y.glm.qlt))
topTags(y.glm.qlt)

# subset DGEs
y.dge.qlt<- topTags(y.glm.qlt,n=80, sort.by="PValue")


tiff("WT_B_MD.tiff",units='in',width=5,height=5,res=1200)
plotMD(y.glm.lrt, xlim=c(-2,20),main="GSK343 treated WT B cells")
dev.off()


# RowSideColor vector

colfunc<-colorRampPalette(c("blue","black","red"))


colfunc(8)
# sort DGEs (by P value) now by log FC
y.dge.fcsort.lrt <- y.dge.lrt[order(y.dge.lrt$table$logFC),]
# side column of repeat type
# match rownames of fcsorted dge table with that of the unsorted count table
# outputs a vector of positions in the count table where those rownames are found
# search for the names of repeat types at those positions from the repeat type column from WB1.txt
types<-WD1$V2[match(rownames(y.dge.fcsort.lrt),rownames(counts))]
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

#svg("WT_heatmap.svg", family='sans', width=6, height=7)
pdf(file="WT_heatmap.pdf")
par(font=2, font.axis=2, font.lab=2,cex.lab=1.5, cex.axis=1.5,
    lty=1)
heatmap.2(cpm(y)[rownames(y.dge.fcsort.lrt),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA,
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
          lmat=rbind(c(5,0,4),c(3,1,2)),
          lhei=c(1,4),
          lwid=c(1.75,0.4, 3.5),margins=c(5,12),cexRow=1.2)
dev.off()
write.table(file="rownames.txt",t(rownames(y.dge.fcsort.qlt)), quote=F, sep=" ", row.names = F, col.names = F)

#EnhancedVolcano

res <- topTags(y.glm.lrt, sort.by="logFC", n="Inf")
svg("WT_volcano.svg", family='sans', width=7, height=5)
EnhancedVolcano(data.frame(res),
                lab=rownames(res),
                selectLab=c("7SLRNA","FordPrefect",
                            "7SK","SSU-rRNA_Hsa",
                            "MMVL30-int","Lx3B"
                            ,"RLR4_MM-int",
                            "MMSAT4","MMERGLN-int","L1Md_A"),
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
                xlim=c(-0.5,2),
                ylim=c(0,6),
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
