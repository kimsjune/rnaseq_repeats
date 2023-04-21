# edgeR Repenrich2 analysis
# Author: Seung June Kim
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
# LRT --obsolete
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
install.packages("gplots")
install.packages("RColorBrewer")
library(gplots)
library(RColorBrewer)
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

# heatmap
# LRT --obsolete
#####
# tiff(filename="WT_heatmap_lrt.tiff",units='in',width=7,height=12,res=1200)
# heatmap.2(cpm(y)[rownames(y.dge.fcsort),],
#           col=colfunc2, scale="row", Rowv=NA, Colv=NA,
#           cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
#           density.info = "none",trace="none", dendrogram = "none",
#           symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
#           lmat=rbind(c(5,0,4),c(3,1,2)),
#           lhei=c(1,8), lwid=c(3,0.3, 5),margins=c(5,20),cexRow=1.3)
# dev.off()
#####
#tiff(filename="WT_heatmap_qlt.tiff",units='in',width=7,height=12,res=1200)
#heatmap.2(cpm(y)[rownames(y.dge.fcsort.lrt),],
#          col=colfunc2, scale="row", Rowv=NA, Colv=NA,
#          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
#          density.info = "none",trace="none", dendrogram = "none",
#          symkey=FALSE,symbreaks=TRUE,revC = FALSE, RowSideColors = types,
#          lmat=rbind(c(5,0,4),c(3,1,2)),
#          lhei=c(1,8), lwid=c(3,0.3, 5),margins=c(5,20),cexRow=1)
#dev.off()

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
library(EnhancedVolcano)
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



### RT-qPCR heatmaps OBSOLETE
L1ORF1a <- c(1.122177458,
             0.918290561,
             0.959531981,
             0.684175008,
             0.888232465,
             1.427592527,
             1.082998169,
             0.892381607,
             1.288763088,
             1.314370747,
             1.620535122,
             1.982208957
)
L1ORF1b <- c(1.034249304,
             0.948205158,
             1.017545538,
             0.637538249,
             0.971149498,
             1.391312253,
             1.067347953,
             0.731595899,
             1.226734139,
             1.508073614,
             1.818490755,
             1.653265162
)
L1ORF2a<-c(1.155172796,
           0.908541801,
           0.936285402,
           0.611705393,
           0.827748928,
           1.56054568,
           1.254994058,
           0.7845881,
           1.105958782,
           1.275800167,
           1.363330606,
           2.179661727
)
L1ORF2b<-c(1.104008994,
           1.026809003,
           0.869182003,
           0.613481904,
           0.932661117,
           1.453856979,
           1.118540605,
           0.636893095,
           1.163127733,
           1.179625408,
           1.369228401,
           2.136557565
)
muERVL <- c(1.878682287,
            0.372248657,
            0.749069056,
            1.590977604,
            0.923708241,
            0.485314155,
            0.870879954,
            0.509679202,
            0.723963352,
            0.627777952,
            1.086634208,
            0.839095993
)

musD <- c(1.448736936,
          0.848057427,
          0.703205638,
          0.785236067,
          0.983080034,
          1.231683898,
          0.771832069,
          0.528536326,
          0.790808125,
          1.099452209,
          1.359575189,
          1.549171723
)
IAP_G <- c(1.164048275,
           0.83660903,
           0.999342695,
           1.018443134,
           1.231587685,
           0.749969181,
           1.390771491,
           0.677247444,
           0.83572987,
           1.249986348,
           1.804741548,
           0.78900664
)
IAP_R<- c(1.444099796,
          0.87224132,
          0.683658885,
          1.347851447,
          1.022708966,
          0.629439587,
          1.030533893,
          0.543163425,
          0.78408925,
          0.966949407,
          1.465214495,
          0.696388425
)
MuLV<-c(1.089818663,
        0.967249037,
        0.942932299,
        0.892657526,
        1.00150621,
        1.105836264,
        1.03537379,
        0.809944929,
        1.068098922,
        0.988001848,
        1.041955157,
        1.173355684
)
RLTR4int1 <-c(1.30679259,
              0.869730086,
              0.823477324,
              1.045955473,
              0.731513949,
              1.222530579,
              0.435070298,
              0.682483188,
              0.641807681,
              1.188841165,
              0.772073278,
              1.32692087
)
colfunc1<- colfunc<-colorRampPalette(c("darkblue","red")(n=1000))
colors = c(seq(0.8,1.2,length=8),seq(1.21,2,length=20))
tiff(filename="RT_qpcr.tiff",units='in',width=4,height=6,res=1200)
heatmap.2(t(cbind(L1ORF1a,L1ORF1b,L1ORF2a,L1ORF2b,muERVL,musD)),breaks=colors,
          col=colfunc1, scale="none", Rowv=NA, Colv=NA,
          cexCol=1, labCol=NA,srtCol=0, adjCol=c(0.5,0),
          density.info = "none",trace="none", 
          dendrogram = "none", symkey=FALSE,
          symbreaks=FALSE,revC = FALSE,
          lmat=rbind(c(0,3,4),c(0,1,2),c(0,0,0)),
          lhei=c(1,2,2), lwid=c(0.5,3,2),cexRow=1
          )
dev.off()
