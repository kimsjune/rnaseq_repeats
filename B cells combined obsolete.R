# cbind all counts
all.counts<- data.frame(
  
  WD1=WD1[,4],WD2=WD2[,4],WD3=WD3[,4],WG1=WG1[,4],WG2=WG2[,4],WG3=WG3[,4],
  RD1=RD1[,4],RD2=RD2[,4],RD3=RD3[,4],RG1=RG1[,4],RG2=RG2[,4],RG3=RG3[,4]
)
# information about columns
all.meta <-data.frame(
  row.names=colnames(all.counts),
  treatment=c("DMSO","DMSO","DMSO","GSK","GSK","GSK","DMSO","DMSO","DMSO","GSK","GSK","GSK"
             ),
 # libsize=c(26516586,27145908,24603751,24236315,23238498,27951010,
 #           27780203,30044784,26603562,27200176,26838724,23460174),
  replicate=c("1","2","3","1","2","3","4","5","6","4","5","6"),
  genotype=c(rep("WT",6),rep("RLR",6))
  
)

#group<- factor(paste(all.meta$treatment,all.meta$genotype,sep="."))
#cbind(all.meta,group=group)
#libsize=z.meta$libsize
treatment<-factor(all.meta$treatment)
replicate<-factor(all.meta$replicate)
genotype<-factor(all.meta$genotype)
# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~replicate)
WT.GSK <- genotype == "WT" & treatment== "GSK"
RLR.GSK <- genotype == "RLR" & treatment == "GSK"
design <- cbind(design, WT.GSK,RLR.GSK)



#colnames(design)<-levels(group)



a<-DGEList(counts=all.counts, group=group)
#keep<- filterByExpr(a)
#a <- a[keep, , keep.lib.sizes=F]
a<- calcNormFactors(a, method="TMMwsp")
a$samples
# PCA plot equivalent
tiff("PCA.tiff")
plotMDS(a)
dev.off()
# general linear model
a.glm<-estimateDisp(a,design)
# how well does the model fit
plotBCV(a.glm)
# determine differentially expressed genes
a.glm.fit<-glmFit(a.glm,design,robust=F)
a.glm.fit.ql<-glmQLFit(a.glm,design,robust=F)
test <- glmQLFTest(a.glm.fit.ql, coef="WT.GSK")
topTags(test)
summary(decideTests(test))
test2 <- glmQLFTest(a.glm.fit, coef="RLR.GSK")
topTags(test2)
summary(decideTests(test2))

my.contrasts<- makeContrasts(
  WT.GvD= GSK.WT-DMSO.WT,
  RLR.GvD= GSK.RLR-DMSO.RLR,
  GSK.WTvRLR=GSK.WT-GSK.RLR,
  DMSO.WTvRLR=DMSO.WT-DMSO.RLR,
  levels=design
)

# DGE from a given comparison

a.gsk.wt.rlr.lrt <- glmLRT(a.glm.fit, contrast=my.contrasts[,"GSK.WTvRLR"])
a.dmso.wt.rlr.lrt <- glmLRT(a.glm.fit, contrast=my.contrasts[,"DMSO.WTvRLR"])

a.wt.gskvdmso.lrt <- glmLRT(a.glm.fit,contrast=my.contrasts[,"WT.GvD"])
a.rlr.gskvdmso.lrt <- glmLRT(a.glm.fit,contrast=my.contrasts[,"RLR.GvD"])

a.wt.gskvdmso.qlt<- glmQLFTest(a.glm.fit.ql,contrast=my.contrasts[,"WT.GvD"])
a.rlr.gskvdmso.qlt <- glmQLFTest(a.glm.fit.ql,contrast=my.contrasts[,"RLR.GvD"])

summary(decideTests(a.wt.gskvdmso.lrt))
summary(decideTests(a.rlr.gskvdmso.lrt))
summary(decideTests(a.wt.gskvdmso.qlt))
summary(decideTests(a.rlr.gskvdmso.qlt))
a.wt.gskvdmso.dge<- topTags(a.wt.gskvdmso.lrt,n=55, sort.by="PValue")
library(gplots)
library(RColorBrewer)
# coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(8))
colfunc<-colorRampPalette(c("darkblue","red"))
colfunc2<-colorRampPalette(c("black","darkblue","red"))
colfunc(10)
# sort DGEs (by P value) now by log FC
a.wt.gskvdmso.dge.fcsort <- a.wt.gskvdmso.dge[order(a.wt.gskvdmso.dge$table$logFC),]

tiff(filename="test.tiff",units='in',width=7,height=6,res=600)
heatmap.2(cpm(a[,c(1,2,3,4,5,6)])[rownames(a.wt.gskvdmso.dge.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c( 4, 10),margins=c(5,10),cexRow=1)

dev.off()

summary(decideTests(a.rlr.gskvdmso.lrt))
a.rlr.gskvdmso.dge<- topTags(a.rlr.gskvdmso.lrt,n=55, sort.by="PValue")
library(gplots)
library(RColorBrewer)
# coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(8))
colfunc<-colorRampPalette(c("darkblue","red"))
colfunc2<-colorRampPalette(c("black","darkblue","red"))
colfunc(10)
# sort DGEs (by P value) now by log FC
a.rlr.gskvdmso.dge.fcsort <- a.rlr.gskvdmso.dge[order(a.rlr.gskvdmso.dge$table$logFC),]

tiff(filename="test.tiff",units='in',width=7,height=6,res=600)
heatmap.2(cpm(a[,c(7,8,9,10,11,12)])[rownames(a.rlr.gskvdmso.dge.fcsort),],
          col=colfunc(10), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c( 4, 10),margins=c(5,10),cexRow=1)

dev.off()