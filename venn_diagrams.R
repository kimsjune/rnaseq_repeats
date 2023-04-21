# Follows from WT_B and RLR_B.R
# venn diagram
# up and downregulated repeat list
deg.rlr<- topTags(x.glm.lrt,n=Inf, p=0.05)$table
up.rlr<- row.names(deg.rlr[deg.rlr$logFC>0,])
down.rlr <- row.names(deg.rlr[deg.rlr$logFC<0,])

# up and downregulated repeat list
deg.wt<- topTags(y.glm.lrt ,n=Inf, p=0.05)$table
up.wt<- row.names(deg.wt[deg.wt$logFC>0,])
down.wt <- row.names(deg.wt[deg.wt$logFC<0,])

install.packages("VennDiagram")

library(VennDiagram)

venn.up <-venn.diagram(
  x=list(up.wt,up.rlr),
  category.names=c("",""),
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
  inverted=T,
  
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily='sans'
)
ggsave(venn.up, file="venn_up.svg",device="svg",height=4,width=4) 
dev.off()
venn.dn <- venn.diagram(
  x=list(down.wt,down.rlr),
  category.names=c("",""),
  resolution=600,
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
  inverted=F,
  
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily='sans'
)
ggsave(venn.dn, file="venn_dn.svg",device="svg",height=4,width=4) 
dev.off()
  