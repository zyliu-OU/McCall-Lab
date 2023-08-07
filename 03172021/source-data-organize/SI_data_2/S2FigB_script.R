# Script to generate PCA plot of samples
# Please ensure that the input files are in the same folder as the script or provide the correct path to the files.
# Please ensure that the following libraries are installed in your R environment
    # DESeq2
    # ggplot2

# File description: DESeqObject.Rdata; DESeq summarized experiment containing gene counts, correction factors, and sample metadata. 
# File output: S2FigB_plot.pdf; PCA plot.

library("DESeq2")
library("ggplot2")

load("DESeqObject.Rdata")
rld<-rlog(ddsTxi,blind=TRUE)
PCAData<-as.data.frame(plotPCA(rld,intgroup="Treatment",returnData=TRUE))
pdf("S2FigB_plot.pdf")
ggplot(data=PCAData,aes(x=PC1,y=PC2, color=Treatment))+geom_point()
dev.off()
