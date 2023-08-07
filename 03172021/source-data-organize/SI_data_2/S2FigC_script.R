# This script generates a heatmap of differentially expressed genes across the three conditions ("vehicle","uninfected","benznidazole")
# Please ensure that the input files are in the same folder as the script or provide the correct path to the files.
# Please ensure that the library "pheatmap" is installed in your R environment
# File description: S2FigC_data.txt; Counts data for differentially expressed genes. Genes are annotated with Gene Symbol, Chromosome, MGI ID, and GO Biological Process
# File description: Metadata.txt; Sample label and Treatment condition.  

library("pheatmap")
dataIn<-read.table("S2FigC_data.txt",header=TRUE,row.names=1,sep="\t")
metaIn<-read.table("Metadata.txt",header=TRUE,row.names=1,sep="\t")
annot_row<-data.frame(GO_Annot=factor(dataIn$GO_Annot))
annot_col<-data.frame(Treatment=factor(metaIn$Treatment))
row.names(annot_row)<-row.names(dataIn)
row.names(annot_col)<-row.names(metaIn)
pdf("S2FigC_plot.pdf")
pheatmap(log(dataIn[,1:15]+1),scale="row",annotation_row = annot_row, annotation_col=annot_col, show_rownames =F, show_colnames = F)
dev.off()