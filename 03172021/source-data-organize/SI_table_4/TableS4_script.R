# Comparisons of KEGG Pathway Module counts between the three conditions ("vehicle","uninfected","benznidazole")
# Please ensure that the input files are in the same folder as the script or provide the correct path to the files.
# Please ensure that the following libraries are installed in your R environment
    # DESeq2
    # tximport
    # apeglm
# File description: TableS4_data.txt; Counts data aggregated by KEGG Pathway Modules. 
# File description: Metadata.txt; Sample label and Treatment condition.  

#Load required R packages 
library("tximport")
library("DESeq2")
library("apeglm")

#load input data: KEGG Module counts, Sample metadata
countsFile<-as.matrix(read.csv("TableS4_data.txt",sep="\t",row.names=1,header=TRUE))
metaData<-read.csv("Metadata.txt",header=TRUE,row.names=1,sep="\t")
metaData$Treatment<-factor(metaData$Treatment)

#Create DESeq Object from imported data and filter modules with low counts (total counts across all samples < 10)
dds<-DESeqDataSetFromMatrix(countData = countsFile, colData = metaData, design = ~ Treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds <- DESeq(dds, test="LRT", reduced = ~1)
resTable<-results(dds)
dds$Treatment<-relevel(dds$Treatment,ref="uninfected")
dds_pair<-DESeq(dds, test="Wald")

res_Group1<-lfcShrink(dds_pair, coef="Treatment_benznidazole_vs_uninfected",type="apeglm")
res_Group2<-lfcShrink(dds_pair, coef="Treatment_vehicle_vs_uninfected",type="apeglm")
dds$Treatment<-relevel(dds$Treatment,ref="vehicle")
dds_pair<-DESeq(dds, test="Wald")
res_Group3<-lfcShrink(dds_pair, coef="Treatment_benznidazole_vs_vehicle",type="apeglm")
Pairwise_out<-merge(as.data.frame(res_Group1),as.data.frame(res_Group2),by="row.names",all=TRUE)
Pairwise_out <- data.frame(Pairwise_out, row.names = 1)
temp<-merge(Pairwise_out,as.data.frame(res_Group3),by="row.names",all=TRUE)
temp <- data.frame(temp, row.names = 1)
Pairwise_out<-temp
temp<-Pairwise_out[,c(1,2,4,5,7,9,10,12,14,15)]
names(temp)<-c("BaseMean","LFC_BvsU","Pval_BvsU","Padj_BvsU","LFC_VvsU","Pval_VvsU","Padj_VvsU","LFC_BvsV","Pval_BvsV","Padj_BvsV")
Pairwise_out<-merge(as.data.frame(resTable),temp,by="row.names",all=TRUE)
Pairwise_out<-data.frame(Pairwise_out,row.names=1)
Pairwise_out<-Pairwise_out[-c(1:4)]
names(Pairwise_out)[1]<-"Pval_LRT"
names(Pairwise_out)[2]<-"Padj_LRT"

write.table(Pairwise_out,file="TableS4_results.txt",sep="\t",quote=FALSE)
