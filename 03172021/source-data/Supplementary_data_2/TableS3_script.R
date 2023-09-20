# Comparisons of Gene Expression Counts between the three conditions ("vehicle","uninfected","benznidazole")
# Please ensure that the input files are in the same folder as the script or provide the correct path to the files.
# Please ensure that the following libraries are installed in your R environment
    # DESeq2
    # tximport
    # apeglm

# File description: DESeqObject.Rdata; DESeq summarized experiment containing gene counts, correction factors, and sample metadata. 
# File output: S3_table.txt; Import this file in Excel (tab-separated) and process with filtering criteria (LRT P <= 0.05, LRT FC >= 2 or <= -2)
# Use the ENSEMBL Gene IDs from the filtered results as input for GO analysis (https://david.ncifcrf.gov/tools.jsp)

load("DESeqObject.Rdata")
ddsTxi <- DESeq(ddsTxi, test="LRT", reduced = ~1)
resTable<-results(ddsTxi)

#Pairwise comparisons
ddsTxi$Treatment<-relevel(ddsTxi$Treatment,ref="uninfected")
dds_pair<-DESeq(ddsTxi, test="Wald")
res_Group1<-lfcShrink(dds_pair, coef="Treatment_benznidazole_vs_uninfected",type="apeglm")
res_Group2<-lfcShrink(dds_pair, coef="Treatment_vehicle_vs_uninfected",type="apeglm")

ddsTxi$Treatment<-relevel(ddsTxi$Treatment,ref="vehicle")
dds_pair<-DESeq(ddsTxi, test="Wald")
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
names(Pairwise_out)[1] <- "baseMean"
Pairwise_out<-Pairwise_out[-c(1,3:4)]
names(Pairwise_out)[1]<-"Log2FC"
names(Pairwise_out)[2]<-"Pval_LRT"
names(Pairwise_out)[3]<-"Padj_LRT"

write.table(Pairwise_out,file="S3_table.txt",sep="\t",quote=FALSE)

