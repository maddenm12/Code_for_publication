a=getwd()
a
setwd("C:/Users/mmoor/OneDrive/Desktop/GDC_Data/luad_mirnaseq")
dir()
a=data.matrix(read.table("luad_raw_reads_matrix.txt",header=T,row.names=1,stringsAsFactors=F))
a[is.na(a)]=0
str(a)

b=read.table("hsa_miR_accessionTOname_copy.txt",header=T,stringsAsFactors=F)
miR_names=vector()
for(i in 1:length(rownames(a)))
{
  temp=rownames(a)[i]
  temp2=which(b[,1]==temp)
  if(length(temp2)==0)
  {
    miR_name="Empty"
  } else miR_name=b[temp2,2]
  miR_names=c(miR_names,miR_name)
}
miR_names
rownames(a)=miR_names
remove=which(rownames(a)=="Empty")
a=a[-remove,]
dim(a)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----quickStart, eval=FALSE----------------------------------------------
data=a
data=as.matrix(data)
dim(data)
length(rownames(data)) 
miR_data=a
tumor_samples=c(which(grepl(".01A.",colnames(miR_data))),which(grepl(".01B.",colnames(miR_data))),which(grepl(".02A.",colnames(miR_data))),which(grepl(".06A.",colnames(miR_data))))
length(tumor_samples)
control_samples=c(which(grepl(".11A.",colnames(miR_data))),which(grepl(".11B.",colnames(miR_data))))
length(control_samples)
length(tumor_samples)+length(control_samples)
length(colnames(a))
x=vector()
x=rep("",length(colnames(miR_data)))
x[tumor_samples]="Tumor"
x[control_samples]="Normal"
conditions2=as.data.frame(x)
rownames(conditions2)=colnames(a)
colnames(conditions2)=c("condition")
colnames(data)=gsub("\\.","-",colnames(data))
rownames(conditions2)=gsub("\\.","-",rownames(conditions2))
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = conditions2,
                              design= ~ condition)
dds <- DESeq(dds,parallel=TRUE)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_Tumor_vs_Normal")
res2=as.matrix(res)
sig_miRs=rownames(res2)[which(res2[,6]<10^-30)]
length(which(res2[,6]<0.1))
length(res2[,6])
sig_miRs
head(res2)
write.table(sig_miRs,file="sig_miRs_deseq2.txt",sep="\t") 





