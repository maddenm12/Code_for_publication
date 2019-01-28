setwd("C:/Users/mmoor/OneDrive/Desktop/GDC_Data/luad_mirnaseq/Tissue_Results_RPM/")
sig_miRs_deseq=unlist(read.table("sig_miRs_deseq2.txt",header=T,stringsAsFactors=F))
a=read.table("raw_reads_luad_na_replaced.txt",header=T,row.names=1)
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
remove=which(miR_names=="Empty")
a=a[-remove,]
miR_names=miR_names[-remove]
rownames(a)=miR_names
a=a[sig_miRs_deseq,]
data2=t(a)
#cor(data2,method="pearson")
data2=as.matrix(data2)
correlation_data=cor(data2,method="pearson")
corrplot(correlation_data)

setwd("C:/Users/mmoor/OneDrive/Desktop/GDC_Data/luad_mirnaseq/Tissue_Results_RPM/")
sig_miRs_rf=unlist(read.table("most_predictive_miRs_10outof10.txt",header=T,stringsAsFactors=F))
sig_miRs_rf=gsub("\\.","-",sig_miRs_rf)
a=read.table("raw_reads_luad_na_replaced.txt",header=T,row.names=1)
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
remove=which(miR_names=="Empty")
a=a[-remove,]
miR_names=miR_names[-remove]
rownames(a)=miR_names
a=a[sig_miRs_rf,]
data2=t(a)
#cor(data2,method="pearson")
data2=as.matrix(data2)
correlation_data=cor(data2,method="pearson")
corrplot(correlation_data)

cor.test(data2[,"hsa-miR-331-3p"],data2[,"hsa-miR-144-5p"])

top3=c("hsa-miR-144-5p","hsa-miR-331-3p","hsa-miR-423-3p")
data3=data2[,top3]
correlation_data2=cor(data3,method="pearson")
corrplot(correlation_data2)

