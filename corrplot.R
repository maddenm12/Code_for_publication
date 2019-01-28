library(corrplot)
setwd("C:/Users/mmoor/OneDrive/Desktop/GDC_Data/luad_mirnaseq/Blood_Results")
data=read.table("GSE61741table.txt",header=T,row.names=1)
sig_miRs=read.table("miRs_both_15and20_blood_rf.txt",header=T,stringsAsFactors=F)
sig_miRs=sig_miRs[,1]
l=nchar(sig_miRs[1])
substr(sig_miRs[1],l,l)
substr(sig_miRs[1],1,l-1)
for(i in 1:length(sig_miRs))
{
  l=nchar(sig_miRs[i])
  if(substr(sig_miRs[i],l,l)==".")
  {
    temp=paste(substr(sig_miRs[i],1,l-1),"*",sep="")
    sig_miRs[i]=temp
  }
}

sig_miRs=gsub("\\.","-",sig_miRs)
indicies=vector()
for(i in 1:length(sig_miRs))
{
  temp=which(rownames(data)==sig_miRs[i])
  indicies=c(indicies,temp)
}
indicies
data=data[indicies,]
new_miR_names=read.table("15_to_20_10outof10_miRs.txt",header=T,stringsAsFactors=F)[,1]
rownames(data)=new_miR_names
#transpose
data2=t(data)
cor(data2,method="pearson")
data2=as.matrix(data2)
correlation_data=cor(data2,method="pearson")
corrplot(correlation_data)
