x=read.table("random_forest_preprocessed_data_lungcancertissue_vs_prostatecancertissue.txt",header=T)
dim(x)
#hsa.let.7i.5p hsa.miR.181a.5p
which(colnames(x)=="hsa.let.7i.5p")
which(colnames(x)=="hsa.miR.181a.5p")
miR_indicies=c(19,252)
which(x[,1]=="PRAD")
prad_indicies=1:499
which(x[,1]=="LUAD")
luad_indicies=500:1020
foldchanges=vector()
foldchanges=c(foldchanges,mean(x[prad_indicies,19])/mean(x[luad_indicies,19]))
foldchanges=c(foldchanges,mean(x[prad_indicies,252])/mean(x[luad_indicies,252]))
foldchanges
table=data.frame()
miRs=c("hsa.let.7i.5p/hsa.let.7i", "hsa.miR.181a.5p/hsa.miR.181a")
foldchanges_tissue=foldchanges
table=cbind(miRs,foldchanges_tissue)
table
x=read.table("random_forest_preprocessed_data_pradVSluad_blood.txt",header=T)
which(colnames(x)=="hsa.let.7i")
which(colnames(x)=="hsa.miR.181a")
miR_indicies=c(17,228)
foldchanges=vector()
which(x[,1]=="PRAD")
which(x[,1]=="LUAD")
prad_indicies=1:73
luad_indicies=74:138
foldchanges=c(foldchanges,mean(x[prad_indicies,17])-mean(x[luad_indicies,17]))
foldchanges=c(foldchanges,mean(x[prad_indicies,228])-mean(x[luad_indicies,228]))
foldchanges_blood=foldchanges
table=cbind(table,foldchanges_blood)
colnames(table)
rownames(table)
table
str(table)
table
colnames(table)=c("miRNA","fold change: mean(prostate)/mean(lung) tissue", "log fold change: mean(prostate)-mean(lung) blood")
write.table(table,file="fold_changes.txt",sep="\t")
foldchanges
2^foldchanges






















setwd("C:\\Users\\mmoor\\OneDrive\\Desktop\\GDC_Data\\Counts_and_ColData/BRCA")
miR_data=read.table("miR_counts_matrix.txt",sep="\t",row.names=1)
col_data=read.table("col_data.txt",header=T,row.names=1)
x=col_data[,1]
which(col_data[,1]=="Normal")
ids=rownames(col_data)[which(col_data[,1]=="Normal")]
ids
substring_ids=substring(ids,1,12)
substring_ids
indicies=vector()
c="TCGA-E2-A153"
which(grepl(c,rownames(col_data)))
for(a in substring_ids)
{
    indicies=c(indicies,which(grepl(a,rownames(col_data))))
}
indicies
length(indicies)
x=miR_data
miR_data=miR_data[,indicies]
dim(miR_data)
y=x
miR_data=miR_data+1
x=col_data[indicies,1]
str(x)
row_1=as.vector(c(x[1],miR_data[,1]))
row_2=as.vector(c(x[2],miR_data[,2]))
table=rbind(row_1,row_2)
#combine the rest of the rows
for(i in 3:length(colnames(miR_data)))
{
  
  z=as.vector(c(x[i],miR_data[,i]))
  table=rbind(table,z)
}
miR_data=y
str(miR_data)
miR_data=as.numeric(miR_data)
