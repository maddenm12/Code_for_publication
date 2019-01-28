setwd("C:/Users/mmoor/OneDrive/Desktop/GDC_Data/luad_mirnaseq/Blood_Results/")
a=data.matrix(read.table("GSE61741table.txt",header=T,row.names=1))
tumor_samples=1:73
normal_samples=74:167
pvals=vector()
fold_changes=vector()
log_foldchanges=vector()
a=a[rowSums(a)>0,]
dim(a)
is.numeric(a)
for(i in 1:length(rownames(a)))
{
  pval=t.test(a[i,tumor_samples],a[i,normal_samples])$p.val
  fold_change=mean(a[i,tumor_samples])/mean(a[i,normal_samples])
  log_foldchange=log2(fold_change)
  pvals=c(pvals,pval)
  fold_changes=c(fold_changes,fold_change)
  log_foldchanges=c(log_foldchanges,log_foldchange)
}

pvals
fold_changes
log_foldchanges
table=cbind(pvals,fold_changes,log_foldchanges)
colnames(table)=c("P-value","Fold-Change","Log-Fold-Change")
rownames(table)=rownames(a)
old_miRs=rownames(a)[order(pvals)[1:25]]
write.table(old_miRs,file="old_miRs.txt",sep="\t")
b=read.table("miRBase_15_miR_accession.txt",fill=TRUE,stringsAsFactors=F)
dim(b)
accessions=vector()
for(rna in rownames(a))
{
  temp=which(b[,1]==paste(">",rna,sep=""))
  accession=b[temp,2]
  accessions=c(accessions,accession)
}
length(accessions)
miR_accessions=read.table("hsa_miR_accessionTOname_copy.txt",header=T,stringsAsFactors=F)
head(miR_accessions)
miR_names=vector()
for(i in 1:length(accessions))
{
  temp=accessions[i]
  temp2=which(miR_accessions[,1]==temp)
  if(length(temp2)==0)
  {
    miR_name="Empty"
  } else miR_name=miR_accessions[temp2,2]
  miR_names=c(miR_names,miR_name)
}
accessions[1]
accessions[5]
length(miR_names)
head(cbind(rownames(a),miR_names))
miR_names
empty=which(miR_names=="Empty")
table=table[-empty,]
miR_names=miR_names[-empty]
length(miR_names)
rownames(table)=miR_names
colnames(table)=c("P-value","Fold-Change","Log-Fold-Change")
which(grepl("miR-423",miR_names[order(pvals)[1:50]]))
top_25=miR_names[order(table[,1])[1:25]]
which(grepl("423-3p",top_25))
length(which(pvals<4.68*10^-7))
4.68*10^-7<10^-7
rf_miRs=read.table("luad_blood_miRs_new.txt",header=T)
dim(rf_miRs)
rf_miRs=rf_miRs[,1]
rf_miRs
miR_names2=vector()
for(i in 1:length(rf_miRs))
{
  temp=rf_miRs[i]
  temp2=which(miR_accessions[,1]==temp)
  if(length(temp2)==0)
  {
    miR_name="Empty"
  } else miR_name=miR_accessions[temp2,2]
  miR_names2=c(miR_names2,miR_name)
}
miR_names2
miR_names2=miR_names2[-which(miR_names2=="Empty")]
miR_names2
which(miR_names2=="Empty")
miR_names2
write.table(miR_names2,file="15_to_20_10outof10_miRs.txt",sep="\t")
orig_miR_names=read.table("most_predictive_miRs_10outof10.txt")[,1]
top_miRs_old=orig_miR_names[-c(6,22)]
write.table(top_miRs_old,file="miRs_both_15and20_blood_rf.txt",sep="\t")
getwd()
write.table(top_25,file="top_25_miRs.txt",sep="\t")



write.table(table,file="t-test_results.txt",sep="\t")

