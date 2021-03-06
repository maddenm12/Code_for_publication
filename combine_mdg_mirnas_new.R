full_list=vector()
setwd("C://Users//mmoor//OneDrive//Desktop//GDC_Data//prad_miRnaseq/Tissue_Results_RPM/")
trial_1=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_1.txt",header=T,row.names=1)
miRs_1=rownames(trial_1)
full_list=miRs_1
trial_2=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_2.txt",header=T,row.names=1)
miRs_2=rownames(trial_2)
full_list=c(full_list,setdiff(miRs_2,full_list))
trial_3=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_3.txt",header=T,row.names=1)
miRs_3=rownames(trial_3)
full_list=c(full_list,setdiff(miRs_3,full_list))
#length(full_list)
trial_4=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_4.txt",header=T,row.names=1)
miRs_4=rownames(trial_4)
full_list=c(full_list,setdiff(miRs_4,full_list))
trial_5=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_5.txt",header=T,row.names=1)
miRs_5=rownames(trial_5)
full_list=c(full_list,setdiff(miRs_5,full_list))
trial_6=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_6.txt",header=T,row.names=1)
miRs_6=rownames(trial_6)
full_list=c(full_list,setdiff(miRs_6,full_list))
trial_7=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_7.txt",header=T,row.names=1)
miRs_7=rownames(trial_7)
full_list=c(full_list,setdiff(miRs_7,full_list))
trial_8=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_8.txt",header=T,row.names=1)
miRs_8=rownames(trial_8)
full_list=c(full_list,setdiff(miRs_8,full_list))
trial_9=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_9.txt",header=T,row.names=1)
miRs_9=rownames(trial_9)
full_list=c(full_list,setdiff(miRs_9,full_list))
trial_10=read.table("top_100_miRs_MeanDecreaseGini_prad_trial_10.txt",header=T,row.names=1)
miRs_10=rownames(trial_10)
full_list=c(full_list,setdiff(miRs_10,full_list))
length(full_list)
head(full_list)
miRs=read.table("preprocessed_data_firehose_prad_rpm.txt",header=T)
all_miRs=colnames(miRs)[2:length(colnames(miRs))]
length(all_miRs)
head(all_miRs)
#all_miRs=gsub("-",".",all_miRs)
all_significant_miRs=cbind(trial_1,trial_2,trial_3,trial_4,trial_5,trial_6,trial_7,trial_8,trial_9,trial_10)
dim(all_significant_miRs)
scores=vector()
for(i in 1:length(all_miRs))
{
  score=0
  for(j in 1:length(colnames(all_significant_miRs)))
  {
    string=paste("^",all_miRs[i],"$",sep="")
    score=score+length(grep(string,all_significant_miRs[,j]))
    scores[i]=score
  }
test=scores[which(scores>6)]
test[5]
test[6]
test[18]
all_miRs[which(scores>4)]
}
length(which(scores==10))
miRs_appear_all_times=all_miRs[which(scores==10)]
#lengths=vector(
length(miRs_appear_all_times)
library(stringi)
new_list=vector()
for(i in 1:length(miRs_appear_all_times))
{
  
  a=stri_locate_all(pattern = 'MIMAT', miRs_appear_all_times[i], fixed = TRUE)
  a=unlist(a)[1]
  temp=substr(miRs_appear_all_times[i],1,a-2)
  new_list=c(new_list,temp)
}
write.table(new_list,file="most_predictive_miRs_10outof10.txt",sep="\t")
#length(all_miRs)
#length(full_list)
a=as.matrix(test_network_rf)
write.table(a,file="test_network.txt",sep="\t")
getwd()
