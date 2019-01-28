setwd("C:\\Users\\mmoor\\OneDrive\\Desktop\\GDC_Data\\prad_mirnaseq")
dir()
x=read.table("prad_tissue_rpm.txt",header=T,row.names=1)
dim(x)
x[is.na(x)]=0
which(is.na(x))
miR_data=x
dim(miR_data)
tumor_samples=c(which(grepl(".01A.",colnames(miR_data))),which(grepl(".01B.",colnames(miR_data))),which(grepl(".02A.",colnames(miR_data))),which(grepl(".06A.",colnames(miR_data))))
length(tumor_samples)
#ratio=498 tumor 46 control
control_samples=c(which(grepl(".11A.",colnames(miR_data))),which(grepl(".11B.",colnames(miR_data))))
tumor_samples
control_samples
#length(tumor_samples)+length(control_samples)==length(colnames(miR_data))
x=vector()
x=rep("",length(colnames(miR_data)))
x[tumor_samples]="Tumor"
x[control_samples]="Normal"
x=as.factor(x)
row_1=as.vector(c(x[1],miR_data[,1]))
row_2=as.vector(c(x[2],miR_data[,2]))
table=rbind(row_1,row_2)
str(table)
#combine the rest of the rows
for(i in 3:length(colnames(miR_data)))
{
  
  z=as.vector(c(x[i],miR_data[,i]))
  table=rbind(table,z)
}
#Adjust the column names and row names
colnames(table)=c("Condition",rownames(miR_data))
rownames(table)=c()
table=as.data.frame(table)
table[,1]=x
str(table)
write.table(table,file="preprocessed_data_firehose_prad_rpm.txt",sep="\t",row.names=F)
