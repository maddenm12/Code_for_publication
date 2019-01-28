b=read.table("miRBase_15_miR_accession.txt",fill=TRUE,stringsAsFactors=F)
dim(b)
accessions=vector()
old_miRs=read.table("old_miRs.txt",header=T,stringsAsFactors=F)[,1]

for(rna in old_miRs)
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
miR_names
old_miRs=old_miRs[-which(miR_names=="Empty")]

miR_names=miR_names[-which(miR_names=="Empty")]
miR_names
old_miRs
write.table(old_miRs,file="miRs_both_15and20_blood_t-test.txt",sep="\t")
write.table(miR_names,file="15_to_20_t-test_miRs.txt",sep="\t")
