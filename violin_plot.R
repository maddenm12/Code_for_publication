library(ggplot2)

all_miRs_results=c(.998,1,1,.9945,.9983,.9997,.9988,1,.9959,1)
miRs_results_22=c(1,.9987,1,.993,1,.999,1,1,.9984,1)
miRs_results_3=c(.9995,.9899,.994,.997,.9861,.9515,.9948,.9994,.9512,.987)

ovr_results=c(all_miRs_results,miRs_results_22,miRs_results_3)
  miRs=c(rep("All miRNAs",10),rep("top 51 miRNAs",10),rep("z (overlap)",10))
results=as.data.frame(cbind(as.numeric(ovr_results),miRs))
results[,1]=as.numeric(ovr_results)
str(results)
colnames(results)=c("Area_under_ROC_Curve","Number_of_miRNAs")
results
p <- ggplot(results, aes(x=Number_of_miRNAs, y=Area_under_ROC_Curve,fill=Number_of_miRNAs) ) + geom_violin(trim=TRUE) + coord_cartesian(ylim = c(0.95, 1))+ geom_boxplot(width=0.05) + theme(legend.position="none")+labs(title="Lung Adenocarcinoma Blood Data",x="Number of miRNAs", y = "Area Under ROC Curve")
p

#BRCA,HNSC,LUAD,UCEC
#KIRC,KIRP,PRAD,THCA
ovr_results=c(KIRC_Results,KIRP_Results,PRAD_Results,THCA_Results)
cancer=c(rep("KIRC",10),rep("KIRP",10),rep("PRAD",10),rep("THCA",10))
results=as.data.frame(cbind(as.numeric(ovr_results),cancer))
results[,1]=as.numeric(ovr_results)
colnames(results)=c("Accuracy","Cancer")
results
p <- ggplot(results, aes(x=Cancer, y=Accuracy,fill=Cancer) ) + geom_violin(trim=TRUE) + coord_cartesian(ylim = c(0.75, 1))+ geom_boxplot(width=0.1) + theme(legend.position="none")+labs(title="Random Forest Classification Accuracy",x="Cancer", y = "Prediction Accuracy")
p
help(mean_sdl)
