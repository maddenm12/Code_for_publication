#Library required libraries and set seed
library(randomForest)
library(e1071)
library(caret)
library(ggplot2)
library(lattice)
library(chron)
library(ROCR)
#If you want the results to be the same every time (mainy for debugging), uncomment
#this line to have a set seed each time
 
#Set a random seed each time the program is run, comment out if you want a set seed every time
time <- as.numeric(Sys.time())
seed2 <- 1e8 * (time - floor(time))
set.seed(seed2)
#read in the data
#setwd("C:\\users\\mmoor\\OneDrive\\Desktop\\GDC_Data\\PRAD_Blood/")
setwd("C:\\Users\\mmoor\\OneDrive\\Desktop\\GDC_Data/luad_miRnaseq/")
dir()
dir()
data=as.data.frame(read.table("preprocessed_data_firehose_luad_rpm.txt",header=T))
#data[1:5,1]
#dir()
#Create training and testing datasets
num_tumor=length(which(data[,1]=="Tumor"))
num_control = length(which(data[,1]=="Normal"))
#tumor=data[which(data[,1]=="Tumor"),]
#control=data[which(data[,1]=="Normal"),]
#difference = num_tumor-num_control
#adjusted_control=control
#sample.control=sample(1:num_control,difference,replace=T)
#adjusted_control=rbind(adjusted_control,control[sample.control,])
sample.ind = sample(2, 
                    nrow(data),
                    replace = T,
                    prob = c(0.7,0.3))
#train_tumor = tumor[sample.ind==1,]
#test_tumor = tumor[sample.ind==2,]
#train_control= adjusted_control[sample.ind==1,]
#test_control= adjusted_control[sample.ind==2,]
#train=rbind(train_tumor,train_control)
#test=rbind(test_tumor,test_control)
train=data[sample.ind==1,]
test=data[sample.ind==2,]
#Run the random forest
rf=randomForest(Condition ~ ., 
                ntree = 500,
                data = train
)
#scale = false: Don't divide by standard error
var.imp = data.frame(importance(rf,
                                type=2))


 #make row names as columns
var.imp$Variables <- row.names(var.imp)
#order(var.imp$MeanDecreaseGini,decreasing=T)
important_miRs=var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
setwd("C:\\Users\\mmoor\\OneDrive\\Desktop\\GDC_Data\\luad_mirnaseq/Tissue_Results_RPM/")
write.table(important_miRs[1:100,],file="top_100_miRs_MeanDecreaseGini_luad_trial_10.txt",sep="\t")
#see how accurate the predictions are with the trained data set (theoretically should be 100%)
train$predicted.response <- predict(rf , train)
#confusionMatrix(data = train$predicted.response,
#               reference = train$Condition,
#                positive = 'Normal')
#See how accurate the predictions are with the test dataset, write to a file
test$predicted.response <- predict(rf ,test)
test.predicted=predict(rf,test)
a=confusionMatrix(data = test$predicted.response,
                reference = test$Condition,
                positive = 'Tumor')

#predictions=as.numeric(predictions)
#lables=as.numeric(labels)
rf.predict.2=predict(rf,newdata=test,type="prob")
rf.predict=prediction(rf.predict.2[,2],test$Condition)
rf.performance = performance(rf.predict,"tpr","fpr")
plot(rf.performance,main="ROC Curve for Lung Adenocarcinoma",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
auc.perf=performance(rf.predict,measure="auc")
capture.output(a,file="luad_cm_trial_10.txt")
a
auc.perf@y.values
