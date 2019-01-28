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
#set.seed(123) 
#Set a random seed each time the program is run, comment out if you want a set seed every time
time <- as.numeric(Sys.time())
seed2 <- 1e8 * (time - floor(time))
set.seed(seed2)
#read in the data
setwd("C:\\users\\mmoor\\OneDrive\\Desktop\\GDC_Data\\Counts_and_ColData\\brca")
data=as.data.frame(read.table("random_forest_preprocessed_data_brca2.txt",header=T))
#Create training and testing datasets
sample.ind=sample(2,
                  nrow(data),
                  replace=T,
                  prob=c(0.9,0.1))
train=data[sample.ind==1,]
test=data[sample.ind==2,]
#Run the random forest
rf=randomForest(Condition ~ ., 
                ntree = 500,
                data = train,
                nodesize=100
                #classwt = c(0.1, 0.9)
)
#train=data[sample.ind==1,]
#test=data[sample.ind==2,]
#Run the random forest
#rf=randomForest(Condition ~ ., 
#                ntree = 100,
#                data = train
#)
#scale = false: Don't divide by standard error
var.imp = data.frame(importance(rf,
                                type=2))


# make row names as columns
var.imp$Variables <- row.names(var.imp)
#order(var.imp$MeanDecreaseGini,decreasing=T)
important_miRs=var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
write.table(important_miRs[1:100,],file="top_100_miRs_MeanDecreaseGini_ucec_trial_20.txt",sep="\t")
#see how accurate the predictions are with the trained data set (theoretically should be 100%)
train$predicted.response <- predict(rf , train)
confusionMatrix(data = train$predicted.response,
              reference = train$Condition,
                positive = 'Normal')
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
plot(rf.performance,main="ROC Curve for Random Forest Breast Cancer",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
capture.output(a,file="confusion_matrix_ucec_trial_20.txt")
dev.copy(png,'ucec_roc_10.png')
dev.off()
a
options(scipen=5)

