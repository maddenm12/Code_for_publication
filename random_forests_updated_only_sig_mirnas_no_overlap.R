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
  data=as.data.frame(read.table("random_forest_preprocessed_data_brca.txt",header=T))
  top_miRs = as.vector(read.table("combined_mdg_miRs_7.txt"))[,1]
  columns=vector()
  for(i in 1:length(colnames(data)))
  {
    for(j in 1:length(top_miRs))
    {
      if(colnames(data)[i]==top_miRs[j])
      {
        columns=c(columns,i)
      }
    }
  }
  data=data[,c(1,columns)]
  #Create training and testing datasets
  #num_control = length(which(data[,1]=="Normal"))
  #tumor=data[which(data[,1]=="Tumor"),]
  #control=data[which(data[,1]=="Normal"),]
  sample.ind=sample(2,
                    nrow(data),
                    replace=T,
                    prob=c(0.60,0.40))
 # train_tumor=tumor[sample.ind==1,]
#  test_tumor=tumor[sample.ind==2,]
#  sample.ind=sample(2,
 #                   nrow(control),
  #                  replace=T,
   #                 prob=c(0.60,0.40))
  #train_control=control[sample.ind==1,]
  #test_control=control[sample.ind==2,]
  #num_control=length(rownames(train_control))
  #num_tumor=length(rownames(train_tumor))
  #difference = num_tumor-num_control
  #adjusted_control=train_control
  #sample.control=sample(1:num_control,difference,replace=T)
  #adjusted_control=rbind(adjusted_control,train_control[sample.control,])
  #train=rbind(train_tumor,adjusted_control)
  #test=rbind(test_tumor,test_control)
  train=data[sample.ind==1,]
  test=data[sample.ind==2,]
  #Run the random forest
  rf=randomForest(Condition ~ ., 
                  ntree = 100,
                  data = train#,
                  #classwt = c(0.1, 0.9)
  )
  #scale = false: Don't divide by standard error
  var.imp = data.frame(importance(rf,
                                  type=2))
  
  
  # make row names as columns
  var.imp$Variables <- row.names(var.imp)
  #order(var.imp$MeanDecreaseGini,decreasing=T)
  important_miRs=var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
  #write.table(important_miRs[1:100,],file="top_100_miRs_MeanDecreaseGini_ucec_model_1.txt",sep="\t")
  #see how accurate the predictions are with the trained data set (theoretically should be 100%)
  train$predicted.response <- predict(rf , train)
  #confusionMatrix(data = train$predicted.response,
  #               reference = train$Condition,
  #              positive = 'Tumor')
  #See how accurate the predictions are with the test dataset, write to a file
  test$predicted.response <- predict(rf ,test)
  test.predicted=predict(rf,test)
  a=confusionMatrix(data = test$predicted.response,
                    reference = test$Condition,
                    positive = 'Tumor'
                  )
  
  #predictions=as.numeric(predictions)
  #lables=as.numeric(labels)
  rf.predict.2=predict(rf,newdata=test,type="prob")
  rf.predict=prediction(rf.predict.2[,2],test$Condition)
  rf.performance = performance(rf.predict,"tpr","fpr")
  plot(rf.performance,main="ROC Curve for Random Forest Uterine Corpus Endometrial Carcinoma",col=2,lwd=2)
  abline(a=0,b=1,lwd=2,lty=2,col="gray")
  capture.output(a,file="confusion_matrix_ucec_top_2_trial_10.txt")
  dev.copy(png,'ucec_top_2_roc_trial_10.png')
  dev.off()
  a
predictive_miRs[1:20,2]
x=vector()
for(i in 1:length(test$predicted.response))
{
  if(test$predicted.response[i]=="Tumor")
  {
    x=c(x,1)
  } else x=c(x,0)
}
y=vector()
for(i in 1:length(test$Condition))
{
  if(test$Condition[i]=="Tumor")
  {
    y=c(y,1)
  } else y=c(y,0)
}
act=y
pred=x
mccr(act,pred)
