setwd("D:/thesis/proposal/codes/R_Code/machine_learning/data")
options(stringsAsFactors = FALSE);

library(randomForest)
library(tidyverse)
library(DMwR)
library(ROSE)
library(doParallel)
#library(glmnet)
#library(RSNNS)
library(caret)
library(dplyr)


### plot result of model test based on input test data and original disease type.
### By default the last column is refClass 
### returns confusion matrix
testResultPlot<-function(testData,model,origin){
  a <- predict(model,testData)
  colnames(testData)
  refClass <- testData[,ncol(testData)]
  
  b <- confusionMatrix(a,as.factor(refClass))
  c <- which(a != refClass)
  
  misClass <- table(origin[c])
  # plt1 <- barplot(misClass, xaxt="n",main = 'Mis classified', ylab = "Frequency")
  # text(plt1, par("usr")[3], labels = row.names(misClass),
  #      srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6)
  
  riClass <- table(origin[-c])
  # plt2 <- barplot(riClass, xaxt="n", main = 'Correctly classified', ylab = "Frequency")
  # text(plt2, par("usr")[3], labels = row.names(riClass),
  #      srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6)
  
  a <- merge(riClass,misClass, by = 'Var1', all = TRUE)
  ## important: replace NA values with zero
  a[is.na(a)] <- 0
  #b <- merge(riClass,misClass, all = TRUE)
  xp <- cbind(a$Freq.x,a$Freq.y)
  rownames(xp) <- a$Var1
  colnames(xp) <- c('correct','mis')
  yy <- a$Freq.y+a$Freq.x
  plt3 <- barplot(t(xp), xaxt="n", main = "Model test Result", ylim = range(0,0,1.25*max(yy)),
                  col=c("darkolivegreen3","firebrick3"), ylab = "Frequency", legend.text = TRUE,
                  args.legend = list(x ='topright', bty='n', inset=c(-0.10,-0.15),cex=0.7) )
  text(plt3, par("usr")[3], labels = row.names(xp),
       srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6)
  text(plt3, yy+80, signif(a$Freq.x/(yy),2),cex=0.6, srt = 90)
  
  b
}

ll <- load("finalData.RData")
l2 <- load("testData.RData")
l3 <- load('../result/model.Rdata')

## predict based on original data
rf <- model.rf
a <- predict(rf,finalData)
b <- confusionMatrix(a,as.factor(finalData$classes))
b

### attribute reduction for test data
### Test with same attributes as model data 
i2 <- which(colnames(testData) %in% colnames(finalData))

### remove unreliable data like nagative prostate biopsy
x <- testData$classes != "Negative prostate biopsy"
testData <- testData[x,]


canList = c("non-Cancer")
testD <- testData[,i2]
r2 <- which(testD[,ncol(testD)] %in% canList)
testD$classes[r2] <- "cancer"
testD$classes[-r2] <- "other"
colnames(testD)

## predict based on test data and plot result
testResultPlot(testD, model = rf, origin = testData$classes)

