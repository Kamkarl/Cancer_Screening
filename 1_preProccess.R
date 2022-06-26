setwd("D:/thesis/proposal/codes/R_Code/machine_learning/data")

stringsAsFactors = FALSE

library(randomForest)
library(tidyverse)
library(caret)
library(DMwR)
library(ROSE)
library(doParallel)

#### -------------------------------------------------
#### Load the data
#### -------------------------------------------------

# Inspect the data
# canList = c("Gastric Cancer", "Colorectal Cancer", "Esophageal Cancer")
canList = c("Gastric Cancer")
data1 <- read.csv("..\\result\\other_gastric.csv")
id = which(data1$classes %in% canList)
b <- data1[-id,]

# 40 is number of true in the list above.
index <- createDataPartition(b$classes, p = length(id)/nrow(b), list = FALSE)
finalData <- rbind(b[index,],data1[id,])

# Remove classes and replace it with other lable
c <- which(finalData$classes %in% canList)
finalData$classes <- as.character(finalData$classes)
finalData$classes[c] <- "cancer"
finalData$classes[-c] <- "other"
rownames(finalData) <- finalData[,1]
finalData <- finalData[,2:ncol(finalData)]
dim(finalData)
finalData$classes <- as.factor(finalData$classes)

#### -------------------------------------------------
#### Preproccess Data
#### -------------------------------------------------
# standard deviation
fdev <- sapply(finalData[,1:(ncol(finalData))-1], sd)
plot(fdev)
hist(fdev, main = "miRNAs standard deviation",xlab = "stDev")
### --------------------------------------------------

# create histograms for each attribute
par(mfrow=c(3,5))
for(i in 1:15) {
  hist(finalData[,i], main = "", xlab = names(finalData)[i])
}
title(main = "histogram of miRNAs")

# create density plot for each attribute
library(lattice)
par(mfrow=c(3,5))
for(i in 1:15) {
  plot(density(finalData[,i]), main=names(finalData)[i])
}

# Create separate boxplots for each attribute
par(mfrow=c(1,4))
for(i in 1:4) {
  boxplot(finalData[,i], main=names(finalData)[i])
}

# create a bar plot of each categorical attribute
par(mfrow=c(1,1))
for(i in 604:604) {
  counts <- table(finalData[,i])
  name <- names(finalData)[i]
  barplot(counts, main=name)
}

# create a missing map
library(Amelia)
library(mlbench)

missmap(finalData, col=c("black", "grey"), legend=FALSE)

### skewness -----------------------------------------
library(mlbench)
library(e1071)
skew <- apply(finalData[,1:(ncol(finalData))-1], 2, skewness)
# display skewness, larger/smaller deviations from 0 show more skew
plot(skew)

colnames(finalData[match(max(skew),skew)])
colnames(finalData[match(min(skew),skew)])
idp <- (skew < -0.6)
colnames(finalData[match(skew[id],skew)])
write.csv(colnames(finalData[match(skew[idp],skew)]),"../result/skewp.csv", quote = TRUE, row.names = FALSE)
write.csv(colnames(finalData[match(skew[idn],skew)]),"../result/skewn.csv", quote = TRUE, row.names = FALSE)
### --------------------------------------------------


# calculate correlations

correlations <- cor(finalData[,1:4])

### calc correlation matrix for numeric vars ---------
correlations <- cor(finalData[,1:(ncol(finalData))-1])
# display the correlation matrix
print(correlations)
plot(correlations,main = "corelation between miRNAs", xlab = 'miRNA', ylab = 'corelation')

# create correlation plot
library(corrplot)
corrplot(correlations[1:10,1:10], method="circle")
pid <- which(correlations[,] > 0.8 & correlations[,] != 1, arr.ind = TRUE)
nid <- which(correlations < -0.8, arr.ind = TRUE)

write_csv(as.data.frame(rownames(pid)),"../result/corr_p.csv") 
write_csv(as.data.frame(rownames(nid)),"../result/corr_n.csv")

### --------------------------------------------------
# pair-wise scatterplots of all 4 attributes
pairs(finalData[1:5])
# use class color to discriminate between samples
pairs(classes~., data=finalData[,600:604], col=finalData$classes)

### --------------------------------------------------
# density plots for each attribute by class value
x <- finalData[,1:6]
y <- finalData[,604]
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot="density", scales=scales)


### --------------------------------------------------
# box and whisker plots for each attribute by class value
x <- finalData[,1:6]
y <- finalData[,604]
featurePlot(x=x, y=y, plot="box")



### summarize data
summary(finalData)
# calculate the pre-process parameters from the dataset
preprocessParams <- preProcess(finalData, method=c("nzv"))
# summarize transform parameters
print(preprocessParams)
# transform the dataset using the parameters
transformed <- predict(preprocessParams, finalData)
# summarize the transformed dataset
summary(transformed)
### ----------------------------------------------------