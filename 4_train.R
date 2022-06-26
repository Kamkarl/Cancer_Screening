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
#canList = c("Non-cancer control")

#fname <- "other_gastric.csv"
#fname <- "otherWithNormal_gastric.csv" # Or blue_yellow.csv or other_gastric.csv
fname <- "impMir_CN.csv"

data1 <- read.csv(paste0("..\\result\\",fname))
id = which(data1$classes %in% canList)
b <- data1[-id,]

# 40 is number of true in the list above.
#index <- createDataPartition(b$classes, p = length(id)/nrow(b), list = FALSE)

# finalD <- rbind(b[index,],data1[id,])
# finalData <-finalD

finalData <- data1
finalData$classes

# create binary classification and remove classes and replace it with other lable
c <- which(finalData$classes %in% canList)
finalData$classes <- as.character(finalData$classes)
finalData$classes[c] <- "cancer"
finalData$classes[-c] <- "other"
rownames(finalData) <- finalData[,1]
finalData <- finalData[,2:ncol(finalData)]
dim(finalData)
finalData$classes <- as.factor(finalData$classes)
finalData$classes
#### -------------------------------------------------
#### attribute reduction or use all attributes
#### -------------------------------------------------
## All atributes
i1 <- 1:ncol(finalData) 

## List of attributes along with class column. count will be attrNum+1
pp <- c("hsa.miR.8073","hsa.miR.614","hsa.miR.548ah.5p","hsa.miR.1258")##,"hsa.miR.4536.3p","hsa.miR.4706","hsa.miR.6802.5p","hsa.miR.4705")##,"hsa.miR.4477b","hsa.miR.4464","hsa.miR.4719","hsa.miR.5591.5p","hsa.miR.4703.5p","hsa.miR.3118","hsa.miR.8076","hsa.miR.3607.5p","hsa.miR.208b.5p","hsa.miR.3686","hsa.miR.605.3p","hsa.miR.1277.3p","hsa.miR.3924","hsa.miR.548ai..hsa.miR.570.5p","hsa.miR.651.3p","hsa.miR.4696","hsa.miR.875.5p")
i1 <- which(colnames(finalData) %in% pp)
i1 <- c(i1,ncol(finalData)) ## Add class column
length(i1)
#### -------------------------------------------------
#### save final data based on selected attributes
#### -------------------------------------------------
finalData <- finalData[,i1]
save(finalData, file = "finalData.RData")

#### -------------------------------------------------
#### shuffle the data
#### -------------------------------------------------
set.seed(42)
shuffd <- sample_n(finalData, nrow(finalData))
shuffd$classes
dim(shuffd)


#### Split data
# index <- createDataPartition(bc_data$classes, p = 0.8, list = FALSE)
# train_data <- bc_data[index, ]
# test_data  <- bc_data[-index, ]

#### -------------------------------------------------
#### Train the model with different methods with K-Fold Cross Validation Approach
#### (method=rf, mlp, ... -- sampling = down, up, smote, rose)
#### -------------------------------------------------

metric <- "Accuracy" 
## OR
metric <-   "ROC" # or "Sens" or "Spec" whichever you desire

## Define training control
## classProbs=TRUE, summaryFunction=twoClassSummary are for ROC
## in accuracy case summaryFunction=twoClassSummary should be omitted from train.control

#### Repeated K-fold cross-validation
train.control <- trainControl(method = "repeatedcv", number = 5, repeats = 2)

#### K-fold cross-validation
set.seed(123) 
train.control <- trainControl(method = "cv", number = 5, classProbs=TRUE,
                        summaryFunction=twoClassSummary)

# RF: Random Forest
set.seed(7, sample.kind = "Rejection")
model.rf <- train(classes ~., data = shuffd,
               method = "rf", metric=metric, trControl = train.control)
# save model
save(model.rf,file = '../result/model.Rdata')
model.rf

# MLP: Multi Layer perceptron
set.seed(10, sample.kind = "Rejection")
model.mlp <- train(classes ~., data = shuffd[,], 
                  method = "mlp", metric=metric, trControl = train.control)
# LG: Logistic Regression
set.seed(10, sample.kind = "Rejection")
model.glm <- train(classes ~., data = shuffd[,], method="glm", metric=metric,
                 trControl=train.control)
# LDA: Linear Discriminate Analysis
set.seed(10, sample.kind = "Rejection")
model.lda <- train(classes ~., data = shuffd[,], method="lda", metric=metric,
                 trControl=train.control)

# GLMNET: Regularized Logistic Regression
set.seed(10, sample.kind = "Rejection")
model.glmnet <- train(classes ~., data = shuffd[,], method="glmnet", metric=metric,
                     trControl=train.control)
# KNN: k-Nearest Neighbors
set.seed(10, sample.kind = "Rejection")
model.knn <- train(classes ~., data = shuffd[,],  method="knn", metric=metric, 
                 trControl=train.control)
# CART: Classication and Regression Trees
set.seed(10, sample.kind = "Rejection")
model.cart <- train(classes ~., data = shuffd[,],  method="rpart", metric=metric,
                  trControl=train.control)
# NB: Naive Bayes
set.seed(10, sample.kind = "Rejection")
model.nb <- train(classes ~., data = shuffd[,],  method="nb", metric=metric,
                trControl=train.control)
# SVM: Support Vector Machines with Radial Basis Functions
set.seed(10, sample.kind = "Rejection")
model.svm <- train(classes ~., data = shuffd[,],  method="svmRadial", metric=metric,
                 trControl=train.control)
# Compare algorithms
mlist <- list(RF=model.rf, MLP=model.mlp , LR=model.glm, LDA=model.lda, GLMNET=model.glmnet, KNN=model.knn,
     CART=model.cart, SVM=model.svm)
transformResults <- resamples(mlist)
summary(transformResults)
dotplot(transformResults)
save(mlist,file = '../result/Allmodels.Rdata')

mlist

# Model$results

# Train the model with oob error
set.seed(10, sample.kind = "Rejection")
model <- train(classes ~., data = shuffd[,], 
               method = "rf", trControl = train.control, importance = TRUE)
model

#Tune the model parameters
tunegrid <- expand.grid(.mtry=20:50)
#### OR
x <- shuffd[,1:ncol(shuffd)-1]
y <- shuffd[,ncol(shuffd)]
set.seed(1)
bestMtry <- tuneRF(x, as.factor(y), stepFactor = 1.5, improve = 1e-5, ntree = 500) 
tunegrid <- expand.grid(.mtry=bestMtry[match(min(bestMtry[,2]),bestMtry[,2]),1])
model <- train(classes ~., data = shuffd[,]
                , method = "rf", trControl = train.control, tuneGrid=tunegrid, importance=TRUE)
model
#### dimension reduction and retrain

cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)

selected <- varImp(model, type=2, scale=FALSE)$importance
modellist <- list()

#train with different ntree parameters
for (j in c(5:50)){
  set.seed(123)
  s1 <- order(selected, decreasing = TRUE)[1:j]
  finalmir <- cbind(rownames(selected)[s1],selected$Overall[s1])
  model <- train(classes ~., data = shuffd[,c(s1,ncol(shuffd))]
                 , metric = 'Accuracy', method = "rf", trControl = train.control)
  key <- toString(j)
  modellist[[key]] <- model
}

#Compare results
save(modellist,file = "../result/modellist.Rdata")
results <- resamples(modellist)
summary(results)
plot(modellist)
model <- modellist$`39`

print(model)
plot(model)

## roc plot
library(ROCR)
pred <- prediction(model$finalModel$votes[,2],shuffd$classes)
perf <- performance(pred,"tpr","fpr")
plot(perf)#,colorize=TRUE)

## My summary report
model
model$results$Accuracy
confmat = confusionMatrix(model)
modelSummary = list(confmat$table,c(sub("%",max(model$results$Accuracy),"Accuracy: %"),
                                    sub("%", sensitivity(confmat$table),"Sensitivity: %"),
                                    sub("%", specificity(confmat$table),"Specificity: %")))
modelSummary

#### other code that may have some error
## run MLeval
library(MLeval)
res <- evalm(model)

## get ROC
res$roc

## get calibration curve
res$cc

## get precision recall gain curve
res$prg

#### Test many models

cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
#Manual search by create 10 folds and repeat 3 times
control <- trainControl(method = 'repeatedcv',
                        number = 5,
                        repeats = 2,
                        search = 'grid')
#create tunegrid
tunegrid <- expand.grid(.mtry = c(sqrt(ncol(shuffd))))
modellist <- list()

#train with different ntree parameters
for (ntree in c(200,300)){
  set.seed(123)
  fit <- train(classes ~.,
               data = shuffd,
               method = 'rf',
               metric = 'Accuracy',
               tuneGrid = tunegrid,
               trControl = control,
               ntree = ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

#Compare results
results <- resamples(modellist)
summary(results)
plot(results)

### --------------------------
##  Tuning glmnet
### --------------------------
myGrid <- expand.grid(
  alpha = seq(0, 1, length = 3),
  lambda = seq(0.0001, 0.1, length = 10)
)

myGrid <- expand.grid(
  alpha = 1,
  lambda = 0.01
)


set.seed(7)
model.glmnet <- train(classes ~., data = shuffd[,], method="glmnet", metric=metric,
                      trControl=train.control, tuneGrid=myGrid)
plot(model.glmnet)
model.glmnet[["results"]]
model.glmnet[["results"]][["ROC"]]

