setwd("D:/thesis/proposal/codes/R_Code/machine_learning/data")
stringsAsFactors=FALSE

library(tidyverse)
library(caret)

bc_data <- read.csv("..\\result\\test.csv")
row.names(bc_data) <- bc_data[,1]
bc_data <- bc_data[,2:dim(bc_data)[2]]
summary(bc_data$classes)

set.seed(42)
shuffd <- sample_n(data1, dim(data1)[1])

index <- createDataPartition(bc_data$classes, p = 0.8, list = FALSE)

train_data <- bc_data[index, ]
test_data  <- bc_data[-index, ]

#### ---------------------------------------
#### Modeling the original unbalanced data
#### ---------------------------------------

set.seed(42)
# model_rf1 <- caret::train(classes ~ .,
#                          data = train_data,
#                          method = "mlp",
#                          preProcess = c("scale", "center"),
#                          trControl = trainControl(method = "repeatedcv", 
#                          number = 10, 
#                          repeats = 10, 
#                          verboseIter = FALSE))

model_rf <- caret::train(classes ~ .,
      data = train_data,
      method = "mlp",
      preProcess = c("scale", "center"),
      trControl = trainControl(method = "repeatedcv", 
                                          number = 5, 
                                          repeats = 1, 
                                          verboseIter = FALSE))


#### ----------------------------------------------------
#### Under-sampling or Up-sampling =up, down, smote, rose  
#### ---------------------------------------
ctrl <- trainControl(method = "repeatedcv", 
                     number = 5, 
                     repeats = 1, 
                     verboseIter = FALSE,
                     sampling = "down")

set.seed(42)
model_rf <- caret::train(classes ~ .,
                               data = train_data,
                               method = "rf",
                               preProcess = c("scale", "center"),
                               trControl = ctrl)


#### ---------------------------------------
#### Test Data
#### ---------------------------------------

final <- data.frame(actual = test_data$classes,
                    predict(model_rf, newdata = test_data, type = "prob"))
final$predict <- ifelse(final$Gastric.Cancer > 0.5, "Gastric Cancer", "Non-cancer control")
cm_original <- confusionMatrix(as.factor(final$predict), test_data$classes)
cm_original
