getwd()
setwd("D:/thesis/proposal/codes/R_Code/machine_learning/data")
options(stringsAsFactors = FALSE)
####-----------------------------------------------
#### To sample from each class type equally to be balanced with selected type
####-----------------------------------------------
balanceDataClass<-function(data,refClass,class){
  ## Maybe it is better to replace == with %in%
  i1 <- which(data$class%in%refClass)
  d2 <- data[i1,]
  num <- length(i1)
  print(num)
  #i2 <- sample(i1,40)
  others <- unique(data$class)
  n <- abs(num/(length(others)-1))+1
  print(n)
  for (i in others) {
    print(i)
    if(!(i %in% refClass)){
      i2 <- which(data$class==i)
      i2 <- sample(i2,n)
      d2 <- rbind(d2,f1[i2,])
    }
  }
  d2
}




ll1 <- load("cancer-healthy-01-dataInput.RData")
ll2 <- load("cancer-healthy-02-dataInput.RData")
ll3 <- load('ex-gr-01-dataInput.RData')
#### ------------------------------------------
#### Use selected miRNAs
#### ------------------------------------------
mirList<- read.table("impMir.txt")
mod1 <- which(colnames(datExpr) %in% mirList$miRNA)
mcolor <- "impMir" ## for file name

#### ------------------------------------------
#### or use selected modules
#### ------------------------------------------
mcolor <- c("blue","yellow")
mod1 <- which(moduleColors %in% mcolor)

#### ------------------------------------------
####  write selected miRNA exp in a file along with disease
#### ------------------------------------------
id <- rownames(datExpr)
tmp2 <- t(rbind(colnames(ex),disease))
dtype <- match(id,tmp2[,1])
classes <- disease[dtype]
f1 <- cbind(id,datExpr[,mod1],classes)
####-------------------------------------------
#### balance classes
#### in case of multiple refs two or more ref classes if size of ref classes is not the same 
#### by calling the function on reduced data and the one ref with min size it can be balanced first.
#### then all of them become balanced with other classes. since we are doing binary classification.
#### in current problem it s not needed.
####-------------------------------------------
ref <- c("Non-cancer control")
bb <- balanceDataClass(f1,ref,classes)

####Write the file with id added
write.csv(bb,paste0("..\\result\\",paste(mcolor,collapse = "_"),"_CN.csv"),row.names = FALSE)

#### ------------------------------------------
#### Write the comparison file for selected module and cancer via other cancers
#### ------------------------------------------
# "Biliary Tract Cancer"     "Bladder Cancer"           "Breast Cancer"           
# "Colorectal Cancer"        "Esophageal Cancer"        "Gastric Cancer"          
# "Glioma"                   "Hepatocellular Carcinoma" "Lung Cancer"             
# "Non-cancer control"       "Ovarian Cancer"           "Pancreatic Cancer"       
# "Prostate Cancer"          "Sarcoma"  
#temp1<- gr=="BiliaryTractCancer" | gr=="ColorectalCancer" | gr== "EsophagealCancer" | gr=="GastricCancer" | gr=="HepatocellularCarcinoma" | gr=="PancreaticCancer"

# Unbalanced classes are Bladder and Non-cancer control
temp2 <- which(disease=="Bladder Cancer")
a <- sample(temp2,40)

temp1 <- which(disease=="Breast Cancer" | disease== "Glioma" |
                    disease=="Lung Cancer" | disease=="Ovarian Cancer" |
                    disease=="Prostate Cancer" | disease=="Sarcoma")
temp1 <- c(temp1,a)
d1 <- disease[temp1]
f2 <- rbind(t(ex[mod1,temp1]))
#class2 <- rep("others",dim(f2)[1])
f2<- cbind(rownames(f2),f2)
colnames(f2)[1] <- "id"
f2<- cbind(f2,d1)
colnames(f2)[dim(f2)[2]] <- "classes"
# write control data along with other cancers into a file
write.csv(rbind(f1,f2),"..\\result\\otherWithNormal_gastric.csv",row.names = FALSE)

# Remove control data from first file
sel1 <- which(f1$classes != "Non-cancer control")
write.csv(rbind(f1[sel1,],f2),"..\\result\\other_gastric.csv",row.names = FALSE)
