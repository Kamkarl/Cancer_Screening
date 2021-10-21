###############################################
###################Installation################
###############################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
BiocManager::install(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
BiocManager::install("WGCNA")

# Load the package
library(WGCNA)
library(flashClust)
###############################################
########Loading expression data################
###############################################

getwd()
workingDir = "D:/thesis/proposal/codes/R_Code/wgcna/data/"
setwd(workingDir)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Read in the cancer data set , ex contains expression data
lnames = load(file = "ex-gr-01-dataInput.RData")
levels(as.factor(disease))



#### for module preservation analysis
otherCancer <- cbind(ex[1:2565,1:392],ex[1:2565,493:532],ex[1:2565,693:732],ex[1:2565,773:812],ex[1:2565,813:852],ex[1:2565,893:932],ex[1:2565,933:972])
write.csv(otherCancer, "othercancer.csv")

#### Define groups for which network will be constructed.

# temp1 <- disease =="Esophageal Cancer" | disease=="Gastric Cancer" | disease=="Colorectal Cancer"
# temp1 <- disease!="Non-cancer control"
### All other cancer types
temp1 <- disease != "Gastric Cancer" & disease != "Non-cancer control"
g1 <- ex[,temp1]

# temp2 <- disease=="Non-cancer control"
temp2 <- disease=="Gastric Cancer"
g2 <- ex[,temp2]
all <- cbind(g2,g1)
datExpr0 = as.data.frame(t(all))
#### CheckUp of dimension and normalization
dim(g1)
names(datExpr0)
rownames(datExpr0)
boxplot(t(all[1:20,1:20]),las=2)
#############



gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(50,20)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
cut1 <- 250
abline(h = cut1, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut1, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
gr2 <- c(disease[temp1],disease[temp2])
gr2 = gr2[keepSamples]

###############################################
########Loading clinical trait data############
###############################################

traitData = read.csv("ClinicalTraits2.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(23:36)];
##allTraits = allTraits[, c(2, 13:16) ]; #"Series_sample_id","Sex","Age","cancer","GasterCancers"  
allTraits = allTraits[, c(2, 13:14, 18)]
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
cancerSamples = rownames(datExpr);
traitRows = match(cancerSamples, allTraits$Series_sample_id);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
# repalce gr2 labeles for better representation
gr2 <- gsub(" Cancer", "", gr2)
gr2 <- gsub(" Carcinoma", "", gr2)
label1 <-levels(as.factor(gr2))
traitColors = numbers2colors(as.numeric(as.factor(gr2)), signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf()
plotDendroAndColors(sampleTree2, traitColors, groupLabels=length(label1),
                    rowText = gr2, rowTextAlignment = "left", addTextGuide = TRUE,
                    main = c("Dendrogram of trait heatmap of", length(label1)-1, "cancers"))
dev.off()

save(datExpr, datTraits, file = "..\\result\\cancer-healthy-01-dataInput.RData")
save(datExpr, file = "..\\result\\datExpr.RData")
collectGarbage()


