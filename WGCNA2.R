getwd();
setwd("D:/thesis/proposal/codes/R_Code/wgcna/data");
# Load the WGCNA package
library(WGCNA)
library(flashClust)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "..\\result\\cancer-healthy-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

###############################################
#####Choosing the soft-thresholding power######
###############################################
# Choose a set of soft-thresholding powers
powers = c(c(1:15), seq(from = 16, to=20, by=2))
netType <- "unsigned" # signed or unsigned or signed hybrid
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = netType, verbose = 5)

# Plot the results:
sizeGrWindow(10,8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab=c("Scale Free Topology Model Fit, R^2",netType),type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$fitIndices[,2]>0.9

###############################################
#####Co-expression similarity and adjacency####
###############################################
tmp1<--sign(sft$fitIndices[,3])*sft$fitIndices[,2]>0.9
softPower = powers[match(TRUE,tmp1)];
adjacency = adjacency(datExpr, type = netType, power = softPower);


###############################################
############Topological Overlap Matrix#########
###############################################
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

###############################################
############Clustering using TOM###############
###############################################
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "miRNA clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#table(dynamicMods)
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "miRNA dendrogram and module colors")

###############################################
##############Merging of modules###############
###############################################
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.02
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, TOM, file = "..\\result\\cancer-healthy-02-dataInput.RData")
#save requires data for module preservation step like:module colors and other info
mPreservationFile <- cbind(t(datExpr),moduleColors)
write.csv(mPreservationFile,file= "../result/3cancerModules.csv")
