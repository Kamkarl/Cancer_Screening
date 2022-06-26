getwd();
workingDir = "D:/thesis/proposal/codes/R_Code/wgcna/data";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "..\\result\\cancer-healthy-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "..\\result\\cancer-healthy-02-dataInput.RData");
lnames


###############################################
####Quantifying module-trait associations######
###############################################
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


###############################################
####Gene Significance and Module Membership####
###############################################
# Define variable cancer containing the cancer column of datTrait
#cancer = as.data.frame(datTraits$cancer);
cancer = as.data.frame(datTraits$GastricCan);
names(cancer) = "cancer"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, cancer, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(cancer), sep="");
names(GSPvalue) = paste("p.GS.", names(cancer), sep="");

###############################################
###########Intramodular analysis###############
###############################################
sizeGrWindow(15, 100);
par(mfrow = c(3,ceiling(length(modNames)/3)));
#### my code
impMir <- matrix(nrow = 0, ncol = 4)
colnames(impMir) <- c("module","miRNA","geneModuleMembership","geneTraitSignificance")
#### end of my code

for (module in  modNames){
  print(module)
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
#### my code
  miRNames = colnames(datExpr)[moduleGenes]
  myGMM = abs(geneModuleMembership[moduleGenes, column])
  myGTS = abs(geneTraitSignificance[moduleGenes, 1])
  myIndex = myGMM > 0.7 & myGTS > 0.5
  a <- matrix(nrow = 0, ncol = 4)
  a <- cbind(module,miRNames[myIndex],myGMM[myIndex],myGTS[myIndex])
  if(ncol(a)>1)
    impMir <- rbind(impMir,a)
  
#### end of my code  
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste(" Membership in", module),
                   ylab = "Gene sig for cancer",
                   main = paste(" \n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
write.table(impMir,"impMir.txt",sep = "\t", row.names = TRUE, col.names = TRUE)
###############################################
###########Summary output######################
###############################################
names(datExpr)
#names(datExpr)[moduleColors=="red"]
#annot = read.csv(file = "GeneAnnotation.csv");
# dim(annot)
# names(annot)
probes = names(datExpr)
# probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
# geneInfo0 = data.frame(substanceBXH = probes,
#                        geneSymbol = annot$gene_symbol[probes2annot],
#                        LocusLinkID = annot$LocusLinkID[probes2annot],
#                        moduleColor = moduleColors,
#                        geneTraitSignificance,
#                        GSPvalue)
geneInfo0 = data.frame(substanceBXH = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for cancer
modOrder = order(-abs(cor(MEs, cancer, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.cancer));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "..\\result\\geneInfo.csv")
