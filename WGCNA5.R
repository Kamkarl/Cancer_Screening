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
###########Exporting to Cytoscape##############
###############################################
# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("blue", "turquoise");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
#                                weighted = TRUE,
#                                threshold = 0.02,
#                                nodeNames = modProbes,
#                                altNodeNames = modGenes,
#                                nodeAttr = moduleColors[inModule]);
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("..\\result\\CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("..\\result\\CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])
