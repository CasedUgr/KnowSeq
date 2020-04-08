library(KnowSeq)

# Downloading one series from NCBI/GEO and one series from ArrayExpress
downloadPublicSeries(c("GSE74251","GSE81593"))

# Using read.csv for NCBI/GEO files (read.csv2 for ArrayExpress files)
GSE74251csv <- read.csv("ReferenceFiles/GSE74251.csv")
GSE81593csv <- read.csv("ReferenceFiles/GSE81593.csv")

# Loading the transcripts to genes converter variable
dir <- system.file("extdata", package="tximportData")
tx2gene <- read.csv(file.path(dir, "tx2gene.ensembl.v87.csv"))

# Performing the alignment of the samples by using kallisto aligner
rawAlignment(GSE74251csv,seq = "kallisto", downloadRef = TRUE,downloadSamples = TRUE, createIndex = TRUE, referenceGenome = 37, tx2Counts = tx2gene)
rawAlignment(GSE81593csv,seq = "kallisto", downloadRef = FALSE,downloadSamples = TRUE, createIndex = FALSE, referenceGenome = 37, tx2Counts = tx2gene)

# Creating the csv file with the information about the counts files location and the labels
Run <- GSE74251csv$Run
Path <- paste("ReferenceFiles/Samples/RNAseq/CountFiles/",GSE74251csv$Run,sep = "")
Class <- rep("Tumor", length(GSE74251csv$Run))
GSE74251CountsInfo <-  data.frame(Run = Run, Path = Path, Class = Class)

Run <- GSE81593csv$Run
Path <- paste("ReferenceFiles/Samples/RNAseq/CountFiles/",GSE81593csv$Run,sep = "")
Class <- rep("Control", length(GSE81593csv$Run))
GSE81593CountsInfo <-  data.frame(Run = Run, Path = Path, Class = Class)

mergedCountsInfo <- rbind(GSE74251CountsInfo, GSE81593CountsInfo)

write.csv(mergedCountsInfo, file = "ReferenceFiles/mergedCountsInfo.csv")

# Merging in one matrix all the count files indicated inside the CSV file
countsInformation <- countsToMatrix("ReferenceFiles/mergedCountsInfo.csv")

# Exporting to independent variables the counts matrix and the labels
countsMatrix <- countsInformation$countsMatrix
labels <- countsInformation$labels

# Downloading human annotation
myAnnotation <- getGenesAnnotation(rownames(countsMatrix),referenceGenome=37)

# Calculating gene expression values matrix using the counts matrix
expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation,genesNames = TRUE)

# Plotting the boxplot of the expression of each samples for all the genes
dataPlot(expressionMatrix,labels,mode = "boxplot", colours = c("blue", "red"),toPNG = TRUE, toPDF = TRUE)

# Performing the quality analysis of the samples
RNAseqQA(expressionMatrix)

# Applying sva batch effect removal method in order to try to correct unknown batch effect
svaMod <- batchEffectRemoval(expressionMatrix, as.factor(labels), method = "sva")

# Extracting DEGs that pass the imposed restrictions
DEGsInformation <- DEGsExtraction(expressionMatrixCorrected, as.factor(labels), lfc = 1.0, pvalue = 0.01,
                                       number = 100, svaCorrection = TRUE, svaMod = svaMod)
topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

# Plotting the expression of the first 12 DEGs for each of the samples in an ordered way
dataPlot(DEGsMatrix[1:12,],labels,mode = "orderedBoxplot", colours = c("blue", "red"),toPNG = TRUE, toPDF = TRUE)

# Plotting the expression of the first 12 DEGs separatelly for all the samples
dataPlot(DEGsMatrix[1:12,],labels,mode = "genesBoxplot",toPNG = TRUE, toPDF = TRUE)

# Plotting the heatmap of the first 12 DEGs separatelly for all the samples
dataPlot(DEGsMatrix[1:12,],labels,mode = "heatmap",toPNG = TRUE, toPDF = TRUE)

# Preparing matrix for Machine Learning process (ML)
DEGsMatrixML <- t(DEGsMatrix)

# Feature selection process with mRMR
mrmrRanking <- featureSelection(DEGsMatrixML,as.factor(labels),colnames(DEGsMatrixML), mode = "mrmr")

# CV function with k-NN
results_cv_knn <- knn_CV(DEGsMatrixML,as.factor(labels),names(mrmrRanking),10)

# Plotting the accuracy of all the folds evaluated in the CV process
dataPlot(results_cv_knn$accMatrix,mode = "classResults", main = "Accuracy for each fold", xlab = "Genes", ylab = "Accuracy",toPNG = TRUE, toPDF = TRUE)

# Plotting the sensitivity of all the folds evaluated in the CV process
dataPlot(results_cv_knn$sensMatrix,mode = "classResults", main = "Sensitivity for each fold", xlab = "Genes", ylab = "Sensitivity",toPNG = TRUE, toPDF = TRUE)

# Plotting the specificity of all the folds evaluated in the CV process
dataPlot(results_cv_knn$specMatrix,mode = "classResults", main = "Specificity for each fold", xlab = "Genes", ylab = "Specificity",toPNG = TRUE, toPDF = TRUE)

# Plotting the confusion matrix with the sum of the confusion matrices of each folds evaluated in the CV process
allCfMats <- results_cv_knn$cfMats[[1]]$table + results_cv_knn$cfMats[[2]]$table +
  results_cv_knn$cfMats[[3]]$table + results_cv_knn$cfMats[[4]]$table +
  results_cv_knn$cfMats[[5]]$table + results_cv_knn$cfMats[[6]]$table +
  results_cv_knn$cfMats[[7]]$table + results_cv_knn$cfMats[[8]]$table +
  results_cv_knn$cfMats[[9]]$table + results_cv_knn$cfMats[[10]]$table

dataPlot(allCfMats,labels,mode = "confusionMatrix",toPNG = TRUE, toPDF = TRUE)

# Test function with SVM
distribution <- sample(1:33, 33, replace=FALSE)
trainingDataset <- DEGsMatrixML[distribution[1:25],]
trainingLabels <- labels[distribution[1:25]]
testDataset <- DEGsMatrixML[distribution[26:33],]
testLabels <- labels[distribution[26:33]]

results_test_svm <- KnowSeq::svm_test(trainingDataset,as.factor(trainingLabels),testDataset,as.factor(testLabels),names(mrmrRanking))

# Plotting the accuracy for all the genes evaluated
dataPlot(results_test_svm$accVector,mode = "classResults", main = "Accuracy for all the genes evaluated", xlab = "Genes", ylab = "Accuracy",toPNG = TRUE, toPDF = TRUE)

# Plotting the sensitivity for all the genes evaluated
dataPlot(results_test_svm$sensVector,mode = "classResults", main = "Sensitivity for all the genes evaluated", xlab = "Genes", ylab = "Sensitivity",toPNG = TRUE, toPDF =TRUE)

# Plotting the specificity for all the genes evaluated
dataPlot(results_test_svm$specVector,mode = "classResults", main = "Specificity for all the genes evaluated", xlab = "Genes", ylab = "Specificity",toPNG = TRUE, toPDF = TRUE)

# Plotting the confusion matrix by using 100 genes to classify and assess the model
dataPlot(results_test_svm$cfMats[[100]]$table,testLabels,mode = "confusionMatrix",toPNG = TRUE, toPDF = TRUE)

# Retrieving the GO information from the three different ontologies
GOsInfo <- geneOntologyEnrichment(DEGsMatrix,labelsGo,nGOs = 20)

# Downloading and filling with expression the pathways of the DEGs
myDEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix)[1:10],
                                             referenceGenome=38,attributes = c("external_gene_name","entrezgene_id","gene_biotype")
                                             , filters = "external_gene_name")
allMyAnnotation <- getAnnotationFromEnsembl(rownames(expressionMatrix),
                                            referenceGenome=38,attributes = c("external_gene_name","entrezgene_id","gene_biotype")
                                            , filters = "external_gene_name")

DEGsPathwayVisualization(DEGsMatrix[1:10,], myDEGsAnnotation, expressionMatrix, allMyAnnotation, labels)

# Downloading the information about the DEGs related diseases
diseases <- DEGsToDiseases(rownames(DEGsMatrix))
