#' getReport creates a report for a given set of genes which their label
#'
#' getReport creates a report for a given set of genes which their label. This provide an html and/or pdf file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process)
#' @param data A matrix that contains the gene expression or counts values.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param outdir The output directory to store the report
#' @param baseline A string that indicates the start point. This will be 'expression' if data contains genes expression values or 'counts' if data contains genes counts values.
#' @param clasifAlgs A vector with including algorithms names that will be used in training cv.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' getReport(expressionMatrix,labels,'pdf-report',clasifAlgs=c('rf'))


getReport <- function(data,labels,outdir,baseline='expression',
                      toHTML = TRUE, toPDF = TRUE,
                      featureSelectionAlg = 'mrmr',
                      clasifAlgs=c('knn','rf','svm'),
                      metrics=c('accuracy','specificity','sensitivity')){
  # --- Check params --- #
  if(baseline == 'expression'){
    if(!is.data.frame(data) && !is.matrix(data)){
      
      stop("The data argument must be a dataframe or a matrix.")
      
    }
    if (dim(data)[2] != length(labels)){
      
      stop("The length of the columns of the argument data must be the same than the length of the lables. Please, ensures that the rows are the samples and the columns are the variables.")
      
    }
    expressionMatrix <- data
  }
  else if (baseline == 'count'){

    myAnnotation <- getAnnotationFromEnsembl(rownames(data),notHSapiens = FALSE)
    expressionMatrix <- calculateGeneExpressionValues(data,myAnnotation, genesNames = TRUE)

  }
  # Create output's directory if it doesn't exists
  if (! dir.exists(outdir)) dir.create(outdir)
  if (substr(outdir,nchar(outdir),nchar(outdir)) != '/') outdir <- paste(outdir,'/',sep='')
  act.folder <- getwd()

  # markobj contain the text that will be displayed in html or pdf file
  markobj <- c()

  #cat("Performing the quality analysis of the samples\n")
  #RNAseqQA(expressionMatrix)
  
  # --- Differencia Expressed Genes --- #
  markobj <- c(markobj,'## Differential Expressed Genes extraction and visualization',
               '### Extracting DEGs\n')

  #DEGsInformation <- limmaDEGsExtraction(t(expressionMatrixCorrected), labels, lfc = 1.0, pvalue = 0.01, number = 100)
  svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
  DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = 2.0, pvalue = 0.01, number = Inf, svaCorrection = TRUE, svaMod = svaMod)
  topTable <- DEGsInformation$Table 
  DEGsMatrix <- DEGsInformation$DEGsMatrix
  
  markobj <- c(markobj,'\nPlotting the expression of the first 12 DEGs for each of the samples in an ordered way\n',
               '```{r echo=FALSE}',
               "dataPlot(DEGsMatrix[1:12,],labels,mode = 'orderedBoxplot',toPNG = FALSE,toPDF = FALSE)",
               '```\n')

  markobj <- c(markobj,"Plotting the expression of the first 12 DEGs separatelly for all the samples.\n")
  
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               "dataPlot(DEGsMatrix[1:12,],labels,mode = 'genesBoxplot',toPNG = FALSE,toPDF = FALSE)",
               '```\n')
  
  markobj <- c(markobj,'Plotting the heatmap of the first 12 DEGs separatelly for all the samples\n')
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               'dataPlot(DEGsMatrix[1:12,],labels,mode = "heatmap",toPNG = FALSE,toPDF = FALSE)',
               '```\n')
  
  
  # --- Machine learning --- #
  # --- ---  Feature Selection --- --- #
  DEGsMatrixML <- t(DEGsMatrix)
  mrmrRanking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = "mrmr")
  mrmrRanking <- names(sort(mrmrRanking,decreasing = FALSE))

  markobj <- c(markobj,'### MRMR Feature Selection',
               '```{r echo=FALSE}',
               'mrmrRanking[1:10]',
               '```\n')
  
  # --- ---  Training --- --- #
  markobj <- c(markobj,'### Training \n')

  for (clasifAlg in clasifAlgs){
    if (clasifAlg == 'knn') results_cv <- knn_CV(DEGsMatrixML,labels,mrmrRanking[1:10],5)
    else if (clasifAlg == 'rf') results_cv <- rf_CV(DEGsMatrixML,labels,mrmrRanking[1:10],5)
    else if (clasifAlg == 'svm') results_cv <- svm_CV(DEGsMatrixML,labels,mrmrRanking[1:10],5)
    
    markobj <- c(markobj,paste('####',clasifAlg),'\n')
    
    for (metric in metrics){
      if (metric == 'accuracy') act.metric = 'accMatrix'
      else if (metric == 'specificity') act.metric = 'specMatrix'
      else if (metric == 'sensitivity') act.metric = 'sensMatrix'
      markobj <- c(markobj,'```{r echo = FALSE}',
                  paste('dataPlot(results_cv[["',act.metric,'"]],
                       mode = "classResults",
                       main = "',metric,' for each fold with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',metric,'")',sep=''),'```')
    }
    #allCfMats <- results_cv$cfMats[[1]]$table + results_cv$cfMats[[2]]$table + results_cv$cfMats[[3]]$table + results_cv$cfMats[[4]]$table + results_cv$cfMats[[5]]$table
    #dataPlot(allCfMats,labels,mode = "confusionMatrix")
  }
  # --- DEGs enrichment methodology --- #
  # --- Gene Ontology --- #
  markobj <- c(markobj,'## DEGs enrichment',
               '### Gene Ontology')
  
  labelsGo <- labels
  for (i in seq(length(unique(labels)))){
    labelsGo <- gsub(unique(labels)[i],i-1,labelsGo) 
  }
  
  GOsMatrix <- geneOntologyEnrichment(DEGsMatrix,labelsGo,nGOs = 20)
  GOsMatrix$`BP Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`BP Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  markobj <- c(markobj,'#### BP Ontology GOs\n','```{r}','knitr::kable(data.frame(GOsMatrix$`BP Ontology GOs`),longtable = T)','```\n')
  GOsMatrix$`MF Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`MF Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  markobj <- c(markobj,'#### MF Ontology GOs\n','```{r}','knitr::kable(data.frame(GOsMatrix$`MF Ontology GOs`),longtable = T)','```\n')
  GOsMatrix$`CC Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`CC Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  markobj <- c(markobj,'#### CC Ontology GOs\n','```{r}','knitr::kable(data.frame(GOsMatrix$`CC Ontology GOs`),longtable = T)','```\n')
  
  # --- Pathways Visualization --- #
  DEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix),notHSapiens=FALSE)
  genomeAnnotation <- getAnnotationFromEnsembl('allGenome',notHSapiens=FALSE)
  
  markobj <- c(markobj,'\n### Pathways visualization\n',       
               '```{r echo=FALSE}','DEGsPathwayVisualization(DEGsMatrix,DEGsAnnotation,expressionMatrix,genomeAnnotation)',
               '```\n')
  
  # --- Related Diseases --- #
  markobj <- c(markobj,'### Related diseases')
  diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 5)
  
  for (gene in names(diseases)){
    markobj <- c(markobj,paste('\n\t- **',gene,'**.',sep=''))
    for (disease in diseases[[gene]]$summary[,1])
      markobj <- c(markobj,disease,', ')
    markobj <- markobj[-length(markobj)]
  }
  
  # --- Save Report --- #
  if ( toHTML ){
    mark.header.html <- c('---',
                          'title: "Genes Report"',
                          'output: html_document',
                          '---',
                          '')
    setwd(outdir)
    markdown::markdownToHTML(text = knitr::knit(text = c(mark.header.html,markobj)), output ='report.html')
    browseURL('report.html')
  }
  if ( toPDF ){
    mark.header.pdf <- c('---',
                         'title: "Genes Report"',
                         'output: pdf_document',
                         '---',
                         '')
    
    file.create('report.rmd')
    fileConn<-file('report.rmd')
    writeLines(c(mark.header.pdf,markobj), fileConn)
    close(fileConn)
    render('report.rmd')
    file.remove('report.rmd')
    file.remove('report.log')
  }
  setwd(act.folder)
}


