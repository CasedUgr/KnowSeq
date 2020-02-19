#' getReport creates a report for a given set of genes which their label
#'
#' getReport creates a report for a given set of genes which their label. This provide an html and/or pdf file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process)
#' @param data A matrix that contains the gene expression or counts values.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param outdir The output directory to store the report
#' @param baseline A string that indicates the start point. This will be 'expression' if data contains genes expression values or 'counts' if data contains genes counts values.
#' @param outputFormat String to indicate in which formar will be the report saved. There are two ossible values: pdf or hmtl.
#' @param featureSelectionMode String that indicates which feature selection algorithm is going to be used. Possible values are: mrmr, rf or da.
#' @param maxGenes Integer that indicates the maximun number of genes which information will be shown and that will be used to train models
#' @param clasifAlgs A vector with including algorithms names that will be used in training cv.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' getReport(expressionMatrix,labels,'pdf-report',clasifAlgs=c('rf'))



getReport <- function(data,labels,outdir,baseline='expression',
                      outputFormat='pdf',
                      featureSelectionMode = 'mrmr',
                      maxGenes = 12,
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
               '### Remove Batch Effect\n')
  
  #DEGsInformation <- limmaDEGsExtraction(t(expressionMatrixCorrected), labels, lfc = 1.0, pvalue = 0.01, number = 100)
  svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
  DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = 2.0, pvalue = 0.01, number = Inf, svaCorrection = TRUE, svaMod = svaMod)
  topTable <- DEGsInformation$Table
  DEGsMatrix <- DEGsInformation$DEGsMatrix
  
  markobj <- c(markobj,paste(dim(expressionMatrix)[1] - dim(DEGsMatrix)[1],'genes have been removed using *sva* method. \n'))
  
  
  # --- Feature Selection --- #
  DEGsMatrixML <- t(DEGsMatrix)
  ranking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = featureSelectionMode)
  
  if (featureSelectionMode == 'mrmr') ranking <- names(sort(ranking,decreasing = FALSE))
  else if (featureSelectionMode == 'da') ranking <- names(ranking)
  
  genes <- ''
  for (gene in ranking[1:12]) genes <- paste(genes,gene,sep=', ')
  markobj <- c(markobj,'### Feature Selection\n','First 12 selected genes are:',sub(".","",genes))
  
  markobj <- c(markobj,'\nThe expression of those 12 selected DEGs for each of the samples in an ordered way is plotted below\n',
               '```{r echo=FALSE}',
               "dataPlot(DEGsMatrix[ranking[1:12],],labels,mode = 'orderedBoxplot',toPNG = FALSE,toPDF = FALSE)",
               '```\n')
  
  markobj <- c(markobj,"The expression of those 12 DEGs separatelly for all the samples is plotted below.\n")
  
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               "dataPlot(DEGsMatrix[ranking[1:12],],labels,mode = 'genesBoxplot',toPNG = FALSE,toPDF = FALSE)",
               '```\n')
  
  markobj <- c(markobj,'The heatmap of those 12 DEGs separatelly for all the samples is plotted below\n')
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               'dataPlot(DEGsMatrix[ranking[1:12],],labels,mode = "heatmap",toPNG = FALSE,toPDF = FALSE)',
               '```\n')
  
  
  # --- Machine learning --- #
  # --- ---  Training --- --- #
  markobj <- c(markobj,'### Training \n')
  
  for (clasifAlg in clasifAlgs){
    if (clasifAlg == 'knn') results_cv <- knn_CV(DEGsMatrixML,labels,ranking[1:12],5)
    else if (clasifAlg == 'rf') results_cv <- rf_CV(DEGsMatrixML,labels,ranking[1:12],5)
    else if (clasifAlg == 'svm') results_cv <- svm_CV(DEGsMatrixML,labels,ranking[1:12],5)
    
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
  markobj <- c(markobj,'\n## DEGs enrichment\n',
               '### Gene Ontology\n')
  
  labelsGo <- labels
  for (i in seq(length(unique(labels)))){
    labelsGo <- gsub(unique(labels)[i],i-1,labelsGo) 
  }
  
  GOsMatrix <- geneOntologyEnrichment(DEGsMatrix,labelsGo,nGOs = 20)
  GOsMatrix$`BP Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`BP Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  bp.frame <- data.frame(GOsMatrix$`BP Ontology GOs`)
  bp.frame <- bp.frame[,c("GO.ID","GO_Genes","Term","Description")]
  markobj <- c(markobj,'#### BP Ontology GOs\n','```{r}',paste('knitr::kable(bp.frame)',sep=''),'```\n')
  
  #GOsMatrix$`MF Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`MF Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  #markobj <- c(markobj,'#### MF Ontology GOs\n','```{r echo=FALSE}','knitr::kable(data.frame(GOsMatrix$`MF Ontology GOs`))','```\n')
  #GOsMatrix$`CC Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`CC Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  #markobj <- c(markobj,'#### CC Ontology GOs\n','```{r echo=FALSE}','knitr::kable(data.frame(GOsMatrix$`CC Ontology GOs`),longtable = T)','```\n')
  
  # --- Pathways Visualization --- #
  #DEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix),notHSapiens=FALSE)
  #genomeAnnotation <- getAnnotationFromEnsembl('allGenome',notHSapiens=FALSE)
  
  #markobj <- c(markobj,'\n### Pathways visualization\n','```{r echo=FALSE}','DEGsPathwayVisualization(DEGsMatrix,DEGsAnnotation,expressionMatrix,genomeAnnotation)','```\n')
  
  # --- Related Diseases --- #
  markobj <- c(markobj,'### Related diseases')
  diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 5, getEvidences = TRUE)
  
  for (gene in names(diseases)){
    markobj <- c(markobj,paste('\n####',gene,sep=' '))

    for (act.disease in names(diseases[[gene]]$evidences)){
      if ( class(diseases[[gene]]$evidences[[act.disease]]) == 'list' ){
        markobj <- c(markobj,paste('\n#####',act.disease,sep=' '))
        for ( evidence.type in names(diseases[[gene]]$evidences[[act.disease]]) ){
          act.evidences.frame <- c()
          for (act.evidence in diseases[[gene]]$evidences[[act.disease]][[evidence.type]]){
            act.evidences.frame <- rbind(act.evidences.frame,act.evidence$evidence)
          }
          
          markobj <- c(markobj,
                       '```{r echo=FALSE}',paste('knitr::kable(data.frame(act.evidences.frame),caption="',evidence.type,' evidences for ',act.disease,'")',sep=''),'```')
        }
      }
    }
  }
  
  # --- Save Report --- #
  if ( outputFormat == 'html' ){
    mark.header.html <- c('---',
                          'title: "Genes Report"',
                          'output: html_document',
                          '---',
                          '')
    markdown::markdownToHTML(text = knitr::knit(text = c(mark.header.html,markobj)), output =paste(outdir,'report.html',sep='/'))
    browseURL('report.html')
  }
  else{
    mark.header.pdf <- c('---',
                         'title: "Genes Report"',
                         'output: pdf_document',
                         '---',
                         '')
    file.create(paste(outdir,'report.rmd',sep='/'))
    fileConn<-file(paste(outdir,'report.rmd',sep='/'))
    writeLines(c(mark.header.pdf,markobj), fileConn)
    close(fileConn)
    render(paste(outdir,'report.rmd',sep='/'))
    file.remove(paste(outdir,'report.rmd',sep='/'))
    file.remove(paste(outdir,'report.log',sep='/'))
  }
}


