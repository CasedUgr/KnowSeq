#' getReport creates a report for a given set of genes which their label
#'
#' getReport creates a report for a given set of genes which their label. This provide an html and/or pdf file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process)
#' @param data A matrix that contains the gene expression or counts values.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param outdir The output directory to store the report
#' @param baseline A string that indicates the start point. This will be 'expression' if data contains genes expression values or 'counts' if data contains genes counts values.
#' @param outputFormat String to indicate in which formar will be the report saved. There are two ossible values: pdf or hmtl.
#' @param featureSelectionMode String that indicates which feature selection algorithm is going to be used. Possible values are: mrmr, rf or da.
#' @param disease String that indicates from which disease wants the user to know if selected genes are related to. Found evidences will be shown. Default empty, this means that all related diseases, and found evidences, will be shown.
#' @param maxGenes Integer that indicates the maximun number of genes which information will be shown and that will be used to train models
#' @param clasifAlgs A vector with including algorithms names that will be used in training cv.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' getReport(expressionMatrix,labels,'pdf-report',outputFormat='pdf',clasifAlgs=c('rf'),disease='cancer',maxGenes = 9)


getReport <- function(data,labels,outdir,baseline='expression',
                      outputFormat='pdf',
                      featureSelectionMode = 'mrmr',
                      disease = '',
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
  if ( outputFormat == 'html') table.format <- 'html'
  else table.format = 'pandoc'
  
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
               '### Remove Batch Effect\n',
              'It is widely known that this is a crucial step in the omics data processing due to the intrinsic 
              deviations that the data can present due to its origin, sequencing design, etc...\n',
              'In this section, batch effect will be removed by using surrogate variable analysis or sva.\n',
              'This method does not return a matrix with the batch effect corrected, instead of this, 
              the function returns a model that has to be used as single input parameter of the function limmaDEGsExtraction.\n')
  
  #DEGsInformation <- limmaDEGsExtraction(t(expressionMatrixCorrected), labels, lfc = 1.0, pvalue = 0.01, number = 100)
  svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
  DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = 2.0, pvalue = 0.01, number = Inf, svaCorrection = TRUE, svaMod = svaMod)
  topTable <- DEGsInformation$Table
  DEGsMatrix <- DEGsInformation$DEGsMatrix
  
  markobj <- c(markobj,'### Differential Expressed Genes extraction using *Limma*',
               'limmaDEGsExtraction performs the analysis to extract the Differentially Expressed Genes (DEGs) 
               among the classes to compare. This method is used setting LFC param to 2 and p-value to 0.01.\n',
               'Finnally',paste(dim(DEGsMatrix)[1],'genes have been keeped after using limma DEGs extraction with *sva* method. \n'))
    
  
  # --- Feature Selection --- #
  markobj <- c(markobj,'### Feature Selection',
               paste('With the purpose of ﬁnding the best DEGs order to assess the data, the',featureSelectionMode,'method
                     will be used in order to select the',maxGenes,'most relevant genes for the machine learning process.\n'))
  DEGsMatrixML <- t(DEGsMatrix)
  ranking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = featureSelectionMode)
  
  if (featureSelectionMode == 'mrmr') ranking <- names(sort(ranking,decreasing = FALSE))
  else if (featureSelectionMode == 'da') ranking <- names(ranking)
  
  genes <- ''
  for (gene in ranking[1:maxGenes]) genes <- paste(genes,gene,sep=', ')
  markobj <- c(markobj,paste('First',maxGenes,'selected genes are:'),sub(".","",genes),'.\n')
  
  markobj <- c(markobj,'### Visualization\n',
                'DEGs are genes that have a truly diﬀerent expression among the studied classes, 
               for that it is important to try to see graphically if those DEGs comply with this requirement. 
               In order to provide a tool to perform this task, the function dataPlot encapsulate a set of 
               graphs that allows plotting in different ways the expression of the DEGs.\n')

  markobj <- c(markobj,paste('\nIn the next boxplot, the expression of the first',maxGenes,'DEGs for each sample its showed.\n'),
               '```{r echo=FALSE}',
               paste("dataPlot(DEGsMatrix[ranking[1:",maxGenes,"],],labels,mode = 'orderedBoxplot',toPNG = FALSE,toPDF = FALSE)",sep=''),
               '```\n')

  markobj <- c(markobj,paste("However it is interesting to see the differentiation at gene expression level for each of the top",
                             maxGenes,"genes used before separately. This information is plotted below.\n"))
  
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               paste("dataPlot(DEGsMatrix[ranking[1:",maxGenes,"],],labels,mode = 'genesBoxplot',toPNG = FALSE,toPDF = FALSE)",sep=''),
               '```\n')
  
  markobj <- c(markobj,paste('Finally, the heatmap of those',maxGenes,'DEGs separatelly for all the samples is plotted below\n'))
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               paste('dataPlot(DEGsMatrix[ranking[1:',maxGenes,'],],labels,mode = "heatmap",toPNG = FALSE,toPDF = FALSE)',sep=''),
               '```\n')
  
  
  # --- Machine learning --- #
  # --- ---  Training --- --- #
  markobj <- c(markobj,'### Training \n',
               'With the purpose of evaluate the robustness of the DEGs in the discernment among the studied pathologies.\n')
  clasifNames <- paste(clasifAlgs,collapse = ',')
  clasifNames <- str_replace(clasifNames,'knn','K-Nearest Neighbors (K-NN)')
  clasifNames <- str_replace(clasifNames,'rf','Random Forest')
  clasifNames <- str_replace(clasifNames,'svm','Support Vector Machine (SVM)')
  
  markobj <- c(markobj,paste('To this effect,',clasifNames,'classiﬁcation algorithms will be trained using 5-Fold Cross Validation.
                    To evaluate obtained results the this metrics will be shown in the following plots:\n'))
  
  #if ('accuracy' %in% metrics) markobj <- c(markobj,)
  
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
  markobj <- c(markobj,'\n## DEGs enrichment\n',
               'The main goal of the this process is the extraction of biological relevant information from the DEGs,
               and this enrichment has three different points of view:\n
               \t- The gene ontology information.\n
               \t- The pathway visualization.\n 
               \t- The relationship between the DEGs and diseases related to the studied pathologies.\n')
  
  # --- Gene Ontology --- #
  markobj <- c(markobj,'### Gene Ontology\n',
               'Gene ontology (GO) provides information about the biological functions of the genes. 
              The following, information from the three different ontologies (BP, MF and CC) will be shown.\n')
  labelsGo <- labels
  for (i in seq(length(unique(labels)))){
    labelsGo <- gsub(unique(labels)[i],i-1,labelsGo) 
  }
  
  GOsMatrix <- geneOntologyEnrichment(DEGsMatrix,labelsGo,nGOs = 20)
  GOsMatrix$`BP Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`BP Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  bp.frame <- data.frame(GOsMatrix$`BP Ontology GOs`)
  bp.frame <- bp.frame[,c("GO.ID","Term","Description","GO_Genes")]
  markobj <- c(markobj,'#### BP Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(bp.frame,"',table.format,'")',sep=''),'```\n')
  
  GOsMatrix$`MF Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`MF Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  mf.frame <- data.frame(GOsMatrix$`MF Ontology GOs`)
  mf.frame <- mf.frame[,c("GO.ID","Term","Description","GO_Genes")]
  markobj <- c(markobj,'#### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(mf.frame,"',table.format,'")',sep=''),'```\n')

  GOsMatrix$`CC Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`CC Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
  cc.frame <- data.frame(GOsMatrix$`CC Ontology GOs`)
  cc.frame <- cc.frame[,c("GO.ID","Term","Description","GO_Genes")]
  markobj <- c(markobj,'#### CC Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(cc.frame,"',table.format,'")',sep=''),'```\n')
  
  # --- Pathways Visualization --- #
  #DEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix),notHSapiens=FALSE)
  #genomeAnnotation <- getAnnotationFromEnsembl('allGenome',notHSapiens=FALSE)
  #markobj <- c(markobj,'\n### Pathways visualization\n','```{r echo=FALSE}','DEGsPathwayVisualization(DEGsMatrix,DEGsAnnotation,expressionMatrix,genomeAnnotation)','```\n')
  
  # --- Related Diseases --- #
  markobj <- c(markobj,'### Related diseases\n',
            'Finally, the related diseases enrichment is displayed. DEGs related diseases are searched 
            from *targetValidation* plastform.\n')

  diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 5, getEvidences = TRUE)
  
  evidences.frame <- list()
  for (gene in names(diseases)){
    act.markobj <- c()
    # If user want to see evidences for all diseases or this diseases match with solicitated disease
    if (disease == ''){ check.diseases <- names(diseases[[gene]]$evidences)
    }else if (!is.na(grep(disease,diseases[[gene]]$summary[,1])[1])) {
      check.diseases <- diseases[[1]]$summary[,1][grep(disease,diseases[[1]]$summary[,1])[1]]
    }else check.diseases <- c()
    
    for (act.disease in check.diseases){
      if ( class(diseases[[gene]]$evidences[[act.disease]]) == 'list' ){
        evidences.frame[[gene]] <- list()
        act.markobj <- c(act.markobj,paste('Found evidences of',gene,'related with',act.disease,sep=' '))
        for ( evidence.type in names(diseases[[gene]]$evidences[[act.disease]]) ){
          act.evidences.frame <- c()
          for (act.evidence in diseases[[gene]]$evidences[[act.disease]][[evidence.type]]){
            act.evidences.frame <- rbind(act.evidences.frame,act.evidence$evidence)
          }
          evidences.frame[[gene]][[evidence.type]] <- act.evidences.frame
        }
        for ( evidence.type in names(evidences.frame[[gene]]))
          act.markobj <- c(act.markobj,'```{r echo=FALSE}',
                       paste('knitr::kable(data.frame(evidences.frame[["',gene,'"]][["',evidence.type,'"]]),"',table.format,'",caption="',evidence.type,' evidences for ',act.disease,'")',sep=''),'```')
      }
    }
    if (length(act.markobj) > 0)
      markobj <- c(markobj,paste('\n####',gene,sep=' '),act.markobj)
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
                         '---','')#, '```{r}','library(knitr)','```')
    file.create(paste(outdir,'report.rmd',sep='/'))
    fileConn<-file(paste(outdir,'report.rmd',sep='/'))
    writeLines(c(mark.header.pdf,markobj), fileConn)
    close(fileConn)
    render(paste(outdir,'report.rmd',sep='/'))
    file.remove(paste(outdir,'report.rmd',sep='/'))
    file.remove(paste(outdir,'report.tex',sep='/'))
    file.remove(paste(outdir,'report.log',sep='/'))
  }
}
