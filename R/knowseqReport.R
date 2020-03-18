#' knowseqReport creates a report for a given set of genes which their label.
#'
#' knowseqReport creates a report for a given set of genes which their label. This provide an html file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process).
#' @param data A matrix that contains the gene expression or counts values.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param outdir The output directory to store the report.
#' @param baseline A string that indicates the start point. This will be 'expression' if data contains genes expression values or 'counts' if data contains genes counts values.
#' @param featureSelectionMode String that indicates which feature selection algorithm is going to be used. Possible values are: mrmr, rf or da.
#' @param disease String that indicates from which disease wants the user to know if selected genes are related to. Found evidences will be shown. Default empty, this means that all related diseases, and found evidences, will be shown.
#' @param maxGenes Integer that indicates the maximun number of genes which information will be shown and that will be used to train models.
#' @param clasifAlgs A vector with including algorithms names that will be used in training cv.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' knowseqReport(expressionMatrix,labels,'knowSeq-report',clasifAlgs=c('rf'),disease='lung-cancer',maxGenes = 9)


knowseqReport <- function(data,labels,outdir="knowSeq-report",baseline='expression', lfc=2.0, pvalue=0.01, qualityAnalysis = TRUE, batchEffectTreatment =  TRUE,
                          geneOntology = TRUE, getPathways = TRUE, getDiseases = TRUE,
                          featureSelectionMode = 'mrmr', disease = '', maxGenes = Inf, clasifAlgs=c('knn','rf','svm'),
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
  
  table.format <- 'html'
  
  # Create output's directory if it doesn't exists
  if (! dir.exists(outdir)) dir.create(outdir)
  if (substr(outdir,nchar(outdir),nchar(outdir)) != '/') outdir <- paste(outdir,'/',sep='')
  act.folder <- getwd()
  # markobj contain the text that will be displayed in html or pdf file
  markobj <- c()
  
  markobj <- c(markobj,'<img src="https://github.com/CasedUgr/KnowSeq/blob/master/logoKnow.png?raw=true" style="position:absolute;top:0px;right:0px;" width=265px height=200px />\n')
  
  table.labels <- data.frame(table(labels))
  colnames(table.labels) <- c("Labels", "Freq.")
  
  markobj <- c(markobj,'# Summary',
               paste('This experiment is performed over ',dim(expressionMatrix)[2], ' samples, taking into account an 
                     inicial number of genes equal to ', dim(expressionMatrix)[1], '.\n','Furthermore, 
                     information about the different classes and samples per classes is shown herein: \n', sep = ""))
  
  markobj <- c(markobj,
               '```{r, echo=FALSE, fig.align="center"}',
               paste('knitr::kable(table.labels,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
               '```\n')
  
  
  
  if(qualityAnalysis){
    #cat("Performing the quality analysis of the samples\n")
    RNAseqQA(expressionMatrix)
  }
  
  # --- Differencia Expressed Genes --- #
  markobj <- c(markobj,'# Differential Expressed Genes extraction')
  
  if(batchEffectTreatment){
    
    markobj <- c(markobj,'## Treating Batch Effect\n',
                 'It is widely known that this is a crucial step in the omics data processing due to the intrinsic 
              deviations that the data can present due to its origin, sequencing design, etc...\n',
                 'Batch effect will be removed by using surrogate variable analysis or sva.\n')
    
    svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
    
    DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = lfc, pvalue = pvalue, number = Inf, svaCorrection = TRUE, svaMod = svaMod)
  }else{
    
    DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = lfc, pvalue = pvalue, number = Inf)
  }
  
  topTable <- DEGsInformation$Table
  
  if(dim(topTable)[1] == 0){
    stop("There is no any DEGs for this combination of LFC and P-value. Please, impose less restrictive thressholds.")
  }
  
  if(is.infinite(maxGenes)){
    maxGenes <- dim(topTable)[1]
  }
  
  if(dim(topTable)[1] < maxGenes){
    maxGenes <- dim(topTable)[1]
  }
  
  DEGsMatrix <- DEGsInformation$DEGsMatrix
  
  topTable.dataframe <- data.frame(GeneSymbol=rownames(topTable),logFC=topTable$logFC,AveExpr=topTable$AveExpr,t=topTable$t,
                                   P.Value=formatC(topTable$P.Value, format = "e", digits = 2),
                                   adj.P.Val=formatC(topTable$adj.P.Val, format = "e", digits = 2),B=topTable$B)
  
  colnames(topTable.dataframe) <- c("Gene Symbol","logFC","AveExpr", "t", "P-Value","adj. P-Value","B")
  
  markobj <- c(markobj,'## Searching for DEGs\n',
               paste('The search and extraction of Differential Expressed Genes is the main challenge for this type of analysis. In this 
                 sense, to achieve this extraction, a LFC greater or equal than ', lfc,' along with a P-Value lower or equal than ',pvalue,' are imposed.\n', sep = ""))
  
  markobj <- c(markobj,'Finally',paste(dim(DEGsMatrix)[1],'genes have been keeped after using DEGs extraction and can be seen in the table below: \n'))
  
  markobj <- c(markobj,
               '```{r, echo=FALSE, fig.align="center"}',
               paste('knitr::kable(topTable.dataframe,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
               '```\n')
  
  # --- Feature Selection --- #
  markobj <- c(markobj,'## Feature Selection',
               paste('With the purpose of finding the best DEGs order to assess the data, the',featureSelectionMode,'method
                     will be used in order to select the',maxGenes,'most relevant genes for the machine learning process.\n'))
  DEGsMatrixML <- t(DEGsMatrix)
  ranking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = featureSelectionMode)
  
  if (featureSelectionMode == 'mrmr') ranking <- names(sort(ranking,decreasing = FALSE))
  else if (featureSelectionMode == 'rf') ranking <- names(ranking)
  else if (featureSelectionMode == 'da') ranking <- names(ranking)
  
  genes <- ''
  for (gene in ranking[1:maxGenes]) genes <- paste(genes,gene,sep=', ')
  markobj <- c(markobj,paste('First',maxGenes,'selected genes are:'),sub(".","",genes),'.\n')
  
  markobj <- c(markobj,'## Visualization\n',
               'DEGs are genes that have a truly different expression among the studied classes, 
               for that it is important to try to see graphically if those DEGs comply with this requirement. 
               In order to provide a tool to perform this task, the function dataPlot encapsulate a set of 
               graphs that allows plotting in different ways the expression of the DEGs.\n')
  
  if(maxGenes > 12)
    boxplotGenes <- 12
  else
    boxplotGenes <- maxGenes
  
  markobj <- c(markobj,paste('\nIn the next boxplot, the expression of the first',maxGenes,'DEGs for each sample ordered by classes its showed.\n'),
               '```{r echo=FALSE}',
               paste("dataPlot(DEGsMatrix[ranking[1:",maxGenes,"],],labels,mode = 'orderedBoxplot',toPNG = FALSE,toPDF = FALSE)",sep=''),
               '```\n')
  
  markobj <- c(markobj,paste("However it is interesting to see the differentiation at gene expression level for each of the top",
                             boxplotGenes,"genes used before separately. This information is plotted below.\n"))
  
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               paste("dataPlot(DEGsMatrix[ranking[1:",boxplotGenes,"],],labels,mode = 'genesBoxplot',toPNG = FALSE,toPDF = FALSE)",sep=''),
               '```\n')
  
  markobj <- c(markobj,paste('Finally, the heatmap of those',maxGenes,'DEGs separatelly for all the samples is plotted below\n'))
  markobj <- c(markobj,
               '```{r echo=FALSE}',
               paste('dataPlot(DEGsMatrix[ranking[1:',maxGenes,'],],labels,mode = "heatmap",toPNG = FALSE,toPDF = FALSE)',sep=''),
               '```\n')
  
  
  # --- Machine learning --- #
  # --- ---  Training --- --- #
  markobj <- c(markobj,'# Machine Learning Assessment \n',
               'With the purpose of evaluate the robustness of the DEGs in the discernment among the studied pathologies.\n')
  clasifNames <- paste(clasifAlgs,collapse = ',')
  clasifNames <- str_replace(clasifNames,'knn','K-Nearest Neighbors (K-NN)')
  clasifNames <- str_replace(clasifNames,'rf','Random Forest')
  clasifNames <- str_replace(clasifNames,'svm','Support Vector Machine (SVM)')
  
  markobj <- c(markobj,paste('To this effect,',clasifNames,'classification algorithms will be trained using 5-Fold Cross Validation.
                    To evaluate obtained results the this metrics will be shown in the following plots:\n'))
  
  #if ('accuracy' %in% metrics) markobj <- c(markobj,)
  
  for (clasifAlg in clasifAlgs){
    if (clasifAlg == 'knn') results_cv <- knn_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
    else if (clasifAlg == 'rf') results_cv <- rf_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
    else if (clasifAlg == 'svm') results_cv <- svm_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
    
    markobj <- c(markobj,paste('## CV Results implementing ',clasifAlg),'\n')
    
    for (metric in metrics){
      if (metric == 'accuracy') act.metric = 'accMatrix'
      else if (metric == 'specificity') act.metric = 'specMatrix'
      else if (metric == 'sensitivity') act.metric = 'sensMatrix'
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(results_cv[["',act.metric,'"]],
                       mode = "classResults",
                       main = "',metric,' for each fold with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',metric,'")',sep=''),'```\n')
    }
    allCfMats <- results_cv$cfMats[[1]]$table + results_cv$cfMats[[2]]$table + results_cv$cfMats[[3]]$table + results_cv$cfMats[[4]]$table + results_cv$cfMats[[5]]$table
    markobj <- c(markobj,'```{r echo = FALSE}',
                 paste('dataPlot(allCfMats, labels,
                       mode = "confusionMatrix",
                       main = "Confusion Matrix for CV results using ',clasifAlg,'",
                       )',sep=''),'```\n')
    #dataPlot(allCfMats,labels,mode = "confusionMatrix")
  }
  # --- DEGs enrichment methodology --- #
  markobj <- c(markobj,'\n# DEGs enrichment\n',
               'The main goal of the this process is the extraction of biological relevant information from the DEGs,
               and this enrichment has three different points of view:\n
               \t- The gene ontology information.\n
               \t- The pathway visualization.\n 
               \t- The relationship between the DEGs and diseases related to the studied pathologies.\n')
  
  # --- Gene Ontology --- #
  
  if(geneOntology){
    markobj <- c(markobj,'## Gene Ontology\n',
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
    markobj <- c(markobj,'### BP Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(bp.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
    
    GOsMatrix$`MF Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`MF Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
    mf.frame <- data.frame(GOsMatrix$`MF Ontology GOs`)
    mf.frame <- mf.frame[,c("GO.ID","Term","Description","GO_Genes")]
    markobj <- c(markobj,'### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(mf.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
    
    GOsMatrix$`CC Ontology GOs`[,10] <- as.character(lapply(GOsMatrix$`CC Ontology GOs`[,10], function(x) {gsub(",", ", ", x)}))
    cc.frame <- data.frame(GOsMatrix$`CC Ontology GOs`)
    cc.frame <- cc.frame[,c("GO.ID","Term","Description","GO_Genes")]
    markobj <- c(markobj,'### CC Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(cc.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
  }
  
  # --- Pathways Visualization --- #
  if(getPathways){
    DEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix),notHSapiens=FALSE)
    genomeAnnotation <- getAnnotationFromEnsembl('allGenome',notHSapiens=FALSE)
    markobj <- c(markobj,'\n## Pathways visualization\n','```{r echo=FALSE}','DEGsPathwayVisualization(DEGsMatrix,DEGsAnnotation,expressionMatrix,genomeAnnotation)','```\n')
  }
  
  # --- Related Diseases --- #
  if(getDiseases){
    if (disease == ''){
      diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 5, getEvidences = TRUE)
      
      evidences.frame <- list()
      for (gene in names(diseases)){
        act.markobj <- c()
        # If user want to see evidences for all diseases or this diseases match with solicitated disease
        check.diseases <- names(diseases[[gene]]$evidences)

        for (act.disease in check.diseases){
          if ( class(diseases[[gene]]$evidences[[act.disease]]) == 'list' ){
            if (!gene %in% names(evidences.frame)) evidences.frame[[gene]] <- list()

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
                               paste('knitr::kable(data.frame(evidences.frame[["',gene,'"]][["',evidence.type,'"]]),"',table.format,'",caption="',evidence.type,' evidences for ',disease,'")',sep=''),'```')
          }
        }
        if (length(act.markobj) > 0)
          markobj <- c(markobj,paste('\n###',gene,sep=' '),act.markobj)
      }
    }else{
      markobj <- c(markobj,'## Related diseases\n',
                   'Finally, the related diseases enrichment is displayed. DEGs related diseases are searched 
                  from *targetValidation* plastform.\n')

      r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease,"&size=1&filter=disease",sep = ""))
      respon <- httr::content(r_Ensembl)
      
      if ( 'size' %in% names(respon) && respon$size == 0){
        markobj  <-  c(markobj,'Disease not found.')
      }else{
        disease.id <- respon$data[[1]]$id
        url  <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=",disease.id,"&size=10000",sep='')
        response <- httr::GET(url)
        response <- httr::content(response)
        found.symbols <- unlist(list.map(r$data,target$gene_info$symbol))
        
        found.symbols <- intersect(found.symbols,rownames(DEGsMatrix))

        evidences <- DEGsEvidences(found.symbols,disease)
        
        evidences.frame <- list()
        for (gene in names(evidences)){
          act.markobj <- c()
          if (class(evidences[[gene]])=='list'){
            evidences.frame[[gene]] = list()
            for ( evidence.type in names(evidences[[gene]]) ){
              act.evidences.frame <- c()
              for (act.evidence in evidences[[gene]][[evidence.type]]){
                act.evidences.frame <- rbind(act.evidences.frame,act.evidence$evidence)
              }
              remove.cols <- c()
              for (col in seq(dim(act.evidences.frame)[2])){
                if (all(act.evidences.frame[,col]=='*')  || all(act.evidences.frame[,col]=='')){
                  remove.cols <- c(remove.cols,col)
                }
              }
              if (length(remove.cols)>0){
                act.evidences.frame <-  act.evidences.frame[,-remove.cols]
              }
              evidences.frame[[gene]][[evidence.type]] <- act.evidences.frame
            }
            for ( evidence.type in names(evidences.frame[[gene]]))
              act.markobj <- c(act.markobj,'```{r echo=FALSE}',
                               paste('knitr::kable(data.frame(evidences.frame[["',gene,'"]][["',evidence.type,'"]]),"',table.format,'",caption="',evidence.type,' evidences for ',disease,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```')
          }
          if (length(act.markobj) > 0)
            markobj <- c(markobj,paste('\n###',gene,sep=' '),act.markobj)
        }
      }
    }
  }
  
  # --- Save Report --- #
  mark.header.html <- c('---',
                        'subtitle: Powered by KnowSeq R/Bioc package from University of Granada',
                        'title: "Gene Expression Intelligent Pipeline Report"',
                        paste('date: "Date of the Report: ', Sys.time(),'"',sep = ""),
                        'output:',
                        ' html_document:',
                        '  number_sections: yes',
                        '  toc: yes',
                        '---')
  
  write(file = "report.Rmd", c(mark.header.html,markobj))
  
  dir <- system.file("extdata", package="KnowSeq")

  rmarkdown::render(input = "report.Rmd", output_file = paste(outdir,'report.html',sep='/'),output_format = rmarkdown::html_document(
    theme = "default",
    mathjax = NULL,
    highlight = NULL,
    toc = TRUE,
    toc_depth = 3,
    toc_float = TRUE,
    number_sections = TRUE,
    self_contained = TRUE,
    css = paste(dir,"/report_style.css",sep = "")
  ))
  file.remove("report.Rmd")
  browseURL(paste(outdir,'report.html',sep='/'))
}