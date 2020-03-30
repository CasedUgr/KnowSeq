#' knowseqReport creates a report for a given set of genes which their label.
#'
#' knowseqReport creates a report for a given set of genes which their label. This provide an html file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process).
#' @param data A matrix that contains the gene expression.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param MLTest This parameter enables the classification process for a test dataset. 
#' @param testData A matrix that contains the unseen samples for the test process.
#' @param testLabels A vector or factor that contains the labels for the unseen samples for the test process.
#' @param outdir The output directory to store the report.
#' @param qualityAnalysis A logical parameter that indicates if the user wants to perform the quality anaylisis or not.
#' @param batchEffectTreatment A logical parameter that indicates if the user wants to perform the batch effect treatment or not.
#' @param geneOntology A logical parameter that indicates if the user wants to show genes ontologies or not.
#' @param getPathways A logical parameter that indicates if the user wants to show genes pathways or not.
#' @param getDiseases A logical parameter that indicates if the user wants to show genes related diseases or not.
#' @param pvalue The value of the p-value which determines the DEGs. If one or more genes have a p-value lower or equal to the selected p-value, they would be considered as DEGs.
#' @param lfc The value of the LFC which determines the DEGs. If one or more genes have a LFC greater or equal to the selected LFC, they would be considered as DEGs.
#' @param cov This value only works when there are more than two classes in the labels. This parameter stablishs a minimum number of pair of classes combination in which exists differential expression to consider a genes as expressed genes.
#' @param featureSelectionMode String that indicates which feature selection algorithm is going to be used. Possible values are: mrmr, rf or da. By default, no feature selection algorithm will be applied.
#' @param disease String that indicates from which disease wants the user wants to know if selected genes are related to. Found evidences will be shown for each subdiseases. Default empty, this means that all related diseases, and found evidences, will be shown.
#' @param subdiseases String that indicates the name of a particular subtype from disease, which  the  user to know if selected genes are related to. Found evidences will be shown. Default empty, this means that there are not subtypes of disease to look for, all found evidences for disease will be shown.
#' @param maxGenes Integer that indicates the maximun number of genes which information will be shown and that will be used to train models.
#' @param clasifAlgs A vector with including algorithms names that will be used in training cv.
#' @param metrics A list with metrics that the user wants to be shown in machine learning process. Metrics could be accuracy, specificity and/or sensitivity.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' \dontrun{knowseqReport(expressionMatrix,labels,'knowSeq-report',clasifAlgs=c('rf'),disease='lung-cancer',maxGenes = 9)}
#' \dontrun{knowseqReport(expressionMatrix,labels,'knowSeq-report',clasifAlgs=c('rf'),disease='lung-cancer',subdiseases=c('squamous cell lung carcinoma','lung adenocarcinoma'),maxGenes = 9)}


knowseqReport <- function(data, labels, MLTest = FALSE, testData="", testLabels="",outdir="knowSeq-report", qualityAnalysis = TRUE, batchEffectTreatment =  TRUE,
                          geneOntology = TRUE, getPathways = TRUE, getDiseases = TRUE,
                          lfc=2.0, pvalue=0.01, cov=2, 
                          featureSelectionMode = 'nofs', disease = '',subdiseases=c(''), maxGenes = Inf, clasifAlgs=c('knn','rf','svm'),
                          metrics=c('accuracy','specificity','sensitivity')){
  
  removeEmptyColumns <- function(data){
    remove.cols <- c()
    for (col in seq(dim(data)[2])){
      if (all(data[,col]=='*')  || all(data[,col]=='')){
        remove.cols <- c(remove.cols,col)
      }
    }
    if (length(remove.cols)>0){
      col.names <- colnames(data)
      data <-  data[,-remove.cols]
      col.names <- col.names[-remove.cols]
      
      if (is(data) != 'matrix')
        data <- as.matrix(data)
      if (dim(data)[2] != length(col.names))
        data <- t(data)
      colnames(data) <- col.names
    }
    return(data)
  }
  
  # --- Check params --- #
  if(!is.data.frame(data) && !is.matrix(data)){
    
    stop("The data argument must be a dataframe or a matrix.")
    
  }
  if (dim(data)[2] != length(labels)){
    
    stop("The length of the columns of the argument data must be the same than the length of the labels. Please, ensures that the rows are the samples and the columns are the variables.")
    
  }
  
  if (MLTest == TRUE){
    
    if(!is.data.frame(testData) && !is.matrix(testData)){
      
      stop("The data for test must be a dataframe or a matrix.")
      
    }
    if (dim(testData)[2] != length(testLabels)){
      
      stop("The length of the columns of the argument data for test must be the same than the length of the labels for test. Please, ensures that the rows are the samples and the columns are the variables.")
      
    }
  }
  
  expressionMatrix <- data
  if (geneOntology || getPathways)
    myAnnotation <- getGenesAnnotation(rownames(expressionMatrix),attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"),filter="external_gene_name",notHSapiens = FALSE)
    

  disease <- gsub(' ','-',disease)
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
               paste('This experiment is performed over ',dim(expressionMatrix)[2], ' with an initial number of ', dim(expressionMatrix)[1], ' genes.\n','A summary 
                     of the classes labels along with the samples per class can be observed below: \n', sep = ""))
  
  markobj <- c(markobj,
               '```{r, echo=FALSE, fig.align="center"}',
               paste('knitr::kable(table.labels,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
               '```\n')
  
  
  
  if(qualityAnalysis){
    cat("Performing the quality analysis of the samples\n")
    RNAseqQA(expressionMatrix,outdir=paste(outdir,'RNAseqQA',sep=''))
    markobj <- c(markobj,'# Quality analysis\n',
                 'A quality analysis must be performed in order to detect and remove any possible outlier contained within the samples. ',
                 'Outliers are numerically different samples when compared with the rest of the samples, and this numerical difference might affect study results. ', 
                 'For that purpose, arrayQualityMetrics bioc package is used to carry out different statistical tests, detecting possible outliers. ',
                 'arrayQualityMetrics generates a report containing all this information that can be seen [**HERE**](RNAseqQA/index.html).\n')
  }
  
  # --- Differencia Expressed Genes --- #
  markobj <- c(markobj,'# Differential Expressed Genes extraction\n')
  
  if(batchEffectTreatment){
    
    markobj <- c(markobj,'## Treating Batch Effect\n',
                 'Batch effect is produced due to intrinsic deviations inside the data due to its origin, sequencing design, lab, technician, etc... Therefore, it is crucial to perform a correct treatment of 
                 it in this type of analysis.','Taking into account that the different Batches are unknown, the effect will be treated by using surrogate variable analysis or sva algorithm.\n')
    
    svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
    
    DEGsInformation <- DEGsExtraction(expressionMatrix, labels, lfc = lfc, pvalue = pvalue, cov = cov, number = Inf, svaCorrection = TRUE, svaMod = svaMod)
  }else{
    
    DEGsInformation <- DEGsExtraction(expressionMatrix, labels, lfc = lfc, pvalue = pvalue, cov = cov, number = Inf)
  }
  
  topTable <- DEGsInformation$Table
  
  if(dim(topTable)[1] == 0){
    stop("There is no any DEGs for this combination of LFC and P-value. Please, impose less restrictive thressholds.")
  }
  
  
  if(length(levels(as.factor(labels))) == 2){
    
    DEGsMatrix <- DEGsInformation$DEGsMatrix
    
    topTable.dataframe <- data.frame(GeneSymbol=rownames(topTable),logFC=topTable$logFC,AveExpr=topTable$AveExpr,t=topTable$t,
                                     P.Value=formatC(topTable$P.Value, format = "e", digits = 2),
                                     adj.P.Val=formatC(topTable$adj.P.Val, format = "e", digits = 2),B=topTable$B)
    
    colnames(topTable.dataframe) <- c("Gene Symbol","logFC","AveExpr", "t", "P-Value","adj. P-Value","B")
    
    markobj <- c(markobj,'## Searching for DEGs\n',
                 paste('The search and extraction of Differential Expressed Genes is the main challenge for this type of study. Thus, to achieve a set 
                       of possible biomarkers the following thresholds will be selected: LFC greater or equal than ', lfc,', P-Value lower or equal than ',pvalue,'.\n', sep = ""))
    
      markobj <- c(markobj,'Finally',paste(dim(DEGsMatrix)[1],'DEGs have been kept after using DEGs extraction and they can be seen in the table below: \n'))
    
    markobj <- c(markobj,
                 '```{r, echo=FALSE, fig.align="center"}',
                 paste('knitr::kable(topTable.dataframe,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
                 '```\n')
    
  }else if(length(levels(as.factor(labels))) > 2){
    
    DEGsMatrix <- DEGsInformation$DEGsMatrix
    genes.selected <- rownames(DEGsInformation$MulticlassLFC)
    
    topTable.dataframe <- data.frame(GeneSymbol=genes.selected,logFC=rowMeans(abs(DEGsInformation$MulticlassLFC)),
                                     t=rowMeans(topTable[genes.selected,]$t),
                                     P.Value=formatC(rowMeans(topTable[genes.selected,]$p.value), format = "e", digits = 2),
                                     F=topTable[genes.selected,]$F,B=rowMeans(topTable[genes.selected,]$lods))
    colnames(topTable.dataframe) <- c("Gene Symbol","logFC ($\\mu$)","t ($\\mu$)", "P-Value ($\\mu$)","F","B ($\\mu$)")
    topTable.dataframe <- topTable.dataframe[order(topTable.dataframe[,2],decreasing=TRUE),]
    rownames(topTable.dataframe) <- NULL
    
    markobj <- c(markobj,'## Searching for Multiclass DEGs\n',
                 paste('The search and extraction of Differential Expressed Genes is the main challenge for this type of study. Thus, to achieve a 
                 set of possible biomarkers the following thresholds will be selected: LFC greater or equal than ', lfc,', P-Value lower or equal than ',pvalue,'. Furthermore, for multiclass 
                       assessment a coverage equal or greater than ', cov, ' will be used.\n', sep = ""))
    
    markobj <- c(markobj,'Finally',paste(dim(DEGsMatrix)[1],'multiclass DEGs have been kept after using DEGs extraction and they can be seen in the table below: \n'))
    
    markobj <- c(markobj,
                 '```{r, echo=FALSE, fig.align="center"}',
                 paste('knitr::kable(topTable.dataframe,"',table.format,'", table.attr = "class=\'paleBlueRows\'", escape=FALSE)',sep=''),
                 '```\n')
    
  }
  
  if(geneOntology || getPathways){
    myAnnotation <- myAnnotation[myAnnotation$external_gene_name %in% rownames(DEGsMatrix),]
    myAnnotation <- myAnnotation[complete.cases(myAnnotation), ]
  }
  
  if(is.infinite(maxGenes)){
    maxGenes <- dim(DEGsMatrix)[1]
  }
  
  if(dim(DEGsMatrix)[1] < maxGenes){
    maxGenes <- dim(DEGsMatrix)[1]
  }
  
  
  # --- Feature Selection --- #
  DEGsMatrixML <- t(DEGsMatrix)
  
  if(featureSelectionMode != 'nofs'){
    
    markobj <- c(markobj,'## Feature Selection',
                 paste('With the purpose of finding the best combination of DEGs to assess the data, the',featureSelectionMode,'method
                     will be used in order to select the',maxGenes,'most relevant genes for the classification process.\n'))
    
    ranking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = featureSelectionMode)
    if (featureSelectionMode == 'mrmr') ranking <- names(sort(ranking,decreasing = FALSE))
    else if (featureSelectionMode == 'rf') ranking <- names(ranking)
    else if (featureSelectionMode == 'da') ranking <- names(ranking)
    
    markobj <- c(markobj,paste('First',maxGenes,'selected genes by ',featureSelectionMode,' algorithm/method are:'),ranking[seq_len(maxGenes)],'.\n')
    
    
  } else{ranking <- rownames(DEGsMatrix)}
  
  
  genes <- ''
  for (gene in ranking[seq_len(maxGenes)]) genes <- paste(genes,gene,sep=', ')
  
  markobj <- c(markobj,'## Visualization\n',
               'DEGs are genes that have a truly different expression among the studied classes, 
               so it is important to see graphically if those DEGs comply with this requirement. 
               In order to provide a tool to perform this task, the function dataPlot encapsulates a set of 
               graphs that allows plotting in different ways the expression of the DEGs.\n')
  
  if(maxGenes > 12)
    boxplotGenes <- 12
  else
    boxplotGenes <- maxGenes
  
  markobj <- c(markobj,paste('\nIn the next boxplot, the expression of the first',maxGenes,'DEGs for each sample ordered by classes its showed.\n'),
               '```{r echo=FALSE}',
               paste("dataPlot(DEGsMatrix[ranking[1:",maxGenes,"],],labels,mode = 'orderedBoxplot',toPNG = FALSE,toPDF = FALSE)",sep=''),
               '```\n')
  
  markobj <- c(markobj,paste("However it is of interest to observe the differentiation at gene expression level for each of the top",
                             boxplotGenes,"genes previously used. This information is plotted below.\n"))
  
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
  markobj <- c(markobj,'# Machine Learning Assessment \n')
  clasifNames <- paste(clasifAlgs,collapse = ',')
  clasifNames <- str_replace(clasifNames,'knn','K-Nearest Neighbors (K-NN)')
  clasifNames <- str_replace(clasifNames,'rf','Random Forest')
  clasifNames <- str_replace(clasifNames,'svm','Support Vector Machine (SVM)')
  
  if(MLTest == FALSE){
  markobj <- c(markobj,paste('With the purpose of evaluating the robustness of the DEGs for the classification task between the studied pathologies, 
  a supervised classification step will be performed. 
  For it,',clasifNames,'classification algorithms will be trained using 10-Fold Cross Validation.
                    To evaluate obtained results the Mean Accuracy,Mean Specificity, Mean Sensitivity and the Confusion Matrix will be shown in the following plots:\n'))
  }else{
    markobj <- c(markobj,paste('With the purpose of evaluating the robustness of the DEGs for the classification task between the studied pathologies, 
  a supervised classification step will be performed. 
  For it,',clasifNames,'classification algorithms will be trained using 10-Fold Cross Validation and then, tested in a separated classification step using unseen samples.
                    To evaluate obtained results the Mean Accuracy,Mean Specificity, Mean Sensitivity and the Confusion Matrix will be shown in the following plots:\n'))
  }
  
  for (clasifAlg in clasifAlgs){
    if (clasifAlg == 'knn'){ 
      results_cv_knn <- knn_CV(DEGsMatrixML,labels,ranking[seq_len(maxGenes)],10)
      markobj <- c(markobj,paste('## Results for 10-CV implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy'){
          act.metric = 'accMatrix'
          colour = "red"
        }else if (metric == 'specificity'){
          act.metric = 'specMatrix'
          colour = "blue"
        }else if (metric == 'sensitivity'){
          act.metric = 'sensMatrix'
          colour = "green"
        }
        
        s <- strsplit(metric, " ")[[1]]
        s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                   sep = "", collapse = " ")
        
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(colMeans(results_cv_knn[["',act.metric,'"]]),
                       mode = "classResults",
                       main = "Mean ',s,' results for 10-CV with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
        
      }
      
      allCfMats_knn <- results_cv_knn$cfMats[[1]]$table + results_cv_knn$cfMats[[2]]$table + results_cv_knn$cfMats[[3]]$table + results_cv_knn$cfMats[[4]]$table + results_cv_knn$cfMats[[5]]$table + results_cv_knn$cfMats[[6]]$table + results_cv_knn$cfMats[[7]]$table + results_cv_knn$cfMats[[8]]$table + results_cv_knn$cfMats[[9]]$table + results_cv_knn$cfMats[[10]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_knn, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
      
      if(MLTest == TRUE){
        
        results_test_knn <- knn_test(DEGsMatrixML,labels,t(testData),testLabels,ranking[seq_len(maxGenes)],bestK = results_cv_knn$bestK)
        markobj <- c(markobj,paste('## Test Results implementing ',clasifAlg),'\n')
        
        for (metric in metrics){
          if (metric == 'accuracy'){
            test.metric = 'accVector'
            colour = "red"
          }else if (metric == 'specificity'){
            test.metric = 'specVector'
            colour = "blue"
          }else if (metric == 'sensitivity'){
            test.metric = 'sensVector'
            colour = "green"
          }
          
          s <- strsplit(metric, " ")[[1]]
          s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                     sep = "", collapse = " ")
          
          markobj <- c(markobj,'```{r echo = FALSE}',
                       paste('dataPlot(results_test_knn[["',test.metric,'"]],
                       mode = "classResults",
                       main = "',s,' test results with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
        }
      
      
        testConfMatrixknn <- results_test_knn$cfMats[[maxGenes]]$table
        markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(testConfMatrixknn, testLabels, mode = "confusionMatrix")',sep=''),'```\n')
      }
      
    }else if (clasifAlg == 'rf'){
      results_cv_rf <- rf_CV(DEGsMatrixML,labels,ranking[seq_len(maxGenes)],10)
      markobj <- c(markobj,paste('## Results for 10-CV implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy'){
          act.metric = 'accMatrix'
          colour = "red"
        }else if (metric == 'specificity'){
          act.metric = 'specMatrix'
          colour = "blue"
        }else if (metric == 'sensitivity'){
          act.metric = 'sensMatrix'
          colour = "green"
        }
        
        s <- strsplit(metric, " ")[[1]]
        s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                   sep = "", collapse = " ")
        
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(colMeans(results_cv_rf[["',act.metric,'"]]),
                       mode = "classResults",
                       main = "Mean ',s,' results for 10-CV with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
      }
      
      allCfMats_rf <- results_cv_rf$cfMats[[1]]$table + results_cv_rf$cfMats[[2]]$table + results_cv_rf$cfMats[[3]]$table + results_cv_rf$cfMats[[4]]$table + results_cv_rf$cfMats[[5]]$table + results_cv_rf$cfMats[[6]]$table + results_cv_rf$cfMats[[7]]$table + results_cv_rf$cfMats[[8]]$table + results_cv_rf$cfMats[[9]]$table + results_cv_rf$cfMats[[10]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_rf, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
      
      if(MLTest == TRUE){
        
        results_test_rf <- rf_test(DEGsMatrixML,labels,t(testData),testLabels,ranking[seq_len(maxGenes)])
        markobj <- c(markobj,paste('## Test Results implementing ',clasifAlg),'\n')
        
        for (metric in metrics){
          if (metric == 'accuracy'){
            test.metric = 'accVector'
            colour = "red"
          }else if (metric == 'specificity'){
            test.metric = 'specVector'
            colour = "blue"
          }else if (metric == 'sensitivity'){
            test.metric = 'sensVector'
            colour = "green"
          }
          
          s <- strsplit(metric, " ")[[1]]
          s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                     sep = "", collapse = " ")
          
          markobj <- c(markobj,'```{r echo = FALSE}',
                       paste('dataPlot(results_test_rf[["',test.metric,'"]],
                       mode = "classResults",
                       main = "',s,' test results with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
        }
        
        
        testConfMatrixrf <- results_test_rf$cfMats[[maxGenes]]$table
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(testConfMatrixrf, testLabels, mode = "confusionMatrix")',sep=''),'```\n')
      }
      
    }else if (clasifAlg == 'svm'){
      results_cv_svm <- svm_CV(DEGsMatrixML,labels,ranking[seq_len(maxGenes)],10)
      markobj <- c(markobj,paste('## Results for 10-CV implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy'){
          act.metric = 'accMatrix'
          colour = "red"
        }else if (metric == 'specificity'){
          act.metric = 'specMatrix'
          colour = "blue"
        }else if (metric == 'sensitivity'){
          act.metric = 'sensMatrix'
          colour = "green"
        }
        
        s <- strsplit(metric, " ")[[1]]
        s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
              sep = "", collapse = " ")
        
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(colMeans(results_cv_svm[["',act.metric,'"]]),
                       mode = "classResults",
                       main = "Mean ',s,' results for 10-CV with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
      }
      
      allCfMats_svm <- results_cv_svm$cfMats[[1]]$table + results_cv_svm$cfMats[[2]]$table + results_cv_svm$cfMats[[3]]$table + results_cv_svm$cfMats[[4]]$table + results_cv_svm$cfMats[[5]]$table + results_cv_svm$cfMats[[6]]$table + results_cv_svm$cfMats[[7]]$table + results_cv_svm$cfMats[[8]]$table + results_cv_svm$cfMats[[9]]$table + results_cv_svm$cfMats[[10]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_svm, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
      
      if(MLTest == TRUE){
        
        results_test_svm <- svm_test(DEGsMatrixML,labels,t(testData),testLabels,ranking[seq_len(maxGenes)], bestParameters = results_cv_svm$bestParameters)
        markobj <- c(markobj,paste('## Test Results implementing ',clasifAlg),'\n')
        
        for (metric in metrics){
          if (metric == 'accuracy'){
            test.metric = 'accVector'
            colour = "red"
          }else if (metric == 'specificity'){
            test.metric = 'specVector'
            colour = "blue"
          }else if (metric == 'sensitivity'){
            test.metric = 'sensVector'
            colour = "green"
          }
          
          s <- strsplit(metric, " ")[[1]]
          s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                     sep = "", collapse = " ")
          
          markobj <- c(markobj,'```{r echo = FALSE}',
                       paste('dataPlot(results_test_svm[["',test.metric,'"]],
                       mode = "classResults",
                       main = "',s,' test results with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',s,'", colours = "',colour,'", xgrid=TRUE, ygrid=TRUE)',sep=''),'```\n')
        }
        
        
        testConfMatrixsvm <- results_test_svm$cfMats[[maxGenes]]$table
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(testConfMatrixsvm, testLabels, mode = "confusionMatrix")',sep=''),'```\n')
      }
    }
    
  }
  
  if(geneOntology | getPathways | getDiseases){
    # --- DEGs enrichment methodology --- #
    markobj <- c(markobj,'\n# DEGs enrichment\n',
                 'The main goal of the this process is the extraction of biological relevant information from the DEGs.
                   The enrichment process has three different approaches:\n
                   \t- The gene ontology information.\n
                   \t- The pathway visualization.\n 
                   \t- The relationship between the DEGs and diseases related to the studied pathologies.\n')
    
    # --- Gene Ontology --- #
    
    if(geneOntology){
      markobj <- c(markobj,'## Gene Ontology\n',
                   'Gene ontology (GO) provides information about the biological functions of the genes. 
                      Information from the three different ontologies (BP, MF and CC) will be shown.\n')
      amigo.url <- 'http://amigo.geneontology.org/amigo/term/'
      GOsMatrix <- geneOntologyEnrichment(as.character(myAnnotation$entrezgene_id),geneType='ENTREZ_GENE_ID',pvalCutOff=0.1,returnGeneSymbols = FALSE)
      
      if(dim(GOsMatrix$`BP Ontology GOs`)[1] != 0){
        gene.names <- as.character(GOsMatrix$`BP Ontology GOs`$Genes)
        for ( gen in seq(dim(myAnnotation)[1])){
          gene.names <- str_replace(gene.names,as.character(myAnnotation[gen,'entrezgene_id']),
                                    as.character(myAnnotation[gen,'external_gene_name']))
        }

        GOsMatrix$`BP Ontology GOs`['Gene Symbols'] <- gene.names
        GOsMatrix$`BP Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`BP Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
        bp.frame <- GOsMatrix$`BP Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
        rownames(bp.frame) <- NULL
        bp.frame$GO.ID  <- paste('[',bp.frame$GO.ID,'](',amigo.url,bp.frame$GO.ID,')',sep='')
        
        markobj <- c(markobj,'### BP Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(bp.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
        
      }else{markobj <- c(markobj,'### BP Ontology GOs\n There are no GO terms related to this set of DEGs for BP ontology\n')}
      
      if(dim(GOsMatrix$`MF Ontology GOs`)[1] != 0){
        gene.names <- GOsMatrix$`MF Ontology GOs`$Genes
        for ( gen in seq(dim(myAnnotation)[1])){
          gene.names <- str_replace(gene.names,as.character(myAnnotation[gen,'entrezgene_id']),
                                    as.character(myAnnotation[gen,'external_gene_name']))
        }
        GOsMatrix$`MF Ontology GOs`['Gene Symbols'] <- gene.names
        
        GOsMatrix$`MF Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`MF Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
        mf.frame <- GOsMatrix$`MF Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
        rownames(mf.frame) <- NULL
        mf.frame$GO.ID  <- paste('[',mf.frame$GO.ID,'](',amigo.url,mf.frame$GO.ID,')',sep='')
        
        markobj <- c(markobj,'### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(mf.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
        
      }else{markobj <- c(markobj,'### MF Ontology GOs\n There are no GO terms related to this set of DEGs for MF ontology\n')}
      
      if(dim(GOsMatrix$`CC Ontology GOs`)[1] != 0){
        gene.names <- GOsMatrix$`CC Ontology GOs`$Genes
        for ( gen in seq(dim(myAnnotation)[1])){
          gene.names <- str_replace(gene.names,as.character(myAnnotation[gen,'entrezgene_id']),
                                    as.character(myAnnotation[gen,'external_gene_name']))
        }
        GOsMatrix$`CC Ontology GOs`['Gene Symbols'] <- gene.names
        
        GOsMatrix$`CC Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`CC Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
        cc.frame <- GOsMatrix$`CC Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
        rownames(cc.frame) <- NULL
        cc.frame$GO.ID  <- paste('[',cc.frame$GO.ID,'](',amigo.url,cc.frame$GO.ID,')',sep='')
        markobj <- c(markobj,'### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(cc.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
      }else{markobj <- c(markobj,'### CC Ontology GOs\n There are no GO terms related to this set of DEGs for CC ontology\n')}
      
    }
    
    # --- Pathways Visualization --- #
    if(getPathways){
      
      myAnnotation <- myAnnotation[which(!is.na(myAnnotation$entrezgene_id) == TRUE),]
      
      commonDEGs <- intersect(rownames(DEGsMatrix),unique(myAnnotation$external_gene_name))
      posCommonDEGs <- match(rownames(DEGsMatrix[commonDEGs,]),myAnnotation$external_gene_name)
      pathMatrix <- DEGsMatrix[commonDEGs,]
      rownames(pathMatrix) <- myAnnotation$entrezgene_id[posCommonDEGs]
      
      paths.data <- matrix(ncol = 3)
      
      cat("Retrieving DEGs associated pathways...\n")
      
      for(gene.index in seq(dim(myAnnotation)[1])){
        gene <- myAnnotation[gene.index,'entrezgene_id']
        gene.name <- myAnnotation[gene.index,'external_gene_name']
        if(!is.na(gene)){
          get_GO <- GET(paste("http://rest.kegg.jp/get/hsa:",gene,sep = ""))
          get_GO_text <- content(get_GO, "text")
          pathway_start <- str_locate_all(pattern = "PATHWAY", get_GO_text)[[1]][2]
          pathway_end <- str_locate_all(pattern = "BRITE", get_GO_text)[[1]][1]
          pathways <- substr(get_GO_text,pathway_start+1,pathway_end-1)
          pathways <- strsplit(pathways,split = "\n")
          
          index <- grep("hsa", unlist(pathways))
          
          for(i in index){
            pathway <- as.character(unlist(str_extract_all(unlist(pathways)[i],"hsa[a-zA-Z0-9]{5}")))
            if (length(pathway) == 0) break
            start <- str_locate_all(pattern = "hsa", unlist(pathways)[i])[[1]][2]
            name <- substr(unlist(pathways)[i],start+8,nchar(unlist(pathways)[i]))
            
            paths.data <- rbind(paths.data,c(as.character(pathway),as.character(name),as.character(gene.name)))
          }
          
        }
      }
      paths.data <- paths.data[-1,]
      paths.data <- as.data.frame(paths.data,stringsAsFactors=FALSE)
      names(paths.data) <- c("KEGG_hsa","Name","Genes")
      naPos <- which(is.na(paths.data) == TRUE)
      
      # Collapse repeated pathways
      remove.index <- c()
      repeated.pathways <- names(which(table(paths.data$KEGG_hsa) > 1))
      for (pathway in repeated.pathways){
        index <- which(paths.data$KEGG_hsa == pathway,)
        genes <- paths.data[index,'Genes']
        paths.data[index[1],'Genes'] <- paste(genes,collapse=', ')
        remove.index <- c(remove.index,index[-1])
      }
      paths.data <- paths.data[-remove.index,]
      rownames(paths.data) <- NULL
      
      paths.data$KEGG_hsa <- paste('[',paths.data$KEGG_hsa,'](https://www.genome.jp/dbget-bin/www_bget?pathway:',paths.data$KEGG_hsa,')',sep='')
      markobj <- c(markobj,'\n## Pathways Extraction\n',
                   'In this step the pathways in which inserted genes appear are shown.\n')
      
      markobj <- c(markobj,'```{r echo=FALSE}',paste('knitr::kable(paths.data,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
      
    }
    
    # --- Related Diseases --- #
    if(getDiseases){
      if (disease == ''){
        
        markobj <- c(markobj,'## DEGs Related diseases\n',
                     'Finally, the related diseases enrichment is displayed. Related diseases to the final set of DEGs are searched 
                      from *targetValidation* plastform. For each disease it is also showed a list of evidences. \n')
        
        diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 10, getEvidences = TRUE)
        evidences.frame <- list()
        for (gene in names(diseases)){
          act.markobj <- c()
          # If user want to see evidences for all diseases or this diseases match with solicitated disease
          check.diseases <- names(diseases[[gene]]$evidences)
          
          for (act.disease in check.diseases){
            if ( is(diseases[[gene]]$evidences[[act.disease]]) == 'list' ){
              if (!gene %in% names(evidences.frame)) evidences.frame[[gene]] <- list()
              
              act.markobj <- c(act.markobj,paste('####',act.disease,sep=' '))
              for ( evidence.type in names(diseases[[gene]]$evidences[[act.disease]]) ){
                act.evidences.frame <- c()
                for (act.evidence in diseases[[gene]]$evidences[[act.disease]][[evidence.type]]){
                  act.evidences.frame <- rbind(act.evidences.frame,act.evidence$evidence)
                }
                # Remove empty columns
                act.evidences.frame <- removeEmptyColumns(act.evidences.frame)
                evidences.frame[[gene]][[evidence.type]] <- act.evidences.frame
              }
              for ( evidence.type in names(evidences.frame[[gene]]))
                act.markobj <- c(act.markobj,'```{r echo=FALSE}',
                                 paste('knitr::kable(data.frame(evidences.frame[["',gene,'"]][["',evidence.type,'"]]),"',table.format,'",caption="',evidence.type,' evidences for ',disease,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```')
            }
          }
          if (length(act.markobj) > 0)
            markobj <- c(markobj,paste('\n###',gene,sep=' '),act.markobj)
        }
      }else{
        
        dis <- gsub("-"," ",disease)
        dis <- strsplit(dis, " ")[[1]]
        dis <- paste(toupper(substring(dis, 1, 1)), substring(dis, 2),
                   sep = "", collapse = " ")
        
        markobj <- c(markobj,paste('## DEGs Evidences for ',dis,'\n',sep = ""),
                     'Finally, the ',dis,' related evidences enrichment is displayed. DEGs related evidences are searched 
                      from *targetValidation* platform.\n')
        
        r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease,"&size=1&filter=disease",sep = ""))
        respon <- content(r_Ensembl)
        
        if ( 'size' %in% names(respon) && respon$size == 0){
          markobj  <-  c(markobj,'\nDisease not found.')
        }else{
          disease.id <- respon$data[[1]]$id
          url  <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=",disease.id,"&size=10000",sep='')
          response <- GET(url)
          response <- content(response)
          found.symbols <- unlist(list.map(response$data,target$gene_info$symbol))
          found.symbols <- intersect(found.symbols,rownames(DEGsMatrix))
          
          if(length(found.symbols) > 0){
            evidences_ <- c()
            for (subdisease in subdiseases){
              evidences_ <- rbind(evidences_,DEGsEvidences(found.symbols,disease,subdisease))
            }
            evidences.frame <- list()
            for (gene in found.symbols){
              act.markobj <- c()
              for (ev.index in seq(dim(evidences_)[1])){
                if (!as.character(ev.index)  %in% names(evidences.frame))
                  evidences.frame[[as.character(ev.index)]]<-list()
                evidences <- evidences_[ev.index,]
                if (is(evidences[[gene]])=='list'){
                  evidences.frame[[as.character(ev.index)]][[gene]] = list()
                  for ( evidence.type in names(evidences[[gene]]) ){
                    act.evidences.frame <- c()
                    for (act.evidence in evidences[[gene]][[evidence.type]]){
                      act.evidences.frame <- rbind(act.evidences.frame,act.evidence$evidence)
                    }
                    # Remove empty columns
                    act.evidences.frame <- removeEmptyColumns(act.evidences.frame)
                    evidences.frame[[as.character(ev.index)]][[gene]][[evidence.type]] <- act.evidences.frame
                  }
                  if (subdisease != '' && length(evidences.frame[[as.character(ev.index)]])>0)
                    act.markobj <- c(act.markobj,paste('####',subdiseases[ev.index]))
                  for ( evidence.type in names(evidences.frame[[as.character(ev.index)]][[gene]]))
                    act.markobj <- c(act.markobj,'```{r echo=FALSE}',
                                     paste('knitr::kable(data.frame(evidences.frame[["',as.character(ev.index),'"]][["',gene,'"]][["',evidence.type,'"]]),"',table.format,'",caption="',evidence.type,' evidences for ',disease,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```')
                }
              }
              if (length(act.markobj) > 0)
                markobj <- c(markobj,paste('\n###',gene,sep=' '),act.markobj)
            }
          }
          else{
            markobj <- c(markobj,paste('\nNo introduced gene is related with',disease,sep=' '))
          }
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

  render(input = "report.Rmd", output_file = paste(outdir,'report.html',sep='/'),output_format = rmarkdown::html_document(
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

}


