#' knowseqReport creates a report for a given set of genes which their label.
#'
#' knowseqReport creates a report for a given set of genes which their label. This provide an html file with all the information that can be obtained for a certain set of genes (as GO, pathway visualization, associated diseases) and their labels (machine learning process).
#' @param data A matrix that contains the gene expression or counts values.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param outdir The output directory to store the report.
#' @param baseline A string that indicates the start point. This will be 'expression' if data contains genes expression values or 'counts' if data contains genes counts values.
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
#' knowseqReport(expressionMatrix,labels,'knowSeq-report',clasifAlgs=c('rf'),disease='lung-cancer',maxGenes = 9)
#' knowseqReport(expressionMatrix,labels,'knowSeq-report',clasifAlgs=c('rf'),disease='lung-cancer',subdiseases=c('squamous cell lung carcinoma','lung adenocarcinoma'),maxGenes = 9)


knowseqReport <- function(data,labels,outdir="knowSeq-report",baseline='expression', qualityAnalysis = TRUE, batchEffectTreatment =  TRUE,
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
      
      if (class(data) != 'matrix')
        data <- as.matrix(data)
      if (dim(data)[2] != length(col.names))
        data <- t(data)
      colnames(data) <- col.names
    }
    return(data)
  }
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
                 paste('The search and extraction of Differential Expressed Genes is the main challenge for this type of analysis. In this 
                   sense, to achieve this extraction, a LFC greater or equal than ', lfc,' along with a P-Value lower or equal than ',pvalue,' are imposed.\n', sep = ""))
    
    markobj <- c(markobj,'Finally',paste(dim(DEGsMatrix)[1],'genes have been keeped after using DEGs extraction and can be seen in the table below: \n'))
    
    markobj <- c(markobj,
                 '```{r, echo=FALSE, fig.align="center"}',
                 paste('knitr::kable(topTable.dataframe,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
                 '```\n')
    
  }else if(length(levels(as.factor(labels))) > 2){
    
    DEGsMatrix <- DEGsInformation$DEGsMatrix
    
    markobj <- c(markobj,'## Searching for Multiclass DEGs\n',
                 paste('The search and extraction of Differential Expressed Genes is the main challenge for this type of analysis. In this 
                   sense, to achieve this extraction, a LFC greater or equal than ', lfc,' along with a P-Value lower or equal than ',pvalue,' are imposed. Furthermore, for multiclass 
                       assessment a coverage equal or greater than ', cov, ' will be used\n', sep = ""))
    
    markobj <- c(markobj,'Finally',paste(dim(DEGsMatrix)[1],'genes have been keeped after using DEGs extraction and can be seen in the table below: \n'))
    
    markobj <- c(markobj,
                 '```{r, echo=FALSE, fig.align="center"}',
                 paste('knitr::kable(DEGsInformation$MulticlassIndex,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),
                 '```\n')
    
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
                 paste('With the purpose of finding the best DEGs order to assess the data, the',featureSelectionMode,'method
                     will be used in order to select the',maxGenes,'most relevant genes for the machine learning process.\n'))
    
    ranking <- featureSelection(DEGsMatrixML,labels,colnames(DEGsMatrixML), mode = featureSelectionMode)
    if (featureSelectionMode == 'mrmr') ranking <- names(sort(ranking,decreasing = FALSE))
    else if (featureSelectionMode == 'rf') ranking <- names(ranking)
    else if (featureSelectionMode == 'da') ranking <- names(ranking)
    
  } else{ranking <- rownames(DEGsMatrix)}
  
  
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
  
  for (clasifAlg in clasifAlgs){
    if (clasifAlg == 'knn'){ 
      results_cv_knn <- knn_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
      markobj <- c(markobj,paste('## CV Results implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy') act.metric = 'accMatrix'
        else if (metric == 'specificity') act.metric = 'specMatrix'
        else if (metric == 'sensitivity') act.metric = 'sensMatrix'
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(results_cv_knn[["',act.metric,'"]],
                       mode = "classResults",
                       main = "',metric,' for each fold with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',metric,'")',sep=''),'```\n')
      }
      allCfMats_knn <- results_cv_knn$cfMats[[1]]$table + results_cv_knn$cfMats[[2]]$table + results_cv_knn$cfMats[[3]]$table + results_cv_knn$cfMats[[4]]$table + results_cv_knn$cfMats[[5]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_knn, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
      
    }else if (clasifAlg == 'rf'){
      results_cv_rf <- rf_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
      markobj <- c(markobj,paste('## CV Results implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy') act.metric = 'accMatrix'
        else if (metric == 'specificity') act.metric = 'specMatrix'
        else if (metric == 'sensitivity') act.metric = 'sensMatrix'
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(results_cv_rf[["',act.metric,'"]],
                       mode = "classResults",
                       main = "',metric,' for each fold with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',metric,'")',sep=''),'```\n')
      }
      allCfMats_rf <- results_cv_rf$cfMats[[1]]$table + results_cv_rf$cfMats[[2]]$table + results_cv_rf$cfMats[[3]]$table + results_cv_rf$cfMats[[4]]$table + results_cv_rf$cfMats[[5]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_rf, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
      
    }else if (clasifAlg == 'svm'){
      results_cv_svm <- svm_CV(DEGsMatrixML,labels,ranking[1:maxGenes],5)
      markobj <- c(markobj,paste('## CV Results implementing ',clasifAlg),'\n')
      
      for (metric in metrics){
        if (metric == 'accuracy') act.metric = 'accMatrix'
        else if (metric == 'specificity') act.metric = 'specMatrix'
        else if (metric == 'sensitivity') act.metric = 'sensMatrix'
        markobj <- c(markobj,'```{r echo = FALSE}',
                     paste('dataPlot(results_cv_svm[["',act.metric,'"]],
                       mode = "classResults",
                       main = "',metric,' for each fold with ',clasifAlg,'",
                       xlab = "Genes", ylab ="',metric,'")',sep=''),'```\n')
      }
      allCfMats_svm <- results_cv_svm$cfMats[[1]]$table + results_cv_svm$cfMats[[2]]$table + results_cv_svm$cfMats[[3]]$table + results_cv_svm$cfMats[[4]]$table + results_cv_svm$cfMats[[5]]$table
      markobj <- c(markobj,'```{r echo = FALSE}',
                   paste('dataPlot(allCfMats_svm, labels,
                       mode = "confusionMatrix")',sep=''),'```\n')
    }

  }
  
  if(geneOntology | getPathways | getDiseases){
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
        amigo.url <- 'http://amigo.geneontology.org/amigo/term/'
        data <- getAnnotationFromEnsembl(rownames(DEGsMatrix),attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"),filter='external_gene_name')
        GOsMatrix <- geneOntologyEnrichment(as.character(data$entrezgene_id),geneType='ENTREZ_GENE_ID',pvalCutOff=0.1,returnGeneSymbols = TRUE)
        
        if(dim(GOsMatrix$`BP Ontology GOs`)[1] != 0){
          GOsMatrix$`BP Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`BP Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
          bp.frame <- GOsMatrix$`BP Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
          rownames(bp.frame) <- NULL
          bp.frame$GO.ID  <- paste('[',bp.frame$GO.ID,'](',amigo.url,bp.frame$GO.ID,')',sep='')
  
          markobj <- c(markobj,'### BP Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(bp.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
        
        }else{markobj <- c(markobj,'### BP Ontology GOs\n There are no GO terms related to this set of DEGs for BP ontology\n')}
        
        if(dim(GOsMatrix$`MF Ontology GOs`)[1] != 0){
          GOsMatrix$`MF Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`MF Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
          mf.frame <- GOsMatrix$`MF Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
          rownames(mf.frame) <- NULL
          mf.frame$GO.ID  <- paste('[',mf.frame$GO.ID,'](',amigo.url,mf.frame$GO.ID,')',sep='')
          
          markobj <- c(markobj,'### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(mf.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
  
        }else{markobj <- c(markobj,'### MF Ontology GOs\n There are no GO terms related to this set of DEGs for MF ontology\n')}
        
        if(dim(GOsMatrix$`CC Ontology GOs`)[1] != 0){
          GOsMatrix$`CC Ontology GOs`$`Gene Symbols` <- as.character(lapply(GOsMatrix$`CC Ontology GOs`$`Gene Symbols`, function(x) {gsub(",", ", ", x)}))
          cc.frame <- GOsMatrix$`CC Ontology GOs`[,c('GO.ID','Term','Description','Gene Symbols')]
          rownames(cc.frame) <- NULL
          cc.frame$GO.ID  <- paste('[',cc.frame$GO.ID,'](',amigo.url,cc.frame$GO.ID,')',sep='')
          markobj <- c(markobj,'### MF Ontology GOs\n','```{r echo=FALSE}',paste('knitr::kable(cc.frame,"',table.format,'", table.attr = "class=\'paleBlueRows\'")',sep=''),'```\n')
        }else{markobj <- c(markobj,'### CC Ontology GOs\n There are no GO terms related to this set of DEGs for CC ontology\n')}
        
      }

      # --- Pathways Visualization --- #
      if(getPathways){
        
        DEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix),notHSapiens=FALSE, attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"), filter = "external_gene_name")

        DEGsAnnotation <- DEGsAnnotation[which(!is.na(DEGsAnnotation$entrezgene_id) == TRUE),]
        
        commonDEGs <- intersect(rownames(DEGsMatrix),unique(DEGsAnnotation$external_gene_name))
        posCommonDEGs <- match(rownames(DEGsMatrix[commonDEGs,]),DEGsAnnotation$external_gene_name)
        pathMatrix <- DEGsMatrix[commonDEGs,]
        rownames(pathMatrix) <- DEGsAnnotation$entrezgene_id[posCommonDEGs]
        
        paths.data <- matrix(ncol = 3)
        
        cat("Retrieving DEGs associated pathways...\n")
        
        for(gene.index in seq(dim(DEGsAnnotation)[1])){
          gene <- DEGsAnnotation[gene.index,'entrezgene_id']
          gene.name <- DEGsAnnotation[gene.index,'external_gene_name']
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
          diseases <- DEGsToDiseases(rownames(DEGsMatrix), size = 5, getEvidences = TRUE)
          evidences.frame <- list()
          for (gene in names(diseases)){
            act.markobj <- c()
            # If user want to see evidences for all diseases or this diseases match with solicitated disease
            check.diseases <- names(diseases[[gene]]$evidences)
    
            for (act.disease in check.diseases){
              if ( class(diseases[[gene]]$evidences[[act.disease]]) == 'list' ){
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
          markobj <- c(markobj,'## Related diseases\n',
                       'Finally, the related diseases enrichment is displayed. DEGs related diseases are searched 
                      from *targetValidation* plastform.\n')
    
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
                  if (class(evidences[[gene]])=='list'){
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

