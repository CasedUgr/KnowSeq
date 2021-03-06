#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function generates different plots over expression data in order to detect possible outliers.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @param toPNG Boolean variable to indicate if a plot would be save to PNG.
#' @param toPDF Boolean variable to indicate if a plot would be save to PDF.
#' @param toRemoval Boolean variable to indicate if detected outliers will be removed from original data
#' @return A list containing found outliers for each realized test or corrected data if toRemoval is TRUE. 
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' outliers <- RNAseqQA(expressionMatrix)

RNAseqQA <- function(expressionMatrix, outdir = "SamplesQualityAnalysis", toPNG = TRUE, toPDF = TRUE, toRemoval = FALSE){

  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  if(!is.logical(toPNG)){stop("toPNG parameter can only take the values TRUE or FALSE.")}
  if(!is.logical(toPDF)){stop("toPDF parameter can only take the values TRUE or FALSE.")}
  
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  cat("Performing samples quality analysis...\n")
  
  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]
  outliers <- list()
  found.outliers <- c(-1)
  removed.outliers <- c()

  outlierBarPlot <- function(data,title,limit,xlab){
    yticks=data$y
    outlier.bar.plot <-  ggplot(data,aes(x=x,y=y)) + geom_point(color = "#56B4E9") +
      geom_segment(aes(yend=y),xend=0,color = "#56B4E9") + ylab('Samples') + xlab(xlab) +
      ggtitle(title) + theme_classic() + geom_vline(xintercept = limit) +
      scale_y_discrete(breaks = yticks[seq(1,length(yticks),by=3)],labels=yticks[seq(1,length(yticks),by=3)])
    return(outlier.bar.plot)
  }

  num.samples <- ncol(expressionMatrix)
  
  cat("Running Distances Outliers Detection test...")
  # --- --- SAMPLES DISTANCES --- --- #
  # Array distances matrix
  distance.matrix <- as.matrix(dist(t(expressionMatrix),'manhattan'))/nrow(expressionMatrix)
  distance.sum  <- colSums(distance.matrix)
  dist.data <-  data.frame('x'=distance.sum,'y'=colnames(distance.matrix))

  distance.matrix <- melt(distance.matrix)
  num.data <- ncol(distance.matrix)

  distance.plot <- ggplot(data = distance.matrix, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + ggtitle('Distances between arrays') + 
    scale_y_discrete(breaks = levels(seq(num.data))[floor(seq(1, nlevels(seq(num.data)),length.out = 10))]) +
    scale_x_discrete(breaks = levels(seq(num.data))[floor(seq(1, nlevels(seq(num.data)),length.out = 10))]) +
    xlab('') + ylab('')

  if (toPNG) ggsave(paste(outdir,'distance-plot.png',sep='/'),distance.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF)  ggsave(paste(outdir,'distance-plot.pdf',sep='/'),distance.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)

  q3 <- quantile(dist.data$x)[4]
  dist.limit <- q3 + 1.5 * IQR(dist.data$x)

  # Distance outlier detection
  dist.outlier.plot <- outlierBarPlot(dist.data,'Distance based Outliers',dist.limit,'Sum-Distance')
  
  if (toPNG) ggsave(paste(outdir,'distance-outlier-plot.png',sep='/'),dist.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF)  ggsave(paste(outdir,'distance-outlier-plot.pdf',sep='/'),dist.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  
  outliers[['Distance']] <- list('limit'=dist.limit,outliers=distance.sum[which(distance.sum > dist.limit)])
  
  cat("Done!\n")
  
  cat("Running Kolmogorov-Smirnov test...")
  # --- --- KS --- --- #
  # KS
  # Empirical cumulative distribution function
  fx = ecdf(as.vector(expressionMatrix))
  ks <- suppressWarnings(apply(expressionMatrix, 2, function(v)
    ks.test(v, y = fx, alternative='two.sided')$statistic))
  
  ks.data <- data.frame('x'=ks,'y'=seq_len(length(ks)))
  q3 <- quantile(ks.data$x)[4]
  KS.limit <- q3 + 1.5 * IQR(ks.data$x)
  ks.plot <- outlierBarPlot(ks.data,'KS - Outliers',KS.limit,'KS')
  
  ks.outliers.index <- which(ks > KS.limit)
  
  quantiles <-  apply(expressionMatrix,2, quantile)
  min.x <- max(quantiles[1,])
  max.x <- min(quantiles[5,])
  
  boxplot.data <- melt(expressionMatrix)
  colnames(boxplot.data) <- c('Tags','Samples','Expression')
  
  box.plot <- ggplot(boxplot.data, aes(x=Expression , y=Samples)) + 
    scale_y_discrete(breaks = levels(seq(nrow(boxplot.data)))[floor(seq(1, nlevels(seq(nrow(boxplot.data))),length.out = 10))]) +
    xlim(min.x,max.x) + geom_boxplot(outlier.shape=NA,fill='#56B4E9') 
  
  for ( i in seq(length(ks.outliers.index)))  {
    box.plot <-  box.plot + annotate(geom="point", x = min.x+0.1, y = ks.outliers.index[i], colour="red",size=3)
  }
  

  options(warn=-1)
  if (toPNG) {
    ggsave(paste(outdir,'box-plot.png',sep='/'),box.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
    ggsave(paste(outdir,'ks-plot.png',sep='/'),ks.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  }
  if (toPDF) {
    ggsave(paste(outdir,'box-plot.pdf',sep='/'),box.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
    ggsave(paste(outdir,'ks-plot.pdf',sep='/'),ks.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  }
  options(warn=0)

  outliers[['KS']] <- list('limit'=KS.limit,'outliers'=ks[which(ks > KS.limit)])
  
  cat("Done!\n")
  
  cat("Running MAD Outliers Detection test...")
  
  # --- --- MA --- --- #
  outliersMA <- list()
  rowExpression <- colMeans(expressionMatrix)
  
  for (i in seq(length(rowExpression))){
    exprMatrix <- rowExpression[-i]
    upperBound <- median(exprMatrix) + 3 * mad(exprMatrix, constant = 1)
    lowerBound <- median(exprMatrix) - 3 * mad(exprMatrix, constant = 1)
    
    outliersMA[[i]] <- which(rowExpression[i] < lowerBound | rowExpression[i] > upperBound)
  }
  

  outliers[['MAD']] <- list('limit'="-",'outliers'=unlist(outliersMA))

  cat("Done!\n")
  
  if (! toRemoval) return(outliers) 

  found.outliers <- union(intersect(names(outliers[[1]]$outliers),names(outliers[[2]]$outliers)),
                  union(intersect(names(outliers[[1]]$outliers),names(outliers[[3]]$outliers)),
                  intersect(names(outliers[[2]]$outliers),names(outliers[[3]]$outliers))))
  expressionMatrix <- expressionMatrix[, ! colnames(expressionMatrix) %in% found.outliers]
  removed.outliers <- c(removed.outliers,found.outliers)

  return(list('matrix'=expressionMatrix,'outliers'=removed.outliers))
}

