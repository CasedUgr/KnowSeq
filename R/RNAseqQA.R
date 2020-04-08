#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function generates different plots over expression data in order to detect possible outliers.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @param toPNG Boolean variable to indicate if a plot would be save to PNG.
#' @param toPDF Boolean variable to indicate if a plot would be save to PDF.
#' @param D.limit Numeric variable to indicate Hoeffding's statistic, D, limit to a sample to be outlier
#' @param KS.limit Numeric variable to indicate Kolmogorov-Smirnov statistic, KS, limit to a sample to be outlier
#' @return A list containing found outliers for each realized test. 
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' outliers <- RNAseqQA(expressionMatrix)

RNAseqQA <- function(expressionMatrix, outdir = "SamplesQualityAnalysis", toPNG = TRUE, toPDF = TRUE, D.limit = 0.1 , KS.limit = 0.1){
  
  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  if(!is.logical(toPNG)){stop("toPNG parameter can only take the values TRUE or FALSE.")}
  if(!is.logical(toPDF)){stop("toPDF parameter can only take the values TRUE or FALSE.")}
  
  if(! dir.exists(outdir)) system(paste("mkdir",outdir))
  
  cat("Performing samples quality analysis...\n")
  
  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]
  num.samples <- ncol(expressionMatrix)
  outliers <- list()
  
  outlierBarPlot <- function(data,title,limit,xlab){
    yticks=data$y
    outlier.bar.plot <-  ggplot(data,aes(x=x,y=y)) + geom_point(color = "#56B4E9") +
      geom_segment(aes(yend=y),xend=0,color = "#56B4E9") + ylab('Samples') + xlab(xlab) +
      ggtitle(title) + theme_classic() + geom_vline(xintercept = limit) +
      scale_y_discrete(breaks = yticks[seq(1,length(yticks),by=3)],labels=yticks[seq(1,length(yticks),by=3)])
    return(outlier.bar.plot)
  }

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

  # Distance outlier detection
  dist.outlier.plot <- outlierBarPlot(dist.data,'Distance based Outliers',min(dist.data$x)*2,'Sum-Distance')
  
  if (toPNG) ggsave(paste(outdir,'distance-outlier-plot.png',sep='/'),dist.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF)  ggsave(paste(outdir,'distance-outlier-plot.pdf',sep='/'),dist.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  
  outliers[['Distance']] <- distance.sum[which(distance.sum > min(dist.data$x)*2)]
  
  # --- --- BOXPLOT --- --- #
  quantiles <-  apply(expressionMatrix,2, quantile)
  min.x <- max(quantiles[1,])
  max.x <- min(quantiles[5,])

  boxplot.data <- melt(expressionMatrix)
  colnames(boxplot.data) <- c('Tags','Samples','Expression')

  box.plot <- ggplot(boxplot.data, aes(x=Expression , y=Samples)) + 
    scale_y_discrete(breaks = levels(seq(nrow(boxplot.data)))[floor(seq(1, nlevels(seq(nrow(boxplot.data))),length.out = 10))]) +
    xlim(min.x,max.x) + geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,fill='#56B4E9') 
  #box.plot
  
  
  if (toPNG) ggsave(paste(outdir,'box-plot.png',sep='/'),box.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF) ggsave(paste(outdir,'box-plot.pdf',sep='/'),box.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)


  # KS
  # Empirical cumulative distribution function
  # AÑADIR ecdf DE stats
  fx = ecdf(as.vector(expressionMatrix))
  ks <- suppressWarnings(apply(expressionMatrix, 2, function(v)
    ks.test(v, y = fx, alternative='two.sided')$statistic))

  ks.data <- data.frame('x'=ks,'y'=c(1:length(ks)))
  ks.plot <- outlierBarPlot(ks.data,'KS - Outliers',KS.limit,'KS')

  #ks.plot

  if (toPNG) ggsave(paste(outdir,'ks-plot.png',sep='/'),ks.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF) ggsave(paste(outdir,'ks-plot.pdf',sep='/'),ks.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)

  outliers[['KS']] <- ks[which(ks > KS.limit)]
  
  # --- --- MA --- --- #
  M <- list()
  A <- list()
  Da <- list()
  
  for ( i in seq(num.samples)){
    act.mean <- rowMeans(expressionMatrix[,setdiff(seq(num.samples),i)])
    M[[i]] <- expressionMatrix[,i] - act.mean
    A[[i]] = rowMeans(cbind(expressionMatrix[,i],act.mean))
    Da[[i]] <- round(hoeffd(cbind(A[[i]],M[[i]]))$D[1,2],3)
  }
  
  Da <- unlist(Da)
  names(Da) <- colnames(expressionMatrix)
  sort.da <- sort(Da)

  ma.plot <- function(ma.data){
    par(mfrow=c(3,3))
    par(mar=c(1,2,1,1),mgp=c(3,0.3,0))
    for ( i in seq(length(ma.data))){
      act.data  <- data.frame('M'=M,'A'=A)
      smoothScatter(A[[i]],M[[i]],main=paste(names(ma.data)[i],'.D =',ma.data[i]),xlab='',ylab='')
      title(xlab='A',line=0)
      title(ylab='M',line=1)
    }
    par(mfrow=c(1,1))
  }
  
  if (toPNG){
    png(paste(outdir,'MA-plot.png',sep='/'),units="in", width=5, height=5, res=300)
    ma.plot(sort.da[1:9])
    dev.off()
  }
  if (toPDF){
    pdf(paste(outdir,'MA-plot.pdf',sep='/'))
    ma.plot(sort.da[1:9])
    dev.off()
  }

  # MA - Outliers
  ma.data <- data.frame('x'=unlist(Da),'y'=seq(length(Da)))
  ma.outlier.plot <- outlierBarPlot(ma.data,'MA - Outliers',D.limit,'D')
  
  
  if (toPNG) ggsave(paste(outdir,'MA-outlier-plot.png',sep='/'),ma.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)
  if (toPDF) ggsave(paste(outdir,'MA-outlier-plot.pdf',sep='/'),ma.outlier.plot,width=5, height=ceiling(ncol(expressionMatrix)/5),limitsize=FALSE,units = "in", dpi = 300)

  outliers[['MA-D']] <- Da[which(Da > D.limit)]

  return(outliers)

}

