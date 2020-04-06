#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function generates different plots over expression data in order to detect possible outliers.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @param toPNG Boolean variable to indicate if a plot would be save to PNG.
#' @param toPDF Boolean variable to indicate if a plot would be save to PDF.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' RNAseqQA(expressionMatrix)


RNAseqQA <- function(expressionMatrix, outdir = "myPlots", toPNG = TRUE, toPDF = TRUE){
  
  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  if(!is.logical(toPNG)){stop("toPNG parameter can only take the values TRUE or FALSE.")}
  if(!is.logical(toPDF)){stop("toPDF parameter can only take the values TRUE or FALSE.")}
  
  if(! dir.exists(outdir)) system(paste("mkdir",outdir))
  
  cat("Performing samples quality analysis...\n")
  
  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]
  num.samples <- ncol(expressionMatrix)
  
  # Array distances matrix
  distance.matrix <- as.matrix(dist(t(expressionMatrix),'manhattan'))/nrow(expressionMatrix)
  distance.matrix <- melt(distance.matrix)

  distance.plot <- ggplot(data = distance.matrix, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + ggtitle('Distances between arrays') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('')
  #distance.plot
  if (toPNG) ggsave(paste(outdir,'distance-plot.png',sep='/'),distance.plot)
  if (toPDF)  ggsave(paste(outdir,'distance-plot.pdf',sep='/'),distance.plot)

  # BOXPLOT
  boxplot.data <- melt(expressionMatrix)
  colnames(boxplot.data) <- c('Tags','Samples','Expression')
  box.plot <- ggplot(boxplot.data, aes(x=Expression , y=Samples)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,fill='#56B4E9') 
  #box.plot
  if (toPNG) ggsave(paste(outdir,'box-plot.png',sep='/'),box.plot)
  if (toPDF) ggsave(paste(outdir,'box-plot.pdf',sep='/'),box.plot)
  
  # KS
  # Empirical cumulative distribution function
  # AÑADIR ecdf DE stats
  fx = ecdf(as.vector(expressionMatrix))
  r <- suppressWarnings(apply(expressionMatrix, 2, function(v)
    ks.test(v, y = fx, alternative='two.sided')$statistic))
  ks.data <- data.frame('ks'=as.vector(r),'x'=c(1:length(r)))
  ks.plot <- ggplot(ks.data,aes(x=ks,y=x)) + geom_point(color = "#56B4E9") + 
    geom_segment(aes(yend=x),xend=0,color = "#56B4E9") +
    scale_y_discrete(limits=c(1:num.samples),labels=colnames(expressionMatrix)) + ylab('') +
    theme_classic()
  #ks.plot
  if (toPNG) ggsave(paste(outdir,'ks-plot.png',sep='/'),ks.plot)
  if (toPDF) ggsave(paste(outdir,'ks-plot.pdf',sep='/'),ks.plot)
  
  # MA-PLOT
  ma.plot <- function(){
    par(mfrow=c(ceiling(num.samples/3),3))
    par(mar=c(1,2,1,1),mgp=c(3,0.3,0))
    for ( i in seq(num.samples)){
      act.mean <- rowMeans(expressionMatrix[,setdiff(seq(num.samples),i)])
      M <- expressionMatrix[,i] - act.mean
      A = rowMeans(cbind(expressionMatrix[,i],act.mean))
      act.data  <- data.frame('M'=M,'A'=A)
      d <- hoeffd(cbind(A,M))$D[1,2]
      d <- round(d,3)
      smoothScatter(A,M,main=paste(colnames(expressionMatrix)[i],'.D =',d),xlab='',ylab='')
      title(xlab='A',line=0)
      title(ylab='M',line=1)
    }
    par(mfrow=c(1,1))
  }
  if (toPNG){
    png("MA-plot.png",units="in", width=5, height=5, res=300)
    ma.plot()
    dev.off()
  }
  if (toPDF){
    pdf("MA-plot.pdf")
    ma.plot()
    dev.off()
  }

}
