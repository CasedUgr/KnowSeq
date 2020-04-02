#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function adapts the RNA-seq data in order to allows using arrayQualityMetrics expression analysis.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
RNAseqQA(expressionMatrix)

RNAseqQA <- function(expressionMatrix, outdir = "RNAseqQA"){
  
  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  
  cat("Performing samples quality analysis...\n")
  
  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]
  num.samples <- num.samples
  
  # Array distances matrix
  melted_cormat <- melt(cor(expressionMatrix))
  distance.plot <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + ggtitle('Distances between arrays') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('')
  #distance.plot
  ggsave(paste(outdir,'distance-plot.png',sep='/'),distance.plot)


  # BOXPLOT
  boxplot.data <- melt(expressionMatrix)
  colnames(boxplot.data) <- c('Tags','Samples','Expression')
  box.plot <- ggplot(boxplot.data, aes(x=Samples, y=Expression)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,fill='#56B4E9') 
  #box.plot
  ggsave(paste(outdir,'box-plot.png',sep='/'),box.plot)
  
  
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
  ggsave(paste(outdir,'ks-plot.png',sep='/'),ks.plot)
  
  # MA-PLOT
  plots <- list()
  for ( i in seq(dim(expressionMatrix)[2])){
    act.mean <- rowMeans(expressionMatrix[,setdiff(seq(num.samples),i)])
    M <- expressionMatrix[,i] - act.mean
    A = rowMeans(cbind(expressionMatrix[,i],act.mean))
    act.data  <- data.frame('M'=M,'A'=A)
    d <- hoeffd(cbind(A,M))$D[1,2]
    plots[[i]] <- ( ggplot(act.data,aes(A,M)) +  geom_point(shape=1) +  
                      ggtitle(paste(colnames(expressionMatrix)[i],'. D =',d)) )
  }
  ma.plots <- do.call(grid.arrange,c(plots,ncol=2))
  #ma.plots
  ggsave(paste(outdir,'MA-plot.png',sep='/'),ma.plots)
}
