#' Plot different graphs depending on the current step of KnowSeq pipeline.
#'
#' This function allows to plot different charts only by changing the parameters, for the different KnowSeq pipeline steps. Furthermore, the chosen plot can be saved to PNG and PDF.
#'
#' @param data Normally, the data parameter is an expression matrix or data.frame, however for the confusionMatrix plot, the data are a confussion matrix that can be achieved by using the output of any of the machine learning functions of this package.
#' @param labels A vector or factor that contains the labels for each of the samples in the data parameter.
#' @param colours A vector that contains the desired colours to plot the different charts. Example: c("red","green","blue").
#' @param main The title for the plot.
#' @param ylab The description for the y axis.
#' @param xlab The description for the x axis.
#' @param xgrid Shows the x grid into the plot
#' @param ygrid Shows the y grid into the plot
#' @param legend A vector with the elements in the legend of the plot.
#' @param mode The different plots supported by this package. The possibilities are boxplot, orderedBoxplot, genesBoxplot, heatmap, confusionMatrix and classResults.
#' @param toPNG Boolean variable to indicate if a plot would be save to PNG.
#' @param toPDF Boolean variable to indicate if a plot would be save to PDF.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' dataPlot(expressionMatrix,labels,mode = "boxplot",toPNG = TRUE,toPDF = TRUE)
#' dataPlot(DEGsMatrix[1:12,],labels,mode = "orderedBoxplot",toPNG = TRUE,toPDF = TRUE)
#' dataPlot(DEGsMatrix[1:12,],labels,mode = "genesBoxplot",toPNG = TRUE,toPDF = FALSE)
#' dataPlot(DEGsMatrix[1:12,],labels,mode = "heatmap",toPNG = TRUE,toPDF = TRUE)


dataPlot <- function(data, labels, colours = c("green", "red"), main = "", ylab = "Expression", xlab = "Samples", xgrid = FALSE, ygrid = FALSE, legend = "", mode="boxplot", toPNG = FALSE, toPDF = FALSE){
  
  if(!is.logical(toPNG)){stop("toPNG parameter can only take the values TRUE or FALSE.")}
  if(!is.logical(toPDF)){stop("toPDF parameter can only take the values TRUE or FALSE.")}
  
  # dev.new()
  
  if(mode == "boxplot"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
    }else{
      
      coloursPalette <- colours
      
    }
    
    coloursAux <- labels
    
    for(i in seq_len(length(levels(as.factor(labels))))){
      
      coloursAux <- gsub(levels(as.factor(labels))[i],coloursPalette[i], coloursAux)
      
    }
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("boxplot.png",width = 1024, height = 720)
      boxplot(data, col=coloursAux, ylab=ylab, xlab=xlab, main=main)
      
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("boxplot.pdf")
      boxplot(data, col=coloursAux, ylab=ylab, xlab=xlab, main=main)
      
      dev.off()
    }
    
    boxplot(data, col=coloursAux, ylab=ylab, xlab=xlab, main=main)
    
    
  }else if(mode == "orderedBoxplot"){
    
    sortedLabels <- order(labels)
    sortedLabelsNames <- labels[sortedLabels]
    
    if(length(levels(as.factor(sortedLabelsNames))) != 2 && length(levels(as.factor(sortedLabelsNames))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],length(levels(as.factor(sortedLabelsNames))) - 2))
      
    }else{
      
      coloursPalette <- colours
      
    }
    
    coloursAux <- sortedLabelsNames
    
    for(i in seq_len(length(levels(as.factor(sortedLabelsNames))))){
      
      coloursAux <- gsub(levels(as.factor(sortedLabelsNames))[i],coloursPalette[i], coloursAux)
      
    }
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("orderedBoxplot.png",width = 1024, height = 720)
      boxplot(data[,sortedLabels] , col=coloursAux, ylab=ylab , xlab=xlab, main=main)
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("orderedBoxplot.pdf")
      boxplot(data[,sortedLabels] , col=coloursAux, ylab=ylab , xlab=xlab, main=main)
      
      dev.off()
    }
    
    boxplot(data[,sortedLabels] , col=coloursAux, ylab=ylab , xlab=xlab, main=main)
    
    
  }else if(mode == "genesBoxplot"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
    }else{
      
      coloursPalette <- colours
      
    }
    
    
    meltMatrix <- t(data)
    rownames(meltMatrix) <- labels
    
    if(ncol(meltMatrix) < 24){
      col = ncol(meltMatrix)
    }else{
      col = 24
    }
    
    xx <- melt(meltMatrix[,seq_len(col)])
    names(xx) <- c("Classes", "Gen", "Value")
    
    print(ggplot(xx, aes(x=as.factor(Classes),y=Value,fill=as.factor(Classes))) + geom_boxplot() + facet_wrap(~Gen, ncol = 3) 
          + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylab(ylab) + labs(fill = "Classes"))
    
    if(toPNG){
      cat("Creating PNG...\n")
      ggplot(xx, aes(x=as.factor(Classes),y=Value,fill=as.factor(Classes))) + geom_boxplot() + facet_wrap(~Gen, ncol = 3) + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylab(ylab) + labs(fill = "Classes")
      
      ggsave("genesBoxplot.png", width = 15, height = 10)
      
    }
    if(toPDF){
      cat("Creating PDF...\n")
      ggplot(xx, aes(x=as.factor(Classes),y=Value,fill=as.factor(Classes))) + geom_boxplot() + facet_wrap(~Gen, ncol = 3) + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylab(ylab) + labs(fill = "Classes")
      
      ggsave("genesBoxplot.pdf", width = 15, height = 10)
      
    }
    
  }else if(mode == "heatmap"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
    }else{
      
      coloursPalette <- colours
      
    }
    
    coloursAux <- labels
    
    for(i in seq_len(length(levels(as.factor(labels))))){
      
      coloursAux <- gsub(levels(as.factor(labels))[i],coloursPalette[i], coloursAux)
      
    }
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("heatmap.png",width = 1024, height = 720)
      heatmap.2(t(data), col=redgreen(75), scale="row", RowSideColors=coloursAux,
                key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1,
                cexCol=1, margins=c(6,11),srtCol=45)
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("heatmap.pdf")
      heatmap.2(t(data), col=redgreen(75), scale="row", RowSideColors=coloursAux,
                key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1,
                cexCol=1, margins=c(6,11),srtCol=45)
      dev.off()
    }
    
    heatmap.2(t(data), col=redgreen(75), scale="row", RowSideColors=coloursAux ,
              key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1,
              cexCol=1, margins=c(6,11),srtCol=45)
    
  }else if(mode == "confusionMatrix"){
    
    plotConfMatrix(data)
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("confusionMatrix.png",width = 1024, height = 720)
      
      plotConfMatrix(data)
      
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("confusionMatrix.pdf")
      
      plotConfMatrix(data)
      
      dev.off()
    }
    
  }else if(mode == "classResults"){
    
    if(!is.matrix(data)){
      
      plot(data,type='l',col=colours[2], main=main,
           xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
      
      if(xgrid){
        grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
      }
      if(ygrid){
        grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
      }
      
      if(legend != ""){
        legend("bottomright", legend=legend, col=colours[2],lwd = 2.5,cex=0.8,lty = "solid")
      }
      
    }else if(is.matrix(data)){
      
      if(length(colours) != dim(data)[1]){ colours = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],dim(data)[1])}
      
      plot(data[1,],type='l',col=colours[1], main=main,
           xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5, ylim = c(min(data),max(data)))
      
      for(i in c(2:dim(data)[1])){
        
        lines(data[i,],col=colours[i], lty="solid",lwd = 2.5)
        
      }
      
      if(xgrid){
        grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
      }
      if(ygrid){
        grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
      }
      
      if(legend == ""){legend = rownames(data)}
      legend("bottomright", legend=legend,
             col=colours, lty=rep("solid",dim(data)[2]), lwd = rep(2.5,dim(data)[2]),cex=0.8)
      
    }
    
    
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("classResults.png",width = 1024, height = 720)
      
      if(!is.matrix(data)){
        
        plot(data,type='l',col=colours[2], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        if(legend != ""){
          legend("bottomright", legend=legend, col=colours[2],lwd = 2.5,cex=0.8,lty = "solid")
        }
        
      }else if(is.matrix(data)){
        
        if(length(colours) != dim(data)[1]){ colours = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],dim(data)[1])}
        
        plot(data[1,],type='l',col=colours[1], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5, ylim = c(min(data),max(data)))
        
        for(i in c(2:dim(data)[1])){
          
          lines(data[i,],col=colours[i], lty="solid",lwd = 2.5)
          
        }
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        if(legend == ""){legend = rownames(data)}
        legend("bottomright", legend=legend,
               col=colours, lty=rep("solid",dim(data)[2]), lwd = rep(2.5,dim(data)[2]),cex=0.8)
        
      }
      
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("classResults.pdf")
      
      if(!is.matrix(data)){
        
        plot(data,type='l',col=colours[2], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        
        if(legend != ""){
          legend("bottomright", legend=legend, col=colours[2],lwd = 2.5,cex=0.8,lty = "solid")
        }
        
      }else if(is.matrix(data)){
        
        if(length(colours) != dim(data)[1]){ colours = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)],dim(data)[1])}
        
        plot(data[1,],type='l',col=colours[1], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5, ylim = c(min(data),max(data)))
        
        for(i in c(2:dim(data)[1])){
          
          lines(data[i,],col=colours[i], lty="solid",lwd = 2.5)
          
        }
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        if(legend == ""){legend = rownames(data)}
        legend("bottomright", legend=legend,
               col=colours, lty=rep("solid",dim(data)[2]), lwd = rep(2.5,dim(data)[2]),cex=0.8)
        
      }
      
      dev.off()
    }
    
  }else{
    stop("Mode unrecognized. Please, use a mode listed in the documentation.")
  }
}
