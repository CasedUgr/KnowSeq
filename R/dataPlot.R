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
#' @param mode The different plots supported by this package. The possibilities are boxplot, orderedBoxplot, genesBoxplot, heatmap, confusionMatrix, classResults and heatmapResults.
#' @param heatmapResultsN Number of genes to show if mode is heatmapResults.
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

dataPlot <- function(data, labels, colours = c("red", "green"), main = "", ylab = "Expression", xlab = "Samples", xgrid = FALSE, ygrid = FALSE, legend = "", mode="boxplot", heatmapResultsN = 0, toPNG = FALSE, toPDF = FALSE){
  
  if(!is.logical(toPNG)){stop("toPNG parameter can only take the values TRUE or FALSE.")}
  if(!is.logical(toPDF)){stop("toPDF parameter can only take the values TRUE or FALSE.")}
  
  # dev.new()
  
  if(mode == "boxplot"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
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
      
      coloursPalette <- c("green","red",sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],length(levels(as.factor(sortedLabelsNames))) - 2))
      
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
    
  }else if(mode == "heatmapResults"){
    
    # data is the output of svm_trn, knn_trn or rf_trn
    acc <- data$accuracyInfo$meanAccuracy
    sen <- data$sensitivityInfo$meanSensitivity
    spe <- data$specificityInfo$meanSpecificity
    f1 <- data$F1Info$meanF1
    
    data2 <- data.frame(acc, sen, spe, f1)
    
    
    # Only showing first heatmapResultsN genes
    heatmapResultsN <- ifelse(heatmapResultsN != 0, heatmapResultsN, nrow(data2))
    data2 <- data2[1:heatmapResultsN, ]
    names(data2) <- c("Accuracy", "Sensitivity", "Specificity", "F1-Score")
    data2 <- melt(t(data2))
    
    graph <- ggplot(data2, aes(Var1, Var2, fill = value)) +
      geom_tile(colour = "black") + 
      scale_fill_gradient(low = colours[1], high = colours[2]) +
      labs(x = ifelse(xlab != "Samples", xlab, ""),
           y = ifelse(ylab != "Expression", ylab, "Number of genes used"),
           fill = "",
           title = main,
           subtitle = "") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = .5),
            panel.grid = element_blank(),
            panel.spacing = unit(c(0, 0, 0, 0), "null"),
            axis.text.x = element_text(size = 11),
            panel.grid.major.x = element_blank(),
            axis.text = element_text(color = "black"))
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("heatmapResults.png",width = 1024, height = 720)
      print(graph)
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("heatmapResults.pdf")
      print(graph)
      dev.off()
    }
    
    print(graph)
    
  }else if(mode == "genesBoxplot"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
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
          + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylim(min(xx$Value),max(xx$Value)) + ylab(ylab) + labs(fill = "Classes") + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90)))
    
    if(toPNG){
      cat("Creating PNG...\n")
      ggplot(xx, aes(x=as.factor(Classes),y=Value,fill=as.factor(Classes))) + geom_boxplot() + facet_wrap(~Gen, ncol = 3) + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylim(min(xx$Value),max(xx$Value)) + ylab(ylab) + labs(fill = "Classes") + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90))      
      ggsave("genesBoxplot.png", width = 15, height = 10)
      
    }
    if(toPDF){
      cat("Creating PDF...\n")
      ggplot(xx, aes(x=as.factor(Classes),y=Value,fill=as.factor(Classes))) + geom_boxplot() + facet_wrap(~Gen, ncol = 3) + scale_fill_manual(values=coloursPalette) + ggtitle(main) + xlab(xlab) + ylim(min(xx$Value),max(xx$Value)) + ylab(ylab) + labs(fill = "Classes") + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90))      
      ggsave("genesBoxplot.pdf", width = 15, height = 10)
      
    }
    
  }else if(mode == "heatmap"){
    
    if(length(levels(as.factor(labels))) != 2 && length(levels(as.factor(labels))) != length(colours)){
      
      coloursPalette <- c("green","red",sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],length(levels(as.factor(labels))) - 2))
      
    }else{
      
      coloursPalette <- colours
      
    }
    
    # Shorten long labels if one of them is longer than 20 characters
    limit_characters <- 20
    if(sum(nchar(labels) > limit_characters) > 0){
      labels <- substr(labels, 1, limit_characters)
      cat(paste0("Warning! The labels were too long (>", limit_characters, " characters) and were shortened for the plot. \n"))
    }
    
    # Reorder data and labels
    labels_levels <- levels(as.factor(labels))
    for(i in 1:length(labels_levels)){
      if(i == 1) {
        data_heatmap <- data[, which(labels == labels_levels[1])]
        labels_heatmap <- labels[which(labels == labels_levels[1])]
      }
      else {
        data_heatmap <- cbind(data_heatmap, data[, which(labels == labels_levels[i])])
        labels_heatmap <- c(labels_heatmap, labels[which(labels == labels_levels[i])])
      }
    }
    
    # Prepare data for ggplot2::geom_tile
    colnames(data_heatmap) <- seq_len(ncol(data_heatmap))
    data_heatmap <- melt(data_heatmap)
    
    # Add labels_heatmap to data_heatmap
    data_heatmap$labels <- NA
    j <- 1
    for(i in 1:nrow(data_heatmap)){
      if(i == 1) data_heatmap$labels[i] <- labels_heatmap[1]
      else {
        if(data_heatmap$Var2[i] == data_heatmap$Var2[i - 1]) data_heatmap$labels[i] <- data_heatmap$labels[i - 1]
        else{
          j <- j + 1
          data_heatmap$labels[i] <- labels_heatmap[j]
        }
      }
    }
    
    # First graph: labels
    g1 <- ggplot(data = data_heatmap, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill = labels), size = 0, linetype = "blank") +
      scale_fill_manual(values = coloursPalette) +
      ggtitle("") + 
      xlab("") + 
      ylab("") + 
      coord_flip() + 
      theme_minimal() + 
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "white"),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust = .5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "left",
            legend.margin=margin(t = 0, unit='cm'))
    
    # Second graph: data
    g2 <- ggplot(data = data_heatmap, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill = value), size = 0, linetype = "blank") +
      scale_fill_gradientn(colors = c("red", "black", "green")) +
      ggtitle(main) + 
      xlab(xlab) + 
      ylab(ylab) + 
      coord_flip() + 
      theme_minimal() + 
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = .5))
    
    # Join the two graphs
    grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(1, 1.5))
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("heatmap.png", width = 1024, height = 720)
      grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(1, 1.5))
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("heatmap.pdf")
      grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(1, 1.5))
      dev.off()
    }
    
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
    c_palette <- c('red', 'blue', 'green', 'orange', 'yellow', 'brown', 'purple', 'pink', 'tan', 'sienna')
    if(!is.matrix(data)){
      
      plot(data,type='l',col=colours[1], main=main,
           xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
      
      if(xgrid){
        grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
      }
      if(ygrid){
        grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
      }
      
      if(length(legend) != 0){
        legend("bottomright", legend=legend, col=colours[1],lwd = 2.5,cex=0.8,lty = "solid")
      }
      
    }else if(is.matrix(data)){
      
      if(length(colours) == dim(data)[1]){
        colours = colours
      }else if(length(colours) != dim(data)[1] && dim(data)[1] <= length(c_palette)){
        colours = c_palette[1:dim(data)[1]]
      }else {
        colours = sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],dim(data)[1])
      }
      
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
      
      if(length(legend) == 0){legend = rownames(data)}
      legend("bottomright", legend=legend,
             col=colours, lty=rep("solid",dim(data)[2]), lwd = rep(2.5,dim(data)[2]),cex=0.8)
      
    }
    
    
    
    if(toPNG){
      cat("Creating PNG...\n")
      png("classResults.png",width = 1024, height = 720)
      
      if(!is.matrix(data)){
        
        plot(data,type='l',col=colours[1], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        if(length(legend) != 0){
          legend("bottomright", legend=legend, col=colours[2],lwd = 2.5,cex=0.8,lty = "solid")
        }
        
      }else if(is.matrix(data)){
        
        if(length(colours) == dim(data)[1]){
          colours = colours
        }else if(length(colours) != dim(data)[1] && dim(data)[1] <= length(c_palette)){
          colours = c_palette[1:dim(data)[1]]
        }else {
          colours = sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],dim(data)[1])
        }
        
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
        
        if(length(legend) != 0){legend = rownames(data)}
        legend("bottomright", legend=legend,
               col=colours, lty=rep("solid",dim(data)[2]), lwd = rep(2.5,dim(data)[2]),cex=0.8)
        
      }
      
      dev.off()
    }
    
    if(toPDF){
      cat("Creating PDF...\n")
      pdf("classResults.pdf")
      
      if(!is.matrix(data)){
        
        plot(data,type='l',col=colours[1], main=main,
             xlab=xlab,ylab=ylab,axes=TRUE,frame.plot=TRUE,lwd = 2.5)
        
        if(xgrid){
          grid(nx = TRUE, ny = NULL, col = "gray", lty = "dashed")
        }
        if(ygrid){
          grid(nx = NULL, ny = TRUE, col = "gray", lty = "dashed")
        }
        
        
        if(legend != ""){
          legend("bottomright", legend=legend, col=colours[1],lwd = 2.5,cex=0.8,lty = "solid")
        }
        
      }else if(is.matrix(data)){
        
        if(length(colours) == dim(data)[1]){
          colours = colours
        }else if(length(colours) != dim(data)[1] && dim(data)[1] <= length(c_palette)){
          colours = c_palette[1:dim(data)[1]]
        }else {
          colours = sample(colors()[grep('gr(a|e)y', colors(), invert = TRUE)],dim(data)[1])
        }
        
        plot(data[1,],type='l',col=colours[2], main=main,
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