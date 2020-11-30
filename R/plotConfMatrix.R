#' plotConfMatrix plots a confusion matrix with some statistics.
#'
#' The function plots a confusion matrix with some statistics. The function is used internally by \code{\link{dataPlot}} but it can be used separatelly.
#' @param data A table which contains a confusion matrix.
#' @return Nothing to return.
#' @examples
#' data <- table(as.factor(c(1,2,4,2,4,5)),as.factor(c(1,2,5,4,5,2)))
#' plotConfMatrix(data)

plotConfMatrix <- function(data){

  res <- confusionMatrix(data)
  
  layout(matrix(c(1,1,2)))

  # The above rescales the confusion matrix such that columns sum to 100.
  opar <- par(mar=c(6.1, 9.1, 1, 2))
  x <- x.orig <- unclass(data)
  x <- log(x + 0.5) * 2.33
  x[x < 0] <- NA
  x[x > 10] <- 10
  diag(x) <- -diag(x)
  image(seq_len(ncol(x)), seq_len(ncol(x)),
        -(x[, nrow(x):1]), xlab='', ylab='',
        col=colorRampPalette(c(hsv(h = 0, s = 0.9, v = 0.9, alpha = 1),
                               hsv(h = 0, s = 0, v = 0.9, alpha = 1),
                               hsv(h = 2/6, s = 0.9, v = 0.9, alpha = 1)))(41),
        xaxt='n', yaxt='n', zlim=c(-10, 10))
  axis(1, at=seq_len(ncol(x)), labels=FALSE, cex.axis=1.2, font = 2)
  title(xlab='Prediction', line=4.5, cex = 1.2)
  text(seq_len(ncol(x)), par("usr")[3] - 0.1, labels = colnames(x), pos = 1 ,font = 2, xpd = TRUE)
  
  axis(2, at=ncol(x):1, labels=colnames(x), las=1, cex.axis=1.2,font = 2)
  title(ylab='Reference', line=7.5, cex = 1.2)
  abline(h = 0:ncol(x) + 0.5, col = 'gray')
  abline(v = 0:ncol(x) + 0.5, col = 'gray')
  text(seq_len(ncol(x)), rep(ncol(x):1, each=ncol(x)),
       labels = c(x.orig),cex=1.2, font=2)
  box(lwd=2)
  par(opar) # reset par

  # add in the specifics
  plot(c(100, 0),c(0, 50),  type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(20, 35, names(res$overall[1]), cex=1.2, font=2)
  text(20, 15, round(as.numeric(res$overall[1])*100, 3), cex=1.4)

  if(ncol(x) > 2){

      text(40, 35, colnames(res$byClass)[7], cex=1.2, font=2)
      res$byClass[is.na(res$byClass[,7]),7] <- 0
      text(40, 15, round(mean(as.numeric(res$byClass[,7]))*100, 3), cex=1.4)
      text(60, 35, colnames(res$byClass)[1], cex=1.2, font=2)
      text(60, 15, round(mean(as.numeric(res$byClass[,1]))*100, 3), cex=1.4)
      text(80, 35, colnames(res$byClass)[2], cex=1.2, font=2)
      text(80, 15, round(mean(as.numeric(res$byClass[,2]))*100, 3), cex=1.4)

  }else{

    text(50, 35, names(res$byClass)[1], cex=1.2, font=2)
    text(50, 15, round(mean(as.numeric(res$byClass[1]))*100, 3), cex=1.4)
    text(80, 35, names(res$byClass)[2], cex=1.2, font=2)
    text(80, 15, round(mean(as.numeric(res$byClass[2]))*100, 3), cex=1.4)

  }

}
