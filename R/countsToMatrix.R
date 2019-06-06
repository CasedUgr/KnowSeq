#' countsToMatrix merges in a matrix the information in the count files.
#'
#' The function merges in a matrix the information in the count files. It can be used from 1 to N count files. These count files can be created by using the function \code{\link{rawAlignment}} with the raw files of RNA-seq.
#' @param csvFile The csv that contains the name and the path to each of the count files. The column of the name of the file must be named Run and the column that contains the paths must be named Path. Furthermore, to facilitate the posterior steps, a column named Class that contains the classes for the samples must be required.
#' @param sep The separator character of the csvFile or tsvFile.
#' @return A matrix with the ensembl ID in the rows and all the samples of each count files in the columns.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' countsInfo <- read.csv(paste(dir,"/countFiles/mergedCountsInfo.csv",sep = ""))
#' 
#' countsInfo$Path <- paste(dir,"/countFiles/",countsInfo$Run,sep = "")
#' 
#' write.csv(countsInfo, file = "countsInfo.csv")
#'
#' countsInformation <- countsToMatrix("countsInfo.csv")
#'
#' countsMatrix <- countsInformation$countsMatrix
#' labels <- countsInformation$labels
#'
#' file.remove("countsInfo.csv")

countsToMatrix <- function(csvFile,sep=','){

    if(!file.exists(csvFile)){
      stop("Unable to find the CSV file. Please check the path to the file.")
    }

    countsData <- read.csv(csvFile, stringsAsFactors = FALSE, sep = sep)

    if(sum(c("Run","Path","Class") %in% colnames(countsData), na.rm = TRUE) != 3){
        stop("The CSV file has to be the following three columns: Run, Path,Class.")
    }

    countf <- vector(mode="character", length=0)
    for(i in seq_len(dim(countsData)[1])){
      dir = paste(countsData$Path[i], countsData$Run[i], sep = "/")
      countf[i] = paste(dir, "count", sep = ".")
      cat(paste("\n",countf[i],sep = ""))
    }
    
    cat(paste("\nMerging ",dim(countsData)[1]," counts files...\n",sep = ""))

    counts = readDGE(countf)$counts
    colnames(counts) = countsData$Run
    noint  =  rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
    cpms  =  cpm(counts)
    keep  =  rowSums(cpms  > 1) >= 1 & !noint
    counts  =  counts[keep,]
    colnames(counts) <- countsData$Run

    ensemblId <- vector(mode="character", length=0)

    for(i in seq_len(length(rownames(counts)))){

      aux <- unlist(strsplit(rownames(counts)[i],"[.]"))[1]
      ensemblId <- c(ensemblId, aux)

    }

    ensemblCounts <- counts
    rownames(ensemblCounts) <- ensemblId

    results <- list(ensemblCounts,countsData$Class)
    names(results) <- c("countsMatrix","labels")

    return(results)

}
