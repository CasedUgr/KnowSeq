#' Download automatically samples from NCBI/GEO and ArrayExpress public databases.
#'
#' Download automatically samples from series of either microarray and RNA-seq. Furthermore, both NCBI/GEO and ArrayExpress public databases are supported. In the case of Microarray, the raw file are downloaded, if they are available, but for RNA-seq a csv is created with the necessary information to download the samples with the function \code{\link{rawAlignment}}.
#'
#' @param samplesVector A vector which contains the different IDs of the wanted series. These IDs are the IDs of the series from NCBI/GEO or ArrayExpress.
#' @return Nothing to return.
#' @examples
#' downloadPublicSeries(c("GSE74251"))

downloadPublicSeries <- function(samplesVector) {

  if(!is.character(samplesVector)){

    stop("The parameter passed is wrong. Please, check that the type of the vector be character.")

  }else if(length(samplesVector)[1] == 0 || is.null(samplesVector)){

    stop("The samplesVector parameter is empty! Please, provide a right geneList.")

  }else{

    if(dir.exists("ReferenceFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/")}
    if(dir.exists("ReferenceFiles/Samples/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/BAMFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/BAMFiles/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/SAMFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/SAMFiles/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/CountFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/CountFiles/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/SRAFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/SRAFiles/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/FASTQFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/FASTQFiles/")}
    if(dir.exists("ReferenceFiles/Samples/RNAseq/QuantFiles/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/RNAseq/QuantFiles/")}

    for(serie in samplesVector){
      if(grepl("GSE",serie)){

        cat("\n******************************** NCBI/GEO format detected ********************************\n\n")

        webpage <- GET(paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",serie,sep = ""))
        webpage_content <- content(webpage, as = "text", encoding = "ISO-8859-1")
        webpage <- readLines(tc <- textConnection(webpage_content)); close(tc)
        pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)

        # parse the tree by tables
        x <- xpathSApply(pagetree, "//*/table", xmlValue)
        # do some clean up with regular expressions
        x <- unlist(strsplit(x, "\n"))
        x <- gsub("\t","",x)
        x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
        x <- x[!(x %in% c("", "|"))]

        if(any(grepl("Expression profiling by array", x))){

          cat("******************************** MICROARRAY series detected *******************************\n\n")
          cat(paste("Searching the RAW data from the microarray series ",serie,"\n",sep = ""))
          if(url.exists(paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substring(serie,1,nchar(serie)-3),"nnn/",serie,"/suppl/",serie,"_RAW.tar",sep = ""))){

            cat("The RAW data have been found correctly...\n")
            POST(paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",substring(serie,1,nchar(serie)-3),"nnn/",serie,"/suppl/",serie,"_RAW.tar",sep = ""), write_disk(paste(serie,"_RAW.tar",sep = ""), overwrite=TRUE))
            cat(paste("Decompressing ", serie,"_RAW.tar...\n",sep=""))

            if(dir.exists("ReferenceFiles/Samples")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/")}
            if(dir.exists("ReferenceFiles/Samples/Microarray/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/Microarray/")}

            fileMove(from = paste(serie,"_RAW.tar",sep = ""),
                      to =  paste("ReferenceFiles/Samples/Microarray/",serie,"/",serie,"_RAW.tar",sep = ""))

            if(dir.exists(paste("ReferenceFiles/Samples/Microarray/",serie,sep=""))){}else{ system2("mkdir", args = paste("ReferenceFiles/Samples/Microarray/",serie,sep=""))}

            system2("tar", args = paste("-xf ReferenceFiles/Samples/Microarray/", serie,"/",serie,"_RAW.tar -C ", "ReferenceFiles/Samples/Microarray/",serie,sep = ""))
            system2("rm", args = paste("ReferenceFiles/Samples/Microarray/", serie,"/",serie,"_RAW.tar",sep = ""))

            CELFiles <- list.files(paste("ReferenceFiles/Samples/Microarray/",serie,"/",sep=""))



            for(CEL in CELFiles){

              cat(paste("Decompressing ", CEL, " file...\n",sep=""))

              system2("gzip", args = paste("-d ReferenceFiles/Samples/Microarray/",serie,"/",CEL,sep = ""))

            }



          }else{

            cat("There aren't RAW data available for this series or the series doesn't exist\n")

          }



        }else if(any(grepl("Expression profiling by high throughput sequencing", x))){

          cat("******************************** RNA-SEQ series detected ********************************\n\n")
          cat(paste("Building the CSV file for RNA-seq samples for series ",serie,"\n",sep = ""))

          gsmList <- unique(x[grep("\\bGSM.*\\d", x)])

          Run = character()
          download_path = character()
          LibraryLayout = character()
          i = 1

          for(gsm in gsmList){

            webpage <- GET(paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gsm,sep = ""))
            webpage_content <- content(webpage, as = "text", encoding = "ISO-8859-1")
            webpage <- readLines(tc <- textConnection(webpage_content)); close(tc)
            pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)

            # parse the tree by tables
            xx <- xpathSApply(pagetree, "//*/table", xmlValue)
            # do some clean up with regular expressions
            xx <- unlist(strsplit(xx, "\n"))
            xx <- gsub("\t","",xx)
            xx <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", xx, perl=TRUE)
            xx <- xx[!(xx %in% c("", "|"))]

            srxID <- unique(unlist(xx[grep("\\bSRX.*\\d", xx)])[1])

            webpage <- GET(paste("https://www.ncbi.nlm.nih.gov/sra?term=",srxID,sep = ""))
            webpage_content <- content(webpage, as = "text", encoding = "ISO-8859-1")
            webpage.treated <- readLines(tc <- textConnection(webpage_content)); close(tc)
            pagetree <- htmlTreeParse(webpage.treated, error=function(...){}, useInternalNodes = TRUE)

            # parse the tree by tables
            xxx <- xpathSApply(pagetree, "//*/table", xmlValue)
            # do some clean up with regular expressions
            xxx <- unlist(strsplit(xxx, "\n"))
            xxx <- gsub("\t","",xxx)
            xxx <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", xxx, perl=TRUE)
            xxx <- xxx[!(xxx %in% c("", "|"))]

            LibraryLayoutElement = ""
            RunElement = ""

            if(grepl("<span>SINGLE</span>",webpage)){
              LibraryLayoutElement = "SINGLE"
            }else if(grepl("<span>PAIRED</span>",webpage)){

              LibraryLayoutElement = "PAIRED"

            }

            m <- gregexpr("\\brun=SRR(.*?)\"", webpage_content)
            regx <- regmatches(webpage_content,m)
            RunElement <- as.character(unlist(regx))

            m <- gregexpr("\\bSRP/.*\\d", x)
            regx <- regmatches(x,m)

            if(length(RunElement)==1){

              RunElement <- substr(RunElement,5,nchar(RunElement))
              RunElement <- substr(RunElement,1,nchar(RunElement)-1)


              Run[i] <- RunElement
              downloadpath = paste("ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr(RunElement,1,6),"/",RunElement,"/",RunElement,".sra",sep = "")
              download_path[i] <- downloadpath
              LibraryLayout[i] <- LibraryLayoutElement

              i = i + 1

            }else{


                RunElement[1] <- substr(RunElement[1],5,nchar(RunElement[1]))
                RunElement[1] <- substr(RunElement[1],1,nchar(RunElement[1])-1)

                Run[i] <- RunElement[1]
                downloadpath =  paste("ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr(RunElement[1],1,6),"/",RunElement[1],"/",RunElement[1],".sra",sep = "")
                download_path[i] <- downloadpath
                LibraryLayout[i] <- LibraryLayoutElement

                RunElement[2] <- substr(RunElement[2],5,nchar(RunElement[2]))
                RunElement[2] <- substr(RunElement[2],1,nchar(RunElement[2])-1)

                Run[i+1] <- RunElement[2]
                downloadpath =  paste("ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr(RunElement[2],1,6),"/",RunElement[2],"/",RunElement[2],".sra",sep = "")
                download_path[i+1] <- downloadpath
                LibraryLayout[i+1] <- LibraryLayoutElement


                i = i + 2

            }

          }

          rnaseq.data <- data.frame(Run, LibraryLayout, download_path)
          cat("Exporting CSV to ReferenceFiles folder...")
          write.csv(rnaseq.data,paste("ReferenceFiles/",serie,".csv",sep = ""))

        }


      }else if(grepl("E-",serie)){

        cat("\n******************************** ArrayExpress format detected ********************************\n\n")

        webpage <- GET(paste("https://www.ebi.ac.uk/arrayexpress/experiments/",serie,"/",sep = ""))
        webpage_content <- content(webpage, as = "text", encoding = "ISO-8859-1")
        webpage <- readLines(tc <- textConnection(webpage_content)); close(tc)
        pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)

        # parse the tree by tables
        x <- xpathSApply(pagetree, "//*/table", xmlValue)
        # do some clean up with regular expressions
        x <- unlist(strsplit(x, "\n"))
        x <- gsub("\t","",x)
        x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
        x <- x[!(x %in% c("", "|"))]

        if(any(grepl("transcription profiling by array", x))){

          cat("******************************** MICROARRAY series detected *******************************\n\n")

          cat(paste("Searching the RAW data from the microarray series ",serie,"\n",sep = ""))

          if(any(grepl("\\braw\\b",x))){

            cat("The RAW data have been found correctly...\n")

            m <- gregexpr("\\bE-MTAB(.*?)zip", x)
            regx <- regmatches(x,m)
            regx <- unlist(regx)

            regx <- unique(regx)

            for(raw in regx){
              cat(paste("https://www.ebi.ac.uk/arrayexpress/files/",serie,"/",raw,"\n",sep = ""))
              POST(paste("https://www.ebi.ac.uk/arrayexpress/files/",serie,"/",raw,sep = ""), write_disk(raw), overwrite=TRUE)
              
              cat(paste("Decompressing ",raw,"\n",sep=""))

              if(dir.exists("ReferenceFiles/Samples")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/")}
              if(dir.exists("ReferenceFiles/Samples/Microarray/")){}else{ system2("mkdir", args = "ReferenceFiles/Samples/Microarray/")}

              if(dir.exists(paste("ReferenceFiles/Samples/Microarray/",serie,sep=""))){}else{ system2("mkdir", args = paste("ReferenceFiles/Samples/Microarray/",serie,sep=""))}

              fileMove(from = paste(raw,sep = ""),
                       to =  paste("ReferenceFiles/Samples/Microarray/",serie,"/",raw,sep = ""))
              
            }

          }else{

            cat("There aren't RAW data available for this series or the series doesn't exist\n")

          }


        }else if(any(grepl("RNA-seq of coding RNA", x))){

          cat("******************************** RNA-SEQ series detected ********************************\n\n")

          if(grepl(paste(serie,".sdrf.txt",sep = ""),x[1])){

            cat("Downloading the sdrf.txt and converting to TSV...\n")

            POST(paste("https://www.ebi.ac.uk/arrayexpress/files/",serie,"/",serie,".sdrf.txt",sep = ""), write_disk(paste(serie,".sdrf.txt",sep = "")), overwrite=TRUE)
            
            fileMove(from = paste(serie,".sdrf.txt",sep = ""),
                      to =  paste("ReferenceFiles/",serie,".tsv",sep = ""))

          }else{

            cat("There isn't any sdrf file available for this series or the series doesn't exist\n")

          }

        }

      }else{

        stop("The IDs of the series doesn't exist or it's not supported. Please, try again with NCBI/GEO and ArrayExpress IDs\n")

      }
    }
  }
}
