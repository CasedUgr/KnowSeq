#' sraToFastq downloads and converts the sra files to fastq files. The function admits both gz and sra formats.
#'
#' This function downloads and converts the sra files to fastq files by using the URLs indicated through the urlsVector argument. The function admits both gz and sra formats. This function is used internally by \code{\link{rawAlignment}} but it can be used separatelly.
#' @param urlsVector A vector that contains a list with the URLs requested.
#' @return Nothing.

sraToFastq <- function(urlsVector){

  if(!is.character(urlsVector)){

    stop("The argument has to be a character vector. Please format your URLs variable")

  }else{
    if(version$os == "linux-gnu"){

      fastqDump = "unixUtils/sratoolkit/bin/fastq-dump.2.8.0"

    }else{
      stop("This function is only supported by GNU/Linux distributions due to the external pre-compiled tools are designed for these type of operating system. The version of MAC-OS will be added in next releases.")

    }


    splittedUrl <- strsplit(urlsVector,"/")
    splittedUrl <- unlist(splittedUrl)
    file <- splittedUrl[length(splittedUrl)]
    cat(paste("\nDownloading the file ",file," from url ", urlsVector,"...\n\n",sep = ""))

    if(url.exists(urlsVector)){

      f = CFILE(file, mode = "wb")
      curlPerform(url = urlsVector, writedata = f@ref, noprogress=FALSE)
      close(f)

    }else{stop("The URL doesn't exist or is unreachable.")}

    if(grepl(".gz",file)){

      cat(paste("Decompressing ", file," file...\n",sep = ""))

      gunzip(file)

      file = substr(file,1,nchar(file)-3)

      if(grepl(".sra",file)){

        cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
        system(paste(fastqDump," --split-3 ", file,sep = ""))
        fileMove(file,paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,sep = ""))
        fileMove(paste(substr(file,1,nchar(file) - 4),".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",substr(file,1,nchar(file) - 4),".fastq",sep = ""))

      }else if(grepl(".fastq",file)){

        cat(paste(file," is already in fastq format\n",sep = ""))
        fileMove(file,paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",file,sep = ""))

      }

    }else if(grepl(".sra",file)){

      cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
      system(paste(fastqDump," --split-3 ", file,sep = ""))
      fileMove(file,paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,sep = ""))
      fileMove(paste(substr(file,1,nchar(file) - 4),".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",substr(file,1,nchar(file) - 4),".fastq",sep = ""))

    }else if(grepl(".fastq",file)){

      cat(paste(file," is already in fastq format\n",sep = ""))
      fileMove(file,paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",file,sep = ""))


    }else{

      fileMove(file,paste(file,".sra",sep = ""))
      cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
      system(paste(fastqDump," --split-3 ", file,sep = ""))
      fileMove(paste(file,".sra",sep = ""),paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,".sra",sep = ""))
      fileMove(paste(file,".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",file,".fastq",sep = ""))


    }
  }
}
