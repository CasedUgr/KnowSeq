#' sraToFastq downloads and converts the sra files to fastq files. The function admits both gz and sra formats.
#'
#' This function downloads and converts the sra files to fastq files by using the URLs indicated through the identifier argument. The function admits both gz and sra formats. This function is used internally by \code{\link{rawAlignment}} but it can be used separatelly.
#' @param identifier A vector that contains a list with the URLs requested.
#' @return Nothing.
#' @examples
#' # This function needs the download of the pre-compiled tools supplied by KnowSeq.
#' 
#' \dontrun{sraToFastq("SRA1")}

sraToFastq <- function(identifier){

  if(!is.character(identifier)){

    stop("The argument has to be a character vector.")

  }else{
    if(version$os == "linux-gnu"){
      
      cat("\nGNU/Linux OS detected. All the tools and aligners are available for these type of operative systems.\n")
      fastqDump = "unixUtils/sratoolkit/bin/fastq-dump"
      prefecth = "unixUtils/sratoolkit/bin/prefetch"
      
      if(dir.exists("unixUtils/")){
        cat("Directory unixUtils found. Checking the tools...\n")
        if(file.exists("unixUtils/hisat2/hisat2")){cat("Hisat2 found!\n")}else{stop("Hisat2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
        if(file.exists("unixUtils/bowtie2/bowtie2")){cat("Bowtie2 found!\n")}else{stop("Bowtie2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
        if(file.exists("unixUtils/samtools/samtools")){cat("Samtools found!\n")}else{stop("Samtools not found, please remove unixUtils folder and re-run the function to download it.\n")}
        if(file.exists("unixUtils/sratoolkit/bin/fastq-dump")){cat("Sratoolkit found!\n")}else{stop("Sratoolkit not found, please remove unixUtils folder and re-run the function to download it.\n")}
        if(file.exists("unixUtils/htseq/scripts-2.7/htseq-count")){cat("Htseq found!\n")}else{stop("Htseq not found, please remove unixUtils folder and re-run the function to download it.\n")}
        if(file.exists("unixUtils/gdcClient/gdc-client")){cat("GDC client found!\n")}else{stop("GDC client not found, please remove unixUtils folder and re-run the function to download it.\n")}
        
      }else{
        
        download.sucess = FALSE
        
        while(!download.sucess){
          
          decission <- readline(prompt="In order to use the aligners, it is necessary to download a pre-compiled version of them. \nThe file has the following aligners and tools: tophat2, hisat2, salmon, bowtie2, samtools, sratoolkit, htseq and gdc-client. \nDo you accept the download? (Y/N): ")
          
          if(decission == 'Y' || decission == 'y'){
            
            cat("Downloading the pre-compiled version of the aligners...\n")
            download.file(url = "http://iwbbio.ugr.es/utils/unixUtils.tar.gz",destfile="unixUtils.tar.gz",method="libcurl")
            cat("Decompressing unixUtils folders with the aligners...\n")
            untar("unixUtils.tar.gz")
            file.remove("rm unixUtils.tar.gz")
            download.sucess = TRUE
            
          } else if(decission == 'N' || decission == 'n'){
            
            stop("Alignment aborted because the pre-compiled aligners are necessaries.")
            
          }
          
        }
        
      }
      
    }else{
      stop("This function is only supported by GNU/Linux distributions due to the external pre-compiled tools are designed for these type of operating system. The version of MAC-OS will be added in next releases.")
    }

    cat(paste("\nDownloading the file ",identifier,".sra...\n\n",sep = ""))

    system2(prefecth, args = identifier)
    
    if(file.exists(paste(identifier,"/",identifier,".sra",sep = ""))){
      
      file = paste(identifier,"/",identifier,".sra",sep = "")
      
    }else{
      
      stop("The download has failed due to a connection error, please check your internet connection or try again.")
      
    }
    
    if(grepl(".gz",file)){

      cat(paste("Decompressing ", file," file...\n",sep = ""))

      gunzip(file)

      file = substr(file,1,nchar(file)-3)

      if(grepl(".sra",file)){

        cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
        system2(fastqDump, args = paste(" --split-e ", file,sep = ""))
        fileMove(file,paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,sep = ""))
        fileMove(paste(identifier,".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",identifier,".fastq",sep = ""))

      }else if(grepl(".fastq",file)){

        cat(paste(file," is already in fastq format\n",sep = ""))
        fileMove(file,paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",file,sep = ""))

      }

    }else if(grepl(".sra",file)){

      cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
      system2(fastqDump, args = paste(" --split-e ", file,sep = ""))
      fileMove(paste(identifier,".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",identifier,".fastq",sep = ""))
      fileMove(file,paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,sep = ""))

    }else if(grepl(".fastq",file)){

      cat(paste(file," is already in fastq format\n",sep = ""))
      fileMove(file,paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",file,sep = ""))


    }else{

      fileMove(file,paste(file,".sra",sep = ""))
      cat(paste(file, " downloaded successfully, converting file to fastq format...\n",sep = ""))
      system2(fastqDump, args = paste(" --split-e ", file,sep = ""))
      fileMove(paste(file,".sra",sep = ""),paste("ReferenceFiles/Samples/RNAseq/SRAFiles/",file,".sra",sep = ""))
      fileMove(paste(identifier,".fastq",sep = ""),paste("ReferenceFiles/Samples/RNAseq/FASTQFiles/",identifier,".fastq",sep = ""))


    }
  }
}
