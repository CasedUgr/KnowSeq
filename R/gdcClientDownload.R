#' This function downloads a list of controlled files from GDC Portal with the user token and the manifest with the information about the desired controlled files.
#'
#' @param manifestPath Path to the samples manifest
#' @param controlled Parameter that indicates if data to download are controlled or not
#' @param tokenPath Path to the GDC token if data are controlled
#' @return Nothing to return.
#' @examples 
#' # This function needs the download of the pre-compiled tools supplied by KnowSeq.
#' \dontrun{gdcClientDownload("PathToTheToken", "PathToTheFileWithDownloadInfo", dataMatrix)}

gdcClientDownload <- function(manifestPath, controlled = FALSE, tokenPath = "") {

  if(version$os == "linux-gnu"){
    
    cat("\nGNU/Linux OS detected. All the tools and aligners are available for these type of operative systems.\n")
    
    if(dir.exists("unixUtils/")){
      cat("Directory unixUtils found. Checking the tools...\n")
      if(file.exists("unixUtils/salmon/bin/salmon")){cat("Salmon found!\n")}else{stop("Salmon not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/hisat2/hisat2")){cat("Hisat2 found!\n")}else{stop("Hisat2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/bowtie2/bowtie2")){cat("Bowtie2 found!\n")}else{stop("Bowtie2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/samtools/samtools")){cat("Samtools found!\n")}else{stop("Samtools not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/sratoolkit/bin/fastq-dump.2.8.0")){cat("Sratoolkit found!\n")}else{stop("Sratoolkit not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/htseq/scripts-2.7/htseq-count")){cat("Htseq found!\n")}else{stop("Htseq not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/kallisto/kallisto")){cat("Kallisto found!\n")}else{stop("Kallisto not found, please remove unixUtils folder and re-run the function to download it.\n")}
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
  
  if(file.exists(manifestPath)){cat("Manifest found!\n")}else{stop("Manifest not found, please revise the path to the manifest.\n")}
  if(controlled){
    if(file.exists(tokenPath)){cat("Token found!\n")}else{stop("Token not found, please revise the path to the token.\n")}
    system2("unixUtils/gdcClient/gdc-client download", args = paste("-t ",tokenPath," -m ", manifestPath," -n 64 --no-segment-md5sums --no-file-md5sum"))
  }else{
    system2("unixUtils/gdcClient/gdc-client download", args = paste(" -m ", manifestPath," -n 64 --no-segment-md5sums --no-file-md5sum"))
  }


  cat("Moving the downloaded files to ReferenceFiles/Samples/RNAseq/BAMFiles/ \n")

}
