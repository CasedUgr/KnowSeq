#' rawAlignment allows downloading and processing the fastq samples in a CSV file.
#'
#' This function allows downloading and processing the fastq samples in a CSV file. Also, samples can be aligned by using hisat2. Finally, the function can downloads the reference files required: FASTA Reference Genome and GTF file.
#' @param data The ID of the variable which contains the samples. Our recommendation is to load this variable from a CSV file.
#' @param downloadRef A logical parameter that represents if the reference files will be downloaded or not.
#' @param downloadSamples A logical parameter that represents if the samples of the CSV file will be downloaded or not.
#' @param createIndex A logical parameter that represents if the index of the aligner would be created or not.
#' @param BAMfiles A logical parameter that represents if the you want the BAM files or not.
#' @param SAMfiles A logical parameter that represents if the you want the SAM files or not.
#' @param countFiles A logical parameter that represents if the you want the Count files or not.
#' @param referenceGenome This parameter allows choosing the reference genome that will be used for the alignment. The options are 37,38 or custom. The two first are human genomes, but with the third option you can choose any genome stored in the computer.
#' @param customFA The path to the custom FASTA file of the reference genome.
#' @param customGTF The path to the custom GTF file.
#' @param fromGDC A logical parameter that allows processing BAM files from GDC portal by using the custom reference genome from GDC.
#' @param tokenPath The path to the GDC portal user token. It is required to downloads the controlled BAM files.
#' @param manifestPath The path to the manifest with the information required to downloads the controlled BAM files selected in GDC Portal.
#' @param hisatParameters Parameter that allow to modify the default configuration for the Hisat2 aligner.
#' @return Nothing to return.
#' @examples
#' # Due to the high computational cost, we strongly recommend it to see the offical documentation and the complete example included in this package:
#'
#' dir <- system.file("extdata", package="KnowSeq")
#' 
#' #Using read.csv for NCBI/GEO files (read.csv2 for ArrayExpress files)
#' GSE74251csv <- read.csv(paste(dir,"/GSE74251.csv",sep = ""))
#' 
#' \dontrun{rawAlignment(GSE74251csv,downloadRef=FALSE,downloadSamples=FALSE, createIndex = TRUE, BAMfiles = TRUE, SAMfiles = TRUE, countFiles = TRUE, referenceGenome = 38, customFA = "", customGTF = "", fromGDC = FALSE, tokenPath = "", manifestPath = "")}


rawAlignment <- function(data,downloadRef=FALSE,downloadSamples=FALSE, createIndex = TRUE, BAMfiles = TRUE, SAMfiles = TRUE, countFiles = TRUE, referenceGenome = 38, customFA = "", customGTF = "", fromGDC = FALSE, tokenPath = "", manifestPath = "", hisatParameters="-p 8 --dta-cufflinks"){

  if(version$os == "linux-gnu"){

    cat("\nGNU/Linux OS detected. All the tools are available for these type of operative systems.\n")

    if(dir.exists("unixUtils/")){
      cat("Directory unixUtils found. Checking the tools...\n")
      if(file.exists("unixUtils/hisat2/hisat2")){cat("Hisat2 found!\n")}else{stop("Hisat2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/bowtie2/bowtie2")){cat("Bowtie2 found!\n")}else{stop("Bowtie2 not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/samtools/samtools")){cat("Samtools found!\n")}else{stop("Samtools not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/sratoolkit/bin/fastq-dump")){cat("Sratoolkit found!\n")}else{stop("Sratoolkit not found, please remove unixUtils folder and re-run the function to download it.\n")}
      if(file.exists("unixUtils/gdcClient/gdc-client")){cat("GDC client found!\n")}else{stop("GDC client not found, please remove unixUtils folder and re-run the function to download it.\n")}

    }else{

      download.sucess = FALSE

      while(!download.sucess){

        decission <- readline(prompt="In order to use the aligners, it is necessary to download a pre-compiled version of them. \nThe file has the following aligners and tools: hisat2, salmon, bowtie2, samtools, sratoolkit, htseq and gdc-client. \nDo you accept the download? (Y/N): ")

        if(decission == 'Y' || decission == 'y'){

          cat("Downloading the pre-compiled version of the tools...\n")
          download.file(url = "http://iwbbio.ugr.es/utils/unixUtils.tar.gz",destfile="unixUtils.tar.gz",method="libcurl")
          cat("Decompressing unixUtils.tar.gz...\n")
          untar("unixUtils.tar.gz")
          system2("rm", args="unixUtils.tar.gz")
          download.sucess = TRUE

        } else if(decission == 'N' || decission == 'n'){

          stop("Alignment aborted because the pre-compiled tools are necessaries.")

        }

      }

    }

  }else{
    stop("This function is only supported by GNU/Linux distributions due to the external pre-compiled tools are designed for these type of operating system. The version of MAC-OS will be added in next releases.")
  }

  if(!is.data.frame(data)){

    stop("Please, use a dataframe for the data parameter.")

  }else if(!is.logical(downloadRef)){

    stop("Please, use a TRUE or FALSE for the downloadRef parameter.")

  }else if(!is.logical(downloadSamples)){

    stop("Please, use a TRUE or FALSE for the downloadSamples parameter.")

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

    if(file.exists("/usr/bin/python3")){
      cat("Python3 found!\n")
      if(file.exists("/usr/bin/pip3")){
        
        cat("Pip3 found!\n")
        output = system2("python3", args = "-c 'import HTSeq;'")
        
        if(output == 0){
          cat("HTSeq found!\n")
        }else{
          stop("HTSeq is required to calculate the counts. To install it just the next command:\n sudo pip3 install HTSeq\n")
        }
      }else{
        stop("Pip3 is required to install HTSeq. To install it just the next command:\n sudo apt install python3-pip\n")
      }
    }else{
      stop("Python3 is required to carry out this step. Please, install it and rerun the function.\n")
    }
    
    if(fromGDC){

      if(downloadRef){

        cat ("Downloading reference genome of GDC from url https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f\n")
        download.file("https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f",destfile="gencode.v22.annotation.gtf.gz",method="libcurl")
        cat ("Decompressing reference genome...\n")
        gunzip("gencode.v22.annotation.gtf.gz")

      }

      if(downloadSamples){

        if(file.exists(tokenPath)){cat("Token found!\n")}else{stop("Token not found, please revise the path to the token.\n")}
        if(file.exists(manifestPath)){cat("Manifest found!\n")}else{stop("Manifest not found, please revise the path to the manifest.\n")}

        cat ("Downloading samples from GDC by using the selected token and manifest...\n")
        gdcClientDownload(tokenPath = tokenPath, manifestPath = manifestPath)

      }
      gf = "gencode.v22.annotation.gtf"

      for(i in seq_len(dim(data)[1])){

        cat("Converting to count file...\n")

        filePathBam = paste("ReferenceFiles/Samples/RNAseq/BAMFiles", data$File.ID[i], data$File.Name[i], sep = "/")
        filePathSam = paste("ReferenceFiles/Samples/RNAseq/SAMFiles/", data$Sample.ID[i],".sam", sep = "")
        params = paste(filePathSam, filePathBam)

        mkdirs(paste("ReferenceFiles/Samples/RNAseq/CountFiles/", data$Sample.ID[i],sep = ""))
        countFile = paste("ReferenceFiles/Samples/RNAseq/CountFiles/", data$Sample.ID[i],"/", data$Sample.ID[i], ".count", sep = "")
        gtf = paste(gf,"  > ", countFile)

        filePathBamSorted = paste("ReferenceFiles/Samples/RNAseq/BAMFiles/", data$Sample.ID[i],"Sorted.bam", sep = "")
        system2("unixUtils/samtools/samtools", args = paste("sort -n ", filePathBam , " -o ",filePathBamSorted,sep = ""))

        system2("unixUtils/samtools/samtools", args = paste("view -f 0x0002 ", filePathBamSorted ," |  awk '!/\t\\*\t/' - | htseq-count -s no -a 10 - ", gtf, sep=""))

      }

    }else{
        hisatAlignment(data,downloadRef=downloadRef,downloadSamples=downloadSamples,createIndex=createIndex,BAMfiles=BAMfiles,SAMfiles=SAMfiles,countFiles=countFiles,referenceGenome=referenceGenome,customFA = customFA,customGTF = customGTF, hisatParameters = hisatParameters)
    }
  }
}
