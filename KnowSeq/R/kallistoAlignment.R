#' kallistoAlignment allows downloading and processing the fastq samples in a CSV file by using kallisto aligner.
#'
#' This function allows downloading and processing the fastq samples in a CSV file by using kallisto aligner. This funtion is used internally by \code{\link{rawAlignment}} but it can be used separatelly. Furthermore, the function can download the reference files required: FASTA Reference Genome and GTF file.
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
#' @param tx2Counts A matrix with two columns that contains the conversion of transcripts IDs to genes IDs. There is more information in the function \code{\link{tximport}}.
#' @return Nothing to return.
#' @examples
#' # Due to the high computational cost, we strongly recommend it to see the offical documentation and the complete example included in this package:
#'
#' dir <- system.file("extdata", package="KnowSeq")
#' 
#' #Using read.csv for NCBI/GEO files (read.csv2 for ArrayExpress files)
#' GSE74251csv <- read.csv(paste(dir,"/GSE74251.csv",sep = ""))
#' 
#' \dontrun{kallistoAlignment(GSE74251csv,downloadRef=FALSE,downloadSamples=FALSE, createIndex = TRUE, BAMfiles = TRUE, SAMfiles = TRUE, countFiles = TRUE, referenceGenome = 38, customFA = "", customGTF = "",tx2Counts = tx2Counts)}

kallistoAlignment <- function(data,downloadRef=FALSE,downloadSamples=FALSE, createIndex = TRUE, BAMfiles = TRUE, SAMfiles = TRUE, countFiles = TRUE,referenceGenome = 38, customFA = "", customGTF = "",tx2Counts = tx2Counts){

  if(version$os != "linux-gnu"){stop("This function is only supported by GNU/Linux distributions due to the external pre-compiled tools are designed for these type of operating system. The version of MAC-OS will be added in the future if all the tools are available for it.")}

  FastqPath = "ReferenceFiles/Samples/RNAseq/FASTQFiles/"

  cat("Using kallisto sequencing tool...\n")


  if(referenceGenome == 38){

    cat("Using Reference Genome 38\n")
    faurl = "ftp://ftp.ensembl.org/pub/release-90//fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    fatargz = "Homo_sapiens.GRCh38.cdna.toplevel.fa.gz"
    fa = "Homo_sapiens.GRCh38.cdna.toplevel.fa"
    gtfurl = "ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
    gtftargz = "Homo_sapiens.GRCh38.90.gtf.gz"
    gtf = "Homo_sapiens.GRCh38.90.gtf"

    genomeIndexCommand = paste("index ReferenceFiles/", fa, " -i ReferenceFiles/kallisto_Homo_sapiens.index", sep = "")

  }else if(referenceGenome == 37){

    cat("Using Reference Genome 37\n")
    faurl = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz"
    fatargz = "Homo_sapiens.GRCh37.75.cdna.toplevel.fa.gz"
    fa = "Homo_sapiens.GRCh37.75.cdna.toplevel.fa"
    gtfurl = "ftp://ftp.ensembl.org/pub/release-75//gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
    gtftargz = "Homo_sapiens.GRCh37.75.gtf.gz"
    gtf = "Homo_sapiens.GRCh37.75.gtf"

    genomeIndexCommand = paste("index ReferenceFiles/", fa, " -i ReferenceFiles/kallisto_Homo_sapiens.index", sep = "")

  }else if(referenceGenome == "custom"){

    cat("Using custom reference genome and GTF\n")
    if(file.exists(customFA)){cat("Custom FA file found!\n")}else{stop("Custom FA file not found, please check the path to the file.\n")}
    if(file.exists(customGTF)){cat("Custom GTF file found!\n")}else{stop("Custom GTF file not found, please check the path to the file.\n")}

    fa = customFA
    gtf = customGTF
    downloadRef = FALSE

    genomeIndexCommand = paste("index ReferenceFiles/", fa, " -i ReferenceFiles/kallisto_Homo_sapiens.index", sep = "")
    gtf <- gtf
    bowind <- fa

  }else{
    stop("The version of the reference genome is wrong. Please use 37, 38 or custom.")
  }


  if(!is.null(data$Run)){

    cat("###################  GEO FORMAT DETECTED  ###################\n\n")
    samples = data[,c("Run","LibraryLayout")]


    for(i in seq_len(nrow(samples))){
      if(samples[i,]$LibraryLayout == "PAIRED"){
        samples$fastq1[i] = paste0(FastqPath,data$Run[i],"_1.fastq", collapse = ",")
        samples$fastq2[i] = paste0(FastqPath,data$Run[i],"_2.fastq", collapse = ",")
      }else{
        samples$fastq1[i] = paste0(FastqPath,data$Run[i],".fastq", collapse = ",")
        samples$fastq2[i] = ""
      }
    }


    if(downloadSamples){

      cat("Downloading the requested samples...\n")
      urls <- as.character(data$download_path)
      lapply(urls,sraToFastq)
      if(length(list.files(pattern = "*.fastq")) != 0){
        system2("mv", args = "*.fastq ReferenceFiles/Samples/RNAseq/FASTQFiles/")

      }

    }

    if(downloadRef){

      cat(paste("Downloading reference human genome from url ", faurl,"\n",sep = ""))
      f = CFILE(fatargz, mode = "wb")
      curlPerform(url = faurl, writedata = f@ref, noprogress=FALSE)
      cat ("Decompressing reference genome...\n")
      gunzip(fatargz, overwrite = TRUE)
      cat("Downloading GTF file from url\n")
      download.file(url = gtfurl,destfile=gtftargz,method="libcurl")
      cat ("Decompressing GTF file genome...\n")
      gunzip(gtftargz, overwrite = TRUE)

      fileMove(from = fa,
                to = paste("ReferenceFiles/", fa, sep = ""))
      fileMove(from = gtf,
               to = paste("ReferenceFiles/", gtf, sep = ""))

    }


    if(createIndex){

      cat("Building index file for kallisto...\n")
      system2("unixUtils/kallisto/kallisto", args = genomeIndexCommand)

    }

    indexName = "ReferenceFiles/kallisto_Homo_sapiens.index"
    gtf = paste("ReferenceFiles/", gtf, sep = "")

    for(i in seq_len(dim(samples)[1])){

      cat("Calculating genes count...\n")

      fastq = paste(samples[i,]$fastq1, samples[i,]$fastq2)
      filePath = paste("ReferenceFiles/Samples/RNAseq/QuantFiles/", samples[i,]$Run, sep = "")
      if(samples[i,]$LibraryLayout == "PAIRED"){
        kallistoCommand = paste("quant -i ", indexName, " --gtf ", gtf," -o ", filePath," -b 100 ", samples[i,]$fastq1 , " ", samples[i,]$fastq2)
      }else{
        kallistoCommand = paste("quant -i ", indexName, " --gtf ", gtf," -o ", filePath," -b 100 --single -l 180 -s 20 ", samples[i,]$fastq1)
      }
      system2("unixUtils/kallisto/kallisto", args = kallistoCommand)

      txi <- tximport(paste("ReferenceFiles/Samples/RNAseq/QuantFiles/", samples[i,]$Run,"/abundance.h5",sep = ""), type = "kallisto", tx2gene = tx2Counts)

      mkdirs(paste("ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Run,"/",sep = ""))

      write.table(as.matrix(data.frame(rownames(txi$counts),txi$counts)),file = paste("ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Run,"/",samples[i,]$Run, ".count",sep = ""), row.names=FALSE, col.names=FALSE,sep = "\t")

      cat(paste("File succesfully stored at ",paste("ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Run,"/",samples[i,]$Run, ".count\n",sep = "")))

    }



  }else if(!is.null(data$Comment.ENA_RUN.)){

    cat("###################  ARRAY EXPRESS FORMAT DETECTED  ###################\n\n")

    samples = data[,c("Comment.ENA_RUN.","Comment.LIBRARY_LAYOUT.")]


    for(i in seq_len(nrow(samples))){
      if(samples[i,]$Comment.LIBRARY_LAYOUT. == "PAIRED"){
        samples$fastq1[i] = paste0(FastqPath,data$Comment.ENA_RUN.[i],"_1.fastq", collapse = ",")
        samples$fastq2[i] = paste0(FastqPath,data$Comment.ENA_RUN.[i],"_2.fastq", collapse = ",")
      }else{
        samples$fastq1[i] = paste0(FastqPath,data$Comment.ENA_RUN.[i],".fastq", collapse = ",")
        samples$fastq2[i] = ""
      }
    }

    if(downloadSamples){

      cat("Downloading the requested samples...\n")
      urls <- as.character(levels(data$Comment.FASTQ_URI.)[as.integer(data$Comment.FASTQ_URI.)])
      lapply(urls,sraToFastq)
      if(length(list.files(pattern = "*.fastq")) != 0){
        system2("mv", args = "*.fastq ReferenceFiles/Samples/RNAseq/FASTQFiles/")

      }

    }

    if(downloadRef){

      cat(paste("Downloading reference human genome from url ", faurl,"\n",sep = ""))
      f = CFILE(fatargz, mode = "wb")
      curlPerform(url = faurl, writedata = f@ref, noprogress=FALSE)
      cat ("Decompressing reference genome...\n")
      gunzip(fatargz, overwrite = TRUE)
      cat("Downloading GTF file from url\n")
      download.file(url = gtfurl,destfile=gtftargz,method="libcurl")
      cat ("Decompressing GTF file genome...\n")
      gunzip(gtftargz, overwrite = TRUE)

      fileMove(from = fa,
               to = paste("ReferenceFiles/", fa, sep = ""))
      fileMove(from = gtf,
               to = paste("ReferenceFiles/", gtf, sep = ""))

    }


    if(createIndex){

      cat("Building index file for kallisto...\n")
      system2("unixUtils/kallisto/kallisto", args = genomeIndexCommand)

    }

    indexName = "ReferenceFiles/kallisto_Homo_sapiens.index"
    gtf = paste("ReferenceFiles/", gtf, sep = "")

    for(i in seq_len(dim(samples)[1])){

      cat("Calculating genes count...\n")

      fastq = paste(samples[i,]$fastq1, samples[i,]$fastq2)
      filePath = paste("ReferenceFiles/Samples/RNAseq/QuantFiles/", samples[i,]$Comment.ENA_RUN., sep = "")
      if(samples[i,]$Comment.LIBRARY_LAYOUT. == "PAIRED"){
        kallistoCommand = paste("quant -i ", indexName, " --gtf ", gtf," -o ", filePath," -b 100 ", samples[i,]$fastq1 , " ", samples[i,]$fastq2)
      }else{
        kallistoCommand = paste("quant -i ", indexName, " --gtf ", gtf," -o ", filePath," -b 100 --single -l 180 -s 20 ", samples[i,]$fastq1)
      }
      system2("unixUtils/kallisto/kallisto", args = kallistoCommand)

      txi <- tximport(paste("ReferenceFiles/Samples/RNAseq/QuantFiles/", samples[i,]$Comment.ENA_RUN.,"/abundance.h5",sep = ""), type = "kallisto", tx2gene = tx2Counts)

      mkdirs(paste("ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Comment.ENA_RUN.,"/",sep = ""))

      write.table(as.matrix(data.frame(rownames(txi$counts),txi$counts)),file = paste("ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Comment.ENA_RUN.,"/",samples[i,]$Comment.ENA_RUN., ".count",sep = ""), row.names=FALSE, col.names=FALSE,sep = "\t")

      cat(paste("File succesfully stored at ",paste(" ReferenceFiles/Samples/RNAseq/CountFiles/", samples[i,]$Comment.ENA_RUN.,"/",samples[i,]$Comment.ENA_RUN., ".count\n",sep = "")))

    }


  }

}
