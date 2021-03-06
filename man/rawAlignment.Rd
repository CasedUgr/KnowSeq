% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawAlignment.R
\name{rawAlignment}
\alias{rawAlignment}
\title{rawAlignment allows downloading and processing the fastq samples in a CSV file.}
\usage{
rawAlignment(
  data,
  downloadRef = FALSE,
  downloadSamples = FALSE,
  createIndex = TRUE,
  BAMfiles = TRUE,
  SAMfiles = TRUE,
  countFiles = TRUE,
  referenceGenome = 38,
  customFA = "",
  customGTF = "",
  fromGDC = FALSE,
  tokenPath = "",
  manifestPath = "",
  hisatParameters = "-p 8 --dta-cufflinks"
)
}
\arguments{
\item{data}{The ID of the variable which contains the samples. Our recommendation is to load this variable from a CSV file.}

\item{downloadRef}{A logical parameter that represents if the reference files will be downloaded or not.}

\item{downloadSamples}{A logical parameter that represents if the samples of the CSV file will be downloaded or not.}

\item{createIndex}{A logical parameter that represents if the index of the aligner would be created or not.}

\item{BAMfiles}{A logical parameter that represents if the you want the BAM files or not.}

\item{SAMfiles}{A logical parameter that represents if the you want the SAM files or not.}

\item{countFiles}{A logical parameter that represents if the you want the Count files or not.}

\item{referenceGenome}{This parameter allows choosing the reference genome that will be used for the alignment. The options are 37,38 or custom. The two first are human genomes, but with the third option you can choose any genome stored in the computer.}

\item{customFA}{The path to the custom FASTA file of the reference genome.}

\item{customGTF}{The path to the custom GTF file.}

\item{fromGDC}{A logical parameter that allows processing BAM files from GDC portal by using the custom reference genome from GDC.}

\item{tokenPath}{The path to the GDC portal user token. It is required to downloads the controlled BAM files.}

\item{manifestPath}{The path to the manifest with the information required to downloads the controlled BAM files selected in GDC Portal.}

\item{hisatParameters}{Parameter that allow to modify the default configuration for the Hisat2 aligner.}
}
\value{
Nothing to return.
}
\description{
This function allows downloading and processing the fastq samples in a CSV file. Also, samples can be aligned by using hisat2. Finally, the function can downloads the reference files required: FASTA Reference Genome and GTF file.
}
\examples{
# Due to the high computational cost, we strongly recommend it to see the offical documentation and the complete example included in this package:

dir <- system.file("extdata", package="KnowSeq")

#Using read.csv for NCBI/GEO files (read.csv2 for ArrayExpress files)
GSE74251csv <- read.csv(paste(dir,"/GSE74251.csv",sep = ""))

\dontrun{rawAlignment(GSE74251csv,downloadRef=FALSE,downloadSamples=FALSE, createIndex = TRUE, BAMfiles = TRUE, SAMfiles = TRUE, countFiles = TRUE, referenceGenome = 38, customFA = "", customGTF = "", fromGDC = FALSE, tokenPath = "", manifestPath = "")}
}
