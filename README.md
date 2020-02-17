# KnowSeq R/bioc package: Beyond the traditional RNA-seq pipeline.

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/KnowSeq.svg)](https://bioconductor.org/packages/release/bioc/html/KnowSeq "How long has been KnowSeq in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/KnowSeq.svg)](https://bioconductor.org/packages/stats/bioc/KnowSeq "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months")
[![Support posts](https://bioconductor.org/shields/posts/KnowSeq.svg)](https://support.bioconductor.org/t/KnowSeq/ "Support site activity on KnowSeq, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")


**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/3.10/KnowSeq.svg)](https://bioconductor.org/packages/release/bioc/html/KnowSeq.html#archives "Whether KnowSeq release is available on all platforms") 
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/KnowSeq.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/KnowSeq "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/KnowSeq.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/KnowSeq/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/3.11/KnowSeq.svg)](https://bioconductor.org/packages/devel/bioc/html/KnowSeq.html#archives "Whether KnowSeq devel is available on all platforms") 
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/KnowSeq.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/KnowSeq "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/KnowSeq.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/KnowSeq/ "Bioconductor devel build")

# Introduction

[KnowSeq](https://bioconductor.org/packages/release/bioc/html/KnowSeq.html) proposes a whole pipeline that comprises the most relevant steps in the RNA-seq gene expression analysis, with the main goal of extracting biological knowledge from raw data (Differential Expressed Genes, Gene Ontology enrichment, pathway visualization and related diseases). In this sense, KnowSeq allows aligning raw data from the original fastq or sra files, by using the most renowned aligners such as tophat2, hisat2, salmon and kallisto. Nowadays, there is no package that only from the information of the samples to align -included in a text file-, automatically performs the download and alignment of all of the samples. Furthermore, the package includes functions to: calculate the gene expression values; remove batch effect; calculate the Differentially Expressed Genes (DEGs); plot different graphs; and perform the DEGs enrichment with the GO information, pathways visualization and related diseases information retrieval. Moreover, KnowSeq is the only package that allows applying both a machine learning and DEGs enrichment processes just after the DEGs extraction. To achieve these objectives, there are functions that allows performing a feature selection process as well as a machine learning process using k-NN, RF or SVM algorithms. Similarly, there are functions allowing the retrieval of biological knowledge of the DEGs candidates. This idea emerged with the aim of proposing a complete tool to the research community containing all the necessary steps to carry out complete studies in a simple and fast way. To achieve this goal, the package uses the most relevant and widespread tools in the scientific community for the aforementioned tasks. The current version of the aligner functions works under Unix, but further version will be extended to MAC_OS and to Windows (if the tools were available). This pipeline has been used in our previous publications for processing raw RNA-seq data and to perform the DEGs extraction and the machine learning classifier design steps, also for their integration with microarray data [1,2,3].

![](https://github.com/CasedUgr/KnowSeq/blob/master/vignettes/KnowSeqPipeline.png)

The whole pipeline included in KnowSeq has been designed carefully with the purpose of achieving a great quality and robustness in each of the steps that conform the pipeline. For that, the pipeline has three fundamental processes:

- RNA-seq RAW data processing.
- Biomarkers identification & assessment.
- DEGs enrichment methodology.

The first process is focused on the RNA-seq RAW data treatment. This step has the purpose of extracting a set of count files from raw files stored in the repositories supported by our package (NCBI/GEO [4] ArrayExpress [5] and GDC-Portal). The second one comprises the Differential Expressed Genes (DEGs) identification and extraction, and the assessment of those DEGs by applying advanced machine learning techniques (feature selection process and supervised classification). The last process, once the DEGs were assessed, is the DEGs enrichment methodology which allows retrieving biological information from the DEGs. In this process, relevant information (such as related diseases, biological processes associated and pathways) about the DEGs is retrieved by using very well-known tools and databases. The three types of enrichment are the Gene Ontology (GO) study, the pathways visualization taking into account the gene expression, and the related diseases to the DEGs.
With the pipeline designed and addressed by KnowSeq, researchers can convert the RAW data of RNA-seq into real knowledge on the identification of possible gene signatures about the studied diseases.

# Installation
To install and load KnowSeq package in R, it is necessary the previous installation of BiocManager from Bioconductor. The next code shows how this install can be performed:

```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
BiocManager::install("KnowSeq")

library(KnowSeq)
```

KnowSeq is now available on Docker by running the next command:

```
Docker run -it casedugr/knowseq
```

# Contributors

Contributors are listed in the [DESCRIPTION](https://github.com/CasedUgr/KnowSeq/blob/master/DESCRIPTION) file.

# Contributions

Contributions are welcomed and each PR would be reviewed by KnowSeq maintainers. Please, provide a useful description for the new contribution. PR guidelines will be added soon.

# Citation

If you find KnowSeq useful and you use it in your work, please cite it as follows:

```
Castillo-Secilla D, Galvez JM, Ortuno FM, Herrera LJ, Rojas. I (2019). KnowSeq: A R package to extract knowledge by using RNA-seq raw files. R package version 1.0.0. 
```
# References

1. Castillo, D., Gálvez, J. M., Herrera, L. J., San Román, B., Rojas, F., & Rojas, I. (2017). Integration of RNA-Seq data with heterogeneous microarray data for breast cancer profiling. BMC bioinformatics, 18(1), 506.

2. Gálvez, J. M., Castillo, D., Herrera, L. J., San Roman, B., Valenzuela, O., Ortuno, F. M., & Rojas, I. (2018). Multiclass classification for skin cancer profiling based on the integration of heterogeneous gene expression series. PloS one, 13(5), e0196836.

3. Castillo, D., Galvez, J. M., Herrera, L. J., Rojas, F., Valenzuela, O., Caba, O., ... & Rojas, I. (2019). Leukemia multiclass assessment and classification from Microarray and RNA-seq technologies integration at gene expression level. PloS one, 14(2), e0212127.
