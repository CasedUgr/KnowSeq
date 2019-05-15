#' This function downloads a list of controlled files from GDC Portal with the user token and the manifest with the information about the desired controlled files.
#'
#' @param tokenPath Path to the GDC token
#' @param manifestPath Path to the samples manifest
#' @param data The matrix or data.frame with the information from the Samples Sheet downloaded from GDC Portal.
#' @return Nothing to return.
#' @examples
#' gdcClientDownload("~/user/token.txt","~/user/manifest.txt", data)

gdcClientDownload <- function(tokenPath, manifestPath, data) {

  if(file.exists(tokenPath)){cat("Token found!\n")}else{stop("Token not found, please revise the path to the token.\n")}
  if(file.exists(manifestPath)){cat("Manifest found!\n")}else{stop("Manifest not found, please revise the path to the manifest.\n")}

  system(paste("unixUtils/gdcClient/gdc-client download -t ",tokenPath," -m ", manifestPath," -n 64 --no-segment-md5sums --no-file-md5sum"))

  cat("Moving the downloaded files to ReferenceFiles/Samples/RNAseq/BAMFiles/ \n")

  for(i in seq_len(dim(data)[1])){
    system(paste("mv ",data$File.ID[i], " ReferenceFiles/Samples/RNAseq/BAMFiles/",sep = ""))

  }

}
