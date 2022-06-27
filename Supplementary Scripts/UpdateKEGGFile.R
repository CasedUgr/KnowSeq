require(httr)
require(stringr)

get_GO <- GET("https://rest.kegg.jp/list/pathway")
get_GO_text <- content(get_GO, "text")
write.table(get_GO_text[1], "KEGGPathways.tsv", sep = "\t",col.names = FALSE, row.names = FALSE)

pathways <- read.table("KEGGPathways.tsv", sep = "\t")

KEGG.data <- matrix(ncol = 4)

for(pathway in pathways$V1){
  
  path_endpoint <- GET(paste("https://rest.kegg.jp/get/", substr(pathway,6,nchar(pathway)), sep = ""))
  path_info <- content(path_endpoint,"text")

  start <- str_locate_all(pattern = "ENTRY", path_info)[[1]][2]
  final <- str_locate_all(pattern = "Pathway", path_info)[[1]][2]
  path.index <- substr(path_info,start,final)
  path.index <- strsplit(path.index,split = " ")[[1]]
  path <- path.index[path.index != ""]
  path <- path[2]
  
  start <- str_locate_all(pattern = "\nNAME ", path_info)[[1]][2]
  final <- str_locate_all(pattern = "\nDESCRIPTION", path_info)[[1]][2]
  name <- substr(path_info,start + 8,final - 12)
  
  start <- str_locate_all(pattern = "\nDESCRIPTION ", path_info)[[1]][2]
  final <- str_locate_all(pattern = "\nCLASS", path_info)[[1]][2]
  description <- substr(path_info,start + 1,final - 6)
  description <- gsub("\"","",description)
  
  start <- str_locate_all(pattern = "\nCLASS ", path_info)[[1]][2]
  final <- str_locate_all(pattern = "\nPATHWAY_MAP", path_info)[[1]][2]
  class <- substr(path_info,start + 7,final - 12)
  
  KEGG.data <- rbind(KEGG.data,c(as.character(path),as.character(name), as.character(description),as.character(class)))
  
}

colnames(KEGG.data) <- c("KEGG_Id", "Name", "Description", "Class")
KEGG.data = KEGG.data[-1,]
write.csv(KEGG.data, "KEGGPathsDB.csv", col.names = T, row.names = F, na = "nothing")
