#' @param evidences A list which names are genes and which contains found evidences for each gene and a disease. getEvidences output
#' @return A Matrix that contains the redundance between genes evidences ([i,j] contains the redundance between the evidences from gene i with evidences from gene j)
#' @example 
#' test <- list('gen1'= matrix(c('literature','literature','animal_model','code','code','code','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/21441245','TRUE','*','*','mouse'),nrow=3),
#' 'gen2'= matrix(c('literature','literature','animal_model','known_drug','code','code','code','code','score','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/182355','TRUE','http://europepmc.org/abstract/MED/26872','*','*','mouse','pancreatic cancer vs normal'),nrow=4),
#' 'gen3'= matrix(c('animal_model','code','score','TRUE','human'),nrow=1))
#' test
#' evs.red <- getRedundance(test)
#' print(evs.red)
 evidencesToRedundance <- function(evidences){
  result = matrix(0, nrow = length(evidences), ncol = length(evidences))
  unused <- 3
  type.cols <- list('literature'=1,'known_drug'=2,'affected_pathway'=3,'animal_model'=2,
                    'rna_expression'=2,'genetic_association'=2,'somatic_mutation'=3)
  
  for(gen1 in seq(length(evidences))){
    types1 <- unique(evidences[[gen1]][,1])
    for(gen2 in c(1:length(evidences))){
      if (gen1 != gen2){
        for (type in types1){
          rows1 <- which( evidences[[gen1]][,1] == type)
          rows2 <- which( evidences[[gen2]][,1] == type)
          type.total  <- 0
          for (row1 in rows1){
            total <- 0
            nadd <- 0
            for (row2 in rows2){
              ncol <- type.cols[[type]]
              
              act.total <- 0
              for (col in c((unused+1):(unused+ncol))){
                if (!is.na(evidences[[gen1]][row1,col]) && !is.na(evidences[[gen2]][row2,col]) && evidences[[gen1]][row1,col]!='*' && evidences[[gen2]][row2,col]!='*'){
                  if (evidences[[gen1]][row1,col] == evidences[[gen2]][row2,col]){
                    act.total = act.total + 1
                  }
                }
              }
              total = total + act.total / ncol
            }
            type.total = type.total + total 
          }
          result[gen1,gen2] = result[gen1,gen2] + type.total
        }
        result[gen1,gen2] = result[gen1,gen2] / dim(evidences[[gen1]])[1]
      }
      else result[gen1,gen2] = 1
    }
  }
  colnames(result) <- names(evidences)
  rownames(result) <- names(evidences)
  return(result)
}