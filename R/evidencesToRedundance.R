#' @param data A list which names are genes and which contains found data for each gene and a disease. getdata output
#' @return A Square Matrix that contains the redundance between genes data ([i,j] contains the redundance between the data from gene i with data from gene j)
#' @example 
#' test <- list('gen1'= matrix(c('literature','literature','animal_model','code','code','code','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/21441245','TRUE','*','*','mouse'),nrow=3),'gen2'= matrix(c('literature','literature','animal_model','known_drug','code','code','code','code','score','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/182355','TRUE','http://europepmc.org/abstract/MED/26872','*','*','mouse','pancreatic cancer vs normal'),nrow=4),'gen3'= matrix(c('animal_model','code','score','TRUE','human'),nrow=1))
# evs.red <- evidencesToRedundance(test)
#' evs.red

evidencesToRedundance <- function(data){
   cat('Calculating redundance between genes...')

   # Create empty output
   result = matrix(0, nrow = length(data), ncol = length(data))

   # Iterate in genes
   for(gen1 in seq(length(data))){
     for(gen2 in seq(length(data))){
       if (gen1 != gen2){
         # If evidence from gen1 and gen2 are not empty
         if ( length(data[[gen1]]$evidences) > 0 && length(data[[gen2]]$evidences) > 0){
           # Iter in types
            gen1.nevs <- 0
            for (type in names(data[[gen1]]$evidences)){
               if (type %in% names(data[[gen1]]$evidences)){
                  # type.total contains found coincidences for actual type of data
                  type.total  <- 0
                  
                   # Iter on gen1 data
                   for (row1 in data[[gen1]]$evidences[[type]]){
                     gen1.nevs <- gen1.nevs + 1
                     # Boolean matrix that contains coincidences between gen1 and gen2
                     print(length(data[[gen2]]$evidences[[type]]))
                     print(data[[gen2]]$evidences[[type]])
                     act.total.matrix <- matrix(0,nrow=length(data[[gen2]]$evidences[[type]]),ncol=ncol)
                     # Iter on gen1 data
                     for (row2 in data[[gen2]]$evidences[[type]]){
                        # Add row to act.total.matrix with boolean values
                        act.total.matrix <- rbind(act.total.matrix, row1$evidence == row2$evidence)
                        # This row fully coincide with actual evidencie, so we stop searching
                        if (rowSums(tail(act.total.matrix,1) == ncol))  break
                     }
               
                     # If there are any row that fully coincide add 1 (this data is fully contained in gen2 data)
                     if ( any(rowSums(act.total.matrix) == ncol)) type.total = type.total + 1
                     else{
                       # Check which row is the most similar to the actual evidence
                       # Add the percentage of coincidence ( num. coincidences / ncol)
                       act.total.field <- colSums(act.total.matrix)
                       type.total = type.total + length(which(colSums(act.total.matrix) >= 1))/ncol
                     }
                  }
                  # Add score of type data
                  result[gen1,gen2] = result[gen1,gen2] + type.total
               }
            }
           # Normalize score dividing by the number of evidences
           result[gen1,gen2] = result[gen1,gen2] / gen1.nevs
         }else result[gen1,gen2] = 0 # No data = no redundance
       }
       else result[gen1,gen2] = 1 # If gen1 = gen2
     }
   }
   colnames(result) <- names(data)
   rownames(result) <- names(data)
   return(result)
}