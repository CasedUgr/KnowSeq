#' @param evidences A list which names are genes and which contains found evidences for each gene and a disease. getEvidences output
#' @return A Square Matrix that contains the redundance between genes evidences ([i,j] contains the redundance between the evidences from gene i with evidences from gene j)
#' @example 
#' test <- list('gen1'= matrix(c('literature','literature','animal_model','code','code','code','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/21441245','TRUE','*','*','mouse'),nrow=3),'gen2'= matrix(c('literature','literature','animal_model','known_drug','code','code','code','code','score','score','score','score','http://europepmc.org/abstract/MED/6783302','http://europepmc.org/abstract/MED/182355','TRUE','http://europepmc.org/abstract/MED/26872','*','*','mouse','pancreatic cancer vs normal'),nrow=4),'gen3'= matrix(c('animal_model','code','score','TRUE','human'),nrow=1))
# evs.red <- evidencesToRedundance(test)
#' evs.red
#' print(evs.red)
evidencesToRedundance <- function(evidences){
   cat('Calculating redundance between genes...')
   
   # Create empty output
   result = matrix(0, nrow = length(evidences), ncol = length(evidences))
   # There are 3 unused rows for redundance: type, code and score
   unused <- 3 
   # Each type of evidence has a fixed numer of usefull columns
   type.cols <- list('literature'=1,'known_drug'=2,'affected_pathway'=3,'animal_model'=2,
                     'rna_expression'=2,'genetic_association'=2,'somatic_mutation'=3)
   
   # Iterate in genes
   for(gen1 in seq(length(evidences))){
     for(gen2 in seq(length(evidences))){
       if (gen1 != gen2){
         # If evidence from gen1 and gen2 are not empty
         if (class(evidences[[gen1]]) == 'matrix' && class(evidences[[gen2]]) == 'matrix' && dim(evidences[[gen1]])[1]>0){
           # Get evidence types in gen1
           types1 <- unique(evidences[[gen1]][,1])
           # Iter in types
           for (type in types1){
             # Rows from gen1 and gen2 which contain evidence of the actual type
             rows1 <- which( evidences[[gen1]][,1] == type)
             rows2 <- which( evidences[[gen2]][,1] == type)
             # Usefull columns of actual type
             ncol <- type.cols[[type]]
             
             # type.total contains found coincidences for actual type of evidences
             type.total  <- 0
             
             # Iter on gen1 evidences
             for (row1 in rows1){
               # Boolean matrix that contains coincidences between gen1 and gen2
               act.total.matrix <- matrix(0,nrow=length(rows2),ncol=ncol)
               # Iter on gen1 evidences
               for (row2 in rows2){
                 # Add row to act.total.matrix with boolean values
                 act.total.matrix <- rbind(act.total.matrix, evidences[[gen1]][row1,(unused+1):(unused+ncol)] == evidences[[gen2]][row2,(unused+1):(unused+ncol)] )
               }
               
               # If there are any row that fully coincide add 1 (this evidences is fully contained in gen2 evidences)
               if ( any(rowSums(act.total.matrix) == ncol)) type.total = type.total + 1
               else{
                 # Check which row is the most similar to the actual evidence
                 # Add the percentage of coincidence ( num. coincidences / ncol)
                 act.total.field <- colSums(act.total.matrix)
                 type.total = type.total + length(which(colSums(act.total.matrix) >= 1))/ncol
               }
             }
             # Add score of type evidences
             result[gen1,gen2] = result[gen1,gen2] + type.total
           }
           # Normalize score dividing by the number of evidences
           result[gen1,gen2] = result[gen1,gen2] / dim(evidences[[gen1]])[1]
         }else result[gen1,gen2] = 0 # No evidences = no redundance
       }
       else result[gen1,gen2] = 1 # If gen1 = gen2
     }
   }
   colnames(result) <- names(evidences)
   rownames(result) <- names(evidences)
   return(result)
}