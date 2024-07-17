
library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

ssn.dir <- "SSN_SignifEdges_BothSets"
dir.create(ssn.dir)
gcn.dir <- "GCN_SignifEdges_BothSets"
dir.create(gcn.dir)

Cancer_ExprM <- readRDS("ExprM/BALL_NBM_norm.RDS")

n_samples <- ncol(Cancer_ExprM)

ssn.files <- list.files(ssn.dir)

Lioness_PerSource <- function(i) {
  
  Source <- rownames(Cancer_ExprM)[i]
  
  gene.name <- gsub("-", "_", Source)
  
  file.name <- paste0(gene.name, "_SSN.csv")
  
  if(!(file.name %in% ssn.files)) {
    
    targets <- which(Source < rownames(Cancer_ExprM))
    
    if(length(targets) > 0){
      
      SSN_df <- list()
      GCN_df <- list()
      
      for(j in 1:length(targets)) {
        #print(j)
        #Check if ID is relevant
        Target <- rownames(Cancer_ExprM)[targets[j]]
        
        ID <- paste(pmin(Source,Target),pmax(Source,Target),sep="_")
        
        #Calculate GCN
        alpha <- Hmisc::rcorr(Cancer_ExprM[rownames(Cancer_ExprM) == Source,],
                              Cancer_ExprM[rownames(Cancer_ExprM) == Target,], type = "spearman")
        p <- alpha$P[1,2]
        
        if(p < 1e-08){
          
          alpha <- alpha$r[1,2]
          
          SSN <- numeric(n_samples)
          
          #Calculate SSN
          for(x in 1:n_samples){
            
            j_cor <- Hmisc::rcorr(Cancer_ExprM[i,-x],
                                  Cancer_ExprM[targets[j],-x], type = "spearman")$r[1,2]
            
            SSN[x] <- round(n_samples*(alpha - j_cor) + j_cor, 4)
            
          }
          
          CancerMedian_SSN <- median(SSN[1:1416])
          Cancer_absolute_deviation <- abs(SSN[1:1416] - CancerMedian_SSN)
          CancerMAD_SSN <- median(Cancer_absolute_deviation)
          
          #Prepare data to save
          
          GCN <- tibble(ID = ID,
                        Rho = round(alpha, 10),
                        p_val = round(p, 10),
                        CancerVar = round(var(SSN[1:1416]), 10),
                        CancerMAD = round(CancerMAD_SSN, 10))
          
          position <- length(SSN_df) + 1
          
          SSN_df[[position]] <- SSN
          GCN_df[[position]] <- GCN
          
          names(SSN_df)[position] <- ID
          names(GCN_df)[position] <- ID
          
        }
      }
      #Save if worked
      if(length(SSN_df) > 0){
        SSN_df <- bind_rows(SSN_df)
        GCN_df <- bind_rows(GCN_df)
        
        Source <- gsub("-", "_", Source)
        
        write.csv(SSN_df, file.path(ssn.dir, paste0(Source, "_SSN.csv")), row.names = F, quote = F)
        write.csv(GCN_df, file.path(gcn.dir, paste0(Source, "_GCN.csv")), row.names = F, quote = F)
      }
    }
  }
}
