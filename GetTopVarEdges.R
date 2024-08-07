library(dplyr)
library(tidyr)
library(future.apply)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)

gcn.dir <- "GCN_SignifEdges_BothSets/"

gcns <- list.files(gcn.dir)

get_varEdges <- function(i){
  
  gcni <- read.csv(paste0(gcn.dir, gcns[i]))
  
  if(nrow(gcni) >= 100){
    
    gcni.var <- gcni %>% arrange(desc(CancerVar))
    
    #gcni.MAD <- gcni %>% arrange(desc(CancerMAD))
    
    gcni.var <- gcni.var[1:100,]
    
    #gcni.MAD <- gcni.MAD[1:100,]
    
  } else {
    
    gcni.var <- gcni
    #gcni.MAD <- gcni
  
    }
  
  return(gcni.var)
  
}

future::plan(multisession, workers = 10)
results <- future_lapply(1:length(gcns), FUN = get_varEdges, future.seed = TRUE)
plan(sequential)

varEdges <- bind_rows(results)

dim(varEdges)

varEdges <- varEdges %>% arrange(desc(CancerVar))

varEdges <- varEdges[1:100,]

saveRDS(varEdges, "Top100varEdges.RDS")

nodes <- varEdges$ID %>% stringr::str_split("_") %>% unlist

sources <- nodes[rep(c(T,F), 100)] %>% unique
sources.f <- gsub("-", "_", sources)

ssn.dir <- "SSN_SignifEdges_BothSets/"
ssns <- list.files(ssn.dir)

ids <- varEdges$ID
ids <- gsub("-", ".", ids)

for(i in 1:length(sources)){
  
  ssni <- read.csv(paste0(ssn.dir, paste0(sources.f[i], "_SSN.csv")))
  
  edges <- colnames(ssni)[which(colnames(ssni) %in% ids)]
  
  if(length(edges) > 1){
    
    ssni <- ssni[,colnames(ssni) %in% edges]
    
  } else {
    
    ssni <- data.frame(weigths = ssni[,colnames(ssni) %in% edges])
    
    colnames(ssni) <- edges
    
  }
  
  if(i == 1){
    
    SSNs <- ssni
    
  } else {
    SSNs <- bind_cols(SSNs, ssni)
  }

}

saveRDS(SSNs, "Top100_SSNs_BothDatasets.RDS")
