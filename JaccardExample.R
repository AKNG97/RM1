library(dplyr)
library(tidyr)
library(M3C)
library(ggsurvfit)
library(survival)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(parallel)
library(future.apply)

Lioness_GOYA <- read.csv("/STORAGE/csbig/anakamura/DLBCL_A4/A5/SSNsCommonGenes/SSN_GOYA.csv", row.names = 1)

Lioness_GOYA.scale <- as.data.frame(t(scale(t(Lioness_GOYA))))

Lioness_GOYA.scale <- Lioness_GOYA.scale %>% mutate(Edge = rownames(Lioness_GOYA.scale)) %>%
  pivot_longer(!Edge, values_to = "Z_score", names_to = "Sample") %>% 
  group_by(Edge) %>%
  mutate(EdgeClass = case_when(Z_score > 1 ~ paste0("+", Edge),
                               Z_score < -1 ~ paste0("-", Edge),
                               .default = paste0("Neutral", Edge))) %>% ungroup %>%
  dplyr::select(Sample, Z_score, EdgeClass)

samples <- unique(Lioness_GOYA.scale$Sample)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

dir.create("jaccard_summary_GOYA3/")

options(future.globals.maxSize = 8 * 1024^3)

jaccard_par <- function(x, TopEdges){
  
  sample1 <- samples[x]
  
  jaccard_list <- list()
  
  s1_edges = Lioness_GOYA.scale[Lioness_GOYA.scale$Sample == sample1,] %>%
    arrange(desc(abs(Z_score))) %>% head(TopEdges) %>% pull(EdgeClass)
  
  for(y in 1:length(samples)){
    
    print(y)
    
    sample2 <- samples[y]
    
    if(y != x && sample1 < sample2){
      
      s2_edges = Lioness_GOYA.scale[Lioness_GOYA.scale$Sample == sample2,] %>%
        arrange(desc(abs(Z_score))) %>% head(TopEdges) %>% pull(EdgeClass)
      
      jaccard_i <- jaccard(s1_edges, s2_edges)
      
      jaccard_list[[length(jaccard_list) + 1]] <- c(sample1, sample2, jaccard_i)
      
    }
  }
  
  if(length(jaccard_list) > 0){
    jaccard_summary <- as.data.frame(do.call(rbind, jaccard_list))
    colnames(jaccard_summary) <- c("Sample1", "Sample2", "JaccardIndex")  
    write.csv(jaccard_summary, 
              paste0("jaccard_summary_GOYA3/Jacard_", sample1, ".csv"),
              quote = F,
              row.names = F)
    
  }
  
  
}
