

library(dplyr)
library(tidyr)
library(M3C) #en caso de que vayas a hacer clustering con M3C, si no, lo puedes hacer en el heatmap (puse instrucciones abajo)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#Filtra genes más variables

 filtered_results <- featurefilter(Cancer_Degrees,
                                    percentile=15, #porcetaje de features mas variables a elegir
                                    method='MAD', #puede ser MAD o var
                                    topN=10)
  
  Cancer_Degrees.scaled <- as.data.frame(t(scale(t(filtered_results$filtered_data)))) #Escala los datos


#Clustetering analysis
  r.hc_Discovery <- M3C(Cancer_Degrees.scaled, method=1, 
                        clusteralg = "spectral", cores = 30, 
                        maxK = 20, seed = 123) #Va a imprimir un número óptimo de clusters cuando termine

    OptimalK <- k #k = número optimo
    assigments.Discovery <- r.hc_Discovery$realdataresults[[OptimalK]]$ordered_annotation
    assigments.Discovery <- assigments.Discovery %>% mutate(Sample = rownames(assigments.Discovery))
    assigments.Discovery <- assigments.Discovery %>% arrange(consensuscluster)

    data <- r.hc_Discovery$realdataresults[[OptimalK]]$ordered_data

    #Guarda UMAP y PCA del OptimalK
    umap <- umap(data,labels=assigments.Discovery$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)
    ggsave(paste0(Clustering.dir,  "/UMAP_K", k, "_Top", top, ".png") , plot = umap)
    pca <- pca(data,labels=assigments.Discovery$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)
    ggsave(paste0(Clustering.dir,  "/PCA_K",k , "_Top", top, ".png"), plot = pca)

#Guarda también estas figuras
r.hc_Discovery$plots[[4]] #RCSI
r.hc_Discovery$plots[[3]] #Pvalues

#Crea una paleta de colores para tus cluters
clusters_unique <- names(table(assigments.Discovery$consensuscluster))
    clusters_palette <- colorRampPalette(brewer.pal(9, "Paired"))(length(clusters_unique))
    clusters_palette <- setNames(clusters_palette, clusters_unique)
    
    column_ha_NCCI = HeatmapAnnotation(Cluster = assigments.Discovery$consensuscluster,
                                       col = list(Cluster = clusters_palette))

#Tu matriz y tu assigments deben tener el mismo orden de muestras
Cancer_Degrees.scaled <- Cancer_Degrees.scaled[,match(assigments.Discovery$Sample, colnames(Cancer_Degrees.scaled))]

#Heatmap
set.seed(123)
    NCCI <- ComplexHeatmap::Heatmap(as.matrix(Cancer_Degrees.scaled), #El objeto tiene que se runa matriz
                                    cluster_columns = F, cluster_rows = T, #Cluster rows 
                                    #column_km = 3, column_km_repeats = 100, #Esto hace un clusterin jerarquico, que no haremos porque ya hicimos uno arriba
                                    #Sin embargo, si no hiciste clustering, puedes hacer uno con esos argumentos, km_repeats no lo cambies,
                                    #solo cambia column_km para el numero de clusters que desees
                                    show_row_names = F, show_column_names = F,
                                    #row_km = 4, row_km_repeats = 100,
                                    use_raster = FALSE,
                                    top_annotation = column_ha_NCCI, #Esto es opcional
                                    #left_annotation = row_ha_NCCI,
                                    name = "Scaled Gene Expression",
                                    #row_names_gp = grid::gpar(fontsize = 5),
                                    col = colorRamp2(c(-1,  0 , 1), #Cambia el intervalo que quieras, recuerda que tus datos los transformamos a Z-scores
                                                     c("black", "yellow", "red"))) #Cambia los colores que quieras, pueden ser Verde, Negro y Rojo
    
    png("NombreDeArchivoQueQuieras.png", units="in", width=10, height=15, res=300) #Cambia width y height como desees, pero creo que 10 y 15 sería adecuado res dejalo en 300
    #Mientras más grande sea tu figura, más memoría usara

    draw(NCCI)
    
    dev.off()
#despues de correr 57, 60 y 62 tendrás el archivo png para que descargues

