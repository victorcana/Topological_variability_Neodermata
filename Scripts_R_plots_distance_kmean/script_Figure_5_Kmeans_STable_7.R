# Set working directory
setwd("Kmeans_Figure_5/")

# Set seed for reproducibility
set.seed(123)

# Load required libraries
library(readxl)
library(ggfortify)
library(factoextra)
library(cluster)

# Load the data
Todos_arboles_listo_nwk <- read_excel("Todos_arboles_listo.nwk.xlsx")

# ------------------------------------------------------
# Determine optimal number of clusters
# ------------------------------------------------------

# Silhouette method
fviz_nbclust(Todos_arboles_listo_nwk[,7:111], kmeans, method = "silhouette", k.max = 20)+
  geom_vline(xintercept = 10, linetype = 2)
#ggsave("plot_nbclustsilhouete.png", width = 30, height = 20,units = "cm",limitsize = FALSE)

# Elbow method (WSS)
fviz_nbclust(Todos_arboles_listo_nwk[,7:111], kmeans, method = "wss", k.max = 40) +
  geom_vline(xintercept = 10, linetype = 2)
#ggsave("plot_nbclustwss.png", width = 30, height = 20,units = "cm",limitsize = FALSE)

# Gap statistic method
fviz_nbclust(Todos_arboles_listo_nwk[,7:111], kmeans, method = "gap", k.max = 20)


# ------------------------------------------------------
# Perform K-means clustering (k = 10)
# ------------------------------------------------------
set.seed(123)
km.result <- kmeans(Todos_arboles_listo_nwk[,7:111], 10)  ###Se seleccionaron 4 cluster porque fue el segundo mejor evaluado
km.result

# ------------------------------------------------------
# Visualize silhouette widths for each observation
# ------------------------------------------------------
set.seed(123)
km.sil <-cluster::silhouette(km.result$cluster, dist(Todos_arboles_listo_nwk[,7:111]), ordered = FALSE)
row.names(km.sil) <- Todos_arboles_listo_nwk$ID # Needed to use label option
fviz_silhouette(km.sil,palette = "jco", label = TRUE) +
  theme_ggstatsplot()+
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
#ggsave("plot_silhouete.png",width = 35, height = 17,units = "cm", limitsize = FALSE)
#ggsave("plot_silhouete.pdf",width = 35, height = 17,units = "cm", limitsize = FALSE)
km.result$cluster

# ------------------------------------------------------
# Save clustering result as table
# ------------------------------------------------------
datos_ncikm <-cbind(Todos_arboles_listo_nwk[,7:111], Cluster = km.result$cluster)
write.table(datos_ncikm, file = "S_Table_7_datos_nci_clusteringkm.tsv", sep = "\t", row.names = FALSE)


# ------------------------------------------------------
# Visualize K-means clusters using PCA
# ------------------------------------------------------
Todos_arboles_listo_nwk2 = as.data.frame(Todos_arboles_listo_nwk)
rownames(Todos_arboles_listo_nwk2) <- Todos_arboles_listo_nwk2$ID

# PCA plot with cluster coloring
fviz_cluster(km.result, data = Todos_arboles_listo_nwk2[,7:111],labelsize = 10, repel = TRUE ,
             star.plot = TRUE, star.plot.lty = 1, star.plot.lwd = 0.07,alpha=0.0,main =" ", obs.scale = 3, var.scale = 1,
             pointsize =8,
             show.clust.cent = FALSE,
             ellipse.type = "convex", palette = "jco"
             ) +
  geom_point( alpha = 0.4, size = 1)+ 
 scale_x_continuous(breaks = seq(-10, 20, 10), limits = (c(-10, 23))) +
  theme_classic()+theme(
    axis.title.x = element_text(size = 13),  
    axis.title.y = element_text(size = 13),  
    axis.text.x = element_text(size = 11),  
    axis.text.y = element_text(size = 11)
    )    

# Save PCA cluster plot
ggsave("Figure_5_plot_clustered_kmean.png",width = 35, height = 20,units = "cm", limitsize = FALSE)
ggsave("Figure_5_plot_clustered_kmean.pdf",width = 35, height = 20,units = "cm", limitsize = FALSE)
