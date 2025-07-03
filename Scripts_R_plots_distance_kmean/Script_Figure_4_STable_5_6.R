# Load necessary libraries
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggstatsplot)
library(statsExpressions)
# Set working directory containing .mldist files
setwd("Distance_Figure_4")

# Initialize list to store data
datos <- list()    
i <- 1
# Read each file and convert to long-format distance matrix
for (file in list.files()) {  
  datos[[i]] <- read.table(file, row.names = 1); 
  colnames(datos[[i]]) <- rownames(datos[[i]]); 
  datos[[i]][datos[[i]] == "0"] <- NA; 
  datos[[i]] <- cbind(datos[[i]], rownames(datos[[i]]));
  datos[[i]]<- datos[[i]]%>%pivot_longer(cols = 1:18, names_to = "Species", values_to = "Distance");
  datos[[i]]<-unite(datos[[i]], Pairwise,c(1:2),  sep = "-", remove = FALSE)
  colnames(datos[[i]]) <- c('Pairwise','SpeciesA','SpeciesB', file)
  if (i == 1)
    A <- (datos[[i]][c(1,4)])
  else
    A <- merge(x = A, y = datos[[i]][c(1,4)], all = TRUE, by = "Pairwise" )
  i <- i + 1
}


# Remove rows with more than 20% missing values
B <- A[-which(rowMeans(is.na(A)) >= 0.20),]
B

# Convert to long format for plotting and stats
C <- B %>% pivot_longer(cols = 2:32, names_to = "Genes", values_to = "Distance")
C

# Extract gene identity and replicate info from column name
C <- cbind(C, C[2])
colnames(C) <- c("Pairwise","Genes", "Distance", "Genes_2")
C <- separate(C, Genes_2, c("Genes_2","Sequences"))
C

# ----------------------------------------------
# Generate plots for nucleotide and amino acid distances
# ----------------------------------------------

# Plot for nucleotide distances
library("ggstatsplot")
set.seed(123)
P1<- ggbetweenstats(
  data  = C%>% dplyr::filter( str_detect(Genes, "NT")), ggsignif.args = list(textsize = 3, tip_length = 0.005, vjust = 0.7),
  x     = Genes_2,
  y     = Distance,   type            = "np", 
  ylab = "Phylogenetic distance",
   pairwise.comparisons= FALSE,# conf.level = 0.50,
  title = "Nucleotides", xlab = "Genes",
  plot.type = "boxviolin", 
  ggtheme = theme_classic(),
  pairwise.display ="none", package = "pals", palette = "alphabet"#package = "yarrr", palette = "info2",
  #outlier.tagging = TRUE, outlier.label = Pairwise
)  +
  scale_x_discrete(limits=c("MitorRNA", "Mitogenome","atp6","cox1","cox2", "cox3" ,"cytb", "nad1","nad2","nad3","nad4","nad4L" ,"nad5","nad6","12S","16S","18S","28S")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 


# Plot for amino acid distances
set.seed(123)
P2 <- ggbetweenstats(
  data  = C%>% filter( str_detect(Genes, "AA")), ggsignif.args = list(textsize = 3, tip_length = 0.005, vjust = 0.7),
  x     = Genes_2,
  y     = Distance,   type            = "np", 
  ylab = "",
   pairwise.comparisons= FALSE,# conf.level = 0.50,
  title = "Amino acids", xlab = "Genes",
  plot.type = "boxviolin", 
  ggtheme = theme_classic(),
  pairwise.display ="none", package = "pals", palette = "alphabet"#package = "yarrr", palette = "info2",
  #outlier.tagging = TRUE, outlier.label = Pairwise
)  +
  scale_x_discrete(limits=c("Mitogenome","atp6","cox1","cox2", "cox3" ,"cytb", "nad1","nad2","nad3","nad4","nad4L" ,"nad5","nad6"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Combine the two plots into one figure
P3<- combine_plots(
  list(P1, P2),
  plotgrid.args = list(nrow = 1))

P3

# Save plots to file
ggsave(plot = P3, "ZFigure_4_plot_distance.png", width = 42.6, height = 13.3,units = "cm",limitsize = FALSE)
ggsave(plot = P3, "ZFigure_4_plot_distance.pdf",width = 42.6, height = 13.3,units = "cm", limitsize = FALSE)


# ----------------------------------------------
# Non-parametric statistical tests
# ----------------------------------------------

# Kruskal-Wallis test for amino acid distances
set.seed(123)
oneway_anova(
 data = C%>%filter( str_detect(Genes, "AA")),
  x     = Genes_2,
  y     = Distance,
  type      = "np"
)

# Pairwise comparisons – amino acids
Statis_AA <-pairwise_comparisons(
  data = C%>%filter( str_detect(Genes, "AA")),
  x     = Genes_2,
  y     = Distance,
  type            = "nonparametric",
  paired          = FALSE,
  p.adjust.method =  "hommel"
)

# Pairwise comparisons – nucleotides
Statis_NT <-pairwise_comparisons(
  data = C%>%filter( str_detect(Genes, "NT")),
  x     = Genes_2,
  y     = Distance,
  type            = "nonparametric",
  paired          = FALSE,
  p.adjust.method =  "hommel"
)

# Convert to character format and save results as tables
Statis_NT <- apply(Statis_NT,2,as.character)
write.table(Statis_NT,"S_Table_5_Stadistis_NT", row.names = FALSE, sep = "\t")
Statis_AA <- apply(Statis_AA,2,as.character)
write.table(Statis_AA,"S_Table_6_Stadistis_AA", row.names = FALSE, sep = "\t")

