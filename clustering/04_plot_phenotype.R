# Step4.Phenotype boxplot

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/clustering/")
library(plot3D)
library(ggplot2)

data = read.csv("data/clustering_result.csv")

{
  # Read and prepare "landscapes_overview.xlsx"
  library(readxl)
  table = read_excel("../Data/20240207_landscapes_overview.xlsx")
  colnames(table)[3:ncol(table)] = table[1,3:ncol(table)]
  table = table[-1,]
  table[table == "-"] = NA
  table = data.frame(table)
  for (i in 3:ncol(table))
    table[,i] = as.numeric(table[,i])
  for (i in 1:nrow(table)){
    g = table$PAH.genotype[i]
    g = strsplit(g, "-")[[1]][1]
    g = strsplit(g, "\\+")[[1]][1]
    g = strsplit(g, "/")[[1]]
    if (g[1] == "WT")
      g = "WT"
    else
      g = paste0(g[1], "-", g[2])
    table$PAH.genotype[i] = g
  }
  colnames(table)[1] = "genotype"
  
  # pheno = table$Classic_PKU...1200µmol.L.Phe....
  # pheno = table$Mild_PKU..600.1200µmol.L.Phe....
  pheno = table$MHPA...600µmol.L.Phe....
  
  pheno[is.na(pheno)] = -100
  pheno = pheno / 100
  table$pheno = pheno
  table$pheno[pheno == -1] = NA
  table$pheno[table$Phenotype..Nr.of..tested.patients < 10] = NA

  # Merge
  # Merge the two tables
  table$genotype_exp = paste0(table$genotype, "|", table$EXPERIMENT)
  data$genotype_exp = paste0(data$genotype, "|", data$experiment)
  table = merge(table, data, by = "genotype_exp", suffixes=c("", ".y"))
}

X = table[,c("genotype", "experiment", "Max_theory", "Max_x", "Max_y", "pheno", "clusters")]
X_pheno = X[!is.na(X$pheno),]

order = c(0, 2, 5, 1, 3, 4)
col.pal = c("#2C1453", "#A48AD3", "#1CC5FE", "#6FC7CF", "#FBA27D", "#FB7D80")
col.pal = c("red", "orange", "green", "skyblue", "blue", "magenta")

# pheno box plot
X_pheno_no6 = X_pheno[X_pheno$clusters < 6,]
X_pheno_no6$clusters_factor = factor(X_pheno_no6$clusters, 
                                        levels = order)

pdf(file = "Phenotype_MHPA.pdf", width = 5, height = 5)

ggplot(X_pheno_no6, aes(x=clusters_factor, y=pheno, color=clusters_factor)) + 
  geom_boxplot(fill="gray95", size = .6) + 
  geom_point() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values = col.pal[order+1]) + 
  theme_classic() +
  theme(legend.position="none")
dev.off()



# col.pal = c(
#   rgb(44/255, 20/255, 83/255),
#   rgb(164/255, 138/255, 211/255),
#   rgb(28/255, 197/255, 254/255),
#   rgb(111/255, 199/255, 207/255),
#   rgb(251/255, 162/255, 125/255),
#   rgb(251/255, 125/255, 128/255)
# )
