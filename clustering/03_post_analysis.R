# Step3. treatment response

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
  
  response_rate_1 = table$Yes... + table$Slow...
  response_rate_1[is.na(response_rate_1)] = -100
  response_rate_1 = response_rate_1 / 100
  table$response_rate_1 = response_rate_1
  table$response_rate_1[response_rate_1 == -1] = NA
  table$response_rate_1[table$BH4.response..Nr.of..tested.patients < 10] = NA
  exclude_genotypes = c("F39L-IVS12", "I65T-I65T", "I65T-IVS12", "I65-R261Q", "I65T-R408W", "F39L-R408W", "I65T-IVS10")
  table$response_rate_1[table$genotype %in% exclude_genotypes] = NA
  
  # Merge
  # Merge the two tables
  table$genotype_exp = paste0(table$genotype, "|", table$EXPERIMENT)
  data$genotype_exp = paste0(data$genotype, "|", data$experiment)
  table = merge(table, data, by = "genotype_exp", suffixes=c("", ".y"))
}

X = table[,c("genotype", "experiment", "Max_theory", "Max_x", "Max_y", "response_rate_1", "clusters")]
X_response = X[!is.na(X$response_rate_1),]

# Scatter plot
order = c(0, 2, 5, 1, 3, 4)
col.pal = c("red", "orange", "green", "skyblue", "blue", "magenta")
# col.pal = grDevices::rainbow(length(unique(X$clusters)))
X$clusters[X$clusters == 6] = 2
# X = X[X$clusters > 0,]


pdf(file = "Clusters.pdf", width = 5, height = 5)

scatter3D(X$Max_x, X$Max_y, log(X$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 0, box = TRUE, colvar = X$clusters, col = col.pal)

scatter3D(X$Max_x, X$Max_y, log(X$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = X$clusters, col = col.pal)

XX = X[X$clusters > 0,]
scatter3D(XX$Max_x, XX$Max_y, log(XX$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 0, box = TRUE, colvar = XX$clusters, col = col.pal)

scatter3D(XX$Max_x, XX$Max_y, log(XX$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = XX$clusters, col = col.pal)

dev.off()

# Response box plot
c = "clusters"
boxplot(
  X_response$response_rate_1[X_response[, c] == 0],
  X_response$response_rate_1[X_response[, c] == 2],
  X_response$response_rate_1[X_response[, c] == 5],
  X_response$response_rate_1[X_response[, c] == 1],
  X_response$response_rate_1[X_response[, c] == 3],
  X_response$response_rate_1[X_response[, c] == 4],
  names = order
)

X_response_no6 = X_response[X_response$clusters < 6,]
X_response_no6$clusters_factor = factor(X_response_no6$clusters, 
                                        levels = order)

pdf(file = "Response.pdf", width = 5, height = 5)

ggplot(X_response_no6, aes(x=clusters_factor, y=response_rate_1, color=clusters_factor)) + 
  geom_boxplot(fill="gray95", size = .6) + 
  geom_point() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values = col.pal[order+1]) + 
  theme_classic() +
  theme(legend.position="none")
dev.off()






