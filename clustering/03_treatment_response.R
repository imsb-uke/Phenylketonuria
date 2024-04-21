# Step3. treatment response

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/")
library(plot3D)

data = read.csv("Data/Clustering_Result_final_v5.csv")

{
  # Read and prepare "landscapes_overview.xlsx"
  library(readxl)
  table = read_excel("Data/20240207_landscapes_overview.xlsx")
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
  
  response_rate_1 = table$Yes...
  response_rate_1[is.na(response_rate_1)] = -100
  response_rate_1 = response_rate_1 / 100
  table$response_rate_1 = response_rate_1
  table$response_rate_1[response_rate_1 == -1] = NA
  
  # Merge
  # Merge the two tables
  table$genotype_exp = paste0(table$genotype, "|", table$EXPERIMENT)
  data$genotype_exp = paste0(data$genotype, "|", data$experiment)
  table = merge(table, data, by = "genotype_exp", suffixes=c("", ".y"))
}

X = table[,c("genotype", "experiment", "Max_theory", "Max_x", "Max_y", "response_rate_1", "clusters")]
X_response = X[!is.na(X$response_rate_1),]

col.pal = grDevices::rainbow(length(unique(X$clusters)))
scatter3D(X$Max_x, X$Max_y, log(X$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 0, box = TRUE, colvar = X$clusters, col = col.pal)

scatter3D(X$Max_x, X$Max_y, log(X$Max_theory) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = X$clusters, col = col.pal)


c = "clusters"
boxplot(
  X_response$response_rate_1[X_response[, c] == 0],
  X_response$response_rate_1[X_response[, c] == 2],
  X_response$response_rate_1[X_response[, c] == 1],
  X_response$response_rate_1[X_response[, c] == 3],
  X_response$response_rate_1[X_response[, c] == 4],
  names = c(0, 2, 1, 3, 4),
  main = paste0("Response in clusters:", c)
)

c = "clusters"
boxplot(
  X_response$response_rate_1[X_response[, c] == 0],
  X_response$response_rate_1[X_response[, c] == 2],
  X_response$response_rate_1[X_response[, c] == 4],
  X_response$response_rate_1[X_response[, c] == 1],
  X_response$response_rate_1[X_response[, c] == 3],
  names = c(0, 2, 4, 1, 3),
  main = paste0("Response in clusters:", c)
)



