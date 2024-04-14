# Step3. Clustering the landscapes using Consensus Clustering

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/")

library(ConsensusClustering)
library(plot3D)

data = read.csv("Data/data_processed_2.csv")

# Remove WT and non-responders
data = data[data$genotype != "WT",]
data = data[data$response != 0,]
data = data[data$genotype_exp != "L348V-R408W|REP5", ]
dim(data)

# Feature extraction
X = data.frame(
  genotype_exp = data$genotype_exp,
  x = data$Max_x,
  y = data$Max_y,
  z = log(data$Max),
  sx = log(data$s_x),
  sy = log(data$s_y),
  rmse = data$rmse
)

# Clustering
X = X[(X$x < 1400) & (X$y < 160),]
X_scaled = data.frame(scale(X[,-1]))
X_scaled = X_scaled[,c("x", "y", "z")]

Adj = adj_mat(X_scaled, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.8, max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt));  # Kopt = 6
pheatmap::pheatmap(CM[[Kopt]])

clusters = clusters = pam_clust_from_adj_mat(Adj, k = Kopt, alpha = 1, adj.conv = FALSE)
X$clusters = clusters

col.pal = grDevices::rainbow(Kopt)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = clusters, col = col.pal)


# Save
data_with_cluster = merge(data, X[,c("genotype_exp", "clusters")], by = "genotype_exp", all.x = TRUE)
data_with_cluster$clusters[is.na(data_with_cluster$clusters)] = 0
write.csv(data_with_cluster[,-c(2)], "Data/Clustering_Result_final_v2.csv")

col.pal = grDevices::rainbow(Kopt+1)
scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max, 
          pch = 19, cex = 1, theta = 10, phi = 10, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max, 
          pch = 19, cex = 1, theta = 0, phi = 90, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)


###############
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
XX$genotype_exp = paste0(XX$genotype, "|", XX$experiment)
table = merge(table, XX, by = "genotype_exp", suffixes=c("", ".y"))
}

XX = table[,c("genotype", "experiment", "x", "y", "z", "response_rate_1", "clusters")]
c = "clusters"
boxplot(
  XX$response_rate_1[XX[, c] == 1],
  XX$response_rate_1[XX[, c] == 2],
  XX$response_rate_1[XX[, c] == 3],
  XX$response_rate_1[XX[, c] == 4],
  # XX$response_rate_1[XX[, c] == 5],
  # XX$response_rate_1[XX[, c] == 6],
  names = c(1, 2, 3, 4),
  main = paste0("Response in clusters:", c)
)


