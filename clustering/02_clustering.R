# Step2. Clustering the landscapes using Consensus Clustering

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/")

library(ConsensusClustering)
library(plot3D)

data_all = read.csv("Data/data_processed_6.csv")

# Remove WT and non-responders
data_all = data_all[data_all$genotype != "WT",]
data = data_all[data_all$response != 0,]
# data = data[data$genotype_exp != "L348V-R408W|REP5", ]
dim(data)

# Feature extraction
X = data.frame(
  genotype_exp = data$genotype_exp,
  x = data$Max_x,
  y = data$Max_y,
  z = log(data$Max_theory),
  sx = log(data$s_x),
  sy = log(data$s_y),
  rmse = data$rmse
)

# Clustering
filter = ((X$x > 1400) | (X$y > 170)) & (X$z < log(.5))
plot(X[,2], X[,3], pch = 19, col = ifelse(filter, "green", "red"))
X = X[!filter,]

# X$z = exp(X$z)
X_scaled = data.frame(scale(X[,-1]))
X_scaled = X_scaled[,c("x", "y", "z")]
Adj = adj_mat(X_scaled, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b", pch = 19)

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
pheatmap::pheatmap(CM[[Kopt]])

clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = TRUE)
X$clusters = clusters

col.pal = grDevices::rainbow(Kopt + 1)
scatter3D(X$x, X$y, (X$z) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 0, box = TRUE, colvar = clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = clusters, col = col.pal)

# Save
data_with_cluster = merge(data_all, X[,c("genotype_exp", "clusters")], by = "genotype_exp", all.x = TRUE)
data_with_cluster$clusters[is.na(data_with_cluster$clusters)] = 0

col.pal = grDevices::rainbow(Kopt+1)
scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 10, phi = 10, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 0, phi = 90, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

write.csv(data_with_cluster[,-c(2)], "Data/Clustering_Result_final_v6.csv")
