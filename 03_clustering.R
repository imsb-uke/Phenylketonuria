# Step3. Clustering the landscapes using Consensus Clustering

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/")

library(ConsensusClustering)
library(plot3D)

data = read.csv("Data/data_processed.csv")

# Remove WT and non-responders
data = data[data$genotype != "WT",]
data = data[data$response != 0,]
dim(data)

# Feature extraction
X = data.frame(
  x = exp(data$Max_x),
  y = exp(data$Max_y),
  z = log(data$Max),
  sx = log(data$s_x),
  sy = log(data$s_y),
  rmse = data$rmse
)

# Clustering
## Clustering without scaling
Adj = adj_mat(X, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt)); # Kopt = 6
pheatmap::pheatmap(CM[[Kopt]])

clusters = clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
Clusters = clusters

col.pal = grDevices::rainbow(Kopt)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = clusters, col = col.pal)


## Clustering with scaling
X_scaled = scale(X)

Adj = adj_mat(X_scaled, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt)); # Kopt = 6
pheatmap::pheatmap(CM[[Kopt]])

clusters = clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
Clusters = cbind(Clusters, clusters)

Clusters = data.frame(Clusters)
colnames(Clusters) = c("not_scaled", "scaled")

col.pal = grDevices::rainbow(Kopt)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = clusters, col = col.pal)


# Save
XX = cbind(data[,c("genotype", "experiment")], X, Clusters)
write.csv(XX, "Data/Clustering_Result_v5.csv")

