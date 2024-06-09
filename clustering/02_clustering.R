# Step2. Clustering the landscapes using Consensus Clustering

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/clustering/")

library(ConsensusClustering)
library(plot3D)

data_all = read.csv("data/data_processed.csv")

# Remove WT and non-responders
data_all = data_all[data_all$genotype != "WT",]
data = data_all[data_all$response != 0,]
data = data_all
# data = data[data$genotype_exp != "L348V-R408W|REP5", ]
dim(data)

# Feature extraction
X = data.frame(
  genotype_exp = data$genotype_exp,
  x = data$Max_x,
  y = data$Max_y,
  z = data$Max_theory,
  sx = log(data$s_x),
  sy = log(data$s_y),
  rmse = data$rmse
)

# Clustering
filter = ((X$x > 1400) | (X$y > 190)) & (X$z < .5)
plot(X[,2], X[,3], pch = 19, col = ifelse(filter, "green", "red"))
X = X[!filter,]

# X$z = log(X$z)
X_scaled = data.frame(scale(X[,-1]))
X_scaled = X_scaled[,c("x", "y", "z")]
Adj = adj_mat(X_scaled, method = "euclidian")
# CM = consensus_matrix(Adj, max.cluster = 10, resample.ratio = 0.7, max.itter = 500, clustering.method = "pam")
# saveRDS(CM, "data/CM.rds")
CM = readRDS("data/CM.rds")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b", pch = 19)

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
pheatmap::pheatmap(CM[[Kopt]], border_color = NA, show_rownames = FALSE, show_colnames = FALSE)

clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = TRUE)
X$clusters = clusters

col.pal = grDevices::rainbow(Kopt + 2)

scatter3D(X$x, X$y, exp(X$z) , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 0, box = TRUE, colvar = clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = clusters, col = col.pal)

# Merge
data_with_cluster = merge(data_all, X[,c("genotype_exp", "clusters")], by = "genotype_exp", all.x = TRUE)
data_with_cluster$clusters[is.na(data_with_cluster$clusters)] = 0

# Get averages of clusters
for (i in 1:5){
  avg_x_i = mean(data_with_cluster$Max_x[data_with_cluster$clusters == i])
  avg_y_i = mean(data_with_cluster$Max_y[data_with_cluster$clusters == i])
  print(paste0("cluster ", i, ": aveage_x = ", round(avg_x_i, 2), ", aveage_y = ", round(avg_y_i, 2)))
}


# Plot
col.pal = grDevices::rainbow(Kopt+1)
scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 10, phi = 10, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 0, phi = 90, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

###
# Put those three sampples as cluster 6
data_with_cluster$clusters[data_with_cluster$Max_y < 190 & data_with_cluster$clusters == 0] = 6
data_with_cluster$clusters[data_with_cluster$Max_x < 1400 & data_with_cluster$clusters == 0] = 6

col.pal = grDevices::rainbow(Kopt+2)
scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 10, phi = 10, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

scatter3D(data_with_cluster$Max_x, data_with_cluster$Max_y, data_with_cluster$Max_theory, 
          pch = 19, cex = 1, theta = 0, phi = 90, box = TRUE, 
          colvar = data_with_cluster$clusters, col = col.pal)

###
write.csv(data_with_cluster[,-c(2)], "data/clustering_result.csv")



