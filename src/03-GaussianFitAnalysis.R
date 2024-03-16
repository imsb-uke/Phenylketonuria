# Analysis of Gaussian fitted features

# Modelling_FeatureExtraction script has been used to extract Gaussian fitted features.
# Compared to the old "01-InitialAnalysis", we have new contributions:
# 1. x, y, z are now based on Gaussian modellig
# 2. Three new features are added
# 3. Complete responders (WT samples) and complete non-responders (NA values) have also been added

rm(list = ls())
setwd("~/Desktop/R_Root/Polina/")

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

# Read Gaussian fitted features
features = read.csv("Data/extracted_features_v1.csv")

# Merge the two tables
table = merge(table, features, by="genotype")

# Add response information
# Remove WT
table = table[table$genotype != "WT",]


# Normalization
hist(table$Max, 20, main = "Histogram of max residual activity")
hist(table$Max_x, 20, main = "Histogram of max _ x")
hist(table$Max_y, 20, main = "Histogram of max _ y")
hist(table$s_x, 20, main = "Histogram of max _ y")
hist(table$s_y, 20, main = "Histogram of max _ y")
hist(table$mse, 20, main = "Histogram of max _ y")
hist(table$s_x * table$s_y, 20, main = "Histogram of max _ y")

hist(log(table$Max), 20, main = "Histogram of log max residual activity")
# hist(log(table$Max_x), 20, main = "Histogram of log max _ x")
# hist(log(table$Max_y), 20, main = "Histogram of log max _ y")
hist(log(table$s_x), 20, main = "Histogram of max _ y")
hist(log(table$s_y), 20, main = "Histogram of max _ y")
hist(log(table$mse), 20, main = "Histogram of max _ y")

data = data.frame(
  genotype = table$genotype,
  x = table$Max_x,
  y = table$Max_y,
  z = log(table$Max),
  sx = table$s_x,
  sy = table$s_y,
  s = table$s_x * table$s_y,
  sxl = log(table$s_x),
  syl = log(table$s_y),
  sl = log(table$s_x * table$s_y),
  ss = sqrt(table$s_x * table$s_y),
  mse = table$mse,
  msel = log(table$mse),
  mses = sqrt(table$mse)
)

for (i in 2:ncol(data)){
  mu = mean(data[,i])
  sd = sd(data[,i])
  print(sd)
  data[,i] = (data[,i] - mu) / sd
}

# Add responses
response_rate_1 = table$Yes...
response_rate_1[is.na(response_rate_1)] = -100
response_rate_1 = response_rate_1 / 100
data$response_rate_1 = response_rate_1
data$response_rate_1[response_rate_1 == -1] = NA

response_rate_2 = table$No...
response_rate_2[is.na(response_rate_2)] = 200
response_rate_2 = 1 - response_rate_2 / 100
data$response_rate_2 = response_rate_2
data$response_rate_2[response_rate_2 == -1] = NA

# 3D plot
library(plot3D)

col.pal = colorRamps::green2red(10)
scatter3D(data$x, data$y, data$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 30, phi = 10, box = TRUE, colvar = data$mse, col = col.pal)

data_naOmit = na.omit(data)

color1 = ifelse(data_naOmit$response_rate_1 > .5, "blue", "red")
color2 = ifelse(data_naOmit$response_rate_2 > .5, "blue", "red")

col.pal = colorRamps::green2red(100)[100:1]

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 0, phi = 90, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 0, phi = 0, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 90, phi = 0, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$msel, data_naOmit$ss, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)


# Clustering
library(ConsensusClustering)

# X = data[,c(2,3,4,8)]
X = data[,c("x", "y", "z")]  # -> 4
X = data[,c("x", "y", "z", "sl")] # -> 5
X = data[,c("x", "y", "z", "s", "mse")]  # -> 4
X = data[,c("x", "y", "z", "sl", "mse")] # -> 6
X = data[,c("x", "y", "z", "ss", "mse")] # -> 4

X = data[,c("x", "y", "z", "ss", "msel")] # -> hc6, pam5
X = data[,c("x", "y", "z", "ss", "mses")] # -> hc5, pam6
X = data[,c("x", "y", "z", "sx", "sy", "mses")] # -> hc5, pam6

Adj = adj_mat(X, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 10, resample.ratio = 0.7, max.itter = 50, clustering.method = "hclust")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
Kopt = 5
pheatmap::pheatmap(CM[[Kopt]])

clusters = hir_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
# clusters = hir_clust_from_adj_mat(Adj, k = Kopt, alpha = 1, adj.conv = FALSE)
clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)

data$clusters = clusters

col.pal = grDevices::rainbow(Kopt)
scatter3D(data$x, data$y, data$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = data$clusters, col = col.pal)

scatter3D(data$x, data$y, data$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = data$clusters, col = col.pal)

scatter3D(data$msel, data$ss, data$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = data$clusters, col = col.pal)

scatter3D(data$sx, data$sy, data$z , pch = 19, cex = 1, main = "sx-sy", 
          theta = 0, phi = 90, box = TRUE, colvar = data$clusters, col = col.pal)


#
# Clusters = data.frame(genotype = data$genotype)
Clusters$C3h = clusters

clusters_hand = read.table("Data/Clusters_hand.csv", sep=";", header = TRUE)
d = merge(Clusters, clusters_hand, by="genotype")

col.pal = grDevices::rainbow(6)

scatter3D(d$x, d$y, d$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = d$cluster, col = col.pal)

scatter3D(d$x, d$y, d$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = d$cluster, col = col.pal)


# Save
# write.csv(Clusters, "Data/Clustering_Result_v2_all.csv")
# write.csv(d, "Data/Clustering_Result_v2_annotated.csv")

# Response

c = "C3h"
data$clusters = Clusters[,c]
boxplot(
  data$response_rate_1[data$cluster==1],
  data$response_rate_1[data$cluster==2],
  data$response_rate_1[data$cluster==3],
  data$response_rate_1[data$cluster==4],
  data$response_rate_1[data$cluster==5],
  data$response_rate_1[data$cluster==6],
  main = paste0("Response in clusters:", c)
)

library(aricode)
ARI(d$cluster, d$C1h)
ARI(d$cluster, d$C1p)
ARI(d$cluster, d$C2h)
ARI(d$cluster, d$C2p)
ARI(d$cluster, d$C3p)
ARI(d$cluster, d$C3h)




