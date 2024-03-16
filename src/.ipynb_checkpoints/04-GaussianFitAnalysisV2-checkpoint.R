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
features = read.csv("Data/extracted_features_v2.csv")
# features$genotype_rev = sapply(features$genotype, function(x){
#   if (x != "WT"){
#     a = strsplit(x, "-")[[1]]
#     y = paste0(a[2], "-", a[1])
#   }else{
#     y = "WT"
#   }
#   return(y)
#   })

# Merge the two tables
# table$genotype_rev = table$genotype
table = merge(table, features, by="genotype")
# table$genotype_rev = table$genotype_rev.x
# table = merge(table, features, by="genotype_rev")

# Add response information
# Remove WT
table = table[table$genotype != "WT",]

# Normalization
data = data.frame(
  genotype = table$genotype,
  xi = table$Phe.at.max..µmol.L,
  yi = table$BH4.at.max..µmol.L,
  xl = table$Max_x,
  yl = table$Max_y,
  x = exp(table$Max_x),
  y = exp(table$Max_y),
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

data$xi[is.na(data$xi)] = 2500
data$yi[is.na(data$yi)] = 250

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


# Hand-made clusters
HandClust = read_excel("Data/20240221_Clustering_Result_v2_all_PG.XLSX")
data = merge(data, HandClust, by = "genotype")

col.pal = grDevices::rainbow(6)
scatter3D(data$x, data$y, data$z , pch = 19, cex = 1, main="ByHand clusters",
          theta = 30, phi = 10, box = TRUE, colvar = data$byHand, col = col.pal)
scatter3D(data$x, data$y, data$z , pch = 19, cex = 1, main="ByHand clusters",
          theta = 0, phi = 90, box = TRUE, colvar = data$byHand, col = col.pal)
scatter3D(data$sx, data$sy, data$mses , pch = 19, cex = 1, main="ByHand clusters",
          theta = 30, phi = 10, box = TRUE, colvar = data$byHand, col = col.pal)

# Clustering
library(ConsensusClustering)

X = data[,c("genotype", "x", "y", "z", "sx", "sy", "mses")]
X = data[,c("genotype", "x", "y", "z", "sx", "sy")]
X = data[,c("genotype", "x", "y", "z", "sxl", "syl", "mses")]
X = data[,c("genotype", "x", "y", "z", "sxl", "syl")]
X = data[,c("genotype", "x", "y", "z", "sxl", "syl", "mses")]

X = X[X$x < max(X$x),]

Adj = adj_mat(X[,-1], method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, 
                      max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
# Kopt = 6
pheatmap::pheatmap(CM[[Kopt]])

clusters = hir_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
# clusters = hir_clust_from_adj_mat(Adj, k = Kopt, alpha = 1, adj.conv = FALSE)
#clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)

X$clusters = clusters

col.pal = grDevices::rainbow(Kopt)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 10, phi = 10, box = TRUE, colvar = X$clusters, col = col.pal)

scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters of all samples", 
          theta = 0, phi = 90, box = TRUE, colvar = X$clusters, col = col.pal)

##
HandClust = HandClust[,c("genotype", "byHand")]
# XX = X
XX$clusters5 = clusters
toSave = merge(XX, HandClust, by = "genotype")

write.csv(toSave, "Data/Clustering_Result_v3.csv")


###
# Read the clustering results and then do prediction for each
Clusters = read.csv("Data/Clustering_Result_v3.csv")[,c(2, 9:13)]
HandClust = read_excel("Data/20240221_Clustering_Result_v2_all_PG.XLSX")
HandClust = HandClust[,c("genotype", "byHand")]
Clusters = merge(HandClust, Clusters, by = "genotype")


X = data[,c("genotype", "x", "y", "z", "sx", "sy")]
x0 = X[X$x >= max(X$x),1]
X0 = data.frame(genotype = x0, 
                clusters1 = rep(0, length(x0)),
                clusters2 = rep(0, length(x0)),
                clusters3 = rep(0, length(x0)),
                clusters4 = rep(0, length(x0)),
                clusters5 = rep(0, length(x0)),
                byHand = rep(0, length(x0))
                )

Clusters = rbind(Clusters, X0)

Clusters_with_response = merge(Clusters, data[,c("genotype", "response_rate_1")])

Clusters_with_response$experiment = sapply(Clusters_with_response$genotype, function(x){strsplit(x, "-")[[1]][1]})
Clusters_with_response$response_rate_rm_I65T = Clusters_with_response$response_rate_1
Clusters_with_response$response_rate_rm_I65T[Clusters_with_response$experiment == "I65T"] = NA

write.csv(Clusters_with_response, "Data/Clustering_Result_v4.csv")

c = "clusters1"
boxplot(
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 0],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 3],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 1],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 4],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 2],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 5],
  names = c(0, 3, 1, 4, 2, 5),
  main = paste0("Response in clusters:", c)
)

c = "clusters2"
boxplot(
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 0],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 3],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 2],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 4],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 1],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 5],
  names = c(0, 3, 2, 4, 1, 5),
  main = paste0("Response in clusters:", c)
)

c = "clusters3"
boxplot(
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 0],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 3],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 4],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 2],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 1],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 5],
  names = c(0, 3, 4, 2, 1, 5),
  main = paste0("Response in clusters:", c)
)

c = "clusters4"
boxplot(
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 0],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 4],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 3],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 2],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 1],
  Clusters_with_response$response_rate_rm_I65T[Clusters_with_response[, c] == 5],
  names = c(0, 4, 3, 2, 1, 5),
  main = paste0("Response in clusters:", c)
)

c = "byHand"
boxplot(
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 6],
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 4],
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 3],
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 2],
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 1],
  Clusters_with_response$response_rate_1[Clusters_with_response[, c] == 5],
  names = c(0, 4, 3, 2, 1, 5),
  main = paste0("Response in clusters:", c)
)

X = data[,c("genotype", "x", "y", "z", "sx", "sy")]
X = merge(X, Clusters_with_response, by = "genotype")
col.pal = grDevices::rainbow(6)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters 4", 
          theta = 10, phi = 10, box = TRUE, colvar = X$clusters4, col = col.pal)
scatter3D(X$x, X$y, X$z , pch = 19, cex = 1, main = "Clusters 4", 
          theta = 0, phi = 90, box = TRUE, colvar = X$clusters4, col = col.pal)

