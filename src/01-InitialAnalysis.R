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


# Deal with NA values in landscape summary
# 1. Remove NA in max residual activity (z) value
dim(table)
table = table[!is.na(table[,3]),]
table = table[!is.na(table[,4]),]
table = table[!is.na(table[,5]),]
dim(table)
# 2. Cap NA in x,y value
# max_phe_at_max = max(table$Phe.at.max..µmol.L, na.rm = TRUE)
# max_bh4_at_max = max(table$BH4.at.max..µmol.L, na.rm = TRUE)
# table$Phe.at.max..µmol.L[is.na(table$Phe.at.max..µmol.L)] = max_phe_at_max + 1000
# table$BH4.at.max..µmol.L[is.na(table$BH4.at.max..µmol.L)] = max_bh4_at_max + 100

# Normalization
hist(table[,3], 20, main = "Histogram of max residual activity")
hist(table[,4], 20, main = "Histogram of max _ x")
hist(table[,5], 20, main = "Histogram of max _ y")

hist(log(table[,3]), 20, main = "Histogram of log max residual activity")
hist(log(table[,4]), 20, main = "Histogram of log max _ x")
hist(log(table[,5]), 20, main = "Histogram of log max _ y")

data = data.frame(
  genotype = table[,1],
  x = log(table[,4]),
  y = log(table[,5]),
  z = log(table[,3])
)

plot(data$x, data$y)
plot(data$x, data$z)
plot(data$y, data$z)

# Response visualization
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

# Phe level
severity_level_1 = table$Classic_PKU...1200µmol.L.Phe....
severity_level_1[is.na(severity_level_1)] = -100
severity_level_1 = severity_level_1 / 100
data$severity_level_1 = severity_level_1
data$severity_level_1[severity_level_1 == -1] = NA

severity_level_2 = table$MHPA...600µmol.L.Phe....
severity_level_2[is.na(severity_level_2)] = 200
severity_level_2 = 1 - severity_level_2 / 100
data$severity_level_2 = severity_level_2
data$severity_level_2[severity_level_2 == -1] = NA

# Some visualization
data_naOmit = na.omit(data)

color1 = ifelse(data_naOmit$response_rate_1 > .5, "blue", "red")
color2 = ifelse(data_naOmit$response_rate_2 > .5, "blue", "red")

plot(data_naOmit$x, data_naOmit$y, col = color1)
plot(data_naOmit$x, data_naOmit$z, col = color1)
plot(data_naOmit$y, data_naOmit$z, col = color1)

# saveRDS(data, file = "Data/table1_v2.rds")

# Plot 3D
# library(plot3D)

col.pal = colorRamps::green2red(100)[100:1]

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 0, phi = 90, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 0, phi = 0, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Response to treatment 1",
          theta = 90, phi = 0, box = TRUE, colvar = data_naOmit$response_rate_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Severity level 1",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$severity_level_1, col = col.pal)

scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main="Severity level 2",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$severity_level_2, col = col.pal)

# Clustering
library(ConsensusClustering)

X = data[,1:3]
Adj = adj_mat(X, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
Kopt = 6
pheatmap::pheatmap(CM[[Kopt]])

clusters = hir_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
data$clustersK6 = clusters



# remove NAs in response to treatment
data_naOmit = na.omit(data)

col.pal = grDevices::rainbow(Kopt)
scatter3D(data_naOmit$x, data_naOmit$y, data_naOmit$z , pch = 19, cex = 1, main = "Clusters for samples with \n response to treatment",
          theta = 30, phi = 10, box = TRUE, colvar = data_naOmit$clusters, col = col.pal)

# Boxplot
boxplot(
  data_naOmit$response_rate_1[data_naOmit$cluster==1],
  data_naOmit$response_rate_1[data_naOmit$cluster==2],
  data_naOmit$response_rate_1[data_naOmit$cluster==3],
  main = "Clusters VS Response to treatment 1"
  )

boxplot(
  data_naOmit$response_rate_2[data_naOmit$cluster==1],
  data_naOmit$response_rate_2[data_naOmit$cluster==2],
  data_naOmit$response_rate_2[data_naOmit$cluster==3],
  main = "Clusters VS Response to treatment 2"
)

boxplot(
  data_naOmit$severity_level_1[data_naOmit$cluster==1],
  data_naOmit$severity_level_1[data_naOmit$cluster==2],
  data_naOmit$severity_level_1[data_naOmit$cluster==3],
  main = "Clusters VS Severity level 1"
)

boxplot(
  data_naOmit$severity_level_2[data_naOmit$cluster==1],
  data_naOmit$severity_level_2[data_naOmit$cluster==2],
  data_naOmit$severity_level_2[data_naOmit$cluster==3],
  main = "Clusters VS Severity level 2"
)


table$clustersK3 = data$clustersK3
table$clustersK6 = data$clustersK6

col.pal = grDevices::rainbow(Kopt)
scatter3D(log(table[,4]), log(table[,5]), log(table[,3]) , pch = 19, cex = 1, main = "Clusters for samples with \n response to treatment",
          theta = 30, phi = 10, box = TRUE, colvar = table$clustersK6, col = col.pal)

write.csv(table, "Data/Table_clustering_v1.csv")
