## Gaussian fit analysis v3
rm(list = ls())
setwd("~/Desktop/R_Root/Phenylketonuria/")


data = read.csv("Data/extracted_features_v3.csv")
data$name = paste0(data$genotype, "_", data$experiment)
data$rmse = sqrt(data$mse + .000001)

qc = read.csv("Data/qc_var_v1.csv")
data$var1 = qc$var1
data$var2 = qc$var2

bad_samples1 = c(
  "A300S-R408W_REP8",
  "A300S-IVS10_REP10",
  "D222X-R261Q_RUS7",
  "I65T-IVS10_REP4",
  "I65T-IVS12_REP4",
  "I65T-R408W_REP9_",
  "IVS10-E390G_RUS7",
  "P281L-A300S_REP6",
  "R252W-E390G_RUS7",
  "R261Q-E280K_RUS5",
  "R261Q-I306V_RUS5",
  "R261Q-R408W_REP3",
  "V245A-R408W_REP6"
)

bad_samples2 = c(
  "R408W-R408W_REP3", 
  "R408W-Y414C_REP4",  
  "I65T-IVS12_REP4", 
  "L48S-L48S_REP6", 
  "IVS10-E390G_REP11",
  "V388M-R408W_REP12", 
  "R261Q-E280K_REP12", 
  "E280K-I306V_RUS6", 
  "F39L-IVS12_UK1"
)

non_responders = c(
  "IVS12-1G",     #>A/IVS12+1G>A
  "ex5del-R408W",
  "R158Q-R158Q",
  "R261X-R408W",
  "R158Q-P281L",
  "R243X-R243X",
  "R111X-R408W",
  "R243X-R408W",
  "R158Q-R408W",
  "R252W-R252W",
  "R252W-R408W",
  "R408W-R408W",
  "R408W-IVS12",  #+1G>A
  "IVS10-11G",    #>A/R408W
  "IVS10-11G",    #>A/IVS10-11G>A
  "P281L-P281L",
  "P281L-R408W",
  "P281L-IVS10" #-11G>A
)


col = rep("green", nrow(data))
col[data$name %in% bad_samples1] = "orange"  # bad looking
col[data$name %in% bad_samples2] = "yellow"  # kicked-out reps
col[data$genotype %in% non_responders] = NA
col[data$genotype %in% "WT"] = NA  # WT

plot(data$rmse, data$n_peaks, pch = 19, col = col)
plotext = ((data$rmse > .2) | (data$n_peaks > 8)) & !(data$genotype %in% non_responders) & !(data$genotype %in% "WT")
text(data$rmse[plotext], 
     data$n_peaks[plotext] + runif(sum(plotext),-.5, .5), 
     data$genotype[plotext],
     cex = .5)

plot(data$rmse, data$var1, pch = 19, col = col, ylim = c(0, .6))
plotext = (data$var1 > .19) & !(data$genotype %in% non_responders) & !(data$genotype %in% "WT")
text(data$rmse[plotext], 
     data$var1[plotext] + runif(sum(plotext),-.01, .01), 
     data$genotype[plotext],
     cex = .7)

plot(data$rmse, data$var2, pch = 19, col = col)

plot(data$Max_x, data$Max_y, pch = 19, col = col)
plot(data$rmse, data$Max, pch = 19, col = col)
plot(data$rmse, data$var2, pch = 19, col = col)

library(plot3D)
col_num = rep(1, nrow(data))
col_num[data$genotype %in% non_responders] = 0
col_num[data$genotype %in% "WT"] = 2
col.pal = grDevices::rainbow(length(unique(col_num)))
scatter3D(exp(data$Max_x), exp(data$Max_y), data$Max, pch = 19, cex = .8,
          theta = 30, phi = 10, box = TRUE, colvar = col_num, col = col.pal)



library(ConsensusClustering)
X = data[, c("Max_x", "Max_y", "Max", "rmse")]
X$Max_x = exp(X$Max_x)
X$Max_y = exp(X$Max_y)
Adj = adj_mat(X, method = "euclidian")
clusters = pam_clust_from_adj_mat(Adj, k = 3, alpha = 1, adj.conv = FALSE)
col.pal = grDevices::rainbow(length(unique(clusters)))
scatter3D(data$Max_x, data$Max_y, data$rmse, pch = 19, cex = .8,
          theta = 30, phi = 10, box = TRUE, colvar = clusters, col = col.pal)



CM = consensus_matrix(Adj, max.cluster = 8, resample.ratio = 0.7, 
                      max.itter = 100, clustering.method = "pam")

Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, type = "b")

Kopt = Scores[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))
Kopt = 7
pheatmap::pheatmap(CM[[Kopt]])

clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
col.pal = grDevices::rainbow(length(unique(clusters)))
scatter3D(data$Max_x, data$Max_y, data$rmse, pch = 19, cex = .8,
          theta = 30, phi = 10, box = TRUE, colvar = clusters, col = col.pal)

col.pal = grDevices::rainbow(length(unique(clusters)))
scatter3D(data$Max_x, data$Max_y, data$rmse, pch = 19, cex = .8,
          theta = 30, phi = 10, box = TRUE, colvar = clusters, col = col.pal)


non_responders = (data$Max_x == max(data$Max_x)) & (data$Max_y == max(data$Max_y))
sum(non_responders)
data$non_responders = non_responders



boxplot(data$rmse)
boxplot(data$avg_range)
boxplot(data$n_peaks)




# mse
R408W-IVS12_REP4
I65T-IVS10_REP4

# var
H107P-R261Q_REP4
I65T-IVS12_REP11
I65T-IVS12_REP4

A104D-K320N_REP7
IVS10-V388M_REP11
R261Q-R261Q_REP7
P281L-IVS10_REP7


