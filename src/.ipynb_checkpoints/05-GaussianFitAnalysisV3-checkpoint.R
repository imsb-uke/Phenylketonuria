## Gaussian fit analysis v3
rm(list = ls())
setwd("~/Desktop/R_Root/Phenylketonuria/")


qc = read.csv("Data/QC_range_by_rmse.csv")

bad_samples = c(
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

col = rep("green", 135)
col[qc$X %in% bad_samples] = "red"
plot(qc$rmse, qc$avg_range, pch = 19, col = col)

boxplot(qc$rmse)
boxplot(qc$avg_range)



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


