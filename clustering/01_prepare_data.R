# Step 1. Read and process the tables

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria")

library(readxl)
library(MASS)

# Read the gm.py output
data = read.csv("gm_output/features/extracted_features.csv")

# Read final qc description excel file
qc_table = read_excel("Data/20240321_features_and_qc_PG.XLSX")
qc_table = data.frame(qc_table[, c("genotype", "experiment", "final decision")])

# Read the table with response to treatment
#...

# Merge the two tables
data$genotype_exp = paste0(data$genotype, "|", data$experiment)
qc_table$genotype_exp = paste0(qc_table$genotype, "|", qc_table$experiment)
data = merge(data, qc_table, by = "genotype_exp", suffixes=c("", ".y"))

# RMSE≥0.25, n_peaks≥8, var≥0.25
hist(data$rmse, 20)
hist(data$n_peaks)
hist(data$variation, 20)

2 / fitdistr(data$rmse, "exponential")$estimate
2 / fitdistr(data$n_peaks, "exponential")$estimate
2 / fitdistr(data$variation, "exponential")$estimate

# Selected QC passed samples
idx = !is.na(data$final.decision)
print(paste0("Number of samples passed QC: ", sum(idx)))
data = data[idx, ]

# Exclude "WTs" and "Complete non-responders"
## Hint: In the current version of the qc_table, all "WT" samples are marked to be omitted.
# data_wt = data[data$genotype == "WT", ]

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
  "P281L-IVS10", #-11G>A
  "IVS10-IS10"
)



data$response_available = rep(1, nrow(data))
data$response_available[data$genotype %in% non_responders] = 0

# Remove extra cols
data = data[, -c(18, 19)]

# Save as csv
write.csv(data, "Data/data_processed_6.csv")

