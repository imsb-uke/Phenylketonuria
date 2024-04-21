# Step2. treatment response

rm(list = ls())
setwd("~/Desktop/My_Codes/Phenylketonuria/")

library(plot3D)

data = read.csv("Data/data_processed_6.csv")

col.pal = colorRamps::green2red(10)[10:1]
scatter3D(data$Max_x, data$Max_y, log(data$Max_theory) , pch = 19, cex = 1, main="Response to treatment",
          theta = 30, phi = 10, box = TRUE, colvar = data$response, col = col.pal)

plot(data$Max_x, data$Max_y, pch = 19, cex = 1, 
     col = ifelse(data$response == 0, "red", "green"),
     main="Response to treatment", xlab = "Max x", ylab = "Max y")
abline(v = 1400, h = 170)
