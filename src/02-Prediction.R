# Prediction
rm(list = ls())
setwd("~/Desktop/R_Root/Polina/")

classification_metrics = function(y_hat, y){
  
  tab = table(y_hat, y)
  
  tp = tab[1,1]
  fp = tab[1,2]
  fn = tab[2,1]
  tn = tab[2,2]
  
  acc = (tn + tp) / (tn + fp +fn + tp)
  sen = tp / (tp + fn)
  spc = tn / (tn + fp)
  prc = tp / (tp + fp)
  
  r = c(acc, sen, spc, prc)
  names(r) = c("acc", "sen", "spc", "prc")
  
  return(r)
}


### Read data
data = readRDS("Data/table1_v2.rds")
data_naOmit = na.omit(data)

X = data_naOmit[,1:3]
y = data_naOmit$response_rate_1
y[y > .5] = 1
y[y < .5] = 0
y = as.factor(y)

## Cross validation loop
N_data = nrow(X)
train_rate = 0.8
N_rep = 100

Result = c()
for (rep in 1:N_rep){
  
  index = sample(N_data)
  X_train = X[index[1:floor(train_rate * N_data)], ]
  X_test = X[-index[1:floor(train_rate * N_data)], ]
  y_train = y[index[1:floor(train_rate * N_data)]]
  y_test = y[-index[1:floor(train_rate * N_data)]]
  
  # m = apply(X_train, 2, mean)
  # s = apply(X_train, 2, sd)
  m = apply(X, 2, mean)
  s = apply(X, 2, sd)
  
  for (i in 1:ncol(X)){
    X_train[,i] = (X_train[,i] - m[i]) / s[i]
    X_test[,i] = (X_test[,i] - m[i]) / s[i]
  }
  
  data_train = cbind(X_train, y_train)
  model = glm(y_train ~ ., data = data_train, family = "binomial")
  y_hat = predict(model, X_test,  type = "response")
  y_hat = ifelse(y_hat > 0.5, 1, 0)
  
  Result = rbind(Result, classification_metrics(y_hat, y_test))

}
apply(Result, 2, mean)


## Prediction for NA samples
X = data[,1:3]
y = data$response_rate_1

m = apply(X, 2, mean)
s = apply(X, 2, sd)
for (i in 1:ncol(X))
  X[,i] = (X[,i] - m[i]) / s[i]

X_train = X[!is.na(y),]
X_test = X[is.na(y),]
y_train = y[!is.na(y)]

y_train[y_train > .5] = 1
y_train[y_train < .5] = 0
y_train = as.factor(y_train)

data_train = cbind(X_train, y_train)
model = glm(y_train ~ ., data = data_train, family = "binomial")

y_hat = predict(model, X_test,  type = "response")
y_hat = as.numeric(ifelse(y_hat > 0.5, 1, 0))

data_test = cbind(X_test, y_hat)

colnames(data_train) = colnames(data_test) = c("x", "y", "z", "response")
data_total = rbind(data_train, data_test)
data_total$response = as.numeric(data_total$response)

library(plot3D)

col.pal = colorRamps::green2red(100)[100:1]
scatter3D(data_total$x, data_total$y, data_total$z , pch = 19, cex = 1, main="All",
          theta = 30, phi = 10, box = TRUE, colvar = as.numeric(data_total$response), col = col.pal)

scatter3D(data_train$x, data_train$y, data_train$z , pch = 19, cex = 1, main="Train",
          theta = 30, phi = 10, box = TRUE, colvar = as.numeric(data_train$response), col = col.pal)

scatter3D(data_test$x, data_test$y, data_test$z , pch = 19, cex = 1, main="Test",
          theta = 30, phi = 10, box = TRUE, colvar = as.numeric(data_test$response), col = col.pal)


#####
data$response_prediction = -1
data$response_prediction[is.na(data$response_rate_1)] = data_test$response
data$response_clustering = ifelse(data$clusters == 2, 0, 1)
saveRDS(data, file = "Data/table1_v3.rds")

data_NA = data[is.na(data$response_rate_1), ]
sum(data_NA$response_prediction == data_NA$response_clustering)

write.csv(data_NA, file = "Data/Predition_for_missing_responses.csv")
