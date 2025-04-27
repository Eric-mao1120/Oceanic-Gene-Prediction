# Random Forest Model for phnJ Gene Prediction
# Author: [Shihai Mao]
# Date: [2024-11-22]

# ============================================================
# 1. Environment Setup
# ============================================================

setwd("D:/R language/R")

# Load libraries
library(randomForest)
library(tidyverse)
library(caret)
library(ggplot2)
library(RColorBrewer)
library(Metrics)

# ============================================================
# 2. Data Preprocessing
# ============================================================

# Load dataset
MMOX_unnormalize <- read.csv("D:/Desktop/phnJ DATA for model.csv")

# Partition data into training (80%) and testing (20%)
set.seed(21212)
train_indices <- createDataPartition(y = MMOX_unnormalize$phnJ, p = 0.8, list = FALSE)
traindata_raw <- MMOX_unnormalize[train_indices, ]
testdata_raw  <- MMOX_unnormalize[-train_indices, ]

phnJ_mean <- mean(traindata_raw$phnJ)
phnJ_sd <- sd(traindata_raw$phnJ)

# Standardize all features
traindata_raw$phnJ <- scale(traindata_raw$phnJ, center = phnJ_mean, scale = phnJ_sd)
testdata_raw$phnJ  <- scale(testdata_raw$phnJ, center = phnJ_mean, scale = phnJ_sd)

traindata <- traindata_raw
testdata <- testdata_raw

# Define the random forest model formula
form_reg <- as.formula("phnJ ~ Chla + SSS + POC + SST + Pi + N + DOP")

# ============================================================
# 3. Model Training and Cross-Validation
# ============================================================

# 10-fold cross-validation to optimize mtry
set.seed(3342)
train_control <- trainControl(method = "cv", number = 10)
tune_grid <- expand.grid(mtry = 2:6)

rf_cv_model <- train(
  form = form_reg,
  data = traindata,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 125
)

# Extract the best mtry value
best_mtry <- rf_cv_model$bestTune$mtry
cat("Best mtry value:", best_mtry, "\n")

# Hyperparameter sensitivity analysis (ntree vs mtry)
set.seed(3152)
ntrees_seq <- seq(25, 300, by = 25)  
rf_acc <- numeric(length(ntrees_seq))

for (i in seq_along(ntrees_seq)) {
  model_temp <- randomForest(
    form_reg,
    data = traindata,
    ntree = ntrees_seq[i],
    mtry = best_mtry,
    importance = TRUE
  )
  prediction_temp <- predict(model_temp, testdata)
  rf_acc[i] <- cor(testdata$phnJ, prediction_temp)
}

rf_ntree_df <- data.frame(
  ntree = ntrees_seq,
  accuracy = rf_acc
)

# Plot ntree vs accuracy
p2 <- ggplot(rf_ntree_df, aes(x = ntree, y = accuracy)) +
  geom_line(size = 1.2) +
  geom_point(size = 3,color="red3") +ylim(0.8,0.95)+
  labs(x = "Number of Trees", y = "Accuracy (Correlation)", title = "ntree Sensitivity Analysis") +
  theme_bw()
print(p2)

# ============================================================
# 4. Final Model Fitting and Variable Importance
# ============================================================

# Train the final model with optimal mtry
set.seed(13452)
RRF_model <- randomForest(
  form_reg,
  data = traindata,
  ntree = 125,
  mtry = best_mtry,
  importance = TRUE
)

# Plot variable importance
varImpPlot(RRF_model, type = 1, main = "Variable Importance - IncMSE",color="BLUE",pch=9)
varImpPlot(RRF_model, type = 2, main = "Variable Importance - NodePurity",color="BLUE",pch=9)

# ============================================================
# 5. Partial Dependence Analysis
# ============================================================

variables <- c("Chla", "SSS", "POC", "SST", "Pi", "N", "DOP")
partialPlot(x = RRF_model, 
            pred.data = traindata, 
            x.var = "SST", 
            col = "red", 
            lwd = 4, 
            xlab = "DIP (nmol/L)", 
            ylab = "scale phnJ (%)")

# ============================================================
# 6. Model Evaluation on Training and Testing Sets
# ============================================================

train_pred <- predict(RRF_model, newdata = traindata)
test_pred  <- predict(RRF_model, newdata = testdata)

train_obs_pred <- data.frame(obs = traindata$phnJ, pred = train_pred, Type = "Training")
test_obs_pred  <- data.frame(obs = testdata$phnJ, pred = test_pred, Type = "Testing")

# Performance metrics function
calc_metrics <- function(df) {
  list(
    n = nrow(df),
    R2 = round(R2(df$pred, df$obs)[1], 2),
    RMSE = round(rmse(df$obs, df$pred), 2),
    MAE = round(mae(df$obs, df$pred), 2)
  )
}

metrics_train <- calc_metrics(train_obs_pred)
metrics_test  <- calc_metrics(test_obs_pred)

# Print performance summary to console
cat("Model Performance on Training dataSet \n")
cat("Samples:", metrics_train$n, "\n")
cat("R2     :", metrics_train$R2, "\n")
cat("RMSE   :", metrics_train$RMSE, "\n")
cat("MAE    :", metrics_train$MAE, "\n\n")

cat("Model Performance on Testing dataSet \n")
cat("Samples:", metrics_test$n, "\n")
cat("R2     :", metrics_test$R2, "\n")
cat("RMSE   :", metrics_test$RMSE, "\n")
cat("MAE    :", metrics_test$MAE, "\n")

# ============================================================
# 7. Monte Carlo Prediction for External Dataset
# ============================================================

# Load external dataset for prediction
Yuce <- read.csv("D:/Desktop/data fet for prediction 1бу.csv")
Yuce <- na.omit(Yuce)

Yuce_normalized <- data.frame(Yuce)

# Set number of iterations
n_iterations <- 10000

result_matrix <- matrix(NA, nrow = nrow(Yuce_normalized), ncol = n_iterations)

set.seed(114053)

all_seeds <- sample(10000:99999, size = n_iterations, replace = FALSE)

for (i in 1:n_iterations) {
  
  set.seed(all_seeds[i])
  
  trains_idx <- createDataPartition(y = MMOX_unnormalize$phnJ, p = 0.8, list = FALSE)
  traindata_r <- MMOX_unnormalize[trains_idx, ]
  
  traindata_r$phnJ <- scale(traindata_r$phnJ, center = TRUE, scale = TRUE)
  
  # Train random forest model
  RRF_model <- randomForest(
    form_reg,
    data = traindata_r,
    ntree = 125,
    mtry = 2,
    importance = TRUE
  )
  
  # Predict on normalized external data
  RRF_result <- predict(RRF_model, newdata = Yuce_normalized)
  
  mean_phnJ <- mean(MMOX_unnormalize$phnJ)
  sd_phnJ   <- sd(MMOX_unnormalize$phnJ)
  RRF_predicted_values_unnormalized <- RRF_result * sd_phnJ + mean_phnJ
  
  result_matrix[, i] <- RRF_predicted_values_unnormalized
}

# Combine longitude, latitude, and prediction results
mean_pred <- rowMeans(result_matrix)
sd_pred <- apply(result_matrix, 1, sd)

final_result <- data.frame(
  Lon = Yuce$Lon,
  Lat = Yuce$Lat,
  phnJ_mean = mean_pred,
  phnJ_sd = sd_pred
)

write.csv(final_result, file = "D:/Desktop/phnJ_prediction_result.csv", row.names = FALSE)

