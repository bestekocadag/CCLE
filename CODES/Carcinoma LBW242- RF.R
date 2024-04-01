input_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE//OUTPUT/carcinoma/LBW242"
setwd(input_dir)
#install.packages("randomForest")
library(randomForest)
#install.packages("ranger")
library(ranger)
#install.packages("lifecycle")
library(lifecycle)
#install.packages("dplyr")
library(dplyr)
#install.packages("magrittr")
library(magrittr)
#install.packages("lifecycle")
library(lifecycle)
#install.packages("tidyverse")
library(tidyverse)

#.........................................................

LBW242_Drug <- read.delim("LBW242 - Expression Data (IC50) - Normalized.tsv")
dim(LBW242_Drug)
head(LBW242_Drug)
colnames(LBW242_Drug)
LBW242_Drug_without_sample_name <- subset(LBW242_Drug, select = c("IC50","CCNL1","CD58","FAM220A","GBP5","LINC01425",   
                                                                  "LINC02804","LOC101928326","LOC105371934","LOC105378146","MFAP2","NPB","PERM1",       
                                                                  "PITPNM3","RFLNB","RNF128","RPLP2","SNHG7","SPATA3.AS1"))
head(LBW242_Drug_without_sample_name)
#..........................................................

ind <- sample(2, nrow(LBW242_Drug_without_sample_name), replace=TRUE, prob=c(0.7, 0.3))
train_LBW242_Drug <- LBW242_Drug_without_sample_name[ind==1,]
head(train_LBW242_Drug)
dim(train_LBW242_Drug)

test_LBW242_Drug <- LBW242_Drug_without_sample_name[ind==2,]
head(test_LBW242_Drug)
dim(test_LBW242_Drug)

#........................................................

n_features <- length(setdiff(names(train_LBW242_Drug), "IC50"))

#........................................................

# Burada sadece random forest yapıldı ,test datası üzerinde henüz bir oynama yapmadan.
install.packages("ranger")
library(ranger)

LBW242_rf1 <- ranger(
  IC50 ~ ., 
  data = train_LBW242_Drug,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order",
  seed = 123
)

LBW242_rf1


# get OOB RMSE

(default_rmse <- sqrt(LBW242_rf1$prediction.error))


hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)

print(hyper_grid)

for(i in seq_len(nrow(hyper_grid))) {
  
  fit <- ranger(
    formula         = IC50 ~ ., 
    data            = train_LBW242_Drug, 
    num.trees       = n_features * 10,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  
  print(fit)
  
  hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

print(hyper_grid$rmse[i] )

install.packages("lifecycle")
library(lifecycle)
install.packages("tidyverse")
library(tidyverse)

hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) 
  head(10) 

hyper_grid_new <- hyper_grid[order(hyper_grid$rmse, decreasing = FALSE),] 
hyper_grid_new
head(hyper_grid_new, 10)
best_model_from_rf <- hyper_grid_new[1,]
best_model_from_rf
class(best_model_from_rf)

#...............

fit_best_model_rf <- ranger(
  formula         = IC50 ~ .,
  data            = train_LBW242_Drug,
  num.trees       = n_features * 10,
  mtry            = 6,
  min.node.size   = 3,
  replace         = FALSE,
  sample.fraction = 0.63,
  verbose         = FALSE,
  seed            = 123,
  respect.unordered.factors = 'order',
)
fit_best_model_rf

#....................................................

train_prediction <- predict(fit_best_model_rf, train_LBW242_Drug)
train_prediction

R_Square_Train <- cor(train_LBW242_Drug$IC50, train_prediction$prediction)^2
R_Square_Train

Pearson_Correlation_Train <- cor(train_LBW242_Drug$IC50,train_prediction$prediction)
Pearson_Correlation_Train

Pearson_Correlation_P_Value_Train <- cor.test(train_LBW242_Drug$IC50,train_prediction$prediction)$p.value
Pearson_Correlation_P_Value_Train

Spearman_Correlation_Train <- cor(train_LBW242_Drug$IC50, train_prediction$prediction, method = "spearman")
Spearman_Correlation_Train

Spearman_Correlation_P_Value_Train <- cor.test(train_LBW242_Drug$IC50, train_prediction$prediction, method = "spearman")$p.value
Spearman_Correlation_P_Value_Train

Square_Train <- (train_prediction$prediction - train_LBW242_Drug$IC50)^2
Mean_Square_Train <- mean(Square_Train)
RMSE_Train <- sqrt(Mean_Square_Train)
RMSE_Train

List_result_table_IC50_Train <- list()
IC50_Train_and_Prediction <- cbind(train_LBW242_Drug$IC50, train_prediction$prediction)
colnames(IC50_Train_and_Prediction) <- c("Actual_IC50_Train_Data","Predicted_IC50_Test_Data")
IC50_Train_and_Prediction
List_result_table_IC50_Train[["LBW242"]] <- IC50_Train_and_Prediction




test_prediction <- predict(fit_best_model_rf, test_LBW242_Drug)
test_prediction

R_Square_Test <- cor(test_LBW242_Drug$IC50, test_prediction$prediction)^2
R_Square_Test

Pearson_Correlation_Test <- cor(test_LBW242_Drug$IC50,test_prediction$prediction)
Pearson_Correlation_Test

Pearson_Correlation_P_Value_Test <- cor.test(test_LBW242_Drug$IC50,test_prediction$prediction)$p.value
Pearson_Correlation_P_Value_Test

Spearman_Correlation_Test <- cor(test_LBW242_Drug$IC50, test_prediction$prediction, method = "spearman")
Spearman_Correlation_Test

Spearman_Correlation_P_Value_Test <- cor.test(test_LBW242_Drug$IC50, test_prediction$prediction, method = "spearman")$p.value
Spearman_Correlation_P_Value_Test

Square_Test <- (test_prediction$prediction - test_LBW242_Drug$IC50)^2
Mean_Square_Test <- mean(Square_Test)
RMSE_Test <- sqrt(Mean_Square_Test)
RMSE_Test


List_result_table_IC50_Test <- list()
IC50_Test_and_Prediction <- cbind(test_LBW242_Drug$IC50, test_prediction$prediction)
colnames(IC50_Test_and_Prediction) <- c("Actual_IC50_Test_Data","Predicted_IC50_Test_Data")
IC50_Test_and_Prediction
List_result_table_IC50_Test[["LBW242"]] <- IC50_Test_and_Prediction





Result_table_drug <- data.frame(
  Tissue = c("Carcinoma"),
  Drug = c("LBW242"),
  R_Square_Train,
  R_Square_Test,
  
  Pearson_Correlation_Train,
  Pearson_Correlation_Test,
  
  Pearson_Correlation_P_Value_Train,
  Pearson_Correlation_P_Value_Test,
  
  Spearman_Correlation_Train,
  Spearman_Correlation_Test,
  
  Spearman_Correlation_P_Value_Train,
  Spearman_Correlation_P_Value_Test,
  
  RMSE_Train,
  RMSE_Test
  
)
