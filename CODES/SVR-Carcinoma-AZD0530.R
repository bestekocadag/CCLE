install.packages("e1071")
library(e1071)

install.packages("caret")
library(caret)

input_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/INPUT"
output_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/OUTPUT"

setwd(paste0(output_dir,"/","carcinoma","/","AZD0530"))

data2 <- read.delim("AZD0530 - Expression Data (IC50) - Normalized.tsv")
head(data2)
dim(data2)
data2_without_sample_name <- subset(data2, select = - Sample_Name)
head(data2_without_sample_name)
dim(data2_without_sample_name)

# Splitting The Data

set.seed(1234)
trainIndex2 <- createDataPartition(data2_without_sample_name$IC50, p = 0.7, list = FALSE)
train_Data2 <- data2_without_sample_name[trainIndex2, ]
test_Data2 <- data2_without_sample_name[-trainIndex2, ]


dim(train_Data2)
dim(test_Data2)


# Tuning Parameters

tuneResult_all_parameters <- tune(svm, IC50 ~. ,  data = train_Data2 ,
                                  ranges = list(epsilon = seq(0,1,0.1), cost = 2^(seq(0.5,8,.5)))
)


tuneResult2 <- tune(svm, IC50 ~. ,  data = train_Data2,
                    ranges = list(epsilon = seq(tuneResult_all_parameters$best.model$epsilon-.15,
                                                tuneResult_all_parameters$best.model$epsilon+.15,
                                                0.01), 
                                  cost = seq(2^(log2(tuneResult_all_parameters$best.model$cost)-1),
                                             2^(log2(tuneResult_all_parameters$best.model$cost)+1),
                                             length=6))
)

tunedVals2 <-tuneResult2$best.model


# Train Prediction

predict_IC50_train_svm2 <- predict(tunedVals2, train_Data2)

predict_IC50_train_svm2

Predict_IC50_Train_SVM <- cbind(train_Data2$IC50, predict_IC50_train_svm2)
colnames(Predict_IC50_Train_SVM) <- c("Train_Data_IC50", "Predict_Train_Data_IC50")


R_Square_Train2 <- cor(train_Data2$IC50, predict_IC50_train_svm2)^2
R_Square_Train2

Pearson_Correlation_Train2 <- cor(train_Data2$IC50, predict_IC50_train_svm2)
Pearson_Correlation_Train2

Pearson_Correlation_P_Value_Train2 <- cor.test(train_Data2$IC50, predict_IC50_train_svm2)$p.value
Pearson_Correlation_P_Value_Train2

Spearman_Correlation_Train2 <- cor(train_Data2$IC50, predict_IC50_train_svm2, method = "spearman")
Spearman_Correlation_Train2

Spearman_Correlation_P_Value_Train2 <- cor.test(train_Data2$IC50, predict_IC50_train_svm2, method = "spearman")$p.value
Spearman_Correlation_P_Value_Train2

Square_Train2 <- (predict_IC50_train_svm2 - train_Data2$IC50)^2
Mean_Square_Train2 <- mean(Square_Train2)
RMSE_Train2 <- sqrt(Mean_Square_Train2)
RMSE_Train2


write.table(Predict_IC50_Train_SVM,
            "C:/Users/MONSTER/Desktop/TEZ/CCLE/DENEME/SVM-Carcinoma-AZD0530-Train-IC50-Predict.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")




# Test Prediction

predict_IC50_test_svm2 <- predict(tunedVals2, test_Data2)

predict_IC50_test_svm2

Predict_IC50_Test_SVM <- cbind(test_Data2$IC50, predict_IC50_test_svm2)
colnames(Predict_IC50_Test_SVM) <- c("Test_Data_IC50", "Predict_Test_Data_IC50")


R_Square_Test2 <- cor(test_Data2$IC50, predict_IC50_test_svm2)^2
R_Square_Test2

Pearson_Correlation_Test2 <- cor(test_Data2$IC50, predict_IC50_test_svm2)
Pearson_Correlation_Test2

Pearson_Correlation_P_Value_Test2 <- cor.test(test_Data2$IC50, predict_IC50_test_svm2)$p.value
Pearson_Correlation_P_Value_Test2

Spearman_Correlation_Test2 <- cor(test_Data2$IC50, predict_IC50_test_svm2, method = "spearman")
Spearman_Correlation_Test2

Spearman_Correlation_P_Value_Test2 <- cor.test(test_Data2$IC50, predict_IC50_test_svm2, method = "spearman")$p.value
Spearman_Correlation_P_Value_Test2

Square_Test2 <- (predict_IC50_test_svm2 - test_Data2$IC50)^2
Mean_Square_Test2 <- mean(Square_Test2)
RMSE_Test2 <- sqrt(Mean_Square_Test2)
RMSE_Test2


SVM_RESULT_TABLE <- data.frame(
  Tissue = "Carcinoma",
  Drug = "AZD0530",
  R_Square_Train2,
  R_Square_Test2,
  Pearson_Correlation_Train2,
  Pearson_Correlation_Test2,
  Pearson_Correlation_P_Value_Train2,
  Pearson_Correlation_P_Value_Test2,
  Spearman_Correlation_Train2,
  Spearman_Correlation_Test2,
  Spearman_Correlation_P_Value_Train2,
  Spearman_Correlation_P_Value_Test2,
  RMSE_Train2,
  RMSE_Test2)



write.table(SVM_RESULT_TABLE,
            "C:/Users/MONSTER/Desktop/TEZ/CCLE/DENEME/SVM-Carcinoma-AZD0530-Result.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")




write.table(Predict_IC50_Test_SVM,
            "C:/Users/MONSTER/Desktop/TEZ/CCLE/DENEME/SVM-Carcinoma-AZD0530-Test-IC50-Predict.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

