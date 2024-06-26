#install.packages("e1071")
library(e1071)

#install.packages("caret")
library(caret)

input_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/INPUT"
output_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/OUTPUT"

setwd(output_dir)
Histologies <- list.dirs(full.name = FALSE, recursive = FALSE)

Drug_SVR_Result <- list()
Histology_SVR_Result <- list()

for(histology in Histologies){
  setwd(paste0(output_dir,"/",histology))
  if(length(list.dirs(full.name = FALSE, recursive = FALSE)) > 5){
    Drugs <- list.dirs(full.name = FALSE, recursive = FALSE)
  for(drug in Drugs){
    setwd(paste0(output_dir,"/",histology,"/",drug))
    drug_expression_data_normalized <- read.delim(paste(drug, "- Expression Data (IC50) - Normalized.tsv", sep = " "))
    if(length(drug_expression_data_normalized$Sample_Name) >=20){
       drug_expression_data_normalized_without_sample_name <- subset(drug_expression_data_normalized, 
                                                                     select = - Sample_Name)
      trainIndex <- createDataPartition(drug_expression_data_normalized_without_sample_name$IC50, 
                    p = 0.7, list = FALSE )
      train_Data <- drug_expression_data_normalized_without_sample_name[trainIndex, ]
      test_Data <- drug_expression_data_normalized_without_sample_name[-trainIndex, ]
    
      if(length(setdiff(names(train_Data), "IC50")) >= 5){
        cat(paste0("\nSupport Vector Regression is running for ", histology,"-",drug))
        flush.console()
        
        tuneResult_all_parameters <- tune(svm, IC50 ~. ,  data = train_Data ,
                                          ranges = list(epsilon = seq(0,1,0.1), cost = 2^(seq(0.5,8,.5))),
                                          type = "eps-regression",
                                          kernel = "radial"
        )
        
        
        tuneResult <- tune(svm, IC50 ~. ,  data = train_Data,
                            ranges = list(epsilon = seq(tuneResult_all_parameters$best.model$epsilon-.15,
                                                        tuneResult_all_parameters$best.model$epsilon+.15,
                                                        0.01), 
                                          cost = seq(2^(log2(tuneResult_all_parameters$best.model$cost)-1),
                                                     2^(log2(tuneResult_all_parameters$best.model$cost)+1),
                                                     length=6)),
                           type = "eps-regression",
                           kernel = "radial"
        )
        
        tunedVals <-tuneResult$best.model
        
        predict_IC50_train_svm <- predict(tunedVals, train_Data)
        
        Predict_IC50_Train_SVM <- cbind(train_Data$IC50, predict_IC50_train_svm)
        colnames(Predict_IC50_Train_SVM) <- c("Train_Data_IC50", "Predict_Train_Data_IC50")
        
        R_Square_Train <- cor(train_Data$IC50, predict_IC50_train_svm)^2
        
        Pearson_Correlation_Train <- cor(train_Data$IC50, predict_IC50_train_svm)
        
        Pearson_Correlation_P_Value_Train <- cor.test(train_Data$IC50, predict_IC50_train_svm)$p.value
        
        Spearman_Correlation_Train <- cor(train_Data$IC50, predict_IC50_train_svm, method = "spearman")
        
        Spearman_Correlation_P_Value_Train <- cor.test(train_Data$IC50, predict_IC50_train_svm, method = "spearman")$p.value
        
        Square_Train <- (predict_IC50_train_svm - train_Data$IC50)^2
        Mean_Square_Train <- mean(Square_Train)
        RMSE_Train <- sqrt(Mean_Square_Train)
        
        
        
        
        predict_IC50_test_svm <- predict(tunedVals, test_Data)
        
        Predict_IC50_Test_SVM <- cbind(test_Data$IC50, predict_IC50_test_svm)
        colnames(Predict_IC50_Test_SVM) <- c("Test_Data_IC50", "Predict_Test_Data_IC50")
        
        
        R_Square_Test <- cor(test_Data$IC50, predict_IC50_test_svm)^2
   
        Pearson_Correlation_Test <- cor(test_Data$IC50, predict_IC50_test_svm)
        
        Pearson_Correlation_P_Value_Test <- cor.test(test_Data$IC50, predict_IC50_test_svm)$p.value
        
        Spearman_Correlation_Test <- cor(test_Data$IC50, predict_IC50_test_svm, method = "spearman")
        
        Spearman_Correlation_P_Value_Test <- cor.test(test_Data$IC50, predict_IC50_test_svm, method = "spearman")$p.value
        
        Square_Test <- (predict_IC50_test_svm - test_Data$IC50)^2
        Mean_Square_Test <- mean(Square_Test)
        RMSE_Test <- sqrt(Mean_Square_Test)
       
        
        SVR_Drug_RESULT_TABLE <- data.frame(
                            Tissue = histology,
                            Drug = drug,
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
                            RMSE_Test)
        
        Drug_SVR_Result[[drug]] <- SVR_Drug_RESULT_TABLE
        dir.create(paste0(output_dir,"/",histology,"/",drug,"/","Support Vector Regression"))
        write.table(Drug_SVR_Result,
                    paste0(output_dir,"/",histology,"/",drug,"/","Support Vector Regression","/","Drug SVR Result.tsv"),
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                    )
        write.table(Predict_IC50_Train_SVM,
                    paste0(output_dir,"/",histology,"/",drug,"/","Support Vector Regression","/","Predict IC50 Train Data.tsv"),
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                    )
        
        write.table(Predict_IC50_Test_SVM,
                    paste0(output_dir,"/",histology,"/",drug,"/","Support Vector Regression","/","Predict IC50 Test Data.tsv"),
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
        )
        
        
      }
      
    }
  }
  
  Histology_SVR_Result[[histology]] <- ldply(Drug_SVR_Result, "rbind")
  write.table(Histology_SVR_Result,
              paste0(output_dir,"/",histology,"/",histology," - SVR Result.tsv"),
              row.names = FALSE,
              quote = FALSE,
              sep = "\t"
              
    
    )
  }
}

All_SVR_RESULT <- ldply(Histology_SVR_Result, "rbind")
write.table(All_SVR_RESULT,
            "C:/Users/MONSTER/Desktop/TEZ/CCLE/ML RESULT/SVR Result.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
