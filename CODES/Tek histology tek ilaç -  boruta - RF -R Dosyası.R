suppressPackageStartupMessages(library(Boruta))

#install.packages("randomForest")
library(randomForest)

#install.packages("ranger")
library(ranger)

install.packages("dplyr")
library(dplyr)

#install.packages("magrittr")
library(magrittr)

#install.packages("tidyverse")
library(tidyverse)

#install.packages("lifecycle")
library(lifecycle)

output_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/OUTPUT"
input_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/INPUT"
deneme_output_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/Deneme Tek İlaç Boruta RF"

setwd(paste0(output_dir,"/","carcinoma"))

Histology_data <- read.delim("carcinoma - Histology.tsv")

indice_drug <- which(Histology_data$Drug == "PD-0332991")
Drug_data <- Histology_data[indice_drug,]


#BORUTA

setwd(input_dir)

ccle_rnaseq_data <- read.delim("GSE36133.BrainArray.RMAlog2Average.SYMBOL.Expr.tsv")
ccle_rnaseq_data <- t(ccle_rnaseq_data)

ccle_rnaseq_annotation_data <- read.delim("CCLE_expressions.annotations.tsv")
colnames(ccle_rnaseq_annotation_data)

head(Drug_data)
merge_histology_annotation <- merge(ccle_rnaseq_annotation_data,
                                    Drug_data,
                                    by.x = "title",
                                    by.y = "Cell_Line")

colnames(merge_histology_annotation)
filtered_histology_annotation <- merge_histology_annotation[,c("GSM","title","Drug","IC50","Activity_Area","histology")]
colnames(filtered_histology_annotation) <- c("Sample_Name","Cell_Line","Drug","IC50","Activity_Area","Histology")
head(filtered_histology_annotation)

if(length(unique(filtered_histology_annotation$Sample_Name)) >= 20){
  setwd(deneme_output_dir)
  drug <- "PD-0332991"
  dir.create(paste0(deneme_output_dir,"/",drug))
  setwd(paste0(deneme_output_dir,"/",drug))
  selected_drug_sample_info <- filtered_histology_annotation[filtered_histology_annotation$Drug %in% drug , ]
  selected_drug_samples_rnaseq_data <- ccle_rnaseq_data[rownames(ccle_rnaseq_data) %in% selected_drug_sample_info$Sample_Name,]
  selected_drug_samples_rnaseq_data_final <- cbind(rownames(selected_drug_samples_rnaseq_data),
                                                   selected_drug_samples_rnaseq_data)
  colnames(selected_drug_samples_rnaseq_data_final)[1] <- "Sample_Name"
  rownames(selected_drug_samples_rnaseq_data_final) <- NULL
  boruta_ready_data <- merge(selected_drug_sample_info,
                             selected_drug_samples_rnaseq_data_final,
                             by = "Sample_Name")
  
  IC50_boruta_ready_data <- boruta_ready_data[,c(4, 7:ncol(boruta_ready_data))]
  rownames(IC50_boruta_ready_data) <- boruta_ready_data$Sample_Name
  IC50_boruta_ready_data <- as.data.frame(sapply(IC50_boruta_ready_data, as.numeric))
  boruta <- Boruta(IC50 ~ ., data = IC50_boruta_ready_data, doTrace = 2, maxRuns = 100)
  boruta_check <- TentativeRoughFix(boruta)
  boruta_stats <- attStats(boruta_check)
  head(boruta_stats[order(boruta_stats$meanImp, decreasing = TRUE),], 20) 
  boruta_confirmed_features <- boruta_stats[boruta_stats$decision == "Confirmed",]  
  selected_drug_features <- rownames(boruta_confirmed_features)
  selected_drug_after_boruta_data <- boruta_ready_data[,c("Sample_Name","IC50",selected_drug_features)] 
  setwd(paste0(deneme_output_dir,"/",drug))
  write.table(selected_drug_after_boruta_data,
              paste(drug, "- Expression Data (IC50) Deneme1.tsv" , sep = " "),
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
  
}


#Normalization............

setwd(paste0(deneme_output_dir,"/",drug))
Drug_data2 <- read.delim(paste(drug, "- Expression Data (IC50) Deneme1.tsv", sep = " "))
indice_IC50_8 <- which(Drug_data2$IC50 == "8")
Drug_data_without_IC50_8<- Drug_data2[-indice_IC50_8,]
Drug_Data2_New <- Drug_data_without_IC50_8

write.table(Drug_Data2_New,
            paste(drug, "- Expression Data (IC50) Deneme1 - Normalized.tsv"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

#Random Forest....................................

drug <- "PD-0332991"
setwd(paste0(deneme_output_dir,"/",drug))
Drug_data3 <- read.delim(paste(drug, "- Expression Data (IC50) Deneme1 - Normalized.tsv" , sep = " "))

if(length(Drug_data3$Sample_Name) >= 20){
  Drug_data_without_sample_name <- subset(Drug_data3, select = - Sample_Name)
  ind <- sample(2, nrow(Drug_data_without_sample_name), replace=TRUE, prob=c(0.7, 0.3))
  train_Drug_data <- Drug_data_without_sample_name[ind==1,]
  test_Drug_data <- Drug_data_without_sample_name[ind==2,]
  
  
  if(length(setdiff(names(train_Drug_data), "IC50")) >= 5){
    
    n_features <- length(setdiff(names(train_Drug_data), "IC50"))  
    
    Drug_data_rf1 <- ranger(
      IC50 ~ ., 
      data = train_Drug_data,
      mtry = floor(n_features / 3),
      respect.unordered.factors = "order",
      seed = 123
    )
    
    (default_rmse <- sqrt(Drug_data_rf1$prediction.error)) 
    
    hyper_grid <- expand.grid(mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
                              min.node.size = c(1, 3, 5, 10), 
                              replace = c(TRUE, FALSE),                               
                              sample.fraction = c(.5, .63, .8),                       
                              rmse = NA                                               
    )
    
    for(i in seq_len(nrow(hyper_grid))) {
      
      fit <- ranger(
        formula         = IC50 ~ ., 
        data            = train_Drug_data, 
        num.trees       = n_features * 10,
        mtry            = hyper_grid$mtry[i],
        min.node.size   = hyper_grid$min.node.size[i],
        replace         = hyper_grid$replace[i],
        sample.fraction = hyper_grid$sample.fraction[i],
        verbose         = FALSE,
        seed            = 123,
        respect.unordered.factors = 'order',
      )
      hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
    }
    
    
    hyper_grid %>%
      arrange(rmse) %>%
      mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) 
    
    
    
    hyper_grid_new <- hyper_grid[order(hyper_grid$rmse, decreasing = FALSE),] 
    hyper_grid_new
    head(hyper_grid_new, 10)
    best_model_from_rf <- hyper_grid_new[1,]
    best_model_from_rf
    class(best_model_from_rf)     
    
    fit_best_model_rf <- ranger(
      formula         = IC50 ~ .,
      data            = train_Drug_data,
      num.trees       = n_features * 10,
      mtry            = best_model_from_rf[1,1],
      min.node.size   = best_model_from_rf[1,2],
      replace         = best_model_from_rf[1,3],
      sample.fraction = best_model_from_rf[1,4],
      verbose         = FALSE,
      seed            = 123,
      respect.unordered.factors = 'order',
    )     
    
    # For Train Data
    
    train_prediction <- predict(fit_best_model_rf, train_Drug_data)
    train_prediction
    
    R_Square_Train <- cor(train_Drug_data$IC50, train_prediction$prediction)^2
    R_Square_Train
    
    Pearson_Correlation_Train <- cor(train_Drug_data$IC50,train_prediction$prediction)
    Pearson_Correlation_Train
    
    Pearson_Correlation_P_Value_Train <- cor.test(train_Drug_data$IC50,train_prediction$prediction)$p.value
    Pearson_Correlation_P_Value_Train
    
    Spearman_Correlation_Train <- cor(train_Drug_data$IC50, train_prediction$prediction, method = "spearman")
    Spearman_Correlation_Train
    
    Spearman_Correlation_P_Value_Train <- cor.test(train_Drug_data$IC50, train_prediction$prediction, method = "spearman")$p.value
    Spearman_Correlation_P_Value_Train
    
    Square_Train <- (train_prediction$prediction - train_Drug_data$IC50)^2
    Mean_Square_Train <- mean(Square_Train)
    RMSE_Train <- sqrt(Mean_Square_Train)
    RMSE_Train
    
    IC50_Train_and_Prediction <- cbind(train_Drug_data$IC50, train_prediction$prediction)
    colnames(IC50_Train_and_Prediction) <- c("Actual_IC50_Train_Data","Predicted_IC50_Test_Data")
    IC50_Train_and_Prediction
    
    
    # For Test Data
    
    test_prediction <- predict(fit_best_model_rf, test_Drug_data)
    test_prediction
    
    R_Square_Test <- cor(test_Drug_data$IC50, test_prediction$prediction)^2
    R_Square_Test
    
    Pearson_Correlation_Test <- cor(test_Drug_data$IC50,test_prediction$prediction)
    Pearson_Correlation_Test
    
    Pearson_Correlation_P_Value_Test <- cor.test(test_Drug_data$IC50,test_prediction$prediction)$p.value
    Pearson_Correlation_P_Value_Test
    
    Spearman_Correlation_Test <- cor(test_Drug_data$IC50, test_prediction$prediction, method = "spearman")
    Spearman_Correlation_Test
    
    Spearman_Correlation_P_Value_Test <- cor.test(test_Drug_data$IC50, test_prediction$prediction, method = "spearman")$p.value
    Spearman_Correlation_P_Value_Test
    
    Square_Test <- (test_prediction$prediction - test_Drug_data$IC50)^2
    Mean_Square_Test <- mean(Square_Test)
    RMSE_Test <- sqrt(Mean_Square_Test)
    RMSE_Test
    
    
    IC50_Test_and_Prediction <- cbind(test_Drug_data$IC50, test_prediction$prediction)
    colnames(IC50_Test_and_Prediction) <- c("Actual_IC50_Test_Data","Predicted_IC50_Test_Data")
    IC50_Test_and_Prediction
    
    # Result Table
    
    Result_table_drug <- data.frame(
      Tissue = "Carcinoma",
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
      RMSE_Test
      
    )
    
    dir.create(paste0(deneme_output_dir,"/",drug,"/","Random Forest Deneme1"))
    setwd(paste0(deneme_output_dir,"/",drug,"/","Random Forest Deneme1"))
    
    # Writing 
    
    write.table(Result_table_drug,
                paste(drug, "- Result_RF_deneme1.tsv" , sep = " "),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    write.table(IC50_Train_and_Prediction,
                paste(drug, "- IC50_Train_Result_deneme1.tsv" , sep = " "),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    write.table(IC50_Test_and_Prediction,
                paste(drug, "- IC50_Test_Result_deneme1.tsv" , sep = " "),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    
    
    
    
  }
  
}
