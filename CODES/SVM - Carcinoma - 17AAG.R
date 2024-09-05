install.packages("e1071")
library(e1071)

install.packages("caret")
library(caret)

input_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/INPUT"
output_dir <- "C:/Users/MONSTER/Desktop/TEZ/CCLE/OUTPUT"

setwd(paste0(output_dir,"/","carcinoma","/","17-AAG"))

data <- read.delim("17-AAG - Expression Data (IC50) - Normalized.tsv")
head(data)
dim(data)
data_without_sample_name <- subset(data, select = - Sample_Name)
head(data_without_sample_name)
dim(data_without_sample_name)

set.seed(1234)
trainIndex <- createDataPartition(data_without_sample_name$IC50, p = 0.7, list = FALSE)
train_Data <- data_without_sample_name[trainIndex, ]
test_Data <- data_without_sample_name[-trainIndex, ]


dim(train_Data)
dim(test_Data)

tuneResult1 <- tune(svm, IC50 ~. ,  data = train_Data ,
                    ranges = list(epsilon = seq(0,1,0.1), cost = 2^(seq(0.5,8,.5))),
                    type = "eps-regression",
                    kernel = "radial"
)


tuneResult <- tune(svm, IC50 ~. ,  data = train_Data,
                   ranges = list(epsilon = seq(tuneResult1$best.model$epsilon-.15,
                                               tuneResult1$best.model$epsilon+.15,
                                               0.01), 
                                 cost = seq(2^(log2(tuneResult1$best.model$cost)-1),
                                            2^(log2(tuneResult1$best.model$cost)+1),
                                            length=6)),
                   type = "eps-regression",
                   kernel = "radial"
)

tunedVals <-tuneResult$best.model
predictYsvm2 <- predict(tunedVals, train_Data)

predictYsvm2
