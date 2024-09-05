install.packages("writexl")
library(writexl)

setwd("C:/Users/MONSTER/Desktop/TEZ/CCLE/OUTPUT3")
ALL_HISTOLOGIES_RF <- read.delim("All_Drugs_RF_Result_3.tsv")
ALL_HISTOLOGIES_RF

filter <- order(ALL_HISTOLOGIES_RF$R_Square_Test, decreasing = TRUE)
ALL_HISTOLOGIES_RF <- ALL_HISTOLOGIES_RF[filter,]

ALL_HISTOLOGIES_RF

getwd()

write_xlsx(ALL_HISTOLOGIES_RF, 
           "ALL HISTOLOGIES RF 3.xlsx")


