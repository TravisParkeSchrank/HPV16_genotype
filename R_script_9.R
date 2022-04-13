library(here)
library(stringr)
files <- read.table("hla_manifest.txt")
files$V1 <- as.character(files$V1)

hla_data <- read.table(here("../25B_HLA_calls_OptiType_FR/",files$V1[1]), header = TRUE, sep = "\t")

for(i in 2:length(files$V1)){
  holder <- read.table(here("../25B_HLA_calls_OptiType_FR/", files$V1[i]), header = TRUE, sep = "\t")
  hla_data <- rbind(hla_data, holder)
}

hla_data$file_origin <- files$V1

library(stringr)
sample_ID <- c()
sample_TN <- c()
#sample_FR <- c()
for(i in 1:length(hla_data$A1)){
  sample_ID <- c(sample_ID, str_split(hla_data$file_origin[i], pattern = "[.]")[[1]][1])
   #sample_TN <- c(sample_TN, str_split(hla_data$file_origin[i], pattern = "[.]")[[1]][2])
   #sample_FR <- c(sample_FR, str_split(hla_data$file_origin[i], pattern = "[.]")[[1]][3])
}

hla_data$sample_ID <- sample_ID
hla_data$sample_TN <- str_sub(hla_data$sample_ID, start = 12L, end =12L)
hla_data$study_id <- str_sub(hla_data$sample_ID, start = 1L, end =10L)
#hla_data$sample_FR <- sample_FR

hla_data_tumor <- subset(hla_data, hla_data$sample_TN == "T")
hpv_pos <- read.csv("104_hpv16_pos_cases.csv")


A <- c(hla_data_tumor$A1, hla_data_tumor$A2)
B <- c(hla_data_tumor$B1, hla_data_tumor$B2)
C <- c(hla_data_tumor$C1, hla_data_tumor$C2)
study_id <- c(hla_data_tumor$study_id, hla_data_tumor$study_id)

hla_data_tumor_simple <- data.frame(A =A, B = B, C= C, study_id = study_id)
hla_data_tumor_simple$hpv_status <- hla_data_tumor_simple$study_id %in% hpv_pos$study_id

A_test_data <- data.frame(table(hla_data_tumor_simple$hpv_status, hla_data_tumor_simple$A))
A_test_data_pos <- subset(A_test_data, A_test_data$Var1 == TRUE)
A_test_data_neg <- subset(A_test_data, A_test_data$Var1 == FALSE)
chisq.test(A_test_data_pos$Freq, A_test_data_neg$Freq)

B_test_data <- data.frame(table(hla_data_tumor_simple$hpv_status, hla_data_tumor_simple$B))
B_test_data_pos <- subset(B_test_data, B_test_data$Var1 == TRUE)
B_test_data_neg <- subset(B_test_data, B_test_data$Var1 == FALSE)
chisq.test(B_test_data_pos$Freq, B_test_data_neg$Freq)

C_test_data <- data.frame(table(hla_data_tumor_simple$hpv_status, hla_data_tumor_simple$C))
C_test_data_pos <- subset(C_test_data, C_test_data$Var1 == TRUE)
C_test_data_neg <- subset(C_test_data, C_test_data$Var1 == FALSE)
chisq.test(C_test_data_pos$Freq, C_test_data_neg$Freq)

library(ggplot2)

s <- ggplot(hla_data_tumor_simple, aes(hpv_status, fill = B))
s+geom_bar(position = "fill")


out_hla_data <- hla_data[, names(hla_data) %in% c("A1", "A2", "B1", "B2", "C1", "C2", "study_id", "sample_TN")]
write.table(out_hla_data, file = "out_hla_data.71*6csv", row.names = FALSE, sep = ",")

normal_hla_data <- subset(hla_data, hla_data$sample_TN == "N")
tumor_hla_data <- subset(hla_data, hla_data$sample_TN == "T")
tumor_hla_data_with_normal <- subset(tumor_hla_data, tumor_hla_data$study_id %in% normal_hla_data$study_id)

dim(normal_hla_data)
dim(tumor_hla_data_with_normal)

identical_calls <- c()
for(i in 1:71){
  identical_calls <- c(identical_calls, 
  identical(sort(c(normal_hla_data$A1[i],normal_hla_data$A2[i],normal_hla_data$B1[i],normal_hla_data$B2[i],normal_hla_data$C1[i],normal_hla_data$C2[i])), 
            sort(c(tumor_hla_data_with_normal$A1[i],tumor_hla_data_with_normal$A2[i],tumor_hla_data_with_normal$B1[i],tumor_hla_data_with_normal$B2[i],tumor_hla_data_with_normal$C1[i],tumor_hla_data_with_normal$C2[i])))
  )
}

call_identical_data <- data.frame(identical_TN = identical_calls, study_id = normal_hla_data$study_id)
not_identical_call_id_data <- subset(call_identical_data, call_identical_data$identical_TN == FALSE)

discordant_hla_call_examples <- subset(out_hla_data, out_hla_data$study_id %in% not_identical_call_id_data$study_id)
write.table(discordant_hla_call_examples, file = "discordant_calls_out_hla_data.csv", row.names = FALSE, sep = ",")


