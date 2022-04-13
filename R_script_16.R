
### see https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#4_overview_of_the_package

#require(TCGAbiolinksGUI)
require(DT)
require(dplyr)
require(SummarizedExperiment)
library(here)
library(stringr)
library(seqinr)

##### truncate a sample name to make it in normal TCGA format
truncate <- function(x){
  return(substr(x, 1, 12))
}
##### truncate a sample name to make it in normal TCGA format
truncate_new <- function(x){
  return(str_sub(x, start = -12L, end = -3L))
}
##### import hpv pos

hpv_pos_sample_data <- read.csv("104_hpv16_pos_cases.csv")

#group_coverage <- read.table(file = paste("UNCseq0733", ".hpv16.coverage.txt", sep = "" ))
##
##
##  Extract the coverage data
##

#blank_df_initial <- blank_df_initial[, !(names(blank_df_initial) %in% c("REF"))]
##########
##########
blank_df_initial <- data.frame("pos" = 1:7906, "UNCseqTEST" = rep(0, times = 7906))

for(i in 1:length(hpv_pos_sample_data$study_id)){
  blank_df_initial[as.character(hpv_pos_sample_data$study_id[i])] <- rep(0, times = 7906)
}


test <- blank_df_initial
##
##
##  Extract the coverage data
##
test_x <- test

samples <- hpv_pos_sample_data$study_id
for(i in samples){
  coverage <- read.table(file = paste("/home/parke/Desktop/HPV16_genotype_paper_revision_CLEAN/24_deep_loss_analysis/old_depth_files/", i, ".hpv16.coverage.txt", sep = "" ))
  for(ii in 1:length(coverage$V2)){
      test_x[[i]][coverage$V2[ii]] <- coverage$V3[ii]
    }
  }


grouped_coverage <- test_x[, names(test_x) %in% hpv_pos_sample_data$study_id]

grouped_coverage_remove_3150_3351 <- grouped_coverage[!(row.names(grouped_coverage) %in% as.character(3150:3351)), ]
grouped_coverage_3150_3351 <- grouped_coverage[(row.names(grouped_coverage) %in% as.character(3150:3351)), ]


grouped_deep_loss <- data.frame(grouped_coverage)


for(i in samples){
  grouped_deep_loss[[i]] <- grouped_deep_loss[[i]] < (mean(grouped_deep_loss[[i]][7000:7906])/100)
  grouped_deep_loss[[i]] <- as.numeric(grouped_deep_loss[[i]]) 
}

grouped_deep_loss <- grouped_deep_loss[ , order(colSums(grouped_deep_loss))]
pdf(file = "deep_loss_heatmap.pdf", width = 10, height = 10)
heatmap(t(as.matrix(grouped_deep_loss)), Colv = NA, Rowv = NA, scale = "none")
dev.off()



#pdf(file = "hist_mean_coverage.pdf", width = 10, height = 3)
pdf(file = "hist_mean_coverage.pdf", width = 10, height = 3)
plot(1:7906, rowMeans(grouped_coverage), cex = 0.25)
lines(spline(1:7906, rowMeans(grouped_coverage)), lwd = 3)
dev.off()


deep_loss_remove_3150_3351 <- grouped_deep_loss[!(row.names(grouped_deep_loss) %in% as.character(3150:3351)), ]

E2_deep_loss <- grouped_deep_loss[(row.names(grouped_deep_loss) %in% as.character(1892:2989)), ]


deep_loss_remove_3150_3351 <- grouped_deep_loss[!(row.names(grouped_deep_loss) %in% as.character(3150:3351)), ]

E2_deep_loss <- grouped_deep_loss[(row.names(grouped_deep_loss) %in% as.character(1892:2989)), ]

copy_loss_data_out <- data.frame(colSums(deep_loss_remove_3150_3351) > 7704/10)
names(copy_loss_data_out) <- c("quarter_loss")
copy_loss_data_out$study_id <-row.names(copy_loss_data_out)

copy_loss_data_E2 <- data.frame(colSums(E2_deep_loss) > 0)
names(copy_loss_data_E2) <- c("E2_loss")
copy_loss_data_E2$study_id <-row.names(copy_loss_data_E2)


copy_loss_data_out_op <- merge(copy_loss_data_out, copy_loss_data_E2, by = "study_id")

save(copy_loss_data_out_op, file = "copy_loss_data_out_op.RData")
