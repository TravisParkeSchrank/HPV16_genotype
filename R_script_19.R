#library(Gviz)
#library(GenomicRanges)
library(Biostrings)
library(stringr)
#library(GenomicRanges)
library(janitor)
#library(rtracklayer)
#library(trackViewer)
library(stringr)
library(ggsignif)
library(ggplot2)
library("wesanderson")
library("ggpubr")


library(survival);
library(stringr)
library(survminer)
library(wesanderson)
surv_pal <- c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2])
surv_pal_sublin <- c(wes_palette(n=4, name="Royal1")[1],wes_palette(n=5, name="IsleofDogs2")[3])
## Needed For KAP MEIRE analysis stuff
#library(plotrix); 
#Needed for Multi Hist


data <- read.csv("clinical_data_update_7_21_21_ed_last_contact.csv", header = TRUE)
load("n_clonal_missense.RData")
merged <- merge(data, n_clonal_missense_vars, by = "study_id")

merged$stage_advanced_8ed <- as.numeric(as.character(merged$Sum_AJCC8) %in% c("IV", "III"))
merged$gt_10_pack_year <- as.numeric(as.numeric(merged$smoking_exposure) > 10)


merged <- subset(merged, merged$OP_subsite %in% c("BOT", "tonsil"))
###REMOVE METASTATIC CASES
merged <- subset(merged, !(merged$Sum_AJCC8 %in% c("IV")))
merged$smoking_exposure_zero <- merged$smoking_exposure
merged$smoking_exposure_zero[is.na(merged$smoking_exposure_zero)] <-0 
#merged <- subset(merged, smoking_exposure_zero < 50)

rfs_time_days <- c()


for(i in 1:length(merged$study_id)){
  if(merged$recurrence_y_n[i] == 1){
    #print("rec")
    rfs_time_days <- c(rfs_time_days, as.Date(merged$recurr_date[i], "%m/%d/%y") - as.Date(merged$date_of_diagnosis[i], "%m/%d/%y"))
    
  } else {
    #print("no rec")
    rfs_time_days <- c(rfs_time_days, as.Date(merged$date_last_contact[i], "%m/%d/%y") - as.Date(merged$date_of_diagnosis[i], "%m/%d/%y"))
  }
}
fu_time_days <- c()
for(i in 1:length(merged$study_id)){
  fu_time_days <- c(fu_time_days, as.Date(merged$date_last_contact[i], "%m/%d/%y") - as.Date(merged$date_of_diagnosis[i], "%m/%d/%y"))	
}

merged$rfs_time_days <- rfs_time_days
merged$fu_time_days <- fu_time_days


#merged$race_binary <- merged$Race %in% c(1)

load("far_clade_parsimony.RData")

merged$high_nonsyn_poly <- as.numeric(merged$study_id %in% far_clade_parsimony)
## OLD SECTION
merged$factor_n_poly <- as.factor(!(as.logical(merged$high_nonsyn_poly)))

#HA_peptide_out_data$study_id <- HA_peptide_out_data$sample_id
#merged <- merge(merged, HA_peptide_out_data, by = "study_id")


#### IMPORT CN DATA
cn_data <- read.table(file = "cn_stats.txt", header = TRUE, sep = "\t")

merged <- merge(merged, cn_data, by = "study_id", all=T)




merged$pST <- merged$Sum_AJCC8

merged$plt_stg <- merged$Sum_AJCC8
merged$plt_smo <- as.numeric(merged$smoking_exposure_zero)

merged$plt_cnv <- as.numeric(merged$n_events)

dictionary_n_poly <- c("BFar", "ANear")
names(dictionary_n_poly) <- c("TRUE", "FALSE")
merged$np_plt <- dictionary_n_poly[merged$factor_n_poly]

merged <- subset(merged, merged$np_plt %in% c("BFar", "ANear"))

#### bring in netMHCpan Data
load("n_HA_pept_data.RData")
merged <- merge(merged, n_HA_pept_data, by="study_id", all=T)
merged <- subset(merged, merged$np_plt %in% c("BFar", "ANear"))
merged$plt_ant <- as.numeric(merged$n_HA_pept)


high_n_poly <- subset(merged, merged$high_nonsyn_poly  == TRUE)
low_n_poly <- subset(merged, merged$high_nonsyn_poly == FALSE)

test <-ggplot(merged, aes(x=np_plt, fill=pST))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise3")[1],wes_palette("Moonrise3")[3:5]))+
  labs(y="AJCC8 Stage (Proportion)", x = "Clade")

pdf(file = "stage_vs_np_coding.pdf", width = 3.5, height = 4)
mybar
dev.off()


pik3ca_data <- data.frame("viral_clade" = c(rep("anear", times = 18), rep("far", times = 19)), "pik3ca_status" = c(rep("zvar", times = 7), rep("wt", times = 11), rep("zvar", times = 8), rep("wt", times = 11)))


test_pik3ca <-ggplot(pik3ca_data, aes(x=viral_clade, fill=pik3ca_status))
mybar_pik <- test_pik3ca + geom_bar(position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = c("grey75", "grey25"))+
  labs(y="PIK3CA Mutation (Proportion)", x = "Viral Clade")
  
  
pdf(file = "clade_vs_pik3ca.pdf", width = 3.5, height = 4)
mybar_pik
dev.off()

                                                                            
                                                                                 


bar_table <- t(table(merged$factor_n_poly, merged$pST))
print(chisq.test(bar_table))

my_box_pal <- c(wes_palette("Moonrise3")[5], wes_palette("Moonrise3")[1])

smoking_panel <-ggplot(merged, aes(x=np_plt, fill=np_plt, y = plt_smo)) + 
  geom_boxplot(aes(y=plt_smo)) + 
  theme_classic() + 
 #theme(legend.position = "none") +
  scale_fill_manual(values = my_box_pal) + 
  coord_cartesian(ylim = c(  min(merged$plt_smo[!is.na(merged$plt_smo)])*1 , max(merged$plt_smo[!is.na(merged$plt_smo)])*1.1)) +
  geom_signif(y_position = c(max(merged$plt_smo[!is.na(merged$plt_smo)])*1.05), xmin = c(1), xmax = c(2), annotation = c("NS"), tip_length = 0.03)+
  labs(y="Smoking Exposure (pack-years)", x = "Clade")
  

ant_panel <-ggplot(merged, aes(x=np_plt, fill=np_plt, y = plt_ant)) + 
  geom_boxplot(aes(y=plt_ant)) + 
  theme_classic() + 
  #theme(legend.position = "none") +
  scale_fill_manual(values = my_box_pal) + 
  coord_cartesian(ylim = c(  min(merged$plt_ant[!is.na(merged$plt_ant)])*1 , max(merged$plt_ant[!is.na(merged$plt_ant)])*0.9)) +
  geom_signif(y_position = c(max(merged$plt_ant[!is.na(merged$plt_ant)])*0.85), xmin = c(1), xmax = c(2), annotation = c("NS"), tip_length = 0.03)+
  labs(y="Antigenic HPV16 Peptides", x = "Clade")



cnv_panel <-ggplot(merged, aes(x=np_plt, fill=np_plt, y = plt_cnv)) + 
  geom_boxplot(aes(y=plt_cnv)) + 
  theme_classic() + 
  #theme(legend.position = "none") +
  scale_fill_manual(values = my_box_pal) + 
  coord_cartesian(ylim = c(  min(merged$plt_cnv[!is.na(merged$plt_cnv)])*1 , max(merged$plt_cnv[!is.na(merged$plt_cnv)])*1.1)) +
  geom_signif(y_position = c(max(merged$plt_cnv[!is.na(merged$plt_cnv)])*1.05), xmin = c(1), xmax = c(2), annotation = c("*"), tip_length = 0.03)+
  labs(y="Frequency Genomic CNV ", x = "Clade")


sublineage_figure <- ggarrange(mybar, smoking_panel,ant_panel,  cnv_panel,
                    labels = c("", "", "", "", ""),
                    ncol = 4, nrow = 1)

pdf(file = "phylo_box_figure.pdf", height = 2, width = 10)
sublineage_figure
dev.off()

print(wilcox.test(high_n_poly$n_cn_events, low_n_poly$n_cn_events))
print(wilcox.test(high_n_poly$n_events, low_n_poly$n_events))
print(wilcox.test(high_n_poly$smoking_exposure_zero, low_n_poly$smoking_exposure_zero))
print(wilcox.test(high_n_poly$n_HA_pept, low_n_poly$n_HA_pept))