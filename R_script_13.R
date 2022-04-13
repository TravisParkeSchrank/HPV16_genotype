load("out_tumor_hpv_integr_data_lisle.RData")
load("far_clade_parsimony.RData")


data <- read.csv("clinical_data_update_7_21_21_ed_last_contact.csv", header = TRUE)
merged <- data
merged <- subset(merged, merged$OP_subsite %in% c("BOT", "tonsil"))

lisle <- out_tumor_hpv_integr_data
library("stringr")
lisle$study_id <- str_sub(lisle$sample, start = 1L, end = 10L)

lisle <- subset(lisle, lisle$study_id %in% merged$study_id)

lisle$far_clade <- lisle$study_id %in% far_clade_parsimony



### my integration calls
load("int_call_data.RData")
load("vars_with_apobec_status.RData")
vcf_data_with_effect_all <- apobec_data
subclonal_vars <- vcf_data_with_effect_all
subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf < 90)
subclonal <- subset(subclonal_vars, subclonal_vars$vaf > 5)

n_sub <- table(subclonal$study_id)
n_sub <- data.frame(n_sub)
names(n_sub) <- c("tumors", "n_subclonal")
merge_int_subclonal <- merge(int_call_data, n_sub, by = "tumors", all.x = TRUE)
merge_int_subclonal$n_subclonal[is.na(merge_int_subclonal$n_subclonal)] <- 0
load("far_clade_parsimony.RData")
merge_int_subclonal$far_clade <- merge_int_subclonal$tumors %in% far_clade_parsimony

my_bar_pal <- c(wes_palette(n=5, name="IsleofDogs2")[1], wes_palette(n=5, name="IsleofDogs2")[3])
test_int <-ggplot(merge_int_subclonal, aes(x=far_clade, fill=highest_conf_int))
mybar_int <- test_int+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = my_bar_pal)+
  labs(y="Integrated", x = "Far Clade")
print(chisq.test(table(merge_int_subclonal$far_clade, merge_int_subclonal$highest_conf_int)))

pdf(file = "phylo_integration_bar.pdf", height = 2, width = 3)
mybar_int
dev.off()

##
##
## Subclonal 
##

data <- read.csv("clinical_data_update_7_21_21_ed_last_contact.csv", header = TRUE)
merged <- data
merged <- subset(merged, merged$OP_subsite %in% c("BOT", "tonsil"))
merged$far_clade <- merged$study_id %in% far_clade_parsimony

load("all_vars_lt95_apobec_annotated.RData")

vcf_data_with_effect_all <- apobec_data_lt95
subclonal_vars <- vcf_data_with_effect_all
subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf < 95)
subclonal <- subclonal_vars#subset(subclonal_vars, subclonal_vars$vaf > 5)
subclonal$depth <- as.numeric(subclonal$deth)
subclonal <- subclonal[subclonal$depth > 450, ]

n_sub <- table(subclonal$study_id)
n_sub <- data.frame(n_sub)
names(n_sub) <- c("study_id", "n_subclonal")

merged_subclonal <- merge(merged, n_sub, by = "study_id", all.x = TRUE)
merged_subclonal$n_subclonal[is.na(merged_subclonal$n_subclonal)] <- 0
merged_subclonal$log_subclonal <- log(merged_subclonal$n_subclonal + 1, base = 10)

load("out_sublineage_data.RData")

A1_data <- subset(out_sublineage_data, out_sublineage_data$sublineage_JC69 %in% c("A1"))

merged_subclonal$A1_status <- merged_subclonal$study_id %in% A1_data$study_id
merged_subclonal$has_subclone <- merged_subclonal$n_subclonal > 0

#merged_subclonal <- subset(merged_subclonal, merged_subclonal$smoking_exposure <= 50) 

plot(as.factor(merged_subclonal$far_clade), merged_subclonal$log_subclonal)
points(as.factor(merged_subclonal$far_clade), merged_subclonal$log_subclonal)
print(wilcox.test(n_subclonal ~ far_clade, data = merged_subclonal))


library(stringr)
library(ggsignif)
library(ggplot2)
library("wesanderson")
library("ggpubr")
library("ggbeeswarm")


library(survival);
library(stringr)
library(survminer)
library(wesanderson)
surv_pal <- c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2])
surv_pal_sublin <- c(wes_palette(n=4, name="Royal1")[1],wes_palette(n=5, name="IsleofDogs2")[3])

my_box_pal <- c(wes_palette("Moonrise3")[5], wes_palette("Moonrise3")[1])

subclonal_panel <-ggplot(merged_subclonal, aes(x=far_clade, y = n_subclonal)) + 
  geom_boxplot(aes(y=n_subclonal), outlier.alpha = 0.0) +
  geom_jitter(height = 0.2, width = 0.15, aes(color=far_clade), cex = 3, alpha = 0.5) +
  #geom_dotplot(binaxis = "y", stackdir = "center",  aes(y=log_subclonal, fill = far_clade), position=position_dodge(0.05), binwidth = 0.02255)+
  theme_classic() +
  scale_color_manual(values = surv_pal) + 
  #coord_cartesian(ylim = c(  min(merged$plt_smo[!is.na(merged$plt_smo)])*1 , max(merged$plt_smo[!is.na(merged$plt_smo)])*1.1)) +
  #geom_signif(y_position = c(max(merged$plt_smo[!is.na(merged$plt_smo)])*1.05), xmin = c(1), xmax = c(2), annotation = c("NS"), tip_length = 0.03)+
  labs(y="N. Subclonal Poly.", x = "Clade")
print(wilcox.test(n_subclonal ~ far_clade, merged_subclonal))

my_bar_pal <- c(wes_palette(n=5, name="IsleofDogs2")[1], wes_palette(n=5, name="IsleofDogs2")[3])
test <-ggplot(merged_subclonal, aes(x=far_clade, fill=has_subclone))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = my_bar_pal)+
  labs(y="Has Subclonal", x = "Far Clade")
print(chisq.test(table(merged_subclonal$far_clade, merged_subclonal$has_subclone)))

pdf(file = "phylo_subclone_bar.pdf", height = 2, width = 3)
mybar
dev.off()

subclonal <- data.frame(subclonal)
subclonal_gt5 <-subclonal[c(subclonal$vaf > 5), ]
freq_subclonal_gt5  <- data.frame(table(subclonal_gt5$study_id))
multip_subclonal_gt5  <- subset(freq_subclonal_gt5 , freq_subclonal_gt5$Freq > 3)

multip_subclonal_vars_gt5  <- subclonal_gt5[c(subclonal_gt5$study_id %in% multip_subclonal_gt5$Var1), ]
multip_subclonal_vars_gt5$far_clade <- multip_subclonal_vars_gt5$study_id %in% far_clade_parsimony

## PLOT THIS OUT BELOW TO ORGANIZE by N Clones

#subclonal_clone_groups <-ggplot(multip_subclonal_vars_gt5 , aes(x=study_id, y = vaf, fill = POS, shape = apobec_status)) + 
#  geom_beeswarm(size= 2, alpha = 0.85, aes(fill = POS), cex = 1.25) +
#  scale_fill_continuous(type = "viridis")+
#  scale_shape_manual(values = c(21,22,23))+
#  theme_classic()+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#
#pdf(file = "raw_subclone_vars_annotated.pdf", height = 3.5, width = 8)
#subclonal_clone_groups
#dev.off()


library(Canopy) ## I think this package will be adequate for
subclonal_gt1 <- subclonal[subclonal$vaf > 1.5, ]
freq_subclonal <- data.frame(table(subclonal_gt1$study_id))
multip_subclonal <- subset(freq_subclonal, freq_subclonal$Freq > 3)
multip_subclonal_vars <- subclonal[c(subclonal$study_id %in% multip_subclonal$Var1), ]
multip_subclonal_vars$far_clade <- multip_subclonal_vars$study_id %in% far_clade_parsimony
multip_subclonal_vars$depth_var <- as.integer(multip_subclonal_vars$depth  * (0.01 * multip_subclonal_vars$vaf))
multip_subclonal_vars$var_id <- paste(multip_subclonal_vars$REF, multip_subclonal_vars$POS, multip_subclonal_vars$ALT, sep = "")



n_hpv16_clones <- c()
for(i in multip_subclonal$Var1){

canopy_data <- subset(multip_subclonal_vars, multip_subclonal_vars$study_id %in% c(i))
canopy_data <- subset(canopy_data, canopy_data$vaf > 1.5)

my_X <- as.numeric(canopy_data$depth)
names(my_X) <- canopy_data$var_id
my_X <- data.frame(my_X)
names(my_X) <- "tumor1"
my_X <- as.matrix(my_X)

my_R <- as.numeric(canopy_data$depth_var)
names(my_R) <- canopy_data$var_id
my_R <- data.frame(my_R)
names(my_R) <- "tumor1"
my_R <- as.matrix(my_R)

num_cluster=c()
for(i in 2:7){
  if(i < dim(my_X)[1]){
    num_cluster <- c(num_cluster, i)
  }
}


num_run=10
my_canopy_cluster = canopy.cluster(R= my_R, X = my_X, num_cluster = num_cluster, num_run = num_run)

n_hpv16_clones <- c(n_hpv16_clones,  sum(as.numeric(my_canopy_cluster$Mu > 0.05)))

#bic_output=my_canopy_cluster$bic_output
#plot(num_cluster,bic_output,xlab='Number of mutation clsuters',ylab='BIC',type='b',main='BIC for model selection')
#abline(v=num_cluster[which.max(bic_output)],lty=2)

}

num_hpv16_subclones_data <- data.frame("study_id" = multip_subclonal$Var1, "n_hpv_clones" = n_hpv16_clones)


save(num_hpv16_subclones_data, file = "num_hpv16_subclones_data.RData")






## ORGANIZE by N Clones

library(tidyr)
library(tidyverse)
n_clone_dictionary <- num_hpv16_subclones_data$n_hpv_clones
names(n_clone_dictionary) <- as.character(num_hpv16_subclones_data$study_id)

n_clones <- c()
for(i in multip_subclonal_vars_gt5$study_id){
  n_clones <- c(n_clones, n_clone_dictionary[[i]])
}

multip_subclonal_vars_gt5$n_clones <- n_clones

multip_subclonal_vars_gt5 <- multip_subclonal_vars_gt5[rev(order(multip_subclonal_vars_gt5$n_clones)), ]

multip_subclonal_vars_gt5$study_id <- as.factor(multip_subclonal_vars_gt5$study_id)
multip_subclonal_vars_gt5$study_id <- factor(multip_subclonal_vars_gt5$study_id , levels=unique(as.character(multip_subclonal_vars_gt5$study_id)))



subclonal_clone_groups <-ggplot(multip_subclonal_vars_gt5 , aes(x=study_id, y = vaf, fill = POS, shape = apobec_status)) + 
  geom_beeswarm(size= 2, alpha = 0.85, aes(fill = POS), cex = 1.15) +
  scale_fill_continuous(type = "viridis")+
  scale_shape_manual(values = c(21,22,23))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "raw_subclone_vars_annotated.pdf", height = 3.5, width = 9)
subclonal_clone_groups
dev.off()


n_clones_plot_data <- multip_subclonal_vars_gt5[, names(multip_subclonal_vars_gt5) %in% c("study_id", "n_clones")]
n_clones_plot_data <- unique(n_clones_plot_data)

n_subclones_bar <- ggplot(n_clones_plot_data, aes(x=study_id, y = n_clones)) +
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "n_subclone_for_vars_annotated.pdf", height = 1.75, width = 9)
n_subclones_bar
dev.off()


####
####
### merge clinical n clones
###
names(num_hpv16_subclones_data) <- c("study_id", "n_hpv_clones")
merged_clin_clones <- merge(merged, num_hpv16_subclones_data, by = "study_id", all.x = TRUE)
merged_clin_clones$n_hpv_clones[is.na(merged_clin_clones$n_hpv_clones)] <- 0



wilcox.test(n_hpv_clones ~ far_clade, merged_clin_clones)

all_samples_subclones_bar <- ggplot(merged_clin_clones, aes(x = n_hpv_clones)) +
  geom_histogram()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "n_subclone_all_hit.pdf", height = 1.75, width = 1.75)
all_samples_subclones_bar
dev.off()


####
###
rfs_time_days <- c()

for(i in 1:length(merged$study_id)){
  if(merged$recurrence_y_n[i] == 1){
    #print("rec")
    rfs_time_days <- c(rfs_time_days, as.Date(merged_clin_clones$recurr_date[i], "%m/%d/%y") - as.Date(merged_clin_clones$date_of_diagnosis[i], "%m/%d/%y"))
    
  } else {
    #print("no rec")
    rfs_time_days <- c(rfs_time_days, as.Date(merged_clin_clones$date_last_contact[i], "%m/%d/%y") - as.Date(merged_clin_clones$date_of_diagnosis[i], "%m/%d/%y"))
  }
}
fu_time_days <- c()
for(i in 1:length(merged$study_id)){
  fu_time_days <- c(fu_time_days, as.Date(merged_clin_clones$date_last_contact[i], "%m/%d/%y") - as.Date(merged_clin_clones$date_of_diagnosis[i], "%m/%d/%y"))	
}

merged_clin_clones$rfs_time_days <- rfs_time_days
merged_clin_clones$fu_time_days <- fu_time_days

surv_data <- subset(merged_clin_clones, !(merged_clin_clones$Sum_AJCC8 %in% c("IV")))
surv_data$poly_clonal_HPV <- surv_data$n_hpv_clones > 0
#surv_data <- subset(surv_data, surv_data$smoking_exposure <= 20)

CoxFit <- coxph(Surv(rfs_time_days, recurrence_y_n) ~ poly_clonal_HPV, data = surv_data);
print(summary(CoxFit))

