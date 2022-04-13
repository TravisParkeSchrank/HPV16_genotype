#library(Gviz)
#library(GenomicRanges)
library(Biostrings)
library(stringr)
#library(GenomicRanges)
library(janitor)
#library(rtracklayer)
#library(trackViewer)

library(ggplot2)


load("common_poly_data_cervix.RData")
load("common_poly_data_oropharynx.RData")

common_poly_data_cervix$origin <- rep("UC", times = length(common_poly_data_cervix$var_id))


common_poly_data_op$origin <- rep("OP", times = length(common_poly_data_op$var_id))


common_op <-head(common_poly_data_op, n = 15)

common_uc <- subset(common_poly_data_cervix, common_poly_data_cervix$var_id %in% common_op$var_id)


merged <- rbind(common_op, common_uc)



myplt <- ggplot(data = merged, aes(x = var_id, y = prop, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  geom_errorbar(aes(x= var_id, y= prop, ymax = high_ci, ymin = low_ci, fill= origin), position = "dodge", width = 0.5)+
  theme_classic() +
  scale_fill_manual(values = c("grey25", "grey75"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
pdf(file = "compare_common_poly_OP_UC.pdf", height = 2.5, width = 6)
print(myplt)
dev.off()





load("copy_loss_data_out_op.RData")

n_tumors <- 104
prop <- c()
low_ci <- c()
high_ci <- c()
var_id <- c("quarter_loss", "E2_loss")
origin <- c("OP", "OP")

for(i in var_id){
  prop <- c(prop, 1 - prop.test(table(copy_loss_data_out_op[[i]]))[[4]][[1]])
  high_ci <- c(high_ci, 1 - prop.test(table(copy_loss_data_out_op[[i]]))[[6]][1])
  low_ci <- c(low_ci, 1 - prop.test(table(copy_loss_data_out_op[[i]]))[[6]][2])
}
plot_CL_data_op <- data.frame(prop = prop, low_ci = low_ci, high_ci = high_ci, var_id, origin  = origin )



load("copy_loss_data_out_cervix.RData")

n_tumors <- 44
prop <- c()
low_ci <- c()
high_ci <- c()
var_id <- c("quarter_loss", "E2_loss")
origin <- c("UC", "UC")

for(i in var_id){
  prop <- c(prop, 1 - prop.test(table(copy_loss_data_out_cervix[[i]]))[[4]][[1]])
  high_ci <- c(high_ci, 1 - prop.test(table(copy_loss_data_out_cervix[[i]]))[[6]][1])
  low_ci <- c(low_ci, 1 - prop.test(table(copy_loss_data_out_cervix[[i]]))[[6]][2])
}
plot_CL_data_cervix <- data.frame(prop = prop, low_ci = low_ci, high_ci = high_ci, var_id, origin  = origin )

plot_CL_data <- rbind(plot_CL_data_op, plot_CL_data_cervix)


plot_CL_data$var_id <- c("tenth_loss", "E2_loss", "tenth_loss", "E2_loss")

myplt_CN_loss <- ggplot(data = plot_CL_data, aes(x = var_id, y = prop, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  geom_errorbar(aes(x= var_id, y= prop, ymax = high_ci, ymin = low_ci), position = "dodge", width = 0.5)+
  theme_classic() +
  scale_fill_manual(values = c("grey25", "grey75"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = "compare_viral_CN_loss_OP_UC.pdf", height = 2.5, width = 2)
print(myplt_CN_loss)
dev.off()















