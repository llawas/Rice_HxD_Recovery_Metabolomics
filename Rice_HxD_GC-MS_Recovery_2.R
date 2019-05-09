#load required packages and functions
#NOTE: Refer to "Rice_HxD_Metabolomics" repository for the functions
library(VennDiagram)
library(grid)
library(RColorBrewer)
library(pcaMethods)
library(gplots)
library(outliers)
library(stats)
library(nortest)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(pheatmap)
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/RemoveFactors_function.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_normalize.R")
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_find_one_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_replace_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_hist_outlier.R')
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log_transform.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log2_median_transform.R")


#NOTE: early stress means mild stress; late stress means severe stress
#NOTE: for initial data pre-processing (e.g. log and median transformation) and related data, refer to the other GitHub repository


###################################################
#FLOWERING SPIKELET#
###################################################



##################################################
#PRINCIPAL COMPONENT ANALYSIS#
##################################################


#PCA - control, late stress, recovery

#exclude early stress
remove_earlystress = c("FAES", "FDES", "FNES")

data_combined_log2_mean_C_LS_R_flowering_spikelet = data_combined_log2_mean_flowering_spikelet[!(data_combined_log2_mean_flowering_spikelet$Code) %in% remove_earlystress, ]    #data without early stress
data_combined_log2_mean_C_LS_R_flowering_spikelet_2 = data_combined_log2_mean_flowering_spikelet_2[!rownames(data_combined_log2_mean_flowering_spikelet_2) %in% remove_earlystress, ]    #data without early stress

#PCA
pca_res_combined_C_LS_R_flowering_spikelet = pca(data_combined_log2_mean_C_LS_R_flowering_spikelet_2, method = "ppca", nPcs = 5, scale = "pareto", center = TRUE)   
pca_res_scores_combined_C_LS_R_flowering_spikelet = scores(pca_res_combined_C_LS_R_flowering_spikelet)    #scores
pca_res_loadings_combined_C_LS_R_flowering_spikelet = loadings(pca_res_combined_C_LS_R_flowering_spikelet)  #loadings

#scores + additional info as data frame to use in score plot - ggplot
scores_combined_C_LS_R_flowering_spikelet = as.data.frame(pca_res_scores_combined_C_LS_R_flowering_spikelet)
scores_combined_C_LS_R_flowering_spikelet$Timepoint = data_combined_log2_mean_C_LS_R_flowering_spikelet$Timepoint
scores_combined_C_LS_R_flowering_spikelet$Cultivar = data_combined_log2_mean_C_LS_R_flowering_spikelet$Cultivar
scores_combined_C_LS_R_flowering_spikelet$Label = sub(".", "", rownames(scores_combined_C_LS_R_flowering_spikelet))   #remove first character from rowname


#biplot
pdf("biplot_pca_combined_C_LS_R_flowering_spikelet.pdf")
biplot(pca_res_combined_C_LS_R_flowering_spikelet, cex = 0.8, main = "Combined_Control, Late stress, Recovery \nFlowering spikelet")
dev.off()

#screeplot
pdf("screeplot_pca_combined_C_LS_R_flowering_spikelet.pdf")
text(barplot(pca_res_combined_C_LS_R_flowering_spikelet@R2*100, names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), main = "Combined_Control, Late stress, Recovery \nFlowering spikelet", ylab = "Variance (%)", ylim = c(0, 40)), 0, round(pca_res_combined_C_LS_R_flowering_spikelet@R2*100, 3), pos = 3)
box()
dev.off()


#score plot

#symbols
#in assigning colors and shapes to timepoints and cultivars, check the order of factor levels first
#different colors
#with labels below symbols
png("score_plot_pca_combined_C_LS_R_flowering_spikelet.png", width = 8*300, height = 8*300, res = 300)
ggplot(scores_combined_C_LS_R_flowering_spikelet, 
       aes(x = PC1, 
           y = PC2, 
           color = Timepoint, 
           fill = Timepoint, 
           shape = Cultivar, 
           label = Label)) +
  geom_point(size = 7) +
  geom_text(vjust = 2, color = "black") +
  scale_color_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                     values = c("#ffeda0", "#bd0026", "#e31a1c", "#fc4e2a", "#800026"),
                     labels = c("Control", "Late stress", "12h rewatering", "36h rewatering", "60h rewatering")) +
  scale_fill_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                    values = c("#ffeda0", "#bd0026", "#e31a1c", "#fc4e2a", "#800026"),
                    labels = c("Control", "Late stress", "12h rewatering", "36h rewatering", "60h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22", "Dular", "Anjali")) +
  scale_x_continuous(limits = c(-8, 6), breaks = seq(-8, 6, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  xlab(paste("PC1 (", round(pca_res_combined_C_LS_R_flowering_spikelet@R2[1], 4)*100, "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca_res_combined_C_LS_R_flowering_spikelet@R2[2], 4)*100, "%)", sep = "")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.box.just = "top",
        legend.position = c(0.995, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.spacing = unit(25.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)))
dev.off()


#export as table then upload in flag leaf R script to make a multipanel figure of PCAs of all organs
write.table(scores_combined_C_LS_R_flowering_spikelet, file = "PCA_scores_combined_flowering_spikelet.txt", sep = "\t", quote = F)



###################################################
#analysis of rewatering time points
###################################################


#late stress vs. recovery - NOTE: consider only the stress-responsive metabolites (i.e. metabolites which have significant difference in LS/C comparison)


#log2-fold change - late stress vs. recovery

stress_recovery_order = c("FNLS", "FNR12", "FNR36", "FNR60", "FDLS", "FDR12", "FDR36", "FDR60", "FALS", "FAR12", "FAR36", "FAR60")
data_combined_latestress_recovery_log2_mean_flowering_spikelet = data_combined_log2_mean_flowering_spikelet_2[stress_recovery_order, ]  #late stress and recovery data in preferred order
data_combined_latestress_recovery_log2_mean_flowering_spikelet = t(data_combined_latestress_recovery_log2_mean_flowering_spikelet)

log2fc_combined_N22_latestress_recovery_flowering_spikelet = data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 2:4] - data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 1]
log2fc_combined_Dular_latestress_recovery_flowering_spikelet = data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 6:8] - data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 5]
log2fc_combined_Anjali_latestress_recovery_flowering_spikelet = data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 10:12] - data_combined_latestress_recovery_log2_mean_flowering_spikelet[ , 9]

log2fc_combined_latestress_recovery_flowering_spikelet = cbind(log2fc_combined_N22_latestress_recovery_flowering_spikelet, 
                                                               log2fc_combined_Dular_latestress_recovery_flowering_spikelet, 
                                                               log2fc_combined_Anjali_latestress_recovery_flowering_spikelet)   #combine log2-fc of all cultivars in one dataset
colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[1:3] = c(paste(colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[1:3], "- FNLS"))
colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[4:6] = c(paste(colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[4:6], "- FDLS"))
colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[7:9] = c(paste(colnames(log2fc_combined_latestress_recovery_flowering_spikelet)[7:9], "- FALS"))

#overview
heatmap.2(log2fc_combined_latestress_recovery_flowering_spikelet, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-1.5, -0.01, length.out = 300), seq(0.01, 1.5, length.out = 300)), 
          margins = c(7, 20),  
          main = "Flowering spikelet - Combined - Log2 FC, Recovery/Late Stress", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


#late stress vs. recovery

#N22

data_combined_N22_log2median_stress_responsive_flowering_spikelet = data_combined_N22_log2median_flowering_spikelet
colnames(data_combined_N22_log2median_stress_responsive_flowering_spikelet)[3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)] = reduced_metabolite_list_final_flowering_spikelet_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_N22_log2median_stress_responsive_flowering_spikelet = data_combined_N22_log2median_stress_responsive_flowering_spikelet[ , c(1, 2, (which(colnames(data_combined_N22_log2median_stress_responsive_flowering_spikelet) %in% rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined))))]

#data distribution 
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[which(data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet)
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet$Distribution == "Normal")     #37 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_flowering_spikelet$Distribution == "Non-normal") #16 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance == "*") #5 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance == "**") #1 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance == "***") #0 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance != "ns") #6 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance == "ns") #47 metabolites 


#Dular

data_combined_Dular_log2median_stress_responsive_flowering_spikelet = data_combined_Dular_log2median_flowering_spikelet
colnames(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)[3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)] = reduced_metabolite_list_final_flowering_spikelet_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Dular_log2median_stress_responsive_flowering_spikelet = data_combined_Dular_log2median_stress_responsive_flowering_spikelet[ , c(1, 2, (which(colnames(data_combined_Dular_log2median_stress_responsive_flowering_spikelet) %in% rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined))))]

#data distribution 
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[which(data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet$Distribution == "Normal")     #35 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_flowering_spikelet$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined,
                                                                                          ifelse(p.value <= 0.001, "***", ifelse(
                                                                                            p.value <= 0.01, "**", ifelse(
                                                                                              p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance == "*") #11 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance == "**") #3 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance == "***") #2 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance != "ns") #16 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance == "ns") #37 metabolites


#Anjali

data_combined_Anjali_log2median_stress_responsive_flowering_spikelet = data_combined_Anjali_log2median_flowering_spikelet
colnames(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)[3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)] = reduced_metabolite_list_final_flowering_spikelet_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Anjali_log2median_stress_responsive_flowering_spikelet = data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[ , c(1, 2, (which(colnames(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet) %in% rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined))))]

#data distribution 
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[which(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet$Distribution == "Normal")     #33 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_flowering_spikelet$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined,
                                                                                           ifelse(p.value <= 0.001, "***", ifelse(
                                                                                             p.value <= 0.01, "**", ifelse(
                                                                                               p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance == "*") #9 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance == "**") #4 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance == "***") #8 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance == "ns") #32 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined = cbind(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined, 
                                                                        wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance, 
                                                                        wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance)
wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance"
wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined)
wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R12_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined[-which(wilcox_sig_stress_responsive_LS_R12_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 30 metabolites

wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined = log2fc_combined_latestress_recovery_flowering_spikelet[which(rownames(log2fc_combined_latestress_recovery_flowering_spikelet) %in% rownames(wilcox_sig_stress_responsive_final_LS_R12_flowering_spikelet_combined)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R12_flowering_spikelet_combined = wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R12_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$`FNR12 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$Up.Down == "Down")  #4 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R12_flowering_spikelet_combined = wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R12_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$`FDR12 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$Up.Down == "Down")  #14 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined = wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_flowering_spikelet_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$`FAR12 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$Up.Down == "Down")  #19 metabolites


#venn - increased - 12h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R12_flowering_spikelet_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R12_latestress_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_stress_responsive_wilcox_LS_R12_venn_up_flowering_spikelet_combined = attr(venn_up_stress_responsive_wilcox_LS_R12_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined[-1, ]


#venn - decreased - 12h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R12_flowering_spikelet_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R12_latestress_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_stress_responsive_wilcox_LS_R12_venn_down_flowering_spikelet_combined = attr(venn_down_stress_responsive_wilcox_LS_R12_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flowering_spikelet_R12_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flowering_spikelet_R12_latestress_combined.txt", sep = "\t", quote = F)


####################


#late stress vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[which(data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet)
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet$Distribution == "Normal")     #41 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_flowering_spikelet$Distribution == "Non-normal") #12 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance == "*") #6 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance == "**") #1 metabolite
sum(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance == "***") #3 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance != "ns") #10 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance == "ns") #43 metabolites 


#Dular

#data distribution
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[which(data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet$Distribution == "Normal")     #34 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_flowering_spikelet$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined,
                                                                                          ifelse(p.value <= 0.001, "***", ifelse(
                                                                                            p.value <= 0.01, "**", ifelse(
                                                                                              p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance == "*") #7 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance == "**") #5 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance == "***") #7 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance != "ns") #19 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance == "ns") #34 metabolites


#Anjali

#data distribution
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[which(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet$Distribution == "Normal")     #37 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_flowering_spikelet$Distribution == "Non-normal") #16 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined,
                                                                                           ifelse(p.value <= 0.001, "***", ifelse(
                                                                                             p.value <= 0.01, "**", ifelse(
                                                                                               p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance == "*") #6 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance == "**") #2 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance == "***") #7 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance != "ns") #15 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance == "ns") #38 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined = cbind(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined, 
                                                                        wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance, 
                                                                        wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance)
wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance"
wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined)
wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R36_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined[-which(wilcox_sig_stress_responsive_LS_R36_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 22 metabolites

wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined = log2fc_combined_latestress_recovery_flowering_spikelet[which(rownames(log2fc_combined_latestress_recovery_flowering_spikelet) %in% rownames(wilcox_sig_stress_responsive_final_LS_R36_flowering_spikelet_combined)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R36_flowering_spikelet_combined = wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R36_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$`FNR36 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$Up.Down == "Down")  #9 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R36_flowering_spikelet_combined = wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R36_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$`FDR36 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$Up.Down == "Down")  #17 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined = wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_flowering_spikelet_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$`FAR36 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$Up.Down == "Up")    #3 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$Up.Down == "Down")  #12 metabolites


#venn - increased - 36h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R36_flowering_spikelet_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R36_latestress_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_stress_responsive_wilcox_LS_R36_venn_up_flowering_spikelet_combined = attr(venn_up_stress_responsive_wilcox_LS_R36_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined[-1, ]


#venn - decreased - 36h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R36_flowering_spikelet_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R36_latestress_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_stress_responsive_wilcox_LS_R36_venn_down_flowering_spikelet_combined = attr(venn_down_stress_responsive_wilcox_LS_R36_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flowering_spikelet_R36_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flowering_spikelet_R36_latestress_combined.txt", sep = "\t", quote = F)


####################


#late stress vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[which(data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet)
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet$Distribution == "Normal")     #43 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_flowering_spikelet$Distribution == "Non-normal") #10 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined)
wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance == "*") #8 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance == "**") #2 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance == "***") #4 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance != "ns") #14 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance == "ns") #39 metabolites 


#Dular

#data distribution
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[which(data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet$Distribution == "Normal")     #40 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_flowering_spikelet$Distribution == "Non-normal") #13 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined)
wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined,
                                                                                          ifelse(p.value <= 0.001, "***", ifelse(
                                                                                            p.value <= 0.01, "**", ifelse(
                                                                                              p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance == "*") #4 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance == "**") #3 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance == "***") #10 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance != "ns") #17 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance == "ns") #36 metabolites


#Anjali

#data distribution
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[which(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet$Distribution == "Normal")     #39 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_flowering_spikelet$Distribution == "Non-normal") #14 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_flowering_spikelet$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined) = rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined)
wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined,
                                                                                           ifelse(p.value <= 0.001, "***", ifelse(
                                                                                             p.value <= 0.01, "**", ifelse(
                                                                                               p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance == "*") #5 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance == "**") #9 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance == "***") #10 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance != "ns") #24 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance == "ns") #29 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined = cbind(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined, 
                                                                        wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance, 
                                                                        wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance)
wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance"
wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined)
wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R60_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined[-which(wilcox_sig_stress_responsive_LS_R60_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 30 metabolites

wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined = log2fc_combined_latestress_recovery_flowering_spikelet[which(rownames(log2fc_combined_latestress_recovery_flowering_spikelet) %in% rownames(wilcox_sig_stress_responsive_final_LS_R60_flowering_spikelet_combined)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R60_flowering_spikelet_combined = wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R60_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$`FNR60 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$Up.Down == "Down")  #14 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R60_flowering_spikelet_combined = wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R60_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$`FDR60 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$Up.Down == "Down")  #15 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined = wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_flowering_spikelet_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$`FAR60 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$Up.Down == "Down")  #24 metabolites


#venn - increased - 60h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R60_flowering_spikelet_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R60_latestress_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_stress_responsive_wilcox_LS_R60_venn_up_flowering_spikelet_combined = attr(venn_up_stress_responsive_wilcox_LS_R60_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined[-1, ]


#venn - decreased - 60h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R60_flowering_spikelet_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R60_latestress_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_stress_responsive_wilcox_LS_R60_venn_down_flowering_spikelet_combined = attr(venn_down_stress_responsive_wilcox_LS_R60_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flowering_spikelet_R60_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flowering_spikelet_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flowering_spikelet_R60_latestress_combined.txt", sep = "\t", quote = F)


####################


#significance of wilcoxon test - all cultivars, all three RW timepoints
wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined = cbind(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined,
                                                                                wilcox_res_stress_responsive_N22_LS_R36_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_N22_LS_R60_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Dular_LS_R12_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Dular_LS_R36_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Dular_LS_R60_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Anjali_LS_R12_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Anjali_LS_R36_flowering_spikelet_combined$Significance,
                                                                                wilcox_res_stress_responsive_Anjali_LS_R60_flowering_spikelet_combined$Significance)
wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined$Significance"
wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R12_flowering_spikelet_combined)
wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R12_R36_R60_flowering_spikelet_combined = wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined[-which(wilcox_sig_stress_responsive_LS_R12_R36_R60_flowering_spikelet_combined$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 41 metabolites

#log2fc values of metabolites with significant difference in any one of the comparisons - 12, 36, 60h recovery vs. late stress in any cultivar
wilcox_sig_stress_responsive_log2fc_LS_R12_R36_R60_flowering_spikelet_combined = log2fc_combined_latestress_recovery_flowering_spikelet[which(rownames(log2fc_combined_latestress_recovery_flowering_spikelet) %in% rownames(wilcox_sig_stress_responsive_final_LS_R12_R36_R60_flowering_spikelet_combined)), ]


#visualization
png("wilcox_sig_stress_responsive_log2fc_LS_RW_flowering_spikelet_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_stress_responsive_log2fc_LS_R12_R36_R60_flowering_spikelet_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-1.5, -0.01, length.out = 300), seq(0.01, 1.5, length.out = 300)), 
          margins = c(7.5, 23),  
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("12 h RW/Severe stress", "36 h RW/Severe stress", "60 h RW/Severe stress"), 3),
          cexCol = 1, 
          cexRow = 1.2, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 1.5),
          srtCol = 45, 
          cellnote = wilcox_sig_stress_responsive_final_LS_R12_R36_R60_flowering_spikelet_combined, 
          notecol = "black", 
          notecex = 1.8, 
          sepcolor = "black", 
          colsep = c(seq(0, 9)), 
          rowsep = c(seq(0, 41)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 3), rep("#fc8d62", 3), rep("#8da0cb", 3)),
          hclustfun = function(x) hclust(x, method = "average"),
          key.xtickfun = function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 1
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("-1.5", side=side, at=-0.15, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("0", side=side, at=0.5, adj=0.5,
                  line=line, cex=cex, col=col, font=font)
            mtext("1.5", side=side, at=1.1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=TRUE))
          })
legend(x = -0.009, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.105, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.22, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.865, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
legend(x = -0.15, y = 1.13, xpd = TRUE, legend = c("G"), bty = "n", cex = 1.05, text.font = 2)
dev.off()


####################

#all venn diagrams

#multiplot
png("venn_up_down_flowering_spikelet_stress_responsive_RW_latestress_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_stress_responsive_R12_latestress_flowering_spikelet_combined), 
             gTree(children = venn_down_stress_responsive_R12_latestress_flowering_spikelet_combined), 
             gTree(children = venn_up_stress_responsive_R36_latestress_flowering_spikelet_combined),
             gTree(children = venn_down_stress_responsive_R36_latestress_flowering_spikelet_combined),
             gTree(children = venn_up_stress_responsive_R60_latestress_flowering_spikelet_combined),
             gTree(children = venn_down_stress_responsive_R60_latestress_flowering_spikelet_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#recovery vs. corresponding control - NOTE: consider all metabolites, not just the stress-responsive


#log2-fold change - control vs. rewatering
control_recovery_order = c("FNC", "FNR12", "FNR36", "FNR60", "FDC", "FDR12", "FDR36", "FDR60", "FAC", "FAR12", "FAR36", "FAR60") #line command already ran previously
data_combined_control_recovery_log2_mean_flowering_spikelet = data_combined_log2_mean_flowering_spikelet_2[control_recovery_order, ]  #control and recovery data in preferred order
data_combined_control_recovery_log2_mean_flowering_spikelet = t(data_combined_control_recovery_log2_mean_flowering_spikelet)

log2fc_combined_N22_control_recovery_flowering_spikelet = data_combined_control_recovery_log2_mean_flowering_spikelet[ , 2:4] - data_combined_control_recovery_log2_mean_flowering_spikelet[ , 1]
log2fc_combined_Dular_control_recovery_flowering_spikelet = data_combined_control_recovery_log2_mean_flowering_spikelet[ , 6:8] - data_combined_control_recovery_log2_mean_flowering_spikelet[ , 5]
log2fc_combined_Anjali_control_recovery_flowering_spikelet = data_combined_control_recovery_log2_mean_flowering_spikelet[ , 10:12] - data_combined_control_recovery_log2_mean_flowering_spikelet[ , 9]

#combine in one data set
log2fc_combined_control_recovery_flowering_spikelet = cbind(log2fc_combined_N22_control_recovery_flowering_spikelet, 
                                                            log2fc_combined_Dular_control_recovery_flowering_spikelet, 
                                                            log2fc_combined_Anjali_control_recovery_flowering_spikelet)   #combine log2-fc of all cultivars in one dataset
#rename columns
colnames(log2fc_combined_control_recovery_flowering_spikelet)[1:3] = c(paste(colnames(log2fc_combined_control_recovery_flowering_spikelet)[1:3], "- FNC"))
colnames(log2fc_combined_control_recovery_flowering_spikelet)[4:6] = c(paste(colnames(log2fc_combined_control_recovery_flowering_spikelet)[4:6], "- FDC"))
colnames(log2fc_combined_control_recovery_flowering_spikelet)[7:9] = c(paste(colnames(log2fc_combined_control_recovery_flowering_spikelet)[7:9], "- FAC"))

#overview
#heatmap including all metabolites
heatmap.2(log2fc_combined_control_recovery_flowering_spikelet, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(6, 20), 
          main = "Flowering spikelet - Combined - Log2 FC, Recovery/Control",
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA,
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


####################

#control vs. 12h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[which(data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_N22_control_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_N22_control_recovery12_flowering_spikelet)
shapiro_res_combined_N22_control_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_N22_control_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery12_flowering_spikelet$Distribution == "Normal")     #62 metabolites
sum(shapiro_res_combined_N22_control_recovery12_flowering_spikelet$Distribution == "Non-normal") #26 metabolites

#wilcox test
wilcox_res_N22_C_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_N22_C_R12_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_N22_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_N22_C_R12_flowering_spikelet_combined)
wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance = with(wilcox_res_N22_C_R12_flowering_spikelet_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance == "*") #8 metabolites
sum(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance == "**") #4 metabolites
sum(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance == "***") #9 metabolites
sum(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance == "ns") #67 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[which(data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_Dular_control_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_Dular_control_recovery12_flowering_spikelet)
shapiro_res_combined_Dular_control_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery12_flowering_spikelet$Distribution == "Normal")     #62 metabolites
sum(shapiro_res_combined_Dular_control_recovery12_flowering_spikelet$Distribution == "Non-normal") #26 metabolites

#wilcox test
wilcox_res_Dular_C_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Dular_C_R12_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Dular_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_Dular_C_R12_flowering_spikelet_combined)
wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance = with(wilcox_res_Dular_C_R12_flowering_spikelet_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance == "*") #14 metabolites
sum(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance == "***") #20 metabolites
sum(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance != "ns") #42 metabolites
sum(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance == "ns") #46 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[which(data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet) = "p.value"
shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet = as.data.frame(shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet)
shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet$Distribution == "Normal")     #46 metabolites
sum(shapiro_res_combined_Anjali_control_recovery12_flowering_spikelet$Distribution == "Non-normal") #42 metabolites

#wilcox test
wilcox_res_Anjali_C_R12_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R12_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R12_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Anjali_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_res_Anjali_C_R12_flowering_spikelet_combined)
wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance = with(wilcox_res_Anjali_C_R12_flowering_spikelet_combined,
                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                          p.value <= 0.01, "**", ifelse(
                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance == "***") #19 metabolites
sum(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance != "ns") #34 metabolites
sum(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance == "ns") #54 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R12_flowering_spikelet_combined = cbind(wilcox_res_N22_C_R12_flowering_spikelet_combined, 
                                                     wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance, 
                                                     wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance)
wilcox_sig_C_R12_flowering_spikelet_combined = wilcox_sig_C_R12_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_C_R12_flowering_spikelet_combined)[1] = "wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance"
wilcox_sig_C_R12_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_C_R12_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R12_flowering_spikelet_combined) = rownames(wilcox_res_N22_C_R12_flowering_spikelet_combined)
wilcox_sig_C_R12_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_C_R12_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R12_flowering_spikelet_combined = wilcox_sig_C_R12_flowering_spikelet_combined[-which(wilcox_sig_C_R12_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 62 metabolites

wilcox_sig_log2fc_C_R12_flowering_spikelet_combined = log2fc_combined_control_recovery_flowering_spikelet[which(rownames(log2fc_combined_control_recovery_flowering_spikelet) %in% rownames(wilcox_sig_final_C_R12_flowering_spikelet_combined)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(2, 20),  
          main = "Flowering spikelet-Combined-Wilcoxon \nLog2 FC, FL, 12h RW/Control", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_final_C_R12_flowering_spikelet_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 62)), 
          sepwidth = c(0.01, 0.01),
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R12_flowering_spikelet_combined = wilcox_res_N22_C_R12_flowering_spikelet_combined[-which(wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_N22_C_R12_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)[1]
wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$`FNR12 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$Up.Down == "Up")    #15 metabolites
sum(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$Up.Down == "Down")  #6 metabolites

#Dular
wilcox_sig_Dular_C_R12_flowering_spikelet_combined = wilcox_res_Dular_C_R12_flowering_spikelet_combined[-which(wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_Dular_C_R12_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)[2]
wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$`FDR12 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$Up.Down == "Up")    #24 metabolites
sum(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$Up.Down == "Down")  #18 metabolites

#Anjali
wilcox_sig_Anjali_C_R12_flowering_spikelet_combined = wilcox_res_Anjali_C_R12_flowering_spikelet_combined[-which(wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined[match(rownames(wilcox_sig_Anjali_C_R12_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R12_flowering_spikelet_combined)[3]
wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$`FAR12 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$Up.Down == "Up")    #16 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$Up.Down == "Down")  #18 metabolites


#venn - increased - 12h rewatering/control
vennlist_up_wilcox_C_R12_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R12_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R12_flowering_spikelet_combined = venn(vennlist_up_wilcox_C_R12_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R12_control_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_C_R12_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R12_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_C_R12_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R12_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_C_R12_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R12_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_C_R12_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_wilcox_C_R12_venn_up_flowering_spikelet_combined = attr(venn_up_wilcox_C_R12_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R12_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined) = venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined[1, ]
venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined = venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined[-1, ]


#venn - decreased - 12h rewatering/control
vennlist_down_wilcox_C_R12_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R12_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R12_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R12_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R12_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R12_flowering_spikelet_combined = venn(vennlist_down_wilcox_C_R12_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R12_control_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_C_R12_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R12_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_C_R12_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R12_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_C_R12_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R12_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_C_R12_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_wilcox_C_R12_venn_down_flowering_spikelet_combined = attr(venn_down_wilcox_C_R12_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R12_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined) = venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined[1, ]
venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined = venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R12_flowering_spikelet_combined, file = "wilcox_sig_up_metabolites_flowering_spikelet_R12_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R12_flowering_spikelet_combined, file = "wilcox_sig_down_metabolites_flowering_spikelet_R12_control_combined.txt", sep = "\t", quote = F)


####################


#control vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[which(data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_N22_control_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_N22_control_recovery36_flowering_spikelet)
shapiro_res_combined_N22_control_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_N22_control_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery36_flowering_spikelet$Distribution == "Normal")     #63 metabolites
sum(shapiro_res_combined_N22_control_recovery36_flowering_spikelet$Distribution == "Non-normal") #25 metabolites

#wilcox test
wilcox_res_N22_C_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_N22_C_R36_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_N22_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_N22_C_R36_flowering_spikelet_combined)
wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance = with(wilcox_res_N22_C_R36_flowering_spikelet_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance == "*") #9 metabolites
sum(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance == "***") #10 metabolites
sum(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance != "ns") #24 metabolites
sum(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance == "ns") #64 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[which(data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_Dular_control_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_Dular_control_recovery36_flowering_spikelet)
shapiro_res_combined_Dular_control_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery36_flowering_spikelet$Distribution == "Normal")     #53 metabolites
sum(shapiro_res_combined_Dular_control_recovery36_flowering_spikelet$Distribution == "Non-normal") #35 metabolites

#wilcox test
wilcox_res_Dular_C_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Dular_C_R36_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Dular_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_Dular_C_R36_flowering_spikelet_combined)
wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance = with(wilcox_res_Dular_C_R36_flowering_spikelet_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance == "*") #13 metabolites
sum(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance == "***") #16 metabolites
sum(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance != "ns") #37 metabolites
sum(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance == "ns") #51 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[which(data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet) = "p.value"
shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet = as.data.frame(shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet)
shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet$Distribution == "Normal")     #59 metabolites
sum(shapiro_res_combined_Anjali_control_recovery36_flowering_spikelet$Distribution == "Non-normal") #29 metabolites

#wilcox test
wilcox_res_Anjali_C_R36_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R36_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R36_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Anjali_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_res_Anjali_C_R36_flowering_spikelet_combined)
wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance = with(wilcox_res_Anjali_C_R36_flowering_spikelet_combined,
                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                          p.value <= 0.01, "**", ifelse(
                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance == "*") #10 metabolites
sum(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance == "**") #14 metabolites
sum(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance == "***") #16 metabolites
sum(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance != "ns") #40 metabolites
sum(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance == "ns") #48 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R36_flowering_spikelet_combined = cbind(wilcox_res_N22_C_R36_flowering_spikelet_combined, 
                                                     wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance, 
                                                     wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance)
wilcox_sig_C_R36_flowering_spikelet_combined = wilcox_sig_C_R36_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_C_R36_flowering_spikelet_combined)[1] = "wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance"
wilcox_sig_C_R36_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_C_R36_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R36_flowering_spikelet_combined) = rownames(wilcox_res_N22_C_R36_flowering_spikelet_combined)
wilcox_sig_C_R36_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_C_R36_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R36_flowering_spikelet_combined = wilcox_sig_C_R36_flowering_spikelet_combined[-which(wilcox_sig_C_R36_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 60 metabolites

wilcox_sig_log2fc_C_R36_flowering_spikelet_combined = log2fc_combined_control_recovery_flowering_spikelet[which(rownames(log2fc_combined_control_recovery_flowering_spikelet) %in% rownames(wilcox_sig_final_C_R36_flowering_spikelet_combined)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined,
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(2, 20), 
          main = "Flowering spikelet-Combined-Wilcoxon \nLog2 FC, FL, 36h RW/Control", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          cellnote = wilcox_sig_final_C_R36_flowering_spikelet_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 60)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R36_flowering_spikelet_combined = wilcox_res_N22_C_R36_flowering_spikelet_combined[-which(wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_N22_C_R36_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)[1]
wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$`FNR36 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$Up.Down == "Up")    #18 metabolites
sum(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$Up.Down == "Down")  #6 metabolites

#Dular
wilcox_sig_Dular_C_R36_flowering_spikelet_combined = wilcox_res_Dular_C_R36_flowering_spikelet_combined[-which(wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_Dular_C_R36_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)[2]
wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$`FDR36 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$Up.Down == "Up")    #19 metabolites
sum(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$Up.Down == "Down")  #18 metabolites

#Anjali
wilcox_sig_Anjali_C_R36_flowering_spikelet_combined = wilcox_res_Anjali_C_R36_flowering_spikelet_combined[-which(wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined[match(rownames(wilcox_sig_Anjali_C_R36_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R36_flowering_spikelet_combined)[3]
wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$`FAR36 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$Up.Down == "Up")    #26 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$Up.Down == "Down")  #14 metabolites


#venn - increased - 36h rewatering/control
vennlist_up_wilcox_C_R36_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R36_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R36_flowering_spikelet_combined = venn(vennlist_up_wilcox_C_R36_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R36_control_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_C_R36_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R36_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_C_R36_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R36_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_C_R36_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R36_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_C_R36_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_wilcox_C_R36_venn_up_flowering_spikelet_combined = attr(venn_up_wilcox_C_R36_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R36_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined) = venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined[1, ]
venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined = venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined[-1, ]


#venn - decreased - 36h rewatering/control
vennlist_down_wilcox_C_R36_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R36_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R36_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R36_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R36_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R36_flowering_spikelet_combined = venn(vennlist_down_wilcox_C_R36_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R36_control_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_C_R36_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R36_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_C_R36_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R36_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_C_R36_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R36_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_C_R36_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_wilcox_C_R36_venn_down_flowering_spikelet_combined = attr(venn_down_wilcox_C_R36_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R36_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined) = venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined[1, ]
venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined = venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R36_flowering_spikelet_combined, file = "wilcox_sig_up_metabolites_flowering_spikelet_R36_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R36_flowering_spikelet_combined, file = "wilcox_sig_down_metabolites_flowering_spikelet_R36_control_combined.txt", sep = "\t", quote = F)


####################


#Control vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[which(data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_N22_control_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_N22_control_recovery60_flowering_spikelet)
shapiro_res_combined_N22_control_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_N22_control_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery60_flowering_spikelet$Distribution == "Normal")     #61 metabolites
sum(shapiro_res_combined_N22_control_recovery60_flowering_spikelet$Distribution == "Non-normal") #27 metabolites

#wilcox test
wilcox_res_N22_C_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_N22_log2median_flowering_spikelet[ , 3:ncol(data_combined_N22_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_N22_log2median_flowering_spikelet$Timepoint, subset = data_combined_N22_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_N22_C_R60_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_N22_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_N22_C_R60_flowering_spikelet_combined)
wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance = with(wilcox_res_N22_C_R60_flowering_spikelet_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance == "*") #16 metabolites
sum(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance == "**") #4 metabolites
sum(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance == "***") #6 metabolites
sum(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance != "ns") #26 metabolites
sum(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance == "ns") #62 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[which(data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_Dular_control_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_Dular_control_recovery60_flowering_spikelet)
shapiro_res_combined_Dular_control_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery60_flowering_spikelet$Distribution == "Normal")     #56 metabolites
sum(shapiro_res_combined_Dular_control_recovery60_flowering_spikelet$Distribution == "Non-normal") #32 metabolites

#wilcox test
wilcox_res_Dular_C_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Dular_log2median_flowering_spikelet[ , 3:ncol(data_combined_Dular_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_flowering_spikelet$Timepoint, subset = data_combined_Dular_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Dular_C_R60_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Dular_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_Dular_C_R60_flowering_spikelet_combined)
wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance = with(wilcox_res_Dular_C_R60_flowering_spikelet_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance == "*") #13 metabolites
sum(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance == "**") #14 metabolites
sum(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance == "***") #13 metabolites
sum(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance != "ns") #40 metabolites
sum(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance == "ns") #48 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[which(data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet) = "p.value"
shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet = as.data.frame(shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet)
shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet$Distribution == "Normal")     #47 metabolites
sum(shapiro_res_combined_Anjali_control_recovery60_flowering_spikelet$Distribution == "Non-normal") #41 metabolites

#wilcox test
wilcox_res_Anjali_C_R60_flowering_spikelet_combined = t(data.frame(lapply(data_combined_Anjali_log2median_flowering_spikelet[ , 3:ncol(data_combined_Anjali_log2median_flowering_spikelet)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_flowering_spikelet$Timepoint, subset = data_combined_Anjali_log2median_flowering_spikelet$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R60_flowering_spikelet_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R60_flowering_spikelet_combined) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
wilcox_res_Anjali_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_res_Anjali_C_R60_flowering_spikelet_combined)
wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance = with(wilcox_res_Anjali_C_R60_flowering_spikelet_combined,
                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                          p.value <= 0.01, "**", ifelse(
                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance == "*") #19 metabolites
sum(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance == "**") #13 metabolites
sum(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance == "***") #14 metabolites
sum(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance != "ns") #46 metabolites
sum(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance == "ns") #42 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R60_flowering_spikelet_combined = cbind(wilcox_res_N22_C_R60_flowering_spikelet_combined,
                                                     wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance, 
                                                     wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance)
wilcox_sig_C_R60_flowering_spikelet_combined = wilcox_sig_C_R60_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_C_R60_flowering_spikelet_combined)[1] = "wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance"
wilcox_sig_C_R60_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_C_R60_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R60_flowering_spikelet_combined) = rownames(wilcox_res_N22_C_R60_flowering_spikelet_combined)
wilcox_sig_C_R60_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_C_R60_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R60_flowering_spikelet_combined = wilcox_sig_C_R60_flowering_spikelet_combined[-which(wilcox_sig_C_R60_flowering_spikelet_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 65 metabolites

wilcox_sig_log2fc_C_R60_flowering_spikelet_combined = log2fc_combined_control_recovery_flowering_spikelet[which(rownames(log2fc_combined_control_recovery_flowering_spikelet) %in% rownames(wilcox_sig_final_C_R60_flowering_spikelet_combined)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(2, 20), 
          main = "Flowering spikelet-Combined-Wilcoxon \nLog2 FC, FL, 60h RW/Control", 
          cexCol = 1.5,
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          cellnote = wilcox_sig_final_C_R60_flowering_spikelet_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)),
          rowsep = c(seq(0, 65)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R60_flowering_spikelet_combined = wilcox_res_N22_C_R60_flowering_spikelet_combined[-which(wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_N22_C_R60_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)[1]
wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$`FNR60 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$Up.Down == "Up")    #12 metabolites
sum(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$Up.Down == "Down")  #14 metabolites

#Dular
wilcox_sig_Dular_C_R60_flowering_spikelet_combined = wilcox_res_Dular_C_R60_flowering_spikelet_combined[-which(wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_Dular_C_R60_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)[2]
wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$`FDR60 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$Up.Down == "Up")    #18 metabolites
sum(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$Up.Down == "Down")  #22 metabolites

#Anjali
wilcox_sig_Anjali_C_R60_flowering_spikelet_combined = wilcox_res_Anjali_C_R60_flowering_spikelet_combined[-which(wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined[match(rownames(wilcox_sig_Anjali_C_R60_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined) = colnames(wilcox_sig_log2fc_C_R60_flowering_spikelet_combined)[3]
wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$`FAR60 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$Up.Down == "Up")    #11 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$Up.Down == "Down")  #35 metabolites


#venn - increased - 60h rewatering/control
vennlist_up_wilcox_C_R60_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R60_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R60_flowering_spikelet_combined = venn(vennlist_up_wilcox_C_R60_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R60_control_flowering_spikelet_combined = venn.diagram(vennlist_up_wilcox_C_R60_flowering_spikelet_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R60_flowering_spikelet_combined[1]), " (", length(vennlist_up_wilcox_C_R60_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R60_flowering_spikelet_combined[2]), " (", length(vennlist_up_wilcox_C_R60_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R60_flowering_spikelet_combined[3]), " (", length(vennlist_up_wilcox_C_R60_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_wilcox_C_R60_venn_up_flowering_spikelet_combined = attr(venn_up_wilcox_C_R60_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R60_venn_up_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined) = venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined[1, ]
venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined = venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined[-1, ]


#venn - decreased - 60h rewatering/control
vennlist_down_wilcox_C_R60_flowering_spikelet_combined = list(rownames(wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_N22_C_R60_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_Dular_C_R60_flowering_spikelet_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined)[wilcox_sig_log2fc_Anjali_C_R60_flowering_spikelet_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R60_flowering_spikelet_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R60_flowering_spikelet_combined = venn(vennlist_down_wilcox_C_R60_flowering_spikelet_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R60_control_flowering_spikelet_combined = venn.diagram(vennlist_down_wilcox_C_R60_flowering_spikelet_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R60_flowering_spikelet_combined[1]), " (", length(vennlist_down_wilcox_C_R60_flowering_spikelet_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R60_flowering_spikelet_combined[2]), " (", length(vennlist_down_wilcox_C_R60_flowering_spikelet_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R60_flowering_spikelet_combined[3]), " (", length(vennlist_down_wilcox_C_R60_flowering_spikelet_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_wilcox_C_R60_venn_down_flowering_spikelet_combined = attr(venn_down_wilcox_C_R60_flowering_spikelet_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined = t(ldply(intersect_wilcox_C_R60_venn_down_flowering_spikelet_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined) = venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined[1, ]
venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined = venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R60_flowering_spikelet_combined, file = "wilcox_sig_up_metabolites_flowering_spikelet_R60_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R60_flowering_spikelet_combined, file = "wilcox_sig_down_metabolites_flowering_spikelet_R60_control_combined.txt", sep = "\t", quote = F)


####################


#rownames were checked to be the same for the data files to be combined
#significance of wilcoxon test - all cultivars, all three RW timepoints and LS/C (refer to previous paper and associated GitHub files)
wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined = cbind(wilcox_res_N22_C_LS_flowering_spikelet_combined,
                                                                wilcox_res_N22_C_R12_flowering_spikelet_combined$Significance,
                                                                wilcox_res_N22_C_R36_flowering_spikelet_combined$Significance,
                                                                wilcox_res_N22_C_R60_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Dular_C_LS_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Dular_C_R12_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Dular_C_R36_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Dular_C_R60_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Anjali_C_LS_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Anjali_C_R12_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Anjali_C_R36_flowering_spikelet_combined$Significance,
                                                                wilcox_res_Anjali_C_R60_flowering_spikelet_combined$Significance)
wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined = wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined[ , -1]
colnames(wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined)[1] = "wilcox_res_N22_C_LS_flowering_spikelet_combined$Significance"
wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined = data.frame(lapply(wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined) = rownames(wilcox_res_N22_C_LS_flowering_spikelet_combined)
wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined$Nonsig = rowSums(wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_LS_R12_R36_R60_flowering_spikelet_combined = wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined[-which(wilcox_sig_C_LS_R12_R36_R60_flowering_spikelet_combined$Nonsig == "12"), ] #metabolites significant in at least one of the comparisons - 82 metabolites

#log2fc values of late stress and recovery relative to control (with reference to stress paper and associated GitHub files)
log2fc_combined_control_latestress_recovery_flowering_spikelet = cbind(log2fc_combined_control_latestress_flowering_spikelet, 
                                                                       log2fc_combined_control_recovery_flowering_spikelet)
log2fc_combined_control_latestress_recovery_flowering_spikelet = log2fc_combined_control_latestress_recovery_flowering_spikelet[ , c(1, 4:6, 2, 7:9, 3, 10:12)]   #rearrange in preferred order

#log2fc values of metabolites with significant difference in at least one of the comparisons - late stress, 12, 36, 60h recovery vs. control in any cultivar
wilcox_sig_log2fc_C_LS_R12_R36_R60_flowering_spikelet_combined = log2fc_combined_control_latestress_recovery_flowering_spikelet[which(rownames(log2fc_combined_control_latestress_recovery_flowering_spikelet) %in% rownames(wilcox_sig_final_C_LS_R12_R36_R60_flowering_spikelet_combined)), ]


#rewatering-responsive metabolites only - indicate as red font in heat map
setdiff(rownames(wilcox_sig_log2fc_C_LS_R12_R36_R60_flowering_spikelet_combined), rownames(wilcox_sig_log2fc_C_LS_flowering_spikelet_combined))

#visualization
png("wilcox_sig_log2fc_C_LS_RW_flowering_spikelet_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_LS_R12_R36_R60_flowering_spikelet_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 18.5),
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("Severe stress/Control", "12 h RW/Control", "36 h RW/Control", "60 h RW/Control"), 3),
          cexCol = 1, 
          cexRow = 0.8, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(1, 5, 1.6), 
          srtCol = 45, 
          cellnote = wilcox_sig_final_C_LS_R12_R36_R60_flowering_spikelet_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 82)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 4), rep("#fc8d62", 4), rep("#8da0cb", 4)),
          hclustfun = function(x) hclust(x, method = "average"),
          colRow = c("black", rep("red", 3), rep("black", 2), rep("red", 4), "black", rep("red", 2), rep("black", 11), rep("red", 3),
                     "black", rep("red", 2), "black", "red", rep("black", 2), "red", rep("black", 5), "red", rep("black", 4), "red",
                     rep("black", 4), "red", rep("black", 4), rep("red", 2), "black", "red", rep("black", 2), "red", rep("black", 2),
                     rep("red", 3), "black", rep("red", 2), rep("black", 2), "red", rep("black", 9)))
legend(x = 0.026, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.155, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.295, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.863, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
dev.off()


####################

#venn diagrams

#multiplot
png("venn_up_down_flowering_spikelet_RW_control_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_R12_control_flowering_spikelet_combined), 
             gTree(children = venn_down_R12_control_flowering_spikelet_combined), 
             gTree(children = venn_up_R36_control_flowering_spikelet_combined),
             gTree(children = venn_down_R36_control_flowering_spikelet_combined),
             gTree(children = venn_up_R60_control_flowering_spikelet_combined),
             gTree(children = venn_down_R60_control_flowering_spikelet_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#correlation analysis

#correlation between change in yield and metabolite changes during 60h RW
#correlation between change in proportion of chalky grains (>50% chalk content) and metabolite changes during 60h RW

#only 60h RW was considered since it is the latest sample collection time point


#yield and chalkiness - mean per cultivar per treatment per year already ran previously


#recovery vs. late stress - log2FC

#2013
stress_recovery_order = c("FNLS", "FNR12", "FNR36", "FNR60", "FDLS", "FDR12", "FDR36", "FDR60", "FALS", "FAR12", "FAR36", "FAR60")
data_stress_recovery_log2_mean_flowering_spikelet_2013 = data_log2_mean_flowering_spikelet_2013_2[stress_recovery_order, ]  #data with stress and recovery values only
data_stress_recovery_log2_mean_flowering_spikelet_2013 = t(data_stress_recovery_log2_mean_flowering_spikelet_2013)  #samples in columns
log2fc_N22_stress_recovery_flowering_spikelet_2013 = data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 2:4] - data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 1]
log2fc_Dular_stress_recovery_flowering_spikelet_2013 = data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 6:8] - data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 5]
log2fc_Anjali_stress_recovery_flowering_spikelet_2013 = data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 10:12] - data_stress_recovery_log2_mean_flowering_spikelet_2013[ , 9]
log2fc_stress_recovery_flowering_spikelet_2013 = cbind(log2fc_N22_stress_recovery_flowering_spikelet_2013, log2fc_Dular_stress_recovery_flowering_spikelet_2013, log2fc_Anjali_stress_recovery_flowering_spikelet_2013)
colnames(log2fc_stress_recovery_flowering_spikelet_2013) = c("FNR12 - FNLS", "FNR36 - FNLS", "FNR60 - FNLS", "FDR12 - FDLS", "FDR36 - FDLS", "FDR60 - FDLS", "FAR12 - FALS", "FAR36 - FALS", "FAR60 - FALS")

#2014
#recovery vs. late stress
data_stress_recovery_log2_mean_flowering_spikelet_2014 = data_log2_mean_flowering_spikelet_2014_2[stress_recovery_order, ]  #data with stress and recovery values only
data_stress_recovery_log2_mean_flowering_spikelet_2014 = t(data_stress_recovery_log2_mean_flowering_spikelet_2014)  #samples in columns
log2fc_N22_stress_recovery_flowering_spikelet_2014 = data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 2:4] - data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 1]
log2fc_Dular_stress_recovery_flowering_spikelet_2014 = data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 6:8] - data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 5]
log2fc_Anjali_stress_recovery_flowering_spikelet_2014 = data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 10:12] - data_stress_recovery_log2_mean_flowering_spikelet_2014[ , 9]
log2fc_stress_recovery_flowering_spikelet_2014 = cbind(log2fc_N22_stress_recovery_flowering_spikelet_2014, log2fc_Dular_stress_recovery_flowering_spikelet_2014, log2fc_Anjali_stress_recovery_flowering_spikelet_2014)
colnames(log2fc_stress_recovery_flowering_spikelet_2014) = c("FNR12 - FNLS", "FNR36 - FNLS", "FNR60 - FNLS", "FDR12 - FDLS", "FDR36 - FDLS", "FDR60 - FDLS", "FAR12 - FALS", "FAR36 - FALS", "FAR60 - FALS")

#2015
#recovery vs. late stress
data_stress_recovery_log2_mean_flowering_spikelet_2015 = data_log2_mean_flowering_spikelet_2015_2[stress_recovery_order, ]  #data with stress and recovery values only
data_stress_recovery_log2_mean_flowering_spikelet_2015 = t(data_stress_recovery_log2_mean_flowering_spikelet_2015)  #samples in columns
log2fc_N22_stress_recovery_flowering_spikelet_2015 = data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 2:4] - data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 1]
log2fc_Dular_stress_recovery_flowering_spikelet_2015 = data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 6:8] - data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 5]
log2fc_Anjali_stress_recovery_flowering_spikelet_2015 = data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 10:12] - data_stress_recovery_log2_mean_flowering_spikelet_2015[ , 9]
log2fc_stress_recovery_flowering_spikelet_2015 = cbind(log2fc_N22_stress_recovery_flowering_spikelet_2015, log2fc_Dular_stress_recovery_flowering_spikelet_2015, log2fc_Anjali_stress_recovery_flowering_spikelet_2015)
colnames(log2fc_stress_recovery_flowering_spikelet_2015) = c("FNR12 - FNLS", "FNR36 - FNLS", "FNR60 - FNLS", "FDR12 - FDLS", "FDR36 - FDLS", "FDR60 - FDLS", "FAR12 - FALS", "FAR36 - FALS", "FAR60 - FALS")


#check first if rownames are the same before combining data
rownames(log2fc_stress_recovery_flowering_spikelet_2013) == rownames(log2fc_stress_recovery_flowering_spikelet_2014)
rownames(log2fc_stress_recovery_flowering_spikelet_2013) == rownames(log2fc_stress_recovery_flowering_spikelet_2015)

#log2fc - 60h RW / LS
log2fc_flowering_spikelet_R60_latestress_20131415 = cbind(log2fc_stress_recovery_flowering_spikelet_2013[ , c(9, 6, 3)], 
                                                          log2fc_stress_recovery_flowering_spikelet_2014[ , c(9, 6, 3)], 
                                                          log2fc_stress_recovery_flowering_spikelet_2015[ , c(9, 6, 3)])    #log2fc of 60h RW / late stress - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(log2fc_flowering_spikelet_R60_latestress_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", "2014 Anjali", "2014 Dular", "2014 N22", "2015 Anjali", "2015 Dular", "2015 N22")
log2fc_flowering_spikelet_R60_latestress_20131415 = log2fc_flowering_spikelet_R60_latestress_20131415[ , match(colnames(yield_chalk_change_HxD_FL), colnames(log2fc_flowering_spikelet_R60_latestress_20131415))]     #reorder columns to match column order of yield and chalkiness data



#yield, chalkiness, and log2FC in one data table

#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_FL) == colnames(log2fc_flowering_spikelet_R60_latestress_20131415)

yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD = rbind(yield_chalk_change_HxD_FL, log2fc_flowering_spikelet_R60_latestress_20131415)
yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD = t(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD)
yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD = as.data.frame(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD)

#test for normality of data
shapiro_yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD = func_shapiro_test(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD)
shapiro_yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$Distribution = ifelse(shapiro_yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$Distribution == "Normal")   #86 variables
sum(shapiro_yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$Distribution == "Non-normal")   #5 variables


#correlation test - spearman

#correlation coefficient

#accounts for change in yield/chalkiness
cor_res_spearman_rho_R60_LS_flowering_spikelet_HxD = cor(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD, method = "spearman")


#correlation with change in yield and chalkiness

#yield
yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD[4:ncol(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD)], function(x) cor.test(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD = as.data.frame(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD)
colnames(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD) = "p.value"
rownames(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD$Significance = with(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD,
                                                                              ifelse(p.value <= 0.001, "***", ifelse(
                                                                                p.value <= 0.01, "**", ifelse(
                                                                                  p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD$p.value < 0.05))   #7 metabolites
yield_cor_sig_metabolites_R60_LS_flowering_spikelet_HxD = rownames(yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD)[yield_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = data.frame(cor_res_spearman_rho_R60_LS_flowering_spikelet_HxD[1, which(colnames(cor_res_spearman_rho_R60_LS_flowering_spikelet_HxD) %in% yield_cor_sig_metabolites_R60_LS_flowering_spikelet_HxD)])     #metabolites with significant correlation with change in yield and rho values
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = round(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD, 2)   #round to 2 decimals
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = as.matrix(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD)  #convert to matrix
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = cbind(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD, yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD[order(rownames(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD)), ]   #sort metabolites alphabetially
yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD = as.data.frame(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD[ , -2])
colnames(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD) = "rho value"
#export table
write.table(yield_sig_cor_rho_spearman_R60_LS_flowering_spikelet_HxD, file = "cor_yield_sig_R60_LS_flowering_spikelet_HxD.txt", sep = "\t", quote = F)


#chalky grains
chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD[4:ncol(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD)], function(x) cor.test(yield_chalk_log2fc_R60_LS_flowering_spikelet_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD) = reduced_metabolite_list_final_flowering_spikelet_2013$Name
chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD,
                                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                                      p.value <= 0.01, "**", ifelse(
                                                                                        p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_R60_LS_flowering_spikelet_HxD$p.value < 0.05))   #0 metabolite
