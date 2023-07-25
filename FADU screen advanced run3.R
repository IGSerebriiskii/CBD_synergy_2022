setwd("~/Work/CBD_synergy_2022")
library(tidyverse)
library(readxl)
library(pheatmap)

fadu_data = read_excel("CBD screening _3 hours and 30 min reading_R3_121222_column format.xlsx",
                       na = "NA")

fadu_data %>% nrow() #1824 = 19*96


cbd_plates = c(10:18)
ctrl_plates = c(1:9)
fadu_data  = fadu_data %>% mutate(cbd = ifelse(Plate %in% cbd_plates, "cbd","control"))
# fadu_data  = fadu_data %>% mutate(comp = ifelse(Plate %in% c(1:3,7:9), 2,200))
fadu_data  = fadu_data %>% filter(Plate<19) %>%  mutate(comp = case_when(Plate %in% c(1:3,10:12) ~ 200,
                                                   Plate %in% c(4:6,13:15) ~ 2000,
                                                   TRUE ~ 5000))


# Defining control wells
wells = fadu_data$Well %>% unique()
fadu_dmso_wells1 = c(paste0(LETTERS[1:4],"01"))
fadu_dmso_wells2 = c(paste0(LETTERS[5:8],12))
fadu_staur_wells1 = c(paste0(LETTERS[5:8],"01"))
fadu_staur_wells2 = c(paste0(LETTERS[1:4],12))
fadu_dmso_wells = c(fadu_dmso_wells1,fadu_dmso_wells2)
fadu_staur_wells = c(fadu_staur_wells1,fadu_staur_wells2)
fadu_cbd_wells = c(paste0(LETTERS[1:5],"02"))
# - - - - - - - - - - - - - - - - - - -

fadu_data_annotated  = fadu_data %>% mutate( well_type = case_when(Well %in% fadu_dmso_wells ~ "dmso", 
                                                        Well %in% fadu_staur_wells ~ "PC",
                                                        Well %in% fadu_cbd_wells ~ "cbd",
                                                        TRUE ~ "test_compounds"))

fadu_data_summary = fadu_data_annotated%>% group_by(well_type,cbd,comp) %>% 
  summarise(Avg = mean(Signal, na.rm=T), SD = sd(Signal, na.rm=T)) %>% arrange(-comp, desc(cbd),well_type)
fadu_data_summary %>% write_csv("fadu_data_summary.csv")
fadu_data_annotated = fadu_data_annotated %>% mutate(annot = paste(cbd,comp,well_type, sep  = "_"))
fadu_data_annotated %>% ggplot(aes(x = as.factor(annot),y = Signal)) +
  geom_violin(trim = F)+ 
  geom_boxplot(width=0.1, color="red") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("summary_violin run3.png",
       width = 30, #note it was 20 for 2 concentrations
       height = 10,
       units =  "cm",
       dpi = 300)

# did not calculate
# comp_matrix = as.data.frame(matrix(0, nrow = length(sort(unique(fadu_data_annotated$annot))), ncol  = length(sort(unique(fadu_data_annotated$annot)))))
# # colnames(my_results) = rownames(my_results) = paste0("h-",unique(PTEN_variants$htsp_MTL))
# colnames(comp_matrix) = rownames(comp_matrix) = sort(unique(fadu_data_annotated$annot))
# 
# comp_matrix_ks_test = comp_matrix_ttest = comp_matrix
# 
# for (i in 1:nrow(comp_matrix)) {
#   for (j in 1:nrow(comp_matrix)) {
#     comp_matrix_ttest[i,j] = t.test(filter(fadu_data_annotated, annot == colnames(comp_matrix)[i])$Signal, filter(fadu_data_annotated, annot == colnames(comp_matrix)[j])$Signal)$p.value
#     comp_matrix_ks_test[i,j] = ks.test(filter(fadu_data_annotated, annot == colnames(comp_matrix)[i])$Signal, filter(fadu_data_annotated, annot == colnames(comp_matrix)[j])$Signal)$p.value
#     
#   }
# }

# diag(comp_matrix) = NA
# comp_matrix[upper.tri(comp_matrix)] <- comp_matrix_ks_test[upper.tri(comp_matrix)]
# comp_matrix[lower.tri(comp_matrix)] <- comp_matrix_ttest[lower.tri(comp_matrix)]
# comp_matrix %>% write.csv("comp_matrix run2.csv")

# - - - - - - - - - - - - - - - - - - -

fadu_data_mean = fadu_data %>% group_by(Well,cbd,comp) %>% 
  summarise(Avg = mean(Signal, na.rm=T), SD = sd(Signal, na.rm=T))

# - - - - - - - - - - - - - - - - - - -
fadu_2000 = matrix(filter(fadu_data_mean,comp == 2000, cbd == "control")$Avg,8,12, byrow = T)
pheatmap(fadu_2000,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_2000,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 control comp_2000nm avg.pdf")
fadu_200 = matrix(filter(fadu_data_mean,comp == 200, cbd == "control")$Avg,8,12, byrow = T)
pheatmap(fadu_200,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_200,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 control comp_200nm avg.pdf")

fadu_5000 = matrix(filter(fadu_data_mean,comp == 5000, cbd == "control")$Avg,8,12, byrow = T)
pheatmap(fadu_5000,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 control comp_5000nm avg.pdf")

fadu_cbd2000 = matrix(filter(fadu_data_mean,comp == 2000, cbd == "cbd")$Avg,8,12, byrow = T)
pheatmap(fadu_cbd2000,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_cbd2000,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd comp_2000nm avg.pdf")

fadu_cbd200 = matrix(filter(fadu_data_mean,comp == 200, cbd == "cbd")$Avg,8,12, byrow = T)
pheatmap(fadu_cbd200,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_cbd200,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd comp_200nm avg.pdf")

fadu_cbd5000 = matrix(filter(fadu_data_mean,comp == 5000, cbd == "cbd")$Avg,8,12, byrow = T)
pheatmap(fadu_cbd5000,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_cbd5000,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd comp_5000nm avg.pdf")

 # note the "breaks..." part comes from the inspection of the unconstrained results

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cbd to control ratios
fadu_data_by_cbd_plot = fadu_data_mean %>% 
  pivot_wider(names_from = cbd, values_from = c(Avg,SD)) %>% 
  mutate(ratio = Avg_cbd/Avg_control, log2ratio = log2(ratio), 
         ratioSD = ratio*sqrt((SD_control/Avg_control)^2 + (SD_cbd/Avg_cbd)^2))


fadu_data_by_cbd_plot = fadu_data_by_cbd_plot %>% 
  mutate( well_type = case_when(Well %in% fadu_dmso_wells ~ "dmso", 
                                Well %in% fadu_staur_wells ~ "PC",
                                Well %in% fadu_cbd_wells ~ "cbd",
                                TRUE ~ "test_compounds"))

fadu_data_by_cbd_plot %>% filter(comp == 2000)  %>% 
  ggplot(aes(x = Avg_control, y = Avg_cbd, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_cbd - SD_cbd, ymax = Avg_cbd + SD_cbd), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control - SD_control, xmax = Avg_control + SD_control), height = 0.1) +
  coord_cartesian(xlim = c(0,250000),ylim = c(0,250000))

ggsave("control vs cbd run3 conc2000 detailed with errorbars.png")

fadu_data_by_cbd_plot %>% filter(comp == 200)  %>% 
  ggplot(aes(x = Avg_control, y = Avg_cbd, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_cbd - SD_cbd, ymax = Avg_cbd + SD_cbd), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control - SD_control, xmax = Avg_control + SD_control), height = 0.1) +
  coord_cartesian(xlim = c(0,250000),ylim = c(0,250000))

ggsave("control vs cbd run3 conc200 detailed with errorbars.png")

fadu_data_by_cbd_plot %>% filter(comp == 5000)  %>% 
  ggplot(aes(x = Avg_control, y = Avg_cbd, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_cbd - SD_cbd, ymax = Avg_cbd + SD_cbd), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control - SD_control, xmax = Avg_control + SD_control), height = 0.1) +
  coord_cartesian(xlim = c(0,250000),ylim = c(0,250000))

ggsave("control vs cbd run3 conc5000 detailed with errorbars.png")

fadu_data_by_cbd = fadu_data_by_cbd_plot

fadu_data_by_cbd200 = matrix(filter(fadu_data_by_cbd,comp == 200)$ratio,8,12, byrow = T) 
# pheatmap(fadu_data_by_cbd2,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F, breaks = seq(0.9, 1.5, length.out = 100), 
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd_to_control comp 200nm avg.pdf")
fadu_data_by_cbd2000 = matrix(filter(fadu_data_by_cbd,comp == 2000)$ratio,8,12, byrow = T) 
pheatmap(fadu_data_by_cbd2000,cluster_rows = F,cluster_cols = F, breaks = seq(0.9, 1.5, length.out = 100), 
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd_to_control log2ratio 2000nm avg.pdf")

fadu_data_by_cbd5000 = matrix(filter(fadu_data_by_cbd,comp == 5000)$ratio,8,12, byrow = T) 
# pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_cbd5000,cluster_rows = F,cluster_cols = F, breaks = seq(0.9, 1.5, length.out = 100), 
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd_to_control comp 5000nm avg.pdf")

# fadu_data_by_cbd200 = matrix(filter(fadu_data_by_cbd,comp == 200)$log2ratio,8,12, byrow = T) 
# pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F, breaks = seq(-1,1, length.out = 100), 
#          cellwidth = 15, cellheight = 15, fontsize = 8,  
#          filename = "FADU run2 cbd_to_control log2ratio 200nm avg.pdf")


fadu_data_by_comp = fadu_data_mean %>% 
  pivot_wider(names_from = comp, values_from = c(Avg,SD))

fadu_data_by_comp = fadu_data_by_comp %>% mutate(ratio2 = Avg_2000/Avg_200, 
                                                 log2ratio2 = log2(ratio2), 
                                                 ratio5 = Avg_5000/Avg_200, 
                                                 log2ratio5 = log2(ratio5))

fadu_data_by_comp_control = matrix(filter(fadu_data_by_comp,cbd == "control")$ratio2,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
# pheatmap(fadu_data_by_comp_control, cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_comp_control, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
         filename = "FADU run3 no_cbd 2000_to_low comp avg.pdf")

fadu_data_by_comp_control = matrix(filter(fadu_data_by_comp,cbd == "control")$log2ratio2,8,12, byrow = T) 
pheatmap(fadu_data_by_comp_control, cluster_rows = F, cluster_cols = F, cellwidth = 15, 
         cellheight = 15, breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run3 no_cbd 2000_to_low log2ratio avg.pdf")

fadu_data_by_comp_control = matrix(filter(fadu_data_by_comp,cbd == "control")$ratio5,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
pheatmap(fadu_data_by_comp_control, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
         filename = "FADU run3 no_cbd 5000_to_low comp avg.pdf")

matrix(filter(fadu_data_by_comp,cbd == "control")$log2ratio5,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, 
         cellheight = 15, breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run3 no_cbd 5000_to_low log2ratio avg.pdf")


matrix(filter(fadu_data_by_comp,cbd == "cbd")$log2ratio5,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, 
         cellwidth = 15, cellheight = 15,  breaks = seq(-1,1, length.out = 100),
         filename = "FADU run3 cbd 5000_to_low comp log2ratio avg.pdf")

matrix(filter(fadu_data_by_comp,cbd == "cbd")$log2ratio2,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, 
            cellwidth = 15, cellheight = 15,  breaks = seq(-1,1, length.out = 100),
            filename = "FADU run3 cbd 2000_to_low comp log2ratio avg.pdf")





norm_sig = c()
for (i in 1:18) {
  norm_sig = c(norm_sig,filter(fadu_data_annotated,Plate == i)$Signal/mean(filter(fadu_data_annotated,Plate == i, well_type == "cbd")$Signal))
}

fadu_data_normalized = fadu_data_annotated %>% filter(Plate<19) %>% 
  mutate(sig_norm = norm_sig)

fadu_data_norm_mean = fadu_data_normalized %>% group_by(Well,cbd,comp) %>% 
  summarise(Avg = mean(sig_norm, na.rm = T), SD = sd(sig_norm, na.rm = T))

# - - - - - - - - - - - - - - - - - - -
matrix(filter(fadu_data_norm_mean,comp == 2000, cbd == "control")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU norm run3 control comp_2000nm avg.pdf")

matrix(filter(fadu_data_norm_mean,comp == 5000, cbd == "control")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU norm run3 control comp_5000nm avg.pdf")

matrix(filter(fadu_data_norm_mean,comp == 200, cbd == "control")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU norm run3 control comp_200nm avg.pdf")

matrix(filter(fadu_data_norm_mean,comp == 2000, cbd == "cbd")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU norm run3 cbd comp_2000nm avg.pdf")

matrix(filter(fadu_data_norm_mean,comp == 5000, cbd == "cbd")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU norm run3 cbd comp_5000nm avg.pdf")

matrix(filter(fadu_data_norm_mean,comp == 200, cbd == "cbd")$Avg,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(0.1,2, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU norm run3 cbd comp_200nm avg.pdf")

fadu_data_norm_by_cbd = fadu_data_norm_mean %>% 
  pivot_wider(names_from = cbd, values_from = c(Avg,SD)) %>% 
  mutate(ratio = Avg_cbd/Avg_control, log2ratio = log2(ratio), 
         ratioSD = ratio*sqrt((SD_control/Avg_control)^2 + (SD_cbd/Avg_cbd)^2))

fadu_data_norm_by_cbd = fadu_data_norm_by_cbd %>% 
  mutate( well_type = case_when(Well %in% fadu_dmso_wells ~ "dmso", 
                                Well %in% fadu_staur_wells ~ "PC",
                                Well %in% fadu_cbd_wells ~ "cbd",
                                TRUE ~ "test_compounds"))

fadu_data_norm_by_cbd = fadu_data_norm_by_cbd %>% 
  mutate( correction = ifelse(ratio > 1, ratio - 2*ratioSD, ratio + 2*ratioSD) ,
          significance = sign(log2ratio)*sign(log2(correction)) + 1,
          sign_values = log2ratio*significance/2)
# after this step, column "sign_values" contains the same values as in log2ratio, 
# but all non-significant values are converted to zeroes

fadu_data_norm_by_cbd %>% 
  arrange(-comp, well_type) %>%  write_csv("fadu_data run3_norm_by_cbd.csv")

matrix(filter(fadu_data_norm_by_cbd,comp == 2000)$sign_values,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(-1,1, length.out = 100),  
         filename = "FADU run3 cbd 2000 comp avg.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 200)$sign_values,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 200 comp avg.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 5000)$sign_values,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 5000 comp avg.pdf")


matrix(filter(fadu_data_norm_by_cbd,comp == 2000)$sign_values,8,12, byrow = T) %>% 
  na_if(0) %>%   pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 2000 sign diff.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 200)$sign_values,8,12, byrow = T) %>%  
  na_if(0) %>%   pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 200 sign diff.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 5000)$sign_values,8,12, byrow = T) %>%  
  na_if(0) %>%   pheatmap( cluster_rows = F, cluster_cols = F, 
                           cellwidth = 15, cellheight = 15, 
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 5000 sign diff.pdf")




fadu_data_by_comp_cbd = matrix(filter(fadu_data_by_comp,cbd == "cbd")$log2ratio,8,12, byrow = T) 
pheatmap(fadu_data_norm_by_cbd_viz_comp_2, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
         breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run2 cbd_to_nocbd comp2000 log2ratio normalized .pdf")
pheatmap(na_if(fadu_data_norm_by_cbd_viz_comp_2,0), cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
         breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run2 cbd_to_nocbd comp2000 log2ratio normalized withNAs.pdf")

fadu_data_norm_by_cbd_viz_comp_200 = matrix(filter(fadu_data_norm_by_cbd,comp == 200)$sign_values,8,12, byrow = T) 
pheatmap(fadu_data_norm_by_cbd_viz_comp_200, cluster_rows = F, cluster_cols = F, 
         cellwidth = 15, cellheight = 15, breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run2 cbd_to_nocbd comp200 log2ratio normalized .pdf")
pheatmap(na_if(fadu_data_norm_by_cbd_viz_comp_200,0), cluster_rows = F, cluster_cols = F, 
         cellwidth = 15, cellheight = 15, breaks = seq(-1,1, length.out = 100), 
         filename = "FADU run2 cbd_to_nocbd comp200 log2ratio normalized withNAs.pdf")

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
matrix(filter(fadu_data_norm_by_cbd,comp == 2000)$Avg_control,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(0.5,1.5, length.out = 100),  
            filename = "FADU run3 comp 2000 normalized avg.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 200)$Avg_control,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(0.5,1.5, length.out = 100),  
            filename = "FADU run3 comp 200 normalized avg.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 5000)$Avg_control,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15,
            breaks = seq(0.5,1.5, length.out = 100),  
            filename = "FADU run3 comp 5000 normalized avg.pdf")
matrix(filter(fadu_data_norm_by_cbd,comp == 5000)$sign_values,8,12, byrow = T) %>%  
  pheatmap( cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, 
            breaks = seq(-1,1, length.out = 100),  
            filename = "FADU run3 cbd 5000 comp avg.pdf")
# - - - - - - - - - - - - - - - - - - - - - - - - - - 
fadu_data_normalized_controls_run3 = fadu_data_norm_by_cbd %>% select(Well,well_type, comp, Avg_control,SD_control) %>% 
  pivot_wider(names_from = comp, values_from = c(Avg_control, SD_control))
fadu_data_normalized_controls_run3 %>% write_csv("fadu_data_normalized_controls_run3.csv")
fadu_data_normalized_controls_run3 %>% ggplot(aes(x = Avg_control_200, y = Avg_control_2000, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_control_2000 - SD_control_2000, ymax = Avg_control_2000 + SD_control_2000), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control_200 - SD_control_200, xmax = Avg_control_200 + SD_control_200), height = 0.1) +
  coord_cartesian(xlim = c(0,1.5),ylim = c(0,1.5))
ggsave("run3 200 vs 2000.png") 

fadu_data_normalized_controls_run3 %>% ggplot(aes(x = Avg_control_200, y = Avg_control_5000, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_control_5000 - SD_control_5000, ymax = Avg_control_5000 + SD_control_5000), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control_200 - SD_control_200, xmax = Avg_control_200 + SD_control_200), height = 0.1) +
  coord_cartesian(xlim = c(0,1.5),ylim = c(0,1.5))
ggsave("run3 200 vs 5000.png") 

fadu_data_normalized_controls_run3 %>% ggplot(aes(x = Avg_control_2000, y = Avg_control_5000, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_control_5000 - SD_control_5000, ymax = Avg_control_5000 + SD_control_5000), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control_2000 - SD_control_2000, xmax = Avg_control_2000 + SD_control_2000), height = 0.1) +
  coord_cartesian(xlim = c(0,1.5),ylim = c(0,1.5))
ggsave("run3 2000 vs 5000.png") 


matrix(filter(fadu_data,Plate == 16)$Signal,8,12, byrow = T) %>% 
pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(30000, 150000, length.out = 100),  
         cellwidth = 15, cellheight = 15, fontsize = 8,  
         filename = "FADU run3 cbd plate16.pdf")

fadu_data_annotated %>% filter(Plate == 16, well_type != "test_compounds") %>% 
  ggplot(aes(x = as.factor(well_type),y = Signal)) + geom_dotplot()
fadu_data_annotated %>% filter(Plate == 16, well_type != "test_compounds") %>% 
  ggplot(aes(x = as.factor(well_type),y = Signal)) + geom_boxplot()
ggsave("run3 plate16 controls boxplot.png")
fadu_data_annotated %>% filter(Plate == 16, well_type != "test_compounds") %>% 
  ggplot(aes(x = well_type,y = Signal,color = well_type)) + geom_point()
ggsave("run3 plate16 controls dotplot.png")

matrix(filter(fadu_data,Plate == 16)$Signal,8,12, byrow = T) %>% 
  pheatmap(cluster_rows = F,cluster_cols = F, breaks = seq(15000, 75000, length.out = 100),  
           cellwidth = 15, cellheight = 15, fontsize = 8,  
           filename = "FADU run3 cbd plate16 nucl area.pdf")
