setwd("~/Work/Flavi CBD screen/CBD_synergy_2022")
setwd("~/Work/Flavi_screen")

library(tidyverse)
library(readxl)
library(pheatmap)
library(plotly)
 # see the original version of this script for sources of some lines and for the playground with the visuals
# this version gets all four plates
fadu_data = read_excel("plate_map_1_1.xlsx",
                       na = "NA")
fadu_drug_map = read_excel("plate_map_1_drugs.xlsx",
                       na = "NA")


# creating annotation template
#####

fadu_data_ann = bind_cols(rep(1:16,each=24),rep(1:24,16),fadu_data)
colnames(fadu_data_ann) = c("row","col",colnames(fadu_data))
fadu_data_ann = fadu_data_ann %>%  left_join(fadu_drug_map)   %>%  mutate(drugs = case_when(drugs == "DMSO + dmso" ~ "bckg",
                                                                                        drugs == "STAUR + dmso" ~ "staur",
                                                                                        drugs == "dmso wellmate only/unpinned" ~ "dmso",
                                                                                        TRUE ~ drugs))
fadu_data_ann = fadu_data_ann  %>%  
          mutate(group_code = case_when(groups == "DMSO + dmso" ~ 0,
                                    groups == "STAUR + dmso" ~ 10,
                                    groups == "dmso wellmate only/unpinned" ~ 1,
                                    TRUE ~ 2)) %>% 
          mutate(groups = case_when(groups == "DMSO + dmso" ~ "bckg",
                                         groups == "STAUR + dmso" ~ "staur",
                                         groups == "dmso wellmate only/unpinned" ~ "dmso",
                                         TRUE ~ "lib"))

fadu_data_ann %>% write_csv("fadu_data_annotated_plate.csv")

fadu_data_ann = read_csv("fadu_data_annotated_plate.csv")

# creating a non-redundant list of drugs (NOT includind staurosporin, which is a positive control)
drug_list = fadu_data_ann %>% filter(groups == "lib") %>% pull(drugs) %>% unique()
# length(drug_list)
# a total of 75 drugs



# matching DMSO wells to drugs and visuzlizing the hits
#####
#  create pattern for DMSO wells, juxtaposed with averages/SD for the drugs 
# then plotting in 3D

fadu_data_ann_dmso_wells_matching = fadu_data_ann %>% 
          filter(groups == "lib") %>% select(row,col, well,drugs) %>% arrange(drugs, row, col) %>% 
          mutate(col = col + 1) %>%  filter(row_number() %% 3 == 0) %>% select(-well)
# these are coordinate of the wells,
# which are the  in the rightmost bottom corner of 2x2 square of the corresponding drug

fadu_data_ann_dmso_matching = fadu_data_ann %>% select(-drugs) %>% 
          filter(group_code < 2) %>% left_join(fadu_data_ann_dmso_wells_matching) 
fadu_data_ann_dmso_matching = fadu_data_ann_dmso_matching %>%  left_join(fadu_data_ann_calc)



# creating template for all wells with negative controls, with the avg/SD of the drugs in the corresponding positions

fadu_data_controls_left = fadu_data_ann %>% filter(groups != "lib", col<3) %>% 
  mutate(rand = sample(32)) %>% arrange(drugs,rand) %>% 
  mutate(comp_ctrl = paste0(groups, rep(1:5,c(3,3,3,3,4))))
fadu_data_controls_right = fadu_data_ann %>% filter(groups != "lib", col>22) %>% 
  mutate(rand = sample(32)) %>% arrange(drugs,rand) %>% 
  mutate(comp_ctrl = paste0(groups, rep(6:10,c(3,3,3,3,4))))
fadu_data_controls = bind_rows(fadu_data_controls_left,fadu_data_controls_right)
# we will only use randomly assigned status for each well
rm(fadu_data_controls_left)
rm(fadu_data_controls_right)
fadu_data_controls %>% write_csv("fadu_data_controls.csv")

# now creating 20 sets  of 4  randomly selected dmso wells 
# (with replacement, from the block of 20 wells in cols 3 and 4)
dmso_contr_groups = fadu_data_ann %>% filter(groups == "dmso", col==3 | col == 4, row < 11) %>% 
  select(well) %>% slice_sample(n=80,replace=T)  %>%  
  bind_cols(drugs =c(rep(paste0("bckg",1:10),4), rep(paste0("staur",1:10),4)))


fadu_data_controls_summary = tibble(plate = rep(1:4,each=20), drugs = rep(c(paste0("bckg",1:10),paste0("staur",1:10)),4), 
                                      drug_avg = NA, drug_sd = NA, ctrl_count = NA, ctrl_avg = NA, ctrl_sd = NA,
                                      ttest_pval = NA, wilcox_pval = NA)


for (i in 1:nrow(fadu_data_controls_summary)) {
  ctrl = fadu_data_all_ann %>% 
    filter(plate == fadu_data_controls_summary$plate[i], 
           well %in% (fadu_data_controls %>% filter(comp_ctrl == fadu_data_controls_summary$drugs[i]) %>% pull(well))) %>% pull(signal)
  dmso = fadu_data_all_ann %>% 
    filter(plate == fadu_data_controls_summary$plate[i], 
           well %in% (dmso_contr_groups %>% filter(drugs == fadu_data_controls_summary$drugs[i]) %>% pull(well) %>% unique())) %>% pull(signal)
  fadu_data_controls_summary$drug_avg[i] = mean(ctrl)
  fadu_data_controls_summary$drug_sd[i] = sd(ctrl)
  fadu_data_controls_summary$ctrl_count[i] = length(dmso)
  fadu_data_controls_summary$ctrl_avg[i] = mean(dmso)
  fadu_data_controls_summary$ctrl_sd[i] = sd(dmso)
  fadu_data_controls_summary$ttest_pval[i] = t.test(ctrl,dmso)$p.value
  fadu_data_controls_summary$wilcox_pval[i] = wilcox.test(ctrl,dmso)$p.value
    
}

fadu_data_controls_summary = fadu_data_controls_summary %>%
  mutate(ttest_color = case_when(ttest_pval <= 0.05 ~ "orange",
                                 ttest_pval <= 0.005 ~ "red",
                                 TRUE ~ "gray") ,
         wilcox_color = case_when(wilcox_pval <= 0.05 ~ "orange",
                                  wilcox_pval <= 0.005 ~ "red",
                                  TRUE ~ "gray"),
         drug_to_ctlr = drug_avg/ctrl_avg, d_t_c_SD = drug_to_ctlr*sqrt((ctrl_sd/ctrl_avg)^2 + (drug_sd/drug_avg)^2 ))  
  


pl <- plot_ly(
          fadu_data_ann_dmso_matching, x= ~row, y= ~col, z= ~signal,
          type='mesh3d', intensity = ~signal,
          colors=  colorRamp(gray.colors(5))
)
# surface plot for negative controls
pl2 = pl %>% layout(scene = list(zaxis=list(
          range = c(0,650000))))
pl2
#####

# creating template for averaging neg controls surrounding the drugs
#####

# 
drug_coord = fadu_data_ann %>% 
  filter(groups == "lib") %>% select(row,col, well,drugs) %>% arrange(drugs, row, col)
dmso_wells_around_drugs = tibble(drugs = drug_list, row_max = NA,row_min = NA, col_max = NA,col_min = NA)

# extracting the quadrant, where each drug is, and expanding it by 1 well in each direction
for (i in 1:length(drug_list)) {
  dmso_wells_around_drugs$row_max[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(row) %>% max()+1
  dmso_wells_around_drugs$row_min[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(row) %>% min()-1
  dmso_wells_around_drugs$col_max[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(col) %>% max()+1
  dmso_wells_around_drugs$col_min[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(col) %>% min()-1
  
}


dmso_wells_around_drugs %>% write_csv("dmso_wells_around_drugs.csv")
dmso_wells_around_drugs = read_csv("dmso_wells_around_drugs.csv")



# importing all data
#####
# reading all image-based data (note the Excel file was manually purged of all comments!!!)
fadu_data_all = read_excel("Golemis_FS_08112023_reanalysis_screen2_NC_SUM_edIS.xlsx") %>% select(1,2,6)
colnames(fadu_data_all) = c("plate","well","signal")
fadu_data_all = fadu_data_all %>% mutate(plate = rep(1:5, each = 384))
tail(fadu_data_all)
# visuzlizing unperturbed plate
#####
unperturbed = bind_cols(rep(1:16,each=24),rep(1:24,16),fadu_data_all %>% filter(plate==5) ) 
colnames(unperturbed) = c("row","col","plate","well","signal")
matrix(unperturbed$signal,16,24, byrow = T) %>% pheatmap(cluster_rows = F,cluster_cols = F)
# image saved as "unperturbed plate run 2 total area heatmap.png"

summary(unperturbed$signal)
sd(unperturbed$signal)

# plu <- plot_ly(
#   unperturbed, x= ~row, y= ~col, z= ~signal,
#   type='mesh3d', intensity = ~signal,
#   colors=  colorRamp(gray.colors(5))
# )
# plu2 = plu %>% layout(scene = list(zaxis=list(
#   range = c(0,600000))))
# plu2
# # image saved as "uperturbed view1.png", "uperturbed view2.png"


# annotating all data, calculating average and SD, plotting comparison (violin plots)

# creating rough summary and violin plots
#####
fadu_data_all_ann = fadu_data_all %>% filter(plate<5) %>% left_join(fadu_data_ann %>% select(-signal))
# rough summary by library/controls
fadu_data_all_ann_summary  = fadu_data_all_ann %>% group_by(plate,groups) %>% 
  summarise(Avg = mean(signal, na.rm=T), SD = sd(signal, na.rm=T))
fadu_data_all_ann_summary %>% pivot_wider(names_from = plate, values_from = c(Avg,SD)) %>% 
  write_csv("384plate run 2 total area summary.csv")

fadu_data_all_ann %>% mutate(annotation = paste(plate, groups)) %>%  
  ggplot(aes(x = as.factor(annotation),y = signal)) +
  geom_violin(trim = F)+ 
  geom_boxplot(width=0.1, color="red") +
  theme(axis.text.x = element_text(angle = 45))
ggsave("384plate run 2 total area summary.png")


#####

# calculating average and SD by drugs
fadu_data_all_drugs_vs_ctrl = tibble(plate = rep(1:4,each=75), drugs = rep(drug_list,4), 
                                     drug_avg = NA, drug_sd = NA, ctrl_count = NA, ctrl_avg = NA, ctrl_sd = NA,
                                     ttest_pval = NA, wilcox_pval = NA)



## temporarily joining col and row values to define control matrix
fadu_data_all_drugs_vs_ctrl =fadu_data_all_drugs_vs_ctrl %>% left_join(dmso_wells_around_drugs) 

# calculating averages, standard deviations, 
# and p-values for comparison between drugs and corresponding controls

for (i in 1:nrow(fadu_data_all_drugs_vs_ctrl)) {
  control_matrix = fadu_data_all_ann %>% 
    filter(plate == fadu_data_all_drugs_vs_ctrl$plate[i],row <= fadu_data_all_drugs_vs_ctrl$row_max[i] & row >= fadu_data_all_drugs_vs_ctrl$row_min[i] ) %>% 
    filter(plate == fadu_data_all_drugs_vs_ctrl$plate[i], col <= fadu_data_all_drugs_vs_ctrl$col_max[i] & col >= fadu_data_all_drugs_vs_ctrl$col_min[i] ) %>% 
    filter(drugs == "dmso")
  # print(control_matrix)
  control  = control_matrix %>% pull(signal)
  # print(mean(control))
  fadu_data_all_drugs_vs_ctrl$ctrl_count[i] =control_matrix %>% nrow()
  fadu_data_all_drugs_vs_ctrl$ctrl_avg[i] = mean(control)
  fadu_data_all_drugs_vs_ctrl$ctrl_sd[i] = sd(control)
  
  drug = fadu_data_all_ann %>% filter(plate == fadu_data_all_drugs_vs_ctrl$plate[i],drugs ==  fadu_data_all_drugs_vs_ctrl$drugs[i]) %>% 
    pull(signal)
  fadu_data_all_drugs_vs_ctrl$drug_avg[i] = mean(drug)
  fadu_data_all_drugs_vs_ctrl$drug_sd[i] = sd(drug)
  
  fadu_data_all_drugs_vs_ctrl$ttest_pval[i] =t.test(drug,control)$p.value
  fadu_data_all_drugs_vs_ctrl$wilcox_pval[i] =wilcox.test(drug,control)$p.value
  
  
}



fadu_data_all_drugs_vs_ctrl = fadu_data_all_drugs_vs_ctrl %>%
  mutate(ttest_color = case_when(ttest_pval <= 0.05 ~ "orange",
                                ttest_pval <= 0.005 ~ "red",
                               TRUE ~ "gray") ,
         wilcox_color = case_when(wilcox_pval <= 0.05 ~ "orange",
                                  wilcox_pval <= 0.005 ~ "red",
                                 TRUE ~ "gray"),
         drug_to_ctlr = drug_avg/ctrl_avg, d_t_c_SD = drug_to_ctlr*sqrt((ctrl_sd/ctrl_avg)^2 + (drug_sd/drug_avg)^2 )) %>% 
  select(-row_max, -row_min,-col_max, -col_min)


fadu_data_all_drugs_vs_ctrl %>% write_csv("run2 fadu_data_all_drugs_vs_ctrl total area.csv")



fadu_data_all_drugs_vs_ctrl_plot = fadu_data_all_drugs_vs_ctrl %>% bind_rows(fadu_data_controls_summary) %>% 
  select(plate, drugs, drug_avg,drug_sd,ctrl_avg, ctrl_sd, drug_to_ctlr, d_t_c_SD,ttest_pval,wilcox_pval) %>% 
  pivot_wider(names_from = plate, values_from = c(drug_avg,drug_sd,ctrl_avg, ctrl_sd, drug_to_ctlr, d_t_c_SD,ttest_pval,wilcox_pval))

fadu_data_all_drugs_vs_ctrl_plot %>% write_csv("run2 fadu_data_all_drugs_vs_ctrl_plot  total area.csv")

#  compare run 1. vs run 2
run_comparison = fadu_data_all_drugs_vs_ctrl_plot %>% select(drugs, starts_with("drug_to"), starts_with("d_t_c")) 
colnames(run_comparison)
colnames(run_comparison)[2:9] = paste0(colnames(run_comparison)[2:9], "run1")
qq = read_csv("run2 fadu_data_all_drugs_vs_ctrl_plot.csv") %>% select(drugs, starts_with("drug_to"), starts_with("d_t_c"))
colnames(qq)[2:9] = paste0(colnames(qq)[2:9], "run2")
run_comparison = run_comparison %>% left_join(qq)

run_comparison   %>% 
  ggplot(aes(x = drug_to_ctlr_1run1, y = drug_to_ctlr_1run2))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  # geom_errorbar(aes(ymin = drug_to_ctlr_1run2 - 2*d_t_c_SD_1run1, ymax = drug_to_ctlr_1run2 + 2*d_t_c_SD_1run2), width = 0.1, alpha = 0.3) +
  # geom_errorbarh(aes(xmin = drug_to_ctlr_1run1 - 2*d_t_c_SD_1run1, xmax = drug_to_ctlr_1run1 + 2*d_t_c_SD_1run1), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(0,1.5),ylim = c(0,1.5))
ggsave("plate1 drugs run1 vs run2.png")
run_comparison %>% select(drug_to_ctlr_1run1, drug_to_ctlr_1run2) %>% cor()
run_comparison %>% select(-drugs, -starts_with("d_t_c")) %>% cor() %>% as.data.frame() %>% write.csv("correlation between runs.csv")

fadu_data_all_no_drugs = fadu_data_all_ann  %>%   filter(group_code !=2) %>% 
  pivot_wider(names_from = plate, values_from = c(signal)) %>% 
  select(-c(well:col, group_code,  drugs))

colnames(fadu_data_all_no_drugs)  = c("groups",paste0("drug_avg_",1:4))

fadu_data_all_for_correl = fadu_data_all_drugs_vs_ctrl %>% 
  select(plate,drugs,drug_avg,drug_sd) %>% pivot_wider(names_from = plate, values_from = c(drug_avg,drug_sd)) %>% 
  mutate(groups = "lib") %>% select(-drugs) %>% relocate(groups,) %>% 
  bind_rows(fadu_data_all_no_drugs)


fadu_data_all_for_correl   %>% 
  ggplot(aes(x = drug_avg_1, y = drug_avg_2, color = groups))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  geom_errorbar(aes(ymin = drug_avg_2 - 2*drug_sd_2, ymax = drug_avg_2 + 2*drug_sd_2), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(xmin = drug_avg_1 - 2*drug_sd_1, xmax = drug_avg_1 + 2*drug_sd_1), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(0,500000),ylim = c(0,500000))
# note errorbars are now close to 95% confidence intervals
ggsave("run1 plate1 vs plate2.png", width = 2400, height = 2100, units = "px")


fadu_data_all_for_correl   %>% 
  ggplot(aes(x = drug_avg_1, y = drug_avg_3, color = groups))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  geom_errorbar(aes(ymin = drug_avg_3 - 2*drug_sd_3, ymax = drug_avg_3 + 2*drug_sd_3), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(xmin = drug_avg_1 - 2*drug_sd_1, xmax = drug_avg_1 + 2*drug_sd_1), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(0,500000),ylim = c(0,500000))
# note errorbars are now close to 95% confidence intervals
ggsave("run1 plate1 vs plate3.png", width = 2400, height = 2100, units = "px")


fadu_data_all_for_correl   %>% 
  ggplot(aes(x = drug_avg_4, y = drug_avg_3, color = groups))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  geom_errorbar(aes(ymin = drug_avg_3 - 2*drug_sd_3, ymax = drug_avg_3 + 2*drug_sd_3), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(xmin = drug_avg_4 - 2*drug_sd_4, xmax = drug_avg_4 + 2*drug_sd_4), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(0,500000),ylim = c(0,500000))
# note errorbars are now close to 95% confidence intervals
ggsave("run1 plate3 vs plate4.png", width = 2400, height = 2100, units = "px")


pl <- plot_ly(
  (fadu_data_all_ann %>% filter(plate == 1, group_code <2)), x= ~row, y= ~col, z= ~signal,
  type='mesh3d', intensity = ~signal,
  colors=  colorRamp(gray.colors(5))
)
# surface plot for negative controls + drugs colored by significance
pl1 = pl %>% layout(scene = list(zaxis=list(
  range = c(0,450000))))
pl1

plate1_hits = fadu_data_ann_dmso_wells_matching %>% 
  left_join( fadu_data_all_drugs_vs_ctrl %>% filter(plate ==1) %>%  select(drugs, drug_avg, ttest_color))
add_markers(pl1, x = plate1_hits$row, y = plate1_hits$col, z = plate1_hits$drug_avg, 
            marker = list(color = ~plate1_hits$ttest_color))



fadu_data_all_for_plot = fadu_data_all_for_plot %>% 
  mutate(ratio_2_1 = Avg_2/Avg_1, ratio_3_1 = Avg_3/Avg_1)
fadu_data_all_for_plot %>% filter(ratio_2_1 < 0.75) %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_3_1 < 0.75) %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_2_1 > 1.25) %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_3_1 > 1.25) %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_2_1 < 0.75, groups == "lib") %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_3_1 < 0.75, groups == "lib") %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_2_1 > 1.25, groups == "lib") %>% nrow()
fadu_data_all_for_plot %>% filter(ratio_3_1 > 1.25, groups == "lib") %>% nrow()

fadu_data_all_for_plot %>% write_csv("fadu_data_all_for_plot.csv")
fadu_data_all_ann_calc %>% left_join()



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
