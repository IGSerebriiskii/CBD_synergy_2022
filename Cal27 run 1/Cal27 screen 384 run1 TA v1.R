setwd("~/Work/Flavi CBD screen/CBD_synergy_2022/Cal27 run 1")
setwd("~/Work/Flavi_screen/Cal27 run 1") # at home

library(tidyverse)
library(readxl)
library(pheatmap)
library(plotly)

# note that the names of variables are still saying "fadu" while the data are really from Cal27
#   !!!

 # see the original version of this script for sources of some lines and for the playground with the visuals
# this version gets all four plates

# reading all data
cal27_data = read_excel("Cal27_TA  run 1.xlsx",
                       na = "NA")
drug_map = read_excel("../plate_maps/plate_map_1_drugs.xlsx",
                           na = "NA")
drug_map_ann = bind_cols(row = rep(1:16,each=24),col = rep(1:24,16),drug_map) 
drug_map_ann %>% write_csv("../plate_maps/drug_map_ann.csv")

# creating annotation template
#####

# cal27_data = cal27_data %>%  left_join(drug_map_ann)   %>%  mutate(drugs = case_when(drugs == "DMSO + dmso" ~ "bckg",
                                                                                        # drugs == "STAUR + dmso" ~ "staur",
                                                                                        # drugs == "dmso wellmate only/unpinned" ~ "dmso",
                                                                                        # TRUE ~ drugs))
cal27_data_ann = cal27_data %>%  left_join(drug_map_ann) %>%  
          mutate(group_code = case_when(drugs == "DMSO + dmso" ~ 0,
                                        drugs == "STAUR + dmso" ~ 10,
                                        drugs == "dmso wellmate only/unpinned" ~ 1,
                                    TRUE ~ 2)) %>% 
          mutate(groups = case_when(drugs == "DMSO + dmso" ~ "bckg",
                                    drugs == "STAUR + dmso" ~ "staur",
                                    drugs == "dmso wellmate only/unpinned" ~ "dmso",
                                         TRUE ~ "lib"))

cal27_data_ann %>% write_csv("cal27_data_ann_plate TA run1.csv")


cal27_data_ann_calc = cal27_data_ann %>% group_by(drugs) %>% 
  summarise(Avg = mean(signal, na.rm=T), SD = sd(signal, na.rm=T))




# matching DMSO wells to drugs and visuzlizing the hits
#####
#  create pattern for DMSO wells, juxtaposed with averages/SD for the drugs 
# then plotting in 3D

cal27_data_ann_dmso_wells_matching = cal27_data_ann %>% filter(plate == 1) %>% 
          filter(groups == "lib") %>% select(row,col, well,drugs) %>% arrange(drugs, row, col) %>% 
          mutate(col = col + 1) %>%  filter(row_number() %% 3 == 0) %>% select(-well)
# these are coordinate of the wells,
# which are the  in the rightmost bottom corner of 2x2 square of the corresponding drug

cal27_data_ann_dmso_matching = cal27_data_ann %>% select(-drugs) %>% 
          filter(group_code < 2) %>% left_join(cal27_data_ann_dmso_wells_matching) 
cal27_data_ann_dmso_matching = cal27_data_ann_dmso_matching %>%  left_join(cal27_data_ann_calc)

#  this is how data controls were created
#####
# creating template for all wells with negative controls, with the avg/SD of the drugs in the corresponding positions

# fadu_data_controls_left = fadu_data_ann %>% filter(groups != "lib", col<3) %>% 
#   mutate(rand = sample(32)) %>% arrange(drugs,rand) %>% 
#   mutate(comp_ctrl = paste0(groups, rep(1:5,c(3,3,3,3,4))))
# fadu_data_controls_right = fadu_data_ann %>% filter(groups != "lib", col>22) %>% 
#   mutate(rand = sample(32)) %>% arrange(drugs,rand) %>% 
#   mutate(comp_ctrl = paste0(groups, rep(6:10,c(3,3,3,3,4))))
# fadu_data_controls = bind_rows(fadu_data_controls_left,fadu_data_controls_right)
# # we will only use randomly assigned status for each well
#####
data_controls = read_csv("../plate_maps/fadu_data_controls.csv",
                           na = "NA")

# now creating 20 sets  of 4  randomly selected dmso wells 
# (with replacement, from the block of 20 wells in cols 3 and 4)
dmso_contr_groups = cal27_data_ann %>% filter(plate == 1, groups == "dmso", col==3 | col == 4, row < 11) %>% 
  select(well) %>% slice_sample(n=80,replace=T)  %>%  
  bind_cols(drugs =c(rep(paste0("bckg",1:10),4), rep(paste0("staur",1:10),4)))


cal27_data_controls_summary = tibble(plate = rep(1:4,each=20), drugs = rep(c(paste0("bckg",1:10),paste0("staur",1:10)),4), 
                                      drug_avg = NA, drug_sd = NA, ctrl_count = NA, ctrl_avg = NA, ctrl_sd = NA,
                                      ttest_pval = NA, wilcox_pval = NA)


for (i in 1:nrow(cal27_data_controls_summary)) {
  ctrl = cal27_data_ann %>% 
    filter(plate == cal27_data_controls_summary$plate[i], 
           well %in% (data_controls %>% filter(comp_ctrl == cal27_data_controls_summary$drugs[i]) %>% pull(well))) %>% pull(signal)
  dmso = cal27_data_ann %>% 
    filter(plate == cal27_data_controls_summary$plate[i], 
           well %in% (dmso_contr_groups %>% filter(drugs == cal27_data_controls_summary$drugs[i]) %>% pull(well) %>% unique())) %>% pull(signal)
  cal27_data_controls_summary$drug_avg[i] = mean(ctrl)
  cal27_data_controls_summary$drug_sd[i] = sd(ctrl)
  cal27_data_controls_summary$ctrl_count[i] = length(dmso)
  cal27_data_controls_summary$ctrl_avg[i] = mean(dmso)
  cal27_data_controls_summary$ctrl_sd[i] = sd(dmso)
  cal27_data_controls_summary$ttest_pval[i] = t.test(ctrl,dmso)$p.value
  cal27_data_controls_summary$wilcox_pval[i] = wilcox.test(ctrl,dmso)$p.value
    
}

cal27_data_controls_summary = cal27_data_controls_summary %>%
  mutate(ttest_color = case_when(ttest_pval <= 0.05 ~ "orange",
                                 ttest_pval <= 0.005 ~ "red",
                                 TRUE ~ "gray") ,
         wilcox_color = case_when(wilcox_pval <= 0.05 ~ "orange",
                                  wilcox_pval <= 0.005 ~ "red",
                                  TRUE ~ "gray"),
         drug_to_ctlr = drug_avg/ctrl_avg, d_t_c_SD = drug_to_ctlr*sqrt((ctrl_sd/ctrl_avg)^2 + (drug_sd/drug_avg)^2 ))  
  
# surface plot
#####
pl <- plot_ly(
          fadu_data_ann_dmso_matching, x= ~row, y= ~col, z= ~signal,
          type='mesh3d', intensity = ~signal,
          colors=  colorRamp(gray.colors(5))
)
# surface plot for negative controls
pl2 = pl %>% layout(scene = list(zaxis=list(
          range = c(0,450000))))
pl2
fadu_data_drugs = fadu_data_ann_dmso_matching %>% 
  filter(!is.na(Avg)) %>% filter(drugs !="dmso"|drugs !="bckg") %>% 
   mutate(dot_color = case_when(abs(signal-Avg) > 2*SD ~ "orange",
                                                 abs(signal-Avg) > 3*SD ~ "red",
                                                 TRUE ~ "gray") )
# this adds an estimate of significance 
# (> 2 SD between drug and dmso in the samer "corner" as drug)

add_markers(pl2, x = fadu_data_drugs$row, y = fadu_data_drugs$col, z = fadu_data_drugs$Avg, 
            marker = list(color = ~fadu_data_drugs$dot_color))
#####

# creating template for averaging neg controls surrounding the drugs
#####

# 
drug_coord = cal27_data_ann %>% 
  filter(plate == 1, groups == "lib") %>% select(row,col, well,drugs) %>% arrange(drugs, row, col)
dmso_wells_around_drugs = tibble(drugs = drug_list, row_max = NA,row_min = NA, col_max = NA,col_min = NA)

# # extracting the quadrant, where each drug is, and expanding it by 1 well in each direction
# for (i in 1:length(drug_list)) {
#   dmso_wells_around_drugs$row_max[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(row) %>% max()+1
#   dmso_wells_around_drugs$row_min[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(row) %>% min()-1
#   dmso_wells_around_drugs$col_max[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(col) %>% max()+1
#   dmso_wells_around_drugs$col_min[i] = drug_coord %>% filter(drugs == drug_list[i]) %>% select(col) %>% min()-1
#   
# }
# 
##### 


# dmso_wells_around_drugs %>% write_csv("dmso_wells_around_drugs.csv")
dmso_wells_around_drugs = read_csv("../plate_maps/dmso_wells_around_drugs.csv")

# color-coding for significance
cal27_data_drugs_vs_ctrl = cal27_data_ann_dmso_wells_matching %>% select(row, col, drugs,ctrl_avg,ctrl_sd,p.val) %>% 
  left_join(cal27_data_ann_calc) %>% rename(drug_avg = Avg, drug_sd = SD) %>% 
  mutate(dot_color = case_when(p.val <= 0.05 ~ "orange",
                               p.val <= 0.005 ~ "red",
                               TRUE ~ "gray") )


# importing all data
#####

# visuzlizing unperturbed plate
#####
unperturbed = bind_cols(rep(1:16,each=24),rep(1:24,16),cal27_data_ann %>% filter(plate==5) ) 
matrix(unperturbed$signal,16,24, byrow = T) %>% pheatmap(cluster_rows = F,cluster_cols = F)
# image saved as "unperturbed plate run 1 TA heatmap.png"

summary(unperturbed$signal)
sd(unperturbed$signal)

plu <- plot_ly(
  unperturbed, x= ~row, y= ~col, z= ~signal,
  type='mesh3d', intensity = ~signal,
  colors=  colorRamp(gray.colors(5))
)
plu2 = plu %>% layout(scene = list(zaxis=list(
  range = c(0,600000))))
plu2
# image saved as "uperturbed view1.png", "uperturbed view2.png"


# annotating all data, calculating average and SD, plotting comparison (violin plots)

# creating rough summary and violin plots
#####
# rough summary by library/controls
cal27_data_all_ann_summary  = cal27_data_ann %>% filter(plate<5) %>% group_by(plate,groups) %>% 
  summarise(Avg = mean(signal, na.rm=T), SD = sd(signal, na.rm=T))
cal27_data_all_ann_summary %>% pivot_wider(names_from = plate, values_from = c(Avg,SD)) %>% 
  write_csv("cal27 run 1 TA summary.csv")

cal27_data_ann %>% mutate(annotation = paste(plate, groups)) %>%  
  ggplot(aes(x = as.factor(annotation),y = signal)) +
  geom_violin(trim = F)+ 
  geom_boxplot(width=0.1, color="red") +
  theme(axis.text.x = element_text(angle = 45))
ggsave("cal27 run 1 TA summary.png")


#####

# calculating average and SD by drugs
cal27_drugs_vs_ctrl = tibble(plate = rep(1:4,each=75), drugs = rep(drug_list,4), 
                                     drug_avg = NA, drug_sd = NA, ctrl_count = NA, ctrl_avg = NA, ctrl_sd = NA,
                                     ttest_pval = NA, wilcox_pval = NA)



## temporarily joining col and row values to define control matrix
cal27_drugs_vs_ctrl =cal27_drugs_vs_ctrl %>% left_join(dmso_wells_around_drugs) 

# calculating averages, standard deviations, 
# and p-values for comparison between drugs and corresponding controls

for (i in 1:nrow(cal27_drugs_vs_ctrl)) {
  control_matrix = cal27_data_ann %>% 
    filter(plate == cal27_drugs_vs_ctrl$plate[i],row <= cal27_drugs_vs_ctrl$row_max[i] & row >= cal27_drugs_vs_ctrl$row_min[i] ) %>% 
    filter(plate == cal27_drugs_vs_ctrl$plate[i], col <= cal27_drugs_vs_ctrl$col_max[i] & col >= cal27_drugs_vs_ctrl$col_min[i] ) %>% 
    filter(groups == "dmso")
  # print(control_matrix)
  control  = control_matrix %>% pull(signal)
  # print(mean(control))
  cal27_drugs_vs_ctrl$ctrl_count[i] =control_matrix %>% nrow()
  cal27_drugs_vs_ctrl$ctrl_avg[i] = mean(control)
  cal27_drugs_vs_ctrl$ctrl_sd[i] = sd(control)
  
  drug = cal27_data_ann %>% filter(plate == cal27_drugs_vs_ctrl$plate[i],
                                      drugs ==  cal27_drugs_vs_ctrl$drugs[i]) %>% 
    pull(signal)
  cal27_drugs_vs_ctrl$drug_avg[i] = mean(drug)
  cal27_drugs_vs_ctrl$drug_sd[i] = sd(drug)
  
  cal27_drugs_vs_ctrl$ttest_pval[i] =t.test(drug,control)$p.value
  cal27_drugs_vs_ctrl$wilcox_pval[i] =wilcox.test(drug,control)$p.value
  
  
}



cal27_drugs_vs_ctrl = cal27_drugs_vs_ctrl %>%
  mutate(ttest_color = case_when(ttest_pval <= 0.05 ~ "orange",
                                ttest_pval <= 0.005 ~ "red",
                               TRUE ~ "gray") ,
         wilcox_color = case_when(wilcox_pval <= 0.05 ~ "orange",
                                  wilcox_pval <= 0.005 ~ "red",
                                 TRUE ~ "gray"),
         drug_to_ctlr = drug_avg/ctrl_avg, d_t_c_SD = drug_to_ctlr*sqrt((ctrl_sd/ctrl_avg)^2 + (drug_sd/drug_avg)^2 )) %>% 
  select(-row_max, -row_min,-col_max, -col_min)


cal27_drugs_vs_ctrl %>% write_csv("cal27_drugs_vs_ctrl TA.csv")



cal27_drugs_vs_ctrl_plot = cal27_drugs_vs_ctrl %>% bind_rows(cal27_data_controls_summary) %>% 
  select(plate, drugs, drug_avg,drug_sd,ctrl_avg, ctrl_sd, drug_to_ctlr, d_t_c_SD,ttest_pval,wilcox_pval) %>% 
  pivot_wider(names_from = plate, values_from = c(drug_avg,drug_sd,ctrl_avg, ctrl_sd, drug_to_ctlr, d_t_c_SD,ttest_pval,wilcox_pval))

cal27_drugs_vs_ctrl_plot %>% write_csv("cal27_drugs_vs_ctrl_plot TA.csv")


fadu_data_all_ann_calc %>% 
  pivot_wider(names_from = plate, values_from = c(Avg,SD)) %>% 
  mutate(groups = "lib") %>% select(-drugs) %>% relocate(groups,)


fadu_data_all_no_drugs = fadu_data_all_ann  %>%   filter(group_code !=2) %>% 
  pivot_wider(names_from = plate, values_from = c(signal)) %>% 
  select(-c(well:col, group_code))

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
  coord_cartesian(xlim = c(100000,600000),ylim = c(100000,600000))
# note errorbars are now close to 95% confidence intervals
ggsave("plate1 vs plate2.png", width = 2400, height = 2100, units = "px")


fadu_data_all_for_correl   %>% 
  ggplot(aes(x = drug_avg_1, y = drug_avg_3, color = groups))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  geom_errorbar(aes(ymin = drug_avg_3 - 2*drug_sd_3, ymax = drug_avg_3 + 2*drug_sd_3), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(xmin = drug_avg_1 - 2*drug_sd_1, xmax = drug_avg_1 + 2*drug_sd_1), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(100000,600000),ylim = c(100000,600000))
# note errorbars are now close to 95% confidence intervals
ggsave("plate1 vs plate3.png", width = 2400, height = 2100, units = "px")


fadu_data_all_for_correl   %>% 
  ggplot(aes(x = drug_avg_4, y = drug_avg_3, color = groups))+
  geom_point(aes(alpha = 0.5)) + geom_smooth(method = 'lm',linetype = "dashed", size = 0.7) +
  geom_errorbar(aes(ymin = drug_avg_3 - 2*drug_sd_3, ymax = drug_avg_3 + 2*drug_sd_3), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(xmin = drug_avg_4 - 2*drug_sd_4, xmax = drug_avg_4 + 2*drug_sd_4), height = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(100000,600000),ylim = c(100000,600000))
# note errorbars are now close to 95% confidence intervals
ggsave("plate1 vs plate4.png", width = 2400, height = 2100, units = "px")


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



