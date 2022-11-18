setwd("~/Work/CBD_synergy_2022")
library(tidyverse)
library(readxl)
library(pheatmap)
# fadu_data = read_csv("CBD screening _3 hours reading_R1_111422_column format.csv")
fadu_data = read_excel("CBD screening _3 hours reading_R1_111422_column format.xlsx")

fadu_data %>% nrow() #1152 = 12*96

fadu_data %>% mutate()

cbd_plates = c(7:12)
ctrl_plates = c(1:6)
fadu_data  = fadu_data %>% mutate(cbd = ifelse(Plate %in% cbd_plates, "cbd","control"))
fadu_data  = fadu_data %>% mutate(comp = ifelse(Plate %in% c(1:3,7:9), 2,200))


# Defining control wells

fadu_dmso_wells1 = c(paste0(LETTERS[1:4],"01"))
fadu_dmso_wells2 = c(paste0(LETTERS[5:8],12))
fadu_staur_wells1 = c(paste0(LETTERS[5:8],"01"))
fadu_staur_wells2 = c(paste0(LETTERS[1:4],12))
fadu_dmso_wells = c(fadu_dmso_wells1,fadu_dmso_wells2)
fadu_staur_wells = c(fadu_staur_wells1,fadu_staur_wells2)


fadu_data_mean = fadu_data %>% group_by(Well,cbd,comp) %>% 
  summarise(Avg = mean(Signal), SD = sd(Signal))

# - - - - - - - - - - - - - - - - - - -
fadu_v1 = matrix(filter(fadu_data_mean,comp == 2, cbd == "control")$Avg,8,12, byrow = T)
pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 control comp_2000nm avg.pdf")
fadu_v2 = matrix(filter(fadu_data_mean,comp == 200, cbd == "control")$Avg,8,12, byrow = T)
pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 control comp_2000nm avg.pdf")

fadu_cbd1 = matrix(filter(fadu_data_mean,comp == 2, cbd == "cbd")$Avg,8,12, byrow = T)
pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd comp_2000nm avg.pdf")

fadu_cbd2 = matrix(filter(fadu_data_mean,comp == 200, cbd == "cbd")$Avg,8,12, byrow = T)
pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F, breaks = seq(50000, 250000, length.out = 100),  cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd comp_200nm avg.pdf")
 # note the "breaks..." part comes from the inspection of the unconstrained results

fadu_data_by_cbd = fadu_data_mean %>% 
  pivot_wider(names_from = cbd, values_from = Avg) %>% 
  mutate(ratio = cbd/control, log2ratio = log2(cbd/control))

fadu_data_by_cbd2 = matrix(filter(fadu_data_by_cbd,comp == 2)$ratio,8,12, byrow = T) 
pheatmap(fadu_data_by_cbd2,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_cbd2,cluster_rows = F,cluster_cols = F, breaks = seq(0.9, 1.5, length.out = 100), cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd_to_control comp 2000nm avg.pdf")
fadu_data_by_cbd2000 = matrix(filter(fadu_data_by_cbd,comp == 2)$log2ratio,8,12, byrow = T) 
pheatmap(fadu_data_by_cbd2000,cluster_rows = F,cluster_cols = F, breaks = seq(-1,1, length.out = 100), cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd_to_control log2ratio 2000nm avg.pdf")

fadu_data_by_cbd200 = matrix(filter(fadu_data_by_cbd,comp == 200)$ratio,8,12, byrow = T) 
pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F, breaks = seq(0.9, 1.5, length.out = 100), cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd_to_control comp 200nm avg.pdf")
fadu_data_by_cbd200 = matrix(filter(fadu_data_by_cbd,comp == 200)$log2ratio,8,12, byrow = T) 
pheatmap(fadu_data_by_cbd200,cluster_rows = F,cluster_cols = F, breaks = seq(-1,1, length.out = 100), cellwidth = 15, cellheight = 15, fontsize = 8,  filename = "FADU run1 cbd_to_control log2ratio 200nm avg.pdf")


fadu_data_by_comp = fadu_data_mean %>% 
  pivot_wider(names_from = comp, values_from = Avg)
colnames(fadu_data_by_comp) = c("Well", "cbd",  "comp_2000",    "comp_200")
fadu_data_by_comp = fadu_data_by_comp %>% mutate(ratio = comp_2000/comp_200, log2ratio = log2(comp_2000/comp_200))

fadu_data_by_comp_control = matrix(filter(fadu_data_by_comp,cbd == "control")$ratio,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
pheatmap(fadu_data_by_comp_control, cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_comp_control, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, filename = "FADU run1 no_cbd high_to_low comp avg.pdf")

fadu_data_by_comp_control = matrix(filter(fadu_data_by_comp,cbd == "control")$log2ratio,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
pheatmap(fadu_data_by_comp_control, cluster_rows = F,cluster_cols = F)
pheatmap(fadu_data_by_comp_control, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, breaks = seq(-1,1, length.out = 100), filename = "FADU run1 no_cbd high_to_low log2ratio avg.pdf")

fadu_data_by_comp_cbd = matrix(filter(fadu_data_by_comp,cbd == "cbd")$ratio,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
pheatmap(fadu_data_by_comp_cbd, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, filename = "FADU run1 cbd high_to_low comp avg.pdf")
fadu_data_by_comp_cbd = matrix(filter(fadu_data_by_comp,cbd == "cbd")$log2ratio,8,12, byrow = T) 
# running the line below breaks pheatmap permanently from showing the results, but files can be still saved
pheatmap(fadu_data_by_comp_cbd, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15, breaks = seq(-1,1, length.out = 100), filename = "FADU run1 cbd high_to_low log2ratio avg.pdf")

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
fadu_data_by_cbd_plot = fadu_data_mean %>% 
  pivot_wider(names_from = cbd, values_from = c(Avg,SD))

fadu_data_by_cbd_plot = fadu_data_by_cbd_plot %>% 
  mutate( well_type = case_when(Well %in% fadu_dmso_wells ~ "dmso", 
                                Well %in% fadu_staur_wells ~ "PC",
                                TRUE ~ "test_compounds"))

fadu_data_by_cbd_plot %>% filter(comp == 2)  %>% 
  ggplot(aes(x = Avg_control, y = Avg_cbd, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_cbd - SD_cbd, ymax = Avg_cbd + SD_cbd), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control - SD_control, xmax = Avg_control + SD_control), height = 0.1) +
  coord_cartesian(xlim = c(0,250000),ylim = c(0,250000))
 
ggsave("control vs cbd run1 conc2000 detailed with errorbars.png")

fadu_data_by_cbd_plot %>% filter(comp == 200)  %>% 
  ggplot(aes(x = Avg_control, y = Avg_cbd, color = well_type))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = Avg_cbd - SD_cbd, ymax = Avg_cbd + SD_cbd), width = 0.1) +
  geom_errorbarh(aes(xmin = Avg_control - SD_control, xmax = Avg_control + SD_control), height = 0.1) +
  coord_cartesian(xlim = c(0,250000),ylim = c(0,250000))

ggsave("control vs cbd run1 conc200 detailed with errorbars.png")

# - - - - - - - - - - - - - - - - - - - - - - - - - - -

ggplot(fadu_data, aes(x = V_pl_1)) +
          geom_density()
fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
          ggplot(aes(x = CTB,fill = plate, col = plate)) +
          geom_density(alpha = 0.5)
fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
          ggplot(aes(x = CTB,fill = plate, col = plate)) +
          geom_density(alpha = 0.2)
ggsave("FADU CTB density by plates.png")

fadu_data_calc = fadu_data %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# same as above, a bit shorter
fadu_data_calc = fadu_data %>% mutate(veh = rowMeans(.[2:3]) , cbd = rowMeans(.[4:5]), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))

fadu_data_calc %>% write_csv("fadu_data_calc.csv")
fr_change = matrix(fadu_data_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change,cluster_rows = F,cluster_cols = F)
log_fr_change = matrix(fadu_data_calc$log_fr_change,16,24, byrow = T)
pheatmap(log_fr_change,cluster_rows = F,cluster_cols = F)

fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change> 1.28) %>% nrow() # 2
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change < 0.76) %>% nrow() # 1
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change> 0.5) %>% nrow()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change < -0.5) %>% nrow()

fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max()
# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# and it is NOT a normal distribution although looks not too distorted
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = fr_change)) + geom_density()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()


fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% max()

fadu_hits = fadu_data_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > 1.27844 | fr_change < 0.7594916)) 
fadu_hits %>% write_csv("fadu_hits.csv")


fadu_data %>% filter(Well %in% fadu_controls) %>% select(V_pl_1,V_pl_2) %>% cor() #0.3645661
fadu_data %>% filter(Well %in% fadu_controls) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.3401433
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(V_pl_1,V_pl_2) %>% cor() #0.7407735
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.6920442
fadu_data %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()

fadu_data %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = V_pl_1, y = CBD_pl_1))+
  geom_point()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point() + geom_smooth()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          ggplot(aes(x = veh, y = cbd))+
          geom_point() + geom_smooth(method = 'lm')
fadu_data_calc_for_plot_comp = fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd) %>% mutate(well = "compounds")
fadu_data_calc_for_plot_control = fadu_data_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd) %>% mutate(well = "controls")
fadu_data_calc_for_plot_gem = fadu_data_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd) %>% mutate(well = "gem")
fadu_data_calc_for_plot_doc = fadu_data_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd) %>% mutate(well = "doc")

fadu_data_calc_for_plot = bind_rows(fadu_data_calc_for_plot_control,
                                    fadu_data_calc_for_plot_comp,
                                    fadu_data_calc_for_plot_gem,
                                    fadu_data_calc_for_plot_doc)
fadu_data_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')

#  * * *
fadu_data2 = read_excel("FADU_6_plates_run_2.xlsx")
fadu_data2 %>% nrow() # 384 which is correct



# fadu_v1 = matrix(fadu_data$V_pl_1,16,24, byrow = T)
# pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
# fadu_v2 = matrix(fadu_data$V_pl_2,16,24, byrow = T)
# pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
# fadu_cbd1 = matrix(fadu_data$CBD_pl_1,16,24, byrow = T)
# pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
# fadu_cbd2 = matrix(fadu_data$CBD_pl_2,16,24, byrow = T)
# pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)
# 
# ggplot(fadu_data, aes(x = V_pl_1)) +
#   geom_density()
# fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
# fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.5)
# fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.2)
# ggsave("FADU CTB density by plates.png")

# fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
#                                       fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# same as above, a bit shorter
fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(.[2:4]) , cbd = rowMeans(.[5:7]), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
fadu_data2_calc$cbd_sd = apply(fadu_data2_calc[5:7],1,sd)
fadu_data2_calc$veh_sd = apply(fadu_data2_calc[2:4],1,sd)

fadu_data2_calc %>% write_csv("fadu_data2_calc.csv")
fr_change2 = matrix(fadu_data2_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change2,cluster_rows = F,cluster_cols = F)
# log_fr_change = matrix(fadu_data_calc$log_fr_change,16,24, byrow = T)
# pheatmap(log_fr_change,cluster_rows = F,cluster_cols = F)

fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min() # 0.6357736
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max() # 1.136606
fadu_data2_calc %>% filter(Well %in% fadu_controls, fr_change> 1.28) %>% nrow() # 2
fadu_data2_calc %>% filter(Well %in% fadu_controls, fr_change < 0.76) %>% nrow() # 1
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% min() # -0.653
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% max() #0.185
fadu_data2_calc %>% filter(Well %in% fadu_controls, log_fr_change> 0.5) %>% nrow()
fadu_data2_calc %>% filter(Well %in% fadu_controls, log_fr_change < -0.5) %>% nrow()
# fadu_controls %>% length() #96
# fadu_inner_compounds = c()
# for (i in c(2,4,6,8,10,12,14)) {
#   fadu_inner_compounds = c(fadu_inner_compounds, fadu_wells[i,3:22])
#   
# }
# fadu_inner_compounds %>% length() #140
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min() # -1.132
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max()
# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# normally distributed
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = fr_change)) + geom_density()
fadu2_controls_mean = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu2_controls_sd = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
# fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()
fadu2_controls_mean # 0.8594187
fadu2_controls_sd # 0.08428393
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% min() 0.456178
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% max() 1.544866

fadu_hits2 = fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > (fadu2_controls_mean + 2* fadu2_controls_sd)| fr_change < (fadu2_controls_mean - 2* fadu2_controls_sd))) 
fadu_hits2 %>% write_csv("fadu_hits2.csv")

fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point()
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point() + geom_smooth(method = "lm")
ggsave("veh vs cbd run2.png")

fadu_data2_calc_for_plot_comp = fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd,veh_sd,cbd_sd) %>% mutate(well = "compounds")
fadu_data2_calc_for_plot_control = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd,veh_sd,cbd_sd) %>% mutate(well = "controls")
fadu_data2_calc_for_plot_gem = fadu_data2_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd,veh_sd,cbd_sd) %>% mutate(well = "gem")
fadu_data2_calc_for_plot_doc = fadu_data2_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd,veh_sd,cbd_sd) %>% mutate(well = "doc")

fadu_data2_calc_for_plot = bind_rows(fadu_data2_calc_for_plot_control,
                                    fadu_data2_calc_for_plot_comp,
                                    fadu_data2_calc_for_plot_gem,
                                    fadu_data2_calc_for_plot_doc)
fadu_data2_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')
ggsave("veh vs cbd run2 detailed.png")
fadu_data2_calc_for_plot  %>% 
  ggplot(aes(x = veh, y = cbd, color = well))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = cbd - cbd_sd, ymax = cbd + cbd_sd), width = 0.2) +
  geom_errorbarh(aes(xmin = veh - veh_sd, xmax = veh + veh_sd), height = 0.2)

fadu_data2_calc_for_plot  %>% 
  ggplot(aes(x = veh, y = cbd, color = well))+
  geom_point() + geom_smooth(method = 'lm') +
  geom_errorbar(aes(ymin = cbd - cbd_sd, ymax = cbd + cbd_sd), width = 0.1) +
  geom_errorbarh(aes(xmin = veh - veh_sd, xmax = veh + veh_sd), height = 0.1) +
  geom_point(data = fadu_hits2, aes(x = veh, y = cbd, color = "hits"),color = "black") +
  coord_cartesian(xlim = c(0,500000),ylim = c(0,500000))

ggsave("veh vs cbd run2 detailed with errorbars.png")

# - - - - - - - - - - - - - -
#  * * *
fadu_data2n = read_excel("FADU run2 nuclei count.xlsx")
fadu_data2n %>% nrow() # 2304  which is 384 * 6 which is correct
fadu_data2n = fadu_data2n %>% pivot_wider(names_from = Plate, values_from = Nuclei)


# fadu_v1 = matrix(fadu_data$V_pl_1,16,24, byrow = T)
# pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
# fadu_v2 = matrix(fadu_data$V_pl_2,16,24, byrow = T)
# pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
# fadu_cbd1 = matrix(fadu_data$CBD_pl_1,16,24, byrow = T)
# pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
# fadu_cbd2 = matrix(fadu_data$CBD_pl_2,16,24, byrow = T)
# pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)
# 
# ggplot(fadu_data, aes(x = V_pl_1)) +
#   geom_density()
# fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
# fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.5)
# fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.2)
# ggsave("FADU CTB density by plates.png")

# fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
#                                       fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# same as above, a bit shorter
fadu_data2n_calc = fadu_data2n %>% mutate(veh = rowMeans(.[5:7]) , cbd = rowMeans(.[2:4]), 
fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# note the column order is different
# 
# fadu_data2n_calc = fadu_data2n %>% rowwise() %>% mutate(veh = mean(c_across(5:7)) , cbd = mean(c_across(2:4)),
#                                                         veh_sd = sd(c_across(5:7)) , cbd_sd = sd(c_across(2:4)),
#                        fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# yet a third way

# maybe the easiest way - on the top
fadu_data2n_calc$veh_sd = apply(fadu_data2n_calc[5:7],1,sd)
fadu_data2n_calc$cbd_sd = apply(fadu_data2n_calc[2:4],1,sd)
fadu_data2n_calc = fadu_data2n_calc %>% mutate(cbd_min = cbd - 2*cbd_sd, cbd_max = cbd + 2*cbd_sd)
fadu_data2n_calc %>% head()

fadu_data2n_calc %>% write_csv("fadu_data2_nuclei_calc.csv")
fr_change2n = matrix(fadu_data2n_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change2n,cluster_rows = F,cluster_cols = F)
# saved as "FADU run 2 nuclei count fr_change.png"

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min() # 0.3801222
fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max() # 0.634187
fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min() # -1.617
fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max() # 0.035

# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# NOT normally distributed but looks just fine

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% 
          ggplot(aes(x = fr_change)) + geom_density()
ggsave("FADU run2 nuclei count controls fr_change distribution.png")


fadu2n_controls_mean = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu2n_controls_sd = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
# fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()
fadu2n_controls_mean # 0.5182667
fadu2n_controls_sd # .04414386

fadu_hits2n = fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > (fadu2n_controls_mean + 2* fadu2n_controls_sd)| fr_change < (fadu2n_controls_mean - 2* fadu2n_controls_sd))) 
fadu_hits2n %>% write_csv("fadu_hits2_nuclei.csv")


# fadu_data %>% filter(Well %in% fadu_controls) %>% select(V_pl_1,V_pl_2) %>% cor() #0.3645661
# fadu_data %>% filter(Well %in% fadu_controls) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.3401433
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(V_pl_1,V_pl_2) %>% cor() #0.7407735
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.6920442
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# 
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = CBD_pl_1))+
#   geom_point()

fadu_data2n_calc_for_plot_comp = fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd ) %>% mutate(well = "compounds")
fadu_data2n_calc_for_plot_control = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "controls")
fadu_data2n_calc_for_plot_gem = fadu_data2n_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "gem")
fadu_data2n_calc_for_plot_doc = fadu_data2n_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "doc")

fadu_data2n_calc_for_plot = bind_rows(fadu_data2n_calc_for_plot_control,
                                     fadu_data2n_calc_for_plot_comp,
                                     fadu_data2n_calc_for_plot_gem,
                                     fadu_data2n_calc_for_plot_doc)
fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')

fadu_data2n_calc_for_plot  %>%  
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd_min, ymax = cbd_max), width = 0.2)

fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd - 2*cbd_sd, ymax = cbd + 2*cbd_sd), width = 0.2)
# same result
ggsave("veh vs cbd run2 nuclei detailed_2.png")

fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd - 2*cbd_sd, ymax = cbd + 2*cbd_sd), width = 0.2) +
          geom_errorbarh(aes(xmin = veh - 2*veh_sd, xmax = veh + 2*veh_sd), height = 0.2) +
          coord_cartesian(xlim = c(2000,16000))
ggsave("veh vs cbd run2 nuclei detailed_3.png")

         
# -------------------


ggsave("veh vs cbd run2 nuclei detailed_2.png")

fadu_hits2_comb_wells = append(fadu_hits2$Well,fadu_hits2n$Well) %>% unique()
fadu_hits2_extr = fadu_data2_calc %>% filter(Well %in% fadu_hits2_comb_wells) %>% 
          select(Well, fr_change) %>% rename(fr_change_CTB = fr_change)
fadu_hits2n_extr = fadu_data2n_calc %>% filter(Well %in% fadu_hits2_comb_wells) %>% 
          select(Well, fr_change) %>% rename(fr_change_Nuc = fr_change)
fadu_hits2_comparison = left_join(fadu_hits2_extr,fadu_hits2n_extr)
pheatmap(as.matrix(select(fadu_hits2_comparison, -Well)),cluster_rows = F,cluster_cols = F)
fadu_hits2_comparison %>% write_csv("fadu_hits2_comparison.csv")
