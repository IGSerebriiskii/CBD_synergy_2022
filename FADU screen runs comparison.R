setwd("~/Work/Flavi CBD screen/CBD_synergy_2022")
setwd("~/Work/Flavi_screen")

library(tidyverse)
library(readxl)
library(pheatmap)
library(plotly)

# calculating Z' factors

run1_CTB_summary = read_csv("384plate run 1 CTB summary.csv")
run1_CTB_summary %>% t()
run1_CTB_summary %>% t() %>% rename_with(c("bckg" ,     "dmso" ,     "lib"  ,     "staur"), everything())
run1_CTB_summary_z %>% mutate(Zprime_dmso = )

# as_tibble(cbind(nms = names(run1_CTB_summary), t(run1_CTB_summary)))
run1_CTB_summary_z = as_tibble(cbind(nms = names(run1_CTB_summary), t(run1_CTB_summary))) %>%
          .[-1, ]  %>%  `colnames<-`(c("nms" ,  run1_CTB_summary$groups)) %>% 
          mutate(across(2:5, as.numeric)) %>% select(-4)


run1_CTB_summary %>% t() %>% as_tibble() %>% select(-3) %>%
          `colnames<-`(c("bckg" ,  "dmso" ,   "staur")) %>% slice(-1) %>% 
          mutate(across(everything(), as.numeric))

run1_CTB_summary %>% t() %>% as_tibble() %>% select(-3) %>% set_names(., nm = c("bckg" ,  "dmso" ,   "staur")) %>%
          .[-1, ]


run1_CTB_summary_z = cbind(plate = paste0("plate",1:4), t(run1_CTB_summary)[2:5,1:4],t(run1_CTB_summary)[6:9,1:4]) %>% 
          as_tibble() %>%
          `colnames<-`(c("plate" ,  paste0(run1_CTB_summary$groups,"_avg"),paste0(run1_CTB_summary$groups,"_sd"))) %>% 
          mutate(across(2:9, as.numeric)) %>% select(-4, -8)

# cbind(plate = paste0(1:4), t(run1_CTB_summary)[2:5,1:4],t(run1_CTB_summary)[6:9,1:4]) %>% 
#           as_tibble() %>%
#           `colnames<-`(c("plate" ,  paste0(run1_CTB_summary$groups,"_avg"),paste0(run1_CTB_summary$groups,"_sd"))) %>% 
#           mutate(across(everything(), as.numeric))
run1_CTB_summary_z = run1_CTB_summary_z %>% 
          mutate(dmso_diff = dmso_avg - staur_avg, 
                 dmso_sd_sum = 3*(dmso_sd + staur_sd), 
                 Zprime_dmso = 1-dmso_sd_sum/dmso_diff,
                 bckg_diff = bckg_avg - staur_avg, 
                 bckg_sd_sum = 3*(bckg_sd + staur_sd), 
                 Zprime_bckg = 1-bckg_sd_sum/bckg_diff)

run1_CTB_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.matrix.data.frame() %>% write_csv("run1_CTB_summary_z.csv")
run1_CTB_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.data.frame.matrix()%>% write.csv("run1_CTB_summary_z.csv")


run1_TA_summary = read_csv("384plate run 1 total area summary.csv")
run1_TA_summary_z = cbind(plate = paste0("plate",1:4), t(run1_TA_summary)[2:5,1:4],t(run1_TA_summary)[6:9,1:4]) %>% 
  as_tibble() %>%
  `colnames<-`(c("plate" ,  paste0(run1_TA_summary$groups,"_avg"),paste0(run1_TA_summary$groups,"_sd"))) %>% 
  mutate(across(2:9, as.numeric)) %>% select(-4, -8)

run1_TA_summary_z = run1_TA_summary_z %>% 
  mutate(dmso_diff = dmso_avg - staur_avg, 
         dmso_sd_sum = 3*(dmso_sd + staur_sd), 
         Zprime_dmso = 1-dmso_sd_sum/dmso_diff,
         bckg_diff = bckg_avg - staur_avg, 
         bckg_sd_sum = 3*(bckg_sd + staur_sd), 
         Zprime_bckg = 1-bckg_sd_sum/bckg_diff)

run1_TA_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.data.frame.matrix()%>% write.csv("run1_TA_summary_z.csv")


run2_CTB_summary = read_csv("384plate run 2 CTB summary.csv")
run2_CTB_summary_z = cbind(plate = paste0("plate",1:4), t(run2_CTB_summary)[2:5,1:4],t(run2_CTB_summary)[6:9,1:4]) %>% 
  as_tibble() %>%
  `colnames<-`(c("plate" ,  paste0(run2_CTB_summary$groups,"_avg"),paste0(run2_CTB_summary$groups,"_sd"))) %>% 
  mutate(across(2:9, as.numeric)) %>% select(-4, -8)
run2_CTB_summary_z = run2_CTB_summary_z %>% 
  mutate(dmso_diff = dmso_avg - staur_avg, 
         dmso_sd_sum = 3*(dmso_sd + staur_sd), 
         Zprime_dmso = 1-dmso_sd_sum/dmso_diff,
         bckg_diff = bckg_avg - staur_avg, 
         bckg_sd_sum = 3*(bckg_sd + staur_sd), 
         Zprime_bckg = 1-bckg_sd_sum/bckg_diff)

run2_CTB_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.data.frame.matrix()%>% write.csv("run2_CTB_summary_z.csv")





run2_TA_summary = read_csv("384plate run 2 total area summary.csv")
run2_TA_summary_z = cbind(plate = paste0("plate",1:4), t(run2_TA_summary)[2:5,1:4],t(run2_TA_summary)[6:9,1:4]) %>% 
  as_tibble() %>%
  `colnames<-`(c("plate" ,  paste0(run2_TA_summary$groups,"_avg"),paste0(run2_TA_summary$groups,"_sd"))) %>% 
  mutate(across(2:9, as.numeric)) %>% select(-4, -8)

run2_TA_summary_z = run2_TA_summary_z %>% 
  mutate(dmso_diff = dmso_avg - staur_avg, 
         dmso_sd_sum = 3*(dmso_sd + staur_sd), 
         Zprime_dmso = 1-dmso_sd_sum/dmso_diff,
         bckg_diff = bckg_avg - staur_avg, 
         bckg_sd_sum = 3*(bckg_sd + staur_sd), 
         Zprime_bckg = 1-bckg_sd_sum/bckg_diff)

run2_TA_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.data.frame.matrix()%>% write.csv("run2_TA_summary_z.csv")


run3_96_CTB_summary = read_csv("96 well fadu_ run3 CTB_summary.csv")
run3_96_CTB_summary = run3_96_CTB_summary %>% filter(comp < 5000)  %>% mutate(plate = rep(c(3,4,1,2), each =4))  %>% arrange(plate) %>% 
  select(well_type, Avg, SD,plate) %>% pivot_wider(names_from = plate, values_from = c(Avg,SD))

run3_96_CTB_summary_z = cbind(plate = paste0("plate",1:4), t(run3_96_CTB_summary)[2:5,1:4],t(run3_96_CTB_summary)[6:9,1:4]) %>% 
  as_tibble() %>%
  `colnames<-`(c("plate" ,  paste0(c("bckg",  "dmso", "staur", "lib" ),"_avg"),paste0(c("bckg",  "dmso", "staur", "lib" ),"_sd"))) %>% 
  mutate(across(2:9, as.numeric)) %>% select(-5, -5)
run3_96_CTB_summary_z = run3_96_CTB_summary_z %>% 
  mutate(dmso_diff = dmso_avg - staur_avg, 
         dmso_sd_sum = 3*(dmso_sd + staur_sd), 
         Zprime_dmso = 1-dmso_sd_sum/dmso_diff,
         bckg_diff = bckg_avg - staur_avg, 
         bckg_sd_sum = 3*(bckg_sd + staur_sd), 
         Zprime_bckg = 1-bckg_sd_sum/bckg_diff)

run3_96_CTB_summary_z %>% select(plate,Zprime_dmso,Zprime_bckg) %>% t() %>% as.data.frame.matrix()%>% write.csv("run3_96_CTB_summary_z.csv")




