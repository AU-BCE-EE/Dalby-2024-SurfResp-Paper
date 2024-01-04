# surface respiration estimate fresh manure
rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')

dat.org <- data.frame(read_excel("../data/dat_herald.xlsx"))

R <- 0.082057
T <- 298
f_CO2 <- 1/1000000 /(R * T) * 60 * 24 * 1000# mmol/day
f_CH4 <- 1/1000000 /(R * T) * 60 * 24 * 1000 # mmol/day

dat.mod <- dat.org %>% mutate(CO2_emis = (CO2 - 450) * flow * f_CO2, CH4_emis = (CH4 - 2) * flow * f_CH4, ratio = CH4_emis / (CO2_emis + CH4_emis)) %>%
 filter(CO2_emis > 0, CH4_emis > 0) %>% group_by(id) %>% summarise(ratio = mean(ratio, na.rm = T), CO2 = mean(CO2, na.rm = T), CH4 = mean(CH4, na.rm = T), 
                            CO2_emis = mean(CO2_emis, na.rm = T), CH4_emis = mean(CH4_emis, na.rm = T))

dat.mod.long <- dat.org %>% mutate(CO2_emis = (CO2 - 450) * flow * f_CO2, CH4_emis = (CH4 - 2) * flow * f_CH4, ratio = CH4_emis / (CO2_emis + CH4_emis)) %>%
  filter(CO2_emis > 0, CH4_emis > 0) %>%  pivot_longer(cols = c('CH4_emis', 'CO2_emis', 'ratio'), names_to = 'comp', values_to = 'value') %>% filter(id != 7)

dat.mod.long$id <- as.factor(dat.mod.long$id) 

ggplot(dat.mod.long, aes(x = date.time, y = value, col = id)) + geom_point() + facet_wrap(~comp, scales = "free") 

