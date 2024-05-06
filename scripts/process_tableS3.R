rm(list = ls())

library("readxl")
library('openxlsx')
library('ggplot2')
library('tidyr')
library('dplyr')
library('broom')

#prep
dat.inf <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))

conc_cols <- c('dm','vs','lig','cel','hem','ndf', 'adf', 'lip','tan','tn','vfa')

dat.inf.f <- dat.inf %>% select(c('reactor', 'temp','gas','datetime','day', 'wet_weight', 'dm','vs','lig','cel','hem', 'ndf', 'adf', 'lip','tan','tn','vfa','pH')) %>%
  mutate(dm = dm * 10, vs = vs * dm, ndf = ndf/100 * dm, lig = lig/100 * dm, cel = cel/100 * dm, hem = hem/100 * dm, adf = adf/100 * dm, lip = lip/100 * dm) %>% 
  mutate(across(all_of(conc_cols), ~ . / 1000 * wet_weight, .names = "{.col}")) %>%
  mutate(CP = (tn-tan)*6.25)

stat_cols <- c(conc_cols, 'CP', 'pH', 'wet_weight')
#start and end concentrations
table1 <- dat.inf.f %>% filter(!is.na(wet_weight), temp != 'none') %>% 
          group_by(day, temp, gas) %>% 
          summarise(across(all_of(stat_cols), .fns = list(mean = mean, sd = sd)))

write.xlsx(table1, '../output/table1.xlsx')

#  group_by(comp) %>% 
table1_stats <- dat.inf.f %>% filter(!is.na(wet_weight), temp != 'none') %>% 
          pivot_longer(all_of(stat_cols), names_to = 'comp', values_to = 'value') %>%
  filter(!is.na(value)) %>%
  group_by(comp) %>% 
  do({
    fit <- aov(value ~ gas * temp * day, data = .)
    posthoc <- TukeyHSD(fit)
    bind_rows(tidy(fit), tidy(posthoc))
  })

write.xlsx(table1_stats, '../output/table1_stats.xlsx')
####