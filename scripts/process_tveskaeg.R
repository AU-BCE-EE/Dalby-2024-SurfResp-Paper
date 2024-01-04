rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')
library("openxlsx")


dat.org <- data.frame(read_excel("../data/TveskaegData_22.3.2023.xlsx")) #%>% 
  select(id...12, chemname...16, uivalue...18, uivalue...28,uivalue...38) %>%
  mutate(reactor = id...12, names = chemname...16, 
         rep1 =uivalue...18, rep2 = uivalue...28, rep3 = uivalue...38) %>% 
  select(reactor, names, rep1, rep2, rep3) %>% 
  pivot_longer(cols = c('rep1','rep2', 'rep3'), values_to = 'value', names_to = 'rep') %>%
  pivot_wider(names_from = 'names', values_from = 'value') %>% 
  group_by(reactor) %>% 
  summarise(across(-c('rep'), .fns = list(mean = ~mean(., na.rm =TRUE), sd = ~sd(., na.rm = TRUE))))

write.xlsx(dat.org, '../data/sorted_tveskaeg.xlsx')

dat.org <- data.frame(read_excel("../data/TveskaegData18.09.2023.xlsx")) %>% 
  select(ID...12, chemname...6, uivalue...8, uivalue...19,uivalue...30) %>%
  mutate(reactor = ID...12, names = chemname...6, 
         rep1 =uivalue...8, rep2 = uivalue...19, rep3 = uivalue...30) %>% 
  select(reactor, names, rep1, rep2, rep3) %>% 
  mutate(rep1 = as.numeric(rep1), rep2 = as.numeric(rep2), rep3 = as.numeric(rep3)) %>%
  pivot_longer(cols = c('rep1','rep2', 'rep3'), values_to = 'value', names_to = 'rep') %>%
  pivot_wider(names_from = 'names', values_from = 'value') %>% 
  group_by(reactor) %>% 
  summarise(across(-c('rep'), .fns = list(mean = ~mean(., na.rm =TRUE), sd = ~sd(., na.rm = TRUE))))

write.xlsx(dat.org, '../data/sorted_tveskaeg_biogas.xlsx')

