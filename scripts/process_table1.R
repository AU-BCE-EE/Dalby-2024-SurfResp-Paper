#cumulative emission
rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')
library('broom')

#prep
dat.inf <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))

dat.inf.f <- dat.inf %>% select(c('reactor', 'temp','gas','datetime','day', 'wet_weight', 'dm','vs')) %>%
  mutate(dm = dm * 10, vs = vs * dm) 

dat.org <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "all"))
mass <- dat.inf.f$wet_weight[dat.inf.f$day == 0]/1000 # kg
mass <- mass[!is.na(mass)]
dat.org <- dat.org[dat.org$reactor != 'bg',]
dat.org$mass <- mass

f_CH4 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 16.04 *1000  # mg pr day
f_CO2 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 44.01 *1000  # mg pr day

CO2_bg <- 430
CH4_bg <- 2 

ifelse(dat.org$gas == "air", dat.org$CO2_emis <- (as.numeric(dat.org$co2) - CO2_bg) * f_CO2, dat.org$CO2_emis <- dat.org$co2 * f_CO2)  
ifelse(dat.org$gas == "air", dat.org$CH4_emis <- (as.numeric(dat.org$ch4) - CH4_bg) * f_CH4, dat.org$CH4_emis <- dat.org$ch4 * f_CH4)  

day_spaced <- seq(from = min(dat.org$date), to = max(dat.org$date), 
                  length.out = max(dat.org$date) - min(dat.org$date) +1)

dat.l <- dat.org %>% group_by(temp) %>% mutate(CH4_emis_C = CH4_emis * 12/16.04/0.089675, CO2_emis_C = CO2_emis * 12/44.01/0.089675) %>%
pivot_longer(cols = c('CH4_emis_C', 'CO2_emis_C'), 
               names_to = 'comp', values_to = 'value')

emis.stat.storage <- dat.l %>% group_by(reactor, comp, temp, gas) %>%
  mutate(cum = sum(approx(x = date, y = value, xout = day_spaced)$y)) %>% 
  mutate(cum = cum/1000) %>% 
  summarise(cum = cum[day == 283]) %>% mutate(incubation = 'storage') 

write.csv(emis.stat.storage, '../output/emis_stat_storage.csv', row.names = F)

## biogas 

source('biogas.R')

emis.stat.biogas <- GDout %>% pivot_longer(cols = c('mCH4_C_cum', 'mCO2_C_cum'), values_to = 'value', names_to = 'comp') %>% 
  mutate(comp = gsub('mCH4_C_cum', 'CH4_emis_C', comp)) %>% 
  mutate(comp = gsub('mCO2_C_cum', 'CO2_emis_C', comp)) %>%
  group_by(reactor, comp, temp, gas) %>% summarise(cum = max(value), incubation = 'AD') 

write.csv(emis.stat.biogas, '../output/emis_stat_biogas.csv', row.names = F)
