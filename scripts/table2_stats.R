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

dat.l <- dat.org %>% group_by(temp) %>% mutate(CH4_emis_C = CH4_emis * 12/16.04, CO2_emis_C = CO2_emis * 12/44.01) %>%
pivot_longer(cols = c('CH4_emis_C', 'CO2_emis_C'), 
               names_to = 'comp', values_to = 'value')

#add surface respiration rates g day m2
emis.cum.storage <- dat.l %>% group_by(reactor, comp, temp, gas) %>%
  mutate(cum = sum(approx(x = date, y = value, xout = day_spaced)$y)) %>% 
  distinct(cum, keep_all = TRUE) %>% 
  mutate(cum = cum/1000) %>%
  group_by(temp, gas, comp) %>% 
  summarise(across(cum, .fns = list(mean = mean, sd = sd)))
end.storage <- max(dat.l$day) 

emis.stat.storage <- dat.l %>% group_by(reactor, comp, temp, gas) %>%
  mutate(cum = sum(approx(x = date, y = value, xout = day_spaced)$y)) %>% 
  mutate(cum = cum/1000) %>% 
  summarise(cum = cum[day == end.storage]) %>%
  group_by(comp, temp) %>% 
  do({
    fit <- aov(cum ~ gas, data = .)
     posthoc <- TukeyHSD(fit)
    bind_rows(tidy(fit), tidy(posthoc))
  })
  

## biogas 
rho_CH4 <- 0.665 # at 22 deg C, kg/m3
rho_CO2 <- 1.802

dat.info <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))
weights_scale <- dat.info %>% select(wet_weight, day, reactor) %>% filter(day >= 283, day< 283.1) %>% 
  mutate(weights_scale = wet_weight)

end <- max(dat.info$day)

emis.stat.biogas <- dat.info %>% 
  select(day, datetime, bio_wet_weight, mass_vented, vol_vented, mass_leaked, x_CH4, reactor, temp, gas) %>% 
  filter(day > 0) %>% mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>%
  mutate(mass_leaked = ifelse(mass_leaked < 0, 0, mass_leaked), x_CH4 = ifelse(x_CH4 < 0, 0, ifelse(x_CH4 > 1, 1, x_CH4))) %>% # correct impossible mol fractions and negative leakage
  mutate(vol_leak = vol_vented/mass_vented * mass_leaked, vol_vent_leak = vol_leak + vol_vented, 
         vol_ch4 = vol_vent_leak * x_CH4, vol_co2 = vol_vent_leak * (1 - x_CH4)) %>% 
  mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>% 
  mutate(weights_scale = rep(weights_scale$weights_scale, length = nrow(dat.info[dat.info$day >= 283,]))) %>% 
  group_by(reactor) %>%  
  mutate(cum_CH4 = cumsum(vol_ch4)/bio_wet_weight * weights_scale * rho_CH4/1000 * 12/16.04, 
         cum_CO2 = cumsum(vol_co2)/bio_wet_weight * weights_scale * rho_CO2/1000 * 12/44.01) %>%
  pivot_longer(c('cum_CH4', 'cum_CO2'), names_to = "comp", values_to = "cum") %>%
  group_by(reactor, temp, gas, comp) %>% summarise(cum = cum[day == end]) %>% 
  group_by(comp, temp) %>%
  do({
    fit <- aov(cum ~ gas, data = .)
    posthoc <- TukeyHSD(fit)
    bind_rows(tidy(fit), tidy(posthoc))
  })


## combined
emis.cum.storage <- dat.l %>% group_by(reactor, comp, temp, gas) %>%
  mutate(cum = sum(approx(x = date, y = value, xout = day_spaced)$y)) %>% 
  distinct(cum, keep_all = TRUE) %>% 
  mutate(cum = cum/1000, experiment = "storage") %>% select(-keep_all)


emis.cum.biogas <- dat.info %>% 
  select(day, datetime, bio_wet_weight, mass_vented, vol_vented, mass_leaked, x_CH4, reactor, temp, gas) %>% 
  filter(day > 0) %>% mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>%
  mutate(mass_leaked = ifelse(mass_leaked < 0, 0, mass_leaked), x_CH4 = ifelse(x_CH4 < 0, 0, ifelse(x_CH4 > 1, 1, x_CH4))) %>% # correct impossible mol fractions and negative leakage
  mutate(vol_leak = vol_vented/mass_vented * mass_leaked, vol_vent_leak = vol_leak + vol_vented, 
         vol_ch4 = vol_vent_leak * x_CH4, vol_co2 = vol_vent_leak * (1 - x_CH4)) %>% 
  mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>% 
  mutate(weights_scale = rep(weights_scale$weights_scale, length = nrow(dat.info[dat.info$day >= 283,]))) %>% 
  group_by(reactor) %>%  
  mutate(CH4_emis_C = cumsum(vol_ch4)/bio_wet_weight * weights_scale * rho_CH4/1000 * 12/16, 
         CO2_emis_C = cumsum(vol_co2)/bio_wet_weight * weights_scale * rho_CO2/1000 * 12/44) %>%
  pivot_longer(c('CH4_emis_C', 'CO2_emis_C'), names_to = "comp", values_to = "cum") %>% 
  mutate(experiment = "biogas") %>% 
  filter(day == end) %>% select(colnames(emis.cum.storage)) %>% mutate(temp = as.numeric(temp))

emis.cum.combined <- bind_rows(emis.cum.storage, emis.cum.biogas) %>% 
  group_by(reactor, temp, gas, comp) %>% summarise(cum = sum(cum)) %>% 
  mutate(experiment = "combined")

emis.stat.combined <- emis.cum.combined %>% 
  group_by(comp, temp) %>%
  do({
    fit <- aov(cum ~ gas, data = .)
    posthoc <- TukeyHSD(fit)
    bind_rows(tidy(fit), tidy(posthoc))
  })


emis_all <- rbind(emis.cum.storage, emis.cum.biogas, emis.cum.combined)
stat_all <- rbind(emis.stat.storage, emis.stat.biogas, emis.stat.combined)

write.csv(emis_all, "../output/emis_cum_dat.csv", row.names = F)
write.csv(stat_all, "../output/emis_cum_stat.csv", row.names = F)


table2.storage <- emis.cum.storage %>% group_by(temp, gas, comp) %>% summarise(across(cum, .fns = list(mean = mean, std = sd)))
table2.biogas <- emis.cum.biogas %>% group_by(temp, gas, comp) %>% summarise(across(cum, .fns = list(mean = mean, std = sd)))
table2.combined <- emis.cum.combined %>% group_by(temp, gas, comp) %>% summarise(across(cum, .fns = list(mean = mean, std = sd)))

write.csv(rbind(table2.storage, table2.biogas, table2.combined), '../output/table2.csv', row.names = F)
