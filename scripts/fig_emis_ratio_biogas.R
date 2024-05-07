# surface respiration estimate
rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')
library('gridExtra')
library('lubridate')

#prep
dat_gas <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "all"))
dat.info <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))

dat.info.f <- dat.info %>% select(c('reactor', 'temp','gas','datetime','day', 'wet_weight', 'dm','vs','lig','cel','hem','lip','tan','tn','vfa','pH')) %>%
  mutate(dm = dm * 10, vs = vs * dm)
mass <- dat.info.f$wet_weight[dat.info.f$day == 0]/1000 # kg
mass <- mass[!is.na(mass)]
vs_mass <- dat.info.f$vs[dat.info.f$day == 0 & dat.info.f$gas != 'none']/1000
vs_mass <- vs_mass[!is.na(vs_mass)] 

# kg
dat_gas <- dat_gas[dat_gas$reactor != 'bg',]
dat_gas$mass <- mass
dat_gas$vs_mass <- vs_mass


f_CH4 <- 1/1000000 * dat_gas$cor_flow/(0.082057 * 293) * 60 * 24 * 16.04  / dat_gas$vs_mass # mg pr. kg slurry pr day
f_CO2 <- 1/1000000 * dat_gas$cor_flow/(0.082057 * 293) * 60 * 24 * 44.01  / dat_gas$vs_mass # mg pr. kg slurry pr day
  
CO2_bg <- 430
CH4_bg <- 2 
  
ifelse(dat_gas$gas == "air", dat_gas$CO2_emis <- (as.numeric(dat_gas$co2) - CO2_bg) * f_CO2, dat_gas$CO2_emis <- dat_gas$co2 * f_CO2)  
ifelse(dat_gas$gas == "air", dat_gas$CH4_emis <- (as.numeric(dat_gas$ch4) - CH4_bg) * f_CH4, dat_gas$CH4_emis <- dat_gas$ch4 * f_CH4)  

dat.mod <- dat_gas %>% group_by(temp, gas, day, date) %>% 
  mutate(ratio = (CH4_emis * 12.01/16.04)/(CH4_emis * 12.01/16.04 + CO2_emis * 12.01/44.01)) %>% mutate(CH4_emis_C = CH4_emis * 12.01/16.04, CO2_emis_C = CO2_emis * 12.01/44.01) %>% 
  group_by(temp, gas, day, date) %>% summarise(across(c('CO2_emis_C','CH4_emis_C', 'CO2_emis','CH4_emis', 'ratio'), .fns = list(mean = mean, sd = sd), na.rm = TRUE)) 

dat.mod.long <- dat.mod %>% pivot_longer(cols = c('CH4_emis_C_mean', 'CO2_emis_C_mean', 'CH4_emis_C_sd', 'CO2_emis_C_sd'), names_to = 'comp', values_to = 'value')

dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 
dat.mod.both$temp <- as.factor(dat.mod.both$temp)

# plot emission rates
new.lab <- as_labeller(c(air = "Air", n2 = "N[2]", CH4_emis_C_mean = "CH[4]", CO2_emis_C_mean = "CO[2]"), label_parsed)

fig_emis <- ggplot(dat.mod.both, aes(x = day, y = value, col = temp)) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(comp~gas, labeller = new.lab, scales = "free_y") + 
  theme_bw() +
  theme() + 
  theme(legend.position = 'top', axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 10)) +
  labs(y = expression('Emission rate (g C kg'^{-1}~VS~'d'^{-1}*')'), 
       x = "", col = expression('Temp (\u00b0C)'), tag = 'a') +
  scale_color_manual(values = c('blue', 'red'))

# plot ratios
dat.mod.long <- dat.mod %>% pivot_longer(cols = c('ratio_mean', 'ratio_sd'), names_to = 'comp', values_to = 'value')

dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 

dat.mod.both$temp <- as.factor(dat.mod.both$temp)

new.lab <- as_labeller(c(air = "Air", n2 = "N[2]", ratio_mean = "C[CH4]/(C[CH4]+C[CO2])"), label_parsed)

fig_ratio <- ggplot(dat.mod.both, aes(x = day, y = value, col = temp)) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(comp~gas, labeller = new.lab, scales = "free_y") + 
  scale_x_continuous(breaks = seq(0, 250, by = 50)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.position = 'none') + 
  labs(y = expression('Molar fraction'), x = "Time (d)", col = expression('Temp (\u00b0C)'))+
  scale_color_manual(values = c('blue', 'red'))

# plot biogas
source('biogas.R')

GDout_tb <- GDout %>% group_by(gas, temp, time) %>% 
  summarize(across(c(mCH4_C, mCO2_C, mCH4_C_cum, mCO2_C_cum), .fns = list(mean = mean, sd = sd))) %>%  
  filter(gas != 0)

GDout_tb$compound <- 'CH4'

new.lab = as_labeller(c(air = 'Air', n2 = 'N[2]', CH4 = 'CH[4]'), label_parsed)

fig_biogas <- ggplot(GDout_tb, aes(time, mCH4_C_mean, col = temp)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin = mCH4_C_mean - mCH4_C_sd, ymax = mCH4_C_mean + mCH4_C_sd, x = time)) +
  facet_grid(compound~gas, labeller = new.lab) + 
  labs(y = expression('Emission rate (g C kg'^{-1}~VS~'d'^{-1}*')'), x = 'Time (d)', col  = expression('Temp (\u00b0C)'), tag = 'b') +
  coord_cartesian(ylim = c(0,2.5)) + 
  theme_bw() + theme(legend.position = '', axis.title.y = element_text(size = 10)) + scale_color_manual(values = c("blue", "red"))

png('../figures/fig_emis_ratio_bio.png',  width = 18/2.54, height = 20/2.54, units = 'in', res = 600)
grid::grid.draw(rbind(ggplotGrob(fig_emis), ggplotGrob(fig_ratio), ggplotGrob(fig_biogas)))
dev.off()

