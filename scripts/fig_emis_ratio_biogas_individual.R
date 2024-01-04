# surface respiration estimate
rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')
library('gridExtra')
library('lubridate')
library("gg.gap")

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


f_CH4 <- 1/1000000 * dat_gas$cor_flow/(0.082057 * 293) * 60 * 24 * 16  / dat_gas$vs_mass # mg pr. kg slurry pr day
f_CO2 <- 1/1000000 * dat_gas$cor_flow/(0.082057 * 293) * 60 * 24 * 44  / dat_gas$vs_mass # mg pr. kg slurry pr day

CO2_bg <- 430
CH4_bg <- 2 

surf_area <- (9.5/2)^2 * pi / 10000# m^2 

ifelse(dat_gas$gas == "air", dat_gas$CO2_emis <- (as.numeric(dat_gas$co2) - CO2_bg) * f_CO2, dat_gas$CO2_emis <- dat_gas$co2 * f_CO2)  
ifelse(dat_gas$gas == "air", dat_gas$CH4_emis <- (as.numeric(dat_gas$ch4) - CH4_bg) * f_CH4, dat_gas$CH4_emis <- dat_gas$ch4 * f_CH4)  

surf_conv <- mean(mass/surf_area) # convert mg / kg slurry / day to mg / m2 / day

dat.mod <- dat_gas %>% mutate(across(contains('CH4_emis'), function(x) x * 12/ 16), across(contains('CO2_emis'), function(x) x * 12/44))
dat.mod.long <- dat.mod %>% pivot_longer(cols = c('CH4_emis', 'CO2_emis'), names_to = 'comp', values_to = 'value')
dat.mod.long$temp <- as.factor(dat.mod.long$temp)

# plot emission rates
new.lab <- as_labeller(c(air = "Air", n2 = "N[2]", CH4_emis = "CH[4]", CO2_emis = "CO[2]"), label_parsed)

fig_emis <- ggplot(dat.mod.long, aes(x = day, y = value, col = temp, group = reactor)) + 
  geom_line(data = dat.mod.long[!is.na(dat.mod.long$value),]) +
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

dat.mod <- dat_gas %>% group_by(temp, gas, day, date) %>% 
  mutate(ratio = (CH4_emis * 12.01/16.04)/(CH4_emis * 12.01/16.04 + CO2_emis * 12.01/44.01)) %>% mutate(CH4_emis_C = CH4_emis * 12/16, CO2_emis_C = CO2_emis * 12/44) %>% 
  group_by(temp, gas, day, date) %>% summarise(across(c('CO2_emis_C','CH4_emis_C', 'CO2_emis','CH4_emis', 'ratio'), .fns = list(mean = mean, sd = sd), na.rm = TRUE)) 

dat.mod.long <- dat.mod %>% pivot_longer(cols = c('ratio_mean', 'ratio_sd'), names_to = 'comp', values_to = 'value')
dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 

dat.mod.both$temp <- as.factor(dat.mod.both$temp)

new.lab <- as_labeller(c(air = "Air", n2 = "N[2]", ratio_mean = "C[CH4]/(C[CH4]+C[CO2])"), label_parsed)

fig_ratio <- ggplot(dat.mod.both, aes(x = day, y = value, col = temp)) + 
  geom_line() +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(comp~gas, labeller = new.lab, scales = "free_y") + 
  scale_x_continuous(breaks = seq(0, 250, by = 50)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10)) +
  theme(legend.position = 'none') + 
  labs(y = expression('Molar fraction'), x = "Time (d)", col = expression('Temp (\u00b0C)'))+
  scale_color_manual(values = c('blue', 'red'))

rho_CH4 <- 0.665 # at 22 deg C, kg/m3
rho_CO2 <- 1.802

dat_biogas <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))
weights_scale <- dat_biogas %>% select(wet_weight, day, reactor) %>% filter(day >= 283, day< 283.1) %>% 
  mutate(weights_scale = wet_weight)
vs_mass <- dat_biogas %>% filter(day == 283) %>% mutate(vs_mass_load = dm/100 * vs * bio_wet_weight/1000) %>% select(vs_mass_load)
dat_biogas$vs_mass_load <- NULL
rep_vs <- c(0,0,0, c(rep(t(vs_mass), length.out = nrow(dat_biogas)-3)))
dat_biogas$vs_mass_load <- rep_vs
dat.s <- dat_biogas %>% 
  
  select(day, datetime, bio_wet_weight, mass_vented, vol_vented, vs_mass_load, mass_leaked, x_CH4, reactor, temp, gas) %>% 
  filter(day > 0) %>% mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>%
  mutate(mass_leaked = ifelse(mass_leaked < 0, 0, mass_leaked), x_CH4 = ifelse(x_CH4 < 0, 0, ifelse(x_CH4 > 1, 1, x_CH4))) %>% # correct impossible mol fractions and negative leakage
  mutate(vol_leak = vol_vented/mass_vented * mass_leaked, vol_vent_leak = vol_leak + vol_vented, 
         vol_ch4 = vol_vent_leak * x_CH4, vol_co2 = vol_vent_leak * (1 - x_CH4)) %>% 
  mutate_all(~ ifelse(is.na(.), 0, .), ~ ifelse(.<0, 0, .)) %>% 
  mutate(weights_scale = rep(weights_scale$weights_scale, length = nrow(dat_biogas[dat_biogas$day >= 283,]))) %>% 
  group_by(reactor) %>%  mutate(diff = diff(c(283, day))) %>% 
  mutate(CH4_rate = (vol_ch4)/vs_mass_load * rho_CH4/1000/diff, 
         CO2_rate = (vol_co2)/vs_mass_load * rho_CO2/1000/diff,
         cum_CH4 = cumsum(vol_ch4)/bio_wet_weight * weights_scale * rho_CH4/1000, 
         cum_CH4_vs = cumsum(vol_ch4)/vs_mass_load * rho_CH4/1000,
         cum_CO2 = cumsum(vol_co2)/bio_wet_weight * weights_scale * rho_CO2/1000, 
         cum_CO2_vs = cumsum(vol_co2)/vs_mass_load * rho_CO2/1000) %>% 
  filter(gas != 0, datetime != 0) %>% 
  mutate(datetime = as_datetime(datetime),  
         time_diff = datetime - min(datetime),  
         days = as.numeric(time_diff, unit = "days")) %>% 
  mutate(across(contains('CH4'), function(x) x*12/16), across(contains('CO2'), function(x) x * 12/44))#, 

#plot cumulative emission
dat.s$comp <- 'CH4'
new.lab = as_labeller(c(air = 'Air', n2 = 'N[2]', CH4 = 'CH[4]'), label_parsed)

fig_biogas<- ggplot(dat.s, aes(x = days, y = CH4_rate, col = temp, group = reactor)) + geom_line() + 
  facet_grid(comp~gas, labeller = new.lab) +
  labs(y = expression('Emission rate (g C kg'^{-1}~VS~'d'^{-1}*')'), x = 'Time (d)', col  = expression('Temp (\u00b0C)'), tag = 'b') +
  theme_bw() + theme(legend.position = '', axis.title.y = element_text(size = 10)) + scale_color_manual(values = c("blue", "red"))

fig_biogas
#gg.gap(plot=fig_biogas,
#       segments=list(c(4.5,9.5)),
#       tick_width = c(1,0.5,10),
#       rel_heights=c(0.5,0,0.2,0,0.05),
#       ylim=c(0,10))



png('../figures/fig_emis_ratio_bio_individual.png',  width = 18/2.54, height = 20/2.54, units = 'in', res = 600)
grid::grid.draw(rbind(ggplotGrob(fig_emis), ggplotGrob(fig_ratio), ggplotGrob(fig_biogas)))
dev.off()

mean_ratios <- dat.mod.both %>% group_by(temp, gas) %>% summarise(mean_ratio = mean(value))



