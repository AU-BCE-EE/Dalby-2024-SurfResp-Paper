rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(marelac)
library(ReacTran)
library(readxl)
library(viridis)

source('pde_fun.R')

dat <- read_excel('../data/O2_profiles_corr.xlsx')
optim_pars <- read.csv('../output/optim_pars.csv')
pars <- cbind(optim_pars, D = 2.2 * 10^-5, kH = 0.0013/1000*10^6)
input_dat <- dat %>% mutate(id = as.factor(id)) %>% filter(exclude == 0, depth < 1000)

plot_dat <- data.frame()

for(i in unique(input_dat$id)){
  parms <- pars[pars$id == i,]
  dat <- input_dat[input_dat$id == i,] %>% select(-manure, -exclude, -mV) %>% filter(umol > 0) %>% 
  mutate(r = parms$r, KO2 = parms$KO2, type = "data", flux = pars$D[pars$id == i] * max(abs(diff(umol/1000))/abs(diff(depth/10000))))
  x <- max(dat$depth)
  model1 <- cbind(pde_fun(parms = parms, x = x), id = i, type = "model", man_name = dat$man_name[1], time = dat$time[1], r = parms$r[1], KO2 = parms$KO2[1]) %>% 
  rename('umol' = 'O2')
    
  
  plot_dat <- rbind(plot_dat, model1, dat)
  
}

plot_dat <- plot_dat[!grepl('Cattle', plot_dat$man_name),]

fig_profiles <- ggplot(plot_dat, aes(umol, y = depth, group = id, col = time)) +
  geom_path(data = subset(plot_dat, type == "data")) +  # Add points for "data" category
  scale_y_continuous(trans = "reverse") +
  scale_x_continuous(position = "top") +
  labs(y = "Depth (µm)", x = expression("Oxygen concentration (µmol"~L^{-1}*")"), 
       col = "Time since mixing (d)") + 
  facet_wrap(~man_name) + 
  scale_color_viridis(option = "D")  +
  theme_bw()

png('../figures/fig_profiles.png',  width = 18/2.54, height = 12/2.54, units = 'in', res = 600)
grid::grid.draw(fig_profiles)
dev.off()

svg('../figures/fig_profiles.svg',  width = 18/2.54, height = 12/2.54)
grid::grid.draw(fig_profiles)
dev.off()

pdf('../figures/fig_profiles.pdf',   width = 17/2.54, height = 12/2.54)
grid::grid.draw(fig_profiles)
dev.off()

flux_stat <- plot_dat %>% group_by(type, id) %>% 
  summarise(Vmax = mean(r), KO2 = mean(KO2), time = mean(time), slurry = man_name[1], flux_g_O2_m2_day = mean(flux) * 3600*24*10^-6*32*10000, flux_g_CO2_m2_day = flux_g_O2_m2_day * 0.33)

write.csv(flux_stat, '../output/flux_stat.csv', row.names = F)

# remove strange profiles, id c(30, 31, 34, 35 ,36, 56, 70, 71)
flux_dat <- read.csv('../output/flux_stat.csv')
flux_dat$animal <- 'Pig'
flux_dat$animal[grepl('Cattle',flux_dat$slurry)] <- 'Cattle'

flux_dat_summary <- flux_dat %>% 
  filter(!id %in% c(30, 31, 34, 35 ,36, 56, 70, 71), time > 0.7, type == 'model') %>% 
  filter(animal == 'Pig') %>%
  group_by(slurry, animal) %>%
  summarise(across(c('Vmax', 'KO2','flux_g_O2_m2_day','flux_g_CO2_m2_day'),
                   .fns = list(mean = mean, sd = sd)))
flux_dat_summary[is.na(flux_dat_summary)] <- 0

# uncertainties are propagated

flux_dat_prop <- flux_dat_summary %>% 
  mutate(Vmax_e = ((Vmax_sd/Vmax_mean)*100)^2, 
         KO2_e = ((KO2_sd/KO2_mean)*100)^2,
         flux_O2_e = ((flux_g_O2_m2_day_sd/flux_g_O2_m2_day_mean)*100)^2,
         flux_CO2_e = ((flux_g_CO2_m2_day_sd/flux_g_CO2_m2_day_mean)*100)^2) %>% 
  group_by(animal) %>% summarise(across(c('Vmax_e', 'KO2_e','flux_O2_e','flux_CO2_e'),
                   function(x) sqrt(sum(x))/100),
                   across(ends_with('mean'), function(x) mean(x))) %>% 
  mutate(Vmax_abs_e = Vmax_e * Vmax_mean, 
         KO2_abs_e = KO2_e * KO2_mean, 
         flux_O2_abs_e = flux_O2_e * flux_g_O2_m2_day_mean, 
         flux_CO2_abs_e = flux_CO2_e * flux_g_CO2_m2_day_mean)
                
 


