rm(list=ls())

### Get latest biogas package (v1.42)
##remove.packages('biogas')
##devtools::install_github('sashahafner/biogas')
##packageVersion('biogas')

library(biogas)
library(readxl)
library(dplyr)
library(ggplot2)

rho_CH4 <- 0.702 # 0 deg C, 1 ATM, kg/m3
rho_CO2 <- 1.802 # 0 deg C, 1 ATM, kg/m3
C_CH4 <- 12.0107/16.042
C_CO2 <- 12.0107/44.009

dat_biogas <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info")) %>% 
  filter(day >= 283) %>% select(reactor, gas, temp, day, bio_wet_weight, GD_weight_before, GD_weight_after, mass_vented, vol_vented, room_temp, temp_HS, p_amb) %>%
  mutate(time = day - 283) %>% select(-day)

dat <- as.data.frame(dat_biogas)

GDout <- calcBgGD(dat, temp.vol = 22, temp.grav = 34, pres.vol = 'p_amb', pres.grav = 1500, 
         id.name = 'reactor', vol.name = 'vol_vented', m.pre.name = 'GD_weight_before', m.post.name = 'GD_weight_after', 
         time.name = 'time', comp.name = 'xCH4', 
         vented.mass = TRUE, averaging = 'final', vmethod = 'grav',
         extrap = TRUE,
         addt0 = TRUE, showt0 = TRUE, comp.lim = c(0.1, 0.7), comp.sub = 'lim',
         unit.pres = 'mbar')

GDout$mCH4_C <- GDout$vCH4 * rho_CH4 * C_CH4
GDout$mCO2_C <- (GDout$vBg - GDout$vCH4) * rho_CO2 * C_CO2
GDout$mCH4_C_cum <- GDout$cvCH4 * rho_CH4 * C_CH4
GDout$mCO2_C_cum <- (GDout$cvBg - GDout$cvCH4) * rho_CO2 * C_CO2

GDout[, c('reactor', 'time', 'mass.tot', 'mass.leak')]
GDout$frac.leak <- GDout$cmass.leak / GDout$cmass.tot
GDout[GDout$time > 170, ]

GDout_tb <- GDout %>% group_by(gas, temp, time) %>% 
  summarize(across(c(mCH4_C, mCO2_C, mCH4_C_cum, mCO2_C_cum), .fns = list(mean = mean, sd = sd))) %>%  
  filter(gas != 0)

new.lab = as_labeller(c(air = 'Air', n2 = 'N[2]', CH4 = 'CH[4]'), label_parsed)

# Looks like sub value is used for 3 of 4 20 C air reactors
ggplot(subset(GDout, time > 170), aes(gas, xCH4, colour = temp)) + geom_jitter(height = 0)

ggplot(GDout_tb, aes(time, mCH4_C_mean, col = temp)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin = mCH4_C_mean - mCH4_C_sd, ymax = mCH4_C_mean + mCH4_C_sd, x = time)) +
  facet_grid(comp~gas) + 
  labs(y = expression('Emission rate (g C kg'^{-1}~VS~'d'^{-1}*')'), x = 'Time (d)', col  = expression('Temp (\u00b0C)'), tag = 'b') +
  coord_cartesian(ylim = c(0,70)) + 
  theme_bw() 

#plot cumulative emission
#png('../figures/fig_emis_ratio_bio.png',  width = 18/2.54, height = 20/2.54, units = 'in', res = 600)
#grid::grid.draw(rbind(ggplotGrob(fig_emis), ggplotGrob(fig_ratio), ggplotGrob(fig_biogas)))
#dev.off()




