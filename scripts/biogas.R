devtools::install_github('sashahafner/biogas')

library(biogas)
library(readxl)
library(dplyr)
library(ggplot2)

rho_CH4 <- 0.702 # 0 deg C, 1 ATM, kg/m3
rho_CO2 <- 1.977 # 0 deg C, 1 ATM, kg/m3
C_CH4 <- 12.0107/16.042
C_CO2 <- 12.0107/44.009

dat_biogas <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info")) %>% 
  filter(day >= 283) %>% select(reactor, gas, temp, day, dm, vs, bio_wet_weight, GD_weight_before, GD_weight_after, mass_vented, vol_vented, room_temp, temp_HS, p_amb) %>%
  mutate(time = day - 283)
vs_mass <- dat_biogas %>% filter(time == 0) %>% mutate(vs_mass_load = dm/100 * vs * bio_wet_weight/1000) %>% select(vs_mass_load)
dat_biogas$vs_mass_load <- vs_mass$vs_mass_load

dat <- as.data.frame(dat_biogas) %>% select(-vs, -dm)

GDout <- calcBgGD(dat, temp.vol = 22, temp.grav = 34, pres.vol = 'p_amb', pres.grav = 1500, 
         id.name = 'reactor', vol.name = 'vol_vented', m.pre.name = 'GD_weight_before', m.post.name = 'GD_weight_after', 
         time.name = 'time', comp.name = 'xCH4', 
         vented.mass = TRUE, averaging = 'final', vmethod = 'vol',
         extrap = TRUE,
         addt0 = TRUE, showt0 = TRUE, comp.lim = c(0.05, 0.7), comp.sub = 'lim',
         unit.pres = 'mbar')

GDout$mCH4_C <- GDout$rvCH4 * rho_CH4/1000 * C_CH4/GDout$vs_mass_load
GDout$mCO2_C <- (GDout$rvBg - GDout$rvCH4) * rho_CO2/1000 * C_CO2/GDout$vs_mass_load
GDout$mCH4_C_cum <- GDout$cvCH4 * rho_CH4/1000 * C_CH4/GDout$vs_mass_load
GDout$mCO2_C_cum <- (GDout$cvBg - GDout$cvCH4) * rho_CO2/1000 * C_CO2/GDout$vs_mass_load