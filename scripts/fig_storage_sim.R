rm(list=ls())

library(ABM)
library(dplyr)
library(tidyr)
library(ggplot2)

source('../R/slurry_app.R')
source('../R/doy.R')

# barn inputs: 
wthr_pars_barn = list(temp_air_C = 20, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars_barn = list(evap = 0.5 * ABM:::et(temp_C = wthr_pars_barn$temp_air_C, pres_kpa = wthr_pars_barn$pres_kpa, rs = wthr_pars_barn$rs))
mng_pars_barn = list(slurry_prod_rate = 2000 * 5.7,  
                 slurry_mass = 0,          
                 storage_depth = 0.6,        
                 resid_depth = 0.05,         
                 floor_area = 1300,
                 area = 1430, 
                 empty_int = 28,
                 temp_C = 18.6,
                 wash_water = 75 * 2000,
                 wash_int = 84,
                 rest_d = 5,
                 RH = 90, 
                 cover = NA,
                 resid_enrich = 0.9,
                 slopes = c(urea = NA, slurry_prod_rate = NA),
                 graze = c(start = 'May', duration = 0, hours_day = 0),
                 scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1))

man_pars_barn = list(conc_fresh = list(sulfide = 0.01, urea = 3.45, sulfate = 0.2, TAN = 0, starch = 5.25, 
                                  VFA = 1.71, xa_dead = 0, Cfat = 27.58, CP = 21.13, RFd = 25.43, iNDF = 11.31, 
                                  VSd = 0, VSd_A = 54.78, VSnd_A = 23.48, ash = 15), pH = 7, dens = 1000)

days <- 3*365

# call abm for barn simulation
barn_dat <- cbind(abm(days, 1, wthr_pars = wthr_pars_barn, evap_pars = evap_pars_barn, 
    mng_pars = mng_pars_barn, man_pars = man_pars_barn), source = 'barn')

# storage inputs
wthr_pars_storage = list(temp_air_C = 15, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars_storage = list(evap = 0.5 * ABM:::et(temp_C = wthr_pars_storage$temp_air_C, pres_kpa = wthr_pars_storage$pres_kpa, rs = wthr_pars_storage$rs))

# get fresh concentrations as the effluent from the barn (averages)
conc_fresh <- barn_dat[barn_dat$slurry_mass_eff != 0, grepl('eff_conc', names(barn_dat))] %>% 
  summarise(across(everything(), mean)) %>% rename_with(~gsub('_eff_conc', '', .), everything())

xa_fresh <- barn_dat[barn_dat$slurry_mass_eff != 0, grepl('eff_conc', names(barn_dat))] %>% 
  select(starts_with('xa') & -contains('dead')) %>% 
  summarise(across(everything(), mean)) %>% rename_with(~gsub('_eff_conc', '', .), everything()) %>%
  rename_with(~gsub('xa_', '', .), everything()) 

conc_fresh_storage  <- as.list(conc_fresh[names(man_pars_barn$conc_fresh)])
xa_fresh_storage <- setNames(as.numeric(xa_fresh), names(xa_fresh))

# get slurry mass effluents from barn
slurry_cum <- cumsum(barn_dat$slurry_mass_eff[!duplicated(barn_dat$time)])
time <- seq(from = 0, to = max(barn_dat$time), 1)
slurry_mass_dat <- data.frame(time = time, slurry_mass = slurry_cum)

# slurry application pattern
specs <- data.frame(app_t1s = "March", app_t2s = "April", app_t3s = "May", app_t4s = "June",
                    app_t5s = "July", app_t6s = "August", app_t7s = "September", 
                    app1s = 0.3, app2s = 0.89, app3s = 0.65, app4s = 0.07, 
                    app5s = 0.05, app6s = 0.11, app7s = 0.12)

# adjust for slurry application  
slurry_mass_dat <- slurry_app(days = days, begin = 'January', specs = specs, from = 'storage', slurry = slurry_mass_dat)

# get temp dat from Vechi 2023
temp_dat <- data.frame(time = c(outside_slurry_temp_vechi$time, outside_slurry_temp_vechi$time+365,
                                 outside_slurry_temp_vechi$time + 365+365), temp_C = rep(outside_slurry_temp_vechi$temp_C, 3))

temp_dat <- temp_dat[!duplicated(temp_dat),]

# combine inputs for storage
mng_pars_storage = list(slurry_prod_rate = 0,  
                     slurry_mass = slurry_mass_dat,          
                     storage_depth = 4.6,        
                     resid_depth = 0,         
                     floor_area = 0,
                     area = 1000, 
                     empty_int = 28,
                     temp_C = temp_dat,
                     wash_water = NA,
                     wash_int = 84,
                     rest_d = 5,
                     RH = 90, 
                     cover = 'tent',
                     resid_enrich = 0.9,
                     slopes = c(urea = NA, slurry_prod_rate = NA),
                     graze = c(start = 'May', duration = 0, hours_day = 0),
                     scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1))

man_pars_storage = list(conc_fresh = conc_fresh_storage, pH = 7.3, dens = 1000)

storage_dat <- cbind(abm(days, 1, mng_pars = mng_pars_storage, man_pars = man_pars_storage,
                         add_pars = list(xa_fresh = xa_fresh_storage)), source = 'storage')

# prepare for figure 5
# remove slurry rem peaks
rem_times <- storage_dat$time[storage_dat$CH4_emis_rate>1.5e05]

dat <- bind_rows(barn_dat, storage_dat) %>% filter(!time %in% rem_times)


dat$n_anim = 2000

dat$mol.CH4 <- dat$CH4_emis_rate * 12.01/16.04
dat$mol.CO2 <- dat$CO2_emis_rate * 12.01/44.01
dat$mol.ratio <- dat$mol.CH4/(dat$mol.CH4 + dat$mol.CO2)

dat.plot_CO2 <- dat %>% filter(time > 2*365) %>% mutate(time = time - 2*365, 
                                                        CH4_emis_rate = CH4_emis_rate * 12/16/n_anim, 
                                                        CO2_emis_rate = CO2_emis_rate * 12/44/n_anim,
                                                        respiration = respiration * 0.436 * 12/44/n_anim,
                                                        ureolysis = rut_urea * 1.57*12/44/n_anim,
                                                        fermentation = CO2_ferm_CO2*12/44/n_anim,
                                                        methanogenesis = (CO2_ferm_meth_sr_CO2 - CO2_ferm_CO2)*12/44/n_anim,
                                                        slurry_mass = slurry_mass / 1000
) %>%
  pivot_longer(c('CO2_emis_rate', 'respiration', 'ureolysis', 'fermentation', 'methanogenesis','CH4_emis_rate', 'mol.ratio', 'temp_C', 'slurry_mass'), names_to = "comp", values_to = "value")

source.lab <- as_labeller(c(barn = "Barn", storage = "Storage", 
                            "1" = "CO[2]~(g~C~pig^{-1}~d^{-1})",
                            "2" = "CH[4]~(g~C~pig^{-1}~d^{-1})",
                            "3" = "C[CH4]/(C[CH4]+C[CO2])",
                            "4" = "Temperature~(degree*C)",
                            "5" = 'Slurry~volume~(m^{3})'), label_parsed)

custom_colors <- c(
  "CO2_emis_rate" = "gray",
  "respiration" = "red",
  "ureolysis" = "palegreen4",
  "fermentation" = "purple",
  "methanogenesis" = "orange"
)

custom_labels <- c(
  "CO2_emis_rate" = expression(paste("CO"[2], " total")),
  "respiration" = expression(paste("CO"[2], " respiration")),
  "ureolysis" = expression(paste("CO"[2], " ureolysis")),
  "fermentation" = expression(paste("CO"[2], " fermentation")),
  "methanogenesis" = expression(paste("CO"[2], " ", methanogenesis + sulfate~reduction))
)


dat.plot_CO2$plot <- 1
dat.plot_CO2$plot[dat.plot_CO2$comp == "CH4_emis_rate"] <- 2
dat.plot_CO2$plot[dat.plot_CO2$comp == "mol.ratio"] <- 3
dat.plot_CO2$plot[dat.plot_CO2$comp == "temp_C"] <- 4
dat.plot_CO2$plot[dat.plot_CO2$comp == "slurry_mass"] <- 5

dat.plot_CO2$plot <- as.factor(dat.plot_CO2$plot)




# Plot
CO2_plot_full <- ggplot(dat.plot_CO2, aes(x = time, y = value, col = comp)) + 
  geom_line(size = 1) +
  theme_bw() +
  theme(text = element_text(size = 11), axis.text.x = element_text(angle = 60, hjust = 1),
        legend.text.align = 0, legend.position = "top" ) +
  scale_x_continuous(
    breaks = c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336),
    labels = month.abb[c(3:12, 1, 2)]
  ) +
  labs(y = "", x = "", col = "") +
  facet_grid(plot~source, labeller = source.lab, scales = "free_y") + 
  scale_color_manual(
    values = custom_colors,
    labels = custom_labels
  ) + 
  guides(col=guide_legend(nrow=2, byrow=TRUE))


png('../figures/fig_storage_sim.png',  width = 15/2.54, height = 20/2.54, units = 'in', res = 600)
grid::grid.draw(CO2_plot_full)
dev.off()

dat.plot_red <- dat.plot_CO2[dat.plot_CO2$plot %in% c(1,2,3),]

CO2_plot_paper <- ggplot(dat.plot_red, aes(x = time, y = value, col = comp)) + 
  geom_line(size = 1) +
  theme_bw() +
  theme(text = element_text(size = 11), axis.text.x = element_text(angle = 60, hjust = 1),
        legend.text.align = 0, legend.position = "top" ) +
  scale_x_continuous(
    breaks = c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336),
    labels = month.abb[c(3:12, 1, 2)]
  ) +
  labs(y = "", x = "", col = "") +
  facet_grid(plot~source, labeller = source.lab, scales = "free_y") + 
  scale_color_manual(
    values = custom_colors,
    labels = custom_labels
  ) + 
  guides(col=guide_legend(nrow=2, byrow=TRUE))

png('../figures/fig_storage_sim_paper.png',  width = 6, height = 6, units = 'in', res = 600)
grid::grid.draw(CO2_plot_paper)
dev.off()

pdf('../figures/fig_storage_sim.pdf',  width = 15/2.54, height = 18/2.54)
grid::grid.draw(CO2_plot_paper)
dev.off()

dat.stat <- dat %>% filter(time > 2*365) %>% mutate(time = time - 2*365) %>% 
  group_by(source) %>% 
  summarise(avg_ratio = mean(mol.ratio, na.rm = T), cum_ratio = sum(mol.CH4, na.rm = T)/(sum(mol.CH4, na.rm = T) + sum(mol.CO2, na.rm = T)),
            CH4_kg_m3_year = mean(CH4_emis_rate/1000, na.rm = T)/mean(slurry_mass/1000, na.rm = T)*365,
            CH4_A_kg_m3_year = mean(CH4_A_emis_rate/1000, na.rm = T)/mean(slurry_mass/1000, na.rm = T)*365,
            CH4_g_pig_day = mean(CH4_emis_rate, na.rm = T)/mean(n_anim),
            CH4_kg_m3_excreted = mean(CH4_emis_rate/1000, na.rm = T)/mean((n_anim * 0.47 * 365 / 89 / 1000)), 
            CO2_kg_m3_year = mean(CO2_emis_rate/1000, na.rm = T)/mean(slurry_mass/1000, na.rm = T)*365,
            CO2_kg_m2_year = mean(CO2_emis_rate/1000, na.rm = T)/area[1] * 365, 
            C_loss_resp = mean(respiration * 0.436 * 12.01/44.01/1000, na.rm = T)/(mean(CO2_emis_rate * 12.01/44.01/1000, na.rm = T) + mean(CH4_emis_rate* 12.01/16.04/1000, na.rm = T)),
            C_loss_urea = mean(rut_urea * 1.57 * 12.01/44.01/1000, na.rm = T)/(mean(CO2_emis_rate * 12.01/44.01/1000, na.rm = T) + mean(CH4_emis_rate* 12.01/16.04/1000, na.rm = T)))

write.csv(dat.stat, '../output/storage_sim_stat.csv', row.names = F)



