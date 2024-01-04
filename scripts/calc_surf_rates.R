rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')
library('gridExtra')

#prep
dat.org <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "all"))
dat.info <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "info"))

dat.info.f <- dat.info %>% select(c('reactor', 'temp','gas','datetime','day', 'wet_weight', 'dm','vs','lig','cel','hem','lip','tan','tn','vfa','pH')) %>%
  mutate(dm = dm * 10, vs = vs * dm)
mass <- dat.info.f$wet_weight[dat.info.f$day == 0]/1000 # kg
mass <- mass[!is.na(mass)]

dat.org <- dat.org[dat.org$reactor != 'bg',]
dat.org$mass <- mass

f_CH4 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 16 *1000 / dat.org$mass # mg pr. kg slurry pr day
f_CO2 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 44 *1000 / dat.org$mass # mg pr. kg slurry pr day

CO2_bg <- 430
CH4_bg <- 2 

surf_area <- (9.5/2)^2 * pi / 10000# m^2 

ifelse(dat.org$gas == "air", dat.org$CO2_emis <- (as.numeric(dat.org$co2) - CO2_bg) * f_CO2, dat.org$CO2_emis <- dat.org$co2 * f_CO2)  
ifelse(dat.org$gas == "air", dat.org$CH4_emis <- (as.numeric(dat.org$ch4) - CH4_bg) * f_CH4, dat.org$CH4_emis <- dat.org$ch4 * f_CH4)  

surf_conv <- mean(mass/surf_area) # convert mg / kg slurry / day to mg / m2 / day

dat.mod <- dat.org %>% group_by(reactor, temp, gas, day, date) %>% 
  mutate(ratio = (CH4_emis * 12.01/16.04)/(CH4_emis * 12.01/16.04 + CO2_emis * 12.01/44.01)) %>% 
  group_by(temp, gas, day, date) %>% summarise(across(c('CO2_emis','CH4_emis','ratio'), .fns = list(mean = mean, sd = sd), na.rm = TRUE)) 
dat.mod$CO2_surf_emis_mean <- 0
dat.mod$CO2_surf_emis_mean[dat.mod$temp == 10 & dat.mod$gas == 'air'] <- dat.mod$CO2_emis_mean[dat.mod$temp == 10 & dat.mod$gas == 'air'] - dat.mod$CO2_emis_mean[dat.mod$temp == 10 & dat.mod$gas == 'n2']
dat.mod$CO2_surf_emis_mean[dat.mod$temp == 20 & dat.mod$gas == 'air'] <- dat.mod$CO2_emis_mean[dat.mod$temp == 20 & dat.mod$gas == 'air'] - dat.mod$CO2_emis_mean[dat.mod$temp == 20 & dat.mod$gas == 'n2']
dat.mod$CO2_surf_emis_sd <- 0
dat.mod$CO2_surf_emis_sd[dat.mod$temp == 10 & dat.mod$gas == 'air'] <- sqrt((dat.mod$CO2_emis_sd[dat.mod$temp == 10 & dat.mod$gas == 'air'])^2 - (dat.mod$CO2_emis_sd[dat.mod$temp == 10 & dat.mod$gas == 'n2'])^2)
dat.mod$CO2_surf_emis_sd[dat.mod$temp == 20 & dat.mod$gas == 'air'] <- sqrt((dat.mod$CO2_emis_sd[dat.mod$temp == 20 & dat.mod$gas == 'air'])^2 + (dat.mod$CO2_emis_sd[dat.mod$temp == 20 & dat.mod$gas == 'n2'])^2)

dat.mod$CO2_surf_emis_mean <- dat.mod$CO2_surf_emis_mean /1000 * surf_conv
dat.mod$CO2_surf_emis_sd <- dat.mod$CO2_surf_emis_sd /1000 * surf_conv

# CO2_surf_emis has units of g CO2/m2/d. Omit first 12 days, to ensure stable bacteria in surface. 
# Omit above 100 days, to ensure CO2 resp has not depleted carbon so much that we cannot compare 
# air and n2 reactors anymore and assume that the difference is only due to surface resp.
# this is also where we see the change in delta values. 

dat.mod.surf <- dat.mod %>% filter(day < 100) %>% group_by(temp, gas) %>% 
 summarize(m_surf = mean(CO2_surf_emis_mean), sd_rel_surf = sqrt(sum((CO2_surf_emis_sd/CO2_surf_emis_mean*100)^2))/100) %>%
  mutate(abs_sd = sd_rel_surf * m_surf)


# calculate kl_O2

temp_standard <- 298
temp_K <- dat.mod.surf$temp[dat.mod.surf$m_surf != 0] + 273.15

kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
E_surf <- dat.mod.surf$m_surf[dat.mod.surf$m_surf != 0]
area <- surf_area
temp_C <- c(10, 20)
pars <- data.frame(kH_oxygen = kH_oxygen, E_surf = E_surf, area = area, temp_C = temp_C) 
CO2_aer <- 1.12 # amount of CO2 produced pr gCOD consumed by aerobic bacteria. if there was no bacteria growth it would be 44gCO2/32gO2
pars$kl.oxygen <- pars$E_surf/CO2_aer / (pars$kH_oxygen * 0.208)

# Step 3: Fit the linear regression model
model <- lm(log(kl.oxygen) ~ temp_C, data = pars)

plot(pars$temp_C, pars$kl.oxygen, ylim = c(0,7), 
     xlim = c(5, 25), ylab = expression('kl'[O2]*'(m day'^{-1}*')'), xlab = 'Temperature (\u00B0 C)')
temps <- c(5:25)
pred <- exp(model$coefficients[1] + model$coefficients[2] * temps)

pred_klo2 <- data.frame(cbind(temp = as.numeric(temps), klo2 = as.numeric(pred), type = 'model'))

plot_klo2 <- data.frame(rbind(pred_klo2, data.frame(temp = pars$temp_C, klo2 = pars$kl.oxygen, type = 'data')))

ggplot(plot_klo2, aes(as.numeric(temp), as.numeric(klo2))) +
  geom_line() +
  geom_point(data = subset(plot_klo2, type == 'data'), aes(as.numeric(temp), as.numeric(klo2)), col = "red", size = 1) + theme_bw() +
  labs(y = expression('kl'[O2]*'(m day'^{-1}*')'), x = 'Temperature (\u00B0 C)')
ggsave('../figures/fig_klo2.png', height = 3, width = 3)
write.csv(pars, '../output/surf_resp_dat.csv', row.names = F)
