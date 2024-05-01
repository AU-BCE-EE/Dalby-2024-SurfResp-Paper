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

dat.surf_resp <- dat.org %>% group_by(reactor, temp, gas, day, date) %>% 
  mutate(ratio = (CH4_emis * 12.01/16.04)/(CH4_emis * 12.01/16.04 + CO2_emis * 12.01/44.01)) %>% filter(day < 100) %>%
  group_by(temp, gas) %>% summarise(across(c('CO2_emis','CH4_emis','ratio'), .fns = list(mean = mean, sd = sd), na.rm = TRUE)) 

surf_resp <- dat.surf_resp %>% mutate(CO2_emis_mean = CO2_emis_mean/ 1000 * surf_conv,
                                CO2_emis_sd = CO2_emis_sd/1000 * surf_conv) %>% ungroup()

means <- surf_resp %>% mutate(CO2_surf_mean_20 = CO2_emis_mean[gas == 'air' & temp == 20] - CO2_emis_mean[gas == 'n2' & temp == 20],
                              CO2_surf_mean_10 = CO2_emis_mean[gas == 'air' & temp == 10] - CO2_emis_mean[gas == 'n2' & temp == 10])

sds <- surf_resp %>% mutate(CO2_emis_sd_20 = sqrt(CO2_emis_sd[gas == 'air' & temp == 20]^2 + CO2_emis_sd[gas == 'n2' & temp == 20]^2),
                            CO2_emis_sd_10 = sqrt(CO2_emis_sd[gas == 'air' & temp == 10]^2 + CO2_emis_sd[gas == 'n2' & temp == 10]^2))

# calculate kl_O2
temp_standard <- 298
temp_K <- c(10+273.15, 20+273.15)

kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
E_surf <- c(unique(means$CO2_surf_mean_10), unique(means$CO2_surf_mean_20))
area <- surf_area
temp_C <- temp_K - 273.15
pars <- data.frame(kH_oxygen = kH_oxygen, E_surf = E_surf, area = area, temp_C = temp_C, temp_K = temp_K) 

# calculate CO2 productivity coeff:
source('stoich.R')

y <- list()
y$CP <- 17.37/0.6541602 # gCOD CP/kg manure
y$starch <- 0 # no starch in the residual pig slurry
y$Cfat <- 10.27/0.3117844 # gCOD Cfat/kg manure
y$RFd <- 51.82 * (25.4/36.7) /0.84447 # g RFd in residual slurry. Calculated from NDF x typical fraction of RFd/NDF divided by gCOD/gRFd
y$VSd <- y$CP + y$Cfat + y$RFd

conc_fresh <- list()
conc_fresh[['VSd']] <- 0

sub_resp <- y$Cfat + y$CP + y$RFd + y$starch 

alpha <- 0 # hydrolysis rate, set to 0 as it does not matter but function requires it. 
pCO2 <- stoich(alpha, y, conc_fresh, sub_resp, respiration = 1)$CO2_resp/1 # respiration is just 1, the size of it does not influence the pCO2 coef (try yourself)
pars$kl.oxygen <- pars$E_surf/pCO2 / (pars$kH_oxygen * 0.208)

#Fit the linear regression model
model <- lm(log(kl.oxygen) ~ temp_C, data = pars)
kl.oxygen <- exp(model$coefficients[1] + model$coefficients[2] * temp_C) # from own lab experiments (Dalby et al. 2023..unpublished) 

plot(pars$temp_C, pars$kl.oxygen, 
     xlim = c(5, 25), ylab = expression('kl'[O2]*'(m day'^{-1}*')'), xlab = 'Temperature (\u00B0 C)')

temps <- c(5:25)

pred <- exp(model$coefficients[1] + model$coefficients[2] * temps)
pred_klo2 <- data.frame(cbind(temp = as.numeric(temps), klo2 = as.numeric(pred), type = 'model'))
plot_klo2 <- data.frame(rbind(pred_klo2, data.frame(temp = pars$temp_C, klo2 = pars$kl.oxygen, type = 'data')))

# uncertainty on kl.o2: 

sds_plot <- c(sds$CO2_emis_sd_10[1], sds$CO2_emis_sd_20[1])/pCO2 / (pars$kH_oxygen * 0.208)

ggplot(plot_klo2, aes(as.numeric(temp), as.numeric(klo2))) +
  geom_line() + 
  geom_linerange(data = subset(plot_klo2, type == 'data'), 
                aes(x = as.numeric(temp), y = as.numeric(klo2), ymin = as.numeric(klo2) - sds_plot, ymax = as.numeric(klo2) + sds_plot)) +
  geom_point(data = subset(plot_klo2, type == 'data'), aes(as.numeric(temp), as.numeric(klo2)), col = "red", size = 1) + theme_bw() +
  labs(y = expression('kl'[O2]*'(m day'^{-1}*')'), x = 'Temperature (\u00B0 C)')
ggsave('../figures/fig_klo2.png', height = 3, width = 3)
write.csv(pars, '../output/surf_resp_dat.csv', row.names = F)
