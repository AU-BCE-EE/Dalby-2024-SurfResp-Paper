# isotope ratios

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

CO2_bg <- 430
CH4_bg <- 2 

dCO2_bg <- -8
dCH4_bg <- -47.2

ifelse(dat.org$gas == "air", dat.org$CO2_emis <- (as.numeric(dat.org$co2) - CO2_bg) , dat.org$CO2_emis <- dat.org$co2)  
ifelse(dat.org$gas == "air", dat.org$CH4_emis <- (as.numeric(dat.org$ch4) - CH4_bg) , dat.org$CH4_emis <- dat.org$ch4)  

dat.mod <- dat.org %>% mutate(delta_ch4_corr = (delta_ch4 * (CH4_emis + CH4_bg) -  (dCH4_bg * CH4_bg))/CH4_emis, 
                              delta_co2_corr =  (delta_co2 * (CO2_emis + CO2_bg) -  (dCO2_bg * CO2_bg))/CO2_emis) %>% group_by(temp, gas, day, date) %>% 
  summarise(across(c('delta_ch4','delta_co2'), .fns = list(mean = mean, sd = sd), na.rm = TRUE)) 

dat.mod.long <- dat.mod %>% 
pivot_longer(cols = c('delta_ch4_mean', 'delta_co2_mean', 'delta_ch4_sd', 'delta_co2_sd'), names_to = 'comp', values_to = 'value')

dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 

dat.mod.both$temp <- as.factor(dat.mod.both$temp)
dat.mod.both <- dat.mod.both %>% mutate(group = ifelse(gas == 'n2' & temp == 10, "10 \u00b0C, N2", 
                                                       ifelse(gas == 'n2' & temp == 20, "20 \u00b0C, N2",
                                                              ifelse(gas == 'air' & temp == 10, "10 \u00b0C, Air", "20 \u00b0C, Air"))))
dat.mod.both$group <- as.factor(dat.mod.both$group)

# plot emission rates
new.lab <- as_labeller(c(delta_ch4_mean = "CH[4]", delta_co2_mean = "CO[2]"), label_parsed)

fig_delta <- ggplot(dat.mod.both, aes(x = day, y = value, col = group)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(~comp, labeller = new.lab, scales = "free_y") + 
  scale_x_continuous(breaks = seq(0, 250, by = 50)) +
  theme_bw() +
  theme(text = element_text(size = 10)) + 
  theme(legend.position = 'top') +
  labs(y = expression(delta^{13}*'C (\u2030)'), x = "", col = expression('Temp (\u00b0C)'))

# save figures
png('../figures/fig_delta.png',  width = 18/2.54, height = 8/2.54, units = 'in', res = 600)
grid::grid.draw(ggplotGrob(fig_delta))
dev.off()

pdf('../figures/fig_delta.pdf',  width = 18/2.54, height = 8/2.54)
grid::grid.draw(ggplotGrob(fig_delta))
dev.off()
