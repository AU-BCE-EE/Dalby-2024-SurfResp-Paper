# surface respiration estimate
rm(list = ls())

library("readxl")
library('ggplot2')
library('tidyr')
library('dplyr')

dat.org <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "all"))


f_CH4 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 16 * 1000 # mg/day
f_CO2 <- 1/1000000 * dat.org$cor_flow/(0.082057 * 293) * 60 * 24 * 44 * 1000 # mg/day
  
CO2_bg <- 430
CH4_bg <- 2 
  
ifelse(dat.org$gas == "air", dat.org$CO2_emis <- (as.numeric(dat.org$co2) - CO2_bg) * f_CO2, dat.org$CO2_emis <- dat.org$co2 * f_CO2)  
ifelse(dat.org$gas == "air", dat.org$CH4_emis <- (as.numeric(dat.org$ch4) - CH4_bg) * f_CH4, dat.org$CH4_emis <- dat.org$ch4 * f_CH4)  

dat.mod <- dat.org %>% filter(reactor != 'bg') %>% group_by(temp, gas, day, date) %>% mutate(ratio = (CH4_emis * 12/16)/(CH4_emis * 12/16 + CO2_emis * 12/44)) %>% 
  group_by(temp, gas, day, date) %>% summarise(across(c('CO2_emis','CH4_emis','ratio'), .fns = list(mean = mean, sd = sd), na.rm = TRUE))

dat.mod.long <- dat.mod %>% pivot_longer(cols = c('CH4_emis_mean', 'CO2_emis_mean', 'CH4_emis_sd', 'CO2_emis_sd'), names_to = 'comp', values_to = 'value')

dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 

dat.mod.both$temp <- as.factor(dat.mod.both$temp)

write.csv(dat.mod.both, '../output/respiration_emis.csv', row.names = F)

# plot emission rates

new.lab <- as_labeller(c(air = "Air", n2 = "N[2]", CH4_emis_mean = "CH[4]", CO2_emis_mean = "CO[2]"), label_parsed)

fig_emis <- ggplot(dat.mod.both, aes(x = day, y = value, col = temp)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(comp~gas, labeller = new.lab, scales = "free_y") + 
  theme_bw() +
  theme(text = element_text(size = 16)) + 
  labs(y = expression('emission rate (mg day'^{-1}*')'), x = "", col = expression('Temp (\u00b0C)'))

## calc cumulative emission
day_spaced <- seq(from = min(dat.mod.both$date), to = max(dat.mod.both$date), 
                  length.out = max(dat.mod.both$date) - min(dat.mod.both$date) +1)

end <- max(dat.org$day)

dat_emis_table_mean <- dat.org %>% filter(reactor != 'bg') %>% 
  group_by(reactor) %>% 
  mutate(ratio = (CH4_emis * 12/16)/(CH4_emis * 12/16 + CO2_emis * 12/44)) %>% 
  mutate(cum_CH4 = sum(approx(x = date, y = CH4_emis, xout = day_spaced)$y),
         cum_CO2 = sum(approx(x = date, y = CO2_emis, xout = day_spaced)$y)) %>% 
  group_by(temp, gas) %>%
  summarise(cum_CH4 = mean(cum_CH4[day == end])/1000, 
            cum_CO2 = mean(cum_CO2[day == end])/1000
            )

dat_emis_table_sd <- dat.org %>% filter(reactor != 'bg') %>% 
    group_by(reactor) %>% 
    mutate(ratio = (CH4_emis * 12/16)/(CH4_emis * 12/16 + CO2_emis * 12/44)) %>% 
    mutate(cum_CH4 = sum(approx(x = date, y = CH4_emis, xout = day_spaced)$y),
           cum_CO2 = sum(approx(x = date, y = CO2_emis, xout = day_spaced)$y)) %>% 
  group_by(temp, gas) %>%
    summarise(sd_CH4 = sd(cum_CH4[day == end])/1000, 
              sd_CO2 = sd(cum_CO2[day == end])/1000
    )
  
stat <- cbind(dat_emis_table_mean, dat_emis_table_sd) 
write.csv(stat, '../output/emis_table2_1.csv')
  
dat.mod.both_cum <- dat.mod.both %>% group_by(temp, gas, comp) %>% 
  mutate(cum = sum(approx(x = date, y = value, xout = day_spaced)$y)) #%>% 
  distinct(cum, keep_all = TRUE) %>% mutate(cum = cum/1000) %>% ungroup() %>%
  filter(comp != 'ratio_mean')

dat.mod.long <- dat.mod %>% pivot_longer(cols = c('ratio_mean', 'ratio_sd'), names_to = 'comp', values_to = 'value')

dat.mod.mean <- dat.mod.long[grepl('mean$', dat.mod.long$comp),]
dat.mod.sd <- dat.mod.long[grepl('sd$', dat.mod.long$comp),]
dat.mod.both <- cbind(dat.mod.mean, sd = dat.mod.sd$value) 

new.lab <- as_labeller(c(air = "Air", n2 = "N[2]"), label_parsed)

fig_ratio <- ggplot(dat.mod.both, aes(x = day, y = value, col = as.factor(temp))) + 
  geom_point() + 
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd)) + 
  facet_grid(~gas, labeller = new.lab, scales = "free_y") + 
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(y = expression('C-CH'[4]*'/ (C-CH'[4]*' + C-CO'[2]*')'), x = "", col = expression('Temp (\u00b0C)'))

p_all <- grid.arrange(fig_emis, fig_ratio, layout_matrix = rbind(c(1), c(2)))

png('../figures/fig_emis.png',  width = 18/2.54, height = 18/2.54, units = 'in', res = 600)
grid::grid.draw(p_all)
dev.off()

pdf('../figures/fig_emis.pdf',  width = 18/2.54, height = 18/2.54)
grid::grid.draw(p_all)
dev.off()