rm(list = ls())

library('ggplot2')
library('tidyr')
library('dplyr')
library(readxl)
library("gridExtra")

# fig_loss
dat <- read_excel('../output/table1.xlsx')

dat.loss <- dat %>% group_by(temp, gas) %>% 
  summarise(across(contains('mean'), function(.) ./.[day == 0] * 100),
            day  = day) %>% 
  filter(day != 0) 

dat.sd <- dat %>% pivot_longer(cols = 4:ncol(.), names_to = 'comp', values_to = 'value') %>% 
  summarise(values_sd = value[grepl('sd', comp)], values_mean = value[grepl('mean', comp)],
            names = comp[grepl('sd', comp)], temp = temp[grepl('sd', comp)], gas = gas[grepl('sd', comp)],
            day = day[grepl('sd', comp)])
dat.sd$names <- gsub('_sd','',dat.sd$names)

dat.rsd <- dat.sd %>% group_by(temp, gas, names) %>% 
  mutate(rsd = sqrt((values_sd/values_mean*100)^2 + (values_sd[day == 0]/values_mean[day == 0]*100)^2)) %>% 
  filter(day != 0) %>% mutate(abssd = rsd/100 * values_mean)

dat.rsd_wide <- dat.rsd %>% select(names, abssd, temp, gas, day) %>% 
  pivot_wider(names_from = names, values_from = abssd)

joined <- full_join(dat.loss, dat.rsd_wide, join_by('temp','gas','day'))

joined.long.mean <- joined %>% select(day, temp, gas, contains('mean')) %>% 
  pivot_longer(cols = contains('mean'), names_to = 'names', values_to = 'values_mean')
joined.long.mean$names <- gsub('_mean','',joined.long.mean$names)

joined.long.sd <- joined %>% select(day, temp, gas, !contains('mean')) %>% 
  pivot_longer(cols = -c(temp, day, gas), names_to = 'names', values_to = 'values_sd') 

joined.long <- full_join(joined.long.mean, joined.long.sd, join_by('day', 'temp', 'gas','names'))

plot_dat <- joined.long %>% filter(names %in% c('vs','cel','hem','lig','CP','lip'))
plot_dat$names <- factor(plot_dat$names, levels = c('vs','cel','hem','lig','CP','lip'))

plot_dat$day <- round(plot_dat$day, 0)
plot_dat$day <- gsub('283','After storage',plot_dat$day)
plot_dat$day <- gsub('455','After AD',plot_dat$day)
plot_dat$day <- factor(plot_dat$day, levels = c('After storage', 'After AD'))

plot_dat$temp <- as.factor(plot_dat$temp)

new.lab <- as_labeller(c('vs' = 'VS', 'cel' = 'Cellulose', 'hem' = 'Hemicellulose', 'lig' = 'Lignin', 'CP' = 'CP', 'lip' = 'Lipids',
                         'air' = 'Air', 'n2' = 'N[2]'), label_parsed)

p_loss <- ggplot(plot_dat, aes(x = as.numeric(day), y = values_mean, col = temp, shape = gas)) + 
  geom_point(position = position_dodge(width = 0.2), size = 3) + geom_line(position = position_dodge(width = 0.2))  +
  geom_linerange(aes(y = values_mean, ymin = values_mean - values_sd, ymax = values_mean + values_sd), position = position_dodge(width = 0.2)) +
  labs(y = 'Remaining mass, %', x = "", col = 'Temp, \u00B0 C', shape = 'Headspace gas', tag = 'a') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60,  hjust = 1, vjust = 1)) + facet_grid(~names, labeller = new.lab) +
  scale_color_manual(values = c("blue", "red")) + scale_shape_manual(values = c(1, 2), labels = c("Air", expression('N'[2]))) +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("After storage", 'After AD')) + 
  geom_hline(yintercept = 100, color = 'gray', lty = 'dashed')

#png('../figures/fig_loss.png', height = 4, width =7, units = 'in', res = 600)
#grid::grid.draw(p_loss)
#dev.off()

dat.org <- data.frame(read_excel("../data/dat_resp.xlsx", sheet = "vfa"))
temp <- c(rep('start', 3), rep(c(10, 20, 10, 20), each = 4, times = 2))
gas <- c(rep('start', 3), rep(c('n2', 'air'), each = 8, times = 2))

dat.org$temp <- temp
dat.org$gas <- gas
dat.org$days 

dat.plot <- dat.org %>% select(-sum, -...11, -reactor) %>% 
  pivot_longer(cols = c(1:8), names_to = 'vfa', values_to = 'value') %>% group_by(temp, gas, vfa, days) %>% 
  summarise(mean = mean(value)/1000, sd = sd(value)/1000)

start.dat <- do.call(rbind, replicate(8, dat.plot[dat.plot$temp == 'start',], simplify = F))
start.dat$temp <- c(rep(10, 16), rep(20, 16), rep(10, 16), rep(20, 16))
start.dat$gas <- c(rep(c('air','n2'), each = 8, times = 2), rep(c('air','n2'), each = 8, times = 2))


dat.plot.mod <- dat.plot[!grepl('start', dat.plot$gas),]
dat.plot.mod$temp <- as.numeric(dat.plot.mod$temp)
dat.plot.all <- rbind(dat.plot.mod, start.dat)

new.lab = as_labeller(c("10"= "10~degree~C", "n2"= "N[2]", "20"= "20~degree~C", "air"= "Air", 'vs' = 'VS', 'cel' = 'Cellulose', 'hem' = 'Hemicellulose', 'lig' = 'Lignin', 'CP' = 'CP', 'lip' = 'Lipids'
                        ), label_parsed)

dat.plot.all$days <- as.factor(dat.plot.all$days )
dat.plot.all$days <- recode(dat.plot.all$days, '0' = "Start", '283' = 'After storage', '455' = 'After AD')

dat.plot.all$vfa <- factor(dat.plot.all$vfa, levels = c('acetic' = 'acetic', 'propanoic'= 'propanoic',
                                                        'butanoic' = 'butanoic', 'iso.butanoic' = 'iso.butanoic',
                                                        'meth.butanoic'='meth.butyric', 'pentanoic' = 'pentanoic',
                                                        'hexanoic' = 'hexanoic','iso.caproic' = 'iso.caproic'))

p_vfa <- 
  ggplot(dat.plot.all, aes(vfa, mean, fill = days)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.8), width = 0.2) +
  facet_grid(temp~gas, labeller = new.lab) + theme_bw() +
  labs(y = expression('VFA (g kg'^{-1}*')'), x = '', fill = '', tag = 'b') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.0, hjust = 1))

comb <- grid.arrange(p_loss, p_vfa, nrow = 2)

png('../figures/fig_loss_vfa.png',  width = 7, height = 8, units = 'in', res = 600)
grid::grid.draw(comb)
dev.off()
