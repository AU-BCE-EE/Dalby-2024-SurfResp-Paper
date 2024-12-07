rm(list = ls())

library('dplyr')
library('tidyr')
library('broom')

storage_dat <- read.csv('../output/emis_stat_storage.csv')
biogas_dat <- read.csv('../output/emis_stat_biogas.csv')

combined <- rbind(storage_dat, biogas_dat) %>% 
  group_by(reactor, comp, temp, gas) %>% 
  summarise(cum = sum(cum)) %>% 
  mutate(incubation = 'combined')

dat_stat <- rbind(storage_dat, biogas_dat, combined) %>% mutate(temp = factor(temp))

# statistics
model1 <- dat_stat %>% group_by(comp, incubation) %>% 
  do({
    fit <- aov(cum ~ gas * temp, data = .)
    posthoc <- TukeyHSD(fit)
    bind_rows(tidy(fit), tidy(posthoc))
  })

model2 <- dat_stat %>% group_by(comp, incubation) %>% 
  do({
    fit <- aov(cum ~ gas + temp, data = .)
    bind_rows(tidy(fit))
  })



# output for table 1
table1 <- dat_stat %>% group_by(gas, temp, comp, incubation) %>% 
  summarise(mean = mean(cum), std = sd(cum))
