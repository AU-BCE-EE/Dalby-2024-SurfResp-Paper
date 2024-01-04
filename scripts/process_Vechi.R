rm(list = ls())

library(readxl)
library(lubridate)
library(dplyr)

#import Vechi 2023 et al. data
dat <- read_excel('../data/Vechi_dat.xlsx')
dat$date <- ymd(dat$date)
dat$month <- month(dat$date)
dat$doy <- as.numeric(as.Date(dat$date) - as.Date("2023-01-01"))
dat$fill <- dat$m3/dat$capacity

dat.mod <- dat %>% filter(grepl('P',Tank)) %>% group_by(month) %>% 
  summarise(doy = mean(doy, na.rm = T), 
            temp = mean(temp, na.rm = T), 
            fill = mean(fill, na.rm = T)
            )

dat.mod$temp <- approx(dat.mod$doy, dat.mod$temp, xout = dat.mod$doy)$y

temp_dat <- data.frame(time = c(dat.mod$doy, dat.mod$doy + 365, dat.mod$doy + 2*365, dat.mod$doy + 3*365), temp_C = rep(dat.mod$temp, 4))
write.csv(temp_dat, "../data/temp_dat.csv", row.names = F)







