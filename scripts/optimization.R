# Author: Frederik Dalby, 21-06-2023

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(marelac)
library(ReacTran)
library(readxl)

source('cal_functions.R')

dat <- read_excel('../data/O2_profiles_corr.xlsx')

input_dat <- dat %>% mutate(id = as.factor(id)) %>% filter(exclude == 0)

model_dat <- NULL

for(i in unique(input_dat$id)){

dat <- input_dat[input_dat$id == i,] %>% filter(umol > 0)
id <- i
manure <- unique(dat$man_name)
time <- unique(dat$time)
parms <- c(D = 2.2*10^-5, kH = 0.0013/1000*10^6, Cair = max(dat$umol)/1000/1.3)
pars.cal <- log10(data.frame(r = 0.1, KO2 = 0.006))

  # run calibration function
  cal <- optim(par = pars.cal, 
               fn = resCalc,
               dat = dat,
               to = c('umol'), 
               parms = parms,
               method = "L-BFGS-B",
               lower = log10(c(0.0001, 0.0001)),
               upper = log10(c(10, 0.1)),
               control = list(reltol = 0.01),
               hessian = TRUE
  )
  
r <- 10^cal$par

model_dat1 <- c(r , id = id, slurry = manure, time = time, Cair = max(dat$umol)/1000/1.3)
model_dat <- rbind(model_dat, model_dat1)

}

model_dat <- as.data.frame(model_dat) %>% arrange(slurry, time)

write.csv(model_dat, '../output/optim_pars.csv', row.names = F)
