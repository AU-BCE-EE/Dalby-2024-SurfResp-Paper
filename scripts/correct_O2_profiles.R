# Author: Frederik Dalby, 21-06-2023

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(marelac)
library(ReacTran)
library(readxl)
library(openxlsx)

dat <- read_excel('../data/O2_profiles.xlsx')

dat.m <- dat %>% group_by(manure, time) %>% filter(mV > 0) %>%
  mutate(max = max(mV, na.rm = T), min = min(mV, na.rm = T)) %>% 
         mutate(slope = (273)/(max-min)) %>% 
  mutate(int = 273 - (slope * max)) %>%
  mutate(umol_corr = mV * slope + int)

dat$max <- 0
dat$min <- 0
dat$slope <- 0
dat$int <- 0
dat$umol_corr <- dat$umol

dat <- dat[dat$manure %in% c('man4', 'man3', 'man_old'), ]

dat.corr <- rbind(dat, dat.m) %>% mutate(umol = umol_corr) %>% select(-max, -min, -int, -slope, -umol_corr) %>% 
  filter(exclude == 0)

dat.corr$man_name[dat.corr$manure == 'man6'] <- 'Pig 3'
dat.corr$man_name[dat.corr$manure == 'man7'] <- 'Cattle 2'
dat.corr$man_name[dat.corr$manure == 'man8'] <- 'Pig 4'
dat.corr$man_name[dat.corr$manure == 'man9'] <- 'Pig 5'
dat.corr$man_name[dat.corr$manure == 'man10'] <- 'Cattle 3'


write.xlsx(dat.corr,'../data/O2_profiles_corr.xlsx')
