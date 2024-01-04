validation_script_prime <- function(variant){

path_fun <- '/home/frd/ABM_carbon_accounting/ABM_functions/'
path_dat <- '/home/frd/ABM_carbon_accounting/paper/data/dat_comp.xlsx'

#path_fun <- 'C:/Program Files/GitHub/ABM_carbon_accounting/ABM_functions/'
#path_dat <- 'C:/Program Files/GitHub/ABM_carbon_accounting/paper/data/dat_comp.xlsx'

ff <- list.files(path = path_fun) 
ff <- ff[ff != 'x.R']
for (i in ff) source(paste0(path_fun,i))
  
abm_packages()

# import data
dat <- read_excel(path_dat, sheet = "emis")
dat_ent <- read_excel(path_dat, sheet = "enteric")
dat_temp <- read_excel(path_dat, sheet = "temp")
dat_mass <- read_excel(path_dat, sheet = "mass")
dat_NH3 <- read_excel(path_dat, sheet = "emis NH3")
dat_analysis <- read_excel(path_dat, sheet = "analysis")
dat_weights <- read_excel(path_dat, sheet = "weights")
dat_vfa <- read_excel(path_dat, sheet = "vfa")
dat_odor <- read_excel(path_dat, sheet = "odor")

# add derived output to dat
dat$CH4_rate <- dat$CH4E * dat$pigs
dat$doy <- yday(dat$date)
dat$airflow <- dat$airflow * 24 # from m3/h to m3/day

C <- dat[dat$treatment == "control",]
FF <- dat[dat$treatment == "frequentflushing",]
ST <- dat[dat$treatment == "slurrytrays",]
SF <- dat[dat$treatment == "slurryfunnels",]

# interpolate enteric methane prod, temperature, slurry mass to periods and treatment data
dat_ent_C <- dat_ent[which(dat_ent$treatment == "control"),]
C$enteric <- as.data.frame(approx(y = dat_ent_C$CH4_enterisk * dat_ent_C$No_pigs, x = dat_ent_C$date, xout = C$date))[,2]
dat_temp_C <- dat_temp[which(dat_temp$treatment == "control"), c("date", "temp")]
C <- merge(C, dat_temp_C, all = T)
dat_mass_C <- dat_mass[which(dat_mass$treatment == "control"), c("date", "mass")]
C <- merge(C, dat_mass_C, all = T)
dat_NH3_C <- dat_NH3[which(dat_NH3$treatment == "control"),]
C$NH3 <- as.data.frame(approx(y = dat_NH3_C$NH3, x = dat_NH3_C$date, xout = C$date))[,2]
dat_weights_C <- dat_weights[which(dat_weights$treatment == "control"),]
C$weights <- as.data.frame(approx(y = dat_weights_C$weights, x = dat_weights_C$date, xout = C$date))[,2]
dat_vfa_C <- dat_vfa[which(dat_vfa$treatment == 'control'),]
C <- merge(C, dat_vfa_C, all = T)
dat_H2S_C <- dat_odor[which(dat_odor$treatment == 'Control'),]
C$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_C$M35E, x = dat_H2S_C$date, xout = C$date))[,2] / 24 # g H2S /day

dat_ent_FF <- dat_ent[which(dat_ent$treatment == "frequentflushing"),]
FF$enteric <- as.data.frame(approx(y = dat_ent_FF$CH4_enterisk * dat_ent_FF$No_pigs, x = dat_ent_FF$date, xout = FF$date))[,2]
dat_temp_FF <- dat_temp[which(dat_temp$treatment == "frequentflushing"), c("date", "temp")]
FF <- merge(FF, dat_temp_FF, all = T)
dat_mass_FF <- dat_mass[which(dat_mass$treatment == "frequentflushing"), c("date", "mass")]
FF <- merge(FF, dat_mass_FF, all = T)
dat_NH3_FF <- dat_NH3[which(dat_NH3$treatment == "frequentflushing"),]
FF$NH3 <- as.data.frame(approx(y = dat_NH3_FF$NH3, x = dat_NH3_FF$date, xout = FF$date))[,2]
dat_weights_FF <- dat_weights[which(dat_weights$treatment == "frequentflushing"),]
FF$weights <- as.data.frame(approx(y = dat_weights_FF$weights, x = dat_weights_FF$date, xout = FF$date))[,2]
dat_vfa_FF <- dat_vfa[which(dat_vfa$treatment == 'frequentflushing'),]
FF <- merge(FF, dat_vfa_FF, all = T)
dat_H2S_FF <- dat_odor[which(dat_odor$treatment == 'Frequentflushing'),]
FF$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_FF$M35E, x = dat_H2S_FF$date, xout = FF$date))[,2] / 24 # g H2S /day

dat_ent_SF <- dat_ent[which(dat_ent$treatment == "slurryfunnels"),]
SF$enteric <- as.data.frame(approx(y = dat_ent_SF$CH4_enterisk * dat_ent_SF$No_pigs, x = dat_ent_SF$date, xout = SF$date))[,2]
dat_temp_SF <- dat_temp[which(dat_temp$treatment == "slurryfunnels"), c("date", "temp")]
SF <- merge(SF, dat_temp_SF, all = T)
dat_mass_SF <- dat_mass[which(dat_mass$treatment == "slurryfunnels"), c("date", "mass")]
SF <- merge(SF, dat_mass_SF, all = T)
dat_NH3_SF <- dat_NH3[which(dat_NH3$treatment == "slurryfunnels"),]
SF$NH3 <- as.data.frame(approx(y = dat_NH3_SF$NH3, x = dat_NH3_SF$date, xout = SF$date))[,2]
dat_weights_SF <- dat_weights[which(dat_weights$treatment == "slurryfunnels"),]
SF$weights <- as.data.frame(approx(y = dat_weights_SF$weights, x = dat_weights_SF$date, xout = SF$date))[,2]
dat_vfa_SF <- dat_vfa[which(dat_vfa$treatment == 'slurryfunnels'),]
SF <- merge(SF, dat_vfa_SF, all = T)
dat_H2S_SF <- dat_odor[which(dat_odor$treatment == 'Slurryfunnels'),]
SF$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_SF$M35E, x = dat_H2S_SF$date, xout = SF$date))[,2] / 24 # g H2S /day

dat_ent_ST <- dat_ent[which(dat_ent$treatment == "slurrytrays"),]
ST$enteric <- as.data.frame(approx(y = dat_ent_ST$CH4_enterisk * dat_ent_ST$No_pigs, x = dat_ent_ST$date, xout = ST$date))[,2]
dat_temp_ST <- dat_temp[which(dat_temp$treatment == "slurrytrays"), c("date", "temp")]
ST <- merge(ST, dat_temp_ST, all = T)
dat_mass_ST <- dat_mass[which(dat_mass$treatment == "slurrytrays"), c("date", "mass")]
ST <- merge(ST, dat_mass_ST, all = T)
dat_NH3_ST <- dat_NH3[which(dat_NH3$treatment == "slurrytrays"),]
ST$NH3 <- as.data.frame(approx(y = dat_NH3_ST$NH3, x = dat_NH3_ST$date, xout = ST$date))[,2]
dat_weights_ST <- dat_weights[which(dat_weights$treatment == "slurrytrays"),]
ST$weights <- as.data.frame(approx(y = dat_weights_ST$weights, x = dat_weights_ST$date, xout = ST$date))[,2]
dat_vfa_ST <- dat_vfa[which(dat_vfa$treatment == 'slurrytrays'),]
ST <- merge(ST, dat_vfa_ST, all = T)
dat_H2S_ST <- dat_odor[which(dat_odor$treatment == 'Slurrytrays'),]
ST$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_ST$M35E, x = dat_H2S_ST$date, xout = ST$date))[,2] / 24 # g H2S /day


C$CO2_enteric <- (0.136 * C$weights^0.573 ) * 1000
FF$CO2_enteric <- (0.136 * FF$weights^0.573 ) * 1000
ST$CO2_enteric <- (0.136 * ST$weights^0.573 ) * 1000
SF$CO2_enteric <- (0.136 * SF$weights^0.573 ) * 1000

# correct for enteric CH4
C$CH4_emis_rate <- C$CH4_rate - C$enteric
FF$CH4_emis_rate <- FF$CH4_rate - FF$enteric
SF$CH4_emis_rate <- SF$CH4_rate - SF$enteric
ST$CH4_emis_rate <- ST$CH4_rate - ST$enteric

#rename NH3
C$NH3_emis_rate <- C$NH3
FF$NH3_emis_rate <- FF$NH3
SF$NH3_emis_rate <- SF$NH3
ST$NH3_emis_rate <- ST$NH3

# correct for enteric CO2
C$CO2_emis_rate <- C$CO2E - C$CO2_enteric
FF$CO2_emis_rate <- FF$CO2E - FF$CO2_enteric
SF$CO2_emis_rate <- SF$CO2E - SF$CO2_enteric
ST$CO2_emis_rate <- ST$CO2E - ST$CO2_enteric


C$days <- as.double(difftime(C$date, C$date[1], units = "days"))
FF$days <- as.double(difftime(FF$date, FF$date[1], units = "days"))
SF$days <- as.double(difftime(SF$date, SF$date[1], units = "days"))
ST$days <- as.double(difftime(ST$date, ST$date[1], units = "days"))

################# setup optimization problem #################
meas <- C

# Create slurry_mass data frame
meas$time <- meas$days
meas$slurry_mass <- meas$mass

slurry_mass_dat <- meas[!is.na(meas$mass), c('date', 'time', 'slurry_mass')]
wash_dates <-ymd_hms(c("2020-08-13 07:00:00 UTC", "2020-11-05 02:00:00 UTC", "2021-01-28 02:00:00 UTC"))
slurry_mass_dat$wash_water <- 0
slurry_mass_dat$wash_water[slurry_mass_dat$date %in% wash_dates] <- 3500
slurry_mass_dat$date <- NULL

temp_dat <- meas[!is.na(meas$temp), c('time', 'temp')]

########################## Setup model parameters ##############################

wthr_pars = list(temp_C = 18.6, temp_air_C = 20, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs))
mng_pars = list(slurry_prod_rate = 0,  
                slurry_mass = slurry_mass_dat,          
                storage_depth = 0.6,        
                resid_depth = 0.035,         
                floor_area = 22.09,
                area = 22.09, 
                empty_int = 35,
                temp_C = temp_dat,
                temp_air_C = 18.6,
                wash_water = 0,
                wash_int = NA,
                rest_d = 7,
                RH = 90, 
                cover = NA,
                resid_enrich = 0.9,
                slopes = c(urea = NA, slurry_prod_rate = NA),
                scale = c(ks_coefficient = 1, qhat_opt = 0.4, xa_fresh = 4, yield = 1, alpha_opt = 1))
grz_pars = list(graze_start = "may",
                graze_days = 0,
                graze_hours = 0)
man_pars = list(conc_fresh = list(S2 = 0.01, urea = 2.4, SO4 = 0.2, TAN = 0.63, starch = 0, 
                                  VFA = 2.83, xa_dead = 0, CF = 0, CP = 0, NDF = 0, iNDF = 15, VSd = 73.8, VSd_A = 44.4), pH = 6.88, dens = 1000)

grp_pars = list(grps = c('m1','m2', 'sr1'),
                yield = c(default = 0.05, sr1 = 0.065),
                xa_fresh = c(default = 0.02),
                xa_init = c(all = 0.0001),
                decay_rate = c(all = 0.02),
                ks_coefficient = c(default = 1, sr1 = 0.4),
                qhat_opt = c(m1 = 3.6, m2 = 5.6 , m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8),
                T_opt = c(m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                T_min = c(m1 = 0, m2 = 8, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                T_max = c(m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                ki_NH3_min = c(all = 0.015),
                ki_NH3_max = c(all = 0.13),
                ki_NH4_min = c(all = 2.7),
                ki_NH4_max = c(all = 4.8),
                pH_upr = c(all = 8.0),
                pH_lwr = c(default = 6.0))
mic_pars = list(ks_SO4 = 0.0067,
                ki_H2S_meth = 0.23,
                ki_H2S_sr = 0.25,
                alpha_opt = c(xa_dead= 0.02, starch = 0.2, CF = 0.02, CP = 0.02, urea = 70, NDF = 0.02, iNDF = 0, VSd = 0.02),
                alpha_T_min = c(xa_dead= 0, starch = 0, CF = 0, CP = 0, urea = 0, NDF = 0, iNDF = 0, VSd = 0),
                alpha_T_opt = c(xa_dead= 50, starch = 50, CF = 50, CP = 50, urea = 50, NDF = 50, iNDF = 50, VSd = 50),
                alpha_T_max = c(xa_dead= 60, starch = 60, CF = 60, CP = 60, urea = 60, NDF = 60, iNDF = 60, VSd = 60))

chem_pars = list(COD_conv = c(CH4 = 0.2507, xa_dead = 0.73, NDF = 0.84, iNDF = 0.65, starch = 0.85, 
                              CF = 0.35, CP = 0.65, VFA = 0.93, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2, CO2_ureo = 1.57,
                              N_CP = 0.1014, C_xa_dead = 0.358, C_NDF = 0.376, C_iNDF = 0.358
                              , C_starch = 0.377, C_CF = 0.265, C_CP = 0.359 , C_VFA = 0.374, C_VSd = 0.344, C_N_urea = 0.429), 
                 kl = c(NH3 = 52, NH3_floor = 22, H2S = 0.02)) 

################### adjust for optimizer and exclude period 4 ###################
# Get emission data frame with average daily emission CH4 rate
days_sim <- max(meas$days[meas$period == 3], na.rm = T)

# NTS: dplyr alternative
meas$day <- ceiling(meas$days)
dat <- as.data.frame(summarise(group_by(meas, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE)))
# FRD: Need to change line below to include the last slurry_mass point at day 244.3...
# before it was days_sim (no ceiling)
dat <- dat[dat$day <= ceiling(days_sim), ] 
plot(CH4_emis_rate ~ day, data = dat, type = 'b', col = 'red')

## data.table package alternative is probably simpler
#dat <- data.table::as.data.table(C)[day <= days_sim, .(CH4_emis_rate1 = mean(CH4_emis_rate, na.rm = TRUE)), by = .(day = ceiling(day))]
#dat

# NTS: drop missing rates?

################### run optimization calculation ###################

## start with LGBS ##
new_pars <- data.frame(scale.qhat_opt = 0.4, scale.xa_fresh = 4, resid_enrich = 0.9, scale.alpha_opt = 1)
res <- 1 # make one so nrow of res is = nrow new_pars
loops <- 3 # number of times the parameters should be optimized in the optim loop
pars.cal <- log10(new_pars) # initial pars.cal

for (i in 1:loops){
  
  weights <- 1 * !is.na(dat$CH4_emis_rate)
  
  # run calibration function
  cal <- optim(par = pars.cal, 
               fn = resCalc,
               dat = dat,
               weights = weights,
               to = 'CH4_emis_rate', 
               mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, wthr_pars = wthr_pars, evap_pars = evap_pars,
               plot = F,
               method = 'L-BFGS-B',
               lower = log10(c(0.2, 1, 0.5, 1)), upper = log10(c(1, 50, 1, 10)),
               control = list(reltol = 0.01),
               hessian = TRUE
  )
  res1 <- cal$value
  res <- rbind(res, res1)
  new_pars1 <- 10^cal$par
  new_pars1 <- new_pars1[c("scale.qhat_opt", "scale.xa_fresh", "resid_enrich", "scale.alpha_opt")] # set to original order
  new_pars <- rbind(new_pars, new_pars1)
  pars.cal <- log10(new_pars[i + 1, sample(1:ncol(new_pars))]) # resample in preparation for next call. 
}

new_pars <- cbind(new_pars, res)
rmin <- which(new_pars$res == min(new_pars$res[-c(1)]))[1]


new_pars2 <- data.frame(scale.qhat_opt =  new_pars$scale.qhat_opt[rmin], scale.xa_fresh = new_pars$scale.xa_fresh[rmin],
                        resid_enrich = new_pars$resid_enrich[rmin], scale.alpha_opt = new_pars$scale.alpha_opt[rmin])
res <- 1 # make one so nrow of res is = nrow new_pars
loops <- 3 # number of times the parameters should be optimized in the optim loop
pars.cal <- log10(new_pars2) # initial pars.cal

for (i in 1:loops){
  
  # take previous new_pars values and randomize the order in which they are passed to optim.R
  
  weights <- 1 * !is.na(dat$CH4_emis_rate)
  
  # run calibration function
  cal <- optim(par = pars.cal, 
               fn = resCalc,
               dat = dat,
               weights = weights,
               to = 'CH4_emis_rate', 
               mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, wthr_pars = wthr_pars, evap_pars = evap_pars,
               plot = F,
               method = 'Nelder-Mead',
               control = list(reltol = 0.01),
               hessian = TRUE
  )
  res1 <- cal$value
  res <- rbind(res, res1)
  new_pars21 <- 10^cal$par
  new_pars21 <- new_pars21[c("scale.qhat_opt", "scale.xa_fresh", "resid_enrich", "scale.alpha_opt")] # set to original order
  new_pars2 <- rbind(new_pars2, new_pars21)
  pars.cal <- log10(new_pars2[i + 1, sample(1:ncol(new_pars2))]) # resample in preparation for next call. 
}

new_pars2 <- cbind(new_pars2, res)

rmin2 <- which(new_pars2$res == min(new_pars2$res[-c(1)]))[1]

new_pars3 <- data.frame(scale.qhat_opt = new_pars2$scale.qhat_opt[rmin2], scale.xa_fresh = new_pars2$scale.xa_fresh[rmin2], 
                        resid_enrich = new_pars2$resid_enrich[rmin2], scale.alpha_opt = new_pars2$scale.alpha_opt[rmin2])
res <- 1 # make one so nrow of res is = nrow new_pars
loops <- 3 # number of times the parameters should be optimized in the optim loop
pars.cal <- log10(new_pars3) # initial pars.cal

for (i in 1:loops){
  
  # take previous new_pars values and randomize the order in which they are passed to optim.R
  
  weights <- 1 * !is.na(dat$CH4_emis_rate)
  
  # run calibration function
  cal <- optim(par = pars.cal, 
               fn = resCalc,
               dat = dat,
               weights = weights,
               to = 'CH4_emis_rate', 
               mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, wthr_pars = wthr_pars, evap_pars = evap_pars,
               plot = F,
               method = 'L-BFGS-B',
               lower = log10(c(0.2, 1, 0.5)), upper = log10(c(1, 50, 1)),
               control = list(reltol = 0.01),
               hessian = TRUE
  )
  res1 <- cal$value
  res <- rbind(res, res1)
  new_pars31 <- 10^cal$par
  new_pars31 <- new_pars31[c("scale.qhat_opt", "scale.xa_fresh", "resid_enrich", "scale.alpha_opt")] # set to original order
  new_pars3 <- rbind(new_pars3, new_pars31)
  pars.cal <- log10(new_pars3[i + 1, sample(1:ncol(new_pars3))]) # resample in preparation for next call. 
}

new_pars3 <- cbind(new_pars3, res)

all_pars <- data.frame(new_pars, new_pars2, new_pars3)

}
