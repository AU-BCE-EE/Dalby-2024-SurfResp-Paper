abm_storage <- function(days, years, rem_dat, conc_fresh, xa_fresh, temps, doy = doy){
  
  #data frames for holding results. 
  stor_dat <- NULL
  storage_dat <- NULL
  digestate_dat <- NULL
  
  if(exists('slurry_mass')) {rm('slurry_mass')}
  
  # repeat this loop for number of storages simulated
  for (i in 1:length(unique(rem_dat$storage_ID))){
    
    if(any(rem_dat$slurry_mass) > 0){
      
      # take out slurry mass to the first storage and sort with respect to time
      stor_dat <- rem_dat %>% filter(storage_ID == i) %>% arrange(time) 
      specs <- stor_dat[1,]
      # calculated slurry weighted average from different sections of conc_fresh and xa_fresh for input to storage 
      inflows <- as.data.frame(stor_dat) %>% select(names(conc_fresh), names(xa_fresh), slurry_mass, section_ID) %>% 
        group_by(section_ID) %>%
        summarise(across(c(names(conc_fresh), names(xa_fresh)), mean), across(slurry_mass, sum)) %>%
        mutate(across(-c(slurry_mass, section_ID), function(x) x*slurry_mass/sum(slurry_mass))) %>% 
        summarise_all(sum) 

      storage_conc_fresh <- as.list(inflows[names(conc_fresh)])
      storage_xa_fresh <- setNames(as.numeric(inflows[names(xa_fresh)]), names(xa_fresh))

      # sort out slurry_mass data.frame for storage
      slurry_mass_dat <- rbind(c(0,0), stor_dat[, c("time", "slurry_mass")] %>% mutate_at(vars(slurry_mass), cumsum))
      
      # interpolate slurry mass with "constant" method to ensure slurry mass is added in portions, (not continously)
      xout <- sort(c(slurry_mass_dat$time, 1:slurry_mass_dat$time[length(slurry_mass_dat$time)], decreasing = F))
      slurry_mass_dat <- data.frame(time = xout, 
                                          slurry_mass = approx(slurry_mass_dat$time, slurry_mass_dat$slurry_mass, 
                                                               xout = xout, method = "constant")$y)  %>% distinct()
      
      ifelse(specs$cover_s == "tent", rain <- 0, rain <- specs$ex_storage_rain)
      ifelse(specs$cover_s == "tent", rs <- 0, rs <- specs$ex_storage_radiation)
      rain_dat <- rain * slurry_mass_dat$time # assumed constant
      evap_dat <- rep(et(cmakk = 0.7, temp_C = mean(temps$temp_C_dat$temp_C[2:12]), pres_kpa = 101, rs = rs)) * slurry_mass_dat$time
      
      # Subtracting evap and adding rain from the slurry mass (cumulatives) 
      slurry_mass_dat$slurry_mass <- slurry_mass_dat$slurry_mass - evap_dat * specs$ex_storage_area + rain_dat * specs$ex_storage_area
      
      # handling emptying events / application times. The year starts in March here
      slurry_mass_dat <- slurry_app(days = days, begin = 'March', specs = specs, from = 'storage', slurry = slurry_mass_dat)
      
      if(any(slurry_mass_dat$slurry_mass > specs$ex_storage_area * specs$ex_storage_depth * 1000)){
        stop("slurry mass exceeded storage capacity") 
      }
      
      # if storage acidification is applied then correct SO4_conc_fresh passed to storage. 
      if(specs$storage_acid == 'continuous H2SO4 acidification'){
        storage_conc_fresh$sulfate <- specs$sulfate_conc_storage
      }  

      wthr_pars <- list(temp_C = mean(temps$temp_C_dat$temp_C[2:12]), temp_air_C = mean(temps$temp_air_C_dat$temp_air_C[2:12]), RH = 90, rain = rain, pres_kpa = 101, rs = rs)  
      
      # add arrhenius specs
      storage_conc_fresh$VSd_A <- storage_conc_fresh$VSd_A + storage_conc_fresh$VSnd_A

      # run model for storage
      storage <- cbind(abm(days = days, 1, wthr_pars = wthr_pars, add_pars = list(slurry_mass = slurry_mass_dat, 
                                    conc_fresh = storage_conc_fresh, xa_fresh = storage_xa_fresh, 
                                    pH = specs$pH_storage, cover = specs$cover_s, storage_depth = specs$storage_depth, 
                                    area = specs$ex_storage_area, floor_area = 0, rain = rain, 
                                    temp_C = temps$temp_C_dat, lnA.VSd_A = specs$lnA., E_CH4.VSd_A = specs$E_CH4, VS_CH4 = specs$VS_CH4,
                                    kl.NH3 = specs$kl_NH3_storage)), source = 'storage', storage_ID = specs$storage_ID, 
                                    prod_area = specs$prod_area, n_anim = specs$n_anim, batch_time =  specs$batch_time, wash_water = 0)
      
      storage$temp_air_C <- approx(temps$temp_air_C_dat$time, temps$temp_air_C_dat$temp_air_C, storage$time)$y
      
      if(specs$venting_s == TRUE | specs$flaring_s == TRUE){
        storage <- mitigate_storage(specs = specs, from = 'storage', storage = storage)
      }
      
      storage_dat <- as.data.frame(rbind(storage_dat, storage))
    }
  }
  
  # if biogas is used then repeat everything as was done for storage. Some exemption see below.  
  for (i in 1:length(unique(rem_dat$digestate_ID[rem_dat$f_biogas > 0]))){  
    
    if (any(rem_dat$digestate.slurry_mass) > 0){
      
      stor_dat <- rem_dat %>% filter(digestate_ID == i) %>% arrange(time) 
      specs <- stor_dat[1,]
      
      # calculated slurry weighted average from different sections of conc_fresh and xa_fresh for input to storage 
      inflows <- as.data.frame(stor_dat) %>% select(names(conc_fresh), names(xa_fresh), digestate.slurry_mass, section_ID) %>% 
        group_by(section_ID) %>%
        summarise(across(c(names(conc_fresh), names(xa_fresh)), mean), across(digestate.slurry_mass, sum)) %>%
        mutate(across(-c(digestate.slurry_mass, section_ID), function(x) x*digestate.slurry_mass/sum(digestate.slurry_mass))) %>% 
        summarise_all(sum) 
      
      digestate_conc_fresh <- as.list(inflows[names(conc_fresh)])
      
      # Sulfate is completely reduced to H2S in anaerobic digestion
      digestate_conc_fresh$sulfate <- 0
      
      digestate_xa_fresh <- setNames(as.numeric(inflows[names(xa_fresh)]), names(xa_fresh))
      
      # sort out slurry_mass for digestate storage
      slurry_digestate <- rbind(c(0,0), stor_dat[, c("digestate.time", "digestate.slurry_mass")] %>% mutate_at(vars(digestate.slurry_mass), cumsum))
      xout <- sort(c(slurry_digestate$digestate.time, 1:slurry_digestate$digestate.time[length(slurry_digestate$digestate.time)], decreasing = F))
      slurry_digestate_dat <- data.frame(time = xout, slurry_mass = approx(slurry_digestate$digestate.time, slurry_digestate$digestate.slurry_mass, 
                                                                                 xout = xout, method = "constant")$y) %>% distinct()
      # calculate rain and evap data
      ifelse(specs$cover_d == "tent", rain <- 0, rain <- specs$ex_storage_rain)
      ifelse(specs$cover_d == "tent", rs <- 0, rs <- specs$ex_storage_radiation)
      rain_dat <- rain * slurry_digestate_dat$time # assumed constant
      evap_dat <- rep(et(cmakk = 0.7, temp_C = mean(temps$d.temp_C_dat$temp_[2:12]), pres_kpa = 101, rs = rs)) * slurry_digestate_dat$time
      
      # Subtracting evap and adding rain from the slurry mass (cumulatives) 
      slurry_digestate_dat$slurry_mass <- slurry_digestate_dat$slurry_mass - evap_dat * specs$digestate_area + rain_dat * specs$digestate_area
      
      # slurry removal due to soil application 
      slurry_digestate_dat <- slurry_app(days = days, begin = 'March', specs = specs, from = 'digestate', slurry = slurry_digestate_dat)
     
      if(any(slurry_digestate_dat$slurry_mass > specs$digestate_area * specs$digestate_depth * 1000)){
        stop("slurry mass exceeded storage capacity") 
      }
      
      # !OBS! Here we subtract slurry organic matter which was consumed in the digester
      # resulting in a lower conc_fresh state variables passed to the digestate storage
      # iNDF is not consumed in digester and therefore the fractoin of degradable OM is calculted first
      C_before <- digestate_conc_fresh$C
      N_before <- digestate_conc_fresh$N
      
      f_COD <- (digestate_conc_fresh[['iNDF']] + digestate_conc_fresh[['VFA']] + digestate_conc_fresh[['starch']] + digestate_conc_fresh[['RFd']] + 
                  digestate_conc_fresh[['CP']] + digestate_conc_fresh[['CF']] + digestate_conc_fresh[['VSd']])/
        (digestate_conc_fresh[['VFA']] + digestate_conc_fresh[['starch']] + digestate_conc_fresh[['RFd']] + digestate_conc_fresh[['CP']] + 
           digestate_conc_fresh[['CF']] + digestate_conc_fresh[['VSd']])
      
      dCOD <- c('VFA', 'starch', 'RFd', 'CP', 'CF', 'VSd')
      
      digestate_conc_fresh <- as.data.frame(digestate_conc_fresh)
      
      # simulate digestion by reducing OM in state variables according to "stor_dat$removal_biogas"
      digestate_conc_fresh[dCOD] <- digestate_conc_fresh[dCOD] * (1 - f_COD * specs$removal_biogas/100)
      digestate_conc_fresh <- as.list(digestate_conc_fresh)
      
      if(specs$digestate_acid == 'continuous H2SO4 acidification'){
        digestate_conc_fresh$sulfate <- specs$sulfate_conc_digestate
      }  
      
      wthr_pars <- list(temp_C = mean(temps$d.temp_C_dat$temp_C[2:12]), temp_air_C = mean(temps$temp_air_C_dat$temp_air_C[2:12]), RH = 90, rain = rain, pres_kpa = 101, rs = 20)  
      
      digestate_conc_fresh$VSd_A <- digestate_conc_fresh$VSd_A + digestate_conc_fresh$VSnd_A
      
      digestate <- cbind(abm(days = days, 1, wthr_pars = wthr_pars, add_pars = list(slurry_mass = slurry_digestate_dat, conc_fresh = digestate_conc_fresh, xa_fresh = digestate_xa_fresh, 
                                                                              pH = specs$pH_digestate, cover = specs$cover_d, storage_depth = specs$digestate_depth, rain = rain,
                                                                              area = specs$digestate_area, floor_area = 0, temp_C = temps$d.temp_C_dat, 
                                                                              lnA.VSd_A = specs$lnA., E_CH4.VSd_A = specs$E_CH4, VS_CH4 = 16.67)), source = 'digestate', 
                                                                              digestate_ID = specs$digestate_ID, prod_area = specs$prod_area, n_anim = specs$n_anim, batch_time = 0, wash_water = 0,
                                                                              C_before = C_before, N_before = N_before)
      
      if(specs$venting_s == TRUE | specs$flaring_s == TRUE){
        digestate <- mitigate_storage(specs = specs, from = 'digestate', storage = digestate)
      }
      
      digestate_dat <- as.data.frame(rbind(digestate_dat, digestate))
    }
  }

  return(list(storage_dat = storage_dat, digestate_dat = digestate_dat))
}

