get_farm <- function(w, days, dat){
  
  COD_conv <- c(CH4 = 0.2507, xa_dead = 0.73, RFd = 0.8444792, iNDF = 0.8444792, starch = 0.8444792, CF = 0.3117844, CP = 0.6541602,
                VFA = 0.9383125, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2, 
                n_CP = 0.1014)

  #get information about the different sections 
  sections <- read_excel(path = dat, sheet = "in-barn", skip = 1, col_names = TRUE)
  col_names <- c(t(sections[, (grepl(paste("Section ", w, sep=""), colnames(sections)))]))
  sections <- as.data.frame(t(sections[paste("Value ", w, sep="")]))
  colnames(sections) <- gsub("[ ]", "_", col_names)
  sections[!grepl("[a-zA-Z]", sections)] <- as.numeric(sections[!grepl("[a-zA-Z]", sections)]) 

  section_parms <- get_sections(sections)
  attach(section_parms)
  
  # modify parms
  vent_air <-       vent_air * n_anim/prod_area * 24
  wash_water <-     wash_water * n_anim
  
  # if natural ventilation occurs calculate temperatures from the outside temperature (variable)
  if(vent == "natural" & temp_air_C == "variable"){
    temp <- read_excel(dat, sheet = "outside temp", skip = 0, col_names = TRUE)
    temp_air_C_dat <- data.frame(time = temp$day, temp_air_C = temp$Air_temp_C + 2.4)
    t_add <- c(rep(seq(0, days - 365, 365), each = nrow(temp)))
    assign("temp_air_C_dat", do.call("rbind", replicate(days/365, eval(parse(text = "temp_air_C_dat")) , simplify = FALSE)))
    temp_air_C_dat$time <- temp_air_C_dat$time + t_add
    temp_air_C <- mean(temp_air_C_dat$temp_air_C) + 2.4
    temp_C_dat <- temp_air_C_dat
  }
  
  ### check some arguments ####
  if(biogas == "none" && f_biogas > 0) stop("slurry is exported to biogas, but biogas is selected as \"none\"")
  if(sum_feed < 99 | sum_feed > 101) stop("feed componenets does not add up to 100%")
  if(biogas == "yes" && between(removal_biogas, 40, 80) == FALSE) warning("the removal efficiency in your biogas plant is unusually high or low")
  if(isTRUE(sections < 0)) stop("some section values are negative")

  storages <-               read_excel(path = dat, sheet = "out-of-barn", skip = 1, col_names = TRUE)
  storages <-               storages[, !grepl("barn", colnames(storages))]
  col_names <-              c(t(storages[, (grepl(paste("Storage ", storage_tank_ID, sep=""), colnames(storages)))]))
  storages <-               as.data.frame(t(storages[,grepl(paste("Value ", storage_tank_ID, sep=""), colnames(storages))]))
  colnames(storages) <-     gsub("[ ]", "_", col_names)
  storages[!grepl("[a-zA-Z]", storages)] <- as.numeric(storages[!grepl("[a-zA-Z]", storages)]) 
  
  storage_parms <- get_storages(storages)
  attach(storage_parms)
  
  # get information about composition of excreted manure
  if (class_anim == "pig" & type_anim == "finisher" | type_anim == "piglet") type <- "Grow"
  if (class_anim == "pig" & type_anim != "finisher" & type_anim != "piglet") type <- "Adult"
  if (class_anim == "cattle") type <- ""
  

  feed_dat <- read_excel(dat, sheet = paste(class_anim, "_feed", sep = ""), 
                         skip = 0, col_names = TRUE) 
  names <- c("Feed_name", "DM", "CP", "Ex_CF", paste("Ex_CP", type, sep = ""), "Ex_Starch", paste("Ex_deg_ResFib", type, sep = ""), paste("Ex_nondeg_ResFib", type, sep = ""), 
    paste("Ex_OM", type, sep = ''), 'OM', 'N', 'CF', 'NDF', 'Sugar', 'Starch', 'iNDF', 'UndigStarch', 'ResFib', 'Potassium', 'Sodium')
  
  feed_dat <- feed_dat[which(feed_dat$Feed_name %in% colnames(sections)), which(colnames(feed_dat) %in% names)] %>% distinct(Feed_name, .keep_all = TRUE)
  # handle non-numerics 
  feed_dat_num <- sections[, which(colnames(sections) %in% feed_dat$Feed_name)]
  feed_dat_num[sapply(feed_dat_num, is.character)] <- as.numeric(feed_dat_num[sapply(feed_dat_num, is.character)])
  feed_comp <- data.frame(t(feed_dat_num/100))
  feed_comp <- data.frame(cbind(feed_comp, Feed_name = rownames(feed_comp)))
  feed_dat <- merge(feed_dat, feed_comp[], by = 'Feed_name')
  colnames(feed_dat)[ncol(feed_dat)] <- 'proportion'
  feed_dat <- filter(feed_dat, proportion > 0)
  feed_dat[, -1] <- lapply(feed_dat[, -1], as.numeric)

  excreta_dat <- feedFun(class_anim, type_anim, type, feed_dat, feed_intake, milk_prod, body_weight, batch_time)
  
  attach(excreta_dat)
  
  feed.names <- c('CP', 'CF', 'iNDF', 'ResFib', 'Starch', 'Sugar')
  feed_spill <- colSums(feed_dat[, feed.names] * feed_dat$proportion, na.rm = T)/1000 * feed_DM_conc * feed_intake * (feed_spill_frac/100) * 1000 # g of feed components spilled during batch
  excreta_dat$feed_spill <- as.data.frame(feed_spill)
  CH4_ent <- CH4_ent * 365 / (batch_time + rest_d)
  
  # using feed composition and functions for urine and feces volumes and masses 
  if(use_comp == 'yes'){
    slurry <- (feces + urine + bedding + sum(feed_spill)/feed_DM_conc/1000)
    slurry_tot <- slurry * n_anim
    slurry_tot_day <- slurry_tot/batch_time
  }
  
  # using normative system values
  if(use_comp == 'no'){
    if(class_anim == 'pig'){
      urine <- feed_intake * feed_DM * spec_urine
      feces <- feed_intake * feed_DM * (1 - feed_digest)/ (TS_feces/100)
      VS_anim <- feed_intake * feed_DM * (1 - feed_digest)
      slurry <- urine + feces + bedding
      VSd <- VS_anim / slurry / COD_conv[['VS']] * 1000 * VSd_VS
      slurry_tot <- slurry * n_anim
      slurry_tot_day <- slurry_tot / batch_time
      ash <- VS_anim / slurry * 0.176 *1000 # g/kg slurry
      urea <- 3.17
      TAN <- 0
    } else if(class_anim == 'cattle'){
      VS_anim <- feed_intake * (1 - feed_digest)
      feces <- feed_intake * (1 - feed_digest)/(TS_feces/100)
      urine <- feces * spec_urine
      slurry <- feces + urine + bedding
      VSd <- VS_anim / slurry / COD_conv[['VS']] * 1000 * VSd_VS
      slurry_tot <- slurry * n_anim
      slurry_tot_day <- slurry_tot / batch_time
      ash <- VS_anim / slurry * 0.176 *1000 # g/kg slurry
      urea <- 2.5
      TAN <- 0
    }
  }

  #if (class_anim == "pig" && type_anim == "finisher"){
  #  slopes.slurry_prod_rate <- 0.0668 * n_anim
  #} else {
  #  slopes.slurry_prod_rate <- NA
  #}
  slopes.slurry_prod_rate <- NA
  
  # convert to COD units in slurry (g COD or g element/kg slurry)
  # bedding material adds a little to CP
  
  CP <- (feces_comp[['CP_feces']] + feed_spill[['CP']] + (bedding * bedding_TS * 0.005 * 6.25)/1000) / COD_conv[['CP']]/slurry
  CF <- ((feces_comp[['CF_feces']] + feed_spill[['CF']]) / COD_conv[['CF']])/slurry
  starch <- ((feces_comp[['starch_feces']] + feed_spill[['Starch']] + feed_spill[['Sugar']]) / COD_conv[['starch']])/slurry
  RFd <- ((feces_comp[['RFd_feces']] + feed_spill[['ResFib']] - feed_spill[['iNDF']])/ COD_conv[['RFd']])/slurry
  iNDF <- (feces_comp[['iNDF_feces']] + feed_spill[['iNDF']] + (bedding * bedding_TS - bedding * bedding_TS * 0.005 * 6.25)/1000) / COD_conv[['iNDF']]/slurry
  VFA <- (feces_comp[['VFA_feces']] / COD_conv[['VFA']])/slurry
  
  if(use_comp == 'yes'){
    if(class_anim == 'pig'){
      urea <- urea_N/slurry
      TAN <- 0
      ash <- sum(CP * COD_conv[['CP']], CF * COD_conv[['CF']], 
                 starch * COD_conv[['starch']], RFd * COD_conv[['RFd']],
                 iNDF * COD_conv[['iNDF']]) * 0.176
    } else if(class_anim == 'cattle'){
      urea <- urea_N/slurry
      CP <- CP + (urine_N - urea_N)/slurry * 4 * COD_conv[['VFA']] # 4 is g molecule per g N. the molecule is assumed to be average of hippuric acid, allantoin, creatinine = C5.67H7.33N2.67O2.33, with COD of 1.08
      TAN <- 0
      ash <- sum(CP * COD_conv[['CP']], CF * COD_conv[['CF']], 
                 starch * COD_conv[['starch']], RFd * COD_conv[['RFd']],
                 iNDF * COD_conv[['iNDF']]) * 0.176
    }
  }
  
  sulfate <- 0.2
  
  # when individual components are not used we need VSd
  
  
  if (barn_acid == "none"){
    if (class_anim == "pig") pH <- 7
    if (class_anim == "cattle") pH <- 7
    sulfate_conc <- sulfate
  } else if (barn_acid == "continuous H2SO4 acidification"){
    pH <- H2SO4_titrat(conc_SO4 = barn_acid_dose * 32.06/98.08, class_anim = class_anim)$pH  
    sulfate_conc <- barn_acid_dose * 32.06/98.08
  }

  if (storage_acid == "none"){
    if (class_anim == "pig") pH_storage <- 7.2
    if (class_anim == "cattle") pH_storage <- 7
    pH_digestate <- 7.66
    sulfate_conc_storage <- sulfate
    sulfate_conc_digestate <- 0
  } 
  if (storage_acid == "continuous H2SO4 acidification"){
    pH_storage <- H2SO4_titrat(conc_SO4 = storage_acid_dose * 32.06/98.08, class_anim = class_anim)$pH
    pH_digestate <- H2SO4_titrat(conc_SO4 = digestate_acid_dose * 32.06/98.08, class_anim = "degassed")$pH
    sulfate_conc_storage <- storage_acid_dose * 32.06/98.08 
    sulfate_conc_digestate <- digestate_acid_dose * 32.06/98.08 
  }
  if (storage_acid == "single dose acid"){

  }
  
  if (biogas == 'none'){
    f_storage <- 1
    f_biogas <- 0
  }
  
  # height of slurry
  height <- resid_depth + (slurry_tot_day/1000/area) * empty_int
  
  # warn about overflowing the storage
  if (height < resid_depth | height > storage_depth){
    stop('In get_farm(), slurry height is lower than minimum slurry height or larger than maximum slurry heigh, problem with slurry production versus dimensions og slurry pit - check data input')
  }
  
  average_height <- resid_depth + (height - resid_depth)/2
  # The following is a rough estimate of the average slurry temp reduction. Data from Holm et al. 2017 is used. 
  # It is assumed that cooling effect is proportional to reduction in temperature (which is most likely not correct)
  # there is no temperature gradient from the 0.1 m slurry height to the surface.
  
  if (barn_cool != "none" & average_height <= 0.1){
    cool_temp <- (average_height * 21 - 4.5) * cool_eff/26 
  } else if(barn_cool != "none" & average_height > 0.1){
    cool_temp <- -(2.4 * cool_eff/26 )
  } else {
    cool_temp <- 0
  }
  
  # sort out temperature if it is variable. 
  if(exists("temp_air_C_dat")){
  temp_C_dat$temp_C <- temp_air_C_dat$temp_air_C + cool_temp * cool_days/365
  temp_C_dat <- as.data.frame(temp_C_dat)
  temp_air_C_wthr <- mean(temp_air_C_dat$temp_air_C)
  temp_C_wthr <- mean(temp_C_dat$temp_C)
  } else {
  temp_C_dat <- temp_air_C + cool_temp * cool_days/365  
  temp_air_C_dat <- temp_air_C
  temp_air_C_wthr <- temp_air_C_dat
  temp_C_wthr <- temp_C_dat
  }
  
  slurry_prod_rate <- slurry_tot_day
  floor_area <- prod_area
  section_ID <- w

  slopes <- c(urea = NA, slurry_prod_rate = slopes.slurry_prod_rate)
  
  # Arrhenius VS
  
  if(use_comp == 'no'){
    VSd_A <- VSd * COD_conv[['VS']]
    VSnd_A <- VSd_A/VSd_VS - VSd_A
    } else{
    VSd <- CP + CF + starch + RFd 
    VSd_A <- VSd * COD_conv[['VS']]
    VSnd_A <- VSd_A/VSd_VS - VSd_A
  }


  if(use_comp == 'yes'){
  conc_fresh <- list(sulfide = 0, urea = urea, sulfate = sulfate_conc, TAN = TAN, CF = CF, 
                     CP = CP, starch = starch, RFd = RFd, iNDF = iNDF, VFA = VFA, 
                     xa_dead = 0, VSd = 1e-10, VSd_A = VSd_A, VSnd_A = VSnd_A, ash = ash)
  } else if(use_comp == 'no'){
  conc_fresh <- list(sulfide = 0, urea = urea, sulfate = sulfate_conc, TAN = TAN, CF = 0, 
                     CP = 0, starch = 0, RFd = 0, iNDF = iNDF, VFA = VFA, 
                     xa_dead = 0, VSd = VSd, VSd_A = VSd_A, VSnd_A = VSnd_A, ash = ash)  
  }
    
  lnA <- c('VSd_A' = lnA)
  E_CH4 <- c('VSd_A' = E_CH4)

  # correct slurry_prod_rate if grazing occurs (for cattle). slurry excreted from the cattle would be the same but the
  # emission from the barn would be affected and moved to outside. But we cannot have another would be the same, but emission from the barn should be corrected
  # since 

  # compile farm specific pars for abm simulation, correct radiation and rain because it is inside a building.
  wthr_pars <- list(temp_C = temp_C_wthr, temp_air_C = temp_air_C_wthr, RH = 90, rain = 0, pres_kpa = 101, rs = 20 * 0.5)
  
  farm_dat <- list(slopes = slopes, conc_fresh = conc_fresh, slurry_prod_rate = slurry_prod_rate,
                   storage_depth = storage_depth, resid_depth = resid_depth, rain = 0,
                   empty_int = empty_int, area = area, floor_area = floor_area, 
                   wash_water = wash_water, wash_int = wash_int, rest_d = rest_d, 
                   temp_air_C = temp_air_C_dat, pH = pH, temp_C = temp_C_dat,
                   graze_days = graze_days, graze_hours = graze_hours, graze_start = graze_start,
                   lnA = lnA, VS_CH4 = VS_CH4, E_CH4 = E_CH4, kl.NH3 = kl_NH3_pit, kl.NH3_floor = kl_NH3_floor)
 
  extra_pars <- list(type_anim = type_anim, n_anim = n_anim, batch_time = batch_time, 
                       prod_area = prod_area, barn_acid = barn_acid, barn_acid_dose = barn_acid_dose,
                       wash_water = wash_water, f_ex_storage = f_ex_storage, f_biogas = f_biogas, removal_biogas = removal_biogas, 
                       section_ID = section_ID, rest_d = rest_d, ex_storage_rain = ex_storage_rain, ex_storage_radiation = ex_storage_radiation,
                     ex_storage_area = ex_storage_area, ex_storage_depth = ex_storage_depth,  
                       storage_ID = storage_ID, storage_acid = storage_acid, storage_acid_dose = storage_acid_dose,
                       cover_s = cover_s, venting_s = venting_s, venting_eff_s = venting_eff_s, 
                       flaring_s = flaring_s, flaring_eff_s = flaring_eff_s,  
                       pH_storage = pH_storage, sulfate_conc_storage = sulfate_conc_storage,
                     digestate_area = digestate_area, digestate_depth = digestate_depth, 
                       digestate_ID = digestate_ID, digestate_acid = digestate_acid, digestate_acid_dose = digestate_acid_dose,
                       cover_d = cover_d, venting_d = venting_d, venting_eff_d = venting_eff_d, 
                       flaring_d = flaring_d, flaring_eff_d = flaring_eff_d,
                     pH_digestate = pH_digestate, sulfate_conc_digestate = sulfate_conc_digestate,
                     app1s = app1s, app2s = app2s, app3s = app3s, app4s = app4s, app5s = app5s, app6s = app6s, app7s = app7s,
                     app_t1s = app_t1s, app_t2s = app_t2s, app_t3s = app_t3s, app_t4s = app_t4s, app_t5s = app_t5s, app_t6s = app_t6s, app_t7s = app_t7s,
                     app1d = app1d, app2d = app2d, app3d = app3d, app4d = app4d, app5d = app5d, app6d = app6d, app7d = app7d,
                     app_t1d = app_t1d, app_t2d = app_t2d, app_t3d = app_t3d, app_t4d = app_t4d, app_t5d = app_t5d, app_t6d = app_t6d, app_t7d = app_t7d,
                     lnA. = lnA., E_CH4 = E_CH4, VS_CH4 = VS_CH4, kl_NH3_storage = kl_NH3_storage, use_comp = use_comp, 
                     CH4_ent = CH4_ent
                     )
  prod_dat <- data.frame(batch_time = batch_time, bedding = bedding, n_anim = n_anim, type_anim = type_anim, 
                            empty_int = empty_int, rest_d = rest_d,
                   area = area, floor_area = floor_area, storage_depth = storage_depth)
  excreta_dat$prod_dat <- prod_dat
  
  detach(section_parms)
  detach(storage_parms)
  
  
  return(c(list(farm_dat = farm_dat, extra_pars = extra_pars, wthr_pars = wthr_pars, excreta_dat = excreta_dat)))

}

