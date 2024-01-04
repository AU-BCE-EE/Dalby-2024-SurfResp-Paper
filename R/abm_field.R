field <- function(dat, years){
  
  library(ALFAM2)

  dat$year <- ceiling(dat$time/365)
  dat$dm_eff_conc <- (dat$ash_conc + dat$VS_conc) * !is.na(dat$ash_eff_conc)
  
  field_input <- dat %>% group_by(storage_ID, year) %>% 
    summarise(air.temp = sum(slurry_mass_eff * temp_air_C) / sum(slurry_mass_eff),
              man.tan = sum(slurry_mass_eff * TAN_eff_conc, na.rm = T) / sum(slurry_mass_eff), 
              man.dm = sum(slurry_mass_eff * dm_eff_conc, na.rm = T) / sum(slurry_mass_eff)/10,
              slurry_out = sum(slurry_mass_eff), 
              man.ph = pH) %>% distinct() %>% 
    mutate(ctime = 168, app.mthd = 'trailing hose', TAN.app = 100, wind.2m = 3, tan_tot = man.tan/1000 * slurry_out) %>%
    ungroup() %>% filter(year == 3)
  
  NH3_kg_yr_field <- alfam2(as.data.frame(field_input), pars = alfam2pars02, 
                  app.name = 'TAN.app', time.name = 'ctime')$er * field_input$tan_tot

  return(NH3_kg_yr_field)
}


