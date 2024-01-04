# Author: Frederik Dalby
# Last modified 09-11-2022
abm_farm <- function(dat, storage_mode = TRUE, update_feed = FALSE, years = 3){

  days <- 365 * years
  
  # merge_rawfeed takes the raw data for pig digestibilities and writes it into the demo sheet. 
  if (isTRUE(update_feed)){
    merge_rawfeed('../farm_dat/pigs/')
  }

  # read demo and determine number of abm runs that has to be done for the barn (due to different sections) 
  s <- read_excel(dat, sheet = "in-barn", skip = 1, col_names = TRUE)
  sim <- as.numeric(sum(grepl("Section ", colnames(s))))
 
  # simulate all sections in the barn
  abm_barn_out <- abm_barn(dat = dat, sim = sim, days = days)

  # pull out results, both barn emission (used in calcNorm.R) and removed slurry (used in abm_storage.R)
  rem_dat <- abm_barn_out$rem_dat
  barn_dat <- abm_barn_out$barn_dat
  conc_fresh <- abm_barn_out$conc_fresh
  xa_fresh <- abm_barn_out$xa_fresh

  if (isTRUE(storage_mode)){
  
    # take out temperatures and repeat for years.
    temps <- get_temp(days = days, dat = dat, sheet = "outside temp")
    
    # simulate storages
    abm_storage_out <- abm_storage(days = days, years = years, rem_dat = rem_dat, conc_fresh = conc_fresh, xa_fresh = xa_fresh, temps = temps, doy = doy) 
  
  }
  
  # pull out storage data for calculating emission factors in calcNorm
  storage_dat <- abm_storage_out$storage_dat
  digestate_dat <- abm_storage_out$digestate_dat
  
  #if (isTRUE(field_mode))
  
  norm <- calcNorm(barn_dat, storage_dat, digestate_dat, ave, days, years)

  
  #norm$emission_dat$NH3_kg_yr[4] <- abm_field(storage_dat, years)
  return(norm)
}


