abm_barn <- function(dat, sim, days){
  
  # need empty dat frames for holding removed slurry and emission dat from the different sections. 
  rem_dat <- NULL
  barn_dat <- NULL
  
  # loop through sections
  for (w in 1:sim){
    
    # get all the section and storage specific data
    farm_dat <- get_farm(w, days, dat)
    wthr_pars <- farm_dat$wthr_pars
  
    barn <- abm(days, 1, wthr_pars = wthr_pars, add_pars = farm_dat$farm_dat)
 
    # add farm_dat needed for storage simulation and more
    barn <- cbind(barn, data.frame(farm_dat$extra_pars), source = 'barn')
    
    # slurry removed to a storage
    rem.rows <- which(barn$slurry_mass_eff > 0)
    slurry_mass_dat <- data.frame(time = c(barn$time[rem.rows]), slurry_mass = c(barn$slurry_mass_eff * barn$f_ex_storage)[rem.rows])
    slurry_digestate <- data.frame(time = c(barn$time[rem.rows]), slurry_mass = c(barn$slurry_mass_eff * barn$f_biogas)[rem.rows])
    
    xa_fresh <- barn[rem.rows,grepl("m[0-9]_eff_conc|sr[0-9]_eff_conc", names(barn))]
    names(xa_fresh) <- mgsub(names(xa_fresh), c("xa_", "_eff_conc"), c("", ""))
    xa_fresh <- as.list(xa_fresh)
    
    # NTS need to add time to conc_fresh, if variable fresh conc shall work in abm_storage
    conc_fresh <- list(sulfide = barn[rem.rows, "sulfide_eff_conc"], urea = barn[rem.rows, "urea_eff_conc"], 
                       sulfate = barn[rem.rows, "sulfate_eff_conc"], TAN = barn[rem.rows, "TAN_eff_conc"],
                       starch = barn[rem.rows, "starch_eff_conc"], VFA = barn[rem.rows, "VFA_eff_conc"], 
                       xa_dead = barn[rem.rows, "xa_dead_eff_conc"], CF = barn[rem.rows, "CF_eff_conc"], 
                       CP = barn[rem.rows, "CP_eff_conc"], RFd = barn[rem.rows, "RFd_eff_conc"], 
                       iNDF = barn[rem.rows, "iNDF_eff_conc"], VSd = barn[rem.rows, "VSd_eff_conc"], 
                       VSd_A = barn[rem.rows, 'VSd_A_eff_conc'], VSnd_A = barn[rem.rows, 'VSnd_A_eff_conc'],
                       ash = barn[rem.rows, 'ash_eff_conc'],
                       C = barn[rem.rows, "C_eff_conc"],
                       N = barn[rem.rows, "N_eff_conc"])
    
    ifelse(barn$use_comp[1] == 'no', conc_fresh[c('xa_dead', 'starch', 'CF', 'CP', 'RFd')] <- 0, conc_fresh[c('VSd')] <- 0) 
    
    rem_dat1 <- cbind(slurry_mass_dat, digestate = slurry_digestate, conc_fresh, xa_fresh, farm_dat$extra_pars)
    rem_dat <- rbind(rem_dat, rem_dat1)
    barn_dat <- as.data.frame(rbind(barn_dat, assign(paste("section_type", w, sep = ""), barn)))
    
  }
  
  return(list(barn_dat = barn_dat, rem_dat = rem_dat, conc_fresh = conc_fresh, xa_fresh = xa_fresh))
  
}
