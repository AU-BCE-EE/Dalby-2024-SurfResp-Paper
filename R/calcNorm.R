calcNorm <- function(barn_dat, storage_dat, digestate_dat, ave, days, years){
  
  norm_vars <- c('time', 'CH4_emis_rate', 'CO2_emis_rate', 'NH3_emis_rate', 'N2O_emis_rate', 'CH4_A_emis_rate',
                 'N_load_cum', 'C_load_cum',
                 'n_anim', 'area', 'prod_area', 'slurry_prod_rate', 'evap', 'rain',
                 'slurry_mass', 'slurry_mass_eff', 'N_conc_fresh', 'C_conc_fresh', 
                 'source', 'section_ID', 'storage_ID', 'digestate_ID', 'C_before', 'N_before',
                 'conc_fresh_TAN', 'conc_fresh_urea', 'C_eff_conc','N_eff_conc',
                 'TAN_eff_conc','urea_eff_conc','wash_water', 'batch_time', 'CH4_ent')
  
  # make separate frame with enteric data
  ent_dat <- barn_dat[, names(barn_dat) %in% norm_vars] %>% 
    mutate(source = 'enteric', CH4_emis_rate = CH4_ent/365 * n_anim, 
           CH4_A_emis_rate = CH4_ent/365 * n_anim, NH3_emis_rate = 0, CO2_emis_rate = 0)
  
  barn_dat$storage_ID <- NA
  barn_dat$digestate_ID <- NA 

  norm_dat <- bind_rows(ent_dat[, names(ent_dat) %in% norm_vars], barn_dat[, names(barn_dat) %in% norm_vars], storage_dat[, names(storage_dat) %in% norm_vars], digestate_dat[, names(digestate_dat) %in% norm_vars]) %>% 
    mutate(year = ceiling(time/365))
  
  .rain <- expression(rain[1] * area[1]/1000 * 365) # rain in kg per year
  .evap <- expression(evap[1] * area[1]/1000 * 365) # evaporation in kg per year
  .wash_water <- expression(sum(wash_water/batch_time/1000))
  .slurry_produced <- expression(sum(slurry_mass_eff, na.rm = T)/1000)
  
  C_load <- expression((max(C_load_cum) - min(C_load_cum))/1000) # kg/yr
  N_load <- expression((max(N_load_cum) - min(N_load_cum))/1000) # kg/yr

  emission_dat <- norm_dat %>% 
                  group_by(source, year, section_ID, storage_ID, digestate_ID) %>%
                  summarise(CH4_kg_yr = mean(CH4_emis_rate, na.rm = T)/1000 * 365, CH4_kg_anim_yr = mean(CH4_emis_rate, na.rm = T)/1000 * 365 /n_anim, CH4_kg_parea_yr = mean(CH4_emis_rate, na.rm = T)/1000 * 365 /prod_area, 
                            CH4_kg_m3_ex_yr = mean(CH4_emis_rate, na.rm = T)/1000 * 365 /eval(.slurry_produced), CH4_kg_m3_yr = mean(CH4_emis_rate, na.rm = T)/mean(slurry_mass, na.rm = T) * 365,
                            
                            CH4_A_kg_yr = mean(CH4_A_emis_rate, na.rm = T)/1000 * 365, CH4_A_kg_anim_yr = mean(CH4_A_emis_rate, na.rm = T)/1000 * 365 /n_anim, CH4_A_kg_parea_yr = mean(CH4_A_emis_rate[source != 'enteric'], na.rm = T)/1000 * 365 /prod_area, 
                            CH4_A_kg_m3_ex_yr = mean(CH4_A_emis_rate, na.rm = T)/1000 * 365 /eval(.slurry_produced), CH4_A_kg_m3_yr = mean(CH4_A_emis_rate, na.rm = T)/mean(slurry_mass, na.rm = T) * 365,
                            
                            NH3_kg_yr = mean(NH3_emis_rate, na.rm = T)/1000 * 365,  NH3_kg_anim_yr = mean(NH3_emis_rate, na.rm = T)/1000 * 365 /n_anim, NH3_kg_parea_yr = mean(NH3_emis_rate, na.rm = T)/1000 * 365 /prod_area,
                            NH3_kg_m3_ex_yr = mean(NH3_emis_rate, na.rm = T)/1000 * 365 /eval(.slurry_produced), NH3_kg_m3_year = mean(NH3_emis_rate, na.rm = T)/mean(slurry_mass, na.rm = T) * 365,
                            
                            N2O_kg_yr = mean(N2O_emis_rate, na.rm = T)/1000 * 365,  N2O_kg_anim_yr = mean(N2O_emis_rate, na.rm = T)/1000 * 365 /n_anim, N2O_kg_parea_yr = mean(N2O_emis_rate, na.rm = T)/1000 * 365 /prod_area,
                            N2O_kg_m3_ex_yr = mean(N2O_emis_rate, na.rm = T)/1000 * 365 /eval(.slurry_produced), N2O_kg_m3_year = mean(N2O_emis_rate, na.rm = T)/mean(slurry_mass, na.rm = T) * 365,
                            
                            
                            CO2_kg_yr = mean(CO2_emis_rate, na.rm = T)/1000 * 365, CO2_kg_anim_yr = mean(CO2_emis_rate, na.rm = T)/1000 * 365 /n_anim, CO2_kg_parea_yr = mean(CO2_emis_rate, na.rm = T)/1000 * 365 /prod_area,
                            CO2_kg_m3_ex_yr = mean(CO2_emis_rate, na.rm = T)/1000 * 365 /eval(.slurry_produced), CO2_kg_m3_year = mean(CO2_emis_rate, na.rm = T)/mean(slurry_mass, na.rm = T) * 365) %>% 
                            
                            filter (year == years-1) %>% distinct() %>% ungroup() %>% 
                            mutate(across(starts_with(c('CH4','NH3', 'N2O', 'CO2','slurry_')), function(x) round(x, digits = 2))) %>% 
    
                            bind_rows(summarise(select_if(., is.numeric), across(.fns = sum, na.rm = TRUE)))
  
  emission_dat$source[nrow(emission_dat)]  <- 'farm'                        


  nutrient_dat <- norm_dat %>% filter(source != 'enteric') %>% group_by(source, year, section_ID, storage_ID, digestate_ID) %>%
                  summarise(C_conc_in = sum(C_conc_fresh * slurry_prod_rate/1000)/eval(.slurry_produced), C_conc_out = mean(C_eff_conc, na.rm = T), C_conc_out = sum(C_eff_conc * slurry_mass_eff/1000, na.rm = T)/eval(.slurry_produced),
                            N_conc_in = sum(N_conc_fresh * slurry_prod_rate/1000)/eval(.slurry_produced), N_conc_out = sum(N_eff_conc * slurry_mass_eff/1000, na.rm = T)/eval(.slurry_produced),
                            TAN_conc_in = sum(conc_fresh_TAN * slurry_prod_rate/1000)/eval(.slurry_produced), TAN_conc_out = sum(TAN_eff_conc * slurry_mass_eff/1000, na.rm = T)/eval(.slurry_produced),
                            urea_conc_in = sum(conc_fresh_urea * slurry_prod_rate/1000)/eval(.slurry_produced), urea_conc_out = sum(urea_eff_conc * slurry_mass_eff/1000, na.rm = T)/eval(.slurry_produced)) %>%
                  mutate(across(everything(), function(x) round(x,2))) %>% filter(year == years-1)
  
  
  production_dat <- norm_dat %>%
    filter(source != 'enteric') %>%
    group_by(source, year, section_ID, storage_ID, digestate_ID) %>%
    summarise(slurry_produced = eval(.slurry_produced), rain = eval(.rain), 
              evaporation = eval(.evap), wash_water = eval(.wash_water)) %>%
    mutate(across(everything(), function(x) round(x, 1))) %>%
    filter(year == years - 1)
  
  
  explore_vars <- c('source', 'year', 'section_ID', 'storage_ID', 'digestate_ID',
                    'cum_inhib_m0', 'HAC_inhib_m0', 'H2S_inhib_m0', 'H2S_inhib_sr1', 'NH4_inhib_m0', 'NH3_inhib_m0', 'H2SO4_inhib_m0', 
                    'xa_fresh_m0', 'xa_fresh_m1', 'xa_fresh_m2', 'xa_fresh_sr1', 
                    'm0_conc', 'm1_conc', 'm2_conc',
                    'conc_fresh_urea', 'conc_fresh_sulfate', 'conc_fresh_sulfide', 
                    'conc_fresh_TAN', 'conc_fresh_VFA', 'conc_fresh_CP', 
                    'conc_fresh_CF', 'conc_fresh_RFd', 'conc_fresh_starch',
                    'COD_conc_fresh', 'dCOD_conc_fresh', 
                    'CH4_emis_rate_slurry', 'time')
  
  explore_dat <- bind_rows(barn_dat[, names(barn_dat) %in% explore_vars], storage_dat[, names(storage_dat) %in% explore_vars], digestate_dat[, names(digestate_dat) %in% explore_vars]) %>% 
    mutate(year = ceiling(time/365))
  
  expl_dat <- explore_dat %>% 
    group_by(source, year, section_ID, storage_ID, digestate_ID) %>%
    summarise(across(, .fns = ~ mean(.x, na.rm = TRUE))) %>% 
    mutate(across(everything(), function(x) round(x,2))) %>% 
    select(-time) %>% filter(year == years-1)
  
  
  return(list(emission_dat = emission_dat, nutrient_dat = nutrient_dat, production_dat = production_dat, expl_dat = expl_dat))
}  

