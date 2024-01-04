get_storages <- function(storages = storages){
  
  storage_parms <- list(ex_storage_depth =       storages[['storage_depth,_m']],
                        ex_storage_diam =        storages[['storage_diameter,_m']],
                        ex_storage_vol =         storages[['storage_volume,_m3']],
                        ex_storage_area =        (storages[['storage_diameter,_m']]/2)^2 * pi, # m2
                        storage_ID =             storages[['storageID']],
                        
                        storage_acid =           storages[["acidification_s"]],
                        storage_acid_dose =      storages[["acid_dose_s,_kg_m-3"]],
                        cover_s =                storages[["cover_s"]],
                        venting_s =              ifelse(storages[["controlled_ventilation_s"]] == 'none', FALSE, TRUE),
                        venting_eff_s =          storages[["methane_oxidized_s,_%_of_produced"]],
                        flaring_s =              ifelse(storages[["flaring_s"]] == 'none', FALSE, TRUE),
                        flaring_eff_s =          storages[["methane_flared_s,_%_of_produced"]],
                        app1s =                  storages[["fraction_used1_s"]],
                        app2s =                  storages[["fraction_used2_s"]],
                        app3s =                  storages[["fraction_used3_s"]],
                        app4s =                  storages[["fraction_used4_s"]],
                        app5s =                  storages[["fraction_used5_s"]],
                        app6s =                  storages[["fraction_used6_s"]],
                        app7s =                  storages[["fraction_used7_s"]],
                        app_t1s =                storages[["app_time1_s"]],
                        app_t2s =                storages[["app_time2_s"]],
                        app_t3s =                storages[["app_time3_s"]],
                        app_t4s =                storages[["app_time4_s"]],
                        app_t5s =                storages[["app_time5_s"]],
                        app_t6s =                storages[["app_time6_s"]],
                        app_t7s =                storages[["app_time7_s"]],
                        
                        # digestate for degassed manure
                        digestate_depth =        storages[['digestate_depth,_m']],
                        digestate_diam =         storages[['digestate_diameter,_m']],
                        digestate_vol =          storages[['digestate_volume,_m3']],
                        digestate_area =         (storages[['digestate_diameter,_m']]/2)^2 * pi, # m2
                        digestate_ID =           storages[['digestateID']],
                        
                        digestate_acid =         storages[["acidification_d"]],
                        digestate_acid_dose =    storages[["acid_dose_d,_kg_m-3"]],
                        cover_d =                storages[["cover_d"]],
                        venting_d =              ifelse(storages[["controlled_ventilation_d"]] == 'none', FALSE, TRUE),
                        venting_eff_d =          storages[["methane_oxidized_d,_%_of_produced"]],
                        flaring_d =              ifelse(storages[["flaring_d"]] == 'none', FALSE, TRUE),
                        flaring_eff_d =          storages[["methane_flared_d,_%_of_produced"]],
                        app1d =                  storages[["fraction_used1_d"]],
                        app2d =                  storages[["fraction_used2_d"]],
                        app3d =                  storages[["fraction_used3_d"]],
                        app4d =                  storages[["fraction_used4_d"]],
                        app5d =                  storages[["fraction_used5_d"]],
                        app6d =                  storages[["fraction_used6_d"]],
                        app7d =                  storages[["fraction_used7_d"]],
                        app_t1d =                storages[["app_time1_d"]],
                        app_t2d =                storages[["app_time2_d"]],
                        app_t3d =                storages[["app_time3_d"]],
                        app_t4d =                storages[["app_time4_d"]],
                        app_t5d =                storages[["app_time5_d"]],
                        app_t6d =                storages[["app_time6_d"]],
                        app_t7d =                storages[["app_time7_d"]],
                        
                        # Wheater
                        ex_storage_rain =        storages[["rain,_mm_day-1"]],
                        ex_storage_radiation =   storages[["radiation,_MJ_m-2_day-1"]],
                        
                        kl_NH3_storage =         storages[['MFC_NH3_storage']])
  
  return(storage_parms)
  
}