feedFun <- function(class_anim, type_anim, type, feed_dat, feed_intake, milk_prod, body_weight, batch_time){
  
 names.feces <- c('CP_feces', 'CF_feces', 'RFd_feces', 'iNDF_feces', 'starch_feces')
 
 if(class_anim == 'pig'){
 
   feed_DM_conc <- sum(feed_dat$DM * feed_dat$proportion)/1000
   DM_intake <- feed_DM_conc * feed_intake
   feed_dat$CP_feces <- feed_dat$DM * eval(parse(text = paste("feed_dat$Ex_CP", type, sep = ""))) / 1000 * feed_intake * feed_dat$proportion
   feed_dat$CF_feces <- feed_dat$DM * as.numeric(feed_dat$Ex_CF) / 1000 * feed_intake * feed_dat$proportion
   feed_dat$RFd_feces <- feed_dat$DM * eval(parse(text = paste("feed_dat$Ex_deg_ResFib",type, sep = ""))) / 1000 * feed_intake * feed_dat$proportion
   feed_dat$iNDF_feces <- feed_dat$DM * eval(parse(text = paste("feed_dat$Ex_nondeg_ResFib",type, sep = ""))) / 1000 * feed_intake * feed_dat$proportion
   feed_dat$starch_feces <- feed_dat$DM * feed_dat$Ex_Starch / 1000 * feed_intake * feed_dat$proportion
   feed_dat$OM_feces <- feed_dat$DM * eval(parse(text = paste("feed_dat$Ex_OM",type, sep = ""))) / 1000 * feed_intake * feed_dat$proportion
   feed_dat$OM_intake <- feed_dat$DM * as.numeric(feed_dat$OM) / 1000 * feed_intake * feed_dat$proportion
   
   feces_comp <- as.data.frame(t(colSums(feed_dat[, names.feces])))
   feces_comp$ash_feces <- sum(feces_comp)/0.90 - sum(feces_comp)
   feces_comp$VFA_feces <- sum(feces_comp) * 0.028 # VFA is 2.8 % of DM (n = 83). Source Saman 
   
   OM_dig <- (1 - sum(feed_dat$OM_feces)/sum(feed_dat$OM_intake))
   
   feces <- (5.405 - 6.31 * OM_dig + 0.505 * DM_intake/batch_time) * batch_time # kg feces per production period
   urine <- 2 * DM_intake # kg urine per production period
   urine_N <- (-21.20 + 0.134 * sum(feed_dat$CP * feed_dat$proportion) + 10.15 * (feed_intake/batch_time * feed_DM_conc)) * batch_time # g urine_N per production period
   urea_N <- urine_N * 0.75 # g urea_N per production period
   
   # H. Jørgensen 2007
   ent <- list(const = c(Grow = 15.2, Adult = 10.4), slope = c(Grow = 0.68, Adult = 0.91))
  
   CH4_ent <- (ent$const[[type]] + ent$slope[[type]] * # g CH4 per production period, 55.65 KJ / g CH4 
            (feces_comp$RFd_feces + feces_comp$iNDF_feces) / batch_time) * batch_time / 55.65 
   
  }
 
 if(class_anim == 'cattle'){
   
   DM_intake <- feed_intake  # kg year, since feed composition is given directly as DM. 
   feed_DM_conc <- 1/sum(feed_dat$proportion/(feed_dat$DM/10) * 100) # kg DM/ kg feed 
   
   #below is g/year
   
   feed_dat.mod <- feed_dat[, ! colnames(feed_dat) %in% c('DM', 'Feed_name', 'proportion')] *  DM_intake * feed_dat$proportion # g pr year
   feed_dat.mod <- data.frame(t(colSums(feed_dat.mod)))
   feed_dat.mod$OM_dig <- feed_dat.mod$OM * 0.73 
   feed_dat.mod$CF_dig <- (0.767 * feed_dat.mod$CF/batch_time - 6.6 * DM_intake/batch_time) * batch_time
   feed_dat.mod$CP_dig <- ((feed_dat.mod$N/batch_time * 0.96 - 8 * DM_intake/batch_time) * 6.25) * batch_time
   feed_dat.mod$starch_dig <- feed_dat.mod$Starch
   feed_dat.mod$sugar_dig <- feed_dat.mod$Sugar
   feed_dat.mod$ResFib_dig <- feed_dat.mod$OM_dig - feed_dat.mod$CF_dig - feed_dat.mod$CP_dig - feed_dat.mod$starch_dig - feed_dat.mod$sugar_dig
   
   feed_dat.mod$OM_feces <- feed_dat.mod$OM - feed_dat.mod$OM_dig
   feed_dat.mod$CP_feces <- (feed_dat.mod$N * 6.25) - feed_dat.mod$CP_dig
   feed_dat.mod$CF_feces <- feed_dat.mod$CF - feed_dat.mod$CF_dig
   feed_dat.mod$starch_feces <- feed_dat.mod$Starch - feed_dat.mod$starch_dig # sugar is assumed to be completely digested
   feed_dat.mod$RFd_feces <- feed_dat.mod$ResFib - feed_dat.mod$ResFib_dig - feed_dat.mod$iNDF
   feed_dat.mod$iNDF_feces  <- feed_dat.mod$iNDF
   
   feces_comp <- feed_dat.mod[,  colnames(feed_dat.mod) %in% names.feces]
   
   CH4_ent <- (76 + 13.5 * DM_intake/batch_time - 9.55 * feed_dat.mod$CF/DM_intake + 2.24 * feed_dat.mod$NDF/DM_intake) * batch_time
   
   if(grepl('heifer|cattle 0-6|cattle, tung race, dry|cattle, jersey, dry', type_anim)){
     feces <- (-0.3314 + 2.1819 * DM_intake/batch_time) * batch_time # kg/year
     dm_feces <- (0.1096 + 0.2849 * DM_intake/batch_time) * batch_time # kg/year
    # Johansen et al. Livestock Science 264 (2022) 105058
     urine <- (-16.4644 + 1.0469 * DM_intake/batch_time + 0.6972 * # L/year
       (feed_dat.mod$N / (DM_intake * 1000) * 100) + 27.7501 * 
       (feed_dat.mod$Sodium / (DM_intake * 1000) * 100) + 5.1286 *
       (feed_dat.mod$Potassium / (DM_intake * 1000) * 100)) * batch_time 
     # Johansen et al. Livestock Science 264 (2022) 105058
     urine_N <- (-181.56 + 13.9797 * DM_intake/batch_time + 75.9430 * # g/year
       (feed_dat.mod$N / (DM_intake * 1000) * 100)) * batch_time 
     # Johansen et al. Livestock Science 264 (2022) 105058
     urea_N <- 0.692 * urine_N # g/year
     # Møller
     feces_comp$VFA_feces <-  0.0287 * dm_feces * 1000
     
     #urine <- (3.3591 + 0.03001 * feed_dat.mod$N/batch_time) * batch_time 
     #urine_N <- (12 + 0.333 * feed_dat.mod$N/batch_time) * batch_time
     #urea_N <- (-54.4484 + 0.6640 * feed_dat.mod$N/batch_time) * batch_time
   }
   
   if(grepl('bull calf, 6-slaugter|bull calf, 0-6', type_anim)){
     feces <- (-0.6001 + 1.7689 * DM_intake/batch_time) * batch_time # kg/year
     dm_feces <- (-0.03894 + 0.2817 * DM_intake/batch_time) * batch_time # kg/year
     
     #P.A. Madsen et al. Livestock Science 267 (2023) 105139
     urine <- (-7.2850 + 1.3881 * DM_intake/batch_time - 8.6929 * # L/year
       (feed_dat.mod$N / (DM_intake * 1000) * 100) - 0.7668 * 
       (feed_dat.mod$Sodium / (DM_intake * 1000) * 100) + 14.4688 * 
       (feed_dat.mod$Potassium / (DM_intake * 1000) * 100)) * batch_time
     #Alternative equation, so urine_N is larger than urea_N
     urine_N <- (-95.69 + 10.87 * DM_intake/batch_time + 38.53 * # g/year
       (feed_dat.mod$N / (DM_intake * 1000) * 100)) * batch_time
     #P.A. Madsen et al. Livestock Science 267 (2023) 105139 data
     urea_N <- 0.694 * urine_N 
     #Møller et al. 
     feces_comp$VFA_feces <-  0.0287 * dm_feces * 1000
     
     #urine <- (3.8777 + 0.03446 * feed_dat.mod$N/batch_time) * batch_time ## old
     #urine_N <- (-95.69 + 10.87 * DM_intake/batch_time + 38.53 * # g/year
     #  (feed_dat.mod$N / (DM_intake * 1000) * 100)) * batch_time 
     #urine_N <- (12 + 0.333 * feed_dat.mod$N/batch_time) * batch_time
     #urea_N <- (-30.3995 + 0.5562 * feed_dat.mod$N/batch_time) * batch_time
   }
   
   if(grepl('cattle, tung race|cattle, jersey', type_anim) & !grepl("cattle, tung race, dry|cattle, jersey, dry", type_anim)){
     feces <- (-1.3548 + 2.2475 * DM_intake/batch_time) * batch_time # kg/year
     dm_feces <- (-0.1143 + 0.2805 * DM_intake/batch_time) * batch_time # kg/year
     
     # BANNINK ET AL. Journal of Dairy Science Vol. 82, No. 5, 1999
     urine <- (1.3441 * DM_intake/batch_time + (1.079 * # L/year
      (feed_dat.mod$Sodium / (DM_intake * 1000) * 100) + 0.5380 * 
      (feed_dat.mod$Potassium / (DM_intake * 1000) * 100) + 0.1266 *
      (feed_dat.mod$N / (DM_intake * 1000) * 100)) - milk_prod/batch_time * 
      (0.1216 + 0.0275 * milk_pro/10)) * batch_time # not working. gives negative 
     # NRC (2021) 
     urine_N <- (12 + 0.333 * feed_dat.mod$N/batch_time) * batch_time # g/year
     # from Dijkstra et al., 2013 Animal (2013), 7:s2, pp 292–302
     urea_N <- 0.761 * urine_N # g/year
     feces_comp$VFA_feces <-  0.0287 * dm_feces * 1000
     
     #urea_N <- (-33.2666 + 0.7891 * feed_dat.mod$N/batch_time -2.9549 * milk_prod/batch_time - 9.5565 * DM_intake/batch_time) * batch_time # g/year
     #urine <- (-10.6752 + 0.07641 * feed_dat.mod$N/batch_time) * batch_time # kg/year
     
   }
 
  }
 
 
  return(list(feces_comp = feces_comp, feces = feces, urine = urine, 
              feed_DM_conc = feed_DM_conc, urine_N = urine_N, urea_N = urea_N,
              CH4_ent = CH4_ent))
}