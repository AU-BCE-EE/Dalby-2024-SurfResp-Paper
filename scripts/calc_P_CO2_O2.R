library(biogas)

fpd <- 'C16 H27 O8.7 N' # average manure
molMass('CO2') * 16 / (calcCOD(fpd) * molMass(fpd))

fs0 <- 0.65 # from rittman

P_CO2_O2 <- (1-fs0) * molMass('CO2') * 16 / (calcCOD(fpd) * molMass(fpd))
