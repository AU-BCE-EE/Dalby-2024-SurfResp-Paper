masstransfer <- function(temp_C, temp_air_C, temp_K, temp_air_K, vent_air, vel_pit, RH){
  
# air velocities (m/s). ref Aarnink 1998
vel.pit <- vel_pit                                 
vel.floor <- 4.62 * vent_air/(24*60*60) + 0.106
  
vel <- c(vel.pit = vel.pit, vel.floor = vel.floor)   

# std constants
P <- 1            # pressure. atm
R.atm <- 0.082057 # gas constant. L*ATM/Kelvin/mol
R.bar <- 0.083145 # gas constant. L*bar/Kelvin/mol

# mol weight (g/mol)
MW1 <- c(H2S = 34.08, NH3 = 17.03, CO2 = 44.01, CH4 = 16.04, O2 = 32.0, HAC = 60.052, H2O = 18.016)
MW_air <- 28.964

# diffusion volumes (cm3/mol)
Vi1 <- c(H2S = 20.96, NH3 = 14.9, CO2 = 26.9, CH4 = 24.42, O2 = 16.6, HAC = 51.88, H2O = 12.7)
Vi_air <- 20.1

# height above surface of measured air velocity (m)
z <- 0.0005

# air velocities special (m/s)
U10 <- (10.4/(log(z) + 8.1)) * vel      # effective air velocity at 10 m height. ref Schwarsenbach 2003
U. <- 0.02 * U10 ^ 1.5                  # friction velocity. ref Rotz et al. 2012 IFSM

# dry matter NTS import from manure composition. 
dm <- 7 

# densities g/cm3
pw <- ((dm + 279)/0.28)/1000                                                                              # ref Thygesen 2012
pa <- 0.001 * (353/(273.16 + temp_air_C)) * ((760 - 0.3783 * RH * exp(0.0596 * temp_air_C - 1.6662))/760) # moist air

# dynamic (u in g/cm/s) and kinematic (v in cm2/s) viscosities
ua <- 1.8325*10^-4
uw <- exp(-52.843 + 3703.6/temp_K + 5.866*log(temp_K) -5.879*10^-(29) * temp_K^10) * 10 
va <- ua/pa 
vw <- uw/pw

# diffusion coefficient (cm2/s)
Dw <- 13.26*10^-5/((uw * 100)^1.14 * Vi1^0.589)
Da <- 10^-3 * ((273 + temp_C)^1.75 * (1/MW1 + 1/MW_air)^0.5)/(P * (Vi1^(1/3) + Vi_air^(1/3))^2)

# Schmidt number (dimless)
Sca <- va / Da 
Scw <- vw / Dw 

# mass transfer velocities (cm/s). ref blunden 2006
if(U.[1] < 0.3){
  kl.pit <- (1*10^-6 + 144*10^-4 * (U.[1])^2.2 * Scw^(-0.5)) * 100
} else{
  kl.pit <- (1*10^-6 + 34.1 * 10^-4 * (U.[1]) * Scw^(-0.5)) * 100 
}
ka.pit <- (10^-3 + 46.2 * 10^-3 * U.[1]* Sca^-0.67) * 100 

if(U.[2] < 0.3){
  kl.floor <- (1*10^-6 + 144*10^-4 * (U.[2])^2.2 * Scw^(-0.5)) * 100 
} else{
  kl.floor <- (1*10^-6 + 34.1 * 10^-4 * (U.[2]) * Scw^(-0.5)) * 100 
}
ka.floor <- (10^-3 + 46.2 * 10^-3 * U.[2]* Sca^-0.67) * 100

# Henrys constants ( mol/m3/Pa) 
H.cp <- c(H2S = 1*10^-3, NH3 = 5.9*10^-1, CO2 = 3.3*10^-4, CH4 = 1.4*10^-5, O2 = 1.3*10^-5, HAC = 40)
dH.cp <- c(H2S = 2100, NH3 = 4200, CO2 = 2400, CH4 = 1600, O2 = 1500, HAC = 6200)

# Henrys constant (dimless liquid/gas). ref beutier Renon 1978, Sander 2015, Blanes 2009, Aarnink 1998
H.H2S <- exp(403.658 - 7056.07/temp_K - 74.6926 * log(temp_K) + 0.14529* temp_K) * R.atm * temp_K
H.NH3 <- 1431 * 1.053^(293 - temp_K)
H.CO2 <- exp(-1082.37 + 34417.2/temp_K + 182.28 * log(temp_K) - 0.25159* temp_K) * R.atm * temp_K
H.CH4 <- H.cp[['CH4']] * 100 * exp(dH.cp[['CH4']] * (1/temp_K - 1/(298.16))) * R.bar * temp_K 
H.O2 <- H.cp[['O2']] * 100 * exp(dH.cp[['O2']] * (1/temp_K - 1/(298.16))) * R.bar * temp_K
H.HAC <- 10^(3.65 + 2596/temp_K) * R.atm * temp_K 

H <- c(H2S = H.H2S, NH3 = H.NH3, CO2 = H.CO2, CH4 = H.CH4, O2 = H.O2, HAC = H.HAC, H2O = NA)

# overall mass transfer coefficients (m/d)
conv <- 1/100 * 60 * 60 * 24                                  
K.ol.pit <- 1/(1/kl.pit + H/ka.pit) * conv  
K.ol.floor <- 1/(1/kl.floor + H/ka.floor) * conv  

return(c(list(K.ol.pit = K.ol.pit,  K.ol.floor = K.ol.floor, H = H, vel = vel)))
}
