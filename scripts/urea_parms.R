## constants for 25 deg C for urea hydrolysis in pig manure

Vmax <- 2.06 # mmol urea/kg wet feces/min
Km <- 32.59 # mmol urea/L
M_urea <- 60.06 # molar mass of urea
gN_gUrea <- 2* 14.01/M_urea # g N per g urea
# changing units
Vmax_gN_day <- Vmax/1000 * M_urea * gN_gUrea * 60* 24 * 0.25
Km_gN_L <- Km/1000 * M_urea * gN_gUrea 

## temperature sensitivty is modelled with the CTM model. 
# We have to define Vmax at the optimum temperature which we set to 50 deg C

tmax <- 60
tmin <- 0
topt <- 50
yopt_25 <- Vmax_gN_day
yopt <- 60 # guess this one such that yopt_25 == y
tt <- 25


y <- yopt * ((tt - tmax) * (tt - tmin)^2) / 
  ((topt - tmin) * ((topt - tmin) * (tt - topt) - 
                      (topt - tmax) * (topt + tmin - 2*tt))
  )


Km_gN_L

