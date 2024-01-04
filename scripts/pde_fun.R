pde_fun <- function(parms, x){
  
Length = x/10000 # cm
N = 200 
Grid <- setup.grid.1D(L = Length, N = N, dx.1 = 0.001)

diffO2_fun <- function(N){
  diffO2 <- parms[['D']] # cm^2/s
  return(diffO2)
}

diffO2_grid <- setup.prop.1D(diffO2_fun, grid = Grid)$int

names <- c('O2')
nspec <- length(names)
O2.ini <- rep(0, length = N)

state <- c(O2.ini)

model <- function(t, state, pars){
  
  with(as.list(pars), {
    O2 <- state[1:N] # extract state variables
    
    tran.O2 <- tran.1D(C = O2, C.up = kH * Cair, dx = Grid, D = diffO2_grid, v = 0) # calculate transport
    O2_resp <- r*O2/(KO2 + O2) # 
    dO2.dt <- tran.O2$dC - O2_resp # define the derivative as a combination of transport and reaction. 
    
    return(list(c(dO2.dt), O2_resp, O2_flux.up = tran.O2$flux.up))
  })
  
}

# get steady state solution
std <- steady.1D(y = state, func = model, parms = parms, nspec = nspec, 
                 dimens = N, names = names)

dat <- data.frame(depth = Grid$x.mid, conc = std$y, flux = std$O2_flux.up) %>% mutate(depth = depth * 10000, O2 = O2 * 1000) 
# above line converts cm to um for depth abd umol/cm^3 to umol/L

return(dat = dat)

}
