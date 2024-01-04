resCalc <- function(p, dat, weights = 1, to, plot = FALSE, ...){
  # Use ... to pass fixed parameters
  
  # Extract observations, keeping missing values (fixed interval needed to compare to measurements)
  obs <- obs.orig <- dat[, to] # umol/L 
  depth <- dat$depth # um
  # Get depth
  
  if (length(unique(diff(depth))) > 1) {
    warning('Measurement data frame dat does not have fixed interval.')
  }
  
  p <- 10^p
  pars <- c(parms, p)
  
  source('pde_fun.R')
  
  pred <- pde_fun(parms = pars, x = max(depth, na.rm = T)) %>% transmute(umol = O2, depth = depth)
  
  int_pred <- data.frame(approx(x = pred$depth, y = pred$umol, xout = depth)) %>% 
    transmute(umol = y, depth = x) %>% select(umol)
  
  rout <- which(is.na(int_pred))
  
  # Center and scale measurements and predicted values by *measured* mean sd
  #ss <- apply(obs, 2, sd, na.rm = TRUE)
  #sm <- apply(obs, 2, mean, na.rm = TRUE)
  #obs <- scale(obs, center = sm, scale = ss)
  #pred <- scale(pred, center = sm, scale = ss)
  
  res <- int_pred[-rout,] - obs[-rout,]
  
  obj <- sum(abs(res))
  
  cat('Parameter estimates: ', signif(as.numeric(p), 3), '\n')
  cat('Objective: ', signif(obj, 3), '\n')
  
  return(obj)
  
}