
pre.marginals.copula <- function(data= data) {
  
  # pseudo obs
  var_a <- pobs(data)[,2]
  var_b <- pobs(data)[,1]
  
  selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
  
  familyparam <- c(selectedCopula$family, selectedCopula$par, selectedCopula$par2)
}

generator <- function(draws=draws){
  
  w <- rMvdc(draws, my_dist_bb7)
  x <- rMvdc(draws, my_dist_60)
  y <- rMvdc(draws, my_dist_70_bb7)
  z <- rMvdc(draws, my_dist_80)
  
  x.draws <- list(w,x,y,z)
  
  return(x.draws)
}

nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}
