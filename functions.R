
pre.marginals.copula <- function(data= data) {
  
  # pseudo obs
  var_a <- pobs(data)[,2]
  var_b <- pobs(data)[,1]
  
  selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
  
  familyparam <- c(selectedCopula$family, selectedCopula$par, selectedCopula$par2)
}

