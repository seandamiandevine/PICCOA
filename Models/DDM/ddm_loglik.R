ddm_loglik <- function(pars, rt, response, drift, ndrift = 5) {
  
  if(any(is.na(names(pars)))) names(pars) = c("a" ,"t0",  "sv",  "sz",  "z",   "v_1", "v_2", "v_3", "v_4", "v_5")  # avoid some errors
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in 1:ndrift) {
    densities[drift == levels(drift)[i]] <- 
      ddiffusion(rt[drift == levels(drift)[i]], 
                 response = response[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"],  
                 sv=pars["sv"],
                 sz=if ("sz" %in% non_v_pars) pars["sz"] else 0.1,
                 z=if ("z" %in% non_v_pars) pars["z"]*pars["a"] else 0.5*pars["a"],
                 st0=if ("st0" %in% non_v_pars) pars["st0"] else 0, 
                 v=pars[base_par+i])
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
  
}