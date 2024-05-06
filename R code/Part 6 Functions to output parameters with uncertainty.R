#######################################################
##
## Using this process to estimate uncertainty...
##
#######################################################
D1load = readRDS("parameters/pyrethroid_uncertainty_LancetGH2024")
D2load = readRDS("parameters/pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024")
D3load = readRDS("parameters/pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024")

itn_standard_params = function(how_many_draws,  ## First choose the number of uncertainty runs you need
                               resistance_level## Next select the relevant level of resistance you need
){
  
  rng = sample(1:1000,how_many_draws,replace = FALSE)
  
  uncertainty_draws = D1load[which(D1load$resistance == resistance_level),]
  
  ## And translate mn_dur to gamman
  uncertainty_draws$gamman = uncertainty_draws$mean_duration/log(2)
  
  key_data = data.frame(dn0    = uncertainty_draws$dn0,
                        rn0    = uncertainty_draws$rn_pyr,
                        gamman = uncertainty_draws$gamman)
  
  return(key_data[rng,])
}

itn_pbo_params = function(how_many_draws,  ## First choose the number of uncertainty runs you need
                          resistance_level## Next select the relevant level of resistance you need
){
  
  rng = sample(1:1000,how_many_draws,replace = FALSE)
  
  relevant_data = D2load[which(D2load$resistance.x == resistance_level),]
  uncertainty_draws = relevant_data[rng,c(1,2,3,4,10)]
  
  ## And translate mn_dur to gamman
  uncertainty_draws$gamman_pbo = uncertainty_draws$mn_dur/log(2)
  
  key_data = data.frame(dn0    = uncertainty_draws$dn0,
                        rn0    = uncertainty_draws$rn_pbo,
                        gamman = uncertainty_draws$gamman_pbo)
  
  return(key_data)
}


itn_pyrrole_params = function(how_many_draws,  ## First choose the number of uncertainty runs you need
                          resistance_level## Next select the relevant level of resistance you need
){
  
  rng = sample(1:1000,how_many_draws,replace = FALSE)
  
  relevant_data = D3load[which(D3load$resistance.x == resistance_level),]
  uncertainty_draws = relevant_data[rng,c(1,2,3,4,10)]
  
  ## And translate mn_dur to gamman
  uncertainty_draws$gamman_pbo = uncertainty_draws$mn_dur/log(2)
  
  key_data = data.frame(dn0    = uncertainty_draws$dn0,
                        rn0    = uncertainty_draws$rn_pbo,
                        gamman = uncertainty_draws$gamman_pbo)
  
  return(key_data)
}

itn_standard_params(50,0.08)
itn_pbo_params(50,0.08)
itn_pyrrole_params(50,0.08)

###########################################################
##
## For median parameter estimates
##
#############################################################

median_ITN_params = function(which_ITNs,
                             resistance_level){
  
  median_ITN_parameters = expand.grid(nets = c("standard","PBO","pyrrole"))
  
  median_ITN_parameters$dn0 = c(median(itn_standard_params(1000,resistance_level)$dn0),
                                median(itn_pbo_params(1000,resistance_level)$dn0),
                                median(itn_pyrrole_params(1000,resistance_level)$dn0))
  
  median_ITN_parameters$rn0 = c(median(itn_standard_params(1000,resistance_level)$rn0),
                                median(itn_pbo_params(1000,resistance_level)$rn0),
                                median(itn_pyrrole_params(1000,resistance_level)$rn0))
  
  median_ITN_parameters$gamman = c(median(itn_standard_params(1000,resistance_level)$gamman),
                                   median(itn_pbo_params(1000,resistance_level)$gamman),
                                   median(itn_pyrrole_params(1000,resistance_level)$gamman))
  
  if(which_ITNs == c("standard"))
    return(c(median_ITN_parameters[1,])) else if(which_ITNs == c("PBO"))
      return(c(median_ITN_parameters[2,])) else if(which_ITNs == c("pyrrole"))
        return(c(median_ITN_parameters[3,])) else if (which_ITNs == "all")
          return(median_ITN_parameters)
}

median_ITN_params("all",)
