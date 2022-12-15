#' @title Steady state solution to the  Q decomposition function
#' @description
#'
#' @param inputs the C inputs for the simulated cohorts
#' @param fc the C ratio of the input material
#' @param e0 the efficiency of the microbial biomass
#' @param eta_11
#' @param u0 the growth rate of the microbial biomass
#' @param q0 the initial quality of the C inputs coming to the cohort
#' @param tmax the maximum diameter of the material particles or residuals
#'
#' @details The function is used to initialize the decomposing pools in the C stocks simulation
#'
#' @example
#' pars<-as.data.frame(t(i_parsQlitter[,c(4)]))
#' colnames(pars)<-t(i_parsQlitter[,c(1)])
#'
#'  #testing
#' Q_ss(fc=0.5,
#'      b=pars$beta_fol,
#'      e0=pars$e0_fol,
#'      q0=0.87, #pine material, Menichetti et al., 2021
#'      u0=u0_calc(45),
#'      eta11=pars$eta_11_fol,
#'      Tmax=3,
#'      I=3)
#'
#' @references
#' Menichetti, L., Mäkinen, H., Stendahl, J., Ågren, G.I., Hyvönen, R., 2021. Modeling persistence of coarse woody debris residuals in boreal forests as an ecological property. Ecosphere 12. https://doi.org/10.1002/ecs2.3792

Q_ss= function(fc, inputs, e0, eta11, b, u0, q0, Tmax){

  I=inputs
  a = fc*b*eta11*u0*q0^b
  z = (1-e0)/(b*eta11*e0) #zeta is not variable over time, it is edaphic

  C_ss=(1/a*1/(z-1)+Tmax/3)*I
  return(C_ss)
}


#' @title Scaling of the microbial activity term
#' @description
#'
#' @param lat latitude, self explanatory
#'
#' @details The function is used to rescale the microbial activity (and therefore the kinetic of decomposition)
#' based on latitude
#'
#' @references
#' Menichetti, L., Mäkinen, H., Stendahl, J., Ågren, G.I., Hyvönen, R., 2021. Modeling persistence of coarse woody debris residuals in boreal forests as an ecological property. Ecosphere 12. https://doi.org/10.1002/ecs2.3792


u0_calc<-function(lat){
  u0_calc=(0.0855+0.0157*(50.6-0.768*lat))
  return(u0_calc)
}



#' @title Initializing the steady states of the decomposing pools
#' @description
#'
#' @param latitude self explanatory
#' @param site,
#' @param species,
#' @param climate,
#' @param thinning,
#' @param parameters,
#' @param parsQlitter,
#' @param soil,
#' @param size_dist,
#' @param settings
#'
#' @details

C_init<-function(lat,
                 site,
                 species,
                 climate,
                 thinning,
                 parameters,
                 parsQlitter,
                 size_dist,
                 settings = list(light_model = 2, transp_model = 2, phys_model = 2,
                                 height_model = 1, correct_bias = 0, calculate_d13c = 0)){

  soilInit<-mat.or.vec(dim(my_species)[1], 5)
  for(i in 1:dim(soilInit)[1]){
    soilInit[i,]<-c(3,3,3,3,0)
  }


  init_out = run_3PG(
    site        = site,
    species     = species,
    climate     = climate,
    thinning    = thinning,
    parameters  = parameters,
    parsQlitter = parsQlitter,
    soil        = soilInit,
    size_dist   = size_dist,
    settings    = settings,
    check_input = TRUE, df_out = F)

  soilInit_out<-mat.or.vec(dim(my_species)[1], 5)
  Qpars<-as.data.frame(t(my_parsQlitter)[-1,])
  colnames(Qpars)<-(t(my_parsQlitter)[1,])

  for(i in 1:dim(soilInit)[1]){
    soilInit_out[i,1]<-Q_ss(b=as.numeric(Qpars[i,]$beta_fol),
                            e0=as.numeric(Qpars[i,]$e0_fol),
                            q0=as.numeric(Qpars[i,]$q0_fol),
                            u0=u0_calc(lat),
                            fc=as.numeric(Qpars[i,]$fc_fol),
                            eta11=as.numeric(Qpars[i,]$eta_11_fol),
                            Tmax=3,
                            inputs=mean(init_out[,1,4,7]))

    soilInit_out[i,2]<-Q_ss(b=as.numeric(Qpars[i,]$beta_root),
                            e0=as.numeric(Qpars[i,]$e0_root),
                            q0=as.numeric(Qpars[i,]$q0_root),
                            u0=u0_calc(lat),
                            fc=as.numeric(Qpars[i,]$fc_root),
                            eta11=as.numeric(Qpars[i,]$eta_11_root),
                            Tmax=3,
                            inputs=mean(init_out[,1,4,8]))

    soilInit_out[i,3]<-Q_ss(b=as.numeric(Qpars[i,]$beta_stem),
                            e0=as.numeric(Qpars[i,]$e0_stem),
                            q0=as.numeric(Qpars[i,]$q0_stem),
                            u0=u0_calc(lat),
                            fc=as.numeric(Qpars[i,]$fc_stem),
                            eta11=as.numeric(Qpars[i,]$eta_11_stem),
                            Tmax=mean(init_out[,i,2,5]),
                            inputs=mean(init_out[,i,4,12]))

    soilInit_out[i,4]<-Q_ss(b=as.numeric(Qpars[i,]$beta_bran),
                            e0=as.numeric(Qpars[i,]$e0_bran),
                            q0=as.numeric(Qpars[i,]$q0_bran),
                            u0=u0_calc(lat),
                            fc=as.numeric(Qpars[i,]$fc_bran),
                            eta11=as.numeric(Qpars[i,]$eta_11_bran),
                            Tmax=3,
                            inputs=mean(init_out[,i,4,13]))

    soilInit_out[i,5]<-0

  }

  return(soilInit_out)

}

