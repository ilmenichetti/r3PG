#' @title Q decomposition function
#' @description A new implementation of the Q model included in this package.
#' This implementation is derived from Menichetti et al. (2021)
#'
#' @param inputs the C inputs for the simulated cohorts
#' @param fc the C ratio of the input material
#' @param e0 the efficiency of the microbial biomass
#' @param eta_11
#' @param u0 the growth rate of the microbial biomass
#' @param q0 the initial quality of the C inputs coming to the cohort
#' @param delay the C inputs for the simulated cohorts
#' @param tmax the maximum diameter of the material particles or residuals
#' @param nSim time steps of the simulation
#'
#' @details this function is just a demo version of what was then implemented in the FORTRAN code run by the
#' wrapper \code{\link{run_3PG}}, and its purpose is to simulate the decomposition of the outputs from the growth module 3PG.
#'
#'
#' @example example here
#'
#' @references
#' Menichetti, L., Mäkinen, H., Stendahl, J., Ågren, G.I., Hyvönen, R., 2021. Modeling persistence of coarse woody debris residuals in boreal forests as an ecological property. Ecosphere 12. https://doi.org/10.1002/ecs2.3792

Q_dec<-function(inputs=inputs,
                fc,
                beta,
                e0,
                eta_11,
                u0,
                q0,
                delay,
                tmax,
                nSim){

    zeta = (1-e0)/(beta*eta_11*e0) #zeta is not variable over time, it is edaphic

  #alpha is variable over time, in theory, since u0 could later depend on climate too
  alpha = fc*beta*eta_11*u0*q0^beta #in this version alpha is calculated inside the loop because I want to have a different u0 for each time step for later on in the project, even if at this point it is not needed

  mass_left<-c()
  for(t in 1:nSim){ # loop running for each time step, calculating the ratio of decomposed inputs left
    #decomposition proportion function, runs again each timestep
    if(t<(tmax[t])){
      mass_left[t] = ((2/(tmax[t]))*(1/(alpha*(1-zeta)))*((1+alpha*t)^(1-zeta)-
                                                            (1-(t/(tmax[t]))))+
                        ((2/(tmax[t])^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*t)^(2-zeta)))+
                        (1-(t/(tmax[t])))^2)
    } else {
      mass_left[t] = (2/(tmax[t]))*(1/(alpha*(1-zeta)))*(1+alpha*t)^(1-zeta)+
        ((2/((tmax[t])^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(t-(tmax[t])))^(2-zeta))-
                                                              ((1+alpha*t)^(2-zeta))))
    }# end of the ifelse

  }

  # decomposing_matrix[j:length(inputs),j]<-inputs[j]*mass_left
  # time_dec_vec_sim=head(time_dec_vec_sim, -1) #here I am dropping the last element of the simulation vector. You can achieve the same in other ways...
  #

  Ct <-mass_left*inputs
  return(Ct)
}

u0_calc<-function(latitude){
  u0_calc=(0.0855+0.0157*(50.6-0.768*latitude))
}

Q_dec_m<-function(inputs,
                  fc,
                  beta,
                  e0,
                  eta_11,
                  u0,
                  q0,
                  delay,
                  tmax,
                  nSim){

  eta_11 = eta_11
  u0 = u0
  tm <- (1:nSim)/12 ####months as fraction of years
  #the matrix is needed to store each decomposing input flux, one per each time step. It is then summed at the end.

  zeta = (1-e0)/(beta*eta_11*e0) #zeta is not variable over time, it is edaphic

  # for(j in 1:length(inputs)){

  #alpha is variable over time, in theory, since u0 could later depend on climate too
  alpha = fc*beta*eta_11*u0*q0^beta #in this version alpha is calculated inside the loop because I want to have a different u0 for each time step for later on in the project, even if at this point it is not needed

  mass_left<-c()
  for(i in 1:nSim){ # loop running for each time step, calculating the ratio of decomposed inputs left
    #decomposition proportion function, runs again each timestep
    if(tm[i]<(tmax[i])){
      mass_left[i] = ((2/(tmax[i]))*(1/(alpha*(1-zeta)))*((1+alpha*tm[i])^(1-zeta)-
                                                            (1-(tm[i]/(tmax[i]))))+
                        ((2/(tmax[i])^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*tm[i])^(2-zeta)))+
                        (1-(tm[i]/(tmax[i])))^2)
    } else {
      mass_left[i] = (2/(tmax[i]))*(1/(alpha*(1-zeta)))*(1+alpha*tm[i])^(1-zeta)+
        ((2/((tmax[i])^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(tm[i]-(tmax[i])))^(2-zeta))-
                                                              ((1+alpha*tm[i])^(2-zeta))))
    }# end of the ifelse

  }

  # decomposing_matrix[j:length(inputs),j]<-inputs[j]*mass_left
  # time_dec_vec_sim=head(time_dec_vec_sim, -1) #here I am dropping the last element of the simulation vector. You can achieve the same in other ways...
  #

  Ct <-mass_left*inputs
  return(Ct)
}


Q_soc<-function(SOC_inputs,
                fc,
                beta,
                e0,
                eta_11,
                u0,
                q0,
                initSoil=0){

  nSim <- length(SOC_inputs)
  time_dec_vec=1:nSim
  #function describing the variation of SOC quality over time

  q_t = SOCt = rep(0,nSim)
  # SOC_inputs[1] <- SOC_inputs[1] + initSoil
  zeta = (1-e0)/(eta_11*e0) #zeta is not variable over time, it is edaphic
  if(initSoil>0){
    q_t <- 1/(1+beta*fc*eta_11*u0*q0^beta*time_dec_vec)^(1/beta)
    SOCt <- initSoil * q_t ^(zeta-beta)
  }
  # SOCt<-c()
  for(i in 1:length(time_dec_vec)){
    #function describing the SOC accumulation over time
    q_t[1:(nSim-i+1)] <- 1/(1+beta*fc*eta_11*u0*q0^beta*time_dec_vec[1:(nSim-i+1)])^(1/beta)
    SOCt[i:nSim]= SOCt[i:nSim] +
      (SOC_inputs[i]*(q_t[1:(nSim-i+1)]^zeta))
  }

  #calculating also the steady state while we're at it...
  SOC_ss= (mean(SOC_inputs)*e0)/(fc*u0*q0^beta)*
    1/(1-e0-beta*eta_11*e0)

  return(list(SOCt=SOCt, SOCss=SOC_ss))

}

Q_soc_m<-function(SOC_inputs,
                  fc,
                  beta,
                  e0,
                  eta_11,
                  u0,
                  q0,
                  initSoil=0){

  nSim <- length(SOC_inputs)
  tm <- (1:nSim)/12 ####months as fraction of years
  #the matrix is needed to store each decomposing input flux, one per each time step. It is then summed at the end.

  q_t = SOCt = rep(0,nSim)
  # SOC_inputs[1] <- SOC_inputs[1] + initSoil
  zeta = (1-e0)/(eta_11*e0) #zeta is not variable over time, it is edaphic

  if(initSoil>0){
    q_t <- 1/(1+beta*fc*eta_11*u0*q0^beta*tm)^(1/beta)
    SOCt <- initSoil * q_t ^(zeta-beta)
  }

  for(i in 1:nSim){
    #function describing the SOC accumulation over time
    q_t[1:(nSim-i+1)] <- 1/(1+beta*fc*eta_11*u0*q0^beta*tm[1:(nSim-i+1)])^(1/beta)
    SOCt[i:nSim]= SOCt[i:nSim] +
      (SOC_inputs[i]*(q_t[1:(nSim-i+1)]^zeta))
  }

  #calculating also the steady state while we're at it...
  SOC_ss= ((mean(SOC_inputs)*12)*e0)/(fc*u0*q0^beta)*
    1/(1-e0-beta*eta_11*e0)

  return(list(SOCt=SOCt, SOCss=SOC_ss))
}



#
# # Examples
# nSim=10
# nMsim=nSim*12
#
# ###inputs
# inputs = seq(from=3, to=5, length.out=nSim)
# inputsM <- rep(inputs/12,each=12)
# tmax = seq(from=1, to=2, length.out=nSim)
# tmaxM = rep(tmax,each=12)#seq(from=1, to=2, length.out=nMsim)
# inSOC = seq(from=3, to=3, length.out=nSim)
#
# ##parameters
# beta=7
# beta_SOC=7
# eta_11=0.36
# e0=0.25
# fc=0.5
# delay=0
# u0_om=u0_calc(55)
# u0_soc=0.8*u0_calc(55)
# q0_om=1.1
# q0_soc=0.7
#
#
#
# ###Run q_dec and monthly version (q_dec_m)
# Ct <- rep(0,nSim)
# for(i in 1:nSim){
#   Ct[i:nSim] <- Ct[i:nSim] + Q_dec(inputs = inputs[i],
#                                      tmax = tmax[1:(nSim-i+1)],
#                                      beta=7,
#                                      eta_11=0.36,
#                                      e0=0.25,
#                                      fc=0.5,
#                                      delay=0,
#                                      u0=u0_calc(55),
#                                      q0=1.1,
#                                      nSim=(nSim-i+1))
#
# }


# ###Run monthly version of q_dec (q_dec_m)
# Ct_m <- rep(0,nMsim)
# for(i in 1:nMsim){
#   Ct_m[i:nMsim] <- Ct_m[i:nMsim] + Q_dec_m(inputs = inputsM[i],
#                                            tmax = tmaxM[1:(nMsim-i+1)],
#                                            beta=7,
#                                            eta_11=0.36,
#                                            e0=0.25,
#                                            fc=0.5,
#                                            delay=0,
#                                            u0=u0_calc(55),
#                                            q0=1.1,
#                                            nSim=(nMsim-i+1))
#
# }
#
# plot(Ct_m)
# points((((1:nSim)*12)-6),Ct,pch=20,col=2)
#
#
# #####Soc function example
# nMsim=12
# nRep = 30
# inSOC2 <- rep(inSOC,nRep)
# inSOC_m <- rep(0,each=100)
# nSim <- length(inSOC2)
# # inSOC2[1] <- 1150
# inSOC_m[1] <- 3
# initSoil=0
# kk <- Q_soc(inSOC2, fc=fc, beta=beta_SOC, e0=e0, eta_11=eta_11,
#              u0=u0_soc, q0=q0_soc,initSoil = initSoil)
# kk_m <- Q_soc_m(inSOC_m, fc=fc, beta=beta_SOC, e0=e0, eta_11=eta_11,
#              u0=u0_soc, q0=q0_soc,initSoil = initSoil)
#
# ylim=range(kk$SOCt,kk_m$SOCt,kk$SOCss,kk_m$SOCss)
# plot(kk_m$SOCt,ylim=ylim,type="l")
# lines(seq(6, (length(inSOC2)*12), by=12),kk$SOCt,col=2,pch=20)
# abline(h=kk$SOCss,col=2)
# kk$SOCss
# kk_m$SOCss
#
#
