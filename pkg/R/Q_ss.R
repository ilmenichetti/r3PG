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
#'
#' @seealso \code{\link{C_init}}

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
#'
#' @seealso \code{\link{Q_ss}}
#'

u0_calc<-function(lat){
  u0_calc=(0.0855+0.0157*(50.6-0.768*lat))
  return(u0_calc)
}



#' @title Initializing the steady states of the decomposing pools (based on a simulated Picea abies stand)
#' @description This function utilizes the steady states solutions for the Q decomposition model  \code{\link{Q_ss}} to initialize the decomposing pools assuming that the simulation starts at their steady states, which are determined once the Q model parameters are known and we can assume steady inputs.
#'  It assumes, with some approximation, as steady inputs of C the average inputs from a whole management cycle.
#'  When initializing multiple layers, the function assigns the initialized C masses calculated from the initialization simulation proportionally to the number of stems of each layer.
#'
#' @param latitude self explanatory, it is used to initialize the Q model throught the function  refer to \code{\link{u0_calc}}
#' @param SI the site index
#' @param site refer to \code{\link{run_3PG}}
#' @param species  refer to \code{\link{run_3PG}}
#' @param climate  refer to \code{\link{run_3PG}}
#' @param thinning  refer to \code{\link{run_3PG}}
#' @param parameters  refer to \code{\link{run_3PG}}
#' @param parsQlitter  refer to \code{\link{run_3PG}}
#' @param soil  refer to \code{\link{run_3PG}}
#' @param size_dist  refer to \code{\link{run_3PG}}
#' @param settings  refer to \code{\link{run_3PG}}
#'
#' @details the function assumes the management plan from the function \code{\link{management_plan}}, which contains the references for standard Swedish clear-cut forestry practices
#' The function considers as simulation length the whole cycle indicated by \code{\link{management_plan}} plus one year, and utilizes the function \code{\link{extend_time}} to repeat the provided climatic time series over the needed time.
#' Because of this, the function assumes a pure stand of Picea abies as former management.
#'
#' @seealso \code{\link{management_plan}}, \code{\link{extend_time}}, \code{\link{u0_calc}}, , \code{\link{Q_ss}}


C_init<-function(lat,
                 SI,
                 site,
                 species,
                 climate,
                 thinning,
                 parameters,
                 parsQlitter,
                 size_dist,
                 settings = list(light_model = 2, transp_model = 2, phys_model = 2,
                                 height_model = 1, correct_bias = 0, calculate_d13c = 0)){

  #simulating a standard management plan (Picea abies)
  management<-management_plan(SI,lat)
  thinning_boolean<-management[8:11]!=1
  management[9]<-management[8]*management[9]
  management[10]<-management[9]*management[10]
  management[11]<-management[10]*management[11]

  main_species<-species[which.max(species$stems_n),] #legacy, when the simulation was based on the most represented species. It is there to initialize the matrix
  main_species$species<- "Picea abies"
  main_species$stems_n<-management$N
  main_thinning<-thinning[0,]
  for(i in 1:length(which(thinning_boolean))){
    if(thinning_boolean[i]==T){
      main_thinning_line<-thinning[1,][1,]
      main_thinning_line$species<- "Picea abies"
      main_thinning_line$age=management[2+i]
      main_thinning_line$stems_n=round(management[7+i]*management$N)
      main_thinning_line$stem=1
      main_thinning_line$root=1
      main_thinning_line$foliage=1
      main_thinning<-rbind(main_thinning, main_thinning_line)
    }
  }

  #adding the final felling to the management plan
  main_thinning_last_line<-thinning[thinning$species==as.character(main_species[1]),][1,]
  main_thinning_last_line$species<- "Picea abies"
  main_thinning_last_line$age<-management$gm
  main_thinning_last_line$stems_n<-0 #no stems are left
  main_thinning_last_line$stem=1
  main_thinning_last_line$root=1
  main_thinning_last_line$foliage=1
  main_thinning<-rbind(main_thinning, main_thinning_last_line)


  #selecting one dummy initialization value for simulating the first time series (to calculate the average inputs)
  soilInit<-c(3,3,3,3,0)


  #extending the cycle to one full length (management$gm+1 length, final felling)
  main_site<-site
  main_site$from<-main_species$planted
  end_cycle<-as.Date(paste(main_site$from, "-01", sep=""), format = "%Y-%m-%d") %m+% years(management$gm+1)
  main_site$to<-format(end_cycle, "%Y-%m")


  #extending the climate time series
  climate_long<-extend_time(climate=climate,from=main_site$from, length=management$gm+1)

  #the initialization assumes spruce parameterization
  Qpars<- list(7.00, 0.36, 0.25, 0.50, 0.00, 1.10, 0.50, 7.00, 0.36,
               0.25, 0.50, 0.00, 1.10, 0.50, 7.00, 0.36, 0.25, 0.50,
               0.00, 1.10, 0.50, 7.00, 0.36, 0.25, 0.50, 0.00, 1.10, 0.30)
  names(Qpars)<-(t(my_parsQlitter)[1,])

  my_parsQlitter<-data.frame(parameter=names(Qpars), "Picea abies"=as.vector(unlist(Qpars)))
  colnames(my_parsQlitter)[2]="Picea abies"

  init_out = run_3PG(
    site        = main_site,
    species     = main_species,
    climate     = climate_long,
    thinning    = main_thinning,
    parameters  = parameters,
    parsQlitter = my_parsQlitter,
    soil        = soilInit,
    size_dist   = size_dist,
    settings    = settings,
    check_input = TRUE,
    df_out = F)



  soilInit_out<-as.data.frame(mat.or.vec(1, 5))
  soilInit_out[1]<-Q_ss(b=as.numeric(Qpars$beta_fol),
                        e0=as.numeric(Qpars$e0_fol),
                        q0=as.numeric(Qpars$q0_fol),
                        u0=u0_calc(lat),
                        fc=as.numeric(Qpars$fc_fol),
                        eta11=as.numeric(Qpars$eta_11_fol),
                        Tmax=3,
                        inputs=mean(init_out[,1,4,7], na.rm = T))

  soilInit_out[2]<-Q_ss(b=as.numeric(Qpars$beta_root),
                        e0=as.numeric(Qpars$e0_root),
                        q0=as.numeric(Qpars$q0_root),
                        u0=u0_calc(lat),
                        fc=as.numeric(Qpars$fc_root),
                        eta11=as.numeric(Qpars$eta_11_root),
                        Tmax=3,
                        inputs=mean(init_out[,1,4,8], na.rm = T))

  soilInit_out[3]<-Q_ss(b=as.numeric(Qpars$beta_stem),
                        e0=as.numeric(Qpars$e0_stem),
                        q0=as.numeric(Qpars$q0_stem),
                        u0=u0_calc(lat),
                        fc=as.numeric(Qpars$fc_stem),
                        eta11=as.numeric(Qpars$eta_11_stem),
                        Tmax=mean(init_out[,1,2,5], na.rm=T),
                        inputs=mean(init_out[,1,4,12], na.rm = T))

  soilInit_out[4]<-Q_ss(b=as.numeric(Qpars$beta_bran),
                        e0=as.numeric(Qpars$e0_bran),
                        q0=as.numeric(Qpars$q0_bran),
                        u0=u0_calc(lat),
                        fc=as.numeric(Qpars$fc_bran),
                        eta11=as.numeric(Qpars$eta_11_bran),
                        Tmax=3,
                        inputs=mean(init_out[,1,4,13], na.rm = T))

  #this pool is not in use in the current model version
  soilInit_out[5]<-0



  #calculating the proportions to distribute the initialization
  species_proportions<-species$stems_n/sum(species$stems_n)
  species_soilInit_out<-mat.or.vec(dim(species)[1], 5)

  for(i in 1:dim(species_soilInit_out)[1]){
    species_soilInit_out[i,]<-as.numeric(species_proportions[i]*soilInit_out)
  }

  return(species_soilInit_out)

}




#main function for management plan selection
#' @title The function for determining a simulated management plan for a Picea abies pure stand.
#'
#' @description
#' The management_plan function contains a series of vectors with management data from Sweden, for a pure Picea abies stand.
#' It selects the most appropriate values based on site index and latitude.
#'
#' @param SI The site index.
#' @param Lat The latitude of the simulated site.
#'
#' @return a data frame with :
#' @H0 the start age, when the planted trees are 1.3 m
#' @N the density of the initial plantation
#' @g a series of numbers showing at which age the stand is thinned, including the last final felling
#' @s what ratio is left at the current stand at each thinning
#'
#' @seealso \code{\link{C_init}}

management_plan<-function(SI, Lat){

  ### Self-contained declaration of matrices for the management plans
  # this is in order to skip the reading of the specific csv and make the function more general

  N0<-c(1250,1250,1250,1250,1400,1600,1800,2000,2200,2300,2500,2700,2900,3000)

  gm<-c(141, 141, 141, 141, 126, 121, 111, 106, 101, 95, 85, 80, 75, 70)

  s1_vec<-c(0.75, 0.75, 0.75, 0.75, 0.72, 0.69, 0.71, 0.73, 0.76, 0.75, 0.72, 0.73, 0.76, 0.76)
  s2_vec<-c(0.71, 0.71, 0.71, 0.71, 0.74, 0.72, 0.73, 0.73, 0.73, 0.74, 0.73, 0.67, 0.71, 0.69)
  s3_vec<-c(1, 1, 1, 1, 1, 0.86, 0.84, 0.71, 0.74, 0.74, 0.76, 0.72, 0.74, 0.74)
  s4_vec<-c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.75, 0.72, 0.76, 0.77, 0.74)
  s<-as.data.frame(cbind(s1_vec,s2_vec,s3_vec,s4_vec))

  g1_vec<-c(75, 75, 75, 75, 62, 55, 47, 37, 33, 30, 27, 25, 24, 21)
  g2_vec<-c(100, 100, 100, 100, 82, 76, 65, 53, 48, 43, 36, 35, 35, 30)
  g3_vec<-c(140, 140, 140, 140, 125, 91, 85, 71, 65, 55, 45, 46, 45, 39)
  g4_vec<-c(140, 140, 140, 140, 125, 120, 110, 105, 100, 70, 57, 59, 57, 47)
  g<-as.data.frame(cbind(g1_vec,g2_vec,g3_vec,g4_vec))

  H0<-c(12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11)
  gm<-c(141, 141, 141, 141, 126, 121, 111, 106, 101, 95, 85, 80, 75, 70)


  #selection id
  ID<-((SI-10)/2)+1

  ### Select the values from the specific site

  #Year when start diameter at breast hight is measured
  T.=H0[ID]#was T


  #Year of final felling
  gm=gm[ID]
  #Thinnings
  g1=as.numeric(g[ID,1])
  g2=as.numeric(g[ID,2])
  g3=as.numeric(g[ID,3])
  g4=as.numeric(g[ID,4])
  s1=as.numeric(s[ID,1])
  s2=as.numeric(s[ID,2])
  s3=as.numeric(s[ID,3])
  s4=as.numeric(s[ID,4])

  N=N0[ID]

  thinnings<-data.frame(g1, g2, g3, g4, gm, s1, s2, s3, s4)
  results<-data.frame(H0=T., N, thinnings)

  return(results)
}




# function to repeat the climate time series
#' @title This function repeats the climate time series for a specified length of time. It is used for initialization.
#'
#' @description
#'
#' @param climate the original climatic time series.
#' @param from the start date of the desired new time series
#' @param length length in years of the desired new time series.
#'
#' @seealso \code{\link{C_init}}

extend_time<-function(climate, from, length){

  require(lubridate)
  length<-length+1

  from_date<-as.Date(paste(from, "-01", sep=""), format = "%Y-%m-%d")
  seq_date<-seq(from_date, by="month", length.out=length*12)

  repetitions<-ceiling(length/(diff(range(climate$year))+1))
  climate_long<-do.call(rbind, replicate(repetitions, climate, simplify=FALSE)) # where m is your matrix
  climate_long_cut<-climate_long[1:(length*12),]
  climate_long_cut$year<-year(seq_date)
  climate_long_cut$month<-month(seq_date)
  return(climate_long_cut)

}
