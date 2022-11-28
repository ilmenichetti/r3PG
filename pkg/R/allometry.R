###define list of parameters
params_allometry <- list()
##
params_allometry$Marklund$mleaf$pine <- list(a=7.7681, b=7, c=3.7983)
params_allometry$Marklund$mleaf$spruce <- list(a = 7.8171, b = 12, c = 1.9602)
params_allometry$Marklund$mleaf$birch <- list(a=7.7681, b=7, c=3.7983)

###define list of functions
eqs_allometry <- list()
###function inputs: d=diameter at breast height; pars= parameters of the equation
eqs_allometry$Marklund$mleaf <- function(d,pars){
  a=pars$a
  b=pars$b
  c=pars$c
  mleaf <- exp(a*d/(d+b) - c)
  return(mleaf=mleaf)
}

eqType = "Marklund"
treeOrg <- "mleaf"
species <- "spruce"


###extract parameters
parX <- params_allometry[[eqType]][[treeOrg]][[species]]
###extract function
funX <- eqs_allometry[[eqType]][[treeOrg]]

wX <- funX(d=1:100,parX)
plot(wX,type='l')



