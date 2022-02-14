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




params_allometry = {
  "Marklund": {
    "mleaf": {
      'pine': {'a': 7.7681, 'b': 7, 'c': 3.7983},
      'spruce': {'a': 7.8171, 'b': 12, 'c': 1.9602},
      'birch': {'a': 8.058, 'b': 8, 'c': 3.9823} # Kellomäki et al. 2001 Atm. Env.
    },

    "mstemwood": {
      'pine': {'a': 11.4219, 'b': 14, 'c': 2.2184},
      'spruce': {'a': 11.4873, 'b': 14, 'c': 2.2471},
      'birch': {'a': 10.8109, 'b': 11, 'c': 2.3327}
    },

    "mstembark": {
      'pine': {'a': 8.8489, 'b': 16, 'c': 2.9748},
      'spruce': {'a': 9.8364, 'b': 15, 'c': 3.3912},
      'birch': {'a': 10.3876, 'b': 14, 'c': 3.2518}
    },

    "mlivingbranches": {
      'pine': {'a': 9.1015, 'b': 10, 'c': 2.8604}, # including needles
      'spruce': {'a': 8.5242, 'b': 13, 'c': 1.2804}, # including needles
      'birch': {'a': 10.2806, 'b': 10, 'c': 3.3633} # excluding leaves
    },

    "mstump": {
      'pine': {'a': 11.0481, 'b': 15, 'c': 3.9657},
      'spruce': {'a': 10.6686, 'b': 17, 'c': 3.3645},
    },

    "mroots": {
      'pine': {
        'a1': 13.2902, 'b1': 9, 'c1': 6.3413,
        'a2': 8.8795, 'b2': 10, 'c2': 3.8375
      },
      'spruce': {
        'a1': 13.3703, 'b1': 8, 'c1': 6.3851,
        'a2': 7.6283, 'b2': 12, 'c2': 2.5706
      }
    },
  },