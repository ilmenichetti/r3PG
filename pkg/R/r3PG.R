#' @name r3PG
#'
#' @title Simulating Forest Growth using the 3-PG Process-Based Vegetation Model
#'
#' @docType package
#'
#' @description The r3PG package provides a flexible and easy-to-use interface for Fortran implementations of the 3-PGpjs (monospecific, evenaged and evergreen forests) or 3-PGmix (deciduous, uneven-aged or mixed-species forests) forest growth models. The user can flexibly switch between various options and submodules, to use the original 3-PGpjs model version for monospecific, even-aged and evergreen forests and the 3-PGmix model, which can also simulate multi-cohort stands (e.g. mixtures, uneven-aged) that contain deciduous species. The core function to run the model is \code{\link{run_3PG}}. For more background, please consult the vignette via vignette(package = "r3PG")
#'
#' @return  None
#'
#' @seealso \code{\link{run_3PG}}
#'
#' @example inst/examples/run_3PG-help.R
#'
#' @references
#' Forrester, D. I., 2020. 3-PG User Manual. Swiss Federal Institute for Forest, Snow and Landscape Research WSL, Birmensdorf, Switzerland. 70 p. Available at the following web site: \url{http://sites.google.com/site/davidforresterssite/home/projects/3PGmix/3pgmixdownload}
#'
#' Forrester, D. I., & Tang, X. (2016). Analysing the spatial and temporal dynamics of species interactions in mixed-species forests and the effects of stand density using the 3-PG model. Ecological Modelling, 319, 233–254. \doi{10.1016/j.ecolmodel.2015.07.010}
#'
#'Landsberg, J. J., & Waring, R. H., 1997. A generalised model of forest productivity using simplified concepts of radiation-use efficiency, carbon balance and partitioning. Forest Ecology and Management, 95(3), 209–228. \doi{10.1016/S0378-1127(97)00026-1}
#'
#'Sands, P. J., 2010. 3PGpjs user manual. Available at the following web site: \url{http://3pg.sites.olt.ubc.ca/files/2014/04/3PGpjs_UserManual.pdf}
NULL
