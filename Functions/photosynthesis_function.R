#######################################################################################
# GROWTH - PHOTOSYNTHESIS - LIGHT function with self shading
# Reference: Burd & Dunton 2001 - MEPS
# Species: Halodule wrightii
#######################################################################################
#' Jassby-Platt parameterization (Jassby & Platt 1976): gross production (P) mmol O2 g dry wt-1 h-1 at irradiance (I) mmol m-2 s-1
#' @param I. 
#' 
#' 
Photo.Light.Shade <- function(I., PL.max, Ik, Ca, B.max){ 
  PL.max*tanh(I./Ik)*(1-Ca/B.max) 
}


######################################################################################
# GROWTH - TEMPERATURE function: Adams et al. 2017 - Scientific Reports: sp. C. serrulata
#####################################################################################
# Yan and Hunt model (Yan & Hunt 1999): gross photosynthesis (P) in mg C g-1 DW h-1 at time T
Photo.Temp <- function(T., PT.max, T.opt, T.max){ 
  PT.max*((T.max-T.)/(T.max-T.opt)) * ((T./T.opt) ^ (T.opt/(T.max-T.opt))) 
}


#######################################################################################
# PHOTOSYNTHESIS MODEL - Combine light and temp into one model 
#######################################################################################
photosynthesis.T.I <- function(T., I., PT.max, T.opt, T.max, Ik, Ca, B.max){
  PL.max <- Photo.Temp(T., PT.max, T.opt, T.max)
  Photo.Light.Shade(I., PL.max, Ik, Ca, B.max)
}

