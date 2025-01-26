require(reticulate)
require(pracma)
require(fields)
require(plotrix)
require(stats)
require(keras)
require(nloptr)
require(FITSio)

py_path <- "/bin/python3.8"
# py_path <- "/bin/python3.10"
home_folder <- "/media/joaomfras/GalaxyPol/Pol-Gal/"

libs_folder <- paste0(home_folder, "Polarimetric-Imaging-Reduction-Scripts-(FO",
                      "RS2)/Commit/")

use_python(py_path, required = T)
rm(py_path)

mypy_lib <- paste0(libs_folder, "sexyDarce_2.py")
load_mod_cmd <- paste0("sexD = SourceFileLoader('sexyDarce_2', '", mypy_lib, 
                       "').load_module()")
py_run_string("from importlib.machinery import SourceFileLoader")
py_run_string("import numpy as np")
py_run_string("import sep")
py_run_string("import scipy")
py_run_string("from photutils import MedianBackground")
py_run_string("from photutils.segmentation import make_source_mask")
py_run_string("from photutils.isophote import EllipseGeometry")
py_run_string("from photutils.isophote import Ellipse")
py_run_string("from astropy.io import fits")
py_run_string("from astropy.stats import SigmaClip")
py_run_string("sig_Clip = SigmaClip(sigma=3.0)")
py_run_string("bkg_Est = MedianBackground()")
py_run_string(load_mod_cmd)

################################################################################
# Given an array of 3d "arr" (2d map with 2 layers, one for data and a second  #
# for uncertainties), and an array of binning steps "bin_l", bin_code_map will #
# return "bin_arr" an array holding the median and the propagated uncertainty  #
# of the binned elements of "arr". "bin_l" must be of length 1 or 2, "bin_arr" #
######################### is a compressed version of "arr". ####################
################################################################################
bin_code_map <- function(arr, bin_l, ref_pix = NULL, min_r = 0.67,             #
                         func = "median"){                                     #
                                                                               #
  dims <- dim(arr)                                                             #
  dnames <- dimnames(arr)                                                      #
  Ndim <- length(dims)                                                         #
  lref <- length(ref_pix)                                                      #
  Nbin <- length(bin_l)                                                        #
                                                                               #
  if(length(min_r) != 1 || min_r < 0 || min_r > 1){                            #
    print(paste0("ERROR: 'min_p' must be a length 1 numeric between 0 and 1. ",#
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(Ndim != 3 && dims[3] != 2){                                               #
    print(paste0("ERROR: 'arr' must be 3d array holding two 2d maps. Returnin",#
                 "g NULL."))                                                   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_pix' must be a length 2 or 3 vector. Proceedi",#
                 "ng with usage of 'arr' central pixel as 'ref_pix'."))        #
    ref_pix <- floor(dims[1:2] / 2)                                            #
  }else{                                                                       #
    ref_pix <- floor(ref_pix[1:2])                                             #
  }                                                                            #
                                                                               #
  if(Nbin == 1){                                                               #
    bin_l <- c(bin_l, bin_l)                                                   #
  }                                                                            #
  if(Nbin != 1 && Nbin != 2){                                                  #
    print(paste0("ERROR: Length of 'bin_l' must be either 1 or 2. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
  if(func != "median" & func != "mean"){                                       #
    print(paste0("ERROR: 'func' must be either 'median' or 'mean'. Returning ",#
                 "NULL."))                                                     #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  rows_f <- c(sort(seq(ref_pix[1] + floor(bin_l[1] / 2), 1, -bin_l[1])),       #
              seq(ref_pix[1] + floor(bin_l[1] / 2) + bin_l[1], dims[1],        #
                  bin_l[1]))                                                   #
  l_r <- length(rows_f)                                                        #
  rows_i <- c(1, rows_f[1:(l_r - 1)] + 1)                                      #
                                                                               #
  if(rows_f[l_r] > dims[1]){                                                   #
    rows_f[l_r] <- dims[1]                                                     #
  }                                                                            #
  if(rows_f[l_r] <  dims[1]){                                                  #
    rows_f <- append(rows_f, dims[1])                                          #
    rows_i <- append(rows_i, rows_f[l_r] + 1)                                  #
  }                                                                            #
  l_r <- length(rows_i)                                                        #
                                                                               #
  cols_f <- c(sort(seq(ref_pix[2] + floor(bin_l[2] / 2), 1, -bin_l[2])),       #
              seq(ref_pix[2] + floor(bin_l[2] / 2) + bin_l[2], dims[2],        #
                  bin_l[2]))                                                   #
  l_c <- length(cols_f)                                                        #
  cols_i <- c(1, cols_f[1:(l_c - 1)] + 1)                                      #
                                                                               #
  if(cols_f[l_c] > dims[2]){                                                   #
    cols_f[l_c] <- dims[2]                                                     #
  }                                                                            #
  if(cols_f[l_c] < dims[2]){                                                   #
    cols_f <- append(cols_f, dims[2])                                          #
    cols_i <- append(cols_i, cols_f[l_c] + 1)                                  #
  }                                                                            #
  l_c <- length(cols_i)                                                        #
                                                                               #
  bin_arr <- array(NA, dim = c(l_r, l_c, 2), dimnames = dnames)                #
  n_bin <- bin_l[1] * bin_l[2]                                                 #
                                                                               #
  for(r in 1:l_r){                                                             #
    ri <- rows_i[r]                                                            #
    rf <- rows_f[r]                                                            #
                                                                               #
    for(c in 1:l_c){                                                           #
      ci <- cols_i[c]                                                          #
      cf <- cols_f[c]                                                          #
                                                                               #
      temp_bin <- arr[ri:rf, ci:cf, 1]                                         #
      temp_bin_unc <- arr[ri:rf, ci:cf, 2]                                     #
                                                                               #
      temp_bin_NA <- which(is.na(temp_bin), arr.ind = T)                       #
      temp_bin_unc_NA <- which(is.na(temp_bin_unc), arr.ind = T)               #
                                                                               #
      temp_bin[temp_bin_unc_NA] <- NA                                          #
      temp_bin_unc[temp_bin_NA] <- NA                                          #
                                                                               #
      n_temp <- length(which(!is.na(temp_bin), arr.ind = T)) / 2               #
                                                                               #
      if(n_temp >= (min_r * n_bin)){                                           #
        if(func == "median"){                                                  #
          bin_arr[r, c, 1] <- median(temp_bin, na.rm = TRUE)                   #
          bin_arr[r, c, 2] <- unc_median(temp_bin, temp_bin_unc, 0)            #
        }else{                                                                 #
          bin_arr[r, c, 1] <- mean(temp_bin, na.rm = TRUE)                     #
          bin_arr[r, c, 2] <- unc_mean(temp_bin, temp_bin_unc, 0)              #
        }                                                                      #
      }                                                                        #
    }                                                                          #
  }                                                                            #
  return(bin_arr)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a binned 3d array "arr" (2 binned maps, data and uncertainty), and an ##
#### array of binning steps "bin_l", bin_shrink_map will return "bin_arr" a ####
# smaller array holding the same information "arr". "bin_l" must have length 1 #
########## or 2, "bin_arr" is the smaller compressed version of "arr". #########
################################################################################
bin_shrink_map <- function(arr, bin_l, ref_pix = NULL){                        #
                                                                               #
  dims <- dim(arr)                                                             #
  dnames <- dimnames(arr)                                                      #
  Ndim <- length(dims)                                                         #
  lref <- length(ref_pix)                                                      #
  Nbin <- length(bin_l)                                                        #
                                                                               #
  if(Ndim != 3 && dims[3] != 2){                                               #
    print(paste0("ERROR: 'arr' must 3d array holding two 2d maps. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_pix' must be a length 2 or 3 vector. Proceedi",#
                 "ng with usage of 'arr' central pixel as 'ref_pix'."))        #
    ref_pix <- floor(dims[1:2] / 2)                                            #
  }else{                                                                       #
    ref_pix <- floor(ref_pix[1:2])                                             #
  }                                                                            #
                                                                               #
  if(Nbin == 1){                                                               #
    bin_l <- c(bin_l, bin_l)                                                   #
  }                                                                            #
  if(Nbin != 1 && Nbin != 2){                                                  #
    print(paste0("ERROR: Length of 'bin_l' must be either 1 or 2. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  rows_f <- c(sort(seq(ref_pix[1] + floor(bin_l[1] / 2), 1, -bin_l[1])),       #
              seq(ref_pix[1] + floor(bin_l[1] / 2) + bin_l[1], dims[1],        #
                  bin_l[1]))                                                   #
  l_r <- length(rows_f)                                                        #
  rows_i <- c(1, rows_f[1:(l_r - 1)] + 1)                                      #
                                                                               #
  if(rows_f[l_r] > dims[1]){                                                   #
    rows_f[l_r] <- dims[1]                                                     #
  }                                                                            #
  if(rows_f[l_r] <  dims[1]){                                                  #
    rows_f <- append(rows_f, dims[1])                                          #
    rows_i <- append(rows_i, rows_f[l_r] + 1)                                  #
  }                                                                            #
  l_r <- length(rows_i)                                                        #
                                                                               #
  cols_f <- c(sort(seq(ref_pix[2] + floor(bin_l[2] / 2), 1, -bin_l[2])),       #
              seq(ref_pix[2] + floor(bin_l[2] / 2) + bin_l[2], dims[2],        #
                  bin_l[2]))                                                   #
  l_c <- length(cols_f)                                                        #
  cols_i <- c(1, cols_f[1:(l_c - 1)] + 1)                                      #
                                                                               #
  if(cols_f[l_c] > dims[2]){                                                   #
    cols_f[l_c] <- dims[2]                                                     #
  }                                                                            #
  if(cols_f[l_c] < dims[2]){                                                   #
    cols_f <- append(cols_f, dims[2])                                          #
    cols_i <- append(cols_i, cols_f[l_c] + 1)                                  #
  }                                                                            #
  l_c <- length(cols_i)                                                        #
                                                                               #
  bin_arr <- array(NA, dim = c(l_r, l_c, 2), dimnames = dnames)                #
                                                                               #
  for(r in 1:l_r){                                                             #
    rm <- (rows_i[r] + rows_f[r]) / 2                                          #
                                                                               #
    for(c in 1:l_c){                                                           #
      cm <- (cols_i[c] + cols_f[c]) / 2                                        #
                                                                               #
      bin_arr[r, c, 1] <- arr[rm, cm, 1]                                       #
      bin_arr[r, c, 2] <- arr[rm, cm, 2]                                       #
    }                                                                          #
  }                                                                            #
  return(bin_arr)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given a binned 3d array "arr" (2d map with 2 layers, one for data and  ####
# a second for uncertainties), an array of binning steps "bin_l", and a vector #
### with the dimensions of a target array "dims", bin_expand_map will return ###
# "big_arr" where each pixel of "arr" is expanded by "bin_l" to match "dims". ##
## If "dense" is TRUE, then all pixels within the bin will be filled, if it is #
##### FALSE, then only the central pixel of the bin will be filled and the #####
## remaining will be NAs. "bin_l" must be of length 1 or 2, "dims" must be of ##
################################## length 3. ###################################
################################################################################
bin_expand_map <- function(arr, bin_l, dims, ref_pix = NULL, dense = TRUE){    #
                                                                               #
  Ndim <- length(dims)                                                         #
  Nbin <- length(bin_l)                                                        #
  lref <- length(ref_pix)                                                      #
  dimbin <- dim(arr)                                                           #
  dnames <- dimnames(arr)                                                      #
                                                                               #
  if(Ndim != 3 && dims[3] != 2){                                               #
    print(paste0("ERROR: 'arr' must 3d array holding two 2d maps. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_pix' must be a length 2 or 3 vector. Proceedi",#
                 "ng with determination of 'ref_pix' as the central pixel of ",#
                 "'dims'."))                                                   #
    ref_pix <- floor(dims[1:2] / 2)                                            #
  }else{                                                                       #
    ref_pix <- floor(ref_pix[1:2])                                             #
  }                                                                            #
                                                                               #
  if(Nbin == 1){                                                               #
    bin_l <- c(bin_l, bin_l)                                                   #
  }                                                                            #
  if(Nbin != 1 && Nbin != 2){                                                  #
    print(paste0("ERROR: Length of 'bin_l' must be either 1 or 2. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  big_arr <- array(NA, dim = dims, dimnames = dnames)                          #
                                                                               #
  rows_f <- c(sort(seq(ref_pix[1] + floor(bin_l[1] / 2), 1, -bin_l[1])),       #
              seq(ref_pix[1] + floor(bin_l[1] / 2) + bin_l[1], dims[1],        #
                  bin_l[1]))                                                   #
  l_r <- length(rows_f)                                                        #
  rows_i <- c(1, rows_f[1:(l_r - 1)] + 1)                                      #
                                                                               #
  if(rows_f[l_r] > dims[1]){                                                   #
    rows_f[l_r] <- dims[1]                                                     #
  }                                                                            #
  if(rows_f[l_r] <  dims[1]){                                                  #
    rows_f <- append(rows_f, dims[1])                                          #
    rows_i <- append(rows_i, rows_f[l_r] + 1)                                  #
  }                                                                            #
                                                                               #
  cols_f <- c(sort(seq(ref_pix[2] + floor(bin_l[2] / 2), 1, -bin_l[2])),       #
              seq(ref_pix[2] + floor(bin_l[2] / 2) + bin_l[2], dims[2],        #
                  bin_l[2]))                                                   #
  l_c <- length(cols_f)                                                        #
  cols_i <- c(1, cols_f[1:(l_c - 1)] + 1)                                      #
                                                                               #
  if(cols_f[l_c] > dims[2]){                                                   #
    cols_f[l_c] <- dims[2]                                                     #
  }                                                                            #
  if(cols_f[l_c] < dims[2]){                                                   #
    cols_f <- append(cols_f, dims[2])                                          #
    cols_i <- append(cols_i, cols_f[l_c] + 1)                                  #
  }                                                                            #
                                                                               #
  if(dense){                                                                   #
    for(r in 1:dimbin[1]){                                                     #
      ri <- rows_i[r]                                                          #
      rf <- rows_f[r]                                                          #
                                                                               #
      for(c in 1:dimbin[2]){                                                   #
        ci <- cols_i[c]                                                        #
        cf <- cols_f[c]                                                        #
                                                                               #
        big_arr[ri:rf, ci:cf, 1] <- arr[r, c, 1]                               #
        big_arr[ri:rf, ci:cf, 2] <- arr[r, c, 2]                               #
      }                                                                        #
    }                                                                          #
  }else{                                                                       #
    for(r in 1:dimbin[1]){                                                     #
      rm <- floor(rows_i[r] + rows_f[r]) / 2                                   #
                                                                               #
      for(c in 1:dimbin[2]){                                                   #
        cm <- floor(cols_i[c] + cols_f[c]) / 2                                 #
                                                                               #
        big_arr[rm, cm, 1] <- arr[r, c, 1]                                     #
        big_arr[rm, cm, 2] <- arr[r, c, 2]                                     #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(big_arr)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects all variables to either have the same dimensions if  ##
## they are N-Dimensional or, if one is not N-Dimension, to have length 1 if  ##
###### they are 1-Dimensional. The function calculates the uncertainty of ######
####### dividing "idend" by "isor" given their respective uncertainties. #######
############### PS: Tests on variables length to be implemented. ###############
################################################################################
unc_div <- function(idend, isor, unc_idend, unc_isor, quot = idend / isor){    #
                                                                               #
  return(abs(quot) * sqrt((unc_idend / idend)^2 + (unc_isor / isor)^2))        #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects all variables to either have the same dimensions if  ##
#### they are N-Dimensional, or to have length 1 if they are 1-Dimensional. ####
## The function calculates the uncertainty of multiplying "m1" by "m2" given  ##
######################## their respective uncertainties. #######################
############### PS: Tests on variables length to be implemented. ###############
################################################################################
unc_mult <- function(m1, m2, unc_m1, unc_m2){                                  #
                                                                               #
  return(sqrt((m2 * unc_m1)^2 + (m1 * unc_m2)^2))                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects all variables to either have the same dimensions if  ##
#### they are N-Dimensional, or to have length 1 if they are 1-Dimensional. ####
### The function propagates the uncertainty of an addition between two terms ###
##################### given their respective uncertainties. ####################
################################################################################
unc_add <- function(unc_add1, unc_add2){                                       #
                                                                               #
  #Check if dims are compatible                                                #
  duadd1 <- dim(unc_add1)                                                      #
  luadd1 <- length(duadd1)                                                     #
  if(is.null(duadd1)){                                                         #
    duadd1 <- length(unc_add1)                                                 #
    luadd1 <- 1                                                                #
  }                                                                            #
                                                                               #
  duadd2 <- dim(unc_add2)                                                      #
  luadd2 <- length(duadd2)                                                     #
  if(is.null(duadd2)){                                                         #
    duadd2 <- length(unc_add2)                                                 #
    luadd2 <- 1                                                                #
  }                                                                            #
                                                                               #
  bool_d <- F                                                                  #
                                                                               #
  if((luadd1 == 1 & duadd1[1] == 1) | (luadd2 == 1 & duadd2[1] == 1)){         #
    bool_d <- T                                                                #
  }else{                                                                       #
    if(luadd2 == luadd1 & prod(duadd2 == duadd1)){                             #
      bool_d <- T                                                              #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(!bool_d){                                                                 #
    print(paste0("ERROR: the uncertainties of ADD1 and of ADD2 must either ha",#
                 "ve 1 dimension or length 1 or the same as each other."))     #
    print(paste0("Uncertainty of ADD1 has ", luadd1, " dims, uncertainty of A",#
                 "DD2 has ", luadd2, " dims."))                                #
    print(paste0("Uncertainty of ADD1 has dims ", paste(duadd1, collapse=" "), #
                 ", ", "uncertainty of ADD2 has dims ",                        #
                 paste(duadd2, collapse=" "), "."))                            #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  return(sqrt(unc_add1^2 + unc_add2^2))                                        #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Expects a numerical array "data". Outputs the medians standard deviation ###
########################### estimation of that array. ##########################
################################################################################
sd_median <- function(data){                                                   #
                                                                               #
  N <- length(which(!is.na(data), arr.ind = T))                                #
                                                                               #
  if(N > 1){                                                                   #
    return(mad(data, na.rm = T) * sqrt(pi / (2 * (N - 1))))                    #
  }else{                                                                       #
    if(N == 1){                                                                #
      return(0)                                                                #
    }else{                                                                     #
      return(NaN)                                                              #
    }                                                                          #
  }                                                                            #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
###### Expects a numerical array "data" holding values, a numerical array ######
#### "unc_data" holding the uncertainty of the previous arrays values, and #####
### another numerical array "dims2med" with the margins along which to work ####
## (see the documentation on the 'apply()' function). Outputs the uncertainty ##
######## of the median of "data" calculate along the margins "dims2med". #######
################################################################################
unc_median <- function(data, unc_data, dims2med = c(1,2)){                     #
                                                                               #
  #Check if dims are compatible                                                #
  dd <- dim(data)                                                              #
  if(is.null(dd)){                                                             #
    dd <- length(data)                                                         #
  }                                                                            #
  ld <- length(dd)                                                             #
                                                                               #
  du <- dim(unc_data)                                                          #
  if(is.null(du)){                                                             #
    du <- length(unc_data)                                                     #
  }                                                                            #
  lu <- length(du)                                                             #
                                                                               #
  ldm <- length(dims2med)                                                      #
                                                                               #
  all_flag <- ldm == 0 || (ldm == 1 && dims2med == 0)                          #
                                                                               #
  bool_d <- 0                                                                  #
  bool_l <- 0                                                                  #
  bool_ldm <- 0                                                                #
                                                                               #
  if(lu == ld | lu == 1){                                                      #
    bool_l <- 1                                                                #
  }                                                                            #
                                                                               #
  if(lu == 1){                                                                 #
    if(du == 1){                                                               #
      bool_d <- 1                                                              #
    }else{                                                                     #
      if(prod(du == dd))                                                       #
        bool_d <- 1                                                            #
    }                                                                          #
  }else{                                                                       #
    if(prod(du == dd))                                                         #
      bool_d <- 1                                                              #
  }                                                                            #
                                                                               #
  if((ldm > 0 && ldm <= ld && prod(0 < dims2med) && prod(ld >= dims2med)) ||   #
     (ldm == 1 & dims2med == 0))                                               #
    bool_ldm <- 1                                                              #
                                                                               #
  if(!bool_l){                                                                 #
    print(paste0("ERROR: the uncertainty must either have 1 dimension or the ",#
                 "same as the data."))                                         #
    print(paste0("Uncertainty has ", lu, " dims, data has ", ld, " dims."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!bool_d){                                                                 #
    print(paste0("ERROR: the uncertainty must have the same length as the dat",#
                 "a in all dimensions, or length 1 if it is one dimensional."))#
    print(paste0("Uncertainty has length ", paste(du, collapse = " "), ", dat",#
                 "a has length ", paste(dd, collapse = " "), "."))             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!bool_ldm){                                                               #
    print(paste0("ERROR: the dimensions across which you want to estimate the",#
                 " uncertainty of the median must match the dimensions along ",#
                 "which the median was calculated, and as such must be compat",#
                 "ible with the dimensions of the data."))                     #
    print(paste0("You state you want to estimate the uncertainty along dimens",#
                 "ions ", paste(dims2med, collapse = " "), ", but the data ha",#
                 "s ", ld, " dims them being: ", paste(dd, collapse = " "),    #
                 "."))                                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(all_flag){                                                                #
    medsd <- sd_median(data)                                                   #
  }else{                                                                       #
    medsd <- apply(data, dims2med, sd_median)                                  #
  }                                                                            #
                                                                               #
  if(lu == ld){                                                                #
    if(all_flag){                                                              #
      NU <- length(which(!is.na(unc_data), arr.ind = T))                       #
                                                                               #
      if(NU > 1){                                                              #
        unc_rep <- median(unc_data, na.rm=T) * sqrt(pi / (2 * (NU - 1)))       #
      }else{                                                                   #
        if(NU == 1){                                                           #
          unc_rep <- unc_data[which(!is.na(unc_data))]                         #
        }else{                                                                 #
          unc_rep <- NaN                                                       #
        }                                                                      #
      }                                                                        #
    }else{                                                                     #
      NU <- apply(unc_data, dims2med,                                          #
                  function(x) length(which(!is.na(x), arr.ind = T)))           #
                                                                               #
      unc_rep <- apply(unc_data, dims2med,                                     #
                       function(x){                                            #
                         NU <- length(which(!is.na(x), arr.ind = T))           #
                         if(NU > 1){                                           #
                           return(median(x, na.rm = T) *                       #
                                    sqrt(pi / (2 * (NU - 1))))                 #
                         }else{                                                #
                           if(NU == 1){                                        #
                             return(x[which(!is.na(x), arr.ind = T)])          #
                           }else{                                              #
                             return(NA)                                        #
                           }                                                   #
                         }                                                     #
                       })                                                      #
    }                                                                          #
  }else{                                                                       #
    unc_rep <- unc_data                                                        #
  }                                                                            #
                                                                               #
  return(sqrt(unc_rep^2 + medsd^2))                                            #
}                                                                              #
################################################################################
#------------------------------------------------------------------------------#

################################################################################
###### Expects a numerical array "data" holding values, a numerical array ######
#### "unc_data" holding the uncertainty of the previous arrays values, and #####
### another numerical array "dims2med" with the margins along which to work ####
## (see the documentation on the 'apply()' function). Outputs the uncertainty ##
######## of the mean of "data" calculated along the margins "dims2med". ########
################################################################################
unc_mean <- function(data, unc_data, dims2med = c(1,2)){                       #
                                                                               #
  #Check if dims are compatible                                                #
  dd <- dim(data)                                                              #
  if(is.null(dd)){                                                             #
    dd <- length(data)                                                         #
  }                                                                            #
  ld <- length(dd)                                                             #
                                                                               #
  du <- dim(unc_data)                                                          #
  if(is.null(du)){                                                             #
    du <- length(unc_data)                                                     #
  }                                                                            #
  lu <- length(du)                                                             #
                                                                               #
  ldm <- length(dims2med)                                                      #
                                                                               #
  all_flag <- ldm == 0 || (ldm == 1 && dims2med == 0)                          #
                                                                               #
  bool_d <- 0                                                                  #
  bool_l <- 0                                                                  #
  bool_ldm <- 0                                                                #
                                                                               #
  if(lu == ld | lu == 1){                                                      #
    bool_l <- 1                                                                #
  }                                                                            #
                                                                               #
  if(lu == 1){                                                                 #
    if(du == 1){                                                               #
      bool_d <- 1                                                              #
    }else{                                                                     #
      if(prod(du == dd))                                                       #
        bool_d <- 1                                                            #
    }                                                                          #
  }else{                                                                       #
    if(prod(du == dd))                                                         #
      bool_d <- 1                                                              #
  }                                                                            #
                                                                               #
  if((ldm > 0 && ldm <= ld && prod(0 < dims2med) && prod(ld >= dims2med)) ||   #
     (ldm == 1 & dims2med == 0))                                               #
    bool_ldm <- 1                                                              #
                                                                               #
  if(!bool_l){                                                                 #
    print(paste0("ERROR: the uncertainty must either have 1 dimension or the ",#
                 "same as the data."))                                         #
    print(paste0("Uncertainty has ", lu, " dims, data has ", ld, " dims."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!bool_d){                                                                 #
    print(paste0("ERROR: the uncertainty must have the same length as the dat",#
                 "a in all dimensions, or length 1 if it is one dimensional."))#
    print(paste0("Uncertainty has length ", paste(du, collapse = " "), ", dat",#
                 "a has length ", paste(dd, collapse = " "), "."))             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!bool_ldm){                                                               #
    print(paste0("ERROR: the dimensions across which you want to estimate the",#
                 " uncertainty of the median must match the dimensions along ",#
                 "which the median was calculated, and as such must be compat",#
                 "ible with the dimensions of the data."))                     #
    print(paste0("You state you want to estimate the uncertainty along dimens",#
                 "ions ", paste(dims2med, collapse = " "), ", but the data ha",#
                 "s ", ld, " dims them being: ", paste(dd, collapse = " "),    #
                 "."))                                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lu == ld){                                                                #
    if(all_flag){                                                              #
      NU <- length(which(!is.na(unc_data), arr.ind = T))                       #
                                                                               #
      if(NU > 1){                                                              #
        unc_rep <- sqrt(sum(unc_data^2, na.rm = T)) / NU                       #
      }else{                                                                   #
        if(NU == 1){                                                           #
          unc_rep <- unc_data[which(!is.na(unc_data))]                         #
        }else{                                                                 #
          unc_rep <- NaN                                                       #
        }                                                                      #
      }                                                                        #
    }else{                                                                     #
      unc_rep <- apply(unc_data, dims2med,                                     #
                       function(x){                                            #
                         NU <- length(which(!is.na(x), arr.ind = T))           #
                         if(NU > 1){                                           #
                           return(sqrt(sum(x^2, na.rm = T)) / NU)              #
                         }else{                                                #
                           if(NU == 1){                                        #
                             return(x[which(!is.na(x), arr.ind = T)])          #
                           }else{                                              #
                             return(NA)                                        #
                           }                                                   #
                         }                                                     #
                       })                                                      #
    }                                                                          #
  }else{                                                                       #
    unc_rep <- unc_data                                                        #
  }                                                                            #
                                                                               #
  return(unc_rep)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# This function requires that the variables have the same dimensions, it also  #
#### requires that the last dimension of each variable has length 2. First  ####
# position refers to HWP angle 0º and the second to HWP angle 45º if you wish  #
# to calculate Q; if you wish to calculate U, then the first position of each  #
# variables last dimension should refer to HWP angle 22.5º, and the second to  #
##### HWP angle 67.5º. The function returns an array with the same spatial #####
##### dimensions and an extra dimension, first position holding the Stokes #####
############### parameter value, and the second its uncertainty. ###############
################################################################################
stokesPar_from_dualBeams_Bagnulo <- function(Os, Es, unc_Os, unc_Es){          #
                                                                               #
  dOs <- dim(Os)                                                               #
  ldOs <- length(dOs)                                                          #
                                                                               #
  dEs <- dim(Es)                                                               #
  duOs <- dim(unc_Os)                                                          #
  duEs <- dim(unc_Es)                                                          #
                                                                               #
  null_flag <- is.null(dOs)                                                    #
                                                                               #
  if(null_flag & ldOs == 0){                                                   #
    ldOs <- length(Os)                                                         #
  }                                                                            #
                                                                               #
  if(null_flag & ldOs != 0){                                                   #
    if(ldOs == 3){                                                             #
      dOs <- c(1,1,2)                                                          #
      dEs <- c(1,1,2)                                                          #
      duOs <- c(1,1,2)                                                         #
      duEs <- c(1,1,2)                                                         #
    }else{                                                                     #
      dOs <- c(1,2)                                                            #
      dEs <- c(1,2)                                                            #
      duOs <- c(1,2)                                                           #
      duEs <- c(1,2)                                                           #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(dOs[ldOs] != 2 || !prod(dOs == dEs) || !prod(dOs == duOs) ||              #
     !prod(dOs == duEs) || ldOs > 3 || ldOs < 2){                              #
    cat("All variables must have the same dimensions (either 2 or 3), the las",#
        "t dimension must have length 2. Os is (")                             #
    for(l in ldOs){                                                            #
      cat(dOs[l])                                                              #
      if(l < ldOs){                                                            #
        cat(paste0(", "))                                                      #
      }                                                                        #
    }                                                                          #
    cat("), Es is (")                                                          #
    for(l in ldOs){                                                            #
      cat(dEs[l])                                                              #
      if(l < ldOs){                                                            #
        cat(paste0(", "))                                                      #
      }                                                                        #
    }                                                                          #
    cat("), unc_Os is (")                                                      #
    for(l in ldOs){                                                            #
      cat(duOs[l])                                                             #
      if(l < ldOs){                                                            #
        cat(paste0(", "))                                                      #
      }                                                                        #
    }                                                                          #
    cat(") and unc_Es is (")                                                   #
    for(l in ldOs){                                                            #
      cat(duEs[l])                                                             #
      if(l < ldOs){                                                            #
        cat(paste0(", "))                                                      #
      }                                                                        #
    }                                                                          #
    cat("). Returning NULL.")                                                  #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dnames <- list(NULL)                                                         #
                                                                               #
  if(ldOs == 3){                                                               #
      dnames <- append(dnames, list(NULL))                                     #
  }                                                                            #
                                                                               #
  dnames <- append(dnames, list(c("Data", "Unc")))                             #
                                                                               #
  stokesPar <- array(NA, dim = c(dOs[-ldOs], 2), dimnames = dnames)            #
                                                                               #
  OiE <- Os / Es                                                               #
  unc_OiE <- unc_div(Os, Es, unc_Os, unc_Es, OiE)                              #
                                                                               #
  if(ldOs == 3){                                                               #
    if(null_flag){                                                             #
      t1 <- OiE[1] / OiE[2]                                                    #
      unc_t1 <- unc_div(OiE[1], OiE[2], unc_OiE[1], unc_OiE[2], t1)            #
                                                                               #
      stokesPar[1,1,"Data"] <- (sqrt(t1) - 1) / (sqrt(t1) + 1)                 #
      stokesPar[1,1,"Unc"] <- abs(unc_t1 / ((sqrt(t1) + 1)^2 * sqrt(t1)))      #
    }else{                                                                     #
      t1 <- OiE[,,1] / OiE[,,2]                                                #
      unc_t1 <- unc_div(OiE[,,1], OiE[,,2], unc_OiE[,,1], unc_OiE[,,2], t1)    #
                                                                               #
      stokesPar[,,"Data"] <- (sqrt(t1) - 1) / (sqrt(t1) + 1)                   #
      stokesPar[,,"Unc"] <- abs(unc_t1 / ((sqrt(t1) + 1)^2 * sqrt(t1)))        #
    }                                                                          #
  }else{                                                                       #
    if(null_flag){                                                             #
      t1 <- OiE[1] / OiE[2]                                                    #
      unc_t1 <- unc_div(OiE[1], OiE[2], unc_OiE[1], unc_OiE[2], t1)            #
                                                                               #
      stokesPar[1,"Data"] <- (sqrt(t1) - 1) / (sqrt(t1) + 1)                   #
      stokesPar[1,"Unc"] <- abs(unc_t1 / ((sqrt(t1) + 1)^2 * sqrt(t1)))        #
    }else{                                                                     #
      t1 <- OiE[,1] / OiE[,2]                                                  #
      unc_t1 <- unc_div(OiE[,1], OiE[,2], unc_OiE[,1], unc_OiE[,2], t1)        #
                                                                               #
      stokesPar[,"Data"] <- (sqrt(t1) - 1) / (sqrt(t1) + 1)                    #
      stokesPar[,"Unc"] <- abs(unc_t1 / ((sqrt(t1) + 1)^2 * sqrt(t1)))         #
    }                                                                          #
  }                                                                            #
                                                                               #
  stokesPar[which(is.infinite(stokesPar), arr.ind = T)] <- NA                  #
  return(stokesPar)                                                            #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
##### This function assumes that each variable has the same dimensions. It #####
# calculates the polarization degree P from the Stokes parameters Q and U. The #
## function returns a datafrane with the same spatial dimensions as the input ##
# and 2 parameters "Data" and "Unc"; the first holds the Data, and the second  #
######################## the respective uncertainties.  ########################
################################################################################
P_from_QU <- function(Q, U, unc_Q, unc_U){                                     #
                                                                               #
  dQ <- dim(Q)                                                                 #
  dU <- dim(U)                                                                 #
  duQ <- dim(unc_Q)                                                            #
  duU <- dim(unc_U)                                                            #
                                                                               #
  if(is.null(dQ)){                                                             #
    dQ <- length(Q)                                                            #
  }                                                                            #
  if(is.null(dU)){                                                             #
    dU <- length(U)                                                            #
  }                                                                            #
  if(is.null(duQ)){                                                            #
    duQ <- length(unc_Q)                                                       #
  }                                                                            #
  if(is.null(duU)){                                                            #
    duU <- length(unc_U)                                                       #
  }                                                                            #
                                                                               #
  ldQ <- length(dQ)                                                            #
                                                                               #
  if(!prod(dQ == dU) || !prod(dQ == duQ) || !prod(dQ == duU)){                 #
    print(paste0("All variables must have the same dimensions. Q is ", dQ, ",",#
                 " U is ", dU, ", unc_Q is ", duQ, ", and unc_U is ", duU,     #
                 "."))                                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  P <- NULL                                                                    #
  P$Data <- sqrt(Q^2 + U^2)                                                    #
  P$Unc <- sqrt((Q * unc_Q)^2 + (U * unc_U)^2) / P$Data                        #
                                                                               #
  return(P)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
##### This function assumes that each variable has the same dimensions. It #####
#### calculates the debiased polarization degree estimator, P_bias, from the ###
### Stokes parameters Q and U and their uncertainties. The function returns a ##
### datafrane with the same spatial dimensions as the input and 2 parameters ###
###### "Data" and "Unc"; the first holds the debiased Data, the second the #####
################################ uncertainties. ################################
################################################################################
debiased_P_from_QU <- function(Q, U, unc_Q, unc_U){                            #
                                                                               #
  dQ <- dim(Q)                                                                 #
  dU <- dim(U)                                                                 #
  duQ <- dim(unc_Q)                                                            #
  duU <- dim(unc_U)                                                            #
                                                                               #
  if(is.null(dQ)){                                                             #
    dQ <- length(Q)                                                            #
  }                                                                            #
  if(is.null(dU)){                                                             #
    dU <- length(U)                                                            #
  }                                                                            #
  if(is.null(duQ)){                                                            #
    duQ <- length(unc_Q)                                                       #
  }                                                                            #
  if(is.null(duU)){                                                            #
    duU <- length(unc_U)                                                       #
  }                                                                            #
                                                                               #
  ldQ <- length(dQ)                                                            #
                                                                               #
  if(!prod(dQ == dU) || !prod(dQ == duQ) || !prod(dQ == duU)){                 #
    print(paste0("All variables must have the same dimensions. Q is ", dQ, ",",#
                 " U is ", dU, ", unc_Q is ", duQ, ", and unc_U is ", duU,     #
                 "."))                                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  p2 <- Q^2 + U^2                                                              #
  b2 <- ((Q * unc_U)^2 + (U * unc_Q)^2) / p2                                   #
  p <- sqrt(p2)                                                                #
                                                                               #
  P <- NULL                                                                    #
  P$Data <- p - b2 * (1 - exp(- p2 / b2)) / (2 * p)                            #
  P$Unc <- sqrt((Q * unc_Q)^2 + (U * unc_U)^2) / p                             #
                                                                               #
  return(P)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
##### This function assumes that each variable has the same dimensions. It #####
# calculates the polarization angle X from the Stokes parameters Q and U. The  #
## function returns a dataframe with the same spatial dimensions as the input ##
# and 2 parameters "Data" and "Unc"; the first holds the Data, and the second  #
######################## the respective uncertainties.  ########################
################################################################################
X_from_QU <- function(Q, U, unc_Q, unc_U){                                     #
                                                                               #
  dQ <- dim(Q)                                                                 #
  dU <- dim(U)                                                                 #
  duQ <- dim(unc_Q)                                                            #
  duU <- dim(unc_U)                                                            #
                                                                               #
  if(is.null(dQ)){                                                             #
    dQ <- length(Q)                                                            #
  }                                                                            #
  if(is.null(dU)){                                                             #
    dU <- length(U)                                                            #
  }                                                                            #
  if(is.null(duQ)){                                                            #
    duQ <- length(unc_Q)                                                       #
  }                                                                            #
  if(is.null(duU)){                                                            #
    duU <- length(unc_U)                                                       #
  }                                                                            #
                                                                               #
  ldQ <- length(dQ)                                                            #
                                                                               #
  if(!prod(dQ == dU) || !prod(dQ == duQ) || !prod(dQ == duU)){                 #
    print(paste0("All variables must have the same dimensions. Q is ", dQ, ",",#
                 " U is ", dU, ", unc_Q is ", duQ, ", unc_U is ", duU, "."))   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  P <- sqrt(Q^2 + U^2)                                                         #
                                                                               #
  X <- NULL                                                                    #
  X$Data <- atan2(U, Q) * 90 / pi                                              #
  X$Unc <- 90 / (pi * P^2) * sqrt((U * unc_Q)^2 + (Q * unc_U)^2)               #
                                                                               #
  return(X)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Expects a numerical array "data". Outputs the Z transform of that array. ###
################################################################################
norm_data <- function(data){                                                   #
  return((data - mean(data, na.rm = T)) / sd(data, na.rm = T)^2)               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
###### Expects a numerical array "data". Outputs the log10 of that array. ######
################################################################################
log_data <- function(data){                                                    #
                                                                               #
  data[which(data <= 0, arr.ind = T)] <- NA                                    #
                                                                               #
  return(log10(data))                                                          #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
###### Expects a 1 dimensional array 'x'. Outputs the indexes of 'x' where #####
################### transitions of 1 to 0 and 0 to 1 occur. ####################
################################################################################
find_0to1_and_1to0 <- function(x){                                             #
                                                                               #
  if (!is.vector(x)){                                                          #
    print("ERROR: 'x' must be a 1 dimensional vector. Returning NULL.")        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  lim_bol <- NULL                                                              #
  for(i in 1:length(x)){                                                       #
    if((i == 1 || i == length(x)) && !x[i]){                                   #
      lim_bol <- append(lim_bol, i)                                            #
    }else{                                                                     #
      if(!x[i] && (x[i+1] || x[i-1])){                                         #
        lim_bol <- append(lim_bol, i)                                          #
      }                                                                        #
    }                                                                          #
  }                                                                            #
  return(lim_bol)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects the variables "a" and "b" to have 3 dimensions, the  ##
# third being of length 2, with the first position having the label "Data" and #
### the second position having the label "Unc"; it also expects the variable ###
#### "op" to be one of four possible strings ("+", "-", "*", "/"); it also  ####
## expects "unbin" to either TRUE or FALSE. If these expectations are not met ##
### the function will return NULL. If "bin_l" is NA, the function checks if ####
#### the first two dimensions of "a" and "b" have the same length, if so it ####
### returns the result of "a" "op" "b", else it checks if the length of each ###
# dimension of "a" is either a multiple or a divisor of the same that of "b",  #
#### if it isn't the function will return NULL; if "bin_l" is specified, the ###
### function checks which of "a" and "b" has the larger dimensions; then, if ###
## "unbin" is TRUE the variable with the shorter dimensions will be parsed by ##
### "bin_expand_map()" and if "unbin" is FALSE the variable with the longer  ###
### dimensions will be parsed by "bin_code_map()" and the function will then ###
# return of "op" on the now dimensionally matched variables, else the function #
################################ will return NULL. #############################
################################################################################
dimMatch_and_operate <- function(a, b, op, bin_l = NA, unbin = TRUE,           #
                                 ref_px = NULL){                               #
                                                                               #
  da <- dim(a)                                                                 #
  db <- dim(b)                                                                 #
  lda <- length(da)                                                            #
  ldb <- length(db)                                                            #
  lref <- length(ref_px)                                                       #
                                                                               #
  if(lda != 3 && ldb != 3){                                                    #
    print(paste0("ERROR: this function expects 'a' and 'b' to have 3 dimensio",#
                 "ns, but instead 'a' has ", lda, " and 'b' has ", ldb, ". Re",#
                 "turning NULL."))                                             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(da[3] != 2 && db[3] != 2){                                                #
    print(paste0("ERROR: this function expects the third dimension of 'a' and",#
                 " 'b' to have length 2, but instead they have lengths ",      #
                 da[3], " and ", db[3], ", respectively. Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_px' must be a length 2 or 3 vector. Proceedi", #
                 "ng with usage of 'arr' central pixel as 'ref_px'."))         #
    ref_px <- floor(dims[1:2] / 2)                                             #
  }else{                                                                       #
    ref_px <- floor(ref_px[1:2])                                               #
  }                                                                            #
                                                                               #
  na <- unlist(dimnames(a)[3])                                                 #
  nb <- unlist(dimnames(b)[3])                                                 #
                                                                               #
  if(na[1] != "Data" && na[1] != nb[1] && na[2] != "Unc" && na[2] != nb[2]){   #
    print(paste0("ERROR: this function expects the first position of the thir",#
                 "d dimension of 'a' and 'b' to be labeled 'Data' and the sec",#
                 "ond position to be labeled 'Unc', but instead the labels of",#
                 " the third dimension of 'a' are labeled: ", na, "; and the ",#
                 "those of 'b' are labeled: ", nb, ". Returning NULL."))       #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!is.character(op) || op != "+" && op != "-" && op != "*" && op != "/"){   #
    print(paste0("ERROR: this function expects the variable 'op' to be one of",#
                 " four possible characters ('+','-','*','/'), but instead is",#
                 " ", op, ". Returning NULL."))                                #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!is.logical(unbin)){                                                      #
    print(paste0("ERROR: this function expects the variable 'unbin' to be eit",#
                 "her TRUE or FALSE, but instead is ", unbin, ". Returning NU",#
                 "LL."))                                                       #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!is.na(bin_l) && ((bin_l %% 1 != 0) || bin_l < 1)){                       #
    print(paste0("ERROR: this function expects the variable 'bin_l' to be ",   #
                 "either NA or a integer larger than 1, but instead is ",      #
                 bin_l, ". Returning NULL."))                                  #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  #Determining which of the input arrays will be either compressed/unbinned    #
  #If the dimensions lengths are consistent with expectations                  #
  #And if the dimensions lengths of one array are a multiple of the other array#
  #it also determines bin_l                                                    #
  if(da[1] > db[1] && da[2] > db[2]){                                          #
    a_is_larger <- TRUE                                                        #
    b_is_larger <- FALSE                                                       #
    m1 <- floor(da[1] / db[1])                                                 #
    m2 <- floor(da[2] / db[2])                                                 #
  }else{                                                                       #
    a_is_larger <- FALSE                                                       #
                                                                               #
    if(db[1] > da[1] && db[2] > da[2]){                                        #
      b_is_larger <- TRUE                                                      #
      m1 <- floor(db[1] / da[1])                                               #
      m2 <- floor(db[2] / da[2])                                               #
    }else{                                                                     #
      b_is_larger <- FALSE                                                     #
      m1 <- NULL                                                               #
    }                                                                          #
  }                                                                            #
                                                                               #
  if((da[1] >= db[1] && da[2] < db[2]) || (da[1] <= db[1] && da[2] > db[2])){  #
    print(paste0("ERROR: this function expects the length of the dimensions o",#
                 "f 'a' and 'b' to either match or those of one to be consist",#
                 "ently larger than the others, but instead 'a' has dimension",#
                 "s of length [", da, "] and 'b' of length [", db, "]. Return",#
                 "ing NULL."))                                                 #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.null(m1)){                                                            #
    if(m1 == m2){                                                              #
      if(is.na(bin_l)){                                                        #
        bin_l <- m1                                                            #
      }                                                                        #
    }else{                                                                     #
      print(paste0("ERROR: this function expects the ratio between the length",#
                   " of the dimensions of 'a' and 'b' to be constant, but ins",#
                   "tead the ratio between 'a' and 'b' first dimensions lengt",#
                   "h is ", m1, " and between their second dimensions length ",#
                   "is ", m2, ". Returning NULL."))                            #
      return(NULL)                                                             #
    }                                                                          #
  }                                                                            #
                                                                               #
  #Setting up output array                                                     #
  if(unbin){                                                                   #
    if(a_is_larger){                                                           #
      c <- array(NA, dim = da, dimnames = dimnames(a))                         #
    }else{                                                                     #
      c <- array(NA, dim = db, dimnames = dimnames(b))                         #
    }                                                                          #
  }else{                                                                       #
    if(a_is_larger){                                                           #
      c <- array(NA, dim = db, dimnames = dimnames(b))                         #
    }else{                                                                     #
      rows_f <- c(sort(seq(ref_px[1] + floor(bin_l / 2), 1, -bin_l)),          #
                  seq(ref_px[1] + floor(bin_l / 2) + bin_l, db[1], bin_l))     #
      l_r <- length(rows_f)                                                    #
      rows_i <- c(1, rows_f[1:(l_r - 1)] + 1)                                  #
                                                                               #
      if(rows_f[l_r] > db[1]){                                                 #
        rows_f[l_r] <- db[1]                                                   #
      }                                                                        #
      if(rows_f[l_r] <  db[1]){                                                #
        rows_f <- append(rows_f, db[1])                                        #
        rows_i <- append(rows_i, rows_f[l_r] + 1)                              #
      }                                                                        #
      l_r <- length(rows_i)                                                    #
                                                                               #
      cols_f <- c(sort(seq(ref_px[2] + floor(bin_l / 2), 1, -bin_l)),          #
                  seq(ref_px[2] + floor(bin_l / 2) + bin_l, db[2], bin_l))     #
      l_c <- length(cols_f)                                                    #
      cols_i <- c(1, cols_f[1:(l_c - 1)] + 1)                                  #
                                                                               #
      if(cols_f[l_c] > db[2]){                                                 #
        cols_f[l_c] <- db[2]                                                   #
      }                                                                        #
      if(cols_f[l_c] < db[2]){                                                 #
        cols_f <- append(cols_f, db[2])                                        #
        cols_i <- append(cols_i, cols_f[l_c] + 1)                              #
      }                                                                        #
      l_c <- length(cols_i)                                                    #
                                                                               #
      c <- array(NA, dim = c(l_r, l_c, 2), dimnames = dimnames(a))             #
    }                                                                          #
  }                                                                            #
                                                                               #
  #Calculating output                                                          #
  if((bin_l %% 1 == 0) && bin_l >= 1){                                         #
    if(bin_l > 1){                                                             #
      if(unbin){                                                               #
        if(a_is_larger){                                                       #
          b <- bin_expand_map(b, bin_l, da, ref_pix = ref_px)                  #
        }else{                                                                 #
          if(b_is_larger){                                                     #
            a <- bin_expand_map(a, bin_l, db, ref_pix = ref_px)                #
          }else{                                                               #
            b <- bin_expand_map(b, bin_l, da, ref_pix = ref_px)                #
            a <- bin_expand_map(a, bin_l, da, ref_pix = ref_px)                #
          }                                                                    #
        }                                                                      #
      }else{                                                                   #
        if(a_is_larger){                                                       #
          a <- bin_code_map(a, bin_l, ref_px)                                  #
        }else{                                                                 #
          b <- bin_code_map(b, bin_l, ref_px)                                  #
          if(!b_is_larger){                                                    #
            a <- bin_code_map(a, bin_l, ref_px)                                #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
    c[,,1] <- oper_data(a, b, op)                                              #
    c[,,2] <- oper_unc(a, b, op)                                               #
  }else{                                                                       #
    print(paste0("ERROR: this function expects 'bin_l' to be an integer large",#
                 "r than 1, but instead is ", bin_l, ". Returning NULL."))     #
    return(NULL)                                                               #
  }                                                                            #
  return(c)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects the variables "a" and "b" to have 3 dimensions, and  ##
# dimensions of matching lengths, the third being of length 2, with the first  #
## position having the label "Data" and the second position having the label  ##
# "Unc"; it also expects the variable "op" to be one of four possible strings  #
## ("+", "-", "*", "/"). If these expectations are not met the function will  ##
##### return NULL, if they are the function will then return "a" "op" "b". #####
################################################################################
oper_data <- function(a, b, op){                                               #
                                                                               #
  da <- dim(a)                                                                 #
  db <- dim(b)                                                                 #
  lda <- length(da)                                                            #
  ldb <- length(db)                                                            #
                                                                               #
  if(lda != 3 && ldb != 3){                                                    #
    print(paste0("ERROR: this function expects 'a' and 'b' to have 3 dimensio",#
                 "ns, but instead 'a' has ", lda, " and 'b' has ", ldb, ". Re",#
                 "turning NULL."))                                             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(da[3] != 2 && db[3] != 2){                                                #
    print(paste0("ERROR: this function expects the third dimension of 'a' and",#
                 " 'b' to have length 2, but instead they have lengths ",      #
                 da[3], " and ", db[3], ", respectively. Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(da[1] != db[1] || da[2] != db[2]){                                        #
    print(paste0("ERROR: this function expects the all dimensions of 'a' and ",#
                 "'b' to have same length, but instead they have lengths ", da,#
                 " and ", db, ", respectively. Returning NULL."))              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  na <- unlist(dimnames(a)[3])                                                 #
  nb <- unlist(dimnames(b)[3])                                                 #
                                                                               #
  if(na[1] != "Data" && na[1] != nb[1] && na[2] != "Unc" && na[2] != nb[2]){   #
    print(paste0("ERROR: this function expects the first position of the thir",#
                 "d dimension of 'a' and 'b' to be labeled 'Data' and the sec",#
                 "ond position to be labeled 'Unc', but instead the labels of",#
                 " the third dimension of 'a' are labeled: ", na, "; and the ",#
                 "those of 'b' are labeled: ", nb, ". Returning NULL."))       #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!is.character(op) || op != "+" && op != "-" && op != "*" && op != "/"){   #
    print(paste0("ERROR: this function expects the variable 'op' to be one of",#
                 " four possible characters ('+','-','*','/'), but instead is",#
                 " ", op, ". Returning NULL."))                                #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(op == '+'){                                                               #
    op <- 'sum'                                                                #
  }                                                                            #
  if(op == '-'){                                                               #
    op <- 'dif'                                                                #
  }                                                                            #
  if(op == '*'){                                                               #
    op <- 'mult'                                                               #
  }                                                                            #
  if(op == '/'){                                                               #
    op <- 'div'                                                                #
  }                                                                            #
                                                                               #
  switch(op,                                                                   #
         sum = a[,,1] + b[,,1], dif = a[,,1] - b[,,1],                         #
         mult = a[,,1] * b[,,1], div = a[,,1] / b[,,1])                        #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## This function expects the variables "a" and "b" to have 3 dimensions, and  ##
# dimensions of matching lengths, the third being of length 2, with the first  #
## position having the label "Data" and the second position having the label  ##
# "Unc"; it also expects the variable "op" to be one of four possible strings  #
## ("+", "-", "*", "/"). If these expectations are not met the function will  ##
## return NULL, if they are the function will then return the uncertainty of  ##
################################# "a" "op" "b". ################################
################################################################################
oper_unc <- function(a, b, op){                                                #
                                                                               #
  da <- dim(a)                                                                 #
  db <- dim(b)                                                                 #
  lda <- length(da)                                                            #
  ldb <- length(db)                                                            #
                                                                               #
  if(lda != 3 && ldb != 3){                                                    #
    print(paste0("ERROR: this function expects 'a' and 'b' to have 3 dimensio",#
                 "ns, but instead 'a' has ", lda, " and 'b' has ", ldb, ". Re",#
                 "turning NULL."))                                             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(da[3] != 2 && db[3] != 2){                                                #
    print(paste0("ERROR: this function expects the third dimension of 'a' and",#
                 " 'b' to have length 2, but instead they have lengths ",      #
                 da[3], " and ", db[3], ", respectively. Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(da[1] != db[1] || da[2] != db[2]){                                        #
    print(paste0("ERROR: this function expects all dimensions of 'a' and 'b' ",#
                 "to have same length, but instead they have lengths ", da,    #
                 " and ", db, ", respectively. Returning NULL."))              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  na <- unlist(dimnames(a)[3])                                                 #
  nb <- unlist(dimnames(b)[3])                                                 #
                                                                               #
  if(na[1] != "Data" && na[1] != nb[1] && na[2] != "Unc" && na[2] != nb[2]){   #
    print(paste0("ERROR: this function expects the first position of the thir",#
                 "d dimension of 'a' and 'b' to be labeled 'Data' and the sec",#
                 "ond position to be labeled 'Unc', but instead the labels of",#
                 " the third dimension of 'a' are labeled: ", na, "; and the ",#
                 "those of 'b' are labeled: ", nb, ". Returning NULL."))       #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!is.character(op) || op != "+" && op != "-" && op != "*" && op != "/"){   #
    print(paste0("ERROR: this function expects the variable 'op' to be one of",#
                 " four possible characters ('+','-','*','/'), but instead is",#
                 " ", op, ". Returning NULL."))                                #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(op == '+'){                                                               #
    op <- 'sum'                                                                #
  }                                                                            #
  if(op == '-'){                                                               #
    op <- 'dif'                                                                #
  }                                                                            #
  if(op == '*'){                                                               #
    op <- 'mult'                                                               #
  }                                                                            #
  if(op == '/'){                                                               #
    op <- 'div'                                                                #
  }                                                                            #
                                                                               #
  switch(op,                                                                   #
         sum = unc_add(a[,,2], b[,,2]),                                        #
         dif = unc_add(a[,,2], b[,,2]),                                        #
         mult = unc_mult(a[,,1], b[,,1], a[,,2], b[,,2]),                      #
         div = unc_div(a[,,1], b[,,1], a[,,2], b[,,2]))                        #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
###### Expects a numerical array "data", a replacement number "rep_val", a #####
## limiting number "num_to_rep" and relation symbol "rel". The function then ###
# The function then determines the positions within "data" that hold +-Inf and #
######################## replaces those with "rep_val". ########################
################################################################################
replace_qtts <- function(data, num_to_rep, rel = "=", rep_val = NA){           #
                                                                               #
  if(!is.finite(num_to_rep) || length(num_to_rep) != 1){                       #
    print(paste0("ERROR: 'num_to_rep' must be a length 1 finite variable. Ret",#
                 "urning NULL."))                                              #
    return(NULL)                                                               #
  }                                                                            #
  if(length(rep_val) != 1){                                                    #
    print(paste0("ERROR: 'rep_val' must be a length 1 variable. Returning NUL",#
                 "L."))                                                        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  possible_rels <- c("=", "<", ">", "<=", ">=")                                #
                                                                               #
  if(switch(rel, "=" = res <- 'eq', "<" = res <- 'sml', ">" = res <- 'lrg',    #
            "<=" = res <- 'seq', ">=" = res <- 'leq', "Err") == "Err"){        #
    print(paste0("ERROR: Invalid 'rel', it should be one of {'=','<','>','<='",#
                 ",'>='}. Returning NULL."))                                   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  switch(res,                                                                  #
         eq = c(inds <- which(data == num_to_rep, arr.ind = T),                #
                data[inds] <- rep_val),                                        #
         sml = c(inds <- which(data < num_to_rep, arr.ind = T),                #
                 data[inds] <- rep_val),                                       #
         lrg = c(inds <- which(data > num_to_rep, arr.ind = T),                #
                 data[inds] <- rep_val),                                       #
         seq = c(inds <- which(data <= num_to_rep, arr.ind = T),               #
                 data[inds] <- rep_val),                                       #
         leq = c(inds <- which(data >= num_to_rep, arr.ind = T),               #
                 data[inds] <- rep_val),                                       #
  )                                                                            #
                                                                               #
  return(data)                                                                 #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
##### Expects a numerical array "data" and a replacement number "rep_val". #####
# The function then determines the positions within "data" that hold +-Inf and #
######################## replaces those with "rep_val". ########################
################################################################################
remove_infs <- function(data, rep_val = NA){                                   #
                                                                               #
  if(length(rep_val) != 1){                                                    #
    print(paste0("ERROR: 'rep_val' must be a length 1 variable. Returning NUL",#
                 "L."))                                                        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  inds <- which(is.infinite(data), arr.ind = T)                                #
  data[inds] <- rep_val                                                        #
                                                                               #
  return(data)                                                                 #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Expects two numerical data arrays, with the same dimensions, "data" and ###
## "target", a number to be replaced, "num_to_rep", and a replacement number, ##
##### "rep_val". The function then seeks in "data" the positions with value ####
## "num_to_rep", save the indexes of those positions and will place "rep_val" ##
##################### in those positions within "target". ######################
################################################################################
bij_match_val <- function(data, target, num_to_rep, rep_val = NA){             #
                                                                               #
  dimd <- dim(data)                                                            #
  dimt <- dim(target)                                                          #
                                                                               #
  if(!prod(dimd == dimt)){                                                     #
    print(paste0("ERROR: 'data' and 'target' dimensions must be equal. Return",#
    "ing NULL."))                                                              #
  }                                                                            #
                                                                               #
  if(length(rep_val) != 1){                                                    #
    print(paste0("ERROR: 'rep_val' must be a length 1 variable. Returning NUL",#
                 "L."))                                                        #
    return(NULL)                                                               #
  }                                                                            #
  if(length(num_to_rep) != 1){                                                 #
    print(paste0("ERROR: 'num_to_rep' must be a length 1 variable. Returning ",#
                 "NULL."))                                                     #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(is.na(num_to_rep)){                                                       #
    inds <- which(is.na(data), arr.ind = T)                                    #
  }else{                                                                       #
    if(is.infinite(num_to_rep)){                                               #
      inds <- which(is.infinite(data), arr.ind = T)                            #
    }else{                                                                     #
      inds <- which(data == num_to_rep, arr.ind = T)                           #
    }                                                                          #
  }                                                                            #
  target[inds] <- rep_val                                                      #
                                                                               #
  return(target)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## Expects a numerical array 'data' with 2 or 3 dimensions, a numerical vector #
# 'cpix' holding the 2d coordinates of the pixel to the center of the cropping #
### and another numerical vector 'delta' with the same length as 'cpix' that ###
# will define the width of the box to be croped within the first 2 dimensions ##
#### of 'data'. Outputs and new array with the same number of dimensions as ####
#### 'data' but holding only the data within the previously defined region. ####
################################################################################
crop_array <- function(data, cpix, delta){                                     #
                                                                               #
  if(!(length(cpix) == 2 && length(delta) == 2 &&                              #
       is.vector(cpix) && is.vector(delta))){                                  #
    print(paste0("ERROR: 'cpix' and 'delta' must be 1-dimensional vectors of ",#
                 "length 2. Returning NULL"))                                  #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!(all(delta %% 1 == 0) && all(delta > 0) &&                               #
       all(cpix %% 1 == 0) && all(cpix > 0))){                                 #
    print(paste0("ERROR: the components of 'cpix' and 'delta' must be positiv",#
                 "e integers. Returning NULL."))                               #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dimd <- dim(data)                                                            #
                                                                               #
  if(is.null(dimd) || (length(dimd) != 2 && length(dimd) != 3)){               #
    print(paste0("ERROR: 'data' must be an array with either 2 or 3 dimension",#
                 "s. Returning NULL."))                                        #
  }                                                                            #
                                                                               #
  if(cpix[1] < 1 || cpix[1] > dimd[1] || cpix[2] < 1 || cpix[2] > dimd[2]){    #
    print("ERROR: 'cpix' is outside the boundaries of 'data'. Returning NULL.")#
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  ri <- cpix[1] - delta[1]                                                     #
  rf <- cpix[1] + delta[1]                                                     #
  ci <- cpix[2] - delta[2]                                                     #
  cf <- cpix[2] + delta[2]                                                     #
                                                                               #
  if(ri < 1 || rf > dimd[1] || ci < 1 || cf > dimd[2]){                        #
    print(paste0("WARNING: the region to be cropped surpasses the limits of t",#
                 "he input map. The output will coerced to those limits."))    #
                                                                               #
    if(ri < 1){                                                                #
      ri <- 1                                                                  #
    }                                                                          #
    if(ci < 1){                                                                #
      ci <- 1                                                                  #
    }                                                                          #
    if(rf < 1){                                                                #
      rf <- dimd[1]                                                            #
    }                                                                          #
    if(cf < 1){                                                                #
      cf <- dimd[2]                                                            #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(length(dimd) == 2){                                                       #
    cropd_data <- data[ri:rf, ci:cf]                                           #
  }else{                                                                       #
    cropd_data <- data[ri:rf, ci:cf,]                                          #
  }                                                                            #
                                                                               #
  return(cropd_data)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Expects three numerical 2d data arrays, with the same dimensions, "mag" ####
### "ang" and "bkg", which are to hold respectively the magnitudes and angles ##
# (in degrees) of the vectors to be plotted and a background image for them to #
### sit on; it also expects a binning size "bin", in case the vector spatial ###
## density is too large; lastly it expects an output folder, name and title, ###
## "outpath" "outname" and "mtitle". The function then produces an arrow plot ##
################ sitting on an image and saves it to a PDF file. ###############
######### If 'relative' is TRUE, 'mag' will be converted to percentage. ########
################################################################################
create_arrow_plot_pdf <- function(mag, u_mag = sqrt(mag), ang,                 #
                                  u_ang = sqrt(ang), bkg = mag,                #
                                  binsize = 1, outpath = "~/",                 #
                                  outname = "arrow_plot", relative = TRUE,     #
                                  mtitle = "Vector map", a_col = "green",      #
                                  points = NULL, mask = 1, x_ticks = NULL,     #
                                  y_ticks = NULL, x_labs = NULL, y_labs = NULL,#
                                  ref_pix = NULL, stat_flag = TRUE){           #
                                                                               #
  chars <- c(3, 4, 1, 8)                                                       #
  colors <- c("cyan", "yellow", "pink", "orange")                              #
                                                                               #
  std_dims <- c(660, 660)                                                      #
                                                                               #
  dimd <- dim(mag)                                                             #
  dimud <- dim(u_mag)                                                          #
  dima <- dim(ang)                                                             #
  dimad <- dim(u_ang)                                                          #
  lref <- length(ref_pix)                                                      #
  dimb <- dim(bkg)                                                             #
  delt <- dimd[1] * .8 / dimd[2]                                               #
                                                                               #
  if(is.null(points)){                                                         #
    lpoints <- 0                                                               #
  }else{                                                                       #
    dimpoints <- dim(points)                                                   #
    if(length(dimpoints) != 2 || dimpoints[2] != 2 || dimpoints[1] > 4){       #
      print(paste0("WARNING: 'points' is expected to be a two-column matrix w",#
                   "ith a maximum 4 rows. Ignoring points."))                  #
      lpoints <- 0                                                             #
    }else{                                                                     #
      lpoints <- dimpoints[1]                                                  #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(is.null(dimd) || is.null(dimud) || is.null(dima) ||  is.null(dimad) ||    #
     is.null(dimb) || length(dimd) != 2 || length(dimud) != 2 ||               #
     length(dima) != 2 || length(dimad) != 2 || length(dimb) != 2 ||           #
     !prod(dimd == dimud) || !prod(dimd == dima) || !prod(dimd == dimad) ||    #
     !prod(dimd == dimb)){                                                     #
    print(paste0("ERROR: 'mag', 'u_mag', 'ang', 'u_ang' and 'bkg' must all be",#
                 " 2d arrays and have the same dimensions. Returning NULL."))  #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(length(mask) != 1 && prod(dim(mask) != dimd)){                            #
    print(paste0("WARNING: 'mask' must either be '1' or a 2d array with the ", #
                 "same dimensions as 'mag'. Proceeding with 'mask' = 1."))     #
  }                                                                            #
                                                                               #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_pix' must be a length 2 or 3 vector. Proceedi",#
                 "ng with usage of 'arr' central pixel as 'ref_pix'."))        #
    ref_pix <- floor(dimd[1:2] / 2)                                            #
  }else{                                                                       #
    ref_pix <- floor(ref_pix[1:2])                                             #
  }                                                                            #
                                                                               #
  if(length(outpath) == 0){                                                    #
    print("WARNING: Output path is empty, coerced to '~/'.")                   #
    outpath <- "~/"                                                            #
  }                                                                            #
  if(!dir.exists(outpath)){                                                    #
    mkdir(outpath)                                                             #
  }                                                                            #
                                                                               #
  if(length(outname) == 0){                                                    #
    outname <- "arrow_plot"                                                    #
    print("WARNING: Output name is empty, coercing it to 'arrow_plot'.")       #
  }                                                                            #
                                                                               #
  i <- 1                                                                       #
  while(file.exists(paste0(outpath, outname, ".pdf"))){                        #
    sep_name <- strsplit(outname, "[(]")[[1]]                                  #
    outname <- paste0(sep_name[1], "(", i, ")")                                #
                                                                               #
    print(paste0("WARNING: There is already a file with that name in the spec",#
                 "ified folder, coercing it to ", outname, " ."))              #
    i <- i + 1                                                                 #
  }                                                                            #
                                                                               #
  out_pdf <- paste0(outpath, outname, ".pdf")                                  #
                                                                               #
  delt <- .8                                                                   #
  rat <- dimd[2] / dimd[1]                                                     #
  rx <- std_dims[1] / dimd[1]                                                  #
  ry <- std_dims[2] / dimd[2]                                                  #
                                                                               #
  rows <- rep(1:dimd[1], dimd[2])                                              #
  cols <- rep(1:dimd[2], each = dimd[1])                                       #
  xs <- (rows - .5) / dimd[1]                                                  #
  ys <- (cols - .5) / dimd[2]                                                  #
                                                                               #
  temp_ang <- as.vector((ang * mask))                                          #
  temp_u_ang <- as.vector((u_ang * mask))                                      #
  temp_mag <- as.vector((mag * mask))                                          #
  temp_u_mag <- as.vector((u_mag * mask))                                      #
                                                                               #
  NAsMag <- which(is.na(temp_mag), arr.ind = T)                                #
  NAsAng <- which(is.na(temp_ang), arr.ind = T)                                #
  nonNAsMag <- which(!is.na(temp_mag), arr.ind = T)                            #
                                                                               #
  if(relative){                                                                #
    temp_mag <- temp_mag * 100                                                 #
  }                                                                            #
                                                                               #
  min_mag <- round(min(temp_mag, na.rm = T), 2)                                #
  med_mag <- round(median(temp_mag, na.rm = T), 2)                             #
  unc_med_mag <- round(unc_median(temp_mag, temp_u_mag, 0), 2)                 #
  max_mag <- round(max(temp_mag, na.rm = T), 2)                                #
                                                                               #
  if(relative){                                                                #
    mag_stat_str <- paste0("min mag = ", min_mag, "%, median mag = ", med_mag, #
                            "% +- ", unc_med_mag, "%, max mag = ", max_mag,"%")#
  }else{                                                                       #
    mag_stat_str <- paste0("min mag = ", min_mag, ", median mag = ", med_mag,  #
                           " +- ", unc_med_mag, ", max mag = ", max_mag)       #
  }                                                                            #
                                                                               #
  min_ang <- round(min(temp_ang, na.rm = T), 2)                                #
  med_ang <- round(median(temp_ang, na.rm = T), 2)                             #
  unc_med_ang <- round(unc_median(temp_ang, temp_u_ang, 0), 2)                 #
  max_ang <- round(max(temp_ang, na.rm = T), 2)                                #
                                                                               #
  ang_stat_str <- paste0("min ang = ", min_ang, "º, median ang = ", med_ang,   #
                         "º +- ", unc_med_ang, "º, max ang = ", max_ang, "º")  #
                                                                               #
  na_inds <- union(which(is.na(temp_ang), arr.ind = TRUE),                     #
                   which(is.na(temp_mag), arr.ind = TRUE))                     #
                                                                               #
  if(length(na_inds) != 0){                                                    #
    xs <- xs[-na_inds]                                                         #
    ys <- ys[-na_inds]                                                         #
    temp_ang <- temp_ang[-na_inds]                                             #
    temp_u_ang <- temp_u_ang[-na_inds]                                         #
    temp_mag <- temp_mag[-na_inds]                                             #
    temp_u_mag <- temp_u_mag[-na_inds]                                         #
  }                                                                            #
                                                                               #
  # The angle is given in relation to the optical axis of the WP               #
  # So for proper arrow representation a rotation is needed                    #
  temp_ang <- temp_ang + 90                                                    #
                                                                               #
  bkg[which(is.na(bkg) | is.infinite(bkg), arr.ind = T)] <- 0                  #
                                                                               #
  legticks <- round(seq(min(bkg[which(!is.infinite(bkg))], na.rm = T),         #
                        max(bkg[which(!is.infinite(bkg))], na.rm = T),         #
                        length.out = 5))                                       #
  leglabs <- paste0(legticks)                                                  #
                                                                               #
  mainleg_cex <- 4 * ry                                                        #
  labax_cex <- 2 * ry                                                          #
  tit_cex <- 3 * ry                                                            #
  ptsStat_cex <- 2 * ry                                                        #
  ref_xlen <- 0.025 * rat * rx                                                 #
  ref_ylen <- 0.025 * rat * ry                                                 #
  tip_xlen <- 0.125 * rat * rx                                                 #
  tip_ylen <- 0.125 * rat * ry                                                 #
  ref_xlwd <- 4 * ry                                                           #
  ref_ylwd <- 4 * rx                                                           #
  arr_lwd <- 3 * mean(c(rx, ry))                                               #
  ax_lwd <- 2 * mean(c(rx, ry))                                                #
  scl_b <- 0.01 * mean(c(rx, ry))                                              #
                                                                               #
  tryCatch({                                                                   #
    pdf(out_pdf, width = round(30 / rat), height = 30, bg = "black")           #
    par(cex.main = mainleg_cex)                                                #
    par(cex.axis = labax_cex)                                                  #
    par(col.axis = "white")                                                    #
    par(cex.lab = labax_cex)                                                   #
                                                                               #
    #Main Plot                                                                 #
    image.plot(bkg, xlab =' ', ylab =' ', xaxt ="n", yaxt ="n",                #
               col = gray.colors(96, 0, 1, 2.5),                               #
               bigplot = c(.13, .13 + delt, .13, .93),                         #
               smallplot = c(.14 + delt, .15 + delt, .13, .93),                #
               legend.lab = "Counts", legend.cex = mainleg_cex, legend.line=15,#
               axis.args = list(at = legticks, labels = leglabs, col ="white"))#
    #Arrow Field                                                               #
    # Getting uv projections, scaling to max_mag and dividing by two           #
    # because the projections will be centered on a given coordinate           #
    u <- 0.5 * temp_mag * cos(temp_ang * pi / 180) / max_mag                   #
    v <- 0.5 * temp_mag * sin(temp_ang * pi / 180) / max_mag                   #
                                                                               #
    # Scaling projections to plot grid                                         #
    us <- u * binsize / dimd[1]                                                #
    vs <- v * binsize / dimd[2]                                                #
                                                                               #
    scl_mag <- 0.5 * binsize * floor(max_mag) / max_mag / dimd[1]              #
                                                                               #
    arrows(xs - us, ys - vs, xs + us, ys + vs, length = 0, lwd = arr_lwd,      #
           col = a_col)                                                        #
                                                                               #
    #Points                                                                    #
    if(lpoints > 0){                                                           #
      for(p in 1:lpoints){                                                     #
        points(points[p, 1], points[p, 2], pch = chars[p], col = colors[p],    #
               cex = ptsStat_cex)                                              #
      }                                                                        #
    }                                                                          #
    # Main Title                                                               #
    mtext(mtitle, cex = tit_cex, col = "white", line = 5)                      #
    #N/E referencial                                                           #
    arrows(x0 = .15, y0 = 0.05, x1 = .15, y1 = 0.05 + ref_ylen,                #
           length = tip_ylen, col = "red", lwd = ref_ylwd)                     #
    arrows(x0 = .15, y0 = 0.05, x1 = .15 - ref_xlen, y1 = 0.05,                #
           length = tip_xlen, col = "red", lwd = ref_xlwd)                     #
    text(x = .15, y = 0.05 + ref_ylen * 1.05, labels = "N", pos = 3,           #
         col = "red", cex = labax_cex)                                         #
    text(x = .15 - ref_xlen * 1.05, y = 0.05, labels = "E", pos = 2,           #
         col = "red", cex = labax_cex)                                         #
                                                                               #
    line_offs <- 0                                                             #
    if(stat_flag){                                                             #
      # Scale legend                                                           #
      arrows(x0 = .1 + delt - scl_mag, y0 = 0.05,                              #
             x1 = .1 + delt + scl_mag, y1 = 0.05, length = 0, col = "red",     #
             lwd = arr_lwd)                                                    #
      text(x = .1 + delt + scl_mag + scl_b, y = 0.05,                          #
           labels = paste0(floor(max_mag), "%"), pos = 4, col = "red",         #
           cex = labax_cex)                                                    #
      text(x = .1 + delt - scl_mag - scl_b, y = 0.05, labels = "scale",        #
           pos = 2, col = "red", cex = labax_cex)                              #
                                                                               #
      #Statistical Info                                                        #
      mtext(mag_stat_str, side = 1, line = 3, cex = ptsStat_cex, col = "white")#
      mtext(ang_stat_str, side = 1, line = 5, cex = ptsStat_cex, col = "white")#
      line_offs <- 7                                                           #
    }                                                                          #
    #X-axis                                                                    #
    axis(1, at = x_ticks, labels = paste0(x_labs, "°"), line = 1 + line_offs,  #
         lwd = ax_lwd, padj = .5, col = "white")                               #
    mtext("Ra (deg)", side = 1, cex = tit_cex, line = 7 + line_offs,           #
          col = "white")                                                       #
    #Y-axis                                                                    #
    axis(2, at = y_ticks, labels = paste0(y_labs, "°"), lwd = ax_lwd, line = 1,#
         col = "white")                                                        #
    mtext("Dec (deg)", side = 2, cex = tit_cex, line = 6, col = "white")       #
    dev.off()                                                                  #
    }, error = function(e) {                                                   #
      print("ERROR: Could not create an arrow plot for the given input.")      #
      print(paste("Error message:", e$message))                                #
      return(e)                                                                #
    }                                                                          #
  )                                                                            #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Expects three numerical 2d data arrays, with the same dimensions, "mag" ####
### "ang" and "bkg", which are to hold respectively the magnitudes and angles ##
# (in degrees) of the vectors to be plotted and a background image for them to #
### sit on; it also expects a binning size "bin", in case the vector spatial ###
## density is too large; lastly it expects an output folder, name and title, ###
## "outpath" "outname" and "mtitle". The function then produces an arrow plot ##
################ sitting on an image and saves it to a PDF file. ###############
######### If 'relative' is TRUE, 'mag' will be converted to percentage. ########
################################################################################
create_arrow_iso_plot_pdf <- function(mag, u_mag = sqrt(mag), ang,             #
                                      u_ang = sqrt(ang), bkg = mag,            #
                                      binsize = 1, outpath = "~/",             #
                                      outname = "arrow_iso_plot", relative = T,#
                                      mtitle = "Vector map", a_col = "black",  #
                                      points = NULL, mask = 1, x_ticks = NULL, #
                                      y_ticks = NULL, x_labs = NULL,           #
                                      y_labs = NULL, ref_pix = NULL,           #
                                      stat_flag = TRUE){                       #
                                                                               #
  chars <- c(3, 4, 1, 8)                                                       #
  colors <- c("blue", "red", "purple", "orange")                               #
                                                                               #
  std_dims <- c(660, 660)                                                      #
                                                                               #
  dimd <- dim(mag)                                                             #
  dimud <- dim(u_mag)                                                          #
  dima <- dim(ang)                                                             #
  dimad <- dim(u_ang)                                                          #
  lref <- length(ref_pix)                                                      #
  dimb <- dim(bkg)                                                             #
                                                                               #
  if(is.null(points)){                                                         #
    lpoints <- 0                                                               #
  }else{                                                                       #
    dimpoints <- dim(points)                                                   #
    if(length(dimpoints) != 2 || dimpoints[2] != 2 || dimpoints[1] > 4){       #
      print(paste0("WARNING: 'points' is expected to be a two-column matrix w",#
                   "ith a maximum 4 rows. Ignoring points."))                  #
      lpoints <- 0                                                             #
    }else{                                                                     #
      lpoints <- dimpoints[1]                                                  #
    }                                                                          #
  }                                                                            #
  if(is.null(dimd) || is.null(dimud) || is.null(dima) ||  is.null(dimad) ||    #
     is.null(dimb) || length(dimd) != 2 || length(dimud) != 2 ||               #
     length(dima) != 2 || length(dimad) != 2 || length(dimb) != 2 ||           #
     !prod(dimd == dimud) || !prod(dimd == dima) || !prod(dimd == dimad) ||    #
     !prod(dimd == dimb)){                                                     #
    print(paste0("ERROR: 'mag', 'u_mag', 'ang', 'u_ang' and 'bkg' must all be",#
                 " 2d arrays and have the same dimensions. Returning NULL."))  #
    return(NULL)                                                               #
  }                                                                            #
  if(length(mask) != 1 && prod(dim(mask) != dimd)){                            #
    print(paste0("WARNING: 'mask' must either be '1' or a 2d array with the ", #
                 "same dimensions as 'mag'. Proceeding with 'mask' = 1."))     #
  }                                                                            #
  if(lref != 2 && lref != 3){                                                  #
    print(paste0("WARNING: 'ref_pix' must be a length 2 or 3 vector. Proceedi",#
                 "ng with usage of 'arr' central pixel as 'ref_pix'."))        #
    ref_pix <- floor(dimd[1:2] / 2)                                            #
  }else{                                                                       #
    ref_pix <- floor(ref_pix[1:2])                                             #
  }                                                                            #
  if(length(outpath) == 0){                                                    #
    print("WARNING: Output path is empty, coerced to '~/'.")                   #
    outpath <- "~/"                                                            #
  }                                                                            #
  if(!dir.exists(outpath)){                                                    #
    mkdir(outpath)                                                             #
  }                                                                            #
  if(length(outname) == 0){                                                    #
    outname <- "arrow_plot"                                                    #
    print("WARNING: Output name is empty, coercing it to 'arrow_plot'.")       #
  }                                                                            #
                                                                               #
  i <- 1                                                                       #
  while(file.exists(paste0(outpath, outname, ".pdf"))){                        #
    sep_name <- strsplit(outname, "[(]")[[1]]                                  #
    outname <- paste0(sep_name[1], "(", i, ")")                                #
                                                                               #
    print(paste0("WARNING: There is already a file with that name in the spec",#
                 "ified folder, coercing it to ", outname, " ."))              #
    i <- i + 1                                                                 #
  }                                                                            #
                                                                               #
  out_pdf <- paste0(outpath, outname, ".pdf")                                  #
                                                                               #
  delt <- .8                                                                   #
  rat <- dimd[2] / dimd[1]                                                     #
  rx <- std_dims[1] / dimd[1]                                                  #
  ry <- std_dims[2] / dimd[2]                                                  #
                                                                               #
  rows <- rep(1:dimd[1], dimd[2])                                              #
  cols <- rep(1:dimd[2], each = dimd[1])                                       #
  xs <- (rows - .5) / dimd[1]                                                  #
  ys <- (cols - .5) / dimd[2]                                                  #
                                                                               #
  temp_ang <- as.vector((ang * mask))                                          #
  temp_u_ang <- as.vector((u_ang * mask))                                      #
  temp_mag <- as.vector((mag * mask))                                          #
  temp_u_mag <- as.vector((u_mag * mask))                                      #
                                                                               #
  NAsMag <- which(is.na(temp_mag), arr.ind = T)                                #
  NAsAng <- which(is.na(temp_ang), arr.ind = T)                                #
  nonNAsMag <- which(!is.na(temp_mag), arr.ind = T)                            #
                                                                               #
  if(relative){                                                                #
    temp_mag <- temp_mag * 100                                                 #
  }                                                                            #
                                                                               #
  min_mag <- round(min(temp_mag, na.rm = T), 2)                                #
  med_mag <- round(median(temp_mag, na.rm = T), 2)                             #
  unc_med_mag <- round(unc_median(temp_mag, temp_u_mag, 0), 2)                 #
  max_mag <- round(max(temp_mag, na.rm = T), 2)                                #
                                                                               #
  if(relative){                                                                #
    mag_stat_str <- paste0("min mag = ", min_mag, "%, median mag = ", med_mag, #
                           "% +- ", unc_med_mag, "%, max mag = ", max_mag,"%") #
  }else{                                                                       #
    mag_stat_str <- paste0("min mag = ", min_mag, ", median mag = ", med_mag,  #
                           " +- ", unc_med_mag, ", max mag = ", max_mag)       #
  }                                                                            #
                                                                               #
  min_ang <- round(min(temp_ang, na.rm = T), 2)                                #
  med_ang <- round(median(temp_ang, na.rm = T), 2)                             #
  unc_med_ang <- round(unc_median(temp_ang, temp_u_ang, 0), 2)                 #
  max_ang <- round(max(temp_ang, na.rm = T), 2)                                #
                                                                               #
  ang_stat_str <- paste0("min ang = ", min_ang, "º, median ang = ", med_ang,   #
                         "º +- ", unc_med_ang, "º, max ang = ", max_ang, "º")  #
                                                                               #
  na_inds <- union(which(is.na(temp_ang), arr.ind = TRUE),                     #
                   which(is.na(temp_mag), arr.ind = TRUE))                     #
                                                                               #
  if(length(na_inds) != 0){                                                    #
    xs <- xs[-na_inds]                                                         #
    ys <- ys[-na_inds]                                                         #
    temp_ang <- temp_ang[-na_inds]                                             #
    temp_u_ang <- temp_u_ang[-na_inds]                                         #
    temp_mag <- temp_mag[-na_inds]                                             #
    temp_u_mag <- temp_u_mag[-na_inds]                                         #
  }                                                                            #
                                                                               #
  # The angle is given in relation to the optical axis of the WP               #
  # So for proper arrow representation a rotation is needed                    #
  temp_ang <- temp_ang + 90                                                    #
                                                                               #
  bkg[which(is.na(bkg) | is.infinite(bkg), arr.ind = T)] <- 0                  #
                                                                               #
  legticks <- round(seq(min(bkg[which(!is.infinite(bkg))], na.rm = T),         #
                        max(bkg[which(!is.infinite(bkg))], na.rm = T),         #
                        length.out = 5))                                       #
  leglabs <- paste0(legticks)                                                  #
                                                                               #
  mainleg_cex <- 4 * ry                                                        #
  labax_cex <- 2 * ry                                                          #
  tit_cex <- 3 * ry                                                            #
  ptsStat_cex <- 2 * ry                                                        #
  ref_xlen <- 0.025 * rat * rx                                                 #
  ref_ylen <- 0.025 * rat * ry                                                 #
  tip_xlen <- 0.125 * rat * rx                                                 #
  tip_ylen <- 0.125 * rat * ry                                                 #
  ref_xlwd <- 4 * ry                                                           #
  ref_ylwd <- 4 * rx                                                           #
  arr_lwd <- 3 * mean(c(rx, ry))                                               #
  ax_lwd <- 2 * mean(c(rx, ry))                                                #
  scl_b <- 0.01 * mean(c(rx, ry))                                              #
                                                                               #
  tryCatch({                                                                   #
    pdf(out_pdf, width = round(30 / rat), height = 30, bg = "white")           #
    par(cex.main = mainleg_cex)                                                #
    par(cex.axis = labax_cex)                                                  #
    par(col.axis = "black")                                                    #
    par(cex.lab = labax_cex)                                                   #
                                                                               #
    #Main Plot                                                                 #
    image.plot(bkg, xlab =' ', ylab =' ', xaxt ="n", yaxt ="n",                #
               col = c("white", "#BBBBBB"),                                    #
               bigplot = c(.13, .13 + delt, .13, .93),                         #
               lab.breaks = rep(" ", 3), legend.shrink = 0)                    #
    #Arrow Field                                                               #
    # Getting uv projections, scaling to max_mag and dividing by two           #
    # because the projections will be centered on a given coordinate           #
    u <- 0.5 * temp_mag * cos(temp_ang * pi / 180) / max_mag                   #
    v <- 0.5 * temp_mag * sin(temp_ang * pi / 180) / max_mag                   #
                                                                               #
    # Scaling projections to plot grid                                         #
    us <- u * binsize / dimd[1]                                                #
    vs <- v * binsize / dimd[2]                                                #
                                                                               #
    scl_mag <- 0.5 * binsize * floor(max_mag) / max_mag / dimd[1]              #
                                                                               #
    arrows(xs - us, ys - vs, xs + us, ys + vs, length = 0, lwd = arr_lwd,      #
           col = a_col)                                                        #
                                                                               #
    #Points                                                                    #
    if(lpoints > 0){                                                           #
      for(p in 1:lpoints){                                                     #
        points(points[p, 1], points[p, 2], pch = chars[p], col = colors[p],    #
               cex = ptsStat_cex)                                              #
      }                                                                        #
    }                                                                          #
    #Main Title                                                                #
    mtext(mtitle, cex = tit_cex, col = "black", line = 5)                      #
    #N/E referencial                                                           #
    arrows(x0 = .15, y0 = 0.05, x1 = .15, y1 = 0.05 + ref_ylen,                #
           length = tip_ylen, col = "black", lwd = ref_ylwd)                   #
    arrows(x0 = .15, y0 = 0.05, x1 = .15 - ref_xlen, y1 = 0.05,                #
           length = tip_xlen, col = "black", lwd = ref_xlwd)                   #
    text(x = .15, y = 0.05 + ref_ylen * 1.05, labels = "N", pos = 3,           #
         col = "black", cex = labax_cex)                                       #
    text(x = .15 - ref_xlen * 1.05, y = 0.05, labels = "E", pos = 2,           #
         col = "black", cex = labax_cex)                                       #
                                                                               #
    line_offs <- 0                                                             #
    if(stat_flag){                                                             #
      # Scale legend                                                           #
      arrows(x0 = .1 + delt - scl_mag, y0 = 0.05,                              #
             x1 = .1 + delt + scl_mag, y1 = 0.05, length = 0, col = "black",   #
             lwd = arr_lwd)                                                    #
      text(x = .1 + delt + scl_mag + scl_b, y = 0.05,                          #
           labels = paste0(floor(max_mag), "%"), pos = 4, col = "black",       #
           cex = labax_cex)                                                    #
      text(x = .1 + delt - scl_mag - scl_b, y = 0.05, labels = "scale",        #
           pos = 2, col = "black", cex = labax_cex)                            #
                                                                               #
      #Statistical Info                                                        #
      mtext(mag_stat_str, side = 1, line = 3, cex = ptsStat_cex, col = "black")#
      mtext(ang_stat_str, side = 1, line = 6, cex = ptsStat_cex, col = "black")#
      line_offs <- 7                                                           #
    }                                                                          #
    #X-axis                                                                    #
    axis(1, at = x_ticks, labels = paste0(x_labs, "°"), line = 1 + line_offs,  #
         lwd = ax_lwd, padj = .5)                                              #
    mtext("Ra (deg)", side = 1, cex = tit_cex, line = 7 + line_offs,           #
          col = "black")                                                       #
    #Y-axis                                                                    #
    axis(2, at = y_ticks, labels = paste0(y_labs, "°"), lwd = ax_lwd, line = 1)#
    mtext("Dec (deg)", side = 2, cex = tit_cex, line = 6, col = "black")       #
    dev.off()                                                                  #
  }, error = function(e) {                                                     #
    print("ERROR: Could not create an arrow plot for the given input.")        #
    print(paste("Error message:", e$message))                                  #
    return(e)                                                                  #
  }                                                                            #
  )                                                                            #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------# 

################################################################################
get_NED_coords <- function(obj, N){                                            #
                                                                               #
  py_run_string(paste0("ra, dec = sexD.get_NED_coord('", obj,"',", N,")"))     #
                                                                               #
  coords <- c(py$ra, py$dec)                                                   #
                                                                               #
  return(coords)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_largest_diam <- function(obj){                                         #
                                                                               #
  t_diam <- 0                                                                  #
  n <- 1                                                                       #
  r_diam <- NA                                                                 #
                                                                               #
  while(!is.na(t_diam)){                                                       #
                                                                               #
    py_run_string(paste0("diam = sexD.get_NED_major_axis('", obj, "',", n,")"))#
    t_diam <- py$diam                                                          #
    n <- n + 1                                                                 #
                                                                               #
    if(!is.na(t_diam) && (is.na(r_diam) || t_diam > r_diam)){                  #
      r_diam <- t_diam                                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(r_diam)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_largest_axis_ratio <- function(obj){                                   #
                                                                               #
  t_axRat <- 0                                                                 #
  n <- 1                                                                       #
  axRat <- NA                                                                  #
                                                                               #
  while(!is.na(t_axRat)){                                                      #
                                                                               #
    py_run_string(paste0("aRat = sexD.get_NED_axis_ratio('", obj, "',", n,")"))#
    t_axRat <- py$aRat                                                         #
    n <- n + 1                                                                 #
                                                                               #
    if(!is.na(t_axRat) && (is.na(axRat) || t_axRat > axRat) && t_axRat < 1){   #
      axRat <- t_axRat                                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(axRat)                                                                #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_smallest_diam <- function(obj){                                        #
                                                                               #
  t_diam <- 0                                                                  #
  n <- 1                                                                       #
  r_diam <- NA                                                                 #
                                                                               #
  while(!is.na(t_diam)){                                                       #
                                                                               #
    py_run_string(paste0("diam = sexD.get_NED_major_axis('", obj, "',", n,")"))#
                                                                               #
    t_diam <- py$diam                                                          #
    n <- n + 1                                                                 #
                                                                               #
    if(!is.na(t_diam) && (is.na(r_diam) || t_diam < r_diam)){                  #
      r_diam <- t_diam                                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(r_diam)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_newest_diam <- function(obj){                                          #
                                                                               #
  py_run_string(paste0("diam = sexD.get_NED_major_axis('", obj, "',", 1, ")")) #
                                                                               #
  r_diam <- py$diam                                                            #
                                                                               #
  return(r_diam)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_newest_posang <- function(obj){                                        #
                                                                               #
  PA <- 0                                                                      #
  i <- 0                                                                       #
  while((PA == 0 | is.na(PA)) & i < 100){                                      #
    py_run_string(paste0("PA = sexD.get_NED_position_angle('", obj,"',",i,")"))#
    PA <- py$PA                                                                #
    i <- i + 1                                                                 #
  }                                                                            #
                                                                               #
  if(i == 100){                                                                #
    PA <- NA                                                                   #
  }                                                                            #
                                                                               #
  return(PA)                                                                   #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
get_NED_median_diam <- function(obj){                                          #
                                                                               #
  t_diam <- 0                                                                  #
  n <- 1                                                                       #
  array_diam <- NULL                                                           #
                                                                               #
  while(!is.na(t_diam)){                                                       #
                                                                               #
    py_run_string(paste0("diam = sexD.get_NED_major_axis('", obj, "',", n,")"))#
    t_diam <- py$diam                                                          #
    n <- n + 1                                                                 #
                                                                               #
    if(!is.na(t_diam)){                                                        #
      array_diam <- append(array_diam, t_diam)                                 #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(median(array_diam))                                                   #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given as input the coordinates "ra" and "dec" in 'deg', and the width of the #
# field "Wra" and "Wdec" also in 'deg'; Returns an array of Gaia DR 3 stars    #
# within the input field. The array includes:                                  #
### * Gaia 'source_id'                                                         #
### * reference year in 'yr'                                                   #
### * RA and DEC in 'deg'                                                      #
### * Proper motions in RA and DEC in 'deg/yr'                                 #
### * Distance and its uncertainty in 'pc'                                     #
################################################################################
get_GAIA_stars_coords <- function(ra, dec, Wra, Wdec){                         #
                                                                               #
  lra <- length(ra)                                                            #
  ldec <- length(dec)                                                          #
  lWra <- length(Wra)                                                          #
  lWdec <- length(Wdec)                                                        #
  mas_2_deg <- 3600000                                                         #
                                                                               #
  if(lra != 1 | ldec != 1 | lWra != 1 | lWdec != 1){                           #
    print(paste0("ERROR: 'ra' has length ", lra, ", 'dec' has length ", ldec,  #
                 ", 'Wra' has length ", lWra, " and 'Wdec' has length ", lWdec,#
                 ". All must have length 1. Returning NULL."))                 #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  py_run_string("import astropy.units as u")                                   #
  py_run_string("from astropy.coordinates import SkyCoord")                    #
  py_run_string("from astroquery.gaia import Gaia")                            #
  py_run_string("Gaia.ROW_LIMIT = -1")                                         #
                                                                               #
  set_coords <- paste0("coord = SkyCoord(ra=", ra, ", dec=", dec, ", unit=(u.",#
                       "degree, u.degree), frame='icrs')")                     #
  py_run_string(set_coords)                                                    #
                                                                               #
  set_grid_Wra <- paste0("width = u.Quantity(", Wra,", u.deg)")                #
  set_grid_Wdec <- paste0("height = u.Quantity(", Wdec,", u.deg)")             #
  py_run_string(set_grid_Wra)                                                  #
  py_run_string(set_grid_Wdec)                                                 #
                                                                               #
  get_objs <- paste0("objs = Gaia.query_object_async(coordinate = coord, widt",#
                     "h = width, height = height)")                            #
  py_run_string(get_objs)                                                      #
                                                                               #
  py_run_string("objs_N = len(objs['ra'])")                                    #
                                                                               #
  N <- py$objs_N                                                               #
                                                                               #
  objs_params <- array(0, dim = c(N, 12),                                      #
                       dimnames = list(NULL, c("source_id", "ref_yr", "ra",    #
                                               "dec", "u_ra", "u_dec", "pm_ra",#
                                               "pm_dec", "u_pm_ra", "u_pm_dec",#
                                               "dist", "dist_unc")))           #
                                                                               #
  Plx <- py$objs['parallax'] / 1000                                            #
  uPlx <- py$objs['parallax_error'] / 1000                                     #
                                                                               #
  objs_params[,"source_id"] <- py$objs['SOURCE_ID']                            #
  objs_params[,"ref_yr"] <- py$objs['ref_epoch']                               #
  objs_params[,"ra"] <- py$objs['ra']                                          #
  objs_params[,"dec"] <- py$objs['dec']                                        #
  objs_params[,"u_ra"] <- py$objs['ra_error'] / mas_2_deg                      #
  objs_params[,"u_dec"] <- py$objs['dec_error'] / mas_2_deg                    #
  objs_params[,"pm_ra"] <- py$objs['pmra'] / mas_2_deg                         #
  objs_params[,"pm_dec"] <- py$objs['pmdec'] / mas_2_deg                       #
  objs_params[,"u_pm_ra"] <- py$objs['pmra_error'] / mas_2_deg                 #
  objs_params[,"u_pm_dec"] <- py$objs['pmdec_error'] / mas_2_deg               #
  objs_params[,"dist"] <- 1 / Plx                                              #
  objs_params[,"dist_unc"] <- 1 / Plx^2 * uPlx                                 #
                                                                               #
  rej_ind_ra <- which(is.na(objs_params[,"pm_ra"]), arr.ind = T)               #
  rej_ind_dec <- which(is.na(objs_params[,"pm_dec"]), arr.ind = T)             #
  rej_ind <- union(rej_ind_ra, rej_ind_dec)                                    #
                                                                               #
  stars_coords <- objs_params[-rej_ind,]                                       #
                                                                               #
  py_run_string("del(coord, width, height, objs, objs_N)")                     #
                                                                               #
  return(stars_coords)                                                         #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given an input 2d array, "stars", with parameters "x" and "y" in 'pix', and  #
# a path. "astroPath", to an astrometry '.new' file;                           #
# Assumes parameters "x" and "y" are in positions 1 and 2 of the second        #
# dimension of "stars", with the positions of the first dimension holding      #
# different objects.                                                           #
# Returns a 2d array, "coords", with star coordinates "ra" and "dec" in 'deg'  #
# reporting to J2000.                                                          #
################################################################################
convert_pixel_to_wcs <- function(stars, astroPath){                            #
  ldim <- length(dim(stars))                                                   #
                                                                               #
  if(ldim != 2){                                                               #
    print(paste0("ERROR: 'stars' has ", ldim, " dimensions instead of 2. ",    #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dims <- dim(stars)                                                           #
  Nstars <- dims[1]                                                            #
                                                                               #
  if(dims[2] != 2){                                                            #
    print(paste0("ERROR: 'stars' has ", dims[2], " columns instead of 2 (X, Y",#
                 "). Returning NULL."))                                        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  py_run_string("import astropy.units as u")                                   #
  py_run_string("from astropy.wcs import WCS")                                 #
  py_run_string("from astropy.io import fits")                                 #
                                                                               #
  load_cmd <- paste0("temp_fits = fits.open('", astroPath, "')")               #
  py_run_string(load_cmd)                                                      #
  py_run_string("wcs_header = WCS(temp_fits[0].header)")                       #
  py_run_string("temp_fits.close()")                                           #
                                                                               #
  coords <- array(NA, dim = c(Nstars, 2), dimnames = list(NULL, c("ra","dec")))#
                                                                               #
  py$xs <- stars[,1] - 1                                                       #
  py$ys <- stars[,2] - 1                                                       #
  py_run_string("coords = wcs_header.wcs_pix2world(xs, ys, 0)")                #
                                                                               #
  coords[,"ra"] <- py$coords[[1]]                                              #
  coords[,"dec"] <- py$coords[[2]]                                             #
                                                                               #
  py_run_string("del(coords, wcs_header, xs, ys)")                             #
                                                                               #
  return(coords)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given an input 2d array, "stars", with parameters "ra" and "dec" in 'deg',   #
# and a path. "astroPath", to an astrometry '.new' file;                       #
# Assumes parameters "ra" and "dec" are in positions 1 and 2 of the second     #
# dimension of "stars", with the positions of the first dimension holding      #
# different objects.                                                           #
# Returns a 2d array, "coords", with star coordinates "x" and "y" in 'pix'.    #
################################################################################
convert_wcs_to_pixel <- function(stars, astroPath){                            #
  ldim <- length(dim(stars))                                                   #
                                                                               #
  if(ldim != 2){                                                               #
    print(paste0("ERROR: 'stars' has ", ldim, " dimensions instead of 2. ",    #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dims <- dim(stars)                                                           #
  Nstars <- dims[1]                                                            #
                                                                               #
  if(dims[2] != 2){                                                            #
    print(paste0("ERROR: 'stars' has ", dims[2], " columns instead of 2 (RA, ",#
                 "DEC). Returning NULL."))                                     #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  py_run_string("import astropy.units as u")                                   #
  py_run_string("from astropy.wcs import WCS")                                 #
  py_run_string("from astropy.io import fits")                                 #
                                                                               #
  load_cmd <- paste0("temp_fits = fits.open('", astroPath, "')")               #
  py_run_string(load_cmd)                                                      #
  py_run_string("wcs_header = WCS(temp_fits[0].header)")                       #
  py_run_string("temp_fits.close()")                                           #
                                                                               #
  coords <- array(NA, dim = c(Nstars, 2), dimnames = list(NULL, c("x", "y")))  #
                                                                               #
  py$ras <- stars[,1]                                                          #
  py$decs <- stars[,2]                                                         #
  py_run_string("coords = wcs_header.wcs_world2pix(ras, decs, 0)")             #
                                                                               #
  coords[,"x"] <- py$coords[[1]] + 1                                           #
  coords[,"y"] <- py$coords[[2]] + 1                                           #
                                                                               #
  py_run_string("del(coords, wcs_header, ras, decs)")                          #
                                                                               #
  return(coords)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2d input array, "stars", with parameters "ra" and "dec" in 'deg',    #
# "pm_ra" and "pm_dec" in 'deg/yr', and the reference year for those data;     #
# Assumes parameters "ra", "dec", "u_ra", "u_dec", "pm_ra", "pm_dec",          #
# "u_pm_ra", "u_pm_dec" and "ref_yr" are in that same order in the second      #
# dimension, with the first dimension positions holding different objects;     #
# Returns a 2d array, "coords", with star coordinates "ra", "dec", "u_ra" and  #
# "u_dec" in 'deg' reporting to J2000.                                         #
################################################################################
convert_yr_to_J2000 <- function(stars){                                        #
  ldim <- length(dim(stars))                                                   #
                                                                               #
  if(ldim != 2){                                                               #
    print(paste0("ERROR: 'stars' has ", ldim, " dimensions instead of 2. ",    #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dims <- dim(stars)                                                           #
  Nstars <- dims[1]                                                            #
                                                                               #
  if(dims[2] != 9){                                                            #
    print(paste0("ERROR: 'stars' has ", dims[2], " columns instead of 9 (ref_",#
                 "yr, ra, dec, u_ra, u_dec, pm_ra, pm_dec, u_pm_ra and u_pm_d",#
                 "ec). Returning NULL."))                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  coords <- array(NA, dim = c(Nstars, 4),                                      #
                  dimnames = list(NULL, c("ra", "dec", "u_ra", "u_dec")))      #
                                                                               #
  for(n in 1:Nstars){                                                          #
    cd_yr <- as.POSIXct(paste0(round(stars[n, 1]), "-01-01"), tz = "UTC")      #
    elled_days <- as.numeric(difftime(cd_yr, "2000-01-01", units = "days"))    #
    elled_yrs <- elled_days / 365.25                                           #
                                                                               #
    coords[n, "ra"] <- stars[n, 2] - stars[n, 6] * elled_yrs                   #
    coords[n, "dec"] <- stars[n, 3] - stars[n, 7] * elled_yrs                  #
    coords[n, "u_ra"] <- sqrt(stars[n, 4]^2 + (stars[n, 8] * elled_yrs)^2)     #
    coords[n, "u_dec"] <- sqrt(stars[n, 5]^2 + (stars[n, 9] * elled_yrs)^2)    #
  }                                                                            #
  return(coords)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Creates a solid ellipse within "matrix", centered on 'x0,y0', with semi-axis #
# 'a' and 'b', and tilted by 'ang' degrees counter clockwise with relation to  #
############################# the horizontal axis. #############################
################################################################################
drawEllipse <- function(matrix, x0, y0, a, b, ang, fill = 1){                  #
                                                                               #
  dimM <- dim(matrix)                                                          #
  ldimM <- length(dimM)                                                        #
                                                                               #
  if(ldimM != 2){                                                              #
    print(paste0("ERROR: 'matrix' must be a 2-dimensional. Returning NULL."))  #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(length(x0) != 1 || length(y0) != 1 || length(a) != 1 ||                   #
     length(b) != 1 || length(ang) != 1 || length(fill) != 1){                 #
    print(paste0("ERROR: 'x0', 'y0', 'a', 'b', 'ang' and 'fill' must have len",#
                 "gth 1. Returning NULL."))                                    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  #x0, y0, a and b must be positive                                            #
  if(x0 <= 0 || y0 <= 0 || a <= 0 || b <= 0){                                  #
    print(paste0("ERROR: 'x0', 'y0', 'a' and 'b' must be positive. Returning ",#
                 "NULL."))                                                     #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  angpi <- ang / 180 * pi                                                      #
                                                                               #
  for(x in 1:dimM[1]){                                                         #
    for(y in 1:dimM[2]){                                                       #
      if((((x - x0) * cos(angpi) + (y - y0) * sin(angpi))^2 / a^2 +            #
         ((x0 - x) * sin(angpi) + (y - y0) * cos(angpi))^2 / b^2) <= 1){       #
        matrix[x, y] <- fill                                                   #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  return(matrix)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given the parameters p_max, wl_max and k, and an input vector x, returns the #
###### Serkowski law output y = par[1] * exp(-par[2] * log(x / par[3])^2) ######
################################################################################
Serkowski_law <- function(x, par){                                             #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 3){                                      #
    print("ERROR: par is not a length 3 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] * exp(-par[2] * log(x / par[3])^2)                               #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given the parameters A, B and P, and an input vector x, returns the sinosoid #
######################## output y = A + B * sin(C * x + P) #####################
################################################################################
sinosoid_model <- function(x, par){                                            #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 4){                                      #
    print("ERROR: par is not a length 4 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] + par[2] * sin(par[3] * x + par[4])                              #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given the parameters 1, 2, 3 and 4 and an input vector x, returns the #####
########## Power law output y = par[1] + par[2] * (x + par[3])^par[4] ##########
################################################################################
Power_law <- function(x, par){                                                 #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 4){                                      #
    print("ERROR: par is not a length 4 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] + par[2] * (x + par[3])^par[4]                                   #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given the parameters 1, 2, 3 and 4 and an input vector x, returns the #####
####### Exp law output y = par[1] + par[2] * exp^(x / par[3] + par[4])   #######
################################################################################
Exp_law <- function(x, par){                                                   #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 4){                                      #
    print("ERROR: par is not a length 4 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] + par[2] * exp(x / par[3] + par[4])                              #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
## Given the parameters 1, 2, 3, 4 and 5, and an input vector x, returns the ###
# Sig law output y = par[1] / (par[2] + par[3] * exp^(-(x / par[4]) + par[5])) #
################################################################################
Sig_law <- function(x, par){                                                   #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 5){                                      #
    print("ERROR: par is not a length 5 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] / (par[2] + par[3] * exp(-(x / par[4]) + par[5]))                #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given the parameters 1 and 2, and an input vector x, returns the linear ###
##################### model output y = par[1] * x + par[2] #####################
################################################################################
linear_model <- function(x, par){                                              #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 2){                                      #
    print("ERROR: par is not a length 2 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] * x + par[2]                                                     #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Given the parameters 1, 2, 3, 4 and 5 and an input vector x, returns the ###
### Sigmoid output y = par[1] + par[2] / (par[3] + par[4] * exp(x / par[5])) ###
################################################################################
Sigmoid_model <- function(x, par){                                             #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 5){                                      #
    print("ERROR: par is not a length 5 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] + par[2] / (par[3] + par[4] * exp(x / par[5]))                   #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given the parameters 1, 2, 3 and 4 and an input vector x, returns the #####
######## rational law output y = par[1] + par[2] / (par[3] + x)^par[4] #########
################################################################################
Rational_model <- function(x, par){                                            #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 4){                                      #
    print("ERROR: par is not a length 4 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] + par[2] / (par[3] + x)^par[4]                                   #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#### Given the parameters 1, 2, 3 and 4 and an input vector x, returns the #####
######### vanish law output y = par[1] / (par[2] * x^par[3] + par[4]) ##########
################################################################################
Vanish_model <- function(x, par){                                              #
  if(!is.vector(x)){                                                           #
    print("ERROR: x is not a vector. Returning NULL.")                         #
    return(NULL)                                                               #
  }                                                                            #
  if(!is.vector(par) | length(par) != 4){                                      #
    print("ERROR: par is not a length 4 vector. Returning NULL.")              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  y <- par[1] / (par[2] * x^par[3] + par[4])                                   #
  return(y)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Given a formula 'form', a data dependent vector 'y', a data feature 'x', ###
### their respective uncertainties 'unc_y' and 'unc_x', and a set of control ###
## parameters, this function returns the mean of 'fit_N' models obtained with ##
#################### nls() and random draws from the data, #####################
################## and the 68% and 95% confidence boundaries. ##################
################################################################################
bootstrap_nls <- function(form, y, x, unc_y = 0, unc_x = 0, fit_N = 100,       #
                          fits_per_draw_N = 1, forceN = F, rnd_method = 'norm',#
                          par_up_lim = NA, par_lo_lim = NA, manual_ini = F,    #
                          par_ini = NA, par_ini_mean = NA, par_ini_width = NA, #
                          nls_ctrl = nls.control(maxiter = 500, tol = 1e-08,   #
                                                 warnOnly = T),                #
                          nls_method = 'port', use_all = F, max_j = 100,       #
                          minN_corr = .85, maxN_res = .15,                     #
                          CI_1sd = T, CI_2sd = T, verbose = F){                #
                                                                               #
  # Test 'form' class                                                          #
  if(class(form) != "formula"){                                                #
    print(paste0("ERROR: 'form' must be of class 'formula'. Instead, it is of",#
                 " class ", class(form), ". Returning NULL."))                 #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  form_str <- paste0(paste0(form)[2]," ", paste0(form)[1]," ", paste0(form)[3])#
  form_vars <- all.vars(form)                                                  #
                                                                               #
  cmd <- tail(as.character(form), 1)                                           #
  exp_cmd <- parse(text = cmd)                                                 #
  model_call <- function(...) eval(exp_cmd, list(...))                         #
                                                                               #
  # Test name of dependent variable                                            #
  if(form_vars[1] != "Y"){                                                     #
    print(paste0("ERROR: 'form' must include a dependent variable named 'Y'.", #
                 " Returning NULL."))                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test name of independent variable                                          #
  if(length(which(form_vars == "X")) != 1){                                    #
    print(paste0("ERROR: 'form' must include one independent variable named ", #
                 "'X'. Returning NULL."))                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  par_list <- form_vars[-c(which(form_vars == "X"), which(form_vars == "Y"))]  #
  len_par <- length(par_list)                                                  #
                                                                               #
  model_par_list <- c(par_list, "data_corr", "new_corr", "data_res(%)",        #
                      "new_res(%)", "data_chisq", "new_chisq", "data_p-val",   #
                      "new_p-val", "data_BIC", "new_BIC")                      #
  model_par_len <- length(model_par_list)                                      #
                                                                               #
  # Test length of y = x                                                       #
  len_data <- length(y)                                                        #
  len_x <- length(x)                                                           #
                                                                               #
  if(len_data != len_x){                                                       #
    print(paste0("ERROR: length of 'x' must be equal to the length of 'y'. ",  #
                 "'y' has length ", len_data, " and 'x' has length ", len_x,   #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test length of unc_y = y or 1                                              #
  len_unc_y <- length(unc_y)                                                   #
                                                                               #
  if(len_data != len_unc_y & len_unc_y != 1){                                  #
    print(paste0("ERROR: length of 'unc_y' must be either 1 or equal to the ", #
                 "length of 'y'. 'y' has length ", len_data, " and 'unc_y' ",  #
                 "has length ", len_unc_y, ". Returning NULL."))               #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(len_unc_y == 1){                                                          #
    unc_y <- rep(unc_y, len_data)                                              #
  }                                                                            #
                                                                               #
  # Test length of unc_x = x or 1                                              #
  len_unc_x <- length(unc_x)                                                   #
                                                                               #
  if(len_data != len_unc_x & len_unc_x != 1){                                  #
    print(paste0("ERROR: length of 'unc_x' must be either 1 or equal to the ", #
                 "length of 'x'. 'x' has length ", len_data, " and 'unc_x' ",  #
                 "has length ", len_unc_x, ". Returning NULL."))               #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(len_unc_x == 1){                                                          #
    unc_x <- rep(unc_x, len_data)                                              #
  }                                                                            #
                                                                               #
  # Filtering NAs from 'y' and 'x'                                             #
  x_NAs <- which(is.na(x), arr.ind = T)                                        #
  x_NAs <- c(x_NAs, which(is.infinite(x), arr.ind = T))                        #
  x_NAs <- c(x_NAs, which(is.null(x), arr.ind = T))                            #
                                                                               #
  y_NAs <- which(is.na(y), arr.ind = T)                                        #
  y_NAs <- c(y_NAs, which(is.infinite(y), arr.ind = T))                        #
  y_NAs <- c(y_NAs, which(is.null(y), arr.ind = T))                            #
                                                                               #
  inds_to_rm <- union(x_NAs, y_NAs)                                            #
  if(length(inds_to_rm) != 0){                                                 #
    x <- x[-inds_to_rm]                                                        #
    y <- y[-inds_to_rm]                                                        #
    unc_x <- unc_x[-inds_to_rm]                                                #
    unc_y <- unc_y[-inds_to_rm]                                                #
  }                                                                            #
                                                                               #
  len_data <- length(y)                                                        #
                                                                               #
  # Replacing NAs from 'unc_y' and 'unc_x'                                     #
  x_NAs <- which(is.na(unc_x), arr.ind = T)                                    #
  x_NAs <- c(x_NAs, which(is.infinite(unc_x), arr.ind = T))                    #
  x_NAs <- c(x_NAs, which(is.null(unc_x), arr.ind = T))                        #
  unc_x[x_NAs] <- 0                                                            #
                                                                               #
  y_NAs <- which(is.na(unc_y), arr.ind = T)                                    #
  y_NAs <- c(y_NAs, which(is.infinite(unc_y), arr.ind = T))                    #
  y_NAs <- c(y_NAs, which(is.null(unc_y), arr.ind = T))                        #
  unc_y[y_NAs] <- 0                                                            #
                                                                               #
  # Test degrees of freedom                                                    #
  df <- len_data - len_par                                                     #
                                                                               #
  if(df < 0){                                                                  #
    print(paste0("ERROR: negative degrees of freedom. Increase number of data",#
                 " points or choose a model with less parameters. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'fit_N'                                                               #
  if(length(fit_N) != 1){                                                      #
    print(paste0("WARNING: 'fit_N' is length ", length(fit_N), ". Coercing ",  #
                 "use of first element."))                                     #
    fit_N <- fit_N[1]                                                          #
  }                                                                            #
  if((fit_N - floor(fit_N)) != 0 | fit_N == 0){                                #
    print("ERROR: 'fit_N' must be a positive integer. Returning NULL.")        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'max_j'                                                               #
  if(length(max_j) != 1){                                                      #
    print(paste0("WARNING: 'max_j' is length ", length(max_j), ". Coercing ",  #
                 "use of first element."))                                     #
    max_j <- max_j[1]                                                          #
  }                                                                            #
  if((max_j - floor(max_j)) != 0 | max_j == 0){                                #
    print("ERROR: 'max_j' must be a positive integer. Returning NULL.")        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'fits_per_draw_N'                                                     #
  if(length(fits_per_draw_N) != 1){                                            #
    print(paste0("WARNING: 'fits_per_draw_N' is length ",                      #
                 length(fits_per_draw_N), ". Coercing use of first element.")) #
    fits_per_draw_N <- fits_per_draw_N[1]                                      #
  }                                                                            #
  if((fits_per_draw_N - floor(fits_per_draw_N)) != 0 | fits_per_draw_N == 0){  #
    print(paste0("ERROR: 'fits_per_draw_N' must be a positive integer. Return",#
                 "ing NULL."))                                                 #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  while(manual_ini & fits_per_draw_N != 1){                                    #
    print(paste0("WARNING: with manual initialization of parameters ",         #
                 "'fits_per_draw_N' must be set to 1. Otherwise you'd be ",    #
                 "multiple equal fits which would hurt the statistics."))      #
    print("Change or confirm 'manual_ini': ")                                  #
    manual_ini <- edit(manual_ini)                                             #
                                                                               #
    print("Change or confirm 'fits_per_draw_N': ")                             #
    fits_per_draw_N <- edit(fits_per_draw_N)                                   #
                                                                               #
    cancel <- FALSE                                                            #
    print("If you'd rather cancel the execution enter 'TRUE':")                #
    cancel <- edit(cancel)                                                     #
    if(cancel){                                                                #
      return(NULL)                                                             #
    }                                                                          #
  }                                                                            #
  draw_N <- fit_N %/% fits_per_draw_N                                          #
                                                                               #
  # Test 'rnd_method'                                                          #
  valid_rnd <- c('unif', 'norm', 'pois', 'exp')                                #
  if(sum(rnd_method == valid_rnd) != 1){                                       #
    print(paste0("ERROR: invalid choice of distribution from which to draw ",  #
                 "data samples. Valid options are 'unif', 'norm', 'pois' and", #
                 " 'exp'. Returning NULL."))                                   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'nls_method'                                                          #
  valid_nls <- c('default', 'plinear', 'port')                                 #
  if(sum(nls_method == valid_nls) != 1){                                       #
    print(paste0("ERROR: invalid choice of fitting method. Valid options are ",#
                 "'default', 'plinear' and 'port'. Returning NULL."))          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test length of par_up_lim and par_lo_lim = length par_list                 #
  # And initialize limits if so required                                       #
  if(prod(!is.na(par_up_lim))){                                                #
    if(length(par_up_lim) != len_par){                                         #
      if(length(par_up_lim) > len_par){                                        #
        print(paste0("WARNING: length of upper limits for model parameters is",#
                     "larger than the number of model parameters. Ignoring ",  #
                     "limits in excess."))                                     #
        par_up_lim <- par_up_lim[1:len_par]                                    #
      }else{                                                                   #
        print(paste0("WARNING: length of upper limits for model parameters is",#
                     "smaller than the number of model parameters. Remaining ",#
                     "limits will be set to Inf."))                            #
        par_dif_len <- len_par - length(par_up_lim)                            #
        par_up_lim <- c(par_up_lim, rep(Inf, par_dif_len))                     #
      }                                                                        #
    }                                                                          #
  }else{                                                                       #
    if(nls_method == 'port'){                                                  #
      print(paste0("WARNING: 'port' method requires the definition of upper ", #
                   "limits for the model parameters. Coercing these limits to",#
                   " Inf."))                                                   #
      par_up_lim <- rep(Inf, len_par)                                          #
    }                                                                          #
  }                                                                            #
  if(prod(!is.na(par_lo_lim))){                                                #
    if(length(par_lo_lim) != len_par){                                         #
      if(length(par_lo_lim) > len_par){                                        #
        print(paste0("WARNING: length of lower limits for model parameters is",#
                     "larger than the number of model parameters. Ignoring ",  #
                     "limits in excess."))                                     #
        par_lo_lim <- par_lo_lim[1:len_par]                                    #
      }else{                                                                   #
        print(paste0("WARNING: length of lower limits for model parameters is",#
                     "smaller than the number of model parameters. Remaining ",#
                     "limits will be set to -Inf."))                           #
        par_dif_len <- len_par - length(par_lo_lim)                            #
        par_lo_lim <- c(par_lo_lim, rep(-Inf, par_dif_len))                    #
      }                                                                        #
    }                                                                          #
  }else{                                                                       #
    if(nls_method == 'port'){                                                  #
      print(paste0("WARNING: 'port' method requires the definition of lower ", #
                   "limits for the model parameters. Coercing these limits to",#
                   " -Inf."))                                                  #
      par_lo_lim <- rep(-Inf, len_par)                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Test parameter initialization                                              #
  if(manual_ini){                                                              #
    if(length(par_ini) > len_par){                                             #
      print(paste0("WARNING: length of 'par_ini' is larger than the number of",#
                   " model parameters. Ignoring those in excess."))            #
      par_ini <- par_ini[1:len_par]                                            #
    }                                                                          #
    if(length(par_ini) < len_par){                                             #
      print(paste0("ERROR: length of 'par_ini' is smaller than the number of ",#
                   "model parameters. Returning NULL."))                       #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    start_list <- list()                                                       #
                                                                               #
    for(p in 1:len_par){                                                       #
      start_list[[p]] <- NA                                                    #
    }                                                                          #
                                                                               #
    names(start_list) <- par_list                                              #
                                                                               #
    for(p in 1:len_par){                                                       #
      start_list[[p]] <- par_ini[p]                                            #
    }                                                                          #
  }else{                                                                       #
    if(length(par_ini_mean) > len_par){                                        #
      print(paste0("WARNING: length of 'par_ini_mean' is larger than the numb",#
                   "er of model parameters. Ignoring those in excess."))       #
      par_ini_mean <- par_ini_mean[1:len_par]                                  #
    }                                                                          #
    if(length(par_ini_mean) < len_par){                                        #
      print(paste0("ERROR: length of 'par_ini_mean' is smaller than the ",     #
                   "number of model parameters. Returning NULL."))             #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    if(prod(is.na(par_ini_width))){                                            #
      if(fits_per_draw_N != 1){                                                #
        print(paste0("WARNING: 'par_ini_width' was not provided while the ",   #
                     "number of fits per sample requested is greater than 1. ",#
                     "Coercing width of parameter initialization to ",         #
                     "sqrt(abs('par_ini_mean'))."))                            #
        par_ini_width <- sqrt(abs(par_ini_mean))                               #
      }                                                                        #
    }else{                                                                     #
      if(length(par_ini_width) > len_par){                                     #
        print(paste0("WARNING: length of 'par_ini_width' is larger than the",  #
                     " number of model parameters. Ignoring those in excess."))#
        par_ini_width <- par_ini_width[1:len_par]                              #
      }                                                                        #
      if(length(par_ini_width) < len_par){                                     #
        print(paste0("WARNING: length of 'par_ini_width' is smaller than the", #
                     " number of model parameters. Remaining values will be ", #
                     "set to sqrt(abs('par_ini_mean'))."))                     #
        par_ini_width <- c(par_ini_width,                                      #
                           sqrt(abs(par_ini_mean[-(1:len_par)])))              #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Check if user is sure to force all fits to be successful                   #
  if(forceN){                                                                  #
    print(paste0("WARNING: forcing all samples fit to succeed may result in a",#
                 " longer execution time, it may also lead to an infinite ",   #
                 "loop. If 'manual_ini' was activated there is also a risk of",#
                 " error in execution to prevent such an infinity loop, due ", #
                 "to not being possible to vary the initial values."))         #
  }                                                                            #
                                                                               #
                                                                               #
  if(!use_all){                                                                #
    # Test min_corr                                                            #
    if(length(minN_corr) != 1){                                                #
      print(paste0("WARNING: 'minN_corr' has length greater than 1. Using ",   #
                   "only the first value."))                                   #
      minN_corr <- minN_corr[1]                                                #
    }                                                                          #
    if(minN_corr <= -1 | minN_corr >= 1){                                      #
      print(paste0("WARNING: 'minN_corr' must be greater than -1 and smaller ",#
                   "than 1. Coercing to 0.85."))                               #
      minN_corr <- 0.85                                                        #
    }                                                                          #
                                                                               #
    # Test max_res                                                             #
    if(length(maxN_res) != 1){                                                 #
      print(paste0("WARNING: 'maxN_res' has length greater than 1. Using only",#
                   " the first value."))                                       #
      maxN_res <- maxN_res[1]                                                  #
    }                                                                          #
    if(maxN_res <= 0 | maxN_res >= 1){                                         #
      print(paste0("WARNING: 'maxnN_res' must be greater than 0 and smaller ", #
                   "than 1. Coercing to 0.15."))                               #
      maxN_res <- 0.15                                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Initialize list of models                                                  #
  models <- list()                                                             #
  models_par <- array(NA, dim = c(len_par, fit_N),                             #
                      dimnames = list(par_list, NULL))                         #
  model_par_stat <- array(NA, dim = c(model_par_len, 2),                       #
                          dimnames = list(model_par_list, c("mean", "sd")))    #
                                                                               #
  sampled_Y <- array(NA, dim = c(len_data, draw_N))                            #
  sampled_X <- array(NA, dim = c(len_data, draw_N))                            #
  for(d in 1:len_data){                                                        #
    if(draw_N != 1){                                                           #
      switch(rnd_method,                                                       #
             'unif' = c(sampled_Y[d,] <- runif(draw_N, y[d] - unc_y[d],        #
                                               y[d] + unc_y[d]),               #
                        sampled_X[d,] <- runif(draw_N, x[d] - unc_x[d],        #
                                               x[d] + unc_x[d])),              #
             'norm' = c(sampled_Y[d,] <- rnorm(draw_N, y[d], unc_y[d]),        #
                        sampled_X[d,] <- rnorm(draw_N, x[d], unc_x[d])),       #
             'pois' = c(sampled_Y[d,] <- rpois(draw_N,                         #
                                               rnorm(1, y[d], unc_y[d])),      #
                        sampled_X[d,] <- rpois(draw_N,                         #
                                               rnorm(1, x[d],unc_x[d]))),      #
             'exp' = c(sampled_Y[d,] <- rexp(draw_N,                           #
                                             1 / rnorm(1,y[d],unc_y[d])),      #
                       sampled_X[d,] <- rexp(draw_N,                           #
                                             1 / rnorm(1,x[d],unc_x[d])))      #
      )                                                                        #
    }else{                                                                     #
      sampled_Y[d,] <- y[d]                                                    #
      sampled_X[d,] <- x[d]                                                    #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Initialize loop of sampling and fitting                                    #
  for(nd in 1:draw_N){                                                         #
                                                                               #
    # Draw data sample to be fitted                                            #
    data_fr <- NULL                                                            #
    data_fr$Y <- sampled_Y[,nd]                                                #
    data_fr$X <- sampled_X[,nd]                                                #
                                                                               #
    for(nf in 1:fits_per_draw_N){                                              #
                                                                               #
      ns <- fits_per_draw_N * (nd - 1) + nf                                    #
                                                                               #
      # Fitting model                                                          #
      print(paste0("Fitting ", form_str, " model ", ns, "/", fit_N, ":"))      #
                                                                               #
      EorW_flag <- TRUE                                                        #
      j <- 1                                                                   #
                                                                               #
      keep_going <- FALSE                                                      #
      if(forceN){                                                              #
        keep_going <- TRUE                                                     #
      }                                                                        #
                                                                               #
      while((EorW_flag && j < max_j) | keep_going){                            #
        EorW_flag <- FALSE                                                     #
                                                                               #
        # Initialize model parameters                                          #
        if(!manual_ini){                                                       #
          start_list <- list()                                                 #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- NA                                              #
          }                                                                    #
                                                                               #
          names(start_list) <- par_list                                        #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- runif(1,                                        #
                                     par_ini_mean[p] - par_ini_width[p] / 2,   #
                                     par_ini_mean[p] + par_ini_width[p] / 2)   #
          }                                                                    #
        }                                                                      #
                                                                               #
        models[[ns]] <- tryCatch(                                              #
          nls(form, data = data_fr, start = start_list, control = nls_ctrl,    #
              lower = par_lo_lim, upper = par_up_lim,                          #
              algorithm = nls_method),                                         #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
                                                                               #
        if(is.logical(models[[ns]])){                                          #
          EorW_flag <- TRUE                                                    #
          j <- j + 1                                                           #
          if(verbose){                                                         #
            cat(".")                                                           #
          }                                                                    #
        }else{                                                                 #
          keep_going <- FALSE                                                  #
        }                                                                      #
      }                                                                        #
      if(!EorW_flag){                                                          #
        cat("!!Success!!\n")                                                   #
        print(paste0("--> ", form_str, " model ", ns, "/", fit_N, " has bee",  #
                     "n fitted..."))                                           #
      }else{                                                                   #
        cat("!!Failure!!\n")                                                   #
        print(paste0("--> After ", max_j, " attempts, ", form_str, " model ",  #
                     ns, "/", fit_N, " was not fitted..."))                    #
      }                                                                        #
                                                                               #
      # Saving model parameters                                                #
      if(!is.logical(models[[ns]])){                                           #
        model_par_list <- list()                                               #
        for(p in 1:len_par){                                                   #
          model_par_list[[p]] <- assign(par_list[p],                           #
                                        models[[ns]]$m$getPars()[p])           #
        }                                                                      #
        Y_tilda <- model_call(X = x, model_par_list)                           #
                                                                               #
        temp_cor <- cor(data_fr$Y, Y_tilda, use ="pairwise.complete.obs")      #
        temp_res <- median(abs((data_fr$Y - Y_tilda) / data_fr$Y))             #
                                                                               #
        if(is.na(temp_cor)){                                                   #
          temp_cor <- 0                                                        #
        }                                                                      #
        if(is.na(temp_res)){                                                   #
          temp_res <- 1                                                        #
        }                                                                      #
                                                                               #
        if(!use_all & temp_cor >= minN_corr & temp_res <= maxN_res){           #
          models_par[,ns] <- models[[ns]]$m$getPars()                          #
        }else{                                                                 #
          if(use_all){                                                         #
            models_par[,ns] <- models[[ns]]$m$getPars()                        #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(fit_N > 1){                                                               #
    # Estimate average separation between a pair of 'x'                        #
    delta_x <- max(x + unc_x) - min(x - unc_x)                                 #
    avg_Sep <- delta_x / len_data                                              #
                                                                               #
    # Create a new 'x' sequence                                                #
    new_avg_Sep <- avg_Sep / 100                                               #
    new_x <- seq(min(x - 3 * unc_x), max(x + 3 * unc_x), by = new_avg_Sep)     #
    len_new_data <- length(new_x)                                              #
    new_df <- len_new_data - len_par                                           #
                                                                               #
    models_pred <- array(NA, dim = c(fit_N, len_new_data))                     #
                                                                               #
    # Make predictions using every successfully fitted model                   #
    for(ns in 1:fit_N){                                                        #
      if(!is.atomic(models[[ns]])){                                            #
                                                                               #
        model_par_list <- list()                                               #
        for(p in 1:len_par){                                                   #
          model_par_list[[p]] <- assign(par_list[p],                           #
                                        models[[ns]]$m$getPars()[p])           #
        }                                                                      #
        models_pred[ns,] <- model_call(X = new_x, model_par_list)              #
      }                                                                        #
    }                                                                          #
                                                                               #
    # Estimate median and sd of model parameters to initialize new models      #
    avg_model_par <- apply(models_par, 1, median, na.rm = T)                   #
    sd_model_par <- apply(models_par, 1, sd, na.rm = T)                        #
                                                                               #
    # Estimate median and sd of predictions                                    #
    avg_models_pred <- apply(models_pred, 2, median, na.rm = T)                #
    sd_models_pred <- apply(models_pred, 2, sd)                                #
                                                                               #
    # Fit model to median of predictions to get median model                   #
    EorW_flag <- TRUE                                                          #
    j <- 1                                                                     #
                                                                               #
    data_fr$Y <- avg_models_pred                                               #
    data_fr$X <- new_x                                                         #
                                                                               #
    print(paste0("Fitting median ", form_str, " model:"))                      #
                                                                               #
    while(EorW_flag && j < max_j){                                             #
      EorW_flag <- FALSE                                                       #
                                                                               #
      # Initialize model parameters                                            #
      start_list <- list()                                                     #
                                                                               #
      for(p in 1:len_par){                                                     #
        start_list[[p]] <- NA                                                  #
      }                                                                        #
                                                                               #
      names(start_list) <- par_list                                            #
                                                                               #
      for(p in 1:len_par){                                                     #
        start_list[[p]] <- rnorm(1, avg_model_par[p], sd_model_par[p])         #
      }                                                                        #
                                                                               #
      # Fit model                                                              #
      avg_model <- tryCatch(                                                   #
        nls(form, data = data_fr, start = start_list, control = nls_ctrl,      #
            lower = par_lo_lim, upper = par_up_lim, algorithm = nls_method),   #
        error = function(e){                                                   #
          return(TRUE)                                                         #
        },                                                                     #
        warning = function(w){                                                 #
          return(TRUE)                                                         #
        }                                                                      #
      )                                                                        #
      if(is.logical(avg_model)){                                               #
        EorW_flag <- TRUE                                                      #
        j <- j + 1                                                             #
        if(verbose){                                                           #
          cat(".")                                                             #
        }                                                                      #
      }                                                                        #
    }                                                                          #
    if(!EorW_flag){                                                            #
      cat("!!Success!!\n")                                                     #
      print(paste0("--> Median ", form_str, " model has been fitted..."))      #
    }else{                                                                     #
      cat("!!Failure!!\n")                                                     #
      print(paste0("--> After ", max_j, " attempts, median ", form_str,        #
                   " model was not fitted. Skipping Conf. Int. fits."))        #
    }                                                                          #
                                                                               #
    # Saving model parameters                                                  #
    if(!is.logical(avg_model)){                                                #
                                                                               #
      model_par_list <- list()                                                 #
      for(p in 1:len_par){                                                     #
        model_par_list[[p]] <- assign(par_list[p], avg_model$m$getPars()[p])   #
      }                                                                        #
                                                                               #
      avg_model_pred <- model_call(X = new_x, model_par_list)                  #
      avg_data_pred <- model_call(X = x, model_par_list)                       #
                                                                               #
      new_corr <- cor(data_fr$Y, avg_model_pred, use = "pairwise.complete.obs")#
      data_corr <- cor(y, avg_data_pred, use = "pairwise.complete.obs")        #
                                                                               #
      new_res <- median((abs((data_fr$Y - avg_model_pred) / data_fr$Y) * 100), #
                        na.rm = T)                                             #
      data_res <- median((abs((y - avg_data_pred) / y) * 100), na.rm = T)      #
                                                                               #
      new_chisq <- sum((data_fr$Y - avg_model_pred)^2 / data_fr$Y, na.rm=T)    #
      data_chisq <- sum((y - avg_data_pred)^2 / y, na.rm=T)                    #
                                                                               #
      new_BIC <- log(sum((data_fr$Y - avg_model_pred)^2) / new_df) + len_par * #
        log(len_new_data)                                                      #
      data_BIC <- log(sum((y - avg_data_pred)^2) / df) + len_par *             #
        log(len_data)                                                          #
                                                                               #
      model_par_stat[1:len_par, "mean"] <- avg_model$m$getPars()               #
                                                                               #
      vcov_mat <- tryCatch(                                                    #
        vcov(avg_model),                                                       #
        error = function(e){                                                   #
          return(TRUE)                                                         #
        }                                                                      #
      )                                                                        #
                                                                               #
      if(!is.logical(vcov_mat)){                                               #
        se <- sqrt(diag(vcov_mat))                                             #
        alpha <- 1 - 0.682                                                     #
        z_value <- qnorm(1 - alpha / 2)                                        #
                                                                               #
        model_par_stat[1:len_par, "sd"] <- z_value * se                        #
        rm(vcov_mat, se, alpha, z_value)                                       #
      }else{                                                                   #
        model_par_stat[1:len_par, "sd"] <- NA                                  #
      }                                                                        #
                                                                               #
      model_par_stat["data_corr", "mean"] <- data_corr                         #
      model_par_stat["new_corr", "mean"] <- new_corr                           #
      model_par_stat["data_res(%)", "mean"] <- data_res                        #
      model_par_stat["new_res(%)", "mean"] <- new_res                          #
      model_par_stat["data_chisq", "mean"] <- data_chisq                       #
      model_par_stat["new_chisq", "mean"] <- new_chisq                         #
      model_par_stat["data_p-val", "mean"] <- 1 - pchisq(data_chisq, df)       #
      model_par_stat["new_p-val", "mean"] <- 1 - pchisq(new_chisq, new_df)     #
      model_par_stat["data_BIC", "mean"] <- data_BIC                           #
      model_par_stat["new_BIC", "mean"] <- new_BIC                             #
    }                                                                          #
                                                                               #
    # Fit model to median +- sd of predictions to get CI 68% of avg_model      #
    if(!is.atomic(avg_model)){                                                 #
      if(CI_1sd){                                                              #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        data_fr$Y <- avg_models_pred + sd_models_pred                          #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          start_list <- list()                                                 #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- NA                                              #
          }                                                                    #
                                                                               #
          names(start_list) <- par_list                                        #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- rnorm(1, avg_model_par[p], sd_model_par[p])     #
          }                                                                    #
                                                                               #
          # Fit upper bound model                                              #
          sd1u_model <- tryCatch(                                              #
            nls(form, data = data_fr, start = start_list, control = nls_ctrl,  #
                lower = par_lo_lim, upper = par_up_lim,                        #
                algorithm = nls_method),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
          if(is.logical(sd1u_model)){                                          #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Upper bound of 68% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, upper bound of 68% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        data_fr$Y <- avg_models_pred - sd_models_pred                          #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          start_list <- list()                                                 #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- NA                                              #
          }                                                                    #
                                                                               #
          names(start_list) <- par_list                                        #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- rnorm(1, avg_model_par[p], sd_model_par[p])     #
          }                                                                    #
                                                                               #
          # Fit lower bound model                                              #
          sd1l_model <- tryCatch(                                              #
            nls(form, data = data_fr, start = start_list, control = nls_ctrl,  #
                lower = par_lo_lim, upper = par_up_lim,                        #
                algorithm = nls_method),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
          if(is.logical(sd1l_model)){                                          #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Lower bound of 68% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, lower bound of 68% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
      }                                                                        #
                                                                               #
      if(CI_2sd){                                                              #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        data_fr$Y <- avg_models_pred + 2 * sd_models_pred                      #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          start_list <- list()                                                 #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- NA                                              #
          }                                                                    #
                                                                               #
          names(start_list) <- par_list                                        #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- rnorm(1, avg_model_par[p], sd_model_par[p])     #
          }                                                                    #
                                                                               #
          # Fit upper bound model                                              #
          sd2u_model <- tryCatch(                                              #
            nls(form, data = data_fr, start = start_list, control = nls_ctrl,  #
                lower = par_lo_lim, upper = par_up_lim,                        #
                algorithm = nls_method),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
          if(is.logical(sd2u_model)){                                          #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Upper bound of 95% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, upper bound of 95% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        data_fr$Y <- avg_models_pred - 2 * sd_models_pred                      #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          start_list <- list()                                                 #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- NA                                              #
          }                                                                    #
                                                                               #
          names(start_list) <- par_list                                        #
                                                                               #
          for(p in 1:len_par){                                                 #
            start_list[[p]] <- rnorm(1, avg_model_par[p], sd_model_par[p])     #
          }                                                                    #
                                                                               #
          # Fit lower bound model                                              #
          sd2l_model <- tryCatch(                                              #
            nls(form, data = data_fr, start = start_list, control = nls_ctrl,  #
                lower = par_lo_lim, upper = par_up_lim,                        #
                algorithm = nls_method),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
          if(is.logical(sd2l_model)){                                          #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Lower bound of 95% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, lower bound of 95% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    output_list <- list()                                                      #
    output_list[[1]] <- avg_model                                              #
    output_list[[2]] <- model_par_stat                                         #
                                                                               #
    if(CI_1sd && !is.atomic(avg_model)){                                       #
      output_list[[length(output_list) + 1]] <- sd1u_model                     #
      output_list[[length(output_list) + 1]] <- sd1l_model                     #
    }                                                                          #
    if(CI_2sd && !is.atomic(avg_model)){                                       #
      output_list[[length(output_list) + 1]] <- sd2u_model                     #
      output_list[[length(output_list) + 1]] <- sd2l_model                     #
    }                                                                          #
  }else{                                                                       #
    if(!is.logical(models[[1]])){                                              #
                                                                               #
      model_par_stat[1:len_par, "mean"] <- models[[1]]$m$getPars()             #
      model_par_list <- list()                                                 #
      for(p in 1:len_par){                                                     #
        model_par_list[[p]] <- assign(par_list[p], models[[1]]$m$getPars()[p]) #
      }                                                                        #
                                                                               #
      data_pred <- model_call(X = x, model_par_list)                           #
      model_par_stat["data_corr", "mean"] <-                                   #
        cor(y, data_pred, use = "pairwise.complete.obs")                       #
      model_par_stat["data_res(%)", "mean"] <-                                 #
        median((abs((y - data_pred) / y) * 100), na.rm = T)                    #
      model_par_stat["data_chisq", "mean"] <-                                  #
        sum((y - data_pred)^2 / y, na.rm = T)                                  #
      model_par_stat["data_p-val", "mean"] <-                                  #
        1 - pchisq(model_par_stat["data_chisq", "mean"], df)                   #
      model_par_stat["data_BIC", "mean"] <- log(sum((y - data_pred)^2) / df) + #
        len_par * log(len_data)                                                #
    }                                                                          #
                                                                               #
    output_list <- list()                                                      #
    output_list[[1]] <- models[[1]]                                            #
    output_list[[2]] <- model_par_stat                                         #
  }                                                                            #
                                                                               #
  return(output_list)                                                          #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
### Given a formula 'form', a data dependent vector 'y', a data feature 'x', ###
### their respective uncertainties 'unc_y' and 'unc_y', and a set of control ###
# parameters, this function returns the mean of N models obtained with optim() #
## random draws from the data, and the 68% and 95% confidence in boundaries. ###
################################################################################
bootstrap_optim <- function(form, y, x, unc_y = 0, unc_x = 0, fit_N = 100,     #
                            fits_per_draw_N = 1, forceN = F, max_j = 100,      #
                            rnd_method = 'norm', manual_ini = F, par_ini = NA, #
                            par_ini_mean = NA, par_ini_width = NA,             #
                            optim_ctrl = list(maxit = 750, abstol = 1e-8,      #
                                              reltol = 1e-5),                  #
                            use_all = F, minN_corr = .85, maxN_res = .15,      #
                            CI_1sd = T, CI_2sd = T, verbose = F){              #
                                                                               #
  # Test 'form' class                                                          #
  if(class(form) != "formula"){                                                #
    print(paste0("ERROR: 'form' must be of class 'formula'. Instead, it is of",#
                 " class ", class(form), ". Returning NULL."))                 #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  form_str <- paste0(paste0(form)[2]," ", paste0(form)[1]," ", paste0(form)[3])#
  form_vars <- all.vars(form)                                                  #
                                                                               #
  cmd <- tail(as.character(form), 1)                                           #
  exp <- parse(text = cmd)                                                     #
  model_call <- function(...) eval(exp, list(...))                             #
                                                                               #
  # Test name of dependent variable                                            #
  if(form_vars[1] != "Y"){                                                     #
    print(paste0("ERROR: 'form' must include a dependent variable named 'Y'.", #
                 " Returning NULL."))                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test name of independent variable                                          #
  if(length(which(form_vars == "X")) != 1){                                    #
    print(paste0("ERROR: 'form' must include one independent variable named ", #
                 "'X'. Returning NULL."))                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  par_list <- form_vars[-c(which(form_vars == "X"), which(form_vars == "Y"))]  #
  len_par <- length(par_list)                                                  #
                                                                               #
  model_par_list <- c(par_list, "data_corr", "new_corr", "data_res(%)",        #
                      "new_res(%)", "data_chisq", "new_chisq", "data_p-val",   #
                      "new_p-val", "data_BIC", "new_BIC")                      #
  model_par_len <- length(model_par_list)                                      #
                                                                               #
  # Test length of y = x                                                       #
  len_data <- length(y)                                                        #
  len_x <- length(x)                                                           #
                                                                               #
  if(len_data != len_x){                                                       #
    print(paste0("ERROR: length of 'x' must be equal to the length of 'y'. ",  #
                 "'y' has length ", len_data, " and 'x' has length ", len_x,   #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test length of unc_y = y or 1                                              #
  len_unc_y <- length(unc_y)                                                   #
                                                                               #
  if(len_data != len_unc_y & len_unc_y != 1){                                  #
    print(paste0("ERROR: length of 'unc_y' must be either 1 or equal to the ", #
                 "length of 'y'. 'y' has length ", len_data, " and 'unc_y' ",  #
                 "has length ", len_unc_y, ". Returning NULL."))               #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(len_unc_y == 1){                                                          #
    unc_y <- rep(unc_y, len_data)                                              #
  }                                                                            #
                                                                               #
  # Test length of unc_x = x or 1                                              #
  len_unc_x <- length(unc_x)                                                   #
                                                                               #
  if(len_data != len_unc_x & len_unc_x != 1){                                  #
    print(paste0("ERROR: length of 'unc_x' must be either 1 or equal to the ", #
                 "length of 'x'. 'x' has length ", len_data, " and 'unc_x' ",  #
                 "has length ", len_unc_x, ". Returning NULL."))               #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(len_unc_x == 1){                                                          #
    unc_x <- rep(unc_x, len_data)                                              #
  }                                                                            #
                                                                               #
  # Filtering NAs from 'y' and 'x'                                             #
  x_NAs <- which(is.na(x), arr.ind = T)                                        #
  x_NAs <- c(x_NAs, which(is.infinite(x), arr.ind = T))                        #
  x_NAs <- c(x_NAs, which(is.null(x), arr.ind = T))                            #
                                                                               #
  y_NAs <- which(is.na(y), arr.ind = T)                                        #
  y_NAs <- c(y_NAs, which(is.infinite(y), arr.ind = T))                        #
  y_NAs <- c(y_NAs, which(is.null(y), arr.ind = T))                            #
                                                                               #
  inds_to_rm <- union(x_NAs, y_NAs)                                            #
  if(length(inds_to_rm) != 0){                                                 #
    x <- x[-inds_to_rm]                                                        #
    y <- y[-inds_to_rm]                                                        #
    unc_x <- unc_x[-inds_to_rm]                                                #
    unc_y <- unc_y[-inds_to_rm]                                                #
  }                                                                            #
                                                                               #
  len_data <- length(y)                                                        #
                                                                               #
  # Replacing NAs from 'unc_y' and 'unc_x'                                     #
  x_NAs <- which(is.na(unc_x), arr.ind = T)                                    #
  x_NAs <- c(x_NAs, which(is.infinite(unc_x), arr.ind = T))                    #
  x_NAs <- c(x_NAs, which(is.null(unc_x), arr.ind = T))                        #
  unc_x[x_NAs] <- 0                                                            #
                                                                               #
  y_NAs <- which(is.na(unc_y), arr.ind = T)                                    #
  y_NAs <- c(y_NAs, which(is.infinite(unc_y), arr.ind = T))                    #
  y_NAs <- c(y_NAs, which(is.null(unc_y), arr.ind = T))                        #
  unc_y[y_NAs] <- 0                                                            #
                                                                               #
  # Test degrees of freedom                                                    #
  df <- len_data - len_par                                                     #
                                                                               #
  if(df < 0){                                                                  #
    print(paste0("ERROR: negative degrees of freedom. Increase number of data",#
                 " points or choose a model with less parameters. Returning N",#
                 "ULL."))                                                      #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Building objective function                                                #
  obj_func <- function(in_params, xx, yy, in_form){                            #
                                                                               #
    in_form_vars <- all.vars(in_form)                                          #
    in_par_list <- in_form_vars[-c(which(in_form_vars == "X"),                 #
                                   which(in_form_vars == "Y"))]                #
    in_len_par <- length(in_par_list)                                          #
                                                                               #
    in_cmd <- tail(as.character(in_form), 1)                                   #
    in_exp <- parse(text = in_cmd)                                             #
    in_model_call <- function(...) eval(in_exp, list(...))                     #
                                                                               #
    # Assign parameters                                                        #
    params_ls <- list()                                                        #
    for(p in 1:in_len_par){                                                    #
      params_ls[[p]] <- assign(in_par_list[p], in_params[p])                   #
    }                                                                          #
                                                                               #
    # Calculate predicted values                                               #
    Y_pred <- in_model_call(X = xx, params_ls)                                 #
                                                                               #
    res <- 0                                                                   #
    for(i in 1:length(yy)){                                                    #
      if(yy[i] != 0){                                                          #
        res <- res + ((yy[i] - Y_pred[i]) / yy[i])^2                           #
      }else{                                                                   #
        res <- res + (yy[i] - Y_pred[i])^2                                     #
      }                                                                        #
    }                                                                          #
    return(res)                                                                #
  }                                                                            #
                                                                               #
  # Test 'fit_N'                                                               #
  if(length(fit_N) != 1){                                                      #
    print(paste0("WARNING: 'fit_N' is length ", length(fit_N), ". Coercing ",  #
                 "use of first element."))                                     #
    fit_N <- fit_N[1]                                                          #
  }                                                                            #
  if((fit_N - floor(fit_N)) != 0 | fit_N == 0){                                #
    print("ERROR: 'fit_N' must be a positive integer. Returning NULL.")        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'max_j'                                                               #
  if(length(max_j) != 1){                                                      #
    print(paste0("WARNING: 'max_j' is length ", length(max_j), ". Coercing ",  #
                 "use of first element."))                                     #
    max_j <- max_j[1]                                                          #
  }                                                                            #
  if((max_j - floor(max_j)) != 0 | max_j == 0){                                #
    print("ERROR: 'max_j' must be a positive integer. Returning NULL.")        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test 'fits_per_draw_N'                                                     #
  if(length(fits_per_draw_N) != 1){                                            #
    print(paste0("WARNING: 'fits_per_draw_N' is length ",                      #
                 length(fits_per_draw_N), ". Coercing use of first element.")) #
    fits_per_draw_N <- fits_per_draw_N[1]                                      #
  }                                                                            #
  if((fits_per_draw_N - floor(fits_per_draw_N)) != 0 | fits_per_draw_N == 0){  #
    print(paste0("ERROR: 'fits_per_draw_N' must be a positive integer. Return",#
                 "ing NULL."))                                                 #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  while(manual_ini & fits_per_draw_N != 1){                                    #
    print(paste0("WARNING: with manual initialization of parameters ",         #
                 "'fits_per_draw_N' must be set to 1. Otherwise you'd be ",    #
                 "multiple equal fits which would hurt the statistics."))      #
    print("Change or confirm 'manual_ini': ")                                  #
    manual_ini <- edit(manual_ini)                                             #
                                                                               #
    print("Change or confirm 'fits_per_draw_N': ")                             #
    fits_per_draw_N <- edit(fits_per_draw_N)                                   #
                                                                               #
    cancel <- FALSE                                                            #
    print("If you'd rather cancel the execution enter 'TRUE':")                #
    cancel <- edit(cancel)                                                     #
    if(cancel){                                                                #
      return(NULL)                                                             #
    }                                                                          #
  }                                                                            #
  draw_N <- fit_N %/% fits_per_draw_N                                          #
                                                                               #
  # Test 'rnd_method'                                                          #
  valid_rnd <- c('unif', 'norm', 'pois', 'exp')                                #
  if(sum(rnd_method == valid_rnd) != 1){                                       #
    print(paste0("ERROR: invalid choice of distribution from which to draw ",  #
                 "data samples. Valid options are 'unif', 'norm', 'pois' and", #
                 " 'exp'. Returning NULL."))                                   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  # Test parameter initialization                                              #
  if(manual_ini){                                                              #
    if(length(par_ini) > len_par){                                             #
      print(paste0("WARNING: length of 'par_ini' is larger than the number of",#
                   " model parameters. Ignoring those in excess."))            #
      par_ini <- par_ini[1:len_par]                                            #
    }                                                                          #
    if(length(par_ini) < len_par){                                             #
      print(paste0("ERROR: length of 'par_ini' is smaller than the number of ",#
                   "model parameters. Returning NULL."))                       #
      return(NULL)                                                             #
    }                                                                          #
  }else{                                                                       #
    if(length(par_ini_mean) > len_par){                                        #
      print(paste0("WARNING: length of 'par_ini_mean' is larger than the numb",#
                   "er of model parameters. Ignoring those in excess."))       #
      par_ini_mean <- par_ini_mean[1:len_par]                                  #
    }                                                                          #
    if(length(par_ini_mean) < len_par){                                        #
      print(paste0("ERROR: length of 'par_ini_mean' is smaller than the ",     #
                   "number of model parameters. Returning NULL."))             #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    if(prod(is.na(par_ini_width))){                                            #
      if(fits_per_draw_N != 1){                                                #
        print(paste0("WARNING: 'par_ini_width' was not provided while the ",   #
                     "number of fits per sample requested is greater than 1. ",#
                     "Coercing width of parameter initialization to ",         #
                     "sqrt(abs('par_ini_mean'))."))                            #
        par_ini_width <- sqrt(abs(par_ini_mean))                               #
      }                                                                        #
    }else{                                                                     #
      if(length(par_ini_width) > len_par){                                     #
        print(paste0("WARNING: length of 'par_ini_width' is larger than the",  #
                     " number of model parameters. Ignoring those in excess."))#
        par_ini_width <- par_ini_width[1:len_par]                              #
      }                                                                        #
      if(length(par_ini_width) < len_par){                                     #
        print(paste0("WARNING: length of 'par_ini_width' is smaller than the", #
                     " number of model parameters. Remaining values will be ", #
                     "set to sqrt(abs('par_ini_mean'))."))                     #
        par_ini_width <- c(par_ini_width,                                      #
                           sqrt(abs(par_ini_mean[-(1:len_par)])))              #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Check if user is sure to force all fits to be successful                   #
  if(forceN){                                                                  #
    print(paste0("WARNING: forcing all samples fit to succeed may result in a",#
                 " longer execution time, it may also lead to an infinite ",   #
                 "loop. If 'manual_ini' was activated there is also a risk of",#
                 " error in execution to prevent such an infinity loop, due ", #
                 "to not being possible to vary the initial values."))         #
  }                                                                            #
                                                                               #
  if(!use_all){                                                                #
    # Test min_corr                                                            #
    if(length(minN_corr) != 1){                                                #
      print(paste0("WARNING: 'minN_corr' has length greater than 1. Using ",   #
                   "only the first value."))                                   #
      minN_corr <- minN_corr[1]                                                #
    }                                                                          #
    if(minN_corr <= -1 | minN_corr >= 1){                                      #
      print(paste0("WARNING: 'minN_corr' must be greater than -1 and smaller ",#
                   "than 1. Coercing to 0.85."))                               #
      minN_corr <- 0.85                                                        #
    }                                                                          #
                                                                               #
    # Test max_res                                                             #
    if(length(maxN_res) != 1){                                                 #
      print(paste0("WARNING: 'maxN_res' has length greater than 1. Using only",#
                   " the first value."))                                       #
      maxN_res <- maxN_res[1]                                                  #
    }                                                                          #
    if(maxN_res <= 0 | maxN_res >= 1){                                         #
      print(paste0("WARNING: 'maxnN_res' must be greater than 0 and smaller ", #
                   "than 1. Coercing to 0.15."))                               #
      maxN_res <- 0.15                                                         #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Initialize list of models                                                  #
  models <- list()                                                             #
  models_par <- array(NA, dim = c(len_par, fit_N),                             #
                      dimnames = list(par_list, NULL))                         #
  model_par_stat <- array(NA, dim = c(model_par_len, 5),                       #
                          dimnames = list(model_par_list,                      #
                                          c("mean", "ci_68+", "ci_68-",        #
                                            "ci_95+", "ci_95-")))              #
                                                                               #
  sampled_Y <- array(NA, dim = c(len_data, draw_N))                            #
  sampled_X <- array(NA, dim = c(len_data, draw_N))                            #
  for(d in 1:len_data){                                                        #
    if(draw_N != 1){                                                           #
      switch(rnd_method,                                                       #
             'unif' = c(sampled_Y[d,] <- runif(draw_N, y[d] - unc_y[d],        #
                                               y[d] + unc_y[d]),               #
                        sampled_X[d,] <- runif(draw_N, x[d] - unc_x[d],        #
                                               x[d] + unc_x[d])),              #
             'norm' = c(sampled_Y[d,] <- rnorm(draw_N, y[d], unc_y[d]),        #
                        sampled_X[d,] <- rnorm(draw_N, x[d], unc_x[d])),       #
             'pois' = c(sampled_Y[d,] <- rpois(draw_N,                         #
                                               rnorm(1, y[d], unc_y[d])),      #
                        sampled_X[d,] <- rpois(draw_N,                         #
                                               rnorm(1, x[d],unc_x[d]))),      #
             'exp' = c(sampled_Y[d,] <- rexp(draw_N,                           #
                                             1 / rnorm(1,y[d],unc_y[d])),      #
                       sampled_X[d,] <- rexp(draw_N,                           #
                                             1 / rnorm(1,x[d],unc_x[d])))      #
      )                                                                        #
    }else{                                                                     #
      sampled_Y[d,] <- y[d]                                                    #
      sampled_X[d,] <- x[d]                                                    #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Initialize loop of sampling and fitting                                    #
  for(nd in 1:draw_N){                                                         #
                                                                               #
    # Draw data sample to be fitted                                            #
    ys <- sampled_Y[,nd]                                                       #
    xs <- sampled_X[,nd]                                                       #
                                                                               #
    for(nf in 1:fits_per_draw_N){                                              #
                                                                               #
      ns <- fits_per_draw_N * (nd - 1) + nf                                    #
                                                                               #
      # Fitting model                                                          #
      print(paste0("Fitting ", form_str, " model ", ns, "/", fit_N, ":"))      #
                                                                               #
      EorW_flag <- TRUE                                                        #
      j <- 1                                                                   #
                                                                               #
      keep_going <- FALSE                                                      #
      if(forceN){                                                              #
        keep_going <- TRUE                                                     #
      }                                                                        #
                                                                               #
      while((EorW_flag && j < max_j) | keep_going){                            #
        EorW_flag <- FALSE                                                     #
                                                                               #
        # Initialize model parameters                                          #
        if(!manual_ini){                                                       #
          par_ini <- NULL                                                      #
          for(p in 1:len_par){                                                 #
            par_ini <- c(par_ini,                                              #
                         runif(1, par_ini_mean[p] - par_ini_width[p] / 2,      #
                               par_ini_mean[p] + par_ini_width[p] /2))         #
          }                                                                    #
        }                                                                      #
                                                                               #
        prelim <- tryCatch(                                                    #
          optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,           #
                control = optim_ctrl),                                         #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
                                                                               #
        if(is.atomic(prelim)){                                                 #
          EorW_flag <- TRUE                                                    #
          j <- j + 1                                                           #
          if(verbose){                                                         #
            cat(".")                                                           #
          }                                                                    #
        }else{                                                                 #
          if(sum(prelim$convergence == c(1, 10)) != 0){                        #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }else{                                                               #
            keep_going <- FALSE                                                #
          }                                                                    #
        }                                                                      #
      }                                                                        #
                                                                               #
      if(!EorW_flag){                                                          #
        models[[ns]] <- tryCatch(                                              #
          optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,        #
                control = optim_ctrl, method = "BFGS"),                        #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
      }else{                                                                   #
        models[[ns]] <- tryCatch(                                              #
          optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,           #
                control = optim_ctrl, method = "BFGS"),                        #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
      }                                                                        #
                                                                               #
      if(!is.atomic(models[[ns]])){                                            #
        cat("!!Success!!\n")                                                   #
        print(paste0("--> ", form_str, " model ", ns, "/", fit_N," has been ", #
                     "fitted..."))                                             #
      }else{                                                                   #
        cat("!!Failure!!\n")                                                   #
        print(paste0("--> After ", max_j, " attempts, ", form_str, " model ",  #
                     ns, "/", fit_N, " was not fitted..."))                    #
      }                                                                        #
                                                                               #
      # Saving model parameters                                                #
      if(!is.atomic(models[[ns]])){                                            #
        if(models[[ns]]$convergence == 0){                                     #
                                                                               #
          pars_ls <- list()                                                    #
          for(p in 1:len_par){                                                 #
            pars_ls[[p]] <- assign(par_list[p], models[[ns]]$par[p])           #
          }                                                                    #
                                                                               #
          ys_tilda <- model_call(X = xs, pars_ls)                              #
          temp_cor <- cor(ys, ys_tilda, use ="pairwise.complete.obs")          #
          temp_res <- median(abs((ys - ys_tilda) / ys))                        #
                                                                               #
          if(is.na(temp_cor)){                                                 #
            temp_cor <- 0                                                      #
          }                                                                    #
          if(is.na(temp_res)){                                                 #
            temp_res <- 1                                                      #
          }                                                                    #
                                                                               #
          if(!use_all & temp_cor >= minN_corr & temp_res <= maxN_res){         #
            models_par[,ns] <- models[[ns]]$par                                #
          }else{                                                               #
            if(use_all){                                                       #
              models_par[,ns] <- models[[ns]]$par                              #
            }                                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
  }                                                                            #
                                                                               #
  # Estimate average separation between a pair of 'x'                          #
  delta_x <- max(x + unc_x) - min(x - unc_x)                                   #
  avg_Sep <- delta_x / len_data                                                #
                                                                               #
  # Create a new 'x' sequence                                                  #
  new_avg_Sep <- avg_Sep / 100                                                 #
  new_x <- seq(min(x - 3 * unc_x), max(x + 3 * unc_x), by = new_avg_Sep)       #
  len_new_data <- length(new_x)                                                #
  new_df <- len_new_data - len_par                                             #
                                                                               #
  models_pred <- array(NA, dim = c(fit_N, len_new_data))                       #
                                                                               #
  # Make predictions using every successfully fitted model                     #
  if(fit_N > 1){                                                               #
    for(ns in 1:fit_N){                                                        #
      if(!is.atomic(models[[ns]])){                                            #
        if(models[[ns]]$convergence == 0){                                     #
                                                                               #
          model_par_list <- list()                                             #
          for(p in 1:len_par){                                                 #
            model_par_list[[p]] <- assign(par_list[p], models[[ns]]$par[p])    #
          }                                                                    #
          models_pred[ns,] <- model_call(X = new_x, model_par_list)            #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    # Estimate median and sd of model parameters to initialize new models      #
    avg_model_par <- apply(models_par, 1, median, na.rm = T)                   #
    avg_model_par[which(is.na(avg_model_par), arr.ind = T)] <- 0               #
    sd_model_par <- apply(models_par, 1, sd, na.rm = T)                        #
    sd_model_par[which(is.na(sd_model_par), arr.ind = T)] <- 0                 #
                                                                               #
    # Estimate median and sd of predictions                                    #
    avg_models_pred <- apply(models_pred, 2, median, na.rm = T)                #
    sd_models_pred <- apply(models_pred, 2, sd)                                #
                                                                               #
    # Fit model to median of predictions to get median model                   #
    EorW_flag <- TRUE                                                          #
    j <- 1                                                                     #
    ys <- avg_models_pred                                                      #
    xs <- new_x                                                                #
                                                                               #
    print(paste0("Fitting median ", form_str, " model:"))                      #
                                                                               #
    while(EorW_flag && j < max_j){                                             #
      EorW_flag <- FALSE                                                       #
                                                                               #
      # Initialize model parameters                                            #
      par_ini <- NULL                                                          #
      for(p in 1:len_par){                                                     #
        par_ini <- c(par_ini, rnorm(1, avg_model_par[p], sd_model_par[p]))     #
      }                                                                        #
                                                                               #
      # Fit model                                                              #
      prelim <- tryCatch(                                                      #
        optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,             #
              control = optim_ctrl),                                           #
        error = function(e){                                                   #
          return(TRUE)                                                         #
        },                                                                     #
        warning = function(w){                                                 #
          return(TRUE)                                                         #
        }                                                                      #
      )                                                                        #
                                                                               #
      if(is.atomic(prelim)){                                                   #
        EorW_flag <- TRUE                                                      #
        j <- j + 1                                                             #
        if(verbose){                                                           #
          cat(".")                                                             #
        }                                                                      #
      }else{                                                                   #
        if(sum(prelim$convergence == c(1, 10)) != 0){                          #
          EorW_flag <- TRUE                                                    #
          j <- j + 1                                                           #
          if(verbose){                                                         #
            cat(".")                                                           #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    j <- 1                                                                     #
                                                                               #
    if(!EorW_flag){                                                            #
      EorW_flag <- TRUE                                                        #
                                                                               #
      while(EorW_flag && j < max_j){                                           #
        EorW_flag <- FALSE                                                     #
                                                                               #
        avg_model <- tryCatch(                                                 #
          optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,        #
                control = optim_ctrl, method = "BFGS"),                        #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
                                                                               #
        if(is.atomic(avg_model)){                                              #
          EorW_flag <- TRUE                                                    #
          j <- j + 1                                                           #
          if(verbose){                                                         #
            cat(".")                                                           #
          }                                                                    #
        }else{                                                                 #
          if(sum(avg_model$convergence == c(1, 10)) != 0){                     #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }else{                                                                     #
      while(EorW_flag && j < max_j){                                           #
        EorW_flag <- FALSE                                                     #
                                                                               #
        avg_model <- tryCatch(                                                 #
          optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,           #
                control = optim_ctrl, method = "BFGS"),                        #
          error = function(e){                                                 #
            return(TRUE)                                                       #
          },                                                                   #
          warning = function(w){                                               #
            return(TRUE)                                                       #
          }                                                                    #
        )                                                                      #
                                                                               #
        if(is.atomic(avg_model)){                                              #
          EorW_flag <- TRUE                                                    #
          j <- j + 1                                                           #
          if(verbose){                                                         #
            cat(".")                                                           #
          }                                                                    #
        }else{                                                                 #
          if(sum(avg_model$convergence == c(1, 10)) != 0){                     #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    if(!EorW_flag){                                                            #
      cat("!!Success!!\n")                                                     #
      print(paste0("--> Median ", form_str, " model has been fitted..."))      #
    }else{                                                                     #
      cat("!!Failure!!\n")                                                     #
      print(paste0("--> After ", max_j, " attempts, median ", form_str,        #
                   " model was not fitted. Skipping Conf. Int. fits."))        #
    }                                                                          #
                                                                               #
    # Saving model parameters                                                  #
    if(!is.atomic(avg_model)){                                                 #
      if(avg_model$convergence == 0){                                          #
                                                                               #
        model_par_list <- list()                                               #
        for(p in 1:len_par){                                                   #
          model_par_list[[p]] <- assign(par_list[p], avg_model$par[p])         #
        }                                                                      #
                                                                               #
        avg_model_pred <- model_call(X = new_x, model_par_list)                #
        avg_data_pred <- model_call(X = x, model_par_list)                     #
                                                                               #
        new_corr <- cor(ys, avg_model_pred, use = "pairwise.complete.obs")     #
        data_corr <- cor(y, avg_data_pred, use = "pairwise.complete.obs")      #
                                                                               #
        new_res <- median((abs((ys - avg_model_pred) / ys) * 100), na.rm = T)  #
        data_res <- median((abs((y - avg_data_pred) / y) * 100), na.rm = T)    #
                                                                               #
        new_chisq <- sum((ys - avg_model_pred)^2 / ys, na.rm=T)                #
        data_chisq <- sum((y - avg_data_pred)^2 / y, na.rm=T)                  #
                                                                               #
        new_BIC <- log(sum((ys - avg_model_pred)^2) / new_df) + len_par *      #
          log(len_new_data)                                                    #
        data_BIC <- log(sum((y - avg_data_pred)^2) / df) + len_par *           #
          log(len_data)                                                        #
                                                                               #
        model_par_stat[1:len_par, "mean"] <- avg_model$par                     #
        model_par_stat["data_corr", "mean"] <- data_corr                       #
        model_par_stat["new_corr", "mean"] <- new_corr                         #
        model_par_stat["data_res(%)", "mean"] <- data_res                      #
        model_par_stat["new_res(%)", "mean"] <- new_res                        #
        model_par_stat["data_chisq", "mean"] <- data_chisq                     #
        model_par_stat["new_chisq", "mean"] <- new_chisq                       #
        model_par_stat["data_p-val", "mean"] <- 1 - pchisq(data_chisq, df)     #
        model_par_stat["new_p-val", "mean"] <- 1 - pchisq(new_chisq, new_df)   #
        model_par_stat["data_BIC", "mean"] <- data_BIC                         #
        model_par_stat["new_BIC", "mean"] <- new_BIC                           #
      }                                                                        #
    }                                                                          #
                                                                               #
    # Fit model to median +- sd of predictions to get CI 68% of avg_model      #
    if(!is.atomic(avg_model)){                                                 #
      if(CI_1sd){                                                              #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        ys <- avg_models_pred + sd_models_pred                                 #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          par_ini <- NULL                                                      #
          for(p in 1:len_par){                                                 #
            par_ini <- c(par_ini, rnorm(1, avg_model_par[p], sd_model_par[p])) #
          }                                                                    #
                                                                               #
          # Fit upper bound model                                              #
          prelim <- tryCatch(                                                  #
            optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,         #
                  control = optim_ctrl),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
                                                                               #
          if(is.atomic(prelim)){                                               #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }else{                                                               #
            if(sum(prelim$convergence == c(1, 10)) != 0){                      #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
                                                                               #
        j <- 1                                                                 #
                                                                               #
        if(!EorW_flag){                                                        #
          EorW_flag <- TRUE                                                    #
                                                                               #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd1u_model <- tryCatch(                                            #
              optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,    #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd1u_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd1u_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }else{                                                                 #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd1u_model <- tryCatch(                                            #
              optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,       #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd1u_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd1u_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Upper bound of 68% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, upper bound of 68% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        # Saving model parameters                                              #
        if(!is.atomic(sd1u_model)){                                            #
          if(sd1u_model$convergence == 0){                                     #
                                                                               #
            model_par_list <- list()                                           #
            for(p in 1:len_par){                                               #
              model_par_list[[p]] <- assign(par_list[p], sd1u_model$par[p])    #
            }                                                                  #
                                                                               #
            model_pred <- model_call(X = new_x, model_par_list)                #
            data_pred <- model_call(X = x, model_par_list)                     #
                                                                               #
            new_corr <- cor(ys, model_pred, use = "pairwise.complete.obs")     #
            data_corr <- cor(y, data_pred, use = "pairwise.complete.obs")      #
                                                                               #
            new_res <- median((abs((ys - model_pred) / ys) * 100), na.rm = T)  #
            data_res <- median((abs((y - data_pred) / y) * 100), na.rm = T)    #
                                                                               #
            new_chisq <- sum((ys - model_pred)^2 / ys, na.rm=T)                #
            data_chisq <- sum((y - data_pred)^2 / y, na.rm=T)                  #
                                                                               #
            new_BIC <- log(sum((ys - model_pred)^2) / new_df) + len_par *      #
              log(len_new_data)                                                #
            data_BIC <- log(sum((y - data_pred)^2) / df) + len_par *           #
              log(len_data)                                                    #
                                                                               #
            model_par_stat[1:len_par, "ci_68+"] <- sd1u_model$par              #
            model_par_stat["data_corr", "ci_68+"] <- data_corr                 #
            model_par_stat["new_corr", "ci_68+"] <- new_corr                   #
            model_par_stat["data_res(%)", "ci_68+"] <- data_res                #
            model_par_stat["new_res(%)", "ci_68+"] <- new_res                  #
            model_par_stat["data_chisq", "ci_68+"] <- data_chisq               #
            model_par_stat["new_chisq", "ci_68+"] <- new_chisq                 #
            model_par_stat["data_p-val", "ci_68+"] <- 1 - pchisq(data_chisq,   #
                                                                 df)           #
            model_par_stat["new_p-val", "ci_68+"] <- 1 - pchisq(new_chisq,     #
                                                                new_df)        #
            model_par_stat["data_BIC", "ci_68+"] <- data_BIC                   #
            model_par_stat["new_BIC", "ci_68+"] <- new_BIC                     #
          }                                                                    #
        }                                                                      #
                                                                               #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        ys <- avg_models_pred - sd_models_pred                                 #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          par_ini <- NULL                                                      #
          for(p in 1:len_par){                                                 #
            par_ini <- c(par_ini, rnorm(1, avg_model_par[p], sd_model_par[p])) #
          }                                                                    #
                                                                               #
          # Fit lower bound model                                              #
          prelim <- tryCatch(                                                  #
            optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,         #
                  control = optim_ctrl),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
                                                                               #
          if(is.atomic(prelim)){                                               #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }else{                                                               #
            if(sum(prelim$convergence == c(1, 10)) != 0){                      #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
                                                                               #
        j <- 1                                                                 #
                                                                               #
        if(!EorW_flag){                                                        #
          EorW_flag <- TRUE                                                    #
                                                                               #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd1l_model <- tryCatch(                                            #
              optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,    #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd1l_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd1l_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }else{                                                                 #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd1l_model <- tryCatch(                                            #
              optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,       #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd1l_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd1l_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Lower bound of 68% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, lower bound of 68% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        # Saving model parameters                                              #
        if(!is.atomic(sd1l_model)){                                            #
          if(sd1l_model$convergence == 0){                                     #
                                                                               #
            model_par_list <- list()                                           #
            for(p in 1:len_par){                                               #
              model_par_list[[p]] <- assign(par_list[p], sd1l_model$par[p])    #
            }                                                                  #
                                                                               #
            model_pred <- model_call(X = new_x, model_par_list)                #
            data_pred <- model_call(X = x, model_par_list)                     #
                                                                               #
            new_corr <- cor(ys, model_pred, use = "pairwise.complete.obs")     #
            data_corr <- cor(y, data_pred, use = "pairwise.complete.obs")      #
                                                                               #
            new_res <- median((abs((ys - model_pred) / ys) * 100), na.rm = T)  #
            data_res <- median((abs((y - data_pred) / y) * 100), na.rm = T)    #
                                                                               #
            new_chisq <- sum((ys - model_pred)^2 / ys, na.rm=T)                #
            data_chisq <- sum((y - data_pred)^2 / y, na.rm=T)                  #
                                                                               #
            new_BIC <- log(sum((ys - model_pred)^2) / new_df) + len_par *      #
              log(len_new_data)                                                #
            data_BIC <- log(sum((y - data_pred)^2) / df) + len_par *           #
              log(len_data)                                                    #
                                                                               #
            model_par_stat[1:len_par, "ci_68-"] <- sd1l_model$par              #
            model_par_stat["data_corr", "ci_68-"] <- data_corr                 #
            model_par_stat["new_corr", "ci_68-"] <- new_corr                   #
            model_par_stat["data_res(%)", "ci_68-"] <- data_res                #
            model_par_stat["new_res(%)", "ci_68-"] <- new_res                  #
            model_par_stat["data_chisq", "ci_68-"] <- data_chisq               #
            model_par_stat["new_chisq", "ci_68-"] <- new_chisq                 #
            model_par_stat["data_p-val", "ci_68-"] <- 1 - pchisq(data_chisq,   #
                                                                 df)           #
            model_par_stat["new_p-val", "ci_68-"] <- 1 - pchisq(new_chisq,     #
                                                                new_df)        #
            model_par_stat["data_BIC", "ci_68-"] <- data_BIC                   #
            model_par_stat["new_BIC", "ci_68-"] <- new_BIC                     #
          }                                                                    #
        }                                                                      #
      }                                                                        #
                                                                               #
      if(CI_2sd){                                                              #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        ys <- avg_models_pred + 2 * sd_models_pred                             #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          par_ini <- NULL                                                      #
          for(p in 1:len_par){                                                 #
            par_ini <- c(par_ini, rnorm(1, avg_model_par[p], sd_model_par[p])) #
          }                                                                    #
                                                                               #
          # Fit upper bound model                                              #
          prelim <- tryCatch(                                                  #
            optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,         #
                  control = optim_ctrl),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
                                                                               #
          if(is.atomic(prelim)){                                               #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }else{                                                               #
            if(sum(prelim$convergence == c(1, 10)) != 0){                      #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
                                                                               #
        j <- 1                                                                 #
                                                                               #
        if(!EorW_flag){                                                        #
          EorW_flag <- TRUE                                                    #
                                                                               #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd2u_model <- tryCatch(                                            #
              optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,    #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd2u_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd2u_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }else{                                                                 #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd2u_model <- tryCatch(                                            #
              optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,       #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd2u_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd2u_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Upper bound of 95% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, upper bound of 95% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        # Saving model parameters                                              #
        if(!is.atomic(sd2u_model)){                                            #
          if(sd2u_model$convergence == 0){                                     #
                                                                               #
            model_par_list <- list()                                           #
            for(p in 1:len_par){                                               #
              model_par_list[[p]] <- assign(par_list[p], sd2u_model$par[p])    #
            }                                                                  #
                                                                               #
            model_pred <- model_call(X = new_x, model_par_list)                #
            data_pred <- model_call(X = x, model_par_list)                     #
                                                                               #
            new_corr <- cor(ys, model_pred, use = "pairwise.complete.obs")     #
            data_corr <- cor(y, data_pred, use = "pairwise.complete.obs")      #
                                                                               #
            new_res <- median((abs((ys - model_pred) / ys) * 100), na.rm = T)  #
            data_res <- median((abs((y - data_pred) / y) * 100), na.rm = T)    #
                                                                               #
            new_chisq <- sum((ys - model_pred)^2 / ys, na.rm=T)                #
            data_chisq <- sum((y - data_pred)^2 / y, na.rm=T)                  #
                                                                               #
            new_BIC <- log(sum((ys - model_pred)^2) / new_df) + len_par *      #
              log(len_new_data)                                                #
            data_BIC <- log(sum((y - data_pred)^2) / df) + len_par *           #
              log(len_data)                                                    #
                                                                               #
            model_par_stat[1:len_par, "ci_95+"] <- sd2u_model$par              #
            model_par_stat["data_corr", "ci_95+"] <- data_corr                 #
            model_par_stat["new_corr", "ci_95+"] <- new_corr                   #
            model_par_stat["data_res(%)", "ci_95+"] <- data_res                #
            model_par_stat["new_res(%)", "ci_95+"] <- new_res                  #
            model_par_stat["data_chisq", "ci_95+"] <- data_chisq               #
            model_par_stat["new_chisq", "ci_95+"] <- new_chisq                 #
            model_par_stat["data_p-val", "ci_95+"] <- 1 - pchisq(data_chisq,   #
                                                                 df)           #
            model_par_stat["new_p-val", "ci_95+"] <- 1 - pchisq(new_chisq,     #
                                                                new_df)        #
            model_par_stat["data_BIC", "ci_95+"] <- data_BIC                   #
            model_par_stat["new_BIC", "ci_95+"] <- new_BIC                     #
          }                                                                    #
        }                                                                      #
                                                                               #
        EorW_flag <- TRUE                                                      #
        j <- 1                                                                 #
        ys <- avg_models_pred - 2 * sd_models_pred                             #
                                                                               #
        while(EorW_flag && j < max_j){                                         #
          EorW_flag <- FALSE                                                   #
                                                                               #
          # Initialize model parameters                                        #
          par_ini <- NULL                                                      #
          for(p in 1:len_par){                                                 #
            par_ini <- c(par_ini, rnorm(1, avg_model_par[p], sd_model_par[p])) #
          }                                                                    #
                                                                               #
          # Fit lower bound model                                              #
          prelim <- tryCatch(                                                  #
            optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,         #
                  control = optim_ctrl),                                       #
            error = function(e){                                               #
              return(TRUE)                                                     #
            },                                                                 #
            warning = function(w){                                             #
              return(TRUE)                                                     #
            }                                                                  #
          )                                                                    #
                                                                               #
          if(is.atomic(prelim)){                                               #
            EorW_flag <- TRUE                                                  #
            j <- j + 1                                                         #
            if(verbose){                                                       #
              cat(".")                                                         #
            }                                                                  #
          }else{                                                               #
            if(sum(prelim$convergence == c(1, 10)) != 0){                      #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
                                                                               #
        j <- 1                                                                 #
                                                                               #
        if(!EorW_flag){                                                        #
          EorW_flag <- TRUE                                                    #
                                                                               #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd2l_model <- tryCatch(                                            #
              optim(prelim$par, obj_func, xx = xs, yy = ys, in_form = form,    #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd2l_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd2l_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }else{                                                                 #
          while(EorW_flag && j < max_j){                                       #
            EorW_flag <- FALSE                                                 #
                                                                               #
            sd2l_model <- tryCatch(                                            #
              optim(par_ini, obj_func, xx = xs, yy = ys, in_form = form,       #
                    control = optim_ctrl, method = "BFGS"),                    #
              error = function(e){                                             #
                return(TRUE)                                                   #
              },                                                               #
              warning = function(w){                                           #
                return(TRUE)                                                   #
              }                                                                #
            )                                                                  #
                                                                               #
            if(is.atomic(sd2l_model)){                                         #
              EorW_flag <- TRUE                                                #
              j <- j + 1                                                       #
              if(verbose){                                                     #
                cat(".")                                                       #
              }                                                                #
            }else{                                                             #
              if(sum(sd2l_model$convergence == c(1, 10)) != 0){                #
                EorW_flag <- TRUE                                              #
                j <- j + 1                                                     #
                if(verbose){                                                   #
                  cat(".")                                                     #
                }                                                              #
              }                                                                #
            }                                                                  #
          }                                                                    #
        }                                                                      #
        if(!EorW_flag){                                                        #
          cat("!!Success!!\n")                                                 #
          print(paste0("--> Lower bound of 95% Confidence Interval for ",      #
                       form_str, " model has been fitted..."))                 #
        }else{                                                                 #
          cat("!!Failure!!\n")                                                 #
          print(paste0("--> After ", max_j, " attempts, lower bound of 95% Co",#
                       "nf. Interval for ", form_str, " model was not fitted.",#
                       ".."))                                                  #
        }                                                                      #
                                                                               #
        # Saving model parameters                                              #
        if(!is.atomic(sd2l_model)){                                            #
          if(sd2l_model$convergence == 0){                                     #
                                                                               #
            model_par_list <- list()                                           #
            for(p in 1:len_par){                                               #
              model_par_list[[p]] <- assign(par_list[p], sd2l_model$par[p])    #
            }                                                                  #
                                                                               #
            model_pred <- model_call(X = new_x, model_par_list)                #
            data_pred <- model_call(X = x, model_par_list)                     #
                                                                               #
            new_corr <- cor(ys, model_pred, use = "pairwise.complete.obs")     #
            data_corr <- cor(y, data_pred, use = "pairwise.complete.obs")      #
                                                                               #
            new_res <- median((abs((ys - model_pred) / ys) * 100), na.rm = T)  #
            data_res <- median((abs((y - data_pred) / y) * 100), na.rm = T)    #
                                                                               #
            new_chisq <- sum((ys - model_pred)^2 / ys, na.rm=T)                #
            data_chisq <- sum((y - data_pred)^2 / y, na.rm=T)                  #
                                                                               #
            new_BIC <- log(sum((ys - model_pred)^2) / new_df) + len_par *      #
              log(len_new_data)                                                #
            data_BIC <- log(sum((y - data_pred)^2) / df) + len_par *           #
              log(len_data)                                                    #
                                                                               #
            model_par_stat[1:len_par, "ci_95-"] <- sd2l_model$par              #
            model_par_stat["data_corr", "ci_95-"] <- data_corr                 #
            model_par_stat["new_corr", "ci_95-"] <- new_corr                   #
            model_par_stat["data_res(%)", "ci_95-"] <- data_res                #
            model_par_stat["new_res(%)", "ci_95-"] <- new_res                  #
            model_par_stat["data_chisq", "ci_95-"] <- data_chisq               #
            model_par_stat["new_chisq", "ci_95-"] <- new_chisq                 #
            model_par_stat["data_p-val", "ci_95-"] <- 1 - pchisq(data_chisq,   #
                                                                 df)           #
            model_par_stat["new_p-val", "ci_95-"] <- 1 - pchisq(new_chisq,     #
                                                                new_df)        #
            model_par_stat["data_BIC", "ci_95-"] <- data_BIC                   #
            model_par_stat["new_BIC", "ci_95-"] <- new_BIC                     #
          }                                                                    #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    output_list <- list()                                                      #
    output_list[[1]] <- avg_model                                              #
    output_list[[2]] <- model_par_stat                                         #
                                                                               #
    if(CI_1sd && !is.atomic(avg_model)){                                       #
      output_list[[length(output_list) + 1]] <- sd1u_model                     #
      output_list[[length(output_list) + 1]] <- sd1l_model                     #
    }                                                                          #
    if(CI_2sd && !is.atomic(avg_model)){                                       #
      output_list[[length(output_list) + 1]] <- sd2u_model                     #
      output_list[[length(output_list) + 1]] <- sd2l_model                     #
    }                                                                          #
  }else{                                                                       #
    if(!is.logical(models[[1]])){                                              #
                                                                               #
      model_par_stat[1:len_par, "mean"] <- models[[1]]$m$getPars()             #
      model_par_list <- list()                                                 #
      for(p in 1:len_par){                                                     #
        model_par_list[[p]] <- assign(par_list[p], models[[1]]$m$getPars()[p]) #
      }                                                                        #
                                                                               #
      data_pred <- model_call(X = x, model_par_list)                           #
      model_par_stat["data_corr", "mean"] <-                                   #
        cor(y, data_pred, use = "pairwise.complete.obs")                       #
      model_par_stat["data_res(%)", "mean"] <-                                 #
        median((abs((y - data_pred) / y) * 100), na.rm = T)                    #
      model_par_stat["data_chisq", "mean"] <-                                  #
        sum((y - data_pred)^2 / y, na.rm = T)                                  #
      model_par_stat["data_p-val", "mean"] <-                                  #
        1 - pchisq(model_par_stat["data_chisq", "mean"], df)                   #
      model_par_stat["data_BIC", "mean"] <- log(sum((y - data_pred)^2) / df) + #
        len_par * log(len_data)                                                #
    }                                                                          #
                                                                               #
    output_list <- list()                                                      #
    output_list[[1]] <- models[[1]]                                            #
    output_list[[2]] <- model_par_stat                                         #
  }                                                                            #
                                                                               #
  return(output_list)                                                          #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
#          Helper function to calculate the normalization factor C_N,          #
#           used in perturbation function for spiral arm definition            #
#                       N -> index for inter-arm ratio                         #
################################################################################
calc_C_N <- function(N){                                                       #
  return(sqrt(pi) * gamma(N + 1) / gamma(N + 0.5))                             #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Perturbation function to calculate xi(R, phi) for general N                  #
# R -> polar radius; phi -> polar angle; m -> number of arms, p -> pitch;      #
# R0 -> reference radius; phi0 -> reference angle;                             #
# w -> perturbation weight factor; N -> index for inter-arm ratio              #
################################################################################
xi_general <- function(R, phi, m, p, R0, phi0, w, N){                          #
  C_N <- calc_C_N(N)                                                           #
  ang <- (m / 2) * (log(R / R0) / tan(p) - (phi - phi0)) + pi / 4              #
                                                                               #
  return((1 - w) + w * C_N * (sin(ang))^(2 * N))                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Perturbation function to calculate xi(R, phi) for N = 1                      #
# R -> polar radius; phi -> polar angle; m -> number of arms, p -> pitch;      #
# R0 -> reference radius; phi0 -> reference angle;                             #
# w -> perturbation weight factor; N -> index for inter-arm ratio              #
################################################################################
xi_special_case <- function(R, phi, m, p, R0, phi0, w){                        #
  ang <- m * (log(R / R0) / tan(p) - (phi - phi0))                             #
                                                                               #
  return(1 + w * sin(ang))                                                     #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Density function for spiral arm structure                                    #
# R -> polar radius; R_s -> scale radius; z -> height; z0 -> scale height;     #
# rho0 -> central density                                                      #
################################################################################
rho_ax <- function(R, R_s, rho0, z = 0, z0 = 1){                               #
  return(rho0 * exp(-R / R_s) * exp(-abs(z) / z0))                             #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Calculate the angle of a line connecting a point to a reference point        #
# This angle is measured in regard to the North (positive y-axis)              #
################################################################################
point_angle <- function(x, y, ref_x, ref_y){                                   #
  return(atan2(x - ref_x, ref_y - y))                                          #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Calculate the euclidian distance between a point and a reference point       #
################################################################################
point_dist <- function(x, y, ref_x, ref_y){                                    #
  sqrt((x - ref_x)^2 + (y - ref_y)^2)                                          #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Normalize an angle to the interval [-pi/2, pi/2]                             #
################################################################################
norm_angle <- function(ang){                                                   #
  # Normalize to [0, 2*pi]                                                     #
  ang <- ang %% (2 * pi)                                                       #
                                                                               #
  # Normalize to [-pi, pi]                                                     #
  if(ang > pi){                                                                #
    ang <- ang - 2 * pi                                                        #
  }                                                                            #
                                                                               #
  # Normalize to [-pi/2, pi/2]                                                 #
  if(ang > pi / 2){                                                            #
    ang <- ang - pi                                                            #
  }                                                                            #
  if(ang <= (-pi / 2)){                                                        #
    ang <- ang + pi                                                            #
  }                                                                            #
  return(ang)                                                                  #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Calculates tangential direction of a spiral (centered at ref_x, ref_y) with  # 
# pitch parameter "b" (b = tan(pitch)), starting angle "start_ang",            # 
# and reference radius "R", at an angle defined by the direction of            #
# a point x,y relative to the center of the spiral.                            #
# Returns both the tangential direction of the spiral at that angle and the    #
# dist from the point to the spiral.                                           #
################################################################################
spiral_direction <- function(x, y, ref_x, ref_y, start_ang, R, b){             #
  # Angle and radial distance of the point                                     #
  ang <- point_angle(x, y, ref_x, ref_y)                                       #
  dist <- point_dist(x, y, ref_x, ref_y)                                       #
                                                                               #
  # Expected radius for this angle on the spiral                               #
  expt_radius <- R * exp(b * (ang - start_ang))                                #
                                                                               #
  # Calculate tangential direction at this point on the spiral                 #
  tan_ang <- ang + atan(1 / b)                                                 #
                                                                               #
  # Normalize the tangent angle to [-pi/2, pi/2]                               #
  tan_ang <- norm_angle(tan_ang)                                               #
                                                                               #
  # Return the difference in distance and the tangent angle                    #
  return(list(dist_diff = abs(dist - expt_radius), direction = tan_ang))       #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'map' and a reference pixel coordinate pair 'cpix'                #
# Returns 'output_map', a 3D array where the 'data' layer holds the distance   #
# of each pixel to the reference pixel, and the 'unc' layer holds the          #
# respective uncertainty.                                                      #
################################################################################
create_dist_map <- function(map, cpix = dim(map) %/% 2 + 1, narm = F){         #
                                                                               #
  dims <- dim(map)                                                             #
  ldims <- length(dims)                                                        #
  lcpix <- length(cpix)                                                        #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims, " dimension. Returning NULL."))#
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lcpix != 2){                                                              #
    print(paste0("ERROR: This function expects 'cpix' to have length 2. ",     #
                 "'cpix' has instead length ", lcpix, ". Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  output_map <- array(NA, dim = c(dims, 2),                                    #
                      dimnames = list(NULL, NULL, c("data", "unc")))           #
  # corrected uncertainty of each pixel coordinate is 1 / sqrt(2) pix          #
  # uncertainty of each distance, via propagation, will be 1                   #
                                                                               #
  pix_unc <- 1 / sqrt(2)                                                       #
                                                                               #
  cx <- cpix[1]                                                                #
  cy <- cpix[2]                                                                #
                                                                               #
  for(x in 1:dims[1]){                                                         #
    for(y in 1:dims[2]){                                                       #
                                                                               #
      if(narm == T & is.na(map[x,y])){                                         #
        next                                                                   #
      }                                                                        #
                                                                               #
      dx <- (cx - x)                                                           #
      dy <- (cy - y)                                                           #
      unc_dif <- 2 * pix_unc                                                   #
                                                                               #
      output_map[x, y, "data"] <- sqrt(dx^2 + dy^2)                            #
      output_map[x, y, "unc"] <- unc_dif / sqrt(output_map[x, y, "data"])      #
    }                                                                          #
  }                                                                            #
  return(output_map)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'map' and a reference pixel coordinate pair 'cpix'                #
# Returns 'output_map', a 3D array where the 'data' layer holds the angle,     #
# in degrees, made between lines, made by each pixel and the reference pixel,  #
# and the North # (Y-axis); the 'unc' layer holds the respective uncertainty.  #
################################################################################
create_ang_map <- function(map, cpix = dim(map) %/% 2 + 1, narm = F){          #
                                                                               #
  dims <- dim(map)                                                             #
  ldims <- length(dims)                                                        #
  lcpix <- length(cpix)                                                        #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims, " dimension. Returning NULL."))#
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(lcpix != 2){                                                              #
    print(paste0("ERROR: This function expects 'cpix' to have length 2. ",     #
                 "'cpix' has instead length ", lcpix, ". Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  output_map <- array(NA, dim = c(dims, 2),                                    #
                      dimnames = list(NULL, NULL, c("data", "unc")))           #
                                                                               #
  pix_unc <- 1 / sqrt(2)                                                       #
                                                                               #
  cx <- cpix[1]                                                                #
  cy <- cpix[2]                                                                #
                                                                               #
  for(x in 1:dims[1]){                                                         #
    for(y in 1:dims[2]){                                                       #
      if(narm == T & is.na(map[x,y])){                                         #
        next                                                                   #
      }                                                                        #
                                                                               #
      dx <- (x - cx)                                                           #
      dy <- (y - cy)                                                           #
      unc_dif <- 2 * pix_unc                                                   #
                                                                               #
      ang <- atan2(dy, dx) * 180 / pi                                          #
      mag <- sqrt(dx^2 + dy^2)                                                 #
      unc_ang <- 90 * unc_dif / (pi * mag)                                     #
                                                                               #
      if(ang < 0){                                                             #
        ang <- 360 + ang                                                       #
      }                                                                        #
                                                                               #
      output_map[x, y, "data"] <- ang                                          #
      output_map[x, y, "unc"] <- unc_ang                                       #
    }                                                                          #
  }                                                                            #
  return(output_map)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'map' and one of two sets of inputs,                              #
# Returns 'ra_map', a 3D array where the 'data' layer holds                    #
# the right ascension of each pixel, and the 'unc' layer holds the respective  #
# uncertainty.                                                                 #
# Input set 1: reference pixel coordinate pair 'cpix'; that pixels RA 'cra' &  #
#              DEC 'cdec' in degrees;                                          #
#              the pixel scale 'pix_scl' in degrees per pixel;                 #
#              and the uncertainties of RA and DEC                             #
#              'unc_ra' and 'unc_dec' in degrees.                              #
#                                                                              #
# Input set 2: the uncertainty of RA 'unc_ra' in degrees;                      #
#              the path to a WCS fits file 'path_to_wcs_file'.                 #
#                                                                              #
# To chose between input sets there is the variable 'method' which accepts     #
# the values "pixel_scale" or "wcs_file".                                      #
################################################################################
create_RA_map <- function(map, cpix = dim(map) %/% 2 + 1,                      #
                          cra, cdec, pix_scl, unc_ra = pix_scl,                #
                          unc_dec = pix_scl,                                   #
                          path_to_wcs_file, method = "pixel_scale"){           #
                                                                               #
  if(length(method) != 1 | (method != "pixel_scale" & method != "wcs_file")){  #
    print(paste0("ERROR: invalid 'method' value. 'method' must be either ",    #
                 " 'pixel_scale' or 'wcs_file'. Returning NULL."))             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dims <- dim(map)                                                             #
                                                                               #
  if(is.null(dims)){                                                           #
    print(paste0("ERROR: 'map' is not an array and this function expects it",  #
                 " to be a 2D array. Returning NULL."))                        #
  }                                                                            #
                                                                               #
  ldims <- length(dims)                                                        #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimension. Returning NULL."))                             #
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(method == "pixel_scale"){                                                 #
    lcpix <- length(cpix)                                                      #
                                                                               #
    if(lcpix != 2){                                                            #
      print(paste0("ERROR: This function expects 'cpix' to have length 2. ",   #
                   "'cpix' has instead length ", lcpix, ". Returning NULL."))  #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    if((length(cra) != 1) | (length(cdec) != 1) |                              #
       (length(unc_ra) != 1) | (length(unc_dec) != 1) |                        #
       (length(pix_scl != 1))){                                                #
      print(paste0("ERROR: This function expects 'cra', 'cdec', 'unc_ra', ",   #
                   "'unc_dec' and 'pix_scl' to be length 1. Returning NULL.")) #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    ra_map <- array(NA, dim = c(dims, 2),                                      #
                    dimnames = list(NULL, NULL, c("data", "unc")))             #
                                                                               #
    unc_mdec <- 0.5 * sqrt(2 * unc_dec^2)^2 * pi / 180 # to rads               #
                                                                               #
    for(x in 1:dims[1]){                                                       #
      for(y in 1:dims[2]){                                                     #
                                                                               #
        if(narm == T & is.na(map[x,y])){                                       #
          next                                                                 #
        }                                                                      #
                                                                               #
        difx <- cpix[1] - x                                                    #
        dify <- y - cpix[2]                                                    #
                                                                               #
        dec <- cdec + dify * pix_scl                                           #
        mdec <- mean(c(dec, cdec)) * pi / 180 # to rads                        #
        cosd <- cos(mdec)                                                      #
        unc_cosd <- sqrt((sin(mdec) * unc_mdec)^2)                             #
                                                                               #
        ra_map[x, y, "data"] <- cra + difx * pix_scl / cosd                    #
        ra_map[x, y, "unc"] <- sqrt(2 * unc_ra^2 +                             #
                                      (difx * pix_scl / cosd^2 * unc_cosd)^2)  #
      }                                                                        #
    }                                                                          #
  }else{                                                                       #
    if(length(path_to_wcs_file) != 1){                                         #
      print(paste0("ERROR: 'path_to_wcs_file' must be a single string. ",      #
                   "Returning NULL."))                                         #
      return(NULL)                                                             #
    }                                                                          #
    if(!file.exists(path_to_wcs_file)){                                        #
      print(paste0("ERROR: 'path_to_wcs_file' does not exist. ",               #
                   "Returning NULL."))                                         #
      return(NULL)                                                             #
    }                                                                          #
    if(length(unc_ra) != 1){                                                   #
      print(paste0("ERROR: This function expects 'unc_ra' to be a length 1 ",  #
                   "constants. Returning NULL."))                              #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    ra_map <- array(NA, dim = c(dims, 2),                                      #
                    dimnames = list(NULL, NULL, c("data", "unc")))             #
                                                                               #
    for(x in 1:dims[1]){                                                       #
      t_coords <- c(rep(x, times = dims[2]), 1:dims[2])                        #
      ra_map[x,,"data"] <- convert_pixel_to_wcs(array(t_coords,                #
                                                      dim = c(dims[2], 2)),    #
                                                path_to_wcs_file)[,"ra"]       #
      ra_map[x,,"unc"] <- unc_ra                                               #
    }                                                                          #
    rm(t_coords)                                                               #
  }                                                                            #
  return(ra_map)                                                               #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'map' and one of two sets of inputs,                              #
# Returns 'dec_map', a 3D array where the 'data' layer holds                   #
# the declination of each pixel, and the 'unc' layer holds the respective      #
# uncertainty.                                                                 #
# Input set 1: reference pixel coordinate pair 'cpix'; that pixels DEC 'cdec'  #
#              in degrees; the pixel scale 'pix_scl' in degrees per pixel;     #
#              and the uncertainties of DEC 'unc_dec' in degrees.              #
#                                                                              #
# Input set 2: the uncertainty of DEC 'unc_dec' in degrees;                    #
#              the path to a WCS fits file 'path_to_wcs_file'.                 #
#                                                                              #
# To chose between input sets there is the variable 'method' which accepts     #
# the values "pixel_scale" or "wcs_file".                                      #
################################################################################
create_DEC_map <- function(map, cpix = dim(map) %/% 2 + 1, cdec, pix_scl,      #
                           unc_dec = pix_scl, path_to_wcs_file,                #
                           method = "pixel_scale"){                            #
                                                                               #
  if(length(method) != 1 | (method != "pixel_scale" & method != "wcs_file")){  #
    print(paste0("ERROR: invalid 'method' value. 'method' must be either ",    #
                 " 'pixel_scale' or 'wcs_file'. Returning NULL."))             #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  dims <- dim(map)                                                             #
                                                                               #
  if(is.null(dims)){                                                           #
    print(paste0("ERROR: 'map' is not an array and this function expects it",  #
                 " to be a 2D array. Returning NULL."))                        #
  }                                                                            #
                                                                               #
  ldims <- length(dims)                                                        #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimension. Returning NULL."))                             #
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'map' to be a 2D array. ",    #
                   "'map' has instead ", ldims,                                #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(method == "pixel_scale"){                                                 #
    lcpix <- length(cpix)                                                      #
                                                                               #
    if(lcpix != 2){                                                            #
      print(paste0("ERROR: This function expects 'cpix' to have length 2. ",   #
                   "'cpix' has instead length ", lcpix, ". Returning NULL."))  #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    if((length(cdec) != 1) | (length(unc_dec) != 1) | (length(pix_scl != 1))){ #
      print(paste0("ERROR: This function expects 'cdec', 'unc_dec' and 'pix_s",#
                   "cl' to be length 1 constants. Returning NULL."))           #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    dec_map <- array(NA, dim = c(dims, 2),                                     #
                     dimnames = list(NULL, NULL, c("data", "unc")))            #
                                                                               #
    for(x in 1:dims[1]){                                                       #
      for(y in 1:dims[2]){                                                     #
                                                                               #
        if(narm == T & is.na(map[x,y])){                                       #
          next                                                                 #
        }                                                                      #
                                                                               #
        dify <- y - cpix[2]                                                    #
                                                                               #
        dec_map[x, y, "data"] <- cdec + dify * pix_scl                         #
        dec_map[x, y, "unc"] <- sqrt(2 * unc_dec^2)                            #
      }                                                                        #
    }                                                                          #
  }else{                                                                       #
    if(length(path_to_wcs_file) != 1){                                         #
      print(paste0("ERROR: 'path_to_wcs_file' must be a single string. ",      #
                   "Returning NULL."))                                         #
      return(NULL)                                                             #
    }                                                                          #
    if(!file.exists(path_to_wcs_file)){                                        #
      print(paste0("ERROR: 'path_to_wcs_file' does not exist. ",               #
                   "Returning NULL."))                                         #
      return(NULL)                                                             #
    }                                                                          #
    if(length(unc_dec) != 1){                                                  #
      print(paste0("ERROR: This function expects 'unc_dec' to be 1. ",         #
                   " Returning NULL."))                                        #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    dec_map <- array(NA, dim = c(dims, 2),                                     #
                     dimnames = list(NULL, NULL, c("data", "unc")))            #
                                                                               #
    for(x in 1:dims[1]){                                                       #
      t_coords <- c(rep(x, times = dims[2]), 1:dims[2])                        #
      dec_map[x,,"data"] <- convert_pixel_to_wcs(array(t_coords,               #
                                                       dim = c(dims[2], 2)),   #
                                                 path_to_wcs_file)[,"dec"]     #
      dec_map[x,,"unc"] <- unc_dec                                             #
    }                                                                          #
    rm(t_coords)                                                               #
  }                                                                            #
                                                                               #
  return(dec_map)                                                              #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'ra_map' of right ascensions; a 2D 'dec_map' of declinations;     #
# 'unc_ra' and 'unc_dec' which can either be 2D maps, or length 1 constants,   #
# holding the respective uncertainties of the previous maps;                   #
# and a reference coordinate pair 'cpix',                                      #
# Returns an 'output_map', a 3D array where the 'data' layer holds the         # 
# difference in the right ascension of each pixel compared to the reference    #
# pixel, and the 'unc' layer holding the respective uncertainties.             #
# If 'narm' is set to TRUE, pixels without numerical values will be skipped,   #
# and will hold NA in the 'output_map' as well.                                #
# All input angles are expected to be given in degrees.                        #
################################################################################
create_dRA_map <- function(ra_map, dec_map, unc_ra = 0, unc_dec = 0,           #
                           cpix = dim(map) %/% 2 + 1, narm = F){               #
                                                                               #
  ra_dims <- dim(ra_map)                                                       #
  dec_dims <- dim(dec_map)                                                     #
                                                                               #
  ura_dims <- dim(unc_ra)                                                      #
  udec_dims <- dim(unc_dec)                                                    #
                                                                               #
  if(is.null(ra_dims)){                                                        #
    print(paste0("ERROR: 'ra_map' is not an array and this function expects i",#
                 "t to be a 2D array. Returning NULL."))                       #
  }                                                                            #
  if(is.null(dec_dims)){                                                       #
    print(paste0("ERROR: 'dec_map' is not an array and this function expects ",#
                 "it to be a 2D array. Returning NULL."))                      #
  }                                                                            #
                                                                               #
  ldims <- length(ra_dims)                                                     #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'ra_map' to be a 2D array. ", #
                   "'ra_map' has instead ", ldims,                             #
                   " dimension. Returning NULL."))                             #
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'ra_map' to be a 2D array. ", #
                   "'ra_map' has instead ", ldims,                             #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
  if(!prod(ra_dims == dec_dims)){                                              #
    print(paste0("ERROR: This function expects 'dec_map' to be a 2D array wit",#
                 "h the same dimensions as 'ra_map'. 'dec_map' has instead ",  #
                 "dimensions = ", dec_dims, ". Returning NULL."))              #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  lcpix <- length(cpix)                                                        #
                                                                               #
  if(lcpix != 2){                                                              #
    print(paste0("ERROR: This function expects 'cpix' to have length 2. ",     #
                 "'cpix' has instead length ", lcpix, ". Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  cra <- ra_map[cpix[1], cpix[2]]                                              #
  cdec <- dec_map[cpix[1], cpix[2]]                                            #
                                                                               #
  if((is.null(ura_dims) & (length(unc_ra) != 1)) |                             #
     (!prod(ura_dims == ra_dims))){                                            #
    print(paste0("ERROR: This function expects 'unc_ra' to either have the ",  #
                 "same dimensions as 'ra_map', or to have length 1. ",         #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
  if((is.null(udec_dims) & (length(unc_dec) != 1)) |                           #
     (!prod(udec_dims == ra_dims))){                                           #
    print(paste0("ERROR: This function expects 'unc_dec' to either have the ", #
                 "same dimensions as 'dec_map', or to have length 1. ",        #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(is.null(ura_dims) & (length(unc_ra) == 1)){                               #
    unc_ra <- array(unc_ra, dim = ra_dims)                                     #
  }                                                                            #
  if(is.null(udec_dims) & (length(unc_dec) == 1)){                             #
    unc_dec <- array(unc_dec, dim = ra_dims)                                   #
  }                                                                            #
                                                                               #
  unc_cra <- unc_ra[cpix[1], cpix[2]]                                          #
  unc_cdec <- unc_dec[cpix[1], cpix[2]]                                        #
                                                                               #
  output_map <- array(NA, dim = c(ra_dims, 2),                                 #
                      dimnames = list(NULL, NULL, c("data", "unc")))           #
                                                                               #
  if(length(unc_dec) == 1){                                                    #
    unc_mdec <- 0.5 * sqrt(unc_dec^2 + unc_cdec^2)^2                           #
  }                                                                            #
                                                                               #
  for(x in 1:ra_dims[1]){                                                      #
    for(y in 1:ra_dims[2]){                                                    #
                                                                               #
      if(narm == T & is.na(ra_map[x,y])){                                      #
        next                                                                   #
      }                                                                        #
                                                                               #
      mdec <- 0.5 * (dec_map[x, y] + cdec) * pi / 180 # to rads                #
      if(length(unc_dec) != 1){                                                #
        unc_mdec <- 0.5 * sqrt(unc_dec[x, y] + unc_cdec^2)^2 * pi / 180        #
      }                                                                        #
                                                                               #
      output_map[x, y, "data"] <- (ra_map[x, y] - cra) * cos(mdec)             #
      output_map[x, y, "unc"] <- sqrt(cos(mdec)^2 *                            #
                                        (unc_ra[x, y]^2 + unc_cra^2) +         #
                                        ((ra_map[x, y] - cra) *                #
                                           sin(mdec) * unc_mdec)^2)            #
    }                                                                          #
  }                                                                            #
  return(output_map)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 2D 'dec_map' of declinations; 'unc_dec' which can either be a 2D map,#
# or length 1 constant, holding the respective uncertainties of the previous   #
# map; and a reference coordinate pair 'cpix',                                 #
# Returns an 'output_map', a 3D array where the 'data' layer holds the         #
# difference in the declination of each pixel compared to the reference pixel, #
# and the 'unc' layer holding the respective uncertainties.                    #
# If 'narm' is set to TRUE, pixels without numerical values will be skipped,   #
# and will hold NA in the 'output_map' as well.                                #
# All input angles are expected to be given in degrees.                        #
################################################################################
create_dDEC_map <- function(dec_map, unc_dec = 0, cpix = dim(map) %/% 2 + 1,   #
                            narm = F){                                         #
                                                                               #
  dec_dims <- dim(dec_map)                                                     #
                                                                               #
  udec_dims <- dim(unc_dec)                                                    #
                                                                               #
  if(is.null(dec_dims)){                                                       #
    print(paste0("ERROR: 'dec_map' is not an array and this function expects ",#
                 "it to be a 2D array. Returning NULL."))                      #
  }                                                                            #
                                                                               #
  ldims <- length(dec_dims)                                                    #
                                                                               #
  if(ldims != 2){                                                              #
    if(ldims == 1){                                                            #
      print(paste0("ERROR: This function expects 'dec_map' to be a 2D array. ",#
                   "'dec_map' has instead ", ldims,                            #
                   " dimension. Returning NULL."))                             #
    }else{                                                                     #
      print(paste0("ERROR: This function expects 'dec_map' to be a 2D array. ",#
                   "'dec_map' has instead ", ldims,                            #
                   " dimensions. Returning NULL."))                            #
    }                                                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  lcpix <- length(cpix)                                                        #
                                                                               #
  if(lcpix != 2){                                                              #
    print(paste0("ERROR: This function expects 'cpix' to have length 2. ",     #
                 "'cpix' has instead length ", lcpix, ". Returning NULL."))    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  cdec <- dec_map[cpix[1], cpix[2]]                                            #
                                                                               #
  if((is.null(udec_dims) & (length(unc_dec) != 1)) |                           #
     (!prod(udec_dims == dec_dims))){                                          #
    print(paste0("ERROR: This function expects 'unc_dec' to either have the ", #
                 "same dimensions as 'dec_map', or to have length 1. ",        #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(is.null(udec_dims) & (length(unc_dec) == 1)){                             #
    unc_dec <- array(unc_dec, dim = dec_dims)                                  #
  }                                                                            #
                                                                               #
  unc_cdec <- unc_dec[cpix[1], cpix[2]]                                        #
                                                                               #
  output_map <- array(NA, dim = c(dec_dims, 2),                                #
                      dimnames = list(NULL, NULL, c("data", "unc")))           #
                                                                               #
  for(x in 1:dec_dims[1]){                                                     #
    for(y in 1:dec_dims[2]){                                                   #
                                                                               #
      if(narm == T & is.na(dec_map[x,y])){                                     #
        next                                                                   #
      }                                                                        #
                                                                               #
      output_map[x, y, "data"] <- (dec_map[x, y] - cdec)                       #
      output_map[x, y, "unc"] <- sqrt(unc_dec[x, y]^2 + unc_cdec^2)            #
    }                                                                          #
  }                                                                            #
  return(output_map)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given an array / constant 'dRA' with expected format either being a 2D or 3D #
# array, with two layers in the last dimension, holding right ascension        #
# differences (towards a reference pixel); 3 other arrays (with the same       #
# dimensions as 'dRa') / constants (with length 1) 'dDec', holding declination #
# difference (expected, but not tested, to be towards the same reference       #
# pixel), and the uncertainties 'unc_dRa' and 'unc_dDec'; and length 1         #
# constants 'phi' and 'unc_phi', which hold the inclination and respective     #
# uncertainty of a galaxy, in radians,                                         #
# Returns a 3D array 'beta' holding the beta values in the 'data' layer,       #
# and the respective uncertainties in the 'unc' layer.                         #
# this beta is based on Gupta et al 2016 and Gagliano et al 2021.              #
################################################################################
create_beta_array <- function(dRa, dDec, unc_dRa, unc_dDec, phi, unc_phi){     #
                                                                               #
  dims_dRa <- dim(dRa)                                                         #
                                                                               #
  if(!is.null(dims_dRa)){                                                      #
    dims_dDec <- dim(dDec)                                                     #
    dims_unc_dRa <-dim(unc_dRa)                                                #
    dims_unc_dDec <- dim(unc_dDec)                                             #
                                                                               #
    ldims_dRa <- length(dims_dRa)                                              #
    ldims_unc_dRa <- length(dims_unc_dRa)                                      #
    ldims_dDec <- length(dims_dDec)                                            #
    ldims_unc_dDec <- length(dims_unc_dDec)                                    #
                                                                               #
    ldRa <- length(dRa)                                                        #
    ldDec <- length(dDec)                                                      #
    lunc_dRa <- length(unc_dRa)                                                #
    lunc_dDec <- length(unc_dDec)                                              #
                                                                               #
    if((ldims_dDec != ldims_dRa) & (ldDec != 1)){                              #
      print(paste0("ERROR: This function expects 'dDec' to have the same numb",#
                   "er of dimensions as 'dRa' or to have length 1. 'dRa' has ",#
                   ldims_dRa, " dimensions and 'dDec' has ", ldims_dDec,       #
                   " and length ", ldDec, ". Returning NULL."))                #
      return(NULL)                                                             #
    }                                                                          #
    if((ldims_unc_dRa != ldims_dRa) & (lunc_dRa != 1)){                        #
      print(paste0("ERROR: This function expects 'unc_dRa' to have the same n",#
                   "umber of dimensions as 'dRa' or to have length 1. ",       #
                   "'dRa' has ", ldims_dRa, " dimensions and 'unc_dRa' has ",  #
                   ldims_unc_dRa, " and length ", lunc_dRa,                    #
                   ". Returning NULL."))                                       #
      return(NULL)                                                             #
    }                                                                          #
    if((ldims_unc_dDec != ldims_dRa) & (lunc_dDec != 1)){                      #
      print(paste0("ERROR: This function expects 'unc_dDec' to have the same ",#
                   "number of dimensions as 'dRa' or to have length 1.",       #
                   " 'dRa' has ", ldims_dRa, " dimensions and 'unc_dRa' has ", #
                   ldims_unc_dDec, " and length ", lunc_dDec,                  #
                   ". Returning NULL."))                                       #
      return(NULL)                                                             #
    }                                                                          #
                                                                               #
    if((ldDec != 1) & !prod(dims_dDec == dims_dRa)){                           #
      print(paste0("ERROR: This function expects 'dDec' to have the same ",    #
                   "dimensions as 'dRa' or to have length 1. ",                #
                   "'dRa' has dimensions = {", toString(dims_dRa),             #
                   "} and 'dDec' has dimensions = {", toString(dims_dDec),     #
                   "}. Returning NULL."))                                      #
      return(NULL)                                                             #
    }                                                                          #
    if((lunc_dRa != 1) & !prod(dims_unc_dRa == dims_dRa)){                     #
      print(paste0("ERROR: This function expects 'unc_dRa' to have the same ", #
                   "dimensions as 'dRa' or to have length 1. ",                #
                   "'dRa' has dimensions = {", toString(dims_dRa),             #
                   "} and 'unc_dRa' has dimensions = {",                       #
                   toString(dims_unc_dRa), "}. Returning NULL."))              #
      return(NULL)                                                             #
    }                                                                          #
    if((lunc_dDec != 1) & !prod(dims_unc_dDec == dims_dRa)){                   #
      print(paste0("ERROR: This function expects 'unc_dDec' to have the same ",#
                   "dimensions as 'dRa' or to have length 1. ",                #
                   "'dRa' has dimensions = {", toString(dims_dRa),             #
                   "} and 'unc_dDec' has dimensions = {",                      #
                   toString(dims_unc_dDec), "}. Returning NULL."))             #
      return(NULL)                                                             #
    }                                                                          #
  }else{                                                                       #
    ldRa <- length(dRa)                                                        #
    ldDec <- length(dDec)                                                      #
    lunc_dRa <- length(unc_dRa)                                                #
    lunc_dDec <- length(unc_dDec)                                              #
                                                                               #
    if((ldDec != ldRa) & (ldDec != 1)){                                        #
      print(paste0("ERROR: This function expects 'dDec' to have the same ",    #
                   "length as 'dRa' or to have length 1. 'dRa' has length ",   #
                   ldRa, " and 'dDec' has length ", ldDec,                     #
                   ". Returning NULL."))                                       #
      return(NULL)                                                             #
    }                                                                          #
    if((lunc_dRa != ldRa) & (lunc_dRa != 1)){                                  #
      print(paste0("ERROR: This function expects 'unc_dRa' to have the same ", #
                   "length as 'dRa' or to have length 1. 'dRa' has length ",   #
                   ldRa, " and 'unc_dRa' has length ", lunc_dRa,               #
                   ". Returning NULL."))                                       #
      return(NULL)                                                             #
    }                                                                          #
    if((lunc_dDec != ldRa) & (lunc_dDec != 1)){                                #
      print(paste0("ERROR: This function expects 'unc_dDec' to have the same ",#
                   "length as 'dRa' or to have length 1. 'dRa' has length ",   #
                   ldRa, " and 'unc_dDec' has length ", lunc_dDec,             #
                   ". Returning NULL."))                                       #
      return(NULL)                                                             #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(length(phi) != 1){                                                        #
    print(paste0("ERROR: This function expects 'phi' to have length 1",        #
                 ". Instead 'phi' has length ", length(phi),                   #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(unc_phi) != 1){                                                    #
    print(paste0("ERROR: This function expects 'unc_phi' to be a length 1 co", #
                 "nstant. Instead 'unc_phi' has length ", length(unc_phi),     #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  beta_data <- phi - atan(dDec/dRa)                                            #
  mag <- sqrt(dDec^2 + dRa^2)                                                  #
  beta_unc <- sqrt(unc_phi^2 + (dRa / mag * unc_dDec)^2 +                      #
                     (dDec / mag * unc_dRa)^2)                                 #
  rm(mag)                                                                      #
                                                                               #
  if(!is.null(dims_dRa)){                                                      #
    beta <- array(NA, dim = c(dims_dRa, 2),                                    #
                  dimnames = list(NULL, NULL, c("data", "unc")))               #
                                                                               #
    beta[,,"data"] <- beta_data                                                #
    beta[,,"unc"] <- beta_unc                                                  #
  }else{                                                                       #
    beta <- array(NA, dim = c(ldRa, 2),                                        #
                  dimnames = list(NULL, c("data", "unc")))                     #
                                                                               #
    beta[,"data"] <- beta_data                                                 #
    beta[,"unc"] <- beta_unc                                                   #
  }                                                                            #
                                                                               #
  return(beta)                                                                 #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given a 3D array 'beta_map', with two 2D maps of beta values, a given galaxy,#
# (first layer of third dimension) and uncertainties (second layer), and       #
# major axis 'a', minor axis 'b', and respective uncertainties 'unc_a' and     #
# 'unc_b' (these 4 should be length 1 constants and share the same units), of  #
# said galaxy,                                                                 #
# Returns a 3D map 'dlr' holding the direction light radius of each pixel in   #
# the 'data' layer, and the respective uncertainties in the 'unc' layer.       #
# This DLR is based on Gupta et al 2016 and Gagliano et al 2021.               #
################################################################################
create_DLR_map <- function(beta_map, a, b, unc_a, unc_b){                      #
                                                                               #
  dims_beta <- dim(beta_map)                                                   #
  ldims_beta <- length(dims_beta)                                              #
                                                                               #
  if(ldims_beta != 3){                                                         #
    print(paste0("ERROR: This function expects 'beta_map' to be a 3 D array. ",#
                 "'beta_map' has ", ldims_beta, " dimensions instead. ",       #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
  if(dims_beta[3] != 2){                                                       #
    print(paste0("ERROR: This function expects 'beta_map' third dimension to ",#
                 "have length 2. Instead 'beta_map' third dimension has ",     #
                 "length ", dims_beta[3], ". Returning NULL."))                #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(length(a) != 1){                                                          #
    print(paste0("ERROR: This function expects 'a' to have length 1",          #
                 ". Instead 'a' has length ", length(a),                       #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(unc_a) != 1){                                                      #
    print(paste0("ERROR: This function expects 'unc_a' to be length 1. ",      #
                 "Instead 'unc_a' has length ", length(unc_a),                 #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(b) != 1){                                                          #
    print(paste0("ERROR: This function expects 'b' to have length 1",          #
                 ". Instead 'b' has length ", length(b),                       #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(unc_b) != 1){                                                      #
    print(paste0("ERROR: This function expects 'unc_b' to be a length 1. ",    #
                 "Instead 'unc_b' has length ", length(unc_b),                 #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  ab_r <- a / b                                                                #
                                                                               #
  dlr <- array(NA, dim = dims_beta, dimnames = dimnames(beta_map))             #
  dlr[,, "data"] <- sqrt((ab_r * sin(beta_map[,,1]))^2 + cos(beta_map[,,1])^2) #
  dlr[,,"unc"] <- sqrt((a * sin(beta_map[,,1])^2 * unc_a /                     #
                          (b^2 * dlr[,,"data"]))^2 +                           #
                         (a^2 * sin(beta_map[,,1])^2 * unc_b /                 #
                            (b^3 * dlr[,,"data"]))^2 + ((a^2 - b^2) *          #
                                                          sin(beta_map[,,1]) * #
                            cos(beta_map[,,1]) * beta_map[,,2] /               #
                              (b^2 * dlr[,,"data"]))^2)                        #
                                                                               #
  return(dlr)                                                                  #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given two 3D arrays 'dRa_map' and 'dDec_map', holding 2D maps of RA and DEC, #
# differences relative to some pixel (reference pixel is expected, but not     #
# tested, to be shared between maps), in the first layer of the third dimension#
# and the respective uncertainties in the second layer,                        #
# Returns a 3D array 'offset_map' holding the RADEC offset (relative to said   #
# pixel) in the # 'data' layer, and the respective uncertainties in the 'unc'  #
# layer. No particular units are expected.                                     #
################################################################################
create_radec_offset_map <- function(dRa_map, dDec_map){                        #
                                                                               #
  dims_dRa <- dim(dRa_map)                                                     #
  dims_dDec <- dim(dDec_map)                                                   #
  ldims_dRa <- length(dims_dRa)                                                #
  ldims_dDec <- length(dims_dDec)                                              #
                                                                               #
  if((ldims_dRa != 3) | (ldims_dDec != 3)){                                    #
    print(paste0("ERROR: This function expects 'dRa_map' and 'dDec_map' to be",#
                 " 3D arrays. Instead 'dRa_map' is a ", ldims_dRa, "D array a",#
                 "nd 'dDec_map' is a ", ldims_dDec, "D array. Returning NULL.")#
          )                                                                    #
    return(NULL)                                                               #
  }                                                                            #
  if((dims_dRa[3] != 2) | (dims_dDec[3] != 2)){                                #
    print(paste0("ERROR: This function expects 'dRa_map' and 'dDec_map' third",#
                 " dimension to have length 2. Instead 'dRa_map' third ",      #
                 "dimension has length ", dims_dRa[3], " and 'dDec_map' third",#
                 " dimension has length ", dims_dDec[3], ". Returning NULL.")) #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!prod(dims_dRa == dims_dDec)){                                            #
    print(paste0("ERROR: This function expects 'dRa_map' and 'dDec_map' to ",  #
                 "have the same dimensions. Instead 'dRa_map' has ",           #
                 "dimensions = {", toString(dims_dRa), "}, and 'dDec_map' has",#
                 " dimensions = {", toString(dims_dDec), "}. Returning NULL."))#
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  offset_map <- array(NA, dim = dims_dRa,                                      #
                      dimnames = list(NULL, NULL, c("data", "unc")))           #
                                                                               #
  offset_map[,, "data"] <- sqrt(dRa_map[,,1]^2 + dDec_map[,,2]^2)              #
  offset_map[,, "unc"] <- 1 / offset_map[,, "data"] *                          #
    sqrt((dRa_map[,,1] * dRa_map[,,2])^2 + (dDec_map[,,1] * dDec_map[,,2])^2)  #
                                                                               #
  return(offset_map)                                                           #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Given two 3D arrays 'offset_map' and 'dlr_map', holding 2D maps of the       #
# RADEC offset and directional light radius, respectively, in the first layer  #
# of the third dimension and the respective uncertainties in the second layer, #
# the major axis 'a', and its uncertainty 'unc_a', of the galaxy to which those#
# maps refer to, 'offset_map', 'a' and 'unc_a' should share the same units.    #
# 'a' and 'unc_a' are expected to be length 1 constants.                       #
# Returns a 3D map 'norm_offset_map' holding the same offsets now normalized   #
# (and as such, unit-less) by the respective pixels DLR in the 'data' layer,   #
# and the respective uncertainties in the 'unc' layer.                         #
# Normalization is based on Gupta et al 2016 and Gagliano et al 2021.          #
################################################################################
create_radec_offset_norm_by_dlr_map <- function(offset_map, dlr_map, a, unc_a){#
                                                                               #
  dims_offset <- dim(offset_map)                                               #
  dims_dlr <- dim(dlr_map)                                                     #
  ldims_offset <- length(dims_offset)                                          #
  ldims_dlr <- length(dims_dlr)                                                #
                                                                               #
  if((ldims_offset != 3) | (ldims_dlr != 3)){                                  #
    print(paste0("ERROR: This function expects 'offset_map' and 'dlr_map' ",   #
                 "to be 3D arrays. Instead 'offset_map' is a ", ldims_offset,  #
                 "D array and 'dlr_map' is a ", ldims_dlr, "D array. ",        #
                 "Returning NULL."))                                           #
    return(NULL)                                                               #
  }                                                                            #
  if((dims_offset[3] != 2) | (dims_dlr[3] != 2)){                              #
    print(paste0("ERROR: This function expects 'offset_map' and 'dlr_map' ",   #
                 "third dimension to have length 2. Instead 'offset_map' ",    #
                 "third dimension has length ", dims_offset[3],                #
                 " and 'dlr_map' third dimension has length ", dims_dlr[3],    #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(!prod(dims_offset == dims_dlr)){                                          #
    print(paste0("ERROR: This function expects 'offset_map' and 'dlr_map' to ",#
                 "have the same dimensions. Instead 'offset_map' has ",        #
                 "dimensions = {", toString(dims_offset), "}, and 'dlr_map' h",#
                 "as dimensions = {", toString(dims_dlr), "}. Returning NULL.")#
          )                                                                    #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(length(a) != 1){                                                          #
    print(paste0("ERROR: This function expects 'a' to have length 1",          #
                 ". Instead 'a' has length ", length(a),                       #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(unc_a) != 1){                                                      #
    print(paste0("ERROR: This function expects 'unc_a' to be a length 1 co",   #
                 "nstant. Instead 'unc_a' has length ", length(unc_a),         #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  norm_offset_map <- array(NA, dim = dims_offset,                              #
                           dimnames = list(NULL, NULL, c("data", "unc")))      #
                                                                               #
  d_dlr <- a / dlr_map[,,1]                                                    #
  unc_d_dlr <- sqrt((unc_a / dlr_map[,,1])^2 +                                 #
                      (a * dlr_map[,,2] / dlr_map[,,1]^2)^2)                   #
                                                                               #
  norm_offset_map[,, "data"] <- offset_map[,,1] / d_dlr                        #
  norm_offset_map[,, "unc"] <- 1 / d_dlr * sqrt(offset_map[,,2]^2 +            #
                                                  norm_offset_map[,, "data"] * #
                                                  unc_d_dlr^2)                 #
                                                                               #
  return(norm_offset_map)                                                      #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# UV Bump Noll et al 2009                                                      #
# Based on python function of Duarte et al 2023                                #
# 'n' is the dust index                                                        #
# 'lambda' is wavelength in nanometers                                         #
################################################################################
Drude_profile <- function(n, lambda){                                          #
                                                                               #
  if(length(lambda) != 1){                                                     #
    print(paste0("ERROR: this function expects 'lambda' to be length 1. Ins",  #
                 "tead it is length ", length(lambda), ". Returning NULL."))   #
    return(NULL)                                                               #
  }                                                                            #
  if(length(n) != 3 & length(n) != 2){                                         #
    print(paste0("ERROR: this function expects 'n' to be at length 2 or 3:",   #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(n),               #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  l <- lambda * 10 # to angstroms                                              #
                                                                               #
  n_up <- n[1] + n[2]                                                          #
                                                                               #
  if(length(n) == 2){                                                          #
    n_down <- n[1] - n[2]                                                      #
  }else{                                                                       #
    n_down <- n[1] - n[3]                                                      #
  }                                                                            #
                                                                               #
  d <- c(NA, NA, NA)                                                           #
  d[1] <- (0.85 - 1.9 * n[1]) * (l * 350)^2 /((l^2-2175^2)^2 + (l * 350)^2)    #
  d[2] <- abs((0.85 - 1.9 * n_up) *                                            #
                (l * 350)^2 /((l^2-2175^2)^2 + (l * 350)^2) - d[1])            #
  d[3] <- abs((0.85 - 1.9 * n_down) *                                          #
                (l * 350)^2 /((l^2-2175^2)^2 + (l * 350)^2) - d[1])            #
                                                                               #
  return(d)                                                                    #
}                                                                              #
################################################################################

################################################################################
# UV Bump Noll et al 2009, derivative for uncertainty propagation              #
# 'lambda' is wavelength in nanometers                                         #
################################################################################
Drude_derivative <- function(lambda){                                          #
                                                                               #
  if(length(lambda) != 1){                                                     #
    print(paste0("ERROR: this function expects 'lambda' to be length 1. Ins",  #
                 "tead it is length ", length(lambda), ". Returning NULL."))   #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  l <- lambda * 10 # to angstroms                                              #
                                                                               #
  d <- -1.9 * (l * 350)^2 /((l^2-2175^2)^2 + (l * 350)^2)                      #
                                                                               #
  return(d)                                                                    #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Calzetti Attenutation Curve, from Calzetti et al 2000                        #
# Based on python function of Duarte et al 2023                                #
# 'lambda' is wavelength in nanometers                                         #
################################################################################
Calzetti_attenuation_curve <- function(lambda){                                #
  l <- lambda * 1e-3 # to microns                                              #
                                                                               #
  if(l > 0.63){                                                                #
    kl <- 2.659 * (-1.857 + 1.040 / l) + 4.05                                  #
  }else{                                                                       #
    kl <- 2.659 * (-2.156 + 1.509 / l - 0.198 / (l^2) + 0.011 / (l^3)) + 4.05  #
  }                                                                            #
                                                                               #
  return(kl)                                                                   #
}                                                                              #
################################################################################
  
#------------------------------------------------------------------------------#

################################################################################
# Modified Calzetti Attenutation Curve, Noll et al 2009 / Leja 2017            #
# Based on python function of Duarte et al 2023                                #
# 'lambda' is wavelength in nanometers                                         #
# 'n' is the dust index                                                        #
# param_in can be either 'tau' or 'A'                                          #
################################################################################
modified_Calzetti_attenuation_curve <- function(param_in, n, lambda){          #
                                                                               #
  if(length(lambda) != 1){                                                     #
    print(paste0("ERROR: this function expects 'lambda' to be length 1. Ins",  #
                 "tead it is length ", length(lambda), ". Returning NULL."))   #
    return(NULL)                                                               #
  }                                                                            #
  if(length(param_in) != 3 & length(param_in) != 2){                           #
    print(paste0("ERROR: this function expects 'param_in' to be at length ",   #
                 "2 or 3: position 1 - data, position 2 - upper unc, positio", #
                 "n 3 -  lower unc. Instead it is length ", length(param_in),  #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(n) != 3 & length(n) != 2){                                         #
    print(paste0("ERROR: this function expects 'n' to be at length 2 or 3:",   #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(n),               #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  l <- lambda * 10 # to angstroms                                              #
                                                                               #
  n_up <- n[2]                                                                 #
  par_in_up <- param_in[2]                                                     #
                                                                               #
  if(length(n) == 2){                                                          #
    n_down <- n[2]                                                             #
  }else{                                                                       #
    n_down <- n[3]                                                             #
  }                                                                            #
  if(length(param_in) == 2){                                                   #
    par_in_down <- param_in[2]                                                 #
  }else{                                                                       #
    par_in_down <- param_in[3]                                                 #
  }                                                                            #
                                                                               #
  param_out <- c(NA, NA, NA)                                                   #
  param_out[1] <- param_in[1] * (l / 5510)^n[1] / 4.05 *                       #
    (Calzetti_attenuation_curve(lambda) + Drude_profile(n, lambda)[1])         #
                                                                               #
  der_n <- modified_Calzetti_derivative_n(param_in, n, lambda)                 #
  der_in <- param_out[1] / param_in[1]                                         #
                                                                               #
  param_out[2] <- sqrt((der_n * n_up)^2 + (der_in * par_in_up)^2)              #
  param_out[3] <- sqrt((der_n * n_down)^2 + (der_in * par_in_down)^2)          #
                                                                               #
  return(param_out)                                                            #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Modified Calzetti Attenutation Curve, Noll et al 2009 / Leja 2017, derivative#
# relative to 'n' for uncertainty propagation                                  #
# 'lambda' is wavelength in nanometers                                         #
# 'n' is the dust index                                                        #
# param_in can be either 'tau' or 'A'                                          #
################################################################################
modified_Calzetti_derivative_n <- function(param_in, n, lambda){               #
                                                                               #
  if(length(lambda) != 1){                                                     #
    print(paste0("ERROR: this function expects 'lambda' to be length 1. Ins",  #
                 "tead it is length ", length(lambda), ". Returning NULL."))   #
    return(NULL)                                                               #
  }                                                                            #
  if(length(param_in) != 3 & length(param_in) != 2){                           #
    print(paste0("ERROR: this function expects 'param_in' to be at length ",   #
                 "2 or 3: position 1 - data, position 2 - upper unc, positio", #
                 "n 3 -  lower unc. Instead it is length ", length(param_in),  #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(n) != 3 & length(n) != 2){                                         #
    print(paste0("ERROR: this function expects 'n' to be at length 2 or 3:",   #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(n),               #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  l <- lambda * 10 # to angstroms                                              #
                                                                               #
  der_n <- param_in[1] * (l / 5510)^n[1] / 4.05 *                              #
    (log(l / 5510) * (Calzetti_attenuation_curve(lambda) +                     #
                        Drude_profile(n, lambda)[1]) +                         #
       Drude_derivative(lambda))                                               #
                                                                               #
  return(der_n)                                                                #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Based on python function of Duarte et al 2023                                #
# converts the optical depth 'tau' and dust index 'n' to Rv                    #
################################################################################
get_Rv_from_tau_and_n <- function(tau, n){                                     #
                                                                               #
  if(length(tau) != 3 & length(tau) != 2){                                     #
    print(paste0("ERROR: this function expects 'tau' to be at length 2 or 3:", #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(tau),             #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  if(length(n) != 3 & length(n) != 2){                                         #
    print(paste0("ERROR: this function expects 'n' to be at length 2 or 3:",   #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(n),               #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  n_up <- n[2]                                                                 #
  tau_up <- tau[2]                                                             #
                                                                               #
  if(length(n) == 2){                                                          #
    n_down <- n[2]                                                             #
  }else{                                                                       #
    n_down <- n[3]                                                             #
  }                                                                            #
  if(length(tau) == 2){                                                        #
    tau_down <- tau[2]                                                         #
  }else{                                                                       #
    tau_down <- tau[3]                                                         #
  }                                                                            #
                                                                               #
  Av <- modified_Calzetti_attenuation_curve(tau, n, 551)                       #
  Ab <- modified_Calzetti_attenuation_curve(tau, n, 445)                       #
                                                                               #
  Rv <- c(NA, NA, NA)                                                          #
  Rv[1] <- Av[1] / (Ab[1] - Av[1])                                             #
                                                                               #
  der_Av_tau <- Av[1] / tau[1]                                                 #
  der_Ab_tau <- Ab[1] / tau[1]                                                 #
  der_Av_n <- modified_Calzetti_derivative_n(tau, n, 551)[1]                   #
  der_Ab_n <- modified_Calzetti_derivative_n(tau, n, 445)[1]                   #
                                                                               #
  der_rv_tau <- (Ab[1] * der_Av_tau - Av[1] * der_Ab_tau) / ((Ab[1] - Av[1])^2)#
  der_rv_n <- (Ab[1] * der_Av_n - Av[1] * der_Ab_n) / ((Ab[1] - Av[1])^2)      #
                                                                               #
  Rv[2] <- sqrt((der_rv_tau * tau_up)^2 + (der_rv_n * n_up)^2)                 #
  Rv[3] <- sqrt((der_rv_tau * tau_down)^2 + (der_rv_n * n_down)^2)             #
                                                                               #
  return(Rv)                                                                   #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Based on python function of Duarte et al 2023                                #
# Converts the optical depth 'tau' to optical attenuation Av                   #
################################################################################
get_Av_from_tau <- function(tau){                                              #
                                                                               #
  if(length(tau) != 3 & length(tau) != 2){                                     #
    print(paste0("ERROR: this function expects 'tau' to be at length 2 or 3:", #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(tau),             #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  tau_up <- tau[1] + tau[2]                                                    #
                                                                               #
  if(length(tau) == 2){                                                        #
    tau_down <- tau[1] - tau[2]                                                #
  }else{                                                                       #
    tau_down <- tau[1] - tau[3]                                                #
  }                                                                            #
                                                                               #
  Av <- c(NA, NA, NA)                                                          #
  Av[1] <- -2.5 * log10(exp(-tau[1]))                                          #
  Av[2] <- abs(-2.5 * log10(exp(-tau_up)) - Av[1])                             #
  Av[3] <- abs(-2.5 * log10(exp(-tau_down)) - Av[1])                           #
                                                                               #
  return(Av)                                                                   #
}                                                                              #
################################################################################

#------------------------------------------------------------------------------#

################################################################################
# Converts the optical attenuation Av to optical depth 'tau'                   #
################################################################################
get_tau_from_Av <- function(Av){                                               #
                                                                               #
  if(length(Av) != 3 & length(Av) != 2){                                       #
    print(paste0("ERROR: this function expects 'Av' to be at length 2 or 3:",  #
                 " position 1 - data, position 2 - upper unc, position 3 - ",  #
                 " lower unc. Instead it is length ", length(Av),              #
                 ". Returning NULL."))                                         #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  Av_up <- Av[1] + Av[2]                                                       #
                                                                               #
  if(length(Av) == 2){                                                         #
    Av_down <- Av[1] - Av[2]                                                   #
  }else{                                                                       #
    Av_down <- Av[1] - Av[3]                                                   #
  }                                                                            #
                                                                               #
  tau <- c(NA, NA, NA)                                                         #
  tau[1] <- -log(10^(-Av[1] / 2.5))                                            #
  tau[2] <- abs(-log(10^(-(Av_up) / 2.5)) - tau[1])                            #
  tau[3] <- abs(-log(10^(-(Av_down) / 2.5)) - tau[1])                          #
                                                                               #
  return(tau)                                                                  #
}                                                                              #
################################################################################