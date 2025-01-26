rm(list=ls())
# my_lib <- "/media/joaomfras/GalaxyPol/Pol-Gal/Polarimetric-Imaging-Reduction-Scripts-(FORS2)/Commit/process_lib.R"
my_lib <- "/home/joaomfras/Downloads/process_lib.R"
source(my_lib)
require(FITSio)
require(stringr)
require(jjb)
require(stats)
require(colorspace)
require(ggplot2)
require(plot3D)
require(spatstat)
require(data.table)
require(bigsnpr)
require(latex2exp)

# Estimate of distance from observer to target
# This should be queried in a future implementation
dist2Target <- 1 * 1e3
unc_dist2Target <- 0

axAng <- 0 * pi / 180 - pi / 2

if(axAng < 0){
  axAng <- axAng + pi
}
if(axAng > pi){
  axAng <- axAng - pi
}

axRatio <- .9
pxscl <- 0.1

# Define main input folder
# main_folder <- "/home/joaomfras/Desktop/New_Grain_Tests"
main_folder <- "/home/joaomfras/Desktop/Isotropic-Homogeneous_New_Round"
# Si_folder <- paste0(main_folder, "/Sil")##
# C_folder <- paste0(main_folder, "/Carb")##
# MB_folder <- paste0(main_folder, "/Zubko/Mix")##
# Zubko_folder <- paste0(main_folder, "/Zubko")##
# out_main <- paste0(main_folder, "/Param_vs_WL")
out_main <- "/Param_vs_WL/"
# MB_folder <- list.dirs(Zubko_folder, recursive = F)[2]
# MB_files <- list.files(MB_folder, full.names = T)

# dust_folders <- c(Si_folder, C_folder)##
# dust_folders <- Zubko_folder##

dust_folders <- list.dirs(main_folder, recursive = F, full.names = T)

# dust_outs <- paste0(out_main, "/Dist_Grain_(ZubkoMixManual)")#,  c("/Si", "/C"))##

wl_data <- round(seq_log(100, 2000, 30) / 1000, 2)
log_wl_data <- log10(wl_data)
wls <- paste0(wl_data, " ", TeX("$\\mu$"), "m")

for(dustit in 1:length(dust_folders)){
  dust_folder <- dust_folders[dustit]
  MB_folder <- dust_folder
  MB_files <- list.files(MB_folder, full.names = T)
  print(paste0("* Entering ", dust_folder, " *"))
  
  # size_folders <- list.dirs(dust_folder, recursive = F, full.names = T)[-2]##
  # sizes <- list.dirs(dust_folder, recursive = F, full.names = F)[-2]##
  # N_Sizes <- length(size_folders)
  
  output_folder <- paste0(dust_folder, out_main)
  
  # output_folder <- dust_outs[dustit]
  if(!dir.exists(output_folder)){
    mkdir(output_folder)
  }
  
  param_set1 <- c("polDeg", "polFlux", "polAng", "total")
  param_tit1 <- c("Pol. Degree", "Pol. Flux", "Pol. Angle", "Total Flux")
  param_set1_lab <- c("p(%)", "P", "X(º)_section", "I")
  
  param_set2 <- c("primaryscattered.fits", "primaryscattered1.fits", 
                  "primaryscattered2.fits", "multipleScattering")
  param_tit2 <- c("Scattered Flux (% of Total)", 
                  "Single Scattered Flux (% of scattered)", 
                  "Double Scattered Flux (% of scattered)", 
                  "Multiple Scattered Flux (% of scattered)")
  param_set2_lab <- c("Is (%)", "I1s (%)", "I2s (%)", "Ims (%)")
  
  N_pars <- length(param_set1)
  
  param_set1_pdf <- paste0(output_folder, "/Polarization_vs_WL.pdf")
  param_set2_pdf <- paste0(output_folder, "/Scattered-Flux_vs_WL.pdf")
  
  param_set1_paths <- NULL
  param_set1_paths$p <- list()
  param_set1_paths$P <- list()
  param_set1_paths$X <- list()
  param_set1_paths$I <- list()
  param_set2_paths <- NULL
  param_set2_paths$Is <- list()
  param_set2_paths$I1s <- list()
  param_set2_paths$I2s <- list()
  param_set2_paths$Ims <- list()
  
  print("** Listing polarization parameter file paths across all grain size folders **")
  for(p in 1:N_pars){
    par1 <- param_set1[p]
    par2 <- param_set2[p]
    
    param_set1_paths[[p]][[1]] <- MB_files[grep(par1, MB_files)]
    param_set2_paths[[p]][[1]] <- MB_files[grep(par2, MB_files)]
    # 
    # for(s in 1:N_Sizes){
    #   size_files <- list.files(size_folders[s], full.names = T)
    #   
    #   param_set1_paths[[p]][[1 + s]] <- size_files[grep(par1, size_files)]
    #   param_set2_paths[[p]][[1 + s]] <- size_files[grep(par2, size_files)]
    # }
  }
  
  temp <- readFITS(param_set1_paths[[1]][[1]])$imDat
  
  dimsObj <- dim(temp)
  rm(temp)
  
  cx <- dimsObj[1] / 2 + 1
  cy <- dimsObj[2] / 2 + 1
  
  nN <- dimsObj[1] * dimsObj[2]
  # types <- c("MB", sizes)
  types <- "MB"
  N_types <- length(types)
  
  nCoords <- array(NA, dim = c(nN, 2), dimnames = list(NULL, c("row", "col")))
  nMajaxs <- array(NA, dim = nN)
  nAngs <- array(NA, dim = nN)
  
  for(nx in 1:dimsObj[1]){
    for(ny in 1:dimsObj[2]){
      n <- ny + dimsObj[1] * (nx - 1)
          
      d <- sqrt((nx - cx)^2 + (ny - cy)^2)
      
      # Calculate major axis of ellipse that inscribes that pix
      t_ang <- atan2(ny - cy, nx - cx)
          
      if(t_ang < 0){
        t_ang <- t_ang + pi
      }
      if(t_ang > pi){
        t_ang <- t_ang - pi
      }
      nAngs[n] <- t_ang
      nMajaxs[n] <- sqrt(d^2 * (cos(t_ang - axAng)^2 + 
                                  (sin(t_ang - axAng) / axRatio)^2))
      nCoords[n, "row"] <- nx    
      nCoords[n, "col"] <- ny
    }
  }
  
  Xinds <- which(nAngs <= (pi/4) & nAngs >= 0, arr.ind = T)
  rm(nAngs)
  
  # Preparing labels for param stats arrays
  arr1_labs <- NULL
  arr2_labs <- NULL
  for(p in 1:N_pars){
    arr1_labs <- c(arr1_labs, param_set1_lab[p], 
                   paste0("unc_", param_set1_lab[p]))
    arr2_labs <- c(arr2_labs, param_set2_lab[p], 
                   paste0("unc_", param_set2_lab[p]))
  }
  
  len_pars <- length(arr1_labs)
  
  par_inds <- seq(1, len_pars, 2)
  
  # Holder arrays for whole map parameter averages
  pwl_arr1 <- array(NA, dim = c(N_types, len_pars, dimsObj[3]), 
                    dimnames = list(types, arr1_labs, wls))
  pwl_arr2 <- array(NA, dim = c(N_types, len_pars, dimsObj[3]), 
                    dimnames = list(types, arr2_labs, wls))
  
  # Holder arrays for random pixels
  samp_N <- 4
  pwl_rp1 <- array(NA, dim = c(samp_N, N_types, len_pars, dimsObj[3]), 
                   dimnames = list(NULL, types, arr1_labs, wls))
  pwl_rp2 <- array(NA, dim = c(samp_N, N_types, len_pars, dimsObj[3]), 
                   dimnames = list(NULL, types, arr2_labs, wls))
  
  print("** Sampling random pixels from maps of all parameters **")
  for(t in 1:N_types){
    Itot <- readFITS(param_set1_paths$I[[t]])$imDat
    
    for(p in 1:N_pars){
      data1 <- readFITS(param_set1_paths[[p]][[t]])$imDat
      data2 <- readFITS(param_set2_paths[[p]][[t]])$imDat
      
      if(p == 2 | p == 3){
        data2 <- data2 / readFITS(param_set2_paths[[1]][[t]])$imDat
      }
      
      par <- par_inds[p]
      unc_par <- par + 1
      
      tempI <- apply(Itot, 1:2, prod)
      null_inds <- which(as.vector(tempI) == 0, arr.ind = T)
      if(length(null_inds) != 0){
        tempI <- as.vector(tempI)[-null_inds]
      }
      
      if(length(tempI) == 0){
        next
      }
      
      if(t == 1 & p == 1){
        rnd_n <- sample.int(length(tempI), samp_N)
      }
      
      for(rn in 1:samp_N){
      
        for(wl in 1:dimsObj[3]){
          tempIt <- as.vector(Itot[,,wl])
          
          if(length(null_inds) != 0){
            tempIt <- tempIt[-null_inds]
          }
          
          d1 <- as.vector(data1[,,wl])
          d2 <- as.vector(data2[,,wl])
          d1 <- d1[-null_inds]
          d2 <- d2[-null_inds]
          
          if(p == 1){
            d2 <- d2 / tempIt
          }
          
          pwl_rp1[rn, t, par, wl] <- d1[rnd_n[rn]]
          pwl_rp2[rn, t, par, wl] <- d2[rnd_n[rn]] * 100
        }
      }
    }
  }
  
  print("** Getting whole map averages of all parameters **")
  for(t in 1:N_types){
    Itot <- readFITS(param_set1_paths$I[[t]])$imDat
    
    for(p in 1:N_pars){
      data1 <- readFITS(param_set1_paths[[p]][[t]])$imDat
      data2 <- readFITS(param_set2_paths[[p]][[t]])$imDat
      
      if(p == 2 | p == 3){
        data2 <- data2 / readFITS(param_set2_paths[[1]][[t]])$imDat
      }
      
      par <- par_inds[p]
      unc_par <- par + 1
      
      for(wl in 1:dimsObj[3]){
        tempI <- as.vector(Itot[,,wl])
        d1 <- as.vector(data1[,,wl])
        d2 <- as.vector(data2[,,wl]) 
        
        if(p == 1){
          d2 <- d2 / tempI
        }
        
        # Checking which pixels have flux, other params stats only make sense
        # where there actually is flux
        null_inds <- which(tempI == 0, arr.ind = T)
        
        if(length(null_inds) != 0){
          if(p != 4){
            d1 <- d1[-null_inds]
          }
          d2 <- d2[-null_inds]
        }
        rm(tempI)
        
        if(p == 3){
          tempI <- Itot[,,wl]
          d1 <- data1[,,wl]
          
          tempI <- tempI[Xinds]
          d1 <- d1[Xinds]
          
          null_inds <- which(tempI == 0, arr.ind = T)
          d1 <- d1[-null_inds] 
          rm(tempI, null_inds)
        }
        
        pwl_arr1[t, par, wl] <- median(d1, na.rm = T)
        pwl_arr1[t, unc_par, wl] <- sd(d1, na.rm = T)
        pwl_arr2[t, par, wl] <- mean(d2, na.rm = T) * 100
        pwl_arr2[t, unc_par, wl] <- sd(d2, na.rm = T) * 100
      }
    }
  }
  rm(data1, data2, Itot, d1, d2, null_inds)
  
  # Calculate minimum and maximum distance to center for pix in gPix
  min_M <- 0
  max_M <- max(nMajaxs)
  
  # Set a distance sequence
  print("** Defining distance sequence **")
  d_step <- 20
  unc_dist <- d_step / 2
  
  maj_list <- seq(min_M, max_M - unc_dist, d_step)
  unc_maj <- unc_dist
  
  kpc_min_M <- min_M * pxscl
  kpc_max_M <- max_M * pxscl
  
  kpc_maj_list <- maj_list * pxscl
  mN <- length(kpc_maj_list)
  kpc_tol_maj <- unc_maj * pxscl
  
  # Updating param stats arrays labels
  arr1_labs <- c(arr1_labs, "d", "unc_d")
  arr2_labs <- c(arr2_labs, "d", "unc_d")
  
  len_pars <- length(arr1_labs)
  # Holder for ring averages
  pd_arr1 <- array(NA, dim = c(N_types, mN, len_pars, dimsObj[3]), 
                   dimnames = list(types, NULL, arr1_labs, wls))
  pd_arr2 <- array(NA, dim = c(N_types, mN, len_pars, dimsObj[3]), 
                   dimnames = list(types, NULL, arr2_labs, wls))
  
  par_inds <- seq(1, len_pars - 2, 2)
  
  print("** Getting ring averages of all parameters **")
  for(t in 1:N_types){
    Itot <- readFITS(param_set1_paths$I[[t]])$imDat
    
    for(m in 1:mN){
      maj <- maj_list[m]
      
      if(m == 1 | m == mN){
        m_pixs <- which(nMajaxs == maj, arr.ind = T)
      }else{
        m_pixs <- which(nMajaxs < (maj + unc_maj) & 
                          nMajaxs >= (maj - unc_maj), arr.ind = T )
      }
      
      if(t == 1){
        pd_arr1[, m, "d",] <- mean(nMajaxs[m_pixs]) * pxscl
        pd_arr1[, m, "unc_d",] <- sd(nMajaxs[m_pixs]) * pxscl
        
        if(pd_arr1[1, m, "d",1] == 0){
          pd_arr1[, m, "unc_d",] <- 0
        }
        
        pd_arr2[, m, "d",] <- pd_arr1[, m, "d",]
        pd_arr2[, m, "unc_d",] <- pd_arr1[, m, "unc_d",]
      }
      
      coords <- nCoords[m_pixs,]
      
      # Selecting ring of flux
      ringItot <- array(NA, dim = c(length(m_pixs), dimsObj[3]))
      for(wl in 1:dimsObj[3]){
        tempI <- Itot[,,wl]
        
        if(length(m_pixs) == 1){
          ringItot[,wl] <- tempI[coords[1], coords[2]]
        }else{
          ringItot[,wl] <- tempI[coords]
        }
      }
      rm(tempI)
    
      for(p in 1:N_pars){
        data1 <- readFITS(param_set1_paths[[p]][[t]])$imDat
        data2 <- readFITS(param_set2_paths[[p]][[t]])$imDat
        
        if(p == 2 | p == 3){
          data2 <- data2 / readFITS(param_set2_paths[[1]][[t]])$imDat
        }
        
        par <- par_inds[p]
        unc_par <- par + 1
        
        # Selecting ring of data
        ringData1 <- array(NA, dim = c(length(m_pixs), dimsObj[3]))
        ringData2 <- array(NA, dim = c(length(m_pixs), dimsObj[3]))
        for(wl in 1:dimsObj[3]){
          tempD1 <- data1[,,wl]
          tempD2 <- data2[,,wl]
          
          if(length(m_pixs) == 1){
            ringData1[,wl] <- tempD1[coords[1], coords[2]]
            ringData2[,wl] <- tempD2[coords[1], coords[2]]
          }else{
            ringData1[,wl] <- tempD1[coords]
            ringData2[,wl] <- tempD2[coords]
          }
        }
        
        if(p == 3){
          mX_pixs <- intersect(m_pixs, Xinds)
          coordsX <- nCoords[mX_pixs,]
          
          # Intersect ring of data with double cone for X
          ringData1 <- array(NA, dim = c(length(mX_pixs), dimsObj[3]))
          ringIX <- array(NA, dim = c(length(mX_pixs), dimsObj[3]))
          for(wl in 1:dimsObj[3]){
            tempD1 <- data1[,,wl]
            tempI <- Itot[,,wl]
            
            if(length(mX_pixs) != 0){
              ringData1[,wl] <- tempD1[coordsX[1], coordsX[2]]
              ringIX[,wl] <- tempI[coordsX[1], coordsX[2]]
            }
          }
          rm(tempI)
        }
        rm(tempD1, tempD2, data1, data2)
        
        for(wl in 1:dimsObj[3]){
          tempI <- ringItot[,wl]
          null_inds <- which(tempI == 0, arr.ind = T)
          
          rd1 <- ringData1[,wl]
          rd2 <- ringData2[,wl]
          
          if(p == 1){
            rd2 <- rd2 / tempI
          }
          
          if(length(null_inds) != 0){
            if(p != 4){
              rd1 <- rd1[-null_inds]
            }
            rd2 <- ringData2[,wl]
            rd2 <- rd2[-null_inds]
          }
          rm(tempI)
          
          if(p == 3){
            tempI <- ringIX[,wl]
            rd1 <- ringData1[,wl]
            
            null_inds <- which(tempI == 0, arr.ind = T)
            if(length(null_inds) != 0){
              rd1 <- rd1[-null_inds]
            }
            rm(tempI, null_inds)
          }
          
          pd_arr1[t, m, par, wl] <- median(rd1, na.rm = T)
          pd_arr1[t, m, unc_par, wl] <- sd(rd1, na.rm = T)
          pd_arr2[t, m, par, wl] <- mean(rd2, na.rm = T) * 100
          pd_arr2[t, m, unc_par, wl] <- sd(rd2, na.rm = T) * 100
          rm(rd1, rd2, null_inds)
          
          if(is.na(pd_arr1[t, m, unc_par, wl])){
            pd_arr1[t, m, unc_par, wl] <- 0
          }
          if(is.na(pd_arr2[t, m, unc_par, wl])){
            pd_arr2[t, m, unc_par, wl] <- 0
          }
        }
      }
      # pd_arr1 <- pd_arr1[,-mN,,]
      # mN <- dim(pd_arr1)[2]
    }
    rm(ringItot, ringData1, ringData2)
  }
  
  # if(dust_folder == Si_folder){##
  #   legend_labs <- c(types[1], paste0("Si: ", types[-1], "um"))##
  # }else{##
  #   legend_labs <- c(types[1], paste0("C: ", types[-1], "um"))##
  # }##
  
  legend_labs <- c("Mix", types[-1])##
  
  print("** Creating plots of paramenters vs wavelength **")
  # Plot settings
  x_lims <- c(0, 1000)
  pdf_width <- 33
  pdf_height <- 33
  avg_col <- c('black', 'cyan', 'steelblue2', 'blue')
  dat_ch <- c(8, 21, 22, 24)
  ax_ind <- seq(1, length(wl_data), length.out = 10)
  x_lims <- c(min(log_wl_data), max(log_wl_data))
  
  pdf(param_set1_pdf, width = pdf_width, height = pdf_height)
  
  par(mfrow = c(2, 2))
  par(cex = 2)
  
  for(p in 1:N_pars){
    
    par <- par_inds[p]
    unc_par <- par + 1
    
    data <- pwl_arr1[, par,]
    unc <- pwl_arr1[, unc_par,]
    
    if(p == 1){
      y_max <- min(c(max(data, na.rm=T), 100))
    }else{
      y_max <- max(data, na.rm=T)
    }
    if(p == 3){
      y_max <- 90
      y_min <- -90
    }else{
      y_min <- max(c(min(data, na.rm=T) * .75, 0))
    }
    y_lims <- c(y_min, y_max)
    y_med <- mean(y_lims)
    
    matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
            ylab = param_set1_lab[p], xlim = x_lims, ylim = y_lims, 
            pch = dat_ch[1], col = avg_col[1], xaxt = 'n',
            main = paste0(param_tit1[p], " vs wavelength (nm) for the whole field"))
    axis(1, at = log_wl_data[ax_ind], labels = wls[ax_ind])
    
    # for(t in 2:N_types){
    #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
    # }
    legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
  }
  
  # for(rn in 1:samp_N){
  #   for(p in 1:N_pars){
  #     
  #     par <- par_inds[p]
  #     unc_par <- par + 1
  #     
  #     data <- pwl_rp1[rn,, par,]
  #     
  #     if(p == 1){
  #       y_max <- min(c(max(data, na.rm=T), 100))
  #     }else{
  #       y_max <- max(data, na.rm=T)
  #     }
  #     if(p == 3){
  #       y_max <- 90
  #       y_min <- -90
  #     }else{
  #       y_min <- max(c(min(data, na.rm=T) * .75, 0))
  #     }
  #     y_lims <- c(y_min, y_max)
  #     
  #     matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
  #             ylab = param_set1_lab[p], xlim = x_lims, ylim = y_lims, 
  #             pch = dat_ch[1], col = avg_col[1],
  #             main = paste0(param_tit1[p], " vs wavelength (nm) for random pixel #", rn))
  #     
  #     # for(t in 2:N_types){
  #     #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
  #     # }
  #     legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
  #   }
  # }
  
  
  for(m in 1:mN){
    for(p in 1:N_pars){
      par <- par_inds[p]
      unc_par <- par + 1
      
      data <- pd_arr1[, m, par,]
      unc <- pd_arr1[, m, unc_par,]
      
      if(p == 1){
        y_max <- min(c(max(data, na.rm=T), 100))
      }else{
        y_max <- max(data, na.rm=T)
      }
      if(p == 3){
        y_max <- 90
        y_min <- -90
      }else{
        y_min <- max(c(min(data, na.rm=T) * .75, 0))
      }
      y_lims <- c(y_min, y_max)
      
      if(is.infinite(y_lims[1]) | is.infinite(y_lims[2])){
        next
      }
      
      matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
              ylab = param_set1_lab[p], xlim = x_lims, ylim = y_lims, 
              pch = dat_ch[1], col = avg_col[1], xaxt = 'n',
              main = paste0(param_tit1[p], " vs wavelength (nm) for iso with majAx ",
                            round(pd_arr1[1, m, "d", 1], 2), "+-",
                            round(pd_arr1[1, m, "unc_d", 1], 2), " kpc"))
      
      axis(1, at = log_wl_data[ax_ind], labels = wls[ax_ind])
      
      # for(t in 2:N_types){
      #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
      # }
      legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
    }
  }
  dev.off()
  
  pdf(param_set2_pdf, width = pdf_width, height = pdf_height)
  
  par(mfrow = c(2, 2))
  par(cex = 2)
  
  for(p in 1:N_pars){
    
    par <- par_inds[p]
    unc_par <- par + 1
    
    data <- pwl_arr2[, par,]
    unc <- pwl_arr2[, unc_par,]
    
    y_max <- max(data, na.rm=T)
    y_min <- max(c(min(data, na.rm=T) * .75, 0))
    y_lims <- c(y_min, y_max)
    
    if(is.infinite(y_lims[1]) | is.infinite(y_lims[2])){
      next
    }
    matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
            ylab = param_set2_lab[p], xlim = x_lims, ylim = y_lims, 
            pch = dat_ch[1], col = avg_col[1], xaxt = 'n',
            main = paste0(param_tit2[p], " vs wavelength (nm) for the whole field"))
    
    axis(1, at = log_wl_data[ax_ind], labels = wls[ax_ind])
    
    # for(t in 2:N_types){
    #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
    # }
    legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
  }
  
  # for(rn in 1:samp_N){
  #   for(p in 1:N_pars){
  #     
  #     par <- par_inds[p]
  #     unc_par <- par + 1
  #     
  #     data <- pwl_rp2[rn,, par,]
  #     
  #     y_max <- max(data, na.rm=T)
  #     y_min <- max(c(min(data, na.rm=T) * .75, 0))
  #     y_lims <- c(y_min, y_max)
  #     
  #     if(is.infinite(y_lims[1]) | is.infinite(y_lims[2])){
  #       next
  #     }
  #     matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
  #             ylab = param_set2_lab[p], xlim = x_lims, ylim = y_lims, 
  #             pch = dat_ch[1], col = avg_col[1],
  #             main = paste0(param_tit2[p], " vs wavelength (nm) for random pixel #", rn))
  #     
  #     # for(t in 2:N_types){
  #     #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
  #     # }
  #     legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
  #   }
  # }
  
  
  for(m in 1:mN){
    for(p in 1:N_pars){
      par <- par_inds[p]
      unc_par <- par + 1
      
      data <- pd_arr2[, m, par,]
      unc <- pd_arr2[, m, unc_par,]
      
      y_max <- max(data, na.rm=T)
      y_min <- max(c(min(data, na.rm=T) * .75, 0))
      y_lims <- c(y_min, y_max)
      
      if(is.infinite(y_lims[1]) | is.infinite(y_lims[2])){
        next
      }
      matplot(log_wl_data, data, xlab = 'Wavelength (nm)', 
              ylab = param_set2_lab[p], xlim = x_lims, ylim = y_lims, 
              pch = dat_ch[1], col = avg_col[1], xaxt = "n",
              main = paste0(param_tit2[p], " vs wavelength (nm) for iso with majAx ",
                            round(pd_arr1[1, m, "d", 1], 2), "+-",
                            round(pd_arr1[1, m, "unc_d", 1], 2), " kpc"))
      
      axis(1, at = log_wl_data[ax_ind], labels = wls[ax_ind])
      # for(t in 2:N_types){
      #   points(wl_data, data[t, ], pch = dat_ch[t], col = avg_col[t])
      # }
      legend(2650, y_med, legend = legend_labs, pch = dat_ch, col = avg_col)
    }
  }
  dev.off()
}
