rm(list=ls())
# my_lib <- "/media/joaomfras/GalaxyPol/Pol-Gal/Polarimetric-Imaging-Reduction-Scripts-(FORS2)/Commit/process_lib.R"
my_lib <- "/home/joaomfras/Downloads/process_lib.R"
source(my_lib)
require(FITSio)
main_folder <- "/home/joaomfras/Desktop/Isotropic-Homogeneous_New_Round"

mat_folders <- list.dirs(main_folder, recursive = F, full.names = T)

for(m_folder in mat_folders){
  # input_folders <- list.dirs(m_folder, recursive = F, full.names = T)
  input_folder <- m_folder
  # for(input_folder in input_folders){
    
    print(paste0("* Entering ", input_folder, " *"))
    
    all_files <- list.files(input_folder, full.names = T)
    
    is_scat <- grep("primaryscattered", all_files)
    
    Q_file <- all_files[grep("_stokesQ", all_files)]
    U_file <- all_files[grep("_stokesU", all_files)]
    I_file <- all_files[grep("_total", all_files)]
    Is_file <- all_files[is_scat[1]]
    Is1_file <- all_files[is_scat[2]]
    Is2_file <- all_files[is_scat[3]]
    
    Q <- readFITS(Q_file)$imDat
    U <- readFITS(U_file)$imDat
    I <- readFITS(I_file)$imDat
    Is <- readFITS(Is_file)$imDat
    Is1 <- readFITS(Is1_file)$imDat / Is
    Is2 <- readFITS(Is2_file)$imDat / Is
    
    dimsObj <- dim(Q)
    
    unc <- array(1, dim = dimsObj[1:2])
    P <- array(NA, dim = dimsObj)
    IP <- array(NA, dim = dimsObj)
    IsR <- array(NA, dim = dimsObj)
    X <- array(NA, dim = dimsObj)
    
    print("** Calculating P, IP, X and IsR **")
    for(wl in 1:dimsObj[3]){
      IP[,,wl] <- P_from_QU(Q[,,wl], U[,,wl], unc, unc)$Data
      P[,,wl] <- IP[,,wl] / I[,,wl] * 100
      IsR[,,wl] <- Is[,,wl] / I[,,wl] * 100
      X[,,wl] <- X_from_QU(Q[,,wl], U[,,wl], unc, unc)$Data
    }
    
    file_base_name <- strsplit(Q_file, "stokesQ")[[1]][1]
    
    P_file <- paste0(file_base_name, "polDeg(%).fits")
    IP_file <- paste0(file_base_name, "polFlux.fits")
    IsR_file <- paste0(file_base_name, "scatFluxRatio(%).fits")
    X_file <- paste0(file_base_name, "polAng(deg).fits")
    mIs_file <- paste0(file_base_name, "multipleScattering.fits")
    
    writeFITSim(P, P_file)
    writeFITSim(IP, IP_file)
    writeFITSim(IsR, IsR_file)
    writeFITSim(X, X_file)
    
    mult_Scat <- array(NA, dim = dimsObj)
    
    print("** Estimating multiple scattering **")
    for(wl in 1:dimsObj[3]){
      mult_Scat[,,wl] <- 1 - (Is1[,,wl] + Is2[,,wl])
    }
    
    writeFITSim(mult_Scat, mIs_file)
  # }
}