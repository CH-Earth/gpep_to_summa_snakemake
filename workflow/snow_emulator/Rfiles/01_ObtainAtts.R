

################################################################################
#
# Script to obtain the attributes for the snow stations
#
# nicolas.vasquez@ucalgary.ca
#
################################################################################

# On Anvil
salloc -A ees250064 -N 1 -n 1 -c 21 --mem=30G -t 05:00:00 --partition=shared

module --force purge
module load r/4.4.1
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3
R

library("stringr")
library("lubridate")
library("dplyr")
library("airGR")
library("SoilHyP")
library("ncdf4")
library("parallel")

################################################################################
# Variables required for the script to run
################################################################################

basin      <- "bow"
dirForcing <- "/anvil/projects/x-ees240082/users/dcasson/gpep/bow/ensemble_generation/rf_best_regression_static_dynamic_backup/summa/forcing" 
dirSetting <- "/anvil/projects/x-ees240082/users/dcasson/gpep/bow/summa/settings"

basin      <- "chena"
dirForcing <- "/anvil/projects/x-ees240082/users/dcasson/gpep/chena/ensemble_generation/rf_best_regression_static_dynamic/summa/forcing" 
dirSetting <- "/anvil/projects/x-ees240082/users/dcasson/gpep/chena/summa/settings"

basin      <- "tuolumne"
dirForcing <- "/anvil/projects/x-ees240082/users/dcasson/gpep/tuolumne/ensemble_generation/rf_best_regression_static_dynamic/summa/forcing" 
dirSetting <- "/anvil/projects/x-ees240082/users/dcasson/gpep/tuolumne/summa/settings"


dirSave    <- "/home/x-nvasquez/ObtainSnowParams"

dir.create(file.path( dirSave, basin), recursive = TRUE)

nCores <- 20

source(paste0( dirSave, "/fn_getMeanForcing.R"))

# Load forcing data to understand structure
nc     <- nc_open( paste0( dirForcing, "/summa_mean_forcing.nc"))
ncAtts <- nc_open( paste0( dirSetting, "/attributes.nc")) 

forcing <- data.frame(hruID   = ncvar_get(nc = ncAtts, varid = "hruId"     ),
                      lat     = ncvar_get(nc = ncAtts, varid = "latitude"  ),
                      lon     = ncvar_get(nc = ncAtts, varid = "longitude" ),
                      elev_m  = ncvar_get(nc = ncAtts, varid = "elevation" ), 
                      area_m2 = ncvar_get(nc = ncAtts, varid = "HRUarea"   ))

ncTi <- nc$dim$time$units
ncTi <- str_remove(string = ncTi, pattern = "minutes since ")
ncTi <- as.POSIXlt(x = paste0(ncTi, " 00:00"))

# Compare hruIds
hruIds <- data.frame(forcing   = ncvar_get(nc = nc,     varid = "hruId" ),
                     setttings = ncvar_get(nc = ncAtts, varid = "hruId" ))


save(hruIds,
     file = file.path( dirSave, basin, "CheckHruIds.RData" ) )

################################################################################
# Attributes that will be computed
################################################################################

# Precipitation attributes
forcing$Pmean           <- 0 * NA
forcing$Pwinter         <- 0 * NA
forcing$Pspring         <- 0 * NA
forcing$highPrecDurMean <- 0 * NA
forcing$highPrecFreq    <- 0 * NA
forcing$lowPrecDurMean  <- 0 * NA
forcing$lowPrecFreq     <- 0 * NA
# Temperature attributes
forcing$Tmean           <- 0 * NA
forcing$Twinter         <- 0 * NA
forcing$Tspring         <- 0 * NA
# Prec. + Temp. attributes
forcing$PETmean         <- 0 * NA
forcing$Aridity         <- 0 * NA
forcing$seasonality     <- 0 * NA
forcing$snowFractMean   <- 0 * NA
forcing$snowFractWinter <- 0 * NA
forcing$snowFractSpring <- 0 * NA
# Wind
forcing$WINDmean        <- 0 * NA
forcing$WINDwinter      <- 0 * NA
forcing$WINDspring      <- 0 * NA
# Shortwave
forcing$SHORTWAVEmean   <- 0 * NA
forcing$SHORTWAVEwinter <- 0 * NA
forcing$SHORTWAVEspring <- 0 * NA
# Rel. Humidity
forcing$SpecHumMean     <- 0 * NA
forcing$SpecHumWinter   <- 0 * NA
forcing$SpecHumspring   <- 0 * NA

################################################################################
# Function to obtain seasonality
################################################################################

fn_pSeason <- function(params  = c(1, 1, 1, 1),
                       pSeason = NA,
                       flagOpt = FALSE){
  
  dP  <- params[1]
  sP  <- params[2]
  dT  <- params[3]
  sT  <- params[4]
  
  Tau <- 365.25
  
  pSeason$pEst   <- mean(pSeason$pr)  * (1 + dP * sin( 2 * pi * (pSeason$jd - sP) / Tau) )
  pSeason$tasEst <- mean(pSeason$tas) +      dT * sin( 2 * pi * (pSeason$jd - sT) / Tau)
  
  pSeasonality   <- dP * sign(dT) * cos(2 * pi * (sP - sT) / Tau)
    
  # Compute errors
  pr_rmse  <- sqrt( mean( (pSeason$pr  - pSeason$pEst  )^2) )
  tas_rmse <- sqrt( mean( (pSeason$tas - pSeason$tasEst)^2) )
  
  if(flagOpt){
    y <- sqrt( pr_rmse^2 + tas_rmse^2  )
  }else{
    y <- pSeasonality
  }
  
  return(y)
}

################################################################################
# Load data and obtain attributes
################################################################################

ncTime <- ncTi + 60 * ncvar_get(nc = nc, varid = "time")

ti     <- as.Date("2010-10-01", format = "%Y-%m-%d")
tf     <- as.Date("2018-09-30", format = "%Y-%m-%d")
nYears <- as.numeric( year(tf) - year(ti) )

avgForcing <- mclapply(X              = 1:nrow(forcing),
                       FUN            = fn_getMeanForcing,
                       mc.cores       = nCores,
                       mc.preschedule = FALSE,
                       #
                       forcing = forcing,
                       ncTime  = ncTime,
                       ti      = ti,
                       tf      = tf,
                       nYears  = nYears)

avgForcing <- do.call(rbind, avgForcing)


nc_close(nc)
nc_close(ncAtts)

# # #
save(avgForcing,  
     file = file.path( dirSave, basin,  paste0(basin, "_avgForcing.RData") ) )

