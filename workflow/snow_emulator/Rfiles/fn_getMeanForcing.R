
fn_getMeanForcing <- function(iHru    = 1,
                              forcing = NA,
                              ncTime  = NA,
                              ti      = NA,
                              tf      = NA,
                              nYears  = NA){
  
  # ============================================================================
  # Load HRU data
  # ============================================================================
  hruForcing <- data.frame(time     = ncTime,
                           pptrate  = ncvar_get(nc = nc, varid = "pptrate",  start = c(iHru, 1), count = c(1, -1)),
                           LWRadAtm = ncvar_get(nc = nc, varid = "LWRadAtm", start = c(iHru, 1), count = c(1, -1)),
                           SWRadAtm = ncvar_get(nc = nc, varid = "SWRadAtm", start = c(iHru, 1), count = c(1, -1)),
                           airpres  = ncvar_get(nc = nc, varid = "airpres",  start = c(iHru, 1), count = c(1, -1)), 
                           airtemp  = ncvar_get(nc = nc, varid = "airtemp",  start = c(iHru, 1), count = c(1, -1)), 
                           spechum  = ncvar_get(nc = nc, varid = "spechum",  start = c(iHru, 1), count = c(1, -1)),
                           windspd  = ncvar_get(nc = nc, varid = "windspd",  start = c(iHru, 1), count = c(1, -1)))
  
  hruForcing$yyyy <- year(  hruForcing$time )
  hruForcing$mm   <- month( hruForcing$time )
  hruForcing$dd   <- day(   hruForcing$time )
  hruForcing$date <- as.Date( paste0(hruForcing$yyyy, "-", hruForcing$mm, "-", hruForcing$dd ))
  
  # Compute PET
  hruForcing$pet <- PE_Oudin(JD          = yday(hruForcing$date),
                             Temp        = hruForcing$airtemp - 273.15,
                             Lat         = forcing$lat[iHru],
                             LatUnit     = "deg",
                             TimeStepIn  = "hourly",
                             TimeStepOut = "hourly")
  
  # Obtain corresponding rows for all period, spring and winter
  rowsAll    <- which( hruForcing$date >= ti & hruForcing$time <= tf )
  rowsWinter <- which( hruForcing$date >= ti & hruForcing$time <= tf & (hruForcing$mm == 12 | hruForcing$mm == 1 | hruForcing$mm == 2))
  rowsSpring <- which( hruForcing$date >= ti & hruForcing$time <= tf & (hruForcing$mm == 3  | hruForcing$mm == 4 | hruForcing$mm == 5))
  
  # ============================================================================
  # Precipitation attributes
  # ============================================================================
  forcing <- forcing[iHru, ]
  iHru    <- 1 # because now there is only one row
  
  forcing$Pmean[iHru]   <- 3600 * sum( hruForcing$pptrate[rowsAll]   ) / nYears
  forcing$Pwinter[iHru] <- 3600 * sum( hruForcing$pptrate[rowsWinter]) / nYears
  forcing$Pspring[iHru] <- 3600 * sum( hruForcing$pptrate[rowsSpring]) / nYears 
  
  # ============================================================================
  # Temperature attributes
  # ============================================================================
  forcing$Tmean[iHru]   <- mean( hruForcing$airtemp[rowsAll   ] ) - 273.15
  forcing$Twinter[iHru] <- mean( hruForcing$airtemp[rowsWinter] ) - 273.15
  forcing$Tspring[iHru] <- mean( hruForcing$airtemp[rowsSpring] ) - 273.15
  
  # ============================================================================
  # Wind
  # ============================================================================
  forcing$WINDmean[iHru]    <- mean( hruForcing$windspd[rowsAll    ] )
  forcing$WINDwinter[iHru]  <- mean( hruForcing$windspd[rowsWinter ] )
  forcing$WINDspring[iHru]  <- mean( hruForcing$windspd[rowsSpring ] )
  
  # ============================================================================
  # Shortwave
  # ============================================================================
  forcing$SHORTWAVEmean[iHru]   <- mean( hruForcing$SWRadAtm[rowsAll    ] )
  forcing$SHORTWAVEwinter[iHru] <- mean( hruForcing$SWRadAtm[rowsWinter ] )
  forcing$SHORTWAVEspring[iHru] <- mean( hruForcing$SWRadAtm[rowsSpring ] )
  
  # ============================================================================
  # Rel. Humidity
  # ============================================================================
  forcing$SpecHumMean[iHru]     <- mean( hruForcing$spechum[rowsAll    ] )
  forcing$SpecHumWinter[iHru]   <- mean( hruForcing$spechum[rowsWinter ] )
  forcing$SpecHumspring[iHru]   <- mean( hruForcing$spechum[rowsSpring ] )
  
  # ============================================================================
  # Precipitation and temperature 
  # ============================================================================
  forcing$PETmean[iHru] <- sum( hruForcing$pet[rowsAll] ) / nYears
  forcing$Aridity[iHru] <- sum( hruForcing$pet[rowsAll] ) / sum( 3600 * hruForcing$pptrate[rowsAll] )
  
  # Snow fractions
  rowsAllBelowZero    <- which( hruForcing$airtemp <= (0+273.15) & hruForcing$date >= ti   &   hruForcing$time <= tf )
  rowsWinterBelowZero <- which( hruForcing$airtemp <= (0+273.15) & hruForcing$date >= ti   &   hruForcing$time <= tf   &   (hruForcing$mm == 12 | hruForcing$mm == 1 | hruForcing$mm == 2))
  rowsSpringBelowZero <- which( hruForcing$airtemp <= (0+273.15) & hruForcing$date >= ti   &   hruForcing$time <= tf   &   (hruForcing$mm == 3  | hruForcing$mm == 4 | hruForcing$mm == 5))
  
  forcing$snowFractMean[iHru]   <- sum( hruForcing$pptrate[rowsAllBelowZero] )    / sum( hruForcing$pptrate[rowsAll] )
  forcing$snowFractWinter[iHru] <- sum( hruForcing$pptrate[rowsWinterBelowZero] ) / sum( hruForcing$pptrate[rowsWinter] )
  forcing$snowFractSpring[iHru] <- sum( hruForcing$pptrate[rowsSpringBelowZero] ) / sum( hruForcing$pptrate[rowsSpring] )
  
  # p-seasonality
  pSeason    <- hruForcing[rowsAll, ]
  pSeason$jd <- yday( pSeason$date ) 
  pSeason    <- group_by(pSeason, jd)
  pSeason    <- summarise(.data = pSeason, pr = 3600 * sum(pptrate)/nYears, tas = mean(airtemp) - 273.15)
  
  # Optimiztion p-seasonality
  res_sce <- SCEoptim(FUN     = fn_pSeason,            
                      par     = runif(n = 4),          
                      lower   = c(0.1,  0.1, 0.1, 0.1),           
                      upper   = c(100,  365, 7,  365),          
                      control = list(ncomplex = 10),
                      #
                      pSeason = pSeason,
                      flagOpt = TRUE)
  
  forcing$seasonality[iHru] <- fn_pSeason(params  = res_sce$par,
                                          pSeason = pSeason,
                                          flagOpt = FALSE)
  
  # ============================================================================
  # High and low precipitation
  # ============================================================================
  forcingDaily <- hruForcing[rowsAll, ]
  forcingDaily <- group_by(forcingDaily, date)
  forcingDaily <- summarise(.data = forcingDaily, 
                            pr    = round(x = 3600 * sum(pptrate), digits = 2), 
                            tas   = mean(airtemp) - 273.15)
  
  forcing$highPrecFreq[iHru] <- sum(forcingDaily$pr >= (5 * mean(forcingDaily$pr))) / nYears # days / year
  forcing$lowPrecFreq[iHru]  <- sum(forcingDaily$pr <= 1                          ) / nYears # days / year
  
  forcingDaily$cum_pr    <- 0
  forcingDaily$cum_pr[1] <- forcingDaily$pr[1]
  forcingDaily$cumDays   <- 1
  
  prHigh                  <- 5 * mean(forcingDaily$pr) # mm/d
  prLow                   <- 1                         # mm/d
  forcingDaily$lowNevent  <- 0
  forcingDaily$highNevent <- 0
  
  lowNevent  <- 0
  highNevent <- 0
  
  for (i in 2:nrow(forcingDaily)) {
    forcingDaily$cum_pr[i]  <- forcingDaily$pr[i] + forcingDaily$cum_pr[i-1]
    forcingDaily$cumDays[i] <- forcingDaily$cumDays[i-1] + 1
    
    if(forcingDaily$pr[i] <= 1){ # 1 mm/d
      forcingDaily$cum_pr[i]  <- 0
      forcingDaily$cumDays[i] <- 0
    }
    
    # Check low events
    if(forcingDaily$pr[i] <= prLow){
      
      if(forcingDaily$pr[i] <= prLow & forcingDaily$pr[i-1] > prLow){
        lowNevent <- lowNevent + 1 
      }
      forcingDaily$lowNevent[i] <- lowNevent
    }else{
      forcingDaily$lowNevent[i] <- NA
    }
    
    # Check high events
    if(forcingDaily$cum_pr[i] >= prHigh){
      
      if(forcingDaily$cum_pr[i] >= prHigh & forcingDaily$cum_pr[i-1] < prHigh){
        highNevent <- highNevent + 1 
      }
      forcingDaily$highNevent[i] <- highNevent
    }else{
      forcingDaily$highNevent[i] <- NA
    }
  }
  
  # ============================================================================
  # Obtain the duration of the high/low events
  # ============================================================================
  rows      <- which(is.na(forcingDaily$lowNevent) == FALSE)
  lowEvents <- forcingDaily[rows, ] 
  lowEvents <- group_by(.data = lowEvents, lowNevent)
  lowEvents <- summarize(.data = lowEvents, consecDays = n())
  
  forcing$lowPrecDurMean[iHru] <- mean(lowEvents$consecDays)
  
  rows       <- which(is.na(forcingDaily$highNevent) == FALSE)
  highEvents <- forcingDaily[rows, ] 
  highEvents <- group_by(.data = highEvents, highNevent)
  highEvents <- summarize(.data = highEvents, consecDays = max(cumDays))
  
  forcing$highPrecDurMean[iHru] <- mean(highEvents$consecDays)
  
  # ============================================================================
  
  return(forcing)
}
