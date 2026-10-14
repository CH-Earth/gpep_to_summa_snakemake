
fn_paramToMetric <- function(x             = NA, # LHS params
                             Atts          = NA,
                             paramsNames   = NA,
                             colAtts       = NA,
                             rf            = NULL,
                             emulator      = "RF",
                             emulatorType  = "SSE",
                             f             = 1){
  
  # Check if there are RF model
  if( is.null(rf) ){
    rf <- get("emulatorModel", envir = .GlobalEnv)
  }
  
  ##############################################################################
  # Order data
  ##############################################################################
  basinParam           <- as.data.frame(  matrix(data = x, nrow = 1)  )      # Obtain the parameter set for the i_basin
  colnames(basinParam) <- paramsNames
  
  ##############################################################################
  # Add attributes or not? CHeck if SSE or LSE will be considered.
  ##############################################################################
    
  # Prepare parameter set for emulator
  if(emulatorType == "SSE" ){
    param       <- basinParam 
    cols_train  <- paramsNames
  }else{
    basinAtt   <- Atts[, colAtts]
    param      <- cbind(basinAtt, basinParam) 
    cols_train <- c(paramsNames, colAtts)
  }

  ##############################################################################
  # Obtain estimate based on emulator
  ##############################################################################
  
  # ==========================================================================
  if(emulator == "RF"){
  # ==========================================================================
      
    rf.value <- predict(object = rf, data = param[, cols_train], num.threads = 1)
    y.est    <- rf.value$predictions   
  
  # ==========================================================================  
  }else if(emulator == "GBM"){
  # ==========================================================================
      
    param$targetMetric  <- 0 # Target Metric!!!!! Dummy variable.
    df_encoded          <- model.matrix(~ . - 1, data = param[, cols_train])  
    dtrain              <- xgb.DMatrix(data = df_encoded, label = param[, c("targetMetric")])
    rf.value            <- predict(rf, newdata = dtrain, nthread = 1)
    y.est               <- rf.value[1]     
    
  # ==========================================================================  
  }else{
  # ==========================================================================
    # Do nothing
  }
    
  ##############################################################################
  # Return only one variable for single objective calibration
  ##############################################################################
  if( is.na(       median(y.est, na.rm = TRUE)  ) == TRUE  |
      is.nan(      median(y.est, na.rm = TRUE)  ) == TRUE  |
      is.infinite( median(y.est, na.rm = TRUE)  ) == TRUE  |
      is.numeric(  median(y.est, na.rm = TRUE)  ) == FALSE ){
    y <- f * -1
  }else{
    y <- f * median(y.est, na.rm = TRUE)
  }
  
  return(y)
}