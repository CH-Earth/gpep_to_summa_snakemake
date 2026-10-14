#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

################################################################################
#
# Script to compute the SSE
#
# nicolas.vasquez@ucalgary.ca
#
################################################################################

# For debugging, use these lines:
# dirYamlFile  <- "/Users/nvasquez/Library/CloudStorage/Dropbox/UofC/10_SnowStations/lseParameters/Rfiles"
# yamlFileName <- "config_LSE.yaml" 
# nCores       <- 5
# basin        <- "chena"
# targetMetric <- "TSSv1"

library("lhs")      
library("parallel")
library("doParallel")
library("ncdf4")
library("stringr")
library("dplyr")
library("readxl")
library("yaml")
library("ranger")
library("xgboost")

dirYamlFile  <- args[1]
yamlFileName <- args[2]
nCores       <- args[3]
basin        <- args[4]
targetMetric <- args[5]

################################################################################
# Assign all required variables
################################################################################

# Read yaml file
config <- read_yaml( file.path(dirYamlFile, yamlFileName) )

# Directories
dirEmulators <- config$dirEmulators
dirFunctions <- config$dirFunctions
dirSave      <- config$dirSave
dirUseful    <- config$dirUseful
dirOutput    <- file.path(dirSave, basin)
dirAtts      <- file.path(config$dirAtts, basin)
# Emulator
iStep        <- config$iStep 
# Genetic Algorithm
gaPop        <- config$gaPop     # Population of the GA
gaIters      <- config$gaIters   # Number of iterations for the GA
# New parameter set
nNewParams   <- config$nNewParams
# ParamVersion
paramVersion <- config$paramVersion

# Create folder where results will be saved.
dir.create(dirOutput, recursive = TRUE)

# Load required function to use emulators
source( paste0( dirFunctions, "/fn_paramToMetric.R"))

################################################################################
# Load more info: files, etc.
################################################################################

# Read attributes file
attFileName <- paste0(basin, "_avgForcing_", config$versionAtts, ".RData")

load( file.path(dirAtts, attFileName) )

# Load parameters
snowParams <- read_xlsx(path = paste0(dirUseful, "/snow_parameters_", paramVersion, ".xlsx"))

################################################################################
# Fix the name of columns
################################################################################
fixCols <- data.frame(from = c("elev_m"),
                      to   = c("elev"  ))

for(iCol in 1:nrow(fixCols)){
  cols                       <- which( colnames(avgForcing) == fixCols$from[iCol] )
  colnames(avgForcing)[cols] <- fixCols$to[iCol]
}

################################################################################
################################################################################
#
#                            LSE iteration
#
################################################################################
################################################################################

AttsNames <- c("elev", "Pmean", "Pwinter", "Pspring", "highPrecDurMean", "highPrecFreq",
               "lowPrecDurMean", "lowPrecFreq", "Tmean", "Twinter", 
               "Tspring", "PETmean", "Aridity", "seasonality", "snowFractMean", 
               "snowFractWinter", "snowFractSpring", "WINDmean", "WINDwinter", 
               "WINDspring", "SHORTWAVEmean", "SHORTWAVEwinter", "SHORTWAVEspring", 
               "SpecHumMean", "SpecHumWinter", "SpecHumspring")

################################################################################
# Variables required by the function used to search the parameter space
################################################################################

# Maximize or minimize.
if(targetMetric == "rmse" | targetMetric == "nrmse"){
  f <- -1 # RMSE, minimize
}else{
  f <- 1  # Spearman correlation coefficient or Taylor Skill Score, maximize.
}

gaNcores <- FALSE # The parameter search won't be parallelized
emulator <- "RF"  # Is it random forest or GBM?

################################################################################
#                       L O A D   E M U L A T O R
################################################################################

load( paste0(dirEmulators, "/", emulator, "_step_", iStep, "_", targetMetric, ".RData") )
emulatorModel <- emulatorSummary$emulator
rm(emulatorSummary)

gc()

################################################################################
#       F U N C T I O N   T O   O B T A I N   P A R A M E T E R S
################################################################################

# Iterate for each HRU
fn <- function(iHru = 1){

  # ==========================================================================
  # Obtain attributes
  # ==========================================================================
  Atts  <- avgForcing[iHru, AttsNames]
  
  dirTo <- paste0(dirOutput, "/hruID_", avgForcing$hruID[iHru])
  dir.create(dirTo)
    
  if( file.exists( paste0(dirTo, "/", emulator, "_BestParamSet_targetMetric_", targetMetric, ".txt") ) == FALSE){
      
    #==========================================================================
    #              S E A R C H   P A R A M E T E R    S P A C E
    #==========================================================================
      
    # Monitor GA
    monitor_ga <- function(obj){
      iter            <- obj@iter
      population      <- obj@population                      # Get entire population
      fitness         <- obj@fitness                         # Get corresponding fitness values
      ga_log[[iter]] <<- data.frame(iteration  = iter,       # Store as a data frame in the global list
                                    population = population,
                                    fitness    = fitness )
    }
      
    # Search param Space
    ga_log <- list()
      
    # Genetic Algorithm (maximization)
    nsga_optim <- GA::ga(type          = "real-valued",
                         fitness       = fn_paramToMetric,
                         lower         = as.numeric(rep(0, times = nrow(snowParams) )),
                         upper         = as.numeric(rep(1, times = nrow(snowParams) )),
                         popSize       = gaPop,
                         maxiter       = gaIters,
                         suggestions   = randomLHS(n = gaPop, k = nrow(snowParams)),
                         run           = ceiling(2 * gaIters / 3),
                         optim         = FALSE,
                         parallel      = gaNcores,
                         monitor       = monitor_ga,
                         # Arguments passed to function
                         rf            = NULL,  # if "NULL", fn_paramToMetric reads it from global environment
                         emulator      = emulator,
                         emulatorType  = "LSE",
                         f             = f,               
                         paramsNames   = snowParams$name,
                         Atts          = Atts,     
                         colAtts       = AttsNames) 
      
    #--------------------------  Obtain best params --------------------------
    ga_newParams <- nsga_optim@population                              # Final population
    pop_fitness  <- nsga_optim@fitness                                 # Fitness(final population)
    rows         <- order(pop_fitness, decreasing = TRUE)              # Order parameters
    pop_fitness  <- pop_fitness[rows   ]                               # Reorder fitnes
    ga_newParams <- ga_newParams[rows, ]                               # Reorder from best to worse
      
    # This will be the new parameter set
    nParams   <- min(nrow(ga_newParams), nNewParams)
    lhsParams <- ga_newParams[1:nParams,  ]  # Leave only ga_params parameters
      
    #==========================================================================
    #               O B T A I N   P A R A M E T E R   V A L U E
    #==========================================================================
      
    tentativeParams <- data.frame(emulatedMetric = pop_fitness[1:nParams],
                                  params         = NA,
                                  targetMetric   = targetMetric)
      
    for(iParam in 1:nrow(tentativeParams)){
        
      snowParams$lhs   <- lhsParams[iParam, ]
      snowParams$toRun <- NA 
        
      for(i in 1:nrow(snowParams)){
          
        if( is.na(snowParams$transformation[i]) ){
          snowParams$toRun[i] <- snowParams$min[i] + snowParams$lhs[i] * (snowParams$max[i] - snowParams$min[i])
            
        }else if(snowParams$transformation[i] == "log10"){
          snowParams$toRun[i] <- 10^ (log10(snowParams$min[i]) + snowParams$lhs[i] * (log10(snowParams$max[i]) - log10(snowParams$min[i])))
            
        }else{
            # Do nothing
        }
      }
        
      tentativeParams$params[iParam] <- list(snowParams = snowParams)
    }
      
    #==========================================================================
    #                         S A V E   D A T A
    #==========================================================================
    
    # Save GA search
    save(ga_log,
         file = paste0(dirTo, "/", emulator, "_paramSearch_targetMetric_", targetMetric, ".RData"))
      
    # Save parameters
    save(tentativeParams, 
         file = paste0(dirTo, "/", emulator, "_allParams_targetMetric", targetMetric, ".RData"))
      
    # Isolate the best parameter set
    if(targetMetric == "rmse" | targetMetric == "nrmse"){
      rows <- which.min(tentativeParams$emulatedMetric) # Best value is the minimum
    }else{
      rows <- which.max(tentativeParams$emulatedMetric) # Best value is the maximum
    }
    
    paramsToSave <- tentativeParams$params[[rows]]
      
    write.table(x         = paramsToSave,
                file      = paste0(dirTo, "/", emulator, "_BestParamSet_targetMetric_", targetMetric, ".txt"),
                quote     = FALSE,
                sep       = "|",
                row.names = FALSE)
      
  } # Check if file already exists
  
  return(iHru)
  
} # Each HRU


################################################################################
# Parallelize the search
################################################################################

obtainParams <- mclapply(X              = 1:nrow(avgForcing),
                         FUN            = fn,
                         mc.cores       = nCores,
                         mc.preschedule = FALSE)


################################################################################
#                         The end. Good luck!
################################################################################
