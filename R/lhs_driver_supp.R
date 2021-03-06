##==============================================================================
## lhs_driver_supp.R
##
## Lists the experiments and creates conditions to re-run the lhs.R script
## for each experiment with the different settings. These include:
## 1) %outbound permitted = .3, .35, .4, .45, .5 (precalibration constraint)
## 2) ... for each of using (2a) CO2 data only, (2b) temperature data only,
##    and (2c) both CO2 and temperature data together. That is, we enforce the
##    %outbound limit for both CO2 and temperature in experiment 2c, for
##    example.
## Runs the supplemental experiment for gymnosperm/angiosperm temperature control
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

## Clear workspace
rm(list=ls())

# needed packages
library(lhs)
library(Hmisc)
library(CholWishart)
library(abind)
library(MASS)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB-precalibration/R')
} else if(Sys.info()['user']=='aewsma') {
  machine <- 'office'
  setwd('/Users/aewsma/codes/GEOCARB-precalibration/R')
} else {
  # assume on another cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/GEOCARB-precalibration/R')
}

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 20000000
n_sample_per_chunk <- 10000 # maximum number of time series samples to consider at once
n_sample_min <- 10000 # minimum number of samples we'd be happy with; stop after this to avoid overrunning RAM
# if you don't want to use this, just set it greater than n_sample
n_node <- 1 # distribute the chunks across multiple nodes?

param_choice <- 'all'   # Calibrate all 68 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
fSR_choice <- 'DT2019'     # Which fSR time series? ("PR2011", "LENTON", "DT2019")
plot.dir <- '../figures/'

threshold_choices <- c(.5)
experiments <- expand.grid(threshold = threshold_choices,
                           use_co2 = c(FALSE, TRUE),
                           use_temperature = c(FALSE, TRUE))
# remove the experiments with both use_co2 = use_temperatures = FALSE
experiments <- experiments[-which(experiments$use_co2==FALSE & experiments$use_temperature==FALSE),]

# only doing the CO2+temperature experiment
for (idx_experiment in c(3)) {

  set.seed(1234*idx_experiment) # for reproducibility, but also for different samples for each experiment
  prcout_threshold <- experiments[idx_experiment,"threshold"]  # only save simulations with percent_outbound < this
  use_temperature <- experiments[idx_experiment,"use_temperature"]
  use_co2 <- experiments[idx_experiment,"use_co2"]

  # supplemental experiment modifying both GYM and timing
  appen <- "gym+timing"; supp_experiment_parameters <- c(110,60,0.25,-1)
  source("lhs_supp.R")
  Sys.sleep(180) # cool down for a few minutes

  # supplemental experiment modifying GYM but leaving timing alone
  appen <- "gym"; supp_experiment_parameters <- c(130,80,0.25,-1)
  source("lhs_supp.R")
  Sys.sleep(180)

  # supplemental experiment modifying both timing but leaving GYM alone
  appen <- "timing"; supp_experiment_parameters <- c(110,60,1,-1)
  source("lhs_supp.R")
  Sys.sleep(180)

  # supplemental experiment modifying deltaT2X for 130-40 Myr ago
  appen <- "dT2X"; supp_experiment_parameters <- c(130,80,1,1)
  source("lhs_supp.R")
  Sys.sleep(180)

  # supplemental experiment modifying deltaT2X for 130-40 Myr ago, and GYM+timing modifications
  appen <- "gym+timing+dT2X"; supp_experiment_parameters <- c(110,60,0.25,1)
  source("lhs_supp.R")
  Sys.sleep(180)
}


##==============================================================================
## End
##==============================================================================
