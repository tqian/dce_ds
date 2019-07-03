###################################
# Data analysis of morogoro (n=1773) double-sampling data set
#     For paper at Statistics of Biosciences
#
#      Tianchen Qian
#        2019.06.11
###################################

# Computes deductive estimator for the data analysis.

# This code is intended to run on HarvardRC cluster (https://www.rc.fas.harvard.edu/) using array job submission.
# It can be easily modified to run on a local computer. The entire simulation would take a long time (~13,000 hours) on a single core CPU.

# Note: the data set is not publicly available. The code is just to showcase how the data analysis was conducted.

rm(list = ls())

# Parallel setup ----------------------------------------------------------

on_HarvardRC <- TRUE
on_JHPCE <- !on_HarvardRC

if (on_JHPCE) {
    parallel_within_R <- FALSE
    version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))
    setwd("~/git_deductive/data_analysis20190611/")
} else if (on_HarvardRC) {
    parallel_within_R <- TRUE
    version_fut <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    setwd("~/data_analysis_morogoro/")
} else {
    # on local computer
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.11_data_analysis_morogoro_n1773")
}

taskID <- version_fut
# parallel_seeds <- data.matrix(read.csv("parallel_seeds.csv"))
# .Random.seed <- parallel_seeds[taskID, ]


library(survival)
# library(tidyverse)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(numDeriv)

# deductive estimation
source("fcn_double-sampling_survival-prob.R") 

# Kaplan-Meier estimator
source("fcn_weighted-KM.R")

# parametric estimator by An et al. (2014)
setwd("functions_from_An/")
source("functions.parametric.R")
source("functions.getMLE.R")
source("functions.getData.R")
source("routines.r")
setwd("..")



dt <- readRDS("data_morogoro_n1773/cymorogororicd_lancet_analysis_recentCD4WeightIncluded.RDS")
dt <- dt[dt$std > dt$ccpstartd, ] # restrict analysis to patients with ART start date > start of study
dt <- dt[, c("ptidno", "age10", "prxcd4v", "robs", "l", "x", "delta", "s", "c", "recent_CD4", "recent_weight")]
names(dt)[1] <- "PatientID"

# impute missing baseline and missing longitudinal
imputeByMean <- function(c) {
    num_missing <- sum(is.na(c))
    print(paste0("Among ", length(c), " obs, there are ", num_missing, " (", round(num_missing / length(c) * 100, 1), "%) missing obs."))
    mean_c <- mean(c, na.rm = TRUE)
    c[is.na(c)] <- mean_c
    return(c)
}

dt$prxcd4v <- imputeByMean(dt$prxcd4v) # "Among 1773 obs, there are 316 (17.8%) missing obs."
dt$recent_CD4 <- imputeByMean(dt$recent_CD4) # "Among 1773 obs, there are 113 (6.4%) missing obs."
dt$recent_weight <- imputeByMean(dt$recent_weight) # "Among 1773 obs, there are 3 (0.2%) missing obs."


# set.seed(123)
# dt$prxcd4v <- dt$prxcd4v + runif(nrow(dt)) * 0.3
# set.seed(123)
# dt$x <- dt$x + runif(nrow(dt), 0, .3)

dta <- dt
dta <- as_tibble(dta)

# dta$age <- as.integer(dta$age10 * 10)
# dta$prxcd4v_int <- as.integer(dta$prxcd4v)
# dta$recent_CD4_int <- as.integer(dta$recent_CD4)
# dta$x <- as.integer(dta$x)
# dta$delta <- as.integer(dta$delta)
# dta$s <- as.integer(dta$s)
# dta$l <- as.integer(dta$l)
# dta$robs <- as.integer(dta$robs)

dta$x[dta$x == 0] <- 0.5 # so that lognormal can be fitted

info_0 <- c("age10", "prxcd4v")
info_loss <- c("l", "recent_CD4")

eps <- 1e-4
parallel_within_R <- TRUE

c_minus_l_max <- 1061

# design_parallel <- expand.grid(t0_year = seq(from = 0.05, to = 2, by = 0.05), xdelta_working_model = c("coxph", "lognormal"), 
#                                gamma = c(365*2, floor(1.5 * 365), 365))

design_parallel <- expand.grid(t0_year = seq(from = 0.05, to = 2, by = 0.05), xdelta_working_model = c("coxph", "lognormal"), 
                               gamma = c(365*2, floor(1.5 * 365), 365))

t0 <- 365 * design_parallel$t0_year[taskID]
xdelta_working_model <- as.character(design_parallel$xdelta_working_model[taskID])
gamma <- design_parallel$gamma[taskID]


if (xdelta_working_model == "coxph") {
    uniroot_interval <- c(0.1, 4)
} else {
    uniroot_interval <- c(-2, 2)
}

result_filename <- paste0("t0year=", design_parallel$t0_year[taskID], ",xdelta=", xdelta_working_model, ",gamma=", gamma, ".csv")
output_filename <- paste0("result/", result_filename)
dir.create("result", showWarnings = FALSE)

parallel_within_R_log_filename <- paste0("log_parallel/log_parallel-", taskID, ".txt")
dir.create("log_parallel", showWarnings = FALSE)


id <- which( (dta$s == 1) & (dta$c - dta$l > gamma) ) # id to mask
dta$s[id] <- 0
dta$x[id] <- NA
dta$delta[id] <- NA

root <- try(
    {
        uniroot(function(alpha) {
            tmp <- sumEifSurvProbInDblSampling(dta, info_0 = info_0, info_loss = info_loss, time_question = t0, tuning_param_value = alpha,
                                               eps = eps, parallel_within_R = parallel_within_R, 
                                               xdelta_working_model = xdelta_working_model, survfit_method = "KM", parallel_within_R_log_filename = parallel_within_R_log_filename)
            cat("alpha =", alpha, "; sumEIF =", tmp$sumEif, "; survProb =", tmp$survProb_unperturbed, "; se =", tmp$se, "\n")
            return(tmp$sumEif)
        },
        interval = uniroot_interval,
        extendInt = "yes")
    })

if (!(class(root) == "try-error")) {
    result <- sumEifSurvProbInDblSampling(dta, info_0 = info_0, info_loss = info_loss, time_question = t0, tuning_param_value = root$root,
                                          eps = eps, parallel_within_R = parallel_within_R,
                                          xdelta_working_model = xdelta_working_model, survfit_method = "KM", parallel_within_R_log_filename = parallel_within_R_log_filename)
    result$data <- NULL
    output <- data.frame(survProb = result$survProb_unperturbed, se = result$se, sumEif = result$sumEif, tuning_param_value = result$tuning_param_value,
                         xdelta_working_model = xdelta_working_model, t0_year = design_parallel$t0_year[taskID], gamma = gamma)
} else {
    output <- data.frame(survProb = NA, se = NA, sumEif = NA, tuning_param_value = NA,
                         xdelta_working_model = xdelta_working_model, t0_year = design_parallel$t0_year[taskID], gamma = gamma)
}

write.csv(output, file = output_filename, row.names = FALSE)


if (0) {
    
    # collect results and save as a single CSV
    
    rm(list = ls())
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.11_data_analysis_morogoro_n1773")
    
    design_parallel <- expand.grid(t0_year = c(1, 2), xdelta_working_model = c("coxph", "lognormal"), gamma = c(1061, 566))
    design_parallel_extra <- expand.grid(t0_year = c(0.25, 0.5, 0.75, 1.25, 1.5, 1.75), xdelta_working_model = c("coxph", "lognormal"), gamma = c(1061))
    design_parallel <- rbind(design_parallel, design_parallel_extra)
    
    design_parallel$xdelta_working_model <- as.character(design_parallel$xdelta_working_model)
    
    
    
    have_read_files <- FALSE
    for (i in 1:nrow(design_parallel)) {
        
        result_filename <- paste0("t0year=", design_parallel$t0_year[i], ",xdelta=", design_parallel$xdelta_working_model[i], ",gamma=", design_parallel$gamma[i], ".csv")
        output_filename <- paste0("result/", result_filename)
        if (file.exists(output_filename)) {
            result <- read.csv(output_filename)
            if (!have_read_files) {
                result_gathered <- result
                have_read_files <- TRUE
            } else {
                result_gathered <- rbind(result_gathered, result)
            }
        }
    }
    if (names(result_gathered)[1] == "X") {
        result_gathered$X <- NULL
    }
    
    write.csv(result_gathered, file = "result/DE.csv")
    
}