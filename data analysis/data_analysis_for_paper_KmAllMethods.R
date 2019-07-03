###################################
# Data analysis of morogoro (n=1773) double-sampling data set
#     For paper at Statistics of Biosciences
#
#      Tianchen Qian
#        2019.06.11
###################################


# Computes Kaplan-Meier estimators for the data analysis.

# Note: the data set is not publicly available. The code is just to showcase how the data analysis was conducted.

rm(list = ls())



setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.11_data_analysis_morogoro_n1773")

library(survival)
# library(tidyverse)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(numDeriv)
library(boot)

# Kaplan-Meier estimator
source("fcn_weighted-KM.R")


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

dta_original <- dta

info_0 <- c("age10", "prxcd4v")
info_loss <- c("l", "recent_CD4")

c_minus_l_max <- 1061

design_parallel <- expand.grid(t0_year = seq(from = 0.05, to = 2, by = 0.05), 
                               estimator = c("KmWeightedByPropScore", "KmWeightedByRobs", "KmCompleteCase", "KmStratifiedByRobs"), 
                               gamma = c(c_minus_l_max, 365, floor(365 * 1.5), 365 * 2))
design_parallel$estimator <- as.character(design_parallel$estimator)

varnames <- c("estimator", "t0_year", "gamma", "survProb", "se", "varlog", "ci_left", "ci_right")
output <- data.frame(matrix(NA, nrow = nrow(design_parallel), ncol = length(varnames)))
names(output) <- varnames


for (taskID in 1:nrow(design_parallel)) {
    cat(taskID, "")
    t0 <- 365 * design_parallel$t0_year[taskID]
    estimator <- design_parallel$estimator[taskID]
    gamma <- design_parallel$gamma[taskID]
    
    dta <- dta_original
    
    id <- which( (dta$s == 1) & (dta$c - dta$l > gamma) ) # id to mask
    dta$s[id] <- 0
    dta$x[id] <- NA
    dta$delta[id] <- NA
    
    if (estimator == "KmWeightedByPropScore") {
        result <- try(KmWeightedByPropScore(dta, info_0, info_loss, t0))    # survProb, varlog, ci
    } else if (estimator == "KmWeightedByRobs") {
        result <- try(KmWeightedByRobs(dta, t0))                            # survProb, varlog, ci
    } else if (estimator == "KmCompleteCase") {
        result <- try(KmCompleteCase(dta, t0))                              # survProb, se, ci
    } else if (estimator == "KmStratifiedByRobs") {
        set.seed(20190617)
        result <- try(KmStratifiedByRobs(dta, t0, n_boot = 2000))           # survProb, se, ci
    }
    
    output[taskID, "estimator"] <- estimator
    output[taskID, "t0_year"] <- design_parallel$t0_year[taskID]
    output[taskID, "gamma"] <- design_parallel$gamma[taskID]
    if (!(class(result) == "try-error")) {
        output[taskID, "survProb"] <- result$survProb
        output[taskID, "left_ci"] <- result$ci[1]
        output[taskID, "right_ci"] <- result$ci[2]
        output[taskID, "varlog"] <- ifelse(is.null(result$varlog), NA, result$varlog)
        output[taskID, "se"] <- ifelse(is.null(result$se), NA, result$se)
    }
    
    if (taskID %% 10 == 0) {
        write.csv(output, file = "result/KmAllMethods.csv", row.names = FALSE)
    }
}

write.csv(output, file = "result/KmAllMethods.csv", row.names = FALSE)

