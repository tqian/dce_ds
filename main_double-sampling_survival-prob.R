###################################
# Deductive estimation of survival probability in double sampling design
# Main function
#      Tianchen Qian
#        2016.09.27
###################################

rm(list = ls())

library(survival)
library(foreach)
# library(doSNOW)
library(doMC)

setwd("~/Dropbox/Research/git_deductive/double_sampling/2016.09.22_rewrite_function/")
source("fcn_double-sampling_survival-prob.R")

dt <- readRDS("YannuData/cymorogororicd_lancet_analysis_recentCD4WeightIncluded.RDS")
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

dt$prxcd4v <- imputeByMean(dt$prxcd4v)
dt$recent_CD4 <- imputeByMean(dt$recent_CD4)
dt$recent_weight <- imputeByMean(dt$recent_weight)


set.seed(123)
dt$prxcd4v <- dt$prxcd4v + runif(nrow(dt)) * 0.3
set.seed(123)
dt$x <- dt$x + runif(nrow(dt), 0, .3)

info_0 <- c("age10", "prxcd4v")
info_loss <- c("l", "recent_CD4")

tmp <- sumEifSurvProbInDblSampling(dt, info_0, info_loss, time_question = 365, tuning_param_value = 1, eps = 1e-4, parallel_within_R = TRUE)
# 
# root <- uniroot(function(alpha) {
#     tmp <- sumEifSurvProbInDblSampling(dt, info_0, info_loss, time_question = 365, tuning_param_value = alpha, eps = 1e-4, parallel_within_R = TRUE)
#     print(alpha)
#     print(tmp$sumEif)
#     print(tmp$survProb_unperturbed)
#     return(tmp$sumEif)
# }, interval = c(1.4, 1.6))


c_minus_l_cutoffs <- c(1061)

for (cutoff in c_minus_l_cutoffs) {
    print(cutoff)
    dta <- dt
    id <- which( (dta$s == 1) & (dta$c - dta$l > cutoff) ) # id to mask
    dta$s[id] <- 0
    dta$x[id] <- NA
    dta$delta[id] <- NA
    
    root <- uniroot(function(alpha) {
        tmp <- sumEifSurvProbInDblSampling(dta, info_0, info_loss, time_question = 365, tuning_param_value = alpha, eps = 1e-4, parallel_within_R = TRUE)
        print(alpha)
        print(tmp$sumEif)
        print(tmp$survProb_unperturbed)
        filename <- paste0("result/1773pts_cutoff=", cutoff, "_extend-ver2_eps=1e-4_info0=age10,prxcd4v_infoloss=l,recentCD4.rds")
        saveRDS(tmp, filename)
        return(tmp$sumEif)
    }, interval = c(1.2, 2))
}



# gather results ----------------------------------------------------------

output <- data.frame(cutoff = 262, alpha = NA, sumEif = NA, survProb = NA, se = NA)

files <- list.files("result/")
files <- files[grep("recentCD4", files)]
setwd("result")
for (file in files) {
    tmp <- readRDS(file)
    cutoff <- sub("1773pts_cutoff=", "", file)
    cutoff <- sub("_extend-ver2_eps=1e-4_info0=age10,prxcd4v_infoloss=l,recentCD4.rds", "", cutoff)
    
    var_eif <- var(tmp$data$eif) / nrow(tmp$data)
    se_from_eif <- sqrt(var_eif)
    
    output <- rbind(output, c(as.numeric(cutoff), tmp$tuning_param_value, tmp$sumEif, tmp$survProb_unperturbed, se_from_eif))
}

names(output) <- c("cutoff", "alpha", "sumEif", "survProb", "se")
output <- output[order(output$cutof), ]

output$pctg_masked <- seq(from = 80, to = 0, by = -10)
output$survProb <- 1 - output$survProb

library(xtable)
rownames(output) <- NULL
output$survProb <- output$survProb * 100
output$se <- output$se * 100
latex_table <- xtable(output[, c("cutoff", "pctg_masked", "alpha", "survProb", "se")], digits = c(0, 0, 0, 2, 1, 2))
print(latex_table, include.rownames=FALSE)


