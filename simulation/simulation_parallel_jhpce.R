###################################
# Simulation study for the double-sampling design
#      Tianchen Qian
#        2019.06.09
###################################

# This code is intended to run on JHPCE cluster (https://jhpce.jhu.edu/) using array job submission.
# It can be easily modified to run on a local computer. The entire simulation would take a long time (~34,000 hours) on a single core CPU.

# The index of estimators in this code that correspond to those in the paper is:
#     DE.Cox: 1
#     DE.LN: 3
#     PAR: 7
#     KM.S: 12
#     KM.C: 13
#     DE.Cox.WrongS: 8
#     DE.LN.WrongS: 9
#     DE.Cox(alpha=0): 2
#     DE.LN(alpha=0): 4

# Note: the PAR estimator (estimator #7 in this code) cannot be executed without the source code from Dr. Ming-Wen An.
#       I am not posting Dr. An's source code online as I do not have her permission to do so.


rm(list = ls())

# Parallel setup ----------------------------------------------------------

on_JHPCE <- TRUE
on_HarvardRC <- !on_JHPCE

nsim <- 1000
npara <- nsim

taskID_increase <- 0

if (on_JHPCE) {
    parallel_within_R <- FALSE
    version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))
    setwd("~/git_deductive/simulation20190608/")
} else if (on_HarvardRC) {
    parallel_within_R <- TRUE
    version_fut <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    setwd("~/deductive_simulation/")
} else {
    # on local computer
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.09_simulation")
}

parallel_seeds <- data.matrix(read.csv("parallel_seeds.csv"))
taskID <- ((version_fut - 1) %% npara) + 1 + taskID_increase
case_id <- (version_fut - 1) %/% npara + 1
.Random.seed <- parallel_seeds[taskID, ]


library(survival)
# library(tidyverse)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(numDeriv)
library(boot)

# deductive estimation
source("fcn_double-sampling_survival-prob.R") 
source("fcn_DE.R")

# Kaplan-Meier estimator
source("fcn_weighted-KM.R")

# parametric estimator by An et al. (2014)
setwd("functions_from_An/")
source("functions.parametric.R")
source("functions.getMLE.R")
source("functions.getData.R")
source("routines.r")
setwd("..")

# generative models
source("dgm.R")

design_simulation <- expand.grid(n = c(50, 100, 200, 500), dgm_type = 1)

design_simulation_extra <- expand.grid(n = c(50, 100, 200, 500), dgm_type = 5)
design_simulation <- rbind(design_simulation, design_simulation_extra)

design_simulation_extra <- expand.grid(n = c(50, 100, 200, 500), dgm_type = 6)
design_simulation <- rbind(design_simulation, design_simulation_extra)


info_0 <- "z"
info_loss <- c("l")
info_loss_cox <- info_loss
info_loss_lognormal <- info_loss

n <- design_simulation$n[case_id]
dgm_type <- design_simulation$dgm_type[case_id]

if (dgm_type == 1) {
    dgm_this_simulation <- dgm_start_from_R
    t0 <- 0.7
} else if (dgm_type == 2) {
    dgm_this_simulation <- dgm_discrete_var
    t0 <- 3.2
} else if (dgm_type == 4) {
    dgm_this_simulation <- dgm_Weilbull_LCT
    t0 <- 0.7
} else if (dgm_type == 5) {
    dgm_this_simulation <- dgm_trunclognormal_LCT
    t0 <- 0.7
    info_loss_cox <- "logL"
} else if (dgm_type == 6) {
    dgm_this_simulation <- dgm_start_from_R_trunclognormal
    t0 <- 0.7
    info_loss_cox <- "logL"
}

jobname <- paste0("n=", n, ",dgm_type=", dgm_type)
filenames <- paste0("result/", jobname, "/job_", 1:(npara + taskID_increase), ".csv")
dir.create(paste0("result/", jobname, "/"), showWarnings = FALSE)

print(paste0("--- case_id = ", case_id, " ---"))
print(paste0("--- taskID = ", taskID, " ---"))


# Simulation ------------------------------------------------------

dta <- dgm_this_simulation(n)

dta$logL <- log(dta$l) # for DE.Cox with dgm_lognormal_LCT()

# ordering of estimators:
# 1) coxph
# 2) coxph, alpha = 0
# 3) lognormal
# 4) lognormal, alpha = 0
# 5) weighted-KM by propensity score
# 6) weighted-KM by Robs
# 7) parametric estimator by An et al. (2014)
# 8) coxph, incorrect model for S (not including Z, L)
# 9) coxph, incorrect model for S (not including Z, L), alpha = 0
# 10) lognormal, incorrect model for S (not including Z, L)
# 11) lognormal, incorrect model for S (not including Z, L), alpha = 0
# 12) stratified-KM by Robs with bootstrap standard error
# 13) complete-case KM

varnames <- c("estimator", "survProb", "se", "varlog", "sumEif", "tuning_param_value", "ci_left", "ci_right")
output <- data.frame(matrix(NA, nrow = 13, ncol = length(varnames)))
names(output) <- varnames
output$estimator <- 1:13

eps <- 1e-4

write.csv(output, filenames[taskID], row.names = FALSE)


###### 1) coxph #####

result <- try(DE.Cox(dta = dta, info_0 = info_0, info_loss = info_loss_cox, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[1, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 2) coxph, alpha = 0 ######

result <- try(DE.Cox_fixAlpha(alpha = 0, dta = dta, info_0 = info_0, info_loss = info_loss_cox, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[2, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 3) lognormal #####

result <- try(DE.LN(dta = dta, info_0 = info_0, info_loss = info_loss, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[3, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 4) lognormal, alpha = 0 #####

result <- try(DE.LN_fixAlpha(alpha = 0, dta = dta, info_0 = info_0, info_loss = info_loss_cox, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[4, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 5) weighted-KM by propensity score #####

result <- try(KmWeightedByPropScore(dta, info_0, info_loss, t0))
if (!(class(result) == "try-error")) {
    output[5, c("survProb", "ci_left", "ci_right", "varlog")] <- c(result$survProb, result$ci[1], result$ci[2], result$varlog)
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 6) stratefied-KM by Robs #####

result <- try(KmWeightedByRobs(dta, t0))
if (!(class(result) == "try-error")) {
    output[6, c("survProb", "ci_left", "ci_right", "varlog")] <- c(result$survProb, result$ci[1], result$ci[2], result$varlog)
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 7) parametric estimator by An et al. (2014) #####

dta_forEM <- dta
names(dta_forEM) <- c("Robs", "cov1", "cov2", "L", "probS", "S", "X", "delta", "r", "t", "C", "logL")
dta_forEM <- dta_forEM[, c("L", "C", "X", "S", "Robs", "delta", "cov1", "cov2")]
dta_forEM$id <- 1:nrow(dta_forEM)

dta_forEM$Y <- NA
dta_forEM$Y[!is.na(dta_forEM$delta) & dta_forEM$delta == 1] <- dta_forEM$X[!is.na(dta_forEM$delta) & dta_forEM$delta == 1]

dir.create("em_tmpfile/", showWarnings = FALSE)
em_tmpfiles <- c(paste0("em_tmpfile/em-", version_fut, ".txt"),
                 paste0("em_tmpfile/mle-", version_fut, ".txt"),
                 paste0("em_tmpfile/gs-", version_fut, ".txt"))
if (any(file.exists(em_tmpfiles))) {
    file.remove(em_tmpfiles[file.exists(em_tmpfiles)])
}

result <- try({
    analyze(data = dta_forEM,
            EM = TRUE, MAXITER.EM = 5000,
            em.filename = em_tmpfiles[1],
            mle.filename = em_tmpfiles[2],
            tol.EM = 1e-7,
            gs.theta.t0 = TRUE,
            MAXITER.GS = 200,
            gs.thetat0.filename = em_tmpfiles[3],
            tol.GS = 1e-7,
            myTvals = t0, variance = TRUE, theta.init = NULL)
})
if (!(class(result) == "try-error")) {
    output[7, c("survProb", "se", "ci_left", "ci_right")] <- 
        c(1 - result$pT, sqrt(result$var), 1 - result$pT - 1.96 * sqrt(result$var), 1 - result$pT + 1.96 * sqrt(result$var))
}
if (any(file.exists(em_tmpfiles))) {
    file.remove(em_tmpfiles[file.exists(em_tmpfiles)])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 8) coxph, incorrect model for S (not including Z, L) #####

models <- list(model_t.robs0 = c(info_0, info_loss_cox),
               model_c.robs0 = c(info_0, info_loss_cox),
               model_s.robs0 = NULL,
               model_t.robs1 = c(info_0, info_loss_cox),
               model_c.robs1 = c(info_0, info_loss_cox))
result <- try(DE.Cox(dta = dta, info_0 = info_0, info_loss = info_loss_cox, models = models, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[8, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 9) coxph, incorrect model for S (not including Z, L), alpha = 0 #####

models <- list(model_t.robs0 = c(info_0, info_loss_cox),
               model_c.robs0 = c(info_0, info_loss_cox),
               model_s.robs0 = NULL,
               model_t.robs1 = c(info_0, info_loss_cox),
               model_c.robs1 = c(info_0, info_loss_cox))
result <- try(DE.Cox_fixAlpha(alpha = 0, dta = dta, info_0 = info_0, info_loss = info_loss_cox, models = models, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[9, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 10) lognormal, incorrect model for S (not including Z, L) #####

models <- list(model_t.robs0 = c(info_0, info_loss_lognormal),
               model_c.robs0 = c(info_0, info_loss_lognormal),
               model_s.robs0 = NULL,
               model_t.robs1 = c(info_0, info_loss_lognormal),
               model_c.robs1 = c(info_0, info_loss_lognormal))
result <- try(DE.LN(dta = dta, info_0 = info_0, info_loss = info_loss, models = models, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[10, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 11) lognormal, incorrect model for S (not including Z, L), alpha = 0 #####

models <- list(model_t.robs0 = c(info_0, info_loss_lognormal),
               model_c.robs0 = c(info_0, info_loss_lognormal),
               model_s.robs0 = NULL,
               model_t.robs1 = c(info_0, info_loss_lognormal),
               model_c.robs1 = c(info_0, info_loss_lognormal))
result <- try(DE.LN_fixAlpha(alpha = 0, dta = dta, info_0 = info_0, info_loss = info_loss, models = models, time_question = t0, eps = eps, parallel_within_R = parallel_within_R, survfit_method = "KM"))
if (!(class(result) == "try-error")) {
    output[11, c("survProb", "se", "sumEif", "tuning_param_value", "ci_left", "ci_right")] <-
        c(result$survProb, result$se, result$sumEif, result$tuning_param_value, result$ci[1], result$ci[2])
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 12) stratified-KM by Robs with bootstrap standard error #####

result <- try(KmStratifiedByRobs(dta, t0))
if (!(class(result) == "try-error")) {
    output[12, c("survProb", "ci_left", "ci_right", "se")] <- c(result$survProb, result$ci[1], result$ci[2], result$se)
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")


###### 13) complete-case KM #####

result <- try(KmCompleteCase(dta, t0))
if (!(class(result) == "try-error")) {
    output[13, c("survProb", "ci_left", "ci_right", "se")] <- c(result$survProb, result$ci[1], result$ci[2], result$se)
}
print(output)
write.csv(output, filenames[taskID], row.names = FALSE)
rm("result")



# analyze results ---------------------------------------------------------

if (0) {
    
    rm(list = ls())
    
    nsim <- 1000
    npara <- nsim
    
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.09_simulation")
    # setwd("~/git_deductive/simulation20190608/")
    # setwd("~/deductive_simulation/")
    
    
    library(dplyr)
    
    
    design_simulation <- expand.grid(n = c(50, 100, 200, 500), dgm_type = 1)
    
    # design_simulation_extra <- expand.grid(n = c(50, 100, 200, 500), dgm_type = 5)
    design_simulation_extra <- expand.grid(n = c(50, 100, 200), dgm_type = 5)
    design_simulation <- rbind(design_simulation, design_simulation_extra)
    design_simulation_extra <- expand.grid(n = c(50, 100, 200), dgm_type = 6)
    design_simulation <- rbind(design_simulation, design_simulation_extra)
    
    info_0 <- "z"
    info_loss <- c("l")
    
    for (case_id in 1:nrow(design_simulation)) {
        n <- design_simulation$n[case_id]
        dgm_type <- design_simulation$dgm_type[case_id]
        
        if (dgm_type == 1) {
            t0 <- 0.7
            tau_true <- 0.7698933
        } else if (dgm_type == 2) {
            t0 <- 3.2
            tau_true <- 0.6123776  
        } else if (dgm_type == 5) {
            t0 <- 0.7
            tau_true <- 0.589393
        }
        
        jobname <- paste0("n=", n, ",dgm_type=", dgm_type)
        filenames <- paste0("result/", jobname, "/job_", 1:npara, ".csv")
        
        result <- data.frame()
        for (ifile in filenames) {
            if (file.exists(ifile)) {
                tmp <- read.csv(ifile)
                result <- rbind(result, tmp)
            }
        }
        
        result <- as_tibble(result)
        
        if (nrow(result) >= 1) {
            performance <- result %>% group_by(estimator) %>%
                do(data.frame(bias = mean(.$survProb, na.rm = TRUE) - tau_true,
                              sd = sd(.$survProb, na.rm = TRUE),
                              cp = mean((.$ci_right >= tau_true) & (.$ci_left <= tau_true), na.rm = TRUE),
                              n_nonNA = sum(!is.na(.$survProb)),
                              n_NA = sum(is.na(.$survProb)))
                )
            performance$n <- n
            performance$dgm_type <- dgm_type
            
            if (case_id == 1) {
                performance_gathered <- performance
            } else {
                performance_gathered <- rbind(performance_gathered, performance)
            }
        }

    }
    
    performance_gathered <- performance_gathered[order(performance_gathered$estimator, performance_gathered$dgm_type, performance_gathered$n), ]
    performance_gathered <- performance_gathered[, c(8,7,1:6)]
    print(performance_gathered, n=Inf)
    print(filter(performance_gathered, dgm_type == 6), n=Inf)
    
    # ordering of estimators:
    # 1) coxph
    # 2) coxph, alpha = 0
    # 3) lognormal
    # 4) lognormal, alpha = 0
    # 5) weighted-KM by propensity score
    # 6) weighted-KM by Robs
    # 7) parametric estimator by An et al. (2014)
    # 8) coxph, incorrect model for S (not including Z, L)
    # 9) coxph, incorrect model for S (not including Z, L), alpha = 0
    # 10) lognormal, incorrect model for S (not including Z, L)
    # 11) lognormal, incorrect model for S (not including Z, L), alpha = 0
    # 12) stratified-KM by Robs with bootstrap standard error
    # 13) complete-case KM
    
    # # A tibble: 33 x 7
    # # Groups:   estimator [11]
    #    estimator      bias  n_NA     sd    cp     n dgm_type
    #        <int>     <dbl> <int>  <dbl> <dbl> <dbl>    <dbl>
    #  1         1  0.00272     44 0.0619 0.917    50        1
    #  2         1 -0.00718      2 0.0465 0.939   100        1
    #  3         1 -0.00142      1 0.0335 0.930   200        1
    #  4         2  0.00786      6 0.0615 0.933    50        1
    #  5         2 -0.00288      0 0.0450 0.945   100        1
    #  6         2 -0.000384     0 0.0332 0.94    200        1
    #  7         3 -0.00568      7 0.0648 0.922    50        1
    #  8         3 -0.00821      0 0.0465 0.935   100        1
    #  9         3 -0.00207      0 0.0337 0.93    200        1
    # 10         4 -0.0364       3 0.0823 0.827    50        1
    # 11         4 -0.0380       0 0.0573 0.805   100        1
    # 12         4 -0.0346       0 0.0408 0.74    200        1
    # 13         5 -0.0211       0 0.0704 0.965    50        1
    # 14         5 -0.0326       0 0.0511 0.93    100        1
    # 15         5 -0.0287       0 0.0380 0.9     200        1
    # 16         6 -0.0196       0 0.0686 0.965    50        1
    # 17         6 -0.0292       0 0.0498 0.935   100        1
    # 18         6 -0.0238       0 0.0370 0.925   200        1
    # 19         7 -0.0213       0 0.0615 0.889    50        1
    # 20         7 -0.0217       0 0.0428 0.895   100        1
    # 21         7 -0.0231       0 0.0308 0.84    200        1
    # 22         8 -0.000850    33 0.0629 0.922    50        1
    # 23         8 -0.00584      1 0.0457 0.940   100        1
    # 24         8 -0.000581     0 0.0332 0.95    200        1
    # 25         9  0.00786      6 0.0615 0.938    50        1
    # 26         9 -0.00288      0 0.0450 0.94    100        1
    # 27         9 -0.000384     0 0.0332 0.945   200        1
    # 28        10 -0.00458      7 0.0641 0.927    50        1
    # 29        10 -0.00604      0 0.0458 0.935   100        1
    # 30        10 -0.00110      0 0.0331 0.935   200        1
    # 31        11 -0.0363       3 0.0822 0.827    50        1
    # 32        11 -0.0380       0 0.0573 0.81    100        1
    # 33        11 -0.0346       0 0.0408 0.75    200        1
    
}
