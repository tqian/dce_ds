###################################
# Simulation study for the double-sampling design
#      Tianchen Qian
#        2016.11.17
###################################

rm(list = ls())

npara <- 100

version_fut <- as.integer(Sys.getenv("SGE_TASK_ID"))
on_cluster <- !is.na(version_fut)
if (on_cluster) {
    setwd("~/git_deductive/double_sampling/2016.09.22_rewrite_function/")
} else {
    version_fut <- npara
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2016.09.22_rewrite_function/")
}

# parallel job id and random seeds
parallel_seeds <- data.matrix(read.csv("misc/parallel_seeds.csv"))
taskID <- (version_fut %% npara) + 1
case_id <- (version_fut - 1) %/% npara + 1
.Random.seed <- parallel_seeds[taskID, ]


library(survival)
library(foreach)
# library(doSNOW)
library(doMC)

source("fcn_double-sampling_survival-prob.R")

expit <- function(x) { 1 / (1 + exp(-x)) }

# simulate data

parallel <- FALSE
n <- 500
ndb <- 50
t0 <- 0.15
eps <- 1e-4

if (case_id == 1) { # correct model
    info_0 <- c("z1", "z2")
    info_loss <- c("z3", "l")
} else if (case_id == 2) { # completely incorrect model
    info_0 <- c("z6")
    info_loss <- c("z7")
} else if (case_id == 3) { # incorrect T model, correct S model
    info_0 <- "z6"
    info_loss <- c("z3", "l")
} else if (case_id == 4) {
    info_0 <- c("z1", "z2")
    info_loss <- "z7"
}


jobname <- paste0("n=", n, "_ndb=", ndb, "_t0=", t0, "_eps=", eps,"_info0=", paste0(info_0, collapse = ','),
                  "_infoloss=", paste0(info_loss, collapse = ','))

filenames <- paste0("result/5 - simulation_R depends on zr, n=500, 4 scenarios/", jobname, "/job_", 1:npara, ".rds")

dir.create(paste0("result/simulation/", jobname, "/"), showWarnings = FALSE)

print(paste0("--- taskID = ", taskID, " ---"))


z1 <- runif(n) # info0
z2 <- runif(n) # info0
z3 <- runif(n) # infoloss
z4 <- runif(n)
z5 <- runif(n)
z6 <- runif(n) # exo
z7 <- runif(n) # exo
zr <- runif(n) # confounder associated with r and t

r <- rbinom(n, 1, (z1 + zr) / 2) # large z1, zr -- less likely to dropout

# larger z1, z2 -- live longer
ut <- runif(n)
t <- - log(ut) / ((1.5 - r) * exp(1 - (z1 + z2) / 2)) 

l_ <- runif(n, min = t/3, max = t)
l_[r == 1] <- NA

uc <- runif(n)
c <- - log(uc)

l <- x <- delta <- robs <- rep(as.numeric(NA), n)

# phase 1
for (i in 1:n) {
    if (r[i] == 1) {
        if (c[i] >= t[i]) { # (a)
            robs[i] <- 1; x[i] <- t[i]; delta[i] <- 1
        } else if (c[i] < t[i]) { # (b2)
            robs[i] <- 1; x[i] <- c[i]; delta[i] <- 0
        }
    } else if (r[i] == 0) {
        if (c[i] < l_[i]) { # (b1)
            robs[i] <- 1; x[i] <- c[i]; delta[i] <- 0
        } else if (c[i] >= t[i]) { # (c) or (e1)
            robs[i] <- 0; x[i] <- NA; delta[i] <- NA; l[i] <- l_[i]
        } else if (c[i] < t[i]) { # (d) or (e2)
            robs[i] <- 0; x[i] <- NA; delta[i] <- NA; l[i] <- l_[i]
        }
    }
}

# phase 2

e <- expit(-1 + z3 + 2 * l)
s <- rbinom(n, 1, e)
s[robs == 1] <- 0
# s[order(e, decreasing = T)[1:ndb]] <- 1

for (i in 1:n) {
    if (robs[i] == 0 & s[i] == 1) {
        if (c[i] >= t[i]) {
            x[i] <- t[i]; delta[i] <- 1
        } else {
            x[i] <- c[i]; delta[i] <- 0
        }
    }
}

dt <- data.frame(z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5, z6 = z6, z7 = z7, zr = zr,
                 robs = robs, l = l, c = c, x = x, delta = delta, s = s)

dt <- dt[order(dt$x), ]

root <- uniroot(function(alpha) {
    tmp <- sumEifSurvProbInDblSampling(dt, info_0, info_loss, time_question = t0, tuning_param_value = alpha, eps = eps, parallel_within_R = parallel)
    print(alpha)
    print(tmp$sumEif)
    print(tmp$survProb_unperturbed)
    return(tmp$sumEif)
}, interval = c(0.1, 4), extendInt = "yes")

result <- sumEifSurvProbInDblSampling(dt, info_0, info_loss, time_question = t0, tuning_param_value = root$root, eps = eps, parallel_within_R = parallel)
result$t <- t
result$l_ <- l_

saveRDS(result, filenames[taskID])




# analyze results ---------------------------------------------------------

if (0) {
    setwd("~/Dropbox/Research/git_deductive/double_sampling/2016.09.22_rewrite_function/")
    
    setwd("~/git_deductive/double_sampling/2016.09.22_rewrite_function/")
    
    tau0 <- tauhat <- sumeif <- se <- alpha <- numeric(npara)
    for (i in 1:npara) {
        if (file.exists(filenames[i])) {
            tmp <- readRDS(filenames[i])
        }
        tau0[i] <- mean(tmp$t >= t0)
        tauhat[i] <- tmp$survProb_unperturbed
        sumeif[i] <- tmp$sumEif
        se[i] <- sqrt(var(tmp$data$eif) / n)
        alpha[i] <- tmp$tuning
    }
    
    summary(tau0)
    summary(tauhat)
    summary(tau0 - tauhat)
    
    sqrt(mean((tauhat - mean(tau0))^2)) # RMSE
}
