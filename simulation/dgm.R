###################################
# Generative models for simulation
#      Tianchen Qian
#        2019.06.04
###################################

## In the final paper, the dgm used are:

# GM 1 is dgm_start_from_R()
# GM 2 is dgm_start_from_R_trunclognormal()


expit <- function(x) { 1 / (1 + exp(-x)) }

rtlnorm <- function(n, meanlog, sdlog, truncpoint, maxiter = 100, ...) {
    # truncated lognormal using rejection sampling
    
    if (length(meanlog) == 1) {
        meanlog <- rep(meanlog, n)
    }
    if (length(sdlog) == 1) {
        sdlog <- rep(sdlog, n)
    }
    if (length(truncpoint) == 1) {
        truncpoint <- rep(truncpoint, n)
    }
    
    samples <- rlnorm(n = n, meanlog = meanlog, sdlog = sdlog, ...)
    if (all(samples <= truncpoint)) {
        all_truncated <- TRUE
    } else {
        all_truncated <- FALSE
    }
    
    iter <- 0
    while (iter < maxiter & !all_truncated) {
        id_to_truncate <- which(samples > truncpoint)
        samples[id_to_truncate] <- rlnorm(n = length(id_to_truncate), meanlog = meanlog[id_to_truncate], sdlog = sdlog[id_to_truncate], ...)
        if (all(samples <= truncpoint)) {
            all_truncated <- TRUE
        }
    }
    if (iter == maxiter & !all_truncated) {
        warning("rtlnorm: maxiter reached and not all truncated yet.")
    }
    return(samples)
}


dgm_discrete_var <- function(n) {
    # dgm_type = 2
    
    # A generative model with discrete values of z, l, t, c
    
    # This is to evaluate the performance of the estimate with discrete generative models.
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Baseline covariate
    dta$z <- base::sample((-2):2, n, replace = TRUE)
    
    # True dropout status
    dta$r <- rbinom(n, 1, (dta$z + 3) / 6)
    
    # Weilbull distributed survival time, then discretized
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- 1
    dta$t <- floor((- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu) * 5)
    
    # Uniform censoring time (constant accrual rate)
    dta$c <- floor(runif(n, min = 2, max = 8)) + 0.5
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # Uniform dropout time
    dta$l <- floor(runif(n, min = dta$t / 3, max = dta$t)) + 0.3
    dta$l[dta$r == 1] <- NA
    dta$l[dta$l > dta$c] <- NA
    
    # Observed dropout status
    dta$robs <- as.numeric(dta$l >= dta$x)
    dta$robs[is.na(dta$robs)] <- 1
    
    # Some other covariate "w" (to feed the function in place of "l")
    dta$w <- base::sample((-2):2, n, replace = TRUE)
    
    # Selection for double-sampling
    dta$probS <- expit(0.5 * dta$z + 0.5 * dta$l - 0.5)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA
    
    return(dta)
}

dgm_start_from_R <- function(n) {
    # dgm_type = 1
    
    # A generative model that starts from R (true dropout status).
    # This is a more plausible generative model.
    
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Baseline covariate
    dta$z <- runif(n, min = -2, max = 2)
    
    # True dropout status
    dta$r <- rbinom(n, 1, (dta$z + 3) / 6)
    
    # Weilbull distributed survival time
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- 1
    dta$t <- (- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu)
    
    # Uniform censoring time (constant accrual rate)
    dta$c <- runif(n, min = 0.5, max = 2)
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # Uniform dropout time
    dta$l <- runif(n, min = dta$t / 3, max = dta$t)
    dta$l[dta$r == 1] <- NA
    dta$l[dta$l > dta$c] <- NA
    
    # Observed dropout status
    dta$robs <- as.numeric(dta$l >= dta$x)
    dta$robs[is.na(dta$robs)] <- 1
    
    # Some other covariate "w" (to feed the function in place of "l")
    dta$w <- runif(n, min = 0, max = 1)
    
    # Selection for double-sampling
    dta$probS <- expit(0.5 * dta$z + 0.5 * dta$l + 0.5)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA

    return(dta)
}

if (0) {
    # summary statistics for dgm_start_from_R()
    set.seed(123)
    
    n <- 1000000
    
    dta <- dgm_start_from_R(n)
    
    1 - mean(dta$robs) # 0.42638
    
    mean(dta$s[dta$robs == 0]) # 0.6454266
    
    mean(dta$delta, na.rm = TRUE) # 0.6961395
    
    quantile(dta$x, probs = c(0.1, 0.9), na.rm = TRUE) # 0.5334267 1.1839490 
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.1, 0.9), na.rm = TRUE) # 0.5027582 0.7955256
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.25, 0.75), na.rm = TRUE) # 0.5576875 0.7336320 
    
    mean(dta$c <= dta$t) # 0.298438
    
    tau_true <- mean(dta$t > 0.7) # 0.770895
    
    # partial correlation conditional on Robs and Z
    library(ppcor)
    dta_test <- dta[, c("t", "c", "z", "robs")]
    pcor(dta_test, method = "pearson")
    
    # partial correlation, further conditional on L
    dta_test2 <- dta[dta$robs == 0, c("t", "c", "z", "l")]
    pcor(dta_test2, method = "pearson")
}


dgm_start_from_Robs <- function(n) {
    # dgm_type = 3
    
    # A generative model that starts from Robs (observed dropout status).
    # This is a less plausible generative model.
    # The purpose of this generative model is to evaluate the performance of the estimator
    # when the working models are correct.
    
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Observed dropout status
    dta$robs <- rbinom(n, 1, 0.6)
    
    # Baseline covariate
    dta$z <- runif(n, min = -2 + dta$robs, max = 2 + dta$robs)
    
    # Some other covariate "w" (to feed the function in place of "l")
    dta$w <- runif(n, min = 0, max = 1)
    
    # Weilbull distributed survival time
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- 1
    dta$t <- (- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu)
    
    # Uniform censoring time (constant accrual rate)
    dta$c <- runif(n, min = 0.5, max = 2)
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # # Uniform dropout time
    # dta$l <- runif(n, min = dta$t / 3, max = dta$t)
    # dta$l[dta$r == 1] <- NA
    
    # Selection for double-sampling
    dta$probS <- expit(0.5 * dta$z + 0.5 * dta$w)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA
    
    return(dta)
}


dgm_Weilbull_LCT <- function(n) {
    # dgm_type = 4
    
    # A generative model that starts from L,C,T for all patients
    # L,C,T all follow Weilbull distribution
    
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Baseline covariate
    dta$z <- runif(n, min = -2, max = 2)
    
    # Weilbull distributed survival time
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- 1
    dta$t <- (- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu)
    
    # Weilbull distributed censoring time
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- 1
    dta$c <- (- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu)
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # Weilbull distributed dropout time
    u <- runif(n)
    lambda <- 1 # scale parameter
    nu <- 5 # shape parameter
    beta <- -1
    dta$l <- (- log(u) / (lambda * exp(dta$z * beta))) ^ (1/nu)
    
    # True dropout status
    dta$r <- as.numeric(dta$l >= dta$t)
    
    # Observed dropout status
    dta$robs <- as.numeric(dta$l >= dta$x)
    
    # L missing for robs = 1
    dta$l[dta$robs == 1] <- NA

    
    # Some other covariate "w" (not used except as a place holder in the EM algorithm by An et al. (2014))
    dta$w <- runif(n, min = 0, max = 1)
    
    # Selection for double-sampling
    dta$probS <- expit(0.5 * dta$z + 0.5 * dta$l + 0.5)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA
    
    return(dta)
}

if (0) {
    # summary statistics for dgm_Weilbull_LCT()
    set.seed(123)
    
    n <- 1000000
    
    dta <- dgm_Weilbull_LCT(n)
    
    1 - mean(dta$robs) # 0.33321
    
    mean(dta$s[dta$robs == 0]) # 0.6961736
    
    mean(dta$delta, na.rm = TRUE) # 0.5000901
    
    quantile(dta$x, probs = c(0.1, 0.9), na.rm = TRUE) # 0.4873241 1.1618135
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.1, 0.9), na.rm = TRUE) # 0.5518518 0.8276378 
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.25, 0.75), na.rm = TRUE) # 0.6118453 0.7868786 
    
    mean(dta$c <= dta$t) # 0.500001
    
    tau_true <- mean(dta$t > 0.7) # 0.769462
    
    # set.seed(123)
    # 
    # n <- 500
    # dta <- dgm_Weilbull_LCT(n)
    
    # partial correlation conditional on Robs and Z
    library(ppcor)
    dta_test <- dta[, c("t", "c", "z", "robs")]
    pcor(dta_test, method = "pearson") # -0.08218737
    
    # partial correlation, further conditional on L
    dta_test2 <- dta[dta$robs == 0, c("t", "c", "z", "l")]
    pcor(dta_test2, method = "pearson") # 0.03011994
}



dgm_trunclognormal_LCT <- function(n, truncpoint = 7) {
    # dgm_type = 5
    
    # A generative model that starts from L,C,T for all patients
    # L,C,T all follow truncated lognormal distribution
    
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Baseline covariate
    dta$z <- runif(n, min = -2, max = 2)
    
    # lognormal distributed survival time
    beta <- 1
    sigma <- 0.5
    dta$t <- rtlnorm(n, meanlog = dta$z * beta, sdlog = sigma, truncpoint = truncpoint)
    
    # lognormal distributed censoring time
    beta <- 1
    sigma <- 0.5
    dta$c <- rtlnorm(n, meanlog = dta$z * beta, sdlog = sigma, truncpoint = truncpoint)
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # lognormal distributed dropout time
    beta <- -1
    sigma <- 0.5
    dta$l <- rtlnorm(n, meanlog = dta$z * beta, sdlog = sigma, truncpoint = truncpoint)
    
    # True dropout status
    dta$r <- as.numeric(dta$l >= dta$t)
    
    # Observed dropout status
    dta$robs <- as.numeric(dta$l >= dta$x)
    
    # L missing for robs = 1
    dta$l[dta$robs == 1] <- NA

    
    # Some other covariate "w" (not used except as a place holder in the EM algorithm by An et al. (2014))
    dta$w <- runif(n, min = 0, max = 1)
    
    # Selection for double-sampling
    dta$probS <- expit(-0.5 * dta$z + 0.5 * dta$l + 0.5)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA
    
    return(dta)
}

if (0) {
    # summary statistics for dgm_trunclognormal_LCT()
    set.seed(123)
    
    n <- 1000000
    
    dta <- dgm_trunclognormal_LCT(n, 7)
    
    1 - mean(dta$robs) # 0.46544
    
    mean(dta$s[dta$robs == 0]) # 0.5474605
    
    mean(dta$delta, na.rm = TRUE) # 0.5003864
    
    quantile(dta$x, probs = c(0.1, 0.9), na.rm = TRUE) # 0.1305549 2.7231948
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.1, 0.9), na.rm = TRUE) # 0.4215289 0.6753728
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.25, 0.75), na.rm = TRUE) # 0.4630679 0.6187772 
    
    mean(dta$c <= dta$t) # 0.49965
    
    hist(dta$t)
    tau_true <- mean(dta$t >= 0.7) # 0.589393
    print(tau_true)
    
    # set.seed(123)
    # 
    # n <- 500
    # dta <- dgm_lognormal_LCT(n)
    
    # partial correlation conditional on Robs and Z
    library(ppcor)
    dta_test <- dta[, c("t", "c", "z", "robs")]
    pcor(dta_test, method = "pearson") # 0.2104279
    
    # partial correlation, further conditional on L
    dta_test2 <- dta[dta$robs == 0, c("t", "c", "z", "l")]
    pcor(dta_test2, method = "pearson") # about 0.01485925
}



dgm_start_from_R_trunclognormal <- function(n, truncpoint = 7) {
    # dgm_type = 6
    
    # n: sample size
    
    var_names <- c("robs", "z", "w", "l", "probS", "s", "x", "delta", "r", "t", "c")
    
    dta <- data.frame(matrix(NA, nrow = n, ncol = length(var_names)))
    names(dta) <- var_names
    dta <- as_tibble(dta)
    
    # Baseline covariate
    dta$z <- runif(n, min = -2, max = 2)
    
    # True dropout status
    dta$r <- rbinom(n, 1, (dta$z + 3) / 6)
    
    # lognormal distributed survival time
    beta <- 1
    sigma <- 0.5
    dta$t <- rtlnorm(n, meanlog = dta$z * beta, sdlog = sigma, truncpoint = truncpoint)
    
    # lognormal distributed censoring time
    beta <- 1
    sigma <- 0.5
    dta$c <- rtlnorm(n, meanlog = dta$z * beta, sdlog = sigma, truncpoint = truncpoint)
    
    # Survival information (not affected by dropout yet)
    dta$x <- pmin(dta$t, dta$c)
    dta$delta <- as.numeric(dta$t <= dta$c)
    
    # Uniform dropout time
    dta$l <- runif(n, min = dta$t / 3, max = dta$t)
    dta$l[dta$r == 1] <- NA
    dta$l[dta$l > dta$c] <- NA
    
    # Observed dropout status
    dta$robs <- as.numeric(dta$l >= dta$x)
    dta$robs[is.na(dta$robs)] <- 1
    
    # Some other covariate "w" (not used except as a place holder in the EM algorithm by An et al. (2014))
    dta$w <- runif(n, min = 0, max = 1)
    
    # Selection for double-sampling
    dta$probS <- expit(-0.5 * dta$z + 0.5 * dta$l + 0.5)
    dta$probS[is.na(dta$probS)] <- 0
    dta$s <- rbinom(n, 1, dta$probS)
    dta$s[dta$robs == 1] <- 0
    
    # (x, delta) missing for those with robs = 0 and s = 0
    dta$x[dta$robs == 0 & dta$s == 0] <- NA
    dta$delta[dta$robs == 0 & dta$s == 0] <- NA
    
    return(dta)
}

if (0) {
    # summary statistics for dgm_start_from_R_trunclognormal()
    set.seed(123)
    
    n <- 1000000
    
    dta <- dgm_start_from_R_trunclognormal(n)
    
    1 - mean(dta$robs) # 0.361514
    
    mean(dta$s[dta$robs == 0]) # 0.7335705
    
    mean(dta$delta, na.rm = TRUE) # 0.4786285
    
    
    quantile(dta$x, probs = c(0.1, 0.9), na.rm = TRUE) # 0.1468549 3.5869629
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.1, 0.9), na.rm = TRUE) # 0.6462089 0.8104520
    
    quantile(dta$probS[dta$robs == 0], probs = c(0.25, 0.75), na.rm = TRUE) # 0.6908630 0.7850612
    
    mean(dta$c <= dta$t) # 0.499805
    
    hist(dta$t)
    tau_true <- mean(dta$t >= 0.7) # 0.588411
    print(tau_true)
    
    # set.seed(123)
    # 
    # n <- 500
    # dta <- dgm_lognormal_LCT(n)
    
    # partial correlation conditional on Robs and Z
    library(ppcor)
    dta_test <- dta[, c("t", "c", "z", "robs")]
    pcor(dta_test, method = "pearson") # 0.2436019
    
    # partial correlation, further conditional on L
    dta_test2 <- dta[dta$robs == 0, c("t", "c", "z", "l")]
    pcor(dta_test2, method = "pearson") # 0.1221957
    
}
