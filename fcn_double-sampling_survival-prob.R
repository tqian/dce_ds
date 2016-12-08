###################################
# Deductive estimation of survival probability in double sampling design
# Source code of functions
#      Tianchen Qian
#        2016.09.22
###################################

# In variable names, x.y means x | y ("." means conditional)

# Possible issues:
# 1) Code breaks down when the true estimand is small (e.g. ~ 0.5).
#    The issue is in getSurvProb.info_robsFromPxdelta.info_robs(), where dLambda has too much NA's,
#    and hence this line
#        time_to_return <- max(c(1:n_s1)[x_value <= time_question])
#    returns NA.
#
# 2) Sometimes coxph doesn't converge.
#    Don't know the exact reason for this, but probably too little variation in the survival data.


sumEifSurvProbInDblSampling <- function(dta, # unique x
                                        info_0 = c("age10", "prxcd4v"),
                                        info_loss = c("l"),
                                        time_question = 365,
                                        eps = 1e-4,
                                        tuning_param_value = 1,
                                        parallel_within_R = TRUE,
                                        max_cores = 10){
    
    ### check ###
    # unique x values
    stopifnot(length(na.omit(dta$x)) == length(unique(na.omit(dta$x))))
    
    ### sort dta ###
    # (1) first robs = 0, then robs = 1
    # (2) among each robs, x is increasing
    # (3) among robs = 0, first s = 1, then s = 0 (x = NA)
    dta <- dta[order(dta$robs, dta$x, decreasing = FALSE), ]
    
    # for internal tracking
    dta$.rowid <- 1:nrow(dta) 
    
    masAndProb_extended <- getExtendedMasAndProbFromData(dta, info_0, info_loss, tuning_param_value)
    estimand_unperturbed <- getSurvProbFromMasAndProb(masAndProb_extended, time_question)

    ### check ###
    # in masAndProb_extended, rowid is unchanged
    n_robs0 <- sum(dta$robs == 0)
    n_obs <- nrow(dta)
    stopifnot(identical(masAndProb_extended$robs0$data$.rowid, 1:n_robs0))
    stopifnot(identical(masAndProb_extended$robs1$data$.rowid, (n_robs0 + 1):n_obs))

    ### perturb to get sumEif ###
    if (parallel_within_R) {
        
        # cl <- makeSOCKcluster(min(detectCores() - 1, max_cores))
        registerDoMC(min(detectCores() - 1, max_cores))
        # registerDoSNOW(cl)
        
        # progress bar using tcltk
        # library(tcltk)
        # pb <- tkProgressBar(max = n_obs)
        # progress <- function(n) setTkProgressBar(pb, n)
        # browser()
        # progress <- function(n) cat(n, '')
        # opts <- list(progress = progress) # pass .options.snow = opts into foreach
        
        tmp <- foreach(iobs = 1:n_obs, .combine = rbind) %dopar% {
            message(iobs)
            masAndProb_perturbed <- perturbPmas(masAndProb_extended, iobs, eps)
            estimand_perturbed_i <- getSurvProbFromMasAndProb(masAndProb_perturbed, time_question)
            return(c(iobs, estimand_perturbed_i))
        }
        # stopCluster(cl)
        estimand_perturbed <- tmp[order(tmp[, 1]), 2]
        
    } else {
        estimand_perturbed <- numeric(n_obs)
        for (iobs in 1:n_obs) {
            cat(iobs, '')
            masAndProb_perturbed <- perturbPmas(masAndProb_extended, iobs, eps)
            estimand_perturbed[iobs] <- getSurvProbFromMasAndProb(masAndProb_perturbed, time_question)
        }
    }
    
    dta$survProb_perturbed <- estimand_perturbed
    dta$eif <- (dta$survProb_perturbed - estimand_unperturbed) / eps
    sumEif <- sum(dta$eif)
    
    return(list(survProb_unperturbed = estimand_unperturbed,
                sumEif = sumEif,
                data = dta,
                tuning_param_value = tuning_param_value,
                eps = eps))
}


getSurvProbFromMasAndProb <- function(masAndProb,
                                      time_question
) {
    ##########
    # Return prob(T > time_question).
    #
    # Keyword Arguments:
    # masAndProb -- a list of two lists: robs0 and robs1
    #
    # Method:
    ##########
    
    ### 1) decompose into Probs, Pinfo_robs, Ps.info_robs, Pxdelta.s_info_robs ###
    
    Pmas_robs0 <- masAndProb$robs0$prob
    n_s1.robs0 <- masAndProb$robs0$mas$n_s1
    
    Probs0 <- sum(Pmas_robs0)
    Pinfo_robs0 <- rowSums(Pmas_robs0)
    Ps.info_robs0 <- 1 - Pmas_robs0[, ncol(Pmas_robs0)] / Pinfo_robs0   # Pmas$robs0[, ncol(Pmas$robs0)] = prob(s=0 | info, robs=0)
    Pxdelta.s_info_robs0 <- Pmas_robs0[, 1:n_s1.robs0] / (Ps.info_robs0 * Pinfo_robs0)
    
    Pmas_robs1 <- masAndProb$robs1$prob
    n_s1.robs1 <- masAndProb$robs1$mas$n_s1
    
    Probs1 <- sum(Pmas_robs1)
    Pinfo_robs1 <- rowSums(Pmas_robs1)
    Pxdelta.info_robs1 <- Pmas_robs1 / Pinfo_robs1
    
    ### 2) calculate survProb.info_robs0 and survProb.info_robs1 ###
    
    survProb.info_robs0 <- getSurvProb.info_robsFromPxdelta.info_robs(
        Pxdelta.s_info_robs0, masAndProb$robs0$mas$x_value[1:n_s1.robs0], masAndProb$robs0$mas$delta_value[1:n_s1.robs0], time_question)
    
    survProb.info_robs1 <- getSurvProb.info_robsFromPxdelta.info_robs(
        Pxdelta.info_robs1, masAndProb$robs1$mas$x_value[1:n_s1.robs1], masAndProb$robs1$mas$delta_value[1:n_s1.robs1], time_question)
    
    ### 3) survProb ###
    
    survProb <- sum(c(survProb.info_robs0 * Pinfo_robs0, survProb.info_robs1 * Pinfo_robs1))
    return(survProb)
    
}

getSurvProb.info_robsFromPxdelta.info_robs <- function(Pxdelta.info_robs,
                                                       x_value,
                                                       delta_value,
                                                       time_question
) {
    ### get Nelson-Aalen estimator for cumulative hazard ###

    n_obs <- nrow(Pxdelta.info_robs)
    n_s1 <- length(x_value)
    
    stopifnot(ncol(Pxdelta.info_robs) == n_s1)
    
    Y <- rowCumsum(Pxdelta.info_robs[, n_s1:1])[, n_s1:1]
    # Y can be 0 due to numerical issues and the fact that Y can be very close to zero
    dN <- Pxdelta.info_robs * rep.row(delta_value, n_obs)
    dLambda <- dN / Y # cumulative hazard
    
    # 20160920: debug
    # dLambda[is.na(dLambda)] <- 0 # NA is from 0/0 in dN/Y
    # if (any(is.na(dLambda))) {
    #     browser()
    # }
    
    ### get survival probability from cumulative hazard ###
    
    S <- exp(- rowCumsum(dLambda))
    time_to_return <- max(c(1:n_s1)[x_value <= time_question])
    return(S[, time_to_return])
}


perturbPmas <- function(masAndProb, iobs, eps) {
    
    n_robs0 <- nrow(masAndProb$robs0$prob)
    n_s1 <- masAndProb$robs0$mas$n_s1
    
    masAndProb$robs0$prob <- (1 - eps) * masAndProb$robs0$prob
    masAndProb$robs1$prob <- (1 - eps) * masAndProb$robs1$prob
    
    if (iobs <= n_s1) {
        # being perturbed: robs = 0, s = 1
        masAndProb$robs0$prob[iobs, iobs] <- masAndProb$robs0$prob[iobs, iobs] + eps
    } else if (iobs <= n_robs0) {
        # being perturbed: robs = 0, s = 0
        masAndProb$robs0$prob[iobs, n_s1 + 1] <- masAndProb$robs0$prob[iobs, n_s1 + 1] + eps
    } else {
        # being perturbed: robs = 1
        masAndProb$robs1$prob[iobs - n_robs0, iobs - n_robs0] <- masAndProb$robs1$prob[iobs - n_robs0, iobs - n_robs0] + eps
    }
    return(masAndProb)
}



getExtendedMasAndProbFromData <- function(dta, # unique x
                                    info_0,
                                    info_loss,
                                    tuning_param_value) {
    ##########
    # Return a list of two matrices.
    # Matrix robs0 stands for Pmas.and.robs0, with nrow = #{unique info_0*info_loss: robs=0}, ncol = #{unique x: robs=0} + 1 (1 for NA).
    # Matrix robs1 stands for Pmas.and.robs1, with nrow = #{unique info_0: robs=1}, ncol = #{unique x: robs=1}.
    #
    # Keyword Arguments:
    # dta         -- entire data set
    # info_0      -- covariate names of info_0
    # info_loss   -- covariate names of info_loss
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs) for both robs0 and robs1
    #
    # Method:
    #   See getExtendedPxdelta_s_info.robs0() and getExtendedPxdelta_info.robs1().
    ##########
    
    dta.robs0 <- subset(dta, robs == 0)
    mas_and_prob.robs0 <- getExtendedPxdelta_s_info.robs0(dta.robs0, info_0, info_loss, tuning_param_value)
    Probs0 <- nrow(dta.robs0) / nrow(dta)
    
    
    dta.robs1 <- subset(dta, robs == 1)
    mas_and_prob.robs1 <- getExtendedPxdelta_info.robs1(dta.robs1, info_0, tuning_param_value)
    Probs1 <- nrow(dta.robs1) / nrow(dta)
    
    return(list(robs0 = list(data = mas_and_prob.robs0$data,
                             mas = mas_and_prob.robs0$mas, 
                             prob = mas_and_prob.robs0$prob * Probs0),
                robs1 = list(data = mas_and_prob.robs1$data,
                             mas = mas_and_prob.robs1$mas, 
                             prob = mas_and_prob.robs1$prob * Probs1)))
}


getExtendedPxdelta_s_info.robs0 <- function(dta.robs0, # unique x, unique (info_0, info_loss)
                                            info_0,
                                            info_loss,
                                            tuning_param_value
) {
    ##########
    # Return a matrix of prob((x, delta) * s, s, info | robs = 0).
    # nrow = #{unique info_0*info_loss: robs=0}, ncol = #{unique x: robs=0} + 1 (1 for NA).
    #
    # Keyword Arguments:
    # dta.robs0   -- data set with robs=0
    # info_0      -- covariate names of info_0
    # info_loss   -- covariate names of info_loss
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    #
    # Method:
    #   1) Get extended Pxdelta.info_robs0: prob(x, delta | info, robs = 0)
    #   2) Get Ps.info_robs0: prob(s = 1 | info, robs = 0)
    #   3) Get Pinfo.robs0: prob(info | robs = 0)
    #   4) Multiply 1) - 3) to get prob((x, delta) * s, s, info | robs = 0)
    ##########
    
    ### order dta.robs0: first s=1, then s=0. Among s=1, x is increasing. ###
    
    # dta.robs0 <- dta.robs0[order(dta.robs0$s, decreasing = TRUE), ]
    # dta.robs0 <- dta.robs0[order(dta.robs0$x, decreasing = FALSE), ]
    stopifnot(!is.unsorted(rev(dta.robs0$s))) # require s to be decreasing
    stopifnot(!is.unsorted(dta.robs0$x, na.rm = TRUE)) # require x to be increasing
    
    ### create mas ###
    
    x_values <- na.omit(dta.robs0$x)
    delta_values <- na.omit(dta.robs0$delta)
    
    mas <- list(x_value = c(x_values, NA),
                delta_value = c(delta_values, NA),
                n_s1 = length(x_values))
    
    ### 1) create extended Pxdelta.info_robs0 ###
    ### prob(x, delta | info, robs = 0) ###
    
    Pxdelta.info_robs0 <- getPxdelta.info_robs(dta.robs0, c(info_0, info_loss), mas, robs_type = 0)
    # Pxdelta.info_robs0 <- extendPxdelta.info_robs(Pxdelta.info_robs0, tuning_param_value)
    Pxdelta.info_robs0 <- extendPxdelta.info_robs_ver2(Pxdelta.info_robs0, x_values, max(dta.robs0$c), tuning_param_value)
    
    ### 2) create Ps.info_robs0 ###
    ### prob(s = 1 | info, robs = 0) ###
    
    Ps.info_robs0 <- getPs.info(dta.robs0, c(info_0, info_loss))
    
    ### 3) create Pinfo.robs0 ###
    ### prob(info | robs = 0) ###
    
    Pinfo.robs0 <- rep(1 / nrow(dta.robs0), nrow(dta.robs0))
    
    ### 4) combine the above three ###
    ### prob((x, delta) * s, s, info | robs = 0) ###
    
    Pxdelta_s_info.robs0 <- cbind(Pxdelta.info_robs0 * Ps.info_robs0, 1 - Ps.info_robs0) * Pinfo.robs0
    
    ### 5) remove (x_max + 1, 0) and (x_max + 2, 1) from mas
    return(list(data = dta.robs0, mas = mas, prob = Pxdelta_s_info.robs0))
}


getExtendedPxdelta_info.robs1 <- function(dta.robs1, # unique x, unique info_0
                                          info_0,
                                          tuning_param_value
) {
    ##########
    # Return a matrix of prob(x, delta, info | robs = 1).
    # nrow = #{unique info_0: robs=1}, ncol = #{unique x: robs=1}.
    #
    # Keyword Arguments:
    # dta.robs1   -- data set with robs=1
    # info_0      -- covariate names of info_0
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    #
    # Method:
    #   1) Get extended Pxdelta.info_robs1: prob(x, delta | info, robs = 1)
    #   2) Get Pinfo.robs1: prob(info | robs = 1)
    #   3) Multiply 1) - 2) to get prob(x, delta, info | robs = 1)
    ##########
    
    ### order dta.robs1: x is increasing. ###
    
    # dta.robs1 <- dta.robs1[order(dta.robs1$x, decreasing = FALSE), ]
    stopifnot(!is.unsorted(dta.robs1$x, na.rm = TRUE)) # require x to be increasing
    
    ### create mas ###
    
    x_values <- dta.robs1$x
    delta_values <- dta.robs1$delta
    
    mas <- list(x_value = x_values,
                delta_value = delta_values,
                n_s1 = length(x_values))
    
    ### 1) create extended Pxdelta.info_robs1 ###
    ### prob(x, delta | info, robs = 1) ###
    
    Pxdelta.info_robs1 <- getPxdelta.info_robs(dta.robs1, info_0, mas, robs_type = 1)
    # Pxdelta.info_robs1 <- extendPxdelta.info_robs(Pxdelta.info_robs1, tuning_param_value)
    Pxdelta.info_robs1 <- extendPxdelta.info_robs_ver2(Pxdelta.info_robs1, x_values, max(dta.robs1$c), tuning_param_value)
    
    ### 2) create Pinfo.robs1 ###
    ### prob(info | robs = 1) ###
    
    Pinfo.robs1 <- rep(1 / nrow(dta.robs1), nrow(dta.robs1))
    
    ### 3) combine the above two ###
    ### prob(x, delta, info | robs = 1) ###
    
    Pxdelta_info.robs1 <- Pxdelta.info_robs1 * Pinfo.robs1
    
    return(list(data = dta.robs1, mas = mas, prob = Pxdelta_info.robs1))
}


getPxdelta.info_robs <- function(dta.robs,
                                 covariates,
                                 mas,
                                 robs_type = c(0, 1)
                                 ) {
    ##########
    # Return a matrix of prob(x, delta | info, robs) for the corresponding robs,
    # with nrow = # unique(info), ncol = n_s1
    #
    # Keyword Arguments:
    # dta.robs    -- data set for the corresponding robs
    # covariates  -- covariates to be used in the cox model, same for T and C
    # mas         -- mas for the corresponding robs
    #
    # Method:
    #   1) Fit Coxph for T and C, then get probability mass prob(T = x | info, robs)
    #   and prob(C = x | info, robs) on each grid (x | info, robs).
    #   2) Calculate uncalibrated prob(x, delta | info, robs) using probs in 1).
    #   3) Calibrate using 2015.12.20_calibration of discrete probability on observed data.pdf.
    ##########
    
    n_robs <- nrow(dta.robs)
    n_s1 <- mas$n_s1
    
    if (robs_type == 1) {
        dta.robs$s <- 1
    }
    dta.robs_s1 <- subset(dta.robs, s == 1)
    dta.robs_s1_deltaflipped <- dta.robs_s1
    dta.robs_s1_deltaflipped$delta <- 1 - dta.robs_s1_deltaflipped$delta
    
    ### Cox model for T ###
    
    formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates, collapse = "+")))
    model_t.info_robs <- coxph(formula, data = dta.robs_s1)
    st <- surv.from.cox(dta.robs_s1, model_t.info_robs, dta.robs, covariates)[, dta.robs_s1$delta==1]
    jumps.t <- cbind(1, st) - cbind(st, 0)
    
    ### Cox model for C ###
    formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates, collapse = "+")))
    model_c.info_robs <- coxph(formula, data = dta.robs_s1_deltaflipped)
    sc <- surv.from.cox(dta.robs_s1_deltaflipped, model_c.info_robs, dta.robs, covariates)[, dta.robs_s1_deltaflipped$delta==1]
    jumps.c <- cbind(1, sc) - cbind(sc, 0)
    
    ### Get uncalibrated prob(x, delta | info, robs) ###
    
    x_max <- max(mas$x_value, na.rm = TRUE)
    x_forcalc <- c(mas$x_value[1:n_s1], x_max + 1, x_max + 2)
    delta_forcalc <- c(mas$delta_value[1:n_s1], 0, 1)
    
    # Pt: prob(T = x | info, robs)
    # Pc: prob(C = x | info, robs)
    Pt <- matrix(0, nrow = n_robs, ncol = length(x_forcalc))
    t_location_in_x_forcalc <- which(delta_forcalc == 1)
    Pc <- Pt # except for the NA from s=0
    c_location_in_x_forcalc <- which(delta_forcalc == 0)
    for(i in 1:nrow(dta.robs)){
        Pt[i, t_location_in_x_forcalc] <- jumps.t[i, ]
        Pc[i, c_location_in_x_forcalc] <- jumps.c[i, ]
    }
    
    # St: prob(T > x | info, robs)
    # Sc: prob(C > x | info, robs)
    St <- 1 - rowCumsum(Pt)
    Sc <- 1 - rowCumsum(Pc)
    
    # calculation of prob(x, delta | info, robs)
    Pxdelta <- Pt * Sc * rep.row(delta_forcalc, n_robs) + Pc * St * (1 - rep.row(delta_forcalc, n_robs))
    
    # removing prob on the added (xmax+1, 0) and (xmax+2, 1)
    Pxdelta <- normalizeByRow(Pxdelta[, 1:n_s1])
    
    return(normalizeByRow(Pxdelta))
    
    ### Calibrate prob(x, delta | info, robs) ###
    # see 2015.12.20_calibration of discrete probability on observed data.pdf
    # not needed for this estimand, since we are already working on CDF not PDF
    
    # source.cov <- c(rep(1, n_s1), rep(0, n_robs - n_s1))
    # denominator <- as.vector(t(source.cov) %*% Pxdelta)
    # qhat_ij <- Pxdelta / rep.row(denominator, nrow(Pxdelta))
    # 
    # return(normalizeByRow(qhat_ij))
}

extendPxdelta.info_robs <- function(Pxdelta.info_robs, tuning_param_value) {
    factor <- seq(from = 1, to = tuning_param_value, length.out = ncol(Pxdelta.info_robs))
    Pxdelta.info_robs <- Pxdelta.info_robs * rep.row(factor, nrow(Pxdelta.info_robs))
    return(normalizeByRow(Pxdelta.info_robs))
}

extendPxdelta.info_robs_ver2 <- function(Pxdelta.info_robs, x_value_withoutNA, cmax, alpha) {
    factor <- 1 + (alpha - 1) * x_value_withoutNA / cmax
    if (length(factor) < ncol(Pxdelta.info_robs)) {
        factor <- c(factor, 1)
    }
    stopifnot(length(factor) == ncol(Pxdelta.info_robs))
    Pxdelta.info_robs <- Pxdelta.info_robs * rep.row(factor, nrow(Pxdelta.info_robs))
    return(normalizeByRow(Pxdelta.info_robs))
}

getPs.info <- function(dta.robs0, covariates) {
    ##########
    # Return a vector of estimated propensity score: prob(s=1 | covariates, robs=0)
    #
    # Keyword Arguments:
    # dta.robs0   -- data set with robs=0
    # covariates  -- covariates to be used in propensity score model
    ##########
    
    formula <- as.formula(paste0("s ~ ", paste(covariates, collapse = "+")))
    model_s <- glm(formula, data = dta.robs0, family = binomial(link = logit))
    return(predict(model_s, data = dta.robs0, type = "response"))
}



# Utility functions -------------------------------------------------------

cumsum_from_last <- function(x) {
    n <- length(x)
    return(cumsum(x[n:1])[n:1])
}

rowCumsum <- function(mat) {
    return(t(apply(mat, 1, cumsum)))
    
    # Example:
    # > x
    # [,1] [,2] [,3] [,4]
    # [1,]    1    4    7   10
    # [2,]    2    5    8   11
    # [3,]    3    6    9   12
    # 
    # > t(apply(x, 1, cumsum))
    # [,1] [,2] [,3] [,4]
    # [1,]    1    5   12   22
    # [2,]    2    7   15   26
    # [3,]    3    9   18   30
}

normalizeByRow <- function(mat) {
    return(mat * (1 / rowSums(mat)))
    # > zz
    # [,1] [,2] [,3] [,4]
    # [1,]    1    4    7   10
    # [2,]    2    5    8   11
    # [3,]    3    6    9   12
    # > normalizeByRow(zz)
    # [,1]      [,2]      [,3]      [,4]
    # [1,] 0.04545455 0.1818182 0.3181818 0.4545455
    # [2,] 0.07692308 0.1923077 0.3076923 0.4230769
    # [3,] 0.10000000 0.2000000 0.3000000 0.4000000
}

rep.row <- function(x, n) {
    return(matrix(rep(x, each = n), nrow = n))
    # > rep.row(1:3,5)
    # [,1] [,2] [,3]
    # [1,]    1    2    3
    # [2,]    1    2    3
    # [3,]    1    2    3
    # [4,]    1    2    3
    # [5,]    1    2    3
}
rep.col <- function(x, n) {
    return(matrix(rep(x, each = n), ncol = n, byrow = TRUE))
    # > rep.col(1:3,5)
    # [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    2    2    3
    # [2,]    1    1    2    3    3
    # [3,]    1    2    2    3    3
}



#--- begin: needs simplification ---
# Breslow's estimator for baseline hazard
# see p.333 eq. (12.3.2) of Lee - statistical methods for survival data analysis
surv.from.cox <- function(datain, model, datanew, covname){
    ## uses cox model on datain to predict survival on datanew
    ## checked nov 7 2015 (assume data are sorted by X)
    ndata <- nrow(datain) 
    Y <- exp(as.matrix(datain[, covname]) %*% model$coef)
    Y <- cumsum_from_last(Y)
    Lambda <- cumsum(datain$delta / Y) # vector of length ndata, cumulative baseline hazard function
    expfactor <- exp(as.matrix(datanew[, covname]) %*% model$coef)
    surv <- exp(-expfactor %*% t(Lambda))
    # nrow(surv) = nrow(datanew); ncol(surv) = length(unique(datain$x))
    stopifnot(all(surv <= 1))
    return(surv)
    
    # Note: 
    # If error message:
    # """
    # Error: all(surv <= 1) is not TRUE
    # In addition: Warning message:
    #     In fitter(X, Y, strats, offset, init, control, weights = weights,  :
    #                   Ran out of iterations and did not converge
    # """
    # Then it is likely that the coxph routine in parent function did not converge.
}
#--- end: needs simplification ---