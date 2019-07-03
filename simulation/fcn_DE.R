###################################
# Deductive estimation of survival probability in double sampling design
# Wrapper function to obtain the estimator
#      Tianchen Qian
#       created: 2019.06.13
###################################


DE.Cox <- function(...) {
    # try a series of expanding intervals for uniroot() to find alpha that makes sum(Eif) = 0, using Coxph working model
    # see DE.wrapper() for the arguments
    
    uniroot_intervals <- rbind(c(1.5, 3), c(1.1, 5), c(-1, 5), c(-3, 5), c(-6, 6))
    return(DE.wrapper(xdelta_working_model = "coxph", uniroot_intervals = uniroot_intervals, ...))
}

DE.LN <- function(...) {
    # try a series of expanding intervals for uniroot() to find alpha that makes sum(Eif) = 0, using lognormal regression working model
    # see DE.wrapper() for the arguments
    
    uniroot_intervals <- rbind(c(-0.5, 0.5), c(-1, 1), c(-2, 2), c(-4, 4))
    return(DE.wrapper(xdelta_working_model = "lognormal", uniroot_intervals = uniroot_intervals, ...))
}



DE.wrapper <- function(dta, # unique x
                       info_0, # baseline covariates
                       info_loss, # longitudinal measurement (including dropout time)
                       models = NULL,
                       time_question, # t in P(T > t)
                       eps = 1e-4,
                       parallel_within_R = TRUE,
                       xdelta_working_model = c("coxph", "lognormal"),
                       survfit_method = c("NA", "Nelson-Aalen", "KM", "Kaplan-Meier"),
                       parallel_within_R_log_folder = NULL,
                       parallel_within_R_log_filename = NULL,
                       uniroot_intervals,
                       print_uniroot_progress = TRUE) {
    
    # A wrapper that combines sumEifSurvProbInDblSampling() and uniroot() to compute the deductive estimator.
    
    # uniroot_intervals is a matrix with 2 columns and K rows.
    # uniroot() will start from interval = uniroot_intervals[1, ],
    # and if it gets error that the two ends have the same sign,
    # it will use interval = uniroot_intervals[2, ], and so on.
    # Until it uses interval = uniroot_intervals[K, ],
    # and if an error is still observed, then it will return NA.
    
    # For all other arguments, see the documentation of sumEifSurvProbInDblSampling()

    if (class(uniroot_intervals) == "numeric") {
        uniroot_intervals <- as.matrix(uniroot_intervals, nrow = 1)
    }
    stopifnot(class(uniroot_intervals) == "matrix")
    
    output <- list(survProb = NA, se = NA, sumEif = NA, tuning_param_value = NA, ci = c(NA, NA))
    for (k in 1:nrow(uniroot_intervals)) {
        root <- try(
            {
                uniroot(function(alpha) {
                    tmp <- sumEifSurvProbInDblSampling(dta, info_0 = info_0, info_loss = info_loss, models = models, time_question = t0, tuning_param_value = alpha,
                                                       eps = eps, parallel_within_R = parallel_within_R, 
                                                       xdelta_working_model = xdelta_working_model, survfit_method = survfit_method)
                    cat("alpha =", alpha, "; sumEIF =", tmp$sumEif, "; survProb =", tmp$survProb_unperturbed, "; se =", tmp$se, "\n")
                    return(tmp$sumEif)
                },
                interval = uniroot_intervals[k, ],
                extendInt = "no")
            })
        if (class(root) != "try-error") {
            result <- sumEifSurvProbInDblSampling(dta, info_0 = info_0, info_loss = info_loss, models = models, time_question = t0, tuning_param_value = root$root,
                                                  eps = eps, parallel_within_R = parallel_within_R, 
                                                  xdelta_working_model = xdelta_working_model, survfit_method = survfit_method)
            output$survProb <- result$survProb_unperturbed
            output$se <- result$se
            output$sumEif <- result$sumEif
            output$tuning_param_value <- result$tuning_param_value
            output$ci <- c(output$survProb - 1.96 * output$se, output$survProb + 1.96 * output$se)
            break
        }
    }
    
    return(output)
    
}


DE.Cox_fixAlpha <- function(alpha, ...) {
    # compute the deductive estimator with fixed alpha value, using Coxph working model
    
    return(DE.wrapper_fixAlpha(alpha = alpha, xdelta_working_model = "coxph", ...))
}

DE.LN_fixAlpha <- function(alpha, ...) {
    # compute the deductive estimator with fixed alpha value, using lognormal regression working model
    
    return(DE.wrapper_fixAlpha(alpha = alpha, xdelta_working_model = "lognormal", ...))
}

DE.wrapper_fixAlpha <- function(alpha, ...) {
    # compute the deductive estimator with fixed alpha value.
    # See the documentation of sumEifSurvProbInDblSampling() for other arguments.
    
    output <- list(survProb = NA, se = NA, sumEif = NA, tuning_param_value = NA, ci = c(NA, NA))
    result <- try(sumEifSurvProbInDblSampling(tuning_param_value = alpha, ...))
    if (class(result) != "try-error") {
        output$survProb <- result$survProb_unperturbed
        output$se <- result$se
        output$sumEif <- result$sumEif
        output$tuning_param_value <- result$tuning_param_value
        output$ci <- c(output$survProb - 1.96 * output$se, output$survProb + 1.96 * output$se)
    }
    return(output)
}