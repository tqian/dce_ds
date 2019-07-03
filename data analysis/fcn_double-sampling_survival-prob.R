###################################
# Deductive estimation of survival probability in double sampling design
# Source code of functions
#      Tianchen Qian
#       created: 2016.09.22
#       most recent update: 2019.06.05
###################################

# In variable names, x.y means x | y ("." means conditional), and _ means comma.
# For example, Pxdelta.s_info_robs means P(x, delta | s, info, robs).

# Updates:
#   2017.04.08: add small positive mass to Pxdelta so that no cell has prob 0 (this is removed at the 2019.06.05 update)
#   2019.06.05: revise the infrastructure so that the covariates and the failure time can have duplicate values


library(survival)
library(dplyr)
library(doMC) # for parallel_within_R in sumEifSurvProbInDblSampling()


sumEifSurvProbInDblSampling <- function(dta, # unique x
                                        info_0, # baseline covariates
                                        info_loss, # longitudinal measurement (including dropout time)
                                        models = NULL,
                                        time_question, # t in P(T > t)
                                        eps = 1e-4,
                                        tuning_param_value = 1,
                                        parallel_within_R = TRUE,
                                        xdelta_working_model = c("coxph", "lognormal"),
                                        survfit_method = c("KM", "Kaplan-Meier", "NA", "Nelson-Aalen"), # KM seems to be correct, while NA seems to be giving wrong answers?
                                        parallel_within_R_log_folder = NULL,
                                        parallel_within_R_log_filename = NULL){
    ##########
    # Return a list of:
    # survProb_unperturbed  -- estimand at tuning parameter (i.e., the deductive estimator)
    # se                    -- standard error of the estimator
    # sumEif                -- value of sum(eif) at tuning parameter (i.e., alpha in the paper)
    # data                  -- a copy of the original data
    # tuning_param_value    -- value of tuning parameter (i.e., alpha in the paper)
    # eps                   -- epsilon used in numeric Gateaux derivative
    # xdelta_working_model  -- the working model used for constructing P(x, delta | s, info, robs), one of c("coxph", "lognormal")
    # survfit_method        -- two methods to calculate the survival probability P(T > t) given P(x, delta), KM and NA.
    #                          They should be equivalent, but currently it seems that only KM gives correct results.
    # info_0                -- covariate names of info_0 (baseline variables)
    # info_loss             -- covariate names of info_loss (longitudinal measurements up to loss-to-followup)
    #
    #
    # Keyword Arguments:
    # dta                   -- data set
    # info_0                -- covariate names of info_0 (baseline variables)
    # info_loss             -- covariate names of info_loss (longitudinal measurements up to loss-to-followup)
    # models                -- a list of model_t.robs0, model_c.robs0, model_s.robs0, model_t.robs1, model_c.robs1.
    #                          if models = NULL, use info_0 and info_loss to construct all models
    # time_question         -- t in the estimand P(T > t)
    # eps                   -- epsilon used in numeric Gateaux derivative
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    # parallel_within_R     -- whether to parallel within R (using doMC package); need to require for multiple cores each job if running on cluster.
    # xdelta_working_model  -- the working model used for constructing P(x, delta | s, info, robs), one of c("coxph", "lognormal")
    # survfit_method        -- two methods to calculate the survival probability P(T > t) given P(x, delta), KM and NA.
    #                          They should be equivalent, but currently it seems that only KM gives correct results.
    # parallel_within_R_log_folder, parallel_within_R_log_filename -- folder and filename to print the parallel progress
    #
    #
    # Method:
    # See paper.
    ##########
    
    # if (0) {
    #     # debugging for sumEifSurvProbInDblSampling()
    #     set.seed(123)
    #     n <- 50
    #     dta <- dgm_start_from_R(n)
    #     info_0 <- "z"
    #     info_loss <- c("w", "l")
    #     models <- NULL
    #     time_question <- 0.7
    #     eps <- 1e-4
    #     tuning_param_value <- 1
    #     parallel_within_R <- FALSE
    #     xdelta_working_model <- "coxph"
    #     survfit_method <- "NA"
    #     
    #     # create duplicate rows
    #     dta[1:2, ] <- dta[1,] # robs=1, delta=0
    #     dta[3:4, ] <- dta[3,] # robs=1, delta=1
    #     dta[6:7, ] <- dta[6,] # robs=0, s=1, delta=1
    #     dta[13:14, ] <- dta[13,] # robs=0, s=1, delta=0
    #     dta[10:11, ] <- dta[10,] #robs=0, s=0, x=NA, delta=NA
    # 
    #     dta$tmp_rowid <- 1:n
    # 
    #     dta$tmp_twin_rowid <- NA
    #     dta$tmp_twin_rowid[1] <- 2
    #     dta$tmp_twin_rowid[2] <- 1
    #     dta$tmp_twin_rowid[3] <- 4
    #     dta$tmp_twin_rowid[4] <- 3
    #     dta$tmp_twin_rowid[6] <- 7
    #     dta$tmp_twin_rowid[7] <- 6
    #     dta$tmp_twin_rowid[13] <- 14
    #     dta$tmp_twin_rowid[14] <- 13
    #     dta$tmp_twin_rowid[10] <- 11
    #     dta$tmp_twin_rowid[11] <- 10
    #     
    #     print(dta, n = 20)
    # }
    # 
    # if (0) {
    #     # debugging for sumEifSurvProbInDblSampling()
    #     set.seed(123)
    #     n <- 200
    #     dta <- dgm_discrete_var(n)
    #     info_0 <- "z"
    #     info_loss <- c("l")
    #     models <- NULL
    #     time_question <- 3.2
    #     eps <- 1e-8
    #     tuning_param_value <- 1
    #     parallel_within_R <- FALSE
    #     xdelta_working_model <- "coxph"
    #     survfit_method <- "NA"
    #     
    #     print(dta, n = 20)
    #     print(dta, n = Inf)
    # }
    
    dta <- as_tibble(dta)
    
    # type of method used in calculating the survival probability given info_0
    survfit_method <- match.arg(survfit_method)

    # type of working model used in constructing P(x, delta | s, info, robs)
    xdelta_working_model <- match.arg(xdelta_working_model)    
    
    # models (regressors) used in modeling T, C, S, for robs=0 and for robs=1
    
    if (is.null(models)) {
        models <- list(model_t.robs0 = c(info_0, info_loss),
                       model_c.robs0 = c(info_0, info_loss),
                       model_s.robs0 = c(info_0, info_loss),
                       model_t.robs1 = c(info_0, info_loss),
                       model_c.robs1 = c(info_0, info_loss))
        # If info_loss contains dropout time (which will be NA for robs1 patients),
        # the model fitting function will check and remove those covariates before proceeding.
    }
    stopifnot(names(models) == c("model_t.robs0", "model_c.robs0", "model_s.robs0", "model_t.robs1", "model_c.robs1"))
    
    
    ### sort dta ###
    # (1) first robs = 0, then robs = 1
    # (2) among each robs, x is increasing
    # (3) among robs = 0, first s = 1, then s = 0 (x = NA)
    dta <- dta[order(dta$robs, dta$x, decreasing = FALSE), ]
    
    # for internal tracking
    dta$row_id <- 1:nrow(dta) 
    
    
    ### make hash id's for the unique values of (info_0) ###
    # Independent censoring T \perp C only holds conditional on info_0
    dta <- dta[do.call(order, dta[, info_0]), ]
    dta$info0_hashid <- as.integer(NA) # integer hashid
    current_hashid <- 0L
    for (irow in 1:nrow(dta)) {
        if (irow == 1 | !identical(dta[irow, info_0], dta[irow - 1, info_0]) ) {
            current_hashid <- current_hashid + 1L
        }
        dta$info0_hashid[irow] <- current_hashid
    }
    dta$info0_hashid <- as.integer(dta$info0_hashid)
    
    
    ### make hash id's for the unique values of (info_0, info_loss) ###
    info_varname <- unique(do.call(c, unname(models)))
    dta <- dta[do.call(order, dta[, info_varname]), ]
    dta$info_hashid <- as.integer(NA) # integer hashid
    current_hashid <- 0L
    for (irow in 1:nrow(dta)) {
        if (irow == 1 | !identical(dta[irow, info_varname], dta[irow - 1, info_varname]) ) {
            current_hashid <- current_hashid + 1L
        }
        dta$info_hashid[irow] <- current_hashid
    }
    dta$info_hashid <- as.integer(dta$info_hashid)
    
    ### make hash id's for the unique values of (x, delta) ###
    dta <- dta[do.call(order, dta[, c("x", "delta")]), ]
    dta$xdelta_hashid <- as.integer(NA) # integer hashid
    current_hashid <- 0L
    for (irow in 1:nrow(dta)) {
        if (irow == 1 | !identical(dta[irow, c("x", "delta")], dta[irow - 1, c("x", "delta")]) ) {
            current_hashid <- current_hashid + 1L
        }
        dta$xdelta_hashid[irow] <- current_hashid
    }
    dta$xdelta_hashid <- as.integer(dta$xdelta_hashid)
    
    ### sort dta ###
    # (1) first robs = 0, then robs = 1
    # (2) among each robs, x is increasing
    # (3) among robs = 0, first s = 1, then s = 0 (x = NA)
    dta <- dta[order(dta$robs, dta$x, decreasing = FALSE), ]
    
    
    masAndProb_extended <- getExtendedMasAndProbFromData(dta, info_varname, models, tuning_param_value, xdelta_working_model)
    estimand_unperturbed <- getSurvProbFromMasAndProb(dta, masAndProb_extended, time_question, survfit_method)

    ### perturb to get sumEif ###
    if (parallel_within_R) {
        
        max_cores <- 64
        n_using_cores <- max(min(detectCores() - 1, max_cores), 1)
        registerDoMC(n_using_cores)
        
        if (!is.null(parallel_within_R_log_filename)) {
            if (!is.null(parallel_within_R_log_folder)) {
                dir.create(parallel_within_R_log_folder, showWarnings = FALSE)
                if (substr(parallel_within_R_log_folder, nchar(parallel_within_R_log_folder)-1, nchar(parallel_within_R_log_folder)) != "/") {
                    parallel_within_R_log_folder <- paste0(parallel_within_R_log_folder, "/")
                }
                log_filename <- paste0(parallel_within_R_log_folder, parallel_within_R_log_filename)
            } else {
                log_filename <- parallel_within_R_log_filename
            }
            writeLines(c(""), log_filename)
            sink(log_filename, append=FALSE)
        }
        tmp <- foreach(iobs = 1:nrow(dta), .combine = rbind) %dopar% {
            cat(iobs, "")
            masAndProb_perturbed <- perturbPmas(masAndProb_extended, dta[iobs, ], eps)
            estimand_perturbed_i <- getSurvProbFromMasAndProb(dta, masAndProb_perturbed, time_question, survfit_method)
            return(c(iobs, estimand_perturbed_i))
        }
        if (!is.null(parallel_within_R_log_filename)) {
            sink()
        }
        estimand_perturbed <- tmp[order(tmp[, 1]), 2]
        
    } else {
        estimand_perturbed <- numeric(nrow(dta))
        for (iobs in 1:nrow(dta)) {
            # cat(iobs, '')
            masAndProb_perturbed <- perturbPmas(masAndProb_extended, dta[iobs, ], eps)
            estimand_perturbed[iobs] <- getSurvProbFromMasAndProb(dta, masAndProb_perturbed, time_question, survfit_method)
        }
    }
    
    dta$survProb_perturbed <- estimand_perturbed
    dta$eif <- (dta$survProb_perturbed - estimand_unperturbed) / eps
    sumEif <- sum(dta$eif)
    
    se <- sqrt(sum(dta$eif^2)) / nrow(dta)
    
    return(list(survProb_unperturbed = estimand_unperturbed,
                se = se,
                sumEif = sumEif,
                data = dta,
                tuning_param_value = tuning_param_value,
                eps = eps,
                xdelta_working_model = xdelta_working_model,
                survfit_method = survfit_method,
                info_0 = info_0,
                info_loss = info_loss))
}


getSurvProbFromMasAndProb <- function(dta,
                                      masAndProb,
                                      time_question,
                                      survfit_method = c("NA", "Nelson-Aalen", "KM", "Kaplan-Meier")
) {
    ##########
    # Return prob(T > time_question).
    #
    # Keyword Arguments:
    # dta                   -- data set
    # masAndProb            -- a list of two lists: ds0 and ds1; see getExtendedMasAndProbFromData() for details
    # time_question         -- t in the estimand P(T > t)
    # survfit_method        -- two methods to calculate the survival probability P(T > t) given P(x, delta), KM and NA.
    #                          They should be equivalent, but currently it seems that only KM gives correct results.
    #
    ##########
    
    ds0 <- masAndProb$ds0
    ds1 <- masAndProb$ds1
    
    
    ### 1) construct prob_table_xdeltainfo0, which is the pmf of unique(info0) * unique(xdelta) ###
    
    # Expand ds0 and ds1 so that they share the same unique(info) and unique(x,delta),
    # and the newly added points to ds0 and ds1 have 0 probability.
    
    stopifnot(!is.unsorted(ds0$row_index$info_hashid))
    stopifnot(!is.unsorted(ds0$col_index$xdelta_hashid))
    stopifnot(!is.unsorted(ds1$row_index$info_hashid))
    stopifnot(!is.unsorted(ds1$col_index$xdelta_hashid))
    
    info_hashid_combined <- sort(unique(c(ds0$row_index$info_hashid, ds1$row_index$info_hashid)))
    xdelta_hashid_combined <- sort(unique(c(ds0$col_index$xdelta_hashid, ds1$col_index$xdelta_hashid)))
    
    # info_hashid_combined and info_hashid_combined need to be consecutive integers for the following
    # construction of prob_table_expanded to work.
    stopifnot(identical(info_hashid_combined, 1:length(info_hashid_combined)))
    stopifnot(identical(xdelta_hashid_combined, 1:length(xdelta_hashid_combined)))
    
    ds0$prob_table_expanded <- matrix(0, nrow = length(info_hashid_combined), ncol = length(xdelta_hashid_combined))
    ds0$prob_table_expanded[ds0$row_index$info_hashid, ds0$col_index$xdelta_hashid] <- ds0$prob_table
    ds1$prob_table_expanded <- matrix(0, nrow = length(info_hashid_combined), ncol = length(xdelta_hashid_combined))
    ds1$prob_table_expanded[ds1$row_index$info_hashid, ds1$col_index$xdelta_hashid] <- ds1$prob_table
    
    prob_table_xdeltainfo <- ds0$prob_table_expanded + ds1$prob_table_expanded
    # prob_table_xdeltainfo is the pmf table with size unique(info) * unique(xdelta) = length(info_hashid_combined) * length(xdelta_hashid_combined)
    # and the (i,j)-th entry is the probability of P(info_hashid_combined[i], xdelta_hashid_combined[j])
    
    # collapse prob_table_xdeltainfo so that info_loss is integrated out
    
    info_info0_mapping_table <- unique(dta[, c("info_hashid", "info0_hashid")])
    info_info0_mapping_table <- info_info0_mapping_table[order(info_info0_mapping_table$info_hashid), ]
    df_tmp <- as_tibble(cbind(info_info0_mapping_table, prob_table_xdeltainfo))
    prob_table_xdeltainfo0 <- group_by(df_tmp, info0_hashid) %>% summarize_all(sum)
    prob_table_xdeltainfo0 <- as.matrix(prob_table_xdeltainfo0[, 3:ncol(prob_table_xdeltainfo0)])
    
    if (0) {
        # test the above collapse functionality
        mapping_table <- data.frame(id_fine = 1:6, id_coarse = rep(1:3, each = 2))
        num_matrix <- matrix(1:12, nrow = 6, ncol = 2)
        df_tmp <- as_tibble(cbind(mapping_table, num_matrix))
        num_matrix_collapsed <- group_by(df_tmp, id_coarse) %>% summarize_all(sum)
        num_matrix_collapsed <- as.matrix(num_matrix_collapsed[, 3:ncol(num_matrix_collapsed)])
    }
    
    ### 2) compute Pinfo0 ###
    # a vector of marginal probability P(info0) from the prob_table_xdeltainfo0
    Pinfo0 <- rowSums(prob_table_xdeltainfo0)
    
    
    ### 3) compute survProb.info0 ###
    
    # construct (x_uniquevalues, delta_uniquevalues) to calculate the survival probability
    tmp <- unique(dta[, c("x", "delta", "xdelta_hashid")])
    tmp <- tmp[order(tmp$xdelta_hashid), ]
    x_has_na <- any(is.na(tmp$x))
    x_uniquevalues <- tmp$x[!is.na(tmp$x)]
    delta_uniquevalues <- tmp$delta[!is.na(tmp$delta)]
    if (x_has_na) {
        x_non_na_colindex <- 1:(ncol(prob_table_xdeltainfo0) - 1)
    } else {
        x_non_na_colindex <- 1:ncol(prob_table_xdeltainfo0)
    }
    Pxdelta.info0_noNA <- normalizeByRow(prob_table_xdeltainfo0[, x_non_na_colindex])
    
    # We have two ways to compute survProb.info0:
    # (1) NA: using getSurvProb.info0_vectorized() based on Nelson-Aalen estimator.
    # (2) KM: using survfit() based on Kaplan-Meier estimator
    
    # (1) is faster.
    # (1) and (2) should be equal if there are no numerical approximation error.
    # Hence, when n is large and all values are discrete, (1) and (2) should be equal.
    # In reality with continuous info and continuous failure time, (1) and (2) can be different for certain info and time_question values.
    # Caution: there might be some bug in NA so that currently only KM gives correct answers.
    
    if (survfit_method %in% c("NA", "Nelson-Aalen")) {
        # Method (1):
        survProb.info0 <- getSurvProb.info0_vectorized(Pxdelta.info0_noNA, x_uniquevalues, delta_uniquevalues, time_question)
    } else if (survfit_method %in% c("KM", "Kaplan-Meier")) {
        # Method (2):
        # estimand is P(T > time_question)
        index_in_KMfit_for_time_question <- rev(which(x_uniquevalues < time_question))[1]
        # Below, the KM_fit$surv gives P(T > KM_fit$time), hence this defintion of index_in_KMfit_for_time_question
        survProb.info0 <- rep(NA, length(Pinfo0))
        for (i_info0 in 1:nrow(Pxdelta.info0_noNA)) {
            # Fit Kaplan-Meier with weights being the probabilities
            # From survfit documentation:
            # "The weights must be nonnegative and it is strongly recommended that they be strictly positive,
            # since zero weights are ambiguous, compared to use of the subset argument."
            index_nonzero_weight <- which(Pxdelta.info0_noNA[i_info0, ] > 0)
            KM_fit <- survfit(Surv(x_uniquevalues[index_nonzero_weight], delta_uniquevalues[index_nonzero_weight]) ~ 1,
                              weights = Pxdelta.info0_noNA[i_info0, index_nonzero_weight])
            index_smaller_x <- which(x_uniquevalues[index_nonzero_weight] < time_question)
            if (length(index_smaller_x) >= 1) {
                index_in_KMfit_for_time_question <- index_smaller_x[length(index_smaller_x)]
                survProb.info0[i_info0] <- KM_fit$surv[index_in_KMfit_for_time_question]
            } else {
                survProb.info0[i_info0] <- 1
            }
        }
    } 

    ### 3) compute survProb ###
    
    survProb <- sum(Pinfo0 * survProb.info0)
    
    return(survProb)
}


getSurvProb.info0_vectorized <- function(Pxdelta.info0_noNA,
                                        x_uniquevalues,
                                        delta_uniquevalues,
                                        time_question
) {
    ##########
    # This function is used in getSurvProbFromMasAndProb() when survfit_method = "NA".
    # I think this function needs to be revised. Currently the main function always uses survfit_method = "KM".
    #
    ##########
    
    # estimand is P(T > time_question)
    index_smaller_x <- which(x_uniquevalues < time_question)
    if (length(index_smaller_x) >= 1) {
        # Pxdelta.info0_noNA is a probability matrix for (x,delta) | info0, excluding the (x,delta)=(NA,NA)
        
        ### get Nelson-Aalen estimator for cumulative hazard ###
        
        R <- rowCumsum_from_last(Pxdelta.info0_noNA) # R: size of risk set
        d <- Pxdelta.info0_noNA * rep.row(delta_uniquevalues, nrow(Pxdelta.info0_noNA)) # d: deaths occurring at each time
        cumulative_hazard <- d / R # cumulative hazard
        
        cumulative_hazard[is.nan(cumulative_hazard)] <- 0 # define 0 / 0 = 0
        # cumulative_hazard[R < .Machine$double.eps ] <- 0
        
        # 20160920: debug
        # dLambda[is.na(dLambda)] <- 0 # NA is from 0/0 in dN/Y
        # if (any(is.na(dLambda))) {
        #     browser()
        # }
        
        ### get survival probability from cumulative hazard ###
        
        survProb_mat <- exp(- rowCumsum(cumulative_hazard)) # P(T > time)
        # estimand is P(T > time_question)
        index_for_time_question <- index_smaller_x[length(index_smaller_x)]
        return(survProb_mat[, index_for_time_question])
    } else {
        return(rep(1, nrow(Pxdelta.info0_noNA)))
    }
}


perturbPmas <- function(masAndProb, obs, eps) {
    ##########
    # Return the masAndProb perturbed by obs (with probility mass equals eps).
    #
    # Keyword Arguments:
    # masAndProb            -- a list of two lists: ds0 and ds1; see getExtendedMasAndProbFromData() for details
    # obs                   -- a row in the dta data.frame
    # eps                   -- epsilon used in numeric Gateaux derivative
    #
    ##########
    
    masAndProb$ds0$prob_table <- (1 - eps) * masAndProb$ds0$prob_table
    masAndProb$ds1$prob_table <- (1 - eps) * masAndProb$ds1$prob_table
    
    # add eps to the corresponding cell in the prob_table
    if (obs$robs == 0) {
        irow_prob_table <- which(masAndProb$ds0$row_index$info_hashid == obs$info_hashid)
        icol_prob_table <- which(masAndProb$ds0$col_index$xdelta_hashid == obs$xdelta_hashid)
        masAndProb$ds0$prob_table[irow_prob_table, icol_prob_table] <- masAndProb$ds0$prob_table[irow_prob_table, icol_prob_table] + eps
    } else if (obs$robs == 1) {
        irow_prob_table <- which(masAndProb$ds1$row_index$info_hashid == obs$info_hashid)
        icol_prob_table <- which(masAndProb$ds1$col_index$xdelta_hashid == obs$xdelta_hashid)
        masAndProb$ds1$prob_table[irow_prob_table, icol_prob_table] <- masAndProb$ds1$prob_table[irow_prob_table, icol_prob_table] + eps
    }
    
    return(masAndProb)
}



getExtendedMasAndProbFromData <- function(dta,
                                          info_varname, 
                                          models,
                                          tuning_param_value,
                                          xdelta_working_model) {
    ##########
    # Return a list of two objects:
    #     ds0 is the discretized support \Omega_0 and its probability table that reflects the value of alpha.
    #     ds1 is the discretized support \Omega_1 and its probability table that reflects the value of alpha.
    #
    # Keyword Arguments:
    # dta                   -- entire data set
    # info_varname          -- a vector of all variable names in info_0 and info_loss
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs) for both robs0 and robs1
    # xdelta_working_model  -- the working model used for constructing P(x, delta | s, info, robs), one of c("coxph", "lognormal")
    #
    ##########
    
    dta.robs0 <- subset(dta, robs == 0)
    ds0 <- getExtendedPxdelta_s_info.robs(dta.robs0, info_varname, models, tuning_param_value, xdelta_working_model, robs_type = 0)
    
    dta.robs1 <- subset(dta, robs == 1)
    ds1 <- getExtendedPxdelta_s_info.robs(dta.robs1, info_varname, models, tuning_param_value, xdelta_working_model, robs_type = 1)
    
    Probs0 <- nrow(dta.robs0) / nrow(dta)
    Probs1 <- 1 - Probs0
    ds0$prob_table <- ds0$prob_table * Probs0
    ds1$prob_table <- ds1$prob_table * Probs1
    
    return(list(ds0 = ds0, ds1 = ds1))
}


getExtendedPxdelta_s_info.robs <- function(dta.robs,
                                           info_varname, 
                                           models,
                                           tuning_param_value,
                                           xdelta_working_model,
                                           robs_type
) {
    
    ##########
    # Return the ds (either ds0 or ds1) with probability table computed.
    #
    # Keyword Arguments:
    # dta.robs              -- data set with robs = robs_type
    # info_varname          -- a vector of all variable names in info_0 and info_loss
    # models                -- a list of model_t.robs0, model_c.robs0, model_s.robs0, model_t.robs1, model_c.robs1.
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    # xdelta_working_model  -- the working model used for constructing P(x, delta | s, info, robs), one of c("coxph", "lognormal")
    # robs_type             -- 0 or 1
    #
    # Method:
    #   1) create extended Pxdelta.info_robs: prob(x, delta | info, robs)
    #   2) create Ps.info_robs: prob(s = 1 | info, robs)
    #   3) create Pinfo.robs: prob(info | robs)
    #   4) combine the above three to get prob((x, delta), s, info | robs = 0)
    ##########
    
    ### order dta.robs: first s=1, then s=0. Among s=1, x is increasing. ###
    # This is needed in calculating the survival probability for Cox model (i.e., Pxdelta.info_robs).
    
    stopifnot(!is.unsorted(rev(dta.robs$s))) # require s to be decreasing
    stopifnot(!is.unsorted(dta.robs$x, na.rm = TRUE)) # require x to be increasing
    
    ### create discretized support ds for robs = robs_type ###
    
    row_index <- unique(dta.robs[, c("info_hashid", info_varname)])
    row_index <- row_index[order(row_index$info_hashid), ]
    row_index$prob_table_rowid <- 1:nrow(row_index)
    
    col_index <- unique(dta.robs[, c("xdelta_hashid", "x", "delta")])
    col_index <- col_index[order(col_index$xdelta_hashid), ]
    col_index$prob_table_colid <- 1:nrow(col_index)
    
    prob_table <- matrix(NA, nrow = nrow(row_index), ncol = nrow(col_index))
    
    ds <- list(row_index = row_index, col_index = col_index, prob_table = prob_table)
    
    ### 1) create extended Pxdelta.info_robs ###
    ### prob(x, delta | info, robs) ###
    
    # Pxdelta.info_robs is a probability table with size unique(info) * (# unique(x,delta) without NA).
    # I.e., with size nrow(ds$row_index) * nrow(ds$col_index - 1)    (if there are (x,delta)=(NA,NA) in the original data with this robs).
    
    if (robs_type == 0) {
        covariates_t <- models$model_t.robs0
        covariates_c <- models$model_c.robs0
    } else if (robs_type == 1) {
        covariates_t <- models$model_t.robs1
        covariates_c <- models$model_c.robs1
    }
    
    if (xdelta_working_model == "coxph") {
        Pxdelta.info_robs <- getPxdelta.info_robs_coxph(dta.robs, covariates_t, covariates_c, ds, robs_type = robs_type)
        # Pxdelta.info_robs <- extendPxdelta.info_robs(Pxdelta.info_robs, tuning_param_value)
        Pxdelta.info_robs <- extendPxdelta.info_robs_ver2(Pxdelta.info_robs, ds$col_index$x[!is.na(ds$col_index$x)], max(dta.robs$c), tuning_param_value)
    } else if  (xdelta_working_model == "lognormal") {
        Pxdelta.info_robs <- getPxdelta.info_robs_lognormal_with_extension(dta.robs, covariates_t, covariates_c, ds, robs_type = robs_type, tuning_param_value)
    }
    
    
    
    ### 2) create Ps.info_robs ###
    ### prob(s = 1 | info, robs) ###
    
    # Ps.info_robs is a probability vector with length unique(info) = nrow(ds$row_index).
    # It denotes prob(s = 1 | info, robs = ).
    
    if (robs_type == 0) {
        Ps.info_robs <- getPs.info(dta.robs, models$model_s.robs, ds)
    } else if (robs_type == 1) {
        Ps.info_robs <- rep(1, nrow(ds$row_index))
    }
    
    ### 3) create Pinfo.robs ###
    ### prob(info | robs) ###
    
    # Pinfo.robs is a probability vector with length unique(info) = nrow(ds$row_index)
    
    Pinfo.robs <- getPinfo.robs(dta.robs, ds)
    
    ### 4) combine the above three ###
    ### prob((x, delta) * s, s, info | robs) ###
    
    tmp <- multiply_vector_to_each_col_of_matrix(vec = Ps.info_robs, mat = Pxdelta.info_robs)
    # probability (conditional on info) for the part of ds with robs=0, s=1 or for the part of ds with robs=1
    if (robs_type == 0 & any(is.na(ds$col_index$x))) {
        tmp <- cbind(tmp, 1 - Ps.info_robs)
        # combine with probability (conditional on info) for the part of ds with s =  for robs = 
    }
    
    Pxdelta_s_info.robs <- multiply_vector_to_each_col_of_matrix(vec = Pinfo.robs, mat = tmp)
    
    stopifnot(identical(dim(ds$prob_table), dim(Pxdelta_s_info.robs)))
    ds$prob_table <- Pxdelta_s_info.robs
    
    return(ds)
}


getPinfo.robs <- function(dta.robs, ds) {
    ##########
    # Return Pinfo.robs: prob(info | robs) as a vector.
    #
    # Keyword Arguments:
    # dta.robs  -- data set with a fixed robs value (either 0 or 1)
    # ds        -- the corresponding discretized support, to compute Pinfo.robs for
    #
    ##########
    
    # ds$row_index should have been sorted by prob_table_rowid
    stopifnot(!is.unsorted(ds$row_index$prob_table_rowid))
    
    # Sort info_hashid in both ds$row_index and dta.robs, then count the frequency.
    # (This is for speed; alternatively one could have used a for loop with which() function.)
    
    row_index <- ds$row_index
    row_index <- row_index[order(row_index$info_hashid), ]
    
    tabular_info <- table(dta.robs$info_hashid)
    if (!identical(as.integer(names(tabular_info)), row_index$info_hashid)) {
        warning("In getPinfo.robs: identical(as.integer(names(tabular_info)), row_index$info_hashid) is FALE. Resorting tabular_info.")
        tabular_info <- tabular_info[order(as.integer(names(tabular_info)))]
    }
    stopifnot((identical(as.integer(names(tabular_info)), row_index$info_hashid)))
    row_index$freq <- as.vector(tabular_info)
    row_index$prob <- row_index$freq / sum(row_index$freq)
    
    # sort row_index back accordingly to prob_table_rowid
    row_index <- row_index[order(row_index$prob_table_rowid), ]
    return(row_index$prob)
}


getPxdelta.info_robs_coxph <- function(dta.robs,
                                       covariates_t,
                                       covariates_c,
                                       ds, # a list of row_index, col_index, prob_table
                                       robs_type
) {

    ##########
    # Return a matrix of prob(x, delta | info, robs) for the corresponding robs,
    # with nrow = # unique(info), ncol = # unique(x,delta) (after excluding NA's).
    # Using Cox proportional hazards model.
    #
    # Keyword Arguments:
    # dta.robs      -- data set with robs = robs_type
    # covariates_t  -- covariates to be used in the cox model for T
    # covariates_c  -- covariates to be used in the cox model for C
    # ds            -- the corresponding discretized support, to compute Pxdelta.info_robs for
    # robs_type     -- 0 or 1
    #
    # Method:
    #   See Section 3.3 of the paper.
    #   Note: this function doesn't include tuning_param_value as an argument
    #         (unlike getPxdelta.info_robs_lognormal_with_extension() below),
    #         because for the coxph working model, tuning parameter alpha is added
    #         onto Pxdelta directly.
    ##########
    

    # Check that if there is NA in ds$col_index$x and ds$col_index$delta, then NA should be the last element.
    # (Otherwise the following definition of delta_forcalc is incorrect.)
    stopifnot(!any(is.na(ds$col_index$x)) | is.na(ds$col_index$x[length(ds$col_index$x)]))
    stopifnot(!any(is.na(ds$col_index$delta)) | is.na(ds$col_index$delta[length(ds$col_index$delta)]))
    
    #### Note:
    #### In the rest of this function,
    #### all comments that refer to unique(x,delta) or unique(x)
    #### mean the unique values after excluding NA's.
    
    n_unique_info <- nrow(ds$row_index)
    n_unique_xdelta <- sum(!is.na(ds$col_index$x)) # number of unique (x,delta) pairs
    
    if (robs_type == 0) {
        dta.robs_s1 <- subset(dta.robs, s == 1)
    } else {
        dta.robs_s1 <- dta.robs
    }
    
    dta.robs_s1_deltaflipped <- dta.robs_s1
    dta.robs_s1_deltaflipped$delta <- 1 - dta.robs_s1_deltaflipped$delta
    
    # Handling observations with same x but opposite delta.
    # This needs special handling because the following survival_t and survival_c from Coxph
    # has ncol = # unique(x), but Pxdelta has ncol = # unique(x,delta).
    x_has_duplicates <- anyDuplicated(ds$col_index$x)
    if (x_has_duplicates) {
        index_in_col_index_with_equal_x_but_opposite_delta <- sort(unique(c(which(duplicated(ds$col_index$x)), which(duplicated(ds$col_index$x, fromLast = TRUE)))))
        index_col_to_copy_sequentially_in_pmf_and_S <- index_in_col_index_with_equal_x_but_opposite_delta[seq(from = 1, to = length(index_in_col_index_with_equal_x_but_opposite_delta), by = 2)]
        # If x is a duplicate in unique(x,delta), then the same x will appear exactly twice.
        # index_col_to_copy_sequentially_in_pmf_and_S are the first index of each duplicated x in ds$col_index$x.
    }
        
    # Below I added (x_max+1, 0) and (x_max+2, 1) to the unique values of (x,delta) to handle the issue that
    # the estimated survival probability doesn't equal to 0 at the last observed survival time (if the last time is delta = 1)
    # or the estimated survival probability doesn't equal to 0 at the last observed censoring time (if the last time is delta = 0).
    # This issue is due to the fact that the estimated cumulative hazard cannot equal to infinity.
    
    if (all(dta.robs_s1$delta == 0)) {
        # all delta = 0: no observed surrival time, all observations are censored
        pmf_t <- matrix(0, nrow = nrow(ds$prob_table), ncol = n_unique_xdelta)
        pmf_t <- cbind(pmf_t, 0, 1) # only 1 probability at x_max + 2
        St <- 1 - rowCumsum(pmf_t)
        
    }  else if (sum(ds$col_index$delta == 1, na.rm = TRUE) <= 2) {
        # use empirical frequency for pmf_t for the few x values with delta = 1
        pmf_t <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
        id_in_col_index <- which(ds$col_index$delta == 1)
        prob_t_nonzero <- 1 / length(id_in_col_index)
        
        pmf_t[, ds$col_index$prob_table_colid[id_in_col_index]] <- prob_t_nonzero
        pmf_t <- cbind(pmf_t, 0, 0) # 0 probability at x_max+1 and x_max+2
        St <- 1 - rowCumsum(pmf_t)
        
    } else {
        
        ### Cox model for T ###
        
        # remove variables with any NA values from covariates_t
        covariates_t_iter <- covariates_t
        for (covar in covariates_t_iter) {
            if (any(is.na(dta.robs_s1[, covar]))) {
                covariates_t <- setdiff(covariates_t, covar)
            }
        }
        
        formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates_t, collapse = "+")))
        model_t.info_robs <- coxph(formula, data = dta.robs_s1)
        
        baseline_cumulative_hazard <- basehaz(model_t.info_robs, centered = FALSE)
        # above line: Breslow's estimator for baseline hazard
        # see p.320 eq. (12.3.2) of Lee - statistical methods for survival data analysis.
        log_survival <- outer(as.vector(exp(as.matrix(ds$row_index[, covariates_t]) %*% coef(model_t.info_robs))),
                              - baseline_cumulative_hazard$hazard)
        # above line is based on log S(t) = exp(Z * beta) * {- H_0(t)}
        survival_t <- exp(log_survival) 
        # survival_t: a table of survival (cumulative) probabilities for T,
        # size is nrow(ds$prob_table) * # unique(x)
        # Note that each row of survival_t does not start at 1 and does not end at 0;
        # this means that there are >0 probability that T > x_max.
        # To fix this, we introduce x_max + 1, x_max + 2 in the following x_forcalc.
        
        pmf_t <- cbind(1, survival_t) - cbind(survival_t, 0)
        pmf_t <- cbind(pmf_t[, 1:(ncol(pmf_t)-1)], 0, pmf_t[, ncol(pmf_t)])
        # probability mass function for P(T = x | info, robs)
        # size is nrow(ds$prob_table) * (# unique(x) + 2)
        # This +2 is due to the introduction of x_max + 1, x_max + 2 in the following x_forcalc.
        # The second to the last column is all 0 for pmf_t, because that is for an artificial point of (x_max+1, 0).
        
        if (0) {
            # check equality with st, when all (x, delta) and all info are distinct
            # survival[, dta.robs_s1$delta==1] is a differently-ordered version of st
            # This is very strange! I don't understand why. Maybe there is something wrong with surv.from.cox().
            st <- surv.from.cox(dta.robs_s1, model_t.info_robs, dta.robs, covariates_t)[, dta.robs_s1$delta==1]
            tmp1 <- survival_t[, dta.robs_s1$delta==1]
            sort(st)==sort(tmp1)
            summary(sort(st)-sort(tmp1))
        }
        
        # St: prob(T > x | info, robs), of size nrow(ds$prob_table) * (# unique(x) + 2)
        St <- 1 - rowCumsum(pmf_t)
        
        if (x_has_duplicates) {
            # handling duplicate x with opposite delta in unique(x,delta)
            for (col_to_copy in index_col_to_copy_sequentially_in_pmf_and_S) {
                pmf_t <- copy_icol(pmf_t, col_to_copy)
                St <- copy_icol(St, col_to_copy)
            }
        }
        # pmf_t and St are now of size nrow(ds$prob_table) * (# unique(x,delta) + 2)
    }
    
    if (all(dta.robs_s1$delta == 1)) {
        # all delta = 1: no censoring, all survival times are observed
        pmf_c <- matrix(0, nrow = nrow(ds$prob_table), ncol = n_unique_xdelta)
        pmf_c <- cbind(pmf_c, 1, 0) # only 1 probability at x_max + 1
        Sc <- 1 - rowCumsum(pmf_c)
        
    } else if (sum(ds$col_index$delta == 0, na.rm = TRUE) <= 2) {
        # use empirical frequency for pmf_c for the few x values with delta = 0
        pmf_c <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
        id_in_col_index <- which(ds$col_index$delta == 0)
        prob_c_nonzero <- 1 / length(id_in_col_index)
        
        pmf_c[, ds$col_index$prob_table_colid[id_in_col_index]] <- prob_c_nonzero
        
        pmf_c <- cbind(pmf_c, 0, 0)# 0 probability at x_max+1 and x_max+2
        Sc <- 1 - rowCumsum(pmf_c)
        
    } else {
        
        ### Cox model for C ###
        
        # remove variables with any NA values from covariates_c
        covariates_c_iter <- covariates_c
        for (covar in covariates_c_iter) {
            if (any(is.na(dta.robs_s1_deltaflipped[, covar]))) {
                covariates_c <- setdiff(covariates_c, covar)
            }
        }
        
        formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates_c, collapse = "+")))
        model_c.info_robs <- coxph(formula, data = dta.robs_s1_deltaflipped)
        
        baseline_cumulative_hazard <- basehaz(model_c.info_robs, centered = FALSE)
        # above line: Breslow's estimator for baseline hazard
        # see p.320 eq. (12.3.2) of Lee - statistical methods for survival data analysis.
        log_survival <- outer(as.vector(exp(as.matrix(ds$row_index[, covariates_c]) %*% coef(model_c.info_robs))),
                              - baseline_cumulative_hazard$hazard)
        # above line is based on log S(t) = exp(Z * beta) * {- H_0(t)}
        survival_c <- exp(log_survival) 
        # survival_c: a table of survival (cumulative) probabilities for C,
        # size is nrow(ds$prob_table) * # unique(x)
        # Note that each row of survival_c does not start at 1 and does not end at 0;
        # this means that there are >0 probability that C > x_max.
        # To fix this, we introduce x_max + 1, x_max + 2 in the following x_forcalc.
        
        pmf_c <- cbind(1, survival_c) - cbind(survival_c, 0)
        pmf_c <- cbind(pmf_c, 0)
        # probability mass function for P(C = x | info, robs)
        # size is nrow(ds$prob_table) * (# unique(x) + 2).
        # This +2 is due to the introduction of x_max + 1, x_max + 2 in the following x_forcalc.
        # The last column is all 0 for pmf_c, because that is for an artificial point of (x_max+2, 1).
        
        # Sc: prob(C > x | info, robs), of size nrow(ds$prob_table) * (# unique(x) + 2)
        Sc <- 1 - rowCumsum(pmf_c)
        
        if (x_has_duplicates) {
            # handling duplicate x with opposite delta in unique(x,delta)
            for (col_to_copy in index_col_to_copy_sequentially_in_pmf_and_S) {
                pmf_c <- copy_icol(pmf_c, col_to_copy)
                Sc <- copy_icol(Sc, col_to_copy)
            }
        }
        # pmf_c and Sc are now of size nrow(ds$prob_table) * (# unique(x,delta) + 2)
    }
    
    
    # calculation of prob(x, delta | info, robs)
    delta_forcalc <- c(ds$col_index$delta[1:n_unique_xdelta], 0, 1)
    Pxdelta <- pmf_t * Sc * rep.row(delta_forcalc, n_unique_info) + pmf_c * St * (1 - rep.row(delta_forcalc, n_unique_info))
    

    
    # make it so that Pxdelta too close to 0 be some small number (added 2017.04.08)
    # I don't remember why I added this...
    # Pxdelta[Pxdelta <= 1e-8] <- 1e-8
    
    # removing prob on the added (xmax+1, 0) and (xmax+2, 1), and renormalize
    Pxdelta <- normalizeByRow(Pxdelta[, 1:n_unique_xdelta])
    
    return(Pxdelta)
    
    ### Calibrate prob(x, delta | info, robs) ###
    # see 2015.12.20_calibration of discrete probability on observed data.pdf
    # not needed for this estimand, since we are already working on CDF not PDF
    
    # source.cov <- c(rep(1, n_s1), rep(0, n_robs - n_s1))
    # denominator <- as.vector(t(source.cov) %*% Pxdelta)
    # qhat_ij <- Pxdelta / rep.row(denominator, nrow(Pxdelta))
    # 
    # return(normalizeByRow(qhat_ij))
}


getPxdelta.info_robs_lognormal_with_extension <- function(dta.robs,
                                                          covariates_t,
                                                          covariates_c,
                                                          ds, # a list of row_index, col_index, prob_table
                                                          robs_type,
                                                          tuning_parameter_value
) {
    ##########
    # Return a matrix of prob(x, delta | info, robs) for the corresponding robs,
    # with nrow = # unique(info), ncol = # unique(x,delta) (after excluding NA's).
    # Using Cox proportional hazards model.
    #
    # Keyword Arguments:
    # dta.robs              -- data set with robs = robs_type
    # covariates_t          -- covariates to be used in the cox model for T
    # covariates_c          -- covariates to be used in the cox model for C
    # ds                    -- the corresponding discretized support, to compute Pxdelta.info_robs for
    # robs_type             -- 0 or 1
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    #
    # Method:
    #   See Section 3.3 of the paper.
    #   Note: this function includes tuning_param_value as an argument
    #         (unlike getPxdelta.info_robs_coxph() above),
    #         because for the lognormal working model, tuning parameter alpha is added
    #         as a regression coefficient in the lognormal regression fit.
    ##########
    
    stopifnot(!any(is.na(ds$col_index$x)) | is.na(ds$col_index$x[length(ds$col_index$x)]))
    stopifnot(!any(is.na(ds$col_index$delta)) | is.na(ds$col_index$delta[length(ds$col_index$delta)]))
    
    n_unique_info <- nrow(ds$row_index)
    n_unique_xdelta <- sum(!is.na(ds$col_index$x)) # number of unique (x,delta) pairs (excluding NA)
    
    if (robs_type == 0) {
        dta.robs_s1 <- subset(dta.robs, s == 1)
    } else {
        dta.robs_s1 <- dta.robs
    }
    
    dta.robs_s1_deltaflipped <- dta.robs_s1
    dta.robs_s1_deltaflipped$delta <- 1 - dta.robs_s1_deltaflipped$delta
    
    x_unique_values <- ds$col_index$x[!is.na(ds$col_index$x)]
    stopifnot(!is.unsorted(x_unique_values))
    midpoint_times <- c(0,
                        (x_unique_values[1:(length(x_unique_values)-1)] + x_unique_values[2:length(x_unique_values)]) / 2,
                        Inf)
    
    if (all(dta.robs_s1$delta == 0)) {
        # all delta = 0: no observed surrival time, all observations are censored
        pmf_t <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
    } else if (sum(ds$col_index$delta == 1, na.rm = TRUE) <= 2) {
        # use empirical frequency for pmf_t for the few x values with delta = 1
        pmf_t <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
        id_in_col_index <- which(ds$col_index$delta == 1)
        prob_t_nonzero <- 1 / length(id_in_col_index)
        
        pmf_t[, ds$col_index$prob_table_colid[id_in_col_index]] <- prob_t_nonzero
        
    } else {
        ### lognormal regression for T ###
        
        # remove variables with any NA values from covariates_t
        covariates_t_iter <- covariates_t
        for (covar in covariates_t_iter) {
            if (any(is.na(dta.robs_s1[, covar]))) {
                covariates_t <- setdiff(covariates_t, covar)
            }
        }
        
        formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates_t, collapse = "+")))
        model_t.info_robs <- survreg(formula, data = dta.robs_s1, dist = "lognormal")
        
        # tuning parameter is added to the intercept of the fitted lognormal regression for T
        coefs <- coef(model_t.info_robs)
        coefs[1] <- coefs[1] + tuning_parameter_value
        
        mean_lognormal <- as.matrix(cbind(1, ds$row_index[, covariates_t])) %*% coefs
        scale_lognormal <- model_t.info_robs$scale
        
        # discretize the fitted lognormal to observed x values
        pmf_t <- matrix(NA, nrow = n_unique_info, ncol = n_unique_xdelta)
        for (i_unique_x in 1:length(x_unique_values)) {
            left_midpoint <- midpoint_times[i_unique_x]
            right_midpoint <- midpoint_times[i_unique_x + 1]
            pmf_t[, i_unique_x] <- psurvreg(right_midpoint, mean = mean_lognormal, scale = scale_lognormal, dist = "lognormal") -
                psurvreg(left_midpoint, mean = mean_lognormal, scale = scale_lognormal, dist = "lognormal")
        }
    }
    
    if (all(dta.robs_s1$delta == 1)) {
        # all delta = 1: no censoring, all survival times are observed
        pmf_c <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
    } else if (sum(ds$col_index$delta == 0, na.rm = TRUE) <= 2) {
        # use empirical frequency for pmf_c for the few x values with delta = 0
        pmf_c <- matrix(0, nrow = n_unique_info, ncol = n_unique_xdelta)
        
        id_in_col_index <- which(ds$col_index$delta == 0)
        prob_c_nonzero <- 1 / length(id_in_col_index)
        
        pmf_c[, ds$col_index$prob_table_colid[id_in_col_index]] <- prob_c_nonzero

    } else {
        
        ### lognormal regression for T ###
        
        # remove variables with any NA values from covariates_c
        covariates_c_iter <- covariates_c
        for (covar in covariates_c_iter) {
            if (any(is.na(dta.robs_s1_deltaflipped[, covar]))) {
                covariates_c <- setdiff(covariates_c, covar)
            }
        }
        
        formula <- as.formula(paste0("Surv(x, delta) ~ ", paste(covariates_c, collapse = "+")))
        model_c.info_robs <- survreg(formula, data = dta.robs_s1_deltaflipped, dist = "lognormal")
        
        mean_lognormal <- as.matrix(cbind(1, ds$row_index[, covariates_c])) %*% coef(model_c.info_robs)
        scale_lognormal <- model_c.info_robs$scale
        
        # discretize the fitted lognormal to observed x values
        pmf_c <- matrix(NA, nrow = n_unique_info, ncol = n_unique_xdelta)
        for (i_unique_x in 1:length(x_unique_values)) {
            left_midpoint <- midpoint_times[i_unique_x]
            right_midpoint <- midpoint_times[i_unique_x + 1]
            pmf_c[, i_unique_x] <- psurvreg(right_midpoint, mean = mean_lognormal, scale = scale_lognormal, dist = "lognormal") -
                psurvreg(left_midpoint, mean = mean_lognormal, scale = scale_lognormal, dist = "lognormal")
        }
    }
    
    # St: prob(T > x | info, robs), of size nrow(ds$prob_table) * (ncol(ds$prob_table) + 1)
    # Sc: prob(C > x | info, robs), of size nrow(ds$prob_table) * (ncol(ds$prob_table) + 1)
    St <- 1 - rowCumsum(pmf_t)
    Sc <- 1 - rowCumsum(pmf_c)
    
    # calculation of prob(x, delta | info, robs)
    delta_forcalc <- ds$col_index$delta[1:n_unique_xdelta]
    Pxdelta <- pmf_t * Sc * rep.row(delta_forcalc, n_unique_info) + pmf_c * St * (1 - rep.row(delta_forcalc, n_unique_info))
    
    Pxdelta <- normalizeByRow(Pxdelta)
    
    return(Pxdelta)
}


extendPxdelta.info_robs <- function(Pxdelta.info_robs, tuning_param_value) {
    ##########
    # Return Pxdelta.info_robs that incorporates the tuning parameter.
    #
    # Keyword Arguments:
    # Pxdelta.info_robs     -- Probablity table (a matrix) for Pxdelta.info_robs
    # tuning_param_value    -- value of tuning parameter, used in extending prob(x, delta | info, robs)
    #
    # Note: this function is obselete, because this does not allow sum(Eif) to explore values around 0.
    #       We are currently using extendPxdelta.info_robs_ver2() instead.
    #
    ##########
    
    factor <- seq(from = 1, to = tuning_param_value, length.out = ncol(Pxdelta.info_robs))
    Pxdelta.info_robs <- Pxdelta.info_robs * rep.row(factor, nrow(Pxdelta.info_robs))
    return(normalizeByRow(Pxdelta.info_robs))
}

extendPxdelta.info_robs_ver2 <- function(Pxdelta.info_robs, x_unique_values_withoutNA, cmax, alpha) {
    ##########
    # Return Pxdelta.info_robs that incorporates the tuning parameter.
    #
    # Keyword Arguments:
    # Pxdelta.info_robs            -- Probablity table (a matrix) for Pxdelta.info_robs
    # x_unique_values_withoutNA    -- The vector of unique values of x, sorted increasingly
    # cmax                         -- max(C) for the corresponding robs
    # alpha                        -- tuning parameter value
    #
    # Method:
    #    See Section 3.3 of the paper.
    #
    ##########
    
    factor <- 1 + alpha * x_unique_values_withoutNA / cmax
    # I don't remember why I wrote the following line.
    # if (length(factor) < ncol(Pxdelta.info_robs)) {
    #     factor <- c(factor, 1)
    # }
    stopifnot(length(factor) == ncol(Pxdelta.info_robs))
    Pxdelta.info_robs <- Pxdelta.info_robs * rep.row(factor, nrow(Pxdelta.info_robs))
    Pxdelta.info_robs[Pxdelta.info_robs < 0] <- 0
    Pxdelta.info_robs <- normalizeByRow(Pxdelta.info_robs)
    return(Pxdelta.info_robs)
}

getPs.info <- function(dta.robs0, covariates, ds0) {
    ##########
    # Return a vector of estimated propensity score: prob(s=1 | covariates, robs=0)
    #
    # Keyword Arguments:
    # dta.robs0   -- data set with robs = 0
    # covariates  -- covariates to be used in propensity score model
    # ds0         -- the discretized support for robs = 0
    ##########
    
    # If all S=1, then return a vector of 1
    if (all(dta.robs0$s == 1)) {
        return(rep(1, nrow(ds0$row_index)))
    } else {
        if (is.null(covariates)) {
            covariates_term <- 1
        } else {
            covariates_term <- paste(covariates, collapse = "+")
        }
        formula <- as.formula(paste0("s ~ ", covariates_term))
        model_s <- glm(formula, data = dta.robs0, family = binomial(link = logit))
        prob_s_equal_1 <- predict(model_s, newdata = ds0$row_index, type = "response")
        return(prob_s_equal_1)
    }

}



# Utility functions -------------------------------------------------------

cumsum_from_last <- function(x) {
    # calculate cumulative sum of a vector from the last element instead of the first
    # a reversed version of cumsum() in basic R
    
    n <- length(x)
    return(cumsum(x[n:1])[n:1])
}

rowCumsum_from_last <- function(mat) {
    # calculate cumulative sum of each row for a mattrix from the last element instead of the first
    # a reversed version of rowCumsum() defined below
    
    n <- ncol(mat)
    return(rowCumsum(mat[, n:1])[, n:1])
    
    # Example:
    # > x
    #      [,1] [,2] [,3] [,4]
    # [1,]    1    4    7   10
    # [2,]    2    5    8   11
    # [3,]    3    6    9   12
    # 
    # > rowCumsum_from_last(x)
    #      [,1] [,2] [,3] [,4]
    # [1,]   22   21   17   10
    # [2,]   26   24   19   11
    # [3,]   30   27   21   12
}

rowCumsum <- function(mat) {
    # cumsum for each row in a matrix
    
    return(t(apply(mat, 1, cumsum)))
    
    # Example:
    # > x
    #      [,1] [,2] [,3] [,4]
    # [1,]    1    4    7   10
    # [2,]    2    5    8   11
    # [3,]    3    6    9   12
    # 
    # > t(apply(x, 1, cumsum))
    #      [,1] [,2] [,3] [,4]
    # [1,]    1    5   12   22
    # [2,]    2    7   15   26
    # [3,]    3    9   18   30
}

normalizeByRow <- function(mat, add_equal_mass_for_zero_rows = TRUE) {
    # normalize a nonnegative matrix by row so that each row adds up to 1
    # If add_equal_mass_for_zero_rows is TRUE, an equal weight is added to each cell of the rows that has sum 0.
    
    if (add_equal_mass_for_zero_rows) {
        which_row_to_add_mass <- which(rowSums(mat) < .Machine$double.eps^0.5)
        mat[which_row_to_add_mass, ] <- 1
    }
    return(mat * (1 / rowSums(mat)))
    # > zz
    #      [,1] [,2] [,3] [,4]
    # [1,]    1    4    7   10
    # [2,]    2    5    8   11
    # [3,]    3    6    9   12
    # > normalizeByRow(zz)
    #            [,1]      [,2]      [,3]      [,4]
    # [1,] 0.04545455 0.1818182 0.3181818 0.4545455
    # [2,] 0.07692308 0.1923077 0.3076923 0.4230769
    # [3,] 0.10000000 0.2000000 0.3000000 0.4000000
    
    # > zz
    #      [,1] [,2] [,3] [,4]
    # [1,]    0    0    0    0
    # [2,]    1    2    3    4
    # > normalizeByRow(zz)
    #      [,1] [,2] [,3] [,4]
    # [1,] 0.25 0.25 0.25 0.25
    # [2,] 0.10 0.20 0.30 0.40
}

rep.row <- function(x, n) {
    # repeat a row vector x by n times
    
    return(matrix(rep(x, each = n), nrow = n))
    # > rep.row(1:3,5)
    #      [,1] [,2] [,3]
    # [1,]    1    2    3
    # [2,]    1    2    3
    # [3,]    1    2    3
    # [4,]    1    2    3
    # [5,]    1    2    3
}

rep.col <- function(x, n) {
    # repeat a column vector x by n times
    
    return(matrix(rep(x, each = n), ncol = n, byrow = TRUE))
    # > rep.col(1:3,5)
    #      [,1] [,2] [,3] [,4] [,5]
    # [1,]    1    1    2    2    3
    # [2,]    1    1    2    3    3
    # [3,]    1    2    2    3    3
}

multiply_vector_to_each_col_of_matrix <- function(vec, mat) {
    # multiply vec to each column of mat
    
    stopifnot(length(vec) == nrow(mat))
    return(rep.col(vec, ncol(mat)) * mat)
    # multiply_vector_to_each_col_of_matrix(c(1,10), matrix(1:6, nrow = 2))
    #      [,1] [,2] [,3]
    # [1,]    1    3    5
    # [2,]   20   40   60
}

copy_icol <- function(mat, icol) {
    # make a copy the i-th column of mat, and insert it between the i-th and (i+1)-th column in mat 
    if (icol != ncol(mat)) {
        return(cbind(mat[, 1:icol], mat[, icol], mat[, (icol+1):ncol(mat)]))
    } else {
        return(cbind(mat, mat[, ncol(mat)]))
    }
    # copy_icol(matrix(1:6, nrow = 2), icol = 2)
    #      [,1] [,2] [,3] [,4]
    # [1,]    1    3    3    5
    # [2,]    2    4    4    6
}


# The following is obselete and is not used in the current code.
# --- Tianchen Qian, 2019.06.08

#--- begin: needs simplification ---
# Get survival probability from Cox model fits
# This is using Breslow's estimator for baseline hazard
# see p.320 eq. (12.3.2) of Lee - statistical methods for survival data analysis
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
    if (!all(surv <= 1)) {
        stop("all(surv <= 1) is not TRUE.
    If also have Warning message:
        In fitter(X, Y, strats, offset, init, control, weights = weights,  :
        Ran out of iterations and did not converge
    Then it is likely that the coxph routine in parent function did not converge.")
    }
    return(surv)
}
#--- end: needs simplification ---