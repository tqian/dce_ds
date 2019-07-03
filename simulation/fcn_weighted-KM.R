###################################
# Weighted Kaplan-Meier estimator of survival probability in double sampling design
# Source code of functions
#      Tianchen Qian
#        2016.11.19
###################################

# In variable names, x.y means x | y ("." means conditional)

# Updates:
# 2019.06.08: corrected the time index of interest in KM fit
# 2019.06.09: added standard error to the return value

library(survey)

KmWeightedByPropScore <- function(dta,
                                  info_0 = c("age10", "prxcd4v"),
                                  info_loss = c("l"),
                                  time_question = 365) {
    ##########
    # Return a list of:
    # SurvProb  -- estimated survival probablity
    # ci        -- 95% confidence interval (output from suvkm())
    # varlog    -- variance on the log-scale (output from suvkm())
    #
    #
    # Keyword Arguments:
    # dta                   -- data set
    # info_0                -- covariate names of info_0 (baseline variables)
    # info_loss             -- covariate names of info_loss (longitudinal measurements up to loss-to-followup)
    # time_question         -- t in the estimand P(T > t)
    #
    # Method:
    #     A weighted Kaplan-Meier estimator, where
    #     the double-sampled participants are weighted by inverse of the estimated propensity score of being selected.
    ##########
    
    dta.robs0 <- subset(dta, robs == 0)
    dta.robs1 <- subset(dta, robs == 1)
    
    covariates <- c(info_0, info_loss)
    formula <- as.formula(paste0("s ~ ", paste(covariates, collapse = "+")))
    model_s <- glm(formula, data = dta.robs0, family = binomial(link = logit))
    
    dta.robs0$propen_score <- predict(model_s, newdata = dta.robs0, type = "response")
    dta.robs0$weight <- 1 / dta.robs0$propen_score
    dta.robs1$propen_score <- 1
    dta.robs1$weight <- 1
    dt <- rbind(dta.robs0, dta.robs1)
    dt <- subset(dt, !is.na(x))
    
    dt_svy <- svydesign(id = ~1, prob = dt$weight, data = dt)
    
    surv_fit <- svykm(Surv(x, delta) ~ 1, design = dt_svy, se = TRUE)
    
    index_smaller_time <- which(surv_fit$time < time_question)
    if (length(index_smaller_time) >= 1) {
        time_index <- index_smaller_time[length(index_smaller_time)]
        survProb <- surv_fit$surv[time_index]
        ci <- confint(surv_fit,surv_fit$time[time_index],level=0.95)
        varlog <- surv_fit$varlog[time_index]
        return(list(survProb = survProb, ci = ci, varlog = varlog))
    } else {
        return(list(survProb = 1, ci = c(1,1), varlog = NA))
    }
}


KmWeightedByRobs <- function(dta, # unique x
                             time_question = 365) {
    ##########
    # Return a list of:
    # SurvProb  -- estimated survival probablity
    # ci        -- 95% confidence interval (output from suvkm())
    # varlog    -- variance on the log-scale (output from suvkm())
    #
    #
    # Keyword Arguments:
    # dta                   -- data set
    # time_question         -- t in the estimand P(T > t)
    #
    # Method:
    #     A weighted Kaplan-Meier estimator, where
    #     the double-sampled participants are weighted by sum(dta$robs == 0) / sum(dta_noNA$robs == 0)
    ##########
    
    
    dta_noNA <- dta[!is.na(dta$x), ]
    dta_noNA$weight <- 1
    dta_noNA$weight[dta_noNA$robs == 0] <- sum(dta$robs == 0) / sum(dta_noNA$robs == 0)
    
    # # Use survfit with weights argument: the estimator seems very off.
    # survfit_stratified <- survfit(Surv(x, delta) ~ robs, data = dta_noNA, weights = dta_noNA$weight)
    # index_smaller_time <- which(survfit_stratified$time < time_question)
    # 
    # if (length(index_smaller_time) >= 1) {
    #     time_index <- index_smaller_time[length(index_smaller_time)]
    #     survProb <- survfit_stratified$surv[time_index]
    #     se <- survfit_stratified$std.err[time_index]
    #     ci <- c(survProb - 1.96 * se, survProb + 1.96 * se)
    #     return(list(survProb = survProb, se = se, ci = ci))
    # } else {
    #     return(list(survProb = 1, se = NA, ci = c(1,1)))
    # }
    
    
    # Use svykm
    dt_svy <- svydesign(id = ~1, prob = dta_noNA$weight, data = dta_noNA)
    
    surv_fit <- svykm(Surv(x, delta) ~ 1, design = dt_svy, se = TRUE)
    
    index_smaller_time <- which(surv_fit$time < time_question)
    if (length(index_smaller_time) >= 1) {
        time_index <- index_smaller_time[length(index_smaller_time)]
        survProb <- surv_fit$surv[time_index]
        ci <- confint(surv_fit,surv_fit$time[time_index],level=0.95)
        varlog <- surv_fit$varlog[time_index]
        return(list(survProb = survProb, ci = ci, varlog = varlog))
    } else {
        return(list(survProb = 1, ci = c(1,1), varlog = NA))
    }
    
    
    # # Directly compute stratified KM, but don't have a way to calculate the standard errors yet
    # dta.robs1 <- subset(dta, robs == 1)
    # km_robs1 <- survfit(Surv(dta.robs1$x, dta.robs1$delta) ~ 1)
    # 
    # dta.robs0s1 <- subset(dta, (robs == 0) & (s == 1))
    # km_robs0s1 <- survfit(Surv(dta.robs0s1$x, dta.robs0s1$delta) ~ 1)
    # 
    # pr_robs0 <- mean(dta$robs == 0)
    # pr_robs1 <- mean(dta$robs == 1)
    # 
    # index_smaller_time_1 <- which(km_robs1$time < time_question)
    # if (length(index_smaller_time_1) >= 1) {
    #     time_index <- index_smaller_time_1[length(index_smaller_time_1)]
    #     prob_1 <- km_robs1$surv[time_index]
    # } else {
    #     prob_1 <- 1
    # }
    # 
    # index_smaller_time_0 <- which(km_robs0s1$time < time_question)
    # if (length(index_smaller_time_0) >= 1) {
    #     time_index <- index_smaller_time_0[length(index_smaller_time_0)]
    #     prob_0 <- km_robs0s1$surv[time_index]
    # } else {
    #     prob_0 <- 1
    # }
    # 
    # return(pr_robs1 * prob_1 + pr_robs0 * prob_0)
}


KmCompleteCase <- function(dta, # unique x
                           time_question = 365) {
    ##########
    # Return a list of:
    # SurvProb  -- estimated survival probablity
    # se        -- standard error (output from survfit())
    # ci        -- 95% confidence interval (output from survfit())
    #
    #
    # Keyword Arguments:
    # dta                   -- data set
    # time_question         -- t in the estimand P(T > t)
    #
    # Method:
    #     A Kaplan-Meier estimator that only uses those with Robs = 1 or S = 1.
    ##########
    
    dta <- subset(dta, !is.na(x))
    
    surv_fit <- survfit(Surv(dta$x, dta$delta) ~ 1)
    
    index_smaller_time <- which(surv_fit$time < time_question)
    if (length(index_smaller_time) >= 1) {
        time_index <- index_smaller_time[length(index_smaller_time)]
        survProb <- surv_fit$surv[time_index]
        se <- surv_fit$std.err[time_index]
        ci <- c(survProb - 1.96 * se, survProb + 1.96 * se)
        return(list(survProb = survProb, se = se, ci = ci))
    } else {
        return(list(survProb = 1, se = NA, ci = c(1,1)))
    }
}



KmStratifiedByRobs <- function(dta, # unique x
                               time_question = 365,
                               n_boot = 1000) {
    ##########
    # Return a list of:
    # SurvProb  -- estimated survival probablity
    # se        -- standard error (computed by bootstrap)
    # ci        -- 95% confidence interval (computed by bootstrap)
    #
    #
    # Keyword Arguments:
    # dta                   -- data set
    # time_question         -- t in the estimand P(T > t)
    # n_boot                -- the number of bootstrap replicates for calculating se and ci
    #
    # Method:
    #     A Kaplan-Meier estimator using robs as stratification variable.
    #     First it computes the K-M survival probability separately for robs = 0 and robs = 1,
    #     then it takes weighted average of the two with weights being P(robs = 0) and P(robs = 1).
    #     Standard error and confidence interval is calculated by bootstrapping the original data set.
    ##########
    
    survProb <- KmStratifiedByRobs_only_estimate(dta, time_question)
    # KmStratifiedByRobs_only_estimate_for_boot <- function(dta){
        # return(KmStratifiedByRobs_only_estimate(dta, time_question))
    # }
    boot_results <- boot(data = dta, statistic = KmStratifiedByRobs_only_estimate_for_boot, R = n_boot, time_question = time_question)
    se <- sd(boot_results$t)
    ci <- boot.ci(boot_results, type="bca")
    ci_left <- ci$bca[length(ci$bca) - 1]
    ci_right <- ci$bca[length(ci$bca)]
    return(list(survProb = survProb, se = se, ci = c(ci_left, ci_right)))
}


KmStratifiedByRobs_only_estimate_for_boot <- function(dta,
                                                      time_question,
                                                      indices) {
    # Used in KmStratifiedByRobs()
    return(KmStratifiedByRobs_only_estimate(dta[indices, ], time_question))
}

KmStratifiedByRobs_only_estimate <- function(dta, # unique x
                               time_question = 365
                               ) {
    # Used in KmStratifiedByRobs()
    
    if (any(dta$robs == 1)) {
        dta.robs1 <- subset(dta, robs == 1)
        km_robs1 <- survfit(Surv(dta.robs1$x, dta.robs1$delta) ~ 1)
        pr_robs1 <- mean(dta$robs == 1)
        index_smaller_time_1 <- which(km_robs1$time < time_question)
        if (length(index_smaller_time_1) >= 1) {
            time_index <- index_smaller_time_1[length(index_smaller_time_1)]
            prob_1 <- km_robs1$surv[time_index]
        } else {
            prob_1 <- 1
        }
    } else {
        prob_1 <- 0
        pr_robs1 <- 0
    }
    
    if (any(dta$robs == 0)) {
        if (any(dta$s[dta$robs == 0] == 1)) {
            dta.robs0s1 <- subset(dta, (robs == 0) & (s == 1))
            km_robs0s1 <- survfit(Surv(dta.robs0s1$x, dta.robs0s1$delta) ~ 1)
            pr_robs0 <- mean(dta$robs == 0)
            
            index_smaller_time_0 <- which(km_robs0s1$time < time_question)
            if (length(index_smaller_time_0) >= 1) {
                time_index <- index_smaller_time_0[length(index_smaller_time_0)]
                prob_0 <- km_robs0s1$surv[time_index]
            } else {
                prob_0 <- 1
            }
        } else {
            # if there are robs = 0 but none of the double-sampled, use only robs = 1 information to estimate
            prob_0 <- 0
            pr_robs0 <- 0
            pr_robs1 <- 1
        }
    } else {
        prob_0 <- 0
        pr_robs0 <- 0
    }

    
    return(pr_robs1 * prob_1 + pr_robs0 * prob_0)
}


