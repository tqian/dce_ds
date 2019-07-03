###################################
# Make tables for simulation section
#      Tianchen Qian
#        2019.06.16
###################################


rm(list = ls())

nsim <- 1000
npara <- nsim

setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.09_simulation")
# setwd("~/git_deductive/simulation20190608/")
# setwd("~/deductive_simulation/")


library(dplyr)


# Collect simulation results ---------------------------------------------------------

design_simulation <- expand.grid(n = c(50, 100, 200, 500), dgm_type = c(1,6))

info_0 <- "z"
info_loss <- c("l")

for (case_id in 1:nrow(design_simulation)) {
    n <- design_simulation$n[case_id]
    dgm_type <- design_simulation$dgm_type[case_id]
    
    if (dgm_type == 1) {
        t0 <- 0.7
        tau_true <- 0.7698933
    } else if (dgm_type == 6) {
        t0 <- 0.7
        tau_true <- 0.589842  
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

performance_gathered <- performance_gathered[order(performance_gathered$estimator, performance_gathered$dgm_type, performance_gathered$n), ]
performance_gathered <- performance_gathered[, c(8,7,1:6)]
print(performance_gathered, n=Inf)

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



library(reshape)
library(kableExtra)
library(knitr)
library(xtable)



# Table 1: consistency ----------------------------------------------------

df1 <- filter(performance_gathered, estimator %in% c(1, 3, 7, 12, 13))[, c("dgm_type", "n", "estimator", "bias", "sd", "cp")]

df1$bias <- round(df1$bias * 100, 1)
df1$sd <- round(df1$sd * 100, 1)
df1$cp <- round(df1$cp * 100, 1)

table1_names <- c("dgm_type", "n", "estimator", "measure", "value")
table1 <- data.frame(matrix(NA, nrow = 3 * nrow(df1), ncol = length(table1_names))) # 3 is the number of different measurements
names(table1) <- table1_names
table1 <- as_tibble(table1)
measures <- c("bias", "sd", "cp")

i_table <- 1
for (irow in 1:nrow(df1)) {
    for (measure in measures) {
        table1[i_table, "dgm_type"] <- df1$dgm_type[irow]
        table1[i_table, "n"] <- df1$n[irow]
        table1[i_table, "estimator"] <- df1$estimator[irow]
        table1[i_table, "measure"] <- measure
        table1[i_table, "value"] <- df1[irow, measure]
        i_table <- i_table + 1
    }
}

casted <- cast(table1, dgm_type + n ~ estimator + measure, value = "value")

sink("table_generation/simulation_1.txt", append=FALSE)
mycaption <- "caption for simulation 1"
latex_code <- kable(casted, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
    add_header_above(c(" ", " ", "DE.Cox" = 3, "DE.LN" = 3, "PAR" = 3, "KM.S" = 3, "KM.C" = 3)) %>%
    # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
    # column_spec(1, bold=T) %>%
    collapse_rows(columns = 1, latex_hline = "major") %>%
    kable_styling(latex_options = c("scale_down"))
print(latex_code)
sink()





# Table 2: incorrect model for S ----------------------------------------------------

df2 <- filter(performance_gathered, estimator %in% c(1, 8, 3, 10))[, c("dgm_type", "n", "estimator", "bias", "sd", "cp")]

df2$bias <- round(df2$bias * 100, 1)
df2$sd <- round(df2$sd * 100, 1)
df2$cp <- round(df2$cp * 100, 1)

df2$estimator[which(df2$estimator == 1)] <- "1"
df2$estimator[which(df2$estimator == 3)] <- "3"
df2$estimator[which(df2$estimator == 8)] <- "1_1"
df2$estimator[which(df2$estimator == 10)] <- "3_1"

table2_names <- c("dgm_type", "n", "estimator", "measure", "value")
table2 <- data.frame(matrix(NA, nrow = 3 * nrow(df2), ncol = length(table2_names))) # 3 is the number of different measurements
names(table2) <- table2_names
table2 <- as_tibble(table2)
measures <- c("bias", "sd", "cp")

i_table <- 1
for (irow in 1:nrow(df2)) {
    for (measure in measures) {
        table2[i_table, "dgm_type"] <- df2$dgm_type[irow]
        table2[i_table, "n"] <- df2$n[irow]
        table2[i_table, "estimator"] <- df2$estimator[irow]
        table2[i_table, "measure"] <- measure
        table2[i_table, "value"] <- df2[irow, measure]
        i_table <- i_table + 1
    }
}

casted <- cast(table2, dgm_type + n ~ estimator + measure, value = "value")

sink("table_generation/simulation_2.txt", append=FALSE)
mycaption <- "caption for simulation 2"
latex_code <- kable(casted, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
    add_header_above(c(" ", " ", "DE.Cox" = 3, "DE.Cox.WrongS" = 3, "DE.LN" = 3, "DE.LN.WrongS" = 3)) %>%
    # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
    # column_spec(1, bold=T) %>%
    collapse_rows(columns = 1, latex_hline = "major") %>%
    kable_styling(latex_options = c("scale_down"))
print(latex_code)
sink()



# Table 3: impact of alpha ----------------------------------------------------

df3 <- filter(performance_gathered, estimator %in% c(1, 2, 3, 4))[, c("dgm_type", "n", "estimator", "bias", "sd", "cp")]

df3$bias <- round(df3$bias * 100, 1)
df3$sd <- round(df3$sd * 100, 1)
df3$cp <- round(df3$cp * 100, 1)

table3_names <- c("dgm_type", "n", "estimator", "measure", "value")
table3 <- data.frame(matrix(NA, nrow = 3 * nrow(df3), ncol = length(table3_names))) # 3 is the number of different measurements
names(table3) <- table3_names
table3 <- as_tibble(table3)
measures <- c("bias", "sd", "cp")

i_table <- 1
for (irow in 1:nrow(df3)) {
    for (measure in measures) {
        table3[i_table, "dgm_type"] <- df3$dgm_type[irow]
        table3[i_table, "n"] <- df3$n[irow]
        table3[i_table, "estimator"] <- df3$estimator[irow]
        table3[i_table, "measure"] <- measure
        table3[i_table, "value"] <- df3[irow, measure]
        i_table <- i_table + 1
    }
}

casted <- cast(table3, dgm_type + n ~ estimator + measure, value = "value")

sink("table_generation/simulation_3.txt", append=FALSE)
mycaption <- "caption for simulation 3"
latex_code <- kable(casted, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
    add_header_above(c(" ", " ", "DE.Cox" = 3, "DE.Cox(alpha=0)" = 3, "DE.LN" = 3, "DE.LN(alpha=0)" = 3)) %>%
    # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
    # column_spec(1, bold=T) %>%
    collapse_rows(columns = 1, latex_hline = "major") %>%
    kable_styling(latex_options = c("scale_down"))
print(latex_code)
sink()
