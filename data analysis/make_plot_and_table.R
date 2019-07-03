###################################
# Deductive estimation of survival probability in double sampling design
# Making plots
#      Tianchen Qian
#        2017.03.29
#        2017.04.07: remake plot for mortality curve with different cutoff (now cutoff = 1 year and 2 years)
#        2017.05.25: put two plots on the same page
#        2019.06.12: use the updated code with coxph and lognormal methods
###################################

rm(list = ls())

library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#F0E442", "#0072B2", "#CC79A7")

setwd("~/Dropbox/Research/git_deductive/double_sampling/2019.06.11_data_analysis_morogoro_n1773/")

DE_result <- read.csv("result/DE.csv", as.is = TRUE)
if (names(DE_result)[1]=="X") {
    DE_result$X <- NULL
}
KM_result <- read.csv("result/KmAllMethods.csv", as.is = TRUE)
if (names(KM_result)[1]=="X") {
    KM_result$X <- NULL
}

varnames <- c("estimator", "year", "gamma", "mortality", "left_ci", "right_ci")
df_plot <- data.frame(matrix(NA, nrow = nrow(DE_result) + nrow(KM_result), ncol = length(varnames)))
names(df_plot) <- varnames
irow <- 1

for (i in 1:nrow(DE_result)) {
    df_plot[irow, "estimator"] <- DE_result[i, "xdelta_working_model"]
    df_plot[irow, "year"] <- DE_result[i, "t0_year"]
    df_plot[irow, "gamma"] <- DE_result[i, "gamma"]
    mortality <- 1 - DE_result[i, "survProb"]
    left_ci <- mortality - 1.96 * DE_result[i, "se"]
    right_ci <- mortality + 1.96 * DE_result[i, "se"]
    df_plot[irow, "mortality"] <- mortality
    df_plot[irow, "left_ci"] <- left_ci
    df_plot[irow, "right_ci"] <- right_ci

    irow <- irow + 1
}

for (i in 1:nrow(KM_result)) {
    df_plot[irow, "estimator"] <- KM_result[i, "estimator"]
    df_plot[irow, "year"] <- KM_result[i, "t0_year"]
    df_plot[irow, "gamma"] <- KM_result[i, "gamma"]
    df_plot[irow, "mortality"] <- 1 - KM_result[i, "survProb"]
    df_plot[irow, "left_ci"] <- 1 - KM_result[i, "right_ci"]
    df_plot[irow, "right_ci"] <- 1 - KM_result[i, "left_ci"]
    
    irow <- irow + 1
}

df_plot$gamma[df_plot$gamma == 513.5] <- 514 # some are 513.5, some are 514


# some of the confidence intervals are swapped
for (irow in 1:nrow(df_plot)) {
    if (df_plot$left_ci[irow] > df_plot$right_ci[irow]) {
        cat(irow, "")
        tmp <- df_plot$left_ci[irow]
        df_plot$left_ci[irow] <- df_plot$right_ci[irow]
        df_plot$right_ci[irow] <- tmp
    }
}


# Plot 1: mortality curve with no restriction of double-sampling ----------

df_plot1 <- subset(df_plot, gamma == 1061)
estimators <- unique(df_plot1$estimator)
df_plot1 <- rbind(df_plot1,
                  data.frame(estimator = estimators,
                             year = rep(0, length(estimators)),
                             gamma = rep(1061, length(estimators)),
                             mortality = rep(0, length(estimators)),
                             left_ci = rep(0, length(estimators)),
                             right_ci = rep(0, length(estimators))))

df_plot1$mortality <- df_plot1$mortality * 100
df_plot1$left_ci <- df_plot1$left_ci * 100
df_plot1$right_ci <- df_plot1$right_ci * 100

df_plot1 <- subset(df_plot1, estimator %in% c("coxph", "lognormal", "KmCompleteCase", "KmStratifiedByRobs"))
df_plot1$estimator[df_plot1$estimator == "coxph"] <- 1
df_plot1$estimator[df_plot1$estimator == "lognormal"] <- 2
df_plot1$estimator[df_plot1$estimator == "KmStratifiedByRobs"] <- 3
df_plot1$estimator[df_plot1$estimator == "KmCompleteCase"] <- 4


ggplot(df_plot1) +
    geom_line(aes(x = year, y = mortality, color = estimator), size = 1.5) +
    xlab("Time on ART (years)") +
    ylab("Mortality rate (%)") +
    scale_color_manual(name = "Estimator",
                       labels = c("DE.Cox", "DE.LN", "KM.S", "KM.C"),
                       values = cbbPalette[1:4]) +
    geom_line(aes(x = year, y = left_ci, color = estimator),
              linetype = 2) +
    geom_line(aes(x = year, y = right_ci, color = estimator),
              linetype = 2) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave("plot/morogoro1733-mortality-DE,KM.pdf", width = 4, height = 5)

# ggplot(df_plot1) +
#     geom_line(aes(x = year, y = mortality, color = estimator, linetype = estimator), size = 2) +
#     xlab("Time on ART (years)") +
#     ylab("Mortality rate (%)") +
#     geom_line(aes(x = year, y = left_ci, color = estimator, linetype = estimator)) +
#     geom_line(aes(x = year, y = right_ci, color = estimator, linetype = estimator)) +
#     scale_color_manual(name = "Estimator",
#                        labels = c("DE.Cox", "DE.LN", "KM.S", "KM.C"),
#                        values = cbbPalette[1:4]) +
#     scale_linetype_manual(name = "Estimator",
#                        labels = c("DE.Cox", "DE.LN", "KM.S", "KM.C"),
#                        values = c(1,2,4,3)) +
#     coord_cartesian(ylim = c(0, 30)) +
#     theme_bw() +
#     theme(legend.position="bottom")



# Plot 2: mortality curve with different cutoffs (gamma) ----------

df_plot2 <- df_plot
estimators <- unique(df_plot2$estimator)
# unique(df_plot2$gamma)
gammas_to_use <- c(1061, 730, 547, 365)
df_plot2 <- subset(df_plot2, gamma %in% gammas_to_use)
df_plot2 <- rbind(df_plot2,
                  data.frame(estimator = rep(estimators, length(gammas_to_use)),
                             year = rep(0, length(estimators) * length(gammas_to_use)),
                             gamma = rep(gammas_to_use, each = length(estimators)),
                             mortality = rep(0, length(estimators) * length(gammas_to_use)),
                             left_ci = rep(0, length(estimators) * length(gammas_to_use)),
                             right_ci = rep(0, length(estimators) * length(gammas_to_use))))

df_plot2$gamma <- as.character(df_plot2$gamma)

df_plot2$mortality <- df_plot2$mortality * 100
df_plot2$left_ci <- df_plot2$left_ci * 100
df_plot2$right_ci <- df_plot2$right_ci * 100

# df_plot2 <- subset(df_plot2, estimator %in% c("coxph", "lognormal", "KmCompleteCase", "KmStratifiedByRobs"))
# df_plot2$estimator[df_plot2$estimator == "coxph"] <- 1
# df_plot2$estimator[df_plot2$estimator == "lognormal"] <- 2
# df_plot2$estimator[df_plot2$estimator == "KmStratifiedByRobs"] <- 3
# df_plot2$estimator[df_plot2$estimator == "KmCompleteCase"] <- 4

df_plot2$gamma_plot[df_plot2$gamma == 1061] <- 1
df_plot2$gamma_plot[df_plot2$gamma == 730] <- 2
df_plot2$gamma_plot[df_plot2$gamma == 547] <- 3
df_plot2$gamma_plot[df_plot2$gamma == 365] <- 4
df_plot2$gamma_plot <- as.character(df_plot2$gamma_plot)



ggplot(subset(df_plot2, estimator %in% c("coxph"))) +
    geom_line(aes(x = year, y = mortality, color = gamma_plot), size = 1.5) +
    xlab("Time on ART (years)") +
    ylab("Mortality rate (%) estimated by DE.Cox") +
    scale_color_manual(name = expression(gamma),
                       labels = c(expression(infinity), "2", "1.5", "1"),
                       # values = cbbPalette[1:4]) +
                       values = cbbPalette[1:4]) +
    geom_line(aes(x = year, y = left_ci, color = gamma_plot), linetype = 2) +
    geom_line(aes(x = year, y = right_ci, color = gamma_plot), linetype = 2) +
    # geom_ribbon(aes(x = year, ymin = left_ci, ymax = right_ci, fill = gamma_plot), alpha = 0.2) +
    # scale_fill_manual(name = expression(gamma),
    #                   labels = c(expression(infinity), "2", "1.5", "1"),
    #                   values = cbbPalette[1:4]) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave("plot/morogoro1733-mortality-DECox-gamma.pdf", width = 4, height = 5)

ggplot(subset(df_plot2, estimator %in% c("lognormal"))) +
    geom_line(aes(x = year, y = mortality, color = gamma_plot), size = 1.5) +
    xlab("Time on ART (years)") +
    ylab("Mortality rate (%) estimated by DE.LN") +
    scale_color_manual(name = expression(gamma),
                       labels = c(expression(infinity), "2", "1.5", "1"),
                       values = cbbPalette[1:4]) +
    geom_line(aes(x = year, y = left_ci, color = gamma_plot), linetype = 2) +
    geom_line(aes(x = year, y = right_ci, color = gamma_plot), linetype = 2) +
    # geom_ribbon(aes(x = year, ymin = left_ci, ymax = right_ci, fill = gamma_plot), alpha = 0.2) +
    # scale_fill_manual(name = expression(gamma),
    #                   labels = c(expression(infinity), "2", "1.5", "1"),
    #                   values = cbbPalette[1:4]) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave("plot/morogoro1733-mortality-DELN-gamma.pdf", width = 4, height = 5)

ggplot(subset(df_plot2, estimator %in% c("KmStratifiedByRobs"))) +
    geom_line(aes(x = year, y = mortality, color = gamma_plot), size = 1.5) +
    xlab("Time on ART (years)") +
    ylab("Mortality rate (%) estimated by KM.S") +
    scale_color_manual(name = expression(gamma),
                       labels = c(expression(infinity), "2", "1.5", "1"),
                       values = cbbPalette[1:4]) +
    geom_line(aes(x = year, y = left_ci, color = gamma_plot), linetype = 2) +
    geom_line(aes(x = year, y = right_ci, color = gamma_plot), linetype = 2) +
    # geom_ribbon(aes(x = year, ymin = left_ci, ymax = right_ci, fill = gamma_plot), alpha = 0.2) +
    # scale_fill_manual(name = expression(gamma),
    #                   labels = c(expression(infinity), "2", "1.5", "1"),
    #                   values = cbbPalette[1:4]) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave("plot/morogoro1733-mortality-KMS-gamma.pdf", width = 4, height = 5)

ggplot(subset(df_plot2, estimator %in% c("KmCompleteCase"))) +
    geom_line(aes(x = year, y = mortality, color = gamma_plot), size = 1.5) +
    xlab("Time on ART (years)") +
    ylab("Mortality rate (%) estimated by KM.C") +
    scale_color_manual(name = expression(gamma),
                       labels = c(expression(infinity), "2", "1.5", "1"),
                       values = cbbPalette[1:4]) +
    geom_line(aes(x = year, y = left_ci, color = gamma_plot), linetype = 2) +
    geom_line(aes(x = year, y = right_ci, color = gamma_plot), linetype = 2) +
    # geom_ribbon(aes(x = year, ymin = left_ci, ymax = right_ci, fill = gamma_plot), alpha = 0.2) +
    # scale_fill_manual(name = expression(gamma),
    #                   labels = c(expression(infinity), "2", "1.5", "1"),
    #                   values = cbbPalette[1:4]) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave("plot/morogoro1733-mortality-KMC-gamma.pdf", width = 4, height = 5)




# ggplot(df_plot2) +
#     geom_line(aes(x = year, y = mortality, color = gamma_plot), size = 1.5) +
#     xlab("Time on ART (years)") +
#     ylab("Mortality rate (%)") +
#     scale_color_manual(name = expression(gamma),
#                        labels = c(expression(infinity), "2", "1.5", "1"),
#                        values = cbbPalette[1:4]) +
#     geom_line(aes(x = year, y = left_ci, color = gamma_plot), 
#               linetype = 2) +
#     geom_line(aes(x = year, y = right_ci, color = gamma_plot), 
#               linetype = 2) +
#     theme_bw()

# ggplot(df_plot2) +
#     geom_line(aes(x = year, y = mortality, color = gamma_plot, linetype = gamma_plot), size = 1.5) +
#     geom_line(aes(x = year, y = left_ci, color = gamma_plot, linetype = gamma_plot)) +
#     geom_line(aes(x = year, y = right_ci, color = gamma_plot, linetype = gamma_plot)) +
#     xlab("Time on ART (years)") +
#     ylab("Mortality rate (%)") +
#     scale_color_manual(name = expression(gamma),
#                        labels = c(expression(infinity), "2", "1.5", "1"),
#                        values = cbbPalette[1:4]) +
#     scale_linetype_manual(name = expression(gamma),
#                           labels = c(expression(infinity), "2", "1.5", "1"),
#                           values = c(1,2,4,3)) +
#     coord_cartesian(ylim = c(0, 30)) +
#     theme_bw() +
#     theme(legend.position="bottom")





# Table 1: estimated mortality rates with various gamma ----------------------------------------

library(reshape)
library(kableExtra)
library(knitr)
library(xtable)
library(dplyr)

# df_table1 <- subset(df_plot, gamma == 1061)
df_table1 <- df_plot
gammas_to_use <- c(1061, 730, 547, 365)
df_table1 <- subset(df_table1, year %in% c(0.5, 1, 1.5, 2) & 
                        estimator %in% c("coxph", "lognormal", "KmStratifiedByRobs", "KmCompleteCase") &
                        gamma %in% gammas_to_use)
df_table1$mortality <- sprintf(fmt = "%.1f", df_table1$mortality * 100)
df_table1$left_ci <- sprintf(fmt = "%.1f", df_table1$left_ci * 100)
df_table1$right_ci <- sprintf(fmt = "%.1f", df_table1$right_ci * 100)
df_table1$ci <- paste0("(", df_table1$left_ci, ", ", df_table1$right_ci, ")")


df_table1$estimator[df_table1$estimator == "coxph"] <- "est1"
df_table1$estimator[df_table1$estimator == "lognormal"] <- "est2"
df_table1$estimator[df_table1$estimator == "KmStratifiedByRobs"] <- "est3"
df_table1$estimator[df_table1$estimator == "KmCompleteCase"] <- "est4"

table1_names <- c("gamma", "estimator", "year", "measure", "value")
table1 <- data.frame(matrix(NA, nrow = 2 * nrow(df_table1), ncol = length(table1_names))) # 2 is the number of different measurements
names(table1) <- table1_names
measures <- c("mortality", "ci")

i_table <- 1
for (irow in 1:nrow(df_table1)) {
    for (measure in measures) {
        table1[i_table, "gamma"] <- df_table1$gamma[irow]
        table1[i_table, "estimator"] <- df_table1$estimator[irow]
        table1[i_table, "year"] <- df_table1$year[irow]
        table1[i_table, "measure"] <- measure
        table1[i_table, "value"] <- df_table1[irow, measure]
        i_table <- i_table + 1
    }
}

table1$measure[table1$measure == "mortality"] <- "1.mortality"
table1$measure[table1$measure == "ci"] <- "2.ci"

table1$gamma <- as.character(table1$gamma)
table1$gamma[table1$gamma == "1061"] <- "1.1061"
table1$gamma[table1$gamma == "730"] <- "2.730"
table1$gamma[table1$gamma == "547"] <- "3.547"
table1$gamma[table1$gamma == "365"] <- "4.365"

casted <- cast(table1, gamma + year ~ estimator + measure, value = "value")

sink("table_generation/data_analysis_1.txt", append=FALSE)
mycaption <- "caption for data analysis 1"
latex_code <- kable(casted, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
    add_header_above(c(" ", " ", "DE.Cox" = 2, "DE.LN" = 2, "KM.S" = 2, "KM.C" = 2)) %>%
    # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
    # column_spec(1, bold=T) %>%
    collapse_rows(columns = 1, latex_hline = "major") %>%
    kable_styling(latex_options = c("scale_down"))
print(latex_code)
sink()















