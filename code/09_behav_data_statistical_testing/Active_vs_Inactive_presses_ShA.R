
library(here)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(cowplot)
library(sessioninfo)


####################   Statistical testing for behavioral data    ######################


## ShA behavioral data:
#  - Number of Active and Inactive lever presses x rat x ShA session
Short_Acess_Active_LPs <- read.delim(file = here("raw-data/raw_behavioral_data/Short_Access_Active_LPs.txt"))
colnames(Short_Acess_Active_LPs) <- c("session", "hr",
                                      paste0("rat_Saline_Inactive_", 1:11),
                                      paste0("rat_Saline_Active_", 1:11),
                                      paste0("rat_Fentanyl_Inactive_", 1:11),
                                      paste0("rat_Fentanyl_Active_", 1:11))

Short_Acess_Active_LPs[, c("rat_Fentanyl_Inactive_9", "rat_Fentanyl_Inactive_10", "rat_Fentanyl_Inactive_11",
                           "rat_Fentanyl_Active_9", "rat_Fentanyl_Active_10", "rat_Fentanyl_Active_11")] <- NULL


#______________________________________________________________________________
#                 1. Number of Active vs Inactive lever presses
#______________________________________________________________________________
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
#                            Short Access sessions
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                   All rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
Short_Acess_Active_LPs_longer <-  pivot_longer(Short_Acess_Active_LPs, setdiff(colnames(Short_Acess_Active_LPs), c("session", "hr")))
colnames(Short_Acess_Active_LPs_longer) <- c("session", "hr", "rat", "number")
Short_Acess_Active_LPs_longer$substance <- str_split_i(Short_Acess_Active_LPs_longer$rat, "_", 2)
Short_Acess_Active_LPs_longer$type <- str_split_i(Short_Acess_Active_LPs_longer$rat, "_", 3)
Short_Acess_Active_LPs_longer$rat <- str_split_i(Short_Acess_Active_LPs_longer$rat, "_", 4)

Short_Acess_Active_LPs_longer$type <- factor(Short_Acess_Active_LPs_longer$type, levels = c("Inactive", "Active"))
Short_Acess_Active_LPs_longer$session <- Short_Acess_Active_LPs_longer$session %>% replace(., . == "Sha4", "ShA4")

## Two-sample, two-tailed paired t-test
# ------------------------------------------------------------------------------
# Used to determine whether the mean num of Active vs Inactive lever presses is
# significantly different, with paired rat measurements across lever types.
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality of the differences: not satisfied (substance and session number must be considered)
## Note: observations are not independent within nor between lever types

differences_x_fent_rat_x_session <- Short_Acess_Active_LPs[, paste0("rat_Fentanyl_Inactive_", 1:8)] - Short_Acess_Active_LPs[, paste0("rat_Fentanyl_Active_", 1:8)]
differences_x_saline_rat_x_session <- Short_Acess_Active_LPs[, paste0("rat_Saline_Inactive_", 1:11)] - Short_Acess_Active_LPs[, paste0("rat_Saline_Active_", 1:11)]
differences_x_rat_x_session_all <- c(unlist(differences_x_fent_rat_x_session), unlist(differences_x_saline_rat_x_session)) %>% as.data.frame()

shapiro.test(differences_x_rat_x_session_all$.)

# Shapiro-Wilk normality test
#
# data:  differences_x_rat_x_session_all$.
# W = 0.24577, p-value < 2.2e-16

qqplot <- ggplot(data = differences_x_rat_x_session_all, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_all_rats.pdf", height = 8, width = 8)
# ------------------------------------------------------------------------------

## Plot:
ggplot(Short_Acess_Active_LPs_longer,
       aes(x = type, y = number,
           colour = substance)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width = 0.2, color = 'gray50') +
    scale_colour_manual(values = c('Fentanyl'='turquoise3', 'Saline'='yellow3')) +
    labs(y = "Number of presses", x = "Lever type", title = "Fentanyl + Saline rats", color = "Substance") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_all_rats.pdf", height = 5, width = 6)

## Analysis of variance:
lm_lever_presses <- lm(number ~ type,
                       data = Short_Acess_Active_LPs_longer)
anova(lm_lever_presses)
# Analysis of Variance Table
# Response: number
#            Df Sum Sq Mean Sq F value  Pr(>F)
# type        1   7130  7129.9    4.79 0.02965 *
# Residuals 226 336401  1488.5
# ---

## t-test
t_test(number ~ type,
       alternative = "two.sided",
       paired = TRUE,
       data = Short_Acess_Active_LPs_longer, ref.group = "Active")
# .y.    group1   group2    n1    n2 statistic    df      p        (F-stat is t-stat**2)
#  <chr>  <chr>    <chr>  <int> <int>     <dbl> <dbl>  <dbl>
#  number Active Inactive   114   114      2.21   113 0.0294


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                Fentanyl rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
fentanyl_rats_lever_presses <- Short_Acess_Active_LPs[, c("session",
                                                        paste0("rat_Fentanyl_Inactive_", 1:8),
                                                        paste0("rat_Fentanyl_Active_", 1:8))]

fentanyl_rats_lever_presses_longer <- pivot_longer(fentanyl_rats_lever_presses, colnames(fentanyl_rats_lever_presses)[-1])
colnames(fentanyl_rats_lever_presses_longer) <- c("session", "rat", "number")
fentanyl_rats_lever_presses_longer$type <- str_split_i(fentanyl_rats_lever_presses_longer$rat, "_", 3)
fentanyl_rats_lever_presses_longer$rat <- str_split_i(fentanyl_rats_lever_presses_longer$rat, "_", 4)

fentanyl_rats_lever_presses_longer$type <- factor(fentanyl_rats_lever_presses_longer$type, levels = c("Inactive", "Active"))
fentanyl_rats_lever_presses_longer$session <- fentanyl_rats_lever_presses_longer$session %>% replace(., . == "Sha4", "ShA4")


# * * * * * * * * * * * *  Differences by lever type  * * * * * * * * * * * * *

## Two-sample, two-tailed paired t-test
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality of the differences: not satisfied (session number must be considered)
## Note: observations are not independent within nor between lever types

differences_x_rat_x_session <- fentanyl_rats_lever_presses[, paste0("rat_Fentanyl_Inactive_", 1:8)] - fentanyl_rats_lever_presses[, paste0("rat_Fentanyl_Active_", 1:8)]
differences_x_rat_x_session_all <- unlist(differences_x_rat_x_session) %>% as.data.frame()

shapiro.test(differences_x_rat_x_session_all$.)
# Shapiro-Wilk normality test
#
# data:  differences_x_rat_x_session_all$.
# W = 0.3362, p-value = 2.047e-13

qqplot <- ggplot(data = differences_x_rat_x_session_all, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_fentanyl_rats.pdf", height = 8, width = 8)
# ------------------------------------------------------------------------------

## Plot:
ggplot(fentanyl_rats_lever_presses_longer,
       aes(x = type, y = number,
           colour = type)) +
    geom_jitter(width = 0.05, alpha = 0.8) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width=0.2, color='gray50') +
    guides(color = "none") +
    scale_colour_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Lever type", title = "Fentanyl rats") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_fentanyl_rats.pdf", height = 5, width = 6)

## Test:
t_test(number ~ type,
       alternative = "two.sided",
       paired = TRUE,
       data = fentanyl_rats_lever_presses_longer, ref.group = "Active")
#  .y.    group1   group2    n1    n2 statistic    df      p
#  <chr>  <chr>    <chr>  <int> <int>     <dbl> <dbl>  <dbl>
#  number Inactive Active    48    48      2.17    47 0.0348

lm_lever_presses <- lm(number ~ type, fentanyl_rats_lever_presses_longer)
anova(lm_lever_presses)
# Analysis of Variance Table
#
# Response: number
#           Df Sum Sq Mean Sq F value  Pr(>F)
# type       1  15632 15631.5  4.6445 0.03371 *
# Residuals 94 316365  3365.6


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                Saline rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
saline_rats_lever_presses <- Short_Acess_Active_LPs[, c("session",
                                                          paste0("rat_Saline_Inactive_", 1:11),
                                                          paste0("rat_Saline_Active_", 1:11))]

saline_rats_lever_presses_longer <- pivot_longer(saline_rats_lever_presses, colnames(saline_rats_lever_presses)[-1])
colnames(saline_rats_lever_presses_longer) <- c("session", "rat", "number")
saline_rats_lever_presses_longer$type <- str_split_i(saline_rats_lever_presses_longer$rat, "_", 3)
saline_rats_lever_presses_longer$rat <- str_split_i(saline_rats_lever_presses_longer$rat, "_", 4)

saline_rats_lever_presses_longer$type <- factor(saline_rats_lever_presses_longer$type, levels = c("Inactive", "Active"))
saline_rats_lever_presses_longer$session <- saline_rats_lever_presses_longer$session %>% replace(., . == "Sha4", "ShA4")


# * * * * * * * * * * * *  Differences by lever type  * * * * * * * * * * * * *

## Two-sample, two-tailed paired t-test
# ------------------------------------------------------------------------------
#  A) Normality of the differences: not satisfied (session number must be considered)
## Note: observations are not independent within nor between lever types

differences_x_rat_x_session <- saline_rats_lever_presses[, paste0("rat_Saline_Inactive_", 1:11)] - saline_rats_lever_presses[, paste0("rat_Saline_Active_", 1:11)]
differences_x_rat_x_session_all <- unlist(differences_x_rat_x_session) %>% as.data.frame()

shapiro.test(differences_x_rat_x_session_all$.)
# Shapiro-Wilk normality test
#
# data:  differences_x_rat_x_session_all$.
# W = 0.93552, p-value = 0.001961

qqplot <- ggplot(data = differences_x_rat_x_session_all, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_saline_rats.pdf", height = 8, width = 8)
# ------------------------------------------------------------------------------

## Plot:
ggplot(saline_rats_lever_presses_longer,
       aes(x = type, y = number,
           colour = type)) +
    geom_jitter(width = 0.05, alpha = 0.8) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width=0.2, color='gray50') +
    guides(color = "none") +
    scale_colour_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Lever type", title = "Saline rats") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_saline_rats.pdf", height = 5, width = 6)

## Test:
t_test(number ~ type,
       alternative = "two.sided",
       paired = TRUE,
       data = saline_rats_lever_presses_longer, ref.group = "Active")
#  .y.    group1   group2    n1    n2 statistic    df      p
#  <chr>  <chr>    <chr>  <int> <int>     <dbl> <dbl>  <dbl>
#   number Active Inactive    66    66     0.951    65 0.345

lm_lever_presses <- lm(number ~ type, saline_rats_lever_presses_longer)
anova(lm_lever_presses)
# Analysis of Variance Table
#
# Response: number
#            Df Sum Sq Mean Sq F value Pr(>F)
# type        1   18.9  18.939  0.6342 0.4273
# Residuals 130 3882.4  29.865



#______________________________________________________________________________
#                 2. Number of lever presses across sessions
#______________________________________________________________________________
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
#                            Short Access sessions
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                   All rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## One-way ANOVA test
# ------------------------------------------------------------------------------
# Used to determine whether the mean nums of lever presses are significantly
# different across sessions.
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality of the residuals: not satisfied (substance and lever type must be considered)
lm_lever_presses <- lm(number ~ session,
                       data = Short_Acess_Active_LPs_longer)
residuals_lever_presses <- residuals(lm_lever_presses)
shapiro.test(residuals_lever_presses)
#        Shapiro-Wilk normality test
#
# data:  residuals_lever_presses
# W = 0.20263, p-value < 2.2e-16

residuals_lever_presses <- residuals_lever_presses %>% as.data.frame()
qqplot <- ggplot(data = residuals_lever_presses, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_residuals_all_rats.pdf", height = 8, width = 8)


#  B) Variance equality (Levene test for non-normal data: 31% confident about it):
levene_test(data = Short_Acess_Active_LPs_longer,
            formula = number ~ session)
#   df1   df2 statistic      p
# <int> <int>     <dbl>  <dbl>
#     5   222      1.19  0.316

#  C) Independent observations: not satisfied
# ------------------------------------------------------------------------------

## Plot:
ggplot(Short_Acess_Active_LPs_longer,
       aes(x = session, y = number,
           colour = substance)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width = 0.2, color = 'gray50') +
    scale_colour_manual(values = c('Fentanyl'='turquoise3', 'Saline'='yellow3')) +
    labs(y = "Number of presses (Active + Inactive)", x = "Session", title = "Fentanyl + Saline rats", color = "Substance") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/sessions_all_rats.pdf", height = 5, width = 6)


## Analysis of variance:
lm_session <- lm(number ~ session,
                       data = Short_Acess_Active_LPs_longer)
anova(lm_session)
# Analysis of Variance Table
# Response: number
#            Df Sum Sq Mean Sq F value Pr(>F)
# session     5   9741  1948.3  1.2958 0.2667
# Residuals 222 333789  1503.6
# ---
anova_test(Short_Acess_Active_LPs_longer, dv = number, between = session)
#    ANOVA Table (type II tests)
#
#   Effect DFn DFd     F     p p<.05   ges
#  session   5 222 1.296 0.267       0.028



## One-way ANOVA paired test
# ------------------------------------------------------------------------------
# Used to determine whether the mean nums of lever presses are significantly
# different across sessions, with measurements in each session coming from the
# same rats (and lever type).
# ------------------------------------------------------------------------------
## Assumptions:
# - Normality: not satisfied for any session (substance and lever type must be considered)
Short_Acess_Active_LPs_longer %>%
    group_by(session) %>%
    shapiro_test(number)
# session    variable statistic        p
#    ShA1    number       0.812 1.82e- 5
#    ShA2    number       0.851 1.38e- 4
#    ShA3    number       0.856 1.79e- 4
#    ShA4    number       0.918 8.78e- 3
#    ShA5    number       0.290 2.39e-12
#    ShA6    number       0.252 1.10e-12
# ------------------------------------------------------------------------------

## Analysis of variance (repeated measures given by rat x substance x lever type):
Short_Acess_Active_LPs_longer$rat_sub_type <- paste(Short_Acess_Active_LPs_longer$rat,
                                                    Short_Acess_Active_LPs_longer$substance,
                                                    Short_Acess_Active_LPs_longer$type, sep = "_")
anova_test(Short_Acess_Active_LPs_longer, dv = number, within = session, wid = rat_sub_type)
#        ANOVA Table (type III tests)
#
# $ANOVA
#    Effect DFn DFd     F     p p<.05   ges
#   session   5 185 1.614 0.158       0.028
#
# $`Mauchly's Test for Sphericity`
#     Effect        W         p p<.05
#   session 3.44e-07 9.83e-102     *
#
# $`Sphericity Corrections`
#    Effect   GGe      DF[GG] p[GG] p[GG]<.05   HFe      DF[HF] p[HF] p[HF]<.05
#   session 0.206 1.03, 38.05 0.212           0.206 1.03, 38.14 0.212
# ------------------------------------------------------------------------------


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                Fentanyl rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## One-way ANOVA test
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality of the residuals: not satisfied (lever type must be considered)
lm_lever_presses <- lm(number ~ session,
                       data = fentanyl_rats_lever_presses_longer)
residuals_lever_presses <- residuals(lm_lever_presses)
shapiro.test(residuals_lever_presses)
#        Shapiro-Wilk normality test
#
# data:  residuals_lever_presses
# W = 0.35035, p-value < 2.2e-16

residuals_lever_presses <- residuals_lever_presses %>% as.data.frame()
qqplot <- ggplot(data = residuals_lever_presses, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_residuals_fentanyl_rats.pdf", height = 8, width = 8)


#  B) Variance equality (Levene test for non-normal data: 30% confident about it):
levene_test(data = fentanyl_rats_lever_presses_longer,
            formula = number ~ session)
#   df1   df2 statistic      p
# <int> <int>     <dbl>  <dbl>
#     5    90      1.21 0.309

#  C) Independent observations: not satisfied
# ------------------------------------------------------------------------------

## Plot:
ggplot(fentanyl_rats_lever_presses_longer,
       aes(x = session, y = number,
           colour = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width = 0.2, color = 'gray50') +
    scale_colour_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses (Active + Inactive)", x = "Session", title = "Fentanyl", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/sessions_fentanyl_rats.pdf", height = 5, width = 6)


## Analysis of variance:
lm_session <- lm(number ~ session,
                 data = fentanyl_rats_lever_presses_longer)
anova(lm_session)
# Analysis of Variance Table
# Response: number
#           Df Sum Sq Mean Sq F value Pr(>F)
# session    5  21315    4263  1.2349 0.2994
# Residuals 90 310682    3452
# ---


## One-way ANOVA paired test
# ------------------------------------------------------------------------------
## Assumptions:
# - Normality: only satisfied for session 4.
fentanyl_rats_lever_presses_longer %>%
    group_by(session) %>%
    shapiro_test(number)
# session    variable statistic        p
#    ShA1    number       0.812 0.00394
#    ShA2    number       0.835 0.00839
#    ShA3    number       0.829 0.00682
#    ShA4    number       0.942 0.377
#    ShA5    number       0.411 0.000000421
#    ShA6    number       0.387 0.000000280
# ------------------------------------------------------------------------------

## Analysis of variance (repeated measures given by rat x lever type):
fentanyl_rats_lever_presses_longer$rat_type <- paste(fentanyl_rats_lever_presses_longer$rat,
                                                     fentanyl_rats_lever_presses_longer$type, sep = "_")
anova_test(fentanyl_rats_lever_presses_longer, dv = number, within = session, wid = rat_type)
#       ANOVA Table (type III tests)
#
# $ANOVA
#    Effect DFn DFd     F     p p<.05   ges
# 1 session   5  75 1.525 0.192       0.064
#
# $`Mauchly's Test for Sphericity`
#    Effect        W        p p<.05
# 1 session 2.57e-08 4.98e-40     *
#
#     $`Sphericity Corrections`
#    Effect   GGe      DF[GG] p[GG] p[GG]<.05   HFe      DF[HF] p[HF] p[HF]<.05
# 1 session 0.203 1.02, 15.24 0.236           0.204 1.02, 15.29 0.236
# ------------------------------------------------------------------------------


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                Saline rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## One-way ANOVA test
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality of the residuals: not satisfied (lever type must be considered)
lm_lever_presses <- lm(number ~ session,
                       data = saline_rats_lever_presses_longer)
residuals_lever_presses <- residuals(lm_lever_presses)
shapiro.test(residuals_lever_presses)
#        Shapiro-Wilk normality test
#
# data:  residuals_lever_presses
# W = 0.87017, p-value = 2.201e-09

residuals_lever_presses <- residuals_lever_presses %>% as.data.frame()
qqplot <- ggplot(data = residuals_lever_presses, aes(sample = .)) +
    stat_qq() +
    stat_qq_line(colour = "blue")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/qqplot_residuals_saline_rats.pdf", height = 8, width = 8)


#  B) Variance equality (Levene test for non-normal data: 90% confident about it):
levene_test(data = saline_rats_lever_presses_longer,
            formula = number ~ session)
#   df1   df2 statistic      p
# <int> <int>     <dbl>  <dbl>
#    5   126     0.306  0.909

#  C) Independent observations: not satisfied
# ------------------------------------------------------------------------------

## Plot:
ggplot(saline_rats_lever_presses_longer,
       aes(x = session, y = number,
           colour = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    geom_boxplot(alpha = 0, size = 0.6, width = 0.2, color = 'gray50') +
    scale_colour_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses (Active + Inactive)", x = "Session", title = "Saline", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/sessions_saline_rats.pdf", height = 5, width = 6)


## Analysis of variance:
lm_session <- lm(number ~ session,
                 data = saline_rats_lever_presses_longer)
anova(lm_session)
# Analysis of Variance Table
# Response: number
#            Df Sum Sq Mean Sq F value Pr(>F)
# session     5   34.0  6.7939  0.2213 0.9528
# Residuals 126 3867.4 30.6934
# ---


## One-way ANOVA paired test
# ------------------------------------------------------------------------------
## Assumptions:
# - Normality: not satisfied for any session.
saline_rats_lever_presses_longer %>%
    group_by(session) %>%
    shapiro_test(number)
# session    variable statistic        p
#   ShA1    number       0.770 0.000178
#   ShA2    number       0.878 0.0111
#   ShA3    number       0.891 0.0199
#   ShA4    number       0.896 0.0245
#   ShA5    number       0.892 0.0202
#   ShA6    number       0.743 0.0000725
# ------------------------------------------------------------------------------

## Analysis of variance (repeated measures given by rat x lever type):
saline_rats_lever_presses_longer$rat_type <- paste(saline_rats_lever_presses_longer$rat,
                                                   saline_rats_lever_presses_longer$type, sep = "_")
anova_test(saline_rats_lever_presses_longer, dv = number, within = session, wid = rat_type)
#       ANOVA Table (type III tests)
#
# $ANOVA
#    Effect DFn DFd     F     p p<.05   ges
# 1 session   5 105 0.362 0.874       0.009
#
# $`Mauchly's Test for Sphericity`
#    Effect     W     p p<.05
# 1 session 0.177 0.003     *
#
# $`Sphericity Corrections`
#    Effect   GGe      DF[GG] p[GG] p[GG]<.05  HFe      DF[HF] p[HF] p[HF]<.05
# 1 session 0.613 3.06, 64.36 0.785           0.73 3.65, 76.61 0.818
# ------------------------------------------------------------------------------



#______________________________________________________________________________
#         3. Number of Active vs Inactive lever presses across sessions
#______________________________________________________________________________
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
#                            Short Access sessions
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                   All rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## Three-way mixed ANOVA test                        <- *** Correct test ***
# ------------------------------------------------------------------------------
# Used to investigate if there are differences in the num of lever presses by
# lever type over time (sessions), i.e. if there's a signif interaction between
# lever type and session on the lever presses, differing between fentanyl and
# saline rats (between-subjects variable). This, with repeated measures of the
# same rats across lever types and sessions (within-subjects variables).
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality: not satisfied for most samples.

Short_Acess_Active_LPs_longer %>%
    group_by(type, session, substance) %>%
    shapiro_test(number)
#    session substance     type variable statistic            p
#       ShA1  Fentanyl Inactive   number 0.6018055 1.669266e-04
#       ShA1    Saline Inactive   number 0.9017693 1.942591e-01
#       ShA2  Fentanyl Inactive   number 0.8078586 3.472143e-02
#       ShA2    Saline Inactive   number 0.8284270 2.226057e-02
#       ShA3  Fentanyl Inactive   number 0.8766270 1.747996e-01
#       ShA3    Saline Inactive   number 0.8565220 5.189650e-02
#       ShA4  Fentanyl Inactive   number 0.8871411 2.200697e-01
#       ShA4    Saline Inactive   number 0.7383837 1.467053e-03
#       ShA5  Fentanyl Inactive   number 0.8121579 3.856780e-02
#       ShA5    Saline Inactive   number 0.8539725 4.808088e-02
#       ShA6  Fentanyl Inactive   number 0.7653255 1.204717e-02
#       ShA6    Saline Inactive   number 0.7899918 6.947008e-03
#       ShA1  Fentanyl   Active   number 0.8962911 2.674433e-01
#       ShA1    Saline   Active   number 0.8785557 9.973906e-02
#       ShA2  Fentanyl   Active   number 0.9089372 3.466401e-01
#       ShA2    Saline   Active   number 0.8855543 1.223119e-01
#       ShA3  Fentanyl   Active   number 0.8583053 1.154880e-01
#       ShA3    Saline   Active   number 0.9184426 3.058812e-01
#       ShA4  Fentanyl   Active   number 0.9108799 3.603014e-01
#       ShA4    Saline   Active   number 0.9064510 2.212623e-01
#       ShA5  Fentanyl   Active   number 0.5507110 4.173957e-05
#       ShA5    Saline   Active   number 0.9679332 8.649351e-01
#       ShA6  Fentanyl   Active   number 0.5281744 2.249989e-05
#       ShA6    Saline   Active   number 0.7639903 3.166671e-03
# ------------------------------------------------------------------------------

## Differences by lever type across sessions:
p1 <- ggplot(Short_Acess_Active_LPs_longer,
       aes(x = session, y = number)) +
    theme_light() +
    geom_boxplot(aes(color = type), size = 0.5, outlier.alpha = 0, show.legend = F) +
    scale_colour_manual(values = c("Inactive" = "gray50", "Active" = "gray50")) +
    ggnewscale::new_scale(new_aes = "color") +
    geom_point(aes(color = type), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl + Saline rats", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

## Differences by lever type across sessions in fentanyl and saline separately:
p1_sub <- ggplot(Short_Acess_Active_LPs_longer,
             aes(x = session, y = number)) +
    theme_light() +
    geom_boxplot(aes(color = type), size = 0.5, outlier.alpha = 0, show.legend = F) +
    scale_colour_manual(values = c("Inactive" = "gray50", "Active" = "gray50")) +
    ggnewscale::new_scale(new_aes = "color") +
    geom_point(aes(color = type), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl | Saline rats", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          strip.background = element_rect(
              color = "black", fill = "white", size = 0.3, linetype = "solid", linewidth = 0.4),
          strip.text.x = element_text(size = 10, color = "black")
          ) +
    facet_grid(cols = vars(substance))


## Lever type as focal variable and session as moderator variable
p2 <- ggplot(Short_Acess_Active_LPs_longer,
       aes(x = type, y = number,
           colour = session, group = session)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    labs(y = "Number of presses", x = "Lever type", title = "Fentanyl + Saline rats", color = "Session") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_colour_brewer(palette = "Dark2") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

p2_sub <- ggplot(Short_Acess_Active_LPs_longer,
             aes(x = type, y = number,
                 colour = session, group = session)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    labs(y = "Number of presses", x = "Lever type", title = "Fentanyl | Saline rats", color = "Session") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_colour_brewer(palette = "Dark2") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
    strip.background = element_rect(
              color = "black", fill = "white", size = 0.3, linetype = "solid", linewidth = 0.4),
    strip.text.x = element_text(
              size = 10, color = "black")) +
    facet_grid(cols = vars(substance))

## Session as focal variable and lever type as moderator variable
p3 <- ggplot(Short_Acess_Active_LPs_longer,
       aes(x = session, y = number,
           colour = type, group = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl + Saline rats", color = "Lever type") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

p3_sub <- ggplot(Short_Acess_Active_LPs_longer,
             aes(x = session, y = number,
                 colour = type, group = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl | Saline rats", color = "Lever type") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          strip.background = element_rect(
              color = "black", fill = "white", size = 0.3, linetype = "solid", linewidth = 0.4),
          strip.text.x = element_text(
              size = 10, color = "black")) +
    facet_grid(cols = vars(substance))

plots <- plot_grid(p1, p1_sub, p2, p2_sub, p3, p3_sub, nrow = 3,
                  rel_widths = c(1, 1.5), align = "vh")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_x_sessions_interaction_all_rats.pdf", height = 8, width = 11)


## Test: substance is between-subjects and type and session are within-subjects variables
res.aov <- anova_test(
    data = Short_Acess_Active_LPs_longer, dv = number, wid = rat_sub,
    within = c(type, session), between = substance
)
get_anova_table(res.aov)
# ANOVA Table (type III tests)
#
#                   Effect  DFn   DFd     F     p p<.05   ges
# 1              substance 1.00 17.00 2.842 0.110       0.027
# 2                   type 1.00 17.00 3.922 0.064       0.033
# 3                session 1.03 17.51 2.312 0.146       0.045
# 4         substance:type 1.00 17.00 3.483 0.079       0.030
# 5      substance:session 1.03 17.51 2.049 0.170       0.040
# 6           type:session 1.02 17.41 1.749 0.203       0.034
# 7 substance:type:session 1.02 17.41 2.341 0.144       0.045


# ? ? ? ? ? ? ? ? ? ? ? ?    My questions   ? ? ? ? ? ? ? ? ? ? ? ?
## 1. Meaning of coeffs in linear model - done
summary(lm(number ~ type * session * substance, data = Short_Acess_Active_LPs_longer))

# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)
# (Intercept)                 11.526      2.480   4.647 6.04e-06 ***
# type1                       -6.570      2.480  -2.649  0.00871 **
# session1                    -6.097      5.546  -1.099  0.27291
# session2                    -5.111      5.546  -0.922  0.35783
# session3                    -4.912      5.546  -0.886  0.37680
# session4                    -3.591      5.546  -0.648  0.51800
# session5                     4.426      5.546   0.798  0.42581
# substance1                   5.859      2.480   2.362  0.01910 *
# type1:session1               3.669      5.546   0.662  0.50901
# type1:session2               3.956      5.546   0.713  0.47648
# type1:session3               4.007      5.546   0.723  0.47081
# type1:session4               5.158      5.546   0.930  0.35349
# type1:session5              -3.280      5.546  -0.591  0.55492
# type1:substance1            -6.191      2.480  -2.496  0.01335 *
# session1:substance1         -5.476      5.546  -0.987  0.32465
# session2:substance1         -4.899      5.546  -0.883  0.37809
# session3:substance1         -4.973      5.546  -0.897  0.37095
# session4:substance1         -3.107      5.546  -0.560  0.57601
# session5:substance1          4.001      5.546   0.721  0.47144
# type1:session1:substance1    6.154      5.546   1.110  0.26848
# type1:session2:substance1    4.304      5.546   0.776  0.43858
# type1:session3:substance1    3.628      5.546   0.654  0.51371
# type1:session4:substance1    5.415      5.546   0.976  0.33002
# type1:session5:substance1   -4.022      5.546  -0.725  0.46914
---

## 2. Meaning of interaction terms - done

## 3. F-stat for interaction terms with 2+ predictors, each with 2+ levels

# ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ? ?


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                 Fentanyl rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## Two-way ANOVA paired test
# ------------------------------------------------------------------------------
# Used to investigate if there are differences in the num of lever presses by
# lever type over time (sessions). This, with repeated measures of the same rats
# in the fentanyl group across lever types and sessions (two within-subjects variables).
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality: satisfied for most samples.

fentanyl_rats_lever_presses_longer %>%
    group_by(type, session) %>%
    shapiro_test(number)
#  session type     variable statistic         p
# 1 ShA1    Inactive number       0.602 0.000167
# 2 ShA2    Inactive number       0.808 0.0347
# 3 ShA3    Inactive number       0.877 0.175
# 4 ShA4    Inactive number       0.887 0.220
# 5 ShA5    Inactive number       0.812 0.0386
# 6 ShA6    Inactive number       0.765 0.0120
# 7 ShA1     Active   number       0.896 0.267
# 8 ShA2     Active   number       0.909 0.347
# 9 ShA3     Active   number       0.858 0.115
# 10 ShA4    Active   number       0.911 0.360
# 11 ShA5    Active   number       0.551 0.0000417
# 12 ShA6    Active   number       0.528 0.0000225
# ------------------------------------------------------------------------------

## Differences by lever type across sessions:
p1 <- ggplot(fentanyl_rats_lever_presses_longer,
             aes(x = session, y = number)) +
    theme_light() +
    geom_boxplot(aes(color = type), size = 0.5, outlier.alpha = 0, show.legend = F) +
    scale_colour_manual(values = c("Inactive" = "gray50", "Active" = "gray50")) +
    ggnewscale::new_scale(new_aes = "color") +
    geom_point(aes(color = type), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl rats", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

## Lever type as focal variable and session as moderator variable
p2 <- ggplot(fentanyl_rats_lever_presses_longer,
             aes(x = type, y = number,
                 colour = session, group = session)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    labs(y = "Number of presses", x = "Lever type", title = "Fentanyl rats", color = "Session") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_colour_brewer(palette = "Dark2") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

## Session as focal variable and lever type as moderator variable
p3 <- ggplot(fentanyl_rats_lever_presses_longer,
             aes(x = session, y = number,
                 colour = type, group = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Fentanyl rats", color = "Lever type") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

plots <- plot_grid(p1, p2, p3, align = "vh")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_x_sessions_interaction_fentanyl_rats.pdf", height = 8, width = 11)


## Test
res.aov <- anova_test(
    data = fentanyl_rats_lever_presses_longer, dv = number, wid = rat,
    within = c(type, session)
)
get_anova_table(res.aov)
# ANOVA Table (type III tests)
#
#         Effect  DFn DFd     F     p p<.05   ges
# 1         type 1.00 7.0 2.655 0.147       0.054
# 2      session 1.01 7.1 1.565 0.251       0.072
# 3 type:session 1.01 7.1 1.453 0.267       0.067


# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
#                                 Saline rats
# ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
## Two-way ANOVA paired test
# ------------------------------------------------------------------------------
## Assumptions:
#  A) Normality: satisfied for most samples.

saline_rats_lever_presses_longer %>%
    group_by(type, session) %>%
    shapiro_test(number)
# session type     variable statistic       p
#   <chr>   <fct>    <chr>        <dbl>   <dbl>
# 1  ShA1    Inactive number       0.902 0.194
# 2  ShA2    Inactive number       0.828 0.0223
# 3  ShA3    Inactive number       0.857 0.0519
# 4  ShA4    Inactive number       0.738 0.00147
# 5  ShA5    Inactive number       0.854 0.0481
# 6  ShA6    Inactive number       0.790 0.00695
# 7  ShA1    Active   number       0.879 0.0997
# 8  ShA2    Active   number       0.886 0.122
# 9  ShA3    Active   number       0.918 0.306
# 10 ShA4    Active   number       0.906 0.221
# 11 ShA5    Active   number       0.968 0.865
# 12 ShA6    Active   number       0.764 0.00317
# ------------------------------------------------------------------------------

## Differences by lever type across sessions:
p1 <- ggplot(saline_rats_lever_presses_longer,
             aes(x = session, y = number)) +
    theme_light() +
    geom_boxplot(aes(color = type), size = 0.5, outlier.alpha = 0, show.legend = F) +
    scale_colour_manual(values = c("Inactive" = "gray50", "Active" = "gray50")) +
    ggnewscale::new_scale(new_aes = "color") +
    geom_point(aes(color = type), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Saline rats", color = "Lever type") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

## Lever type as focal variable and session as moderator variable
p2 <- ggplot(saline_rats_lever_presses_longer,
             aes(x = type, y = number,
                 colour = session, group = session)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    labs(y = "Number of presses", x = "Lever type", title = "Saline rats", color = "Session") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    scale_colour_brewer(palette = "Dark2") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

## Session as focal variable and lever type as moderator variable
p3 <- ggplot(saline_rats_lever_presses_longer,
             aes(x = session, y = number,
                 colour = type, group = type)) +
    geom_jitter(width = 0.05, alpha = 0.5) +
    theme_light() +
    scale_color_manual(values = c("Inactive" = "thistle", "Active" = "mediumorchid")) +
    labs(y = "Number of presses", x = "Session", title = "Saline rats", color = "Lever type") +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line") +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11))

plots <- plot_grid(p1, p2, p3, align = "vh")
ggsave(filename = "plots/09_behav_data_statistical_testing/01_Active_vs_Inactive_presses_ShA/lever_type_x_sessions_interaction_saline_rats.pdf", height = 8, width = 11)


## Test
res.aov <- anova_test(
    data = saline_rats_lever_presses_longer, dv = number, wid = rat,
    within = c(type, session)
)
get_anova_table(res.aov)
# ANOVA Table (type III tests)
#
#         Effect DFn DFd     F     p p<.05   ges
# 1         type   1  10 0.513 0.490       0.005
# 2      session   5  50 0.347 0.882       0.009
# 3 type:session   5  50 3.203 0.014     * 0.063






## TODOs: remove outliers and check more carefully assumptions for two-way mixed and paired ANOVA tests




#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
#                            Long Access sessions
#  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------












