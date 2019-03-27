#!/usr/bin/env ./Rscriptee
## Use the project-specific ./Rscriptee

################################################################################
### Assess aggregated results
##
## Created on: 2018-05-25
## Author: Kazuki Yoshida
################################################################################

## When running non-interactively
.script_name. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(.script_name.) == 1) {
    cat("### Running:", paste(commandArgs()), "\n")
    options(width = 100)
}

###
### Capture data filename argument
################################################################################

## Specify data file as the first argument
data_file_name <- commandArgs(trailingOnly = TRUE)[1]
## Specify the core count as the second argument
n_cores <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

## Execution not allowed without data file
stopifnot(!is.na(data_file_name))
## Execution not allowed without n_cores
stopifnot(!is.na(n_cores))

## Check it is a scenario file
if (!grepl("all_analysis_summary", data_file_name)) {
    stop("Not a summary result file")
}


###
### Prepare environment
################################################################################

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")

## Configure parallelization
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)
## Used by parallel::mclapply() as default
options(mc.cores = n_cores)
## Used by doParallel as default
options(cores = n_cores)
## Register doParallel as the parallel backend for foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = n_cores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Load packages
library(tidyverse)
library(grid)
library(gtable)



cat("
###
### Define functions
################################################################################\n")

## Grid label manipulation
## https://stackoverflow.com/questions/40732543/seeking-workaround-for-gtable-add-grob-code-broken-by-ggplot-2-2-0/40838827#40838827
## Need to overwrite the strip. Thus, make strip background white.
clean_contrast_trim_method_columns <- function(gg) {
    ## Generate a ggplot2 plot grob.
    pg <- ggplotGrob(gg)
    ## Get a list of strips from the original plot
    strip <- lapply(grep("strip-t", pg$layout$name), function(x) {pg$grobs[[x]]})
    ## Construct gtable to contain the new strip
    newStrip <- gtable(widths = unit(rep(1, length(strip)), "null"), heights = strip[[1]]$heights)
    ## Top row
    cols <- seq(1, by = 3, length.out = length(strip)/3)
    newStrip <- gtable_add_grob(newStrip,
                                lapply(strip[cols], `[`, 1),
                                t = 1, l = cols, r = cols + 2)
    ## Bottom row
    newStrip <- gtable_add_grob(newStrip,
                                lapply(strip, `[`, 2),
                                t = 2, l = seq_along(strip))
    ## Put the strip into the plot
    pgNew <- gtable_add_grob(pg, newStrip, t = 6, l = 4, r = 20)
    ## Vertical lines
    pgNew <- gtable_add_grob(x = pgNew,
                             grobs = segmentsGrob(x0 = unit(0.5, "npc"), x1 = unit(0.5, "npc"),
                                                  y0 = unit(0, "npc"), y1 = unit(1, "npc"),
                                                  gp = gpar(lty = 2)),
                             t = 6, l = 9, b = 17, r = 9)
    pgNew <- gtable_add_grob(x = pgNew,
                             grobs = segmentsGrob(x0 = unit(0.5, "npc"), x1 = unit(0.5, "npc"),
                                                  y0 = unit(0, "npc"), y1 = unit(1, "npc"),
                                                  gp = gpar(lty = 2)),
                             t = 6, l = 15, b = 17, r = 15)
    pgNew
}


cat("
###
### Load summary
################################################################################\n")

cat("###  Loading", data_file_name, "\n")
load(data_file_name)

pryr::object_size(df_scenarios)
df_scenarios %>%
    rename(`#` = scenario_count) %>%
    print(n = Inf, width = Inf)

pryr::object_size(df_summary)
df_summary %>%
    print(n = 20)


cat("
###
### Split description for easier use
################################################################################\n")

df_summary <- df_summary %>%
    ## Only retain untrimmed and trimmed at Walker with threshold 0.2
    filter(thres %in% c(0, 0.2)) %>%
    separate(col = scenario_description,
             into = c("n","p","correlation","lambda","prev","size","unmeasured","equipoise","incidence","main_effect","x_on_y","modification"),
             sep = ";") %>%
    mutate(n = as.numeric(gsub(".*=","", n)),
           p = as.numeric(gsub(".*=","", p)),
           correlation = as.numeric(gsub(".*=","", correlation)),
           lambda = as.numeric(gsub(".*=","", lambda)),
           prev = as.numeric(gsub(".*=","", prev)),
           size = gsub(".*=","", size),
           unmeasured = as.numeric(gsub(".*=","", unmeasured)),
           equipoise = as.numeric(gsub(".*=","", equipoise)),
           incidence = as.numeric(gsub(".*=","", incidence)),
           main_effect = gsub(".*=","", main_effect),
           x_on_y = as.numeric(gsub(".*=","", x_on_y)),
           modification = as.numeric(gsub(".*=","", modification))) %>%
    ## Ordering only
    mutate(adjustment = factor(adjustment,
                               levels = c("unadj","iptw","mw","ow"),
                               labels = c("Unadj","IPTW","MW","OW"))) %>%
    mutate(p_kept = mean_n_kept / n)


cat("
###
### Sanity check results for the correctness of true results
################################################################################\n")

cat("###  Examine scenario 1 results \n")
df_scenarios %>%
    filter(scenario_count  == 1)

cat("###   No trimming \n")
df_summary %>%
    filter(scenario_count  == 1,
           thres == 0,
           measure == "coef") %>%
    mutate(mean = round(mean, 3),
           var = round(var, 3)) %>%
    print(n = 100)


cat("
###
### Sample size calculations
################################################################################\n")

cat("
###  Full cohort\n")
df_n <- df_summary %>%
    ## remove redundancies
    filter(measure == "coef",
           contrast == "1vs0",
           adjustment == "Unadj",
           !true,
           !reest)
df_n %>%
    select(scenario_count, trim_method_name, thres,
           mean_n_kept,
           p_kept)


cat("
###  Group-wise\n")
df_n_group <- df_n %>%
    ## Change to long format by the treatment groups
    ## Three times longer than df_n
    gather(key = group, value = mean_n_group_kept,
           mean_n0_kept, mean_n1_kept, mean_n2_kept) %>%
    ## Extract group names
    mutate(group = group %>%
               gsub(pattern = "mean_n", replacement = "") %>%
               gsub(pattern = "_kept", replacement = ""))
## Show
df_n_group %>%
    arrange(scenario_count, thres, group) %>%
    select(scenario_count, trim_method_name, thres,
           mean_n_kept,
           group,
           mean_n_group_kept) %>%
    print(n = 200)

## Calculate group-wise proportions with respect group-wise full size
df_p_mean_n_group_kept <- df_n_group %>%
    select(scenario_count, thres, group, mean_n_group_kept) %>%
    filter(thres %in% c(0, 0.2)) %>%
    spread(key = thres, value = mean_n_group_kept) %>%
    ## Divide group-wise sample size when trimmed by group-wise sample size when untrimmed.
    mutate(p_mean_n_group_kept = `0.2` / `0`) %>%
    select(-`0.2`, -`0`) %>%
    spread(key = group, value = p_mean_n_group_kept) %>%
    rename(p_mean_n0_kept = `0`,
           p_mean_n1_kept = `1`,
           p_mean_n2_kept = `2`) %>%
    mutate(p_mean_n_group_kept_min = pmin(p_mean_n0_kept, p_mean_n1_kept, p_mean_n2_kept))
df_p_mean_n_group_kept

cat("###  Any duplicated scenarios?\n")
any(duplicated(names(table(df_p_mean_n_group_kept$scenario_count))))


cat("
###
### Bias calculation for mean coef (compare to corresponding true value)
################################################################################\n")

df_bias <- df_summary %>%
    filter(measure == "coef") %>%
    select(-var, -sd) %>%
    spread(key = true, value = mean) %>%
    ## mean of estimates - mean of calculated truth
    mutate(bias = `FALSE` - `TRUE`) %>%
    ## Drop these variables
    select(-`FALSE`, -`TRUE`)


cat("
###
### Assign proportion kept at different threshold to the full cohort result
################################################################################\n")
## Bias is obtained from the full cohort analysis because we are assessing the appropriateness of the full-cohort analysis.

## Pull p_kept from threshold 0.20 and add to full cohort bias results
##                                # Full cohort analysis except for p_kept in 0.20 region
df_bias_full_cohort <- inner_join(df_bias %>%
                                  filter(thres == 0) %>%
                                  select(-p_kept),
                                  ##
                                  ## Trimmed cohort analysis for p_kept in 0.20 region
                                  df_bias %>%
                                  filter(thres == 0.20,
                                         reest == FALSE) %>%
                                  select(scenario_count,
                                         adjustment, measure, contrast,
                                         p_kept))
cat("###  df_bias_full_cohort\n")
df_bias_full_cohort

cat("###  df_bias_full_cohort with p_mean_n_group_kept_min\n")
## Pull p_mean_n_group_kept_min from df_p_mean_n_group_kept
df_bias_full_cohort <- left_join(df_bias_full_cohort,
                                 df_p_mean_n_group_kept)
df_bias_full_cohort

df_coef_full_cohort <- inner_join(df_summary %>%
                                  filter(measure == "coef",
                                         thres == 0) %>%
                                  select(-p_kept),
                                  ##
                                  ## Trimmed cohort analysis for p_kept in 0.20 region
                                  df_summary %>%
                                  filter(measure == "coef",
                                         thres == 0.20,
                                         reest == FALSE) %>%
                                  select(scenario_count,
                                         adjustment, measure, contrast,
                                         p_kept))


cat("
###
### Preliminary visualizations
################################################################################\n")

df_bias_full_cohort <- df_bias_full_cohort %>%
    arrange(equipoise)

if (FALSE) {

cat("
###  Coefficients vs p_kept (overall)\n")
pdf(file = "./out/exp_mean_coef_vs_p_kept.pdf", width = 10, height = 7, family = "sans")

## Plot coefficients
df_coef_full_cohort %>%
    ## Need to restrict to full cohort analyses only
    filter(trim_method_name == "none",
           thres == 0,
           ##
           measure == "coef",
           !reest,
           !true) %>%
    group_by(size, x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(size, x_on_y, data), function(size, x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(mean), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_line() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ adjustment + contrast, scales = "free_x") +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("exp(mean(coef)); Group sizes:%s; Xu on Y: %s", size, x_on_y),
                 x = "exp(mean(coef))", y = "Proportion in Walker region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
        magrittr::extract2("gg")

dev.off()


cat("
###  Bias vs p_kept (overall)\n")
pdf(file = "./out/bias_vs_p_kept.pdf", width = 10, height = 7, family = "sans")

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(unmeasured ~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(correlation ~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    geom_smooth() +
    facet_grid(correlation ~ adjustment) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept, color = factor(unmeasured), shape = factor(equipoise))) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept, color = factor(unmeasured), shape = factor(equipoise))) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "33:33:33") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "33:33:33", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "10:45:45") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "10:45:45", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "10:10:80") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "10:10:80", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

## Automate panels
df_bias_full_cohort %>%
    group_by(size, x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(size, x_on_y, data), function(size, x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ adjustment + contrast, scales = "free_x") +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("Group sizes:%s; Xu on Y: %s", size, x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(size, x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(size, x_on_y, data), function(size, x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ adjustment + contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("Group sizes:%s; Xu on Y: %s", size, x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            filter(contrast == "1vs0") %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment, scales = "free_x") +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            filter(contrast == "1vs0") %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(contrast == "1vs0",
           size == "33:33:33") %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 and 33:33:33 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(contrast == "1vs0",
           size == "33:33:33",
           equipoise != 0.25,
           correlation == 0) %>%
    ggplot(mapping = aes(x = exp(bias), y = p_kept,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(x_on_y ~ size + adjustment) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = sprintf("1vs0 and 33:33:33 only"),
         x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

dev.off()


cat("
###  Bias vs p_mean_n_group_kept_min (group-wise minimum)\n")
pdf(file = "./out/bias_vs_p_mean_n_group_kept_min.pdf", width = 10, height = 7, family = "sans")

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(unmeasured ~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    facet_grid(correlation ~ adjustment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    geom_smooth() +
    facet_grid(correlation ~ adjustment) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min, color = factor(unmeasured), shape = factor(equipoise))) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min, color = factor(unmeasured), shape = factor(equipoise))) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    ggplot(mapping = aes(x = bias, y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast, labeller = label_both) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Bias in log rate ratio", y = "Proportion in Walker region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "33:33:33") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "33:33:33", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "10:45:45") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "10:45:45", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

df_bias_full_cohort %>%
    filter(size == "10:10:80") %>%
    ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(correlation ~ adjustment + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "10:10:80", x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

## Automate panels
df_bias_full_cohort %>%
    group_by(size, x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(size, x_on_y, data), function(size, x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ adjustment + contrast, scales = "free_x") +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("Group sizes:%s; Xu on Y: %s", size, x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(size, x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(size, x_on_y, data), function(size, x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ adjustment + contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("Group sizes:%s; Xu on Y: %s", size, x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            filter(contrast == "1vs0") %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment, scales = "free_x") +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            filter(contrast == "1vs0") %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(contrast == "1vs0",
           size == "33:33:33") %>%
    group_by(x_on_y) %>%
    nest() %>%
    mutate(gg = pmap(list(x_on_y, data), function(x_on_y, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 color = factor(unmeasured), group = factor(unmeasured),
                                 shape = factor(equipoise))) +
            geom_point() +
            ## This geom gets the ordering correct (order in the data set?).
            geom_path() +
            geom_hline(yintercept = 0.5) +
            ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
            facet_grid(correlation ~ size + adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("1vs0 and 33:33:33 only; Xu on Y: %s", x_on_y),
                 x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(contrast == "1vs0",
           size == "33:33:33",
           equipoise != 0.25,
           correlation == 0) %>%
    ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                         color = factor(unmeasured), group = factor(unmeasured),
                         shape = factor(equipoise))) +
    geom_point() +
    ## This geom gets the ordering correct (order in the data set?).
    geom_path() +
    geom_hline(yintercept = 0.5) +
    ## geom_smooth(se = FALSE, alpha = 0.5, size = 0.5) +
    facet_grid(x_on_y ~ size + adjustment) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = sprintf("1vs0 and 33:33:33 only"),
         x = "Multiplicative Bias in Rate Ratio", y = "Proportion in Walker Region") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())

dev.off()


}
cat("
###
### Production figures
################################################################################\n")

equipoise_levels <- c("Perfect", "Good", "Moderate", "Poor")
unmeasured_levels <- c("Zero", "Half", "Same", "Twice")


pdf(file = "./out/bias_vs_p_kept_no_corr.pdf", width = 10, height = 7, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0) %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size, contrast) %>%
    nest() %>%
    mutate(gg = pmap(list(size, contrast, data), function(size, contrast, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s; %s", size, contrast),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./out/bias_vs_p_kept_no_corr_slides.pdf", width = 10, height = 7, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./out/bias_vs_p_kept_no_corr_slides_focus.pdf", width = 5, height = 3.5, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5) %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./out/bias_vs_p_kept_no_corr_slides_focus_build.pdf", width = 5, height = 3.5, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point(alpha = 0) +
            geom_path(alpha = 0) +
            geom_hline(yintercept = 0.5, alpha = 0) +
            geom_segment(x = 1, y = 0, xend = 2.1, yend = 0, arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            annotate(geom = "text", x = 1.5, y = 0.05, label = "Residual Bias", size = 5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point(alpha = 0) +
            geom_path(alpha = 0) +
            geom_hline(yintercept = 0.5, alpha = 0) +
            geom_segment(x = 1, y = 0, xend = 1, yend = 1, arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            annotate(geom = "text", x = 1, y = 0.8, label = "Proportion in Region\nof Empirical Equipoise", size = 5, hjust = 0) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point(alpha = 0) +
            geom_path(alpha = 0) +
            geom_segment(x = 1, y = 0, xend = 2.1, yend = 0, arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            annotate(geom = "text", x = 1.5, y = 0.05, label = "Residual Bias", size = 5) +
            geom_segment(x = 1, y = 0, xend = 1, yend = 1, arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            annotate(geom = "text", x = 1, y = 0.8, label = "Proportion in Region\nof Empirical Equipoise", size = 5, hjust = 0) +
            geom_hline(yintercept = 0.5) +
            annotate(geom = "text", x = 1.5, y = 0.5, label = "50% threshold", size = 5, vjust = 0) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            geom_segment(x = 1.8, y = 0, xend = 1.8, yend = 0.15,
                         arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            geom_point(x = 1.37, y = 0.5, shape = 4, size = 5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            geom_segment(x = 1.2, y = 0, xend = 1.01, yend = 0.15,
                         arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            geom_point(x = 1, y = 0.5, shape = 4, size = 5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0,
           adjustment == "MW",
           contrast == "2vs0",
           x_on_y == 1.5,
           size == "33:33:33") %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size) %>%
    nest() %>%
    mutate(gg = pmap(list(size, data), function(size, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            geom_segment(x = 1.9, y = 0.25, xend = 2.0, yend = 0.4,
                         arrow = arrow(length = unit(0.5, "cm")), show.legend = FALSE) +
            geom_point(x = 1.825, y = 0.5, shape = 4, size = 5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ contrast) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s", size),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./out/bias_vs_p_mean_n_group_kept_min_no_corr.pdf", width = 10, height = 7, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           correlation == 0) %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size, contrast) %>%
    nest() %>%
    mutate(gg = pmap(list(size, contrast, data), function(size, contrast, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_mean_n_group_kept_min,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(x_on_y ~ adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s; %s", size, contrast),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Minimum group-wise proportion in empirical equipoise region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./out/bias_vs_p_kept_corr.pdf", width = 10, height = 7, family = "sans")

## Automate panels
df_bias_full_cohort %>%
    filter(equipoise != 0.25,
           x_on_y == 1.5) %>%
    mutate(equipoise = factor(equipoise, labels = equipoise_levels),
           unmeasured = factor(unmeasured, labels = unmeasured_levels)) %>%
    group_by(size, contrast) %>%
    nest() %>%
    mutate(gg = pmap(list(size, contrast, data), function(size, contrast, data) {
        data %>%
            ggplot(mapping = aes(x = exp(bias), y = p_kept,
                                 linetype = unmeasured, group = unmeasured)) +
            geom_point() +
            geom_path() +
            geom_hline(yintercept = 0.5) +
            scale_shape() +
            scale_linetype() +
            facet_grid(correlation ~ adjustment) +
            coord_cartesian(ylim = c(0, 1)) +
            labs(title = sprintf("%s; %s", size, contrast),
                 x = "Multiplicative Bias in Rate Ratio",
                 y = "Proportion in Empirical Equipoise Region",
                 linetype = "Relative Treatment\nAssociation of X7") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  strip.background = element_blank())
    })) %>%
    magrittr::extract2("gg")

dev.off()


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("\n### Started  ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
