#!/usr/bin/env Rscriptee
## Rscriptee: Rscript $@ 2>&1 | tee ${OUTPUT_FILE}

################################################################################
### Generate data for different scenarios
##
## Created on: 2018-05-24
## Author: Kazuki Yoshida
################################################################################

## When running non-interactively
.script_name. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(.script_name.) == 1) {
    cat("### Running:", paste(commandArgs()), "\n")
    options(width = 100)
}

## Specify the core count as the second argument
n_cores <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
## Execution not allowed without n_cores
stopifnot(!is.na(n_cores))

cat("
###
### Prepare environment
################################################################################\n")

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")

## Configure parallelization
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)
## Detect core count
## n_cores <- min(parallel::detectCores(), 8)
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
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
library(empeq3)


cat("
###
### Load configurations
################################################################################\n")

## Effect of covariates on log probability ratio of treatments
base_alphas_unmeasured0 <- c(1, 1, 1, 1, 1, 1, 0)
base_alphas_unmeasured0.5 <- c(1, 1, 1, 1, 1, 1, 0.5)
base_alphas_unmeasured1.0 <- c(1, 1, 1, 1, 1, 1, 1.0)
base_alphas_unmeasured2.0 <- c(1, 1, 1, 1, 1, 1, 2.0)

## Multiplier for alphas
sigma0 <- 0.0
sigma0.25 <- 0.25
sigma0.5 <- 0.5
sigma1.0 <- 1.0
sigma2.0 <- 2.0

## Effect of covariates on rate of outcome
XonY1.2 <- rep(log(1.2), 7)
XonY1.2_unmeasured1.5 <- c(rep(log(1.2), 6), log(1.5))
XonY1.2_unmeasured2.0 <- c(rep(log(1.2), 6), log(2.0))


cat("
###
### Generate scenarios
################################################################################\n")

lst_lst_possible_values <-
    list(n = list("n=6000" = 6000),
         p = list("p=7" = 7),
         rho = list("Corr=0" = 0,
                    "Corr=0.1" = 0.1,
                    "Corr=0.3" = 0.3,
                    "Corr=0.5" = 0.5,
                    "Corr=0.7" = 0.7,
                    "Corr=0.9" = 0.9),
         lambda = list("lambda=1" = 1),
         prev = list("prev=0.2" = c(0.2, 0.2, 0.2, 0.2)),
         alphas = list(
             ## 33:33:33
             "size=33:33:33;unmeasured=0;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured0,
                                                          0, 1.0 * sigma0 * base_alphas_unmeasured0),
             "size=33:33:33;unmeasured=0;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured0,
                                                             -0.5, 1.0 * sigma0.25 * base_alphas_unmeasured0),
             "size=33:33:33;unmeasured=0;Equipoise=0.5" = c(-0.34, 0.5 * sigma0.5 * base_alphas_unmeasured0,
                                                            -1.0, 1.0 * sigma0.5 * base_alphas_unmeasured0),
             "size=33:33:33;unmeasured=0;Equipoise=1.0" = c(-0.4, 0.5 * sigma1.0 * base_alphas_unmeasured0,
                                                            -1.7, 1.0 * sigma1.0 * base_alphas_unmeasured0),
             "size=33:33:33;unmeasured=0;Equipoise=2.0" = c(-0.4, 0.5 * sigma2.0 * base_alphas_unmeasured0,
                                                            -3.1, 1.0 * sigma2.0 * base_alphas_unmeasured0),
             ##
             "size=33:33:33;unmeasured=0.5;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured0.5,
                                                            0, 1.0 * sigma0 * base_alphas_unmeasured0.5),
             "size=33:33:33;unmeasured=0.5;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured0.5,
                                                               -0.5, 1.0 * sigma0.25 * base_alphas_unmeasured0.5),
             "size=33:33:33;unmeasured=0.5;Equipoise=0.5" = c(-0.34, 0.5 * sigma0.5 * base_alphas_unmeasured0.5,
                                                              -1.0, 1.0 * sigma0.5 * base_alphas_unmeasured0.5),
             "size=33:33:33;unmeasured=0.5;Equipoise=1.0" = c(-0.4, 0.5 * sigma1.0 * base_alphas_unmeasured0.5,
                                                              -1.7, 1.0 * sigma1.0 * base_alphas_unmeasured0.5),
             "size=33:33:33;unmeasured=0.5;Equipoise=2.0" = c(-0.4, 0.5 * sigma2.0 * base_alphas_unmeasured0.5,
                                                              -3.1, 1.0 * sigma2.0 * base_alphas_unmeasured0.5),
             ##
             "size=33:33:33;unmeasured=1;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured1.0,
                                                          0, 1.0 * sigma0 * base_alphas_unmeasured1.0),
             "size=33:33:33;unmeasured=1;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured1.0,
                                                             -0.5, 1.0 * sigma0.25 * base_alphas_unmeasured1.0),
             "size=33:33:33;unmeasured=1;Equipoise=0.5" = c(-0.34, 0.5 * sigma0.5 * base_alphas_unmeasured1.0,
                                                            -1.0, 1.0 * sigma0.5 * base_alphas_unmeasured1.0),
             "size=33:33:33;unmeasured=1;Equipoise=1.0" = c(-0.4, 0.5 * sigma1.0 * base_alphas_unmeasured1.0,
                                                            -1.7, 1.0 * sigma1.0 * base_alphas_unmeasured1.0),
             "size=33:33:33;unmeasured=1;Equipoise=2.0" = c(-0.4, 0.5 * sigma2.0 * base_alphas_unmeasured1.0,
                                                            -3.1, 1.0 * sigma2.0 * base_alphas_unmeasured1.0),
             ##
             "size=33:33:33;unmeasured=2;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured2.0,
                                                          0, 1.0 * sigma0 * base_alphas_unmeasured2.0),
             "size=33:33:33;unmeasured=2;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured2.0,
                                                             -0.5, 1.0 * sigma0.25 * base_alphas_unmeasured2.0),
             "size=33:33:33;unmeasured=2;Equipoise=0.5" = c(-0.34, 0.5 * sigma0.5 * base_alphas_unmeasured2.0,
                                                            -1.0, 1.0 * sigma0.5 * base_alphas_unmeasured2.0),
             "size=33:33:33;unmeasured=2;Equipoise=1.0" = c(-0.4, 0.5 * sigma1.0 * base_alphas_unmeasured2.0,
                                                            -1.7, 1.0 * sigma1.0 * base_alphas_unmeasured2.0),
             "size=33:33:33;unmeasured=2;Equipoise=2.0" = c(-0.4, 0.5 * sigma2.0 * base_alphas_unmeasured2.0,
                                                            -3.1, 1.0 * sigma2.0 * base_alphas_unmeasured2.0),
             ##
             ## 10:45:45
             "size=10:45:45;unmeasured=0;Equipoise=0" = c(1.5, 0.5 * sigma0 * base_alphas_unmeasured0,
                                                          1.5, 1.0 * sigma0 * base_alphas_unmeasured0),
             "size=10:45:45;unmeasured=0;Equipoise=0.25" = c(1.2, 0.5 * sigma0.25 * base_alphas_unmeasured0,
                                                             1.0, 1.0 * sigma0.25 * base_alphas_unmeasured0),
             "size=10:45:45;unmeasured=0;Equipoise=0.5" = c(0.9, 0.5 * sigma0.5 * base_alphas_unmeasured0,
                                                            0.4, 1.0 * sigma0.5 * base_alphas_unmeasured0),
             "size=10:45:45;unmeasured=0;Equipoise=1.0" = c(1.5, 0.5 * sigma1.0 * base_alphas_unmeasured0,
                                                            0.5, 1.0 * sigma1.0 * base_alphas_unmeasured0),
             "size=10:45:45;unmeasured=0;Equipoise=2.0" = c(1.7, 0.5 * sigma2.0 * base_alphas_unmeasured0,
                                                            -0.3, 1.0 * sigma2.0 * base_alphas_unmeasured0),
             ##
             "size=10:45:45;unmeasured=0.5;Equipoise=0" = c(1.5, 0.5 * sigma0 * base_alphas_unmeasured0.5,
                                                            1.5, 1.0 * sigma0 * base_alphas_unmeasured0.5),
             "size=10:45:45;unmeasured=0.5;Equipoise=0.25" = c(1.2, 0.5 * sigma0.25 * base_alphas_unmeasured0.5,
                                                               1.0, 1.0 * sigma0.25 * base_alphas_unmeasured0.5),
             "size=10:45:45;unmeasured=0.5;Equipoise=0.5" = c(0.9, 0.5 * sigma0.5 * base_alphas_unmeasured0.5,
                                                              0.4, 1.0 * sigma0.5 * base_alphas_unmeasured0.5),
             "size=10:45:45;unmeasured=0.5;Equipoise=1.0" = c(1.5, 0.5 * sigma1.0 * base_alphas_unmeasured0.5,
                                                              0.5, 1.0 * sigma1.0 * base_alphas_unmeasured0.5),
             "size=10:45:45;unmeasured=0.5;Equipoise=2.0" = c(1.7, 0.5 * sigma2.0 * base_alphas_unmeasured0.5,
                                                              -0.3, 1.0 * sigma2.0 * base_alphas_unmeasured0.5),
             ##
             "size=10:45:45;unmeasured=1;Equipoise=0" = c(1.5, 0.5 * sigma0 * base_alphas_unmeasured1.0,
                                                          1.5, 1.0 * sigma0 * base_alphas_unmeasured1.0),
             "size=10:45:45;unmeasured=1;Equipoise=0.25" = c(1.2, 0.5 * sigma0.25 * base_alphas_unmeasured1.0,
                                                             1.0, 1.0 * sigma0.25 * base_alphas_unmeasured1.0),
             "size=10:45:45;unmeasured=1;Equipoise=0.5" = c(0.9, 0.5 * sigma0.5 * base_alphas_unmeasured1.0,
                                                            0.4, 1.0 * sigma0.5 * base_alphas_unmeasured1.0),
             "size=10:45:45;unmeasured=1;Equipoise=1.0" = c(1.5, 0.5 * sigma1.0 * base_alphas_unmeasured1.0,
                                                            0.5, 1.0 * sigma1.0 * base_alphas_unmeasured1.0),
             "size=10:45:45;unmeasured=1;Equipoise=2.0" = c(1.7, 0.5 * sigma2.0 * base_alphas_unmeasured1.0,
                                                            -0.3, 1.0 * sigma2.0 * base_alphas_unmeasured1.0),
             ##
             "size=10:45:45;unmeasured=2;Equipoise=0" = c(1.5, 0.5 * sigma0 * base_alphas_unmeasured2.0,
                                                          1.5, 1.0 * sigma0 * base_alphas_unmeasured2.0),
             "size=10:45:45;unmeasured=2;Equipoise=0.25" = c(1.2, 0.5 * sigma0.25 * base_alphas_unmeasured2.0,
                                                             1.0, 1.0 * sigma0.25 * base_alphas_unmeasured2.0),
             "size=10:45:45;unmeasured=2;Equipoise=0.5" = c(0.9, 0.5 * sigma0.5 * base_alphas_unmeasured2.0,
                                                            0.4, 1.0 * sigma0.5 * base_alphas_unmeasured2.0),
             "size=10:45:45;unmeasured=2;Equipoise=1.0" = c(1.5, 0.5 * sigma1.0 * base_alphas_unmeasured2.0,
                                                            0.5, 1.0 * sigma1.0 * base_alphas_unmeasured2.0),
             "size=10:45:45;unmeasured=2;Equipoise=2.0" = c(1.7, 0.5 * sigma2.0 * base_alphas_unmeasured2.0,
                                                            -0.3, 1.0 * sigma2.0 * base_alphas_unmeasured2.0),
             ##
             ## 10:10:80
             "size=10:10:80;unmeasured=0;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured0,
                                                          2, 1.0 * sigma0 * base_alphas_unmeasured0),
             "size=10:10:80;unmeasured=0;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured0,
                                                             1.7, 1.0 * sigma0.25 * base_alphas_unmeasured0),
             "size=10:10:80;unmeasured=0;Equipoise=0.5" = c(0, 0.5 * sigma0.5 * base_alphas_unmeasured0,
                                                            1.7, 1.0 * sigma0.5 * base_alphas_unmeasured0),
             "size=10:10:80;unmeasured=0;Equipoise=1.0" = c(0.1, 0.5 * sigma1.0 * base_alphas_unmeasured0,
                                                            1.6, 1.0 * sigma1.0 * base_alphas_unmeasured0),
             "size=10:10:80;unmeasured=0;Equipoise=2.0" = c(0.8, 0.5 * sigma2.0 * base_alphas_unmeasured0,
                                                            2.0, 1.0 * sigma2.0 * base_alphas_unmeasured0),
             ##
             "size=10:10:80;unmeasured=0.5;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured0.5,
                                                            2, 1.0 * sigma0 * base_alphas_unmeasured0.5),
             "size=10:10:80;unmeasured=0.5;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured0.5,
                                                               1.7, 1.0 * sigma0.25 * base_alphas_unmeasured0.5),
             "size=10:10:80;unmeasured=0.5;Equipoise=0.5" = c(0, 0.5 * sigma0.5 * base_alphas_unmeasured0.5,
                                                              1.7, 1.0 * sigma0.5 * base_alphas_unmeasured0.5),
             "size=10:10:80;unmeasured=0.5;Equipoise=1.0" = c(0.1, 0.5 * sigma1.0 * base_alphas_unmeasured0.5,
                                                              1.6, 1.0 * sigma1.0 * base_alphas_unmeasured0.5),
             "size=10:10:80;unmeasured=0.5;Equipoise=2.0" = c(0.8, 0.5 * sigma2.0 * base_alphas_unmeasured0.5,
                                                              2.0, 1.0 * sigma2.0 * base_alphas_unmeasured0.5),
             ##
             "size=10:10:80;unmeasured=1;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured1.0,
                                                          2, 1.0 * sigma0 * base_alphas_unmeasured1.0),
             "size=10:10:80;unmeasured=1;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured1.0,
                                                             1.7, 1.0 * sigma0.25 * base_alphas_unmeasured1.0),
             "size=10:10:80;unmeasured=1;Equipoise=0.5" = c(0, 0.5 * sigma0.5 * base_alphas_unmeasured1.0,
                                                            1.7, 1.0 * sigma0.5 * base_alphas_unmeasured1.0),
             "size=10:10:80;unmeasured=1;Equipoise=1.0" = c(0.1, 0.5 * sigma1.0 * base_alphas_unmeasured1.0,
                                                            1.6, 1.0 * sigma1.0 * base_alphas_unmeasured1.0),
             "size=10:10:80;unmeasured=1;Equipoise=2.0" = c(0.8, 0.5 * sigma2.0 * base_alphas_unmeasured1.0,
                                                            2.0, 1.0 * sigma2.0 * base_alphas_unmeasured1.0),
             ##
             "size=10:10:80;unmeasured=2;Equipoise=0" = c(0, 0.5 * sigma0 * base_alphas_unmeasured2.0,
                                                          2, 1.0 * sigma0 * base_alphas_unmeasured2.0),
             "size=10:10:80;unmeasured=2;Equipoise=0.25" = c(-0.2, 0.5 * sigma0.25 * base_alphas_unmeasured2.0,
                                                             1.7, 1.0 * sigma0.25 * base_alphas_unmeasured2.0),
             "size=10:10:80;unmeasured=2;Equipoise=0.5" = c(0, 0.5 * sigma0.5 * base_alphas_unmeasured2.0,
                                                            1.7, 1.0 * sigma0.5 * base_alphas_unmeasured2.0),
             "size=10:10:80;unmeasured=2;Equipoise=1.0" = c(0.1, 0.5 * sigma1.0 * base_alphas_unmeasured2.0,
                                                            1.6, 1.0 * sigma1.0 * base_alphas_unmeasured2.0),
             "size=10:10:80;unmeasured=2;Equipoise=2.0" = c(0.8, 0.5 * sigma2.0 * base_alphas_unmeasured2.0,
                                                            2.0, 1.0 * sigma2.0 * base_alphas_unmeasured2.0)
         ),
         beta0 = list("Outcome=0.20" = log(0.20)),
         betaA = list("Main=0,0" = c(0,0)),
         betaX = list("Xu_on_Y=1.2" =  XonY1.2,
                      "Xu_on_Y=1.5" =  XonY1.2_unmeasured1.5,
                      "Xu_on_Y=2.0" =  XonY1.2_unmeasured2.0),
         betaXA = list("Modification=0" = c(0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0)))

## Split alphas, prev_params, and contraindication
scenarios <- datagen3::generate_scenario_data_frame(lst_lst_possible_values)
scenarios
scenarios$description


## Move some new scenarios to the end.
## unmeasured0 scenarios and missing scenario.
scenarios <- scenarios %>%
    mutate(unmeasured0 = grepl("unmeasured=0;", description),
           missing_scenario = grepl("size=10:45:45;unmeasured=0.5;Equipoise=0.5", description)) %>%
    ## FALSE comes before TRUE
    arrange(unmeasured0, missing_scenario)

scenarios %>%
    group_by(unmeasured0, missing_scenario) %>%
    summarize(n = n())

## Drop arrange variables
scenarios <- scenarios %>%
    select(-unmeasured0, -missing_scenario)

## Add back the scenario class
class(scenarios) <- c("scenarios", class(scenarios))


cat("
###
### Generate datasets
################################################################################\n")

n_parts <- 1
R <- 200

datagen3::generate_data_for_all_scenarios(fun = generate_p_norm_count_bin_data_count,
                                          scenarios = scenarios,
                                          n_parts = n_parts,
                                          R = R,
                                          ## Skip these rows
                                          skip = 1:792,
                                          path = "./data/")

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
