# Code to run the 10-fold cross-validation for all the four alternative models (for both datas).
# There are intructions on the comments how to run this for different models. The current setup
# considers the main model.

library(rstan)
library(spdep)
library(tidyverse)
library(loo)

# fixed version of loo::kfold_split_stratified which handles singletons properly
kfold_split_stratified_fixed <- function (K = 10, x = NULL) {
  stopifnot(!is.null(x), K == as.integer(K), length(K) == 1, 
            K > 1, K <= length(x))
  x <- as.integer(as.factor(x))
  Nlev <- length(unique(x))
  N <- length(x)
  xids <- numeric()
  for (l in 1:Nlev) {
    idx <- which(x == l)
    if (length(idx) > 1) {
      xids <- c(xids, sample(idx))
    } else {
      xids <- c(xids, idx)
    }
  }
  bins <- rep(NA, N)
  bins[xids] <- rep(1:K, ceiling(N/K))[1:N]
  return(bins)
}


# The data describe utilised paternal leave quotas and other demographic information on 
# Finnish municipalities in 2009-2017. The political information (main parties, shares of
# votes for each party, familialism and defamilialism) are same for all years. The party
# shares are values between 0-1 as is the variable shareOver18. Other shares are multiplied with
# 100 and thus straightforwardly interpreted as percentages.
#
# The data include variables:
# site = municipality number used by Statistics Finland
# name = name of municipality
# main_party = the main party in the municipality, according to the data omitting SPP
# CD = proprotion of votes for Christian Democrats, according to the data omitting SPP
# Centre = proprotion of votes for Centre Party, according to the data omitting SPP
# NCP = proprotion of votes for National Coalition Party, according to the data omitting SPP
# FP = proprotion of votes for Finns Party, according to the data omitting SPP
# SDP = proprotion of votes for Social Democratic Party of Finland, according to the data omitting SPP
# Left = proprotion of votes for Left Alliance, according to the data omitting SPP
# Green = proprotion of votes for Green League, according to the data omitting SPP
# year = year of the observation
# fin = proportion of Finnish and Sami speakers
# swe = proportion of Swedish speakers
# other = proportion of foreign language speakers
# low_income = proportion of people at risk of poverty rate
# education = education indicator value
# unemp_prop = proportion of unemployed job seekers
# entrepreneurs = proportion of entrepreneurs
# shareOver18 = proportion of fathers using their quota
# n_father = number of fathers
# n_over18 = number of fathers using their quota
# fam = familialisation index, according to the data omitting SPP
# defam = defamilialisation index, according to the data omitting SPP
# main_party_spp = the main party in the municipality, according to the data including SPP
# CD_spp = proprotion of votes for Christian Democrats, according to the data including SPP
# Centre_spp = proprotion of votes for Centre Party, according to the data including SPP
# NCP_spp = proprotion of votes for National Coalition Party, according to the data including SPP
# FP_spp = proprotion of votes for Finns Party, according to the data including SPP
# SDP_spp = proprotion of votes for Social Democratic Party of Finland, according to the data including SPP
# Left_spp = proprotion of votes for Left Alliance, according to the data including SPP
# Green_spp = proprotion of votes for Green League, according to the data including SPP
# SPP_spp = proprotion of votes for Swedish Peoples' Party, according to the data including SPP
# fam_spp = familialisation index, according to the data including SPP
# defam_spp = defamilialisation index, according to the data including SPP
# geom = geometry object for the municipality, retrieved via geofi::get_municipalities

# read data, select variables, and add variables for region and time indeces
# after this the data are ordered primarily by year and secondarily by site
# select appropriate data: variables depending on SPP are denoted with postfix "_spp"
# when they include SPP, other exclude SPP (main model)
# the example code handles the main model, the SPP is included in the comments
leaves <- readRDS("data_leaves.rds") %>% 
  select(site, name, year, swe, other, low_income, education, unemp_prop, entrepreneurs,
         fam, defam, main_party,
 #        fam_spp, defam_spp, main_party_spp, #replace the above line with this to include SPP
         shareOver18, n_father, n_over18, geom) %>% 
  arrange(factor(year)) %>% 
  mutate(region_id = rep(1:length(unique(site)), length(unique(year))),
         time = year - (min(year) - 1))

# number of municipalities
N <- length(unique(leaves$site))

# get the neighbours
tempne <- poly2nb(select(leaves[1:N, ], -geom), queen = TRUE)
mat <- t(sapply(tempne, "length<-", max(lengths(tempne))))

# filter to include only sites with at least one neighbour
#leaves[(mat[, 1] == 0), ]
# this would exclude Hailuoto, Kustavi, Kemiönsaari, and Parainen,
# but we fix this below manually

# add neigbours manually
mat[18, 1] <- 171 # Hailuoto -> Oulu
mat[111, 1] <- 252 # Kustavi -> Taivassalo
mat[116, 1] <- 224 # Kemiönsaari -> Salo
mat[116, 2] <- 225 # Kemiönsaari -> Sauvo
mat[143, 1] <- 64 # Parainen -> Kaarina
# and the other way around
mat[171, ] <- c(sort(c(18, mat[171, ])), rep(NA, 4)) # Oulu -> Hailuoto
mat[252, ] <- c(sort(c(111, mat[252, ])), rep(NA, 9)) # Taivassalo -> Kustavi
mat[224, ] <- c(sort(c(116, mat[224, ])), rep(NA, 4)) # Salo -> Kemiönsaari
mat[225, ] <- c(sort(c(116, mat[225, ])), rep(NA, 8)) # Sauvo -> Kemiönsaari
mat[64, ] <- c(sort(c(143, mat[64, ])), rep(NA, 7)) # Kaarina -> Parainen

# total number of neighbours for each site
n_mat <- rowSums(!is.na(mat))

# neighbourhood structure as pairs for ICAR
edges <- vector("list", N)
for(i in 1:N) {
  n <- mat[i, 1:n_mat[i]]
  n <- n[n>i]
  if(length(n) > 0) {
    edges[[i]] <- cbind(i, n)
  } else {
    edges[[i]] <- NULL
  }
}
edges <- do.call("rbind",edges)

# response = n_over18
response_n <- leaves %>% 
  as.data.frame() %>% 
  pull(matches("n_over"))

# total numbers of fathers
N_all <- leaves %>% 
  as.data.frame() %>% 
  pull(matches("n_father"))

# indeces of non-missing observations
ind <- which(!is.na(response_n))

# number of non-missing observations
n_obs <- sum(!is.na(response_n))

# replace NAs with zeros, these are not used for modelling but to create permissible input data for stan
response_n[which(is.na(response_n))] <- N_all[which(is.na(response_n))] <- 0

# set types as integers for stan
#mode(response_n) <- mode(N_all) <- "integer" # this rounds down
response_n <- as.integer(round(response_n))
N_all <- as.integer(round(N_all))

# indices of regions
region_ind <- leaves %>%
  pull(region_id)

# principal component analysis for family values
pcs <- leaves %>% 
  as.data.frame() %>% 
  select(fam, defam) %>% 
  prcomp(center = TRUE, scale = TRUE)

# time dependent covariates, these are standardised for each year separately
covariates_t <- leaves %>% 
  as.data.frame() %>% 
  select(c(
    "swe",
    "other", 
    "low_income",
    "education",
    "unemp_prop",
    "entrepreneurs",
    "year"
  )) %>%
  group_by(year) %>% 
  mutate_at(vars(-year), scale) %>% 
  ungroup() %>% 
  select(-year) %>% 
  cbind(pcs$x[, 1] %>% scale()) %>% 
  cbind(pcs$x[, 2] %>% scale())

# for the k-fold cross-validation set the number of folds, go through values 1-10
k <- 1 # change this to get other folds

# indices for the K-fold cross-validation, use fixed version to handle the singletons
set.seed(222)
ids <- kfold_split_stratified_fixed(K = 10, x = region_ind[ind])

# collect the parameters for the model, exclude the current fold
pars <- list('T' = length(unique(leaves$year)),
             N = N,
             n_obs = n_obs - sum(ids == k),
             response = response_n,
             N_edges = nrow(edges),
             node1 = edges[, 1],
             node2 = edges[, 2],
             N_all = N_all,
             ind = ind[ids != k],
             region_ind = region_ind[ind[ids != k]],
             n_covariates_t = ncol(covariates_t),
             covariates_t = covariates_t,
             time = (leaves$time)
)

# set the initial values for sampling
# comment in or out appropriate lines for other models
inits <- replicate(4, list(#kappa = 100, # this parameter is needed for the beta-binomial models
  phi = matrix(0, nrow = N, ncol = pars$T), # comment out for models without phi
  sigma_phi = 1, # comment out for models without phi
  c = rep(1, pars$T),
  sigma_c = 1),
  simplify = FALSE)

# construct the stan model, choose appropriate row
#model <- stan_model("model_bin_nophi.stan", allow_optimizations = TRUE)
model <- stan_model("model_bin_phi.stan", allow_optimizations = TRUE)
#model <- stan_model("model_beta_bin_nophi.stan", allow_optimizations = TRUE)
#model <- stan_model("model_beta_bin_phi.stan", allow_optimizations = TRUE)

# fit the model
fit <- sampling(model, data = pars, init = inits,
                iter = 8000, warmup = 1500,
                chains = 4, cores = 4, refresh = 50,
                save_warmup = FALSE)

# gather test parameters to get log likelihood values for the observations excluded from the fit
pars_test <- list('T' = length(unique(leaves$year)),
                  N = N,
                  K = 0,
                  n_obs = sum(ids == k),
                  response = response_n,
                  N_edges = nrow(edges),
                  node1 = edges[, 1],
                  node2 = edges[, 2],
                  N_all = N_all,
                  ind = ind[ids == k],
                  region_ind = region_ind[ind[ids == k]],
                  n_covariates_t = ncol(covariates_t),
                  covariates_t = covariates_t,
                  time = (leaves$time))

# run generated quantities
gen_test <- gqs(model, draws = as.matrix(fit), data = pars_test)

# get the log likelihoods
log_pd <- extract_log_lik(gen_test)

# save the samples, these are needed in the file comparison.R
#saveRDS(log_pd, paste0("log_pd_bin_phi_", k, "_data1.rds")) # full data with SPP
saveRDS(log_pd, paste0("log_pd_bin_phi_", k, "_data2.rds")) # data without SPP
