#main
library(cmdstanr)
library(parallel)

#step 1: load data--------
source('data_load.R')

#Step 2: run stan model--------
source('stan_code.R')
n.sim=5e3
model=cmdstan_model(stan_file = write_stan_file(stan_code))
set.seed(123); n.chains=5
fit=model$sample(data = data_list, chains = n.chains, iter_sampling = n.sim, iter_warmup = 2e3, parallel_chains = n.chains)

#Step 3: derive estimates-------
source('analysis.R')