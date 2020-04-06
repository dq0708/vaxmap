rm(list=ls())

library(data.table)
library(survey)
library(sp)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(raster)
library(maptools)
library(latticeExtra)
library(xtable)
library(Hmisc)
library(INLA)
library(spdep)

options(survey.adjust.domain.lonely = T)
options(survey.lonely.psu = "adjust")

##################################
# load previous results
##################################

load("data_cleaned_1k.RData")

inladata_dt <- cluster_data_dt[, c("i_id", "j_id", "c_id", "IsUrban", "long", "lat", "y", "n", 
                                   "aridity10k", "poverty", "evi1k", "log_travel_time", "log_night"), with = F]
inladata_dt[, "beta0"] <- 1

range0 <- min(c(diff(range(nga_mesh$loc[, 1])), diff(range(nga_mesh$loc[, 2]))))/5
nga_spde <- inla.spde2.pcmatern(nga_mesh, prior.range = c(range0, 0.5), prior.sigma = c(1, 0.01))

A_obs <- inla.spde.make.A(mesh = nga_mesh,
                          loc = cbind(inladata_dt$long, inladata_dt$lat))

stack_obs <- inla.stack(tag = "obs", 
                        data = list(y = inladata_dt$y,
                                    n = inladata_dt$n),
                        A = list(A_obs, 1),
                        effects = list(u = 1:nga_mesh$n,
                                       data.frame(beta0 = inladata_dt$beta0, 
                                                  IsUrban = inladata_dt$IsUrban,
                                                  aridity10k = inladata_dt$aridity10k,
                                                  poverty = inladata_dt$poverty,
                                                  evi1k = inladata_dt$evi1k,
                                                  log_travel_time = inladata_dt$log_travel_time,
                                                  log_night = inladata_dt$log_night,
                                                  i_id = inladata_dt$i_id,
                                                  j_id = inladata_dt$j_id,
                                                  c_id = inladata_dt$c_id)))

##################################
# cov_nostrat_Binom_SPDE_NE, i.e. Utazi et al 2018
##################################

cov_nostrat_Binom_SPDE_NE <- y ~ -1 + beta0 +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde)

cov_nostrat_Binom_SPDE_NE_fit <- inla(formula = cov_nostrat_Binom_SPDE_NE,
                                      family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                      data = inla.stack.data(stack_obs),
                                      control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                      control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                      verbose = F)

save(cov_nostrat_Binom_SPDE_NE_fit, nga_spde, A_obs, stack_obs,
     file = "cov_nostrat_Binom_SPDE_NE_fit.RData")

##################################
# cov_strat_Binom_SPDE_NE
##################################

cov_strat_Binom_SPDE_NE <- y ~ -1 + beta0 +
  IsUrban +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde)

cov_strat_Binom_SPDE_NE_fit <- inla(formula = cov_strat_Binom_SPDE_NE,
                                    family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                    data = inla.stack.data(stack_obs),
                                    control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                    control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                    verbose = F)

save(cov_strat_Binom_SPDE_NE_fit, nga_spde, A_obs, stack_obs,
     file = "cov_strat_Binom_SPDE_NE_fit.RData")

##################################
# cov_nostrat_LonoB_SPDE_OD
##################################

cov_nostrat_LonoB_SPDE_OD <- y ~ -1 + beta0 +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde) +
  f(c_id, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

cov_nostrat_LonoB_SPDE_OD_fit <- inla(formula = cov_nostrat_LonoB_SPDE_OD,
                                      family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                      data = inla.stack.data(stack_obs),
                                      control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                      control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                      verbose = F)

save(cov_nostrat_LonoB_SPDE_OD_fit, nga_spde, A_obs, stack_obs,
     file = "cov_nostrat_LonoB_SPDE_OD_fit.RData")

##################################
# cov_strat_LonoB_SPDE_OD
##################################

cov_strat_LonoB_SPDE_OD <- y ~ -1 + beta0 +
  IsUrban +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde) +
  f(c_id, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

cov_strat_LonoB_SPDE_OD_fit <- inla(formula = cov_strat_LonoB_SPDE_OD,
                                    family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                    data = inla.stack.data(stack_obs),
                                    control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                    control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                    verbose = F)

save(cov_strat_LonoB_SPDE_OD_fit, nga_spde, A_obs, stack_obs,
     file = "cov_strat_LonoB_SPDE_OD_fit.RData")

##################################
# cov_nostrat_Binom_SPDE_TS
##################################

cov_nostrat_Binom_SPDE_TS <- y ~ -1 + beta0 +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde) +
  f(c_id, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

cov_nostrat_Binom_SPDE_TS_fit<- inla(formula = cov_nostrat_Binom_SPDE_TS,
                                      family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                      data = inla.stack.data(stack_obs),
                                      control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                      control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                      verbose = F)

save(cov_nostrat_Binom_SPDE_TS_fit, nga_spde, A_obs, stack_obs,
     file = "cov_nostrat_Binom_SPDE_TS_fit.RData")

##################################
# cov_strat_Binom_SPDE_TS
##################################

cov_strat_Binom_SPDE_TS <- y ~ -1 + beta0 +
  IsUrban +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde) +
  f(c_id, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

cov_strat_Binom_SPDE_TS_fit <- inla(formula = cov_strat_Binom_SPDE_TS,
                                    family = "binomial", Ntrials = inla.stack.data(stack_obs)$n,
                                    data = inla.stack.data(stack_obs),
                                    control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                    control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                    verbose = F)

save(cov_strat_Binom_SPDE_TS_fit, nga_spde, A_obs, stack_obs,
     file = "cov_strat_Binom_SPDE_TS_fit.RData")

##################################
# cov_nostrat_BetaB_SPDE_OD
##################################

cov_nostrat_BetaB_SPDE_OD <- y ~ -1 + beta0 +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde)

cov_nostrat_BetaB_SPDE_OD_fit <- inla(formula = cov_nostrat_BetaB_SPDE_OD,
                                      family = "betabinomial", Ntrials = inla.stack.data(stack_obs)$n,
                                      data = inla.stack.data(stack_obs),
                                      control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                      control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                      verbose = F)

save(cov_nostrat_BetaB_SPDE_OD_fit, nga_spde, A_obs, stack_obs,
     file = "cov_nostrat_BetaB_SPDE_OD_fit.RData")


##################################
# cov_strat_BetaB_SPDE_OD
##################################

cov_strat_BetaB_SPDE_OD <- y ~ -1 + beta0 +
  IsUrban +
  aridity10k + poverty + evi1k + log_travel_time + log_night +
  f(u, model = nga_spde)

cov_strat_BetaB_SPDE_OD_fit <- inla(formula = cov_strat_BetaB_SPDE_OD,
                                      family = "betabinomial", Ntrials = inla.stack.data(stack_obs)$n,
                                      data = inla.stack.data(stack_obs),
                                      control.predictor = list(A = inla.stack.A(stack_obs), compute = T, link = 1),
                                      control.compute = list(cpo = T, dic = T, config = T, waic = T),
                                      verbose = F)

save(cov_strat_BetaB_SPDE_OD_fit, nga_spde, A_obs, stack_obs,
     file = "cov_strat_BetaB_SPDE_OD_fit.RData")
