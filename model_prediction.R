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

nSamp <- 1000

##################################
# load previous results
##################################

load("data_cleaned_1k.RData")
load("cov_strat_LonoB_SPDE_OD_fit.RData")
inlafit <- cov_strat_LonoB_SPDE_OD_fit

# for lga and state level estimates
pred_dt_comp <- pred_dt_wp

pred_dt_comp[,  c("pop_nation_u5") := list(sum(pop_u5))]
pred_dt_comp[,  c("pop_state_u5") := list(sum(pop_u5)), by = StateID]
pred_dt_comp[,  c("pop_lga_u5") := list(sum(pop_u5)), by = LGAID]

pred_dt_comp[,  c("q_nation_u5") := list(pop_u5/pop_nation_u5)]
pred_dt_comp[,  c("q_state_u5") := list(pop_u5/pop_state_u5)]
pred_dt_comp[,  c("q_lga_u5") := list(pop_u5/pop_lga_u5)]

pred_dt_comp <- pred_dt_comp[order(StateID, LGAID, long, lat), ]

n_pred <- nrow(pred_dt_comp)

# Indices for different effects
idx_u <- inlafit$misc$configs$contents$start[3]
len_u <- inlafit$misc$configs$contents$length[3]
idx_X <- inlafit$misc$configs$contents$start[5]
X <- as.matrix(cbind(rep(1, n_pred), 
                     pred_dt_comp[, inlafit$misc$configs$contents$tag[6:11], with = F]))
len_X <- ncol(X)

# for pixel level estimates
pred_dt_pixel <- pred_dt_grid
pred_dt_pixel[, "IsUrban" := IsUrban_wp]
pred_dt_pixel <- pred_dt_pixel[order(StateID, LGAID, long, lat), ]
n_pred_pixel <- nrow(pred_dt_pixel)
X_pixel <- as.matrix(cbind(rep(1, n_pred_pixel), 
                           pred_dt_pixel[, inlafit$misc$configs$contents$tag[6:11], with = F]))

# Draw all samples
set.seed(12345)
inla.seed <- as.integer(runif(1)*.Machine$integer.max)
allSamples <- inla.posterior.sample(n = nSamp, inlafit, seed = inla.seed)

# Extract samples
xLatent <- matrix(0, nrow = length(allSamples[[1]]$latent), ncol = nSamp)
xHyperpar <- matrix(0, nrow = length(allSamples[[1]]$hyperpar), ncol = nSamp)
for(i in 1:nSamp){
  xLatent[, i] <- allSamples[[i]]$latent
  xHyperpar[, i] <- allSamples[[i]]$hyperpar
}

xX <- xLatent[idx_X:(idx_X + len_X - 1), ]
xU <- xLatent[idx_u:(idx_u + len_u - 1), ]
xSig <- sqrt(1/xHyperpar[3, ])
k <- 16*sqrt(3)/15/pi

# Iterate through samples and construct predictions for each sample of latent field
nation_dt <- data.table(mean = 0, sd = 0, med = 0,
                        low_95 = 0, up_95 = 0,
                        low_90 = 0, up_90 = 0,
                        low_80 = 0, up_80 = 0,
                        low_50 = 0, up_50 = 0)

##
state_dt <- unique(pred_dt_comp[, c("StateID"), with = F])
state_dt <- state_dt[order(StateID)]
state_dt[, c("mean", "sd", "med", 
             "low_95", "up_95",
             "low_90", "up_90",
             "low_80", "up_80",
             "low_50", "up_50")] <- 0
state_dt[, "row_num"] <- 1:nrow(state_dt)


##
state_lga_dt <- unique(pred_dt_comp[, c("StateID", "LGAID"), with = F])
state_lga_dt <- state_lga_dt[order(StateID, LGAID)]
state_lga_dt[, c("mean", "sd", "med",  
                 "low_95", "up_95",
                 "low_90", "up_90",
                 "low_80", "up_80",
                 "low_50", "up_50")] <- 0
state_lga_dt[, "row_num"] <- 1:nrow(state_lga_dt)

##
pixel_p_quantile_mt <- matrix(0, nrow = n_pred_pixel, ncol = 11)

##
q_nation_vt_u5 <- pred_dt_comp$q_nation_u5
q_state_vt_u5 <- pred_dt_comp$q_state_u5
q_lga_vt_u5 <- pred_dt_comp$q_lga_u5

##
A_pred <- inla.spde.make.A(mesh = nga_mesh,
                           loc = cbind(pred_dt_comp$long, pred_dt_comp$lat))
A_pred_pixel <- inla.spde.make.A(mesh = nga_mesh,
                           loc = cbind(pred_dt_pixel$long, pred_dt_pixel$lat))

postSamples_p_weighted <-  NULL

pixel_postsamp_mt <- matrix(0, nrow = nrow(pred_dt_pixel), ncol = nSamp)
lga_postsamp_mt <- matrix(0, nrow = nrow(state_lga_dt), ncol = nSamp)
state_postsamp_mt <- matrix(0, nrow = nrow(state_dt), ncol = nSamp)

for (i in unique(state_dt[, StateID])){
  # i <- unique(state_dt[, StateID])[2]
  
  postSamples_p_i_weighted <-  NULL
  
  for (j in unique(state_lga_dt[StateID == i, LGAID])){
    # j <- unique(state_lga_dt[StateID == i, LGAID])[2]
    
    # pixel level
    idx_pred_j_pixel <- which(pred_dt_pixel[, LGAID] == j)
    
    postSamples_p_j_pixel <- expit((as.matrix(X_pixel[idx_pred_j_pixel, ] %*% xX) + 
                                as.matrix(A_pred_pixel[idx_pred_j_pixel, ] %*% xU))/
                               (matrix(rep(sqrt(1+k^2*xSig^2), length(idx_pred_j_pixel)), 
                                      nrow = length(idx_pred_j_pixel), byrow = T)))
    
    pixel_p_quantile_mt[idx_pred_j_pixel, 1] <- t(apply(X = postSamples_p_j_pixel, 
                                                  MARGIN = 1, 
                                                  FUN = mean, na.rm = T))
    pixel_p_quantile_mt[idx_pred_j_pixel, 2] <- t(apply(X = postSamples_p_j_pixel,
                                                  MARGIN = 1, 
                                                  FUN = sd, na.rm = T))
    pixel_p_quantile_mt[idx_pred_j_pixel, 3:11] <- t(apply(X = postSamples_p_j_pixel, 
                                                     MARGIN = 1, 
                                                     FUN = quantile, c(0.5, 
                                                                       0.025, 0.975,
                                                                       0.05, 0.95,
                                                                       0.1, 0.9, 
                                                                       0.25, 0.75), 
                                                     na.rm = T))
    pixel_postsamp_mt[idx_pred_j_pixel, ] <- postSamples_p_j_pixel
    
    # areal level
    idx_pred_j <- which(pred_dt_comp[, LGAID] == j)
    
    postSamples_p_j <- expit((as.matrix(X[idx_pred_j, ] %*% xX) + 
                                as.matrix(A_pred[idx_pred_j, ] %*% xU))/
                               (matrix(rep(sqrt(1+k^2*xSig^2), length(idx_pred_j)), 
                                       nrow = length(idx_pred_j), byrow = T)))
    
    # lga level
    q_lga_mt_j <- matrix(rep(q_lga_vt_u5[idx_pred_j], nSamp), ncol = nSamp, byrow = F)
    
    postSamples_pq_j <- postSamples_p_j * q_lga_mt_j
    postSamples_p_j_weighted <- apply(postSamples_pq_j, MARGIN = 2, sum, na.rm = T)
    
    state_lga_dt[LGAID == j, "mean"] <- mean(postSamples_p_j_weighted, na.rm = T)
    state_lga_dt[LGAID == j, "sd"] <- sd(postSamples_p_j_weighted, na.rm = T)
    state_lga_dt[LGAID == j, "med"] <- quantile(postSamples_p_j_weighted, 0.5, na.rm = T)
    
    state_lga_dt[LGAID == j, "low_95"] <- quantile(postSamples_p_j_weighted, 0.025, na.rm = T)
    state_lga_dt[LGAID == j, "up_95"] <- quantile(postSamples_p_j_weighted, 0.975, na.rm = T)
    
    state_lga_dt[LGAID == j, "low_90"] <- quantile(postSamples_p_j_weighted, 0.05, na.rm = T)
    state_lga_dt[LGAID == j, "up_90"] <- quantile(postSamples_p_j_weighted, 0.95, na.rm = T)
    
    state_lga_dt[LGAID == j, "low_80"] <- quantile(postSamples_p_j_weighted, 0.1, na.rm = T)
    state_lga_dt[LGAID == j, "up_80"] <- quantile(postSamples_p_j_weighted, 0.9, na.rm = T)
    
    state_lga_dt[LGAID == j, "low_50"] <- quantile(postSamples_p_j_weighted, 0.25, na.rm = T)
    state_lga_dt[LGAID == j, "up_50"] <- quantile(postSamples_p_j_weighted, 0.75, na.rm = T)
    
    row_num_j <- state_lga_dt[LGAID == j, row_num]
    lga_postsamp_mt[row_num_j, ] <- postSamples_p_j_weighted
    
    # state level
    q_state_mt_j <- matrix(rep(q_state_vt_u5[idx_pred_j], nSamp), ncol = nSamp, byrow = F)
    
    postSamples_pq_j <- postSamples_p_j * q_state_mt_j
    postSamples_p_j_weighted <- apply(postSamples_pq_j, MARGIN = 2, sum, na.rm = T)
    postSamples_p_i_weighted <- rbind(postSamples_p_i_weighted, 
                                      postSamples_p_j_weighted)
    
    # nation level
    q_nation_mt_j <- matrix(rep(q_nation_vt_u5[idx_pred_j], nSamp), ncol = nSamp, byrow = F)
    
    postSamples_pq_j <- postSamples_p_j * q_nation_mt_j
    postSamples_p_j_weighted <- apply(postSamples_pq_j, MARGIN = 2, sum, na.rm = T)
    postSamples_p_weighted <- rbind(postSamples_p_weighted, 
                                    postSamples_p_j_weighted)
    
    print(j)
  }
  
  # state level
  postSamples_p_i_weighted <- apply(postSamples_p_i_weighted, MARGIN = 2, sum, na.rm = T)
  
  state_dt[StateID == i, "mean"] <- mean(postSamples_p_i_weighted, na.rm = T)
  state_dt[StateID == i, "sd"] <- sd(postSamples_p_i_weighted, na.rm = T)
  state_dt[StateID == i, "med"] <- quantile(postSamples_p_i_weighted, 0.5, na.rm = T)
  
  state_dt[StateID == i, "low_95"] <- quantile(postSamples_p_i_weighted, 0.025, na.rm = T)
  state_dt[StateID == i, "up_95"] <- quantile(postSamples_p_i_weighted, 0.975, na.rm = T)
  
  state_dt[StateID == i, "low_90"] <- quantile(postSamples_p_i_weighted, 0.05, na.rm = T)
  state_dt[StateID == i, "up_90"] <- quantile(postSamples_p_i_weighted, 0.95, na.rm = T)
  
  state_dt[StateID == i, "low_80"] <- quantile(postSamples_p_i_weighted, 0.1, na.rm = T)
  state_dt[StateID == i, "up_80"] <- quantile(postSamples_p_i_weighted, 0.9, na.rm = T)
  
  state_dt[StateID == i, "low_50"] <- quantile(postSamples_p_i_weighted, 0.25, na.rm = T)
  state_dt[StateID == i, "up_50"] <- quantile(postSamples_p_i_weighted, 0.75, na.rm = T)
  
  row_num_i <- state_dt[StateID == i, row_num]
  state_postsamp_mt[row_num_i, ] <- postSamples_p_i_weighted
}

# nation level
postSamples_p_weighted <- apply(postSamples_p_weighted, MARGIN = 2, sum, na.rm = T)

nation_dt[, "mean"] <- mean(postSamples_p_weighted, na.rm = T)
nation_dt[, "sd"] <- sd(postSamples_p_weighted, na.rm = T)
nation_dt[, "med"] <- quantile(postSamples_p_weighted, 0.5, na.rm = T)

nation_dt[, "low_95"] <- quantile(postSamples_p_weighted, 0.025, na.rm = T)
nation_dt[, "up_95"] <- quantile(postSamples_p_weighted, 0.975, na.rm = T)

nation_dt[, "low_90"] <- quantile(postSamples_p_weighted, 0.05, na.rm = T)
nation_dt[, "up_90"] <- quantile(postSamples_p_weighted, 0.95, na.rm = T)

nation_dt[, "low_80"] <- quantile(postSamples_p_weighted, 0.1, na.rm = T)
nation_dt[, "up_80"] <- quantile(postSamples_p_weighted, 0.9, na.rm = T)

nation_dt[, "low_50"] <- quantile(postSamples_p_weighted, 0.25, na.rm = T)
nation_dt[, "up_50"] <- quantile(postSamples_p_weighted, 0.75, na.rm = T)

# pixel level
pred_dt_pixel[, c("mean")] <- pixel_p_quantile_mt[, 1]
pred_dt_pixel[, c("sd")] <- pixel_p_quantile_mt[, 2]
pred_dt_pixel[, c("med")] <- pixel_p_quantile_mt[, 3]

pred_dt_pixel[, c("low_95")] <- pixel_p_quantile_mt[, 4]
pred_dt_pixel[, c("up_95")] <- pixel_p_quantile_mt[, 5]

pred_dt_pixel[, c("low_90")] <- pixel_p_quantile_mt[, 6]
pred_dt_pixel[, c("up_90")] <- pixel_p_quantile_mt[, 7]

pred_dt_pixel[, c("low_80")] <- pixel_p_quantile_mt[, 8]
pred_dt_pixel[, c("up_80")] <- pixel_p_quantile_mt[, 9]

pred_dt_pixel[, c("low_50")] <- pixel_p_quantile_mt[, 10]
pred_dt_pixel[, c("up_50")] <- pixel_p_quantile_mt[, 11]

# save results
pred_dt_nation <- nation_dt
pred_dt_state <- state_dt
pred_dt_lga <- state_lga_dt
pred_dt_pixel <- pred_dt_pixel

save(pred_dt_nation, pred_dt_state, pred_dt_lga, pred_dt_pixel,
     state_postsamp_mt, lga_postsamp_mt, pixel_postsamp_mt, 
     file = "cov_strat_LonoB_SPDE_OD_pred.RData")
