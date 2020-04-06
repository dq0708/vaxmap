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

##################################
# load previous results
##################################

load("cov_strat_LonoB_SPDE_OD_pred.RData")

Nclass=4

##################################
# define function to get measure
##################################

get_measure_dt <- function(postsamp_mt, pred_dt, grp_thresh){
  # postsamp_mt <- state_postsamp_mt
  # pred_dt <- pred_dt_state
  # grp_thresh <- c(0, 0.5, 1)
  
  # create group lookup table based on thresholds
  n_grp <- length(grp_thresh)-1
  grp_dt <- data.table(grp = 1:n_grp,
                       lower = grp_thresh[1:(n_grp)],
                       upper = grp_thresh[2:(n_grp+1)])
  
  n_area <- nrow(postsamp_mt)
  n_postsamp <- ncol(postsamp_mt)
  
  ##################################
  # assign group based on model max post prob
  ##################################
  
  grp_cnt_mt <- matrix(0, nrow = n_area, ncol = n_grp)
  for (i in 1:n_grp){
    grp_cnt_mt[, i] <- apply(postsamp_mt, 1, 
                             function(x){sum(x > grp_dt[i, lower] & x <= grp_dt[i, upper])})
  }
  
  DF <- data.frame(grp_cnt_mt)
  DT <- data.table(value = unlist(DF, use.names=FALSE), 
                   colid = 1:nrow(DF), 
                   rowid = rep(names(DF), each=nrow(DF)))
  setkey(DT, colid, value)
  
  grp_cnt_dt_max <- as.data.table(DF)
  grp_cnt_dt_max[, "grp"] <- DT[J(unique(colid), DT[J(unique(colid)), value, mult="last"]), rowid, mult="first"]
  grp_cnt_dt_max[, "TCP" := 0]
  for (i in 1:n_grp){
    idx_grp <- which(grp_cnt_dt_max[, grp] == paste0("X", i))
    grp_cnt_dt_max[idx_grp, "TCP"] <- grp_cnt_dt_max[idx_grp, paste0("X", i), with = F]/n_postsamp
  }
  
  pred_dt[, "grp"] <- grp_cnt_dt_max[, grp]
  pred_dt[, "TCP"] <- grp_cnt_dt_max[, TCP]

  return(pred_dt)
}

##################################
# analysis for Nclass
##################################

L_vt_state <- c(0, 0.2, 0.5, 0.8, 1)
L_val_vt_state <- c(0.1, 0.35, 0.65, 0.9)
L_dt_state <- data.table(grp = paste0("X", 1:Nclass),
                         grp_low = L_vt_state[1:Nclass],
                         grp_up = L_vt_state[2:(Nclass+1)],
                         grp_val = L_val_vt_state)

L_vt_lga <- c(0, 0.2, 0.5, 0.8, 1)
L_val_vt_lga <- c(0.1, 0.35, 0.65, 0.9)
L_dt_lga <- data.table(grp = paste0("X", 1:Nclass),
                       grp_low = L_vt_lga[1:Nclass],
                       grp_up = L_vt_lga[2:(Nclass+1)],
                       grp_val = L_val_vt_lga)

L_vt_pixel <- c(0, 0.2, 0.5, 0.8, 1)
L_val_vt_pixel <- c(0.1, 0.35, 0.65, 0.9)
L_dt_pixel <- data.table(grp = paste0("X", 1:Nclass),
                         grp_low = L_vt_pixel[1:Nclass],
                         grp_up = L_vt_pixel[2:(Nclass+1)],
                         grp_val = L_val_vt_pixel)

##################################
# state
##################################

t1 <- Sys.time()
measure_state_dt <- get_measure_dt(postsamp_mt = state_postsamp_mt,
                                   pred_dt = pred_dt_state,
                                   grp_thresh = L_vt_state)
measure_state_dt <- merge(measure_state_dt, L_dt_state, by = "grp", all.x = T)
measure_state_dt[, "Nclass"] <- Nclass
t2 <- Sys.time()
print(paste0("State ", round(t2-t1, 3), " sec"))

save(measure_state_dt, L_dt_state,
     file = paste0("cov_strat_LonoB_SPDE_OD_measure258_state_Nclass", Nclass, ".RData"))

##################################
# lga
##################################
t1 <- Sys.time()
measure_lga_dt <- get_measure_dt(postsamp_mt = lga_postsamp_mt,
                                   pred_dt = pred_dt_lga,
                                   grp_thresh = L_vt_lga)
measure_lga_dt <- merge(measure_lga_dt, L_dt_lga, by = "grp", all.x = T)
measure_lga_dt[, "Nclass"] <- Nclass
t2 <- Sys.time()
print(paste0("LGA ", round(t2-t1, 3), " sec"))

save(measure_lga_dt, L_dt_lga,
     file = paste0("cov_strat_LonoB_SPDE_OD_measure258_lga_Nclass", Nclass, ".RData"))

##################################
# pixel
##################################
t1 <- Sys.time()
measure_pixel_dt <- get_measure_dt(postsamp_mt = pixel_postsamp_mt,
                                   pred_dt = pred_dt_pixel,
                                   grp_thresh = L_vt_pixel)
measure_pixel_dt <- merge(measure_pixel_dt, L_dt_pixel, by = "grp", all.x = T)
measure_pixel_dt[, "Nclass"] <- Nclass
t2 <- Sys.time()
print(paste0("pixel ", round(t2-t1, 3), " sec"))

save(measure_pixel_dt, L_dt_pixel,
     file = paste0("cov_strat_LonoB_SPDE_OD_measure258_pixel_Nclass", Nclass, ".RData"))

