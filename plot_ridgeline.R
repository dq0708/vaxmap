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
library(rasterVis)
library(plotrix)
library(ggplot2)
library(ggridges)

##################################
# load previous results
##################################

load("data_cleaned_1k.RData")

##################################
# maps for models
##################################

model_vt <- c("cov_strat_LonoB_SPDE_OD",
              "cov_strat_Binom_SPDE_TS",
              "cov_strat_BetaB_SPDE_OD",
              "cov_strat_Binom_SPDE_NE")

for (m in model_vt){
  
  load(paste0(m, "_pred.RData"))
  
  ##########################################
  # ridge State post med
  ##########################################
  
  pred_dt_state_order <- merge(pred_dt_state, map_nga_state_shp@data, by = "StateID")
  pred_dt_state_order <- pred_dt_state_order[order(-med)]

  data_plot_dt <- data.table(StateID = rep(map_nga_state_shp$StateID, ncol(state_postsamp_mt)),
                             StateName = rep(map_nga_state_shp$StateName, ncol(state_postsamp_mt)),
                             MCV1 = 100*as.numeric(state_postsamp_mt))
  data_plot_dt[, "State"] <- factor(data_plot_dt[, StateName], levels = rev(pred_dt_state_order[, StateName]))

  # plot
  bitmap(paste0("plots_ridge/", m, "_all_state_med.jpeg"),
         width = 8, height = 8, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  p <- ggplot(data_plot_dt, aes(x = MCV1, y = State)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
    scale_fill_gradientn(colours = rev(c("#1a9641", "#a6d96a", "#fdae61", "#d7191c")), guide = F) + 
    labs(x = "MCV1 Coverage (%)", y = "")
  print(p)
  
  dev.off()
  
  ##########################################
  # ridge State post mean
  ##########################################
  
  pred_dt_state_order <- merge(pred_dt_state, map_nga_state_shp@data, by = "StateID")
  pred_dt_state_order <- pred_dt_state_order[order(-mean)]
  
  data_plot_dt <- data.table(StateID = rep(map_nga_state_shp$StateID, ncol(state_postsamp_mt)),
                             StateName = rep(map_nga_state_shp$StateName, ncol(state_postsamp_mt)),
                             MCV1 = 100*as.numeric(state_postsamp_mt))
  data_plot_dt[, "State"] <- factor(data_plot_dt[, StateName], levels = rev(pred_dt_state_order[, StateName]))
  
  # plot
  bitmap(paste0("plots_ridge/", m, "_all_state_mean.jpeg"),
         width = 8, height = 8, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  p <- ggplot(data_plot_dt, aes(x = MCV1, y = State)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
    scale_fill_gradientn(colours = rev(c("#1a9641", "#a6d96a", "#fdae61", "#d7191c")), guide = F) + 
    labs(x = "MCV1 Coverage (%)", y = "")
  print(p)
  
  dev.off()
  
  print(m)
}
