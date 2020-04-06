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

##################################
# load previous results
##################################

load("data_cleaned_1k.RData")

##################################
# maps for models
##################################

pred_dt_pixel_all <- pred_dt_lga_all <- pred_dt_state_all <- pred_dt_nation_all<- NULL

model_vt <- c("cov_strat_LonoB_SPDE_OD",
              "cov_strat_Binom_SPDE_TS",
              "cov_strat_BetaB_SPDE_OD", 
              "cov_strat_Binom_SPDE_NE")

for (m in model_vt){
  
  load(paste0(m, "_pred.RData"))
  
  # med color scheme
  col_regions_med <- colorRampPalette(brewer.pal(8, "RdYlGn"))(1000)
  at_med <- seq(0, 1, 0.001)
  labelat_med <- seq(0, 1, 0.2)
  labeltext_med <- format(round(labelat_med*100, 0), nsmall = 0)

  # width color scheme
  col_regions_wid1 <- colorRampPalette(rev(brewer.pal(8, "Blues")))(201)
  col_regions_wid2 <- colorRampPalette(brewer.pal(8, "Reds"))(401)
  col_regions_wid <- c(col_regions_wid1, col_regions_wid2)
  at_wid <- seq(-0.001, 0.601, 0.001)
  labelat_wid <- seq(0, 0.6, 0.1)
  labeltext_wid <- c(format(round(labelat_wid*100, 0), nsmall = 0)[1:6], ">=60")
  
  # distance scale
  scale_x <- min(nga_boundary$loc[, 1]) + 0.75*max(nga_boundary$loc[, 1])-min(nga_boundary$loc[, 1])
  scale_y <- min(nga_boundary$loc[, 2]) + 0.4*max(nga_boundary$loc[, 2])-min(nga_boundary$loc[, 2])
  scale_bar <- list("SpatialPolygonsRescale", layout.scale.bar(),
                    offset = c(scale_x, scale_y), scale = 2/111*100, 
                    fill = c("transparent","black"))
  text_x <- list("sp.text", c(scale_x, scale_y-0.4), "0", cex = 0.75)
  text_y <- list("sp.text", c(scale_x+2, scale_y-0.4), "200 km", cex = 0.75)
  
  #################################
  # state maps
  #################################
  
  # posterior med
  shp_plot <- merge(map_nga_state_shp, pred_dt_state, by = "StateID")

  bitmap(paste0("plots_pred/", m, "_med_state.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "med",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_med,
                    at = at_med,
                    colorkey = list(col = col_regions_med,
                                    at = at_med,
                                    labels = list(at = labelat_med, labels = labeltext_med)))
  print(sp_plot)
  grid.text("% vaccinated", x = unit(0.88, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()

  # posterior wid90
  shp_plot$wid90 <- shp_plot$up_90 - shp_plot$low_90
  shp_plot$wid90 <- ifelse(shp_plot$wid90 < 0.6, shp_plot$wid90, 0.6)

  bitmap(paste0("plots_pred/", m, "_wid90_state.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "wid90",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("width (%)", x = unit(0.91, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()

  #################################
  # lga maps
  #################################

  # posterior med
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")

  bitmap(paste0("plots_pred/", m, "_med_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "med",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_med,
                    at = at_med,
                    colorkey = list(col = col_regions_med,
                                    at = at_med,
                                    labels = list(at = labelat_med, labels = labeltext_med)))
  print(sp_plot)
  grid.text("% vaccinated", x = unit(0.88, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()

  # posterior wid90
  shp_plot$wid90 <- shp_plot$up_90 - shp_plot$low_90
  shp_plot$wid90 <- ifelse(shp_plot$wid90 < 0.6, shp_plot$wid90, 0.6)

  bitmap(paste0("plots_pred/", m, "_wid90_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "wid90",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("width (%)", x = unit(0.91, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  #################################
  # pixel maps
  #################################
  
  # posterior med
  shp_plot <- SpatialPixelsDataFrame(points = cbind(pred_dt_pixel[pop_all_wp > 0, long], 
                                                    pred_dt_pixel[pop_all_wp > 0, lat]), 
                                     data = pred_dt_pixel[pop_all_wp > 0, ],
                                     tolerance = 0.56)
  
  bitmap(paste0("plots_pred/", m, "_med_pixel.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "med",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_med,
                    at = at_med,
                    colorkey = list(col = col_regions_med,
                                    at = at_med,
                                    labels = list(at = labelat_med, labels = labeltext_med)))
  print(sp_plot)
  grid.text("% vaccinated", x = unit(0.88, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  # posterior wid90
  shp_plot$wid90 <- shp_plot$up_90 - shp_plot$low_90
  shp_plot$wid90 <- ifelse(shp_plot$wid90 < 0.6, shp_plot$wid90, 0.6)
  
  bitmap(paste0("plots_pred/", m, "_wid90_pixel.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "wid90",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("width (%)", x = unit(0.91, "npc"), y = unit(0.95, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  print(m)
}


