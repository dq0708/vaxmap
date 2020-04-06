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

  # width color scheme
  col_regions_wid <- colorRampPalette(brewer.pal(8, "RdYlGn"))(1002)
  at_wid <- seq(-0.001, 1.001, 0.001)
  labelat_wid <- seq(0, 1, 0.2)
  labeltext_wid <- format(round(labelat_wid, 1), nsmall = 1)
  
  # distance scale
  scale_x <- min(nga_boundary$loc[, 1]) + 0.75*max(nga_boundary$loc[, 1])-min(nga_boundary$loc[, 1])
  scale_y <- min(nga_boundary$loc[, 2]) + 0.4*max(nga_boundary$loc[, 2])-min(nga_boundary$loc[, 2])
  scale_bar <- list("SpatialPolygonsRescale", layout.scale.bar(),
                    offset = c(scale_x, scale_y), scale = 2/111*100, 
                    fill = c("transparent","black"))
  text_x <- list("sp.text", c(scale_x, scale_y-0.4), "0", cex = 0.75)
  text_y <- list("sp.text", c(scale_x+2, scale_y-0.4), "200 km", cex = 0.75)

  #################################
  # lga maps
  #################################

  pred_dt_lga[, "above80"] <- apply(lga_postsamp_mt, 1, function(x){mean(x>=0.8)})
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")

  bitmap(paste0("plots_pred80/", m, "_above80_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "above80",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("Pr(MCV1>=80%)", x = unit(0.85, "npc"), y = unit(0.97, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  ######################
  
  pred_dt_lga[, "above60"] <- apply(lga_postsamp_mt, 1, function(x){mean(x>=0.6)})
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")
  
  bitmap(paste0("plots_pred80/", m, "_above60_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "above60",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("Pr(MCV1>=60%)", x = unit(0.85, "npc"), y = unit(0.97, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  ######################
  
  pred_dt_lga[, "above50"] <- apply(lga_postsamp_mt, 1, function(x){mean(x>=0.5)})
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")
  
  bitmap(paste0("plots_pred80/", m, "_above50_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "above50",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("Pr(MCV1>=50%)", x = unit(0.85, "npc"), y = unit(0.97, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  ######################
  
  pred_dt_lga[, "above40"] <- apply(lga_postsamp_mt, 1, function(x){mean(x>=0.4)})
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")
  
  bitmap(paste0("plots_pred80/", m, "_above40_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "above40",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("Pr(MCV1>=40%)", x = unit(0.85, "npc"), y = unit(0.97, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  ######################
  
  pred_dt_lga[, "above20"] <- apply(lga_postsamp_mt, 1, function(x){mean(x>=0.2)})
  shp_plot <- merge(map_nga_lga_shp, pred_dt_lga, by = "LGAID")
  
  # sum(pred_dt_lga[, above20] < 0.5)/774 # 28, 0.06847545
  
  bitmap(paste0("plots_pred80/", m, "_above20_lga.jpeg"),
         width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
  sp_plot <- spplot(shp_plot, "above20",
                    main = "", xlab = "", ylab = "",
                    sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions_wid,
                    at = at_wid,
                    colorkey = list(col = col_regions_wid,
                                    at = at_wid,
                                    labels = list(at = labelat_wid, labels = labeltext_wid)))
  print(sp_plot)
  grid.text("Pr(MCV1>=20%)", x = unit(0.85, "npc"), y = unit(0.97, "npc"),
            gp = gpar(fontsize = 10))
  dev.off()
  
  print(m)
}