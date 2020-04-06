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

model_vt <- c("cov_strat_LonoB_SPDE_OD",
              "cov_strat_Binom_SPDE_TS",
              "cov_strat_BetaB_SPDE_OD",
              "cov_strat_Binom_SPDE_NE")

for (m in model_vt){
  
  ##################################
  # maps and histograms
  ##################################
  
  Nclass_vt <- c(2,3,4,5)
  
  for (K in Nclass_vt){
    # K <- Nclass_vt[2]
    load(paste0(m, "_measure1_state_Nclass", K, ".RData"))
    
    shp_plot <- merge(map_nga_state_shp, measure_state_dt, by = "StateID")
    
    manual.col <- colorRampPalette(brewer.pal(8, "RdYlGn"))(1000)
    color.match <- manual.col[round(sort(L_dt_state$grp_val)*1000)]
    lookup_dt <- data.table(grp_val = sort(L_dt_state$grp_val),
                            col = color.match)
    shp_plot <- merge(shp_plot, lookup_dt, by = "grp_val")
    shp_plot$grp_val <- as.factor(shp_plot$grp_val)
    col_regions <- as.vector(lookup_dt[grp_val %in% shp_plot$grp_val, col])
    
    labelat <- sort(unique(c(L_dt_state$grp_low, 
                             L_dt_state$grp_up)))
    labeltext <- format(round(labelat*100, 0), nsmall = 0)
    
    # distance scale
    scale_x <- min(nga_boundary$loc[, 1]) + 0.75*max(nga_boundary$loc[, 1])-min(nga_boundary$loc[, 1])
    scale_y <- min(nga_boundary$loc[, 2]) + 0.4*max(nga_boundary$loc[, 2])-min(nga_boundary$loc[, 2])
    scale_bar <- list("SpatialPolygonsRescale", layout.scale.bar(),
                      offset = c(scale_x, scale_y), scale = 2/111*100, 
                      fill = c("transparent","black"))
    text_x <- list("sp.text", c(scale_x, scale_y-0.4), "0", cex = 0.75)
    text_y <- list("sp.text", c(scale_x+2, scale_y-0.4), "200 km", cex = 0.75)
    
    #########
    # map
    #########
    
    bitmap(paste0("plots_measure1/plots_state/", m, "_map_Nclass", K, "_state.jpeg"),
           width = 4, height = 3.2, units = "in", type = "jpeg", res = 1000)
    sp_plot <- spplot(shp_plot, zcol = "grp_val",
                      main = paste0("State Map K = ", K, ", ATCP = ", format(round(mean(measure_state_dt$TCP), 2), nsmall = 2)),
                      xlab = "", ylab = "",
                      sp.layout = list(scale_bar, text_x, text_y),
                      col.regions = col_regions,
                      colorkey = list(col = color.match,
                                      at = labelat,
                                      labels = list(at = labelat, labels = labeltext)))
    print(sp_plot)
    grid.text("% vaccinated", x = unit(0.88, "npc"), y = unit(0.05, "npc"),
              gp = gpar(fontsize = 10))
    dev.off()
    
    #########
    # stack hist
    #########
    L_dt_state[, "grp_name" := paste0("(", format(round(grp_low*100, 0), nsmall = 0), "%, ", format(round(grp_up*100, 0), nsmall = 0), "%]")]
    measure_state_dt[, "grp_name" := paste0("(", format(round(grp_low*100, 0), nsmall = 0), "%, ", format(round(grp_up*100, 0), nsmall = 0), "%]")]
    measure_state_dt$grp_name <- factor(measure_state_dt$grp_name, rev(L_dt_state$grp_name))
    
    bitmap(paste0("plots_measure1/plots_state/", m, "_stack_hist_Nclass", K, "_state.jpeg"),
           width = 2.5, height = 2, units = "in", type = "jpeg", res = 1000)
    par(mar = c(4, 4, 2, 2))
    histStack(TCP ~ grp_name, 
              data = measure_state_dt,
              col = rev(color.match), 
              # main = paste0("State Map K = ", K, ", ATCP = ", format(round(mean(measure_state_dt$TCP), 2), nsmall = 2)),
              xlab = "TCP",
              legend.pos = "topleft", 
              breaks = seq(0,1, 0.05))
    abline(v = mean(measure_state_dt$TCP), col = "blue", lwd = 2)
    dev.off()
  }
  
  print(m)
}

