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
  
  # add hatching 
  get_odds <- function(p){p/(1-p)}
  get_cv <- function(x_vt){sqrt(var(x_vt))/mean(x_vt)}
  pred_dt_state[, "odds_cv"] <- apply(get_odds(state_postsamp_mt), 1, get_cv)
  
  #################################
  # state maps
  #################################

  # posterior med

  bitmap(paste0("plots_hatch/", m, "_med_state_hatch.jpeg"),
         width = 6, height = 3.2, units = "in", type = "jpeg", res = 1000)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(3, 2), heights = c(3))
  par(mai = c(0.25, 0.1, 0.3, 0.1), oma = c(0.5, 0.1, 0.1, 0.1))
  # map
  med.palette <- col_regions_med
  zlim <- range(pred_dt_state$med, na.rm = TRUE)
  pred_dt_state$value.col <- med.palette[floor((pred_dt_state$med - zlim[1])/(zlim[2] - zlim[1]) * (1000 - 1)) + 1]
  tmp <- pred_dt_state
  tmp <- tmp[match(map_nga_state_shp[["StateID"]], tmp[, StateID]), ]
  sp::plot(map_nga_state_shp, col = tmp$value.col, border = "black", 
           main = "", lwd = 1)
  map_nga_state_shp$odds_cv <- tmp[, odds_cv]
  rrt1 <- map_nga_state_shp@data$odds_cv
  breaks_vt <- seq(min(pred_dt_state[, odds_cv], na.rm = TRUE), 
                   max(pred_dt_state[, odds_cv], na.rm = TRUE), 
                   length.out = 6)
  nbrks <- length(breaks_vt)
  brklabels <- paste(signif(breaks_vt[1:(nbrks - 1)], 2), signif(breaks_vt[2:nbrks], 2), sep = " - ")
  dens <- (2:nbrks) * 3
  sp::plot(map_nga_state_shp, 
           density = dens[findInterval(pred_dt_state[, odds_cv], breaks_vt, all.inside = TRUE)], 
           add = T, col = "black", border = FALSE)
  
  # legend
  legend.col <- function(col, lev, hadj = -1.5, title = NULL) {
    opar <- par
    n <- length(col)
    bx <- graphics::par("usr")
    box.cx <- c(bx[1] + (bx[2] - bx[1])/1000, bx[1] + (bx[2] - bx[1])/1000 + (bx[2] - bx[1])/30)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3])/n
    xx <- rep(box.cx, each = 2)
    graphics::par(xpd = TRUE)
    for (i in 1:n) {
      yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + 
                (box.sy * (i)), box.cy[1] + (box.sy * (i)), box.cy[1] + 
                (box.sy * (i - 1)))
      graphics::polygon(xx, yy, col = col[i], border = col[i])
    }
    if (!is.null(title)) 
    graphics::text(box.cx[1] - 0.05, box.cy[1] + box.sy * (n + 2) + 0.05, title, pos = 4, cex = 1.5)
    graphics::par(new = TRUE)
    graphics::plot(0, 0, type = "n", ylim = c(min(lev), max(lev)), 
                   yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
    graphics::axis(side = 2, las = 2, tick = FALSE, line = 0.25, cex.axis = 1.5,
                   hadj = hadj)
    par <- opar
  }
  
  graphics::plot(1, type = "n", axes = FALSE, xlab = "",  ylab = "")
  legend.col(col = med.palette, lev = c(0, 1), hadj = -2, title = "MCV1 Posterior Median")
  graphics::legend(x = "right", inset = 0, legend = brklabels, 
                   col = rep("black", 2), cex = 2.2, density = dens, bty = "n", 
                   title = "CV of MCV1 Odds")
  
  dev.off()

  print(m)
}



