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

model_vt <- c("cov_strat_LonoB_SPDE_OD",
              "cov_strat_Binom_SPDE_TS",
              "cov_strat_BetaB_SPDE_OD",
              "cov_strat_Binom_SPDE_NE")

for (m in model_vt){
  
  load(paste0(m, "_pred.RData"))
  
  ##########################################
  # rank State post samples
  ##########################################
  
  state_rank_mt <- apply(-state_postsamp_mt, 2, rank, )
  pred_dt_state[, "avg_rank"] <- apply(state_rank_mt, 1, mean)
  pred_dt_state[, "low_rank"] <- apply(state_rank_mt, 1, min)
  pred_dt_state[, "up_rank"] <- apply(state_rank_mt, 1, max)
  
  # top 5 state
  bitmap(paste0("plots_rank/", m, "_top5_state.jpeg"),
         width = 2, height = 4, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(5, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_state_order <- pred_dt_state[order(avg_rank)]
  
  for (i in 1:5){
    # i <- 1
    
    state_id <- pred_dt_state_order[i, StateID]
    state_name <- map_nga_state_shp$StateName[state_id]
    
    state_rank_vt <- state_rank_mt[state_id, ]
    state_rank_vt <- ifelse(state_rank_vt <= 10, state_rank_vt, 10)
    
    avg_rank <- pred_dt_state_order[i, avg_rank]
    
    state_ranktable <- as.data.table(table(state_rank_vt))
    state_ranktable <- merge(data.table(rank = as.character(1:10)), state_ranktable, 
                             by.x = "rank", by.y = "state_rank_vt", all.x = T)
    state_ranktable[, "rank" := as.integer(rank)]
    state_ranktable <- state_ranktable[order(rank)]
    state_ranktable[is.na(N), "N"] <- 0
    
    barplot(state_ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(state_name, ", ER = ", format(round(avg_rank, 1), nsmall = 1)),
            xaxt = "n", yaxt = "n", col = "#31a354", border = F)
    axis(1, at = 10:1-0.5, labels = c("10+", as.character(9:1)), tick = F)
  }
  dev.off()
  
  # bottom 5 state
  bitmap(paste0("plots_rank/", m, "_bottom5_state.jpeg"),
         width = 2, height = 4, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(5, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_state_order <- pred_dt_state[order(-avg_rank)]
  
  for (i in 1:5){
    # i <- 1
    
    state_id <- pred_dt_state_order[i, StateID]
    state_name <- map_nga_state_shp$StateName[state_id]
    
    avg_rank <- pred_dt_state_order[i, avg_rank]
    
    state_rank_vt <- state_rank_mt[state_id, ]
    state_rank_vt <- ifelse(state_rank_vt >= (37-9), state_rank_vt, (37-9))
    
    state_ranktable <- as.data.table(table(state_rank_vt))
    state_ranktable <- merge(data.table(rank = as.character(37:(37-9))), state_ranktable, 
                             by.x = "rank", by.y = "state_rank_vt", all.x = T)
    state_ranktable[, "rank" := as.integer(rank)]
    state_ranktable <- state_ranktable[order(rank)]
    state_ranktable[is.na(N), "N"] <- 0
    
    barplot(state_ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(state_name, ", ER = ", format(round(avg_rank, 1), nsmall = 1)),
            xaxt = "n", yaxt = "n", col = "#fc4e2a", border = F)
    axis(1, at = 10:1-0.5, labels = c(as.character(37:(37-8)), "28-"), tick = F)
  }
  dev.off()
  
  # all states hist
  bitmap(paste0("plots_rank/", m, "_all_state.jpeg"),
         width = 15, height = 20, units = "in", type = "jpeg", res = 400)
  par(mar = c(2.5, 1, 2, 1), mfcol = c(13, 3))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_state_order <- pred_dt_state[order(-avg_rank)]
  
  for (i in 1:37){
    # i <- 1
    
    state_id <- pred_dt_state_order[i, StateID]
    state_name <- map_nga_state_shp$StateName[state_id]
    
    state_rank_vt <- state_rank_mt[state_id, ]
    
    avg_rank <- pred_dt_state_order[i, avg_rank]
    
    state_ranktable <- as.data.table(table(state_rank_vt))
    state_ranktable <- merge(data.table(rank = as.character(1:37)), state_ranktable, 
                             by.x = "rank", by.y = "state_rank_vt", all.x = T)
    state_ranktable[, "rank" := as.integer(rank)]
    state_ranktable <- state_ranktable[order(rank)]
    state_ranktable[is.na(N), "N"] <- 0
    
    barplot(state_ranktable$N, width = 0.825, 
            xlim = c(37, 0), xlab = "", ylab = "",
            main = paste0(state_name, ", ER = ", format(round(avg_rank, 1), nsmall = 1)),
            xaxt = "n", yaxt = "n", col = "#08519c", border = F)
    axis(1, at = 37:1-0.5, labels = as.character(37:1), tick = F)
  }
  dev.off()
  
  # all states line plots
  bitmap(paste0("plots_rank/", m, "_all_state_line.jpeg"),
         width = 4, height = 3, units = "in", type = "jpeg", res = 400)
  par(mar=c(5.5, 4, 1.5, 1), mfrow = c(1, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_state_order <- merge(pred_dt_state, 
                               map_nga_state_shp@data)
  pred_dt_state_order <- pred_dt_state_order[order(-avg_rank)]
  pred_dt_state_order[, "row_num"] <- 1:37
  
  plot(x = -100, y = -100, 
       xlim = c(0, 37), ylim = c(0, 37), 
       ann = FALSE, xaxt = "n", yaxt = "n")
  title(main = "Range of State Posterior Rankings", font.main = 2, cex.main = 0.9)
  mtext(side = 2, text = "rank", line = 2.7, font = 2, cex = 0.9)
  axis(1, at = 1:37, labels = pred_dt_state_order$StateName, 
       las = 2, cex.axis = 0.9, font = 2)
  axis(2, at = 1:37, labels = as.character(37:1), 
       las = 1, cex.axis = 0.9, font = 2)
  
  # 2015 estimates
  points(x = pred_dt_state_order[, row_num],
         y = 38-pred_dt_state_order[, avg_rank], 
         col = "#08519c", pch = 16, cex = 0.6)
  
  segments(pred_dt_state_order[, row_num], 
           38-pred_dt_state_order[, low_rank],
           pred_dt_state_order[, row_num], 
           38-pred_dt_state_order[, up_rank],
           col = "#08519c", lwd = 0.65)
  dev.off()
  
  ##########################################
  # rank LGA post samples
  ##########################################
  
  lga_rank_mt <- apply(-lga_postsamp_mt, 2, rank, )
  pred_dt_lga[, "avg_rank"] <- apply(lga_rank_mt, 1, mean)
  pred_dt_lga[, "low_rank"] <- apply(lga_rank_mt, 1, min)
  pred_dt_lga[, "up_rank"] <- apply(lga_rank_mt, 1, max)
  
  # top 5 LGA
  bitmap(paste0("plots_rank/", m, "_top5_lga.jpeg"),
         width = 2, height = 4, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(5, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_lga_order <- pred_dt_lga[order(avg_rank)]
  
  for (i in 1:5){
    # i <- 5
    lga_id <- pred_dt_lga_order[i, LGAID]
    lga_name <- map_nga_lga_shp$LGAName[lga_id]
    
    state_id <- pred_dt_lga_order[i, StateID]
    state_name <- map_nga_lga_shp$StateName[lga_id]
    
    avg_rank <- pred_dt_lga_order[i, avg_rank]
    
    lga_rank_vt <- lga_rank_mt[lga_id, ]
    lga_rank_vt <- ifelse(lga_rank_vt <= 10, lga_rank_vt, 10)
    
    lga_ranktable <- as.data.table(table(lga_rank_vt))
    lga_ranktable <- merge(data.table(rank = as.character(1:10)), lga_ranktable, 
                             by.x = "rank", by.y = "lga_rank_vt", all.x = T)
    lga_ranktable[, "rank" := as.integer(rank)]
    lga_ranktable <- lga_ranktable[order(rank)]
    lga_ranktable[is.na(N), "N"] <- 0
    
    barplot(lga_ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(lga_name, ", ", state_name, ", ER = ", format(round(avg_rank, 1), nsmall = 1)),
            xaxt = "n", yaxt = "n", col = "#31a354", border = F)
    axis(1, at = 10:1-0.5, labels = c("10+", as.character(9:1)), tick = F)
  }
  dev.off()
  
  # bottom 5 LGA
  bitmap(paste0("plots_rank/", m, "_bottom5_lga.jpeg"),
         width = 2, height = 4, units = "in", type = "jpeg", res = 1000)
  par(mar = c(2.5, 1, 2, 1), mfrow = c(5, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_lga_order <- pred_dt_lga[order(-avg_rank)]
  
  for (i in 1:5){
    # i <- 1
    lga_id <- pred_dt_lga_order[i, LGAID]
    lga_name <- map_nga_lga_shp$LGAName[lga_id]
    
    state_id <- pred_dt_lga_order[i, StateID]
    state_name <- map_nga_lga_shp$StateName[lga_id]
    
    avg_rank <- pred_dt_lga_order[i, avg_rank]
    
    lga_rank_vt <- lga_rank_mt[lga_id, ]
    lga_rank_vt <- ifelse(lga_rank_vt >= (774-9), lga_rank_vt, (774-9))
    
    lga_ranktable <- as.data.table(table(lga_rank_vt))
    lga_ranktable <- merge(data.table(rank = as.character(774:(774-9))), lga_ranktable, 
                           by.x = "rank", by.y = "lga_rank_vt", all.x = T)
    lga_ranktable[, "rank" := as.integer(rank)]
    lga_ranktable <- lga_ranktable[order(rank)]
    lga_ranktable[is.na(N), "N"] <- 0
    
    barplot(lga_ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(lga_name, ", ", state_name, ", ER = ", format(round(avg_rank, 1), nsmall = 1)),
            xaxt = "n", yaxt = "n", col = "#fc4e2a", border = F)
    axis(1, at = 10:1-0.5, labels = c(as.character(774:(774-8)), "765-"), tick = F)
  }
  dev.off()
  
  # all states line plots
  bitmap(paste0("plots_rank/", m, "_all_lga_line.jpeg"),
         width = 30, height = 15, units = "in", type = "jpeg", res = 400)
  par(mar=c(4, 4, 1.5, 1), mfrow = c(1, 1))
  # par(mar = c(2.5, 1, 2, 1), mfrow = c(1, 1))
  
  pred_dt_lga_order <- merge(pred_dt_lga, 
                               map_nga_lga_shp@data)
  pred_dt_lga_order <- pred_dt_lga_order[order(-avg_rank)]
  pred_dt_lga_order[, "row_num"] <- 1:774
  
  plot(x = -100, y = -100, 
       xlim = c(0, 774), ylim = c(0, 774), 
       ann = FALSE, xaxt = "n", yaxt = "n")
  title(main = "Range of LGA Posterior Rankings", font.main = 2, cex.main = 0.9)
  mtext(side = 2, text = "rank", line = 2.7, font = 2, cex = 0.9)
  axis(1, at = 1:774, labels = pred_dt_lga_order$LGAName, 
       las = 2, cex.axis = 0.1, font = 1)
  axis(2, at = c(1, seq(50, 774, 50), 774), labels = rev(c(1, seq(50, 774, 50), 774)), 
       las = 1, cex.axis = 0.9, font = 1)
  
  points(x = pred_dt_lga_order[, row_num],
         y = 775-pred_dt_lga_order[, avg_rank], 
         col = "#08519c", pch = 16, cex = 0.6)
  
  segments(pred_dt_lga_order[, row_num], 
           775-pred_dt_lga_order[, low_rank],
           pred_dt_lga_order[, row_num], 
           775-pred_dt_lga_order[, up_rank],
           col = "#08519c", lwd = 0.65)
  dev.off()
  
  print(m)
}
