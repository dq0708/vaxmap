library(foreign)
library(INLA)
library(Hmisc)
library(maptools)
library(raster)
library(xtable)
library(data.table)
library(broom)
library(rgdal)
library(spdep)
library(rgeos)
library(rasterVis)
library(survey)
library(sp)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(latticeExtra)
library(stringi)
library(Biograph)

##############################
# useful functions
##############################

logit<-function(x){log(x/(1-x))}

expit<-function(x){exp(x)/(1+exp(x))}

match_loc_spdfvar <- function(loc_x, loc_y, spdf, varname){
  loc_df <- data.frame(x = loc_x, y = loc_y)
  coordinates(loc_df) <- ~x+y
  proj4string(loc_df) <- proj4string(spdf)
  return(over(loc_df, spdf)[, varname])
}

match_loc_raster <- function(loc_x, loc_y, raster){
  loc_df <- data.frame(x = loc_x, y = loc_y)
  return(extract(x = raster, y = loc_df))
}

match_loc_raster_buffer <- function(loc_x, loc_y, raster, radius){
  loc_df <- data.frame(x = loc_x, y = loc_y)
  return(extract(x = raster, y = loc_df, buffer = radius, fun = mean, na.rm = T))
}

##################################
# Boundary files
##################################

map_nga_state_shp <- readOGR(dsn = "folder", 
                             layer = "boundary", stringsAsFactors = F)
map_nga_state_shp <- map_nga_state_shp[order(map_nga_state_shp$StateName), ]
map_nga_state_shp$StateID <- 1:length(map_nga_state_shp$StateName)

###

map_nga_lga_shp <- readOGR(dsn = "folder", 
                             layer = "boundary", stringsAsFactors = F)
map_nga_lga_shp <- merge(map_nga_lga_shp, 
                         data.table(StateCode = map_nga_state_shp$StateCode,
                                    StateID = map_nga_state_shp$StateID,
                                    StateName = map_nga_state_shp$StateName),
                         by = "StateCode")
map_nga_lga_shp <- map_nga_lga_shp[order(map_nga_lga_shp$StateName, map_nga_lga_shp$LGAName), ]
map_nga_lga_shp$LGAID <- 1:length(map_nga_lga_shp$LGAName)

###
map_nga_state_shp <- map_nga_state_shp[, c("StateID", "StateName")]
map_nga_lga_shp <- map_nga_lga_shp[, c("StateID", "StateName", "LGAID", "LGAName")]

############################################
# adjacency matrix
############################################

nb_state <- poly2nb(map_nga_state_shp, queen = F, row.names = map_nga_state_shp$StateID)
Amat_state <- nb2mat(nb_state, style = "B", zero.policy = TRUE)

nb_lga <- poly2nb(map_nga_lga_shp, queen = F, row.names = map_nga_lga_shp$LGAID)
Amat_lga <- nb2mat(nb_lga, style = "B", zero.policy = TRUE)

############################################
# create mesh using lga map
############################################

nga_border <- gUnaryUnion(map_nga_lga_shp, map_nga_lga_shp@data$LGAID)
nga_boundary <- inla.sp2segment(nga_border)
nga_mesh <- inla.mesh.2d(boundary = nga_boundary, cutoff = 0.08, max.edge = c(0.08, 0.8))
# plot(nga_mesh , asp = 1, main = "")

##################################
# Nigeria 2018 DHS data
##################################

# survey data
dhs2018_chil_raw <- read.dta("folder/data.dta", 
                             convert.factors = F)

col_vt <- c("caseid", "v001", "v002", "v003", "v005", "v006", "v007", "v008", "v016", 
            "v021", "v022", "v023", "v024", "v025", 
            "h1", "h2", "h2d", "h2m", "h2y", "h3", "h3d", "h3m", "h3y", 
            "h4", "h4d", "h4m", "h4y", "h5", "h5d", "h5m", "h5y", 
            "h6", "h6d", "h6m", "h6y", "h7", "h7d", "h7m", "h7y", 
            "h8", "h8d", "h8m", "h8y", "h9", "h9d", "h9m", "h9y", 
            "h0", "h0d", "h0m", "h0y", "h10",
            "b5", "b8", "hw1", "hw16", "sstate")

col_name_vt <- c("caseid", "cluster", "hh", "respond", "hhsampwgt", "interv_m", "interv_y", "interv_date", "interv_d", 
                 "psu", "strata", "domain", "region", "type", 
                 "hascard", "bcg", "bcg_d", "bcg_m", "bcg_y", "dpt1", "dpt1_d", "dpt1_m", "dpt1_y",
                 "polio1", "polio1_d", "polio1_m", "polio1_y", "dpt2", "dpt2_d", "dpt2_m", "dpt2_y",
                 "polio2", "polio2_d", "polio2_m", "polio2_y", "dpt3", "dpt3_d", "dpt3_m", "dpt3_y",
                 "polio3", "polio3_d", "polio3_m", "polio3_y", "measles", "measles_d", "measles_m", "measles_y",
                 "polio0", "polio0_d", "polio0_m", "polio0_y", "evervacc",
                 "alive", "age_y", "age_m", "birth_d", "state")

dhs2018_chil_dt <- as.data.table(dhs2018_chil_raw[, col_vt])
setnames(dhs2018_chil_dt, col_name_vt)

state_lookup <- data.table(state = seq(10, 370, 10),
                           statename = c("SOKOTO", "ZAMFARA", "KATSINA", "JIGAWA", "YOBE", "BORNO", "ADAMAWA", 
                                         "GOMBE", "BAUCHI", "KANO", "KADUNA", "KEBBI", "NIGER", "FCT-ABUJA", "NASARAWA", 
                                         "PLATEAU", "TARABA", "BENUE", "KOGI", "KWARA", "OYO", "OSUN", "EKITI", "ONDO", 
                                         "EDO", "ANAMBRA", "ENUGU", "EBONYI", "CROSS RIVER", "AKWA IBOM", "ABIA", "IMO", 
                                         "RIVERS", "BAYELSA", "DELTA", "LAGOS", "OGUN"))
state_lookup <- state_lookup[order(statename)]
state_lookup$StateID <- sort(map_nga_state_shp$StateID)
state_lookup$StateName <- sort(map_nga_state_shp$StateName)
dhs2018_chil_dt <- merge(dhs2018_chil_dt, state_lookup[, .(state, StateID, StateName)], by = "state", all.x = T)

dhs2018_chil_dt[, "sampwgt" := hhsampwgt/1e6] 
dhs2018_chil_dt[, "IsUrban" := ifelse(type == 1, 1, 0)]
dhs2018_chil_dt[, "MCV" := ifelse(is.na(measles), 0, 
                                  ifelse(measles %in% c(1, 2, 3), 1, 0))]

# GPS data
dhs2018_shp <- readOGR("folder",
                       layer = "cluster_locations")

dhs2018_chil_dt <- merge(dhs2018_chil_dt, 
                         data.table(cluster = dhs2018_shp$DHSCLUST,
                                    long = dhs2018_shp$LONGNUM,
                                    lat = dhs2018_shp$LATNUM),
                         by = "cluster", all.x = T)

# match LGA
dhs2018_chil_dt[, "LGAID"] <- match_loc_spdfvar(loc_x = dhs2018_chil_dt$long,
                                                 loc_y = dhs2018_chil_dt$lat,
                                                 spdf = map_nga_lga_shp,
                                                 varname = "LGAID")
dhs2018_chil_dt[, "LGAName"] <- match_loc_spdfvar(loc_x = dhs2018_chil_dt$long,
                                                 loc_y = dhs2018_chil_dt$lat,
                                                 spdf = map_nga_lga_shp,
                                                 varname = "LGAName")

svy_data_dt <- dhs2018_chil_dt[long > 0 & lat > 0 & alive == 1 & age_y == 1, 
                               c("StateID", "StateName", "LGAID", "LGAName",
                                 "IsUrban", "cluster", "hh", "sampwgt", 
                                 "MCV", "long", "lat")]

LGAID_miss <- order(gDistance(SpatialPoints(coords = data.frame(c(3.405828), c(6.56706))), 
                              map_nga_lga_shp, byid = T))[1]
LGAName_miss <- map_nga_lga_shp$LGAName[LGAID_miss]
svy_data_dt[is.na(LGAID), "LGAID"] <- LGAID_miss
svy_data_dt[is.na(LGAName), "LGAName"] <- LGAName_miss

svy_data_dt <- svy_data_dt[order(StateID, LGAID, long, lat)]

##########

cluster_data_dt <- svy_data_dt[, .(y = sum(MCV), n = length(MCV)),
                               by = .(StateID, StateName, LGAID, LGAName, cluster, long, lat, IsUrban)]

cluster_data_dt[, "i_id" := StateID]
cluster_data_dt[, "j_id" := LGAID]
cluster_data_dt <- cluster_data_dt[order(StateID, LGAID, long, lat)]
cluster_data_dt[, "c_id"] <- 1:nrow(cluster_data_dt)

##################################
# create prediction grid
##################################

# create grid for predictions
border_tidy <- tidy(nga_border)
x_min <- min(border_tidy$long)
x_max <- max(border_tidy$long)
y_min <- min(border_tidy$lat)
y_max <- max(border_tidy$lat)

xy_min <- floor(min(c(x_min, y_min)))
xy_max <- ceiling(max(c(x_max, y_max)))

limit_points <- SpatialPoints(cbind(x = c(xy_min, xy_min, xy_max, xy_max), 
                                    y = c(xy_min, xy_max, xy_min, xy_max)), 
                              proj4string =  CRS("+init=epsg:4326"))

cell_size <- 0.00833333
x_n_cell <- ceiling((coordinates(limit_points)[3, 1] - coordinates(limit_points)[2, 1])/cell_size)
y_n_cell <- ceiling((coordinates(limit_points)[2, 2] - coordinates(limit_points)[3, 2])/cell_size)
nxy <- c(x_n_cell, y_n_cell)

pred_grid <- inla.mesh.projector(nga_mesh, xlim = c(xy_min, xy_max),  ylim = c(xy_min, xy_max), dims = nxy)

pred_dt_grid <- data.table(long = pred_grid$lattice$loc[, 1],
                           lat = pred_grid$lattice$loc[, 2])

pred_dt_grid <- cbind(pred_dt_grid, 
                      match_loc_spdfvar(loc_x = pred_dt_grid$long,
                                        loc_y = pred_dt_grid$lat,
                                        spdf = map_nga_lga_shp,
                                        varname = c("StateID", "StateName", "LGAID", "LGAName")))
pred_dt <- pred_dt_grid[!is.na(LGAID), ]

pred_dt[, "i_id" := StateID]
pred_dt[, "j_id" := LGAID]
pred_dt[, "c_id"] <- 0

##################################
# read in aridity layer
##################################

aridity_dat <- raster("folder/data")
xylim <- extent(x_min-0.2, x_max+0.2, y_min-0.2, y_max+0.2)
aridity_raster_1km <- crop(aridity_dat, xylim)

pred_dt[, "aridity"] <- match_loc_raster(loc_x = pred_dt$long,
                                         loc_y = pred_dt$lat,
                                         raster = aridity_raster_1km)
pred_dt[is.na(aridity), "aridity"] <- median(pred_dt[!is.na(aridity), aridity])

cluster_data_dt[, "aridity"] <- match_loc_raster_buffer(loc_x = cluster_data_dt$long,
                                                        loc_y = cluster_data_dt$lat,
                                                        raster = aridity_raster_1km, 
                                                        radius = ifelse(cluster_data_dt$IsUrban== 1, 2000, 5000))

##################################
# read in poverty layer
##################################

poverty_dat <- raster("folder/data")
xylim <- extent(x_min-0.2, x_max+0.2, y_min-0.2, y_max+0.2)
poverty_raster_1km <- crop(poverty_dat, xylim)

pred_dt[, "poverty"] <- match_loc_raster(loc_x = pred_dt$long,
                                         loc_y = pred_dt$lat,
                                         raster = poverty_raster_1km)
pred_dt[is.na(poverty), "poverty"] <- median(pred_dt[!is.na(poverty), poverty])

cluster_data_dt[, "poverty"] <- match_loc_raster_buffer(loc_x = cluster_data_dt$long,
                                                        loc_y = cluster_data_dt$lat,
                                                        raster = poverty_raster_1km, 
                                                        radius = ifelse(cluster_data_dt$IsUrban== 1, 2000, 5000))

##################################
# read in EVI layer
##################################

evi_dat <- raster("folder/data")
xylim <- extent(x_min-0.2, x_max+0.2, y_min-0.2, y_max+0.2)
evi_raster_1km <- crop(evi_dat, xylim)

pred_dt[, "evi"] <- match_loc_raster(loc_x = pred_dt$long,
                                     loc_y = pred_dt$lat,
                                     raster = evi_raster_1km)

pred_dt[is.na(evi), "evi"] <- median(pred_dt[!is.na(evi), evi])

cluster_data_dt[, "evi"] <- match_loc_raster_buffer(loc_x = cluster_data_dt$long,
                                                    loc_y = cluster_data_dt$lat,
                                                    raster = evi_raster_1km, 
                                                    radius = ifelse(cluster_data_dt$IsUrban== 1, 2000, 5000))

##################################
# read in travel time layer
##################################

travel_time_dat <- raster("folder/data")
xylim <- extent(x_min-0.2, x_max+0.2, y_min-0.2, y_max+0.2)
travel_time_raster_1km <- crop(travel_time_dat, xylim)

pred_dt[, "travel_time"] <- match_loc_raster(loc_x = pred_dt$long,
                                             loc_y = pred_dt$lat,
                                             raster = travel_time_raster_1km)

cluster_data_dt[, "travel_time"] <- match_loc_raster_buffer(loc_x = cluster_data_dt$long,
                                                            loc_y = cluster_data_dt$lat,
                                                            raster = travel_time_raster_1km, 
                                                            radius = ifelse(cluster_data_dt$IsUrban== 1, 2000, 5000))
pred_dt[, "log_travel_time" := ifelse(is.na(travel_time), NA, 
                                      ifelse(travel_time == 0, log(1), log(travel_time)))]
pred_dt[is.na(log_travel_time), "log_travel_time"] <- median(pred_dt[!is.na(log_travel_time), log_travel_time])

cluster_data_dt[, "log_travel_time" := ifelse(travel_time == 0, 0, log(travel_time))]

##################################
# read in night light layer
##################################

night_dat <- raster("folder/data")
xylim <- extent(x_min-0.2, x_max+0.2, y_min-0.2, y_max+0.2)
night_raster_1km <- crop(night_dat, xylim)

pred_dt[, "night"] <- match_loc_raster(loc_x = pred_dt$long,
                                       loc_y = pred_dt$lat,
                                       raster = night_raster_1km)
pred_dt[, "night" := night + 0.5]

cluster_data_dt[, "night"] <- match_loc_raster_buffer(loc_x = cluster_data_dt$long,
                                                      loc_y = cluster_data_dt$lat,
                                                      raster = night_raster_1km, 
                                                      radius = ifelse(cluster_data_dt$IsUrban== 1, 2000, 5000))
cluster_data_dt[, "night":= night + 0.5] 
pred_dt[, "log_night" := ifelse(night <= 0, log(0.01), log(night))]

cluster_data_dt[, "log_night" := ifelse(night <= 0, log(0.01), log(night))]

##################################
# population data - WP, 2018, all age
##################################

wp_pop_2018_raw  <- raster("folder/data")

wp_pop_2018_raster_1km <- raster::aggregate(wp_pop_2018_raw, fact = 10, fun = sum, na.rm = T)

pred_dt[, "pop_all_wp"] <- match_loc_raster(loc_x = pred_dt$long,
                                            loc_y = pred_dt$lat,
                                            raster = wp_pop_2018_raster_1km)

pred_dt[, "pop_all_wp" := ifelse(is.na(pop_all_wp), 0, pop_all_wp)]

##################################
# population data - WP, 2014, under 5
##################################

wp_pop_2014_raw <- raster("folder/data")
wp_pop_2014_raster_1km <- raster::aggregate(wp_pop_2014_raw, fact = 10, fun = sum, na.rm = T)
pred_dt[, "pop_u5_wp"] <- match_loc_raster(loc_x = pred_dt$long,
                                           loc_y = pred_dt$lat,
                                           raster = wp_pop_2014_raster_1km)
pred_dt[, "pop_u5_wp" := ifelse(is.na(pop_u5_wp), 0, pop_u5_wp)]

##################################
# create urban-rural layer based on pred grid and DHS 
##################################

dhs_2018_pop_dt <- fread("folder/data")
dhs_2018_pop_dt <- dhs_2018_pop_dt[order(state)]
dhs_2018_pop_dt[, "StateName"] <- map_nga_state_shp$StateName
dhs_2018_pop_dt[, "StateID"] <- map_nga_state_shp$StateID
dhs_2018_pop_dt[, "urban_frac" := pop_urban/(pop_urban+pop_rural)]

####

pred_dt_temp <- NULL

for(i in 1:37){
  # i <- 1
  pred_dt_i <- pred_dt[StateID == i, ]
  
  # using worldpop pop
  pred_dt_i <- pred_dt_i[order(-pop_all_wp)]
  pred_dt_i[, "cum_pop_frac_wp" := cumsum(pop_all_wp)/sum(pop_all_wp)]
  pred_dt_i[, "IsUrban_wp" := ifelse(cum_pop_frac_wp <= dhs_2018_pop_dt[StateID == i, urban_frac]*1.01, 1, 0)]
  
  pred_dt_temp <- rbind(pred_dt_temp, pred_dt_i)
  
  print(i)
}

pred_dt_grid <- pred_dt_temp[, c("long", "lat", "StateID", "StateName", "LGAID", "LGAName",
                                 "i_id", "j_id", "c_id", "aridity", "poverty", "evi", 
                                 "log_travel_time", "log_night", "pop_all_wp", "pop_u5_wp",
                                 "IsUrban_wp")]

##################################
# create urban-rural pop and EA summary by state and lga based on pred grid and DHS 
##################################

dhs_2018_ea_dt <- fread("folder/data")
dhs_2018_ea_dt <- dhs_2018_ea_dt[order(state)]
dhs_2018_ea_dt[, "StateName"] <- map_nga_state_shp$StateName
dhs_2018_ea_dt[, "StateID"] <- map_nga_state_shp$StateID

####
# using worldpop pop
####

pop_ea_state_dt_wp <- data.table(StateID = map_nga_state_shp$StateID,
                                 StateName = map_nga_state_shp$StateName)

pop_ea_lga_dt_wp <- data.table(StateID = map_nga_lga_shp$StateID,
                               StateName = map_nga_lga_shp$StateName,
                               LGAID = map_nga_lga_shp$LGAID,
                               LGAName = map_nga_lga_shp$LGAName)

for(i in 1:37){
  # i <- 1
  pred_dt_i <- pred_dt_grid[StateID == i, ]
  
  total_pop_i <- sum(pred_dt_i[, pop_all_wp])
  urban_pop_i <- sum(pred_dt_i[IsUrban_wp == 1, pop_all_wp])
  rural_pop_i <- total_pop_i - urban_pop_i
  
  pop_ea_state_dt_wp[StateID == i, "pop_all"] <- total_pop_i
  pop_ea_state_dt_wp[StateID == i, "urban_pop_all"] <- urban_pop_i
  pop_ea_state_dt_wp[StateID == i, "rural_pop_all"] <- rural_pop_i
  
  pop_ea_state_dt_wp[StateID == i, "urban_popfrac_all"] <- urban_pop_i/total_pop_i
  pop_ea_state_dt_wp[StateID == i, "rural_popfrac_all"] <- 1 - urban_pop_i/total_pop_i
  
  ###
  
  total_u5pop_i <- sum(pred_dt_i[, pop_u5_wp])
  urban_u5pop_i <- sum(pred_dt_i[IsUrban_wp == 1, pop_u5_wp])
  rural_u5pop_i <- total_u5pop_i - urban_u5pop_i
  
  pop_ea_state_dt_wp[StateID == i, "pop_u5"] <- total_u5pop_i
  pop_ea_state_dt_wp[StateID == i, "urban_pop_u5"] <- urban_u5pop_i
  pop_ea_state_dt_wp[StateID == i, "rural_pop_u5"] <- rural_u5pop_i
  
  pop_ea_state_dt_wp[StateID == i, "urban_popfrac_u5"] <- urban_u5pop_i/total_u5pop_i
  pop_ea_state_dt_wp[StateID == i, "rural_popfrac_u5"] <- 1 - urban_u5pop_i/total_u5pop_i
  
  ###
  
  total_ea_i <- dhs_2018_ea_dt[StateID == i, ea_total]
  urban_ea_i <- dhs_2018_ea_dt[StateID == i, ea_urban]
  rural_ea_i <- total_ea_i - urban_ea_i
  
  pop_ea_state_dt_wp[StateID == i, "ea"] <- total_ea_i
  pop_ea_state_dt_wp[StateID == i, "urban_ea"] <- urban_ea_i
  pop_ea_state_dt_wp[StateID == i, "rural_ea"] <- rural_ea_i
  
  pop_ea_state_dt_wp[StateID == i, "urban_eafrac"] <- urban_ea_i/total_ea_i
  pop_ea_state_dt_wp[StateID == i, "rural_eafrac"] <- 1 - urban_ea_i/total_ea_i
  
  ###
  
  urban_avgeapop_i <- urban_pop_i/urban_ea_i
  rural_avgeapop_i <- rural_pop_i/rural_ea_i
  
  ###
  
  for (j in pop_ea_lga_dt_wp[StateID == i, LGAID]){
    # j <- 1
    pred_dt_j <- pred_dt_grid[LGAID == j, ]
    
    total_pop_j <- sum(pred_dt_j[, pop_all_wp])
    urban_pop_j <- sum(pred_dt_j[IsUrban_wp == 1, pop_all_wp])
    rural_pop_j <- total_pop_j - urban_pop_j
    
    pop_ea_lga_dt_wp[LGAID == j, "pop_all"] <- total_pop_j
    pop_ea_lga_dt_wp[LGAID == j, "urban_pop_all"] <- urban_pop_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_pop_all"] <- rural_pop_j
    
    pop_ea_lga_dt_wp[LGAID == j, "urban_popfrac_all"] <- urban_pop_j/total_pop_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_popfrac_all"] <- 1 - urban_pop_j/total_pop_j
    
    ###
    
    total_u5pop_j <- sum(pred_dt_j[, pop_u5_wp])
    urban_u5pop_j <- sum(pred_dt_j[IsUrban_wp == 1, pop_u5_wp])
    rural_u5pop_j <- total_u5pop_j - urban_u5pop_j
    
    pop_ea_lga_dt_wp[LGAID == j, "pop_u5"] <- total_u5pop_j
    pop_ea_lga_dt_wp[LGAID == j, "urban_pop_u5"] <- urban_u5pop_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_pop_u5"] <- rural_u5pop_j
    
    pop_ea_lga_dt_wp[LGAID == j, "urban_popfrac_u5"] <- urban_u5pop_j/total_u5pop_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_popfrac_u5"] <- 1 - urban_u5pop_j/total_u5pop_j
    
    ###
    
    urban_ea_j <- round(urban_pop_j/urban_avgeapop_i)
    rural_ea_j <- round(rural_pop_j/rural_avgeapop_i)
    total_ea_j <- urban_ea_j + rural_ea_j
    
    pop_ea_lga_dt_wp[LGAID == j, "ea"] <- total_ea_j
    pop_ea_lga_dt_wp[LGAID == j, "urban_ea"] <- urban_ea_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_ea"] <- rural_ea_j
    
    pop_ea_lga_dt_wp[LGAID == j, "urban_eafrac"] <- urban_ea_j/total_ea_j
    pop_ea_lga_dt_wp[LGAID == j, "rural_eafrac"] <- 1 - urban_ea_j/total_ea_j

    
    print(j)
    
  }
}

##################################
# simulate prediction clusters
##################################

pred_dt_wp <- NULL

set.seed(12345)

for(j in 1:774){
  # j <- 1
  
  # urban wp
  pred_dt_j <- pred_dt_grid[LGAID == j & IsUrban_wp == 1, ]
  N_pix_j <- nrow(pred_dt_j)
  N_ea_j <- pop_ea_lga_dt_wp[LGAID == j, urban_ea]
  if (N_pix_j > 0 & N_ea_j > 0){
    idx_samp <- sample(x = N_pix_j, size = N_ea_j, replace = T, 
                       prob = pred_dt_j$pop_all_wp/sum(pred_dt_j$pop_all_wp))
    pred_dt_wp_j <- pred_dt_j[idx_samp, c("long", "lat", "StateID", "StateName", "LGAID", "LGAName", "IsUrban_wp",
                                          "i_id", "j_id", "c_id", "aridity", "poverty", "evi", "log_travel_time", "log_night")]
    setnames(pred_dt_wp_j, c("long", "lat", "StateID", "StateName", "LGAID", "LGAName", "IsUrban",
                             "i_id", "j_id", "c_id", "aridity", "poverty", "evi", "log_travel_time", "log_night"))
    pred_dt_wp_j[, "pop_all"] <- pop_ea_lga_dt_wp[LGAID == j, urban_pop_all]/N_ea_j
    pred_dt_wp_j[, "pop_u5"] <- pop_ea_lga_dt_wp[LGAID == j, urban_pop_u5]/N_ea_j
    pred_dt_wp <- rbind(pred_dt_wp, pred_dt_wp_j)
  }
  
  # rural wp
  pred_dt_j <- pred_dt_grid[LGAID == j & IsUrban_wp == 0, ]
  N_pix_j <- nrow(pred_dt_j)
  N_ea_j <- pop_ea_lga_dt_wp[LGAID == j, rural_ea]
  if (N_pix_j > 0 & N_ea_j > 0){
    idx_samp <- sample(x = N_pix_j, size = N_ea_j, replace = T, 
                       prob = pred_dt_j$pop_all_wp/sum(pred_dt_j$pop_all_wp))
    pred_dt_wp_j <- pred_dt_j[idx_samp, c("long", "lat", "StateID", "StateName", "LGAID", "LGAName", "IsUrban_wp",
                                          "i_id", "j_id", "c_id", "aridity", "poverty", "evi", "log_travel_time", "log_night")]
    setnames(pred_dt_wp_j, c("long", "lat", "StateID", "StateName", "LGAID", "LGAName", "IsUrban",
                             "i_id", "j_id", "c_id", "aridity", "poverty", "evi", "log_travel_time", "log_night"))
    pred_dt_wp_j[, "pop_all"] <- pop_ea_lga_dt_wp[LGAID == j, rural_pop_all]/N_ea_j
    pred_dt_wp_j[, "pop_u5"] <- pop_ea_lga_dt_wp[LGAID == j, rural_pop_u5]/N_ea_j
    pred_dt_wp <- rbind(pred_dt_wp, pred_dt_wp_j)
  }
  
  print(j)
    
}

#######

describe(cluster_data_dt$aridity)
describe(cluster_data_dt$poverty)
describe(cluster_data_dt$evi)
describe(cluster_data_dt$log_travel_time)
describe(cluster_data_dt$log_night)

cluster_data_dt$aridity10k <- cluster_data_dt$aridity/10000
cluster_data_dt$evi1k <- cluster_data_dt$evi/1000

pred_dt_grid$aridity10k <- pred_dt_grid$aridity/10000
pred_dt_grid$evi1k <- pred_dt_grid$evi/1000

pred_dt_wp$aridity10k <- pred_dt_wp$aridity/10000
pred_dt_wp$evi1k <- pred_dt_wp$evi/1000


#######

save(expit, logit, match_loc_raster, match_loc_spdfvar, match_loc_raster_buffer,
     map_nga_state_shp, map_nga_lga_shp, nb_state, Amat_state, nb_lga, Amat_lga,
     nga_border, nga_boundary, nga_mesh, dhs_2018_ea_dt, dhs_2018_pop_dt,
     svy_data_dt, cluster_data_dt, pred_dt_grid, pred_dt_wp, 
     pop_ea_state_dt_wp, pop_ea_lga_dt_wp, 
     file = "data_cleaned_1k.RData")
