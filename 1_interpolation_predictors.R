# ================================================================
# Process and Interpolate Global Data for Predictions
# Example: phosphate Data
# Target: Interpolate depth data to a unified 1бу б┴ 1бу resolution
# Author: [Shihai Mao]
# Date: [2024-11-22]
# ================================================================
# Data Sources
# The following environmental predictors are used for the modeling process:
#- **Sea Surface Temperature (SST)**: 0.25бу resolution, sourced from the World Ocean Atlas 2018 (WOA18).
#- **Sea Surface Salinity (SSS)**: 0.25бу resolution, sourced from the World Ocean Atlas 2018 (WOA18).
#- **Phosphate Concentration ([PO43иC])**: 1бу resolution, sourced from the World Ocean Atlas 2018 (WOA18).
#- **Nitrate Concentration ([NO3иC])**: 1бу resolution, sourced from the World Ocean Atlas 2018 (WOA18).
#- **Seasonal mean surface seawater concentration of total ammonium ([NHx]), sourced from the previous publication (Paulot et al., 2020; https://doi.org/10.1029/2019MS002026).
#- **Dissolved Organic Phosphorus (DOP)**: 2бу resolution, sourced from previous publications (Liang et al., 2022; https://doi.org/10.1038/s41561-022-00988-1).
#- **Surface Chlorophyll-a (Chla)**: 0.083бу resolution, derived from NASA Aqua MODIS climatology.
#- **Particulate Organic Carbon (POC)**: 0.083бу resolution, derived from NASA Aqua MODIS climatology.
#- **Surface N* Values**: Computed as the difference between DIN and 16б┴DIP.
#- **[DIN] is the sum of surface seawater concentration of [NO3иC] and [NHx]. 
#- **Global total inorganic nitrogen deposition fluxes sourced from https://conservancy.umn.edu/items/1b5b51bd-5167-4407-aac3-232795284584

# ================================================================
# Part 1: Setup and Data Loading
# ================================================================
setwd("D:/R_language/R")

# Load libraries
library(sf)
library(sp)
library(raster)
library(gstat)
library(tidyverse)
library(maptools)
library(ncdf4)

# Load data
pivalue_sf <- na.omit(read.csv("D:/R language/R/phosphate.csv", header = T))

# Convert to spatial data frame
pivalue_sf <- st_as_sf(pivalue_sf, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
pivalue_sf <- st_transform(pivalue_sf, crs = 4326)

# Load boundary shapefile
bound <- st_read("ne_10m_ocean.shp")
bound_sf <- st_as_sf(bound)
bound_sf <- st_transform(bound_sf, CRS = 4326)
bound_sf <- st_make_valid(bound_sf)

# ================================================================
# Part 2: Grid Creation and Projection Setup
# ================================================================
# Convert boundary to Spatial object
bound_sp <- as_Spatial(bound_sf)

pivalue_sp <- pivalue_sf %>% as_Spatial()
pivalue_sp@bbox <- pivalue_sp@bbox

# Create a 1бу x 1бу grid
grd <- expand.grid(seq(-179.5, 179.5, 1), seq(-89.5, 89.5, 1))
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE
fullgrid(grd) <- TRUE

# Set the same projection for the grid
proj4string(pivalue_sp) <- proj4string(pivalue_sp)
proj4string(grd) <- proj4string(pivalue_sp)

# ================================================================
# Part 3: Interpolation and Output
# ================================================================
# Perform interpolation
vgm_exp <- variogram(Pi ~ 1, data = pivalue_sp)

vgm_model <- fit.variogram(vgm_exp, model = vgm("Gau"))

plot(vgm_exp, vgm_model, main = "Variogram Fit (Kriging)")

grd_pivalue <- gstat(id = "Pi", formula = Pi ~ 1, data = pivalue_sp, model = vgm_model, nmax = 10)

grd_pivalue_pred <- predict(grd_pivalue, newdata = grd)

# Convert prediction to raster
r <- raster(grd_pivalue_pred)

# Mask the raster to the ocean region
r.m <- mask(r, bound, inverse = FALSE)

# Convert to data frame and include coordinates
r.m_df <- as.data.frame(r.m)
r.m_df_xy <- cbind(r.m_df, coordinates(r.m))

# Save results
write.csv(r.m_df_xy, "D:/C/Desktop/Pi-1x1.csv")

# ================================================================
# Part 4: Interpolate Global Depth Data from NetCDF (.nc)
# ================================================================

# Open NetCDF file
ncfile <- nc_open("D:/R language/R/etopo2.nc")

# Extract coordinates and depth variable
lon <- ncvar_get(ncfile, "lon")
lat <- ncvar_get(ncfile, "lat")
depth <- ncvar_get(ncfile, "btdata")

# Create raster
r <- raster(t(depth), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +datum=WGS84"))
r <- flip(r, direction = "y")

# Extract non-NA points
r_points <- rasterToPoints(r, spatial = TRUE)
r_df <- as.data.frame(r_points)
colnames(r_df) <- c("depth", "lon", "lat")

# Convert to spatial points
coordinates(r_df) <- ~lon + lat
proj4string(r_df) <- CRS("+proj=longlat +datum=WGS84")

# Create target grid
grd <- expand.grid(seq(-179.5, 179.5, 1), seq(-89.5, 89.5, 1))
names(grd) <- c("lon", "lat")
coordinates(grd) <- ~lon + lat
gridded(grd) <- TRUE
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")

# Sample 5000 points for Variogram
r_df_sample <- r_df[sample(1:nrow(r_df), 5000), ]
coordinates(r_df_sample) <- ~lon + lat
proj4string(r_df_sample) <- CRS("+proj=longlat +datum=WGS84")

vgm_exp <- variogram(depth ~ 1, data = r_df_sample)
vgm_model <- fit.variogram(vgm_exp, model = vgm("Sph"))

plot(vgm_exp, vgm_model, main = "Variogram Model Fit")

krig_model <- gstat(id = "depth", formula = depth ~ 1, data = r_df, model = vgm_model, nmax = 10)
krig_result <- predict(krig_model, newdata = grd)

# Export result
r_krig <- raster(krig_result["depth.pred"])

r_masked <- mask(r_krig, bound_sp, inverse = FALSE)

masked_df <- as.data.frame(r_masked, xy = TRUE)

write.csv(masked_df, "D:/R language/R/depth_kriging_1x1.csv", row.names = FALSE)

nc_close(ncfile)
