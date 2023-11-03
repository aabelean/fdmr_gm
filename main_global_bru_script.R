###
# load packages
###
library(rgdal)
library(sp)
library(GEOmap)
library(INLA)
library(inlabru)
library(ncdf4)
library(data.table)

###
# Set paths
###
mesh_path <- NULL
obs_path <- NULL

###
# Load datasets
###
# Load in GRACE mascons
load(base::file.path(obs_path, "EWH_JPL_grace_global_GMGIA_to2015.RData"))
grace = base::data.frame(grace)

# Load in Greenland alt
load(base::file.path(obs_path,"Greenland_alt_annual_diff.RData"))
gris = base::data.frame(gris)

# Load in Antarctica alt
load(base::file.path(obs_path, "Antarctica_alt_annual_diff.RData"))
ais = base::data.frame(ais)

# Load in GIC
load(base::file.path(obs_path, "Hugonnet_GIC_annual_diff.RData"))
gic = base::data.frame(gic)
# replace NAs
gic$SERate[base::is.na(gic$SERate)] <- base::mean(gic$SERate,na.rm=TRUE)
gic$SERate[(gic$SERate==0)] <- base::mean(gic$SERate,na.rm=TRUE)

# Load in ocean altimetry
load(base::file.path(obs_path, "Altimetry_filtered_GMGIA_arctic.RData"))
altimetry = base::data.frame(altimetry)

# Load in Argo annual rates
load(base::file.path(obs_path, "Argo_flagged_annual_diff.RData"))
argo = base::data.frame(argo)

###
# Set number of year slices
###
n_years = base::as.integer(base::length(base::unique(grace$Year))) #10 in this case

###
# Set a function to create polygon blocks (resampling taking into account the disassociation distance/range)
###
poly_block <- function(i, df_sp, dis = 100000) {
  sp_i <- sp::SpatialPolygons(list(df_sp@polygons[[i]]), proj4string = sp::CRS("+proj=longlat"))
  area_i <- df_sp$area[i]
  n <- base::round(area_i / dis^2)

  if (n < 2) {
    grid_i <- sp::SpatialPoints(sp_i)
    ngrid_i <- 1
  } else {
    grid_i <- sp::spsample(sp_i, n = base::round(area_i / dis^2), type = "regular", offset = c(0.5, 0.5)) #need to convert to sf?
    ngrid_i <- base::length(grid_i)
  }

  grid_xyz <- base::do.call(cbind, GEOmap::Lll2xyz(lat = grid_i@coords[, 2], lon = grid_i@coords[, 1]))
  block_i <- base::rep(i, ngrid_i)
  weights <- base::rep(area_i / ngrid_i, ngrid_i)

  return(list(grid_xyz = grid_xyz, block = block_i, weights = weights, ngrid = ngrid_i))
}

###
# Prepare GRACE polygons (assuming information resolution of GRACE to be one mascon cell)
###
df_sp <- grace_sp
dis_n = 300000

grace_polygon_block <- base::lapply(1:base::nrow(grace_sp), poly_block, df_sp = df_sp, dis = dis_n)
grace_xyz <- base::do.call(rbind, base::lapply(grace_polygon_block, "[[", "grid_xyz"))
grace_block <- base::do.call(c, base::lapply(grace_polygon_block, "[[", "block"))
grace_weights <- base::do.call(c, base::lapply(grace_polygon_block, "[[", "weights"))
message(base::paste0('Number of new points (by blocks): ', base::length(grace_block)))

grace_points_dat<-base::cbind.data.frame(grace_xyz,grace[grace_block,],weight=grace_weights)
grace_points_dat$hydrology.group<-grace_points_dat$time
grace_points_dat$glaciers.group<-grace_points_dat$time
grace_points_dat$omass.group<-grace_points_dat$time
grace_points_dat$steric.group<-NA
grace_points_dat$icesheet.group<-grace_points_dat$time

###
# Similarly apply polygon blocks to GIC. Disassociation distance is expected to be ~10 km
###
message(c('Number of Hugonnet tile observations (polygons * timesteps): ', base::nrow(Hugonnet_sp)))
df_sp <- Hugonnet_sp

Hugonnet_polygon_block <- base::lapply(1:base::nrow(Hugonnet_sp), poly_block, df_sp = df_sp, dis = 25000)
Hugonnet_xyz <- base::do.call(rbind, base::lapply(Hugonnet_polygon_block, "[[", "grid_xyz"))
Hugonnet_block <- base::do.call(c, base::lapply(Hugonnet_polygon_block, "[[", "block"))
Hugonnet_weights <- base::do.call(c, base::lapply(Hugonnet_polygon_block, "[[", "weights"))
message(base::paste0('Number of new points (by blocks): ', base::length(Hugonnet_block)))

gic_points_dat<-base::cbind.data.frame(Hugonnet_xyz,gic[Hugonnet_block,],weight=Hugonnet_weights)
gic_points_dat$hydrology.group<-NA
gic_points_dat$glaciers.group<-gic_points_dat$time
gic_points_dat$omass.group<-NA
gic_points_dat$steric.group<-NA
gic_points_dat$icesheet.group<-NA

###
# Altimetry and steric data - just convert to Cartesian coords
###

gris_xyz<-base::do.call(cbind, GEOmap::Lll2xyz(lat = gris$Lat,
                                               lon = gris$Long))
gris$x<-gris_xyz[,1]
gris$y<-gris_xyz[,2]
gris$z<-gris_xyz[,3]
gris$hydrology.group<-NA
gris$glaciers.group<-NA
gris$omass.group<-NA
gris$steric.group<-NA
gris$icesheet.group<-gris$time

ais_xyz<-base::do.call(cbind, GEOmap::Lll2xyz(lat = ais$Lat,
                                              lon = ais$Long))
ais$x<-ais_xyz[,1]
ais$y<-ais_xyz[,2]
ais$z<-ais_xyz[,3]
ais$hydrology.group<-NA
ais$glaciers.group<-NA
ais$omass.group<-NA
ais$steric.group<-NA
ais$icesheet.group<-ais$time

altimetry2<- altimetry[!base::is.na(altimetry$Lat),]
altimetry_xyz<-base::do.call(cbind, GEOmap::Lll2xyz(lat = altimetry2$Lat,
                                                    lon = altimetry2$Long))
altimetry2$x<-altimetry_xyz[,1]
altimetry2$y<-altimetry_xyz[,2]
altimetry2$z<-altimetry_xyz[,3]
altimetry2$hydrology.group<-NA
altimetry2$glaciers.group<-NA
altimetry2$omass.group<-altimetry2$time
altimetry2$steric.group<-altimetry2$time
altimetry2$icesheet.group<-NA

argo2<- argo[!base::is.na(argo$Lat),]
argo_xyz<-base::do.call(cbind, GEOmap::Lll2xyz(lat = argo2$Lat,
                                               lon = argo2$Long))
argo2$x<-argo_xyz[,1]
argo2$y<-argo_xyz[,2]
argo2$z<-argo_xyz[,3]
argo2$hydrology.group<-NA
argo2$glaciers.group<-NA
argo2$omass.group<-NA
argo2$steric.group<-argo2$time
argo2$icesheet.group<-NA

###
# Bind together
###
comb_dat<-base::rbind(grace_points_dat[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                                          "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")],
                      gic_points_dat[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                                        "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")],
                      gris[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                              "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")],
                      ais[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                             "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")],
                      altimetry2[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                                    "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")],
                      argo2[,c("x","y","z","Rate","SERate","Year","time","Long","Lat",
                               "hydrology.group","glaciers.group","omass.group","steric.group","icesheet.group")])
sp::coordinates(comb_dat)<- c("x","y","z")

###
# Define groups
###
hydrology.group<-comb_dat$hydrology.group
glaciers.group<-comb_dat$glaciers.group
omass.group<-comb_dat$omass.group
steric.group<-comb_dat$steric.group
icesheet.group<-comb_dat$icesheet.group

###
# Set SPDEs
###
prior_rangeH <- 20
prior_rangeG <- 10
prior_rangeI <- 15
prior_rangeS <- 25
prior_rangeO <- 65

prior_sigmaH <- 100
prior_sigmaG <- 100
prior_sigmaI <- 100
prior_sigmaS <- 30
prior_sigmaO <- 10

# setting NA in pcmatern P fixes the hyperpar

SPDE_Hyd <- INLA::inla.spde2.pcmatern(mesh_land,
                                      prior.range=c(prior_rangeH * (pi/180),NA),prior.sigma=c(prior_sigmaH,NA))
SPDE_Glaciers <- INLA::inla.spde2.pcmatern(mesh_GIC,
                                           prior.range=c(prior_rangeG * (pi/180),NA),prior.sigma=c(prior_sigmaG,NA))
SPDE_IceSheet <- INLA::inla.spde2.pcmatern(mesh_IceSheets,
                                           prior.range=c(prior_rangeI * (pi/180),NA),prior.sigma=c(prior_sigmaI,NA))
SPDE_Steric <- INLA::inla.spde2.pcmatern(mesh_ocean,
                                         prior.range=c(prior_rangeS * (pi/180),NA),prior.sigma=c(prior_sigmaS,NA))
SPDE_omass <- INLA::inla.spde2.pcmatern(mesh_ocean,
                                        prior.range=c(prior_rangeO * (pi/180),NA),prior.sigma=c(prior_sigmaO,NA))
###
# Set fixed observation errors
###
hyper <- list(prec = list(fixed = TRUE, initial = 0))
prec_scale <- c(1/(grace_points_dat$SERate^2),1/(gic_points_dat$SERate^2),1/(gris$SERate^2),
                1/(ais$SERate^2),1/(altimetry2$SERate^2),1/(argo2$SERate^2))

rhoprior <- list(theta = list(prior = 'pccor1',
                              param = c(0, 0.9)))
###
# Define the model
###
cmp<- ~ -1 + hyd_field(main=coordinates,
                       model=SPDE_Hyd,
                       group=hydrology.group,
                       ngroup=n_years,
                       control.group=list(model="ar1",hyper=rhoprior))+
  GIC_field(main=coordinates,
            model=SPDE_Glaciers,
            group=glaciers.group,
            ngroup=n_years,
            control.group=list(model="ar1",hyper=rhoprior))+
  ice_field(main=coordinates,
            model=SPDE_IceSheet,
            group=icesheet.group,
            ngroup=n_years,
            control.group=list(model="ar1",hyper=rhoprior))+
  steric_field(main=coordinates,
               model=SPDE_Steric,
               group=steric.group,
               ngroup=n_years,
               control.group=list(model="ar1",hyper=rhoprior))+
  omass_field(main=coordinates,
              model=SPDE_omass,
              group=omass.group,
              ngroup=n_years,
              control.group=list(model="ar1",hyper=rhoprior))

model_formula<- Rate ~ .

###
# Run the model
###
res_bru<- inlabru::bru(cmp,
                       inlabru::like(
                         formula = model_formula,
                         family = "gaussian",
                         control.family = list(hyper = hyper),
                         data = comb_dat
                       ),
                       options = list(
                         scale=prec_scale,
                         inla.mode = 'classic',
                         verbose = TRUE))

###
# Analyse the output
###

## Give lat-lon coords
ocean_coords <- base::data.frame(GEOmap::Lxyz2ll(list(x = mesh_ocean$loc[,1], y = mesh_ocean$loc[,2], z = mesh_ocean$loc[,3])))
ocean_coords$index = base::seq(1:base::length(ocean_coords$lat))
ocean_coords_rep = base::data.frame(base::sapply(ocean_coords,base::rep.int,time=n_years))
ocean_coords_rep$time = base::rep(1:n_years, each = base::length(ocean_coords$lon))
ocean_coords_rep$lat <- base::round(ocean_coords_rep$lat, digits = 3)
ocean_coords_rep$lon <- base::round(ocean_coords_rep$lon, digits = 3)

land_coords <- base::data.frame(GEOmap::Lxyz2ll(list(x = mesh_land$loc[,1], y = mesh_land$loc[,2], z = mesh_land$loc[,3])))
land_coords$index = base::seq(1:base::length(land_coords$lat))
land_coords_rep = base::data.frame(base::sapply(land_coords,base::rep.int,time=n_years))
land_coords_rep$time = base::rep(1:n_years, each = base::length(land_coords$lon))
land_coords_rep$lat <- base::round(land_coords_rep$lat, digits = 3)
land_coords_rep$lon <- base::round(land_coords_rep$lon, digits = 3)

gic_coords <- base::data.frame(GEOmap::Lxyz2ll(list(x = mesh_GIC$loc[,1], y = mesh_GIC$loc[,2], z = mesh_GIC$loc[,3])))
gic_coords$index = base::seq(1:base::length(gic_coords$lat))
gic_coords_rep = base::data.frame(base::sapply(gic_coords,base::rep.int,time=n_years))
gic_coords_rep$time = base::rep(1:n_years, each = base::length(gic_coords$lon))
gic_coords_rep$lat <- base::round(gic_coords_rep$lat, digits = 3)
gic_coords_rep$lon <- base::round(gic_coords_rep$lon, digits = 3)

IceSheets_coords <- base::data.frame(GEOmap::Lxyz2ll(list(x = mesh_IceSheets$loc[,1], y = mesh_IceSheets$loc[,2], z = mesh_IceSheets$loc[,3])))
IceSheets_coords$index = base::seq(1:base::length(IceSheets_coords$lat))
IceSheets_coords_rep = base::data.frame(base::sapply(IceSheets_coords,base::rep.int,time=n_years))
IceSheets_coords_rep$time = base::rep(1:n_years, each = base::length(IceSheets_coords$lon))
IceSheets_coords_rep$lat <- base::round(IceSheets_coords_rep$lat, digits = 3)
IceSheets_coords_rep$lon <- base::round(IceSheets_coords_rep$lon, digits = 3)

