switch(Sys.info()[['sysname']],
       Windows = {base_path = '~/../Box Sync/NIFA'},
       Linux   = {base_path = '~/Documents/NIFA'})
CASA_path = paste0(base_path, '/global_data/CASA_TESTBED')
R_path = paste0(base_path, '/R')
time_step = 'day'
clm.fit = TRUE

setwd(R_path)
source('SCAM packages.r')
source('SCAM Functions 1.02.r')
sourceCpp('SCAM1.01.cpp')

rasterOptions(maxmemory = 1e+09, progress = "text")
setwd(CASA_path)

##################################################################
# coastline map to overplot on raster plots
##################################################################
data(world2MapEnv)
world = map("world", fill=TRUE, col="transparent", plot=FALSE)
world = map2SpatialPolygons(world, world$names, CRS("+proj=longlat +ellps=WGS84"))
world2 = map("world2", fill=TRUE, col="transparent", plot=FALSE)
world2 = map2SpatialPolygons(world2, world2$names, CRS("+proj=longlat +ellps=WGS84"))

wmap = readOGR(dsn="world_outlines/ne_110m_land.shp", layer="ne_110m_land")
bbox = readOGR("world_outlines/ne_110m_wgs84_bounding_box.shp", layer="ne_110m_wgs84_bounding_box")
bbox = spTransform(bbox, CRS("+proj=robin"))
bbox_df = fortify(bbox)
wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df <- fortify(wmap_robin)

##################################################################
# Load litter layers and determine which raster cells contain data
##################################################################
cflux = ncdf4::nc_open("casaclm_pool_flux_1901_daily.nc")
litter.dpm = brick("casaclm_pool_flux_1901_daily.nc", lvar = 3, varname = "litInptMet") #gC m-2 day-1
litter.rpm = brick("casaclm_pool_flux_1901_daily.nc", lvar = 3, varname = "litInptStruc") #gC m-2 day-1
litter.dpm = litter.dpm * 100*100 / 1e6 # Mg C / ha / day
litter.rpm = litter.rpm * 100*100 / 1e6 # Mg C / ha / day
litter.tot = litter.dpm + litter.rpm
Nlat  = dim(litter.dpm)[1]
Nlon  = dim(litter.dpm)[2]
Nday  = dim(litter.dpm)[3]
#use this commented version to include cells that have only zeros
# has.dpm = sort(which(getValues(all(!is.na(litter.dpm)))))
# has.rpm = sort(which(getValues(all(!is.na(litter.rpm)))))
has.dpm = sort(which(getValues(max(litter.dpm, na.rm=T)) > 0))
has.rpm = sort(which(getValues(max(litter.rpm, na.rm=T)) > 0))
has.litter = sort(union(has.dpm, has.rpm))

##################################################################
# Load Soil layers
##################################################################
soil = ncdf4::nc_open('HWSD_SOIL_CLM_RES.nc4')
PCT_SAND = unrotate(raster('HWSD_SOIL_CLM_RES.nc4', varname='PCT_SAND')) # %
PCT_CLAY = unrotate(raster('HWSD_SOIL_CLM_RES.nc4', varname='PCT_CLAY')) # %
# AWT_SOC = unrotate(raster('HWSD_SOIL_CLM_RES.nc4', varname='AWT_SOC'))  # area weighted soc to 1 m from HWSD
BD       = unrotate(raster('HWSD_SOIL_CLM_RES.nc4', varname='BULK_DEN')) # Mg/m3

# aggregate and resample soil rasters to same dimensions as atmosphere
PCT_SAND = aggregate(PCT_SAND, fact = round(dim(PCT_SAND)[1:2] / dim(litter.tot)[1:2]))
PCT_SAND = resample(PCT_SAND, litter.tot)
PCT_CLAY = aggregate(PCT_CLAY, fact = round(dim(PCT_CLAY)[1:2] / dim(litter.tot)[1:2]))
PCT_CLAY = resample(PCT_CLAY, litter.tot)
BD = aggregate(BD, fact = round(dim(BD)[1:2] / dim(litter.tot)[1:2]))
BD = resample(BD, litter.tot)
PCT_SILT = 100 - PCT_SAND - PCT_CLAY

clm_file = 'clm4_5_12_r191_CLM45spHIST_CRU.clm2.h1.1901-01-01-00000.nc'
clm.soil.depths = c(1.75, 4.51, 9.06, 16.56, 28.91, 49.29, 82.89,138.28, 229.61, 380.19) #http://onlinelibrary.wiley.com/doi/10.1002/2013WR014586/pdf
clm_max_depth = 5

#H2OSOI volumetric soil water (vegetated only) (m3/m3). Use mean of top 6 soil laters (weighted by layer thickness)
H2OSOI = lapply(1:clm_max_depth, function(lev) brick(clm_file, lvar = 3, varname = "H2OSOI", level = lev))
H2OSOI = wmean.bricks(H2OSOI, c(clm.soil.depths[1], diff(clm.soil.depths[1:clm_max_depth])))

# TSOI: soil temperature (vegetated landunits only) (K)
TSOI = lapply(1:clm_max_depth, function(lev) brick(clm_file, lvar = 3, varname = "TSOI", level = lev))
TSOI = wmean.bricks(TSOI, c(clm.soil.depths[1], diff(clm.soil.depths[1:clm_max_depth])))
TSOI = TSOI - 273.13

# Calculate saturation water capacity
# Scheinost 3 eqn: http://www.soil.tu-bs.de/pubs/reprints/02_StudienarbeitLipsius.pdf
# theta_sat = 0.86 - 0.34*BD + 0.14*PCT_CLAY/100 # BD in g/cm3, clay in %, theta_sat in m3/m3
# theta_sat = max(H2OSOI, theta_sat, na.rm=TRUE)
theta_sat = wsat(PCT_SAND, PCT_SILT, PCT_CLAY)
theta_field = wfield(PCT_SAND, PCT_SILT, PCT_CLAY)
theta_wilt = wwilt(PCT_SAND, PCT_SILT, PCT_CLAY)

#########################################
# soc map to compare model predictions to
#########################################
# use two levels for 1m, only one level for topsoil
# soc.hwsd = raster('HWSD_SOIL_CLM_RES.nc4', varname='AWT_SOC', level=1) + raster('HWSD_SOIL_CLM_RES.nc4', varname='AWT_SOC', level=2) 
soc.hwsd = raster('HWSD_SOIL_CLM_RES.nc4', varname='AWT_SOC', level=1)
soc.hwsd.clm = resample(unrotate(soc.hwsd), BD)
soc.hwsd.plot = projectRaster(soc.hwsd, crs=CRS("+proj=robin +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
soc.hwsd.plot = mask(soc.hwsd.plot, wmap_robin)
g.soc = g.map(soc.hwsd.plot); g.soc
# ggsave(gsub(':','_', paste0('HWSD SOC 30cm ', Sys.time(), '.pdf')), height=5, width=8)


###############################
# experiment data creation
###############################

#number months and years through experiment
SpinUpYears = 250
clm.litter.quality = TRUE
expt.data = data.table(day = 1:(365 * SpinUpYears))
expt.data[, month :=  as.integer(cut(day %% 365, 12))]
expt.data[, precip := NA] # set precip and pet to NA, because we are using soil hydrology for water, rather than climate
expt.data[, pet := NA]
expt.data[, soc.d13c := -25]
expt.data[, added.d13c := -25]

Nlitter = length(has.litter)
soc = litter.tot[[1]]
soc[] = NA
soc.init = unrotate(soc.hwsd) * 10 # guess value to start spin up (Mg / ha)
soc.init = aggregate(soc.init, max(round(dim(soc.hwsd)[1:2]/ dim(litter.tot)[1:2])))
soc.init = resample(soc.init, litter.tot)
soc.init[is.na(soc.init)] = 40



##################################################################################
# loop to generate data structure for each gridcell and run PrimC model on it
##################################################################################
pb <- tkProgressBar(title = "Calculating Global SOC", min = 0, max = Nlitter, width = 400)
for (i in 1:Nlitter) {
  cell = has.litter[i]
  #  climate
  expt.data[, temp := extract(TSOI, cell)[1,]]
  expt.data[, sat := extract(theta_sat, cell)]
  expt.data[, h2o := extract(H2OSOI, cell)[1,]]
  #  carbon inputs
  expt.data[, fym := 0] # should set fym as fraction of litter, depending on land cover (rangeland or not)
  expt.data[, litter := extract(litter.tot, cell)[1,]]
  expt.data[, added.bio := 0.0]

  # use either CLM or rothc litter quality
  if (clm.litter.quality){
    expt.data[, added.dpm := extract(litter.dpm, cell)[1,] + fym*0.49]
    expt.data[, added.rpm := extract(litter.rpm, cell)[1,] + fym*0.49]
    expt.data[, added.hum := fym*0.02]
  } else{
    expt.data[, added.dpm := litter * 0.59 + fym*0.49]
    expt.data[, added.rpm := litter * 0.41 + fym*0.49]
    expt.data[, added.hum := fym*0.02]
  }

  #  soil
  clay = PCT_CLAY[cellFromXY(PCT_CLAY, xyFromCell(litter.dpm, cell))] 
  max_tsmd = -(20.0 + 1.3*clay -0.01*clay^2)
  expt.data[, soc := soc.init[cell]] # guess value to start spin up

  #  Initialise pools, flows, and rate modifers
  with.dt(expt.data, expression(
    cover = 1,
    dpm = soc * 0.001,
    rpm = soc * 0.081,
    doc = 0.0,
    bio = soc * 0.0331,
    sta = soc * (1-0.001-0.081-0.0331),
    atsmd = 0.0,
    a = fT.PrimC(temp, method = 'RothC'),
    c = 0.6,
    added.doc = 0.0
  ))

  scamc = as.data.table(
    scam(expt.data, 
         mic_fact = 1.0, 
         mic_offset = fit.par['mic_offset'],
         kdissolution = fit.par['kdissolution'],
         kdepoly = fit.par['kdepoly'],
         kdeath_and_exudates = fit.par['kdeath_and_exudates'],
         kdesorb = fit.par['kdesorb'],
         ksorb = fit.par['ksorb'],
         kmicrobial_uptake = fit.par['kmicrobial_uptake'],
         cue_0 = fit.par['cue_0'],
         mcue = fit.par['mcue'], 
         mclay = fit.par['mclay'],       
         clay, 
         max_tsmd, 
         use_atsmd = FALSE, 
         monthly = FALSE)
    )

  soc[cell] = scamc[.N, soc]
  setTkProgressBar(pb, i, label=paste(round(i/Nlitter*100, 0),"% done"))
}
writeRaster(soc, gsub(':','_', paste0('PrimC output', Sys.time(), '.tif')))
close(pb) # shut the progress bar


##################################################################
# Draw map
##################################################################
#soc = raster(file.choose())
soc.plot = soc / 10 # kg/m2
soc.plot = rotate(soc.plot)
soc.plot = disaggregate(soc.plot, max(round(dim(soc.hwsd.plot)[1:2]/ dim(soc.plot)[1:2])), 'bilinear')
soc.plot = projectRaster(soc.plot, crs=CRS("+proj=robin +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
soc.plot = mask(soc.plot, wmap_robin)
soc.plot = resample(soc.plot, soc.hwsd.plot)

g.soc = g.map(soc.plot)
g.soc
ggsave(gsub(':','_', paste0('PrimC output ', Sys.time(), '.pdf')), height=5, width=8)


##################################################################
# Difference from HWSD
##################################################################
soc.diff =  (soc.plot - soc.hwsd.plot) #/ soc.hwsd.plot
# g.map(soc.hwsd.plot)
# g.map(soc.plot)
# rng = cellStats(soc.diff,range)
# cuts = matrix (c(-Inf,-20,-10,-5,-1,0,1, 5,10,20,
#                  -20 ,-10, -5,-1, 0,1,5,10,20,Inf,
#                  1   , 2,   3, 4, 5,6,7,8,9,10), ncol=3)
# soc.diff = reclassify(soc.diff, cuts)
# cols = rev(brewer.pal(10,"RdBu"))
# names(cols) = cuts[,3]
library(scales)

gplot(soc.diff) +
  geom_raster(aes(fill=value)) +
  geom_path(data = wmap_df, mapping = aes(long, lat, group=group),size=0.4) +
  geom_path(data=bbox_df, mapping=aes(long,lat, group=group), size=0.4) +
  scale_fill_gradient2(na.value = 'white', 
                       low = muted("blue"), mid = 'grey98', high = muted("red"),
                       limits = c(-20,20), 
    name = expression("difference in SOC")) + 
  theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position='right',
            legend.title.align=1,
            legend.title = element_text(size=10),
            legend.direction="vertical",
            plot.margin = unit(c(0,0,0,0), "cm"))





soc.diff =  (soc/10 - resample(soc.hwsd, soc)) / resample(soc.hwsd, soc)
soc.diff
PCT_CLAY
ggplot(data.frame(clay=values(PCT_CLAY), cdiff=values(soc.diff)), aes(clay,cdiff))+
  geom_point() +
  geom_smooth(method = 'lm')
