#default values for model parameters
mic.fact = 1.0
# this version tuned on long term data
fit.par = c(mic_offset = 0.27179, kdissolution = 2.62852, kdepoly = 0.25887, 
            kdeath_and_exudates = 0.50878, kdesorb = 0.00432, ksorb = 0.23789, 
            kmicrobial_uptake = 1.3858, cue_0 = 0.33668, mcue = 0.008, mclay = 5e-05)

# this version tuned on clm and hwsd
fit.par = c(mic_offset = 0.76519, kdissolution = 4.77183, kdepoly = 1.17511, kdeath_and_exudates = 0.03843, kdesorb = 5e-04, ksorb = 0.1, kmicrobial_uptake = 1.60289, cue_0 = 0.15019, mcue = 0, mclay = 0.03268)
  
fit.lower = c(mic_offset = 0.0, kdissolution = 0.3, kdepoly = 0.008, 
              kdeath_and_exudates = 0.02, kdesorb = 0.0005, ksorb = 0.1, kmicrobial_uptake = 0.1, 
              cue_0 = 0.15, mcue = 0, mclay = 0)
fit.upper = c(mic_offset = 1.0, kdissolution = 6, kdepoly = 1.5, 
              kdeath_and_exudates = 1.0, kdesorb = 0.01, ksorb = 1.0, kmicrobial_uptake = 2.0, 
              cue_0 = 0.6, mcue = 0.008, mclay = 1/23)

  
if (time_step == 'day') {
  fit.par[2:7] = fit.par[2:7] * 12/365
  fit.lower[2:7] = fit.lower[2:7] * 12/365
  fit.upper[2:7] = fit.upper[2:7] * 12/365
}
if (clm.fit) fit.par = c(mic_offset = 0.76519, kdissolution = 0.15688, kdepoly = 0.03863, kdeath_and_exudates = 0.00126, kdesorb = 2e-05, ksorb = 0.00329, kmicrobial_uptake = 0.0527, cue_0 = 0.15019, mcue = 0, mclay = 0.03268)

# function to evaluate list of expressions on data.table, sequentially
with.dt = function(dt, expr){
  for (j in 1:length(expr)) set(dt, NULL, names(expr)[j], dt[, eval(expr[[j]])])
}

# function to calculate temperature dependence of decomposition rate
fT.PrimC = function(temp, method = 'Century2', t.ref = 30){
  switch(method,
    RothC = fT.RothC(temp),
    Century1 = fT.Century1(temp) * fT.RothC(t.ref) / fT.Century1(t.ref),    
    Century2 = fT.Century2(temp) * fT.RothC(t.ref) / fT.Century2(t.ref),
    stop('Unrecognised temperature function'))
}

wsat = function(sand, silt, clay) 0.6658*silt + 0.1567*sand - 0.0079*silt^2 - 12.31121/sand - 
  6.4756*log(sand) - 0.0038*clay*silt + 0.0038*clay*sand - 0.0042*silt*sand + 52.7526

wfield = function(sand, silt, clay) 118.932*clay + 119.0866*silt + 119.1104*sand + 162.31731/clay - 
  46.21921/silt-5.12991/sand  + 18.1733*log(clay) + 0.0013*clay*silt + 0.0022*silt*sand - 11939.3493

wwilt = function(sand, silt, clay) 92.3851 - 1.5722*silt - 0.5423*sand - 0.0072*clay^2 + 0.0072*silt^2 - 
  0.0059*sand^2 + 160.14591/clay  +  6.60011/sand + 0.0022*clay*silt - 0.0039*clay*sand


unrotate = function(x) { # inverse of rotate: convert to 0 to 360 longitudes
  raster::shift(rotate(raster::shift(x, 180)), 180)
}

pl = function(x){
  xname = deparse(substitute(x))
  plot(rotate(x), ylim=c(-90,90), main = xname)
  lines(world, lwd=0.3)
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(rev(brewer.pal(9, "GnBu")))

sf = scale_fill_gradientn(
  colours = myPalette(50),
  breaks = (c(0,2,4,6,8,10,15,20)),
  limits = (c(0,25)),
  na.value = 'white',
  name = expression("SOC (kg m"^-2*")"))

t <-theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='bottom',
          legend.title.align=1,
          legend.title = element_text(size=10),
          legend.direction="horizontal",
          plot.margin = unit(c(0,0,0,0), "cm"))

g.map = function (r, map_theme=t, map_cols=sf){
  gplot((r)) +
    geom_raster(aes(fill=value)) +
    geom_path(data = wmap_df,
              mapping = aes(long, lat, group=group),
              size=0.4) +
    geom_path(data=bbox_df,
              mapping=aes(long,lat, group=group), size=0.4) +
    guides(fill= guide_colorbar(barwidth=15)) +
    map_cols + map_theme
}


elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

get_data = function(data_path) {
  files <- list.files(path=data_path, pattern="*.socdat", full.names=T, recursive=FALSE)
  experiments <<- basename(file_path_sans_ext(files))
  all.data <- setDT(ldply(files, function(fn)  data.frame(read.table(fn, header = T, sep = ","),exp=basename(file_path_sans_ext(fn)))))
  all.data$added.bio <- as.double(all.data$added.bio)
  all.data <- merge (all.data, exp.const, by="exp")
  setnames(all.data, "init.soc", "soc")
  all.data
}

# function to intialise values in soc dataframe
initialise.dayprimc.data <- function(soc.data) {
  t.ref =30
  roth.cent.ratio = fT.RothC(t.ref) / fT.Century2(t.ref)
  
  #convert monthly to daily data
  start.date = as.Date('1800-01-01')
  soc.data[, date := start.date %m+% months(month-1)]
  daily.soc.data = soc.data[, .(date = seq(date[1], date[.N], by='1 day')), by = exp]
  daily.soc.data = merge(soc.data, daily.soc.data, by=c('exp', 'date'), all.y = T)
  setorder(daily.soc.data, exp, date)

  #interpolate missing daily values
  daily.soc.data[, month := elapsed_months(date, start.date) + 1]
  daily.soc.data[, day := 1:.N, by = exp]
  daily.soc.data[, temp := na.spline(temp), by = exp]
  daily.soc.data[, precip := na.spline(precip)*12/365, by = exp]
  daily.soc.data[, pet := na.spline(pet)*12/365, by = exp]
  daily.soc.data[, cover := na.locf(cover), by = exp]
  daily.soc.data[, clay := na.locf(clay), by = exp]
  daily.soc.data[, depth := na.locf(depth), by = exp]
  daily.soc.data[, soc := na.locf(soc), by = exp]
  daily.soc.data[, soc.d13c := na.locf(soc.d13c), by = exp]
  daily.soc.data[, max_tsmd := na.locf(max_tsmd), by = exp]
  daily.soc.data[, max_tsmd := na.locf(max_tsmd), by = exp]
  daily.soc.data[is.na(added.dpm), added.dpm := 0.0]
  daily.soc.data[is.na(added.rpm), added.rpm := 0.0]
  daily.soc.data[is.na(added.bio), added.bio := 0.0]
  daily.soc.data[is.na(added.sta), added.sta := 0.0]
  daily.soc.data[, added.d13c := na.locf(added.d13c), by = exp]

  with.dt(daily.soc.data, expression(
    dpm = soc * 0.001,
    rpm = soc * 0.081,
    doc = 0.0,
    bio = soc * 0.0331,
    sta = soc * (1-0.001-0.081-0.0331),
    atsmd = 0.0,
    sat = 0.0,
    #field = 0.0,
    h2o = 0.0,
    a = fT.PrimC(temp),
    b = 1,
    c = ifelse (cover==0, 1, 0.6),
    mic = 1,
    x = 1.67*(1.85+1.6*exp(-0.0786*clay)),
    added.doc = 0.0,
    decomp.dpm = dpm * (1 - exp(-k1 * a * b * c * mic)),
    decomp.rpm = rpm * (1 - exp(-k2 * a * b * c * mic)),
    decomp.doc = doc * (1 - exp(-kdoc * a * b * c * mic)),
    decomp.bio = bio * (1 - exp(-k3 * a * b * c * mic)),
    decomp.sta = sta * (1 - exp(-k4 * a * b * c * mic)),
    decomp.tot = decomp.dpm + decomp.rpm + decomp.bio + decomp.sta,
    mineralise.doc = decomp.doc * frac.co2,
    mineralise.tot = mineralise.doc,
    mineralise.cum = 0.0,
    dpm.d13c = soc.d13c,
    rpm.d13c = soc.d13c,
    doc.d13c = soc.d13c,
    bio.d13c = soc.d13c,
    sta.d13c = soc.d13c,
    co2.d13c = 0.0
  ))
  daily.soc.data
}

Mode <- function(x, na.rm = FALSE) {
  if (na.rm) x = x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# function to intialise values in soc dataframe
initialise.monthly.dayprimc.data <- function(soc.data) {
  with.dt(soc.data, expression(
    day = 1L,
    dpm = soc * 0.001,
    rpm = soc * 0.081,
    doc = 0.0,
    bio = soc * 0.0331,
    sta = soc * (1-0.001-0.081-0.0331),
    atsmd = 0.0,
    sat = 0.0,
    #field = 0.0,
    h2o = 0.0,
    # a = fT.PrimC(temp, method = 'Century1', t.ref = 28),
    a = fT.PrimC(temp),
    c = ifelse (cover==0, 1, 0.6),
    mic = 1,
    added.doc = 0.0,
    added.rpm = added.rpm + added.hum,
    added.hum = 0.0,
    dpm.d13c = soc.d13c,
    rpm.d13c = soc.d13c,
    doc.d13c = soc.d13c,
    bio.d13c = soc.d13c,
    sta.d13c = soc.d13c,
    co2.d13c = 0.0
  ))
  soc.data
}

wmean.bricks = function(b.list, w) {
  b.wmean = b.list[[1]] * w[1]
  for (i in seq_along(w)[-1]){
    b.wmean = b.wmean + b.list[[i]] * w[i]
  }
  b.wmean / sum(w)
}
