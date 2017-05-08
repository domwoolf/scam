my_packages = c('data.table', 'ggplot2', 'cowplot', 'plyr', 'tools', 
                'Rcpp', 'lubridate', 'tcltk2', 'maps', 'maptools', 'rgdal', 
                'raster', 'rasterVis', 'ncdf4', 'SoilR', 'zoo')

for (package in my_packages)  {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
