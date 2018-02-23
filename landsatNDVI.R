# AUTHOR: Manuel Gimond
# PURPOSE: Automate the computation of NDVI values from Landsat images 
# retrieved from # http://earthexplorer.usgs.gov.
# Radiance and reflectance formulas provided by the Landsat Handbook
# (http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html)
#
# Before running this script, make sure that all input bands reside in a folder
# called input which needs to reside at the same level as this script. You will
# also need to create a new folder called output (this is where the script will save
# the new NDVI rasters).
#

library(raster)

### Function computes fraction of mean distance to sun (adopted from 'Landsat' R package)
Sdist = function (adate) 
{
  edist <- julian(as.Date(adate), origin = as.Date(paste(substring(adate,
                                                                   1, 4), "12", "31", sep = "-")))[[1]]
  edist <- 1 - 0.016729 * cos((2 * pi) * (0.9856 * (edist - 4)/360))
  edist
}
### End function

# Some constants
E.red    = 1551  # Red solar spectral irradiance
E.nir    = 1044  # NIR solar spectral irradiance

# Get list of landsat metadata files in directory
flist <- list.files(path="./input",pattern = "\\MTL\\.txt$",full.names=T)

# Loop through each file in list
for (i in 1:length(flist) ){
  
  # Open metadata file
  inFile = read.delim(flist[i],header=F,sep="=", stringsAsFactors=F)
  inFile$V1 = gsub(' +', '', inFile$V1) # Remove spaces from 1st column
  
  # Identifying Satellite
  sat.id = gsub(' +', '',inFile$V2[inFile$V1 == "SPACECRAFT_ID"]) # get value and remove spaces
  
  # Extract parameters from metadata file
  if (sat.id == "LANDSAT_8") {
    sun.elev = as.numeric(inFile$V2[inFile$V1 == "SUN_ELEVATION"]) * pi / 180 # Solar elevation in radians
    lmax.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_4"])
    lmin.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_4"])
    lmax.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_5"])
    lmin.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_5"])
    Qmax.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_4"])
    Qmin.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_4"])
    Qmax.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_5"])
    Qmin.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_5"])
    L.date   = gsub(' +', '',inFile$V2[inFile$V1 == "DATE_ACQUIRED"]) # get value and remove spaces
    file.red = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_4"]) # get value and remove spaces
    file.nir = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_5"]) # get value and remove spaces
  }else if (sat.id ==  "LANDSAT_2"){
    sun.elev = as.numeric(inFile$V2[inFile$V1 == "SUN_ELEVATION"]) * pi / 180 # Solar elevation in radians
    lmax.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_5"])
    lmin.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_5"])
    lmax.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_6"])
    lmin.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_6"])
    Qmax.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_5"])
    Qmin.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_5"])
    Qmax.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_6"])
    Qmin.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_6"])
    L.date   = gsub(' +', '',inFile$V2[inFile$V1 == "DATE_ACQUIRED"]) # get value and remove spaces
    file.red = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_5"]) # get value and remove spaces
    file.nir = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_6"]) # get value and remove spaces
  }else{
    sun.elev = as.numeric(inFile$V2[inFile$V1 == "SUN_ELEVATION"]) * pi / 180 # Solar elevation in radians
    lmax.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_3"])
    lmin.red = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_3"])
    lmax.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MAXIMUM_BAND_4"])
    lmin.nir = as.numeric(inFile$V2[inFile$V1 == "RADIANCE_MINIMUM_BAND_4"])
    Qmax.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_3"])
    Qmin.red = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_3"])
    Qmax.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MAX_BAND_4"])
    Qmin.nir = as.numeric(inFile$V2[inFile$V1 == "QUANTIZE_CAL_MIN_BAND_4"])
    L.date   = gsub(' +', '',inFile$V2[inFile$V1 == "DATE_ACQUIRED"]) # get value and remove spaces
    file.red = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_3"]) # get value and remove spaces
    file.nir = gsub(' +', '',inFile$V2[inFile$V1 == "FILE_NAME_BAND_4"]) # get value and remove spaces
  }
  
  
  # Open raster layers
  red = raster(paste("./input/",file.red,sep=""))
  nir = raster(paste("./input/",file.nir,sep=""))
  
  # Set 0 to NA
  red[red == 0] = NA
  nir[nir == 0] = NA
  
  # Compute Radiance (L)
  red.L = (lmax.red - lmin.red)/(Qmax.red - Qmin.red) * (red - Qmin.red) + lmin.red
  nir.L = (lmax.nir - lmin.nir)/(Qmax.nir - Qmin.nir) * (nir - Qmin.nir) + lmin.nir
  
  # Compute Reflectance (R)
  red.R = (pi * red.L * Sdist(L.date)^2) / (E.red * cos(sun.elev))
  nir.R = (pi * nir.L * Sdist(L.date)^2) / (E.nir * cos(sun.elev))
  
  # Compute NDVI
  NDVI = (nir.R - red.R) / (nir.R + red.R)
  
  # Convert NDVI to integer values after multiplying by 1000
  NDVI_int = NDVI * 1000
  NDVI_int[] = round(NDVI_int[])
  # Save NDVI raster to file (as an ERDAS Imagine file)
  outFile = paste("./output/",unlist(strsplit(flist[i], "\\.|/"))[4],".img",sep="") # Set file name to landsat filename
  writeRaster(NDVI_int,filename=outFile,format="HFA", overwrite=T, datatype='INT2S')
  ## Save NDVI raster to file (as a GeoTIFF file)
  # outFile = paste("./output/",unlist(strsplit(flist[i], "\\.|/"))[4],".img",sep="") # Set file name to landsat filename
  # writeRaster(NDVI_int,filename=outFile,format="GTiff", overwrite=T, datatype='INT2S')
}

  