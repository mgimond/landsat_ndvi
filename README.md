# landsat_ndvi

This repo houses an R-script that converts multiple red and near-infrared Landsat bands to normalized difference vegetation index (NDVI) rasters. This script is currently used in a couple of undergraduate courses. Since this script has not been fully vetted, caution should be taken if used for research.

## Input files stored in the ./input folder:

* Separate red and near-infrared Landsat raster files (note that the band numbers associated with these files are not the same across Landsat platforms). Note that you do **not** need to store all other Landsat bands in this folder.
* Landsat metadata files (usually ending in *_MTL.txt*). One file per Landsat scene.

## Output files created in the ./output folder:

* NDVI rasters (for each input Landsat scene) rescaled to 1000. The rescaling allows the raster to be stored as an integer thus significantly reducing the file size.

## Sample dataset
A sample set of landsat files (subsets of their original extents) are available in the sample.zip file. Simply unzip this file in the folder that houses the landsatNDVI.R script then run the script. Output will be saved in the ./output folder which is initially empty. 



