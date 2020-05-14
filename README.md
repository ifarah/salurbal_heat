# UNDER CONSTRUCTION



# Temperature data from 2001 - 2015 from ERA5Land
## Sample of 371 Latin American Cities 

Many researchers involved in SALURBAL are interested in using historical temperature reanalysis data. ERA5Land data neglects pixels that have more than 50% water. However, many cities across the world are situated next to the ocean. Since we were losing information, our team interpolated data from ERA5 (at a 30 km x 30 km resolution) and imputed it at the ERA5Land 9km x 9km resolution, filling the gaps from those "missing" pixels. (Here include Yang's workflow). 

Once we got the "complete" data for all of our 371 cities, we decided to estimate temperatures at an AD (level 2) level. In order to get temperature data at the AD city level, we weight the temperature pixels by population data (using 100m x 100m WorldPop data for 2010).

Note: Since population data is not as accuradte for Panama and Peru, we weight temperature by urban footprint data instead for cities in those countries.

The workflow is as follows:  
1. Create a vector from the temperature raster (in 9km x 9km cells). 
2. Import the vector and the AD boundaries into Google Earth Engine (GEE).
3. Process WorldPop data in GEE, estimating the number of people by the 9km x 9km.
4. Export the 9km x 9km pixels with population data into R.
5. Carry out the population weight: Since both the temeprature and the population data are at the same resolution, we can carry out the estimation. 

For data processing, please view the Jupyter notebook and the folder containing the [temperature data here](https://github.com/ifarah/salurbal_era5_ncdc/blob/master/data_processing.ipynb) and [raw table 1 here](https://github.com/ifarah/salurbal_era5_ncdc/blob/master/output/cor_table.xlsx).

### Access to raw data:
- [ERA5 hourly data on single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
- [WorldPop data](https://www.worldpop.org/project/categories?id=3)
- [Global Urban Footprint data]()

### Access to imputed data:
- [ERA5Land imputed data]()





