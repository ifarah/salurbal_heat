from osgeo import gdal
import numpy as np
import numpy.ma as ma
import pandas as pd
from netCDF4 import Dataset
from datetime import date, timedelta
import multiprocessing as mp
import glob

def per_nc(era5_land_nc):
    final = pd.DataFrame()
    idx = 0
    era5land = Dataset(era5_land_nc, mode='r')
    t2m = era5land.variables['t2m']
    year = int(era5_land_nc.split('\\')[-1].split('_')[0])
    Q = int(era5_land_nc.split('\\')[-1].split('_')[-1][1])

    # read in the area-weighted mask file, specifiy the exact location where you have these files stored
    masksfile = glob.glob(r'..\processing\irene_ft_raster\UX\by_area_weight\*_area_weight.tif') #change here  # choose between AD UX L2

    # read in worldpop raster file,, specifiy the exact location where you have this stored
    worldpop = gdal.Open(r'salurbal_heat-master\WorldPop\worldpoppers_bueno_convec.tif').ReadAsArray()
    worldpop[worldpop < 0] = 0 # set no data area with 0 weight

    # open the GUF urban footprint file, sepcify the exact location where you have this file stored
    guf = gdal.Open(r'salurbal_heat-master\GUF\guf_panperu_vectorized_complete.tif').ReadAsArray()
    guf[guf < 0] = 0 # set no data area with 0 weight
    f_date = date(year, 3 * Q - 2, 1)

    for maskfile in masksfile: # for each area unit
        SALID = maskfile.split('\\')[-1].split('_')[0]
        mask = gdal.Open(maskfile).ReadAsArray() # read in the area-weight mask
        mask = ma.masked_array(mask, mask < 0)  # generate a mask for the area unit
        days = t2m.shape[0]

        # get population or GUF weight
        if SALID.startswith('105') or SALID.startswith('206'): # decide if the area unit is in Peru or Panama
            guf_wp = guf # if yes, use GUF as the "population" weight
        else:
            guf_wp = worldpop # otherwise, use world pop as the weight

        mask1d = mask[mask.mask == False]  # get area weight = how much grid cells intersecting with area units
        guf_wp1d = guf_wp[mask.mask == False]*mask1d  # get area weight * population (GUF) weight

        # get the total weights
        total_area = np.sum(mask)
        total_guf_wp1d = np.sum(guf_wp1d)

        # normalize individual weight by the total weights
        mask1d = mask1d / total_area
        guf_wp1d = guf_wp1d/total_guf_wp1d

        for day in range(0, days): # iterate through the days
            t2m_day = t2m[day] - 273.15
            t2m_day = t2m_day[mask.mask == False]
            # apply weights to temperature
            weighted_t2m = t2m_day * mask1d
            weighted_t2m_gufwp = t2m_day * guf_wp1d
            date2 = f_date + timedelta(days=day)
            final.loc[idx, 'SALID'] = SALID
            final.loc[idx, 'date'] = date2
            # get the area-weighted temperature (weight = how much a grid cell intersects with boundary)
            final.loc[idx, 'UXtemp_pw'] = np.sum(weighted_t2m_gufwp) #change here # choose between AD UX L2
            # get the area and population -weighted temperature (weight = area weight * population weight)
            final.loc[idx, 'UXtemp_x'] = np.sum(weighted_t2m) #change here
            idx = idx + 1
    # write the results to a CSV per year-quarter
    final.to_csv(r'..\40Metrics_calculation\ERA5land_temperature\v2\L1UX_%s_%s.csv'%(year,Q)) #change here # choose between L1AD L1UX L2

# run the calculations in parellel
if __name__ == "__main__":
    era5_land_ncs = glob.glob(r'..\processing\ERA5land_fill\Results_v2\*_Q*.nc')

    pool = mp.Pool(14)
    results = pool.map(per_nc,era5_land_ncs)
    pool.close()






