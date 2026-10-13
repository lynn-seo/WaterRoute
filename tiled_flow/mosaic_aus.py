# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\rea095\.spyder2\.temp.py
"""

import numpy as np
import os
from osgeo import gdal

halfcell = 1.0/ 3600
halfcell = cell / 2
top = -10 - halfcell
left = 113 - halfcell
geo = (left, cell, 0.0, top, 0.0, -cell)

proj = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]'

# dir = r'S:\rea095\dem_h_masked'
dir = r'\\Wron\gis\national\dem\SRTM_derived\1sec\V1_0_unofficial_or_prerelease\DEM_H\Tiles\\'
outArray = np.empty(((-10 - -44)*3600,(154-113)*3600),np.float32)
for i,x in enumerate(xrange(113,154)):
    print x
    for j,y in enumerate(reversed(range(-44,-10))):
        print y
        p = os.path.join(dir, 'e' + str(x) + r'\s' + str(abs(y)) + r'\e' + str(x) + '_s' + str(abs(y)))
        #print p
        if os.path.exists(p):
            print 'Got', p
            g = gdal.Open(p, gdal.gdalconst.GA_ReadOnly)
            outArray[j*3600:j*3600 + 3600,i*3600:i*3600 + 3600] = g.GetRasterBand(1).ReadAsArray()
        else:
            print 'No tile'
            outArray[j*3600:j*3600 + 3600,i*3600:i*3600 + 3600] = -9999
        
format = "GTiff"
driver = gdal.GetDriverByName( format )
outDataset = driver.Create(r'H:\JohnGallant\demh.tif',outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Float32)
outDataset.SetGeoTransform(geo)
outDataset.SetProjection(proj)
outBand = outDataset.GetRasterBand(1)
outBand.SetNoDataValue(-9999)
outBand.WriteArray(outArray)
del(outDataset)
