# -*- coding: utf-8 -*-
"""
Created on Fri Jul 08 15:48:01 2011

@author: rea095
"""
import flow1_c
import gdal, gdalconst
import numpy as np
from itertools import product

#TODO: to be calculated in site finding script (?)
dam_top_elev = 310.0         
damx = 126.020827
damy = -17.629144
#TODO:calculate the area from the gdal transform data
cell_area = 30*30 #in metres

#open datasets
g = gdal.Open(r'\\wron\project\WINA\3_work\3_trev\3_f_demh_1v0\f_filh1v0_clp', gdalconst.GA_ReadOnly)
trans = g.GetGeoTransform()
dem = g.GetRasterBand(1).ReadAsArray(0,0,g.RasterXSize,g.RasterYSize)
nodata1 = g.GetRasterBand(1).GetNoDataValue()
g = gdal.Open(r'\\wron\project\WINA\3_work\3_trev\3_f_demh_1v0\f_fd_fil', gdalconst.GA_ReadOnly)
#trans = g.GetGeoTransform()
fd = g.GetRasterBand(1).ReadAsArray(0,0,g.RasterXSize,g.RasterYSize)
fd = fd.astype(numpy.int)
nodata = g.GetRasterBand(1).GetNoDataValue()

#In a north up image, trans[1] is the pixel width, and trans[5] is the pixel height.
#The upper left corner of the upper left pixel is at position (trans[0],trans[3]).
#calc row and col of dam
d_col = int((damx - trans[0]) / trans[1])
# 10875.477199864492 to 10875
d_row = int((damy - trans[3]) / trans[5])
# 5864.418400133739 to 5864
dam_base_elev = dem[d_row, d_col]
#the catchments code has a hard coded nodata value (oops)
#quick work around till that is fixed properly
fd[fd == nodata] = -9999
#now find the catchment
regions = np.zeros(fd.shape, dtype=int)
regions[d_row, d_col] = 1
catch = flow1_c.catchments(fd,regions)

#depth of water at full water level for cells within the catchment
d = dam_top_elev - dem[catch == 1]
#limit to only those cells that are flooded
d_depth = d[d > 0]

#cell height above base of dam
d_flood = (dam_top_elev - dam_base_elev) - d_depth

d_flood.sort()
#count of cells flooded proxy for area
area = np.arange(1,len(d_flood)+1,dtype='float64')  * cell_area
#get the change in height as the dam floods another cell
d_increase = ediff1d(d_flood, to_begin=0)
volume_step = d_increase * area
volume = np.cumsum(volume_step)
a = np.rec.array([d_flood,area,volume])
#np.savetxt('test.csv',a,delimiter=',')
#a.dump('output.dmp') #just output a pickle of the data

#someone else can work out what format they want it in



dam = np.where(catch ==1,dem < dam_top_elev,False)
#note: havent tested this bit of code outside australia.
xlist = []    
for x in [x * 0.05 for x in xrange(-180*20, 180*20)]:
    xmin = round((x - 0.05/2 - trans[0]) / trans[1] -0.5)
    xmax = round((x + 0.05/2 - trans[0]) / trans[1] -0.5)
    xmin = max(0,xmin)
    if xmax > 0 and xmin < dem.shape[0]:
        xlist.append ((xmin,xmax,x))

ylist = []
for y in [y * 0.05 for y in xrange(-90*20, 90*20)]:
    ymax = round((y - 0.05/2 - trans[3]) / trans[5] +0.5)
    ymin = round((y + 0.05/2 - trans[3]) / trans[5] +0.5)
    ymin = max(0,ymin)
    if ymax > 0 and ymin < dem.shape[1]:
        ylist.append ((ymin,ymax,y))

full_count = round(0.05 / trans[1])**2
resultlist = []
for b in product(xlist,ylist):
    #print b[0][0]:b[0][1],b[1][0]:b[1][1]
    dam_block = dam[b[0][0]:b[0][1],b[1][0]:b[1][1]]
    block = dem[b[0][0]:b[0][1],b[1][0]:b[1][1]]
    block = block[dam_block == True]
    block = block[block != nodata1]
    if len(block) > 0:
        med = np.median(block)
        ratio = len(block) / full_count
        name = ("%.2f" % b[0][2] + "%.2f" % b[1][2]).replace('-','_').replace('.','')
        resultlist.append((name,med,ratio))
