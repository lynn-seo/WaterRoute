# -*- coding: utf-8 -*-
#TODO: need to buffer the tile with cells from neighbouring tiles (2 extra pixels exactly)
import flow1_c
import gdal, gdalconst
import numpy as np
import scipy.ndimage as ndimage

import networkx as nx


def pixelToXY(geoTransform,row,col):
    #from: http://www.gdal.org/classGDALDataset.html#f9593cc241e7d140f5f3c4798a43a668
    X = geoTransform[0] + col*geoTransform[1] + row*geoTransform[2];
    Y = geoTransform[3] + col*geoTransform[4] + row*geoTransform[5];
    return (X,Y)

def run():
#open data
g = gdal.Open(r'H:\rea095\479\e124s27', gdalconst.GA_ReadOnly)
trans = g.GetGeoTransform()
dem = g.GetRasterBand(1).ReadAsArray(0,0,g.RasterXSize,g.RasterYSize)
#flow direction
fd = flow1_c.d8(dem)
#group up neighbouring sinks and label them (like regiongroup)
regions = np.zeros(fd.shape, dtype=int)
connected = np.array([[1, 1, 1],[1, 1, 1],[1, 1, 1]], dtype=bool)
num_sinks = ndimage.label(fd[2:-2,2:-2]==0, connected,regions[2:-2,2:-2])
#cells that recieve flow from the tile (outlets)
outlet_list = flow1_c.flow_band(fd,2,1)
# add region labels to outflowing points.
for i,coord in enumerate(outlet_list,start=num_sinks + 1):
    regions[coord] = i
    
#calculate catchments for outflowing points
catch = flow1_c.catchments(fd,regions)

#HINT: could process sinks that are fully internal to the tile now
#but leave that for now and implement later if optimisation is required
#sink processing is required for the merged tiles network and can deal with the lot

#gather together all the info that is required to make the merged tiles network.
#need
#outlet_list - using lat/long coords and labeled with catchment
#inlet_list - using lat/long coords and labeled with catchment
#lowest saddles between catchments using lat/long coords and labeled with both catchments
#sinks or flat areas on the edge of the tile

#note: could use some sort of tileTheCellFallsOn index rather then lat/long

inlet_list = flow1_c.flow_band(fd,1,2)
inlet_nodes = [(pixelToXY(trans,*x),catch[x]) for x in inlet_list]
outlet_nodes = [(pixelToXY(trans,*x),catch[x]) for x in outlet_list]
saddles_dict = flow1_c.saddles(dem,catch)
#The saddles_dict should have catchment 0 for any saddles on edge of tile


DG=nx.DiGraph()
out_dict = {}
for coord,id in outlet_nodes:
    out_dict[id] = coord
    if not coord in DG.node:
        DG.add_node(coord)
            
for coord,id in inlet_nodes:
    if not coord in DG.node:
        DG.add_node(coord)
    DG.add_edge(key,out_dict[id]) #doesnt work for sinks!! or flats
    
    
    
    