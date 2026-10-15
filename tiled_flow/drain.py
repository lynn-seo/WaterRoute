#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      rea095
#
# Created:     25/04/2012
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import pickle
import numpy as np
import scipy
from scipy import ndimage
import networkx as nx
from waterroute.tiled_flow import flow1_c
from osgeo import gdal, gdalconst
import sys

"""
a way to deal with large flat areas in a dem

large flat areas are problematic for finding a flow route over. there are many
techniques in the literture with varying levels of niceness.
This is just another one... there is no right answer just as long as it looks ok
and produces correct d8 flow.
It is aimed to produce a nice enough stream network

The real difference is it doesnt care if it produces sinks as the later
processing steps are designed to deal nicely with sinks so there is no use
avoiding them in this step.

This process is quite likely to still leave some small (a few pixels) flat areas
but these are considered ok as each flat cell can be considered as a seperate sink
in the later processing steps

one benefit is that if there are no flat areas left then
d8 only depends on its neighbouring cells

If working with data in a tiles form this process should be applied with buffers
to ensure flat areas are treated consistently.

"""
def flat2sinks(dem):
    ##useful reference http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    ##http://pythonvision.org/basic-tutorial#id6
    neighborhood = ndimage.morphology.generate_binary_structure(2,2)
    sinksNflats = ndimage.filters.minimum_filter(dem, footprint=neighborhood)==dem
    #introduce a slope away from the edges of a flat area forming sinks towards
    #the middle of the flat area
    sunk_dem = dem - ndimage.distance_transform_cdt(sinksNflats)
    return sunk_dem



coordOffsets = {1:(0,1), 2:(1,1), 4:(1,0), 8:(1,-1), 16:(0,-1), 32:(-1,-1), 64:(-1,0), 128:(-1,1)}
coordOffsets_dir = {(0,1):1, (1,1):2, (1,0):4, (1,-1):8, (0,-1):16, (-1,-1):32, (-1,0):64, (-1,1):128}

"""
using the saddle information change flow directions to make it flow out of this from catchment.
"""
def reverse_fd(fd,d,from_catch):
    if d['high_cat'] == from_catch:
        from_row,from_col = d['high_row'],d['high_col']
        to_row,to_col = d['low_row'],d['low_col']
    else:
        from_row,from_col = d['low_row'],d['low_col']
        to_row,to_col = d['high_row'],d['high_col']
    from_fd = fd[from_row,from_col]
    fd[from_row,from_col] = coordOffsets_dir[(to_row - from_row,to_col - from_col)]
    #print from_row,from_col
    to_row,to_col = from_row,from_col
    #process a sink and no further and stop if the flow is already in the right direction (flow between catchments is already correct like an outlet of a catchment)
    while from_fd != 0 and fd[from_row,from_col] != from_fd:
        #progresivly flip the flow
        row_off,col_off = coordOffsets[from_fd]
        from_row = to_row+row_off
        from_col = to_col+col_off
        from_fd = fd[from_row,from_col]
        fd[from_row,from_col] = coordOffsets_dir[(to_row - from_row,to_col - from_col)]
        #print from_row,from_col
        to_row,to_col = from_row,from_col

def gdalWrite(outArray, filename, trans):
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    if outArray.dtype in (np.int8, np.byte):
        outDataset = driver.Create(filename,outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Byte)
    elif outArray.dtype == np.int16:
        outDataset = driver.Create(filename,outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Int16)
    elif outArray.dtype == np.int32:
        outDataset = driver.Create(filename,outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Int32)
    elif outArray.dtype == np.float32:
        outDataset = driver.Create(filename,outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Float32)
    else:
        outDataset = driver.Create(filename,outArray.shape[1],outArray.shape[0],1,gdal.gdalconst.GDT_Float64)
    outDataset.SetGeoTransform(trans)
    outBand = outDataset.GetRasterBand(1)
    outBand.SetNoDataValue(-9999)
    outBand.WriteArray(outArray)
    del(outDataset)

def test():
    import os,datetime
    os.chdir(r's:\rea095')
    print (datetime.datetime.now().time())
    #testing code
    g = gdal.Open(r'S:\rea095\dem_h_mosaic_c.tif', gdalconst.GA_ReadOnly)
    trans = g.GetGeoTransform()
    dem = g.GetRasterBand(1).ReadAsArray(0,0,g.RasterXSize,g.RasterYSize)
    print (datetime.datetime.now().time())

    sunk_dem = flat2sinks(dem)
    print (datetime.datetime.now().time())
    sunk_dem = sunk_dem.astype(np.float32)
    print (datetime.datetime.now().time())
    #d8 flow direction
    fd = flow1_c.d8(sunk_dem)
    print (datetime.datetime.now().time())
    gdalWrite(fd,'fd.tif',trans)
    print (datetime.datetime.now().time())
    #unique label all catchments to sinks/outlets (remember d8 sets all edge cells as sinks)
    catch = flow1_c.sink_catchments(fd)  #use sink catchments in flow1_c
    print (datetime.datetime.now().time())
    #better yet rewrite catchments to use the flow order yielding code
    gdalWrite(catch,'catch.tif',trans)
    print (datetime.datetime.now().time())
    #bit of a hack so all outflowing catchments flow to the same catchment make all
    #the edges flow to a catchment labeled 0
    #coastal nodata and legit sinks should be set to zero as well.
    catch[0] = 0
    catch[-1] = 0
    catch[:,0] = 0
    catch[:,-1] = 0
    print (datetime.datetime.now().time())
    saddles_dict = flow1_c.saddles(sunk_dem,catch)
    print (datetime.datetime.now().time())
    #graph of catchments as nodes and saddles as edges
    G=nx.Graph()
    G.add_edges_from((u[0],u[1],d) for u,d in saddles_dict.items() )
    print (datetime.datetime.now().time())
    nx.write_gpickle(G, 'graph.gz')
    print (datetime.datetime.now().time())
    mst=nx.minimum_spanning_tree(G,weight='elev')
    P = nx.predecessor(mst,0)
    print (datetime.datetime.now().time())
    reversal = [reverse_fd(fd,mst[f][t[0]],f) for f,t in P.items() if len(t)==1]
    print (datetime.datetime.now().time())
    gdalWrite(fd,'fd_new.tif',trans)
    print (datetime.datetime.now().time())
    print ('done')
    # if you dont need a flow direction grid just do any processing in a network
    #generate stream lines if you need a vector representation.

"""
sinks is a numpy array where non zero values indicate the location of a sink
"""
def fixedFD(dem,sinks = None):
    # print ('flat2sink')
    sys.stdout.flush()
    sunk_dem = flat2sinks(dem)
    sunk_dem = sunk_dem.astype(np.float32)
    #np.save('sunk_dem.npy',sunk_dem)
    
    #d8 flow direction
    # print ('d8')
    sys.stdout.flush()
    fd = flow1_c.d8(sunk_dem)
    np.save('fd.npy', fd)

    #unique label all catchments to sinks/outlets (remember d8 sets all edge cells as sinks)
    # print ('catchments')
    sys.stdout.flush()
    catch = flow1_c.sink_catchments(fd)  #use sink catchments in flow1_c
    np.save('catch.npy', catch)

    #bit of a hack so all outflowing catchments flow to the same catchment make all
    #the edges flow to a catchment labeled 0
    #user defined sinks ie coastal nodata and legit sinks are also set to zero.
    catch[0] = 0
    catch[-1] = 0
    catch[:,0] = 0
    catch[:,-1] = 0
    if sinks is not None:
        catch[sinks != 0] = 0
        fd[sinks != 0] = 0
    # print ('saddles')
    sys.stdout.flush()
    saddles_dict = flow1_c.saddles(sunk_dem,catch)
    
    # print ('graph')
    sys.stdout.flush()
    #graph of catchments as nodes and saddles as edges
    G=nx.Graph()
    G.add_edges_from((u[0],u[1],d) for u,d in saddles_dict.items() )
    # nx.write_gpickle(G, 'graph.pkl')
    # with open('graph.pkl', 'wb') as f:
    #     pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    
    # print ('mst')
    sys.stdout.flush()
    mst=nx.minimum_spanning_tree(G,weight='elev')
    # nx.write_gpickle(mst, 'mst.pkl')
    # with open('mst.pkl', 'wb') as f:
    #     pickle.dump(mst, f, pickle.HIGHEST_PROTOCOL)

    # print ('predecessor')
    sys.stdout.flush()
    P = nx.predecessor(mst,0)
    fd_new = fd.copy()
    # print ('reversal')
    sys.stdout.flush()
    reversal = [reverse_fd(fd_new,mst[f][t[0]],f) for f,t in P.items() if len(t)==1]
    return sunk_dem, np.array(fd_new, dtype = np.int32)


if __name__ == '__main__':
    pass


