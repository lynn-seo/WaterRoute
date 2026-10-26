import os
import sys

from waterroute import utils
if sys.version_info < (3,8):
    raise ImportError("You are running on Python version below 3.8 which hasnt been tested.")

from scipy.ndimage import gaussian_filter
from waterroute.channel import Losses
import datetime
import requests
import re
import pathlib
import json
import pysheds  
import numpy as np
import time
import dijkstra3d

import pyproj

import ipyleaflet
from shapely.geometry import mapping,LineString
from shapely.geometry import shape
import geopandas as gpd
import pandas as pd
import pickle
import rasterio
from rasterio import features
from scipy.optimize import root_scalar
from scipy.signal import find_peaks
from scipy import interpolate
import matplotlib.pyplot as plt

from scipy import ndimage

#import livemaptiles
#from matplotlib import cm

import skfmm
from skimage.transform import resize
from skimage.graph import MCP,MCP_Geometric

from types import SimpleNamespace

# import fiona
#force geopandas to support kml
# fiona.supported_drivers['KML'] = 'rw'
# gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'


import site
# site.addsitedir(r'C:\Users\rea095\Documents\repos\ta_svn_various_tools\tiled_flow')
from tiled_flow import drain 

from waterroute.channel import earthworks


#probably dont need all the above imports need to cull

def trapezoidal_discharge(width, depth, S, n = 0.03, z=2):
    """
    Calculate discharge for a trapezoidal channel using Manning's equation.

    Parameters
    ----------
    width : float
        Bottom width (m)
    depth : float
        Flow depth (m)
    S : float
        Channel slope (m/m)
    n : float
        Manning roughness coefficient
    z : float, optional
        Side slope (horizontal : vertical). Default = 2

    Returns
    -------
    Q : float
        Discharge (m3/s)
    V : float
        Velocity (m/s)
    """

    # Cross-sectional area
    A = depth * (width + z * depth)

    # Wetted perimeter
    P = width + 2 * depth * np.sqrt(1 + z**2)

    # Hydraulic radius
    R = A / P

    # Velocity (Manning)
    V = (1 / n) * (R ** (2/3)) * (S ** 0.5)

    # Discharge
    Q = A * V

    return Q, V, A


def route_through_array(array, start, end, fully_connected=True,
                        geometric=True):
    """based on skimage.graph route_through_array
    modified to be faster"""
    start, end = tuple(start), tuple(end)
    if geometric:
        mcp_class = MCP_Geometric
    else:
        mcp_class = MCP
    m = mcp_class(array, fully_connected=fully_connected,sampling=(30,30))
    
    costs, traceback_array = m.find_costs([start], [end],find_all_ends=False)
    #print (costs[end]/1000000)
    return m.traceback(end), costs[end]

# #no long use
# def route_through_array_length_penalty_v2(array, start, end,
#                                        lenghth_penalty_per_km,
#                                        fully_connected=True,
#                                        geometric=True):
#     """
#     length_penalty : annualised cost per km
#     this uses buffer along the path
#     """
#     # print(f'length_penalty ($/cell)= {length_penalty}')
#     temp_array = array.copy()
#     paths = []
#     costs_list = []

#     for i in range(30):

#         m = mcp_class(temp_array, fully_connected=fully_connected, sampling=(30,30))
#         costs, _ = m.find_costs([start], [end], find_all_ends=False)

#         cost = costs[end]
#         if np.isinf(cost):
#             break

#         path = m.traceback(end)

#         paths.append(path)
#         costs_list.append(cost)

#         # apply penalty to get next route by buffering the path otherwise each of route is too same which is computationally inefficient 
#         for r, c in path:
#             r0 = max(0, r - 2)
#             r1 = min(temp_array.shape[0], r + 3)
#             c0 = max(0, c - 2)
#             c1 = min(temp_array.shape[1], c + 3)

#             temp_array[r0:r1, c0:c1] += length_penalty

#     # ✅ return cheapest
#     if paths:
#         min_idx = np.argmin(costs_list)
#         best_path = paths[min_idx]
#         best_cost = costs_list[min_idx]
#     else:
#         best_path, best_cost = None, np.inf

#     return best_path, best_cost


#this add linear cost on the earthwork cost array
def route_through_array_length_penalty(array, start, end,
                                    #    lenghth_penalty_per_km=None,
                                       fully_connected=True,
                                       geometric=True):
    """
    length_penalty : annualised cost per km 
    """
    
    start, end = tuple(start), tuple(end)

    
    # 🔹 Add length penalty to cost surface
    # penalised_array = array + length_penalty/buffer_size
    # length_penalty = 0
    # length_penalty = annualised_cost(577045/1000)
    
    if lenghth_penalty_per_km is None:
        lenghth_penalty_per_km = 0

    # print(f'length_penalty ($/cell)= {length_penalty}')
    # print(f'length_penalty ($/m)= {lenghth_penalty_per_km/1000}')
    # penalised_array = array + lenghth_penalty_per_km/1000
    penalised_array = array
    if geometric:
        mcp_class = MCP_Geometric
    else:
        mcp_class = MCP

    m = mcp_class(
        penalised_array,
        fully_connected=fully_connected,
        sampling=(30, 30)   # your grid spacing
    )

    costs, traceback_array = m.find_costs(
        [start],
        [end],
        find_all_ends=False
    )

    return m.traceback(end), costs[end]


_geod = pyproj.Geod(ellps="WGS84")


def geodeticLineLengthInKm(geom, km=False,total=True):
    """ an accurate calculation of length (in km) of a linestring
    assumes coordinates in longitude and latitude
    assumes geom is a shapely geom (or at least that numpy knows how to represent the coordinates in an array)
    for background reading see:
    https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid
    """
    # factor = 1
    # if km:
    #     factor = 1000
    # arr = np.asarray(geom)
    # segDist = _geod.inv(arr[:-1, 0], arr[:-1, 1], arr[1:, 0], arr[1:, 1])[2] / factor
    # if total:
    #     return segDist.sum()
    # return segDist
    if isinstance(geom, LineString):
        arr = np.array(geom.coords)  # Extract coordinates from LineString
    else:
        arr = np.asarray(geom)  # Ensure it's a NumPy array

    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("Input geometry must be a 2D array with shape (n, 2).")

    # factor = 1 if km else 1000
    factor = 1000 if km else 1
    segDist = _geod.inv(arr[:-1, 0], arr[:-1, 1], arr[1:, 0], arr[1:, 1])[2] / factor
    if total:
        return segDist.sum()
    return segDist



def crossDrainage (dem,aff,crs,path,thresholdAreaKm=1,catchDrainDepth=2, name=None):
    '''returns information on cross drainages
    
    thresholdAreaKm: only catchments above this threshold get a Q40
    
    '''
    #adjust the channel alignment relative to the terrain a small amount.
    #This simulates the redirection of small drainages to larger cross drainages where possible
    

    zzz = np.zeros_like(dem)
    zzz[path[:,0],path[:,1]] = catchDrainDepth
    insetDem = dem - zzz

    sinks = np.zeros_like(dem,dtype=int)
    sinks[path[:,0],path[:,1]] = 1

    #put in sinks at highest flow accumulation crossing points
    dem_sunk_original, fd_original = drain.fixedFD(dem, None)#FD based original dem 
    accum_original = drain.flow1_c.accumulation(fd_original)
    
    #find the cells in accum_original, greater than 1111
    #1111 is the threshold for a 10 cumec flow
    # accum_original_threshold = np.zeros_like(accum_original)
    # accum_original_threshold[accum_original>1111] = -2
    # dem3 = dem + accum_original_threshold
    # #min of two np arrays
    # modified_dem = np.minimum(dem, dem3)

    # dem_sunk, fd = drain.fixedFD(insetDem, None)
    # accum = drain.flow1_c.accumulation(fd)
 

    #Lynn's method
    crossings_overFlowThreshold = path[accum_original[path[:,0],path[:,1]]>1111]   
    #modified DEM based - method3, insetDEM +2 at crossings_overFlowThreshold
    modified_dem = insetDem.copy()
    for i in range(crossings_overFlowThreshold.shape[0]):
        modified_dem[crossings_overFlowThreshold[i,0],crossings_overFlowThreshold[i,1]] = modified_dem[crossings_overFlowThreshold[i,0],crossings_overFlowThreshold[i,1]] +2

    sunk_modified_dem, fd_modified_dem = drain.fixedFD(modified_dem, None)
    accum_modified_dem = drain.flow1_c.accumulation(fd_modified_dem)

    indexInPathToCrossings = find_peaks(accum_modified_dem[path[:,0],path[:,1]])[0]
    # indexInPathToCrossings = find_peaks(accum[path[:,0],path[:,1]])[0]
    
    crossings = path[indexInPathToCrossings]

    # sinkFD = fd.copy()
    sinkFD = fd_modified_dem.copy()
    sinkFD[crossings[:,0],crossings[:,1]] = 0
    sc = drain.flow1_c.sink_catchments(sinkFD)
    crossings_ids = sc[crossings[:,0],crossings[:,1]]
    crossing_geo_x,crossing_geo_y = aff * (crossings[:,1],crossings[:,0])
    crossingsDF = pd.DataFrame({'id': crossings_ids,'x':crossing_geo_x, 'y':crossing_geo_y,'indexInPathToCrossings':indexInPathToCrossings},index = crossings_ids)

    #calculate areas and centroids and polygons
    mask = np.isin(sc,sc[crossings[:,0],crossings[:,1]])
    # shapes = features.shapes(sc, mask=mask,connectivity=8,transform=aff,)
    # Fix: cast sc to int32 because rasterio doesn't support int64
    shapes = features.shapes(sc.astype(np.int32), mask=mask,connectivity=8,transform=aff,)
    crossingsGDF = gpd.GeoDataFrame(
        [
            {
                "id": int(ident),
                "geom": shape(poly),
                "centroid": shape(poly).centroid,
                "area_km": abs(_geod.geometry_area_perimeter(shape(poly))[0]) / 1_000_000,
            }
            for poly, ident in shapes
        ]
    )
    crossingsGDF.set_index('id',inplace=True)
    full_crossingsDF = crossingsGDF.join(crossingsDF)
    full_crossingsDF['q40'] = 0
    # Fix: Ensure q40 column is float to accept calculated float values
    full_crossingsDF['q40'] = full_crossingsDF['q40'].astype(float)

    subsetDF = full_crossingsDF[full_crossingsDF['area_km'] > thresholdAreaKm]
    import time
    # t_start = time.perf_counter()
    for r in subsetDF.itertuples():
        # q40 = q40_rffe(r.y,r.x,r.centroid.y,r.centroid.x,r.area_km)
        q40 = 16.3857 * r.area_km**0.8493 # nonlinear fit between catchment area and q40 from section, 
        full_crossingsDF.loc[r.Index,'q40'] = q40
        # print (f'extracting rffe for drainage {r.Index}. area {r.area_km}. q40 = {q40}')
    # print(f'Total RFFE time for {len(subsetDF)} catchments: {time.perf_counter()-t_start:.2f}s')

    return full_crossingsDF


def writeCrossDrainage (name,crossDF):
    gdf = gpd.GeoDataFrame(crossDF[['area_km','q40']], geometry=gpd.points_from_xy(x=crossDF.x, y=crossDF.y))
    #save gdf as pickle 
    gdf.to_pickle(f'{name}_cross_drainage.pkl')

    #only crossDF[['area_km','q40']] where q40>0
    # gdf_valid = gdf[gdf['q40'] > 0]
    gdf.crs = "EPSG:4326"
    gdf.to_file(f'{name}_cross_drainage.shp')
    gdf.to_file(f'{name}_cross_drainage.kml',driver='KML')
    # gdf_valid.to_file(f'{name}_cross_drainage_valid.shp')
    # gdf_valid.to_file(f'{name}_cross_drainage_valid.kml',driver='KML')

    crossDF[['area_km','q40']].to_csv(f'{name}_cross_drainage.csv')

    gdf = gpd.GeoDataFrame(crossDF[['geom','area_km','q40']], geometry='geom')
    # gdf_valid = gdf[gdf['q40'] > 0]
    gdf.crs = "EPSG:4326"
    gdf.to_file(f'{name}_cross_drainage_catchment.shp')
    gdf.to_file(f'{name}_cross_drainage_catchment.kml',driver='KML')
    # gdf_valid.to_file(f'{name}_cross_drainage_catchment_valid.shp')
    # gdf_valid.to_file(f'{name}_cross_drainage_catchment_valid.kml',driver='KML')


def q40_rffe(lato,lono,latc,lonc,areakm,savePath=None,region=0):
    '''extracts q40 flow from the arr rffe site, regional flood frequency estimation
    
    NOTE: ignores warning and errors that are sometimes reported with the model results web presentation.
    
    returns Q40 (2.5% AEP) Discharge (m3/s).
    
    Usually used through a user interface.
    Dont use heavily. It is unknown if the site is designed for heavy loads.
    You should manually check the results here https://rffe.arr-software.org/
    '''
    
    payload = {
        'catchment_name': 'Catchment1',
        'lato': f'{lato:.3f}',
        'lono': f'{lono:.3f}',
        'latc': f'{latc:.3f}',
        'lonc': f'{lonc:.3f}',
        'area': f'{areakm:.0f}',
        'region': f'{region:.0f}',
        'i2_6h': '0',
        'i50_6h': '0',
        }

    baseURL = 'https://rffe.arr-software.org/'
    response = requests.post(baseURL, data=payload)
    match = re.search(r'results = (\[{(.*?)}\]);', response.text)
    
    if not match:
        print(f"Warning: Could not find Q40 results for coordinates {latc}, {lonc}. The service may be unavailable or data may not exist for this location. Assuming Q40 is 0.")
        return 0
    
    result_string = re.search(r'results = (\[{(.*?)}\]);', response.text).group(1)
    arr_data = json.loads(result_string.replace("'", '"'))
    
    #get the Q20(=5 AEP) and Q50(2 AEP) (its not stored as %) and log interpolate the Q40 flow
    Q40 = np.interp(np.log10(40), np.log10([20,50]), (float(arr_data[3]['flow']),float(arr_data[4]['flow'])))
    # Q40 =  float(arr_data[5]['flow']) #1% AEP
    # Q40 =  float(arr_data[3]['flow']) #10% AEP
    if savePath:
        with open(savePath, "wb") as f:
            f.write(response.content)
    
    return Q40
    
def q40_rffe_old(lato,lono,latc,lonc,areakm,savePath=None):
    '''not used any more gets results directly from the first page

    extracts q40 flow from the arr rffe site
    
    returns Q40 (2.5% AEP) Discharge (m3/s).
    
    Usually used through a user interface.
    Dont use heavily. It is unknown if the site is designed for heavy loads.
    You should manually check the results here https://rffe.arr-software.org/
    '''

    payload = {
        'catchment_name': 'Catchment1',
        'lato': f'{lato:.3f}',
        'lono': f'{lono:.3f}',
        'latc': f'{latc:.3f}',
        'lonc': f'{lonc:.3f}',
        'area': f'{areakm:.0f}',
        'region': '0',
        'i2_6h': '0',
        'i50_6h': '0',
        }

    baseURL = 'https://rffe.arr-software.org/'
    downloadURL = 'https://rffe.arr-software.org/download-json'
    response = requests.post(baseURL, data=payload)

    docID = re.search(r'reef/out/json/(.*?)\.json', response.text).group(1)
    payload2 = {'fpath': f'/var/www/reef/reef/out/json/{docID}.json'}
    response2 = requests.post(downloadURL, data=payload2)
    arr_data = response2.json()['weighted_resp_base']

    #get the 20 and 50 AEP (its not stored as %) and log interpolate the Q40 flow
    Q40 = np.interp(np.log10(40), np.log10([20,50]), (arr_data[3]['flow'],arr_data[4]['flow']))
    
    if savePath:
        with open(savePath, "wb") as f:
            f.write(response.content)
    
    return Q40

def make_specific_cost_array(relDEM, bed_width,water_depth,freeboard,additional_per_m_cost):
    '''applies a cost of fill/cut function to a surface 
    
    relDEM is an elevation surface. Elevation are relative to the target natural surface gradeline
    ---
    rock: depth of rock layer
    '''

    results = []
    
    for dfl in range(-100,100):
        #negative dfl is excavation in getCostEarthworks but relDEM is positive in these locations
        works = earthworks.getCostEarthworks(DFL = -1 * dfl,
                                                    bedWidth = bed_width,
                                                    depthFlow = water_depth,
                                                    freeboard = freeboard,
                                                    crossGradient = 0,
                                                    sides = 2,
                                                    batterSides = 1.5,
                                                    bermSlope = 1/20,
                                                    batterHeight = 5,
                                                    bermWidth = 4,
                                                    bankWidth = 5,
                                                    overExcavateBed = 1,
                                                    overExcavateWall = 3,
                                                    channelLined = False, #To update this to pass initial parameter into this function
                                                    rockDepth = 5,)
        results.append({'dfl':dfl,'cost':(works['Total Cost'] )
                        +(additional_per_m_cost )})
        # results.append({'dfl':dfl,'cost':annualised_cost(works['Total Cost'] )
        #                 +annualised_cost(additional_per_m_cost )})

    df = pd.DataFrame(results)
    df.index = df.dfl
    #write pickle of df
    # with open('df.pkl', 'wb') as f:
    #      pickle.dump(df, f)
    # with open('realDEM.pkl', 'wb') as f:
    #      pickle.dump(relDEM, f)

    return np.interp(relDEM, df.dfl, df.cost)


def annualised_cost(cost, discount_rate = 0.07, lifespan= 50, QandM_percent=0.02):
        '''
        cost = capital cost
        discount_rate = discount rate (e.g. 0.05 for 5%)
        lifespan = lifespan of the asset in years
        QandM = ratio of capital cost for annual operation and maintenance
        '''
        return cost * (discount_rate * (1 + discount_rate) ** lifespan) / ((1 + discount_rate) ** lifespan - 1) + cost * QandM_percent


def make_specific_cost_array_withRockDepth(earthWorkCost_curve, relDEM,bed_width, water_depth, freeboard, additional_per_m_cost, rockDepth_raster,
                                           earthwork_cost_weight, linear_cost_factor):
    '''applies a cost of fill/cut function to a surface 
    
    relDEM is an elevation surface. Elevation are relative to the target natural surface gradeline
    rock: depth of rock layer
    '''
    
    # ...existing code...
    # import numpy as np

    # Expected df domain from your generation loop:
    # rock_depth: 0..49, dfl: -50..49
    ROCK_MIN, ROCK_MAX = 0, 49
    DFL_MIN, DFL_MAX = -50, 49

    # Build LUT once (very fast lookup)
    rock_idx = np.rint(earthWorkCost_curve["rock_depth"].to_numpy()).astype(np.int16)
    dfl_idx  = np.rint(earthWorkCost_curve["dfl"].to_numpy()).astype(np.int16)
    cost_val = earthWorkCost_curve["cost"].to_numpy(dtype=np.float64)

    cost_lut = np.full((ROCK_MAX - ROCK_MIN + 1, DFL_MAX - DFL_MIN + 1), np.nan, dtype=np.float64)
    cost_lut[rock_idx - ROCK_MIN, dfl_idx - DFL_MIN] = cost_val


    # Raster -> integer bins; leave invalid cells out of the cast.
    rock_finite = np.isfinite(rockDepth_raster)
    dfl_finite = np.isfinite(relDEM)
    rock_i = np.zeros(rockDepth_raster.shape, dtype=np.int16)
    dfl_i = np.zeros(relDEM.shape, dtype=np.int16)
    rock_i[rock_finite] = np.clip(
        np.rint(rockDepth_raster[rock_finite]), ROCK_MIN, ROCK_MAX
    ).astype(np.int16)
    dfl_i[dfl_finite] = np.clip(
        np.rint(relDEM[dfl_finite]), DFL_MIN, DFL_MAX
    ).astype(np.int16)
    valid = rock_finite & dfl_finite

    # Vectorized lookup
    specific_cost_array2 = np.full(relDEM.shape, np.nan, dtype=np.float64)
    specific_cost_array2[valid] = cost_lut[
        rock_i[valid] - ROCK_MIN,
        dfl_i[valid] - DFL_MIN
    ]

    # # print("Done:", specific_cost_array2.shape, "NaN fraction:", np.isnan(specific_cost_array).mean())
    # #write pickle of df
    # with open('df.pkl', 'wb') as f:
    #      pickle.dump(earthWorkCost_curve, f)
    # with open('realDEM.pkl', 'wb') as f:
    #      pickle.dump(relDEM, f)
    print(f"Earthwork cost weight: {earthwork_cost_weight}")

    return specific_cost_array2*earthwork_cost_weight




def makeCostArray(relDEM,slope = '20k_195cumec',additional_per_m_cost = 0):
    '''applies a cost of fill/cut function to a surface 
    This is original Bradfield that calculated from the excel spreadsheet by Kev, 
    This has been replaced by make_specific_cost_array
    
    relDEM is an elevation surface. Elevation are relative to the target natural surface gradeline
    
    test like this
    r = np.arange(20,-36,-1)
    c = makeCostArray(r)
    pd.DataFrame({'depth':r,'cost':c})
    '''
    
    depth = relDEM *-1
    if slope=='20k_195cumec':
    #Earth lined channel cost from Kev D 21/10/2020
    #modelled as 10k costed with channels at 20k
    #note: this is a 195cumec flow
        cutCost = 15.681*depth**2 - 396.14*depth + 2825.7 + additional_per_m_cost
        fillCost = 17.03*depth**2 + 637.00*depth - 1464.38 + additional_per_m_cost
    elif slope=='20k_195cumec_deepRock':    
        rockCutCost = 129.273*depth**2 - 839.2*depth + 4857.2 + additional_per_m_cost
        soilCutCost = 15.681*depth**2 - 396.14*depth + 2825.7 + additional_per_m_cost
        fillCost = 17.03*depth**2 + 637.00*depth - 1464.38 + additional_per_m_cost
        cutCost = np.where(depth<0,rockCutCost,soilCutCost)
    elif slope=='15k_156cumec':      
        #note: this is a 156cumec flow
        cutCost = 15.784*depth**2 - 349.52*depth + 2368.4 + additional_per_m_cost
        fillCost = 16.001*depth**2 + 625.53*depth - 1242.8 + additional_per_m_cost
    elif slope=='10k_156cumec':
        #note: this is a 156cumec flow
        cutCost = 15.051*depth**2 - 354.39*depth + 2116.4   + additional_per_m_cost    
        fillCost = 15.982*depth**2 + 612.13*depth - 875.95 + additional_per_m_cost
    elif slope=='5k_130cumec':
        #note: this is a 130cumec flow
        cutCost = 14.514*depth**2 - 344.96*depth + 1696.6 + additional_per_m_cost
        fillCost = 15.338*depth**2 + 613.45*depth - 485.03 + additional_per_m_cost
    elif slope=='10k_130cumec':
        #note: this is a 130cumec flow
        cutCost = 15.072*depth**2 - 321.51*depth + 1911.9  + additional_per_m_cost     
        fillCost = 15.766*depth**2 + 578.85*depth - 762.55 + additional_per_m_cost
    elif slope=='5k_156cumec':
        #note: this is a 156cumec flow
        cutCost = 15.394*depth**2 - 385.19*depth + 2047.9 + additional_per_m_cost
        fillCost = 15.041*depth**2 + 706.23*depth - 797.53 + additional_per_m_cost
    elif slope == 'pipe_195cumec':
        #79023 is the cost of 1m of pipe line on a 1:10000 slope carrying 195cumecs (it is pumping)
        cutCost = 79023 + (15619769*depth) #this isnt correct. I need to know the gradient that im going up ie or the elev change from the previous pixel
        fillCost = 79023
    return np.maximum(cutCost,fillCost)

def makeCostArray3D(dem,elevations):
    '''applies a cost of fill/cut function to a surface 
    
    relDEM is an elevation surface. Elevation are relative to the target natural surface gradeline
    
    test like this
    r = np.arange(20,-36,-1)
    c = makeCostArray(r)
    pd.DataFrame({'depth':r,'cost':c})
    '''
    costs = []
    for e in elevations:
        depth = (dem - e) *-1
        #Earth lined channel cost from Kev D 21/10/2020
        cutCost = 15.681*depth**2 - 396.14*depth + 2825.7
        fillCost = 17.03*depth**2 + 637.00*depth - 1464.38
        costs.append(np.maximum(cutCost,fillCost))
        
    return np.stack(costs)


def extrapolateSmoothFast(dem,path,cumDist,sampleF = 1000):
    #surface grid, slope grid and slope surface, 28 July 2025 commentedLynn
    grid_x, grid_y =np.mgrid[int(sampleF/2):dem.shape[0]:sampleF,int(sampleF/2):dem.shape[1]:sampleF]
    slopeSurface = interpolate.griddata(path[::100], cumDist[::100], (grid_x, grid_y), method='nearest')
    return resize(slopeSurface, dem.shape,order=1)

def rle(inarray):
        """ run length encoding. Partial credit to R rle function. 
            Multi datatype arrays catered for including non Numpy
            returns: tuple (runlengths, startpositions, values) """
        #from https://stackoverflow.com/a/32681075/604456
        ia = np.asarray(inarray)                # force numpy
        n = len(ia)
        if n == 0: 
            return (None, None, None)
        else:
            y = np.array(ia[1:] != ia[:-1])     # pairwise unequal (string safe)
            i = np.append(np.where(y), n - 1)   # must include last element posi
            z = np.diff(np.append(-1, i))       # run lengths
            p = np.cumsum(np.append(0, z))[:-1] # positions
            return(z, p, ia[i])
        
def add_control_structures_V3(profile):
    """
    Insert 'control_structure_id' where:
      - use_siphon_v2 == True
      - first row
      - last row
    """

    profile = profile.copy()
    profile['control_structure_id'] = np.nan

    # Identify rows that must contain a control structure
    control_mask = profile['use_siphon_v2'].fillna(False)

    # Ensure first and last rows included
    control_mask.iloc[0] = True
    control_mask.iloc[-1] = True

    # Assign sequential IDs
    profile.loc[control_mask, 'control_structure_id'] = \
        range(1, control_mask.sum() + 1)

    return profile


def add_control_structures_V2(profile,maxControlSpacing=28000,channelSlope=1/20000,minWaterDepth=1.5):
    '''
    split up the channel between siphons to have control structures at no further then 28km spacing.
    
    structures must occur at start and end of channel and at each siphon.
    each structure is given an ID and added to the profile.
    
    
    Can you suggest a generic ruleset for how control structure spacing is determined at various channel gradients and channel depths? 
    The maintaining 2m of water depth rule works in some cases but not all.  The 2m (and we used 3m in cases where this is possible) is about two things. 
      Firstly, that flow depth limits the growth of aquatic weeds to the banks only.  This is extremely important in any tropical or sub-tropical areas. 
        Secondly, it guarantees a minimum water depth for energy dissipation at part flow conditions at control points.  We just located them at points where the critical depth was reached.
            In the real world, you would select your control locations with a few other points in mind (such as doubling as access points, being close to power, avoiding high banks or fill,etc)
    '''
    #base all this on the lowest water depth in the channel to be conservative and avoid the complexity of varying it with length
    channelDepth = profile['water_depth'].min()
    controlSpacing=(channelDepth-minWaterDepth)/channelSlope
    if maxControlSpacing:
        controlSpacing=min(controlSpacing,maxControlSpacing)
    
    dist = np.insert(profile[profile['use_siphon']]['short_distance(m)'].array.copy(), 0, 0)
    dist = np.append(dist,profile['short_distance(m)'].iat[-1])
    NoControlStructures = (np.floor(np.diff(dist)/controlSpacing)+1).astype(int)
    controlChainage = np.cumsum(np.repeat(np.diff(dist)/(NoControlStructures),(NoControlStructures).astype(int)))
    controlChainage = np.insert(controlChainage, 0, 0)

    controlChainageDF = pd.DataFrame(controlChainage,columns=['short_distance(m)'],)
    controlChainageDF['control_structure_id'] = controlChainageDF.index + 1

    pc = profile.copy()
    pc['uid'] = pc.index
    control_join = pd.merge_asof(controlChainageDF,pc,on='short_distance(m)',direction='nearest')[['uid','control_structure_id']]
    outDF = profile.merge(control_join,how='left',right_on='uid',left_index=True)
    outDF.index = profile.index
    outDF.drop(columns=['uid'],inplace=True)
    return outDF

def add_control_structures(profile,q40SiphonThreshold=400,maxControlSpacing=28000,channelSlope=1/20000,channelDepth=4.4,minWaterDepth=3):
    '''
    split up the channel between siphons to have control structures at no further then 28km spacing.
    
    structures must occur at start and end of channel and at each siphon.
    each structure is given an ID and added to the profile.
    
    
    Can you suggest a generic ruleset for how control structure spacing is determined at various channel gradients and channel depths? 
    The maintaining 2m of water depth rule works in some cases but not all.  The 2m (and we used 3m in cases where this is possible) is about two things. 
      Firstly, that flow depth limits the growth of aquatic weeds to the banks only.  This is extremely important in any tropical or sub-tropical areas.  
      Secondly, it guarantees a minimum water depth for energy dissipation at part flow conditions at control points.  We just located them at points where the critical depth was reached.
          In the real world, you would select your control locations with a few other points in mind (such as doubling as access points, being close to power, avoiding high banks or fill,etc)
    '''
    controlSpacing=(channelDepth-minWaterDepth)/channelSlope
    if maxControlSpacing:
        controlSpacing=min(controlSpacing,maxControlSpacing)
    
    dist = np.insert(profile[profile['q40']>q40SiphonThreshold]['short_distance(m)'].array.copy(), 0, 0)
    dist = np.append(dist,profile['short_distance(m)'].iat[-1])
    NoControlStructures = (np.floor(np.diff(dist)/controlSpacing)+1).astype(int)
    controlChainage = np.cumsum(np.repeat(np.diff(dist)/(NoControlStructures),(NoControlStructures).astype(int)))
    controlChainage = np.insert(controlChainage, 0, 0)

    controlChainageDF = pd.DataFrame(controlChainage,columns=['short_distance(m)'],)
    controlChainageDF['control_structure_id'] = controlChainageDF.index + 1

    pc = profile.copy()
    pc['uid'] = pc.index
    control_join = pd.merge_asof(controlChainageDF,pc,on='short_distance(m)',direction='nearest')[['uid','control_structure_id']]
    outDF = profile.merge(control_join,how='left',right_on='uid',left_index=True)
    outDF.index = profile.index
    outDF.drop(columns=['uid'],inplace=True)
    return outDF


def get_min_max_indices(indices, sunwater_raster):
    '''
    this is to get cooridnate indices based on sunwater alginment raster
    this shouldn't be used for oridinary modelling
    '''
    import numpy as np
    max_x = np.max(indices[:,0])
    min_x = np.min(indices[:,0])
    max_y = np.max(indices[:,1])
    min_y = np.min(indices[:,1])
    #get the location of where sunwaterPath = 0
    zero_locations = np.where(sunwater_raster == 0)
    # print("Locations where sunwaterPath is 0:", zero_locations)
    #generate array of x by 2 with the x and y coordinates of the zero locations
    zero_coordinates = np.column_stack(zero_locations)
    #constraint zero_coordinates[:,0] to be between min_x and max_x, and zero_coordinates[:,1] to be between min_y and max_y
    zero_coordinates = zero_coordinates[(zero_coordinates[:,0] >= min_x) & (zero_coordinates[:,0] <= max_x) & (zero_coordinates[:,1] >= min_y) & (zero_coordinates[:,1] <= max_y)]
    # print("Zero coordinates within bounds:", zero_coordinates)
    resized_sunwater_indices = []
    for i in range(len(indices)):
        distances = np.sqrt((zero_coordinates[:,0] - indices[i,0])**2 + (zero_coordinates[:,1] - indices[i,1])**2)
        closest_index = np.argmin(distances)
        resized_sunwater_indices.append(zero_coordinates[closest_index])
    indices = np.array(resized_sunwater_indices)
    return indices


def calc_path(output_path, weight, target_elev, startCell,finishCell,aff, loss_value_dollar_per_ML, total_linear_cost, total_loss_ML,
              withLengthPenalty = True, is_sunwater_raster = False, approximate_length = 1,
              discount_rate = 0.07, lifespan = 50, linear_cost_factor = 1.0):
    
    # indices, totalWeight = route_through_array(weight, startCell, finishCell,geometric=True)
    # farmCrossing_cost = 4081275
    # easement_lining_cost = 9714603
    # catchDrainage_cost = 792007
    # SCADA_cost = 224402
    
    if total_linear_cost is None:
        lenghth_penalty = 0
    else:
        total_annualised_linearCost_per_km = annualised_cost(
            total_linear_cost / approximate_length,
            discount_rate=discount_rate,
            lifespan=lifespan,
        ) #this is per km value not per cell
        #total_loss_ML is annual value; no need to annualise it. 
        lenghth_penalty_per_km = total_annualised_linearCost_per_km + total_loss_ML*loss_value_dollar_per_ML/approximate_length
    
    # legnth_penalty = 0
    # 7487, dollar per cell, from linear cost
    # 12579, is value of loss per cell 
    # print(f'length_penalty ($/cell)= {legnth_penalty}')
    # print('approximated length (km) =',approximate_length)
    # print('linear cost per km =',lenghth_penalty_per_km)
    weight = weight + lenghth_penalty_per_km/1000*linear_cost_factor
    # penalised_array = array + lenghth_penalty_per_km/1000
    if (withLengthPenalty):
        indices, totalWeight = route_through_array_length_penalty(weight, startCell, finishCell,geometric=True)
    else:
        indices, totalWeight = route_through_array(weight, startCell, finishCell,geometric=True)
    
    #this overwrite the existing path indices
    if is_sunwater_raster:
        sunwaterPath = r"\\fs1-cbr.nexus.csiro.au\{lw-bradfield}\work\3_New_Bradfield_2020\1_Project_management\8_Special_issue\2_WaterRoute\3_Model\0_DEM\rasterize_sunwater_route.tif"
        left = 147
        top = -19
        right = 149
        bottom = -21
        routing_bounds = (left, bottom, right, top)
        sunwater_raster, _, _ = utils.load_raster(sunwaterPath,routing_bounds,1,1,utils.Resampling.average)
        path = np.array(indices)
        path = get_min_max_indices(path, sunwater_raster) 
    else:
        path = np.array(indices)

    
    # l = asLineString(np.column_stack(aff * (path[:,1],path[:,0])))
    l = LineString(np.column_stack(aff * (path[:,1],path[:,0])))
    lDist = np.insert(geodeticLineLengthInKm(l,False,False), 0, 0)
    cumDist = np.cumsum(lDist)
    l = l.simplify(0.0005, preserve_topology=False)
    return (path,l,lDist,cumDist)

def path_outputs(dem,aff,crs,name,weight,target_elev,targetSlope,path,l,lDist,cumDist,water_depth,max_fill_siphon_ends = 5,doCrossDrainage=False, rockDepth = None):
    leng = geodeticLineLengthInKm(l,False)
    pathDF = gpd.GeoDataFrame([l], columns=['geometry'])
    path_layer = ipyleaflet.GeoData(geo_dataframe=pathDF, style={'color':'red', 'weight':2})
    path_layer.name = str(target_elev)
    #removeLayer(m,path_layer.name)
    #m.add_layer(path_layer)
    pathDF.crs = "EPSG:4326"

    pathDF.to_file(f'{name}.shp')
    pathDF.to_file(f'{name}.kml',driver='KML')
    
    pathElev = dem[path[:,0],path[:,1]]
    #adjust to scale to the slightly short route length of the simplified line
    # unit_earthWorkCost = weight[path[:,0],path[:,1]]
    
    pathCost = lDist * weight[path[:,0],path[:,1]]*(leng/cumDist[-1])
    rockDepth = rockDepth[path[:,0],path[:,1]]
    adjustedDist = cumDist*(leng/cumDist[-1])
    
    cumPathCost = np.cumsum(pathCost)
    # print('targetSlope =', targetSlope)
    targetElev = target_elev + adjustedDist*targetSlope
    coords = (aff * (path[:,1],path[:,0]))
    profile = pd.DataFrame({'x':coords[0],
                            'y':coords[1],
                            #'long_distance(m)':cumDist,
                            'short_distance(m)':adjustedDist,
                            'flow_level_elevation(m)':targetElev,
                            'natural_surface_elevation(m)':pathElev,
                            'smoothedSurface(m)' : smooth(pathElev),
                            'diff_flowLevel_naturalSurface(m)':pathElev-targetElev,
                            'rock_depth(m)':rockDepth,
                            'routing cost (2020 AUD)':pathCost,
                            'routing cost to this point (2020 AUD)':cumPathCost,
                            
                           })
    
    # smoothedSurface = pd.DataFrame({'x':coords[0],
    #                         'y':coords[1],
    #                         #'long_distance(m)':cumDist,
    #                         'short_distance(m)':adjustedDist,
    #                         'flow_level_elevation(m)':targetElev,
    #                         'natural_surface_elevation(m)':pathElev,
    #                         'diff_flowLevel_naturalSurface(m)':pathElev-targetElev,                                               'smoothedSurface(m)' : pathElev
    #                         'routing cost (2020 AUD)':pathCost,
    #                         'routing cost to this point (2020 AUD)':cumPathCost,
    #                        })
    #length of channel through each cell
    s = profile['short_distance(m)']
    profile['cell_distance(m)'] = (s.rolling(3,center=True,min_periods=1).max() - s.rolling(3,center=True,min_periods=1).min())/2
    profile['gradeline_to_NS_(m)'] = profile['natural_surface_elevation(m)'] - profile['flow_level_elevation(m)'] 
    
    crossDF = None
    if doCrossDrainage:
        #add crossings information.
        # print(class(dem))
        crossDF = crossDrainage (dem,aff,crs,path,thresholdAreaKm=1, catchDrainDepth=2, name = name)
        profileWithCrossings = profile.merge(crossDF[['area_km','id','indexInPathToCrossings','q40']],how='left',left_index=True,right_on='indexInPathToCrossings')
        profileWithCrossings.index = profile.index
        profileWithCrossings.drop(columns=['indexInPathToCrossings'],inplace=True)
        profile = profileWithCrossings
        #profile = add_control_structures(profileWithCrossings,q40SiphonThreshold=400,maxControlSpacing=None,channelSlope=0.00005,
        # channelDepth=4.4,minWaterDepth=3)
        
        #work out the length of sections of channel below the threshold where fill can be used prior to a siphon
        e = profile['natural_surface_elevation(m)'] - (profile['flow_level_elevation(m)'] - water_depth - max_fill_siphon_ends)
        l,pos,v = rle(e<0)
        segLength = np.repeat(np.diff(np.append(profile['short_distance(m)'][pos],profile['short_distance(m)'].iat[-1])),l)
        #if NS doesnt drop low enough use a default siphon length.
        #q40 / 2 is a rule of thumb i developed and tested with the following code.
        #bed_width = max(80,Losses.channel_width_for_discharge(q40,flow_depth=2,depth_scale = None, S= 1/2000, n=0.035))
        #could just run the above code I guess.
        #My logic behind this is: 
        #a)	If the natural channel is much wider than this its going to be slow sluggish flood flows and its seems to be
        #  common practice for roads to be built up on such flood plain forcing the flow into a narrower but still quite wide flow path. 
        # Seems reasonable that a channel can do the same.
        #b)	If the natural channel in much narrower than this then there is a good chance its more incised and
        #  it passes the 5m of fill criteria and the siphon length will be determined from the terrain
        
        # bed_width = max(80,Losses.channel_width_for_discharge(q40,flow_depth=2,depth_scale = None, S= 1/1000, n=0.025))
        # segLength[e>=0]=profile[e>=0]['q40']/2
        
        
        #for loop for profile[e >= 0]['q40'] only noNan
        # subset_mask = e >= 0
        # valid_q40_mask = profile['q40'].notna() & (profile['q40'] > 0)
        
        # # Combine with the subset mask
        # update_mask = subset_mask & valid_q40_mask

        # # Calculate widths only for the valid rows
        # if update_mask.any():
        #     widths = [Losses.channel_width_for_discharge(d, max_flow_depth=2, depth_scale=None, S=1/1000, n=0.025)
        #               for d in profile.loc[update_mask, 'q40']]
            
        #     # Update segLength only at the valid locations
        #     segLength[update_mask] = widths

        segLength[e>=0] = 19.606*profile[e>=0]['area_km']**0.2995

        #hard limit of siphons must be over 80m long
        segLength = segLength.clip(min = 80) 

        profile['length_of_siphon(m)'] = segLength

    profile.to_csv(f'{name}+{pathCost.sum()/1000000:.2f}.csv')
    #these bits are useful but slow sometimes randomly over the network
    #gdf = gpd.GeoDataFrame(profile, geometry=gpd.points_from_xy(profile['x'], profile['y']))
    #gdf.to_file(f'{name}_profile_points.shp')
    #gdf.to_file(f'{name}_profile_points.kml',driver='KML')
    
    print(f'length of {target_elev} is {leng/1000:.2f}km. Cost is (2020 AUD) M${pathCost.sum()/1000000:.2f}')
    return (profile,crossDF)


def smooth(y, box_pts=30):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


# def elevRangePathsWithJunction(output_path,dem,aff,fr,to,targetSlope,start,finish,max_fill_siphon_ends, contour,doCrossDrainage,sampleF,
#                    bed_width,water_depth,freeboard,additional_per_m_cost, channel_capacity, rockDepth_raster, 
#                    loss_value_dollar_per_ML, total_linear_cost, sunWater_route):
#     #junction point considered where the slope changes
#     paths = {}
#     for elev in range(fr,to,1):
#     #for endElev in range(fr, to, 1):
#         print (f'target elevation is {elev}')
#         path = elevPath(output_path,dem,aff,targetSlope,start,finish,max_fill_siphon_ends,elev, contour,doCrossDrainage,sampleF,
#                         bed_width,water_depth,freeboard,additional_per_m_cost, channel_capacity, rockDepth_raster, 
#                         loss_value_dollar_per_ML, total_linear_cost, sunWater_route)
#         paths[elev] = path
#     return paths


def elevRangePaths(output_path,dem,aff,crs,fr,to,targetSlope,start,finish,max_fill_siphon_ends, contour,doCrossDrainage,sampleF,
                   bed_width,water_depth,freeboard,additional_per_m_cost, channel_capacity, rockDepth_raster, 
                   loss_value_dollar_per_ML, total_linear_cost, withLengthPenalty, cost_array_with_rock_depth, 
                   sunWater_route,earthwork_cost_weight, approximate_length, discount_rate=0.07, lifespan=50,
                   linear_cost_factor=1.0):
    paths = {}
    for elev in range(fr,to,1):
    # for i in range(0,1): #modified by Lynn to handle float start_elevation
        # elev = fr + i
        # print (f'target elevation is {elev}')
        path = elevPath(output_path,dem,aff,crs,targetSlope,start,finish,max_fill_siphon_ends,elev, contour,doCrossDrainage,sampleF,
                        bed_width,water_depth,freeboard,additional_per_m_cost, channel_capacity, rockDepth_raster, 
                        loss_value_dollar_per_ML, total_linear_cost,  withLengthPenalty, cost_array_with_rock_depth,
                        sunWater_route,earthwork_cost_weight, approximate_length, discount_rate, lifespan,
                        linear_cost_factor)
        paths[elev] = path
    return paths

        
def elevPath(output_path,dem,aff,crs,targetSlope,start,finish,max_fill_siphon_ends,target_elev,
             contour,doCrossDrainage,sampleF,bed_width,water_depth,freeboard,additional_per_m_cost, channel_capacity, rockDepth_raster, 
             loss_value_dollar_per_ML, total_linear_cost, withLengthPenalty, cost_array_with_rock_depth, sunWater_route,earthwork_cost_weight,
               approximate_length, discount_rate=0.07, lifespan=50, linear_cost_factor=1):
    
    # startCell = (~aff * start)
    # finishCell = (~aff * finish)
    col_startCell, row_startCell = (~aff * start)
    col_finishCell, row_finishCell  = (~aff * finish)
    # startCell = (np.floor(startCell[1]).astype(int),np.floor(startCell[0]).astype(int))
    # finishCell = (np.floor(finishCell[1]).astype(int),np.floor(finishCell[0]).astype(int))
    startCell = (np.floor(row_startCell).astype(int),np.floor(col_startCell).astype(int))
    finishCell = (np.floor(row_finishCell).astype(int),np.floor(col_finishCell).astype(int))
    
    relDEM = dem - contour


    results = []
    for rock_depth in range(0, 50, 1):
        for dfl in range(-50, 50):

            works = earthworks.getCostEarthworks(
                DFL = - 1 * dfl,  ####IMPORTANT: should put negative sign : DFL is negative of the depth to flow level, this is because the cost curves are based on depth of cut or fill which is the negative of the DFL
                bedWidth = bed_width,
                depthFlow = water_depth,
                freeboard = freeboard,
                crossGradient = 0,
                sides = 2,
                batterSides = 1.5,
                bermSlope = 1/20,
                batterHeight = 5,
                bermWidth = 4,
                bankWidth = 5,
                overExcavateBed = 1,
                overExcavateWall = 3,
                channelLined = False,
                rockDepth = rock_depth
            )
            results.append({
                'rock_depth': rock_depth,
                'dfl': dfl,
                'cost': annualised_cost(
                    works['Total Cost'],
                    discount_rate=discount_rate,
                    lifespan=lifespan,
                ) + annualised_cost(
                    additional_per_m_cost,
                    discount_rate=discount_rate,
                    lifespan=lifespan,
                ),
                'channel_lined': False
            })

    earthWorkCost_curve = pd.DataFrame(results)
    #index by dfl and rock depth for interpolation
    earthWorkCost_curve.index = earthWorkCost_curve.apply(lambda row: f"{row.rock_depth}_{row.dfl}", axis=1)

    # weight_original = make_specific_cost_array(dem - contour,bed_width,water_depth,freeboard,additional_per_m_cost)
    if cost_array_with_rock_depth:
        weight = make_specific_cost_array_withRockDepth(earthWorkCost_curve, dem - contour,bed_width,water_depth,freeboard,additional_per_m_cost, 
                                                        rockDepth_raster, earthwork_cost_weight, linear_cost_factor)
    else:
        weight = make_specific_cost_array(dem - contour,bed_width,water_depth,freeboard,additional_per_m_cost)
        
    # weight = make_specific_cost_array_withRockDepth(earthWorkCost_curve, dem - contour,bed_width,water_depth,freeboard,additional_per_m_cost, rockDepth_raster)
    
    # with open(f'{output_path}_{target_elev}_weight_original.pkl', 'wb') as f:
    #         pickle.dump(weight_original, f)
    # with open(f'{output_path}_{target_elev}_weight_with_rock_depth.pkl', 'wb') as f:
    #         pickle.dump(weight, f)
    
    properties_earthworks = earthworks.getCostEarthworks(DFL = 0,bedWidth = bed_width, depthFlow = water_depth, freeboard = freeboard,channelLined = False , rockDepth = 5)
    water_section_area = properties_earthworks['Geometries']['waterSection'].area
    wetted_perimeter = properties_earthworks['Geometries']['wettedPerimeterSection'].length
    water_surface_width = properties_earthworks['Geometries']['waterSection'].bounds[2] - properties_earthworks['Geometries']['waterSection'].bounds[0]

    annual_loss_ML_per_annum = Losses.getAnnualLossesLinear(annualDuration = 330, channelLength = 1, surfaceWidth = water_surface_width, evapPanAnnual = 2,
                        evapPeakDaily = 6, bedPerimeter = bed_width, k = 0.002, 
                        lossEachStockDomestic = 2, spacingStockDomestic = 2, leakageRate = 0, channelXSection = water_section_area,rolling_shutdown_factor = 1, 
                        waterDepth= water_depth)
    
    pd.Series(annual_loss_ML_per_annum, name='ML_per_annum').to_csv(
        f'{output_path}_{target_elev}_annual_loss.csv',
        header=True
    )

    total_loss_ML = annual_loss_ML_per_annum['Total'] 

    # #addLayer(m,'contour_weight',1/weight,aff,crs,1,10,cm.viridis)
    
    path,l,lDist,cumDist = calc_path(output_path, weight,target_elev,startCell,finishCell,aff, loss_value_dollar_per_ML,
                                    total_linear_cost=total_linear_cost, total_loss_ML=total_loss_ML, 
                                    withLengthPenalty=False, is_sunwater_raster = False,approximate_length=approximate_length,
                                    discount_rate=discount_rate, lifespan=lifespan,
                                    linear_cost_factor=linear_cost_factor)
    
    # with open(f'{output_path}_{target_elev}_weigth0.pkl', 'wb') as f:
    #         pickle.dump(weight, f)
        
    # path,l,cumDist,totalCost,profile,crossDF = calc_path(f'{output_path}_{target_elev}_contour_path',weight,target_elev,startCell,finishCell)
    leng = geodeticLineLengthInKm(l,False)
    # with open(f'{output_path}_{target_elev}_leng.pkl', 'wb') as f:
    #     pickle.dump(leng, f)
    

    for iterN in range(1,5):
        print (f'iteration number: {iterN}')
        if targetSlope is None:
            #estimate a slope
            # print(f'iterN {iterN} leng {leng}')
            # est_targetSlope = (dem[finishCell] - dem[startCell])/leng
            # print(f'target slope is {est_targetSlope}')
            adj_targetSlope = (dem[finishCell] - dem[startCell])/cumDist[-1]
            #always negative adj_targetSlope, Lynn edited
            adj_targetSlope = -(abs(adj_targetSlope))
        else:
            adj_targetSlope = leng/cumDist[-1]*targetSlope
        #print (f'target slope = {adj_targetSlope}')
        
        slopeSurface = extrapolateSmoothFast(dem,path,cumDist*adj_targetSlope,sampleF = sampleF)
        flattened = dem - target_elev - slopeSurface

        if cost_array_with_rock_depth:
            weight = make_specific_cost_array_withRockDepth(earthWorkCost_curve, flattened,bed_width,water_depth,freeboard,additional_per_m_cost, 
                                                            rockDepth_raster, earthwork_cost_weight, linear_cost_factor)
        else:
            weight = make_specific_cost_array(flattened,bed_width,water_depth,freeboard,additional_per_m_cost)
                      
        # with open(f'{output_path}_{target_elev}_weigth{iterN}.pkl', 'wb') as f:
        #     pickle.dump(weight, f)
        
        #addLayer(m,f'weight',weight,aff,crs,1,50000,cm.viridis)
        
        if sunWater_route and iterN == 4: #not used anymore
            path,l,lDist,cumDist = calc_path(output_path, weight,target_elev,startCell,finishCell,aff, loss_value_dollar_per_ML,
                                          total_linear_cost=total_linear_cost, total_loss_ML=total_loss_ML, 
                                        withLengthPenalty=withLengthPenalty, is_sunwater_raster = True,approximate_length=approximate_length,
                                        discount_rate=discount_rate, lifespan=lifespan,
                                        linear_cost_factor=linear_cost_factor)
        else:
            if iterN == 4:
                path,l,lDist,cumDist = calc_path(output_path, weight,target_elev,startCell,finishCell,aff, loss_value_dollar_per_ML,
                                          total_linear_cost=total_linear_cost, total_loss_ML=total_loss_ML,
                                            withLengthPenalty=withLengthPenalty, is_sunwater_raster = False, approximate_length=approximate_length,
                                            discount_rate=discount_rate, lifespan=lifespan,
                                            linear_cost_factor=linear_cost_factor)
            else:                
                path,l,lDist,cumDist = calc_path(output_path, weight,target_elev,startCell,finishCell,aff, loss_value_dollar_per_ML,
                                          total_linear_cost=total_linear_cost, total_loss_ML=total_loss_ML,
                                            withLengthPenalty=withLengthPenalty, is_sunwater_raster = False, approximate_length=approximate_length,
                                            discount_rate=discount_rate, lifespan=lifespan,
                                            linear_cost_factor=linear_cost_factor)
            
        
        leng = geodeticLineLengthInKm(l,False)
        # print(f'iterN {iterN} after adjusting leng {leng}')
        #pathDF = gpd.GeoDataFrame([l], columns=['geometry'])


    # with open(f'{output_path}_{target_elev}_path.pkl', 'wb') as f:
    #     pickle.dump(path, f)
    # 
    # #when targetSlope is None    
    targetSlope = adj_targetSlope
    #after last iteration
    # water_depth = max_flow_depth
    #pickle rockDepth_raster
    # with open(f'{output_path}_{target_elev}_rockDepth_raster.pkl', 'wb') as f:
    #         pickle.dump(rockDepth_raster, f)
    
    

    profile,crossDF = path_outputs(dem,aff,crs,f'{output_path}_{target_elev}',weight, target_elev,targetSlope,path,l,lDist,cumDist,water_depth,
                                   max_fill_siphon_ends = 5,doCrossDrainage=doCrossDrainage, rockDepth = rockDepth_raster)
    
    if crossDF is not None:
        writeCrossDrainage (f'{output_path}_{target_elev}',crossDF)
        crossDF.to_pickle(f'{output_path}_{target_elev}_crossDF.pkl')
    
    profile.to_pickle(f'{output_path}_{target_elev}_profile.pkl')
    
    return (path,profile,l)
