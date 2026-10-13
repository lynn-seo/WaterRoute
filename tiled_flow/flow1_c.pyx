'''
Flowgrid implements a generator for traversing a grid of single (D8) flow directions
Two flavours
 - one that traverses the flow network in an oriented manner starting at outflows across edges and reports
   each cell as it is encountered; used for stream function mapping
 - one that starts at a local top (no inflow) and yields cells in order so that no cell is yielded before
   all its inflows are resolved (all upslope cells have been processed)
Input is a numpy array of flow directions using the ArcInfo convention:
    32 64 128
    16     1
     8  4  2

'''
#you build this by using python setup.py build_ext --inplace

#NB numpy arrays use [row,col] indexing which is the opposite way around to (x,y)
#and [0,0] is top left corner

import numpy as np
import math

cimport numpy as np
DTYPE = np.float32
#ctypedef np.float32_t DTYPE_t
#DTYPE_int = np.int_
#ctypedef np.int_t DTYPE_int_t
#DTYPE_int32 = np.int32
#ctypedef np.int32_t DTYPE_int32_t

# floating point
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

# 32-bit integer (e.g. flow direction grids)
DTYPE_int32 = np.int32
ctypedef np.int32_t DTYPE_int32_t

# index / label / catchment IDs (portable across OS)
DTYPE_int = np.int_
ctypedef np.npy_intp DTYPE_int_t

cdef float NodataValue = -9999.0

coordOffsets = {1:(0,1), 2:(1,1), 4:(1,0), 8:(1,-1), 16:(0,-1), 32:(-1,-1), 64:(-1,0), 128:(-1,1)}
coordOffsets_dir = {(0,1):1, (1,1):2, (1,0):4, (1,-1):8, (0,-1):16, (-1,-1):32, (-1,0):64, (-1,1):128}
coordOffsets_list = [(0,1), (0,-1), (1,0), (-1,0), (1,1), (1,-1), (-1,-1), (-1,1)] #order can be important
half_coordOffsets_list = [(0,1), (1,1), (1,0), (1,-1)]
coordOffsets_flows_to = {(0,1):16, (1,1):32, (1,0):64, (1,-1):128, (0,-1):1, (-1,-1):2, (-1,0):4, (-1,1):8}

def InGrid(np.ndarray[DTYPE_t,ndim=2] grid, int row, int col):
    return row >= 0 and row < grid.shape[0] and col >= 0 and col < grid.shape[1]

def NoData(np.ndarray[DTYPE_t,ndim=2] grid, int row, int col):
    return not InGrid(grid,row,col) or grid[row,col] == NodataValue

## if there is no neighbour lower its a sink - includes flat areas - value of 0
##flows towards nodata neighbours
##nodata as sinks - value of 0
##edge cells  as sinks - value of 0
def d8 (np.ndarray[DTYPE_t,ndim=2] dem not None):
    cdef int row_len = dem.shape[0]
    cdef int col_len = dem.shape[1]
    cdef np.ndarray fd = np.zeros([row_len,col_len], dtype=DTYPE_int32)
    cdef int row,col,n_row,n_col
    cdef int row_off,col_off
    cdef float steepest_slope,slope,sqrt2
    sqrt2 =  math.sqrt(2)

    #skip on the edge cells - they can just stay set to 0
    for row in range(1,row_len-1):
        for col in range(1,col_len-1):
            steepest_slope = 0
            #no flow default #its already all zeros.
#            if dem[row,col] == NodataValue:
#                fd[row,col] = 0
#            else:
            if dem[row,col] != NodataValue:
                for row_off,col_off in coordOffsets_list:
                    n_row = row+row_off
                    n_col = col+col_off
                    #if NoData(dem,n_row,n_col): #works but is a bit slower due to continual access to grid.shape
                    #if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len):  #just not looking at the edge cells to avoid this check
                    if dem[n_row,n_col] == NodataValue:
                        fd[row,col] = coordOffsets_dir[(row_off,col_off)]
                        break
                    else:
                        slope = dem[row,col] - dem[n_row,n_col]
                        if 0 not in [row_off,col_off]:
                            slope = slope / sqrt2
                        if slope > steepest_slope:
                            steepest_slope = slope
                            fd[row,col] = coordOffsets_dir[(row_off,col_off)]
    return fd


#assumes single flow direction
def catchments(np.ndarray[DTYPE_int32_t,ndim=2] fd not None,np.ndarray[DTYPE_int_t,ndim=2] catch_ends not None):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef np.ndarray catch = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int l_row,l_col,row,col,n_row,n_col
    cdef int row_off,col_off

    todo = []
    for row in range(row_len):
        for col in range(col_len):
            if catch_ends[row,col] > 0:
                catch[row,col] = catch_ends[row,col]
                todo.append((row,col))
            while todo:
                l_row,l_col = todo.pop()
                #search upstream
                for row_off,col_off in coordOffsets_list:
                    n_row = l_row+row_off
                    n_col = l_col+col_off
                    if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len) or fd[n_row,n_col] == NodataValue:
                        pass #not interested in no data or off edges data??
                    elif fd[n_row,n_col] == coordOffsets_flows_to[row_off,col_off] and catch_ends[n_row,n_col] == 0:
                        catch[n_row,n_col] = catch[l_row,l_col]
                        todo.append((n_row,n_col))
    return catch

#assumes single flow direction
#all edge pixels flow direction needs to be zero
#todo could try a variation where it doesnt keep a list but just tracks its last step.
def accumulation(np.ndarray[DTYPE_int32_t,ndim=2] fd not None):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef np.ndarray accum = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int l_row,l_col,row,col,n_row,n_col
    cdef int row_off,col_off
    cdef int dir_i,init_dir_i
#    cdef int flag  #used as a bool but simpler then declaring a bool in cython
#    accum = accum - 1
    init_dir_i = 0
    todo = []
    #skip on the edge cells - they can just stay set to 0
    for row in range(0,row_len):
        for col in range(0,col_len):
#            if accum[row,col] == -1: #not processed yet
            if fd[row,col] == 0:
                todo.append([row,col,init_dir_i])
            while todo:
                #search upstream
                l_row,l_col,dir_i = todo[-1]
#                l_row,l_col,dir_i = todo.pop()
#                flag = 0
#                if dir_i == init_dir_i: #first look at cell
#                    if accum[l_row,l_col] == -1:
#                        accum[l_row,l_col] = 0
#                    else:
#                        #already processed
#                        flag = 1

#                if flag == 0 and dir_i < 8:
                if dir_i < 8:
                    #put it back on the list
                    todo[-1][2] = dir_i+1
#                    todo.append((l_row,l_col,dir_i+1))
                    #look at the neighbouring cell
                    row_off,col_off = coordOffsets_list[dir_i]
                    n_row = l_row+row_off
                    n_col = l_col+col_off
                    #if it flows here queue it up
                    if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len) or fd[n_row,n_col] == NodataValue:
                        #print n_row,n_col,dir_i,l_row,l_col,fd[l_row,l_col]
                        pass
                    else:
                        if fd[n_row,n_col] == coordOffsets_flows_to[row_off,col_off]:
                            todo.append([n_row,n_col,init_dir_i])
                else: # all neighbours have been processed - deal with this cell
                    #write my accumulation downstream
                    todo.pop()
                    if fd[l_row,l_col] != 0:
                        row_off,col_off = coordOffsets[fd[l_row,l_col]]
                        n_row = l_row+row_off
                        n_col = l_col+col_off
#                        if accum[n_row,n_col] >= 0:
#                            accum[n_row,n_col] = accum[n_row,n_col] + accum[l_row,l_col] + 1
                        accum[n_row,n_col] = accum[n_row,n_col] + accum[l_row,l_col] + 1
    return accum

def accumulation_method2(np.ndarray[DTYPE_int32_t,ndim=2] fd not None):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef np.ndarray accum = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int l_row,l_col,row,col,n_row,n_col
    cdef int row_off,col_off
    cdef int dir_i,init_dir_i
    init_dir_i = 1
    for row in range(0,row_len):
        for col in range(0,col_len):
            if fd[row,col] == 0:
                l_row,l_col,dir_i = row,col,init_dir_i
            while True:
                if dir_i <= 128:
                    row_off,col_off = coordOffsets[dir_i]
                    n_row = l_row+row_off
                    n_col = l_col+col_off
                    if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len) or fd[n_row,n_col] == NodataValue:
                        l_row,l_col,dir_i = l_row,l_col,dir_i * 2
                    else:
                        if fd[n_row,n_col] == coordOffsets_flows_to[row_off,col_off]:
                            l_row,l_col,dir_i = n_row,n_col,init_dir_i
                        else:
                            l_row,l_col,dir_i = l_row,l_col,dir_i * 2
                else: # all neighbours have been processed
                    #write my accumulation downstream
                    if fd[l_row,l_col] != 0:
                        row_off,col_off = coordOffsets[fd[l_row,l_col]]
                        n_row = l_row+row_off
                        n_col = l_col+col_off
                        accum[n_row,n_col] = accum[n_row,n_col] + accum[l_row,l_col] + 1
                        #queue up the downstream cell
                        l_row,l_col,dir_i = n_row,n_col,coordOffsets_flows_to[row_off,col_off] * 2
                    else:
                        break  #back at the starting sink
    return accum

def flowd_remap(np.ndarray[DTYPE_int32_t,ndim=2] fd not None):
    coordremap = {0:-1, 1:0, 2:1, 4:2, 8:3, 16:4, 32:5, 64:6, 128:7}
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef int row,col
    for row in range(0,row_len):
        for col in range(0,col_len):
            fd[row,col] = coordremap[fd[row,col]]

def accumulation_method3(np.ndarray[DTYPE_int32_t,ndim=2] fd not None):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef np.ndarray accum = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int l_row,l_col,row,col,n_row,n_col
    cdef int row_off,col_off
    cdef int dir_i,init_dir_i
    init_dir_i = 0
    offsets_list = [(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
    for row in range(0,row_len):
        for col in range(0,col_len):
            if fd[row,col] == -1:
                l_row,l_col,dir_i = row,col,init_dir_i
            while True:
                if dir_i != fd[l_row,l_col]:
                    row_off,col_off = offsets_list[dir_i]
                    n_row = l_row+row_off
                    n_col = l_col+col_off
                    if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len):
                        l_row,l_col,dir_i = l_row,l_col,(dir_i + 1) % 8
                    else:
                        if fd[n_row,n_col] == (dir_i + 4) % 8:
                            l_row,l_col,dir_i = n_row,n_col,(fd[n_row,n_col] + 1) % 8
                        else:
                            l_row,l_col,dir_i = l_row,l_col,(dir_i + 1) % 8
                else: # all neighbours have been processed
                    #write my accumulation downstream
                    if fd[l_row,l_col] != -1:
                        row_off,col_off = offsets_list[fd[l_row,l_col]]
                        n_row = l_row+row_off
                        n_col = l_col+col_off
                        accum[n_row,n_col] = accum[n_row,n_col] + accum[l_row,l_col] + 1
                        #queue up the downstream cell
                        l_row,l_col,dir_i = n_row,n_col,(fd[l_row,l_col] + 5) % 8
                    else:
                        break #back at the starting sink
    return accum

#assumes single flow direction
def saddles (np.ndarray[DTYPE_t,ndim=2] dem not None,np.ndarray[DTYPE_int_t,ndim=2] catchments not None):
    cdef int row_len = dem.shape[0]
    cdef int col_len = dem.shape[1]
    cdef np.ndarray catch = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int row,col,n_row,n_col
    cdef int row_off,col_off,high_row,high_col,low_row,low_col,high_cat,low_cat
    cdef float elev,slope,sqrt2
    sqrt2 =  math.sqrt(2)

    saddles = {}
    #skip on the edge cells - saves on bounds checking - and i dont care about them
    #but doing the first row as im only in half the directions
    #NOTE will miss a saddle into the tip left corner.....
    for row in range(0,row_len-1):
        for col in range(1,col_len-1):
            for row_off,col_off in half_coordOffsets_list:
                n_row = row+row_off
                n_col = col+col_off
                if catchments[n_row,n_col] != catchments[row,col]:
                    #put this combo in the dictionary
                    #the key is the pair of catchment ids (sorted)
                    if catchments[row,col] < catchments[n_row,n_col]:
                        ID = (catchments[row,col],catchments[n_row,n_col])
                    else:
                        ID = (catchments[n_row,n_col],catchments[row,col])
                    #highest elevation
                    if dem[row,col] < dem[n_row,n_col]:
                        elev = dem[n_row,n_col]
                        slope = elev - dem[row,col]
                        specs = dict(elev=elev,slope=slope,
                            high_row = n_row,
                            high_col = n_col,
                            high_cat = catchments[n_row,n_col],
                            low_row = row,
                            low_col = col,
                            low_cat = catchments[row,col])
                    else:
                        elev = dem[row,col]
                        slope = elev - dem[n_row,n_col]
                        specs = dict(elev=elev,slope=slope,
                            high_row = row,
                            high_col = col,
                            high_cat = catchments[row,col],
                            low_row = n_row,
                            low_col = n_col,
                            low_cat = catchments[n_row,n_col])
                    if 0 not in [row_off,col_off]:
                        slope = slope / sqrt2
                    #if there is already an entry in the dictionary for that pair
                    #replace the elevation if this one is lower
                    if ID in saddles:
                        if saddles[ID]['elev'] > elev:
                            saddles[ID] = specs
                        elif saddles[ID]['elev'] == elev and slope > saddles[ID]['slope']:
                            saddles[ID] = specs
                    else:
                        saddles[ID] = specs
    return saddles


#ignore outside edge of pixels
#if it flows from 3 cells in to 2 cells in its an outlet
#if it flows from 2 cells in to 3 cells in its an inlet
def flow_band(np.ndarray[DTYPE_int32_t,ndim=2] fd not None,int souce_band, int dest_band):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    flowList = []

    #TODO: souce_band and dest_band should be neighbours - check
    #TODO: error check should be based on and souce_band and dest_band
    if row_len < 6 or col_len < 6:
        #the array is too small
        return
    cdef int s = souce_band
    cdef int d = dest_band
    cdef int d_row_max = row_len-d
    cdef int d_col_max = col_len-d
    cdef int row,col,n_row,n_col
    cdef int row_off,col_off

    for row in [s,row_len-s]:
        for col in range(s,col_len-s):
            if fd[row,col] is not 0:
                row_off,col_off = coordOffsets[fd[row,col]]
                n_row = row+row_off
                n_col = col+col_off
                if n_row in [d, d_row_max] and n_col >=d and n_col <= d_col_max:
                    flowList.append((n_row,n_col))
    for row in range(s+1,row_len-(s+1)): #dont reprocess the corner cells
        for col in (s,col_len-s):
            if fd[row,col] is not 0:
                row_off,col_off = coordOffsets[fd[row,col]]
                n_row = row+row_off
                n_col = col+col_off
                if n_col in [d, d_col_max] and n_row >= d and n_row <= d_row_max:
                    flowList.append((n_row,n_col))
    return flowList



#assumes single flow direction
def sink_catchments(np.ndarray[DTYPE_int32_t,ndim=2] fd not None):
    cdef int row_len = fd.shape[0]
    cdef int col_len = fd.shape[1]
    cdef np.ndarray catch = np.zeros([row_len,col_len], dtype=DTYPE_int)
    cdef int l_row,l_col,row,col,n_row,n_col
    cdef int row_off,col_off
    cdef int sinkID = 0
    todo = []
    for row in range(row_len):
        for col in range(col_len):
            if fd[row,col] == 0:
                sinkID = sinkID + 1
                catch[row,col] = sinkID
                todo.append((row,col))
            while todo:
                l_row,l_col = todo.pop()
                #search upstream
                for row_off,col_off in coordOffsets_list:
                    n_row = l_row+row_off
                    n_col = l_col+col_off
                    if not (n_row >= 0 and n_row < row_len and n_col >= 0 and n_col < col_len) or fd[n_row,n_col] == NodataValue:
                        pass #not interested in no data or off edges data??
                    elif fd[n_row,n_col] == coordOffsets_flows_to[row_off,col_off]:
                        catch[n_row,n_col] = catch[l_row,l_col]
                        todo.append((n_row,n_col))
    return catch

###
###def fd_test():
###    import gdal, gdalconst
###    fdgrid = gdal.Open(r'H:\rea095\brind\brind', gdalconst.GA_ReadOnly)
###    dem = fdgrid.GetRasterBand(1).ReadAsArray(0,0,fdgrid.RasterXSize,fdgrid.RasterYSize)
###    fd = d8(dem)
###    p = sink_catchments(fd)
###    return p
