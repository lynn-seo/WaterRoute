#TODO:
#handle junctions for new streams
#TODO:
#couple of streams try to flow the wrong way down a vector and we end up not deleting a vector and getting duplicates
#not sure on the solution to this one
#could use our same loop analysis to identify loops with most of its length marked for keeping
#TODO:
#have had the wrong stream flow into a sink of a different stream.

#TODO:
#final streams can have some funniny kinks in them where they drop into a sink and then climb back out
#could simplify final lines or
#burn them into the dem again and reproduce them
#it would short cut any bits where it can jump to a neighbouring pixel.
#or
#first time round you could enforce edge conditions like anudem does.

#FIXED
#stop streams at a distance threshold if they flow past a sink -DONE: these have now been flagged
#also stop them if they flow to a different tree in the forest of a very different part of the same tree (maybe)
#that would deal with streams that go the wong side of a major divide. -DONE: these have now been flagged


import collections
import sys
from osgeo import gdal, gdalconst, ogr
#import gdal, gdalconst
import numpy as np
import scipy
sys.path.append(r'E:\temp\svn\TerrainAnalysis\various_tools\Trunk\vector2raster')
import vector2raster
import shapely.wkb
from shapely.geometry import Point,LineString
#import osgeo.ogr as ogr
import networkx as nx
sys.path.append(r'E:\temp\svn\TerrainAnalysis\various_tools\Trunk\tiled_flow')
import drain
import flow1_c
sys.path.append(r'E:\temp\svn\TerrainAnalysis\various_tools\Trunk\NASYdams')
import FlowOrdertest
import FlowUtilities
import ogr

"""
from 
https://github.com/om-henners/ndimage_talk/blob/master/calculate_slope.py
and
http://www.slideshare.net/JHasthorpe/raster-processing-with-scipyndimage-dev-meet-ii

slope = generic_filter(rast_data, calc_slope, size=3, extra_arguments=(x_cellsize, y_cellsize))
"""
# def calc_slope(in_filter, x_cellsize, y_cellsize):
    # #slope calculation here - note need to reshape in array to be 3*3
    # if 0.0 in in_filter:
        # return -9999
    # else:
        # [a, b, c, d, e, f, g, h, i] = in_filter
    
        # dz_dx = ((c + 2*f + i) - (a + 2 * d + g)) / (8 * x_cellsize)
        # dz_dy = ((g + 2*h + i) - (a + 2 * b + c)) / (8 * y_cellsize)
        # slope = np.sqrt(dz_dx ** 2 + dz_dy**2)
    
        # return np.degrees(slope)

def build_graph(streams_path,streams_ID):
    #read streams and build a directed graph
    str_source = ogr.Open(streams_path)
    str_layer = str_source.GetLayer()
    str_layer.ResetReading()
    str_defn = str_layer.GetLayerDefn()
    ID_index = str_defn.GetFieldIndex(streams_ID)
    DG=nx.MultiDiGraph()
    for feat in str_layer:
        ID = feat.GetFieldAsInteger(ID_index)
        stream = shapely.wkb.loads(feat.GetGeometryRef().ExportToWkb())
        #build topology from geometry of line rather then FROM/TO node attributes
        #stream digitising to be in direction of flow with identical stream end coords
        FROM = stream.coords[0]
        TO = stream.coords[-1]
        DG.add_edge(FROM, TO, id=ID,stream=stream)

    return DG
"""
Cuts a line in two at a distance from its starting point

from http://sgillies.net/blog/1040/shapely-recipes/
distance >= line.length this test fails as times due to floating point arithmatic
dont trust it
 """ 
def cut(line, distance):
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]

"""
a variation to cut that accepts a point rather then a distance

The point gets used as the vertex at the cut rather then an interpolated point.
this ensure the point is snaped to and avoid floating point numerical differences
"""
def cutPoint(line, point):
    distance = line.project(point)
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            if i == 0 or i == len(coords) - 1:
                break
            return [
                LineString(coords[:i] + [(point.x, point.y)]),
                LineString([(point.x, point.y)] + coords[i+1:])]
        if pd > distance:
            return [
                LineString(coords[:i] + [(point.x, point.y)]),
                LineString([(point.x, point.y)] + coords[i:])]
    return [LineString(line)]

def cutPoints_old(line, cutPoints):
    if type(cutPoints) is Point:
        cutPoints = [cutpoints]
    cuts = sorted([line.project(p) for p in cutPoints])
    #work out the added distance for each sorted cut
    cuts_inc = cuts[:1]
    cuts_inc.extend( i - cuts[id-1] for id, i in enumerate(cuts[1:],1))
    segs = []
    for c in cuts_inc:
        parts = cut(line,c)
        if len(parts) ==2:
            seg,line = parts
            if seg.length > 0.0:
                segs.append(seg)
        else:
            line = parts[0]
    #keep whats left of the line
    if line.length > 0.0:
        segs.append(line)
    return segs
"""
cuts a line at a series of points

if you know all the points need to cut this line then it saves 
work in working out which half to cut with the second point after you have already cut it
"""        
def cutPoints(line, cutPoints):
    if type(cutPoints) is Point:
        cutPoints = [cutPoints]
    cuts = sorted([(line.project(p),p) for p in cutPoints])
    #work out the added distance for each sorted cut
    segs = []
    for d,p in cuts:
        parts = cutPoint(line,p)
        if len(parts) ==2:
            seg,line = parts
            if seg.length > 0.0:
                segs.append(seg)
        else:
            line = parts[0]
    #keep whats left of the line
    if line.length > 0.0:
        segs.append(line)
    return segs

"""
assumes no cell flows out of the raster
to ensure this all edge cells should be sinks
works with the special case of wanting the flow between features in 'ends'
TODO: a little bit slow - convert to cython
TODO: need to use accumulation to handle junctions correctly between 'ends'
"""
def raster2stream(fd,ends):
    from shapely.geometry import LineString
    coordOffsets = {1:(0,1), 2:(1,1), 4:(1,0), 8:(1,-1), 16:(0,-1), 32:(-1,-1), 64:(-1,0), 128:(-1,1)}
    row_len = fd.shape[0]
    col_len = fd.shape[1]
    links = []
    #skip on the edge cells - assumed to always be sinks
    for row in range(1,row_len-1):
        for col in range(1,col_len-1):
            from_row,from_col = row,col

            if ends[from_row,from_col] != 0:
                startVal = ends[from_row,from_col]
                endVal = 0
                coordList = [FlowUtilities.RCtoXY(trans2,(from_row,from_col))]
                from_fd = fd[from_row,from_col]
                while from_fd != 0:
                    row_off,col_off = coordOffsets[from_fd]
                    to_row = from_row+row_off
                    to_col = from_col+col_off
                    coordList.append(FlowUtilities.RCtoXY(trans2,(to_row,to_col)))
                    if ends[to_row,to_col] != 0:
                        #flow off network
                        endVal = ends[to_row,to_col]
                        break
                    #get ready for next iteration
                    from_row,from_col = to_row,to_col
                    from_fd = fd[from_row,from_col]
                if len(coordList) > 2:
                    links.append ((LineString(coordList),startVal,endVal))
    return links

    
print 'start'
raster_file_source = r'\\wron\gis\national\dem\SRTM_derived\1sec\V1_0_official_release\DEM-S\Mosaic\dems1sv1_0'
shp_file_source = r"E:\temp\drain\tas\streams.shp"
IdField = "AUSHYDRO_I"
cell_size_in_m = 30
g = gdal.Open(raster_file_source, gdalconst.GA_ReadOnly)
trans = g.GetGeoTransform()

DG = build_graph(shp_file_source,IdField)
stream_bC,stream_ras,stream_mapping = vector2raster.Streams2raster(DG,trans[0],trans[3],trans[1])
stream_w_ids = stream_mapping[stream_ras]
id_lookup = {e_data['id']: (e_from,e_to,e_key,e_data) for e_from,e_to,e_key,e_data in DG.edges(keys=True,data=True)}

ulX,ulY,lrX,lrY,pxWidth,pxHeight,minX_snap,minY_snap,maxX_snap,maxY_snap = stream_bC
dem = g.GetRasterBand(1).ReadAsArray(int(ulX),int(ulY),int(pxWidth),int(pxHeight))

#im assuming nodata is just ocean
#todo should maintain this nodata throughout the processing as sinks
r = g.GetRasterBand(1)
dem[dem == r.GetNoDataValue()] = 0

trans2 = list(trans)
trans2[0] = minX_snap
trans2[3] = maxY_snap

#this was a good idea but i had some issues with it 
#to be investigated
## create raster of cells that are the  outlets of lower end streams.
outs = [Point(item[0]) for item in DG.out_degree_iter() if item[1] == 0]
orig_bC,bottoms1 = vector2raster.vectorList2raster(outs,minX_snap,maxY_snap,trans[1])
bottoms = np.zeros(dem.shape,dtype=np.int)
oulX,oulY,olrX,olrY,opxWidth,opxHeight,ominX_snap,ominY_snap,omaxX_snap,omaxY_snap = orig_bC
# this could fail if the new array somehow go outside the edge of the other array
bottoms[oulY:olrY,oulX:olrX] = bottoms1
drain.gdalWrite(bottoms,r"E:\temp\drain\tas\bottoms.tif",trans2)

## create raster of cells that are the  inlets of top end streams.
ins = [Point(item[0]) for item in DG.in_degree_iter() if item[1] == 0]
orig_bC,tops1 = vector2raster.vectorList2raster(ins,minX_snap,maxY_snap,trans[1])
tops = np.zeros(dem.shape,dtype=np.int)
oulX,oulY,olrX,olrY,opxWidth,opxHeight,ominX_snap,ominY_snap,omaxX_snap,omaxY_snap = orig_bC
# this could fail if the new array somehow go outside the edge of the other array
tops[oulY:olrY,oulX:olrX] = tops1
drain.gdalWrite(tops,r"E:\temp\drain\tas\tops.tif",trans2)

print 'tops and bottoms done'

print 'dem burning'
burn_dem = np.copy(dem)
burn_dem[stream_ras > 0] -= 1000
#dont burn in edge cells
burn_dem[0] = dem[0]
burn_dem[-1] = dem[-1]
burn_dem[:,0] = dem[:,0]
burn_dem[:,-1] = dem[:,-1]


print 'dilation and masking'

#mask streams on steep slopes (using local relief in this case)
#todo: develop good method to decide what to throw away.
#small gaps in the network give problamatic results
#small bits remaining shouldnt be an issue

# mask is 10 pixels along stream from any pixel that is over the relief threshold (including neighbouring pixels)
#pixels directly next to the stream and over the threshold still trigger the mask.
#if you dont want this behaviour you need to filter to stream_ras first
dem_max = scipy.ndimage.filters.maximum_filter(dem,3)
dem_min = scipy.ndimage.filters.minimum_filter(dem,3)
dem_range = dem_max-dem_min
mask = np.zeros(dem.shape)
mask[dem_range > 20] = 1
expand_mask = scipy.ndimage.morphology.binary_dilation(mask,scipy.ndimage.generate_binary_structure(2, 2),10,stream_ras)
expand_mask[stream_ras == 0] = 0
drain.gdalWrite(expand_mask,r"E:\temp\drain\tas\maskexpandclip.tif",trans2)

not_mask = np.logical_not(mask)
not_mask[stream_ras == 0] = 0



# #accumulate weighted by masked streams using fd burnt but not masked
# #find accumulations of one
# #save these for seeds for later accumulation
# #these are the tops of the original network except where the top is masked
# #its the first bit of the network that is unmasked.
# burn_fd = drain.fixedFD(burn_dem,bottoms)
# flow1_c.flowd_remap(burn_fd) #need to recode stuff to use the same fd scheme
# ok_accum,ok_count = FlowOrdertest.upstreamAverage(burn_fd,not_mask.astype(np.dtype('float32')),-9999)
# #TODO: seeds occuring along channels where there is a one pixel shorter route
# #also at flow splits
# #TODO: doesnt handle flow splits.

print 'remove burning in masked areas'
burn_dem[expand_mask] = dem[expand_mask]

print 'flow direction'
mask_burn_fd = drain.fixedFD(burn_dem)#,bottoms)

print 'raster stream generation'
#use earlier seeds with masked fd
#this gives paths from top of streams around gaps
mask_burn_fd2 = np.copy(mask_burn_fd)
flow1_c.flowd_remap(mask_burn_fd2) #need to recode stuff to use the same fd scheme
##not_mask[ok_accum != 1] = 0
##seed_accum,seed_count = FlowOrdertest.upstreamAverage(mask_burn_fd2,not_mask.astype(np.dtype('float32')),-9999)
seed_accum,seed_count = FlowOrdertest.upstreamAverage(mask_burn_fd2,tops.astype(np.dtype('float32')),-9999)

#Extract just the streams around the gaps
drain.gdalWrite(seed_accum,r"E:\temp\drain\tas\streams.tif",trans2)
#TODO: need to force it to stop if the stream stops
#could put a nodata sink at the end of such a stream
#would only stop it if flow was in that stream at that stage.
#since that tends to happen on flatter areas that seems reasonable.
#if the goal isnt input into anudem then this doesnt matter as we can just flag stream that are climbing out of a depression.


print 'create vector lines for the new stream segments'
top_ids = np.copy(stream_w_ids)
stream_w_ids[expand_mask] = 0
stream_w_ids[seed_accum == 0] = 0
stream_w_ids[tops != 0] = top_ids[tops != 0]
drain.gdalWrite(stream_w_ids,r"E:\temp\drain\tas\stream_w_ids.tif",trans2)
links = raster2stream(mask_burn_fd,stream_w_ids)

print 'snap the ends of the lines to the original vector linework'
links2 = []
cutSites= collections.defaultdict(list)
for l,f,t in links:
    if f != 0:
        fLine = id_lookup[f][3]['stream']
        cp = fLine.interpolate(fLine.project(Point(l.coords[0])))
        l = LineString([(cp.x, cp.y)] + l.coords[1:])
        cutSites[f].append(cp)
    if t != 0:
        tLine = id_lookup[t][3]['stream']
        cp = tLine.interpolate(tLine.project(Point(l.coords[-1])))
        l = LineString(l.coords[:-1] + [(cp.x, cp.y)])
        cutSites[t].append(cp)
    links2.append(l)
    
print 'split the original linework at the new intersections'
cut_streams = { id: cutPoints(id_lookup[id][3]['stream'],points) for id,points in cutSites.items() }

print 'rebuild the graph'
DG2=nx.MultiDiGraph()
#new lines
for stream in links2:
    FROM = stream.coords[0]
    TO = stream.coords[-1]
    DG2.add_edge(FROM, TO, id=-1,stream=stream,source=100,length=stream.length)
#cut lines
for id,streamList in cut_streams.items():
    for stream in streamList:
        FROM = stream.coords[0]
        TO = stream.coords[-1]
        DG2.add_edge(FROM, TO, id=id,stream=stream,source=10,length=stream.length)

#unmodified lines        
id_lookup2 = {e_data['id']: (e_from,e_to,e_key,e_data) for e_from,e_to,e_key,e_data in DG2.edges(keys=True,data=True)}
gen = ((id,i[3]['stream']) for id,i in id_lookup.items() if id not in id_lookup2)
for id,stream in gen:
        FROM = stream.coords[0]
        TO = stream.coords[-1]
        DG2.add_edge(FROM, TO, id=id,stream=stream,source=1,length=stream.length)


print 'write some output diagnostic streams'
out_streams = [e_data['stream'] for e_from,e_to,e_key,e_data in DG2.edges(keys=True,data=True)]# if e_data['keep'] == True]
driverName = "ESRI Shapefile"
drv = ogr.GetDriverByName(driverName)
ds = drv.CreateDataSource(r"E:\temp\drain\tas")
lyr = ds.CreateLayer('out_streams_diag2',None,ogr.wkbLineString)

fDef = lyr.GetLayerDefn()
for i in out_streams:
    feat = ogr.Feature(fDef)
    feat.SetGeometry(ogr.CreateGeometryFromWkb(i.wkb))
    lyr.CreateFeature(feat)
    feat.Destroy()
ds=None

        
        
print 'hacked topological sort to remove cycle nodes '
#not a good fix but a temporary solution.

seen={}
order_explored=[] # provide order and 
explored={}       # fast search without more general priorityDictionary
cycleNodes = []
nbunch = DG2.nodes_iter() 
for v in nbunch:     # process all vertices in G
    if v in explored: 
        continue
    fringe=[v]   # nodes yet to look at
    while fringe:
        w=fringe[-1]  # depth first search
        if w in explored: # already looked down this branch
            fringe.pop()
            continue
        seen[w]=1     # mark as seen
        # Check successors for cycles and for new nodes
        new_nodes=[]
        for n in DG2[w]:
            if n not in explored:
                if n in seen: #CYCLE !!
                    print n
                    explored[n]=1
                    cycleNodes.append(n)
                    #raise nx.NetworkXUnfeasible("Graph contains a cycle.")
                new_nodes.append(n)
        if new_nodes:   # Add new_nodes to fringe
            fringe.extend(new_nodes)
        else:           # No new nodes so w is fully explored
            explored[w]=1
            order_explored.insert(0,w) # reverse order explored
            fringe.pop()    # done considering this node

for n in cycleNodes:
    DG2.remove_node(n)

if len(cycleNodes):
    print 'exporting cycles for diagnostics'         
    driverName = "ESRI Shapefile"
    drv = ogr.GetDriverByName(driverName)
    ds = drv.CreateDataSource(r"E:\temp\drain\tas")
    lyr = ds.CreateLayer('cycle_locations2',None,ogr.wkbPoint)
    fDef = lyr.GetLayerDefn()
    for i in cycleNodes:
        feat = ogr.Feature(fDef)
        feat.SetGeometry(ogr.CreateGeometryFromWkb(Point(i).wkb))
        lyr.CreateFeature(feat)
        feat.Destroy()
    ds=None



print 'traverse network working out which edges to keep and which are no longer required'
for n in nx.topological_sort(DG2):
    if DG2.out_degree(n,weight='source') >= 100:
            keep = False
    else:
        keepList = [e_data['keep'] for e_from,e_to,e_key,e_data in DG2.in_edges_iter(n,True,True)]
        if len(keepList) == 0 or True in keepList:
            keep = True
        else:
            keep=False
    for e_from,e_to,e_key,e_data in DG2.out_edges_iter(n,True,True):
        e_data['keep'] = keep
        if e_data['source'] >= 100:
            e_data['keep'] = True

            
print 'analyse the path lengths to find out if the streams come back together'
for e_from,e_to,e_key,e_data in DG2.edges_iter(data=True,keys=True):
    e_data['dist'] = 9999
    e_data['keep_dist'] = 9999
    if e_data['keep'] == True:
        e_data['keep_length'] = e_data['length']
    else:
        e_data['keep_length'] = 0 

        #todo there is an issue the path needs to prefer to go donw kept paths
for n,deg in DG2.out_degree_iter(weight='source'):
    if deg >= 100:
        for e_from,e_to,e_key,e_data in DG2.out_edges(n,True,True):
            if e_data['source'] ==100:
                #will only ever get here once per node
                e_from1,e_to1,e_key1,e_data1 = e_from,e_to,e_key,e_data
                pathLengths1 = nx.shortest_path_length(DG2, source=e_to1, target=None, weight='length')
                keep_pathLengths1 = nx.shortest_path_length(DG2, source=e_to1, target=None, weight='keep_length')
                if deg == 100:
                    #new stream coming out of the end of network (ie root of an isolated tree)
                    e_data1['dist'] = 9997
        
        for e_from,e_to,e_key,e_data in DG2.out_edges(n,True,True):
            if e_data['source'] !=100:
                e_from2,e_to2,e_key2,e_data2 = e_from,e_to,e_key,e_data
                pathLengths2 = nx.shortest_path_length(DG2, source=e_to2, target=None, weight='length')
                keep_pathLengths2 = nx.shortest_path_length(DG2, source=e_to2, target=None, weight='keep_length')

                candidates = []
                for key in pathLengths1:
                    if key in pathLengths2:
                        candidates.append([pathLengths1[key] + e_data1['length'] + pathLengths2[key] + e_data2['length'],pathLengths1[key] + e_data1['length'],pathLengths2[key] + e_data2['length'],key,keep_pathLengths1[key] + e_data1['keep_length'] + keep_pathLengths2[key] + e_data2['keep_length']])
                if len (candidates) == 0:
                    # the new and old streams flow to different outlets
                    #print "no path", n
                    #gchek if freshly initialised
                    if e_data1['dist'] == 9999: 
                        e_data1['dist'] = 9998
                else:
                    #the shortest loop found with a minimum distance.
                    dist,dist1,dist2,key,keep_dist = min(candidates)
                    #create a 0 to 1 index as a measure of the path difference 0 is identical 1 is very different
                    #values over about 0.15 are worth checking
                    dist_index = abs(0.5 - e_data1['length']/dist) * 2
                    if dist_index < e_data1['dist']:
                        e_data1['dist'] = dist_index
                        e_data1['keep_dist'] = keep_dist/dist







print 'write some output'
out_streams = [(e_data['stream'],e_data['keep'],e_data['dist'],e_data['keep_dist']) for e_from,e_to,e_key,e_data in DG2.edges(keys=True,data=True)]# if e_data['keep'] == True]
driverName = "ESRI Shapefile"
drv = ogr.GetDriverByName(driverName)
ds = drv.CreateDataSource(r"E:\temp\drain\tas")
lyr = ds.CreateLayer('out_streams_vec',None,ogr.wkbLineString)
defn = ogr.FieldDefn( 'x', ogr.OFTInteger)
lyr.CreateField(defn)
defn = ogr.FieldDefn( 'dist', ogr.OFTReal)
lyr.CreateField(defn)
defn = ogr.FieldDefn( 'keep_dist', ogr.OFTReal)
lyr.CreateField(defn)

fDef = lyr.GetLayerDefn()
for i in out_streams:
    feat = ogr.Feature(fDef)
    feat.SetField('x',int(i[1]))
    feat.SetField('dist',i[2])
    feat.SetField('keep_dist',i[3])
    feat.SetGeometry(ogr.CreateGeometryFromWkb(i[0].wkb))
    lyr.CreateFeature(feat)
    feat.Destroy()
ds=None


"""
ways to get fg streams
filter out bad streams - do it well ie remove short remaining bits
burn good streams - stops multiple lines occuring

rest of fg technique would be ok

seed from top of remaining network instead of from everywhere

network analysis of alternate paths to check sanity

its all a bit pointless
why not burn in the good streams and then create a detailed network then provide some info on sinks in the network and meaningful environmental info to allow users to define the extent of the network.
The position of the blue lines in flat areas is useful information. The extent of the blue line network is arbitrary in many cases. Use the positions where we need that information ignore the extent


dont forget distributaries and anabranches
"""