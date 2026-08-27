# -*- coding: utf-8 -*-
"""
a tool for simplifying the geo fabric to course and courser 'levels'

a good basic example of stream network simplification.
deals with flow bifurcations rather crudely but simpily.

"""
#defaults use test data from the BOM
streams_path = r'S:\rea095\geofabric\AHGFNetworkStreams.shp' #r"bom\AHGFNetworkStream.shp"
streams_ID = "DrainID"
subcatchs_path = r"S:\rea095\geofabric\AHFGCatchment.shp"
subcatchs_ID = "HydroID"
outshp = r'S:\rea095\geofabric\outV2_'

import shapely.wkb
from shapely.geometry import Polygon
import osgeo.ogr
import networkx as nx
import datetime

def build_graph(streams_path,subcatchs_path,streams_ID,subcatchs_ID):
    #read watersheds to a dictionary indexed by their ID
    wshd_source = osgeo.ogr.Open(subcatchs_path)
    wshd_layer = wshd_source.GetLayer()
    wshd_layer.ResetReading()
    wshd_defn = wshd_layer.GetLayerDefn()
    ID_index = wshd_defn.GetFieldIndex(subcatchs_ID)
    subcatch_dict = {}
    for feat in wshd_layer:
        ID = feat.GetFieldAsInteger(ID_index)
        subcatch = shapely.wkb.loads(feat.GetGeometryRef().ExportToWkb())
        if ID in subcatch_dict:
            #is features have the same ID make multipart shapes
            subcatch_dict[ID] = subcatch.union(subcatch_dict[ID])
        else:
            subcatch_dict[ID] = subcatch
    #read streams and build a directed graph and match to the watersheds
    str_source = osgeo.ogr.Open(streams_path)
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
        if ID in subcatch_dict:
            DG.add_edge(FROM, TO, id=ID,stream=stream,subcatch=subcatch_dict[ID])
        else:
            DG.add_edge(FROM, TO, id=ID,stream=stream,subcatch=Polygon())
    return DG


def run_levels():
    starttime = datetime.datetime.now().isoformat()
    print 'building graph 1'
    DG = build_graph(streams_path,subcatchs_path,streams_ID,subcatchs_ID)
    shortest_path_to_furthest_outlet(DG)
    remove_pseudo_nodes(DG)
    write_to_shape(outshp + '000',DG)
    for x in xrange(15):
        simplify(DG)
        remove_pseudo_nodes(DG)
        write_to_shape(outshp + str(x),DG)
    print starttime
    print datetime.datetime.now().isoformat()

def simplify(DG):
    for item in DG.in_degree_iter():
        if item[1] == 0:
            for e1_from,e1_to,e1_key,e1_data in DG.out_edges_iter(item[0],data=True,keys=True):
                out_es =DG.out_edges(e1_to,data=True,keys=True)
                if len(out_es) > 0:
                    e2_from,e2_to,e2_key,e2_data = out_es[0]
                    e2_data['subcatch'] = e2_data['subcatch'].union(e1_data['subcatch'])
                    #dont tab this in if you want to remove streams once the rest of their network is gone
                DG.remove_edge(e1_from,e1_to,e1_key)
                print e1_from,e1_to,e1_key
    #now go through the graph again and remove psuedo nodes and parallel edges

def remove_pseudo_nodes(DG):
    found = True
    while found:
        print 'loop started'
        found = False
        for n in DG:
            out_es =DG.out_edges(n,data=True,keys=True)
            in_es =DG.in_edges(n,data=True,keys=True)
            if len(in_es) == 0 and len(out_es) == 0:
                DG.remove_node(n)
            if len(in_es) == 1 and len(out_es) == 1:
                found = True
                e2_from,e2_to,e2_key,e2_data =out_es[0]
                e1_from,e1_to,e1_key,e1_data =in_es[0]
                DG.add_edge(e1_from, e2_to, id=e2_data['id'],stream=e2_data['stream'].union(e1_data['stream']),subcatch=e2_data['subcatch'].union(e1_data['subcatch']))
                DG.remove_node(n)
                print 'node',n
            #get rid of parellel flow paths
            if len(in_es) > 1:
                for e1_from,e1_to,e1_key,e1_data in in_es:
                    if e1_key > 0:
                        found = True
                        #arbitarily just keep the first one
                        DG[e1_from][e1_to][0]['subcatch'] =  DG[e1_from][e1_to][0]['subcatch'].union(e1_data['subcatch'])
                        DG.remove_edge(e1_from,e1_to,e1_key)
                        print 'edge',e1_from,e1_to,e1_key

#def dup_finder(DG):
#    stream_dict = {}
#    for e1_from,e1_to,e1_data in DG.edges(data=True):
#        if e1_data["id"] in stream_dict:
#            print e1_data["id"]
#        else:
#            stream_dict[e1_data["id"]] = (e1_from,e1_to)

def mean(nums):
    if len(nums):
        return float(sum(nums)) / len(nums)
    else:
        return 0.0

def shortest_path_to_furthest_outlet(DG):
    print 'sorting'
    ts1 = nx.topological_sort(DG)
    print 'finished sorting ... searching....'
    for n1 in reversed(ts1):
        d = {}
        dist = 0
        outlet = n1
        out_es =DG.out_edges(n1,data=True,keys=True)
        for  e1_from,e1_to,e1_key,e1_data in out_es:
            outlet,dist = e1_data['outlet']
            if outlet in d:
                if dist < d[outlet]:
                    d[outlet] = (dist,e1_from,e1_to)
            else:
                d[outlet] = (dist,e1_from,e1_to)
        if len(d) > 0:
            outlet = sorted(d, key=d.get, reverse=True)[0]
            dist,elast_from,elast_to = d[outlet]
        for  e1_from,e1_to,e1_key,e1_data in out_es:
            if (e1_from,e1_to) != (elast_from,elast_to):
                #TODO:need to union the catchment downstream
                DG.remove_edge(e1_from,e1_to,e1_key)
                out_es2 =DG.out_edges(e1_to,data=True,keys=True)
                if len(out_es2) > 0:
                    out_es2[0][3]['subcatch'] = out_es2[0][3]['subcatch'].union(e1_data['subcatch'])
                print 'removed',e1_from,e1_to,e1_key
        for  e1_from,e1_to,e1_key,e1_data in DG.in_edges_iter(n1,data=True,keys=True):
            e1_data['outlet'] = (outlet,dist + e1_data['stream'].length)


def write_to_shape(outFile,DG):
    import sys
    import ogr
    import string

    driverName = "ESRI Shapefile"
    drv = ogr.GetDriverByName( driverName )
    ds = drv.CreateDataSource( outFile )
    lyr = ds.CreateLayer( "streams", None, ogr.wkbMultiLineString )
    field_defn = ogr.FieldDefn( "ID", ogr.OFTInteger )
    lyr.CreateField ( field_defn )
    featureDefinition = lyr.GetLayerDefn()

    lyr2 = ds.CreateLayer( "subcatch", None, ogr.wkbMultiPolygon )
    field_defn2 = ogr.FieldDefn( "ID", ogr.OFTInteger )
    lyr2.CreateField ( field_defn2 )
    featureDefinition2 = lyr2.GetLayerDefn()

    for e1_from,e1_to,e1_key,e1_data in DG.edges_iter(data=True,keys=True):
        feat = ogr.Feature( featureDefinition)
        feat.SetField( "ID", e1_data['id'] )
        feat.SetGeometry(ogr.CreateGeometryFromWkb(e1_data['stream'].wkb))
        if lyr.CreateFeature(feat) != 0:
            print "Failed to create feature in shapefile.\n"
        feat.Destroy()

        feat2 = ogr.Feature( featureDefinition2)
        feat2.SetField( "ID", e1_data['id'] )
        feat2.SetGeometry(ogr.CreateGeometryFromWkb(e1_data['subcatch'].wkb))
        if lyr2.CreateFeature(feat2) != 0:
            print "Failed to create feature in shapefile2.\n"
        feat2.Destroy()


    ds = None



