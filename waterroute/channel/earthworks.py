# seline.ng@csiro.au, start: 9/dec/2020, clayton vic
# major additions by arthur.read@csiro.au, 17/dec/2020, black mountain act
from math import sqrt
from shapely.affinity import translate
from shapely.ops import unary_union

from shapely.geometry import LineString, Polygon
from shapely.affinity import scale
import math

def getChannelFlowAndVelocity(depthFlow = 2.45, slope = 1/10000, n = 0.03, sides = 2, bedWidth = 15):

    '''
    Standalone code to check channel dimensions and flow capacity matche (Lynn, May2025)
    Function to compute the volumetric flow rate and flow velocity
    of an open channel using the Manning's Equation

    Inputs:
    -----------
    depthFlow, float = depth of the flow, m, >= 0
    slope, float = slope of the channel, m, >= 0
    n, float = Manning's roughness coefficient, >= 0
    sides, float = inverse of the bank gradient, m/m, >= 0
    bedWidth, float = width of the channel bed, m, >= 0

    Output:
    -----------
    res, dictionary
    - res['Flow'], float = channel flow, m3/s
    - res['Velocity'], float = flow velocity, m/s

    Corresponding Excel sheets:
    -----------------------------
    Channel Design Unlined
    Channel lined rock exc
    3-Channel Design Lined Buc-Th 5
    4-Channel Lined Cor to Aram 10
    1 Channel Lined Buc-Big Drop 15
    Channel Lined HGD to Buc
    '''

    # check inputs
    li = [depthFlow, slope, n, sides, bedWidth]
    if any(l < 0 for l in li):
        raise ValueError('Negative input detected. All inputs must be 0 or positive.')

    # compute channel characteristics
    widthOfFlow = bedWidth + depthFlow * (sides + sides) #m
    areaH2O = depthFlow / 2 * (bedWidth + widthOfFlow) #m2
    wettedPerimeter = bedWidth + depthFlow * (sqrt(1 + sides**2)* 2) #m
    aOverP = areaH2O / wettedPerimeter #dimensionless
    flow = areaH2O * slope**0.5 * aOverP**0.6666 / n #m3/s, mannings
    vel = flow / areaH2O #m/s

    # function output
    res = {}
    res['Flow'] = flow #m3/s
    res['Velocity'] = vel #m/s
    return res


def getCostEarthworks(DFL = 5, crossGradient = 0, bedWidth = 50, depthFlow = 4.4, freeboard = 1, sides = 2, \
    batterSides = 1.5, bermSlope = 1/20, batterHeight = 5, bermWidth = 4, bankWidth = 5, \
    overExcavateBed = 1, overExcavateWall = 3, \
    channelLined = True, rockDepth = 5):

    '''
    Function to compute the cost of earthworks as a function of the required cut and fill volumes

    Inputs:
    -----------
    DFL, float = difference between FSL and NSL, FSL-NSL, m
    crossGradient, float = cross-sectional gradient of channel
    bedWidth, float = width of the channel bed, m, >= 0
    depthFlow, float = depth of the flow, m, >= 0
    freeboard, float = vertical distance between the top of the channel and the water surface, m, >= 0
    sides, float = inverse of the gradient of the water-side of the bank, m/m, >= 0
    batterSides, float = sides of batters for deep cut 1 in, >= 0
    bermSlope, float = slope of berm, m/m
    batterHeight, float = maximum height of each batter, m, >= 0
    bermWidth, float = width of intermediate berms, m, >= 0
    bankWidth, float = width of the bank, m, >= 0
    overExcavateWall = width of over excavation parallel to side of channel, m, >= 0
    overExcavateBed = depth of over excavation below depth of cut, m, >= 0
    channelLined, boolean = True if the channel is lined, False otherwise
    rockDepth, float = 5 depth at which costs of excevation increase
   
    Output:
    -----------
    res, dictionary
    - res['Excav'], float = volume of soil excavation, m3
    - res['Fill-Normal'], float = volume of normal fill, m3
    - res['Fill-High Bank'], float = volume of high bank fill, m3
    - res['Rock Excav'], float = volume of rock excavation, m3
    - res['Total Cost'], float = total cost of earthworks, $
    - res['DFL'], float = difference between FSL and NSL, FSL-NSL, m

    Corresponding Excel sheets:
    -----------------------------
    Channel Design Unlined
    Channel lined rock exc
    3-Channel Design Lined Buc-Th 5
    4-Channel Lined Cor to Aram 10
    1 Channel Lined Buc-Big Drop 15
    Channel Lined HGD to Buc
    '''

    # ==============
    # check inputs
    # ==============

    li = [depthFlow, bedWidth, depthFlow, freeboard, sides]
    li += [batterSides, bermSlope, batterHeight, bermWidth, bankWidth]
    li += [overExcavateBed, overExcavateWall]
    if any(l < 0 for l in li):
        raise ValueError('Invalid negative input or inputs detected.')

    if type(channelLined) != type(True):
        raise ValueError('Invalid input or inputs detected.')

    # ============================================
    # over excavation figures
    # ============================================

    if channelLined == False: overExcavateWall = overExcavateBed = 0.0000000001
    overExcavateCompactFactor = 1.15

    # ==========================
    # cut and fill geometries
    # ==========================

    # get settings in right format
    settings = dict(DFLwrtNS = DFL, crossGradient = crossGradient, bedWidth = bedWidth, depthFlow = depthFlow, freeboard = freeboard, sides = sides, \
        batterSides = batterSides, bermSlope = bermSlope, batterHeight = batterHeight, bermWidth = bermWidth, bankWidth = bankWidth, \
        overExcavateBed = overExcavateBed, overExcavateWall = overExcavateWall, rockDepth=rockDepth)

    # get geometries
    res = {}
    res['Geometries'] = {}
    
    (shallow_maskedCut,maskedCut,maskedFill,maskedBorrowFill,
     shallow_maskedoverExCut,maskedOverExCut,maskedOverExFill,
     waterSection,wettedPerimeterSection,full_width) =  channel_geometry(**settings)

    res['Geometries']['shallow_maskedCut'] = shallow_maskedCut
    res['Geometries']['maskedCut'] = maskedCut
    res['Geometries']['maskedFill'] = maskedFill
    res['Geometries']['maskedBorrowFill'] = maskedBorrowFill
    res['Geometries']['shallow_maskedoverExCut'] = shallow_maskedoverExCut
    res['Geometries']['maskedOverExCut'] = maskedOverExCut
    res['Geometries']['maskedOverExFill'] = maskedOverExFill
    res['Geometries']['waterSection'] = waterSection
    res['Geometries']['wettedPerimeterSection'] = wettedPerimeterSection
    res['full_width'] = full_width
    # show geometries graphically
    #
    #sss = shapelySVG()
    #sss.add_geom(maskedCut, 0.1, 'forestgreen')
    #sss.add_geom(maskedFill, 0.1, 'red')
    #sss.add_geom(maskedBorrowFill, 0.1, 'orange')
    #sss.add_geom(wettedPerimeterSection, 0.3, 'black')
    #if maskedOverExCut: sss.add_geom(maskedOverExCut, 0.1, 'greenyellow')
    #if maskedOverExFill: sss.add_geom(maskedOverExFill, 0.1, 'hotpink')
    #sss.add_geom(waterSection, 0.1, 'blue')
    #res['Geometries'] = sss

    # ===============================================================
    # get cut and fill volumes, on a per unit length basis, m3/m
    # ===============================================================

    # channel and overexcavation cuts
    res['Rock Excav'] = maskedCut.area + maskedOverExCut.area
    res['Excav'] = shallow_maskedCut.area + shallow_maskedoverExCut.area

    # overexcavation and normal fills
    res['Fill-Normal'] = (shallow_maskedoverExCut.area + maskedOverExCut.area + maskedOverExFill.area) * overExcavateCompactFactor + maskedFill.area

    # borrow fill
    res['Fill-Borrow'] = maskedBorrowFill.area

    # ================================
    # get individual and total costs
    # ================================

    # channel and overexcavation cuts
    # scraping unit cost, $6/m3 if normal excavation, $12/m3 if rock excavation
    # costCuts = 12 * res['Rock Excav'] + 6 * res['Excav']
    costCuts = 60 * res['Rock Excav'] + 10 * res['Excav']#updated on 22/09/2025
    # overexcavation and normal fill
    # scraping unit cost, $6/m3
    # costFillNormal = 6 * res['Fill-Normal']
    costFillNormal = 12 * res['Fill-Normal'] #updated on 22/09/2025

    # borrow fill
    # scraping unit cost, $8/m3
    # costFillBorrow = 8 * res['Fill-Borrow']
    costFillBorrow = 12 * res['Fill-Borrow'] #updated on 22/09/2025

    # total cost
    res['Total Cost'] = costCuts + costFillNormal + costFillBorrow

    # =================
    # return output
    # =================

    res['DFL'] = DFL
    return res


# addition by arthur read
class shapelySVG():

    def __init__(self):
        self._geom_SVG = ''
        self._list_geom = []
    def add_geom(self,geom,scale=1.0,colour=None):
        self._list_geom.append (geom)
        self._geom_SVG += f' {geom.svg(scale,colour)}'
    def svg(self):
        from shapely.geometry import GeometryCollection
        geom = GeometryCollection(self._list_geom)
        xmin, ymin, xmax, ymax = geom.bounds
        svgNamespace = 'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"'
        svgStyle = 'style="background-color:cornsilk"'
        #transform = 'transform = "matrix(1,0,0,-1,0,0)"'
        transform = 'transform = "scale(1, -1) scale(0.8)"'
        width = 'width="100%"'
        #height = 'height="100%"'
        AR = 'preserveAspectRatio="xMidYMid meet"'
        header = f'<svg {width} {AR} viewBox="{xmin} {ymin} {xmax - xmin} {ymax - ymin}" {svgStyle} {transform} {svgNamespace}>\n'
        #header = f'<svg {width} {AR} {svgStyle} {transform} {svgNamespace}>\n'
        #header = f'{header} <g {transform}>'
        footer = f'</svg>'
        #footer = f'</g> {footer}'
        out = f'{header} {self._geom_SVG} {footer}'
        return out
    def _repr_svg_(self):
        return self.svg()


# addition by arthur read
def coord_offset(a, b):
    return tuple(p + q for p, q in zip(a, b))


# addition by arthur read
def channel_geometry(DFLwrtNS , crossGradient , bedWidth , depthFlow , freeboard , sides ,
    batterSides , bermSlope , batterHeight , bermWidth , bankWidth ,
    overExcavateBed , overExcavateWall ,rockDepth):

    '''
    (DFLwrtNS = 0, crossGradient = 0, bedWidth = 50, depthFlow = 4.4, freeboard = 1, sides = 2,
    batterSides = 1.5, bermSlope = 1/20, batterHeight = 5, bermWidth = 4, bankWidth = 5,
    overExcavateBed = 1, overExcavateWall = 3,rockDepth = 5)
    returns shapely geometry of channel features
    DFLwrtNS Design Flow Level with regard to Natural Surface. Positive values indicate the water surface is above natural ground level.
    for percent slope conversions see https://www.engineeringtoolbox.com/slope-degrees-gradient-grade-d_1562.html
    '''

    # make the design flow level above the center of the bed as origin 0,0
    origin = (0, 0)
    ##create a maximum cut section
    channelDepth = depthFlow + freeboard
    maxCut = 1000  # make it silly big - from DFL
    maxFill = maxCut  # from DFL
    maskSize = maxCut * 10  # just a big number
    bedCentre = (0, -depthFlow)
    bedEdge = coord_offset(bedCentre, (bedWidth / 2, 0))
    channelSide = coord_offset(bedEdge, ((channelDepth) * sides, channelDepth))
    channelPerimeterPoints = [bedCentre, bedEdge, channelSide]
    cutPoints = [*channelPerimeterPoints]
    cutPoints.append(coord_offset(cutPoints[-1], (bankWidth, -bankWidth * bermSlope)))
    for x in range(math.ceil(maxCut / batterHeight)):
        cutPoints.append(coord_offset(cutPoints[-1], (batterHeight * batterSides, batterHeight)))
        cutPoints.append(coord_offset(cutPoints[-1], (bermWidth, -bermWidth * bermSlope)))
    # dont need the last horizontal berm
    cutPoints = cutPoints[:-1]
    # finish off to a polygon representing the cut
    cutPoints.append((0, cutPoints[-1][1]))
    cutPoly = Polygon(cutPoints)
    cutSection = cutPoly.union(scale(cutPoly, xfact=-1, origin=origin))

    # create a maximum fill section
    bankTop = coord_offset(channelSide, (bankWidth, -bankWidth * bermSlope))
    bankToe = coord_offset(bankTop, (maxFill * sides, -maxFill))
    fillPoints = [*channelPerimeterPoints, bankTop, bankToe]
    fillPoints.append((0, fillPoints[-1][1]))
    fillPoly = Polygon(fillPoints)
    fillSection = fillPoly.union(scale(fillPoly, xfact=-1, origin=origin))

    # make a large polygon to select the required parts of the sections and cut them with the Natural Surface
    # rather then shift all the parts to the DFL shift the NS in the opposite direction
    NSline = [coord_offset((-maskSize,-maskSize*crossGradient),(0,-DFLwrtNS)),
              coord_offset((maskSize,maskSize*crossGradient),(0,-DFLwrtNS))]
    shallow_cut_line = [coord_offset((-maskSize,-maskSize*crossGradient),(0,-DFLwrtNS-rockDepth)),
                        coord_offset((maskSize,maskSize*crossGradient),(0,-DFLwrtNS-rockDepth))]
    fillMask = Polygon((*NSline,(maskSize,maskSize),(-maskSize,maskSize)))
    cutMask = Polygon((*shallow_cut_line,(maskSize,-maskSize),(-maskSize,-maskSize)))
    shallow_cutMask = Polygon((*NSline,*reversed(shallow_cut_line)))
    maskedFill = fillMask.intersection(fillSection)
    maskedCut = cutMask.intersection(cutSection)
    shallow_maskedCut = shallow_cutMask.intersection(cutSection)

    # if over excevating create sections for this.
    maskedoverExFill = maskedoverExCut = shallow_maskedoverExCut = None
    if overExcavateBed or overExcavateWall:
        bedCentreOverEx = coord_offset(bedCentre,(0,-overExcavateBed))
        bedEdgeOverEx = coord_offset(bedEdge,(overExcavateWall/sides,-overExcavateBed))
        channelSideOverEx = coord_offset(channelSide,(overExcavateWall,-overExcavateWall*bermSlope))
        overExPoints = [*channelPerimeterPoints,channelSideOverEx,bedEdgeOverEx,bedCentreOverEx]
        overExPoly = Polygon(overExPoints)
        overExSection = overExPoly.union(scale(overExPoly, xfact = -1, origin = origin))
        maskedoverExFill = fillMask.intersection(overExSection)
        maskedoverExCut = cutMask.intersection(overExSection)
        shallow_maskedoverExCut = shallow_cutMask.intersection(overExSection)
        maskedFill = maskedFill.difference(maskedoverExFill)

    borrowMask = Polygon(((maskSize, -depthFlow - overExcavateBed),
                          (maskSize, -maskSize),
                          (-maskSize, -maskSize),
                          (-maskSize, -depthFlow - overExcavateBed)))
    maskedBorrowFill = borrowMask.intersection(maskedFill)
    maskedFill = maskedFill.difference(borrowMask)

    # create a water section
    waterSide = coord_offset(bedEdge, ((depthFlow) * sides, depthFlow))
    waterPoints = [bedCentre, bedEdge, waterSide]
    waterPoly = Polygon([origin, *waterPoints])
    waterSection = waterPoly.union(scale(waterPoly, xfact=-1, origin=origin))

    wettedPerimeterLine = LineString(waterPoints)
    wettedPerimeterSection = wettedPerimeterLine.union(scale(wettedPerimeterLine, xfact=-1, origin=origin))
    
    everything = unary_union([shallow_maskedCut,maskedCut,maskedFill,maskedBorrowFill,shallow_maskedoverExCut,maskedoverExCut,maskedoverExFill,waterSection,wettedPerimeterSection])
    full_width = everything.bounds[2] - everything.bounds[0]
    
    return (shallow_maskedCut,maskedCut,maskedFill,maskedBorrowFill,shallow_maskedoverExCut,maskedoverExCut,maskedoverExFill,waterSection,wettedPerimeterSection,full_width)


def plot_channel_geometry(DFLwrtNS=0, crossGradient=0, bedWidth=50, depthFlow=4.4, freeboard=1, sides=2,
    batterSides=1.5, bermSlope=1/20, batterHeight=5, bermWidth=4, bankWidth=5,
    overExcavateBed=1, overExcavateWall=3, rockDepth=5, ax=None):
    """Plot the cross-section returned by ``channel_geometry`` and return its axes."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Polygon as PolygonPatch

    geometries = channel_geometry(
        DFLwrtNS=DFLwrtNS,
        crossGradient=crossGradient,
        bedWidth=bedWidth,
        depthFlow=depthFlow,
        freeboard=freeboard,
        sides=sides,
        batterSides=batterSides,
        bermSlope=bermSlope,
        batterHeight=batterHeight,
        bermWidth=bermWidth,
        bankWidth=bankWidth,
        overExcavateBed=overExcavateBed,
        overExcavateWall=overExcavateWall,
        rockDepth=rockDepth,
    )
    labels = (
        'Soil excavation',
        'Rock excavation',
        'Normal fill',
        'Borrow fill',
        'Soil over-excavation',
        'Rock over-excavation',
        'Over-excavation fill',
        'Water',
        'Wetted perimeter',
    )
    colors = ('#d9a066', '#777777', '#83b96b', '#9a6fc0', '#efb27b',
              '#555555', '#a6cf8b', '#6bb9e8', '#174a6e')
    if ax is None:
        _, ax = plt.subplots()

    def draw_geometry(geometry, color, label):
        if geometry is None or geometry.is_empty:
            return
        if geometry.geom_type == 'Polygon':
            ax.add_patch(PolygonPatch(
                geometry.exterior.coords,
                closed=True,
                facecolor=color,
                edgecolor='black',
                alpha=0.55,
                label=label,
            ))
        elif geometry.geom_type in ('LineString', 'LinearRing'):
            x, y = geometry.xy
            ax.plot(x, y, color=color, linewidth=2, label=label)
        elif hasattr(geometry, 'geoms'):
            for part in geometry.geoms:
                draw_geometry(part, color, label)

    for geometry, color, label in zip(geometries[:9], colors, labels):
        draw_geometry(geometry, color, label)

    full_width = geometries[9]
    ax.set_xlim(-full_width / 2, full_width / 2)
    ax.set_xlabel('Width (m)')
    ax.set_ylabel('Elevation relative to channel bed (m)')
    ax.set_title('Channel cross-section')
    ax.set_aspect('equal', adjustable='datalim')
    handles = [
        Patch(facecolor=color, edgecolor='black', alpha=0.55, label=label)
        for geometry, color, label in zip(geometries[:9], colors, labels)
        if geometry is not None and not geometry.is_empty
    ]
    if handles:
        ax.legend(handles=handles, loc='best')
    ax.grid(True, alpha=0.25)
    return ax

