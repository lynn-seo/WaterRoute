# seline.ng@csiro.au, start: 9/feb/2021, clayton vic
# corresponding Excel sheet: 'Losses'
import numpy as np
import fluids
from waterroute.channel import earthworks
from scipy.optimize import root

import math

def channel_v (flow_area, wetted_perimeter, S= 1/20000, n=0.0225):
    ''' returns flow velocity in m/s
    requires Flow Area [m2] and Wetted perimeter [m]
    S is slope in m/m. default is 1/20000 as the lowest practical slope for a channel.
    n is manning n. default is 0.0225 from Kev for constructed channels.
    '''
    #Hydraulic radius of the channel, Flow Area/Wetted perimeter [m]
    Rh = flow_area/wetted_perimeter
    v = fluids.open_flow.V_Manning(Rh, S, n)
    # assert v<=1 #TODO ask Arthur 28/Aug/2025
    return v

def max_channel_velocity(flow_depth,C = 0.84,upper_v_limit = 1):
    '''There is also a check for non-scouring velocities, to further refine the 1m/sec limit noted above.
        The formula is 
        Vs =0.5CD0.64 where,
        Vs is non scouring limit
        D is depth of flow
        C is a coefficient for soil conditions.  Use 0.84 in the absence of other info.  
        This is for fine, light sandy soil.  
        Note these numbers are based on US research, so relate to sandy soils.  
        Not a whole lot of use in clay profiles, where sodicity is a much better predictor of erodeability.  
        
        also applies an upper limit 
    '''
    Vs =0.5*C*flow_depth*0.64
    return min(Vs,upper_v_limit)

def channel_discharge_and_v (flow_area, wetted_perimeter, S= 1/20000, n=0.0225):
    ''' returns discharge in cumecs
    requires Flow Area [m2] and Wetted perimeter [m]
    S is slope in m/m. default is 1/20000 as the lowest practical slope for a channel.
    n is manning n. default is 0.0225 from Kev for constructed channels.
    '''
    v = channel_v (flow_area, wetted_perimeter, S= S, n=n)
    return (flow_area * v,v)

def scaled_depth(bed_width,max_flow_depth=4.4,depth_scale=None):
    '''
    if depth_scale is provided the bed width is multiplied by this factor to determine the flow_depth. 1/8 might be a useful value.
    in this case max_flow_depth is used to set a maximum water depth.
    If depth_scale is None and relationship between width and depth is used
    '''
    if depth_scale is not None:
        flow_depth = bed_width * depth_scale
    if depth_scale is None:
        # Changes to the width to depth ratio. KEV -  I looked a bit further into this, and while I don’t think this will help too much, I did find a reference on a design chart that gave the following:
        # BW>15m, then BW/d=8
        # 15>BW>7.5, then BW/d=5
        # BW<7.5, then BW/d=4
        #AR. simplified this down to the below. A linear change between the thresholds better then step changes.
        depth_scale = 1/sorted([4, bed_width * 0.533333, 8])[1]
        flow_depth = bed_width*depth_scale

    return min(max_flow_depth,flow_depth)

def channel_discharge_and_v_for_width (bed_width,max_flow_depth,depth_scale, S= 1/20000, n=0.0225):
    '''
    if depth_scale is provided the bed width is multiplied by this factor to determine the flow_depth. 1/8 might be a useful value.
    in this case flow_depth is used to set a maximum water depth.
    IF depth_scale is None (default) then flow_depth is used (4.4 default)
    '''
    #if bed_width <= 0:
    #    return 0
    flow_depth = scaled_depth(bed_width,max_flow_depth,depth_scale)
    #print (bed_width,flow_depth)
    ew = earthworks.channel_geometry(DFLwrtNS = 0, crossGradient = 0, bedWidth = bed_width, depthFlow = flow_depth, freeboard = 1, sides = 2, \
    batterSides = 1.5, bermSlope = 1/20, batterHeight = 5, bermWidth = 4, bankWidth = 5, \
    overExcavateBed = 1, overExcavateWall = 3,rockDepth = 5)
    (shallow_maskedCut,maskedCut,maskedFill,maskedBorrowFill,
     shallow_maskedoverExCut,maskedoverExCut,maskedoverExFill,
     waterSection,wettedPerimeterSection,full_width) = ew
    dis,v = channel_discharge_and_v(waterSection.area,wettedPerimeterSection.length, S= S, n=n)
    return dis,v,flow_depth

def channel_width_for_discharge(discharge,max_flow_depth,depth_scale, S= 1/20000, n=0.0225):
    f =  lambda bed_width: channel_discharge_and_v_for_width(bed_width,max_flow_depth,depth_scale, S= S, n=n)[0]-discharge
    r = (root(f,1))
    if not r.success:
        raise RuntimeError (r.message)
    return r.x[0]

def exploring_non_scouring_limits(bed_width,max_flow_depth,depth_scale, S= 1/20000, n=0.0225):
    '''end results was unless channel is steeper then 1:13000 or narower then 7m these velocity limits dont come into play.'''
    dis2,v,flow_depth = channel_discharge_and_v_for_width(bed_width,max_flow_depth,depth_scale, S= S, n=n)
    print (dis2,v,flow_depth)
    max_v = max_channel_velocity(flow_depth)
    print (max_v)
    if v>max_v:
        for d in range(int(flow_depth*10), 0, -1):
            dis2,v,flow_depth = channel_discharge_and_v_for_width(bed_width,d/10,depth_scale, S= S, n=n)
            max_v = max_channel_velocity(flow_depth)
            print (dis2,v,flow_depth,max_v)

            if v<= max_v:
                break
    print ('')

def GL_annum_to_cumecs(GL,num_days):
    return GL * 1E9 / num_days /24/60/60/1000
def ML_day_to_cumecs(ML):
    return ML * 1E6 /24/60/60/1000
def cumecs_to_GL_day (cumecs):
    return cumecs*1000*60*60*24/1e9

#used for some testing
#def channel_discharge_for_width2 (bed_width,flow_depth = 4.4,depth_scale = None):
#    '''
#    if depth_scale is provided the bed width is multiplied by this factor to determine the flow_depth. 1/8 might be a useful value.
#    in this case flow_depth is used to set a maximum water depth.
#    IF depth_scale is None (default) then flow_depth is used (4.4 default)
#    '''
#    flow_depth = scaled_depth(bed_width,flow_depth,depth_scale)#

#    ew = earthworks.channel_geometry(DFLwrtNS = 0, crossGradient = 0, bedWidth = bed_width, depthFlow = flow_depth, freeboard = 1, sides = 2, \
#    batterSides = 1.5, bermSlope = 1/20, batterHeight = 5, bermWidth = 4, bankWidth = 5, \
#    overExcavateBed = 1, overExcavateWall = 3)
#    (shallow_maskedCut,maskedCut,maskedFill,maskedBorrowFill,
#    shallow_maskedoverExCut,maskedoverExCut,maskedoverExFill,
#    waterSection,wettedPerimeterSection,full_width) = ew
#    dis = channel_discharge(waterSection.area,wettedPerimeterSection.length)

#    S= 1/20000
#    n=0.0225
#    v = channel_v (waterSection.area,wettedPerimeterSection.length)
#    print (f'discharge {dis:.2f}, velocity {v:.2f}, flow_area {waterSection.area:.2f}, wetted_perimeter {wettedPerimeterSection.length:.2f}, S {S:.5f}, n {n} ')



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# INDIVIDUAL FUNCTIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def getLossEvap(evapPanAnnual, evapPeakDaily, channelLength, surfaceWidth):
    '''
    ------------
    Inputs:
    ------------
    evapPanAnnual: m/annum
    evapPeakDaily: mm/day
    channelLength: km
    surfaceWidth: m
    ------------
    Outputs:
    ------------
    lossEvap['Annual']: ML/annum
    lossEvap['Peak daily']: ML/day
    '''
    lossEvap = {}
    lossEvap['Annual'] = evapPanAnnual * channelLength * surfaceWidth #ML/annum
    lossEvap['Peak daily'] = evapPeakDaily * channelLength * surfaceWidth / 1000 #ML/daily
    return lossEvap


def getLossSeepage(waterDepth, bedPerimeter, channelLength, k, i, annualDuration):
    '''
    ------------
    Inputs:
    ------------
    bedPerimeter: m - use bed width as a proxy for perimeter.  The actual perimeter has the side slopes, but the "i" is a lot less there, so for all practical purposes, they can be ignored, and BW can be used.
    channelLength: km 
    k: m/day 
    i: unitless -  "i" is gradient, expressed as water depth over lining thickness.  So, given that we are using 1m lining, I suppose you can think of it as water depth.
    annualDuration: number of days
    waterDepth : m
    ------------
    Outputs:
    ------------
    lossSeepage['Daily']: ML/day
    lossSeepage['Annual']: ML/annum 
    '''
    lossSeepage = {}
    area = bedPerimeter * channelLength * 1000 #m2
    i = 3 #fix gradient as 3
    lossSeepage['Daily'] = waterDepth * k * i * area / 1000 #ML/day #we incorporate the water depth16/03/2026
    lossSeepage['Annual'] = lossSeepage['Daily'] * annualDuration #ML/annum
    return lossSeepage


def getLossStockAndDomestic(lossEachStockDomestic, spacingStockDomestic, channelLength, annualDuration):
    '''
    ------------
    Inputs:
    ------------
    lossEachStockDomestic: ML/annum
    spacingStockDomestic: km
    channelLength: km
    annualDuration: number of days
    ------------
    Outputs:
    ------------
    lossStockDomestic['Annual']: ML/annum
    lossStockDomestic['Daily']: ML/day
    '''
    lossStockDomestic = {}
    lossStockDomestic['Annual'] = channelLength / spacingStockDomestic * lossEachStockDomestic #ML/annnum
    lossStockDomestic['Daily'] = lossStockDomestic['Annual'] / annualDuration #ML/day
    return lossStockDomestic


def getLossLeakage(leakageRate, channelLength, annualDuration):
    '''
    ------------
    Inputs:
    ------------
    leakageRate: L/sec/km
    channelLength: km
    annualDuration: number of days
    ------------
    Outputs:
    ------------
    lossLeakage['Annual']: ML/annum
    lossLeakage['Daily']: ML/day
    '''
    lossLeakage = {}
    leakage = leakageRate * channelLength #L/sec
    lossLeakage['Daily'] = leakage * 60 * 60 * 24 /1000000 #ML/day
    lossLeakage['Annual'] = lossLeakage['Daily'] * annualDuration #ML/annum
    return lossLeakage


def getLossMeteringInaccuracy(diversionFromDam, diversionEfficiency):
    '''
    ------------
    Inputs:
    ------------
    diversionFromDam: GL/annum
    diversionEfficiency: %
    ------------
    Output:
    ------------
    lossMeteringInaccuracy['Annual']: ML/annum
    '''
    # standard AS4747 says + / - 2.5 % lab, and 5 % field
    # allow 1.5% of gross diversions
    lossMeteringInaccuracy = {}
    lossMeteringInaccuracy['Annual'] = diversionFromDam * (diversionEfficiency / 100) * 0.015 * 1000 #ML/annum
    return lossMeteringInaccuracy


def getLossInitialAnnualRefillAndDrain(channelLength, channelXSection, rolling_shutdown_factor=1):
    '''
    ------------
    Inputs:
    ------------
    channelLength: km
    channelXSection: m2
    rolling_shutdown_factor: default of 1. proportion of the channel that is drained during shutdown
    ------------
    Output:
    ------------
    lossInitialRefillAndDrain['Annual']: ML/annum
    '''
    # assume pipe remains full
    lossInitialRefillAndDrain = {}
    lossInitialRefillAndDrain['Annual'] = rolling_shutdown_factor * channelLength * channelXSection * 1000 / 1000 #ML/annum
    return lossInitialRefillAndDrain


def getLossTheft():
    '''
    ------------
    Outputs:
    ------------
    lossTheft['Annual']: 0 (negligible)
    lossTheft['Daily']: 0 (negligible)
    '''
    # If we go with the assumption of TCC, we can assume this is negligible
    lossTheft = {}
    lossTheft['Annual'] = 0
    lossTheft['Daily'] = 0
    return lossTheft


def getLossOperationalOverflow(diversionFromDam, distributionEfficiency):
    '''
    ------------
    Inputs:
    ------------
    diversionFromDam: GL/annum
    distributionEfficiency: %
    ------------
    Output:
    ------------
    lossOperationalOverflow['Annual']: ML/annum
    '''
    # This is a tricky one to estimate, as it will be low with TCC; however, some is inevitable
    # In the absence of other info, we will assume 2% of the total flow (based purely on the fact that this has been measured at up to 10% in manual systems)
    lossOperationalOverflow = {}
    lossOperationalOverflow['Annual'] = diversionFromDam * (distributionEfficiency / 100) * 0.02 * 1000 #ML/annum
    return lossOperationalOverflow


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FINAL FUNCTIONS WRAPPED AROUND INDIVIDUAL FUNCTIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def getPeakDailyLossesLinear(annualDuration = 330, channelLength = 426, surfaceWidth = 68, evapPanAnnual = 2, \
    evapPeakDaily = 10, bedPerimeter = 50, k = 0.002, i = 4.4, lossEachStockDomestic = 2, spacingStockDomestic = 2, leakageRate = 2, waterDepth = 3):
    '''
    ------------
    Inputs:
    ------------
    annualDuration: number of days
    channelLength: km
    surfaceWidth: m
    evapPanAnnual: m/annum
    evapPeakDaily: mm/day
    bedPerimeter: m
    k: m/day
    i: unitless
    lossEachStockDomestic: ML/annum
    spacingStockDomestic: km
    leakageRate: L/sec/km
    ------------
    Outputs:
    ------------
    peakDailyLosses['Evap']: ML/day
    peakDailyLosses['Seepage']: ML/day
    peakDailyLosses['Stock and domestic']: ML/day
    peakDailyLosses['Leakage']: ML/day
    peakDailyLosses['Metering inaccuracy']: ML/day
    peakDailyLosses['Total']: ML/day
    '''
    peakDailyLosses = {}
    peakDailyLosses['Evap'] = getLossEvap(evapPanAnnual, evapPeakDaily, channelLength, surfaceWidth)['Peak daily']
    peakDailyLosses['Seepage'] = getLossSeepage(waterDepth, bedPerimeter, channelLength, k, i, annualDuration)['Daily']
    peakDailyLosses['Stock and domestic'] = getLossStockAndDomestic(lossEachStockDomestic, spacingStockDomestic, channelLength, annualDuration)['Daily']
    peakDailyLosses['Leakage'] = getLossLeakage(leakageRate, channelLength, annualDuration)['Daily']
    peakDailyLosses['Metering inaccuracy'] = 1.5 / 100 * sum(peakDailyLosses.values())
    peakDailyLosses['Total'] = sum(peakDailyLosses.values())
    return peakDailyLosses

def getAnnualLossesDamToChannel(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87):
    '''
    ------------
    Inputs:
    ------------
    diversionFromDam: GL/annum
    diversionEfficiency: %
    distributionEfficiency: %
    ------------
    Outputs:
    ------------
    annualLosses['Metering inaccuracy']: ML/annum
    annualLosses['Operational overflow']: ML/annum
    annualLosses['Total']: ML/annum
    '''
    annualLosses = {}
    annualLosses['Metering inaccuracy'] = getLossMeteringInaccuracy(diversionFromDam, diversionEfficiency)['Annual']
    annualLosses['Operational overflow'] = getLossOperationalOverflow(diversionFromDam, distributionEfficiency)['Annual']
    annualLosses['Total'] = sum(annualLosses.values())
    return annualLosses


def getAnnualLossesLinear(annualDuration = 330, channelLength = 426, surfaceWidth = 68, evapPanAnnual = 2, \
    evapPeakDaily = 10, bedPerimeter = 50, k = 0.002, i = 4.4, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
    leakageRate = 2, channelXSection = 258.72,rolling_shutdown_factor=1, waterDepth= 3):
    '''
    ------------
    Inputs:
    ------------
    annualDuration: number of days
    channelLength: km
    surfaceWidth: m
    evapPanAnnual: m/annum
    evapPeakDaily: mm/day
    bedPerimeter: m
    k: m/day
    i: unitless
    lossEachStockDomestic: ML/annum
    spacingStockDomestic: km
    leakageRate: L/sec/km
    channelXSection: m2
    waterDepth : m used for seepage loss calculation, 16/03/2026
    ------------
    Outputs:
    ------------
    annualLosses['Evap']: ML/annum
    annualLosses['Seepage']: ML/annum
    annualLosses['Stock and domestic']: ML/annum
    annualLosses['Leakage']: ML/annum
    annualLosses['Initial refill and drain']: ML/annum
    annualLosses['Theft']: ML/annum
    annualLosses['Total']: ML/annum
    '''
    annualLosses = {}
    annualLosses['Evap'] = getLossEvap(evapPanAnnual, evapPeakDaily, channelLength, surfaceWidth)['Annual']
    annualLosses['Seepage'] = getLossSeepage(waterDepth, bedPerimeter, channelLength, k, i, annualDuration)['Annual']
    annualLosses['Stock and domestic'] = getLossStockAndDomestic(lossEachStockDomestic, spacingStockDomestic, channelLength, annualDuration)['Annual']
    annualLosses['Leakage'] = getLossLeakage(leakageRate, channelLength, annualDuration)['Annual']
    annualLosses['Initial refill and drain'] = getLossInitialAnnualRefillAndDrain(channelLength, channelXSection,rolling_shutdown_factor)['Annual']
    annualLosses['Theft'] = getLossTheft()['Annual']
    annualLosses['Total'] = sum(annualLosses.values())
    return annualLosses


# ++++++++++++++++++++++++
# ++++++++++++++++++++++++
# EXAMPLE APPLICATIONS
# ++++++++++++++++++++++++
# ++++++++++++++++++++++++

if __name__ == "__main__":

    # Losses as far as Buchanan
    # -----------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 426, surfaceWidth = 68, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 50, k = 0.002, i = 4.4, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 258.72)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamA = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearA = getAnnualLossesLinear(**settings2)
    print(annualLossesLinearA, annualLossesDamA, annualLossesLinearA['Total'] + annualLossesDamA['Total'])

    peakDailyLossesLinearA = getPeakDailyLossesLinear(**settings3)
    print(peakDailyLossesLinearA)
    print()


    # Losses for Buchanan to Big Drop
    # ------------------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 69.301, surfaceWidth = 60.2, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 45, k = 0.002, i = 4.4, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 199.88)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamB = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearB = getAnnualLossesLinear(**settings2)
    for i in annualLossesLinearB.keys(): annualLossesLinearB[i] += annualLossesLinearA[i]
    print(annualLossesLinearB, annualLossesDamB, annualLossesLinearB['Total'] + annualLossesDamB['Total'])

    peakDailyLossesLinearB = getPeakDailyLossesLinear(**settings3)
    for i in peakDailyLossesLinearB.keys(): peakDailyLossesLinearB[i] += peakDailyLossesLinearA[i]
    print(peakDailyLossesLinearB)
    print()


    # Losses as far as Section 3-Big Drop to Cornish Ck
    # -----------------------------------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 49, surfaceWidth = 64.8, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 55, k = 0.002, i = 2.45, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 146.76)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamC = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearC = getAnnualLossesLinear(**settings2)
    for i in annualLossesLinearC.keys(): annualLossesLinearC[i] += annualLossesLinearB[i]
    print(annualLossesLinearC, annualLossesDamC, annualLossesLinearC['Total'] + annualLossesDamC['Total'])

    peakDailyLossesLinearC = getPeakDailyLossesLinear(**settings3)
    for i in peakDailyLossesLinearC.keys(): peakDailyLossesLinearC[i] += peakDailyLossesLinearB[i]
    print(peakDailyLossesLinearC)
    print()


    # Losses as far as Section 4-to Aramac Ck
    # ------------------------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 121, surfaceWidth = 53, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 40, k = 0.002, i = 3.25, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 151.13)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamD = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearD = getAnnualLossesLinear(**settings2)
    for i in annualLossesLinearD.keys(): annualLossesLinearD[i] += annualLossesLinearC[i]
    print(annualLossesLinearD, annualLossesDamD, annualLossesLinearD['Total'] + annualLossesDamD['Total'])

    peakDailyLossesLinearD = getPeakDailyLossesLinear(**settings3)
    for i in peakDailyLossesLinearD.keys(): peakDailyLossesLinearD[i] += peakDailyLossesLinearC[i]
    print(peakDailyLossesLinearD)


    #  TODO
    #  Losses for Elliot, section1
    # ------------------------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 70, surfaceWidth = 53, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 40, k = 0.002, i = 3.25, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 151.13)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamD = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearD = getAnnualLossesLinear(**settings2)
    for i in annualLossesLinearD.keys(): annualLossesLinearD[i] += annualLossesLinearC[i]
    print(annualLossesLinearD, annualLossesDamD, annualLossesLinearD['Total'] + annualLossesDamD['Total'])

    peakDailyLossesLinearD = getPeakDailyLossesLinear(**settings3)
    for i in peakDailyLossesLinearD.keys(): peakDailyLossesLinearD[i] += peakDailyLossesLinearC[i]
    print(peakDailyLossesLinearD)

    # Losses for Elliot, section2, to update
    # ------------------------------------------
    settings1 = dict(diversionFromDam = 2520, diversionEfficiency = 87, distributionEfficiency = 87)

    settings2 = dict(annualDuration = 330, channelLength = 70, surfaceWidth = 53, evapPanAnnual = 2, \
        evapPeakDaily = 10, bedPerimeter = 40, k = 0.002, i = 3.25, lossEachStockDomestic = 2, spacingStockDomestic = 2, \
        leakageRate = 2, channelXSection = 151.13)

    settings3 = settings2.copy()
    settings3.pop('channelXSection')

    annualLossesDamD = getAnnualLossesDamToChannel(**settings1)
    annualLossesLinearD = getAnnualLossesLinear(**settings2)
    for i in annualLossesLinearD.keys(): annualLossesLinearD[i] += annualLossesLinearC[i]
    print(annualLossesLinearD, annualLossesDamD, annualLossesLinearD['Total'] + annualLossesDamD['Total'])

    peakDailyLossesLinearD = getPeakDailyLossesLinear(**settings3)
    for i in peakDailyLossesLinearD.keys(): peakDailyLossesLinearD[i] += peakDailyLossesLinearC[i]
    print(peakDailyLossesLinearD)

