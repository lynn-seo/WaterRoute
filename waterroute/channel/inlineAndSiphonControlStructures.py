# seline.ng@csiro.au, start: 21/dec/2020, clayton vic
# corresponding Excel sheets: 'Siphons', 'Control Structures'
import math
import numpy as np
from waterroute import piperoute

def siphon_inletOutlet_loss(V):
    '''calculates head loss at inlet and outlet of siphon
    V: velocity in m/s
    returns head loss in metres'''
    
    lossInletoutlet = 1.2 * V**2 / (2 * 9.81)
    return lossInletoutlet


def siphon_headLoss(L, n, V, D, n_barrels=1):
    '''calculates head loss in a siphon
    
    L: length of siphon in metres
    n: manning coefficient
    V: velocity in m/s
    D: diameter in metres
    
    returns head loss in metres'''
    #https://www.calculatoratoz.com/en/head-loss-by-manning-formula-calculator/Calc-25959
    ''''''
    headLoss = n_barrels*(L*(n*V)**2)/(0.157*D**(4/3 ))
    return headLoss

def siphon_headLoss_v2(L, n, Q, D, n_barrels=1):
    '''calculates head loss in a siphon
    
    L: length of siphon in metres
    n: manning coefficient
    Q: flow in m^3/s
    D: diameter in metres
    
    returns head loss in metres'''
    #https://www.calculatoratoz.com/en/head-loss-by-manning-formula-calculator/Calc-25959
    ''''''
    headLoss = 16*L*n**2*Q**2/(0.157* math.pi**2 * D**(16/3) * n_barrels**2) 
    return headLoss


# structureType, 'INLINE', 'SIPHON-INLET', 'SIPHON-OUTLET'
def getCostSteelworkSansControlGates(structureType, noPiers, t, spacingBarrels, D):
    from math import ceil

    # check inputs
    # ---------------
    if structureType.lower() not in ['inline', 'siphon-inlet', 'siphon-outlet']:
        raise ValueError('Invalid structureType. structureType must be: "INLINE", "SIPHON-INLET" or "SIPHON-OUTLET".')

    li = [noPiers, t, spacingBarrels, D]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # set number of units of each item
    # ---------------------------------
    # For inline control structures, you may want to take a gate out of service without shutting down the channel.
    # So each bay has facility for a bulkhead each side of each gate. The situation with siphons is a bit
    # different. Since you can dewater a barrel, the normal practice would be to only have one set of guides per barrel
    # upstream on the inlet, and one set downstream on the outlet.
    noBulkheadGuides = noPiers + 1
    if structureType.lower() == 'inline': noBulkheadGuides = (noPiers + 1) * 2

    # The number of bulkheads is an operational question. Operators like to have bulkheads for all occasions,
    # bean counters only want one for the entire scheme. If we have 4 for a 8 barrel structure, I think you can
    # safely proportion for smaller no.of barrels. Don’t go above 4 though.
    noBulkheads = min(4, ceil(0.5 * (noPiers+1)))
    if structureType.lower() == 'siphon-outlet': noBulkheads = ceil(0.5 * noBulkheads)

    # This covers the protective handrails around the structure to meet code requirements, and the necessary
    # walkways to get to the inner gates and bulkheads. We probably should have a proportioning for smaller barrel
    # diameters as well. It is not a major item, and this number is just experience based.
    noWalkwaysHandrails = noPiers + 1

    if structureType.lower() in ['inline', 'siphon-outlet']: noFloatingWeedDeflector = 0
    if structureType.lower() == 'siphon-inlet': noFloatingWeedDeflector = 1

    # find individual and total costs
    # ----------------------------------------
    # bulkhead guides, unit cost 3000 $/unit * scaling factor
    res = {}
    scale = (2 * D + spacingBarrels - t) / (2 * 5.4 + 5 - 0.4)
    res['Bulkhead guides'] = 3000 * scale * noBulkheadGuides

    # bulkheads, unit cost 5000 $/unit * scaling factor
    scale = D * (spacingBarrels - t) / (5.4 * (5 - 0.4))
    res['Bulkheads'] = 5000 * scale * noBulkheads

    # walkways and handrails, unit cost 6000 $/unit * scaling factor
    scale = spacingBarrels * (noPiers + 1) / (5 * (7 + 1))
    res['Walkways and handrails'] = 6000 * scale * noWalkwaysHandrails

    # floating weed deflector, unit cost 30000 $/unit
    res['Floating weed deflector'] = 30000 * noFloatingWeedDeflector

    # total cost
    res['Total'] = res['Bulkhead guides'] + res['Bulkheads'] + res['Walkways and handrails'] + res['Floating weed deflector']
    return res


def getCostControlGates(waterDepth, noPiers, num_barrels, pipeOD):
    # check inputs
    li = [waterDepth, noPiers]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # control (overshot) gates
    # as per Rubicon quote of 5/11/2020
    # scaled for area of gate
    noControlGates = noPiers + 1
    #replace 4 with pipe diamtre 
    #replace noControlGates with no barrels
    # cost = ((200000 * (waterDepth * 4) / (4.4 * 4)) + 8225.58 + 19372) * noControlGates
    cost = ((200000 * (waterDepth * pipeOD) / (4.4 * 4)) + 8225.58 + 19372) * num_barrels
    return cost


# lengthCarrWalls = L3
def getCostDownstreamRockwork(bedWidth, D, lengthCarrWalls, rockworkBeyondCarrWalls):
    # check inputs
    li = [bedWidth, D, lengthCarrWalls, rockworkBeyondCarrWalls]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # bed area t=0.5, d/s, rate 100 $/m2
    res = {}
    res['Bed area t=0.5'] = bedWidth * (lengthCarrWalls + rockworkBeyondCarrWalls) * 100

    # batters, d/s, rate 100 $/m2, assumes a 1:2 gradient
    res['Batters'] = (lengthCarrWalls + rockworkBeyondCarrWalls) * (D * 5 ** 0.5) * 2 * 100

    # total cost
    res['Total'] = res['Bed area t=0.5'] + res['Batters']
    return res


def getCostUpstreamRockwork(bedWidth, D, rockworkUpstream):
    # check inputs
    li = [bedWidth, D, rockworkUpstream]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # bed area t=0.5, u/s, rate 100 $/m2
    res = {}
    res['Bed area t=0.5'] = bedWidth * rockworkUpstream * 100

    # batters, u/s, rate 100 $/m2, assumes a 1:2 gradient
    res['Batters'] = rockworkUpstream * (D * 5**0.5) * 2 * 100

    # total cost
    res['Total'] = res['Bed area t=0.5'] + res['Batters']
    return res


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# --------------------------------------------------------------
# INDIVIDUAL FUNCTIONS TO COMPUTE CONCRETE VOLUMES AND COSTS
# --------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# lengthGateSection = L1, lengthDissipater = L2, lengthCarrWalls = L3
def getCostConcrete_inlineControlStructure(lengthGateSection, lengthDissipater, lengthCarrWalls, D, H, \
    cutoffMain, cutoffT, t, spacingBarrels, noPiers, frontCutoff, concRate):

    # check inputs
    li = [lengthGateSection, lengthDissipater, lengthCarrWalls, D, H, cutoffMain, cutoffT, t, spacingBarrels, noPiers, frontCutoff, concRate]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # get costs
    res = {}
    number = 2
    res['Side walls'] = (2 * D * D / 2 + D * lengthGateSection + lengthDissipater * (D + H) + (D + H) * lengthCarrWalls) * t * number * concRate

    number = noPiers
    res['Piers'] = (lengthGateSection + 1) * D * t * number * concRate

    number = 1
    res['Floor'] = (2 * D + lengthGateSection + lengthDissipater) * ((noPiers + 1) * spacingBarrels) * t * number * concRate

    number = 1
    res['Front cutoff'] = frontCutoff * (noPiers + 1) * spacingBarrels * cutoffT * number * concRate

    number = 2
    res['Gate cutoff, sides'] = (D+H) * cutoffMain * cutoffT * number * concRate

    number = 1
    res['Gate cutoff, base'] = (noPiers+1) * spacingBarrels * cutoffMain * cutoffT * number * concRate

    number = 1
    res['Vertical section of floor'] = (noPiers+1) * spacingBarrels * H * t * number * concRate

    number = 2
    res['Main wall side cutoff at jump'] = (D+H) * cutoffMain * cutoffT * number * concRate

    res['Total'] = res['Side walls'] + res['Piers'] + res['Floor'] + res['Front cutoff'] + res['Gate cutoff, sides']  \
        + res['Gate cutoff, base'] + res['Vertical section of floor'] + res['Main wall side cutoff at jump']

    return res


# lengthMainSection = L1, lengthInletSection = L2
def getCostConcrete_siphonInlet(lengthMainSection, lengthInletSection, \
    D, t, cutoffMain, noPiers, spacingBarrels, cutoffT, pipeOD, concRate):

    # check inputs
    li = [lengthMainSection, lengthInletSection, D, t, cutoffMain, noPiers, spacingBarrels, cutoffT, pipeOD, concRate]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # get costs
    res = {}
    number = 2
    res['Side walls'] = D * (lengthMainSection + lengthInletSection/2) * t * number * concRate

    number = noPiers
    res['Piers'] = D * (lengthMainSection+0.5) * t * number * concRate

    number = 1
    res['Floor'] = (lengthMainSection+lengthInletSection) * ((noPiers+1)*spacingBarrels) * t * number * concRate

    number = 1
    res['Front cutoff'] = cutoffMain * (noPiers+1) * spacingBarrels * cutoffT * number * concRate

    number = 1
    res['Back wall, main'] = ((noPiers+1)*spacingBarrels*D - (noPiers+1)*3.14*pipeOD**2/4) * t * number * concRate

    number = 1
    res['Main wall, bottom cutoff'] = (noPiers+1)*spacingBarrels * cutoffMain * cutoffT * number * concRate

    number = 2
    res['Main wall, side cutoff'] = (D + cutoffMain) * cutoffMain * cutoffT * number * concRate

    res['Total'] = res['Side walls'] + res['Piers'] + res['Floor'] + res['Front cutoff'] \
        + res['Back wall, main'] + res['Main wall, bottom cutoff'] + res['Main wall, side cutoff']

    return res


# lengthFirstSection = L1, lengthMiddleSection = L2, lengthCarrWalls = L3
def getCostConcrete_siphonOutlet(lengthFirstSection, lengthMidSection, lengthCarrWalls, \
    D, H, t, cutoffMain, cutoffT, noPiers, spacingBarrels, pipeOD, concRate):

    # --- ADD THESE 3 LINES to prevent negative geometry crashes ---
    lengthFirstSection = max(0, lengthFirstSection)
    lengthMidSection = max(0, lengthMidSection)
    lengthCarrWalls = max(0, lengthCarrWalls)
    
    #check inputs
    li = [lengthFirstSection, lengthMidSection, lengthCarrWalls, D, H, t, cutoffMain, cutoffT, noPiers, spacingBarrels, pipeOD, concRate]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')
   
    # get costs
    res = {}
    number = 2
    res['Side walls'] = (D+H) * (lengthFirstSection + lengthMidSection + lengthCarrWalls) * t * number * concRate

    number = noPiers
    res['Piers'] = (D+H) * lengthFirstSection * t * number * concRate

    number = 1
    res['Floor'] = (lengthFirstSection + lengthMidSection) * ((noPiers + 1) * spacingBarrels) * t * number * concRate

    number = 1
    res['Front cutoff'] = H * (noPiers+1) * spacingBarrels * cutoffT * number * concRate

    number = 1
    res['Back wall, main'] = ((noPiers+1)*spacingBarrels*(D+H) - (noPiers+1)*3.14*pipeOD**2/4) * t * number * concRate

    number = 1
    res['Main wall, bottom cutoff'] = (noPiers+1)*spacingBarrels * cutoffMain * cutoffT * number * concRate

    number = 2
    res['Main wall, side cutoff'] = (D + H + cutoffMain) * cutoffMain * cutoffT * number * concRate

    res['Total'] = res['Side walls'] + res['Piers'] + res['Floor'] + res['Front cutoff'] \
        + res['Back wall, main'] + res['Main wall, bottom cutoff'] + res['Main wall, side cutoff']

    return res


# lengthFirstSection = L1, lengthMiddleSection = L2, lengthCarrWalls = L3
def getCostConcrete_siphonBaffledOutlet(lengthFirstSection, lengthMidSection, lengthCarrWalls, \
    Dup, H, t, pipeOD, spacingBarrels, Ddown, noPiers, lipWidth, concRate):

    # check inputs
    li = [lengthFirstSection, lengthMidSection, lengthCarrWalls, Dup]
    li += [H, t, pipeOD, spacingBarrels, Ddown, noPiers, lipWidth, concRate]
    if any(l < 0 for l in li): raise ValueError('Invalid negative input or inputs.')

    # get costs
    res = {}
    number = 2
    res['Side walls'] = Dup * (lengthFirstSection + lengthMidSection) #first and middle sections
    res['Side walls'] += -1 * (Dup - H - Ddown) * (2 * lengthMidSection - 1.5 * (Dup - H - Ddown)) / 2 #first and middle sections
    res['Side walls'] += lengthCarrWalls * (Ddown + H) #Carr walls
    res['Side walls'] = res['Side walls'] * t * number * concRate

    number = noPiers
    res['Piers'] = Dup * (lengthFirstSection + lengthMidSection)  # first and middle sections
    res['Piers'] += -1 * (Dup - H - Ddown) * (2 * lengthMidSection - 1.5 * (Dup - H - Ddown)) / 2 # first and middle sections
    res['Piers'] = res['Piers'] * t * number * concRate

    number = 1
    res['Floor'] = (lengthFirstSection + lengthMidSection) * spacingBarrels * (noPiers + 1) * t * number * concRate

    number = 1
    res['Front cutoff'] = H * spacingBarrels * (noPiers + 1) * t * number * concRate

    number = 1
    res['Main back wall'] = Dup * spacingBarrels * (noPiers + 1) - (noPiers + 1) * 3.14 * pipeOD ** 2 / 4
    res['Main back wall'] = res['Main back wall'] * t * number * concRate

    number = 2
    res['Top lip, side'] = lengthFirstSection * lipWidth * t * number * concRate

    number = 1
    res['Back lip'] = lipWidth * spacingBarrels * (noPiers + 1) * t * number * concRate

    res['Total'] = res['Side walls'] + res['Piers'] + res['Floor'] + res['Front cutoff'] \
        + res['Main back wall'] + res['Top lip, side'] + res['Back lip']

    return res


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# -----------------------------------------------------------
# FINAL FUNCTIONS WRAPPING INDIVIDUAL FUNCTIONS
# -----------------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def select_siphon(bedWidth,waterDepth,freeboard,capacity_m3_s,siphon_length):
    '''
    returns a tuple of direct_cost,pipe_ID,num_barrels,V
    '''
    max_ID = piperoute.standard_GRP_Ds.max()
    siphon_results = []
    for num_barrels in range(1,300):
        barrel_Q = capacity_m3_s/num_barrels
        #assume 1:2000 slope
        D = piperoute.diameter_for_Q(1,barrel_Q,2000)
        if D < max_ID:
            #pick the next standard D above
            actual_D = piperoute.standard_GRP_Ds.iloc[np.searchsorted(piperoute.standard_GRP_Ds,D)]
            actual_OD = piperoute.standard_GRP_Ds.index[np.searchsorted(piperoute.standard_GRP_Ds,D)]/1000
            direct_cost = sum(piperoute.pipe_diameter_cost(actual_D,siphon_length,num_barrels))
            V = piperoute.volume_to_velocity(barrel_Q,D=actual_D)
            direct_cost += simplified_siphon_structures_cost(bedWidth,waterDepth,freeboard,num_barrels,actual_OD)
            siphon_results.append([direct_cost,actual_D,num_barrels,V])
            #print (siphon_results[-1])
        if len(siphon_results) >= 5:
            #probably only need the smallest number of barrels but lets test a few
            break
    return (min(siphon_results))

def simplified_siphon_structures_cost(bedWidth,waterDepth,freeboard,pipe_OD,num_barrels):
    D = waterDepth + freeboard
    gate_height = waterDepth * 1.2
    settings = dict(bedWidth = bedWidth, waterDepth = waterDepth, lengthMainSection = gate_height*2, D = D , 
    cutoffMain = 1, cutoffT = 0.3, t = D/12, spacingBarrels = pipe_OD + 0.6, noPiers = num_barrels-1, pipeOD = pipe_OD, 
    rockworkUpstream = 10, concRate = 1200)

    inlet_cost = getTotalCost_siphonInlet(**settings)
    
    # settings = dict(bedWidth = bedWidth, D = D, H = 1, lengthCarrWalls = waterDepth*1.5, 
    # t = D/12, cutoffMain = 1.75, cutoffT = 0.3, noPiers = num_barrels-1, spacingBarrels = pipe_OD + 0.6, 
    # pipeOD = pipe_OD, rockworkBeyondCarrWalls=15, concRate=1200)
    #Lynn changed noPiers to max(0,num_barrels-1) not to assign a negative value 
    settings = dict(bedWidth = bedWidth, D = D, H = 1, lengthCarrWalls = waterDepth*1.5, 
    t = D/12, cutoffMain = 1.75, cutoffT = 0.3, noPiers = max(0,num_barrels-1), spacingBarrels = pipe_OD + 0.6, 
    pipeOD = pipe_OD, rockworkBeyondCarrWalls=15, concRate=1200)

    outlet_cost = getTotalCost_siphonOutlet(**settings)
    return inlet_cost['Total'] + outlet_cost['Total']

def simplified_control_structures_cost(bedWidth,waterDepth,freeboard, pipe_OD=None, num_barrels=None):

    height_of_gate = 1.2*waterDepth
    D = waterDepth + freeboard
    t = max(0.2,D/12)
    
    #a couple of guessses while waiting for Kev
    pipe_OD = 0.9*waterDepth #ie  roughly square gates
    num_barrels = math.ceil(bedWidth/(pipe_OD + t)*7/10)
    

    settings = dict(bedWidth = bedWidth, 
                    waterDepth = waterDepth, 
                    lengthGateSection = 1.1*height_of_gate, 
                    lengthDissipater = 3*waterDepth, 
                    lengthCarrWalls = 1.5*waterDepth, 
                    D = D, 
                    H = pipe_OD/2,
                    cutoffMain = 2,
                    cutoffT = 0.3,
                    t = t,
                    spacingBarrels = pipe_OD + t,
                    noPiers = num_barrels-1,
                    frontCutoff = 1,
                    rockworkBeyondCarrWalls = 15,
                    rockworkUpstream = 10,
                    concRate = 1200, 
                    num_barrels = num_barrels,
                    pipeOD = pipe_OD)

    cost = getTotalCost_inlineControlStructure(**settings)
    return cost['Total']

def getTotalCost_inlineControlStructure(bedWidth = 50, waterDepth = 4.4, \
    lengthGateSection = 5, lengthDissipater = 14, lengthCarrWalls = 7, \
    D = 5.4, H = 2, cutoffMain = 2, cutoffT = 0.3, t = 0.4, spacingBarrels = 5, noPiers = 7, frontCutoff = 1, \
    rockworkBeyondCarrWalls = 15, rockworkUpstream = 10, concRate = 1200, num_barrels = 7, pipeOD = 4):
    '''
    bedWidth: Width of the channel bed, m
    waterDepth: Water depth, m
    lengthGateSection: L1, m, this has to be longer than the gate length, which will be about 1.2 * waterDepth
    lengthDissipater: L2, m, this has to be about 3 times the incoming water depth, since it must be fully contain the resulting hydraulic jump at low flows
    lengthCarrWalls: L3, m, this is driven by the dimension of the Carr wall, a mechanism to prevent erosion from eddies forming on the outlet, 
    minimum requirement is 1.5 * water depth
    D: Water depth plus freeboard, m, freeboard is normally 0.2 * waterDepth + 0.3
    H: m, this one is a bit of a guess, but the idea is to absolutely guarantee the formation of a hydraulic jump, pipe diameter/2 will do that
    cutoffMain: This dimension is critical to percolation check for the structure; basically, vertical walls count, horizontals not so much, 
    this gives us 4m of vertical flow path, we need it
    cutoffT: Cutoff thickness, m, absolute minimum thickness for a deep trench, can be excavated by machine with narrow bucket
    t: Wall thickness, m, based roughly on D/12
    spacingBarrels: Spacing of barrels, m, this is a bit conservative, 2 drivers are adequate space to ensure compaction of pipe embedment material, 
    and to line up with the above ‘t’, and the necessary offsets for the gates; final design may reduce this, but I wanted to be conservative; 
    it may be down to 4.7 in final design
    noPiers: Number of intermediate piers, equals the number of barrels minus 1
    frontCutoff: m, this is another point that helps to percolation protection, but it really a guess
    rockworkBeyondCarrWalls: Extension of rockwork beyond Carr Walls, m
    rockworkUpstream: Extension of rockwork upstream, m
    concRate: Unit cost of concrete, $/m3
    '''

    # get individual and total costs
    res = {}
    res['Concrete'] = getCostConcrete_inlineControlStructure(lengthGateSection, lengthDissipater, lengthCarrWalls, D, H, cutoffMain, cutoffT, t, spacingBarrels, noPiers, frontCutoff, concRate)
    res['Steelwork'] = getCostSteelworkSansControlGates('INLINE', noPiers, t, spacingBarrels, D)
    res['Rockwork, d/s'] = getCostDownstreamRockwork(bedWidth, D, lengthCarrWalls, rockworkBeyondCarrWalls)
    res['Rockwork, u/s'] = getCostUpstreamRockwork(bedWidth, D, rockworkUpstream)
    res['Control gates'] = getCostControlGates(waterDepth, noPiers, num_barrels, pipeOD)
    res['Total'] = res['Concrete']['Total'] + res['Steelwork']['Total'] + res['Rockwork, d/s']['Total'] + res['Rockwork, u/s']['Total'] + res['Control gates']
    return res


def getTotalCost_siphonInlet(bedWidth = 50, waterDepth = 4.4, D = 5.4, \
    lengthMainSection = 10, t = 0.4, cutoffMain = 1, cutoffT = 0.2, noPiers = 7, spacingBarrels = 5, \
    pipeOD = 4, rockworkUpstream = 10, concRate = 1200, num_barrels = 1):
    '''
    bedWidth: Width of the channel bed, m
    waterDepth: d, water depth, m
    D: Water depth plus freeboard, m, freeboard is normally 0.2 * waterDepth + 0.3
    lengthMainSection: L1, m, normally twice of the gate length, to allow the gate to fully open, but flow to settle
    t: Wall thickness, m, based roughly on D/12
    cutoffMain: To achieve percolation resistance for the total structure
    cutoffT: Cutoff thickness, m, based on the practical minimum for machine excavation for 1 m deep
    noPiers: Number of intermediate piers, equals the number of barrels minus 1
    spacingBarrels: Spacing of barrels, m, this is conservative and could be as low as D + 0.6
    pipeOD: Pipe outer diameter, m, from pipe selection
    rockworkUpstream: Extension of rockwork upstream, m
    concRate: Unit cost of concrete, $/m3
    '''

    # define additional parameters
    lengthInletSection = 2 * D #L2, m, based on the earth batter at 1:2

    # get individual and total costs
    res = {}
    res['Concrete'] = getCostConcrete_siphonInlet(lengthMainSection, lengthInletSection, D, t, cutoffMain, noPiers, spacingBarrels, cutoffT, pipeOD, concRate)
    res['Steelwork'] = getCostSteelworkSansControlGates('SIPHON-INLET', noPiers, t, spacingBarrels, D)
    res['Rockwork, u/s'] = getCostUpstreamRockwork(bedWidth, D, rockworkUpstream)
    res['Control gates'] = getCostControlGates(waterDepth, noPiers, num_barrels, pipeOD)
    res['Total'] = res['Concrete']['Total'] + res['Steelwork']['Total'] + res['Rockwork, u/s']['Total'] + res['Control gates']
    return res


def getTotalCost_siphonOutlet(bedWidth = 50, D = 5.4, H = 1, lengthCarrWalls = 7, \
    t = 0.4, cutoffMain = 2, cutoffT = 0.3, noPiers = 7, spacingBarrels = 5, \
    pipeOD = 4, rockworkBeyondCarrWalls = 15, concRate = 1200):
    '''
    bedWidth: Width of the channel bed, m
    D: Water depth plus freeboard, m, freeboard is normally 0.2 * waterDepth + 0.3
    H: m, based on the minimum to stop shooting flow along the base
    lengthCarrWalls: L3, m, based on 1.5 times the downstream water depth
    t: Wall thickness, m, based roughly on D/12
    cutoffMain: Co, m, based on the percolation path from the downstream channel to the pipe bedding
    cutoffT: Cutoff thickness, m, based on the minimum for construction and hence, a function of cutoffMain
    noPiers: Number of intermediate piers, equals the number of barrels minus 1
    spacingBarrels: Spacing of barrels, m, this is conservative and could be as low as 4.6 m without introducing too much construction difficulty
    pipeOD: Pipe outer diameter, m, from pipe selection
    rockworkBeyondCarrWalls: Extension of rockwork beyond Carr Walls, m
    concRate: Unit cost of concrete, $/m3
    '''

    # define additional parameters
    lengthFirstSection = D #L1, m, length of intermediate piers, needs to be long enough to allow bulkheads to be fitted and to allow for flow training; suggest D as a minimum
    lengthMidSection = 2 * spacingBarrels - lengthFirstSection #L2, m, based on L1 + L2 being 2 * spacingBarrels

    # get individual and total costs
    res = {}
    res['Concrete'] = getCostConcrete_siphonOutlet(lengthFirstSection, lengthMidSection, lengthCarrWalls, D, H, t, cutoffMain, cutoffT, noPiers, spacingBarrels, pipeOD, concRate)
    res['Steelwork'] = getCostSteelworkSansControlGates('SIPHON-OUTLET', noPiers, t, spacingBarrels, D)
    res['Rockwork, d/s'] = getCostDownstreamRockwork(bedWidth, D, lengthCarrWalls, rockworkBeyondCarrWalls)
    res['Total'] = res['Concrete']['Total'] + res['Steelwork']['Total'] + res['Rockwork, d/s']['Total']
    return res


def getTotalCost_siphonBaffledOutlet(bedWidth = 50, lengthCarrWalls = 7, t = 0.5, pipeOD = 4, \
    spacingBarrels = 10.5, Ddown = 5.4, noPiers = 1, rockworkBeyondCarrWalls = 15, concRate = 1200):
    '''
    bedWidth: Width of the channel bed, m
    lengthCarrWalls: L3, m, based on 1.5 times the downstream water depth
    H: m, called f in USBR, is W/6
    t: Wall thickness, m, based roughly on D/12 and structure projecting 2 m above n/s
    pipeOD: Pipe outer diameter, m, from pipe selection
    spacingBarrels: W, width of the structure per barrel, m, based on Froude Number calculations (and not on the pipe spacing)
    Ddown: Downstream height, m, based on the water depth + freeboard, corresponding to D for the other control structures
    noPiers: Number of intermediate piers, equals the number of barrels minus 1
    rockworkBeyondCarrWalls: Extension of rockwork beyond Carr Walls, m
    concRate: Unit cost of concrete, $/m3
    '''

    from math import ceil

    # define additional parameters
    # NOTE: W = spacingBarrels, based on Froude Number calculations (and not on the pipe spacing)
    lengthFirstSection = 0.5 * spacingBarrels #L1, m, based on 0.5*W
    lengthMidSection = 4 / 3 * spacingBarrels - lengthFirstSection #L2, m, based on L1 + L2 = 4/3*W
    Dup = 3 / 4 * spacingBarrels #Depth of the structure at the upstream end, m, called H in USBR, is 3/4*W
    H = spacingBarrels / 6 #m, called f in USBR, is W/6
    lipWidth = ceil(spacingBarrels / 12 * 10) / 10 #m, based on e in USBR, is W/12 (round up to 1 decimal)

    # get individual and total costs
    res = {}
    res['Concrete'] = getCostConcrete_siphonBaffledOutlet(lengthFirstSection, lengthMidSection, lengthCarrWalls, Dup, H, t, pipeOD, spacingBarrels, Ddown, noPiers, lipWidth, concRate)
    res['Steelwork'] = getCostSteelworkSansControlGates('SIPHON-OUTLET', noPiers, t, spacingBarrels, Ddown)
    res['Rockwork, d/s'] = getCostDownstreamRockwork(bedWidth, Ddown, lengthCarrWalls, rockworkBeyondCarrWalls)
    res['Total'] = res['Concrete']['Total'] + res['Steelwork']['Total'] + res['Rockwork, d/s']['Total']
    return res
