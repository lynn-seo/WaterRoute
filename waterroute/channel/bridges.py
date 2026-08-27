# seline.ng@csiro.au, start: 8/dec/2020, clayton vic
import math

def getTotalCostPerCrossing(crossingType = 'farm access', BW = 50,
    D = 5.4, concRate = 1200, deckRate = 800):

    '''
    Function purpose:
    -------------------
    To compute the cost of a farm access, highway or railway crossing

    Inputs:
    ---------
    crossingType, string = type of crossing, 'farm access', 'highway' or 'railway'
    BW, float = parameter, >= 0
    D, float = parameter, >= 0
    concRate, float = concrete rate for structures in $/m3, >= 0
    deckRate, float = rate for decks in $/m2, >= 0
   
    Output:
    ---------
    costFinal, dictionary:
    - costFinal['Concrete'], float = concrete cost, $
    - costFinal['Deck'], float = deck cost, $
    - costFinal['Fencing'], float = fencing or guardrails cost, $
    - costFinal['Total'], float = total cost of crossing, $

    Corresponding Excel sheet:
    -----------------------------
    Bridges
    '''

    # ==============
    # check inputs
    # ==============

    if not(crossingType.lower() == 'farm access' or crossingType.lower() == 'highway' or crossingType.lower() == 'railway'):
        raise ValueError('crossingType must be one of the following: "farm access", "highway", "railway".')

    if BW < 0 or D < 0 or concRate < 0 or deckRate < 0:
        raise ValueError('One of the following is negative when it should be positive: BW, D, concRate, deckRate.')

    # ============================================
    # assign parameters according to bridgeType
    # ============================================

    if crossingType.lower() == 'farm access':
        # Will span from top of Batter to Batter
        # Will be pre-stressed, post tensioned deck units, with spans < 25m.
        # Will cost on a m2 basis
        # Cost Intermediate piers as normal reinforced conc
        # Assume 6 m to inside of kerb units
        # Assume end headstocks are same as intermediate headstocks
        # Cost kerb units as decks (add 1.2 m to width)
        widthBetweenKerbs = 6; W = 3; DF = 0.4; t = 0.5; H = 5.1; HS = 0.5; HW = 1; costFencing = 20000

    elif crossingType.lower() == 'highway':
        # Width assumptions: (i) lane and verge, 7.6 m, (ii) minimum width, 19.2 m, (iii) use, 20 m
        # Increase sections compared to farm access, as dynamic loads will apply
        widthBetweenKerbs = 20; W = 4; DF = 0.6; t = 0.75; H = 5.1; HS = 0.6; HW = 1.5; costFencing = 100000

    elif crossingType.lower() == 'railway':
        # Width assumptions: (i) track and ballast and access, 10 m, (ii) use, 10 m
        widthBetweenKerbs = 10; W = 4; DF = 0.6; t = 0.75; H = 5.1; HS = 0.6; HW = 1.5; costFencing = 20000

    else:
        raise valueError('Unrecognized crossingType.')
        
    # ======================
    # compute cost of deck
    # ======================

    batterToBatter = BW + 2 * D * 2 #default 68.8 for farm crossing, therefore 2 intermediate piers required
    deckArea = batterToBatter * (widthBetweenKerbs + 1.2) #m2
    costDeck = deckArea * deckRate #$
    
    # ============================================================
    # compute concrete volumes for piers and total concrete cost
    # ============================================================

    widthOfPiers = widthBetweenKerbs + (2 * 0.3) + (2 * 0.25)

    concVolPerPier = {}
    concVolPerPier['Base'] = W * DF * widthOfPiers
    concVolPerPier['Riser'] = widthOfPiers * (H - HS) * t
    concVolPerPier['Headstock'] = HS * HW * widthOfPiers
    concVolPerPier['Total'] = concVolPerPier['Base'] + concVolPerPier['Riser'] + concVolPerPier['Headstock']

    noPiersPerSite = math.floor(batterToBatter/25) #max 25m span
    concVolPiers = concVolPerPier['Total'] * noPiersPerSite #m3
    costConcVolPiers = concVolPiers * concRate #$



    # =========================================================================
    # compute cost of fencing
    # farm access crossing: fencing at and across bridge (steel rails)
    # highway crossing: guardrails across the bridge, about 90 m per side
    # railway crossing: fencing at and across bridge (steel rails)
    #==========================================================================

    costFencing = costFencing #$

    # ==============================
    # compute total cost
    # ==============================

    costFinal = {}
    costFinal['Concrete'] = costConcVolPiers #concrete cost, $
    costFinal['Deck'] = costDeck #deck cost, $
    costFinal['Fencing'] = costFencing #fencing or guardrails cost, $
    costFinal['Total'] = costConcVolPiers + costDeck + costFencing #total cost of crossing, $
    return costFinal


###############################
# TESTING
###############################

if __name__ == "__main__":

    # From HGD to Buchannan
    print('From HGD to Buchannan')
    BW = 50; D = 5.4
    print('Cost of a farm access bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='farm access')))
    print('Cost of a highway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='highway')))
    print('Cost of a railway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='railway')))
    print('')

    # Part 1-From Buchannan to Big Drop
    print('Part 1-From Buchannan to Big Drop')
    BW = 45; D = 4.8
    print('Cost of a farm access bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='farm access')))
    print('Cost of a highway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='highway')))
    print('Cost of a railway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='railway')))
    print('')

    # Part 3-From Big Drop to Cornish
    print('Part 3-From Big Drop to Cornish')
    BW = 55; D = 3.45
    print('Cost of a farm access bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='farm access')))
    print('Cost of a highway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='highway')))
    print('Cost of a railway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='railway')))
    print('')

    # Sect 4-Cornish Ck to Aramac Ck
    print('Sect 4-Cornish Ck to Aramac Ck')
    BW = 40; D = 4.25
    print('Cost of a farm access bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='farm access')))
    print('Cost of a highway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='highway')))
    print('Cost of a railway bridge in $: ' + str(getTotalCostPerCrossing(BW = BW, D = D, crossingType='railway')))
    print('')




