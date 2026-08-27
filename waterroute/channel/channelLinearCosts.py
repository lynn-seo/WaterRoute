# seline.ng@csiro.au, start: 24/dec/2020, clayton vic

def getTotalChannelLinearCost(reserveWidth = 120, noGates = 1):
    '''
    Function to compute the total channel linear cost

    Inputs:
    --------
    reserveWidth, float = reserve width, m
    noGates, int = number of gates per km of channel

    Output:
    --------
     res, dictionary
    - res['Easement'], float = total cost of easement excluding catch drainage, $/km
    - res['Catch drainage'], float = total cost of catch drainage, $/km
    - res['Total'], float = total linear cost of channel, $/km

    Corresponding Excel sheet:
    -----------------------------
    channel Linear costs
    '''

    import pandas as pd
    import numpy as np

    # easement costs per km of channel
    # -------------------------------------
    # this includes:
    # -clearing and grubbing
    # -stock fencing of easement
    # -provision of temporary arrangement and reinstatement of existing facilities
    # -assume cleared 10 m outside easement for fence access

    # assume
    clearedArea = (reserveWidth + 20) * 1000 / 10000 #ha/km
    fencingRequirements = 2 #km per km of channel
    gates = noGates #number per km of channel
    topsoil = 1000 * reserveWidth * 0.15 #m3 per km of channel
    replaceTopsoil = 15 #m of slope avg

    # lists of items and unit rates
    # warning: do not change the order of the elements in items
    items = ['Clear, grub, burn', 'Fencing', 'Gates, strainers', 'Reinstate', 'Strip topsoil', 'Replace topsoil', 'Access road on berm']
    rate = [1300, 2000, 800, 4000, 4, 1.5, 10]

    # compute quantity of each individual item
    qty = np.empty(len(items)) * np.nan
    qty[0] = clearedArea #clear, grub, burn [ha/km]
    qty[1] = fencingRequirements #fencing [km/km]
    qty[2] = gates #gates, strainers [number/km]
    qty[3] = 1 #reinstate

    # CHECK: typo in excel sheet, G21, 1800 or 18000?
    qty[4] = topsoil #strip topsoil, m3/km

    qty[5] = 1000 * replaceTopsoil #replace topsoil
    qty[6] = 5 * 1000 #access road on berm

    # create pandas dataframe
    di = {'Items': items, 'Qty': qty, 'Rate': rate, 'Cost': np.zeros(len(items))}
    df = pd.DataFrame(di)

    # compute individual and total costs
    df['Cost'] = df['Qty'] * df['Rate']
    res = {}; res['Easement'] = np.sum(df['Cost'])

    # catch drainage allowance per km of channel
    # ----------------------------------------------
    # assumptions:
    # -since the banks are fairly high, formed catch drains will only be necessary for a portion of the route
    # -The primary purpose will be to divert small adjacent catchments to the ones featuring cross drains
    # -use 4 m2 per m for 50 % of length, gives $1/m

    effectiveCatchDrainageCost = 12000 #per km
    res['Catch drainage'] = effectiveCatchDrainageCost

    # total linear cost
    # ----------------------------------------------
    res['Total'] = res['Easement'] + res['Catch drainage']
    return res



###############################
# TESTING
###############################

if __name__ == "__main__":

    print(getTotalChannelLinearCost(reserveWidth = 120, noGates = 1))