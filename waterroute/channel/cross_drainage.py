# seline.ng@csiro.au, start: 7/dec/2020, clayton vic

import math

# key = diameter, value = supplyAndInstallRate
_supplyAndInstallRate = {1500: 1280, 1800: 1810, 2100: 2400}

def select_cross_drain(cross_drain_discharge, cross_drain_length):
    '''
    return tuple of cross_drain_cost,number_of_barrels,diameter 
    '''
    r = []
    for diameter in _supplyAndInstallRate.keys():
        number_of_barrels = math.ceil(cross_drain_discharge / getFlowCapacity(diameter, 1))
        r.append((getCrossDrainageUnitCost(diameter, cross_drain_length, number_of_barrels,),number_of_barrels,diameter))
    cross_drain_cost,number_of_barrels,diameter = min(r)
    return (cross_drain_cost,number_of_barrels,diameter)

def getCrossDrainageUnitCost(diameter, length, noBarrels, concRate = 1200):

    '''
    Function purpose:
    -------------------
    To compute the unit cost of cross drainage

    Inputs:
    -----------------------
    diameter, int: 1500, 1800 or 2100
    length, float: >= 0
    nbarrels, float: number of barrels, >= 0
    concRate, float: concrete rate for structures in $/m3, >= 0

    Output:
    -----------------------
    cost, float: cost in $

    Assumptions:
    -----------------------
    1. We will have a tabulation of locations, design capacity, elevations wrt n/s, and cross drain length
    2. We can assume there is sufficient cross fall for all cross drains to be inlet controlled
    3. We will standardise the diameter, and then use number of barrels to cover the varying capacities
    4. We will use 4 m/sec as this will be the limit without baffled outlets

    Corresponding Excel sheet:
    -----------------------------
    cross drainage_unit costs
    '''

    # ==============
    # check inputs
    # ==============

    if diameter not in _supplyAndInstallRate.keys():
        raise ValueError(f'Value of diameter must be one of {list(_supplyAndInstallRate.keys())}.')

    if length < 0 and noBarrels < 0 and concRate < 0:
        raise ValueError('One of the following is negative when it should be positive: length, noBarrels and concRate.')

    # ==================
    # compute cost ($)
    # ==================

    # supply and install rate
    # source: Sunwater pipe supply and install sheet

    # concrete requirements on a per barrel basis
    # weep collars, inlet structure, outlet structure
    concReqWeepCollars = ((1.2 + diameter/1000) * (1.2 + diameter/1000) * 0.25 - 3.14 * ((diameter/1000)**2 / 4 * 0.25))*2
    concReqInletStruct = (2.5 + 2.5 + 0.6) * 0.2 * ((diameter/1000) + 1.2)
    concReqOutletStruct = (2.5 + 5.5 + 0.6) * 0.2 * ((diameter/1000) + 1.2)
    concReqTotal = concReqWeepCollars + concReqInletStruct + concReqOutletStruct

    # cost ($)
    cost = (_supplyAndInstallRate[diameter] * length * noBarrels) + (noBarrels * concReqTotal * concRate)

    return cost



def getFlowCapacity(diameter, noBarrels):

    '''
    Function purpose:
    -------------------
    To compute flow capacity

    Inputs:
    ---------
    diameter, int: 1500, 1800 or 2100
    nbarrels, float: number of barrels, >= 0

    Output:
    ---------
    flowCapacity, float: flow capacity in cumec

    Corresponding Excel sheet:
    -----------------------------
    cross drainage_unit costs
    '''

    if diameter not in [1500, 1800, 2100]:
        raise ValueError('Value of diameter must be 1500, 1800 or 2100.')

    if noBarrels < 0:
        raise ValueError('noBarrels cannot be negative.')

    flowAt4MPerSec = 3.14 * (diameter/1000)**2 / 4 * 4
    flowCapacity = noBarrels * flowAt4MPerSec #cumec

    return flowCapacity

