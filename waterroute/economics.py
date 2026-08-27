import numpy as np
import numpy_financial as npf
from math import inf

def indirect_costs(direct_cost:float,overhead_rate=0.3,insurance_rate=0.0,contingency_rate = 0.4):
    #overheads = overhead_rate * direct_cost
    #insurance = insurance_rate * (direct_cost + overheads)
    #contingencies = contingency_rate * (direct_cost + overheads + insurance)
    #total_cost = direct_cost + overheads + insurance + contingencies
    #return (total_cost,overheads,insurance,contingencies)
    return inflate_cost(direct_cost,overhead_rate,insurance_rate,contingency_rate)

def inflate_cost(cost,*ratios):
    result = []
    for r in ratios:
        indirect_cost = cost * r
        cost+=indirect_cost
        result.append(indirect_cost)
    return (cost,*result)

def lcoe(annual_output, capital_cost, annual_operating_cost, discount_rate, lifetime,capital_cost_years = 1):
    """Compute levelised cost of electricity
    Arguments
    ---------
    annual_output : float
    capital_cost : float
    annual_operating_cost : float
    discount_rate : float
    lifetime : int
    capital_cost_years: int default of 1, min of 1. spreads the capital cost over the first x years. Delays annual output by the same amount of time. lifetime starts after x.
    Returns
    -------
    float
    adapted from https://github.com/willu47/lcoe/blob/master/src/lcoe/lcoe.py
    MIT License
    """
    #annual_cost_capital = npf.pmt(discount_rate, lifetime, -capital_cost)
    #total_annual_cost = annual_cost_capital + annual_operating_cost
    #lcoe = total_annual_cost / annual_output if annual_output>0 else inf
        
    annual_expenses = [*np.repeat(capital_cost/capital_cost_years,capital_cost_years),*np.repeat(annual_operating_cost,lifetime)]
    npv_expenses = npf.npv(discount_rate,annual_expenses)
    annual_generation = [*np.repeat(0,capital_cost_years),*np.repeat(annual_output,lifetime)]
    npv_generation = npf.npv(discount_rate,annual_generation)
    lcoe = npv_expenses/npv_generation if npv_generation>0 else inf
    
    return lcoe