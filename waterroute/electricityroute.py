import pandas as pd
from dataclasses import dataclass, field, asdict
from operator import attrgetter

from waterroute import economics

#powerline price model

#more then you ever want to know about netowrk standards
#https://www.ausgrid.com.au/-/media/Documents/Technical-Documentation/NS/NS220.pdf

#inflation factors
#https://www.abs.gov.au/statistics/economy/price-indexes-and-inflation/producer-price-indexes-australia/jun-2020#data-downloads-data-cubes
inflation_2015_2020 = 117.5/107.6 #from  heavy construction index
inflation_2019_2020 = 117.5/115.8 #from  heavy construction index
#inflation_2015_2020 = 1.08498 #jennys number

power_factor = 0.8

#didnt included 1000kV as this seems to be the cutting edge of AC power lines and none in use in Aus
#i dont know if the rest of the model holds for these though.
standard_voltages_kV = [11,22,33,66,132,275,330,500]#,1000]

#raw data for the line limit model
_cols = ['km','kV','MW']
_data = [[50,66,22],
[50,132,120],
[50,220,400],
[50,400,1500],
[100,66,11],
[100,132,60],
[100,220,200],
[100,400,750],
[200,66,5.5],
[200,132,30],
[200,220,100],
[200,400,375],
[300,66,3.6],
[300,132,20],
[300,220,67],
[300,400,250],
]
_df = pd.DataFrame(_data,columns=_cols)

#def line_MW_limit_v1(length_km:float,line_V_kV:float):
#    '''MW load limits on lines by transmission distance
#    from table 2 in http://www.iitk.ac.in/npsc/Papers/NPSC1998/p2.pdf'''
#    assert line_V_kV > 50, 'this model doesn not apply for voltages this low'
#    limit = (0.5542*line_V_kV**2 - 37.217*line_V_kV + 1205.2)/length_km
#    return limit

#def _prep_line_kV_limit ():
#    '''rearranging the formula'''
#    MW, line_V_kV, length_km = sympy.symbols(['MW', 'line_V_kV', 'length_km'],real=True)
#    eq = sympy.Eq(MW,(0.5542*line_V_kV**2 - 37.217*line_V_kV + 1205.2)/length_km)
#    eq2 = sympy.solve(eq, line_V_kV)[1]
#    fnc = sympy.lambdify([length_km,MW],eq2)
#    return fnc

#_line_kV_limit_func = _prep_line_kV_limit ()

#def line_kV_limit (length_km:float,MW:float):
#    '''returns the minimum voltage (kV) transmission line that should be used
#    from this you should round up to the next standard size
#    
#    the model isnt applicable to low voltage lines so if the result is lower then 50kV then 50 is returned'''
#    kv_limit = _line_kV_limit_func(length_km,MW)
#    if kv_limit < 50:
#        return 50
#    if np.isnan(kv_limit):
#        return 50
#    return kv_limit

def line_kv_lower_limit(MW:float):
    '''notes on current limits
    indian doc that specifies the standard conductor for each line kv rating and the current carrying capacity
    https://upenergy.in/upptcl/en/article/current-carrying
    current capacities of many different lines
    https://www.prysmiancable.com.au/wp-content/uploads/2018/05/PGAU_BOHC2018-05_20pg_lowres.pdf
    lines used for varous voltages in QLD - unfortunatly not very perscriptive
    https://www.ergon.com.au/__data/assets/pdf_file/0010/146584/STNW3355_V4-Subtransmission-Line-Design-Manual.pdf
    it looks like the indian current limits are broadly in line with the QLd values.
    Lets use the indian data as broadly indicative of the current limits
    These are then converted to 3 phase MW limits using 0.8 power factor
    The result is the following relationship for a new conductor at 75C limit'''
    #MW = 0.8331*line_V_kV - 38.191
    line_V_kV = 1.1963*MW + 46.374
    return line_V_kV


#def line_MW_limit(length_km:float,line_V_kV:float):
#    '''MW load limits on lines by transmission distance
#    from table 2 in http://www.iitk.ac.in/npsc/Papers/NPSC1998/p2.pdf'''
#    MW = (0.0684344*line_V_kV**2.3260363)/length_km
#    return MW


def line_kV_limit (length_km:float,MW:float,allowed_undersize:float = 0.0):
    '''returns the minimum voltage (kV) transmission line that should be used
    from this it is round up to the next standard size'''
    #print (MW,length_km)
    #formula rearranged with https://www.wolframalpha.com/widgets/view.jsp?id=b833b9ec052dba3d28359e9de6532441
    #i couldnt get sympy to do the rearranging
    min_line_V_kV = 3.16763*(length_km * MW)**0.4299159045798210
    #The above relationship ends up with current s that are too high for very short lines so lets put a floor in the kv
    min_line_V_kV = max([min_line_V_kV,line_kv_lower_limit(MW)])
    #if the model ends up over the max in standard_voltages_kV need to use the max.
    line_V_kV = min([max(standard_voltages_kV),*[v for v in standard_voltages_kV if v>(1-allowed_undersize)*min_line_V_kV]])
    return line_V_kV

def transformer_cost(kV1:int,kV2:int):
    ''' from analysis of Jenny's spreadsheet
     APGT (2015) - I think this is the report https://anlecrd.com.au/projects/australian-power-technology-generation-technology-apgt-study/
     return is in dollars AUD 2020'''
    direct_cost = 61305*abs(kV1-kV2)
    return direct_cost * inflation_2015_2020

def transmission_line_cost_old (length_km:float,kV:int,double_circuit:bool=False):
    ''' from analysis of Jenny's spreadsheet
     APGT (2015) - I think this is the report https://anlecrd.com.au/projects/australian-power-technology-generation-technology-apgt-study/
     return is in dollars AUD 2020'''
    direct_cost_per_km = 3297.8 * kV + 80652

    #extra circuts are on the same towers. and a significant cost saving
    circuit_factor = 1
    if double_circuit:
        circuit_factor = 1.3
    #no economies of scale assumed for long lengths
    return direct_cost_per_km * length_km * circuit_factor * inflation_2015_2020

def transmission_line_cost (length_km:float,kV:int,double_circuit:bool=False):
    ''' from analysis of AEMO assumptions and inputs workbook 2019
     return is in dollars AUD 2020
     1. For transmission lines longer than 150 km, an economy of scale factor of 0.95 is applied to $/km cost estimate.
     '''
    #derived from https://aemo.com.au/-/media/files/electricity/nem/planning_and_forecasting/inputs-assumptions-methodologies/2020/2019-input-and-assumptions-workbook-v1-5-jul-20.xlsx?la=en
    if double_circuit:
        #extra circuts are on the same towers. and a significant cost saving
        direct_cost_per_km = 2886.92 * kV + 1008387
    else:
        direct_cost_per_km = 2360.6 * kV + 848463

    #economies of scale assumed for long lengths
    economies_scale = 1
    if length_km >= 150:
        economies_scale = 0.95
        
    return direct_cost_per_km * length_km * inflation_2019_2020 * economies_scale

def annual_tuos (peak_usage_kW,annual_consumption_kWh,locational_price=4):
    '''assumes that peak usage occurs in every month
    for simplicity you can probably equate peak usage to installed capacity
    based on
    https://www.powerlink.com.au/sites/default/files/2020-03/2020-2021%20Schedule%20of%20Transmission%20Shared%20Network%20Prices.pdf
    Prescribed TUOS Service
    Non-Locational
    Demand Based Price
    ($/kW/month) = 1.6148 
    
    Prescribed TUOS Service
    Non-Locational
    Energy Based Price
    (c/kWh) = 0.4719
    
    Prescribed Common
    Transmission Service
    Demand Based Price
    ($/kW/month) = 1.8667
    
    Prescribed Common
    Transmission Service
    Energy Based Price
    (c/kWh) = 0.5456
    
    Prescribed TUOS Service
    Locational Price
    ($/kW/month)
    This varies quite widely from under $ 0.3 to over $8
    since we are generally considering remote locations $4 seems reasonable.
    '''
    locational = peak_usage_kW * locational_price * 12
    nl_demand = peak_usage_kW * 1.6148 * 12
    nl_energy = annual_consumption_kWh * (0.4719/100)
    nl = min(nl_demand,nl_energy)
    ts_demand = peak_usage_kW * 1.8667 * 12
    ts_energy = annual_consumption_kWh * (0.5456/100)
    ts = min(ts_demand,ts_energy)
    #print (locational,nl_demand,nl_energy,ts_demand,ts_energy)
    return sum ([locational,nl,ts])


@dataclass
class GridCost:
    overhead_rate: float = 0.3
    contingency_rate: float = 0.4
    length_km: float = None
    MW: float = None
    connection_V_kV: int = None
    double_circuit:bool = None
    line_V_kV: int = None
    transformer_cost: float = None
    transmission_line_cost: float = None
    direct_cost: float = None
    total_cost: float = None
    overheads: float = None
    contingencies: float = None
        
    def asdict(self,prefix=''):
        return {f'{prefix}{key}': val for key, val in asdict(self).items()} 
        
        
def grid_cost(length_km:float,MW:float,connection_V_kV:int = 11,line_V_kV:int = None,include_transformer:bool=True):
    #im going to assum that all generators and pumps are operating at 11kV
    #The data that we found indicates some of the bigger ones may operate at 
    #maybe double that but in large schemes like that the cost different from this assumption is small
    def line_cost(double_circuit:bool=False):
        g = GridCost()
        g.length_km = length_km
        g.MW = MW
        g.double_circuit = double_circuit
        g.connection_V_kV = connection_V_kV
        d = 1
        if double_circuit:
            d = 2
        g.line_V_kV = line_V_kV
        if line_V_kV is None:
            g.line_V_kV = line_kV_limit(g.length_km,g.MW/d)
        g.transformer_cost = 0
        if include_transformer:
            g.transformer_cost = transformer_cost(g.connection_V_kV,g.line_V_kV)
        g.transmission_line_cost = transmission_line_cost(g.length_km,g.line_V_kV,g.double_circuit)
        if MW == 0:
            g.transformer_cost = 0
            g.transmission_line_cost = 0
        
        g.direct_cost = g.transformer_cost + g.transmission_line_cost
        g.total_cost,g.overheads,_insurance,g.contingencies = economics.indirect_costs(g.direct_cost,g.overhead_rate,0,g.contingency_rate)
        return g
    results = [line_cost(False),line_cost(True)]
    return min(results, key=attrgetter('total_cost'))
