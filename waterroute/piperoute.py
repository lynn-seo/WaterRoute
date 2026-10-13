from operator import attrgetter
from functools import total_ordering

import itertools
from math import pi,log10,inf
import math
import copy
from dataclasses import dataclass, field, asdict
from scipy import optimize
from scipy.optimize import root,minimize_scalar,minimize,differential_evolution,brute
import fluids
import numpy_financial as npf
import numpy as np

import pandas as pd
import seaborn as sns


seconds_per_day = 60*60*24
gravity = 9.80665 #standard gravity is, by definition, 9.80665 m/s2

from waterroute import economics
from waterroute import electricityroute

def clover_GRP_IDs():
    '''scrape all the standard GRP pipe diameters from clover
    
    however lets not limit ourselves to these dimensions other sizes may be available
    eg http://grandpipe.com/en/product-information (have 0.2m increments up to 4m and say they can make any size)
    https://www.amiblu.com/pressure-pipes/
    and FLOWTITE
    also in the real world pipe diameters can be mixed along a barrel to give a desired friction loss
    This mixing probably doesnt end up much different in cost then using a cost for the right intermediate size.
    
    this is the pandas series that is returned
    300     0.333
    375     0.411
    450     0.489
    525     0.566
    600     0.644
    675     0.721
    750     0.797
    900     0.891
    1000    0.990
    1100    1.085
    1200    1.185
    1300    1.284
    1400    1.382
    1500    1.482
    1600    1.581
    1700    1.678
    1800    1.777
    1900    1.875
    2000    1.973
    3000    2.956
    4000    3.941
    Name: PN1, PN6 & PN10, dtype: float64
    https://www.cloverpipe.com.au/product/grp-pipes/grp-pipe-and-fittings/#size-chart
    '''
    pipe_IDs = pd.read_html('https://www.cloverpipe.com.au/product/grp-pipe-and-fittings/',header=[0,1],index_col=0)[0]['ID']['PN1, PN6 & PN10']/1000
    return pipe_IDs

# Hardcoded copy of clover_GRP_IDs() so importing doesn't hit the network (checked 13/10/2026)
standard_GRP_Ds = pd.Series(
    [0.333, 0.411, 0.489, 0.566, 0.644, 0.721, 0.797, 0.891, 0.990, 1.085, 1.185,
     1.284, 1.382, 1.482, 1.581, 1.678, 1.777, 1.875, 1.973, 2.956, 3.941],
    index=[300, 375, 450, 525, 600, 675, 750, 900, 1000, 1100, 1200,
           1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 3000, 4000],
    name='PN1, PN6 & PN10',
)

def velocity_head(V):
    #(V**2)/(gravity*2)
    return fluids.head_from_K(1, V)



def head_loss(V,D,L,k_mm = 0.029 ):
    #k = Nikuradse's equivalent sand grain diameter
    #http://www.wrc.org.za/wp-content/uploads/mdocs/579-1-96.pdf
    #TBM k = 5.53mm
    #https://www.sciencedirect.com/science/article/pii/S0886779809000716?via%3Dihub#fig13
    #TBM k = designed at 10mm   measured average 7mm
    #a fairly clear summary of this calculation https://me.queensu.ca/People/Sellens/LossesinPipes.html
    # a full list of roughnesses can be found in these 2 spots
    #fluids.friction.roughness_clean_dict
    #fluids.friction.HHR_roughness
    #k can be obtained from
    #fluids.material_roughness("Fiberglass lining")
    #or getting fancy  fluids.roughness_Farshad('Fiberglass lining', D=D)
    #https://www.cloverpipe.com.au/wp-content/uploads/2016/03/grp-technical-design-manual.pdf
    #for concrete including some fouling k is 0.3mm pers comm kev
    #Equivalent roughness k = 0.02mm and 0.06mm ± 0.015mm from https://www.iplex.com.au/assets/Uploads/5188d1e29b/FLOWTITE-Pipeline-Design-Hydraulic.pdf
    k = k_mm/1000 #Colebrook-White k = 0.029 mm from https://www.cloverpipe.com.au/wp-content/uploads/2019/02/Superlit_GRP_Pipe_and_Fittings_v5-1.pdf
    rho = 998       #density of water @ 20C
    #from https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html?vA=20&units=C#
    #mu is dynamic viscosity Pa*s @ 20C
    mu = 0.0010005
    #nu is kinematic viscosity m^2/s @ 20C
    nu = fluids.nu_mu_converter(rho=rho, mu=mu)
    Re = fluids.Reynolds(V=V, D=D, rho=rho, mu=mu,nu=nu)
    fd = fluids.friction_factor(Re, eD=k/D)
    K = fluids.K_from_f(fd, L, D)
    #pressure drop [Pa]
    #print (fluids.dP_from_K(K, rho, V))
    #head loss
    return fluids.head_from_K(K, V)

def residual_head(total_head,V,D,L):
    loss = head_loss(V,D,L)
    return (total_head - loss)

def velocity_for_diameter(total_head,D,L):
    '''finds the velocity that would fully scrub the available head at the given diameter
    
    returns velocity in m/s'''
    f =  lambda V: residual_head(total_head,V,D,L)
    r = (root(f,1))
    if not r.success:
        raise RuntimeError (r.message)
    return r.x[0]
    
def diameter_for_velocity(total_head,V,L):
    '''finds the diameter that would fully scrub the available head at the given velocity
    
    returns diameter in metres'''
    f =  lambda D: residual_head(total_head,V,D,L)
    r = (root(f,1))
    if not r.success:
        raise RuntimeError (r.message)
    return r.x[0]

def diameter_for_Q(total_head,Q,L):
    '''finds the diameter that would fully scrub the available head at the given discharge
    
    returns diameter in metres'''
    def fcn(D):
        V = volume_to_velocity(Q,D=D)
        return residual_head(total_head,V,D,L)

    r = (root(fcn,1))
    if not r.success:
        raise RuntimeError (r.message)
    return r.x[0]

def cost_of_pipes(V,D,L,diversion_m3_per_s):
    '''not really useful as it returns fractions of a barrel'''
    vol_m3_per_s = velocity_to_volume(V,D=D)#*1000* seconds_per_day
    #num_barrels = math.ceil(diversion_m3_per_s/vol_m3_per_s)
    num_barrels = diversion_m3_per_s/vol_m3_per_s
    cost_of_pipe = sum(pipe_diameter_cost(D,L,num_barrels))
    return (cost_of_pipe,num_barrels) 

def pipe_diameter_cost(D,L=1,num_barrels=1):
    '''cost of pipe supply and install
    
    D: pipe diameter in metres
    
    returns $AUD 2020 cost of supply and install of GRP pipes'''
    #todo increase supply cost for higher pressure rated pipes(Bradfield 2020)
    # supply = 378.97174*D**2 + 46.60519*D + 99.03634
    # install = 25.38122*D**2 + 198.15220*D - 12.69834
    #Updated supply and install 23/09/2025, relationships were developed from costs provided by RLB 2024, $/mm
    supply = 1.1533*(D*1000) + 98
    install = 342.29 * math.exp(0.0008*(D*1000))

    supply = supply * num_barrels * L
    install = install * num_barrels * L
    return (supply,install)

def pipe_o_and_m_cost(direct_pipe_cost):
    return 0.01 * direct_pipe_cost

def operating_cost(W):
    ''' operating cost per watt of installed capacity reported in AUD2020
    based on this
    https://www.ifc.org/wps/wcm/connect/topics_ext_content/ifc_external_corporate_site/sustainability-at-ifc/publications/hydroelectric_power_a_guide_for_developers_and_investors
    published 2015 
    QUOTE:
    Annual O&M costs range from 1.0–4.0 percent of
    investment costs. The IEA assumes 2.2 percent for large HPP
    and 2.2–3 percent for small HPP. However, these figures do not
    include major electro-mechanical equipment replacement, which
    would be required a couple of times during the HPP lifespan and
    would raise average O&M costs to US$45/kW/year for large
    HPP; US$52/kW/year for small HPP. '''
    USD_to_AUD_2015 = 1.3175 #https://www.poundsterlinglive.com/best-exchange-rates/best-us-dollar-to-australian-dollar-history-2015
    #not using construction index as it doesnt seem applicable to O&M
    AUD_2015_to_AUD_2020_cpi = 1.078 #https://www.rba.gov.au/calculator/annualDecimal.html
    annual_O_and_M_per_kW = 50 * USD_to_AUD_2015 * AUD_2015_to_AUD_2020_cpi
    return (W/1000)* annual_O_and_M_per_kW
    

def pump_station_cost(W):
    '''make this better'''
    return hydro_power_station_cost(W)

def hydro_power_station_cost(W):
    '''capital cost of a hydro power station
    
    does not include connection to the grid or penstock?
    
    W: is the installed capacity of the station in watts
    
    based on Entura hydro costs
    https://www.csiro.au/en/Research/Major-initiatives/Northern-Australia/Current-work/NAWRA/Mitchell-report
    
    Note: this function needs to be updated with more/improved sources of data.
    
    returns $AUD 2020? cost'''
    AUD2018_to_AUD2020 = 1.046263345 #heavy construction index

    if W == 0:
        #no power/pump station
        capex = 0
    elif W > 0:
        #generation model
        #todo:refine model with more data
        capex = 2.538638 * W + 889000
    else:
        #pumping model
        #todo: replace with a pump specific model
        capex = 2.538638 * -W + 889000
    return capex * AUD2018_to_AUD2020


def hydro_power(Q,head,efficiency=0.9,rho=998,gravity=9.80665):
    '''returns hydro power in watts for a given discharge Q(m3/s) and head (m)'''
    assert head >= 0
    return efficiency*rho*Q*gravity*head
    
def pump_power(Q,head,efficiency=0.9,rho=998,gravity=9.80665):
    '''returns pumping power in watts for a given discharge Q(m3/s) and head (m)'''
    assert head >= 0
    return (1/efficiency)*rho*Q*gravity*head

def velocity_to_volume(V,A=None,D=None):
    '''converts m/s into m^3/s given either a cross sectional area or a pipe diameter.
    '''
    assert A is None, 'A has been giving me issues by accidently using it instead of D - disabling this functionality'
    if (A is not None and D is not None) or V is None or (A is None and D is None):
        raise ValueError('Inputs must be V and one of A and D.')
    if D is not None:
        A = pi * (D/2)**2
    return V * A
    
def volume_to_velocity(Q,A=None,D=None):
    '''converts m^3/s into m/s given either a cross sectional area or a pipe diameter.
    '''
    assert A is None, 'A has been giving me issues by accidently using it instead of D - disabling this functionality'
    if (A is not None and D is not None) or Q is None or (A is None and D is None):
        raise ValueError('Inputs must be Q and one of A and D.')
    if D is not None:
        A = pi * (D/2)**2
    return Q/A

def velocity_to_D(V,Q):
    '''finds the pipe diameter that can take a given dischage at a given velocity
    '''
    A = Q/V
    D = math.sqrt(A/pi)*2
    return D

def tunnel_cost(D,L,hard_rock=0.5,soft_rock=0.5):
    '''
    D is tunnel diameter in meters
    L is length in meters
    returns direct cost in $AUD2020
    hard_rock and soft_rock are proportions of the tunnel modelled as each type. Note that the proportion isnt applied to length as the scalling of length relates to the tunnel as a whole.
    
    based on the water tunnel relationships in the following paper.
    Planning level tunnel cost estimation based on statistical analysis of historical data
    Jamal Rostami, Mahmoud Sepehrmanesh, Ehsan Alavi Gharahbagh, Navid Mojtabai
    '''
    if D <= 0: return 0
    assert(hard_rock>=0 and soft_rock>=0 and hard_rock + soft_rock == 1)
    L_km = L/1000
    USD2008_to_AUD2008 = 1.1765 #exchange rate
    AUD2008_to_AUD2020 = 1.307 #heavy construction index
    
    #table 5 Water Hard rock mechanized Cost (M$) = 10^(0.553 + 0.866 log of length (km) + 1.23 log of diameter (m))
    hard_cost = 10**(0.553 + 0.866 * log10(L_km) + 1.23 * log10(D))
    #table 5 Water Soft ground mechanized  Cost(M$) = 10^(1.07 + 0.725 log of length (km) + 1.02 log of diameter (m))
    soft_cost = 10**(1.07 + 0.725 * log10(L_km) + 1.02 * log10(D))
    #mix the 2 costs - im not sure what the original authors would think of this but it lets us take an average of the 2 models when geology is uncertain
    cost = hard_cost * hard_rock + soft_cost * soft_rock
    return cost * 1E6 * USD2008_to_AUD2008 * AUD2008_to_AUD2020

@dataclass
class Pipe:
    pipe_type: str = ''
    hard_rock_ratio: float = 0.5 #only used for tunnels
    material_k_mm: float = None
    overhead_rate: float = 0.3
    contingency_rate: float = 0.4
    start_elev: float = None
    end_elev: float = None
    initial_head: float = None
    V: float = field(default=None,metadata={'units':'m/s'})
    D: float = field(default=2) #just a random default
    L: float = field(default=None)
    min_V: float = field(default=0.3)#default is 0.3m/s to entrain any sediments?
    #http://www.cwc.gov.in/sites/default/files/16052019.pdf a lit review of max velocity in concrete lined tunnels
    #6m/s is a commun upper max for water with no sediment
    max_V: float = field(default=6.5)#default is 6.5m/s for GRP pipe
    min_D: float = field(default=0.3)#default is 0.3m smallest practical pipe
    max_D: float = field(default=3.941)# default is 3.941 nominal 4m Internal diameter

    direct_cost:float = None
    overheads:float = None
    contingencies :float = None
    total_cost:float = None

    num_barrels: int = field(default=1)
    residual_head: float = field(default=None)
        
    def asdict(self,prefix=''):
        return {f'{prefix}{key}': val for key, val in asdict(self).items()} 

@total_ordering
@dataclass
class Pipe_Route:
    description: str = ''
    optimise_on: str = 'marginal_LCOE'
    optimise_minimise: bool = True
    optimise_abs: bool = False    
    minimum_head_tolerance: float = 0
    highpoint_storage: bool = False
    min_num_barrels:int = 1
    max_num_barrels: int = field(default=5)  #all parts of the pipeline have the same barrel limit
    
    
    diversion_m3_per_s: float = field(default=None)
    annual_diversions_l: float = field(default=None)
    water_source_cost_per_ML: float = 0
    hydro_efficiency: float = 0.9
    hydro_overhead_rate: float = 0 #the model is already using total costs
    hydro_contingency_rate: float = 0 #the model is already using total costs
    pump_efficiency: float = 0.9
    pump_overhead_rate: float = 0 #the model is already using total costs
    pump_contingency_rate: float = 0 #the model is already using total costs
    
    behind_the_meter_hours_per_day: float = 0
    behind_the_meter_LCOE: float = 0
    use_average_day_power_prices: bool = True
    power_sell_price: float = 53.41 #$/MWh 53.41 is the QLD average in the 2019/20 financial year.
    power_buy_price: float = 53.41 #$/MWh
    annual_connection_tariff:float = 0
    use_tuos: bool = True
    pump_electricity_retailer_margin: float = 0.1
    tuos_locational_price: float = 4
    include_transformer: bool = True
    _MW_per_ML: float = 0
    
    hydro_grid_distance_km:float= 0
    pump_grid_distance_km:float= 0

    base_scenario = None
    #base_case_npv: float = field(default=0)
    #base_case_hydro_annual_output: float = field(default=0)
    
    diversion_days_per_year: float = 365.25
    flow_annual_hours: float = field(default=None) #hours per year
    flow_hours_per_day: float = field(default=None) #hours   
    
    pipes: list = None
    residual_head: float = 0
    pumping_head: float = 0
    generation_head: float = 0
    total_head: float = 0

    hydro_direct_cost: float = 0
    hydro_overheads: float = 0
    hydro_contingencies: float = 0
    hydro_total_cost: float = 0
    hydro_installed_capacity: float = 0 #watts
    hydro_annual_output: float = 0 #Wh per year
    hydro_annual_gross_income:float = 0
    hydro_grid=None

    pump_direct_cost: float = 0
    pump_overheads: float = 0 
    pump_contingencies: float = 0
    pump_total_cost: float = 0
    pump_installed_capacity: float = 0
    pump_annual_consumption: float = 0 #Wh per year
    pump_metered_annual_consumption: float = 0
    pump_behind_the_meter_annual_consumption: float = 0        
    pump_annual_tuos_cost: float = 0
    pump_annual_metered_cost: float = 0
    pump_behind_the_meter_annual_cost: float = 0
    pump_grid=None

    operating_cost: float  = 0
    water_purchase_cost: float  = 0
    other_costs: float  = 0
    other_operating_costs: float  = 0
    discount_rate: float  = 0.07
    lifetime: float  = 100

    #direct_cost:float = None
    #overheads:float = None
    #contingencies :float = None

    total_cost:float = None
    total_npv:float = None
    total_irr:float = None
    total_mirr:float = None
    total_LCOE:float = None
    annualised_cost_per_ML:float = None

    marginal_total_cost:float = None
    marginal_npv:float = None
    marginal_irr:float = None
    marginal_mirr:float = None
    marginal_LCOE:float = None

    #allow_pumping:bool = True
    
    def asdict(self,prefix=''):
        d = {f'{prefix}{key}': val for key, val in asdict(self).items() if key not in  ['pipes','base_scenario','pump_grid','hydro_grid']}
        for i,p in enumerate(self.pipes):
            d.update(p.asdict(f'{prefix}pipes_{i}_'))
        if self.pump_grid is not None:
            d.update(self.pump_grid.asdict(f'{prefix}pump_grid_'))
        if self.hydro_grid is not None:
            d.update(self.hydro_grid.asdict(f'{prefix}hydro_grid_'))
        return d
    
    def __float__(self):
        return self._comp_val()
    
    def __sub__(self, other):
        assert self.optimise_on == other.optimise_on
        assert self.optimise_minimise == other.optimise_minimise
        assert self.optimise_abs == other.optimise_abs
        return (self._comp_val() - other._comp_val())
    
    def __eq__(self, other):
        assert self.optimise_on == other.optimise_on
        assert self.optimise_minimise == other.optimise_minimise
        assert self.optimise_abs == other.optimise_abs
  
        return (self._comp_val() == other._comp_val())

    def __lt__(self, other):
        assert self.optimise_on == other.optimise_on
        assert self.optimise_minimise == other.optimise_minimise
        assert self.optimise_abs == other.optimise_abs
  
        return (self._comp_val() < other._comp_val())
    
    def _comp_val(self):
        if self.optimise_abs:
            opt_val = abs(getattr(self, self.optimise_on))
        else:
            opt_val = getattr(self, self.optimise_on)
        if self.optimise_minimise:
            return opt_val
        else:
            return -opt_val

def flow_time_calcs(route: Pipe_Route):
    route.flow_annual_hours = route.annual_diversions_l/1000/route.diversion_m3_per_s/seconds_per_day*24
    route.flow_hours_per_day = route.flow_annual_hours/route.diversion_days_per_year
    assert route.flow_hours_per_day <= 24.01 #added a little margin for error

def route_calcs(route: Pipe_Route):
    '''before this need to set each pipes D and num_barrels'''
    flow_time_calcs(route)
    
    pipe_head_calcs(route)

    if route.minimum_head_tolerance > 0: #avoids a hydro station on very small residual head that is either just error in the calculations or can be handled by a simple control structure
        if route.generation_head <= route.minimum_head_tolerance: 
            route.generation_head = 0
        if route.pumping_head <= route.minimum_head_tolerance: 
            route.pumping_head = 0 

    route.hydro_installed_capacity = hydro_power(route.diversion_m3_per_s,route.generation_head,route.hydro_efficiency)
    route.hydro_grid = electricityroute.grid_cost(route.hydro_grid_distance_km,route.hydro_installed_capacity/1e6,include_transformer = route.include_transformer)
    route.pump_installed_capacity = pump_power(route.diversion_m3_per_s,route.pumping_head,route.pump_efficiency)
    route.pump_grid = electricityroute.grid_cost(route.pump_grid_distance_km,route.pump_installed_capacity/1e6,include_transformer = route.include_transformer)
    
    if not route.highpoint_storage:
        if route.hydro_installed_capacity < route.pump_installed_capacity:
            #pumps and hydro have to work at same time since there is more pumping then hydro buy the cheapest power
            route.power_sell_price =  route.power_buy_price
        else:
            #pumps and hydro have to work at same time since there is more hydro then pumping sell the expensive power
            route.power_buy_price =  route.power_sell_price
            
    route.hydro_direct_cost = hydro_power_station_cost(route.hydro_installed_capacity)
    route.hydro_annual_output = route.hydro_installed_capacity * route.flow_annual_hours
    route.hydro_annual_gross_income = route.hydro_annual_output * route.power_sell_price/1E6
    route.hydro_total_cost,route.hydro_overheads,_insurance,route.hydro_contingencies = economics.indirect_costs(route.hydro_direct_cost,route.hydro_overhead_rate,0,route.hydro_contingency_rate)
    
    route.pump_direct_cost = pump_station_cost(route.pump_installed_capacity)
    
    route.pump_annual_consumption = route.pump_installed_capacity * route.flow_annual_hours
    _proportion_behind_the_meter = min(1,route.behind_the_meter_hours_per_day/route.flow_hours_per_day)
    route.pump_behind_the_meter_annual_consumption =route.pump_annual_consumption * _proportion_behind_the_meter
    route.pump_metered_annual_consumption = route.pump_annual_consumption * (1-_proportion_behind_the_meter)
    route.pump_annual_tuos_cost = 0
    if route.use_tuos:
        route.pump_annual_tuos_cost = electricityroute.annual_tuos(route.pump_installed_capacity/1000,
                                                               route.pump_metered_annual_consumption/1000,route.tuos_locational_price)
    route.pump_annual_metered_cost = route.annual_connection_tariff + (route.pump_annual_tuos_cost + (route.pump_metered_annual_consumption * route.power_buy_price/1E6))*(1+route.pump_electricity_retailer_margin)
    route.pump_behind_the_meter_annual_cost = route.pump_behind_the_meter_annual_consumption * route.behind_the_meter_LCOE/1E6
    route.pump_total_cost,route.pump_overheads,_insurance,route.pump_contingencies = economics.indirect_costs(route.pump_direct_cost,route.pump_overhead_rate,0,route.pump_contingency_rate)
    
    route.water_purchase_cost = (route.water_source_cost_per_ML/1E6)* route.annual_diversions_l
    _pipe_cost = sum(p.total_cost for p in route.pipes)
    _pump_hydro_operating = operating_cost(route.hydro_installed_capacity+route.pump_installed_capacity)
    _pipe_operating_cost = pipe_o_and_m_cost(_pipe_cost)
    route.operating_cost = _pump_hydro_operating + _pipe_operating_cost
    route.total_cost = _pipe_cost + route.hydro_total_cost + route.pump_total_cost + route.pump_grid.total_cost + route.hydro_grid.total_cost+route.other_costs
    #route.total_cost,route.overheads,_insurance,route.contingencies = economics.indirect_costs(route.direct_cost)
    _total_annual_expenses = route.pump_behind_the_meter_annual_cost + route.pump_annual_metered_cost + route.operating_cost + route.water_purchase_cost + route.other_operating_costs
    route._annual_cash_flows = [-route.total_cost,*np.repeat(route.hydro_annual_gross_income - _total_annual_expenses, route.lifetime)]
    route.total_npv = npf.npv(route.discount_rate,route._annual_cash_flows)
    #route.total_irr = npf.irr(annual_cash_flows) #commented out as this is slow to calculate
    route.total_mirr = npf.mirr(route._annual_cash_flows,route.discount_rate,route.discount_rate)
    route.annualised_cost_per_ML = economics.lcoe(route.annual_diversions_l/1E6, -route.total_npv, 0 , route.discount_rate, route.lifetime)

    #LCOE gets tricky when there is also pumping in the system
    #need to include the NPV of the pumping into the costs for the LCOE calculation
    route._annual_cash_flows_sans_hydro = [-route.total_cost,*np.repeat(0 - _total_annual_expenses,route.lifetime)]
    route._total_npv_sans_hydro = npf.npv(route.discount_rate,route._annual_cash_flows_sans_hydro)
    route.total_LCOE = economics.lcoe(route.hydro_annual_output, -route._total_npv_sans_hydro, 0 , route.discount_rate, route.lifetime) #operating costs are 0 as already included in the npv calculation
    
    #marginal values compared to a base case
    if route.base_scenario is None:
        route.marginal_total_cost = None
        route.marginal_npv = None
        route.marginal_mirr = None
        route.marginal_LCOE = None
    else:
        route.marginal_total_cost = route.total_cost - route.base_scenario.total_cost
        route.marginal_npv = route.total_npv - route.base_scenario.total_npv
        marginal_annual_cash_flows = [x - y for x, y in zip(route._annual_cash_flows, route.base_scenario._annual_cash_flows)]
        #route.marginal_irr = npf.irr(marginal_annual_cash_flows) #commented out as this is slow to calculate
        route.marginal_mirr = npf.mirr(marginal_annual_cash_flows,route.discount_rate,route.discount_rate)

        #LCOE gets tricky when there is also pumping in the system
        #need to include the NPV of the pumping into the costs for the LCOE calculation
        marginal_annual_cash_flows_sans_hydro = [x - y for x, y in zip(route._annual_cash_flows_sans_hydro, route.base_scenario._annual_cash_flows_sans_hydro)]
        marginal_npv_sans_hydro = route._total_npv_sans_hydro - route.base_scenario._total_npv_sans_hydro
        route.marginal_LCOE = economics.lcoe(route.hydro_annual_output - route.base_scenario.hydro_annual_output, -marginal_npv_sans_hydro, 0, route.discount_rate, route.lifetime) #operating costs are 0 as already included in the npv calculation

def pipe_head_calcs(route:Pipe_Route):
    #make sure some things are zero

    route.total_head = 0
    #route.residual_head = 0
    route.pumping_head = 0
    #route.generation_head = 0
    #work thought the pipes sequentially
    available_head = 0
    for p in route.pipes:
        #print (p.pipe_type)
        local_head = p.start_elev - p.end_elev
        route.total_head += local_head
        p.initial_head = available_head + local_head
        if p.pipe_type == 'tunnel':
            p.direct_cost = tunnel_cost(p.D,p.L,p.hard_rock_ratio,1-p.hard_rock_ratio)
        else: #assume GRP
            p.direct_cost = sum(pipe_diameter_cost(p.D,p.L,p.num_barrels))
        p.total_cost,p.overheads,_insurance,p.contingencies = economics.indirect_costs(p.direct_cost,p.overhead_rate,0,p.contingency_rate)

        p.V = volume_to_velocity(route.diversion_m3_per_s/p.num_barrels,D=p.D)
        p.residual_head = residual_head(p.initial_head,p.V,p.D,p.L)
        if p.residual_head < 0: #todo: think about the velocity head
            route.pumping_head += abs(p.residual_head)
            available_head = 0 #if the residual is negative pumping will bring it back to zero
        elif p.residual_head > 0:
            #any residual head can be scrubbed in the next pipes
            available_head = p.residual_head
    route.generation_head = max(0,route.pipes[-1].residual_head)
    route.residual_head = route.generation_head - route.pumping_head
    
def pipefcn(args,pr,return_obj = False):
    #print (args)
    pr2 = copy.deepcopy(pr)
    for k,D in enumerate(args):
        pr2.pipes[k].D = D
    route_calcs(pr2)
    if return_obj:
        return pr2
    #flip the sign is asked for so the result is maximised by the minimise function
    #print (getattr(pr, pr.optimise_on))
    if pr2.optimise_abs:
        opt_val = abs(getattr(pr2, pr2.optimise_on))
    else:
        opt_val = getattr(pr2, pr2.optimise_on)
    if pr2.optimise_minimise:
        return opt_val
    else:
        return -opt_val

def pipe_optimise(pr: Pipe_Route):
    '''a simple wrapper to set up for the scipy optimise function
    '''
    pr = copy.deepcopy(pr)

    guess = [p.max_D for p in pr.pipes]
    bounds = [(p.min_D,p.max_D) for p in pr.pipes]
    args = pr
    #return (guess,bounds,args)
    #print (bounds)
    #r = minimize(pipefcn,guess,bounds=bounds)
    r = minimize(pipefcn,guess,bounds=bounds,args=args,method='Powell')#,options={'ftol':0.01,'xtol':0.001})
    #r = differential_evolution(pipefcn,bounds=bounds)
    #print (guess,pr.pipes[0].num_barrels)
    if not r.success:
        #print (r,guess,bounds,args)
        raise RuntimeError (r.message)
    return pipefcn(r.x,pr,return_obj = True)

def best_pipe(pr: Pipe_Route,return_all_barrels = False): 
    def barrel_setup():
        num_barrels_list = list(num_barrels)
        for p in barrel_pr.pipes:
            if p.pipe_type=='pipe':
                p.num_barrels = num_barrels_list.pop()
                p.max_D = min(p.max_D, velocity_to_D(p.min_V,barrel_pr.diversion_m3_per_s/p.num_barrels))
                p.min_D = max(p.min_D, velocity_to_D(p.max_V,barrel_pr.diversion_m3_per_s/p.num_barrels))
                if p.min_D > p.max_D:
                    #no pipe diameters that meet the velocity and diameter restrictions
                    return
            #    #then there are pipe diameters that meet the velocity and diameter restrictions
            #
            #    #start with the fully scrubbed case ( the optimiser sometimes doesnt pick this up)
            #    pipe_scrubbed = copy.copy(barrel_p)
            #    pipe_scrubbed.D = diameter_for_Q(pipe_scrubbed.total_head,pipe_scrubbed.diversion_m3_per_s/pipe_scrubbed.num_barrels,pipe_scrubbed.L)
            #    print (f'{num_barrels} {pipe_scrubbed.D} scrubbed' )
            #    if pipe_scrubbed.min_D <= pipe_scrubbed.D <= pipe_scrubbed.max_D:
            #        pipe_hydro_cost(pipe_scrubbed,clamp_residual_head=True)
            #        results.append(pipe_scrubbed)
            #        barrel_p.min_D = pipe_scrubbed.D # use this line to avoid pumping
            #        print (f'{num_barrels} {pipe_scrubbed.D} scrubbed' )
            #    #then optimise diameter for the least cost case

        optimum_pr = pipe_optimise(barrel_pr)
        results.append(optimum_pr)
        #print (f'{num_barrels} {optimum_p.D} optimum' )

    results = []
    #optimise doesnt handle integers so doing barrels manually
    num_of_pipes = len([p for p in pr.pipes if p.pipe_type=='pipe'])
    c = list(itertools.product(*itertools.repeat(range(pr.min_num_barrels,pr.max_num_barrels+1),num_of_pipes)))
    for num_barrels in c:
        #print (num_barrels)
        barrel_pr = copy.deepcopy(pr)
        barrel_setup()
    if return_all_barrels:
        return results
    #now pick the best barrel configuration
    if pr.optimise_minimise:
        return min(results, key=attrgetter(pr.optimise_on))
    else:
        return max(results, key=attrgetter(pr.optimise_on))

def optimise_route(route,best_buy_value,best_sell_value):
    route = copy.deepcopy(route)
    try:
        flow_time_calcs(route)
    except:
        print ('24hour limit')
        return route
    
    if route.use_average_day_power_prices:
        route.power_buy_price = np.interp(route.flow_hours_per_day, best_buy_value.index, best_buy_value)
        route.power_sell_price = np.interp(route.flow_hours_per_day, best_sell_value.index, best_sell_value)
        #used for tully
    
    route.water_source_cost_per_ML = route.power_buy_price * route._MW_per_ML
    
    route.optimise_on = 'total_npv'
    route.optimise_abs = False
    route.optimise_minimise=False
    route.base_scenario = None
    best_npv = best_pipe(route)
    #results_best_npv.append(best_npv)
    return best_npv

    #route.optimise_on = 'marginal_LCOE'
    #route.optimise_abs = False
    #route.optimise_minimise=True
    ##route.highpoint_storage = True
    #route.base_scenario = best_npv
    ##route.base_case_npv = best_npv.total_npv
    ##route.base_case_hydro_annual_output = best_npv.hydro_annual_output
    #marginal_best_lcoe = best_pipe(route)
    ##results_marginal_best_lcoe.append(marginal_best_lcoe)
    #return marginal_best_lcoe
    ##return (best_npv,marginal_best_lcoe)