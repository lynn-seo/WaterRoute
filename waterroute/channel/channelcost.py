import math
from dataclasses import dataclass, field, asdict
import pickle
from scipy.optimize import root
import pandas as pd
import numpy as np
import numpy_financial as npf

from waterroute.channel import bridges
from waterroute.channel import cross_drainage
from waterroute.channel import inlineAndSiphonControlStructures
from waterroute.channel import earthworks
from waterroute.channel import Losses
from waterroute.channel import channelLinearCosts
from waterroute import piperoute
from waterroute import utils
from waterroute import channelroute
from waterroute import economics


@dataclass
class Channel_Reach:
    start_diversion_m3_per_s: float = 0
    end_diversion_m3_per_s: float = 0
    start_diversions_annual_l: float = None
    end_diversions_annual_l: float = None
    start_location_id: int = 0
    end_location_id: int = 0
    start_distance: float = 0
    end_distance: float = 0
    length_km: float = 0
    bed_width: float = 0
    earthworks = None
    water_section_area: float = None
    wetted_perimeter: float = None
    water_surface_width: float = None
    peak_loss = None
    annual_loss = None
    access_road_cost:float=0
    mean_toe_to_toe_width:float=None
    reserve_width_m:float=None
    
@dataclass
class Channel_Route:
    reaches:list  = field(default_factory=list)
    description: str = ''
    linear_cost_factor: float = 1.0
    loss_factor: float = 1.0
    hydraulic_loss_factor: float = 1
    cross_drainage_cost_factor: float = 1
    mannings_n: float = 0.0225

    min_slope: float = 1/20000
    
    path = None
    profile = None
    lined = True
    rockDepth = 5
    max_reach_length_m = 200E3
    rolling_shutdown_factor: float = 2/3
    access_road_direct_cost_per_km:float = 0.4*1E6 #high speed gravel road
    #from Kev: Reserve Width.  The only fault with toe to toe plus clearance is that it won’t allow for catch drainage requirements.  So, if you say toe to toe plus 5m one side, and say 15 the other, I think you will be safe.  Will be conservative probably though, as catch drainage won’t always be necessary.  
    additional_reserve_width_m:float= 15 + 5 #this is in addition to the average toe to toe width
    land_resumptions_and_relocations_per_ha:float = 5000 #a very generous valuation on land to cover land value,construction disruption and ongoing inconvenience
    railway_crossing_count:int = 0
    highway_crossings_count:int = 0
    distance_between_farm_crossing_km:float = 4 #default was 2km but this is excessive in remote sections of the channel
    large_terminal_storage:bool = True # used to switch off some losses
    depth_scale:float = None
    max_water_depth:float  = 4.4
    freeboard:float=None
    OandM_rate_per_km:float = 60000
    OandM_cost:float = 0
    construction_time_years:float = 0
    earthworks_volume_m3_per_year:float = 65 *1E6 #m3 of material moved/placed each year
    discount_rate:float = 0.07
    lifetime:float = 100
    annualised_cost_per_saleable_ML:float = 0
        
def costs_calc(CR):
    ts_factor = 1
    if CR.large_terminal_storage:
        ts_factor = 0

    CR.total_peak_loss = {}
    CR.total_annual_loss = {}
    CR.farm_crossing_costs = 0

    CR.profile['reach_id'] = np.floor(CR.profile['short_distance(m)']/(CR.max_reach_length_m)).astype(int)
    for step in range(0,CR.profile['reach_id'].max()+1):
        print (step)
        CR.reaches.append(Channel_Reach())
        reach = CR.reaches[-1]
        if len(CR.reaches) == 1:
            reach.start_location_id = 0
            reach.start_distance = 0
            reach.start_diversion_m3_per_s = CR.gross_diversions_annual_m3_per_s
            reach.start_diversions_annual_l = CR.gross_diversions_annual_l
        else:
            reach.start_location_id = CR.reaches[-2].end_location_id + 1
            reach.start_distance = CR.reaches[-2].end_distance
            reach.start_diversion_m3_per_s = CR.reaches[-2].end_diversion_m3_per_s
            reach.start_diversions_annual_l = CR.reaches[-2].end_diversions_annual_l

        #size the channel to the peak discharge
        # reach.bed_width = Losses.channel_width_for_discharge(reach.start_diversion_m3_per_s,CR.max_water_depth,CR.depth_scale)
        # print('water_width',CR.routing_bed_width, reach.bed_width)
        reach.bed_width = CR.routing_bed_width
        # reach.water_depth = Losses.scaled_depth(reach.bed_width,max_flow_depth = CR.max_water_depth,depth_scale = CR.depth_scale)
        # print('water_depth',CR.routing_water_depth, reach.water_depth)
        reach.water_depth = CR.routing_water_depth
        
        #change above to refer the depth and width dected in channelroute.py, solver, def find_channel_dimension

        reach.freeboard = CR.freeboard
        if CR.freeboard is None:
            reach.freeboard = 0.2 * reach.water_depth + 0.3

        #lots of default values here
        reach.earthworks = earthworks.getCostEarthworks(DFL = 0,bedWidth = reach.bed_width, depthFlow = reach.water_depth, 
                                                        freeboard = reach.freeboard,channelLined = CR.lined , rockDepth = CR.rockDepth, 
                                                        sides = CR.side_slope)
        reach.water_section_area = reach.earthworks['Geometries']['waterSection'].area
        reach.wetted_perimeter = reach.earthworks['Geometries']['wettedPerimeterSection'].length
        reach.water_surface_width = reach.earthworks['Geometries']['waterSection'].bounds[2] - reach.earthworks['Geometries']['waterSection'].bounds[0]

        #work out where to end the reach - at a control point.
        CR.water_section_area = reach.water_section_area
        CR.wetted_perimeter = reach.wetted_perimeter
        CR.water_surface_width = reach.water_surface_width

        reach.end_location_id = CR.profile.index[CR.profile['reach_id']==step][-1]
        reach.end_distance = CR.profile.iloc[reach.end_location_id,:]['short_distance(m)']

        reach.length_km = (reach.end_distance - reach.start_distance)/1000
        #calculate losses.
        reach.peak_loss = Losses.getPeakDailyLossesLinear(annualDuration = 330, channelLength = reach.length_km, surfaceWidth = reach.water_surface_width, evapPanAnnual = 2,
                        evapPeakDaily = 6, bedPerimeter = reach.bed_width, k = 0.002, 
                        lossEachStockDomestic = 2, spacingStockDomestic = 2, leakageRate = 0, waterDepth=reach.water_depth)
        
        reach.annual_loss = Losses.getAnnualLossesLinear(annualDuration = 330, channelLength = reach.length_km, surfaceWidth = reach.water_surface_width, evapPanAnnual = 2,
                        evapPeakDaily = 6, bedPerimeter = reach.bed_width, k = 0.002, 
                        lossEachStockDomestic = 2, spacingStockDomestic = 2, leakageRate = 0, channelXSection = reach.water_section_area,rolling_shutdown_factor = CR.rolling_shutdown_factor, 
                        waterDepth=reach.water_depth)
        
        for k in reach.peak_loss.keys():
            CR.total_peak_loss[k] = CR.total_peak_loss.get(k,0) + reach.peak_loss[k]
        for k in reach.annual_loss.keys():
            CR.total_annual_loss[k] = CR.total_annual_loss.get(k,0) + reach.annual_loss[k]

        reach.end_diversion_m3_per_s = reach.start_diversion_m3_per_s - Losses.ML_day_to_cumecs(reach.peak_loss['Total'])
        reach.end_diversions_annual_l = reach.start_diversions_annual_l - reach.annual_loss['Total']*1E6

        #earthworks costs for reach
        p =  CR.profile.loc[reach.start_location_id:reach.end_location_id]
        # print(CR.profile.loc[reach.start_location_id:reach.end_location_id,'rock_depth(m)'])
        # print(CR.rockDepth)
        rd = p['rock_depth(m)'].mean()
        s = p['flow_level_elevation(m)'] - p['natural_surface_elevation(m)']
        ew = s.apply(earthworks.getCostEarthworks,bedWidth = reach.bed_width, depthFlow = reach.water_depth, 
                     freeboard = reach.freeboard,channelLined = CR.lined , rockDepth = rd, sides = CR.side_slope)
        ew = pd.DataFrame(ew.tolist())
        ew.index = s.index
        
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'full_width'] = ew['full_width']
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'earthworks_direct_cost_per_m'] = ew['Total Cost'] * CR.earthwork_cost_weight
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'earthworks_rock_excav_vol_per_m'] = ew['Rock Excav']
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'earthworks_soil_excav_vol_per_m'] = ew['Excav']
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'earthworks_fill_vol_per_m'] = ew['Fill-Normal']
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'earthworks_borrow_fill_vol_per_m'] = ew['Fill-Borrow']
        
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'bed_width'] = reach.bed_width
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'water_depth'] = reach.water_depth
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'freeboard'] = reach.freeboard
        CR.profile.loc[reach.start_location_id:reach.end_location_id,'capacity_m3_per_s'] = reach.start_diversion_m3_per_s

        reach.farm_crossing_costs = (reach.length_km/CR.distance_between_farm_crossing_km) * bridges.getTotalCostPerCrossing(BW = reach.bed_width, D = reach.water_depth+reach.freeboard, crossingType='farm access')['Total']
        CR.farm_crossing_costs += reach.farm_crossing_costs

    #this is the last reach
    #finalise losses
    CR.end_losses_Ml = Losses.getAnnualLossesDamToChannel(diversionFromDam = reach.end_diversions_annual_l/1E9 , diversionEfficiency = ts_factor*100, distributionEfficiency = ts_factor*100)
    CR.total_loss_GL = CR.total_annual_loss['Total']/1000 + CR.end_losses_Ml['Total']/1000
    CR.saleable_allocation_GL = (CR.gross_diversions_annual_l/1E9) - CR.total_loss_GL
    CR.losses_percent = CR.total_loss_GL/(CR.gross_diversions_annual_l/1E9)*100


    CR.profile['earthworks_direct_cost'] = CR.profile['earthworks_direct_cost_per_m'] * CR.profile['cell_distance(m)']
    CR.profile['earthworks_rock_excav_vol'] = CR.profile['earthworks_rock_excav_vol_per_m'] * CR.profile['cell_distance(m)']
    CR.profile['earthworks_soil_excav_vol'] = CR.profile['earthworks_soil_excav_vol_per_m'] * CR.profile['cell_distance(m)']
    CR.profile['earthworks_fill_vol'] = CR.profile['earthworks_fill_vol_per_m'] * CR.profile['cell_distance(m)']
    CR.profile['earthworks_borrow_fill_vol'] = CR.profile['earthworks_borrow_fill_vol_per_m'] * CR.profile['cell_distance(m)']
    
    
    CR.unit_farm_crossing_costs = bridges.getTotalCostPerCrossing(BW = CR.reaches[0].bed_width, D = CR.reaches[0].water_depth+CR.reaches[0].freeboard,
                                                                   crossingType='farm access')['Total'] 
    #major bridges costs - dont have a location for the bridges so use the dimensions for the widest part of the channel
    CR.railway_crossing_costs = CR.railway_crossing_count * bridges.getTotalCostPerCrossing(BW = CR.reaches[0].bed_width, D = CR.reaches[0].water_depth+CR.reaches[0].freeboard, crossingType='railway')['Total']
    CR.highway_crossing_costs = CR.highway_crossings_count * bridges.getTotalCostPerCrossing(BW = CR.reaches[0].bed_width, D = CR.reaches[0].water_depth+CR.reaches[0].freeboard, crossingType='highway')['Total']

    CR.earthworks_direct_cost = CR.profile['earthworks_direct_cost'].sum()
    CR.earthworks_rock_excav_vol = CR.profile['earthworks_rock_excav_vol'].sum()
    CR.earthworks_soil_excav_vol = CR.profile['earthworks_soil_excav_vol'].sum()
    CR.earthworks_fill_vol = CR.profile['earthworks_fill_vol'].sum()
    CR.earthworks_borrow_fill_vol = CR.profile['earthworks_borrow_fill_vol'].sum()
    #costs
    CR.length_km = CR.profile['short_distance(m)'].max()/1000
    
    CR.unit_access_road_cost = CR.access_road_direct_cost_per_km * CR.linear_cost_factor
    CR.access_road_cost = CR.access_road_direct_cost_per_km * CR.length_km * CR.linear_cost_factor
    CR.mean_toe_to_toe_width =  CR.profile['full_width'].mean() 
    
    CR.unit_farm_crossing_cost = bridges.getTotalCostPerCrossing(BW = CR.reaches[0].bed_width, D = CR.reaches[0].water_depth+CR.reaches[0].freeboard, 
                                                                 crossingType='farm access')['Total'] 
    CR.reserve_width_m = CR.mean_toe_to_toe_width + CR.additional_reserve_width_m
    _linear_cost_per_km = channelLinearCosts.getTotalChannelLinearCost(reserveWidth = CR.reserve_width_m, noGates = 1)

    CR.unit_easement_linear_costs_km = _linear_cost_per_km['Easement'] * CR.linear_cost_factor
    CR.easement_linear_costs =_linear_cost_per_km['Easement']  * CR.length_km * CR.linear_cost_factor
    
    CR.unit_catch_drainage_linear_costs_km = _linear_cost_per_km['Catch drainage'] * CR.linear_cost_factor
    CR.catch_drainage_costs =_linear_cost_per_km['Catch drainage']  * CR.length_km * CR.linear_cost_factor
    CR.land_resumptions_and_relocations_cost = CR.land_resumptions_and_relocations_per_ha * CR.reserve_width_m * CR.length_km * 1000 / 10000
    CR.SCADA_network_cost = CR.length_km * 3400 * CR.linear_cost_factor #from Kev, 2025

    CR.siphon_cost = 0
    CR.siphon_earthworks_cost_saving = 0
    CR.total_siphon_length = 0
    CR.cross_drain_cost = 0
    CR.profile['use_siphon'] = False
    CR.profile['use_cross_drain'] = False
        
    #insert to code to debug relavant to siphon_cost
    # with open ('CR.profile.pkl', 'wb') as f:
    #     f.write(pickle.dumps(CR.profile))

    if 'q40' in CR.profile.columns and (CR.profile['q40'] > 0).any():
        #costs for cross drains at all streams - some of these might not be used as they will be siphons
        #cross_drain_cost,number_of_barrels,diameter = cross_drainage.select_cross_drain(cross_drain_discharge, cross_drain_length)
        df = CR.profile[CR.profile['q40']>0].apply(lambda x: cross_drainage.select_cross_drain(x['q40'],x['full_width']),axis='columns', result_type='expand')
        df.rename(columns={0:'cross_drain_cost', 1:'cross_drain_number_of_barrels', 2:'cross_drain_diameter'},inplace = True)
        CR.profile =  CR.profile.merge(df,left_index=True, right_index=True,how='left')

        #costs for siphons at all streams - some of these might not be used as they will be cross drains
        #direct_cost,actual_D,num_barrels,V = inlineAndSiphonControlStructures.select_siphon(bedWidth,waterDepth,freeboard,capacity_m3_s,siphon_length)
        df = CR.profile[CR.profile['q40']>0].apply(lambda x: inlineAndSiphonControlStructures.select_siphon(x['bed_width'],x['water_depth'],x['freeboard'],x['capacity_m3_per_s'],x['length_of_siphon(m)']),axis='columns', result_type='expand')
        df.rename(columns={0:'siphon_cost', 1:'siphon_pipe_D', 2:'siphon_number_of_barrels', 3:'siphon_velocity'},inplace = True)
        CR.profile =  CR.profile.merge(df,left_index=True, right_index=True,how='left') 
        
        #If cross drain capacity larger than 6 times of channel flow, use siphon else culvert 
        # CR.profile['type_cross_drain'] = CR.profile['type_cross_drain'].astype(object)
        # CR.profile.loc[CR.profile['q40'] > (6 * CR.profile['capacity_m3_per_s']), 'type_cross_drain'] = 'Siphon'
        
        #not entierly accurate but close enough i think.
        CR.profile['siphon_earthworks_cost_saving'] = CR.profile['earthworks_direct_cost_per_m'] * CR.profile['length_of_siphon(m)']

        # #siphon structure cost 
        df = CR.profile[CR.profile['q40']>0].apply(lambda x: inlineAndSiphonControlStructures.simplified_control_structures_cost(x['bed_width'],x['water_depth'],x['freeboard'], x['siphon_pipe_D'], x['siphon_number_of_barrels']),axis='columns', result_type='expand')
        # df.rename(columns={0:'siphon_control_structure_cost'},inplace = True)
        # CR.profile =  CR.profile.merge(df,left_index=True, right_index=True,how='left') 
        if isinstance(df, pd.Series):
            df = df.rename('siphon_control_structure_cost').to_frame()
        else:
            # df is already a DataFrame
            if 0 in df.columns:
                df = df.rename(columns={0: 'siphon_control_structure_cost'})
            elif len(df.columns) == 1:
                df.columns = ['siphon_control_structure_cost']
        
        CR.profile =  CR.profile.merge(df,left_index=True, right_index=True,how='left')
        CR.profile['siphon_cost'] = CR.profile['siphon_cost']*CR.cross_drainage_cost_factor
        CR.profile['cross_drain_cost'] = CR.profile['cross_drain_cost'] * CR.cross_drainage_cost_factor
        
        _net_siphon_cost = CR.profile['siphon_cost'] - CR.profile['siphon_earthworks_cost_saving'] + CR.profile['siphon_control_structure_cost']
        # siphon_cost already considered siphon structure cost
        # _net_siphon_cost = CR.profile['siphon_cost'] - CR.profile['siphon_earthworks_cost_saving'] 
                
        CR.profile['use_siphon'] = _net_siphon_cost < CR.profile['cross_drain_cost']
        CR.profile['use_cross_drain'] = _net_siphon_cost >= CR.profile['cross_drain_cost']
        CR.siphon_cost = CR.profile[CR.profile['use_siphon']]['siphon_cost'].sum()
        # CR.siphon_cost = _net_siphon_cost #TODO delete
        CR.siphon_earthworks_cost_saving = CR.profile[CR.profile['use_siphon']]['siphon_earthworks_cost_saving'].sum()
        CR.cross_drain_cost = CR.profile[CR.profile['use_cross_drain']]['cross_drain_cost'].sum()
        CR.total_siphon_length = CR.profile[CR.profile['use_siphon']]['length_of_siphon(m)'].sum()
    
    #now control_structure refers to actual slope the model decided
    temp_slope = (CR.profile['flow_level_elevation(m)'].iloc[-1] - CR.profile['flow_level_elevation(m)'].iloc[0]) / (CR.profile['short_distance(m)'].iloc[-1] - CR.profile['short_distance(m)'].iloc[0])
    CR.profile = channelroute.add_control_structures_V2(CR.profile,maxControlSpacing=CR.maxControlSpacing,channelSlope= abs(temp_slope),minWaterDepth=CR.bay_minWaterDepth)
    
    #use_siphon_v2 is the final decision on whether to use a siphon or not, 2/10/2025 Lynn
    #this is culvert 
    CR.profile['control_structure_cost']  = CR.profile[(CR.profile['control_structure_id']>0) & ~CR.profile['use_siphon']].apply(lambda x: inlineAndSiphonControlStructures.simplified_control_structures_cost(x['bed_width'],x['water_depth'],x['freeboard']),axis='columns')
    CR.control_structure_cost = CR.profile['control_structure_cost'].sum()
    
    siphon_rows = CR.profile[CR.profile['use_siphon']]
    df_siphonLosses = pd.DataFrame(index=siphon_rows.index, columns=['siphon_head_loss_single'], dtype=float)
    if not siphon_rows.empty:
        df_siphonLosses['siphon_head_loss_single'] = siphon_rows.apply(lambda x: inlineAndSiphonControlStructures.siphon_headLoss(x['length_of_siphon(m)'], 0.013, x['siphon_velocity'], x['siphon_pipe_D']),axis='columns')
    CR.profile =  CR.profile.merge(df_siphonLosses,left_index=True, right_index=True,how='left')
    
    df_siphonLosses_v2 = pd.DataFrame(index=siphon_rows.index, columns=['siphon_head_loss_allBarrels'], dtype=float)
    if not siphon_rows.empty:
        df_siphonLosses_v2['siphon_head_loss_allBarrels'] = siphon_rows.apply(lambda x: inlineAndSiphonControlStructures.siphon_headLoss_v2(x['length_of_siphon(m)'], 0.013, x['capacity_m3_per_s'], x['siphon_pipe_D'], x['siphon_number_of_barrels']),axis='columns')
    CR.profile =  CR.profile.merge(df_siphonLosses_v2,left_index=True, right_index=True,how='left')

    df_siphon_inletOutlet_Losses = pd.DataFrame(index=siphon_rows.index, columns=['siphon_inletOutlet_loss'], dtype=float)
    if not siphon_rows.empty:
        df_siphon_inletOutlet_Losses['siphon_inletOutlet_loss'] = siphon_rows.apply(lambda x: inlineAndSiphonControlStructures.siphon_inletOutlet_loss(x['siphon_velocity']),axis='columns')
    CR.profile =  CR.profile.merge(df_siphon_inletOutlet_Losses,left_index=True, right_index=True,how='left')
    
    #add extra costs for bays (skip first control structure)
    CR.profile.loc[CR.profile['control_structure_id']>1,'weed_excluder_cost'] = 20000
    CR.profile.loc[CR.profile['control_structure_id']>1,'drainage_structures_cost'] = 100000
    CR.profile.loc[CR.profile['control_structure_id']>1,'bank_overflows_cost'] = 81730
    
    CR.weed_excluder_cost = CR.profile['weed_excluder_cost'].sum()
    CR.drainage_structures_cost = CR.profile['drainage_structures_cost'].sum()
    CR.bank_overflows_cost = CR.profile['bank_overflows_cost'].sum()

    #sum up items
    # CR.semi_linear_cost = sum([CR.farm_crossing_costs,
    #                      CR.easement_linear_costs,
    #                      CR.catch_drainage_costs,
    #                      CR.SCADA_network_cost,
    #                      CR.control_structure_cost,
    #                      CR.weed_excluder_cost,
    #                      CR.drainage_structures_cost,
    #                      CR.bank_overflows_cost,
    #                     ])
    #Lynn splitted semin_linear_cost into linear cost and semi linear cost
    #potentially add access road cost in linear cost
    CR.linear_cost = sum([CR.access_road_cost,
                         CR.farm_crossing_costs,
                         CR.easement_linear_costs,
                         CR.catch_drainage_costs,
                         CR.SCADA_network_cost])*CR.linear_cost_factor
    
    #within_channel_infrastructure_cost
    # CR.semi_linear_cost = sum([CR.control_structure_cost,
    #                      CR.weed_excluder_cost,
    #                      CR.drainage_structures_cost,
    #                      CR.bank_overflows_cost,
    #                     ])
    CR.semi_linear_cost = sum([
                         CR.weed_excluder_cost,
                         CR.drainage_structures_cost,
                         CR.bank_overflows_cost,
                        ])

    CR.non_linear_costs = sum([CR.earthworks_direct_cost,
                        CR.railway_crossing_costs,
                        CR.highway_crossing_costs,
                        CR.siphon_cost,
                        -CR.siphon_earthworks_cost_saving,
                        CR.cross_drain_cost,
                        CR.control_structure_cost,
                        ])
    
    CR.total_major_items_cost = sum([CR.linear_cost,
                                     CR.semi_linear_cost,                                    
                                     
                                     CR.non_linear_costs,
                                    ])
    
    CR.minor_items_cost = 0.1 * CR.total_major_items_cost

    CR.total_direct_costs = CR.minor_items_cost + CR.total_major_items_cost

    # CR.overheads = 0.3 * CR.total_direct_costs
    #onsite cost
    CR.total_onsite_overheads = 0.12 * CR.total_direct_costs #changed from 0.3 to 0.12, 30/09/2025
    CR.profit = 0.1 * (CR.total_direct_costs + CR.total_onsite_overheads)
    CR.total_outturn_costs = CR.total_direct_costs + CR.total_onsite_overheads + CR.profit
    # CR.contingencies = 0.4 * sum([CR.land_resumptions_and_relocations_cost,
    #                         CR.overheads,
    #                         CR.accom_camps_costs,
    #                         CR.total_direct_costs,
    #                        ])
    #contigency
    CR.contingencies = 0.4 * sum([CR.total_onsite_overheads,CR.total_direct_costs]) #changged on 30/09/2025 to match overheads change
    #Owners costs
    CR.investigation_and_design = 0.05 * CR.total_direct_costs
    CR.accom_camps_costs = 0.15 * CR.total_direct_costs
    CR.insurance = 0.0125 * CR.total_outturn_costs
    CR.approval = 0.05 * CR.total_direct_costs
    CR.total_owner_costs = sum([CR.investigation_and_design,CR.accom_camps_costs,CR.insurance,CR.approval])

    CR.total_cost = sum([CR.total_outturn_costs, CR.contingencies, CR.total_owner_costs])
    
    #calculate a head budget
    #the assumptions we worked with last time on Control Structures and Access Crossings. 
    #We used:
    #Control Structures -  0.15m - related to flow measurement
    #Access Crossings  -  0.025m - as small as you can make i
    #Control Structures -  0.075m - related to flow measurement, based on Kev's advice, 2025
    #Access Crossings  -  0.025m - as small as you can make i

    _siphon_slope = 1/2000
    _length_of_channel_m = CR.length_km*1000 - CR.total_siphon_length
    _channel_head_loss_m = CR.real_channelSlope * _length_of_channel_m

    
    _siphon_head_loss_m = df_siphonLosses['siphon_head_loss_single'].sum()
    _siphon_head_loss_m_v2 = df_siphonLosses_v2['siphon_head_loss_allBarrels'].sum()
    _siphon_inlet_outlet_loss = df_siphon_inletOutlet_Losses['siphon_inletOutlet_loss'].sum()
    
    
    _control_head_loss_m = 0.075 * (CR.profile['control_structure_id']>0).sum() #a count

    #not including  (CR.railway_crossing_count  + CR.highway_crossings_count) as they would probably replace a farm bridge
    _crossing_head_loss_m = 0.025 *((_length_of_channel_m/1000)/CR.distance_between_farm_crossing_km)
    CR.calculated_total_head_loss_m = sum([_channel_head_loss_m,
                             _siphon_head_loss_m,
                             _control_head_loss_m,
                             _crossing_head_loss_m,])
    CR.channel_head_loss_m = _channel_head_loss_m
    CR.siphon_head_loss_m = _siphon_head_loss_m
    # CR.siphon_head_loss_m_v2 = _siphon_head_loss_m_v2
    CR.siphon_inlet_outlet_loss_m = _siphon_inlet_outlet_loss
    CR.control_head_loss_m = _control_head_loss_m
    CR.crossing_head_loss_m = _crossing_head_loss_m
    CR.average_slope_head_loss_m = CR.profile['flow_level_elevation(m)'].iloc[0] - CR.profile['flow_level_elevation(m)'].iloc[-1]
    
    CR.OandM_cost = CR.OandM_rate_per_km *CR.length_km
    
    CR.construction_time_years =     sum([CR.earthworks_rock_excav_vol,
                                          CR.earthworks_soil_excav_vol, #tempted to leave this out so not double counted with fill
                                          CR.earthworks_fill_vol,
                                          CR.earthworks_borrow_fill_vol,
                                         ]) / CR.earthworks_volume_m3_per_year
    
    #economics
    # _time = math.ceil(CR.construction_time_years)
    #was getting to tricky.
    #_cash_flows = [*np.repeat(-CR.total_cost/_time,_time),*np.repeat(CR.OandM_cost,CR.lifetime)]
    # _cash_flows = [-CR.total_cost,*np.repeat(CR.OandM_cost,CR.lifetime)]
    # CR.total_npv = npf.npv(CR.discount_rate,_cash_flows)
    # CR.total_mirr = npf.mirr(_cash_flows,CR.discount_rate,CR.discount_rate)
    #now the model uses constant channel dimensions; it refers to CR
    _total_loss_v2 = Losses.getAnnualLossesLinear(annualDuration = 330, channelLength = CR.length_km, surfaceWidth = CR.water_surface_width,
                                                    evapPanAnnual = 2, evapPeakDaily = 6, bedPerimeter = CR.wetted_perimeter, k = 0.002, 
                        lossEachStockDomestic = 2, spacingStockDomestic = 2, leakageRate = 0, channelXSection = CR.water_section_area, 
                        rolling_shutdown_factor = CR.rolling_shutdown_factor, waterDepth=CR.routing_water_depth)
    
    CR.total_loss_GL_v2 = _total_loss_v2['Total']/1000*CR.loss_factor
    CR.total_elevation_change_m = CR.profile['flow_level_elevation(m)'].iloc[0] - CR.profile['flow_level_elevation(m)'].iloc[-1]
    CR.total_length_m = CR.profile['short_distance(m)'].iloc[-1] - CR.profile['short_distance(m)'].iloc[0]
    CR.final_slope = (CR.profile['flow_level_elevation(m)'].iloc[-1] - CR.profile['flow_level_elevation(m)'].iloc[0]) / (CR.profile['short_distance(m)'].iloc[-1] - CR.profile['short_distance(m)'].iloc[0])

    

    CR.total_controlled_head_loss_m = CR.siphon_head_loss_m + CR.control_head_loss_m + CR.crossing_head_loss_m + CR.siphon_inlet_outlet_loss_m
    CR.total_controlled_head_loss_m = CR.total_controlled_head_loss_m * CR.hydraulic_loss_factor
    CR.adjusted_final_slope = (CR.profile['flow_level_elevation(m)'].iloc[-1]+CR.total_controlled_head_loss_m - CR.profile['flow_level_elevation(m)'].iloc[0]) / (CR.profile['short_distance(m)'].iloc[-1] - CR.profile['short_distance(m)'].iloc[0])
    _, V, A = channelroute.trapezoidal_discharge(width=CR.water_surface_width, depth=CR.routing_water_depth, S=abs(CR.adjusted_final_slope), n=0.03, z=2)
    CR.channel_velocity_m_per_s = V
    # CR.total_annual_loss_GL = CR.total_annual_loss['Total']/1E3
    CR.GL_delivered_per_year = CR.gross_diversions_annual_m3_per_s * 86400 / 1E6 *120 - CR.total_loss_GL_v2
    CR.dollar_per_GL_delivered = (CR.total_cost/1E6) / (CR.GL_delivered_per_year)
    annualised_capital_cost = -npf.pmt(CR.discount_rate, CR.lifetime, CR.total_cost)
    CR.annualised_cost_per_saleable_ML = (
        (annualised_capital_cost + CR.OandM_cost) / (CR.GL_delivered_per_year * 1000)
    )
    
