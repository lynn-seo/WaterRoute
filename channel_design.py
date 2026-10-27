"""
Channel design and routing functions for the WaterRoute notebook.

Import this before running the notebook:

    from channel_design import start_map, route_with_head_loss_iteration, candidate_designs

Importing it also enables KML output in fiona and silences the "Lossy conversion"
warning. The livemaptiles server is loaded only when adding a raster map layer.

Paths, hydraulics and per-section inputs are defined in the notebook and passed in
with `configure(...)` before running anything. Functions read them when called.
"""

import os
import warnings

import fiona
import geopandas as gpd
import ipyleaflet
import numpy as np
import pandas as pd
from matplotlib import cm
from pyproj import Transformer
from scipy.optimize import root_scalar

from waterroute import channelroute, utils
from waterroute.channel import channelcost

# Let geopandas/fiona write KML
fiona.supported_drivers['KML'] = 'rw'

warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Lossy conversion occurred.*")



def _array_tile_init(self, in_array, in_affine, in_crs):
    self.in_array = in_array
    self.in_affine = in_affine
    self.in_crs = in_crs
    self.tile_size = (256, 256)
    self.colourMap = cm.autumn
    self.scale_min = self.in_array.min()
    self.scale_max = self.in_array.max()
    self.alpha = None
    # EPSG:3857 = spherical mercator -> the array's CRS
    self.transformer = Transformer.from_crs(3857, self.in_crs, always_xy=True)

# ---- Settings (set from the notebook with configure) ----
SETTINGS = (
    'ROUTES_DIR', 'SEQUENTIAL_RUN_TEMPLATE', 'DEM_PATH', 'ROCK_DEPTH_PATH', 'ROUTING_BOUNDS',
    'MANNING_N', 'SIDE_SLOPE', 'SAMPLE_F', 'MAX_FILL_SIPHON_ENDS',
    'MAX_RUNS', 'RESIDUAL_HEAD_LOSS_TOL_M', 'CROSS_DRAINAGE_COST_FACTOR',
    'EARTHWORK_COST_FACTOR', 'HYDRAULIC_LOSS_FACTOR', 'LINEAR_COST_FACTOR','LOSS_FACTOR',
    'DISCOUNT_RATE', 'SECTIONS',
)
ROUTES_DIR = SEQUENTIAL_RUN_TEMPLATE = DEM_PATH = ROCK_DEPTH_PATH = ROUTING_BOUNDS = None
MANNING_N = SIDE_SLOPE = SAMPLE_F = MAX_FILL_SIPHON_ENDS = None
MAX_RUNS = RESIDUAL_HEAD_LOSS_TOL_M = CROSS_DRAINAGE_COST_FACTOR = SECTIONS = None
EARTHWORK_COST_FACTOR = 1.0
LINEAR_COST_FACTOR = 1.0
DISCOUNT_RATE = 0.07
HYDRAULIC_LOSS_FACTOR = 1.0
LOSS_FACTOR = 1.0

def configure(**settings):
    """Set the module settings listed in SETTINGS, e.g. configure(MANNING_N=0.03, ...)."""
    unknown = set(settings) - set(SETTINGS)
    if unknown:
        raise ValueError(f"Unknown settings {sorted(unknown)}; expected names from {SETTINGS}")
    globals().update(settings)


# ---- Map helpers ----
_map = None  # set by start_map(); layers are only drawn when a map exists


def remove_layer(the_map, name):
    for layer in the_map.layers:
        if layer.name == name:
            the_map.remove_layer(layer)


def add_layer(the_map, name, array, aff, crs, scale_min, scale_max, colour_map=cm.autumn):
    import livemaptiles

    # Pass a pyproj Transformer without skip_equivalent.
    livemaptiles.arrayTile.__init__ = _array_tile_init
    remove_layer(the_map, name)
    tile = livemaptiles.arrayTile(array, aff, crs)
    tile.scale_min = scale_min
    tile.scale_max = scale_max
    tile.colourMap = colour_map
    tile.colourMap.set_under(alpha=0)
    tile.colourMap.set_over(alpha=0)
    livemaptiles.tileServer.layers[name] = tile.tileMaker

    tile_layer = ipyleaflet.TileLayer(
        name=name,
        url=f'http://localhost:8080/livemaptiles/{name}/{{z}}/{{x}}/{{y}}.png',
    )
    the_map.add_layer(tile_layer)


def start_map(start=(145.829917, -24.693563)):
    """Create the interactive map (also used by runner). `start` is (lon, lat)."""
    global _map
    m = ipyleaflet.Map(center=start[::-1], zoom=8)
    m.scroll_wheel_zoom = True
    m.add_control(ipyleaflet.SearchControl(
        position='topleft',
        url='https://nominatim.openstreetmap.org/search?format=json&q={s}',
        property_name='display_name',
        marker=ipyleaflet.Marker(),
    ))
    m.add_control(ipyleaflet.FullScreenControl())
    m.add_control(ipyleaflet.LayersControl())
    m.add_control(ipyleaflet.MeasureControl(
        position='topleft',
        active_color='orange',
        primary_length_unit='kilometers',
        secondary_length_unit='meters',
        primary_area_unit='sqmeters',
    ))
    base = ipyleaflet.basemap_to_tiles(ipyleaflet.basemaps.OpenTopoMap)
    base.base = True
    m.add_layer(base)
    _map = m
    return m


def output_location(run_name):
    return utils.out_location(ROUTES_DIR, run_name)


# ---- Channel hydraulics (Manning, trapezoidal section) ----
def trapezoidal_discharge(width, depth, n=None, S=None, z=None):
    """
    Discharge for a trapezoidal channel using Manning's equation.

    Parameters
    ----------
    width : bottom width (m)
    depth : flow depth (m)
    S : channel slope (m/m)
    n : Manning roughness coefficient (default MANNING_N)
    z : side slope, horizontal : vertical (default SIDE_SLOPE)

    Returns
    -------
    Q : discharge (m3/s)
    V : velocity (m/s)
    A : cross-sectional area (m2)
    """
    n = MANNING_N if n is None else n
    S = 1/15000 if S is None else S
    z = SIDE_SLOPE if z is None else z
    A = depth * (width + z * depth)            # cross-sectional area
    P = width + 2 * depth * np.sqrt(1 + z**2)  # wetted perimeter
    R = A / P                                  # hydraulic radius
    V = (1 / n) * R ** (2 / 3) * S ** 0.5      # Manning velocity
    Q = A * V
    return Q, V, A


def solve_depth_ratio(Q, S, n, ratio, z=2):
    """Water depth (m) that carries Q at slope S when bed width = ratio * depth."""
    S = abs(S)

    def residual(depth):
        Q_calc, _, _ = trapezoidal_discharge(width=ratio * depth, depth=depth, S=S, n=n, z=z)
        return Q_calc - Q

    return root_scalar(residual, bracket=[0.1, 10], method='brentq').root



import numpy as np

def solve_S(Q, V, z=None, n=None, ratio=5):
    """Calculate the uniform-flow slope and trapezoidal channel dimensions.

    Parameters
    ----------
    Q : float
        Design discharge (m³/s).
    V : float
        Target mean flow velocity (m/s).
    z : float, optional
        Side slope, horizontal-to-vertical. Defaults to SIDE_SLOPE.
    n : float, optional
        Manning roughness coefficient. Defaults to MANNING_N.
    ratio : float, optional
        Bed-width-to-water-depth ratio. Defaults to 5.

    Returns
    -------
    tuple
        (S, depth, bed_width), with slope in m/m and dimensions in m.
    """
    n = MANNING_N if n is None else n
    z = SIDE_SLOPE if z is None else z

    if Q <= 0 or V <= 0 or n <= 0 or z < 0 or ratio <= 0:
        raise ValueError(
            "Q, V, n and ratio must be positive; z must be non-negative."
        )

    # Channel geometry from Q = V * A
    d = np.sqrt(Q / (V * (ratio + z)))
    b = ratio * d

    # Trapezoidal cross-section
    A = d * (b + z * d)
    P = b + 2 * d * np.sqrt(1 + z**2)
    R = A / P

    # Manning equation
    S = (V * n / R**(2 / 3))**2

    return S, d, b


def annualised_cost(cost, discount_rate=0.07, lifespan=50, om_fraction=0.02):
    """
    cost : capital cost
    discount_rate : e.g. 0.05 for 5%
    lifespan : asset life in years
    om_fraction : annual operation and maintenance as a fraction of capital cost
    """
    growth = (1 + discount_rate) ** lifespan
    capital_recovery = cost * discount_rate * growth / (growth - 1)
    return capital_recovery + cost * om_fraction


# ---- Run one route ----
def runner(CR, dem_ratio=1, dem=None, aff=None, crs=None):
    """Route a channel for the configured Channel_Route CR; returns the output path prefix."""
    if dem is None or aff is None or crs is None:
        dem, aff, crs = utils.load_raster(DEM_PATH, CR.routing_bounds, 1, dem_ratio, utils.Resampling.average)
    if _map is not None:
        add_layer(_map, 'dem', dem, aff, crs, 200, 1000, cm.viridis)
    rock_depth, _, _ = utils.load_raster(ROCK_DEPTH_PATH, CR.routing_bounds, 1, dem_ratio, utils.Resampling.average)

    output_path = output_location(CR.routing_out_name)


    #always insert/change parameters here to update the entire model
    paths = channelroute.elevRangePaths(
        output_path,
        dem,
        aff,
        crs,
        CR.routing_start_elevation,
        CR.routing_start_elevation + 1,
        targetSlope=CR.routing_targetSlope,
        start=CR.routing_start,
        finish=CR.routing_finish,
        max_fill_siphon_ends=MAX_FILL_SIPHON_ENDS,
        contour=CR.routing_contour,
        doCrossDrainage=CR.routing_doCrossDrainage,
        sampleF=SAMPLE_F,
        bed_width=CR.routing_bed_width,
        water_depth=CR.routing_water_depth,
        freeboard=CR.routing_freeboard,
        additional_per_m_cost=CR.routing_additional_per_m_cost,
        channel_capacity=CR.gross_diversions_annual_m3_per_s,
        rockDepth_raster=rock_depth,
        loss_value_dollar_per_ML=CR.loss_value_dollar_per_ML,
        total_linear_cost=CR.total_linear_cost,
        # total_loss_GL=None, #now this is calulcated inside the cahnnelroute.py, considering configured channel dimension based on potental channel length
        # no_cells=CR.no_cells,
        withLengthPenalty=CR.withLengthPenalty,
        cost_array_with_rock_depth=CR.cost_array_with_rock_depth,
        sunWater_route=CR.sunWater_route,
        earthwork_cost_weight = CR.earthwork_cost_weight,
        approximate_length=CR.approximate_length,
        discount_rate=CR.discount_rate,
        lifespan=CR.lifetime,
        linear_cost_factor=CR.linear_cost_factor,
        loss_factor=CR.loss_factor,
        # cross_drainage_cost_factor=CR.cross_drainage_cost_factor,
        # hydraulic_loss_factor=CR.hydraulic_loss_factor,
    )

    _, profile, line = paths[list(paths.keys())[0]]
    CR.profile = profile.copy()

    # Detailed cost, adjusting the distance of each reach
    channelcost.costs_calc(CR)
    
    CR.profile.to_csv(f'{output_path}_full_profile.csv')

    costs = CR.__dict__.copy()
    del costs['reaches']
    del costs['profile']
    pd.DataFrame([costs]).T.to_csv(f'{output_path}_costs.csv')

    if _map is not None:
        path_layer = ipyleaflet.GeoData(
            geo_dataframe=gpd.GeoDataFrame([line], columns=['geometry']),
            style={'color': 'blue', 'weight': 4},
        )
        path_layer.name = str(output_path)
        remove_layer(_map, path_layer.name)
        _map.add_layer(path_layer)

    return output_path


def run_channel(section='section1', ratio_width_depth=5, run_name=None, target_slope=None,
                bed_width=None, water_depth=None, linear_cost=None,
                loss_value_dollar_per_ML=75, is_length_penalty=True,
                cost_array_with_rock_depth=True, sunWater_route=False, 
                # earthwork_cost_weight=None, linear_cost_factor=None, loss_factor=None
                ):
    """Configure a Channel_Route for the given section and run it."""
    if section not in SECTIONS:
        raise ValueError(f"Unknown section {section!r}; expected one of {list(SECTIONS)}")
    cfg = SECTIONS[section]

    CR = channelcost.Channel_Route()

    CR.routing_bounds = ROUTING_BOUNDS
    CR.routing_out_name = run_name
    if target_slope is not None:
        CR.routing_targetSlope = -abs(target_slope)

    # Section-specific inputs
    CR.routing_start = cfg['start']
    CR.routing_finish = cfg['finish']
    CR.routing_start_elevation = cfg['start_elevation']
    CR.gross_diversions_annual_m3_per_s = cfg['gross_diversions_m3_s']
    CR.linear_cost_factor = LINEAR_COST_FACTOR 
    
    base_linear_cost = cfg['default_linear_cost'] if linear_cost is None else linear_cost
    CR.total_linear_cost = base_linear_cost * CR.linear_cost_factor
    # print(f"Linear cost: {base_linear_cost} x {CR.linear_cost_factor} = {CR.total_linear_cost}")
    # CR.total_loss_GL = cfg['total_loss_GL']
    CR.approximate_length = cfg['approximate_length_km']
    # CR.no_cells = cfg['no_cells']

    CR.loss_value_dollar_per_ML = loss_value_dollar_per_ML
    CR.cost_array_with_rock_depth = cost_array_with_rock_depth
    CR.withLengthPenalty = is_length_penalty
    CR.sunWater_route = sunWater_route
    CR.earthwork_cost_weight = EARTHWORK_COST_FACTOR
    CR.side_slope = SIDE_SLOPE
    

    # Channel dimensions
    CR.mannings_n = MANNING_N
    CR.ratio_width_depth = ratio_width_depth
    CR.routing_water_depth = water_depth
    CR.routing_bed_width = bed_width
    print(f"Calculated channel dimensions: slope = {CR.routing_targetSlope}, "
          f"mannings_n = {CR.mannings_n}, "
          f"side_slope = {CR.side_slope}, "
          f"ratio_width_depth = {CR.ratio_width_depth}, "
          f"water depth = {CR.routing_water_depth:.2f} m, bed width = {CR.routing_bed_width:.2f} m")

    # Routing options
    CR.routing_contour = 27
    CR.routing_doCrossDrainage = True
    CR.routing_freeboard = 1
    CR.routing_additional_per_m_cost = 2000
    CR.lined = False

    # Cost model
    CR.cross_drainage_cost_factor = CROSS_DRAINAGE_COST_FACTOR
    CR.hydraulic_loss_factor = HYDRAULIC_LOSS_FACTOR
    CR.loss_factor = LOSS_FACTOR
    CR.lifetime = 50
    CR.discount_rate = DISCOUNT_RATE
    CR.rockDepth = 5  # no longer used; kept in case channelcost still reads it
    # Based on 19.4 m3/s; may be superseded by the capacity-based delivery calculation
    CR.gross_diversions_annual_l = 608.645 * 1E9

    CR.railway_crossing_count = 0
    CR.highway_crossings_count = 0
    CR.distance_between_farm_crossing_km = 4
    CR.rolling_shutdown_factor = 2 / 3
    CR.large_terminal_storage = True
    CR.max_reach_length_m = 200000

    CR.maxControlSpacing = 28000
    CR.real_channelSlope = 1 / 12000
    CR.bay_minWaterDepth = 1

    return runner(CR=CR, dem_ratio=1)


# ---- Sequential runs ----
# Each run re-sizes the channel for the slope left over after head losses in the
# previous run, and re-routes. Results for every run are collected in a CSV.
def get_cost_scalar(df, key):
    v = df.loc[key]
    if isinstance(v, pd.Series):
        v = v.dropna().iloc[0] if not v.dropna().empty else np.nan
    return pd.to_numeric(v, errors="coerce")


def route_with_head_loss_iteration(number_of_runs=3, width_to_depth_ratio=None, water_depth=None, bed_width=None,
                                   run_name=None, section='section1', slope_for_channel_dimension=1/15000,
                                   linear_cost=None, loss_value_dollar_per_ML=75, is_length_penalty=True,
                                   cost_array_with_rock_depth=True, earthwork_cost_weight=None,
                                   linear_cost_factor=None, loss_factor=None):
    channel_capacity_m3_s = SECTIONS[section]['gross_diversions_m3_s']

    # 1 or 2 stop early; any larger value allows up to MAX_RUNS runs
    max_runs = number_of_runs if number_of_runs in (1, 2) else MAX_RUNS

    seq_df = pd.read_csv(SEQUENTIAL_RUN_TEMPLATE)
    seq_df['cost_file_name'] = seq_df['cost_file_name'].astype('string')

    prev_adjusted_slope = None
    prev_head_loss = None

    for run_id in range(1, max_runs + 1):
        row = seq_df['run_ID'].astype(str) == str(run_id)

        # Channel dimensions for this run
        if run_id == 1:
            slope, depth, width = slope_for_channel_dimension, water_depth, bed_width
        else:
            slope = prev_adjusted_slope
            # print(MANNING_N, "delete")
            depth = solve_depth_ratio(Q=channel_capacity_m3_s, S=slope, n=MANNING_N, ratio=width_to_depth_ratio, z = SIDE_SLOPE)
            width = width_to_depth_ratio * depth
        seq_df.loc[row, 'width_to_depth_ratio'] = width_to_depth_ratio
        seq_df.loc[row, 'slope'] = slope
        seq_df.loc[row, 'water_depth'] = depth
        seq_df.loc[row, 'bed_width'] = width

        # Route
        output_path = run_channel(
            section=section, ratio_width_depth=width_to_depth_ratio, run_name=run_name,
            target_slope=slope_for_channel_dimension, bed_width=width, water_depth=depth,
            linear_cost=linear_cost, 
            is_length_penalty=is_length_penalty, cost_array_with_rock_depth=cost_array_with_rock_depth,
            earthwork_cost_weight=earthwork_cost_weight,
            linear_cost_factor=linear_cost_factor,
            loss_factor=  loss_factor
        )

        # Collect results
        cost_file = f'{output_path}_costs.csv'
        cost_df = pd.read_csv(cost_file, index_col=0)
        surface_width = abs(get_cost_scalar(cost_df, 'water_surface_width'))
        _, velocity, _ = channelroute.trapezoidal_discharge(
            width=surface_width, depth=depth, S=slope, n=MANNING_N, z=SIDE_SLOPE)
        final_slope = abs(get_cost_scalar(cost_df, 'final_slope'))
        total_length = get_cost_scalar(cost_df, 'total_length_m')
        elevation_difference = get_cost_scalar(cost_df, 'total_elevation_change_m')
        head_loss = get_cost_scalar(cost_df, 'total_controlled_head_loss_m')

        if run_id == 1:
            residual = None
            adjusted_slope = abs(get_cost_scalar(cost_df, 'adjusted_final_slope'))
        else:
            # Positive residual head loss -> steeper slope; negative -> flatter slope
            residual = abs(prev_head_loss) - abs(head_loss)
            adjusted_slope = (elevation_difference - head_loss + residual) / total_length

        seq_df.loc[row, 'water_surface_width'] = surface_width
        seq_df.loc[row, 'velocity'] = velocity
        seq_df.loc[row, 'final_slope'] = final_slope
        seq_df.loc[row, '1inX_final_slope'] = 1 / final_slope
        seq_df.loc[row, 'total_length'] = total_length
        seq_df.loc[row, 'elevation_difference'] = elevation_difference
        seq_df.loc[row, 'head_loss_to_consider'] = head_loss
        if residual is not None:
            seq_df.loc[row, 'residual_head_loss_to_consider'] = residual
        seq_df.loc[row, 'adjusted_slope'] = adjusted_slope
        seq_df.loc[row, '1inX_adjusted_slope'] = 1 / adjusted_slope
        seq_df.loc[row, 'cost_file_name'] = cost_file
        seq_df.to_csv(f'{output_path}_seq{run_id}.csv', index=False)

        # Decide whether to continue
        if 1 < run_id < MAX_RUNS and residual < RESIDUAL_HEAD_LOSS_TOL_M:
            print(f"No need for run {run_id + 1} as the residual head loss is less than "
                  f"{RESIDUAL_HEAD_LOSS_TOL_M} m.")
            return
        if run_id == max_runs:
            return

        prev_adjusted_slope = adjusted_slope
        prev_head_loss = head_loss


# ---- Candidate channel designs ----
def candidate_designs(section, ratios, velocities, min_gap=500, manning_n=None, z=None):
    """
    Size a channel for each velocity. Returns tuples of
    (section, ratio, depth, width, 1-in-X slope, velocity), dropping designs whose
    1-in-X slope is within `min_gap` of the previous one. Uses `manning_n` and
    `z` when supplied; otherwise uses the configured `MANNING_N` and `SIDE_SLOPE`.
    """
    Q = SECTIONS[section]['gross_diversions_m3_s']
    designs = []
    for ratio in ratios:
        for V in velocities:
            slope, depth, width = solve_S(Q=Q, V=V, ratio=ratio, n=manning_n, z=z)
            designs.append((section, ratio, depth, width, int(round(1 / slope)), round(V, 3)))

        i = 0
        while i < len(designs) - 1:
            if abs(designs[i][4] - designs[i + 1][4]) < min_gap:
                designs.pop(i + 1)
            else:
                i += 1
    return designs
