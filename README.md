# waterRoute

waterRoute identifies least-cost gravity-flow canal alignments using iterative least-cost path optimisation and hydraulic design principles.
The framework explicitly evaluates and costs earthworks, water losses, linear and hydraulic infrastructure, and uniquely in the publish literature incorporates cross-drainage infrastructure identification and sizing.
Benchmarking against an independent real-world canal design produced comparable engineering outcomes, lower costs, and substantially faster assessment.

## Model layout

| Path | Contents |
|---|---|
| `waterroute/` | |
| `waterroute/channelroute.py` | identifies least-cost path routing over a DEM/cost arrays, locate cross drainage, control structures, and calculate route outputs |
| `waterroute/channelcost.py/` | channel cost components: earthworks, bridges, cross drainage, control structures and siphons, seepage/evaporation losses |
| `waterroute/piperoute.py` | pipe hydraulics (head loss, diameter, velocity), pump and hydro station costs, pipe route optimisation |
| `waterroute/economics.py` | all the cost terms |
| `waterroute/utils.py` | raster loading and other utility functions |
| `tiled_flow/` | flow-direction / drainage code (`drain` is imported by `channelroute.py`) |
| `channel_design.py` | model configuration : map display, channel sizing, routing runs with head-loss iteration |
| `elliot_waterroute.ipynb` | main workflow notebook |
| `inputs/` | Input data (`Sequential_Run.csv` template; Note DEM rasters are not in the repo) |

## Installation
Install the packages using mamba or conda
```
mamba env create -f waterroute/environment.yml
mamba activate channelroute
```



## Input data

Place it locally at:

```
inputs/DEM/dem.tif           # digital elevation map, raster
inputs/DEM/rock_depth.tif    # rock depth raster
```

`inputs/Sequential_Run.csv` is the template that records each iteration of a sequential run: slope, channel dimensions, velocity, length, head loss and adjusted slope. This is used to generate summart statistic table at the end of model run for the configuration you set.

## Usage

The model has a flexible architecture and is modular so that users can easily adjust parameters or swap a trapezoid cross section for a different type of cross section. Open `elliot_waterroute.ipynb` and run the cells in order.

1. **Settings.** Set the paths, routing bounds, Manning's n, side slope and the per-section inputs (start/finish coordinates, start elevation, diversion flow, losses), then pass them to `channel_design.configure(...)`.
2. **Candidate designs.** `candidate_designs(section, ratios, velocities)` produces channel geometries (width:depth ratio, depth, bed width, slope) for a range of flow velocities. 
3. **Routing.** For each design, `route_with_head_loss_iteration(...)` routes the channel across the DEM. It then recalculates the head loss and adjusts the slope. It repeats until the residual head loss is below `RESIDUAL_HEAD_LOSS_TOL_M` or `MAX_RUNS` is reached.

Results are written to `run1/<run_name>/`. Run names encode the design, for example `run1-r5-d2.10-w10.50-S20000-v0.490` for width:depth ratio 5, depth 2.1m, width 10.5m, channel slope is 1:20,000, channel velocity is 0.490 m/s. 

Minimal example:

```python
import channel_design
from channel_design import candidate_designs, route_with_head_loss_iteration

channel_design.configure(ROUTES_DIR="run1/1_Routes", DEM_PATH="inputs/DEM/merged_dem.tif", ...)

for section, ratio, depth, width, slope_1_in_x, velocity in candidate_designs("section1", ratios=[5], velocities=[0.49]):
    route_with_head_loss_iteration(
        number_of_runs=3,
        width_to_depth_ratio=ratio, water_depth=depth, bed_width=width,
        run_name=f"test-{section}-v{velocity:.3f}",
        section=section,
        slope_for_channel_dimension=1 / slope_1_in_x,
    )
```
