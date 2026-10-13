# waterRoute

Least-cost gravity-flow canal alignment and design.

waterRoute finds the cheapest route for an open channel across a DEM between a start and a finish point, then sizes and costs the channel along that route. It also includes tools for costing pipe and electricity-transmission routes.

## Repository layout

| Path | Contents |
|---|---|
| `waterroute/` | The Python package |
| `waterroute/channelroute.py` | Least-cost path routing over a DEM (Dijkstra), cost arrays, cross drainage, control structures, route outputs |
| `waterroute/channel/` | Channel cost components: earthworks, bridges, cross drainage, control structures and siphons, seepage/evaporation losses |
| `waterroute/piperoute.py` | Pipe hydraulics (head loss, diameter, velocity), pump and hydro station costs, pipe route optimisation |
| `waterroute/electricityroute.py` | Transmission line and transformer cost model |
| `waterroute/economics.py` | Indirect costs, inflation, levelised cost |
| `waterroute/utils.py` | Raster loading and output folders |
| `tiled_flow/` | Cython flow-direction / drainage code (`drain` is imported by `channelroute.py`) |
| `channel_design.py` | Helpers used by the notebook: map display, channel sizing, routing runs with head-loss iteration |
| `elliot_waterroute.ipynb` | Main workflow notebook |
| `inputs/` | Input data (`Sequential_Run.csv` template; DEM rasters are not in the repo) |
| `tests/` | Module test notebook used durintg development |

## Installation



```
TODO
mamba env create -f waterroute/environment.yml
mamba activate channelroute
```



## Input data

The DEM is too large for GitHub (>100 MB) and isn't committed. Place it locally at:

```
inputs/DEM/dem.tif           # digital elevation map, raster
inputs/DEM/rock_depth.tif    # rock depth raster
```

`inputs/Sequential_Run.csv` is the template that records each iteration of a sequential run: slope, channel dimensions, velocity, length, head loss and adjusted slope. This is used to generate summart statistic table at the end of model run for the configuration you set.

## Usage

Open `elliot_waterroute.ipynb` and run the cells in order.

1. **Settings.** Set the paths, routing bounds, Manning's n, side slope and the per-section inputs (start/finish coordinates, start elevation, diversion flow, losses), then pass them to `channel_design.configure(...)`.
2. **Candidate designs.** `candidate_designs(section, ratios, velocities)` produces channel geometries (width:depth ratio, depth, bed width, slope) for a range of flow velocities.
3. **Routing.** For each design, `route_with_head_loss_iteration(...)` routes the channel across the DEM. It then recalculates the head loss and adjusts the slope. It repeats until the residual head loss is below `RESIDUAL_HEAD_LOSS_TOL_M` or `MAX_RUNS` is reached.

Results are written to `run1/1_Routes/<run_name>/`. Run names encode the design, for example `v22-section1-r5-d2.10-w10.50-S20000-v0.490`. A run is skipped if its folder already exists.

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
