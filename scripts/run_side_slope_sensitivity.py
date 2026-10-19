"""Run the section 2 side-slope sensitivity sweep outside the notebook."""

import os
import sys
import argparse
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import channel_design
from channel_design import candidate_designs, route_with_head_loss_iteration


PROJECT_ROOT = REPO_ROOT.parent / "elliot_main_channel"
ROUTES_DIR = Path(
    os.environ.get("WATERROUTE_ROUTES_DIR", PROJECT_ROOT / "0_run1")
).expanduser().resolve()
SEQUENTIAL_RUN_TEMPLATE = Path(
    os.environ.get("WATERROUTE_SEQUENTIAL_RUN_TEMPLATE", REPO_ROOT / "inputs/Sequential_Run.csv")
).expanduser().resolve()
DEM_PATH = Path(
    os.environ.get("WATERROUTE_DEM_PATH", REPO_ROOT / "inputs/DEM/merged_dem.tif")
).expanduser().resolve()
ROCK_DEPTH_PATH = Path(
    os.environ.get("WATERROUTE_ROCK_DEPTH_PATH", REPO_ROOT / "inputs/DEM/low_res_resampled.tif")
).expanduser().resolve()

SIDE_SLOPES = (1.5, 2.0, 2.5)
VELOCITY_RANGES = {
    "section1": np.arange(0.490, 0.491, 0.01),
    "section2": np.arange(0.425, 0.444, 0.005),
}
RATIOS = (5,)
MANNING_N = 0.03
MAX_RUNS = 4
RESIDUAL_HEAD_LOSS_TOL_M = 0.1
SECTIONS = {
    "section1": {
        "start": (147.31, -19.81),
        "finish": (147.857, -19.978),
        "start_elevation": 36,
        "gross_diversions_m3_s": 19.4,
        "default_linear_cost": 48227516,
        "approximate_length_km": 67,
    },
    "section2": {
        "start": (147.866, -19.974),
        "finish": (148.079, -20.086),
        "start_elevation": 40,
        "gross_diversions_m3_s": 10,
        "default_linear_cost": 17772388,
        "approximate_length_km": 40,
    },
}


TASKS = [
    (section, side_slope, float(velocity))
    for side_slope in SIDE_SLOPES
    for section, velocities in VELOCITY_RANGES.items()
    for velocity in velocities
]


def main(task_index):
    for input_path in (SEQUENTIAL_RUN_TEMPLATE, DEM_PATH, ROCK_DEPTH_PATH):
        if not input_path.is_file():
            raise FileNotFoundError(f"Required input file not found: {input_path}")

    if not 0 <= task_index < len(TASKS):
        raise ValueError(f"Task index {task_index} is out of range 0..{len(TASKS) - 1}")

    section, side_slope, velocity = TASKS[task_index]
    ROUTES_DIR.mkdir(parents=True, exist_ok=True)
    channel_design.configure(
        ROUTES_DIR=str(ROUTES_DIR),
        SEQUENTIAL_RUN_TEMPLATE=str(SEQUENTIAL_RUN_TEMPLATE),
        DEM_PATH=str(DEM_PATH),
        ROCK_DEPTH_PATH=str(ROCK_DEPTH_PATH),
        ROUTING_BOUNDS=(147, -21, 149, -19),
        MANNING_N=MANNING_N,
        SIDE_SLOPE=side_slope,
        SAMPLE_F=250,
        MAX_FILL_SIPHON_ENDS=5,
        MAX_RUNS=MAX_RUNS,
        RESIDUAL_HEAD_LOSS_TOL_M=RESIDUAL_HEAD_LOSS_TOL_M,
        CROSS_DRAINAGE_COST_FACTOR=1,
        SECTIONS=SECTIONS,
    )

    designs = candidate_designs(
        section,
        ratios=RATIOS,
        velocities=(velocity,),
        z=side_slope,
    )
    print(
        f"Task {task_index}: section={section}, velocity={velocity:.3f}, "
        f"side slope={side_slope}; designs={designs}",
        flush=True,
    )

    for section, ratio, depth, width, slope_1_in_x, velocity in designs:
        run_name = (
            f"SS-SA-v1-{section}-r{ratio}-d{depth}-w{width}"
            f"-S{slope_1_in_x}-v{velocity:.3f}-SS{side_slope:.1f}"
        )
        existing_run_dir = ROUTES_DIR / run_name
        if existing_run_dir.exists():
            print(f"Exists, skipping: {existing_run_dir}", flush=True)
            continue

        print(f"Running {run_name}; side slope={side_slope}", flush=True)
        route_with_head_loss_iteration(
            number_of_runs=3,
            width_to_depth_ratio=ratio,
            water_depth=depth,
            bed_width=width,
            run_name=run_name,
            section=section,
            slope_for_channel_dimension=1 / slope_1_in_x,
            earthwork_cost_weight=1,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-index",
        type=int,
        default=os.environ.get("SLURM_ARRAY_TASK_ID"),
        help="Index into the section/side-slope/velocity combinations.",
    )
    args = parser.parse_args()
    if args.task_index is None:
        parser.error("provide --task-index or run as a SLURM array task")
    main(int(args.task_index))
