"""Run the section 2 linear-cost sensitivity sweep outside the notebook."""

import os
import sys
import argparse
from pathlib import Path

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

VERSION = "LC-SA-v2"
SECTION = "section2"
VELOCITY = 0.437
SIDE_SLOPE = 2.0
LINEAR_COST_FACTORS = (0.7, 1.3)
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


TASKS = LINEAR_COST_FACTORS


def main(task_index):
    for input_path in (SEQUENTIAL_RUN_TEMPLATE, DEM_PATH, ROCK_DEPTH_PATH):
        if not input_path.is_file():
            raise FileNotFoundError(f"Required input file not found: {input_path}")

    if not 0 <= task_index < len(TASKS):
        raise ValueError(f"Task index {task_index} is out of range 0..{len(TASKS) - 1}")

    linear_cost_factor = TASKS[task_index]
    ROUTES_DIR.mkdir(parents=True, exist_ok=True)
    channel_design.configure(
        ROUTES_DIR=str(ROUTES_DIR),
        SEQUENTIAL_RUN_TEMPLATE=str(SEQUENTIAL_RUN_TEMPLATE),
        DEM_PATH=str(DEM_PATH),
        ROCK_DEPTH_PATH=str(ROCK_DEPTH_PATH),
        ROUTING_BOUNDS=(147, -21, 149, -19),
        MANNING_N=MANNING_N,
        SIDE_SLOPE=SIDE_SLOPE,
        SAMPLE_F=250,
        MAX_FILL_SIPHON_ENDS=5,
        MAX_RUNS=MAX_RUNS,
        RESIDUAL_HEAD_LOSS_TOL_M=RESIDUAL_HEAD_LOSS_TOL_M,
        CROSS_DRAINAGE_COST_FACTOR=1,
        SECTIONS=SECTIONS,
    )

    designs = candidate_designs(
        SECTION,
        ratios=RATIOS,
        velocities=(VELOCITY,),
        z=SIDE_SLOPE,
    )
    print(
        f"Task {task_index}: section={SECTION}, velocity={VELOCITY:.3f}, "
        f"linear cost factor={linear_cost_factor}; designs={designs}",
        flush=True,
    )

    for section, ratio, depth, width, slope_1_in_x, velocity in designs:
        run_name = (
            f"{VERSION}-{section}-r{ratio}-d{depth}-w{width}"
            f"-S{slope_1_in_x}-v{velocity:.3f}"
            f"-c{int(round(linear_cost_factor * 100))}"
        )
        existing_run_dir = ROUTES_DIR / run_name
        if existing_run_dir.exists():
            print(f"Exists, skipping: {existing_run_dir}", flush=True)
            continue

        print(
            f"Running {run_name}; linear cost factor={linear_cost_factor}",
            flush=True,
        )
        route_with_head_loss_iteration(
            number_of_runs=3,
            width_to_depth_ratio=ratio,
            water_depth=depth,
            bed_width=width,
            run_name=run_name,
            section=section,
            slope_for_channel_dimension=1 / slope_1_in_x,
            earthwork_cost_weight=1,
            linear_cost_factor=linear_cost_factor,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-index",
        type=int,
        default=os.environ.get("SLURM_ARRAY_TASK_ID"),
        help="Index into the linear-cost-factor sweep.",
    )
    args = parser.parse_args()
    if args.task_index is None:
        parser.error("provide --task-index or run as a SLURM array task")
    main(int(args.task_index))
