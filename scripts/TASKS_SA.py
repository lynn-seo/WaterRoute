"""Run selected section 2 sensitivity scenarios outside the notebook."""

import os
import sys
import argparse
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import channel_design
from channel_design import candidate_designs, route_with_head_loss_iteration


PROJECT_ROOT = REPO_ROOT.parent / "elliot_main_channel/run11"
ROUTES_DIR = Path(
    os.environ.get("WATERROUTE_ROUTES_DIR", PROJECT_ROOT)
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

SECTION = "section2"
BASE_VELOCITY = 0.437


class VELOCITIES(tuple):
    """Tuple-like container for the channel velocities used in sensitivity runs."""

    def __new__(cls, start: float = 0.421, stop: float = 0.451, step: float = 0.002):
        values = tuple(float(v) for v in np.arange(start, stop, step))
        return super().__new__(cls, values)

    @property
    def min_value(self) -> float:
        return float(min(self)) if self else np.nan

    @property
    def max_value(self) -> float:
        return float(max(self)) if self else np.nan

    def contains(self, velocity: float, *, tol: float = 1e-12) -> bool:
        velocity = float(velocity)
        return any(np.isclose(v, velocity, atol=tol, rtol=0.0) for v in self)

    def as_list(self):
        return list(self)

    def __repr__(self):
        return f"{self.__class__.__name__}({list(self)!r})"


VELOCITIES = VELOCITIES(0.301, 0.701, 0.002)
# VELOCITIES = (BASE_VELOCITY,)
default_SIDE_SLOPE = 2
default_MANNING_N = 0.031

default_DISCOUNT_RATE = 0.07
default_EARTHWORK_COST_WEIGHT = 1
default_LINEAR_COST_WEIGHT = 1
default_CROSS_DRAINAGE_COST_WEIGHT = 1
default_HYDRAULIC_LOSS_WEIGHT = 1
default_LOSS_WEIGHT = 1
default_FREEBOARD_WEIGHT = 1


EARTHWORK_COST_WEIGHTS = (0.7, 1, 1.3)
LINEAR_COST_WEIGHTS = (0.7, 1, 1.3)
SIDE_SLOPES = (1.5, 2.0, 2.5)

CROSS_DRAINAGE_COST_WEIGHTS = (0.5, 1, 1.5)
MANNING_N_VALUES = (0.025, 0.031, 0.037)
DISCOUNT_RATES = (0.03, 0.07, 0.1)

HYDRAULIC_LOSS_FACTORS = (0.7, 1, 1.3)
LOSS_FACTORS = (0.5, 1, 1.5)
FREEBOARD_FACTORS = (0.8, 1, 1.2)

RATIOS = (5,)

MAX_RUNS = 5
RESIDUAL_HEAD_LOSS_TOL_M = 0.1  # head difference tolerance to finish iterative runs
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
        "approximate_length_km": 30,
    },
}


TASKS = (
    [
            ("EW", velocity, "earthwork_cost_weight", factor)
            for factor in EARTHWORK_COST_WEIGHTS
            for velocity in VELOCITIES
       ]
    + [
        ("LC", velocity, "linear_cost_weight", factor)
        for factor in LINEAR_COST_WEIGHTS
        for velocity in VELOCITIES
    ]
    
    + [
        ("SS", velocity, "side_slope", side_slope)
        for side_slope in SIDE_SLOPES
        for velocity in VELOCITIES
    ]
    + [
        ("CD", velocity, "cross_drainage_cost_weight", factor)
        for factor in CROSS_DRAINAGE_COST_WEIGHTS
        for velocity in VELOCITIES
    ]
    + [
        ("MN", velocity, "manning_n", manning_n) #5
        for manning_n in MANNING_N_VALUES
        for velocity in VELOCITIES
    ]
    + [
        ("DR", velocity, "discount_rate", rate)
        for rate in DISCOUNT_RATES
        for velocity in VELOCITIES
    ]
    + [
        ("HL", velocity, "hydraulic_loss_weight", factor)
        for factor in HYDRAULIC_LOSS_FACTORS
        for velocity in VELOCITIES
    ]
    + [
        ("LF", velocity, "loss_weight", factor)
        for factor in LOSS_FACTORS
        for velocity in VELOCITIES
    ]
    + [
        ("FB", velocity, "freeboard_weight", factor)
        for factor in FREEBOARD_FACTORS
        for velocity in VELOCITIES
    ]
)

TASK_GROUPS = {
    "EW": "EW",
    "LC": "LC",
    "SS": "SS",
    "CD": "CD",
    "MN": "MN",
    "DR": "DR",
    "HL": "HL", 
    "LF": "LF", 
    "FB": "FB" 
}


def main(task_index, task_group):
    task_prefix = TASK_GROUPS[task_group]
    selected_tasks = [task for task in TASKS if task[0] == task_prefix]

    for input_path in (SEQUENTIAL_RUN_TEMPLATE, DEM_PATH, ROCK_DEPTH_PATH):
        if not input_path.is_file():
            raise FileNotFoundError(f"Required input file not found: {input_path}")

    if not 0 <= task_index < len(selected_tasks):
        raise ValueError(
            f"Task index {task_index} is out of range for {task_group}; "
            f"expected 0..{len(selected_tasks) - 1}"
        )

    version, velocity, parameter, value = selected_tasks[task_index]

    # Sensitivity parameters
    earthwork_cost_weight = value if parameter == "earthwork_cost_weight" else default_EARTHWORK_COST_WEIGHT
    linear_cost_weight = value if parameter == "linear_cost_weight" else default_LINEAR_COST_WEIGHT

    side_slope = value if parameter == "side_slope" else default_SIDE_SLOPE
    cross_drainage_cost_weight = value if parameter == "cross_drainage_cost_weight" else default_CROSS_DRAINAGE_COST_WEIGHT

    manning_n = value if parameter == "manning_n" else default_MANNING_N
    discount_rate = value if parameter == "discount_rate" else default_DISCOUNT_RATE

    hydraulic_loss_weight = value if parameter == "hydraulic_loss_weight" else default_HYDRAULIC_LOSS_WEIGHT
    loss_weight = value if parameter == "loss_weight" else default_LOSS_WEIGHT

    freeboard_weight = value if parameter == "freeboard_weight" else default_FREEBOARD_WEIGHT


    task_routes_dir = ROUTES_DIR / task_group
    task_routes_dir.mkdir(parents=True, exist_ok=True)
    #this defines modelling configurations
    channel_design.configure(
        ROUTES_DIR=str(task_routes_dir),
        SEQUENTIAL_RUN_TEMPLATE=str(SEQUENTIAL_RUN_TEMPLATE),
        DEM_PATH=str(DEM_PATH),
        ROCK_DEPTH_PATH=str(ROCK_DEPTH_PATH),
        ROUTING_BOUNDS=(147, -21, 149, -19),
        SAMPLE_F=250,
        MAX_FILL_SIPHON_ENDS=5,
        MAX_RUNS=MAX_RUNS,
        RESIDUAL_HEAD_LOSS_TOL_M=RESIDUAL_HEAD_LOSS_TOL_M,

        EARTHWORK_COST_WEIGHT=earthwork_cost_weight,
        LINEAR_COST_WEIGHT=linear_cost_weight,
        SIDE_SLOPE=side_slope,
        CROSS_DRAINAGE_COST_WEIGHT=cross_drainage_cost_weight,

        MANNING_N=manning_n,
        DISCOUNT_RATE=discount_rate,

        HYDRAULIC_LOSS_WEIGHT=hydraulic_loss_weight,
        LOSS_WEIGHT=loss_weight,
        FREEBOARD_WEIGHT=freeboard_weight,

        SECTIONS=SECTIONS,
    )
    #this desides intiial channel dimension, slope, velocity
    designs = candidate_designs(
        SECTION,
        ratios=RATIOS,
        velocities=(velocity,),
        manning_n=manning_n,
        z=side_slope,
    )
    print(
        f"Task {task_index}: {version}, section={SECTION}, "
        f"velocity={velocity:.3f}, {parameter}={value}; designs={designs}",
        flush=True,
    )

    for section, ratio, depth, width, slope_1_in_x, velocity in designs:
        if parameter == "linear_cost_weight":
            suffix = f"-LC{float(round(value, 3))}"
        elif parameter == "earthwork_cost_weight":
            suffix = f"-EW{float(round(value, 3))}"
        elif parameter == "side_slope":
            suffix = f"-SS{float(round(value, 3))}"
        elif parameter == "cross_drainage_cost_weight":
            suffix = f"-CD{float(round(value, 3))}"
        elif parameter == "manning_n":
            suffix = f"-MN{float(round(value, 3))}"
        elif parameter == "discount_rate":
            suffix = f"-DR{float(round(value, 3))}"
        elif parameter == "hydraulic_loss_weight":
            suffix = f"-HL{float(round(value, 3))}"
        elif parameter == "loss_weight":
            suffix = f"-LF{float(round(value, 3))}"
        elif parameter == "freeboard_weight":
            suffix = f"-FB{float(round(value, 3))}"
        

        run_name = (
            f"{version}-{section}-r{ratio}-d{depth}-w{width}"
            f"-S{slope_1_in_x}-v{velocity:.3f}{suffix}"
        )
        existing_run_dir = task_routes_dir / run_name
        if existing_run_dir.exists():
            print(f"Exists, skipping: {existing_run_dir}", flush=True)
            continue

        print(f"Running {run_name}; {parameter}={value}", flush=True)
        route_kwargs = {
            "earthwork_cost_weight": earthwork_cost_weight,
            "linear_cost_factor": linear_cost_weight,
        }
        route_with_head_loss_iteration(
            number_of_runs=3,
            width_to_depth_ratio=ratio,
            water_depth=depth,
            bed_width=width,
            run_name=run_name,
            section=section,
            slope_for_channel_dimension=1 / slope_1_in_x,
            **route_kwargs,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-index",
        type=int,
        default=os.environ.get("SLURM_ARRAY_TASK_ID"),
        help="Index within the selected sensitivity task group.",
    )
    parser.add_argument(
        "--task-group",
        required=True,
        choices=TASK_GROUPS,
        help="Sensitivity task group to run.",
    )
    args = parser.parse_args()
    if args.task_index is None:
        parser.error("provide --task-index or run as a SLURM array task")
    main(int(args.task_index), args.task_group)
