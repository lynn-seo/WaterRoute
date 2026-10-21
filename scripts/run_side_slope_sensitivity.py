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


PROJECT_ROOT = REPO_ROOT.parent / "elliot_main_channel"
ROUTES_DIR = Path(
    os.environ.get("WATERROUTE_ROUTES_DIR", PROJECT_ROOT / "0_run2")
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


VELOCITIES = VELOCITIES(0.421, 0.451, 0.002)
default_SIDE_SLOPE = 2
default_MANNING_N = 0.03
default_DISCOUNT_RATE = 0.07

SIDE_SLOPES = (1.5, 2.0, 2.5)
LINEAR_COST_FACTORS = (0.7, 1, 1.3)
EARTHWORK_COST_WEIGHTS = (0.7, 1, 1.3)
CROSS_DRAINAGE_COST_FACTORS = (0.5, 1, 1.5)
MANNING_N_VALUES = (0.022, 0.031, 0.036)
DISCOUNT_RATES = (0.03, 0.07, 0.1)
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
        ("LC-SA-v2", velocity, "linear_cost_factor", factor)
        for factor in LINEAR_COST_FACTORS
        for velocity in VELOCITIES
    ]
    + [
        ("earthwork-SA-v1", velocity, "earthwork_cost_weight", weight)
        for weight in EARTHWORK_COST_WEIGHTS
        for velocity in VELOCITIES
    ]
    + [
        ("SS-SA-v1", velocity, "side_slope", side_slope)
        for side_slope in SIDE_SLOPES
        for velocity in VELOCITIES
    ]
    + [
        ("CD-SA-v1", velocity, "cross_drainage_cost_factor", factor)
        for factor in CROSS_DRAINAGE_COST_FACTORS
        for velocity in VELOCITIES
    ]
    + [
        ("MN-SA-v1", velocity, "manning_n", manning_n)
        for manning_n in MANNING_N_VALUES
        for velocity in VELOCITIES
    ]
    + [
        ("DR-SA-v1", velocity, "discount_rate", rate)
        for rate in DISCOUNT_RATES
        for velocity in VELOCITIES
    ]
)


def main(task_index):
    for input_path in (SEQUENTIAL_RUN_TEMPLATE, DEM_PATH, ROCK_DEPTH_PATH):
        if not input_path.is_file():
            raise FileNotFoundError(f"Required input file not found: {input_path}")

    if not 0 <= task_index < len(TASKS):
        raise ValueError(f"Task index {task_index} is out of range 0..{len(TASKS) - 1}")

    version, velocity, parameter, value = TASKS[task_index]
    side_slope = value if parameter == "side_slope" else default_SIDE_SLOPE
    manning_n = value if parameter == "manning_n" else default_MANNING_N
    cross_drainage_cost_factor = (
        value if parameter == "cross_drainage_cost_factor" else 1
    )
    discount_rate = value if parameter == "discount_rate" else default_DISCOUNT_RATE
    ROUTES_DIR.mkdir(parents=True, exist_ok=True)
    channel_design.configure(
        ROUTES_DIR=str(ROUTES_DIR),
        SEQUENTIAL_RUN_TEMPLATE=str(SEQUENTIAL_RUN_TEMPLATE),
        DEM_PATH=str(DEM_PATH),
        ROCK_DEPTH_PATH=str(ROCK_DEPTH_PATH),
        ROUTING_BOUNDS=(147, -21, 149, -19),
        MANNING_N=manning_n,
        SIDE_SLOPE=side_slope,
        SAMPLE_F=250,
        MAX_FILL_SIPHON_ENDS=5,
        MAX_RUNS=MAX_RUNS,
        RESIDUAL_HEAD_LOSS_TOL_M=RESIDUAL_HEAD_LOSS_TOL_M,
        CROSS_DRAINAGE_COST_FACTOR=cross_drainage_cost_factor,
        DISCOUNT_RATE=discount_rate,
        SECTIONS=SECTIONS,
    )

    designs = candidate_designs(
        SECTION,
        ratios=RATIOS,
        velocities=(velocity,),
        manning_n=manning_n if parameter == "manning_n" else None,
        z=side_slope if parameter == "side_slope" else None,
    )
    print(
        f"Task {task_index}: {version}, section={SECTION}, "
        f"velocity={velocity:.3f}, {parameter}={value}; designs={designs}",
        flush=True,
    )

    for section, ratio, depth, width, slope_1_in_x, velocity in designs:
        if parameter in ("linear_cost_factor", "earthwork_cost_weight"):
            suffix = (
                f"-{'c' if parameter == 'linear_cost_factor' else 'EW'}"
                f"{float(round(value, 3))}"
            )
        elif parameter == "side_slope":
            suffix = f"-SS{float(round(value, 3))}"
        elif parameter == "cross_drainage_cost_factor":
            suffix = f"-CD{float(round(value, 3))}"
        elif parameter == "manning_n":
            suffix = f"-ME{float(round(value, 3))}"
        else:
            suffix = f"-DR{value:.3f}"
        run_name = (
            f"{version}-{section}-r{ratio}-d{depth}-w{width}"
            f"-S{slope_1_in_x}-v{velocity:.3f}{suffix}"
        )
        existing_run_dir = ROUTES_DIR / run_name
        if existing_run_dir.exists():
            print(f"Exists, skipping: {existing_run_dir}", flush=True)
            continue

        print(f"Running {run_name}; {parameter}={value}", flush=True)
        route_kwargs = {
            "earthwork_cost_weight": (
                value if parameter == "earthwork_cost_weight" else 1
            ),
        }
        if parameter == "linear_cost_factor":
            route_kwargs["linear_cost_factor"] = value
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
        help="Index into the combined notebook sensitivity scenarios.",
    )
    args = parser.parse_args()
    if args.task_index is None:
        parser.error("provide --task-index or run as a SLURM array task")
    main(int(args.task_index))
