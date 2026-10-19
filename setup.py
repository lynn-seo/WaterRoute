"""Install the WaterRoute package so `import waterroute` works from any folder.

Put this file next to the `waterroute` folder (not inside it), then run:

    pip install -e .

(`-e` links to your files, so edits take effect without reinstalling.)
"""
from setuptools import find_packages, setup

setup(
    name="waterroute",
    version="1.0.0",
    description="Channel, pipe and electricity routing and costing tools",
    # waterroute and waterroute.channel. The tiled_flow folder has no __init__.py,
    # so it is not packaged (it is installed separately as tiled-flow).
    packages=find_packages(include=["waterroute", "waterroute.*"]),
    python_requires=">=3.11",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "geopandas",
        "shapely",
        "pyproj",
        "rasterio",
        "matplotlib",
        "seaborn",
        "requests",
        "ipyleaflet",
        "numpy-financial",
        "fluids",
        "dijkstra3d",
        "scikit-fmm",
        "scikit-image",
        "pysheds",
    ],
    extras_require={
        # used by the notebook rather than the package itself
        "notebook": [
            "livemaptiles @ git+https://github.com/artttt/livemaptiles.git@7d9c7ac2e44399f9b071ad3fdadb6534a2ea66c3",
            "fiona",
            "ipykernel",
            "mercantile",
        ],
    },
)
