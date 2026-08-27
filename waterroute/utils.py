import rasterio
from rasterio.enums import Resampling

import datetime
import pathlib


def load_raster(file_path,bounds=None,band_indexes = None ,upscale_factor=None,resampling=Resampling.nearest):
    '''
    bounds default is None (loads the full raster)
    otherwise a tuple of spatial coordinates in this order (left, bottom, right, top)
    '''
    with rasterio.open(file_path) as src:
        if bounds is None:
            bounds = src.bounds
        #convert from bounds to a window
        win = rasterio.windows.from_bounds(*bounds, src.transform).round_offsets().round_lengths()
        
        if upscale_factor is None:
            upscale_factor = 1
        if band_indexes is None:
            band_indexes = src.indexes
        
        if type(band_indexes) is int or len(band_indexes) == 1:
            out_shape = (int(win.height * upscale_factor), int(win.width * upscale_factor))
        else:
            out_shape = (len(band_indexes), int(win.height * upscale_factor), int(win.width * upscale_factor))
        
        arr = src.read(band_indexes, window=win,out_shape=out_shape,resampling=resampling)
        
        aff = src.window_transform(win)
        # scale image transform
        aff = aff * aff.scale(
            (win.width / arr.shape[-1]),
            (win.height / arr.shape[-2])
        )

        crs = src.crs
    return (arr,aff,crs)


def out_location(base_dir,run_name):
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = pathlib.Path(base_dir)
    output_path = output_path.joinpath(run_name,stamp,run_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path