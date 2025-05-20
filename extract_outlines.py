import os
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from skimage import measure
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob

input_dir = "/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/2025_trajectories"


# Main loop
for folder in os.listdir(input_dir):
    print(f"Processing event {folder}...")


    filepaths = sorted(glob(os.path.join(input_dir, folder, 'footprints', "*_foot.nc")))

    before_files = filepaths[:4]
    after_files = filepaths[-4:]
    during_files = [f for f in filepaths if f not in before_files + after_files]

    def average_footprints(files):
        datasets = []
        for f in files:
            if not os.path.exists(f):
                continue
            ds = xr.open_dataset(f)
            arr = ds[list(ds.data_vars)[0]].mean(dim='time', skipna=True)
            datasets.append(arr)
        if datasets:
            stacked = xr.concat(datasets, dim='z').mean(dim='z')
            return stacked
        return None

    during = average_footprints(during_files)
    if during is None:
        print("No valid data for 'during' period.")
        continue

    # Extract outlines from log10 footprint (ignoring zeros)
    data = during.values
    data[data == 0] = np.nan
    with np.errstate(divide='ignore'):
        log_data = np.log10(data)
    log_data = np.nan_to_num(log_data, nan=-10)

    # Threshold to define footprint region (adjust if needed)
    threshold = -8
    binary_mask = log_data > threshold

    contours = measure.find_contours(binary_mask.astype(float), 0.5)
    lon = during['lon'].values
    lat = during['lat'].values

    outlines = []
    for contour in contours:
        row_idx, col_idx = contour[:, 0], contour[:, 1]
        outline_lon = np.interp(col_idx, np.arange(len(lon)), lon)
        outline_lat = np.interp(row_idx, np.arange(len(lat)), lat)
        outlines.append(np.stack([outline_lon, outline_lat], axis=1))

    # Compute bounding box from all outlines
    all_coords = np.vstack(outlines)
    min_lon, max_lon = np.min(all_coords[:, 0]), np.max(all_coords[:, 0])
    min_lat, max_lat = np.min(all_coords[:, 1]), np.max(all_coords[:, 1])

    rectangle_coords = np.array([
        [min_lon, min_lat],
        [min_lon, max_lat],
        [max_lon, max_lat],
        [max_lon, min_lat],
        [min_lon, min_lat]  # close the box
    ])

    # Pad outlines to same length
    max_len = max(len(outline) for outline in outlines)
    padded_outlines = np.full((len(outlines), max_len, 2), np.nan)
    for j, outline in enumerate(outlines):
        padded_outlines[j, :len(outline), :] = outline
    
    # Create Dataset with two variables: 'outlines' and 'bounding_box'
    outline_ds = xr.Dataset({
        "outlines": xr.DataArray(
            padded_outlines,
            dims=["outline", "point", "coord"],
            coords={
                "outline": np.arange(len(outlines)),
                "point": np.arange(max_len),
                "coord": ["lon", "lat"]
            }
        ),
        "bounding_box": xr.DataArray(
            rectangle_coords,
            dims=["corner", "coord"],
            coords={
                "corner": np.arange(5),
                "coord": ["lon", "lat"]
            }
        )
    })

    output_path = os.path.join(input_dir, folder, f'footprint_outlines/{folder}_outlines.nc')
    outline_ds.to_netcdf(output_path)
    print(f"Saved outline and bounding box NetCDF to: {output_path}")

    # Plot the outlines and bounding box
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Add map features
    ax.set_extent([-125, -100, 30, 50], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':', linewidth=0.5)

    # Plot the bounding box
    ax.plot(rectangle_coords[:, 0], rectangle_coords[:, 1], transform=ccrs.PlateCarree(),
            color='blue', linewidth=2, linestyle='--')

    # Plot the outlines
    for outline in outlines:
        outline = np.array(outline)
        ax.plot(outline[:, 0], outline[:, 1], transform=ccrs.PlateCarree(),
                color='red', linewidth=1.5)

    # Save the plot
    outline_plot_file = output_path.replace('.nc', '.png')
    plt.savefig(outline_plot_file, dpi=300, bbox_inches='tight')
    plt.close()
