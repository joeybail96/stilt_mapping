import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
from glob import glob


input_dir = "/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/2025_trajectories"
# Output directory
output_dir = "/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/scripts/figures/raw_footprints"


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

    before = average_footprints(before_files)
    during = average_footprints(during_files)
    after = average_footprints(after_files)
    
    during_save = os.path.join(input_dir, folder, 'footprints',"during.nc")
    during.to_netcdf(during_save)
    

    def plot_footprint(ax, data, title, vmax):
        if data is None:
            return
        lon = data['lon'].values
        lat = data['lat'].values
        val = data.values
        val[val == 0] = np.nan
        with np.errstate(divide='ignore'):
            log_val = np.log10(val)
        alpha = np.where(np.isnan(log_val), 0.0, 0.7)

        mesh = ax.pcolormesh(lon, lat, log_val, cmap='rainbow', vmin=-8, vmax=vmax, transform=ccrs.PlateCarree())
        ax.set_extent([-125, -100, 30, 50], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)
        ax.set_title(title)
        return mesh

    # Determine colorbar limits
    vmax = np.nanmax([np.nanmax(np.log10(d.values)) for d in [before, during, after] if d is not None])

    # Plot
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), subplot_kw={'projection': ccrs.PlateCarree()})
    mesh_b = plot_footprint(axs[0], before, "Before", vmax)
    mesh_d = plot_footprint(axs[1], during, "During", vmax)
    mesh_a = plot_footprint(axs[2], after, "After", vmax)

    # Shared colorbar
    cbar = fig.colorbar(mesh_d, ax=axs.ravel().tolist(), shrink=0.7, orientation='horizontal', pad=0.05)
    cbar.set_label("log10(Footprint)")

    # Save plot
    output_file = os.path.join(output_dir, folder)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
