import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os
from matplotlib.path import Path
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader, natural_earth
from glob import glob



# ==== User Inputs ====
input_dir = "/uufs/chpc.utah.edu/common/home/hallar-group2/climatology/stilt/dust_spl/out/2025_trajectories"

# prompt user for which dataset to plot
choice = input("Which dataset would you like to plot? (percent/aridity): ").strip().lower()

if choice == 'percent':
    nc_file = "aridity_data/percent_change_west.nc"
    var_name = 'percent_change'
    color_label = 'Percent Change (%)'
    cmap = 'RdBu'
    vmin, vmax = -100, 100
elif choice == 'aridity':
    nc_file = "aridity_data/2025_west.nc"
    var_name = 'VegDRI'
    color_label = 'Aridity Index'
    cmap = 'YlOrRd'
    vmin, vmax = 0, 2  # adjust as appropriate for your data
else:
    raise ValueError("Invalid selection. Choose 'percent' or 'aridity'.")



# ==== Load Datasets ====

# iterate through stilt trajectory outputs
for folder in os.listdir(input_dir):
    
    # define path to the outlined footprint
    outline_nc = glob(os.path.join(input_dir, folder, 'footprint_outlines', "*_outlines.nc"))
    
    # define the save name for the figure
    output_plot = f"figures/{folder}_{var_name}_with_outline.png"
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)
    
    # open the aridity or % difference map
    data_ds = xr.open_dataset(nc_file)
    if choice == 'aridity':
        data_ds = data_ds.mean(dim='time')
        
    # open the outline nc file
    outline_ds = xr.open_dataset(outline_nc[0])
    
    # pull out variables, lon, and lat values of aridity or % diff data
    data_var = data_ds[var_name]
    lon = data_ds['lon'].values
    lat = data_ds['lat'].values
    
    # generate meshgrid for spatial masking
    lon2d, lat2d = np.meshgrid(lon, lat)
    
    
    
    # ==== set up plot ====
    crs_proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Set extent based on data bounds
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    
    
    # ==== add features and border info to the map ====
    # mask out Canada and Mexico
    shapefile = shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
    reader = shpreader.Reader(shapefile)
    for record in reader.records():
        if record.attributes['NAME'] in ['Mexico', 'Canada']:
            geometry = record.geometry
            feature = ShapelyFeature([geometry], crs=crs_proj, facecolor='gray', edgecolor='none')
            ax.add_feature(feature, zorder=10)
    
    # add ocean and lakes with consistent colors
    ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#3b9b9b', zorder=1)
    ax.add_feature(cfeature.LAKES.with_scale('10m'), facecolor='#3b9b9b', zorder=1)
    
    # add border information
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=4)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.25, alpha=0.3, linestyle='--', edgecolor='black', zorder=4)
    ax.coastlines(resolution='10m', linewidth=0.25, zorder=4)
    
    # add high-quality water features on top
    ocean_shp = natural_earth(resolution='10m', category='physical', name='ocean')
    ocean_geom = Reader(ocean_shp).geometries()
    ocean_feature = ShapelyFeature(ocean_geom, crs_proj, facecolor='#3b9b9b', edgecolor='none')
    ax.add_feature(ocean_feature, linewidth=0.25, zorder=12)
    lakes_shp = natural_earth(resolution='10m', category='physical', name='lakes')
    lakes_geom = Reader(lakes_shp).geometries()
    lakes_feature = ShapelyFeature(lakes_geom, crs_proj, facecolor='#3b9b9b', edgecolor='black')
    ax.add_feature(lakes_feature, linewidth=0.25, zorder=12)
    
   
    
    # ==== plot percent difference or aridity map based on choice ====
    if choice == 'percent':
        c = ax.pcolormesh(lon, lat, data_var, cmap='RdBu', vmin=-100, vmax=100, transform=ccrs.PlateCarree(), zorder=2)
        plt.colorbar(c, ax=ax, label=f'{var_name}')
    
    elif choice == 'aridity':
        boundaries = [0, 64, 81, 97, 113, 161, 178, 193, 253, 254, 255, 256]
        colors = [
            '#800000', '#FF0000', '#FFA500', '#FFFF00', '#ADD8E6',
            '#90EE90', '#008000', '#006400', '#00008B', '#D3D3D3', '#FFFFFF'
        ]
        cmap = ListedColormap(colors)
        cmap.set_bad(color='#FFFFFF')
        norm = BoundaryNorm(boundaries, ncolors=len(colors))   
        c = ax.pcolormesh(lon, lat, data_var, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), zorder=2)
        cbar = plt.colorbar(c, ax=ax, orientation='vertical', pad=0.05, boundaries=boundaries)
        cbar.set_label(var_name)
        cbar.set_ticks(boundaries[:-3])  
        
        try:
            mask = np.ma.masked_invalid(data_var)
            outline = ax.contour(
                lon, lat, mask.mask.astype(int), levels=[0.5],
                colors='black', linewidths=0.25, transform=crs_proj, zorder=3.5
            )
        except Exception as e:
            print(f"Failed to draw vegdri outline: {e}")
        
    # ==== plot outlines ====
    outlines = outline_ds['outlines'].values  # shape: (n_outlines, n_points, 2)
    for outline in outlines:
        outline = outline[~np.isnan(outline[:, 0])]
        if len(outline) > 1:
            ax.plot(outline[:, 0], outline[:, 1], color='black', linewidth=1.0, transform=ccrs.PlateCarree(), zorder=9)
    
    
    
    # ==== plot bounding box ====
    bbox = outline_ds['bounding_box'].values  # shape: (5, 2)
    if bbox.shape[0] == 5:
        bbox_closed = np.vstack([bbox, bbox[0]])  # Close loop for plotting
        ax.plot(bbox_closed[:, 0], bbox_closed[:, 1], color='black', linewidth=0.5, linestyle='--', transform=ccrs.PlateCarree(), label='Bounding Box', zorder=9)
    

        

    
    
    # ==== Save Figure ====
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to: {output_plot}")
