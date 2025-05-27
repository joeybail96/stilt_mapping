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

# ==== Function to Compute Mean Within Outlines and BBox ====
def compute_mean_within_regions(lon2d, lat2d, data2d, outlines, bbox):
    # Assume the first outline is the primary trajectory
    trajectory = outlines[0]
    trajectory = trajectory[~np.isnan(trajectory[:, 0])]

    if len(trajectory) < 2 or bbox.shape[0] < 4:
        raise ValueError("Invalid trajectory or bounding box for constructing closed polygon.")

    # Construct a closed polygon by appending the bounding box to trajectory
    # Find the nearest points on the bbox to the start and end of trajectory
    start_pt = trajectory[0]
    end_pt = trajectory[-1]

    # Remove potential duplicate closing point in bbox
    bbox_unique = bbox[:-1] if np.allclose(bbox[0], bbox[-1]) else bbox

    # Compute distances from bbox points to start and end
    dists_to_start = np.linalg.norm(bbox_unique - start_pt, axis=1)
    dists_to_end = np.linalg.norm(bbox_unique - end_pt, axis=1)

    start_idx = np.argmin(dists_to_start)
    end_idx = np.argmin(dists_to_end)

    # Get the bbox segment connecting end → start (in bbox order)
    if start_idx < end_idx:
        bbox_segment = bbox_unique[start_idx:end_idx + 1][::-1]  # reverse for CCW
    else:
        bbox_segment = np.concatenate([bbox_unique[start_idx:], bbox_unique[:end_idx + 1]])[::-1]

    # Combine trajectory + bbox segment to form closed polygon
    closed_polygon = np.vstack([trajectory, bbox_segment, trajectory[0]])

    # Masking
    points = np.vstack((lon2d.ravel(), lat2d.ravel())).T
    path = Path(closed_polygon)
    mask = path.contains_points(points).reshape(data2d.shape)
    mean_value = np.nanmean(np.where(mask, data2d, np.nan))
    return mean_value, mask



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
    ax.add_feature(cfeature.BORDERS, linewidth=1, zorder=4)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.5, edgecolor='black', zorder=4)
    ax.coastlines(resolution='10m', linewidth=0.5, zorder=4)
    
    # add high-quality water features on top
    ocean_shp = natural_earth(resolution='10m', category='physical', name='ocean')
    ocean_geom = Reader(ocean_shp).geometries()
    ocean_feature = ShapelyFeature(ocean_geom, crs_proj, facecolor='#3b9b9b', edgecolor='none')
    ax.add_feature(ocean_feature, zorder=12)
    lakes_shp = natural_earth(resolution='10m', category='physical', name='lakes')
    lakes_geom = Reader(lakes_shp).geometries()
    lakes_feature = ShapelyFeature(lakes_geom, crs_proj, facecolor='#3b9b9b', edgecolor='black')
    ax.add_feature(lakes_feature, zorder=12)
    
   
    
    # ==== plot percent difference or aridity map based on choice ====
    if choice == 'percent':
        c = ax.pcolormesh(lon, lat, data_var, cmap='RdBu', vmin=-100, vmax=100, transform=ccrs.PlateCarree(), zorder=9)
        plt.colorbar(c, ax=ax, label=f'{var_name}')
    
    elif choice == 'aridity':
        boundaries = [0, 64, 81, 97, 113, 161, 178, 193, 253, 254, 255, 256]
        colors = [
            '#800000', '#FF0000', '#FFA500', '#FFFF00', '#FFFFFF',
            '#90EE90', '#008000', '#006400', '#00008B', '#D3D3D3', '#FFFFFF'
        ]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(boundaries, ncolors=len(colors))   
        c = ax.pcolormesh(lon, lat, data_var, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), zorder=9)
        cbar = plt.colorbar(c, ax=ax, orientation='vertical', pad=0.05, boundaries=boundaries)
        cbar.set_label(var_name)
        cbar.set_ticks(boundaries[:-3])  
        
        
        
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
    
    
    # # add high-quality water features on top
    # ocean_shp = natural_earth(resolution='10m', category='physical', name='ocean')
    # ocean_geom = Reader(ocean_shp).geometries()
    # ocean_feature = ShapelyFeature(ocean_geom, crs_proj, facecolor='#3b9b9b', edgecolor='none')
    # ax.add_feature(ocean_feature, zorder=12)
    
    # lakes_shp = natural_earth(resolution='10m', category='physical', name='lakes')
    # lakes_geom = Reader(lakes_shp).geometries()
    # lakes_feature = ShapelyFeature(lakes_geom, crs_proj, facecolor='#3b9b9b', edgecolor='black')
    # ax.add_feature(lakes_feature, zorder=12)
    
    
    
    # ==== compute mean and mask ====
    # this code doesn't always work. Might need to scorch this and start over
    mean_within_region, mask_2d = compute_mean_within_regions(lon2d, lat2d, data_var.values, outlines, bbox)
    print(f"Average {var_name} within outlines and bounding box: {mean_within_region:.2f}")
    
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    masked_values = np.where(mask_2d, data_var.values, np.nan)
    
    if choice == 'percent':
        pc = ax.pcolormesh(lon2d, lat2d, masked_values, cmap='RdBu', vmin=-100, vmax=100, shading='auto', transform=ccrs.PlateCarree(), zorder=9)
        plt.colorbar(pc, ax=ax, label=f'{var_name}')
        
    elif choice == 'aridity':
        boundaries = [0, 64, 81, 97, 113, 161, 178, 193, 253, 254, 255, 256]
        colors = [
            '#800000', '#FF0000', '#FFA500', '#FFFF00', '#FFFFFF',
            '#90EE90', '#008000', '#006400', '#00008B', '#D3D3D3', '#FFFFFF'
        ]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(boundaries, ncolors=len(colors))
        
        pc = ax.pcolormesh(lon2d, lat2d, masked_values, cmap=cmap, norm=norm, shading='auto', transform=ccrs.PlateCarree(), zorder=9)
        cbar = plt.colorbar(pc, ax=ax, orientation='vertical', pad=0.05, boundaries=boundaries)
        cbar.set_label(var_name)
        cbar.set_ticks(boundaries[:-3])  
        
    # Draw black outlines around the trajectory outlines
    for outline in outlines:
        outline = outline[~np.isnan(outline[:, 0])]  # remove NaNs if any
        if len(outline) > 1:
            ax.plot(outline[:, 0], outline[:, 1], color='black', linewidth=1.0, transform=ccrs.PlateCarree(), zorder=9)
    
    # Draw bounding box outline as well (red dashed line if you want, or black)
    if bbox.shape[0] == 5:
        bbox_closed = np.vstack([bbox, bbox[0]])  # close the loop
        ax.plot(bbox_closed[:, 0], bbox_closed[:, 1], color='black', linewidth=0.5, linestyle='--', transform=ccrs.PlateCarree(), zorder=9)
        
    # Set the map extent to match your data bounds
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Mask out Canada and Mexico
    shapefile = shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
    reader = shpreader.Reader(shapefile)
    for record in reader.records():
        if record.attributes['NAME'] in ['Mexico', 'Canada']:
            geometry = record.geometry
            feature = ShapelyFeature([geometry], crs=crs_proj, facecolor='gray', edgecolor='none')
            ax.add_feature(feature, zorder=10)
    
    # Optional: add geographic features for context
    # Add ocean and lakes with consistent colors
    ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#3b9b9b', zorder=1)
    ax.add_feature(cfeature.LAKES.with_scale('10m'), facecolor='#3b9b9b', zorder=1)
    
    ax.add_feature(cfeature.BORDERS, linewidth=1, zorder=4)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.5, edgecolor='black', zorder=4)
    ax.coastlines(resolution='10m', linewidth=0.5, zorder=4)
    
    # Add high-quality water features on top
    ocean_shp = natural_earth(resolution='10m', category='physical', name='ocean')
    ocean_geom = Reader(ocean_shp).geometries()
    ocean_feature = ShapelyFeature(ocean_geom, crs_proj, facecolor='#3b9b9b', edgecolor='none')
    ax.add_feature(ocean_feature, zorder=12)
    
    lakes_shp = natural_earth(resolution='10m', category='physical', name='lakes')
    lakes_geom = Reader(lakes_shp).geometries()
    lakes_feature = ShapelyFeature(lakes_geom, crs_proj, facecolor='#3b9b9b', edgecolor='black')
    ax.add_feature(lakes_feature, zorder=12)
    
    ax.set_title("Values Included in Regional Mean")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    
    
    masked_plot = output_plot.replace('.png', '_masked.png')
    plt.savefig(masked_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Masked {var_name} plot saved to: figures/{var_name}_with_outline_masked.png")
    
    
    
    
    # ==== Save Figure ====
    plt.title(f"Mean = {mean_within_region:.2f}")
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to: {output_plot}")
