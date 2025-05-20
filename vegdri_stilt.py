import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os
from matplotlib.path import Path

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
percent_diff_nc = "/aridity_data/percent_change_CONUS.nc"
outline_nc = "../2025_trajectories/footprint_outlines/04022025_outline.nc"
output_plot = "figures/04022025_percent_change_with_outline.png"
os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# ==== Load Datasets ====
percent_ds = xr.open_dataset(percent_diff_nc)
outline_ds = xr.open_dataset(outline_nc)

# Assume variable name is 'percent_change'
percent_change = percent_ds['percent_change']
lon = percent_ds['lon'].values
lat = percent_ds['lat'].values

# Generate meshgrid for spatial masking
lon2d, lat2d = np.meshgrid(lon, lat)

# ==== Set up Plot ====
fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Set extent based on data bounds
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

# Add geographic features
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)

# ==== Plot Percent Difference ====
c = ax.pcolormesh(lon, lat, percent_change, cmap='RdBu', vmin=-100, vmax=100, transform=ccrs.PlateCarree())
plt.colorbar(c, ax=ax, label='Percent Change (%)')

# ==== Plot Outlines ====
outlines = outline_ds['outlines'].values  # shape: (n_outlines, n_points, 2)

for outline in outlines:
    outline = outline[~np.isnan(outline[:, 0])]
    if len(outline) > 1:
        ax.plot(outline[:, 0], outline[:, 1], color='black', linewidth=1.0, transform=ccrs.PlateCarree())

# ==== Plot Bounding Box ====
bbox = outline_ds['bounding_box'].values  # shape: (5, 2)
if bbox.shape[0] == 5:
    bbox_closed = np.vstack([bbox, bbox[0]])  # Close loop for plotting
    ax.plot(bbox_closed[:, 0], bbox_closed[:, 1], color='black', linewidth=0.5, linestyle='--', transform=ccrs.PlateCarree(), label='Bounding Box')


# ==== Compute Mean and Mask ====
mean_within_region, mask_2d = compute_mean_within_regions(lon2d, lat2d, percent_change.values, outlines, bbox)
print(f"Average percent change within outlines and bounding box: {mean_within_region:.2f}%")

fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})

masked_values = np.where(mask_2d, percent_change.values, np.nan)

pc = ax.pcolormesh(lon2d, lat2d, masked_values, cmap='RdBu', vmin=-100, vmax=100, shading='auto', transform=ccrs.PlateCarree())

# Draw black outlines around the trajectory outlines
for outline in outlines:
    outline = outline[~np.isnan(outline[:, 0])]  # remove NaNs if any
    if len(outline) > 1:
        ax.plot(outline[:, 0], outline[:, 1], color='black', linewidth=1.0, transform=ccrs.PlateCarree())

# Draw bounding box outline as well (red dashed line if you want, or black)
if bbox.shape[0] == 5:
    bbox_closed = np.vstack([bbox, bbox[0]])  # close the loop
    ax.plot(bbox_closed[:, 0], bbox_closed[:, 1], color='black', linewidth=0.5, linestyle='--', transform=ccrs.PlateCarree())
    
# Set the map extent to match your data bounds
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

# Optional: add geographic features for context
ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)

ax.set_title("Values Included in Regional Mean")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(pc, ax=ax, label='Percent Change (%)')

plt.savefig("figures/percent_change_within_region.png", dpi=300, bbox_inches='tight')
plt.close()

print("Masked percent change plot saved to: figures/percent_change_within_region.png")




# ==== Save Figure ====
plt.title(f"Percent Change with STILT Footprint Outlines\nMean = {mean_within_region:.2f}%")
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved to: {output_plot}")
