#
# Download a digital elevation model (DEM) of Mount Rainier, Washington
#
# The DEM comes from the 3D Elevation Program (3DEP) of the U.S. Geological Survey
# (USGS), whose data are in the public domain. The 3DEP image service resamples its
# elevation mosaic to a grid of 1 arc-second (about 30 m) in geographic coordinates
# (WGS 84) and returns a GeoTIFF file. The example includes this file,
# dem/rainier.tif, so this script needs to run only to download the DEM again.
#
import requests
import rasterio
from pathlib import Path

# Longitude and latitude bounds (degrees), about 17.5 km by 17 km around the summit;
# they cover the model, its absorbing layers and a margin
lon = [-121.875, -121.645]
lat = [46.775, 46.930]
res = 1.0/3600.0

url = 'https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage'
params = {'bbox': f'{lon[0]},{lat[0]},{lon[1]},{lat[1]}',
          'bboxSR': 4326,
          'imageSR': 4326,
          'size': f'{round((lon[1] - lon[0])/res)},{round((lat[1] - lat[0])/res)}',
          'format': 'tiff',
          'pixelType': 'F32',
          'interpolation': 'RSP_BilinearInterpolation',
          'f': 'image'}
r = requests.get(url, params=params, timeout=300)
r.raise_for_status()

Path('dem').mkdir(exist_ok=True)
Path('dem/download.tif').write_bytes(r.content)

# Save the DEM as a compressed GeoTIFF file
with rasterio.open('dem/download.tif') as src:
    profile = src.profile
    profile.update(compress='deflate', predictor=3)
    h = src.read(1)
    with rasterio.open('dem/rainier.tif', 'w', **profile) as dst:
        dst.write(h, 1)
Path('dem/download.tif').unlink()

print(f'DEM: {h.shape[1]} x {h.shape[0]} points, elevation {h.min():.0f} to {h.max():.0f} m')
