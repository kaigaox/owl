# 3-D elastic FWI of Mount Rainier

This example sets up 3-D elastic full-waveform inversion (FWI) with the topography of Mount Rainier, Washington. It goes through the entire process: a digital elevation model (DEM), the topography of the model, the models, the modeling of the observed data and the FWI. The subsurface model is artificial. It follows the main geologic units of Mount Rainier, but their shapes and velocities are not taken from a published model. The FWI updates the P- and S-wave velocities from smooth initial models.

## DEM

`dem/rainier.tif` is a GeoTIFF file of the DEM. It comes from the 3D Elevation Program (3DEP) of the U.S. Geological Survey, whose data are in the public domain. The DEM covers longitudes from 121.875°W to 121.645°W and latitudes from 46.775°N to 46.930°N on a grid of 1 arc-second (about 30 m) in geographic coordinates (WGS 84), with 828 by 558 points. Its elevation ranges from 818 to 4391 m above sea level. `get_dem.py` downloads the DEM again from the 3DEP image service; another source of DEMs works as well if the file covers the model.

## Topography

`dem_to_topo.py` converts the DEM into the topography of the model with rasterio, a Python interface to GDAL. It reprojects the DEM to UTM zone 10N (EPSG:32610) on a 20 m grid. It then smooths the DEM with a Gaussian of 120 m, twice the grid spacing of the model, to avoid aliasing. The S waves of the simulation have wavelengths of about 1 km, so they do not sense the removed relief. Last, it samples the DEM at the grid points of the model and of its absorbing layers.

The model grid has 201 by 201 points with a spacing of 60 m along x (east) and y (north). Its origin lies at 588,500 m E and 5,183,500 m N, so the model covers 12 km by 12 km around the summit. The elevation of the topography ranges from 1388 to 4365 m, with the summit at x = 5.94 km and y = 6.00 km. The slope has a mean of 23°, a 95th percentile of 38° and a maximum of 51°.

`model/topo.bin` stores the elevation above sea level on the model grid. `model/topo.txt` lists x, y and elevation on the grid extended by the absorbing layers, as OWL requires. In this file, the elevation is measured from the lowest point, 1241 m above sea level. OWL maps each column between the surface and the bottom of the model to a mesh whose height is at most the depth of the bottom below the zero elevation. With the elevation measured from the lowest point, this depth is 2.9 km instead of 1.6 km, and the mesh beneath the summit is stretched less.

## Model

The model extends 6 km downward from the summit to an elevation of −1.6 km, with 101 grid points in depth. `create_test.f90` builds it with libflit and librgm.

The volcanic edifice is a stack of 20 layers of andesite lava flows and breccias. The layers follow the surface smoothed over 600 m, so they are parallel to the flanks of the cone. Their P-wave velocity increases with depth from 3000 to 4200 m/s. The edifice overlies an unconformity whose elevation is 1900 m plus 15% of the difference between the topography, smoothed over 2 km, and 1900 m. The unconformity thus lies between 1.87 and 2.11 km, and the basement crops out in the deep valleys, over about a fifth of the surface. The basement consists of 12 layers of Tertiary volcanic rocks, cut by 4 faults whose slip varies along the faults. Its P-wave velocity increases with depth from 4200 to 5600 m/s. The ratio of the P- and S-wave velocities varies from layer to layer, between 1.72 and 1.85 in the edifice and between 1.70 and 1.80 in the basement. The density follows the Gardner relation, ρ = 310·vp^0.25, with ρ in kg/m³ and vp in m/s.

Weathering lowers the velocities near the surface. At the surface, the P- and S-wave velocities are 80% and 76% of those of the fresh rock, and they increase linearly to the base of the weathered zone. The zone is 60 m thick on the ridges and up to 180 m thick in the valleys, where the topography is concave. The grid points above the surface take the values just below the surface, and `model/mask.bin` marks them with zero.

Two bodies are the targets of the FWI. A granodiorite pluton, like the Tatoosh pluton, lies beneath the south flank. It is centered 3 km south of the summit, with semi-axes of 4.5 km east-west and 3.5 km north-south and irregular flanks. Its domed top lies at an elevation of 1300 m and descends below the model toward the flanks. The pluton has a P-wave velocity of 5900 m/s, an S-wave velocity of 3400 m/s and a density of 2700 kg/m³. Hydrothermal alteration weakens the upper edifice within 900 m of the axis of the summit, from 150 m below the surface down to an elevation of 2800 m. There, the P- and S-wave velocities decrease by up to 15% and 25%, with smooth edges. In the whole model, the P-wave velocity ranges from 2400 to 5900 m/s, and the S-wave velocity from 1335 to 3400 m/s.

The initial models come from the model before the pluton and the alteration are inserted. `create_test.f90` smooths their slownesses in coordinates that follow the surface, with Gaussian widths of 180 m in depth and 1.2 km horizontally. The initial models thus keep a smoothed weathered zone and the smoothed contrast between the edifice and the basement, but they contain neither target. The density stays fixed during the inversion.

## Acquisition and modeling

Thirty-six explosive sources lie 25 m below the surface on a 2 km grid, from 1 to 11 km in x and y. They have a Ricker wavelet with a center frequency of 1.5 Hz. Three-component receivers record on the surface on a 300 m grid, from 0.3 to 11.7 km in x and y (1521 receivers).

The modeling uses the elastic solver with a free surface that follows the topography. OWL refines the vertical grid near the surface by `free_surface_dz_refine = 4`, which increases the number of vertical grid points from 101 to 173 for the true model and 183 for the initial models. With this refinement, the frequency that OWL considers free of dispersion is 1.62 Hz for the true model, above the 1.5 Hz center frequency. The refined cells on the steepest slopes limit the stable time step to 0.65 ms. The parameter files set 0.625 ms, so each simulation of 10 s takes 16,001 time steps. The records of 10 s include the surface waves at the largest offsets of about 15 km. The seismograms are saved every 10 ms.

## Inversion

`param_fwi.rb` inverts the waveforms for the P- and S-wave velocities in 40 iterations. The gradients are smoothed with Gaussian widths of 180 m horizontally and 120 m in depth, masked above the surface, and preconditioned with the source energy. The P- and S-wave velocities stay within 2000–6500 m/s and 1100–3800 m/s, and each update changes them by at most 150 and 100 m/s.

## Running the example

`run_test.sh` converts the DEM, builds and plots the models, computes the observed seismograms and runs the FWI. `get_dem.py` needs `requests` and `rasterio`, `dem_to_topo.py` needs `rasterio` and `scipy`, `x_runf90` compiles `create_test.f90` with libflit and librgm, and `plot_model.py` needs `matplotlib`. The `rasterio` wheels include GDAL, so `pip install rasterio` suffices. The example is meant for a cluster. The parameter files assign each source to one of 36 groups, each of which decomposes the model into 4 by 4 by 2 subdomains, so they need 1152 MPI ranks. For fewer ranks, reduce `ngroup` (the groups then process several sources in turn) or `rankx`, `ranky` and `rankz`.

`plot_model.py` saves the figure as `figures/model.pdf` and `.png`, and `dem/mt_rainer_artificial_model.png` is a copy of it. The figure shows (a) the topography with the sources (stars) and receivers (dots), (b) a horizontal slice of the true P-wave velocity at an elevation of 1 km, and (c–f) north-south slices through the summit of the true and initial P- and S-wave velocities.
