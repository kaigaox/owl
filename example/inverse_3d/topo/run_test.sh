
set -e

################################################################################
# DEM: the example includes dem/rainier.tif; get_dem.py downloads it again from
# the 3D Elevation Program of the USGS

# python get_dem.py


################################################################################
# Convert the DEM into the topography, make the models and the geometry, and plot them

python dem_to_topo.py
x_runf90 create_test.f90
python plot_model.py


################################################################################
# Modeling and FWI, which need 1152 MPI ranks on a cluster (see README.md)

export OMP_NUM_THREADS=2
mpirun -np 1152 owl_modeling3 ./param_modeling.rb

export OMP_NUM_THREADS=2
mpirun -np 1152 owl_fwi3 ./param_fwi.rb
