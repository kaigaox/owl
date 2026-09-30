#!/bin/bash

set -e

export OMP_NUM_THREADS=1

# Create the models and the geometry
python3 create_model_and_geometry.py

# Compute the seismograms in the true model, then add noise to them
mpirun -np 20 owl_modeling2 param_modeling.rb
python3 add_noise.py

# FWI without regularization
mpirun -np 20 owl_fwi2 param_fwi.rb

# FWI with TGpV regularization
mpirun -np 20 owl_fwi2 param_fwi_tgpv.rb

# FWI with TV regularization
mpirun -np 20 owl_fwi2 param_fwi_tv.rb

# FWI with TGpV regularization from iteration 11
mpirun -np 20 owl_fwi2 param_fwi_tgpv_delayed.rb

# Plot the models and the model errors
python3 plot.py
