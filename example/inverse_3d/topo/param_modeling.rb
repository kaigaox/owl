
nx = 201
ny = 201
nz = 101

dx = 60
dy = 60
dz = 60

dt = 6.25e-4
tmax = 10.0
data_dt = 1.0e-2

ns = 36
file_geometry = ./geometry/geometry.txt

which_medium = elastic-tti
anisotropy_type = iso
model_name = vp, vs, rho
file_vp = model/vp.bin
file_vs = model/vs.bin
file_rho = model/rho.bin

yn_free_surface = y
free_surface_dz_refine = 4
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
file_topo = model/topo.txt

dir_synthetic = data

ngroup = 36
rankx = 4
ranky = 4
rankz = 2
