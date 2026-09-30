
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
model_update = vp, vs
model_aux = rho
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_rho = model/rho.bin

min_vp = 2000
max_vp = 6500
min_vs = 1100
max_vs = 3800
step_max_vp = 150
step_max_vs = 100

yn_free_surface = y
free_surface_dz_refine = 4
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
file_topo = model/topo.txt

process_grad = smooth, mask
grad_smooth_x = 180
grad_smooth_y = 180
grad_smooth_z = 120
grad_mask = model/mask.bin

yn_energy_precond = y

dir_record = data

niter_max = 40
dir_working = test_fwi

ngroup = 36
rankx = 4
ranky = 4
rankz = 2
