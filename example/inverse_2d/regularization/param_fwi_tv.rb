
nx = 101
nz = 101
dx = 20
dz = 20

dt = 2.5e-3
tmax = 1.3

ns = 20
file_geometry = geometry/geometry.txt

model_update = vp
file_vp = model/vp_init.bin

dir_record = data_noisy

process_grad = smooth
grad_smooth_x = 40
grad_smooth_z = 40

yn_energy_precond = y

niter_max = 60
step_max_vp = 50
min_vp = 2000
max_vp = 4000
jumpout_factor = 1.0

model_regularization_method = tgpv
reg_tv_lambda2 = 0
reg_scale_vp = 0.5

dir_working = test_tv
