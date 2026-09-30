# FWI with model regularization

This example tests the model regularization of 2-D acoustic full-waveform inversion (FWI) on a checkerboard model. The observed seismograms contain random noise. Without regularization, the inversion fits the noise once it has recovered the checkerboard: the model error stops decreasing, and noise artifacts cover the cells. TGpV regularization removes most of these artifacts. The model error changes less, because the smooth edges of the cells, which FWI cannot sharpen at 8 Hz, dominate it. The example mirrors the LATTE example of traveltime tomography with regularization, with the same model.

## Model and data

The model is 2 km wide and 2 km deep, with a grid spacing of 20 m (101 by 101 grid points). The true model is a checkerboard of 5 by 5 cells, each 400 m wide, with P-wave velocities of 2850 and 3150 m/s. The density is constant, and the initial model has a constant velocity of 3000 m/s. Twenty sources lie on the four sides of the model, 5 on each side. They are explosions with a Ricker wavelet of 8 Hz center frequency. One hundred receivers lie on the four sides, 80 m apart. The modeling uses a time step of 2.5 ms and records 1.3 s.

We add Gaussian noise to the seismograms. The noise has the amplitude spectrum of the Ricker wavelet, so it lies in the frequency band of the data, where FWI can fit it. Its root-mean-square (RMS) amplitude is half that of each shot gather. The noise is therefore stronger than the signal of the checkerboard: the RMS difference between the data of the true and initial models is about 0.3 of the RMS amplitude of the data.

## Inversions

All four inversions run 60 iterations and differ only in the regularization:

| Parameter file | Regularization | Output directory |
|---|---|---|
| `param_fwi.rb` | None | `test_noreg` |
| `param_fwi_tgpv.rb` | TGpV from iteration 1 | `test_tgpv` |
| `param_fwi_tv.rb` | TV from iteration 1 | `test_tv` |
| `param_fwi_tgpv_delayed.rb` | TGpV from iteration 11 | `test_tgpv_delayed` |

TGpV denotes total generalized _p_-variation. It has first- and second-order terms, and their power _p_ is 0.5 by default. Without the second-order term, TGpV becomes TV (total variation), which `param_fwi_tv.rb` selects with `reg_tv_lambda2 = 0`. In each iteration, OWL denoises the updated model with the listed method, and the next iteration adds the difference between the model and its denoised version to the gradient. With `reg_scale_vp = 0.5`, this term has half the RMS amplitude of the data-misfit gradient.

In `param_fwi_tgpv_delayed.rb`, the schedule `reg_scale_vp = 1~10:0, 11:0.5` sets the strength to zero in iterations 1 to 10. These iterations have no regularization term, so they equal those of `test_noreg`. OWL still denoises the model in these iterations, so the denoised model is up to date when the strength becomes positive.

All inversions smooth the gradient with a Gaussian of 40 m and precondition it with the source energy (`yn_energy_precond = y`). They also set `jumpout_factor = 1`. By default, when the misfit stops decreasing, OWL accepts steps that increase the misfit by up to 5% to leave a local minimum. Such steps would change the model errors of some inversions and obscure the comparison. With `jumpout_factor = 1`, an inversion whose step-size search fails keeps its model.

## Running the example

`run.sh` creates the models and the geometry, computes the seismograms, adds the noise, runs the four inversions and plots the results. It needs `owl_modeling2` and `owl_fwi2` in the search path, and `numpy` and `matplotlib` for the Python scripts. It runs 20 MPI ranks, one for each source. Each inversion writes about 0.6 GB, mostly the synthetic seismograms and the adjoint sources of each iteration.

`plot.py` saves the figure as `figures/checkerboard_regularization.pdf` and `.png`. The figure shows (a) the true model, (b–e) the models after 60 iterations, and (f) the model error in each iteration. The model error is the RMS difference between the inverted and true models. `plot.py` also prints the final error, the minimum error and the iteration of the minimum for each inversion.
