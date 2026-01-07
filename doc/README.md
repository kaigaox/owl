
# Short Tutorial for OWL

This is short tutorial on how to install and use `OWL` for forward wavefield modeling and FWI.

## Table of Contents
- [Installation](#installation)
- [Parameters](#parameters)
- [Examples](#examples)


## Installation

Installation of `OWL` is straightforward:

```bash
git clone https://github.com/lanl/owl.git
cd owl
ruby install.rb clean
```

You need to install [`FLIT`](https://github.com/lanl/flit) before installing `OWL`. 


## Geometry

To use `OWL` for forward modeling or FWI, a source-receiver geometry file must be provided. The format of `OWL`'s geometry file is:

### Overall geometry file

`OWL` requires an overall geometry file (say, `./geometry/geometry.txt`) in the following form:
 
```ruby
	shot_1_geometry.txt
	shot_2_geometry.txt
	shot_3_geometry.txt
	....
```

where each row is the name of each common-shot gather's geometry file. The names of the individual geometry files can be arbitrary, e.g., 

```ruby
	shot_1_geometry.txt
	source_222_geometry.txt
	event_33_sr.txt
	...
```

However, these geometry files should be distinct (otherwise, it does not make much sense to include two identical common-shot gathers in one survey) and should present in the same directory with the overall geometry file, i.e., the path to individual geometry files is `./geometry/shot_1_geometry.txt`, ...

### Individual geometry file

An individual geometry file is to describe the source-receiver distribution as well as source parameters. The form is:

```ruby
shot id

number of point sources
point source 1 location
point source 1 mechanism
point source 1 time function, frequency, amplitude, origin time
point source 1 time function processing (frequency filtering)
point source 2 location
point source 2 mechanism
point source 2 time function, frequency, amplitude, origin time
point source 2 time function processing (frequency filtering)
...

number of receivers
receiver 1 location, weight
receiver 2 location, weight
...
```

For example:
```ruby
1								# The unique id of this shot is 1

1								# The shot contains 1 point source
100.0 500.0 10.0				# The x, y, z location of this point source is (100, 500, 10) meters
force 45.0 45.0					# The source is a force vector with (polar, azimuth) = (45, 45) degrees
ricker 20.0 1.0e4 0.0			# The source time function is a Ricker with f0 = 20, A0 = 1e4, and t0 = 0
0.0 0.0							# No processing (frequency filtering) is applied to this source

2								# The shot contains 2 receivers
200.0 400.0 10.0 1.0			# The 1st receiver locates at (200, 400, 10) meters, with a weight of 1
1200.0 1400.0 100.0 1.0			# The 2nd receiver locates at (1200, 1400, 100) meters, with a weight of 1
```

The source can be an explosion source (i.e., isotropic moment tensor): 
```ruby
	...
	explosion					# The source is an explosion source
	...
```

Or a general moment tensor source:
```ruby
	...
	mt 1.0 1.0 -1.0 -0.5 -0.2 0.1		# The notation convention is (Mxx, Myy, Mzz, Mxy, Mxz, Myz) 
	...
```

For source time function, the valid options are:
- `gaussian`
- `gaussian_deriv` (the 1st-order derivative of Gaussian)
- `gaussian_deriv_deriv`, `ricker` (the 2nd-order derivative of Gaussian)
- `gaussian_deriv_deriv_deriv`, `ricker_deriv` (the 3rd-order derivative of Gaussian)
- `ormsby` (approximation to sinc with 4 corner frqeuencies)
- `custom` (user-provided custom source time function)

To use custom stf, the user must use the following form:
```ruby
	...
	custom 10.0 1.0 0.0			# The f0 here is to accurately using ADE-CFS-MPML
	custom_wavelet.txt			# Name of the custom stf file, in ASCII format
	... 
```

The custom stf file must be in the following format:
```ruby
	t1 amp1
	t2 amp2
	t3 amp3
	...
	tn ampn
```

that is, each row contains `time amplitude` of the stf, where the time is in second. For example:

```ruby
	0.0			0.0
	1.0e-3		2.0e-5
	2.0e-3		2.5e-5
	....
	1.0e-1		1.0e0
	...
	2.0e-1		0.0
```
The `time` can be irregularly sampled, although in practice they are usually regularly sampled. The provided stf will be resampled during modeling or FWI to be consistent with the solver. Therefore, the sampling interval can be arbitrary, as long as it is consistent with physics of the survey. 

For FWI applications, due to historical limitations of SU header, here we don't use SU to store source-receiver geometry infomration. Therefore, in the provided observed data, the SU files can have null headers (except time-related headers like sampling interval and number of samples), but the number of receivers must be consistent with that in the geometry file. Future versions of `OWL` may seek using other more flexible data format for I/O.  

## Model


## Dimension

## Other parameters

<!-- > **`n1` (integer)** 
- **Description**: Number of grid points along axis-1 of the generated random model.
- **Default**: `128` -->