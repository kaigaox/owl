#
# Plot the topography, the acquisition geometry, and slices of the true and initial models
#
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LightSource, LinearSegmentedColormap, Normalize
from pathlib import Path

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 11, 'axes.titlesize': 11, 'axes.labelsize': 11,
                     'xtick.labelsize': 10, 'ytick.labelsize': 10})

nx = ny = 201
nz = 101
d = 0.06
# Origin of the model in UTM zone 10N (km)
x0 = 588.5
y0 = 5183.5


def load(path, shape):
    """Read an OWL array, whose first axis is the fastest, into a numpy array indexed as in Fortran."""
    return np.fromfile(path, np.float32).reshape(shape[::-1]).T


def minus(v):
    """Format a number with a typographic minus sign."""
    return f'{v:.1f}'.replace('-', '−')


topo = load('model/topo.bin', (ny, nx))/1000.0
mask = load('model/mask.bin', (nz, ny, nx))
vp = load('model/vp.bin', (nz, ny, nx))
iy, ix = np.unravel_index(np.argmax(topo), topo.shape)
hmax = topo.max()

src = np.array([[x0 + 1.0 + 2.0*j, y0 + 1.0 + 2.0*i] for j in range(6) for i in range(6)])
rec = np.array([[x0 + 0.3*j, y0 + 0.3*i] for j in range(1, 40) for i in range(1, 40)])

fig = plt.figure(figsize=(7.0, 7.0), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1.75, 1, 1])
extent_map = [x0 - 0.5*d, x0 + (nx - 0.5)*d, y0 - 0.5*d, y0 + (ny - 0.5)*d]
xticks = [590, 595, 600]
yticks = [5185, 5190, 5195]

# Topography with the sources and receivers; the land part of the terrain color map
land = LinearSegmentedColormap.from_list('land', plt.cm.terrain(np.linspace(0.25, 1.0, 256)))
ax = fig.add_subplot(gs[0, 0])
ls = LightSource(azdeg=315, altdeg=35)
ax.imshow(ls.shade(topo, cmap=land, vert_exag=1, dx=d, dy=d, blend_mode='soft', vmin=1.2, vmax=4.4),
          origin='lower', extent=extent_map)
ax.contour(x0 + np.arange(nx)*d, y0 + np.arange(ny)*d, topo, levels=np.arange(1.5, 4.5, 0.5), colors='k', linewidths=0.4)
ax.plot(rec[:, 0], rec[:, 1], '.', color='0.15', ms=0.6, alpha=0.6)
ax.plot(src[:, 0], src[:, 1], '*', color='r', ms=6, mec='none')
ax.axvline(x0 + ix*d, color='w', lw=0.8, ls='--')
ax.set_xlabel('Easting (km)')
ax.set_ylabel('Northing (km)')
ax.set_title('(a) Topography', loc='left')
ax.set_xticks(xticks)
ax.set_yticks(yticks)
cb = fig.colorbar(ScalarMappable(Normalize(1.2, 4.4), land), ax=ax, shrink=0.7, aspect=15, pad=0.01)
cb.set_label('Elevation (km)')

# Horizontal slice of the true P-wave velocity at an elevation of 1 km
k = int(round((hmax - 1.0)/d))
ax = fig.add_subplot(gs[0, 1])
im = ax.imshow(vp[k, :, :], origin='lower', extent=extent_map, cmap='jet', vmin=2000, vmax=6000, interpolation='none')
ax.axvline(x0 + ix*d, color='w', lw=0.8, ls='--')
ax.set_xlabel('Easting (km)')
ax.set_title(f'(b) P-wave velocity at {minus(hmax - k*d)} km', loc='left')
ax.set_xticks(xticks)
ax.set_yticks(yticks)
cb = fig.colorbar(im, ax=ax, shrink=0.7, aspect=15, pad=0.01)
cb.set_label('m/s')

# North-south slices of the true and initial models through the summit
extent = [y0 - 0.5*d, y0 + (ny - 0.5)*d, hmax - (nz - 0.5)*d, hmax + 0.5*d]
panels = [('vp.bin', '(c) True P-wave velocity', 2000, 6000),
          ('vs.bin', '(d) True S-wave velocity', 1100, 3400),
          ('vp_init.bin', '(e) Initial P-wave velocity', 2000, 6000),
          ('vs_init.bin', '(f) Initial S-wave velocity', 1100, 3400)]
for n, (f, title, vmin, vmax) in enumerate(panels):
    ax = fig.add_subplot(gs[1 + n//2, n % 2])
    v = load(f'model/{f}', (nz, ny, nx))[:, :, ix]
    v = np.where(mask[:, :, ix] > 0, v, np.nan)
    im = ax.imshow(v, extent=extent, cmap='jet', vmin=vmin, vmax=vmax, interpolation='none')
    ax.set_title(title, loc='left')
    ax.set_xticks(yticks)
    ax.set_yticks([-1, 0, 1, 2, 3, 4])
    ax.set_yticklabels([minus(t).replace('.0', '') for t in [-1, 0, 1, 2, 3, 4]])
    if n % 2 == 0:
        ax.set_ylabel('Elevation (km)')
    if n >= 2:
        ax.set_xlabel('Northing (km)')
    cb = fig.colorbar(im, ax=ax, shrink=0.9, aspect=10, pad=0.01)
    cb.set_label('m/s')

Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/model.pdf')
fig.savefig('figures/model.png', dpi=200)
