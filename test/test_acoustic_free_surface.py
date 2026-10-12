"""
Test 10: acoustic free surface against the analytical half-space solution.

owl_modeling2 with acoustic-iso and a flat free surface (p = 0 on the first grid row) in a
homogeneous water half-space is compared in absolute amplitude with the exact solution given by
the method of images: the 2D whole-space line-source field of the source minus (explosion) or
plus (vertical force) that of its mirror image at (xs, -zs). The whole-space fields are evaluated
in the frequency domain (numpy FFT convention, e^{+i w t} time dependence, outgoing waves are
Hankel functions of the second kind):

    explosion (OWL adds dt*amp*s(t) to p)
        p^ = dx dz amp (i w) S^(w) (-i/(4 c^2)) H0^(2)(w r/c)
    vertical force (OWL adds dt*amp*s(t)/rho to vz)
        p^ = -(i w/(4 c)) dx dz amp S^(w) H1^(2)(w r/c) (z - zs)/r

Marine-streamer-like acquisition: explosions at 2.5 m (half a cell), 6 m (airgun depth, off the
grid) and 50 m, a vertical force at 6 m, and receivers at 0, 8, 15, 100, 500 and 800 m depth with 100-500 m
offset. The sources and the shallow receivers lie inside the 4-point interpolation window of the
surface, so the image-source injection, the mirrored receiver sampling and the mirrored vz rows
of the free-surface update all matter; the 500 m and 800 m lines check the far field, where the
dispersion of the interior scheme accumulates. The 6 m explosion is repeated on a twice finer
grid and time step.

Time alignment: OWL records p before the update of step t, updates p then v, and injects the
source sample s(t) (nominal time (t-1) dt) afterwards, so the v update never sees the current
injection. Eliminating v from the recurrence of the recorded pressure R_t gives

    (R_{t+1} - 2 R_t + R_{t-1})/dt^2 + c^2 k^2 R_t = (q_t - q_{t-1})/dt + dt c^2 k^2 q_{t-1},

whose right-hand side on the numerical dispersion surface c^2 k^2 = (2/dt)^2 sin^2(w dt/2) is
exactly i Omega q^ e^{-1.5 i w dt}: the explosion acts exactly 1.5 dt after its nominal time. The
same algebra gives exactly 1.0 dt for a force injected into v. The Ricker is lengthened to 3/f0
(f0_factor = 1.5): OWL's default truncation at +-1/f0, where the Ricker is still 1e-3 of its peak,
turns the two jumps into impulses of the source derivative that excite the 1/sqrt(t - r/c)
singularity of the 2D Green's function in the exact reference but not on the grid, and that
artifact dominates the ghost-dominated shallow traces.

Pass criteria: relative L2 error per trace below 4% on the 5 m grid for both mechanisms, the
surface receivers below 1e-6 of the peak pressure, the fine-grid error below 0.7 of the
coarse-grid error (the residual is the interior scheme's dispersion and must vanish with
refinement), and the snapshot difference below 3%. Observed: 0.2-2.3% (explosion) and 0.2-1.6%
(force), growing linearly with the propagation distance, fine/coarse 0.34, exactly 0 at the
surface, 1.05% for the snapshot.
"""

import os
import sys
import shutil
import numpy as np
from scipy.special import hankel2
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import write_model, read_model, read_su, write_param, run_owl, report

WORK = os.path.join(HERE, 'work', 'acoustic_free_surface')
PLOT = os.path.join(HERE, 'plots')

# Water half-space and marine-like acquisition
VP, RHO = 1500.0, 1000.0
F0, F0_FACTOR, AMP = 15.0, 1.5, 1.0e6
D = 5.0                      # dx = dz (m)
NX, NZ = 401, 241            # 2000 m x 1200 m
DT, TMAX = 5.0e-4, 1.0
PML = 30                     # cells (150 m)
XS = 750.0
OFFSETS = np.arange(100.0, 501.0, 100.0)
REC_Z = [0.0, 8.0, 15.0, 100.0, 500.0, 800.0]
SRC_Z = [2.5, 6.0, 50.0]
SNAP_T = 0.4
TOL, TOL_SURFACE, TOL_REFINE, TOL_SNAP = 0.04, 1.0e-6, 0.7, 0.03


def ricker_owl(nt, dt):
    """OWL's Ricker: ceil(2 f0_factor/(f0 dt)) samples at (t-1) dt, peak at f0_factor/f0, unit maximum."""
    nstf = int(np.ceil(2.0*F0_FACTOR/F0/dt))
    tt = np.arange(nstf)*dt
    a = (np.pi*F0*(tt - F0_FACTOR/F0))**2
    s = (1.0 - 2.0*a)*np.exp(-a)
    s /= np.max(np.abs(s))
    out = np.zeros(nt)
    out[:nstf] = s[:nt]
    return out


def delay_of(mech, dt):
    """Delay of OWL's source relative to its nominal time (see the module docstring)."""
    return 1.5*dt if mech == 'explosion' else dt


def whole_space_spectrum(xr, zr, xs, zs, w, S, mech, d):
    """Pressure spectrum (nf, nrec) of a line source at (xs, zs) in a whole space."""
    r = np.hypot(xr - xs, zr - zs)[None, :]
    k = (w/VP)[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        if mech == 'explosion':
            P = d*d*AMP*(1j*w)[:, None]*S[:, None]*(-1j/(4.0*VP**2))*hankel2(0, k*r)
        else:
            P = -(1j*k/4.0)*d*d*AMP*S[:, None]*hankel2(1, k*r)*((zr - zs)[None, :]/r)
    P[0, :] = 0.0
    return P


def half_space_spectrum(xr, zr, zs, mech, w, S, d):
    """Source at (XS, zs) plus its image: p is odd about the surface, so the explosion image is
    opposite to the source and the vertical-force image (vz even) equals it."""
    sign = -1.0 if mech == 'explosion' else 1.0
    return (whole_space_spectrum(xr, zr, XS, zs, w, S, mech, d)
            + sign*whole_space_spectrum(xr, zr, XS, -zs, w, S, mech, d))


def analytic(xr, zr, zs, mech, nt, dt, d, delay=None, pad=8):
    """Half-space pressure (nt, nrec) at the receivers (xr, zr) from a source at (XS, zs)."""
    npad = pad*nt
    f = np.fft.rfftfreq(npad, dt)
    w = 2.0*np.pi*f
    S = np.fft.rfft(ricker_owl(npad, dt))*dt
    delay = delay_of(mech, dt) if delay is None else delay
    P = half_space_spectrum(xr, zr, zs, mech, w, S, d)*np.exp(-1j*w*delay)[:, None]
    return np.fft.irfft(P, n=npad, axis=0)[:nt]/dt


def receivers():
    xr = np.concatenate([XS + OFFSETS for _ in REC_Z])
    zr = np.concatenate([np.full(len(OFFSETS), z) for z in REC_Z])
    return xr, zr


def cols(iz):
    n = len(OFFSETS)
    return slice(iz*n, (iz + 1)*n)


def rel_l2(a, b):
    return np.linalg.norm(a - b)/np.linalg.norm(b)


def owl_run(tag, zs, mech, d, nx, nz, dt, pml, snap):
    """One owl_modeling2 shot; returns (data (nt, nrec), pressure snapshot at SNAP_T or None)."""
    wd = os.path.join(WORK, tag)
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(os.path.join(wd, 'geometry'))
    write_model(os.path.join(wd, 'model', 'vp.bin'), np.full((nz, nx), VP))
    write_model(os.path.join(wd, 'model', 'rho.bin'), np.full((nz, nx), RHO))
    with open(os.path.join(wd, 'geometry', 'geometry.txt'), 'w') as f:
        f.write('shot_1_geometry.txt\n')
    with open(os.path.join(wd, 'geometry', 'shot_1_geometry.txt'), 'w') as f:
        f.write(f'1\n1\n{XS} 0.0 {zs}\n{"explosion" if mech == "explosion" else "force 0 0"}\n'
                f'ricker {F0} {AMP} 0.0\n0 0\n{len(OFFSETS)*len(REC_Z)}\n')
        for rz in REC_Z:
            for r in OFFSETS:
                f.write(f'{XS + r} 0.0 {rz} 1.0\n')
    par = {'nx': nx, 'nz': nz, 'dx': d, 'dz': d, 'dt': dt, 'data_dt': dt, 'tmax': TMAX, 'ns': 1,
           'file_geometry': './geometry/geometry.txt', 'which_medium': 'acoustic-iso',
           'model_name': 'vp, rho', 'file_vp': './model/vp.bin', 'file_rho': './model/rho.bin',
           'npml': pml, 'yn_free_surface': 'y', 'f0_factor': F0_FACTOR, 'dir_synthetic': 'data',
           'dir_snapshot': 'snapshot', 'verbose': 'n'}
    if snap:
        # snapshot times as begin, interval, end; the second one is at SNAP_T
        par['snaps'] = f'{0.5*SNAP_T:g}, {0.5*SNAP_T:g}, {SNAP_T:g}'
    write_param(os.path.join(wd, 'param.rb'), par)
    run_owl('owl_modeling2', 'param.rb', wd)
    data = read_su(os.path.join(wd, 'data', 'shot_1_seismogram_p.su'))[0].astype(np.float64)
    field = None
    if snap:
        field = read_model(os.path.join(wd, 'snapshot', 'shot_1_forward_wavefield_p_2.bin'),
                           (nz, nx)).astype(np.float64)
    return data, field


def compare(mech, zs, d=D, nx=NX, nz=NZ, dt=DT, pml=PML, snap=False):
    """Run OWL for one source and compare with the reference; returns a results dict."""
    xr, zr = receivers()
    owl, field = owl_run(f'{mech}_{zs:g}m_dx{d:g}m', zs, mech, d, nx, nz, dt, pml, snap)
    nt = owl.shape[0]
    ref = analytic(xr, zr, zs, mech, nt, dt, d)
    errs = np.full(len(xr), np.nan)
    for ir in range(len(xr)):
        if zr[ir] > 0.0:
            errs[ir] = rel_l2(owl[:, ir], ref[:, ir])
    surface = np.max(np.abs(owl[:, cols(0)]))/np.max(np.abs(owl))
    # diagnostic: the best-fitting source delay over the receivers between the surface and 100 m,
    # where the accumulated dispersion is smallest
    sub = (zr > 0.0) & (zr <= 100.0)
    taus = np.linspace(0.0, 2.0*dt, 41)
    e = [rel_l2(owl[:, sub], analytic(xr[sub], zr[sub], zs, mech, nt, dt, d, tau)) for tau in taus]
    print(f'  acoustic-iso, {mech} at {zs:g} m, dx = {d:g} m: surface |p| / max |p| = {surface:.1e}, '
          f'best-fit source delay {taus[int(np.argmin(e))]/dt:.2f} dt (theory {delay_of(mech, dt)/dt:.1f} dt)')
    for iz, rz in enumerate(REC_Z):
        if rz > 0.0:
            print(f'    receivers at {rz:5.1f} m: relative L2 error per offset '
                  + ' '.join(f'{v*100:5.2f}%' for v in errs[cols(iz)]))
    return dict(mech=mech, zs=zs, d=d, nx=nx, nz=nz, dt=dt, nt=nt, owl=owl, ref=ref, errs=errs,
                surface=surface, snap=field)


def snapshot_reference(r, X, Z):
    """Analytical pressure at t = SNAP_T on the points (X, Z), from frequencies up to 80 Hz."""
    nt, dt, d = r['nt'], r['dt'], r['d']
    npad = 8*nt
    f = np.fft.rfftfreq(npad, dt)
    w = 2.0*np.pi*f
    S = np.fft.rfft(ricker_owl(npad, dt))*dt
    nf = int(np.searchsorted(f, 80.0))
    P = half_space_spectrum(X.ravel(), Z.ravel(), r['zs'], r['mech'], w[:nf], S[:nf], d)
    P *= np.exp(-1j*w[:nf]*(delay_of(r['mech'], dt) - SNAP_T))[:, None]
    return ((2.0*np.real(np.sum(P[1:], axis=0)) + np.real(P[0]))/(npad*dt)).reshape(X.shape)


def plot_traces(res_list, passed):
    rows = [8.0, 800.0]
    fig, axs = plt.subplots(len(rows), len(res_list), figsize=(3.5*len(res_list) + 0.5, 3.3*len(rows) + 0.6),
                            sharex=True, squeeze=False)
    for j, r in enumerate(res_list):
        t = np.arange(r['nt'])*r['dt']
        for i, rz in enumerate(rows):
            ax = axs[i, j]
            iz = REC_Z.index(rz)
            for k in range(len(OFFSETS)):
                ir = iz*len(OFFSETS) + k
                a = np.max(np.abs(r['ref'][:, ir]))
                ax.plot(t, 2.4*k + r['owl'][:, ir]/a, '-', color='blue', lw=1.3)
                ax.plot(t, 2.4*k + r['ref'][:, ir]/a, '--', color='red', lw=1.1)
                ax.plot(t, 2.4*k + 10.0*(r['owl'][:, ir] - r['ref'][:, ir])/a, '-', color='black', lw=0.7)
            ax.set_yticks(2.4*np.arange(len(OFFSETS)))
            ax.set_yticklabels([f'{o:g}' for o in OFFSETS])
            ax.set_xlim(0.0, TMAX)
            ax.set_ylim(-1.6, 2.4*(len(OFFSETS) - 1) + 1.6)
            ax.set_title(f'{r["mech"]} at {r["zs"]:g} m, receivers at {rz:g} m', fontsize=9)
            if j == 0:
                ax.set_ylabel('Offset (m)')
            if i == len(rows) - 1:
                ax.set_xlabel('Time (s)')
    axs[0, 0].plot([], [], '-', color='blue', lw=1.3, label='OWL')
    axs[0, 0].plot([], [], '--', color='red', lw=1.1, label='Analytical')
    axs[0, 0].plot([], [], '-', color='black', lw=0.7, label='Difference x 10')
    axs[0, 0].legend(loc='upper left', fontsize=8)
    fig.suptitle(f'Test 10: acoustic free surface vs analytical half-space '
                 f'(water, {F0:g} Hz Ricker, dx = {D:g} m) [{"PASSED" if passed else "FAILED"}]')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test10_acoustic_free_surface.png'), dpi=150)
    plt.close(fig)


def plot_errors(res_list, fine):
    panels = res_list + [fine]
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 7.0), sharey=True)
    colors = ['blue', 'red', 'black', 'green', 'magenta']
    for ax, r in zip(axs.ravel(), panels):
        for i, rz in enumerate([z for z in REC_Z if z > 0.0]):
            ax.plot(OFFSETS, 100.0*r['errs'][cols(REC_Z.index(rz))], '-o', color=colors[i], ms=4,
                    label=f'receivers {rz:g} m')
        ax.axhline(100.0*TOL, color='0.5', ls=':')
        ax.set_title(f'{r["mech"]} at {r["zs"]:g} m, dx = {r["d"]:g} m, dt = {1e3*r["dt"]:g} ms', fontsize=10)
        ax.grid(True, alpha=0.3)
    for ax in axs.ravel()[len(panels):]:
        ax.axis('off')
    for ax in axs[:, 0]:
        ax.set_ylabel('Relative L2 error (%)')
    for ax in axs[1, :len(panels) - 3]:
        ax.set_xlabel('Offset (m)')
    for ax in axs[0, :]:
        ax.set_xlabel('Offset (m)')
    axs[0, 0].legend(fontsize=8)
    fig.suptitle('Test 10: relative L2 error per trace (dotted: tolerance)')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test10_acoustic_free_surface_errors.png'), dpi=150)
    plt.close(fig)


def plot_snapshot(r, zmax=700.0, xlim=(XS - 600.0, XS + 600.0)):
    """OWL snapshot, analytical field and their difference at t = SNAP_T; returns the relative L2
    difference outside the source vicinity."""
    d, nx, nz = r['d'], r['nx'], r['nz']
    x = np.arange(nx)*d
    z = np.arange(nz)*d
    iz = np.where(z <= zmax)[0]
    ix = np.where((x >= xlim[0]) & (x <= xlim[1]))[0]
    X, Z = np.meshgrid(x[ix], z[iz])
    ana = snapshot_reference(r, X, Z)
    ana[np.hypot(X - XS, Z - r['zs']) < 4.0*d] = np.nan
    owl = r['snap'][np.ix_(iz, ix)]
    mask = np.isfinite(ana)
    err = np.linalg.norm((owl - ana)[mask])/np.linalg.norm(ana[mask])
    vmax = np.nanmax(np.abs(ana))
    fig, axs = plt.subplots(3, 1, figsize=(8.0, 8.6), sharex=True)
    ext = [x[ix][0], x[ix][-1], z[iz][-1], z[iz][0]]
    for ax, fld, lab, scale in zip(axs, [owl, ana, owl - ana], ['OWL', 'Analytical', 'Difference x 20'],
                                   [1.0, 1.0, 20.0]):
        ax.imshow(scale*fld, extent=ext, cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal',
                  interpolation='nearest')
        ax.plot([XS], [r['zs']], 'o', color='black', ms=4)
        ax.set_ylabel('z (m)')
        ax.set_title(lab, fontsize=10)
    axs[-1].set_xlabel('x (m)')
    fig.suptitle(f'Test 10: {r["mech"]} at {r["zs"]:g} m, pressure at t = {SNAP_T:g} s\n'
                 f'relative L2 difference {err*100:.2f}% (source vicinity excluded)')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test10_acoustic_free_surface_snapshot.png'), dpi=150)
    plt.close(fig)
    return err


def run():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(PLOT, exist_ok=True)

    res_list = [compare('explosion', zs, snap=(zs == 6.0)) for zs in SRC_Z]
    res_list.append(compare('force', 6.0))
    coarse = res_list[SRC_Z.index(6.0)]
    fine = compare('explosion', 6.0, d=D/2, nx=2*(NX - 1) + 1, nz=2*(NZ - 1) + 1, dt=DT/2, pml=2*PML)

    worst = max(np.nanmax(r['errs']) for r in res_list)
    surface = max(r['surface'] for r in res_list)
    ratio = np.nanmean(fine['errs'])/np.nanmean(coarse['errs'])
    passed_so_far = worst <= TOL and surface <= TOL_SURFACE and ratio <= TOL_REFINE
    snap_err = plot_snapshot(coarse)
    passed = passed_so_far and snap_err <= TOL_SNAP
    print(f'  refinement: mean error fine/coarse = {ratio:.2f}; snapshot difference {snap_err*100:.2f}%')

    plot_traces(res_list, passed)
    plot_errors(res_list, fine)
    report('Acoustic free surface vs analytical half-space solution', passed,
           f'max trace error {worst*100:.2f}% (tol {TOL*100:g}%), surface |p| ratio {surface:.1e} '
           f'(tol {TOL_SURFACE:g}), fine/coarse error {ratio:.2f} (tol {TOL_REFINE:g}), '
           f'snapshot difference {snap_err*100:.2f}% (tol {TOL_SNAP*100:g}%)')
    return passed, dict(worst=worst, surface=surface, ratio=ratio, snapshot=snap_err)


if __name__ == '__main__':
    sys.exit(0 if run()[0] else 1)
