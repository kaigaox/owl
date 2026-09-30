#
# Add band-limited Gaussian noise to the seismograms
#
import numpy as np
from pathlib import Path

ns = 20
f0 = 8.0

# Root-mean-square (RMS) amplitude of the noise relative to that of each shot gather
level = 0.5

rng = np.random.default_rng(2026)
Path('data_noisy').mkdir(exist_ok=True)
for i in range(1, ns + 1):

    # A SU file is a sequence of traces, each with a 240-byte header followed by
    # float32 samples; the header stores the number of samples and the sampling
    # interval (in microseconds) at bytes 114 and 116
    raw = np.fromfile(f'data/shot_{i}_seismogram_p.su', np.uint8)
    nt = int(raw[114:116].view(np.int16)[0])
    dt = int(raw[116:118].view(np.int16)[0])*1.0e-6
    traces = raw.reshape(-1, 240 + 4*nt)
    d = traces[:, 240:].copy().view(np.float32).astype(np.float64)

    # Shape white noise with the amplitude spectrum of the Ricker wavelet,
    # so that the noise lies in the frequency band of the data
    f = np.fft.rfftfreq(nt, dt)
    w = (f/f0)**2*np.exp(-(f/f0)**2)
    n = np.fft.irfft(np.fft.rfft(rng.standard_normal(d.shape), axis=1)*w, n=nt, axis=1)
    n *= level*np.sqrt(np.mean(d**2)/np.mean(n**2))

    traces[:, 240:] = (d + n).astype(np.float32).view(np.uint8)
    traces.tofile(f'data_noisy/shot_{i}_seismogram_p.su')
