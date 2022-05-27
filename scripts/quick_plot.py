import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen
import physics_util as pu

rng = np.random.default_rng()

base = '/Users/manifestation/Stanford/mossbauer/darkbox_clone/'

# data_dir = base + '20220510_CAENTest/data'
# data_dir = base + '20220513_preampTest/data'
# data_dir = base + '20220517_preampCLC/data'
data_dir = base + '20220520/data'

data_files = \
    [
     # 'ba133_1600V_ortecVT120.dat', \
     # 'ba133_1600V_CLC144.dat', \
     # 'background_1600V.dat', \
     'test2.dat', \
     # '../wave0.dat', \
    ]


n_to_plot = 15

rng = np.random.default_rng()

# waveforms_to_plot = range(5)
# waveforms_to_plot = [11870, 44187, 58880]

data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[0]), \
                              header=True)

data.compute_baseline(pulse_start_ind=50)
# data.plot_baseline(amp_scale='b', mean_subtract=False)

waveforms_to_plot = rng.integers(low=0, high=data.n_waveform, size=n_to_plot)
print(waveforms_to_plot)
for i in waveforms_to_plot:
    data.plot_waveform(i, amp_scale='bits', baseline=True)




