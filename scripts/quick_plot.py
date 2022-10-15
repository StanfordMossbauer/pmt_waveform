import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen

rng = np.random.default_rng()

# base = '/Users/manifestation/Stanford/mossbauer/darkbox_clone/'
base = '/Users/manifestation/tmp/'
# base = '/home/cblakemore/caen_data/'

# data_dir = base + '20220510_CAENTest/data'
# data_dir = base + '20220513_preampTest/data'
# data_dir = base + '20220517_preampCLC/data'
# data_dir = base + '20220520/data'
# data_dir = base + '20220527/data'
# data_dir = base + '20220605/data'
# data_dir = base + '20220607/data'

data_dir = base

data_files = \
    [
     # 'background_1600V.dat', \
     # 'am241_800V_13300trig.dat', \
     # 'baseline_test.dat', \
     # '../wave1.dat', \
     'cs137_20dB.dat', \
    ]


n_to_plot = 5

rng = np.random.default_rng()

# waveforms_to_plot = range(5)
# waveforms_to_plot = [11870, 44187, 58880]

data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[0]), \
                              filetype='compass')

data.compute_baseline(pulse_start_ind=50)

data.plot_baseline(amp_scale='b', mean_subtract=False)

waveforms_to_plot = rng.integers(low=0, high=data.n_waveform, size=n_to_plot)
print(waveforms_to_plot)
for i in waveforms_to_plot:
    data.plot_waveform(i, amp_scale='volts', baseline=False)




