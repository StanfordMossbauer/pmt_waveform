import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen
import physics_util as pu

plt.rcParams.update({'font.size': 14})

rng = np.random.default_rng()

base = '/Users/manifestation/Stanford/mossbauer/darkbox_clone/'

# data_dir = base + '20220510_CAENTest/data'
data_dir = base + '20220513_preampTest/data'
# data_dir = base + '20220517_preampCLC/data'
# data_dir = base + '20220519/data'

data_files = \
    [
     # 'ba133_1600V_ortecVT120.dat', \
     'ba133_1600V_CLC144.dat', \
     # 'background_1600V.dat', \
     # 'test.dat', \
     # 'test2.dat', \
     # '../wave0.dat', \
    ]


data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[0]), \
                              header=True)

data.compute_baseline(pulse_start_ind=50, plot=False, amp_scale='v')

# for i in range(100):
#     ind = int(222 + 2*i)
#     data.plot_waveform(ind, baseline=True, amp_scale='b')
#     input()

test_index = 573
data.integrate_waveforms(integration_start=1.5e-6, integration_window=0.9e-6, \
                         baseline=False, plot=True, plot_index=222, \
                         polarity=-1, adaptive_window=True, asymmetry=0.2)

# a, b = data._load_caen_binary(first_index=test_index, last_index=test_index+1)
# index = a[0]
# waveform = b[0]

# print(index)
# # data.plot_waveform(index, baseline=False, amp_scale='v')

# plt.figure()
# plt.loglog(np.abs(np.fft.rfft(waveform - data.baseline_arr[index])))

# plt.show()

data.plot_pulse_spectra(hist_range=(0,5000), filled=False)


