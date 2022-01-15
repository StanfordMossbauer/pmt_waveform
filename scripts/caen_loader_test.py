import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen
import physics_util as pu

from bead_util import get_color_map

rng = np.random.default_rng()

# data_dir = '/home/cblakemore/caen_data/PlasticScintillatorTestForChas/'
data_dir = '/Users/manifestation/tmp/'

data_files = \
    [
     'no_source_background.txt', \
     'no_source_background_dark.txt', \
     'cs137_dark.txt', \
     'ba133_dark.txt'
    ]

first_n_waveform = 0
# first_n_waveform = 70000

cs137_data = caen.WaveformContainer(os.path.join(data_dir, data_files[3]), \
                                    first_n_waveform=first_n_waveform)
cs137_data.subtract_baseline(pulse_start_ind=50)

cs137_data.find_pulse_maxima_mean(sample_window=10, presmooth=True)


fig, ax = plt.subplots(1,1)

vals, _, _ = ax.hist(cs137_data.pulse_maxima, bins=100, range=(0,100))

ax.set_yscale('log')
ax.set_ylim(0.5, 2*np.max(vals))
ax.set_ylabel('Counts [abs.]')
ax.set_xlabel('Pulse Maximum [ADC bins]')

plt.show()


# cs137_data.find_pulse_maxima_gauss(sample_window=30)

# for i in range(10):
#     plt.plot(cs137_data.waveform_arr_nomean[i,:], color='C{:d}'.format(i))
#     plt.axhline(cs137_data.pulse_maxima[i], color='C{:d}'.format(i))
# plt.show()



# test_ind = rng.choice(np.arange(cs137_data.waveform_arr.shape[0]))

# popt = pu.fitting.generate_histogram_and_fit_gaussian( \
#                         cs137_data.waveform_arr[:,:50], bins=15, \
#                         range=(8276,8290), plot=True, print_level=0)

