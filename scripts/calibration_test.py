import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen
import physics_util as pu

plt.rcParams.update({'font.size': 14})

rng = np.random.default_rng()

base = '/Users/manifestation/Stanford/mossbauer/darkbox_clone/'


# data_dir = base + '20220519/data'

# data_files = \
#     [
#      'background_1600V.dat', \
#      'ba133_1600V_far.dat', \
#      # 'ba133_1600V_close.dat', \
#      'cs137_1600V_far.dat', \
#      # 'cs137_1600V_close.dat', \
#     ]

# colors = ['C0', 'C1', 'C2']
# labels = ['background', 'Ba-133', 'Cs-137']



data_dir = base + '20220513_preampTest/data'

data_files = \
    [
     'ba133_1600V_ortecVT120.dat', \
     'ba133_1600V_CLC144.dat', \
    ]

colors = ['C0', 'C1']
labels = ['Ortec VT120', 'CLC144']


n_spectra = len(data_files)

show = False
fig, ax = plt.subplots(1,1)

for i in range(n_spectra):
    if i == n_spectra - 1:
        show = True

    data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[i]), \
                                  header=True)

    data.compute_baseline(pulse_start_ind=50)

    data.integrate_waveforms(integration_start=1.5e-6, integration_window=0.6e-6, \
                             baseline=False, plot=False, plot_index=111, \
                             polarity=-1, adaptive_window=True, asymmetry=0.3)

    fig, ax = data.plot_pulse_spectra(hist_range=(0,8000), filled=False, \
                                      fig=fig, ax=ax, show=show, color=colors[i], 
                                      label=labels[i])


    del data


