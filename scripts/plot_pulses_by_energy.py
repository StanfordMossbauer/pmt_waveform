import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen

plt.rcParams.update({'font.size': 14})

rng = np.random.default_rng()

base = '/Users/manifestation/Stanford/mossbauer/darkbox_clone/'

load_class = True

data_dir = base + '20220527/data'
data_files = \
    [
     'am241_600V.dat', \
     # 'cs137_600V.dat', \
     # 'ba133_600V.dat', \
     # 'background_600V.dat', \
    ]

colors = ['C0', 'C1', 'C2', 'C3']
labels = ['Am-241', 'Cs-137', 'Ba-133', 'background']






n_spectra = len(data_files)


for i in range(n_spectra):

    data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[i]), \
                                  header=True)

    if load_class:
        data.load(verbose=True)

    else:
        data.compute_baseline(pulse_start_ind=50)

        data.find_pulse_extrema(polarity=-1)

        data.integrate_waveforms(integration_start=1.5e-6, integration_window=3.0e-6, \
                                 baseline=False, plot=False, plot_index=111, \
                                 polarity=-1, adaptive_window=True, asymmetry=0.1)

        data.save(verbose=True)


    fig, ax = plt.subplots(1,1)
    fig, ax = data.plot_pulse_spectra(hist_range=(0,40000), nbin=2000, filled=False, \
                                      fig=fig, ax=ax, show=True, color=colors[i], \
                                      label=labels[i], log_scale=True)

    energy_to_plot = float(input('Energy to plot: '))
    pulses_to_plot = int(input('Pulses to plot: '))

    energy_sorter = np.argsort(np.abs(data.waveform_integrals - energy_to_plot))
    colors = caen.get_color_map(pulses_to_plot, cmap='plasma')

    fig2, ax2 = plt.subplots(1,1)
    for i in range(pulses_to_plot):
        color = colors[i]
        fig2, ax2 = data.plot_waveform(energy_sorter[i], baseline=False, \
                                       amp_scale='v', tscale='u', show=False, \
                                       fig=fig2, ax=ax2, color=color)
    plt.show()

    del data

plt.show()


