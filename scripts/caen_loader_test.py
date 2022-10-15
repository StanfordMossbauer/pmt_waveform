import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen

plt.rcParams.update({'font.size': 14})

rng = np.random.default_rng()

base = '/Users/manifestation/tmp/'

# data_dir = base + '20220510_CAENTest/data'
# data_dir = base + '20220513_preampTest/data'
# data_dir = base + '20220517_preampCLC/data'
# data_dir = base + '20220519/data'

data_dir = base

data_files = \
    [
     'cs137_20dB.dat', \
    ]


data = caen.WaveformContainer(fname=os.path.join(data_dir, data_files[0]), \
                              filetype='compass')

