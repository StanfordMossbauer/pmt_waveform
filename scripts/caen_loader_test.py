import os, timeit

import numpy as np
import matplotlib.pyplot as plt

import caen_loader as caen

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
                              filetype='compass')

