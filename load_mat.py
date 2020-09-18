# -*- coding: utf-8 -*-
"""
Created on Sat May  2 21:46:13 2020

@author: lapi2
"""

import scipy.io as sio
import os
import numpy as np


def load_mat(osdir,matname):
    os.chdir(osdir)
    mat = sio.loadmat(matname)

    data = mat['exp1'];
    channels = mat['canales']; # canales
    nsamples = mat['nsamples']; # nsamples
    sampling_frequency = mat['sampling_frequency']; # samplingfreq

    n_electrodes = int(mat['n_electrodes']);
    elec = np.arange(n_electrodes,dtype='int')
    return data,channels,nsamples,sampling_frequency,n_electrodes,elec

