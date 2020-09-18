# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:26:02 2020

@author: lapi2
"""
import scipy.io as sio
from scipy import signal
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('YOURDRIVE:\\YOURPATH') # change directory
filename = 'YOURDRIVE:\\YOURPATH\\YOURMATFILE.mat';
mat = sio.loadmat(filename);

a = mat['EEG'];

canales = a.item(0)[20] # canales
nsamples = int(a.item(0)[10]) # nsamples
sampling_frequency = int(a.item(0)[11]) # samplingfreq

n_electrodes = len(a.item(0)[15])
elec = np.arange(n_electrodes,dtype='int')


exp1 = a.item(0)[15][elec[0:n_electrodes],0:nsamples] # data

