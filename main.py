# -*- coding: utf-8 -*-
"""
Created on Sat May  2 21:55:41 2020

@author: lapi2
"""

from load_mat import *
from Im_Coherency import coherency

def main():
    data,channels,nsamples,sampling_frequency,n_electrodes,elec = load_mat('D:\\APP_BrainConn\\Ezequiel_ImgCoh','EXPERIMENTS_EEG_FOR_GITHUB.mat');
    Coh = coherency(int(sampling_frequency),data)
    return Coh

Coh = main()

