# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:31:23 2018

@author: lapi2
"""
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

"""
500 Hz sampling rate
100 samples/2=50 FFT bins
250 Hz/50 FFT bins≃5 Hz/bin
"""

def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n


def coherency(sampling_frequency,data):
    nfft = 2^nextpow2(len(data))
    timestep = 0.1
    Fs = sampling_frequency 
    dt = 1/Fs

    nchan, segleng = data.shape # row or nchan = 64 channels, columns or segleng = all sampling points
    t = np.arange(0,segleng/Fs,dt) 

    channels_number = np.arange(0,nchan)

    count = 0
    incount = 1
    j = 0
    k = 250
    window = signal.hann(250)
    frq = np.arange(0,Fs/2)*0.5
    N = len(window)
    positions = np.arange(len(window))
    T = len(window)/Fs

    frq = np.arange(0,Fs/4) # one side frequency range
    Sxy = np.zeros(250, dtype = complex)
    Sxx = np.zeros(250, dtype = complex)
    Syy = np.zeros(250, dtype = complex)
    Coh = np.zeros((nchan,250,nchan), dtype = complex)


    electrode = 0
    TRIO = nchan - 1 
    
    while TRIO > 0:
        
        for index in range(electrode+1,nchan):
            
            while (k < segleng):
                
                s1 = signal.detrend(data[electrode,j:k])
                s2 = signal.detrend(data[index,j:k])
                temp = (s1)*window
                Fx = np.fft.fft(temp)/len(temp)
                temp = (s2)*window        
                Fy = np.fft.fft(temp)/len(temp)
                # Fx = Fx ./ abs(Fx) for PLV
                # Fy = Fy ./ abs(Fy) for PLV
                Sxy += Fx * np.conjugate(Fy)
                Syy += Fy * np.conjugate(Fy)
                Sxx += Fx * np.conjugate(Fx)
                incount += 1
                j += 250
                k += 250
                
        # PLV = Sxy /  # Sxx*Syy
        Sxy_mean = Sxy / incount
        Sxx_mean = Sxx / incount
        Syy_mean = Syy / incount
        coherency_matrix = Sxy_mean / np.sqrt(Sxx_mean*Syy_mean)
        Coh[electrode,:,index] = coherency_matrix;
            
        Sxy = np.zeros(250, dtype = complex)
        Sxx = np.zeros(250, dtype = complex)
        Syy = np.zeros(250, dtype = complex)
        coherency_matrix = np.zeros((250,16), dtype = complex)
        j = 0
        k = 250
        incount = 1
        
    electrode += 1
    TRIO -=  1
    print('TRIO',TRIO)
    return Coh
