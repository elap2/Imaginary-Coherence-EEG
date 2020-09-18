# Imaginary-Coherence-EEG

## 1. Background
### Paper selected: Identifying true brain interaction from EEG data using the imaginary part of coherency
Guido Nolte 1, Ou Bai, Lewis Wheaton, Zoltan Mari, Sherry Vorbach, Mark Hallett
Affiliations collapse
Affiliation 1Human Motor Control Section, NINDS, NIH, 10 Center Drive MSC 1428, Bldg 10, Room 5N226, Bethesda, MD 20892-1428, USA. nolteg@ninds.nih.gov

#### https://pubmed.ncbi.nlm.nih.gov/15351371/

## 2. This is the Motor-Code

```python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:31:23 2017

@author: lapi2
"""
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

nfft = 256 # number of bins
timestep = 0.1
Fs = 500 # Frequency
dt = 1/Fs # sampled at 

nchan, segleng = data.shape # row or nchan = 64 channels, columns or segleng = all sampling points
t = np.arange(0,segleng/Fs,dt) 

#pylab.plot(t,data[1,:],'g')

channels_number = np.arange(0,nchan)

"""
500 Hz sampling rate
100 samples/2=50 FFT bins
250 Hz/50 FFT bins≃5 Hz/bin
"""

count = 0
incount = 1
j = 0
k = 250
window = signal.hann(250)
N = len(window)
positions = np.arange(len(window))
T = len(window)/Fs

#frq = np.arange(0,Fs/4) # one side frequency range
Sxy = np.zeros(250, dtype = complex)
Sxx = np.zeros(250, dtype = complex)
Syy = np.zeros(250, dtype = complex)
Coh = np.zeros((nchan,250,nchan), dtype = complex)


#freq = np.fft.fftfreq(s1.size, d=timestep)
frq = np.arange(0,Fs/2)*0.5



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
            #print(k)
            
        # PLV = Sxy /  # Sxx*Syy
        #print(index)
        Sxy_mean = Sxy / incount
        Sxx_mean = Sxx / incount
        Syy_mean = Syy / incount
        coherency = Sxy_mean / np.sqrt(Sxx_mean*Syy_mean)
            #xf = np.linspace(0.0, 1.0/(2.0*T), len(window)//2)
            #plt.plot(xf, 2.0/N * np.abs(Coherence[0:N//2].imag))
        Coh[electrode,:,index] = coherency;
            
        Sxy = np.zeros(250, dtype = complex)
        Sxx = np.zeros(250, dtype = complex)
        Syy = np.zeros(250, dtype = complex)
        coherency = np.zeros((250,16), dtype = complex)
        j = 0
        k = 250
        incount = 1
        
    electrode += 1
    TRIO -=  1
    print('TRIO',TRIO)
    #next_electrode += 1

# Plot
#r = Coh[0:nchan-1,0:int(Fs/4),0:nchan-1].imag.T * 2 #imag and transpose
#output = np.where(r>0.3,r,0) # output of imaginary coherence >0.5
r = abs(Coh.imag) #imag and transpose
output = np.where(r>0.25,r,0)

plt.imshow(r[:,4,:])

fig,axs = plt.subplots(2,1)
axs[0].plot(t,data[0,:],t,data[4,:])
axs[1].plot(frq,output)

n = np.arange(0,window.size)
```
