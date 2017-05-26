# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 15:46:54 2017

@author: G4-101b
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter


#==============================================================================
# Function definition
#==============================================================================

def prewrap(Td,w):
    return (2/Td)*np.tan(w/2)
    
def cot(x):
    return np.cos(x)/np.sin(x)

#==============================================================================
# Filter specification
#==============================================================================

pi = np.pi
fs = 16000 
T = 1/fs 

WC = 5000 
wc = T*2*pi*WC
wc2 = T*2*pi*5000 
#w400 = T*2*pi*200
w1 = T*2*pi*4600
w2 = T*2*pi*5400

d1 = 10**(((-1)/(20))) 
d2 = 10**(((-10)/(20)))
dc = 10**(((-3)/(20)))
d5 = 10**(((-3)/(20)))

""" Prewrap """
wcn = prewrap(T,wc)
w1n = prewrap(T,w1)
w2n = prewrap(T,w2)


#==============================================================================
# Finding the order N
#==============================================================================

N1 = 0.5*np.log10(1/(d1**2) - 1)/np.log10(w1n/wcn)
N2 = 0.5*np.log10(1/(d2**2) - 1)/np.log10(w2n/wcn)

#==============================================================================
# Get poles
#==============================================================================

def pole(k,N):
    return WC*np.exp(((1j*np.pi)/(2*N))*(2*k - 1))*1j #Redueced pole plot

N = 7
poles = [pole(k,N) for k in range(0,int(2*N))]

#==============================================================================
# Butterworth-lowpass filter 
#==============================================================================

def butter_lowpass(cutoff, fs, order=2):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
   
b, a = butter_lowpass(WC,fs, order = 7)

def H(z,a,b):
    M = len(b)
    N = len(a)
    if M != N:
        z1 = np.array([z**(-k) for k in range(M)])
        z2 = np.array([z**(-k) for k in range(N)])
        return b.dot(z1)/(a.dot(z2))
    else:
        z = np.array([z**(-k) for k in range(N)])
        return b.dot(z)/(a.dot(z))
        

w = np.linspace(0,pi,1000)
Hw = [H(np.exp(-1j*i),a,b) for i in w]


#==============================================================================
# Butterwoth-bandpass filter
#==============================================================================

""" Specifications """
Wc1 = 100 
wc1 = T * 2 * pi * Wc1 
Wc2 = 5000 
wc2 = T * 2 * pi * Wc2 

wstop1 = T * 2 * pi
wstop2 = T * 2 * pi * 5400
wpass1 = T * 2 * pi * 500 
wpass2 = T * 2 * pi * 4600 

""" Defining filter """
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
    
b, a = butter_bandpass(Wc1,Wc2,fs,order = 7)

w = np.linspace(0,pi,1000) # x-axis in frequency
Hw = [H(np.exp(-1j*i),a,b) for i in w] # Frequency respons

""" Plot filter """
plt.plot(w, np.abs(Hw))
plt.plot([wc1],[1/np.sqrt(2)],"r.", label = r"$\omega_{c1}$")
plt.plot([wc2],[1/np.sqrt(2)],"r.",label = r"$\omega_{c2}$")
plt.plot([wpass1,wpass2],[d1,d1],"k-")
plt.plot([0,wstop1],[d2,d2],"k-")
plt.plot([wstop2,pi],[d2,d2],"k-")
plt.text(.7,.92,"Passband",fontsize = 15)
plt.text(2.2,.35,"Stopband", fontsize = 15)
plt.text(2.,.7,r"$\omega_{c2}$", fontsize = 15)
plt.text(.1,.7,r"$\omega_{c1}$", fontsize = 15)
#plt.axis([0,0.5,0,2])
#plt.legend()
plt.ylabel(r"$|H(e^{j\omega})|$", fontsize = 15)
plt.xlabel(r"$\omega$ [rad]", fontsize = 15)
plt.savefig("bandpass.png", dpi = 500)
plt.show()
