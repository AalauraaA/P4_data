# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:05:50 2017

@author: G4-101b
"""
from __future__ import division
import sys
import os

"""import liberaies and files"""
local_path = os.getcwd()[:-36]
lib_path = "Python_code\librares"
sys.path.insert(0, local_path + lib_path)

import numpy as np
import matplotlib.pyplot as plt
import P4

pi = np.pi
cos = np.cos
sin = np.sin
tan = np.tan
wc = (5*pi)/18

"""Define constans and frequencies""" 
fs = 16000
T = 1/fs
dc = 1/np.sqrt(2)
d1 = 0.89125
d2 = 0.31623
w1 = T*2*pi*4600
w2 = T*2*pi*5400
wc = 5000*2*pi*T
WC = tan(wc/2)
N = 7

#==============================================================================
# Functions for calculating filter 
#==============================================================================
"""Calculate poles"""
def s_pol(k,N):
    return WC*np.exp(((1j*pi)/(2*N))*(2*k + N - 1)) #Reduced pole plot

"""Calculate continuous filter"""
def H(s):
    product = 1
    for i in range(1,N+1):
        product = product * (s - s_pol(i,N))
    return (WC**N)/product

"""Bilinear transformation on continuous function"""
def Hbil(z):
    return H((z-1.)/(z+1.))

#==============================================================================
# Scipy lowpass
#==============================================================================
W_digital = 5000
fs = 16000
Nr = 100

wlist = np.linspace(0,pi,Nr)
b, a = P4.butter_lowpass(W_digital,fs, order = 7) # Get coefficients
Hw = [P4.H(i,a,b, mode = "ejw") for i in wlist] # Calculate filter on unit circle

#==============================================================================
# Our lowpass
#==============================================================================
Hz = [Hbil(np.exp(1j*w)) for w in wlist]

#==============================================================================
# Plot
#==============================================================================
plt.plot(wlist,np.abs(Hz),'r*', label='Lowpass design')
plt.plot(wlist,np.abs(Hw), label='Scipy lowpass')
plt.plot([wc],[dc],"ko")
plt.plot([0,w1],[d1,d1],"k-")
plt.plot([w2,pi],[d2,d2],"k-")
plt.xlabel('$\omega$ [rad]', fontsize = 15)
plt.ylabel('$|H(e^{j\omega})|$, (Gain)',fontsize = 15)
plt.legend()
plt.savefig('IIRlowcorrect.png', dpi = 500)
plt.show()
