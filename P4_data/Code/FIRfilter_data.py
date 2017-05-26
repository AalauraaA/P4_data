# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 10:00:38 2017

@author: G4-101b
"""
from __future__ import division
import sys
import os

"""import liberaies and files"""
local_path = os.getcwd()[:-4]
git_path = 'Code\librares'
sys.path.insert(0, local_path + git_path)
import P4
import numpy as np
from scipy import signal
from scipy.io import wavfile


data_path = "Sound_files\Chapter_8\\"
data_name = "Raw_word"
fs, data = wavfile.read(local_path + data_path + data_name + '.wav',mmap=True)

""" Assign parameters """
nyq = fs/2 #Nyquist
Lowcut = 100/nyq
Highcut = 5000/nyq
f = [Lowcut, Highcut]
K = 8
N = 2**K

#==============================================================================
# FIR filter
#==============================================================================
b = signal.firwin(N, f, window='boxcar', pass_zero=False)
freq = P4.FreqBins(N,fs)
x_filter = P4.use_diff(data,[],b)

#==============================================================================
# Save file as output
#==============================================================================
gain = 1

""" Save new audio file """
out_path = data_path
wavfile.write(local_path + out_path + "FIRfiltered_word.wav", fs, np.int16(x_filter*gain))