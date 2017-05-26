# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:57:06 2017

@author: G4-101b
"""

from __future__ import division
import sys
import os

"""import liberaies and files"""
local_path = os.getcwd()[:-4]
lib_path = "Code\librares"
sys.path.insert(0, local_path + lib_path)
data_path = "Sound_files\Chapter_8\\"

import P4
import numpy as np
from scipy.io import wavfile # get the api

data_name = 'Raw_word'
fs, data = wavfile.read(local_path + data_path + data_name + ".wav",mmap=True) # load the data

""" Assign parameters """
T = 1/fs
Wc1 = 100
Wc2 = 5000

N = 7 #result from get_butterworth

b, a = P4.butter_bandpass(Wc1,Wc2,fs,order = N) #Get filter coefficients

#==============================================================================
# Use the Bandpass-Butterworth filter on the signal
#==============================================================================
x_filter = P4.use_diff(data,a,b) # Filter with coefficents 

gain = 1

""" Save new audio file """
out_path = data_path
wavfile.write(local_path + out_path + "IIRfiltered_word.wav", fs, np.int16(x_filter*gain))

print "Job's done :D" 