# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:15:56 2017

@author: G4-101b
"""

from __future__ import division
import sys
import os
from scipy.io import wavfile 
import numpy as np 

"""import liberaies and files"""
local_path = os.getcwd()[:-4]
lib_path = "Code\librares"
sys.path.insert(0, local_path + lib_path)
data_path = "Sound_files\Chapter_8\\"
import P4
import fcompres as fc


data_name = "Raw_sentence2" #21 year old female saying: "Michael ejer tyve pæne ringe"
fs, data = wavfile.read(local_path + data_path + data_name + ".wav",mmap=True) # load the data

#==============================================================================
# Filter signal
#==============================================================================

Wc1 = 100 #Lower cutooff frequency
Wc2 = 5000 #Higher cutooff frequency
N = 7 #Filter order

b, a = P4.butter_bandpass(Wc1,Wc2,fs,order = N) #Get filter coefficients
data_filter = P4.use_diff(data,a,b) # Filter with coefficents 

#==============================================================================
# Frequency lower signal 
#==============================================================================
"""Assign parameters"""
K = 2.2 #Compression factor
a = 0.4 #Wraping factor
w_length = 0.02 #Window length of compression
mode = "nonlin" #The method used in the f compression. Here: Non linear

data_fcompressed = fc.compress_signal(data_filter,w_length,fs,mode,(K,a))

#==============================================================================
# Save file as output
#==============================================================================
gain = 1

out_path = "Sound_files\output_for_scripts\\"

""" Save new audio file """
wavfile.write(local_path + out_path + data_name + "_compressed.wav", fs, np.int16(data_fcompressed*gain))