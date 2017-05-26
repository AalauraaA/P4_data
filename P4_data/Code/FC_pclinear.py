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
import fcompres as fc

data_path = "Sound_files\Chapter_9\\"
data_name = "IIRfiltered_sentence2" 
fs, data = wavfile.read(local_path + data_path + data_name + ".wav",mmap=True) # load the data

#==============================================================================
# Frequency lower signal 
#==============================================================================
"""Assign parameters"""
b =  9.76744 #y axis cuts
a = 0.3023 #gradient
f1 = 500
w_length = 0.02 #Window length of compression
mode = "pclin" #The method used in the f compression. Here: pc linear
data_fcompressed = fc.compress_signal(data,w_length,fs,mode,(a,b,f1))

#==============================================================================
# Save file as output
#==============================================================================
gain = 1

out_path = data_path

""" Save new audio file """
wavfile.write(local_path + out_path + data_name + "_FCpclin.wav", fs, np.int16(data_fcompressed*gain))