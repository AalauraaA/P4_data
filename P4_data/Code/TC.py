# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:22:22 2017

@author: G4-101b
"""
from __future__ import division

import numpy as np
from scipy.io import wavfile 
import os
import sys

"""import liberaies and files"""
local_path = os.getcwd()[:-15]
data_path = 'P4\P4_data\Sound_files\Chapter_9\\'
git_path = 'P4_data\Sound_files\Chapter_9'
sys.path.insert(0, local_path + git_path)

import fcompres as f

data_name = 'IIRfiltered_sentence1'
fs, x = wavfile.read(local_path + data_path + data_name + ".wav",mmap=True) # load the data


#==============================================================================
# Perform compression
#==============================================================================
ts =[]
ks = []
K = 1.3 
w_length = 0.125 #Window length  
gain = 0.7

xcomp = f.t_comp_gen(x,K,w_length,fs) #Time compression function

output_path = data_path

#==============================================================================
# Save file as output
#==============================================================================
""" Save new audio file """
wavfile.write(local_path + output_path + data_name + "_TC" + ".wav", fs, np.int16(gain*xcomp))
