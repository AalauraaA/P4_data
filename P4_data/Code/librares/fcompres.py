# -*- coding: utf-8 -*-
"""
Created on Thu May 11 10:09:54 2017

@author: Tobias
"""

from __future__ import division
import numpy as np

pi = np.pi

#==============================================================================
# Methods for time compression
#==============================================================================

def lin_interpolate(x,array):
    """
    Given some potentially non integer value x, the value of the array 
    interpolated between the nearest two points with a linear line is obtained.
    
    x    : Is the value interpolated at. This can be an integer value or not. Only positive numbers allowed
    array: Is the "function" evaluated.       
    """
    a = int(np.floor(x)) #nearest lower index in array
    b = int(np.ceil(x)) #nearest highter index in array
    return array[b]*(x + 1 - b) + array[a]*(b-x) #interpolate with linear line

def t_comp_window(signal,Np,I):
    """
    Streches window by a potentially non integer value.
    
    signal: Is the values in the window as array type
    I     : Is the index list 
    """
    signal_strech = np.zeros(Np)
    for i in range(Np):
        signal_strech[i] = lin_interpolate(I[i],signal)
    return signal_strech
  

def t_comp_gen(signal,c,wlength,fs):
    """
    Streches a signal by a factor of K. Every second window will be thrown away.
    
    signal : Is the sigal as array type
    c      : Is the strech factor 
    wlength: Is the length of the window in seconds
    """
    Nw = int(np.floor(fs*wlength))       #samples in one window
    N = len(signal)                 #samles in signal
    Np = int(Nw*c)             #samles in streched window
    I = np.linspace(0,Nw-1,Np) #Index of samples in streched window
    operations = int(N/(Nw*c)) #number of windows that have to be streched
    signal = np.array(signal)
    
    for n in range(operations):
        x_wi = signal[int(n*Np):int(n*Np) + Nw] #window slice
        x_wi_strech = t_comp_window(x_wi,Np,I) #Streched window
        signal[int(n*Np):int((n+1)*Np)] = x_wi_strech
        
    return signal

def timecomp(x,K,wlength,fs):
    #Streches window of length "wlength" with factor "K"
    x = np.array(x)
    N = len(x)
    
    Nw = int(np.floor(fs*wlength))
    Np = int(Nw*K) 
    operations = int(np.floor(N/(Np)))
    
    x = np.array(x)
    for n in range(operations):
        x[int(n*Np):int((n+1)*Np)] =  x[int(n*Np):int(n*Np) + Nw].repeat(K)
    return x

#==============================================================================
# Methods for frequency compression  
#==============================================================================

def k_bin(fk):
    N = f_w_length*f_fs
    return int(np.round(N*fk/f_fs))

def f_bin(k):
    N = f_w_length*f_fs
    return int(k*f_fs/N)

def dialation_nonlin(k,K,a):
    return k_bin(((f_fs/(pi*K))*np.arctan(((1+a)/(1-a))*np.tan(f_bin(k)*pi/f_fs))))

def dialation_pclin(k,a,b,f_1,flag = False):
    if k > k_bin(f_fs/2):    
        return dialation_pclin(int(2*(k_bin(f_fs/2)) - k),a,b,f_1,True)
    else:
        if k < k_bin(f_1):
            k_new = k
        else:
            k_new = int(round(k*a + b))
        if flag == True:
            return -k_new
        else:
            return k_new

def dialation_lin(k,a,b,flag = False):
    if k > k_bin(f_fs/2):    
        return dialation_lin(int(2*(k_bin(f_fs/2)) - k),a,b,True)
    else:
        k_new = int(round(k*a + b))
        if flag == True:
            return -k_new
        else:
            return k_new
            
def dialation_badlin(k,a,b):
        return int(round(k*a + b))

def compress_window(X,g_list):
    N = len(X)
    X_compressed = np.zeros(N, dtype = "complex")
    for k in range(N):
        X_compressed[g_list[k]] += X[k]
    return X_compressed

def compress_signal(x,w_length,fs,mode,parameters):
    """
    x         : Is the signal
    w_length  : Is the window length in seconds
    fs        : Is the sample rate
    mode      : Is the type of compression wanted
    parameters: Are are the parameters for the compression as tuple type 
    The modes are: "lin", "pclin" and "nonlin" for linenar, piecewise linear 
    and nonlinear
    
    The parameters must be of the form:
    lin: SKRIV NOGET
    pclin: SKRIV NOGET
    nonlin: (K,a) where K is the compression factor and a is the wraping factor
    """
    global f_fs, f_w_length #Set global variable to avoid passing to functions
    f_fs = fs; f_w_length = w_length        
    
    N = len(x)
    samples = int(w_length*fs)
    n_windows = int(N/samples)
    x_comp = np.zeros(N)

    #Calculate compress functions in advance

    if mode == "badlin":
        a , b  = parameters
        g_list = [dialation_badlin(i,a,b) for i in range(samples)]
    
    #Linear
    if mode == "lin":
        a , b  = parameters
        g_list = [dialation_lin(i,a,b) for i in range(samples)]
       
    
    #Piecewise linear
    elif mode == "pclin":
        a , b , f1 = parameters
        g_list = [dialation_pclin(i,a,b,f1) for i in range(samples)]
        
    #Nonlinear method
    elif mode == "nonlin":     
        K , a = parameters
        #list with values in dialation function
        g_list = [dialation_nonlin(i,K,a) for i in range(samples)] 

    for n in range(n_windows):
        x_w = x[n*samples:(n+1)*samples]
        X_w = np.fft.fft(x_w)
        X_com = compress_window(X_w,g_list)
        x_com = np.real(np.fft.ifft(X_com))
        x_comp[n*samples:(n+1)*samples] = x_com
    return x_comp
