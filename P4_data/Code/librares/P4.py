# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 13:14:39 2017

@author: Jonas og Tobias 
"""

import numpy as np
from scipy.signal import butter

def FreqBins(N, fs):
    f_bin = [k*fs/float(N) for k in range(N)]
    return f_bin

def timeBins(N, fs):
    t_bin = np.zeros(N)
    for i in range(N):
        t_bin[i] = i * (1./fs)
    return t_bin
    
def Hamming(n,L):
    if 0 <= n <= L:
        return 0.54 - 0.46*np.cos((2*np.pi*n)/L)
    else:
        return 0
def Hann(n,M):
    while(0<=n<=M):
        return 0.5-0.5*np.cos((2*np.pi*n)/M)
    else:
        return 0
def spectrogram(x,ts,window_length,window = Hann,fs = 1):
    """
    x : Is discrete signal as array type
    ts: Is sampling of time signal
    window_length: is... length of window
    window: is window used to truncate signal, set to Hann as default
    fs: It the sampling frequency the signal was obtained by
    
    Returns: time bin, frequency bin and specrogram in matrix 
    
    The first loop takes about 16 times longer to excecute than the second
    
    """
    N = len(x)
    L2 = int(np.ceil(window_length/2)) #Length of FFT and sampled xr
     #gives the highest value of r that can be evaluated
    Nr = int(np.floor(N/float(ts)))
    T = 1./fs #Sampling period
    xr = np.zeros((Nr,window_length)) #Matric for x_r
    #Get array with window
    w_array = [window(m,window_length) for m in range(window_length)] 
    
    #Compute sampled signal with window. This takes 16 longer than next loop
    for r in range(Nr):
        for m in range(window_length):
            if r*ts + m <= N -1:
                xr[r,m] = x[r*ts + m]*w_array[m]
    Xr = np.zeros((Nr,L2),dtype = "complex") #Matric for specgrogram
    #Compute the spectrogram using numpys fft    
    for r in range(Nr):
        Xr[r] = np.fft.fft(xr[r])[:L2]
    Xr = Xr.T
    t = [T*r*ts for r in range(Nr)] #Time bin
    f = [(k*fs)/float(window_length) for k in range(L2)] #Frequency bin
    return np.abs(Xr)**2, t,f #Returns maginitude squared and the bins
    
def butter_bandpass(lowcut, highcut, fs, order=7):
    """
    Conputes coefficients of bandpass filter using scipys butterworth function
    
    
    lowcut : Low cut frequency in Hz 
    highcut: High cut frequency in Hz 
    fs     : is sample frequency
    order  : is order for filter (coefficients will be of lenght ceil(order*2))
    """
    nyq = 0.5 * fs #Nyquist frequency
    low = lowcut / nyq #Lowcut in relation to normalized frequency
    high = highcut / nyq #Same for highcut
    b, a = butter(order, [low, high], btype='band') #Calls scipys butterworth function
    return b, a
    
    
def butter_lowpass(cutoff, fs, order=5):
    """
    Conputes coefficients of lowpass filter using scipys butterworth function
    
    
    lowcut : Low cut frequency in Hz 
    fs     : is sample frequency
    """
    nyq = 0.5 * fs #Nyquist frequency
    normal_cutoff = cutoff / nyq #Lowcut in relation to normalized frequency
    b, a = butter(order, normal_cutoff, btype='low', analog=False)  #Calls scipys butterworth function
    return b, a

def difference_eq(n,i_signal,o_signal,a_coef,b_coef):
    """
    Filters the n'th entry of an input signal and returns output
    
    i_signal: input signal as array type
    o_signal: output signal as array type (valued for up to o[n-1])
    a_coef  : coefficeints for output signal
    b_coef  : coefficeints for input signal 
    n       : index of array x
    """    
    N = len(a_coef)
    M = len(b_coef)
    val = 0 #variable for adding
    if N == M:
       for k in range(M):
           if n - k >= 0:
               val += b_coef[k]*i_signal[n-k] #Output in realation ot input
           if k > 0:
               val -= a_coef[k]*o_signal[n-k] #Input in realtiion to output
       return val #value for y[n] = (h*x)[n]
    elif N == 0:  
        for k in range(M):
            if n - k >= 0:
                val += b_coef[k]*i_signal[n-k] 
        return val
    else:
        print("Implements M != N later")

def use_diff(i_signal,a_coef,b_coef):
    """
    Filters an entire sequence with a LTI filter. Returns output signal
    
    i_signal: input signal as array type
    a_coef  : coefficeints for output signal
    b_coef  : coefficeints for input signal 
    """
    N = len(i_signal)
    y = np.zeros(N, dtype = "float64")
    for n in range(N):
        y[n] = difference_eq(n,i_signal,y,a_coef,b_coef)
    return y
    
def H(z,a,b, mode = "z"): #Discrete frequency response of LTI filter from coefficients
    """
    z   : is complex variable
    a   : is filter coefficients for denominator
    b   : is filter coefficicents for numerator 
    mode: either z or ejw. 
          If "z" evaluate in complex plane
          If "ejw2 evaluate in unit circle
    """
    M = len(b)
    N = len(a)
    if mode == "ejw": #calculatex exponent to the radian
        z = np.exp(1j*z)
    if M != N:
        z1 = np.array([z**(-k) for k in range(M)])
        z2 = np.array([z**(-k) for k in range(N)])
        return b.dot(z1)/(a.dot(z2))
    else:
        z = np.array([z**(-k) for k in range(N)])
        return b.dot(z)/(a.dot(z))

#==============================================================================
# FFT and DFT algorithm
#==============================================================================

def tw(k,n,N):
    """ Calculates twilde factors with parameters k, n and N """
    return np.exp((-1j*np.pi*2*n*k)/N)

def FFT_core(x,X):
    """
    Calulates Discrete Fourier Transform coefficients from signal in time domain
    Returns signal in wrong order so bitrevsing is needed after. 

    x: Is the signal in the timedomain
    X: Is the list the coefficients are saved in (in wrong order)
    s: Is global counting variable  
    """
    N = len(x)
    if N == 2: #Perform length DFT for signals of length 2
        global s
        X[s] = x[0] + x[1]
        X[s + 1] = x[0] - x[1]
        s += 2
    else: 
        """ Devide signal x in to even and odd parts and call this function 
        until all sub sequences are length 2 and can be computed by DFT"""
        N2 = int(N/2)
        a = [x[n] + x[n + N2] for n in range(N2)] #Even coeficients
        b = [(x[n] - x[n + N2])*tw(1,n,N) for n in range(N2)] #Odd coefficients
        FFT_core(a,X)
        FFT_core(b,X)

def bitreverse_step(X): 
    """
    Take 1 sequence and reorder into odd and even pairs by the following 
    algorithm
    
    X: Is input signal
    
    """
    N2 = int(len(X)/2)
    X_even = X[:N2] ; X_odd = X[N2:]
    X_new = np.zeros(len(X),dtype = "complex")    
    for i in range(N2):
        X_new[2*i] = X_even[i]
        X_new[2*i + 1] = X_odd[i]
    return X_new

def bitreverse_sequence(X):
    """
    Takes a sequence in wrong order from FFT and bitreverses it 
    (decimation in frequence)

    X: Is input signal    
    """
    X = np.array(X)
    N = len(X)
    power = int(np.log(N)/np.log(2)) #Get length og signal as power of 2
    powerlist = range(power)
    powerlist.reverse()
    for p in powerlist: #Split signal and rearrange "power" times
        if p > 0:
            #Split sequence into 2**(p-1) sequences that need to be rearranged
            Xsplit = np.array(np.split(X,2**(p-1))) 
            for l in range(len(Xsplit)): #Reverse every sub sequence
                Xsplit[l] = bitreverse_step(Xsplit[l]) 
            X = Xsplit.flatten() #Collaps Xsplit into a 1-dim list
    return X

def fft(x):
    """
    Performs fft of sequences of length 2^i, i is natural number.
    
    x: Is input signal
    
    """
    x = np.array(x)
    N = len(x)
    #Check if N is a pure power of 2
    if str(np.log(N)/np.log(2))[-1] == '0':
        X = np.zeros(N,dtype = "complex") #List for storing FFT
        
        #Calculate fft of X and return uorderd
        global s
        s = 0
        FFT_core(x,X)
        
        #Bitreverse sequence
        X = bitreverse_sequence(X)
        return X
    else:
        print "Input a sequence with lengt radix-2"
    
def DFT(x,k):
    """
    Performs DFT on a signal for a single value of k

    x: Is input sequence
    k: Is desired value evaluated
    
    """
    val = 0
    N = len(x)
    for n in range(N):
        val += x[n]*np.exp((-1j*2*np.pi*k*n)/N)
    return val

def DFT_signal(x):
    """
    Performs DFT on a signal for all values of k
    x: is input sequence
    """
    N = len(x)
    X = np.zeros(N,dtype = "complex")
    for k in range(N):
        X[k] = DFT(x,k) #Performs DFT on a singel value of k
    return X 