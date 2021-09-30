# /******************************************************************************
# *                        ETSI TS 103 634 V1.3.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/

import math as m
import numpy as np

float_epsilon = np.finfo(float).eps

def generate_sinusoid(freq, fs, duration=1, amplitude=0):
    alpha = 2*m.pi*freq/fs
    Nsamples = int(duration * fs)
    fac = 10**(amplitude/20)
    
    sig = [fac * m.cos(alpha*n) for n in range(Nsamples)]
    
    return sig



def _get_notch_iir(alpha, r):
    scf = (1 + 2*r*m.cos(alpha) + r**2)/(2 + 2*m.cos(alpha))

    #print(scf)
    
    A = [1, -2*m.cos(alpha), 1]
    B = [1, -2*r*m.cos(alpha), r**2]

    return A, B

def thd_plus_n(x, x2, freq, fs):
    """thd_plus_n calculates the total harmonic distortion plus noise added in decibel

    x      :   distorted sinusoidal signal

    freq   :   frequency of sinusoid in Hz

    fs     :   sampling frequency in Hz

    return :   ratio of nodge filtered signal energy to signal energy in decibel
    """
    alpha = freq*2*m.pi/fs
    A, B = _get_notch_iir(alpha, .95)
    L = len(x)

    n0 = fs//50
    
    sig_energy = 0
    noise_energy = 0

    ymem = [0, 0]
    
    for i in range(2,L):
        y = A[0]*x[i] + A[1]*x[i-1] + A[2]*x[i-2]  \
            - B[1]*ymem[0] - B[2]*ymem[1]

        if i > n0:
            sig_energy += x2[i]**2
            noise_energy += y**2

        ymem[1] = ymem[0]
        ymem[0] = y
            
    thd_plus_n_ratio = noise_energy/sig_energy

    return float(10*m.log10(thd_plus_n_ratio))


def snr(x, y):
    N = min(len(x),len(y))

    signal_energy, noise_energy = 0, 0

    for i in range(N):
        signal_energy += x[i]**2
        noise_energy += (x[i] - y[i])**2

    signal_to_noise_ratio = signal_energy / noise_energy

    return 10*m.log10(signal_to_noise_ratio)
