import numpy as np
import scipy.signal as sig
import scipy.optimize as opt
import math
import matplotlib.pyplot as plt
# import h5py
from functools import partial
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import os
print(os.getcwd())

def peak_search(f, z, fwindow=5e-4, start_f=None, stop_f=None, nsig=3, max_N_peaks=10):
    '''
    ## Perform a search for resonance peaks using a high-pass filter that will
    ## identify points that are above some number of sigma (RMS about the mean)
    ## ARGUMENTS
    ## -              f: array of frequency values [GHz]
    ## -              z: array of corresponding z values
    ## -        fwindow: half the window cut around each resonator before fitting [GHz]
    ## -           nsig: nsig*sigma is the threshold for resonator identification
    ## -    max_N_peaks: if more than this many peaks are identified, call the only peak the minimum of the transmission spectrum
    ## -        start_f: lower bound of resonance identification region [GHz]
    ## -         stop_f: upper bound of resonance identification region [GHz]
    ## RETURNS:
    ## -       peaklist: array of indeces corresponding to the peaks
    ## -            mfz: filtered transmission spectrum
    '''

    ## Extract Nyquist frequency [s]
    nfreq = 1/(2*(abs(f[-1]-f[0])/(len(f)-1)))

    ## The frequency corresponding to the expected window size [s]
    evfreq = 1/(2*fwindow) 
    print('evfreq', evfreq)

    ## Butter some bread?
    print('evfreq/nfreq', evfreq/nfreq)
    print('nfreq', nfreq)
    b, a = sig.butter(2, evfreq/nfreq, btype='highpass')

    ## The magnitude of filtered z, The filtfilt part calls a deprication warning for unknown reasons
    mfz = np.sqrt(sig.filtfilt(b, a, z.real)**2 + sig.filtfilt(b, a, z.imag)**2)  

    ## Do some averaging
    mfz = (mfz+np.append(0,mfz[:-1])+np.append(mfz[1:],0)+np.append([0,0],mfz[:-2])+np.append(mfz[2:],[0,0]))/5
    mfz = (mfz+np.append(0,mfz[:-1])+np.append(mfz[1:],0))/3

    ## Record the standard deviation of mfz
    bstd = np.std(mfz)

    ## initialize peaklist
    peaklist  = np.array([], dtype = int)

    ## initialize mx below min
    mx = -np.inf
    peak_pos = 0
    mx_pos = np.nan
    lookformax = False
    delta = nsig*bstd
    gamma = 3*np.mean(mfz[mfz<delta])

    ## Definte the frequency space to search if none is provided
    if (start_f is None) or (start_f<f[0]):
        start_f = f[0]
    if (stop_f  is None) or (stop_f >f[-1]):
        stop_f  = f[-1]

    ## find peaks and add them to peaklist
    for i in range(len(mfz)):
        if (f[i] >= start_f)*(f[i] <= stop_f):
            cp = mfz[i]
            if cp >= mx:
                mx = cp
                mx_pos = i
            if lookformax == True:
                if cp < gamma:
                    peak_pos  = mx_pos
                    peaklist  = np.append(peaklist , peak_pos)
                    lookformax = False
            else:
                # if cp > delta and f[i] > (min(f)+2*fwindow):
                if cp > delta:
                    mx = cp
                    mx_pos = i
                    lookformax = True

    ## Handle too many peaks by picking the minimum of the unfiltered transmission
    if (len(peaklist) > max_N_peaks):
        # peaklist = np.array([ np.argmax(20*np.log10(abs(np.array(z)))) ])
        # peaklist = np.array([ np.argmax(mfz) ])
        peaklist = np.array([ np.argmin(20*np.log10(np.array(z))) ])

    return peaklist, mfz

def plot_filtered_trace_with_peaks(f, z, plot_range, fwindow=5e-4, start_f=None, stop_f=None, nsig=3, max_N_peaks=10):
    peaklist, mfz = peak_search(f, z, fwindow=fwindow, start_f=start_f, stop_f=stop_f, nsig=nsig, max_N_peaks=max_N_peaks)

    # Convert to GHz if needed (assuming input is already in GHz)
    f_plot = f*1e3 # mhz
    
    # Print peak frequencies
    peak_freqs = f_plot[peaklist]
    print("Identified peak frequencies (MHz):", peak_freqs)

    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(f_plot, mfz, label='Filtered magnitude (mfz)')
    plt.plot(f_plot[peaklist], mfz[peaklist], 'ro', label='Identified Peaks')
    for i in peaklist:
        plt.axvline(f_plot[i], color='red', linestyle='--', alpha=0.5)

    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Filtered |S21|')
    plt.xlim(plot_range[0], plot_range[1])
    plt.title('Filtered Transmission with Detected Peaks')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()