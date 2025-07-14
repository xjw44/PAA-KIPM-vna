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

def peak_search(f, z, save_path, fwindow=5e-4, start_f=None, stop_f=None, nsig=3, max_N_peaks=10):
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
    # print('evfreq', evfreq)

    # ## Butter some bread?
    print('evfreq/nfreq', evfreq/nfreq)
    # print('nfreq', nfreq)
    b, a = sig.butter(2, evfreq/nfreq, btype='highpass')

    w, h = sig.freqz(b, a)
    plt.plot(w, 20 * np.log10(abs(h)))
    plt.title('Butterworth filter frequency response')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Amplitude [dB]')
    plt.margins(0, 0.1)
    plt.grid(which='both', axis='both')
    save_dir = os.path.dirname(save_path)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path+'filter_resp.pdf', dpi=300)
    plt.savefig(save_path+'filter_resp.png', dpi=300)
    plt.show()
    plt.close()

    ## The magnitude of filtered z, The filtfilt part calls a deprication warning for unknown reasons
    # mfz = z
    # mfz = sig.filtfilt(b, a, z.real)
    mfz = np.sqrt(z.real**2 + z.imag**2)  
    # mfz = np.sqrt(sig.filtfilt(b, a, z.real)**2 + sig.filtfilt(b, a, z.imag)**2)  
    print('before ave', len(mfz))

    # ## Do some averaging
    mfz = (mfz+np.append(0,mfz[:-1])+np.append(mfz[1:],0)+np.append([0,0],mfz[:-2])+np.append(mfz[2:],[0,0]))/5
    mfz = (mfz+np.append(0,mfz[:-1])+np.append(mfz[1:],0))/3
    print('after ave', len(mfz))

    ## Record the standard deviation of mfz
    bstd = np.std(mfz)
    print('bstd', bstd)

    ## initialize peaklist
    peaklist  = np.array([], dtype = int)

    ## initialize mx below min
    mx = -np.inf
    peak_pos = 0
    mx_pos = np.nan
    lookformax = False
    delta = nsig*bstd
    print('delta', delta)
    # gamma = 3*np.mean(mfz[mfz<delta])
    # print('gamma', gamma)
    # print('gamma', mfz<delta)

    # ## Definte the frequency space to search if none is provided
    # if (start_f is None) or (start_f<f[0]):
    #     start_f = f[0]
    # if (stop_f  is None) or (stop_f >f[-1]):
    #     stop_f  = f[-1]

    # # ## find peaks and add them to peaklist
    # # for i in range(len(mfz)):
    # #     if (f[i] >= start_f)*(f[i] <= stop_f):
    # #         cp = mfz[i]
    # #         if cp >= mx:
    # #             mx = cp
    # #             mx_pos = i
    # #         if lookformax == True:
    # #             # if cp < gamma:
    # #             peak_pos  = mx_pos
    # #             peaklist  = np.append(peaklist , peak_pos)
    # #             lookformax = False
    # #         else:
    # #             # if cp > delta and f[i] > (min(f)+2*fwindow):
    # #             if cp > delta:
    # #                 mx = cp
    # #                 mx_pos = i
    # #                 lookformax = True
    # # print('len_peaklist', len(peaklist))

    # ## Handle too many peaks by picking the minimum of the unfiltered transmission
    # if (len(peaklist) > max_N_peaks):
    #     # peaklist = np.array([ np.argmax(20*np.log10(abs(np.array(z)))) ])
    #     # peaklist = np.array([ np.argmax(mfz) ])
    #     peaklist = np.array([ np.argmin(20*np.log10(np.array(z))) ])

    peaklist, properties = sig.find_peaks(mfz, width=[10, 200])
    print('peaklist', len(peaklist))
    print('peaklist', peaklist)
    print('properties', properties)

    return peaklist, mfz

def plot_filtered_trace_with_peaks(f, z, plot_range, plot_range_y, save_path, fwindow=5e-4, start_f=None, stop_f=None, nsig=3, max_N_peaks=10):
    peaklist, mfz = peak_search(f, z, save_path, fwindow=fwindow, start_f=start_f, stop_f=stop_f, nsig=nsig, max_N_peaks=max_N_peaks)

    # Convert to GHz if needed (assuming input is already in GHz)
    f_plot = f*1e3 # mhz
    
    # Print peak frequencies
    peak_freqs = f_plot[peaklist]
    print("Identified peak frequencies (MHz):", peak_freqs)

    z = np.polyfit(f_plot, mfz, 100)
    z_50 = np.polyfit(f_plot, mfz, 50)
    p = np.poly1d(z)
    p_50 = np.poly1d(z_50)

    # Compute the fitted values at original x points
    y_fit = p(f_plot)
    y_fit_50 = p_50(f_plot)

    # Compute residuals
    residuals = mfz - y_fit
    residuals_50 = mfz - y_fit_50

    # residual peaks 
    peaklist_res, properties_res = sig.find_peaks(residuals, width=[10, 200], prominence=[0, 0.02])
    print('res prop', properties_res)
    print('len peaklist', len(peaklist_res))

    # Create subplots: top for trace, bottom for residuals
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6), gridspec_kw={'height_ratios': [1, 1]})

    # --- Top panel: Original trace + fit + peaks ---
    ax1.plot(f_plot, mfz, label='Filtered magnitude (mfz)')
    ax1.plot(f_plot, y_fit, 'r--', label='High Degree Poly Fit')
    ax1.plot(f_plot, y_fit_50, 'r--', label='High Degree Poly Fit')
    ax1.plot(f_plot[peaklist], mfz[peaklist], 'ro', label='Identified Peaks')

    for i in peaklist:
        ax1.axvline(f_plot[i], color='red', linestyle='--', alpha=0.5)

    ax1.set_ylabel('|S21| (dB)')
    ax1.set_xlim(plot_range[0], plot_range[1])
    if plot_range_y: 
        ax1.set_ylim(plot_range_y[0], plot_range_y[1])
    ax1.set_title('Filtered Transmission with Detected Peaks')
    ax1.legend()
    ax1.grid(True)

    # --- Bottom panel: Residuals ---
    ax2.plot(f_plot, residuals, label='Residuals (mfz - polyfit)')
    ax2.plot(f_plot, residuals_50, label='Residuals (mfz - polyfit)')
    ax2.axhline(0, color='k', linestyle='--', linewidth=1)
    ax2.set_ylim(-0.02, 0.02)

    for i in peaklist_res:
        ax2.axvline(f_plot[i], color='red', linestyle='--', alpha=0.5)

    ax2.set_xlabel('Frequency (MHz)')
    ax2.set_ylabel('Residual')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()

    # Save
    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    plt.savefig(save_path + '.pdf', dpi=300)
    plt.savefig(save_path + '.png', dpi=300)
    plt.show()
    plt.close()

    return peaklist_res

