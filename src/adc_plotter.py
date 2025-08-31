#!/usr/bin/python

"""
Author: Sravan Munagavalasa
Date: December 2024
Description: A set of tools to quickly view the ADC data comming from the ACM

Update: Sravan (June 2025). Works with new (proper version) of file type with explicit RU and RD
"""

import sys, os, csv
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram


def adc_plotter(fname, trace=False, psd=False, V0=0, t0=0, delta_t=500):
    df = pd.read_csv(fname)
    df["t"] = df.index / 15 # 15MHz Clock
    
    if trace:
        plot_trace(df, V0, t0, delta_t, title_=fname)
    
    if psd:
        plot_psd(df, npsd=5000, nsampint=150, fs=15.E6, title_=fname)


def adc_plotter_compare(fnames, trace=False, psd=False, V0=0, t0=0, delta_t=500):

    # Read the first two files in the list as DataFrames and add the file name as a column
    df_list = [
        pd.read_csv(fname).assign(
            t=lambda df: df.index / 15,  # Add 't' column based on index
            filename=fname  # Add 'filename' column with the file name
        )
        for fname in fnames  # Read only the first two files
    ]


    if trace:
        plot_trace(df_list, t0, delta_t, label_=fnames)
    
    if psd:
        plot_psd(df_list, npsd=5000, nsampint=150, fs=15.E6, label_=fnames)



def plot_trace(df_list, V0_list=0, t0_list=0, delta_t=300, title_="ADC Trace", label_=None):
    """ plot the adc trace from the start point to end point. Pedestal and signal as well"""
    if not isinstance(df_list, list):
        df_list = [df_list]
    
    if not isinstance(t0_list, list):
        t0_list = len(df_list)*[t0_list]
        
    if not isinstance(V0_list, list):
        V0_list = len(df_list)*[V0_list]

    if not label_:
        label_ = len(df_list) * ["Trace"]
    elif not isinstance(label_, list):
        label_ = len(df_list) * [label_]


    ped_colors = ['peru', 'sandybrown', 'chocolate', 'sienna']
    sig_colors = ['forestgreen', 'seagreen', 'springgreen', 'darkseagreen']
    adc_colors = ['dimgrey', 'slategray' , 'CornflowerBlue', 'dodgerblue']
    
    plt.figure(figsize=(10, 5))

    for i, (df, t0, V0) in enumerate(zip(df_list, t0_list, V0_list)):
        # clipping time is convulted but is need to properly line up
        df.t = df.t - t0
        df = df[(df.t >= 0) & (df.t <= delta_t)]
        df.val = df.val - V0

        if "mode" in df:
            df_ped = df[df["mode"]==2]
            df_sig = df[df["mode"]==0]
        elif ("RD" in df) and ("RU" in df):
            df_ped = df[df["RD"]==1]
            df_sig = df[df["RU"]==1]



        # Plotting the adc column within the specified range
        plt.scatter(df_ped.t, df_ped.val, color=ped_colors[i],
                    label='Pedestal' if i==0 else None, zorder=3)
        plt.scatter(df_sig.t, df_sig.val, color=sig_colors[i],
                     label='Signal'   if i==0 else None, zorder=3)
        plt.plot(df.t,     df.val,     color=adc_colors[i], label=label_[i])

            
    # Adding labels and legend
    plt.xlabel('time ($\mu$s)')
    plt.ylabel('ADC Value')
    plt.title(title_)
    plt.legend()


    #plt.show()



def plot_psd(df_list, npsd=5000, nsampint=150, fs=15.E6, title_="Power Spectral Distribution", label_=None):

    if not isinstance(df_list, list):
        df_list = [df_list]

    if not label_:
        label_ = len(df_list) * [""]
    elif not isinstance(label_, list):
        label_ = len(df_list) * [label_+" "]
    else:
        label_ = [label__ + " " for label__ in label_]

    # Make figure and super title
    fig, axs = plt.subplots(2,2,figsize=(14,8.5))
    if title_:
        fig.suptitle(title_, fontsize=16, fontweight='bold')

    # Plot 0,0 is the trace
    axs[0,0].set_title('ADC Trace')
    axs[0,0].set_xlabel('N Sample')
    axs[0,0].set_ylabel('Sample Value (ADU)')

    # Plot 0, 1 is 
    axs[0,1].set_ylabel(r"ADU/$\sqrt{Hz}$")
    axs[0,1].set_xlabel("Frequency (Hz)")
    
    # Plot 1, 0 is the log log PSD
    axs[1,0].set_ylabel(r"ADU/$\sqrt{Hz}$")
    axs[1,0].set_xlabel("Frequency (Hz)")

    # Plot 1, 1 is the linear PSD
    axs[1,1].set_ylabel(r"ADU/$\sqrt{Hz}$")
    axs[1,1].set_xlabel("Frequency (Hz)")

        
    for i, (df, label__) in enumerate(zip(df_list, label_)):
        # adctrace = df.val[800000:1400000]
        adctrace = df.val
        
        # Calculate PSD
        f0, Vxx0 = make_psd(adctrace, fs)
        f,  Vxx  = make_avg_psd(adctrace, fs)

        # calculate CDS noise
        sig_CDS = calculate_CDS(adctrace, nsampint)
        sig_ADU = np.std(adctrace)


        # label2__  = "$\sigma_\\text{trace}$="+f"{sig_ADU:.2e}   "
        # label2__ += "$\sigma_\\text{trace}$="+f"{sig_CDS:.2e}"

        label2__  = "$\sigma_{ADU}$="+f"{sig_ADU:.2e}   "
        label2__ += "$\sigma_{CDS}$="+f"{sig_CDS:.2e}"
        
        axs[0,0].plot(adctrace, '.', label=label__)      # Plot the ADC Trace
        axs[0,1].loglog(f0[1:], Vxx0[1:], label=label__) # Plot non averaged psd
        axs[1,1].plot(f[1:], Vxx[1:], label=label2__)     # Plot PSD in lineaer scale
        axs[1,0].loglog(f[1:], Vxx[1:], label=label2__)   # Plot PSD in log log scale


    plt.tight_layout()
    axs[0,0].legend()
    axs[1,1].legend()
    #plt.show()


def calculate_CDS(adctrace, nsampint):
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/(2*nsampint))-1)
    valCDS = np.zeros(nblocks, dtype=np.float64)
    for isampl in range (0,nblocks):
        valCDS[isampl] = np.mean(adctrace[isampl+nsampint:isampl+2*nsampint])-np.mean(adctrace[isampl:isampl+nsampint])
    sigmaCDS = np.std(valCDS)
    return sigmaCDS


def make_psd(adctrace, samplingFreq=15.e6):
    f, Pxx = periodogram(adctrace-np.mean(adctrace), fs=samplingFreq, scaling='density');
    Vxx=np.sqrt(Pxx)
    return f,Vxx


def make_avg_psd(adctrace, samplingFreq=15.e6, npsd=10000):
    Pxx_avg = np.zeros(int(npsd/2)+1, dtype=np.float64)
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/npsd)-1)
    for isampl in range (0,nblocks):
        istart = isampl*npsd
        iend = istart+npsd
        y_unsigned = adctrace[istart:iend]
        f, Pxx = periodogram(y_unsigned-np.mean(y_unsigned), fs=samplingFreq, scaling='density');
        Pxx_avg = Pxx_avg + Pxx/nblocks

    Vxx=np.sqrt(Pxx_avg)
    return f,Vxx

def calculate_CDS(adctrace,nsampint):
    nsamples = len(adctrace)
    nblocks = int(np.ceil(nsamples/(2*nsampint))-1)
    valCDS = np.zeros(nblocks, dtype=np.float64)
    for isampl in range (0,nblocks):
        valCDS[isampl] = np.mean(adctrace[isampl+nsampint:isampl+2*nsampint])-np.mean(adctrace[isampl:isampl+nsampint])
    sigmaCDS = np.std(valCDS)
    return sigmaCDS


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process some parameters.")
    parser.add_argument('--fname', nargs="+", type=str, help='The file name as a string.')
    parser.add_argument('--compare', action='store_true', help='compare 2 different traces')
    parser.add_argument('--psd', action='store_true', help='Enable PSD (boolean).')
    parser.add_argument('--trace', action='store_true', help='Enable trace (boolean).')
    parser.add_argument('--t0', nargs="+", type=float, default=0, help='Start time as a us')
    parser.add_argument('--V0', nargs="+", type=float, default=0, help='voltage offset')
    parser.add_argument('--delta_t', type=float, default=500, help='Time interval as a number.')

    args = parser.parse_args()
    
    if not args.compare:
        for f in args.fname:
            adc_plotter(f, args.trace, args.psd, args.V0, args.t0, args.delta_t)
            # different file extensions just for the work flow...
            if args.trace==True and args.psd==False:
                outname = os.path.splitext(f)[0] + ".pdf" 
            if args.trace==False and args.psd==True:
                outname = os.path.splitext(f)[0] + ".png"
            plt.savefig(outname, dpi=300, bbox_inches="tight")

            plt.close()   # prevents figures from piling up
    else:
        adc_plotter_compare(args.fname, args.trace, args.psd, args.V0, args.t0, args.delta_t)

    plt.show()

