import scipy.signal as sp
import numpy as np


def highpass(signal, fs):
    # filter cut-offs for PPG and BP
    lpf_cutoff = 0.7 # Hz
    hpf_cutoff = 10 # Hz
    sos_filter = sp.butter(10, [lpf_cutoff, hpf_cutoff],
                       btype = 'bp',
                       analog = False,
                       output = 'sos',
                       fs = fs)

    signal_filt = sp.sosfiltfilt(sos_filter, signal)
    return signal_filt

def u3_transform(signal, fs):
    signal_filt = highpass(signal, fs)
    N = int(0.1*fs)
    u3_all = []
    qrs_location_all = []
    
    u30 = 0
    y = signal_filt
    for j in range(2,int(N+1)):
        u30 += (y[j] - y[j-2])**2
    u3 = [u30]
    for i in range(1,len(y)-N//2):
        u3.append(u3[i-1] + (y[i+N//2] - y[i+N//2-2])**2 - (y[i-N//2] - y[i-N//2+2])**2)
    R_peaks = panPeakDetect(u3, fs)
        
    return u3, R_peaks, signal_filt

#adapted from ecgdetectors library but using theirs lead to errors with low quality ECGs
def panPeakDetect(detection, fs): 
    min_distance = int(0.25*fs)

    signal_peaks = [0]
    noise_peaks = []

    SPKI = 0.0
    NPKI = 0.0

    threshold_I1 = 0.0
    threshold_I2 = 0.0

    RR_missed = 0
    index = 0
    indexes = []

    missed_peaks = []
    peaks = []

    for i in range(1,len(detection)-1):
        if detection[i-1]<detection[i] and detection[i+1]<detection[i]:
            peak = i
            peaks.append(i)

            if detection[peak]>threshold_I1 and (peak-signal_peaks[-1])>0.3*fs:
                    
                signal_peaks.append(peak)
                indexes.append(index)
                SPKI = 0.125*detection[signal_peaks[-1]] + 0.875*SPKI
                if RR_missed!=0:
                    if signal_peaks[-1]-signal_peaks[-2]>RR_missed:
                        missed_section_peaks = peaks[indexes[-2]+1:indexes[-1]]
                        missed_section_peaks2 = []
                        for missed_peak in missed_section_peaks:
                            if missed_peak-signal_peaks[-2]>min_distance and signal_peaks[-1]-missed_peak>min_distance and detection[missed_peak]>threshold_I2:
                                missed_section_peaks2.append(missed_peak)

                        if len(missed_section_peaks2)>0:
                            signal_missed = [detection[i] for i in missed_section_peaks2]
                            index_max = np.argmax(signal_missed)
                            missed_peak = missed_section_peaks2[index_max]
                            missed_peaks.append(missed_peak)
                            signal_peaks.append(signal_peaks[-1])
                            signal_peaks[-2] = missed_peak   

            else:
                noise_peaks.append(peak)
                NPKI = 0.125*detection[noise_peaks[-1]] + 0.875*NPKI

            threshold_I1 = NPKI + 0.25*(SPKI-NPKI)
            threshold_I2 = 0.5*threshold_I1

            if len(signal_peaks)>8:
                RR = np.diff(signal_peaks[-9:])
                RR_ave = int(np.mean(RR))
                RR_missed = int(1.66*RR_ave)

            index = index+1      

    signal_peaks.pop(0)
    return signal_peaks

