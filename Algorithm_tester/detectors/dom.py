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
def split(list_a, chunk_size):
    for i in range(0, len(list_a), chunk_size):
        yield list_a[i:i + chunk_size]
def dom(signal, fs):
    '''DOM algorithm'''
    signal_filt = highpass(signal, fs)


    #difference operation
    x_d = []
    for n in range(len(signal_filt)):
        x_d.append(signal_filt[n] - signal_filt[n-1])
    
    #low pass filter
    sos_filter_low = sp.butter(10, 100,
                    btype = 'lowpass',
                    analog = False,
                    output = 'sos',
                    fs = fs)
    x_df = sp.sosfiltfilt(sos_filter_low, x_d)

    
    #step 4 thresholding
    signal_pos = []
    signal_neg = []
    for i in x_d: #should use original ECG but doesn't work - this is where our threshold changes
        if i >= 0:
            signal_pos.append(i)
        if i < 0:
            signal_neg.append(i)
    MVp = np.sum(signal_pos)/len(signal_pos)
    MVn = np.sum(signal_neg)/len(signal_neg)
    T1 = 2*MVp
    T2 = 2*MVn
    x_df_hat = []
    for x in x_df:
        if np.any((x > 0) & (x < T1)) or np.any((x > T2) & (x < 0)):
            x_df_hat.append(0)
        else:
            x_df_hat.append(x)

    #Find point R
    #step 1
    x_df_hat_pos = []
    x_df_hat_neg = []
    for x in x_df_hat:
        if x >= 0:
            x_df_hat_pos.append(x)
            x_df_hat_neg.append(0)
        else:
            x_df_hat_neg.append(x)
            x_df_hat_pos.append(0)

    #step 2

    chunk_size = 50
    x_df_hat_pos_chunks = list(split(x_df_hat_pos, chunk_size))
    x_df_hat_neg_chunks = list(split(x_df_hat_neg, chunk_size))
    max_pos = []
    max_pos_index = []
    for i in range(len(x_df_hat_pos_chunks)):
        max_pos.append(np.max(x_df_hat_pos_chunks[i]))
        max_pos_index.append(np.argmax(x_df_hat_pos_chunks[i])+chunk_size*i)
    max_neg = []
    max_neg_index = []
    for i in range(len(x_df_hat_neg_chunks)):
        max_neg.append(np.min(x_df_hat_neg_chunks[i]))
        max_neg_index.append(np.argmax(x_df_hat_neg_chunks[i])+chunk_size*i)

    max_pos_real = max_pos.copy()
    for i in range(len(max_pos_index)-1):
        if not any(max_pos[i:i+2]) == 0:
            prev_index = max_pos_index[i]
            after_index = max_pos_index[i+1]
            if after_index - prev_index < 50:
                if max_pos_real[i] > max_pos_real[i+1]:
                    max_pos_real[i+1] = 0
                else:
                    max_pos_real[i] = 0
    max_neg_real = max_neg.copy()
    for i in range(len(max_neg_index)-1):
        if not any(max_neg[i:i+2]) == 0:
            prev_index = max_neg_index[i]
            after_index = max_neg_index[i+1]
            if after_index - prev_index < 50:
                if max_neg_real[i] < max_neg_real[i+1]:
                    max_neg_real[i+1] = 0
                else:
                    max_neg_real[i] = 0

    max_pos_collated = np.column_stack([max_pos_real, max_pos_index])
    max_neg_collated = np.column_stack([max_neg_real, max_neg_index])
    pos = []
    neg = []
    for i in range(len(max_pos_collated)):
        if max_pos_collated[i][0] != 0:
            pos.append(max_pos_collated[i][1])
        if max_neg_collated[i][0] != 0:
            neg.append(max_neg_collated[i][1])

    correct_peak = []       
    for i in range(len(pos)):
        upper_limit = pos[i]+50
        lower_limit = pos[i]-50
        if any(y > lower_limit and y < upper_limit for y in neg):
            correct_peak.append(int(pos[i]))
    
    return correct_peak




    