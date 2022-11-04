import numpy as np
from scipy import signal
from ts2vg import NaturalVG

class VisGraphDetector:
    """
    "Fast and Sample Accurate R-Peak Detection for Noisy ECG Using
    Visibility Graphs"
    License: 
    Authors: Taulant Koka, Michael Muma
    """

    def __init__(self, sampling_frequency=250):
        """
        Takes in the sampling rate used for acquisition of the ECG-signal.
        """

        self.fs = sampling_frequency
        

    def ts2vg_adjacency(self,data):
        """
        Computes the adjacency matrix to the directed visibility graph of the input sequence. 
        The direction of edges is taken with respect to the sample amplitude, i.e. the maximum is the root and the minima are the sinks of the graph.
        Uses a fast implementation of the divide and conquer algorithm for visibility graph construction from [1].
        
        [1] https://github.com/CarlosBergillos/ts2vg

        Args:
            data (array): Input sequence.

        Returns:
            adjacency (array): Adjacency matrix of shape N x N of the directed visibility graph.
        """

        size = len(data)
        nvg = NaturalVG()
        nvg.build(data)
        edgelist = nvg.edges
        adjacency = np.zeros((size,size))
        for edge in edgelist: 
            adjacency[edge[0]][edge[1]] = 1
        return adjacency

    def panPeakDetect(self, detection, thresh1=0,thresh2=0):    
        """
        This function implements the thresholding introduced by Pan and Tompkins in [2]. Full credit for this implementation goes to [3].
        [2] J. Pan and W. J. Tompkins, “A Real-Time QRS Detection Algorithm,” IEEE Transactions on Biomedical Engineering, vol. BME-32, Mar. 1985. pp. 230-236.
        [3] https://github.com/berndporr/py-ecg-detectors

        """
        min_distance = int(0.25*self.fs)

        signal_peaks = [0]
        noise_peaks = []

        SPKI = 0.0
        NPKI = 0.0

        threshold_I1 = thresh1
        threshold_I2 = thresh2

        RR_missed = 0
        index = 0
        indexes = []

        missed_peaks = []
        peaks = []

        for i in range(len(detection)):

            if i>0 and i<len(detection)-1:
                if detection[i-1]<detection[i] and detection[i+1]<detection[i]:
                    peak = i
                    peaks.append(i)

                    if detection[peak]>threshold_I1 and (peak-signal_peaks[-1])>0.3*self.fs:
                            
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
                                    missed_peak = missed_section_peaks2[np.argmax(detection[missed_section_peaks2])]
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

    def highpass(self, lowcut = 3, order=2):
        ''' 
            Implements a highpass Butterworth filter. Takes in cutoff frequency 'lowcut', sampling frequency 'fs' and the order.
            Returns numerator (b) and denominator (a) polynomials of the IIR filter. 
        '''
        nyq = 0.5 * self.fs
        high = lowcut / nyq
        b, a = signal.butter(order, high,btype='highpass')
        return b, a

    def calc_weight(self, s,beta):
        """
        This function computes the weights for the input signal s using its visibility graph transformation.
        The weights are calculated by iterative multiplication of the inital weights (1,1,..,1) with the adjacacency matrix of the visibility graph 
        and normalization at each iteration. Terminates when a predetermined amount of weights are equal to zero and returns the weights.

        Args:
            s (array): Signal for which to compute the weights.
            beta (float): Sparsity parameter between 0 and 1. Defines the termination criterion.

        Returns:
            w (array): Weights corresponding to Signal s.
        """
        nvg = self.ts2vg_adjacency(s)
        w = np.ones(len(s))
        while np.count_nonzero(w)> beta*len(s):
            Av = nvg @ w
            w_new =  Av / np.linalg.norm(Av)
            if np.any(np.isnan(w_new)):
                break
            w = w_new
        return w

    def visgraphdetect(self, sig, beta=0.55, gamma=0.5, lowcut=4.0, M = 500):
        """
        This function implements a R-peak detector using the directed natural visibility graph.
        Takes in an ECG-Signal and returns the R-peak indices, the weights and the weighted signal.
        
        Args:
            sig <float>(array): The ECG-signal as a numpy array of length N.
            beta (float, optional): Sparsity parameter for the compuation of the weights. Defaults to 0.55.
            gamma (float, optional): Overlap between consecutive segments in the interval (0,1). Defaults to 0.5.
            lowcut (float, optional): Cutoff frequency of the highpass filter in Hz. Defaults to 4.
            M (int, optional): Segment size. Defaults to 500.

        Returns:
            R_peaks <int>(list): List of the R-peak indices. 
            weights <float>(array): Array of length N containing the weights for the full signal.
            weighted_signal <float>(array): Array of length N containing the weighted filtered signal.
        """
        
        ### filter the signal with a highpass butterworth filter of order 2 ###
        b, a = self.highpass(lowcut)
        sig = signal.filtfilt(b,a,sig)
        N = len(sig)

        ### Initialize some variables ###
        weights = np.zeros(N) # Empty array to store the weights
        l = 0 # Left boundary
        r = M # Right boundary
        dM = int(np.ceil(gamma*M)) # Size of overlap
        L = int(np.ceil(((N-r)/(M-dM)+1))) # Number of segments

        ### Compute the weights for the filtered signal ###
        # for loop is faster
        for jj in range(L):#while right <= N and left<=right:
            s = sig[l:r]
            w = self.calc_weight(s,beta)

            ### Update full weight vector ###
            if l == 0: 
                weights[l:r] = w

            elif N-dM+1 <= l and l+1 <= N:
                weights[l:] = 0.5*(w + weights[l:])

            else: 
                weights[l:l+dM] = 0.5*(w[:dM] + weights[l:l+dM])
                weights[l+dM:r] =  w[dM:] 

            ### break loop, if end of signal is reached ###        
            if r-l < M:
                break

            ### Update segment boundaries ###
            l = l + (M-dM)
            if r +(M-dM) <= N:
                r = r+(M-dM)
            else: 
                r = N

        ### weight the signal and use thresholding algorithm for the peak detection ###
        weighted_signal = sig*weights
        R_peaks = self.panPeakDetect(weighted_signal)
        return R_peaks, weights, weighted_signal