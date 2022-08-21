function rr = simPAF_gen_SR_RR_intervals(N)
%
% rr = simPAF_gen_SR_RR_intervals() returns synthetic SR RR intervals. 
% Ventricular rhythm during SR is simulated using RR interval generator 
% proposed by McSharry et al. 2003 in which both the impact of parasympathetic 
% stimulation (respiratory sinus arrhythmia) and baroreflex regulation
% (Mayer waves) is modeled by a bimodal power spectrum.
%
% Original code:
% Copyright (c) 2003 by Patrick McSharry & Gari Clifford, All Rights Reserved  
% Contact P. McSharry (patrick@mcsharry.net) or G. Clifford (gari@mit.edu)
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details. You should have received a copy of the 
% GNU General Public License along with this program; if not, write to the 
% Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
% MA 02111-1307  USA
%
% Original code is freely availble from Physionet at 
% https://physionet.org/content/ecgsyn/1.0.0/
%
% Modified by: Andrius Petrenas, Biomedical Engineering Institute, 
% Kaunas University of Technology, September 2017
%

if N < 2  % does not work when N = 1
    N = 2;
end
% 

hrMean = randi([50,90]);         % Generate heart rate from an interval of [50-80] bpm
hrStd = randi([5,30])/10;               % Generate SD of heart rate from an interval of [0.5-3] bpm
lfhfRatio = randi([5,20])/10;           % Generate LF/HF ratio [0.5 - 2]
respRate = randi([20,50])/100;          % Respiratory rate [0.2 - 0.5]
MayerFreq = 0.1;                        % Mayer waves 
% Define frequency parameters for RR process 
floStd = 0.01;
fhiStd = 0.01;
% Calculate time scales for RR and total output
rrMean = (60/hrMean);	 
% nRR = 2^(ceil(log2(N*rrMean)));
nRR = 2^(ceil(log2(N)));
% Compute rr process
rr = simPAF_rrprocess(MayerFreq,respRate,floStd,fhiStd,lfhfRatio,hrMean,hrStd,1,nRR);

% Slow changing RR component
y = simPAF_AR_modeling(rr);
rr = rr + y;
end