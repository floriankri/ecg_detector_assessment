close all
clear 
clc

%% Paroxysmal atrial fibrillation (PAF) generator
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% Related publication: Petrenas, A., Marozas, V., Solosenko, A., Kubilius, R., 
% Skibarkiene, J., Oster, J., & Sornmo, L. (2017). Electrocardiogram modeling 
% during paroxysmal atrial fibrillation: Application to the detection of brief 
% episodes. Physiological Measurement, 38(11), 2058–2080. 
% http://doi.org/10.1088/1361-6579/aa9153
% For more information see help simPAF_ECG_generator!!!

% ---Initial parameters---
rrLength = 50;     % A desired ECG signal length (the number of RR intervals) 
APBrate = 0.10;    % Rate of atrial premature beats (APB). A number between 0 and 0.5
onlyRR = 0;        % 1 - only RR intervals are generated, 0 - multilead ECG is generated

medEpis = 15;       % Median duration of an atrial fibrillation (AF) episode
stayInAF = 1-log(2)/medEpis;   % Probability to stay in AF state
AFburden = 0.8;     % AF burden. 0 - the entire signal is sinus rhythm (SR), 1 - the entire signal is AF

noiseType = 0;      % Type of noise. A number from 0 to 4. 0 - no noise added (noise RMS = 0 mV), 
                    % 1 - motion artefacts, 2 - electrode movement artefacts, 3 - baseline wander, 
                    % 4 - mixture of type 1, type 2 and type 3 noises
noiseRMS = 0.02;    % Noise level in milivolts 

realRRon = 1;       % 1 - real RR series are used, 0 - synthetic
realVAon = 1;       % 1 - real ventricular activity is used, 0 - synthetic
realAAon = 1;       % 1 - real atrial activity is used, 0 - synthetic
% Note: cannot select real atrial activity and synthetic ventricular activity

noiseType_array = [0 1 2 3 4];



load('DATA_PQRST_real');
load('DATA_f_waves_real');
DATAnoises = load('DATA_noises_real'); 


for n = 1 : length(noiseType_array)
% ---PAF generator--- 
    [simPAFdata, initialParameters] = Adapted_simPAF_ECG_generator(rrLength, realRRon, realVAon, realAAon, AFburden, stayInAF, APBrate, noiseType_array(n), noiseRMS, onlyRR, DATApqrst, DATAfWaves, DATAnoises);
    outputdata(n) = simPAFdata;
end

clear DATA_f_waves_real
clear DATApqrst
clear DATAnoises

% ---Returned data and initial parameters---
% Initial parameters
fibFreqz = initialParameters.fibFreqz;      % Dominant fibrillatory frequency
% Data
rr = simPAFdata.rr;                         % RR intervals in miliseconds
multileadECG = simPAFdata.multileadECG;     % 15-lead ECG
multileadVA = simPAFdata.multileadVA;       % 15-lead ventricular activity
multileadAA = simPAFdata.multileadAA;       % 15-lead atrial activity
multileadNoise = simPAFdata.multileadNoise; % 15-lead physiological noise
QRSindex = simPAFdata.QRSindex;             % QRS index. Shows sample number when R peak occurs
targets_SR_AF = simPAFdata.targets_SR_AF;   % Marks SR and AF beats (1 - AF, 0 - SR)
targets_APB = simPAFdata.targets_APB;       % Marks atrial premature beats (1 - APB, 0 - not APB)
pafBoundaries = simPAFdata.pafBoundaries;   % Start and the end of each PAF episode
pafEpisLength = simPAFdata.pafEpisLength;   % Length (in beats) of each PAF episode
ecgLength = simPAFdata.ecgLength;           % ECG length in samples


%% version 1.3.1 of April 29, 2022