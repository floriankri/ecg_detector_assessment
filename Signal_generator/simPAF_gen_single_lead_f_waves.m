function fWaves = simPAF_gen_single_lead_f_waves(fWaveLength, fibFreqz, fWavesRMS)
%
% fWaves = simPAF_gen_single_lead_f_waves() returns single lead f-waves,
% modeled using an extended sawtooth f wave model first proposed in 
% Stridh and Sörnmo 2001, and later extended to include a stochastic 
% signal component(Petrėnas et al. 2012).
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%

Fs = 1000;
fWaveLength = fWaveLength+1000;      %Prolong f wave signal
df = 0.25;
ff = 0.2;
I = 3;
al = fWavesRMS*10^(-3);
dal = (fWavesRMS/3)*10^(-3);
fa = 0.2;

for n = 1:fWaveLength
    theta = (2*pi*fibFreqz*n)/Fs + (df/ff)*sin((2*pi*ff*n)/Fs);
    yl_temp = 0;
    for i = 1:I
        al_temp = (2/(i*pi))*(al+dal*sin((2*pi*fa*n)/Fs));
        yl_temp = yl_temp + al_temp.*sin(i*theta);
    end
    fWavesModel(n) = -yl_temp; 
end

%% Add noise component
noise = randn(1, fWaveLength);

Fstop1 = fibFreqz*0.6;                % First Stopband Frequency
Fpass1 = fibFreqz*0.7;                % First Passband Frequency
Fpass2 = fibFreqz*0.9;                % Second Passband Frequency
Fstop2 = fibFreqz;                    % Second Stopband Frequency
Astop1 = 50;                          % First Stopband Attenuation (dB)
Apass  = 1;                           % Passband Ripple (dB)
Astop2 = 50;                          % Second Stopband Attenuation (dB)
match  = 'both';                      % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

Fstop1 = fibFreqz;                    % First Stopband Frequency
Fpass1 = fibFreqz+fibFreqz*0.1;       % First Passband Frequency
Fpass2 = fibFreqz+fibFreqz*0.3;       % Second Passband Frequency
Fstop2 = fibFreqz+fibFreqz*0.4;       % Second Stopband Frequency
Astop1 = 50;                          % First Stopband Attenuation (dB)
Apass  = 1;                           % Passband Ripple (dB)
Astop2 = 50;                          % Second Stopband Attenuation (dB)
match  = 'both';                      % Band to match exactly

h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
Hd2 = design(h, 'ellip', 'MatchExactly', match);

noise1 = filtfilt(Hd.sosMatrix,Hd.ScaleValues, noise);
noise2 = filtfilt(Hd2.sosMatrix,Hd2.ScaleValues, noise);
fWaves = 0.5*fWavesModel + 0.25*rms(fWavesModel)*(noise1/rms(noise1))+ 0.25*rms(fWavesModel)*(noise2/rms(noise2));
fWaves = fWaves(1001:end);            % Remove first second of the produced signal
end
