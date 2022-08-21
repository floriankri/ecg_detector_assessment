function [signal] = gen_PPG_artifacts(duration_sampl,prob,dur_mu,RMS_shape,RMS_scale,slope,Fd)
%
% [] = gen_PPG_artifacts() generates PPG artifacts.
%
% Inputs: 
% duration_sampl - duration of signal to generate in samples
% prob - probabilities of artifacts, [1x4] vector
% dur_mu - mean duration of artifact-free and artifact intervals in seconds, [1x5] vector
% RMS_shape - shape parameter for gamma distribution of normalized RMS artifact amplitudes, [1x4] vector
% RMS_scale - scale parameter for gamma distribution of normalized RMS artifact amplitudes, [1x4] vector
% slope - generated artifact PSD slope, [1x4] vector
% Fd - sampling frequency in Hz
% artifacts in vectors are in the folowing order:
% 1 - device displacement
% 2 - forearm motion
% 3 - hand motion
% 4 - poor contact
%
% Copyright (C) 2020  Birute Paliakaite
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
signal = [];	% empty vector to place signal with artifacts
states_vec = [];	% empty vector to save states
dur_mu = dur_mu*Fd;     % convert mean inerval duration to samples 

fv = linspace(0, 1, 100);	% frequency vector for filter design
for n=1:4
    a = [0 0 log10(fv(3:end)*(Fd/2))*slope(n)];
    a = sqrt(10.^(a/10));
    b(n,:) = firls(250, fv, a); % create filter with predefined slope
end

% Model transitions between artifact-free intervals and artifacts. The transition from the artifact-free interval to all four artifact
% types is possible; however, only the transition to the artifact-free interval is allowed from the artifact.
P = [0 prob; 1 zeros(1,4); 1 zeros(1,4); 1 zeros(1,4); 1 zeros(1,4)]; % stochastic transition matrix
X = 1;	% starts with an initial state corresponding to artifact-free interval
% States are: 1 - artifact-free interval, 2 - device displacement, 3 - forearm motion, 4 - hand motion, 5 - poor contact

while length(states_vec) < duration_sampl % generate for defined signal duration
    T = round(exprnd(dur_mu(X)),0); % duration of the current interval
    if T < Fd
        T = Fd; % the shortest interval allowed is 1 second
    end
    if duration_sampl-(length(states_vec)+T) < Fd
        T = duration_sampl-length(states_vec); % lengthen the current interval if less than 1 second is left to achieve duration
    end
    if X > 1 % if state is with an artifact
        tr = filter(b(X-1,:), 1, randn(1, T+200)); % generate an artifact by filtering white noise
        tr = tr(201:end);
        tr = (tr-mean(tr))/std(tr); % standardize
        R = gamrnd(RMS_shape(X-1),RMS_scale(X-1)); % generate corresponding RMS level 
        tr = R*tr; % tune amplitude
        signal = [signal tr]; % add to already generated artifact signal 
    else
        signal = [signal zeros(1,T)]; % if state is artifact-free, add 0's
    end
    states_vec = [states_vec (X-1)*ones(1,T)]; % add the current state to already generated states
    u = rand; % generate next state
    X = find(u < cumsum(P(X,:)),1);
end
    signal = signal(1:duration_sampl);
    states_vec = states_vec(1:duration_sampl);
end