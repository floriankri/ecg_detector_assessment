function [y] = simPAF_AR_modeling(rr)
%
% [] = simPAF_AR_modeling() generates slow changing RR component
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%

% AR filter coefficients
b = 3.02093251536985e-07;
a = [1 -2.92818743236790 2.85647129690703 -0.928283828650815];

% Resampling RR series to 2 Hz
N = length(rr); 
fs = 2; 
% Create indexes
ind = rr(1);
for i = 2:N
    ind = [ind ind(i-1)+rr(i)];
end
 
% Linear interpolation
ind2 = ind(1):1/fs:ind(end);
rr2  = interp1(ind,rr,ind2); 
rr2 = rr2 - mean(rr2);

noise = randn(1, length(rr2));
y = filter(b,a,noise);

y = 0.2*(2*(y - min(y))/(max(y) - min(y)) - 1);
y = resample(y, length(rr), length(rr2))';

end