function [signal,check] = filtersignal(signal,fs,type,order,hpf,lpf,stopf)

%  INPUT:   signal  - Signal
%           fs   - Sampling frequency
%           type - This is a vector of ones and zeros. 
%           order - This is a cell with the high pass and low pass filter orders
%                  
%           hpf   - cutoff frequency for highpass filter
%           lpf   - cutoff frequency for lowpass filter
%           stopf - cutoff frequency for band stop filter
%           other - 1 or 0, define a new filter when other ==1
%  OUTPUT:  signal - Filtered signal
%
%
% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 06/05/2019   JM      Initial version
%
% Copyright (c) 2019,  Jonathan Moeyersons, KULeuven-ESAT-STADIUS 
%
% This software is made available for non commercial research purposes only 
% under the GNU General Public License. However, notwithstanding any 
% provision of the GNU General Public License, this software may not be 
% used for commercial purposes without explicit written permission after 
% contacting jonathan.moeyersons@esat.kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-define check
check = 1;

% Set the signal in the correct form (collumn vectors)
[R_temp,K_temp] = size(signal);
if R_temp < K_temp
    signal = signal';
end

% Get the Nyquist frequency
Nf = fs/2;  

% High pass filter
if type(1)
    try
        Cf      = hpf/Nf;        % Cutoff freq. Passband corner freq. 0.5Hz

        [bh,ah] = butter(order{1},Cf,'high');   % High pass filter
        signal  = filtfilt(bh,ah,signal);
    catch ME
        % Show an error dialog
        errordlg(ME.message)
        
        % Change check
        check=0;
    end    
end

% Low pass filter
if type(2)
    try
        Cf      = lpf/Nf;        % Cutoff freq. Passband corner freq. 100Hz

        [bl,al] = butter(order{2},Cf,'low');    % Low pass filter
        signal  = filtfilt(bl,al,signal);
    catch ME
        % Show an error dialog
        errordlg(ME.message)
        
        % Change check
        check=0;
    end 
end

% Notch filter
if type(3)
    try
        % Define zeros
        omega   = stopf/Nf*pi;
        zero(1) = cos(omega)+1i*sin(omega);
        zero(2) = cos(omega)-1i*sin(omega);

        % Define poles
        pole    = 0.99*zero;

        % Define coefficients
        bs      = [1 -zero(1)-zero(2) zero(1)*zero(2)] ;
        as      = [1 -pole(1)-pole(2) pole(1)*pole(2)] ;

        signal  = filtfilt(bs,as,signal);
    catch ME
        % Show an error dialog
        errordlg(ME.message)
        
        % Change check
        check = 0;
    end 
end 

% Self designed filter
if type(4)
    try
        % Design filter
        Filt = designfilt;

        if ~isempty(Filt)
            if strcmp(Filt.FrequencyResponse,'differentiator') || strcmp(Filt.FrequencyResponse,'hilbert')
                 % Filter the signal
                 signal = filter(Filt,signal);
            else
                % Filter the signal
                signal = filtfilt(Filt,signal);   
            end
        else
            % Change check
            check = 0;
        end
    catch ME
        % Show an error dialog
        errordlg(ME.message)

        % Change check
        check = 0;
    end
end


   
    
