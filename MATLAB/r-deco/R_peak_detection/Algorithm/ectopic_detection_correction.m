function [Rpos, RR, Pod] = ectopic_detection_correction(fs,RR,Rpos)
% Input
% fs   - sampling frequency
% RR   - RR signal
% Rpos - Position of the R peaks
%
% Output
% Rpos - New position of the peaks
% RR   - Corrected RR interval 
% Pod  - Position of the ectopic peaks. Used for respiration analysis
%        The two beats are replaced, even if the last one is ok.
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

% Pre-allocate
Pod = [];  

% Remove ectopics first 
RRmed = medfilt1(RR,5);   % Fifth order median filter
for ii = 6:length(RR)-1
   if((RR(ii) < RRmed(ii)*0.9) && (RR(ii+1) > RRmed(ii)*1.05)) ...
       || ((RR(ii) > RRmed(ii)*1.05) && (RR(ii+1) < RRmed(ii)*0.9))% one way down, next one way up and viceversa

        % Adjust the investigated RR-interval
        RR(ii) = (RR(ii)+RR(ii+1))/2;
        Rpos(ii+1) = Rpos(ii)+round(RR(ii)*fs/1000);
        
        % Adjust the next RR-interval
        RR(ii+1) = RR(ii);
        Rpos(ii+2) = Rpos(ii+1)+round(RR(ii+1)*fs/1000);
             
        % Store the index of the ectopic beats
        Pod = [Pod; Rpos(ii+1); Rpos(ii+2)];  %#ok 
    end
end

% Correction of ectopic beats. 20% filter is used but for cases when the
% actual interval is small and the next one is large and viceversa
RRmed = medfilt1(RR,5);   % Fifth order median filter
for ii = 6:length(RR)-1
   if((RR(ii) < RRmed(ii)*0.85) && (RR(ii+1) > RRmed(ii)*1.15)) ...
       || ((RR(ii) > RRmed(ii)*1.15) && (RR(ii+1) < RRmed(ii)*0.85))% one way down, next one way up and viceversa
   
        % Adjust the investigated RR-interval
        RR(ii) = (RR(ii)+RR(ii+1))/2;
        Rpos(ii+1) = Rpos(ii)+round(RR(ii)*fs/1000);
        
        % Adjust the next RR-interval
        RR(ii+1) = RR(ii);
        Rpos(ii+2) = Rpos(ii+1)+round(RR(ii+1)*fs/1000);
        
        % Store the index of the ectopic beats
        Pod = [Pod; Rpos(ii+1); Rpos(ii+2)];  %#ok 
    end
end

% Clean and sort Pod
Pod = sort(unique(Pod));



