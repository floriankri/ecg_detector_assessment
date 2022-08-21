function [ect_loc] = ectopic_detection(Rpos)
% Input
% Rpos - Position of the R peaks
%
% Output
% ect_loc  - Position of the ectopic peaks. Used for respiration analysis
%
% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 17/02/2019   JM      Initial version
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

% Set basic variables
a = 0.1; % 10% difference
b = 0.05; % 5% difference
max_k = 5; % Maximum ammount of previous RR-intervals

% Get the RR-intervals
RR = diff(Rpos);

% Get the amount of RR-intervals
len = length(RR);

% Pre-allocate
RRref_all = zeros(1,len); % Reference RR intervals
ectopic = zeros(1,len); % Index of ectopics

% Set first RR-reference
RR_ref = median(RR(1:max_k));
RRref_all(1) = RR_ref;

% Set loop start
counter = 1;

% Loop over the different RR-intervals
while counter < len-1
    % Define the RR-ref if you are further than the first RR-interval
    if counter > 2
        
        % Set variables
        k = 1; % Number of previous RR-intervals
        sum_ref = 0; % Sum of the good previous RRI's
        num = 0; % Number of good previous RRI's
        wght = [0.5 0.3 0.2]; % Weights: the farther away, the lower the weight
        
        while counter-k > 0 && num < 3 && k < max_k
            
            % Check if the previous RRI is an ectopic
            if ectopic(counter-k) == 1
                prob_check = 1;
            elseif ectopic(counter-k) == 0
                prob_check = 0;
            end
            
            % If this is not the case, use it for reference
            if ~prob_check
                num = num + 1;
                sum_ref = sum_ref + wght(num)*RR(counter-k);
            end
            k = k+1;
        end
        
        % If no good previous RR-intervals are found, this probably means that
        % this is the new normal, so just take the three previous
        if num <= 1
            k = 1;
            while counter-k > 0 && num < 3
                num = num + 1;
                sum_ref = sum_ref + wght(num)*RR(counter-k);
                k = k+1;
            end
        end
        
        % Define the reference variables
        RR_ref = sum_ref/sum(wght(1:num));
        no_RR_ref = nnz(RRref_all);
        RRref_all(no_RR_ref+1) = RR_ref;
    end
    
    % Check for ectopics
    if ((RR(counter) < (1-a)*RR_ref) && (RR(counter+1) > (1+b)*RR_ref)) ||...
            ((RR(counter) > (1+b)*RR_ref) && (RR(counter+1) < (1-a)*RR_ref))
        
        % Indicate the presence of an ectopic
        ectopic(counter) = 1;
    end
    
    % Add to the counter
    counter = counter + 1;
end

% Define the ectopic locations
ect_loc = find(ectopic == 1)+1;
